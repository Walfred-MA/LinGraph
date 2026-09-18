"""Piecewise base coordinates for insertion representatives encoded on others."""
from bisect import bisect_right
from collections import OrderedDict
import json
import re

from graph_cigar_payloads import query_only_graph_cigar
from tools.grvcf_to_vcf import OP, parse_info, unescape


def cigar_chunks(text):
    text = query_only_graph_cigar(unescape(text))
    chunks, end = [], 0
    for match in re.finditer(r'([<>])([^<>]+)', text):
        if match.start() != end:
            raise ValueError('malformed graph-CIGAR chunk')
        end = match.end()
        direction, body = match.groups()
        name, separator, rest = body.partition(':')
        if separator:
            body = rest
        else:
            name = None
        operations, cursor = [], 0
        for op in OP.finditer(body):
            if op.start() != cursor or int(op[1]) <= 0:
                raise ValueError('malformed graph-CIGAR operation')
            cursor = op.end()
            operations.append((int(op[1]), op[2], op[3].upper()))
        if cursor != len(body) or not operations:
            raise ValueError('malformed graph-CIGAR suffix')
        chunks.append((direction, name, operations))
    if end != len(text) or not chunks:
        raise ValueError('malformed graph CIGAR')
    return chunks


def append_run(runs, start, end, name, target, step):
    if start == end:
        return
    if runs and runs[-1][1] == start and runs[-1][2] == name and runs[-1][4] == step:
        previous = runs[-1]
        if previous[3] + (previous[1] - previous[0]) * step == target:
            previous[1] = end
            return
    runs.append([start, end, name, target, step])


class Coordinates:
    """Flatten transitive alignments into nonoverlapping runs, stored on disk.

    A run [a,b,name,t,step] maps source base x to t + (x-a)*step.
    Query-only I/S bases stay on the source insertion; they have no target base.
    """
    def __init__(self, database, catalog, *, build=False):
        self.db, self.catalog, self.build = database, catalog, build
        self.active, self.cache, self.sequences = set(), OrderedDict(), OrderedDict()
        if build:
            self.db.executescript('CREATE TABLE projections (id TEXT PRIMARY KEY, data TEXT); '
                                  'CREATE TABLE loci (id TEXT PRIMARY KEY, length INTEGER);')

    def length(self, name):
        row = self.db.execute('SELECT length(sequence) FROM resolved WHERE id=?', (name,)).fetchone()
        if row:
            return row[0]
        if name in self.catalog.readers:
            return self.catalog.length(name)
        raise ValueError(f'coordinate target {name!r} has no resolved sequence')

    def reference(self, name, start, end):
        if name in self.catalog.readers:
            return self.catalog.fetch(name, start, end)
        if name not in self.sequences:
            row = self.db.execute('SELECT sequence FROM resolved WHERE id=?', (name,)).fetchone()
            if row is None:
                raise ValueError(f'coordinate target {name!r} has no resolved sequence')
            self.sequences[name] = row[0]
            # Bound cached sequence data; one long current template is allowed.
            while len(self.sequences) > 1 and sum(map(len, self.sequences.values())) > 32 * 1024 * 1024:
                self.sequences.popitem(last=False)
        self.sequences.move_to_end(name)
        sequence = self.sequences[name]
        if not 0 <= start <= end <= len(sequence):
            raise ValueError(f'coordinate interval outside insertion {name!r}')
        return sequence[start:end]

    def segments(self, name):
        if name not in self.cache:
            row = self.db.execute('SELECT data FROM projections WHERE id=?', (name,)).fetchone()
            if row is None:
                runs = self._build(name)
                if self.build:
                    self.db.execute('INSERT INTO projections VALUES (?,?)', (name, json.dumps(runs)))
                    for _, _, destination, _, _ in runs:
                        self.db.execute('INSERT OR IGNORE INTO loci VALUES (?,?)', (destination, self.length(destination)))
            else:
                runs = json.loads(row[0])
            if len(self.cache) >= 64:
                self.cache.popitem(last=False)
            self.cache[name] = runs, [r[1] for r in runs]
        self.cache.move_to_end(name)
        return self.cache[name]

    def project(self, name, start, end):
        if not 0 <= start <= end <= self.length(name):
            raise ValueError(f'projection interval outside {name!r}: {start}-{end}')
        runs, ends = self.segments(name)
        for index in range(bisect_right(ends, start), len(runs)):
            a, b, target, t, step = runs[index]
            if a >= end:
                break
            left, right = max(start, a), min(end, b)
            yield left, right, target, t + (left-a)*step, step

    def oriented(self, name, start, count, reverse=False):
        length = self.length(name)
        lo, hi = (length-start-count, length-start) if reverse else (start, start+count)
        pieces = list(self.project(name, lo, hi))
        if reverse:
            pieces.reverse()
        for a, b, destination, t, step in pieces:
            offset = hi-b if reverse else a-lo
            if reverse:
                t += (b-a-1)*step
                step = -step
            yield offset, b-a, destination, t, step

    def _build(self, name):
        if name in self.active:
            raise ValueError(f'circular insertion coordinate dependency at {name!r}')
        self.active.add(name)
        try:
            length = self.length(name)
            row = self.db.execute('SELECT info FROM alleles WHERE id=?', (name,)).fetchone()
            encoded = None
            if row:
                info = parse_info(row[0])
                for key in ('SEQ', 'SVINSSEQ', 'EXTENDGRAPHCIGAR', 'CIGAR'):
                    value = info.get(key)
                    if value not in (None, '', '.'):
                        encoded = value if value.startswith(('<', '>')) else None
                        break
                # Some files retain literal SEQ plus the partial alignment in
                # INFO/CIGAR. A self-target CIGAR is only the usual row label.
                if encoded is None:
                    for key in ('EXTENDGRAPHCIGAR', 'CIGAR'):
                        value = info.get(key)
                        if value and value.startswith(('<', '>')):
                            if any(target is not None and target != name
                                   for _, target, _ in cigar_chunks(value)):
                                encoded = value
                                break
            if encoded is None:
                return [[0, length, name, 0, 1]]
            runs, query = [], 0
            for direction, target, ops in cigar_chunks(encoded):
                cursor = 0
                for n, op, payload in ops:
                    if op in '=MXIS':
                        source_bases = self.reference(name, query, query+n)
                        if payload and payload != source_bases:
                            raise ValueError(f'coordinate CIGAR payload disagrees with insertion {name!r}')
                    if op in '=MX':
                        if target is None:
                            raise ValueError('encoded representative has an unnamed aligned segment')
                        if op == '=':
                            tlen = self.length(target)
                            lo, hi = (tlen-cursor-n, tlen-cursor) if direction == '<' else (cursor, cursor+n)
                            reference = self.reference(target, lo, hi)
                            if direction == '<':
                                reference = reference.translate(str.maketrans('ACGTN', 'TGCAN'))[::-1]
                            if reference != source_bases and any(a != b and a in 'ACGT' and b in 'ACGT'
                                                                 for a, b in zip(reference, source_bases)):
                                raise ValueError(f'coordinate = alignment disagrees with insertion {name!r}')
                        for offset, size, dest, t, step in self.oriented(target, cursor, n, direction == '<'):
                            append_run(runs, query+offset, query+offset+size, dest, t, step)
                        query += n
                        cursor += n
                    elif op in 'IS':
                        append_run(runs, query, query+n, name, query, 1)
                        query += n
                    elif op in 'DH':
                        cursor += n
            if query != length or not runs or runs[0][0] != 0 or runs[-1][1] != length:
                raise ValueError(f'coordinate map does not cover insertion {name!r}')
            if any(a[1] != b[0] for a, b in zip(runs, runs[1:])):
                raise ValueError(f'noncontiguous coordinate map for {name!r}')
            return runs
        finally:
            self.active.remove(name)
