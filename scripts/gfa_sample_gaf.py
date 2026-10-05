#!/usr/bin/env python3
"""Write one GAF per sample from its pseudo-linear map and merged-VCF alleles.

Inputs are the rGFA written by ``merged_vcf_to_gfa.py`` (with its
``OUTPUT.variants.tsv`` sidecar), the same merged VCFs in the same order, and
each sample's own VCF, whose ``##pseudoLinearMapping`` header lines map query
intervals onto reference paths.

Steps (RAM is bounded by the graph arrays and one sample's variant shard):

1. Load GFA segment lengths, P-line walks (CSR arrays), and links.
2. Load the sidecar as arrays indexed by VCF data-line number, and record each
   data line's byte offset as one uint64 array across all merged VCFs
   (``variant_offsets.npy``; file 2 offsets continue after file 1's end).
3. Stream the merged VCFs once. For every exported variant and every carrying
   sample observation, buffer (contig, query start/end, strand, variant). Every
   ``--shard-records`` records, append each sample's buffer to its shard file
   and close it again, so the number of open files never grows.
4. In parallel per sample: load its shard, sort by contig and query start,
   chain its pseudo-linear intervals into runs (query and reference both
   contiguous, as in tools/check_vcf_lossless.py), and splice each variant's
   allele walk into the run's reference path at its merged breakpoints
   (nested variants are spliced into their parent allele the same way; a
   deletion's allele is the rows nested in it). Consecutive runs of a contig
   are chained into one GAF record when they are contiguous on the query; a
   record breaks only where the walk cannot continue (a query gap, unmapped
   query, a run without graph support or whose walk does not have the
   query's length, a missing link, or a mid-node discontinuity).

With ``--add-links FILE``, a join between two steps that the GFA does not
link is kept (the two are adjacent in this sample's haplotype) and the link
is written to FILE as an L line; the final graph is the GFA plus FILE.

GAF columns 10/11 (matches/block) assume the walked graph sequence equals the
query; no base-level CIGAR is computed.
"""
import argparse
from array import array
from bisect import bisect_left
from collections import Counter, defaultdict
import gzip
import json
import multiprocessing as mp
import os
from pathlib import Path
import re
import sys
import tempfile

import numpy as np

KINDS = {'insertion': 0, 'snp': 1, 'substitution': 2, 'deletion': 3}
KIND_NAMES = {value: key for key, value in KINDS.items()}
# CHROM of rows nested in a merged row (graphvcfmerge row IDs) or on a shared
# full-locus-dup template (gfa_interval_metadata.TEMPLATE_PATH).
NESTED_CHROM = re.compile(r'(?:INS|DEL|SUB|DUP)_')
SHARD_DTYPE = np.dtype([('contig', '<u4'), ('qs', '<i8'), ('qe', '<i8'),
                        ('var', '<u8'), ('strand', 'u1')])
_META = re.compile(r'(\w+)="((?:[^"\\]|\\.)*)"')
_INTERVAL = re.compile(r'(.+):(\d+)-(\d+)([+-])')

# Graph state is loaded once in the parent and shared with fork workers.
G = {}


def _open(path, mode='rt'):
    return gzip.open(path, mode) if str(path).endswith('.gz') else open(path, mode)


def _log(message):
    print(f'[gfa-sample-gaf] {message}', file=sys.stderr, flush=True)


# ---------------------------------------------------------------- graph ----

def load_gfa(path):
    lengths = {}
    names, offsets, steps = [], array('q', [0]), array('q')
    edges = array('Q')
    with _open(path) as handle:
        for line in handle:
            kind = line[:1]
            if kind == 'S':
                fields = line.rstrip('\n').split('\t')
                length = next((int(tag[5:]) for tag in fields[3:] if tag.startswith('LN:i:')),
                              len(fields[2]))
                lengths[int(fields[1])] = length
            elif kind == 'L':
                fields = line.split('\t', 6)
                a = int(fields[1]) * 2 + (fields[2] == '-')
                b = int(fields[3]) * 2 + (fields[4] == '-')
                edges.append(_edge_key(a, b))
            elif kind == 'P':
                fields = line.split('\t', 3)
                names.append(fields[1])
                for token in fields[2].split(','):
                    steps.append(int(token[:-1]) * 2 + (token[-1] == '-'))
                offsets.append(len(steps))
    size = max(lengths, default=0) + 1
    segment_length = np.zeros(size, dtype=np.int64)
    for identifier, length in lengths.items():
        segment_length[identifier] = length
    G.update(
        segment_length=segment_length,
        path_names=names,
        path_index={name: index for index, name in enumerate(names)},
        offsets=np.frombuffer(offsets, dtype=np.int64),
        steps=np.frombuffer(steps, dtype=np.int64),
        edges=np.unique(np.frombuffer(edges, dtype=np.uint64)),
    )
    _log(f'GFA: {len(lengths)} segments, {len(names)} paths, {len(G["edges"])} links')


def _edge_key(a, b):
    """Canonical key of a directed link a->b (encoded id*2+reverse)."""
    ra, rb = b ^ 1, a ^ 1
    if (ra, rb) < (a, b):
        a, b = ra, rb
    return (a << 32) | b


def _has_edge(a, b):
    key = np.uint64(_edge_key(a, b))
    edges = G['edges']
    position = np.searchsorted(edges, key)
    return position < len(edges) and edges[position] == key


def load_variant_index(path):
    rows = []
    with open(path) as handle:
        for line in handle:
            if line.startswith('#'):
                continue
            rows.append(line.rstrip('\n').split('\t'))
    count = max((int(row[0]) for row in rows), default=-1) + 1
    path_of = np.full(count, -1, dtype=np.int64)
    parent = np.full(count, -1, dtype=np.int64)
    start = np.zeros(count, dtype=np.int64)
    end = np.zeros(count, dtype=np.int64)
    kind = np.full(count, -1, dtype=np.int8)
    reverse = np.zeros(count, dtype=np.int8)
    index = G['path_index']
    missing = 0
    for vcf_index, _identifier, name, variant_kind, parent_name, low, high, strand in rows:
        i = int(vcf_index)
        if i < 0 or parent_name == '.':
            continue
        if parent_name not in index or (variant_kind != 'deletion' and name not in index):
            missing += 1
            continue
        # A deletion's P path (its flanks) names it as the parent of rows
        # nested in it; its allele is those rows (Sample.allele).
        path_of[i] = index.get(name, -1)
        parent[i] = index[parent_name]
        start[i], end[i] = int(low), int(high)
        kind[i] = KINDS.get(variant_kind, 0)
        reverse[i] = strand == '-'
    if missing:
        _log(f'warning: {missing} indexed variants name a path absent from the GFA')
    G.update(var_path=path_of, var_parent=parent, var_start=start, var_end=end,
             var_kind=kind, var_reverse=reverse)
    _log(f'variant index: {int((parent >= 0).sum())} exported of {count} VCF records')


def load_local_paths(path):
    aliases = {}
    if path and Path(path).is_file():
        with open(path) as handle:
            next(handle, None)
            for line in handle:
                name, target, low, high, strand = line.rstrip('\n').split('\t')
                aliases[name] = (target, int(low), int(high), strand)
    G['aliases'] = aliases


# ----------------------------------------------------- merged VCF shards ----

def build_shards(vcfs, samples, shard_dir, shard_records):
    """Stream merged VCFs once into per-sample (contig, query, variant) shards."""
    contigs = {}
    buffers = {sample: [array('I'), array('q'), array('q'), array('Q'), array('B')]
               for sample in samples}
    buffered = 0
    offsets = array('Q')
    shard_paths = {sample: shard_dir / f'{index}.shard'
                   for index, sample in enumerate(samples)}
    for path in shard_paths.values():
        path.write_bytes(b'')

    def flush():
        nonlocal buffered
        for sample, columns in buffers.items():
            if not columns[0]:
                continue
            records = np.empty(len(columns[0]), dtype=SHARD_DTYPE)
            for name, column in zip(SHARD_DTYPE.names, columns):
                records[name] = np.frombuffer(column, dtype=records.dtype[name])
            with open(shard_paths[sample], 'ab') as out:
                out.write(records.tobytes())
            buffers[sample] = [array('I'), array('q'), array('q'), array('Q'), array('B')]
        buffered = 0

    exported = G['var_parent']
    base = 0
    variant = 0
    for vcf in vcfs:
        columns = {}
        position = 0
        with _open(vcf, 'rb') as handle:
            for raw in handle:
                start = position
                position += len(raw)
                if raw.startswith(b'#'):
                    if raw.startswith(b'#CHROM'):
                        header = raw.decode().rstrip('\r\n').split('\t')
                        columns = {index: name for index, name in enumerate(header)
                                   if index >= 9 and name in buffers}
                    continue
                if not raw.strip():
                    continue
                offsets.append(base + start)
                current = variant
                variant += 1
                if current >= len(exported) or exported[current] < 0 or not columns:
                    continue
                fields = raw.decode().rstrip('\r\n').split('\t')
                if len(fields) < 10:
                    continue
                keys = fields[8].split(':')
                try:
                    gt, ctg, qc = keys.index('GT'), keys.index('ASSEMBLYCONTIG'), keys.index('QUERYCOORD')
                except ValueError:
                    continue
                snp = G['var_kind'][current] == KINDS['snp']
                kind_key = keys.index('TYPE') if 'TYPE' in keys else None
                # A SNP row placed directly on its backbone path is top level:
                # its INS_SNP observations restate an insertion's bases against
                # that insertion's template and are skipped
                # (tools/check_vcf_lossless.py, check_merge_lossless.py). A row
                # on a merged row's path (an insertion or deletion row, or a
                # shared DUP_ template) is nested: its INS_SNPs are the copy's
                # own differences and are kept.
                top_level = (snp and G['path_names'][int(exported[current])] == fields[0]
                             and not NESTED_CHROM.match(fields[0]))
                for column, sample in columns.items():
                    text = fields[column]
                    if not text.startswith('1'):
                        continue
                    values = text.split(':')
                    if len(values) != len(keys) or values[gt] != '1':
                        continue
                    names, coords = values[ctg].split(','), values[qc].split(',')
                    kinds = values[kind_key].split(',') if kind_key is not None else []
                    for index, coord in enumerate(coords):
                        if top_level and len(kinds) == len(coords) and kinds[index] != 'SNP':
                            continue
                        match = re.fullmatch(r'(\d+)(?:-(\d+))?([+-]?)', coord)
                        if match is None:
                            continue
                        name = names[index if len(names) > 1 else 0]
                        low = int(match[1])
                        if snp and not match[2] and match[3] == '-':
                            low -= 1  # '-' points are traversal boundaries
                        high = int(match[2]) if match[2] else low + (1 if snp else 0)
                        low, high = min(low, high), max(low, high)
                        target = buffers[sample]
                        target[0].append(contigs.setdefault(name, len(contigs)))
                        target[1].append(low)
                        target[2].append(high)
                        target[3].append(current)
                        target[4].append(match[3] == '-')
                        buffered += 1
                if buffered >= shard_records:
                    flush()
        base += position
    flush()
    return shard_paths, [name for name, _index in sorted(contigs.items(), key=lambda item: item[1])], \
        np.frombuffer(offsets, dtype=np.uint64)


# ------------------------------------------------------------- walking ----

class Walk:
    """One GAF record under construction: encoded steps plus edge trims."""

    def __init__(self, qstart):
        self.qstart = qstart
        self.qend = qstart
        self.steps = []
        self.start_offset = 0
        self.end_trim = 0
        self.variants = 0
        self.matches = 0
        self.pending = []      # links added by the extend() in progress

    def append(self, piece, stats):
        """Add (steps, start_offset, end_trim); False if it cannot continue."""
        steps, offset, trim = piece
        if not steps:
            return True
        if not self.steps:
            self.steps = list(steps)
            self.start_offset, self.end_trim = offset, trim
            return True
        last, first = self.steps[-1], steps[0]
        length = int(G['segment_length'][last >> 1])
        if last == first and length - self.end_trim == offset:
            # Continue inside the same node (adjacent pseudo-linear pieces).
            self.steps.extend(steps[1:])
            self.end_trim = trim
            return True
        if self.end_trim or offset:
            stats['break_mid_node'] += 1
            return False
        if not _has_edge(last, first):
            added = G.get('added_links')
            if added is None:
                stats['break_missing_link'] += 1
                return False
            if _edge_key(last, first) not in added:
                added.add(_edge_key(last, first))
                self.pending.append(_edge_key(last, first))
                stats['link_added'] += 1
        self.steps.extend(steps)
        self.end_trim = trim
        return True

    def extend(self, pieces, stats):
        """Append all pieces or none (rolls back on a failed join)."""
        state = len(self.steps), self.start_offset, self.end_trim
        self.pending = []
        for piece in pieces:
            if not self.append(piece, stats):
                del self.steps[state[0]:]
                _count, self.start_offset, self.end_trim = state
                # Links added for this failed join are not used by any walk.
                for key in self.pending:
                    G['added_links'].discard(key)
                stats['link_added'] -= len(self.pending)
                self.pending = []
                return False
        return True


def _path_cumulative(path, cache):
    cumulative = cache.get(path)
    if cumulative is None:
        steps = G['steps'][G['offsets'][path]:G['offsets'][path + 1]]
        cumulative = np.zeros(len(steps) + 1, dtype=np.int64)
        np.cumsum(G['segment_length'][steps >> 1], out=cumulative[1:])
        cache[path] = cumulative
    return cumulative


def span(path, low, high, cache):
    """Forward piece covering [low, high) of a P path."""
    if high <= low:
        return None
    cumulative = _path_cumulative(path, cache)
    if low < 0 or high > cumulative[-1]:
        raise ValueError(f'{G["path_names"][path]}: interval {low}-{high} outside 0-{cumulative[-1]}')
    first = int(np.searchsorted(cumulative, low, side='right')) - 1
    last = int(np.searchsorted(cumulative, high, side='left')) - 1
    base = G['offsets'][path]
    steps = G['steps'][base + first:base + last + 1].tolist()
    return steps, int(low - cumulative[first]), int(cumulative[last + 1] - high)


def _edge_owned(sample, nested, low, high, descending, length, path):
    """Nested candidates of one allele, without the edge points it does not
    own (allele_edges.py)."""
    if length is None or low == high:
        return nested
    from allele_edges import edge_owned
    first = high if descending else low
    kept, starts, ends, spans = [], [], [], []
    current = length
    for candidate in nested:
        point = int(sample.qs[candidate])
        variant = int(sample.var[candidate])
        child = G['var_parent'][variant] == path
        a, b = int(G['var_start'][variant]), int(G['var_end'][variant])
        if child and point == int(sample.qe[candidate]) and point in (low, high):
            (starts if point == first else ends).append((a, b, candidate))
            continue
        kept.append(candidate)
        if child:
            spans.append((a, b))
            current += int(sample.qe[candidate]) - point - (b - a)
    return kept + edge_owned(starts, ends, spans, current, high - low, length)


def reverse_pieces(pieces):
    return [([step ^ 1 for step in reversed(steps)], trim, offset)
            for steps, offset, trim in reversed(pieces)]


class Sample:
    """Per-sample variant table sorted by contig and query start."""

    def __init__(self, records):
        order = np.lexsort((records['qe'], records['qs'], records['contig']))
        self.records = records[order]
        # Contiguous int64 columns: searchsorted on a structured-array field
        # view, or on uint32 with a Python int key, copies the whole column on
        # every call (quadratic over a sample).
        for name in ('contig', 'qs', 'qe', 'var', 'strand'):
            setattr(self, name, np.ascontiguousarray(self.records[name], dtype=np.int64))
        self.used = np.zeros(len(self.records), dtype=bool)
        self.skipped = set()     # positions skipped for overlapping another variant
        self.cache = {}

    def snapshot(self, contig, low, high):
        """The used flags of records starting in [low, high] (restore())."""
        if contig is None:
            return None
        left = int(np.searchsorted(self.contig, contig, side='left'))
        right = int(np.searchsorted(self.contig, contig, side='right'))
        first = left + int(np.searchsorted(self.qs[left:right], low, side='left'))
        last = left + int(np.searchsorted(self.qs[left:right], high, side='right'))
        return first, self.used[first:last].copy()

    def restore(self, saved):
        if saved is not None:
            first, values = saved
            self.used[first:first + len(values)] = values

    def between(self, contig, low, high):
        """Record positions with query span inside [low, high] on contig."""
        left = int(np.searchsorted(self.contig, contig, side='left'))
        right = int(np.searchsorted(self.contig, contig, side='right'))
        start = left + int(np.searchsorted(self.qs[left:right], low, side='left'))
        output = []
        qs, qe, used = self.qs, self.qe, self.used
        for position in range(start, right):
            if qs[position] > high:
                break
            if qe[position] <= high and not used[position]:
                output.append(position)
        return output

    def splice(self, path, low, high, candidates, stats, descending=False):
        """Pieces for [low, high) of path with candidate variants spliced.

        Variants at one breakpoint follow the query; ``descending`` when the
        path runs against the query."""
        chosen = []
        for position in candidates:
            variant = int(self.var[position])
            if G['var_parent'][variant] != path:
                continue
            a, b = int(G['var_start'][variant]), int(G['var_end'][variant])
            if low <= a <= b <= high:
                chosen.append((a, b, -position if descending else position, variant))
        chosen.sort()
        # A SNP inside one of this sample's SVs here (a deletion or
        # substitution covering it) is redundant: the SV replaces that base.
        # Placing it first would make the SV look overlapping and drop it
        # (tools/check_vcf_lossless.py skips such SNPs the same way).
        sv_starts, sv_reach = [], []
        for a, b, _position, variant in chosen:
            if b > a and G['var_kind'][variant] != KINDS['snp']:
                sv_starts.append(a)
                sv_reach.append(max(b, sv_reach[-1] if sv_reach else b))
        pieces = []
        cursor = low
        for a, b, position, variant in chosen:
            position = abs(position)
            if G['var_kind'][variant] == KINDS['snp'] and sv_starts:
                covering = bisect_left(sv_starts, a + 1) - 1
                if covering >= 0 and sv_reach[covering] >= b:
                    self.used[position] = True
                    stats['snp_inside_sv'] += 1
                    continue
            if a < cursor:
                stats['variant_overlap_skipped'] += 1
                self.skipped.add(position)
                continue
            pieces.append(span(path, cursor, a, self.cache))
            pieces.extend(self.allele(position, variant, stats))
            self.used[position] = True
            stats['variants_spliced'] += 1
            cursor = b
        pieces.append(span(path, cursor, high, self.cache))
        return [piece for piece in pieces if piece]

    def allele(self, position, variant, stats):
        """Allele walk in its parent's orientation, with nested variants."""
        path = int(G['var_path'][variant])
        low, high = int(self.qs[position]), int(self.qe[position])
        # Nested rows follow the query along the parent path: backwards when
        # this allele is carried on the query's reverse strand.
        descending = bool(self.strand[position])
        length = (int(_path_cumulative(path, self.cache)[-1])
                  if path >= 0 and G['var_kind'][variant] != KINDS['deletion'] else None)
        nested = _edge_owned(self, [candidate for candidate in self.between(
            int(self.contig[position]), low, high) if candidate != position],
            low, high, descending, length, path)
        self.used[position] = True
        if G['var_kind'][variant] == KINDS['deletion']:
            # A deletion has no bases; rows nested in it (--exact: a member's
            # kept bases) are its allele, by offset and then query order.
            if path < 0:
                return []
            chosen = sorted(
                (int(G['var_start'][int(self.var[candidate])]),
                 -candidate if descending else candidate)
                for candidate in nested if G['var_parent'][int(self.var[candidate])] == path)
            pieces = []
            for _offset, candidate in chosen:
                candidate = abs(candidate)
                if self.used[candidate]:
                    continue
                pieces.extend(self.allele(candidate, int(self.var[candidate]), stats))
                stats['variants_spliced'] += 1
            return reverse_pieces(pieces) if G['var_reverse'][variant] else pieces
        cumulative = _path_cumulative(path, self.cache)
        pieces = self.splice(path, 0, int(cumulative[-1]), nested, stats, descending)
        return reverse_pieces(pieces) if G['var_reverse'][variant] else pieces


def read_pseudolinear(vcf, only=None):
    """``only``: keep that Sample's lines (a combined mapping header file)."""
    intervals = []
    seen = set()
    sample = None
    with _open(vcf) as handle:
        for line in handle:
            if line.startswith('#CHROM'):
                columns = line.rstrip('\n').split('\t')
                sample = columns[9] if len(columns) > 9 else None
                break
            if not line.startswith('##pseudoLinearMapping=<'):
                continue
            values = {key: bytes(value, 'utf-8').decode('unicode_escape')
                      for key, value in _META.findall(line)}
            if only is not None and values.get('Sample') != only:
                continue
            if values.get('Path') == 'alt':
                continue          # sequence source: places no bases
            query = _INTERVAL.fullmatch(values.get('Query', ''))
            if query is None:
                continue
            key = (values['Query'], values.get('Category', '.'))
            if query[2] == query[3]:
                # A deletion (a query point) claims no query bases: several
                # at one point are distinct reference spans, all kept (as in
                # tools/check_vcf_lossless.py).
                key += (values.get('Reference'),)
            if key in seen:
                # Locus duplications repeat the query with their source
                # interval; keep the first (insertion-point) line.
                continue
            seen.add(key)
            reference = _INTERVAL.fullmatch(values.get('Reference', ''))
            intervals.append((query[1], int(query[2]), int(query[3]), query[4],
                              values.get('Category', '.'),
                              (reference[1], int(reference[2]), int(reference[3]), reference[4])
                              if reference else None))
    # Ties at one query point (deletions) in walking order along the
    # reference, so the one adjacent to the run extends it first.
    intervals.sort(key=lambda value: (value[0], value[1], value[2], 0 if value[5] is None else (
        value[5][1] if value[3] == value[5][3] else -value[5][2])))
    return sample, intervals


def _canonical_reference(reference):
    """Translate a local reference name onto its exported P path."""
    name, low, high, strand = reference
    index = G['path_index']
    if name in index:
        return index[name], low, high, strand
    alias = G['aliases'].get(name)
    if alias is None or alias[0] not in index:
        return None
    target, start, end, alias_strand = alias
    if alias_strand == '+':
        low, high = start + low, start + high
    else:
        low, high = end - high, end - low
    return index[target], low, high, strand if alias_strand == '+' else ('-' if strand == '+' else '+')


def _runs(intervals, stats):
    """Chain contiguous pseudo-linear intervals on one graph path into runs.

    Yields (contig, qs, qe, qstrand, (path, low, high, rstrand) or None,
    category, member intervals as (qs, qe, canonical, category)). As in tools/check_vcf_lossless.py, a run continues while the
    query and the reference are both contiguous on the same strands, so a
    variant is placed wherever its merged breakpoints fall inside the run.
    Intervals without a graph anchor come out on their own (run None).
    """
    current = None
    for contig, qs, qe, qstrand, category, reference in intervals:
        stats['intervals'] += 1
        canonical = None
        if category == 'UNMAPPED':
            stats['break_unmapped'] += 1
        elif reference is not None:
            canonical = _canonical_reference(reference)
            if canonical is None:
                stats['break_reference_not_in_gfa'] += 1
        elif category not in ('INSERTION', 'DELETION'):
            stats['break_reference_not_in_gfa'] += 1
        member = (qs, qe, canonical, category)
        if canonical is not None and current is not None:
            path, low, high, rstrand = current[4]
            other, start, end, strand = canonical
            if (contig, qstrand, other, strand) == (current[0], current[3], path, rstrand) \
                    and qs == current[2]:
                if qstrand == strand and start == high:
                    current = (contig, current[1], qe, qstrand, (path, low, end, strand), 'RUN',
                               current[6])
                    current[6].append(member)
                    continue
                if qstrand != strand and end == low:
                    current = (contig, current[1], qe, qstrand, (path, start, high, strand), 'RUN',
                               current[6])
                    current[6].append(member)
                    continue
        if current is not None:
            stats['runs'] += 1
            yield current
            current = None
        if canonical is None:
            yield contig, qs, qe, qstrand, None, category, [member]
        else:
            current = (contig, qs, qe, qstrand, canonical, category, [member])
    if current is not None:
        stats['runs'] += 1
        yield current


def _walked_length(pieces):
    return sum(int(G['segment_length'][np.asarray(steps) >> 1].sum()) - offset - trim
               for steps, offset, trim in pieces)


def write_sample_gaf(task):
    sample_name, vcf, shard, contig_names, lengths, output, add_links, only = task
    # Links this sample's walks need that the GFA lacks (--add-links).
    G['added_links'] = set() if add_links else None
    records = np.fromfile(shard, dtype=SHARD_DTYPE)
    sample = Sample(records)
    contig_ids = {name: index for index, name in enumerate(contig_names)}
    _vcf_sample, intervals = read_pseudolinear(vcf, only)
    stats = Counter()
    written = 0
    breaks = open(output + '.breaks.tsv.tmp', 'w')
    debug = open(output + '.lengthdebug.tsv.tmp', 'w')
    debug.write('contig\tquery_start\tquery_end\tcategory\trun_path\twalk_minus_query\t'
                'candidates\tvcf_index:kind:start-end:query:state\n')
    breaks.write('contig\tquery_start\tquery_end\tbases\tcategory\toutcome\treason\n')

    def report(contig, qs, qe, category, outcome, before):
        """One line per split run or dropped interval: the break counters it
        raised (the reason)."""
        reason = ','.join(sorted(key for key in stats
                                 if key.startswith('break_') and stats[key] > before.get(key, 0)))
        breaks.write(f'{contig}\t{qs}\t{qe}\t{qe - qs}\t{category}\t{outcome}\t{reason or "."}\n')

    with open(output + '.tmp', 'w') as out:
        walk = None
        current_contig = None

        def close():
            nonlocal walk, written
            if walk is not None and walk.steps and walk.qend > walk.qstart:
                node_lengths = G['segment_length'][np.asarray(walk.steps) >> 1]
                path_length = int(node_lengths.sum())
                path_start = walk.start_offset
                path_end = path_length - walk.end_trim
                query_span = walk.qend - walk.qstart
                block = max(query_span, path_end - path_start)
                matches = min(query_span, path_end - path_start)
                path = ''.join(('<' if step & 1 else '>') + str(step >> 1) for step in walk.steps)
                length = lengths.get(current_contig, walk.qend)
                out.write(f'{current_contig}\t{length}\t{walk.qstart}\t{walk.qend}\t+\t'
                          f'{path}\t{path_length}\t{path_start}\t{path_end}\t{matches}\t'
                          f'{block}\t255\tnv:i:{walk.variants}\n')
                written += 1
            walk = None

        def place(qs, qe, qstrand, run, category, contig_id):
            """The walk pieces of one run (or anchorless interval), or None."""
            candidates = sample.between(contig_id, qs, qe) if contig_id is not None else []
            pieces = None
            if run is not None:
                # Contiguous mappings: the reference walk with the sample's
                # variants spliced at their merged breakpoints, which may lie
                # anywhere in the run (--exact moves shifted breakpoints).
                path, low, high, rstrand = run
                pieces = sample.splice(path, low, high, candidates, stats, rstrand != qstrand)
                if rstrand != qstrand:
                    pieces = reverse_pieces(pieces)
            elif category in ('INSERTION', 'DELETION'):
                # No usable anchor (e.g. a locus duplication): use the
                # variants carried inside this query interval directly.
                pieces = []
                variants = [int(sample.var[position]) for position in candidates]
                inner = {int(G['var_path'][variant]) for variant in variants} - {-1}
                for position, variant in zip(candidates, variants):
                    # Nested variants are spliced inside their parent allele.
                    if sample.used[position] or int(G['var_parent'][variant]) in inner:
                        continue
                    # A deletion adds no bases here; it belongs to the
                    # anchored run next to this interval.
                    if G['var_kind'][variant] == KINDS['deletion']:
                        continue
                    allele = sample.allele(position, variant, stats)
                    if sample.strand[position]:
                        allele = reverse_pieces(allele)
                    pieces.extend(allele)
                    stats['variants_spliced'] += 1
                if qe > qs and not pieces:
                    stats['break_insertion_without_variant'] += 1
                    pieces = None
            if pieces is not None and _walked_length(pieces) != qe - qs:
                # A variant was not placed here: never write a walk that does
                # not spell the query.
                stats['break_length_mismatch'] += 1
                length_report(qs, qe, category, run, candidates, _walked_length(pieces) - (qe - qs))
                pieces = None
            return pieces

        def length_report(qs, qe, category, run, candidates, delta):
            """Every candidate variant of a walk whose length is off: placed
            or not, and why not (another path, outside the interval, overlap)."""
            path = run[0] if run is not None else -1
            low, high = (run[1], run[2]) if run is not None else (0, 0)
            items = []
            for position in candidates:
                variant = int(sample.var[position])
                parent = int(G['var_parent'][variant])
                a, b = int(G['var_start'][variant]), int(G['var_end'][variant])
                if sample.used[position]:
                    state = 'placed'
                elif position in sample.skipped:
                    state = 'overlap_skipped'
                elif parent != path:
                    state = 'other_path:' + (G['path_names'][parent] if parent >= 0 else '.')
                elif not low <= a <= b <= high:
                    state = 'outside'
                else:
                    state = 'unused'
                items.append(f'{variant}:{KIND_NAMES.get(int(G["var_kind"][variant]), "?")}:'
                             f'{a}-{b}:q{int(sample.qs[position])}-{int(sample.qe[position])}:{state}')
            debug.write(f'{current_contig}\t{qs}\t{qe}\t{category}\t'
                        f'{G["path_names"][path] if path >= 0 else "."}:{low}-{high}\t{delta}\t'
                        f'{len(items)}\t{";".join(items)}\n')

        def walk_on(qs, qe, pieces, spliced):
            """Add pieces to the record, or start a new one; False (nothing
            changed) when the pieces are not a walk of the graph by themselves."""
            nonlocal walk
            # A record continues only where the query is contiguous.
            if walk is None or qs != walk.qend or not walk.extend(pieces, stats):
                fresh = Walk(qs)
                if not fresh.extend(pieces, stats):
                    return False
                close()
                walk = fresh
            walk.qend = qe
            walk.variants += spliced
            stats['runs_walked'] += 1
            return True

        for contig, qs, qe, qstrand, run, category, members in _runs(intervals, stats):
            if contig != current_contig:
                close()
                current_contig = contig
            contig_id = contig_ids.get(contig)
            saved = Counter(stats), sample.snapshot(contig_id, qs, qe)
            spliced_start = stats['variants_spliced']
            pieces = place(qs, qe, qstrand, run, category, contig_id)
            if pieces is not None and walk_on(qs, qe, pieces,
                                              stats['variants_spliced'] - spliced_start):
                continue
            report(contig, qs, qe, category, 'run_failed', saved[0])
            # Variants consumed by the failed attempt stay available to the
            # neighbouring runs.
            sample.restore(saved[1])
            if len(members) > 1:
                # Retry this run interval by interval (the pre-run granularity),
                # so only intervals that cannot be walked are lost.
                stats.clear()
                stats.update(saved[0])
                sample.restore(saved[1])
                stats['run_split'] += 1
                for low, high, canonical, member_category in members:
                    spliced_start = stats['variants_spliced']
                    before = Counter(stats)
                    snapshot = sample.snapshot(contig_id, low, high)
                    pieces = place(low, high, qstrand, canonical, member_category, contig_id)
                    if pieces is None or not walk_on(
                            low, high, pieces, stats['variants_spliced'] - spliced_start):
                        sample.restore(snapshot)
                        stats['run_dropped'] += pieces is not None
                        report(contig, low, high, member_category,
                               'interval_dropped' if pieces is not None else 'interval_unwalked',
                               before)
                        close()
                continue
            stats['run_dropped'] += pieces is not None
            close()
        close()
    breaks.close()
    debug.close()
    os.replace(output + '.breaks.tsv.tmp', output + '.breaks.tsv')
    os.replace(output + '.lengthdebug.tsv.tmp', output + '.lengthdebug.tsv')
    os.replace(output + '.tmp', output)
    stats['variants_unplaced'] = int((~sample.used).sum())
    stats['gaf_records'] = written
    # The sample's variants no walk placed: VCF index (variants.tsv column 1),
    # query interval, graph placement, and whether one was skipped for
    # overlapping another of the sample's variants.
    with open(output + '.unplaced.tsv.tmp', 'w') as out:
        out.write('contig\tquery_start\tquery_end\tstrand\tvcf_index\tparent\t'
                  'parent_start\tparent_end\treason\n')
        for position in np.flatnonzero(~sample.used).tolist():
            variant = int(sample.var[position])
            parent = int(G['var_parent'][variant]) if variant < len(G['var_parent']) else -1
            out.write(f'{contig_names[int(sample.contig[position])]}\t{int(sample.qs[position])}\t'
                      f'{int(sample.qe[position])}\t{"-" if sample.strand[position] else "+"}\t'
                      f'{variant}\t{G["path_names"][parent] if parent >= 0 else "."}\t'
                      f'{int(G["var_start"][variant])}\t{int(G["var_end"][variant])}\t'
                      f'{"overlap_skipped" if position in sample.skipped else "unused"}\n')
    os.replace(output + '.unplaced.tsv.tmp', output + '.unplaced.tsv')
    return sample_name, dict(stats), sorted(G['added_links'] or ())


# ----------------------------------------------------------------- main ----

def _fai_lengths(query_list, samples):
    lengths = {}
    if not query_list:
        return lengths
    base = Path(query_list).resolve().parent
    with open(query_list) as handle:
        for line in handle:
            fields = line.split()
            if len(fields) < 2 or fields[0] not in samples:
                continue
            fasta = Path(fields[1]) if Path(fields[1]).is_absolute() else base / fields[1]
            fai = Path(fields[2]) if len(fields) > 2 else Path(str(fasta) + '.fai')
            if not fai.is_absolute():
                fai = base / fai
            if fai.is_file():
                table = lengths.setdefault(fields[0], {})
                with open(fai) as index:
                    for row in index:
                        name, length = row.split('\t')[:2]
                        table[name] = int(length)
    return lengths


def build_parser():
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('-g', '--gfa', required=True, help='rGFA from merged_vcf_to_gfa.py')
    parser.add_argument('--variant-index', help='default: GFA.variants.tsv')
    parser.add_argument('--local-paths', help='default: GFA.anchors.bed.local-paths.tsv if present')
    parser.add_argument('-v', '--vcf', required=True, action='append', nargs='+',
                        help='merged VCFs, in the same order given to merged_vcf_to_gfa.py')
    parser.add_argument('-s', '--sample-vcf', required=True, action='append', nargs='+',
                        help='per-sample VCFs carrying ##pseudoLinearMapping headers, or one '
                             'merged header of all samples\' mapping lines (tools/rebuild_assemblies.py '
                             'headers, or a VCF header whose #CHROM line lists several samples): '
                             'its samples are the lines\' Sample values; no per-sample VCF is read')
    parser.add_argument('-q', '--query-fasta-list',
                        help='NAME FASTA [FAI] per line, for GAF query lengths')
    parser.add_argument('-o', '--output-folder', required=True, help='writes SAMPLE.gaf here')
    parser.add_argument('-t', '--processes', type=int, default=1)
    parser.add_argument('--shard-records', type=int, default=100_000,
                        help='buffered records before flushing shards (default: 100000)')
    parser.add_argument('--tmpdir', help='shard folder parent (default: output folder)')
    parser.add_argument('--add-links', metavar='FILE',
                        help='keep joins the GFA does not link (adjacent in a sample) and write '
                             'those links to FILE as L lines; the final graph is the GFA plus FILE')
    return parser


def main(argv=None):
    args = build_parser().parse_args(argv)
    vcfs = [path for group in args.vcf for path in group]
    sample_vcfs = [path for group in args.sample_vcf for path in group]
    if args.processes < 1 or args.shard_records < 1:
        raise SystemExit('--processes and --shard-records must be positive')
    output = Path(args.output_folder)
    output.mkdir(parents=True, exist_ok=True)

    load_gfa(args.gfa)
    load_variant_index(args.variant_index or args.gfa + '.variants.tsv')
    load_local_paths(args.local_paths or args.gfa + '.anchors.bed.local-paths.tsv')

    samples = {}        # name -> (file, Sample filter or None)
    for vcf in sample_vcfs:
        name, listed = None, []
        with _open(vcf) as handle:
            for line in handle:
                if line.startswith('#CHROM'):
                    columns = line.rstrip('\n').split('\t')
                    name = columns[9] if len(columns) == 10 else None
                    break
                if line.startswith('##pseudoLinearMapping=<'):
                    value = dict(_META.findall(line)).get('Sample')
                    if value and value not in listed:
                        listed.append(value)
        # A per-sample VCF names its one sample column.  A merged header
        # (no #CHROM line, or a #CHROM line with several samples) holds every
        # sample's lines, each named by its Sample value, so the original
        # per-sample VCFs are not needed.
        entries = [(name, None)] if name else [(value, value) for value in listed]
        if not entries:
            raise SystemExit(f'{vcf}: no sample column or ##pseudoLinearMapping samples')
        for value, only in entries:
            if value in samples:
                raise SystemExit(f'sample {value!r} is given twice')
            samples[value] = (vcf, only)
    lengths = _fai_lengths(args.query_fasta_list, samples)

    with tempfile.TemporaryDirectory(prefix='gaf-shards-', dir=args.tmpdir or output) as directory:
        _log(f'streaming {len(vcfs)} merged VCF(s) into shards for {len(samples)} sample(s)')
        shards, contig_names, offsets = build_shards(
            vcfs, list(samples), Path(directory), args.shard_records)
        np.save(output / 'variant_offsets.npy', offsets)
        _log(f'{len(offsets)} VCF records indexed; shards in {directory}')
        tasks = [(name, vcf, str(shards[name]), contig_names, lengths.get(name, {}),
                  str(output / f'{name}.gaf'), bool(args.add_links), only)
                 for name, (vcf, only) in samples.items()]
        if args.processes > 1 and len(tasks) > 1:
            context = mp.get_context('fork')
            with context.Pool(min(args.processes, len(tasks))) as pool:
                results = list(pool.imap_unordered(write_sample_gaf, tasks))
        else:
            results = [write_sample_gaf(task) for task in tasks]
    if args.add_links:
        links = sorted({key for _name, _stats, added in results for key in added})
        with open(args.add_links + '.tmp', 'w') as out:
            for key in links:
                a, b = key >> 32, key & 0xFFFFFFFF
                out.write(f'L\t{a >> 1}\t{"-" if a & 1 else "+"}\t'
                          f'{b >> 1}\t{"-" if b & 1 else "+"}\t0M\n')
        os.replace(args.add_links + '.tmp', args.add_links)
        _log(f'{len(links)} link(s) added by sample walks -> {args.add_links}')
    results = [(name, stats) for name, stats, _added in results]
    with open(output / 'gaf_stats.tsv', 'w') as out:
        keys = sorted({key for _name, stats in results for key in stats})
        out.write('sample\t' + '\t'.join(keys) + '\n')
        for name, stats in sorted(results):
            out.write(name + '\t' + '\t'.join(str(stats.get(key, 0)) for key in keys) + '\n')
            _log(f'{name}: ' + ', '.join(f'{key}={stats.get(key, 0)}' for key in keys))
    return 0


if __name__ == '__main__':
    raise SystemExit(main())
