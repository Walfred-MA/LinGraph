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
   walk its pseudo-linear intervals in query order, and splice each variant's
   allele walk into the reference path at its breakpoints (nested variants
   are spliced into their parent insertion the same way). Consecutive
   intervals of a contig are chained into one GAF record; a record breaks
   only where the walk cannot continue (unmapped query, an interval without
   graph support, a missing link, or a mid-node discontinuity).

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
        path_of[i] = index[name] if variant_kind != 'deletion' else -1
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
                for column, sample in columns.items():
                    text = fields[column]
                    if not text.startswith('1'):
                        continue
                    values = text.split(':')
                    if len(values) != len(keys) or values[gt] != '1':
                        continue
                    names, coords = values[ctg].split(','), values[qc].split(',')
                    for index, coord in enumerate(coords):
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
            stats['break_missing_link'] += 1
            return False
        self.steps.extend(steps)
        self.end_trim = trim
        return True

    def extend(self, pieces, stats):
        """Append all pieces or none (rolls back on a failed join)."""
        state = len(self.steps), self.start_offset, self.end_trim
        for piece in pieces:
            if not self.append(piece, stats):
                del self.steps[state[0]:]
                _count, self.start_offset, self.end_trim = state
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


def reverse_pieces(pieces):
    return [([step ^ 1 for step in reversed(steps)], trim, offset)
            for steps, offset, trim in reversed(pieces)]


class Sample:
    """Per-sample variant table sorted by contig and query start."""

    def __init__(self, records):
        order = np.lexsort((records['qe'], records['qs'], records['contig']))
        self.records = records[order]
        self.used = np.zeros(len(self.records), dtype=bool)
        self.cache = {}

    def between(self, contig, low, high):
        """Record positions with query span inside [low, high] on contig."""
        records = self.records
        left = int(np.searchsorted(records['contig'], contig, side='left'))
        right = int(np.searchsorted(records['contig'], contig, side='right'))
        start = left + int(np.searchsorted(records['qs'][left:right], low, side='left'))
        output = []
        for position in range(start, right):
            if records['qs'][position] > high:
                break
            if records['qe'][position] <= high and not self.used[position]:
                output.append(position)
        return output

    def splice(self, path, low, high, candidates, stats):
        """Pieces for [low, high) of path with candidate variants spliced."""
        chosen = []
        for position in candidates:
            variant = int(self.records['var'][position])
            if G['var_parent'][variant] != path:
                continue
            a, b = int(G['var_start'][variant]), int(G['var_end'][variant])
            if low <= a <= b <= high:
                chosen.append((a, b, position, variant))
        chosen.sort()
        pieces = []
        cursor = low
        for a, b, position, variant in chosen:
            if a < cursor:
                stats['variant_overlap_skipped'] += 1
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
        if G['var_kind'][variant] == KINDS['deletion']:
            return []
        path = int(G['var_path'][variant])
        cumulative = _path_cumulative(path, self.cache)
        records = self.records
        nested = [candidate for candidate in self.between(
            int(records['contig'][position]), int(records['qs'][position]),
            int(records['qe'][position])) if candidate != position]
        self.used[position] = True
        pieces = self.splice(path, 0, int(cumulative[-1]), nested, stats)
        return reverse_pieces(pieces) if G['var_reverse'][variant] else pieces


def read_pseudolinear(vcf):
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
            query = _INTERVAL.fullmatch(values.get('Query', ''))
            if query is None:
                continue
            key = (values['Query'], values.get('Category', '.'))
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
    intervals.sort(key=lambda value: (value[0], value[1], value[2]))
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


def write_sample_gaf(task):
    sample_name, vcf, shard, contig_names, lengths, output = task
    records = np.fromfile(shard, dtype=SHARD_DTYPE)
    sample = Sample(records)
    contig_ids = {name: index for index, name in enumerate(contig_names)}
    _vcf_sample, intervals = read_pseudolinear(vcf)
    stats = Counter()
    written = 0
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

        for contig, qs, qe, qstrand, category, reference in intervals:
            if contig != current_contig:
                close()
                current_contig = contig
            stats['intervals'] += 1
            spliced_start = stats['variants_spliced']
            contig_id = contig_ids.get(contig)
            candidates = sample.between(contig_id, qs, qe) if contig_id is not None else []
            pieces = None
            if category == 'UNMAPPED':
                stats['break_unmapped'] += 1
            elif category == 'DELETION' and reference is not None and (
                    canonical := _canonical_reference(reference)):
                # Deleted reference bases are skipped, never walked. The next
                # piece must then be reachable by the deletion's graph link.
                path, low, high, _strand = canonical
                pieces = []
                for position in candidates:
                    variant = int(sample.records['var'][position])
                    if (G['var_kind'][variant] == KINDS['deletion'] and G['var_parent'][variant] == path
                            and low <= G['var_start'][variant] <= G['var_end'][variant] <= high):
                        sample.used[position] = True
                        stats['variants_spliced'] += 1
            elif reference is not None and (canonical := _canonical_reference(reference)):
                path, low, high, rstrand = canonical
                spliced_before = stats['variants_spliced']
                pieces = sample.splice(path, low, high, candidates, stats)
                if rstrand != qstrand:
                    pieces = reverse_pieces(pieces)
                if category == 'INSERTION' and qe > qs and stats['variants_spliced'] == spliced_before:
                    stats['break_insertion_without_variant'] += 1
                    pieces = None
            elif category in ('INSERTION', 'DELETION'):
                # No usable anchor (e.g. a locus duplication): use the
                # variants carried inside this query interval directly.
                pieces = []
                records = sample.records
                variants = [int(records['var'][position]) for position in candidates]
                inner = {int(G['var_path'][variant]) for variant in variants} - {-1}
                for position, variant in zip(candidates, variants):
                    # Nested variants are spliced inside their parent allele.
                    if sample.used[position] or int(G['var_parent'][variant]) in inner:
                        continue
                    allele = sample.allele(position, variant, stats)
                    if records['strand'][position]:
                        allele = reverse_pieces(allele)
                    pieces.extend(allele)
                    stats['variants_spliced'] += 1
                if qe > qs and not pieces:
                    stats['break_insertion_without_variant'] += 1
                    pieces = None
            else:
                stats['break_reference_not_in_gfa'] += 1
            if pieces is None:
                close()
                continue
            spliced = stats['variants_spliced'] - spliced_start
            if walk is None:
                walk = Walk(qs)
            if not walk.extend(pieces, stats):
                # Cannot chain onto the previous interval: start a new record.
                close()
                walk = Walk(qs)
                if not walk.extend(pieces, stats):
                    stats['interval_dropped'] += 1
                    walk = None
                    continue
            walk.qend = qe
            walk.variants += spliced
            stats['intervals_walked'] += 1
        close()
    os.replace(output + '.tmp', output)
    stats['variants_unplaced'] = int((~sample.used).sum())
    stats['gaf_records'] = written
    return sample_name, dict(stats)


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
                        help='per-sample VCFs carrying ##pseudoLinearMapping headers')
    parser.add_argument('-q', '--query-fasta-list',
                        help='NAME FASTA [FAI] per line, for GAF query lengths')
    parser.add_argument('-o', '--output-folder', required=True, help='writes SAMPLE.gaf here')
    parser.add_argument('-t', '--processes', type=int, default=1)
    parser.add_argument('--shard-records', type=int, default=100_000,
                        help='buffered records before flushing shards (default: 100000)')
    parser.add_argument('--tmpdir', help='shard folder parent (default: output folder)')
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

    samples = {}
    for vcf in sample_vcfs:
        name, _intervals = None, None
        with _open(vcf) as handle:
            for line in handle:
                if line.startswith('#CHROM'):
                    columns = line.rstrip('\n').split('\t')
                    name = columns[9] if len(columns) > 9 else None
                    break
        if not name:
            raise SystemExit(f'{vcf}: no sample column')
        if name in samples:
            raise SystemExit(f'sample {name!r} is given twice')
        samples[name] = vcf
    lengths = _fai_lengths(args.query_fasta_list, samples)

    with tempfile.TemporaryDirectory(prefix='gaf-shards-', dir=args.tmpdir or output) as directory:
        _log(f'streaming {len(vcfs)} merged VCF(s) into shards for {len(samples)} sample(s)')
        shards, contig_names, offsets = build_shards(
            vcfs, list(samples), Path(directory), args.shard_records)
        np.save(output / 'variant_offsets.npy', offsets)
        _log(f'{len(offsets)} VCF records indexed; shards in {directory}')
        tasks = [(name, vcf, str(shards[name]), contig_names, lengths.get(name, {}),
                  str(output / f'{name}.gaf')) for name, vcf in samples.items()]
        if args.processes > 1 and len(tasks) > 1:
            context = mp.get_context('fork')
            with context.Pool(min(args.processes, len(tasks))) as pool:
                results = list(pool.imap_unordered(write_sample_gaf, tasks))
        else:
            results = [write_sample_gaf(task) for task in tasks]
    with open(output / 'gaf_stats.tsv', 'w') as out:
        keys = sorted({key for _name, stats in results for key in stats})
        out.write('sample\t' + '\t'.join(keys) + '\n')
        for name, stats in sorted(results):
            out.write(name + '\t' + '\t'.join(str(stats.get(key, 0)) for key in keys) + '\n')
            _log(f'{name}: ' + ', '.join(f'{key}={stats.get(key, 0)}' for key in keys))
    return 0


if __name__ == '__main__':
    raise SystemExit(main())
