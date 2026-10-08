#!/usr/bin/env python3
"""Check that per-sample GAF walks spell the sample's assembly exactly.

Every record must also be a walk of the graph: consecutive steps need an L
line (records_missing_link).

Read-only. For every record of SAMPLE.gaf (from gfa_sample_gaf.py) the path
is spelled from the GFA segment sequences, cut to [path_start, path_end), and
compared with the assembly bases [query_start, query_end) of that record
(reverse-complemented for strand '-'), case-insensitively.

With -s/--sample-vcf, each sample's ##pseudoLinearMapping regions (anchored,
non-UNMAPPED; same selection as check_vcf_lossless.py) are compared with the
query bases that exact GAF records cover, so bases lost by breaks or dropped
intervals are counted too. The mapping category at the first mismatching
base is written to the report to separate insertion collapse from reference
spans.

Requires the GFA S lines to carry sequences. Only the segments and links the
GAF paths use are loaded; with -t N the GFA is read in N parallel chunks.
numpy, when present, locates mismatches.
"""
import argparse
from bisect import bisect_right
from collections import Counter
import csv
import multiprocessing as mp
import os
from pathlib import Path
import re
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from check_vcf_lossless import FastaIndex, _open, read_header, revcomp  # noqa: E402


_STEP = re.compile(r'([<>])(\d+)')
REPORT_FIELDS = ('status', 'query', 'path_start', 'path_end', 'query_length',
                 'path_length', 'variants', 'first_mismatch', 'mismatch_category',
                 'hamming', 'expected_context', 'walked_context')
try:
    import numpy as np
except ImportError:          # pragma: no cover - plain loops instead
    np = None

SEGMENTS = {}
LINKS = set()   # canonical link keys (_link) of the L lines the GAF paths need
NEEDED_SEGMENTS = set()
NEEDED_LINKS = set()
GFA_PATH = None


def _link(a, b):
    """Canonical key of the directed link a -> b, steps encoded id*2+reverse:
    a link and its reverse share one key."""
    ra, rb = b ^ 1, a ^ 1
    if (ra, rb) < (a, b):
        a, b = ra, rb
    return (a << 32) | b


def _steps(path):
    return [int(identifier) * 2 + (direction == '<')
            for direction, identifier in _STEP.findall(path)]


def missing_links(path):
    """Consecutive steps of a GAF path with no L line between them."""
    steps = _steps(path)
    return sum(_link(a, b) not in LINKS for a, b in zip(steps, steps[1:]))


def _needed(gafs):
    """Segments and links the GAF paths use."""
    for gaf in gafs:
        with open(gaf) as handle:
            for line in handle:
                fields = line.split('\t', 6)
                if len(fields) < 6:
                    continue
                steps = _steps(fields[5])
                NEEDED_SEGMENTS.update(step >> 1 for step in steps)
                NEEDED_LINKS.update(_link(a, b) for a, b in zip(steps, steps[1:]))


def _load_chunk(span):
    """Needed S sequences and L keys of the GFA bytes [start, end)."""
    start, end = span
    segments, links = [], []
    position = start
    with open(GFA_PATH, 'rb') as handle:
        handle.seek(start)
        for line in handle:
            if position >= end:
                break
            position += len(line)
            if line.startswith(b'L\t'):
                fields = line.split(b'\t', 5)
                key = _link(int(fields[1]) * 2 + (fields[2] == b'-'),
                            int(fields[3]) * 2 + (fields[4] == b'-'))
                if key in NEEDED_LINKS:
                    links.append(key)
            elif line.startswith(b'S\t'):
                fields = line.split(b'\t', 3)
                identifier = int(fields[1])
                if identifier in NEEDED_SEGMENTS:
                    sequence = fields[2].rstrip(b'\n')
                    if sequence == b'*':
                        raise SystemExit(f'{GFA_PATH}: segment {identifier} has no sequence')
                    segments.append((identifier, sequence.decode().upper()))
    return segments, links


def load_segments(gfa, gafs, processes=1):
    """Load the S sequences and L lines that the GAF paths need."""
    global GFA_PATH
    _needed(gafs)
    if str(gfa).endswith('.gz') or processes <= 1:
        with _open(gfa) as handle:
            for line in handle:
                if line.startswith('L\t'):
                    fields = line.split('\t', 5)
                    key = _link(int(fields[1]) * 2 + (fields[2] == '-'),
                                int(fields[3]) * 2 + (fields[4] == '-'))
                    if key in NEEDED_LINKS:
                        LINKS.add(key)
                elif line.startswith('S\t'):
                    fields = line.split('\t', 3)
                    if int(fields[1]) in NEEDED_SEGMENTS:
                        sequence = fields[2].rstrip('\n')
                        if sequence == '*':
                            raise SystemExit(f'{gfa}: segment {fields[1]} has no sequence')
                        SEGMENTS[int(fields[1])] = sequence.upper()
        return
    GFA_PATH = str(gfa)
    size = os.path.getsize(GFA_PATH)
    bounds = [size * part // processes for part in range(processes + 1)]
    with open(GFA_PATH, 'rb') as handle:
        for part in range(1, processes):
            handle.seek(bounds[part])
            handle.readline()           # chunks start at line starts
            bounds[part] = max(bounds[part - 1], handle.tell())
    with mp.get_context('fork').Pool(processes) as pool:
        for segments, links in pool.imap_unordered(_load_chunk, list(zip(bounds[:-1], bounds[1:]))):
            SEGMENTS.update(segments)
            LINKS.update(links)


def _first_mismatch(expected, walked):
    """(index of the first differing base, Hamming distance or None)."""
    same = len(expected) == len(walked)
    if np is not None:
        a = np.frombuffer(expected.encode(), dtype=np.uint8)
        b = np.frombuffer(walked.encode(), dtype=np.uint8)
        n = min(len(a), len(b))
        diff = np.flatnonzero(a[:n] != b[:n])
        first = int(diff[0]) if len(diff) else n
        return first, (int(len(diff)) if same else None)
    first = next((i for i, (x, y) in enumerate(zip(expected, walked)) if x != y),
                 min(len(expected), len(walked)))
    return first, (sum(x != y for x, y in zip(expected, walked)) if same else None)


def spell(path, start, end):
    pieces, offset = [], 0
    for direction, identifier in _STEP.findall(path):
        if offset >= end:
            break
        sequence = SEGMENTS[int(identifier)]
        length = len(sequence)
        if offset + length > start:
            if direction == '<':
                sequence = revcomp(sequence)
            pieces.append(sequence[max(0, start - offset):end - offset])
        offset += length
    return ''.join(pieces)


def read_queries(query_list):
    fastas = {}
    base = Path(query_list).resolve().parent
    with open(query_list) as handle:
        for line in handle:
            fields = line.split()
            if len(fields) >= 2:
                fasta = Path(fields[1])
                fastas[fields[0]] = str(fasta if fasta.is_absolute() else base / fasta)
    return fastas


def merge(intervals):
    """Interval union per contig: {contig: [[low, high], ...]}."""
    merged = {}
    for contig, low, high in sorted(intervals):
        spans = merged.setdefault(contig, [])
        if spans and low <= spans[-1][1]:
            spans[-1][1] = max(spans[-1][1], high)
        else:
            spans.append([low, high])
    return merged


def union_length(intervals):
    return sum(high - low for spans in merge(intervals).values() for low, high in spans)


def covered_by(intervals, others):
    """Bases of the union of `intervals` that also lie in the union of `others`."""
    cover = merge(others)
    total = 0
    for contig, spans in merge(intervals).items():
        targets, index = cover.get(contig, []), 0
        for low, high in spans:
            while index < len(targets) and targets[index][1] <= low:
                index += 1
            probe = index
            while probe < len(targets) and targets[probe][0] < high:
                a, b = targets[probe]
                total += min(b, high) - max(a, low)
                probe += 1
    return total


def mapping_category(mappings, contig, position):
    starts = mappings['starts'].get(contig)
    if not starts:
        return '.'
    index = bisect_right(starts, position) - 1
    rows = mappings['rows'][contig]
    for probe in (index, index - 1):
        if 0 <= probe < len(rows):
            low, high, category = rows[probe]
            if low <= position < high or (low == high == position):
                return category
    return 'outside_mappings'


def check_sample(task):
    sample, gaf, fasta, vcf, report, max_report = task
    stats = Counter()
    query = FastaIndex([fasta])
    mappings = {'starts': {}, 'rows': {}}
    mapped = []
    if vcf:
        _name, rows, _flags = read_header(vcf)
        for (contig, low, high, _strand), reference, category in rows:
            mappings['rows'].setdefault(contig, []).append((low, high, category))
            if category != 'UNMAPPED' and reference is not None:
                mapped.append((contig, low, high))
        for contig, values in mappings['rows'].items():
            values.sort()
            mappings['starts'][contig] = [low for low, _high, _category in values]
    exact_spans = []
    reported = 0
    with open(gaf) as handle, open(report, 'w', newline='') as out:
        writer = csv.writer(out, delimiter='\t', lineterminator='\n')
        writer.writerow(REPORT_FIELDS)
        for line in handle:
            fields = line.rstrip('\n').split('\t')
            if len(fields) < 12:
                continue
            contig, qs, qe, strand = fields[0], int(fields[2]), int(fields[3]), fields[4]
            path, ps, pe = fields[5], int(fields[7]), int(fields[8])
            variants = next((tag[5:] for tag in fields[12:] if tag.startswith('nv:i:')), '.')
            expected = query.fetch(contig, qs, qe).upper()
            if strand == '-':
                expected = revcomp(expected)
            walked = spell(path, ps, pe)
            stats['records'] += 1
            gaps = missing_links(path)
            if gaps:
                # The walk is not a path of the graph.
                stats['records_missing_link'] += 1
                stats['missing_links'] += gaps
            stats['record_query_bases'] += qe - qs
            if walked == expected:
                stats['records_exact'] += 1
                stats['record_query_bases_exact'] += qe - qs
                exact_spans.append((contig, qs, qe))
                continue
            first, hamming = _first_mismatch(expected, walked)
            same = len(expected) == len(walked)
            kind = 'substitutions_only' if same else 'length_differs'
            stats['records_mismatch_' + kind] += 1
            position = qs + first if strand == '+' else qe - 1 - first
            category = mapping_category(mappings, contig, position) if vcf else '.'
            stats[f'first_mismatch_in_{category}'] += 1
            if reported < max_report:
                reported += 1
                writer.writerow((
                    'MISMATCH_' + kind, f'{contig}:{qs}-{qe}{strand}', ps, pe,
                    qe - qs, pe - ps, variants, position, category,
                    hamming if same else '.',
                    expected[max(0, first - 10):first + 10],
                    walked[max(0, first - 10):first + 10]))
    if vcf:
        stats['mapped_bases'] = union_length(mapped)
        stats['mapped_bases_in_exact_gaf'] = covered_by(mapped, exact_spans)
        stats['mapped_bases_not_exact'] = stats['mapped_bases'] - stats['mapped_bases_in_exact_gaf']
    return sample, dict(stats)


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('-g', '--gfa', required=True, help='rGFA from merged_vcf_to_gfa.py')
    parser.add_argument('-a', '--gaf', required=True, nargs='+',
                        help='SAMPLE.gaf files from gfa_sample_gaf.py')
    parser.add_argument('-q', '--query-fasta-list', required=True,
                        help='NAME FASTA [FAI] per line (the gfa_sample_gaf.py list)')
    parser.add_argument('-s', '--sample-vcf', nargs='+', default=[],
                        help='per-sample VCFs, for mapped-region coverage')
    parser.add_argument('-o', '--output-folder', required=True,
                        help='writes SAMPLE.gaf_check.tsv and summary.tsv here')
    parser.add_argument('-t', '--processes', type=int, default=1)
    parser.add_argument('--max-report', type=int, default=100000)
    args = parser.parse_args(argv)

    fastas = read_queries(args.query_fasta_list)
    vcfs = {}
    for vcf in args.sample_vcf:
        name, _rows, _flags = read_header(vcf)
        vcfs[name] = vcf
    output = Path(args.output_folder)
    output.mkdir(parents=True, exist_ok=True)
    tasks = []
    for gaf in args.gaf:
        sample = Path(gaf).name.removesuffix('.gaf')
        if sample not in fastas:
            raise SystemExit(f'{gaf}: sample {sample!r} is not in {args.query_fasta_list}')
        tasks.append((sample, gaf, fastas[sample], vcfs.get(sample),
                      str(output / f'{sample}.gaf_check.tsv'), args.max_report))
    load_segments(args.gfa, args.gaf, args.processes)
    if args.processes > 1 and len(tasks) > 1:
        with mp.get_context('fork').Pool(args.processes) as pool:
            results = pool.map(check_sample, tasks, chunksize=1)
    else:
        results = [check_sample(task) for task in tasks]
    keys = sorted({key for _sample, stats in results for key in stats})
    with open(output / 'summary.tsv', 'w') as handle:
        handle.write('sample\t' + '\t'.join(keys) + '\tlossless\n')
        failed = 0
        for sample, stats in results:
            lossless = (stats.get('records_exact', 0) == stats.get('records', 0)
                        and not stats.get('mapped_bases_not_exact', 0)
                        and not stats.get('records_missing_link', 0))
            failed += not lossless
            handle.write(sample + '\t' + '\t'.join(str(stats.get(key, 0)) for key in keys)
                         + '\t' + ('yes' if lossless else 'no') + '\n')
    print((output / 'summary.tsv').read_text(), end='')
    return 1 if failed else 0


if __name__ == '__main__':
    sys.exit(main())
