#!/usr/bin/env python3
"""Check that a per-sample VCF reconstructs its ##pseudoLinearMapping regions.

Read-only. For every mapped region in the header, the reference interval is
edited with the sample's own VCF records and compared with the assembly:

    orient(edit(REF[rs:re], records), strand) == ASSEMBLY[qs:qe]

Mappings that are contiguous on both axes (same query contig, same reference
path and strands) are checked as one run, because one event may cross a
mapping boundary. Conventions follow graphreftovcf.py:

* SV rows (TYPE INS/DEL/SUB): 0-based reference [POS, END) is replaced by
  INFO/SEQ. SEQ is the assembly slice at QUERYCOORD in the QUERYCOORD strand;
  it is compared with the assembly on its own (seq_vs_querycoord).
* SNP rows (TYPE SNP): 0-based reference POS-1 becomes ALT (forward reference
  strand). A '-' QUERYCOORD point is a traversal boundary: its base is at
  QUERYCOORD-1. INS_SNP observations restate bases already in their parent
  insertion and are ignored. SNPs inside an SV's replaced interval are counted
  (snp_inside_sv) and not applied twice.
* An observation belongs to a run when its assembly contig and QUERYCOORD
  lie inside the run's query interval and its CHROM/[POS,END) lie inside the
  run's reference interval.

Mappings whose Reference is '.' (no anchor) and UNMAPPED rows have no VCF
representation by construction; their bases are reported, not checked.
The second header line of a locus duplication (same Query and Category,
Reference = source interval) is ignored, as in gfa_sample_gaf.py.

Only the plain SEQ encoding (no --seqcompress) is supported.
"""
import argparse
from bisect import bisect_right
from collections import Counter, defaultdict
import csv
import gzip
import os
import re
import sys


_META = re.compile(r'(\w+)="((?:[^"\\]|\\.)*)"')
_INTERVAL = re.compile(r'(.+):(\d+)-(\d+)([+-])')
_COORD = re.compile(r'(\d+)(?:-(\d+))?([+-]?)')
_COMPLEMENT = str.maketrans('ACGTNacgtnRYKMSWBDHVrykmswbdhv',
                            'TGCANtgcanYRMKSWVHDByrmkswvhdb')
SV_TYPES = {'INS', 'DEL', 'SUB'}
REPORT_FIELDS = ('status', 'categories', 'query', 'reference', 'mappings',
                 'query_length', 'predicted_length', 'events', 'snps',
                 'first_mismatch', 'hamming', 'expected_context',
                 'predicted_context')


def revcomp(sequence):
    return sequence.translate(_COMPLEMENT)[::-1]


def orient(sequence, strand):
    return revcomp(sequence) if strand == '-' else sequence


def vcf_unescape(text):
    return (text.replace('@0A', '\n').replace('@09', '\t').replace('@2C', ',')
            .replace('@3B', ';').replace('@3A', ':').replace('@40', '@'))


def _open(path):
    return gzip.open(path, 'rt') if str(path).endswith('.gz') else open(path)


class FastaIndex:
    """Random access to one or more plain FASTA files through their .fai."""

    def __init__(self, paths):
        self.records = {}
        self.handles = {}
        for path in paths:
            index = str(path) + '.fai'
            if not os.path.isfile(index):
                raise SystemExit(f'{path}: missing .fai (run samtools faidx)')
            with open(index) as handle:
                for line in handle:
                    fields = line.rstrip('\n').split('\t')
                    if len(fields) < 5:
                        continue
                    name = fields[0]
                    if name in self.records and self.records[name][0] != path:
                        raise SystemExit(f'sequence {name!r} is in both '
                                         f'{self.records[name][0]} and {path}')
                    self.records[name] = (path, *map(int, fields[1:5]))

    def __contains__(self, name):
        return name in self.records

    def length(self, name):
        return self.records[name][1]

    def fetch(self, name, start, end):
        path, length, offset, line_bases, line_width = self.records[name]
        if start < 0 or end > length or end < start:
            raise ValueError(f'{name}:{start}-{end} is outside [0,{length})')
        if end == start:
            return ''
        handle = self.handles.get(path)
        if handle is None:
            handle = self.handles[path] = open(path, 'rb')
        first = offset + (start // line_bases) * line_width + start % line_bases
        last = offset + ((end - 1) // line_bases) * line_width + (end - 1) % line_bases
        handle.seek(first)
        data = handle.read(last - first + 1).decode('ascii')
        sequence = data.replace('\n', '').replace('\r', '')
        if len(sequence) != end - start:
            raise ValueError(f'{name}:{start}-{end}: short read from {path}')
        return sequence


def parse_meta(line):
    return {key: re.sub(r'\\(.)', r'\1', value) for key, value in _META.findall(line)}


def parse_interval(text):
    match = _INTERVAL.fullmatch(text or '')
    if match is None:
        return None
    return match[1], int(match[2]), int(match[3]), match[4]


def read_header(vcf):
    """Mappings (first line per Query/Category; per Query/Reference for a
    deletion, a query point), sample name, header flags."""
    mappings, seen, flags = [], set(), {}
    sample = None
    with _open(vcf) as handle:
        for line in handle:
            if line.startswith('#CHROM'):
                columns = line.rstrip('\n').split('\t')
                sample = columns[9] if len(columns) > 9 else None
                break
            if line.startswith('##minsetrefSequenceEncoding='):
                flags['encoding'] = line.strip().split('=', 1)[1]
            elif line.startswith('##querySequenceAlphabetFilter='):
                flags['alphabet'] = line.strip().split('=', 1)[1]
            elif line.startswith('##pseudoLinearMapping=<'):
                values = parse_meta(line)
                if values.get('Path') == 'alt':
                    continue          # sequence source: places no bases
                query = parse_interval(values.get('Query'))
                if query is None:
                    continue
                category = values.get('Category', '.')
                key = (values['Query'], category)
                if query[1] == query[2]:
                    key += (values.get('Reference'),)
                if key in seen:
                    continue
                seen.add(key)
                mappings.append((query, parse_interval(values.get('Reference')), category))
    return sample, mappings, flags


class Run:
    __slots__ = ('contig', 'qs', 'qe', 'qstrand', 'path', 'rs', 're', 'rstrand',
                 'categories', 'mappings', 'events', 'last')

    def __init__(self, query, reference, category):
        self.contig, self.qs, self.qe, self.qstrand = query
        self.path, self.rs, self.re, self.rstrand = reference
        self.categories = Counter([category])
        self.mappings = 1
        self.events = []
        self.last = (category, self.rs, self.re)

    def extend(self, query, reference, category):
        contig, qs, qe, qstrand = query
        path, rs, re_, rstrand = reference
        if (contig, qstrand, path, rstrand) != (self.contig, self.qstrand, self.path, self.rstrand):
            return False
        if qs != self.qe:
            return False
        if category == 'INSERTION' and self.last == ('DELETION', rs, re_):
            # A replacement listed as DELETION + INSERTION of one reference
            # span: the deletion already took the reference bases, the
            # insertion adds only query bases (its row spans both lines).
            self.qe = qe
            self.categories[category] += 1
            self.mappings += 1
            self.last = (category, rs, re_)
            return True
        if qstrand == rstrand:
            if rs != self.re:
                return False
            self.re = re_
        else:
            if re_ != self.rs:
                return False
            self.rs = rs
        self.qe = qe
        self.categories[category] += 1
        self.mappings += 1
        self.last = (category, rs, re_)
        return True


def build_runs(mappings, stats):
    # One claim per query interval: a full-locus-dup copy is listed both as
    # PRIMARY (replacing its reference span) and as INSERTION (a reference
    # point); both would put the same bases in two overlapping runs. Keep
    # the mapping with the longest reference interval. A deletion (a query
    # point) claims no query bases: several at one point are distinct
    # reference spans, all kept.
    best = {}
    for mapping in mappings:
        query, reference, category = mapping
        if category == 'UNMAPPED' or reference is None:
            best[(query, category)] = mapping
            continue
        if query[1] == query[2]:
            best[(query, reference)] = mapping
            continue
        span = reference[2] - reference[1]
        kept = best.get(query)
        if kept is None or span > kept[1][2] - kept[1][1]:
            best[query] = mapping
        if kept is not None:
            stats['duplicate_query_mappings'] += 1
    mappings = list(best.values())
    runs = []
    # Ties at one query point (deletions) in walking order along the
    # reference, so the one adjacent to the run extends it first.
    def walk_order(mapping):
        query, reference, _category = mapping
        if reference is None:
            return (query[0], query[1], query[2], 0)
        return (query[0], query[1], query[2],
                reference[1] if query[3] == reference[3] else -reference[2])

    for query, reference, category in sorted(mappings, key=walk_order):
        if category == 'UNMAPPED' or reference is None:
            key = 'unmapped' if category == 'UNMAPPED' else 'unanchored_' + category.lower()
            stats[key + '_mappings'] += 1
            stats[key + '_bases'] += query[2] - query[1]
            continue
        if runs and runs[-1].extend(query, reference, category):
            continue
        runs.append(Run(query, reference, category))
    return runs


def parse_coord(text):
    match = _COORD.fullmatch(text)
    if match is None:
        return None
    low = int(match[1])
    high = int(match[2]) if match[2] else low
    return low, high, match[3] or '+'


DUP_LIFT_CATEGORIES = frozenset({'Dup', 'DivergeDup'})


def read_events(vcf, stats, restated=None):
    """Yield (chrom, r0, r1, bases, contig, q0, q1, qstrand, kind, row).

    SV bases are in assembly (+) orientation; SNP bases are forward reference.
    restated: if a set, collects the IDs of SV rows lifted through a dup PA
    (LIFTCATEGORY Dup/DivergeDup), the copies' source-relative restatements.
    Records made only from Path=alt header entries (INFO/PAMAP; graphreftovcf's
    sequence-source layer: a copy's source or calls along it; older files:
    INFO PAPATH=alt) add no bases to the haplotype; they are counted
    (alternative_path_observations) and not yielded.
    """
    alt_entries = set()
    with _open(vcf) as handle:
        for line in handle:
            if line.startswith('#'):
                if line.startswith('##pseudoLinearMapping=<'):
                    values = parse_meta(line)
                    if values.get('Path') == 'alt' and values.get('ID', '').isdigit():
                        alt_entries.add(int(values['ID']))
                continue
            fields = line.rstrip('\r\n').split('\t')
            if len(fields) < 10:
                continue
            keys = fields[8].split(':')
            values = fields[9].split(':')
            if len(values) != len(keys):
                stats['malformed_rows'] += 1
                continue
            sample = dict(zip(keys, values))
            if sample.get('GT') != '1':
                continue
            chrom, pos = fields[0], int(fields[1])
            types = sample.get('TYPE', '.').split(',')
            info_items = fields[7].split(';')
            entries = next((item[6:] for item in info_items if item.startswith('PAMAP=')), '')
            if 'PAPATH=alt' in info_items or (entries and all(
                    int(entry) in alt_entries for entry in entries.split(','))):
                stats['alternative_path_observations'] += len(types)
                continue
            contigs = [vcf_unescape(name) for name in sample.get('ASSEMBLYCONTIG', '.').split(',')]
            coords = sample.get('QUERYCOORD', '.').split(',')
            if fields[4] in ('<INS>', '<DEL>', '<SUB>'):
                info = dict(item.split('=', 1) if '=' in item else (item, '')
                            for item in fields[7].split(';'))
                end = int(info['END'])
                sequence = info.get('SEQ', '.')
                for index, kind in enumerate(types):
                    coord = parse_coord(coords[index])
                    if coord is None or kind not in SV_TYPES:
                        stats['unparsed_sv_observations'] += 1
                        continue
                    q0, q1, qstrand = coord
                    contig = contigs[index if len(contigs) > 1 else 0]
                    if q1 > q0 and sequence in ('', '.'):
                        raise SystemExit(f'{fields[2]}: INS/SUB without INFO/SEQ; '
                                         'the --seqcompress encoding is not supported')
                    bases = '' if q1 == q0 else orient(sequence, qstrand)
                    if len(bases) != q1 - q0:
                        stats['sv_seq_length_mismatch'] += 1
                    # 'dup': a full-locus-dup copy; other calls inside its
                    # query span restate it relative to the source locus.
                    kind_name = 'dup' if info.get('PACLASS') == 'fulllocusdup' else 'sv'
                    if restated is not None and sample.get('LIFTCATEGORY') in DUP_LIFT_CATEGORIES:
                        restated.add(fields[2])
                    yield chrom, pos, end, bases, contig, q0, q1, qstrand, kind_name, fields[2]
            else:
                bases = sample.get('BASE', fields[4]).split(',')
                for index, kind in enumerate(types):
                    if kind != 'SNP':
                        stats['skipped_' + kind.lower() + '_observations'] += 1
                        continue
                    coord = parse_coord(coords[index])
                    if coord is None:
                        stats['unparsed_snp_observations'] += 1
                        continue
                    q0, _q1, qstrand = coord
                    if qstrand == '-':
                        q0 -= 1
                    contig = contigs[index if len(contigs) > 1 else 0]
                    yield (chrom, pos - 1, pos, fields[4].upper(), contig, q0, q0 + 1,
                           qstrand, 'snp', fields[2])


def assign(runs, events, query, stats, unassigned):
    by_contig = defaultdict(list)
    for run in runs:
        by_contig[run.contig].append(run)
    starts = {contig: [run.qs for run in values] for contig, values in by_contig.items()}
    for event in events:
        chrom, r0, r1, bases, contig, q0, q1, qstrand, kind, row = event
        if kind == 'sv' and q1 > q0 and contig in query:
            actual = query.fetch(contig, q0, q1)
            if actual.upper() != bases.upper():
                stats['seq_vs_querycoord_mismatch'] += 1
        candidates = by_contig.get(contig, [])
        index = bisect_right(starts.get(contig, []), q0) - 1
        chosen = None
        # Runs are disjoint on the query, but a boundary point (e.g. a pure
        # deletion) can touch two; the reference interval decides, then the
        # QUERYCOORD strand (at a strand switch both runs can hold the span).
        for probe in (index, index - 1, index + 1):
            if 0 <= probe < len(candidates):
                run = candidates[probe]
                if (run.qs <= q0 and q1 <= run.qe and run.path == chrom
                        and run.rs <= r0 and r1 <= run.re):
                    if chosen is None:
                        chosen = run
                    if run.qstrand == qstrand:
                        chosen = run
                        break
        if chosen is None:
            stats[f'unassigned_{kind}_observations'] += 1
            unassigned.append(event)
            continue
        chosen.events.append(event)
        stats[f'assigned_{kind}_observations'] += 1


def reconstruct(run, reference, stats):
    same = run.qstrand == run.rstrand
    # Ties (insertions at one point with one query start, as when a VCF
    # describes overlapping query spans twice) break on the query end and
    # bases, so the result does not depend on the order events were read in.
    events = sorted(run.events, key=lambda e: (e[1], e[2], e[8] == 'snp',
                                               e[5] if same else -e[5],
                                               e[6] if same else -e[6],
                                               (e[3] or '').upper()))
    # SV replaced intervals sorted by start, with the running maximum end, so
    # "some SV covers this SNP" is one binary search.
    sv_starts, sv_max_end = [], []
    for e in events:
        if e[8] == 'sv' and e[2] > e[1]:
            sv_starts.append(e[1])
            sv_max_end.append(max(e[2], sv_max_end[-1] if sv_max_end else e[2]))
    pieces, cursor, applied, snps = [], run.rs, 0, 0
    for chrom, r0, r1, bases, contig, q0, q1, qstrand, kind, row in events:
        covering = bisect_right(sv_starts, r0) - 1
        if kind == 'snp' and covering >= 0 and sv_max_end[covering] >= r1:
            stats['snp_inside_sv'] += 1
            continue
        if r0 < cursor:
            stats['overlapping_events'] += 1
            continue
        if kind == 'sv':
            # bases are in assembly orientation; flip them onto the forward
            # reference when the query maps to the opposite strand.
            forward = bases if same else revcomp(bases)
        else:
            forward = bases
            snps += 1
        pieces.append(reference.fetch(run.path, cursor, r0))
        pieces.append(forward)
        cursor = r1
        applied += 1
    pieces.append(reference.fetch(run.path, cursor, run.re))
    edited = ''.join(pieces)
    return (edited if same else revcomp(edited)), applied, snps


def compare(expected, predicted):
    first = next((i for i, (a, b) in enumerate(zip(expected, predicted)) if a != b),
                 None)
    if first is None and len(expected) != len(predicted):
        first = min(len(expected), len(predicted))
    hamming = (sum(a != b for a, b in zip(expected, predicted))
               if len(expected) == len(predicted) else '.')
    return first, hamming


def inside_copy(q0, q1, spans, restated=False):
    """Query interval [q0, q1) inside a full-locus-dup copy's span. A point
    (pure deletion) strictly inside is a restatement; on the copy's edge only
    when lifted through the dup PA (restated), otherwise it is the
    neighbouring primary call (graphvcfmerge._inside_dup_span)."""
    if q0 == q1:
        return any(a < q0 < b or (restated and a <= q0 <= b) for a, b in spans)
    return any(a <= q0 and q1 <= b for a, b in spans)


def check(vcf, reference_paths, query_path, report_path, max_report,
          events=None, stats=None):
    """events: replace the VCF's own records (check_merge_lossless.py)."""
    stats = Counter() if stats is None else stats
    sample, mappings, flags = read_header(vcf)
    if 'SEQ=omitted' in flags.get('encoding', ''):
        raise SystemExit(f'{vcf}: --seqcompress output (SEQ omitted) is not supported')
    if not mappings:
        raise SystemExit(f'{vcf}: no ##pseudoLinearMapping header lines')
    reference = FastaIndex(reference_paths)
    query = FastaIndex([query_path])
    runs = build_runs(mappings, stats)
    missing = {run.path for run in runs if run.path not in reference}
    missing |= {run.contig for run in runs if run.contig not in query}
    if missing:
        raise SystemExit(f'sequences not in the FASTA inputs: {sorted(missing)[:10]}')
    unassigned = []
    restated = set()
    events = list(read_events(vcf, stats, restated) if events is None else events)
    copies = defaultdict(list)
    for event in events:
        if event[8] == 'dup':
            copies[event[4]].append((event[5], event[6]))
    kept = []
    for event in events:
        if event[8] == 'dup':
            kept.append(event[:8] + ('sv',) + event[9:])
        elif inside_copy(event[5], event[6], copies.get(event[4], ()), event[9] in restated):
            stats['source_relative_in_dup'] += 1
        else:
            kept.append(event)
    assign(runs, kept, query, stats, unassigned)
    reported = 0
    with open(report_path, 'w', newline='') as handle:
        writer = csv.writer(handle, delimiter='\t', lineterminator='\n')
        writer.writerow(REPORT_FIELDS)
        for run in runs:
            expected = query.fetch(run.contig, run.qs, run.qe).upper()
            predicted, applied, snps = reconstruct(run, reference, stats)
            predicted = predicted.upper()
            bases = run.qe - run.qs
            stats['runs'] += 1
            stats['run_query_bases'] += bases
            if predicted == expected:
                stats['runs_exact'] += 1
                stats['run_query_bases_exact'] += bases
                continue
            first, hamming = compare(expected, predicted)
            kind = 'substitutions_only' if hamming != '.' else 'length_differs'
            stats['runs_mismatch_' + kind] += 1
            stats['run_query_bases_mismatch'] += bases
            if reported >= max_report:
                continue
            reported += 1
            writer.writerow((
                'MISMATCH_' + kind,
                ','.join(f'{name}:{count}' for name, count in sorted(run.categories.items())),
                f'{run.contig}:{run.qs}-{run.qe}{run.qstrand}',
                f'{run.path}:{run.rs}-{run.re}{run.rstrand}',
                run.mappings, bases, len(predicted), applied, snps,
                run.qs + first, hamming,
                expected[max(0, first - 10):first + 10],
                predicted[max(0, first - 10):first + 10],
            ))
        for chrom, r0, r1, _bases, contig, q0, q1, qstrand, kind, row in unassigned[:max_report]:
            writer.writerow(('UNASSIGNED_' + kind, '.', f'{contig}:{q0}-{q1}{qstrand}',
                             f'{chrom}:{r0}-{r1}', '.', q1 - q0, '.', '.', '.', '.', '.',
                             row, '.'))
    return sample, stats


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('-v', '--vcf', required=True, help='per-sample VCF (graphreftovcf)')
    parser.add_argument('-r', '--reference', required=True, action='append',
                        help='reference FASTA; repeat for the local reference templates')
    parser.add_argument('-q', '--query', required=True, help="the sample's assembly FASTA")
    parser.add_argument('-o', '--report', required=True, help='mismatch report TSV')
    parser.add_argument('--max-report', type=int, default=100000,
                        help='maximum rows written per category (default: 100000)')
    args = parser.parse_args(argv)
    sample, stats = check(args.vcf, args.reference, args.query, args.report, args.max_report)
    print(f'sample\t{sample}')
    for key in sorted(stats):
        print(f'{key}\t{stats[key]}')
    exact = stats['runs_exact'] == stats['runs'] and not any(
        key.startswith('unassigned_') or key == 'seq_vs_querycoord_mismatch'
        for key in stats if stats[key])
    print('lossless\t' + ('yes' if exact else 'no'))
    return 0 if exact else 1


if __name__ == '__main__':
    sys.exit(main())
