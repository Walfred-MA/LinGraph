#!/usr/bin/env python3
"""Check that merged VCFs still reconstruct every sample's mapped regions.

Read-only. Each sample's alleles are rebuilt from the merged VCFs alone
(sv, indel and snp outputs of graphvcfmerge.py / cohort_vcf_merge.py):

    top-level row carried by the sample
      -> representative INFO/SEQ (INS) or [POS, POS+|SVLEN|) (DEL)
      -> nested rows the sample carries inside that observation
         (CHROM = parent row ID, matched by QUERYCOORD), recursively
      -> INS_SNP rows (CHROM = insertion ID)

and then checked exactly like check_vcf_lossless.py: the reference over each
##pseudoLinearMapping run of the sample's own VCF, edited with these
alleles, must equal the assembly. The merged row POS is used, not
TEMPLATEOFFSET, because that is where the merged record places the allele.

Nested coordinates follow the merger: a nested insertion goes at offset
POS-1 of its parent sequence, a nested deletion removes
[POS-1, POS-1+|SVLEN|), and an INS_SNP replaces offset POS-1.

The per-sample VCFs give the mapped regions (-s); their own records are not
used. Standard library only.

Line checks (unification): every carrier of every row, top-level or nested,
must sit at the row's breakpoint (TEMPLATEOFFSET 0), and on every line its
allele - the row's representative edited by the nested rows and INS_SNPs it
carries - must equal its own assembly bases at QUERYCOORD. Failures are
listed in OUTPUT/line_check.tsv; the summary has a `unified` column.
"""
import argparse
from collections import Counter, defaultdict
import multiprocessing as mp
import os
from pathlib import Path
import re
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(1, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from allele_edges import edge_owned  # noqa: E402  (pipeline scripts)
from check_vcf_lossless import (  # noqa: E402
    FastaIndex, _open, check, orient, parse_coord, read_header, vcf_unescape,
)


SAMPLES = {}   # sample -> list of projected events (filled before forking)
LINE_ISSUES = []   # (sample, row, check, detail) from the line checks
QUERIES = {}       # sample -> assembly FASTA path (line checks)


class Record:
    __slots__ = ('id', 'chrom', 'pos', 'end', 'svtype', 'svlen', 'seq', 'alt', 'path')

    def __init__(self, fields, info):
        self.id, self.chrom, self.pos = fields[2], fields[0], int(fields[1])
        # Nested rows name their parent path: the row itself, or a shared
        # template path (full-locus duplications) given by >PATH:N=.
        named = re.fullmatch(r'[<>]([^:<>]+):\d+=', vcf_unescape(info.get('EXTENDGRAPHCIGAR', '')))
        self.path = named.group(1) if named else self.id
        self.alt = fields[4].upper()
        self.svtype = info.get('SVTYPE', 'SNP')
        self.end = int(info['END']) if 'END' in info else self.pos
        self.svlen = int(info['SVLEN']) if 'SVLEN' in info else 0
        seq = info.get('SEQ', '.')
        self.seq = '' if seq in ('', '.') else vcf_unescape(seq)


def read_queries(path):
    base = Path(path).resolve().parent
    output = {}
    with open(path) as handle:
        for line in handle:
            fields = line.split()
            if len(fields) >= 2 and not fields[0].startswith('#'):
                output[fields[0]] = str(base / fields[1])
    return output


def read_merged(vcfs, samples, stats):
    """Records carried by the target samples: [(record, sample, obs)]."""
    carried = []
    for vcf in vcfs:
        columns = {}
        with _open(vcf) as handle:
            for line in handle:
                if line.startswith('##'):
                    continue
                if line.startswith('#CHROM'):
                    header = line.rstrip('\n').split('\t')
                    columns = {index: name for index, name in enumerate(header)
                               if index >= 9 and name in samples}
                    continue
                fields = line.rstrip('\r\n').split('\t')
                picked = [(index, fields[index]) for index in columns
                          if index < len(fields) and fields[index].startswith('1')]
                if not picked:
                    continue
                info = dict(item.split('=', 1) if '=' in item else (item, '')
                            for item in fields[7].split(';'))
                record = Record(fields, info)
                keys = fields[8].split(':')
                for index, text in picked:
                    values = dict(zip(keys, text.split(':')))
                    if values.get('GT') != '1':
                        continue
                    types = values.get('TYPE', '.').split(',')
                    contigs = [vcf_unescape(name) for name in
                               values.get('ASSEMBLYCONTIG', '.').split(',')]
                    offsets = values.get('TEMPLATEOFFSET', '0').split(',')
                    allele_names = values.get('ALLELENAME', '.').split(',')
                    for number, coord in enumerate(values.get('QUERYCOORD', '.').split(',')):
                        parsed = parse_coord(coord)
                        if parsed is None:
                            stats['unparsed_observations'] += 1
                            continue
                        contig = contigs[number if len(contigs) > 1 else 0]
                        kind = types[number if len(types) > 1 else 0]
                        offset = offsets[number if len(offsets) > 1 else 0]
                        if offset not in ('0', '.', ''):
                            LINE_ISSUES.append((columns[index], record.id, 'breakpoint',
                                                f'TEMPLATEOFFSET={offset} at {coord}'))
                        low, high, strand = parsed
                        if record.svtype == 'SNP' and low == high:
                            # A SNP point on '-' is the traversal boundary
                            # (graphreftovcf): the base is the one before it.
                            low = low - 1 if strand == '-' else low
                            high = low + 1
                        # The ALLELENAME keeps two observations of one sample on
                        # the same assembly bases apart (a full-locus-dup copy
                        # lying inside another copy of the same template).
                        name = allele_names[number if len(allele_names) > 1 else 0]
                        carried.append((record, columns[index],
                                        (kind, contig, low, high, strand, name)))
    return carried


def project(carried, reference, stats_by_sample):
    """Per-sample events in check_vcf_lossless.py's tuple layout."""
    children = defaultdict(list)      # (sample, parent id) -> [(record, obs)]
    top = []
    for record, sample, obs in carried:
        if record.chrom in reference:
            top.append((record, sample, obs))
        else:
            children[(sample, record.chrom)].append((record, obs))
    used = set()

    def inside(obs, parent):
        _kind, contig, q0, q1, _strand = obs[:5]
        return contig == parent[1] and parent[2] <= q0 and q1 <= parent[3]

    # A nested row inside two observations of one sample on its path (e.g. a
    # full-locus-dup copy lying inside another copy of the same template on
    # the assembly) belongs to the one with its ALLELENAME.
    parents = defaultdict(list)       # (sample, path) -> [(record, obs)]
    for record, sample, obs in carried:
        parents[(sample, record.path)].append((record, obs))
    owner = {}
    for (sample, path), items in children.items():
        for child, child_obs in items:
            candidates = [(parent, parent_obs) for parent, parent_obs in parents.get((sample, path), ())
                          if inside(child_obs, parent_obs)]
            if len(candidates) < 2:
                continue
            name = child_obs[5]
            named = [(parent.id, parent_obs) for parent, parent_obs in candidates
                     if name not in ('', '.') and parent_obs[5] == name]
            if len(named) == 1:
                owner[(sample, child.id, child_obs)] = named[0]

    readers = {}

    def verify_line(record, obs, sample, bases, stats):
        """The carrier's allele on this line equals its assembly bases."""
        kind, contig, q0, q1, strand = obs[:5]
        if q1 <= q0 or sample not in QUERIES:
            return
        if sample not in readers:
            readers[sample] = FastaIndex([QUERIES[sample]])
        reader = readers[sample]
        stats['line_alleles_checked'] += 1
        if contig not in reader or q1 > reader.length(contig):
            stats['line_allele_mismatch'] += 1
            LINE_ISSUES.append((sample, record.id, 'allele', f'{contig}:{q0}-{q1} not in assembly'))
            return
        expected = reader.fetch(contig, q0, q1).upper()
        if orient(bases, strand).upper() != expected:
            stats['line_allele_mismatch'] += 1
            LINE_ISSUES.append((sample, record.id, 'allele',
                                f'{contig}:{q0}-{q1}{strand} len {len(bases)} vs {q1 - q0}'))

    def allele(record, obs, sample, stats, depth=0):
        """Parent sequence edited by the sample's nested rows, recursively."""
        sequence = record.seq
        edits = []
        rows, starts, ends, kept_spans = [], [], [], []
        current = len(record.seq)

        def extent(child):
            first = child.pos - 1
            return first, (first + abs(child.svlen) if child.svtype == 'DEL' else
                           first + 1 if child.svtype == 'SNP' else
                           child.end - 1 if child.end > child.pos else first)
        for child, child_obs in children.get((sample, record.path), ()):
            if not inside(child_obs, obs) or (sample, child.id, child_obs) in used:
                continue
            if owner.get((sample, child.id, child_obs), (record.id, obs)) != (record.id, obs):
                continue                 # another copy of this sample owns it
            first, last = extent(child)
            if (child_obs[2] == child_obs[3] and obs[2] < obs[3] and record.svtype != 'DEL'
                    and child_obs[2] in (obs[2], obs[3])):
                edge = starts if child_obs[2] == (obs[3] if obs[4] == '-' else obs[2]) else ends
                edge.append((first, last, (child, child_obs)))
            else:
                rows.append((child, child_obs))
                kept_spans.append((first, last))
                current += (child_obs[3] - child_obs[2]) - (last - first)
        if record.svtype != 'DEL':
            # Zero-length rows on this allele's query edge may belong to an
            # abutting allele of the same path (allele_edges.py).
            rows.extend(edge_owned(starts, ends, kept_spans, current, obs[3] - obs[2],
                                   len(record.seq)))
        for child, child_obs in rows:
            used.add((sample, child.id, child_obs))
            offset = child.pos - 1
            # Rows at one offset (e.g. insertions of different lines) follow
            # the sample's query order in the parent's orientation.
            order = child_obs[2] if child_obs[4] != '-' else -child_obs[3]
            if child_obs[0] == 'INS_SNP' or child.svtype == 'SNP':
                edits.append((offset, offset + 1, child.alt, order))
            elif child.svtype == 'DEL':
                # A deletion keeps only its own nested insertions (if any).
                edits.append((offset, offset + abs(child.svlen),
                              allele(child, child_obs, sample, stats, depth + 1), order))
            else:
                end = child.end - 1 if child.end > child.pos else offset
                edits.append((offset, end, allele(child, child_obs, sample, stats, depth + 1),
                              order))
            stats['nested_applied'] += 1
        # A deletion's nested rows are insertions on its deleted interval
        # (offsets 0..|SVLEN|); the deletion itself contributes no bases.
        limit = abs(record.svlen) if record.svtype == 'DEL' else len(sequence)
        pieces, cursor = [], 0
        for start, end, bases, _order in sorted(edits, key=lambda e: (e[0], e[1], e[3])):
            if start < cursor or end > limit:
                stats['nested_overlap_or_out_of_bounds'] += 1
                continue
            pieces.append(sequence[cursor:start])
            pieces.append(bases)
            cursor = end
        pieces.append(sequence[cursor:])
        result = ''.join(pieces)
        verify_line(record, obs, sample, result, stats)
        return result

    events = defaultdict(list)
    for record, sample, obs in top:
        stats = stats_by_sample[sample]
        kind, contig, q0, q1, qstrand = obs[:5]
        if record.svtype == 'SNP':
            if kind != 'SNP':
                stats[f'skipped_top_level_{kind.lower()}'] += 1
                continue
            point = q0   # already the SNP base (read_merged)
            events[sample].append((record.chrom, record.pos - 1, record.pos, record.alt,
                                   contig, point, point + 1, qstrand, 'snp', record.id))
        elif record.svtype == 'DEL':
            bases = allele(record, obs, sample, stats)
            events[sample].append((record.chrom, record.pos, record.pos + abs(record.svlen),
                                   orient(bases, qstrand), contig, q0, q1, qstrand, 'sv',
                                   record.id))
        else:
            bases = allele(record, obs, sample, stats)
            end = max(record.pos, record.end)
            events[sample].append((record.chrom, record.pos, end, orient(bases, qstrand),
                                   contig, q0, q1, qstrand, 'sv', record.id))
    for (sample, _parent), rows in children.items():
        for child, child_obs in rows:
            if (sample, child.id, child_obs) not in used:
                stats_by_sample[sample]['nested_without_carried_parent'] += 1
    return events


def check_sample(task):
    sample, vcf, reference_paths, query_path, report, max_report, extra = task
    stats = Counter(extra)
    _sample, stats = check(vcf, reference_paths, query_path, report, max_report,
                           events=SAMPLES.get(sample, []), stats=stats)
    return sample, dict(stats)


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('-v', '--merged-vcf', required=True, nargs='+',
                        help='merged VCFs (cohort sv, indel and snp outputs)')
    parser.add_argument('-s', '--sample-vcf', required=True, nargs='+',
                        help='per-sample VCFs, for their ##pseudoLinearMapping regions')
    parser.add_argument('-r', '--reference', required=True, action='append',
                        help='reference FASTA; repeat for the local reference templates')
    parser.add_argument('-q', '--query-paths', required=True,
                        help='NAME FASTA [FAI] per line (query_paths.normalized.txt)')
    parser.add_argument('-o', '--output-folder', required=True)
    parser.add_argument('-t', '--processes', type=int, default=1)
    parser.add_argument('--max-report', type=int, default=100000)
    args = parser.parse_args(argv)

    queries = read_queries(args.query_paths)
    vcfs = {}
    for path in args.sample_vcf:
        name, _mappings, _flags = read_header(path)
        if name not in queries:
            raise SystemExit(f'{path}: sample {name!r} is not in {args.query_paths}')
        vcfs[name] = path
    reference = FastaIndex(args.reference)
    QUERIES.update({name: queries[name] for name in vcfs})
    read_stats = Counter()
    carried = read_merged(args.merged_vcf, set(vcfs), read_stats)
    stats_by_sample = defaultdict(Counter)
    SAMPLES.update(project(carried, reference, stats_by_sample))
    output = Path(args.output_folder)
    output.mkdir(parents=True, exist_ok=True)
    tasks = [(name, path, args.reference, queries[name],
              str(output / f'{name}.merge_check.tsv'), args.max_report,
              dict(stats_by_sample[name])) for name, path in vcfs.items()]
    if args.processes > 1 and len(tasks) > 1:
        with mp.get_context('fork').Pool(args.processes) as pool:
            results = pool.map(check_sample, tasks, chunksize=1)
    else:
        results = [check_sample(task) for task in tasks]
    with open(output / 'line_check.tsv', 'w') as handle:
        handle.write('sample\trow\tcheck\tdetail\n')
        for issue in LINE_ISSUES:
            handle.write('\t'.join(map(str, issue)) + '\n')
    issues_by_sample = defaultdict(int)
    for issue in LINE_ISSUES:
        issues_by_sample[issue[0]] += 1
    results = [(name, dict(stats, line_issues=issues_by_sample.get(name, 0)))
               for name, stats in results]
    keys = sorted({key for _name, stats in results for key in stats})
    failed = 0
    with open(output / 'summary.tsv', 'w') as handle:
        handle.write('sample\t' + '\t'.join(keys) + '\tlossless\tunified\n')
        for name, stats in sorted(results):
            lossless = (stats.get('runs_exact', 0) == stats.get('runs', 0) and not any(
                stats.get(key) for key in keys
                if key.startswith(('unassigned_', 'nested_without', 'nested_overlap'))
                or key == 'seq_vs_querycoord_mismatch'))
            unified = not stats.get('line_issues')
            failed += not (lossless and unified)
            handle.write(name + '\t' + '\t'.join(str(stats.get(key, 0)) for key in keys)
                         + '\t' + ('yes' if lossless else 'no')
                         + '\t' + ('yes' if unified else 'no') + '\n')
    if read_stats:
        print('\n'.join(f'{key}\t{value}' for key, value in sorted(read_stats.items())))
    print((output / 'summary.tsv').read_text(), end='')
    return 1 if failed else 0


if __name__ == '__main__':
    sys.exit(main())
