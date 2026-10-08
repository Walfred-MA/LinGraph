#!/usr/bin/env python3
"""Rebuild each haplotype's assembly from the merged VCFs (no assemblies needed).

    rebuild_assemblies.py headers -s SAMPLE.vcf ... -o mappings.header
        One file with every ##pseudoLinearMapping line of the per-sample VCFs
        (each line names its Sample). It replaces the per-sample VCFs for
        `rebuild` and for gfa_sample_gaf.py -s.

    rebuild_assemblies.py rebuild -v MERGED.vcf ... -m mappings.header
                                  -r REFERENCE.fa [-r TEMPLATES.fa] -o OUTDIR
        OUTDIR/SAMPLE.fa (+ .fai) per sample and OUTDIR/query_paths.txt
        (NAME FASTA), to give merged_vcf_to_gfa.py / gfa_sample_gaf.py -q in
        place of the real assemblies.

    rebuild_assemblies.py rebuild -s SAMPLE.vcf ... -r REFERENCE.fa ... -o OUTDIR
        The same from the per-sample VCFs (each sample's own records and
        mapping lines), e.g. as the query-path list of cohort_vcf_merge.py
        --exact when the assemblies are not at hand.

Each contig is a run of N as long as the largest query coordinate its
mapping lines name; every mapped run is then filled with the reference over
the run, edited with the sample's merged alleles - the same reconstruction
tools/check_merge_lossless.py verifies against the real assemblies. Bases
outside runs (unmapped, insertions without a reference anchor, the
sequence-source layer) stay N; the GFA/GAF never use them. A contig longer
than its last mapped base gets that shorter length.
"""
import argparse
from collections import Counter, defaultdict
import multiprocessing as mp
import os
from pathlib import Path
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from check_vcf_lossless import (  # noqa: E402
    FastaIndex, _open, assign, build_runs, inside_copy, parse_interval, parse_meta,
    read_events, reconstruct,
)
import check_merge_lossless as merged  # noqa: E402

LINE_BASES = 80
SAMPLES = {}      # sample -> (mappings, contig lengths, events); filled before forking


def write_headers(vcfs, output):
    """Copy every ##pseudoLinearMapping line (and the coordinate-system line)."""
    written, seen = 0, set()
    with open(output + '.tmp', 'w') as out:
        for vcf in vcfs:
            with _open(vcf) as handle:
                for line in handle:
                    if line.startswith('#CHROM'):
                        break
                    if line.startswith('##pseudoLinearMappingCoordinateSystem'):
                        if line not in seen:
                            seen.add(line)
                            out.write(line)
                    elif line.startswith('##pseudoLinearMapping=<'):
                        if 'Sample=' not in line:
                            raise SystemExit(f'{vcf}: mapping line without Sample: {line[:120]}')
                        out.write(line)
                        written += 1
    os.replace(output + '.tmp', output)
    print(f'{written} mapping line(s) from {len(vcfs)} VCF(s) -> {output}')


def read_mappings(path):
    """Per sample: mappings as check_vcf_lossless.read_header keeps them, and
    each contig's length (the largest query coordinate of any of its lines).
    A line without Sample belongs to the file's #CHROM sample column."""
    mappings, seen = defaultdict(list), set()
    lengths = defaultdict(dict)
    with _open(path) as handle:
        for line in handle:
            if line.startswith('#CHROM'):
                columns = line.rstrip('\n').split('\t')
                if None in lengths and len(columns) > 9:
                    lengths[columns[9]] = lengths.pop(None)
                    if None in mappings:
                        mappings[columns[9]] = mappings.pop(None)
                break
            if not line.startswith('##pseudoLinearMapping=<'):
                continue
            values = parse_meta(line)
            sample = values.get('Sample') or None
            query = parse_interval(values.get('Query'))
            if query is None:
                continue
            contig = lengths[sample]
            contig[query[0]] = max(contig.get(query[0], 0), query[2])
            if values.get('Path') == 'alt':
                continue          # sequence source: places no bases
            category = values.get('Category', '.')
            key = (sample, values['Query'], category)
            if query[1] == query[2]:
                key += (values.get('Reference'),)
            if key in seen:
                continue
            seen.add(key)
            mappings[sample].append((query, parse_interval(values.get('Reference')), category))
    if None in lengths:
        raise SystemExit(f'{path}: mapping lines without Sample and no #CHROM sample column')
    return mappings, lengths


def sample_events(vcf, stats):
    """A per-sample VCF's own records, kept as check_vcf_lossless.check keeps them."""
    restated = set()
    events = list(read_events(vcf, stats, restated))
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
    return kept


def rebuild_sample(task):
    sample, references, output = task
    mappings, lengths, events = SAMPLES[sample]
    stats = Counter()
    reference = FastaIndex(references)
    runs = build_runs(mappings, stats)
    unassigned = []
    assign(runs, events, {}, stats, unassigned)
    by_contig = defaultdict(list)
    for run in runs:
        by_contig[run.contig].append(run)
    fasta = Path(output) / f'{sample}.fa'
    with open(str(fasta) + '.tmp', 'wb') as out, open(str(fasta) + '.fai.tmp', 'w') as fai:
        for contig in sorted(lengths):
            length = lengths[contig]
            bases = bytearray(b'N') * length
            for run in by_contig.get(contig, ()):
                predicted, _applied, _snps = reconstruct(run, reference, stats)
                if len(predicted) != run.qe - run.qs:
                    stats['runs_length_differs'] += 1
                    continue
                bases[run.qs:run.qe] = predicted.encode('ascii')
                stats['runs'] += 1
                stats['run_bases'] += run.qe - run.qs
            name = contig.encode()
            out.write(b'>' + name + b'\n')
            fai.write(f'{contig}\t{length}\t{out.tell()}\t{LINE_BASES}\t{LINE_BASES + 1}\n')
            view = memoryview(bases)
            for start in range(0, length, LINE_BASES):
                out.write(view[start:start + LINE_BASES])
                out.write(b'\n')
            stats['contigs'] += 1
            stats['contig_bases'] += length
    os.replace(str(fasta) + '.tmp', fasta)
    os.replace(str(fasta) + '.fai.tmp', str(fasta) + '.fai')
    stats['unassigned_observations'] = len(unassigned)
    return sample, dict(stats)


def rebuild(args):
    stats_by_sample = defaultdict(Counter)
    if args.sample_vcf:
        mappings, lengths, events = {}, {}, {}
        for vcf in args.sample_vcf:
            found, found_lengths = read_mappings(vcf)
            for sample in found:
                if sample in mappings:
                    raise SystemExit(f'sample {sample!r} is given twice')
                mappings[sample], lengths[sample] = found[sample], found_lengths[sample]
                events[sample] = sample_events(vcf, stats_by_sample[sample])
    else:
        if not args.merged_vcf or not args.mappings:
            raise SystemExit('rebuild needs -v and -m, or -s')
        mappings, lengths = read_mappings(args.mappings)
        reference = FastaIndex(args.reference)
        read_stats = Counter()
        carried = merged.read_merged(args.merged_vcf, set(mappings), read_stats)
        events = merged.project(carried, reference, stats_by_sample)
        del carried
    if not mappings:
        raise SystemExit('no ##pseudoLinearMapping lines')
    samples = sorted(mappings)
    for sample in samples:
        SAMPLES[sample] = (mappings[sample], lengths[sample], events.get(sample, []))
    output = Path(args.output_folder)
    output.mkdir(parents=True, exist_ok=True)
    tasks = [(sample, args.reference, str(output)) for sample in samples]
    if args.processes > 1 and len(tasks) > 1:
        with mp.get_context('fork').Pool(min(args.processes, len(tasks))) as pool:
            results = pool.map(rebuild_sample, tasks, chunksize=1)
    else:
        results = [rebuild_sample(task) for task in tasks]
    with open(output / 'query_paths.txt', 'w') as handle:
        for sample in samples:
            handle.write(f'{sample}\t{sample}.fa\n')
    keys = sorted({key for _sample, stats in results for key in stats}
                  | {key for stats in stats_by_sample.values() for key in stats})
    failed = 0
    with open(output / 'rebuild_summary.tsv', 'w') as handle:
        handle.write('sample\t' + '\t'.join(keys) + '\n')
        for sample, stats in sorted(results):
            stats = Counter(stats) + stats_by_sample[sample]
            failed += bool(stats.get('runs_length_differs') or stats.get('unassigned_observations')
                           or stats.get('overlapping_events'))
            handle.write(sample + '\t' + '\t'.join(str(stats.get(key, 0)) for key in keys) + '\n')
    print((output / 'rebuild_summary.tsv').read_text(), end='')
    return 1 if failed else 0


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    commands = parser.add_subparsers(dest='command', required=True)
    headers = commands.add_parser('headers', help='collect the mapping lines into one file')
    headers.add_argument('-s', '--sample-vcf', required=True, nargs='+')
    headers.add_argument('-o', '--output', required=True)
    build = commands.add_parser('rebuild', help='write the rebuilt assemblies')
    build.add_argument('-v', '--merged-vcf', nargs='+',
                       help='merged VCFs (cohort sv, indel and snp outputs)')
    build.add_argument('-m', '--mappings',
                       help='mapping header file (`headers`) for -v')
    build.add_argument('-s', '--sample-vcf', nargs='+',
                       help='instead of -v/-m: per-sample VCFs (their own records)')
    build.add_argument('-r', '--reference', required=True, action='append',
                       help='reference FASTA; repeat for the local reference templates')
    build.add_argument('-o', '--output-folder', required=True)
    build.add_argument('-t', '--processes', type=int, default=1)
    args = parser.parse_args(argv)
    if args.command == 'headers':
        write_headers(args.sample_vcf, args.output)
        return 0
    return rebuild(args)


if __name__ == '__main__':
    sys.exit(main())
