#!/usr/bin/env python3
"""Merge existing cohort VCFs and publish the requested variant categories."""
from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path
import pickle
import shutil

import graphvcfmerge as sv
import graphvcfmerge_snp as snp


def add_modes(parser):
    group = parser.add_mutually_exclusive_group()
    for name, help_text in (
        ('all', 'merge into separate cohort.sv.vcf, cohort.indel.vcf and cohort.snp.vcf'),
        ('svonly', 'merge SVs at or above --svcutoff into cohort.sv.vcf'),
        ('snp', 'merge only SNPs into cohort.snp.vcf'),
        ('svindel', 'merge SVs and split by size into cohort.sv.vcf and cohort.indel.vcf'),
    ):
        group.add_argument('--' + name, dest='merge_mode', action='store_const', const=name, help=help_text)
    return group


def resolve_mode(args):
    return getattr(args, 'merge_mode', None) or ('all' if getattr(args, 'mc_graph', False) else 'svonly')


def output_paths(output, mode):
    kinds = {'all': ('sv', 'indel', 'snp'), 'svindel': ('sv', 'indel'), 'svonly': ('sv',), 'snp': ('snp',)}[mode]
    return [Path(output) / f'cohort.{kind}.vcf' for kind in kinds]


def input_paths(output, listing=None):
    output = Path(output)
    if listing:
        paths = sv.expand_vcf_inputs(sv.read_vcf_input_lists([str(listing)]))
    else:
        listing = next((p for p in (output/'inputs/mergevcf.list', output/'lingraph/vcfs.list') if p.is_file()), None)
        if listing:
            return input_paths(output, listing)
        run_config = output/'inputs/cohort_call.run.json'
        if run_config.is_file():
            config = json.loads(run_config.read_text())
            paths = [str(output/'samples'/sample/(sample + '.vcf')) for sample in config['selected_samples']]
        else:
            paths = [str(p) for p in sorted((output/'samples').glob('*/*.vcf')) if p.name == p.parent.name + '.vcf']
    if not paths:
        raise FileNotFoundError(f'no existing per-sample VCFs or mergevcf.list found in {output}')
    for path in paths:
        if not Path(path).is_file():
            raise FileNotFoundError(f'missing merge input: {path}')
    return paths


def combine(paths, output):
    """Union compatible headers/samples and coordinate-sort VCF records on disk."""
    paths = [str(path) for path in paths]
    samples, metadata, seen, file_samples, definitions = [], [], set(), [], {}
    for path in paths:
        meta, names = sv.collect_vcf_header_info(path)
        if not names or len(names) != len(set(names)):
            raise ValueError(f'{path}: expected unique VCF sample columns')
        file_samples.append(names)
        for name in names:
            if name not in samples:
                samples.append(name)
        for line in meta:
            if line.startswith('##source=merge_locus_vcfs_'):
                continue
            if line.startswith('##contig=<'):
                values = sv._parse_structured_meta(line, 'contig')
                name, length = values['ID'], values.get('length')
                if name in definitions and definitions[name] != length:
                    raise ValueError(f'conflicting contig lengths for {name!r}')
                if name in definitions:
                    continue
                definitions[name] = length
            if line not in seen:
                metadata.append(line)
                seen.add(line)
    output = Path(output)
    output.parent.mkdir(parents=True, exist_ok=True)
    import tempfile
    with tempfile.TemporaryDirectory(prefix='combine-vcf-', dir=output.parent) as directory:
        body, ordered = Path(directory)/'rows.tsv', Path(directory)/'sorted.tsv'
        contig_order = sv._contig_order(metadata)
        with body.open('w') as handle:
            for path, names in zip(paths, file_samples):
                indexes = {name: index for index, name in enumerate(names)}
                with sv.open_text(path) as source:
                    for number, raw in enumerate(source, 1):
                        if raw.startswith('#') or not raw.strip():
                            continue
                        row = sv.parse_vcf_record(raw, number)
                        if len(row.samples) != len(names):
                            raise ValueError(f'{path}:{number}: sample column count differs from header')
                        row.samples = [row.samples[indexes[name]] if name in indexes else
                                       sv._state_sample_field('.', row.fmt) for name in samples]
                        handle.write(f'{contig_order.get(row.chrom, len(contig_order))}\t' + sv.format_vcf_record(row) + '\n')
        sv.run_sort(str(body), str(ordered), ['-k1,1n', '-k2,2', '-k3,3n', '-k5,5', '-k6,6', '-k4,4'])
        with snp.atomic_output(output) as handle:
            sv.write_vcf_header(handle, metadata, samples, 'small')
            with ordered.open() as source:
                for raw in source:
                    handle.write(raw.split('\t', 1)[1])


def cleanup_merge_directories(output, directories, *, protected=()):
    """Remove completed-merge scratch, never the output root or its parents."""
    output = Path(output).resolve()
    protected = [output, *(Path(path).resolve() for path in protected)]
    candidates = []
    for value in directories:
        path = Path(value)
        if path.is_symlink():
            raise ValueError(f'refusing to clean a symlinked merge directory: {path}')
        path = path.resolve()
        if any(path == target or path in target.parents for target in protected):
            raise ValueError(f'merge cleanup would remove a protected path: {path}')
        if path.exists() and not path.is_dir():
            raise ValueError(f'merge cleanup target is not a directory: {path}')
        if path not in candidates:
            candidates.append(path)
    # Validate every target before deleting any of them.
    for path in candidates:
        if path.exists():
            shutil.rmtree(path)
            print(f'[cohort-merge] removed temporary directory: {path}', flush=True)


def publish(output, mode, *, sv_input=None, snp_input=None, indel_input=None,
            cleanup_tmpdirs=()):
    targets = output_paths(output, mode)
    sources = {'all': [sv_input, indel_input, snp_input], 'svonly': [sv_input],
               'snp': [snp_input], 'svindel': [sv_input, indel_input]}[mode]
    if any(path is None or not Path(path).is_file() for path in sources):
        raise ValueError(f'{mode}: missing merged input VCFs: {sources}')
    if mode in ('svindel', 'all'):
        for source, target in zip(sources, targets):
            combine([source], target)
    else:
        combine(sources, targets[0])
    cleanup_merge_directories(output, cleanup_tmpdirs, protected=targets)
    return targets


def run(output, mode='svonly', *, listing=None, paths=None, processes=1, cutoff=20,
        merge_distance=500, size_similarity=.7, sequence_similarity=.7, var_in_insert=100,
        kmermatch=sv.DEFAULT_KMERMATCH, dry_run=False, keep_merge_tmpdir=False):
    paths = list(paths) if paths is not None else input_paths(output, listing)
    if not paths or processes < 1 or cutoff < 1:
        raise ValueError('merge requires input VCFs and positive process count/cutoff')
    print(f'[cohort-merge] mode={mode}; {len(paths)} input VCFs; {processes} workers', flush=True)
    if dry_run:
        print('[cohort-merge] outputs: ' + ', '.join(map(str, output_paths(output, mode))), flush=True)
        return output_paths(output, mode)
    stamps = [(str(Path(path).resolve()), Path(path).stat().st_size, Path(path).stat().st_mtime_ns) for path in paths]
    settings = dict(minsvsize=cutoff, merge_distance=merge_distance, size_similarity=size_similarity,
                    sequence_similarity=sequence_similarity, var_in_insert=var_in_insert,
                    emit_small=mode in ('all', 'svindel'))
    identity = hashlib.sha256(json.dumps([stamps, settings, mode == "all", "worker-snp-v3"]).encode()).hexdigest()[:20]
    root = Path(output)/'tmp'/'merge_only'/identity
    root.mkdir(parents=True, exist_ok=True)
    sv_output, snp_output = root/'sv.vcf', root/'snp.vcf'
    insertion_snps = str(root/'insertion_snps') if mode == 'all' else None
    snp_root = root/'snp_shards'
    if mode != 'snp':
        shards = str(root/'sv_shards')
        if not (Path(shards)/'chroms.txt').is_file():
            sv._fast_stage_scan(paths, shards, processes, snp_root if mode == 'all' else None)
        with open(Path(shards)/'manifest.pkl', 'rb') as handle:
            manifest = pickle.load(handle)
        sv._fast_stage_chrom(shards, sorted(manifest['chrom_sections']), processes=processes,
                             svcutoff=None, winnowmap=None, meryl='meryl', kmermatch=kmermatch,
                             insertion_snps=insertion_snps, **settings)
        sv._fast_stage_concat(shards, str(sv_output), output_small=settings['emit_small'],
                              processes=processes, insertion_snps=insertion_snps)
    if mode in ('all', 'snp'):
        from graphvcfmerge_snp_compact import scan, prepare, merge_chrom
        if mode == 'all':
            manifest = prepare(snp_root, snp_root/'raw.manifest.json', insertion_snps=insertion_snps)
        elif (snp_root/'manifest.json').is_file():
            manifest = json.loads((snp_root/'manifest.json').read_text())
        else:
            manifest = scan(paths, snp_root, processes)
        for key in manifest['chroms']:
            merge_chrom(snp_root, key, processes)
        snp.concat(snp_root, snp_output)
    targets = publish(output, mode, sv_input=sv_output, snp_input=snp_output,
                      indel_input=str(sv_output) + '.small.vcf')
    if not keep_merge_tmpdir:
        cleanup_merge_directories(output, [root], protected=[*paths, *targets])
    return targets


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('action', choices=('run', 'publish'))
    add_modes(parser)
    parser.add_argument('-O', '--output', required=True)
    parser.add_argument('-I', '--vcf-list')
    parser.add_argument('-t', '--processes', type=int, default=1)
    parser.add_argument('--svcutoff', type=int, default=20)
    parser.add_argument('--merge-distance', type=int, default=500)
    parser.add_argument('--size-similarity', type=float, default=.7)
    parser.add_argument('--sequence-similarity', type=float, default=.7)
    parser.add_argument('--var-in-insert', type=int, default=100)
    parser.add_argument('--kmermatch', default=sv.DEFAULT_KMERMATCH)
    parser.add_argument('--sv-input')
    parser.add_argument('--snp-input')
    parser.add_argument('--indel-input')
    parser.add_argument('--dry-run', action='store_true')
    parser.add_argument('--keep-merge-tmpdir', action='store_true',
                        help='retain merge intermediates after successful publication')
    parser.add_argument('--cleanup-tmpdir', action='append', default=[],
                        help='publish: remove this merge directory after all final VCFs are written')
    args = parser.parse_args(argv)
    mode = resolve_mode(args)
    if args.action == 'publish':
        publish(args.output, mode, sv_input=args.sv_input, snp_input=args.snp_input, indel_input=args.indel_input,
                cleanup_tmpdirs=() if args.keep_merge_tmpdir else args.cleanup_tmpdir)
    else:
        run(args.output, mode, listing=args.vcf_list, processes=args.processes, cutoff=args.svcutoff,
            merge_distance=args.merge_distance, size_similarity=args.size_similarity,
            sequence_similarity=args.sequence_similarity, var_in_insert=args.var_in_insert,
            kmermatch=args.kmermatch, dry_run=args.dry_run, keep_merge_tmpdir=args.keep_merge_tmpdir)
    return 0


if __name__ == '__main__':
    raise SystemExit(main())
