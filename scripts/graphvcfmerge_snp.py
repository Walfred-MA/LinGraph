#!/usr/bin/env python3
"""Exact, chromosome-wise cohort merge for graphreftovcf SNPs and small indels.

Scan each input once into compact shards. Each chromosome worker loads all its
observations into RAM, sorts with Polars, and unions identical coordinates and
alleles. No alignment, distance clustering, or disk-backed sort rounds.
"""
from __future__ import annotations

import argparse
from bisect import bisect_right
from collections import OrderedDict, defaultdict
from contextlib import contextmanager
import hashlib
import json
import os
from pathlib import Path
import shutil
import sys
import tempfile

import graphvcfmerge as vcf


def output_paths(prefix):
    prefix = str(prefix)
    if prefix.endswith('.vcf.gz'):
        prefix, extension = prefix[:-7], '.vcf.gz'
    else:
        prefix, extension = prefix.removesuffix('.vcf'), '.vcf'
    return tuple(prefix + suffix + extension for suffix in ('.snp', '.indel', '.sv'))


@contextmanager
def atomic_output(path):
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    suffix = '.gz' if str(path).endswith('.gz') else ''
    fd, temporary = tempfile.mkstemp(prefix=path.name + '.tmp.', suffix=suffix, dir=path.parent)
    os.close(fd)
    try:
        with vcf.open_output_text(temporary) as handle:
            yield handle
        os.replace(temporary, path)
    finally:
        if os.path.exists(temporary):
            os.unlink(temporary)


def _kind_and_key(row, cutoff):
    info = vcf.parse_info_field(row.info)
    if vcf.is_snp_format(row.fmt):
        return 'snp', (row.pos, row.ref, row.alt)
    if not vcf.is_sv_format(row.fmt):
        raise ValueError(f'unsupported graphreftovcf FORMAT {row.fmt!r}')
    size = abs(int(info.get('MAXSIZE', info.get('SVLEN', '0'))))
    if size >= cutoff:
        return None, None
    if size <= 0:
        raise ValueError('indel needs positive MAXSIZE or nonzero SVLEN')
    end = int(info.get('END', row.pos))
    kind = info.get('SVTYPE', row.alt.strip('<>'))
    sequence = vcf.vcf_unescape(info.get('SEQ', info.get('SVINSSEQ', '')))
    if kind in ('INS', 'SUB') and (not sequence or not sequence.isalpha()):
        raise ValueError('small INS/SUB requires plain INFO/SEQ for exact allele matching')
    # Symbolic <INS>/<SUB> alone does not identify an allele. Include its
    # actual bases and reference span; different alleles stay separate.
    return 'indel', (row.pos, row.ref, row.alt, end, kind, sequence.upper())


def _scan_one(task):
    index, path, directory, cutoff, snp_only = task
    directory = Path(directory) / str(index)
    directory.mkdir()
    handles, chroms, meta, samples = OrderedDict(), {}, [], []

    def write(chrom, text):
        key = hashlib.sha256(chrom.encode()).hexdigest()
        chroms[key] = chrom
        if key not in handles:
            if len(handles) >= 32:
                handles.popitem(last=False)[1].close()
            handles[key] = (directory / (key + '.vcf')).open('at')
        handles.move_to_end(key)
        handles[key].write(text.rstrip('\r\n') + '\n')

    try:
        with vcf.open_text(path) as handle:
            for number, raw in enumerate(handle, 1):
                try:
                    if raw.startswith('##referenceCoverage=<'):
                        values = vcf._parse_structured_meta(raw.strip(), 'referenceCoverage')
                        if values.get('Sample'):
                            write(values['Chrom'], raw)
                    elif raw.startswith('##'):
                        if not raw.startswith('##referenceCoverage'):
                            meta.append(raw.strip())
                    elif raw.startswith('#CHROM\t'):
                        samples = raw.rstrip('\r\n').split('\t')[9:]
                        if not samples or len(samples) != len(set(samples)):
                            raise ValueError('expected unique sample names in the VCF header')
                    elif raw.strip() and not raw.startswith('#'):
                        if not samples:
                            raise ValueError('missing #CHROM sample header')
                        row = vcf.parse_vcf_record(raw, number)
                        if len(row.samples) != len(samples):
                            raise ValueError('sample column count differs from header')
                        if snp_only and not vcf.is_snp_format(row.fmt):
                            continue
                        kind, _key = _kind_and_key(row, cutoff)
                        if kind:
                            write(row.chrom, raw)
                except (ValueError, KeyError) as error:
                    raise ValueError(f'{path}:{number}: {error}') from error
        if not samples:
            raise ValueError(f'{path}: missing #CHROM sample header')
    finally:
        for handle in handles.values():
            handle.close()
    return dict(index=index, samples=samples, chroms=chroms, meta=meta)


def scan(paths, shards_dir, cutoff=20, processes=1, snp_only=False, manifest_output=None, insertion_snps=None):
    if snp_only:
        from graphvcfmerge_snp_compact import scan as compact_scan
        return compact_scan(paths, shards_dir, processes, manifest_output, insertion_snps)
    if insertion_snps:
        raise ValueError('insertion SNP spools require snp_only=True')
    if not paths or cutoff < 1 or processes < 1:
        raise ValueError('scan needs inputs, positive --svcutoff and positive --processes')
    # Preserve the caller's path spelling (including symlink prefixes), since
    # Snakemake matches dynamic input paths against output patterns literally.
    root = Path(os.path.abspath(shards_dir))
    root.mkdir(parents=True, exist_ok=True)
    directory = tempfile.mkdtemp(prefix='scan-', dir=root)
    tasks = [(i, str(path), directory, cutoff, snp_only) for i, path in enumerate(paths)]
    if processes > 1 and len(tasks) > 1:
        with vcf.mp_context().Pool(min(processes, len(tasks))) as pool:
            results = pool.map(_scan_one, tasks)
    else:
        results = [_scan_one(task) for task in tasks]
    samples, metadata, seen, chroms, definitions = [], [], set(), {}, {}
    for result in results:
        samples.extend(name for name in result['samples'] if name not in samples)
        chroms.update(result['chroms'])
        for line in result['meta']:
            if line.startswith('##contig=<'):
                values = vcf._parse_structured_meta(line, 'contig')
                name, length = values.get('ID'), values.get('length')
                if name in definitions:
                    if definitions[name] != length:
                        raise ValueError(f'conflicting contig lengths for {name!r}')
                    continue
                definitions[name] = length
            if line not in seen:
                metadata.append(line)
                seen.add(line)
    order = vcf._contig_order(metadata)
    chroms = dict(sorted(chroms.items(), key=lambda item: (order.get(item[1], len(order)), item[1])))
    manifest = dict(protocol=1, directory=directory, cutoff=cutoff, samples=samples,
                    metadata=metadata, chroms=chroms,
                    files=[dict(index=r['index'], samples=r['samples'], chroms=list(r['chroms']))
                           for r in results])
    with atomic_output(root / 'manifest.json') as handle:
        json.dump(manifest, handle)
    if manifest_output:
        # Keep the small DAG-discovery record outside disposable shard data.
        with atomic_output(manifest_output) as handle:
            json.dump(manifest, handle)
    print(f'[merge:small] scanned {len(paths)} VCFs; {len(samples)} samples, {len(chroms)} chromosomes', file=sys.stderr)
    return manifest


def load_manifest(root):
    with open(Path(root) / 'manifest.json') as handle:
        manifest = json.load(handle)
    if manifest.get('protocol') in (2, 3, 4):
        manifest['chroms'] = {key: manifest['chroms'][key] for key in manifest['chrom_order']}
    return manifest


def part_path(manifest, key, kind):
    return Path(manifest['directory']) / 'parts' / (key + '.' + kind + '.part')


def merge_chrom(shards_dir, chrom, processes=1):
    manifest = load_manifest(shards_dir)
    if manifest.get("protocol") in (2, 3, 4):
        from graphvcfmerge_snp_compact import merge_chrom as compact_chrom
        return compact_chrom(shards_dir, chrom, processes)
    key = chrom if chrom in manifest['chroms'] else next(
        (key for key, name in manifest['chroms'].items() if name == chrom), None)
    if key is None:
        raise ValueError(f'unknown chromosome {chrom!r}')
    chrom = manifest['chroms'][key]
    sites, coverage = {}, defaultdict(list)
    for source in manifest['files']:
        if key not in source['chroms']:
            continue
        path = Path(manifest['directory']) / str(source['index']) / (key + '.vcf')
        with path.open() as handle:
            for number, raw in enumerate(handle, 1):
                if raw.startswith('##referenceCoverage=<'):
                    values = vcf._parse_structured_meta(raw.strip(), 'referenceCoverage')
                    coverage[values['Sample']].append((int(values['Start']), int(values['End'])))
                    continue
                row = vcf.parse_vcf_record(raw, number)
                kind, allele_key = _kind_and_key(row, manifest['cutoff'])
                site = sites.setdefault((kind, allele_key), dict(row=row, alleles=defaultdict(dict), states=defaultdict(set)))
                if row.filt == 'PASS':
                    site['row'].filt = 'PASS'
                for sample, field in zip(source['samples'], row.samples):
                    gt = field.split(':', 1)[0]
                    alleles = vcf.split_sample_alleles(field, row.fmt)
                    if not alleles and gt not in ('', '.', './.', '0', '0/0'):
                        raise ValueError(f'{path}:{number}: malformed allele for sample {sample!r}')
                    for allele in alleles:
                        values = vcf.parse_hsv_allele(allele)
                        if values is None:
                            raise ValueError(f'{path}:{number}: malformed observation for {sample!r}')
                        if kind == 'indel' and values[8] == '.':
                            shift, values[3] = vcf._fast_split_row_position_shift(values[3])
                            values[8] = str(shift)
                        # Exact duplicate observations are kept once; separate
                        # assembly locations and PA provenance remain intact.
                        site['alleles'][sample][vcf.format_hsv_allele(values)] = None
                    site['states'][sample].add(gt)
    coverage = {sample: vcf._merge_intervals(intervals) for sample, intervals in coverage.items()}

    def covered(sample, start, end):
        intervals = coverage.get(sample, ())
        index = bisect_right(intervals, (start, float('inf'))) - 1
        return index >= 0 and intervals[index][0] <= start and end <= intervals[index][1]

    for kind in ('snp', 'indel'):
        count = 0
        with atomic_output(part_path(manifest, key, kind)) as handle:
            for (site_kind, allele_key), site in sorted(sites.items()):
                if site_kind != kind or not site['alleles']:
                    continue
                row = site['row']
                fmt = vcf.SNP_FORMAT if kind == 'snp' else vcf.SV_FORMAT
                info = vcf.parse_info_field(row.info)
                start = max(0, row.pos - 1) if kind == 'snp' else row.pos
                end = start + 1 if kind == 'snp' else max(start + 1, int(info.get('END', row.pos)))
                fields = []
                for sample in manifest['samples']:
                    alleles = site['alleles'].get(sample, {})
                    if alleles:
                        fields.append(vcf.format_sample_alleles(sorted(alleles), fmt))
                    else:
                        states = site['states'].get(sample, set())
                        # Explicit missing calls take priority over coverage.
                        genotype = ('.' if states & {'.', './.', ''} else
                                    '0' if states & {'0', '0/0'} or covered(sample, start, end) else '.')
                        fields.append(vcf._state_sample_field(genotype, fmt))
                digest = hashlib.sha256(json.dumps([chrom, allele_key]).encode()).hexdigest()[:20]
                row.id = f'{kind.upper()}_{vcf._safe_variant_token(chrom)}_{row.pos}_{digest}'
                row.samples, row.fmt = fields, fmt
                row.info = vcf.update_info_nsup(row.info, sum(bool(a) for a in site['alleles'].values()))
                handle.write(vcf.format_vcf_record(row) + '\n')
                count += 1
        print(f'[merge:small] {chrom}: {count} {kind} rows', file=sys.stderr)
    return key


def _chrom_task(task):
    return merge_chrom(*task)


def concat(shards_dir, snp_output, indel_output=None, *, insertion_snps=None):
    manifest = load_manifest(shards_dir)
    from graphvcfmerge_snp_compact import concat_insertions, append_insertions
    metadata, insertions = concat_insertions(manifest, insertion_snps)
    targets = [(snp_output, 'snp')] + ([(indel_output, 'indel')] if indel_output else [])
    for target, kind in targets:
        parts = [part_path(manifest, key, kind) for key in manifest['chroms']]
        if any(not path.is_file() for path in parts):
            raise ValueError(f'incomplete chromosome parts for {kind}')
        with atomic_output(target) as handle:
            # Supply both public schemas even for an empty input category.
            vcf.write_vcf_header(handle, metadata, manifest['samples'], 'small')
            for part in parts:
                with part.open() as source:
                    shutil.copyfileobj(source, handle)
            if kind == 'snp':
                append_insertions(handle, manifest, insertions)


def merge(paths, snp_output, indel_output=None, *, cutoff=20, processes=1, tmpdir=None, snp_only=False, insertion_snps=None):
    with tempfile.TemporaryDirectory(prefix='graphvcfmerge-small-', dir=tmpdir) as root:
        manifest = scan(paths, root, cutoff, processes, snp_only=snp_only, insertion_snps=insertion_snps)
        tasks = [(root, key) for key in manifest['chroms']]
        if manifest.get("protocol") in (2, 3, 4):
            for _, key in tasks:
                merge_chrom(root, key, processes)
        elif processes > 1 and len(tasks) > 1:
            with vcf.mp_context().Pool(min(processes, len(tasks))) as pool:
                pool.map(_chrom_task, tasks)
        else:
            for task in tasks:
                _chrom_task(task)
        concat(root, snp_output, indel_output)


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('-i', '--input', nargs='+', default=[])
    parser.add_argument('-I', '--input-list', action='append', default=[])
    parser.add_argument('-o', '--output', help='prefix for .snp.vcf and .indel.vcf')
    parser.add_argument('--svcutoff', type=int, default=20, help='small indels have MAXSIZE < cutoff [20]')
    parser.add_argument('-t', '--processes', type=int, default=1,
                        help='scan workers or chromosome loading/sorting threads [1]')
    parser.add_argument('--stage', choices=('scan', 'prepare', 'chrom', 'concat'))
    parser.add_argument('--shards-dir')
    parser.add_argument('--manifest-output', help='scan: durable copy of the chromosome manifest')
    parser.add_argument('--raw-manifest', help='prepare: reuse saved original SNP shards')
    parser.add_argument('--insertion-snps', help='saved insertion SNP directory to sort and append at concat; scan/prepare defer loading it')
    parser.add_argument('--chrom', help='chromosome name or manifest key')
    parser.add_argument('--snp-only', action='store_true', help='ignore SV/indel input rows; at concat, write only the SNP VCF')
    args = parser.parse_args(argv)
    if args.svcutoff < 1 or args.processes < 1:
        parser.error('cutoff and process count must be positive')
    if args.stage and not args.shards_dir:
        parser.error('--stage requires --shards-dir')
    if args.stage == 'chrom':
        if not args.chrom:
            parser.error('--stage chrom requires --chrom')
        merge_chrom(args.shards_dir, args.chrom, args.processes)
    elif args.stage == 'prepare':
        from graphvcfmerge_snp_compact import prepare
        if not args.raw_manifest:
            parser.error('--stage prepare requires --raw-manifest')
        prepare(args.shards_dir, args.raw_manifest, args.manifest_output, args.insertion_snps)
    elif args.stage == 'concat':
        if not args.output:
            parser.error('--stage concat requires --output')
        snp_output, indel_output, _ = output_paths(args.output)
        concat(args.shards_dir, snp_output, None if args.snp_only else indel_output,
               insertion_snps=args.insertion_snps)
    else:
        paths = vcf.expand_vcf_inputs(args.input + vcf.read_vcf_input_lists(args.input_list))
        if not paths:
            parser.error('need input VCFs via -i or -I')
        if args.stage == 'scan':
            scan(paths, args.shards_dir, args.svcutoff, args.processes, snp_only=args.snp_only,
                 manifest_output=args.manifest_output, insertion_snps=args.insertion_snps)
        else:
            if not args.output:
                parser.error('--output is required')
            merge(paths, *output_paths(args.output)[:2], cutoff=args.svcutoff, processes=args.processes, snp_only=args.snp_only, insertion_snps=args.insertion_snps)
    return 0


if __name__ == '__main__':
    raise SystemExit(main())
