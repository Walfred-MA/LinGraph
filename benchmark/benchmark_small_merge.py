#!/usr/bin/env python3
"""Compare native/Python SNP or indel merging on identical chromosome shards.

Each backend runs in a fresh process; timings include loading, sorting, grouping
and writing the dense VCF body. Input generation and native compilation are
excluded. The result is accepted only when output bytes match.
"""
import argparse
import hashlib
import json
import os
from pathlib import Path
import resource
import subprocess
import sys
import tempfile
import time

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / 'scripts'))


def worker(root, backend, threads, native_binary=None, kind='snp'):
    os.environ['GRAPHVCFMERGE_BACKEND'] = backend
    if native_binary:
        import graphvcfmerge_native as native
        native._BINARY = str(Path(native_binary).resolve())
    import graphvcfmerge_snp_compact as compact
    import graphvcfmerge_snp as small
    start = time.monotonic()
    small.merge_chrom(root, 'chr1', threads)
    elapsed = time.monotonic() - start
    path = small.part_path(small.load_manifest(root), compact.chrom_key('chr1'), kind)
    digest = hashlib.sha256()
    with path.open('rb') as handle:
        while block := handle.read(1 << 20):
            digest.update(block)
    rss = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss
    children = resource.getrusage(resource.RUSAGE_CHILDREN).ru_maxrss
    # Native worker and small Python parent overlap, so use their sum as a
    # conservative peak; Linux reports KiB and macOS reports bytes.
    scale = 2**20 if sys.platform == 'darwin' else 1024
    print(json.dumps(dict(backend=backend, seconds=elapsed,
                          peak_mib=(rss + children) / scale,
                          output_bytes=path.stat().st_size, sha256=digest.hexdigest())))


def indel_shards(root, args, names, rng):
    import graphvcfmerge as vcf
    import graphvcfmerge_snp_compact as compact
    key = compact.chrom_key('chr1')
    chroms = {key: 'chr1'}
    files = []
    for sample, name in enumerate(names):
        directory = root / str(sample)
        directory.mkdir()
        count = args.observations // args.samples + (sample < args.observations % args.samples)
        with (directory / (key + '.vcf')).open('w') as handle:
            handle.write(f'##referenceCoverage=<Chrom=chr1,Sample={name},Start=0,End={args.sites+4}>\n')
            for pos in rng.integers(1, args.sites + 1, count):
                kind = 'INS' if pos % 2 else 'DEL'
                end = pos if kind == 'INS' else pos + 3
                seq = 'ACG' if kind == 'INS' else '.'
                field = f'1:{kind}:3:>2H<1H>3I:assembly:{pos+100}+:.:.:allele:0H0H'
                handle.write(f'chr1\t{pos}\toriginal\tN\t<{kind}>\t.\tPASS\t'
                             f'SVTYPE={kind};END={end};MAXSIZE=3;SEQ={seq}\t{vcf.SV_FORMAT}\t{field}\n')
        files.append(dict(index=sample, samples=[name], chroms=[key]))
    manifest = dict(protocol=1, directory=str(root), cutoff=20, samples=names,
                    metadata=[], chroms=chroms, files=files)
    (root / 'manifest.json').write_text(json.dumps(manifest))


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--samples', type=int, default=2000)
    parser.add_argument('--observations', type=int, default=1_000_000)
    parser.add_argument('--sites', type=int, default=5000)
    parser.add_argument('--threads', type=int, default=4)
    parser.add_argument('--kind', choices=('snp', 'indel'), default='snp')
    parser.add_argument('--native-binary', type=Path,
                        help='use this native worker for before/after comparisons')
    parser.add_argument('--worker', nargs=2, metavar=('ROOT', 'BACKEND'))
    args = parser.parse_args()
    if args.worker:
        return worker(*args.worker, args.threads, args.native_binary, args.kind)
    if min(args.samples, args.observations, args.sites, args.threads) < 1:
        parser.error('sizes and threads must be positive')
    import numpy as np
    import graphvcfmerge_snp_compact as compact
    import graphvcfmerge_native as native
    if not args.native_binary:
        native.executable()
    with tempfile.TemporaryDirectory(prefix='small-merge-benchmark-') as directory:
        root = Path(directory)
        source = root / 'source'
        source.mkdir()
        names = [f's{i}' for i in range(args.samples)]
        rng = np.random.default_rng(42)
        if args.kind == 'indel':
            indel_shards(root, args, names, rng)
        else:
            records = np.empty(args.observations, dtype=compact.DTYPE)
            records['pos'] = rng.integers(1, args.sites + 1, args.observations, dtype=np.uint32)
            records['bases'] = 1  # A>C
            records['query_pos'] = records['pos'] + 100
            records['query'] = rng.integers(0, args.samples, args.observations, dtype=np.uint32)
            records['allele'] = 0
            records['label_h'] = 0
            records.tofile(source / 'records.bin')
            del records
            key = compact.chrom_key('chr1')
            data = dict(protocol=4, directory=str(source), queries=[
                [s, 'assembly', '+', 'SNP', '1', True, 'PASS'] for s in names],
                alleles=['allele'], label_groups=[], chroms={key: 'chr1'},
                coverage={'chr1': [[s, 0, args.sites] for s in names]},
                counts={key: args.observations}, sections={key: [[0, args.observations]]})
            (source / 'source.json').write_text(json.dumps(data))
            compact.manifest_from_sources([dict(source=str(source / 'source.json'),
                                               samples=names, metadata=[], chroms=data['chroms'])], root, str(root))
        results = []
        for backend in ('python', 'native'):
            command = [sys.executable, __file__, '--worker', str(root), backend,
                       '--threads', str(args.threads), '--kind', args.kind]
            if args.native_binary:
                command += ['--native-binary', str(args.native_binary)]
            output = subprocess.check_output(command, text=True)
            results.append(json.loads(output))
        if results[0]['sha256'] != results[1]['sha256']:
            raise RuntimeError('native/Python outputs differ')
        print(json.dumps(dict(kind=args.kind, samples=args.samples, observations=args.observations,
                              sites=args.sites, threads=args.threads, results=results,
                              speedup=results[0]['seconds'] / results[1]['seconds']), indent=2))


if __name__ == '__main__':
    main()
