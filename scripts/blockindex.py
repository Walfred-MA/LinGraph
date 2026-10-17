#!/usr/bin/env python3
"""Build reference-coordinate grouped blocks using the established 50-kb pipeline."""
from __future__ import annotations

import argparse
import fcntl
import json
import os
from pathlib import Path
import re
import shlex
import shutil
import subprocess
import sys


HERE = Path(__file__).resolve().parent
PROTOCOL = 1


def atomic_write(path, content):
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_name(path.name + f'.tmp.{os.getpid()}')
    temporary.write_text(content)
    temporary.replace(path)


def stamps(paths):
    return [[str(p), p.stat().st_size, p.stat().st_mtime_ns] for p in map(Path, paths)]


class Stages:
    """Reuse only completed stages whose inputs and output inventory still match."""
    def __init__(self, root, dry_run=False):
        self.root, self.dry_run = root, dry_run

    def run(self, name, command, inputs, outputs, action=None):
        command = list(map(str, command))
        print(f'[blockindex] {name}: {shlex.join(command)}', flush=True)
        if self.dry_run:
            return
        marker = self.root / 'checkpoints' / (name + '.json')
        signature = dict(protocol=PROTOCOL, command=command, inputs=stamps(inputs),
                         driver=stamps([Path(__file__)]))
        try:
            saved = json.loads(marker.read_text())
            if (saved['signature'] == signature and saved['outputs']
                    == stamps(row[0] for row in saved['outputs'])):
                print(f'[blockindex] Reuse {name}', flush=True)
                return
        except (OSError, ValueError, KeyError, TypeError):
            pass
        marker.unlink(missing_ok=True)
        if action is None:
            subprocess.run(command, check=True)
        else:
            action()
        inventory = list(outputs() if callable(outputs) else outputs)
        if not inventory:
            raise ValueError(f'{name}: no outputs')
        atomic_write(marker, json.dumps(dict(signature=signature, outputs=stamps(inventory)), indent=2) + '\n')


def read_lengths(fasta):
    lengths = {}
    file_size = fasta.stat().st_size
    with fasta.open('rb') as source:
        first_header = source.readline()
        first_offset = source.tell()
    if not first_header.startswith(b'>') or not first_header[1:].split():
        raise ValueError('input must start with a nonempty FASTA header')
    first_name = first_header[1:].split()[0].decode('ascii')
    with Path(str(fasta) + '.fai').open() as handle:
        for line in handle:
            fields = line.split()
            if len(fields) < 5:
                raise ValueError('invalid FASTA index; regenerate it with samtools faidx')
            name = fields[0]
            length, offset, width, stride = map(int, fields[1:5])
            if name in lengths or length <= 0 or offset < 0 or width <= 0 or stride < width:
                raise ValueError(f'invalid FASTA index entry: {name}')
            if not lengths and (name != first_name or offset != first_offset):
                raise ValueError('FASTA/index first record mismatch; regenerate with samtools faidx')
            last_base = offset + ((length - 1) // width) * stride + (length - 1) % width
            if last_base >= file_size:
                raise ValueError(f'FASTA index entry exceeds the input file: {name}')
            lengths[name] = length
    if not lengths:
        raise ValueError('empty FASTA index')
    return lengths


def export_blocks(localblocks, output, lengths, sample):
    """Publish retained templates, including source loci with no filtered lift."""
    targets = sorted(localblocks.glob('*/*.fasta'),
                     key=lambda p: [int(x) if x.isdigit() else x
                                    for x in re.split(r'(\d+)', p.parent.name)])
    targets = [p for p in targets if p.name == p.parent.name + '.fasta']
    rows = set()
    for target in targets:
        with target.with_suffix('.bed').open() as handle:
            for line in handle:
                fields = line.split()
                if len(fields) < 7 or fields[6] != sample:
                    raise ValueError(f'{target}: template is not from reference {sample}')
                contig, start, end, group = fields[:4]
                start, end = int(start), int(end)
                if (contig not in lengths or not 0 <= start < end <= lengths[contig]
                        or group != target.parent.name):
                    raise ValueError(f'{target}: invalid source interval or group: {line.strip()}')
                rows.add((contig, start, end, group))
    if not targets or not rows:
        raise ValueError('No blocks retained; inspect the reference masking and --unmask threshold')
    order = {name: i for i, name in enumerate(lengths)}
    rows = sorted(rows, key=lambda row: (order[row[0]], *row[1:]))
    atomic_write(output / 'blocks.bed', ''.join('\t'.join(map(str, row)) + '\n' for row in rows))
    atomic_write(output / 'partition_targets.list', ''.join(str(p.resolve()) + '\n' for p in targets))
    covered = {name: [] for name in lengths}
    for contig, start, end, _group in rows:
        covered[contig].append((start, end))
    gaps = []
    for contig, length in lengths.items():
        cursor = 0
        for start, end in sorted(covered[contig]):
            if start > cursor:
                gaps.append((contig, cursor, start))
            cursor = max(cursor, end)
        if cursor < length:
            gaps.append((contig, cursor, length))
    atomic_write(output / 'uncovered.bed', ''.join('\t'.join(map(str, row)) + '\n' for row in gaps))
    atomic_write(output / 'blocks.stats.json', json.dumps(dict(
        partitions=len(targets), intervals=len(rows), reference_bases=sum(lengths.values()),
        uncovered_bases=sum(end - start for _, start, end in gaps)), indent=2) + '\n')


def build_parser():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('-i', '--input', required=True, help='one uncompressed reference FASTA with adjacent .fai')
    parser.add_argument('-s', '--sample', required=True, help='reference name, e.g. Species_h1; native contig names are preserved')
    parser.add_argument('-o', '-O', '--output', required=True, help='dedicated block-index directory')
    parser.add_argument('-t', '--threads', type=int, default=16)
    parser.add_argument('--block-size', type=int, default=50000)
    parser.add_argument('--minimum-block-size', type=int, default=20000)
    parser.add_argument('--unmask', type=int, default=1000, help='minimum uppercase bases per retained partition')
    parser.add_argument('--resume', action='store_true', help='accepted for compatibility; completed matching stages are always reused')
    parser.add_argument('--dry-run', action='store_true', help='validate inputs and print stages without writing or running commands')
    return parser


def pipeline(args, source, output, lengths):
    work = output / 'work'
    clean = work / 'selfclean' / 'main_chroms.fa'
    kmers = work / 'reference.kmers.bin'
    raw, named = work / 'blocks.raw.bed', work / 'blocks.named.bed'
    blocks = work / 'blocks.fasta'
    aligned, filtered = work / 'aligned.bed', work / 'aligned.filtered.bed'
    merged = work / 'merged.bed'
    groups = Path(str(merged) + '.groups.tsv')
    local = output / 'localblocks'
    listing = output / 'query_paths.txt'
    targets = output / 'partition_targets.list'
    index = Path(str(targets) + '.bin')
    marker = Path(str(index) + '.fasta_target_unique_v2')
    dummy, dummy_list = work / 'index_dummy.fa', work / 'index_dummy.list'
    stages = Stages(output, args.dry_run)
    python = sys.executable

    def run(name, command, inputs, outputs, action=None):
        scripts = [Path(c) for c in command if str(c).endswith('.py') and Path(c).is_file()]
        stages.run(name, command, [*inputs, *scripts], outputs, action)

    if not args.dry_run:
        work.mkdir(parents=True, exist_ok=True)
        content = f'{args.sample}\t{source}\n'
        if not listing.exists() or listing.read_text() != content:
            atomic_write(listing, content)
    run('self-clean', [python, HERE / 'self_clean_reference.py', '-i', source,
        '-o', clean.parent, '-s', args.sample, '-t', args.threads, '--resume'],
        [source, Path(str(source) + '.fai'), HERE / 'minsetref.py', HERE / 'minsetref_selfclean.py',
         HERE / 'minsetref_align.py', HERE / 'minsetref_core.py', HERE / 'minsetref_segments.py'], [clean])
    run('index-clean-fasta', ['samtools', 'faidx', clean], [clean], [Path(str(clean) + '.fai')])
    run('kmer-index', [HERE / 'kmerindex', '-f', clean, '-o', kmers, '-t', args.threads],
        [clean, Path(str(clean) + '.fai'), HERE / 'kmerindex'], [kmers])
    run('blocking', [HERE / 'kmerblocking', '-i', kmers, '-f', str(clean) + '.fai',
        '-o', raw, '-t', args.threads, '--block-size', args.block_size],
        [kmers, Path(str(clean) + '.fai'), HERE / 'kmerblocking'], [raw])

    def name_blocks():
        with raw.open() as handle:
            rows = [line.split()[:3] for line in handle if line.strip() and not line.startswith('#')]
        if not rows or any(len(row) != 3 for row in rows):
            raise ValueError('blocking produced no valid BED3 intervals')
        atomic_write(named, ''.join('\t'.join([*row, f'g{i}']) + '\n' for i, row in enumerate(rows, 1)))

    run('name-blocks', ['name-blocks', raw, named], [raw], [named], name_blocks)
    run('extract-blocks', [python, HERE / 'extract_named_blocks.py', '-f', clean,
        '-b', named, '-o', blocks, '--sample', args.sample], [clean, named], [blocks])
    run('index-block-fasta', ['samtools', 'faidx', blocks], [blocks], [Path(str(blocks) + '.fai')])
    run('map-reference', [python, HERE / 'map_assemblies_to_block_reference.py',
        '-q', source, '-t', blocks, '-o', aligned, '--threads', args.threads],
        [source, Path(str(source) + '.fai'), blocks, Path(str(blocks) + '.fai')], [aligned])
    run('filter-mappings', [python, HERE / 'filter_block_mapping_bed.py', aligned,
        '-o', filtered, '--minimum-block-size', args.minimum_block_size],
        [aligned, HERE / 'map_assemblies_to_block_reference.py'], [filtered])
    run('merge-partitions', [HERE / 'kmerpartitionmerger', '-f', source, '-T', blocks,
        '-b', filtered, '-o', merged, '-m', '10', '-t', args.threads, '--only-kmer'],
        [source, blocks, filtered, HERE / 'kmerpartitionmerger'], [merged, groups])
    command = [python, HERE / 'build_localblocks_from_filtered_bed.py', '-b', filtered,
               '--merged-bed', merged, '-q', listing, '-t', blocks,
               '-o', local, '--unmask', args.unmask]

    def localblocks():
        # This directory is owned by this dedicated blockindex run. Build into
        # scratch first so removed groups cannot survive as stale targets.
        scratch = work / 'localblocks.next'
        if scratch.exists():
            shutil.rmtree(scratch)
        actual = list(map(str, command))
        actual[actual.index('-o') + 1] = str(scratch)
        subprocess.run(actual, check=True)
        if local.exists():
            shutil.rmtree(local)
        scratch.replace(local)

    run('local-blocks', command, [filtered, merged, groups, listing, blocks, source],
        lambda: sorted(local.glob('*/*.bed')) + sorted(local.glob('*/*.fasta')), localblocks)
    local_files = [] if args.dry_run else sorted(local.glob('*/*.bed')) + sorted(local.glob('*/*.fasta'))
    run('export-blocks', ['export-blocks', local, output], local_files,
        [output / 'blocks.bed', targets, output / 'uncovered.bed', output / 'blocks.stats.json'],
        lambda: export_blocks(local, output, lengths, args.sample))

    def index_targets():
        atomic_write(dummy, '>index_dummy\n' + 'N' * 32 + '\n')
        atomic_write(dummy_list, f'index_dummy\t{dummy}\n')
        index.unlink(missing_ok=True)
        marker.unlink(missing_ok=True)
        subprocess.run(list(map(str, index_command)), check=True)
        if not index.is_file() or index.stat().st_size == 0:
            raise ValueError('kmer_searcher did not produce a nonempty target index')

    index_command = [HERE / 'kmer_searcher', '-T', targets, '-I', dummy_list,
                     '-o', str(work) + '/', '-c', '500', '-n', args.threads]
    run('index-partitions', index_command, [targets, *local_files, HERE / 'kmer_searcher'],
        [index, marker], index_targets)
    print(f'[blockindex] {"Planned" if args.dry_run else "Ready"}: {output / "blocks.bed"}', flush=True)


def main(argv=None):
    parser = build_parser()
    args = parser.parse_args(argv)
    if args.threads < 1 or args.block_size < 1 or min(args.minimum_block_size, args.unmask) < 0:
        parser.error('threads/block-size must be positive; minimum-block-size/unmask must be nonnegative')
    if not re.fullmatch(r'[A-Za-z0-9.-]+_h[1-9][0-9]*', args.sample):
        parser.error('--sample must have the form Species_h1')
    source, output = Path(args.input).expanduser().resolve(), Path(args.output).expanduser().resolve()
    if (not source.is_file() or source.suffix.lower() in {'.gz', '.bgz', '.bgzf'}
            or any(c.isspace() for c in str(source) + str(output))):
        parser.error('input must be an uncompressed FASTA; input/output paths cannot contain whitespace')
    if output == source or output in source.parents:
        parser.error('use a dedicated output directory that does not contain the input reference')
    lengths = read_lengths(source)
    state = dict(protocol=PROTOCOL, reference=stamps([source, Path(str(source) + '.fai')]),
                 sample=args.sample, block_size=args.block_size,
                 minimum_block_size=args.minimum_block_size, unmask=args.unmask)
    config = output / 'blockindex.json'
    if config.exists():
        if json.loads(config.read_text()) != state:
            parser.error('reference or blocking parameters changed; use a new output directory')
    elif output.exists() and any(output.iterdir()):
        parser.error('output is not an existing blockindex run or an empty directory')
    if args.dry_run:
        pipeline(args, source, output, lengths)
        return 0
    for executable in ['samtools', 'minimap2', 'winnowmap', 'meryl', 'blastn', 'makeblastdb', 'sort',
                       *(str(HERE / name) for name in ('kmerindex', 'kmerblocking', 'kmerpartitionmerger', 'kmer_searcher'))]:
        if shutil.which(executable) is None:
            raise FileNotFoundError(f'required executable not found: {executable}')
    output.mkdir(parents=True, exist_ok=True)
    with (output / '.blockindex.lock').open('a') as lock:
        try:
            fcntl.flock(lock, fcntl.LOCK_EX | fcntl.LOCK_NB)
        except BlockingIOError:
            raise RuntimeError(f'another blockindex run is active in {output}') from None
        if config.exists() and json.loads(config.read_text()) != state:
            raise ValueError('blockindex inputs changed while waiting for the output directory')
        if not config.exists():
            atomic_write(config, json.dumps(state, indent=2) + '\n')
        pipeline(args, source, output, lengths)
    return 0


if __name__ == '__main__':
    try:
        raise SystemExit(main())
    except (OSError, ValueError, RuntimeError, subprocess.CalledProcessError) as error:
        print(f'[blockindex] ERROR: {error}', file=sys.stderr)
        raise SystemExit(1)
