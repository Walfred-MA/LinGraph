#!/usr/bin/env python
"""Count merged cohort rows and their bases per nesting level: INS, DEL,
indel_INS, indel_DEL, SNP.

Level 0 rows sit on a reference chromosome; a row whose CHROM is another
row's ID (a nested row on an insertion path, or an insertion SNP) is one level
deeper than that row (a CHROM that looks like a row ID, I_/D_/INS_/DEL_/...,
but is no nested row is taken as a level-0 row). SV and indel rows are split by size (INFO/MAXSIZE, else
|SVLEN|): >= CUTOFF is INS/DEL, smaller is indel_INS/indel_DEL; every row of
the SNP file is a SNP. Bases: the row's size (MAXSIZE, else |SVLEN|; the
representative's inserted or deleted bases), 1 per SNP; columns <category>_bp.

Usage:
  python tools/count_variants_by_level.py MERGE_OUTPUT_DIR
  python tools/count_variants_by_level.py --sv cohort.sv.vcf --indel cohort.indel.vcf --snp cohort.snp.vcf
"""
import argparse
from collections import Counter
from concurrent.futures import ProcessPoolExecutor
from pathlib import Path
import re
import subprocess
import sys

CATEGORIES = ('INS', 'DEL', 'indel_INS', 'indel_DEL', 'SNP')
NESTED_ID = re.compile(r'^(I|D|INS|DEL|SUB|DUP)_')


def columns(path, fields):
    """Yield the tab-separated FIELDS (cut -f syntax) of each body line; the
    sample columns are dropped by cut before Python sees them."""
    if str(path).endswith(('.gz', '.bgz')):
        source = subprocess.Popen(['gzip', '-cd', '--', str(path)], stdout=subprocess.PIPE)
        cut = subprocess.Popen(['cut', '-f', fields], stdin=source.stdout, stdout=subprocess.PIPE)
        source.stdout.close()
        processes = [source, cut]
    else:
        cut = subprocess.Popen(['cut', '-f', fields, '--', str(path)], stdout=subprocess.PIPE)
        processes = [cut]
    with cut.stdout as handle:
        for raw in handle:
            if raw[:1] == b'#':
                continue
            yield raw.decode().rstrip('\n').split('\t')
    for process in processes:
        if process.wait() != 0:
            raise RuntimeError(f'{path}: {process.args[0]} exited with {process.returncode}')


def info_value(info, key):
    match = re.search(r'(?:^|;)' + key + r'=([^;]*)', info)
    return match.group(1) if match else None


def scan_sv(path, cutoff):
    """Counter[(CHROM, category)] of rows and of bases, and {row ID: CHROM}
    of the nested rows; top-level rows are level 0, so they need no entry."""
    counts, bases, parents = Counter(), Counter(), {}
    for chrom, row_id, info in columns(path, '1,3,8'):
        svtype = info_value(info, 'SVTYPE') or '.'
        size = info_value(info, 'MAXSIZE')
        if size in (None, '', '.'):
            size = (info_value(info, 'SVLEN') or '0').lstrip('-')
        size = int(size) if size.isdigit() else 0
        if svtype in ('INS', 'DEL'):
            category = svtype if size >= cutoff else 'indel_' + svtype
        else:
            category = 'other_' + svtype
        counts[chrom, category] += 1
        bases[chrom, category] += size
        if NESTED_ID.match(chrom):
            parents[row_id] = chrom
    return counts, bases, parents


def scan_snp(path):
    counts = Counter()
    for (chrom,) in columns(path, '1'):
        counts[chrom, 'SNP'] += 1
    return counts, counts, {}


def find(directory, kind):
    for name in (f'cohort.{kind}.vcf', f'cohort.{kind}.vcf.gz'):
        path = Path(directory) / name
        if path.is_file():
            return path
    return None


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('directory', nargs='?', help='merge output folder holding cohort.{sv,indel,snp}.vcf[.gz]')
    parser.add_argument('--sv', help='cohort.sv.vcf[.gz]')
    parser.add_argument('--indel', help='cohort.indel.vcf[.gz]')
    parser.add_argument('--snp', help='cohort.snp.vcf[.gz]')
    parser.add_argument('-m', '--cutoff', type=int, default=20,
                        help='SV size cutoff of the merge (default 20, the pipeline --svcutoff)')
    parser.add_argument('-o', '--output', help='TSV to write (default stdout)')
    args = parser.parse_args(argv)

    files = {kind: getattr(args, kind) or (find(args.directory, kind) if args.directory else None)
             for kind in ('sv', 'indel', 'snp')}
    if not any(files.values()):
        parser.error('no cohort.{sv,indel,snp}.vcf found; give a folder or --sv/--indel/--snp')
    for kind, path in files.items():
        print(f'[count-levels] {kind}: {path or "missing, skipped"}', file=sys.stderr, flush=True)

    jobs = [(scan_snp, (path,)) if kind == 'snp' else (scan_sv, (path, args.cutoff))
            for kind, path in files.items() if path]
    counts, bases, parents = Counter(), Counter(), {}
    with ProcessPoolExecutor(max_workers=len(jobs)) as pool:
        for part_counts, part_bases, part_parents in [future.result() for future in
                                                      [pool.submit(function, *call) for function, call in jobs]]:
            counts.update(part_counts)
            bases.update(part_bases)
            parents.update(part_parents)

    levels = {}

    def level(chrom):
        """Level of the rows on CHROM: 0 on a reference chromosome, else one
        more than the level of the row CHROM names."""
        chain, current = [], chrom
        while current not in levels and NESTED_ID.match(current) and current in parents:
            if current in chain:
                raise ValueError(f'cyclic nesting through {current}')
            chain.append(current)
            current = parents[current]
        if current not in levels:
            # A reference chromosome, or a nested-looking ID of a top-level row.
            levels[current] = 1 if NESTED_ID.match(current) else 0
        value = levels[current]
        for row_id in reversed(chain):
            value += 1
            levels[row_id] = value
        return levels[chrom]

    table, table_bp = {}, {}
    for (chrom, category), number in counts.items():
        key = level(chrom)
        table.setdefault(key, Counter())[category] += number
        table_bp.setdefault(key, Counter())[category] += bases[chrom, category]
    extra = sorted({category for row in table.values() for category in row} - set(CATEGORIES))
    names = [*CATEGORIES, *extra]
    header = ['level', *names, 'total', *(name + '_bp' for name in names), 'total_bp']

    def line(key, row, row_bp):
        return '\t'.join([str(key), *(str(row[c]) for c in names), str(sum(row.values())),
                          *(str(row_bp[c]) for c in names), str(sum(row_bp.values()))])

    total, total_bp = Counter(), Counter()
    lines = ['\t'.join(header)]
    for key in sorted(table):
        total.update(table[key])
        total_bp.update(table_bp[key])
        lines.append(line(key, table[key], table_bp[key]))
    lines.append(line('all', total, total_bp))
    text = '\n'.join(lines) + '\n'
    if args.output:
        Path(args.output).write_text(text)
    sys.stdout.write(text)


if __name__ == '__main__':
    main()
