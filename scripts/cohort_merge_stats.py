#!/usr/bin/env python
"""Summary statistics of the published cohort merge (cohort.summary.tsv).

cohort_vcf_merge.py publish counts every row while it writes cohort.sv.vcf,
cohort.sv.10Ksingleton.vcf, cohort.indel.vcf and cohort.snp.vcf (no extra
pass) and writes <output>/cohort.summary.tsv:

  file  kind  level  rows  bases  templates  template_bases

kind: SNP (cohort.snp.vcf), INS / DEL / SUB / ... (SV files: INFO SVTYPE),
INDEL_INS / INDEL_DEL / ... (cohort.indel.vcf).
level: top (CHROM is a reference sequence) or nested (CHROM is a row or
template ID: a variant inside an insertion).
bases: INS / SUB / DUP = length of INFO SEQ (an @XX escape is one base), or
|SVLEN| without SEQ; DEL = |SVLEN| deleted reference bases; SNP = one per
ALT allele.
templates / template_bases: the sequence a row is, counted once (first path
of EXTENDGRAPHCIGAR when it is a variant/template ID, e.g. the shared
DUP_chrom_start_end of full-locus-dup site rows, else the row ID); SNP rows
are their own templates.
Rows with file = ALL add up the files per kind and level.

The same table for files already published:
  python cohort_merge_stats.py OUTPUT_DIR            (writes OUTPUT_DIR/cohort.summary.tsv)
"""
import argparse
import os
import re
import sys
from collections import defaultdict
from pathlib import Path

SUMMARY = 'cohort.summary.tsv'
NESTED_PREFIXES = (b'INS_', b'DEL_', b'I_', b'D_', b'SNP_', b'SUB_', b'DUP_')
INSERTED = (b'INS', b'SUB', b'DUP')
TEMPLATE = re.compile(rb'^[<>]([^<>:]+):')
HEADER = 'file\tkind\tlevel\trows\tbases\ttemplates\ttemplate_bases'


def _info_field(info, key):
    for item in info.split(b';'):
        name, _, value = item.partition(b'=')
        if name == key:
            return value
    return None


def _template(row_id, info):
    cigar = _info_field(info, b'EXTENDGRAPHCIGAR')
    if cigar:
        match = TEMPLATE.match(cigar.replace(b'@3A', b':'))
        if match and match.group(1).startswith(NESTED_PREFIXES):
            return match.group(1)
    return row_id


def file_label(path):
    """INDEL_ for cohort.indel.vcf, SNP for cohort.snp.vcf, '' for SV files."""
    name = os.path.basename(str(path))
    if '.indel.' in name:
        return 'INDEL_'
    if '.snp.' in name:
        return 'SNP'
    return ''


class RowStats:
    """Counts the rows written to one published VCF; add() takes one row (str or bytes)."""

    def __init__(self, path):
        self.path = str(path)
        self.label = file_label(path)
        self.reset()

    def reset(self):
        """Forget everything (the rows are written again from the start)."""
        self.counts = defaultdict(lambda: [0, 0])
        self.templates = defaultdict(dict)

    def add(self, line):
        if isinstance(line, str):
            line = line.encode()
        if line.startswith(b'#'):
            return
        cols = line.split(b'\t', 8)
        if len(cols) < 8:
            return
        chrom, row_id, ref, alt, info = cols[0], cols[2], cols[3], cols[4], cols[7]
        level = 'nested' if chrom.startswith(NESTED_PREFIXES) else 'top'
        svtype = _info_field(info, b'SVTYPE')
        if svtype is None and not alt.startswith(b'<'):
            alts = alt.split(b',')
            if len(ref) == 1 and all(len(a) == 1 for a in alts):
                c = self.counts[('SNP', level)]
                c[0] += 1
                c[1] += len(alts)
                return
            svtype = b'INS' if max(len(a) for a in alts) > len(ref) else b'DEL'
        svtype = svtype or alt.strip(b'<>') or b'.'
        kind = ('INDEL_' if self.label == 'INDEL_' else '') + svtype.decode()
        bases = 0
        if svtype in INSERTED:
            seq = _info_field(info, b'SEQ')
            if seq and seq != b'.':
                bases = len(seq) - 2 * seq.count(b'@')
            else:
                bases = _abs_int(_info_field(info, b'SVLEN'))
        elif svtype == b'DEL':
            bases = _abs_int(_info_field(info, b'SVLEN'))
        c = self.counts[(kind, level)]
        c[0] += 1
        c[1] += bases
        self.templates[(kind, level)].setdefault(_template(row_id, info), bases)

    def table(self):
        """[(kind, level, rows, bases, templates, template_bases)] (SNP rows are their own templates)."""
        result = []
        for key in sorted(self.counts, key=_order):
            rows, bases = self.counts[key]
            tpl = self.templates.get(key)
            if tpl is None:
                result.append((*key, rows, bases, rows, bases))
            else:
                result.append((*key, rows, bases, len(tpl), sum(tpl.values())))
        return result


def _abs_int(value):
    try:
        return abs(int(value))
    except (TypeError, ValueError):
        return 0


def _order(key):
    return key[1] != 'top', key[0]


def write_summary(output, tables):
    """TABLES = [(path, RowStats.table())]; writes OUTPUT/cohort.summary.tsv and logs it."""
    lines = [HEADER]
    totals = defaultdict(lambda: [0, 0, 0, 0])
    for path, table in tables:
        for kind, level, *values in table:
            lines.append('\t'.join([os.path.basename(str(path)), kind, level, *map(str, values)]))
            for i, value in enumerate(values):
                totals[(kind, level)][i] += value
    for key in sorted(totals, key=_order):
        lines.append('\t'.join(['ALL', *key, *map(str, totals[key])]))
    target = Path(output) / SUMMARY
    temporary = target.with_name(target.name + '.tmp')
    temporary.write_text('\n'.join(lines) + '\n')
    os.replace(temporary, target)
    print(f'[cohort-merge] summary statistics ({target}):', flush=True)
    width = [max(len(line.split('\t')[i]) for line in lines) for i in range(len(HEADER.split('\t')))]
    for line in lines:
        print('  ' + '  '.join(f'{v:<{w}}' for v, w in zip(line.split('\t'), width)), flush=True)
    return target


def scan(path):
    stats = RowStats(path)
    with open(path, 'rb', buffering=1 << 22) as handle:
        for line in handle:
            stats.add(line.rstrip(b'\r\n'))
    return path, stats.table()


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('output', help='merge output folder with the published cohort.*.vcf')
    args = parser.parse_args(argv)
    output = Path(args.output)
    def rank(path):    # publish order: SV, its singleton file(s), indel, SNP
        name = path.name
        return (0 if name == 'cohort.sv.vcf' else 1 if name.startswith('cohort.sv.') else
                2 if name == 'cohort.indel.vcf' else 3 if name == 'cohort.snp.vcf' else 4), name

    paths = sorted(output.glob('cohort.*.vcf'), key=rank)
    paths = [p for p in paths if not p.name.endswith('.small.vcf')]
    if not paths:
        raise SystemExit(f'no cohort.*.vcf in {output}')
    from concurrent.futures import ProcessPoolExecutor
    with ProcessPoolExecutor(max_workers=len(paths)) as pool:
        tables = list(pool.map(scan, paths))
    write_summary(output, tables)


if __name__ == '__main__':
    sys.exit(main())
