#!/usr/bin/env python3
"""Per-line trio consistency of merged LinGraph VCFs.

A merged row puts every carrier of one allele on the same line, so each row
can be checked directly: for every child haplotype that carries the row
(GT 1), is it also carried by the parent that haplotype came from?

Classes, per child haplotype carrying a row:
  inherited      the expected parent carries it on at least one haplotype
  other_parent   only the other parent carries it (phase swap, or an error)
  child_only     neither parent carries it, and the expected parent is
                 callable (GT 0) on at least one haplotype
  parent_nocall  neither parent carries it and the expected parent has no
                 call (GT .) on both haplotypes: not assessable
and, per child haplotype NOT carrying a row:
  not_transmitted  the expected parent carries it on both haplotypes
                   (homozygous), so the child's haplotype should carry it

By default the child's h1 comes from the father and h2 from the mother
(LinGraph/HPRC naming: h1 paternal, h2 maternal). --unphased accepts either
parent for both haplotypes. Nested rows (CHROM = a parent row's ID) are
checked like top-level rows and marked nested.

Writes one TSV line per row a child haplotype carries (all rows with
--all-rows) to stdout and a summary by class and variant kind to stderr,
with consistent (carried by either parent), inconsistent (by neither) and
inconsistency_rate = inconsistent / (consistent + inconsistent) and
consistency_percent = 100 * consistent / (consistent + inconsistent),
whatever the phasing.

usage: python trio_line_check.py cohort.sv.vcf [cohort.indel.vcf cohort.snp.vcf] \\
           --child NA19240 --mother NA19238 --father NA19239 > lines.tsv
Standard library only.
"""
import argparse
from collections import Counter
import gzip
import re
import sys

NESTED_PREFIXES = ('I_', 'D_', 'S_', 'SUB_', 'DUP_', 'INS_', 'DEL_')
_SVTYPE = re.compile(r'(?:^|;)SVTYPE=([^;]+)')
_SVLEN = re.compile(r'(?:^|;)SVLEN=(-?\d+)')


def open_text(path):
    return gzip.open(path, 'rt') if path.endswith('.gz') else open(path)


def haplotypes(columns, sample):
    """Column indexes of SAMPLE_h1 and SAMPLE_h2 (or SAMPLE itself)."""
    found = {}
    for index, name in enumerate(columns):
        if name == f'{sample}_h1':
            found[1] = index
        elif name == f'{sample}_h2':
            found[2] = index
    if not found and sample in columns:
        found[1] = columns.index(sample)
    if not found:
        raise SystemExit(f'sample {sample!r} not among VCF columns: {", ".join(columns)}')
    return found


def gt(value):
    first = value.split(':', 1)[0]
    return first if first in ('0', '1') else '.'


def kind_of(info, min_sv):
    svtype = _SVTYPE.search(info)
    if svtype is None:
        return 'SNP', 1
    length = _SVLEN.search(info)
    size = abs(int(length.group(1))) if length else 0
    kind = svtype.group(1)
    return (f'{kind}_SV' if size >= min_sv else f'{kind}_indel'), size


def classify(carries, expected, other):
    """Class of one child haplotype for one row."""
    expected_carries = '1' in expected.values()
    other_carries = '1' in other.values()
    if carries == '1':
        if expected_carries:
            return 'inherited'
        if other_carries:
            return 'other_parent'
        if '0' in expected.values():
            return 'child_only'
        return 'parent_nocall'
    if carries == '0' and expected and all(value == '1' for value in expected.values()):
        return 'not_transmitted'
    return None


def main():
    parser = argparse.ArgumentParser(
        description='Per-line trio consistency of merged LinGraph VCFs (see module docstring).')
    parser.add_argument('vcf', nargs='+', help='merged VCFs (sv, indel, snp; .gz accepted)')
    parser.add_argument('--child', required=True)
    parser.add_argument('--mother', required=True)
    parser.add_argument('--father', required=True)
    parser.add_argument('--unphased', action='store_true',
                        help='accept either parent for both child haplotypes')
    parser.add_argument('--min-sv', type=int, default=50,
                        help='size at or above which an INS/DEL/SUB counts as SV (default: 50)')
    parser.add_argument('--min-size', type=int, default=0,
                        help='skip rows smaller than this (SNPs have size 1; default: 0, all rows)')
    parser.add_argument('--no-nested', action='store_true', help='skip nested rows')
    parser.add_argument('--all-rows', action='store_true',
                        help='write every row, not only rows a child haplotype carries')
    args = parser.parse_args()

    out = sys.stdout
    out.write('file\tchrom\tpos\tid\tkind\tsize\tnested\t'
              'child_h1\tchild_h2\tfather_h1\tfather_h2\tmother_h1\tmother_h2\t'
              'class_h1\tclass_h2\n')
    summary = Counter()
    rows = Counter()
    for path in args.vcf:
        with open_text(path) as handle:
            child = father = mother = None
            for line in handle:
                if line.startswith('##'):
                    continue
                fields = line.rstrip('\n').split('\t')
                if line.startswith('#'):
                    columns = fields[9:]
                    child = haplotypes(columns, args.child)
                    father = haplotypes(columns, args.father)
                    mother = haplotypes(columns, args.mother)
                    continue
                nested = fields[0].startswith(NESTED_PREFIXES)
                if nested and args.no_nested:
                    continue
                kind, size = kind_of(fields[7], args.min_sv)
                if size < args.min_size:
                    continue
                samples = fields[9:]
                c = {hap: gt(samples[i]) for hap, i in child.items()}
                f = {hap: gt(samples[i]) for hap, i in father.items()}
                m = {hap: gt(samples[i]) for hap, i in mother.items()}
                classes = {}
                for hap, carries in c.items():
                    if args.unphased:
                        both = {('f', k): v for k, v in f.items()}
                        both.update({('m', k): v for k, v in m.items()})
                        classes[hap] = classify(carries, both, {})
                    else:
                        expected, other = (f, m) if hap == 1 else (m, f)
                        classes[hap] = classify(carries, expected, other)
                group = kind + (' nested' if nested else '')
                rows[group] += 1
                for hap, value in classes.items():
                    if value is not None:
                        summary[(group, f'h{hap}', value)] += 1
                if args.all_rows or '1' in c.values() or 'not_transmitted' in classes.values():
                    out.write('\t'.join([
                        path, fields[0], fields[1], fields[2], kind, str(size),
                        'yes' if nested else 'no',
                        c.get(1, '-'), c.get(2, '-'), f.get(1, '-'), f.get(2, '-'),
                        m.get(1, '-'), m.get(2, '-'),
                        classes.get(1) or '-', classes.get(2) or '-',
                    ]) + '\n')

    err = sys.stderr
    err.write('kind\thaplotype\trows\tconsistent\tinconsistent\ttotal\tinconsistency_rate\tconsistency_percent\n')
    for group in sorted(rows):
        for hap in ('h1', 'h2'):
            # Either parent: Strand-seq-phased hap1/hap2 are not parent-of-origin.
            in_parents = summary[(group, hap, 'inherited')] + summary[(group, hap, 'other_parent')]
            child_only = summary[(group, hap, 'child_only')]
            total = in_parents + child_only
            rate = f'{child_only / total:.4f}' if total else 'NA'
            percent = f'{100 * in_parents / total:.2f}' if total else 'NA'
            err.write(f'{group}\t{hap}\t{rows[group]}\t{in_parents}\t{child_only}\t{total}\t{rate}\t{percent}\n')
    err.write('consistent: carried by either parent; inconsistent: by neither (a parent callable); '
              'inconsistency_rate = inconsistent / total; consistency_percent = 100 * consistent / total\n')

if __name__ == '__main__':
    main()
