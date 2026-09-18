#!/usr/bin/env python3
"""Copy a saved Snakefile with the current 32-CPU/64G SNP chromosome rule.

Only that rule changes. Original workflow/configuration, shards, SV jobs and
completed chromosome parts remain in place. Stop the old workflow before
running Snakemake with the returned file; running jobs retain their allocation.
"""
import argparse
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1] / "scripts"


def rule_span(text):
    start = text.index('    rule graphvcfmerge_small_chrom:\n')
    end = text.index('\n    def small_merge_parts(', start)
    return start, end


def upgrade(work):
    work = Path(work).resolve()
    source = work / 'Snakefile'
    text = source.read_text()
    start, end = rule_span(text)
    current = (ROOT / 'cohort_call_snakemake/Snakefile').read_text()
    first, last = rule_span(current)
    target = work / 'Snakefile.snp-ram'
    # Exclusive creation avoids overwriting a user's further edits or a file
    # currently being used by another launcher.
    with target.open('x') as handle:
        handle.write(text[:start] + current[first:last] + text[end:])
    return target


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('work', type=Path, help='existing unfinished-graphcigar directory')
    args = parser.parse_args(argv)
    print(upgrade(args.work))


if __name__ == '__main__':
    main()
