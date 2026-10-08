#!/usr/bin/env python3
"""Separate merge-introduced losses from per-sample calling losses.

Compares, per haplotype, the run reports of tools/check_vcf_lossless.py (the
per-sample VCF against its assembly) and tools/check_merge_lossless.py (the
merged VCF's projection of that haplotype). Both list one row per mismatching
run, keyed by the run's query interval, plus UNASSIGNED_* rows.

A run that mismatches only in the merge check was exact in the per-sample VCF:
the merge changed it. Those runs are written, with their first mismatch, to
OUTPUT (default: MERGE_CHECK/merge_introduced.tsv). Unassigned merged
observations whose query interval is not unassigned in the per-sample report
go to MERGE_CHECK/merge_unassigned.tsv (the row column names the merged row;
a split replacement can list a part whose interval changed).

usage: compare_lossless_reports.py -m MERGE_CHECK_DIR -c VCF_CHECK_DIR [-o OUTPUT]
  VCF_CHECK_DIR holds SAMPLE.vcf_check.tsv, MERGE_CHECK_DIR SAMPLE.merge_check.tsv.
Standard library only.
"""
import argparse
import csv
from pathlib import Path
import sys


def read_report(path):
    """Mismatching runs by query interval; unassigned rows as a list."""
    runs, unassigned = {}, []
    with open(path, newline='') as handle:
        for row in csv.DictReader(handle, delimiter='\t'):
            if row['status'].startswith('MISMATCH'):
                runs[row['query']] = row
            elif row['status'].startswith('UNASSIGNED'):
                unassigned.append(row)
    return runs, unassigned


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('-m', '--merge-check', required=True,
                        help='folder of SAMPLE.merge_check.tsv (check_merge_lossless.py -o)')
    parser.add_argument('-c', '--vcf-check', required=True,
                        help='folder of SAMPLE.vcf_check.tsv (check_vcf_lossless.py -o)')
    parser.add_argument('-o', '--output', help='merge-introduced runs TSV')
    args = parser.parse_args(argv)
    merge_dir, vcf_dir = Path(args.merge_check), Path(args.vcf_check)
    output = Path(args.output) if args.output else merge_dir / 'merge_introduced.tsv'
    fields = ('sample', 'query', 'reference', 'status', 'first_mismatch',
              'query_length', 'predicted_length', 'expected_context', 'predicted_context')
    print('sample\truns_mismatch_per_sample\truns_mismatch_merged\tboth\t'
          'merge_introduced\tfixed_by_merge\tunassigned_per_sample\tunassigned_merged')
    missing = []
    with open(output, 'w', newline='') as out, \
            open(merge_dir / 'merge_unassigned.tsv', 'w', newline='') as extra:
        writer = csv.writer(out, delimiter='\t', lineterminator='\n')
        writer.writerow(fields)
        extra_writer = csv.writer(extra, delimiter='\t', lineterminator='\n')
        extra_writer.writerow(('sample', 'status', 'query', 'reference', 'query_length', 'row'))
        for merged_report in sorted(merge_dir.glob('*.merge_check.tsv')):
            sample = merged_report.name[:-len('.merge_check.tsv')]
            vcf_report = vcf_dir / f'{sample}.vcf_check.tsv'
            if not vcf_report.is_file():
                missing.append(str(vcf_report))
                continue
            merged, merged_unassigned = read_report(merged_report)
            own, own_unassigned = read_report(vcf_report)
            introduced = sorted(set(merged) - set(own))
            print('\t'.join(map(str, (
                sample, len(own), len(merged), len(set(own) & set(merged)), len(introduced),
                len(set(own) - set(merged)), len(own_unassigned), len(merged_unassigned)))))
            known = {(row['status'], row['query']) for row in own_unassigned}
            for row in merged_unassigned:
                if (row['status'], row['query']) not in known:
                    extra_writer.writerow((sample, row['status'], row['query'], row['reference'],
                                           row['query_length'], row['expected_context']))
            for query in introduced:
                row = merged[query]
                writer.writerow((sample, query, row['reference'], row['status'],
                                 row['first_mismatch'], row['query_length'],
                                 row['predicted_length'], row['expected_context'],
                                 row['predicted_context']))
    if missing:
        print('missing per-sample reports: ' + ', '.join(missing), file=sys.stderr)
    print(f'merge-introduced runs: {output}; new unassigned: '
          f'{merge_dir / "merge_unassigned.tsv"}', file=sys.stderr)
    return 1 if missing else 0


if __name__ == '__main__':
    sys.exit(main())
