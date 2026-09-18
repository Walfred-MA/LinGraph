"""CLI contracts for combined and per-haplotype trio reports."""
from pathlib import Path
import subprocess
import sys
import tempfile
import unittest

SCRIPT = Path(__file__).resolve().parents[1] / 'benchmark' / 'QuickTriocheck.py'


class ChildHaplotypeSummaryTests(unittest.TestCase):
    def setUp(self):
        self.temp = tempfile.TemporaryDirectory()
        self.addCleanup(self.temp.cleanup)
        self.root = Path(self.temp.name)
        self.paths = {}
        self.write('child_h1', [self.sv(5001, 'edge'), self.sv(20001, 'maternal'),
                              self.sv(30001, 'both', 'INS', 80),
                              self.sv(40001, 'childonly', size=-100),
                              self.sv(70001, 'mother_missing')])
        self.write('child_h2', [self.sv(160001, 'paternal', size=-90),
                              self.sv(120001, 'father_missing'),
                              self.sv(90001, 'both_missing')])
        self.write('mother_h1', [self.sv(20003, 'mom'), self.sv(30003, 'mom_both', 'INS', 80)])
        self.write('mother_h2', [], [(0, 60000), (100000, 200000)])
        self.write('father_h1', [self.sv(160003, 'dad', size=-90), self.sv(30003, 'dad_both', 'INS', 80)])
        self.write('father_h2', [], [(0, 90000), (140000, 200000)])

    @staticmethod
    def sv(pos, identifier, kind='DEL', size=-60, original=False):
        end = pos + abs(size) if kind == 'DEL' else pos
        fmt, sample = 'GT', '1'
        if original:
            fmt = 'GT:TYPE:SIZE:EXTENDGRAPHCIGAR'
            sample = f'1:{kind}:{size}:>{abs(size)}' + ('D' if kind == 'DEL' else 'I')
        return f'chr1\t{pos}\t{identifier}\tN\t<{kind}>\t.\tPASS\tSVTYPE={kind};END={end};SVLEN={size}\t{fmt}\t{sample}\n'

    def write(self, name, records, coverage=None):
        path = self.root / f'{name}.vcf'
        path.parent.mkdir(parents=True, exist_ok=True)
        coverage = [(0, 200000)] if coverage is None else coverage
        header = '##fileformat=VCFv4.2\n##contig=<ID=chr1,length=200000>\n'
        header += '##referenceCoverageCoordinateSystem=0-based-half-open\n'
        header += ''.join(f'##referenceCoverage=<Sample="{name}",Chrom="chr1",Start={start},End={end}>\n'
                          for start, end in coverage)
        header += f'#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t{name}\n'
        path.write_text(header + ''.join(records))
        self.paths[name] = path
        return path

    def run_cli(self, *extra, children=None):
        children = children or [self.paths['child_h1'], self.paths['child_h2']]
        command = [sys.executable, '-B', str(SCRIPT), '--child', *map(str, children)]
        for role in ('mother', 'father'):
            command += ['--' + role, str(self.paths[role + '_h1']), str(self.paths[role + '_h2'])]
        result = subprocess.run([*command, *map(str, extra)], capture_output=True, text=True)
        self.assertEqual(result.returncode, 0, result.stderr)
        values = dict(line.split('\t', 1) for line in result.stdout.splitlines()
                      if line.count('\t') == 1)
        return values, result.stdout

    def test_haplotype_counts_filters_fractions_and_outputs(self):
        values, output = self.run_cli('--fp', self.root / 'fp.vcf', '--tp', self.root / 'tp.vcf')
        expected = {
            'child_h1_vcf': str(self.paths['child_h1']),
            'child_h2_vcf': str(self.paths['child_h2']),
            'child_h1_vcf_records_scanned': '5', 'child_h2_vcf_records_scanned': '3',
            'child_h1_records_excluded_edge_trim': '1', 'child_h2_records_excluded_edge_trim': '0',
            'child_h1_records_excluded_parent_missing': '1', 'child_h2_records_excluded_parent_missing': '2',
            'child_h1_records_excluded_mother_missing': '1', 'child_h2_records_excluded_mother_missing': '1',
            'child_h1_records_excluded_father_missing': '0', 'child_h2_records_excluded_father_missing': '2',
            'child_h1_records_excluded_both_parents_missing': '0', 'child_h2_records_excluded_both_parents_missing': '1',
            'child_h1_total_sv': '3', 'child_h2_total_sv': '1', 'child_total_sv': '4',
            'child_h1_found_in_mother': '2', 'child_h2_found_in_mother': '0',
            'child_h1_found_in_father': '1', 'child_h2_found_in_father': '1',
            'child_h1_found_in_both_parents': '1', 'child_h2_found_in_both_parents': '0',
            'child_h1_found_in_either_parent': '2', 'child_h2_found_in_either_parent': '1',
            'child_h1_not_found_in_parents': '1', 'child_h2_not_found_in_parents': '0',
            'child_h1_found_in_either_parent_fraction': '0.666667',
            'child_h2_found_in_either_parent_fraction': '1.000000',
            'child_found_in_either_parent_fraction': '0.750000',
        }
        for key, expected_value in expected.items():
            with self.subTest(key=key):
                self.assertEqual(values[key], expected_value)
        for suffix in ('vcf_records_scanned', 'records_loaded', 'records_excluded_edge_trim',
                       'records_excluded_parent_missing', 'records_excluded_mother_missing',
                       'records_excluded_father_missing', 'records_excluded_both_parents_missing',
                       'total_sv', 'found_in_mother', 'found_in_father', 'found_in_either_parent',
                       'found_in_both_parents', 'not_found_in_parents'):
            self.assertEqual(int(values['child_' + suffix]),
                             int(values['child_h1_' + suffix]) + int(values['child_h2_' + suffix]))
        self.assertIn('h1\tDEL\t2\t1\t0.500000', output)
        self.assertIn('h1\tINS\t1\t1\t1.000000', output)
        self.assertIn('h2\tDEL\t1\t1\t1.000000', output)
        for name, expected_ids in {'fp.child_h1.vcf': ['childonly'], 'fp.child_h2.vcf': [],
                                   'tp.child_h1.vcf': ['maternal', 'both'],
                                   'tp.child_h2.vcf': ['paternal']}.items():
            rows = [line.split('\t')[2] for line in (self.root / name).read_text().splitlines()
                    if not line.startswith('#')]
            self.assertEqual(rows, expected_ids)

    def test_empty_and_fully_filtered_haplotype_report_na(self):
        for rows in ([], [self.sv(90001, 'both_missing')], [self.sv(5001, 'edge')]):
            with self.subTest(rows=rows):
                self.write('child_h2', rows)
                values, _ = self.run_cli()
                self.assertEqual(values['child_h2_total_sv'], '0')
                self.assertEqual(values['child_h2_found_in_either_parent'], '0')
                self.assertEqual(values['child_h2_found_in_either_parent_fraction'], 'NA')
                self.assertEqual(values['child_h2_not_found_in_parents_fraction'], 'NA')
                self.assertEqual(values['child_h1_total_sv'], '3')

    def test_input_order_and_distinct_paths_with_same_filename(self):
        first = self.write('first/child', [self.sv(20001, 'maternal')])
        second = self.write('second/child', [self.sv(40001, 'childonly'), self.sv(41001, 'childonly2')])
        empty = self.write('third/child', [])
        values, _ = self.run_cli(children=[second, first, empty])
        self.assertEqual(values['child_h1_vcf'], str(second))
        self.assertEqual(values['child_h2_vcf'], str(first))
        self.assertEqual(values['child_h1_total_sv'], '2')
        self.assertEqual(values['child_h1_not_found_in_parents'], '2')
        self.assertEqual(values['child_h2_total_sv'], '1')
        self.assertEqual(values['child_h2_found_in_either_parent'], '1')
        self.assertEqual(values['child_h3_total_sv'], '0')

    def test_all_mode_reports_variant_classes(self):
        snp = 'chr1\t20021\tsnp\tC\tT\t.\tPASS\t.\tGT\t1\n'
        indel = 'chr1\t40021\tindel\tA\tAT\t.\tPASS\t.\tGT\t1\n'
        self.write('child_h1', [snp, indel])
        self.write('child_h2', [self.sv(160001, 'paternal', size=-90)])
        self.write('mother_h1', [snp])
        values, output = self.run_cli('--all')
        self.assertEqual(values['child_h1_total_variants'], '2')
        self.assertEqual(values['child_h2_total_variants'], '1')
        self.assertEqual(values['child_total_variants'], '3')
        self.assertNotIn('child_h1_total_sv', values)
        self.assertIn('h1\tSNP\t1\t1\t1.000000', output)
        self.assertIn('h1\tINDEL\t1\t0\t0.000000', output)
        self.assertIn('h2\tSV\t1\t1\t1.000000', output)

    def test_missing_original_coordinates_are_counted_per_source(self):
        self.write('child_h1', [self.sv(20001, 'maternal', original=True), self.sv(40001, 'unknown')])
        self.write('child_h2', [self.sv(160001, 'paternal', original=True),
                              self.sv(40001, 'unknown'), self.sv(41001, 'unknown2')])
        self.write('mother_h1', [self.sv(20001, 'mom', original=True)])
        self.write('father_h1', [self.sv(160001, 'dad', original=True)])
        values, _ = self.run_cli('--original-distance', 'exact')
        self.assertEqual(values['child_h1_records_excluded_missing_original_coordinates'], '1')
        self.assertEqual(values['child_h2_records_excluded_missing_original_coordinates'], '2')
        self.assertEqual(values['child_records_excluded_missing_original_coordinates'], '3')
        self.assertEqual(values['child_h1_found_in_either_parent'], '1')
        self.assertEqual(values['child_h2_found_in_either_parent'], '1')

    def test_single_child_and_multisample_output_keep_existing_keys(self):
        values, output = self.run_cli(children=[self.paths['child_h1']])
        self.assertEqual(values['child_total_sv'], '3')
        self.assertNotIn('child_h1_total_sv', values)
        self.assertNotIn('child_haplotype\t', output)
        merged = self.root / 'merged.vcf'
        header = '##fileformat=VCFv4.2\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tCHILD\tMOTHER\tFATHER\n'
        merged.write_text(header + self.sv(20001, 'maternal').rstrip() + '\t1\t0\n')
        result = subprocess.run([sys.executable, '-B', str(SCRIPT), '-i', str(merged),
                                 '--child', 'CHILD', '--mother', 'MOTHER', '--father', 'FATHER'],
                                capture_output=True, text=True)
        self.assertEqual(result.returncode, 0, result.stderr)
        self.assertIn('child_found_in_either_parent\t1\n', result.stdout)
        self.assertNotIn('child_h1_', result.stdout)


if __name__ == '__main__':
    unittest.main()
