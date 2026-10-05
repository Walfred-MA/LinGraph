"""tools/grvcf_to_vcf.py: one streaming pass, rows kept as representative alleles."""
import gzip
from pathlib import Path
import shutil
import subprocess
import sys
import tempfile
import unittest

ROOT = Path(__file__).resolve().parents[1]
SCRIPT = ROOT / 'tools/grvcf_to_vcf.py'
FORMAT = 'GT:TYPE:SIZE:QUERYCOORD'


class GrVcfConversionTests(unittest.TestCase):
    def setUp(self):
        self.tmp = tempfile.TemporaryDirectory()
        self.addCleanup(self.tmp.cleanup)
        self.root = Path(self.tmp.name)
        self.source = self.root / 'input.gr.vcf'
        self.output = self.root / 'normal.vcf'
        self.header = ('##fileformat=VCFv4.2\n##source=graphvcfmerge\n'
                       '##contig=<ID=chr1,length=100>\n'
                       '##contig=<ID=INS_chr1_2_1,length=4>\n'
                       '##contig=<ID=locus7,length=50>\n'
                       '##FILTER=<ID=LowQual,Description="Input filter">\n'
                       '##pseudoLinearMapping=<ID="0",Sample="A_h1",Query="q:0-9+">\n'
                       '##INFO=<ID=PAMAP,Number=.,Type=Integer,Description="x">\n'
                       '#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tA_h1\tB_h1\tC_h1\n')

    def row(self, chrom, pos, identifier, ref, alt, info, samples=('1:INS:3:q:2-5+', '0:.:.:.', '.:.:.:.')):
        return '\t'.join([chrom, str(pos), identifier, ref, alt, '.', 'PASS', info, FORMAT, *samples]) + '\n'

    def convert(self, rows, *extra, header=None, success=True):
        self.source.write_text((header or self.header) + ''.join(rows))
        result = subprocess.run([sys.executable, str(SCRIPT), '-i', str(self.source),
                                 '-o', str(self.output), *extra], capture_output=True, text=True)
        self.assertEqual(result.returncode == 0, success, result.stdout + result.stderr)
        return result

    def lines(self, path=None):
        path = path or self.output
        opener = gzip.open if str(path).endswith('.gz') else open
        with opener(path, 'rt') as handle:
            return handle.read().splitlines()

    def records(self, path=None):
        return [line.split('\t') for line in self.lines(path) if not line.startswith('#')]

    def test_alleles_info_and_genotypes(self):
        result = self.convert([
            self.row('chr1', 2, 'INS_chr1_2_1', 'C', '<INS>',
                     'SVTYPE=INS;END=2;NSUP=1;SVLEN=4;MAXSIZE=4;SEQ=aacG;EXTENDGRAPHCIGAR=>4I;PAMAP=0'),
            self.row('chr1', 6, 'DEL_chr1_6_2', 'A', '<DEL>',
                     'SVTYPE=DEL;END=9;NSUP=1;SVLEN=-3;SEQ=.;PACLASS=primary'),
            self.row('chr1', 20, 'SNP_chr1_20_G', 'T', 'G', 'NSUP=1',
                     samples=('1:SNP:1:q:19+', '1:SNP:1:q:20+', '0:.:.:.')),
        ])
        rows = self.records()
        self.assertEqual(rows[0][:8], ['chr1', '2', 'INS_chr1_2_1', 'C', 'CAACG', '.', 'PASS',
                                       'SVTYPE=INS;SVLEN=4;NSUP=1'])
        self.assertEqual(rows[1][3:8], ['A', '<DEL>', '.', 'PASS', 'SVTYPE=DEL;END=9;SVLEN=-3;NSUP=1'])
        self.assertEqual(rows[2][3:5], ['T', 'G'])
        self.assertEqual([r[8:] for r in rows], [['GT', '1', '0', '.'], ['GT', '1', '0', '.'],
                                                 ['GT', '1', '1', '0']])
        self.assertIn('read=3 written=3 nested_removed=0 insertions_with_bases=1', result.stderr)

    def test_nested_rows_and_their_contigs_are_removed(self):
        result = self.convert([
            self.row('chr1', 2, 'INS_chr1_2_1', 'C', '<INS>', 'SVTYPE=INS;END=2;SVLEN=4;SEQ=AAAA'),
            self.row('INS_chr1_2_1', 2, 'SNP_INS_chr1_2_1_2', 'A', 'G', 'NSUP=1'),
            self.row('INS_chr1_2_1', 1, 'INS_INS_chr1_2_1_1', 'N', '<INS>', 'SVTYPE=INS;END=1;SVLEN=2;SEQ=TT'),
            self.row('DEL_chr1_6_3', 1, 'INS_DEL_chr1_6_3_1', 'N', '<INS>', 'SVTYPE=INS;END=1;SVLEN=1;SEQ=C'),
            self.row('DUP_chr1_1_4', 2, 'SNP_DUP_1', 'A', 'C', 'NSUP=1'),
            self.row('locus7', 3, 'SNP_locus7_3_C', 'A', 'C', 'NSUP=1'),
        ])
        self.assertEqual([r[2] for r in self.records()], ['INS_chr1_2_1', 'SNP_locus7_3_C'])
        header = [line for line in self.lines() if line.startswith('##')]
        self.assertEqual([line for line in header if line.startswith('##contig')],
                         ['##contig=<ID=chr1,length=100>', '##contig=<ID=locus7,length=50>'])
        self.assertFalse(any('pseudoLinearMapping' in line or 'PAMAP' in line for line in header))
        self.assertIn('##FILTER=<ID=LowQual,Description="Input filter">', header)
        self.assertEqual(header[0], '##fileformat=VCFv4.2')
        self.assertIn('read=6 written=2 nested_removed=4', result.stderr)

    def test_encoded_insertion_stays_symbolic_and_position_zero_is_removed(self):
        result = self.convert([
            self.row('chr1', 3, 'INS_enc', 'G', '<INS>', 'SVTYPE=INS;END=3;SVLEN=5;SEQ=.;EXTENDGRAPHCIGAR=>x:5='),
            self.row('chr1', 0, 'INS_zero', 'N', '<INS>', 'SVTYPE=INS;END=0;SVLEN=2;SEQ=AC'),
        ])
        self.assertEqual([r[2:5] + [r[7]] for r in self.records()],
                         [['INS_enc', 'G', '<INS>', 'SVTYPE=INS;END=3;SVLEN=5']])
        self.assertIn('position_zero_removed=1', result.stderr)

    def test_sites_only_input(self):
        header = self.header.rsplit('#CHROM', 1)[0] + '#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n'
        self.convert(['chr1\t20\tsnp\tT\tG\t.\tPASS\tNSUP=1\n'], header=header)
        self.assertEqual(self.records(), [['chr1', '20', 'snp', 'T', 'G', '.', 'PASS', 'NSUP=1']])
        self.assertFalse(any(line.startswith('##FORMAT') for line in self.lines()))

    @unittest.skipUnless(shutil.which('bgzip'), 'needs bgzip')
    def test_bgzip_output_and_force(self):
        self.output = self.root / 'normal.vcf.gz'
        rows = [self.row('chr1', 20, 'snp', 'T', 'G', 'NSUP=1')]
        self.convert(rows)
        self.assertEqual(self.records()[0][2], 'snp')
        self.convert(rows, success=False)            # exists: needs --force
        self.convert(rows, '--force')
        with open(self.output, 'rb') as handle:     # BGZF: gzip with the BC extra field
            self.assertEqual(handle.read(16)[12:14], b'BC')


if __name__ == '__main__':
    unittest.main()
