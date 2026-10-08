"""tools/merge_grvcfs.py and tools/convert_merged_grvcf.py: merged grVCFs split
back into individual grVCFs that rebuild every sample's assembly exactly."""
import gzip
from pathlib import Path
import random
import subprocess
import sys

import pytest

ROOT = Path(__file__).resolve().parents[1]
TOOLS, CHECKER = ROOT / 'tools', ROOT / 'tools' / 'check_vcf_lossless.py'
SV = ('GT:TYPE:SIZE:EXTENDGRAPHCIGAR:ASSEMBLYCONTIG:QUERYCOORD:TEMPLATEOFFSET:'
      'LIFTCATEGORY:ALLELENAME:LABEL_H')
SNPF = 'GT:TYPE:SIZE:BASE:ASSEMBLYCONTIG:QUERYCOORD:ALLELENAME:LABEL_H'
HEAD, TAIL = 'A' * 10, 'C' * 10


def revcomp(text):
    return text[::-1].translate(str.maketrans('ACGT', 'TGCA'))


def fasta(path, name, sequence):
    path.write_text(f'>{name}\n{sequence}\n')
    Path(str(path) + '.fai').write_text(
        f'{name}\t{len(sequence)}\t{len(name) + 2}\t{len(sequence)}\t{len(sequence) + 1}\n')


def write_sample(out, sample, ref, hap, rows, reverse=False):
    """Assembly HEAD + haplotype (reverse-complemented on '-') + TAIL, and its grVCF."""
    contig, q0 = 'q' + sample, len(HEAD)
    fasta(out / f'{sample}.fa', contig, HEAD + (revcomp(hap) if reverse else hap) + TAIL)
    (out / f'{sample}.vcf').write_text(
        f'##fileformat=VCFv4.2\n##contig=<ID=c,length={len(ref)}>\n'
        f'##referenceCoverage=<Sample="{sample}",Chrom="c",Start=0,End={len(ref)}>\n'
        '##pseudoLinearMappingCoordinateSystem=0-based-half-open\n'
        f'##pseudoLinearMapping=<ID="0",Sample="{sample}",Query="{contig}:0-{q0}+",'
        'Reference=".",Category="UNMAPPED",Path="pri">\n'
        f'##pseudoLinearMapping=<ID="1",Sample="{sample}",Query="{contig}:{q0}-{q0 + len(hap)}+",'
        f'Reference="c:0-{len(ref)}{"-" if reverse else "+"}",Category="PRIMARY",Path="pri">\n'
        f'##pseudoLinearMapping=<ID="2",Sample="{sample}",'
        f'Query="{contig}:{q0 + len(hap)}-{q0 + len(hap) + len(TAIL)}+",'
        'Reference=".",Category="UNMAPPED",Path="pri">\n'
        f'#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t{sample}\n' + ''.join(rows))


def deletion_cohort(out):
    """Deletions of 40-75 bp at staggered starts, and one sample's SNP."""
    rng = random.Random(7)
    ref = ''.join(rng.choice('ACGT') for _ in range(3000))
    for sample, start, size in (('s1', 1000, 60), ('s2', 1000, 70), ('s3', 995, 70),
                                ('s4', 1002, 75), ('s5', 2000, 40), ('s6', 2003, 52)):
        hap = ref[:start] + ref[start + size:]
        rows = [f'c\t{start}\tDEL_c_{start}_{sample}\t{ref[start - 1]}\t<DEL>\t.\tPASS\t'
                f'SVTYPE=DEL;END={start + size};SVLEN=-{size};SEQ=.;EXTENDGRAPHCIGAR=>{size}D\t{SV}\t'
                f'1:DEL:{size}:>{size}D:q{sample}:{len(HEAD) + start}+:0:.:.:.\n']
        if sample == 's1':
            alt = 'A' if ref[2500] != 'A' else 'C'
            hap = hap[:2500 - size] + alt + hap[2500 - size + 1:]
            rows.append(f'c\t2501\tSNP_c_2501_{alt}\t{ref[2500]}\t{alt}\t.\tPASS\tNSUP=1\t{SNPF}\t'
                        f'1:SNP:1:{alt}:q{sample}:{len(HEAD) + 2500 - size}+:.:.\n')
        write_sample(out, sample, ref, hap, rows)
    fasta(out / 'c.fa', 'c', ref)
    return ref


def insertion_cohort(out):
    """A shared insertion with internal SNPs, deletion and insertion, members
    a few bases off its breakpoint, one on the reverse strand, and SNPs."""
    rng = random.Random(3)
    seq = lambda n: ''.join(rng.choice('ACGT') for _ in range(n))  # noqa: E731
    ref, template = seq(5000), seq(400)

    def mutate(text, at):
        return text[:at] + next(b for b in 'ACGT' if b != text[at]) + text[at + 1:]
    variants = {
        's1': (1500, template), 's2': (1500, mutate(mutate(template, 50), 300)),
        's3': (1503, template[:100] + template[130:]),
        's4': (1498, template[:200] + seq(25) + template[200:]),
        's5': (1500, mutate(template[:100] + template[130:], 250)),
        's6': (3000, seq(80)), 's7': (1501, mutate(template, 10)),
    }
    for sample, (pos, ins) in variants.items():
        reverse = sample == 's7'
        hap = ref[:pos] + ins + ref[pos:]
        q0 = len(HEAD)
        if reverse:
            coord = f'{q0 + len(hap) - pos - len(ins)}-{q0 + len(hap) - pos}-'
        else:
            coord = f'{q0 + pos}-{q0 + pos + len(ins)}+'
        rows = [f'c\t{pos}\tINS_c_{pos}_{sample}\t{ref[pos - 1]}\t<INS>\t.\tPASS\t'
                f'SVTYPE=INS;END={pos};SVLEN={len(ins)};MAXSIZE={len(ins)};SEQ={ins};'
                f'EXTENDGRAPHCIGAR=>{len(ins)}I\t{SV}\t'
                f'1:INS:{len(ins)}:>{len(ins)}I:q{sample}:{coord}:0:.:.:.\n']
        if sample in ('s2', 's4', 's7'):
            alt = next(b for b in 'ACGT' if b != ref[4000])
            at = 4000 + len(ins)
            hap = hap[:at] + alt + hap[at + 1:]
            point = f'{q0 + len(hap) - at}-' if reverse else f'{q0 + at}+'
            rows.append(f'c\t4001\tSNP_c_4001_{alt}\t{ref[4000]}\t{alt}\t.\tPASS\tNSUP=1\t{SNPF}\t'
                        f'1:SNP:1:{alt}:q{sample}:{point}:.:.\n')
        write_sample(out, sample, ref, hap, rows, reverse)
    fasta(out / 'c.fa', 'c', ref)
    return ref


def run(*command):
    result = subprocess.run([sys.executable, *map(str, command)], capture_output=True, text=True)
    assert result.returncode == 0, result.stdout + result.stderr
    return result


def merge(inputs, out, exact=False):
    extra = ['--exact', inputs / 'queries.txt', '--reference-fasta', inputs / 'c.fa'] if exact else []
    run(TOOLS / 'merge_grvcfs.py', '-i', *sorted(inputs.glob('s*.vcf')), '-O', out, *extra)
    return [out / f'cohort.{kind}.vcf' for kind in ('sv', 'indel', 'snp')]


def lossless(vcf, inputs, tmp):
    sample = vcf.stem
    report = run(CHECKER, '-v', vcf, '-r', inputs / 'c.fa', '-q', inputs / f'{sample}.fa',
                 '-o', tmp / f'{sample}.{vcf.parent.name}.report')
    return 'lossless\tyes' in report.stdout


@pytest.fixture(params=['deletions', 'insertions'])
def cohort(request, tmp_path):
    inputs = tmp_path / 'inputs'
    inputs.mkdir()
    (deletion_cohort if request.param == 'deletions' else insertion_cohort)(inputs)
    (inputs / 'queries.txt').write_text(''.join(
        f'{path.stem}\t{path}\n' for path in sorted(inputs.glob('s*.fa'))))
    assert all(lossless(vcf, inputs, tmp_path) for vcf in inputs.glob('s*.vcf'))
    return inputs


@pytest.mark.parametrize('exact', [False, True], ids=['cigar', 'exact'])
def test_split_rebuilds_every_sample_losslessly(cohort, tmp_path, exact):
    merged = merge(cohort, tmp_path / 'merged', exact)
    version = [line for line in merged[0].read_text().splitlines()
               if line.startswith('##graphvcfmergeVersion=')]
    assert version == [f'##graphvcfmergeVersion={"exact" if exact else "cigar"}']
    with gzip.open(tmp_path / 'merged' / 'cohort.samples.headers.gz', 'rt') as handle:
        headers = handle.read()
    assert '##pseudoLinearMapping=' in headers and '##referenceCoverage=' in headers
    split = tmp_path / 'split'
    run(TOOLS / 'convert_merged_grvcf.py', 'split', '-v', *merged, '-r', cohort / 'c.fa', '-o', split)
    vcfs = sorted(split.glob('s*.vcf'))
    assert [vcf.stem for vcf in vcfs] == sorted(path.stem for path in cohort.glob('s*.vcf'))
    assert all(lossless(vcf, cohort, tmp_path) for vcf in vcfs)
    assert (split / 'vcfs.list').read_text().count('\n') == len(vcfs)


def calls(vcf):
    return [line.split('\t')[:2] + [line.split('\t')[7].split(';')[2]]
            for line in vcf.read_text().splitlines() if not line.startswith('#')]


def test_exact_pieces_join_back_into_one_call(tmp_path):
    inputs = tmp_path / 'inputs'
    inputs.mkdir()
    deletion_cohort(inputs)
    (inputs / 'queries.txt').write_text(''.join(
        f'{path.stem}\t{path}\n' for path in sorted(inputs.glob('s*.fa'))))
    merged = merge(inputs, tmp_path / 'merged', exact=True)
    assert '_F1\t' in merged[1].read_text()      # s3: the row plus flank pieces
    run(TOOLS / 'convert_merged_grvcf.py', 'split', '-v', *merged, '-r', inputs / 'c.fa',
        '-o', tmp_path / 'split')
    assert calls(tmp_path / 'split' / 's3.vcf') == [['c', '995', 'SVLEN=-70']]
    assert calls(tmp_path / 'split' / 's6.vcf') == [['c', '2003', 'SVLEN=-52']]  # its SNP joined


def test_no_realignment_keeps_shifted_members_on_their_own_rows(tmp_path):
    inputs = tmp_path / 'inputs'
    inputs.mkdir()
    deletion_cohort(inputs)
    (inputs / 'queries.txt').write_text(''.join(
        f'{path.stem}\t{path}\n' for path in sorted(inputs.glob('s*.fa'))))
    run(TOOLS / 'merge_grvcfs.py', '-i', *sorted(inputs.glob('s*.vcf')), '-O', tmp_path / 'merged',
        '--exact', inputs / 'queries.txt', '--reference-fasta', inputs / 'c.fa', '--no-realignment')
    merged = [tmp_path / 'merged' / f'cohort.{kind}.vcf' for kind in ('sv', 'indel', 'snp')]
    assert '_F1\t' not in merged[0].read_text() + merged[1].read_text()
    # No member is moved onto another's breakpoint: one row per deletion,
    # e.g. s3's 70 bp deletion 5 bp left of s2's.
    assert [call[:2] for call in calls(merged[0])] == [
        ['c', '995'], ['c', '1000'], ['c', '1000'], ['c', '1002'], ['c', '2000'], ['c', '2003']]
    run(TOOLS / 'convert_merged_grvcf.py', 'split', '-v', *merged, '-r', inputs / 'c.fa',
        '-o', tmp_path / 'split')
    assert calls(tmp_path / 'split' / 's3.vcf') == [['c', '995', 'SVLEN=-70']]
    assert all(lossless(vcf, inputs, tmp_path) for vcf in sorted((tmp_path / 'split').glob('s*.vcf')))
    result = subprocess.run([sys.executable, TOOLS / 'merge_grvcfs.py', '-i', *sorted(inputs.glob('s*.vcf')),
                             '-O', tmp_path / 'cigar', '--no-realignment'], capture_output=True, text=True)
    assert result.returncode and '--no-realignment applies to --exact only' in result.stderr


def test_to_cigar_matches_a_direct_merge_up_to_placement(tmp_path):
    inputs = tmp_path / 'inputs'
    inputs.mkdir()
    insertion_cohort(inputs)
    (inputs / 'queries.txt').write_text(''.join(
        f'{path.stem}\t{path}\n' for path in sorted(inputs.glob('s*.fa'))))
    exact = merge(inputs, tmp_path / 'exact', exact=True)
    direct = merge(inputs, tmp_path / 'direct')
    run(TOOLS / 'convert_merged_grvcf.py', 'to-cigar', '-v', *exact, '-r', inputs / 'c.fa',
        '-O', tmp_path / 'converted')
    converted = [tmp_path / 'converted' / path.name for path in direct]
    assert '##graphvcfmergeVersion=cigar' in converted[0].read_text()
    for mine, theirs in zip(converted, direct):
        ids = lambda path: [line.split('\t')[2] for line in path.read_text().splitlines()  # noqa: E731
                            if not line.startswith('#')]
        assert ids(mine) == ids(theirs)
    assert not (tmp_path / 'converted' / 'tmp').exists()
    run(TOOLS / 'convert_merged_grvcf.py', 'split', '-v', *converted, '-r', inputs / 'c.fa',
        '-o', tmp_path / 'split')
    assert all(lossless(vcf, inputs, tmp_path) for vcf in sorted((tmp_path / 'split').glob('s*.vcf')))


@pytest.mark.parametrize('realignment', [True, False], ids=['realign', 'no-realign'])
def test_to_exact_matches_a_direct_exact_merge(cohort, tmp_path, realignment):
    flags = [] if realignment else ['--no-realignment']
    direct = tmp_path / 'direct'
    run(TOOLS / 'merge_grvcfs.py', '-i', *sorted(cohort.glob('s*.vcf')), '-O', direct,
        '--exact', cohort / 'queries.txt', '--reference-fasta', cohort / 'c.fa', *flags)
    cigar = merge(cohort, tmp_path / 'cigar')
    run(TOOLS / 'convert_merged_grvcf.py', 'to-exact', '-v', *cigar, '-r', cohort / 'c.fa',
        '-q', cohort / 'queries.txt', '-O', tmp_path / 'converted', *flags)
    converted = [tmp_path / 'converted' / path.name for path in cigar]
    assert '##graphvcfmergeVersion=exact' in converted[0].read_text()
    for mine, theirs in zip(converted, [direct / path.name for path in cigar]):
        rows = lambda path: [line.split('\t')[:2] for line in path.read_text().splitlines()  # noqa: E731
                             if not line.startswith('#')]
        assert rows(mine) == rows(theirs)
    run(TOOLS / 'convert_merged_grvcf.py', 'split', '-v', *converted, '-r', cohort / 'c.fa',
        '-o', tmp_path / 'split')
    assert all(lossless(vcf, cohort, tmp_path) for vcf in sorted((tmp_path / 'split').glob('s*.vcf')))


def test_split_without_sidecar_warns_and_needs_a_version(tmp_path):
    inputs = tmp_path / 'inputs'
    inputs.mkdir()
    deletion_cohort(inputs)
    merged = merge(inputs, tmp_path / 'merged')
    (tmp_path / 'merged' / 'cohort.samples.headers.gz').unlink()
    for path in merged:            # a merge from before the version line existed
        path.write_text(''.join(line for line in path.read_text().splitlines(keepends=True)
                                if not line.startswith('##graphvcfmergeVersion=')))
    command = [sys.executable, TOOLS / 'convert_merged_grvcf.py', 'split', '-v', *merged,
               '-r', inputs / 'c.fa', '-o', tmp_path / 'split']
    refused = subprocess.run(command, capture_output=True, text=True)
    assert refused.returncode == 1 and '--input-version' in refused.stderr
    result = run(*command[1:], '--input-version', 'cigar')
    assert 'cohort.samples.headers.gz' in result.stderr and 'warning' in result.stderr
    text = (tmp_path / 'split' / 's3.vcf').read_text()
    assert '##referenceCoverage' not in text and '\tD_c_995_s3\t' not in text
    assert calls(tmp_path / 'split' / 's3.vcf') == [['c', '995', 'SVLEN=-70']]
