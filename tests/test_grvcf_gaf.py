"""GrvcfGraph GAFs: per-sample grVCFs -> exact merge -> C++ GFA + GAF, lossless
against the assemblies, with the Python walker's statistics (gfa_sample_gaf.py)."""
import csv
from pathlib import Path
import random
import subprocess
import sys

import pytest

REPO = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(REPO / 'tests'))
import test_merge_tools as merge_tools  # noqa: E402
from test_grvcf_graph import tool  # noqa: E402,F401  (fixture)


def random_cohort(out, seed):
    """Shared pools of SNPs, deletions and insertions; carriers change the shared
    insertions (SNP, small deletion or insertion inside) so the exact merge nests
    rows; some samples are on the reverse strand."""
    rng = random.Random(seed)
    seq = lambda n: ''.join(rng.choice('ACGT') for _ in range(n))  # noqa: E731
    length = rng.choice([4000, 8000])
    ref = seq(length)
    pool = []
    for site in sorted(rng.sample(range(60, length - 120, 40), k=min(40, (length - 200) // 40))):
        kind = rng.choice(['snp', 'snp', 'del', 'ins', 'ins'])
        if kind == 'snp':
            pool.append(('snp', site, next(b for b in 'ACGT' if b != ref[site])))
        elif kind == 'del':
            pool.append(('del', site, rng.randrange(1, 30)))
        else:
            pool.append(('ins', site, seq(rng.choice([5, 20, 60, 150]))))

    def changed(text):
        choice = rng.random()
        if choice < 0.3 or len(text) < 12:
            return text
        at = rng.randrange(2, len(text) - 2)
        if choice < 0.6:
            return text[:at] + next(b for b in 'ACGT' if b != text[at]) + text[at + 1:]
        if choice < 0.8:
            return text[:at] + text[at + rng.randrange(1, min(8, len(text) - at - 1)):]
        return text[:at] + seq(rng.randrange(1, 8)) + text[at:]

    for sample in (f's{i}' for i in range(1, rng.randrange(3, 9))):
        reverse = rng.random() < 0.3
        pieces, placed, position, h = [], [], 0, 0
        for kind, site, data in (v for v in pool if rng.random() < 0.5):
            data = changed(data) if kind == 'ins' else data
            pieces.append(ref[position:site])
            h += site - position
            position = site
            placed.append((kind, site, data, h))
            if kind == 'snp':
                pieces.append(data)
                h, position = h + 1, site + 1
            elif kind == 'del':
                position = site + data
            else:
                pieces.append(data)
                h += len(data)
        pieces.append(ref[position:])
        hap = ''.join(pieces)
        q0, n = len(merge_tools.HEAD), len(hap)
        rows = []
        for kind, site, data, h in placed:
            point = f'{q0 + n - h}-' if reverse else f'{q0 + h}+'
            if kind == 'snp':
                rows.append(f'c\t{site + 1}\tSNP_c_{site + 1}_{data}\t{ref[site]}\t{data}\t.\tPASS\tNSUP=1\t'
                            f'{merge_tools.SNPF}\t1:SNP:1:{data}:q{sample}:{point}:.:.\n')
            elif kind == 'del':
                rows.append(f'c\t{site}\tDEL_c_{site}_{sample}\t{ref[site - 1]}\t<DEL>\t.\tPASS\t'
                            f'SVTYPE=DEL;END={site + data};SVLEN=-{data};SEQ=.;EXTENDGRAPHCIGAR=>{data}D\t'
                            f'{merge_tools.SV}\t1:DEL:{data}:>{data}D:q{sample}:{point}:0:.:.:.\n')
            else:
                coord = (f'{q0 + n - h - len(data)}-{q0 + n - h}-' if reverse
                         else f'{q0 + h}-{q0 + h + len(data)}+')
                rows.append(f'c\t{site}\tINS_c_{site}_{sample}\t{ref[site - 1]}\t<INS>\t.\tPASS\t'
                            f'SVTYPE=INS;END={site};SVLEN={len(data)};MAXSIZE={len(data)};SEQ={data};'
                            f'EXTENDGRAPHCIGAR=>{len(data)}I\t{merge_tools.SV}\t'
                            f'1:INS:{len(data)}:>{len(data)}I:q{sample}:{coord}:0:.:.:.\n')
        merge_tools.write_sample(out, sample, ref, hap, rows, reverse)
    merge_tools.fasta(out / 'c.fa', 'c', ref)


def run(*command):
    result = subprocess.run([str(value) for value in command], capture_output=True, text=True)
    assert result.returncode == 0, result.stdout + result.stderr
    return result


def stats(path):
    """Walker counters per sample; link_added differs by design (the C++ graph
    has no carrier-junction links, the walks add them)."""
    with open(path) as handle:
        return {row['sample']: {key: value for key, value in row.items()
                                if key not in ('sample', 'link_added') and value != '0'}
                for row in csv.DictReader(handle, delimiter='\t')}


def test_gaf_rejects_query_coordinate_overflow(tmp_path, tool):  # noqa: F811
    ref = 'A' * 40
    merge_tools.fasta(tmp_path / 'c.fa', 'c', ref)
    row = ('c\t2\tSNP_c_2_C\tA\tC\t.\tPASS\tNSUP=1\t'
           f'{merge_tools.SNPF}\t1:SNP:1:C:qs1:18446744073709551616+:.:.\n')
    merge_tools.write_sample(tmp_path, 's1', ref, ref, [row])
    vcf = tmp_path / 's1.vcf'
    result = subprocess.run([str(tool), '-v', str(vcf), '-r', str(tmp_path / 'c.fa'),
                             '--gaf', str(tmp_path / 'gaf'), '-s', str(vcf), '-o', str(tmp_path / 'out.gfa')],
                            capture_output=True, text=True)
    assert result.returncode != 0
    assert 'QUERYCOORD beyond 2^32' in result.stderr


@pytest.mark.parametrize('cohort', ['deletions', 'insertions', 3, 17, 28])
def test_gaf_is_lossless_and_matches_the_python_walker(tmp_path, tool, cohort):  # noqa: F811
    inputs = tmp_path / 'inputs'
    inputs.mkdir()
    if cohort == 'deletions':
        merge_tools.deletion_cohort(inputs)
    elif cohort == 'insertions':
        merge_tools.insertion_cohort(inputs)
    else:
        random_cohort(inputs, cohort)
    queries = inputs / 'queries.txt'
    queries.write_text(''.join(f'{path.stem}\t{path}\n' for path in sorted(inputs.glob('s*.fa'))))
    samples = sorted(inputs.glob('s*.vcf'))
    python = sys.executable
    run(python, REPO / 'tools/merge_grvcfs.py', '-i', *samples, '-O', tmp_path / 'merged',
        '--exact', queries, '--reference-fasta', inputs / 'c.fa')
    merged = [tmp_path / 'merged' / f'cohort.{kind}.vcf' for kind in ('sv', 'indel', 'snp')]

    run(tool, '-v', *merged, '-r', inputs / 'c.fa', '-o', tmp_path / 'cpp.gfa', '--gaf', tmp_path / 'cpp',
        '-s', *samples, '-t', 3)
    report = run(python, REPO / 'tools/check_gaf_lossless.py', '-g', tmp_path / 'cpp.gfa',
                 '-a', *sorted((tmp_path / 'cpp').glob('*.gaf')), '-q', queries, '-s', *samples,
                 '-o', tmp_path / 'check').stdout
    lossless = [line.split('\t')[-1] for line in report.splitlines()[1:] if '\t' in line]
    assert len(lossless) == len(samples) and all(value == 'yes' for value in lossless), report

    run(python, REPO / 'scripts/merged_vcf_to_gfa.py', '-v', *merged, '-q', queries, '-r', inputs / 'c.fa',
        '--all', '-o', tmp_path / 'py.gfa')
    run(python, REPO / 'scripts/gfa_sample_gaf.py', '-g', tmp_path / 'py.gfa', '-v', *merged, '-s', *samples,
        '-q', queries, '-o', tmp_path / 'py', '--add-links', tmp_path / 'py' / 'added.gfa')
    assert stats(tmp_path / 'cpp' / 'gaf_stats.tsv') == stats(tmp_path / 'py' / 'gaf_stats.tsv')
