"""GrvcfGraph segments agree with the Python resolver on random exact cohorts."""
from pathlib import Path
import random
import subprocess
import sys

import pytest

REPO = Path(__file__).resolve().parents[1]
SCRIPTS = REPO / 'scripts'
SOURCE = SCRIPTS / 'src' / 'GrvcfGraph' / 'GrvcfGraph.cpp'
COMPLEMENT = str.maketrans('ACGTacgt', 'TGCAtgca')


def fasta(path, records, width=60):
    data, index = bytearray(), []
    for header, sequence in records:
        data.extend(b'>' + header.encode() + b'\n')
        offset = len(data)
        for low in range(0, len(sequence), width):
            data.extend(sequence[low:low+width].encode() + b'\n')
        bases = min(width, len(sequence)) or 1
        index.append(f'{header.split()[0]}\t{len(sequence)}\t{offset}\t{bases}\t{bases+1}\n')
    path.write_bytes(data)
    Path(str(path) + '.fai').write_text(''.join(index))
    return str(path)


@pytest.fixture(scope='module')
def tool(tmp_path_factory):
    built = SCRIPTS / 'src' / 'GrvcfGraph' / 'GrvcfGraph'
    if built.is_file() and built.stat().st_mtime >= SOURCE.stat().st_mtime:
        return built
    output = tmp_path_factory.mktemp('build') / 'GrvcfGraph'
    subprocess.run(['c++', '-O2', '-std=c++17', '-pthread', '-o', str(output), str(SOURCE)], check=True)
    return output


def bases(rng, length, lower=0.2):
    text = ''.join(rng.choice('ACGT') for _ in range(length))
    if rng.random() < lower:
        cut = rng.randrange(length + 1)
        text = text[:cut] + text[cut:].lower()
    return text


def cohort(folder, seed, count=300):
    """Reference, templates and three merged VCFs with nested rows of every kind."""
    rng = random.Random(seed)
    chr1, chr2 = bases(rng, 20000, 1.0), bases(rng, 8000, 1.0)
    reference = fasta(folder / 'ref.fa', [('chr1', chr1), ('chr2', chr2)])
    templates = fasta(folder / 'templates.fa', [
        # Sequence-identical to the backbone: aliases, forward and reverse.
        ('tplA source=REF:chr1:2000-3000+ backbone=REF lift_status=mapped reference=chr1:2000-3000:+',
         chr1[2000:3000]),
        ('tplB source=REF:chr2:1000-1500- backbone=REF lift_status=mapped reference=chr2:1000-1500:-',
         chr2[1000:1500].translate(COMPLEMENT)[::-1]),
        # Their own sequence: unplaced, and lifted onto chr1.
        ('tplC source=Q1:ctg9:0-700+ backbone=REF lift_status=unmapped', bases(rng, 700)),
        ('tplD source=Q1:ctg9:800-1400+ backbone=REF lift_status=mapped reference=chr1:5000-5100:+',
         bases(rng, 600)),
    ])
    # A novel locus exported as its own path, and a lookup-only catalog whose
    # names alias backbone intervals (only the names rows use are resolved).
    alternatives = fasta(folder / 'alternatives.fa', [('novel_1', bases(rng, 900))])
    catalog = fasta(folder / 'catalog.fa', [
        ('loc1 source=REF:chr2:3000-3400+', chr2[3000:3400]),
        ('loc2 source=REF:chr1:12000-12250-', chr1[12000:12250].translate(COMPLEMENT)[::-1]),
        ('unused source=REF:chr1:1-2+', 'N'),
    ])
    lengths = {'chr1': 20000, 'chr2': 8000, 'tplA': 1000, 'tplB': 500, 'tplC': 700,
               'tplD': 600, 'DUP_chr1_7000_7300': 300, 'DUP_tplA_100_400': 300,
               'novel_1': 900, 'loc1': 400, 'loc2': 250}
    nested_roots = {'DUP_chr1_7000_7300', 'DUP_tplA_100_400'}
    kinds = {}
    rows = []
    # Full-locus duplication sites walk the DUP_ template.
    for number in range(3):
        name = f'I_site{number}'
        template = 'DUP_tplA_100_400' if number == 2 else 'DUP_chr1_7000_7300'
        rows.append((name, 'chr1', 9000 + number * 100, f'SVTYPE=INS;SVLEN=300;END={9000 + number * 100};'
                     f'SEQ={bases(rng, 300)};PACLASS=fulllocusdup;EXTENDGRAPHCIGAR=>{template}:300=',
                     'N', '<INS>', rng.random() < 0.8))
        lengths[name] = 300
        kinds[name] = 'site'
    for number in range(count):
        parents = [name for name in lengths if kinds.get(name) != 'deletion' or True]
        parent = rng.choice(parents)
        nested = parent in kinds or parent in nested_roots
        length = lengths[parent]
        passed = rng.random() < 0.85
        if kinds.get(parent) == 'deletion':
            name = f'I_{number}'
            size = rng.randrange(1, 40)
            offset = rng.randrange(0, length + 1)
            rows.append((name, parent, offset + 1, f'SVTYPE=INS;SVLEN={size};END={offset + 1};'
                         f'SEQ={bases(rng, size)}', 'N', '<INS>', passed))
            lengths[name] = size
            kinds[name] = 'insertion'
            continue
        kind = rng.choice(['insertion', 'insertion', 'deletion', 'snp', 'substitution'])
        if kind == 'snp':
            if length < 1:
                continue
            point = rng.randrange(length)
            name = f'S_{number}'
            rows.append((name, parent, point + 1, '.', 'A', rng.choice('ACGT'), passed))
            lengths[name] = 1
            kinds[name] = 'snp'
            continue
        point = rng.randrange(length + 1)
        pos = point + int(nested)
        if kind == 'insertion':
            size = 2500 if rng.random() < 0.03 else rng.randrange(1, 80)
            sequence = bases(rng, size)
            encoded = f'>{size}I{sequence}' if rng.random() < 0.3 else sequence
            name = f'I_{number}'
            rows.append((name, parent, pos, f'SVTYPE=INS;SVLEN={size};END={pos};SEQ={encoded}',
                         'N', '<INS>', passed))
            lengths[name] = size
        elif kind == 'deletion':
            size = rng.randrange(1, max(2, min(60, length - point + 1)))
            if point + size > length:
                continue
            name = f'D_{number}'
            rows.append((name, parent, pos, f'SVTYPE=DEL;SVLEN=-{size};END={pos + size}',
                         'N', '<DEL>', passed))
            lengths[name] = 0
        else:
            size = rng.randrange(1, max(2, min(30, length - point + 1)))
            if point + size > length:
                continue
            sequence = bases(rng, rng.randrange(1, 30))
            name = f'U_{number}'
            rows.append((name, parent, pos, f'SVTYPE=SUB;END={point + size};SEQ={sequence}',
                         'N', '<SUB>', passed))
            lengths[name] = len(sequence)
        kinds[name] = kind
    rng.shuffle(rows)
    header = ('##fileformat=VCFv4.2\n##contig=<ID=chr1,length=20000>\n##contig=<ID=chr2,length=8000>\n'
              '#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tQ1\n')
    files = [folder / name for name in ('cohort.sv.vcf', 'cohort.indel.vcf', 'cohort.snp.vcf')]
    texts = ['', '', '']
    for name, chrom, pos, info, ref, alt, passed in rows:
        line = (f'{chrom}\t{pos}\t{name}\t{ref}\t{alt}\t.\t{"PASS" if passed else "LowQual"}\t{info}'
                f'\tGT:X\t1:{"x" * rng.randrange(0, 50)}\n')
        texts[rng.randrange(3)] += line
    for path, text in zip(files, texts):
        path.write_text(header + text)
    return reference, templates, files, alternatives, catalog


@pytest.mark.parametrize('seed', [1, 2, 3, 5, 8, 13])
@pytest.mark.parametrize('max_node', [1024, 7])
def test_segments_match_python_resolver(tmp_path, tool, seed, max_node):
    reference, templates, files, alternatives, catalog = cohort(tmp_path, seed)
    gfa, placements, roots = tmp_path / 'out.gfa', tmp_path / 'placements.tsv', tmp_path / 'roots.tsv'
    common = ['-r', reference, '-a', alternatives, '--local-reference-templates', templates,
              '--local-path-fasta', catalog, '--reference-haplotype', 'REF']
    result = subprocess.run([str(tool), '-v', *map(str, files), *common, '-o', str(gfa), '-t', '3',
                             '--max-node-length', str(max_node), '--placements', str(placements),
                             '--roots-table', str(roots)], capture_output=True, text=True)
    assert result.returncode == 0, result.stderr
    sys.path.insert(0, str(REPO / 'tools'))
    import check_grvcf_nodes
    assert check_grvcf_nodes.main(['-g', str(gfa), '-p', str(placements), '--roots-table', str(roots),
                                   '-v', *map(str, files), *common,
                                   '--max-node-length', str(max_node)]) == 0


@pytest.mark.parametrize(('info', 'message'), [
    ('SVTYPE=DEL;SVLEN=-9223372036854775808;END=3', 'SVLEN is outside the supported range'),
    ('SVTYPE=INS;SVLEN=3;END=2;SEQ=>18446744073709551616IACG',
     'graph mapping operation length is out of range'),
])
def test_rejects_integer_overflow(tmp_path, tool, info, message):
    reference = fasta(tmp_path / 'ref.fa', [('c', 'A' * 10)])
    vcf = tmp_path / 'bad.vcf'
    vcf.write_text('##fileformat=VCFv4.2\n'
                   '#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n'
                   f'c\t2\tbad\tA\t<SV>\t.\tPASS\t{info}\n')
    result = subprocess.run([str(tool), '-v', str(vcf), '-r', reference, '-o', str(tmp_path / 'out.gfa')],
                            capture_output=True, text=True)
    assert result.returncode != 0
    assert message in result.stderr


sys.path.insert(0, str(SCRIPTS))



def branch_cohort(folder, hg38):
    """Path-resolution branches the random cohorts do not reach."""
    rng = random.Random(99)
    backbone = 'HG38_h1' if hg38 else 'GRAPH_h1'
    chr1, chr3 = bases(rng, 6000, 1.0), bases(rng, 3000, 1.0)
    t1 = bases(rng, 300)
    t3 = chr3[100:400].translate(COMPLEMENT)[::-1]
    # Graph-mode style backbone FASTA: source= headers name the backbone (inferred,
    # no --reference-haplotype) and hold an alternative from another sample and
    # an imported alternative on an HG38 non-primary contig.
    reference = fasta(folder / 'ref.fa', [
        (f'chr1 source={backbone}:chr1:0-6000+', chr1),
        (f'HG38#1#chr3 source={backbone}:HG38#1#chr3:0-3000+', chr3),
        (f'chrUn_x source={backbone}:chrUn_x:0-500+', bases(rng, 500, 1.0)),
        ('altA source=S2_h1:ctg5:100-2100+', bases(rng, 2000)),
        ('imported source=HG38_h1:chr1_KI270706v1_random:0-400+ sequence_role=imported_alternative',
         bases(rng, 400)),
    ])
    templates = fasta(folder / 'templates.fa', [
        # source_haplotype-style header and one-sided lifts in both directions;
        # t2 is a slice of t1's own path (a template matching an earlier one).
        (f't1 source_haplotype=S3_h2 source_contig=c9 source_start=0 source_end=300 source_strand=+ '
         f'backbone={backbone} lift_status=one_sided lift_note=one_side_up reference=chr1:1000-1000:+', t1),
        (f't2 source=S3_h2:c9:50-250+ backbone={backbone} lift_status=one_sided lift_note=one_side_down '
         f'reference=chr1:2000-2000:-', t1[50:250]),
        (f't3 source={backbone}:HG38#1#chr3:100-400- backbone={backbone} lift_status=both', t3),
    ])
    catalog = fasta(folder / 'catalog.fa', [
        (f'chr1copy source={backbone}:chr1:1000-1500+', chr1[1000:1500]),
        ('t3 source=x:y:0-300+', t3),                       # a name that is already a path
    ])
    lengths = {'chr1': 6000, 'HG38#1#chr3': 3000, 'altA': 2000, 't1': 300, 't2': 200, 't3': 300,
               'imported': 400, 'DUP_HG38_1_chr3_10_200': 190, 'chr1copy': 500}
    return reference, templates, catalog, lengths


@pytest.mark.parametrize('hg38', [False, True])
def test_path_resolution_branches(tmp_path, tool, hg38):
    reference, templates, catalog, lengths = branch_cohort(tmp_path, hg38)
    rng = random.Random(5)
    lines = []
    for number, (parent, length) in enumerate((parent, length) for parent, length in lengths.items()
                                               for _ in range(4)):
        nested = parent.startswith('DUP_')
        point, size = rng.randrange(length + 1), rng.randrange(1, 20)
        lines.append(f'{parent}\t{point + int(nested)}\tI_{number}\tN\t<INS>\t.\tPASS\t'
                     f'SVTYPE=INS;SVLEN={size};END={point + int(nested)};SEQ={bases(rng, size)}\n')
    vcf = tmp_path / 'cohort.sv.vcf'
    vcf.write_text('##fileformat=VCFv4.2\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n' + ''.join(lines))
    gfa, placements, roots = tmp_path / 'out.gfa', tmp_path / 'placements.tsv', tmp_path / 'roots.tsv'
    common = ['-r', reference, '--local-reference-templates', templates, '--local-path-fasta', catalog]
    result = subprocess.run([str(tool), '-v', str(vcf), *common, '-o', str(gfa),
                             '--placements', str(placements), '--roots-table', str(roots)],
                            capture_output=True, text=True)
    assert result.returncode == 0, result.stderr
    names = {line.split('\t')[0] for line in roots.read_text().splitlines()}
    assert ('chrUn_x' in names) != hg38      # HG38_h1 keeps primary contigs only
    assert {'imported', 'HG38#1#chr3', 'DUP_HG38_1_chr3_10_200'} <= names
    sys.path.insert(0, str(REPO / 'tools'))
    import check_grvcf_nodes
    assert check_grvcf_nodes.main(['-g', str(gfa), '-p', str(placements), '--roots-table', str(roots),
                                   '-v', str(vcf), *common]) == 0
