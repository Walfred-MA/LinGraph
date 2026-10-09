"""Differential checks against the retained Python small-variant implementation."""
import gzip
import json
from pathlib import Path
import random
import sys

import numpy as np
import pytest

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / 'scripts'))
import graphvcfmerge as vcf
import graphvcfmerge_snp as small
import graphvcfmerge_snp_compact as compact
import graphvcfmerge_native as native


def row(chrom, pos, alt, fields, fmt=None, info='.', filt='PASS', ref='A'):
    return '\t'.join([chrom, str(pos), 'original', ref, alt, '37', filt, info,
                      fmt or vcf.SNP_FORMAT, *fields]) + '\n'


def allele(alt='C', query='ctg', qpos=2, label='0H', name='allele', kind='SNP'):
    return '1:' + ':'.join(vcf.vcf_escape(s) for s in
                          [kind, '1', alt, query, f'{qpos}+', name, label])


def source(path, names, rows, extra=()):
    text = ('##fileformat=VCFv4.2\n##contig=<ID=chr1,length=4294967295>\n'
            '##contig=<ID=chr2,length=1000>\n' + '\n'.join(extra) + '\n'
            '#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t'
            + '\t'.join(names) + '\n' + ''.join(rows))
    if str(path).endswith('.gz'):
        with gzip.open(path, 'wt') as handle:
            handle.write(text)
    else:
        path.write_text(text)
    return path


def run(tmp, inputs, backend, monkeypatch, *, snp_only=True, processes=1):
    monkeypatch.setenv('GRAPHVCFMERGE_BACKEND', backend)
    root = tmp / backend
    manifest = small.scan(inputs, root, processes=processes, snp_only=snp_only)
    for key in manifest['chroms']:
        small.merge_chrom(root, key, processes=processes)
    snp, indel = tmp / f'{backend}.snp.vcf', tmp / f'{backend}.indel.vcf'
    small.concat(root, snp, indel)
    return manifest, (snp.read_bytes(), indel.read_bytes())


@pytest.mark.parametrize('snp_only', [True, False])
def test_exact_outputs_custom_fields_and_coverage(tmp_path, monkeypatch, snp_only):
    rows = [
        row('chr2', 9, 'G', [allele('G', label='-3H&.&12H'), '.']),
        row('chr1', 5, 'C', [allele(query='ctg:one', name='a,;:@\t\n', label='-2147483648H'), '0']),
        row('chr1', 5, 'C', [allele(query='ctg:one', name='a,;:@\t\n', label='-2147483648H'), '.'], filt='q10'),
        row('chr1', 5, 'G', [allele('G', label='2147483647H'), '.'], filt='q20;q10'),
        row('chr1', 11, 'T', ['0', allele('T')]),
        row('chr1', 4294967295, 'N', [allele('N', qpos=4294967295), '.']),
        row('chr1', 7, 'C', [allele(), '.'], info='PAMAP=1,2'),  # annotation
        row('chr1', 8, 'C', [allele(), '.'], info='PAMAP=1,3'),  # real call
        row('chr1', 12, 'C', [allele(), '.'], info='PAPATH=alt'),
    ]
    if not snp_only:
        obs = ['1', 'INS', '3', '>2H<1H>3I', 'ctg', '2-5+', 'a*name', '0H0H', '.', '.']
        field = vcf.format_sample_alleles([vcf.format_hsv_allele(obs)], vcf.SV_FORMAT)
        rows += [row('chr1', 0, '<INS>', [field, '.'], vcf.SV_FORMAT,
                     'SVTYPE=INS;END=0;MAXSIZE=3;SEQ=acg;PAMAP=9;NSUP=9', ref='N'),
                 row('chr1', 10, '<DEL>', [field.replace('INS', 'DEL'), '.'], vcf.SV_FORMAT,
                     'SVTYPE=DEL;END=13;SVLEN=-3', ref='N')]
    extras = ['##referenceCoverage=<Chrom=chr1,Sample=b,Start=0,End=12>',
              '##referenceCoverage=<Chrom=chr1,Sample=b,Start=12,End=20>',
              '##pseudoLinearMapping=<ID=1,Path=alt>',
              '##pseudoLinearMapping=<ID=2,Path="alt">']
    paths = [source(tmp_path / 'one.vcf.gz', ['a', 'b'], rows, extras),
             source(tmp_path / 'two.vcf', ['b', 'c'],
                    [row('chr1', 5, 'C', ['0', allele(name='x')]),
                     row('chr1', 21, 'T', ['0', allele('T')])],
                    ['##referenceCoverage=<Chrom=chr1,Sample=a,Start=0,End=30>'])]
    _, expected = run(tmp_path, paths, 'python', monkeypatch, snp_only=snp_only)
    _, actual = run(tmp_path, paths, 'native', monkeypatch, snp_only=snp_only)
    assert actual == expected


@pytest.mark.parametrize('snp_only', [True, False])
def test_legacy_and_vector_formats(tmp_path, monkeypatch, snp_only):
    vector = vcf.format_sample_alleles([
        '1:SNP:1:C:ctg:10+:z:1H', '1:SNP:1:C:ctg:2+:a:-2H',
        '1:SNP:1:C:ctg:2+:a:-2H'], vcf.SNP_FORMAT)
    legacy = vcf.format_sample_alleles(['1:SNP:1:C:query:8+:x:0H'], 'HSNP')
    paths = [source(tmp_path / 'input.vcf', ['s'], [
        row('chr1', 8, 'C', [vector]), row('chr1', 8, 'C', [legacy], 'HSNP'),
        row('chr1', 9, 'C', ['1:SNP:1:C:query:8+'], vcf.SNP_FORMAT_BASE),
    ])]
    _, expected = run(tmp_path, paths, 'python', monkeypatch, snp_only=snp_only)
    _, actual = run(tmp_path, paths, 'native', monkeypatch, snp_only=snp_only)
    assert actual == expected


def test_saved_compact_shards_realignment_and_spill(tmp_path, monkeypatch):
    # 2,000 output samples, duplicate physical observations, unsorted input,
    # joined labels and a small enough budget to force several disk runs.
    samples = [f's{i}' for i in range(2000)]
    rng = random.Random(42)
    writer = compact.Writer(tmp_path / 'source')
    for _ in range(50000):
        sample = rng.choice(samples)
        pos = rng.randrange(1, 45)
        writer.add('chr1', pos, 'A', 'C', sample,
                   ['1', 'SNP', '1', 'C', 'ctg', f'{pos}+', 'allele', '-3H&1H'])
    writer.add('chr1', 49, 'A', 'C', 's0', ['1', 'SNP', '1', 'C', 'ctg', '50+', '.', '.'])
    writer.coverage['chr1'] = [[s, 0, 20 + i % 20] for i, s in enumerate(samples)]
    data = writer.finish()
    root = tmp_path / 'manifest'
    manifest = compact.manifest_from_sources([dict(source=str(writer.directory / 'source.json'),
                                                   samples=samples, chroms=data['chroms'], metadata=[])],
                                            root, str(tmp_path))
    key = compact.chrom_key('chr1')
    (root / 'realign').mkdir()
    (root / 'realign' / (key + '.json')).write_text(json.dumps(dict(
        chrom='chr1', drops=[['s0', 49, 'C', 'ctg', 50, '+']],
        adds=[['s0', 50, 'A', 'T', 'ctg', 51, '+']])))
    monkeypatch.setenv('GRAPHVCFMERGE_BACKEND', 'python')
    compact.merge_chrom(root, key, 2)
    expected = small.part_path(manifest, key, 'snp').read_bytes()
    monkeypatch.setenv('GRAPHVCFMERGE_BACKEND', 'native')
    compact.merge_chrom(root, key, 2)
    assert small.part_path(manifest, key, 'snp').read_bytes() == expected
    monkeypatch.setenv('GRAPHVCFMERGE_RAM_GB', '0.001')
    compact.merge_chrom(root, key, 2)
    assert small.part_path(manifest, key, 'snp').read_bytes() == expected
    assert not list((tmp_path / 'parts').rglob('*.run.*'))
    # A failed replacement retains the previous valid part.
    realign = root / 'realign' / (key + '.json')
    contents = json.loads(realign.read_text())
    contents['drops'][0][1] = 500
    realign.write_text(json.dumps(contents))
    with pytest.raises(ValueError, match='drop not found'):
        compact.merge_chrom(root, key)
    assert small.part_path(manifest, key, 'snp').read_bytes() == expected


def test_coverage_only_empty_category(tmp_path, monkeypatch):
    paths = [source(tmp_path / 'input.vcf', ['s'], [],
                    ['##referenceCoverage=<Chrom=chr1,Sample=s,Start=0,End=30>'])]
    _, expected = run(tmp_path, paths, 'python', monkeypatch)
    _, actual = run(tmp_path, paths, 'native', monkeypatch)
    assert actual == expected


def test_corrupt_dictionary_does_not_replace_output(tmp_path, monkeypatch):
    paths = [source(tmp_path / 'input.vcf', ['s'], [row('chr1', 5, 'C', [allele()])])]
    monkeypatch.setenv('GRAPHVCFMERGE_BACKEND', 'native')
    root = tmp_path / 'shards'
    manifest = compact.scan(paths, root)
    key = compact.chrom_key('chr1')
    compact.merge_chrom(root, key)
    original = small.part_path(manifest, key, 'snp').read_bytes()
    source_info = json.loads(Path(manifest['sources'][0]).read_text())
    binary = Path(source_info['directory']) / 'records.bin'
    records = np.fromfile(binary, dtype=compact.DTYPE)
    records['query'] = 2**32 - 1
    records.tofile(binary)
    with pytest.raises(ValueError, match='dictionary index'):
        compact.merge_chrom(root, key)
    assert small.part_path(manifest, key, 'snp').read_bytes() == original


def test_ram_defaults(monkeypatch):
    monkeypatch.delenv('GRAPHVCFMERGE_RAM_GB', raising=False)
    monkeypatch.delenv('SLURM_MEM_PER_NODE', raising=False)
    monkeypatch.delenv('SLURM_MEM_PER_CPU', raising=False)
    assert [native.memory_bytes(n) >> 30 for n in (99, 100, 1000, 2000, 10000)] == [64, 128, 128, 128, 512]
    monkeypatch.setenv('SLURM_MEM_PER_NODE', '65536')
    assert native.memory_bytes(2000) == 64 << 30
    assert native.memory_bytes(2000, 4) == 16 << 30


def test_small_indel_spill_keys_and_boundaries(tmp_path, monkeypatch):
    rng = random.Random(13)
    names = ['one', 'two', 'three']
    rows = []
    for i in range(10000):
        pos = rng.randrange(6, 16)
        kind = rng.choice(['INS', 'DEL', 'SUB'])
        size = rng.randrange(1, 4)
        sequence = 'ACG'[:size] if kind != 'DEL' else '.'
        end = pos if kind == 'INS' else pos + size
        observation = ['1', kind, str(size), f'>{size}I', 'ctg', f'{i % 19}+',
                       rng.choice(['a', 'a&b', '*c']), '-1H2H', '.', '.']
        field = vcf.format_sample_alleles([vcf.format_hsv_allele(observation)], vcf.SV_FORMAT)
        fields = [rng.choice(['.', '0']) for _ in names]
        fields[rng.randrange(3)] = field
        rows.append(row('chr雪🚀', pos, '<' + kind + '>', fields, vcf.SV_FORMAT,
                        f'SVTYPE={kind};SVLEN={size};END={end};SEQ={sequence};PAMAP=3;FLAG;NSUP=1',
                        filt='PASS' if i % 4 else 'q20', ref='N'))
    paths = [source(tmp_path / 'indels.vcf', names, rows,
                    [f'##referenceCoverage=<Chrom=chr雪🚀,Sample={s},Start=0,End={10+i}>'
                     for i, s in enumerate(names)])]
    _, expected = run(tmp_path, paths, 'python', monkeypatch, snp_only=False)
    monkeypatch.setenv('GRAPHVCFMERGE_RAM_GB', '0.001')
    _, actual = run(tmp_path, paths, 'native', monkeypatch, snp_only=False)
    assert actual == expected


def test_scan_workers_and_old_shards(tmp_path, monkeypatch):
    paths = [source(tmp_path / f'{i}.vcf', [f's{i}'],
                    [row('chr1', 9, 'C', [allele(query=f'q{i}')]),
                     row('chr2', 6, 'G', [allele('G')]),
                     row('chr1', 3, 'T', [allele('T')])]) for i in range(4)]
    _, expected = run(tmp_path, paths, 'python', monkeypatch, processes=2)
    manifest, actual = run(tmp_path, paths, 'native', monkeypatch, processes=2)
    assert actual == expected
    # Version 2 stores one binary per chromosome, without a section table.
    for source_path in manifest['sources']:
        path = Path(source_path)
        data = json.loads(path.read_text())
        binary = Path(data['directory']) / 'records.bin'
        raw = binary.read_bytes()
        for key, spans in data.pop('sections').items():
            (binary.parent / (key + '.bin')).write_bytes(b''.join(raw[o:o+n*21] for o, n in spans))
        data['protocol'] = 2
        path.write_text(json.dumps(data))
        binary.unlink()
    for key in manifest['chroms']:
        compact.merge_chrom(tmp_path / 'native', key)
    out = tmp_path / 'old-shards.vcf'
    small.concat(tmp_path / 'native', out)
    assert out.read_bytes() == expected[0]


@pytest.mark.parametrize('field,error', [
    ('1:SNP,SNP:1:C:ctg:8+:a:0H', 'observation count'),
    ('1:SNP:1:C:ctg:4294967296+:a:0H', 'uint32'),
    ('1:SNP:1:C:ctg:8+:a:2147483648H', 'signed 32-bit'),
    ('1:SNP:1:C:ctg:8+:a:+-1H', 'integer sign'),
])
def test_invalid_scan_reports_source_line(tmp_path, monkeypatch, field, error):
    path = source(tmp_path / 'bad.vcf', ['s'], [row('chr1', 8, 'C', [field])])
    monkeypatch.setenv('GRAPHVCFMERGE_BACKEND', 'native')
    with pytest.raises(ValueError, match=error) as failure:
        compact.scan([path], tmp_path / 'bad')
    assert 'bad.vcf:' in str(failure.value)


def test_truncated_gzip_is_not_published(tmp_path, monkeypatch):
    path = source(tmp_path / 'bad.vcf.gz', ['s'], [row('chr1', 8, 'C', [allele()])])
    path.write_bytes(path.read_bytes()[:-8])
    monkeypatch.setenv('GRAPHVCFMERGE_BACKEND', 'native')
    root = tmp_path / 'bad'
    with pytest.raises(ValueError, match='compressed input'):
        compact.scan([path], root)
    assert not (root / 'manifest.json').exists()


@pytest.mark.parametrize('snp_only', [True, False])
def test_header_trailing_whitespace(tmp_path, monkeypatch, snp_only):
    paths = [source(tmp_path / 'one.vcf', ['a'], [row('chr1', 5, 'C', [allele()])]),
             source(tmp_path / 'two.vcf', ['b'], [row('chr1', 9, 'C', [allele()], info='PAMAP=1')],
                    ['##referenceCoverage=<Chrom=chr1,Sample=b,Start=0,End=20>  ',
                     '##pseudoLinearMapping=<ID=1,Path=alt>  '])]
    _, expected = run(tmp_path, paths, 'python', monkeypatch, snp_only=snp_only)
    _, actual = run(tmp_path, paths, 'native', monkeypatch, snp_only=snp_only)
    assert actual == expected


def test_long_joined_indel_labels(tmp_path, monkeypatch):
    # Label parsing must be linear and must not recurse once per joined atom.
    label = '&'.join(['-1H_+2H', '.', '+3H-4H'] * 10000)
    values = ['1', 'INS', '3', '>2H<1H>3I', 'ctg', '2-5+', 'allele', label, '.', '.']
    field = vcf.format_sample_alleles([vcf.format_hsv_allele(values)], vcf.SV_FORMAT)
    paths = [source(tmp_path / 'input.vcf', ['a'], [
        row('chr1', 5, '<INS>', [field], vcf.SV_FORMAT, 'SVTYPE=INS;MAXSIZE=3;SEQ=ACG')])]
    _, expected = run(tmp_path, paths, 'python', monkeypatch, snp_only=False)
    _, actual = run(tmp_path, paths, 'native', monkeypatch, snp_only=False)
    assert actual == expected


def test_indel_label_suffix_compatibility(tmp_path, monkeypatch):
    # The legacy parser uses these suffixes to locate CIGAR's right boundary.
    # Include malformed suffixes: they must follow the same fallback as Python.
    labels = ['.', '0H', '-1H_+2H', '+3H-4H', '.&-2H0H&1H', '', '&', '.&',
              '1H&', '&1H', '1H_', '1H__2H', '1H2H3H', '+-1H', '1H&&2H', '1h', 'x']
    offsets = ['.', '0', '-1', '+2', '', '+', '+-1', 'x']
    rows = []
    for i, (label, offset) in enumerate((a, b) for a in labels for b in offsets):
        field = f'1:INS:3:>2H<1H>3I:ctg:2-5+:{offset}:Pri:a:{label}'
        rows.append(row('chr1', 5 + i % 3, '<INS>', [field], vcf.SV_FORMAT,
                        'SVTYPE=INS;MAXSIZE=3;SEQ=ACG'))
    paths = [source(tmp_path / 'input.vcf', ['a'], rows)]
    _, expected = run(tmp_path, paths, 'python', monkeypatch, snp_only=False)
    _, actual = run(tmp_path, paths, 'native', monkeypatch, snp_only=False)
    assert actual == expected


def test_randomized_sparse_emission(tmp_path, monkeypatch):
    """Raw tuple ordering, labels, filters and call states survive buffer reuse."""
    rng = random.Random(20261009)
    samples = [f's{i}' for i in range(127)]
    all_samples = samples + ['absent-from-header']
    writer = compact.Writer(tmp_path / 'source')
    for _ in range(12000):
        sample = rng.choice(all_samples)
        pos = rng.randrange(1, 90)
        alt = rng.choice('CTN')
        filt = rng.choice(['PASS', '.', 'q20;q10', 'q10', 'PASS;q10', ''])
        state = rng.choice(['1', '1', '1', '0', '.'])
        values = None
        if state == '1':
            values = ['1', rng.choice(['SNP', 'INS_SNP']), '1', alt,
                      rng.choice(['query', 'query:one', 'query;two']),
                      str(rng.choice([0, 1, 2, 10, 100, 2**32-1])) + rng.choice('+-'),
                      rng.choice(['.', 'a', 'a,;:@\t\n', '雪']),
                      rng.choice(['.', '0H', '-2H', '10H', '.&-2H&0H', '2147483647H'])]
        writer.add('chr1', pos, 'A', alt, sample, values, state, filt)
    # Adjacent ALT rows at a position can have different states; marker-only
    # sites must not leak their fields into the next emitted site.
    writer.add('chr1', 100, 'A', 'C', 's1', state='.')
    writer.add('chr1', 101, 'A', 'T', 's0', ['1', 'SNP', '1', 'T', 'q', '10+', '.', '.'])
    for sample in samples[:5]:
        for qpos in range(120):
            writer.add('chr1', 200, 'A', 'C', sample,
                       ['1', 'SNP', '1', 'C', 'q', f'{qpos}+', 'long' * 20, '0H'])
    writer.add('chr1', 201, 'A', 'T', 's0', ['1', 'SNP', '1', 'T', 'q', '1+', '.', '.'])
    writer.coverage['chr1'] = [[s, start, end] for s in samples
                               for start, end in [(0, 10), (20, 30), (25, 45), (50, 80), (80, 120)]]
    data = writer.finish()
    root = tmp_path / 'shards'
    manifest = compact.manifest_from_sources([dict(source=str(writer.directory / 'source.json'),
                                                   samples=samples + ['s0'], chroms=data['chroms'], metadata=[])],
                                            root, str(tmp_path))
    key = compact.chrom_key('chr1')
    monkeypatch.setenv('GRAPHVCFMERGE_BACKEND', 'python')
    compact.merge_chrom(root, key)
    expected = small.part_path(manifest, key, 'snp').read_bytes()
    monkeypatch.setenv('GRAPHVCFMERGE_BACKEND', 'native')
    compact.merge_chrom(root, key)
    assert small.part_path(manifest, key, 'snp').read_bytes() == expected


def test_more_than_64_spill_runs(tmp_path, monkeypatch):
    """Exercise an intermediate merge pass before the final 64-way merge."""
    count = 800000
    rng = np.random.default_rng(41)
    writer = compact.Writer(tmp_path / 'source')
    names = [f's{i}' for i in range(128)]
    for sample in names:
        writer.add('chr1', 1, 'A', 'C', sample, ['1', 'SNP', '1', 'C', 'q', '1+', '.', '.'])
    data = writer.finish()
    key = compact.chrom_key('chr1')
    rows = np.empty(count, dtype=compact.DTYPE)
    rows['pos'] = rng.integers(1, 101, count, dtype=np.uint32)
    rows['bases'] = 1
    rows['query_pos'] = rows['pos']
    rows['query'] = rng.integers(0, len(names), count, dtype=np.uint32)
    rows['allele'] = 0
    rows['label_h'] = 0
    rows.tofile(writer.directory / 'records.bin')
    data['counts'][key] = count
    data['sections'][key] = [[0, count]]
    (writer.directory / 'source.json').write_text(json.dumps(data))
    root = tmp_path / 'shards'
    manifest = compact.manifest_from_sources([dict(source=str(writer.directory / 'source.json'),
                                                   samples=names, chroms=data['chroms'], metadata=[])],
                                            root, str(tmp_path))
    monkeypatch.setenv('GRAPHVCFMERGE_BACKEND', 'native')
    compact.merge_chrom(root, key, 2)
    expected = small.part_path(manifest, key, 'snp').read_bytes()
    monkeypatch.setenv('GRAPHVCFMERGE_RAM_GB', '0.001')
    assert count * 96 / native.memory_bytes(len(names)) > 64
    compact.merge_chrom(root, key, 2)
    assert small.part_path(manifest, key, 'snp').read_bytes() == expected
    assert not list((tmp_path / 'parts').rglob('*.run.*'))
