"""FASTA bulk reads and GFA extraction preserve sequences, validation, and offsets."""
import hashlib
import io
from pathlib import Path
import random
import subprocess
import sys

import numpy as np
import pytest

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / 'scripts'))
import gfa_interval_pipeline as pipeline
from gfa_query_anchors import reverse_complement
from gfa_topology import POINT_BITS, Segments
from minsetref_core import IndexedFasta


def fasta(path, records, width=60, newline=b'\n', final_newline=True):
    data, index = bytearray(), []
    for name, sequence in records:
        data.extend(b'>' + name.encode() + newline)
        offset = len(data)
        for low in range(0, len(sequence), width):
            data.extend(sequence[low:low+width].encode() + newline)
        bases = min(width, len(sequence)) or 1
        index.append(f'{name.split()[0]}\t{len(sequence)}\t{offset}\t{bases}\t{bases+len(newline)}\n')
    if not final_newline:
        del data[-len(newline):]
    path.write_bytes(data)
    Path(str(path)+'.fai').write_text(''.join(index))
    return str(path), str(path)+'.fai'


@pytest.mark.parametrize('width', [1, 60, 80, 5000])
@pytest.mark.parametrize('newline', [b'\n', b'\r\n'])
@pytest.mark.parametrize('indexed', [False, True])
def test_fasta_bulk_intervals(tmp_path, width, newline, indexed):
    rng = random.Random(31)
    records = [('empty', ''), ('a', 'ACgtNRY'*100), ('b', 'tGCA'*177)]
    path, fai = fasta(tmp_path/'input.fa', records, width, newline, False)
    if not indexed:
        Path(fai).unlink()
    with IndexedFasta(path) as reader:
        for name, sequence in records:
            assert reader.sequence(name) == sequence
            intervals = [(0, len(sequence)), (0, 0), (0, width),
                         (width-1, width+1), (width, 2*width)]
            intervals.extend((rng.randrange(-10, 900), rng.randrange(-10, 900))
                             for _ in range(100))
            for start, end in intervals:
                low = max(0, min(start, len(sequence)))
                high = max(low, min(end, len(sequence)))
                assert reader.fetch(name, start, end) == sequence[low:high]


def test_bulk_fetch_uses_one_read(tmp_path, monkeypatch):
    source = fasta(tmp_path/'wrapped.fa', [('a', 'ACGT'*30000)])
    with IndexedFasta(*source) as reader:
        calls = []
        original = reader._pread
        def read(size, offset):
            calls.append((size, offset))
            return original(size, offset)
        monkeypatch.setattr(reader, '_pread', read)
        assert reader.fetch('a', 17, 100017) == ('ACGT'*30000)[17:100017]
        assert len(calls) == 1


def request(name, sequence, low, high, start, end, strand, checks=None):
    core = sequence[start:end]
    if strand == '-':
        core = reverse_complement(core)
    if checks is None:
        checks = [(0, len(core), hashlib.sha256(core.upper().encode()).digest())]
    return name, low, high, start, end, strand, checks


def test_extraction_short_streamed_reverse_and_failed_candidates(tmp_path):
    sequence = 'aCGtnACGTTGCA'*20000
    source = fasta(tmp_path/'query.fa', [('a', sequence)])
    requests = [
        request('first', sequence, 0, 55, 0, 5, '+'),
        request('bad', sequence, 10, 111, 60, 61, '+', [(0, 1, b'bad')]),
        request('reverse', sequence, 15, 125, 65, 75, '-'),
        request('long', sequence, 200, 100300, 250, 100250, '-'),
        request('bad-long', sequence, 210, 100310, 260, 100260, '+', [(0, 1, b'bad')]),
        request('later', sequence, 100400, 100501, 100450, 100451, '+'),
        request('unchecked', sequence, 100510, 100620, 100560, 100570, '-', []),
        request('tail', sequence, len(sequence)-55, len(sequence),
                len(sequence)-5, len(sequence), '-'),
    ]
    output = str(tmp_path/'sequences')
    hits, failures = pipeline._extract_batch((source, 'a', requests, output))
    assert [name for name, _ in failures] == ['bad', 'bad-long']
    expected, expected_hits = bytearray(), []
    for name, low, high, start, end, strand, _checks in requests:
        if name.startswith('bad'):
            continue
        data = sequence[low:high]
        if strand == '-':
            data = reverse_complement(data)
        left = start-low if strand == '+' else high-end
        expected_hits.append((name, output, len(expected), len(data), left,
                              len(data)-left-(end-start)))
        expected.extend(data.encode())
    assert hits == expected_hits
    assert Path(output).read_bytes() == expected


def test_nearby_snps_share_reads(tmp_path, monkeypatch):
    sequence = 'ACGT'*10000
    source = fasta(tmp_path/'query.fa', [('a', sequence)])
    calls, original = [], IndexedFasta.fetch
    def fetch(self, name, start, end):
        calls.append((start, end))
        return original(self, name, start, end)
    monkeypatch.setattr(IndexedFasta, 'fetch', fetch)
    requests = [request(str(i), sequence, i, i+101, i+50, i+51, '+')
                for i in range(2048)]
    hits, failures = pipeline._extract_batch((source, 'a', requests, str(tmp_path/'sequences')))
    assert len(hits) == 2048 and not failures
    assert len(calls) == 1


def test_preloaded_roots_share_duplication_record_and_match_streaming(tmp_path, monkeypatch):
    # A segment larger than CHUNK also exercises the bounded output slicing.
    sequence = 'ACGT' * (pipeline.CHUNK // 2 + 7)
    path, fai = fasta(tmp_path/'roots.fa', [('a', sequence), ('alt', 'acgtn'), ('unused', 'TTT')])
    roots = {
        'a': pipeline.Root('a', path, fai, len(sequence), True, 'reference', 0),
        'copy': pipeline.Root('copy', path, fai, 10, True, 'duplication', 1, record='a', base=3),
        'alt': pipeline.Root('alt', path, fai, 5, True, 'alternative', 2),
        # An unused alias need not exist as a record in its backing FASTA.
        'alias': pipeline.Root('alias', path, fai, 5, False, 'reference', -1),
    }
    keys = np.array([0, len(sequence), 1 << POINT_BITS, (1 << POINT_BITS)+10,
                     2 << POINT_BITS, (2 << POINT_BITS)+5], np.int64)
    used = np.array([True, False]*3)
    segments = Segments(keys, used, np.cumsum(used), list(roots), np.zeros(4, bool))
    sequences = pipeline._preload_root_sequences(segments, roots, lambda _message: None)
    assert sequences == {(path, 'a'): sequence, (path, 'alt'): 'acgtn'}
    ranks = dict.fromkeys(roots, 0)
    streamed, cached = io.BytesIO(), io.BytesIO()
    pipeline._write_segments(streamed, segments, {}, roots, {}, '', ranks)
    def unexpected_read(*args):
        raise AssertionError('preloaded writer reopened a FASTA')
    monkeypatch.setattr(pipeline, 'IndexedFasta', unexpected_read)
    pipeline._write_segments(cached, segments, {}, roots, {}, '', ranks, sequences)
    assert streamed.getvalue() == cached.getvalue()


def test_converter_outputs_match_with_preload_workers_and_partitions(tmp_path):
    from graphvcfmerge import SNP_FORMAT, SV_FORMAT_LEGACY

    reference, _ = fasta(tmp_path/'ref.fa', [('chr1', 'A'*100000)])
    query, _ = fasta(tmp_path/'query.fa', [('ctg', 'C'*100000)])
    queries = tmp_path/'queries.list'
    queries.write_text(f'Q {query}\n')
    template, _ = fasta(tmp_path/'template.fa', [
        ('alt source=Q:ctg:500-6500+ backbone=REF lift_status=unmapped', 'C'*6000)])
    header = ('##fileformat=VCFv4.2\n##contig=<ID=chr1,length=100000>\n'
              '#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tQ\n')
    snps = tmp_path/'snp.vcf'
    # Enough records to use multiple extraction and topology batches. Reverse
    # coordinates deliberately need the second candidate to match ALT.
    snps.write_text(header + ''.join(
        f'chr1\t{i+1}\tS_{i}\tA\tC\t.\tPASS\t.\t{SNP_FORMAT}\t'
        f'1:SNP:1:C:ctg:{i+1}-:allele:0H\n' for i in range(1000, 3100)))
    svs = tmp_path/'sv.vcf'
    svs.write_text(header +
        f'alt\t100\tI_1\tN\t<INS>\t.\tPASS\tSVTYPE=INS;SVLEN=30;END=100;SEQ={"C"*30}\t'
        f'{SV_FORMAT_LEGACY}\t1:INS:30:>30=:ctg:600-630+:allele:0H\n' +
        'chr1\t500\tI_2\tN\t<INS>\t.\tPASS\tSVTYPE=INS;SVLEN=30;END=500;SEQ=>I_1:30=\t'
        f'{SV_FORMAT_LEGACY}\t1:INS:30:>30=:ctg:600-630+:allele:0H\n' +
        'I_2\t5\tI_3\tN\t<INS>\t.\tPASS\tSVTYPE=INS;SVLEN=10;END=5;SEQ=>I_1:4H10=\t'
        f'{SV_FORMAT_LEGACY}\t1:INS:10:>10=:ctg:604-614+:allele:0H\n')
    script = Path(pipeline.__file__).with_name('merged_vcf_to_gfa.py')
    common = [sys.executable, str(script), '-v', str(svs), str(snps), '-q', str(queries),
              '-r', reference, '--local-reference-templates', template,
              '--reference-haplotype', 'REF', '--all']
    expected = None
    for number, options in enumerate([
            ['--processes', '1', '--no-preload-reference'],
            ['--processes', '2'],
            ['--processes', '2', '--partition-records', '1000'],
            ['--processes', '2', '--partition-records', '1000', '--no-preload-reference']]):
        folder = tmp_path/str(number)
        folder.mkdir()
        result = subprocess.run(common + options + ['-o', str(folder/'cohort.gfa')],
                                capture_output=True, text=True)
        assert result.returncode == 0, result.stderr
        actual = {path.name: path.read_bytes() for path in folder.glob('cohort*')}
        assert len(actual) == 7
        if expected is None:
            expected = actual
        else:
            assert actual == expected
