"""Hardness eligibility, stage-specific scores, and unchanged adjacent D/I."""
import math
import random
from unittest import mock

import pytest
import graphcigartoref as core
import graphcigartoref_persample as per
from alignment_scoring import MaskedSVScorer, whole_pair_score
from reference_hardness import quantify_reference_hardness, quantify_sequence_hardness
from tests.test_late_di_rescue_score import fixture
from tests.test_tandem_minimap_realign import repeated_stage, stage
from types import SimpleNamespace


@pytest.mark.parametrize('sample,expected', [
    ('father', 7.7714608901416895), ('child', 7.724704687979818),
])
def test_saved_full_query_uses_historical_70_30_hardness(sample, expected):
    import gzip
    import json
    from pathlib import Path
    path = Path(__file__).with_name('fixtures') / 'na19239_na19240_group81264.json.gz'
    query = json.loads(gzip.decompress(path.read_bytes()))['samples'][sample]['sequence']
    report = quantify_sequence_hardness(query)
    summary = report['summary']
    assert report['parameters']['local_weight'] == 0.7
    assert report['score_model'] == 'local-complexity-v2'
    assert summary['hardness_score'] == pytest.approx(expected, abs=1e-10)
    assert summary['hardness_score'] == pytest.approx(
        0.3 * summary['seed_hardness_score'] + 70 * summary['local_complexity_risk'])
    assert core._realignment_hardness_allows(query)


@pytest.mark.parametrize('sequence', [
    'acgt' * 1000, 'ACGT' * 10 + 'NNNN' + 'a' * 200,
    ''.join(random.Random(92).choices('ACGTacgt', k=2000)), 'N' * 30, 'ACG',
])
def test_in_memory_hardness_is_identical_to_file_api(tmp_path, sequence):
    path = tmp_path / 'query.fa'
    path.write_text('>query\n' + sequence + '\n')
    memory = quantify_sequence_hardness(sequence)
    disk = quantify_reference_hardness(path)
    assert memory.pop('reference') is None
    disk.pop('reference')
    assert memory == disk


@pytest.mark.parametrize('score,allowed', [(5, True), (9.999, True), (10, False), (10.001, False), (None, False)])
def test_whole_pair_boundary_is_checked_before_mapping(score, allowed):
    sequence = 'ACGT' * 50
    reader = SimpleNamespace(index={'chr': (200,)}, fetch=lambda *_: sequence)
    original = stage('100D100I' + sequence[:100] + '100=')
    hit = core._tandem_hit_from_cigar(sequence, sequence, 0, 200, 0, 200, '200M')
    with mock.patch.object(core, 'quantify_sequence_hardness',
                           return_value={'summary': {'hardness_score': score}}) as hardness, \
         mock.patch.object(core, '_minimap2_tandem_hits', return_value=[hit]) as mapper:
        result = core._realign_tandem_stage2(repeated_stage(), original,
                                            core.Coord('chr', 0, 200, '+'), sequence, reader)
    hardness.assert_called_once_with(sequence)
    assert mapper.call_count == int(allowed)
    assert (result is not original) == allowed


@pytest.mark.parametrize('strand', ['+', '-'])
@pytest.mark.parametrize('flank_base', ['A', 'N'])
def test_early_hardness_uses_whole_comparison_query(strand, flank_base):
    reference = ''.join(random.Random(68324).choices('ACGT', k=200))
    query = reference[:150] + reference[50:]
    flanks = flank_base * 4000
    contig = flanks + (query if strand == '+' else core.revcomp(query)) + flanks
    qcoord = core.Coord('assembly', len(flanks), len(flanks) + len(query), strand)
    pair = core.PairRow(1, '', 'q', 'r', 'q', 'r', query_name='q', query_coord=qcoord)
    graph_cigar = '>node:150M50H>node:50H150M'
    records = per._query_record_sequences(pair, {'q': graph_cigar}, {'assembly': contig}, {})
    assert records['q'] == query
    expected = quantify_sequence_hardness(query)['summary']['hardness_score']
    assert expected < 10
    if flank_base == 'A':
        assert quantify_sequence_hardness(contig)['summary']['hardness_score'] > 10
    original = stage('150=100I' + reference[50:150] + '50=')
    hit = core._tandem_hit_from_cigar(reference, query, 0, 300, 0, 200, '150M100I50M')
    assert hit is not None
    reader = SimpleNamespace(index={'chr': (200,)}, fetch=lambda *_: reference)
    with mock.patch.object(core, 'quantify_sequence_hardness',
                           wraps=quantify_sequence_hardness) as hardness, \
         mock.patch.object(core, '_minimap2_tandem_hits', return_value=[hit]) as mapper:
        core._realign_tandem_stage2(repeated_stage(), original,
                                    core.Coord('chr', 0, 200, '+'), records['q'], reader)
    hardness.assert_called_once_with(query)
    mapper.assert_called_once_with(reference, query)


@pytest.mark.parametrize('helper,backend', [
    (core._tm_consolidate_sv_ops, '_sv_window_affine_align'),
    (core._tm_realign_score_linked_indel_windows, '_tm_realign_secondary_sv_window'),
])
def test_high_hardness_sv_windows_still_realign_and_improve(helper, backend):
    reference = 'ac' * 2000
    query = reference[:1000] + 'ac' * 40 + reference[1000:]
    assert not core._realignment_hardness_allows(query)
    original = [(n, 'M' if op == '=' else op, payload) for n, op, payload in
                core._tm_ops_from_cigar_on_sequences(reference, query, '1000=40I20=40I2980=')]
    with mock.patch.object(core, '_SV_REALIGNMENT_ENABLED', True), \
         mock.patch.object(core, '_realignment_hardness_allows',
                           side_effect=AssertionError('SV window was gated')), \
         mock.patch.object(core, backend, wraps=getattr(core, backend)) as align:
        result = helper(original, reference, query)
    assert align.call_count > 0
    assert [(n, op) for n, op, _ in result if op in 'ID'] == [(80, 'I')]
    assert MaskedSVScorer(reference, query)(result) > MaskedSVScorer(reference, query)(original)
    assert sum(n for n, op, _ in result if op in 'M=XD') == len(reference)
    assert sum(n for n, op, _ in result if op in 'M=XI') == len(query)


@pytest.mark.parametrize('strand', ['+', '-'])
@pytest.mark.parametrize('improved', [False, True])
def test_composite_sv_uses_score_without_hardness(strand, improved):
    body = '40=100D20=100I' + 'A' * 100 + '40='
    args = fixture(body, 'A' * 200, 'A' * 200, strand)
    assert not core._realignment_hardness_allows('A' * 120)
    candidate = [(120, 'M', '')] if improved else [(100, 'D', ''), (20, 'M', ''), (100, 'I', 'A' * 100)]
    with mock.patch.object(core, '_SV_REALIGNMENT_ENABLED', True), \
         mock.patch.object(core, '_realignment_hardness_allows',
                           side_effect=AssertionError('Composite SV was gated')), \
         mock.patch.object(core, '_realign_long_di_sequences') as di, \
         mock.patch.object(core, '_tm_realign_secondary_sv_window', return_value=candidate) as sv:
        piece, linear = core._rescue_internal_di_windows_in_linear_piece(*args)
    assert piece.cigar == ('200=' if improved else body)
    assert (linear.genome_start, linear.genome_end, linear.strand) == (100, 300, strand)
    di.assert_not_called()
    sv.assert_called_once()


@pytest.mark.parametrize('strand', ['+', '-'])
@pytest.mark.parametrize('backend', ['native', 'M', '='])
def test_successive_composite_repairs_advance_both_sequence_cursors(strand, backend):
    """SV backends emit M; subsequent windows must still use their own bases."""
    reference = ''.join(random.Random(79530).choices('ACGT', k=360))
    cigar = ('40=100D20=100I' + reference[60:160]
             + '40=100D20=100I' + reference[220:320] + '40=')
    args = fixture(cigar, reference, reference, strand)
    native = core._tm_realign_secondary_sv_window
    seen = []

    def align(ref_window, query_window, **kwargs):
        seen.append((ref_window, query_window))
        if backend == 'native':
            return native(ref_window, query_window, **kwargs)
        return [(120, backend, '')]

    with mock.patch.object(core, '_SV_REALIGNMENT_ENABLED', True), \
         mock.patch.object(core, '_realignment_hardness_allows', return_value=True), \
         mock.patch.object(core, '_tm_realign_secondary_sv_window', side_effect=align):
        piece, linear = core._rescue_internal_di_windows_in_linear_piece(*args)
    assert seen == [(reference[40:160], reference[40:160]),
                    (reference[200:320], reference[200:320])]
    assert piece.cigar == '360='
    assert linear.cigar == piece.cigar
    assert (linear.genome_start, linear.genome_end, linear.strand) == (100, 460, strand)


@pytest.mark.parametrize('strand', ['+', '-'])
def test_composite_repair_stores_sequence_verified_matches_and_payloads(strand):
    reference = ''.join(random.Random(79531).choices('ACGT', k=200))
    base = next(base for base in 'ACGT' if base != reference[99]).lower()
    query = reference[:99] + base + reference[100:]
    cigar = '40=100D20=100I' + query[60:160] + '40='
    args = fixture(cigar, reference, query, strand)
    with mock.patch.object(core, '_SV_REALIGNMENT_ENABLED', True), \
         mock.patch.object(core, '_realignment_hardness_allows', return_value=True), \
         mock.patch.object(core, '_tm_realign_secondary_sv_window', return_value=[(120, 'M', '')]):
        piece, linear = core._rescue_internal_di_windows_in_linear_piece(*args)
    # Internal X payloads hold reference followed by query; serialization
    # later emits only the query bases.
    assert piece.cigar == f'99=1X{reference[99]}{base}100='
    assert linear.cigar == piece.cigar


@pytest.mark.parametrize('strand', ['+', '-'])
@pytest.mark.parametrize('helper', [core._rescue_long_di_pairs_in_linear_piece,
                                  core._rescue_internal_di_windows_in_linear_piece])
def test_adjacent_di_repair_never_consults_hardness(helper, strand):
    body = '40=100D100I' + 'a' * 100 + '40='
    args = fixture(body, 'a' * 180, 'a' * 180, strand)
    with mock.patch.object(core, '_realignment_hardness_allows', side_effect=AssertionError('D/I was gated')):
        piece, linear = helper(*args)
    assert piece.cigar == '180='
    assert (linear.genome_start, linear.genome_end, linear.strand) == (100, 280, strand)


def test_masked_gap_formula_axes_offsets_and_payload_case():
    scorer = MaskedSVScorer('AAaaCC', 'TTTTtGG')
    ops = [(2, '=', ''), (2, 'D', 'AA'), (3, 'I', 'TTT'), (2, 'X', '')]
    # Deleted aa: 2 masked; inserted TTt: 1 masked. Payloads are deliberately
    # uppercase to ensure they cannot replace the authoritative FASTA case.
    expected = 8 - 16 - (24 + 16 + 160 * math.log2(3)) - (24 + 24 + 160)
    assert scorer(ops) == pytest.approx(expected)
    assert scorer([(1, 'D', ''), (1, 'D', ''), (3, 'I', '')], 2, 2) == pytest.approx(expected + 8)
    assert scorer([(2, 'D', '')], 0, 0) == -40  # zero masked bases
    assert whole_pair_score(ops) == 8 - 16 - (24 + 16) - (24 + 24)


def test_masked_sv_boundary_joining_uses_exact_full_sequence_offsets():
    reference, query = 'AAA', 'aAAAC'
    original = core.parse_cigar_ops('1IA3=')
    candidate = [(3, '=', ''), (1, 'I', 'C')]
    scorer = MaskedSVScorer(reference, query)
    before = core._tm_ops_from_cigar_on_sequences(reference, query, '2I3=')
    after = core._tm_ops_from_cigar_on_sequences(reference, query, '1I3=1I')
    assert not core._sv_replacement_is_valid(
        original, candidate, reference, query[1:],
        left_context=[core.CigarOp(1, 'I', 'a')],
        scorer=scorer, reference_offset=0, query_offset=1)
    assert scorer(after) < scorer(before)


def test_sv_score_rechecks_matches_and_rejects_span_changes():
    original = core.parse_cigar_ops('1IA3=')
    assert not core._sv_replacement_is_valid(original, [(3, '=', ''), (1, 'I', 'A')], 'AAA', 'CAAA')
    assert not core._sv_replacement_is_valid(original, [(4, '=', '')], 'AAA', 'CAAA')
    with pytest.raises(ValueError):
        MaskedSVScorer('AAA', 'AAA')([(4, 'I', '')])


def test_sv_mask_penalty_does_not_use_neighboring_bases():
    ops = [(3, 'I', '')]
    for flank in ('A' * 600, 'a' * 600):
        scorer = MaskedSVScorer('', flank + 'AcG' + flank)
        assert scorer(ops, 0, 600) == -(24 + 24 + 160)


def test_local_minimap_hit_score_counts_masks_at_hit_offsets():
    reference = 'AAAA' + 'CCaaGG' + 'AAAA'
    query = 'TTTTTT' + 'CCtGG' + 'TTTT'
    scorer = MaskedSVScorer(reference, query)
    hit = core._tandem_hit_from_cigar(
        reference, query, 6, 11, 4, 10, '2M2D1I2M', score=scorer)
    assert hit is not None
    expected = 16 - (24 + 16 + 160 * math.log2(3)) - (24 + 8 + 160)
    assert hit.score == pytest.approx(expected)
    assert hit.score != scorer(hit.ops)
