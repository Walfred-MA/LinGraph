"""Regressions independently reproduced from the September 21 alignment audit."""
import os
import random
import unittest
from unittest import mock

import pytest
import graphcigartoref as c
from tests.test_late_di_rescue_score import Reader as BaseReader, fixture
from tests.test_tandem_real_trio import call_owned, verify_spelling


class Reader(BaseReader):
    def __contains__(self, name):
        return name in self.index


def dna(rng, n):
    return ''.join(rng.choices('ACGT', k=n))


def cigar(ops):
    return ''.join(f'{n}{"=" if op == "M" else op}' for n, op, _ in ops)


def check_spelling(reference, query, ops):
    text = c.query_only_graph_cigar('>chr:' + c._tm_format_pairwise_ops(ops))
    verify_spelling(unittest.TestCase(), text, query, Reader(reference))


def calls(reference, query, text):
    sample = dict(sequence=query, owned_query_slice=(0, len(query)),
                  strand='+', start=0, end=len(query), contig='query_chr')
    return call_owned(dict(chrom='chr', chrom_length=len(reference)), sample, text)


def repeat_pair():
    rng = random.Random(17)
    left, right, motif = dna(rng, 300), dna(rng, 300), dna(rng, 37)
    return left, right, motif, left + motif * 4 + right, left + motif * 7 + right


@pytest.mark.parametrize('body', ['300=37I37=74I111=300=', '300=111=37I37=74I300='])
@pytest.mark.parametrize('second_pass', [False, True])
def test_moving_repeat_gap_can_replace_complete_window(body, second_pass):
    _, _, _, ref, qry = repeat_pair()
    old = c._tm_ops_from_cigar_on_sequences(ref, qry, body)
    if second_pass:
        old = [(n, 'M' if op == '=' else op, p) for n, op, p in old]
        new = c._tm_realign_score_linked_indel_windows(old, ref, qry)
    else:
        new = c._tm_consolidate_sv_ops(old, ref, qry)
    check_spelling(ref, qry, new)
    assert [n for n, op, _ in new if op == 'I'] == [111]
    assert c.alignment_score(new) > c.alignment_score(old)
    events = calls(ref, qry, '>chr:' + c._tm_format_pairwise_ops(new))
    assert events == [('INS', 300, 300, 111, 0)]


def test_upstream_base_alignment_uses_same_orientation_as_downstream():
    left, right, motif, ref, qry = repeat_pair()
    outputs = []
    for upstream in (False, True):
        text = c._align_graph_cigar_extension(
            '>p:300M148I' + motif * 4 + '300M',
            '>p:300M259I' + motif * 7 + '300M',
            'ref', 'query', {'p': left + right}, outward_from_right_edge=upstream)
        text = c.query_only_graph_cigar('>chr:' + text)
        verify_spelling(unittest.TestCase(), text, qry, Reader(ref))
        outputs.append(calls(ref, qry, text))
    assert outputs[0] == outputs[1]
    assert outputs[0][0][3] == 111


@pytest.mark.parametrize('kind', ['INS', 'DEL'])
@pytest.mark.parametrize('strand', ['+', '-'])
def test_core_boundary_preserves_64bp_sv_through_calling(kind, strand):
    rng = random.Random(413)
    if kind == 'INS':
        ref, payload = dna(rng, 720), dna(rng, 64)
        qry = ref[:600] + payload + ref[600:]
        graph_query = f'>node:600M64I{payload}120M'
        rcore, qcore = (100, 620), (100, 684)
    else:
        ref = dna(rng, 784)
        qry = ref[:120] + ref[184:]
        graph_query = '>node:120M64D600M'
        rcore, qcore = (100, 684), (100, 620)
    backbone = c.Coord('chr', 0, len(ref), strand)
    qcoord = c.Coord('query_chr', 0, len(qry), strand)
    pair = c.PairRow(1, '', 'query', 'reference',
                    f'query_chr:0-{len(qry)}{strand}', f'chr:0-{len(ref)}{strand}',
                    query_name='query', ref_name='reference', query_coord=qcoord, ref_coord=backbone)
    pair.query_core_slice, pair.ref_core_slice = qcore, rcore
    reader = Reader(ref if strand == '+' else c.revcomp(ref))
    row = c.build_provisional_row(pair,
        {'reference': f'>node:{len(ref)}M', 'query': graph_query},
        {'node': ref}, {'reference': ref, 'query': qry}, {},
        {'chr': len(ref)}, {}, reader, {'query': qcoord, 'reference': backbone})
    text = row.split('\t')[6]
    verify_spelling(unittest.TestCase(), text, qry, reader)
    events = calls(reader.sequence, qry, text)
    assert len(events) == 1
    assert events[0][3:] == ((64, 0) if kind == 'INS' else (0, 64))


def divergent_repeat(seed):
    rng = random.Random(seed)
    motif = dna(rng, 41)
    def mutate(s, rate):
        return ''.join(rng.choice([x for x in 'ACGT' if x != b])
                       if rng.random() < rate else b for b in s)
    left, right = dna(rng, 100), dna(rng, 100)
    repeat = ''.join(mutate(motif, .035) for _ in range(8))
    return (left + repeat + right,
            left + mutate(repeat[:164], .025) + ''.join(mutate(motif, .035) for _ in range(4))
            + mutate(repeat[164:], .025) + right)


@pytest.mark.parametrize('seed', range(6))
@pytest.mark.parametrize('reverse', [False, True])
def test_short_di_affine_repair_does_not_erase_sv_into_microgaps(seed, reverse):
    ref, qry = divergent_repeat(seed)
    gap = 'I'
    if reverse:
        ref, qry, gap = qry, ref, 'D'
    new = c._realign_long_di_sequences(ref, qry)
    assert new is not None
    check_spelling(ref, qry, new)
    sizes = [n for n, op, _ in new if op == gap]
    assert sum(sizes) == 164 and max(sizes) >= 50
    raw = c.parse_cigar_ops(f'{len(ref)}D{len(qry)}I' + qry)
    assert c._di_rescue_is_valid(raw, new, ref, qry)
    # The specified affine model can still prefer multiple repeat-copy gaps;
    # the regression is losing the entire SV below the calling threshold.


@pytest.mark.parametrize('strand', ['+', '-'])
@pytest.mark.parametrize('kind', ['INS', 'DEL'])
def test_asymmetric_di_pairs_have_sequence_rescue(kind, strand):
    rng = random.Random(826)
    shared, extra = dna(rng, 20), dna(rng, 64)
    ref, qry = (shared + extra, shared) if kind == 'DEL' else (shared, shared + extra)
    args = fixture('40=' + f'{len(ref)}D{len(qry)}I' + qry + '40=',
                   'A' * 40 + ref + 'C' * 40, 'A' * 40 + qry + 'C' * 40, strand)
    piece, linear = c._rescue_long_di_pairs_in_linear_piece(*args)
    assert piece is not args[0]
    assert c.pairwise_ref_span(piece.cigar) == len(ref) + 80
    assert c.pairwise_query_span(piece.cigar) == len(qry) + 80
    assert [(tok.n, tok.op) for tok in c.parse_cigar_ops(piece.cigar) if tok.op in 'ID'] == [(64, 'D' if kind == 'DEL' else 'I')]


def test_late_deletion_consolidation_is_symmetric():
    _, _, _, short, long = repeat_pair()
    text = '>chr:300=37D37=74D111=300='
    result = c._rescue_serialized_di_pairs(text, short, Reader(long))
    verify_spelling(unittest.TestCase(), result, short, Reader(long))
    events = calls(long, short, result)
    assert len(events) == 1 and events[0][3:] == (0, 111)


def test_late_consolidated_insertion_passes_through_encoder():
    _, _, _, ref, qry = repeat_pair()
    stage1 = c.build_stage1('>p:748M', '>p:300M111I' + qry[300:411] + '448M', 'r', 'q')
    text = c.query_only_graph_cigar('>chr:' + c._tm_format_pairwise_ops(
        c._tm_ops_from_cigar_on_sequences(ref, qry, '300=37I37=74I111=300=')))
    with mock.patch.object(c, 'encode_large_insertions_in_piece',
                           wraps=c.encode_large_insertions_in_piece) as encode:
        result = c._rescue_serialized_di_pairs(text, qry, Reader(ref), stage1=stage1,
                                               graph_sequences={'p': ref})
    assert encode.call_count == 1
    assert [t.n for t in c.parse_cigar_ops(encode.call_args.args[0].cigar) if t.op == 'I'] == [111]
    verify_spelling(unittest.TestCase(), c.query_only_graph_cigar(result), qry, Reader(ref))


def test_stage4_consolidates_before_encoding():
    _, _, _, ref, qry = repeat_pair()
    stage1 = c.build_stage1('>p:748M', '>p:300M111I' + qry[300:411] + '448M', 'r', 'q')
    body = c._tm_format_pairwise_ops(c._tm_ops_from_cigar_on_sequences(ref, qry, '300=37I37=74I111=300='))
    piece = c.PairwisePiece('main', 'p', 0, 0, 0, len(ref), 0, len(ref), body)
    with mock.patch.object(c, 'encode_large_insertions_in_piece', wraps=c.encode_large_insertions_in_piece) as encode:
        stage4 = c.build_stage4_linear(stage1, c.Stage3PolishResult([piece], body), {'p': ref}, {},
            {'chr': len(ref)}, {}, Reader(ref), ref_backbone_coord=c.Coord('chr', 0, len(ref), '+'),
            record_sequences={'q': qry}, query_name='q')
    assert encode.call_count == 1
    assert [t.n for t in c.parse_cigar_ops(encode.call_args.args[0].cigar) if t.op == 'I'] == [111]


@pytest.mark.parametrize('backend', ['0', '1'])
def test_native_forward_mapping_cannot_be_suppressed_by_reverse_primary(backend):
    if (backend == '0' and c._mappy is None) or (backend == '1' and not c._minimap2_executable()):
        pytest.skip('native mapper unavailable')
    rng = random.Random(992)
    qry = dna(rng, 3000)
    mutated = ''.join(rng.choice([x for x in 'ACGT' if x != b]) if rng.random() < .02 else b for b in qry)
    ref = c.revcomp(qry) + dna(rng, 1000) + mutated
    with mock.patch.dict(os.environ, {'GRAPH_CIGARTOREF_NO_MAPPY': backend}):
        new = c.minimap2_payload_ops(ref, qry)
    assert sum(n for n, op, _ in new if op == '=') == 2928
    check_spelling(ref, qry, new)


@pytest.mark.parametrize('backend', ['0', '1'])
def test_payload_assembly_retains_both_separated_native_hits(backend):
    if (backend == '0' and c._mappy is None) or (backend == '1' and not c._minimap2_executable()):
        pytest.skip('native mapper unavailable')
    rng = random.Random(6421)
    a, b = dna(rng, 1500), dna(rng, 1500)
    ref, qry = a + dna(rng, 20000) + b, a + dna(rng, 20000) + b
    with mock.patch.dict(os.environ, {'GRAPH_CIGARTOREF_NO_MAPPY': backend}):
        new = c.minimap2_payload_ops(ref, qry)
    assert sum(n for n, op, _ in new if op == '=') == 3000
    check_spelling(ref, qry, new)


def test_early_payload_rejects_worse_global_score():
    rng = random.Random(390)
    shared = dna(rng, 400)
    ref, qry = dna(rng, 5000) + shared + dna(rng, 5000), dna(rng, 5000) + shared + dna(rng, 5000)
    # Candidate span and matches are valid, but its four large gaps cost more
    # than the original two under the unchanged SV logarithmic score.
    candidate = c._tm_ops_from_cigar_on_sequences(ref, qry, '5000D5000I400=5000D5000I')
    with mock.patch.object(c, 'minimap2_payload_ops', return_value=candidate):
        result = []
        c._tm_realign_insertion_payloads(ref, qry, result)
    assert cigar(result) == '10400D10400I'
    check_spelling(ref, qry, result)


def test_trim_accounting_and_reverse_insertion_anchor():
    piece = c.PairwisePiece('main', 'p', 0, 0, 0, 115, 0, 115, '1X3=2X9=5I100=')
    polished = c.polish_match_piece(piece)
    assert c.pairwise_query_span(polished.cigar) == 120
    for sign in ('>', '<'):
        text = c.slice_graphic_mapping_by_query(sign + 'chr:100H10=5IAAAAA10=100H', 10, 15, None)
        segment = c.parse_graphic_segments(text, 'slice')[0]
        assert segment.start == segment.end == 110
    assert c._reverse_template_cigar('2DAC') == '2DGT'


def test_native_and_fallback_affine_objectives_agree():
    if c._parasail is None:
        pytest.skip('native affine comparison requires parasail')
    rng = random.Random(552)
    for match, mismatch, gap_open, gap_extend in [(4, -8, 12, 1), (1, -4, 4, 1)]:
        def score(ops):
            return sum(match * n if op in '=M' else mismatch * n if op == 'X'
                       else -gap_open - (n - 1) * gap_extend for n, op, _ in ops)
        for _ in range(30):
            ref, qry = dna(rng, rng.randrange(5, 50)), dna(rng, rng.randrange(5, 50))
            kwargs = dict(match_score=match, mismatch_score=mismatch, gap_open=gap_open, gap_extend=gap_extend)
            native = c._global_affine_payload_align(ref, qry, **kwargs)
            with mock.patch.object(c, '_parasail', None):
                fallback = c._global_affine_payload_align(ref, qry, **kwargs)
            check_spelling(ref, qry, native)
            check_spelling(ref, qry, fallback)
            assert score(native) == score(fallback)


@pytest.mark.parametrize('is_main', [True, False])
@pytest.mark.parametrize('strand', ['+', '-'])
def test_late_encoding_preserves_main_and_nested_reference_coordinates(is_main, strand):
    _, _, _, ref, qry = repeat_pair()
    payload = qry[300:411]
    chromosome = 'A' * 50 + 'T' * 40 + 'N' * 110 + ref + 'N' * 100 + payload + 'N' * 50
    if strand == '-':
        chromosome = c.revcomp(chromosome)
    reader = Reader(chromosome)
    total = len(chromosome)
    direction = '>' if strand == '+' else '<'
    # Mapping clips are measured in traversal orientation, on either strand.
    mappings = {'p': f'{direction}chr:200H748={total-948}H',
                'dup': f'{direction}chr:1048H111=50H'}
    graph_query = '>p:300M448H>dup:111M>p:300H448M'
    if not is_main:
        graph_query = '>before:50M' + graph_query + '>after:40M'
    stage1 = c.build_stage1('>p:748M', graph_query, 'r', 'q')
    old = c._tm_ops_from_cigar_on_sequences(ref, qry, '300=37I37=74I111=300=')
    body = c.query_only_graph_cigar('>chr:' + c._tm_format_pairwise_ops(old)).split(':', 1)[1]
    segment_text = f'{direction}chr:200H{body}{total-948}H'
    if is_main:
        text, query = segment_text, qry
    else:
        text = f'{direction}chr:50=' + segment_text + f'{direction}40={total-90}H'
        query = 'A' * 50 + qry + 'T' * 40
    result = c._rescue_serialized_di_pairs(text, query, reader, stage1=stage1,
                                           graph_sequences={'p': ref, 'dup': payload,
                                                            'before': 'A' * 50, 'after': 'T' * 40},
                                           graph_mappings=mappings)
    result = c.query_only_graph_cigar(result)
    verify_spelling(unittest.TestCase(), result, query, reader)
    segments = c.parse_reference_cigar_segments(result, 'after')
    mains = [s for s in segments if s.is_main]
    assert sum(s.end - s.start for s in mains) == (748 if is_main else 90)
    encoded = [s for s in segments if not s.is_main]
    if is_main:
        # On '-' normalization moves right in traversal orientation, so the
        # same repeat is encoded from p instead of the original dup chunk.
        expected = (1048, 1159) if strand == '+' else (561, 672)
        assert [(s.start, s.end) for s in encoded] == [expected]
    else:
        # Already encoded segments retain their own nested indel. They must
        # not be split into separate sibling encodings by the late pass.
        assert len(encoded) == 1
        assert [t.n for t in encoded[0].ops if t.op == 'I'] == [111]
        called = calls(chromosome, query, result)
        assert any(event[0] == 'INS_INS' and event[3:] == (111, 0) for event in called)


def test_no_parasail_repeat_consolidation_preserves_one_sv():
    _, _, _, ref, qry = repeat_pair()
    old = c._tm_ops_from_cigar_on_sequences(ref, qry, '300=111=37I37=74I300=')
    with mock.patch.object(c, '_parasail', None):
        new = c._tm_consolidate_sv_ops(old, ref, qry)
    check_spelling(ref, qry, new)
    assert [n for n, op, _ in new if op == 'I'] == [111]
