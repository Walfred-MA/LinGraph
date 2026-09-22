"""Canonical positions and nested VCF calls must not depend on repair route."""
import random
import unittest
from unittest import mock

import pytest
import graphcigartoref as c
from tests.test_confirmed_realign_regressions import (
    Reader, repeat_pair, dna, cigar, calls, check_spelling,
)
from tests.test_tandem_real_trio import verify_spelling


@pytest.mark.parametrize('kind', ['I', 'D'])
@pytest.mark.parametrize('strand', ['+', '-'])
def test_decomposition_and_di_routes_produce_identical_vcf_positions(kind, strand):
    _, _, _, ref, qry = repeat_pair()
    bodies = ['300=37I37=74I111=300=', '411=37I37=74I300=',
              '361=111I387=', '300=111I448=']
    if kind == 'D':
        ref, qry = qry, ref
        bodies = [s.replace('I', 'D') for s in bodies]
    bodies.append(cigar(c._realign_long_di_sequences(ref, qry)))
    reader = Reader('N' * 200 + ref + 'N' * 100)
    signatures = []
    outputs = []
    for body in bodies:
        oriented_ref, oriented_query = ref, qry
        if strand == '-':
            body = ''.join(f'{t.n}{t.op}' for t in reversed(c.parse_cigar_ops(body)))
            oriented_ref, oriented_query = c.revcomp(ref), c.revcomp(qry)
        body = c.query_only_graph_cigar('>chr:' + c._tm_format_pairwise_ops(
            c._tm_ops_from_cigar_on_sequences(oriented_ref, oriented_query, body))).split(':', 1)[1]
        text = ('>chr:200H' + body + '100H' if strand == '+'
                else '<chr:100H' + body + '200H')
        result = c._rescue_serialized_di_pairs(text, oriented_query, reader)
        verify_spelling(unittest.TestCase(), result, oriented_query, reader)
        assert c._rescue_serialized_di_pairs(result, oriented_query, reader) == result
        signatures.append(calls(reader.sequence, oriented_query, result))
        outputs.append(result)
    assert all(result == outputs[0] for result in outputs)
    assert all(signature == signatures[0] for signature in signatures)
    assert len(signatures[0]) == 1
    assert signatures[0][0][1:] == ((500, 500, 111, 0) if kind == 'I'
                                   else (500, 611, 0, 111))


@pytest.mark.parametrize('kind', ['I', 'D'])
@pytest.mark.parametrize('reverse', [False, True])
def test_long_repeat_normalization_rotates_payloads_and_is_idempotent(kind, reverse):
    rng = random.Random(431)
    for trial in range(20):
        motif = dna(rng, rng.randrange(17, 45))
        # Make a tract much longer than the 50-column consolidation flank.
        left = dna(rng, 100) + next(b for b in 'ACGT' if b != motif[-1])
        right = next(b for b in 'ACGT' if b != motif[0]) + dna(rng, 100)
        repeats, delta = 70, 4
        ref, qry = left + motif * repeats + right, left + motif * (repeats + delta) + right
        if kind == 'D':
            ref, qry = qry, ref
        normalized = []
        for offset in (0, 23 * len(motif) + 3, repeats * len(motif)):
            pos = len(left) + offset
            tail = len(ref) - pos - (delta * len(motif) if kind == 'D' else 0)
            parts = [(pos, '='), (delta * len(motif), kind), (tail, '=')]
            r, q = ref, qry
            if reverse:
                r, q = c.revcomp(ref), c.revcomp(qry)
                parts.reverse()
            ops = c._tm_ops_from_cigar_on_sequences(r, q, ''.join(f'{n}{op}' for n, op in parts))
            result = c._normalize_repeat_indels(ops, r, q, reverse=reverse)
            check_spelling(r, q, result)
            assert c.alignment_score(result) == c.alignment_score(ops)
            assert c.di_repair_score(result) == c.di_repair_score(ops)
            assert c._normalize_repeat_indels(result, r, q, reverse=reverse) == result
            normalized.append(result)
        assert normalized[0] == normalized[1] == normalized[2]


def test_normalization_stops_at_mismatch_and_ambiguous_base():
    for middle, op in [('C', '1XG'), ('N', '1=')]:
        ref = 'T' + 'A' * 79 + middle + 'A' * 200
        qry = ref[:80] + ('G' if middle == 'C' else middle) + ref[81:150] + 'A' * 64 + ref[150:]
        ops = c._tm_ops_from_cigar_on_sequences(ref, qry, '80=' + op + '69=64I131=')
        result = c._normalize_repeat_indels(ops, ref, qry)
        check_spelling(ref, qry, result)
        expected = '80=1X64I200=' if middle == 'C' else '81=64I200='
        assert cigar(result) == expected
    # Pure insertion/deletion and segment-edge gaps cannot leave the segment.
    for ref, qry, body in [('', 'A' * 64, '64I'), ('A' * 64, '', '64D'),
                           ('A' * 20, 'A' * 84, '64I20=')]:
        ops = c._tm_ops_from_cigar_on_sequences(ref, qry, body)
        assert c._normalize_repeat_indels(ops, ref, qry) == ops


def test_normalization_does_not_assume_sam_m_means_exact_match():
    ref, qry = 'A' * 100, 'G' + 'A' * 163
    ops = [(50, 'M', ''), (64, 'I', 'A' * 64), (50, 'M', '')]
    result = c._normalize_repeat_indels(ops, ref, qry)
    check_spelling(ref, qry, result)
    assert cigar(result) == '1X64I99='


@pytest.mark.parametrize('kind', ['I', 'D'])
@pytest.mark.parametrize('reverse', [False, True])
def test_adjacent_equivalent_gaps_are_normalized_after_joining(kind, reverse):
    ref, qry = 'T' + 'A' * 100 + 'C', 'T' + 'A' * 220 + 'C'
    if kind == 'D':
        ref, qry = qry, ref
    body = f'31=60{kind}30=60{kind}41='
    expected = f'1=120{kind}101='
    if reverse:
        ref, qry = c.revcomp(ref), c.revcomp(qry)
        body = ''.join(f'{t.n}{t.op}' for t in reversed(c.parse_cigar_ops(body)))
        expected = f'101=120{kind}1='
    ops = c._tm_ops_from_cigar_on_sequences(ref, qry, body)
    result = c._normalize_repeat_indels(ops, ref, qry, reverse=reverse)
    check_spelling(ref, qry, result)
    assert cigar(result) == expected
    assert c._normalize_repeat_indels(result, ref, qry, reverse=reverse) == result


@pytest.mark.parametrize('strand', ['+', '-'])
def test_late_encoded_segment_keeps_nested_111bp_call_and_parent(strand):
    _, _, _, ref, qry = repeat_pair()
    chrom = 'A' * 50 + 'T' * 40 + 'N' * 110 + ref + 'N' * 100 + qry[300:411] + 'N' * 50
    if strand == '-':
        chrom = c.revcomp(chrom)
    direction = '>' if strand == '+' else '<'
    reader = Reader(chrom)
    st = c.build_stage1('>p:748M', '>before:50M>p:300M448H>dup:111M>p:300H448M>after:40M', 'r', 'q')
    body = c.query_only_graph_cigar('>chr:' + c._tm_format_pairwise_ops(
        c._tm_ops_from_cigar_on_sequences(ref, qry, '411=37I37=74I300='))).split(':', 1)[1]
    text = f'{direction}chr:50={direction}chr:200H{body}{len(chrom)-948}H{direction}40={len(chrom)-90}H'
    query = 'A' * 50 + qry + 'T' * 40
    before = c.parse_reference_cigar_segments(text, 'before')[1]
    with mock.patch.object(c, 'encode_large_insertions_in_piece', wraps=c.encode_large_insertions_in_piece) as encode:
        result = c._rescue_serialized_di_pairs(text, query, reader, stage1=st,
            graph_sequences={'p': ref, 'dup': qry[300:411]},
            graph_mappings={'p': f'{direction}chr:200H748={len(chrom)-948}H',
                            'dup': f'{direction}chr:1048H111=50H'})
    encode.assert_not_called()
    verify_spelling(unittest.TestCase(), result, query, reader)
    encoded = [s for s in c.parse_reference_cigar_segments(result, 'after') if not s.is_main]
    assert len(encoded) == 1
    assert (encoded[0].start, encoded[0].end, encoded[0].main_anchor) == (before.start, before.end, before.main_anchor)
    called = calls(chrom, query, result)
    assert sorted((ev[0], ev[3], ev[4]) for ev in called) == [('INS', 859, 0), ('INS_INS', 111, 0)]


def test_late_main_encoding_uses_same_secondary_lcs_as_stage4():
    ref = dna(random.Random(310), 1000)
    qry = ref[:500] + ref[100:400] + ref[500:]
    stage1 = c.build_stage1('>p:1000M', '>p:500M500H>p:100H300M600H>p:500H500M', 'r', 'q')
    backbone = c.Coord('chr', 100, 1100, '+')
    body = '500=100D400I' + qry[500:900] + '400='
    reader = Reader('N' * 100 + ref + 'N' * 100)
    with mock.patch.object(c, 'build_secondary_duplicate_lcs_pieces', wraps=c.build_secondary_duplicate_lcs_pieces) as lcs:
        late = c._rescue_serialized_di_pairs('>chr:100H' + body + '100H', qry, reader,
            stage1=stage1, graph_sequences={'p': ref}, ref_backbone_coord=backbone)
    lcs.assert_called_once_with(stage1, {'p': ref}, backbone)
    piece = c.PairwisePiece('main', 'p', 0, 0, 0, 1000, 0, 1000, body)
    stage4 = c.build_stage4_linear(stage1, c.Stage3PolishResult([piece], body), {'p': ref}, {},
        {'chr': 1200}, {}, reader, ref_backbone_coord=backbone,
        record_sequences={'q': qry}, query_name='q')
    early = c.serialize_stage3_reflalign(stage4, {'chr': 1200}, ref_reader=reader,
                                         main_anchor_coord=backbone, query_sequence=qry)
    assert early == late == '>chr:100H500=>chr:200H300=700H>500=100H'
    verify_spelling(unittest.TestCase(), late, qry, reader)
    assert calls(reader.sequence, qry, early) == calls(reader.sequence, qry, late)
