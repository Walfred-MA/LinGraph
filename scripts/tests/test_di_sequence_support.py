"""D/I rescue requires local homology and keeps its separate affine score."""
import random
import unittest
from unittest import mock

import pytest
import graphcigartoref as c
from tests.test_confirmed_realign_regressions import Reader, calls, dna
from tests.test_late_di_rescue_score import fixture, ops
from tests.test_tandem_real_trio import verify_spelling


SHAPES = [(500, 40), (500, 120), (300, 300), (900, 900),
          (40, 500), (120, 500), (120, 120)]


def uncompressed_supported_matches(operations, width):
    columns = [op == '=' for n, op, _ in operations for _ in range(n)]
    covered = [False] * len(columns)
    for start in range(len(columns) - width + 1):
        if 5 * sum(columns[start:start + width]) >= 4 * width:
            covered[start:start + width] = [True] * width
    return sum(match and seen for match, seen in zip(columns, covered))


def test_compressed_support_matches_literal_alignment_columns():
    rng = random.Random(4152)
    for _ in range(200):
        operations = [(rng.randint(1, 300), rng.choice('=XID'), '')
                      for _ in range(rng.randint(1, 20))]
        width = rng.choice((16, 20, 25, 31))
        assert c._di_supported_match_bases(operations, width) == uncompressed_supported_matches(operations, width)
        split = [(1, op, '') for n, op, _ in operations for _ in range(n)]
        assert c._di_supported_match_bases(split, width) == c._di_supported_match_bases(operations, width)
        assert c._di_supported_match_bases(list(reversed(operations)), width) == c._di_supported_match_bases(operations, width)


def test_support_counts_small_gaps_and_never_double_counts_matches():
    assert c._di_supported_match_bases([(6, '=', ''), (1, 'D', '')] * 8) == 48
    assert c._di_supported_match_bases([(6, '=', ''), (2, 'D', '')] * 8) == 0
    assert c._di_supported_match_bases([(30, '=', '')]) == 0
    assert c._di_supported_match_bases([(31, '=', '')]) == 31
    assert c._di_supported_match_bases([(10**8, '=', ''), (10**8, 'D', ''), (31, '=', '')]) == 10**8 + 31


@pytest.mark.parametrize('order', ['DI', 'ID'])
def test_a_short_matching_island_cannot_justify_shattering_an_unrelated_pair(order):
    rng = random.Random(4081)
    reference, query = dna(rng, 900), dna(rng, 900)
    query = query[:400] + reference[400:440] + query[440:]
    candidate = c._global_affine_payload_align(
        reference, query, match_score=c.DI_MATCH_SCORE, mismatch_score=c.DI_MISMATCH_SCORE,
        gap_open=c.DI_GAP_OPEN, gap_extend=c.DI_GAP_EXTEND)
    assert c._di_supported_match_bases(candidate) >= 40
    assert c._realign_long_di_sequences(reference, query) is None
    raw = {'D': '900D', 'I': '900I' + query}
    original = c.parse_cigar_ops(''.join(raw[axis] for axis in order))
    assert c.di_repair_score(candidate) > c.di_repair_score(original)


@pytest.mark.parametrize('gap', ['I', 'D'])
@pytest.mark.parametrize('position', ['left', 'middle', 'right'])
def test_short_homologous_copy_with_a_long_gap_still_repairs(gap, position):
    rng = random.Random(505)
    shared, extra = dna(rng, 40), dna(rng, 460)
    split = {'left': 0, 'middle': 20, 'right': 40}[position]
    longer = shared[:split] + extra + shared[split:]
    reference, query = (shared, longer) if gap == 'I' else (longer, shared)
    candidate = c._realign_long_di_sequences(reference, query)
    assert candidate is not None
    original = c.parse_cigar_ops(f'{len(reference)}D{len(query)}I' + query)
    assert c._di_rescue_is_valid(original, candidate, reference, query)


def unrelated_pair(nr, nq, trial):
    rng = random.Random(nr * 100000 + nq * 100 + trial)
    r, q = dna(rng, nr), dna(rng, nq)
    left, right = dna(rng, 200), dna(rng, 200)
    return r, q, left, right


@pytest.mark.parametrize('nr,nq', SHAPES)
@pytest.mark.parametrize('order', ['DI', 'ID'])
@pytest.mark.parametrize('strand', ['+', '-'])
def test_unrelated_pairs_keep_variants_across_rescue_paths(nr, nq, order, strand):
    for trial in range(30):
        r, q, left, right = unrelated_pair(nr, nq, trial)
        raw = {'D': f'{nr}D', 'I': f'{nq}I' + q}
        body = '200=' + ''.join(raw[op] for op in order) + '200='
        reference, query = left + r + right, left + q + right
        candidate = c._global_affine_payload_align(
            r, q, match_score=c.DI_MATCH_SCORE, mismatch_score=c.DI_MISMATCH_SCORE,
            gap_open=c.DI_GAP_OPEN, gap_extend=c.DI_GAP_EXTEND)
        original = c.parse_cigar_ops(''.join(raw[op] for op in order))
        # The affine objective alone favors chance matches; the homology
        # filter must reject this candidate without adding an SV-score gate.
        assert candidate is not None
        assert c.di_repair_score(candidate) > c.di_repair_score(original)
        assert c.alignment_score(candidate) < c.alignment_score(original)
        assert c._di_rescue_is_valid(original, candidate, r, q)
        assert c._realign_long_di_sequences(r, q) is None
        args = fixture(body, reference, query, strand)
        piece, linear = c._rescue_long_di_pairs_in_linear_piece(*args)
        assert piece is args[0] and linear is args[1]
        piece, linear = c._rescue_internal_di_windows_multiscale(*args)
        assert piece is args[0] and linear is args[1]
        direction = '>' if strand == '+' else '<'
        text = direction + 'chr:100H' + body + '100H'
        result = c._rescue_serialized_di_pairs(text, query, args[-1])
        verify_spelling(unittest.TestCase(), result, query, args[-1])
        segments = c.parse_reference_cigar_segments(result, 'test')
        tokens = [tok for segment in segments for tok in segment.ops if tok.op != 'H']
        assert c.alignment_score(tokens) >= c.alignment_score(ops(body))
        assert any(tok.op in 'ID' and tok.n > 50 for tok in tokens)
        assert calls(args[-1].sequence, query, result)
        assert sum(c.pairwise_ref_span(''.join(c._op_to_cigar(t) for t in seg.ops
                                             if t.op != 'H')) for seg in segments) == len(reference)
        assert c.pairwise_query_span(''.join(c._op_to_cigar(t) for t in tokens)) == len(query)


@pytest.mark.parametrize('n', [40, 70, 120, 300, 900, 1500])
@pytest.mark.parametrize('strand', ['+', '-'])
def test_ten_percent_diverged_homologous_repairs_still_pass(n, strand):
    for trial in range(10):
        rng = random.Random(n * 100 + trial)
        reference = dna(rng, n)
        query = list(reference)
        for i in rng.sample(range(n), n // 10):
            query[i] = rng.choice([base for base in 'ACGT' if base != query[i]])
        query = ''.join(query)
        original = c.parse_cigar_ops(f'{n}D{n}I' + query)
        candidate = c._realign_long_di_sequences(reference, query)
        assert candidate is not None
        assert c.alignment_score(candidate) > c.alignment_score(original)
        assert c.di_repair_score(candidate) > c.di_repair_score(original)
        assert c._di_rescue_is_valid(original, candidate, reference, query)
        args = fixture('40=' + f'{n}D{n}I' + query + '40=',
                       'A' * 40 + reference + 'C' * 40,
                       'A' * 40 + query + 'C' * 40, strand)
        piece, linear = c._rescue_long_di_pairs_in_linear_piece(*args)
        assert piece is not args[0]
        assert c.pairwise_query_span(piece.cigar) == n + 80
        assert c.pairwise_ref_span(piece.cigar) == n + 80
        assert (linear.genome_start, linear.genome_end, linear.strand) == (100, n + 180, strand)
        text = ('>' if strand == '+' else '<') + 'chr:100H' + piece.cigar + '100H'
        verify_spelling(unittest.TestCase(), c.query_only_graph_cigar(text), args[3][17:], args[-1])


@pytest.mark.parametrize('native', [False, True])
def test_native_and_fallback_random_alignment_cannot_bypass_acceptance(native):
    if native and c._parasail is None:
        pytest.skip('parasail not installed')
    r, q, _, _ = unrelated_pair(120, 120, 0)
    original = c.parse_cigar_ops('120D120I' + q)
    with mock.patch.object(c, '_parasail', c._parasail if native else None):
        candidate = c._realign_long_di_sequences(r, q)
    assert candidate is None


@pytest.mark.parametrize('n', [1000, 5000])
def test_long_mapper_candidates_require_local_homology(n):
    rng = random.Random(n)
    reference, query = dna(rng, n), dna(rng, n)
    # Isolate acceptance from mapper seeding: exercise a valid but fragmented
    # candidate which a mapper could propose for an unrelated long payload.
    candidate = c._global_affine_payload_align(
        reference, query, match_score=c.DI_MATCH_SCORE,
        mismatch_score=c.DI_MISMATCH_SCORE,
        gap_open=c.DI_GAP_OPEN, gap_extend=c.DI_GAP_EXTEND)
    original = c.parse_cigar_ops(f'{n}D{n}I' + query)
    with mock.patch.object(c, 'minimap2_payload_ops', return_value=candidate) as mapper:
        proposed = c._realign_long_di_sequences(reference, query)
    mapper.assert_called_once()
    assert c.di_repair_score(candidate) > c.di_repair_score(original)
    assert proposed is None


@pytest.mark.parametrize('strand', ['+', '-'])
def test_balanced_replacement_survives_stage4_encoding_and_calling(strand):
    r, q, left, right = unrelated_pair(300, 300, 0)
    reference, query = left + r + right, left + q + right
    body = '200=300D300I' + q + '200='
    reader = Reader('N' * 100 + (reference if strand == '+' else c.revcomp(reference))
                    + 'N' * 100)
    backbone = c.Coord('chr', 100, 100 + len(reference), strand)
    stage1 = c.build_stage1(f'>p:{len(reference)}M', '>p:' + body.replace('=', 'M'), 'r', 'q')
    piece = c.PairwisePiece('main', 'p', 0, 0, 0, len(reference), 0, len(reference), body)
    result = c.build_stage4_linear(
        stage1, c.Stage3PolishResult([piece], body), {'p': reference}, {},
        {'chr': len(reader.sequence)}, {}, reader, ref_backbone_coord=backbone,
        record_sequences={'q': query}, query_name='q')
    text = c.serialize_stage3_reflalign(
        result, {'chr': len(reader.sequence)}, ref_reader=reader,
        main_anchor_coord=backbone, query_sequence=query)
    text = c.query_only_graph_cigar(text)
    verify_spelling(unittest.TestCase(), text, query, reader)
    assert calls(reader.sequence, query, text)
