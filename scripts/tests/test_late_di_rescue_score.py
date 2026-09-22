"""Late D/I rescues must not overwrite a better or equally scored alignment."""
import math
import unittest
from unittest import mock

import graphcigartoref as core


class Reader:
    def __init__(self, sequence):
        self.sequence = sequence
        self.index = {'chr': (len(sequence),)}

    def fetch(self, chrom, start, end, strand='+'):
        assert chrom == 'chr'
        sequence = self.sequence[start:end]
        return core.revcomp(sequence) if strand == '-' else sequence


def ops(cigar):
    return [(t.n, t.op, t.payload) for t in core.parse_cigar_ops(cigar)]


def fixture(cigar, reference, query, strand='+'):
    length = len(reference)
    piece = core.PairwisePiece('main', 'node', 0, 0, 0, length, 0, length, cigar)
    linear = core.LinearPiece(piece, 'chr', strand, 100, 100 + length, cigar)
    chromosome = reference if strand == '+' else core.revcomp(reference)
    reader = Reader('N' * 100 + chromosome + 'N' * 100)
    # Exercise the query offset used when a piece follows an encoded insertion.
    return piece, linear, 17, 'C' * 17 + query, reader


class CompositeRescueTests(unittest.TestCase):
    original = '40=100D20=100I' + 'A' * 100 + '40='
    candidates = {
        'tie': [(100, 'I', 'A' * 100), (20, '=', ''), (100, 'D', '')],
        # 40 matches and 80 deleted/inserted bases, but 80 gap openings:
        # supported sequence can still score worse than the original window.
        'worse': [(1, '=', ''), (2, 'D', ''), (2, 'I', 'AA')] * 40,
        'improved': [(120, '=', '')],
        'wrong_span': [(119, '=', '')],
    }

    def setUp(self):
        gate = mock.patch.object(core, '_realignment_hardness_allows', return_value=True)
        gate.start()
        self.addCleanup(gate.stop)

    def test_internal_windows_reject_ties_worse_scores_and_changed_spans(self):
        for strand in ('+', '-'):
            for label in ('tie', 'worse', 'wrong_span'):
                with self.subTest(strand=strand, candidate=label):
                    args = fixture(self.original, 'A' * 200, 'A' * 200, strand)
                    with mock.patch.object(core, '_tm_realign_secondary_sv_window',
                                           return_value=self.candidates[label]) as align:
                        piece, linear = core._rescue_internal_di_windows_in_linear_piece(*args)
                    self.assertEqual(align.call_args.args, ('A' * 120, 'A' * 120))
                    self.assertIs(piece, args[0])
                    self.assertIs(linear, args[1])
                    self.assertEqual(piece.cigar, self.original)

    def test_internal_windows_accept_improvement_on_both_strands(self):
        for strand in ('+', '-'):
            with self.subTest(strand=strand):
                args = fixture(self.original, 'A' * 200, 'A' * 200, strand)
                with mock.patch.object(core, '_tm_realign_secondary_sv_window',
                                       return_value=self.candidates['improved']):
                    piece, linear = core._rescue_internal_di_windows_in_linear_piece(*args)
                self.assertEqual(piece.cigar, '200=')
                self.assertEqual((linear.genome_start, linear.genome_end, linear.strand),
                                 (100, 300, strand))
                self.assertEqual(core.pairwise_query_span(piece.cigar), 200)

    def test_multiscale_pass_cannot_fragment_a_previously_improved_window(self):
        cigar = '1000=100D20=100I' + 'A' * 100 + '1000='
        args = fixture(cigar, 'A' * 2120, 'A' * 2120)
        # A broad pass ties, and the later local pass proposes a worse result.
        with mock.patch.object(core, '_tm_realign_secondary_sv_window', side_effect=[
            self.candidates['tie'], self.candidates['worse'],
        ]) as align:
            piece, linear = core._rescue_internal_di_windows_multiscale(*args)
        self.assertEqual(align.call_count, 2)
        self.assertIs(piece, args[0])
        self.assertIs(linear, args[1])


class AdjacentPairRescueTests(unittest.TestCase):
    def test_both_pair_orders_and_sizes_keep_tied_or_worse_alignments(self):
        for n in (5, 60):
            function = '_realign_small_di_sequences' if n <= 30 else '_realign_long_di_sequences'
            for order in ('DI', 'ID'):
                gap = {'D': f'{n}D', 'I': f'{n}I' + 'A' * n}
                cigar = '40=' + ''.join(gap[op] for op in order) + '40='
                candidates = {
                    'tie': [(n, op, 'A' * n if op == 'I' else '') for op in order[::-1]],
                    'worse': [(k, op, 'A' * k if op == 'I' else '')
                              for k in (n // 2, n - n // 2) for op in 'DI'],
                    'wrong_span': [(n + 1, '=', '')],
                }
                for label, candidate in candidates.items():
                    with self.subTest(n=n, order=order, candidate=label):
                        args = fixture(cigar, 'A' * (80 + n), 'A' * (80 + n))
                        with mock.patch.object(core, function, return_value=candidate) as align:
                            piece, linear = core._rescue_long_di_pairs_in_linear_piece(*args)
                        align.assert_called_once_with('A' * n, 'A' * n)
                        self.assertIs(piece, args[0])
                        self.assertIs(linear, args[1])

    def test_valid_small_and_long_sequence_rescues_still_work(self):
        for n in (3, 60, 1000):
            for strand in ('+', '-'):
                with self.subTest(n=n, strand=strand):
                    reference = 'A' * (80 + n)
                    cigar = '40=' + f'{n}D{n}I' + 'A' * n + '40='
                    args = fixture(cigar, reference, reference, strand)
                    piece, linear = core._rescue_long_di_pairs_in_linear_piece(*args)
                    self.assertEqual(piece.cigar, f'{80 + n}=')
                    self.assertEqual((linear.genome_start, linear.genome_end), (100, 180 + n))

    def test_successful_small_mismatch_rescue_preserves_query_payload(self):
        args = fixture('40=3D3IATG40=', 'A' * 40 + 'ACG' + 'A' * 40,
                       'A' * 40 + 'ATG' + 'A' * 40)
        piece, _ = core._rescue_long_di_pairs_in_linear_piece(*args)
        self.assertEqual(piece.cigar, '41=1XT41=')

    def test_insufficient_sequence_support_is_still_rejected(self):
        cigar = '40=1000D1000I' + 'C' * 1000 + '40='
        args = fixture(cigar, 'A' * 1080, 'A' * 40 + 'C' * 1000 + 'A' * 40)
        with mock.patch.object(core, 'minimap2_payload_ops',
                               return_value=[(1000, 'X', 'C' * 1000)]):
            piece, linear = core._rescue_long_di_pairs_in_linear_piece(*args)
        self.assertIs(piece, args[0])
        self.assertIs(linear, args[1])


class ScoringTests(unittest.TestCase):
    def test_di_affine_score_is_separate_from_sv_score(self):
        for size in (1, 2, 49, 50, 51, 1000):
            for axis in ('I', 'D'):
                with self.subTest(size=size, axis=axis):
                    self.assertEqual(core.di_repair_score([(size, axis, '')]),
                                     -(4 + size - 1))
        self.assertEqual(core.di_repair_score(ops('10=2XAA3IAAA2IAA4D')),
                         10 - 8 - 8 - 7)
        self.assertEqual(core.di_repair_score(ops('10=2XAA5IAAAAA4D')),
                         10 - 8 - 8 - 7)
        # Matches and mismatches are charged separately in D/I acceptance.
        original = core.parse_cigar_ops('20D20I' + 'A' * 20)
        candidate = [(10, '=', ''), (10, 'X', 'C' * 10)]
        self.assertTrue(core._di_rescue_is_valid(
            original, candidate, 'A' * 20, 'A' * 10 + 'C' * 10))
        # Excessive gap openings can still make an exact repair score worse.
        original = core.parse_cigar_ops('40=20D20I' + 'A' * 20 + '40=')
        candidate = [(1, '=', ''), (1, 'D', ''), (1, 'I', 'A')] * 50
        self.assertFalse(core._di_rescue_is_valid(original, candidate, 'A' * 100, 'A' * 100))

    def test_small_di_dp_uses_minus_four_for_mismatches(self):
        candidate = core._realign_small_di_sequences('ACAC', 'CACA')
        # Four mismatches cost 16; shifting recovers three matches for two
        # one-base gap openings, yielding 3 - 8 = -5.
        self.assertEqual(core.di_repair_score(candidate), -5)
        self.assertEqual(sum(n for n, op, _ in candidate if op == '='), 3)

    def test_every_gap_extension_is_logarithmic_with_unchanged_open_cost(self):
        for size in (1, 49, 50, 51, 100, 1000, 100000):
            for axis in ('I', 'D'):
                with self.subTest(size=size, axis=axis):
                    extension = 100 * math.log2(size)
                    self.assertEqual(core._tm_realignment_score([(size, axis, '')]),
                                     -(12 + extension))

    def test_logarithmic_score_coalesces_gaps_before_scoring(self):
        split = [(25, 'I', ''), (26, 'I', ''), (60, 'D', ''), (40, 'D', '')]
        whole = [(51, 'I', ''), (100, 'D', '')]
        expected = -24 - 100 * math.log2(51) - 100 * math.log2(100)
        self.assertAlmostEqual(core._tm_realignment_score(split), expected)
        self.assertEqual(core._tm_realignment_score(split), core._tm_realignment_score(whole))

    def test_equal_scores_are_independent_of_gap_order_and_match_tokenization(self):
        gaps = [(n, 'I', '') for n in (1, 50, 51, 123, 2345, 100000)]
        forward = [op for gap in gaps for op in (gap, (3, '=', ''))]
        reverse = [op for gap in reversed(gaps) for op in
                   (gap, (1, '=', ''), (2, '=', ''))]
        self.assertEqual(core._tm_realignment_score(forward), core._tm_realignment_score(reverse))

    def test_standard_score_counts_runs_not_cigar_tokens(self):
        # Adjacent I tokens are one gap, but an immediately following D is new.
        self.assertAlmostEqual(core._tm_realignment_score(ops('10=2XAA3IAAA2IAA4D')),
                         40 - 16 - 24 - 100 * (math.log2(5) + math.log2(4)))
        self.assertAlmostEqual(core._tm_realignment_score(ops('10=2XAA5IAAAAA4D')),
                         40 - 16 - 24 - 100 * (math.log2(5) + math.log2(4)))

    def test_one_base_gaps_can_outscore_two_long_gaps(self):
        original = core.parse_cigar_ops('40=20D20I' + 'A' * 20 + '40=')
        candidate = [(19, '=', ''), (1, 'D', ''), (1, 'I', 'A')] * 5
        # log2(1) has no extension cost: this requested objective does not
        # guarantee fewer gaps, even though adjacent tokens still coalesce.
        self.assertEqual(core._tm_realignment_score(candidate), 260)
        self.assertAlmostEqual(core._tm_realignment_score(
            [(t.n, t.op, t.payload) for t in original]), 320 - 24 - 200 * math.log2(20))
        self.assertTrue(core._di_rescue_is_valid(
            original, candidate, 'A' * 100, 'A' * 100))

    @unittest.skipIf(core._parasail is None, 'parasail not installed')
    def test_native_affine_candidate_is_rescored_by_standard_objective(self):
        reference, query = 'ACGTTTTTCAGT', 'ACGTCAGT'
        candidate = core._sv_window_affine_align(reference, query)
        self.assertEqual(sum(n for n, op, _ in candidate if op in 'M=X'), 8)
        self.assertAlmostEqual(core._tm_realignment_score(candidate),
                               32 - 12 - 100 * math.log2(4))

    def test_gap_open_cost_includes_the_previous_operation(self):
        original = core.parse_cigar_ops('1IA3=')
        candidate = [(3, '=', ''), (1, 'I', 'C')]
        # Old = annotations are recomputed: AAC aligned to AAA has a mismatch.
        self.assertTrue(core._di_rescue_is_valid(original, candidate, 'AAA', 'AAAC'))
        # Correcting a mismatch gains 5 points, exceeding the additional
        # 3-point cost of splitting a joined gap into two openings.
        self.assertTrue(core._di_rescue_is_valid(
            original, candidate, 'AAA', 'AAAC', left_context=core.parse_cigar_ops('1IA')))

    def test_gap_open_cost_includes_the_following_operation(self):
        original = core.parse_cigar_ops('3=1IA')
        candidate = [(1, 'I', 'C'), (3, '=', '')]
        self.assertTrue(core._di_rescue_is_valid(original, candidate, 'AAA', 'CAAA'))
        self.assertTrue(core._di_rescue_is_valid(
            original, candidate, 'AAA', 'CAAA', right_context=core.parse_cigar_ops('1IA')))

    def test_candidate_match_labels_are_checked_against_actual_bases(self):
        original = core.parse_cigar_ops('1IA3=')
        # Raw op scores tie; actual sequence comparison makes the candidate worse.
        self.assertFalse(core._di_rescue_is_valid(
            original, [(3, '=', ''), (1, 'I', 'A')], 'AAA', 'CAAA'))

    def test_boundary_scoring_uses_complete_gap_length_on_both_sides(self):
        for size in (1, 50, 1000):
            for side in ('left', 'right'):
                with self.subTest(size=size, side=side):
                    left = side == 'left'
                    original = core.parse_cigar_ops('1IA3=' if left else '3=1IA')
                    candidate = ([(3, '=', ''), (1, 'I', 'A')] if left else
                                 [(1, 'I', 'A'), (3, '=', '')])
                    query = 'AAAA'
                    # Split context tokens must still be scored as one full gap.
                    context = [core.CigarOp(1, 'I', 'A')]
                    if size > 1:
                        context.append(core.CigarOp(size - 1, 'I', 'A' * (size - 1)))
                    self.assertFalse(core._di_rescue_is_valid(
                        original, candidate, 'AAA', query,
                        **{side + '_context': context}))

    def test_more_svs_are_allowed_when_di_score_improves(self):
        original = core.parse_cigar_ops('10=980D980I' + 'A' * 980 + '10=')
        candidate = [(60, 'D', ''), (60, 'I', 'A' * 60), (440, '=', '')] * 2
        self.assertEqual(core._tm_count_large_indels(candidate), 4)
        self.assertTrue(core._di_rescue_is_valid(original, candidate, 'A' * 1000, 'A' * 1000))

    def test_invalid_operations_and_spans_are_rejected(self):
        original = core.parse_cigar_ops('200=')
        for candidate in ([(0, 'I', ''), (200, '=', '')], [(200, 'S', '')], [(199, '=', '')]):
            self.assertFalse(core._di_rescue_is_valid(original, candidate, 'A' * 200, 'A' * 200))



class SerializedRescueTests(unittest.TestCase):
    def test_final_sv_consolidation_does_not_use_the_di_affine_gate(self):
        # Isolate replacement acceptance from score-neutral canonical shifts.
        body = '40=100D20=100I' + 'A' * 100 + '40='
        text = '>chr:100H' + body + '100H'
        reader = Reader('N' * 100 + 'A' * 200 + 'N' * 100)
        candidate = [(40, '=', '')] + [(3, '=', ''), (2, 'D', ''),
                     (2, 'I', 'AA')] * 24 + [(40, '=', '')]
        self.assertTrue(core._di_rescue_is_valid(
            core.parse_cigar_ops(body), candidate, 'A' * 200, 'A' * 200))
        self.assertLess(core._tm_realignment_score(candidate),
                        core._tm_realignment_score(ops(body)))
        unchanged = lambda piece, linear, *_: (piece, linear)
        with mock.patch.object(core, '_SV_REALIGNMENT_ENABLED', True), \
             mock.patch.object(core, '_rescue_internal_di_windows_multiscale', side_effect=unchanged), \
             mock.patch.object(core, '_rescue_long_di_pairs_in_linear_piece', side_effect=unchanged), \
             mock.patch.object(core, '_normalize_repeat_indels', side_effect=lambda ops, *_a, **_k: ops), \
             mock.patch.object(core, '_tm_consolidate_sv_ops', return_value=candidate):
            result = core._rescue_serialized_di_pairs(text, 'A' * 200, reader)
        self.assertEqual(result, text)

    def test_serialized_composite_rescue_keeps_tie_and_accepts_improvement(self):
        body = CompositeRescueTests.original
        text = '>chr:100H' + body + '100H'
        reader = Reader('N' * 100 + 'A' * 200 + 'N' * 100)
        for label in ('tie', 'worse', 'improved'):
            with self.subTest(candidate=label), \
                 mock.patch.object(core, '_SV_REALIGNMENT_ENABLED', True), \
                 mock.patch.object(core, '_tm_realign_secondary_sv_window',
                                   return_value=CompositeRescueTests.candidates[label]):
                with mock.patch.object(core, '_realignment_hardness_allows', return_value=True), \
                     mock.patch.object(core, '_tm_consolidate_sv_ops', side_effect=lambda ops, *_: ops), \
                     mock.patch.object(core, '_normalize_repeat_indels', side_effect=lambda ops, *_a, **_k: ops):
                    result = core._rescue_serialized_di_pairs(text, 'A' * 200, reader)
            self.assertEqual(result, '>chr:100H200=100H' if label == 'improved' else text)

    def test_final_consolidation_requires_a_strict_improvement(self):
        # Equal-score placement normalization is tested separately; the
        # realignment candidate itself must still pass the strict score gate.
        body = '40=60I' + 'A' * 60 + '20=60D20=5IAAAAA40='
        text = '>chr:100H' + body + '100H'
        reader = Reader('N' * 100 + 'A' * 180 + 'N' * 100)
        candidates = {
            'tie': ops('40=60D20=60I' + 'A' * 60 + '20=5IAAAAA40='),
            'worse': [(40, '=', '')] + [(1, 'D', ''), (1, 'I', 'A')] * 60
                     + [(80, '=', ''), (5, 'I', 'AAAAA')],
            'improved': [(180, '=', ''), (5, 'I', 'AAAAA')],
        }
        unchanged = lambda piece, linear, *_: (piece, linear)
        for label, candidate in candidates.items():
            with self.subTest(candidate=label), \
                 mock.patch.object(core, '_SV_REALIGNMENT_ENABLED', True), \
                 mock.patch.object(core, '_rescue_internal_di_windows_multiscale', side_effect=unchanged), \
                 mock.patch.object(core, '_rescue_long_di_pairs_in_linear_piece', side_effect=unchanged), \
                 mock.patch.object(core, '_normalize_repeat_indels', side_effect=lambda ops, *_a, **_k: ops), \
                 mock.patch.object(core, '_tm_consolidate_sv_ops', return_value=candidate) as consolidate:
                result = core._rescue_serialized_di_pairs(text, 'A' * 185, reader)
            self.assertEqual(consolidate.call_count, 1)
            self.assertEqual(result, '>chr:100H180=5IAAAAA100H' if label == 'improved' else text)


if __name__ == '__main__':
    unittest.main()
