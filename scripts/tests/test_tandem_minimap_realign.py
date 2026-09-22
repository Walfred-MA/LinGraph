"""Tandem realignment must pass query hardness and improve the linear score."""
import math
import os
from pathlib import Path
import random
import subprocess
import tempfile
from types import SimpleNamespace
import unittest
from unittest import mock

import graphcigartoref as core


def stage(cigar):
    piece = core.PairwisePiece('main', 'node', 0, 0, 0,
                              core.pairwise_ref_span(cigar), 0,
                              core.pairwise_ref_span(cigar), cigar)
    return core.Stage2Result([], [piece], cigar)


def repeated_stage():
    return core.build_stage1('>node:200M', '>node:150M50H>node:50H150M', 'r', 'q')


class SequenceAssertions(unittest.TestCase):
    def assert_spells(self, ops, ref, qry):
        r = 0
        spelled = []
        for n, op, payload in ops:
            self.assertGreater(n, 0)
            if op in 'M=':
                spelled.append(ref[r:r + n])
            elif op == 'X':
                self.assertEqual(len(payload), 2 * n)
                self.assertEqual(payload[:n].upper(), ref[r:r + n].upper())
                spelled.append(payload[n:])
            elif op == 'I':
                self.assertEqual(len(payload), n)
                spelled.append(payload)
            elif op != 'D':
                self.fail(op)
            if op in 'M=XD':
                r += n
            self.assertLessEqual(r, len(ref))
        self.assertEqual(r, len(ref))
        self.assertEqual(''.join(spelled).upper(), qry.upper())

    def hit(self, ref, qry, q0, q1, r0, r1, cigar):
        hit = core._tandem_hit_from_cigar(ref, qry, q0, q1, r0, r1, cigar)
        self.assertIsNotNone(hit)
        return hit


class IntervalSweepTests(SequenceAssertions):
    def test_highest_score_per_interval_and_tie_order(self):
        hit = core.TandemAlignmentHit
        hits = [hit(0, 10, 0, 10, (), 10), hit(4, 15, 4, 15, (), 20),
                hit(18, 20, 18, 20, (), 3)]
        self.assertEqual(core._select_tandem_hit_intervals(hits),
                         [(0, 4, 0), (4, 15, 1), (18, 20, 2)])
        ties = [hit(0, 10, 4, 14, (), 10), hit(0, 10, 0, 10, (), 10)]
        for order in (ties, list(reversed(ties))):
            selected = core._select_tandem_hit_intervals(order)
            self.assertEqual([(a, b, order[i].reference_start) for a, b, i in selected],
                             [(0, 10, 0)])

    def test_deletion_at_query_boundary_belongs_to_right_slice(self):
        ref, qry = 'ACGGTT', 'ACTT'
        hit = self.hit(ref, qry, 0, 4, 0, 6, '2M2D2M')
        self.assertEqual(core._slice_tandem_hit(hit, 0, 2), (0, [(2, '=', '')]))
        self.assertEqual(core._slice_tandem_hit(hit, 2, 4),
                         (2, [(2, 'D', ''), (2, '=', '')]))

    def test_slicing_inside_insertion_and_mismatch_preserves_payload_halves(self):
        hit = self.hit('AAAA', 'TTGGGG', 0, 6, 0, 4, '2I4X')
        self.assertEqual(core._slice_tandem_hit(hit, 1, 4),
                         (0, [(1, 'I', 'T'), (2, 'X', 'AAGG')]))
        self.assertEqual(core._slice_tandem_hit(hit, 4, 6),
                         (2, [(2, 'X', 'AAGG')]))

    def test_reference_overlap_retains_duplicated_query_bases(self):
        ref, qry = 'ACGTTGCA', 'ACGTTTGTGCA'
        hits = [self.hit(ref, qry, 0, 5, 0, 5, '5M'),
                self.hit(ref, qry, 5, 11, 2, 8, '6M')]
        ops = core._assemble_tandem_hit_intervals(ref, qry, hits)
        self.assert_spells(ops, ref, qry)
        self.assertEqual(sum(n for n, op, _ in ops if op == 'I'), 3)

    def test_scattered_hits_include_uncovered_query_and_reference(self):
        ref, qry = 'AACCGGTTAA', 'TACCGTAC'
        hits = [self.hit(ref, qry, 1, 4, 1, 4, '3M'),
                self.hit(ref, qry, 5, 7, 7, 9, '2M')]
        ops = core._assemble_tandem_hit_intervals(ref, qry, hits)
        self.assert_spells(ops, ref, qry)
        self.assertEqual(sum(n for n, op, _ in ops if op == 'D'), 5)
        self.assertEqual(sum(n for n, op, _ in ops if op == 'I'), 3)

    def test_nested_and_fully_backwards_hits_do_not_reconsume_reference(self):
        ref, qry = 'AACCGGTTAA', 'CCGGAATT'
        hits = [self.hit(ref, qry, 0, 4, 4, 8, '4M'),
                self.hit(ref, qry, 4, 8, 0, 4, '4M')]
        ops = core._assemble_tandem_hit_intervals(ref, qry, hits)
        self.assert_spells(ops, ref, qry)
        self.assertEqual(sum(n for n, op, _ in ops if op == 'I'), 4)

    def test_malformed_hits_are_rejected_without_padding(self):
        for args in ((0, 4, 0, 4, '3M'), (0, 4, 0, 4, '4MZ'),
                     (0, 5, 0, 4, '4M1I'), (0, 4, -1, 3, '4M'),
                     (0, 4, 0, 4, '1S3M'), (4, 4, 0, 4, '4D')):
            with self.subTest(args=args):
                self.assertIsNone(core._tandem_hit_from_cigar('ACGT', 'ACGT', *args))

    def test_random_scattered_cigars_preserve_all_bases(self):
        rng = random.Random(81264)
        for case in range(500):
            ref = ''.join(rng.choices('ACGT', k=150))
            qry = ''.join(rng.choices('ACGT', k=170))
            hits = []
            for _ in range(rng.randrange(1, 12)):
                tokens = [(rng.randrange(1, 8), rng.choice('MIDX'))
                          for _ in range(rng.randrange(1, 10))]
                tokens.append((3, 'M'))
                rn = sum(n for n, op in tokens if op in 'MDX')
                qn = sum(n for n, op in tokens if op in 'MIX')
                r0 = rng.randrange(len(ref) - rn + 1)
                q0 = rng.randrange(len(qry) - qn + 1)
                hits.append(self.hit(ref, qry, q0, q0 + qn, r0, r0 + rn,
                                     ''.join(f'{n}{op}' for n, op in tokens)))
            with self.subTest(case=case):
                ops = core._assemble_tandem_hit_intervals(ref, qry, hits)
                self.assert_spells(ops, ref, qry)


class ScoreGateTests(SequenceAssertions):
    def setUp(self):
        self.ref = 'ACGT' * 50
        self.qry = self.ref
        self.reader = SimpleNamespace(index={'chr': (500,)}, fetch=mock.Mock(return_value=self.ref))
        self.backbone = core.Coord('chr', 100, 300, '+')
        self.stage1 = repeated_stage()
        self.enabled = mock.patch.object(core, '_SV_REALIGNMENT_ENABLED', True)
        self.enabled.start()
        self.addCleanup(self.enabled.stop)
        # These tests isolate score acceptance; hardness is exercised separately.
        self.hardness = mock.patch.object(core, '_realignment_hardness_allows', return_value=True)
        self.hardness.start()
        self.addCleanup(self.hardness.stop)

    def run_candidate(self, old, hits):
        with mock.patch.object(core, '_minimap2_tandem_hits', return_value=hits):
            return core._realign_tandem_stage2(self.stage1, old, self.backbone,
                                              self.qry, self.reader)

    def test_strictly_better_alignment_is_accepted(self):
        old = stage('100=60D60I' + self.qry[100:160] + '40=')
        hit = self.hit(self.ref, self.qry, 0, 200, 0, 200, '200M')
        replacement = self.run_candidate(old, [hit])
        self.assertIsNot(replacement, old)
        self.assertEqual(replacement.pairwise_cigar, '200=')
        self.assertEqual((replacement.pieces[0].ref_path_start,
                          replacement.pieces[0].ref_path_end), (0, 200))

    def test_tie_and_short_local_hit_retain_original_object(self):
        old = stage('200=')
        for hit in (self.hit(self.ref, self.qry, 0, 200, 0, 200, '200M'),
                    self.hit(self.ref, self.qry, 80, 120, 80, 120, '40M')):
            self.assertIs(self.run_candidate(old, [hit]), old)

    def test_equal_score_repeat_shift_does_not_move_an_existing_insertion(self):
        self.qry = self.ref + self.ref[:20]
        old = stage('20I' + self.ref[:20] + '200=')
        hit = self.hit(self.ref, self.qry, 0, 220, 0, 200, '200M20I')
        self.assertIs(self.run_candidate(old, [hit]), old)

    def test_comparison_recomputes_mismatches_in_old_equals(self):
        self.qry = self.ref[1:] + self.ref[:1]
        old = stage('60D60I' + self.qry[:60] + '140=')
        hit = self.hit(self.ref, self.qry, 0, 199, 1, 200, '199M')
        new = self.run_candidate(old, [hit])
        self.assertIsNot(new, old)
        ops = [(t.n, t.op, t.payload) for t in core.parse_cigar_ops(new.pairwise_cigar)]
        self.assert_spells(ops, self.ref, self.qry)

    def test_higher_score_accepts_equal_or_more_large_indels(self):
        self.ref = self.qry = 'A' * 1000
        self.reader = SimpleNamespace(index={'chr': (1200,)},
                                      fetch=mock.Mock(return_value=self.ref))
        self.backbone = core.Coord('chr', 100, 1100, '+')
        old = stage('100=800D800I' + 'A' * 800 + '100=')
        original = core._tm_ops_from_cigar_on_sequences(
            self.ref, self.qry, old.pairwise_cigar)
        for cigar in ('400M60D60I540M', '100M60D60I100M60D60I680M'):
            with self.subTest(candidate=cigar):
                hit = self.hit(self.ref, self.qry, 0, 1000, 0, 1000, cigar)
                self.assertGreater(core.whole_pair_score(hit.ops),
                                   core.whole_pair_score(original))
                self.assertIsNot(self.run_candidate(old, [hit]), old)

    def test_fewer_large_indels_still_requires_higher_score(self):
        self.ref, self.qry = 'A' * 500 + 'C' * 500, 'C' * 500 + 'A' * 500
        self.reader.fetch.return_value = self.ref
        self.reader.index = {'chr': (1200,)}
        self.backbone = core.Coord('chr', 100, 1100, '+')
        old = stage('500D500=500I' + 'A' * 500)
        hit = self.hit(self.ref, self.qry, 0, 1000, 0, 1000, '1000M')
        original = core._tm_ops_from_cigar_on_sequences(
            self.ref, self.qry, old.pairwise_cigar)
        self.assertLess(core.whole_pair_score(hit.ops),
                        core.whole_pair_score(original))
        self.assertIs(self.run_candidate(old, [hit]), old)

    def test_better_alignment_is_accepted_on_both_sides_of_sv_cutoff(self):
        hit = self.hit(self.ref, self.qry, 0, 200, 0, 200, '200M')
        for size in (49, 50, 51):
            with self.subTest(size=size):
                old = stage(f'{size}D{size}I' + self.qry[:size]
                            + f'{200-size}=')
                result = self.run_candidate(old, [hit])
                self.assertIsNot(result, old)

    def test_adjacent_gap_tokens_count_as_one_run(self):
        self.assertEqual(core._tm_count_large_indels([
            (25, 'I', 'A' * 25), (26, 'I', 'A' * 26),
            (60, 'D', ''), (40, 'D', ''), (1, '=', ''),
            (50, 'I', 'A' * 50), (1, '=', ''), (51, 'D', ''),
        ]), 3)

    def test_disabled_single_traversal_and_bad_spans_never_invoke_mapper(self):
        old = stage('200=')
        single = core.build_stage1('>node:200M', '>node:80M60D60M', 'r', 'q')
        with mock.patch.object(core, '_minimap2_tandem_hits') as mapper:
            with mock.patch.object(core, '_SV_REALIGNMENT_ENABLED', False):
                self.assertIs(core._realign_tandem_stage2(
                    self.stage1, old, self.backbone, self.qry, self.reader), old)
            for s1, coord, qry in ((single, self.backbone, self.qry),
                                   (self.stage1, core.Coord('chr', 400, 600, '+'), self.qry),
                                   (self.stage1, self.backbone, self.qry[:-1])):
                self.assertIs(core._realign_tandem_stage2(s1, old, coord, qry, self.reader), old)
            mapper.assert_not_called()

    def test_failed_mapper_and_empty_result_retain_original(self):
        old = stage('200=')
        self.assertIs(self.run_candidate(old, []), old)
        for error in (OSError('missing'), ValueError('invalid hit'),
                      subprocess.TimeoutExpired('minimap2', 1)):
            with mock.patch.object(core, '_minimap2_tandem_hits', side_effect=error):
                self.assertIs(core._realign_tandem_stage2(
                    self.stage1, old, self.backbone, self.qry, self.reader), old)

    def test_minus_backbone_uses_oriented_reference(self):
        self.backbone = core.Coord('chr', 100, 300, '-')
        old = stage('100D100I' + self.qry[:100] + '100=')
        new = self.run_candidate(old, [self.hit(self.ref, self.qry, 0, 200, 0, 200, '200M')])
        self.reader.fetch.assert_called_once_with('chr', 100, 300, '-')
        projected = core.project_linear_piece_on_reference_backbone(new.pieces[0], self.backbone, 0)
        self.assertEqual((projected.genome_start, projected.genome_end, projected.strand),
                         (100, 300, '-'))

    def test_trigger_is_over_50_on_either_graph_side(self):
        for overlap in (50, 51):
            repeated = f'>node:150M50H>node:{150-overlap}H{50+overlap}M'
            for ref, qry in (('>node:200M', repeated), (repeated, '>node:200M')):
                s1 = core.build_stage1(ref, qry, 'r', 'q')
                self.assertEqual(core._has_tandem_realign_candidate(s1), overlap > 50)
        # Opposite orientations are not the tandem event targeted by this pass.
        s1 = core.build_stage1('>node:200M', '>node:150M50H<node:150M50H', 'r', 'q')
        self.assertFalse(core._has_tandem_realign_candidate(s1))


class MinimapBackendTests(SequenceAssertions):
    def test_cli_parses_all_valid_forward_hits_and_rejects_invalid_hits(self):
        ref = qry = 'ACGT' * 50
        def paf(q0, q1, r0, r1, cigar, strand='+'):
            return f'query\t200\t{q0}\t{q1}\t{strand}\treference\t200\t{r0}\t{r1}\t80\t80\t60\tcg:Z:{cigar}'
        output = '\n'.join((paf(0, 80, 0, 80, '80M'), paf(100, 200, 100, 200, '100='),
                            paf(0, 80, 0, 80, '80M', '-'), paf(0, 80, 0, 80, '79M')))
        with mock.patch.object(core, '_mappy', None), \
             mock.patch.object(core, '_minimap2_executable', return_value='minimap2'), \
             mock.patch.object(core.subprocess, 'run', return_value=SimpleNamespace(
                 returncode=0, stdout=output)) as run:
            hits = core._minimap2_tandem_hits(ref, qry)
        self.assertEqual([(h.query_start, h.query_end) for h in hits], [(0, 80), (100, 200)])
        self.assertNotIn('-x', run.call_args.args[0])
        self.assertNotIn('asm5', run.call_args.args[0])
        self.assertIn('--secondary=yes', run.call_args.args[0])
        self.assertIn('--for-only', run.call_args.args[0])

    def test_real_minimap2_tandem_expansion_and_contraction(self):
        backends = []
        if core._mappy is not None:
            backends.append('0')
        if core._minimap2_executable():
            backends.append('1')
        if not backends:
            self.skipTest('minimap2/mappy not installed')
        rng = random.Random(20260920)
        dna = lambda n: ''.join(rng.choices('ACGT', k=n))
        left, unit, right = dna(3000), dna(600), dna(3000)
        shorter, longer = left + unit * 2 + right, left + unit * 4 + right
        for backend in backends:
            for ref, qry, gap in ((shorter, longer, 'I'), (longer, shorter, 'D')):
                with self.subTest(cli=backend, gap=gap), \
                     mock.patch.dict(os.environ, {'GRAPH_CIGARTOREF_NO_MAPPY': backend}):
                    hits = core._minimap2_tandem_hits(ref, qry)
                    self.assertTrue(hits)
                    ops = core._assemble_tandem_hit_intervals(ref, qry, hits)
                    self.assert_spells(ops, ref, qry)
                    self.assertEqual(sum(n for n, op, _ in ops if op == gap), 1200)
                    self.assertEqual(core.whole_pair_score(ops),
                                     4 * min(len(ref), len(qry)) - (24 + 8 * 1200))


class ProjectionIntegrationTests(SequenceAssertions):
    def test_replacement_survives_full_projection_on_both_reference_strands(self):
        rng = random.Random(20)
        dna = lambda n: ''.join(rng.choices('ACGT', k=n))
        left, unit, right = dna(2000), dna(200), dna(2000)
        ref, qry = left + unit * 5 + right, left + unit * 11 + right
        # The graph traversal puts 200 reference-flank bases inside the repeat.
        gcigars = {'ref': '>node:5000M', 'query': '>node:3200M1800H>node:2000H3000M'}
        hit = self.hit(ref, qry, 0, 6200, 0, 5000, '3000M1200I2000M')
        # Supply a valid lower-scoring Stage2 alignment. Projection is exercised
        # with real sequences, coordinates, encoding and serialization.
        original = stage('3000=200D1400I' + qry[3000:4400] + '1800=')
        for strand in ('+', '-'):
            with self.subTest(strand=strand), tempfile.TemporaryDirectory() as tmp:
                chromosome = 'N' * 100 + (ref if strand == '+' else core.revcomp(ref)) + 'N' * 100
                path = Path(tmp) / 'reference.fa'
                path.write_text('>chr\n' + chromosome + '\n')
                Path(str(path) + '.fai').write_text('chr\t5200\t5\t5200\t5201\n')
                reader = core.InMemoryFastaRegionReader(str(path))
                try:
                    backbone = core.Coord('chr', 100, 5100, strand)
                    query_coord = core.Coord('ctg', 0, 6200, strand)
                    pair = core.PairRow(1, '', 'query', 'ref', f'ctg:0-6200{strand}',
                                       f'chr:100-5100{strand}', query_name='query', ref_name='ref',
                                       query_coord=query_coord, ref_coord=backbone)
                    with mock.patch.object(core, '_SV_REALIGNMENT_ENABLED', True), \
                         mock.patch.object(core, 'build_stage2', return_value=original), \
                         mock.patch.object(core, '_minimap2_tandem_hits', return_value=[hit]) as mapper, \
                         mock.patch.object(core, '_realign_tandem_stage2', wraps=core._realign_tandem_stage2) as gate:
                        row = core.build_provisional_row(
                            pair, gcigars, {'node': ref}, {'query': qry, 'ref': ref}, {},
                            {'chr': 5200}, {}, reader, {'query': query_coord, 'ref': backbone})
                    mapper.assert_called_once_with(ref, qry)
                    self.assertEqual(gate.call_count, 1)
                    old = gate.call_args.args[1]
                    self.assertLess(core.whole_pair_score(core._tm_ops_from_cigar_on_sequences(
                        ref, qry, old.pairwise_cigar)),
                        core.whole_pair_score(hit.ops))
                    self.assertEqual(core._tm_count_large_indels(hit.ops), 1)
                    self.assertEqual(core._tm_count_large_indels([
                        (t.n, t.op, t.payload) for t in core.parse_cigar_ops(old.pairwise_cigar)
                    ]), 2)
                    segments = core.parse_reference_cigar_segments(row.split('\t')[-1], 'q')
                    spelled = []
                    for segment in segments:
                        self.assertLessEqual(segment.end, 5200)
                        reference = reader.fetch(segment.path, segment.start, segment.end,
                                                 '+' if segment.direction == '>' else '-')
                        r = 0
                        for op in segment.ops:
                            if op.op in '=M':
                                spelled.append(reference[r:r + op.n])
                            elif op.op == 'X':
                                spelled.append(op.payload[-op.n:])
                            elif op.op == 'I':
                                spelled.append(op.payload)
                            if op.op in '=MXD':
                                r += op.n
                    self.assertEqual(''.join(spelled).upper(), qry.upper())
                    self.assertEqual([(s.start, s.end) for s in segments if s.is_main], [(100, 5100)])
                    self.assertEqual(sum(op.n for s in segments for op in s.ops if op.op == 'I'), 1200)
                finally:
                    reader.close()


if __name__ == '__main__':
    unittest.main()
