"""Default mapping for every graph payload caller and the group82000 repair."""
import gzip
import json
import os
from pathlib import Path
from types import SimpleNamespace
import unittest
from unittest import mock

import graphcigartoref as core
from tests.test_late_di_rescue_score import fixture


FIXTURE = Path(__file__).with_name('fixtures') / 'group82000_di_repair.json.gz'


def early_payload(reference, query):
    result = []
    core._tm_realign_insertion_payloads(reference, query, result)
    return result


PAYLOAD_CALLERS = (
    core.minimap2_payload_ops, early_payload,
    core._tm_realign_secondary_sv_window, core._realign_long_di_sequences,
)


class BackendSelectionTests(unittest.TestCase):
    def test_di_hit_ranking_uses_affine_score_on_both_backends(self):
        reference, query = 'A' * 255, 'A' * 254 + 'C'
        hits = [SimpleNamespace(strand=1, is_primary=True, q_st=0, q_en=qend,
                 r_st=0, r_en=qend, mlen=matches, blen=block, cigar_str=cigar)
                for cigar, qend, matches, block in (
                    ('50M1D1I' * 5, 255, 250, 260), ('248M', 248, 248, 248))]
        paf = '\n'.join(
            f'graph_query\t255\t0\t{h.q_en}\t+\tgraph_ref\t255\t0\t{h.r_en}'
            f'\t{h.mlen}\t{h.blen}\t60\tcg:Z:{h.cigar_str}\ttp:A:P' for h in hits)
        module = SimpleNamespace(Aligner=lambda **kw: SimpleNamespace(map=lambda q: hits))
        for backend in ('0', '1'):
            with self.subTest(backend=backend), \
                 mock.patch.dict(os.environ, {'GRAPH_CIGARTOREF_NO_MAPPY': backend}), \
                 mock.patch.object(core, '_mappy', module), \
                 mock.patch.object(core, '_assemble_tandem_hit_intervals', return_value=None), \
                 mock.patch.object(core, '_minimap2_executable', return_value='minimap2'), \
                 mock.patch.object(core.subprocess, 'run', return_value=SimpleNamespace(
                     returncode=0, stdout=paf, stderr='')):
                sv_candidate = core.minimap2_payload_ops(reference, query)
                di_candidate = core.minimap2_payload_ops(reference, query, score=core.di_repair_score)
            # Isolate single-hit ranking when no combined chain is available.
            # The split-hit regression separately tests successful assembly.
            self.assertNotEqual(sv_candidate, di_candidate)
            self.assertEqual(di_candidate, [(248, '=', ''), (7, 'D', ''), (7, 'I', 'AAAAAAC')])
            self.assertGreater(core.di_repair_score(di_candidate), core.di_repair_score(sv_candidate))
            self.assertLess(core._tm_realignment_score(di_candidate), core._tm_realignment_score(sv_candidate))

    def test_long_di_helper_requests_its_affine_ranking_policy(self):
        with mock.patch.object(core, 'minimap2_payload_ops', return_value=[(1500, '=', '')]) as align:
            core._realign_long_di_sequences('A' * 1500, 'A' * 1499 + 'C')
        self.assertIs(align.call_args.kwargs['score'], core.di_repair_score)
        self.assertIsNone(align.call_args.kwargs['preset'])

    def test_di_cli_uses_default_mapping(self):
        reference, query = 'A' * 1500, 'A' * 1499 + 'C'
        paf = ('graph_query\t1500\t0\t1500\t+\tgraph_ref\t1500\t0\t1500'
               '\t1499\t1500\t60\tcg:Z:1500M\ttp:A:P')
        for caller in PAYLOAD_CALLERS:
            with self.subTest(caller=caller.__name__), \
                 mock.patch.dict(os.environ, {'GRAPH_CIGARTOREF_NO_MAPPY': '1'}), \
                 mock.patch.object(core, '_parasail', None), \
                 mock.patch.object(core, '_minimap2_executable', return_value='minimap2'), \
                 mock.patch.object(core.subprocess, 'run', return_value=SimpleNamespace(
                     returncode=0, stdout=paf, stderr='')) as run:
                candidate = caller(reference, query)
                self.assertIsNotNone(candidate)
                run.assert_called_once()
                command = run.call_args.args[0]
                self.assertNotIn('-x', command)
                self.assertNotIn('asm5', command)
                self.assertIn('--for-only', command)
                self.assertIn('--secondary=no', command)
                self.assertEqual(sum(n for n, op, _ in candidate if op in '=M'), 1499)

    def test_di_mappy_uses_default_mapping(self):
        reference, query = 'A' * 1500, 'A' * 1499 + 'C'
        hit = SimpleNamespace(strand=1, is_primary=True, q_st=0, q_en=1500,
                              r_st=0, r_en=1500, mlen=1499, blen=1500,
                              cigar_str='1500M')
        module = SimpleNamespace(Aligner=mock.Mock(return_value=SimpleNamespace(
            map=mock.Mock(return_value=[hit]))))
        for caller in PAYLOAD_CALLERS:
            module.Aligner.reset_mock()
            with self.subTest(caller=caller.__name__), \
                 mock.patch.dict(os.environ, {'GRAPH_CIGARTOREF_NO_MAPPY': '0'}), \
                 mock.patch.object(core, '_parasail', None), \
                 mock.patch.object(core, '_mappy', module):
                candidate = caller(reference, query)
                self.assertIsNotNone(candidate)
                module.Aligner.assert_called_once_with(
                    seq=reference, preset='map-ont', n_threads=1,
                    extra_flags=core._MINIMAP_FOR_ONLY)


class RealRepairTests(unittest.TestCase):
    def test_supported_repair_improves_affine_score_without_changing_spans(self):
        data = json.loads(gzip.decompress(FIXTURE.read_bytes()))
        reference, query = data['reference'], data['query']
        backends = []
        if core._mappy is not None:
            backends.append('0')
        if core._minimap2_executable():
            backends.append('1')
        if not backends:
            self.skipTest('minimap2/mappy not installed')
        old_cigar = '2926D2927I' + query
        old = core.parse_cigar_ops(old_cigar)
        for backend in backends:
            for strand in ('+', '-'):
                with self.subTest(backend=backend, strand=strand), \
                     mock.patch.dict(os.environ, {'GRAPH_CIGARTOREF_NO_MAPPY': backend}):
                    candidate = core._realign_long_di_sequences(reference, query)
                    self.assertIsNotNone(candidate)
                    # The separate affine D/I policy accepts this repair even
                    # though its micro-indels worsen the logarithmic SV score.
                    self.assertTrue(core._di_rescue_is_valid(old, candidate, reference, query))
                    self.assertGreater(core.di_repair_score(candidate),
                                       core.di_repair_score(old))
                    self.assertLess(core._tm_realignment_score(candidate),
                                    core._tm_realignment_score([
                                        (t.n, t.op, t.payload) for t in old]))
                    self.assertGreaterEqual(sum(n for n, op, _ in candidate if op == '='), 2500)
                    self.assertEqual(core._tm_count_large_indels(candidate), 0)
                    args = fixture('40=' + old_cigar + '40=',
                                   'A' * 40 + reference + 'A' * 40,
                                   'A' * 40 + query + 'A' * 40, strand)
                    piece, linear = core._rescue_long_di_pairs_in_linear_piece(*args)
                    self.assertNotEqual(piece.cigar, args[0].cigar)
                    self.assertEqual((linear.genome_start, linear.genome_end, linear.strand),
                                     (100, 3106, strand))
                    ref = 'A' * 40 + reference + 'A' * 40
                    spelled = []
                    r = 0
                    for tok in core.parse_cigar_ops(piece.cigar):
                        if tok.op in '=M':
                            spelled.append(ref[r:r + tok.n])
                        elif tok.op == 'X':
                            spelled.append(tok.payload[-tok.n:])
                        elif tok.op == 'I':
                            spelled.append(tok.payload)
                        if tok.op in '=MXD':
                            r += tok.n
                    self.assertEqual(r, len(ref))
                    self.assertEqual(''.join(spelled).upper(), ('A' * 40 + query + 'A' * 40).upper())
                    self.assertEqual(core.pairwise_query_span(piece.cigar), 3007)
                    self.assertEqual(core._tm_count_large_indels([
                        (t.n, t.op, t.payload) for t in core.parse_cigar_ops(piece.cigar)]), 0)


if __name__ == '__main__':
    unittest.main()
