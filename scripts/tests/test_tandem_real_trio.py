"""Replay separate SV/D-I scoring through real minimap2 and VCF calling."""
import gzip
import json
import os
from pathlib import Path
import unittest
from unittest import mock

import graphcigartoref as core
import graphcigartoref_persample as per
import graphreftovcf as caller


class MemoryReader:
    def __init__(self, sequences):
        self.sequences = dict(sequences)
        self.index = {
            name: (len(sequence), 0, max(1, len(sequence)), max(1, len(sequence)))
            for name, sequence in self.sequences.items()
        }

    def fetch(self, name, start, end, strand="+"):
        sequence = self.sequences[name][int(start):int(end)]
        return caller.revcomp(sequence) if strand == "-" else sequence


FIXTURE = Path(__file__).with_name('fixtures') / 'na19239_na19240_group81264.json.gz'


class LocusReader:
    def __init__(self, data):
        self.data = data
        self.index = {data['chrom']: (data['chrom_length'],)}

    def fetch(self, name, start, end, strand='+'):
        data = self.data
        if name != data['chrom'] or not data['start'] <= start <= end <= data['end']:
            raise ValueError((name, start, end))
        seq = data['reference'][start - data['start']:end - data['start']]
        return core.revcomp(seq) if strand == '-' else seq


def verify_spelling(test, cigar, sequence, reader):
    spelled = []
    for segment in core.parse_reference_cigar_segments(cigar, 'query'):
        reference = reader.fetch(segment.path, segment.start, segment.end,
                                 '+' if segment.direction == '>' else '-')
        cursor = 0
        for op in segment.ops:
            if op.op in 'M=':
                spelled.append(reference[cursor:cursor + op.n])
            elif op.op in 'XI':
                test.assertEqual(len(op.payload), op.n)
                spelled.append(op.payload)
            if op.op != 'H':
                cursor += core.ref_consume_pair(op.op, op.n)
        test.assertEqual(cursor, len(reference))
    test.assertEqual(''.join(spelled).upper(), sequence.upper())


def call_owned(data, sample, cigar):
    lo, hi = sample['owned_query_slice']
    sequence = sample['sequence'][lo:hi]
    if sample['strand'] == '+':
        start, end = sample['start'] + lo, sample['start'] + hi
    else:
        start, end = sample['end'] - hi, sample['end'] - lo
    source = caller.SeqRecord('query', 'query', sample['contig'], start, end,
                              sample['strand'], len(sequence),
                              reader=MemoryReader({'query': sequence}))
    reference = caller.SeqRecord(data['chrom'], data['chrom'], data['chrom'],
                                 0, data['chrom_length'], '+', data['chrom_length'])
    result = caller.parse_row_events(
        dict(line_no=1, schema='persample_v2', query='query',
             query_coord=f"{sample['contig']}:{start}-{end}{sample['strand']}",
             labels='query', label_coords='.', encoded=cigar),
        ref_lookup={data['chrom']: reference}, query_lookup={'query': source},
        label_lookup={'query': source}, sequence_lookup={'query': source},
        max_impute=10000, bp_uncertain=200, gap_break=100, gap_penalty=4,
        var_in_insert_min_bp=100, next_id_start=1, strict=True,
    )
    return sorted((ev.type_prefix + ev.svtype, ev.ref_start, ev.ref_end,
                   ev.ins_size, ev.del_size) for ev in result.svs
                  if max(ev.ins_size, ev.del_size) >= 50)


class RealTrioTests(unittest.TestCase):
    def test_separate_scores_preserve_sequences_and_record_backend_sv_outcomes(self):
        data = json.loads(gzip.decompress(FIXTURE.read_bytes()))
        reader = LocusReader(data)
        backbone = core.Coord(data['chrom'], data['start'], data['end'], '+')
        backends = []
        if core._mappy is not None:
            backends.append('0')
        if core._minimap2_executable():
            backends.append('1')
        if not backends:
            self.skipTest('minimap2/mappy not installed')
        for backend in backends:
            signatures = []
            for name, sample in data['samples'].items():
                with self.subTest(backend=backend, sample=name), \
                     mock.patch.dict(os.environ, {'GRAPH_CIGARTOREF_NO_MAPPY': backend}), \
                     mock.patch.object(core, '_SV_REALIGNMENT_ENABLED', True):
                    query = sample['sequence']
                    qcoord = core.Coord(sample['contig'], sample['start'],
                                        sample['end'], sample['strand'])
                    pair = core.PairRow(1, '', 'query', 'reference',
                        f'{qcoord.chrom}:{qcoord.start}-{qcoord.end}{qcoord.strand}',
                        f'{backbone.chrom}:{backbone.start}-{backbone.end}+',
                        query_name='query', ref_name='reference',
                        query_coord=qcoord, ref_coord=backbone)
                    # Capture the actual gate result without supplying aligner hits.
                    changes = []
                    mapped_hits = []
                    real_gate = core._realign_tandem_stage2
                    real_mapper = core._minimap2_tandem_hits
                    def capture_hits(*args):
                        hits = real_mapper(*args)
                        mapped_hits.extend(hits)
                        return hits
                    def capture(*args):
                        result = real_gate(*args)
                        changes.append((args[1], result))
                        return result
                    with mock.patch.object(core, '_realign_tandem_stage2', side_effect=capture), \
                         mock.patch.object(core, '_minimap2_tandem_hits',
                                           side_effect=capture_hits) as mapper:
                        row = core.build_provisional_row(
                            pair, {'query': sample['graph_cigar'],
                                   'reference': '>' + data['path'] + ':92353M'},
                            {data['path']: data['reference']},
                            {'query': query, 'reference': data['reference']}, {},
                            {data['chrom']: data['chrom_length']}, {}, reader,
                            {'query': qcoord, 'reference': backbone})
                    self.assertEqual(len(changes), 1)
                    old, new = changes[0]
                    hardness = core.quantify_sequence_hardness(
                        query)['summary']['hardness_score']
                    if hardness is not None and hardness < 10:
                        mapper.assert_called_once_with(data['reference'], query)
                        original = core._tm_ops_from_cigar_on_sequences(
                            data['reference'], query, old.pairwise_cigar)
                        candidate = core._assemble_tandem_hit_intervals(
                            data['reference'], query, mapped_hits)
                        self.assertIsNot(old, new)
                        self.assertGreater(core.whole_pair_score(candidate),
                                           core.whole_pair_score(original))
                        self.assertGreaterEqual(core._tm_count_large_indels(candidate),
                                                core._tm_count_large_indels(original))
                    else:
                        mapper.assert_not_called()
                        self.assertIs(old, new)
                    cigar = row.split('\t')[6]
                    verify_spelling(self, cigar, query, reader)
                    lo, hi = sample['owned_query_slice']
                    owned = per.slice_graph_cigar_by_query(cigar, lo, hi,
                                                           reference_alignment=True)
                    verify_spelling(self, owned, query[lo:hi], reader)
                    signatures.append(call_owned(data, sample, owned))
            self.assertEqual(len(signatures), 2, backend)
            # Downstream processing must preserve father/child SV agreement,
            # including when the stricter early gate rejects both queries.
            self.assertTrue(all(signatures), backend)
            self.assertEqual(signatures[0], signatures[1], backend)




if __name__ == '__main__':
    unittest.main()
