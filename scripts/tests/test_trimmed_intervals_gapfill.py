"""Trimming updates effective intervals and exposes gaps to the next stage."""
from types import SimpleNamespace

import pytest

import graphcigartoref as core
import graphcigartoref_persample as per
import fill_graphcigartoref_gaps as gaps


class Reader:
    def __init__(self, name, sequence):
        self.sequence = sequence
        self.index = {name: (len(sequence), 0, 1, 1)}

    def fetch(self, name, start, end, strand='+'):
        text = self.sequence[start:end]
        return core.revcomp(text) if strand == '-' else text


def write_row(name, query, reference, cigar, sequence=None, writer='direct'):
    qcoord, rcoord = core.parse_coord(query), core.parse_coord(reference)
    pair = core.PairRow(1, '', name, 'ref', query, reference, '.',
                        query_name=name, ref_name='ref', query_coord=qcoord,
                        ref_coord=rcoord)
    raw = '\t'.join((name, query, 'ref', reference, name, '.', cigar))
    if writer == 'part':
        pair.output_query_name = name
        pair.output_slice = (0, qcoord.end-qcoord.start)
        pair.output_query_coord = qcoord
    if writer == 'merged':
        return per.format_merged_comparison_row(raw, pair, [pair], sequence, reference_flank_rescue=0)
    return per.format_provisional_row(raw, pair, sequence, reference_flank_rescue=0)


@pytest.mark.parametrize('direction', ['>', '<'])
@pytest.mark.parametrize('side', ['leading', 'trailing'])
@pytest.mark.parametrize('writer', ['direct', 'part', 'merged'])
def test_deletion_only_edge_trim_refreshes_reference_interval(direction, side, writer):
    body = '5D200=' if side == 'leading' else '200=5D'
    reference = 'chr:100-305+' if direction == '>' else 'chr:695-900-'
    output = write_row('q', 'ctg:0-200+', reference, f'{direction}chr:100H{body}695H', 'A'*200, writer)
    fields = output.split('\t')
    assert fields[2] == 'ctg:0-200+'
    assert fields[3] == reference  # The source locus remains provenance.
    expected = {('>', 'leading'): 'chr:105-305+', ('>', 'trailing'): 'chr:100-300+',
                ('<', 'leading'): 'chr:695-895-', ('<', 'trailing'): 'chr:700-900-'}
    assert fields[5] == expected[(direction, side)]
    assert not any(op.op=='D' for s in core.parse_reference_cigar_segments(fields[6], 'test') for op in s.ops)


def lift_map(annotations):
    return {a.query_name: SimpleNamespace(locus=a.query_coord, left_extension='0', right_extension='0')
            for a in annotations}


@pytest.mark.parametrize('strand', ['+', '-'])
def test_removed_query_i_creates_gap_that_is_discovered_and_aligned(strand):
    if strand == '+':
        left = '>chr:100H200=5Iaaaaa700H'
        right = '>chr:305H200=495H'
        left_seq = 'A'*205
    else:
        left = '<chr:700H5Ittttt200=100H'
        right = '<chr:495H200=305H'
        left_seq = 'T'*205
    a = gaps.annotation_from_tsv(write_row('q1', f'ctg:0-205{strand}', f'chr:100-300{strand}', left, left_seq), 1)
    b = gaps.annotation_from_tsv(write_row('q2', f'ctg:205-405{strand}', f'chr:305-505{strand}', right), 2)
    assert (a.query_coord.start, a.query_coord.end) == (0, 200)
    query, reference = Reader('ctg', 'A'*405), Reader('chr', 'A'*1000)
    tasks, stats = gaps.discover_gap_tasks([a,b], lift_map([a,b]), query, reference.index, 'sample')
    task, = tasks
    assert (task.start, task.end, task.mode) == (200,205,'bounded')
    assert (task.left_boundary.position, task.right_boundary.position) == (300,305)
    gap_row = gaps.annotation_from_tsv(gaps.build_gap_row(task, query, reference), 3)
    assert gap_row.query_coord == core.Coord('ctg',200,205,'+')
    assert gap_row.graph_cigar == '>chr:300H5=695H'
    remaining, _ = gaps.discover_gap_tasks([a,gap_row,b], lift_map([a,gap_row,b]), query, reference.index, 'sample')
    assert remaining == []


@pytest.mark.parametrize('strand', ['+', '-'])
def test_trimmed_d_exposes_reference_gap_at_exact_query_breakpoint(strand):
    left = '>chr:100H200=5D695H' if strand == '+' else '<chr:695H5D200=100H'
    right = '>chr:305H200=495H' if strand == '+' else '<chr:495H200=305H'
    a = gaps.annotation_from_tsv(write_row('q1', f'ctg:0-200{strand}', f'chr:100-305{strand}', left), 1)
    b = gaps.annotation_from_tsv(write_row('q2', f'ctg:200-400{strand}', f'chr:305-505{strand}', right), 2)
    deletions, _ = gaps.discover_reference_gap_deletions([a,b], lift_map([a,b]), {'chr':(1000,0,1,1)}, 'sample')
    deletion, = deletions
    assert deletion.ref_coord == core.Coord('chr',300,305,'+')
    assert deletion.query_coord.start == deletion.query_coord.end == 200


def test_unresolved_query_gap_is_not_mistaken_for_a_deletion():
    a = gaps.annotation_from_tsv(write_row('q1','ctg:0-205+','chr:100-300+','>chr:100H200=5Iaaaaa700H'),1)
    b = gaps.annotation_from_tsv(write_row('q2','ctg:205-405+','chr:305-505+','>chr:305H200=495H'),2)
    deletions, stats = gaps.discover_reference_gap_deletions([a,b],lift_map([a,b]),{'chr':(1000,0,1,1)},'sample')
    assert deletions == []
    assert stats['positive_query_gaps_skipped'] == 1


def test_contig_end_exposed_by_trimming_is_not_an_internal_gap():
    a = gaps.annotation_from_tsv(write_row('q1','ctg:0-205+','chr:100-300+','>chr:100H200=5Iaaaaa700H'),1)
    tasks, stats = gaps.discover_gap_tasks([a],lift_map([a]),Reader('ctg','A'*205),{'chr':(1000,0,1,1)},'sample')
    assert tasks == []
    assert stats['internal_gaps'] == 0
