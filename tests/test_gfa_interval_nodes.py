"""Dependency-ordered interval nodes agree with recursive graph expansion."""
from collections import OrderedDict
from pathlib import Path
import random
import sys

import numpy as np
import pytest

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / 'scripts'))
from gfa_interval_nodes import IntervalWalk
from gfa_interval_pipeline import (Event, Hit, PathSpec, Root, Run, _drop_unverified,
                                   _path_leaves, _link_leaves, _path_and_link_leaves)
from gfa_source_catalog import SourceAlias
from gfa_stable_coords import StableResolver
from gfa_junctions import _allele_parent, _involved
from gfa_topology import Leaves, POINT_BITS, boundaries


def graph(seed=34, count=100):
    rng = random.Random(seed)
    roots = {
        'chr': Root('chr', '', None, 10000, True, 'reference', 0),
        'alias': Root('alias', '', None, 100, False, 'alternative', -1),
        'lifted': Root('lifted', '', None, 80, True, 'alternative', 1,
                       lift={'left': ('chr', 600, '-'), 'right': ('chr', 400, '-')}),
        'DUP_chr_10_90': Root('DUP_chr_10_90', '', None, 80, True, 'duplication', 2),
    }
    events, hits = OrderedDict(), {}
    available = list(roots)
    for i in range(count):
        name = f'v{i}'
        parent = rng.choice(available)
        length = events[parent].length if parent in events else roots[parent].length
        size = rng.randrange(1, 50)
        point = rng.randrange(length+1)
        nested = parent in events or parent == 'DUP_chr_10_90'
        pos = point + int(nested)
        runs, cursor = [], 0
        while cursor < size:
            span = rng.randrange(1, size-cursor+1)
            if rng.random() < 0.6:
                target = rng.choice(available)
                target_length = events[target].length if target in events else roots[target].length
                span = min(span, target_length)
                start = rng.randrange(target_length-span+1)
                strand = rng.choice('+-')
                runs.append(Run(cursor, cursor+span, '=', target, start, start+span, strand, 0))
            else:
                runs.append(Run(cursor, cursor+span, 'I', None, 0, 0, '+', 1))
            cursor += span
        events[name] = Event(i, name, parent, pos, pos, size, runs, None, '', True)
        left = rng.randrange(0, 5)
        hits[name] = Hit('', 0, size+left+3, left, 3)
        available.append(name)
    resolver = StableResolver(events, roots, {}, hits, set(events),
                              source_aliases={'alias': SourceAlias('chr', 1000, 1100, '-')})
    return resolver


@pytest.mark.parametrize('seed', [1, 34, 81, 217])
def test_nested_source_and_flank_intervals_match_recursive_resolver(seed):
    resolver = graph(seed)
    rng = random.Random(seed)
    calls, expected = [], []
    for name in [*resolver.roots, *resolver.events]:
        size = resolver.target_length(name)
        intervals = [(0, size), (0, 0), (size, size)]
        intervals.extend(sorted((rng.randrange(size+1), rng.randrange(size+1))) for _ in range(12))
        for start, end in intervals:
            for strand in '+-':
                calls.append(('expand', (name, start, end, strand)))
                expected.append(resolver.expand(name, start, end, strand))
        for point in (0, size//2, size):
            for flank in (0, 1, 10, 50, 75):
                for side in ('left', 'right'):
                    calls.append(('_flank', (name, point, flank, side)))
                    expected.append(resolver._flank(name, point, flank, side))
    resolver.prepare_sources(50, lambda message: None)
    for (method, args), leaves in zip(calls, expected):
        assert getattr(resolver, method)(*args) == leaves, (method, args)


def test_walk_slices_preserve_boundaries_and_reverse_coordinates():
    walk = IntervalWalk([('root', 'a', 10, 20, '+'), ('query', 'b', 30, 35, '-'),
                         ('root', 'a', 20, 30, '+')])
    assert walk.slice(8, 18) == [('root', 'a', 18, 20, '+'), ('query', 'b', 30, 35, '-'),
                                 ('root', 'a', 20, 23, '+')]
    assert walk.slice(11, 14, '-') == [('query', 'b', 31, 34, '+')]
    assert walk.slice(10, 15) == [('query', 'b', 30, 35, '-')]
    assert walk.slice(10, 10) == []


def test_deep_copies_share_walk_and_parent_flanks_without_recursion():
    count = 1500
    roots = {'chr': Root('chr', '', None, 10000, True, 'reference', 0)}
    events = OrderedDict()
    for i in range(count):
        name, parent = f'v{i}', f'v{i-1}' if i else 'chr'
        run = Run(0, 100, '=', parent, 0, 100, '+', 0) if i else Run(0, 100, 'I', None, 0, 0, '+', 0)
        pos = 1 if i else 500
        events[name] = Event(i, name, parent, pos, pos, 100, [run], None, '', True)
    resolver = StableResolver(events, roots, {}, {'v0': Hit('', 0, 100, 0, 0)}, set(events))
    resolver.prepare_sources(50, lambda message: None)
    assert resolver.source_intervals.maximum_depth == count-1
    assert resolver.source_intervals.walks['v1'] is resolver.source_intervals.walks[f'v{count-1}']
    assert resolver.expand(f'v{count-1}', 5, 9) == [('query', 'v0', 5, 9, '+')]
    assert resolver._flank(f'v{count-1}', 0, 50, 'left') == [('root', 'chr', 450, 500, '+')]


def test_variant_coordinate_index_and_deletion_parent_flanks():
    roots = {'chr': Root('chr', '', None, 1000, True, 'reference', 0)}
    events = OrderedDict(
        snp=Event(0, 'snp', 'chr', 11, 11, 1, [Run(0, 1, 'I', None, 0, 0, '+', 0)], None, '', True, 'snp'),
        deletion=Event(1, 'deletion', 'chr', 100, 140, 0, [Run(0, 0, 'D', None, 0, 30, '+', 0)], None, '', True, 'deletion'),
        insertion=Event(2, 'insertion', 'deletion', 5, 5, 2, [Run(0, 2, 'I', None, 0, 0, '+', 0)], None, '', True),
    )
    hits = {'snp': Hit('', 0, 1, 0, 0), 'insertion': Hit('', 0, 2, 0, 0)}
    resolver = StableResolver(events, roots, {}, hits, set(events))
    expected = resolver.local_path('insertion', 0, 2, 50)
    resolver.prepare_sources(50, lambda message: None)
    assert resolver.local_path('insertion', 0, 2, 50) == expected
    assert resolver.index_variants().by_parent == {
        'chr': [(10, 11, 'snp', 'snp'), (100, 130, 'deletion', 'deletion')],
        'deletion': [(0, 0, 'insertion', 'insertion')]}


def test_sorted_junction_index_matches_pairwise_boundaries():
    resolver = graph(count=200)
    events = resolver.events
    spans = {name: resolver.breakpoints(name) for name in resolver.order}
    repeated = set(list(events)[::7])
    expected = set()
    for a in events:
        for b in events:
            if a != b and events[a].chrom == events[b].chrom and spans[a][1] == spans[b][0]:
                expected.update((a, b))
        if a in repeated and spans[a][0] == spans[a][1]:
            expected.add(a)
        parent = _allele_parent(events[a], resolver.roots)
        if parent != a and any(event.chrom == parent for event in events.values()):
            expected.add(a)
    pending = list(expected)
    while pending:
        parent = _allele_parent(events[pending.pop()], resolver.roots)
        for name, event in events.items():
            if event.chrom == parent and name not in expected:
                expected.add(name)
                pending.append(name)
    actual, actual_spans = _involved(events, resolver, repeated)
    assert actual == expected
    assert actual_spans == spans


@pytest.mark.parametrize('anchor', [0, 1, 50])
def test_paired_paths_reuse_core_without_changing_links(anchor):
    resolver = graph(count=30)
    specs = [PathSpec(name, root.kind, name, 0, root.length)
             for name, root in resolver.roots.items()]
    for name, event in resolver.events.items():
        specs.extend([PathSpec(name, event.kind, name, 0, event.length),
                      PathSpec(name+'#unique', 'unique', name, 0, event.length//2)])
    resolver.junctions = [('chr', 0, 2, (('v0', ()), ('v1', ())))]
    specs.append(PathSpec('junction', 'junction', 'v0', 0, 0))
    expected = [(_path_leaves(spec, resolver.roots, resolver, anchor),
                 _link_leaves(spec, resolver.roots, resolver, anchor)) for spec in specs]
    resolver.prepare_sources(anchor, lambda message: None)
    assert [_path_and_link_leaves(spec, resolver.roots, resolver, anchor)
            for spec in specs] == expected


def test_dropped_variants_invalidate_coordinate_and_source_indexes(tmp_path):
    resolver = graph(count=30)
    resolver.prepare_sources(50, lambda message: None)
    _drop_unverified(resolver.events, resolver, {}, [('v0', 'sequence_mismatch')],
                     tmp_path/'unresolved.tsv', {}, lambda message: None)
    assert resolver.variant_intervals is None
    assert resolver.source_intervals is None
    resolver.prepare_sources(50, lambda message: None)
    assert 'v0' not in resolver.index_variants().spans
    assert set(resolver.index_variants().spans) == set(resolver.order)
    assert not (set(resolver.source_intervals.walks) - set(resolver.order))


def test_uncovered_source_intervals_are_rejected():
    resolver = graph(count=1)
    event = resolver.events['v0']
    event.runs = [Run(1, event.length, 'I', None, 0, 0, '+', 1)]
    with pytest.raises(ValueError, match='graph definition has a gap at 0'):
        resolver.prepare_sources(50, lambda message: None)


def make_leaves(source, start, end):
    size = len(source)
    return Leaves(np.array(source, np.int64), np.array(start, np.int64), np.array(end, np.int64),
                  np.zeros(size, np.int8), np.ones(size, np.int64),
                  np.full(size, -1, np.int64), np.zeros(size, np.int64))


@pytest.mark.parametrize('node_length', [1, 17, 1024, 10000])
def test_shared_node_cuts_equal_per_path_chopping(node_length):
    rng = random.Random(141)
    source, start, end = [], [], []
    for _ in range(200):
        source.append(rng.randrange(4))
        start.append(rng.choice([0, 100, 1000, rng.randrange(10000)]))
        end.append(start[-1] + rng.randrange(0, 20000))
    leaves = make_leaves(source, start, end)
    expected = set()
    for name, low, high in zip(source, start, end):
        expected.update((name << POINT_BITS) | cut for cut in [low, high])
        expected.update((name << POINT_BITS) | cut for cut in range(low+node_length, high, node_length))
    assert boundaries(leaves, True, node_length).tolist() == sorted(expected)


def test_many_copies_chop_the_source_once(monkeypatch):
    copies, length, node_length = 10000, 10**9, 10**6
    leaves = make_leaves([0]*copies, [0]*copies, [length]*copies)
    original, expanded = np.repeat, []
    def repeat(values, counts):
        expanded.append(int(counts.sum()))
        return original(values, counts)
    monkeypatch.setattr(np, 'repeat', repeat)
    assert boundaries(leaves, True, node_length).tolist() == list(range(0, length+1, node_length))
    assert expanded == [length//node_length - 1]
