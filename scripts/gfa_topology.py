"""Array-based GFA topology: segments, links, and paths for many variants.

This replaces per-cut Python sets and per-segment tuples with integer arrays.
It produces exactly the same GFA as the original set-based builder:

* Every graph source (FASTA root or query-backed allele) gets an integer ID
  in the final segment order: (rank, kind, name) in rGFA mode, (kind, name)
  otherwise. A cut point is one int64 key, ``source_id << POINT_BITS | point``,
  so sorting the keys orders segments exactly as the old numbering did.
* Segment ``i`` is the pair of consecutive keys ``(K[i], K[i+1])`` of one
  source. A leaf spanning ``[start, end)`` covers pairs ``idx(start)`` up to
  ``idx(end)``; only pairs covered by a link leaf become segments.
* Links are integer pairs (``id * 2 + orientation``), deduplicated with
  their minimum source rank by sorting, instead of a dict of tuples.

Leaves are still computed by the resolver, one path at a time, but in forked
worker processes that exit after each batch (so pages they copy are freed).
"""
from array import array
import gc
import math
import multiprocessing as mp

import numpy as np

POINT_BITS = 34
POINT_MASK = (1 << POINT_BITS) - 1
_ORIENT = {'+': 0, '-': 1}
_STATE = None


def source_order(roots, events, resolver, header_lengths, stable):
    """Return (ids by name, names by ID, query flags) in segment order."""
    if stable:
        ranks = resolver.ranks
        ordered = sorted((name for name in roots), key=lambda name: (ranks[name], name))
        ordered.extend(resolver.order)
        placed = set(ordered)
        ordered.extend(name for name in events if name not in placed)
        placed.update(events)
        ordered.extend(sorted(name for name in header_lengths if name not in placed))
    else:
        ordered = sorted(events)
        others = set(roots) | set(header_lengths)
        ordered.extend(sorted(name for name in others if name not in events))
    if len(ordered) >= 1 << (63 - POINT_BITS):
        raise ValueError(f'too many graph sources ({len(ordered)}) for 64-bit segment keys')
    ids = {name: number for number, name in enumerate(ordered)}
    is_query = np.fromiter((name in events for name in ordered), dtype=bool, count=len(ordered))
    return ids, ordered, is_query


def _leaf_batch(bounds):
    """Leaves of specs[lo:hi] as flat arrays (link leaves, then path leaves)."""
    low, high = bounds
    specs, roots, resolver, anchor, ids, path_leaves, link_leaves, ranks = _STATE
    source, start, end, orient = array('q'), array('q'), array('q'), array('b')
    link_counts, path_counts, spec_ranks = array('q'), array('q'), array('q')

    def add(leaves, expected_kind=None):
        for kind, name, leaf_start, leaf_end, orientation in leaves:
            number = ids.get(name)
            if number is None:
                raise ValueError(f'graph source {name!r} has no coordinate definition')
            source.append(number)
            start.append(leaf_start)
            end.append(leaf_end)
            orient.append(_ORIENT[orientation])
        return len(leaves)

    for spec in specs[low:high]:
        link = link_leaves(spec, roots, resolver, anchor)
        path = path_leaves(spec, roots, resolver, anchor)
        link_counts.append(add(link))
        path_counts.append(-1 if path == link else add(path))
        spec_ranks.append(ranks[spec.source] if ranks is not None else 0)
    as_array = lambda values, dtype: np.frombuffer(values, dtype=dtype).copy() if values else np.zeros(0, dtype)
    return (as_array(source, np.int64), as_array(start, np.int64), as_array(end, np.int64),
            as_array(orient, np.int8), as_array(link_counts, np.int64),
            as_array(path_counts, np.int64), as_array(spec_ranks, np.int64))


def collect_leaves(specs, roots, resolver, anchor, ids, path_leaves, link_leaves,
                   processes, log):
    """Compute every spec's leaves, in parallel when fork is available."""
    global _STATE
    ranks = getattr(resolver, 'ranks', None)
    _STATE = specs, roots, resolver, anchor, ids, path_leaves, link_leaves, ranks
    count = len(specs)
    batch = max(256, math.ceil(count / max(1, processes * 16)))
    bounds = [(low, min(count, low + batch)) for low in range(0, count, batch)]
    try:
        if processes > 1 and len(bounds) > 1 and 'fork' in mp.get_all_start_methods():
            # Freeze the shared structures so cyclic GC does not unshare
            # their pages; each worker exits after one batch, returning any
            # pages it copied.
            gc.freeze()
            try:
                with mp.get_context('fork').Pool(processes, maxtasksperchild=1) as pool:
                    parts = []
                    for done, part in enumerate(pool.imap(_leaf_batch, bounds), 1):
                        parts.append(part)
                        if done % max(1, len(bounds) // 8) == 0:
                            log(f'Topology: leaves for {bounds[done - 1][1]}/{count} paths')
            finally:
                gc.unfreeze()
        else:
            parts = [_leaf_batch(item) for item in bounds]
    finally:
        _STATE = None
    if not parts:
        empty = np.zeros(0, np.int64)
        return Leaves(empty, empty, empty, np.zeros(0, np.int8), empty, empty, empty)
    return Leaves(*(np.concatenate([part[index] for part in parts]) for index in range(7)))


class Leaves:
    """Flat leaf arrays; per spec: link leaves then (unless shared) path leaves."""

    def __init__(self, source, start, end, orient, link_counts, path_counts, ranks):
        self.source, self.start, self.end, self.orient = source, start, end, orient
        self.link_counts, self.path_counts, self.ranks = link_counts, path_counts, ranks
        own = link_counts + np.maximum(path_counts, 0)
        self.offsets = np.zeros(len(own) + 1, np.int64)
        np.cumsum(own, out=self.offsets[1:])
        # Mask of rows that are link leaves (path-only rows are excluded).
        spec_of_row = np.repeat(np.arange(len(own), dtype=np.int64), own)
        within = np.arange(len(spec_of_row), dtype=np.int64) - self.offsets[spec_of_row]
        self.link_row = within < link_counts[spec_of_row]
        self.spec_of_row = spec_of_row

    def path_rows(self, spec):
        low = self.offsets[spec]
        if self.path_counts[spec] < 0:
            return low, low + self.link_counts[spec]
        return low + self.link_counts[spec], self.offsets[spec + 1]


def _keys(source, point):
    if len(point) and int(point.max()) > POINT_MASK:
        raise ValueError('a graph coordinate exceeds the 64-bit segment key range')
    return (source << POINT_BITS) | point


def boundaries(leaves, stable, max_node_length):
    """Sorted unique cut keys from all leaves (plus node chopping in rGFA)."""
    parts = [_keys(leaves.source, leaves.start), _keys(leaves.source, leaves.end)]
    if stable:
        counts = np.maximum(0, (leaves.end - leaves.start - 1) // max_node_length)
        total = int(counts.sum())
        if total:
            row = np.repeat(np.arange(len(counts), dtype=np.int64), counts)
            first = np.zeros(len(counts), np.int64)
            np.cumsum(counts[:-1], out=first[1:])
            step = np.arange(total, dtype=np.int64) - first[row] + 1
            parts.append(_keys(leaves.source[row], leaves.start[row] + step * max_node_length))
    return np.unique(np.concatenate(parts))


def pair_ranges(keys, source, start, end):
    """Segment-pair index range [first, last) of each leaf."""
    first = np.searchsorted(keys, _keys(source, start))
    last = np.searchsorted(keys, _keys(source, end))
    return first, last


def number_segments(keys, leaves, first, last, stable, forbidden):
    """Mark pairs covered by link leaves and number them in key order."""
    rows = leaves.link_row & (last > first)
    marks = (np.bincount(first[rows], minlength=len(keys) + 1)
             - np.bincount(last[rows], minlength=len(keys) + 1))
    used = np.cumsum(marks)[:len(keys)] > 0
    numbers = np.cumsum(used, dtype=np.int64)
    if stable:
        # Skip numbers that equal an existing path name, as before.
        for value in sorted(forbidden):
            numbers[numbers >= value] += 1
    return used, numbers


def check_path_coverage(used, first, last, leaves, specs):
    """Every path leaf must lie on link-leaf segments (as the dict lookup required)."""
    unused = np.concatenate(([0], np.cumsum(~used, dtype=np.int64)))
    bad = (last > first) & ~leaves.link_row & (unused[last] - unused[first] > 0)
    if bad.any():
        spec = specs[int(leaves.spec_of_row[np.flatnonzero(bad)[0]])]
        raise ValueError(f'{spec.name}: path segment is not on the graph')


class Segments:
    """Iterable view with the old ``segment_ids.items()`` shape."""

    def __init__(self, keys, used, numbers, names, is_query):
        self.keys, self.used, self.numbers = keys, used, numbers
        self.names, self.is_query = names, is_query
        self.index = np.flatnonzero(used)
        self.count = len(self.index)

    def __len__(self):
        return self.count

    def part(self, low, high):
        """The same view restricted to segments [low, high) in output order."""
        view = object.__new__(Segments)
        view.__dict__.update(self.__dict__)
        view.index = self.index[low:high]
        view.count = len(view.index)
        return view

    def items(self):
        keys, names, is_query, numbers = self.keys, self.names, self.is_query, self.numbers
        for index in self.index.tolist():
            source = int(keys[index] >> POINT_BITS)
            yield (('query' if is_query[source] else 'root', names[source],
                    int(keys[index] & POINT_MASK), int(keys[index + 1] & POINT_MASK)),
                   int(numbers[index]))


def _reduce_links(parts):
    """Unique (left, right) links with their minimum rank, sorted."""
    left = np.concatenate([part[0] for part in parts])
    right = np.concatenate([part[1] for part in parts])
    rank = np.concatenate([part[2] for part in parts])
    order = np.lexsort((rank, right, left))
    left, right, rank = left[order], right[order], rank[order]
    keep = np.ones(len(left), bool)
    keep[1:] = (left[1:] != left[:-1]) | (right[1:] != right[:-1])
    return left[keep], right[keep], rank[keep]


def links(leaves, first, last, numbers, stable, chunk_steps=2_000_000):
    """All links of the link-leaf walks: (left, right, rank), sorted and unique.

    Steps inside one leaf always join pair ``i`` and ``i + 1`` of one source,
    so they are recorded per pair index with np.minimum.at (in chunks of about
    ``chunk_steps``) instead of being materialized and sorted. Only the
    junctions between consecutive leaves of a walk are sorted.
    """
    rows = np.flatnonzero(leaves.link_row & (last > first))
    a, b = first[rows], last[rows]
    o = leaves.orient[rows].astype(np.int64)
    spec = leaves.spec_of_row[rows]
    rank = leaves.ranks[spec]
    unset = np.iinfo(np.int64).max
    # rGFA links are canonical, so both walk directions give the same edge.
    directions = (0,) if stable else (0, 1)
    internal = {direction: np.full(len(numbers), unset, np.int64) for direction in directions}
    steps = b - a - 1
    ends = np.cumsum(steps)
    total = int(ends[-1]) if len(ends) else 0
    cuts = np.searchsorted(ends, np.arange(chunk_steps, total, chunk_steps), side='right')
    for low, high in zip(np.concatenate(([0], cuts)), np.concatenate((cuts, [len(rows)]))):
        count = steps[low:high]
        size = int(count.sum())
        if not size:
            continue
        row = np.repeat(np.arange(low, high, dtype=np.int64), count)
        begin = np.zeros(high - low, np.int64)
        np.cumsum(count[:-1], out=begin[1:])
        pair = a[row] + np.arange(size, dtype=np.int64) - begin[row - low]
        if stable:
            np.minimum.at(internal[0], pair, rank[row])
        else:
            for direction in directions:
                chosen = o[row] == direction
                np.minimum.at(internal[direction], pair[chosen], rank[row][chosen])
    parts = []
    for direction, values in internal.items():
        pair = np.flatnonzero(values != unset)
        if direction == 0:
            parts.append((numbers[pair] * 2, numbers[pair + 1] * 2, values[pair]))
        else:
            parts.append((numbers[pair + 1] * 2 + 1, numbers[pair] * 2 + 1, values[pair]))
    # Junctions: last token of a leaf -> first token of the next leaf.
    same = spec[1:] == spec[:-1]
    tail = np.where(o == 0, b - 1, a)
    head = np.where(o == 0, a, b - 1)
    left = numbers[tail[:-1][same]] * 2 + o[:-1][same]
    right = numbers[head[1:][same]] * 2 + o[1:][same]
    if stable:
        # Canonical orientation: the smaller of the edge and its reverse.
        r_left, r_right = right ^ 1, left ^ 1
        swap = (r_left < left) | ((r_left == left) & (r_right < right))
        left, right = np.where(swap, r_left, left), np.where(swap, r_right, right)
    parts.append((left, right, rank[1:][same]))
    return _reduce_links(parts)


def path_text(spec_index, leaves, first, last, numbers, prefix):
    """The comma-separated oriented segment list of one P line."""
    low, high = leaves.path_rows(spec_index)
    pieces = []
    for row in range(low, high):
        a, b = int(first[row]), int(last[row])
        if b <= a:
            continue
        orientation = '+' if leaves.orient[row] == 0 else '-'
        values = numbers[a:b] if orientation == '+' else numbers[a:b][::-1]
        separator = f'{orientation},{prefix}'
        pieces.append(prefix + separator.join(map(str, values.tolist())) + orientation)
    return ','.join(pieces)
