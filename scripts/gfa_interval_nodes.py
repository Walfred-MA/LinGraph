"""Resolve variant sequence intervals once, in parent/copy dependency order.

Nodes retain their source interval and strand, not a second copy of the bases.
Keeping run boundaries is intentional: they participate in GFA node chopping.
"""
from array import array
from bisect import bisect_left, bisect_right
from collections import defaultdict
import time


class VariantIntervals:
    """Sorted (start, end, variant ID, type) rows keyed by chromosome/insertion."""

    def __init__(self, resolver):
        self.spans = {}
        self.by_parent = defaultdict(list)
        for name in resolver.order:
            event = resolver.events[name]
            start, end = resolver.breakpoints(name)
            self.spans[name] = start, end
            self.by_parent[event.chrom].append((start, end, name, event.kind))
        for rows in self.by_parent.values():
            # Stable ties retain the established variant dependency order.
            rows.sort(key=lambda row: (row[0], row[1]))


class IntervalWalk:
    """Immutable source intervals with cumulative ends for binary-search slicing."""

    __slots__ = ('leaves', 'ends', 'length')

    def __init__(self, leaves):
        self.leaves = tuple(leaves)
        self.ends = array('q')
        self.length = 0
        for _kind, _name, start, end, _strand in self.leaves:
            self.length += end - start
            self.ends.append(self.length)

    def slice(self, start, end, orientation='+'):
        if start == end:
            return []
        if start == 0 and end == self.length:
            leaves = list(self.leaves)
        else:
            leaves = []
            first = bisect_right(self.ends, start)
            last = bisect_left(self.ends, end) + 1
            for index in range(first, last):
                leaf = self.leaves[index]
                before = self.ends[index-1] if index else 0
                low, high = max(start, before)-before, min(end, self.ends[index])-before
                kind, name, a, b, strand = leaf
                if low == 0 and high == b-a:
                    leaves.append(leaf)
                elif strand == '+':
                    leaves.append((kind, name, a+low, a+high, strand))
                else:
                    leaves.append((kind, name, b-high, b-low, strand))
        if orientation == '-':
            return [(kind, name, a, b, '-' if strand == '+' else '+')
                    for kind, name, a, b, strand in reversed(leaves)]
        return leaves


class SourceIntervals:
    """Materialized interval walks shared by every path and forked topology worker.

    Literal-only variants (including SNPs) need no stored walk: their one source
    interval is already indexed by the verified query hit. Only composite/copy
    paths are flattened. Boundary flanks are stored only for parents of nested
    variants, so flat SNP cohorts do not acquire two extra walks per SNP.
    """

    def __init__(self, resolver, anchor, log):
        self.resolver = resolver
        self.anchor = max(1, anchor)
        self.walks = {}
        self.flanks = {}
        self.maximum_depth = 0
        started = time.monotonic()
        coordinates = resolver.index_variants()
        parents = coordinates.by_parent.keys() & resolver.events.keys()
        depths = {}
        for name in resolver.order:
            event = resolver.events[name]
            dependencies = {event.chrom}
            dependencies.update(run.target for run in event.runs if run.target)
            depth = max((depths.get(target, 0) + 1 for target in dependencies
                         if target in resolver.events), default=0)
            if depth:
                depths[name] = depth
                self.maximum_depth = max(self.maximum_depth, depth)
            if event.length and not self._literal(event):
                self.walks[name] = self._resolve(event)
            if name in parents:
                # Parents are ready before children, including copied targets.
                # Each boundary is expanded once to the maximum requested flank.
                self.flanks[name] = (
                    IntervalWalk(self.flank(name, 0, self.anchor, 'left')),
                    IntervalWalk(self.flank(name, event.length, self.anchor, 'right')))
        log(f'Intervals: resolved {len(resolver.order)} variants through nesting depth '
            f'{self.maximum_depth}; {len(self.walks)} composite walks and '
            f'{len(self.flanks)} parent flank pairs in {time.monotonic()-started:.1f}s')

    @staticmethod
    def _literal(event):
        return (len(event.runs) == 1 and event.runs[0].operation != '='
                and event.runs[0].qstart == 0 and event.runs[0].qend == event.length)

    def _resolve(self, event):
        resolver = self.resolver
        # A complete forward copy shares its source's entire interval table.
        if len(event.runs) == 1:
            run = event.runs[0]
            prior = self.walks.get(run.target)
            if (prior is not None and run.operation == '=' and run.orientation == '+'
                    and run.qstart == run.rstart == 0
                    and run.qend == run.rend == event.length == prior.length):
                return prior
        leaves, covered = [], 0
        for run in event.runs:
            if run.qend <= 0 or run.qstart >= event.length:
                continue
            low, high = max(0, run.qstart), min(event.length, run.qend)
            if low != covered:
                raise ValueError(f'{event.identifier}: graph definition has a gap at {covered}')
            if run.operation == '=':
                a, b = resolver.target_interval(run.target, run.rstart, run.rend, run.orientation)
                if run.orientation == '+':
                    a, b = a+low-run.qstart, a+high-run.qstart
                else:
                    a, b = b-(high-run.qstart), b-(low-run.qstart)
                leaves.extend(self.expand(run.target, a, b, run.orientation))
            else:
                hit = resolver.hits.get(event.identifier)
                if hit is None:
                    raise ValueError(f'{event.identifier}: query sequence is required for {run.operation}')
                leaves.append(('query', event.identifier, hit.left+low, hit.left+high, '+'))
            covered = high
        if covered != event.length:
            raise ValueError(f'{event.identifier}: graph definition ends at {covered}, expected {event.length}')
        return IntervalWalk(leaves)

    def expand(self, name, start, end, orientation='+'):
        resolver = self.resolver
        event = resolver.events.get(name)
        if event is None:
            root = resolver.roots.get(name)
            if root is None:
                raise ValueError(f'graph target {name!r} has no FASTA sequence')
            if not 0 <= start <= end <= root.length:
                raise ValueError(f'{name}: interval {start}-{end} is outside 0-{root.length}')
        if start == end:
            return []
        if event is None:
            name, start, end, orientation = resolver.canonical_interval(name, start, end, orientation)
            return [('root', name, start, end, orientation)]
        if not 0 <= start <= end <= event.length:
            raise ValueError(f'{name}: interval {start}-{end} is outside 0-{event.length}')
        walk = self.walks.get(name)
        if walk is not None:
            return walk.slice(start, end, orientation)
        if not self._literal(event):
            raise ValueError(f'{name}: source intervals requested before its dependencies were resolved')
        hit = resolver.hits.get(name)
        if hit is None:
            raise ValueError(f'{name}: query sequence is required for {event.runs[0].operation}')
        return [('query', name, hit.left+start, hit.left+end, orientation)]

    def flank(self, name, boundary, size, side):
        resolver = self.resolver
        length = resolver.target_length(name)
        if not 0 <= boundary <= length:
            raise ValueError(f'{name}: flank boundary {boundary} is out of bounds')
        if name in resolver.source_aliases and not resolver.roots[name].lift:
            target, point, _end, strand = resolver.canonical_interval(name, boundary, boundary)
            target_side = side if strand == '+' else ('right' if side == 'left' else 'left')
            leaves = self.flank(target, point, size, target_side)
            return leaves if strand == '+' else resolver.reverse(leaves)
        low, high = ((max(0, boundary-size), boundary) if side == 'left'
                     else (boundary, min(length, boundary+size)))
        leaves = self.expand(name, low, high)
        remaining = size - (high-low)
        if not remaining:
            return leaves
        cached = self.flanks.get(name)
        if cached is not None and remaining <= self.anchor:
            walk = cached[0 if side == 'left' else 1]
            take = min(remaining, walk.length)
            outside = (walk.slice(walk.length-take, walk.length) if side == 'left'
                       else walk.slice(0, take))
        elif name in resolver.events:
            start, end = resolver.breakpoints(name)
            outside = self.flank(resolver.events[name].chrom, start if side == 'left' else end,
                                 remaining, side)
        else:
            lift = resolver.roots[name].lift
            attachment = lift.get(side) if lift else None
            if not attachment:
                return leaves
            target, point, strand = attachment
            target_side = side if strand == '+' else ('right' if side == 'left' else 'left')
            outside = self.flank(target, point, remaining, target_side)
            if strand == '-':
                outside = resolver.reverse(outside)
        return outside + leaves if side == 'left' else leaves + outside
