"""Parent-anchored coordinates for the minsetref cohort rGFA exporter."""
from collections import defaultdict
import heapq

from gfa_interval_pipeline import Resolver


class StableResolver(Resolver):
    def __init__(self, events, roots, header_lengths, hits, reachable,
                 nested_pos_base=1, source_aliases=None):
        super().__init__(events, roots, header_lengths, hits)
        self.source_aliases = source_aliases or {}
        self.nested_pos_base = nested_pos_base
        # Link-only walks through carriers' abutting alleles (gfa_junctions).
        self.junctions = []
        self.ranks = {name: 0 if root.kind == 'reference' else 1
                      for name, root in roots.items()}
        self.order = self._dependency_order(reachable)
        next_rank = max(self.ranks.values(), default=0) + 1
        for name in self.order:
            self.ranks[name] = next_rank
            next_rank += 1
            self.breakpoints(name)
            for run in events[name].runs:
                if run.target:
                    self.target_interval(run.target, run.rstart, run.rend,
                                         run.orientation)

    def _dependency_order(self, reachable):
        """Order parents before children, independently of VCF record order."""
        children = defaultdict(list)
        pending = {}
        ready = []
        for name in reachable:
            if name in self.roots:
                raise ValueError(f'insertion ID {name!r} also names a FASTA root')
            event = self.events[name]
            targets = {event.chrom}
            targets.update(run.target for run in event.runs if run.target)
            dependencies = set()
            for target in targets:
                if target in self.events:
                    if target not in reachable:
                        raise ValueError(f'{name}: missing dependency {target!r}')
                    dependencies.add(target)
                elif target not in self.roots:
                    raise ValueError(f'{name}: graph target {target!r} has no FASTA '
                                     'or insertion definition; supply -r/-a, or '
                                     '--local-path-fasta for a source-mapped local target')
            pending[name] = len(dependencies)
            for parent in dependencies:
                children[parent].append(name)
            if not dependencies:
                heapq.heappush(ready, (event.order, name))
        ordered = []
        while ready:
            _ordinal, name = heapq.heappop(ready)
            ordered.append(name)
            for child in children[name]:
                pending[child] -= 1
                if not pending[child]:
                    heapq.heappush(ready, (self.events[child].order, child))
        if len(ordered) != len(reachable):
            names = sorted(name for name, count in pending.items() if count)
            raise ValueError('cyclic insertion/parent graph involving: ' +
                             ', '.join(names[:10]))
        return ordered

    def breakpoints(self, name):
        event = self.events[name]
        if getattr(event, 'kind', 'insertion') == 'snp':
            start, end = event.pos - 1, event.pos
            length = self.target_length(event.chrom)
            if not 0 <= start < end <= length:
                raise ValueError(f'{name}: SNP coordinate is outside parent {event.chrom!r}')
            return start, end
        parent = self.events.get(event.chrom)
        if parent is not None and getattr(parent, 'kind', 'insertion') == 'deletion':
            # A deletion row has no bases: rows nested in it (--exact: a
            # member's kept bases) replace the deleted interval, entering
            # from its left flank and leaving to its right flank.
            if getattr(event, 'kind', 'insertion') != 'insertion':
                raise ValueError(f'{name}: only insertions can be nested in deletion '
                                 f'{event.chrom!r}')
            offset = event.pos - self.nested_pos_base
            deleted = sum(run.rend - run.rstart for run in parent.runs if run.operation == 'D')
            if not 0 <= offset <= max(deleted, parent.ref_end - parent.pos):
                raise ValueError(f'{name}: offset {offset} is outside deletion {event.chrom!r}')
            return 0, 0
        # Rows on a --exact duplication template path use nested POS too.
        nested = parent is not None or (event.chrom in self.roots and
                                         self.roots[event.chrom].kind == 'duplication')
        start = event.pos - (self.nested_pos_base if nested else 0)
        # The cohort writer uses POS = offset + 1 for nested events. For a
        # pure insertion END == POS; a nonempty replacement ends at END.
        end = start if event.ref_end == event.pos else event.ref_end
        if getattr(event, 'kind', 'insertion') == 'deletion':
            # A merged DEL's END spans the cohort, not the representative.
            # A nonzero D span carries abs(SVLEN); old rows without SVLEN
            # retain their END-based interval.
            deleted_length = sum(run.rend - run.rstart for run in event.runs
                                 if run.operation == 'D')
            if deleted_length:
                end = start + deleted_length
        length = self.target_length(event.chrom)
        if not 0 <= start <= end <= length:
            hint = 'check POS/END/SVLEN and the parent length'
            if nested:
                hint += ', including --nested-pos-base'
            raise ValueError(f'{name}: parent interval {event.chrom}:{start}-{end} '
                             f'is outside 0-{length}; {hint}')
        return start, end

    def expand(self, name, start, end, orientation='+', stack=()):
        if name not in self.events:
            if name not in self.roots:
                raise ValueError(f'graph target {name!r} has no FASTA sequence')
            length = self.roots[name].length
            if not 0 <= start <= end <= length:
                raise ValueError(f'{name}: interval {start}-{end} is outside 0-{length}')
        if start == end:
            return []
        if name in self.source_aliases:
            name, start, end, orientation = self.canonical_interval(name, start, end, orientation)
        return super().expand(name, start, end, orientation, stack)

    def canonical_interval(self, name, start, end, orientation='+'):
        alias = self.source_aliases.get(name)
        return (alias.interval(start, end, orientation) if alias else
                (name, start, end, orientation))

    def _flank(self, name, boundary, size, side):
        """Walk beyond an insertion end through that insertion's parent."""
        length = self.target_length(name)
        if not 0 <= boundary <= length:
            raise ValueError(f'{name}: flank boundary {boundary} is out of bounds')
        if name in self.source_aliases and not self.roots[name].lift:
            target, point, _end, strand = self.canonical_interval(name, boundary, boundary)
            target_side = side if strand == '+' else ('right' if side == 'left' else 'left')
            leaves = self._flank(target, point, size, target_side)
            return leaves if strand == '+' else self.reverse(leaves)
        if side == 'left':
            low, high = max(0, boundary-size), boundary
        else:
            low, high = boundary, min(length, boundary+size)
        leaves = self.expand(name, low, high)
        remaining = size - (high-low)
        if remaining and name in self.events:
            start, end = self.breakpoints(name)
            parent = self.events[name].chrom
            outside = self._flank(parent, start if side == 'left' else end,
                                  remaining, side)
            leaves = outside + leaves if side == 'left' else leaves + outside
        elif remaining and name in self.roots and self.roots[name].lift:
            attachment = self.roots[name].lift.get(side)
            if attachment:
                target, point, strand = attachment
                target_side = side if strand == '+' else ('right' if side == 'left' else 'left')
                outside = self._flank(target, point, remaining, target_side)
                if strand == '-':
                    outside = self.reverse(outside)
                leaves = outside + leaves if side == 'left' else leaves + outside
        return leaves

    def junction_leaves(self, number, flank):
        """Link walk of one junction: the parent's flank, the carrier's
        alleles of abutting variants in its order, the parent's flank."""
        parent, start, end, items = self.junctions[number]
        leaves = self._flank(parent, start, flank, 'left')
        for item in items:
            leaves.extend(self.allele_leaves(item))
        return leaves + self._flank(parent, end, flank, 'right')

    def allele_leaves(self, item):
        """Leaves of one carrier's allele: (variant, nested items) with the
        nested items already in the carrier's order."""
        name, nested = item
        event = self.events[name]
        leaves = []
        if getattr(event, 'kind', 'insertion') == 'deletion':
            for child in nested:
                leaves.extend(self.allele_leaves(child))
            return leaves
        cursor = 0
        for child in nested:
            start, end = self.breakpoints(child[0])
            if start < cursor:
                continue
            leaves.extend(self.expand(name, cursor, start))
            leaves.extend(self.allele_leaves(child))
            cursor = end
        return leaves + self.expand(name, cursor, event.length)

    def local_path(self, name, start, end, flank):
        return (self._flank(name, start, flank, 'left') +
                self.expand(name, start, end) +
                self._flank(name, end, flank, 'right'))
