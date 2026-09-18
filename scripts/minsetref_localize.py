#!/usr/bin/env python3
"""Localize a novel insertion to main-chrom coordinates via flank anchoring.

For an insertion occupying query interval [qi_s, qi_e) on an assembly contig,
we look at the main-chrom alignments in the +/- ``flank`` bp windows adjacent to
the insertion.  Each aligned M/=/X block maps query positions to main positions.
Blocks are clustered by nearby (main_chrom, strand) loci. Within each locus,
overlapping query coverage is assigned to the best raw alignment, producing a
non-overlapping mosaic. The mosaic is scored for identity, query gaps, reference
gaps, and path consistency before the best location is selected.

  * up & down anchors agree (same chrom+strand) -> ``mapped``: the insertion is
    the smallest main interval between the upstream anchor end and downstream
    anchor start, tightened inward.  A clean insertion collapses to [X, X].
  * anchors disagree -> ``both``: record each side's placement.
  * neither side anchors within ``flank`` -> ``unmapped``.

Design mirrors graphpathlift_corrected.py's group->score->best-location idea, but
specialized to two single-sided flank anchors instead of whole-query chaining.
"""
from __future__ import annotations

import dataclasses
import heapq
from typing import Dict, List, Optional, Sequence, Tuple


@dataclasses.dataclass
class AnchorBlock:
    """One M/=/X block on a contig, mapping query [q0,q1) to main [t0,t1)."""
    q0: int
    q1: int
    t0: int
    t1: int
    main: str
    strand: str
    identity: float = 100.0
    alignment_score: float = 0.0

    def map_qpos(self, qpos: int) -> int:
        """Map a query position (clamped into the block) to a main position."""
        p = min(max(qpos, self.q0), self.q1)
        if self.strand == "-":
            return self.t0 + (self.q1 - p)  # q1<->t0, q0<->t1
        return self.t0 + (p - self.q0)

    @property
    def priority(self) -> Tuple[float, float, int]:
        """Priority used when raw alignments overlap on the query."""
        return (self.alignment_score, self.identity, self.q1 - self.q0)


@dataclasses.dataclass
class Placement:
    main: str
    start: int
    end: int
    strand: str


@dataclasses.dataclass
class LocalizeResult:
    status: str                     # "mapped" | "one_sided" | "both" | "unmapped"
    placements: List[Placement]
    note: str = ""


def overlapping_blocks_for_windows(
    windows: Sequence[Tuple[str, int, int]],
    blocks_by_contig: Dict[str, List[AnchorBlock]],
) -> List[List[AnchorBlock]]:
    """Join query windows to overlapping blocks with a coordinate sweep.

    The old localization path scanned every block on a contig for every
    insertion.  Here half-open window and block endpoints are sorted once per
    contig.  Work is O((W+B) log(W+B) + K), where K is the number of real
    window/block overlaps that must be returned.

    Identical windows are joined once and share the resulting read-only block
    list.  This matters for placement clones that retain the same assembly
    coordinates.
    """
    output: List[List[AnchorBlock]] = [[] for _ in windows]
    members_by_contig: Dict[str, Dict[Tuple[int, int], List[int]]] = {}
    for index, (contig, raw_start, raw_end) in enumerate(windows):
        start, end = int(raw_start), int(raw_end)
        if end <= start:
            continue
        members_by_contig.setdefault(contig, {}).setdefault((start, end), []).append(index)

    for contig, member_map in members_by_contig.items():
        blocks = [b for b in blocks_by_contig.get(contig, []) if b.q1 > b.q0]
        if not blocks:
            continue
        unique_windows = list(member_map)
        overlap_indices: List[List[int]] = [[] for _ in unique_windows]

        window_starts: Dict[int, List[int]] = {}
        window_ends: Dict[int, List[int]] = {}
        for local_index, (start, end) in enumerate(unique_windows):
            window_starts.setdefault(start, []).append(local_index)
            window_ends.setdefault(end, []).append(local_index)
        block_starts: Dict[int, List[int]] = {}
        block_ends: Dict[int, List[int]] = {}
        for block_index, block in enumerate(blocks):
            block_starts.setdefault(block.q0, []).append(block_index)
            block_ends.setdefault(block.q1, []).append(block_index)

        coordinates = sorted(
            set(window_starts) | set(window_ends) | set(block_starts) | set(block_ends)
        )
        active_windows: set[int] = set()
        active_blocks: set[int] = set()
        for coordinate in coordinates:
            # Half-open intervals: intervals ending here must disappear before
            # intervals starting here are joined.
            for local_index in window_ends.get(coordinate, ()):
                active_windows.discard(local_index)
            for block_index in block_ends.get(coordinate, ()):
                active_blocks.discard(block_index)

            # A newly starting block overlaps windows that started earlier.
            for block_index in block_starts.get(coordinate, ()):
                for local_index in active_windows:
                    overlap_indices[local_index].append(block_index)
                active_blocks.add(block_index)

            # A newly starting window sees prior blocks plus blocks starting at
            # this same coordinate. Every overlap is therefore emitted once.
            for local_index in window_starts.get(coordinate, ()):
                overlap_indices[local_index].extend(active_blocks)
                active_windows.add(local_index)

        for local_index, window in enumerate(unique_windows):
            # Restore the input/PAF order so exact priority ties retain the same
            # deterministic winner as the pre-index localization code.
            selected = [blocks[i] for i in sorted(overlap_indices[local_index])]
            for global_index in member_map[window]:
                output[global_index] = selected
    return output


def localize_windows_from_blocks(
    queries: Sequence[Tuple[str, int, int]],
    blocks_by_contig: Dict[str, List[AnchorBlock]],
    flank: int,
) -> List[LocalizeResult]:
    """Localize many intervals while streaming completed sweep windows.

    Unlike :func:`overlapping_blocks_for_windows`, this does not retain all K
    query/block overlaps for a contig simultaneously. A window is localized
    and released as soon as its right endpoint is reached. Identical source
    intervals are evaluated once.
    """
    results: List[Optional[LocalizeResult]] = [None] * len(queries)
    members_by_contig: Dict[str, Dict[Tuple[int, int], List[int]]] = {}
    for index, (contig, raw_start, raw_end) in enumerate(queries):
        start, end = int(raw_start), int(raw_end)
        members_by_contig.setdefault(contig, {}).setdefault((start, end), []).append(index)

    for contig, member_map in members_by_contig.items():
        blocks = [block for block in blocks_by_contig.get(contig, []) if block.q1 > block.q0]
        unique_queries = list(member_map)
        if not blocks:
            for query in unique_queries:
                start, end = query
                result = localize_insertion(start, end, [], flank)
                for global_index in member_map[query]:
                    results[global_index] = result
            continue

        window_starts: Dict[int, List[int]] = {}
        window_ends: Dict[int, List[int]] = {}
        for local_index, (start, end) in enumerate(unique_queries):
            window_start = max(0, start - flank)
            window_end = end + flank
            window_starts.setdefault(window_start, []).append(local_index)
            window_ends.setdefault(window_end, []).append(local_index)
        block_starts: Dict[int, List[int]] = {}
        block_ends: Dict[int, List[int]] = {}
        for block_index, block in enumerate(blocks):
            block_starts.setdefault(block.q0, []).append(block_index)
            block_ends.setdefault(block.q1, []).append(block_index)

        coordinates = sorted(
            set(window_starts) | set(window_ends) | set(block_starts) | set(block_ends)
        )
        active_blocks: set[int] = set()
        active_windows: set[int] = set()
        overlap_indices: Dict[int, List[int]] = {}
        for coordinate in coordinates:
            # Finalize ending half-open windows before admitting blocks that
            # start at this same coordinate.
            for local_index in window_ends.get(coordinate, ()):
                indices = sorted(overlap_indices.pop(local_index, ()))
                selected = [blocks[index] for index in indices]
                start, end = unique_queries[local_index]
                result = localize_insertion(start, end, selected, flank)
                for global_index in member_map[(start, end)]:
                    results[global_index] = result
                active_windows.discard(local_index)
            for block_index in block_ends.get(coordinate, ()):
                active_blocks.discard(block_index)
            for block_index in block_starts.get(coordinate, ()):
                for local_index in active_windows:
                    overlap_indices[local_index].append(block_index)
                active_blocks.add(block_index)
            for local_index in window_starts.get(coordinate, ()):
                overlap_indices[local_index] = list(active_blocks)
                active_windows.add(local_index)

    return [
        result if result is not None else LocalizeResult("unmapped", [], note="no flank anchor")
        for result in results
    ]


def anchor_blocks_from_hits(hits: Sequence["object"], contig: str) -> List[AnchorBlock]:
    """Flatten main-chrom hits for one contig into AnchorBlocks (contig coords).

    ``hits`` are AlignmentHit records whose query is the contig and whose target
    is a main chrom; ``aligned_pairs`` give (q0,q1,t0,t1) forward-axis blocks.
    """
    blocks: List[AnchorBlock] = []
    for hit in hits:
        if hit.query_id != contig:
            continue
        for q0, q1, t0, t1 in hit.aligned_pairs:
            if q1 > q0:
                blocks.append(AnchorBlock(
                    q0, q1, t0, t1, hit.target_id, hit.strand,
                    float(hit.identity), float(hit.identity) * float(hit.aligned_bases) / 100.0,
                ))
    return blocks


@dataclasses.dataclass
class AnchorCandidate:
    main: str
    strand: str
    blocks: List[AnchorBlock]
    unique_bp: int
    score: float
    ref_start: int
    ref_end: int


def _clip_block(block: AnchorBlock, q0: int, q1: int) -> Optional[AnchorBlock]:
    a, b = max(block.q0, q0), min(block.q1, q1)
    if b <= a:
        return None
    ta, tb = block.map_qpos(a), block.map_qpos(b)
    return AnchorBlock(a, b, min(ta, tb), max(ta, tb), block.main, block.strand,
                       block.identity, block.alignment_score)


def _reference_clusters(blocks: Sequence[AnchorBlock], merge_gap: int) -> List[List[AnchorBlock]]:
    """Cluster one flank by chromosome, strand, and nearby reference locus."""
    grouped: Dict[Tuple[str, str], List[AnchorBlock]] = {}
    for block in blocks:
        grouped.setdefault((block.main, block.strand), []).append(block)
    clusters: List[List[AnchorBlock]] = []
    for vals in grouped.values():
        vals = sorted(vals, key=lambda b: (b.t0, b.t1, b.q0, b.q1))
        current: List[AnchorBlock] = []
        current_end = -1
        for block in vals:
            if not current or block.t0 - current_end <= merge_gap:
                current.append(block)
                current_end = max(current_end, block.t1)
            else:
                clusters.append(current)
                current = [block]
                current_end = block.t1
        if current:
            clusters.append(current)
    return clusters


class _MaxPriority:
    """Heap entry whose smallest item is the largest localization priority."""

    __slots__ = ("key", "index")

    def __init__(self, key: tuple, index: int):
        self.key = key
        self.index = index

    def __lt__(self, other: "_MaxPriority") -> bool:
        if self.key != other.key:
            return self.key > other.key
        return self.index < other.index


def _resolve_query_overlaps(blocks: Sequence[AnchorBlock]) -> List[AnchorBlock]:
    """Build a best-over-all-alignments, non-overlapping query mosaic.

    This is the flank-sized equivalent of graphpathlift_corrected.py's
    resolve_overlap_greedy: every elementary query span is assigned to the
    highest-priority alignment covering it, and adjacent pieces from the same
    alignment are coalesced.
    """
    valid = [(index, block) for index, block in enumerate(blocks) if block.q1 > block.q0]
    if not valid:
        return []
    starts: Dict[int, List[int]] = {}
    ends: Dict[int, List[int]] = {}
    active = [False] * len(blocks)
    by_index = {index: block for index, block in valid}
    for index, block in valid:
        starts.setdefault(block.q0, []).append(index)
        ends.setdefault(block.q1, []).append(index)
    coords = sorted(set(starts) | set(ends))
    heap: List[_MaxPriority] = []
    chosen: List[AnchorBlock] = []
    for a, b in zip(coords, coords[1:]):
        for index in ends.get(a, ()):
            active[index] = False
        for index in starts.get(a, ()):
            block = by_index[index]
            key = (block.priority, -block.t0, block.main, block.strand)
            active[index] = True
            heapq.heappush(heap, _MaxPriority(key, index))
        while heap and not active[heap[0].index]:
            heapq.heappop(heap)
        if b <= a:
            continue
        if not heap:
            continue
        best = by_index[heap[0].index]
        piece = _clip_block(best, a, b)
        if piece is None:
            continue
        if (chosen and chosen[-1].main == piece.main and chosen[-1].strand == piece.strand
                and chosen[-1].q1 == piece.q0 and chosen[-1].identity == piece.identity
                and chosen[-1].alignment_score == piece.alignment_score
                and (chosen[-1].t1 == piece.t0 or piece.t1 == chosen[-1].t0)):
            prev = chosen[-1]
            chosen[-1] = AnchorBlock(prev.q0, piece.q1, min(prev.t0, piece.t0),
                                     max(prev.t1, piece.t1), prev.main, prev.strand,
                                     prev.identity, prev.alignment_score)
        else:
            chosen.append(piece)
    return chosen


def _best_candidate(blocks: Sequence[AnchorBlock], merge_gap: int) -> Optional[AnchorCandidate]:
    candidates: List[AnchorCandidate] = []
    for cluster in _reference_clusters(blocks, merge_gap):
        mosaic = _resolve_query_overlaps(cluster)
        if not mosaic:
            continue
        unique_bp = sum(b.q1 - b.q0 for b in mosaic)
        score = _mosaic_score(mosaic)
        candidates.append(AnchorCandidate(
            mosaic[0].main, mosaic[0].strand, mosaic, unique_bp, score,
            min(b.t0 for b in mosaic), max(b.t1 for b in mosaic),
        ))
    if not candidates:
        return None
    # After overlap resolution, choose the highest-scoring mosaic, as in the
    # final graphpathlift candidate ordering. Unique coverage and interval
    # tightness are deterministic tie-breakers.
    return max(candidates, key=lambda c: (c.score, c.unique_bp,
                                           -(c.ref_end - c.ref_start), c.main, c.strand))


def _mosaic_score(blocks: Sequence[AnchorBlock]) -> float:
    """Graphpathlift-style score for an already non-overlapping query mosaic."""
    ordered = sorted(blocks, key=lambda b: (b.q0, b.q1, b.t0, b.t1))
    score = sum((b.q1 - b.q0) * b.identity / 100.0 for b in ordered)
    for left, right in zip(ordered, ordered[1:]):
        qgap = max(0, right.q0 - left.q1)
        if qgap:
            score -= qgap // 10 + 20
        if left.strand == "+":
            if right.t0 >= left.t1:
                rgap = right.t0 - left.t1
            else:
                # Overlap/backtracking is not a coherent forward path.
                rgap = left.t1 - right.t0
        else:
            if right.t1 <= left.t0:
                rgap = left.t0 - right.t1
            else:
                rgap = right.t1 - left.t0
        if rgap:
            score -= rgap + 20
    return score


def _upstream_anchor(blocks: Sequence[AnchorBlock], qi_s: int,
                     merge_gap: int) -> Optional[Tuple[str, str, int]]:
    """Main (chrom, strand, pos) for the upstream boundary qi_s."""
    win = [b for b in blocks if b.q1 <= qi_s or (b.q0 < qi_s <= b.q1)]
    candidate = _best_candidate(win, merge_gap)
    if candidate is None:
        return None
    # block whose end is closest to (but not past) qi_s -> tightest inward anchor
    anchor = min(candidate.blocks, key=lambda b: (qi_s - min(b.q1, qi_s), -(b.q1 - b.q0)))
    return candidate.main, candidate.strand, anchor.map_qpos(qi_s)


def _downstream_anchor(blocks: Sequence[AnchorBlock], qi_e: int,
                       merge_gap: int) -> Optional[Tuple[str, str, int]]:
    """Main (chrom, strand, pos) for the downstream boundary qi_e."""
    win = [b for b in blocks if b.q0 >= qi_e or (b.q0 <= qi_e < b.q1)]
    candidate = _best_candidate(win, merge_gap)
    if candidate is None:
        return None
    anchor = min(candidate.blocks, key=lambda b: (max(b.q0, qi_e) - qi_e, -(b.q1 - b.q0)))
    return candidate.main, candidate.strand, anchor.map_qpos(qi_e)


def localize_insertion(qi_s: int, qi_e: int, blocks: Sequence[AnchorBlock],
                       flank: int = 5000, require_unique: bool = False) -> LocalizeResult:
    """Localize insertion [qi_s, qi_e) using flank AnchorBlocks (contig coords)."""
    up_win = [clipped for b in blocks
              if (clipped := _clip_block(b, qi_s - flank, qi_s)) is not None]
    dn_win = [clipped for b in blocks
              if (clipped := _clip_block(b, qi_e, qi_e + flank)) is not None]

    if require_unique:
        for window in (up_win, dn_win):
            scores = []
            for cluster in _reference_clusters(window, flank):
                mosaic = _resolve_query_overlaps(cluster)
                if mosaic:
                    scores.append((_mosaic_score(mosaic), sum(b.q1-b.q0 for b in mosaic),
                                   -(max(b.t1 for b in mosaic)-min(b.t0 for b in mosaic))))
            if scores and scores.count(max(scores)) > 1:
                return LocalizeResult('unmapped', [], note='ambiguous flank placements')

    up = _upstream_anchor(up_win, qi_s, flank)
    dn = _downstream_anchor(dn_win, qi_e, flank)

    if require_unique:
        # A nearby flank localizes a region, but is not an exact graph edge.
        # Gaps between an anchor and the core must not delete source bases.
        if up is not None and not any(b.q1 == qi_s and
                (b.main, b.strand, b.map_qpos(qi_s)) == up for b in up_win):
            up = None
        if dn is not None and not any(b.q0 == qi_e and
                (b.main, b.strand, b.map_qpos(qi_e)) == dn for b in dn_win):
            dn = None

    if up is None and dn is None:
        return LocalizeResult("unmapped", [], note="no flank anchor")

    if up is not None and dn is not None:
        up_main, up_strand, up_pos = up
        dn_main, dn_strand, dn_pos = dn
        if up_main == dn_main and up_strand == dn_strand:
            if require_unique and ((up_strand == '+' and up_pos > dn_pos) or
                                   (up_strand == '-' and up_pos < dn_pos)):
                return LocalizeResult('unmapped', [], note='noncollinear flank placements')
            s, e = min(up_pos, dn_pos), max(up_pos, dn_pos)
            return LocalizeResult("mapped", [Placement(up_main, s, e, up_strand)])
        # disagree: record both single-sided placements
        return LocalizeResult("both", [
            Placement(up_main, up_pos, up_pos, up_strand),
            Placement(dn_main, dn_pos, dn_pos, dn_strand),
        ], note="up/down disagree")

    # one-sided anchor: record the single breakpoint we do have
    side = up if up is not None else dn
    main, strand, pos = side
    return LocalizeResult("one_sided", [Placement(main, pos, pos, strand)],
                          note="one_side_up" if up is not None else "one_side_down")
