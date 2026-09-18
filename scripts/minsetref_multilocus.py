#!/usr/bin/env python3
"""Reporting-only multi-locus alignment resolution.

This module deliberately does not decide cleaning residuals.  It mirrors the
query-coordinate sweep used by the older graphcigarlight ``graphicalign``
implementation while retaining exact paired query/reference blocks needed for
graphic-CIGAR construction.
"""
from __future__ import annotations

import dataclasses
import heapq
from collections import defaultdict
from typing import DefaultDict, Dict, Iterable, List, Sequence, Tuple


PairedBlock = Tuple[int, int, int, int]


@dataclasses.dataclass(frozen=True)
class PlacementEvidence:
    query_id: str
    mapped_partition: str
    ref_haplotype: str
    ref_contig: str
    strand: str
    target_record: str
    query_start: int
    query_end: int
    paired_blocks: Tuple[PairedBlock, ...]
    alignment_score: float
    aligned_bases: int
    identity: float
    source: str
    ordinal: int


@dataclasses.dataclass(frozen=True)
class ResolvedPlacementComponent:
    query_id: str
    component_index: int
    mapped_partition: str
    ref_haplotype: str
    ref_contig: str
    ref_start: int
    ref_end: int
    strand: str
    target_record: str
    query_start: int
    query_end: int
    query_blocks: Tuple[Tuple[int, int], ...]
    target_blocks: Tuple[Tuple[int, int], ...]
    paired_blocks: Tuple[PairedBlock, ...]
    alignment_score: float
    aligned_bases: int
    identity: float
    source: str
    evidence_ordinal: int
    liftover_path: str = "."


def merge_intervals(intervals: Iterable[Tuple[int, int]]) -> Tuple[Tuple[int, int], ...]:
    ordered = sorted((int(a), int(b)) for a, b in intervals if int(b) > int(a))
    merged: List[List[int]] = []
    for start, end in ordered:
        if merged and start <= merged[-1][1]:
            merged[-1][1] = max(merged[-1][1], end)
        else:
            merged.append([start, end])
    return tuple((start, end) for start, end in merged)


def _clip_paired_block(
    block: PairedBlock,
    strand: str,
    query_start: int,
    query_end: int,
) -> PairedBlock | None:
    q0, q1, r0, _r1 = block
    left = max(q0, query_start)
    right = min(q1, query_end)
    if right <= left:
        return None
    if strand == "+":
        ref_start = r0 + left - q0
        ref_end = r0 + right - q0
    else:
        ref_start = r0 + q1 - right
        ref_end = r0 + q1 - left
    return left, right, ref_start, ref_end


def _priority(evidence: PlacementEvidence, index: int) -> Tuple:
    return (
        -float(evidence.alignment_score),
        -int(evidence.aligned_bases),
        -float(evidence.identity),
        evidence.ordinal,
        evidence.ref_haplotype,
        evidence.ref_contig,
        evidence.strand,
        evidence.target_record,
        evidence.mapped_partition,
        index,
    )


def _validate_evidence(evidence: PlacementEvidence) -> None:
    if evidence.query_end <= evidence.query_start:
        raise ValueError(
            f"{evidence.query_id}: invalid evidence query span "
            f"{evidence.query_start}-{evidence.query_end}"
        )
    if evidence.strand not in {"+", "-"}:
        raise ValueError(
            f"{evidence.query_id}: invalid evidence strand {evidence.strand!r}"
        )
    for q0, q1, r0, r1 in evidence.paired_blocks:
        if (
            q0 < evidence.query_start
            or q1 > evidence.query_end
            or q1 <= q0
            or r1 <= r0
            or q1 - q0 != r1 - r0
        ):
            raise ValueError(
                f"{evidence.query_id}: invalid paired block "
                f"{q0}-{q1}:{r0}-{r1}"
            )


def resolve_query_evidence(
    evidence_rows: Sequence[PlacementEvidence],
) -> List[ResolvedPlacementComponent]:
    """Select the highest-priority active alignment along one query.

    The sweep runs on full clipped alignment spans, matching the old
    ``graphicalign`` policy. Exact paired blocks are clipped afterwards, so
    insertions remain represented as query gaps rather than mapped bases.
    """
    if not evidence_rows:
        return []
    query_ids = {row.query_id for row in evidence_rows}
    if len(query_ids) != 1:
        raise ValueError("resolve_query_evidence requires exactly one query_id")
    starts: DefaultDict[int, List[int]] = defaultdict(list)
    ends: DefaultDict[int, List[int]] = defaultdict(list)
    coordinates = set()
    for index, evidence in enumerate(evidence_rows):
        _validate_evidence(evidence)
        starts[evidence.query_start].append(index)
        ends[evidence.query_end].append(index)
        coordinates.add(evidence.query_start)
        coordinates.add(evidence.query_end)
    ordered = sorted(coordinates)
    active = set()
    heap: List[Tuple] = []
    selected: List[List[int]] = []
    for coordinate, next_coordinate in zip(ordered, ordered[1:]):
        # Half-open intervals ending here are inactive before intervals that
        # begin here are considered.
        for index in ends.get(coordinate, ()):
            active.discard(index)
        for index in starts.get(coordinate, ()):
            active.add(index)
            heapq.heappush(heap, _priority(evidence_rows[index], index))
        while heap and heap[0][-1] not in active:
            heapq.heappop(heap)
        if not heap or next_coordinate <= coordinate:
            continue
        winner = heap[0][-1]
        if selected and selected[-1][2] == winner and selected[-1][1] == coordinate:
            selected[-1][1] = next_coordinate
        else:
            selected.append([coordinate, next_coordinate, winner])

    components: List[ResolvedPlacementComponent] = []
    for selected_start, selected_end, winner in selected:
        evidence = evidence_rows[winner]
        paired = tuple(
            clipped
            for block in evidence.paired_blocks
            for clipped in (
                _clip_paired_block(
                    block, evidence.strand, selected_start, selected_end,
                ),
            )
            if clipped is not None
        )
        if not paired:
            # A query-only insertion inside an alignment is not mapped
            # reference sequence and must not become a placement component.
            continue
        query_blocks = merge_intervals((q0, q1) for q0, q1, _r0, _r1 in paired)
        target_blocks = merge_intervals((r0, r1) for _q0, _q1, r0, r1 in paired)
        aligned_bases = sum(q1 - q0 for q0, q1 in query_blocks)
        proportional_score = (
            float(evidence.alignment_score)
            * aligned_bases
            / max(1, int(evidence.aligned_bases))
        )
        components.append(ResolvedPlacementComponent(
            evidence.query_id,
            len(components) + 1,
            evidence.mapped_partition,
            evidence.ref_haplotype,
            evidence.ref_contig,
            min(start for start, _end in target_blocks),
            max(end for _start, end in target_blocks),
            evidence.strand,
            evidence.target_record,
            selected_start,
            selected_end,
            query_blocks,
            target_blocks,
            paired,
            proportional_score,
            aligned_bases,
            float(evidence.identity),
            evidence.source,
            evidence.ordinal,
        ))
    # Match graphcigarlight/connectbreaks: when a winner changes but the next
    # selected piece continues monotonically on the same target path, retain
    # one component and encode the intervening query/reference gap in its
    # graphic CIGAR instead of emitting duplicate target markers.
    connected: List[ResolvedPlacementComponent] = []
    for component in components:
        if not connected:
            connected.append(component)
            continue
        prior = connected[-1]
        same_target = (
            prior.mapped_partition == component.mapped_partition
            and prior.ref_haplotype == component.ref_haplotype
            and prior.ref_contig == component.ref_contig
            and prior.strand == component.strand
            and prior.target_record == component.target_record
        )
        prior_blocks = sorted(prior.paired_blocks)
        next_blocks = sorted(component.paired_blocks)
        monotonic = False
        if same_target and prior_blocks and next_blocks:
            if prior.strand == "+":
                monotonic = next_blocks[0][2] >= prior_blocks[-1][3]
            else:
                monotonic = next_blocks[0][3] <= prior_blocks[-1][2]
        if not monotonic:
            connected.append(component)
            continue
        paired = tuple((*prior_blocks, *next_blocks))
        query_blocks = merge_intervals(
            (q0, q1) for q0, q1, _r0, _r1 in paired
        )
        target_blocks = merge_intervals(
            (r0, r1) for _q0, _q1, r0, r1 in paired
        )
        aligned = prior.aligned_bases + component.aligned_bases
        identity = (
            prior.identity * prior.aligned_bases
            + component.identity * component.aligned_bases
        ) / max(1, aligned)
        sources = list(dict.fromkeys((
            *prior.source.split(","), *component.source.split(","),
        )))
        connected[-1] = dataclasses.replace(
            prior,
            query_end=component.query_end,
            ref_start=min(start for start, _end in target_blocks),
            ref_end=max(end for _start, end in target_blocks),
            query_blocks=query_blocks,
            target_blocks=target_blocks,
            paired_blocks=paired,
            alignment_score=(
                prior.alignment_score + component.alignment_score
            ),
            aligned_bases=aligned,
            identity=identity,
            source=",".join(sources),
            evidence_ordinal=min(
                prior.evidence_ordinal, component.evidence_ordinal,
            ),
        )
    return [
        dataclasses.replace(component, component_index=index)
        for index, component in enumerate(connected, 1)
    ]


def resolve_all_evidence(
    evidence_rows: Iterable[PlacementEvidence],
) -> Dict[str, List[ResolvedPlacementComponent]]:
    by_query: DefaultDict[str, List[PlacementEvidence]] = defaultdict(list)
    for row in evidence_rows:
        by_query[row.query_id].append(row)
    return {
        query_id: resolve_query_evidence(rows)
        for query_id, rows in sorted(by_query.items())
    }
