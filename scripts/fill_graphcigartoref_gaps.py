#!/usr/bin/env python3
"""Fill internal query-annotation gaps in per-sample graph-CIGAR output.

The input is the seven-column table emitted by
``graphcigartoref_persample.py``.  Existing rows are immutable: they are copied
to the output byte-for-byte and are never sliced, removed, or reassigned by
this additive stage.

Internal positive gaps on each query contig are selected from the union of the
effective column-3 query intervals.  Finite GenomeLift regions always mask
overlapping ``-inf`` rows; the remaining ``-inf`` rows are assigned by the
same complete-alignment score and >=1000 threshold used by graphreftovcf.py.
This changes only the analysis view: original rows are still copied unchanged.
Every gap shorter than 100 bp is selected; larger gaps require more than 100
uppercase A/C/G/T bases and an uppercase fraction greater than 0.2.  Gaps of
1 Mb or more are never selected.

The neighboring GenomeLift rows select the higher-priority side from columns
13/14.  When the two graph-CIGAR edges describe a positive gap on the same
reference contig whose span is plausible for the query gap, the complete query
and reference gaps are aligned with 100-bp anchors on both sides.  Otherwise
the reference is extended away from the higher-priority block using the same
size rule as the raw-reference flank rescue in
``graphcigartoref_persample.py``.  Anchor sequence participates in the
alignment but is removed from the emitted query interval and graph CIGAR.
Overlapping same-reference breakpoints are emitted as a pure insertion: the
insertion contains the query gap plus the lower-priority query sequence that
projects across the reference overlap, which represents a tandem duplication.

After query-gap alignment, positive same-reference gaps at exact query
breakpoints are subtracted from the union of all reference intervals
represented by the query.  Any uncovered piece is emitted as a pure deletion.
Pairs separated by unresolved query sequence, or overlapping on the query,
cannot establish a deletion and are skipped.  By default the output contains
the exact original input followed by new query-gap and deletion rows.
``--gaps-only`` writes only the additions.

The final step re-elects assembly PA ownership from the score-clipped graph
alignments. Mapped valid PAs rank first; qualifying former ``-inf`` PAs rank by
alignment score and are called resurrected PAs; clipped, insertion-only, and
wholly unmapped PA sequence ranks last. Alignments whose coordinate target is
an alternative/novel template receive half the corresponding ownership score,
and those PAs never receive ownership extensions. Gaps touching such a target
are not filled and do not generate deletion rows. Only one connected core is
retained per PA, then GenomeLift-style neighboring gap ownership is recomputed.
The sparse changed rows are written to ``--genomelift-fix-output`` for use
after the original GenomeLift file in a comma-ordered VCF ``--coord-map``.
"""

from __future__ import annotations

import argparse
import bisect
import dataclasses
import heapq
import multiprocessing as mp
import os
import re
import sys
import tempfile
import time
import traceback
from typing import Dict, Iterable, List, Mapping, Optional, Sequence, Tuple

import graphcigartoref as core
from graph_cigar_payloads import query_only_graph_cigar
from GenomeLift import (
    compute_priority_display_offsets_for_regions,
    parse_allele_name,
)
from graphcigartoref_persample import (
    SyntheticDeletionRow,
    _align_reference_flank_ops,
    read_genomelift,
    render_synthetic_deletion_or_warn,
)
from graphreftovcf import (
    MIN_INFINITE_FALLBACK_ALIGNMENT_SCORE,
    VALID_OWNERSHIP_PRIORITY_SCORE,
    apply_interval_extensions_to_coord,
    infinite_fallback_alignment_score,
    resolve_infinite_fallback_candidates,
)
from minsetref_segments import count_unmasked


DEFAULT_ANCHOR = 100
DEFAULT_MAXIMUM_GAP = 1_000_000
DEFAULT_MINIMUM_UNMASKED = 100
DEFAULT_MINIMUM_UNMASKED_FRACTION = 0.2
DEFAULT_REFERENCE_EXTENSION = 10_000
DEFAULT_MAX_IMPUTE = 10_000
OWNERSHIP_TIER_UNMAPPED = 0
OWNERSHIP_TIER_RESURRECTED = 1
OWNERSHIP_TIER_VALID_MAPPED = 2


@dataclasses.dataclass(frozen=True)
class AnnotationRow:
    line_number: int
    query_name: str
    reference_name: str
    query_coord: core.Coord
    query_alignment: str
    reference_alignment: str
    graph_cigar: str
    reference_coord: Optional[core.Coord] = None
    alternative_intervals: Optional[dict] = None


@dataclasses.dataclass(frozen=True)
class PriorityBlock:
    annotation: AnnotationRow
    allele: str
    locus: core.Coord
    facing_token: str


@dataclasses.dataclass(frozen=True)
class ReferenceBoundary:
    path: str
    position: int
    # Reference direction as physical query coordinates increase.
    strand: str


@dataclasses.dataclass(frozen=True)
class GapTask:
    index: int
    name: str
    contig: str
    start: int
    end: int
    left: AnnotationRow
    right: AnnotationRow
    selected: AnnotationRow
    selected_side: str
    left_boundary: ReferenceBoundary
    right_boundary: ReferenceBoundary
    mode: str
    reference_gap_size: Optional[int]
    # Tandem duplications include query sequence from the lower-priority block
    # in addition to the literal annotation gap.
    insertion_start: Optional[int] = None
    insertion_end: Optional[int] = None


@dataclasses.dataclass(frozen=True)
class GapResult:
    index: int
    row: Optional[str]
    error: str = ""


@dataclasses.dataclass(frozen=True)
class OwnershipLiftRow:
    line_number: int
    fields: Tuple[str, ...]
    allele: str
    locus: core.Coord
    part_index: int
    left_extension: str
    right_extension: str


@dataclasses.dataclass(frozen=True)
class OwnershipCandidate:
    start: int
    end: int
    allele: str
    # 2: mapped valid PA; 1: mapped resurrected PA; 0: unmapped sequence.
    tier: int
    score: int
    line_number: int
    non_primary: bool = False


_WORKER_QUERY_READER: Optional[core.FastaRegionReader] = None
_WORKER_REFERENCE_READER: Optional[core.FastaRegionReader] = None
_WORKER_ANCHOR = DEFAULT_ANCHOR
_WORKER_REFERENCE_EXTENSION = DEFAULT_REFERENCE_EXTENSION


def coord_text(coord: core.Coord) -> str:
    return f"{coord.chrom}:{coord.start}-{coord.end}{coord.strand}"


def _strip_outer_strand(name: str) -> str:
    return name.split("|outer_strand=", 1)[0]


def query_name_parts(name: str) -> Tuple[str, ...]:
    return tuple(
        part.strip()
        for part in _strip_outer_strand(name).split(";")
        if part.strip()
    )


def read_annotations(path: str) -> List[AnnotationRow]:
    rows: List[AnnotationRow] = []
    with open(path, "rt") as handle:
        for line_number, raw in enumerate(handle, 1):
            if not raw.strip() or raw.startswith("#"):
                continue
            fields = raw.rstrip("\r\n").split("\t")
            row = annotation_from_fields(fields, line_number, path)
            if row is None:
                # Explicit DEL rows have a zero-width query breakpoint and do
                # not annotate query sequence.
                continue
            rows.append(row)
    if not rows:
        raise ValueError(f"no positive-width query annotations found: {path}")
    return rows


def annotation_from_fields(
    fields: Sequence[str], line_number: int, source: str,
) -> Optional[AnnotationRow]:
    if len(fields) < 7:
        raise ValueError(
            f"{source}:{line_number}: expected the seven-column "
            f"graphcigartoref_persample layout, found {len(fields)} columns"
        )
    query_coord = core.parse_coord(fields[2])
    if query_coord is None:
        raise ValueError(
            f"{source}:{line_number}: malformed query coordinate {fields[2]!r}"
        )
    if query_coord.end <= query_coord.start:
        return None
    if core.graph_cigar_query_span(fields[6], fields[0]) != (
        query_coord.end - query_coord.start
    ):
        raise ValueError(
            f"{source}:{line_number}: graph CIGAR query span disagrees with "
            f"{fields[2]!r}"
        )
    from alternative_intervals import read_tag
    return AnnotationRow(
        line_number=line_number,
        query_name=fields[0],
        reference_name=fields[1],
        query_coord=query_coord,
        query_alignment=fields[4],
        reference_alignment=fields[5],
        graph_cigar=fields[6],
        reference_coord=core.parse_coord(fields[3]),
        alternative_intervals=read_tag(fields),
    )


def annotation_from_tsv(text: str, line_number: int) -> AnnotationRow:
    row = annotation_from_fields(
        text.rstrip("\r\n").split("\t"), line_number, "generated gap row",
    )
    if row is None:
        raise ValueError("generated gap row has a zero-width query interval")
    return row


def annotation_to_tsv(row: AnnotationRow) -> str:
    reference_coord = row.reference_coord
    if reference_coord is None:
        reference_coord = _primary_reference_coord(
            row.graph_cigar, row.query_name,
        )
    return "\t".join((
        row.query_name,
        row.reference_name,
        coord_text(row.query_coord),
        coord_text(reference_coord),
        row.query_alignment,
        row.reference_alignment,
        row.graph_cigar,
    ))


def infer_query_genome(
    annotations: Sequence[AnnotationRow],
    lift_by_allele: Mapping[str, object],
) -> str:
    genomes = set()
    for annotation in annotations:
        for allele in query_name_parts(annotation.query_name):
            if allele not in lift_by_allele:
                continue
            try:
                genomes.add(parse_allele_name(allele)[4])
            except ValueError:
                continue
    if len(genomes) != 1:
        raise ValueError(
            "could not infer exactly one sample/haplotype from the input and "
            f"GenomeLift rows; found {sorted(genomes)!r}; provide --query-genome"
        )
    return next(iter(genomes))


def _clusters(
    annotations: Sequence[AnnotationRow],
) -> List[Tuple[int, int, Tuple[AnnotationRow, ...]]]:
    ordered = sorted(
        annotations,
        key=lambda row: (
            row.query_coord.start,
            row.query_coord.end,
            row.line_number,
        ),
    )
    output: List[Tuple[int, int, Tuple[AnnotationRow, ...]]] = []
    start = ordered[0].query_coord.start
    end = ordered[0].query_coord.end
    members = [ordered[0]]
    for row in ordered[1:]:
        if row.query_coord.start <= end:
            end = max(end, row.query_coord.end)
            members.append(row)
            continue
        output.append((start, end, tuple(members)))
        start = row.query_coord.start
        end = row.query_coord.end
        members = [row]
    output.append((start, end, tuple(members)))
    return output


def _extension_parts(token: str) -> Tuple[str, str, str]:
    token = (token or "").strip()
    neighbor = ""
    value = token
    if ":" in token:
        neighbor, value = token.rsplit(":", 1)
        neighbor = neighbor.strip()
        value = value.strip()
    if value in {"-inf", "+inf", "inf"}:
        return neighbor, "inf", value
    sign = value[0] if value[:1] in {"+", "-"} else ""
    body = value[1:] if sign else value
    if value in {"", ".", "NA", "None", "none", "null"}:
        return neighbor, "", value
    if not body.isdigit():
        raise ValueError(f"malformed GenomeLift extension token {token!r}")
    return neighbor, sign, value


def _self_priority(token: str) -> int:
    _neighbor, sign, _value = _extension_parts(token)
    if sign == "+":
        return 2
    if sign == "-":
        return -1
    if sign == "inf":
        return -10
    # An unsigned gap says that the adjacent block has first priority.
    return 0


def _priority_blocks(
    annotations: Iterable[AnnotationRow],
    lift_by_allele: Mapping[str, object],
    *,
    facing_side: str,
    boundary_position: int,
) -> List[PriorityBlock]:
    output: List[PriorityBlock] = []
    for annotation in annotations:
        for allele in query_name_parts(annotation.query_name):
            lift = lift_by_allele.get(allele)
            if lift is None or lift.locus.chrom != annotation.query_coord.chrom:
                continue
            token = (
                lift.right_extension
                if facing_side == "right"
                else lift.left_extension
            )
            if _extension_parts(token)[1] == "inf":
                continue
            output.append(PriorityBlock(
                annotation=annotation,
                allele=allele,
                locus=lift.locus,
                facing_token=token,
            ))
    output.sort(key=lambda block: (
        abs(
            (block.locus.end if facing_side == "right" else block.locus.start)
            - boundary_position
        ),
        block.annotation.line_number,
        block.allele,
    ))
    return output


def choose_neighboring_blocks(
    left_annotations: Sequence[AnnotationRow],
    right_annotations: Sequence[AnnotationRow],
    lift_by_allele: Mapping[str, object],
    gap_start: int,
    gap_end: int,
) -> Tuple[PriorityBlock, PriorityBlock, str]:
    left_blocks = _priority_blocks(
        left_annotations,
        lift_by_allele,
        facing_side="right",
        boundary_position=gap_start,
    )
    right_blocks = _priority_blocks(
        right_annotations,
        lift_by_allele,
        facing_side="left",
        boundary_position=gap_end,
    )
    if not left_blocks or not right_blocks:
        raise ValueError("one or both query-gap boundaries lack a GenomeLift block")

    left_neighbors = {
        id(block): _extension_parts(block.facing_token)[0]
        for block in left_blocks
    }
    right_neighbors = {
        id(block): _extension_parts(block.facing_token)[0]
        for block in right_blocks
    }

    def left_key(block: PriorityBlock) -> Tuple[int, int, str]:
        return (
            abs(block.locus.end - gap_start),
            block.annotation.line_number,
            block.allele,
        )

    def right_key(block: PriorityBlock) -> Tuple[int, int, str]:
        return (
            abs(block.locus.start - gap_end),
            block.annotation.line_number,
            block.allele,
        )

    best_left_by_allele: Dict[str, PriorityBlock] = {}
    best_left_by_allele_neighbor: Dict[
        Tuple[str, str], PriorityBlock
    ] = {}
    for block in left_blocks:
        previous = best_left_by_allele.get(block.allele)
        if previous is None or left_key(block) < left_key(previous):
            best_left_by_allele[block.allele] = block
        key = (block.allele, left_neighbors[id(block)])
        previous = best_left_by_allele_neighbor.get(key)
        if previous is None or left_key(block) < left_key(previous):
            best_left_by_allele_neighbor[key] = block

    best_right_by_allele: Dict[str, PriorityBlock] = {}
    best_right_by_allele_neighbor: Dict[
        Tuple[str, str], PriorityBlock
    ] = {}
    for block in right_blocks:
        previous = best_right_by_allele.get(block.allele)
        if previous is None or right_key(block) < right_key(previous):
            best_right_by_allele[block.allele] = block
        key = (block.allele, right_neighbors[id(block)])
        previous = best_right_by_allele_neighbor.get(key)
        if previous is None or right_key(block) < right_key(previous):
            best_right_by_allele_neighbor[key] = block

    def pair_key(
        left: PriorityBlock, right: PriorityBlock,
    ) -> Tuple[int, int, int, int, str, str]:
        reciprocal = int(
            left_neighbors[id(left)] == right.allele
        ) + int(right_neighbors[id(right)] == left.allele)
        return (
            -reciprocal,
            abs(left.locus.end - gap_start)
            + abs(right.locus.start - gap_end),
            left.annotation.line_number,
            right.annotation.line_number,
            left.allele,
            right.allele,
        )

    # The old implementation materialized the complete left x right product.
    # Only linked pairs can have reciprocal score 1/2; the best unlinked pair
    # is the independently closest block on each side. Enumerating these
    # O(L+R) candidates preserves the exact old comparison key.
    candidate_pairs: List[Tuple[PriorityBlock, PriorityBlock]] = [(
        min(left_blocks, key=left_key),
        min(right_blocks, key=right_key),
    )]
    for left_block in left_blocks:
        wanted = left_neighbors[id(left_block)]
        if not wanted:
            continue
        right_block = best_right_by_allele.get(wanted)
        if right_block is not None:
            candidate_pairs.append((left_block, right_block))
        reciprocal_right = best_right_by_allele_neighbor.get((
            wanted, left_block.allele,
        ))
        if reciprocal_right is not None:
            candidate_pairs.append((left_block, reciprocal_right))
    for right_block in right_blocks:
        wanted = right_neighbors[id(right_block)]
        if not wanted:
            continue
        left_block = best_left_by_allele.get(wanted)
        if left_block is not None:
            candidate_pairs.append((left_block, right_block))
        reciprocal_left = best_left_by_allele_neighbor.get((
            wanted, right_block.allele,
        ))
        if reciprocal_left is not None:
            candidate_pairs.append((reciprocal_left, right_block))

    left, right = min(candidate_pairs, key=lambda pair: pair_key(*pair))
    left_priority = _self_priority(left.facing_token)
    right_priority = _self_priority(right.facing_token)
    selected_side = "right" if right_priority > left_priority else "left"
    return left, right, selected_side


def _pair_priority_side(
    left: AnnotationRow,
    right: AnnotationRow,
    lift_by_allele: Mapping[str, object],
) -> Optional[str]:
    """Return the higher-priority side for two physical query neighbors."""
    try:
        _left, _right, selected_side = choose_neighboring_blocks(
            (left,),
            (right,),
            lift_by_allele,
            left.query_coord.end,
            right.query_coord.start,
        )
    except (KeyError, ValueError):
        return None
    return selected_side


def reconcile_query_overlaps(
    annotations: Sequence[AnnotationRow],
    lift_by_allele: Mapping[str, object],
) -> Tuple[List[AnnotationRow], Dict[str, int]]:
    """Assign overlapping query bases to the higher-priority block.

    Rows are swept in physical query order.  When the current row wins, the
    lower-priority suffix already on the stack is trimmed and the comparison
    is repeated against the preceding row.  This keeps the work proportional
    to the number of actual overlaps rather than comparing every row pair.
    """
    by_contig: Dict[str, List[AnnotationRow]] = {}
    for row in annotations:
        by_contig.setdefault(row.query_coord.chrom, []).append(row)
    output: List[AnnotationRow] = []
    stats = {
        "query_overlaps": 0,
        "query_overlap_bases": 0,
        "query_overlap_rows_removed": 0,
        "query_overlap_missing_priority": 0,
    }
    for contig in sorted(by_contig):
        ordered = sorted(by_contig[contig], key=lambda row: (
            row.query_coord.start,
            row.query_coord.end,
            row.line_number,
        ))
        stack: List[AnnotationRow] = []
        for original in ordered:
            current: Optional[AnnotationRow] = original
            while current is not None and stack:
                left = stack[-1]
                overlap_start = max(
                    left.query_coord.start, current.query_coord.start,
                )
                overlap_end = min(
                    left.query_coord.end, current.query_coord.end,
                )
                if overlap_end <= overlap_start:
                    break
                selected_side = _pair_priority_side(
                    left, current, lift_by_allele,
                )
                if selected_side is None:
                    stats["query_overlap_missing_priority"] += 1
                    break
                stats["query_overlaps"] += 1
                stats["query_overlap_bases"] += overlap_end - overlap_start
                if selected_side == "left":
                    # The right/lower-priority block loses its physical prefix.
                    new_start = max(current.query_coord.start, overlap_end)
                    if new_start >= current.query_coord.end:
                        current = None
                        stats["query_overlap_rows_removed"] += 1
                    else:
                        current = slice_annotation_query_interval(
                            current, new_start, current.query_coord.end,
                        )
                    break

                # The current/right block wins. Remove the lower-priority
                # suffix from the previous row, then compare against an older
                # active row in case the overlap chain continues.
                new_end = min(left.query_coord.end, overlap_start)
                stack.pop()
                if new_end > left.query_coord.start:
                    stack.append(slice_annotation_query_interval(
                        left, left.query_coord.start, new_end,
                    ))
                else:
                    stats["query_overlap_rows_removed"] += 1
            if current is not None:
                stack.append(current)
        output.extend(stack)
    output.sort(key=lambda row: (
        row.query_coord.chrom,
        row.query_coord.start,
        row.query_coord.end,
        row.line_number,
    ))
    return output, stats


def _signed_reference_gap(
    left: AnnotationRow,
    right: AnnotationRow,
    reference_index: Mapping[str, Tuple[int, int, int, int]],
) -> Optional[Tuple[int, ReferenceBoundary, ReferenceBoundary]]:
    try:
        left_boundary = reference_boundary(left, "right")
        right_boundary = reference_boundary(right, "left")
    except ValueError:
        return None
    if (
        not _same_reference(left_boundary, right_boundary, reference_index)
        or left_boundary.strand != right_boundary.strand
    ):
        return None
    if left_boundary.strand == "+":
        gap = right_boundary.position - left_boundary.position
    else:
        gap = left_boundary.position - right_boundary.position
    return gap, left_boundary, right_boundary


def reconcile_reference_overlaps(
    annotations: Sequence[AnnotationRow],
    lift_by_allele: Mapping[str, object],
    query_reader: core.FastaRegionReader,
    reference_index: Mapping[str, Tuple[int, int, int, int]],
) -> Tuple[List[AnnotationRow], Dict[str, int]]:
    """Convert lower-priority same-reference overlap into query insertion."""
    by_contig: Dict[str, List[AnnotationRow]] = {}
    for row in annotations:
        by_contig.setdefault(row.query_coord.chrom, []).append(row)
    output: List[AnnotationRow] = []
    stats = {
        "reference_overlaps": 0,
        "reference_overlap_bases": 0,
        "reference_overlap_missing_priority": 0,
        "reference_overlap_unprojectable": 0,
    }
    for contig in sorted(by_contig):
        ordered = sorted(by_contig[contig], key=lambda row: (
            row.query_coord.start,
            row.query_coord.end,
            row.line_number,
        ))
        for index in range(len(ordered) - 1):
            left = ordered[index]
            right = ordered[index + 1]
            relation = _signed_reference_gap(
                left, right, reference_index,
            )
            if relation is None or relation[0] >= 0:
                continue
            selected_side = _pair_priority_side(
                left, right, lift_by_allele,
            )
            if selected_side is None:
                stats["reference_overlap_missing_priority"] += 1
                continue
            gap, left_boundary, right_boundary = relation
            try:
                if selected_side == "left":
                    right = convert_reference_overlap_to_query_insertion(
                        right,
                        left_boundary.path,
                        left_boundary.position,
                        "left",
                        query_reader,
                    )
                    ordered[index + 1] = right
                else:
                    left = convert_reference_overlap_to_query_insertion(
                        left,
                        right_boundary.path,
                        right_boundary.position,
                        "right",
                        query_reader,
                    )
                    ordered[index] = left
            except (KeyError, ValueError):
                stats["reference_overlap_unprojectable"] += 1
                continue
            stats["reference_overlaps"] += 1
            stats["reference_overlap_bases"] += -gap
        output.extend(ordered)
    output.sort(key=lambda row: (
        row.query_coord.chrom,
        row.query_coord.start,
        row.query_coord.end,
        row.line_number,
    ))
    return output, stats


def _segment_query_endpoints(
    graph_cigar: str, row_name: str,
) -> Tuple[Tuple[str, int, str], Tuple[str, int, str]]:
    first: Optional[Tuple[str, int, str]] = None
    last: Optional[Tuple[str, int, str]] = None
    for segment in core.parse_graphic_segments(graph_cigar, row_name):
        rpos = segment.start if segment.direction == ">" else segment.end
        for operation in segment.ops:
            if operation.op == "H":
                continue
            reference_size = core.ref_consume_pair(operation.op, operation.n)
            query_size = core.query_consume(operation.op, operation.n)
            before = rpos
            after = (
                rpos - reference_size
                if segment.direction == "<"
                else rpos + reference_size
            )
            if query_size:
                if first is None:
                    first = (segment.path, before, segment.direction)
                last = (segment.path, after, segment.direction)
            rpos = after
    if first is None or last is None:
        raise ValueError(f"{row_name}: graph CIGAR has no query-consuming operation")
    return first, last


def reference_boundary(
    annotation: AnnotationRow, physical_side: str,
) -> ReferenceBoundary:
    query_start, query_end = _segment_query_endpoints(
        annotation.graph_cigar, annotation.query_name,
    )
    if annotation.query_coord.strand == "+":
        endpoint = query_start if physical_side == "left" else query_end
        physical_strand = "+" if endpoint[2] == ">" else "-"
    else:
        endpoint = query_end if physical_side == "left" else query_start
        physical_strand = "-" if endpoint[2] == ">" else "+"
    return ReferenceBoundary(endpoint[0], endpoint[1], physical_strand)


def project_reference_position_to_query(
    annotation: AnnotationRow,
    reference_path: str,
    reference_position: int,
    physical_side: str,
) -> int:
    """Project one reference breakpoint into a neighboring query block.

    ``physical_side`` selects the occurrence closest to that assembly edge.
    This matters for repeated/multi-segment graph CIGARs.  Insertions at the
    requested reference coordinate contribute both query-axis boundaries as
    candidates; the edge rule deterministically selects the appropriate one.
    """
    wanted_path = core.path_key(reference_path)
    query_position = 0
    candidates: List[int] = []
    for segment in core.parse_graphic_segments(
        annotation.graph_cigar, annotation.query_name,
    ):
        rpos = segment.start if segment.direction == ">" else segment.end
        same_path = core.path_key(segment.path) == wanted_path
        for operation in segment.ops:
            if operation.op == "H":
                continue
            reference_size = core.ref_consume_pair(operation.op, operation.n)
            query_size = core.query_consume(operation.op, operation.n)
            before = rpos
            after = (
                rpos - reference_size
                if segment.direction == "<"
                else rpos + reference_size
            )
            if same_path:
                if reference_size and min(before, after) <= reference_position <= max(
                    before, after,
                ):
                    if query_size:
                        offset = (
                            reference_position - before
                            if segment.direction == ">"
                            else before - reference_position
                        )
                        candidates.append(query_position + offset)
                    else:
                        candidates.append(query_position)
                elif not reference_size and reference_position == before:
                    candidates.extend((
                        query_position,
                        query_position + query_size,
                    ))
            query_position += query_size
            rpos = after
    if not candidates:
        raise ValueError(
            f"{annotation.query_name}: reference breakpoint "
            f"{reference_path}:{reference_position} is not represented by the "
            "neighboring graph CIGAR"
        )

    query_span = annotation.query_coord.end - annotation.query_coord.start
    physical_candidates = []
    for offset in candidates:
        if offset < 0 or offset > query_span:
            continue
        if annotation.query_coord.strand == "+":
            physical = annotation.query_coord.start + offset
        else:
            physical = annotation.query_coord.end - offset
        if annotation.query_coord.start <= physical <= annotation.query_coord.end:
            physical_candidates.append(physical)
    if not physical_candidates:
        raise ValueError(
            f"{annotation.query_name}: projected reference breakpoint lies "
            "outside the query annotation"
        )
    edge = (
        annotation.query_coord.start
        if physical_side == "left"
        else annotation.query_coord.end
    )
    return min(physical_candidates, key=lambda value: (abs(value - edge), value))


def _resolve_reference_path(
    path: str, index: Mapping[str, Tuple[int, int, int, int]],
) -> Optional[str]:
    for candidate in (
        path,
        core.strip_match_name(path),
        core.path_key(path),
    ):
        if candidate in index:
            return candidate
    return None


def annotation_targets_non_primary(
    annotation: AnnotationRow,
    primary_reference_index: Optional[
        Mapping[str, Tuple[int, int, int, int]]
    ],
) -> bool:
    """Return whether an alignment is on an alternative/novel target.

    Column 4 is authoritative. A CIGAR-only legacy row falls back to its
    reference-consuming graph paths. A row spanning any primary path is not
    penalized merely because it also encodes a promoted nested insertion.
    """
    if primary_reference_index is None:
        return False
    if annotation.reference_coord is not None:
        return _resolve_reference_path(
            annotation.reference_coord.chrom, primary_reference_index,
        ) is None
    paths = [
        segment.path
        for segment in core.parse_graphic_segments(
            annotation.graph_cigar, annotation.query_name,
        )
        if segment.end > segment.start
    ]
    return bool(paths) and all(
        _resolve_reference_path(path, primary_reference_index) is None
        for path in paths
    )


def ownership_priority_score(
    annotation: AnnotationRow,
    score: int,
    primary_reference_index: Optional[
        Mapping[str, Tuple[int, int, int, int]]
    ],
) -> int:
    """Apply the half-priority rule to alternative/novel alignments."""
    if annotation_targets_non_primary(
        annotation, primary_reference_index,
    ):
        return score // 2
    return score


def _same_reference(
    left: ReferenceBoundary,
    right: ReferenceBoundary,
    reference_index: Mapping[str, Tuple[int, int, int, int]],
) -> bool:
    left_path = _resolve_reference_path(left.path, reference_index)
    right_path = _resolve_reference_path(right.path, reference_index)
    return left_path is not None and left_path == right_path


def gap_is_eligible(
    gap_size: int,
    unmasked: int,
    maximum_gap: int = DEFAULT_MAXIMUM_GAP,
    minimum_unmasked: int = DEFAULT_MINIMUM_UNMASKED,
    minimum_fraction: float = DEFAULT_MINIMUM_UNMASKED_FRACTION,
    short_gap: int = DEFAULT_ANCHOR,
) -> bool:
    if gap_size <= 0 or gap_size >= maximum_gap:
        return False
    if gap_size < short_gap:
        return True
    return (
        unmasked > minimum_unmasked
        and (unmasked / gap_size) > minimum_fraction
    )


def discover_gap_tasks(
    annotations: Sequence[AnnotationRow],
    lift_by_allele: Mapping[str, object],
    query_reader: core.FastaRegionReader,
    reference_index: Mapping[str, Tuple[int, int, int, int]],
    query_genome: str,
    maximum_gap: int = DEFAULT_MAXIMUM_GAP,
    minimum_unmasked: int = DEFAULT_MINIMUM_UNMASKED,
    minimum_fraction: float = DEFAULT_MINIMUM_UNMASKED_FRACTION,
) -> Tuple[List[GapTask], Dict[str, int]]:
    by_contig: Dict[str, List[AnnotationRow]] = {}
    for row in annotations:
        by_contig.setdefault(row.query_coord.chrom, []).append(row)

    tasks: List[GapTask] = []
    stats = {
        "internal_gaps": 0,
        "eligible_gaps": 0,
        "masked_or_large": 0,
        "missing_priority": 0,
        "unresolvable_reference_boundary": 0,
        "non_primary_reference_boundary": 0,
        "zero_reference": 0,
        "query_insertion": 0,
        "unprojectable_tandem": 0,
        "tandem_duplication": 0,
        "bounded_reference": 0,
        "one_sided_reference": 0,
    }
    for contig in sorted(by_contig):
        if contig not in query_reader.index:
            raise KeyError(f"query contig {contig!r} is absent from the FASTA index")
        contig_clusters = _clusters(by_contig[contig])
        for left_cluster, right_cluster in zip(
            contig_clusters, contig_clusters[1:],
        ):
            _left_start, gap_start, left_members = left_cluster
            gap_end, _right_end, right_members = right_cluster
            if gap_end <= gap_start:
                continue
            stats["internal_gaps"] += 1
            gap_size = gap_end - gap_start
            sequence = query_reader.fetch(contig, gap_start, gap_end, "+")
            unmasked = count_unmasked(sequence)
            if not gap_is_eligible(
                gap_size,
                unmasked,
                maximum_gap,
                minimum_unmasked,
                minimum_fraction,
                # The masking exception is fixed at <100 bp.  Changing the
                # alignment-anchor length must not change gap eligibility.
                DEFAULT_ANCHOR,
            ):
                stats["masked_or_large"] += 1
                continue
            left_edge_rows = tuple(
                row for row in left_members
                if row.query_coord.end == gap_start
            )
            right_edge_rows = tuple(
                row for row in right_members
                if row.query_coord.start == gap_end
            )
            try:
                left_block, right_block, selected_side = choose_neighboring_blocks(
                    left_edge_rows,
                    right_edge_rows,
                    lift_by_allele,
                    gap_start,
                    gap_end,
                )
                left_boundary = reference_boundary(
                    left_block.annotation, "right",
                )
                right_boundary = reference_boundary(
                    right_block.annotation, "left",
                )
            except (KeyError, ValueError):
                stats["missing_priority"] += 1
                continue

            # Alternative and novel loci are immutable templates at this
            # stage. Do not fill or extend across a gap touching one, even if
            # the other (selected) side is primary-reference backed.
            if (
                left_block.annotation.alternative_intervals is not None
                or right_block.annotation.alternative_intervals is not None
                or _resolve_reference_path(
                    left_boundary.path, reference_index,
                ) is None
                or _resolve_reference_path(
                    right_boundary.path, reference_index,
                ) is None
            ):
                stats["non_primary_reference_boundary"] += 1
                continue

            # A neighboring row can carry a self-referential graph CIGAR (a
            # reference-free/novel block whose path is its own assembly
            # contig).  Such a boundary cannot anchor a reference window;
            # skip the gap instead of failing the task at alignment time.
            selected_boundary_path = (
                right_boundary.path
                if selected_side == "right"
                else left_boundary.path
            )
            if _resolve_reference_path(
                selected_boundary_path, reference_index,
            ) is None:
                stats["unresolvable_reference_boundary"] += 1
                continue

            selected = (
                right_block.annotation
                if selected_side == "right"
                else left_block.annotation
            )
            same_reference = _same_reference(
                left_boundary, right_boundary, reference_index,
            )
            reference_gap_size: Optional[int] = None
            insertion_start: Optional[int] = None
            insertion_end: Optional[int] = None
            mode = ""
            if same_reference:
                selected_boundary = (
                    right_boundary if selected_side == "right" else left_boundary
                )
                reference_gap_size = (
                    right_boundary.position - left_boundary.position
                    if selected_boundary.strand == "+"
                    else left_boundary.position - right_boundary.position
                )
                if reference_gap_size == 0:
                    # The query interval is absent from both neighboring rows
                    # but the two reference breakpoints touch. Preserve it as
                    # a pure insertion instead of silently dropping the gap.
                    insertion_start = gap_start
                    insertion_end = gap_end
                    mode = "tandem-duplication"
                    stats["zero_reference"] += 1
                    stats["query_insertion"] += 1
                elif reference_gap_size < 0:
                    try:
                        if selected_side == "left":
                            insertion_start = gap_start
                            insertion_end = project_reference_position_to_query(
                                right_block.annotation,
                                left_boundary.path,
                                left_boundary.position,
                                "left",
                            )
                            if insertion_end <= gap_end:
                                raise ValueError(
                                    "right block does not contribute query "
                                    "sequence across the reference overlap"
                                )
                        else:
                            insertion_start = project_reference_position_to_query(
                                left_block.annotation,
                                right_boundary.path,
                                right_boundary.position,
                                "right",
                            )
                            insertion_end = gap_end
                            if insertion_start >= gap_start:
                                raise ValueError(
                                    "left block does not contribute query "
                                    "sequence across the reference overlap"
                                )
                    except ValueError:
                        stats["unprojectable_tandem"] += 1
                        continue
                    mode = "tandem-duplication"
                    stats["tandem_duplication"] += 1

            # A bounded alignment must include the *complete* interval between
            # the two reference breakpoints.  The former 1-kb cap sent every
            # larger same-reference gap through the one-sided 10-kb rescue
            # window.  For an unmodified reference sample, a 23,667-bp gap was
            # consequently encoded as 13,667I + 10,000=; deletion discovery
            # then emitted the missing 13,667-bp reference prefix as a false
            # DEL.  Keep the existing three-times-query plausibility rule, but
            # allow it across the same <maximum_gap range already accepted for
            # query gaps.
            bounded_limit = min(maximum_gap, 3 * gap_size)
            if mode == "tandem-duplication":
                pass
            elif (
                same_reference
                and reference_gap_size is not None
                and 0 < reference_gap_size < bounded_limit
            ):
                mode = "bounded"
                stats["bounded_reference"] += 1
            else:
                mode = "one-sided"
                stats["one_sided_reference"] += 1

            index = len(tasks) + 1
            tasks.append(GapTask(
                index=index,
                name=f"gap_{query_genome}_{index}",
                contig=contig,
                start=gap_start,
                end=gap_end,
                left=left_block.annotation,
                right=right_block.annotation,
                selected=selected,
                selected_side=selected_side,
                left_boundary=left_boundary,
                right_boundary=right_boundary,
                mode=mode,
                reference_gap_size=reference_gap_size,
                insertion_start=insertion_start,
                insertion_end=insertion_end,
            ))
            stats["eligible_gaps"] += 1
    return tasks, stats


def _reference_window(
    task: GapTask,
    reference_reader: core.FastaRegionReader,
    anchor: int,
    maximum_reference_extension: int,
) -> Tuple[str, int, int, str]:
    if task.mode == "bounded":
        path = _resolve_reference_path(
            task.left_boundary.path, reference_reader.index,
        )
        if path is None:
            raise KeyError(
                f"reference path {task.left_boundary.path!r} is absent from the FAI"
            )
        # The signed-positive check performed during discovery guarantees that
        # traversal from the left query block reaches the right block in the
        # higher-priority block's reference direction.
        selected_boundary = (
            task.right_boundary
            if task.selected_side == "right"
            else task.left_boundary
        )
        strand = selected_boundary.strand
        if strand == "+":
            start = task.left_boundary.position - anchor
            end = task.right_boundary.position + anchor
        else:
            start = task.right_boundary.position - anchor
            end = task.left_boundary.position + anchor
    else:
        boundary = (
            task.left_boundary
            if task.selected_side == "left"
            else task.right_boundary
        )
        path = _resolve_reference_path(boundary.path, reference_reader.index)
        if path is None:
            raise KeyError(f"reference path {boundary.path!r} is absent from the FAI")
        strand = boundary.strand
        gap_size = task.end - task.start
        extension = (
            min(maximum_reference_extension, 3 * gap_size)
            if gap_size < 100
            else maximum_reference_extension
        )
        position = boundary.position
        if task.selected_side == "left":
            if strand == "+":
                start, end = position, position + extension
            else:
                start, end = position - extension, position
        else:
            if strand == "+":
                start, end = position - extension, position
            else:
                start, end = position, position + extension

    contig_size = reference_reader.index[path][0]
    start = max(0, min(start, contig_size))
    end = max(start, min(end, contig_size))
    if end <= start:
        raise ValueError(f"{task.name}: empty reference alignment window")
    return path, start, end, strand


def _query_window(
    task: GapTask,
    query_reader: core.FastaRegionReader,
    anchor: int,
) -> Tuple[int, int, int, int]:
    contig_size = query_reader.index[task.contig][0]
    if task.mode == "bounded":
        start = max(0, task.start - anchor)
        end = min(contig_size, task.end + anchor)
    else:
        start = task.start
        end = task.end
    left_anchor = task.start - start
    right_anchor = end - task.end
    return start, end, left_anchor, right_anchor


def _reference_coord_from_cigar(graph_cigar: str, row_name: str) -> core.Coord:
    segments = core.parse_graphic_segments(graph_cigar, row_name)
    if not segments:
        raise ValueError(f"{row_name}: empty gap graph CIGAR")
    paths = {segment.path for segment in segments}
    directions = {segment.direction for segment in segments}
    if len(paths) != 1 or len(directions) != 1:
        raise ValueError(f"{row_name}: gap graph CIGAR is not one reference interval")
    start = min(segment.start for segment in segments)
    end = max(segment.end for segment in segments)
    direction = next(iter(directions))
    return core.Coord(next(iter(paths)), start, end, "+" if direction == ">" else "-")


def _primary_reference_coord(graph_cigar: str, row_name: str) -> core.Coord:
    """Return the main (first named) reference interval in one graph CIGAR."""
    segments = core.parse_graphic_segments(graph_cigar, row_name)
    if not segments:
        raise ValueError(f"{row_name}: empty graph CIGAR")
    segment = segments[0]
    return core.Coord(
        segment.path,
        min(segment.start, segment.end),
        max(segment.start, segment.end),
        "+" if segment.direction == ">" else "-",
    )


def _reference_alignment_text(graph_cigar: str, row_name: str) -> str:
    """Summarize all reference intervals represented by a graph CIGAR."""
    grouped: Dict[Tuple[str, str], List[Tuple[int, int]]] = {}
    for segment in core.parse_graphic_segments(graph_cigar, row_name):
        if segment.end <= segment.start:
            continue
        strand = "+" if segment.direction == ">" else "-"
        grouped.setdefault((segment.path, strand), []).append((
            min(segment.start, segment.end),
            max(segment.start, segment.end),
        ))
    output: List[str] = []
    for (path, strand), intervals in sorted(grouped.items()):
        intervals.sort()
        start, end = intervals[0]
        for next_start, next_end in intervals[1:]:
            if next_start <= end:
                end = max(end, next_end)
            else:
                output.append(f"{path}:{start}-{end}{strand}")
                start, end = next_start, next_end
        output.append(f"{path}:{start}-{end}{strand}")
    return ";".join(output)


def _refresh_annotation(
    row: AnnotationRow,
    graph_cigar: str,
    query_coord: Optional[core.Coord] = None,
) -> AnnotationRow:
    query_coord = query_coord or row.query_coord
    observed = core.graph_cigar_query_span(graph_cigar, row.query_name)
    expected = query_coord.end - query_coord.start
    if observed != expected:
        raise ValueError(
            f"{row.query_name}: reconciled graph CIGAR spans {observed} query "
            f"bases, expected {expected} from {coord_text(query_coord)}"
        )
    return dataclasses.replace(
        row,
        query_coord=query_coord,
        graph_cigar=graph_cigar,
        reference_coord=_primary_reference_coord(
            graph_cigar, row.query_name,
        ),
        reference_alignment=_reference_alignment_text(
            graph_cigar, row.query_name,
        ),
    )


def _query_axis_slice(
    coord: core.Coord, physical_start: int, physical_end: int,
) -> Tuple[int, int]:
    physical_start = max(coord.start, min(physical_start, coord.end))
    physical_end = max(physical_start, min(physical_end, coord.end))
    if coord.strand == "+":
        return physical_start - coord.start, physical_end - coord.start
    return coord.end - physical_end, coord.end - physical_start


def slice_annotation_query_interval(
    row: AnnotationRow, physical_start: int, physical_end: int,
) -> AnnotationRow:
    """Trim one row on the physical assembly axis and reproject reference."""
    if physical_end <= physical_start:
        raise ValueError(f"{row.query_name}: query reconciliation removed the row")
    qstart, qend = _query_axis_slice(
        row.query_coord, physical_start, physical_end,
    )
    graph_cigar = core.slice_graph_cigar_by_query(
        row.graph_cigar,
        qstart,
        qend,
        row.query_name,
        include_left_boundary_deletions=False,
    )
    return _refresh_annotation(
        row,
        graph_cigar,
        core.Coord(
            row.query_coord.chrom,
            physical_start,
            physical_end,
            row.query_coord.strand,
        ),
    )


def slice_annotation_query_intervals(
    row: AnnotationRow,
    physical_intervals: Sequence[Tuple[int, int]],
) -> Dict[Tuple[int, int], AnnotationRow]:
    """Materialize many disjoint slices without K full CIGAR rescans.

    A balanced recursive partition makes each graph-CIGAR operation appear in
    at most O(log K) intermediate slices. The former per-label loop scanned
    the complete row K times, which became quadratic when both a merged row's
    PA count and its CIGAR operation count grew together.
    """
    intervals = sorted(set(
        (max(row.query_coord.start, int(start)),
         min(row.query_coord.end, int(end)))
        for start, end in physical_intervals
        if int(end) > int(start)
    ))
    if not intervals:
        return {}
    for index, (start, end) in enumerate(intervals):
        if end <= start:
            raise ValueError(f"{row.query_name}: empty ownership slice")
        if index and intervals[index - 1][1] > start:
            raise ValueError(
                f"{row.query_name}: finalized PA ownership slices overlap: "
                f"{intervals[index - 1]} and {(start, end)}"
            )

    output: Dict[Tuple[int, int], AnnotationRow] = {}

    def materialize(
        parent: AnnotationRow,
        wanted: Sequence[Tuple[int, int]],
    ) -> None:
        if len(wanted) == 1:
            interval = wanted[0]
            if interval == (
                parent.query_coord.start, parent.query_coord.end,
            ):
                output[interval] = parent
            else:
                output[interval] = slice_annotation_query_interval(
                    parent, interval[0], interval[1],
                )
            return
        middle = len(wanted) // 2
        left_wanted = wanted[:middle]
        right_wanted = wanted[middle:]
        for subset in (left_wanted, right_wanted):
            subset_start = subset[0][0]
            subset_end = subset[-1][1]
            if (
                subset_start == parent.query_coord.start
                and subset_end == parent.query_coord.end
            ):
                child = parent
            else:
                child = slice_annotation_query_interval(
                    parent, subset_start, subset_end,
                )
            materialize(child, subset)

    materialize(row, intervals)
    return output


def _insert_query_payload_at_edge(
    graph_cigar: str, payload: str, *, at_start: bool,
) -> str:
    if not payload:
        return graph_cigar
    insertion = f"{len(payload)}I{payload}"
    markers = list(re.finditer(r"([<>])([^:<>]+):", graph_cigar))
    if not markers:
        raise ValueError("cannot add an insertion to an empty graph CIGAR")
    if at_start:
        body_start = markers[0].end()
        leading_h = re.match(r"\d+H", graph_cigar[body_start:])
        position = body_start + (leading_h.end() if leading_h else 0)
    else:
        trailing_h = re.search(r"\d+H$", graph_cigar)
        position = trailing_h.start() if trailing_h else len(graph_cigar)
    return graph_cigar[:position] + insertion + graph_cigar[position:]


def convert_reference_overlap_to_query_insertion(
    row: AnnotationRow,
    reference_path: str,
    reference_position: int,
    physical_side: str,
    query_reader: core.FastaRegionReader,
) -> AnnotationRow:
    """Remove a lower-priority reference overlap but retain its query bases.

    The reference breakpoint is projected into the row, the non-owning query
    edge is removed from the aligned CIGAR, and those exact query bases are
    reinserted at the same query-side edge as a pure ``I`` operation.
    """
    split = project_reference_position_to_query(
        row, reference_path, reference_position, physical_side,
    )
    if physical_side == "left":
        insertion_start, insertion_end = row.query_coord.start, split
        kept_start, kept_end = split, row.query_coord.end
    else:
        insertion_start, insertion_end = split, row.query_coord.end
        kept_start, kept_end = row.query_coord.start, split
    if insertion_end <= insertion_start or kept_end <= kept_start:
        raise ValueError(
            f"{row.query_name}: reference-overlap projection produced an "
            "empty insertion or retained alignment"
        )

    qstart, qend = _query_axis_slice(row.query_coord, kept_start, kept_end)
    retained = core.slice_graph_cigar_by_query(
        row.graph_cigar,
        qstart,
        qend,
        row.query_name,
        include_left_boundary_deletions=False,
    )
    payload = query_reader.fetch(
        row.query_coord.chrom, insertion_start, insertion_end, "+",
    )
    # CIGAR query order follows query_coord.strand, not forward assembly order.
    if row.query_coord.strand == "-":
        payload = core.revcomp(payload)
    insertion_at_start = (
        physical_side == "left"
        if row.query_coord.strand == "+"
        else physical_side == "right"
    )
    reconciled = _insert_query_payload_at_edge(
        retained, payload, at_start=insertion_at_start,
    )
    return _refresh_annotation(row, reconciled)


def _merge_intervals(
    intervals: Iterable[Tuple[int, int]],
) -> List[Tuple[int, int]]:
    ordered = sorted((start, end) for start, end in intervals if end > start)
    if not ordered:
        return []
    output = [ordered[0]]
    for start, end in ordered[1:]:
        old_start, old_end = output[-1]
        if start <= old_end:
            output[-1] = (old_start, max(old_end, end))
        else:
            output.append((start, end))
    return output


def _lift_has_infinite_ownership(lift: object) -> bool:
    return any(
        _extension_parts(token)[1] == "inf"
        for token in (lift.left_extension, lift.right_extension)
    )


def select_effective_annotations_for_gap_analysis(
    annotations: Sequence[AnnotationRow],
    lift_by_allele: Mapping[str, object],
    query_genome: str,
    max_impute: int = DEFAULT_MAX_IMPUTE,
) -> Tuple[List[AnnotationRow], Dict[str, int]]:
    """Apply the VCF ownership rule before discovering additive rows.

    Existing input rows remain immutable in the output file.  This function
    builds only the analysis view used for query-gap alignment and synthetic
    deletion discovery.  Every finite GenomeLift interval is a sentinel-score
    mask, including finite blocks for which no graph-CIGAR row survived.  A
    pure ``-inf`` row contributes analysis sequence only where it wins the
    shared score sweep and its complete alignment score is at least 1000.

    Rows containing a finite ownership name remain available as ordinary
    alignment context.  Rows with no matching GenomeLift name are also kept
    for backward compatibility (for example, rows added by an older run).
    """
    valid_regions_by_sample_contig: Dict[
        Tuple[str, str], List[Tuple[int, int]]
    ] = {}
    for lift in lift_by_allele.values():
        if _lift_has_infinite_ownership(lift):
            continue
        start, end = apply_interval_extensions_to_coord(
            lift.locus.start,
            lift.locus.end,
            lift.locus.strand,
            lift.left_extension,
            lift.right_extension,
            max_impute,
        )
        if end > start:
            valid_regions_by_sample_contig.setdefault(
                (query_genome, lift.locus.chrom), []
            ).append((start, end))
    valid_regions_by_sample_contig = {
        key: _merge_intervals(intervals)
        for key, intervals in valid_regions_by_sample_contig.items()
    }

    retained: List[AnnotationRow] = []
    infinite_rows: Dict[str, AnnotationRow] = {}
    candidates_by_sample_contig: Dict[
        Tuple[str, str], List[Tuple[int, int, int, str, int]]
    ] = {}
    stats = {
        "finite_or_legacy_rows": 0,
        "infinite_rows": 0,
        "infinite_rows_below_score": 0,
        "infinite_rows_selected": 0,
        "infinite_rows_excluded": 0,
        "infinite_selected_bases": 0,
        "infinite_competing_bases": 0,
    }

    for row_index, annotation in enumerate(annotations):
        matching_lifts = [
            lift_by_allele[allele]
            for allele in query_name_parts(annotation.query_name)
            if (
                allele in lift_by_allele
                and lift_by_allele[allele].locus.chrom
                == annotation.query_coord.chrom
            )
        ]
        infinite_lifts = [
            lift for lift in matching_lifts
            if _lift_has_infinite_ownership(lift)
        ]
        if not matching_lifts or len(infinite_lifts) != len(matching_lifts):
            retained.append(annotation)
            stats["finite_or_legacy_rows"] += 1
            continue

        stats["infinite_rows"] += 1
        score = infinite_fallback_alignment_score(annotation.graph_cigar)
        if score < MIN_INFINITE_FALLBACK_ALIGNMENT_SCORE:
            stats["infinite_rows_below_score"] += 1
            stats["infinite_rows_excluded"] += 1
            continue

        row_key = f"{annotation.line_number}\x1f{row_index}"
        candidate_intervals = _merge_intervals(
            (
                max(annotation.query_coord.start, lift.locus.start),
                min(annotation.query_coord.end, lift.locus.end),
            )
            for lift in infinite_lifts
        )
        if not candidate_intervals:
            stats["infinite_rows_excluded"] += 1
            continue
        infinite_rows[row_key] = annotation
        bucket = candidates_by_sample_contig.setdefault(
            (query_genome, annotation.query_coord.chrom), []
        )
        for start, end in candidate_intervals:
            bucket.append((
                start,
                end,
                score,
                row_key,
                annotation.line_number,
            ))

    selected_by_row, selected_bases, competing_bases = (
        resolve_infinite_fallback_candidates(
            candidates_by_sample_contig,
            valid_regions_by_sample_contig,
        )
    )
    for row_key, annotation in infinite_rows.items():
        intervals = _merge_intervals(selected_by_row.get(row_key, ()))
        if not intervals:
            stats["infinite_rows_excluded"] += 1
            continue
        stats["infinite_rows_selected"] += 1
        for start, end in intervals:
            retained.append(slice_annotation_query_interval(
                annotation, start, end,
            ))

    retained.sort(key=lambda row: (
        row.query_coord.chrom,
        row.query_coord.start,
        row.query_coord.end,
        row.line_number,
    ))
    stats["infinite_selected_bases"] = selected_bases
    stats["infinite_competing_bases"] = competing_bases
    return retained, stats


def read_ownership_lift_rows(
    path: str, query_genome: str,
) -> List[OwnershipLiftRow]:
    """Read the exact GenomeLift rows that a sparse override may replace."""
    rows: List[OwnershipLiftRow] = []
    seen = set()
    with open(path, "rt") as handle:
        for line_number, raw in enumerate(handle, 1):
            if not raw.strip() or raw.startswith("#"):
                continue
            fields = raw.rstrip("\r\n").split("\t")
            if not fields or fields[0] == "allelename":
                continue
            if len(fields) < 14 or fields[0] in {"DEL", "NA"}:
                continue
            try:
                genome = parse_allele_name(fields[0])[4]
            except ValueError:
                continue
            if genome != query_genome:
                continue
            match = re.fullmatch(r"part_(\d+)", fields[6].strip())
            part_index = int(match.group(1)) if match else 0

            def part_value(value: str) -> str:
                values = value.split(";")
                if part_index:
                    if len(values) < part_index:
                        raise ValueError(
                            f"{path}:{line_number}: part_{part_index} lacks "
                            f"its semicolon-delimited value in {value!r}"
                        )
                    return values[part_index - 1].strip()
                return value.strip()

            locus = core.parse_coord(part_value(fields[2]))
            if locus is None:
                raise ValueError(
                    f"{path}:{line_number}: malformed assembly coordinate "
                    f"{fields[2]!r}"
                )
            allele = fields[0]
            if allele in seen:
                raise ValueError(
                    f"{path}:{line_number}: duplicate query PA {allele!r}"
                )
            seen.add(allele)
            rows.append(OwnershipLiftRow(
                line_number=line_number,
                fields=tuple(fields),
                allele=allele,
                locus=locus,
                part_index=part_index,
                left_extension=part_value(fields[12]),
                right_extension=part_value(fields[13]),
            ))
    if not rows:
        raise ValueError(
            f"no GenomeLift PA rows found for {query_genome!r}: {path}"
        )
    return rows


def annotation_mapped_query_intervals(
    annotation: AnnotationRow,
) -> List[Tuple[int, int]]:
    """Return assembly intervals paired to a graph base by =/M/X."""
    local_position = 0
    physical_intervals: List[Tuple[int, int]] = []
    for segment in core.parse_graphic_segments(
        annotation.graph_cigar, annotation.query_name,
    ):
        for operation in segment.ops:
            query_size = core.query_consume(operation.op, operation.n)
            if operation.op in {"=", "M", "X"} and query_size:
                local_start = local_position
                local_end = local_position + query_size
                if annotation.query_coord.strand == "+":
                    start = annotation.query_coord.start + local_start
                    end = annotation.query_coord.start + local_end
                else:
                    start = annotation.query_coord.end - local_end
                    end = annotation.query_coord.end - local_start
                physical_intervals.append((start, end))
            local_position += query_size
    expected = annotation.query_coord.end - annotation.query_coord.start
    if local_position != expected:
        raise ValueError(
            f"{annotation.query_name}: graph CIGAR consumes {local_position} "
            f"query bases but {coord_text(annotation.query_coord)} spans "
            f"{expected}"
        )
    return _merge_intervals(physical_intervals)


def _intersect_interval_lists(
    left: Sequence[Tuple[int, int]],
    right: Sequence[Tuple[int, int]],
) -> List[Tuple[int, int]]:
    output: List[Tuple[int, int]] = []
    left_index = 0
    right_index = 0
    while left_index < len(left) and right_index < len(right):
        start = max(left[left_index][0], right[right_index][0])
        end = min(left[left_index][1], right[right_index][1])
        if end > start:
            output.append((start, end))
        if left[left_index][1] <= right[right_index][1]:
            left_index += 1
        else:
            right_index += 1
    return output


def _intersect_sorted_intervals_with_one(
    intervals: Sequence[Tuple[int, int]],
    start: int,
    end: int,
) -> List[Tuple[int, int]]:
    """Intersect merged/sorted intervals with one interval in O(log N + K)."""
    if end <= start or not intervals:
        return []
    interval_index = bisect.bisect_left(intervals, (start,))
    if interval_index and intervals[interval_index - 1][1] > start:
        interval_index -= 1
    output: List[Tuple[int, int]] = []
    while interval_index < len(intervals):
        interval_start, interval_end = intervals[interval_index]
        if interval_start >= end:
            break
        overlap_start = max(start, interval_start)
        overlap_end = min(end, interval_end)
        if overlap_end > overlap_start:
            output.append((overlap_start, overlap_end))
        interval_index += 1
    return output


def _uncovered_pieces(
    start: int, end: int, union: Sequence[Tuple[int, int]],
) -> List[Tuple[int, int]]:
    """Subtract merged/sorted coverage from one interval in O(log N + K)."""
    if end <= start:
        return []
    if not union:
        return [(start, end)]
    index = bisect.bisect_left(union, (start,))
    if index and union[index - 1][1] > start:
        index -= 1
    output: List[Tuple[int, int]] = []
    cursor = start
    while index < len(union) and cursor < end:
        covered_start, covered_end = union[index]
        if covered_start >= end:
            break
        if covered_start > cursor:
            output.append((cursor, min(end, covered_start)))
        cursor = max(cursor, covered_end)
        index += 1
    if cursor < end:
        output.append((cursor, end))
    return output


def _score_annotation_cigar(item: Tuple[int, str]) -> Tuple[int, int]:
    key, graph_cigar = item
    return key, infinite_fallback_alignment_score(graph_cigar)


def _score_annotations_parallel(
    rows: Sequence[AnnotationRow], processes: int,
) -> Dict[int, int]:
    """Score surviving resurrection rows, in parallel when worthwhile."""
    items = [(id(row), row.graph_cigar) for row in rows]
    worker_count = min(max(1, processes), len(items))
    if worker_count > 1 and len(items) >= 32:
        chunksize = max(1, len(items) // (worker_count * 4))
        with mp.get_context("fork").Pool(processes=worker_count) as pool:
            return dict(pool.imap_unordered(
                _score_annotation_cigar, items, chunksize,
            ))
    return dict(map(_score_annotation_cigar, items))


def _winning_ownership_segments(
    candidates: Sequence[OwnershipCandidate],
) -> List[Tuple[int, int, str, int, bool]]:
    """Sweep one contig while retaining only its single active PA owner."""
    if not candidates:
        return []
    events: List[Tuple[int, int, int]] = []
    for candidate_id, candidate in enumerate(candidates):
        if candidate.end <= candidate.start:
            continue
        events.append((candidate.start, 1, candidate_id))
        events.append((candidate.end, 0, candidate_id))
    events.sort(key=lambda value: (value[0], value[1], value[2]))

    active = set()
    priority_heap: List[Tuple[int, int, int, str, int, int]] = []
    output: List[Tuple[int, int, str, int, bool]] = []
    previous: Optional[int] = None
    event_index = 0
    while event_index < len(events):
        coordinate = events[event_index][0]
        if previous is not None and coordinate > previous:
            while priority_heap and priority_heap[0][-1] not in active:
                heapq.heappop(priority_heap)
            if priority_heap:
                winner = candidates[priority_heap[0][-1]]
                if (
                    output
                    and output[-1][1] == previous
                    and output[-1][2] == winner.allele
                    and output[-1][3] == winner.tier
                    and output[-1][4] == winner.non_primary
                ):
                    output[-1] = (
                        output[-1][0], coordinate, winner.allele, winner.tier,
                        winner.non_primary,
                    )
                else:
                    output.append((
                        previous, coordinate, winner.allele, winner.tier,
                        winner.non_primary,
                    ))

        while event_index < len(events) and events[event_index][0] == coordinate:
            _coordinate, event_type, candidate_id = events[event_index]
            candidate = candidates[candidate_id]
            if event_type == 1:
                active.add(candidate_id)
                # Highest tier, score, and span win. PA and source line are
                # stable ascending tie-breakers.
                heapq.heappush(priority_heap, (
                    -candidate.tier,
                    -candidate.score,
                    -(candidate.end - candidate.start),
                    candidate.allele,
                    candidate.line_number,
                    candidate_id,
                ))
            else:
                active.discard(candidate_id)
            event_index += 1
        previous = coordinate
    return output


def _interval_lists_overlap(
    left: Sequence[Tuple[int, int]],
    right: Sequence[Tuple[int, int]],
) -> bool:
    """Return whether two sorted/unsorted half-open interval sets overlap."""
    return bool(_intersect_interval_lists(
        _merge_intervals(left), _merge_intervals(right),
    ))


def _suppress_contained_resurrected_pas(
    candidates_by_contig: Mapping[str, Sequence[OwnershipCandidate]],
) -> set:
    """Discard a resurrected PA contained by another resurrected PA.

    Containment is deliberately decided before score competition.  It avoids
    allowing a short inner resurrected PA to split one outer PA into two
    disconnected ownership islands.  Equal spans are not containment; those
    continue to use alignment score and the stable sweep tie-breakers.
    """
    suppressed = set()
    for candidates in candidates_by_contig.values():
        intervals_by_allele: Dict[str, List[Tuple[int, int]]] = {}
        for candidate in candidates:
            if candidate.tier != OWNERSHIP_TIER_RESURRECTED:
                continue
            intervals_by_allele.setdefault(candidate.allele, []).append((
                candidate.start, candidate.end,
            ))
        envelopes = []
        for allele, intervals in intervals_by_allele.items():
            merged = _merge_intervals(intervals)
            if not merged:
                continue
            start = min(value[0] for value in merged)
            end = max(value[1] for value in merged)
            envelopes.append((start, end, allele))
        envelopes.sort(key=lambda value: (
            value[0], -value[1], value[2],
        ))
        active_outer: Optional[Tuple[int, int, str]] = None
        for start, end, allele in envelopes:
            if active_outer is not None:
                outer_start, outer_end, _outer_allele = active_outer
                if (
                    outer_end > start
                    and outer_end >= end
                    and (outer_start < start or outer_end > end)
                ):
                    suppressed.add(allele)
                    continue
                # Starts are nondecreasing, so the previous interval with the
                # farthest end is the only one that can contain this or any
                # later interval. Replacing it keeps the scan linear after
                # sorting even for a dense overlap pileup.
                if end <= outer_end:
                    continue
            active_outer = (start, end, allele)
    return suppressed


def determine_final_pa_ownership(
    annotations: Sequence[AnnotationRow],
    lift_rows: Sequence[OwnershipLiftRow],
    processes: int = 1,
    primary_reference_index: Optional[
        Mapping[str, Tuple[int, int, int, int]]
    ] = None,
) -> Tuple[
    Dict[str, Tuple[int, int]],
    Dict[str, Tuple[str, str]],
    Dict[str, int],
]:
    """Re-elect one connected ownership core per PA from graph evidence.

    Mapped finite PAs have the highest tier. Former ``-inf`` PAs whose full
    alignment score is >=1000 are resurrected at that score. Within either
    tier, alternative/novel-target alignments receive half priority. An outer
    resurrected PA consumes a fully contained resurrected PA before scoring
    unless that would let a non-primary outer suppress a primary inner PA.
    Unmapped-only PAs that overlap mapped evidence are discarded. A partially
    mapped PA may retain only the low-priority sequence connected to its
    winning mapped core.
    """
    lift_by_allele = {row.allele: row for row in lift_rows}
    candidates_by_contig: Dict[str, List[OwnershipCandidate]] = {}
    mapped_rows_by_allele: Dict[str, int] = {}
    resurrectable_rows_by_allele: Dict[str, int] = {}

    # Lowest-priority ownership covers wholly unmapped PAs and the clipped or
    # insertion-only portion of partially mapped PAs.
    for row in lift_rows:
        candidates_by_contig.setdefault(row.locus.chrom, []).append(
            OwnershipCandidate(
                row.locus.start,
                row.locus.end,
                row.allele,
                OWNERSHIP_TIER_UNMAPPED,
                0,
                row.line_number,
            )
        )

    # Column 3 is the block's owned query interval and its span was already
    # validated against the graph CIGAR. Column 4 classifies the coordinate
    # target without parsing the CIGAR in the usual case.
    deferred_resurrections: List[
        Tuple[AnnotationRow, str, object, List[Tuple[int, int]], bool]
    ] = []
    valid_regions_by_contig: Dict[str, List[Tuple[int, int]]] = {}
    for annotation in annotations:
        non_primary = annotation_targets_non_primary(
            annotation, primary_reference_index,
        )
        mapped_intervals = [(
            annotation.query_coord.start, annotation.query_coord.end,
        )]
        for allele in query_name_parts(annotation.query_name):
            lift = lift_by_allele.get(allele)
            if lift is None or lift.locus.chrom != annotation.query_coord.chrom:
                continue
            intersections = _intersect_sorted_intervals_with_one(
                mapped_intervals,
                lift.locus.start,
                lift.locus.end,
            )
            if not intersections:
                continue
            if _lift_has_infinite_ownership(lift):
                # Scored lazily below, only when the candidate reaches far
                # enough outside the valid union to matter.
                deferred_resurrections.append(
                    (annotation, allele, lift, intersections, non_primary)
                )
                continue
            mapped_rows_by_allele[allele] = (
                mapped_rows_by_allele.get(allele, 0) + 1
            )
            bucket = candidates_by_contig.setdefault(lift.locus.chrom, [])
            region_bucket = valid_regions_by_contig.setdefault(
                lift.locus.chrom, [],
            )
            for start, end in intersections:
                bucket.append(OwnershipCandidate(
                    start,
                    end,
                    allele,
                    OWNERSHIP_TIER_VALID_MAPPED,
                    ownership_priority_score(
                        annotation,
                        VALID_OWNERSHIP_PRIORITY_SCORE,
                        primary_reference_index,
                    ),
                    annotation.line_number,
                    non_primary,
                ))
                region_bucket.append((start, end))
            # Preserve GenomeLift's previous valid-owner decision as the
            # tie-breaker when two mapped valid PAs overlap. Positive gap
            # extensions are clipped back to the original PA locus here;
            # only negative overlap trims change this preferred core.
            owned_start, owned_end = apply_interval_extensions_to_coord(
                lift.locus.start,
                lift.locus.end,
                lift.locus.strand,
                lift.left_extension,
                lift.right_extension,
                DEFAULT_MAX_IMPUTE,
            )
            preferred = _intersect_sorted_intervals_with_one(
                intersections,
                max(lift.locus.start, owned_start),
                min(lift.locus.end, owned_end),
            )
            for start, end in preferred:
                bucket.append(OwnershipCandidate(
                    start,
                    end,
                    allele,
                    OWNERSHIP_TIER_VALID_MAPPED,
                    ownership_priority_score(
                        annotation,
                        VALID_OWNERSHIP_PRIORITY_SCORE,
                        primary_reference_index,
                    ) + 1,
                    annotation.line_number,
                    non_primary,
                ))

    valid_union_by_contig = {
        contig: _merge_intervals(intervals)
        for contig, intervals in valid_regions_by_contig.items()
    }

    # Gate 1 (>1000 bp): a resurrected candidate survives only when at least
    # one CONTIGUOUS uncovered stretch outside the valid union exceeds
    # MIN_INFINITE_FALLBACK_ALIGNMENT_SCORE bases; where the union covers it,
    # it can never out-rank valid ownership.
    gated_by_contig: Dict[
        str, List[Tuple[int, int, str, AnnotationRow, bool]]
    ] = {}
    for (
        annotation, allele, lift, intersections, non_primary
    ) in deferred_resurrections:
        union = valid_union_by_contig.get(lift.locus.chrom, ())
        best = 0
        for start, end in intersections:
            for piece_start, piece_end in _uncovered_pieces(
                start, end, union,
            ):
                best = max(best, piece_end - piece_start)
        if best <= MIN_INFINITE_FALLBACK_ALIGNMENT_SCORE:
            continue
        bucket = gated_by_contig.setdefault(lift.locus.chrom, [])
        for start, end in intersections:
            bucket.append((start, end, allele, annotation, non_primary))

    # Gate 2: exclude a resurrected candidate lying fully inside another
    # resurrected candidate on the same contig, before any scoring.
    surviving_resurrections: List[
        Tuple[int, int, str, AnnotationRow, bool]
    ] = []
    removed_alleles: set = set()
    surviving_alleles: set = set()
    for contig in sorted(gated_by_contig):
        entries = sorted(
            gated_by_contig[contig],
            key=lambda entry: (entry[0], entry[1], entry[2]),
        )
        removed = [False] * len(entries)
        for outer_index, entry in enumerate(entries):
            if removed[outer_index]:
                continue
            current_end = entry[1]
            for inner_index in range(outer_index + 1, len(entries)):
                if removed[inner_index]:
                    continue
                if entries[inner_index][0] >= current_end:
                    break
                if (
                    entries[inner_index][1] < current_end
                    and not (entry[4] and not entries[inner_index][4])
                ):
                    removed[inner_index] = True
        for entry_index, entry in enumerate(entries):
            if removed[entry_index]:
                removed_alleles.add(entry[2])
            else:
                surviving_resurrections.append(entry)
                surviving_alleles.add(entry[2])
    contained_resurrected = removed_alleles - surviving_alleles

    # Only the surviving resurrected candidates pay for a CIGAR score parse.
    unique_rows: Dict[int, AnnotationRow] = {}
    for (
        _start, _end, _allele, annotation, _non_primary
    ) in surviving_resurrections:
        unique_rows.setdefault(id(annotation), annotation)
    score_by_annotation = _score_annotations_parallel(
        list(unique_rows.values()), processes,
    )
    for (
        start, end, allele, annotation, non_primary
    ) in surviving_resurrections:
        alignment_score = score_by_annotation[id(annotation)]
        if alignment_score < MIN_INFINITE_FALLBACK_ALIGNMENT_SCORE:
            continue
        priority_score = ownership_priority_score(
            annotation, alignment_score, primary_reference_index,
        )
        mapped_rows_by_allele[allele] = (
            mapped_rows_by_allele.get(allele, 0) + 1
        )
        resurrectable_rows_by_allele[allele] = max(
            priority_score,
            resurrectable_rows_by_allele.get(allele, priority_score),
        )
        lift = lift_by_allele[allele]
        candidates_by_contig.setdefault(lift.locus.chrom, []).append(
            OwnershipCandidate(
                start,
                end,
                allele,
                OWNERSHIP_TIER_RESURRECTED,
                priority_score,
                annotation.line_number,
                non_primary,
            )
        )
    effective_mapped_alleles = (
        set(mapped_rows_by_allele) - contained_resurrected
    )
    mapped_union_by_contig: Dict[str, List[Tuple[int, int]]] = {}
    for contig, candidates in candidates_by_contig.items():
        mapped_union_by_contig[contig] = _merge_intervals(
            (candidate.start, candidate.end)
            for candidate in candidates
            if (
                candidate.tier > OWNERSHIP_TIER_UNMAPPED
                and candidate.allele not in contained_resurrected
            )
        )

    # A wholly unmapped PA has no anchor to which a surviving fragment could
    # remain connected. If any mapped PA overlaps it, remove the entire
    # unmapped PA instead of preserving an arbitrary exposed edge.
    suppressed_unmapped = set()
    unmapped_rows_by_contig: Dict[str, List[OwnershipLiftRow]] = {}
    for row in lift_rows:
        if row.allele not in effective_mapped_alleles:
            unmapped_rows_by_contig.setdefault(row.locus.chrom, []).append(row)
    for contig, unmapped_rows in unmapped_rows_by_contig.items():
        mapped_intervals = mapped_union_by_contig.get(contig, ())
        mapped_index = 0
        for row in sorted(unmapped_rows, key=lambda value: (
            value.locus.start, value.locus.end, value.allele,
        )):
            while (
                mapped_index < len(mapped_intervals)
                and mapped_intervals[mapped_index][1] <= row.locus.start
            ):
                mapped_index += 1
            if (
                mapped_index < len(mapped_intervals)
                and mapped_intervals[mapped_index][0] < row.locus.end
            ):
                suppressed_unmapped.add(row.allele)

    winning_candidates_by_contig: Dict[
        str, List[OwnershipCandidate]
    ] = {}
    for contig, candidates in candidates_by_contig.items():
        winning_candidates_by_contig[contig] = [
            candidate for candidate in candidates
            if (
                candidate.allele not in contained_resurrected
                and not (
                    candidate.tier == OWNERSHIP_TIER_UNMAPPED
                    and candidate.allele in suppressed_unmapped
                )
            )
        ]

    components_by_allele: Dict[str, List[Tuple[int, int]]] = {}
    mapped_wins_by_allele: Dict[str, List[Tuple[int, int]]] = {}
    primary_winning_alleles = set()
    non_primary_winning_alleles = set()
    for contig in sorted(candidates_by_contig):
        for (
            start, end, allele, tier, non_primary
        ) in _winning_ownership_segments(
            winning_candidates_by_contig[contig],
        ):
            components_by_allele.setdefault(allele, []).append((start, end))
            if tier > OWNERSHIP_TIER_UNMAPPED:
                mapped_wins_by_allele.setdefault(allele, []).append((
                    start, end,
                ))
                if non_primary:
                    non_primary_winning_alleles.add(allele)
                else:
                    primary_winning_alleles.add(allele)

    final_core: Dict[str, Tuple[int, int]] = {}
    discarded_components = 0
    ambiguous_mapped_pas = 0
    for row in lift_rows:
        components = _merge_intervals(
            components_by_allele.get(row.allele, ()),
        )
        if not components:
            continue
        if row.allele in effective_mapped_alleles:
            mapped_wins = _merge_intervals(
                mapped_wins_by_allele.get(row.allele, ()),
            )
            anchored = [
                component for component in components
                if _interval_lists_overlap([component], mapped_wins)
            ]
            if len(anchored) != 1:
                # Do not choose a largest scattered enclave. A mapped PA that
                # remains split after containment/priority resolution is
                # ambiguous and is removed as a whole.
                discarded_components += len(components)
                ambiguous_mapped_pas += 1
                continue
            final_core[row.allele] = anchored[0]
            discarded_components += max(0, len(components) - 1)
            continue

        # An unmapped-only PA is already excluded if mapped evidence touched
        # it. It may survive only as one connected low-priority interval.
        if len(components) != 1:
            discarded_components += len(components)
            continue
        final_core[row.allele] = components[0]

    offsets: Dict[str, Tuple[str, str]] = {
        row.allele: ("-inf", "-inf") for row in lift_rows
    }
    rows_by_contig: Dict[str, List[OwnershipLiftRow]] = {}
    for row in lift_rows:
        if row.allele in final_core:
            rows_by_contig.setdefault(row.locus.chrom, []).append(row)
    for contig, contig_rows in rows_by_contig.items():
        regions = [
            (*final_core[row.allele], row.allele)
            for row in contig_rows
        ]
        unique_map = {
            row.allele: (
                final_core[row.allele][0],
                final_core[row.allele][1],
                final_core[row.allele][1] - final_core[row.allele][0],
                0,
                0,
            )
            for row in contig_rows
        }
        offsets.update(compute_priority_display_offsets_for_regions(
            regions, unique_map,
        ))

    # A PA backed only by an alternative/novel target owns its elected core
    # and no adjacent query gap. Classify the PA by the evidence that actually
    # won ownership, rather than by lower-priority alignments that lost.
    non_extendable_alleles = (
        non_primary_winning_alleles - primary_winning_alleles
    )
    for allele in non_extendable_alleles:
        if allele in final_core:
            offsets[allele] = ("0", "0")

    resurrected = sum(
        allele in final_core
        for allele in resurrectable_rows_by_allele
        if allele not in contained_resurrected
    )
    stats = {
        "pas": len(lift_rows),
        "mapped_pas": len(effective_mapped_alleles),
        "resurrected_pas": resurrected,
        "unmapped_only_pas": sum(
            row.allele in final_core and row.allele not in mapped_rows_by_allele
            for row in lift_rows
        ),
        "unowned_pas": len(lift_rows) - len(final_core),
        "discarded_components": discarded_components,
        "contained_resurrected_pas": len(contained_resurrected),
        "suppressed_unmapped_pas": len(suppressed_unmapped),
        "ambiguous_mapped_pas": ambiguous_mapped_pas,
        "non_extendable_pas": sum(
            allele in final_core for allele in non_extendable_alleles
        ),
    }
    return final_core, offsets, stats


def apply_final_ownership_to_lifts(
    lift_by_allele: Mapping[str, object],
    final_core: Mapping[str, Tuple[int, int]],
    offsets: Mapping[str, Tuple[str, str]],
) -> Dict[str, object]:
    """Return an in-memory GenomeLift overlay matching genomeliftfix.tsv."""
    output: Dict[str, object] = {}
    for allele, lift in lift_by_allele.items():
        core_interval = final_core.get(allele)
        if core_interval is None:
            locus = lift.locus
        else:
            locus = core.Coord(
                lift.locus.chrom,
                core_interval[0],
                core_interval[1],
                lift.locus.strand,
            )
        left_extension, right_extension = offsets.get(
            allele, ("-inf", "-inf"),
        )
        replacements = {
            "locus": locus,
            "left_extension": left_extension,
            "right_extension": right_extension,
        }
        if hasattr(lift, "locus_text"):
            replacements["locus_text"] = coord_text(locus)
        output[allele] = dataclasses.replace(lift, **replacements)
    return output


def select_annotations_for_final_ownership(
    annotations: Sequence[AnnotationRow],
    lift_by_allele: Mapping[str, object],
    max_impute: int = DEFAULT_MAX_IMPUTE,
) -> Tuple[List[AnnotationRow], Dict[str, int]]:
    """Filter the analysis view to rows with at least one surviving owner.

    Ownership was already finalized and swapped into the GenomeLift overlay.
    The graph-CIGAR rows were generated per block independently of ownership
    in graphcigartoref_persample.py, so the analysis view keeps every
    surviving row byte-identical: no interval clipping, no per-PA expansion,
    and no CIGAR work.  A row is excluded only when every matching PA lost
    ownership (pure ``-inf`` after finalization); rows with no matching
    GenomeLift name are kept for backward compatibility.
    """
    del max_impute
    retained: List[AnnotationRow] = []
    stats = {
        "input_rows": len(annotations),
        "legacy_rows": 0,
        "logical_pa_rows": 0,
        "excluded_pa_rows": 0,
        "excluded_rows": 0,
        "retained_bases": 0,
    }
    for annotation in annotations:
        matched = [
            lift_by_allele[allele]
            for allele in query_name_parts(annotation.query_name)
            if (
                allele in lift_by_allele
                and lift_by_allele[allele].locus.chrom
                == annotation.query_coord.chrom
            )
        ]
        if matched:
            stats["logical_pa_rows"] += len(matched)
            owned = [
                lift for lift in matched
                if not _lift_has_infinite_ownership(lift)
            ]
            stats["excluded_pa_rows"] += len(matched) - len(owned)
            if not owned:
                stats["excluded_rows"] += 1
                continue
        else:
            stats["legacy_rows"] += 1
        retained.append(annotation)
        stats["retained_bases"] += (
            annotation.query_coord.end - annotation.query_coord.start
        )

    retained.sort(key=lambda row: (
        row.query_coord.chrom,
        row.query_coord.start,
        row.query_coord.end,
        row.line_number,
        row.query_name,
    ))
    return retained, stats


def _replace_part_value(
    fields: List[str], field_index: int, part_index: int, value: str,
) -> None:
    if not part_index:
        fields[field_index] = value
        return
    values = fields[field_index].split(";")
    if len(values) < part_index:
        raise ValueError(
            f"part_{part_index} cannot update field {field_index + 1}: "
            f"{fields[field_index]!r}"
        )
    values[part_index - 1] = value
    fields[field_index] = ";".join(values)


def write_genomeliftfix(
    path: str,
    lift_rows: Sequence[OwnershipLiftRow],
    final_core: Mapping[str, Tuple[int, int]],
    offsets: Mapping[str, Tuple[str, str]],
) -> int:
    """Write only PA rows whose coordinate ownership actually changed."""
    path = os.path.abspath(os.path.expanduser(path))
    os.makedirs(os.path.dirname(path) or ".", exist_ok=True)
    temporary = path + f".tmp.{os.getpid()}"
    updated = 0
    try:
        with open(temporary, "wt") as output:
            output.write(
                "# sparse GenomeLift ownership overrides; later "
                "--coord-map entries replace earlier PA names\n"
            )
            for row in lift_rows:
                fields = list(row.fields)
                core_interval = final_core.get(row.allele)
                if core_interval is None:
                    new_coord = coord_text(row.locus)
                else:
                    new_coord = coord_text(core.Coord(
                        row.locus.chrom,
                        core_interval[0],
                        core_interval[1],
                        row.locus.strand,
                    ))
                new_left, new_right = offsets.get(
                    row.allele, ("-inf", "-inf"),
                )

                old_fields = tuple(fields)
                _replace_part_value(
                    fields, 2, row.part_index, new_coord,
                )
                _replace_part_value(
                    fields, 12, row.part_index, new_left,
                )
                _replace_part_value(
                    fields, 13, row.part_index, new_right,
                )
                if tuple(fields) == old_fields:
                    continue
                output.write("\t".join(fields) + "\n")
                updated += 1
        os.replace(temporary, path)
    finally:
        try:
            os.remove(temporary)
        except FileNotFoundError:
            pass
    return updated


def default_genomeliftfix_path(genomelift_path: str) -> str:
    directory = os.path.dirname(os.path.abspath(genomelift_path))
    basename = os.path.basename(genomelift_path)
    if basename.endswith(".genomelift.tsv"):
        basename = basename[:-len(".genomelift.tsv")] + ".genomeliftfix.tsv"
    elif basename.endswith(".tsv"):
        basename = basename[:-4] + ".fix.tsv"
    else:
        basename += ".fix.tsv"
    return os.path.join(directory, basename)


def _subtract_coverage(
    start: int,
    end: int,
    coverage: Sequence[Tuple[int, int]],
) -> List[Tuple[int, int]]:
    if end <= start:
        return []
    output: List[Tuple[int, int]] = []
    cursor = start
    # ``coverage`` is merged/sorted by reference_coverage_by_contig(). Jump
    # directly to the interval that can touch ``start`` instead of rescanning
    # the contig from its first interval for every candidate deletion gap.
    # This changes G repeated queries from O(G*C) to O(G log C + overlaps).
    coverage_index = bisect.bisect_left(coverage, (start,))
    if coverage_index and coverage[coverage_index - 1][1] > start:
        coverage_index -= 1
    while coverage_index < len(coverage):
        covered_start, covered_end = coverage[coverage_index]
        coverage_index += 1
        if covered_end <= cursor:
            continue
        if covered_start >= end:
            break
        if covered_start > cursor:
            output.append((cursor, min(end, covered_start)))
        cursor = max(cursor, covered_end)
        if cursor >= end:
            break
    if cursor < end:
        output.append((cursor, end))
    return output


def reference_coverage_by_contig(
    annotations: Sequence[AnnotationRow],
    reference_index: Mapping[str, Tuple[int, int, int, int]],
    represented_deletions: Iterable[Tuple[str, int, int]] = (),
) -> Dict[str, List[Tuple[int, int]]]:
    """Return merged reference bases represented anywhere in the query."""
    coverage: Dict[str, List[Tuple[int, int]]] = {}
    for row in annotations:
        for segment in core.parse_graphic_segments(
            row.graph_cigar, row.query_name,
        ):
            path = _resolve_reference_path(segment.path, reference_index)
            if path is None or segment.end <= segment.start:
                continue
            coverage.setdefault(path, []).append((
                min(segment.start, segment.end),
                max(segment.start, segment.end),
            ))
    # Existing explicit/synthetic deletion rows already represent a called
    # absence.  Count them here solely to avoid appending the same call twice.
    for path, start, end in represented_deletions:
        resolved = _resolve_reference_path(path, reference_index)
        if resolved is not None and end > start:
            coverage.setdefault(resolved, []).append((start, end))
    return {
        path: _merge_intervals(intervals)
        for path, intervals in coverage.items()
    }


def read_represented_deletions(path: str) -> List[Tuple[str, int, int]]:
    output: List[Tuple[str, int, int]] = []
    with open(path, "rt") as handle:
        for raw in handle:
            if not raw.strip() or raw.startswith("#"):
                continue
            fields = raw.rstrip("\r\n").split("\t")
            if len(fields) < 7 or not fields[0].startswith("DEL"):
                continue
            coord = core.parse_coord(fields[3])
            if coord is not None and coord.end > coord.start:
                output.append((coord.chrom, coord.start, coord.end))
    return output


def discover_reference_gap_deletions(
    annotations: Sequence[AnnotationRow],
    lift_by_allele: Mapping[str, object],
    reference_index: Mapping[str, Tuple[int, int, int, int]],
    query_genome: str,
    represented_deletions: Iterable[Tuple[str, int, int]] = (),
) -> Tuple[List[SyntheticDeletionRow], Dict[str, int]]:
    """Call uncovered reference gaps only at an observed query breakpoint.

    Missing graph coverage is not, by itself, evidence of a deletion.  The
    neighboring query annotations must meet exactly: if query bases remain
    between them, that interval is unresolved and must not be converted into
    a reference deletion merely because gap alignment was skipped or failed.
    """
    represented_deletions = tuple(represented_deletions)
    coverage = reference_coverage_by_contig(
        annotations, reference_index, represented_deletions,
    )
    by_contig: Dict[str, List[AnnotationRow]] = {}
    for row in annotations:
        by_contig.setdefault(row.query_coord.chrom, []).append(row)
    output: List[SyntheticDeletionRow] = []
    seen: set = set()
    stats = {
        "positive_reference_gaps": 0,
        "reference_gap_bases": 0,
        "reference_gap_covered_bases": 0,
        "reference_deletions": 0,
        "reference_deletion_bases": 0,
        "positive_query_gaps_skipped": 0,
        "positive_query_gap_bases_skipped": 0,
        "query_overlaps_skipped": 0,
    }
    for contig in sorted(by_contig):
        ordered = sorted(by_contig[contig], key=lambda row: (
            row.query_coord.start,
            row.query_coord.end,
            row.line_number,
        ))
        for left, right in zip(ordered, ordered[1:]):
            if left.alternative_intervals is not None or right.alternative_intervals is not None:
                continue
            relation = _signed_reference_gap(
                left, right, reference_index,
            )
            if relation is None or relation[0] <= 0:
                continue
            gap_size, left_boundary, right_boundary = relation
            path = _resolve_reference_path(
                left_boundary.path, reference_index,
            )
            if path is None:
                continue
            gap_start = min(left_boundary.position, right_boundary.position)
            gap_end = max(left_boundary.position, right_boundary.position)
            stats["positive_reference_gaps"] += 1
            stats["reference_gap_bases"] += gap_size
            query_gap = right.query_coord.start - left.query_coord.end
            if query_gap > 0:
                stats["positive_query_gaps_skipped"] += 1
                stats["positive_query_gap_bases_skipped"] += query_gap
                continue
            if query_gap < 0:
                # Overlapping query observations do not define a deletion
                # breakpoint either.  They require ownership reconciliation,
                # not an absence call inferred from missing graph coverage.
                stats["query_overlaps_skipped"] += 1
                continue
            uncovered = _subtract_coverage(
                gap_start, gap_end, coverage.get(path, ()),
            )
            uncovered_bases = sum(end - start for start, end in uncovered)
            stats["reference_gap_covered_bases"] += gap_size - uncovered_bases
            if not uncovered:
                continue
            selected_side = _pair_priority_side(
                left, right, lift_by_allele,
            ) or "left"
            selected = left if selected_side == "left" else right
            query_breakpoint = (
                left.query_coord.end + right.query_coord.start
            ) // 2
            query_strand = (
                left.query_coord.strand
                if left.query_coord.strand == right.query_coord.strand
                else "+"
            )
            for deletion_start, deletion_end in uncovered:
                key = (path, deletion_start, deletion_end, contig, query_breakpoint)
                if key in seen:
                    continue
                seen.add(key)
                safe_path = re.sub(r"[^A-Za-z0-9]+", "p", path).strip("p")
                query_name = (
                    f"DEL{safe_path}p{deletion_start}p{deletion_end}_"
                    f"{query_genome}_{len(output) + 1}"
                )
                output.append(SyntheticDeletionRow(
                    query_name=query_name,
                    query_coord=core.Coord(
                        contig,
                        query_breakpoint,
                        query_breakpoint,
                        query_strand,
                    ),
                    ref_name=selected.reference_name,
                    ref_coord=core.Coord(
                        path, deletion_start, deletion_end, "+",
                    ),
                ))
                stats["reference_deletions"] += 1
                stats["reference_deletion_bases"] += (
                    deletion_end - deletion_start
                )
    return output, stats


def _build_tandem_duplication_row(
    task: GapTask,
    query_reader: core.FastaRegionReader,
    reference_reader: core.FastaRegionReader,
) -> str:
    if task.insertion_start is None or task.insertion_end is None:
        raise ValueError(f"{task.name}: tandem duplication lacks its query interval")
    if task.insertion_end <= task.insertion_start:
        raise ValueError(f"{task.name}: tandem duplication query interval is empty")
    selected_boundary = (
        task.left_boundary
        if task.selected_side == "left"
        else task.right_boundary
    )
    path = _resolve_reference_path(
        selected_boundary.path, reference_reader.index,
    )
    if path is None:
        raise KeyError(
            f"reference path {selected_boundary.path!r} is absent from the FAI"
        )
    position = selected_boundary.position
    contig_size = reference_reader.index[path][0]
    if position < 0 or position > contig_size:
        raise ValueError(
            f"{task.name}: insertion breakpoint {path}:{position} is outside "
            f"the reference contig length {contig_size}"
        )
    insertion_sequence = query_reader.fetch(
        task.contig,
        task.insertion_start,
        task.insertion_end,
        "+",
    )
    direction = ">" if selected_boundary.strand == "+" else "<"
    graph_cigar = core._format_interval_gcigar(
        direction,
        path,
        contig_size,
        position,
        position,
        (core.CigarOp(len(insertion_sequence), "I", insertion_sequence),),
    )
    if core.graph_cigar_query_span(graph_cigar, task.name) != len(
        insertion_sequence
    ):
        raise ValueError(f"{task.name}: pure-insertion graph CIGAR span mismatch")
    overlap_start = min(
        task.left_boundary.position, task.right_boundary.position,
    )
    overlap_end = max(
        task.left_boundary.position, task.right_boundary.position,
    )
    query_coord = core.Coord(
        task.contig, task.insertion_start, task.insertion_end, "+",
    )
    reference_coord = core.Coord(
        path, position, position, selected_boundary.strand,
    )
    overlap_coord = core.Coord(
        path, overlap_start, overlap_end, selected_boundary.strand,
    )
    return "\t".join((
        task.name,
        task.selected.reference_name,
        coord_text(query_coord),
        coord_text(reference_coord),
        coord_text(query_coord),
        coord_text(overlap_coord),
        graph_cigar,
    ))


def build_gap_row(
    task: GapTask,
    query_reader: core.FastaRegionReader,
    reference_reader: core.FastaRegionReader,
    anchor: int = DEFAULT_ANCHOR,
    maximum_reference_extension: int = DEFAULT_REFERENCE_EXTENSION,
) -> str:
    if task.mode == "tandem-duplication":
        return _build_tandem_duplication_row(
            task, query_reader, reference_reader,
        )
    query_start, query_end, left_anchor, _right_anchor = _query_window(
        task, query_reader, anchor,
    )
    query_sequence = query_reader.fetch(
        task.contig, query_start, query_end, "+",
    )
    path, reference_start, reference_end, reference_strand = _reference_window(
        task, reference_reader, anchor, maximum_reference_extension,
    )
    reference_sequence = reference_reader.fetch(
        path, reference_start, reference_end, reference_strand,
    )
    operations = _align_reference_flank_ops(
        reference_sequence,
        query_sequence,
        # The existing upstream rescue aligns core-first, then restores the
        # original far-to-core traversal.  Reuse that tie-breaking behavior
        # for a gap owned by its right-hand neighbor.
        reverse_from_core=(
            task.mode == "one-sided" and task.selected_side == "right"
        ),
    )
    direction = ">" if reference_strand == "+" else "<"
    complete_cigar = core._format_interval_gcigar(
        direction,
        path,
        reference_reader.index[path][0],
        reference_start,
        reference_end,
        operations,
    )
    gap_size = task.end - task.start
    graph_cigar = core.slice_graph_cigar_by_query(
        complete_cigar,
        left_anchor,
        left_anchor + gap_size,
        task.name,
        # A one-sided alignment must remain connected to the selected
        # neighboring block's reference boundary.  Keep a terminal deletion
        # only on that core-facing side; terminal reference sequence on the
        # far/search-window side is not owned by this query gap.  Bounded
        # alignments use query anchors on both sides and retain neither outer
        # boundary deletion after those anchors are removed.
        include_left_boundary_deletions=(
            task.mode == "one-sided" and task.selected_side == "left"
        ),
        include_right_boundary_deletions=(
            task.mode == "one-sided" and task.selected_side == "right"
        ),
    )
    if core.graph_cigar_query_span(graph_cigar, task.name) != gap_size:
        raise ValueError(f"{task.name}: anchor truncation changed the query gap span")
    reference_coord = _reference_coord_from_cigar(graph_cigar, task.name)
    graph_cigar = query_only_graph_cigar(graph_cigar)
    return "\t".join((
        task.name,
        task.selected.reference_name,
        coord_text(core.Coord(task.contig, task.start, task.end, "+")),
        coord_text(reference_coord),
        coord_text(core.Coord(task.contig, query_start, query_end, "+")),
        coord_text(core.Coord(
            path,
            reference_start,
            reference_end,
            reference_strand,
        )),
        graph_cigar,
    ))


def _init_worker(
    query_path: str,
    reference_path: str,
    local_reference_templates: str,
    anchor: int,
    maximum_reference_extension: int,
) -> None:
    global _WORKER_QUERY_READER
    global _WORKER_REFERENCE_READER
    global _WORKER_ANCHOR
    global _WORKER_REFERENCE_EXTENSION
    _WORKER_QUERY_READER = core.FastaRegionReader(query_path)
    _WORKER_REFERENCE_READER = core.FastaRegionReader(
        reference_path,
        template_fasta=(local_reference_templates or None),
    )
    _WORKER_ANCHOR = anchor
    _WORKER_REFERENCE_EXTENSION = maximum_reference_extension


def _worker(task: GapTask) -> GapResult:
    task_started = time.monotonic()
    try:
        if _WORKER_QUERY_READER is None or _WORKER_REFERENCE_READER is None:
            raise RuntimeError("gap worker FASTA readers were not initialized")
        result = GapResult(
            task.index,
            build_gap_row(
                task,
                _WORKER_QUERY_READER,
                _WORKER_REFERENCE_READER,
                _WORKER_ANCHOR,
                _WORKER_REFERENCE_EXTENSION,
            ),
        )
    except Exception:
        result = GapResult(task.index, None, traceback.format_exc())
    task_elapsed = time.monotonic() - task_started
    if task_elapsed > 10:
        sys.stderr.write(
            "[fill_graphcigartoref_gaps] slow gap task "
            f"{task.name} ({task.end - task.start}bp/{task.mode}) "
            f"took {task_elapsed:.0f}s\n"
        )
    return result


def align_gap_tasks(
    tasks: Sequence[GapTask],
    query_path: str,
    reference_path: str,
    local_reference_templates: str,
    processes: int,
    chunksize: int,
    start_method: str,
    anchor: int,
    maximum_reference_extension: int,
) -> Tuple[List[str], List[GapResult]]:
    if not tasks:
        return [], []
    largest = sorted(
        tasks, key=lambda task: task.end - task.start, reverse=True,
    )[:5]
    sys.stderr.write(
        "[fill_graphcigartoref_gaps] aligning "
        f"{len(tasks)} gap task(s); largest: "
        + ", ".join(
            f"{task.name}={task.end - task.start}bp/{task.mode}"
            for task in largest
        )
        + "\n"
    )
    started = time.monotonic()
    completed = 0
    results: List[GapResult] = []

    # Small gaps fan out across worker processes with single-threaded
    # aligners; large gaps run one at a time in this process with minimap2
    # given every requested thread.
    small_tasks = [task for task in tasks if task.end - task.start < 1_000]
    large_tasks = [task for task in tasks if task.end - task.start >= 1_000]
    aligner_threads = max(1, processes)
    sys.stderr.write(
        "[fill_graphcigartoref_gaps] scheduling "
        f"{len(small_tasks)} small (<1000bp) task(s) across "
        f"{min(processes, max(1, len(small_tasks)))} process(es) and "
        f"{len(large_tasks)} large task(s) one-by-one with "
        f"{aligner_threads} aligner thread(s)\n"
    )

    def report_progress(force: bool = False) -> None:
        if force or completed % 25 == 0:
            elapsed = max(1e-6, time.monotonic() - started)
            sys.stderr.write(
                "[fill_graphcigartoref_gaps] aligned "
                f"{completed}/{len(tasks)} gap task(s) "
                f"({elapsed:.0f}s elapsed)\n"
            )

    def run_serial(
        serial_tasks: Sequence[GapTask],
        query_reader: core.FastaRegionReader,
        reference_reader: core.FastaRegionReader,
    ) -> None:
        nonlocal completed
        for task in serial_tasks:
            task_started = time.monotonic()
            try:
                row = build_gap_row(
                    task,
                    query_reader,
                    reference_reader,
                    anchor,
                    maximum_reference_extension,
                )
                results.append(GapResult(task.index, row))
            except Exception:
                results.append(
                    GapResult(task.index, None, traceback.format_exc())
                )
            completed += 1
            task_elapsed = time.monotonic() - task_started
            if task_elapsed > 10:
                sys.stderr.write(
                    "[fill_graphcigartoref_gaps] slow gap task "
                    f"{task.name} ({task.end - task.start}bp/{task.mode}) "
                    f"took {task_elapsed:.0f}s\n"
                )
            report_progress()

    worker_count = min(processes, len(small_tasks))
    if small_tasks and worker_count > 1:
        context = mp.get_context(start_method)
        actual_chunksize = chunksize or max(
            1, min(64, len(small_tasks) // max(1, worker_count * 8)),
        )
        with context.Pool(
            processes=worker_count,
            initializer=_init_worker,
            initargs=(
                query_path,
                reference_path,
                local_reference_templates,
                anchor,
                maximum_reference_extension,
            ),
        ) as pool:
            result_iter = pool.imap_unordered(
                _worker, small_tasks, actual_chunksize,
            )
            pool_pending = len(small_tasks)
            while pool_pending > 0:
                try:
                    result = result_iter.next(timeout=60)
                except mp.TimeoutError:
                    elapsed = max(1e-6, time.monotonic() - started)
                    sys.stderr.write(
                        "[fill_graphcigartoref_gaps] still aligning: "
                        f"{completed}/{len(tasks)} gap task(s) done after "
                        f"{elapsed:.0f}s\n"
                    )
                    continue
                except StopIteration:
                    break
                results.append(result)
                pool_pending -= 1
                completed += 1
                report_progress()

    serial_small = small_tasks if small_tasks and worker_count <= 1 else []
    if serial_small or large_tasks:
        query_reader = core.FastaRegionReader(query_path)
        reference_reader = core.FastaRegionReader(
            reference_path,
            template_fasta=(local_reference_templates or None),
        )
        try:
            if serial_small:
                run_serial(serial_small, query_reader, reference_reader)
            if large_tasks:
                previous_threads = os.environ.get(
                    "GRAPH_CIGARTOREF_MINIMAP2_THREADS"
                )
                os.environ["GRAPH_CIGARTOREF_MINIMAP2_THREADS"] = str(
                    aligner_threads
                )
                try:
                    run_serial(large_tasks, query_reader, reference_reader)
                finally:
                    if previous_threads is None:
                        os.environ.pop(
                            "GRAPH_CIGARTOREF_MINIMAP2_THREADS", None,
                        )
                    else:
                        os.environ["GRAPH_CIGARTOREF_MINIMAP2_THREADS"] = (
                            previous_threads
                        )
        finally:
            query_reader.handle.close()
            reference_reader.close()
    report_progress(force=True)
    results.sort(key=lambda result: result.index)
    return [result.row for result in results if result.row is not None], results


def _all_input_rows(input_path: str) -> List[str]:
    """Return every nonempty input row unchanged and in original order."""
    output: List[str] = []
    with open(input_path, "rt") as source:
        for raw in source:
            text = raw.rstrip("\r\n")
            if not text:
                continue
            output.append(text)
    return output


def _write_output(
    input_path: str,
    output_path: str,
    base_rows: Sequence[str],
    added_rows: Sequence[str],
    gaps_only: bool,
) -> None:
    if not output_path:
        if not gaps_only:
            for row in base_rows:
                sys.stdout.write(row + "\n")
        for row in added_rows:
            sys.stdout.write(row + "\n")
        return

    output_path = os.path.abspath(os.path.expanduser(output_path))
    input_path = os.path.abspath(os.path.expanduser(input_path))
    same_path = output_path == input_path
    if same_path:
        directory = os.path.dirname(output_path) or "."
        fd, temporary = tempfile.mkstemp(
            prefix=os.path.basename(output_path) + ".gapfill.",
            dir=directory,
            text=True,
        )
        os.close(fd)
    else:
        temporary = output_path
    try:
        with open(temporary, "wt") as output:
            if not gaps_only:
                for row in base_rows:
                    output.write(row + "\n")
            for row in added_rows:
                output.write(row + "\n")
        if same_path:
            os.replace(temporary, output_path)
    finally:
        if same_path:
            try:
                os.remove(temporary)
            except FileNotFoundError:
                pass


def parse_args(argv: Optional[Sequence[str]] = None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "-i", "--input", required=True,
        help="seven-column graphcigartoref_persample.py output",
    )
    parser.add_argument("--genomelift", required=True)
    parser.add_argument(
        "--genomelift-fix-output",
        default="",
        help=(
            "sparse corrected GenomeLift ownership TSV; defaults beside "
            "--genomelift with .genomeliftfix.tsv/.fix.tsv suffix"
        ),
    )
    parser.add_argument(
        "-q", "--fasta-query", required=True,
        help="faidxed query assembly FASTA",
    )
    parser.add_argument(
        "-r", "--reference", required=True,
        help="faidxed reference FASTA used by the input graph CIGAR",
    )
    parser.add_argument(
        "--local-reference-templates",
        default="",
        help="optional sparse local-template FASTA for reference-free graphs",
    )
    parser.add_argument("-o", "--output", default="")
    parser.add_argument(
        "--gaps-only", action="store_true",
        help="write only new gap rows instead of copying the input first",
    )
    parser.add_argument(
        "--query-genome", default="",
        help="sample_haplotype for gap names; inferred when omitted",
    )
    parser.add_argument("-t", "--processes", type=int, default=1)
    parser.add_argument("--chunksize", type=int, default=0)
    parser.add_argument(
        "--mp-start-method",
        choices=("fork", "spawn", "forkserver"),
        default="fork",
    )
    parser.add_argument("--anchor", type=int, default=DEFAULT_ANCHOR)
    parser.add_argument("--maximum-gap", type=int, default=DEFAULT_MAXIMUM_GAP)
    parser.add_argument(
        "--minimum-unmasked", type=int, default=DEFAULT_MINIMUM_UNMASKED,
    )
    parser.add_argument(
        "--minimum-unmasked-fraction",
        type=float,
        default=DEFAULT_MINIMUM_UNMASKED_FRACTION,
    )
    parser.add_argument(
        "--maximum-reference-extension",
        type=int,
        default=DEFAULT_REFERENCE_EXTENSION,
    )
    parser.add_argument(
        "--max-impute",
        type=int,
        default=DEFAULT_MAX_IMPUTE,
        help=(
            "ownership extension cap; keep equal to graphreftovcf.py "
            "--max-impute (default: 10000)"
        ),
    )
    parser.add_argument(
        "--strict", action="store_true",
        help="fail when any eligible gap cannot be aligned",
    )
    args = parser.parse_args(argv)
    if args.processes < 1:
        parser.error("--processes must be >= 1")
    if args.chunksize < 0:
        parser.error("--chunksize must be >= 0")
    if args.anchor < 1:
        parser.error("--anchor must be >= 1")
    if args.maximum_gap < 1:
        parser.error("--maximum-gap must be >= 1")
    if args.minimum_unmasked < 0:
        parser.error("--minimum-unmasked must be >= 0")
    if not 0 <= args.minimum_unmasked_fraction <= 1:
        parser.error("--minimum-unmasked-fraction must be between 0 and 1")
    if args.maximum_reference_extension < 1:
        parser.error("--maximum-reference-extension must be >= 1")
    if args.max_impute < 0:
        parser.error("--max-impute must be >= 0")
    return args


def run(args: argparse.Namespace) -> int:
    run_started = time.monotonic()
    input_path = os.path.abspath(os.path.expanduser(args.input))
    genomelift_path = os.path.abspath(os.path.expanduser(args.genomelift))
    genomelift_fix_path = (
        os.path.abspath(os.path.expanduser(args.genomelift_fix_output))
        if args.genomelift_fix_output
        else default_genomeliftfix_path(genomelift_path)
    )
    query_path = os.path.abspath(os.path.expanduser(args.fasta_query))
    reference_path = os.path.abspath(os.path.expanduser(args.reference))
    local_reference_templates = (
        os.path.abspath(os.path.expanduser(args.local_reference_templates))
        if args.local_reference_templates else ""
    )
    for path in (
        input_path, genomelift_path, query_path, reference_path,
        local_reference_templates,
    ):
        if not path:
            continue
        if not os.path.isfile(path):
            raise FileNotFoundError(path)

    # This index deliberately excludes --local-reference-templates. It is the
    # authoritative distinction between primary-reference targets and the
    # alternative/novel loci added by that sparse template catalog.
    primary_reader = core.FastaRegionReader(reference_path)
    try:
        primary_reference_index = dict(primary_reader.index)
    finally:
        primary_reader.close()

    annotations = read_annotations(input_path)
    lift_rows = read_genomelift(genomelift_path)
    lift_by_allele = {row.allele: row for row in lift_rows}
    query_genome = args.query_genome or infer_query_genome(
        annotations, lift_by_allele,
    )
    selected_lift_rows = {}
    for allele, row in lift_by_allele.items():
        try:
            genome = parse_allele_name(allele)[4]
        except ValueError:
            continue
        if genome == query_genome:
            selected_lift_rows[allele] = row
    lift_by_allele = selected_lift_rows

    ownership_started = time.monotonic()
    ownership_lift_rows = read_ownership_lift_rows(
        genomelift_path, query_genome,
    )
    final_core, final_offsets, final_ownership_stats = (
        determine_final_pa_ownership(
            annotations,
            ownership_lift_rows,
            args.processes,
            primary_reference_index,
        )
    )
    updated_pa_count = write_genomeliftfix(
        genomelift_fix_path,
        ownership_lift_rows,
        final_core,
        final_offsets,
    )
    lift_by_allele = apply_final_ownership_to_lifts(
        lift_by_allele,
        final_core,
        final_offsets,
    )
    analysis_annotations, ownership_stats = (
        select_annotations_for_final_ownership(
            annotations,
            lift_by_allele,
            args.max_impute,
        )
    )
    ownership_elapsed = time.monotonic() - ownership_started

    sys.stderr.write(
        "[fill_graphcigartoref_gaps] finalized PA ownership before gap "
        f"analysis: {final_ownership_stats['mapped_pas']} mapped PA(s), "
        f"{final_ownership_stats['resurrected_pas']} resurrected PA(s), "
        f"{final_ownership_stats['contained_resurrected_pas']} contained "
        "resurrected PA(s) removed, "
        f"{final_ownership_stats['suppressed_unmapped_pas']} overlapping "
        "unmapped-only PA(s) removed, "
        f"{final_ownership_stats['ambiguous_mapped_pas']} disconnected "
        "mapped PA(s) removed, "
        f"{final_ownership_stats['unmapped_only_pas']} unmapped-only PA(s) "
        "retained at lowest priority, "
        f"{final_ownership_stats['non_extendable_pas']} alternative/novel "
        "PA(s) retained without extension, "
        f"{final_ownership_stats['unowned_pas']} PA(s) left unowned; wrote "
        f"{updated_pa_count} sparse override row(s) to "
        f"{genomelift_fix_path} in {ownership_elapsed:.2f}s\n"
    )

    discovery_started = time.monotonic()
    represented_deletions = read_represented_deletions(input_path)
    query_reader = core.FastaRegionReader(query_path)
    reference_reader = core.FastaRegionReader(
        reference_path,
        template_fasta=(local_reference_templates or None),
    )
    try:
        reference_lengths = {
            name: values[0] for name, values in reference_reader.index.items()
        }
        tasks, stats = discover_gap_tasks(
            analysis_annotations,
            lift_by_allele,
            query_reader,
            primary_reference_index,
            query_genome,
            args.maximum_gap,
            args.minimum_unmasked,
            args.minimum_unmasked_fraction,
        )
    finally:
        query_reader.handle.close()
        reference_reader.close()
    discovery_elapsed = time.monotonic() - discovery_started

    sys.stderr.write(
        "[fill_graphcigartoref_gaps] additive mode: preserving all "
        f"{len(annotations)} positive-width input row(s); no existing query "
        "or reference alignment will be trimmed or removed\n"
    )
    sys.stderr.write(
        "[fill_graphcigartoref_gaps] corrected ownership analysis view: "
        f"{ownership_stats['logical_pa_rows']} logical PA row(s), "
        f"{ownership_stats['excluded_pa_rows']} disowned PA name(s), "
        f"{ownership_stats['excluded_rows']} row(s) excluded outright, "
        f"{ownership_stats['legacy_rows']} legacy generated row(s), "
        f"{ownership_stats['retained_bases']} retained query bp\n"
    )
    sys.stderr.write(
        "[fill_graphcigartoref_gaps] discovered "
        f"{stats['internal_gaps']} internal gap(s): "
        f"{stats['eligible_gaps']} queued, "
        f"{stats['masked_or_large']} failed size/masking thresholds, "
        f"{stats['missing_priority']} lacked usable neighboring priority, "
        f"{stats['non_primary_reference_boundary']} touched an alternative/"
        "novel target and were skipped, "
        f"{stats['unresolvable_reference_boundary']} had a non-reference-"
        "anchored selected neighbor, "
        f"{stats['zero_reference']} had zero-size same-reference gaps "
        f"({stats['query_insertion']} retained as pure query insertions), "
        f"{stats['unprojectable_tandem']} tandem overlap(s) could not be "
        "projected; "
        f"queued tandem duplications={stats['tandem_duplication']}, "
        f"queued modes bounded={stats['bounded_reference']}, "
        f"one-sided={stats['one_sided_reference']} in "
        f"{discovery_elapsed:.2f}s\n"
    )
    alignment_started = time.monotonic()
    rows, results = align_gap_tasks(
        tasks,
        query_path,
        reference_path,
        local_reference_templates,
        args.processes,
        args.chunksize,
        args.mp_start_method,
        args.anchor,
        args.maximum_reference_extension,
    )
    alignment_elapsed = time.monotonic() - alignment_started
    failures = [result for result in results if result.error]
    for result in failures:
        task = tasks[result.index - 1]
        sys.stderr.write(
            "[fill_graphcigartoref_gaps] warning: failed to align "
            f"{task.name} ({task.contig}:{task.start}-{task.end}):\n"
            f"{result.error}"
        )
    if failures and args.strict:
        raise RuntimeError(f"failed to align {len(failures)} eligible gap(s)")

    deletion_started = time.monotonic()
    generated_annotations = [
        annotation_from_tsv(row, len(annotations) + index + 1)
        for index, row in enumerate(rows)
    ]
    deletion_rows, deletion_stats = discover_reference_gap_deletions(
        tuple(analysis_annotations) + tuple(generated_annotations),
        lift_by_allele,
        primary_reference_index,
        query_genome,
        represented_deletions,
    )
    chrom_lengths = reference_lengths
    rendered_deletions = []
    rendered_deletion_bases = 0
    for row in deletion_rows:
        chrom_length = chrom_lengths[row.ref_coord.chrom]
        rendered = render_synthetic_deletion_or_warn(
            row,
            chrom_length,
            "fill_graphcigartoref_gaps",
        )
        if rendered is not None:
            rendered_deletions.append(rendered)
            clipped_start = max(0, min(row.ref_coord.start, chrom_length))
            clipped_end = max(
                clipped_start, min(row.ref_coord.end, chrom_length),
            )
            rendered_deletion_bases += clipped_end - clipped_start
    base_rows = _all_input_rows(input_path)
    added_rows = list(rows) + rendered_deletions
    _write_output(
        input_path,
        args.output,
        base_rows,
        added_rows,
        args.gaps_only,
    )
    deletion_elapsed = time.monotonic() - deletion_started
    sys.stderr.write(
        f"[fill_graphcigartoref_gaps] emitted {len(rows)}/{len(tasks)} gap row(s)"
        f" with {args.processes} requested process(es); inspected "
        f"{deletion_stats['positive_reference_gaps']} positive reference "
        f"gap(s) ({deletion_stats['reference_gap_bases']} bp), found "
        f"{deletion_stats['reference_gap_covered_bases']} bp represented by "
        "another query alignment; skipped "
        f"{deletion_stats['positive_query_gaps_skipped']} gap(s) "
        f"({deletion_stats['positive_query_gap_bases_skipped']} query bp) "
        "with unresolved query sequence and "
        f"{deletion_stats['query_overlaps_skipped']} pair(s) with query "
        "overlap; emitted "
        f"{len(rendered_deletions)}/"
        f"{deletion_stats['reference_deletions']} uncovered deletion row(s) "
        f"({rendered_deletion_bases} bp); alignment "
        f"{alignment_elapsed:.2f}s, deletion/output "
        f"{deletion_elapsed:.2f}s, total "
        f"{time.monotonic() - run_started:.2f}s\n"
    )
    return 0


def main(argv: Optional[Sequence[str]] = None) -> int:
    try:
        return run(parse_args(argv))
    except Exception as error:
        print(f"fill_graphcigartoref_gaps.py: error: {error}", file=sys.stderr)
        return 1


if __name__ == "__main__":
    raise SystemExit(main())
