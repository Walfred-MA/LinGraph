#!/usr/bin/env python3
"""Build an exclusive pseudo-linear query/reference assignment.

This helper runs at the end of ``fill_graphcigartoref_gaps.py``. It does not alter
the per-line graph alignments. Instead, it elects complete CIGAR mapping
intervals independently on the query and reference axes, then slices each
elected interval from its original CIGAR. Internal matches, mismatches,
insertions, and deletions remain inside that CIGAR and never become election or
report records. Losing query mappings retain their alternative reference path
so duplicated sequence can later be interpreted relative to its source rather
than automatically emitted as a novel insertion. Query gaps and clipped tails
without a supporting alignment CIGAR are reported as UNMAPPED and never gain a
synthetic insertion anchor from neighboring mappings.
Current GenomeLift rows carry separate location and sequence mappings. Their
main and alternative layers are elected independently: a sequence-source
alignment cannot consume query or reference ownership from the location map.
When both cover the same query interval, the sequence layer is reported as a
duplication/source mapping placed on the location interval. The older staged
two-pass policy remains available for graph-CIGAR files without mapping-source
tags. Resurrection remains a query-region tier below extended sequence.
"""

from __future__ import annotations

import argparse
import bisect
import csv
import dataclasses
import functools
import gc
import math
import os
import re
import sys
import time
from collections import defaultdict
from typing import Dict, Iterable, List, Mapping, Optional, Sequence, Tuple

import graphcigartoref as core
from alignment_scoring import alignment_score


# Legacy combined-election order is retained for old unannotated inputs and
# direct callers of ``candidate.priority``. Final GenomeLift annotations use
# the explicit two-pass orders below.
STAGE_PRIORITY = {1: 4, 2: 3, 4: 2, 3: 1, 0: 0}
PRIMARY_STAGE_PRIORITY = {1: 3, 2: 2, 3: 1, 0: 0}
SECONDARY_STAGE_PRIORITY = {2: 3, 3: 2, 4: 1, 0: 0}
REGION_PRIORITY = {
    "owned": 4, "imputed": 3, "extended": 2, "resurrected": 1,
    "trimmed": 0,
}
QUERY_OPS = {"=", "X", "I", "M", "S"}
REFERENCE_OPS = {"=", "X", "D", "M"}
DUP_GROUP_MAX_UNMAPPED = 5_000
DUP_GROUP_MAX_REFERENCE_SPAN = 100_000
OUTER_STRAND_SUFFIX = re.compile(r"\|outer_strand=[+-]$")
LARGE_INSERTION_MIN_QUERY = 500
# Millions of intervals/spans can be live at once; __slots__ removes a per-object
# attribute dict (dataclass slots need Python >= 3.10).
_SLOTS = {"slots": True} if sys.version_info >= (3, 10) else {}

__all__ = (
    "AlignmentInterval", "Lift", "Span", "build_assignments", "elect",
    "interval_assignment_rows", "intervals_from_row", "main",
    "read_intervals", "read_lifts", "write_rows",
)


@dataclasses.dataclass(frozen=True)
class Lift:
    allele: str
    locus: core.Coord
    left: str
    right: str
    stage: int
    secondary: bool
    resurrected: bool = False
    final: bool = False


@dataclasses.dataclass(frozen=True, **_SLOTS)
class _LegacyAtom:
    atom_id: int
    line_number: int
    source_line: str
    query_name: str
    reference_name: str
    query_contig: str
    query_start: int
    query_end: int
    query_strand: str
    reference_path: str
    reference_start: int
    reference_end: int
    reference_strand: str
    op: str
    region: str
    trimmed: bool
    encoded: bool
    lift_stage: int
    secondary_lift: bool
    score: float
    resurrected_lift: bool = False
    anchor_path: str = ""
    anchor_position: int = -1
    # Complete alignment windows from graphcigartoref columns 5/6.  Deletion
    # calls require independent left/right anchors whose extended windows
    # overlap on both axes; the elected atom interval alone cannot establish
    # that the sequence between two blocks was actually assembled.
    query_alignment: str = ""
    reference_alignment: str = ""
    graph_cigar: str = ""
    source_query_interval: str = ""
    alignment_cigar: str = ""
    reference_path_length: int = 0
    # The single CIGAR operation of a main-path atom.  Its anchored CIGAR
    # text is derived only when a partial slice needs it.
    piece_op: Optional[core.CigarOp] = None

    @property
    def region_tier(self) -> str:
        """Return the effective region tier used by the third priority key."""
        return "resurrected" if self.resurrected_lift else self.region

    @property
    def priority(self) -> Tuple[int, int, int, float, int, int]:
        # GenomeLift category is the second key: ordinary assignments retain
        # stage order 1,2,4,3 and dups follow them.  Resurrection belongs to
        # the third (region) key, after extended sequence.
        lift_priority = (
            1 if self.secondary_lift
            else 2 + STAGE_PRIORITY.get(self.lift_stage, 0)
        )
        return (
            0 if self.trimmed else 1,
            lift_priority,
            REGION_PRIORITY[self.region_tier],
            self.score,
            self.query_end - self.query_start + self.reference_end - self.reference_start,
            -self.line_number,
        )


@dataclasses.dataclass(frozen=True, **_SLOTS)
class Span:
    contig: str
    start: int
    end: int
    owner_id: int


@dataclasses.dataclass(frozen=True, **_SLOTS)
class AlignmentInterval:
    """One source alignment interval used by the pseudo-linear election.

    The original graph CIGAR remains intact.  An interval may contain any
    number of match, mismatch, insertion, deletion, or encoded-path
    operations; those operations are never expanded into election records.
    """
    interval_id: int
    line_number: int
    source_line: str
    source_query_name: str
    query_name: str
    reference_name: str
    query_contig: str
    query_start: int
    query_end: int
    query_strand: str
    reference_path: str
    reference_start: int
    reference_end: int
    reference_strand: str
    region: str
    trimmed: bool
    encoded: bool
    lift_stage: int
    secondary_lift: bool
    score: float
    resurrected_lift: bool = False
    anchor_path: str = ""
    anchor_position: int = -1
    query_alignment: str = ""
    reference_alignment: str = ""
    graph_cigar: str = ""
    source_query_interval: str = ""
    reference_path_length: int = 0
    op: str = "ALIGNMENT"
    pa_lifts: Tuple[Lift, ...] = ()
    final_lift_annotation: bool = False
    lift_kind: str = ""
    lift_role: str = ""

    @property
    def region_tier(self) -> str:
        return "resurrected" if self.resurrected_lift else self.region

    @property
    def priority(self) -> Tuple[int, int, int, float, int, int]:
        lift_priority = (
            1 if self.secondary_lift
            else 2 + STAGE_PRIORITY.get(self.lift_stage, 0)
        )
        return (
            0 if self.trimmed else 1,
            lift_priority,
            REGION_PRIORITY[self.region_tier],
            self.score,
            self.query_end - self.query_start
            + self.reference_end - self.reference_start,
            -self.line_number,
        )


def _is_owned_dup_pa(candidate: object) -> bool:
    """Whether this interval is inside a duplicated PA's valid query span.

    A duplicated PA maps an inserted copy back to its source locus.  It must
    remain an alternative-source assignment even when no other row covers
    that source locus; otherwise identical copies can be PRIMARY in one
    sample and full-locus insertions in another solely because their local
    competitors differ.  Extended/imputed sequence outside the PA is not
    included in this rule.
    """
    return bool(
        getattr(candidate, "secondary_lift", False)
        and getattr(candidate, "region", "") == "owned"
    )


def _tag_value(fields: Sequence[str], name: str) -> str:
    prefix = name + ":"
    for field in fields[7:]:
        if not field.startswith(prefix):
            continue
        parts = field.split(":", 2)
        if len(parts) == 3:
            return parts[2].strip()
    return ""


def _row_final_lift_annotation(
    fields: Sequence[str],
) -> Optional[Tuple[int, bool, str, str]]:
    """Return the mapping-specific GenomeLift stage/tier carried by a row."""
    stage_text = _tag_value(fields, "LIFT_STAGE")
    tier = _tag_value(fields, "LIFT_TIER").lower()
    if not stage_text.isdigit() or tier not in {"primary", "secondary"}:
        return None
    kind = _tag_value(fields, "LIFT_KIND").lower()
    role = _tag_value(fields, "LIFT_ROLE").lower()
    if kind not in {"location", "sequence"}:
        kind = ""
    if role not in {"main", "alternative"}:
        role = ""
    return int(stage_text), tier == "secondary", kind, role


def _part_value(value: str, part_index: int) -> str:
    values = value.split(";")
    if part_index and len(values) >= part_index:
        return values[part_index - 1].strip()
    return value.strip()


def _base_query_name(name: str) -> str:
    """Return the GenomeLift allele name for a possibly grouped row member."""
    return OUTER_STRAND_SUFFIX.sub("", name.strip())



def _effective_stage_tier(fields):
    """Columns 15/16 as (stage text, tier): GenomeLift joins its location and
    sequence mappings with ";" (``2;3`` / ``primary;secondary``); the primary
    one is effective, else the first."""
    stages = (fields[14] if len(fields) > 14 else "").strip().split(";")
    tiers = (fields[15] if len(fields) > 15 else "").strip().lower().split(";")
    index = tiers.index("primary") if "primary" in tiers else 0
    return (
        stages[index] if index < len(stages) else stages[0],
        tiers[index] if index < len(tiers) else "",
    )


def read_lifts(paths: Sequence[str]) -> Dict[str, Lift]:
    """Read GenomeLift rows; later sparse files override earlier rows."""
    result: Dict[str, Lift] = {}
    for path in paths:
        if not path:
            continue
        with open(path) as handle:
            for number, raw in enumerate(handle, 1):
                if not raw.strip() or raw.startswith("#"):
                    continue
                fields = raw.rstrip("\r\n").split("\t")
                if fields[0] == "allelename" or fields[0] in {"DEL", "NA"}:
                    continue
                if len(fields) < 14:
                    raise ValueError(f"{path}:{number}: expected at least 14 fields")
                match = re.fullmatch(r"part_(\d+)", fields[6].strip())
                part = int(match.group(1)) if match else 0
                locus = core.parse_coord(_part_value(fields[2], part))
                if locus is None:
                    raise ValueError(f"{path}:{number}: invalid assembly interval")
                previous = result.get(fields[0])
                stage_text, assignment_type = _effective_stage_tier(fields)
                stage = (
                    int(stage_text) if stage_text.isdigit()
                    else previous.stage if previous is not None else 0
                )
                # Grouped GenomeLift rows use part_N in the historical class
                # column.  17-column files retain the PA's original class in
                # column 17; 18-column files hold the location/sequence
                # mappings in 17/18, so the class is column 7.
                class_type = (
                    fields[16].strip().lower()
                    if len(fields) == 17 and fields[16].strip()
                    and "|" not in fields[16]
                    else fields[6].strip().lower()
                )
                # Final stage/tier is authoritative. Historical Pri/Cov/Dup
                # source labels are diagnostic and must not reclassify the
                # mapping after GenomeLift's evidence election.
                final_annotation = assignment_type in {"primary", "secondary"}
                source_is_secondary = (
                    assignment_type == "secondary"
                    if final_annotation else (
                        previous.secondary if previous is not None else False
                    )
                )
                resurrected = (
                    assignment_type == "resurrected"
                    or stage_text.lower() == "resurrected"
                )
                result[fields[0]] = Lift(
                    fields[0], locus, _part_value(fields[12], part),
                    _part_value(fields[13], part), stage,
                    source_is_secondary, resurrected,
                    final_annotation or (
                        previous.final if previous is not None else False
                    ),
                )
    return result


def _offset(token: str, cap: int) -> int:
    token = (token or "0").rsplit(":", 1)[-1].strip()
    if token in {"", "-inf", "+inf", "inf"}:
        return 0
    try:
        value = int(token)
    except ValueError:
        return 0
    return max(-cap, min(cap, value))


def _effective_locus(lift: Lift, cap: int) -> Tuple[int, int]:
    left = _offset(lift.left, cap)
    right = _offset(lift.right, cap)
    # Positive values extend; negative values trim.  Physical left/right are
    # reversed for a minus-strand query interval.
    if lift.locus.strand == "+":
        start = lift.locus.start - left
        end = lift.locus.end + right
    else:
        start = lift.locus.start - right
        end = lift.locus.end + left
    return max(0, min(start, end)), max(start, end)


def _region_for_interval(start: int, end: int, lift: Optional[Lift], cap: int) -> str:
    if lift is None:
        return "extended"
    if start >= lift.locus.start and end <= lift.locus.end:
        return "owned"
    effective_start, effective_end = _effective_locus(lift, cap)
    if start >= effective_start and end <= effective_end:
        return "imputed"
    return "extended"


def _lift_for_interval(
    named_lifts: Sequence[Lift], start: int, end: int,
) -> Optional[Lift]:
    """Choose the grouped PA that owns this physical query interval.

    Exact containment is decisive.  Point operations use half-open point
    containment.  The overlap fallback handles an operation that straddles a
    PA boundary without reverting every atom in the grouped row to its first
    member's metadata.  If there is no overlap at all, the physically nearest
    PA supplies the metadata; PA size is only a deterministic tie-break.
    """
    if not named_lifts:
        return None
    if end > start:
        contained = [
            lift for lift in named_lifts
            if lift.locus.start <= start and lift.locus.end >= end
        ]
    else:
        contained = [
            lift for lift in named_lifts
            if lift.locus.start <= start < lift.locus.end
        ]
        if not contained:
            contained = [
                lift for lift in named_lifts if start == lift.locus.end
            ]
    if contained:
        return min(
            contained,
            key=lambda lift: (lift.locus.end - lift.locus.start, lift.allele),
        )
    overlaps = [
        max(0, min(end, lift.locus.end) - max(start, lift.locus.start))
        for lift in named_lifts
    ]
    if max(overlaps, default=0) > 0:
        return max(
            zip(overlaps, named_lifts),
            key=lambda item: (
                item[0],
                -(item[1].locus.end - item[1].locus.start),
                item[1].allele,
            ),
        )[1]

    def distance(lift: Lift) -> int:
        if end <= lift.locus.start:
            return lift.locus.start - end
        if start >= lift.locus.end:
            return start - lift.locus.end
        return 0

    return min(
        named_lifts,
        key=lambda lift: (
            distance(lift),
            lift.locus.end - lift.locus.start,
            lift.allele,
        ),
    )


def _query_interval_chunks(
    start: int, end: int, named_lifts: Sequence[Lift],
) -> List[Tuple[int, int]]:
    """Split a query-consuming operation at grouped PA boundaries."""
    if end <= start:
        return [(start, end)]
    boundaries = {start, end}
    for lift in named_lifts:
        for boundary in (lift.locus.start, lift.locus.end):
            if start < boundary < end:
                boundaries.add(boundary)
    ordered = sorted(boundaries)
    return list(zip(ordered, ordered[1:]))


def _physical_query(coord: core.Coord, local_start: int, local_end: int) -> Tuple[int, int]:
    if coord.strand == "+":
        return coord.start + local_start, coord.start + local_end
    return coord.end - local_end, coord.end - local_start


def _first_n_gap_query_offset(
    segments: Sequence[core.ReferenceCigarSegment],
) -> Optional[int]:
    """Return the query offset where an N-containing insertion begins.

    Such an operation crosses an assembly gap.  Its anchor and every later
    traversal in the same graph CIGAR depend on unsupported contig adjacency,
    so the pseudo-linear map reports that suffix as UNMAPPED.
    """
    query_cursor = 0
    for segment in segments:
        for operation in segment.ops:
            query_size = operation.n if operation.op in QUERY_OPS else 0
            if (
                operation.op == "I"
                and operation.payload
                and (
                    "N" in operation.payload
                    or "n" in operation.payload
                )
            ):
                return query_cursor
            query_cursor += query_size
    return None


@functools.lru_cache(maxsize=262144)
def _context_intervals(text: str) -> Tuple[Tuple[str, int, int], ...]:
    """Parse a semicolon-separated extended-alignment coordinate column."""
    intervals = []
    for value in text.split(";"):
        coord = core.parse_coord(value.strip())
        if coord is not None and coord.end > coord.start:
            intervals.append((coord.chrom, coord.start, coord.end))
    return tuple(sorted(intervals))


def _context_columns_overlap(left: str, right: str) -> bool:
    """Return whether two coordinate lists have positive same-contig overlap."""
    a = _context_intervals(left)
    b = _context_intervals(right)
    i = j = 0
    while i < len(a) and j < len(b):
        ac, a0, a1 = a[i]
        bc, b0, b1 = b[j]
        if ac < bc:
            i += 1
        elif bc < ac:
            j += 1
        elif a1 <= b0:
            i += 1
        elif b1 <= a0:
            j += 1
        else:
            return True
    return False


def _atoms_have_overlapping_context(left: _LegacyAtom, right: _LegacyAtom) -> bool:
    return (
        _context_columns_overlap(left.query_alignment, right.query_alignment)
        and _context_columns_overlap(
            left.reference_alignment, right.reference_alignment,
        )
    )


def _main_anchor(segment: core.ReferenceCigarSegment) -> Tuple[str, int]:
    direction, path, _length, cursor = segment.main_anchor
    return path, cursor


def _row_score(segments: Sequence[core.ReferenceCigarSegment]) -> float:
    ops = [op for segment in segments for op in segment.ops if op.op != "H"]
    return alignment_score(ops)


def _legacy_atoms_from_row(
    fields: Sequence[str], line_number: int, lifts: Mapping[str, Lift],
    max_impute: int, next_atom_id: int,
) -> Tuple[List[_LegacyAtom], int]:
    if len(fields) < 7:
        raise ValueError(f"line {line_number}: expected seven graphcigartoref columns")
    query_coord = core.parse_coord(fields[2])
    if query_coord is None or query_coord.end <= query_coord.start:
        return [], next_atom_id
    source_line = f"line:{line_number}"
    segments = core.parse_reference_cigar_segments(fields[6], fields[0])
    if not segments:
        return [], next_atom_id
    score = _row_score(segments)
    names = [
        _base_query_name(name) for name in fields[0].split(";")
        if name.strip()
    ]
    named_lifts = [lifts[name] for name in names if name in lifts]
    qcursor = 0
    atoms: List[_LegacyAtom] = []

    for segment in segments:
        rcursor = segment.start if segment.direction == ">" else segment.end
        segment_q_start = qcursor
        segment_query_size = sum(op.n for op in segment.ops if op.op in QUERY_OPS)
        if not segment.is_main:
            q0, q1 = _physical_query(
                query_coord, segment_q_start, segment_q_start + segment_query_size,
            )
            lift = _lift_for_interval(named_lifts, q0, q1)
            anchor_path, anchor_position = _main_anchor(segment)
            atoms.append(_LegacyAtom(
                next_atom_id, line_number, source_line, fields[0], fields[1],
                query_coord.chrom, q0, q1, query_coord.strand,
                segment.path, segment.start, segment.end, segment.direction,
                "ENCODED_INS", _region_for_interval(q0, q1, lift, max_impute),
                False, True, lift.stage if lift else 0,
                lift.secondary if lift else False, score,
                resurrected_lift=lift.resurrected if lift else False,
                anchor_path=anchor_path, anchor_position=anchor_position,
                query_alignment=fields[4], reference_alignment=fields[5],
                graph_cigar=fields[6], source_query_interval=fields[2],
                alignment_cigar=core.format_reference_cigar_segments([segment]),
                reference_path_length=segment.path_len,
            ))
            next_atom_id += 1
            qcursor += segment_query_size
            continue

        for op in segment.ops:
            if op.op == "H" or op.n <= 0:
                continue
            qsize = op.n if op.op in QUERY_OPS else 0
            rsize = op.n if op.op in REFERENCE_OPS else 0
            q0, q1 = _physical_query(query_coord, qcursor, qcursor + qsize)
            for chunk_start, chunk_end in _query_interval_chunks(
                q0, q1, named_lifts,
            ):
                lift = _lift_for_interval(
                    named_lifts, chunk_start, chunk_end,
                )
                if qsize and rsize == qsize:
                    if query_coord.strand == "+":
                        local_start = chunk_start - q0
                        local_end = chunk_end - q0
                    else:
                        local_start = q1 - chunk_end
                        local_end = q1 - chunk_start
                    if segment.direction == ">":
                        r0 = rcursor + local_start
                        r1 = rcursor + local_end
                    else:
                        r0 = rcursor - local_end
                        r1 = rcursor - local_start
                elif rsize:
                    r0, r1 = (
                        (rcursor, rcursor + rsize)
                        if segment.direction == ">"
                        else (rcursor - rsize, rcursor)
                    )
                else:
                    r0 = r1 = rcursor
                if qsize:
                    query_offset = (
                        chunk_start - q0
                        if query_coord.strand == "+" else q1 - chunk_end
                    )
                    piece_size = chunk_end - chunk_start
                else:
                    query_offset = 0
                    piece_size = rsize
                payload = ""
                if op.op == "X" and op.payload:
                    payload = core._slice_x_payload(
                        op.payload, op.n, query_offset, piece_size,
                    )
                elif op.op == "I" and op.payload:
                    payload = op.payload[query_offset:query_offset + piece_size]
                piece_op = core.CigarOp(piece_size, op.op, payload)
                anchor = rcursor
                atoms.append(_LegacyAtom(
                    next_atom_id, line_number, source_line, fields[0], fields[1],
                    query_coord.chrom, chunk_start, chunk_end,
                    query_coord.strand,
                    segment.path, r0, r1, segment.direction, op.op,
                    _region_for_interval(
                        chunk_start, chunk_end, lift, max_impute,
                    ),
                    False, False, lift.stage if lift else 0,
                    lift.secondary if lift else False, score,
                    resurrected_lift=lift.resurrected if lift else False,
                    anchor_path=segment.path, anchor_position=anchor,
                    query_alignment=fields[4], reference_alignment=fields[5],
                    graph_cigar=fields[6], source_query_interval=fields[2],
                    reference_path_length=segment.path_len,
                    piece_op=piece_op,
                ))
                next_atom_id += 1
            qcursor += qsize
            rcursor += rsize if segment.direction == ">" else -rsize

    if qcursor != query_coord.end - query_coord.start:
        raise ValueError(
            f"line {line_number}: CIGAR query span {qcursor} != "
            f"{query_coord.end - query_coord.start}"
        )

    # The complete linear alignment interval may be wider than the final
    # graph-CIGAR row. Preserve those clipped tails as lowest-priority query
    # candidates so the pseudo-linear report marks the query span UNMAPPED.
    # They deliberately carry no CIGAR or reference coordinate and must never
    # acquire a neighboring insertion anchor.
    alignment_coords = [
        core.parse_coord(value.strip())
        for value in fields[4].split(";") if value.strip()
    ]
    for alignment_coord in alignment_coords:
        if alignment_coord is None or alignment_coord.chrom != query_coord.chrom:
            continue
        for start, end in (
            (alignment_coord.start, min(alignment_coord.end, query_coord.start)),
            (max(alignment_coord.start, query_coord.end), alignment_coord.end),
        ):
            if end <= start:
                continue
            lift = _lift_for_interval(named_lifts, start, end)
            atoms.append(_LegacyAtom(
                next_atom_id, line_number, source_line, fields[0], fields[1],
                query_coord.chrom, start, end, alignment_coord.strand,
                "", 0, 0, "+", "TRIMMED", "trimmed", True, False,
                lift.stage if lift else 0, lift.secondary if lift else False,
                score, resurrected_lift=lift.resurrected if lift else False,
                query_alignment=fields[4], reference_alignment=fields[5],
                graph_cigar=fields[6], source_query_interval=fields[2],
                alignment_cigar=".", reference_path_length=0,
            ))
            next_atom_id += 1
    return atoms, next_atom_id


def _legacy_read_atoms(path: str, lifts: Mapping[str, Lift], max_impute: int) -> List[_LegacyAtom]:
    atoms: List[_LegacyAtom] = []
    next_id = 0
    with open(path) as handle:
        for number, raw in enumerate(handle, 1):
            if not raw.strip() or raw.startswith("#"):
                continue
            fields = raw.rstrip("\r\n").split("\t")
            created, next_id = _legacy_atoms_from_row(fields, number, lifts, max_impute, next_id)
            atoms.extend(created)
    return atoms


def _containing_lift(
    named_lifts: Sequence[Lift], start: int, end: int,
) -> Optional[Lift]:
    """Return the PA label containing an interval, or no label outside PAs."""
    if end > start:
        containing = [
            lift for lift in named_lifts
            if lift.locus.start <= start and end <= lift.locus.end
        ]
    else:
        containing = [
            lift for lift in named_lifts
            if lift.locus.start <= start < lift.locus.end
        ]
    if not containing:
        return None
    return min(
        containing,
        key=lambda lift: (lift.locus.end - lift.locus.start, lift.allele),
    )


def _source_query_offsets(
    source: core.Coord, start: int, end: int,
) -> Tuple[int, int]:
    if source.strand == "+":
        return start - source.start, end - source.start
    return source.end - end, source.end - start


def _segments_geometry(
    segments: Sequence[core.ReferenceCigarSegment], row_name: str,
) -> Tuple[str, int, int, str, int, str, int, bool]:
    """Return slice geometry from already parsed reference segments."""
    if not segments:
        raise ValueError(f"{row_name}: empty reference CIGAR interval")
    main = [segment for segment in segments if segment.is_main]
    direction, path, path_length, cursor = segments[0].main_anchor
    if main:
        reference_start = min(segment.start for segment in main)
        reference_end = max(segment.end for segment in main)
        direction = main[0].direction
        path = main[0].path
        path_length = main[0].path_len
        anchor = main[0].start if direction == ">" else main[0].end
    else:
        reference_start = reference_end = cursor
        anchor = cursor
    return (
        path, reference_start, reference_end, direction, path_length,
        path, anchor, any(not segment.is_main for segment in segments),
    )


def intervals_from_row(
    fields: Sequence[str], line_number: int, lifts: Mapping[str, Lift],
    max_impute: int, next_interval_id: int,
) -> Tuple[List[AlignmentInterval], int]:
    """Create interval candidates without expanding the CIGAR operations."""
    if len(fields) < 7:
        raise ValueError(f"line {line_number}: expected seven graphcigartoref columns")
    query_coord = core.parse_coord(fields[2])
    if query_coord is None or query_coord.end <= query_coord.start:
        return [], next_interval_id
    source_segments = core.parse_reference_cigar_segments(fields[6], fields[0])
    if not source_segments:
        return [], next_interval_id
    source_query_span = sum(
        op.n for segment in source_segments for op in segment.ops
        if op.op in QUERY_OPS
    )
    if source_query_span != query_coord.end - query_coord.start:
        raise ValueError(
            f"line {line_number}: CIGAR query span {source_query_span} != "
            f"{query_coord.end - query_coord.start}"
        )
    n_gap_query_offset = _first_n_gap_query_offset(source_segments)
    score = _row_score(source_segments)
    names = [
        _base_query_name(name) for name in fields[0].split(";")
        if name.strip()
    ]
    named_lifts = [lifts[name] for name in names if name in lifts]
    row_lift = _row_final_lift_annotation(fields)
    uniform_named_lift = None
    if named_lifts and len({
        (lift.stage, lift.secondary, lift.final) for lift in named_lifts
    }) == 1:
        uniform_named_lift = named_lifts[0]

    def stage_for(lift: Optional[Lift]) -> int:
        # The current coord map can repair a graph-CIGAR generated from an
        # older GenomeLift file. The row annotation covers extended sequence
        # outside a named PA and is otherwise expected to agree.
        if row_lift is not None and row_lift[2]:
            return row_lift[0]
        if lift is not None:
            return lift.stage
        if uniform_named_lift is not None:
            return uniform_named_lift.stage
        return row_lift[0] if row_lift is not None else 0

    def secondary_for(lift: Optional[Lift]) -> bool:
        stage = stage_for(lift)
        if row_lift is not None and row_lift[3]:
            # The role names the layer; the tier keeps GenomeLift's election
            # within it: a location duplicate is secondary as well.
            return row_lift[3] == "alternative" or row_lift[1]
        if lift is not None:
            secondary = lift.secondary
            final = lift.final
        elif uniform_named_lift is not None:
            secondary = uniform_named_lift.secondary
            final = uniform_named_lift.final
        else:
            secondary = row_lift[1] if row_lift is not None else False
            final = row_lift is not None
        # The pass policy fixes the tiers at the two endpoints.  Stages 2/3
        # use GenomeLift's elected tier, while stage 1 is always primary and
        # stage 4 is always an alternative/duplicate mapping.
        if final and stage == 1:
            return False
        if final and stage == 4:
            return True
        return secondary

    source_line = f"line:{line_number}"
    intervals: List[AlignmentInterval] = []

    qcursor = 0
    for source_segment in source_segments:
        source_segment_query_size = sum(
            op.n for op in source_segment.ops if op.op in QUERY_OPS
        )
        if (
            n_gap_query_offset is not None
            and qcursor >= n_gap_query_offset
        ):
            break
        segment_query_size = source_segment_query_size
        segment = source_segment
        if (
            n_gap_query_offset is not None
            and qcursor + segment_query_size > n_gap_query_offset
        ):
            segment_query_size = n_gap_query_offset - qcursor
            if segment_query_size <= 0:
                break
            sliced_segment = core.slice_query_segment_by_query(
                source_segment, 0, segment_query_size,
                include_left_boundary_deletions=True,
                include_right_boundary_deletions=True,
            )
            parsed_prefix = core.parse_reference_cigar_segments(
                sliced_segment, fields[0],
            )
            if not parsed_prefix:
                break
            segment = parsed_prefix[0]
        segment_query_size = sum(
            op.n for op in segment.ops if op.op in QUERY_OPS
        )
        if segment_query_size <= 0:
            continue
        segment_start, segment_end = _physical_query(
            query_coord, qcursor, qcursor + segment_query_size,
        )
        # A main-path segment has a direct query/reference interval map, so a
        # PA boundary may cut that interval. An encoded-path traversal is
        # already one indivisible mapping interval and is never operation-
        # split merely because its query span crosses a PA boundary.
        chunks = (
            _query_interval_chunks(segment_start, segment_end, named_lifts)
            if segment.is_main else [(segment_start, segment_end)]
        )
        chunk_slices = None
        if segment.is_main and len(chunks) > 1:
            segment_cigar = core.format_reference_cigar_segments([segment])
            temporary_interval = AlignmentInterval(
                interval_id=-1, line_number=line_number,
                source_line=source_line, source_query_name=fields[0],
                query_name="", reference_name=fields[1],
                query_contig=query_coord.chrom,
                query_start=segment_start, query_end=segment_end,
                query_strand=query_coord.strand,
                reference_path=segment.path,
                reference_start=segment.start, reference_end=segment.end,
                reference_strand=segment.direction, region="extended",
                trimmed=False, encoded=False, lift_stage=0,
                secondary_lift=False, score=score,
                graph_cigar=segment_cigar,
                source_query_interval=(
                    f"{query_coord.chrom}:{segment_start}-{segment_end}"
                    f"{query_coord.strand}"
                ),
                reference_path_length=segment.path_len,
            )
            chunk_slices, _unused_deletions, _unused_full_slices = _slice_interval_ranges(
                temporary_interval, chunks, [segment],
            )
        for query_start, query_end in chunks:
            global_local_start, global_local_end = _source_query_offsets(
                query_coord, query_start, query_end,
            )
            if segment.is_main:
                rel_start = global_local_start - qcursor
                rel_end = global_local_end - qcursor
                if chunk_slices is not None:
                    interval_cigar, geometry = chunk_slices[
                        (query_start, query_end)
                    ]
                    (
                        reference_path, reference_start, reference_end,
                        reference_strand, reference_path_length,
                        anchor_path, anchor_position, encoded,
                    ) = geometry
                    interval_segment = None
                elif rel_start == 0 and rel_end == segment_query_size:
                    part = segment
                else:
                    sliced_segment = core.slice_query_segment_by_query(
                        segment, rel_start, rel_end,
                        include_left_boundary_deletions=True,
                        include_right_boundary_deletions=(
                            rel_end == segment_query_size
                        ),
                    )
                    parsed_parts = core.parse_graphic_segments(
                        sliced_segment, fields[0],
                    )
                    if not parsed_parts:
                        raise ValueError(
                            f"line {line_number}: empty main interval slice"
                        )
                    part = parsed_parts[0]
                if chunk_slices is None:
                    reference_path = part.path
                    reference_start, reference_end = part.start, part.end
                    reference_strand = part.direction
                    reference_path_length = part.path_len
                    anchor_path = part.path
                    anchor_position = (
                        part.start if part.direction == ">" else part.end
                    )
                    encoded = False
                    interval_segment = dataclasses.replace(
                        segment, start=part.start, end=part.end,
                        cigar=part.cigar, ops=list(part.ops),
                    )
            else:
                (
                    reference_strand, reference_path,
                    reference_path_length, anchor_position,
                ) = segment.main_anchor
                reference_start = reference_end = anchor_position
                anchor_path = reference_path
                encoded = True
                interval_segment = segment
            if interval_segment is not None:
                interval_cigar = core.format_reference_cigar_segments(
                    [interval_segment],
                )
            source_lift = _lift_for_interval(
                named_lifts, query_start, query_end,
            )
            intervals.append(AlignmentInterval(
                next_interval_id, line_number, source_line, fields[0],
                "", fields[1], query_coord.chrom,
                query_start, query_end, query_coord.strand,
                reference_path, reference_start, reference_end,
                reference_strand,
                _region_for_interval(
                    query_start, query_end, source_lift, max_impute,
                ),
                False, encoded, stage_for(source_lift),
                secondary_for(source_lift), score,
                resurrected_lift=(
                    source_lift.resurrected if source_lift else False
                ),
                anchor_path=anchor_path, anchor_position=anchor_position,
                query_alignment=fields[4], reference_alignment=fields[5],
                graph_cigar=interval_cigar,
                source_query_interval=(
                    f"{query_coord.chrom}:{query_start}-{query_end}"
                    f"{query_coord.strand}"
                ),
                reference_path_length=reference_path_length,
                pa_lifts=tuple(named_lifts),
                final_lift_annotation=(
                    row_lift is not None
                    or bool(source_lift and source_lift.final)
                    or bool(uniform_named_lift and uniform_named_lift.final)
                ),
                lift_kind=(row_lift[2] if row_lift is not None else ""),
                lift_role=(row_lift[3] if row_lift is not None else ""),
            ))
            next_interval_id += 1
        qcursor += segment_query_size

    if n_gap_query_offset is not None:
        gap_start, gap_end = _physical_query(
            query_coord, n_gap_query_offset, source_query_span,
        )
        gap_lift = _lift_for_interval(named_lifts, gap_start, gap_end)
        intervals.append(AlignmentInterval(
            next_interval_id, line_number, source_line, fields[0],
            "", fields[1], query_coord.chrom, gap_start, gap_end,
            query_coord.strand, "", 0, 0, ">", "trimmed",
            True, False, stage_for(gap_lift),
            secondary_for(gap_lift), score,
            resurrected_lift=(
                gap_lift.resurrected if gap_lift else False
            ),
            query_alignment=fields[4], reference_alignment=fields[5],
            graph_cigar="", source_query_interval=(
                f"{query_coord.chrom}:{gap_start}-{gap_end}"
                f"{query_coord.strand}"
            ),
            reference_path_length=0,
            pa_lifts=tuple(named_lifts),
            final_lift_annotation=(
                row_lift is not None
                or bool(gap_lift and gap_lift.final)
                or bool(uniform_named_lift and uniform_named_lift.final)
            ),
            lift_kind=(row_lift[2] if row_lift is not None else ""),
            lift_role=(row_lift[3] if row_lift is not None else ""),
        ))
        next_interval_id += 1

    # Preserve clipped alignment tails as unanchored, lowest-priority query
    # intervals. They become UNMAPPED report rows; no synthetic insertion
    # CIGAR or neighboring reference anchor is created for them.
    alignment_coords = [] if (
        row_lift is not None and row_lift[2] == "sequence"
    ) else [
        core.parse_coord(value.strip())
        for value in fields[4].split(";") if value.strip()
    ]
    for alignment_coord in alignment_coords:
        if alignment_coord is None or alignment_coord.chrom != query_coord.chrom:
            continue
        for tail_start, tail_end in (
            (alignment_coord.start, min(alignment_coord.end, query_coord.start)),
            (max(alignment_coord.start, query_coord.end), alignment_coord.end),
        ):
            if tail_end <= tail_start:
                continue
            tail_lift = _lift_for_interval(named_lifts, tail_start, tail_end)
            intervals.append(AlignmentInterval(
                next_interval_id, line_number, source_line, fields[0],
                "", fields[1], alignment_coord.chrom, tail_start, tail_end,
                alignment_coord.strand, "", 0, 0, ">", "trimmed",
                True, False, stage_for(tail_lift),
                secondary_for(tail_lift), score,
                resurrected_lift=tail_lift.resurrected if tail_lift else False,
                query_alignment=fields[4], reference_alignment=fields[5],
                source_query_interval=fields[2],
                pa_lifts=tuple(named_lifts),
                final_lift_annotation=(
                    row_lift is not None
                    or bool(tail_lift and tail_lift.final)
                    or bool(uniform_named_lift and uniform_named_lift.final)
                ),
                lift_kind=(row_lift[2] if row_lift is not None else ""),
                lift_role=(row_lift[3] if row_lift is not None else ""),
            ))
            next_interval_id += 1
    return intervals, next_interval_id


def read_intervals(
    path: str, lifts: Mapping[str, Lift], max_impute: int,
) -> List[AlignmentInterval]:
    intervals: List[AlignmentInterval] = []
    next_id = 0
    with open(path) as handle:
        for number, raw in enumerate(handle, 1):
            if not raw.strip() or raw.startswith("#"):
                continue
            fields = raw.rstrip("\r\n").split("\t")
            created, next_id = intervals_from_row(
                fields, number, lifts, max_impute, next_id,
            )
            intervals.extend(created)
    return intervals


def _pass_priority(candidate: object, secondary_pass: bool) -> tuple:
    stages = (
        SECONDARY_STAGE_PRIORITY if secondary_pass
        else PRIMARY_STAGE_PRIORITY
    )
    return (
        0 if getattr(candidate, "trimmed", False) else 1,
        stages.get(getattr(candidate, "lift_stage", 0), 0),
        REGION_PRIORITY[getattr(candidate, "region_tier")],
        getattr(candidate, "score", 0.0),
        getattr(candidate, "query_end") - getattr(candidate, "query_start")
        + getattr(candidate, "reference_end")
        - getattr(candidate, "reference_start"),
        -getattr(candidate, "line_number"),
    )


def elect(
    spans: Sequence[Span], candidates: Sequence[object],
    priority_key=None,
) -> List[Span]:
    """Return a non-overlapping winner sweep for one coordinate axis.

    Boundaries are coordinate-compressed once. Spans are considered in winner
    order, and a disjoint-set ``next unassigned`` index visits each elementary
    interval only once. After sorting, election is effectively linear rather
    than scanning every span at every boundary.
    """
    by_contig: Dict[str, List[Span]] = defaultdict(list)
    for span in spans:
        if span.end > span.start:
            by_contig[span.contig].append(span)
    output: List[Span] = []
    for contig, local in sorted(by_contig.items()):
        boundaries = sorted({x for span in local for x in (span.start, span.end)})
        boundary_index = {
            coordinate: index for index, coordinate in enumerate(boundaries)
        }
        cell_count = max(0, len(boundaries) - 1)
        winners: List[Optional[Span]] = [None] * cell_count
        parent = list(range(cell_count + 1))

        def find(index: int) -> int:
            root = index
            while parent[root] != root:
                root = parent[root]
            while parent[index] != index:
                following = parent[index]
                parent[index] = root
                index = following
            return root

        ranked = sorted(
            enumerate(local),
            key=lambda item: (
                (
                    priority_key(candidates[item[1].owner_id])
                    if priority_key is not None
                    else candidates[item[1].owner_id].priority
                ),
                -item[0],
            ),
            reverse=True,
        )
        for _order, span in ranked:
            cell = find(boundary_index[span.start])
            stop = boundary_index[span.end]
            while cell < stop:
                winners[cell] = span
                parent[cell] = find(cell + 1)
                cell = parent[cell]

        for index, winner in enumerate(winners):
            if winner is None:
                continue
            start = boundaries[index]
            end = boundaries[index + 1]
            if output and output[-1].contig == contig and output[-1].end == start and output[-1].owner_id == winner.owner_id:
                output[-1] = Span(contig, output[-1].start, end, winner.owner_id)
            else:
                output.append(Span(contig, start, end, winner.owner_id))
    return output


def _overlaps_for_projected_spans(
    reference_spans: Sequence[Span], requests: Sequence[Tuple[int, str, int, int]],
) -> Dict[int, List[Span]]:
    """Sweep non-overlapping reference winners against projected queries."""
    references_by_contig: Dict[str, List[Span]] = defaultdict(list)
    requests_by_contig: Dict[str, List[Tuple[int, int, int]]] = defaultdict(list)
    for span in reference_spans:
        references_by_contig[span.contig].append(span)
    for request_id, contig, start, end in requests:
        if end > start:
            requests_by_contig[contig].append((start, end, request_id))

    result: Dict[int, List[Span]] = defaultdict(list)
    for contig, local_requests in requests_by_contig.items():
        local_references = references_by_contig.get(contig, ())
        local_requests.sort(key=lambda item: (item[0], item[1], item[2]))
        left = 0
        for start, end, request_id in local_requests:
            while left < len(local_references) and local_references[left].end <= start:
                left += 1
            index = left
            while index < len(local_references) and local_references[index].start < end:
                if local_references[index].end > start:
                    result[request_id].append(local_references[index])
                index += 1
    return result


def _project_query(atom: _LegacyAtom, qstart: int, qend: int) -> Tuple[int, int]:
    if atom.query_end <= atom.query_start or atom.reference_end <= atom.reference_start:
        return atom.anchor_position, atom.anchor_position
    same = _same_mapping_orientation(atom)
    if same:
        return (
            atom.reference_start + qstart - atom.query_start,
            atom.reference_start + qend - atom.query_start,
        )
    return (
        atom.reference_end - (qend - atom.query_start),
        atom.reference_end - (qstart - atom.query_start),
    )


def _project_reference(atom: _LegacyAtom, rstart: int, rend: int) -> Tuple[int, int]:
    if atom.query_end <= atom.query_start or atom.reference_end <= atom.reference_start:
        return atom.query_start, atom.query_start
    same = _same_mapping_orientation(atom)
    if same:
        return (
            atom.query_start + rstart - atom.reference_start,
            atom.query_start + rend - atom.reference_start,
        )
    return (
        atom.query_end - (rend - atom.reference_start),
        atom.query_end - (rstart - atom.reference_start),
    )


def _same_mapping_orientation(atom: object) -> bool:
    """Whether increasing physical query coordinates increase reference coordinates."""
    return (atom.query_strand == "+") == (atom.reference_strand == ">")


def _row_owner(row: Mapping[str, object]):
    """Return the interval owner (or a legacy test owner)."""
    owner = row.get("interval")
    return owner if owner is not None else row["atom"]


def _reference_edge_at_query_side(row: dict, side: str) -> int:
    """Return the reference coordinate at a row's physical query edge."""
    same = _same_mapping_orientation(_row_owner(row))
    if side == "right":
        return row["reference_end"] if same else row["reference_start"]
    if side == "left":
        return row["reference_start"] if same else row["reference_end"]
    raise ValueError(f"invalid physical query side: {side}")


def _bracketed_deletion_intervals(rows: Sequence[dict]) -> Dict[
    Tuple[str, str, str], List[Tuple[int, int]]
]:
    """Find reference gaps supported by two neighboring primary anchors.

    A reference absence is callable only when the physical query has no
    unresolved bases between its anchors and both anchors come from extended
    alignment windows that overlap on the query and reference axes.  This is
    the evidence that distinguishes an assembled deletion from sequence that
    is merely unmapped or cut off at an alignment boundary.
    """
    by_query: Dict[str, List[dict]] = defaultdict(list)
    for row in rows:
        if isinstance(row.get("query_start"), int):
            by_query[row["query_contig"]].append(row)
    supported: Dict[Tuple[str, str, str], List[Tuple[int, int]]] = defaultdict(list)
    for contig, local in by_query.items():
        local.sort(key=lambda row: (
            row["query_start"], row["query_end"], row["kind"],
        ))
        # Walk once along the elected query intervals.  A literal/encoded
        # insertion from the same source alignment accounts for query bases
        # between two anchors and therefore must not hide an adjacent
        # reference deletion (for example 10=5D3I10=).  Unresolved,
        # trimmed, cross-source, or discontinuous query sequence still breaks
        # the bracket.
        left = None
        bridge_end = -1
        bridge_source = ""
        for right in local:
            if right["kind"] == "INSERTION":
                atom = right["atom"]
                accounted = (
                    left is not None
                    and right["query_start"] == bridge_end
                    and atom.source_line == bridge_source
                    and not atom.trimmed
                    and right.get("anchor", ".") != "."
                )
                if accounted:
                    bridge_end = right["query_end"]
                else:
                    left = None
                    bridge_end = -1
                    bridge_source = ""
                continue
            if right["kind"] != "PRIMARY":
                left = None
                bridge_end = -1
                bridge_source = ""
                continue
            if left is None:
                left = right
                bridge_end = right["query_end"]
                bridge_source = right["atom"].source_line
                continue
            if bridge_end != right["query_start"]:
                left = right
                bridge_end = right["query_end"]
                bridge_source = right["atom"].source_line
                continue
            # When an insertion bridged the query interval, it and both
            # anchors must come from the same graph-CIGAR row.  Directly
            # adjacent anchors may continue to come from neighboring rows.
            if (
                bridge_end != left["query_end"]
                and right["atom"].source_line != bridge_source
            ):
                left = right
                bridge_end = right["query_end"]
                bridge_source = right["atom"].source_line
                continue
            if (
                left["query_strand"] != right["query_strand"]
                or left["reference_path"] != right["reference_path"]
                or left["reference_strand"] != right["reference_strand"]
            ):
                left = right
                bridge_end = right["query_end"]
                bridge_source = right["atom"].source_line
                continue
            left_atom = left["atom"]
            right_atom = right["atom"]
            if _atoms_have_overlapping_context(left_atom, right_atom):
                left_edge = _reference_edge_at_query_side(left, "right")
                right_edge = _reference_edge_at_query_side(right, "left")
                same = _same_mapping_orientation(left_atom)
                if (same and right_edge > left_edge) or (
                    not same and right_edge < left_edge
                ):
                    supported[(
                        contig, left["reference_path"], left["reference_strand"],
                    )].append((
                        min(left_edge, right_edge), max(left_edge, right_edge),
                    ))
            left = right
            bridge_end = right["query_end"]
            bridge_source = right["atom"].source_line
    return supported


def _retain_bracketed_deletions(
    rows: Sequence[dict], candidates: Sequence[dict],
) -> List[dict]:
    """Filter deletion candidates by supported intervals in sorted sweeps."""
    supported = _bracketed_deletion_intervals(rows)
    indexed: Dict[Tuple[str, str, str], List[Tuple[int, dict]]] = defaultdict(list)
    for index, row in enumerate(candidates):
        indexed[(
            row["query_contig"], row["reference_path"],
            row["reference_strand"],
        )].append((index, row))
    keep = set()
    for key, local_candidates in indexed.items():
        intervals = sorted(supported.get(key, ()))
        if not intervals:
            continue
        local_candidates.sort(key=lambda item: (
            item[1]["reference_start"], item[1]["reference_end"], item[0],
        ))
        interval_index = 0
        furthest_end = -1
        for original_index, row in local_candidates:
            start = row["reference_start"]
            while (
                interval_index < len(intervals)
                and intervals[interval_index][0] <= start
            ):
                furthest_end = max(furthest_end, intervals[interval_index][1])
                interval_index += 1
            if furthest_end >= row["reference_end"]:
                keep.add(original_index)
    return [row for index, row in enumerate(candidates) if index in keep]


def _same_public_assignment_metadata(left: _LegacyAtom, right: _LegacyAtom) -> bool:
    return (
        left.source_line == right.source_line
        and left.query_name == right.query_name
        and left.reference_name == right.reference_name
        and left.region == right.region
        and left.trimmed == right.trimmed
        and left.encoded == right.encoded
        and left.lift_stage == right.lift_stage
        and left.secondary_lift == right.secondary_lift
        and left.resurrected_lift == right.resurrected_lift
        and left.score == right.score
        and left.graph_cigar == right.graph_cigar
    )


def _legacy_coalesce_primary_rows(rows: Sequence[dict]) -> List[dict]:
    """Collapse adjacent atomic matches into auditable alignment intervals.

    Compatibility path for callers of ``assignment_rows``. The production
    builder elects AlignmentInterval objects directly and does not call this.
    """
    output: List[dict] = []
    old_to_new: Dict[int, int] = {}
    for old_index, original in enumerate(rows):
        row = dict(original)
        if "_alignment_pieces" in row:
            row["_alignment_pieces"] = list(row["_alignment_pieces"])
        if not output:
            output.append(row)
            old_to_new[old_index] = 0
            continue
        previous = output[-1]
        left_atom = previous["atom"]
        right_atom = row["atom"]
        same_orientation = _same_mapping_orientation(left_atom)
        reference_contiguous = (
            previous["reference_end"] == row["reference_start"]
            if same_orientation
            else previous["reference_start"] == row["reference_end"]
        )
        if (
            previous["kind"] == row["kind"] == "PRIMARY"
            and previous["query_contig"] == row["query_contig"]
            and previous["query_end"] == row["query_start"]
            and previous["query_strand"] == row["query_strand"]
            and previous["reference_path"] == row["reference_path"]
            and previous["reference_strand"] == row["reference_strand"]
            and reference_contiguous
            and _same_public_assignment_metadata(left_atom, right_atom)
        ):
            previous["query_end"] = row["query_end"]
            previous["reference_start"] = min(
                previous["reference_start"], row["reference_start"],
            )
            previous["reference_end"] = max(
                previous["reference_end"], row["reference_end"],
            )
            previous_operation = previous.get("operation", left_atom.op)
            previous["operation"] = (
                previous_operation
                if previous_operation == right_atom.op
                else "ALIGNMENT"
            )
            previous.setdefault("_alignment_pieces", [left_atom]).extend(
                row.get("_alignment_pieces", [(right_atom, row["query_start"], row["query_end"])])
            )
            old_to_new[old_index] = len(output) - 1
        else:
            output.append(row)
            old_to_new[old_index] = len(output) - 1
    # Unresolved insertion chains use temporary @<zero-based-row-index>
    # references.  Consolidation only merges PRIMARY rows, but remap every
    # dependency explicitly so the representation remains correct even if a
    # future producer points at a mergeable row.
    for row in output:
        dependency = row.get("dependency", ".")
        if isinstance(dependency, str) and dependency.startswith("@"):
            old_target = int(dependency[1:])
            if old_target not in old_to_new:
                raise ValueError(
                    f"pseudo-linear dependency references missing row {old_target}"
                )
            row["dependency"] = f"@{old_to_new[old_target]}"
    return output


class _ReferenceSpanIndex:
    """Indexed access to the non-overlapping reference-election result."""

    def __init__(self, spans: Sequence[Span]):
        self.by_path = {}
        grouped = defaultdict(list)
        for span in spans:
            grouped[span.contig].append(span)
        for path, local in grouped.items():
            local.sort(key=lambda item: (item.start, item.end))
            self.by_path[path] = (
                local,
                [item.start for item in local],
            )

    def overlaps(self, path: str, start: int, end: int) -> List[Span]:
        indexed = self.by_path.get(path)
        if indexed is None or end <= start:
            return []
        spans, starts = indexed
        index = max(0, bisect.bisect_right(starts, start) - 1)
        while index < len(spans) and spans[index].end <= start:
            index += 1
        output = []
        while index < len(spans) and spans[index].start < end:
            if spans[index].end > start:
                output.append(spans[index])
            index += 1
        return output

    def owned_parts(
        self, path: str, start: int, end: int, owner_id: int,
    ) -> List[Tuple[int, int]]:
        output = []
        for span in self.overlaps(path, start, end):
            if span.owner_id != owner_id:
                continue
            local_start = max(start, span.start)
            local_end = min(end, span.end)
            if local_end <= local_start:
                continue
            if output and output[-1][1] == local_start:
                output[-1] = (output[-1][0], local_end)
            else:
                output.append((local_start, local_end))
        return output

    def owns(
        self, path: str, start: int, end: int, owner_id: int,
    ) -> bool:
        return self.owned_parts(path, start, end, owner_id) == [(start, end)]


def _physical_query_point(source: core.Coord, local: int) -> int:
    return source.start + local if source.strand == "+" else source.end - local


def _slice_interval_ranges(
    interval: AlignmentInterval,
    ranges: Sequence[Tuple[int, int]],
    segments: Sequence[core.ReferenceCigarSegment],
    reference_index: Optional[_ReferenceSpanIndex] = None,
) -> Tuple[
    Dict[Tuple[int, int], Tuple[
        str, Tuple[str, int, int, str, int, str, int, bool]
    ]],
    List[Tuple[int, str, int, int, str, int]],
    Dict[Tuple[int, int], Tuple[
        str, Tuple[str, int, int, str, int, str, int, bool]
    ]],
]:
    """Slice all requested ranges in one CIGAR walk.

    Query-consuming operations are distributed over the sorted, disjoint
    ranges. A deletion at a query cut is assigned once to the following range
    in CIGAR order. If reference election divides a deletion, its owned
    reference pieces are returned as zero-query deletion intervals instead of
    dropping the complete operation.  The third return value retains the full
    source alignment for alternative-source calling.  A deletion that another
    row owns must not be called on the elected main path, but it still consumes
    source-reference coordinates in that alternative alignment.
    """
    source = core.parse_coord(interval.source_query_interval)
    if source is None:
        raise ValueError(
            f"{interval.source_line}: missing source query interval"
        )
    specs = []
    for key in dict.fromkeys(ranges):
        query_start, query_end = key
        if query_end <= query_start:
            continue
        local_start, local_end = _source_query_offsets(
            source, query_start, query_end,
        )
        specs.append({
            "key": key, "start": local_start, "end": local_end,
            "segments": {},
        })
    specs.sort(key=lambda item: (item["start"], item["end"]))
    if not specs:
        return {}, [], {}
    for left, right in zip(specs, specs[1:]):
        if left["end"] > right["start"]:
            raise ValueError(
                f"{interval.source_line}: overlapping interval slice requests"
            )
    starts = [item["start"] for item in specs]
    ends = [item["end"] for item in specs]
    deletion_intervals = []

    def builder(spec, segment_index, segment, anchor):
        result = spec["segments"].get(segment_index)
        if result is None:
            result = {
                "template": segment, "ops": [], "ref_min": None,
                "ref_max": None, "anchor": anchor,
                "full_ops": [], "full_ref_min": None,
                "full_ref_max": None,
            }
            spec["segments"][segment_index] = result
        return result

    def mark_ref(result, start, end):
        if end < start:
            start, end = end, start
        if end <= start:
            return
        result["ref_min"] = (
            start if result["ref_min"] is None
            else min(result["ref_min"], start)
        )
        result["ref_max"] = (
            end if result["ref_max"] is None
            else max(result["ref_max"], end)
        )

    def mark_full_ref(result, start, end):
        if end < start:
            start, end = end, start
        if end <= start:
            return
        result["full_ref_min"] = (
            start if result["full_ref_min"] is None
            else min(result["full_ref_min"], start)
        )
        result["full_ref_max"] = (
            end if result["full_ref_max"] is None
            else max(result["full_ref_max"], end)
        )

    def range_at_point(query_point):
        index = bisect.bisect_left(starts, query_point)
        if (
            index < len(specs)
            and starts[index] == query_point < ends[index]
        ):
            return index
        if index and starts[index - 1] < query_point <= ends[index - 1]:
            return index - 1
        return None

    query_cursor = 0
    for segment_index, segment in enumerate(segments):
        reference_cursor = (
            segment.start if segment.direction == ">" else segment.end
        )
        for operation in segment.ops:
            if operation.op == "H" or operation.n <= 0:
                continue
            query_size = operation.n if operation.op in QUERY_OPS else 0
            reference_size = (
                operation.n if operation.op in REFERENCE_OPS else 0
            )
            operation_start = query_cursor
            operation_end = query_cursor + query_size
            if query_size:
                range_index = bisect.bisect_right(ends, operation_start)
                while (
                    range_index < len(specs)
                    and starts[range_index] < operation_end
                ):
                    overlap_start = max(
                        operation_start, starts[range_index],
                    )
                    overlap_end = min(operation_end, ends[range_index])
                    if overlap_end > overlap_start:
                        offset = overlap_start - operation_start
                        take = overlap_end - overlap_start
                        if operation.op in {"=", "X", "M"}:
                            if segment.direction == "<":
                                local_reference_end = reference_cursor - offset
                                local_reference_start = (
                                    local_reference_end - take
                                )
                            else:
                                local_reference_start = (
                                    reference_cursor + offset
                                )
                                local_reference_end = (
                                    local_reference_start + take
                                )
                        else:
                            local_reference_start = reference_cursor
                            local_reference_end = reference_cursor
                        result = builder(
                            specs[range_index], segment_index, segment,
                            local_reference_start,
                        )
                        mark_ref(
                            result, local_reference_start,
                            local_reference_end,
                        )
                        mark_full_ref(
                            result, local_reference_start,
                            local_reference_end,
                        )
                        payload = ""
                        if (
                            operation.op == "I"
                            and operation.payload
                            and len(operation.payload) == operation.n
                        ):
                            payload = operation.payload[offset:offset + take]
                        elif operation.op == "X" and operation.payload:
                            payload = core._slice_x_payload(
                                operation.payload, operation.n, offset, take,
                            )
                        sliced_operation = core.CigarOp(
                            take, operation.op, payload,
                        )
                        result["ops"].append(sliced_operation)
                        result["full_ops"].append(sliced_operation)
                    range_index += 1
            elif operation.op == "D" and reference_size:
                if segment.direction == "<":
                    deletion_start = reference_cursor - reference_size
                    deletion_end = reference_cursor
                else:
                    deletion_start = reference_cursor
                    deletion_end = reference_cursor + reference_size
                range_index = range_at_point(query_cursor)
                owned = [(deletion_start, deletion_end)]
                if reference_index is not None and segment.is_main:
                    owned = reference_index.owned_parts(
                        segment.path, deletion_start, deletion_end,
                        interval.interval_id,
                    )
                if owned == [(deletion_start, deletion_end)]:
                    if range_index is not None:
                        result = builder(
                            specs[range_index], segment_index, segment,
                            deletion_start,
                        )
                        mark_ref(result, deletion_start, deletion_end)
                        result["ops"].append(operation)
                elif range_index is not None:
                    physical = _physical_query_point(source, query_cursor)
                    for owned_start, owned_end in owned:
                        deletion_intervals.append((
                            physical, segment.path, owned_start, owned_end,
                            segment.direction, segment.path_len,
                        ))
                if range_index is not None:
                    result = builder(
                        specs[range_index], segment_index, segment,
                        deletion_start,
                    )
                    mark_full_ref(result, deletion_start, deletion_end)
                    result["full_ops"].append(operation)
            query_cursor = operation_end
            if segment.direction == "<":
                reference_cursor -= reference_size
            else:
                reference_cursor += reference_size

    output = {}
    full_output = {}
    for spec in specs:
        sliced_segments = []
        for segment_index in sorted(spec["segments"]):
            result = spec["segments"][segment_index]
            if not result["ops"]:
                continue
            start = result["ref_min"]
            end = result["ref_max"]
            if start is None or end is None:
                start = end = result["anchor"]
            sliced_segments.append(dataclasses.replace(
                result["template"], start=start, end=end, cigar="",
                ops=result["ops"],
            ))
        cigar = core.format_reference_cigar_segments(sliced_segments)
        if not cigar:
            raise ValueError(
                f"{interval.source_line}: empty elected interval slice "
                f"{spec['key'][0]}-{spec['key'][1]}"
            )
        output[spec["key"]] = (
            cigar,
            _segments_geometry(
                sliced_segments, interval.source_query_name,
            ),
        )
        full_sliced_segments = []
        for segment_index in sorted(spec["segments"]):
            result = spec["segments"][segment_index]
            if not result["full_ops"]:
                continue
            start = result["full_ref_min"]
            end = result["full_ref_max"]
            if start is None or end is None:
                start = end = result["anchor"]
            full_sliced_segments.append(dataclasses.replace(
                result["template"], start=start, end=end, cigar="",
                ops=result["full_ops"],
            ))
        full_cigar = core.format_reference_cigar_segments(
            full_sliced_segments,
        )
        if not full_cigar:
            raise ValueError(
                f"{interval.source_line}: empty full interval slice "
                f"{spec['key'][0]}-{spec['key'][1]}"
            )
        full_output[spec["key"]] = (
            full_cigar,
            _segments_geometry(
                full_sliced_segments, interval.source_query_name,
            ),
        )
    return output, deletion_intervals, full_output


def _reference_boundary_query_cuts(
    interval: AlignmentInterval, query_start: int, query_end: int,
    boundaries: Iterable[int],
    segments: Optional[Sequence[core.ReferenceCigarSegment]] = None,
) -> List[int]:
    """Map elected reference boundaries through the original source CIGAR."""
    wanted = sorted(set(boundaries))
    if not wanted or interval.trimmed:
        return []
    source_query = core.parse_coord(interval.source_query_interval)
    if source_query is None:
        return []
    cuts = set()
    qcursor = 0
    if segments is None:
        segments = core.parse_reference_cigar_segments(
            interval.graph_cigar, interval.source_query_name,
        )
    for segment in segments:
        rcursor = segment.start if segment.direction == ">" else segment.end
        for op in segment.ops:
            if op.op == "H" or op.n <= 0:
                continue
            qsize = op.n if op.op in QUERY_OPS else 0
            rsize = op.n if segment.is_main and op.op in REFERENCE_OPS else 0
            if rsize:
                rnext = (
                    rcursor + rsize if segment.direction == ">"
                    else rcursor - rsize
                )
                low, high = sorted((rcursor, rnext))
                first = bisect.bisect_left(wanted, low)
                stop = bisect.bisect_right(wanted, high)
                for boundary in wanted[first:stop]:
                    if qsize:
                        offset = (
                            boundary - rcursor
                            if segment.direction == ">"
                            else rcursor - boundary
                        )
                        local_query = qcursor + max(0, min(qsize, offset))
                    else:
                        local_query = qcursor
                    physical = (
                        source_query.start + local_query
                        if source_query.strand == "+"
                        else source_query.end - local_query
                    )
                    if query_start < physical < query_end:
                        cuts.add(physical)
                rcursor = rnext
            qcursor += qsize
    return sorted(cuts)


def _pa_metadata(
    interval: AlignmentInterval, start: int, end: int,
) -> dict:
    lift = _containing_lift(interval.pa_lifts, start, end)
    return {
        "pa_name": lift.allele if lift else "",
        "region_tier": (
            "trimmed" if interval.trimmed
            else _region_for_interval(start, end, lift, 0)
        ),
        # Stage/tier belongs to the mapping row, not merely to the PA label.
        # graphcigartoref carries GenomeLift's final annotation on each row;
        # old inputs still populate these interval fields from the coord map.
        "lift_stage": interval.lift_stage,
        "secondary_lift": interval.secondary_lift,
        "lift_kind": interval.lift_kind,
        "lift_role": interval.lift_role,
        "resurrected_lift": lift.resurrected if lift else False,
    }


def _internal_unmapped_gap_rows(

    query_spans: Sequence[Span], candidates: Sequence[object], owner_key: str,

) -> List[dict]:

    """Report uncovered sequence between elected query intervals as UNMAPPED."""

    by_contig: Dict[str, List[Span]] = defaultdict(list)

    for span in query_spans:

        by_contig[span.contig].append(span)

    candidate_by_id = {

        getattr(candidate, "interval_id", getattr(candidate, "atom_id", -1)):

        candidate for candidate in candidates

    }

    rows: List[dict] = []

    for contig, spans in by_contig.items():

        spans.sort(key=lambda span: (span.start, span.end))

        for left, right in zip(spans, spans[1:]):

            if right.start <= left.end:

                continue

            left_owner = candidate_by_id[left.owner_id]

            right_owner = candidate_by_id[right.owner_id]

            strand = (

                left_owner.query_strand

                if left_owner.query_strand == right_owner.query_strand

                else "+"

            )

            rows.append({

                "kind": "UNMAPPED", "query_contig": contig,

                "query_start": left.end, "query_end": right.start,

                "query_strand": strand,

                "reference_path": ".", "reference_start": ".",

                "reference_end": ".", "reference_strand": ">",

                "alternative_path": ".", "alternative_start": ".",

                "alternative_end": ".", "alternative_strand": ".",

                "anchor": ".", "dependency": ".",

                "operation": "UNMAPPED", owner_key: left_owner,

                "source_line": ".", "pa_name": "",

                "reference_pa": ".", "region_tier": "unmapped",

                "trimmed": 0, "encoded": 0, "lift_stage": 0,

                "secondary_lift": 0, "resurrected_lift": 0,

                "alignment_score": 0.0,

            })

    return rows


def _materialize_interval_assignment_rows(
    intervals: Sequence[AlignmentInterval],
    query_spans: Sequence[Span],
    reference_spans: Sequence[Span],
    *,
    add_unmapped: bool = False,
) -> List[dict]:
    """Slice pre-elected interval spans from their original CIGARs."""
    reference_index = _ReferenceSpanIndex(reference_spans)
    parsed_intervals = {}
    initial_ranges = defaultdict(list)
    for qspan in query_spans:
        interval = intervals[qspan.owner_id]
        if interval.trimmed:
            continue
        if interval.interval_id not in parsed_intervals:
            parsed_intervals[interval.interval_id] = (
                core.parse_reference_cigar_segments(
                    interval.graph_cigar, interval.source_query_name,
                )
            )
        initial_ranges[interval.interval_id].append(
            (qspan.start, qspan.end)
        )

    # First obtain each elected query span's projected reference geometry.
    # Each source interval is parsed once and walked once for this pass.
    initial_slices = {}
    for interval_id, ranges in initial_ranges.items():
        interval = intervals[interval_id]
        sliced, _deletions, _full_sliced = _slice_interval_ranges(
            interval, ranges, parsed_intervals[interval_id],
        )
        initial_slices[interval_id] = sliced

    requests = []
    for request_id, qspan in enumerate(query_spans):
        interval = intervals[qspan.owner_id]
        if interval.trimmed:
            continue
        _cigar, geometry = initial_slices[interval.interval_id][
            (qspan.start, qspan.end)
        ]
        path, start, end = geometry[0], geometry[1], geometry[2]
        if end > start:
            requests.append((request_id, path, start, end))
    overlaps = _overlaps_for_projected_spans(reference_spans, requests)

    # Map all reference-election boundaries through each source CIGAR in one
    # pass. Mapping separately for every query span re-read the same long
    # operation list once per competing interval.
    reference_boundaries = defaultdict(set)
    for request_id, qspan in enumerate(query_spans):
        interval = intervals[qspan.owner_id]
        if interval.trimmed:
            continue
        _whole_cigar, geometry = initial_slices[interval.interval_id][
            (qspan.start, qspan.end)
        ]
        whole_start, whole_end = geometry[1], geometry[2]
        for span in overlaps.get(request_id, ()):
            for coordinate in (span.start, span.end):
                if whole_start < coordinate < whole_end:
                    reference_boundaries[interval.interval_id].add(coordinate)
    mapped_reference_cuts = {}
    for interval_id, boundaries in reference_boundaries.items():
        interval = intervals[interval_id]
        mapped_reference_cuts[interval_id] = (
            _reference_boundary_query_cuts(
                interval, interval.query_start, interval.query_end,
                boundaries, parsed_intervals[interval_id],
            )
        )

    request_ranges = {}
    final_ranges = defaultdict(list)
    for request_id, qspan in enumerate(query_spans):
        interval = intervals[qspan.owner_id]
        label_boundaries = {
            boundary
            for lift in interval.pa_lifts
            for boundary in (lift.locus.start, lift.locus.end)
            if qspan.start < boundary < qspan.end
        }
        cuts = {qspan.start, qspan.end, *label_boundaries}
        if not interval.trimmed:
            cuts.update(
                cut for cut in mapped_reference_cuts.get(
                    interval.interval_id, (),
                )
                if qspan.start < cut < qspan.end
            )
        ordered = sorted(cuts)
        pieces = [
            (start, end) for start, end in zip(ordered, ordered[1:])
            if end > start
        ]
        request_ranges[request_id] = pieces
        if not interval.trimmed:
            final_ranges[interval.interval_id].extend(pieces)

    # Materialize every final CIGAR slice in one more linear walk per source
    # interval. The resulting text is retained on the report row so writing
    # the TSV never re-slices the source CIGAR.
    final_slices = {}
    final_full_slices = {}
    deletion_slices = defaultdict(list)
    for interval_id, ranges in final_ranges.items():
        interval = intervals[interval_id]
        sliced, deletions, full_sliced = _slice_interval_ranges(
            interval, ranges, parsed_intervals[interval_id], reference_index,
        )
        final_slices[interval_id] = sliced
        final_full_slices[interval_id] = full_sliced
        deletion_slices[interval_id].extend(deletions)

    rows: List[dict] = []
    for request_id, qspan in enumerate(query_spans):
        interval = intervals[qspan.owner_id]
        if interval.trimmed:
            for query_start, query_end in request_ranges[request_id]:
                rows.append({
                    "kind": "UNMAPPED", "query_contig": qspan.contig,
                    "query_start": query_start, "query_end": query_end,
                    "query_strand": interval.query_strand,
                    "reference_path": ".", "reference_start": ".",
                    "reference_end": ".", "reference_strand": ">",
                    "alternative_path": ".", "alternative_start": ".",
                    "alternative_end": ".", "alternative_strand": ".",
                    "anchor": ".", "dependency": ".",
                    "operation": "UNMAPPED", "encoded": False,
                    "interval": interval,
                    "pa_name": "", "reference_pa": ".",
                    "region_tier": "unmapped", "lift_stage": 0,
                    "secondary_lift": False, "resurrected_lift": False,
                })
            continue
        for query_start, query_end in request_ranges[request_id]:
            cigar, geometry = final_slices[interval.interval_id][
                (query_start, query_end)
            ]
            full_cigar, full_geometry = final_full_slices[
                interval.interval_id
            ][(query_start, query_end)]
            (
                path, reference_start, reference_end, reference_strand,
                _path_length, anchor_path, anchor_position, encoded,
            ) = geometry
            primary = reference_index.owns(
                path, reference_start, reference_end, interval.interval_id,
            )
            common = {
                "query_contig": qspan.contig,
                "query_start": query_start, "query_end": query_end,
                "query_strand": interval.query_strand,
                "reference_strand": reference_strand,
                "dependency": ".", "operation": "ALIGNMENT",
                "encoded": encoded, "interval": interval,
                "alignment_cigar": cigar,
                **_pa_metadata(interval, query_start, query_end),
            }
            if primary:
                rows.append({
                    **common, "kind": "PRIMARY",
                    "reference_path": path,
                    "reference_start": reference_start,
                    "reference_end": reference_end,
                    "alternative_path": ".", "alternative_start": ".",
                    "alternative_end": ".", "alternative_strand": ".",
                    "anchor": ".",
                })
            else:
                (
                    alternative_path, alternative_start, alternative_end,
                    alternative_strand,
                    _alternative_path_length, _alternative_anchor_path,
                    _alternative_anchor_position, _alternative_encoded,
                ) = full_geometry
                has_alternative = alternative_end > alternative_start
                rows.append({
                    **common, "kind": "INSERTION",
                    "alignment_cigar": full_cigar,
                    "reference_path": (
                        "." if has_alternative else anchor_path or "."
                    ),
                    "reference_start": (
                        "." if has_alternative else anchor_position
                    ),
                    "reference_end": (
                        "." if has_alternative else anchor_position
                    ),
                    "alternative_path": (
                        alternative_path if has_alternative else "."
                    ),
                    "alternative_start": (
                        alternative_start if has_alternative else "."
                    ),
                    "alternative_end": (
                        alternative_end if has_alternative else "."
                    ),
                    "alternative_strand": (
                        alternative_strand if has_alternative else "."
                    ),
                    "anchor": (
                        "." if has_alternative else anchor_position
                    ),
                })

    # A reference winner can divide a zero-query deletion. Preserve every
    # portion still owned by this interval as an explicit reference interval;
    # there is no nonzero query interval on which such a split can be encoded.
    for interval_id, deletions in deletion_slices.items():
        interval = intervals[interval_id]
        for (
            query_point, path, reference_start, reference_end,
            reference_strand, path_length,
        ) in deletions:
            cigar = core._format_interval_gcigar(
                reference_strand, path, max(path_length, reference_end),
                reference_start, reference_end,
                [core.CigarOp(reference_end - reference_start, "D", "")],
            )
            rows.append({
                "kind": "DELETION",
                "query_contig": interval.query_contig,
                "query_start": query_point,
                "query_end": query_point,
                "query_strand": interval.query_strand,
                "reference_path": path,
                "reference_start": reference_start,
                "reference_end": reference_end,
                "reference_strand": reference_strand,
                "alternative_path": ".", "alternative_start": ".",
                "alternative_end": ".", "alternative_strand": ".",
                "anchor": ".", "dependency": ".", "operation": "D",
                "encoded": False, "interval": interval,
                "alignment_cigar": cigar,
                **_pa_metadata(interval, query_point, query_point),
            })

    if add_unmapped:
        rows.extend(_internal_unmapped_gap_rows(
            query_spans, intervals, "interval",
        ))
    return rows


def _merge_coordinate_ranges(
    ranges: Iterable[Tuple[int, int]],
) -> List[Tuple[int, int]]:
    ordered = sorted((start, end) for start, end in ranges if end > start)
    merged: List[Tuple[int, int]] = []
    for start, end in ordered:
        if merged and start <= merged[-1][1]:
            merged[-1] = (merged[-1][0], max(merged[-1][1], end))
        else:
            merged.append((start, end))
    return merged


def _subtract_ranges(
    start: int, end: int, masks: Sequence[Tuple[int, int]],
) -> List[Tuple[int, int]]:
    """Subtract sorted disjoint masks from one half-open interval."""
    if end <= start or not masks:
        return [(start, end)] if end > start else []
    index = bisect.bisect_right(masks, (start, math.inf)) - 1
    if index < 0:
        index = 0
    elif masks[index][1] <= start:
        index += 1
    cursor = start
    output = []
    while index < len(masks) and masks[index][0] < end:
        mask_start, mask_end = masks[index]
        if mask_start > cursor:
            output.append((cursor, min(mask_start, end)))
        cursor = max(cursor, mask_end)
        if cursor >= end:
            break
        index += 1
    if cursor < end:
        output.append((cursor, end))
    return output


def _large_insertion_ranges(
    interval: AlignmentInterval,
) -> List[Tuple[int, int]]:
    """Physical query spans of >=500-bp query-only CIGAR operations."""
    if interval.trimmed or interval.query_end <= interval.query_start:
        return []
    if interval.encoded:
        return (
            [(interval.query_start, interval.query_end)]
            if interval.query_end - interval.query_start
            >= LARGE_INSERTION_MIN_QUERY else []
        )
    source = core.parse_coord(interval.source_query_interval)
    if source is None:
        return []
    ranges = []
    qcursor = 0
    for segment in core.parse_reference_cigar_segments(
        interval.graph_cigar, interval.source_query_name,
    ):
        for operation in segment.ops:
            if operation.op == "H" or operation.n <= 0:
                continue
            qsize = operation.n if operation.op in QUERY_OPS else 0
            if operation.op == "I" and operation.n >= LARGE_INSERTION_MIN_QUERY:
                start, end = _physical_query(
                    source, qcursor, qcursor + operation.n,
                )
                start = max(start, interval.query_start)
                end = min(end, interval.query_end)
                if end > start:
                    ranges.append((start, end))
            qcursor += qsize
    return _merge_coordinate_ranges(ranges)


def _split_spans(
    spans: Sequence[Span], boundaries_by_contig: Mapping[str, Sequence[int]],
) -> List[Span]:
    output = []
    for span in spans:
        boundaries = boundaries_by_contig.get(span.contig, ())
        left = bisect.bisect_right(boundaries, span.start)
        right = bisect.bisect_left(boundaries, span.end)
        cuts = (span.start, *boundaries[left:right], span.end)
        output.extend(
            Span(span.contig, start, end, span.owner_id)
            for start, end in zip(cuts, cuts[1:]) if end > start
        )
    return output


def _covered_by(
    start: int, end: int, masks: Sequence[Tuple[int, int]],
) -> bool:
    if end <= start or not masks:
        return False
    index = bisect.bisect_right(masks, (start, math.inf)) - 1
    return index >= 0 and masks[index][0] <= start and masks[index][1] >= end


def _finish_assignment_rows(
    rows: List[dict], query_spans: Sequence[Span],
    intervals: Sequence[AlignmentInterval],
) -> List[dict]:
    rows.extend(_internal_unmapped_gap_rows(
        query_spans, intervals, "interval",
    ))
    resolve_insertion_anchors(rows)
    for index, row in enumerate(rows, 1):
        row["assignment_id"] = f"A{index}"
    for row in rows:
        dependency = row.get("dependency", ".")
        if isinstance(dependency, str) and dependency.startswith("@"):
            row["dependency"] = f"A{int(dependency[1:]) + 1}"
    return rows


def _two_pass_interval_assignment_rows(
    intervals: Sequence[AlignmentInterval],
) -> List[dict]:
    """Elect primary mappings, then alternatives only on unresolved query."""
    primary_ids = {
        item.interval_id for item in intervals
        if (
            (
                item.lift_stage == 1
                or (
                    not item.secondary_lift
                    and item.lift_stage in {2, 3}
                )
                # New gap-fill alignments have no PA lift of their own. They
                # are primary evidence between already elected anchors.
                or (
                    not item.secondary_lift
                    and not item.final_lift_annotation
                    and item.lift_stage == 0
                )
            )
        )
    }
    secondary_ids = {
        item.interval_id for item in intervals
        if (
            item.secondary_lift
            and item.lift_stage in {2, 3, 4}
        )
    }

    primary_query = elect([
        Span(item.query_contig, item.query_start, item.query_end,
             item.interval_id)
        for item in intervals
        if item.interval_id in primary_ids and item.query_end > item.query_start
    ], intervals, priority_key=lambda item: _pass_priority(item, False))
    primary_reference = elect([
        Span(item.reference_path, item.reference_start, item.reference_end,
             item.interval_id)
        for item in intervals
        if (
            item.interval_id in primary_ids
            and not item.trimmed
            and item.reference_end > item.reference_start
        )
    ], intervals, priority_key=lambda item: _pass_priority(item, False))

    # Query sequence aligned by the primary pass is unavailable to duplicate
    # mappings.  A >=500-bp query-only operation remains available so another
    # alignment can explain it as a duplication source.
    resolved_by_contig = defaultdict(list)
    insertion_boundaries = defaultdict(set)
    # A large query-only operation is deliberately reopened for the
    # alternative pass.  Keep the primary alignment's reference interval with
    # that query window so the same physical alignment cannot explain its own
    # insertion as a duplication.  This exclusion is local to the reopened
    # insertion; it does not reserve that reference against other secondary
    # query intervals.
    insertion_reference_exclusions = defaultdict(list)
    large_insertions_by_interval = {}
    for span in primary_query:
        interval = intervals[span.owner_id]
        if interval.trimmed:
            continue
        unresolved = []
        if interval.interval_id not in large_insertions_by_interval:
            large_insertions_by_interval[interval.interval_id] = (
                _large_insertion_ranges(interval)
            )
        for start, end in large_insertions_by_interval[interval.interval_id]:
            start, end = max(start, span.start), min(end, span.end)
            if end > start:
                unresolved.append((start, end))
                insertion_boundaries[span.contig].update((start, end))
                if (
                    interval.reference_path
                    and interval.reference_end > interval.reference_start
                ):
                    insertion_reference_exclusions[span.contig].append((
                        start, end, interval.reference_path,
                        interval.reference_start, interval.reference_end,
                    ))
        for start, end in _subtract_ranges(
            span.start, span.end, _merge_coordinate_ranges(unresolved),
        ):
            resolved_by_contig[span.contig].append((start, end))
    resolved_by_contig = {
        contig: _merge_coordinate_ranges(values)
        for contig, values in resolved_by_contig.items()
    }
    insertion_reference_exclusions = {
        contig: sorted(values)
        for contig, values in insertion_reference_exclusions.items()
    }

    secondary_candidates = []
    for item in intervals:
        if item.interval_id not in secondary_ids or item.trimmed:
            continue
        for start, end in _subtract_ranges(
            item.query_start, item.query_end,
            resolved_by_contig.get(item.query_contig, ()),
        ):
            forbidden = []
            windows = insertion_reference_exclusions.get(
                item.query_contig, (),
            )
            index = bisect.bisect_right(
                windows, (start, math.inf),
            ) - 1
            if index < 0:
                index = 0
            elif windows[index][1] <= start:
                index += 1
            while index < len(windows) and windows[index][0] < end:
                (
                    query_start, query_end, reference_path,
                    reference_start, reference_end,
                ) = windows[index]
                if (
                    item.reference_path == reference_path
                    and min(item.reference_end, reference_end)
                    > max(item.reference_start, reference_start)
                ):
                    overlap_start = max(start, query_start)
                    overlap_end = min(end, query_end)
                    if overlap_end > overlap_start:
                        forbidden.append((overlap_start, overlap_end))
                index += 1
            for candidate_start, candidate_end in _subtract_ranges(
                start, end, _merge_coordinate_ranges(forbidden),
            ):
                secondary_candidates.append(Span(
                    item.query_contig, candidate_start, candidate_end,
                    item.interval_id,
                ))
    secondary_query = elect(
        secondary_candidates, intervals,
        priority_key=lambda item: _pass_priority(item, True),
    )

    secondary_coverage = defaultdict(list)
    split_boundaries = defaultdict(set)
    for span in secondary_query:
        secondary_coverage[span.contig].append((span.start, span.end))
        split_boundaries[span.contig].update((span.start, span.end))
    for contig, boundaries in insertion_boundaries.items():
        split_boundaries[contig].update(boundaries)
    secondary_coverage = {
        contig: _merge_coordinate_ranges(values)
        for contig, values in secondary_coverage.items()
    }
    split_boundaries = {
        contig: sorted(values) for contig, values in split_boundaries.items()
    }
    primary_query = _split_spans(primary_query, split_boundaries)

    primary_rows = _materialize_interval_assignment_rows(
        intervals, primary_query, primary_reference,
    )
    secondary_rows = _materialize_interval_assignment_rows(
        intervals, secondary_query, (),
    )

    # Secondary mappings replace only unresolved primary pieces.  Primary
    # rows covering aligned sequence cannot overlap secondary_query because
    # those intervals were removed by resolved_by_contig above.
    retained_primary = []
    for row in primary_rows:
        start, end = row.get("query_start"), row.get("query_end")
        if (
            isinstance(start, int) and isinstance(end, int) and end > start
            and _covered_by(
                start, end, secondary_coverage.get(row["query_contig"], ()),
            )
        ):
            continue
        retained_primary.append(row)

    final_spans = []
    for span in primary_query:
        if not _covered_by(
            span.start, span.end,
            secondary_coverage.get(span.contig, ()),
        ):
            final_spans.append(span)
    final_spans.extend(secondary_query)
    final_spans.sort(key=lambda span: (
        span.contig, span.start, span.end, span.owner_id,
    ))
    return _finish_assignment_rows(
        retained_primary + secondary_rows, final_spans, intervals,
    )


def _attach_alternative_targets(
    main_rows: Sequence[dict], alternative_rows: Sequence[dict],
) -> None:
    """Place sequence-source insertions on their overlapping location map.

    Both row lists are already non-overlapping within their own layer.  A
    sorted two-pointer walk therefore finds the location pieces covering each
    alternative without a pairwise scan.
    """
    main_by_contig: Dict[str, List[dict]] = defaultdict(list)
    alternative_by_contig: Dict[str, List[dict]] = defaultdict(list)
    for row in main_rows:
        if (
            (
                row.get("kind") == "PRIMARY"
                or (
                    row.get("kind") == "INSERTION"
                    and row.get("alternative_path") in {None, "", "."}
                )
            )
            and isinstance(row.get("query_start"), int)
            and isinstance(row.get("query_end"), int)
            and row.get("reference_path") not in {None, "", "."}
            and isinstance(row.get("reference_start"), int)
            and isinstance(row.get("reference_end"), int)
        ):
            main_by_contig[row["query_contig"]].append(row)
    for row in alternative_rows:
        if (
            row.get("kind") == "INSERTION"
            and row.get("alternative_path") not in {None, "", "."}
            and isinstance(row.get("query_start"), int)
            and isinstance(row.get("query_end"), int)
        ):
            alternative_by_contig[row["query_contig"]].append(row)

    group_number = 0
    for contig, alternatives in alternative_by_contig.items():
        mains = sorted(main_by_contig.get(contig, ()), key=lambda row: (
            row["query_start"], row["query_end"], row["reference_path"],
        ))
        alternatives.sort(key=lambda row: (
            row["query_start"], row["query_end"], row["alternative_path"],
        ))
        first = 0
        for alternative in alternatives:
            while (
                first < len(mains)
                and mains[first]["query_end"] <= alternative["query_start"]
            ):
                first += 1
            pieces = []
            index = first
            while (
                index < len(mains)
                and mains[index]["query_start"] < alternative["query_end"]
            ):
                pieces.append(mains[index])
                index += 1
            if not pieces:
                continue
            pieces.sort(key=lambda row: (row["query_start"], row["query_end"]))
            if (
                pieces[0]["query_start"] != alternative["query_start"]
                or pieces[-1]["query_end"] != alternative["query_end"]
                or any(
                    left["query_end"] != right["query_start"]
                    for left, right in zip(pieces, pieces[1:])
                )
            ):
                continue
            paths = {row.get("reference_path") for row in pieces}
            strands = {row.get("reference_strand") for row in pieces}
            if len(paths) != 1 or len(strands) != 1 or "." in paths:
                continue
            reference_start = min(row["reference_start"] for row in pieces)
            reference_end = max(row["reference_end"] for row in pieces)
            group_number += 1
            group_id = (
                f"S{group_number}:{contig}:"
                f"{alternative['query_start']}-{alternative['query_end']}"
            )
            alternative["reference_path"] = next(iter(paths))
            alternative["reference_start"] = reference_start
            alternative["reference_end"] = reference_end
            alternative["reference_strand"] = next(iter(strands))
            alternative["anchor"] = reference_start
            alternative["dependency"] = "."
            alternative["dup_group"] = group_id
            alternative["insertion_group"] = group_id
            alternative["insertion_anchor_path"] = next(iter(paths))
            alternative["insertion_anchor_start"] = reference_start
            alternative["insertion_anchor_end"] = reference_end


def _separate_mapping_assignment_rows(
    intervals: Sequence[AlignmentInterval],
) -> List[dict]:
    """Elect location/main and sequence/alternative mappings independently."""
    main_ids = {
        item.interval_id for item in intervals
        if item.lift_role != "alternative"
    }
    alternative_ids = {
        item.interval_id for item in intervals
        if item.lift_role == "alternative"
    }
    main_query = elect([
        Span(item.query_contig, item.query_start, item.query_end,
             item.interval_id)
        for item in intervals
        if item.interval_id in main_ids and item.query_end > item.query_start
    ], intervals, priority_key=lambda item: _pass_priority(item, False))
    alternative_query = elect([
        Span(item.query_contig, item.query_start, item.query_end,
             item.interval_id)
        for item in intervals
        if (
            item.interval_id in alternative_ids
            and item.query_end > item.query_start
        )
    ], intervals, priority_key=lambda item: _pass_priority(item, True))

    # Split the main election at alternative boundaries before CIGAR slicing.
    # This makes each alternative's location target an exact union of main
    # pieces and avoids re-parsing long CIGARs during target attachment.
    boundaries = defaultdict(set)
    for span in alternative_query:
        boundaries[span.contig].update((span.start, span.end))
    main_query = _split_spans(
        main_query,
        {contig: sorted(values) for contig, values in boundaries.items()},
    )
    main_reference = elect([
        Span(item.reference_path, item.reference_start, item.reference_end,
             item.interval_id)
        for item in intervals
        if (
            item.interval_id in main_ids
            and not item.trimmed
            # GenomeLift elected the owner: a location duplicate keeps its
            # query span but never owns the reference.
            and not item.secondary_lift
            and item.reference_end > item.reference_start
        )
    ], intervals, priority_key=lambda item: _pass_priority(item, False))
    main_rows = _materialize_interval_assignment_rows(
        intervals, main_query, main_reference,
    )
    alternative_rows = _materialize_interval_assignment_rows(
        intervals, alternative_query, (),
    )
    _attach_alternative_targets(main_rows, alternative_rows)
    return _finish_assignment_rows(
        main_rows + alternative_rows, main_query, intervals,
    )


def interval_assignment_rows(
    intervals: Sequence[AlignmentInterval],
) -> List[dict]:
    """Run mapping-source election, with legacy staged compatibility."""
    if any(item.lift_kind for item in intervals):
        return _separate_mapping_assignment_rows(intervals)
    if any(item.final_lift_annotation for item in intervals):
        return _two_pass_interval_assignment_rows(intervals)

    query_spans = elect([
        Span(item.query_contig, item.query_start, item.query_end,
             item.interval_id)
        for item in intervals if item.query_end > item.query_start
    ], intervals)
    reference_spans = elect([
        Span(item.reference_path, item.reference_start, item.reference_end,
             item.interval_id)
        for item in intervals
        if (
            not item.trimmed
            and not _is_owned_dup_pa(item)
            and item.reference_end > item.reference_start
        )
    ], intervals)
    rows = _materialize_interval_assignment_rows(
        intervals, query_spans, reference_spans,
    )
    return _finish_assignment_rows(rows, query_spans, intervals)


def _legacy_assignment_rows(atoms: Sequence[_LegacyAtom]) -> List[dict]:
    query_spans = elect([
        Span(a.query_contig, a.query_start, a.query_end, a.atom_id)
        for a in atoms if a.query_end > a.query_start
    ], atoms)
    reference_spans = elect([
        Span(a.reference_path, a.reference_start, a.reference_end, a.atom_id)
        for a in atoms
        if (
            a.reference_end > a.reference_start
            and not a.encoded
            and not _is_owned_dup_pa(a)
        )
    ], atoms)
    projected_requests = []
    for request_id, qspan in enumerate(query_spans):
        atom = atoms[qspan.owner_id]
        if not atom.encoded and atom.reference_end > atom.reference_start:
            r0, r1 = _project_query(atom, qspan.start, qspan.end)
            projected_requests.append(
                (request_id, atom.reference_path, r0, r1),
            )
    projected_overlaps = _overlaps_for_projected_spans(
        reference_spans, projected_requests,
    )
    rows: List[dict] = []

    for request_id, qspan in enumerate(query_spans):
        atom = atoms[qspan.owner_id]
        if atom.trimmed:
            rows.append({
                "kind": "UNMAPPED", "query_contig": qspan.contig,
                "query_start": qspan.start, "query_end": qspan.end,
                "query_strand": atom.query_strand,
                "reference_path": ".", "reference_start": ".",
                "reference_end": ".", "reference_strand": ">",
                "alternative_path": ".", "alternative_start": ".",
                "alternative_end": ".", "alternative_strand": ".",
                "anchor": ".", "dependency": ".",
                "operation": "UNMAPPED", "atom": atom,
                "pa_name": "", "reference_pa": ".",
                "region_tier": "unmapped", "lift_stage": 0,
                "secondary_lift": False, "resurrected_lift": False,
            })
            continue
        paired = projected_overlaps.get(request_id, ())
        boundaries = {qspan.start, qspan.end}
        query_owners = []
        if paired:
            mapped_start, mapped_end = _project_query(atom, qspan.start, qspan.end)
            for ref_span in paired:
                overlap_start = max(mapped_start, ref_span.start)
                overlap_end = min(mapped_end, ref_span.end)
                if overlap_end > overlap_start:
                    query_owner_start, query_owner_end = _project_reference(
                        atom, overlap_start, overlap_end,
                    )
                    boundaries.update((query_owner_start, query_owner_end))
                    query_owners.append(
                        (query_owner_start, query_owner_end, ref_span),
                    )
        if len(query_owners) > 1:
            query_owners.sort(key=lambda item: (item[0], item[1]))
        if len(boundaries) == 2:
            ordered = [qspan.start, qspan.end]
        else:
            ordered = sorted(
                x for x in boundaries if qspan.start <= x <= qspan.end
            )
        owner_index = 0
        for query_start, query_end in zip(ordered, ordered[1:]):
            if query_end <= query_start:
                continue
            r0, r1 = _project_query(atom, query_start, query_end)
            while (
                owner_index < len(query_owners)
                and query_owners[owner_index][1] <= query_start
            ):
                owner_index += 1
            owner = (
                query_owners[owner_index][2]
                if owner_index < len(query_owners)
                and query_owners[owner_index][0] <= query_start
                and query_owners[owner_index][1] >= query_end
                else None
            )
            is_primary = owner is not None and owner.owner_id == atom.atom_id
            if is_primary:
                rows.append({
                    "kind": "PRIMARY", "query_contig": qspan.contig,
                    "query_start": query_start, "query_end": query_end,
                    "query_strand": atom.query_strand,
                    "reference_path": atom.reference_path,
                    "reference_start": r0, "reference_end": r1,
                    "reference_strand": atom.reference_strand,
                    "alternative_path": ".", "alternative_start": ".",
                    "alternative_end": ".", "anchor": ".", "dependency": ".",
                    "_alignment_pieces": [(atom, query_start, query_end)],
                    "atom": atom,
                })
                continue
            anchor_path = atom.anchor_path
            anchor = atom.anchor_position
            # A losing aligned span is sequence inserted at this query
            # location, while r0/r1 describe its alternative source.  Leave
            # its target anchor unresolved so physical query neighbours choose
            # the insertion point.  Literal I and encoded insertion atoms
            # already carry their true main-path anchor.
            if (
                not atom.encoded
                and atom.reference_end > atom.reference_start
                and atom.op in {"=", "X", "M"}
            ):
                anchor_path = ""
                anchor = -1
            has_alternative_alignment = (
                not atom.encoded
                and atom.op in {"=", "X", "M"}
                and atom.reference_end > atom.reference_start
            )
            rows.append({
                "kind": "INSERTION", "query_contig": qspan.contig,
                "query_start": query_start, "query_end": query_end,
                "query_strand": atom.query_strand,
                "reference_path": anchor_path or ".",
                "reference_start": anchor if anchor >= 0 else ".",
                "reference_end": anchor if anchor >= 0 else ".",
                "reference_strand": atom.reference_strand,
                "alternative_path": (
                    atom.reference_path
                    if has_alternative_alignment and atom.reference_path else "."
                ),
                "alternative_start": r0 if has_alternative_alignment else ".",
                "alternative_end": r1 if has_alternative_alignment else ".",
                "alternative_strand": (
                    atom.reference_strand if has_alternative_alignment else "."
                ),
                "anchor": anchor if anchor >= 0 else ".", "dependency": ".",
                "_alignment_pieces": [(atom, query_start, query_end)],
                "atom": atom,
            })

    reference_by_atom: Dict[int, List[Tuple[int, Span]]] = defaultdict(list)
    primary_by_atom: Dict[int, List[dict]] = defaultdict(list)
    for reference_index, rspan in enumerate(reference_spans):
        reference_by_atom[rspan.owner_id].append((reference_index, rspan))
    for row in rows:
        if row["kind"] == "PRIMARY":
            primary_by_atom[row["atom"].atom_id].append(row)

    deletions_by_reference: Dict[int, List[dict]] = defaultdict(list)
    for atom_id, indexed_reference_spans in reference_by_atom.items():
        atom = atoms[atom_id]
        indexed_reference_spans.sort(
            key=lambda item: (item[1].start, item[1].end),
        )
        covered_rows = primary_by_atom.get(atom_id, [])
        covered_rows.sort(
            key=lambda row: (row["reference_start"], row["reference_end"]),
        )
        first_covered = 0
        for reference_index, rspan in indexed_reference_spans:
            while (
                first_covered < len(covered_rows)
                and covered_rows[first_covered]["reference_end"] <= rspan.start
            ):
                first_covered += 1
            cursor = rspan.start
            covered_index = first_covered
            while (
                covered_index < len(covered_rows)
                and covered_rows[covered_index]["reference_start"] < rspan.end
            ):
                row = covered_rows[covered_index]
                covered_start = max(rspan.start, row["reference_start"])
                covered_end = min(rspan.end, row["reference_end"])
                if covered_start > cursor:
                    deletions_by_reference[reference_index].append({
                        "kind": "DELETION", "query_contig": atom.query_contig,
                        "query_start": ".", "query_end": ".",
                        "query_strand": atom.query_strand,
                        "reference_path": rspan.contig,
                        "reference_start": cursor, "reference_end": covered_start,
                        "reference_strand": atom.reference_strand,
                        "alternative_path": ".", "alternative_start": ".",
                        "alternative_end": ".", "anchor": ".", "dependency": ".",
                        "atom": atom,
                    })
                cursor = max(cursor, covered_end)
                if cursor >= rspan.end:
                    break
                covered_index += 1
            if cursor < rspan.end:
                deletions_by_reference[reference_index].append({
                    "kind": "DELETION", "query_contig": atom.query_contig,
                    "query_start": ".", "query_end": ".",
                    "query_strand": atom.query_strand,
                    "reference_path": rspan.contig,
                    "reference_start": cursor, "reference_end": rspan.end,
                    "reference_strand": atom.reference_strand,
                    "alternative_path": ".", "alternative_start": ".",
                    "alternative_end": ".", "anchor": ".", "dependency": ".",
                    "atom": atom,
                })
    deletion_candidates = []
    for reference_index in range(len(reference_spans)):
        deletion_candidates.extend(deletions_by_reference.get(reference_index, ()))
    rows.extend(_retain_bracketed_deletions(rows, deletion_candidates))

    rows.extend(_internal_unmapped_gap_rows(

        query_spans, atoms, "atom",

    ))

    resolve_insertion_anchors(rows)
    rows = _legacy_coalesce_primary_rows(rows)
    for index, row in enumerate(rows, 1):
        row["assignment_id"] = f"A{index}"
    for row in rows:
        dependency = row.get("dependency", ".")
        if isinstance(dependency, str) and dependency.startswith("@"):
            row["dependency"] = f"A{int(dependency[1:]) + 1}"
    return rows


def _truthy_assignment_value(value: object) -> bool:
    return str(value or "0").strip().lower() in {"1", "true", "yes"}


def _is_duplication_assignment_row(row: Mapping[str, object]) -> bool:
    return bool(
        row.get("kind") == "INSERTION"
        and row.get("alternative_path") not in {None, "", "."}
        and row.get("pa_name") not in {None, "", "."}
        and _truthy_assignment_value(row.get("secondary_lift"))
    )


def _assign_duplication_groups(rows: List[dict]) -> None:
    """Group consecutive duplicated PAs and validate their outer anchors.

    A short UNMAPPED interval may be absorbed only when duplicated PAs occur
    on both sides.  The source alignments continue to describe the individual
    duplicated PAs, while ``dup_group`` tells VCF generation to concatenate
    their parent pieces (and the absorbed query gap) once.  This is one sorted
    pass per query contig; no pairwise PA comparison is performed.
    """
    by_contig: Dict[str, List[int]] = defaultdict(list)
    for index, row in enumerate(rows):
        if (
            isinstance(row.get("query_start"), int)
            and isinstance(row.get("query_end"), int)
            and row["query_end"] > row["query_start"]
        ):
            by_contig[row["query_contig"]].append(index)

    group_number = 0
    for contig, indexes in by_contig.items():
        indexes.sort(key=lambda index: (
            rows[index]["query_start"], rows[index]["query_end"], index,
        ))
        position = 0
        while position < len(indexes):
            first_index = indexes[position]
            if not _is_duplication_assignment_row(rows[first_index]):
                position += 1
                continue

            member_positions = [position]
            last_dup_position = position
            cursor = position + 1
            while cursor < len(indexes):
                previous = rows[indexes[member_positions[-1]]]
                candidate = rows[indexes[cursor]]
                if candidate["query_start"] != previous["query_end"]:
                    break
                if _is_duplication_assignment_row(candidate):
                    member_positions.append(cursor)
                    last_dup_position = cursor
                    cursor += 1
                    continue
                if candidate.get("kind") == "UNMAPPED":
                    # Gap reporting can produce more than one adjacent
                    # UNMAPPED row (for example, a trimmed interval beside an
                    # uncovered interval). Treat their combined query span as
                    # one gap while retaining each report row.
                    gap_positions = []
                    gap_start = candidate["query_start"]
                    gap_end = previous["query_end"]
                    gap_cursor = cursor
                    while gap_cursor < len(indexes):
                        gap = rows[indexes[gap_cursor]]
                        if (
                            gap.get("kind") != "UNMAPPED"
                            or gap["query_start"] != gap_end
                        ):
                            break
                        gap_positions.append(gap_cursor)
                        gap_end = gap["query_end"]
                        gap_cursor += 1
                    if (
                        gap_positions
                        and gap_end - gap_start < DUP_GROUP_MAX_UNMAPPED
                        and gap_cursor < len(indexes)
                    ):
                        following = rows[indexes[gap_cursor]]
                        if (
                            following["query_start"] == gap_end
                            and _is_duplication_assignment_row(following)
                        ):
                            member_positions.extend(gap_positions)
                            member_positions.append(gap_cursor)
                            last_dup_position = gap_cursor
                            cursor = gap_cursor + 1
                            continue
                break

            left_position = position - 1
            right_position = last_dup_position + 1
            valid = left_position >= 0 and right_position < len(indexes)
            if valid:
                left = rows[indexes[left_position]]
                right = rows[indexes[right_position]]
                valid = (
                    left.get("kind") == "PRIMARY"
                    and right.get("kind") == "PRIMARY"
                    and left.get("reference_path") not in {None, "", "."}
                    and left.get("reference_path") == right.get("reference_path")
                    and isinstance(left.get("reference_start"), int)
                    and isinstance(right.get("reference_start"), int)
                )
            if valid:
                left_anchor = _reference_edge_at_query_side(left, "right")
                right_anchor = _reference_edge_at_query_side(right, "left")
                valid = (
                    abs(right_anchor - left_anchor)
                    < DUP_GROUP_MAX_REFERENCE_SPAN
                )
            if valid:
                group_number += 1
                group_id = (
                    f"D{group_number}:{contig}:"
                    f"{rows[first_index]['query_start']}-"
                    f"{rows[indexes[last_dup_position]]['query_end']}"
                )
                reference_start = min(left_anchor, right_anchor)
                reference_end = max(left_anchor, right_anchor)
                leader = rows[first_index]
                leader_source_line = leader.get("source_line")
                if leader_source_line in {None, "", "."}:
                    leader_source_line = getattr(
                        _row_owner(leader), "source_line", ".",
                    )
                for member_position in member_positions:
                    member = rows[indexes[member_position]]
                    if member.get("kind") == "UNMAPPED":
                        member["kind"] = "INSERTION"
                        member["operation"] = "GROUPED_UNMAPPED"
                        member["source_line"] = leader_source_line
                    member["dup_group"] = group_id
                    member["reference_path"] = left["reference_path"]
                    member["reference_start"] = reference_start
                    member["reference_end"] = reference_end
                    # Parent sequence is rendered in the common main-path
                    # orientation; alternative source CIGARs retain their own
                    # strand independently.
                    member["query_strand"] = left.get(
                        "query_strand", member.get("query_strand", "+"),
                    )
                    member["reference_strand"] = ">"
                    member["anchor"] = reference_start
                    member["dependency"] = "."

            # Every duplicated PA considered above belongs to this maximal
            # run, valid or not.  Continue after its last duplicated member;
            # an invalid run must not be reconsidered from its second PA.
            position = max(position + 1, last_dup_position + 1)


def _assign_insertion_anchor_intervals(rows: List[dict]) -> None:
    """Record both physical reference flanks for bracketed insertions.

    ``reference_interval`` remains the point used by the materialized source
    CIGAR.  ``insertion_anchor_interval`` is the event location: the interval
    between the nearest PRIMARY rows on the two physical query sides.  A pure
    insertion has equal anchors; differing anchors describe an interval
    replacement.  The rows are sorted once per query contig.
    """
    by_contig: Dict[str, List[int]] = defaultdict(list)
    for index, row in enumerate(rows):
        if (
            isinstance(row.get("query_start"), int)
            and isinstance(row.get("query_end"), int)
            and row["query_end"] > row["query_start"]
        ):
            by_contig[row["query_contig"]].append(index)

    group_number = 0
    for indexes in by_contig.values():
        indexes.sort(key=lambda index: (
            rows[index]["query_start"], rows[index]["query_end"], index,
        ))
        position = 0
        while position < len(indexes):
            first = rows[indexes[position]]
            if first.get("kind") != "INSERTION":
                position += 1
                continue
            # Alternative-source rows are callable on the main path only when
            # _assign_duplication_groups has validated their duplicated PA.
            if (
                first.get("alternative_path") not in {None, "", "."}
                and not first.get("dup_group")
            ):
                position += 1
                continue

            member_positions = [position]
            group_kind = "dup" if first.get("dup_group") else "ordinary"
            cursor = position + 1
            while cursor < len(indexes):
                previous = rows[indexes[member_positions[-1]]]
                candidate = rows[indexes[cursor]]
                candidate_kind = (
                    "dup" if candidate.get("dup_group") else "ordinary"
                )
                if (
                    candidate.get("kind") != "INSERTION"
                    or candidate["query_start"] != previous["query_end"]
                    # Opposite-strand pieces are separate insertion events
                    # (different insertion points); never join them.
                    or candidate.get("query_strand") != previous.get("query_strand")
                    or candidate_kind != group_kind
                    or (
                        candidate.get("alternative_path")
                        not in {None, "", "."}
                        and not candidate.get("dup_group")
                    )
                ):
                    break
                member_positions.append(cursor)
                cursor += 1

            left_position = position - 1
            right_position = cursor
            valid = left_position >= 0 and right_position < len(indexes)
            if valid:
                left = rows[indexes[left_position]]
                right = rows[indexes[right_position]]
                valid = (
                    left.get("kind") == "PRIMARY"
                    and right.get("kind") == "PRIMARY"
                    and left.get("reference_path") not in {None, "", "."}
                    and left.get("reference_path") == right.get("reference_path")
                    and left.get("reference_strand")
                    == right.get("reference_strand")
                    and _same_mapping_orientation(_row_owner(left))
                    == _same_mapping_orientation(_row_owner(right))
                    and left["query_end"]
                    == rows[indexes[position]]["query_start"]
                    and rows[indexes[member_positions[-1]]]["query_end"]
                    == right["query_start"]
                )
            if valid:
                left_anchor = _reference_edge_at_query_side(left, "right")
                right_anchor = _reference_edge_at_query_side(right, "left")
                group_number += 1
                group_id = str(first.get("dup_group") or (
                    f"I{group_number}:{first['query_contig']}:"
                    f"{first['query_start']}-"
                    f"{rows[indexes[member_positions[-1]]]['query_end']}"
                ))
                anchor_start = min(left_anchor, right_anchor)
                anchor_end = max(left_anchor, right_anchor)
                for member_position in member_positions:
                    member = rows[indexes[member_position]]
                    member["insertion_group"] = group_id
                    member["insertion_anchor_path"] = left["reference_path"]
                    member["insertion_anchor_start"] = anchor_start
                    member["insertion_anchor_end"] = anchor_end

            position = max(position + 1, cursor)


def resolve_insertion_anchors(rows: List[dict]) -> None:
    """Anchor alignment-supported insertion mappings through their neighbors.

    UNMAPPED query gaps are barriers and never supply or transmit an anchor.
    Only INSERTION rows backed by an elected alignment CIGAR or alternative
    source mapping receive anchors here.
    """
    _assign_duplication_groups(rows)
    _assign_insertion_anchor_intervals(rows)
    by_contig: Dict[str, List[int]] = defaultdict(list)
    for index, row in enumerate(rows):
        if isinstance(row["query_start"], int):
            by_contig[row["query_contig"]].append(index)
    for indexes in by_contig.values():
        # Keep UNMAPPED rows in coordinate order as barriers.  They cannot be
        # anchors, and an alignment-supported insertion must not search across
        # unsupported query sequence for a distant mapped row.
        ordered = list(indexes)
        ordered.sort(key=lambda i: (
            rows[i]["query_start"], rows[i]["query_end"], i,
        ))
        positions = {
            row_index: offset for offset, row_index in enumerate(ordered)
        }

        # Find the closest row on each side that already has a real reference
        # coordinate.  An unresolved INSERTION is not an anchor: allowing two
        # adjacent pieces of the same losing alignment to select one another
        # creates a dependency cycle and hides the mapped flank just outside
        # the insertion block.
        nearest_mapped_left: List[Optional[int]] = [None] * len(ordered)
        nearest_mapped_right: List[Optional[int]] = [None] * len(ordered)
        mapped = None
        for offset, candidate_index in enumerate(ordered):
            nearest_mapped_left[offset] = mapped
            candidate = rows[candidate_index]
            if candidate.get("kind") == "UNMAPPED":
                mapped = None
            elif (
                candidate.get("reference_path") != "."
                and isinstance(candidate.get("reference_start"), int)
            ):
                mapped = candidate_index
        mapped = None
        for offset in range(len(ordered) - 1, -1, -1):
            nearest_mapped_right[offset] = mapped
            candidate_index = ordered[offset]
            candidate = rows[candidate_index]
            if candidate.get("kind") == "UNMAPPED":
                mapped = None
            elif (
                candidate.get("reference_path") != "."
                and isinstance(candidate.get("reference_start"), int)
            ):
                mapped = candidate_index

        for row_index in ordered:
            row = rows[row_index]
            if row["kind"] != "INSERTION" or row["anchor"] != ".":
                continue
            offset = positions[row_index]
            left_index = nearest_mapped_left[offset]
            right_index = nearest_mapped_right[offset]
            has_alternative_source = row.get("alternative_path") not in {
                None, "", ".",
            }
            if has_alternative_source:
                # Alternative-source rows are anchored only as one grouped
                # duplication replacement by _assign_duplication_groups.
                # Never fall back to a single PA-priority neighbor here.
                continue
            neighbors = [
                index for index in (left_index, right_index)
                if index is not None
            ]
            if not neighbors:
                continue
            chosen = max(neighbors, key=lambda i: _row_owner(rows[i]).priority)
            neighbor = rows[chosen]
            side = (
                "right"
                if neighbor["query_end"] <= row["query_start"]
                else "left"
            )
            anchor = _reference_edge_at_query_side(neighbor, side)
            row["reference_path"] = neighbor["reference_path"]
            row["reference_start"] = anchor
            row["reference_end"] = anchor
            row["reference_strand"] = neighbor["reference_strand"]
            row["anchor"] = anchor
            row["dependency"] = "."

    # Resolve the functional dependency graph once with path compression.
    # The former repeated whole-row relaxation was O(N^2) for a long chain.
    resolved: Dict[int, bool] = {}
    for start_index, start_row in enumerate(rows):
        dependency = start_row.get("dependency", ".")
        if not isinstance(dependency, str) or not dependency.startswith("@"):
            continue
        path = []
        local_positions = {}
        current = start_index
        success = False
        while True:
            if current in resolved:
                success = resolved[current]
                break
            current_row = rows[current]
            if (
                current_row["reference_path"] != "."
                and isinstance(current_row["reference_start"], int)
            ):
                success = True
                resolved[current] = True
                break
            current_dependency = current_row.get("dependency", ".")
            if (
                not isinstance(current_dependency, str)
                or not current_dependency.startswith("@")
            ):
                resolved[current] = False
                break
            if current in local_positions:
                success = False
                break
            local_positions[current] = len(path)
            path.append(current)
            current = int(current_dependency[1:])

        if success:
            for row_index in reversed(path):
                row = rows[row_index]
                target_index = int(row["dependency"][1:])
                target = rows[target_index]
                if (
                    target["reference_path"] == "."
                    or not isinstance(target["reference_start"], int)
                ):
                    success = False
                    break
                side = (
                    "right"
                    if target["query_end"] <= row["query_start"]
                    else "left"
                )
                anchor = _reference_edge_at_query_side(target, side)
                row["reference_path"] = target["reference_path"]
                row["reference_start"] = anchor
                row["reference_end"] = anchor
                row["reference_strand"] = target["reference_strand"]
                row["anchor"] = anchor
                row["dependency"] = "."
                resolved[row_index] = True
        if not success:
            for row_index in path:
                resolved[row_index] = False


def _interval_text(path, start, end, strand) -> str:
    if (
        path in {None, "", "."}
        or not isinstance(start, int)
        or not isinstance(end, int)
    ):
        return "."
    if strand == ">":
        strand = "+"
    elif strand == "<":
        strand = "-"
    if strand not in {"+", "-"}:
        strand = "+"
    return f"{path}:{start}-{end}{strand}"


def _assignment_alignment_cigar(row: Mapping[str, object]) -> str:
    """Slice an elected interval from its original source graph CIGAR.

    The fallback below exists only for private compatibility tests without a
    source interval/CIGAR. Production assignments use one interval slice,
    including internal I/D and encoded-template traversals.
    """
    materialized = row.get("alignment_cigar")
    if isinstance(materialized, str) and materialized not in {"", "."}:
        return materialized
    if (
        row.get("kind") == "UNMAPPED"
        or row.get("operation") == "GROUPED_UNMAPPED"
    ):
        # A duplication group may absorb a short unmapped query gap between
        # its members. That sequence was never aligned; graphreftovcf renders
        # it from the query as part of the group's insertion.
        return "."
    atom = _row_owner(row)
    if row["kind"] == "DELETION":
        start = row.get("reference_start")
        end = row.get("reference_end")
        if not isinstance(start, int) or not isinstance(end, int) or end <= start:
            return "."
        path_length = max(atom.reference_path_length, end)
        return core._format_interval_gcigar(
            row.get("reference_strand", atom.reference_strand),
            row.get("reference_path", atom.reference_path), path_length,
            start, end, [core.CigarOp(end - start, "D", "")],
        )
    if isinstance(atom, AlignmentInterval):
        if atom.trimmed:
            # Clipped query tails are UNMAPPED and have no alignment CIGAR.
            return "."
        raise ValueError(
            f"{atom.source_line}: elected interval row is missing its "
            "materialized CIGAR slice"
        )
    query_start = row.get("query_start")
    query_end = row.get("query_end")
    source_query = core.parse_coord(atom.source_query_interval)
    if (
        atom.graph_cigar not in {None, "", "."}
        and source_query is not None
        and isinstance(query_start, int)
        and isinstance(query_end, int)
        and source_query.chrom == row.get("query_contig")
        and source_query.start <= query_start < query_end <= source_query.end
    ):
        if source_query.strand == "+":
            local_start = query_start - source_query.start
            local_end = query_end - source_query.start
        else:
            local_start = source_query.end - query_end
            local_end = source_query.end - query_start
        return core.slice_reference_cigar_by_query(
            atom.graph_cigar, local_start, local_end,
            getattr(atom, "source_query_name", atom.query_name),
            include_left_boundary_deletions=True,
            include_right_boundary_deletions=(
                local_end == source_query.end - source_query.start
            ),
        )
    fragments = []
    for piece_atom, query_start, query_end in row.get(
        "_alignment_pieces", [(atom, row.get("query_start"), row.get("query_end"))],
    ):
        if (
            not isinstance(query_start, int)
            or not isinstance(query_end, int)
            or query_end <= query_start
        ):
            continue
        whole = (
            query_start == piece_atom.query_start
            and query_end == piece_atom.query_end
        )
        if whole and piece_atom.piece_op is not None:
            # A whole main-path atom is exactly one operation on one
            # reference interval; build its segment without format/parse.
            fragments.append(core.ReferenceCigarSegment(
                piece_atom.query_name, 0, piece_atom.reference_strand,
                piece_atom.reference_path, piece_atom.reference_path_length,
                piece_atom.reference_start, piece_atom.reference_end, "",
                [piece_atom.piece_op], True,
                (
                    piece_atom.reference_strand, piece_atom.reference_path,
                    piece_atom.reference_path_length,
                    piece_atom.reference_start,
                ),
            ))
            continue
        alignment_cigar = _atom_alignment_cigar(piece_atom)
        if alignment_cigar == ".":
            continue
        if whole:
            fragment = alignment_cigar
        else:
            if piece_atom.query_strand == "+":
                local_start = query_start - piece_atom.query_start
                local_end = query_end - piece_atom.query_start
            else:
                local_start = piece_atom.query_end - query_end
                local_end = piece_atom.query_end - query_start
            fragment = core.slice_reference_cigar_by_query(
                alignment_cigar, local_start, local_end,
                piece_atom.query_name,
                include_left_boundary_deletions=True,
                include_right_boundary_deletions=(
                    local_end == piece_atom.query_end - piece_atom.query_start
                ),
            )
        fragments.extend(core.parse_reference_cigar_segments(
            fragment, piece_atom.query_name,
        ))
    merged_fragments = []
    for segment in fragments:
        if merged_fragments:
            previous = merged_fragments[-1]
            contiguous = (
                previous.end == segment.start
                if segment.direction == ">" else previous.start == segment.end
            )
            if (
                previous.is_main and segment.is_main
                and previous.path == segment.path
                and previous.direction == segment.direction
                and previous.path_len == segment.path_len
                and contiguous
            ):
                # Extend in place: copying the op list per piece made long
                # coalesced rows quadratic.
                previous.start = min(previous.start, segment.start)
                previous.end = max(previous.end, segment.end)
                previous.ops.extend(segment.ops)
                continue
        # Every fragment here owns a fresh op list (built above or freshly
        # parsed), so in-place extension cannot alias an atom's state.
        merged_fragments.append(segment)
    return (
        core.format_reference_cigar_segments(merged_fragments)
        if merged_fragments else "."
    )


def _atom_alignment_cigar(atom: _LegacyAtom) -> str:
    """Return one atom's anchored CIGAR text, deriving it from its op."""
    if atom.alignment_cigar:
        return atom.alignment_cigar
    if atom.piece_op is None:
        return "."
    return core._format_interval_gcigar(
        atom.reference_strand, atom.reference_path,
        atom.reference_path_length, atom.reference_start,
        atom.reference_end, [atom.piece_op],
    )


FIELDS = (
    "assignment_id", "pa_name", "query_interval", "type",
    "reference_interval", "reference_alignment_cigar",
    "alternative_interval", "alternative_alignment_cigar",
    "reference_pa", "operation", "source_line",
    "anchor", "dependency", "dup_group", "insertion_group",
    "insertion_anchor_interval", "region_tier", "trimmed", "encoded",
    "lift_stage", "secondary_lift", "lift_kind", "lift_role",
    "alignment_score", "resurrected_lift",
)


def write_rows(path: str, rows: Sequence[dict]) -> None:
    temporary = path + f".tmp.{os.getpid()}"
    with open(temporary, "w", newline="") as handle:
        writer = csv.DictWriter(handle, FIELDS, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        for row in sorted(rows, key=lambda x: (
            x["query_contig"],
            x["query_start"] if isinstance(x["query_start"], int) else math.inf,
            x["reference_path"],
            x["reference_start"] if isinstance(x["reference_start"], int) else math.inf,
            x["kind"],
        )):
            atom = _row_owner(row)
            alignment_cigar = _assignment_alignment_cigar(row)
            has_alternative = row.get("alternative_path") not in {None, "", "."}
            writer.writerow({
                "assignment_id": row["assignment_id"],
                "pa_name": row.get("pa_name", atom.query_name),
                "query_interval": _interval_text(
                    row["query_contig"], row["query_start"], row["query_end"],
                    row["query_strand"],
                ),
                "type": row["kind"],
                "reference_interval": _interval_text(
                    row["reference_path"], row["reference_start"],
                    row["reference_end"], row["reference_strand"],
                ),
                "reference_alignment_cigar": (
                    "." if has_alternative else alignment_cigar
                ),
                "alternative_interval": _interval_text(
                    row.get("alternative_path", "."),
                    row.get("alternative_start", "."),
                    row.get("alternative_end", "."),
                    row.get("alternative_strand", atom.reference_strand),
                ),
                "alternative_alignment_cigar": (
                    alignment_cigar if has_alternative else "."
                ),
                "reference_pa": row.get("reference_pa", atom.reference_name),
                "operation": row.get("operation", atom.op),
                "source_line": row.get("source_line", atom.source_line),
                "anchor": row.get("anchor", "."),
                "dependency": row.get("dependency", "."),
                "dup_group": row.get("dup_group", "."),
                "insertion_group": row.get("insertion_group", "."),
                "insertion_anchor_interval": _interval_text(
                    row.get("insertion_anchor_path", "."),
                    row.get("insertion_anchor_start", "."),
                    row.get("insertion_anchor_end", "."),
                    ">",
                ),
                "region_tier": row.get("region_tier", atom.region_tier),
                "trimmed": int(atom.trimmed),
                "encoded": int(row.get("encoded", atom.encoded)),
                "lift_stage": row.get("lift_stage", atom.lift_stage),
                "secondary_lift": int(row.get(
                    "secondary_lift", atom.secondary_lift,
                )),
                "lift_kind": row.get(
                    "lift_kind", getattr(atom, "lift_kind", ""),
                ),
                "lift_role": row.get(
                    "lift_role", getattr(atom, "lift_role", ""),
                ),
                "resurrected_lift": int(row.get(
                    "resurrected_lift", atom.resurrected_lift,
                )),
                "alignment_score": f"{float(row.get('alignment_score', atom.score)):.6f}",
            })
    os.replace(temporary, path)


def build_assignments(
    input_path: str,
    coord_maps: Sequence[str],
    output_path: str,
    max_impute: int,
) -> Mapping[str, int]:
    """Build and write assignments for gap-fill and the standalone CLI."""
    # Election creates millions of acyclic objects; with the cyclic
    # collector active, each full collection rescans every live interval.
    gc_was_enabled = gc.isenabled()
    gc.disable()
    try:
        return _build_assignments(
            input_path, coord_maps, output_path, max_impute,
        )
    finally:
        if gc_was_enabled:
            gc.enable()


def _build_assignments(
    input_path: str,
    coord_maps: Sequence[str],
    output_path: str,
    max_impute: int,
) -> Mapping[str, int]:
    started = time.monotonic()
    lifts = read_lifts(coord_maps)
    intervals = read_intervals(input_path, lifts, max_impute)
    if any(item.lift_kind for item in intervals):
        policy = (
            "electing location/main and sequence/alternative layers "
            "independently; alternative sources do not reserve main query "
            "or reference intervals"
        )
    else:
        policy = (
            "electing primary stages 1>2>3, then secondary stages 2>3>4 "
            "over unresolved query intervals"
        )
    sys.stderr.write(
        "[pseudolinear] indexed "
        f"{len(lifts)} GenomeLift row(s) and {len(intervals)} alignment "
        f"interval(s); {policy}\n"
    )
    sys.stderr.flush()
    rows = interval_assignment_rows(intervals)
    elected = time.monotonic()
    write_rows(output_path, rows)
    counts = defaultdict(int)
    for row in rows:
        counts[row["kind"]] += 1
    sys.stderr.write(
        "[pseudolinear] wrote "
        f"{len(rows)} assignment row(s) in {time.monotonic() - started:.2f}s "
        f"(election {elected - started:.2f}s)\n"
    )
    sys.stderr.flush()
    return counts


def main(argv: Optional[Sequence[str]] = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("-i", "--input", required=True, help="post-gapfill graphcigartoreffix TSV")
    parser.add_argument("-m", "--coord-map", required=True, action="append", help="GenomeLift TSV; repeat for sparse overrides in replacement order")
    parser.add_argument("-o", "--output", required=True, help="exclusive pseudo-linear assignment TSV")
    parser.add_argument("--max-impute", type=int, default=10000)
    args = parser.parse_args(argv)
    counts = build_assignments(
        args.input, args.coord_map, args.output, args.max_impute,
    )
    print(
        "[pseudolinear] " + ", ".join(f"{key.lower()}={counts[key]}" for key in ("PRIMARY", "INSERTION", "DELETION", "UNMAPPED")),
        file=sys.stderr,
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
