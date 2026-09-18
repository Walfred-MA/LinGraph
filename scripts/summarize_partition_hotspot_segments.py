#!/usr/bin/env python3
"""Project graph source BED intervals through existing graph alignments.

All sample rows and source haplotypes participate by default. Use ``--onlyref``
to retain the historical reference-only selection.
"""

from __future__ import annotations

import argparse
import dataclasses
import heapq
import logging
import os
import re
from collections import defaultdict
from concurrent.futures import ProcessPoolExecutor, as_completed
from pathlib import Path
from typing import (
    DefaultDict, Dict, Iterable, List, Mapping, Optional, Sequence, Tuple,
)

from align_partition_hotspots import (
    OutputRow,
    read_graph_list,
)


LOG = logging.getLogger("summarize_partition_hotspot_segments")
DEFAULT_SEGMENT_CUTOFF = 15_000
DEFAULT_MINIMUM_SEGMENT_SIZE = 1_000
DEFAULT_VALID_ALIGNMENT_SIZE = 100
DEFAULT_VALID_REGION_MERGE_GAP = 1_000
DEFAULT_NEIGHBOR_BAD_SIZE_CUTOFF = 500
DEFAULT_NEIGHBOR_BAD_SIZE_PENALTY = 10
SOURCE_COORDINATE_RE = re.compile(
    r"^([^:]+):([^:]+):(\d+)-(\d+)([+-])$"
)
PATH_TOKEN_RE = re.compile(r"([><])([^><]+)")
CIGAR_TOKEN_RE = re.compile(r"([><])([^:<>]+):([^<>]+)")
CIGAR_OPERATION_RE = re.compile(r"(\d+)([HID=XMN])([A-Za-z]*)")
ALIGNMENT_HEADER_FIELDS = (
    "hotspot_index", "contig", "start", "end", "strand", "queryname",
    "graph_path", "graphcigar", "refpositions", "qpositions",
)
SEGMENT_HEADER = (
    "hotspot_index\tgraph_name\tquery_contig\tquery_start\tquery_end\t"
    "strand\tgraph_start\tgraph_end\tcomponent_count\n"
)


@dataclasses.dataclass(frozen=True)
class OriginalBedInterval:
    haplotype: str
    contig: str
    start: int
    end: int
    strand: str


@dataclasses.dataclass(frozen=True)
class GraphPathRecord:
    name: str
    haplotype: str
    contig: str
    start: int
    end: int
    strand: str
    length: int


@dataclasses.dataclass(frozen=True)
class ProjectedInterval:
    hotspot_index: int
    graph_name: str
    query_contig: str
    query_start: int
    query_end: int
    strand: str
    graph_start: int
    graph_end: int

    @property
    def size(self) -> int:
        return self.graph_end - self.graph_start


@dataclasses.dataclass(frozen=True)
class SegmentMappingPiece:
    query_start: int
    query_end: int
    strand: str
    graph_start: int
    graph_end: int


@dataclasses.dataclass
class FinalSegment:
    hotspot_index: int
    graph_name: str
    query_contig: str
    query_start: int
    query_end: int
    strand: str
    graph_start: int
    graph_end: int
    component_count: int = 1
    gap_duplication_size: int = 0
    bad_size: int = 0
    active: bool = True
    mapping_pieces: List[SegmentMappingPiece] = dataclasses.field(
        default_factory=list, repr=False, compare=False,
    )

    @property
    def size(self) -> int:
        return self.graph_end - self.graph_start

    @property
    def query_size(self) -> int:
        return self.query_end - self.query_start

    def line(self) -> str:
        return "\t".join(map(str, (
            self.hotspot_index, self.graph_name, self.query_contig,
            self.query_start, self.query_end, self.strand,
            self.graph_start, self.graph_end, self.component_count,
        ))) + "\n"


@dataclasses.dataclass
class RedundancyRun:
    start: int
    end: int
    redundancy: int

    @property
    def size(self) -> int:
        return self.end - self.start


@dataclasses.dataclass(frozen=True)
class NeighborMergeCandidate:
    start_index: int
    end_index: int
    score: int
    bad_size: int


@dataclasses.dataclass(frozen=True)
class AlignmentPathComponent:
    marker: str
    path: str
    cigar: str
    query_start: int
    query_end: int
    path_length: int


@dataclasses.dataclass(frozen=True)
class GraphMapping:
    path: str
    graph_start: int
    graph_end: int
    query_start: int
    query_end: int
    component_index: int = 0

    @property
    def size(self) -> int:
        return min(
            self.graph_end - self.graph_start,
            self.query_end - self.query_start,
        )


def read_alignment_output(path: str) -> List[OutputRow]:
    """Read headered or headerless align_partition_hotspots.py output."""
    output: List[OutputRow] = []
    with open(path, "rt") as handle:
        for line_number, raw in enumerate(handle, 1):
            if not raw.strip() or raw.startswith("#"):
                continue
            fields = raw.rstrip("\r\n").split("\t")
            if tuple(fields) == ALIGNMENT_HEADER_FIELDS:
                continue
            if len(fields) != len(ALIGNMENT_HEADER_FIELDS):
                raise ValueError(
                    f"{path}:{line_number}: expected ten alignment columns, "
                    f"found {len(fields)}"
                )
            try:
                hotspot_index, start, end = (
                    int(fields[0]), int(fields[2]), int(fields[3]),
                )
            except ValueError as error:
                raise ValueError(
                    f"{path}:{line_number}: invalid integer alignment field"
                ) from error
            if (
                hotspot_index < 1 or start < 0 or end <= start
                or fields[4] not in {"+", "-"}
            ):
                raise ValueError(
                    f"{path}:{line_number}: invalid alignment coordinates/strand"
                )
            output.append(OutputRow(
                hotspot_index, fields[1], start, end, fields[4], fields[5],
                fields[6], fields[7], fields[8], fields[9],
            ))
    return output


def query_name_matches_haplotype(
    query_name: str, haplotypes: Iterable[str],
) -> bool:
    """Return whether the generated query name contains an exact haplotype.

    Direct partition alignments name records as
    ``PARTITION_HAPLOTYPE_INDEX`` while the standalone hotspot aligner uses
    ``HAPLOTYPE_INDEX``.  Haplotype names themselves may contain underscores,
    so split-on-underscore parsing is not reliable.
    """
    for raw_haplotype in haplotypes:
        haplotype = str(raw_haplotype).strip()
        if not haplotype:
            continue
        if (
            query_name == haplotype
            or query_name.startswith(haplotype + "_")
            or query_name.endswith("_" + haplotype)
            or ("_" + haplotype + "_") in query_name
        ):
            return True
    return False


def alignment_row_matches_haplotype(
    row: OutputRow, haplotypes: Iterable[str],
) -> bool:
    return query_name_matches_haplotype(row.query_name, haplotypes)


def select_alignment_rows(
    rows: Sequence[OutputRow],
    reference_haplotypes: Iterable[str],
    only_reference: bool = False,
) -> List[OutputRow]:
    """Use the complete cohort by default or only named reference rows."""
    selected = list(rows)
    if not only_reference:
        return selected
    accepted = tuple(reference_haplotypes)
    if not accepted:
        raise ValueError("--onlyref requires a reference haplotype")
    return [
        row for row in selected
        if alignment_row_matches_haplotype(row, accepted)
    ]


def graph_name_from_fasta(path: str) -> str:
    name = os.path.basename(path)
    return name[:-3] if name.lower().endswith(".fa") else name


def original_bed_path(graph_fasta: str) -> str:
    return os.path.join(
        os.path.dirname(graph_fasta), graph_name_from_fasta(graph_fasta) + ".bed",
    )


def read_original_bed(
    graph_fasta: str,
    accepted_haplotypes: Optional[Iterable[str]] = None,
) -> List[OriginalBedInterval]:
    """Read original graph intervals, optionally filtering haplotypes."""
    bed_path = original_bed_path(graph_fasta)
    if not os.path.isfile(bed_path):
        LOG.warning(
            "Skipping %s: original BED is missing: %s",
            graph_name_from_fasta(graph_fasta), bed_path,
        )
        return []
    accepted = (
        None if accepted_haplotypes is None
        else set(accepted_haplotypes)
    )
    output: List[OriginalBedInterval] = []
    with open(bed_path, "rt") as handle:
        for line_number, raw in enumerate(handle, 1):
            if not raw.strip() or raw.startswith("#"):
                continue
            fields = raw.rstrip("\n").split("\t")
            if len(fields) < 7:
                fields = raw.split()
            if len(fields) < 7:
                raise ValueError(
                    f"{bed_path}:{line_number}: expected at least seven BED columns"
                )
            try:
                start, end = int(fields[1]), int(fields[2])
            except ValueError as error:
                raise ValueError(
                    f"{bed_path}:{line_number}: invalid BED coordinates"
                ) from error
            strand, haplotype = fields[5], fields[6]
            if start < 0 or end <= start or strand not in {"+", "-"}:
                raise ValueError(
                    f"{bed_path}:{line_number}: invalid interval "
                    f"{fields[0]}:{start}-{end}{strand}"
                )
            if accepted is None or haplotype in accepted:
                output.append(OriginalBedInterval(
                    haplotype, fields[0], start, end, strand,
                ))
    return output


def read_original_reference_bed(
    graph_fasta: str,
    reference_haplotypes: Iterable[str],
) -> List[OriginalBedInterval]:
    """Compatibility wrapper selecting only named reference templates."""
    return read_original_bed(graph_fasta, reference_haplotypes)


def read_graph_path_records(graph_fasta: str) -> Dict[str, GraphPathRecord]:
    """Read graph FASTA path lengths and public source coordinates."""
    output: Dict[str, GraphPathRecord] = {}
    current_name: Optional[str] = None
    current_source: Optional[Tuple[str, str, int, int, str]] = None
    current_length = 0

    def finish() -> None:
        nonlocal current_name, current_source, current_length
        if current_name is None or current_source is None:
            return
        haplotype, contig, start, end, strand = current_source
        if end - start != current_length:
            raise ValueError(
                f"{graph_fasta}: path {current_name!r} has {current_length} "
                f"bases but source interval {contig}:{start}-{end} has "
                f"{end - start} bases"
            )
        if current_name in output:
            raise ValueError(
                f"{graph_fasta}: duplicate graph path {current_name!r}"
            )
        output[current_name] = GraphPathRecord(
            current_name, haplotype, contig, start, end, strand,
            current_length,
        )

    with open(graph_fasta, "rt") as handle:
        for line_number, raw in enumerate(handle, 1):
            if raw.startswith(">"):
                finish()
                fields = raw[1:].strip().split()
                if not fields:
                    raise ValueError(f"{graph_fasta}:{line_number}: empty header")
                current_name = fields[0]
                current_source = None
                current_length = 0
                if len(fields) >= 2:
                    match = SOURCE_COORDINATE_RE.fullmatch(fields[1])
                    if match is not None:
                        haplotype, contig, start, end, strand = match.groups()
                        current_source = (
                            haplotype, contig, int(start), int(end), strand,
                        )
            elif current_name is None:
                if raw.strip():
                    raise ValueError(
                        f"{graph_fasta}:{line_number}: sequence before first header"
                    )
            else:
                current_length += len(raw.strip())
    finish()
    return output


def parse_position_list(
    text: str, context: str, *, allow_empty: bool = False,
) -> List[Tuple[int, int]]:
    if not text or text == ".":
        return []
    output: List[Tuple[int, int]] = []
    for value in text.split(";"):
        fields = value.split("_")
        if len(fields) != 2:
            raise ValueError(f"{context}: invalid position interval {value!r}")
        try:
            start, end = map(int, fields)
        except ValueError as error:
            raise ValueError(
                f"{context}: invalid position interval {value!r}"
            ) from error
        if start < 0 or end < start or (end == start and not allow_empty):
            raise ValueError(f"{context}: invalid position interval {value!r}")
        output.append((start, end))
    return output


def parse_cigar_operations(
    cigar: str, context: str,
) -> List[Tuple[int, str]]:
    matches = CIGAR_OPERATION_RE.findall(cigar)
    if not matches or "".join(
        f"{length}{operation}{payload}"
        for length, operation, payload in matches
    ) != cigar:
        raise ValueError(f"{context}: invalid graph CIGAR {cigar!r}")
    output = []
    for length_text, operation, payload in matches:
        length = int(length_text)
        if payload and operation not in {"I", "X"}:
            raise ValueError(
                f"{context}: unexpected sequence payload on {operation}"
            )
        if payload and len(payload) != length:
            raise ValueError(
                f"{context}: {operation} payload length {len(payload)} "
                f"does not match operation length {length}"
            )
        output.append((length, operation))
    return output


def parse_alignment_path_components(
    row: OutputRow, paths: Mapping[str, GraphPathRecord],
) -> List[AlignmentPathComponent]:
    """Parse synchronized graph-CIGAR components using FASTA path lengths."""
    if not row.graph_path or row.graph_path == "*":
        return []
    path_tokens = PATH_TOKEN_RE.findall(row.graph_path)
    cigar_tokens = CIGAR_TOKEN_RE.findall(row.graph_cigar)
    ref_positions = parse_position_list(
        row.ref_positions, f"hotspot {row.hotspot_index}:refpositions",
        allow_empty=True,
    )
    query_positions = parse_position_list(
        row.query_positions, f"hotspot {row.hotspot_index}:qpositions",
        allow_empty=True,
    )
    if len({
        len(path_tokens), len(cigar_tokens), len(ref_positions),
        len(query_positions),
    }) != 1:
        raise ValueError(
            f"hotspot {row.hotspot_index}: graph path/CIGAR/position column "
            "counts do not match"
        )
    output = []
    for index, (
        (marker, path_name), (cigar_marker, cigar_name, cigar),
        ref_position, query_position,
    ) in enumerate(zip(
        path_tokens, cigar_tokens, ref_positions, query_positions,
    ), 1):
        if (marker, path_name) != (cigar_marker, cigar_name):
            raise ValueError(
                f"hotspot {row.hotspot_index} component {index}: graph path "
                "and graph CIGAR name/orientation disagree"
            )
        path = paths.get(path_name)
        if path is None or ref_position[1] == ref_position[0]:
            continue
        operations = parse_cigar_operations(
            cigar, f"hotspot {row.hotspot_index} component {index}",
        )
        query_consumed = sum(
            length for length, operation in operations
            if operation in {"I", "=", "X", "M"}
        )
        if query_consumed != query_position[1] - query_position[0]:
            raise ValueError(
                f"hotspot {row.hotspot_index} component {index}: CIGAR "
                f"consumes {query_consumed} query bases but qpositions spans "
                f"{query_position[1] - query_position[0]}"
            )
        output.append(AlignmentPathComponent(
            marker, path_name, cigar, query_position[0],
            query_position[1], path.length,
        ))
    return output


def absolute_to_oriented_query(
    row: OutputRow, start: int, end: int,
) -> Tuple[int, int]:
    if row.strand == "+":
        return start - row.start, end - row.start
    return row.end - end, row.end - start


def oriented_to_absolute_query(
    row: OutputRow, start: int, end: int,
) -> Tuple[int, int]:
    if row.strand == "+":
        return row.start + start, row.start + end
    return row.end - end, row.end - start


def query_interval_to_graph_mappings(
    row: OutputRow,
    components: Sequence[AlignmentPathComponent],
    start: int,
    end: int,
) -> List[GraphMapping]:
    """Map an absolute query interval to atomic matched graph-path spans."""
    clipped_start = max(start, row.start)
    clipped_end = min(end, row.end)
    if clipped_end <= clipped_start:
        return []
    target_start, target_end = absolute_to_oriented_query(
        row, clipped_start, clipped_end,
    )
    output = []
    for component_index, component in enumerate(components):
        if (
            component.query_end <= target_start
            or component.query_start >= target_end
        ):
            continue
        graph_cursor = 0 if component.marker == ">" else component.path_length
        query_cursor = component.query_start
        for length, operation in parse_cigar_operations(
            component.cigar, f"graph path {component.path}",
        ):
            consumes_graph = operation in {"H", "D", "=", "X", "M", "N"}
            consumes_query = operation in {"I", "=", "X", "M"}
            if operation in {"=", "X", "M"}:
                operation_start = query_cursor
                operation_end = query_cursor + length
                left = max(target_start, operation_start)
                right = min(target_end, operation_end)
                if right > left:
                    if component.marker == ">":
                        graph_start = graph_cursor + left - operation_start
                        graph_end = graph_cursor + right - operation_start
                    else:
                        graph_start = graph_cursor - (right - operation_start)
                        graph_end = graph_cursor - (left - operation_start)
                    query_start, query_end = oriented_to_absolute_query(
                        row, left, right,
                    )
                    output.append(GraphMapping(
                        component.path, graph_start, graph_end,
                        query_start, query_end, component_index,
                    ))
            if consumes_graph:
                graph_cursor += length if component.marker == ">" else -length
            if consumes_query:
                query_cursor += length
    return output


def graph_interval_to_query_mappings(
    row: OutputRow,
    component: AlignmentPathComponent,
    graph_start: int,
    graph_end: int,
) -> List[GraphMapping]:
    """Map one canonical graph-path span back to atomic query intervals."""
    if graph_start < 0 or graph_end <= graph_start:
        return []
    output = []
    graph_cursor = 0 if component.marker == ">" else component.path_length
    query_cursor = component.query_start
    for length, operation in parse_cigar_operations(
        component.cigar, f"graph path {component.path}",
    ):
        consumes_graph = operation in {"H", "D", "=", "X", "M", "N"}
        consumes_query = operation in {"I", "=", "X", "M"}
        if operation in {"=", "X", "M"}:
            if component.marker == ">":
                operation_start, operation_end = (
                    graph_cursor, graph_cursor + length,
                )
                left = max(graph_start, operation_start)
                right = min(graph_end, operation_end)
                local_start = query_cursor + left - operation_start
                local_end = query_cursor + right - operation_start
            else:
                operation_start, operation_end = (
                    graph_cursor - length, graph_cursor,
                )
                left = max(graph_start, operation_start)
                right = min(graph_end, operation_end)
                local_start = query_cursor + operation_end - right
                local_end = query_cursor + operation_end - left
            if right > left:
                query_start, query_end = oriented_to_absolute_query(
                    row, local_start, local_end,
                )
                output.append(GraphMapping(
                    component.path, left, right, query_start, query_end,
                ))
        if consumes_graph:
            graph_cursor += length if component.marker == ">" else -length
        if consumes_query:
            query_cursor += length
    return output


def path_local_bed_interval(
    path: GraphPathRecord, bed: OriginalBedInterval,
) -> Optional[Tuple[int, int]]:
    if path.haplotype != bed.haplotype or path.contig != bed.contig:
        return None
    source_start = max(path.start, bed.start)
    source_end = min(path.end, bed.end)
    if source_end <= source_start:
        return None
    if path.strand == "+":
        return source_start - path.start, source_end - path.start
    return path.end - source_end, path.end - source_start


def project_cigar_to_path_interval(
    cigar: str,
    marker: str,
    path_length: int,
    query_start: int,
    query_end: int,
    local_start: int,
    local_end: int,
    context: str,
) -> Optional[Tuple[int, int, int, int]]:
    """Project one canonical graph-path interval through a graph CIGAR."""
    operations = parse_cigar_operations(cigar, context)
    ref_cursor = 0 if marker == ">" else path_length
    query_cursor = query_start
    query_pieces: List[Tuple[int, int]] = []
    graph_pieces: List[Tuple[int, int]] = []
    for length, operation in operations:
        consumes_reference = operation in {"H", "D", "=", "X", "M", "N"}
        consumes_query = operation in {"I", "=", "X", "M"}
        if operation in {"=", "X", "M"} and length:
            if marker == ">":
                operation_start, operation_end = ref_cursor, ref_cursor + length
                left = max(local_start, operation_start)
                right = min(local_end, operation_end)
                if right > left:
                    query_pieces.append((
                        query_cursor + left - operation_start,
                        query_cursor + right - operation_start,
                    ))
                    graph_pieces.append((left, right))
            else:
                operation_start, operation_end = ref_cursor - length, ref_cursor
                left = max(local_start, operation_start)
                right = min(local_end, operation_end)
                if right > left:
                    query_pieces.append((
                        query_cursor + operation_end - right,
                        query_cursor + operation_end - left,
                    ))
                    graph_pieces.append((left, right))
        if consumes_reference:
            ref_cursor += length if marker == ">" else -length
        if consumes_query:
            query_cursor += length
    if query_cursor != query_end:
        raise ValueError(
            f"{context}: CIGAR consumes {query_cursor - query_start} query "
            f"bases but qpositions spans {query_end - query_start}"
        )
    if not query_pieces:
        return None
    return (
        min(start for start, _end in query_pieces),
        max(end for _start, end in query_pieces),
        min(start for start, _end in graph_pieces),
        max(end for _start, end in graph_pieces),
    )


def source_interval_from_local(
    path: GraphPathRecord, local_start: int, local_end: int,
) -> Tuple[int, int]:
    if path.strand == "+":
        return path.start + local_start, path.start + local_end
    return path.end - local_end, path.end - local_start


def projected_strand(
    query_strand: str,
    path_marker: str,
    path_source_strand: str,
    bed_strand: str,
) -> str:
    direction = 1
    for value in (
        query_strand, "+" if path_marker == ">" else "-",
        path_source_strand, bed_strand,
    ):
        direction *= 1 if value == "+" else -1
    return "+" if direction == 1 else "-"


def project_alignment_row(
    row: OutputRow,
    graph_name: str,
    paths: Mapping[str, GraphPathRecord],
    beds: Sequence[OriginalBedInterval],
) -> List[ProjectedInterval]:
    if row.graph_path == "*" or not beds:
        return []
    path_tokens = PATH_TOKEN_RE.findall(row.graph_path)
    cigar_tokens = CIGAR_TOKEN_RE.findall(row.graph_cigar)
    ref_positions = parse_position_list(
        row.ref_positions, f"hotspot {row.hotspot_index}:refpositions",
        allow_empty=True,
    )
    query_positions = parse_position_list(
        row.query_positions, f"hotspot {row.hotspot_index}:qpositions",
        allow_empty=True,
    )
    if len({
        len(path_tokens), len(cigar_tokens), len(ref_positions),
        len(query_positions),
    }) != 1:
        raise ValueError(
            f"hotspot {row.hotspot_index}: graph path/CIGAR/position column "
            "counts do not match"
        )
    output: List[ProjectedInterval] = []
    seen = set()
    for component_index, (
        path_token, cigar_token, _ref_position, query_position,
    ) in enumerate(zip(
        path_tokens, cigar_tokens, ref_positions, query_positions,
    ), 1):
        marker, path_name = path_token
        cigar_marker, cigar_name, cigar = cigar_token
        if (marker, path_name) != (cigar_marker, cigar_name):
            raise ValueError(
                f"hotspot {row.hotspot_index} component {component_index}: "
                "graph path and graph CIGAR name/orientation disagree"
            )
        # Older graph-CIGAR files can retain insertion-only or deletion-only
        # components.  Their synchronized coordinate entry is zero-width.
        # Keep accepting the row, but such a component cannot define a block.
        if (
            _ref_position[1] == _ref_position[0]
            or query_position[1] == query_position[0]
        ):
            continue
        path = paths.get(path_name)
        if path is None:
            continue
        for bed in beds:
            local_bed = path_local_bed_interval(path, bed)
            if local_bed is None:
                continue
            projected = project_cigar_to_path_interval(
                cigar, marker, path.length,
                query_position[0], query_position[1],
                local_bed[0], local_bed[1],
                f"hotspot {row.hotspot_index} component {component_index}",
            )
            if projected is None:
                continue
            oriented_start, oriented_end, local_start, local_end = projected
            if row.strand == "+":
                query_start = row.start + oriented_start
                query_end = row.start + oriented_end
            else:
                query_start = row.end - oriented_end
                query_end = row.end - oriented_start
            graph_start, graph_end = source_interval_from_local(
                path, local_start, local_end,
            )
            graph_start = max(graph_start, bed.start)
            graph_end = min(graph_end, bed.end)
            if query_end <= query_start or graph_end <= graph_start:
                continue
            strand = projected_strand(
                row.strand, marker, path.strand, bed.strand,
            )
            key = (query_start, query_end, strand, graph_start, graph_end)
            if key in seen:
                continue
            seen.add(key)
            output.append(ProjectedInterval(
                row.hotspot_index, graph_name, row.contig,
                query_start, query_end, strand, graph_start, graph_end,
            ))
    return output


def interval_overlap(
    left_start: int, left_end: int, right_start: int, right_end: int,
) -> int:
    return max(0, min(left_end, right_end) - max(left_start, right_start))


def merge_intervals_with_gap(
    intervals: Iterable[Tuple[int, int]], gap: int,
) -> List[Tuple[int, int]]:
    """Union half-open intervals, merging positive gaps strictly below GAP."""
    if gap < 0:
        raise ValueError("interval merge gap must be nonnegative")
    ordered = sorted(set(
        (start, end) for start, end in intervals if end > start
    ))
    if not ordered:
        return []
    merged = [list(ordered[0])]
    for start, end in ordered[1:]:
        distance = start - merged[-1][1]
        if distance <= 0 or distance < gap:
            merged[-1][1] = max(merged[-1][1], end)
        else:
            merged.append([start, end])
    return [(start, end) for start, end in merged]


def build_graph_valid_regions(
    rows: Sequence[OutputRow],
    paths: Mapping[str, GraphPathRecord],
    beds: Sequence[OriginalBedInterval],
    minimum_alignment_size: int = DEFAULT_VALID_ALIGNMENT_SIZE,
    merge_gap: int = DEFAULT_VALID_REGION_MERGE_GAP,
) -> Dict[str, List[Tuple[int, int]]]:
    """Build reference-coordinate validity masks through graph spans.

    The original BED intervals are first projected from the reference query
    onto graph paths. A CIGAR component contributes all of its atomic matched
    pieces only when their accumulated query union is strictly longer than
    ``minimum_alignment_size``. Unioned graph spans are projected through every
    alignment row with the same component-level rule. Candidate query intervals
    and original BED intervals are then jointly merged across positive gaps
    strictly smaller than ``merge_gap``; only merged components connected to an
    original BED interval survive.
    """
    if minimum_alignment_size < 0:
        raise ValueError("minimum valid-region alignment size must be nonnegative")
    if merge_gap < 0:
        raise ValueError("valid-region merge gap must be nonnegative")

    defined_by_contig: DefaultDict[str, List[Tuple[int, int]]] = defaultdict(list)
    for bed in beds:
        defined_by_contig[bed.contig].append((bed.start, bed.end))
    defined = {
        contig: merge_intervals_with_gap(intervals, 0)
        for contig, intervals in defined_by_contig.items()
    }
    if not defined:
        return {}

    parsed_rows = [
        (row, parse_alignment_path_components(row, paths))
        for row in rows
    ]
    graph_spans: DefaultDict[str, List[Tuple[int, int]]] = defaultdict(list)
    for row, components in parsed_rows:
        row_mappings: List[GraphMapping] = []
        for start, end in defined.get(row.contig, ()):
            if end <= row.start or start >= row.end:
                continue
            row_mappings.extend(query_interval_to_graph_mappings(
                row, components, start, end,
            ))
        by_component: DefaultDict[int, List[GraphMapping]] = defaultdict(list)
        for mapping in row_mappings:
            by_component[mapping.component_index].append(mapping)
        for mappings in by_component.values():
            valid_query_size = sum(
                end - start for start, end in merge_intervals_with_gap(
                    (
                        (mapping.query_start, mapping.query_end)
                        for mapping in mappings
                    ),
                    0,
                )
            )
            if valid_query_size > minimum_alignment_size:
                for mapping in mappings:
                    graph_spans[mapping.path].append((
                        mapping.graph_start, mapping.graph_end,
                    ))
    graph_spans = defaultdict(list, {
        path: merge_intervals_with_gap(intervals, 0)
        for path, intervals in graph_spans.items()
    })

    candidates: DefaultDict[str, List[Tuple[int, int]]] = defaultdict(list)
    for row, components in parsed_rows:
        for component in components:
            component_mappings: List[GraphMapping] = []
            for graph_start, graph_end in graph_spans.get(component.path, ()):
                component_mappings.extend(graph_interval_to_query_mappings(
                    row, component, graph_start, graph_end,
                ))
            valid_query_size = sum(
                end - start for start, end in merge_intervals_with_gap(
                    (
                        (mapping.query_start, mapping.query_end)
                        for mapping in component_mappings
                    ),
                    0,
                )
            )
            if valid_query_size > minimum_alignment_size:
                for mapping in component_mappings:
                    candidates[row.contig].append((
                        mapping.query_start, mapping.query_end,
                    ))

    valid = {}
    for contig, original_intervals in defined.items():
        merged_with_original = merge_intervals_with_gap(
            [*candidates.get(contig, ()), *original_intervals], merge_gap,
        )
        valid[contig] = [
            interval for interval in merged_with_original
            if any(
                interval_overlap(interval[0], interval[1], start, end) > 0
                for start, end in original_intervals
            )
        ]
    return valid


def clip_segment_to_query_interval(
    segment: FinalSegment, start: int, end: int,
) -> Optional[FinalSegment]:
    """Clip one segment and proportionally retain its graph-coordinate span."""
    query_start = max(segment.query_start, start)
    query_end = min(segment.query_end, end)
    if query_end <= query_start:
        return None
    query_size = segment.query_end - segment.query_start
    graph_size = segment.graph_end - segment.graph_start
    if query_size <= 0 or graph_size <= 0:
        return None
    left = query_start - segment.query_start
    right = query_end - segment.query_start
    if segment.strand == "+":
        graph_start = segment.graph_start + left * graph_size // query_size
        graph_end = (
            segment.graph_start
            + (right * graph_size + query_size - 1) // query_size
        )
    else:
        graph_start = (
            segment.graph_end
            - (right * graph_size + query_size - 1) // query_size
        )
        graph_end = segment.graph_end - left * graph_size // query_size
    if graph_end <= graph_start:
        return None
    return dataclasses.replace(
        segment,
        query_start=query_start,
        query_end=query_end,
        graph_start=graph_start,
        graph_end=graph_end,
    )


def intersect_segments_with_valid_regions(
    segments: Sequence[FinalSegment],
    valid_regions: Mapping[str, Sequence[Tuple[int, int]]],
    minimum_segment_size: int,
) -> List[FinalSegment]:
    """Intersect every final block with its reference-contig validity mask."""
    output = []
    seen = set()
    for segment in segments:
        for start, end in valid_regions.get(segment.query_contig, ()):
            clipped = clip_segment_to_query_interval(segment, start, end)
            if clipped is None:
                continue
            if (
                minimum_segment_size > 0
                and clipped.query_end - clipped.query_start
                < minimum_segment_size
            ):
                continue
            key = (
                clipped.hotspot_index, clipped.graph_name,
                clipped.query_contig, clipped.query_start, clipped.query_end,
                clipped.strand, clipped.graph_start, clipped.graph_end,
                clipped.component_count,
            )
            if key not in seen:
                seen.add(key)
                output.append(clipped)
    output.sort(key=lambda segment: (
        segment.hotspot_index, segment.graph_name, segment.query_contig,
        segment.query_start, segment.query_end, segment.strand,
        segment.graph_start, segment.graph_end,
    ))
    return output


def validate_nonoverlapping_segments(
    segments: Sequence[FinalSegment],
) -> None:
    """Reject overlapping query blocks even when their strands differ."""
    grouped: DefaultDict[
        Tuple[int, str, str], List[FinalSegment]
    ] = defaultdict(list)
    for segment in segments:
        grouped[(
            segment.hotspot_index, segment.graph_name,
            segment.query_contig,
        )].append(segment)
    for key, members in grouped.items():
        ordered = sorted(members, key=lambda segment: (
            segment.query_start, segment.query_end, segment.strand,
        ))
        previous = ordered[0]
        for segment in ordered[1:]:
            if segment.query_start < previous.query_end:
                raise ValueError(
                    f"hotspot {key[0]} graph {key[1]} contig {key[2]}: "
                    f"overlapping final blocks {previous.query_start}-"
                    f"{previous.query_end}{previous.strand} and "
                    f"{segment.query_start}-{segment.query_end}"
                    f"{segment.strand}"
                )
            previous = segment


def merge_projected_intervals(
    intervals: Sequence[ProjectedInterval],
    scutoff: int = DEFAULT_SEGMENT_CUTOFF,
    dcutoff: int = DEFAULT_SEGMENT_CUTOFF,
) -> List[FinalSegment]:
    """Build final segments while retaining retired segments in the result."""
    if scutoff < 0 or dcutoff < 1:
        raise ValueError("scutoff must be nonnegative and dcutoff must be positive")
    grouped: DefaultDict[
        Tuple[int, str, str, str], List[ProjectedInterval]
    ] = defaultdict(list)
    for interval in intervals:
        grouped[(
            interval.hotspot_index, interval.graph_name,
            interval.query_contig, interval.strand,
        )].append(interval)

    final: List[FinalSegment] = []
    for group_key in sorted(grouped):
        active: List[FinalSegment] = []
        ordered = sorted(grouped[group_key], key=lambda item: (
            item.query_start, item.query_end,
            item.graph_start, item.graph_end,
        ))
        for interval in ordered:
            mapping_piece = SegmentMappingPiece(
                interval.query_start, interval.query_end, interval.strand,
                interval.graph_start, interval.graph_end,
            )
            active = [
                segment for segment in active
                if segment.active and segment.gap_duplication_size < dcutoff
            ]
            merged = False
            for active_index, previous in enumerate(tuple(active)):
                overlap = interval_overlap(
                    interval.graph_start, interval.graph_end,
                    previous.graph_start, previous.graph_end,
                )
                if (
                    overlap < scutoff
                    and overlap < 0.5 * interval.size
                    and previous.gap_duplication_size < min(
                        interval.size, previous.size,
                    )
                ):
                    previous.query_start = min(
                        previous.query_start, interval.query_start,
                    )
                    previous.query_end = max(
                        previous.query_end, interval.query_end,
                    )
                    previous.graph_start = min(
                        previous.graph_start, interval.graph_start,
                    )
                    previous.graph_end = max(
                        previous.graph_end, interval.graph_end,
                    )
                    previous.component_count += 1
                    previous.gap_duplication_size = 0
                    previous.mapping_pieces.append(mapping_piece)
                    for later in active[active_index + 1:]:
                        later.active = False
                    active = active[:active_index + 1]
                    merged = True
                    break

                accumulated = previous.gap_duplication_size + overlap
                if accumulated > previous.size or accumulated >= dcutoff:
                    previous.active = False
                else:
                    previous.gap_duplication_size = accumulated

            active = [segment for segment in active if segment.active]
            if not merged:
                segment = FinalSegment(
                    interval.hotspot_index, interval.graph_name,
                    interval.query_contig, interval.query_start,
                    interval.query_end, interval.strand,
                    interval.graph_start, interval.graph_end,
                    mapping_pieces=[mapping_piece],
                )
                final.append(segment)
                active.append(segment)
    final.sort(key=lambda segment: (
        segment.hotspot_index, segment.graph_name, segment.query_contig,
        segment.query_start, segment.query_end, segment.strand,
        segment.graph_start, segment.graph_end,
    ))
    return final


def coalesce_redundancy_runs(
    runs: Sequence[RedundancyRun],
) -> List[RedundancyRun]:
    """Join touching half-open runs that have the same redundancy."""
    output: List[RedundancyRun] = []
    for run in sorted(runs, key=lambda value: (
        value.start, value.end, value.redundancy,
    )):
        if run.end <= run.start:
            continue
        if (
            output and output[-1].end == run.start
            and output[-1].redundancy == run.redundancy
        ):
            output[-1].end = run.end
        else:
            output.append(RedundancyRun(
                run.start, run.end, run.redundancy,
            ))
    return output


def merge_coordinate_intervals(
    intervals: Sequence[Tuple[int, int]],
) -> List[Tuple[int, int]]:
    """Merge overlapping or touching positive half-open intervals."""
    merged: List[List[int]] = []
    for start, end in sorted(intervals):
        if end <= start:
            continue
        if merged and start <= merged[-1][1]:
            merged[-1][1] = max(merged[-1][1], end)
        else:
            merged.append([start, end])
    return [(start, end) for start, end in merged]


def redundancy_runs_from_segments(
    segments: Sequence[FinalSegment],
) -> List[RedundancyRun]:
    """Return constant query-depth runs over a group's complete span."""
    if not segments:
        return []
    normalized_domains = [(
        min(segment.query_start for segment in segments),
        max(segment.query_end for segment in segments),
    )]
    runs: List[RedundancyRun] = []
    for domain_start, domain_end in normalized_domains:
        if domain_end <= domain_start:
            continue
        events: DefaultDict[int, int] = defaultdict(int)
        events[domain_start] += 0
        events[domain_end] += 0
        for segment in segments:
            start = max(domain_start, segment.query_start)
            end = min(domain_end, segment.query_end)
            if end <= start:
                continue
            events[start] += 1
            events[end] -= 1
        coordinates = sorted(events)
        depth = 0
        for index, coordinate in enumerate(coordinates[:-1]):
            depth += events[coordinate]
            next_coordinate = coordinates[index + 1]
            if next_coordinate > coordinate:
                runs.append(RedundancyRun(
                    coordinate, next_coordinate, depth,
                ))
    return coalesce_redundancy_runs(runs)


def demote_short_redundancy_runs(
    runs: Sequence[RedundancyRun], minimum_size: int,
) -> List[RedundancyRun]:
    """Repeatedly demote undersized depth>1 runs and coalesce the result."""
    stabilized = coalesce_redundancy_runs(runs)
    if minimum_size == 0:
        return stabilized
    while True:
        changed = False
        demoted: List[RedundancyRun] = []
        for run in stabilized:
            redundancy = run.redundancy
            if redundancy > 1 and run.size < minimum_size:
                redundancy -= 1
                changed = True
            demoted.append(RedundancyRun(
                run.start, run.end, redundancy,
            ))
        stabilized = coalesce_redundancy_runs(demoted)
        if not changed:
            return stabilized


def recursively_merge_short_runs(
    runs: Sequence[RedundancyRun], minimum_size: int,
) -> List[RedundancyRun]:
    """Absorb the leftmost undersized run into an immediate neighbor.

    A qualifying large neighbor is preferred, followed by the closest
    redundancy level, the longer neighbor, and finally the left neighbor.
    Processing restarts after every merge. A short final run with no neighbor
    is discarded, making ``minimum_size`` a hard output-size floor.
    """
    merged = coalesce_redundancy_runs(runs)
    if minimum_size == 0:
        return merged
    while len(merged) > 1:
        short_index = next((
            index for index, run in enumerate(merged)
            if run.size < minimum_size
        ), None)
        if short_index is None:
            break
        neighbor_indexes = []
        if short_index > 0:
            neighbor_indexes.append(short_index - 1)
        if short_index + 1 < len(merged):
            neighbor_indexes.append(short_index + 1)
        if not neighbor_indexes:
            break
        short = merged[short_index]

        def priority(index: int) -> Tuple[int, int, int, int]:
            neighbor = merged[index]
            return (
                int(neighbor.size >= minimum_size),
                -abs(neighbor.redundancy - short.redundancy),
                neighbor.size,
                int(index < short_index),
            )

        neighbor_index = max(neighbor_indexes, key=priority)
        neighbor = merged[neighbor_index]
        neighbor.start = min(neighbor.start, short.start)
        neighbor.end = max(neighbor.end, short.end)
        del merged[short_index]
        merged = coalesce_redundancy_runs(merged)
    return [run for run in merged if run.size >= minimum_size]


def segment_mapping_pieces(
    segments: Sequence[FinalSegment],
) -> List[SegmentMappingPiece]:
    """Return preserved atomic mappings, with legacy aggregate fallbacks."""
    output = []
    for segment in segments:
        if segment.mapping_pieces:
            output.extend(segment.mapping_pieces)
        else:
            output.append(SegmentMappingPiece(
                segment.query_start, segment.query_end, segment.strand,
                segment.graph_start, segment.graph_end,
            ))
    return [piece for piece in output if (
        piece.query_end > piece.query_start
        and piece.graph_end > piece.graph_start
        and piece.strand in {"+", "-"}
    )]


def neighbor_group_statistics(
    segments: Sequence[FinalSegment],
) -> Tuple[int, bool]:
    """Return cumulative bad size and starting-strand dominance.

    Bad size is the query length covered by at least two atomic blocks plus
    the complete query length contributed by the strand opposite to the
    group's first block. These quantities are intentionally additive: an
    opposite-strand overlap is penalized once as overlap and once as minority
    strand. The first block's strand wins ties but must not be the minority.
    """
    if not segments:
        return 0, False
    pieces = segment_mapping_pieces(segments)
    if not pieces:
        return 0, False

    events: DefaultDict[int, int] = defaultdict(int)
    strand_sizes: DefaultDict[str, int] = defaultdict(int)
    for piece in pieces:
        events[piece.query_start] += 1
        events[piece.query_end] -= 1
        strand_sizes[piece.strand] += piece.query_end - piece.query_start

    overlap_size = 0
    depth = 0
    coordinates = sorted(events)
    for index, coordinate in enumerate(coordinates[:-1]):
        depth += events[coordinate]
        next_coordinate = coordinates[index + 1]
        if depth > 1:
            overlap_size += next_coordinate - coordinate

    starting_strand = segments[0].strand
    starting_size = strand_sizes[starting_strand]
    minority_size = sum(
        size for strand, size in strand_sizes.items()
        if strand != starting_strand
    )
    return overlap_size + minority_size, starting_size >= minority_size


def best_neighbor_merge_candidate(
    segments: Sequence[FinalSegment],
    bad_size_cutoff: int = DEFAULT_NEIGHBOR_BAD_SIZE_CUTOFF,
    bad_size_penalty: int = DEFAULT_NEIGHBOR_BAD_SIZE_PENALTY,
) -> Optional[NeighborMergeCandidate]:
    """Use the requested left-to-right DP to find the best active merge."""
    if bad_size_cutoff < 0 or bad_size_penalty < 0:
        raise ValueError("neighbor bad-size cutoff/penalty must be nonnegative")
    best: Optional[NeighborMergeCandidate] = None
    for start_index in range(len(segments) - 1):
        running_score = 0
        previous_bad_size, _dominant = neighbor_group_statistics(
            segments[start_index:start_index + 1],
        )
        running_start = segments[start_index].query_start
        running_end = segments[start_index].query_end
        for end_index in range(start_index + 1, len(segments)):
            next_segment = segments[end_index]
            smaller_block_size = min(
                running_end - running_start,
                next_segment.query_size,
            )
            candidate_bad_size, starting_strand_is_dominant = (
                neighbor_group_statistics(
                    segments[start_index:end_index + 1],
                )
            )
            # Bad size is monotone as blocks are appended, so no longer group
            # beginning at this position can become eligible after the cutoff.
            if candidate_bad_size > bad_size_cutoff:
                break
            added_bad_size = candidate_bad_size - previous_bad_size
            merging_score = (
                smaller_block_size - bad_size_penalty * added_bad_size
            )
            running_score += merging_score
            running_start = min(running_start, next_segment.query_start)
            running_end = max(running_end, next_segment.query_end)
            previous_bad_size = candidate_bad_size
            if not starting_strand_is_dominant or running_score <= 0:
                continue
            candidate = NeighborMergeCandidate(
                start_index, end_index, running_score, candidate_bad_size,
            )
            # Scanning order supplies the deterministic tie break: earliest
            # starting position, followed by the earliest ending position.
            if best is None or candidate.score > best.score:
                best = candidate
    return best


def merge_neighbor_candidate(
    segments: Sequence[FinalSegment],
    candidate: NeighborMergeCandidate,
) -> List[FinalSegment]:
    """Merge one selected candidate while retaining its atomic mappings."""
    members = list(segments[candidate.start_index:candidate.end_index + 1])
    first = members[0]
    pieces = segment_mapping_pieces(members)
    merged = FinalSegment(
        first.hotspot_index,
        first.graph_name,
        first.query_contig,
        min(segment.query_start for segment in members),
        max(segment.query_end for segment in members),
        first.strand,
        min(segment.graph_start for segment in members),
        max(segment.graph_end for segment in members),
        component_count=sum(segment.component_count for segment in members),
        bad_size=candidate.bad_size,
        mapping_pieces=pieces,
    )
    return [
        *segments[:candidate.start_index],
        merged,
        *segments[candidate.end_index + 1:],
    ]


def neighbor_join_segments(
    segments: Sequence[FinalSegment],
    bad_size_cutoff: int = DEFAULT_NEIGHBOR_BAD_SIZE_CUTOFF,
    bad_size_penalty: int = DEFAULT_NEIGHBOR_BAD_SIZE_PENALTY,
) -> List[FinalSegment]:
    """Iteratively merge the highest-scoring contiguous neighbor group."""
    ordered = sorted(segments, key=lambda segment: (
        segment.query_start, segment.query_end, segment.strand,
        segment.graph_start, segment.graph_end,
    ))
    while True:
        candidate = best_neighbor_merge_candidate(
            ordered, bad_size_cutoff, bad_size_penalty,
        )
        if candidate is None:
            return ordered
        ordered = merge_neighbor_candidate(ordered, candidate)


def project_mapping_piece(
    piece: SegmentMappingPiece, query_start: int, query_end: int,
) -> Tuple[int, int]:
    """Project or locally extrapolate one query piece to graph coordinates."""
    query_size = piece.query_end - piece.query_start
    graph_size = piece.graph_end - piece.graph_start
    if query_size <= 0 or graph_size <= 0 or query_end <= query_start:
        return query_start, query_end
    left = query_start - piece.query_start
    right = query_end - piece.query_start
    if piece.strand == "+":
        graph_start = piece.graph_start + left * graph_size // query_size
        graph_end = (
            piece.graph_start
            + (right * graph_size + query_size - 1) // query_size
        )
    else:
        graph_start = (
            piece.graph_end
            - (right * graph_size + query_size - 1) // query_size
        )
        graph_end = piece.graph_end - left * graph_size // query_size
    graph_start = max(0, graph_start)
    if graph_end <= graph_start:
        return query_start, query_end
    return graph_start, graph_end


def representative_mapping_piece(
    pieces: Sequence[SegmentMappingPiece],
    query_start: int,
    query_end: int,
    run: RedundancyRun,
) -> SegmentMappingPiece:
    """Choose one stable mapping instead of spanning incompatible mappings."""
    supporters = [
        piece for piece in pieces
        if interval_overlap(
            piece.query_start, piece.query_end, query_start, query_end,
        ) > 0
    ]
    pool = supporters or list(pieces)
    if not pool:
        return SegmentMappingPiece(
            query_start, query_end, "+", query_start, query_end,
        )

    def priority(piece: SegmentMappingPiece) -> Tuple[int, int, int, int, int, int]:
        distance = max(
            0, piece.query_start - query_end,
            query_start - piece.query_end,
        )
        run_overlap = interval_overlap(
            piece.query_start, piece.query_end, run.start, run.end,
        )
        query_size = piece.query_end - piece.query_start
        graph_size = piece.graph_end - piece.graph_start
        return (
            -distance,
            run_overlap,
            query_size,
            -abs(graph_size - query_size),
            int(piece.strand == "+"),
            -piece.graph_start,
        )

    return max(pool, key=priority)


def mappings_are_continuous(
    left: FinalSegment, right: FinalSegment,
) -> bool:
    if (
        left.strand != right.strand
        or left.query_end != right.query_start
        or left.component_count != right.component_count
    ):
        return False
    if left.strand == "+":
        return abs(left.graph_end - right.graph_start) <= 1
    return abs(left.graph_start - right.graph_end) <= 1


def coalesce_compatible_final_segments(
    segments: Sequence[FinalSegment],
) -> List[FinalSegment]:
    output: List[FinalSegment] = []
    for segment in sorted(segments, key=lambda value: (
        value.query_start, value.query_end, value.strand,
        value.graph_start, value.graph_end,
    )):
        if output and mappings_are_continuous(output[-1], segment):
            output[-1].query_end = segment.query_end
            output[-1].graph_start = min(
                output[-1].graph_start, segment.graph_start,
            )
            output[-1].graph_end = max(
                output[-1].graph_end, segment.graph_end,
            )
            output[-1].mapping_pieces.extend(segment.mapping_pieces)
        else:
            output.append(segment)
    return output


def merge_or_drop_short_final_segments(
    segments: Sequence[FinalSegment], minimum_size: int,
) -> List[FinalSegment]:
    """Absorb short atomic pieces when possible; otherwise enforce the floor."""
    merged = coalesce_compatible_final_segments(segments)
    if minimum_size == 0:
        return merged
    while len(merged) > 1:
        short_index = next((
            index for index, segment in enumerate(merged)
            if segment.query_end - segment.query_start < minimum_size
        ), None)
        if short_index is None:
            break
        neighbors = []
        if short_index > 0:
            neighbors.append(short_index - 1)
        if short_index + 1 < len(merged):
            neighbors.append(short_index + 1)
        short = merged[short_index]

        def priority(index: int) -> Tuple[int, int, int, int]:
            neighbor = merged[index]
            return (
                int(mappings_are_continuous(
                    neighbor, short,
                ) or mappings_are_continuous(short, neighbor)),
                int(neighbor.strand == short.strand),
                int(
                    neighbor.query_end - neighbor.query_start
                    >= minimum_size
                ),
                neighbor.query_end - neighbor.query_start,
            )

        neighbor_index = max(neighbors, key=priority)
        neighbor = merged[neighbor_index]
        neighbor.query_start = min(neighbor.query_start, short.query_start)
        neighbor.query_end = max(neighbor.query_end, short.query_end)
        # A coordinate block may contain discontinuous graph mappings. Keep
        # their bounding graph span instead of refusing the merge.
        neighbor.graph_start = min(neighbor.graph_start, short.graph_start)
        neighbor.graph_end = max(neighbor.graph_end, short.graph_end)
        neighbor.bad_size += short.bad_size
        neighbor.mapping_pieces.extend(short.mapping_pieces)
        del merged[short_index]
        merged.sort(key=lambda value: (value.query_start, value.query_end))
        merged = coalesce_compatible_final_segments(merged)
    return [
        segment for segment in merged
        if segment.query_end - segment.query_start >= minimum_size
    ]


def project_redundancy_runs(
    runs: Sequence[RedundancyRun],
    segments: Sequence[FinalSegment],
    hotspot_index: int,
    graph_name: str,
    query_contig: str,
    minimum_size: int,
) -> List[FinalSegment]:
    """Project nonoverlapping runs without splitting joined graph mappings."""
    output = []
    for run in runs:
        supporters = [
            segment for segment in segments
            if interval_overlap(
                segment.query_start, segment.query_end, run.start, run.end,
            ) > 0
        ]
        if not supporters:
            continue
        # The earliest supporting block supplies the strand, matching the
        # starting-position dominance rule used by neighbor joining.
        representative = min(supporters, key=lambda segment: (
            segment.query_start,
            -interval_overlap(
                segment.query_start, segment.query_end, run.start, run.end,
            ),
            segment.query_end,
        ))
        aggregate_piece = SegmentMappingPiece(
            representative.query_start,
            representative.query_end,
            representative.strand,
            representative.graph_start,
            representative.graph_end,
        )
        graph_start, graph_end = project_mapping_piece(
            aggregate_piece, run.start, run.end,
        )
        if graph_end <= graph_start:
            continue
        output.append(FinalSegment(
            hotspot_index, graph_name, query_contig,
            run.start, run.end, representative.strand,
            graph_start, graph_end,
            component_count=max(1, run.redundancy),
            bad_size=representative.bad_size,
            mapping_pieces=list(representative.mapping_pieces),
        ))
    return merge_or_drop_short_final_segments(output, minimum_size)


def partition_neighbor_joined_segments(
    segments: Sequence[FinalSegment],
) -> List[FinalSegment]:
    """Resolve residual overlaps without erasing neighbor-join boundaries.

    Every joined segment start/end remains an explicit cut. For a run covered
    by multiple unmerged segments, the earliest starting segment supplies the
    coordinate mapping and strand; the coverage depth remains available as
    ``component_count``. Adjacent runs coalesce only when they retain the same
    representative mapping and depth.
    """
    if not segments:
        return []
    boundaries = sorted({
        coordinate
        for segment in segments
        for coordinate in (segment.query_start, segment.query_end)
    })
    output: List[FinalSegment] = []
    for query_start, query_end in zip(boundaries, boundaries[1:]):
        if query_end <= query_start:
            continue
        supporters = [
            segment for segment in segments
            if segment.query_start < query_end
            and segment.query_end > query_start
        ]
        if not supporters:
            continue
        representative = min(supporters, key=lambda segment: (
            segment.query_start,
            -interval_overlap(
                segment.query_start, segment.query_end,
                query_start, query_end,
            ),
            segment.query_end,
        ))
        aggregate_piece = SegmentMappingPiece(
            representative.query_start,
            representative.query_end,
            representative.strand,
            representative.graph_start,
            representative.graph_end,
        )
        graph_start, graph_end = project_mapping_piece(
            aggregate_piece, query_start, query_end,
        )
        if graph_end <= graph_start:
            continue
        output.append(FinalSegment(
            representative.hotspot_index,
            representative.graph_name,
            representative.query_contig,
            query_start,
            query_end,
            representative.strand,
            graph_start,
            graph_end,
            component_count=len(supporters),
            bad_size=representative.bad_size,
            mapping_pieces=list(representative.mapping_pieces),
        ))
    return coalesce_compatible_final_segments(output)


def stabilize_segment_redundancy(
    segments: Sequence[FinalSegment], minimum_size: int,
) -> List[FinalSegment]:
    """Neighbor-join blocks, then retain the existing downstream cleanup."""
    if minimum_size < 0:
        raise ValueError("minimum segment size must be nonnegative")
    grouped: DefaultDict[
        Tuple[int, str, str], List[FinalSegment]
    ] = defaultdict(list)
    for segment in segments:
        grouped[(
            segment.hotspot_index, segment.graph_name,
            segment.query_contig,
        )].append(segment)

    output: List[FinalSegment] = []
    for group_key in sorted(grouped):
        source = neighbor_join_segments(grouped[group_key])
        partitioned = partition_neighbor_joined_segments(source)
        output.extend(merge_or_drop_short_final_segments(
            partitioned, minimum_size,
        ))
    output.sort(key=lambda segment: (
        segment.hotspot_index, segment.graph_name, segment.query_contig,
        segment.query_start, segment.query_end, segment.strand,
        segment.graph_start, segment.graph_end,
    ))
    return output


def build_final_segments(
    rows: Sequence[OutputRow],
    graph_fastas: Sequence[Optional[str]],
    reference_haplotypes: Iterable[str],
    scutoff: int,
    dcutoff: int,
    jobs: int = 1,
    minimum_segment_size: int = DEFAULT_MINIMUM_SEGMENT_SIZE,
    valid_minimum_alignment: int = DEFAULT_VALID_ALIGNMENT_SIZE,
    valid_merge_gap: int = DEFAULT_VALID_REGION_MERGE_GAP,
    only_reference: bool = False,
) -> List[FinalSegment]:
    """Project and merge rows, processing independent graphs in processes."""
    if jobs < 1:
        raise ValueError("jobs must be positive")
    if minimum_segment_size < 0:
        raise ValueError("minimum segment size must be nonnegative")
    if valid_minimum_alignment < 0 or valid_merge_gap < 0:
        raise ValueError("valid-region alignment size/gap must be nonnegative")
    grouped: DefaultDict[str, List[OutputRow]] = defaultdict(list)
    for row in rows:
        graph_fasta = graph_fastas[row.hotspot_index - 1]
        if graph_fasta is None:
            continue
        grouped[graph_fasta].append(row)

    graph_groups = sorted(grouped.items())
    if not graph_groups:
        return []
    reference_tuple = tuple(reference_haplotypes)
    worker_count = min(jobs, len(graph_groups))
    total_rows = sum(
        len(group_rows) for _path, group_rows in graph_groups
    )
    LOG.info(
        "Projecting %d alignment rows across %d local graphs with %d "
        "worker process%s",
        total_rows,
        len(graph_groups), worker_count, "" if worker_count == 1 else "es",
    )
    if worker_count == 1:
        return _project_graph_batch((
            graph_groups, reference_tuple, scutoff, dcutoff,
            minimum_segment_size, valid_minimum_alignment, valid_merge_gap,
            only_reference,
        ))[2]

    # Use several balanced batches per process. This avoids creating tens of
    # thousands of tiny futures while retaining enough work units to balance
    # graphs with very different numbers of alignment rows.
    batch_count = min(len(graph_groups), worker_count * 4)
    batches: List[List[Tuple[str, List[OutputRow]]]] = [
        [] for _index in range(batch_count)
    ]
    heap = [(0, index) for index in range(batch_count)]
    heapq.heapify(heap)
    for item in sorted(
        graph_groups, key=lambda value: (-len(value[1]), value[0]),
    ):
        load, batch_index = heapq.heappop(heap)
        batches[batch_index].append(item)
        heapq.heappush(heap, (load + len(item[1]), batch_index))

    final: List[FinalSegment] = []
    graphs_done = 0
    rows_done = 0
    report_every = max(1, len(batches) // 20)
    with ProcessPoolExecutor(max_workers=worker_count) as executor:
        futures = [
            executor.submit(
                _project_graph_batch,
                (
                    batch, reference_tuple, scutoff, dcutoff,
                    minimum_segment_size, valid_minimum_alignment,
                    valid_merge_gap,
                    only_reference,
                ),
            )
            for batch in batches if batch
        ]
        for completed, future in enumerate(as_completed(futures), 1):
            graph_count, row_count, segments = future.result()
            graphs_done += graph_count
            rows_done += row_count
            final.extend(segments)
            if completed % report_every == 0 or completed == len(futures):
                LOG.info(
                    "Projection progress: %d/%d graphs; %d/%d rows",
                    graphs_done, len(graph_groups), rows_done,
                    total_rows,
                )
    final.sort(key=lambda segment: (
        segment.hotspot_index, segment.graph_name, segment.query_contig,
        segment.query_start, segment.query_end, segment.strand,
        segment.graph_start, segment.graph_end,
    ))
    return final


def _project_graph_batch(
    task: Tuple[
        Sequence[Tuple[str, Sequence[OutputRow]]],
        Tuple[str, ...], int, int, int, int, int, bool,
    ],
) -> Tuple[int, int, List[FinalSegment]]:
    """Worker entry point for one batch of independent local graphs."""
    (
        graph_groups, reference_haplotypes, scutoff, dcutoff,
        minimum_segment_size, valid_minimum_alignment, valid_merge_gap,
        only_reference,
    ) = task
    final: List[FinalSegment] = []
    row_count = 0
    for graph_fasta, rows in graph_groups:
        graph_name = graph_name_from_fasta(graph_fasta)
        paths = read_graph_path_records(graph_fasta)
        beds = read_original_bed(
            graph_fasta,
            reference_haplotypes if only_reference else None,
        )
        valid_regions = build_graph_valid_regions(
            rows, paths, beds, valid_minimum_alignment, valid_merge_gap,
        )
        projected = [
            interval
            for row in rows
            for interval in project_alignment_row(
                row, graph_name, paths, beds,
            )
        ]
        merged = merge_projected_intervals(
            projected, scutoff, dcutoff,
        )
        stabilized = stabilize_segment_redundancy(
            merged, minimum_segment_size,
        )
        clipped = intersect_segments_with_valid_regions(
            stabilized, valid_regions, minimum_segment_size,
        )
        # Blocks are coordinate partitions, not orientation calls. Keep the
        # required output column for downstream compatibility but normalize it.
        final.extend(
            dataclasses.replace(segment, strand="+")
            for segment in clipped
        )
        row_count += len(rows)
    final.sort(key=lambda segment: (
        segment.hotspot_index, segment.graph_name, segment.query_contig,
        segment.query_start, segment.query_end, segment.strand,
        segment.graph_start, segment.graph_end,
    ))
    validate_nonoverlapping_segments(final)
    return len(graph_groups), row_count, final


def write_final_segments(path: str, segments: Sequence[FinalSegment]) -> None:
    validate_nonoverlapping_segments(segments)
    temporary = path + f".tmp.{os.getpid()}"
    Path(os.path.dirname(os.path.abspath(path)) or ".").mkdir(
        parents=True, exist_ok=True,
    )
    try:
        with open(temporary, "wt") as output:
            output.write(SEGMENT_HEADER)
            for segment in segments:
                output.write(segment.line())
        os.replace(temporary, path)
    finally:
        try:
            os.remove(temporary)
        except FileNotFoundError:
            pass


def parse_args(argv: Optional[Sequence[str]] = None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "-i", "--input", required=True,
        help="existing align_partition_hotspots.py output TSV",
    )
    parser.add_argument(
        "-G", "--graph-folder", required=True,
        help="root containing GRAPH_NAME/GRAPH_NAME.FA and GRAPH_NAME.bed",
    )
    parser.add_argument(
        "-L", "--graph-list", required=True,
        help="the same ordered graph list used for alignment",
    )
    parser.add_argument("-o", "--output", required=True)
    parser.add_argument(
        "--scutoff", type=int, default=DEFAULT_SEGMENT_CUTOFF,
    )
    parser.add_argument(
        "--dcutoff", type=int, default=DEFAULT_SEGMENT_CUTOFF,
    )
    parser.add_argument(
        "--min-segment-size", type=int,
        default=DEFAULT_MINIMUM_SEGMENT_SIZE, metavar="INT",
        help=(
            "demote shorter high-redundancy runs, then recursively absorb "
            "short edge/gap runs into an immediate neighbor, dropping a "
            "short final run that cannot be absorbed "
            "(default: 1000; 0 disables)"
        ),
    )
    parser.add_argument(
        "--valid-min-alignment", type=int,
        default=DEFAULT_VALID_ALIGNMENT_SIZE, metavar="INT",
        help=(
            "minimum matched graph/query span used to extend valid regions; "
            "the comparison is strictly greater than this value "
            "(default: 100)"
        ),
    )
    parser.add_argument(
        "--valid-merge-gap", type=int,
        default=DEFAULT_VALID_REGION_MERGE_GAP, metavar="INT",
        help=(
            "merge graph-derived valid intervals when their positive gap is "
            "strictly smaller than this value (default: 1000)"
        ),
    )
    parser.add_argument(
        "-j", "--jobs", type=int, default=16,
        help="local-graph projection worker processes (default: 16)",
    )
    parser.add_argument(
        "--reference-haplotype", action="append",
        help="reference BED-template haplotype; repeat for aliases (default: CHM13_h1)",
    )
    parser.add_argument(
        "--onlyref", "--only-ref", action="store_true",
        help=(
            "use only alignment rows and source BED intervals belonging to "
            "the named reference haplotype(s); default: use all samples"
        ),
    )
    parser.add_argument("-v", "--verbose", action="store_true")
    args = parser.parse_args(argv)
    if args.scutoff < 0 or args.dcutoff < 1:
        parser.error("--scutoff must be nonnegative and --dcutoff must be positive")
    if args.min_segment_size < 0:
        parser.error("--min-segment-size must be nonnegative")
    if args.valid_min_alignment < 0 or args.valid_merge_gap < 0:
        parser.error(
            "--valid-min-alignment and --valid-merge-gap must be nonnegative"
        )
    if args.jobs < 1:
        parser.error("--jobs must be positive")
    return args


def run(args: argparse.Namespace) -> int:
    graph_fastas = read_graph_list(args.graph_list, args.graph_folder)
    reference_haplotypes = tuple(
        args.reference_haplotype or ["CHM13_h1"]
    )
    input_rows = read_alignment_output(args.input)
    rows = select_alignment_rows(
        input_rows, reference_haplotypes, args.onlyref,
    )
    if args.onlyref and not rows:
        raise ValueError(
            "--onlyref selected no alignment rows for: "
            + ",".join(reference_haplotypes)
        )
    if any(row.hotspot_index > len(graph_fastas) for row in rows):
        bad = next(
            row.hotspot_index for row in rows
            if row.hotspot_index > len(graph_fastas)
        )
        raise IndexError(
            f"alignment hotspot index {bad} exceeds graph-list length "
            f"{len(graph_fastas)}"
        )
    segments = build_final_segments(
        rows, graph_fastas,
        reference_haplotypes,
        args.scutoff, args.dcutoff, args.jobs,
        args.min_segment_size,
        args.valid_min_alignment, args.valid_merge_gap,
        args.onlyref,
    )
    write_final_segments(args.output, segments)
    LOG.info(
        "Projected %d graph-alignment rows into %d final source segments: %s",
        len(rows), len(segments), args.output,
    )
    return 0


def main(argv: Optional[Sequence[str]] = None) -> int:
    args = parse_args(argv)
    logging.basicConfig(
        level=logging.DEBUG if args.verbose else logging.INFO,
        format="[%(asctime)s] %(levelname)s: %(message)s",
    )
    try:
        return run(args)
    except Exception as error:
        LOG.error("%s", error)
        if args.verbose:
            LOG.exception("segment projection failed")
        return 1


if __name__ == "__main__":
    raise SystemExit(main())
