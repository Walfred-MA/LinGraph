#!/usr/bin/env python3
"""Shared parsing and graph/query projection for uniform partition blocks."""

from __future__ import annotations

import dataclasses
import importlib.util
import os
import re
import statistics
from collections import defaultdict
from types import ModuleType
from typing import Dict, Iterable, List, Mapping, Sequence, Tuple

from align_partition_hotspots import OutputRow


PATH_TOKEN_RE = re.compile(r"([><])([^><]+)")
CIGAR_TOKEN_RE = re.compile(r"([><])([^:<>]+):([^<>]+)")
CIGAR_OPERATION_RE = re.compile(r"(\d+)([HID=XMN])([A-Za-z]*)")
PATH_COORDINATE_RE = re.compile(r"_(\d+)_(\d+)$")
INITIAL_BLOCK_HEADER = (
    "hotspot_index", "graph_name", "query_contig", "query_start",
    "query_end", "strand", "graph_start", "graph_end",
    "component_count",
)
ALLGRAPHALIGN_COLUMNS = (
    "gindex", "prefix", "chrom", "strand", "start", "end",
)


@dataclasses.dataclass(frozen=True)
class InitialBlock:
    block_id: int
    hotspot_index: int
    graph_name: str
    query_contig: str
    query_start: int
    query_end: int
    strand: str
    graph_start: int
    graph_end: int
    component_count: int


@dataclasses.dataclass(frozen=True)
class AlignmentComponent:
    marker: str
    path: str
    cigar: str
    ref_start: int
    ref_end: int
    query_start: int
    query_end: int
    path_length: int


@dataclasses.dataclass(frozen=True)
class GraphSegment:
    path: str
    start: int
    end: int
    orientation: str
    path_length: int
    query_start: int
    query_end: int


@dataclasses.dataclass(frozen=True)
class QueryProjection:
    query_start: int
    query_end: int
    graph_start: int
    graph_end: int


def _split_fields(raw: str) -> List[str]:
    fields = raw.rstrip("\r\n").split("\t")
    if len(fields) == 1:
        fields = raw.split()
    return fields


def read_initial_blocks(path: str) -> List[InitialBlock]:
    """Read summarize_partition_hotspot_segments.py's nine-column output."""
    output: List[InitialBlock] = []
    with open(path, "rt") as handle:
        for line_number, raw in enumerate(handle, 1):
            if not raw.strip() or raw.startswith("#"):
                continue
            fields = _split_fields(raw)
            if tuple(fields) == INITIAL_BLOCK_HEADER:
                continue
            if len(fields) != len(INITIAL_BLOCK_HEADER):
                raise ValueError(
                    f"{path}:{line_number}: expected nine initial-block "
                    f"columns, found {len(fields)}"
                )
            try:
                values = [
                    int(fields[index]) for index in (0, 3, 4, 6, 7, 8)
                ]
            except ValueError as error:
                raise ValueError(
                    f"{path}:{line_number}: invalid integer field"
                ) from error
            hotspot, qstart, qend, gstart, gend, count = values
            if (
                hotspot < 1 or qstart < 0 or qend <= qstart
                or gstart < 0 or gend <= gstart
                or count < 1 or fields[5] not in {"+", "-"}
            ):
                raise ValueError(
                    f"{path}:{line_number}: invalid 0-based half-open block"
                )
            output.append(InitialBlock(
                len(output) + 1, hotspot, fields[1], fields[2],
                qstart, qend, fields[5], gstart, gend, count,
            ))
    return output


def parse_position_list(
    text: str, context: str, *, allow_empty: bool = False,
) -> List[Tuple[int, int]]:
    if not text or text == ".":
        return []
    output = []
    for encoded in text.split(";"):
        fields = encoded.split("_")
        if len(fields) != 2:
            raise ValueError(f"{context}: invalid coordinate {encoded!r}")
        try:
            start, end = map(int, fields)
        except ValueError as error:
            raise ValueError(
                f"{context}: invalid coordinate {encoded!r}"
            ) from error
        if start < 0 or end < start or (end == start and not allow_empty):
            raise ValueError(f"{context}: invalid interval {encoded!r}")
        output.append((start, end))
    return output


def path_length_from_name(path: str, cigar: str = "") -> int:
    match = PATH_COORDINATE_RE.search(path)
    if match is not None:
        start, end = map(int, match.groups())
        if end > start:
            return end - start
    if cigar:
        reference_size = sum(
            length for length, operation in cigar_operations(
                cigar, f"graph path {path}",
            )
            if operation in {"H", "D", "=", "X", "M", "N"}
        )
        if reference_size > 0:
            return reference_size
    raise ValueError(
        f"graph path {path!r} does not end in _START_END coordinates and "
        "its CIGAR does not encode a positive reference length"
    )


def cigar_operations(cigar: str, context: str) -> List[Tuple[int, str]]:
    output = []
    cursor = 0
    for match in CIGAR_OPERATION_RE.finditer(cigar):
        if match.start() != cursor:
            raise ValueError(f"{context}: invalid graph CIGAR {cigar!r}")
        length, operation, payload = match.groups()
        size = int(length)
        if size < 1:
            raise ValueError(f"{context}: zero-sized CIGAR operation")
        if payload and operation not in {"I", "X"}:
            raise ValueError(
                f"{context}: unexpected sequence payload on {operation}"
            )
        if payload and len(payload) != size:
            raise ValueError(
                f"{context}: {operation} payload length {len(payload)} "
                f"does not match operation length {size}"
            )
        output.append((size, operation))
        cursor = match.end()
    if not output or cursor != len(cigar):
        raise ValueError(f"{context}: invalid graph CIGAR {cigar!r}")
    return output


def parse_alignment_components(row: OutputRow) -> List[AlignmentComponent]:
    """Parse synchronized path/CIGAR/ref/query component columns."""
    if not row.graph_path or row.graph_path == "*":
        return []
    paths = PATH_TOKEN_RE.findall(row.graph_path)
    cigars = CIGAR_TOKEN_RE.findall(row.graph_cigar)
    ref_positions = parse_position_list(
        row.ref_positions,
        f"hotspot {row.hotspot_index} {row.query_name}:refpositions",
        allow_empty=True,
    )
    query_positions = parse_position_list(
        row.query_positions,
        f"hotspot {row.hotspot_index} {row.query_name}:qpositions",
        allow_empty=True,
    )
    counts = {
        len(paths), len(cigars), len(ref_positions), len(query_positions),
    }
    if len(counts) != 1:
        raise ValueError(
            f"hotspot {row.hotspot_index} {row.query_name}: graph path/"
            "CIGAR/reference/query component counts differ: "
            f"{len(paths)}/{len(cigars)}/{len(ref_positions)}/"
            f"{len(query_positions)}"
        )
    output = []
    for index, (path, cigar, ref_position, query_position) in enumerate(zip(
        paths, cigars, ref_positions, query_positions,
    ), 1):
        marker, path_name = path
        cigar_marker, cigar_name, cigar_text = cigar
        if (marker, path_name) != (cigar_marker, cigar_name):
            raise ValueError(
                f"hotspot {row.hotspot_index} {row.query_name} component "
                f"{index}: path and CIGAR disagree"
            )
        operations = cigar_operations(
            cigar_text,
            f"hotspot {row.hotspot_index} {row.query_name} component {index}",
        )
        query_consumed = sum(
            length for length, operation in operations
            if operation in {"I", "=", "X", "M"}
        )
        if query_consumed != query_position[1] - query_position[0]:
            raise ValueError(
                f"hotspot {row.hotspot_index} {row.query_name} component "
                f"{index}: CIGAR consumes {query_consumed} query bases but "
                f"qpositions spans {query_position[1] - query_position[0]}"
            )
        output.append(AlignmentComponent(
            marker, path_name, cigar_text,
            ref_position[0], ref_position[1],
            query_position[0], query_position[1],
            path_length_from_name(path_name, cigar_text),
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


def query_interval_to_graph_segments(
    row: OutputRow, start: int, end: int,
) -> List[GraphSegment]:
    """Project an absolute query interval onto matched graph-path pieces."""
    clipped_start = max(start, row.start)
    clipped_end = min(end, row.end)
    if clipped_end <= clipped_start:
        return []
    target_start, target_end = absolute_to_oriented_query(
        row, clipped_start, clipped_end,
    )
    output = []
    for component in parse_alignment_components(row):
        if (
            component.query_end <= target_start
            or component.query_start >= target_end
        ):
            continue
        ref_cursor = 0 if component.marker == ">" else component.path_length
        query_cursor = component.query_start
        graph_pieces = []
        query_pieces = []
        for length, operation in cigar_operations(
            component.cigar, f"graph path {component.path}",
        ):
            consumes_graph = operation in {"H", "D", "=", "X", "M", "N"}
            consumes_query = operation in {"I", "=", "X", "M"}
            if operation in {"=", "X", "M"}:
                operation_q_start = query_cursor
                operation_q_end = query_cursor + length
                left = max(target_start, operation_q_start)
                right = min(target_end, operation_q_end)
                if right > left:
                    if component.marker == ">":
                        graph_left = ref_cursor + left - operation_q_start
                        graph_right = ref_cursor + right - operation_q_start
                    else:
                        graph_left = ref_cursor - (right - operation_q_start)
                        graph_right = ref_cursor - (left - operation_q_start)
                    graph_pieces.append((graph_left, graph_right))
                    query_pieces.append((left, right))
            if consumes_graph:
                ref_cursor += length if component.marker == ">" else -length
            if consumes_query:
                query_cursor += length
        if not graph_pieces:
            continue
        local_q_start = min(piece[0] for piece in query_pieces)
        local_q_end = max(piece[1] for piece in query_pieces)
        absolute_start, absolute_end = oriented_to_absolute_query(
            row, local_q_start, local_q_end,
        )
        orientation = component.marker
        if row.strand == "-":
            orientation = ">" if orientation == "<" else "<"
        output.append(GraphSegment(
            component.path,
            max(0, min(piece[0] for piece in graph_pieces)),
            min(component.path_length, max(piece[1] for piece in graph_pieces)),
            orientation, component.path_length,
            absolute_start, absolute_end,
        ))
    return sorted(output, key=lambda segment: (
        segment.query_start, segment.query_end, segment.path,
        segment.start, segment.end,
    ))


def graph_interval_to_query(
    row: OutputRow,
    component: AlignmentComponent,
    graph_start: int,
    graph_end: int,
) -> List[QueryProjection]:
    """Project one canonical graph-path interval through one row component."""
    if graph_start < 0 or graph_end <= graph_start:
        return []
    ref_cursor = 0 if component.marker == ">" else component.path_length
    query_cursor = component.query_start
    output = []
    for length, operation in cigar_operations(
        component.cigar, f"graph path {component.path}",
    ):
        consumes_graph = operation in {"H", "D", "=", "X", "M", "N"}
        consumes_query = operation in {"I", "=", "X", "M"}
        if (
            operation == "I"
            and graph_start <= ref_cursor < graph_end
        ):
            absolute_start, absolute_end = oriented_to_absolute_query(
                row, query_cursor, query_cursor + length,
            )
            output.append(QueryProjection(
                absolute_start, absolute_end, ref_cursor, ref_cursor,
            ))
        if operation in {"=", "X", "M"}:
            if component.marker == ">":
                operation_start, operation_end = ref_cursor, ref_cursor + length
                left = max(graph_start, operation_start)
                right = min(graph_end, operation_end)
                if right > left:
                    local_start = query_cursor + left - operation_start
                    local_end = query_cursor + right - operation_start
                else:
                    local_start = local_end = 0
            else:
                operation_start, operation_end = ref_cursor - length, ref_cursor
                left = max(graph_start, operation_start)
                right = min(graph_end, operation_end)
                if right > left:
                    local_start = query_cursor + operation_end - right
                    local_end = query_cursor + operation_end - left
                else:
                    local_start = local_end = 0
            if local_end > local_start:
                absolute_start, absolute_end = oriented_to_absolute_query(
                    row, local_start, local_end,
                )
                output.append(QueryProjection(
                    absolute_start, absolute_end, left, right,
                ))
        if consumes_graph:
            ref_cursor += length if component.marker == ">" else -length
        if consumes_query:
            query_cursor += length
    return output


def merge_query_intervals(
    intervals: Iterable[Tuple[int, int]], gap: int,
) -> List[Tuple[int, int]]:
    ordered = sorted(set(intervals))
    if not ordered:
        return []
    output = [list(ordered[0])]
    for start, end in ordered[1:]:
        if start - output[-1][1] <= gap:
            output[-1][1] = max(output[-1][1], end)
        else:
            output.append([start, end])
    return [(start, end) for start, end in output if end > start]


def interval_union_size(intervals: Iterable[Tuple[int, int]]) -> int:
    return sum(end - start for start, end in merge_query_intervals(intervals, 0))


def load_graphfixbreaks(path: str) -> ModuleType:
    """Load the requested legacy graphfixbreaks module without sys.path edits."""
    absolute = os.path.abspath(path)
    if not os.path.isfile(absolute):
        raise FileNotFoundError(absolute)
    spec = importlib.util.spec_from_file_location(
        f"minsetref_graphfixbreaks_{abs(hash(absolute))}", absolute,
    )
    if spec is None or spec.loader is None:
        raise ImportError(f"cannot load graphfixbreaks module: {absolute}")
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def uniformize_path_coordinates(
    module: ModuleType, path: str, coordinates: Sequence[int],
) -> Dict[int, int]:
    """Call a graphfixbreaks breakpoint-uniforming function through adapters."""
    values = [int(value) for value in coordinates]
    if not values:
        return {}
    if hasattr(module, "combinebreaks_eachpath"):
        names = [f"reference_{index // 2}" for index in range(len(values))]
        result = module.combinebreaks_eachpath(names, values)
        mapping = result[0] if isinstance(result, tuple) else result
        if isinstance(mapping, Mapping):
            return {int(key): int(value) for key, value in mapping.items()}

    if hasattr(module, "combinebreaks"):
        artificial = {}
        for index in range(0, len(values), 2):
            name = f"reference_{index // 2}"
            endpoints = values[index:index + 2]
            artificial[name] = [[
                [path, coordinate, 0, ">", name]
                for coordinate in endpoints
            ]]
        result = module.combinebreaks(artificial)
        if isinstance(result, Mapping):
            return {
                int(key[1]): int(value)
                for key, value in result.items()
                if isinstance(key, tuple) and len(key) == 2 and key[0] == path
            }
        groups = result[0] if isinstance(result, tuple) and result else result
        if isinstance(groups, Sequence):
            mapping = {}
            for group in groups:
                if not group:
                    continue
                representative = int(statistics.median(group))
                mapping.update({int(value): representative for value in group})
            if mapping:
                return mapping

    if hasattr(module, "coordinate_uniform"):
        representatives = [
            int(value) for value in module.coordinate_uniform(sorted(values))
        ]
        if representatives:
            return {
                value: min(representatives, key=lambda ref: abs(ref - value))
                for value in values
            }
    raise AttributeError(
        "graphfixbreaks module has none of combinebreaks_eachpath, "
        "combinebreaks, or coordinate_uniform"
    )


def allgraphalign_line(
    hotspot_index: int, graph_name: str, contig: str,
    strand: str, start: int, end: int,
) -> str:
    """Return AllGraphAlign.py's non-BED six-column layout."""
    return "\t".join(map(str, (
        hotspot_index, graph_name, contig, strand, start, end,
    ))) + "\n"


def allgraphalign_prefix(graph_name: str) -> str:
    """Match AllGraphAlign.py's target basename ``split('_')[0]`` rule."""
    return graph_name.split("_", 1)[0]
