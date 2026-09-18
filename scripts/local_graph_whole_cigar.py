#!/usr/bin/env python3
"""Align one complete hotspot sequence to one local graph.

This is the alignment-only part of the legacy LocalGraphAlign.py workflow.
It deliberately stops at graphcigarlight's five-column whole-query result and
never loads ``cache.json`` or splits the query into graph-valid subregions.
"""

from __future__ import annotations

import argparse
import dataclasses
import fcntl
import os
import re
import shutil
import subprocess
import sys
from pathlib import Path
from typing import Dict, Mapping, Optional, Sequence, Tuple


REPEAT_CUT = 5000
BLAST_IDENTITY = 90
MIN_PATH_ALIGNMENT_SIZE = 20
PATH_TOKEN_RE = re.compile(r"([><])([^><]+)")
CIGAR_TOKEN_RE = re.compile(r"([><])([^:<>]+):([^<>]+)")
CIGAR_PAYLOAD_OPERATION_RE = re.compile(r"(\d+)([HID=XMN])([A-Za-z]*)")
SAM_CIGAR_OPERATION_RE = re.compile(r"(\d+)([MIDNSHP=X])")
BLAST_DATABASE_COMPONENT_SUFFIXES = (
    ".nal", ".ndb", ".nhr", ".nin", ".njs", ".nog", ".nos",
    ".not", ".nsq", ".ntf", ".nto",
)


@dataclasses.dataclass(frozen=True)
class WholeGraphCigar:
    query_name: str
    graph_path: str
    graph_cigar: str
    ref_positions: str
    query_positions: str

    @property
    def aligned(self) -> bool:
        return bool(self.graph_path and self.graph_path != "*")


def _tool(extension_dir: str, name: str) -> str:
    path = os.path.join(extension_dir, name)
    if os.path.isfile(path):
        return path
    installed = shutil.which(name)
    if installed:
        return installed
    raise FileNotFoundError(
        f"required graph-alignment tool is missing beside {extension_dir} "
        f"and on PATH: {name}"
    )


def graph_alignment_tool(extension_dir: str, name: str) -> str:
    """Resolve one graph-alignment helper for callers using the batch API."""
    return _tool(extension_dir, name)


def resolve_extension_dir() -> str:
    """Locate graphcigarlight; executables may also be resolved from PATH."""
    here = Path(__file__).resolve().parent
    candidates = [
        os.environ.get("MINSETREF_GRAPH_TOOLS", ""),
        str(here),
        str(here.parent / "Graph"),
        str(here.parent / "Ctyper2" / "Newbuild" / "Extension"),
    ]
    directory = next(
        (value for value in candidates if value and os.path.isfile(
            os.path.join(value, "graphcigarlight.py")
        )),
        str(here),
    )
    required = ("KmerStrd", "runmakeblastdb", "runblastn", "graphcigarlight.py")
    missing = []
    for name in required:
        try:
            _tool(directory, name)
        except FileNotFoundError:
            missing.append(name)
    if missing:
        raise FileNotFoundError(
            f"required graph-alignment tools must be beside this script in "
            f"{directory}; missing: {', '.join(missing)}"
        )
    return directory


def _run(
    command: Sequence[str], *, stdout=None, timeout: Optional[int] = None,
) -> None:
    subprocess.run(
        list(command), check=True, stdout=stdout, timeout=timeout,
    )


def run_graph_alignment_command(
    command: Sequence[str], *, stdout=None, timeout: Optional[int] = None,
) -> None:
    """Run one checked graph-alignment command for the batch alignment path."""
    _run(command, stdout=stdout, timeout=timeout)


def fasta_sequence_lengths(path: str) -> Dict[str, int]:
    """Read FASTA record names and lengths without retaining sequences."""
    lengths: Dict[str, int] = {}
    name: Optional[str] = None
    length = 0
    with open(path, "rt") as handle:
        for line_number, raw in enumerate(handle, 1):
            if raw.startswith(">"):
                if name is not None:
                    lengths[name] = length
                fields = raw[1:].strip().split()
                if not fields:
                    raise ValueError(f"{path}:{line_number}: empty FASTA header")
                name = fields[0]
                if name in lengths:
                    raise ValueError(f"{path}:{line_number}: duplicate FASTA name {name!r}")
                length = 0
            elif name is not None:
                length += len(raw.strip())
            elif raw.strip():
                raise ValueError(f"{path}:{line_number}: sequence precedes FASTA header")
    if name is not None:
        lengths[name] = length
    if not lengths:
        raise ValueError(f"{path}: FASTA contains no records")
    return lengths


def _sam_cigar_operations(cigar: str, context: str) -> list[Tuple[int, str]]:
    operations = [
        (int(size), operation)
        for size, operation in SAM_CIGAR_OPERATION_RE.findall(cigar)
    ]
    if not operations or "".join(
        f"{size}{operation}" for size, operation in operations
    ) != cigar:
        raise ValueError(f"{context}: invalid SAM CIGAR {cigar!r}")
    return operations


def _transpose_standard_sam_row(
    fields: list[str], graph_lengths: Mapping[str, int],
    query_lengths: Mapping[str, int], context: str,
) -> list[str]:
    """Convert standard query-first SAM to graphcigarlight's old layout.

    The legacy converter expects graph-name, flag, query-name, query-start,
    MAPQ, and a graph-first CIGAR. Standard SAM instead reports the query
    first, stores a graph coordinate in POS, and uses query-first I/D
    semantics. This conversion changes all three, including reverse-strand
    coordinate orientation.
    """
    query_name, graph_name = fields[0], fields[2]
    try:
        flag = int(fields[1])
        graph_start = int(fields[3]) - 1
    except ValueError as error:
        raise ValueError(f"{context}: invalid FLAG or POS in SAM row") from error
    if graph_start < 0:
        raise ValueError(f"{context}: mapped SAM POS must be positive")

    operations = _sam_cigar_operations(fields[5], context)
    leading_clip = 0
    first = 0
    while first < len(operations) and operations[first][1] in ("S", "H"):
        leading_clip += operations[first][0]
        first += 1
    last = len(operations)
    while last > first and operations[last - 1][1] in ("S", "H"):
        last -= 1
    core = operations[first:last]
    if not core:
        raise ValueError(f"{context}: SAM CIGAR contains only clipping")

    query_aligned = sum(
        size for size, operation in core if operation in "MI=X"
    )
    graph_aligned = sum(
        size for size, operation in core if operation in "MDN=X"
    )
    query_length = query_lengths[query_name]
    graph_length = graph_lengths[graph_name]
    graph_end = graph_start + graph_aligned
    if graph_end > graph_length:
        raise ValueError(
            f"{context}: alignment ends at {graph_end}, beyond graph record "
            f"{graph_name!r} length {graph_length}"
        )

    reverse = bool(flag & 16)
    if reverse:
        query_start = query_length - leading_clip - query_aligned
        graph_oriented_start = graph_length - graph_end
        core = list(reversed(core))
    else:
        query_start = leading_clip
        graph_oriented_start = graph_start
    if query_start < 0 or query_start + query_aligned > query_length:
        raise ValueError(
            f"{context}: CIGAR/query coordinates exceed query {query_name!r} "
            f"length {query_length}"
        )

    # graphcigarlight swaps I and D immediately after reading the row. Store
    # their inverse here so its internal CIGAR has standard graph-reference /
    # sample-query meaning. A leading H records the graph offset because the
    # legacy POS column is repurposed for the sample-query offset.
    legacy_operations: list[Tuple[int, str]] = []
    if graph_oriented_start:
        legacy_operations.append((graph_oriented_start, "H"))
    legacy_operations.extend(
        (
            size,
            "D" if operation == "I" else "I" if operation in ("D", "N") else operation,
        )
        for size, operation in core
        if operation != "P"
    )
    converted = list(fields)
    converted[0] = graph_name
    converted[2] = query_name
    converted[3] = str(query_start + 1)
    converted[5] = "".join(
        f"{size}{operation}" for size, operation in legacy_operations
    )
    return converted


def normalize_graphcigarlight_sam(
    sam_path: str,
    graph_lengths: Mapping[str, int],
    query_lengths: Mapping[str, int],
) -> Tuple[int, int, int]:
    """Make a SAM safe for legacy graphcigarlight, atomically and idempotently.

    Returns ``(legacy_rows, converted_standard_rows, skipped_unmapped_rows)``.
    Mixed BLAST/Winnowmap files are supported.
    """
    temporary = sam_path + f".layout.tmp.{os.getpid()}"
    legacy_rows = 0
    converted_rows = 0
    skipped_rows = 0
    try:
        with open(sam_path, "rt") as source, open(temporary, "wt") as output:
            for line_number, raw in enumerate(source, 1):
                if not raw.strip() or raw.startswith("@"):
                    output.write(raw)
                    continue
                fields = raw.rstrip("\r\n").split("\t")
                if len(fields) < 6:
                    fields = raw.split()
                context = f"{sam_path}:{line_number}"
                if len(fields) < 6:
                    raise ValueError(f"{context}: expected at least six SAM columns")
                try:
                    flag = int(fields[1])
                except ValueError as error:
                    raise ValueError(f"{context}: invalid SAM FLAG {fields[1]!r}") from error
                if fields[2] == "*" or fields[5] == "*" or flag & 4:
                    skipped_rows += 1
                    continue

                first_is_graph = fields[0] in graph_lengths
                third_is_graph = fields[2] in graph_lengths
                legacy = first_is_graph and not third_is_graph
                standard = third_is_graph and not first_is_graph
                if legacy and fields[2] not in query_lengths:
                    if len(query_lengths) != 1:
                        raise ValueError(
                            f"{context}: legacy SAM query {fields[2]!r} is "
                            "not present in the query FASTA"
                        )
                    # Bulk BLAST labels queries Query_1, Query_2, ... in its
                    # legacy graph-first SAM.  Demultiplexing has left exactly
                    # one possible query, so restore its real FASTA name.
                    fields[2] = next(iter(query_lengths))
                if standard and fields[0] not in query_lengths:
                    if len(query_lengths) != 1:
                        raise ValueError(
                            f"{context}: standard SAM query {fields[0]!r} is "
                            "not present in the query FASTA"
                        )
                    # Some BLAST configurations replace the sole FASTA name
                    # with Query_1. The file is already demultiplexed here, so
                    # mapping that alias to its only possible query is safe.
                    fields[0] = next(iter(query_lengths))
                if legacy == standard:
                    raise ValueError(
                        f"{context}: cannot uniquely determine SAM layout from "
                        f"names {fields[0]!r} and {fields[2]!r}"
                    )
                if standard:
                    fields = _transpose_standard_sam_row(
                        fields, graph_lengths, query_lengths, context,
                    )
                    converted_rows += 1
                else:
                    legacy_rows += 1
                output.write("\t".join(fields) + "\n")
        os.replace(temporary, sam_path)
    except BaseException:
        try:
            os.remove(temporary)
        except FileNotFoundError:
            pass
        raise
    return legacy_rows, converted_rows, skipped_rows


def blast_database_is_complete(graph_fasta: str) -> bool:
    """Quickly validate the required nonempty v4 BLAST database trio."""
    prefix = graph_fasta + "_db"
    # A v4 nucleotide database is the core trio and has no LMDB ``.ndb``.
    # ``.ndb`` can survive an interrupted v5 build by itself; never accept it
    # as part of the v4 database required by this pipeline.
    try:
        core_exists = all(
            os.path.isfile(prefix + suffix)
            and os.path.getsize(prefix + suffix) > 0
            for suffix in (".nhr", ".nin", ".nsq")
        )
        return core_exists and not os.path.exists(prefix + ".ndb")
    except OSError:
        return False


def _build_graph_blast_database(
    prefix: str, graph_fasta: str, extension_dir: str,
    timeout: Optional[int],
) -> None:
    # Remove files from an incomplete or v5 database before forcing a v4
    # rebuild. Complete v4 databases never reach this function.
    for suffix in BLAST_DATABASE_COMPONENT_SUFFIXES:
        try:
            os.remove(prefix + suffix)
        except FileNotFoundError:
            pass
    _run(
        [
            "bash", _tool(extension_dir, "runmakeblastdb"),
            "-in", graph_fasta, "-out", prefix,
            "-blastdb_version", "4",
        ],
        timeout=timeout,
    )
    if not blast_database_is_complete(graph_fasta):
        raise RuntimeError(
            f"makeblastdb produced no complete database for {graph_fasta}"
        )


def ensure_graph_blast_database(
    graph_fasta: str, extension_dir: str, timeout: Optional[int],
) -> str:
    """Build the shared graph BLAST database once, with a filesystem lock."""
    prefix = graph_fasta + "_db"
    if blast_database_is_complete(graph_fasta):
        return prefix
    lock_path = prefix + ".build.lock"
    Path(lock_path).parent.mkdir(parents=True, exist_ok=True)
    with open(lock_path, "a+") as lock:
        fcntl.flock(lock.fileno(), fcntl.LOCK_EX)
        if not blast_database_is_complete(graph_fasta):
            _build_graph_blast_database(
                prefix, graph_fasta, extension_dir, timeout,
            )
        fcntl.flock(lock.fileno(), fcntl.LOCK_UN)
    return prefix


def parse_whole_graph_cigar(path: str) -> WholeGraphCigar:
    with open(path, "rt") as handle:
        raw = handle.readline().rstrip("\r\n")
    fields = raw.split("\t")
    if len(fields) != 5:
        raise ValueError(f"{path}: expected five graph-CIGAR columns, found {len(fields)}")
    return WholeGraphCigar(*fields)


def query_name_and_length(path: str) -> tuple[str, int]:
    name = ""
    length = 0
    records = 0
    with open(path, "rt") as handle:
        for raw in handle:
            if raw.startswith(">"):
                records += 1
                if records > 1:
                    raise ValueError(f"{path}: expected exactly one FASTA record")
                fields = raw[1:].strip().split()
                if not fields:
                    raise ValueError(f"{path}: empty FASTA header")
                name = fields[0]
            elif records:
                length += len(raw.strip())
    if records != 1:
        raise ValueError(f"{path}: expected exactly one FASTA record")
    return name, length


def query_strand(path: str) -> str:
    """Read the orientation written into the query header by KmerStrd."""
    with open(path, "rt") as handle:
        header = handle.readline().strip()
    if not header.startswith(">"):
        raise ValueError(f"{path}: missing FASTA header after KmerStrd")
    for token in header[1:].split()[1:]:
        if token.endswith(("+", "-")) and ":" in token and "-" in token:
            return token[-1]
    return "+"


def insertion_only(query_name: str, query_length: int) -> WholeGraphCigar:
    return WholeGraphCigar(
        query_name=query_name,
        graph_path="*",
        graph_cigar=f"{query_length}I",
        ref_positions=".",
        query_positions=f"0_{query_length}",
    )


def add_query_payloads(
    result: WholeGraphCigar, query_sequence: str,
) -> WholeGraphCigar:
    """Attach oriented query bases to every ``I`` and ``X`` operation."""
    if result.graph_path == "*":
        return dataclasses.replace(
            result, graph_cigar=f"{len(query_sequence)}I{query_sequence}",
        )
    cigars = CIGAR_TOKEN_RE.findall(result.graph_cigar)
    query_positions = _coordinate_list(
        result.query_positions, f"{result.query_name}:qpositions",
    )
    if len(cigars) != len(query_positions):
        raise ValueError(
            f"{result.query_name}: graph CIGAR/query position component "
            f"counts differ: {len(cigars)}/{len(query_positions)}"
        )
    output = []
    for (marker, path_name, cigar), (query_start, query_end) in zip(
        cigars, query_positions,
    ):
        operations = CIGAR_PAYLOAD_OPERATION_RE.findall(cigar)
        rebuilt = "".join(
            f"{length}{operation}{payload}"
            for length, operation, payload in operations
        )
        if not operations or rebuilt != cigar:
            raise ValueError(
                f"{result.query_name}: invalid graph CIGAR {cigar!r}"
            )
        cursor = query_start
        encoded = []
        for length_text, operation, _old_payload in operations:
            length = int(length_text)
            payload = ""
            if operation in {"I", "X"}:
                payload = query_sequence[cursor:cursor + length]
                if len(payload) != length:
                    raise ValueError(
                        f"{result.query_name}: CIGAR payload exceeds query "
                        f"sequence at {cursor}+{length}"
                    )
            encoded.append(f"{length}{operation}{payload}")
            if operation in {"I", "=", "X", "M"}:
                cursor += length
        if cursor != query_end:
            raise ValueError(
                f"{result.query_name}: payloaded CIGAR consumes "
                f"{cursor - query_start} bases but qpositions spans "
                f"{query_end - query_start}"
            )
        output.append(f"{marker}{path_name}:{''.join(encoded)}")
    return dataclasses.replace(result, graph_cigar="".join(output))


def add_graph_cigar_payloads(
    graph_cigar: str, query_sequence: str, query_name: str,
) -> str:
    """Restore I/X bases for one complete, oriented comparison interval.

    Prepared comparisons have already been sliced/merged into a continuous
    query traversal, so component query positions follow from CIGAR lengths.
    H and D consume only graph sequence. A '<' path does not reverse the
    query: the caller supplies the sequence in the row's query orientation.
    Reuse the same payload writer as whole-hotspot alignment.
    """
    components = CIGAR_TOKEN_RE.findall(graph_cigar)
    if not components or ''.join(
        f'{marker}{name}:{body}' for marker, name, body in components
    ) != graph_cigar:
        raise ValueError(f'{query_name}: invalid named graph CIGAR')
    positions = []
    cursor = 0
    for _marker, _name, body in components:
        start = cursor
        cursor += sum(
            int(length) for length, operation, _payload
            in CIGAR_PAYLOAD_OPERATION_RE.findall(body)
            if operation in {'I', '=', 'X', 'M'}
        )
        positions.append(f'{start}_{cursor}')
    if cursor != len(query_sequence):
        raise ValueError(
            f'{query_name}: graph CIGAR consumes {cursor} query bases, '
            f'but the local sequence has {len(query_sequence)} bases'
        )
    result = add_query_payloads(WholeGraphCigar(
        query_name=query_name,
        graph_path=''.join(f'{marker}{name}' for marker, name, _ in components),
        graph_cigar=graph_cigar,
        ref_positions='.',
        query_positions=';'.join(positions),
    ), query_sequence)
    return result.graph_cigar


def _coordinate_list(text: str, context: str) -> list[tuple[int, int]]:
    if not text or text == ".":
        return []
    output = []
    for value in text.split(";"):
        fields = value.split("_")
        if len(fields) != 2:
            raise ValueError(f"{context}: invalid interval {value!r}")
        try:
            start, end = map(int, fields)
        except ValueError as error:
            raise ValueError(
                f"{context}: invalid interval {value!r}"
            ) from error
        if start < 0 or end < start:
            raise ValueError(f"{context}: invalid interval {value!r}")
        output.append((start, end))
    return output


def _path_alignment_size(cigar: str) -> int:
    return sum(
        int(length)
        for length, operation in re.findall(r"(\d+)([M=X])", cigar)
    )


def _prepend_insertion(cigar: str, size: int) -> str:
    if size <= 0:
        return cigar
    leading_hard_clip = re.match(r"^\d+H", cigar)
    split = leading_hard_clip.end() if leading_hard_clip else 0
    return f"{cigar[:split]}{size}I{cigar[split:]}"


def _append_insertion(cigar: str, size: int) -> str:
    if size <= 0:
        return cigar
    trailing_hard_clip = re.search(r"\d+H$", cigar)
    split = trailing_hard_clip.start() if trailing_hard_clip else len(cigar)
    return f"{cigar[:split]}{size}I{cigar[split:]}"


def complete_query_coverage(
    result: WholeGraphCigar,
    query_length: int,
) -> WholeGraphCigar:
    """Represent every uncovered query interval as an insertion.

    graphcigarlight can report only the aligned interval, while qpositions
    remains relative to the complete query FASTA.  Preserve leading, internal,
    and trailing unaligned bases by attaching them to the nearest retained
    graph path.  This makes graph_cigar query coordinates cover
    ``0..query_length`` and therefore agree with the assembly interval stored
    by align_partition_hotspots.py.
    """
    query_length = int(query_length)
    if query_length < 0:
        raise ValueError("query length must be nonnegative")
    if not result.aligned:
        return insertion_only(result.query_name, query_length)

    paths = PATH_TOKEN_RE.findall(result.graph_path)
    cigars = CIGAR_TOKEN_RE.findall(result.graph_cigar)
    ref_positions = _coordinate_list(
        result.ref_positions, f"{result.query_name}:refpositions",
    )
    query_positions = _coordinate_list(
        result.query_positions, f"{result.query_name}:qpositions",
    )
    if len({len(paths), len(cigars), len(ref_positions), len(query_positions)}) != 1:
        raise ValueError(
            f"{result.query_name}: graph path/CIGAR/reference/query component "
            f"counts differ: {len(paths)}/{len(cigars)}/"
            f"{len(ref_positions)}/{len(query_positions)}"
        )

    components: list[list[object]] = []
    cursor = 0
    for index, (path, cigar, ref_position, query_position) in enumerate(zip(
        paths, cigars, ref_positions, query_positions,
    ), 1):
        marker, path_name = path
        cigar_marker, cigar_name, cigar_text = cigar
        if (marker, path_name) != (cigar_marker, cigar_name):
            raise ValueError(
                f"{result.query_name} component {index}: graph path and "
                "graph CIGAR disagree"
            )
        query_start, query_end = query_position
        if query_start < cursor:
            raise ValueError(
                f"{result.query_name} component {index}: overlapping/out-of-order "
                f"qpositions {query_start}_{query_end} after query offset {cursor}"
            )
        if query_end > query_length:
            raise ValueError(
                f"{result.query_name} component {index}: qpositions "
                f"{query_start}_{query_end} exceeds query length {query_length}"
            )
        query_consumed = sum(
            int(length)
            for length, operation in re.findall(r"(\d+)([HID=XMN])", cigar_text)
            if operation in {"I", "=", "X", "M"}
        )
        if query_consumed != query_end - query_start:
            raise ValueError(
                f"{result.query_name} component {index}: CIGAR consumes "
                f"{query_consumed} query bases but qpositions spans "
                f"{query_end - query_start}"
            )

        gap = query_start - cursor
        if gap:
            if components:
                components[-1][2] = _append_insertion(
                    str(components[-1][2]), gap,
                )
                previous_query = components[-1][4]
                previous_query[1] = query_start
            else:
                cigar_text = _prepend_insertion(cigar_text, gap)
                query_start = 0
        components.append([
            marker, path_name, cigar_text, ref_position,
            [query_start, query_end],
        ])
        cursor = query_end

    if not components:
        return insertion_only(result.query_name, query_length)
    trailing = query_length - cursor
    if trailing:
        components[-1][2] = _append_insertion(
            str(components[-1][2]), trailing,
        )
        components[-1][4][1] = query_length

    return WholeGraphCigar(
        query_name=result.query_name,
        graph_path="".join(
            f"{marker}{path_name}"
            for marker, path_name, _cigar, _ref, _query in components
        ),
        graph_cigar="".join(
            f"{marker}{path_name}:{cigar}"
            for marker, path_name, cigar, _ref, _query in components
        ),
        ref_positions=";".join(
            f"{ref[0]}_{ref[1]}"
            for _marker, _path, _cigar, ref, _query in components
        ),
        query_positions=";".join(
            f"{query[0]}_{query[1]}"
            for _marker, _path, _cigar, _ref, query in components
        ),
    )


def normalize_whole_graph_cigar(
    result: WholeGraphCigar,
    query_name: str,
    query_length: int,
    minimum_alignment_size: int = MIN_PATH_ALIGNMENT_SIZE,
) -> WholeGraphCigar:
    """Remove unusable path components and preserve their query bases as I."""
    if not result.aligned:
        return insertion_only(result.query_name or query_name, query_length)

    paths = PATH_TOKEN_RE.findall(result.graph_path)
    cigars = CIGAR_TOKEN_RE.findall(result.graph_cigar)
    ref_positions = _coordinate_list(
        result.ref_positions, f"{result.query_name}:refpositions",
    )
    query_positions = _coordinate_list(
        result.query_positions, f"{result.query_name}:qpositions",
    )
    counts = {
        len(paths), len(cigars), len(ref_positions), len(query_positions),
    }
    if len(counts) != 1:
        raise ValueError(
            f"{result.query_name}: graph path/CIGAR/reference/query component "
            f"counts differ: {len(paths)}/{len(cigars)}/"
            f"{len(ref_positions)}/{len(query_positions)}"
        )

    retained: list[list[object]] = []
    leading_query_start: Optional[int] = None
    for index, (path, cigar, ref_position, query_position) in enumerate(zip(
        paths, cigars, ref_positions, query_positions,
    ), 1):
        marker, path_name = path
        cigar_marker, cigar_name, cigar_text = cigar
        if (marker, path_name) != (cigar_marker, cigar_name):
            raise ValueError(
                f"{result.query_name} component {index}: graph path and "
                "graph CIGAR disagree"
            )
        query_start, query_end = query_position
        usable = (
            ref_position[1] > ref_position[0]
            and query_end > query_start
            and _path_alignment_size(cigar_text) >= minimum_alignment_size
        )
        if usable:
            if not retained and leading_query_start is not None:
                inserted = max(0, query_start - leading_query_start)
                cigar_text = _prepend_insertion(cigar_text, inserted)
                query_start = min(query_start, leading_query_start)
            retained.append([
                marker, path_name, cigar_text, ref_position,
                [query_start, query_end],
            ])
            continue

        # A dropped component after a retained path becomes a terminal
        # insertion on that previous path. Leading dropped components are
        # accumulated and become a prefix insertion on the first kept path.
        if retained:
            previous = retained[-1]
            previous_query = previous[4]
            inserted = max(0, query_end - previous_query[1])
            previous[2] = _append_insertion(previous[2], inserted)
            previous_query[1] = max(previous_query[1], query_end)
        elif leading_query_start is None:
            leading_query_start = query_start
        else:
            leading_query_start = min(leading_query_start, query_start)

    if not retained:
        return insertion_only(result.query_name or query_name, query_length)

    normalized = WholeGraphCigar(
        query_name=result.query_name or query_name,
        graph_path="".join(
            f"{marker}{path_name}"
            for marker, path_name, _cigar, _ref, _query in retained
        ),
        graph_cigar="".join(
            f"{marker}{path_name}:{cigar}"
            for marker, path_name, cigar, _ref, _query in retained
        ),
        ref_positions=";".join(
            f"{ref[0]}_{ref[1]}"
            for _marker, _path, _cigar, ref, _query in retained
        ),
        query_positions=";".join(
            f"{query[0]}_{query[1]}"
            for _marker, _path, _cigar, _ref, query in retained
        ),
    )
    return complete_query_coverage(normalized, query_length)


def align_whole_graph_cigar(
    query_fasta: str,
    graph_fasta: str,
    output_path: str,
    threads: int = 1,
    timeout: Optional[int] = None,
    winnowmap_only: bool = False,
    add_payload: bool = False,
) -> tuple[WholeGraphCigar, str]:
    """Run the selected aligners and return one unsplit whole-query graph CIGAR."""
    if threads < 1:
        raise ValueError("threads must be positive")
    tools = resolve_extension_dir()
    query_fasta = os.path.abspath(query_fasta)
    graph_fasta = os.path.abspath(graph_fasta)
    output_path = os.path.abspath(output_path)
    Path(output_path).parent.mkdir(parents=True, exist_ok=True)

    requested_name, requested_length = query_name_and_length(query_fasta)
    oriented = query_fasta + ".kmerstrd.tmp"
    _run([
        _tool(tools, "KmerStrd"), "-i", query_fasta,
        "-r", graph_fasta, "-o", oriented,
        "-t", str(threads),
    ], timeout=timeout)
    if not os.path.isfile(oriented) or os.path.getsize(oriented) == 0:
        raise RuntimeError(f"KmerStrd produced no query FASTA: {oriented}")
    os.replace(oriented, query_fasta)
    oriented_name, oriented_length = query_name_and_length(query_fasta)
    if oriented_length != requested_length:
        raise ValueError(
            f"KmerStrd changed query length from {requested_length} to {oriented_length}"
        )
    strand = query_strand(query_fasta)
    sam_path = output_path + ".sam"
    if winnowmap_only:
        winnow_script = _tool(tools, "runwinnowmaplight.sh")
        with open(sam_path, "wb") as alignment_output:
            _run([
                "bash", winnow_script, query_fasta, graph_fasta, str(threads),
                str(BLAST_IDENTITY), "300",
            ], stdout=alignment_output, timeout=timeout)
    else:
        lowercase = 0
        with open(query_fasta, "rt") as handle:
            for raw in handle:
                if not raw.startswith(">"):
                    lowercase += sum(base.islower() for base in raw.strip())

        db_prefix = ensure_graph_blast_database(graph_fasta, tools, timeout)
        repeat_options = []
        if lowercase > REPEAT_CUT:
            repeat_options = ["-dust", "yes", "-lcase_masking"]
        _run([
            "bash", _tool(tools, "runblastn"),
            "-task", "megablast", "-query", query_fasta,
            "-db", db_prefix, "-gapopen", "10", "-gapextend", "2",
            "-word_size", "30", "-perc_identity", str(BLAST_IDENTITY),
            *repeat_options,
            "-evalue", "1e-200", "-outfmt", "17", "-out", sam_path,
            "-num_threads", str(threads), "-max_target_seqs", "100",
        ], timeout=timeout)

        if lowercase > REPEAT_CUT:
            winnow_script = _tool(tools, "runwinnowmaplight.sh")
            with open(sam_path, "ab") as alignment_output:
                _run([
                    "bash", winnow_script, query_fasta, graph_fasta,
                    str(threads), str(BLAST_IDENTITY), "300",
                ], stdout=alignment_output, timeout=timeout)

    normalize_graphcigarlight_sam(
        sam_path,
        fasta_sequence_lengths(graph_fasta),
        {oriented_name or requested_name: oriented_length},
    )
    raw_result = output_path + ".raw.tsv"
    _run([
        sys.executable, _tool(tools, "graphcigarlight.py"),
        "-i", sam_path, "-q", query_fasta, "-r", graph_fasta,
        "-o", raw_result,
        *(["--addinsert"] if add_payload else []),
    ], timeout=timeout)
    result = parse_whole_graph_cigar(raw_result)
    if not result.graph_path or not result.graph_cigar:
        result = insertion_only(oriented_name or requested_name, oriented_length)
    else:
        result = normalize_whole_graph_cigar(
            result, oriented_name or requested_name, oriented_length,
        )
    if add_payload:
        with open(query_fasta, "rt") as oriented_query:
            query_sequence = "".join(
                raw.strip() for raw in oriented_query
                if raw and not raw.startswith(">")
            )
        result = add_query_payloads(result, query_sequence)

    temporary = output_path + f".tmp.{os.getpid()}"
    with open(temporary, "wt") as output:
        output.write("\t".join(dataclasses.astuple(result)) + "\n")
    os.replace(temporary, output_path)
    return result, strand


def parse_args(argv: Optional[Sequence[str]] = None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Align one complete hotspot FASTA record and emit one unsplit graph CIGAR",
    )
    parser.add_argument("-i", "--input", required=True, help="single-record hotspot FASTA")
    parser.add_argument("-g", "--graph", required=True, help="partition graph FASTA")
    parser.add_argument("-o", "--output", required=True, help="five-column output TSV")
    parser.add_argument("-t", "--threads", type=int, default=1)
    parser.add_argument("--timeout", type=int, default=3600)
    parser.add_argument(
        "--payload", action="store_true",
        help="attach oriented query sequence to I and X CIGAR operations",
    )
    return parser.parse_args(argv)


def main(argv: Optional[Sequence[str]] = None) -> int:
    args = parse_args(argv)
    align_whole_graph_cigar(
        args.input, args.graph, args.output, args.threads, args.timeout,
        add_payload=args.payload,
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
