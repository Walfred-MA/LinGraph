#!/usr/bin/env python3
"""Combine minsetref-light partitions into global- and local-novel BEDs.

For every local block, input FASTA header intervals are converted to source
assembly coordinates. Final reporting anchors are removed from unique BED
intervals, then:

    mapped(partition) = headers(partition) - unique(partition)
    global_mapped      = union(mapped(partition))
    global_novel       = union(headers) - global_mapped
    local_novel        = unique(partition) intersect global_mapped

Non-reference local-novel intervals are realigned only to the selected
reference records of partitions where their source interval was mapped. All
valid placements are saved. The chosen placement maximizes
(-distance_to_that_partition's_reference_unique_intervals, alignment_score).
"""
from __future__ import annotations

import argparse
from bisect import bisect_left
import csv
import dataclasses
import hashlib
import json
import logging
import math
import mmap
import os
import shutil
import sqlite3
import stat
import subprocess
import sys
from collections import defaultdict
from concurrent.futures import (
    FIRST_COMPLETED,
    ProcessPoolExecutor,
    ThreadPoolExecutor,
    wait,
)
from pathlib import Path
from typing import (
    Dict,
    FrozenSet,
    Callable,
    Iterable,
    Iterator,
    List,
    Optional,
    Sequence,
    Set,
    Tuple,
    TypeVar,
)

from minsetref_align import AlignmentHit, collect_alignment_hits, ensure_executable
from minsetref_core import IndexedFasta, mp_context, wrap_fasta
from minsetref_light import (
    LIGHT_WINNOW_PARAMS,
    MASKED_BASE_WEIGHT,
    MINIMAP_PARAMS,
    expand_reference_selectors,
)
from minsetref_segments import (
    complement_intervals,
    merge_intervals,
    novelty_score,
)


LOG = logging.getLogger("summarize_partition_novelty")
Interval = Tuple[int, int]
InputValue = TypeVar("InputValue")
OutputValue = TypeVar("OutputValue")


@dataclasses.dataclass(frozen=True)
class HeaderRecord:
    partition: str
    record_id: str
    haplotype: str
    source_contig: str
    source_start: int
    source_end: int
    record_length: int
    record_order: int


@dataclasses.dataclass(frozen=True)
class PartitionFiles:
    partition: str
    fasta: str
    unique_bed: str
    reference_bed: str
    fingerprints: Tuple[Tuple[str, int, int], ...] = dataclasses.field(
        default_factory=tuple, compare=False, repr=False,
    )


@dataclasses.dataclass
class PartitionCatalog:
    files: PartitionFiles
    headers: List[Tuple]
    unique_regions: List[Tuple]
    mapped_regions: List[Tuple]


def header_path_for_fasta(fasta: str) -> str:
    """Prefer canonical ``_samples.header`` while accepting interim names."""
    stem = os.path.splitext(fasta)[0]
    canonical = (
        stem + ".header"
        if stem.endswith("_samples")
        else stem + "_samples.header"
    )
    legacy = stem + ".header"
    if os.path.isfile(canonical):
        return canonical
    if os.path.isfile(legacy):
        return legacy
    return canonical


def partition_header_path(partition: PartitionFiles) -> str:
    return header_path_for_fasta(partition.fasta)


@dataclasses.dataclass(frozen=True)
class UniqueIntervalWork:
    record_id: str
    haplotype: str
    source_contig: str
    source_start: int
    source_end: int
    local_start: int
    local_end: int
    region_class: str


@dataclasses.dataclass(frozen=True)
class CoverageSegment:
    start: int
    end: int
    partitions: FrozenSet[str]


@dataclasses.dataclass
class LocalCandidateResult:
    candidates: List[Tuple]
    mappings: List[Tuple]


@dataclasses.dataclass(frozen=True)
class QueryCore:
    candidate_id: int
    core_start: int
    core_end: int
    sequence: str


@dataclasses.dataclass(frozen=True)
class Placement:
    candidate_id: int
    mapped_partition: str
    ref_haplotype: str
    ref_contig: str
    ref_start: int
    ref_end: int
    strand: str
    distance: int
    alignment_score: float
    identity: float
    aligned_bases: int
    source: str
    target_record: str
    target_blocks: str


SCHEMA = """
PRAGMA journal_mode=WAL;
PRAGMA synchronous=NORMAL;
CREATE TABLE IF NOT EXISTS metadata (
    key TEXT PRIMARY KEY,
    value TEXT NOT NULL
);
CREATE TABLE IF NOT EXISTS partitions (
    partition TEXT PRIMARY KEY,
    fasta TEXT NOT NULL,
    unique_bed TEXT NOT NULL,
    reference_bed TEXT NOT NULL
);
CREATE TABLE IF NOT EXISTS headers (
    partition TEXT NOT NULL,
    record_id TEXT NOT NULL,
    haplotype TEXT NOT NULL,
    source_contig TEXT NOT NULL,
    source_start INTEGER NOT NULL,
    source_end INTEGER NOT NULL,
    record_length INTEGER NOT NULL,
    record_order INTEGER NOT NULL,
    PRIMARY KEY (partition, record_id)
);
CREATE TABLE IF NOT EXISTS unique_regions (
    unique_id INTEGER PRIMARY KEY AUTOINCREMENT,
    partition TEXT NOT NULL,
    record_id TEXT NOT NULL,
    haplotype TEXT NOT NULL,
    source_contig TEXT NOT NULL,
    source_start INTEGER NOT NULL,
    source_end INTEGER NOT NULL,
    local_start INTEGER NOT NULL,
    local_end INTEGER NOT NULL,
    raw_local_start INTEGER NOT NULL,
    raw_local_end INTEGER NOT NULL,
    region_name TEXT NOT NULL,
    region_class TEXT NOT NULL,
    novelty_score REAL NOT NULL
);
CREATE TABLE IF NOT EXISTS mapped_regions (
    partition TEXT NOT NULL,
    record_id TEXT NOT NULL,
    haplotype TEXT NOT NULL,
    source_contig TEXT NOT NULL,
    source_start INTEGER NOT NULL,
    source_end INTEGER NOT NULL
);
CREATE TABLE IF NOT EXISTS partition_mapped_union (
    partition TEXT NOT NULL,
    haplotype TEXT NOT NULL,
    source_contig TEXT NOT NULL,
    source_start INTEGER NOT NULL,
    source_end INTEGER NOT NULL
);
CREATE TABLE IF NOT EXISTS header_union (
    haplotype TEXT NOT NULL,
    source_contig TEXT NOT NULL,
    source_start INTEGER NOT NULL,
    source_end INTEGER NOT NULL
);
CREATE TABLE IF NOT EXISTS global_mapped (
    haplotype TEXT NOT NULL,
    source_contig TEXT NOT NULL,
    source_start INTEGER NOT NULL,
    source_end INTEGER NOT NULL
);
CREATE TABLE IF NOT EXISTS global_novel (
    global_id INTEGER PRIMARY KEY AUTOINCREMENT,
    haplotype TEXT NOT NULL,
    source_contig TEXT NOT NULL,
    source_start INTEGER NOT NULL,
    source_end INTEGER NOT NULL
);
CREATE TABLE IF NOT EXISTS local_candidates (
    candidate_id INTEGER PRIMARY KEY AUTOINCREMENT,
    novel_partition TEXT NOT NULL,
    record_id TEXT NOT NULL,
    haplotype TEXT NOT NULL,
    source_contig TEXT NOT NULL,
    source_start INTEGER NOT NULL,
    source_end INTEGER NOT NULL,
    local_start INTEGER NOT NULL,
    local_end INTEGER NOT NULL,
    region_class TEXT NOT NULL,
    UNIQUE(
        novel_partition, haplotype, source_contig, source_start, source_end
    )
);
CREATE TABLE IF NOT EXISTS candidate_mapped_partitions (
    candidate_id INTEGER NOT NULL,
    mapped_partition TEXT NOT NULL,
    PRIMARY KEY(candidate_id, mapped_partition)
);
CREATE TABLE IF NOT EXISTS all_placements (
    candidate_id INTEGER NOT NULL,
    mapped_partition TEXT NOT NULL,
    ref_haplotype TEXT NOT NULL,
    ref_contig TEXT NOT NULL,
    ref_start INTEGER NOT NULL,
    ref_end INTEGER NOT NULL,
    strand TEXT NOT NULL,
    distance INTEGER NOT NULL,
    alignment_score REAL NOT NULL,
    identity REAL NOT NULL,
    aligned_bases INTEGER NOT NULL,
    source TEXT NOT NULL,
    target_record TEXT NOT NULL,
    target_blocks TEXT NOT NULL
);
CREATE TABLE IF NOT EXISTS best_placements (
    candidate_id INTEGER PRIMARY KEY,
    mapped_partition TEXT NOT NULL,
    ref_haplotype TEXT NOT NULL,
    ref_contig TEXT NOT NULL,
    ref_start INTEGER NOT NULL,
    ref_end INTEGER NOT NULL,
    strand TEXT NOT NULL,
    distance INTEGER NOT NULL,
    alignment_score REAL NOT NULL,
    identity REAL NOT NULL,
    aligned_bases INTEGER NOT NULL,
    source TEXT NOT NULL,
    target_record TEXT NOT NULL,
    target_blocks TEXT NOT NULL
);
"""

SECONDARY_INDEXES = {
    "headers_source": (
        "CREATE INDEX IF NOT EXISTS headers_source "
        "ON headers(haplotype, source_contig, source_start, source_end)"
    ),
    "unique_source": (
        "CREATE INDEX IF NOT EXISTS unique_source "
        "ON unique_regions(haplotype, source_contig, source_start, source_end)"
    ),
    "unique_record": (
        "CREATE INDEX IF NOT EXISTS unique_record "
        "ON unique_regions(partition, record_id, local_start, local_end)"
    ),
    "mapped_source": (
        "CREATE INDEX IF NOT EXISTS mapped_source "
        "ON mapped_regions(haplotype, source_contig, source_start, source_end)"
    ),
    "partition_mapped_source": (
        "CREATE INDEX IF NOT EXISTS partition_mapped_source "
        "ON partition_mapped_union("
        "haplotype, source_contig, source_start, source_end, partition)"
    ),
    "global_mapped_source": (
        "CREATE INDEX IF NOT EXISTS global_mapped_source "
        "ON global_mapped(haplotype, source_contig, source_start, source_end)"
    ),
    "candidate_partition": (
        "CREATE INDEX IF NOT EXISTS candidate_partition "
        "ON candidate_mapped_partitions(mapped_partition, candidate_id)"
    ),
    "placement_candidate": (
        "CREATE INDEX IF NOT EXISTS placement_candidate "
        "ON all_placements(candidate_id)"
    ),
}


def connect_database(
    path: str,
    create_indexes: bool = True,
) -> sqlite3.Connection:
    connection = sqlite3.connect(path, timeout=120)
    connection.row_factory = sqlite3.Row
    connection.executescript(SCHEMA)
    if create_indexes:
        create_secondary_indexes(connection)
    return connection


def drop_secondary_indexes(connection: sqlite3.Connection) -> None:
    for name in SECONDARY_INDEXES:
        connection.execute(f"DROP INDEX IF EXISTS {name}")


def create_secondary_indexes(connection: sqlite3.Connection) -> None:
    for statement in SECONDARY_INDEXES.values():
        connection.execute(statement)
    connection.commit()


def reset_catalog(connection: sqlite3.Connection) -> None:
    # Bulk deletion and reload are substantially faster without maintaining
    # secondary B-trees row by row. Primary and UNIQUE indexes stay active.
    drop_secondary_indexes(connection)
    for table in (
        "metadata", "partitions", "headers", "unique_regions",
        "mapped_regions", "partition_mapped_union", "header_union",
        "global_mapped", "global_novel", "local_candidates",
        "candidate_mapped_partitions", "all_placements", "best_placements",
    ):
        connection.execute(f"DELETE FROM {table}")
    connection.commit()


def threaded_map_unordered(
    function: Callable[[InputValue], OutputValue],
    items: Iterable[InputValue],
    jobs: int,
) -> Iterator[OutputValue]:
    """Bounded unordered thread map, avoiding 60K queued Future objects."""
    if jobs == 1:
        yield from map(function, items)
        return
    iterator = iter(items)
    with ThreadPoolExecutor(max_workers=jobs) as executor:
        pending = set()
        for _index in range(jobs * 4):
            try:
                item = next(iterator)
            except StopIteration:
                break
            pending.add(executor.submit(function, item))
        while pending:
            finished, pending = wait(pending, return_when=FIRST_COMPLETED)
            for future in finished:
                yield future.result()
                try:
                    item = next(iterator)
                except StopIteration:
                    continue
                pending.add(executor.submit(function, item))


def discover_partitions(root: str, jobs: int = 1) -> List[PartitionFiles]:
    """Discover and fingerprint partition inputs with one stat call per file."""
    root = os.path.abspath(root)
    with os.scandir(root) as entries:
        directories = sorted(
            (
                (entry.name, os.path.abspath(entry.path))
                for entry in entries if entry.is_dir()
            ),
            key=lambda item: item[0],
        )

    def inspect(item: Tuple[str, str]) -> Optional[PartitionFiles]:
        name, directory = item
        paths = (
            os.path.join(directory, f"{name}.fasta"),
            os.path.join(directory, f"{name}.unique.bed"),
            os.path.join(directory, f"{name}.bed"),
        )
        fingerprints: List[Optional[Tuple[str, int, int]]] = []
        for path in paths:
            try:
                file_stat = os.stat(path)
            except FileNotFoundError:
                fingerprints.append(None)
                continue
            if not stat.S_ISREG(file_stat.st_mode):
                fingerprints.append(None)
                continue
            fingerprints.append((path, file_stat.st_size, file_stat.st_mtime_ns))
        present = [fingerprint is not None for fingerprint in fingerprints]
        if not any(present):
            return None
        if not all(present):
            missing = [
                path for path, fingerprint in zip(paths, fingerprints)
                if fingerprint is None
            ]
            raise FileNotFoundError(
                f"partition {name} is incomplete: {', '.join(missing)}"
            )
        return PartitionFiles(
            name, paths[0], paths[1], paths[2],
            tuple(fingerprint for fingerprint in fingerprints if fingerprint),
        )

    inspected = threaded_map_unordered(inspect, directories, jobs)
    partitions = sorted(
        (partition for partition in inspected if partition is not None),
        key=lambda partition: partition.partition,
    )
    if not partitions:
        raise ValueError(f"{root}: no complete PARTITION/PARTITION.* inputs")
    return partitions


def catalog_signature(
    partitions: Sequence[PartitionFiles], input_anchor: int,
) -> str:
    digest = hashlib.sha256()
    digest.update(f"catalog-v1\tanchor={input_anchor}\n".encode())
    for partition in partitions:
        digest.update(f"{partition.partition}\n".encode())
        fingerprints = partition.fingerprints
        if not fingerprints:
            fallback = []
            for path in (
                partition.fasta, partition.unique_bed,
                partition.reference_bed,
            ):
                file_stat = os.stat(path)
                fallback.append((
                    path, file_stat.st_size, file_stat.st_mtime_ns,
                ))
            fingerprints = tuple(fallback)
        for path, size, mtime_ns in fingerprints:
            digest.update(
                f"{path}\t{size}\t{mtime_ns}\n".encode()
            )
    return digest.hexdigest()


def fasta_header_descriptions(path: str) -> List[Tuple[int, str]]:
    """Read extracted .header records first, never touching FASTA when present."""
    extracted = header_path_for_fasta(path)
    if os.path.isfile(extracted):
        headers: List[Tuple[int, str]] = []
        with open(extracted, "rt") as handle:
            for record_number, raw in enumerate(handle, 1):
                line = raw.rstrip("\r\n")
                if not line:
                    raise ValueError(
                        f"{extracted}:{record_number}: empty header line"
                    )
                if not line.startswith(">"):
                    raise ValueError(
                        f"{extracted}:{record_number}: expected leading '>'"
                    )
                headers.append((record_number, line[1:]))
        if not headers:
            raise ValueError(f"{extracted}: no FASTA headers")
        return headers

    # Legacy fallback for callers outside the strict coordinate-analysis mode.
    # The coordinate analyzer validates and requires extracted .header inputs.
    fai = path + ".fai"
    if os.path.isfile(fai) and os.stat(fai).st_mtime_ns >= os.stat(path).st_mtime_ns:
        indexed: List[Tuple[str, int]] = []
        try:
            with open(fai, "rt") as handle:
                for raw in handle:
                    fields = raw.rstrip("\n").split("\t")
                    if len(fields) < 3:
                        indexed = []
                        break
                    indexed.append((fields[0], int(fields[2])))
            if indexed:
                headers: List[Tuple[int, str]] = []
                with open(path, "rb") as handle:
                    for record_number, (indexed_name, offset) in enumerate(
                        indexed, 1,
                    ):
                        window_start = max(0, offset - 1024 * 1024)
                        handle.seek(window_start)
                        prefix = handle.read(offset - window_start)
                        if not prefix.endswith(b"\n"):
                            headers = []
                            break
                        end = len(prefix) - 1
                        prior_newline = prefix.rfind(b"\n", 0, end)
                        if prior_newline < 0 and window_start:
                            headers = []
                            break
                        raw_header = prefix[prior_newline + 1:end].rstrip(b"\r")
                        if not raw_header.startswith(b">"):
                            headers = []
                            break
                        header = raw_header[1:].decode("utf-8")
                        if not header or header.split()[0] != indexed_name:
                            headers = []
                            break
                        headers.append((record_number, header))
                if len(headers) == len(indexed):
                    return headers
        except (OSError, UnicodeDecodeError, ValueError):
            pass

    if os.path.getsize(path) == 0:
        return []
    headers = []
    with open(path, "rb") as handle:
        with mmap.mmap(handle.fileno(), 0, access=mmap.ACCESS_READ) as mapped:
            start = 0 if mapped[:1] == b">" else -1
            record_number = 0
            while start >= 0:
                end = mapped.find(b"\n", start)
                if end < 0:
                    end = len(mapped)
                record_number += 1
                headers.append((
                    record_number,
                    mapped[start + 1:end].rstrip(b"\r").decode("utf-8"),
                ))
                marker = mapped.find(b"\n>", end)
                start = marker + 1 if marker >= 0 else -1
    return headers


def parse_header_records(
    partition: PartitionFiles,
    selected_record_ids: Optional[Set[str]] = None,
    selected_haplotypes: Optional[Set[str]] = None,
) -> List[HeaderRecord]:
    """Parse every header or only selected records/haplotypes."""
    records: List[HeaderRecord] = []
    seen: Set[str] = set()
    selected = (
        None if selected_record_ids is None else set(selected_record_ids)
    )
    selected_haps = (
        None if selected_haplotypes is None else set(selected_haplotypes)
    )
    if selected == set():
        return []

    def selected_descriptions():
        extracted = partition_header_path(partition)
        if (
            selected is None
            and selected_haps is None
        ) or not os.path.isfile(extracted):
            yield from fasta_header_descriptions(partition.fasta)
            return
        with open(extracted, "rt") as handle:
            for record_number, raw in enumerate(handle, 1):
                line = raw.rstrip("\r\n")
                if not line or not line.startswith(">"):
                    raise ValueError(
                        f"{extracted}:{record_number}: expected FASTA header"
                    )
                yield record_number, line[1:]

    for record_number, description in selected_descriptions():
        fields = description.strip().split()
        if len(fields) < 2:
            raise ValueError(
                f"{partition.fasta}:record {record_number}: "
                "header lacks source coordinates"
            )
        record_id = fields[0]
        if selected is not None and record_id not in selected:
            continue
        if record_id in seen:
            raise ValueError(
                f"{partition.fasta}:record {record_number}: "
                f"duplicate record {record_id}"
            )
        id_fields = record_id.split("_")
        if len(id_fields) < 4:
            raise ValueError(
                f"{partition.fasta}:record {record_number}: cannot derive "
                f"haplotype from {record_id!r}"
            )
        haplotype = "_".join(id_fields[1:3])
        if selected_haps is not None and haplotype not in selected_haps:
            continue
        coordinate = None
        for token in fields[1:]:
            if ":" not in token or "-" not in token:
                continue
            contig, interval = token.rsplit(":", 1)
            interval = interval.rstrip("+-")
            try:
                start_text, end_text = interval.split("-", 1)
                start, end = int(start_text), int(end_text)
            except ValueError:
                continue
            if start >= 0 and end > start:
                coordinate = (contig, start, end)
                break
        if coordinate is None:
            raise ValueError(
                f"{partition.fasta}:record {record_number}: "
                "no valid contig:start-end"
            )
        contig, start, end = coordinate
        seen.add(record_id)
        records.append(HeaderRecord(
            partition.partition, record_id, haplotype, contig,
            start, end, end - start, record_number - 1,
        ))
        if selected is not None and seen == selected:
            break
    if not records and selected is None and selected_haps is None:
        raise ValueError(f"{partition.fasta}: no FASTA headers")
    return records


def remove_reporting_anchors(
    start: int, end: int, record_length: int, anchor: int,
) -> Interval:
    """Remove reporting anchors; a contig-clipped side contributes no anchor."""
    if start < 0 or end <= start or end > record_length:
        raise ValueError(
            f"invalid unique BED interval {start}-{end} for length {record_length}"
        )
    core_start = start if start == 0 else start + anchor
    core_end = end if end == record_length else end - anchor
    if core_end <= core_start:
        raise ValueError(
            f"anchor removal collapses interval {start}-{end} "
            f"for record length {record_length}"
        )
    return core_start, core_end


def parse_partition_catalog(
    partition: PartitionFiles,
    anchor: int,
) -> PartitionCatalog:
    # Parsing the reference BED here deliberately invokes minsetref_light's
    # strict numbered-block coordinate-name validation.
    expand_reference_selectors([partition.reference_bed])
    records = parse_header_records(partition)
    by_id = {record.record_id: record for record in records}
    header_rows = [
        (
            record.partition, record.record_id, record.haplotype,
            record.source_contig, record.source_start, record.source_end,
            record.record_length, record.record_order,
        )
        for record in records
    ]
    unique_local: Dict[str, List[Interval]] = defaultdict(list)
    unique_rows = []
    with open(partition.unique_bed, "rt") as handle:
        for line_number, raw in enumerate(handle, 1):
            if not raw.strip() or raw.startswith("#"):
                continue
            fields = raw.rstrip("\n").split("\t")
            if len(fields) < 9:
                fields = raw.split()
            if len(fields) < 9:
                raise ValueError(
                    f"{partition.unique_bed}:{line_number}: expected BED9"
                )
            record = by_id.get(fields[0])
            if record is None:
                raise ValueError(
                    f"{partition.unique_bed}:{line_number}: unknown FASTA "
                    f"record {fields[0]!r}"
                )
            raw_start, raw_end = int(fields[1]), int(fields[2])
            local_start, local_end = remove_reporting_anchors(
                raw_start, raw_end, record.record_length, anchor,
            )
            if fields[6] != record.haplotype:
                raise ValueError(
                    f"{partition.unique_bed}:{line_number}: haplotype "
                    f"{fields[6]!r} disagrees with {record.haplotype!r}"
                )
            source_start = record.source_start + local_start
            source_end = record.source_start + local_end
            unique_local[record.record_id].append((local_start, local_end))
            unique_rows.append((
                partition.partition, record.record_id, record.haplotype,
                record.source_contig, source_start, source_end,
                local_start, local_end, raw_start, raw_end,
                fields[3], fields[7], float(fields[8]),
            ))
    mapped_rows = []
    for record in records:
        for start, end in complement_intervals(
            unique_local.get(record.record_id, ()), 0, record.record_length,
        ):
            mapped_rows.append((
                partition.partition, record.record_id, record.haplotype,
                record.source_contig, record.source_start + start,
                record.source_start + end,
            ))
    return PartitionCatalog(
        partition, header_rows, unique_rows, mapped_rows,
    )


def insert_catalog_batch(
    connection: sqlite3.Connection,
    catalogs: Sequence[PartitionCatalog],
) -> None:
    connection.executemany(
        "INSERT INTO partitions VALUES (?, ?, ?, ?)",
        [
            (
                catalog.files.partition, catalog.files.fasta,
                catalog.files.unique_bed, catalog.files.reference_bed,
            )
            for catalog in catalogs
        ],
    )
    connection.executemany(
        "INSERT INTO headers VALUES (?, ?, ?, ?, ?, ?, ?, ?)",
        [row for catalog in catalogs for row in catalog.headers],
    )
    connection.executemany(
        """
        INSERT INTO unique_regions(
            partition, record_id, haplotype, source_contig,
            source_start, source_end, local_start, local_end,
            raw_local_start, raw_local_end, region_name, region_class,
            novelty_score
        ) VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?)
        """,
        [row for catalog in catalogs for row in catalog.unique_regions],
    )
    connection.executemany(
        "INSERT INTO mapped_regions VALUES (?, ?, ?, ?, ?, ?)",
        [row for catalog in catalogs for row in catalog.mapped_regions],
    )


def load_partition(
    connection: sqlite3.Connection,
    partition: PartitionFiles,
    anchor: int,
) -> None:
    """Compatibility wrapper for loading a single parsed partition."""
    insert_catalog_batch(
        connection, [parse_partition_catalog(partition, anchor)],
    )


def load_partitions_parallel(
    connection: sqlite3.Connection,
    partitions: Sequence[PartitionFiles],
    anchor: int,
    jobs: int,
) -> None:
    """Parse partition files concurrently and batch all SQLite writes."""
    batch_size = max(32, min(512, jobs * 32))
    progress_step = max(1, len(partitions) // 20)
    batch: List[PartitionCatalog] = []

    def parse(partition: PartitionFiles) -> PartitionCatalog:
        return parse_partition_catalog(partition, anchor)

    catalogs = threaded_map_unordered(parse, partitions, jobs)
    for index, catalog in enumerate(catalogs, 1):
        batch.append(catalog)
        if len(batch) >= batch_size:
            insert_catalog_batch(connection, batch)
            connection.commit()
            batch.clear()
        if index == len(partitions) or index % progress_step == 0:
            LOG.info("Loaded %d/%d partitions", index, len(partitions))
    if batch:
        insert_catalog_batch(connection, batch)
        connection.commit()


def grouped_interval_rows(
    rows: Iterable[sqlite3.Row],
    key_fields: Sequence[str],
    start_field: str = "source_start",
    end_field: str = "source_end",
) -> Iterator[Tuple[Tuple[str, ...], List[Interval]]]:
    current_key = None
    intervals: List[Interval] = []
    for row in rows:
        key = tuple(str(row[field]) for field in key_fields)
        if current_key is not None and key != current_key:
            yield current_key, merge_intervals(intervals)
            intervals = []
        current_key = key
        intervals.append((int(row[start_field]), int(row[end_field])))
    if current_key is not None:
        yield current_key, merge_intervals(intervals)


def insert_rows_batched(
    connection: sqlite3.Connection,
    statement: str,
    rows: Iterable[Tuple],
    batch_size: int = 50_000,
) -> int:
    """Insert generated rows without retaining an unbounded Python list."""
    batch: List[Tuple] = []
    count = 0
    for row in rows:
        batch.append(row)
        if len(batch) >= batch_size:
            connection.executemany(statement, batch)
            count += len(batch)
            batch.clear()
    if batch:
        connection.executemany(statement, batch)
        count += len(batch)
    return count


def build_interval_unions(connection: sqlite3.Connection) -> None:
    """Build all global unions with sorted linear sweeps."""
    for table in (
        "partition_mapped_union", "header_union", "global_mapped",
        "global_novel", "local_candidates", "candidate_mapped_partitions",
        "all_placements", "best_placements",
    ):
        connection.execute(f"DELETE FROM {table}")

    partition_groups = grouped_interval_rows(
        connection.execute(
            """
            SELECT partition, haplotype, source_contig,
                   source_start, source_end
            FROM mapped_regions
            ORDER BY partition, haplotype, source_contig,
                     source_start, source_end
            """
        ),
        ("partition", "haplotype", "source_contig"),
    )
    insert_rows_batched(
        connection,
        "INSERT INTO partition_mapped_union VALUES (?, ?, ?, ?, ?)",
        (
            (key[0], key[1], key[2], start, end)
            for key, intervals in partition_groups
            for start, end in intervals
        ),
    )

    header_groups = grouped_interval_rows(
        connection.execute(
            """
            SELECT haplotype, source_contig, source_start, source_end
            FROM headers
            ORDER BY haplotype, source_contig, source_start, source_end
            """
        ),
        ("haplotype", "source_contig"),
    )
    insert_rows_batched(
        connection,
        "INSERT INTO header_union VALUES (?, ?, ?, ?)",
        (
            (key[0], key[1], start, end)
            for key, intervals in header_groups
            for start, end in intervals
        ),
    )

    global_groups = grouped_interval_rows(
        connection.execute(
            """
            SELECT haplotype, source_contig, source_start, source_end
            FROM mapped_regions
            ORDER BY haplotype, source_contig, source_start, source_end
            """
        ),
        ("haplotype", "source_contig"),
    )
    insert_rows_batched(
        connection,
        "INSERT INTO global_mapped VALUES (?, ?, ?, ?)",
        (
            (key[0], key[1], start, end)
            for key, intervals in global_groups
            for start, end in intervals
        ),
    )

    # Walk header and mapped unions together by (haplotype, contig). This
    # replaces one overlap query per header interval with a single sorted pass.
    headers = iter(grouped_interval_rows(
        connection.execute(
            """
            SELECT haplotype, source_contig, source_start, source_end
            FROM header_union
            ORDER BY haplotype, source_contig, source_start, source_end
            """
        ),
        ("haplotype", "source_contig"),
    ))
    mapped = iter(grouped_interval_rows(
        connection.execute(
            """
            SELECT haplotype, source_contig, source_start, source_end
            FROM global_mapped
            ORDER BY haplotype, source_contig, source_start, source_end
            """
        ),
        ("haplotype", "source_contig"),
    ))
    mapped_item = next(mapped, None)

    def novel_rows() -> Iterator[Tuple]:
        nonlocal mapped_item
        for header_key, header_intervals in headers:
            while mapped_item is not None and mapped_item[0] < header_key:
                mapped_item = next(mapped, None)
            masks = (
                mapped_item[1]
                if mapped_item is not None and mapped_item[0] == header_key
                else []
            )
            for header_start, header_end in header_intervals:
                for start, end in complement_intervals(
                    masks, header_start, header_end,
                ):
                    yield header_key[0], header_key[1], start, end

    insert_rows_batched(
        connection,
        """
        INSERT INTO global_novel(
            haplotype, source_contig, source_start, source_end
        ) VALUES (?, ?, ?, ?)
        """,
        novel_rows(),
    )
    connection.commit()


def active_partition_segments(
    start: int,
    end: int,
    mapped_rows: Iterable[sqlite3.Row],
) -> List[Tuple[int, int, Tuple[str, ...]]]:
    """Split [start,end) into pieces with a constant non-empty partition set."""
    starts: Dict[int, List[str]] = defaultdict(list)
    ends: Dict[int, List[str]] = defaultdict(list)
    points = {start, end}
    for row in mapped_rows:
        left = max(start, int(row["source_start"]))
        right = min(end, int(row["source_end"]))
        if right <= left:
            continue
        partition = str(row["partition"])
        starts[left].append(partition)
        ends[right].append(partition)
        points.update((left, right))
    active: Set[str] = set()
    output = []
    ordered = sorted(points)
    for index, position in enumerate(ordered[:-1]):
        for partition in ends.get(position, ()):
            active.discard(partition)
        active.update(starts.get(position, ()))
        next_position = ordered[index + 1]
        if next_position > position and active:
            output.append((position, next_position, tuple(sorted(active))))
    return output


def coverage_segments(
    intervals: Sequence[Tuple[str, int, int]],
) -> List[CoverageSegment]:
    """Convert overlapping partition intervals to disjoint active-set spans."""
    events: List[Tuple[int, str, int]] = []
    for partition, start, end in intervals:
        if end <= start:
            continue
        events.append((start, partition, 1))
        events.append((end, partition, -1))
    events.sort(key=lambda event: event[0])
    counts: Dict[str, int] = {}
    segments: List[CoverageSegment] = []
    index = 0
    while index < len(events):
        position = events[index][0]
        while index < len(events) and events[index][0] == position:
            _event_position, partition, delta = events[index]
            count = counts.get(partition, 0) + delta
            if count:
                counts[partition] = count
            else:
                counts.pop(partition, None)
            index += 1
        if index >= len(events):
            break
        next_position = events[index][0]
        if next_position <= position or not counts:
            continue
        active = frozenset(counts)
        if (
            segments
            and segments[-1].end == position
            and segments[-1].partitions == active
        ):
            previous = segments[-1]
            segments[-1] = CoverageSegment(
                previous.start, next_position, active,
            )
        else:
            segments.append(CoverageSegment(
                position, next_position, active,
            ))
    return segments


CoverageIndex = Dict[Tuple[str, str], List[CoverageSegment]]
_LOCAL_COVERAGE_INDEX: CoverageIndex = {}


def build_mapped_coverage_index(
    connection: sqlite3.Connection,
) -> CoverageIndex:
    """Build a sorted disjoint coverage index once for every source contig."""
    index: CoverageIndex = {}
    current_key: Optional[Tuple[str, str]] = None
    intervals: List[Tuple[str, int, int]] = []
    for row in connection.execute(
        """
        SELECT haplotype, source_contig, partition, source_start, source_end
        FROM partition_mapped_union
        ORDER BY haplotype, source_contig, source_start, source_end, partition
        """
    ):
        key = str(row["haplotype"]), str(row["source_contig"])
        if current_key is not None and key != current_key:
            index[current_key] = coverage_segments(intervals)
            intervals = []
        current_key = key
        intervals.append((
            str(row["partition"]), int(row["source_start"]),
            int(row["source_end"]),
        ))
    if current_key is not None:
        index[current_key] = coverage_segments(intervals)
    return index


def local_candidate_workloads(
    connection: sqlite3.Connection,
) -> List[Tuple[str, List[UniqueIntervalWork]]]:
    workloads: List[Tuple[str, List[UniqueIntervalWork]]] = []
    current_partition: Optional[str] = None
    current: List[UniqueIntervalWork] = []
    for row in connection.execute(
        """
        SELECT partition, record_id, haplotype, source_contig,
               source_start, source_end, local_start, local_end, region_class
        FROM unique_regions
        ORDER BY partition, haplotype, source_contig,
                 source_start, source_end, record_id
        """
    ):
        partition = str(row["partition"])
        if current_partition is not None and partition != current_partition:
            workloads.append((current_partition, current))
            current = []
        current_partition = partition
        current.append(UniqueIntervalWork(
            str(row["record_id"]), str(row["haplotype"]),
            str(row["source_contig"]), int(row["source_start"]),
            int(row["source_end"]), int(row["local_start"]),
            int(row["local_end"]), str(row["region_class"]),
        ))
    if current_partition is not None:
        workloads.append((current_partition, current))
    return workloads


def _initialize_local_candidate_worker(index: CoverageIndex) -> None:
    global _LOCAL_COVERAGE_INDEX
    _LOCAL_COVERAGE_INDEX = index


def process_local_candidate_partition(
    workload: Tuple[str, List[UniqueIntervalWork]],
) -> LocalCandidateResult:
    """Intersect one partition's sorted unique intervals with global coverage."""
    partition, unique_intervals = workload
    candidate_rows: Dict[Tuple, Tuple] = {}
    candidate_mappings: Dict[Tuple, Set[str]] = defaultdict(set)
    group_start = 0
    while group_start < len(unique_intervals):
        first = unique_intervals[group_start]
        key = first.haplotype, first.source_contig
        group_end = group_start + 1
        while (
            group_end < len(unique_intervals)
            and (
                unique_intervals[group_end].haplotype,
                unique_intervals[group_end].source_contig,
            ) == key
        ):
            group_end += 1
        segments = _LOCAL_COVERAGE_INDEX.get(key, ())
        segment_cursor = 0
        for unique in unique_intervals[group_start:group_end]:
            while (
                segment_cursor < len(segments)
                and segments[segment_cursor].end <= unique.source_start
            ):
                segment_cursor += 1
            intersections: List[Tuple[int, int, FrozenSet[str]]] = []
            scan = segment_cursor
            while scan < len(segments) and segments[scan].start < unique.source_end:
                segment = segments[scan]
                start = max(unique.source_start, segment.start)
                end = min(unique.source_end, segment.end)
                mapped_partitions = segment.partitions.difference((partition,))
                if end > start and mapped_partitions:
                    active = frozenset(mapped_partitions)
                    if (
                        intersections
                        and intersections[-1][1] == start
                        and intersections[-1][2] == active
                    ):
                        prior_start, _prior_end, prior_active = intersections[-1]
                        intersections[-1] = prior_start, end, prior_active
                    else:
                        intersections.append((start, end, active))
                scan += 1
            for start, end, mapped_partitions in intersections:
                local_start = (
                    unique.local_start + start - unique.source_start
                )
                local_end = local_start + end - start
                natural_key = (
                    partition, unique.haplotype, unique.source_contig,
                    start, end,
                )
                candidate_rows.setdefault(natural_key, (
                    partition, unique.record_id, unique.haplotype,
                    unique.source_contig, start, end, local_start, local_end,
                    unique.region_class,
                ))
                candidate_mappings[natural_key].update(mapped_partitions)
        group_start = group_end
    candidates = [
        candidate_rows[key] for key in sorted(candidate_rows)
    ]
    mappings = [
        (*key, mapped_partition)
        for key in sorted(candidate_mappings)
        for mapped_partition in sorted(candidate_mappings[key])
    ]
    return LocalCandidateResult(candidates, mappings)


def _create_candidate_staging(connection: sqlite3.Connection) -> None:
    connection.executescript(
        """
        DROP TABLE IF EXISTS temp.local_candidate_stage;
        DROP TABLE IF EXISTS temp.candidate_mapping_stage;
        CREATE TEMP TABLE local_candidate_stage (
            novel_partition TEXT NOT NULL,
            record_id TEXT NOT NULL,
            haplotype TEXT NOT NULL,
            source_contig TEXT NOT NULL,
            source_start INTEGER NOT NULL,
            source_end INTEGER NOT NULL,
            local_start INTEGER NOT NULL,
            local_end INTEGER NOT NULL,
            region_class TEXT NOT NULL,
            PRIMARY KEY(
                novel_partition, haplotype, source_contig,
                source_start, source_end
            )
        ) WITHOUT ROWID;
        CREATE TEMP TABLE candidate_mapping_stage (
            novel_partition TEXT NOT NULL,
            haplotype TEXT NOT NULL,
            source_contig TEXT NOT NULL,
            source_start INTEGER NOT NULL,
            source_end INTEGER NOT NULL,
            mapped_partition TEXT NOT NULL,
            PRIMARY KEY(
                novel_partition, haplotype, source_contig,
                source_start, source_end, mapped_partition
            )
        ) WITHOUT ROWID;
        """
    )


def _insert_candidate_result(
    connection: sqlite3.Connection,
    result: LocalCandidateResult,
) -> None:
    connection.executemany(
        "INSERT OR IGNORE INTO local_candidate_stage VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?)",
        result.candidates,
    )
    connection.executemany(
        "INSERT OR IGNORE INTO candidate_mapping_stage VALUES (?, ?, ?, ?, ?, ?)",
        result.mappings,
    )


def _finalize_candidate_staging(connection: sqlite3.Connection) -> None:
    connection.execute("DELETE FROM local_candidates")
    connection.execute("DELETE FROM candidate_mapped_partitions")
    connection.execute(
        """
        INSERT INTO local_candidates(
            novel_partition, record_id, haplotype, source_contig,
            source_start, source_end, local_start, local_end, region_class
        )
        SELECT novel_partition, record_id, haplotype, source_contig,
               source_start, source_end, local_start, local_end, region_class
        FROM local_candidate_stage
        ORDER BY novel_partition, haplotype, source_contig,
                 source_start, source_end
        """
    )
    connection.execute(
        """
        INSERT OR IGNORE INTO candidate_mapped_partitions(
            candidate_id, mapped_partition
        )
        SELECT candidate.candidate_id, mapping.mapped_partition
        FROM candidate_mapping_stage AS mapping
        JOIN local_candidates AS candidate
          ON candidate.novel_partition=mapping.novel_partition
         AND candidate.haplotype=mapping.haplotype
         AND candidate.source_contig=mapping.source_contig
         AND candidate.source_start=mapping.source_start
         AND candidate.source_end=mapping.source_end
        """
    )
    connection.execute("DROP TABLE local_candidate_stage")
    connection.execute("DROP TABLE candidate_mapping_stage")
    connection.commit()


def build_local_candidates(
    connection: sqlite3.Connection,
    jobs: int = 1,
) -> None:
    """Build local candidates by parallel per-partition sorted sweeps."""
    coverage_index = build_mapped_coverage_index(connection)
    workloads = local_candidate_workloads(connection)
    _create_candidate_staging(connection)
    progress_step = max(1, len(workloads) // 20)
    LOG.info(
        "Local-candidate interval sweeps: %d partitions with %d workers",
        len(workloads), jobs,
    )
    if jobs == 1:
        _initialize_local_candidate_worker(coverage_index)
        results: Iterable[LocalCandidateResult] = map(
            process_local_candidate_partition, workloads,
        )
        for index, result in enumerate(results, 1):
            _insert_candidate_result(connection, result)
            if index == len(workloads) or index % progress_step == 0:
                LOG.info(
                    "Processed local candidates for %d/%d partitions",
                    index, len(workloads),
                )
    else:
        chunksize = max(1, min(64, len(workloads) // max(1, jobs * 8)))
        with ProcessPoolExecutor(
            max_workers=jobs,
            mp_context=mp_context(),
            initializer=_initialize_local_candidate_worker,
            initargs=(coverage_index,),
        ) as executor:
            for index, result in enumerate(executor.map(
                process_local_candidate_partition,
                workloads,
                chunksize=chunksize,
            ), 1):
                _insert_candidate_result(connection, result)
                if index % 100 == 0:
                    connection.commit()
                if index == len(workloads) or index % progress_step == 0:
                    LOG.info(
                        "Processed local candidates for %d/%d partitions",
                        index, len(workloads),
                    )
    _finalize_candidate_staging(connection)


def bed_distance(intervals: Sequence[Interval], start: int, end: int) -> int:
    """Distance to sorted, merged intervals in O(log n)."""
    if not intervals:
        return sys.maxsize
    index = bisect_left(intervals, (end,))
    best = sys.maxsize
    if index:
        _left, right = intervals[index - 1]
        if right > start:
            return 0
        best = start - right
    if index < len(intervals):
        left, _right = intervals[index]
        best = min(best, max(0, left - end))
    return best


def clip_hit_to_query_core(
    hit: AlignmentHit,
    core_start: int,
    core_end: int,
) -> Tuple[List[Interval], List[Interval]]:
    query_blocks: List[Interval] = []
    target_blocks: List[Interval] = []
    for q0, q1, t0, t1 in hit.aligned_pairs:
        left = max(q0, core_start)
        right = min(q1, core_end)
        if right <= left:
            continue
        query_blocks.append((left, right))
        if hit.strand == "-":
            target_blocks.append((
                t0 + (q1 - right),
                t0 + (q1 - left),
            ))
        else:
            target_blocks.append((
                t0 + (left - q0),
                t0 + (right - q0),
            ))
    return merge_intervals(query_blocks), merge_intervals(target_blocks)


def placement_sort_key(placement: Placement) -> Tuple:
    return (
        placement.distance,
        -placement.alignment_score,
        -placement.aligned_bases,
        -placement.identity,
        placement.ref_haplotype,
        placement.ref_contig,
        placement.ref_start,
        placement.ref_end,
        placement.mapped_partition,
        placement.target_record,
        placement.source,
    )


def write_tsv(path: str, fieldnames: Sequence[str], rows: Iterable[Sequence]) -> None:
    Path(path).parent.mkdir(parents=True, exist_ok=True)
    temporary = path + ".tmp"
    try:
        with open(temporary, "wt", newline="") as handle:
            writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
            writer.writerow(fieldnames)
            writer.writerows(rows)
        os.replace(temporary, path)
    finally:
        try:
            os.remove(temporary)
        except FileNotFoundError:
            pass


def write_bed_table(
    path: str, fieldnames: Sequence[str], rows: Iterable[Sequence],
) -> None:
    """Write a BED-compatible table with a comment-prefixed column header."""
    header = list(fieldnames)
    header[0] = "#" + header[0]
    write_tsv(path, header, rows)


def alignment_task(
    database_path: str,
    mapped_partition: str,
    work_root: str,
    options: Dict[str, object],
) -> str:
    workdir = os.path.join(work_root, mapped_partition)
    Path(workdir).mkdir(parents=True, exist_ok=True)
    result_path = os.path.join(workdir, "placements.tsv")
    marker = os.path.join(workdir, "_SUCCESS")
    task_signature = hashlib.sha256(json.dumps(
        {
            key: value for key, value in options.items()
            if key not in {"resume"}
        },
        sort_keys=True,
        separators=(",", ":"),
    ).encode()).hexdigest()
    reusable = False
    if options["resume"] and os.path.isfile(marker) and os.path.isfile(result_path):
        try:
            with open(marker) as handle:
                reusable = json.load(handle).get("signature") == task_signature
        except (OSError, ValueError, AttributeError):
            reusable = False
    if reusable:
        return result_path
    if os.path.isdir(workdir):
        shutil.rmtree(workdir)
    Path(workdir).mkdir(parents=True, exist_ok=True)

    connection = connect_database(database_path)
    partition_row = connection.execute(
        "SELECT * FROM partitions WHERE partition=?", (mapped_partition,),
    ).fetchone()
    if partition_row is None:
        raise RuntimeError(f"unknown mapped partition {mapped_partition}")
    candidates = list(connection.execute(
        """
        SELECT c.*, p.fasta
        FROM candidate_mapped_partitions cmp
        JOIN local_candidates c ON c.candidate_id=cmp.candidate_id
        JOIN partitions p ON p.partition=c.novel_partition
        WHERE cmp.mapped_partition=?
        ORDER BY c.candidate_id
        """,
        (mapped_partition,),
    ))
    reference_records = list(connection.execute(
        """
        SELECT DISTINCT h.*
        FROM headers h
        JOIN unique_regions u
          ON u.partition=h.partition AND u.record_id=h.record_id
        WHERE h.partition=? AND u.region_class='reference'
        ORDER BY h.record_order
        """,
        (mapped_partition,),
    ))
    reference_unique: Dict[str, List[Interval]] = defaultdict(list)
    for row in connection.execute(
        """
        SELECT record_id, local_start, local_end
        FROM unique_regions
        WHERE partition=? AND region_class='reference'
        ORDER BY record_id, local_start
        """,
        (mapped_partition,),
    ):
        reference_unique[str(row["record_id"])].append(
            (int(row["local_start"]), int(row["local_end"]))
        )
    reference_unique = {
        record_id: merge_intervals(intervals)
        for record_id, intervals in reference_unique.items()
    }
    connection.close()
    if not candidates or not reference_records:
        write_tsv(result_path, [field.name for field in dataclasses.fields(Placement)], [])
        with open(marker, "wt") as out:
            json.dump({"signature": task_signature}, out, sort_keys=True)
            out.write("\n")
        return result_path

    target_fasta = os.path.join(workdir, "reference.fa")
    target_meta: Dict[str, sqlite3.Row] = {
        str(row["record_id"]): row for row in reference_records
    }
    with IndexedFasta(str(partition_row["fasta"])) as source, open(target_fasta, "wt") as out:
        for row in reference_records:
            record_id = str(row["record_id"])
            out.write(f">{record_id}\n{wrap_fasta(source.sequence(record_id))}\n")

    query_fasta = os.path.join(workdir, "queries.fa")
    query_meta: Dict[str, QueryCore] = {}
    by_fasta: Dict[str, List[sqlite3.Row]] = defaultdict(list)
    for row in candidates:
        by_fasta[str(row["fasta"])].append(row)
    with open(query_fasta, "wt") as out:
        for fasta_path, members in by_fasta.items():
            with IndexedFasta(fasta_path) as source:
                for row in members:
                    record_id = str(row["record_id"])
                    record_length = source.length(record_id)
                    local_start = int(row["local_start"])
                    local_end = int(row["local_end"])
                    extract_start = (
                        local_start if local_start == 0
                        else max(0, local_start - int(options["anchor"]))
                    )
                    extract_end = (
                        local_end if local_end == record_length
                        else min(record_length, local_end + int(options["anchor"]))
                    )
                    sequence = source.fetch(record_id, extract_start, extract_end)
                    query_id = f"localnovel_{int(row['candidate_id']):09d}"
                    core_start = local_start - extract_start
                    core_end = core_start + local_end - local_start
                    query_meta[query_id] = QueryCore(
                        int(row["candidate_id"]), core_start, core_end, sequence,
                    )
                    out.write(f">{query_id}\n{wrap_fasta(sequence)}\n")

    hits: List[AlignmentHit] = []
    minimap_dir = os.path.join(workdir, "minimap2")
    hits.extend(collect_alignment_hits(
        query_fasta, target_fasta, minimap_dir, int(options["cores"]),
        float(options["min_identity"]), int(options["min_score"]), LOG,
        prefix="local_novel", skip_blastn=True,
        broad_aligner="minimap2", minimap_params=MINIMAP_PARAMS,
        masked_weight=MASKED_BASE_WEIGHT,
        parse_min_segment=min(50, int(options["min_score"])),
        reuse_existing=False,
    ))
    sensitive_dir = os.path.join(workdir, "winnowmap_blastn")
    hits.extend(collect_alignment_hits(
        query_fasta, target_fasta, sensitive_dir, int(options["cores"]),
        float(options["min_identity"]), int(options["min_score"]), LOG,
        winnow_params=LIGHT_WINNOW_PARAMS,
        blast_word_size=int(options["blast_word_size"]),
        blast_evalue=str(options["blast_evalue"]),
        blast_max_target_seqs=int(options["blast_max_target_seqs"]),
        prefix="local_novel", skip_blastn=bool(options["skip_blastn"]),
        blast_mode="independent", broad_aligner="winnowmap",
        masked_weight=MASKED_BASE_WEIGHT,
        parse_min_segment=min(50, int(options["min_score"])),
        reuse_existing=False,
    ))

    placements: Dict[Tuple, Placement] = {}
    for hit in hits:
        query = query_meta.get(hit.query_id)
        target = target_meta.get(hit.target_id)
        if query is None or target is None:
            continue
        query_blocks, target_blocks = clip_hit_to_query_core(
            hit, query.core_start, query.core_end,
        )
        if not query_blocks or not target_blocks:
            continue
        aligned_bases = sum(end - start for start, end in query_blocks)
        aligned_score = sum(
            novelty_score(query.sequence[start:end], MASKED_BASE_WEIGHT)
            for start, end in query_blocks
        )
        if (
            aligned_bases < int(options["min_score"])
            or aligned_score < int(options["min_score"])
        ):
            continue
        target_start = min(start for start, _end in target_blocks)
        target_end = max(end for _start, end in target_blocks)
        distance = bed_distance(
            reference_unique.get(hit.target_id, ()), target_start, target_end,
        )
        source_start = int(target["source_start"])
        source_blocks = [
            (source_start + start, source_start + end)
            for start, end in target_blocks
        ]
        placement = Placement(
            query.candidate_id,
            mapped_partition,
            str(target["haplotype"]),
            str(target["source_contig"]),
            min(start for start, _end in source_blocks),
            max(end for _start, end in source_blocks),
            hit.strand,
            distance,
            float(hit.alignment_score),
            float(hit.identity),
            aligned_bases,
            hit.source,
            hit.target_id,
            json.dumps(source_blocks, separators=(",", ":")),
        )
        key = (
            placement.candidate_id, placement.mapped_partition,
            placement.ref_haplotype, placement.ref_contig,
            placement.ref_start, placement.ref_end, placement.strand,
        )
        prior = placements.get(key)
        if prior is None or placement_sort_key(placement) < placement_sort_key(prior):
            placements[key] = placement
    write_tsv(
        result_path,
        [field.name for field in dataclasses.fields(Placement)],
        [dataclasses.astuple(row) for row in sorted(
            placements.values(), key=lambda row: (
                row.candidate_id, placement_sort_key(row),
            ),
        )],
    )
    with open(marker, "wt") as out:
        json.dump({"signature": task_signature}, out, sort_keys=True)
        out.write("\n")
    return result_path


def load_placement_file(connection: sqlite3.Connection, path: str) -> None:
    with open(path, "rt", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        rows = []
        for row in reader:
            rows.append((
                int(row["candidate_id"]), row["mapped_partition"],
                row["ref_haplotype"], row["ref_contig"],
                int(row["ref_start"]), int(row["ref_end"]), row["strand"],
                int(row["distance"]), float(row["alignment_score"]),
                float(row["identity"]), int(row["aligned_bases"]),
                row["source"], row["target_record"], row["target_blocks"],
            ))
    connection.executemany(
        "INSERT INTO all_placements VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?)",
        rows,
    )


def add_reference_source_placements(connection: sqlite3.Connection) -> None:
    """Add already-reference placements using pre-sorted interval indexes."""
    reference_haplotypes = {
        str(row[0]) for row in connection.execute(
            "SELECT DISTINCT haplotype FROM unique_regions WHERE region_class='reference'"
        )
    }
    eligible_by_candidate: Dict[int, List[str]] = defaultdict(list)
    for row in connection.execute(
        """
        SELECT cmp.candidate_id, cmp.mapped_partition
        FROM candidate_mapped_partitions AS cmp
        JOIN local_candidates AS candidate USING(candidate_id)
        WHERE candidate.haplotype IN (
            SELECT DISTINCT haplotype FROM unique_regions
            WHERE region_class='reference'
        )
        ORDER BY cmp.candidate_id, cmp.mapped_partition
        """
    ):
        eligible_by_candidate[int(row["candidate_id"])].append(
            str(row["mapped_partition"])
        )

    reference_intervals: Dict[
        Tuple[str, str, str], List[Interval]
    ] = {}
    groups = grouped_interval_rows(
        connection.execute(
            """
            SELECT partition, haplotype, source_contig,
                   source_start, source_end
            FROM unique_regions
            WHERE region_class='reference'
            ORDER BY partition, haplotype, source_contig,
                     source_start, source_end
            """
        ),
        ("partition", "haplotype", "source_contig"),
    )
    for key, intervals in groups:
        reference_intervals[(key[0], key[1], key[2])] = intervals

    placement_rows: List[Tuple] = []
    for candidate in connection.execute(
        """
        SELECT * FROM local_candidates
        ORDER BY candidate_id
        """
    ):
        if candidate["haplotype"] not in reference_haplotypes:
            continue
        choices = []
        for partition in eligible_by_candidate.get(
            int(candidate["candidate_id"]), (),
        ):
            intervals = reference_intervals.get(
                (
                    partition, str(candidate["haplotype"]),
                    str(candidate["source_contig"]),
                ),
                (),
            )
            choices.append((
                bed_distance(
                    intervals, int(candidate["source_start"]),
                    int(candidate["source_end"]),
                ),
                partition,
            ))
        distance, partition = min(choices) if choices else (sys.maxsize, ".")
        placement_rows.append((
            candidate["candidate_id"], partition, candidate["haplotype"],
            candidate["source_contig"], candidate["source_start"],
            candidate["source_end"], "+", distance, 0, 100,
            int(candidate["source_end"]) - int(candidate["source_start"]),
            "already_reference", candidate["record_id"],
            json.dumps([[
                int(candidate["source_start"]), int(candidate["source_end"]),
            ]], separators=(",", ":")),
        ))
    connection.executemany(
        "INSERT INTO all_placements VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?)",
        placement_rows,
    )


def select_best_placements(connection: sqlite3.Connection) -> None:
    connection.execute("DELETE FROM best_placements")
    current_candidate = None
    current: List[Placement] = []

    def finish() -> None:
        if not current:
            return
        best = min(current, key=placement_sort_key)
        connection.execute(
            "INSERT INTO best_placements VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?)",
            dataclasses.astuple(best),
        )

    for row in connection.execute(
        "SELECT * FROM all_placements ORDER BY candidate_id",
    ):
        candidate_id = int(row["candidate_id"])
        if current_candidate is not None and candidate_id != current_candidate:
            finish()
            current = []
        current_candidate = candidate_id
        current.append(Placement(
            candidate_id, str(row["mapped_partition"]),
            str(row["ref_haplotype"]), str(row["ref_contig"]),
            int(row["ref_start"]), int(row["ref_end"]), str(row["strand"]),
            int(row["distance"]), float(row["alignment_score"]),
            float(row["identity"]), int(row["aligned_bases"]),
            str(row["source"]), str(row["target_record"]),
            str(row["target_blocks"]),
        ))
    finish()
    connection.commit()


def mapped_partitions_by_candidate(
    connection: sqlite3.Connection,
) -> Dict[int, Tuple[str, ...]]:
    memberships: Dict[int, List[str]] = defaultdict(list)
    for row in connection.execute(
        """
        SELECT candidate_id, mapped_partition
        FROM candidate_mapped_partitions
        ORDER BY candidate_id, mapped_partition
        """
    ):
        memberships[int(row["candidate_id"])].append(
            str(row["mapped_partition"])
        )
    return {
        candidate_id: tuple(partitions)
        for candidate_id, partitions in memberships.items()
    }


def export_catalog(connection: sqlite3.Connection, output_dir: str) -> None:
    candidate_memberships = mapped_partitions_by_candidate(connection)
    write_bed_table(
        os.path.join(output_dir, "all_header_regions.bed"),
        (
            "source_contig", "source_start", "source_end", "record_id",
            "score", "strand", "haplotype", "partition", "record_length",
        ),
        (
            (
                row["source_contig"], row["source_start"], row["source_end"],
                row["record_id"], 0, "+", row["haplotype"],
                row["partition"], row["record_length"],
            )
            for row in connection.execute(
                """
                SELECT * FROM headers
                ORDER BY haplotype, source_contig, source_start, partition
                """
            )
        ),
    )
    write_bed_table(
        os.path.join(output_dir, "all_unique_regions.bed"),
        (
            "source_contig", "source_start", "source_end", "region_name",
            "score", "strand", "haplotype", "partition", "record_id",
            "region_class", "local_start", "local_end",
            "raw_local_start", "raw_local_end", "novelty_score",
        ),
        (
            (
                row["source_contig"], row["source_start"], row["source_end"],
                row["region_name"], 0, "+", row["haplotype"],
                row["partition"], row["record_id"], row["region_class"],
                row["local_start"], row["local_end"],
                row["raw_local_start"], row["raw_local_end"],
                row["novelty_score"],
            )
            for row in connection.execute(
                """
                SELECT * FROM unique_regions
                ORDER BY haplotype, source_contig, source_start, partition
                """
            )
        ),
    )
    write_bed_table(
        os.path.join(output_dir, "all_mapped_regions.bed"),
        (
            "source_contig", "source_start", "source_end", "record_id",
            "score", "strand", "haplotype", "partition",
        ),
        (
            (
                row["source_contig"], row["source_start"], row["source_end"],
                row["record_id"], 0, "+", row["haplotype"], row["partition"],
            )
            for row in connection.execute(
                """
                SELECT * FROM mapped_regions
                ORDER BY haplotype, source_contig, source_start, partition
                """
            )
        ),
    )
    write_bed_table(
        os.path.join(output_dir, "reference_unique_regions.bed"),
        (
            "source_contig", "source_start", "source_end", "region_name",
            "score", "strand", "haplotype", "partition", "record_id",
        ),
        (
            (
                row["source_contig"], row["source_start"], row["source_end"],
                row["region_name"], 0, "+", row["haplotype"],
                row["partition"], row["record_id"],
            )
            for row in connection.execute(
                """
                SELECT * FROM unique_regions WHERE region_class='reference'
                ORDER BY haplotype, source_contig, source_start, partition
                """
            )
        ),
    )
    for filename, table in (
        ("header_union.bed", "header_union"),
        ("global_mapped_union.bed", "global_mapped"),
    ):
        write_bed_table(
            os.path.join(output_dir, filename),
            (
                "source_contig", "source_start", "source_end", "region_id",
                "score", "strand", "haplotype",
            ),
            (
                (
                    row["source_contig"], row["source_start"], row["source_end"],
                    f"{table}_{index:09d}", 0, "+", row["haplotype"],
                )
                for index, row in enumerate(connection.execute(
                    f"""
                    SELECT * FROM {table}
                    ORDER BY haplotype, source_contig, source_start
                    """
                ), 1)
            ),
        )
    write_bed_table(
        os.path.join(output_dir, "local_novel_candidates.bed"),
        (
            "source_contig", "source_start", "source_end", "candidate_id",
            "score", "strand", "haplotype", "novel_partition",
            "mapped_partitions", "record_id", "local_start", "local_end",
            "region_class",
        ),
        (
            (
                row["source_contig"], row["source_start"], row["source_end"],
                row["candidate_id"], 0, "+", row["haplotype"],
                row["novel_partition"],
                ",".join(candidate_memberships.get(
                    int(row["candidate_id"]), (),
                )),
                row["record_id"], row["local_start"], row["local_end"],
                row["region_class"],
            )
            for row in connection.execute(
                """
                SELECT * FROM local_candidates
                ORDER BY haplotype, source_contig, source_start, novel_partition
                """
            )
        ),
    )


def export_global_novel(connection: sqlite3.Connection, output_dir: str) -> None:
    header_coverage: CoverageIndex = {}
    current_key: Optional[Tuple[str, str]] = None
    intervals: List[Tuple[str, int, int]] = []
    for header in connection.execute(
        """
        SELECT haplotype, source_contig, partition, source_start, source_end
        FROM headers
        ORDER BY haplotype, source_contig, source_start, source_end, partition
        """
    ):
        key = str(header["haplotype"]), str(header["source_contig"])
        if current_key is not None and key != current_key:
            header_coverage[current_key] = coverage_segments(intervals)
            intervals = []
        current_key = key
        intervals.append((
            str(header["partition"]), int(header["source_start"]),
            int(header["source_end"]),
        ))
    if current_key is not None:
        header_coverage[current_key] = coverage_segments(intervals)

    path = os.path.join(output_dir, "global_novel.bed")
    rows = []
    segment_cursors: Dict[Tuple[str, str], int] = defaultdict(int)
    for index, row in enumerate(connection.execute(
        """
        SELECT * FROM global_novel
        ORDER BY haplotype, source_contig, source_start
        """
    ), 1):
        key = str(row["haplotype"]), str(row["source_contig"])
        segments = header_coverage.get(key, ())
        cursor = segment_cursors[key]
        start, end = int(row["source_start"]), int(row["source_end"])
        while cursor < len(segments) and segments[cursor].end <= start:
            cursor += 1
        segment_cursors[key] = cursor
        partitions: Set[str] = set()
        scan = cursor
        while scan < len(segments) and segments[scan].start < end:
            partitions.update(segments[scan].partitions)
            scan += 1
        rows.append((
            row["source_contig"], row["source_start"], row["source_end"],
            f"global_novel_{index:09d}", 0, "+", row["haplotype"],
            ",".join(sorted(partitions)),
        ))
    write_bed_table(
        path,
        (
            "source_contig", "source_start", "source_end", "global_novel_id",
            "score", "strand", "haplotype", "partitions",
        ),
        rows,
    )


def export_all_placements(connection: sqlite3.Connection, output_dir: str) -> None:
    fields = [field.name for field in dataclasses.fields(Placement)]
    write_tsv(
        os.path.join(output_dir, "all_local_novel_alignments.tsv"),
        fields,
        (
            tuple(row[field] for field in fields)
            for row in connection.execute(
                """
                SELECT * FROM all_placements
                ORDER BY candidate_id, distance, alignment_score DESC
                """
            )
        ),
    )


def export_local_novel(
    connection: sqlite3.Connection,
    output_dir: str,
    merge_gap: int,
) -> None:
    candidate_memberships = mapped_partitions_by_candidate(connection)
    rows = list(connection.execute(
        """
        SELECT c.*, b.*
        FROM local_candidates c
        JOIN best_placements b USING(candidate_id)
        ORDER BY b.ref_haplotype, b.ref_contig, b.ref_start, b.ref_end,
                 c.candidate_id
        """
    ))
    merged: List[Dict[str, object]] = []
    for row in rows:
        if (
            merged
            and merged[-1]["ref_haplotype"] == row["ref_haplotype"]
            and merged[-1]["ref_contig"] == row["ref_contig"]
            and int(row["ref_start"]) - int(merged[-1]["ref_end"]) <= merge_gap
        ):
            item = merged[-1]
            item["ref_end"] = max(int(item["ref_end"]), int(row["ref_end"]))
            item["alignment_score"] = max(
                float(item["alignment_score"]), float(row["alignment_score"]),
            )
            item["distance"] = min(int(item["distance"]), int(row["distance"]))
            item["candidate_ids"].append(str(row["candidate_id"]))
            item["source_loci"].append(
                f"{row['haplotype']}:{row['source_contig']}:"
                f"{row['source_start']}-{row['source_end']}"
            )
            item["novel_partitions"].add(str(row["novel_partition"]))
            item["mapped_partitions"].add(str(row["mapped_partition"]))
            item["strands"].add(str(row["strand"]))
            item["sources"].add(str(row["source"]))
            continue
        merged.append({
            "ref_haplotype": str(row["ref_haplotype"]),
            "ref_contig": str(row["ref_contig"]),
            "ref_start": int(row["ref_start"]),
            "ref_end": int(row["ref_end"]),
            "alignment_score": float(row["alignment_score"]),
            "distance": int(row["distance"]),
            "candidate_ids": [str(row["candidate_id"])],
            "source_loci": [
                f"{row['haplotype']}:{row['source_contig']}:"
                f"{row['source_start']}-{row['source_end']}"
            ],
            "novel_partitions": {str(row["novel_partition"])},
            "mapped_partitions": {str(row["mapped_partition"])},
            "strands": {str(row["strand"])},
            "sources": {str(row["source"])},
        })
    output_rows = []
    for index, item in enumerate(merged, 1):
        output_rows.append((
            item["ref_contig"], item["ref_start"], item["ref_end"],
            f"local_novel_{index:09d}", f"{item['alignment_score']:.1f}",
            ",".join(sorted(item["strands"])), item["ref_haplotype"],
            ";".join(item["source_loci"]),
            ",".join(sorted(item["novel_partitions"])),
            ",".join(sorted(item["mapped_partitions"])),
            item["distance"], ",".join(item["candidate_ids"]),
            ",".join(sorted(item["sources"])),
        ))
    write_bed_table(
        os.path.join(output_dir, "local_novel.bed"),
        (
            "ref_contig", "ref_start", "ref_end", "local_novel_id",
            "alignment_score", "strand", "ref_haplotype", "source_loci",
            "novel_partitions", "mapped_partitions",
            "distance_to_reference_unique", "candidate_ids",
            "alignment_sources",
        ),
        output_rows,
    )
    write_bed_table(
        os.path.join(output_dir, "unresolved_local_novel.bed"),
        (
            "source_contig", "source_start", "source_end", "candidate_id",
            "score", "strand", "haplotype", "novel_partition",
            "mapped_partitions",
        ),
        (
            (
                row["source_contig"], row["source_start"], row["source_end"],
                row["candidate_id"], 0, "+", row["haplotype"],
                row["novel_partition"],
                ",".join(candidate_memberships.get(
                    int(row["candidate_id"]), (),
                )),
            )
            for row in connection.execute(
                """
                SELECT c.* FROM local_candidates c
                LEFT JOIN best_placements b USING(candidate_id)
                WHERE b.candidate_id IS NULL
                ORDER BY c.haplotype, c.source_contig, c.source_start
                """
            )
        ),
    )


def parse_args(argv: Optional[Sequence[str]] = None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=(
            "Aggregate local minsetref-light results into global-novel and "
            "reference-projected local-novel BED files."
        )
    )
    parser.add_argument(
        "-i", "--localblocks", required=True,
        help="root containing PARTITION/PARTITION.{fasta,bed,unique.bed}",
    )
    parser.add_argument("-o", "--output-dir", required=True)
    parser.add_argument(
        "-j", "--jobs", type=int, default=4,
        help=(
            "parallel partition discovery/loading, local interval sweeps, "
            "and alignment jobs (default: 4)"
        ),
    )
    parser.add_argument("-c", "--cores", type=int, default=4)
    parser.add_argument(
        "--input-anchor", type=int, default=50,
        help="reporting anchor to remove from unique BED sides (default: 50)",
    )
    parser.add_argument(
        "--merge-gap", type=int, default=500,
        help="merge projected Ref intervals within this gap (default: 500)",
    )
    parser.add_argument("--min-score", type=int, default=100)
    parser.add_argument("--min-identity", type=float, default=95.0)
    parser.add_argument("--blast-word-size", type=int, default=50)
    parser.add_argument("--blast-evalue", default="1e-300")
    parser.add_argument("--blast-max-target-seqs", type=int, default=100)
    parser.add_argument("--skip-blastn", action="store_true")
    parser.add_argument("--resume", action="store_true")
    parser.add_argument("--rebuild-catalog", action="store_true")
    parser.add_argument(
        "--log-level", choices=("DEBUG", "INFO", "WARNING", "ERROR"),
        default="INFO",
    )
    args = parser.parse_args(argv)
    for name in ("jobs", "cores", "min_score"):
        if getattr(args, name) < 1:
            parser.error(f"--{name.replace('_', '-')} must be positive")
    if args.input_anchor < 0 or args.merge_gap < 0:
        parser.error("--input-anchor and --merge-gap cannot be negative")
    if not 0 < args.min_identity <= 100:
        parser.error("--min-identity must be in (0,100]")
    return args


def run(args: argparse.Namespace) -> None:
    output_dir = os.path.abspath(args.output_dir)
    Path(output_dir).mkdir(parents=True, exist_ok=True)
    database_path = os.path.join(output_dir, "partition_novelty.sqlite")
    connection = connect_database(database_path, create_indexes=False)
    partitions = discover_partitions(args.localblocks, jobs=args.jobs)
    expected_catalog_signature = catalog_signature(
        partitions, args.input_anchor,
    )
    catalog_complete = connection.execute(
        "SELECT value FROM metadata WHERE key='catalog_complete'"
    ).fetchone()
    stored_signature = connection.execute(
        "SELECT value FROM metadata WHERE key='catalog_signature'"
    ).fetchone()
    if (
        args.resume
        and catalog_complete
        and (
            stored_signature is None
            or stored_signature["value"] != expected_catalog_signature
        )
        and not args.rebuild_catalog
    ):
        raise RuntimeError(
            "saved interval catalog does not match current inputs or "
            "--input-anchor; rerun with --rebuild-catalog"
        )
    if args.rebuild_catalog or not (args.resume and catalog_complete):
        reset_catalog(connection)
        LOG.info(
            "Loading %d local-block partitions with %d workers",
            len(partitions), args.jobs,
        )
        load_partitions_parallel(
            connection, partitions, args.input_anchor, args.jobs,
        )
        LOG.info("Building global interval unions with sorted sweeps")
        build_interval_unions(connection)
        LOG.info("Building per-partition local candidates")
        build_local_candidates(connection, jobs=args.jobs)
        LOG.info("Creating SQLite indexes after bulk catalog loading")
        create_secondary_indexes(connection)
        connection.execute(
            "INSERT OR REPLACE INTO metadata VALUES ('catalog_complete', '1')"
        )
        connection.execute(
            "INSERT OR REPLACE INTO metadata VALUES ('catalog_signature', ?)",
            (expected_catalog_signature,),
        )
        connection.commit()
    else:
        create_secondary_indexes(connection)
        LOG.info("RESUME: using completed interval catalog %s", database_path)

    export_catalog(connection, output_dir)
    export_global_novel(connection, output_dir)

    reference_haplotypes = {
        str(row[0]) for row in connection.execute(
            "SELECT DISTINCT haplotype FROM unique_regions WHERE region_class='reference'"
        )
    }
    mapped_partitions = [
        str(row[0]) for row in connection.execute(
            """
            SELECT DISTINCT cmp.mapped_partition
            FROM candidate_mapped_partitions cmp
            JOIN local_candidates c USING(candidate_id)
            WHERE c.haplotype NOT IN (
                SELECT DISTINCT haplotype FROM unique_regions
                WHERE region_class='reference'
            )
            ORDER BY cmp.mapped_partition
            """
        )
    ]
    LOG.info(
        "Local-novel candidates: %d; Ref haplotypes: %s; alignment partitions: %d",
        connection.execute("SELECT COUNT(*) FROM local_candidates").fetchone()[0],
        ",".join(sorted(reference_haplotypes)) or "none",
        len(mapped_partitions),
    )
    if mapped_partitions:
        for executable in ("minimap2", "winnowmap"):
            ensure_executable(executable)
        if not args.skip_blastn:
            for executable in ("makeblastdb", "blastn"):
                ensure_executable(executable)
        options = {
            "resume": args.resume,
            "cores": args.cores,
            "anchor": args.input_anchor,
            "min_identity": args.min_identity,
            "min_score": args.min_score,
            "blast_word_size": args.blast_word_size,
            "blast_evalue": args.blast_evalue,
            "blast_max_target_seqs": args.blast_max_target_seqs,
            "skip_blastn": args.skip_blastn,
            "catalog_signature": expected_catalog_signature,
        }
        work_root = os.path.join(output_dir, "alignment_work")
        with ThreadPoolExecutor(max_workers=args.jobs) as executor:
            result_paths = list(executor.map(
                lambda partition: alignment_task(
                    database_path, partition, work_root, options,
                ),
                mapped_partitions,
            ))
    else:
        result_paths = []

    connection.execute("DELETE FROM all_placements")
    connection.execute("DELETE FROM best_placements")
    for path in result_paths:
        load_placement_file(connection, path)
    add_reference_source_placements(connection)
    select_best_placements(connection)
    export_all_placements(connection, output_dir)
    export_local_novel(connection, output_dir, args.merge_gap)
    connection.close()
    LOG.info("Wrote global and local novelty outputs under %s", output_dir)


def main(argv: Optional[Sequence[str]] = None) -> int:
    args = parse_args(argv)
    logging.basicConfig(
        level=getattr(logging, args.log_level),
        format="[%(asctime)s] %(levelname)s: %(message)s",
        datefmt="%Y-%m-%d %H:%M:%S",
    )
    try:
        run(args)
    except (OSError, ValueError, RuntimeError, sqlite3.Error,
            subprocess.CalledProcessError) as error:
        LOG.error("%s", error)
        return 1
    return 0


if __name__ == "__main__":
    sys.exit(main())
