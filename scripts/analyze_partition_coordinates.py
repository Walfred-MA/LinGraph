#!/usr/bin/env python3
"""Classify partition novelty with one in-memory coordinate sweep.

For each (haplotype, source contig), the script jointly sweeps:

* H: partitions whose input FASTA headers cover the source interval;
* U: partitions whose anchor-trimmed ``unique.bed`` covers it;
* B: partitions whose numbered ``PARTITION.bed`` uses it as a novel locus.

The complete coordinate analysis is written to flat BED/TSV files.  Local
novel alignment is deliberately handled by ``map_partition_local_novel.py`` so
the large in-memory interval catalog is released before aligners are started.
"""
from __future__ import annotations

import argparse
import csv
import dataclasses
import hashlib
import json
import logging
import os
import sys
from array import array
from collections import defaultdict, deque
from concurrent.futures import ProcessPoolExecutor
from pathlib import Path
from typing import (
    DefaultDict,
    Dict,
    FrozenSet,
    Iterable,
    Iterator,
    List,
    Optional,
    Sequence,
    Set,
    Tuple,
)

from minsetref_core import mp_context
from minsetref_light import _translate_coordinate_reference_interval
from minsetref_segments import complement_intervals
from partition_hotspot_headers import (
    CompactHeaderStore,
    read_blocks,
    read_query_sources,
    remap_blocks_to_target_list,
    reconstruct_headers,
    source_fingerprints,
)
from summarize_partition_novelty import (
    HeaderRecord,
    PartitionFiles,
    discover_partitions,
    parse_header_records,
    partition_header_path,
    remove_reporting_anchors,
    threaded_map_unordered,
    write_bed_table,
    write_tsv,
)


LOG = logging.getLogger("analyze_partition_coordinates")
ContigKey = Tuple[str, str]
Interval = Tuple[int, int]
IGNORED_SOURCE_CONTIGS = frozenset({"chrM"})


@dataclasses.dataclass(frozen=True)
class UniqueRegion:
    partition: str
    record_id: str
    haplotype: str
    source_contig: str
    source_start: int
    source_end: int
    local_start: int
    local_end: int
    raw_local_start: int
    raw_local_end: int
    region_name: str
    region_class: str
    novelty_score: float
    header_index: int = -1


@dataclasses.dataclass(frozen=True)
class UsedNovelRegion:
    partition: str
    haplotype: str
    source_contig: str
    source_start: int
    source_end: int
    encoded_name: str
    local_start: int
    local_end: int
    block_id: str
    block_score: str


@dataclasses.dataclass
class CoordinateCatalog:
    files: PartitionFiles
    headers: List[HeaderRecord]
    unique_regions: List[UniqueRegion]
    used_novel_regions: List[UsedNovelRegion]


@dataclasses.dataclass
class ContigIntervals:
    headers: List[Tuple[int, int, str]]
    unique_regions: List[UniqueRegion]
    used_novel_regions: List[UsedNovelRegion]


@dataclasses.dataclass
class CompactContigIntervals:
    holder_id: int
    header_indices: object
    unique_regions: List[UniqueRegion]
    used_novel_regions: List[UsedNovelRegion]


@dataclasses.dataclass
class HeaderFileContigIntervals:
    source_start: array
    source_end: array
    partition_id: array
    unique_regions: List[UniqueRegion]
    used_novel_regions: List[UsedNovelRegion]


@dataclasses.dataclass
class AnnotationCatalog:
    files: PartitionFiles
    unique_regions: List[UniqueRegion]
    used_novel_regions: List[UsedNovelRegion]


@dataclasses.dataclass(frozen=True)
class StatusSpan:
    start: int
    end: int
    status: str
    header_partitions: FrozenSet[str]
    lifted_partitions: FrozenSet[str]
    unlifted_partitions: FrozenSet[str]
    used_novel_locus: FrozenSet[str]


@dataclasses.dataclass(frozen=True)
class CandidateSpan:
    status_index: int
    start: int
    end: int
    novel_partition: str
    lifted_partitions: FrozenSet[str]
    unique_index: int


@dataclasses.dataclass
class SweepResult:
    key: ContigKey
    statuses: List[StatusSpan]
    candidates: List[CandidateSpan]


_CONTIG_INTERVALS: Dict[ContigKey, ContigIntervals] = {}
_COMPACT_HEADER_STORE: Optional[CompactHeaderStore] = None
_REFERENCE_HAPLOTYPES: FrozenSet[str] = frozenset()
_HEADER_PARTITIONS: Sequence[str] = ()


def parse_reference_haplotypes(text: str) -> FrozenSet[str]:
    values = frozenset(token.strip() for token in text.split(",") if token.strip())
    if not values:
        raise ValueError("--reference-haplotypes cannot be empty")
    return values


def split_coordinate_source(
    source_name: str,
    known_haplotypes: Iterable[str],
    context: str,
) -> Tuple[str, str]:
    del known_haplotypes  # Validation is cohort-wide after all headers load.
    fields = source_name.split("_", 2)
    if len(fields) != 3 or not all(fields):
        raise ValueError(
            f"{context}: coordinate name source {source_name!r} does not match "
            "HAPLOTYPE_FIELD1_HAPLOTYPE_FIELD2_CONTIG"
        )
    haplotype = f"{fields[0]}_{fields[1]}"
    contig = fields[2]
    return haplotype, contig


def parse_used_novel_regions(
    partition: PartitionFiles,
    known_haplotypes: Iterable[str],
) -> List[UsedNovelRegion]:
    rows: List[UsedNovelRegion] = []
    with open(partition.reference_bed, "rt") as handle:
        for line_number, raw in enumerate(handle, 1):
            if not raw.strip() or raw.startswith("#"):
                continue
            fields = raw.split()
            if len(fields) < 3:
                raise ValueError(
                    f"{partition.reference_bed}:{line_number}: expected BED3+"
                )
            try:
                local_start, local_end = int(fields[1]), int(fields[2])
            except ValueError as error:
                raise ValueError(
                    f"{partition.reference_bed}:{line_number}: invalid BED coordinates"
                ) from error
            if (
                len(fields) >= 11
                and fields[10] == "source_coordinates"
            ):
                try:
                    fasta_start, fasta_end = int(fields[8]), int(fields[9])
                except ValueError as error:
                    raise ValueError(
                        f"{partition.reference_bed}:{line_number}: invalid "
                        "anchored FASTA coordinates"
                    ) from error
                if not (
                    fasta_start <= local_start < local_end <= fasta_end
                    and fields[0]
                    and fields[6]
                    and fields[7]
                ):
                    raise ValueError(
                        f"{partition.reference_bed}:{line_number}: core "
                        f"{local_start}-{local_end} is outside anchored FASTA "
                        f"interval {fasta_start}-{fasta_end}"
                    )
                rows.append(UsedNovelRegion(
                    partition.partition,
                    fields[6],
                    fields[0],
                    local_start,
                    local_end,
                    fields[7],
                    local_start - fasta_start,
                    local_end - fasta_start,
                    fields[3],
                    fields[4],
                ))
                continue
            translated = _translate_coordinate_reference_interval(
                fields[0], local_start, local_end,
                f"{partition.reference_bed}:{line_number}",
                (
                    "source"
                    if len(fields) >= 7
                    and fields[6] == "source_coordinates"
                    else None
                ),
            )
            if translated is None:
                raise ValueError(
                    f"{partition.reference_bed}:{line_number}: numbered block BED "
                    f"contig {fields[0]!r} is not coordinate-named; expected "
                    "HAPLOTYPE_CONTIG_SOURCESTART_SOURCEEND. Run "
                    "temp_fixes/fix_existing_block_bed_names.py first."
                )
            (
                source_name,
                source_start,
                source_end,
                local_start,
                local_end,
            ) = translated
            haplotype, source_contig = split_coordinate_source(
                source_name, known_haplotypes,
                f"{partition.reference_bed}:{line_number}",
            )
            rows.append(UsedNovelRegion(
                partition.partition,
                haplotype,
                source_contig,
                source_start,
                source_end,
                fields[0],
                local_start,
                local_end,
                fields[3] if len(fields) > 3 else ".",
                fields[4] if len(fields) > 4 else ".",
            ))
    if not rows:
        raise ValueError(f"{partition.reference_bed}: no BED intervals")
    return rows


def parse_unique_regions(
    partition: PartitionFiles,
    headers: Sequence[HeaderRecord],
    anchor: int,
) -> List[UniqueRegion]:
    by_id = {record.record_id: record for record in headers}
    rows: List[UniqueRegion] = []
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
            try:
                raw_start, raw_end = int(fields[1]), int(fields[2])
                score = float(fields[8])
            except ValueError as error:
                raise ValueError(
                    f"{partition.unique_bed}:{line_number}: invalid numeric field"
                ) from error
            local_start, local_end = remove_reporting_anchors(
                raw_start, raw_end, record.record_length, anchor,
            )
            if fields[6] != record.haplotype:
                raise ValueError(
                    f"{partition.unique_bed}:{line_number}: haplotype "
                    f"{fields[6]!r} disagrees with {record.haplotype!r}"
                )
            rows.append(UniqueRegion(
                partition.partition,
                record.record_id,
                record.haplotype,
                record.source_contig,
                record.source_start + local_start,
                record.source_start + local_end,
                local_start,
                local_end,
                raw_start,
                raw_end,
                fields[3],
                fields[7],
                score,
            ))
    return rows


def parse_coordinate_catalog(
    partition: PartitionFiles,
    anchor: int,
    ignored_source_contigs: FrozenSet[str] = IGNORED_SOURCE_CONTIGS,
) -> CoordinateCatalog:
    all_headers = parse_header_records(partition)
    haplotypes = {record.haplotype for record in all_headers}
    headers = [
        record for record in all_headers
        if record.source_contig not in ignored_source_contigs
    ]
    unique_regions = [
        row for row in parse_unique_regions(
            partition, all_headers, anchor,
        )
        if row.source_contig not in ignored_source_contigs
    ]
    used_novel_regions = [
        row for row in parse_used_novel_regions(partition, haplotypes)
        if row.source_contig not in ignored_source_contigs
    ]
    return CoordinateCatalog(
        partition,
        headers,
        unique_regions,
        used_novel_regions,
    )


def load_catalogs(
    partitions: Sequence[PartitionFiles],
    anchor: int,
    jobs: int,
    ignored_source_contigs: FrozenSet[str] = IGNORED_SOURCE_CONTIGS,
) -> List[CoordinateCatalog]:
    def parse(partition: PartitionFiles) -> CoordinateCatalog:
        return parse_coordinate_catalog(
            partition, anchor, ignored_source_contigs,
        )

    catalogs: List[CoordinateCatalog] = []
    progress_step = max(1, len(partitions) // 20)
    for index, catalog in enumerate(
        threaded_map_unordered(parse, partitions, jobs), 1,
    ):
        catalogs.append(catalog)
        if index == len(partitions) or index % progress_step == 0:
            LOG.info("Loaded %d/%d partitions", index, len(partitions))
    catalogs.sort(key=lambda item: item.files.partition)
    return catalogs


def extracted_header_path(partition: PartitionFiles) -> str:
    return partition_header_path(partition)


def validate_extracted_headers(
    partitions: Sequence[PartitionFiles],
    jobs: int,
) -> Tuple[Dict[str, Tuple[str, int]], List[Tuple[str, int, int]]]:
    """Require nonempty extracted .header files without auxiliary markers."""

    def inspect(
        partition: PartitionFiles,
    ) -> Tuple[str, str, int, Tuple[Tuple[str, int, int], ...]]:
        header = extracted_header_path(partition)
        if not os.path.isfile(header):
            raise FileNotFoundError(header)
        if os.path.getsize(header) <= 0:
            raise ValueError(f"{header}: empty extracted header file")
        file_stat = os.stat(header)
        fingerprints = ((
            os.path.abspath(header),
            file_stat.st_size,
            file_stat.st_mtime_ns,
        ),)
        return (
            partition.partition,
            os.path.abspath(header),
            0,
            fingerprints,
        )

    locations: Dict[str, Tuple[str, int]] = {}
    fingerprints: List[Tuple[str, int, int]] = []
    for partition, header, count, rows in threaded_map_unordered(
        inspect, partitions, jobs,
    ):
        locations[partition] = header, count
        fingerprints.extend(rows)
    return locations, fingerprints


def load_header_file_catalogs(
    partitions: Sequence[PartitionFiles],
    header_locations: Dict[str, Tuple[str, int]],
    anchor: int,
    jobs: int,
    ignored_source_contigs: FrozenSet[str] = IGNORED_SOURCE_CONTIGS,
) -> Tuple[
    List[AnnotationCatalog],
    Dict[ContigKey, HeaderFileContigIntervals],
    Set[str],
    int,
]:
    """Stream authoritative .header rows into compact per-contig arrays."""
    partition_names = [partition.partition for partition in partitions]
    partition_ids = {
        partition: index for index, partition in enumerate(partition_names)
    }
    contig_index: Dict[ContigKey, HeaderFileContigIntervals] = {}
    contig_owner: Dict[str, str] = {}
    observed_haplotypes: Set[str] = set()
    annotations: List[AnnotationCatalog] = []
    header_count = 0

    def parse(partition: PartitionFiles) -> CoordinateCatalog:
        return parse_coordinate_catalog(
            partition, anchor, ignored_source_contigs,
        )

    progress_step = max(1, len(partitions) // 20)
    for completed, catalog in enumerate(
        threaded_map_unordered(parse, partitions, jobs), 1,
    ):
        header_path = header_locations[catalog.files.partition][0]
        header_locations[catalog.files.partition] = (
            header_path,
            len(catalog.headers),
        )
        partition_id = partition_ids[catalog.files.partition]
        for header in catalog.headers:
            owner = contig_owner.get(header.source_contig)
            if owner is None:
                contig_owner[header.source_contig] = header.haplotype
            elif owner != header.haplotype:
                raise ValueError(
                    f"source contig {header.source_contig!r} occurs in both "
                    f"{owner!r} and {header.haplotype!r}; contig names must "
                    "be globally unique across assemblies"
                )
            key = header.haplotype, header.source_contig
            bucket = contig_index.get(key)
            if bucket is None:
                bucket = HeaderFileContigIntervals(
                    array("Q"), array("Q"), array("I"), [], [],
                )
                contig_index[key] = bucket
            bucket.source_start.append(header.source_start)
            bucket.source_end.append(header.source_end)
            bucket.partition_id.append(partition_id)
            observed_haplotypes.add(header.haplotype)
        for row in catalog.unique_regions:
            key = row.haplotype, row.source_contig
            bucket = contig_index.get(key)
            if bucket is None:
                bucket = HeaderFileContigIntervals(
                    array("Q"), array("Q"), array("I"), [], [],
                )
                contig_index[key] = bucket
            bucket.unique_regions.append(row)
        for row in catalog.used_novel_regions:
            key = row.haplotype, row.source_contig
            bucket = contig_index.get(key)
            if bucket is None:
                bucket = HeaderFileContigIntervals(
                    array("Q"), array("Q"), array("I"), [], [],
                )
                contig_index[key] = bucket
            bucket.used_novel_regions.append(row)
        header_count += len(catalog.headers)
        annotations.append(AnnotationCatalog(
            catalog.files,
            catalog.unique_regions,
            catalog.used_novel_regions,
        ))
        if completed == len(partitions) or completed % progress_step == 0:
            LOG.info(
                "Loaded extracted headers for %d/%d partitions: %d headers",
                completed, len(partitions), header_count,
            )
    annotations.sort(key=lambda item: item.files.partition)
    return (
        annotations,
        contig_index,
        observed_haplotypes,
        header_count,
    )


def header_record_from_store(
    store: CompactHeaderStore,
    header_index: int,
) -> HeaderRecord:
    return HeaderRecord(
        store.partition_name(header_index),
        store.record_id(header_index),
        store.haplotype(header_index),
        store.contig(header_index),
        store.start(header_index),
        store.end(header_index),
        int(store.record_length(header_index)),
        0,
    )


def parse_compact_annotation_catalog(
    partition: PartitionFiles,
    store: CompactHeaderStore,
    anchor: int,
    ignored_source_contigs: FrozenSet[str] = IGNORED_SOURCE_CONTIGS,
) -> AnnotationCatalog:
    unique_rows: List[UniqueRegion] = []
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
            try:
                header_index = store.find_record(
                    partition.partition, fields[0], fields[6],
                )
            except KeyError as error:
                raise ValueError(
                    f"{partition.unique_bed}:{line_number}: reconstructed "
                    f"headers contain no record {fields[0]!r}; check hotspot "
                    "reconstruction parameters"
                ) from error
            record_length = int(store.record_length(header_index))
            raw_start, raw_end = int(fields[1]), int(fields[2])
            local_start, local_end = remove_reporting_anchors(
                raw_start, raw_end, record_length, anchor,
            )
            haplotype = store.haplotype(header_index)
            if fields[6] != haplotype:
                raise ValueError(
                    f"{partition.unique_bed}:{line_number}: haplotype "
                    f"{fields[6]!r} disagrees with reconstructed "
                    f"{haplotype!r}"
                )
            source_contig = store.contig(header_index)
            if source_contig in ignored_source_contigs:
                continue
            source_start = store.start(header_index)
            unique_rows.append(UniqueRegion(
                partition.partition,
                fields[0],
                haplotype,
                source_contig,
                source_start + local_start,
                source_start + local_end,
                local_start,
                local_end,
                raw_start,
                raw_end,
                fields[3],
                fields[7],
                float(fields[8]),
                header_index,
            ))
    return AnnotationCatalog(
        partition,
        unique_rows,
        [
            row
            for row in parse_used_novel_regions(
                partition, store.haplotypes,
            )
            if row.source_contig not in ignored_source_contigs
        ],
    )


def load_compact_annotations(
    partitions: Sequence[PartitionFiles],
    store: CompactHeaderStore,
    anchor: int,
    jobs: int,
    ignored_source_contigs: FrozenSet[str] = IGNORED_SOURCE_CONTIGS,
) -> List[AnnotationCatalog]:
    def parse(partition: PartitionFiles) -> AnnotationCatalog:
        return parse_compact_annotation_catalog(
            partition, store, anchor, ignored_source_contigs,
        )

    catalogs: List[AnnotationCatalog] = []
    progress_step = max(1, len(partitions) // 20)
    for index, catalog in enumerate(
        threaded_map_unordered(parse, partitions, jobs), 1,
    ):
        catalogs.append(catalog)
        if index == len(partitions) or index % progress_step == 0:
            LOG.info(
                "Loaded unique/BED annotations for %d/%d partitions",
                index, len(partitions),
            )
    catalogs.sort(key=lambda item: item.files.partition)
    return catalogs


def build_compact_contig_index(
    store: CompactHeaderStore,
    catalogs: Sequence[AnnotationCatalog],
    ignored_source_contigs: FrozenSet[str] = IGNORED_SOURCE_CONTIGS,
) -> Dict[ContigKey, CompactContigIntervals]:
    index: Dict[ContigKey, CompactContigIntervals] = {}
    for key, holder_id, header_indices in store.iter_contigs():
        if key[1] in ignored_source_contigs:
            continue
        index[key] = CompactContigIntervals(
            holder_id, header_indices, [], [],
        )

    def bucket(key: ContigKey) -> CompactContigIntervals:
        value = index.get(key)
        if value is None:
            value = CompactContigIntervals(-1, (), [], [])
            index[key] = value
        return value

    for catalog in catalogs:
        for row in catalog.unique_regions:
            bucket((row.haplotype, row.source_contig)).unique_regions.append(row)
        for row in catalog.used_novel_regions:
            bucket((
                row.haplotype, row.source_contig,
            )).used_novel_regions.append(row)
    return index


def build_contig_index(
    catalogs: Sequence[CoordinateCatalog],
) -> Dict[ContigKey, ContigIntervals]:
    index: Dict[ContigKey, ContigIntervals] = {}

    def bucket(key: ContigKey) -> ContigIntervals:
        value = index.get(key)
        if value is None:
            value = ContigIntervals([], [], [])
            index[key] = value
        return value

    for catalog in catalogs:
        for row in catalog.headers:
            bucket((row.haplotype, row.source_contig)).headers.append(
                (row.source_start, row.source_end, row.partition)
            )
        for row in catalog.unique_regions:
            bucket((row.haplotype, row.source_contig)).unique_regions.append(row)
        for row in catalog.used_novel_regions:
            bucket((row.haplotype, row.source_contig)).used_novel_regions.append(row)
    return index


def _update_count(
    counts: Dict[str, int],
    partition: str,
    delta: int,
    context: str,
) -> None:
    value = counts.get(partition, 0) + delta
    if value < 0:
        raise RuntimeError(f"{context}: negative active count for {partition}")
    if value:
        counts[partition] = value
    else:
        counts.pop(partition, None)


def _classify(
    haplotype: str,
    headers: FrozenSet[str],
    unique: FrozenSet[str],
) -> Tuple[str, FrozenSet[str], FrozenSet[str]]:
    if not unique.issubset(headers):
        raise ValueError(
            f"{haplotype}: unique partitions absent from header coverage: "
            + ",".join(sorted(unique - headers))
        )
    lifted = headers - unique
    unlifted = headers & unique
    if haplotype in _REFERENCE_HAPLOTYPES:
        # The true reference is not called novel against itself.  Report only
        # whether its source coordinate is represented by any partition
        # header; internal min-set uniqueness remains available in the
        # lifted/unlifted membership columns.
        return ("mapped" if headers else "unmapped"), lifted, unlifted
    if not headers:
        return "unmapped", lifted, unlifted
    if unlifted and lifted:
        return "local_novel", lifted, unlifted
    if unlifted:
        return "global_novel", lifted, unlifted
    return "mapped", lifted, unlifted


def _append_status(
    output: List[StatusSpan],
    start: int,
    end: int,
    status: str,
    headers: FrozenSet[str],
    lifted: FrozenSet[str],
    unlifted: FrozenSet[str],
    used: FrozenSet[str],
) -> int:
    if (
        output
        and output[-1].end == start
        and output[-1].status == status
        and output[-1].header_partitions == headers
        and output[-1].lifted_partitions == lifted
        and output[-1].unlifted_partitions == unlifted
        and output[-1].used_novel_locus == used
    ):
        previous = output[-1]
        output[-1] = dataclasses.replace(previous, end=end)
        return len(output) - 1
    output.append(StatusSpan(
        start, end, status, headers, lifted, unlifted, used,
    ))
    return len(output) - 1


def _sweep_materialized(
    key: ContigKey,
    data: ContigIntervals,
) -> SweepResult:
    # event = (position, kind, partition, delta, unique_index)
    events: List[Tuple[int, str, str, int, int]] = []
    for start, end, partition in data.headers:
        events.append((start, "H", partition, 1, -1))
        events.append((end, "H", partition, -1, -1))
    for unique_index, row in enumerate(data.unique_regions):
        events.append((
            row.source_start, "U", row.partition, 1, unique_index,
        ))
        events.append((
            row.source_end, "U", row.partition, -1, unique_index,
        ))
    for row in data.used_novel_regions:
        events.append((row.source_start, "B", row.partition, 1, -1))
        events.append((row.source_end, "B", row.partition, -1, -1))
    if not events:
        return SweepResult(key, [], [])
    events.sort(key=lambda item: item[0])

    header_counts: Dict[str, int] = {}
    unique_counts: Dict[str, int] = {}
    used_counts: Dict[str, int] = {}
    active_unique: DefaultDict[str, Set[int]] = defaultdict(set)
    statuses: List[StatusSpan] = []
    candidates: List[CandidateSpan] = []
    index = 0
    while index < len(events):
        position = events[index][0]
        while index < len(events) and events[index][0] == position:
            _event_position, kind, partition, delta, unique_index = events[index]
            if kind == "H":
                _update_count(header_counts, partition, delta, str(key))
            elif kind == "B":
                _update_count(used_counts, partition, delta, str(key))
            else:
                _update_count(unique_counts, partition, delta, str(key))
                if delta > 0:
                    active_unique[partition].add(unique_index)
                else:
                    active_unique[partition].discard(unique_index)
                    if not active_unique[partition]:
                        active_unique.pop(partition, None)
            index += 1
        if index >= len(events):
            break
        next_position = events[index][0]
        if next_position <= position:
            continue
        headers = frozenset(header_counts)
        unique = frozenset(unique_counts)
        used = frozenset(used_counts)
        status, lifted, unlifted = _classify(key[0], headers, unique)
        status_index = _append_status(
            statuses, position, next_position, status,
            headers, lifted, unlifted, used,
        )
        if status != "local_novel":
            continue
        for novel_partition in sorted(unlifted):
            choices = sorted(
                active_unique.get(novel_partition, ()),
                key=lambda item: (
                    data.unique_regions[item].record_id,
                    data.unique_regions[item].source_start,
                    data.unique_regions[item].source_end,
                ),
            )
            if not choices:
                raise RuntimeError(
                    f"{key}:{position}-{next_position}: no active unique "
                    f"record for partition {novel_partition}"
                )
            unique_index = choices[0]
            if (
                candidates
                and candidates[-1].status_index == status_index
                and candidates[-1].end == position
                and candidates[-1].novel_partition == novel_partition
                and candidates[-1].lifted_partitions == lifted
                and candidates[-1].unique_index == unique_index
            ):
                candidates[-1] = dataclasses.replace(
                    candidates[-1], end=next_position,
                )
            else:
                candidates.append(CandidateSpan(
                    status_index, position, next_position,
                    novel_partition, lifted, unique_index,
                ))
    return SweepResult(key, statuses, candidates)


def sweep_contig(key: ContigKey) -> SweepResult:
    return _sweep_materialized(key, _CONTIG_INTERVALS[key])


def sweep_compact_contig(key: ContigKey) -> SweepResult:
    if _COMPACT_HEADER_STORE is None:
        raise RuntimeError("compact header store is not initialized")
    compact = _CONTIG_INTERVALS[key]
    store = _COMPACT_HEADER_STORE
    holder = (
        store.holders[compact.holder_id]
        if compact.holder_id >= 0 else None
    )
    if compact.holder_id >= 0 and holder is None:
        raise RuntimeError(f"missing assembly holder for {key}")
    headers = [
        (
            int(holder.source_start[header_index]),
            int(holder.source_end[header_index]),
            store.partitions[holder.partition_id[header_index]],
        )
        for header_index in compact.header_indices
    ]
    return _sweep_materialized(
        key,
        ContigIntervals(
            headers,
            compact.unique_regions,
            compact.used_novel_regions,
        ),
    )


def sweep_header_file_contig(key: ContigKey) -> SweepResult:
    data = _CONTIG_INTERVALS[key]
    headers = [
        (
            int(data.source_start[index]),
            int(data.source_end[index]),
            _HEADER_PARTITIONS[data.partition_id[index]],
        )
        for index in range(len(data.partition_id))
    ]
    return _sweep_materialized(
        key,
        ContigIntervals(
            headers,
            data.unique_regions,
            data.used_novel_regions,
        ),
    )


def _initialize_sweep_workers(
    index: Dict[ContigKey, ContigIntervals],
    reference_haplotypes: FrozenSet[str],
) -> None:
    global _CONTIG_INTERVALS, _REFERENCE_HAPLOTYPES
    _CONTIG_INTERVALS = index
    _REFERENCE_HAPLOTYPES = reference_haplotypes


def _initialize_header_file_sweep_workers(
    index,
    partition_names: Sequence[str],
    reference_haplotypes: FrozenSet[str],
) -> None:
    global _CONTIG_INTERVALS, _HEADER_PARTITIONS, _REFERENCE_HAPLOTYPES
    _CONTIG_INTERVALS = index
    _HEADER_PARTITIONS = partition_names
    _REFERENCE_HAPLOTYPES = reference_haplotypes


def _initialize_compact_sweep_workers(
    index,
    store: CompactHeaderStore,
    reference_haplotypes: FrozenSet[str],
) -> None:
    global _CONTIG_INTERVALS, _COMPACT_HEADER_STORE, _REFERENCE_HAPLOTYPES
    _CONTIG_INTERVALS = index
    _COMPACT_HEADER_STORE = store
    _REFERENCE_HAPLOTYPES = reference_haplotypes


def sweep_all_contigs(
    index: Dict[ContigKey, ContigIntervals],
    reference_haplotypes: FrozenSet[str],
    jobs: int,
) -> Iterator[SweepResult]:
    keys = sorted(index)
    _initialize_sweep_workers(index, reference_haplotypes)
    if jobs == 1:
        yield from map(sweep_contig, keys)
        return
    with ProcessPoolExecutor(
        max_workers=jobs,
        mp_context=mp_context(),
        initializer=_initialize_sweep_workers,
        initargs=(index, reference_haplotypes),
    ) as executor:
        # Keep only a small ordered window in flight.  This prevents one slow
        # chromosome from allowing tens of thousands of completed contig
        # results to accumulate behind it inside Executor.map().
        key_iterator = iter(keys)
        pending = deque()
        for _index in range(jobs * 2):
            try:
                key = next(key_iterator)
            except StopIteration:
                break
            pending.append(executor.submit(sweep_contig, key))
        while pending:
            yield pending.popleft().result()
            try:
                key = next(key_iterator)
            except StopIteration:
                continue
            pending.append(executor.submit(sweep_contig, key))


def sweep_all_compact_contigs(
    index: Dict[ContigKey, CompactContigIntervals],
    store: CompactHeaderStore,
    reference_haplotypes: FrozenSet[str],
    jobs: int,
) -> Iterator[SweepResult]:
    keys = sorted(index)
    _initialize_compact_sweep_workers(index, store, reference_haplotypes)
    if jobs == 1:
        yield from map(sweep_compact_contig, keys)
        return
    with ProcessPoolExecutor(
        max_workers=jobs,
        mp_context=mp_context(),
        initializer=_initialize_compact_sweep_workers,
        initargs=(index, store, reference_haplotypes),
    ) as executor:
        key_iterator = iter(keys)
        pending = deque()
        for _index in range(jobs * 2):
            try:
                key = next(key_iterator)
            except StopIteration:
                break
            pending.append(executor.submit(sweep_compact_contig, key))
        while pending:
            yield pending.popleft().result()
            try:
                key = next(key_iterator)
            except StopIteration:
                continue
            pending.append(executor.submit(sweep_compact_contig, key))


def sweep_all_header_file_contigs(
    index: Dict[ContigKey, HeaderFileContigIntervals],
    partition_names: Sequence[str],
    reference_haplotypes: FrozenSet[str],
    jobs: int,
) -> Iterator[SweepResult]:
    keys = sorted(index)
    _initialize_header_file_sweep_workers(
        index, partition_names, reference_haplotypes,
    )
    if jobs == 1:
        yield from map(sweep_header_file_contig, keys)
        return
    with ProcessPoolExecutor(
        max_workers=jobs,
        mp_context=mp_context(),
        initializer=_initialize_header_file_sweep_workers,
        initargs=(index, partition_names, reference_haplotypes),
    ) as executor:
        key_iterator = iter(keys)
        pending = deque()
        for _index in range(jobs * 2):
            try:
                key = next(key_iterator)
            except StopIteration:
                break
            pending.append(executor.submit(sweep_header_file_contig, key))
        while pending:
            yield pending.popleft().result()
            try:
                key = next(key_iterator)
            except StopIteration:
                continue
            pending.append(executor.submit(sweep_header_file_contig, key))


def _csv_list(values: Iterable[str]) -> str:
    ordered = sorted(values)
    return ",".join(ordered) if ordered else "."


def _open_tsv(path: str, fields: Sequence[str], bed: bool = False):
    Path(path).parent.mkdir(parents=True, exist_ok=True)
    temporary = path + ".tmp"
    handle = open(temporary, "wt", newline="")
    writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
    header = list(fields)
    if bed:
        header[0] = "#" + header[0]
    writer.writerow(header)
    return handle, writer, temporary


def export_input_catalogs(
    catalogs: Sequence[CoordinateCatalog],
    output_dir: str,
) -> None:
    write_tsv(
        os.path.join(output_dir, "partition_inputs.tsv"),
        ("partition", "fasta", "unique_bed", "reference_bed"),
        (
            (
                catalog.files.partition,
                os.path.abspath(catalog.files.fasta),
                os.path.abspath(catalog.files.unique_bed),
                os.path.abspath(catalog.files.reference_bed),
            )
            for catalog in catalogs
        ),
    )
    write_bed_table(
        os.path.join(output_dir, "all_header_regions.bed"),
        (
            "source_contig", "source_start", "source_end", "record_id",
            "score", "strand", "haplotype", "partition", "record_length",
        ),
        (
            (
                row.source_contig, row.source_start, row.source_end,
                row.record_id, 0, "+", row.haplotype, row.partition,
                row.record_length,
            )
            for catalog in catalogs for row in catalog.headers
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
                row.source_contig, row.source_start, row.source_end,
                row.region_name, 0, "+", row.haplotype, row.partition,
                row.record_id, row.region_class, row.local_start, row.local_end,
                row.raw_local_start, row.raw_local_end, row.novelty_score,
            )
            for catalog in catalogs for row in catalog.unique_regions
        ),
    )
    write_bed_table(
        os.path.join(output_dir, "all_used_novel_loci.bed"),
        (
            "source_contig", "source_start", "source_end", "encoded_name",
            "score", "strand", "haplotype", "partition", "local_start",
            "local_end", "block_id", "block_score",
        ),
        (
            (
                row.source_contig, row.source_start, row.source_end,
                row.encoded_name, 0, "+", row.haplotype, row.partition,
                row.local_start, row.local_end, row.block_id, row.block_score,
            )
            for catalog in catalogs for row in catalog.used_novel_regions
        ),
    )

    def mapped_rows() -> Iterator[Tuple]:
        for catalog in catalogs:
            by_record: DefaultDict[str, List[Interval]] = defaultdict(list)
            for unique in catalog.unique_regions:
                by_record[unique.record_id].append((
                    unique.local_start, unique.local_end,
                ))
            for header in catalog.headers:
                for start, end in complement_intervals(
                    by_record.get(header.record_id, ()),
                    0, header.record_length,
                ):
                    yield (
                        header.source_contig,
                        header.source_start + start,
                        header.source_start + end,
                        header.record_id,
                        0,
                        "+",
                        header.haplotype,
                        header.partition,
                    )

    write_bed_table(
        os.path.join(output_dir, "all_mapped_regions.bed"),
        (
            "source_contig", "source_start", "source_end", "record_id",
            "score", "strand", "haplotype", "partition",
        ),
        mapped_rows(),
    )


def write_compact_header_catalog(
    store: CompactHeaderStore,
    output_dir: str,
) -> None:
    path = os.path.join(output_dir, "all_header_regions.bed")
    temporary = path + ".tmp"
    index_rows = []
    header = (
        "#source_contig\tsource_start\tsource_end\trecord_id\tscore\t"
        "strand\thaplotype\tpartition\trecord_length\n"
    ).encode()
    try:
        with open(temporary, "wb") as out:
            out.write(header)
            progress_step = max(1, len(store.partitions) // 20)
            for partition_number, partition in enumerate(
                store.partitions, 1,
            ):
                byte_start = out.tell()
                record_count = 0
                for _holder_id, holder, local_index in (
                    store.iter_partition_entries(partition)
                ):
                    start = int(holder.source_start[local_index])
                    end = int(holder.source_end[local_index])
                    values = (
                        holder.contigs[holder.contig_id[local_index]],
                        start,
                        end,
                        (
                            f"g{partition}_{holder.assembly.name}_"
                            f"{holder.sequence_index[local_index]}"
                        ),
                        0,
                        "+",
                        holder.assembly.name,
                        partition,
                        end - start,
                    )
                    out.write(("\t".join(map(str, values)) + "\n").encode())
                    record_count += 1
                index_rows.append((
                    partition, byte_start, out.tell(), record_count,
                ))
                if (
                    partition_number == len(store.partitions)
                    or partition_number % progress_step == 0
                ):
                    LOG.info(
                        "Wrote reconstructed header catalog for %d/%d "
                        "partitions",
                        partition_number, len(store.partitions),
                    )
        os.replace(temporary, path)
    finally:
        try:
            os.remove(temporary)
        except FileNotFoundError:
            pass
    write_tsv(
        os.path.join(output_dir, "all_header_regions.index.tsv"),
        ("partition", "byte_start", "byte_end", "record_count"),
        index_rows,
    )


def export_compact_input_catalogs(
    partitions: Sequence[PartitionFiles],
    store: CompactHeaderStore,
    catalogs: Sequence[AnnotationCatalog],
    output_dir: str,
) -> None:
    write_tsv(
        os.path.join(output_dir, "partition_inputs.tsv"),
        ("partition", "fasta", "unique_bed", "reference_bed"),
        (
            (
                partition.partition,
                os.path.abspath(partition.fasta),
                os.path.abspath(partition.unique_bed),
                os.path.abspath(partition.reference_bed),
            )
            for partition in partitions
        ),
    )
    write_compact_header_catalog(store, output_dir)
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
                row.source_contig, row.source_start, row.source_end,
                row.region_name, 0, "+", row.haplotype, row.partition,
                row.record_id, row.region_class, row.local_start, row.local_end,
                row.raw_local_start, row.raw_local_end, row.novelty_score,
            )
            for catalog in catalogs for row in catalog.unique_regions
        ),
    )
    write_bed_table(
        os.path.join(output_dir, "all_used_novel_loci.bed"),
        (
            "source_contig", "source_start", "source_end", "encoded_name",
            "score", "strand", "haplotype", "partition", "local_start",
            "local_end", "block_id", "block_score",
        ),
        (
            (
                row.source_contig, row.source_start, row.source_end,
                row.encoded_name, 0, "+", row.haplotype, row.partition,
                row.local_start, row.local_end, row.block_id, row.block_score,
            )
            for catalog in catalogs for row in catalog.used_novel_regions
        ),
    )
    unique_by_header: DefaultDict[int, List[Interval]] = defaultdict(list)
    for catalog in catalogs:
        for row in catalog.unique_regions:
            unique_by_header[row.header_index].append((
                row.local_start, row.local_end,
            ))

    def mapped_rows() -> Iterator[Tuple]:
        progress_step = max(1, len(store.partitions) // 20)
        for partition_number, partition in enumerate(store.partitions, 1):
            for holder_id, holder, local_index in (
                store.iter_partition_entries(partition)
            ):
                header_index = store.locator(holder_id, local_index)
                source_start = int(holder.source_start[local_index])
                record_length = int(
                    holder.source_end[local_index] - source_start
                )
                for start, end in complement_intervals(
                    unique_by_header.get(header_index, ()),
                    0,
                    record_length,
                ):
                    yield (
                        holder.contigs[holder.contig_id[local_index]],
                        source_start + start,
                        source_start + end,
                        (
                            f"g{partition}_{holder.assembly.name}_"
                            f"{holder.sequence_index[local_index]}"
                        ),
                        0,
                        "+",
                        holder.assembly.name,
                        partition,
                    )
            if (
                partition_number == len(store.partitions)
                or partition_number % progress_step == 0
            ):
                LOG.info(
                    "Wrote mapped-region catalog for %d/%d partitions",
                    partition_number, len(store.partitions),
                )

    write_bed_table(
        os.path.join(output_dir, "all_mapped_regions.bed"),
        (
            "source_contig", "source_start", "source_end", "record_id",
            "score", "strand", "haplotype", "partition",
        ),
        mapped_rows(),
    )


def export_header_file_input_catalogs(
    partitions: Sequence[PartitionFiles],
    header_locations: Dict[str, Tuple[str, int]],
    catalogs: Sequence[AnnotationCatalog],
    output_dir: str,
) -> None:
    """Save manifests/annotations while retaining .header as source of truth."""
    for obsolete in (
        "all_header_regions.bed",
        "all_header_regions.index.tsv",
        "all_mapped_regions.bed",
    ):
        try:
            os.remove(os.path.join(output_dir, obsolete))
        except FileNotFoundError:
            pass
    write_tsv(
        os.path.join(output_dir, "partition_inputs.tsv"),
        ("partition", "fasta", "header", "unique_bed", "reference_bed"),
        (
            (
                partition.partition,
                os.path.abspath(partition.fasta),
                header_locations[partition.partition][0],
                os.path.abspath(partition.unique_bed),
                os.path.abspath(partition.reference_bed),
            )
            for partition in partitions
        ),
    )
    write_tsv(
        os.path.join(output_dir, "header_inputs.tsv"),
        ("partition", "header", "record_count"),
        (
            (
                partition.partition,
                header_locations[partition.partition][0],
                header_locations[partition.partition][1],
            )
            for partition in partitions
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
                row.source_contig, row.source_start, row.source_end,
                row.region_name, 0, "+", row.haplotype, row.partition,
                row.record_id, row.region_class, row.local_start, row.local_end,
                row.raw_local_start, row.raw_local_end, row.novelty_score,
            )
            for catalog in catalogs for row in catalog.unique_regions
        ),
    )
    write_bed_table(
        os.path.join(output_dir, "all_used_novel_loci.bed"),
        (
            "source_contig", "source_start", "source_end", "encoded_name",
            "score", "strand", "haplotype", "partition", "local_start",
            "local_end", "block_id", "block_score",
        ),
        (
            (
                row.source_contig, row.source_start, row.source_end,
                row.encoded_name, 0, "+", row.haplotype, row.partition,
                row.local_start, row.local_end, row.block_id, row.block_score,
            )
            for catalog in catalogs for row in catalog.used_novel_regions
        ),
    )


def export_sweeps(
    sweeps: Iterable[SweepResult],
    contig_index,
    catalogs,
    output_dir: str,
) -> Tuple[int, int]:
    fasta_by_partition = {
        catalog.files.partition: os.path.abspath(catalog.files.fasta)
        for catalog in catalogs
    }
    status_fields = (
        "source_contig", "source_start", "source_end", "interval_id",
        "haplotype", "status", "header_partitions", "lifted_partitions",
        "unlifted_partitions", "used_novel_locus",
    )
    category_fields = (
        "source_contig", "source_start", "source_end", "interval_id",
        "score", "strand", "haplotype", "header_partitions",
        "lifted_partitions", "unlifted_partitions", "used_novel_locus",
    )
    candidate_fields = (
        "candidate_id", "interval_id", "source_contig", "source_start",
        "source_end", "haplotype", "novel_partition", "lifted_partitions",
        "record_id", "local_start", "local_end", "source_fasta",
        "region_class",
    )
    opened = []
    try:
        status_handle, status_writer, status_tmp = _open_tsv(
            os.path.join(output_dir, "interval_status.tsv"), status_fields,
        )
        opened.append((
            status_handle, status_tmp,
            os.path.join(output_dir, "interval_status.tsv"),
        ))
        candidate_handle, candidate_writer, candidate_tmp = _open_tsv(
            os.path.join(output_dir, "local_novel_candidates.tsv"),
            candidate_fields,
        )
        opened.append((
            candidate_handle, candidate_tmp,
            os.path.join(output_dir, "local_novel_candidates.tsv"),
        ))
        category_writers = {}
        for status in ("mapped", "unmapped", "local_novel", "global_novel"):
            handle, writer, temporary = _open_tsv(
                os.path.join(output_dir, f"{status}.bed"),
                category_fields,
                bed=True,
            )
            opened.append((
                handle, temporary,
                os.path.join(output_dir, f"{status}.bed"),
            ))
            category_writers[status] = writer

        interval_number = 0
        candidate_number = 0
        progress = 0
        for result in sweeps:
            progress += 1
            interval_ids: List[str] = []
            for row in result.statuses:
                interval_number += 1
                interval_id = f"interval_{interval_number:012d}"
                interval_ids.append(interval_id)
                values = (
                    result.key[1], row.start, row.end, interval_id,
                    result.key[0], row.status,
                    _csv_list(row.header_partitions),
                    _csv_list(row.lifted_partitions),
                    _csv_list(row.unlifted_partitions),
                    _csv_list(row.used_novel_locus),
                )
                status_writer.writerow(values)
                category_writers[row.status].writerow((
                    result.key[1], row.start, row.end, interval_id,
                    0, "+", result.key[0],
                    _csv_list(row.header_partitions),
                    _csv_list(row.lifted_partitions),
                    _csv_list(row.unlifted_partitions),
                    _csv_list(row.used_novel_locus),
                ))
            data = contig_index[result.key]
            for row in result.candidates:
                candidate_number += 1
                candidate_id = f"candidate_{candidate_number:012d}"
                unique = data.unique_regions[row.unique_index]
                local_start = unique.local_start + row.start - unique.source_start
                local_end = local_start + row.end - row.start
                candidate_writer.writerow((
                    candidate_id,
                    interval_ids[row.status_index],
                    result.key[1],
                    row.start,
                    row.end,
                    result.key[0],
                    row.novel_partition,
                    _csv_list(row.lifted_partitions),
                    unique.record_id,
                    local_start,
                    local_end,
                    fasta_by_partition[row.novel_partition],
                    unique.region_class,
                ))
            if progress % 1000 == 0:
                LOG.info("Swept %d/%d haplotype-contigs", progress, len(contig_index))
        for handle, temporary, final in opened:
            handle.close()
            os.replace(temporary, final)
        opened.clear()
        return interval_number, candidate_number
    finally:
        for handle, temporary, _final in opened:
            handle.close()
            try:
                os.remove(temporary)
            except FileNotFoundError:
                pass


def input_signature(
    partitions: Sequence[PartitionFiles],
    anchor: int,
    reference_haplotypes: FrozenSet[str],
    extra_fingerprints: Iterable[Tuple[str, int, int]] = (),
    reconstruction_settings: Optional[Dict[str, int]] = None,
    ignored_source_contigs: FrozenSet[str] = IGNORED_SOURCE_CONTIGS,
    allow_missing_reference_haplotypes: bool = False,
    allow_missing_bed_haplotypes: bool = False,
) -> str:
    digest = hashlib.sha256()
    digest.update(
        (
            "coordinate-analysis-v6\n"
            f"anchor={anchor}\n"
            f"reference={','.join(sorted(reference_haplotypes))}\n"
            f"ignored_contigs={','.join(sorted(ignored_source_contigs))}\n"
            "allow_missing_reference_haplotypes="
            f"{int(allow_missing_reference_haplotypes)}\n"
            "allow_missing_bed_haplotypes="
            f"{int(allow_missing_bed_haplotypes)}\n"
        ).encode()
    )
    if reconstruction_settings:
        digest.update(json.dumps(
            reconstruction_settings, sort_keys=True,
        ).encode())
        digest.update(b"\n")
    for partition in partitions:
        digest.update((partition.partition + "\n").encode())
        for path, size, mtime_ns in partition.fingerprints:
            digest.update(f"{path}\t{size}\t{mtime_ns}\n".encode())
    for path, size, mtime_ns in extra_fingerprints:
        digest.update(f"{path}\t{size}\t{mtime_ns}\n".encode())
    return digest.hexdigest()


def parse_args(argv: Optional[Sequence[str]] = None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=(
            "Load all local-block coordinate annotations into RAM, classify "
            "atomic source intervals, and save local-novel mapping inputs."
        )
    )
    parser.add_argument(
        "-i", "--localblocks", required=True,
        help="root containing PARTITION/PARTITION.{fasta,bed,unique.bed}",
    )
    parser.add_argument("-o", "--output-dir", required=True)
    parser.add_argument(
        "-j", "--jobs", type=int, default=16,
        help=(
            "parallel .header/annotation readers; also hotspot processes "
            "with --use-hotspots (default: 16)"
        ),
    )
    parser.add_argument(
        "--sweep-jobs", type=int,
        help=(
            "per-contig sorting processes; defaults to min(--jobs,16) to "
            "limit copy-on-write/process memory"
        ),
    )
    parser.add_argument(
        "-q", "--query-paths",
        help=(
            "original NAME FASTA (with adjacent FASTA.fai) list; used only with --use-hotspots"
        ),
    )
    parser.add_argument(
        "--use-hotspots", action="store_true",
        help=(
            "legacy alternative: reconstruct headers from hotspot files; "
            "default mode strictly requires PARTITION.header files"
        ),
    )
    parser.add_argument(
        "--kmermaps",
        help=(
            "directory containing NAME_hotspot.txt; with --query-paths, "
            "defaults to the kmermaps sibling of --localblocks"
        ),
    )
    parser.add_argument(
        "--blocks-bed",
        help=(
            "original scored block BED; with --query-paths, defaults to "
            "main_chroms_withnovel_blocks_kmerscore.bed beside localblocks"
        ),
    )
    parser.add_argument(
        "--target-list",
        help=(
            "KmerSearcher -T <GraphFolder>.list; remaps hotspot column 2 "
            "from compact target-list row IDs to original partition names"
        ),
    )
    parser.add_argument("--hotspot-kmer-length", type=int, default=31)
    parser.add_argument("--hotspot-merge-distance", type=int, default=30_000)
    parser.add_argument("--hotspot-anchor", type=int, default=15_000)
    parser.add_argument(
        "--hotspot-reference-extension", type=int, default=5_000,
    )
    parser.add_argument(
        "--input-anchor", type=int, default=50,
        help="reporting anchor to remove from unique BED sides (default: 50)",
    )
    parser.add_argument(
        "--reference-haplotypes", default="CHM13_h1",
        help=(
            "comma-separated true reference haplotypes; these receive only "
            "mapped/unmapped statuses (default: CHM13_h1)"
        ),
    )
    parser.add_argument(
        "--allow-missing-reference-haplotypes", action="store_true",
        help=(
            "allow an isolated partition whose extracted headers contain no "
            "true-reference record; retain the configured reference names in "
            "metadata and treat the partition's block BED as its unchanged "
            "original template"
        ),
    )
    parser.add_argument(
        "--allow-missing-bed-haplotypes", action="store_true",
        help=(
            "allow PARTITION.bed source haplotypes supplied by an external "
            "alternative-template FASTA even when those haplotypes are not "
            "members of the cohort query list"
        ),
    )
    parser.add_argument(
        "--ignore-contigs", default="chrM",
        help=(
            "comma-separated source contigs excluded from headers and all "
            "annotations (default: chrM); pass an empty string to disable"
        ),
    )
    parser.add_argument("--resume", action="store_true")
    parser.add_argument(
        "--log-level", choices=("DEBUG", "INFO", "WARNING", "ERROR"),
        default="INFO",
    )
    args = parser.parse_args(argv)
    if args.jobs < 1:
        parser.error("--jobs must be positive")
    if args.sweep_jobs is not None and args.sweep_jobs < 1:
        parser.error("--sweep-jobs must be positive")
    if args.input_anchor < 0:
        parser.error("--input-anchor cannot be negative")
    for name in (
        "hotspot_kmer_length", "hotspot_merge_distance",
        "hotspot_anchor", "hotspot_reference_extension",
    ):
        if getattr(args, name) < 0:
            parser.error(f"--{name.replace('_', '-')} cannot be negative")
    if args.hotspot_kmer_length == 0:
        parser.error("--hotspot-kmer-length must be positive")
    if args.use_hotspots and not args.query_paths:
        parser.error("--use-hotspots requires --query-paths")
    if (args.kmermaps or args.blocks_bed or args.target_list) and not args.use_hotspots:
        parser.error("--kmermaps/--blocks-bed/--target-list require --use-hotspots")
    return args


def run(args: argparse.Namespace) -> None:
    output_dir = os.path.abspath(args.output_dir)
    Path(output_dir).mkdir(parents=True, exist_ok=True)
    reference_haplotypes = parse_reference_haplotypes(
        args.reference_haplotypes
    )
    ignored_source_contigs = frozenset(
        token.strip()
        for token in args.ignore_contigs.split(",")
        if token.strip()
    )
    partitions = discover_partitions(args.localblocks, jobs=args.jobs)
    sweep_jobs = args.sweep_jobs or min(args.jobs, 16)
    hotspot_mode = bool(args.use_hotspots)
    sources = None
    blocks = None
    header_locations = None
    extra_fingerprints = []
    reconstruction_settings = None
    if hotspot_mode:
        cohort_root = os.path.dirname(os.path.abspath(args.localblocks))
        kmermaps = os.path.abspath(
            args.kmermaps or os.path.join(cohort_root, "kmermaps")
        )
        blocks_bed = os.path.abspath(
            args.blocks_bed
            or os.path.join(
                cohort_root,
                "main_chroms_withnovel_blocks_kmerscore.bed",
            )
        )
        sources = read_query_sources(
            os.path.abspath(args.query_paths), kmermaps,
        )
        blocks = read_blocks(blocks_bed)
        target_list = (
            os.path.abspath(args.target_list)
            if args.target_list else None
        )
        if target_list is not None:
            if not os.path.isfile(target_list):
                raise FileNotFoundError(target_list)
            blocks = remap_blocks_to_target_list(blocks, target_list)
        partition_names = {row.partition for row in partitions}
        missing_blocks = partition_names - {
            block.name for block in blocks.values()
        }
        if missing_blocks:
            example = ",".join(sorted(missing_blocks)[:10])
            raise ValueError(
                f"{len(missing_blocks)} local-block partitions are absent "
                f"from {blocks_bed}; examples: {example}"
            )
        extra_fingerprints = list(source_fingerprints(
            os.path.abspath(args.query_paths),
            blocks_bed,
            sources,
            (target_list,) if target_list is not None else (),
        ))
        reconstruction_settings = {
            "kmer_length": args.hotspot_kmer_length,
            "merge_distance": args.hotspot_merge_distance,
            "anchor": args.hotspot_anchor,
            "reference_extension": args.hotspot_reference_extension,
        }
    else:
        header_locations, header_fingerprints = validate_extracted_headers(
            partitions, args.jobs,
        )
        extra_fingerprints = header_fingerprints
    signature = input_signature(
        partitions,
        args.input_anchor,
        reference_haplotypes,
        extra_fingerprints,
        reconstruction_settings,
        ignored_source_contigs,
        args.allow_missing_reference_haplotypes,
        args.allow_missing_bed_haplotypes,
    )
    marker = os.path.join(output_dir, "_COORDINATE_ANALYSIS_SUCCESS.json")
    if args.resume and os.path.isfile(marker):
        try:
            with open(marker, "rt") as handle:
                saved = json.load(handle)
            required = (
                "interval_status.tsv",
                "local_novel_candidates.tsv",
                "partition_inputs.tsv",
                "global_novel.bed",
            )
            if not hotspot_mode:
                required += ("header_inputs.tsv",)
            if (
                saved.get("signature") == signature
                and all(os.path.isfile(os.path.join(output_dir, name))
                        for name in required)
            ):
                LOG.info("RESUME: coordinate analysis is complete in %s", output_dir)
                return
        except (OSError, ValueError, AttributeError):
            pass

    if hotspot_mode:
        assert sources is not None and blocks is not None
        LOG.info(
            "Reconstructing independent per-assembly header holders from %d "
            "hotspot maps with %d processes (local FASTAs will not be scanned)",
            len(sources), args.jobs,
        )
        progress_step = max(1, len(sources) // 20)

        def header_progress(done, total, headers, hotspots):
            if done == total or done % progress_step == 0:
                LOG.info(
                    "Reconstructed %d/%d assemblies: %d headers from "
                    "%d hotspots",
                    done, total, headers, hotspots,
                )

        store, hotspot_count = reconstruct_headers(
            sources,
            blocks,
            [row.partition for row in partitions],
            args.jobs,
            args.hotspot_kmer_length,
            args.hotspot_merge_distance,
            args.hotspot_anchor,
            args.hotspot_reference_extension,
            header_progress,
        )
        LOG.info(
            "Loading unique/BED annotations for %d partitions",
            len(partitions),
        )
        catalogs = load_compact_annotations(
            partitions,
            store,
            args.input_anchor,
            args.jobs,
            ignored_source_contigs,
        )
        LOG.info(
            "Compact header storage before catalog export: %.2f GiB "
            "(numeric columns/indexes only)",
            store.packed_size_bytes() / 1024**3,
        )
        observed_haplotypes = set(store.haplotypes)
    else:
        LOG.info(
            "Loading only extracted .header and annotation files for %d "
            "partitions with %d workers (FASTA sequence bytes are not read)",
            len(partitions), args.jobs,
        )
        assert header_locations is not None
        (
            catalogs,
            contig_index,
            observed_haplotypes,
            header_count,
        ) = load_header_file_catalogs(
            partitions,
            header_locations,
            args.input_anchor,
            args.jobs,
            ignored_source_contigs,
        )
        store = None
        hotspot_count = 0
    bed_haplotypes = {
        row.haplotype
        for catalog in catalogs
        for row in catalog.used_novel_regions
    }
    unknown_bed_haplotypes = bed_haplotypes - observed_haplotypes
    if unknown_bed_haplotypes:
        message = (
            "coordinate-named PARTITION.bed haplotypes absent from all cohort "
            "FASTA headers: " + ",".join(sorted(unknown_bed_haplotypes))
        )
        if not args.allow_missing_bed_haplotypes:
            from fixed_alternatives import unbacked_template_haplotypes
            unbacked = unbacked_template_haplotypes(
                {row.partition: row for row in partitions},
                (row for catalog in catalogs for row in catalog.used_novel_regions),
                observed_haplotypes,
            )
            if unbacked:
                raise ValueError(message)
        LOG.warning(
            "%s; accepting them as external alternative templates and "
            "retaining their original block coordinates",
            message,
        )
    missing_reference = reference_haplotypes - observed_haplotypes
    if missing_reference:
        message = (
            "--reference-haplotypes absent from all FASTA headers: "
            + ",".join(sorted(missing_reference))
        )
        if not args.allow_missing_reference_haplotypes:
            raise ValueError(message)
        LOG.warning(
            "%s; continuing in isolated-partition mode with the original "
            "block template unchanged",
            message,
        )
    LOG.info("Building in-memory interval index")
    if hotspot_mode:
        assert store is not None
        contig_index = build_compact_contig_index(
            store, catalogs, ignored_source_contigs,
        )
        export_compact_input_catalogs(
            partitions, store, catalogs, output_dir,
        )
        packed_before_release = store.packed_size_bytes()
        store.release_for_sweep()
        LOG.info(
            "Released header catalog-only columns: %.2f GiB -> %.2f GiB "
            "of packed coordinate storage before sorting",
            packed_before_release / 1024**3,
            store.packed_size_bytes() / 1024**3,
        )
        sweeps = sweep_all_compact_contigs(
            contig_index, store, reference_haplotypes, sweep_jobs,
        )
    else:
        assert header_locations is not None
        export_header_file_input_catalogs(
            partitions, header_locations, catalogs, output_dir,
        )
        sweeps = sweep_all_header_file_contigs(
            contig_index,
            [partition.partition for partition in partitions],
            reference_haplotypes,
            sweep_jobs,
        )
    LOG.info(
        "Sweeping %d haplotype-contigs with %d workers",
        len(contig_index), sweep_jobs,
    )
    intervals, candidates = export_sweeps(
        sweeps,
        contig_index,
        catalogs,
        output_dir,
    )
    temporary = marker + ".tmp"
    with open(temporary, "wt") as out:
        json.dump(
            {
                "signature": signature,
                "partitions": len(partitions),
                "haplotype_contigs": len(contig_index),
                "intervals": intervals,
                "local_novel_candidates": candidates,
                "reference_haplotypes": sorted(reference_haplotypes),
                "missing_reference_haplotypes": sorted(missing_reference),
                "missing_bed_haplotypes": sorted(unknown_bed_haplotypes),
                "ignored_source_contigs": sorted(ignored_source_contigs),
                "input_anchor": args.input_anchor,
                "header_source": "hotspots" if hotspot_mode else "header_files",
                "headers": len(store) if store is not None else header_count,
                "hotspots": hotspot_count,
                "reconstruction_settings": reconstruction_settings,
            },
            out,
            sort_keys=True,
        )
        out.write("\n")
    os.replace(temporary, marker)
    LOG.info(
        "Wrote %d interval statuses and %d local-novel candidates under %s",
        intervals, candidates, output_dir,
    )


def main(argv: Optional[Sequence[str]] = None) -> int:
    args = parse_args(argv)
    logging.basicConfig(
        level=getattr(logging, args.log_level),
        format="[%(asctime)s] %(levelname)s: %(message)s",
        datefmt="%Y-%m-%d %H:%M:%S",
    )
    try:
        run(args)
    except (OSError, ValueError, RuntimeError) as error:
        LOG.error("%s", error)
        return 1
    return 0


if __name__ == "__main__":
    sys.exit(main())
