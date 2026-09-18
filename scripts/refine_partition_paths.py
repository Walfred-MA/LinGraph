#!/usr/bin/env python3
"""Build adjusted partition paths from completed local/global mappings.

Global-novel alignment is performed separately by ``map_global_novel.py``.
This script reads those completed mappings together with local-novel mappings.
A mapped local or global interval is replaced by cleaned-reference sequence
only when its reference placement belongs to the merged component containing
the original block template. Distant placements remain query alternatives
carrying graphic liftover annotations.
"""
from __future__ import annotations

import argparse
import bisect
import csv
import dataclasses
import hashlib
import json
import logging
import os
import re
import shutil
import sys
from collections import defaultdict
from concurrent.futures import (
    FIRST_COMPLETED,
    ProcessPoolExecutor,
    ThreadPoolExecutor,
    wait,
)
from contextlib import ExitStack
from functools import lru_cache
from pathlib import Path
from typing import (
    DefaultDict,
    Dict,
    FrozenSet,
    Iterable,
    Iterator,
    List,
    Mapping,
    Optional,
    Sequence,
    Set,
    Tuple,
)

from map_partition_local_novel import (
    AssemblyFastaSource,
    BlockTarget,
    CachedCandidateSequence,
    MappingCandidate,
    MappedPlacement,
    PreparedBlockReference,
    REPORTING_RESOLUTION_VERSION,
    load_candidate_header_records,
    load_candidate_sequences,
    map_partition_candidates,
    placement_sort_key,
    read_assembly_sources,
    read_candidates,
    read_partition_inputs,
    read_placements,
    read_reporting_components,
    write_reporting_components,
)
from minsetref_align import AlignmentHit, ensure_executable
from minsetref_core import IndexedFasta, mp_context
from minsetref_graphic import (
    build_graphic_path,
    build_multilocus_graphic_path,
)
from minsetref_multilocus import ResolvedPlacementComponent
from minsetref_light import (
    MASKED_BASE_WEIGHT,
)
from minsetref_segments import (
    count_unmasked,
    merge_intervals,
    novelty_score,
    subtract_intervals,
)
from summarize_partition_novelty import (
    HeaderRecord,
    PartitionFiles,
    clip_hit_to_query_core,
    fasta_header_descriptions,
    parse_header_records,
    write_tsv,
)


LOG = logging.getLogger("refine_partition_paths")
Interval = Tuple[int, int]
SOURCE_COORD_RE = re.compile(r"^(.*):([^:]+):(\d+)-(\d+)([+-]?)$")
_FORKED_PARTITION_BUILDER = None


def _run_forked_partition_builder(partition: str):
    """Run the inherited per-partition closure in a Linux fork worker."""
    if _FORKED_PARTITION_BUILDER is None:
        raise RuntimeError("forked partition builder is not initialized")
    return _FORKED_PARTITION_BUILDER(partition)


def _bounded_thread_map(function, items: Iterable, workers: int):
    """Yield completed results while retaining at most ``workers`` futures."""
    iterator = iter(items)
    with ThreadPoolExecutor(max_workers=workers) as executor:
        pending = set()
        for _index in range(workers):
            try:
                item = next(iterator)
            except StopIteration:
                break
            pending.add(executor.submit(function, item))
        while pending:
            done, pending = wait(pending, return_when=FIRST_COMPLETED)
            for future in done:
                yield future.result()
                try:
                    item = next(iterator)
                except StopIteration:
                    continue
                pending.add(executor.submit(function, item))


def _bounded_process_map(
    function,
    items: Sequence,
    workers: int,
) -> Iterator:
    """Yield process results with bounded task/result memory."""
    if not items:
        return
    if workers <= 1:
        for item in items:
            yield function(item)
        return
    iterator = iter(items)
    worker_count = min(workers, len(items))
    with ProcessPoolExecutor(
        max_workers=worker_count,
        mp_context=mp_context(),
    ) as executor:
        pending = {
            executor.submit(function, next(iterator))
            for _index in range(worker_count)
        }
        while pending:
            done, pending = wait(pending, return_when=FIRST_COMPLETED)
            for future in done:
                yield future.result()
                try:
                    item = next(iterator)
                except StopIteration:
                    continue
                pending.add(executor.submit(function, item))


@dataclasses.dataclass(frozen=True)
class UniqueRow:
    partition: str
    record_id: str
    haplotype: str
    source_contig: str
    source_start: int
    source_end: int
    local_start: int
    local_end: int
    region_name: str
    region_class: str


@dataclasses.dataclass(frozen=True)
class StatusRow:
    haplotype: str
    source_contig: str
    source_start: int
    source_end: int
    interval_id: str
    status: str
    unlifted_partitions: FrozenSet[str]


@dataclasses.dataclass(frozen=True)
class UsedReferenceRow:
    partition: str
    haplotype: str
    source_contig: str
    source_start: int
    source_end: int
    encoded_name: str


@dataclasses.dataclass(frozen=True)
class NovelQuery:
    query_id: str
    origin: str
    partition: str
    interval_id: str
    source_fasta: str
    record_id: str
    haplotype: str
    source_contig: str
    source_start: int
    source_end: int
    local_start: int
    local_end: int


@dataclasses.dataclass(frozen=True)
class QuerySequence:
    query: NovelQuery
    extract_start: int
    extract_end: int
    core_start: int
    core_end: int
    sequence: str


@dataclasses.dataclass(frozen=True)
class ReferenceRecord:
    record_id: str
    haplotype: str
    source_contig: str
    source_start: int
    source_end: int
    strand: str
    length: int


@dataclasses.dataclass(frozen=True)
class GlobalPlacement:
    query_id: str
    partition: str
    ref_haplotype: str
    ref_contig: str
    ref_start: int
    ref_end: int
    strand: str
    alignment_score: float
    aligned_score: float
    aligned_bases: int
    identity: float
    source: str
    target_record: str
    query_blocks: Tuple[Interval, ...]
    target_blocks: Tuple[Interval, ...]
    paired_blocks: Tuple[Tuple[int, int, int, int], ...] = ()
    liftover_path: str = "."


@dataclasses.dataclass
class ReferenceSpan:
    haplotype: str
    contig: str
    start: int
    end: int
    score: float
    origins: Set[str]


@dataclasses.dataclass(frozen=True)
class PathRow:
    partition: str
    path_order: int
    path_id: str
    role: str
    label: str
    source_kind: str
    source_haplotype: str
    source_contig: str
    source_start: int
    source_end: int
    strand: str
    segments: str
    origins: str


def compact_segments(segments: Sequence[Mapping[str, object]]) -> str:
    """Encode materialized local segments without embedding file paths/JSON."""
    if not segments:
        return "."
    encoded: List[str] = []
    for segment in segments:
        record = str(segment["record"])
        start, end = int(segment["start"]), int(segment["end"])
        strand = str(segment["strand"])
        if start < 0 or end <= start or strand not in {"+", "-"}:
            raise ValueError(
                f"invalid compact segment {record}:{start}-{end}{strand}"
            )
        encoded.append(f"{record}_{start}_{end}{strand}")
    return ";".join(encoded)


def _parse_set(text: str) -> FrozenSet[str]:
    if not text or text == ".":
        return frozenset()
    return frozenset(value for value in text.split(",") if value and value != ".")


def _read_dict_rows(path: str) -> Iterable[Dict[str, str]]:
    with open(path, "rt", newline="") as handle:
        first = handle.readline()
        if not first:
            return
        if first.startswith("#"):
            first = first[1:]
        fieldnames = next(csv.reader([first], delimiter="\t"))
        yield from csv.DictReader(
            handle, delimiter="\t", fieldnames=fieldnames,
        )


def _fingerprint(path: str) -> Tuple[str, int, int]:
    stat = os.stat(path)
    return os.path.abspath(path), stat.st_size, stat.st_mtime_ns


def read_novel_queries(path: str) -> List[NovelQuery]:
    output: List[NovelQuery] = []
    seen: Set[str] = set()
    for row in _read_dict_rows(path):
        query_id = row["query_id"]
        if query_id in seen:
            raise ValueError(f"{path}: duplicate query_id {query_id!r}")
        seen.add(query_id)
        output.append(NovelQuery(
            query_id,
            row["origin"],
            row["partition"],
            row["interval_id"],
            row["source_fasta"],
            row["record_id"],
            row["haplotype"],
            row["source_contig"],
            int(row["source_start"]),
            int(row["source_end"]),
            int(row["local_start"]),
            int(row["local_end"]),
        ))
    return output


def read_global_placements(path: str) -> Dict[str, GlobalPlacement]:
    output: Dict[str, GlobalPlacement] = {}
    for row in _read_dict_rows(path):
        query_id = row["query_id"]
        if query_id in output:
            raise ValueError(f"{path}: duplicate query_id {query_id!r}")
        output[query_id] = GlobalPlacement(
            query_id,
            row["partition"],
            row["ref_haplotype"],
            row["ref_contig"],
            int(row["ref_start"]),
            int(row["ref_end"]),
            row["strand"],
            float(row["alignment_score"]),
            float(row["aligned_score"]),
            int(row["aligned_bases"]),
            float(row["identity"]),
            row["source"],
            row["target_record"],
            tuple(_json_intervals(
                row["query_blocks"], f"{path}:{query_id}:query_blocks",
            )),
            tuple(_json_intervals(
                row["target_blocks"], f"{path}:{query_id}:target_blocks",
            )),
            tuple(_json_paired_blocks(
                row.get("paired_blocks", "[]"),
                f"{path}:{query_id}:paired_blocks",
            )),
            row.get("liftover_path") or ".",
        )
    return output


def read_small_novel_annotations(path: str) -> Dict[str, str]:
    """Read one annotation-only large/small novel-minset label per query."""
    output: Dict[str, str] = {}
    seen: Set[str] = set()
    for row in _read_dict_rows(path):
        query_id = row["query_id"]
        if query_id in seen:
            raise ValueError(f"{path}: duplicate query_id {query_id!r}")
        seen.add(query_id)
        label = row.get("label") or "."
        if row.get("status") == "mapped" and label != ".":
            output[query_id] = label
    return output


def read_alignment_summary_annotations(
    paths: Sequence[str],
) -> Dict[str, str]:
    """Load annotations from any number of local/global summary tables.

    A table may identify rows with ``query_id``, ``candidate_id`` or
    ``original_query_id`` and may publish its graph path as
    ``graphic_cigar``, ``liftover_path``, ``label`` or ``labels``. Multiple
    files and multiple rows per query are combined in command-line order.
    """
    output: Dict[str, str] = {}
    for path in paths:
        with open(path, "rt", newline="") as handle:
            first = handle.readline()
        if not first:
            raise ValueError(f"{path}: empty alignment summary")
        if first.startswith("#"):
            first = first[1:]
        fields = next(csv.reader([first], delimiter="\t"))
        if not (
            {"query_id", "candidate_id", "original_query_id"} & set(fields)
        ):
            raise ValueError(
                f"{path}: alignment summary lacks a query identifier"
            )
        for row in _read_dict_rows(path):
            query_id = (
                row.get("query_id") or row.get("candidate_id")
                or row.get("original_query_id")
            )
            if not query_id:
                continue
            label = next((
                row.get(field) for field in (
                    "graphic_cigar", "liftover_path", "label", "labels",
                )
                if row.get(field) not in {None, "", "."}
            ), None)
            if label:
                output[query_id] = merge_annotation_labels(
                    output.get(query_id), label,
                )
    return output


def merge_annotation_labels(*labels: Optional[str]) -> str:
    """Join distinct whitespace-free path annotations in stable order."""
    output: List[str] = []
    seen: Set[str] = set()
    for label in labels:
        if not label or label in {".", "alternative"}:
            continue
        if any(character.isspace() for character in label):
            raise ValueError(f"path annotation contains whitespace: {label!r}")
        if label not in seen:
            seen.add(label)
            output.append(label)
    return ";".join(output) if output else "."


def read_unique_rows(
    analysis_dir: str,
    partitions: Dict[str, PartitionFiles],
) -> List[UniqueRow]:
    path = os.path.join(analysis_dir, "all_unique_regions.bed")
    output: List[UniqueRow] = []
    for row in _read_dict_rows(path):
        partition = row["partition"]
        if partition not in partitions:
            raise ValueError(f"{path}: unknown partition {partition!r}")
        output.append(UniqueRow(
            partition,
            row["record_id"],
            row["haplotype"],
            row["source_contig"],
            int(row["source_start"]),
            int(row["source_end"]),
            int(row["local_start"]),
            int(row["local_end"]),
            row["region_name"],
            row["region_class"],
        ))
    return output


def read_status_rows(analysis_dir: str) -> List[StatusRow]:
    path = os.path.join(analysis_dir, "interval_status.tsv")
    return [
        StatusRow(
            row["haplotype"],
            row["source_contig"],
            int(row["source_start"]),
            int(row["source_end"]),
            row["interval_id"],
            row["status"],
            _parse_set(row["unlifted_partitions"]),
        )
        for row in _read_dict_rows(path)
    ]


def read_used_reference_rows(analysis_dir: str) -> List[UsedReferenceRow]:
    path = os.path.join(analysis_dir, "all_used_novel_loci.bed")
    return [
        UsedReferenceRow(
            row["partition"],
            row["haplotype"],
            row["source_contig"],
            int(row["source_start"]),
            int(row["source_end"]),
            row["encoded_name"],
        )
        for row in _read_dict_rows(path)
    ]


def parse_main_reference(path: str) -> Dict[str, ReferenceRecord]:
    records: Dict[str, ReferenceRecord] = {}
    with IndexedFasta(path) as indexed:
        lengths = {name: indexed.length(name) for name in indexed.names()}
    for line_number, description in fasta_header_descriptions(path):
        fields = description.split()
        if len(fields) < 2:
            raise ValueError(
                f"{path}:record {line_number}: main reference header lacks "
                "source coordinates"
            )
        record_id = fields[0]
        coordinate = None
        match = None
        for field in fields[1:]:
            candidate = field.removeprefix("source=")
            candidate_match = SOURCE_COORD_RE.fullmatch(candidate)
            if candidate_match is not None:
                coordinate = candidate
                match = candidate_match
                break
        if match is None:
            raise ValueError(
                f"{path}:record {line_number}: cannot parse cleaned-reference "
                "source coordinate from its FASTA header"
            )
        haplotype, contig, start_text, end_text, strand = match.groups()
        strand = strand or "+"
        start, end = int(start_text), int(end_text)
        if record_id not in lengths:
            raise ValueError(f"{path}: FASTA index lacks record {record_id!r}")
        if end - start != lengths[record_id]:
            raise ValueError(
                f"{path}:record {line_number}: source span {end - start} "
                f"does not equal sequence length {lengths[record_id]}"
            )
        if record_id in records:
            raise ValueError(f"{path}: duplicate record {record_id!r}")
        records[record_id] = ReferenceRecord(
            record_id, haplotype, contig, start, end, strand, lengths[record_id],
        )
    if not records:
        raise ValueError(f"{path}: no cleaned-reference coordinate records")
    return records


def global_queries_from_status(
    unique_rows: Sequence[UniqueRow],
    statuses: Sequence[StatusRow],
    partitions: Dict[str, PartitionFiles],
) -> List[NovelQuery]:
    """Intersect global status spans with unique records using sorted joins."""
    status_by_key: DefaultDict[Tuple[str, str], List[StatusRow]] = defaultdict(list)
    for row in statuses:
        if row.status == "global_novel":
            status_by_key[(row.haplotype, row.source_contig)].append(row)
    for rows in status_by_key.values():
        rows.sort(key=lambda row: (row.source_start, row.source_end))

    queries: List[NovelQuery] = []
    counter = 0
    seen: Set[Tuple[str, str, int, int, str]] = set()
    unique_by_key: DefaultDict[Tuple[str, str], List[UniqueRow]] = defaultdict(list)
    for row in unique_rows:
        unique_by_key[(row.haplotype, row.source_contig)].append(row)
    for key, members in unique_by_key.items():
        global_rows = status_by_key.get(key, ())
        if not global_rows:
            continue
        members.sort(key=lambda row: (
            row.source_start, row.source_end, row.partition, row.record_id,
        ))
        left = 0
        for unique in members:
            while (
                left < len(global_rows)
                and global_rows[left].source_end <= unique.source_start
            ):
                left += 1
            scan = left
            while (
                scan < len(global_rows)
                and global_rows[scan].source_start < unique.source_end
            ):
                status = global_rows[scan]
                if unique.partition in status.unlifted_partitions:
                    start = max(unique.source_start, status.source_start)
                    end = min(unique.source_end, status.source_end)
                    if end > start:
                        local_start = unique.local_start + start - unique.source_start
                        local_end = local_start + end - start
                        identity = (
                            unique.partition,
                            unique.record_id,
                            local_start,
                            local_end,
                            status.interval_id,
                        )
                        if identity in seen:
                            scan += 1
                            continue
                        seen.add(identity)
                        counter += 1
                        queries.append(NovelQuery(
                            f"global_{counter:012d}",
                            "global_novel",
                            unique.partition,
                            status.interval_id,
                            partitions[unique.partition].fasta,
                            unique.record_id,
                            unique.haplotype,
                            unique.source_contig,
                            start,
                            end,
                            local_start,
                            local_end,
                        ))
                scan += 1
    return queries


def global_queries_from_bed(
    analysis_dir: str,
    partitions: Dict[str, PartitionFiles],
) -> List[NovelQuery]:
    """Stream unique BED rows against the sorted global-novel BED index."""
    global_path = os.path.join(analysis_dir, "global_novel.bed")
    if not os.path.isfile(global_path):
        global_path = os.path.join(analysis_dir, "interval_status.tsv")
    unique_path = os.path.join(analysis_dir, "all_unique_regions.bed")
    status_by_key: DefaultDict[
        Tuple[str, str], List[StatusRow]
    ] = defaultdict(list)
    for row in _read_dict_rows(global_path):
        if row.get("status", "global_novel") != "global_novel":
            continue
        interval_id = row.get("interval_id") or row.get("global_novel_id")
        membership = (
            row.get("unlifted_partitions")
            if "unlifted_partitions" in row
            else row.get("partitions")
        )
        if not interval_id or membership is None:
            raise ValueError(
                f"{global_path}: expected interval/global-novel id and "
                "unlifted_partitions/partitions columns"
            )
        status_by_key[(row["haplotype"], row["source_contig"])].append(
            StatusRow(
                row["haplotype"],
                row["source_contig"],
                int(row["source_start"]),
                int(row["source_end"]),
                interval_id,
                "global_novel",
                _parse_set(membership),
            )
        )
    global_interval_count = sum(
        len(rows) for rows in status_by_key.values()
    )
    LOG.info(
        "Loaded %d global-novel intervals from %s",
        global_interval_count, global_path,
    )
    starts_by_key: Dict[Tuple[str, str], List[int]] = {}
    for key, rows in status_by_key.items():
        rows.sort(key=lambda item: (item.source_start, item.source_end))
        starts_by_key[key] = [item.source_start for item in rows]

    queries: List[NovelQuery] = []
    seen: Set[Tuple[str, str, int, int, str]] = set()
    counter = 0
    scanned = 0
    for row in _read_dict_rows(unique_path):
        scanned += 1
        if scanned % 1_000_000 == 0:
            LOG.info(
                "Global-query BED join: scanned %d unique intervals; "
                "selected %d query intervals",
                scanned, len(queries),
            )
        key = row["haplotype"], row["source_contig"]
        statuses = status_by_key.get(key)
        if not statuses:
            continue
        source_start = int(row["source_start"])
        source_end = int(row["source_end"])
        local_base = int(row["local_start"])
        starts = starts_by_key[key]
        scan = max(0, bisect.bisect_right(starts, source_start) - 1)
        while scan < len(statuses) and statuses[scan].source_end <= source_start:
            scan += 1
        while (
            scan < len(statuses)
            and statuses[scan].source_start < source_end
        ):
            status = statuses[scan]
            partition = row["partition"]
            if partition in status.unlifted_partitions:
                start = max(source_start, status.source_start)
                end = min(source_end, status.source_end)
                if end > start:
                    local_start = local_base + start - source_start
                    local_end = local_start + end - start
                    identity = (
                        partition,
                        row["record_id"],
                        local_start,
                        local_end,
                        status.interval_id,
                    )
                    if identity not in seen:
                        seen.add(identity)
                        counter += 1
                        queries.append(NovelQuery(
                            f"global_{counter:012d}",
                            "global_novel",
                            partition,
                            status.interval_id,
                            partitions[partition].fasta,
                            row["record_id"],
                            row["haplotype"],
                            row["source_contig"],
                            start,
                            end,
                            local_start,
                            local_end,
                        ))
            scan += 1
    LOG.info(
        "Global-query BED join complete: scanned %d unique intervals; "
        "selected %d query intervals",
        scanned, len(queries),
    )
    return queries


def _json_intervals(text: str, context: str) -> List[Interval]:
    try:
        raw = json.loads(text)
        intervals = [(int(start), int(end)) for start, end in raw]
    except (TypeError, ValueError, json.JSONDecodeError) as error:
        raise ValueError(f"{context}: invalid interval JSON") from error
    return merge_intervals(intervals)


def _json_paired_blocks(
    text: str,
    context: str,
) -> List[Tuple[int, int, int, int]]:
    try:
        raw = json.loads(text)
        blocks = [
            (int(q0), int(q1), int(r0), int(r1))
            for q0, q1, r0, r1 in raw
        ]
    except (TypeError, ValueError, json.JSONDecodeError) as error:
        raise ValueError(f"{context}: invalid paired-block JSON") from error
    for q0, q1, r0, r1 in blocks:
        if q0 < 0 or q1 <= q0 or r0 < 0 or r1 <= r0:
            raise ValueError(
                f"{context}: invalid paired block "
                f"{q0}-{q1}:{r0}-{r1}"
            )
        if q1 - q0 != r1 - r0:
            raise ValueError(
                f"{context}: unequal paired block lengths "
                f"{q0}-{q1}:{r0}-{r1}"
            )
    return sorted(blocks, key=lambda row: (row[0], row[1], row[2], row[3]))


def choose_local_placements(
    candidates: Sequence[MappingCandidate],
    placements: Sequence[MappedPlacement],
    reference_records: Optional[Dict[str, ReferenceRecord]] = None,
) -> Tuple[Dict[str, MappedPlacement], Dict[Tuple[str, str], List[Interval]]]:
    """Choose one local placement and return aligned query masks."""
    best: Dict[str, MappedPlacement] = {}
    rejected_candidates: Set[str] = set()
    reference_by_source = (
        reference_catalog_by_source(reference_records)
        if reference_records is not None else None
    )
    empty_mask_placements = 0
    for placement in placements:
        if not placement.query_blocks:
            raise ValueError(
                "local alignment results contain a blank query_blocks value; "
                "rerun "
                "the updated map_partition_local_novel.py"
            )
        if not _json_intervals(
            placement.query_blocks,
            f"{placement.candidate_id}:query_blocks",
        ):
            empty_mask_placements += 1
            continue
        if (
            reference_by_source is not None
            and not reference_span_is_covered(
                reference_by_source,
                placement.ref_haplotype,
                placement.ref_contig,
                placement.ref_start,
                placement.ref_end,
            )
        ):
            rejected_candidates.add(placement.candidate_id)
            continue
        prior = best.get(placement.candidate_id)
        if prior is None or placement_sort_key(placement) < placement_sort_key(prior):
            best[placement.candidate_id] = placement
    if empty_mask_placements:
        LOG.info(
            "Ignored %d local alignment rows that removed no qualifying "
            "query bases; their queries remain available for global mapping",
            empty_mask_placements,
        )
    if reference_by_source is not None:
        for candidate_id in sorted(rejected_candidates - set(best)):
            LOG.warning(
                "%s: all local placements are absent from cleaned "
                "main_chroms.fa; "
                "sending the complete query to global-reference alignment",
                candidate_id,
            )
    masks: DefaultDict[Tuple[str, str], List[Interval]] = defaultdict(list)
    candidate_by_id = {row.candidate_id: row for row in candidates}
    for candidate_id, placement in best.items():
        candidate = candidate_by_id.get(candidate_id)
        if candidate is None:
            raise ValueError(
                f"local placement references unknown candidate {candidate_id!r}"
            )
        for start, end in _json_intervals(
            placement.query_blocks, f"{candidate_id}:query_blocks",
        ):
            start = max(start, candidate.local_start)
            end = min(end, candidate.local_end)
            if end > start:
                masks[(candidate.novel_partition, candidate.record_id)].append(
                    (start, end)
                )
    return best, dict(masks)


def _annotate_local_partition(task):
    partition, members = task
    group_output: Dict[str, MappedPlacement] = {}
    readers: Dict[Tuple[str, str], IndexedFasta] = {}
    try:
        for candidate, placement, query_source, ref_source in sorted(
            members, key=lambda row: row[0].candidate_id,
        ):
            query_key = (
                query_source.fasta_path, query_source.fai_path,
            )
            ref_key = (
                ref_source.fasta_path, ref_source.fai_path,
            )
            query_fasta = readers.get(query_key)
            if query_fasta is None:
                query_fasta = IndexedFasta(*query_key)
                readers[query_key] = query_fasta
            reference = readers.get(ref_key)
            if reference is None:
                reference = IndexedFasta(*ref_key)
                readers[ref_key] = reference
            paired_blocks = _json_paired_blocks(
                placement.paired_blocks,
                f"{candidate.candidate_id}:paired_blocks",
            )
            if not paired_blocks:
                raise ValueError(
                    "local alignment results do not contain paired_blocks; "
                    "rerun the updated map_partition_local_novel.py"
                )
            sequence = query_fasta.fetch(
                candidate.source_contig,
                candidate.source_start,
                candidate.source_end,
            )
            path = build_graphic_path(
                query_id=candidate.candidate_id,
                query_start=candidate.local_start,
                query_end=candidate.local_end,
                query_sequence=sequence,
                ref_contig=placement.ref_contig,
                strand=placement.strand,
                paired_blocks=paired_blocks,
                reference=reference,
            )
            group_output[candidate.candidate_id] = dataclasses.replace(
                placement, liftover_path=path,
            )
    finally:
        for reader in readers.values():
            reader.close()
    return partition, group_output


def annotate_local_placements(
    candidates: Sequence[MappingCandidate],
    best: Dict[str, MappedPlacement],
    assembly_sources: Dict[str, object],
    jobs: int = 1,
) -> Dict[str, MappedPlacement]:
    """Add exact original-reference graphic paths to selected local hits."""
    candidate_by_id = {row.candidate_id: row for row in candidates}
    grouped: DefaultDict[str, List[Tuple[
        MappingCandidate,
        MappedPlacement,
        object,
        object,
    ]]] = defaultdict(list)
    for candidate_id, placement in best.items():
        candidate = candidate_by_id[candidate_id]
        if candidate.haplotype not in assembly_sources:
            raise ValueError(
                f"{candidate_id}: query_paths lacks query haplotype "
                f"{candidate.haplotype!r}"
            )
        if placement.ref_haplotype not in assembly_sources:
            raise ValueError(
                f"{candidate_id}: query_paths lacks reference haplotype "
                f"{placement.ref_haplotype!r}"
            )
        grouped[candidate.novel_partition].append(
            (
                candidate,
                placement,
                assembly_sources[candidate.haplotype],
                assembly_sources[placement.ref_haplotype],
            )
        )

    items = [
        (partition, members)
        for partition, members in sorted(grouped.items())
    ]
    output: Dict[str, MappedPlacement] = {}
    progress_step = max(1, len(items) // 20)
    for completed, (_partition, group_output) in enumerate(
        _bounded_process_map(
            _annotate_local_partition, items, max(1, jobs),
        ),
        1,
    ):
        output.update(group_output)
        if completed == len(items) or completed % progress_step == 0:
            LOG.info(
                "Local-placement annotation progress: %d/%d "
                "partitions with mapped local novels",
                completed, len(items),
            )
    return output


def _annotate_local_reporting_partition(task):
    partition, members = task
    output: Dict[str, str] = {}
    readers: Dict[Tuple[str, str], IndexedFasta] = {}
    try:
        for candidate, components, assembly_sources in members:
            query_source = assembly_sources.get(candidate.haplotype)
            if query_source is None:
                raise ValueError(
                    f"{candidate.candidate_id}: --query-paths lacks query "
                    f"haplotype {candidate.haplotype!r}"
                )
            query_key = (
                query_source.fasta_path, query_source.fai_path,
            )
            query_reader = readers.get(query_key)
            if query_reader is None:
                query_reader = IndexedFasta(*query_key)
                readers[query_key] = query_reader
            query_sequence = query_reader.fetch(
                candidate.source_contig,
                candidate.source_start,
                candidate.source_end,
            )
            reference_map = {}
            for component in components:
                source = assembly_sources.get(component.ref_haplotype)
                if source is None:
                    raise ValueError(
                        f"{candidate.candidate_id}: --query-paths lacks "
                        f"reference haplotype {component.ref_haplotype!r}"
                    )
                key = (source.fasta_path, source.fai_path)
                reader = readers.get(key)
                if reader is None:
                    reader = IndexedFasta(*key)
                    readers[key] = reader
                reference_map[(
                    component.ref_haplotype, component.ref_contig,
                )] = reader
            output[candidate.candidate_id] = build_multilocus_graphic_path(
                query_id=candidate.candidate_id,
                query_start=candidate.local_start,
                query_end=candidate.local_end,
                query_sequence=query_sequence,
                components=components,
                references=reference_map,
            )
    finally:
        for reader in readers.values():
            reader.close()
    return partition, output


def annotate_local_reporting_components(
    candidates: Sequence[MappingCandidate],
    components: Sequence[ResolvedPlacementComponent],
    assembly_sources: Dict[str, AssemblyFastaSource],
    jobs: int,
) -> Dict[str, str]:
    candidate_by_id = {row.candidate_id: row for row in candidates}
    by_candidate: DefaultDict[str, List[ResolvedPlacementComponent]] = (
        defaultdict(list)
    )
    for component in components:
        by_candidate[component.query_id].append(component)
    grouped: DefaultDict[str, List[Tuple]] = defaultdict(list)
    for candidate_id, members in by_candidate.items():
        candidate = candidate_by_id.get(candidate_id)
        if candidate is None:
            raise ValueError(
                "local reporting component references unknown candidate "
                f"{candidate_id!r}"
            )
        members.sort(key=lambda row: (
            row.query_start, row.query_end, row.component_index,
        ))
        grouped[candidate.novel_partition].append((
            candidate, tuple(members), assembly_sources,
        ))
    tasks = sorted(grouped.items())
    output: Dict[str, str] = {}
    for _partition, labels in _bounded_process_map(
        _annotate_local_reporting_partition, tasks, max(1, jobs),
    ):
        output.update(labels)
    return output


def reporting_liftover_labels(
    components: Sequence[ResolvedPlacementComponent],
) -> Dict[str, str]:
    output: Dict[str, str] = {}
    for component in components:
        path = component.liftover_path
        if not path or path == ".":
            continue
        prior = output.get(component.query_id)
        if prior is not None and prior != path:
            raise ValueError(
                f"{component.query_id}: inconsistent multi-locus paths in "
                "reporting table"
            )
        output[component.query_id] = path
    return output


def promoted_local_queries(
    candidates: Sequence[MappingCandidate],
    best: Dict[str, MappedPlacement],
) -> List[NovelQuery]:
    """Return local-query sequence not removed by its single best placement."""
    output: List[NovelQuery] = []
    counter = 0
    for candidate in candidates:
        placement = best.get(candidate.candidate_id)
        aligned = (
            _json_intervals(
                placement.query_blocks,
                f"{candidate.candidate_id}:query_blocks",
            )
            if placement is not None else []
        )
        residual = subtract_intervals(
            [(candidate.local_start, candidate.local_end)], aligned,
        )
        for local_start, local_end in residual:
            counter += 1
            source_start = (
                candidate.source_start + local_start - candidate.local_start
            )
            output.append(NovelQuery(
                f"promoted_local_{counter:012d}",
                "local_residual" if placement is not None else "unresolved_local",
                candidate.novel_partition,
                candidate.interval_id,
                candidate.source_fasta,
                candidate.record_id,
                candidate.haplotype,
                candidate.source_contig,
                source_start,
                source_start + local_end - local_start,
                local_start,
                local_end,
            ))
    return output


def global_mapping_candidates(
    queries: Sequence[NovelQuery],
) -> List[MappingCandidate]:
    """Adapt global queries to the shared staged mapping engine."""
    return [
        MappingCandidate(
            query.query_id,
            query.interval_id,
            query.source_contig,
            query.source_start,
            query.source_end,
            query.haplotype,
            query.partition,
            frozenset({"global_reference"}),
            query.record_id,
            query.local_start,
            query.local_end,
            query.source_fasta,
            query.origin,
        )
        for query in queries
    ]


def filter_queries_by_score(
    queries: Sequence[NovelQuery],
    min_score: float,
    cached_queries: Optional[
        Mapping[str, CachedCandidateSequence]
    ] = None,
) -> List[NovelQuery]:
    if cached_queries is not None:
        retained = []
        for query in queries:
            cache_id = query.query_id.split(".residual.", 1)[0]
            cached = cached_queries[cache_id]
            sequence = cached.sequence[
                query.local_start - cached.local_start:
                query.local_end - cached.local_start
            ]
            if novelty_score(sequence, MASKED_BASE_WEIGHT) >= min_score:
                retained.append(query)
        return retained

    by_fasta: DefaultDict[str, List[NovelQuery]] = defaultdict(list)
    for query in queries:
        by_fasta[query.source_fasta].append(query)
    retained: List[NovelQuery] = []
    for fasta, members in sorted(by_fasta.items()):
        with IndexedFasta(fasta) as source:
            for query in members:
                sequence = source.fetch(
                    query.record_id, query.local_start, query.local_end,
                )
                if novelty_score(sequence, MASKED_BASE_WEIGHT) >= min_score:
                    retained.append(query)
    return retained


def filter_queries_for_global_lift(
    queries: Sequence[NovelQuery],
    min_unmasked: int,
    cached_queries: Mapping[str, CachedCandidateSequence],
) -> List[NovelQuery]:
    """Keep only queries that strictly exceed the unmasked lift threshold."""
    retained: List[NovelQuery] = []
    for query in queries:
        cache_id = query.query_id.split(".residual.", 1)[0]
        cached = cached_queries[cache_id]
        sequence = cached.sequence[
            query.local_start - cached.local_start:
            query.local_end - cached.local_start
        ]
        if count_unmasked(sequence) > min_unmasked:
            retained.append(query)
    return retained


def global_residual_queries(
    queries: Sequence[NovelQuery],
    placements: Dict[str, GlobalPlacement],
    min_score: float,
    cached_queries: Mapping[str, CachedCandidateSequence],
) -> List[NovelQuery]:
    """Return score-passing query pieces not covered by global placements."""
    residuals: List[NovelQuery] = []
    for query in queries:
        placement = placements.get(query.query_id)
        uncovered = subtract_intervals(
            [(query.local_start, query.local_end)],
            placement.query_blocks if placement is not None else (),
        )
        for local_start, local_end in uncovered:
            source_start = (
                query.source_start + local_start - query.local_start
            )
            residuals.append(dataclasses.replace(
                query,
                query_id=(
                    f"{query.query_id}.residual."
                    f"{local_start}-{local_end}"
                ),
                source_start=source_start,
                source_end=source_start + local_end - local_start,
                local_start=local_start,
                local_end=local_end,
            ))
    return filter_queries_by_score(
        residuals, min_score, cached_queries,
    )


def _reference_source_blocks(
    record: ReferenceRecord,
    target_blocks: Sequence[Interval],
) -> List[Interval]:
    output: List[Interval] = []
    for start, end in target_blocks:
        if record.strand == "-":
            output.append((
                record.source_end - end,
                record.source_end - start,
            ))
        else:
            output.append((
                record.source_start + start,
                record.source_start + end,
            ))
    return merge_intervals(output)


def select_global_placements(
    hits: Sequence[AlignmentHit],
    query_meta: Dict[str, QuerySequence],
    reference_records: Dict[str, ReferenceRecord],
    min_score: float,
) -> Dict[str, GlobalPlacement]:
    """Select exactly one highest-scoring reference placement per query."""
    choices: DefaultDict[str, List[GlobalPlacement]] = defaultdict(list)
    for hit in hits:
        meta = query_meta.get(hit.query_id)
        target = reference_records.get(hit.target_id)
        if meta is None or target is None:
            continue
        query_blocks, target_blocks = clip_hit_to_query_core(
            hit, meta.core_start, meta.core_end,
        )
        if not query_blocks or not target_blocks:
            continue
        aligned_score = sum(
            novelty_score(
                meta.sequence[start:end], MASKED_BASE_WEIGHT,
            )
            for start, end in query_blocks
        )
        aligned_bases = sum(end - start for start, end in query_blocks)
        if aligned_score < min_score:
            continue
        local_query_blocks = tuple(merge_intervals([
            (
                meta.query.local_start + start - meta.core_start,
                meta.query.local_start + end - meta.core_start,
            )
            for start, end in query_blocks
        ]))
        source_target_blocks = tuple(
            _reference_source_blocks(target, target_blocks)
        )
        if not source_target_blocks:
            continue
        choices[hit.query_id].append(GlobalPlacement(
            hit.query_id,
            meta.query.partition,
            target.haplotype,
            target.source_contig,
            min(start for start, _end in source_target_blocks),
            max(end for _start, end in source_target_blocks),
            hit.strand,
            float(hit.alignment_score),
            aligned_score,
            aligned_bases,
            float(hit.identity),
            hit.source,
            hit.target_id,
            local_query_blocks,
            source_target_blocks,
        ))
    selected: Dict[str, GlobalPlacement] = {}
    for query_id, rows in choices.items():
        selected[query_id] = min(rows, key=lambda row: (
            -row.alignment_score,
            -row.aligned_score,
            -row.aligned_bases,
            -row.identity,
            row.ref_haplotype,
            row.ref_contig,
            row.ref_start,
            row.ref_end,
            row.target_record,
            row.source,
        ))
    return selected


def align_global_queries(
    queries: Sequence[NovelQuery],
    candidates: Sequence[MappingCandidate],
    cached_queries: Mapping[str, CachedCandidateSequence],
    main_reference: str,
    reference_records: Dict[str, ReferenceRecord],
    workdir: str,
    args: argparse.Namespace,
    mapping_signature: str,
) -> Dict[str, GlobalPlacement]:
    """Run one cohort-wide staged residual alignment against cleaned Ref."""
    if not queries:
        return {}
    reuse_alignments_only = bool(getattr(
        args, "reuse_alignments_only", False,
    ))
    if not reuse_alignments_only:
        for executable in ("minimap2", "winnowmap"):
            ensure_executable(executable)
    if not args.skip_blastn and not reuse_alignments_only:
        for executable in ("makeblastdb", "blastn"):
            ensure_executable(executable)
    LOG.info(
        "Global alignment: %d queries in one combined FASTA against one "
        "cleaned-reference target using %d cores",
        len(queries), args.cores,
    )

    targets = tuple(
        BlockTarget(
            record.record_id,
            record.record_id,
            0,
            record.length,
            record.haplotype,
            record.source_contig,
            record.source_start,
            record.source_end,
            record.strand,
        )
        for record in sorted(
            reference_records.values(),
            key=lambda row: row.record_id,
        )
    )
    reference_signature = hashlib.sha256(json.dumps(
        {
            "main_reference": _fingerprint(main_reference),
            "targets": [dataclasses.asdict(target) for target in targets],
        },
        sort_keys=True,
    ).encode()).hexdigest()
    Path(workdir).mkdir(parents=True, exist_ok=True)
    reference_dir = os.path.join(workdir, "reference")
    target_fasta = os.path.join(reference_dir, "cleaned_reference.fa")
    reference_marker = os.path.join(reference_dir, "_SUCCESS.json")
    reusable_reference = False
    if args.resume and os.path.isfile(reference_marker):
        try:
            with open(reference_marker, "rt") as handle:
                reusable_reference = (
                    json.load(handle).get("signature") == reference_signature
                    and os.path.islink(target_fasta)
                    and os.path.realpath(target_fasta)
                    == os.path.realpath(main_reference)
                )
        except (OSError, ValueError, AttributeError):
            reusable_reference = False
    if not reusable_reference and reuse_alignments_only:
        raise RuntimeError(
            f"report-only replay cannot reuse fixed reference {target_fasta}"
        )
    if not reusable_reference:
        if os.path.isdir(reference_dir):
            shutil.rmtree(reference_dir)
        Path(reference_dir).mkdir(parents=True, exist_ok=True)
        os.symlink(main_reference, target_fasta)
        with open(reference_marker, "wt") as out:
            json.dump({"signature": reference_signature}, out, sort_keys=True)
            out.write("\n")
    reference = PreparedBlockReference(
        "global_reference",
        target_fasta,
        targets,
        reference_signature,
    )
    options = {
        "resume": args.resume,
        "mapping_signature": mapping_signature,
        "cores": args.cores,
        "input_anchor": args.alignment_anchor,
        "min_identity": args.min_identity,
        "min_score": args.min_score,
        "minimap_query_batch": getattr(
            args, "minimap_query_batch", None,
        ),
        "minimap_retry_sigkill": True,
        "minimap_retry_threads": min(args.cores, 32),
        "minimap_retry_query_batch": "100M",
        "blast_word_size": args.blast_word_size,
        "blast_evalue": args.blast_evalue,
        "blast_max_target_seqs": args.blast_max_target_seqs,
        "skip_blastn": args.skip_blastn,
        "converge_aligners": True,
        "reuse_alignments_only": reuse_alignments_only,
    }
    result_path = map_partition_candidates(
        candidates,
        cached_queries,
        reference,
        workdir,
        options,
    )
    first_cycle_queries = os.path.join(
        workdir, "partitions", "global_reference",
        "minimap2", "cycle000", "queries.fa",
    )
    combined_query_fasta = os.path.join(workdir, "global_queries.fa")
    if os.path.lexists(combined_query_fasta):
        os.remove(combined_query_fasta)
    os.symlink(first_cycle_queries, combined_query_fasta)
    query_by_id = {query.query_id: query for query in queries}
    selected: Dict[str, GlobalPlacement] = {}
    for placement in read_placements((result_path,)):
        query = query_by_id[placement.candidate_id]
        query_blocks = tuple(_json_intervals(
            placement.query_blocks,
            f"{placement.candidate_id}:global_query_blocks",
        ))
        target_blocks = tuple(_json_intervals(
            placement.target_blocks,
            f"{placement.candidate_id}:global_target_blocks",
        ))
        paired_blocks = tuple(_json_paired_blocks(
            placement.paired_blocks,
            f"{placement.candidate_id}:global_paired_blocks",
        ))
        cached = cached_queries[placement.candidate_id]
        aligned_score = sum(
            novelty_score(
                cached.sequence[
                    start - cached.local_start:end - cached.local_start
                ],
                MASKED_BASE_WEIGHT,
            )
            for start, end in query_blocks
        )
        selected[placement.candidate_id] = GlobalPlacement(
            placement.candidate_id,
            query.partition,
            placement.ref_haplotype,
            placement.ref_contig,
            placement.ref_start,
            placement.ref_end,
            placement.strand,
            placement.alignment_score,
            aligned_score,
            placement.aligned_bases,
            placement.identity,
            placement.source,
            placement.target_record,
            query_blocks,
            target_blocks,
            paired_blocks,
            ".",
        )
    return selected


def reference_catalog_by_source(
    records: Dict[str, ReferenceRecord],
) -> Dict[Tuple[str, str], List[ReferenceRecord]]:
    output: DefaultDict[Tuple[str, str], List[ReferenceRecord]] = defaultdict(list)
    for record in records.values():
        output[(record.haplotype, record.source_contig)].append(record)
    for members in output.values():
        members.sort(key=lambda row: (
            row.source_start, row.source_end, row.record_id,
        ))
    return dict(output)


def reference_span_is_covered(
    catalog: Dict[Tuple[str, str], List[ReferenceRecord]],
    haplotype: str,
    contig: str,
    start: int,
    end: int,
) -> bool:
    covered = merge_intervals([
        (max(start, row.source_start), min(end, row.source_end))
        for row in catalog.get((haplotype, contig), ())
        if row.source_start < end and start < row.source_end
    ])
    return sum(right - left for left, right in covered) == end - start


def reference_segments(
    catalog: Dict[Tuple[str, str], List[ReferenceRecord]],
    fasta: str,
    haplotype: str,
    contig: str,
    start: int,
    end: int,
) -> List[Dict[str, object]]:
    segments: List[Dict[str, object]] = []
    for record in catalog.get((haplotype, contig), ()):
        if record.source_end <= start:
            continue
        if record.source_start >= end:
            break
        left = max(start, record.source_start)
        right = min(end, record.source_end)
        if right <= left:
            continue
        if record.strand == "-":
            local_start = record.source_end - right
            local_end = record.source_end - left
        else:
            local_start = left - record.source_start
            local_end = right - record.source_start
        segments.append({
            "fasta": os.path.abspath(fasta),
            "record": record.record_id,
            "start": local_start,
            "end": local_end,
            "strand": record.strand,
            "source_start": left,
            "source_end": right,
        })
    return segments


@lru_cache(maxsize=None)
def _fai_lengths(path: str) -> Dict[str, int]:
    with open(path, "rt") as handle:
        return {
            fields[0]: int(fields[1])
            for raw in handle
            if (fields := raw.rstrip("\n").split("\t"))
            and len(fields) >= 2
        }


def complete_reference_segments(
    catalog: Dict[Tuple[str, str], List[ReferenceRecord]],
    main_reference: str,
    assembly_sources: Mapping[str, AssemblyFastaSource],
    haplotype: str,
    contig: str,
    start: int,
    end: int,
) -> List[Dict[str, object]]:
    """Prefer main_chroms records and fill uncovered coordinates from Ref."""
    segments = reference_segments(
        catalog, main_reference, haplotype, contig, start, end,
    )
    covered = merge_intervals([
        (int(segment["source_start"]), int(segment["source_end"]))
        for segment in segments
    ])
    missing = subtract_intervals([(start, end)], covered)
    if missing:
        source = assembly_sources.get(haplotype)
        if source is None:
            raise ValueError(
                f"query_paths lacks reference haplotype {haplotype!r}, "
                f"needed to materialize {contig}:{start}-{end}"
            )
        contig_length = _fai_lengths(source.fai_path).get(contig)
        if contig_length is None:
            raise ValueError(
                f"{source.fai_path}: reference contig {contig!r} is absent"
            )
        if end > contig_length:
            raise ValueError(
                f"{haplotype}:{contig}:{start}-{end} exceeds original "
                f"reference length {contig_length}"
            )
        for left, right in missing:
            segments.append({
                "fasta": os.path.abspath(source.fasta_path),
                "fai": os.path.abspath(source.fai_path),
                "record": contig,
                "start": left,
                "end": right,
                "strand": "+",
                "source_start": left,
                "source_end": right,
            })
    segments.sort(key=lambda row: (
        int(row["source_start"]), int(row["source_end"]),
        str(row["fasta"]), str(row["record"]),
    ))
    return segments


def merge_reference_spans(
    spans: Sequence[ReferenceSpan],
    merge_gap: int,
) -> List[ReferenceSpan]:
    ordered = sorted(spans, key=lambda row: (
        row.haplotype, row.contig, row.start, row.end,
    ))
    output: List[ReferenceSpan] = []
    for row in ordered:
        if (
            output
            and output[-1].haplotype == row.haplotype
            and output[-1].contig == row.contig
            and row.start - output[-1].end <= merge_gap
        ):
            prior = output[-1]
            prior.end = max(prior.end, row.end)
            prior.score = max(prior.score, row.score)
            prior.origins.update(row.origins)
        else:
            output.append(ReferenceSpan(
                row.haplotype, row.contig, row.start, row.end,
                row.score, set(row.origins),
            ))
    return output


def partition_cleaned_reference_chains(
    rows: Sequence[UniqueRow], span: ReferenceSpan,
) -> List[Tuple[int, int, List[Dict[str, object]], Set[str]]]:
    """Build gap-free reference paths from this partition's unique rows.

    Source-coordinate gaps start separate paths. They are never filled from
    the cohort-level cleaned reference or from the original assembly.
    """
    clipped: List[Tuple[int, int, int, int, UniqueRow]] = []
    for row in rows:
        left = max(span.start, row.source_start)
        right = min(span.end, row.source_end)
        if right <= left:
            continue
        local_start = row.local_start + left - row.source_start
        clipped.append((
            left, right, local_start, local_start + right - left, row,
        ))
    clipped.sort(key=lambda item: (
        item[0], item[1], item[4].record_id, item[2], item[3],
    ))

    output: List[Tuple[int, int, List[Dict[str, object]], Set[str]]] = []
    chain_start: Optional[int] = None
    chain_end: Optional[int] = None
    segments: List[Dict[str, object]] = []
    origins: Set[str] = set(span.origins)

    def flush() -> None:
        nonlocal chain_start, chain_end, segments, origins
        if chain_start is not None and chain_end is not None and segments:
            output.append((chain_start, chain_end, segments, origins))
        chain_start = chain_end = None
        segments = []
        origins = set(span.origins)

    for left, right, local_start, local_end, row in clipped:
        if chain_end is not None and left > chain_end:
            flush()
        if chain_end is not None and left < chain_end:
            trim = chain_end - left
            left += trim
            local_start += trim
        if right <= left:
            continue
        if chain_start is None:
            chain_start = left
            chain_end = left
        segments.append({
            "record": row.record_id,
            "start": local_start,
            "end": local_end,
            "strand": "+",
            "source_start": left,
            "source_end": right,
        })
        chain_end = right
        origins.add(f"reference_unique:{row.region_name}")
    flush()
    return output


def _overlaps_original_reference(
    span: ReferenceSpan,
    original: Sequence[UsedReferenceRow],
) -> bool:
    return any(
        row.haplotype == span.haplotype
        and row.source_contig == span.contig
        and row.source_start < span.end
        and span.start < row.source_end
        for row in original
    )


def partition_header_segments(
    files: PartitionFiles,
    headers: Sequence[HeaderRecord],
    row: UsedReferenceRow,
) -> List[Dict[str, object]]:
    sequence_fasta = partition_sample_fasta(files)
    segments: List[Dict[str, object]] = []
    for header in headers:
        if (
            header.haplotype != row.haplotype
            or header.source_contig != row.source_contig
            or header.source_end <= row.source_start
            or header.source_start >= row.source_end
        ):
            continue
        left = max(header.source_start, row.source_start)
        right = min(header.source_end, row.source_end)
        segments.append({
            "fasta": os.path.abspath(sequence_fasta),
            "record": header.record_id,
            "start": left - header.source_start,
            "end": right - header.source_start,
            "strand": "+",
            "source_start": left,
            "source_end": right,
        })
    return segments


def partition_sample_fasta(files: PartitionFiles) -> str:
    """Return the FASTA that owns records described by PARTITION.header.

    Current local graphs keep the block template in ``PARTITION.fasta`` and
    all assembly paths in ``PARTITION_samples.fasta``.  The authoritative
    ``PARTITION.header`` and ``PARTITION.unique.bed`` describe the latter.
    Retain the old single-FASTA fallback for legacy/test partitions.
    """
    candidate = expected_partition_sample_fasta(files)
    return candidate if os.path.isfile(candidate) else files.fasta


def expected_partition_sample_fasta(files: PartitionFiles) -> str:
    """Return the canonical sample path even when it is intentionally absent."""
    stem, extension = os.path.splitext(files.fasta)
    return stem + "_samples" + extension


def _reference_length_for_placement(
    haplotype: str,
    contig: str,
    placement_end: int,
    assembly_sources: Optional[Mapping[str, AssemblyFastaSource]],
) -> int:
    """Return the original target length needed for a semiglobal path."""
    if assembly_sources is not None and haplotype in assembly_sources:
        source = assembly_sources[haplotype]
        length = _fai_lengths(source.fai_path).get(contig)
        if length is not None:
            return length
    return placement_end


def _placement_paired_blocks(placement) -> List[Tuple[int, int, int, int]]:
    value = placement.paired_blocks
    if isinstance(value, str):
        return _json_paired_blocks(
            value, f"{getattr(placement, 'candidate_id', '?')}:paired_blocks",
        )
    return sorted(tuple(value), key=lambda row: (row[0], row[1], row[2], row[3]))


def insertion_graphic_path(
    placement,
    query_start: int,
    query_end: int,
    assembly_sources: Optional[Mapping[str, AssemblyFastaSource]],
) -> str:
    """Place an uncovered query piece as ``I`` beside its accepted target.

    The aligned part of a connected mapped-novel query is removed from the
    alternative FASTA.  A retained flank is therefore an insertion relative
    to the refined reference, not an unmapped record with a ``.`` label.
    """
    length = query_end - query_start
    if length <= 0:
        return "."
    ref_length = _reference_length_for_placement(
        placement.ref_haplotype, placement.ref_contig,
        placement.ref_end, assembly_sources,
    )
    oriented: List[Tuple[int, int, int, int]] = []
    for q0, q1, r0, r1 in _placement_paired_blocks(placement):
        if placement.strand == "+":
            o0, o1 = r0, r1
        else:
            o0, o1 = ref_length - r1, ref_length - r0
        oriented.append((q0, q1, o0, o1))
    oriented.sort(key=lambda row: (row[0], row[1]))
    previous = [row for row in oriented if row[1] <= query_start]
    following = [row for row in oriented if row[0] >= query_end]
    if previous:
        anchor = previous[-1][3]
    elif following:
        anchor = following[0][2]
    elif placement.strand == "+":
        anchor = placement.ref_start
    else:
        anchor = ref_length - placement.ref_end
    anchor = max(0, min(ref_length, anchor))
    marker = ">" if placement.strand == "+" else "<"
    operations = []
    if anchor:
        operations.append(f"{anchor}H")
    operations.append(f"{length}I")
    if ref_length - anchor:
        operations.append(f"{ref_length - anchor}H")
    return f"{marker}{placement.ref_contig}:{''.join(operations)}"


def complete_graphic_query_path(path: str, query_length: int) -> str:
    """Pad a legacy aligned-only graph path with query insertions."""
    if query_length <= 0:
        return path
    if not path or path == ".":
        return f"{query_length}I"
    consumed = sum(
        int(length)
        for length, operation in re.findall(r"(\d+)([MIDNSHP=X])", path)
        if operation in {"M", "I", "=", "X"}
    )
    missing = query_length - consumed
    if missing <= 0:
        return path
    match = re.search(r"(\d+H)$", path)
    if match:
        return path[:match.start()] + f"{missing}I" + match.group(1)
    return path + f"{missing}I"


def build_adjusted_paths(
    partitions: Dict[str, PartitionFiles],
    unique_rows: Sequence[UniqueRow],
    used_rows: Sequence[UsedReferenceRow],
    local_candidates: Sequence[MappingCandidate],
    local_best: Dict[str, MappedPlacement],
    local_masks: Dict[Tuple[str, str], List[Interval]],
    global_queries: Sequence[NovelQuery],
    global_best: Dict[str, GlobalPlacement],
    reference_records: Dict[str, ReferenceRecord],
    main_reference: str,
    merge_gap: int,
    min_score: float,
    assembly_sources: Optional[Mapping[str, AssemblyFastaSource]] = None,
    jobs: int = 1,
    small_novel_annotations: Optional[Mapping[str, str]] = None,
    additional_annotations: Optional[Mapping[str, str]] = None,
    covered_records: Optional[List[Dict[str, object]]] = None,
    precovered_masks: Optional[
        Mapping[Tuple[str, str], Sequence[Interval]]
    ] = None,
    force_reference_rows: bool = False,
    include_reference_unique_spans: bool = True,
    suppress_unique_haplotypes: Optional[Iterable[str]] = None,
    prevalidated_query_ids: Optional[Iterable[str]] = None,
) -> List[PathRow]:
    """Build refined paths and remove query bases promoted to reference.

    Reference connectivity is transitive on one reference contig.  Original
    block spans seed a component; selected local/global placements join it
    when separated by at most ``merge_gap``.  Only placements in a seeded
    component refine the reference and mask their aligned query blocks.
    """
    # Reference paths are now sourced exclusively from partition-local unique
    # rows; the cohort reference inputs remain API metadata for callers.
    del local_masks, reference_records, main_reference, force_reference_rows
    novel_labels = dict(small_novel_annotations or {})
    for query_id, label in (additional_annotations or {}).items():
        novel_labels[query_id] = merge_annotation_labels(
            novel_labels.get(query_id), label,
        )
    audit = covered_records if covered_records is not None else []
    used_by_partition: DefaultDict[str, List[UsedReferenceRow]] = defaultdict(list)
    unique_by_partition: DefaultDict[str, List[UniqueRow]] = defaultdict(list)
    for row in used_rows:
        used_by_partition[row.partition].append(row)
    for row in unique_rows:
        unique_by_partition[row.partition].append(row)
    candidate_by_id = {row.candidate_id: row for row in local_candidates}
    query_by_id = {row.query_id: row for row in global_queries}
    suppressed_unique_haplotypes = set(suppress_unique_haplotypes or ())
    score_validated_queries = set(prevalidated_query_ids or ())
    from fixed_alternatives import FIXED_ROLE, fixed_partitions
    fixed_by_partition = fixed_partitions(partitions)

    def build_partition(partition: str):
        original = used_by_partition.get(partition, ())
        if not original:
            raise ValueError(
                f"partition {partition}: expected at least one block-reference interval"
            )
        fixed = fixed_by_partition.get(partition, {})
        fixed_only = bool(fixed) and all(row.encoded_name in fixed for row in original)
        fixed_sources = {(item['source_haplotype'], item['source_contig'])
                         for item in fixed.values()}
        def may_refine(unique):
            if fixed_only:
                return False
            if (unique.haplotype, unique.source_contig) not in fixed_sources:
                return True
            # A separate novel template may share the imported record's
            # provenance sample/contig. Its policy is still independent.
            return any(row.encoded_name not in fixed and
                       (row.haplotype, row.source_contig) == (unique.haplotype, unique.source_contig) and
                       row.source_start < unique.source_end and unique.source_start < row.source_end
                       for row in original)
        spans: List[ReferenceSpan] = [
            ReferenceSpan(
                row.haplotype, row.source_contig, row.source_start,
                row.source_end, 0.0, {f"original:{row.encoded_name}"},
            )
            for row in original
        ]
        if include_reference_unique_spans and not fixed_only:
            for unique in unique_by_partition.get(partition, ()):
                # Every reference-class row in the partition-cleaned union
                # participates in coordinate connectivity. Other unique rows
                # remain query alternatives.
                if unique.region_class != "reference":
                    continue
                if not may_refine(unique):
                    continue
                spans.append(ReferenceSpan(
                    unique.haplotype, unique.source_contig,
                    unique.source_start, unique.source_end, 0.0,
                    {f"reference_unique:{unique.region_name}"},
                ))
        for candidate_id, placement in local_best.items():
            candidate = candidate_by_id[candidate_id]
            if candidate.novel_partition == partition:
                spans.append(ReferenceSpan(
                    placement.ref_haplotype, placement.ref_contig,
                    placement.ref_start, placement.ref_end,
                    placement.alignment_score, {f"local:{candidate_id}"},
                ))
        for query_id, placement in global_best.items():
            query = query_by_id[query_id]
            if query.partition == partition:
                spans.append(ReferenceSpan(
                    placement.ref_haplotype, placement.ref_contig,
                    placement.ref_start, placement.ref_end,
                    placement.alignment_score, {f"global:{query_id}"},
                ))
        components = merge_reference_spans(spans, merge_gap)
        connected = [
            span for span in components
            if any(origin.startswith("original:") for origin in span.origins)
        ]
        connected_local = {
            origin.split(":", 1)[1]
            for span in connected for origin in span.origins
            if origin.startswith("local:")
        }
        connected_global = {
            origin.split(":", 1)[1]
            for span in connected for origin in span.origins
            if origin.startswith("global:")
        }

        masks: DefaultDict[Tuple[str, str], List[Interval]] = defaultdict(list)
        for key, intervals in (precovered_masks or {}).items():
            if key[0] == partition:
                masks[key].extend(intervals)
        # A selected template is partition-local and may come from any
        # haplotype, including a large novel locus. Its source bases must not
        # be emitted again through an explicit local/global novel group.
        for unique in unique_by_partition.get(partition, ()):
            for span in connected:
                if (
                    unique.haplotype != span.haplotype
                    or unique.source_contig != span.contig
                ):
                    continue
                source_start = max(unique.source_start, span.start)
                source_end = min(unique.source_end, span.end)
                if source_end <= source_start:
                    continue
                masks[(partition, unique.record_id)].append((
                    unique.local_start + source_start - unique.source_start,
                    unique.local_start + source_end - unique.source_start,
                ))
        local_global_labels: DefaultDict[str, List[str]] = defaultdict(list)
        handled_globals: Set[str] = set()

        def record_covered(
            kind: str, query_id: str, record_id: str,
            haplotype: str, source_contig: str, source_base: int,
            local_base: int, blocks: Sequence[Interval], placement,
        ) -> None:
            for start, end in merge_intervals(blocks):
                masks[(partition, record_id)].append((start, end))
                audit.append({
                    "partition": partition,
                    "query_kind": kind,
                    "query_id": query_id,
                    "record_id": record_id,
                    "haplotype": haplotype,
                    "source_contig": source_contig,
                    "source_start": source_base + start - local_base,
                    "source_end": source_base + end - local_base,
                    "local_start": start,
                    "local_end": end,
                    "ref_haplotype": placement.ref_haplotype,
                    "ref_contig": placement.ref_contig,
                    "ref_start": placement.ref_start,
                    "ref_end": placement.ref_end,
                    "strand": placement.strand,
                    "reason": "connected_reference_refinement",
                })

        for candidate_id in sorted(connected_local):
            candidate = candidate_by_id[candidate_id]
            placement = local_best[candidate_id]
            record_covered(
                "local_novel", candidate_id, candidate.record_id,
                candidate.haplotype, candidate.source_contig,
                candidate.source_start, candidate.local_start,
                _json_intervals(
                    placement.query_blocks, f"{candidate_id}:query_blocks",
                ), placement,
            )
        for query_id in sorted(connected_global):
            query = query_by_id[query_id]
            placement = global_best[query_id]
            record_covered(
                "global_novel", query_id, query.record_id,
                query.haplotype, query.source_contig,
                query.source_start, query.local_start,
                placement.query_blocks, placement,
            )

        # A mapped local candidate remains the owning alternative.  Any
        # promoted global pieces inside it contribute masks/labels but are not
        # emitted a second time.
        groups: DefaultDict[Tuple[str, str], List[Dict[str, object]]] = defaultdict(list)
        globals_by_record: DefaultDict[str, List[NovelQuery]] = defaultdict(list)
        for query in global_queries:
            if query.partition == partition:
                globals_by_record[query.record_id].append(query)
        for candidate_id, placement in sorted(local_best.items()):
            candidate = candidate_by_id[candidate_id]
            if candidate.novel_partition != partition:
                continue
            labels = [
                complete_graphic_query_path(
                    placement.liftover_path,
                    candidate.local_end - candidate.local_start,
                ),
                novel_labels.get(candidate_id),
            ]
            for query in globals_by_record.get(candidate.record_id, ()):
                if not (
                    candidate.local_start <= query.local_start
                    and query.local_end <= candidate.local_end
                ):
                    continue
                handled_globals.add(query.query_id)
                global_placement = global_best.get(query.query_id)
                if query.query_id in connected_global and global_placement is not None:
                    for start, end in subtract_intervals(
                        [(query.local_start, query.local_end)],
                        global_placement.query_blocks,
                    ):
                        local_global_labels[candidate_id].append(
                            insertion_graphic_path(
                                global_placement, start, end, assembly_sources,
                            )
                        )
                elif global_placement is not None:
                    local_global_labels[candidate_id].append(
                        global_placement.liftover_path
                    )
                if query.query_id in novel_labels:
                    local_global_labels[candidate_id].append(
                        novel_labels[query.query_id]
                    )
            labels.extend(local_global_labels[candidate_id])
            excluded = merge_intervals(
                masks.get((partition, candidate.record_id), ())
            )
            for start, end in subtract_intervals(
                [(candidate.local_start, candidate.local_end)], excluded,
            ):
                label = merge_annotation_labels(*labels)
                if candidate_id in connected_local:
                    label = merge_annotation_labels(
                        insertion_graphic_path(
                            placement, start, end, assembly_sources,
                        ), *local_global_labels[candidate_id],
                    )
                groups[(partition, candidate.record_id)].append({
                    "group_id": f"local:{candidate_id}:{start}-{end}",
                    "start": start, "end": end, "label": label,
                    "source_kind": "local_novel_annotated",
                    "origins": [f"local_query:{candidate_id}"],
                })

        for query in global_queries:
            if query.partition != partition or query.query_id in handled_globals:
                continue
            placement = global_best.get(query.query_id)
            excluded = list(masks.get((partition, query.record_id), ()))
            if query.query_id in connected_global and placement is not None:
                excluded.extend(placement.query_blocks)
            residuals = subtract_intervals(
                [(query.local_start, query.local_end)],
                merge_intervals(excluded),
            )
            for start, end in residuals:
                label = novel_labels.get(query.query_id, ".")
                if placement is not None:
                    label = merge_annotation_labels(
                        insertion_graphic_path(
                            placement, start, end, assembly_sources,
                        ) if query.query_id in connected_global
                        else complete_graphic_query_path(
                            placement.liftover_path, end - start,
                        ),
                        label,
                    )
                elif label == ".":
                    label = f"{end - start}I"
                groups[(partition, query.record_id)].append({
                    "group_id": f"global:{query.query_id}:{start}-{end}",
                    "start": start, "end": end, "label": label,
                    "source_kind": "global_novel_query",
                    "score_validated": query.query_id in score_validated_queries,
                    "origins": [
                        f"global_query:{query.query_id}",
                        *(
                            [f"global_minset:{query.query_id}"]
                            if query.query_id in novel_labels else []
                        ),
                    ],
                })

        provisional: List[Tuple[int, str, str, str, str, int, int, str, str]] = []
        for row in original:
            provisional.append((
                0, "original", "original", "original_block", row.haplotype,
                row.source_start, row.source_end, row.source_contig,
                json.dumps({"segments": [], "origins": [
                    f"block_reference:{row.encoded_name}"
                ]}, separators=(",", ":")),
            ))
            if row.encoded_name in fixed:
                segment = fixed[row.encoded_name]
                if (row.haplotype, row.source_contig, row.source_start, row.source_end) != (
                        segment['source_haplotype'], segment['source_contig'],
                        segment['source_start'], segment['source_end']):
                    raise ValueError(f'{partition}: imported alternative coordinates changed')
                provisional.append((
                    1, "reference", "reference", FIXED_ROLE, row.haplotype,
                    row.source_start, row.source_end, row.source_contig,
                    json.dumps({"segments": [segment], "origins": [FIXED_ROLE]},
                               separators=(",", ":")),
                ))
        # The reference half of the partition's self-cleaned union is
        # authoritative. Emit every retained reference interval, including
        # pieces outside the original BED core. Coordinate gaps remain
        # separate paths and are never restored from another reference.
        reference_union: DefaultDict[
            Tuple[str, str], List[UniqueRow]
        ] = defaultdict(list)
        if include_reference_unique_spans and not fixed_only:
            for unique in unique_by_partition.get(partition, ()):
                if unique.region_class == "reference" and may_refine(unique):
                    reference_union[
                        (unique.haplotype, unique.source_contig)
                    ].append(unique)
        for (haplotype, contig), reference_rows in sorted(
            reference_union.items()
        ):
            union_span = ReferenceSpan(
                haplotype, contig,
                min(row.source_start for row in reference_rows),
                max(row.source_end for row in reference_rows),
                0.0, set(),
            )
            for left, right, segments, origins in (
                partition_cleaned_reference_chains(
                    reference_rows, union_span,
                )
            ):
                provisional.append((
                    1, "reference", "reference", "refined_reference",
                    haplotype, left, right, contig,
                    json.dumps({
                        "segments": segments,
                        "origins": sorted(origins),
                    }, separators=(",", ":")),
                ))

        partition_fasta = partition_sample_fasta(partitions[partition])
        with ExitStack() as readers:
            assembly_readers: Dict[str, IndexedFasta] = {}
            partition_reader: Optional[IndexedFasta] = None
            if assembly_sources is None:
                partition_reader = readers.enter_context(
                    IndexedFasta(partition_fasta)
                )
            else:
                for haplotype in sorted({
                    row.haplotype
                    for row in unique_by_partition.get(partition, ())
                }):
                    assembly = assembly_sources.get(haplotype)
                    if assembly is None:
                        # Backward compatibility for callers that provide a
                        # partial assembly catalog while retaining the legacy
                        # partition samples FASTA.
                        if partition_reader is None:
                            partition_reader = readers.enter_context(
                                IndexedFasta(partition_fasta)
                            )
                        continue
                    assembly_readers[haplotype] = readers.enter_context(
                        IndexedFasta(assembly.fasta_path, assembly.fai_path)
                    )

            def fetch_unique(
                unique: UniqueRow, local_start: int, local_end: int,
            ) -> str:
                """Fetch a local unique span without requiring samples FASTA."""
                if assembly_sources is None:
                    assert partition_reader is not None
                    return partition_reader.fetch(
                        unique.record_id, local_start, local_end,
                    )
                reader = assembly_readers.get(unique.haplotype)
                if reader is None:
                    assert partition_reader is not None
                    return partition_reader.fetch(
                        unique.record_id, local_start, local_end,
                    )
                source_start = (
                    unique.source_start + local_start - unique.local_start
                )
                source_end = source_start + local_end - local_start
                if unique.source_contig not in reader.index:
                    raise KeyError(
                        f"{partition}: source contig {unique.source_contig!r} "
                        f"is absent from {reader.path}"
                    )
                return reader.fetch(
                    unique.source_contig, source_start, source_end,
                )

            for unique in sorted(
                unique_by_partition.get(partition, ()),
                key=lambda row: (row.record_id, row.local_start, row.local_end),
            ):
                group_rows = [
                    group for group in groups.get((partition, unique.record_id), ())
                    if unique.local_start <= int(group["start"])
                    and int(group["end"]) <= unique.local_end
                ]
                for group in group_rows:
                    start, end = int(group["start"]), int(group["end"])
                    if not bool(group.get("score_validated", False)):
                        sequence = fetch_unique(unique, start, end)
                        if novelty_score(
                            sequence, MASKED_BASE_WEIGHT,
                        ) < min_score:
                            continue
                    source_start = unique.source_start + start - unique.local_start
                    segment = {
                        "record": unique.record_id, "start": start,
                        "end": end, "strand": "+",
                    }
                    role = (
                        "novel" if group["source_kind"] == "global_novel_query"
                        else "alternative"
                    )
                    provisional.append((
                        3 if role == "novel" else 2, role,
                        str(group["label"]), str(group["source_kind"]),
                        unique.haplotype, source_start,
                        source_start + end - start, unique.source_contig,
                        json.dumps({
                            "segments": [segment],
                            "origins": list(group["origins"]),
                        }, separators=(",", ":")),
                    ))
                removed = [
                    (int(group["start"]), int(group["end"]))
                    for group in group_rows
                ] + masks.get((partition, unique.record_id), [])
                if unique.region_class == "reference":
                    removed.append((unique.local_start, unique.local_end))
                if unique.haplotype in suppressed_unique_haplotypes:
                    removed.append((unique.local_start, unique.local_end))
                # True-reference unique sequence is already represented by an
                # expanded reference component and must not reappear.
                for span in connected:
                    if (
                        unique.haplotype == span.haplotype
                        and unique.source_contig == span.contig
                    ):
                        left = max(unique.source_start, span.start)
                        right = min(unique.source_end, span.end)
                        if right > left:
                            removed.append((
                                unique.local_start + left - unique.source_start,
                                unique.local_start + right - unique.source_start,
                            ))
                for start, end in subtract_intervals(
                    [(unique.local_start, unique.local_end)],
                    merge_intervals(removed),
                ):
                    sequence = fetch_unique(unique, start, end)
                    if novelty_score(sequence, MASKED_BASE_WEIGHT) < min_score:
                        continue
                    source_start = unique.source_start + start - unique.local_start
                    provisional.append((
                        2, "alternative", ".", "query_residual",
                        unique.haplotype, source_start,
                        source_start + end - start, unique.source_contig,
                        json.dumps({
                            "segments": [{
                                "record": unique.record_id, "start": start,
                                "end": end, "strand": "+",
                            }],
                            "origins": [f"unique:{unique.region_name}"],
                        }, separators=(",", ":")),
                    ))

        template_spans: DefaultDict[Tuple[str, str], List[Interval]] = (
            defaultdict(list)
        )
        for item in provisional:
            priority, _role, _label, _kind, hap, start, end, contig, _payload = item
            if priority <= 1:
                template_spans[(hap, contig)].append((start, end))
        template_spans = defaultdict(list, {
            key: merge_intervals(intervals)
            for key, intervals in template_spans.items()
        })
        for item in provisional:
            priority, role, _label, kind, hap, start, end, contig, _payload = item
            if priority <= 1:
                continue
            overlap = next((
                (left, right)
                for left, right in template_spans.get((hap, contig), ())
                if start < right and left < end
            ), None)
            if overlap is not None:
                raise ValueError(
                    f"{partition}: {role}/{kind} source interval "
                    f"{hap}:{contig}:{start}-{end} overlaps selected "
                    f"template {overlap[0]}-{overlap[1]}"
                )

        provisional.sort(key=lambda item: (
            item[0], item[4], item[7], item[5], item[6], item[3], item[8],
        ))
        output: List[PathRow] = []
        for order, item in enumerate(provisional, 1):
            _priority, role, label, source_kind, hap, start, end, contig, payload_text = item
            payload = json.loads(payload_text)
            output.append(PathRow(
                partition, order, f"{partition}_{role}_{order:06d}",
                role, label, source_kind, hap, contig, start, end,
                (payload['segments'][0]['source_strand'] if source_kind == FIXED_ROLE else '+'),
                (json.dumps(payload['segments'], separators=(',', ':'))
                 if source_kind == FIXED_ROLE else compact_segments(payload["segments"])),
                ",".join(payload["origins"]),
            ))
        return partition, output

    outputs: Dict[str, List[PathRow]] = {}
    partition_names = sorted(partitions)
    process_context = mp_context()
    use_fork_processes = (
        jobs > 1
        and process_context.get_start_method() == "fork"
        and covered_records is None
    )
    LOG.info(
        "Building adjusted path rows for %d independent partitions with %d "
        "%s workers",
        len(partition_names), max(1, jobs),
        "forked process" if use_fork_processes else "thread",
    )
    progress_step = max(1, len(partition_names) // 20)
    global _FORKED_PARTITION_BUILDER
    if use_fork_processes:
        _FORKED_PARTITION_BUILDER = build_partition
        result_iterator = _bounded_process_map(
            _run_forked_partition_builder,
            partition_names,
            max(1, jobs),
        )
    else:
        result_iterator = _bounded_thread_map(
            build_partition, partition_names, max(1, jobs),
        )
    try:
        for completed, (partition, rows) in enumerate(result_iterator, 1):
            outputs[partition] = rows
            if (
                completed == len(partition_names)
                or completed % progress_step == 0
            ):
                LOG.info(
                    "Adjusted path construction progress: %d/%d partitions",
                    completed, len(partition_names),
                )
    finally:
        if use_fork_processes:
            _FORKED_PARTITION_BUILDER = None
    return [row for partition in sorted(outputs) for row in outputs[partition]]


def _build_adjusted_paths_legacy(
    partitions: Dict[str, PartitionFiles],
    unique_rows: Sequence[UniqueRow],
    used_rows: Sequence[UsedReferenceRow],
    local_candidates: Sequence[MappingCandidate],
    local_best: Dict[str, MappedPlacement],
    local_masks: Dict[Tuple[str, str], List[Interval]],
    global_queries: Sequence[NovelQuery],
    global_best: Dict[str, GlobalPlacement],
    reference_records: Dict[str, ReferenceRecord],
    main_reference: str,
    merge_gap: int,
    min_score: float,
    assembly_sources: Optional[
        Mapping[str, AssemblyFastaSource]
    ] = None,
    jobs: int = 1,
    small_novel_annotations: Optional[Mapping[str, str]] = None,
) -> List[PathRow]:
    """Classify original, reference, alternative, and novel partition paths."""
    del local_masks, reference_records, main_reference, merge_gap
    del assembly_sources
    novel_labels = small_novel_annotations or {}
    used_by_partition: DefaultDict[str, List[UsedReferenceRow]] = defaultdict(list)
    for row in used_rows:
        used_by_partition[row.partition].append(row)
    for partition in partitions:
        blocks = used_by_partition.get(partition, ())
        if not blocks:
            raise ValueError(
                f"partition {partition}: expected at least one "
                "block-reference interval, found 0"
            )

    candidate_by_id = {row.candidate_id: row for row in local_candidates}
    placement_partitions: Set[str] = set()
    for candidate_id, placement in local_best.items():
        candidate = candidate_by_id[candidate_id]
        if candidate.novel_partition not in partitions:
            raise ValueError(
                f"{candidate_id}: unknown novel partition "
                f"{candidate.novel_partition!r}"
        )
        placement_partitions.add(candidate.novel_partition)
    query_by_id = {row.query_id: row for row in global_queries}
    for query_id, placement in global_best.items():
        query = query_by_id[query_id]
        if placement.partition != query.partition:
            raise ValueError(
                f"{query_id}: placement partition {placement.partition!r} "
                f"does not match query partition {query.partition!r}"
            )
        if query.partition not in partitions:
            raise ValueError(
                f"{query_id}: unknown query partition {query.partition!r}"
        )
        placement_partitions.add(query.partition)
    LOG.info(
        "Partition routing: %d/%d partitions contain mapped novels; "
        "%d use the no-mapping pass-through path",
        len(placement_partitions),
        len(partitions),
        len(partitions) - len(placement_partitions),
    )
    LOG.info(
        "Path policy: original blocks are metadata; retained paths that "
        "overlap an original become references, local novels remain "
        "alternatives, and global novels remain novels",
    )

    # Local candidates remain alternatives regardless of their placement.
    # Global queries remain novel paths. A promoted global query contained by
    # one local candidate is represented by that complete local alternative.
    intact_groups: DefaultDict[
        Tuple[str, str], List[Dict[str, object]]
    ] = defaultdict(list)
    for candidate_id, placement in local_best.items():
        candidate = candidate_by_id[candidate_id]
        intact_groups[(
            candidate.novel_partition, candidate.record_id,
        )].append({
            "group_id": f"local:{candidate_id}",
            "start": candidate.local_start,
            "end": candidate.local_end,
            "label": placement.liftover_path or ".",
            "source_kind": "local_novel_annotated",
            "origins": [
                f"local_query:{candidate_id}",
                f"local_liftover:{candidate_id}",
            ],
            "global_ids": set(),
        })

    for query in global_queries:
        key = (query.partition, query.record_id)
        local_container = next((
            group for group in intact_groups.get(key, ())
            if (
                str(group["group_id"]).startswith("local:")
                and int(group["start"]) <= query.local_start
                and query.local_end <= int(group["end"])
            )
        ), None)
        if local_container is not None:
            local_container["origins"].append(  # type: ignore[union-attr]
                f"global_query:{query.query_id}"
            )
            local_container["global_ids"].add(query.query_id)  # type: ignore[union-attr]
            novel_label = novel_labels.get(query.query_id)
            if novel_label:
                local_container["label"] = merge_annotation_labels(
                    str(local_container["label"]), novel_label,
                )
                local_container["origins"].append(  # type: ignore[union-attr]
                    f"global_minset:{query.query_id}"
                )
            continue
        exact = next((
            group for group in intact_groups.get(key, ())
            if (
                int(group["start"]) == query.local_start
                and int(group["end"]) == query.local_end
                and str(group["group_id"]).startswith("global:")
            )
        ), None)
        placement = global_best.get(query.query_id)
        novel_label = novel_labels.get(query.query_id)
        if exact is None:
            exact = {
                "group_id": f"global:{query.query_id}",
                "start": query.local_start,
                "end": query.local_end,
                "label": merge_annotation_labels(
                    (
                        placement.liftover_path
                        if placement is not None else None
                    ),
                    novel_label,
                ),
                "source_kind": "global_novel_query",
                "origins": [],
                "global_ids": set(),
            }
            intact_groups[key].append(exact)
        exact["origins"].append(  # type: ignore[union-attr]
            f"global_query:{query.query_id}"
        )
        exact["global_ids"].add(query.query_id)  # type: ignore[union-attr]
        if (
            placement is not None
            and placement.liftover_path
            and placement.liftover_path != "."
        ):
            exact["label"] = merge_annotation_labels(
                str(exact["label"]), placement.liftover_path,
            )
            exact["origins"].append(  # type: ignore[union-attr]
                f"global_liftover:{query.query_id}"
            )
        if novel_label:
            exact["label"] = merge_annotation_labels(
                str(exact["label"]), novel_label,
            )
            exact["origins"].append(  # type: ignore[union-attr]
                f"global_minset:{query.query_id}"
            )

    for key, groups in intact_groups.items():
        prior: Optional[Dict[str, object]] = None
        for group in sorted(
            groups, key=lambda row: (int(row["start"]), int(row["end"])),
        ):
            if prior is not None and int(group["start"]) < int(prior["end"]):
                raise ValueError(
                    f"{key[0]}:{key[1]}: overlapping intact alternatives "
                    f"{prior['start']}-{prior['end']} and "
                    f"{group['start']}-{group['end']}"
                )
            prior = group

    unique_by_partition: DefaultDict[str, List[UniqueRow]] = defaultdict(list)
    for row in unique_rows:
        unique_by_partition[row.partition].append(row)

    def build_partition(partition: str):
        partition_output: List[PathRow] = []
        partition_emitted_groups: Set[str] = set()
        partition_handled_globals: Set[str] = set()
        provisional: List[Tuple[int, str, str, str, str, int, int, str, str]] = []
        # tuple priority, role, label, source_kind, haplotype, start, end,
        # contig, internal segment/origin payload.
        original = used_by_partition.get(partition, ())
        for row in original:
            payload = json.dumps({
                "segments": [],
                "origins": [f"block_reference:{row.encoded_name}"],
            }, separators=(",", ":"))
            provisional.append((
                0,
                "original",
                "original",
                "original_block",
                row.haplotype,
                row.source_start,
                row.source_end,
                row.source_contig,
                payload,
            ))

        # unique.bed coordinates and PARTITION.header refer to assembly
        # records in PARTITION_samples.fasta, not the block-template FASTA.
        partition_fasta = partition_sample_fasta(partitions[partition])
        with IndexedFasta(partition_fasta) as source:
            for unique in sorted(
                unique_by_partition.get(partition, ()),
                key=lambda row: (
                    row.record_id, row.local_start, row.local_end,
                ),
            ):
                retained = [(unique.local_start, unique.local_end)]
                contained_groups: List[Dict[str, object]] = []
                for group in intact_groups.get(
                    (partition, unique.record_id), (),
                ):
                    group_start = int(group["start"])
                    group_end = int(group["end"])
                    if (
                        group_start < unique.local_start
                        or group_end > unique.local_end
                    ):
                        continue
                    contained_groups.append(group)

                for group in contained_groups:
                    local_start = int(group["start"])
                    local_end = int(group["end"])
                    source_start = (
                        unique.source_start
                        + local_start - unique.local_start
                    )
                    segment = {
                        "fasta": os.path.abspath(partition_fasta),
                        "record": unique.record_id,
                        "start": local_start,
                        "end": local_end,
                        "strand": "+",
                        "source_start": source_start,
                        "source_end": source_start + local_end - local_start,
                    }
                    payload = json.dumps({
                        "segments": [segment],
                        "origins": list(group["origins"]) + [
                            f"unique:{unique.region_name}",
                            f"class:{unique.region_class}",
                        ],
                    }, separators=(",", ":"))
                    source_kind = str(group["source_kind"])
                    role = (
                        "novel"
                        if source_kind == "global_novel_query"
                        else "alternative"
                    )
                    provisional.append((
                        3 if role == "novel" else 2,
                        role, str(group["label"]), source_kind,
                        unique.haplotype, source_start,
                        source_start + local_end - local_start,
                        unique.source_contig, payload,
                    ))
                    partition_emitted_groups.add(str(group["group_id"]))
                    partition_handled_globals.update(group["global_ids"])

                intact_intervals = [
                    (int(group["start"]), int(group["end"]))
                    for group in contained_groups
                ]
                ordinary_retained: List[Interval] = []
                for start, end in retained:
                    ordinary_retained.extend(subtract_intervals(
                        [(start, end)], intact_intervals,
                    ))
                for local_start, local_end in ordinary_retained:
                    sequence = source.fetch(
                        unique.record_id, local_start, local_end,
                    )
                    if novelty_score(
                        sequence, MASKED_BASE_WEIGHT,
                    ) < min_score:
                        continue
                    source_start = (
                        unique.source_start
                        + local_start - unique.local_start
                    )
                    segment = {
                        "fasta": os.path.abspath(partition_fasta),
                        "record": unique.record_id,
                        "start": local_start,
                        "end": local_end,
                        "strand": "+",
                        "source_start": source_start,
                        "source_end": (
                            source_start + local_end - local_start
                        ),
                    }
                    payload = json.dumps({
                        "segments": [segment],
                        "origins": [
                            f"unique:{unique.region_name}",
                            f"class:{unique.region_class}",
                        ],
                    }, separators=(",", ":"))
                    overlaps_original = any(
                        block.haplotype == unique.haplotype
                        and block.source_contig == unique.source_contig
                        and source_start < block.source_end
                        and block.source_start < (
                            source_start + local_end - local_start
                        )
                        for block in original
                    )
                    role = "reference" if overlaps_original else "alternative"
                    provisional.append((
                        1 if role == "reference" else 2,
                        role,
                        role,
                        "partition_reference" if role == "reference"
                        else "query_residual",
                        unique.haplotype, source_start,
                        source_start + local_end - local_start,
                        unique.source_contig, payload,
                    ))

        provisional.sort(key=lambda item: (
            item[0], item[4], item[7], item[5], item[6], item[3], item[8],
        ))
        for order, item in enumerate(provisional, 1):
            (
                _priority, role, label, source_kind, haplotype,
                start, end, contig, payload_text,
            ) = item
            payload = json.loads(payload_text)
            partition_output.append(PathRow(
                partition,
                order,
                f"{partition}_{role}_{order:06d}",
                role,
                label,
                source_kind,
                haplotype,
                contig,
                start,
                end,
                "+",
                compact_segments(payload["segments"]),
                ",".join(payload["origins"]),
            ))
        return (
            partition,
            partition_output,
            partition_emitted_groups,
            partition_handled_globals,
        )

    ordered_partitions = sorted(partitions)
    partition_outputs: Dict[str, List[PathRow]] = {}
    emitted_group_ids: Set[str] = set()
    handled_global_ids: Set[str] = set()
    progress_step = max(1, len(ordered_partitions) // 20)
    LOG.info(
        "Building adjusted paths for %d partitions with %d workers",
        len(ordered_partitions), max(1, jobs),
    )
    for completed, result in enumerate(
        _bounded_thread_map(
            build_partition, ordered_partitions, max(1, jobs),
        ),
        1,
    ):
        (
            partition,
            partition_rows,
            partition_groups,
            partition_globals,
        ) = result
        partition_outputs[partition] = partition_rows
        emitted_group_ids.update(partition_groups)
        handled_global_ids.update(partition_globals)
        if (
            completed == len(ordered_partitions)
            or completed % progress_step == 0
        ):
            LOG.info(
                "Adjusted-path construction progress: %d/%d partitions",
                completed, len(ordered_partitions),
            )

    output = [
        row
        for partition in ordered_partitions
        for row in partition_outputs[partition]
    ]

    expected_groups = {
        str(group["group_id"])
        for groups in intact_groups.values()
        for group in groups
    }
    missing_groups = expected_groups - emitted_group_ids
    if missing_groups:
        raise ValueError(
            "intact alternatives were not contained by a unique interval: "
            + ",".join(sorted(missing_groups)[:20])
            + ("..." if len(missing_groups) > 20 else "")
        )
    missing_global = {
        query.query_id for query in global_queries
    } - handled_global_ids
    if missing_global:
        raise ValueError(
            "global queries were not retained as novel/local paths: "
            + ",".join(sorted(missing_global)[:20])
            + ("..." if len(missing_global) > 20 else "")
        )
    return output


def write_global_alignment_summary(
    queries: Sequence[NovelQuery],
    global_placements: Dict[str, GlobalPlacement],
    output_dir: str,
) -> None:
    fields = (
        "query_id", "query_kind", "partition", "status",
        "ref_haplotype", "ref_contig", "ref_start", "ref_end",
        "strand", "alignment_score", "aligned_bases", "identity",
        "source", "target_record", "query_blocks", "target_blocks",
        "paired_blocks", "graphic_cigar",
    )
    write_tsv(
        os.path.join(output_dir, "global_novel_alignment_summary.tsv"),
        fields,
        (
            (
                query.query_id, "global_novel", query.partition,
                "mapped" if query.query_id in global_placements else "unmapped",
                *(
                    (
                        placement.ref_haplotype, placement.ref_contig,
                        placement.ref_start, placement.ref_end,
                        placement.strand, placement.alignment_score,
                        placement.aligned_bases, placement.identity,
                        placement.source, placement.target_record,
                        json.dumps(placement.query_blocks, separators=(",", ":")),
                        json.dumps(placement.target_blocks, separators=(",", ":")),
                        json.dumps(placement.paired_blocks, separators=(",", ":")),
                        placement.liftover_path,
                    )
                    if (placement := global_placements.get(query.query_id)) is not None
                    else (".", ".", ".", ".", ".", 0, 0, 0, ".", ".",
                          "[]", "[]", "[]", ".")
                ),
            )
            for query in queries
        ),
    )


def write_global_mapping_outputs(
    queries: Sequence[NovelQuery],
    final_global: Sequence[NovelQuery],
    global_placements: Dict[str, GlobalPlacement],
    output_dir: str,
) -> None:
    query_fields = [field.name for field in dataclasses.fields(NovelQuery)]
    write_tsv(
        os.path.join(output_dir, "all_global_queries.tsv"),
        query_fields,
        (dataclasses.astuple(row) for row in queries),
    )
    write_tsv(
        os.path.join(output_dir, "final_global_novel.tsv"),
        query_fields,
        (dataclasses.astuple(row) for row in final_global),
    )
    write_tsv(
        os.path.join(output_dir, "all_global_reference_placements.tsv"),
        [field.name for field in dataclasses.fields(GlobalPlacement)],
        (
            (
                *dataclasses.astuple(row)[:-4],
                json.dumps(row.query_blocks, separators=(",", ":")),
                json.dumps(row.target_blocks, separators=(",", ":")),
                json.dumps(row.paired_blocks, separators=(",", ":")),
                row.liftover_path,
            )
            for row in sorted(
                global_placements.values(),
                key=lambda row: row.query_id,
            )
        ),
    )
    write_global_alignment_summary(queries, global_placements, output_dir)


def write_adjusted_path_outputs(
    rows: Sequence[PathRow],
    output_dir: str,
    jobs: int = 1,
    partition_sources: Optional[Mapping[str, str]] = None,
) -> None:
    fields = [field.name for field in dataclasses.fields(PathRow)]
    write_tsv(
        os.path.join(output_dir, "all_adjusted_paths.tsv"),
        fields,
        (dataclasses.astuple(row) for row in rows),
    )
    if partition_sources is not None:
        write_tsv(
            os.path.join(output_dir, "partition_sources.tsv"),
            ("partition", "samples_fasta"),
            (
                (partition, os.path.abspath(path))
                for partition, path in sorted(partition_sources.items())
            ),
        )
    by_partition: DefaultDict[str, List[PathRow]] = defaultdict(list)
    for row in rows:
        by_partition[row.partition].append(row)
    def write_partition(item):
        partition, members = item
        directory = os.path.join(output_dir, partition)
        Path(directory).mkdir(parents=True, exist_ok=True)
        write_tsv(
            os.path.join(directory, f"{partition}.adjusted.tsv"),
            fields,
            (dataclasses.astuple(row) for row in members),
        )
        return partition

    items = sorted(by_partition.items())
    progress_step = max(1, len(items) // 20)
    LOG.info(
        "Writing %d adjusted partition manifests with %d workers",
        len(items), max(1, jobs),
    )
    for completed, _partition in enumerate(
        _bounded_thread_map(write_partition, items, max(1, jobs)), 1,
    ):
        if completed == len(items) or completed % progress_step == 0:
            LOG.info(
                "Adjusted-manifest write progress: %d/%d partitions",
                completed, len(items),
            )


def write_outputs(
    rows: Sequence[PathRow],
    final_global: Sequence[NovelQuery],
    global_placements: Dict[str, GlobalPlacement],
    output_dir: str,
    global_queries: Sequence[NovelQuery] = (),
    jobs: int = 1,
    partition_sources: Optional[Mapping[str, str]] = None,
) -> None:
    """Compatibility writer retaining the pre-split adjusted-directory files."""
    write_adjusted_path_outputs(
        rows, output_dir, jobs, partition_sources,
    )
    write_global_mapping_outputs(
        global_queries, final_global, global_placements, output_dir,
    )


def parse_args(argv: Optional[Sequence[str]] = None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=(
            "Classify completed partition paths as original metadata, "
            "materialized references overlapping an original, alternatives, "
            "or global novels, while retaining mapping annotations."
        )
    )
    parser.add_argument("-a", "--analysis-dir", required=True)
    parser.add_argument(
        "-l", "--local-mapping-dir", required=True,
        help="output of map_partition_local_novel.py",
    )
    parser.add_argument(
        "-g", "--global-mapping-dir", required=True,
        help="output of map_global_novel.py",
    )
    parser.add_argument(
        "-s", "--small-novel-mapping-dir",
        help=(
            "optional output of map_small_novel_minset.py; its large/small "
            "novel-minset paths are added to alternative FASTA labels"
        ),
    )
    parser.add_argument(
        "--alignment-summary",
        action="append",
        default=[],
        metavar="TSV",
        help=(
            "repeatable local/global alignment-summary TSV; annotations from "
            "any number of supplied files are combined in command-line order"
        ),
    )
    parser.add_argument(
        "-r", "--main-reference", required=True,
        help="cleaned cohort main_chroms.fa",
    )
    parser.add_argument("-o", "--output-dir", required=True)
    parser.add_argument(
        "--merge-gap", type=int, default=5000,
        help="compatibility option retained in the output signature",
    )
    parser.add_argument(
        "--min-score", type=int, default=100,
        help="minimum unmasked + 0.3*masked score for alternative residuals",
    )
    parser.add_argument(
        "-j", "--jobs", type=int, default=8,
        help=(
            "parallel workers for local-placement annotation, per-partition "
            "path construction, and manifest writing (default: 8)"
        ),
    )
    parser.add_argument("--resume", action="store_true")
    parser.add_argument(
        "--log-level", choices=("DEBUG", "INFO", "WARNING", "ERROR"),
        default="INFO",
    )
    args = parser.parse_args(argv)
    if args.min_score < 1:
        parser.error("--min-score must be positive")
    if args.merge_gap < 0:
        parser.error("--merge-gap cannot be negative")
    if args.jobs < 1:
        parser.error("--jobs must be positive")
    return args


def run(args: argparse.Namespace) -> None:
    analysis_dir = os.path.abspath(args.analysis_dir)
    local_mapping_dir = os.path.abspath(args.local_mapping_dir)
    global_mapping_dir = os.path.abspath(args.global_mapping_dir)
    small_novel_mapping_dir = (
        os.path.abspath(args.small_novel_mapping_dir)
        if args.small_novel_mapping_dir else None
    )
    main_reference = os.path.abspath(args.main_reference)
    output_dir = os.path.abspath(args.output_dir)
    alignment_summary_paths = [
        os.path.abspath(path) for path in args.alignment_summary
    ]
    Path(output_dir).mkdir(parents=True, exist_ok=True)
    partition_inputs = os.path.join(analysis_dir, "partition_inputs.tsv")
    local_candidates_path = os.path.join(
        analysis_dir, "local_novel_candidates.tsv",
    )
    local_placements_path = os.path.join(
        local_mapping_dir, "all_local_novel_alignments.tsv",
    )
    local_reporting_path = os.path.join(
        local_mapping_dir, "all_local_novel_reporting_alignments.tsv",
    )
    global_queries_path = os.path.join(
        global_mapping_dir, "all_global_queries.tsv",
    )
    global_placements_path = os.path.join(
        global_mapping_dir, "all_global_reference_placements.tsv",
    )
    global_reporting_path = os.path.join(
        global_mapping_dir,
        "all_global_reference_reporting_alignments.tsv",
    )
    final_global_path = os.path.join(
        global_mapping_dir, "final_global_novel.tsv",
    )
    global_marker = os.path.join(
        global_mapping_dir, "_GLOBAL_NOVEL_MAPPING_SUCCESS.json",
    )
    small_annotations_path = None
    small_marker = None
    if small_novel_mapping_dir is not None:
        small_annotations_path = os.path.join(
            small_novel_mapping_dir, "all_small_novel_annotations.tsv",
        )
        small_marker = os.path.join(
            small_novel_mapping_dir, "_SMALL_NOVEL_MINSET_SUCCESS.json",
        )
    required = (
        partition_inputs,
        os.path.join(analysis_dir, "all_unique_regions.bed"),
        os.path.join(analysis_dir, "all_used_novel_loci.bed"),
        local_candidates_path,
        local_placements_path,
        global_queries_path,
        global_placements_path,
        global_reporting_path,
        final_global_path,
        global_marker,
        main_reference,
        *((small_annotations_path, small_marker)
          if small_annotations_path is not None and small_marker is not None
          else ()),
        *alignment_summary_paths,
    )
    for path in required:
        if not os.path.isfile(path):
            raise FileNotFoundError(path)
    if os.path.isfile(local_reporting_path):
        required = (*required, local_reporting_path)
    with open(global_marker, "rt") as handle:
        global_metadata = json.load(handle)
    expected_global_version = (
        "global-novel-mapping-v5-graphic-annotation"
    )
    if global_metadata.get("version") != expected_global_version:
        raise ValueError(
            "global mapping predates annotation-only graphic liftover; "
            "rerun the updated map_global_novel.py"
        )
    if global_metadata.get("reporting_resolution") != (
        REPORTING_RESOLUTION_VERSION
    ):
        raise ValueError(
            "global mapping has no current multi-locus reporting table; "
            "replay the saved alignments with map_global_novel.py "
            "--reuse-alignments-only"
        )
    if small_marker is not None:
        with open(small_marker, "rt") as handle:
            small_metadata = json.load(handle)
        if small_metadata.get("version") != (
            "small-novel-minset-v3-audit-intermediates"
        ):
            raise ValueError(
                "small-novel mapping predates the supported minset format; "
                "rerun map_small_novel_minset.py"
            )
    query_paths = global_metadata.get("query_paths")
    if not query_paths or not os.path.isfile(query_paths):
        raise ValueError(
            "global mapping does not record a usable query_paths file; "
            "rerun the updated map_global_novel.py"
        )
    query_paths = os.path.abspath(query_paths)
    required = (*required, query_paths)
    if float(global_metadata.get("min_score", -1)) != args.min_score:
        raise ValueError(
            f"--min-score {args.min_score} does not match global mapping "
            f"min_score {global_metadata.get('min_score')!r}"
        )
    signature_payload = {
        "version": 11,
        "files": [_fingerprint(path) for path in required],
        "merge_gap": args.merge_gap,
        "min_score": args.min_score,
    }
    signature = hashlib.sha256(
        json.dumps(signature_payload, sort_keys=True).encode()
    ).hexdigest()
    marker = os.path.join(output_dir, "_ADJUSTED_PATHS_SUCCESS.json")
    if args.resume and os.path.isfile(marker):
        try:
            with open(marker, "rt") as handle:
                saved = json.load(handle)
            if (
                saved.get("signature") == signature
                and os.path.isfile(os.path.join(
                    output_dir, "all_adjusted_paths.tsv",
                ))
                and os.path.isfile(os.path.join(
                    output_dir, "partition_sources.tsv",
                ))
                and os.path.isfile(os.path.join(
                    output_dir, "local_novel_reporting_alignments.tsv",
                ))
                and os.path.isfile(os.path.join(
                    output_dir, "reference_covered_mapped_novel.tsv",
                ))
            ):
                LOG.info("RESUME: adjusted paths are complete in %s", output_dir)
                return
        except (OSError, ValueError, AttributeError):
            pass

    partitions = read_partition_inputs(partition_inputs)
    unique_rows = read_unique_rows(analysis_dir, partitions)
    used_rows = read_used_reference_rows(analysis_dir)
    local_candidates = read_candidates(local_candidates_path)
    local_placements = read_placements([local_placements_path])
    reference_records = parse_main_reference(main_reference)
    assembly_sources = read_assembly_sources(query_paths)
    local_best, local_masks = choose_local_placements(
        local_candidates, local_placements, reference_records,
    )
    if os.path.isfile(local_reporting_path):
        local_reporting_components = read_reporting_components((
            local_reporting_path,
        ))
        local_reporting_labels = annotate_local_reporting_components(
            local_candidates, local_reporting_components,
            assembly_sources, args.jobs,
        )
        write_reporting_components(
            os.path.join(
                output_dir, "local_novel_reporting_alignments.tsv",
            ),
            (
                dataclasses.replace(
                    component,
                    liftover_path=local_reporting_labels.get(
                        component.query_id, ".",
                    ),
                )
                for component in local_reporting_components
            ),
        )
    else:
        local_reporting_components = []
        local_reporting_labels = {}
        LOG.warning(
            "No corrected local multi-locus report at %s; retaining legacy "
            "single-placement annotations. Replay map_partition_local_novel.py "
            "--reuse-alignments-only to correct them.",
            local_reporting_path,
        )
        write_reporting_components(
            os.path.join(
                output_dir, "local_novel_reporting_alignments.tsv",
            ),
            (),
        )
    # Keep the original one-placement cleanup/masking decision.  Only replace
    # its public annotation with the ordered all-target query path.  A legacy
    # annotation is computed solely for an old placement that has no matching
    # reporting component (normally zero rows).
    missing_local_labels = {
        candidate_id: placement
        for candidate_id, placement in local_best.items()
        if candidate_id not in local_reporting_labels
    }
    if missing_local_labels:
        missing_local_labels = annotate_local_placements(
            local_candidates, missing_local_labels,
            assembly_sources, args.jobs,
        )
    local_best = {
        candidate_id: dataclasses.replace(
            placement,
            liftover_path=local_reporting_labels.get(
                candidate_id,
                missing_local_labels.get(candidate_id, placement).liftover_path,
            ),
        )
        for candidate_id, placement in local_best.items()
    }
    global_queries = read_novel_queries(global_queries_path)
    final_global = read_novel_queries(final_global_path)
    global_best = read_global_placements(global_placements_path)
    global_reporting_labels = reporting_liftover_labels(
        read_reporting_components((global_reporting_path,))
    )
    global_best = {
        query_id: dataclasses.replace(
            placement,
            liftover_path=global_reporting_labels.get(
                query_id, placement.liftover_path,
            ),
        )
        for query_id, placement in global_best.items()
    }
    small_novel_annotations = (
        read_small_novel_annotations(small_annotations_path)
        if small_annotations_path is not None else {}
    )
    additional_annotations = read_alignment_summary_annotations(
        alignment_summary_paths
    )
    unknown_placements = set(global_best) - {
        query.query_id for query in global_queries
    }
    if unknown_placements:
        raise ValueError(
            "global placements reference unknown queries: "
            + ",".join(sorted(unknown_placements))
        )
    unknown_annotations = set(small_novel_annotations) - {
        query.query_id for query in global_queries
    }
    if unknown_annotations:
        raise ValueError(
            "small-novel annotations reference unknown global queries: "
            + ",".join(sorted(unknown_annotations)[:20])
        )
    LOG.info(
        "Loaded completed global mappings for %d/%d queries",
        len(global_best), len(global_queries),
    )
    covered_records: List[Dict[str, object]] = []
    rows = build_adjusted_paths(
        partitions,
        unique_rows,
        used_rows,
        local_candidates,
        local_best,
        local_masks,
        global_queries,
        global_best,
        reference_records,
        main_reference,
        args.merge_gap,
        args.min_score,
        assembly_sources,
        args.jobs,
        small_novel_annotations,
        additional_annotations,
        covered_records,
    )
    covered_fields = (
        "partition", "query_kind", "query_id", "record_id",
        "haplotype", "source_contig", "source_start", "source_end",
        "local_start", "local_end", "ref_haplotype", "ref_contig",
        "ref_start", "ref_end", "strand", "reason",
    )
    write_tsv(
        os.path.join(output_dir, "reference_covered_mapped_novel.tsv"),
        covered_fields,
        (
            tuple(record[field] for field in covered_fields)
            for record in sorted(covered_records, key=lambda row: (
                str(row["partition"]), str(row["query_kind"]),
                str(row["query_id"]), int(row["local_start"]),
                int(row["local_end"]),
            ))
        ),
    )
    write_outputs(
        rows, final_global, global_best, output_dir, global_queries,
        args.jobs,
        {
            partition: partition_sample_fasta(files)
            for partition, files in partitions.items()
        },
    )
    temporary = marker + ".tmp"
    with open(temporary, "wt") as out:
        json.dump({
            "version": "adjusted-partition-paths-v11-reference-refinement",
            "signature": signature,
            "partitions": len(partitions),
            "paths": len(rows),
            "global_queries": len(global_queries),
            "global_placements": len(global_best),
            "local_reporting_components": len(local_reporting_components),
            "global_reporting_queries": len(global_reporting_labels),
            "final_global_queries": len(final_global),
            "global_mapping_dir": global_mapping_dir,
            "small_novel_mapping_dir": small_novel_mapping_dir,
            "alignment_summaries": alignment_summary_paths,
            "reference_covered_mapped_novel_intervals": len(covered_records),
            "main_reference": main_reference,
        }, out, sort_keys=True)
        out.write("\n")
    os.replace(temporary, marker)
    LOG.info(
        "Wrote %d adjusted paths for %d partitions under %s",
        len(rows), len(partitions), output_dir,
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
    except Exception as error:
        LOG.error("%s", error)
        if args.log_level == "DEBUG":
            raise
        return 1
    return 0


if __name__ == "__main__":
    sys.exit(main())
