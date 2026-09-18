#!/usr/bin/env python3
"""Clean saved local-novel candidates against selected true-reference blocks.

This is stage two of partition novelty analysis. It consumes the flat files
written by ``analyze_partition_coordinates.py``. ``PARTITION.header`` supplies
all query and Ref coordinates. Query windows are fetched once per assembly via
its existing FAI, required Ref contigs are loaded into shared RAM once, and
small target FASTAs are written from that cache. Each candidate/eligible-
partition assignment is grouped by partition. Each partition is one independent
parallel multi-query job using one Minimap2 pass, one Winnowmap pass, and
finally one BLASTN pass. Local-block FASTAs are never read or indexed.
"""
from __future__ import annotations

import argparse
import csv
import dataclasses
import heapq
import hashlib
import json
import logging
import os
import shutil
import subprocess
import sys
import threading
from collections import defaultdict
from concurrent.futures import FIRST_COMPLETED, ThreadPoolExecutor, wait
from pathlib import Path
from typing import (
    Dict,
    FrozenSet,
    Iterable,
    List,
    Mapping,
    Optional,
    Sequence,
    Tuple,
)

from analyze_partition_coordinates import (
    CoordinateCatalog,
    parse_reference_haplotypes,
    parse_used_novel_regions,
)
from minsetref_align import (
    AlignmentHit,
    build_blast_db,
    collect_alignment_hits,
    collect_blastn_hits,
    ensure_executable,
)
from minsetref_core import (
    DNA_WRAP,
    IndexedFasta,
    iter_n_free_segments,
    mp_context,
    wrap_fasta,
)
from minsetref_light import (
    LIGHT_WINNOW_PARAMS,
    MASKED_BASE_WEIGHT,
    MINIMAP_PARAMS,
)
from minsetref_segments import (
    complement_intervals,
    find_unaligned_segments,
    merge_intervals,
    novelty_score,
    subtract_intervals,
)
from minsetref_multilocus import (
    PlacementEvidence,
    ResolvedPlacementComponent,
    resolve_all_evidence,
)
from summarize_partition_novelty import (
    HeaderRecord,
    PartitionFiles,
    parse_header_records,
    partition_header_path,
    write_bed_table,
    write_tsv,
)


LOG = logging.getLogger("map_partition_local_novel")
Interval = Tuple[int, int]
PLACEMENT_RESOLUTION_VERSION = "query-score-sweep-v3-thread-safe"
REPORTING_RESOLUTION_VERSION = "multilocus-query-sweep-v1"


def _allow_large_csv_fields() -> None:
    """Lift Python's 128-KiB default for JSON-rich placement TSV fields."""
    limit = sys.maxsize
    while True:
        try:
            csv.field_size_limit(limit)
            return
        except OverflowError:
            limit //= 10


_allow_large_csv_fields()


@dataclasses.dataclass(frozen=True)
class MappingCandidate:
    candidate_id: str
    interval_id: str
    source_contig: str
    source_start: int
    source_end: int
    haplotype: str
    novel_partition: str
    lifted_partitions: FrozenSet[str]
    record_id: str
    local_start: int
    local_end: int
    source_fasta: str
    region_class: str


@dataclasses.dataclass(frozen=True)
class MappedPlacement:
    candidate_id: str
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
    query_blocks: str
    target_blocks: str
    paired_blocks: str = "[]"
    liftover_path: str = "."


@dataclasses.dataclass(frozen=True)
class BlockTarget:
    target_id: str
    input_record: str
    local_start: int
    local_end: int
    ref_haplotype: str
    ref_contig: str
    ref_start: int
    ref_end: int
    ref_strand: str = "+"


@dataclasses.dataclass(frozen=True)
class ResidualPiece:
    candidate_id: str
    local_start: int
    local_end: int


@dataclasses.dataclass(frozen=True)
class CycleQuery:
    query_id: str
    piece: ResidualPiece
    core_start: int
    core_end: int
    sequence: Optional[str]
    sequence_length: int


@dataclasses.dataclass
class PlacementAccumulator:
    target_id: str
    strand: str
    query_masks: List[Interval]
    target_blocks: List[Interval]
    aligned_query_blocks: List[Interval]
    paired_blocks: List["ScoredPairedBlock"]
    alignment_score: float
    identity_bases: float
    aligned_bases: int
    sources: set


@dataclasses.dataclass(frozen=True)
class ScoredPairedBlock:
    query_start: int
    query_end: int
    ref_start: int
    ref_end: int
    alignment_score: float
    alignment_bases: int
    identity: float
    source: str
    ordinal: int


REPORTING_FIELDS = (
    "candidate_id", "component_index", "mapped_partition",
    "ref_haplotype", "ref_contig", "ref_start", "ref_end", "strand",
    "target_record", "query_start", "query_end", "query_blocks",
    "target_blocks", "paired_blocks", "alignment_score", "aligned_bases",
    "identity", "source", "evidence_ordinal", "liftover_path",
)


@dataclasses.dataclass(frozen=True)
class PreparedBlockReference:
    mapped_partition: str
    target_fasta: str
    targets: Tuple[BlockTarget, ...]
    signature: str


@dataclasses.dataclass(frozen=True)
class AssemblyFastaSource:
    haplotype: str
    fasta_path: str
    fai_path: str


@dataclasses.dataclass(frozen=True)
class CachedCandidateSequence:
    """One candidate plus its reusable source-side alignment flanks."""

    local_start: int
    local_end: int
    sequence: str


def _bounded_thread_map(function, items: Iterable, workers: int):
    """Yield results while retaining at most ``workers`` submitted futures."""
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


def _parse_set(text: str) -> FrozenSet[str]:
    if not text or text == ".":
        return frozenset()
    return frozenset(item for item in text.split(",") if item and item != ".")


def read_partition_inputs(path: str) -> Dict[str, PartitionFiles]:
    output: Dict[str, PartitionFiles] = {}
    with open(path, "rt", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        expected = {"partition", "fasta", "unique_bed", "reference_bed"}
        if reader.fieldnames is None or not expected.issubset(reader.fieldnames):
            raise ValueError(f"{path}: expected columns {sorted(expected)}")
        for row in reader:
            partition = row["partition"]
            if partition in output:
                raise ValueError(f"{path}: duplicate partition {partition!r}")
            files = PartitionFiles(
                partition,
                row["fasta"],
                row["unique_bed"],
                row["reference_bed"],
            )
            expected_header = partition_header_path(files)
            saved_header = row.get("header") or expected_header
            accepted_headers = {
                os.path.abspath(expected_header),
                os.path.abspath(os.path.splitext(row["fasta"])[0] + ".header"),
                os.path.abspath(
                    os.path.splitext(row["fasta"])[0] + "_samples.header"
                ),
            }
            if os.path.abspath(saved_header) not in accepted_headers:
                raise ValueError(
                    f"{path}: header path for partition {partition!r} "
                    "does not match its FASTA basename"
                )
            if not os.path.isfile(saved_header):
                raise FileNotFoundError(saved_header)
            output[partition] = files
    if not output:
        raise ValueError(f"{path}: no partition inputs")
    return output


def read_assembly_sources(path: str) -> Dict[str, AssemblyFastaSource]:
    """Read NAME FASTA from the working directory, requiring adjacent indexes."""
    output: Dict[str, AssemblyFastaSource] = {}
    with open(path, "rt") as handle:
        for line_number, raw in enumerate(handle, 1):
            if not raw.strip() or raw.startswith("#"):
                continue
            fields = raw.split()
            if len(fields) != 2:
                raise ValueError(
                    f"{path}:{line_number}: expected exactly NAME FASTA; index must be FASTA.fai"
                )
            haplotype = fields[0]
            fasta_path = os.path.abspath(os.path.expanduser(fields[1]))
            fai_path = fasta_path + ".fai"
            if haplotype in output:
                raise ValueError(
                    f"{path}:{line_number}: duplicate haplotype {haplotype!r}"
                )
            for required in (fasta_path, fai_path):
                if not os.path.isfile(required):
                    raise FileNotFoundError(required)
            output[haplotype] = AssemblyFastaSource(
                haplotype, fasta_path, fai_path,
            )
    if not output:
        raise ValueError(f"{path}: no assembly sources")
    return output


def read_candidates(path: str) -> List[MappingCandidate]:
    output: List[MappingCandidate] = []
    seen = set()
    with open(path, "rt", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        required = {
            "candidate_id", "interval_id", "source_contig", "source_start",
            "source_end", "haplotype", "novel_partition",
            "lifted_partitions", "record_id", "local_start", "local_end",
            "source_fasta", "region_class",
        }
        if reader.fieldnames is None or not required.issubset(reader.fieldnames):
            raise ValueError(f"{path}: expected columns {sorted(required)}")
        for line_number, row in enumerate(reader, 2):
            candidate_id = row["candidate_id"]
            if candidate_id in seen:
                raise ValueError(
                    f"{path}:{line_number}: duplicate candidate {candidate_id}"
                )
            seen.add(candidate_id)
            source_start = int(row["source_start"])
            source_end = int(row["source_end"])
            local_start = int(row["local_start"])
            local_end = int(row["local_end"])
            if source_end <= source_start or local_end <= local_start:
                raise ValueError(
                    f"{path}:{line_number}: empty candidate interval"
                )
            if source_end - source_start != local_end - local_start:
                raise ValueError(
                    f"{path}:{line_number}: source/local candidate lengths "
                    "disagree"
                )
            output.append(MappingCandidate(
                candidate_id,
                row["interval_id"],
                row["source_contig"],
                source_start,
                source_end,
                row["haplotype"],
                row["novel_partition"],
                _parse_set(row["lifted_partitions"]),
                row["record_id"],
                local_start,
                local_end,
                row["source_fasta"],
                row["region_class"],
            ))
    return output


def load_candidate_header_records(
    candidates: Sequence[MappingCandidate],
    partition_inputs: Dict[str, PartitionFiles],
    jobs: int,
) -> Dict[Tuple[str, str], HeaderRecord]:
    """Load only candidate-used records from authoritative ``.header`` files."""
    needed: Dict[str, set] = defaultdict(set)
    for candidate in candidates:
        needed[candidate.novel_partition].add(candidate.record_id)

    def load_partition(item):
        partition, record_ids = item
        files = partition_inputs[partition]
        selected = {
            row.record_id: row
            for row in parse_header_records(files, record_ids)
        }
        missing = record_ids - set(selected)
        if missing:
            examples = ",".join(sorted(missing)[:10])
            raise ValueError(
                f"{partition_header_path(files)}: missing "
                f"{len(missing)} candidate records; examples: {examples}"
            )
        return partition, selected

    output: Dict[Tuple[str, str], HeaderRecord] = {}
    items = sorted(needed.items())
    progress_step = max(1, len(items) // 20)
    for completed, result in enumerate(_bounded_thread_map(
        load_partition, items, jobs,
    ), 1):
        partition, selected = result
        output.update(
            ((partition, record_id), record)
            for record_id, record in selected.items()
        )
        if completed == len(items) or completed % progress_step == 0:
            LOG.info(
                "Loaded candidate records from .header: %d/%d partitions",
                completed, len(items),
            )

    for candidate in candidates:
        record = output[(
            candidate.novel_partition, candidate.record_id,
        )]
        expected_start = record.source_start + candidate.local_start
        expected_end = record.source_start + candidate.local_end
        if (
            record.haplotype != candidate.haplotype
            or record.source_contig != candidate.source_contig
            or expected_start != candidate.source_start
            or expected_end != candidate.source_end
        ):
            raise ValueError(
                f"{candidate.candidate_id}: saved source coordinates disagree "
                f"with {candidate.novel_partition}.header record "
                f"{candidate.record_id}: expected "
                f"{record.haplotype}:{record.source_contig}:"
                f"{expected_start}-{expected_end}, observed "
                f"{candidate.haplotype}:{candidate.source_contig}:"
                f"{candidate.source_start}-{candidate.source_end}"
            )
    return output


def load_candidate_sequences(
    candidates: Sequence[MappingCandidate],
    source_records: Dict[Tuple[str, str], HeaderRecord],
    assembly_sources: Dict[str, AssemblyFastaSource],
    anchor: int,
    jobs: int,
) -> Dict[str, CachedCandidateSequence]:
    """Fetch every query window once, opening each assembly only once."""
    by_haplotype: Dict[str, List[MappingCandidate]] = defaultdict(list)
    for candidate in candidates:
        by_haplotype[candidate.haplotype].append(candidate)

    def load_assembly(item):
        haplotype, members = item
        source_info = assembly_sources[haplotype]
        loaded = []
        with IndexedFasta(
            source_info.fasta_path, source_info.fai_path,
        ) as source:
            for candidate in sorted(
                members,
                key=lambda row: (
                    row.source_contig, row.source_start, row.source_end,
                    row.candidate_id,
                ),
            ):
                record = source_records[(
                    candidate.novel_partition, candidate.record_id,
                )]
                local_start = max(0, candidate.local_start - anchor)
                local_end = min(
                    record.record_length, candidate.local_end + anchor,
                )
                sequence = source.fetch(
                    record.source_contig,
                    record.source_start + local_start,
                    record.source_start + local_end,
                )
                if len(sequence) != local_end - local_start:
                    raise ValueError(
                        f"{candidate.candidate_id}: incomplete query extraction"
                    )
                loaded.append((
                    candidate.candidate_id,
                    CachedCandidateSequence(
                        local_start, local_end, sequence,
                    ),
                ))
        return loaded

    output: Dict[str, CachedCandidateSequence] = {}
    items = sorted(by_haplotype.items())
    progress_step = max(1, len(items) // 20)
    total_bases = 0
    for completed, rows in enumerate(_bounded_thread_map(
        load_assembly, items, min(jobs, max(1, len(items))),
    ), 1):
        for candidate_id, sequence in rows:
            output[candidate_id] = sequence
            total_bases += len(sequence.sequence)
        if completed == len(items) or completed % progress_step == 0:
            LOG.info(
                "Loaded query windows from %d/%d assemblies",
                completed, len(items),
            )
    LOG.info(
        "Cached %d query windows (%.2f GiB)",
        len(output), total_bases / (1024 ** 3),
    )
    return output


def load_reference_sequences(
    targets_by_partition: Dict[str, Sequence[BlockTarget]],
    assembly_sources: Dict[str, AssemblyFastaSource],
    jobs: int,
) -> Dict[Tuple[str, str], str]:
    """Load every required true-reference contig into shared RAM once."""
    needed: Dict[str, set] = defaultdict(set)
    for targets in targets_by_partition.values():
        for target in targets:
            needed[target.ref_haplotype].add(target.ref_contig)

    def load_assembly(item):
        haplotype, contigs = item
        source_info = assembly_sources[haplotype]
        rows = []
        with IndexedFasta(
            source_info.fasta_path, source_info.fai_path,
        ) as source:
            for contig in sorted(contigs):
                rows.append(((haplotype, contig), source.sequence(contig)))
        return rows

    output: Dict[Tuple[str, str], str] = {}
    total_bases = 0
    items = sorted(needed.items())
    for rows in _bounded_thread_map(
        load_assembly, items, min(jobs, max(1, len(items))),
    ):
        for key, sequence in rows:
            output[key] = sequence
            total_bases += len(sequence)
    LOG.info(
        "Loaded %d required true-reference contigs into RAM (%.2f GiB)",
        len(output), total_bases / (1024 ** 3),
    )
    return output


def placement_sort_key(row: MappedPlacement) -> Tuple:
    return (
        row.distance,
        -row.alignment_score,
        -row.aligned_bases,
        -row.identity,
        row.ref_haplotype,
        row.ref_contig,
        row.ref_start,
        row.ref_end,
        row.mapped_partition,
        row.target_record,
        row.source,
    )


def _reference_block_targets(
    catalog: CoordinateCatalog,
    reference_haplotypes: FrozenSet[str],
) -> List[BlockTarget]:
    """Return only selected block-reference sequence lying on true Ref."""
    output: List[BlockTarget] = []
    seen = set()
    headers_by_source: Dict[Tuple[str, str], List[HeaderRecord]] = (
        defaultdict(list)
    )
    for header in catalog.headers:
        if header.haplotype in reference_haplotypes:
            headers_by_source[(
                header.haplotype, header.source_contig,
            )].append(header)
    for used in catalog.used_novel_regions:
        if used.haplotype not in reference_haplotypes:
            continue
        for header in headers_by_source.get(
            (used.haplotype, used.source_contig), (),
        ):
            if (
                header.source_end <= used.source_start
                or header.source_start >= used.source_end
            ):
                continue
            source_start = max(header.source_start, used.source_start)
            source_end = min(header.source_end, used.source_end)
            identity = (
                header.record_id, source_start, source_end,
                used.haplotype, used.source_contig,
            )
            if identity in seen:
                continue
            seen.add(identity)
            output.append(BlockTarget(
                "",
                header.record_id,
                source_start - header.source_start,
                source_end - header.source_start,
                used.haplotype,
                used.source_contig,
                source_start,
                source_end,
            ))
    output.sort(key=lambda row: (
        row.ref_haplotype, row.ref_contig, row.ref_start, row.ref_end,
        row.input_record,
    ))
    return [
        dataclasses.replace(row, target_id=f"blockref_{index:06d}")
        for index, row in enumerate(output, 1)
    ]


def _write_block_target_fasta(
    targets: Sequence[BlockTarget],
    reference_sequences: Dict[Tuple[str, str], str],
    path: str,
) -> None:
    with open(path, "wt") as out:
        for target in targets:
            try:
                contig_sequence = reference_sequences[(
                    target.ref_haplotype, target.ref_contig,
                )]
            except KeyError as error:
                raise ValueError(
                    f"reference sequence was not cached for "
                    f"{target.ref_haplotype}:{target.ref_contig}"
                ) from error
            sequence = contig_sequence[target.ref_start:target.ref_end]
            if len(sequence) != target.ref_end - target.ref_start:
                raise ValueError(
                    f"{target.ref_haplotype}:{target.ref_contig}:"
                    f"{target.ref_start}-{target.ref_end}: incomplete "
                    "reference extraction"
                )
            out.write(f">{target.target_id}\n{wrap_fasta(sequence)}\n")


def _write_candidate_queries(
    pieces: Sequence[ResidualPiece],
    candidates: Dict[str, MappingCandidate],
    cached_queries: Mapping[str, CachedCandidateSequence],
    anchor: int,
    path: str,
    reuse_existing: bool = False,
) -> Dict[str, CycleQuery]:
    """Build cycle metadata and write FASTA unless a verified copy exists."""
    metadata: Dict[str, CycleQuery] = {}
    retain_sequences = bool(getattr(
        cached_queries, "retain_cycle_sequences", True,
    ))
    out = None
    fai_out = None
    if not reuse_existing:
        out = open(path, "wb")
        fai_out = open(path + ".fai", "wt")
    try:
        for index, piece in enumerate(pieces):
            try:
                candidate = candidates[piece.candidate_id]
                cached = cached_queries[piece.candidate_id]
            except KeyError as error:
                raise ValueError(
                    f"{piece.candidate_id}: residual lacks cached query data"
                ) from error
            if (
                cached.local_start > candidate.local_start
                or cached.local_end < candidate.local_end
            ):
                raise ValueError(
                    f"{candidate.candidate_id}: cached query window does not "
                    "contain the candidate interval"
                )
            query_id = f"query_{index:08d}"
            if (
                piece.local_start < candidate.local_start
                or piece.local_end <= piece.local_start
                or piece.local_end > candidate.local_end
            ):
                raise ValueError(
                    f"{candidate.candidate_id}: residual local interval "
                    f"{piece.local_start}-{piece.local_end} is outside candidate "
                    f"{candidate.local_start}-{candidate.local_end}"
                )
            extract_start = max(
                cached.local_start, piece.local_start - anchor,
            )
            extract_end = min(
                cached.local_end, piece.local_end + anchor,
            )
            sequence = None
            sequence_length = extract_end - extract_start
            if not reuse_existing:
                sequence = cached.sequence[
                    extract_start - cached.local_start:
                    extract_end - cached.local_start
                ]
                if len(sequence) != sequence_length:
                    raise ValueError(
                        f"{piece.candidate_id}: incomplete cached sequence "
                        f"for {extract_start}-{extract_end}"
                    )
            core_start = piece.local_start - extract_start
            core_end = core_start + piece.local_end - piece.local_start
            metadata[query_id] = CycleQuery(
                query_id,
                piece,
                core_start,
                core_end,
                sequence if retain_sequences else None,
                sequence_length,
            )
            if out is not None and fai_out is not None and sequence is not None:
                out.write(f">{query_id}\n".encode())
                sequence_offset = out.tell()
                out.write(wrap_fasta(sequence).encode())
                out.write(b"\n")
                line_bases = min(DNA_WRAP, sequence_length)
                fai_out.write(
                    f"{query_id}\t{sequence_length}\t{sequence_offset}\t"
                    f"{line_bases}\t{line_bases + 1}\n"
                )
    finally:
        if out is not None:
            out.close()
        if fai_out is not None:
            fai_out.close()
    if reuse_existing:
        existing: List[Tuple[str, int]] = []
        with open(path + ".fai", "rt") as handle:
            for line_number, line in enumerate(handle, 1):
                fields = line.rstrip("\n").split("\t")
                if len(fields) < 2:
                    raise ValueError(
                        f"{path}.fai:{line_number}: malformed FASTA index"
                    )
                existing.append((fields[0], int(fields[1])))
        expected = [
            (query_id, query.sequence_length)
            for query_id, query in metadata.items()
        ]
        if existing != expected:
            raise ValueError(
                f"{path}: reusable cycle FASTA index does not match the "
                "saved residual state"
            )
    return metadata


def _clip_pair_to_query_core(
    pair: Tuple[int, int, int, int],
    strand: str,
    core_start: int,
    core_end: int,
) -> Optional[Tuple[int, int, int, int]]:
    q0, q1, t0, _t1 = pair
    left = max(q0, core_start)
    right = min(q1, core_end)
    if right <= left:
        return None
    if strand == "-":
        target_start = t0 + (q1 - right)
        target_end = t0 + (q1 - left)
    else:
        target_start = t0 + (left - q0)
        target_end = t0 + (right - q0)
    return left, right, target_start, target_end


def _fresh_paired_blocks(
    hits: Sequence[AlignmentHit],
    cycle_query: CycleQuery,
) -> Tuple[List[Interval], List[Tuple[int, int, int, int, AlignmentHit]]]:
    """Clip all blocks to the query core while retaining alignment order."""
    covered_blocks: List[Interval] = []
    accepted: List[Tuple[int, int, int, int, AlignmentHit]] = []
    for hit in hits:
        for pair in hit.aligned_pairs:
            clipped = _clip_pair_to_query_core(
                pair, hit.strand,
                cycle_query.core_start, cycle_query.core_end,
            )
            if clipped is None:
                continue
            query_start, query_end, target_start, target_end = clipped
            covered_blocks.append((query_start, query_end))
            accepted.append((
                query_start,
                query_end,
                target_start,
                target_end,
                hit,
            ))
    return merge_intervals(covered_blocks), accepted


def _collect_reporting_evidence(
    query_meta: Dict[str, CycleQuery],
    hits: Sequence[AlignmentHit],
    target_by_id: Dict[str, BlockTarget],
    evidence_by_candidate: Dict[str, List[PlacementEvidence]],
    mapped_partition: str,
) -> None:
    """Retain all valid targets for reporting without changing cleaning.

    Query identifiers in cycle FASTAs are translated immediately back to the
    candidate's stable local coordinate system.  This is what makes a later
    report-only replay independent of temporary ``query_########`` names.
    """
    for hit in hits:
        cycle_query = query_meta.get(hit.query_id)
        target = target_by_id.get(hit.target_id)
        if cycle_query is None or target is None:
            continue
        clipped_pairs: List[Tuple[int, int, int, int]] = []
        for pair in hit.aligned_pairs:
            clipped = _clip_pair_to_query_core(
                pair,
                hit.strand,
                cycle_query.core_start,
                cycle_query.core_end,
            )
            if clipped is None:
                continue
            q0, q1, t0, t1 = clipped
            query_start = (
                cycle_query.piece.local_start
                + q0 - cycle_query.core_start
            )
            query_end = (
                cycle_query.piece.local_start
                + q1 - cycle_query.core_start
            )
            if target.ref_strand == "-":
                ref_start = target.ref_end - t1
                ref_end = target.ref_end - t0
            else:
                ref_start = target.ref_start + t0
                ref_end = target.ref_start + t1
            clipped_pairs.append((
                query_start, query_end, ref_start, ref_end,
            ))
        if not clipped_pairs:
            continue
        raw_query_start = hit.query_start
        raw_query_end = hit.query_end
        if raw_query_start < 0 or raw_query_end <= raw_query_start:
            raw_query_start = min(row[0] for row in hit.aligned_pairs)
            raw_query_end = max(row[1] for row in hit.aligned_pairs)
        raw_query_start = max(raw_query_start, cycle_query.core_start)
        raw_query_end = min(raw_query_end, cycle_query.core_end)
        if raw_query_end <= raw_query_start:
            continue
        query_start = (
            cycle_query.piece.local_start
            + raw_query_start - cycle_query.core_start
        )
        query_end = (
            cycle_query.piece.local_start
            + raw_query_end - cycle_query.core_start
        )
        strand = hit.strand
        if target.ref_strand == "-":
            strand = "+" if strand == "-" else "-"
        members = evidence_by_candidate.setdefault(
            cycle_query.piece.candidate_id, [],
        )
        members.append(PlacementEvidence(
            cycle_query.piece.candidate_id,
            mapped_partition,
            target.ref_haplotype,
            target.ref_contig,
            strand,
            hit.target_id,
            query_start,
            query_end,
            tuple(clipped_pairs),
            float(hit.alignment_score),
            max(1, int(hit.aligned_bases)),
            float(hit.identity),
            hit.source,
            len(members),
        ))


def _reporting_row(component: ResolvedPlacementComponent) -> Tuple:
    return (
        component.query_id,
        component.component_index,
        component.mapped_partition,
        component.ref_haplotype,
        component.ref_contig,
        component.ref_start,
        component.ref_end,
        component.strand,
        component.target_record,
        component.query_start,
        component.query_end,
        json.dumps(component.query_blocks, separators=(",", ":")),
        json.dumps(component.target_blocks, separators=(",", ":")),
        json.dumps(component.paired_blocks, separators=(",", ":")),
        component.alignment_score,
        component.aligned_bases,
        component.identity,
        component.source,
        component.evidence_ordinal,
        component.liftover_path,
    )


def write_reporting_components(
    path: str,
    components: Iterable[ResolvedPlacementComponent],
) -> None:
    write_tsv(
        path,
        REPORTING_FIELDS,
        (_reporting_row(row) for row in components),
    )


def read_reporting_components(
    paths: Iterable[str],
) -> List[ResolvedPlacementComponent]:
    output: List[ResolvedPlacementComponent] = []
    for path in paths:
        with open(path, "rt", newline="") as handle:
            reader = csv.DictReader(handle, delimiter="\t")
            if reader.fieldnames != list(REPORTING_FIELDS):
                raise ValueError(
                    f"{path}: expected reporting columns "
                    f"{list(REPORTING_FIELDS)}"
                )
            for row in reader:
                output.append(ResolvedPlacementComponent(
                    row["candidate_id"],
                    int(row["component_index"]),
                    row["mapped_partition"],
                    row["ref_haplotype"],
                    row["ref_contig"],
                    int(row["ref_start"]),
                    int(row["ref_end"]),
                    row["strand"],
                    row["target_record"],
                    int(row["query_start"]),
                    int(row["query_end"]),
                    tuple(tuple(map(int, item)) for item in json.loads(
                        row["query_blocks"]
                    )),
                    tuple(tuple(map(int, item)) for item in json.loads(
                        row["target_blocks"]
                    )),
                    tuple(tuple(map(int, item)) for item in json.loads(
                        row["paired_blocks"]
                    )),
                    float(row["alignment_score"]),
                    int(row["aligned_bases"]),
                    float(row["identity"]),
                    row["source"],
                    int(row["evidence_ordinal"]),
                    row.get("liftover_path") or ".",
                ))
    return output


def reporting_components_to_evidence(
    components: Sequence[ResolvedPlacementComponent],
) -> List[PlacementEvidence]:
    return [
        PlacementEvidence(
            row.query_id,
            row.mapped_partition,
            row.ref_haplotype,
            row.ref_contig,
            row.strand,
            row.target_record,
            row.query_start,
            row.query_end,
            row.paired_blocks,
            row.alignment_score,
            row.aligned_bases,
            row.identity,
            row.source,
            index,
        )
        for index, row in enumerate(components)
    ]


def _choose_target(
    hits: Sequence[AlignmentHit],
) -> Optional[Tuple[str, str]]:
    if not hits:
        return None
    best = min(hits, key=lambda hit: (
        -float(hit.alignment_score),
        -int(hit.aligned_bases),
        -float(hit.identity),
        hit.target_id,
        hit.strand,
        hit.target_start,
        hit.target_end,
        hit.source,
    ))
    return best.target_id, best.strand


def _stage_query_thresholds(
    query_meta: Dict[str, CycleQuery],
    min_score: int,
) -> Dict[str, float]:
    """Use the small-query alignment rule from the light cleaner."""
    return {
        query_id: min(float(min_score), 0.5 * query.sequence_length)
        for query_id, query in query_meta.items()
    }


def _subtract_stage_hits(
    current: Sequence[ResidualPiece],
    query_meta: Dict[str, CycleQuery],
    hits: Sequence[AlignmentHit],
    target_by_id: Dict[str, BlockTarget],
    accumulators: Dict[str, PlacementAccumulator],
    min_score: int,
    cached_queries: Optional[
        Mapping[str, CachedCandidateSequence]
    ] = None,
) -> Tuple[List[ResidualPiece], int, float]:
    """Subtract one aligner's accepted coverage and retain score-passing gaps."""
    del current  # query_meta is the authoritative stage-piece ordering.
    hits_by_query: Dict[str, List[AlignmentHit]] = defaultdict(list)
    for hit in hits:
        if hit.target_id in target_by_id:
            hits_by_query[hit.query_id].append(hit)

    next_pieces: List[ResidualPiece] = []
    removed_bp = 0
    remaining_score = 0.0
    for query_id, cycle_query in query_meta.items():
        piece = cycle_query.piece
        accumulator = accumulators.get(piece.candidate_id)
        member_hits = [
            hit
            for hit in hits_by_query.get(query_id, ())
            if any(
                _clip_pair_to_query_core(
                    pair,
                    hit.strand,
                    cycle_query.core_start,
                    cycle_query.core_end,
                ) is not None
                for pair in hit.aligned_pairs
            )
        ]
        selected = (
            (accumulator.target_id, accumulator.strand)
            if accumulator is not None
            else _choose_target(member_hits)
        )
        if selected is None:
            eligible: Sequence[AlignmentHit] = ()
        else:
            target_id, strand = selected
            eligible = [
                hit for hit in member_hits
                if hit.target_id == target_id and hit.strand == strand
            ]
            if accumulator is None and eligible:
                accumulator = PlacementAccumulator(
                    target_id, strand, [], [], [], [], 0.0, 0.0, 0, set(),
                )
                accumulators[piece.candidate_id] = accumulator

        covered_anchor, accepted = _fresh_paired_blocks(
            eligible, cycle_query,
        )
        covered_core = merge_intervals([
            (
                start - cycle_query.core_start,
                end - cycle_query.core_start,
            )
            for start, end in covered_anchor
        ])
        if cycle_query.sequence is not None:
            core_sequence = cycle_query.sequence[
                cycle_query.core_start:cycle_query.core_end
            ]
        elif cached_queries is not None:
            cached = cached_queries[piece.candidate_id]
            core_sequence = cached.sequence[
                piece.local_start - cached.local_start:
                piece.local_end - cached.local_start
            ]
        else:
            raise ValueError(
                f"{piece.candidate_id}: disk-backed cycle lacks query cache"
            )
        residual = find_unaligned_segments(
            covered_core,
            piece.local_end - piece.local_start,
            core_sequence,
            min_score,
            MASKED_BASE_WEIGHT,
        )
        retained_core: List[Interval] = []
        for start, end in residual:
            for nstart, nend, sequence in iter_n_free_segments(
                core_sequence[start:end]
            ):
                if novelty_score(sequence, MASKED_BASE_WEIGHT) >= min_score:
                    retained_core.append((
                        start + nstart, start + nend,
                    ))
        retained_core = merge_intervals(retained_core)
        removed = complement_intervals(
            retained_core, 0, piece.local_end - piece.local_start,
        )
        removed_bp += sum(end - start for start, end in removed)
        if accumulator is not None:
            accumulator.query_masks.extend([
                (
                    piece.local_start + start,
                    piece.local_start + end,
                )
                for start, end in removed
            ])
            target = target_by_id[accumulator.target_id]
            for q0, q1, t0, t1, hit in accepted:
                query_local_start = (
                    piece.local_start + q0 - cycle_query.core_start
                )
                query_local_end = (
                    piece.local_start + q1 - cycle_query.core_start
                )
                accumulator.aligned_query_blocks.append((
                    query_local_start, query_local_end,
                ))
                if target.ref_strand == "-":
                    target_block = (
                        target.ref_end - t1,
                        target.ref_end - t0,
                    )
                else:
                    target_block = (
                        target.ref_start + t0,
                        target.ref_start + t1,
                    )
                accumulator.target_blocks.append(target_block)
                accumulator.paired_blocks.append(ScoredPairedBlock(
                    query_local_start,
                    query_local_end,
                    target_block[0],
                    target_block[1],
                    float(hit.alignment_score),
                    max(1, int(hit.aligned_bases)),
                    float(hit.identity),
                    hit.source,
                    len(accumulator.paired_blocks),
                ))
                length = q1 - q0
                accumulator.aligned_bases += length
                accumulator.identity_bases += length * float(hit.identity)
                denominator = max(1, int(hit.aligned_bases))
                accumulator.alignment_score += (
                    float(hit.alignment_score) * length / denominator
                )
                accumulator.sources.add(hit.source)
        for start, end in retained_core:
            remaining_score += novelty_score(
                core_sequence[start:end], MASKED_BASE_WEIGHT,
            )
            next_pieces.append(ResidualPiece(
                piece.candidate_id,
                piece.local_start + start,
                piece.local_start + end,
            ))
    return next_pieces, removed_bp, remaining_score


def _ensure_reference_blast_db(
    reference: PreparedBlockReference,
) -> str:
    prefix = os.path.join(
        os.path.dirname(reference.target_fasta), "blastdb", "ref",
    )
    if not os.path.isfile(prefix + ".seq"):
        build_blast_db(reference.target_fasta, prefix, LOG)
    return prefix


def _residual_state(
    pieces: Sequence[ResidualPiece],
) -> Tuple[Tuple[str, int, int], ...]:
    """Ordered residual coordinates used for convergence and resume checks.

    Order is significant because cycle FASTA query identifiers are assigned
    from this sequence.  A reordered input must therefore not reuse an older
    PAF/SAM whose query identifiers referred to a different piece.
    """
    return tuple(
        (piece.candidate_id, piece.local_start, piece.local_end)
        for piece in pieces
    )


def _residual_digest(pieces: Sequence[ResidualPiece]) -> str:
    return hashlib.sha256(json.dumps(
        _residual_state(pieces),
        separators=(",", ":"),
    ).encode()).hexdigest()


def _run_reference_clean_stage(
    current: Sequence[ResidualPiece],
    candidates: Dict[str, MappingCandidate],
    cached_queries: Mapping[str, CachedCandidateSequence],
    reference: PreparedBlockReference,
    accumulators: Dict[str, PlacementAccumulator],
    reporting_evidence: Dict[str, List[PlacementEvidence]],
    stage_dir: str,
    options: Dict[str, object],
    aligner: str,
) -> List[ResidualPiece]:
    """Run exactly one aligner stage against one fixed reference block."""
    if not current:
        return []
    input_digest = _residual_digest(current)
    cycle_marker = os.path.join(stage_dir, "_SUCCESS.json")
    reusable_cycle = False
    if bool(options["resume"]) and os.path.isfile(cycle_marker):
        try:
            with open(cycle_marker, "rt") as handle:
                saved_cycle = json.load(handle)
            reusable_cycle = (
                saved_cycle.get("aligner") == aligner
                and saved_cycle.get("input_residual") == input_digest
                and os.path.isfile(os.path.join(stage_dir, "queries.fa"))
                and os.path.isfile(os.path.join(stage_dir, "queries.fa.fai"))
            )
        except (OSError, ValueError, AttributeError):
            reusable_cycle = False
    if (
        not reusable_cycle
        and bool(options.get("reuse_alignments_only", False))
    ):
        raise RuntimeError(
            f"report-only replay cannot reuse alignment cycle {stage_dir}; "
            "the saved query FASTA/marker does not match the frozen "
            "residual state"
        )
    reuse_only = bool(options.get("reuse_alignments_only", False))
    blast_db_prefix = None
    if reuse_only and reusable_cycle:
        expected_alignment = (
            os.path.join(stage_dir, "local_block.blastn.sam")
            if aligner == "blastn"
            else os.path.join(
                stage_dir, f"local_block.{aligner}.paf",
            )
        )
        if not os.path.isfile(expected_alignment):
            raise FileNotFoundError(
                "report-only replay will not create missing alignment "
                f"output: {expected_alignment}"
            )
        if aligner == "blastn":
            blast_db_prefix = os.path.join(
                os.path.dirname(reference.target_fasta), "blastdb", "ref",
            )
            if not os.path.isfile(blast_db_prefix + ".seq"):
                raise FileNotFoundError(
                    "report-only replay will not rebuild missing BLAST "
                    f"database: {blast_db_prefix}.seq"
                )
    if not reusable_cycle and os.path.isdir(stage_dir):
        shutil.rmtree(stage_dir)
    Path(stage_dir).mkdir(parents=True, exist_ok=True)
    query_fasta = os.path.join(stage_dir, "queries.fa")
    query_meta = _write_candidate_queries(
        current,
        candidates,
        cached_queries,
        int(options["input_anchor"]),
        query_fasta,
        reuse_existing=reusable_cycle,
    )
    if reusable_cycle:
        LOG.info("REUSE: cycle query FASTA %s", query_fasta)
    thresholds = _stage_query_thresholds(
        query_meta, int(options["min_score"]),
    )
    if aligner in {"minimap2", "winnowmap"}:
        minimap_params = dict(MINIMAP_PARAMS)
        if options.get("minimap_query_batch"):
            minimap_params["K"] = str(options["minimap_query_batch"])
        if options.get("minimap_retry_sigkill"):
            minimap_params.update({
                "retry_on_sigkill": True,
                "retry_threads": int(options.get(
                    "minimap_retry_threads", 32,
                )),
                "retry_K": str(options.get(
                    "minimap_retry_query_batch", "100M",
                )),
            })
        hits = collect_alignment_hits(
            query_fasta,
            reference.target_fasta,
            stage_dir,
            int(options["cores"]),
            float(options["min_identity"]),
            int(options["min_score"]),
            LOG,
            winnow_params=LIGHT_WINNOW_PARAMS,
            prefix="local_block",
            skip_blastn=True,
            broad_aligner=aligner,
            minimap_params=minimap_params,
            masked_weight=MASKED_BASE_WEIGHT,
            min_score_by_query=thresholds,
            parse_min_segment=min(50, int(options["min_score"])),
            reuse_existing=reusable_cycle,
            query_workers=int(options["cores"]),
        )
    elif aligner == "blastn":
        hits = collect_blastn_hits(
            query_fasta,
            reference.target_fasta,
            stage_dir,
            int(options["cores"]),
            float(options["min_identity"]),
            int(options["min_score"]),
            LOG,
            blast_word_size=int(options["blast_word_size"]),
            blast_evalue=str(options["blast_evalue"]),
            blast_max_target_seqs=int(options["blast_max_target_seqs"]),
            prefix="local_block",
            blast_db_prefix=(
                blast_db_prefix or _ensure_reference_blast_db(reference)
            ),
            masked_weight=MASKED_BASE_WEIGHT,
            min_score_by_query=thresholds,
            parse_min_segment=min(50, int(options["min_score"])),
            reuse_existing=reusable_cycle,
            query_workers=int(options["cores"]),
        )
    else:
        raise ValueError(f"unknown local cleaning aligner {aligner!r}")

    target_by_id = {
        target.target_id: target for target in reference.targets
    }
    _collect_reporting_evidence(
        query_meta,
        hits,
        target_by_id,
        reporting_evidence,
        reference.mapped_partition,
    )
    next_pieces, removed_bp, remaining_score = _subtract_stage_hits(
        current,
        query_meta,
        hits,
        target_by_id,
        accumulators,
        int(options["min_score"]),
        cached_queries,
    )
    LOG.info(
        "%s [%s]: %d -> %d residual pieces; removed %d bp; "
        "remaining novelty score %.1f",
        reference.mapped_partition,
        aligner,
        len(current),
        len(next_pieces),
        removed_bp,
        remaining_score,
    )
    temporary_marker = cycle_marker + ".tmp"
    with open(temporary_marker, "wt") as out:
        json.dump({
            "aligner": aligner,
            "input_residual": input_digest,
            "output_residual": _residual_digest(next_pieces),
            "input_pieces": len(current),
            "output_pieces": len(next_pieces),
            "removed_bp": removed_bp,
        }, out, sort_keys=True)
        out.write("\n")
    os.replace(temporary_marker, cycle_marker)
    return next_pieces


def resolve_scored_paired_blocks(
    blocks: Sequence[ScoredPairedBlock],
    strand: str,
) -> List[ScoredPairedBlock]:
    """Sweep query endpoints and retain the highest-priority active hit."""
    if strand not in {"+", "-"}:
        raise ValueError(f"invalid placement strand {strand!r}")
    starts: Dict[int, List[int]] = defaultdict(list)
    ends: Dict[int, List[int]] = defaultdict(list)
    coordinates = set()
    for index, block in enumerate(blocks):
        if (
            block.query_start < 0
            or block.query_end <= block.query_start
            or block.ref_start < 0
            or block.ref_end <= block.ref_start
            or block.query_end - block.query_start
            != block.ref_end - block.ref_start
        ):
            raise ValueError(
                "invalid scored paired block "
                f"{block.query_start}-{block.query_end}:"
                f"{block.ref_start}-{block.ref_end}"
            )
        starts[block.query_start].append(index)
        ends[block.query_end].append(index)
        coordinates.add(block.query_start)
        coordinates.add(block.query_end)
    ordered_coordinates = sorted(coordinates)
    active = set()
    priority_heap: List[Tuple[float, int, float, int, int]] = []
    selected: List[List[int]] = []
    for position_index, coordinate in enumerate(ordered_coordinates[:-1]):
        for index in ends.get(coordinate, ()):
            active.discard(index)
        for index in starts.get(coordinate, ()):
            active.add(index)
            block = blocks[index]
            heapq.heappush(priority_heap, (
                -block.alignment_score,
                -block.alignment_bases,
                -block.identity,
                block.ordinal,
                index,
            ))
        while priority_heap and priority_heap[0][4] not in active:
            heapq.heappop(priority_heap)
        next_coordinate = ordered_coordinates[position_index + 1]
        if not priority_heap or next_coordinate <= coordinate:
            continue
        winner = priority_heap[0][4]
        if (
            selected
            and selected[-1][2] == winner
            and selected[-1][1] == coordinate
        ):
            selected[-1][1] = next_coordinate
        else:
            selected.append([coordinate, next_coordinate, winner])

    output: List[ScoredPairedBlock] = []
    for query_start, query_end, winner in selected:
        block = blocks[winner]
        if strand == "+":
            ref_start = (
                block.ref_start + query_start - block.query_start
            )
            ref_end = (
                block.ref_start + query_end - block.query_start
            )
        else:
            ref_start = block.ref_start + block.query_end - query_end
            ref_end = block.ref_start + block.query_end - query_start
        output.append(ScoredPairedBlock(
            query_start,
            query_end,
            ref_start,
            ref_end,
            block.alignment_score,
            block.alignment_bases,
            block.identity,
            block.source,
            block.ordinal,
        ))
    return output


_RESOLVE_ACCUMULATORS: Dict[str, PlacementAccumulator] = {}
_RESOLVE_TARGETS: Dict[str, BlockTarget] = {}
_RESOLVE_PARTITION = ""


def _placement_from_accumulator(
    candidate_id: str,
) -> Optional[MappedPlacement]:
    return _build_resolved_placement(
        candidate_id,
        _RESOLVE_PARTITION,
        _RESOLVE_ACCUMULATORS[candidate_id],
        _RESOLVE_TARGETS,
    )


def _build_resolved_placement(
    candidate_id: str,
    mapped_partition: str,
    accumulator: PlacementAccumulator,
    target_by_id: Dict[str, BlockTarget],
) -> Optional[MappedPlacement]:
    query_masks = merge_intervals(accumulator.query_masks)
    target = target_by_id[accumulator.target_id]
    strand = accumulator.strand
    if target.ref_strand == "-":
        strand = "+" if strand == "-" else "-"
    resolved_blocks = resolve_scored_paired_blocks(
        accumulator.paired_blocks, strand,
    )
    paired_blocks = [
        (
            block.query_start,
            block.query_end,
            block.ref_start,
            block.ref_end,
        )
        for block in resolved_blocks
    ]
    aligned_query = merge_intervals([
        (block.query_start, block.query_end)
        for block in resolved_blocks
    ])
    target_blocks = merge_intervals([
        (block.ref_start, block.ref_end)
        for block in resolved_blocks
    ])
    if not aligned_query or not target_blocks or not query_masks:
        return None
    aligned_bases = sum(
        block.query_end - block.query_start
        for block in resolved_blocks
    )
    identity = sum(
        (block.query_end - block.query_start) * block.identity
        for block in resolved_blocks
    ) / aligned_bases
    alignment_score = sum(
        block.alignment_score
        * (block.query_end - block.query_start)
        / block.alignment_bases
        for block in resolved_blocks
    )
    return MappedPlacement(
        candidate_id,
        mapped_partition,
        target.ref_haplotype,
        target.ref_contig,
        min(start for start, _end in target_blocks),
        max(end for _start, end in target_blocks),
        strand,
        0,
        alignment_score,
        identity,
        aligned_bases,
        ",".join(sorted({
            block.source for block in resolved_blocks
        })),
        accumulator.target_id,
        json.dumps(query_masks, separators=(",", ":")),
        json.dumps(target_blocks, separators=(",", ":")),
        json.dumps(paired_blocks, separators=(",", ":")),
        ".",
    )


def _resolve_placement_chunk(
    candidate_ids: Sequence[str],
) -> List[MappedPlacement]:
    return [
        placement
        for candidate_id in candidate_ids
        for placement in (_placement_from_accumulator(candidate_id),)
        if placement is not None
    ]


def _placements_from_accumulators(
    mapped_partition: str,
    accumulators: Dict[str, PlacementAccumulator],
    target_by_id: Dict[str, BlockTarget],
    workers: int = 1,
) -> List[MappedPlacement]:
    candidate_ids = sorted(accumulators)
    worker_count = min(
        max(1, int(workers)), max(1, len(candidate_ids)),
    )
    context = mp_context()
    use_processes = (
        worker_count > 1
        and context.get_start_method() == "fork"
        and len(candidate_ids) >= worker_count * 4
        and threading.current_thread() is threading.main_thread()
    )
    if use_processes:
        global _RESOLVE_ACCUMULATORS
        global _RESOLVE_TARGETS
        global _RESOLVE_PARTITION
        _RESOLVE_ACCUMULATORS = accumulators
        _RESOLVE_TARGETS = target_by_id
        _RESOLVE_PARTITION = mapped_partition
        try:
            buckets: List[Tuple[int, int, List[str]]] = [
                (0, index, []) for index in range(worker_count)
            ]
            heapq.heapify(buckets)
            for candidate_id in sorted(
                candidate_ids,
                key=lambda value: -len(
                    accumulators[value].paired_blocks
                ),
            ):
                load, bucket_index, members = heapq.heappop(buckets)
                members.append(candidate_id)
                heapq.heappush(buckets, (
                    load + len(accumulators[candidate_id].paired_blocks),
                    bucket_index,
                    members,
                ))
            chunks = [
                members
                for _load, _index, members in sorted(
                    buckets, key=lambda row: row[1],
                )
                if members
            ]
            LOG.info(
                "%s: resolving alignment-priority sweeps for %d queries "
                "with %d workers",
                mapped_partition, len(candidate_ids), len(chunks),
            )
            with context.Pool(processes=len(chunks)) as pool:
                output = [
                    placement
                    for rows in pool.map(_resolve_placement_chunk, chunks)
                    for placement in rows
                ]
        finally:
            _RESOLVE_ACCUMULATORS = {}
            _RESOLVE_TARGETS = {}
            _RESOLVE_PARTITION = ""
    else:
        output = [
            placement
            for candidate_id in candidate_ids
            for placement in (_build_resolved_placement(
                candidate_id,
                mapped_partition,
                accumulators[candidate_id],
                target_by_id,
            ),)
            if placement is not None
        ]
    output.sort(key=lambda row: row.candidate_id)
    return output


def read_partition_reference_targets(
    mapped_partition: str,
    partition_files: PartitionFiles,
    reference_haplotypes: FrozenSet[str],
) -> Tuple[str, Tuple[BlockTarget, ...]]:
    """Resolve one partition's exact Ref target from its ``.header`` and BED."""
    headers = parse_header_records(
        partition_files,
        selected_haplotypes=set(reference_haplotypes),
    )
    catalog = CoordinateCatalog(
        partition_files,
        headers,
        [],
        parse_used_novel_regions(
            partition_files,
            {row.haplotype for row in headers},
        ),
    )
    targets = _reference_block_targets(
        catalog, reference_haplotypes,
    )
    return mapped_partition, tuple(targets)


def materialize_partition_reference(
    mapped_partition: str,
    targets: Sequence[BlockTarget],
    reference_sequences: Dict[Tuple[str, str], str],
    work_root: str,
    options: Dict[str, object],
) -> PreparedBlockReference:
    """Write one small reusable target FASTA from the shared RAM cache."""
    workdir = os.path.join(work_root, "references", mapped_partition)
    target_fasta = os.path.join(workdir, "reference.fa")
    marker = os.path.join(workdir, "_SUCCESS.json")
    task_signature = hashlib.sha256(json.dumps(
        {
            "mapping_signature": options["mapping_signature"],
            "mapped_partition": mapped_partition,
            "targets": [
                dataclasses.asdict(target) for target in targets
            ],
        },
        sort_keys=True,
        separators=(",", ":"),
    ).encode()).hexdigest()
    reusable = False
    if bool(options["resume"]) and os.path.isfile(marker):
        try:
            with open(marker, "rt") as handle:
                saved = json.load(handle)
            reusable = (
                saved.get("signature") == task_signature
                and (
                    not targets
                    or (
                        os.path.isfile(target_fasta)
                        and os.path.getsize(target_fasta) > 0
                    )
                )
            )
        except (OSError, ValueError, AttributeError):
            reusable = False
    if not reusable:
        if os.path.isdir(workdir):
            shutil.rmtree(workdir)
        Path(workdir).mkdir(parents=True, exist_ok=True)
        if targets:
            _write_block_target_fasta(
                targets, reference_sequences, target_fasta,
            )
        with open(marker, "wt") as out:
            json.dump(
                {
                    "signature": task_signature,
                    "targets": len(targets),
                },
                out,
                sort_keys=True,
            )
            out.write("\n")
    if not targets:
        LOG.info(
            "%s: selected block is not on a true-reference haplotype; "
            "local mapping will be skipped",
            mapped_partition,
        )
    return PreparedBlockReference(
        mapped_partition,
        target_fasta,
        tuple(targets),
        task_signature,
    )


def map_partition_candidates(
    candidates: Sequence[MappingCandidate],
    cached_queries: Mapping[str, CachedCandidateSequence],
    reference: PreparedBlockReference,
    work_root: str,
    options: Dict[str, object],
) -> str:
    """Map every local-novel query assigned to one partition in one job."""
    workdir = os.path.join(
        work_root, "partitions", reference.mapped_partition,
    )
    result_path = os.path.join(workdir, "placements.tsv")
    reporting_path = os.path.join(workdir, "reporting_placements.tsv")
    marker = os.path.join(workdir, "_SUCCESS.json")
    candidate_by_id = {
        candidate.candidate_id: candidate for candidate in candidates
    }
    task_signature = hashlib.sha256(json.dumps(
        {
            "mapping_signature": options["mapping_signature"],
            "reference": (
                reference.mapped_partition, reference.signature,
            ),
            "candidate_ids": sorted(candidate_by_id),
            "converge_aligners": bool(
                options.get("converge_aligners", False)
            ),
        },
        sort_keys=True,
        separators=(",", ":"),
    ).encode()).hexdigest()
    if bool(options["resume"]) and os.path.isfile(marker) and os.path.isfile(
        result_path
    ):
        try:
            with open(marker, "rt") as handle:
                saved = json.load(handle)
                if (
                    saved.get("signature") == task_signature
                    and saved.get("placement_resolution")
                    == PLACEMENT_RESOLUTION_VERSION
                    and saved.get("reporting_resolution")
                    == REPORTING_RESOLUTION_VERSION
                    and os.path.isfile(reporting_path)
                ):
                    return result_path
        except (OSError, ValueError, AttributeError):
            pass
    task_marker = os.path.join(workdir, "_TASK.json")
    reusable_task = False
    if bool(options["resume"]) and os.path.isfile(task_marker):
        try:
            with open(task_marker, "rt") as handle:
                reusable_task = (
                    json.load(handle).get("signature") == task_signature
                )
        except (OSError, ValueError, AttributeError):
            reusable_task = False
    if os.path.isdir(workdir) and not reusable_task:
        shutil.rmtree(workdir)
    Path(workdir).mkdir(parents=True, exist_ok=True)
    temporary_task_marker = task_marker + ".tmp"
    with open(temporary_task_marker, "wt") as out:
        json.dump({"signature": task_signature}, out, sort_keys=True)
        out.write("\n")
    os.replace(temporary_task_marker, task_marker)

    current = [
        ResidualPiece(
            candidate.candidate_id,
            candidate.local_start,
            candidate.local_end,
        )
        for candidate in candidates
    ]
    accumulators: Dict[str, PlacementAccumulator] = {}
    reporting_evidence: Dict[str, List[PlacementEvidence]] = {}
    converge_aligners = bool(options.get("converge_aligners", False))
    for aligner in ("minimap2", "winnowmap", "blastn"):
        if aligner == "blastn" and bool(options["skip_blastn"]):
            break
        cycle = 0
        while current:
            before = _residual_state(current)
            stage_dir = (
                os.path.join(workdir, aligner, f"cycle{cycle:03d}")
                if converge_aligners
                else os.path.join(workdir, aligner)
            )
            next_pieces = _run_reference_clean_stage(
                current,
                candidate_by_id,
                cached_queries,
                reference,
                accumulators,
                reporting_evidence,
                stage_dir,
                options,
                aligner,
            )
            current = next_pieces
            if not converge_aligners:
                break
            if _residual_state(current) == before:
                LOG.info(
                    "%s [%s] converged after cycle %d: no additional "
                    "query coverage",
                    reference.mapped_partition,
                    aligner,
                    cycle,
                )
                break
            cycle += 1
        if not current:
            break
    target_by_id = {
        target.target_id: target for target in reference.targets
    }
    placements = _placements_from_accumulators(
        reference.mapped_partition,
        accumulators,
        target_by_id,
        int(options["cores"]),
    )
    resolved_reporting = resolve_all_evidence(
        row
        for rows in reporting_evidence.values()
        for row in rows
    )
    reporting_components = [
        row
        for candidate_id in sorted(resolved_reporting)
        for row in resolved_reporting[candidate_id]
    ]
    LOG.info(
        "%s: mapped %d candidates; %s after staged cleaning",
        reference.mapped_partition,
        len(candidates),
        (
            "no score-passing residual"
            if not current
            else f"{len(current)} residual piece(s)"
        ),
    )

    write_tsv(
        result_path,
        [field.name for field in dataclasses.fields(MappedPlacement)],
        (
            dataclasses.astuple(row)
            for row in sorted(placements, key=placement_sort_key)
        ),
    )
    write_reporting_components(reporting_path, reporting_components)
    with open(marker, "wt") as out:
        json.dump({
            "signature": task_signature,
            "placement_resolution": PLACEMENT_RESOLUTION_VERSION,
            "reporting_resolution": REPORTING_RESOLUTION_VERSION,
            "reporting_components": len(reporting_components),
        }, out, sort_keys=True)
        out.write("\n")
    return result_path


def read_placements(paths: Iterable[str]) -> List[MappedPlacement]:
    output: List[MappedPlacement] = []
    for path in paths:
        with open(path, "rt", newline="") as handle:
            reader = csv.DictReader(handle, delimiter="\t")
            if (
                reader.fieldnames is None
                or "query_blocks" not in reader.fieldnames
            ):
                raise ValueError(
                    f"{path}: missing required query_blocks column; rerun "
                    "the updated map_partition_local_novel.py"
                )
            for row in reader:
                output.append(MappedPlacement(
                    row["candidate_id"],
                    row["mapped_partition"],
                    row["ref_haplotype"],
                    row["ref_contig"],
                    int(row["ref_start"]),
                    int(row["ref_end"]),
                    row["strand"],
                    int(row["distance"]),
                    float(row["alignment_score"]),
                    float(row["identity"]),
                    int(row["aligned_bases"]),
                    row["source"],
                    row["target_record"],
                    row["query_blocks"],
                    row["target_blocks"],
                    row.get("paired_blocks", "[]"),
                    row.get("liftover_path") or ".",
                ))
    return output


def write_local_alignment_summary(
    candidates: Sequence[MappingCandidate],
    placements: Sequence[MappedPlacement],
    output_dir: str,
) -> None:
    best: Dict[str, MappedPlacement] = {}
    for row in placements:
        prior = best.get(row.candidate_id)
        if prior is None or placement_sort_key(row) < placement_sort_key(prior):
            best[row.candidate_id] = row
    fields = (
        "query_id", "query_kind", "partition", "status",
        "ref_haplotype", "ref_contig", "ref_start", "ref_end",
        "strand", "alignment_score", "aligned_bases", "identity",
        "source", "target_record", "query_blocks", "target_blocks",
        "paired_blocks", "graphic_cigar",
    )
    write_tsv(
        os.path.join(output_dir, "local_novel_alignment_summary.tsv"),
        fields,
        (
            (
                candidate.candidate_id, "local_novel",
                candidate.novel_partition,
                "mapped" if candidate.candidate_id in best else "unmapped",
                *(
                    (
                        placement.ref_haplotype, placement.ref_contig,
                        placement.ref_start, placement.ref_end,
                        placement.strand, placement.alignment_score,
                        placement.aligned_bases, placement.identity,
                        placement.source, placement.target_record,
                        placement.query_blocks, placement.target_blocks,
                        placement.paired_blocks, placement.liftover_path,
                    )
                    if (placement := best.get(candidate.candidate_id)) is not None
                    else (".", ".", ".", ".", ".", 0, 0, 0, ".", ".",
                          "[]", "[]", "[]", ".")
                ),
            )
            for candidate in candidates
        ),
    )


def export_mapping_results(
    candidates: Sequence[MappingCandidate],
    placements: Sequence[MappedPlacement],
    output_dir: str,
    merge_gap: int,
) -> None:
    placement_fields = [
        field.name for field in dataclasses.fields(MappedPlacement)
    ]
    write_tsv(
        os.path.join(output_dir, "all_local_novel_alignments.tsv"),
        placement_fields,
        (
            dataclasses.astuple(row)
            for row in sorted(
                placements,
                key=lambda row: (row.candidate_id, placement_sort_key(row)),
            )
        ),
    )
    candidate_by_id = {row.candidate_id: row for row in candidates}
    best: Dict[str, MappedPlacement] = {}
    for row in placements:
        prior = best.get(row.candidate_id)
        if prior is None or placement_sort_key(row) < placement_sort_key(prior):
            best[row.candidate_id] = row

    unresolved = [
        row for row in candidates if row.candidate_id not in best
    ]
    write_bed_table(
        os.path.join(output_dir, "unresolved_local_novel.bed"),
        (
            "source_contig", "source_start", "source_end", "candidate_id",
            "score", "strand", "haplotype", "novel_partition",
            "lifted_partitions", "interval_id", "record_id",
        ),
        (
            (
                row.source_contig, row.source_start, row.source_end,
                row.candidate_id, 0, "+", row.haplotype,
                row.novel_partition, ",".join(sorted(row.lifted_partitions)),
                row.interval_id, row.record_id,
            )
            for row in unresolved
        ),
    )

    selected = sorted(
        (
            (
            placement, candidate_by_id[candidate_id]
            )
            for candidate_id, placement in best.items()
        ),
        key=lambda item: (
            item[0].ref_haplotype,
            item[0].ref_contig,
            item[0].ref_start,
            item[0].ref_end,
            item[0].candidate_id,
        ),
    )
    merged: List[Dict[str, object]] = []
    for placement, candidate in selected:
        if (
            merged
            and merged[-1]["ref_haplotype"] == placement.ref_haplotype
            and merged[-1]["ref_contig"] == placement.ref_contig
            and placement.ref_start - int(merged[-1]["ref_end"]) <= merge_gap
        ):
            item = merged[-1]
            item["ref_end"] = max(int(item["ref_end"]), placement.ref_end)
            item["alignment_score"] = max(
                float(item["alignment_score"]), placement.alignment_score,
            )
            item["distance"] = min(int(item["distance"]), placement.distance)
            item["candidate_ids"].append(candidate.candidate_id)
            item["interval_ids"].add(candidate.interval_id)
            item["source_loci"].append(
                f"{candidate.haplotype}:{candidate.source_contig}:"
                f"{candidate.source_start}-{candidate.source_end}"
            )
            item["novel_partitions"].add(candidate.novel_partition)
            item["mapped_partitions"].add(placement.mapped_partition)
            item["strands"].add(placement.strand)
            item["sources"].add(placement.source)
            continue
        merged.append({
            "ref_haplotype": placement.ref_haplotype,
            "ref_contig": placement.ref_contig,
            "ref_start": placement.ref_start,
            "ref_end": placement.ref_end,
            "alignment_score": placement.alignment_score,
            "distance": placement.distance,
            "candidate_ids": [candidate.candidate_id],
            "interval_ids": {candidate.interval_id},
            "source_loci": [
                f"{candidate.haplotype}:{candidate.source_contig}:"
                f"{candidate.source_start}-{candidate.source_end}"
            ],
            "novel_partitions": {candidate.novel_partition},
            "mapped_partitions": {placement.mapped_partition},
            "strands": {placement.strand},
            "sources": {placement.source},
        })
    write_bed_table(
        os.path.join(output_dir, "local_novel.bed"),
        (
            "ref_contig", "ref_start", "ref_end", "local_novel_id",
            "alignment_score", "strand", "ref_haplotype", "source_loci",
            "novel_partitions", "mapped_partitions",
            "distance_to_reference_unique", "candidate_ids", "interval_ids",
            "alignment_sources",
        ),
        (
            (
                item["ref_contig"],
                item["ref_start"],
                item["ref_end"],
                f"local_novel_{index:09d}",
                f"{float(item['alignment_score']):.1f}",
                ",".join(sorted(item["strands"])),
                item["ref_haplotype"],
                ";".join(item["source_loci"]),
                ",".join(sorted(item["novel_partitions"])),
                ",".join(sorted(item["mapped_partitions"])),
                item["distance"],
                ",".join(item["candidate_ids"]),
                ",".join(sorted(item["interval_ids"])),
                ",".join(sorted(item["sources"])),
            )
            for index, item in enumerate(merged, 1)
        ),
    )
    write_local_alignment_summary(candidates, placements, output_dir)


def copy_global_novel(analysis_dir: str, output_dir: str) -> None:
    source = os.path.join(analysis_dir, "global_novel.bed")
    target = os.path.join(output_dir, "global_novel.bed")
    if os.path.abspath(source) == os.path.abspath(target):
        return
    temporary = target + ".tmp"
    try:
        shutil.copyfile(source, temporary)
        os.replace(temporary, target)
    finally:
        try:
            os.remove(temporary)
        except FileNotFoundError:
            pass


def mapping_signature(
    analysis_dir: str,
    candidates_path: str,
    inputs_path: str,
    reference_haplotypes: FrozenSet[str],
    assembly_sources: Dict[str, AssemblyFastaSource],
    args: argparse.Namespace,
) -> str:
    digest = hashlib.sha256()
    digest.update(b"local-novel-mapping-v10-paired-blocks\n")
    for path in (
        os.path.join(analysis_dir, "_COORDINATE_ANALYSIS_SUCCESS.json"),
        candidates_path,
        inputs_path,
        os.path.abspath(args.query_paths),
        *(
            path
            for source in sorted(
                assembly_sources.values(),
                key=lambda row: row.haplotype,
            )
            for path in (source.fasta_path, source.fai_path)
        ),
    ):
        file_stat = os.stat(path)
        digest.update(
            f"{path}\t{file_stat.st_size}\t{file_stat.st_mtime_ns}\n".encode()
        )
    settings = {
        "reference_haplotypes": sorted(reference_haplotypes),
        "input_anchor": args.input_anchor,
        "min_score": args.min_score,
        "min_identity": args.min_identity,
        "blast_word_size": args.blast_word_size,
        "blast_evalue": args.blast_evalue,
        "blast_max_target_seqs": args.blast_max_target_seqs,
        "skip_blastn": args.skip_blastn,
        "merge_gap": args.merge_gap,
    }
    digest.update(json.dumps(settings, sort_keys=True).encode())
    return digest.hexdigest()


def parse_args(argv: Optional[Sequence[str]] = None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=(
            "Clean independent local-novel candidates with staged Minimap2, "
            "Winnowmap, and BLASTN passes against exact selected block "
            "references on true-reference haplotypes."
        )
    )
    parser.add_argument(
        "-a", "--analysis-dir", required=True,
        help="coordinate-analysis output directory",
    )
    parser.add_argument("-o", "--output-dir", required=True)
    parser.add_argument(
        "-q", "--query-paths", required=True,
        help=(
            "original HAPLOTYPE FASTA (with adjacent FASTA.fai) list used for coordinate-based "
            "query and reference extraction"
        ),
    )
    parser.add_argument(
        "-j", "--jobs", type=int, default=16,
        help="concurrent local-novel candidate jobs (default: 16)",
    )
    parser.add_argument(
        "-c", "--cores", type=int, default=4,
        help="aligner threads per job (total potential threads: -j * -c)",
    )
    parser.add_argument(
        "--reference-haplotypes",
        help=(
            "must match the coordinate-analysis reference haplotypes; by "
            "default they are read from its success marker"
        ),
    )
    parser.add_argument(
        "--input-anchor", type=int,
        help=(
            "query anchor size; defaults to and must match the coordinate "
            "analysis setting"
        ),
    )
    parser.add_argument("--merge-gap", type=int, default=500)
    parser.add_argument("--min-score", type=int, default=100)
    parser.add_argument("--min-identity", type=float, default=95.0)
    parser.add_argument("--blast-word-size", type=int, default=50)
    parser.add_argument("--blast-evalue", default="1e-300")
    parser.add_argument("--blast-max-target-seqs", type=int, default=100)
    parser.add_argument(
        "--max-cycles", type=int, default=0,
        help=(
            "deprecated compatibility option; staged local mapping now runs "
            "each aligner exactly once"
        ),
    )
    parser.add_argument("--skip-blastn", action="store_true")
    parser.add_argument("--resume", action="store_true")
    parser.add_argument(
        "--reuse-alignments-only",
        action="store_true",
        help=(
            "rebuild reporting tables only from matching saved query "
            "FASTAs and PAF/SAM files; fail instead of running an aligner"
        ),
    )
    parser.add_argument(
        "--log-level", choices=("DEBUG", "INFO", "WARNING", "ERROR"),
        default="INFO",
    )
    args = parser.parse_args(argv)
    for name in ("jobs", "cores", "min_score"):
        if getattr(args, name) < 1:
            parser.error(f"--{name.replace('_', '-')} must be positive")
    if (
        (args.input_anchor is not None and args.input_anchor < 0)
        or args.merge_gap < 0
        or args.max_cycles < 0
    ):
        parser.error(
            "--input-anchor, --merge-gap, and --max-cycles cannot be negative"
        )
    if not 0 < args.min_identity <= 100:
        parser.error("--min-identity must be in (0,100]")
    return args


def run(args: argparse.Namespace) -> None:
    analysis_dir = os.path.abspath(args.analysis_dir)
    output_dir = os.path.abspath(args.output_dir)
    Path(output_dir).mkdir(parents=True, exist_ok=True)
    analysis_marker = os.path.join(
        analysis_dir, "_COORDINATE_ANALYSIS_SUCCESS.json",
    )
    with open(analysis_marker, "rt") as handle:
        analysis_metadata = json.load(handle)
    saved_anchor = int(analysis_metadata.get("input_anchor", 50))
    if args.input_anchor is None:
        args.input_anchor = saved_anchor
    elif args.input_anchor != saved_anchor:
        raise ValueError(
            f"--input-anchor {args.input_anchor} does not match coordinate "
            f"analysis value {saved_anchor}"
        )
    saved_reference = frozenset(analysis_metadata["reference_haplotypes"])
    if args.reference_haplotypes:
        reference_haplotypes = parse_reference_haplotypes(
            args.reference_haplotypes
        )
        if reference_haplotypes != saved_reference:
            raise ValueError(
                "--reference-haplotypes must match coordinate analysis: "
                + ",".join(sorted(saved_reference))
            )
    else:
        reference_haplotypes = saved_reference

    candidates_path = os.path.join(
        analysis_dir, "local_novel_candidates.tsv",
    )
    inputs_path = os.path.join(analysis_dir, "partition_inputs.tsv")
    assembly_sources = read_assembly_sources(
        os.path.abspath(args.query_paths)
    )
    missing_reference_sources = reference_haplotypes - set(assembly_sources)
    if missing_reference_sources:
        raise ValueError(
            "--query-paths lacks reference haplotypes: "
            + ",".join(sorted(missing_reference_sources))
        )
    signature = mapping_signature(
        analysis_dir,
        candidates_path,
        inputs_path,
        reference_haplotypes,
        assembly_sources,
        args,
    )
    marker = os.path.join(output_dir, "_LOCAL_NOVEL_MAPPING_SUCCESS.json")
    if args.resume and os.path.isfile(marker):
        try:
            with open(marker, "rt") as handle:
                saved = json.load(handle)
            if (
                saved.get("signature") == signature
                and saved.get("placement_resolution")
                == PLACEMENT_RESOLUTION_VERSION
                and os.path.isfile(os.path.join(output_dir, "local_novel.bed"))
                and os.path.isfile(os.path.join(output_dir, "global_novel.bed"))
                and os.path.isfile(os.path.join(
                    output_dir, "all_local_novel_alignments.tsv",
                ))
                and os.path.isfile(os.path.join(
                    output_dir, "local_novel_alignment_summary.tsv",
                ))
                and saved.get("reporting_resolution")
                == REPORTING_RESOLUTION_VERSION
                and os.path.isfile(os.path.join(
                    output_dir,
                    "all_local_novel_reporting_alignments.tsv",
                ))
            ):
                LOG.info("RESUME: local-novel mapping is complete in %s", output_dir)
                return
        except (OSError, ValueError, AttributeError):
            pass

    partition_inputs = read_partition_inputs(inputs_path)
    candidates = read_candidates(candidates_path)
    needed_partition_names = set()
    for candidate in candidates:
        novel_input = partition_inputs.get(candidate.novel_partition)
        if novel_input is None:
            raise ValueError(
                f"{candidate.candidate_id}: unknown novel partition "
                f"{candidate.novel_partition!r}"
            )
        if os.path.abspath(candidate.source_fasta) != os.path.abspath(
            novel_input.fasta
        ):
            raise ValueError(
                f"{candidate.candidate_id}: source FASTA disagrees with "
                f"partition {candidate.novel_partition}"
            )
        if candidate.haplotype not in assembly_sources:
            raise ValueError(
                f"{candidate.candidate_id}: haplotype "
                f"{candidate.haplotype!r} is absent from --query-paths"
            )
        for partition in candidate.lifted_partitions:
            if partition not in partition_inputs:
                raise ValueError(
                    f"{candidate.candidate_id}: unknown lifted partition "
                    f"{partition!r}"
                )
            needed_partition_names.add(partition)
    LOG.info(
        "Mapping %d local-novel candidates through %d lifted partitions "
        "with %d jobs x %d cores",
        len(candidates), len(needed_partition_names), args.jobs, args.cores,
    )
    candidate_header_records = load_candidate_header_records(
        candidates, partition_inputs, args.jobs,
    )
    if needed_partition_names and not getattr(
        args, "reuse_alignments_only", False,
    ):
        for executable in ("minimap2", "winnowmap"):
            ensure_executable(executable)
        if not args.skip_blastn:
            for executable in ("makeblastdb", "blastn"):
                ensure_executable(executable)

    options: Dict[str, object] = {
        "resume": args.resume,
        "analysis_signature": analysis_metadata["signature"],
        "mapping_signature": signature,
        "cores": args.cores,
        "input_anchor": args.input_anchor,
        "min_identity": args.min_identity,
        "min_score": args.min_score,
        "blast_word_size": args.blast_word_size,
        "blast_evalue": args.blast_evalue,
        "blast_max_target_seqs": args.blast_max_target_seqs,
        "skip_blastn": args.skip_blastn,
        "reuse_alignments_only": getattr(
            args, "reuse_alignments_only", False,
        ),
    }
    work_root = os.path.join(output_dir, "alignment_work")
    needed_partitions = sorted(needed_partition_names)
    targets_by_partition: Dict[str, Tuple[BlockTarget, ...]] = {}
    references: Dict[str, PreparedBlockReference] = {}
    if needed_partitions:
        LOG.info(
            "Reading Ref coordinates from %d partition .header/BED pairs "
            "with %d jobs",
            len(needed_partitions), args.jobs,
        )

        def target_task(
            partition: str,
        ) -> Tuple[str, Tuple[BlockTarget, ...]]:
            return read_partition_reference_targets(
                partition,
                partition_inputs[partition],
                reference_haplotypes,
            )

        for partition, targets in _bounded_thread_map(
            target_task, needed_partitions, args.jobs,
        ):
            targets_by_partition[partition] = targets

        if getattr(args, "reuse_alignments_only", False):
            for partition in needed_partitions:
                reference_dir = os.path.join(
                    work_root, "references", partition,
                )
                reference_fasta = os.path.join(
                    reference_dir, "reference.fa",
                )
                reference_marker = os.path.join(
                    reference_dir, "_SUCCESS.json",
                )
                if not os.path.isfile(reference_marker):
                    raise FileNotFoundError(reference_marker)
                with open(reference_marker, "rt") as handle:
                    reference_saved = json.load(handle)
                if targets_by_partition[partition] and not os.path.isfile(
                    reference_fasta
                ):
                    raise FileNotFoundError(reference_fasta)
                references[partition] = PreparedBlockReference(
                    partition,
                    reference_fasta,
                    tuple(targets_by_partition[partition]),
                    str(reference_saved["signature"]),
                )
        else:
            reference_sequences = load_reference_sequences(
                targets_by_partition, assembly_sources, args.jobs,
            )

            def materialize_task(partition: str) -> PreparedBlockReference:
                return materialize_partition_reference(
                    partition,
                    targets_by_partition[partition],
                    reference_sequences,
                    work_root,
                    options,
                )

            for reference in _bounded_thread_map(
                materialize_task, needed_partitions, args.jobs,
            ):
                references[reference.mapped_partition] = reference
            del reference_sequences

    mapping_candidates = [
        candidate
        for candidate in candidates
        if any(
            partition in references and references[partition].targets
            for partition in candidate.lifted_partitions
        )
    ]
    candidate_sequences = load_candidate_sequences(
        mapping_candidates,
        candidate_header_records,
        assembly_sources,
        int(args.input_anchor),
        args.jobs,
    )
    candidates_by_partition: Dict[str, List[MappingCandidate]] = (
        defaultdict(list)
    )
    for candidate in mapping_candidates:
        for partition in sorted(candidate.lifted_partitions):
            reference = references.get(partition)
            if reference is not None and reference.targets:
                candidates_by_partition[partition].append(candidate)
    mapping_partitions = sorted(candidates_by_partition)
    LOG.info(
        "Running %d independent partition jobs for %d candidate/partition "
        "assignments with "
        "%d jobs x %d cores",
        len(mapping_partitions),
        sum(len(rows) for rows in candidates_by_partition.values()),
        args.jobs,
        args.cores,
    )

    def partition_task(partition: str) -> str:
        return map_partition_candidates(
            candidates_by_partition[partition],
            candidate_sequences,
            references[partition],
            work_root,
            options,
        )

    placements: List[MappedPlacement] = []
    reporting_components: List[ResolvedPlacementComponent] = []
    if mapping_partitions:
        progress_step = max(1, len(mapping_partitions) // 20)
        for completed, result_path in enumerate(_bounded_thread_map(
            partition_task, mapping_partitions, args.jobs,
        ), 1):
            placements.extend(read_placements((result_path,)))
            reporting_components.extend(read_reporting_components((
                os.path.join(
                    os.path.dirname(result_path),
                    "reporting_placements.tsv",
                ),
            )))
            if (
                completed == len(mapping_partitions)
                or completed % progress_step == 0
            ):
                LOG.info(
                    "Local-novel mapping progress: %d/%d partitions",
                    completed, len(mapping_partitions),
                )
    if getattr(args, "reuse_alignments_only", False):
        # These tables define the already completed cleaning and local->global
        # promotion.  Preserve their bytes and mtimes so downstream frozen
        # alignment signatures remain reusable.
        for frozen_name in (
            "all_local_novel_alignments.tsv", "local_novel.bed",
            "global_novel.bed",
        ):
            frozen_path = os.path.join(output_dir, frozen_name)
            if not os.path.isfile(frozen_path):
                raise FileNotFoundError(frozen_path)
        write_local_alignment_summary(candidates, placements, output_dir)
    else:
        export_mapping_results(
            candidates, placements, output_dir, args.merge_gap,
        )
    final_reporting = resolve_all_evidence(
        reporting_components_to_evidence(reporting_components)
    )
    flattened_reporting = [
        row
        for candidate_id in sorted(final_reporting)
        for row in final_reporting[candidate_id]
    ]
    write_reporting_components(
        os.path.join(
            output_dir, "all_local_novel_reporting_alignments.tsv",
        ),
        flattened_reporting,
    )
    if not getattr(args, "reuse_alignments_only", False):
        copy_global_novel(analysis_dir, output_dir)
    temporary = marker + ".tmp"
    with open(temporary, "wt") as out:
        json.dump(
            {
                "version": "local-novel-mapping-v10-paired-blocks",
                "placement_resolution": PLACEMENT_RESOLUTION_VERSION,
                "reporting_resolution": REPORTING_RESOLUTION_VERSION,
                "signature": signature,
                "candidates": len(candidates),
                "placements": len(placements),
                "reporting_components": len(flattened_reporting),
                "resolved_candidates": len({
                    row.candidate_id for row in placements
                }),
                "reference_haplotypes": sorted(reference_haplotypes),
            },
            out,
            sort_keys=True,
        )
        out.write("\n")
    os.replace(temporary, marker)
    LOG.info("Wrote local-novel mapping outputs under %s", output_dir)


def main(argv: Optional[Sequence[str]] = None) -> int:
    args = parse_args(argv)
    logging.basicConfig(
        level=getattr(logging, args.log_level),
        format="[%(asctime)s] %(levelname)s: %(message)s",
        datefmt="%Y-%m-%d %H:%M:%S",
    )
    try:
        run(args)
    except (
        OSError,
        ValueError,
        RuntimeError,
        subprocess.CalledProcessError,
    ) as error:
        LOG.error("%s", error)
        return 1
    return 0


if __name__ == "__main__":
    sys.exit(main())
