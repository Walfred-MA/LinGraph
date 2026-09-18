#!/usr/bin/env python3
"""Materialize, summarize, and reversibly package partition graph paths."""
from __future__ import annotations

import argparse
import csv
import heapq
import hashlib
import io
import json
import logging
import mmap
import os
import re
import shutil
import struct
import sys
import tempfile
from collections import OrderedDict, defaultdict
from concurrent.futures import ProcessPoolExecutor, ThreadPoolExecutor, as_completed
from contextlib import nullcontext
from dataclasses import dataclass
from pathlib import Path
from typing import Dict, Iterable, List, Mapping, Optional, Sequence, Set, Tuple

from minsetref_core import IndexedFasta, mp_context, wrap_fasta
from partition_hotspot_headers import (
    AssemblySource,
    merge_and_anchor,
    partition_name_from_target_path,
    read_query_sources,
)


LOG = logging.getLogger("build_local_graphs")
_SHARD_INDEX = struct.Struct("<IIIQQQ")
_MERGE_INDEX = struct.Struct("<IIIIQQQ")
FULL_CONTIG_INTERVAL_CUTOFF = 100
_COMPLEMENT = str.maketrans(
    "ACGTRYKMSWBDHVNacgtrykmswbdhvn",
    "TGCAYRMKSWVHDBNtgcayrmkswvhdbn",
)


def _allow_large_csv_fields() -> None:
    limit = sys.maxsize
    while True:
        try:
            csv.field_size_limit(limit)
            return
        except OverflowError:
            limit //= 10


_allow_large_csv_fields()


def reverse_complement(sequence: str) -> str:
    return sequence.translate(_COMPLEMENT)[::-1]


@dataclass(frozen=True)
class TargetPartition:
    index: int
    name: str


class SharedFastaIndex:
    """FAI metadata shared in RAM; each fetch batch opens then closes FASTA."""

    def __init__(self, source: AssemblySource) -> None:
        if source.fasta.lower().endswith((".gz", ".bgz", ".bgzf")):
            raise ValueError(
                f"hotspot extraction requires uncompressed FASTA: {source.fasta}"
            )
        if not os.path.isfile(source.fasta):
            raise FileNotFoundError(source.fasta)
        self.source = source
        self.records: Dict[str, Tuple[int, int, int, int]] = {}
        with open(source.fai, "rt") as handle:
            for line_number, raw in enumerate(handle, 1):
                if not raw.strip():
                    continue
                fields = raw.rstrip("\n").split("\t")
                if len(fields) < 5:
                    raise ValueError(
                        f"{source.fai}:{line_number}: expected FAI5+"
                    )
                self.records[fields[0]] = tuple(map(int, fields[1:5]))
        if not self.records:
            raise ValueError(f"{source.fai}: empty FAI")

    def length(self, contig: str) -> int:
        try:
            return self.records[contig][0]
        except KeyError as error:
            raise KeyError(
                f"{contig!r} is absent from {self.source.fai}"
            ) from error

    def fetch_many(
        self, requests: Sequence[Tuple[str, int, int]],
    ) -> List[str]:
        descriptor = os.open(self.source.fasta, os.O_RDONLY)
        try:
            output: List[str] = []
            for contig, start, end in requests:
                length, offset, line_bases, line_width = self.records[contig]
                if start < 0 or end <= start or end > length:
                    raise ValueError(
                        f"invalid interval {contig}:{start}-{end}; "
                        f"contig length is {length}"
                    )
                sequence = bytearray()
                position = start
                while position < end:
                    line, column = divmod(position, line_bases)
                    size = min(end - position, line_bases - column)
                    chunk = os.pread(
                        descriptor,
                        size,
                        offset + line * line_width + column,
                    )
                    if len(chunk) != size:
                        raise IOError(
                            f"short read for {contig}:{start}-{end} from "
                            f"{self.source.fasta}"
                        )
                    sequence.extend(chunk)
                    position += size
                output.append(sequence.decode("ascii"))
            return output
        finally:
            os.close(descriptor)


def read_target_partitions(path: str) -> List[TargetPartition]:
    targets: List[TargetPartition] = []
    seen = set()
    with open(path, "rt") as handle:
        for index, raw in enumerate(handle, 1):
            target_path = raw.strip()
            if not target_path:
                raise ValueError(f"{path}:{index}: empty target-list row")
            name = partition_name_from_target_path(target_path)
            if not name or name in {".", ".."} or os.path.basename(name) != name:
                raise ValueError(f"{path}:{index}: unsafe partition name {name!r}")
            if name in seen:
                raise ValueError(f"{path}:{index}: duplicate partition {name!r}")
            seen.add(name)
            targets.append(TargetPartition(index, name))
    if not targets:
        raise ValueError(f"{path}: no target FASTA paths")
    return targets


def _hotspot_record_prefix(partition: str, target_index: int) -> str:
    del target_index
    match = re.match(r"^(g\d+)_", partition)
    if match is not None:
        return match.group(1)
    return partition


def parse_hotspot_regions(
    source: AssemblySource,
    fasta_index: SharedFastaIndex,
    target_count: int,
    kmer_length: int,
    merge_distance: int,
    anchor: int,
) -> Tuple[int, Dict[int, Dict[str, List[Tuple[int, int]]]], int]:
    raw_regions = defaultdict(lambda: defaultdict(list))
    hotspot_count = 0
    with open(source.hotspot, "rt") as handle:
        for line_number, raw in enumerate(handle, 1):
            if not raw.strip() or raw.startswith("#"):
                continue
            fields = raw.split()
            if len(fields) < 4:
                raise ValueError(
                    f"{source.hotspot}:{line_number}: expected "
                    "CONTIG TARGET START END"
                )
            contig = fields[0]
            if contig == "chrM":
                continue
            target_index = int(fields[1])
            hotspot_start, hotspot_end = int(fields[2]), int(fields[3])
            if (
                target_index < 1
                or target_index > target_count
                or hotspot_start <= 0
                or hotspot_end < hotspot_start
            ):
                raise ValueError(
                    f"{source.hotspot}:{line_number}: invalid hotspot row"
                )
            contig_length = fasta_index.length(contig)
            start = max(0, hotspot_start - kmer_length)
            end = min(contig_length, hotspot_end)
            if start < end:
                raw_regions[target_index][contig].append((start, end))
                hotspot_count += 1

    processed: Dict[int, Dict[str, List[Tuple[int, int]]]] = {}
    for target_index, by_contig in raw_regions.items():
        processed[target_index] = {
            contig: merge_and_anchor(
                intervals,
                fasta_index.length(contig),
                merge_distance,
                anchor,
            )
            for contig, intervals in by_contig.items()
        }
    return source.order, processed, hotspot_count


@dataclass(frozen=True)
class HotspotShardResult:
    shard_id: int
    data_path: str
    index_path: str
    samples: int
    hotspots: int
    records: int
    bases: int


@dataclass(frozen=True)
class HotspotMergeResult:
    plan_id: int
    partitions: int
    records: int
    bases: int


def _load_fasta_into_memory(
    source: AssemblySource,
    fasta_index: SharedFastaIndex,
) -> Dict[str, bytearray]:
    """Load one assembly as bytearrays, validating it against its FAI."""
    sequences: Dict[str, bytearray] = {}
    current: Optional[bytearray] = None
    with open(source.fasta, "rb") as handle:
        for line_number, raw in enumerate(handle, 1):
            if raw.startswith(b">"):
                header = raw[1:].strip().split(None, 1)
                if not header:
                    raise ValueError(
                        f"{source.fasta}:{line_number}: empty FASTA header"
                    )
                contig = header[0].decode("utf-8")
                if contig in sequences:
                    raise ValueError(
                        f"{source.fasta}:{line_number}: duplicate FASTA "
                        f"record {contig!r}"
                    )
                current = bytearray()
                sequences[contig] = current
                continue
            if raw in (b"\n", b"\r\n", b""):
                continue
            if current is None:
                raise ValueError(
                    f"{source.fasta}:{line_number}: sequence precedes header"
                )
            current.extend(raw.rstrip(b"\r\n"))

    for contig, expected in fasta_index.records.items():
        sequence = sequences.get(contig)
        if sequence is None:
            raise KeyError(f"{contig!r} is absent from {source.fasta}")
        if len(sequence) != expected[0]:
            raise ValueError(
                f"{source.fasta}: {contig!r} has {len(sequence)} bases but "
                f"{source.fai} reports {expected[0]}"
            )
    return sequences


def _write_wrapped_bytes(handle, sequence: bytearray, width: int = 80) -> None:
    for start in range(0, len(sequence), width):
        handle.write(sequence[start:start + width])
        handle.write(b"\n")


def _build_hotspot_shard(
    shard_id: int,
    sources: Sequence[AssemblySource],
    target_names: Sequence[str],
    workdir: str,
    kmer_length: int,
    merge_distance: int,
    anchor: int,
) -> HotspotShardResult:
    """Process a fixed sample bucket into one FASTA-data/index shard pair."""
    data_path = os.path.join(workdir, f"shard_{shard_id:04d}.fa")
    index_path = os.path.join(workdir, f"shard_{shard_id:04d}.idx")
    hotspots = 0
    records = 0
    bases = 0
    with open(data_path, "wb") as data, open(index_path, "wb") as index:
        for source in sources:
            fasta_index = SharedFastaIndex(source)
            _, parsed, source_hotspots = parse_hotspot_regions(
                source,
                fasta_index,
                len(target_names) - 1,
                kmer_length,
                merge_distance,
                anchor,
            )
            hotspots += source_hotspots
            sequences = _load_fasta_into_memory(source, fasta_index)

            for target_index, by_contig in parsed.items():
                target_name = target_names[target_index]
                prefix = _hotspot_record_prefix(target_name, target_index)
                requests = [
                    (contig, start, end)
                    for contig in sorted(by_contig)
                    for start, end in by_contig[contig]
                ]
                for sequence_index, (contig, start, end) in enumerate(
                    requests, 1,
                ):
                    try:
                        contig_sequence = sequences[contig]
                    except KeyError as error:
                        raise KeyError(
                            f"{contig!r} is absent from {source.fasta}"
                        ) from error
                    sequence = contig_sequence[start:end]
                    if len(sequence) != end - start:
                        raise IOError(
                            f"short in-memory read for {contig}:{start}-{end} "
                            f"from {source.fasta}"
                        )
                    offset = data.tell()
                    header = (
                        f">{prefix}_{source.name}_{sequence_index}\t"
                        f"{contig}:{start}-{end}\t"
                        f"partition={target_name}\n"
                    ).encode("utf-8")
                    data.write(header)
                    _write_wrapped_bytes(data, sequence)
                    byte_length = data.tell() - offset
                    index.write(_SHARD_INDEX.pack(
                        target_index,
                        source.order,
                        sequence_index,
                        offset,
                        byte_length,
                        len(sequence),
                    ))
                    records += 1
                    bases += len(sequence)

            # Release the complete assembly before this worker loads its next
            # sample; peak memory is therefore one assembly per worker.
            del sequences

    return HotspotShardResult(
        shard_id, data_path, index_path, len(sources), hotspots, records, bases,
    )


def _merge_hotspot_plan(
    plan_id: int,
    plan_path: str,
    target_names: Dict[int, str],
    shard_paths: Sequence[str],
    output_root: str,
) -> HotspotMergeResult:
    """Concatenate every partition assigned to one process-owned plan."""
    records_by_target: Dict[
        int, List[Tuple[int, int, int, int, int, int]]
    ] = defaultdict(list)
    with open(plan_path, "rb") as plan:
        while True:
            packed = plan.read(_MERGE_INDEX.size)
            if not packed:
                break
            if len(packed) != _MERGE_INDEX.size:
                raise ValueError(f"truncated hotspot merge plan {plan_path}")
            (
                target_index, source_order, sequence_index, shard_id,
                offset, byte_length, base_count,
            ) = _MERGE_INDEX.unpack(packed)
            if target_index not in target_names:
                raise ValueError(
                    f"unexpected target {target_index} in {plan_path}"
                )
            if shard_id >= len(shard_paths):
                raise ValueError(
                    f"invalid shard {shard_id} in {plan_path}"
                )
            records_by_target[target_index].append((
                source_order, sequence_index, shard_id,
                offset, byte_length, base_count,
            ))

    shard_handles = []
    shard_maps: List[Optional[mmap.mmap]] = []
    merged_records = 0
    merged_bases = 0
    try:
        # Each process opens the completed shards only once, then writes every
        # partition in its plan without returning sequence bytes to the parent.
        for shard_path in shard_paths:
            handle = open(shard_path, "rb")
            shard_handles.append(handle)
            if os.path.getsize(shard_path):
                shard_maps.append(mmap.mmap(
                    handle.fileno(), 0, access=mmap.ACCESS_READ,
                ))
            else:
                shard_maps.append(None)

        target_indices = sorted(records_by_target)
        progress_step = max(250, len(target_indices) // 10)
        for completed, target_index in enumerate(target_indices, 1):
            target_name = target_names[target_index]
            target_records = records_by_target[target_index]
            target_records.sort(key=lambda row: (row[0], row[1]))
            directory = os.path.join(output_root, target_name)
            Path(directory).mkdir(parents=True, exist_ok=True)
            output = os.path.join(directory, f"{target_name}_samples.fasta")
            temporary = output + f".tmp.merge-{plan_id}"
            try:
                with open(temporary, "wb") as out:
                    for (
                        _source_order, _sequence_index, shard_id,
                        offset, byte_length, base_count,
                    ) in target_records:
                        shard = shard_maps[shard_id]
                        if shard is None:
                            raise IOError(
                                f"empty data shard {shard_id} contains an index"
                            )
                        record = shard[offset:offset + byte_length]
                        if len(record) != byte_length:
                            raise IOError(
                                f"short shard read {shard_id}:{offset}+"
                                f"{byte_length}"
                            )
                        out.write(record)
                        merged_records += 1
                        merged_bases += base_count
                os.replace(temporary, output)
            finally:
                try:
                    os.remove(temporary)
                except FileNotFoundError:
                    pass
            if completed % progress_step == 0:
                LOG.info(
                    "Hotspot merge plan %d progress: %d/%d partitions",
                    plan_id + 1, completed, len(target_indices),
                )
    finally:
        for shard in shard_maps:
            if shard is not None:
                shard.close()
        for handle in shard_handles:
            handle.close()

    return HotspotMergeResult(
        plan_id, len(records_by_target), merged_records, merged_bases,
    )


def build_hotspot_graphs(args: argparse.Namespace) -> None:
    target_list = os.path.abspath(args.target_list)
    hotspot_dir = os.path.abspath(args.hotspot_dir)
    query_paths = os.path.abspath(args.query_paths)
    output_root = os.path.abspath(args.output_dir)
    targets = read_target_partitions(target_list)
    sources = read_query_sources(query_paths, hotspot_dir)
    Path(output_root).mkdir(parents=True, exist_ok=True)
    workdir = tempfile.mkdtemp(
        prefix=".build_local_graphs.hotspot.work.", dir=output_root,
    )
    success = False
    try:
        target_names = tuple([""] + [target.name for target in targets])
        source_buckets = [
            tuple(sources[shard_id::args.jobs])
            for shard_id in range(args.jobs)
        ]
        LOG.info(
            "Hotspot shard phase: %d samples -> %d worker shards; each worker "
            "loads at most one complete assembly at a time",
            len(sources), args.jobs,
        )
        shard_results: List[Optional[HotspotShardResult]] = [
            None for _ in range(args.jobs)
        ]
        completed_samples = 0
        with ProcessPoolExecutor(
            max_workers=args.jobs,
            mp_context=mp_context(),
        ) as executor:
            futures = {
                executor.submit(
                    _build_hotspot_shard,
                    shard_id,
                    source_buckets[shard_id],
                    target_names,
                    workdir,
                    args.hotspot_kmer_length,
                    args.hotspot_merge_distance,
                    args.hotspot_anchor,
                ): shard_id
                for shard_id in range(args.jobs)
            }
            for completed_shards, future in enumerate(
                as_completed(futures), 1,
            ):
                result = future.result()
                shard_results[result.shard_id] = result
                completed_samples += result.samples
                LOG.info(
                    "Hotspot shard progress: %d/%d shards, %d/%d samples",
                    completed_shards, args.jobs,
                    completed_samples, len(sources),
                )

        completed_results = [
            result for result in shard_results if result is not None
        ]
        if len(completed_results) != args.jobs:
            raise RuntimeError("not all hotspot shards completed")

        # Route the compact indexes into process-owned merge plans.  This keeps
        # the parent from retaining all records and prevents Python sorting and
        # record assembly from contending in one threaded interpreter.
        plan_paths = [
            os.path.join(workdir, f"merge_plan_{plan_id:04d}.idx")
            for plan_id in range(args.jobs)
        ]
        plan_handles = [open(path, "wb") for path in plan_paths]
        plan_target_indices = [set() for _ in range(args.jobs)]
        total_hotspots = 0
        sharded_records = 0
        sharded_bases = 0
        try:
            for result in completed_results:
                total_hotspots += result.hotspots
                sharded_records += result.records
                sharded_bases += result.bases
                with open(result.index_path, "rb") as index:
                    while True:
                        packed = index.read(_SHARD_INDEX.size)
                        if not packed:
                            break
                        if len(packed) != _SHARD_INDEX.size:
                            raise ValueError(
                                "truncated hotspot shard index "
                                f"{result.index_path}"
                            )
                        (
                            target_index, source_order, sequence_index,
                            offset, byte_length, base_count,
                        ) = _SHARD_INDEX.unpack(packed)
                        if target_index < 1 or target_index > len(targets):
                            raise ValueError(
                                f"invalid target {target_index} in "
                                f"{result.index_path}"
                            )
                        plan_id = (target_index - 1) % args.jobs
                        plan_target_indices[plan_id].add(target_index)
                        plan_handles[plan_id].write(_MERGE_INDEX.pack(
                            target_index, source_order, sequence_index,
                            result.shard_id, offset, byte_length, base_count,
                        ))
        finally:
            for handle in plan_handles:
                handle.close()

        active_target_count = sum(map(len, plan_target_indices))
        shard_paths = tuple(
            result.data_path for result in completed_results
        )
        LOG.info(
            "Hotspot concatenate phase: %d active partitions -> %d "
            "process-owned merge plans",
            active_target_count, args.jobs,
        )
        records = 0
        bases = 0
        completed_partitions = 0
        with ProcessPoolExecutor(
            max_workers=args.jobs,
            mp_context=mp_context(),
        ) as executor:
            futures = [
                executor.submit(
                    _merge_hotspot_plan,
                    plan_id,
                    plan_paths[plan_id],
                    {
                        target_index: target_names[target_index]
                        for target_index in plan_target_indices[plan_id]
                    },
                    shard_paths,
                    output_root,
                )
                for plan_id in range(args.jobs)
            ]
            for future in as_completed(futures):
                result = future.result()
                records += result.records
                bases += result.bases
                completed_partitions += result.partitions
                LOG.info(
                    "Hotspot concatenate progress: %d/%d partitions "
                    "(%d/%d plans)",
                    completed_partitions, active_target_count,
                    sum(done.done() for done in futures), len(futures),
                )

        if records != sharded_records or bases != sharded_bases:
            raise RuntimeError(
                "hotspot shard/merge count mismatch: "
                f"shards={sharded_records}/{sharded_bases}, "
                f"merged={records}/{bases}"
            )
        LOG.info(
            "Hotspot LocalGraphs complete: %d hotspots, %d/%d partitions, "
            "%d sequences, %d bp under %s",
            total_hotspots, active_target_count, len(targets), records, bases,
            output_root,
        )
        success = True
    finally:
        if success:
            LOG.info("Cleaning completed hotspot shard and merge-plan files")
            shutil.rmtree(workdir, ignore_errors=True)
        else:
            LOG.error("Retained failed hotspot shard directory: %s", workdir)


def discover_manifests(root: str) -> List[Tuple[str, str]]:
    output: List[Tuple[str, str]] = []
    with os.scandir(root) as entries:
        for entry in entries:
            if not entry.is_dir(follow_symlinks=False):
                continue
            partition = entry.name
            path = os.path.join(entry.path, f"{partition}.adjusted.tsv")
            if os.path.isfile(path):
                output.append((partition, path))
    output.sort()
    if not output:
        raise ValueError(f"{root}: no PARTITION/PARTITION.adjusted.tsv files")
    return output


def read_partition_sources(
    root: str, allow_missing_fastas: bool = False,
) -> Dict[str, str]:
    path = os.path.join(root, "partition_sources.tsv")
    if not os.path.isfile(path):
        raise FileNotFoundError(
            f"{path}: compact adjusted manifests require the partition "
            "source table from refine_partition_paths.py"
        )
    output: Dict[str, str] = {}
    with open(path, "rt", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        if reader.fieldnames != ["partition", "samples_fasta"]:
            raise ValueError(
                f"{path}: expected partition and samples_fasta columns"
            )
        for row in reader:
            partition = row["partition"]
            source = os.path.abspath(row["samples_fasta"])
            if partition in output:
                raise ValueError(f"{path}: duplicate partition {partition!r}")
            if not allow_missing_fastas and not os.path.isfile(source):
                raise FileNotFoundError(
                    f"{path}: source FASTA for {partition} is missing: {source}"
                )
            output[partition] = source
    return output


def read_manifest(path: str, partition: str) -> List[Dict[str, str]]:
    with open(path, "rt", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        required = {
            "partition", "path_order", "path_id", "role", "label",
            "source_kind", "source_haplotype", "source_contig",
            "source_start", "source_end", "strand", "segments",
            "origins",
        }
        if reader.fieldnames is None or not required.issubset(reader.fieldnames):
            raise ValueError(f"{path}: expected columns {sorted(required)}")
        rows = list(reader)
    if not rows:
        raise ValueError(f"{path}: empty adjusted path manifest")
    for row in rows:
        if row["partition"] != partition:
            raise ValueError(
                f"{path}: row partition {row['partition']!r} does not match "
                f"directory {partition!r}"
            )
    rows.sort(key=lambda row: int(row["path_order"]))
    orders = [int(row["path_order"]) for row in rows]
    if orders != list(range(1, len(rows) + 1)):
        raise ValueError(f"{path}: path_order must be consecutive from 1")
    role_rank = {"original": 0, "reference": 1, "alternative": 2, "novel": 3}
    if rows[0]["role"] != "original":
        raise ValueError(f"{path}: first adjusted path is not original")
    ranks = []
    for row in rows:
        role = row["role"]
        if role not in role_rank:
            raise ValueError(f"{path}: invalid role {row['role']!r}")
        ranks.append(role_rank[role])
    if ranks != sorted(ranks):
        raise ValueError(f"{path}: roles are not in original/reference/alternative/novel order")
    return rows


def public_path_type(row: Mapping[str, str]) -> str:
    """Return stable user-facing path terminology for FASTA/summary output."""
    role = row["role"]
    source_kind = row["source_kind"]
    if role == "original":
        return "original"
    if role == "reference":
        return "reference"
    if source_kind == "local_novel_annotated":
        return "decoy"
    if role == "novel" or source_kind == "global_novel_query":
        return "novel"
    return "alternative"


def public_path_id(row: Mapping[str, str]) -> str:
    """Return ``PARTITION_ID_CONTIG_START_END`` for a graph path."""
    return public_path_id_from_values(
        row["partition"], row["source_contig"],
        row["source_start"], row["source_end"],
    )


def public_path_id_from_values(
    partition: str, source_contig: str, source_start, source_end,
) -> str:
    start, end = int(source_start), int(source_end)
    if start < 0 or end <= start:
        raise ValueError(
            f"{partition}: invalid graph path interval {start}-{end}"
        )
    partition_fields = partition.split("_")
    partition_id = partition_fields[0]
    if partition_id == "merged" and len(partition_fields) > 1:
        partition_id = partition_fields[1]
    if (
        not partition_id
        or not source_contig
        or any(character.isspace() for character in source_contig)
        or "/" in source_contig
    ):
        raise ValueError(
            f"{partition}: cannot construct graph path ID from contig "
            f"{source_contig!r}"
        )
    return f"{partition_id}_{source_contig}_{start}_{end}"


GRAPHIC_SEGMENT_RE = re.compile(r"([><])([^:;< >]+):([^;<> ]+)")
GRAPHIC_OPERATION_RE = re.compile(r"(\d+)([HID=XMN])")


def graphic_mapping_coordinates(label: str) -> str:
    """Convert internal whole-reference graphic paths to explicit intervals."""
    if not label or label == ".":
        return "."
    coordinates: List[str] = []
    for match in GRAPHIC_SEGMENT_RE.finditer(label):
        marker, contig, cigar = match.groups()
        operations = GRAPHIC_OPERATION_RE.findall(cigar)
        if not operations or "".join(
            f"{length}{operation}" for length, operation in operations
        ) != cigar:
            continue
        reference_position = 0
        aligned_start: Optional[int] = None
        aligned_end: Optional[int] = None
        for length_text, operation in operations:
            length = int(length_text)
            if operation in {"=", "X", "M"}:
                if aligned_start is None:
                    aligned_start = reference_position
                reference_position += length
                aligned_end = reference_position
            elif operation in {"H", "D", "N"}:
                reference_position += length
            elif operation == "I":
                continue
        if aligned_start is None or aligned_end is None:
            continue
        if marker == ">":
            start, end, strand = aligned_start, aligned_end, "+"
        else:
            start = reference_position - aligned_end
            end = reference_position - aligned_start
            strand = "-"
        coordinates.append(f"{contig}:{start}-{end}{strand}")
    return ";".join(coordinates) if coordinates else "."


def graphic_mapping_path(label: str) -> str:
    """Return complete, ordered graphic-CIGAR placements from a label.

    Keeping the strand marker and contig with every CIGAR makes multi-locus
    paths self-contained.  The returned components are positionally paired
    with those returned by :func:`graphic_mapping_coordinates`.
    """
    if not label or label == ".":
        return "."
    paths: List[str] = []
    for match in GRAPHIC_SEGMENT_RE.finditer(label):
        marker, contig, cigar = match.groups()
        operations = GRAPHIC_OPERATION_RE.findall(cigar)
        if not operations or "".join(
            f"{length}{operation}" for length, operation in operations
        ) != cigar:
            continue
        paths.append(f"{marker}{contig}:{cigar}")
    return "".join(paths) if paths else "."


def graph_fasta_path(output_root: str, partition: str) -> str:
    return os.path.join(
        output_root, partition, f"{partition}_samples.fasta",
    )


def _has_current_header_format(path: str) -> bool:
    """Avoid resuming FASTAs that still use verbose key/value headers."""
    try:
        with open(path, "rt") as handle:
            for raw in handle:
                if raw.startswith(">"):
                    fields = raw[1:].split()
                    fields = [field for field in fields
                              if field not in {'sequence_role=imported_alternative', 'sequence_role=static_template'}]
                    return (
                        len(fields) == 6
                        and fields[-1] in {
                            "reference", "alternative", "decoy", "novel",
                        }
                        and not any("=" in field for field in fields[1:])
                    )
    except OSError:
        return False
    return False


COMPACT_SEGMENT_RE = re.compile(r"^(.+)_(\d+)_(\d+)([+-])$")


def _parse_segments(
    text: str,
    context: str,
    source_fasta: str,
) -> List[Dict[str, object]]:
    if not text or text == ".":
        return []
    output: List[Dict[str, object]] = []
    for index, encoded in enumerate(text.split(";"), 1):
        match = COMPACT_SEGMENT_RE.fullmatch(encoded)
        if match is None:
            raise ValueError(
                f"{context}: invalid compact segment {encoded!r}; expected "
                "RECORD_START_ENDSTRAND"
            )
        record, start_text, end_text, strand = match.groups()
        start, end = int(start_text), int(end_text)
        if start < 0 or end <= start or strand not in {"+", "-"}:
            raise ValueError(
                f"{context}: invalid segment {index}: {start}-{end}{strand}"
            )
        output.append({
            "fasta": source_fasta,
            "record": record,
            "start": start,
            "end": end,
            "strand": strand,
        })
    return output


def materialize_manifest_row(
    row: Mapping[str, str],
    manifest: str,
    partition: str,
    source_fasta: str,
    reader: IndexedFasta,
) -> Tuple[str, str]:
    """Materialize one non-original row with its complete public header."""
    from fixed_alternatives import STORED_TEMPLATE_ROLES, materialize_fixed
    if row['source_kind'] in STORED_TEMPLATE_ROLES:
        return materialize_fixed(row, partition)
    segments = _parse_segments(
        row["segments"], f"{manifest}:{row['path_id']}", source_fasta,
    )
    if not segments:
        if row["role"] == "reference":
            raise ValueError(
                f"{manifest}:{row['path_id']}: reference path has no "
                "partition-local self-cleaned segments"
            )
        raise ValueError(
            f"{manifest}:{row['path_id']}: materialized path has no segment"
        )
    pieces: List[str] = []
    locals_: List[str] = []
    for segment in segments:
        record = str(segment["record"])
        start, end = int(segment["start"]), int(segment["end"])
        if record not in reader.index:
            if row["role"] == "reference":
                raise ValueError(
                    f"{manifest}:{row['path_id']}: reference segment record "
                    f"{record!r} is absent from this partition's self-cleaned "
                    f"sample FASTA {source_fasta}; refusing to substitute the "
                    "cohort-level cleaned reference or original genome"
                )
            raise ValueError(
                f"{manifest}:{row['path_id']}: compact segment record "
                f"{record!r} is absent from partition sample FASTA "
                f"{source_fasta}; this path must be materialized from its "
                "source assembly. Supply -q/--query-paths while building "
                "the LocalGraphs package. The completed package remains "
                "self-contained and reconstruction does not need -q."
            )
        if end > reader.length(record):
            raise ValueError(
                f"{manifest}:{row['path_id']}: {record}:{start}-{end} "
                "exceeds sequence length"
            )
        sequence = reader.fetch(record, start, end)
        if segment["strand"] == "-":
            sequence = reverse_complement(sequence)
        pieces.append(sequence)
        locals_.append(f"{record}:{start}-{end}{segment['strand']}")
    sequence = "".join(pieces)
    if not sequence:
        raise ValueError(f"{manifest}:{row['path_id']}: empty path sequence")
    expected_length = int(row["source_end"]) - int(row["source_start"])
    if len(sequence) != expected_length:
        raise ValueError(
            f"{manifest}:{row['path_id']}: materialized sequence length "
            f"{len(sequence)} does not match source interval "
            f"{row['source_start']}-{row['source_end']} "
            f"({expected_length} bases)"
        )
    source_coordinate = (
        f"{row['source_haplotype']}:{row['source_contig']}:"
        f"{row['source_start']}-{row['source_end']}{row['strand']}"
    )
    mapped_coordinate = graphic_mapping_coordinates(row["label"])
    mapped_cigar = graphic_mapping_path(row["label"])
    header = (
        f">{public_path_id(row)} {source_coordinate} {mapped_coordinate} "
        f"{mapped_cigar} {';'.join(locals_)} {public_path_type(row)}"
    )
    return header, sequence


NC_TO_CHR = {
    f"NC_{accession:06d}.1": f"chr{chromosome}"
    for accession, chromosome in zip(
        range(60925, 60947), range(1, 23),
    )
}
NC_TO_CHR.update({"NC_060947.1": "chrX", "NC_060948.1": "chrY"})
CHR_TO_NC = {chromosome: accession for accession, chromosome in NC_TO_CHR.items()}


def _read_assembly_paths(path: str, required_haplotypes=None) -> Dict[str, Tuple[str, str]]:
    sources: Dict[str, Tuple[str, str]] = {}
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
            if required_haplotypes is not None and haplotype not in required_haplotypes:
                continue
            fasta = os.path.abspath(os.path.expanduser(fields[1]))
            fai = fasta + ".fai"
            if haplotype in sources:
                raise ValueError(
                    f"{path}:{line_number}: duplicate haplotype {haplotype!r}"
                )
            for required in (fasta, fai):
                if not os.path.isfile(required):
                    raise FileNotFoundError(required)
            sources[haplotype] = (fasta, fai)
    if not sources and required_haplotypes != set():
        raise ValueError(f"{path}: no assembly sources")
    return sources


def _assembly_record(reader: IndexedFasta, table_contig: str) -> str:
    candidates = [table_contig]
    alias = NC_TO_CHR.get(table_contig) or CHR_TO_NC.get(table_contig)
    if alias is not None:
        candidates.append(alias)
    for candidate in candidates:
        if candidate in reader.index:
            return candidate
    raise KeyError(
        f"neither {table_contig!r} nor its chr* alias occurs in {reader.path}"
    )


def materialize_assembly_template_row(
    row: Mapping[str, str],
    partition: str,
    reader: IndexedFasta,
    sequence_override: Optional[str] = None,
) -> Tuple[str, str]:
    """Materialize a selected non-main original template from its assembly."""
    table_contig = row["source_contig"]
    record = _assembly_record(reader, table_contig)
    start, end = int(row["source_start"]), int(row["source_end"])
    length = reader.length(record)
    if start < 0 or end <= start or end > length:
        raise ValueError(
            f"{partition}: {row['source_haplotype']}:{record}:{start}-{end} "
            f"exceeds assembly contig length {length}"
        )
    sequence = (
        reader.fetch(record, start, end)
        if sequence_override is None else sequence_override
    )
    if len(sequence) != end - start:
        raise IOError(
            f"{partition}: short assembly read for "
            f"{row['source_haplotype']}:{record}:{start}-{end}"
        )
    if row["strand"] == "-":
        sequence = reverse_complement(sequence)
    source_coordinate = (
        f"{row['source_haplotype']}:{table_contig}:{start}-{end}"
        f"{row['strand']}"
    )
    local_coordinate = f"{table_contig}:{start}-{end}{row['strand']}"
    header = (
        f">{public_path_id(row)} {source_coordinate} {local_coordinate} "
        f". {local_coordinate} reference"
    )
    return header, sequence


def materialize_assembly_manifest_row(
    row: Mapping[str, str],
    partition: str,
    reader: IndexedFasta,
    sequence_override: Optional[str] = None,
) -> Tuple[str, str]:
    """Materialize any adjusted path directly from source coordinates."""
    table_contig = row["source_contig"]
    record = _assembly_record(reader, table_contig)
    start, end = int(row["source_start"]), int(row["source_end"])
    length = reader.length(record)
    if start < 0 or end <= start or end > length:
        raise ValueError(
            f"{partition}: {row['source_haplotype']}:{record}:{start}-{end} "
            f"exceeds assembly contig length {length}"
        )
    sequence = (
        reader.fetch(record, start, end)
        if sequence_override is None else sequence_override
    )
    if len(sequence) != end - start:
        raise IOError(
            f"{partition}: short assembly read for "
            f"{row['source_haplotype']}:{record}:{start}-{end}"
        )
    if row["strand"] == "-":
        sequence = reverse_complement(sequence)
    source_coordinate = (
        f"{row['source_haplotype']}:{table_contig}:{start}-{end}"
        f"{row['strand']}"
    )
    local_coordinate = f"{table_contig}:{start}-{end}{row['strand']}"
    header = (
        f">{public_path_id(row)} {source_coordinate} "
        f"{graphic_mapping_coordinates(row['label'])} "
        f"{graphic_mapping_path(row['label'])} {local_coordinate} "
        f"{public_path_type(row)}"
    )
    return header, sequence


def build_one(
    item: Tuple[str, str],
    output_root: str,
    resume: bool,
    source_fasta: str,
) -> Tuple[str, int, int]:
    partition, manifest = item
    directory = os.path.join(output_root, partition)
    output = graph_fasta_path(output_root, partition)
    if (
        resume
        and os.path.isfile(output)
        and os.path.getsize(output) > 0
        and os.stat(output).st_mtime_ns >= os.stat(manifest).st_mtime_ns
        and _has_current_header_format(output)
    ):
        return partition, -1, -1
    rows = read_manifest(manifest, partition)
    Path(directory).mkdir(parents=True, exist_ok=True)
    temporary = output + ".tmp"
    sequence_count = 0
    base_count = 0
    try:
        with open(temporary, "wt") as out, IndexedFasta(
            source_fasta,
        ) as reader:
            for row in rows:
                if row["role"] == "original":
                    continue
                header, sequence = materialize_manifest_row(
                    row, manifest, partition, source_fasta, reader,
                )
                out.write(f"{header}\n{wrap_fasta(sequence)}\n")
                sequence_count += 1
                base_count += len(sequence)
        os.replace(temporary, output)
    finally:
        try:
            os.remove(temporary)
        except FileNotFoundError:
            pass
    return partition, sequence_count, base_count


SUMMARY_BASE_FIELDS = (
    "partition",
    "path_order",
    "type",
    "source_haplotype",
    "source_contig",
    "source_start",
    "source_end",
    "strand",
    "segments",
)
SUMMARY_FIELDS = SUMMARY_BASE_FIELDS + (
    "reference_slice",
    "sequence_sha256",
    "fasta_header",
)
PARTITION_CACHE_PACKAGE = "partition_caches.jsonl"


REFERENCE_SLICE_RE = re.compile(r"^(.+):(\d+)-(\d+)([+-])$")


def _data_row_count(path: str, header_prefix: str = "") -> int:
    count = 0
    with open(path, "rt") as handle:
        for raw in handle:
            if not raw.strip() or raw.startswith("#"):
                continue
            if header_prefix and raw.startswith(header_prefix):
                continue
            count += 1
    return count


def validated_partition_cache_text(
    graph_folder: str, partition: str,
) -> Optional[str]:
    """Return one current all-sample cache verbatim, or None if absent."""
    directory = os.path.join(graph_folder, partition)
    cache = os.path.join(directory, f"{partition}cache.json")
    if not os.path.isfile(cache):
        return None
    graph = os.path.join(directory, f"{partition}.FA")
    bed = os.path.join(directory, f"{partition}.bed")
    alignment = os.path.join(directory, f"{partition}_align.txt")
    initial = os.path.join(
        directory, f"{partition}_breakpoint_consistency.initial.tsv",
    )
    blocked = os.path.join(
        directory, f"{partition}_breakpoint_consistency.bed",
    )
    for path in (graph, bed, alignment, initial, blocked):
        if not os.path.isfile(path):
            raise FileNotFoundError(
                f"{partition}: cache exists but dependency is missing: {path}"
            )
    if not all(os.path.getsize(path) > 0 for path in (graph, bed, alignment, initial)):
        raise ValueError(
            f"{partition}: cache dependency graph/BED/alignment/initial is empty"
        )
    with open(cache, "rt", encoding="utf-8") as handle:
        cache_text = handle.read()
    if not cache_text:
        raise ValueError(f"{partition}: cache is empty: {cache}")
    try:
        cache_data = json.loads(cache_text)
    except json.JSONDecodeError as error:
        raise ValueError(f"{partition}: invalid cache JSON: {error}") from error
    metadata = cache_data.get("_minsetref") if isinstance(cache_data, dict) else None
    control = (
        metadata.get("cohort_assemblysmall")
        if isinstance(metadata, dict) else None
    )
    if (
        not isinstance(metadata, dict)
        or metadata.get("alignment_scope") != "all_samples"
        or not isinstance(control, dict)
        or control.get("scope") != "all_samples"
    ):
        raise ValueError(
            f"{partition}: cache was not learned and replayed with all samples"
        )
    alignment_rows = _data_row_count(alignment, "hotspot_index\t")
    blocked_rows = _data_row_count(blocked)
    try:
        cached_alignment_rows = int(metadata["alignment_row_count"])
        cached_blocked_rows = int(metadata["blocked_output_count"])
    except (KeyError, TypeError, ValueError) as error:
        raise ValueError(
            f"{partition}: cache lacks alignment/blocking dependency counts"
        ) from error
    if cached_alignment_rows != alignment_rows:
        raise ValueError(
            f"{partition}: cache records {cached_alignment_rows} alignment rows, "
            f"but {alignment} contains {alignment_rows}"
        )
    if cached_blocked_rows != blocked_rows:
        raise ValueError(
            f"{partition}: cache records {cached_blocked_rows} blocked rows, "
            f"but {blocked} contains {blocked_rows}"
        )
    graph_inputs_time = max(os.stat(path).st_mtime_ns for path in (
        graph, bed, alignment,
    ))
    initial_time = os.stat(initial).st_mtime_ns
    blocked_time = os.stat(blocked).st_mtime_ns
    cache_time = os.stat(cache).st_mtime_ns
    if not (
        initial_time >= graph_inputs_time
        and blocked_time >= max(graph_inputs_time, initial_time)
        and cache_time >= max(graph_inputs_time, initial_time, blocked_time)
    ):
        raise ValueError(
            f"{partition}: cache or breakpoint outputs are older than their "
            "graph/BED/alignment dependencies"
        )
    return cache_text


def write_partition_cache_package(
    manifests: Sequence[Tuple[str, str]],
    graph_folder: Optional[str],
    output_dir: str,
) -> int:
    """Pack current partition caches as partition-keyed exact text records."""
    destination = os.path.join(output_dir, PARTITION_CACHE_PACKAGE)
    temporary = destination + ".tmp"
    count = 0
    try:
        with open(temporary, "wt", encoding="utf-8") as output:
            if graph_folder:
                graph_folder = os.path.abspath(graph_folder)
                if not os.path.isdir(graph_folder):
                    raise NotADirectoryError(graph_folder)
                for partition, _manifest in manifests:
                    cache_text = validated_partition_cache_text(
                        graph_folder, partition,
                    )
                    if cache_text is None:
                        continue
                    output.write(json.dumps({
                        "partition": partition,
                        "cache_text": cache_text,
                    }, ensure_ascii=False, separators=(",", ":")) + "\n")
                    count += 1
        os.replace(temporary, destination)
    finally:
        try:
            os.remove(temporary)
        except FileNotFoundError:
            pass
    return count


@dataclass(frozen=True)
class PackageBatchResult:
    index: int
    summary_shard: str
    alternatives_shard: str
    partitions: int
    paths: int
    materialized: int
    alternatives: int
    bases: int
    reference_slices: int


@dataclass(frozen=True)
class GroupedAssemblyPackageTask:
    order: int
    partition: str
    row: Dict[str, str]
    template: bool


@dataclass(frozen=True)
class GroupedAssemblyPackageResult:
    index: int
    summary_shard: str
    alternatives_shard: str
    alternatives_index: str
    paths: int
    alternatives: int
    bases: int
    reference_slices: int
    full_contigs: int


def _package_summary_row(
    partition: str, row: Mapping[str, str],
) -> Dict[str, str]:
    return {
        "partition": partition,
        "path_order": row["path_order"],
        "type": public_path_type(row),
        "source_haplotype": row["source_haplotype"],
        "source_contig": row["source_contig"],
        "source_start": row["source_start"],
        "source_end": row["source_end"],
        "strand": row["strand"],
        "segments": row["segments"],
        "reference_slice": ".",
        "sequence_sha256": ".",
        "fasta_header": ".",
    }


def _summary_row_text(row: Mapping[str, str]) -> str:
    buffer = io.StringIO()
    csv.DictWriter(
        buffer,
        fieldnames=SUMMARY_FIELDS,
        delimiter="\t",
        lineterminator="\n",
    ).writerow(row)
    return buffer.getvalue()


def _ordered_shard_row(handle) -> Optional[Tuple[int, str]]:
    raw = handle.readline()
    if not raw:
        return None
    order_text, separator, payload = raw.partition("\t")
    if not separator:
        raise ValueError(f"invalid ordered summary shard row: {raw!r}")
    return int(order_text), payload


def _alternative_index_row(handle) -> Optional[Tuple[int, int, int]]:
    raw = handle.readline()
    if not raw:
        return None
    fields = raw.rstrip("\n").split("\t")
    if len(fields) != 3:
        raise ValueError(f"invalid alternatives shard index row: {raw!r}")
    return tuple(map(int, fields))  # type: ignore[return-value]


def _copy_file_range(
    source_descriptor: int, output, offset: int, length: int,
) -> None:
    copied = 0
    while copied < length:
        chunk = os.pread(
            source_descriptor, min(8 * 1024 * 1024, length - copied),
            offset + copied,
        )
        if not chunk:
            raise IOError(
                f"short alternatives shard read at {offset}+{length}"
            )
        output.write(chunk)
        copied += len(chunk)


def write_grouped_assembly_package_indexes(
    manifests: Sequence[Tuple[str, str]],
    assembly_sources: Mapping[str, Tuple[str, str]],
    output_dir: str,
    jobs: int,
    main_haplotype: Optional[str],
) -> Tuple[int, int, int, int, int]:
    """Package deleted partition FASTAs by processing each assembly once.

    Requests are grouped first by haplotype and then by resolved contig. At
    most ``FULL_CONTIG_INTERVAL_CUTOFF`` intervals use indexed reads; a denser
    contig is loaded once and released immediately after all of its requests.
    Worker shards carry global row indexes so final output remains byte-order
    deterministic by partition and path order.
    """
    summary_path = os.path.join(output_dir, "local_graphs.tsv")
    alternatives_path = os.path.join(output_dir, "alternatives.fasta")
    summary_tmp = summary_path + ".tmp"
    alternatives_tmp = alternatives_path + ".tmp"
    workdir = tempfile.mkdtemp(
        prefix=".build_local_graphs.package.work.", dir=output_dir,
    )
    grouped: Dict[str, List[GroupedAssemblyPackageTask]] = defaultdict(list)
    metadata_rows: List[Tuple[int, Dict[str, str]]] = []
    path_count = 0
    success = False
    from fixed_alternatives import STORED_TEMPLATE_ROLES, fixed_segment, materialize_fixed
    assembly_sources = dict(assembly_sources)

    for partition, manifest in manifests:
        rows = read_manifest(manifest, partition)
        has_refined = any(row["role"] == "reference" for row in rows)
        graph_rows = [
            row for row in rows
            if row["role"] != "original" or not has_refined
        ]
        graph_ids: Dict[str, str] = {}
        for row in graph_rows:
            record_id = public_path_id(row)
            previous = graph_ids.get(record_id)
            if previous is not None:
                raise ValueError(
                    f"{partition}: graph paths {previous} and "
                    f"{row['path_id']} both resolve to FASTA ID "
                    f"{record_id!r}"
                )
            graph_ids[record_id] = row["path_id"]

        for row in rows:
            path_count += 1
            store_original = (
                row["role"] == "original"
                and (not has_refined or public_path_id(row) not in graph_ids)
            )
            store_template = (
                row["role"] == "reference" or store_original
            )
            included = row["role"] not in {
                "original", "reference",
            } or store_template
            if not included:
                metadata_rows.append((
                    path_count, _package_summary_row(partition, row),
                ))
                continue
            haplotype = row["source_haplotype"]
            if row['source_kind'] in STORED_TEMPLATE_ROLES:
                segment = fixed_segment(row)
                haplotype = '__fixed_alternative__:' + segment['fasta']
                assembly_sources[haplotype] = (segment['fasta'], None)
            if haplotype not in assembly_sources:
                raise ValueError(
                    f"{partition}: adjusted path {haplotype}:"
                    f"{row['source_contig']}:"
                    f"{row['source_start']}-{row['source_end']} must be "
                    "read from its original assembly because the temporary "
                    "partition FASTA was deleted; supply that assembly "
                    "through -q/--query-paths"
                )
            grouped[haplotype].append(GroupedAssemblyPackageTask(
                path_count, partition, row, store_original,
            ))

    metadata_shard = os.path.join(workdir, "summary.metadata.tsv")
    with open(metadata_shard, "wt", newline="") as output:
        for order, summary_row in metadata_rows:
            output.write(f"{order}\t{_summary_row_text(summary_row)}")

    grouped_items = list(enumerate(sorted(grouped.items())))
    LOG.info(
        "Grouped %d assembly-backed paths by %d haplotypes; contigs with "
        "more than %d requested intervals will be loaded once into RAM",
        sum(len(tasks) for tasks in grouped.values()), len(grouped),
        FULL_CONTIG_INTERVAL_CUTOFF,
    )

    def build_haplotype_shard(item) -> GroupedAssemblyPackageResult:
        group_index, (haplotype, tasks) = item
        summary_shard = os.path.join(
            workdir, f"summary.assembly.{group_index:06d}.tsv",
        )
        alternatives_shard = os.path.join(
            workdir, f"alternatives.assembly.{group_index:06d}.fasta",
        )
        alternatives_index = alternatives_shard + ".idx"
        summary_lines: List[Tuple[int, str]] = []
        alternative_entries: List[Tuple[int, int, int]] = []
        alternative_count = base_count = reference_slice_count = 0
        full_contig_count = 0
        fasta, fai = assembly_sources[haplotype]
        source_reader = IndexedFasta(fasta, fai)
        try:
            by_record: Dict[
                str, List[GroupedAssemblyPackageTask]
            ] = defaultdict(list)
            for task in tasks:
                record = (fixed_segment(task.row)['record']
                          if task.row['source_kind'] in STORED_TEMPLATE_ROLES else
                          _assembly_record(source_reader, task.row["source_contig"]))
                by_record[record].append(task)

            with open(alternatives_shard, "wb") as alternatives:
                records = sorted(
                    by_record.items(),
                    key=lambda item: min(task.order for task in item[1]),
                )
                for record, contig_tasks in records:
                    contig_tasks.sort(key=lambda task: (
                        int(task.row["source_start"]),
                        int(task.row["source_end"]), task.order,
                    ))
                    full_sequence: Optional[str] = None
                    if len(contig_tasks) > FULL_CONTIG_INTERVAL_CUTOFF:
                        full_sequence = source_reader.sequence(record)
                        full_contig_count += 1
                    try:
                        for task in contig_tasks:
                            row = task.row
                            start = int(row["source_start"])
                            end = int(row["source_end"])
                            if row['source_kind'] in STORED_TEMPLATE_ROLES:
                                segment = fixed_segment(row)
                                start, end = segment['start'], segment['end']
                            raw_sequence = (
                                full_sequence[start:end]
                                if full_sequence is not None
                                else source_reader.fetch(record, start, end)
                            )
                            if row['source_kind'] in STORED_TEMPLATE_ROLES:
                                header, sequence = materialize_fixed(row, task.partition, source_reader)
                            elif task.template:
                                header, sequence = materialize_assembly_template_row(
                                    row, task.partition, source_reader,
                                    raw_sequence,
                                )
                            else:
                                header, sequence = materialize_assembly_manifest_row(
                                    row, task.partition, source_reader,
                                    raw_sequence,
                                )
                            summary_row = _package_summary_row(
                                task.partition, row,
                            )
                            encoded = (
                                f"{header}\n{wrap_fasta(sequence)}\n"
                            ).encode("utf-8")
                            offset = alternatives.tell()
                            alternatives.write(encoded)
                            alternative_entries.append((task.order, offset, len(encoded)))
                            alternative_count += 1
                            base_count += len(sequence)
                            summary_lines.append((
                                task.order, _summary_row_text(summary_row),
                            ))
                    finally:
                        # A dense contig is never retained while this worker
                        # advances to another contig or assembly.
                        del full_sequence
        finally:
            source_reader.close()

        summary_lines.sort(key=lambda item: item[0])
        with open(summary_shard, "wt", newline="") as output:
            for order, line in summary_lines:
                output.write(f"{order}\t{line}")
        alternative_entries.sort(key=lambda item: item[0])
        with open(alternatives_index, "wt") as output:
            for order, offset, length in alternative_entries:
                output.write(f"{order}\t{offset}\t{length}\n")
        return GroupedAssemblyPackageResult(
            group_index, summary_shard, alternatives_shard,
            alternatives_index, len(tasks), alternative_count, base_count,
            reference_slice_count, full_contig_count,
        )

    results: List[GroupedAssemblyPackageResult] = []
    try:
        if grouped_items:
            worker_count = min(max(1, jobs), len(grouped_items))
            with ThreadPoolExecutor(max_workers=worker_count) as executor:
                results = list(executor.map(
                    build_haplotype_shard, grouped_items,
                ))
        results.sort(key=lambda result: result.index)

        summary_shards = [metadata_shard] + [
            result.summary_shard for result in results
        ]
        summary_handles = [open(path, "rt") for path in summary_shards]
        try:
            heap: List[Tuple[int, int, str]] = []
            for shard_index, handle in enumerate(summary_handles):
                row = _ordered_shard_row(handle)
                if row is not None:
                    heapq.heappush(heap, (row[0], shard_index, row[1]))
            expected_order = 1
            with open(summary_tmp, "wt", newline="") as summary:
                writer = csv.DictWriter(
                    summary,
                    fieldnames=SUMMARY_FIELDS,
                    delimiter="\t",
                    lineterminator="\n",
                )
                writer.writeheader()
                while heap:
                    order, shard_index, payload = heapq.heappop(heap)
                    if order != expected_order:
                        raise ValueError(
                            "grouped package summary order mismatch: "
                            f"expected {expected_order}, observed {order}"
                        )
                    summary.write(payload)
                    expected_order += 1
                    row = _ordered_shard_row(summary_handles[shard_index])
                    if row is not None:
                        heapq.heappush(
                            heap, (row[0], shard_index, row[1]),
                        )
            if expected_order - 1 != path_count:
                raise ValueError(
                    "grouped package summary row count mismatch: "
                    f"expected {path_count}, wrote {expected_order - 1}"
                )
        finally:
            for handle in summary_handles:
                handle.close()

        index_handles = [
            open(result.alternatives_index, "rt") for result in results
        ]
        data_handles = [
            open(result.alternatives_shard, "rb") for result in results
        ]
        try:
            alternative_heap: List[Tuple[int, int, int, int]] = []
            for shard_index, handle in enumerate(index_handles):
                row = _alternative_index_row(handle)
                if row is not None:
                    heapq.heappush(alternative_heap, (
                        row[0], shard_index, row[1], row[2],
                    ))
            last_order = 0
            with open(alternatives_tmp, "wb") as alternatives:
                while alternative_heap:
                    order, shard_index, offset, length = heapq.heappop(
                        alternative_heap,
                    )
                    if order <= last_order:
                        raise ValueError(
                            "grouped alternatives output order is not unique"
                        )
                    _copy_file_range(
                        data_handles[shard_index].fileno(), alternatives,
                        offset, length,
                    )
                    last_order = order
                    row = _alternative_index_row(
                        index_handles[shard_index],
                    )
                    if row is not None:
                        heapq.heappush(alternative_heap, (
                            row[0], shard_index, row[1], row[2],
                        ))
        finally:
            for handle in index_handles + data_handles:
                handle.close()

        os.replace(summary_tmp, summary_path)
        os.replace(alternatives_tmp, alternatives_path)
        success = True
    finally:
        for temporary in (summary_tmp, alternatives_tmp):
            try:
                os.remove(temporary)
            except FileNotFoundError:
                pass
        if success:
            shutil.rmtree(workdir)
        else:
            LOG.error("Retained failed package work directory: %s", workdir)

    LOG.info(
        "Grouped assembly extraction loaded %d dense contigs in full "
        "(cutoff: more than %d intervals)",
        sum(result.full_contigs for result in results),
        FULL_CONTIG_INTERVAL_CUTOFF,
    )
    return (
        path_count,
        sum(result.paths for result in results),
        sum(result.alternatives for result in results),
        sum(result.bases for result in results),
        sum(result.reference_slices for result in results),
    )


def write_package_indexes(
    manifests: Sequence[Tuple[str, str]],
    partition_sources: Mapping[str, str],
    output_dir: str,
    jobs: int,
    query_paths: Optional[str] = None,
    reference_haplotype: Optional[str] = None,
    embed_all_templates: bool = False,
) -> Tuple[int, int, int, int, int]:
    """Write the consolidated graph-data files with bounded RAM.

    Every selected graph path is materialized and stored in alternatives.fasta,
    including paths from the main reference. Original templates are also
    packed, even when refined graph paths exist,
    so reconstruction can reproduce the original hotspot targets without
    access to any cohort assembly. Shared original/refined IDs are stored once.
    """
    summary_path = os.path.join(output_dir, "local_graphs.tsv")
    alternatives_path = os.path.join(output_dir, "alternatives.fasta")
    summary_tmp = summary_path + ".tmp"
    alternatives_tmp = alternatives_path + ".tmp"
    required_haplotypes = None
    if reference_haplotype is not None:
        from fixed_alternatives import STORED_TEMPLATE_ROLES
        required_haplotypes = set()
        required_haplotypes.update(
            row['source_haplotype'] for partition, manifest in manifests
            for row in read_manifest(manifest, partition)
            if row['source_kind'] not in STORED_TEMPLATE_ROLES
        )
    assembly_sources = (
        _read_assembly_paths(os.path.abspath(query_paths), required_haplotypes)
        if query_paths else {}
    )
    main_haplotype = (
        None if embed_all_templates
        else reference_haplotype or next(iter(assembly_sources), None)
    )
    # Integrated graph construction removes its temporary partition FASTAs.
    # In that common case, schedule extraction by assembly/contig instead of
    # making partition workers repeatedly jump among cohort assemblies.
    if assembly_sources and not any(
        os.path.isfile(source) for source in partition_sources.values()
    ):
        return write_grouped_assembly_package_indexes(
            manifests, assembly_sources, output_dir, jobs, main_haplotype,
        )

    workdir = tempfile.mkdtemp(
        prefix=".build_local_graphs.package.work.", dir=output_dir,
    )
    worker_count = min(max(1, jobs), len(manifests))
    batch_size = (len(manifests) + worker_count - 1) // worker_count
    batches = [
        manifests[start:start + batch_size]
        for start in range(0, len(manifests), batch_size)
    ]
    success = False

    def build_batch(item) -> PackageBatchResult:
        batch_index, batch = item
        summary_shard = os.path.join(
            workdir, f"summary.{batch_index:06d}.tsv",
        )
        alternatives_shard = os.path.join(
            workdir, f"alternatives.{batch_index:06d}.fasta",
        )
        partition_count = path_count = materialized_count = 0
        alternative_count = base_count = reference_slice_count = 0
        assembly_readers: "OrderedDict[str, IndexedFasta]" = OrderedDict()

        def assembly_reader(haplotype: str) -> IndexedFasta:
            if haplotype not in assembly_sources:
                raise KeyError(haplotype)
            source_reader = assembly_readers.get(haplotype)
            if source_reader is None:
                if len(assembly_readers) >= 4:
                    _old_haplotype, old_reader = assembly_readers.popitem(
                        last=False,
                    )
                    old_reader.close()
                fasta, fai = assembly_sources[haplotype]
                source_reader = IndexedFasta(fasta, fai)
                assembly_readers[haplotype] = source_reader
            else:
                assembly_readers.move_to_end(haplotype)
            return source_reader

        try:
            summary = open(summary_shard, "wt", newline="")
            alternatives = open(alternatives_shard, "wt")
            with summary, alternatives:
                writer = csv.DictWriter(
                    summary,
                    fieldnames=SUMMARY_FIELDS,
                    delimiter="\t",
                    lineterminator="\n",
                )
                for partition, manifest in batch:
                    rows = read_manifest(manifest, partition)
                    has_refined = any(
                        row["role"] == "reference" for row in rows
                    )
                    graph_rows = [
                        row for row in rows
                        if row["role"] != "original" or not has_refined
                    ]
                    graph_ids: Dict[str, str] = {}
                    for row in graph_rows:
                        record_id = public_path_id(row)
                        previous = graph_ids.get(record_id)
                        if previous is not None:
                            raise ValueError(
                                f"{partition}: graph paths {previous} and "
                                f"{row['path_id']} both resolve to FASTA ID "
                                f"{record_id!r}"
                            )
                        graph_ids[record_id] = row["path_id"]
                    source_fasta = partition_sources[partition]
                    source_context = (
                        IndexedFasta(source_fasta)
                        if os.path.isfile(source_fasta)
                        else nullcontext(None)
                    )
                    with source_context as reader:
                        for row in rows:
                            store_original = (
                                row["role"] == "original"
                                and (not has_refined or public_path_id(row) not in graph_ids)
                            )
                            store_template = (
                                row["role"] == "reference"
                                or store_original
                            )
                            included = row["role"] not in {
                                "original", "reference",
                            } or store_template
                            summary_row = {
                                "partition": partition,
                                "path_order": row["path_order"],
                                "type": public_path_type(row),
                                "source_haplotype": row["source_haplotype"],
                                "source_contig": row["source_contig"],
                                "source_start": row["source_start"],
                                "source_end": row["source_end"],
                                "strand": row["strand"],
                                "segments": row["segments"],
                                "reference_slice": ".",
                                "sequence_sha256": ".",
                                "fasta_header": ".",
                            }
                            path_count += 1
                            if not included:
                                writer.writerow(summary_row)
                                continue
                            if row['source_kind'] in {'imported_alternative', 'static_template'}:
                                from fixed_alternatives import materialize_fixed
                                header, sequence = materialize_fixed(row, partition)
                            elif store_original or reader is None:
                                haplotype = row["source_haplotype"]
                                if haplotype not in assembly_sources:
                                    raise ValueError(
                                        f"{partition}: adjusted path "
                                        f"{haplotype}:"
                                        f"{row['source_contig']}:"
                                        f"{row['source_start']}-"
                                        f"{row['source_end']} must be read "
                                        "from its original assembly because "
                                        "the temporary partition FASTA was "
                                        "deleted; supply that assembly "
                                        "through -q/--query-paths"
                                    )
                                source_reader = assembly_reader(haplotype)
                                if store_original:
                                    header, sequence = materialize_assembly_template_row(
                                        row, partition, source_reader,
                                    )
                                else:
                                    header, sequence = materialize_assembly_manifest_row(
                                        row, partition, source_reader,
                                    )
                            else:
                                header, sequence = materialize_manifest_row(
                                    row, manifest, partition,
                                    source_fasta, reader,
                                )
                            materialized_count += 1
                            alternatives.write(f"{header}\n{wrap_fasta(sequence)}\n")
                            alternative_count += 1
                            base_count += len(sequence)
                            writer.writerow(summary_row)
                    partition_count += 1
        finally:
            for assembly_reader in assembly_readers.values():
                assembly_reader.close()
        return PackageBatchResult(
            batch_index, summary_shard, alternatives_shard,
            partition_count, path_count, materialized_count,
            alternative_count, base_count, reference_slice_count,
        )

    try:
        with ThreadPoolExecutor(max_workers=worker_count) as executor:
            results = list(executor.map(build_batch, enumerate(batches)))
        results.sort(key=lambda row: row.index)
        with open(summary_tmp, "wt", newline="") as summary:
            writer = csv.DictWriter(
                summary,
                fieldnames=SUMMARY_FIELDS,
                delimiter="\t",
                lineterminator="\n",
            )
            writer.writeheader()
            for result in results:
                with open(result.summary_shard, "rt") as shard:
                    shutil.copyfileobj(shard, summary, 1024 * 1024)
        with open(alternatives_tmp, "wt") as alternatives:
            for result in results:
                with open(result.alternatives_shard, "rt") as shard:
                    shutil.copyfileobj(shard, alternatives, 1024 * 1024)
        os.replace(summary_tmp, summary_path)
        os.replace(alternatives_tmp, alternatives_path)
        success = True
    finally:
        for temporary in (summary_tmp, alternatives_tmp):
            try:
                os.remove(temporary)
            except FileNotFoundError:
                pass
        if success:
            shutil.rmtree(workdir)
        else:
            LOG.error("Retained failed package work directory: %s", workdir)
    return (
        sum(row.paths for row in results),
        sum(row.materialized for row in results),
        sum(row.alternatives for row in results),
        sum(row.bases for row in results),
        sum(row.reference_slices for row in results),
    )


def index_fasta_byte_ranges(path: str) -> Dict[str, Tuple[int, int]]:
    """Index whole-record byte ranges without creating a third package file."""
    ranges: Dict[str, Tuple[int, int]] = {}
    current_id: Optional[str] = None
    current_start = 0
    with open(path, "rb") as handle:
        while True:
            offset = handle.tell()
            raw = handle.readline()
            if not raw:
                if current_id is not None:
                    ranges[current_id] = (current_start, offset)
                break
            if not raw.startswith(b">"):
                continue
            if current_id is not None:
                ranges[current_id] = (current_start, offset)
            current_id = raw[1:].split(None, 1)[0].decode("utf-8")
            if current_id in ranges:
                raise ValueError(f"{path}: duplicate FASTA ID {current_id!r}")
            current_start = offset
    return ranges


def validate_embedded_templates(
    package_dir: str, reference_haplotype: Optional[str],
) -> None:
    """Require graph templates AND original hotspot templates to be packed."""
    summary_path = os.path.join(package_dir, "local_graphs.tsv")
    alternatives_path = os.path.join(package_dir, "alternatives.fasta")
    by_partition: DefaultDict[str, List[Dict[str, str]]] = defaultdict(list)
    with open(summary_path, "rt", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        for row in reader:
            by_partition[row["partition"]].append(row)
    stored = set(index_fasta_byte_ranges(alternatives_path))
    missing: List[str] = []
    for partition, rows in by_partition.items():
        for row in rows:
            if row["type"] not in {'original', 'reference'}:
                continue
            if (
                row['type'] == "original"
                and row["source_haplotype"] == reference_haplotype
            ):
                continue
            record_id = public_path_id_from_values(
                partition, row["source_contig"],
                row["source_start"], row["source_end"],
            )
            if (
                record_id not in stored
                and not any(
                    other.get('reference_slice', '.') not in {'', '.'}
                    and public_path_id_from_values(
                        partition, other['source_contig'], other['source_start'], other['source_end'],
                    ) == record_id
                    for other in rows
                )
            ):
                missing.append(
                    f"{partition}:{row['source_haplotype']}:{record_id}"
                )
    if missing:
        raise ValueError(
            "incomplete self-contained LocalGraphs package: original/selected "
            "partition templates are absent from alternatives.fasta and do "
            "not have a verified reference_slice: "
            + ",".join(missing[:20])
            + (f" (and {len(missing) - 20} more)" if len(missing) > 20 else "")
        )


def expand_two_file_package(
    package_dir: str,
    output_dir: str,
    reference_fasta: str,
    jobs: int,
    resume: bool,
) -> None:
    """Reconstruct annotated per-partition FASTAs from the graph package."""
    summary_path = os.path.join(package_dir, "local_graphs.tsv")
    alternatives_path = os.path.join(package_dir, "alternatives.fasta")
    for path in (summary_path, alternatives_path):
        if not os.path.isfile(path):
            raise FileNotFoundError(path)
    by_partition: Dict[str, List[Dict[str, str]]] = defaultdict(list)
    with open(summary_path, "rt", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        if (
            reader.fieldnames is None
            or not set(SUMMARY_BASE_FIELDS).issubset(reader.fieldnames)
        ):
            raise ValueError(
                f"{summary_path}: expected columns "
                f"{list(SUMMARY_BASE_FIELDS)}"
            )
        for row in reader:
            by_partition[row["partition"]].append(row)
    required_alternatives = {
        public_path_id_from_values(
            partition, row["source_contig"],
            row["source_start"], row["source_end"],
        )
        for partition, rows in by_partition.items()
        for row in rows
        if (
            row["type"] not in {"original", "reference"}
            and row.get("reference_slice", ".") in {"", "."}
        )
    }
    optional_templates = {
        public_path_id_from_values(
            partition, row["source_contig"],
            row["source_start"], row["source_end"],
        )
        for partition, rows in by_partition.items()
        for row in rows if row['type'] in {'original', 'reference'}
    }
    ranges = index_fasta_byte_ranges(alternatives_path)
    allowed_alternatives = required_alternatives | optional_templates
    if (
        not required_alternatives.issubset(ranges)
        or not set(ranges).issubset(allowed_alternatives)
    ):
        missing = required_alternatives - set(ranges)
        extra = set(ranges) - allowed_alternatives
        raise ValueError(
            "alternatives.fasta/TSV mismatch; missing="
            + ",".join(sorted(missing)[:10])
            + "; extra=" + ",".join(sorted(extra)[:10])
        )
    Path(output_dir).mkdir(parents=True, exist_ok=True)
    reference_fasta = os.path.abspath(reference_fasta)
    if not os.path.isfile(reference_fasta):
        raise FileNotFoundError(reference_fasta)

    def expand_partition(item) -> Tuple[str, int]:
        partition, rows = item
        rows.sort(key=lambda row: int(row["path_order"]))
        has_reference = any(row["type"] == "reference" for row in rows)
        materialized_rows = [
            row for row in rows
            if row["type"] != "original" or not has_reference
        ]
        directory = os.path.join(output_dir, partition)
        output = graph_fasta_path(output_dir, partition)
        if not materialized_rows:
            return partition, 0
        if (
            resume and os.path.isfile(output) and os.path.getsize(output) > 0
            and os.stat(output).st_mtime_ns >= os.stat(summary_path).st_mtime_ns
            and os.stat(output).st_mtime_ns >= os.stat(
                alternatives_path
            ).st_mtime_ns
            and _has_current_header_format(output)
        ):
            return partition, -1
        Path(directory).mkdir(parents=True, exist_ok=True)
        temporary = output + ".tmp"
        count = 0
        try:
            with open(temporary, "wt") as output_handle, IndexedFasta(
                reference_fasta,
            ) as reference_reader, open(
                alternatives_path, "rb",
            ) as alternatives:
                for row in materialized_rows:
                    record_id = public_path_id_from_values(
                        partition, row["source_contig"],
                        row["source_start"], row["source_end"],
                    )
                    if record_id in ranges:
                        start, end = ranges[record_id]
                        alternatives.seek(start)
                        remaining = end - start
                        while remaining:
                            chunk = alternatives.read(min(1024 * 1024, remaining))
                            if not chunk:
                                raise IOError(
                                    f"short read from {alternatives_path}"
                                )
                            output_handle.write(chunk.decode("ascii"))
                            remaining -= len(chunk)
                    elif row.get("reference_slice", ".") not in {"", "."}:
                        match = REFERENCE_SLICE_RE.fullmatch(
                            row["reference_slice"],
                        )
                        if match is None:
                            raise ValueError(
                                f"{summary_path}: invalid reference_slice "
                                f"{row['reference_slice']!r}"
                            )
                        table_contig, start_text, end_text, strand = match.groups()
                        start, end = int(start_text), int(end_text)
                        contig = _assembly_record(reference_reader, table_contig)
                        if end > reference_reader.length(contig):
                            raise ValueError(
                                f"{summary_path}: reference_slice "
                                f"{table_contig}:{start}-{end} exceeds supplied "
                                "reference"
                            )
                        sequence = reference_reader.fetch(contig, start, end)
                        if strand == "-":
                            sequence = reverse_complement(sequence)
                        expected_digest = row.get("sequence_sha256", ".")
                        observed_digest = hashlib.sha256(
                            sequence.encode("ascii")
                        ).hexdigest()
                        if (
                            expected_digest in {None, "", "."}
                            or observed_digest != expected_digest
                        ):
                            raise ValueError(
                                f"{summary_path}: supplied reference does not "
                                f"match compacted path {record_id}"
                            )
                        stored_header = row.get("fasta_header", ".")
                        if stored_header in {"", "."}:
                            raise ValueError(
                                f"{summary_path}: compacted path {record_id} "
                                "lacks fasta_header"
                            )
                        if stored_header.split(None, 1)[0] != record_id:
                            raise ValueError(
                                f"{summary_path}: compacted header ID does not "
                                f"match {record_id}"
                            )
                        output_handle.write(
                            f">{stored_header}\n{wrap_fasta(sequence)}\n"
                        )
                    else:
                        contig = row["source_contig"]
                        start = int(row["source_start"])
                        end = int(row["source_end"])
                        strand = row["strand"]
                        if end > reference_reader.length(contig):
                            raise ValueError(
                                f"{summary_path}: {contig}:{start}-{end} "
                                "exceeds supplied reference"
                            )
                        sequence = reference_reader.fetch(contig, start, end)
                        if strand == "-":
                            sequence = reverse_complement(sequence)
                        record_id = public_path_id_from_values(
                            partition, row["source_contig"],
                            row["source_start"], row["source_end"],
                        )
                        contig_length = reference_reader.length(contig)
                        if strand == "+":
                            graphic_cigar = (
                                f">{contig}:{start}H{end - start}="
                                f"{contig_length - end}H"
                            )
                        else:
                            graphic_cigar = (
                                f"<{contig}:{contig_length - end}H"
                                f"{end - start}={start}H"
                            )
                        header = (
                            f">{record_id} {row['source_haplotype']}:"
                            f"{contig}:{start}-{end}{strand} "
                            f"{contig}:{start}-{end}{strand} "
                            f"{graphic_cigar} . reference"
                        )
                        output_handle.write(
                            f"{header}\n{wrap_fasta(sequence)}\n"
                        )
                    count += 1
            os.replace(temporary, output)
        finally:
            try:
                os.remove(temporary)
            except FileNotFoundError:
                pass
        return partition, count

    items = sorted(by_partition.items())
    with ThreadPoolExecutor(max_workers=max(1, jobs)) as executor:
        results = list(executor.map(expand_partition, items))
    LOG.info(
        "Expanded %d partitions (%d newly written paths) into %s",
        len(items), sum(max(0, count) for _partition, count in results),
        output_dir,
    )


def remove_legacy_partition_fastas(
    output_dir: str,
) -> int:
    """Leave the package files, removing only canonical old directories."""
    removed = 0
    allowed = {
        "local_graphs.tsv", "alternatives.fasta", "Graphs.list",
        "alternative_intervals.bed",
        "alternatives.fasta.fai", "partition_targets.list",
        "partition_targets.list.bin",
        PARTITION_CACHE_PACKAGE,
        "summary.complete",
        "reference_covered_mapped_novel.tsv",
    }
    with os.scandir(output_dir) as entries:
        existing = list(entries)
    for entry in existing:
        if entry.name in allowed:
            continue
        if not entry.is_dir(follow_symlinks=False):
            raise ValueError(
                f"package output directory contains an unexpected file: "
                f"{entry.path}"
            )
        partition = entry.name
        directory = entry.path
        expected = f"{partition}_samples.fasta"
        entries = os.listdir(directory)
        unexpected = [name for name in entries if name != expected]
        if unexpected:
            raise ValueError(
                f"refusing to remove legacy graph directory {directory}; "
                "it contains unexpected files: "
                + ",".join(sorted(unexpected)[:20])
            )
        fasta = os.path.join(directory, expected)
        if os.path.isfile(fasta):
            os.remove(fasta)
        os.rmdir(directory)
        removed += 1
    return removed


def remove_stale_package_workdirs(output_dir: str) -> int:
    """Remove only work directories created by failed package builds."""
    prefix = ".build_local_graphs.package.work."
    removed = 0
    with os.scandir(output_dir) as entries:
        stale = [
            entry.path for entry in entries
            if entry.name.startswith(prefix)
            and entry.is_dir(follow_symlinks=False)
        ]
    for path in stale:
        shutil.rmtree(path)
        removed += 1
    return removed


def parse_args(argv: Optional[Sequence[str]] = None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=(
            "Build the adjusted-path LocalGraphs package, "
            "expand it into partition FASTAs, or build directly from "
            "KmerSearcher hotspots."
        )
    )
    mode = parser.add_mutually_exclusive_group(required=True)
    mode.add_argument(
        "-i", "--adjusted-dir",
        help="output directory from refine_partition_paths.py",
    )
    mode.add_argument(
        "--hotspot-dir",
        help="directory containing NAME_hotspot.txt files",
    )
    mode.add_argument(
        "--expand-package",
        help=(
            "restore per-partition FASTAs from a folder containing "
            "local_graphs.tsv and alternatives.fasta; the exclusion audit "
            "file, when present, is ignored during expansion"
        ),
    )
    parser.add_argument(
        "-q", "--query-paths",
        help=(
            "NAME FASTA (with adjacent FASTA.fai) list; used for KmerSearcher queries, or with "
            "--adjusted-dir only as a fallback when partition FASTAs are missing"
        ),
    )
    parser.add_argument(
        "--reference-haplotype",
        help=(
            "main-reference label for calling-interval metadata; all path "
            "sequences are stored. With --query-paths, defaults to its first entry"
        ),
    )
    parser.add_argument(
        "--graph-folder",
        help=(
            "partition root: reuse existing PARTITION.FA and PARTITION.fasta "
            "without reading assemblies; also pack current partition caches. "
            "Defaults to Graphs/ beside the adjusted directory when present"
        ),
    )
    parser.add_argument(
        "--graph-list",
        help="finalized partition list; defaults to <graph-folder>.list when present",
    )
    parser.add_argument(
        "--embed-all-templates", action="store_true",
        help=(
            "treat every original template as a calling interval; use when "
            "partitions have independent references rather than one main one. "
            "All sequence bases are stored regardless of this option"
        ),
    )
    parser.add_argument(
        "-T", "--target-list",
        help="the exact <GraphFolder>.list used with KmerSearcher -T",
    )
    parser.add_argument("--hotspot-kmer-length", type=int, default=31)
    parser.add_argument("--hotspot-merge-distance", type=int, default=30_000)
    parser.add_argument("--hotspot-anchor", type=int, default=15_000)
    parser.add_argument(
        "-o", "--output-dir", default="LocalGraphs",
    )
    parser.add_argument(
        "-r", "--reference-fasta",
        help="reference FASTA required with --expand-package",
    )
    parser.add_argument("-j", "--jobs", type=int, default=16)
    parser.add_argument("--resume", action="store_true")
    parser.add_argument(
        "--log-level", choices=("DEBUG", "INFO", "WARNING", "ERROR"),
        default="INFO",
    )
    args = parser.parse_args(argv)
    if args.jobs < 1:
        parser.error("--jobs must be positive")
    if args.hotspot_dir and (not args.query_paths or not args.target_list):
        parser.error("--hotspot-dir requires --query-paths and --target-list")
    if args.target_list and not args.hotspot_dir:
        parser.error("--target-list requires --hotspot-dir")
    if args.expand_package and args.query_paths:
        parser.error("--query-paths is not used with --expand-package")
    if args.reference_haplotype and not args.adjusted_dir:
        parser.error(
            "--reference-haplotype requires --adjusted-dir"
        )
    if args.graph_folder and not args.adjusted_dir:
        parser.error("--graph-folder requires --adjusted-dir")
    if args.graph_list and not args.adjusted_dir:
        parser.error("--graph-list requires --adjusted-dir")
    if args.embed_all_templates and not args.adjusted_dir:
        parser.error("--embed-all-templates requires --adjusted-dir")
    if args.embed_all_templates and args.reference_haplotype:
        parser.error(
            "--embed-all-templates and --reference-haplotype are mutually exclusive"
        )
    if (
        args.adjusted_dir
        and not args.query_paths
        and not args.reference_haplotype
        and not args.embed_all_templates
    ):
        parser.error(
            "--adjusted-dir without --query-paths requires "
            "--reference-haplotype (for example, CHM13_h1) to identify "
            "the main reference in calling-interval metadata"
        )
    if args.expand_package and not args.reference_fasta:
        parser.error("--expand-package requires --reference-fasta")
    if not args.expand_package and args.reference_fasta:
        parser.error("--reference-fasta requires --expand-package")
    for name in (
        "hotspot_kmer_length", "hotspot_merge_distance", "hotspot_anchor",
    ):
        if getattr(args, name) < 0:
            parser.error(f"--{name.replace('_', '-')} cannot be negative")
    if args.hotspot_kmer_length == 0:
        parser.error("--hotspot-kmer-length must be positive")
    return args


def run(args: argparse.Namespace) -> None:
    if args.expand_package:
        expand_two_file_package(
            os.path.abspath(args.expand_package),
            os.path.abspath(args.output_dir),
            os.path.abspath(args.reference_fasta),
            args.jobs,
            args.resume,
        )
        return
    if args.hotspot_dir:
        build_hotspot_graphs(args)
        return
    adjusted_dir = os.path.abspath(args.adjusted_dir)
    output_dir = os.path.abspath(args.output_dir)
    manifests = discover_manifests(adjusted_dir)
    graph_folder = args.graph_folder
    if graph_folder is None:
        candidate = os.path.join(os.path.dirname(adjusted_dir), "Graphs")
        if os.path.isdir(candidate):
            graph_folder = candidate
    if graph_folder is not None:
        graph_folder = os.path.abspath(os.path.expanduser(graph_folder))
    graph_listing = args.graph_list
    if graph_listing is None and graph_folder:
        graph_listing = next((path for path in (
            graph_folder + '.list',
            os.path.join(os.path.dirname(graph_folder), 'summary', 'Graphs.list'),
        ) if os.path.isfile(path)), None)
    if graph_listing:
        from graph_package_fastas import select_manifests
        manifests = select_manifests(manifests, graph_listing)
    partition_sources_path = os.path.join(
        adjusted_dir, "partition_sources.tsv",
    )
    if os.path.isfile(partition_sources_path):
        partition_sources = read_partition_sources(
            adjusted_dir, allow_missing_fastas=bool(args.query_paths or graph_folder),
        )
    elif args.query_paths or graph_folder:
        # Integrated graph construction deliberately deletes materialized
        # partition FASTAs. With query paths, adjusted rows can be fetched
        # directly from their authoritative assembly coordinates.
        partition_sources = {
            partition: "" for partition, _manifest in manifests
        }
    else:
        raise FileNotFoundError(partition_sources_path)
    missing_sources = {
        partition for partition, _path in manifests
        if partition not in partition_sources
    }
    if missing_sources and not graph_folder:
        raise ValueError(
            "partition_sources.tsv lacks adjusted partitions: "
            + ",".join(sorted(missing_sources)[:20])
        )
    for partition in missing_sources:
        partition_sources[partition] = ""
    Path(output_dir).mkdir(parents=True, exist_ok=True)
    stale_workdirs = remove_stale_package_workdirs(output_dir)
    if stale_workdirs:
        LOG.info(
            "Removed %d stale failed package work directories",
            stale_workdirs,
        )
    package_reference = None if args.embed_all_templates else args.reference_haplotype
    if package_reference is None and not args.embed_all_templates and args.query_paths:
        # Read provenance only. Existing graph FASTAs must remain exportable
        # even after the assembly paths in this list have disappeared.
        with open(args.query_paths) as handle:
            package_reference = next((raw.split()[0] for raw in handle
                                      if raw.strip() and not raw.lstrip().startswith('#')), None)
    writer = write_package_indexes
    writer_options = {}
    if graph_folder:
        from graph_package_fastas import write_graph_package
        writer = write_graph_package
        writer_options['graph_folder'] = graph_folder
    (
        path_count, materialized_count, alternative_count,
        alternative_bases, reference_slice_count,
    ) = writer(
        manifests, partition_sources, output_dir, args.jobs,
        args.query_paths, package_reference,
        args.embed_all_templates, **writer_options,
    )
    validate_embedded_templates(output_dir, package_reference)
    cache_count = write_partition_cache_package(
        manifests, graph_folder, output_dir,
    )
    covered_source = os.path.join(
        adjusted_dir, "reference_covered_mapped_novel.tsv",
    )
    covered_target = os.path.join(
        output_dir, "reference_covered_mapped_novel.tsv",
    )
    temporary_covered = covered_target + ".tmp"
    try:
        if os.path.isfile(covered_source):
            shutil.copyfile(covered_source, temporary_covered)
        else:
            with open(temporary_covered, "wt") as output:
                output.write(
                    "partition\tquery_kind\tquery_id\trecord_id\t"
                    "haplotype\tsource_contig\tsource_start\tsource_end\t"
                    "local_start\tlocal_end\tref_haplotype\tref_contig\t"
                    "ref_start\tref_end\tstrand\treason\n"
                )
        os.replace(temporary_covered, covered_target)
    finally:
        try:
            os.remove(temporary_covered)
        except FileNotFoundError:
            pass
    from alternative_intervals import BED_NAME, write_bed
    bed_count = write_bed(
        os.path.join(output_dir, 'local_graphs.tsv'),
        os.path.join(output_dir, BED_NAME),
        package_reference,
        args.embed_all_templates,
    )
    LOG.info('Wrote %d original alternative/novel calling intervals', bed_count)
    if graph_listing:
        destination = Path(output_dir) / 'Graphs.list'
        temporary = destination.with_suffix('.list.tmp')
        temporary.write_text(''.join(partition + '\n' for partition, _ in manifests))
        os.replace(temporary, destination)
    removed = remove_legacy_partition_fastas(output_dir)
    LOG.info(
        "LocalGraphs package complete: local_graphs.tsv summarizes %d "
        "paths (%d materialized); alternatives.fasta contains %d "
        "stored paths and %d bp; %d exact paths use reference_slice; "
        "partition_caches.jsonl contains %d current "
        "all-sample caches; removed %d legacy partition folders",
        path_count, materialized_count, alternative_count,
        alternative_bases, reference_slice_count, cache_count, removed,
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
