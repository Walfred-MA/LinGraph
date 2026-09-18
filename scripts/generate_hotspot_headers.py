#!/usr/bin/env python3
"""Generate authoritative partition ``.header`` files without FASTA bases.

Each assembly hotspot table is reduced independently against its FAI.  Workers
write compact disk shards, then a second bounded phase merges those shards into
``LOCALBLOCKS/PARTITION/PARTITION_samples.header`` in query-list order.  No assembly
sequence is read and the complete cohort header text is never held in RAM.
"""
from __future__ import annotations

import argparse
import dataclasses
import logging
import os
import shutil
import struct
import tempfile
from collections import defaultdict
from concurrent.futures import ProcessPoolExecutor, as_completed
from pathlib import Path
from typing import Dict, Iterable, List, Optional, Sequence, Tuple

from build_local_graphs import (
    SharedFastaIndex,
    _hotspot_record_prefix,
    parse_hotspot_regions,
    read_target_partitions,
)
from minsetref_core import mp_context
from partition_hotspot_headers import AssemblySource, read_query_sources


LOG = logging.getLogger("generate_hotspot_headers")
HEADER_INDEX = struct.Struct("<IIIQI")
MERGE_INDEX = struct.Struct("<IIIIQI")


@dataclasses.dataclass(frozen=True)
class HeaderShard:
    shard_id: int
    data_path: str
    index_path: str
    assemblies: int
    hotspots: int
    records: int


@dataclasses.dataclass(frozen=True)
class MergeResult:
    plan_id: int
    partitions: Tuple[Tuple[int, str], ...]
    records: int


def source_buckets(
    sources: Sequence[AssemblySource], jobs: int,
) -> List[Tuple[AssemblySource, ...]]:
    count = min(max(1, jobs), len(sources))
    buckets: List[List[AssemblySource]] = [[] for _ in range(count)]
    for index, source in enumerate(sources):
        buckets[index % count].append(source)
    return [tuple(bucket) for bucket in buckets if bucket]


def write_header_shard(
    shard_id: int,
    sources: Sequence[AssemblySource],
    target_names: Sequence[str],
    workdir: str,
    kmer_length: int,
    merge_distance: int,
    anchor: int,
) -> HeaderShard:
    data_path = os.path.join(workdir, f"headers.{shard_id:04d}.data")
    index_path = os.path.join(workdir, f"headers.{shard_id:04d}.index")
    hotspot_count = 0
    record_count = 0
    with open(data_path, "wb") as data, open(index_path, "wb") as index:
        for source in sources:
            fasta_index = SharedFastaIndex(source)
            _order, parsed, source_hotspots = parse_hotspot_regions(
                source,
                fasta_index,
                len(target_names) - 1,
                kmer_length,
                merge_distance,
                anchor,
            )
            hotspot_count += source_hotspots
            for target_index in sorted(parsed):
                partition = target_names[target_index]
                prefix = _hotspot_record_prefix(partition, target_index)
                sequence_index = 0
                for contig in sorted(parsed[target_index]):
                    for start, end in parsed[target_index][contig]:
                        sequence_index += 1
                        header = (
                            f">{prefix}_{source.name}_{sequence_index}\t"
                            f"{contig}:{start}-{end}\tpartition={partition}\t"
                            f"sample={source.name}\n"
                        ).encode("utf-8")
                        offset = data.tell()
                        data.write(header)
                        index.write(HEADER_INDEX.pack(
                            target_index,
                            source.order,
                            sequence_index,
                            offset,
                            len(header),
                        ))
                        record_count += 1
    return HeaderShard(
        shard_id, data_path, index_path, len(sources),
        hotspot_count, record_count,
    )


def iter_header_index(path: str):
    with open(path, "rb") as handle:
        while True:
            payload = handle.read(HEADER_INDEX.size)
            if not payload:
                return
            if len(payload) != HEADER_INDEX.size:
                raise ValueError(f"truncated header index: {path}")
            yield HEADER_INDEX.unpack(payload)


def header_is_current(path: str, prefix: str, newest_input_ns: int) -> bool:
    try:
        stat = os.stat(path)
        if stat.st_size <= 0 or stat.st_mtime_ns < newest_input_ns:
            return False
        with open(path, "rb") as handle:
            return handle.readline().startswith(f">{prefix}_".encode())
    except OSError:
        return False


def merge_plan(
    plan_id: int,
    plan_path: str,
    target_names: Sequence[str],
    shard_paths: Sequence[str],
    output_root: str,
    newest_input_ns: int,
    resume: bool,
) -> MergeResult:
    by_target: Dict[int, List[Tuple[int, int, int, int, int]]] = defaultdict(list)
    with open(plan_path, "rb") as plan:
        while True:
            payload = plan.read(MERGE_INDEX.size)
            if not payload:
                break
            if len(payload) != MERGE_INDEX.size:
                raise ValueError(f"truncated merge plan: {plan_path}")
            target, source, sequence, shard, offset, length = MERGE_INDEX.unpack(
                payload
            )
            by_target[target].append((source, sequence, shard, offset, length))

    descriptors = [os.open(path, os.O_RDONLY) for path in shard_paths]
    completed: List[Tuple[int, str]] = []
    records = 0
    try:
        for target_index in sorted(by_target):
            partition = target_names[target_index]
            directory = os.path.join(output_root, partition)
            output = os.path.join(directory, f"{partition}_samples.header")
            prefix = _hotspot_record_prefix(partition, target_index)
            rows = sorted(by_target[target_index])
            if resume and header_is_current(output, prefix, newest_input_ns):
                completed.append((target_index, output))
                records += len(rows)
                continue
            Path(directory).mkdir(parents=True, exist_ok=True)
            temporary = output + f".tmp.{os.getpid()}"
            try:
                with open(temporary, "wb") as destination:
                    for _source, _sequence, shard, offset, length in rows:
                        payload = os.pread(descriptors[shard], length, offset)
                        if len(payload) != length:
                            raise IOError(
                                f"short header-shard read at {shard}:{offset}+{length}"
                            )
                        destination.write(payload)
                    destination.flush()
                    os.fsync(destination.fileno())
                os.replace(temporary, output)
            finally:
                try:
                    os.remove(temporary)
                except FileNotFoundError:
                    pass
            completed.append((target_index, output))
            records += len(rows)
    finally:
        for descriptor in descriptors:
            os.close(descriptor)
    return MergeResult(plan_id, tuple(completed), records)


def newest_input_mtime(paths: Iterable[str]) -> int:
    return max(os.stat(path).st_mtime_ns for path in paths)


def atomic_write(path: str, text: str) -> None:
    Path(path).parent.mkdir(parents=True, exist_ok=True)
    temporary = path + f".tmp.{os.getpid()}"
    try:
        with open(temporary, "wt") as output:
            output.write(text)
            output.flush()
            os.fsync(output.fileno())
        os.replace(temporary, path)
    finally:
        try:
            os.remove(temporary)
        except FileNotFoundError:
            pass


def parse_args(argv: Optional[Sequence[str]] = None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--hotspot-dir", required=True)
    parser.add_argument("-q", "--query-paths", required=True)
    parser.add_argument("-T", "--target-list", required=True)
    parser.add_argument("-o", "--output-dir", required=True)
    parser.add_argument("--active-list", required=True)
    parser.add_argument("-j", "--jobs", type=int, default=16)
    parser.add_argument("--hotspot-kmer-length", type=int, default=31)
    parser.add_argument("--hotspot-merge-distance", type=int, default=30_000)
    parser.add_argument("--hotspot-anchor", type=int, default=15_000)
    parser.add_argument("--resume", action="store_true")
    parser.add_argument(
        "--log-level", choices=("DEBUG", "INFO", "WARNING", "ERROR"),
        default="INFO",
    )
    args = parser.parse_args(argv)
    if args.jobs < 1 or args.hotspot_kmer_length < 1:
        parser.error("--jobs and --hotspot-kmer-length must be positive")
    if args.hotspot_merge_distance < 0 or args.hotspot_anchor < 0:
        parser.error("merge distance and anchor must be nonnegative")
    return args


def run(args: argparse.Namespace) -> None:
    hotspot_dir = os.path.abspath(args.hotspot_dir)
    query_paths = os.path.abspath(args.query_paths)
    target_list = os.path.abspath(args.target_list)
    output_root = os.path.abspath(args.output_dir)
    active_list = os.path.abspath(args.active_list)
    targets = read_target_partitions(target_list)
    target_names = tuple([""] + [target.name for target in targets])
    sources = read_query_sources(query_paths, hotspot_dir)
    Path(output_root).mkdir(parents=True, exist_ok=True)
    workdir = tempfile.mkdtemp(prefix=".hotspot_headers.work.", dir=output_root)
    success = False
    try:
        buckets = source_buckets(sources, args.jobs)
        LOG.info(
            "Reducing %d assembly hotspot maps into %d header shards",
            len(sources), len(buckets),
        )
        shards: List[HeaderShard] = []
        with ProcessPoolExecutor(
            max_workers=len(buckets), mp_context=mp_context(),
        ) as executor:
            futures = [
                executor.submit(
                    write_header_shard,
                    shard_id,
                    bucket,
                    target_names,
                    workdir,
                    args.hotspot_kmer_length,
                    args.hotspot_merge_distance,
                    args.hotspot_anchor,
                )
                for shard_id, bucket in enumerate(buckets)
            ]
            for completed, future in enumerate(as_completed(futures), 1):
                shards.append(future.result())
                LOG.info("Header shard progress: %d/%d", completed, len(futures))
        shards.sort(key=lambda row: row.shard_id)

        merge_jobs = min(args.jobs, len(targets))
        plan_paths = [
            os.path.join(workdir, f"merge.{index:04d}.index")
            for index in range(merge_jobs)
        ]
        plan_handles = [open(path, "wb") for path in plan_paths]
        try:
            for shard in shards:
                for target, source, sequence, offset, length in iter_header_index(
                    shard.index_path
                ):
                    plan_handles[(target - 1) % merge_jobs].write(
                        MERGE_INDEX.pack(
                            target, source, sequence, shard.shard_id,
                            offset, length,
                        )
                    )
        finally:
            for handle in plan_handles:
                handle.close()

        newest = newest_input_mtime([
            query_paths,
            target_list,
            *(source.fai for source in sources),
            *(source.hotspot for source in sources),
        ])
        shard_paths = tuple(shard.data_path for shard in shards)
        results: List[MergeResult] = []
        with ProcessPoolExecutor(
            max_workers=merge_jobs, mp_context=mp_context(),
        ) as executor:
            futures = [
                executor.submit(
                    merge_plan,
                    plan_id,
                    plan_paths[plan_id],
                    target_names,
                    shard_paths,
                    output_root,
                    newest,
                    args.resume,
                )
                for plan_id in range(merge_jobs)
            ]
            for completed, future in enumerate(as_completed(futures), 1):
                results.append(future.result())
                LOG.info("Header merge progress: %d/%d", completed, len(futures))
        active = sorted(
            (item for result in results for item in result.partitions),
            key=lambda row: row[0],
        )
        if not active:
            raise ValueError("hotspot maps produced no active partition headers")
        atomic_write(active_list, "".join(path + "\n" for _index, path in active))
        LOG.info(
            "Generated %d partition headers containing %d records; "
            "assembly bases read: 0",
            len(active), sum(result.records for result in results),
        )
        success = True
    finally:
        if success:
            shutil.rmtree(workdir, ignore_errors=True)
        else:
            LOG.error("Retained failed header work directory: %s", workdir)


def main(argv: Optional[Sequence[str]] = None) -> int:
    args = parse_args(argv)
    logging.basicConfig(
        level=getattr(logging, args.log_level),
        format="[%(asctime)s] %(levelname)s: %(message)s",
        datefmt="%Y-%m-%d %H:%M:%S",
    )
    run(args)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
