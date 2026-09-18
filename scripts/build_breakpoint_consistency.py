#!/usr/bin/env python3
"""Build cohort breakpoint caches and blocked BEDs for finished partitions.

The input partition list is the same list used by local graph alignment.  Its
paths are used only to recover ordered partition names; every artifact is
resolved below ``--graph-folder`` so a migrated graph root remains resumable.

For each partition with a nonempty ``PARTITION.FA`` and
``PARTITION_align.txt`` this command:

1. derives conservative initial breakpoint segments;
2. constructs an in-memory gfixbreaks graphDB from those segments;
3. learns one Assemblysmall control with all selected rows present;
4. replays that frozen control on every selected row; and
5. atomically writes the blocked BED6 and then ``PARTITIONcache.json`` with
   the learned control, making the cache the final completion artifact.

All samples and all source BED haplotypes participate by default.  ``--onlyref``
restricts both to the named reference haplotype(s). Missing alignment folders
are reported and skipped without deleting anything.
"""

from __future__ import annotations

import argparse
import collections
import concurrent.futures
import dataclasses
import json
import logging
import multiprocessing
import os
from pathlib import Path
from types import ModuleType
from typing import Dict, Iterable, List, Optional, Sequence, Tuple

from block_partition_alignments import (
    COHORT_ASSEMBLYSMALL_KEY,
    COHORT_ASSEMBLYSMALL_VERSION,
    block_one_alignment,
    lossless_projection_results,
    result_within_block_control_paths,
    result_within_block_links,
    write_output,
)
from build_folder_uniformbreaks import (
    CACHE_FORMAT,
    FolderTask,
    uniformbreaks_from_blocks,
    write_graphdb_atomic,
)
from partition_hotspot_segments import build_partition_segments
from summarize_partition_hotspot_segments import (
    read_alignment_output,
    select_alignment_rows,
    write_final_segments,
)
from uniform_graph_blocks import load_graphfixbreaks, read_initial_blocks


LOG = logging.getLogger("build_breakpoint_consistency")
_GFIXBREAKS: Optional[ModuleType] = None
_REFERENCE_HAPLOTYPES: Tuple[str, ...] = ("CHM13_h1",)
_ONLY_REFERENCE = False
_MERGE_SMALL = 20_000


@dataclasses.dataclass(frozen=True)
class BreakpointTask:
    partition: str
    directory: str
    graph: str
    bed: str
    alignment: str
    initial: str
    blocks: str
    cache: str
    graphfixbreaks: str
    resume: bool


@dataclasses.dataclass(frozen=True)
class BreakpointResult:
    partition: str
    status: str
    alignment_rows: int = 0
    segment_count: int = 0
    block_count: int = 0
    cache: str = ""
    blocks: str = ""
    message: str = ""


def _init_worker(
    graphfixbreaks: str,
    reference_haplotypes: Sequence[str],
    only_reference: bool,
    merge_small: int,
) -> None:
    global _GFIXBREAKS, _REFERENCE_HAPLOTYPES, _ONLY_REFERENCE, _MERGE_SMALL
    _GFIXBREAKS = load_graphfixbreaks(graphfixbreaks)
    _REFERENCE_HAPLOTYPES = tuple(reference_haplotypes)
    _ONLY_REFERENCE = bool(only_reference)
    _MERGE_SMALL = int(merge_small)


def _nonempty_file(path: str) -> bool:
    try:
        return os.path.isfile(path) and os.path.getsize(path) > 0
    except OSError:
        return False


def _cache_metadata(path: str) -> Optional[Dict[str, object]]:
    try:
        with open(path, "rt", encoding="utf-8") as handle:
            data = json.load(handle)
        metadata = data.get("_minsetref")
        return metadata if isinstance(metadata, dict) else None
    except (OSError, ValueError, TypeError, json.JSONDecodeError):
        return None


def _outputs_are_current(task: BreakpointTask) -> bool:
    if not task.resume:
        return False
    if not (
        _nonempty_file(task.initial)
        and os.path.isfile(task.blocks)
        and _nonempty_file(task.cache)
    ):
        return False
    metadata = _cache_metadata(task.cache)
    expected_scope = "only_reference" if _ONLY_REFERENCE else "all_samples"
    control = metadata.get(COHORT_ASSEMBLYSMALL_KEY) if metadata else None
    try:
        control_cutoff = int(control.get("cutoff", -1))
    except (AttributeError, TypeError, ValueError):
        return False
    if (
        metadata is None
        or metadata.get("format") != CACHE_FORMAT
        or metadata.get("alignment_scope") != expected_scope
        or not isinstance(control, dict)
        or control.get("version") != COHORT_ASSEMBLYSMALL_VERSION
        or control.get("scope") != expected_scope
        or control_cutoff != _MERGE_SMALL
    ):
        return False
    try:
        inputs_time = max(os.path.getmtime(path) for path in (
            task.graph, task.bed, task.alignment,
        ))
        initial_time = os.path.getmtime(task.initial)
        blocks_time = os.path.getmtime(task.blocks)
        cache_time = os.path.getmtime(task.cache)
    except OSError:
        return False
    return (
        initial_time >= inputs_time
        and blocks_time >= max(inputs_time, initial_time)
        and cache_time >= max(inputs_time, initial_time, blocks_time)
    )


def _base_cache_is_current(
    task: BreakpointTask, metadata: Optional[Dict[str, object]], row_count: int,
) -> bool:
    """Allow old v5 caches to be upgraded without rebuilding graphDB."""
    expected_scope = "only_reference" if _ONLY_REFERENCE else "all_samples"
    try:
        cached_row_count = int(
            metadata.get("alignment_row_count", -1)
            if metadata is not None else -1
        )
    except (TypeError, ValueError):
        return False
    if (
        metadata is None
        or metadata.get("format") != CACHE_FORMAT
        or metadata.get("alignment_scope") != expected_scope
        or cached_row_count != row_count
        or not _nonempty_file(task.initial)
        or not _nonempty_file(task.cache)
    ):
        return False
    try:
        inputs_time = max(os.path.getmtime(path) for path in (
            task.graph, task.bed, task.alignment,
        ))
        initial_time = os.path.getmtime(task.initial)
        cache_time = os.path.getmtime(task.cache)
    except OSError:
        return False
    return initial_time >= inputs_time and cache_time >= initial_time


def learn_cohort_assemblysmall_control(
    module: ModuleType,
    graphdb,
    rows,
    cutoff: int,
):
    """Learn one fixed Assemblysmall control with every row present."""
    cohort_results = {}
    for input_index, row in enumerate(rows):
        if not row.graph_path or row.graph_path == "*":
            continue
        row_results = lossless_projection_results(row, graphdb, module)
        for result_index, regions in enumerate(row_results.values()):
            cohort_results[
                f"{input_index}:{result_index}:{row.query_name}"
            ] = regions
    changed = True
    while changed:
        cohort_results, changed_loci = module.assemblysmall(
            cohort_results, cutoff,
        )
        changed = bool(changed_loci)
    return (
        result_within_block_links(cohort_results),
        result_within_block_control_paths(cohort_results),
    )


def declare_graphdb_coordinate_system(graphdb) -> None:
    """Make fresh and reloaded graphDB objects interpret regions identically."""
    # set_validregions_from_annotations() writes 0-based half-open endpoints.
    # A freshly constructed graphDB has no cache metadata yet, however, and
    # gfixbreaks interprets an undeclared endpoint as legacy inclusive.  The
    # serialized cache is explicitly marked half-open by write_graphdb_atomic,
    # so failing to set the same marker before the first blocking pass can make
    # its committed BED differ by one base when that cache is later replayed.
    graphdb.cache_metadata["coordinate_system"] = "0-based-half-open"


def _process_task(task: BreakpointTask) -> BreakpointResult:
    if not _nonempty_file(task.graph):
        return BreakpointResult(
            task.partition, "skipped_missing_fa", message=task.graph,
        )
    if not _nonempty_file(task.alignment):
        return BreakpointResult(
            task.partition, "skipped_missing_alignment",
            message=task.alignment,
        )
    if not _nonempty_file(task.bed):
        return BreakpointResult(
            task.partition, "skipped_missing_bed", message=task.bed,
        )
    if _outputs_are_current(task):
        metadata = _cache_metadata(task.cache) or {}
        return BreakpointResult(
            task.partition, "reused",
            int(metadata.get("alignment_row_count", 0)),
            int(metadata.get("block_count", 0)),
            int(metadata.get("blocked_output_count", 0)),
            task.cache, task.blocks,
        )
    if _GFIXBREAKS is None:
        raise RuntimeError("breakpoint-consistency worker was not initialized")

    all_rows = read_alignment_output(task.alignment)
    if any(row.hotspot_index != 1 for row in all_rows):
        bad = next(row.hotspot_index for row in all_rows if row.hotspot_index != 1)
        raise ValueError(
            f"{task.partition}: single-partition alignment contains hotspot "
            f"index {bad}, expected 1"
        )
    rows = select_alignment_rows(
        all_rows, _REFERENCE_HAPLOTYPES, _ONLY_REFERENCE,
    )
    if not rows:
        return BreakpointResult(
            task.partition,
            "skipped_no_reference_rows" if _ONLY_REFERENCE else "skipped_no_rows",
        )

    old_metadata = _cache_metadata(task.cache)
    if _base_cache_is_current(task, old_metadata, len(rows)):
        initial_blocks = read_initial_blocks(task.initial)
        graphdb = _GFIXBREAKS.graphDB.load_json(task.cache)
        raw_lengths = old_metadata.get("graph_path_lengths", {})
        if not isinstance(raw_lengths, dict):
            raw_lengths = {}
        path_lengths = {
            str(name): int(length) for name, length in raw_lengths.items()
        }
        usable_block_count = int(old_metadata.get(
            "usable_reference_block_count", len(initial_blocks),
        ))
        skipped_block_count = int(old_metadata.get(
            "skipped_reference_block_count", 0,
        ))
        uniform_region_count = int(old_metadata.get(
            "uniform_reference_region_count", len(initial_blocks),
        ))
        segment_count = len(initial_blocks)
    else:
        segments = build_partition_segments(
            rows,
            task.graph,
            task.bed,
            _REFERENCE_HAPLOTYPES,
            only_reference=_ONLY_REFERENCE,
        )
        write_final_segments(task.initial, segments)
        if not segments:
            return BreakpointResult(
                task.partition, "skipped_no_segments",
                alignment_rows=len(rows), message=task.initial,
            )
        initial_blocks = read_initial_blocks(task.initial)
        (
            graphdb, path_lengths, usable_block_count, skipped_block_count,
            uniform_region_count,
        ) = uniformbreaks_from_blocks(
            _GFIXBREAKS, rows, initial_blocks, task.graph,
        )
        segment_count = len(segments)

    # Blocking must use the same coordinate interpretation before and after
    # the graphDB is serialized.  This is essential for fresh graphDB objects,
    # whose metadata dictionary starts empty.
    declare_graphdb_coordinate_system(graphdb)

    # Pass 1 presents all selected rows to one Assemblysmall call. Its
    # bindrecords can therefore resolve boundary ties using cohort-wide
    # support rather than independently per row. Synthetic keys preserve
    # duplicate query names without affecting the learned graph paths.
    learned_links, learned_paths = learn_cohort_assemblysmall_control(
        _GFIXBREAKS, graphdb, rows, _MERGE_SMALL,
    )

    # Pass 2 applies only the frozen cohort control to every row, including
    # the reference. It must not make any row-specific Assemblysmall choices.
    blocked = []
    for input_index, row in enumerate(rows):
        row_blocks, _links = block_one_alignment(
            _GFIXBREAKS, graphdb, task.partition, input_index, row,
            _MERGE_SMALL,
            reference_merge_links=learned_links,
            reference_block_paths=learned_paths,
            fixed_merge_control=True,
        )
        blocked.extend(row_blocks)
    write_output(task.blocks, blocked, include_header=False)

    # Cache is deliberately committed last. Final target selection treats a
    # current cache as proof that both segment derivation and blocking ended.
    cache_task = FolderTask(
        hotspot_index=1,
        graph_name=task.partition,
        graph_fasta=task.graph,
        row_count=len(rows),
        block_count=len(initial_blocks),
        output=task.cache,
        source_alignment=task.alignment,
        source_blocks=task.initial,
        graphfixbreaks=task.graphfixbreaks,
        resume=False,
        alignment_scope=(
            "only_reference" if _ONLY_REFERENCE else "all_samples"
        ),
        blocked_output=task.blocks,
        blocked_output_count=len(blocked),
        assemblysmall_scope=(
            "only_reference" if _ONLY_REFERENCE else "all_samples"
        ),
        assemblysmall_cutoff=_MERGE_SMALL,
        assemblysmall_links=tuple(sorted(learned_links)),
        assemblysmall_paths=tuple(sorted(learned_paths)),
    )
    write_graphdb_atomic(
        _GFIXBREAKS, graphdb, cache_task, path_lengths,
        usable_block_count, skipped_block_count, uniform_region_count,
    )
    return BreakpointResult(
        task.partition, "built", len(rows), segment_count, len(blocked),
        task.cache, task.blocks,
    )


def partition_name(value: str) -> str:
    token = value.strip().split()[0]
    path = Path(token.rstrip("/"))
    name = path.name
    parent = path.parent.name
    # Graphs.active.list contains a header manifest inside its same-named
    # partition directory.  Recognize only that exact layout so a names-only
    # partition that legitimately ends in ".header" is not rewritten.
    for suffix in ("_samples.header", ".header"):
        if name.endswith(suffix):
            candidate = name[:-len(suffix)]
            if parent and parent == candidate:
                return parent
    suffixes = (
        "_samples.fasta", "_samples.fa", ".fasta", ".fna", ".FA", ".fa",
    )
    stem = name
    for suffix in suffixes:
        if name.endswith(suffix):
            stem = name[:-len(suffix)]
            break
    return parent if parent and parent == stem else stem


def read_partition_names(path: str) -> List[str]:
    names = []
    with open(path, "rt", encoding="utf-8") as handle:
        for line_number, raw in enumerate(handle, 1):
            if not raw.strip() or raw.lstrip().startswith("#"):
                continue
            name = partition_name(raw)
            if not name or name.startswith(".") or "/" in name:
                raise ValueError(
                    f"{path}:{line_number}: invalid partition name {name!r}"
                )
            names.append(name)
    if not names:
        raise ValueError(f"partition list is empty: {path}")
    duplicates = [
        name for name, count in collections.Counter(names).items()
        if count > 1
    ]
    if duplicates:
        raise ValueError(f"duplicate partition in list: {duplicates[0]}")
    return names


def build_tasks(args: argparse.Namespace) -> List[BreakpointTask]:
    root = Path(args.graph_folder).expanduser().resolve()
    if not root.is_dir():
        raise NotADirectoryError(root)
    tasks = []
    for partition in read_partition_names(args.partition_list):
        directory = root / partition
        tasks.append(BreakpointTask(
            partition=partition,
            directory=str(directory),
            graph=str(directory / f"{partition}.FA"),
            bed=str(directory / f"{partition}.bed"),
            alignment=str(directory / f"{partition}_align.txt"),
            initial=str(
                directory / f"{partition}_breakpoint_consistency.initial.tsv"
            ),
            blocks=str(directory / f"{partition}_breakpoint_consistency.bed"),
            cache=str(directory / f"{partition}cache.json"),
            graphfixbreaks=os.path.abspath(args.graphfixbreaks),
            resume=args.resume,
        ))
    return tasks


def write_manifest(path: str, results: Iterable[BreakpointResult]) -> None:
    destination = Path(path).expanduser().resolve()
    destination.parent.mkdir(parents=True, exist_ok=True)
    temporary = destination.with_name(destination.name + f".tmp.{os.getpid()}")
    try:
        with temporary.open("wt", encoding="utf-8") as output:
            output.write(
                "partition\tstatus\talignment_rows\tsegments\tblocks\t"
                "cache\tblocking_output\tmessage\n"
            )
            for result in results:
                message = result.message.replace("\t", " ").replace("\n", " ")
                output.write("\t".join(map(str, (
                    result.partition, result.status, result.alignment_rows,
                    result.segment_count, result.block_count, result.cache,
                    result.blocks, message,
                ))) + "\n")
        os.replace(temporary, destination)
    finally:
        try:
            temporary.unlink()
        except FileNotFoundError:
            pass


def write_success(
    path: str, results: Sequence[BreakpointResult], only_reference: bool,
) -> None:
    counts = collections.Counter(result.status for result in results)
    payload = {
        "format": "minsetref-breakpoint-consistency-v1",
        "alignment_scope": "only_reference" if only_reference else "all_samples",
        "partitions": len(results),
        "status_counts": dict(sorted(counts.items())),
    }
    destination = Path(path).expanduser().resolve()
    destination.parent.mkdir(parents=True, exist_ok=True)
    temporary = destination.with_name(destination.name + f".tmp.{os.getpid()}")
    try:
        temporary.write_text(
            json.dumps(payload, indent=2, sort_keys=True) + "\n",
            encoding="utf-8",
        )
        os.replace(temporary, destination)
    finally:
        try:
            temporary.unlink()
        except FileNotFoundError:
            pass


def run(args: argparse.Namespace) -> int:
    tasks = build_tasks(args)
    worker_count = min(args.jobs, len(tasks))
    reference_haplotypes = tuple(
        args.reference_haplotype or ("CHM13_h1",)
    )
    LOG.info(
        "Building breakpoint consistency for %d partitions with %d spawned "
        "worker process%s (%s)",
        len(tasks), worker_count, "" if worker_count == 1 else "es",
        "reference only" if args.onlyref else "all samples",
    )
    results: List[BreakpointResult] = []
    if worker_count == 1:
        _init_worker(
            args.graphfixbreaks, reference_haplotypes,
            args.onlyref, args.merge_small,
        )
        for task in tasks:
            try:
                results.append(_process_task(task))
            except Exception as error:
                LOG.error("%s: %s", task.partition, error)
                results.append(BreakpointResult(
                    task.partition, "failed", message=str(error),
                ))
    else:
        context = multiprocessing.get_context("spawn")
        with concurrent.futures.ProcessPoolExecutor(
            max_workers=worker_count,
            mp_context=context,
            initializer=_init_worker,
            initargs=(
                args.graphfixbreaks, reference_haplotypes,
                args.onlyref, args.merge_small,
            ),
        ) as executor:
            task_iterator = iter(tasks)
            pending = {}
            for _index in range(min(len(tasks), worker_count * 2)):
                task = next(task_iterator, None)
                if task is None:
                    break
                pending[executor.submit(_process_task, task)] = task
            while pending:
                done, _waiting = concurrent.futures.wait(
                    pending, return_when=concurrent.futures.FIRST_COMPLETED,
                )
                for future in done:
                    task = pending.pop(future)
                    try:
                        results.append(future.result())
                    except Exception as error:
                        LOG.error("%s: %s", task.partition, error)
                        results.append(BreakpointResult(
                            task.partition, "failed", message=str(error),
                        ))
                    next_task = next(task_iterator, None)
                    if next_task is not None:
                        pending[executor.submit(
                            _process_task, next_task,
                        )] = next_task
                completed = len(results)
                if completed % max(1, len(tasks) // 100) == 0 or completed == len(tasks):
                    LOG.info(
                        "Breakpoint consistency progress: %d/%d partitions",
                        completed, len(tasks),
                    )

    order = {task.partition: index for index, task in enumerate(tasks)}
    results.sort(key=lambda result: order[result.partition])
    write_manifest(args.manifest, results)
    counts = collections.Counter(result.status for result in results)
    completed_count = counts["built"] + counts["reused"]
    LOG.info(
        "Breakpoint consistency complete: %d built, %d reused, %d skipped, "
        "%d failed",
        counts["built"], counts["reused"],
        len(results) - completed_count - counts["failed"], counts["failed"],
    )
    skipped_counts = {
        status: count
        for status, count in sorted(counts.items())
        if status.startswith("skipped_") and count
    }
    if skipped_counts:
        LOG.warning(
            "Breakpoint consistency skip reasons: %s",
            ", ".join(
                f"{status}={count}"
                for status, count in skipped_counts.items()
            ),
        )
    if completed_count == 0:
        LOG.error("No partition produced a current breakpoint cache")
        return 1
    if args.strict and counts["failed"]:
        return 1
    write_success(args.success, results, args.onlyref)
    return 0


def parse_args(argv: Optional[Sequence[str]] = None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--partition-list", required=True)
    parser.add_argument("--graph-folder", required=True)
    parser.add_argument("--manifest", required=True)
    parser.add_argument("--success", required=True)
    parser.add_argument(
        "--graphfixbreaks",
        default=str(Path(__file__).resolve().with_name("gfixbreaks.py")),
    )
    parser.add_argument("--reference-haplotype", action="append")
    parser.add_argument(
        "--onlyref", "--only-ref", action="store_true",
        help="use only named reference rows/intervals; default: all samples",
    )
    parser.add_argument("--merge-small", type=int, default=20_000)
    parser.add_argument("-j", "--jobs", type=int, default=16)
    parser.add_argument("--resume", action="store_true")
    parser.add_argument(
        "--strict", action="store_true",
        help="fail if any otherwise eligible partition raises an error",
    )
    parser.add_argument("-v", "--verbose", action="store_true")
    args = parser.parse_args(argv)
    if args.jobs < 1:
        parser.error("--jobs must be positive")
    if args.merge_small < 0:
        parser.error("--merge-small must be nonnegative")
    return args


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
            LOG.exception("breakpoint consistency failed")
        return 1


if __name__ == "__main__":
    raise SystemExit(main())
