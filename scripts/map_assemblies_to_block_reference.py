#!/usr/bin/env python3
"""Map assemblies to a blocked reference plus novel loci with Minimap2 cycles.

Each assembly is aligned with Minimap2 only.  Accepted query coverage is
removed, the exposed residual sequence is aligned again, and cycling continues
until a pass finds no additional alignment (or no residual remains).

All paired alignment blocks from all cycles are then combined by
query-contig/target/strand.  Neighboring query intervals are joined when their
gap is strictly less than ``--merge-gap``.  A group is reported when the union
of target-side uppercase A/C/G/T bases it covers is at least
``min(1000, 0.5 * total target uppercase A/C/G/T)`` by default.
"""
from __future__ import annotations

import argparse
import dataclasses
import hashlib
import heapq
import json
import logging
import os
import shutil
import subprocess
from concurrent.futures import ThreadPoolExecutor, as_completed
from pathlib import Path
from typing import Dict, List, Optional, Sequence, Tuple

from minsetref import AssemblyEntry, read_assembly_list
from minsetref_align import (
    build_minimap2_index,
    parse_winnowmap_paf,
    run_minimap2,
)
from minsetref_core import (
    IndexedFasta,
    iter_n_free_segments,
    open_maybe_gzip,
    sanitize_id,
    wrap_fasta,
)
from minsetref_segments import complement_intervals, count_unmasked, merge_intervals


LOG = logging.getLogger("map_assemblies_to_block_reference")
MINIMAP_PARAMS = {
    "preset": "asm5",
    "N": 1000,
    "p": 0.001,
    "f": 0.001,
    "K": "100M",
    "retry_on_sigkill": True,
    "retry_threads": 32,
    "retry_K": "100M",
}
DEFAULT_THREADS = max(1, min(50, os.cpu_count() or 1))
MIN_IDENTITY = 95.0
MIN_ALIGNMENT_BP = 1
MIN_RESIDUAL_UNMASKED = 1
MERGE_GAP = 500
TARGET_UNMASKED_CAP = 1000.0
MAX_CYCLES = 0
BED_COLUMNS = (
    "query_contig",
    "query_start",
    "query_end",
    "target",
    "score",
    "strand",
    "sample",
    "target_start",
    "target_end",
    "target_unmasked_covered",
    "target_unmasked_total",
    "required_unmasked_exclusive",
    "mean_identity",
    "alignment_blocks",
    "first_cycle",
    "last_cycle",
)


@dataclasses.dataclass(frozen=True)
class QuerySource:
    contig: str
    start: int
    end: int


@dataclasses.dataclass
class AlignmentGroup:
    target: str
    strand: str
    query_contig: str
    query_start: int
    query_end: int
    target_intervals: List[Tuple[int, int]]
    identity_numerator: float
    identity_denominator: int
    block_count: int
    first_cycle: int
    last_cycle: int

    def add(
        self,
        query_start: int,
        query_end: int,
        target_start: int,
        target_end: int,
        identity: float,
        cycle: int,
    ) -> None:
        self.query_start = min(self.query_start, query_start)
        self.query_end = max(self.query_end, query_end)
        self.target_intervals.append((target_start, target_end))
        width = max(0, query_end - query_start)
        self.identity_numerator += identity * width
        self.identity_denominator += width
        self.block_count += 1
        self.first_cycle = min(self.first_cycle, cycle)
        self.last_cycle = max(self.last_cycle, cycle)


@dataclasses.dataclass
class BedMapping:
    fields: List[str]
    row_id: int
    contig: str
    start: int
    end: int
    target: str
    score: int
    strand: str
    target_highest: int = 0

    @property
    def target_key(self) -> Tuple[str, str]:
        return self.target, self.strand


@dataclasses.dataclass
class FilteredMappingRun:
    contig: str
    start: int
    end: int
    target_key: Tuple[str, str]
    best: BedMapping
    effective_score: int


def default_filtered_bed_path(path: str) -> str:
    return path[:-4] + ".filtered.bed" if path.lower().endswith(".bed") else path + ".filtered.bed"


def _filtered_run_row(run: FilteredMappingRun) -> str:
    fields = list(run.best.fields)
    fields[0] = run.contig
    fields[1] = str(run.start)
    fields[2] = str(run.end)
    fields[4] = str(run.effective_score)
    return "\t".join(fields) + "\n"


def merge_short_mapping_runs(
    runs: Sequence[FilteredMappingRun],
    minimum_block_size: int,
) -> List[FilteredMappingRun]:
    """Recursively absorb short query runs into an immediate neighbor.

    Only the absorbing run's query coordinates are enlarged; its target,
    strand, score, and all other BED metadata remain unchanged.  A neighbor is
    eligible only when its query-coordinate gap from the short run is no more
    than half the smaller of their two sizes. Among eligible neighbors,
    selection prefers the same target/strand, then the higher effective score,
    then the longer run, and finally the left neighbor. Processing is repeated
    from the leftmost mergeable short run until no further merge is possible.

    A lone run shorter than ``minimum_block_size`` is retained because it has
    no neighbor into which it can be merged.
    """
    if minimum_block_size < 0:
        raise ValueError("minimum_block_size must be nonnegative")
    merged = list(runs)
    merged.sort(key=lambda run: (
        run.start, run.end, run.target_key[0], run.target_key[1],
        -run.effective_score,
    ))
    if minimum_block_size == 0:
        return merged

    while len(merged) > 1:
        selected_merge = None
        for short_index, short in enumerate(merged):
            short_size = short.end - short.start
            if short_size >= minimum_block_size:
                continue
            neighbor_indexes = []
            if short_index > 0:
                neighbor_indexes.append(short_index - 1)
            if short_index + 1 < len(merged):
                neighbor_indexes.append(short_index + 1)
            eligible_indexes = []
            for index in neighbor_indexes:
                neighbor = merged[index]
                neighbor_size = neighbor.end - neighbor.start
                if index < short_index:
                    gap = max(0, short.start - neighbor.end)
                else:
                    gap = max(0, neighbor.start - short.end)
                if 2 * gap <= min(neighbor_size, short_size):
                    eligible_indexes.append(index)
            if eligible_indexes:
                selected_merge = (short_index, short, eligible_indexes)
                break
        if selected_merge is None:
            break
        short_index, short, neighbor_indexes = selected_merge

        def neighbor_priority(index: int) -> Tuple[int, int, int, int]:
            neighbor = merged[index]
            return (
                int(neighbor.target_key == short.target_key),
                neighbor.effective_score,
                neighbor.end - neighbor.start,
                int(index < short_index),
            )

        neighbor_index = max(neighbor_indexes, key=neighbor_priority)
        neighbor = merged[neighbor_index]
        neighbor.start = min(neighbor.start, short.start)
        neighbor.end = max(neighbor.end, short.end)
        del merged[short_index]
        merged.sort(key=lambda run: (
            run.start, run.end, run.target_key[0], run.target_key[1],
            -run.effective_score,
        ))

    return merged


def _filter_one_contig(
    mappings: Sequence[BedMapping],
    relative_numerator: int,
    relative_denominator: int,
    minimum_score_exclusive: Optional[int],
    scale_partial_scores: bool,
) -> List[FilteredMappingRun]:
    """Sweep endpoints and retain all target/strand paths near the local best.

    Per-target heaps and integer score buckets make updates logarithmic in the
    number of overlapping records. No alignment pair is compared directly.
    """
    if not mappings:
        return []
    active = [False] * len(mappings)
    target_heaps: Dict[Tuple[str, str], List[Tuple[int, int]]] = {}
    best_by_target: Dict[Tuple[str, str], BedMapping] = {}
    score_buckets: List[set] = [set() for _ in range(1001)]
    events: List[Tuple[int, int, int]] = []
    for local_id, mapping in enumerate(mappings):
        if mapping.end <= mapping.start:
            continue
        if not 0 <= mapping.score <= 1000:
            raise ValueError(
                f"BED score must be in [0,1000], got {mapping.score} for "
                f"{mapping.contig}:{mapping.start}-{mapping.end}"
            )
        # End events sort before starts at the same half-open coordinate.
        events.append((mapping.start, 1, local_id))
        events.append((mapping.end, 0, local_id))
    events.sort()
    if not events:
        return []

    open_runs: Dict[Tuple[str, str], FilteredMappingRun] = {}
    finished: List[FilteredMappingRun] = []

    def refresh_target(key: Tuple[str, str]) -> None:
        previous = best_by_target.pop(key, None)
        if previous is not None:
            score_buckets[previous.score].discard(key)
        heap = target_heaps.get(key, [])
        while heap and not active[heap[0][1]]:
            heapq.heappop(heap)
        if heap:
            current = mappings[heap[0][1]]
            best_by_target[key] = current
            score_buckets[current.score].add(key)

    event_index = 0
    while event_index < len(events):
        coordinate = events[event_index][0]
        changed_keys = set()
        while event_index < len(events) and events[event_index][0] == coordinate:
            _position, event_type, local_id = events[event_index]
            mapping = mappings[local_id]
            key = mapping.target_key
            if event_type == 0:
                active[local_id] = False
            else:
                active[local_id] = True
                heapq.heappush(
                    target_heaps.setdefault(key, []),
                    (-mapping.score, local_id),
                )
            changed_keys.add(key)
            event_index += 1
        for key in changed_keys:
            refresh_target(key)

        next_coordinate = (
            events[event_index][0] if event_index < len(events) else coordinate
        )
        selected: Dict[Tuple[str, str], BedMapping] = {}
        if next_coordinate > coordinate:
            highest = next(
                (score for score in range(1000, -1, -1)
                 if score_buckets[score]),
                None,
            )
            if highest is not None:
                for score in range(highest, -1, -1):
                    for key in sorted(score_buckets[score]):
                        mapping = best_by_target[key]
                        # Mutual filtering: remove only when the mapping is
                        # weak relative to both its query-side and target-side
                        # competitors. Equality to 30% is retained.
                        query_weak = (
                            mapping.score * relative_denominator
                            < highest * relative_numerator
                        )
                        target_weak = (
                            mapping.score * relative_denominator
                            < mapping.target_highest * relative_numerator
                        )
                        if not (query_weak and target_weak):
                            selected[key] = mapping

        for key in list(open_runs):
            if key not in selected:
                finished.append(open_runs.pop(key))
        for key, mapping in selected.items():
            run = open_runs.get(key)
            # Do not merge adjacent segments owned by different input rows.
            # Remaining-size score scaling is defined per original mapping.
            mapping_changed = (
                scale_partial_scores
                and run is not None
                and run.best.row_id != mapping.row_id
            )
            if run is None or run.end != coordinate or mapping_changed:
                if run is not None:
                    finished.append(run)
                open_runs[key] = FilteredMappingRun(
                    mapping.contig, coordinate, next_coordinate, key, mapping,
                    mapping.score,
                )
                continue
            run.end = next_coordinate
            if not scale_partial_scores and mapping.score > run.best.score:
                run.best = mapping
                run.effective_score = mapping.score

    finished.extend(open_runs.values())
    # Multi-mapping competition may retain only part of an original interval.
    # Scale its score once by the union of all retained query fragments:
    #     effective = original_score * retained_bp / original_bp
    # Round to the integer BED score first, then reapply the strict cutoff so
    # every emitted row still visibly satisfies column 5 > cutoff and a second
    # filtering pass is idempotent.
    if scale_partial_scores:
        if minimum_score_exclusive is None:
            raise ValueError(
                "scale_partial_scores requires minimum_score_exclusive"
            )
        runs_by_row: Dict[int, List[FilteredMappingRun]] = {}
        for run in finished:
            runs_by_row.setdefault(run.best.row_id, []).append(run)
        scaled: List[FilteredMappingRun] = []
        for runs in runs_by_row.values():
            mapping = runs[0].best
            total_size = mapping.end - mapping.start
            if total_size <= 0:
                continue
            retained_size = sum(
                end - start
                for start, end in merge_intervals([
                    (run.start, run.end) for run in runs
                ])
            )
            score_numerator = mapping.score * retained_size
            effective_score = int(round(score_numerator / total_size))
            if effective_score <= minimum_score_exclusive:
                continue
            for run in runs:
                run.effective_score = effective_score
                scaled.append(run)
        finished = scaled
    finished.sort(key=lambda run: (
        run.start, run.end, run.target_key[0], run.target_key[1], -run.best.score,
    ))
    return finished


def _assign_target_highest(mappings: Sequence[BedMapping]) -> None:
    """Assign each row the best score overlapping its target interval.

    Target records are independent. Within each target, an endpoint sweep
    builds the local-best score on each elementary interval, followed by an
    iterative range-maximum index. Complexity is O(N log N), not all-pairs.
    Both target strands compete because they cover the same target coordinates.
    """
    by_target: Dict[str, List[BedMapping]] = {}
    for mapping in mappings:
        by_target.setdefault(mapping.target, []).append(mapping)

    for target_rows in by_target.values():
        coordinates = sorted({
            coordinate
            for mapping in target_rows
            for coordinate in (
                int(mapping.fields[7]), int(mapping.fields[8]),
            )
        })
        if len(coordinates) < 2:
            for mapping in target_rows:
                mapping.target_highest = mapping.score
            continue
        coordinate_index = {coordinate: index for index, coordinate in enumerate(coordinates)}
        events: List[Tuple[int, int, int]] = []
        active = [False] * len(target_rows)
        for row_id, mapping in enumerate(target_rows):
            target_start = int(mapping.fields[7])
            target_end = int(mapping.fields[8])
            if target_end <= target_start:
                mapping.target_highest = mapping.score
                continue
            events.append((target_start, 1, row_id))
            events.append((target_end, 0, row_id))
        events.sort()
        heap: List[Tuple[int, int]] = []
        best_segments = [0] * (len(coordinates) - 1)
        event_index = 0
        for segment_index, coordinate in enumerate(coordinates[:-1]):
            while event_index < len(events) and events[event_index][0] == coordinate:
                _position, event_type, row_id = events[event_index]
                if event_type == 0:
                    active[row_id] = False
                else:
                    active[row_id] = True
                    heapq.heappush(heap, (-target_rows[row_id].score, row_id))
                event_index += 1
            while heap and not active[heap[0][1]]:
                heapq.heappop(heap)
            if heap:
                best_segments[segment_index] = -heap[0][0]

        tree_size = 1
        while tree_size < len(best_segments):
            tree_size *= 2
        range_tree = [0] * (2 * tree_size)
        range_tree[tree_size:tree_size + len(best_segments)] = best_segments
        for index in range(tree_size - 1, 0, -1):
            range_tree[index] = max(range_tree[2 * index], range_tree[2 * index + 1])

        def range_max(left: int, right: int) -> int:
            left += tree_size
            right += tree_size
            answer = 0
            while left < right:
                if left & 1:
                    answer = max(answer, range_tree[left])
                    left += 1
                if right & 1:
                    right -= 1
                    answer = max(answer, range_tree[right])
                left //= 2
                right //= 2
            return answer

        for mapping in target_rows:
            target_start = int(mapping.fields[7])
            target_end = int(mapping.fields[8])
            if target_end > target_start:
                mapping.target_highest = range_max(
                    coordinate_index[target_start], coordinate_index[target_end],
                )
            if mapping.target_highest == 0:
                mapping.target_highest = mapping.score


def filter_competing_mappings_bed(
    input_bed: str,
    output_bed: str,
    relative_numerator: int = 3,
    relative_denominator: int = 10,
    target_unmasked_cap: float = 1000.0,
    minimum_score_exclusive: Optional[int] = None,
    scale_partial_scores: bool = False,
    minimum_block_size: int = 0,
) -> int:
    """Filter a query-sorted mapper BED with an endpoint sweep.

    A row first needs target_unmasked_covered >= min(target_unmasked_cap,
    0.5 * target_unmasked_total). When ``minimum_score_exclusive`` is given,
    column 5 must also exceed it. ``scale_partial_scores`` enables the separate
    downstream filter's retained-size score update. A row is otherwise removed
    only when weak relative to both query and target local best scores.  When
    ``minimum_block_size`` is positive, surviving shorter query runs are
    recursively absorbed into an immediate neighbor.
    """
    if minimum_block_size < 0:
        raise ValueError("minimum_block_size must be nonnegative")
    if os.path.realpath(input_bed) == os.path.realpath(output_bed):
        raise ValueError("Filtered BED output must differ from its input")
    Path(output_bed).parent.mkdir(parents=True, exist_ok=True)
    temporary = output_bed + ".tmp"
    rows_written = 0
    mappings: List[BedMapping] = []
    current_contig: Optional[str] = None
    previous_start = -1

    try:
        with open(input_bed) as source:
            row_id = 0
            for line_number, raw in enumerate(source, 1):
                if not raw.strip() or raw.startswith("#"):
                    continue
                fields = raw.rstrip("\n").split("\t")
                if len(fields) < 11:
                    raise ValueError(
                        f"{input_bed}:{line_number}: expected the mapper's 16-column BED"
                    )
                try:
                    start, end, score = int(fields[1]), int(fields[2]), int(fields[4])
                    target_start = int(fields[7])
                    target_end = int(fields[8])
                    target_covered = int(fields[9])
                    target_total = int(fields[10])
                except ValueError as error:
                    raise ValueError(
                        f"{input_bed}:{line_number}: invalid numeric BED field"
                    ) from error
                if target_end < target_start:
                    raise ValueError(
                        f"{input_bed}:{line_number}: target_end is before target_start"
                    )
                contig = fields[0]
                if current_contig is None:
                    current_contig = contig
                    previous_start = -1
                elif contig != current_contig:
                    current_contig = contig
                    previous_start = -1
                if start < previous_start:
                    raise ValueError(
                        f"{input_bed} is not query-coordinate sorted at line {line_number}"
                    )
                previous_start = start
                # Column 5 is target unmasked coverage on a 0..1000 scale.
                # The requested cutoff is strict: 200 is removed; 201 passes.
                if (
                    minimum_score_exclusive is not None
                    and score <= minimum_score_exclusive
                ):
                    continue
                required_target_coverage = min(
                    float(target_unmasked_cap), 0.5 * target_total,
                )
                if target_covered < required_target_coverage:
                    continue
                mappings.append(BedMapping(
                    fields, row_id, contig, start, end,
                    fields[3], score, fields[5],
                ))
                row_id += 1

        _assign_target_highest(mappings)
        by_query: Dict[str, List[BedMapping]] = {}
        query_order: List[str] = []
        for mapping in mappings:
            if mapping.contig not in by_query:
                by_query[mapping.contig] = []
                query_order.append(mapping.contig)
            by_query[mapping.contig].append(mapping)
        with open(temporary, "w") as output:
            for contig in query_order:
                rows = by_query[contig]
                rows.sort(key=lambda row: (row.start, row.end, row.target, row.strand))
                runs = _filter_one_contig(
                    rows, relative_numerator, relative_denominator,
                    minimum_score_exclusive, scale_partial_scores,
                )
                runs = merge_short_mapping_runs(runs, minimum_block_size)
                for run in runs:
                    output.write(_filtered_run_row(run))
                    rows_written += 1
        os.replace(temporary, output_bed)
    finally:
        try:
            os.remove(temporary)
        except FileNotFoundError:
            pass
    return rows_written


def atomic_json(path: str, payload: dict) -> None:
    Path(path).parent.mkdir(parents=True, exist_ok=True)
    temporary = path + ".tmp"
    with open(temporary, "w") as out:
        json.dump(payload, out, sort_keys=True)
        out.write("\n")
        out.flush()
        os.fsync(out.fileno())
    os.replace(temporary, path)


def read_json(path: str) -> Optional[dict]:
    if not os.path.isfile(path):
        return None
    try:
        with open(path) as handle:
            value = json.load(handle)
    except (OSError, TypeError, ValueError):
        return None
    return value if isinstance(value, dict) else None


def file_fingerprint(path: str) -> dict:
    stat = os.stat(path)
    return {
        "path": os.path.realpath(path),
        "size": stat.st_size,
        "mtime_ns": stat.st_mtime_ns,
    }


def hash_payload(payload: object) -> str:
    encoded = json.dumps(
        payload, sort_keys=True, separators=(",", ":"),
    ).encode()
    return hashlib.sha256(encoded).hexdigest()


def setup_logging(workdir: str, verbose: bool) -> None:
    Path(workdir).mkdir(parents=True, exist_ok=True)
    LOG.setLevel(logging.DEBUG)
    LOG.handlers.clear()
    formatter = logging.Formatter(
        "[%(asctime)s] %(levelname)s: %(message)s",
        datefmt="%Y-%m-%d %H:%M:%S",
    )
    stream = logging.StreamHandler()
    stream.setLevel(logging.DEBUG if verbose else logging.INFO)
    stream.setFormatter(formatter)
    logfile = logging.FileHandler(
        os.path.join(workdir, "mapping.log"), mode="a",
    )
    logfile.setLevel(logging.DEBUG)
    logfile.setFormatter(formatter)
    LOG.addHandler(stream)
    LOG.addHandler(logfile)


def build_target_bundle(
    input_fastas: Sequence[str],
    target_fasta: str,
    target_table: str,
    marker_path: str,
    resume: bool,
) -> Dict[str, Tuple[int, int]]:
    """Build one duplicate-checked target FASTA and return length/unmasked."""
    signature = hash_payload([file_fingerprint(path) for path in input_fastas])
    marker = read_json(marker_path) if resume else None
    if (
        marker is not None
        and marker.get("signature") == signature
        and os.path.isfile(target_fasta)
        and os.path.getsize(target_fasta) > 0
        and os.path.isfile(target_table)
    ):
        stats: Dict[str, Tuple[int, int]] = {}
        with open(target_table) as table:
            header = next(table, "").rstrip("\n").split("\t")
            if header != ["target", "length", "unmasked"]:
                raise RuntimeError(f"Invalid target table: {target_table}")
            for raw in table:
                target, length, unmasked = raw.rstrip("\n").split("\t")
                stats[target] = (int(length), int(unmasked))
        if len(stats) != int(marker.get("record_count", -1)):
            raise RuntimeError(f"Target checkpoint count mismatch: {target_table}")
        LOG.info("RESUME: reusing combined target with %d records", len(stats))
        return stats

    Path(target_fasta).parent.mkdir(parents=True, exist_ok=True)
    fasta_tmp = target_fasta + ".tmp"
    table_tmp = target_table + ".tmp"
    stats: Dict[str, Tuple[int, int]] = {}
    current_name: Optional[str] = None
    current_length = 0
    current_unmasked = 0

    def finish_record() -> None:
        nonlocal current_name, current_length, current_unmasked
        if current_name is None:
            return
        stats[current_name] = (current_length, current_unmasked)
        current_name = None
        current_length = 0
        current_unmasked = 0

    try:
        with open(fasta_tmp, "w") as fasta_out:
            for source_path in input_fastas:
                LOG.info("Adding target FASTA %s", source_path)
                with open_maybe_gzip(source_path, "rt") as source:
                    for line_number, raw in enumerate(source, 1):
                        if raw.startswith(">"):
                            finish_record()
                            fields = raw[1:].strip().split()
                            if not fields:
                                raise ValueError(
                                    f"Empty FASTA header at {source_path}:{line_number}"
                                )
                            current_name = fields[0]
                            if current_name in stats:
                                raise ValueError(
                                    f"Duplicate target FASTA id {current_name!r} "
                                    f"in {source_path}:{line_number}"
                                )
                            fasta_out.write(raw.rstrip("\r\n") + "\n")
                            continue
                        sequence = raw.strip()
                        if not sequence:
                            continue
                        if current_name is None:
                            raise ValueError(
                                f"Sequence before FASTA header at "
                                f"{source_path}:{line_number}"
                            )
                        current_length += len(sequence)
                        current_unmasked += count_unmasked(sequence)
                        fasta_out.write(sequence + "\n")
                finish_record()
        if not stats:
            raise ValueError("Target FASTAs contain no records")
        with open(table_tmp, "w") as table:
            table.write("target\tlength\tunmasked\n")
            for target, (length, unmasked) in stats.items():
                table.write(f"{target}\t{length}\t{unmasked}\n")
        os.replace(fasta_tmp, target_fasta)
        os.replace(table_tmp, target_table)
        atomic_json(marker_path, {
            "kind": "combined_target",
            "signature": signature,
            "record_count": len(stats),
            "total_bp": sum(length for length, _unmasked in stats.values()),
        })
    finally:
        for temporary in (fasta_tmp, table_tmp):
            try:
                os.remove(temporary)
            except FileNotFoundError:
                pass
    LOG.info(
        "Built combined target: %d records, %.3f Gb",
        len(stats),
        sum(length for length, _unmasked in stats.values()) / 1e9,
    )
    return stats


def scan_target_stats(target_fasta: str) -> Dict[str, Tuple[int, int]]:
    """Count target length/unmasked bp without copying the target FASTA."""
    if target_fasta.endswith(".gz"):
        raise ValueError(
            f"Target FASTA must be uncompressed for base counting: {target_fasta}"
        )
    stats: Dict[str, Tuple[int, int]] = {}
    current_name: Optional[str] = None
    current_length = 0
    current_unmasked = 0

    def finish_record() -> None:
        nonlocal current_name, current_length, current_unmasked
        if current_name is None:
            return
        stats[current_name] = (current_length, current_unmasked)
        current_name = None
        current_length = 0
        current_unmasked = 0

    with open(target_fasta) as source:
        for line_number, raw in enumerate(source, 1):
            if raw.startswith(">"):
                finish_record()
                fields = raw[1:].strip().split()
                if not fields:
                    raise ValueError(
                        f"Empty FASTA header at {target_fasta}:{line_number}"
                    )
                current_name = fields[0]
                if current_name in stats:
                    raise ValueError(
                        f"Duplicate target FASTA id {current_name!r} at "
                        f"{target_fasta}:{line_number}"
                    )
                continue
            sequence = raw.strip()
            if not sequence:
                continue
            if current_name is None:
                raise ValueError(
                    f"Sequence before FASTA header at "
                    f"{target_fasta}:{line_number}"
                )
            current_length += len(sequence)
            current_unmasked += count_unmasked(sequence)
    finish_record()
    if not stats:
        raise ValueError(f"Target FASTA contains no records: {target_fasta}")
    return stats


def external_sort(
    input_path: str,
    output_path: str,
    keys: Sequence[str],
    temporary_dir: str,
) -> None:
    """Use the system external-memory sort with a deterministic locale."""
    Path(output_path).parent.mkdir(parents=True, exist_ok=True)
    if not os.path.isfile(input_path) or os.path.getsize(input_path) == 0:
        Path(output_path).touch()
        return
    command = ["sort", "-t", "\t", "-T", temporary_dir, *keys, input_path]
    environment = os.environ.copy()
    environment["LC_ALL"] = "C"
    temporary = output_path + ".tmp"
    try:
        with open(temporary, "w") as out:
            subprocess.run(
                command, check=True, stdout=out, env=environment,
            )
        os.replace(temporary, output_path)
    finally:
        try:
            os.remove(temporary)
        except FileNotFoundError:
            pass


def initial_query_map(fasta_path: str) -> Dict[str, QuerySource]:
    if fasta_path.endswith(".gz"):
        raise ValueError(
            f"Assembly FASTA must be uncompressed for residual cycling: {fasta_path}"
        )
    with IndexedFasta(fasta_path) as reader:
        return {
            name: QuerySource(name, 0, reader.length(name))
            for name in reader.names()
        }


def parse_cycle_alignments(
    paf_path: str,
    query_map: Dict[str, QuerySource],
    raw_blocks_path: str,
    coverage_path: str,
    min_identity: float,
    min_alignment_bp: int,
    cycle: int,
) -> Tuple[int, int]:
    """Append original-coordinate blocks and write current-query coverage."""
    hit_count = 0
    block_count = 0
    with open(raw_blocks_path, "a") as block_out, open(coverage_path, "w") as coverage:
        for hit in parse_winnowmap_paf(
            paf_path, min_identity, min_alignment_bp,
        ):
            source = query_map.get(hit.query_id)
            if source is None:
                raise RuntimeError(
                    f"Minimap2 returned unknown query {hit.query_id!r} in {paf_path}"
                )
            valid_pairs = [
                (int(q0), int(q1), int(t0), int(t1))
                for q0, q1, t0, t1 in hit.aligned_pairs
                if q1 > q0 and t1 > t0
            ]
            if not valid_pairs:
                continue
            hit_count += 1
            for q0, q1 in hit.query_intervals:
                if q1 > q0:
                    coverage.write(f"{hit.query_id}\t{q0}\t{q1}\n")
            for q0, q1, t0, t1 in valid_pairs:
                block_count += 1
                block_out.write("\t".join([
                    hit.target_id,
                    hit.strand,
                    source.contig,
                    str(source.start + q0),
                    str(source.start + q1),
                    str(t0),
                    str(t1),
                    repr(hit.identity),
                    repr(hit.alignment_score),
                    str(cycle),
                ]) + "\n")
    return hit_count, block_count


def load_merged_coverage(
    coverage_path: str,
    sorted_path: str,
    temporary_dir: str,
) -> Dict[str, List[Tuple[int, int]]]:
    external_sort(
        coverage_path,
        sorted_path,
        ("-k1,1", "-k2,2n", "-k3,3n"),
        temporary_dir,
    )
    coverage: Dict[str, List[Tuple[int, int]]] = {}
    current_query: Optional[str] = None
    current_intervals: List[Tuple[int, int]] = []

    def finish() -> None:
        if current_query is not None:
            coverage[current_query] = merge_intervals(current_intervals)

    with open(sorted_path) as source:
        for raw in source:
            query, start, end = raw.rstrip("\n").split("\t")
            if query != current_query:
                finish()
                current_query = query
                current_intervals = []
            current_intervals.append((int(start), int(end)))
    finish()
    return coverage


def write_residual_queries(
    query_fasta: str,
    query_map: Dict[str, QuerySource],
    coverage: Dict[str, List[Tuple[int, int]]],
    output_fasta: str,
    min_unmasked: int,
    next_cycle: int,
) -> Tuple[Dict[str, QuerySource], int, int]:
    """Write only newly exposed gaps; no-hit records are retired."""
    next_map: Dict[str, QuerySource] = {}
    residual_bp = 0
    serial = 0
    Path(output_fasta).parent.mkdir(parents=True, exist_ok=True)
    with IndexedFasta(query_fasta) as reader, open(output_fasta, "w") as out:
        for query_id in sorted(coverage):
            source = query_map.get(query_id)
            if source is None:
                raise RuntimeError(f"Coverage has unknown query {query_id!r}")
            query_length = reader.length(query_id)
            for gap_start, gap_end in complement_intervals(
                coverage[query_id], 0, query_length,
            ):
                gap_sequence = reader.fetch(query_id, gap_start, gap_end)
                for local_start, local_end, sequence in iter_n_free_segments(
                    gap_sequence,
                ):
                    if count_unmasked(sequence) < min_unmasked:
                        continue
                    serial += 1
                    source_start = source.start + gap_start + local_start
                    source_end = source.start + gap_start + local_end
                    residual_id = f"res_c{next_cycle:03d}_{serial:09d}"
                    next_map[residual_id] = QuerySource(
                        source.contig, source_start, source_end,
                    )
                    residual_bp += len(sequence)
                    out.write(
                        f">{residual_id} source={source.contig}:"
                        f"{source_start}-{source_end}\n{wrap_fasta(sequence)}\n"
                    )
    return next_map, serial, residual_bp


def write_query_map(path: str, query_map: Dict[str, QuerySource]) -> None:
    with open(path, "w") as out:
        out.write("query_id\tsource_contig\tsource_start\tsource_end\n")
        for query_id, source in query_map.items():
            out.write(
                f"{query_id}\t{source.contig}\t{source.start}\t{source.end}\n"
            )


def new_group(
    target: str,
    strand: str,
    query_contig: str,
    query_start: int,
    query_end: int,
    target_start: int,
    target_end: int,
    identity: float,
    cycle: int,
) -> AlignmentGroup:
    width = query_end - query_start
    return AlignmentGroup(
        target=target,
        strand=strand,
        query_contig=query_contig,
        query_start=query_start,
        query_end=query_end,
        target_intervals=[(target_start, target_end)],
        identity_numerator=identity * width,
        identity_denominator=width,
        block_count=1,
        first_cycle=cycle,
        last_cycle=cycle,
    )


def format_group(
    group: AlignmentGroup,
    sample: str,
    target_sequence: str,
    target_unmasked_total: int,
    threshold_cap: float,
) -> Optional[str]:
    target_intervals = merge_intervals(group.target_intervals)
    covered_unmasked = sum(
        count_unmasked(target_sequence[start:end])
        for start, end in target_intervals
    )
    required = min(float(threshold_cap), 0.5 * target_unmasked_total)
    if covered_unmasked < required:
        return None
    target_start = min(start for start, _end in target_intervals)
    target_end = max(end for _start, end in target_intervals)
    coverage_fraction = (
        covered_unmasked / target_unmasked_total
        if target_unmasked_total else 0.0
    )
    score = min(1000, int(round(1000.0 * coverage_fraction)))
    mean_identity = (
        group.identity_numerator / group.identity_denominator
        if group.identity_denominator else 0.0
    )
    return "\t".join([
        group.query_contig,
        str(group.query_start),
        str(group.query_end),
        group.target,
        str(score),
        group.strand,
        sample,
        str(target_start),
        str(target_end),
        str(covered_unmasked),
        str(target_unmasked_total),
        f"{required:.1f}",
        f"{mean_identity:.6f}",
        str(group.block_count),
        str(group.first_cycle),
        str(group.last_cycle),
    ]) + "\n"


def group_blocks_to_bed(
    raw_blocks_path: str,
    sorted_blocks_path: str,
    unsorted_bed_path: str,
    output_bed_path: str,
    target_fasta: str,
    target_stats: Dict[str, Tuple[int, int]],
    sample: str,
    merge_gap: int,
    threshold_cap: float,
    temporary_dir: str,
) -> int:
    """Sorted sweep over blocks; never compares every block with every group."""
    external_sort(
        raw_blocks_path,
        sorted_blocks_path,
        ("-k1,1", "-k2,2", "-k3,3", "-k4,4n", "-k5,5n"),
        temporary_dir,
    )
    kept = 0
    group: Optional[AlignmentGroup] = None
    loaded_target: Optional[str] = None
    target_sequence = ""

    with IndexedFasta(target_fasta) as target_reader, open(
        unsorted_bed_path, "w",
    ) as bed:

        def finish_group() -> None:
            nonlocal kept, group, loaded_target, target_sequence
            if group is None:
                return
            if loaded_target != group.target:
                if group.target not in target_stats:
                    raise RuntimeError(
                        f"Alignment names unknown target {group.target!r}"
                    )
                target_sequence = target_reader.sequence(group.target)
                loaded_target = group.target
            row = format_group(
                group,
                sample,
                target_sequence,
                target_stats[group.target][1],
                threshold_cap,
            )
            if row is not None:
                bed.write(row)
                kept += 1

        with open(sorted_blocks_path) as blocks:
            for raw in blocks:
                (
                    target,
                    strand,
                    query_contig,
                    query_start_text,
                    query_end_text,
                    target_start_text,
                    target_end_text,
                    identity_text,
                    _alignment_score,
                    cycle_text,
                ) = raw.rstrip("\n").split("\t")
                query_start = int(query_start_text)
                query_end = int(query_end_text)
                target_start = int(target_start_text)
                target_end = int(target_end_text)
                identity = float(identity_text)
                cycle = int(cycle_text)
                same_key = (
                    group is not None
                    and group.target == target
                    and group.strand == strand
                    and group.query_contig == query_contig
                )
                # Strictly less than merge_gap, as requested.
                joins = (
                    same_key
                    and query_start - group.query_end < merge_gap
                )
                if not joins:
                    finish_group()
                    group = new_group(
                        target, strand, query_contig,
                        query_start, query_end, target_start, target_end,
                        identity, cycle,
                    )
                else:
                    group.add(
                        query_start, query_end, target_start, target_end,
                        identity, cycle,
                    )
        finish_group()

    external_sort(
        unsorted_bed_path,
        output_bed_path,
        ("-k1,1", "-k2,2n", "-k3,3n", "-k4,4"),
        temporary_dir,
    )
    return kept


def compact_sample_work(sample_dir: str) -> None:
    for name in os.listdir(sample_dir):
        if name in {"aligned.bed", "_SUCCESS.json"}:
            continue
        path = os.path.join(sample_dir, name)
        if os.path.isdir(path):
            shutil.rmtree(path)
        else:
            os.remove(path)


def sample_marker(
    assembly: AssemblyEntry,
    config_signature: str,
) -> dict:
    return {
        "kind": "assembly_block_mapping",
        "sample": assembly.sample,
        "assembly": file_fingerprint(assembly.path),
        "config_signature": config_signature,
    }


def map_one_assembly(
    assembly: AssemblyEntry,
    target_index: str,
    target_fasta: str,
    target_stats: Dict[str, Tuple[int, int]],
    sample_dir: str,
    args,
    config_signature: str,
) -> dict:
    expected = sample_marker(assembly, config_signature)
    marker_path = os.path.join(sample_dir, "_SUCCESS.json")
    shard_bed = os.path.join(sample_dir, "aligned.bed")
    observed = read_json(marker_path) if args.resume else None
    if (
        observed is not None
        and all(observed.get(key) == value for key, value in expected.items())
        and os.path.isfile(shard_bed)
    ):
        LOG.info("CHECKPOINT: skipped finished assembly %s", assembly.sample)
        return observed

    if os.path.isdir(sample_dir):
        shutil.rmtree(sample_dir)
    Path(sample_dir).mkdir(parents=True, exist_ok=True)
    raw_blocks = os.path.join(sample_dir, "all_blocks.unsorted.tsv")
    Path(raw_blocks).touch()
    current_fasta = assembly.path
    current_map = initial_query_map(current_fasta)
    cycle = 0
    total_hits = 0
    total_blocks = 0
    cycles_run = 0

    while current_map:
        if args.max_cycles and cycle >= args.max_cycles:
            LOG.warning(
                "%s reached --max-cycles=%d before a zero-alignment pass",
                assembly.sample, args.max_cycles,
            )
            break
        cycle_dir = os.path.join(sample_dir, f"cycle{cycle:03d}")
        Path(cycle_dir).mkdir(parents=True, exist_ok=True)
        write_query_map(os.path.join(cycle_dir, "queries.map.tsv"), current_map)
        paf = os.path.join(cycle_dir, "align.minimap2.paf")
        run_minimap2(
            current_fasta,
            target_index,
            paf,
            args.threads,
            LOG,
            params=MINIMAP_PARAMS,
        )
        coverage_path = os.path.join(cycle_dir, "coverage.unsorted.tsv")
        hits, blocks = parse_cycle_alignments(
            paf,
            current_map,
            raw_blocks,
            coverage_path,
            args.min_identity,
            args.min_alignment_bp,
            cycle,
        )
        cycles_run += 1
        total_hits += hits
        total_blocks += blocks
        LOG.info(
            "%s cycle %d: %d accepted hits, %d paired blocks",
            assembly.sample, cycle, hits, blocks,
        )
        if hits == 0 or blocks == 0:
            break
        coverage = load_merged_coverage(
            coverage_path,
            os.path.join(cycle_dir, "coverage.sorted.tsv"),
            cycle_dir,
        )
        covered_bp = sum(
            end - start
            for intervals in coverage.values()
            for start, end in intervals
        )
        if covered_bp == 0:
            break
        next_fasta = os.path.join(cycle_dir, "residual.fa")
        next_map, residual_count, residual_bp = write_residual_queries(
            current_fasta,
            current_map,
            coverage,
            next_fasta,
            args.min_residual_unmasked,
            cycle + 1,
        )
        LOG.info(
            "%s cycle %d: removed %d query bp; %d exposed residuals "
            "(%d bp) continue",
            assembly.sample, cycle, covered_bp, residual_count, residual_bp,
        )
        if not next_map:
            break
        current_fasta = next_fasta
        current_map = next_map
        cycle += 1

    rows = group_blocks_to_bed(
        raw_blocks,
        os.path.join(sample_dir, "all_blocks.sorted.tsv"),
        os.path.join(sample_dir, "aligned.unsorted.bed"),
        shard_bed,
        target_fasta,
        target_stats,
        assembly.sample,
        args.merge_gap,
        args.target_unmasked_cap,
        sample_dir,
    )
    marker = {
        **expected,
        "cycles_run": cycles_run,
        "alignment_hits": total_hits,
        "paired_blocks": total_blocks,
        "bed_rows": rows,
    }
    atomic_json(marker_path, marker)
    if not args.keep_work:
        compact_sample_work(sample_dir)
    LOG.info(
        "%s complete: %d BED rows from %d hits across %d cycles",
        assembly.sample, rows, total_hits, cycles_run,
    )
    return marker


def run_with_retries(
    assemblies: Sequence[AssemblyEntry],
    worker,
    jobs: int,
    retries: int,
) -> Dict[str, dict]:
    pending = list(assemblies)
    results: Dict[str, dict] = {}
    concurrency = max(1, min(jobs, len(pending))) if pending else 1
    for attempt in range(retries + 1):
        if not pending:
            break
        LOG.info(
            "Assembly mapping attempt %d/%d: %d jobs in parallel",
            attempt + 1, retries + 1, concurrency,
        )
        failed: List[AssemblyEntry] = []
        failures: List[Tuple[AssemblyEntry, Exception]] = []
        with ThreadPoolExecutor(max_workers=concurrency) as executor:
            futures = {
                executor.submit(worker, assembly): assembly
                for assembly in pending
            }
            for future in as_completed(futures):
                assembly = futures[future]
                try:
                    results[assembly.sample] = future.result()
                except Exception as error:
                    failed.append(assembly)
                    failures.append((assembly, error))
                    LOG.error("%s failed: %s", assembly.sample, error)
        if not failed:
            break
        if attempt >= retries:
            detail = "; ".join(
                f"{assembly.sample}: {error}"
                for assembly, error in failures
            )
            raise RuntimeError(
                f"{len(failed)} assemblies failed after {retries + 1} "
                f"attempts: {detail}"
            )
        if any(
            isinstance(error, subprocess.CalledProcessError)
            and error.returncode in {-9, 137}
            for _assembly, error in failures
        ):
            concurrency = max(1, concurrency // 2)
            LOG.warning(
                "Possible OOM: retrying with %d concurrent jobs", concurrency,
            )
        pending = failed
    return results


def merge_shard_beds(
    assemblies: Sequence[AssemblyEntry],
    workdir: str,
    output_bed: str,
) -> int:
    Path(output_bed).parent.mkdir(parents=True, exist_ok=True)
    temporary = output_bed + ".tmp"
    rows = 0
    try:
        with open(temporary, "w") as output:
            for assembly in assemblies:
                shard = os.path.join(
                    workdir, "assemblies", sanitize_id(assembly.sample),
                    "aligned.bed",
                )
                if not os.path.isfile(shard):
                    raise RuntimeError(f"Missing completed BED shard: {shard}")
                with open(shard) as source:
                    for raw in source:
                        if raw.strip():
                            rows += 1
                            output.write(raw)
        os.replace(temporary, output_bed)
    finally:
        try:
            os.remove(temporary)
        except FileNotFoundError:
            pass
    columns_path = output_bed + ".columns.txt"
    with open(columns_path, "w") as columns:
        for index, name in enumerate(BED_COLUMNS, 1):
            columns.write(f"{index}\t{name}\n")
    return rows


def parse_args(argv: Optional[Sequence[str]] = None):
    parser = argparse.ArgumentParser(
        description=(
            "Cycle one query FASTA against one blocked-reference-plus-novels "
            "target FASTA and emit accepted merged intervals as BED."
        ),
    )
    parser.add_argument("-q", "--query", required=True, help="query assembly FASTA")
    parser.add_argument(
        "-t", "--target", required=True,
        help="target FASTA already containing blocked reference plus novels",
    )
    parser.add_argument("-o", "--output", required=True, help="output BED")
    return parser.parse_args(argv)


def main(argv: Optional[Sequence[str]] = None) -> int:
    args = parse_args(argv)
    if shutil.which("minimap2") is None:
        raise FileNotFoundError("Required executable not found on PATH: minimap2")
    if shutil.which("sort") is None:
        raise FileNotFoundError("Required executable not found on PATH: sort")
    args.query = os.path.abspath(args.query)
    args.target = os.path.abspath(args.target)
    args.output = os.path.abspath(args.output)
    for path in (args.query, args.target):
        if not os.path.isfile(path):
            raise FileNotFoundError(path)
    if os.path.realpath(args.output) in {
        os.path.realpath(args.query), os.path.realpath(args.target),
    }:
        raise ValueError("-o/--output must differ from -q and -t")

    # The command line intentionally exposes only q/t/o. All behavior below is
    # the fixed policy requested for this mapper.
    args.threads = DEFAULT_THREADS
    args.min_identity = MIN_IDENTITY
    args.min_alignment_bp = MIN_ALIGNMENT_BP
    args.min_residual_unmasked = MIN_RESIDUAL_UNMASKED
    args.merge_gap = MERGE_GAP
    args.target_unmasked_cap = TARGET_UNMASKED_CAP
    args.max_cycles = MAX_CYCLES
    args.resume = False
    args.keep_work = True

    workdir = args.output + ".work"
    shutil.rmtree(workdir, ignore_errors=True)
    Path(workdir).mkdir(parents=True, exist_ok=True)
    setup_logging(workdir, verbose=False)
    LOG.info(
        "Mapping one query with %d Minimap2 threads; -N=%d",
        args.threads, MINIMAP_PARAMS["N"],
    )
    success = False
    try:
        target_stats = scan_target_stats(args.target)
        target_index = build_minimap2_index(
            args.target,
            os.path.join(workdir, "target.asm5.mmi"),
            LOG,
            preset=MINIMAP_PARAMS["preset"],
        )
        config_signature = hash_payload({
            "target": file_fingerprint(args.target),
            "minimap": MINIMAP_PARAMS,
            "min_identity": args.min_identity,
            "min_alignment_bp": args.min_alignment_bp,
            "min_residual_unmasked": args.min_residual_unmasked,
            "merge_gap": args.merge_gap,
            "target_unmasked_cap": args.target_unmasked_cap,
            "max_cycles": args.max_cycles,
            "format": 1,
        })
        query_name = os.path.basename(args.query)
        for suffix in (".fasta", ".fa", ".fna"):
            if query_name.lower().endswith(suffix):
                query_name = query_name[:-len(suffix)]
                break
        assembly = AssemblyEntry(sanitize_id(query_name), args.query)
        sample_dir = os.path.join(workdir, "query")
        marker = map_one_assembly(
            assembly,
            target_index,
            args.target,
            target_stats,
            sample_dir,
            args,
            config_signature,
        )
        source_bed = os.path.join(sample_dir, "aligned.bed")
        Path(args.output).parent.mkdir(parents=True, exist_ok=True)
        temporary = args.output + ".tmp"
        shutil.copyfile(source_bed, temporary)
        os.replace(temporary, args.output)
        with open(args.output + ".columns.txt", "w") as columns:
            for index, name in enumerate(BED_COLUMNS, 1):
                columns.write(f"{index}\t{name}\n")
        success = True
        LOG.info(
            "Done: %d accepted merged intervals -> %s. Run "
            "filter_block_mapping_bed.py to apply multi-mapping filters.",
            marker["bed_rows"], args.output,
        )
    finally:
        # Keep every PAF and residual when a job fails. On success the requested
        # BED and its column definition are the only products.
        if success:
            for handler in list(LOG.handlers):
                handler.close()
                LOG.removeHandler(handler)
            shutil.rmtree(workdir, ignore_errors=True)
    return 0


if __name__ == "__main__":
    try:
        raise SystemExit(main())
    except Exception as error:
        LOG.error("%s", error)
        raise
