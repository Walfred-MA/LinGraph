#!/usr/bin/env python3
"""Align one assembly's KmerSearcher hotspots to their partition graphs.

Hotspot column 2 is a 1-based index into ``-L/--graph-list``. Each graph-list
entry supplies a graph name, which is resolved below ``-G/--graph-folder`` as
``GRAPH_NAME/GRAPH_NAME.FA``. Raw hotspot intervals are merged per
partition/contig, extended, extracted from the source assembly, and aligned as
complete sequences. Output graph CIGARs are never split through a graph
``cache.json``. The source contigs needed by the hotspot report are loaded into
RAM once, all initial query FASTAs are materialized once, and only then are the
independent FASTA-to-graph alignments dispatched to worker processes. Each
worker optionally finishes all cycle passes for its current hotspot before
accepting another hotspot job.
"""

from __future__ import annotations

import argparse
import concurrent.futures
import dataclasses
import fcntl
import hashlib
import importlib.util
import logging
import multiprocessing as mp
import os
import re
import shutil
import subprocess
import tempfile
import threading
from collections import Counter, defaultdict
from pathlib import Path
from typing import (
    Any, DefaultDict, Dict, Iterable, Iterator, List, Mapping, Optional,
    Sequence, Tuple,
)

from local_graph_whole_cigar import (
    BLAST_IDENTITY,
    WholeGraphCigar,
    add_query_payloads,
    align_whole_graph_cigar,
    blast_database_is_complete,
    ensure_graph_blast_database,
    fasta_sequence_lengths,
    graph_alignment_tool,
    insertion_only,
    normalize_graphcigarlight_sam,
    normalize_whole_graph_cigar,
    query_strand,
    resolve_extension_dir,
    run_graph_alignment_command,
)
from minsetref_core import IndexedFasta, wrap_fasta
from graph_cigar_payloads import strip_graph_cigar_payloads


LOG = logging.getLogger("align_partition_hotspots")
DEFAULT_EXTENSION = 30_000
DEFAULT_MERGE_GAP = 2 * DEFAULT_EXTENSION
DEFAULT_HOTSPOT_KMER_LENGTH = 31
DEFAULT_CYCLE_MIN_UNALIGNED = 300
FAST_MAX_UNMAPPED_GAP = 50
FAST_MASKED_BASE_THRESHOLD = 1_000
FAST_WINNOWMAP_BATCH_BASES = 256 * 1024 * 1024
BLAST_LOCAL_ID_MAXIMUM = 50
TIMEOUT_EXIT_CODE = 124
DIRECT_WORK_PROTOCOL = "direct-alignment-work-v3"
FAST_COMMAND_COMPLETE_PROTOCOL = "fast-alignment-command-v1\n"
CIGAR_TOKEN_RE = re.compile(r"([><])([^:<>]+):([^<>]+)")
CIGAR_OPERATION_RE = re.compile(r"(\d+)([HID=XMN])([A-Za-z]*)")
GRAPH_NAME_SUFFIXES = (
    "_samples.fasta",
    "_samples.fa",
    ".fasta",
    ".list",
    ".txt",
    ".tsv",
    ".bed",
    ".FA",
    ".fa",
)
FASTA_SUFFIXES = (".fa", ".fasta", ".fna", ".fas", ".FA")
SOURCE_COORDINATE_RE = re.compile(r"^(.+):(\d+)-(\d+)([+-])?$")
COORDINATE_NAME_RE = re.compile(r"^(.+)_(\d+)_(\d+)$")
UNIQUE_REGION_RE = re.compile(r"^unique_region_\d+_(.+)_(\d+)_(\d+)$")
ALIGNMENT_HEADER_FIELDS = (
    "hotspot_index", "contig", "start", "end", "strand", "queryname",
    "graph_path", "graphcigar", "refpositions", "qpositions",
)


def alignment_mode_protocol(args: argparse.Namespace) -> str:
    """Return the final-output mode marker used for safe resume decisions."""
    if bool(getattr(args, "slow_rigorous", False)):
        mode = "slow-rigorous-v1;blastn=all;winnowmap=all"
    elif bool(getattr(args, "fast_mode", False)):
        mode = "fast-mode-v2;max-unmapped-gap=50"
    elif bool(getattr(args, "cycle", False)):
        mode = "cycle-mode-v1;min-unaligned=300"
    else:
        mode = "normal-mode-v1"
    return mode + f";payload={int(bool(getattr(args, 'alignment_payload', False)))}\n"


def atomic_write_text(path: str, text: str) -> None:
    """Atomically replace a small checkpoint file."""
    destination = Path(path)
    destination.parent.mkdir(parents=True, exist_ok=True)
    temporary = destination.with_name(
        destination.name + f".tmp.{os.getpid()}"
    )
    try:
        temporary.write_text(text, encoding="utf-8")
        os.replace(temporary, destination)
    finally:
        try:
            temporary.unlink()
        except FileNotFoundError:
            pass


def completed_command_output(path: str) -> bool:
    """Return true only for an output committed by a finished command."""
    marker = path + ".complete"
    try:
        return (
            os.path.isfile(path)
            and Path(marker).read_text(encoding="utf-8")
            == FAST_COMMAND_COMPLETE_PROTOCOL
        )
    except FileNotFoundError:
        return False


def mark_command_output_complete(path: str) -> None:
    atomic_write_text(path + ".complete", FAST_COMMAND_COMPLETE_PROTOCOL)


@dataclasses.dataclass(frozen=True)
class SourceInterval:
    contig: str
    start: int
    end: int


@dataclasses.dataclass(frozen=True)
class PartitionTask:
    hotspot_index: int
    graph_fasta: str
    intervals: Tuple[SourceInterval, ...]


@dataclasses.dataclass(frozen=True)
class PreparedAlignment:
    hotspot_index: int
    graph_fasta: str
    interval: SourceInterval
    query_name: str
    query_fasta: str
    result_path: str
    source_strand: str = "+"


@dataclasses.dataclass(frozen=True)
class OutputRow:
    hotspot_index: int
    contig: str
    start: int
    end: int
    strand: str
    query_name: str
    graph_path: str
    graph_cigar: str
    ref_positions: str
    query_positions: str

    def line(self) -> str:
        return "\t".join(str(value) for value in dataclasses.astuple(self)) + "\n"


@dataclasses.dataclass(frozen=True)
class AlignmentJobResult:
    output_rows: int
    insertion_rows: int
    additional_rows: int = 0
    cycle_passes: int = 0
    recycled_attempts: int = 0


@dataclasses.dataclass(frozen=True)
class DirectSequence:
    sample: str
    contig: str
    start: int
    end: int
    strand: str
    query_name: str
    sequence: str
    record_class: str = ""


@dataclasses.dataclass(frozen=True)
class FastPreparedAlignment:
    prepared: PreparedAlignment
    sequence: str


_WORKER_THREADS = 1
_WORKER_TIMEOUT: Optional[int] = None
_WORKER_CYCLE = False
_WORKER_OUTPUT_PATH: Optional[str] = None
_WORKER_PAYLOAD = False


def graph_name_from_list_entry(value: str, context: str) -> str:
    """Return the final path component with one recognized suffix removed."""
    final_component = value.rstrip("/").rsplit("/", 1)[-1]
    if not final_component:
        raise ValueError(f"{context}: graph-list entry has no final path component")
    lowered = final_component.lower()
    for suffix in GRAPH_NAME_SUFFIXES:
        if lowered.endswith(suffix.lower()):
            final_component = final_component[:-len(suffix)]
            break
    if not final_component:
        raise ValueError(f"{context}: graph name is empty after removing its suffix")
    return final_component


def read_graph_list(path: str, graph_folder: str) -> List[Optional[str]]:
    graphs: List[Optional[str]] = []
    missing: List[str] = []
    missing_blastdb: List[str] = []
    graph_root = os.path.abspath(os.path.expanduser(graph_folder))
    if not os.path.isdir(graph_root):
        raise NotADirectoryError(f"graph folder not found: {graph_root}")
    with open(path, "rt") as handle:
        for line_number, raw in enumerate(handle, 1):
            content = raw.strip()
            if not content or content.startswith("#"):
                continue
            graph_name = graph_name_from_list_entry(
                content.split()[0], f"{path}:{line_number}",
            )
            graph_fasta = os.path.join(
                graph_root, graph_name, graph_name + ".FA",
            )
            if not os.path.isfile(graph_fasta):
                graphs.append(None)
                missing.append(graph_fasta)
            elif not blast_database_is_complete(graph_fasta):
                graphs.append(None)
                missing_blastdb.append(graph_fasta + "_db")
            else:
                graphs.append(graph_fasta)
    if not graphs:
        raise ValueError(f"graph list is empty: {path}")
    if missing:
        examples = ", ".join(missing[:5])
        remainder = f" (and {len(missing) - 5} more)" if len(missing) > 5 else ""
        LOG.warning(
            "%d/%d graph-list entries have no local graph and will be skipped: %s%s",
            len(missing), len(graphs), examples, remainder,
        )
    if missing_blastdb:
        examples = ", ".join(missing_blastdb[:5])
        remainder = (
            f" (and {len(missing_blastdb) - 5} more)"
            if len(missing_blastdb) > 5 else ""
        )
        LOG.warning(
            "%d/%d graph-list entries have no complete v4 BLAST database; "
            "their alignments will be skipped: %s%s",
            len(missing_blastdb), len(graphs), examples, remainder,
        )
    return graphs


def parse_hotspot_row(
    fields: Sequence[str], context: str,
) -> Tuple[int, str, int, int, bool]:
    """Parse either legacy four-column or explicit six-column hotspot rows."""
    try:
        if len(fields) >= 6 and fields[0].isdigit() and fields[3] in {"+", "-"}:
            hotspot_index = int(fields[0])
            contig = fields[2]
            start, end = int(fields[4]), int(fields[5])
            raw_kmer_coordinates = False
        elif len(fields) >= 4 and fields[1].isdigit():
            contig = fields[0]
            hotspot_index = int(fields[1])
            start, end = int(fields[-2]), int(fields[-1])
            raw_kmer_coordinates = True
        else:
            raise ValueError("unrecognized hotspot layout")
    except (IndexError, ValueError) as error:
        raise ValueError(f"{context}: cannot parse hotspot row: {' '.join(fields)!r}") from error
    if hotspot_index < 1:
        raise ValueError(f"{context}: hotspot index must be 1-based and positive")
    start, end = sorted((start, end))
    if start < 0 or end <= start:
        raise ValueError(f"{context}: invalid hotspot interval {start}-{end}")
    return hotspot_index, contig, start, end, raw_kmer_coordinates


def merge_intervals(
    intervals: Iterable[Tuple[int, int]], gap: int,
) -> List[Tuple[int, int]]:
    ordered = sorted(intervals)
    if not ordered:
        return []
    merged: List[List[int]] = [[ordered[0][0], ordered[0][1]]]
    for start, end in ordered[1:]:
        previous = merged[-1]
        if start <= previous[1] + gap:
            previous[1] = max(previous[1], end)
        else:
            merged.append([start, end])
    return [(start, end) for start, end in merged]


def build_partition_tasks(
    hotspot_path: str,
    graph_fastas: Sequence[Optional[str]],
    fasta: IndexedFasta,
    extension: int,
    merge_gap: int,
    kmer_length: int = DEFAULT_HOTSPOT_KMER_LENGTH,
) -> Tuple[List[PartitionTask], int]:
    raw: DefaultDict[int, DefaultDict[str, List[Tuple[int, int]]]] = defaultdict(
        lambda: defaultdict(list)
    )
    contig_order: Dict[str, int] = {}
    input_rows = 0
    available = set(fasta.names())
    with open(hotspot_path, "rt") as handle:
        for line_number, line in enumerate(handle, 1):
            if not line.strip() or line.startswith("#"):
                continue
            (
                hotspot_index, contig, start, end, raw_kmer_coordinates,
            ) = parse_hotspot_row(
                line.split(), f"{hotspot_path}:{line_number}",
            )
            if hotspot_index > len(graph_fastas):
                raise IndexError(
                    f"{hotspot_path}:{line_number}: hotspot index {hotspot_index} "
                    f"exceeds graph-list length {len(graph_fastas)}"
                )
            input_rows += 1
            if graph_fastas[hotspot_index - 1] is None:
                continue
            if raw_kmer_coordinates:
                # KmerSearcher reports the terminal coordinate of its first
                # k-mer. Match build_local_graphs.py and recover the complete
                # zero-based k-mer span before merging/extension.
                start = max(0, start - kmer_length)
            if contig not in available:
                raise KeyError(
                    f"{hotspot_path}:{line_number}: contig {contig!r} is absent "
                    f"from {fasta.path}"
                )
            contig_length = fasta.length(contig)
            if start >= contig_length:
                raise ValueError(
                    f"{hotspot_path}:{line_number}: start {start} is outside "
                    f"{contig!r} length {contig_length}"
                )
            end = min(end, contig_length)
            raw[hotspot_index][contig].append((start, end))
            contig_order.setdefault(contig, len(contig_order))

    tasks: List[PartitionTask] = []
    for hotspot_index in sorted(raw):
        graph_fasta = graph_fastas[hotspot_index - 1]
        if graph_fasta is None:
            continue
        extracted: List[SourceInterval] = []
        for contig in sorted(raw[hotspot_index], key=lambda name: contig_order[name]):
            contig_length = fasta.length(contig)
            for start, end in merge_intervals(raw[hotspot_index][contig], merge_gap):
                extracted.append(SourceInterval(
                    contig,
                    max(0, start - extension),
                    min(contig_length, end + extension),
                ))
        tasks.append(PartitionTask(
            hotspot_index,
            graph_fasta,
            tuple(extracted),
        ))
    return tasks, input_rows


def load_hotspot_contigs(
    fasta: IndexedFasta, tasks: Sequence[PartitionTask],
) -> Dict[str, str]:
    """Load every assembly contig referenced by retained hotspots exactly once."""
    ordered_names = dict.fromkeys(
        interval.contig
        for task in tasks
        for interval in task.intervals
    )
    return {name: fasta.sequence(name) for name in ordered_names}


def prepare_alignments(
    tasks: Sequence[PartitionTask],
    contigs: Mapping[str, str],
    sample: str,
    work_root: str,
) -> List[PreparedAlignment]:
    """Write all hotspot FASTAs before any graph aligner is started."""
    prepared: List[PreparedAlignment] = []
    for task in tasks:
        graph_name = graph_name_from_list_entry(
            task.graph_fasta, f"hotspot {task.hotspot_index} graph FASTA",
        )
        graph_query_prefix = graph_name.split("_", 1)[0]
        task_dir = os.path.join(work_root, f"hotspot_{task.hotspot_index:08d}")
        Path(task_dir).mkdir(parents=True, exist_ok=True)
        for local_index, interval in enumerate(task.intervals, 1):
            # Downstream PA names begin with the graph's basename-before-first-
            # underscore prefix.  Emit that stable identity at the source so
            # GenomeLift alignment-score priority and graph-CIGAR lookup both
            # work without a post-hoc prefix repair.
            query_name = f"{graph_query_prefix}_{sample}_{local_index}"
            query_fasta = os.path.join(task_dir, f"query_{local_index:06d}.fa")
            result_path = os.path.join(
                task_dir, f"query_{local_index:06d}.graph_cigar.tsv",
            )
            sequence = contigs[interval.contig][interval.start:interval.end]
            expected_length = interval.end - interval.start
            if len(sequence) != expected_length:
                raise ValueError(
                    f"hotspot {task.hotspot_index}: extracted "
                    f"{interval.contig}:{interval.start}-{interval.end} has "
                    f"{len(sequence)} bases instead of {expected_length}"
                )
            with open(query_fasta, "wt") as output:
                output.write(
                    f">{query_name}\t{interval.contig}:"
                    f"{interval.start}-{interval.end}+\n"
                )
                output.write(wrap_fasta(sequence) + "\n")
            prepared.append(PreparedAlignment(
                task.hotspot_index,
                task.graph_fasta,
                interval,
                query_name,
                query_fasta,
                result_path,
            ))
    return prepared


def _worker_init(
    threads: int,
    timeout: Optional[int],
    cycle: bool = False,
    output_path: Optional[str] = None,
    payload: bool = False,
) -> None:
    global _WORKER_THREADS, _WORKER_TIMEOUT, _WORKER_CYCLE
    global _WORKER_OUTPUT_PATH, _WORKER_PAYLOAD
    _WORKER_THREADS = threads
    _WORKER_TIMEOUT = timeout
    _WORKER_CYCLE = cycle
    _WORKER_OUTPUT_PATH = output_path
    _WORKER_PAYLOAD = payload


def _append_completed_job(rows: Sequence[OutputRow]) -> None:
    """Append one fully cycled job to the shared temporary output."""
    if not _WORKER_OUTPUT_PATH:
        raise RuntimeError("worker output path was not initialized")
    payload = "".join(row.line() for row in rows)
    with open(_WORKER_OUTPUT_PATH, "at") as output:
        fcntl.flock(output.fileno(), fcntl.LOCK_EX)
        try:
            output.write(payload)
            output.flush()
        finally:
            fcntl.flock(output.fileno(), fcntl.LOCK_UN)


def _align_prepared(
    prepared: PreparedAlignment, *, winnowmap_only: bool = False,
) -> OutputRow:
    result, strand = align_whole_graph_cigar(
        prepared.query_fasta,
        prepared.graph_fasta,
        prepared.result_path,
        threads=_WORKER_THREADS,
        timeout=_WORKER_TIMEOUT,
        winnowmap_only=winnowmap_only,
        # The legacy converter's --addinsert also resolves M into =/X.
        # Preserve that classification even when saving a compact final row.
        add_payload=True,
    )
    if not _WORKER_PAYLOAD:
        result = dataclasses.replace(
            result, graph_cigar=strip_graph_cigar_payloads(result.graph_cigar),
        )
    if prepared.source_strand == "-":
        strand = "+" if strand == "-" else "-"
    interval = prepared.interval
    return OutputRow(
        prepared.hotspot_index,
        interval.contig,
        interval.start,
        interval.end,
        strand,
        prepared.query_name,
        result.graph_path,
        result.graph_cigar,
        result.ref_positions,
        result.query_positions,
    )


def read_single_fasta_sequence(path: str) -> str:
    """Read one prepared FASTA before KmerStrd reorients it in place."""
    records = 0
    pieces: List[str] = []
    with open(path, "rt") as handle:
        for line_number, raw in enumerate(handle, 1):
            if raw.startswith(">"):
                records += 1
                if records > 1:
                    raise ValueError(f"{path}: expected exactly one FASTA record")
            elif raw.strip():
                if records != 1:
                    raise ValueError(
                        f"{path}:{line_number}: sequence before FASTA header"
                    )
                pieces.append(raw.strip())
    if records != 1:
        raise ValueError(f"{path}: expected exactly one FASTA record")
    return "".join(pieces)


def prepare_worker_cycle_alignments(
    rows: Sequence[OutputRow],
    root: PreparedAlignment,
    root_sequence: str,
    cycle: int,
    minimum_unaligned: int = DEFAULT_CYCLE_MIN_UNALIGNED,
) -> List[PreparedAlignment]:
    """Prepare one original hotspot job's residual chunks inside its worker."""
    expected_length = root.interval.end - root.interval.start
    if len(root_sequence) != expected_length:
        raise ValueError(
            f"hotspot {root.hotspot_index}: original worker sequence has "
            f"{len(root_sequence)} bases instead of {expected_length}"
        )
    cycle_root = root.result_path + ".cycles"
    prepared: List[PreparedAlignment] = []
    for row_index, row in enumerate(rows, 1):
        for chunk_index, interval in enumerate(
            unaligned_query_intervals(row), 1,
        ):
            length = interval.end - interval.start
            if length <= minimum_unaligned:
                continue
            if (
                interval.contig != root.interval.contig
                or interval.start < root.interval.start
                or interval.end > root.interval.end
            ):
                raise ValueError(
                    f"hotspot {root.hotspot_index}: recycled interval "
                    f"{interval.contig}:{interval.start}-{interval.end} is "
                    f"outside original interval {root.interval.contig}:"
                    f"{root.interval.start}-{root.interval.end}"
                )
            if root.source_strand == "+":
                local_start = interval.start - root.interval.start
                local_end = interval.end - root.interval.start
            else:
                local_start = root.interval.end - interval.end
                local_end = root.interval.end - interval.start
            sequence = root_sequence[local_start:local_end]
            if len(sequence) != length:
                raise ValueError(
                    f"cycle {cycle}: extracted {interval.contig}:"
                    f"{interval.start}-{interval.end} has {len(sequence)} "
                    f"bases instead of {length}"
                )
            task_dir = os.path.join(
                cycle_root,
                f"cycle_{cycle:04d}",
                f"row_{row_index:08d}",
            )
            Path(task_dir).mkdir(parents=True, exist_ok=True)
            query_fasta = os.path.join(
                task_dir, f"chunk_{chunk_index:06d}.fa",
            )
            result_path = os.path.join(
                task_dir, f"chunk_{chunk_index:06d}.graph_cigar.tsv",
            )
            with open(query_fasta, "wt") as output:
                output.write(
                    f">{root.query_name}\t{interval.contig}:"
                    f"{interval.start}-{interval.end}+\n"
                )
                output.write(wrap_fasta(sequence) + "\n")
            prepared.append(PreparedAlignment(
                root.hotspot_index,
                root.graph_fasta,
                interval,
                root.query_name,
                query_fasta,
                result_path,
                root.source_strand,
            ))
    return prepared


def _align_prepared_job(prepared: PreparedAlignment) -> AlignmentJobResult:
    """Align one original hotspot and finish all of its cycles in this job."""
    root_sequence = (
        read_single_fasta_sequence(prepared.query_fasta)
        if _WORKER_CYCLE else ""
    )
    initial = _align_prepared(prepared)
    if not _WORKER_CYCLE:
        rows = [initial]
        _append_completed_job(rows)
        return AlignmentJobResult(1, int(initial.graph_path == "*"))

    rows = [initial]
    frontier = [initial]
    cycle = 1
    recycled_attempts = 0
    attempted_passes = 0
    while True:
        recycled = prepare_worker_cycle_alignments(
            frontier, prepared, root_sequence, cycle,
        )
        if not recycled:
            break
        attempted_passes += 1
        recycled_attempts += len(recycled)
        additional = []
        for child in recycled:
            row = _align_prepared(child, winnowmap_only=True)
            if row.graph_path and row.graph_path != "*":
                additional.append(row)
        if not additional:
            break
        rows.extend(additional)
        frontier = additional
        cycle += 1
    _append_completed_job(rows)
    return AlignmentJobResult(
        output_rows=len(rows),
        insertion_rows=sum(row.graph_path == "*" for row in rows),
        additional_rows=len(rows) - 1,
        cycle_passes=attempted_passes,
        recycled_attempts=recycled_attempts,
    )


def align_prepared_batch(
    prepared: Sequence[PreparedAlignment],
    args: argparse.Namespace,
    label: str,
    output_path: str,
) -> Tuple[int, int]:
    """Align one materialized batch within the configured CPU budget."""
    if not prepared:
        return 0, 0
    automatic_jobs = max(1, args.threads // args.threads_per_job)
    jobs = args.jobs if args.jobs is not None else automatic_jobs
    jobs = max(1, min(jobs, len(prepared), automatic_jobs))
    cycle_enabled = getattr(args, "cycle", False)
    cycle_note = (
        "; cycles use Winnowmap only and stay within each job"
        if cycle_enabled else ""
    )
    LOG.info(
        "Aligning %d %s FASTAs with %d parallel jobs x %d aligner threads%s",
        len(prepared), label, jobs, args.threads_per_job, cycle_note,
    )
    initializer_args = (
        args.threads_per_job, args.timeout, cycle_enabled, output_path,
        bool(getattr(args, "alignment_payload", False)),
    )
    output_rows = 0
    insertion_rows = 0
    additional_rows = 0
    cycle_passes = 0
    recycled_attempts = 0

    def collect(result: AlignmentJobResult) -> None:
        nonlocal output_rows, insertion_rows, additional_rows
        nonlocal cycle_passes, recycled_attempts
        output_rows += result.output_rows
        insertion_rows += result.insertion_rows
        additional_rows += result.additional_rows
        cycle_passes += result.cycle_passes
        recycled_attempts += result.recycled_attempts

    def report_progress(completed: int) -> None:
        if (
            completed == len(prepared)
            or completed % max(1, len(prepared) // 20) == 0
        ):
            LOG.info(
                "%s alignment progress: %d/%d FASTAs",
                label.capitalize(), completed, len(prepared),
            )

    if jobs == 1:
        _worker_init(*initializer_args)
        for completed, alignment in enumerate(prepared, 1):
            collect(_align_prepared_job(alignment))
            report_progress(completed)
    else:
        # The parent may retain several GiB of assembly sequence.  Never use
        # Linux's implicit fork context here: workers need clean interpreters,
        # not inherited memory, file descriptors, and lock state.
        with concurrent.futures.ProcessPoolExecutor(
            max_workers=jobs,
            mp_context=mp.get_context("spawn"),
            initializer=_worker_init,
            initargs=initializer_args,
        ) as executor:
            for completed, result in enumerate(
                executor.map(_align_prepared_job, prepared), 1,
            ):
                collect(result)
                report_progress(completed)

    if cycle_enabled:
        LOG.info(
            "Worker-local cycle alignment attempted %d residual chunks "
            "across %d per-job cycle passes and appended %d additional rows",
            recycled_attempts, cycle_passes, additional_rows,
        )
    return output_rows, insertion_rows


def parse_position_list(text: str, context: str) -> List[Tuple[int, int]]:
    if not text or text == ".":
        return []
    output: List[Tuple[int, int]] = []
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


def unaligned_query_intervals(row: OutputRow) -> List[SourceInterval]:
    """Return maximal insertion spans in original assembly coordinates."""
    query_length = row.end - row.start
    if row.graph_path == "*" or not row.graph_path:
        return [SourceInterval(row.contig, row.start, row.end)]

    cigar_tokens = CIGAR_TOKEN_RE.findall(row.graph_cigar)
    query_positions = parse_position_list(
        row.query_positions, f"hotspot {row.hotspot_index}:qpositions",
    )
    if len(cigar_tokens) != len(query_positions):
        raise ValueError(
            f"hotspot {row.hotspot_index}: graph CIGAR/query position "
            f"counts differ: {len(cigar_tokens)}/{len(query_positions)}"
        )

    oriented: List[Tuple[int, int]] = []
    for component_index, (
        (_marker, _path_name, cigar), (query_start, query_end),
    ) in enumerate(zip(cigar_tokens, query_positions), 1):
        operations = CIGAR_OPERATION_RE.findall(cigar)
        if not operations or "".join(
            f"{length}{operation}{payload}"
            for length, operation, payload in operations
        ) != cigar:
            raise ValueError(
                f"hotspot {row.hotspot_index} component {component_index}: "
                f"invalid graph CIGAR {cigar!r}"
            )
        cursor = query_start
        for length_text, operation, payload in operations:
            length = int(length_text)
            if payload and operation not in {"I", "X"}:
                raise ValueError(
                    f"hotspot {row.hotspot_index} component {component_index}: "
                    f"unexpected payload on {operation}"
                )
            if operation in {"I", "X"} and payload and len(payload) != length:
                raise ValueError(
                    f"hotspot {row.hotspot_index} component {component_index}: "
                    f"payload length differs from {length}{operation}"
                )
            if operation == "I":
                oriented.append((cursor, cursor + length))
            if operation in {"I", "=", "X", "M"}:
                cursor += length
        if cursor != query_end:
            raise ValueError(
                f"hotspot {row.hotspot_index} component {component_index}: "
                f"CIGAR consumes {cursor - query_start} query bases but "
                f"qpositions spans {query_end - query_start}"
            )

    merged = merge_intervals(oriented, gap=0)
    output = []
    for start, end in merged:
        if not (0 <= start < end <= query_length):
            raise ValueError(
                f"hotspot {row.hotspot_index}: unaligned query interval "
                f"{start}-{end} exceeds query length {query_length}"
            )
        if row.strand == "+":
            absolute_start, absolute_end = row.start + start, row.start + end
        else:
            absolute_start, absolute_end = row.end - end, row.end - start
        output.append(SourceInterval(
            row.contig, absolute_start, absolute_end,
        ))
    return sorted(output, key=lambda value: (value.contig, value.start, value.end))


def _sample_from_record_id(record_id: str, fallback: str = "") -> str:
    fields = record_id.split("_")
    if len(fields) >= 3 and fields[1] and fields[2]:
        return f"{fields[1]}_{fields[2]}"
    return fallback or "sample"


def read_direct_fasta(path: str, sample: str = "") -> List[DirectSequence]:
    """Read all FASTA records and retain optional source-coordinate headers."""
    records: List[DirectSequence] = []
    description: Optional[str] = None
    pieces: List[str] = []

    def finish() -> None:
        if description is None:
            return
        fields = description.split()
        record_id = fields[0]
        sequence = "".join(pieces)
        if not sequence:
            raise ValueError(f"{path}: FASTA record {record_id!r} is empty")
        record_sample = _sample_from_record_id(record_id, sample)
        record_class = next(
            (token.split("=", 1)[1] for token in fields[1:] if token.startswith("class=")),
            "",
        )
        contig, start, end, source_strand = record_id, 0, len(sequence), "+"
        for token in fields[1:]:
            match = SOURCE_COORDINATE_RE.fullmatch(token)
            if match is None:
                continue
            candidate_contig, start_text, end_text, candidate_strand = match.groups()
            coordinate_sample = ""
            if ":" in candidate_contig:
                coordinate_sample, candidate_contig = candidate_contig.split(":", 1)
            candidate_start, candidate_end = int(start_text), int(end_text)
            if candidate_end - candidate_start == len(sequence):
                contig, start, end = candidate_contig, candidate_start, candidate_end
                source_strand = candidate_strand or "+"
                if coordinate_sample:
                    record_sample = coordinate_sample
                break
        records.append(DirectSequence(
            record_sample, contig, start, end,
            source_strand, record_id, sequence, record_class,
        ))

    with open(path, "rt") as handle:
        for line_number, raw in enumerate(handle, 1):
            if raw.startswith(">"):
                finish()
                description = raw[1:].strip()
                pieces = []
                if not description:
                    raise ValueError(f"{path}:{line_number}: empty FASTA header")
            elif raw.strip():
                if description is None:
                    raise ValueError(f"{path}:{line_number}: sequence before header")
                pieces.append(raw.strip())
    finish()
    if not records:
        raise ValueError(f"{path}: no FASTA records")
    duplicate = next(
        (name for name, count in Counter(
            record.query_name for record in records
        ).items() if count > 1),
        None,
    )
    if duplicate is not None:
        raise ValueError(f"{path}: duplicate FASTA ID {duplicate!r}")
    return records


def read_query_path_list(path: str) -> Tuple[List[str], Dict[str, str]]:
    order: List[str] = []
    sources: Dict[str, str] = {}
    with open(path, "rt") as handle:
        for line_number, raw in enumerate(handle, 1):
            if not raw.strip() or raw.startswith("#"):
                continue
            fields = raw.split()
            if len(fields) != 2:
                raise ValueError(
                    f"{path}:{line_number}: expected exactly NAME FASTA; index must be FASTA.fai"
                )
            sample, fasta_path = fields[0], os.path.abspath(os.path.expanduser(fields[1]))
            if sample in sources:
                raise ValueError(f"{path}:{line_number}: duplicate sample {sample!r}")
            fai_path = os.path.abspath(
                fasta_path + ".fai"
            )
            if not os.path.isfile(fasta_path) or not os.path.isfile(fai_path):
                raise FileNotFoundError(
                    f"{path}:{line_number}: FASTA and faidx are required: "
                    f"{fasta_path}, {fai_path}"
                )
            order.append(sample)
            sources[sample] = fasta_path
    if not order:
        raise ValueError(f"{path}: no query FASTA sources")
    return order, sources


def _coordinate_named_source(
    name: str, samples: Sequence[str], available: Mapping[str, IndexedFasta],
) -> Optional[Tuple[str, str, int, int]]:
    match = COORDINATE_NAME_RE.fullmatch(name)
    if match is None:
        return None
    prefix, start_text, end_text = match.groups()
    start, end = int(start_text), int(end_text)
    if end <= start:
        return None
    for sample in sorted(samples, key=len, reverse=True):
        marker = sample + "_"
        if not prefix.startswith(marker):
            continue
        contig = prefix[len(marker):]
        if contig in available[sample].index:
            return sample, contig, start, end
    return None


def read_direct_bed(path: str, query_paths: str) -> List[DirectSequence]:
    """Extract BED regions from the indexed sample FASTAs in query_paths."""
    sample_order, sources = read_query_path_list(query_paths)
    readers = {sample: IndexedFasta(source) for sample, source in sources.items()}
    output: List[DirectSequence] = []
    try:
        with open(path, "rt") as handle:
            for line_number, raw in enumerate(handle, 1):
                if not raw.strip() or raw.startswith("#"):
                    continue
                fields = raw.rstrip("\r\n").split("\t")
                if len(fields) < 3:
                    fields = raw.split()
                if len(fields) < 3:
                    raise ValueError(f"{path}:{line_number}: expected BED3+")
                try:
                    bed_start, bed_end = int(fields[1]), int(fields[2])
                except ValueError as error:
                    raise ValueError(
                        f"{path}:{line_number}: invalid BED coordinates"
                    ) from error
                if bed_start < 0 or bed_end <= bed_start:
                    raise ValueError(f"{path}:{line_number}: invalid BED interval")
                sample = fields[6] if len(fields) >= 7 and fields[6] in readers else ""
                contig = fields[0]
                start, end = bed_start, bed_end

                coordinate_source = _coordinate_named_source(
                    contig, sample_order, readers,
                )
                if coordinate_source is not None:
                    encoded_sample, encoded_contig, encoded_start, encoded_end = coordinate_source
                    sample = sample or encoded_sample
                    contig = encoded_contig
                    span = encoded_end - encoded_start
                    if 0 <= bed_start < bed_end <= span:
                        start, end = encoded_start + bed_start, encoded_start + bed_end
                    elif encoded_start <= bed_start < bed_end <= encoded_end:
                        start, end = bed_start, bed_end
                    else:
                        raise ValueError(
                            f"{path}:{line_number}: BED interval lies outside "
                            f"coordinate-named source {encoded_start}-{encoded_end}"
                        )

                if not sample:
                    candidates = [
                        value for value in sample_order
                        if contig in readers[value].index
                    ]
                    if len(candidates) == 1:
                        sample = candidates[0]
                    elif len(sample_order) == 1:
                        sample = sample_order[0]

                # minsetref_light BED uses local record coordinates but stores
                # the absolute source interval in its unique_region name.
                if sample in readers and contig not in readers[sample].index and len(fields) >= 4:
                    match = UNIQUE_REGION_RE.fullmatch(fields[3])
                    if match is not None:
                        named_contig, named_start, named_end = (
                            match.group(1), int(match.group(2)), int(match.group(3))
                        )
                        if named_contig in readers[sample].index:
                            contig, start, end = named_contig, named_start, named_end

                if sample not in readers:
                    raise ValueError(
                        f"{path}:{line_number}: cannot select a query-path sample"
                    )
                reader = readers[sample]
                if contig not in reader.index:
                    raise KeyError(
                        f"{path}:{line_number}: {contig!r} is absent from "
                        f"sample {sample} FASTA {reader.path}"
                    )
                if end > reader.length(contig):
                    raise ValueError(
                        f"{path}:{line_number}: {contig}:{start}-{end} exceeds "
                        f"source length {reader.length(contig)}"
                    )
                strand = fields[5] if len(fields) >= 6 and fields[5] in {"+", "-"} else "+"
                query_name = (
                    fields[3] if len(fields) >= 4 and fields[3] not in {"", "."}
                    else f"{sample}_{contig}_{start}_{end}"
                )
                query_name = re.sub(r"[^A-Za-z0-9_.:#-]+", "_", query_name)
                query_name = f"{query_name}__row{line_number}"
                output.append(DirectSequence(
                    sample, contig, start, end, strand, query_name,
                    reader.fetch(contig, start, end),
                    fields[7] if len(fields) >= 8 and fields[7] in {
                        "reference", "unique_query",
                    } else "",
                ))
    finally:
        for reader in readers.values():
            reader.close()
    if not output:
        raise ValueError(f"{path}: no BED regions")
    return output


def prepare_direct_alignments(
    records: Sequence[DirectSequence], graph_fasta: str, work_root: str,
) -> List[PreparedAlignment]:
    """Materialize one ordered _samples.fasta plus per-record graph-CIGAR jobs."""
    combined = os.path.join(work_root, "input_samples.fasta")
    with open(combined, "wt") as output:
        for record in records:
            output.write(
                f">{record.query_name}\t{record.contig}:{record.start}-"
                f"{record.end}{record.strand}\tsample={record.sample}\n"
                f"{wrap_fasta(record.sequence)}\n"
            )
    prepared: List[PreparedAlignment] = []
    jobs_root = os.path.join(work_root, "records")
    Path(jobs_root).mkdir(parents=True, exist_ok=True)
    for index, record in enumerate(records, 1):
        query_fasta = os.path.join(jobs_root, f"record_{index:08d}.fa")
        result_path = os.path.join(jobs_root, f"record_{index:08d}.graph_cigar.tsv")
        with open(query_fasta, "wt") as output:
            output.write(
                f">{record.query_name}\t{record.contig}:{record.start}-"
                f"{record.end}{record.strand}\tsample={record.sample}\n"
                f"{wrap_fasta(record.sequence)}\n"
            )
        prepared.append(PreparedAlignment(
            1, graph_fasta,
            SourceInterval(record.contig, record.start, record.end),
            record.query_name, query_fasta, result_path, record.strand,
        ))
    return prepared


def _write_direct_sequence_fasta(
    records: Sequence[DirectSequence], path: str,
) -> int:
    """Write complete records and return their total sequence length."""
    total = 0
    with open(path, "wt") as output:
        for record in records:
            output.write(
                f">{record.query_name}\t{record.contig}:{record.start}-"
                f"{record.end}{record.strand}\tsample={record.sample}\n"
                f"{wrap_fasta(record.sequence)}\n"
            )
            total += len(record.sequence)
    return total


def alias_long_blast_query_ids(
    source: str,
    destination: str,
    maximum_length: int = BLAST_LOCAL_ID_MAXIMUM,
) -> Tuple[str, Dict[str, str]]:
    """Write a BLAST-safe FASTA and return ``(path, alias -> original)``.

    BLAST+ rejects local FASTA identifiers longer than 50 characters. Other
    aligners and every public output keep the authoritative full query name;
    only the private BLAST invocation sees deterministic ``Q1``, ``Q2``, ...
    aliases.
    """
    identifiers: List[str] = []
    with open(source, "rt") as handle:
        for line_number, raw in enumerate(handle, 1):
            if not raw.startswith(">"):
                continue
            match = re.match(r">(\S+)", raw)
            if match is None:
                raise ValueError(f"{source}:{line_number}: empty FASTA identifier")
            identifiers.append(match.group(1))
    if not identifiers:
        raise ValueError(f"{source}: FASTA contains no records")

    long_names = [
        name for name in identifiers
        if len(name.encode("utf-8")) > maximum_length
    ]
    if not long_names:
        return source, {}

    occupied = set(identifiers)
    alias_by_original: Dict[str, str] = {}
    next_alias = 1
    for original in long_names:
        if original in alias_by_original:
            continue
        while True:
            alias = f"Q{next_alias}"
            next_alias += 1
            if alias not in occupied:
                break
        occupied.add(alias)
        alias_by_original[original] = alias

    with open(source, "rt") as input_handle, open(destination, "wt") as output:
        for raw in input_handle:
            if not raw.startswith(">"):
                output.write(raw)
                continue
            ending_length = len(raw) - len(raw.rstrip("\r\n"))
            ending = raw[-ending_length:] if ending_length else ""
            content = raw[:-ending_length] if ending_length else raw
            match = re.match(r">(\S+)(.*)", content)
            if match is None:
                raise ValueError(f"{source}: empty FASTA identifier")
            original, remainder = match.groups()
            output.write(
                ">" + alias_by_original.get(original, original)
                + remainder + ending
            )

    return destination, {
        alias: original for original, alias in alias_by_original.items()
    }


def restore_blast_query_ids(
    sam_path: str, aliases: Mapping[str, str],
) -> int:
    """Restore private Q-number aliases in BLAST SAM query fields atomically."""
    if not aliases:
        return 0
    temporary = sam_path + f".names.tmp.{os.getpid()}"
    restored_rows = 0
    try:
        with open(sam_path, "rt") as source, open(temporary, "wt") as output:
            for raw in source:
                if raw.startswith("@") or not raw.strip():
                    output.write(raw)
                    continue
                ending = "\n" if raw.endswith("\n") else ""
                fields = raw.rstrip("\r\n").split("\t")
                changed = False
                for index in range(min(3, len(fields))):
                    original = aliases.get(fields[index])
                    if original is not None:
                        fields[index] = original
                        changed = True
                output.write("\t".join(fields) + ending)
                restored_rows += int(changed)
        os.replace(temporary, sam_path)
    finally:
        try:
            os.remove(temporary)
        except FileNotFoundError:
            pass
    return restored_rows


def prepare_fast_alignments(
    records: Sequence[DirectSequence], graph_fasta: str, work_root: str,
    timeout: Optional[int], threads: int, sample: str = "",
    hotspot_index: int = 1,
) -> Tuple[List[FastPreparedAlignment], str]:
    """Orient one multi-record FASTA once, then materialize cheap per-query files."""
    if threads < 1:
        raise ValueError("KmerStrd thread count must be positive")
    combined = os.path.join(work_root, "fast_input_samples.fasta")
    oriented = os.path.join(work_root, "fast_input_samples.oriented.fasta")
    _write_direct_sequence_fasta(records, combined)
    tools = resolve_extension_dir()
    run_graph_alignment_command([
        graph_alignment_tool(tools, "KmerStrd"),
        "-i", combined, "-r", graph_fasta, "-o", oriented,
        "-t", str(threads),
    ], timeout=timeout)
    if not os.path.isfile(oriented) or os.path.getsize(oriented) == 0:
        raise RuntimeError(f"KmerStrd produced no oriented FASTA: {oriented}")

    oriented_records = read_direct_fasta(oriented, sample)
    requested_names = [record.query_name for record in records]
    oriented_names = [record.query_name for record in oriented_records]
    if oriented_names != requested_names:
        raise ValueError(
            "KmerStrd changed fast-mode FASTA record names/order: "
            f"expected {requested_names[:5]!r}, observed {oriented_names[:5]!r}"
        )

    jobs_root = os.path.join(work_root, "fast_records")
    Path(jobs_root).mkdir(parents=True, exist_ok=True)
    prepared: List[FastPreparedAlignment] = []
    for index, (original, corrected) in enumerate(
        zip(records, oriented_records), 1,
    ):
        if len(original.sequence) != len(corrected.sequence):
            raise ValueError(
                f"KmerStrd changed {original.query_name!r} length from "
                f"{len(original.sequence)} to {len(corrected.sequence)}"
            )
        query_fasta = os.path.join(jobs_root, f"record_{index:08d}.fa")
        sam_path = os.path.join(jobs_root, f"record_{index:08d}.sam")
        _write_direct_sequence_fasta([corrected], query_fasta)
        prepared.append(FastPreparedAlignment(
            PreparedAlignment(
                hotspot_index, graph_fasta,
                SourceInterval(original.contig, original.start, original.end),
                original.query_name, query_fasta, sam_path, original.strand,
            ),
            corrected.sequence,
        ))
    return prepared, oriented


def iter_winnowmap_batches(
    prepared: Sequence[FastPreparedAlignment],
    maximum_bases: int = FAST_WINNOWMAP_BATCH_BASES,
) -> Iterator[List[FastPreparedAlignment]]:
    """Pack complete FASTA records without exceeding the sequence-base cap."""
    if maximum_bases < 1:
        raise ValueError("Winnowmap batch size must be positive")
    batch: List[FastPreparedAlignment] = []
    batch_bases = 0
    for item in prepared:
        size = len(item.sequence)
        if size > maximum_bases:
            raise ValueError(
                f"fast-mode query {item.prepared.query_name!r} contains "
                f"{size} bases, exceeding the Winnowmap batch limit of "
                f"{maximum_bases} bases; records are never split"
            )
        if batch and batch_bases + size > maximum_bases:
            yield batch
            batch = []
            batch_bases = 0
        batch.append(item)
        batch_bases += size
    if batch:
        yield batch


def _sam_query_name(
    line: str, query_names: Mapping[str, str], context: str,
) -> str:
    """Identify the query in BLASTN's or Winnowmap's SAM column layout."""
    fields = line.rstrip("\r\n").split("\t")
    if len(fields) < 3:
        fields = line.split()
    matches = list(dict.fromkeys(
        query_names[field] for field in fields[:3] if field in query_names
    ))
    if len(matches) == 1:
        return matches[0]
    if not matches:
        canonical_names = tuple(dict.fromkeys(query_names.values()))
        if len(canonical_names) == 1:
            # BLAST can replace a sole FASTA identifier with Query_1.  This
            # merged SAM has only one possible query, so the same fallback used
            # by normalize_graphcigarlight_sam() is unambiguous here.
            return canonical_names[0]
    expected = list(dict.fromkeys(query_names.values()))
    raise ValueError(
        f"{context}: cannot uniquely identify a query name in SAM row: "
        f"{line.rstrip()!r}; expected one of {expected[:5]!r}"
    )


FastSamSource = Tuple[str, Tuple[int, ...]]


def index_fast_sam_offsets(
    source: str, prepared: Sequence[FastPreparedAlignment],
) -> DefaultDict[str, List[int]]:
    """Scan one merged SAM once and index each query's byte offsets."""
    query_names = {
        item.prepared.query_name: item.prepared.query_name for item in prepared
    }
    # runblastn's legacy SAM layout may label queries by their one-based FASTA
    # order instead of preserving the first header token.  The combined FASTA
    # and ``prepared`` have already been checked to have identical order.
    for index, item in enumerate(prepared, 1):
        query_names.setdefault(f"Query_{index}", item.prepared.query_name)
    offsets: DefaultDict[str, List[int]] = defaultdict(list)
    with open(source, "rb") as input_handle:
        line_number = 0
        while True:
            offset = input_handle.tell()
            raw = input_handle.readline()
            if not raw:
                break
            line_number += 1
            if not raw.strip() or raw.startswith(b"@"):
                continue
            try:
                line = raw.decode("utf-8")
            except UnicodeDecodeError as error:
                raise ValueError(
                    f"{source}:{line_number}: SAM row is not UTF-8 text"
                ) from error
            query_name = _sam_query_name(
                line, query_names, f"{source}:{line_number}",
            )
            offsets[query_name].append(offset)
    return offsets


def materialize_indexed_sam(
    output_path: str, sources: Sequence[FastSamSource],
) -> int:
    """Seek indexed rows and create one query's short-lived SAM view."""
    written = 0
    with open(output_path, "wb") as output:
        for source_path, offsets in sources:
            if not offsets:
                continue
            with open(source_path, "rb") as source:
                for offset in offsets:
                    source.seek(offset)
                    row = source.readline()
                    if not row:
                        raise ValueError(
                            f"{source_path}: saved SAM offset {offset} is past EOF"
                        )
                    output.write(row)
                    written += 1
    return written


_GRAPHCIGAR_MODULES: Dict[str, Any] = {}
_GRAPHCIGAR_MODULES_LOCK = threading.Lock()
_FAST_GRAPHCIGAR_PATH = ""
_FAST_GRAPH_LENGTHS: Mapping[str, int] = {}


def _load_graphcigarlight(path: str) -> Any:
    absolute = os.path.abspath(path)
    module = _GRAPHCIGAR_MODULES.get(absolute)
    if module is not None:
        return module
    with _GRAPHCIGAR_MODULES_LOCK:
        module = _GRAPHCIGAR_MODULES.get(absolute)
        if module is not None:
            return module
        specification = importlib.util.spec_from_file_location(
            f"minsetref_fast_graphcigar_{len(_GRAPHCIGAR_MODULES)}", absolute,
        )
        if specification is None or specification.loader is None:
            raise ImportError(f"cannot import graph-CIGAR converter: {absolute}")
        module = importlib.util.module_from_spec(specification)
        specification.loader.exec_module(module)
        if not callable(getattr(module, "graphcigar", None)):
            raise AttributeError(f"{absolute}: missing callable graphcigar()")
        _GRAPHCIGAR_MODULES[absolute] = module
        return module


def convert_fast_query_alignment(
    prepared: PreparedAlignment, sequence_length: int,
    graphcigar_module: Any,
    graph_lengths: Mapping[str, int],
) -> OutputRow:
    """Convert one indexed query SAM while retaining normal-mode semantics."""
    normalize_graphcigarlight_sam(
        prepared.result_path, graph_lengths,
        {prepared.query_name: sequence_length},
    )
    raw = WholeGraphCigar(*graphcigar_module.graphcigar(
        prepared.graph_fasta, prepared.query_fasta,
        prepared.result_path, False,
    ))
    if not raw.graph_path or not raw.graph_cigar:
        result = insertion_only(prepared.query_name, sequence_length)
    else:
        result = normalize_whole_graph_cigar(
            raw, prepared.query_name, sequence_length,
        )
    strand = query_strand(prepared.query_fasta)
    if prepared.source_strand == "-":
        strand = "+" if strand == "-" else "-"
    interval = prepared.interval
    return OutputRow(
        prepared.hotspot_index, interval.contig, interval.start, interval.end,
        strand, prepared.query_name, result.graph_path, result.graph_cigar,
        result.ref_positions, result.query_positions,
    )


def _initialize_fast_conversion_worker(
    graphcigar_path: str, graph_lengths: Mapping[str, int],
) -> None:
    global _FAST_GRAPHCIGAR_PATH, _FAST_GRAPH_LENGTHS
    _FAST_GRAPHCIGAR_PATH = graphcigar_path
    _FAST_GRAPH_LENGTHS = dict(graph_lengths)


def _convert_indexed_fast_query(
    task: Tuple[PreparedAlignment, int, Tuple[FastSamSource, ...]],
) -> OutputRow:
    prepared, sequence_length, sources = task
    materialize_indexed_sam(prepared.result_path, sources)
    module = _load_graphcigarlight(_FAST_GRAPHCIGAR_PATH)
    result = convert_fast_query_alignment(
        prepared, sequence_length, module, _FAST_GRAPH_LENGTHS,
    )
    os.remove(prepared.result_path)
    return result


def _convert_indexed_fast_query_local(
    task: Tuple[PreparedAlignment, int, Tuple[FastSamSource, ...]],
    graphcigar_path: str, graph_lengths: Mapping[str, int],
) -> OutputRow:
    """Convert without process-global state (safe for graph-level threads)."""
    prepared, sequence_length, sources = task
    materialize_indexed_sam(prepared.result_path, sources)
    module = _load_graphcigarlight(graphcigar_path)
    result = convert_fast_query_alignment(
        prepared, sequence_length, module, graph_lengths,
    )
    os.remove(prepared.result_path)
    return result


def convert_indexed_fast_queries(
    prepared: Sequence[FastPreparedAlignment],
    evidence: Mapping[str, Sequence[FastSamSource]],
    graphcigar_path: str,
    graph_lengths: Mapping[str, int],
    workers: int,
    executor: Optional[concurrent.futures.ProcessPoolExecutor] = None,
) -> Dict[str, OutputRow]:
    """Convert independent indexed queries, in processes when possible."""
    tasks = [
        (
            item.prepared,
            len(item.sequence),
            tuple(evidence.get(item.prepared.query_name, ())),
        )
        for item in prepared
    ]
    if not tasks:
        return {}
    if workers <= 1:
        converted = [
            _convert_indexed_fast_query_local(
                task, graphcigar_path, graph_lengths,
            )
            for task in tasks
        ]
    else:
        if executor is None:
            raise ValueError("parallel fast conversion requires a process executor")
        converted = list(executor.map(_convert_indexed_fast_query, tasks))
    return {row.query_name: row for row in converted}


def add_fast_row_payload(row: OutputRow, sequence: str) -> OutputRow:
    result = add_query_payloads(WholeGraphCigar(
        row.query_name, row.graph_path, row.graph_cigar,
        row.ref_positions, row.query_positions,
    ), sequence)
    return dataclasses.replace(row, graph_cigar=result.graph_cigar)


def fast_row_needs_winnowmap(
    row: OutputRow, maximum_unmapped_gap: int = FAST_MAX_UNMAPPED_GAP,
) -> bool:
    """Return true when any normalized insertion is strictly over the cutoff."""
    return any(
        interval.end - interval.start > maximum_unmapped_gap
        for interval in unaligned_query_intervals(row)
    )


def fast_mask_statistics(
    prepared: Sequence[FastPreparedAlignment],
) -> Tuple[int, int]:
    """Return total lowercase bases and queries strictly over the mask cutoff."""
    masked_counts = [
        sum(base.islower() for base in item.sequence) for item in prepared
    ]
    return sum(masked_counts), sum(
        count > FAST_MASKED_BASE_THRESHOLD for count in masked_counts
    )


def run_fast_direct_alignment(
    records: Sequence[DirectSequence], graph_fasta: str,
    args: argparse.Namespace, work_root: str, output_path: str,
    hotspot_index: int = 1,
) -> Tuple[int, int]:
    """Run bulk BLASTN and bounded Winnowmap batches for one local graph.

    Fast mode sends only incompletely aligned repeat-heavy queries to
    Winnowmap.  Slow-rigorous mode uses the same batched implementation but
    sends every query to Winnowmap regardless of masking or BLASTN coverage.
    """
    slow_rigorous = bool(getattr(args, "slow_rigorous", False))
    prepared, oriented_fasta = prepare_fast_alignments(
        records, graph_fasta, work_root, args.timeout, args.threads,
        args.sample or "", hotspot_index,
    )
    tools = resolve_extension_dir()
    graphcigar_path = graph_alignment_tool(tools, "graphcigarlight.py")
    graph_lengths = fasta_sequence_lengths(graph_fasta)
    masked_bases, repeat_heavy_queries = fast_mask_statistics(prepared)
    use_masked_alignment = slow_rigorous or repeat_heavy_queries > 0
    LOG.info(
        "%s mode oriented %d records once; detected %d lowercase bases "
        "and %d queries above the %d-base mask threshold; Winnowmap "
        "eligibility is %s",
        "Slow-rigorous" if slow_rigorous else "Fast",
        len(prepared), masked_bases, repeat_heavy_queries,
        FAST_MASKED_BASE_THRESHOLD,
        "enabled" if use_masked_alignment else "disabled",
    )

    combined_sam = os.path.join(work_root, "fast_blastn.sam")
    blast_query_fasta, blast_query_aliases = alias_long_blast_query_ids(
        oriented_fasta,
        os.path.join(work_root, "fast_blast_queries.fasta"),
    )
    if blast_query_aliases:
        LOG.info(
            "BLAST local-ID compatibility: replaced %d query ID(s) longer "
            "than %d bytes with Q-number aliases; full names will be restored",
            len(blast_query_aliases), BLAST_LOCAL_ID_MAXIMUM,
        )
    blast_database = ensure_graph_blast_database(
        graph_fasta, tools, args.timeout,
    )
    repeat_options = (
        ["-dust", "yes", "-lcase_masking"] if use_masked_alignment else []
    )
    if completed_command_output(combined_sam):
        LOG.info("RESUME: reusing completed bulk BLASTN output: %s", combined_sam)
    else:
        for stale in (combined_sam, combined_sam + ".complete"):
            try:
                os.remove(stale)
            except FileNotFoundError:
                pass
        run_graph_alignment_command([
            "bash", graph_alignment_tool(tools, "runblastn"),
            "-task", "megablast", "-query", blast_query_fasta,
            "-db", blast_database, "-gapopen", "10", "-gapextend", "2",
            "-word_size", "30", "-perc_identity", str(BLAST_IDENTITY),
            *repeat_options,
            "-evalue", "1e-200", "-outfmt", "17", "-out", combined_sam,
            # Query_N is usable for single-record BLAST but is too fragile for a
            # bulk query. Preserve each FASTA ID in the SAM so byte-offset
            # routing is based on the real query name rather than an inferred
            # ordinal.
            "-parse_deflines",
            "-num_threads", str(args.threads), "-max_target_seqs", "100",
        ], timeout=args.timeout)
        restored_rows = restore_blast_query_ids(
            combined_sam, blast_query_aliases,
        )
        if blast_query_aliases:
            LOG.info(
                "Restored full query IDs in %d BLAST SAM alignment row(s)",
                restored_rows,
            )
        mark_command_output_complete(combined_sam)
    blast_offsets = index_fast_sam_offsets(combined_sam, prepared)
    blast_rows = sum(map(len, blast_offsets.values()))
    evidence: DefaultDict[str, List[FastSamSource]] = defaultdict(list)
    for item in prepared:
        name = item.prepared.query_name
        evidence[name].append((combined_sam, tuple(blast_offsets.get(name, ()))))
    LOG.info(
        "Fast BLASTN indexed %d SAM records for %d/%d queries in one scan",
        blast_rows, len(blast_offsets), len(prepared),
    )

    conversion_workers = min(
        int(getattr(args, "conversion_workers", args.threads)), len(prepared),
    )
    conversion_workers = max(1, conversion_workers)
    LOG.info(
        "Converting indexed graph alignments with %d query process(es)",
        conversion_workers,
    )
    executor: Optional[concurrent.futures.ProcessPoolExecutor] = None
    if conversion_workers > 1:
        # Fast conversion starts after the parent has indexed every query and
        # can be multi-GiB.  An explicit spawn context avoids fork-inheriting
        # that state and the process-pool stalls it can cause.
        executor = concurrent.futures.ProcessPoolExecutor(
            max_workers=conversion_workers,
            mp_context=mp.get_context("spawn"),
            initializer=_initialize_fast_conversion_worker,
            initargs=(graphcigar_path, graph_lengths),
        )
    try:
        rows = convert_indexed_fast_queries(
            prepared, evidence, graphcigar_path, graph_lengths,
            conversion_workers, executor,
        )
        if slow_rigorous:
            unfinished = list(prepared)
            LOG.info(
                "Slow-rigorous mode sends all %d queries to Winnowmap after "
                "bulk BLASTN",
                len(prepared),
            )
        else:
            unfinished = [
                item for item in prepared
                if fast_row_needs_winnowmap(rows[item.prepared.query_name])
            ]
            LOG.info(
                "Fast BLASTN completed %d/%d queries with no unmapped gap "
                "over %d bp",
                len(prepared) - len(unfinished), len(prepared),
                FAST_MAX_UNMAPPED_GAP,
            )

        if use_masked_alignment and unfinished:
            batches = list(iter_winnowmap_batches(unfinished))
            LOG.info(
                "Sending %d unfinished queries to Winnowmap in %d batch(es), "
                "each capped at %d MiB of complete records",
                len(unfinished), len(batches),
                FAST_WINNOWMAP_BATCH_BASES // (1024 * 1024),
            )
            winnow_script = graph_alignment_tool(tools, "runwinnowmaplight.sh")
            for batch_index, batch in enumerate(batches, 1):
                batch_fasta = os.path.join(
                    work_root, f"fast_winnowmap_batch_{batch_index:04d}.fa",
                )
                batch_sam = os.path.join(
                    work_root, f"fast_winnowmap_batch_{batch_index:04d}.sam",
                )
                batch_records = [
                    DirectSequence(
                        item.prepared.query_name,
                        item.prepared.interval.contig,
                        item.prepared.interval.start,
                        item.prepared.interval.end,
                        query_strand(item.prepared.query_fasta),
                        item.prepared.query_name,
                        item.sequence,
                    )
                    for item in batch
                ]
                batch_bases = _write_direct_sequence_fasta(
                    batch_records, batch_fasta,
                )
                if batch_bases > FAST_WINNOWMAP_BATCH_BASES:
                    raise AssertionError("internal Winnowmap batch limit violation")
                if completed_command_output(batch_sam):
                    LOG.info(
                        "RESUME: reusing completed Winnowmap batch %d/%d: %s",
                        batch_index, len(batches), batch_sam,
                    )
                else:
                    for stale in (batch_sam, batch_sam + ".complete"):
                        try:
                            os.remove(stale)
                        except FileNotFoundError:
                            pass
                    with open(batch_sam, "wb") as alignment_output:
                        run_graph_alignment_command([
                            "bash", winnow_script, batch_fasta, graph_fasta,
                            str(args.threads), str(BLAST_IDENTITY), "300",
                        ], stdout=alignment_output, timeout=args.timeout)
                    mark_command_output_complete(batch_sam)
                winnow_offsets = index_fast_sam_offsets(batch_sam, batch)
                added = sum(map(len, winnow_offsets.values()))
                for item in batch:
                    name = item.prepared.query_name
                    evidence[name].append((
                        batch_sam, tuple(winnow_offsets.get(name, ())),
                    ))
                LOG.info(
                    "Fast Winnowmap batch %d/%d: %d queries, %.2f MiB, "
                    "%d indexed SAM records",
                    batch_index, len(batches), len(batch),
                    batch_bases / (1024 * 1024), added,
                )
                rows.update(convert_indexed_fast_queries(
                    batch, evidence, graphcigar_path, graph_lengths,
                    min(conversion_workers, len(batch)), executor,
                ))
        elif unfinished:
            LOG.info(
                "Skipping Winnowmap because no query contains more than %d "
                "lowercase bases",
                FAST_MASKED_BASE_THRESHOLD,
            )
    finally:
        if executor is not None:
            executor.shutdown(wait=True)

    ordered_rows = [rows[item.prepared.query_name] for item in prepared]
    if args.alignment_payload:
        ordered_rows = [
            add_fast_row_payload(row, item.sequence)
            for row, item in zip(ordered_rows, prepared)
        ]
    else:
        ordered_rows = [
            dataclasses.replace(
                row, graph_cigar=strip_graph_cigar_payloads(row.graph_cigar),
            )
            for row in ordered_rows
        ]
    with open(output_path, "wt") as output:
        if args.header:
            output.write("\t".join(ALIGNMENT_HEADER_FIELDS) + "\n")
        for row in ordered_rows:
            output.write(row.line())
    return len(ordered_rows), sum(row.graph_path == "*" for row in ordered_rows)


def align_fast_hotspot_batch(
    prepared: Sequence[PreparedAlignment], args: argparse.Namespace,
    work_root: str, output_path: str,
) -> Tuple[int, int]:
    """Batch fast or slow-rigorous hotspot records by their local graph."""
    grouped: Dict[Tuple[int, str], List[PreparedAlignment]] = {}
    for item in prepared:
        grouped.setdefault(
            (item.hotspot_index, item.graph_fasta), [],
        ).append(item)
    groups = list(grouped.items())
    automatic_jobs = max(1, args.threads // args.threads_per_job)
    jobs = args.jobs if args.jobs is not None else automatic_jobs
    jobs = max(1, min(jobs, len(groups), automatic_jobs))
    slow_rigorous = bool(getattr(args, "slow_rigorous", False))
    LOG.info(
        "%s-aligning %d hotspot FASTAs against %d local graphs with %d "
        "parallel graph jobs x %d aligner threads",
        "Slow-rigorous" if slow_rigorous else "Fast",
        len(prepared), len(groups), jobs, args.threads_per_job,
    )

    def run_group(
        group: Tuple[Tuple[int, str], List[PreparedAlignment]],
    ) -> Tuple[int, str, int, int]:
        (hotspot_index, graph_fasta), items = group
        group_root = os.path.join(
            work_root, f"fast_hotspot_{hotspot_index:08d}",
        )
        Path(group_root).mkdir(parents=True, exist_ok=True)
        group_output = os.path.join(group_root, "alignment.tsv")
        records = [
            DirectSequence(
                args.sample or "", item.interval.contig,
                item.interval.start, item.interval.end, item.source_strand,
                item.query_name, read_single_fasta_sequence(item.query_fasta),
            )
            for item in items
        ]
        group_args = argparse.Namespace(**vars(args))
        group_args.threads = args.threads_per_job
        group_args.header = False
        group_args.slow_rigorous = slow_rigorous
        # Graphs are parallelized here. Avoid nesting a process pool inside
        # every graph job; conversion remains parallel across graph jobs.
        group_args.conversion_workers = 1
        output_rows, insertions = run_fast_direct_alignment(
            records, graph_fasta, group_args, group_root, group_output,
            hotspot_index=hotspot_index,
        )
        return hotspot_index, group_output, output_rows, insertions

    results: Dict[int, Tuple[str, int, int]] = {}
    if jobs == 1:
        completed_groups = map(run_group, groups)
        executor = None
    else:
        executor = concurrent.futures.ThreadPoolExecutor(max_workers=jobs)
        completed_groups = executor.map(run_group, groups)
    try:
        for completed, result in enumerate(completed_groups, 1):
            hotspot_index, group_output, output_rows, insertions = result
            results[hotspot_index] = (group_output, output_rows, insertions)
            if completed == len(groups) or completed % max(1, len(groups) // 20) == 0:
                LOG.info(
                    "%s local-graph alignment progress: %d/%d graphs",
                    "Slow-rigorous" if slow_rigorous else "Fast",
                    completed, len(groups),
                )
    finally:
        if executor is not None:
            executor.shutdown(wait=True)

    output_rows = 0
    insertion_rows = 0
    with open(output_path, "at") as output:
        for (hotspot_index, _graph_fasta), _items in groups:
            group_output, rows, insertions = results[hotspot_index]
            with open(group_output, "rt") as source:
                shutil.copyfileobj(source, output)
            output_rows += rows
            insertion_rows += insertions
    return output_rows, insertion_rows


def direct_alignment_signature(
    records: Sequence[DirectSequence], graph_fasta: str,
    args: argparse.Namespace,
) -> str:
    """Hash all inputs that make retained direct-alignment work reusable."""
    digest = hashlib.sha256()
    digest.update((DIRECT_WORK_PROTOCOL + "\n").encode("utf-8"))
    slow_rigorous = bool(getattr(args, "slow_rigorous", False))
    # Retain the exact signature prefix used by existing fast/legacy work
    # directories so old timed-out commands can still resume their completed
    # alignment batches. Slow-rigorous adds a new discriminator.
    digest.update(
        f"fast={int(bool(getattr(args, 'fast_mode', False)) or slow_rigorous)}\n"
        .encode("utf-8")
    )
    if slow_rigorous:
        digest.update(b"slow-rigorous=1\n")
    with open(graph_fasta, "rb") as graph_handle:
        while True:
            chunk = graph_handle.read(1024 * 1024)
            if not chunk:
                break
            digest.update(chunk)
    for record in records:
        for value in (
            record.sample, record.contig, str(record.start), str(record.end),
            record.strand, record.query_name, record.record_class,
        ):
            encoded = value.encode("utf-8")
            digest.update(len(encoded).to_bytes(8, "little"))
            digest.update(encoded)
        sequence = record.sequence.encode("ascii")
        digest.update(len(sequence).to_bytes(8, "little"))
        digest.update(sequence)
    return digest.hexdigest() + "\n"


def prepare_direct_work_directory(
    path: str, signature: str,
) -> bool:
    """Prepare a stable work directory; return true when it is resumable."""
    signature_path = os.path.join(path, ".resume-signature")
    resumable = False
    if os.path.isdir(path):
        try:
            resumable = (
                Path(signature_path).read_text(encoding="utf-8") == signature
            )
        except FileNotFoundError:
            resumable = False
        if not resumable:
            shutil.rmtree(path)
    elif os.path.exists(path):
        os.remove(path)
    Path(path).mkdir(parents=True, exist_ok=True)
    if not resumable:
        atomic_write_text(signature_path, signature)
    return resumable


def align_direct(args: argparse.Namespace) -> int:
    if args.threads < 1 or args.threads_per_job < 1:
        raise ValueError("thread counts must be positive")
    graph_fasta = os.path.abspath(args.graph)
    if not os.path.isfile(graph_fasta):
        raise FileNotFoundError(graph_fasta)
    input_format = args.input_format
    if input_format == "auto":
        lowered = args.input.lower()
        if lowered.endswith(".bed"):
            input_format = "bed"
        elif lowered.endswith(tuple(suffix.lower() for suffix in FASTA_SUFFIXES)):
            input_format = "fasta"
        else:
            raise ValueError(
                "direct --graph mode needs --input-format bed or fasta for "
                f"unrecognized suffix: {args.input}"
            )
    if input_format == "bed":
        if not args.query:
            raise ValueError("direct BED alignment requires -q/--query query-path list")
        records = read_direct_bed(args.input, args.query)
    elif input_format == "fasta":
        records = read_direct_fasta(args.input, args.sample or "")
    else:
        raise ValueError("--graph direct mode accepts only BED or FASTA input")

    output_path = os.path.abspath(args.output)
    output_parent = os.path.dirname(output_path) or "."
    Path(output_parent).mkdir(parents=True, exist_ok=True)
    stable_prefix = f".{os.path.basename(output_path)}.direct.work"
    work_root = os.path.join(output_parent, stable_prefix)
    work_lock = work_root + ".lock"
    temporary_output = output_path + f".tmp.{os.getpid()}"
    signature = direct_alignment_signature(records, graph_fasta, args)
    with open(work_lock, "a+") as lock_handle:
        fcntl.flock(lock_handle.fileno(), fcntl.LOCK_EX)
        try:
            # Random work directories from older releases cannot be verified
            # or resumed. Remove only this output's precisely scoped legacy
            # directories before installing the stable checkpoint directory.
            for legacy in Path(output_parent).glob(stable_prefix + ".*"):
                if legacy.is_dir():
                    shutil.rmtree(legacy)
            resumed = prepare_direct_work_directory(work_root, signature)
            if resumed:
                LOG.info(
                    "RESUME: reusing direct-alignment checkpoints: %s",
                    work_root,
                )
            try:
                if (
                    getattr(args, "fast_mode", False)
                    or getattr(args, "slow_rigorous", False)
                ):
                    output_rows, insertions = run_fast_direct_alignment(
                        records, graph_fasta, args, work_root, temporary_output,
                    )
                else:
                    prepared = prepare_direct_alignments(
                        records, graph_fasta, work_root,
                    )
                    with open(temporary_output, "wt") as output:
                        if args.header:
                            output.write(
                                "\t".join(ALIGNMENT_HEADER_FIELDS) + "\n"
                            )
                    output_rows, insertions = align_prepared_batch(
                        prepared, args, "sample", temporary_output,
                    )
                os.replace(temporary_output, output_path)
                atomic_write_text(
                    output_path + ".mode", alignment_mode_protocol(args),
                )
                LOG.info(
                    "Aligned %d sample records to one graph; wrote %d rows "
                    "(%d insertion-only) to %s",
                    len(records), output_rows, insertions, output_path,
                )
            except subprocess.TimeoutExpired:
                try:
                    os.remove(temporary_output)
                except FileNotFoundError:
                    pass
                LOG.warning(
                    "Retained timed-out direct-alignment checkpoints for "
                    "resume: %s", work_root,
                )
                raise
            except Exception:
                try:
                    os.remove(temporary_output)
                except FileNotFoundError:
                    pass
                LOG.error(
                    "Retained failed direct-alignment work directory: %s",
                    work_root,
                )
                raise
            else:
                if args.keep_work:
                    LOG.info(
                        "Retained direct-alignment work directory: %s",
                        work_root,
                    )
                else:
                    shutil.rmtree(work_root)
        finally:
            fcntl.flock(lock_handle.fileno(), fcntl.LOCK_UN)
    return 0


def align_all(args: argparse.Namespace) -> int:
    if args.threads < 1 or args.threads_per_job < 1:
        raise ValueError("thread counts must be positive")
    if args.hotspot_kmer_length < 1:
        raise ValueError("hotspot k-mer length must be positive")
    if args.extension < 0 or args.merge_gap < 0:
        raise ValueError("extension and merge gap must be nonnegative")
    graph_fastas = read_graph_list(args.graph_list, args.graph_folder)
    output_parent = os.path.dirname(os.path.abspath(args.output)) or "."
    Path(output_parent).mkdir(parents=True, exist_ok=True)
    work_root = tempfile.mkdtemp(
        prefix=f".{os.path.basename(args.output)}.work.", dir=output_parent,
    )
    temporary_output = args.output + f".tmp.{os.getpid()}"

    try:
        with IndexedFasta(args.query) as fasta:
            tasks, input_rows = build_partition_tasks(
                args.input, graph_fastas, fasta,
                args.extension, args.merge_gap, args.hotspot_kmer_length,
            )
            if not tasks:
                raise ValueError(
                    "no hotspot rows with an available local graph found in "
                    f"{args.input}"
                )
            contigs = load_hotspot_contigs(fasta, tasks)
        loaded_bases = sum(map(len, contigs.values()))
        LOG.info(
            "Loaded %d hotspot rows and %d required assembly contigs "
            "(%.2f GiB) into RAM",
            input_rows, len(contigs), loaded_bases / (1024 ** 3),
        )
        prepared = prepare_alignments(
            tasks, contigs, args.sample, work_root,
        )
        LOG.info(
            "Materialized %d hotspot FASTAs for %d local graphs before "
            "starting alignment",
            len(prepared), len(tasks),
        )
        del contigs
        with open(temporary_output, "wt") as output:
            if args.header:
                output.write(
                    "hotspot_index\tcontig\tstart\tend\tstrand\tqueryname\t"
                    "graph_path\tgraphcigar\trefpositions\tqpositions\n"
                )
        if (
            getattr(args, "fast_mode", False)
            or getattr(args, "slow_rigorous", False)
        ):
            output_rows, insertions = align_fast_hotspot_batch(
                prepared, args, work_root, temporary_output,
            )
        else:
            output_rows, insertions = align_prepared_batch(
                prepared, args, "initial", temporary_output,
            )
        os.replace(temporary_output, args.output)
        atomic_write_text(
            args.output + ".mode", alignment_mode_protocol(args),
        )
        LOG.info(
            "Wrote %d whole-sequence graph CIGAR rows (%d insertion-only) to %s",
            output_rows, insertions, args.output,
        )
    except Exception:
        removed = []
        try:
            if os.path.exists(temporary_output):
                os.remove(temporary_output)
                removed.append(temporary_output)
        except OSError as error:
            LOG.warning(
                "Could not remove failed temporary alignment output %s: %s",
                temporary_output, error,
            )
        try:
            if os.path.isdir(work_root):
                shutil.rmtree(work_root)
                removed.append(work_root)
        except OSError as error:
            LOG.warning(
                "Could not remove failed-job work directory %s: %s",
                work_root, error,
            )
        if removed:
            LOG.error("Removed failed-job artifact(s): %s", ", ".join(removed))
        raise
    else:
        if args.keep_work:
            LOG.info("Retained work directory: %s", work_root)
        else:
            shutil.rmtree(work_root)
    return 0


def parse_args(argv: Optional[Sequence[str]] = None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Align KmerSearcher hotspots to ordered local graph FASTAs",
    )
    parser.add_argument(
        "-i", "--input", required=True,
        help="_hotspot.txt, BED regions, or multi-FASTA",
    )
    parser.add_argument(
        "-q", "--query",
        help="assembly FASTA in hotspot mode; query-path list in direct BED mode",
    )
    parser.add_argument(
        "--input-format", choices=("auto", "hotspot", "bed", "fasta"),
        default="auto",
        help="direct-input format; hotspot mode remains selected when --graph is absent",
    )
    parser.add_argument(
        "--graph",
        help="one graph .FA for direct BED/FASTA alignment",
    )
    parser.add_argument(
        "-G", "--graph-folder",
        help="root containing GRAPH_NAME/GRAPH_NAME.FA local graphs",
    )
    parser.add_argument(
        "-L", "--graph-list",
        help="ordered graph-name/path list; hotspot indices are 1-based",
    )
    parser.add_argument("-s", "--sample", help="sample/haplotype name")
    parser.add_argument("-o", "--output", required=True, help="whole-CIGAR output TSV")
    parser.add_argument("-t", "--threads", type=int, default=4, help="total CPU budget")
    parser.add_argument(
        "-j", "--jobs", type=int, help="parallel hotspot-FASTA alignment jobs",
    )
    parser.add_argument("--threads-per-job", type=int, default=2)
    parser.add_argument(
        "--extension", type=int, default=DEFAULT_EXTENSION,
        help="anchor added to each side after merging (default: 30000)",
    )
    parser.add_argument(
        "--merge-gap", type=int, default=DEFAULT_MERGE_GAP,
        help=(
            "raw-hotspot merge distance before anchor extension "
            "(default: 60000, twice the default anchor)"
        ),
    )
    parser.add_argument(
        "--hotspot-kmer-length", type=int,
        default=DEFAULT_HOTSPOT_KMER_LENGTH,
        help="KmerSearcher k-mer length used to recover raw four-column spans",
    )
    parser.add_argument("--timeout", type=int, default=3600, help="timeout per external command")
    parser.add_argument(
        "--fast", "--fast-mode", dest="fast_mode", action="store_true",
        help=(
            "orient the multi-record query once, run one BLASTN pass, and "
            "send only unfinished queries to Winnowmap in record-safe batches "
            "capped at 256 MiB; --fast-mode is a compatibility alias"
        ),
    )
    parser.add_argument(
        "--rigorous", "--slow-rigorous", dest="slow_rigorous", action="store_true",
        help=(
            "in direct --graph mode, use the same batch engine but run both "
            "BLASTN and Winnowmap for every query"
        ),
    )
    parser.add_argument("--header", action="store_true", help="include a TSV header")
    payload_group = parser.add_mutually_exclusive_group()
    payload_group.add_argument(
        "--alignment-payload", dest="alignment_payload", action="store_true",
        default=False,
        help=(
            "attach oriented query bases to every I and X graph-CIGAR "
            "operation; by default graphcigartoref restores bases per comparison"
        ),
    )
    payload_group.add_argument(
        "--no-alignment-payload", dest="alignment_payload", action="store_false",
        help="save compact graph CIGARs without I/X sequence payloads (default)",
    )
    parser.add_argument(
        "--cycle",
        action="store_true",
        default=False,
        help=(
            "opt in to re-align residual unaligned chunks longer than 300 bp using "
            "Winnowmap only until none remain or a cycle finds no "
            "additional alignments; disabled by default"
        ),
    )
    parser.add_argument("--keep-work", action="store_true")
    parser.add_argument("-v", "--verbose", action="store_true")
    args = parser.parse_args(argv)
    if args.timeout < 1:
        parser.error("--timeout must be positive")
    selected_modes = sum((
        bool(args.fast_mode), bool(args.slow_rigorous), bool(args.cycle),
    ))
    if selected_modes > 1:
        parser.error(
            "--fast, --slow-rigorous, and --cycle are mutually exclusive"
        )
    if args.graph:
        if args.input_format == "hotspot":
            parser.error("--graph direct mode cannot use --input-format hotspot")
    else:
        missing = [
            option for option, value in (
                ("-q/--query", args.query),
                ("-G/--graph-folder", args.graph_folder),
                ("-L/--graph-list", args.graph_list),
                ("-s/--sample", args.sample),
            ) if not value
        ]
        if missing:
            parser.error("hotspot mode requires " + ", ".join(missing))
        if args.input_format not in {"auto", "hotspot"}:
            parser.error("BED/FASTA --input-format requires --graph")
    return args


def main(argv: Optional[Sequence[str]] = None) -> int:
    args = parse_args(argv)
    logging.basicConfig(
        level=logging.DEBUG if args.verbose else logging.INFO,
        format="[%(asctime)s] %(levelname)s: %(message)s",
    )
    try:
        return align_direct(args) if args.graph else align_all(args)
    except subprocess.TimeoutExpired as error:
        LOG.warning(
            "Alignment command timed out after %s seconds; checkpoints were "
            "retained for resume",
            error.timeout,
        )
        return TIMEOUT_EXIT_CODE
    except Exception as error:
        LOG.error("%s", error)
        if args.verbose:
            LOG.exception("graph alignment failed")
        return 1


if __name__ == "__main__":
    raise SystemExit(main())
