#!/usr/bin/env python3
"""Annotate the combined global-novel pool against the cleaned reference.

The input pool contains both stage-one global-novel intervals and residual
pieces left by local-novel mapping. Coordinates are resolved from partition
``.header`` files and sequence is fetched from the original assembly
FASTA/FAI paths in ``query_paths.txt``. All queries share one combined FASTA
and pass through ordered Minimap2, Winnowmap, and BLASTN residual stages.
Accepted blocks are converted back to the original uncleaned-reference
coordinates and encoded as graphic paths. Query sequence is never swapped.
"""
from __future__ import annotations

import argparse
import csv
import dataclasses
import hashlib
import json
import logging
import multiprocessing as mp
import os
import re
import shutil
import sys
from collections import defaultdict
from collections.abc import Mapping as MappingABC
from concurrent.futures import FIRST_COMPLETED, ProcessPoolExecutor, wait
from pathlib import Path
from typing import Dict, Iterator, List, Optional, Sequence, Tuple

from map_partition_local_novel import (
    AssemblyFastaSource,
    CachedCandidateSequence,
    MappingCandidate,
    PLACEMENT_RESOLUTION_VERSION,
    REPORTING_RESOLUTION_VERSION,
    load_candidate_header_records,
    read_assembly_sources,
    read_candidates,
    read_partition_inputs,
    read_placements,
    read_reporting_components,
    write_reporting_components,
)
from refine_partition_paths import (
    GlobalPlacement,
    NovelQuery,
    _fingerprint,
    align_global_queries,
    choose_local_placements,
    filter_queries_by_score,
    filter_queries_for_global_lift,
    global_mapping_candidates,
    global_queries_from_bed,
    promoted_local_queries,
    read_novel_queries,
    write_global_alignment_summary,
    write_global_mapping_outputs,
)
from minsetref_core import IndexedFasta
from minsetref_graphic import (
    build_graphic_path,
    build_multilocus_graphic_path,
)
from minsetref_multilocus import ResolvedPlacementComponent
from summarize_partition_novelty import HeaderRecord
from fixed_alternatives import (
    ANNOTATION_POLICY, annotation_reference_records, annotation_reference_sources,
)


LOG = logging.getLogger("map_global_novel")
_COMPLEMENT = str.maketrans(
    "ACGTRYKMSWBDHVNacgtrykmswbdhvn",
    "TGCAYRMKSWVHDBNtgcayrmkswvhdbn",
)


@dataclasses.dataclass(frozen=True)
class ExtractionInterval:
    candidate_id: str
    source_contig: str
    source_start: int
    source_end: int
    local_start: int
    local_end: int


@dataclasses.dataclass(frozen=True)
class ExtractionTask:
    task_index: int
    haplotype: str
    source: AssemblyFastaSource
    intervals: Tuple[ExtractionInterval, ...]
    shard_fasta: str
    shard_index: str
    shard_marker: str
    signature: str
    resume: bool


def _reverse_complement(sequence: str) -> str:
    return sequence.translate(_COMPLEMENT)[::-1]


def _append_cigar_op(
    operations: List[Tuple[int, str]],
    length: int,
    operation: str,
) -> None:
    if length <= 0:
        return
    if operations and operations[-1][1] == operation:
        prior_length, _prior_operation = operations[-1]
        operations[-1] = (prior_length + length, operation)
    else:
        operations.append((length, operation))


def _exact_base_operations(
    query_sequence: str,
    reference_sequence: str,
) -> List[Tuple[int, str]]:
    if len(query_sequence) != len(reference_sequence):
        raise ValueError(
            "paired alignment block has unequal query/reference sequence "
            f"lengths: {len(query_sequence)} != {len(reference_sequence)}"
        )
    operations: List[Tuple[int, str]] = []
    for query_base, reference_base in zip(
        query_sequence, reference_sequence,
    ):
        _append_cigar_op(
            operations,
            1,
            "="
            if query_base.upper() == reference_base.upper()
            else "X",
        )
    return operations


def build_graphic_liftover_path(
    query: NovelQuery,
    placement: GlobalPlacement,
    cached: CachedCandidateSequence,
    reference: IndexedFasta,
) -> str:
    """Encode current accepted blocks on the original uncleaned reference."""
    query_sequence = cached.sequence[
        query.local_start - cached.local_start:
        query.local_end - cached.local_start
    ]
    if len(query_sequence) != query.local_end - query.local_start:
        raise ValueError(
            f"{query.query_id}: incomplete cached sequence for graphic "
            "liftover"
        )

    return build_graphic_path(
        query_id=query.query_id,
        query_start=query.local_start,
        query_end=query.local_end,
        query_sequence=query_sequence,
        ref_contig=placement.ref_contig,
        strand=placement.strand,
        paired_blocks=placement.paired_blocks,
        reference=reference,
    )


def annotate_global_placements(
    queries: Sequence[NovelQuery],
    placements: Dict[str, GlobalPlacement],
    cached_queries: MappingABC,
    assembly_sources: Dict[str, AssemblyFastaSource],
) -> Dict[str, GlobalPlacement]:
    """Attach original-reference graphic paths without changing query bases."""
    query_by_id = {query.query_id: query for query in queries}
    readers: Dict[str, IndexedFasta] = {}
    output: Dict[str, GlobalPlacement] = {}
    try:
        for query_id, placement in placements.items():
            query = query_by_id[query_id]
            source = assembly_sources.get(placement.ref_haplotype)
            if source is None:
                raise ValueError(
                    "--query-paths lacks the original uncleaned reference "
                    f"haplotype {placement.ref_haplotype!r}, required to "
                    f"annotate {query_id}"
                )
            reader = readers.get(placement.ref_haplotype)
            if reader is None:
                reader = IndexedFasta(
                    source.fasta_path, source.fai_path,
                )
                readers[placement.ref_haplotype] = reader
            path = build_graphic_liftover_path(
                query,
                placement,
                cached_queries[query_id],
                reader,
            )
            output[query_id] = dataclasses.replace(
                placement, liftover_path=path,
            )
    finally:
        for reader in readers.values():
            reader.close()
    return output


def annotate_global_reporting_components(
    queries: Sequence[NovelQuery],
    components: Sequence[ResolvedPlacementComponent],
    cached_queries: MappingABC,
    assembly_sources: Dict[str, AssemblyFastaSource],
) -> List[ResolvedPlacementComponent]:
    """Attach one complete ordered multi-locus path to each query's rows."""
    query_by_id = {query.query_id: query for query in queries}
    grouped: Dict[str, List[ResolvedPlacementComponent]] = defaultdict(list)
    for component in components:
        grouped[component.query_id].append(component)
    readers: Dict[str, IndexedFasta] = {}
    output: List[ResolvedPlacementComponent] = []
    try:
        for query_id in sorted(grouped):
            query = query_by_id.get(query_id)
            if query is None:
                raise ValueError(
                    f"multi-locus report references unknown query {query_id!r}"
                )
            members = sorted(
                grouped[query_id],
                key=lambda row: (
                    row.query_start, row.query_end, row.component_index,
                ),
            )
            reference_map = {}
            for member in members:
                source = assembly_sources.get(member.ref_haplotype)
                if source is None:
                    raise ValueError(
                        "--query-paths lacks original reference haplotype "
                        f"{member.ref_haplotype!r} required by {query_id}"
                    )
                reader = readers.get(member.ref_haplotype)
                if reader is None:
                    reader = IndexedFasta(
                        source.fasta_path, source.fai_path,
                    )
                    readers[member.ref_haplotype] = reader
                reference_map[(member.ref_haplotype, member.ref_contig)] = reader
            cached = cached_queries[query_id]
            query_sequence = cached.sequence[
                query.local_start - cached.local_start:
                query.local_end - cached.local_start
            ]
            path = build_multilocus_graphic_path(
                query_id=query_id,
                query_start=query.local_start,
                query_end=query.local_end,
                query_sequence=query_sequence,
                components=members,
                references=reference_map,
            )
            output.extend(
                dataclasses.replace(member, liftover_path=path)
                for member in members
            )
    finally:
        for reader in readers.values():
            reader.close()
    return output


@dataclasses.dataclass(frozen=True)
class ExtractionResult:
    task_index: int
    haplotype: str
    shard_fasta: str
    shard_index: str
    records: int
    bases: int
    reused: bool


class DiskCandidateSequenceStore(
    MappingABC[str, CachedCandidateSequence]
):
    """Random-access candidate windows without retaining sequence strings."""

    retain_cycle_sequences = False

    def __init__(self, fasta_path: str, index_path: str):
        self.fasta_path = fasta_path
        self.index_path = index_path
        self._index: Dict[str, Tuple[int, int, int, int]] = {}
        with open(index_path, "rt", newline="") as handle:
            reader = csv.DictReader(handle, delimiter="\t")
            required = {
                "candidate_id", "local_start", "local_end",
                "sequence_offset", "sequence_length",
            }
            if reader.fieldnames is None or not required.issubset(
                reader.fieldnames
            ):
                raise ValueError(
                    f"{index_path}: expected columns {sorted(required)}"
                )
            for row in reader:
                candidate_id = row["candidate_id"]
                if candidate_id in self._index:
                    raise ValueError(
                        f"{index_path}: duplicate candidate {candidate_id!r}"
                    )
                self._index[candidate_id] = (
                    int(row["local_start"]),
                    int(row["local_end"]),
                    int(row["sequence_offset"]),
                    int(row["sequence_length"]),
                )
        self._fd = os.open(fasta_path, os.O_RDONLY)

    def __getitem__(self, candidate_id: str) -> CachedCandidateSequence:
        local_start, local_end, offset, length = self._index[candidate_id]
        data = os.pread(self._fd, length, offset)
        if len(data) != length:
            raise ValueError(
                f"{self.fasta_path}: incomplete cached sequence for "
                f"{candidate_id!r}"
            )
        return CachedCandidateSequence(
            local_start, local_end, data.decode("ascii"),
        )

    def __iter__(self) -> Iterator[str]:
        return iter(self._index)

    def __len__(self) -> int:
        return len(self._index)

    def close(self) -> None:
        fd = getattr(self, "_fd", -1)
        if fd >= 0:
            os.close(fd)
            self._fd = -1


def _extract_assembly_shard(task: ExtractionTask) -> ExtractionResult:
    if task.resume and os.path.isfile(task.shard_marker):
        try:
            with open(task.shard_marker, "rt") as handle:
                saved = json.load(handle)
            if (
                saved.get("signature") == task.signature
                and os.path.isfile(task.shard_fasta)
                and os.path.isfile(task.shard_index)
            ):
                return ExtractionResult(
                    task.task_index,
                    task.haplotype,
                    task.shard_fasta,
                    task.shard_index,
                    int(saved["records"]),
                    int(saved["bases"]),
                    True,
                )
        except (OSError, ValueError, KeyError, TypeError):
            pass

    temporary_fasta = task.shard_fasta + ".tmp"
    temporary_index = task.shard_index + ".tmp"
    for path in (temporary_fasta, temporary_index):
        try:
            os.remove(path)
        except FileNotFoundError:
            pass

    total_bases = 0
    with IndexedFasta(
        task.source.fasta_path, task.source.fai_path,
    ) as source, open(temporary_fasta, "wb") as fasta_out, open(
        temporary_index, "wt", newline="",
    ) as index_out:
        writer = csv.writer(
            index_out, delimiter="\t", lineterminator="\n",
        )
        writer.writerow((
            "candidate_id", "local_start", "local_end",
            "sequence_offset", "sequence_length",
        ))
        for interval in task.intervals:
            if interval.source_contig not in source.index:
                raise ValueError(
                    f"{task.haplotype}: contig "
                    f"{interval.source_contig!r} is absent from "
                    f"{task.source.fasta_path}"
                )
            contig_length = source.length(interval.source_contig)
            if (
                interval.source_start < 0
                or interval.source_end <= interval.source_start
                or interval.source_end > contig_length
            ):
                raise ValueError(
                    f"{interval.candidate_id}: extraction "
                    f"{interval.source_contig}:{interval.source_start}-"
                    f"{interval.source_end} is outside assembly contig length "
                    f"{contig_length}"
                )
            # Use the supplied .fai for direct random access. Do not read the
            # complete ~3-Gb assembly when this task needs only a small set of
            # query windows.
            sequence = source.fetch(
                interval.source_contig,
                interval.source_start,
                interval.source_end,
            )
            fasta_out.write(f">{interval.candidate_id}\n".encode())
            offset = fasta_out.tell()
            encoded = sequence.encode("ascii")
            fasta_out.write(encoded)
            fasta_out.write(b"\n")
            writer.writerow((
                interval.candidate_id,
                interval.local_start,
                interval.local_end,
                offset,
                len(encoded),
            ))
            total_bases += len(encoded)
    os.replace(temporary_fasta, task.shard_fasta)
    os.replace(temporary_index, task.shard_index)
    temporary_marker = task.shard_marker + ".tmp"
    with open(temporary_marker, "wt") as out:
        json.dump({
            "signature": task.signature,
            "records": len(task.intervals),
            "bases": total_bases,
        }, out, sort_keys=True)
        out.write("\n")
    os.replace(temporary_marker, task.shard_marker)
    return ExtractionResult(
        task.task_index,
        task.haplotype,
        task.shard_fasta,
        task.shard_index,
        len(task.intervals),
        total_bases,
        False,
    )


def _bounded_process_results(
    tasks: Sequence[ExtractionTask],
    jobs: int,
) -> Iterator[ExtractionResult]:
    iterator = iter(tasks)
    with ProcessPoolExecutor(
        max_workers=jobs, mp_context=mp.get_context("spawn"),
    ) as executor:
        pending = set()
        for _index in range(min(jobs, len(tasks))):
            pending.add(executor.submit(
                _extract_assembly_shard, next(iterator),
            ))
        while pending:
            completed, pending = wait(
                pending, return_when=FIRST_COMPLETED,
            )
            for future in completed:
                yield future.result()
                try:
                    task = next(iterator)
                except StopIteration:
                    continue
                pending.add(executor.submit(
                    _extract_assembly_shard, task,
                ))


def build_disk_query_cache(
    candidates: Sequence[MappingCandidate],
    source_records: Dict[Tuple[str, str], HeaderRecord],
    assembly_sources: Dict[str, AssemblyFastaSource],
    anchor: int,
    jobs: int,
    workdir: str,
    signature: str,
    resume: bool,
    compatible_signatures: Sequence[str] = (),
) -> DiskCandidateSequenceStore:
    """Extract assembly-parallel shards and combine them without sequence RAM."""
    Path(workdir).mkdir(parents=True, exist_ok=True)
    combined_fasta = os.path.join(workdir, "candidate_sequences.fa")
    combined_index = os.path.join(workdir, "candidate_sequences.index.tsv")
    marker = os.path.join(workdir, "_SUCCESS.json")
    if resume and os.path.isfile(marker):
        try:
            with open(marker, "rt") as handle:
                saved = json.load(handle)
            if (
                saved.get("signature") in {
                    signature, *compatible_signatures,
                }
                and os.path.isfile(combined_fasta)
                and os.path.isfile(combined_index)
            ):
                LOG.info(
                    "RESUME: reusing disk-backed query cache (%s)",
                    combined_fasta,
                )
                return DiskCandidateSequenceStore(
                    combined_fasta, combined_index,
                )
        except (OSError, ValueError, AttributeError):
            pass

    shard_dir = os.path.join(workdir, "shards")
    Path(shard_dir).mkdir(parents=True, exist_ok=True)
    by_haplotype: Dict[str, List[ExtractionInterval]] = defaultdict(list)
    for candidate in candidates:
        record = source_records[(
            candidate.novel_partition, candidate.record_id,
        )]
        local_start = max(0, candidate.local_start - anchor)
        local_end = min(record.record_length, candidate.local_end + anchor)
        by_haplotype[candidate.haplotype].append(ExtractionInterval(
            candidate.candidate_id,
            record.source_contig,
            record.source_start + local_start,
            record.source_start + local_end,
            local_start,
            local_end,
        ))

    tasks: List[ExtractionTask] = []
    for task_index, haplotype in enumerate(sorted(by_haplotype)):
        source = assembly_sources[haplotype]
        prefix = os.path.join(shard_dir, f"{task_index:06d}")
        task_signature = hashlib.sha256(json.dumps({
            "signature": signature,
            "haplotype": haplotype,
            "candidate_ids": [
                row.candidate_id for row in by_haplotype[haplotype]
            ],
        }, sort_keys=True).encode()).hexdigest()
        tasks.append(ExtractionTask(
            task_index,
            haplotype,
            source,
            tuple(by_haplotype[haplotype]),
            prefix + ".fa",
            prefix + ".index.tsv",
            prefix + "._SUCCESS.json",
            task_signature,
            resume,
        ))
    LOG.info(
        "Extracting %d query windows from %d assemblies with %d processes "
        "using direct .fai interval fetches (complete assemblies are not "
        "loaded)",
        len(candidates), len(tasks), min(jobs, max(1, len(tasks))),
    )
    results: List[ExtractionResult] = []
    progress_step = max(1, len(tasks) // 20)
    for completed, result in enumerate(
        _bounded_process_results(
            tasks, min(jobs, max(1, len(tasks))),
        ),
        1,
    ):
        results.append(result)
        if completed == len(tasks) or completed % progress_step == 0:
            LOG.info(
                "Assembly query extraction: %d/%d complete",
                completed, len(tasks),
            )

    temporary_fasta = combined_fasta + ".tmp"
    temporary_index = combined_index + ".tmp"
    total_records = 0
    total_bases = 0
    with open(temporary_fasta, "wb") as fasta_out, open(
        temporary_index, "wt", newline="",
    ) as index_out:
        writer = csv.writer(
            index_out, delimiter="\t", lineterminator="\n",
        )
        writer.writerow((
            "candidate_id", "local_start", "local_end",
            "sequence_offset", "sequence_length",
        ))
        for result in sorted(results, key=lambda row: row.task_index):
            base_offset = fasta_out.tell()
            with open(result.shard_fasta, "rb") as shard:
                shutil.copyfileobj(shard, fasta_out, length=16 * 1024 * 1024)
            with open(result.shard_index, "rt", newline="") as handle:
                for row in csv.DictReader(handle, delimiter="\t"):
                    writer.writerow((
                        row["candidate_id"],
                        row["local_start"],
                        row["local_end"],
                        base_offset + int(row["sequence_offset"]),
                        row["sequence_length"],
                    ))
            total_records += result.records
            total_bases += result.bases
    os.replace(temporary_fasta, combined_fasta)
    os.replace(temporary_index, combined_index)
    temporary_marker = marker + ".tmp"
    with open(temporary_marker, "wt") as out:
        json.dump({
            "signature": signature,
            "records": total_records,
            "bases": total_bases,
        }, out, sort_keys=True)
        out.write("\n")
    os.replace(temporary_marker, marker)
    shutil.rmtree(shard_dir)
    LOG.info(
        "Combined %d query windows (%.2f GiB) into disk-backed cache",
        total_records, total_bases / (1024 ** 3),
    )
    return DiskCandidateSequenceStore(combined_fasta, combined_index)


def parse_args(argv: Optional[Sequence[str]] = None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=(
            "Combine original global novelty with local residuals and map the "
            "pool to the complete cleaned reference."
        )
    )
    parser.add_argument("-a", "--analysis-dir", required=True)
    parser.add_argument(
        "-l", "--local-mapping-dir", required=True,
        help="output of map_partition_local_novel.py",
    )
    parser.add_argument(
        "-r", "--main-reference", required=True,
        help="cleaned cohort main_chroms.fa",
    )
    parser.add_argument(
        "-q", "--query-paths", required=True,
        help="original HAPLOTYPE FASTA (with adjacent FASTA.fai) list used to fetch query bases",
    )
    parser.add_argument("-o", "--output-dir", required=True)
    parser.add_argument(
        "-j", "--jobs", type=int, default=16,
        help="parallel original-assembly query loaders (default: 16)",
    )
    parser.add_argument(
        "-c", "--cores", type=int, default=16,
        help="threads for the single combined global alignment (default: 16)",
    )
    parser.add_argument("--alignment-anchor", type=int, default=50)
    parser.add_argument("--min-score", type=int, default=100)
    parser.add_argument(
        "--min-lift-unmasked",
        type=int,
        default=100,
        help=(
            "lift only global queries with more than this many uppercase "
            "A/C/G/T bases; skipped queries remain globally novel "
            "(default: 100)"
        ),
    )
    parser.add_argument("--min-identity", type=float, default=95.0)
    parser.add_argument(
        "--minimap-query-batch", default="50M",
        help=(
            "Minimap2 -K query bases retained per internal mini-batch "
            "(default: 50M; upstream default is 500M)"
        ),
    )
    parser.add_argument("--blast-word-size", type=int, default=50)
    parser.add_argument("--blast-evalue", default="1e-300")
    parser.add_argument("--blast-max-target-seqs", type=int, default=100)
    parser.add_argument("--skip-blastn", action="store_true")
    parser.add_argument("--resume", action="store_true")
    parser.add_argument(
        "--query-cache-dir",
        help=(
            "reuse an existing candidate_sequences.fa and "
            "candidate_sequences.index.tsv (for example "
            "small_novel_mapping/intermediates/query_sequence_cache)"
        ),
    )
    parser.add_argument(
        "--reuse-alignments-only",
        action="store_true",
        help=(
            "rebuild reports from matching saved query FASTAs and PAF/SAM "
            "files; fail instead of running an aligner"
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
    if args.min_lift_unmasked < 0:
        parser.error("--min-lift-unmasked cannot be negative")
    if args.alignment_anchor < 0:
        parser.error("--alignment-anchor cannot be negative")
    if not 0 < args.min_identity <= 100:
        parser.error("--min-identity must be in (0,100]")
    if re.fullmatch(r"[1-9][0-9]*[kKmMgG]?", args.minimap_query_batch) is None:
        parser.error(
            "--minimap-query-batch must be a positive integer with an "
            "optional K/M/G suffix"
        )
    return args


def run(args: argparse.Namespace) -> None:
    analysis_dir = os.path.abspath(args.analysis_dir)
    local_mapping_dir = os.path.abspath(args.local_mapping_dir)
    main_reference = os.path.abspath(args.main_reference)
    query_paths = os.path.abspath(args.query_paths)
    output_dir = os.path.abspath(args.output_dir)
    reuse_alignments_only = bool(
        getattr(args, "reuse_alignments_only", False)
    )
    Path(output_dir).mkdir(parents=True, exist_ok=True)

    partition_inputs = os.path.join(analysis_dir, "partition_inputs.tsv")
    local_candidates_path = os.path.join(
        analysis_dir, "local_novel_candidates.tsv",
    )
    local_placements_path = os.path.join(
        local_mapping_dir, "all_local_novel_alignments.tsv",
    )
    global_catalog_path = os.path.join(
        analysis_dir, "global_novel.bed",
    )
    if not os.path.isfile(global_catalog_path):
        global_catalog_path = os.path.join(
            analysis_dir, "interval_status.tsv",
        )
    required = (
        partition_inputs,
        os.path.join(analysis_dir, "all_unique_regions.bed"),
        global_catalog_path,
        local_candidates_path,
        local_placements_path,
        main_reference,
        query_paths,
    )
    for path in required:
        if not os.path.isfile(path):
            raise FileNotFoundError(path)

    assembly_sources = read_assembly_sources(query_paths)
    signature_files = [
        *required,
        *(
            path
            for source in sorted(
                assembly_sources.values(),
                key=lambda row: row.haplotype,
            )
            for path in (source.fasta_path, source.fai_path)
        ),
    ]
    common_signature_payload = {
        "files": [_fingerprint(path) for path in signature_files],
        "alignment_anchor": args.alignment_anchor,
        "min_score": args.min_score,
        "min_identity": args.min_identity,
        "blast_word_size": args.blast_word_size,
        "blast_evalue": args.blast_evalue,
        "blast_max_target_seqs": args.blast_max_target_seqs,
        "skip_blastn": args.skip_blastn,
    }
    signature_payload = {
        **common_signature_payload,
        "version": "global-novel-mapping-v5-graphic-annotation",
        "annotation_reference_policy": ANNOTATION_POLICY,
        "min_lift_unmasked": args.min_lift_unmasked,
    }
    signature = hashlib.sha256(
        json.dumps(signature_payload, sort_keys=True).encode()
    ).hexdigest()
    # Extraction is independent of the alignment algorithm.  Keep its cache
    # signature separate, while accepting the immediately preceding mapping
    # signature so a cohort extracted by v2 is not read from 1,000+ assemblies
    # again merely because convergence was enabled.
    cache_signature = hashlib.sha256(json.dumps({
        "version": "global-query-cache-v1",
        "files": common_signature_payload["files"],
        "alignment_anchor": args.alignment_anchor,
    }, sort_keys=True).encode()).hexdigest()
    legacy_cache_signature = hashlib.sha256(json.dumps({
        **common_signature_payload,
        "version": "global-novel-mapping-v2-disk-shards",
    }, sort_keys=True).encode()).hexdigest()
    marker = os.path.join(
        output_dir, "_GLOBAL_NOVEL_MAPPING_SUCCESS.json",
    )
    expected_outputs = (
        os.path.join(output_dir, "all_global_queries.tsv"),
        os.path.join(output_dir, "all_global_reference_placements.tsv"),
        os.path.join(output_dir, "global_novel_alignment_summary.tsv"),
        os.path.join(
            output_dir, "all_global_reference_reporting_alignments.tsv",
        ),
        os.path.join(output_dir, "final_global_novel.tsv"),
    )
    if args.resume and os.path.isfile(marker):
        try:
            with open(marker, "rt") as handle:
                saved = json.load(handle)
            if (
                saved.get("signature") == signature
                and saved.get("placement_resolution")
                == PLACEMENT_RESOLUTION_VERSION
                and saved.get("reporting_resolution")
                == REPORTING_RESOLUTION_VERSION
                and all(os.path.isfile(path) for path in expected_outputs)
            ):
                if saved.get("query_paths") != query_paths:
                    saved["query_paths"] = query_paths
                    temporary = marker + ".tmp"
                    with open(temporary, "wt") as out:
                        json.dump(saved, out, sort_keys=True)
                        out.write("\n")
                    os.replace(temporary, marker)
                LOG.info(
                    "RESUME: global-novel mapping is complete in %s",
                    output_dir,
                )
                return
        except (OSError, ValueError, AttributeError):
            pass

    partitions = read_partition_inputs(partition_inputs)
    local_candidates = read_candidates(local_candidates_path)
    local_placements = read_placements([local_placements_path])
    reference_records = annotation_reference_records(main_reference)
    reference_sources = annotation_reference_sources(assembly_sources, main_reference)
    missing_original_references = {
        record.haplotype for record in reference_records.values()
    } - set(reference_sources)
    if missing_original_references:
        raise ValueError(
            "--query-paths lacks original uncleaned reference "
            "haplotype(s) required by main_chroms.fa: "
            + ",".join(sorted(missing_original_references))
        )
    local_best, _local_masks = choose_local_placements(
        local_candidates, local_placements, reference_records,
    )
    original_global = global_queries_from_bed(analysis_dir, partitions)
    promoted = promoted_local_queries(local_candidates, local_best)
    all_queries = [*original_global, *promoted]
    if reuse_alignments_only:
        # Reporting recovery must replay the exact query population used by
        # the completed alignment.  The reconstructed pool above is the
        # pre-score-filter population, whereas all_global_queries.tsv is the
        # frozen, score-passing population whose sequences were placed in the
        # saved alignment query FASTAs and persistent query cache.
        frozen_queries_path = os.path.join(
            output_dir, "all_global_queries.tsv",
        )
        if not os.path.isfile(frozen_queries_path):
            raise FileNotFoundError(frozen_queries_path)
        all_queries = read_novel_queries(frozen_queries_path)
        original_global = [
            query for query in all_queries
            if query.origin == "global_novel"
        ]
        promoted = [
            query for query in all_queries
            if query.origin != "global_novel"
        ]
        LOG.info(
            "REUSE: frozen score-passing global query selection from %s "
            "(%d queries)",
            frozen_queries_path, len(all_queries),
        )
    missing_sources = {
        query.haplotype
        for query in all_queries
        if query.haplotype not in assembly_sources
    }
    if missing_sources:
        raise ValueError(
            "--query-paths lacks global-query haplotypes: "
            + ",".join(sorted(missing_sources))
        )

    all_candidates = global_mapping_candidates(all_queries)
    cache_dir = os.path.join(
        output_dir, "global_query_extraction.work",
    )
    external_cache = (
        os.path.abspath(args.query_cache_dir)
        if getattr(args, "query_cache_dir", None) else None
    )
    if external_cache is not None:
        cache_fasta = os.path.join(
            external_cache, "candidate_sequences.fa",
        )
        cache_index = os.path.join(
            external_cache, "candidate_sequences.index.tsv",
        )
        if not os.path.isfile(cache_fasta):
            raise FileNotFoundError(cache_fasta)
        if not os.path.isfile(cache_index):
            raise FileNotFoundError(cache_index)
        cached_queries = DiskCandidateSequenceStore(
            cache_fasta, cache_index,
        )
        missing_cache = {
            candidate.candidate_id for candidate in all_candidates
        } - set(cached_queries)
        if missing_cache:
            cached_queries.close()
            raise ValueError(
                f"external query cache lacks {len(missing_cache)} global "
                "queries; first: "
                + ",".join(sorted(missing_cache)[:20])
            )
        LOG.info(
            "REUSE: external disk-backed global query cache %s",
            cache_fasta,
        )
    else:
        if reuse_alignments_only:
            raise ValueError(
                "--reuse-alignments-only requires --query-cache-dir so "
                "report recovery never rereads 1,000+ source assemblies"
            )
        header_records = load_candidate_header_records(
            all_candidates, partitions, args.jobs,
        )
        cached_queries = build_disk_query_cache(
            all_candidates,
            header_records,
            assembly_sources,
            args.alignment_anchor,
            args.jobs,
            cache_dir,
            cache_signature,
            args.resume,
            compatible_signatures=(legacy_cache_signature,),
        )
        del header_records
    succeeded = False
    try:
        if reuse_alignments_only:
            # These were already filtered before all_global_queries.tsv was
            # written.  Re-filtering is unnecessary and could silently alter
            # the completed alignment's query identity if thresholds or
            # sequence sources changed after the original run.
            queries = list(all_queries)
        else:
            queries = filter_queries_by_score(
                all_queries, args.min_score, cached_queries,
            )
        lift_queries = filter_queries_for_global_lift(
            queries, args.min_lift_unmasked, cached_queries,
        )
        retained_ids = {query.query_id for query in lift_queries}
        candidates = [
            candidate
            for candidate in all_candidates
            if candidate.candidate_id in retained_ids
        ]
        LOG.info(
            "Global-reference input: %d original global + %d local residual = "
            "%d score-passing queries; %d have >%d unmasked bases and will "
            "be lifted; %d remain global without lifting",
            len(original_global), len(promoted), len(queries),
            len(lift_queries), args.min_lift_unmasked,
            len(queries) - len(lift_queries),
        )
        workdir = os.path.join(output_dir, "global_alignment_work")
        Path(workdir).mkdir(parents=True, exist_ok=True)
        placements = align_global_queries(
            lift_queries,
            candidates,
            cached_queries,
            main_reference,
            reference_records,
            workdir,
            args,
            signature,
        )
        placements = annotate_global_placements(
            lift_queries,
            placements,
            cached_queries,
            reference_sources,
        )
        reporting_path = os.path.join(
            workdir, "partitions", "global_reference",
            "reporting_placements.tsv",
        )
        reporting_components = read_reporting_components((reporting_path,))
        reporting_components = annotate_global_reporting_components(
            lift_queries,
            reporting_components,
            cached_queries,
            reference_sources,
        )
        write_reporting_components(
            os.path.join(
                output_dir,
                "all_global_reference_reporting_alignments.tsv",
            ),
            reporting_components,
        )
        # Global lifting is annotation-only. Even mapped query bases remain
        # intact and available as complete alternative paths.
        final_global = list(queries)
        if reuse_alignments_only:
            # Preserve the exact completed cleaning/query tables.  Their
            # fingerprints are embedded in downstream small-novel cache
            # signatures; only the new reporting table is regenerated.
            for frozen_name in (
                "all_global_queries.tsv",
                "all_global_reference_placements.tsv",
                "final_global_novel.tsv",
            ):
                frozen_path = os.path.join(output_dir, frozen_name)
                if not os.path.isfile(frozen_path):
                    raise FileNotFoundError(frozen_path)
            write_global_alignment_summary(queries, placements, output_dir)
        else:
            write_global_mapping_outputs(
                queries, final_global, placements, output_dir,
            )
        temporary = marker + ".tmp"
        with open(temporary, "wt") as out:
            json.dump({
                "version": "global-novel-mapping-v5-graphic-annotation",
                "annotation_reference_policy": ANNOTATION_POLICY,
                "placement_resolution": PLACEMENT_RESOLUTION_VERSION,
                "reporting_resolution": REPORTING_RESOLUTION_VERSION,
                "signature": signature,
                "global_queries": len(queries),
                "lift_eligible_queries": len(lift_queries),
                "global_placements": len(placements),
                "global_reporting_components": len(reporting_components),
                "final_global_queries": len(final_global),
                "original_global_queries": len(original_global),
                "promoted_local_queries": len(promoted),
                "min_score": args.min_score,
                "main_reference": main_reference,
                "query_paths": query_paths,
            }, out, sort_keys=True)
            out.write("\n")
        os.replace(temporary, marker)
        succeeded = True
        LOG.info(
            "Annotated %d/%d lift-eligible global queries; retained all %d "
            "global query sequences intact in %s",
            len(placements), len(lift_queries), len(final_global), output_dir,
        )
    finally:
        cached_queries.close()
        if succeeded and external_cache is None:
            shutil.rmtree(cache_dir, ignore_errors=True)


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
