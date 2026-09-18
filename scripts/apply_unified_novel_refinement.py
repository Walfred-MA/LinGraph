#!/usr/bin/env python3
"""Apply unified nearby-template refinement and annotate only its residuals.

No aligner is run here.  The refined template intervals and exclusion masks
come exclusively from ``map_all_novels_to_local_templates.py``.  Any number
of old alignment tables may be supplied for reporting; their paired blocks
are clipped to the retained residual intervals and cannot alter templates or
exclude additional sequence.
"""
from __future__ import annotations

import argparse
import csv
import dataclasses
import hashlib
import heapq
import json
import logging
import os
import shutil
from collections import defaultdict
from contextlib import nullcontext
from pathlib import Path
from typing import DefaultDict, Dict, Iterable, List, Mapping, Optional, Sequence, Tuple

from map_partition_local_novel import read_assembly_sources, read_candidates
from minsetref_core import IndexedFasta
from minsetref_graphic import build_multilocus_graphic_path
from minsetref_light import MASKED_BASE_WEIGHT
from minsetref_multilocus import PlacementEvidence, resolve_query_evidence
from minsetref_segments import novelty_score
from refine_partition_paths import (
    NovelQuery,
    UsedReferenceRow,
    _bounded_process_map,
    build_adjusted_paths,
    expected_partition_sample_fasta,
    parse_main_reference,
    partition_sample_fasta,
    read_novel_queries,
    read_partition_inputs,
    read_unique_rows,
    read_used_reference_rows,
    write_adjusted_path_outputs,
)
from summarize_partition_novelty import write_tsv
from fixed_alternatives import ANNOTATION_HAPLOTYPE, ANNOTATION_POLICY


LOG = logging.getLogger("apply_unified_novel_refinement")
Interval = Tuple[int, int]
PairedBlock = Tuple[int, int, int, int]
VERSION = "apply-unified-novel-refinement-v8-verified-outputs"


def select_used_templates(original_used, refined_main, fixed, main_haplotype):
    """Route imported originals per partition and retain refinable novels."""
    fixed_only = {name for name, records in fixed.items()
                  if all(row.encoded_name in records for row in original_used
                         if row.partition == name)}
    refined_main = [row for row in refined_main if row.partition not in fixed_only]
    refined_partitions = {row.partition for row in refined_main}
    used_rows = [row for row in original_used
                 if row.encoded_name in fixed.get(row.partition, {})
                 or row.haplotype != main_haplotype
                 or row.partition not in refined_partitions] + refined_main
    return used_rows, refined_partitions


@dataclasses.dataclass(frozen=True)
class OldQuery:
    query_id: str
    partition: str
    record_id: str
    local_start: int
    local_end: int


@dataclasses.dataclass(frozen=True)
class RawEvidence:
    owner: OldQuery
    mapped_partition: str
    ref_haplotype: str
    ref_contig: str
    strand: str
    target_record: str
    query_start: int
    query_end: int
    paired_blocks: Tuple[PairedBlock, ...]
    alignment_score: float
    aligned_bases: int
    identity: float
    source: str
    ordinal: int


def read_rows(path: str):
    with open(path, "rt", newline="") as handle:
        first = handle.readline()
        if not first:
            return
        if first.startswith("#"):
            first = first[1:]
        fields = next(csv.reader([first.rstrip("\n")], delimiter="\t"))
        reader = csv.DictReader(handle, fieldnames=fields, delimiter="\t")
        for row in reader:
            yield row


def json_intervals(text: str, context: str, width: int) -> Tuple[Tuple[int, ...], ...]:
    if not text or text == ".":
        return ()
    try:
        value = json.loads(text)
    except json.JSONDecodeError as error:
        raise ValueError(f"{context}: invalid JSON intervals") from error
    output = []
    for item in value:
        if not isinstance(item, list) or len(item) != width:
            raise ValueError(f"{context}: expected {width}-column intervals")
        numbers = tuple(int(number) for number in item)
        if numbers[1] <= numbers[0]:
            raise ValueError(f"{context}: invalid interval {numbers}")
        if width == 4 and (
            numbers[3] <= numbers[2]
            or numbers[1] - numbers[0] != numbers[3] - numbers[2]
        ):
            raise ValueError(f"{context}: invalid paired block {numbers}")
        output.append(numbers)
    return tuple(output)


def load_old_query_catalog(
    local_candidates_path: Optional[str],
    global_queries_path: Optional[str],
) -> Dict[str, OldQuery]:
    output: Dict[str, OldQuery] = {}
    if local_candidates_path and os.path.isfile(local_candidates_path):
        for row in read_candidates(local_candidates_path):
            output[row.candidate_id] = OldQuery(
                row.candidate_id, row.novel_partition, row.record_id,
                row.local_start, row.local_end,
            )
    if global_queries_path and os.path.isfile(global_queries_path):
        for row in read_novel_queries(global_queries_path):
            output[row.query_id] = OldQuery(
                row.query_id, row.partition, row.record_id,
                row.local_start, row.local_end,
            )
    return output


def alignment_owner_id(row: Mapping[str, str], catalog: Mapping[str, OldQuery]) -> Optional[str]:
    # Small-novel rows use a residual query_id but preserve the original
    # global ID. Prefer the latter whenever it is known.
    choices = (
        row.get("original_query_id"), row.get("candidate_id"),
        row.get("query_id"),
    )
    return next((value for value in choices if value in catalog), None)


def read_old_evidence(
    paths: Sequence[str], catalog: Mapping[str, OldQuery],
) -> List[RawEvidence]:
    output: List[RawEvidence] = []
    seen = set()
    ordinal = 0
    for path in paths:
        before = len(output)
        LOG.info("Reading old annotation evidence: %s", path)
        for line_number, row in enumerate(read_rows(path), 2):
            if line_number % 250000 == 0:
                LOG.info(
                    "Old evidence scan: %s %d rows read; %d usable "
                    "rows retained from this file",
                    os.path.basename(path), line_number - 1,
                    len(output) - before,
                )
            owner_id = alignment_owner_id(row, catalog)
            if owner_id is None:
                continue
            paired = json_intervals(
                row.get("paired_blocks", "[]"),
                f"{path}:{line_number}:paired_blocks", 4,
            )
            if not paired:
                continue
            ref_haplotype = row.get("ref_haplotype") or row.get("target_haplotype")
            ref_contig = row.get("ref_contig") or row.get("target_contig")
            strand = row.get("strand")
            if not ref_haplotype or not ref_contig or strand not in {"+", "-"}:
                continue
            query_blocks = json_intervals(
                row.get("query_blocks", "[]"),
                f"{path}:{line_number}:query_blocks", 2,
            )
            query_start = int(row["query_start"]) if row.get("query_start") not in {None, "", "."} else (
                min(block[0] for block in query_blocks) if query_blocks
                else min(block[0] for block in paired)
            )
            query_end = int(row["query_end"]) if row.get("query_end") not in {None, "", "."} else (
                max(block[1] for block in query_blocks) if query_blocks
                else max(block[1] for block in paired)
            )
            owner = catalog[owner_id]
            query_start = max(query_start, owner.local_start)
            query_end = min(query_end, owner.local_end)
            paired = tuple(
                block for block in paired
                if block[0] < query_end and query_start < block[1]
            )
            if query_end <= query_start or not paired:
                continue
            identity_key = (
                owner_id, ref_haplotype, ref_contig, strand, paired,
            )
            if identity_key in seen:
                continue
            seen.add(identity_key)
            ordinal += 1
            output.append(RawEvidence(
                owner,
                row.get("mapped_partition") or row.get("partition") or
                row.get("target_kind") or ".",
                ref_haplotype, ref_contig, strand,
                row.get("target_record") or ref_contig,
                query_start, query_end, paired,
                float(row.get("alignment_score") or row.get("aligned_score") or 0),
                int(float(row.get("aligned_bases") or sum(
                    block[1] - block[0] for block in paired
                ))),
                float(row.get("identity") or 0),
                row.get("source") or os.path.basename(path), ordinal,
            ))
        LOG.info(
            "Old evidence loaded: %s contributed %d usable rows "
            "(%d cumulative)",
            os.path.basename(path), len(output) - before, len(output),
        )
    return output


def clip_pair(block: PairedBlock, strand: str, start: int, end: int) -> Optional[PairedBlock]:
    q0, q1, r0, _r1 = block
    left, right = max(q0, start), min(q1, end)
    if right <= left:
        return None
    if strand == "+":
        return left, right, r0 + left - q0, r0 + right - q0
    return left, right, r0 + q1 - right, r0 + q1 - left


def read_target_sources(
    query_paths: str,
    target_fasta_args: Sequence[str],
) -> Tuple[Dict[str, Tuple[str, Optional[str]]], Dict[str, object]]:
    assembly_sources = read_assembly_sources(query_paths)
    sources: Dict[str, Tuple[str, str]] = {
        haplotype: (row.fasta_path, row.fai_path)
        for haplotype, row in assembly_sources.items()
    }
    for value in target_fasta_args:
        if "=" not in value:
            raise ValueError("--target-fasta must be HAPLOTYPE=FASTA")
        haplotype, fasta = value.split("=", 1)
        fasta = os.path.abspath(fasta)
        fai = fasta + ".fai"
        if not os.path.isfile(fasta):
            raise FileNotFoundError(fasta)
        sources[haplotype] = (
            fasta, fai if os.path.isfile(fai) else None,
        )
    return sources, assembly_sources


def annotation_labels(
    queries: Sequence[NovelQuery],
    evidence_rows: Sequence[RawEvidence],
    target_sources: Mapping[str, Tuple[str, Optional[str]]],
    jobs: int = 1,
) -> Tuple[Dict[str, str], List[Tuple[object, ...]]]:
    by_record: DefaultDict[Tuple[str, str], List[RawEvidence]] = defaultdict(list)
    for row in evidence_rows:
        by_record[(row.owner.partition, row.owner.record_id)].append(row)
    readers: Dict[str, IndexedFasta] = {}
    names: Dict[str, set] = {}

    def target_reader(haplotype: str) -> Optional[IndexedFasta]:
        source = target_sources.get(haplotype)
        if source is None:
            return None
        if haplotype not in readers:
            readers[haplotype] = IndexedFasta(*source)
            names[haplotype] = set(readers[haplotype].names())
        return readers[haplotype]

    labels: Dict[str, str] = {}
    audit: List[Tuple[object, ...]] = []
    by_source: DefaultDict[str, List[NovelQuery]] = defaultdict(list)
    for query in queries:
        by_source[query.source_fasta].append(query)
    if jobs > 1 and len(by_source) > 1:
        source_groups = []
        for source_fasta, members in by_source.items():
            keys = {(row.partition, row.record_id) for row in members}
            group_evidence = [
                row for key in keys for row in by_record.get(key, ())
            ]
            cost = sum(row.local_end - row.local_start for row in members)
            cost += sum(len(row.paired_blocks) * 100 for row in group_evidence)
            source_groups.append((
                cost, source_fasta, tuple(members), tuple(group_evidence),
            ))
        bucket_count = min(len(source_groups), max(1, jobs * 4))
        buckets = [(0, index, []) for index in range(bucket_count)]
        heapq.heapify(buckets)
        for cost, source_fasta, members, group_evidence in sorted(
            source_groups, key=lambda row: (-row[0], row[1]),
        ):
            load, bucket_index, groups = heapq.heappop(buckets)
            groups.append((source_fasta, members, group_evidence))
            heapq.heappush(
                buckets, (load + cost, bucket_index, groups),
            )
        tasks = []
        for _load, _bucket_index, groups in sorted(
            buckets, key=lambda row: row[1],
        ):
            bucket_queries = tuple(
                query for _source, members, _evidence in groups
                for query in members
            )
            bucket_evidence = tuple(
                row for _source, _members, evidence in groups
                for row in evidence
            )
            required_haplotypes = {
                row.ref_haplotype for row in bucket_evidence
            }
            tasks.append((
                bucket_queries,
                bucket_evidence,
                {
                    haplotype: target_sources[haplotype]
                    for haplotype in required_haplotypes
                    if haplotype in target_sources
                },
            ))
        LOG.info(
            "Annotating %d residual intervals from %d local graphs using "
            "%d process workers across %d balanced tasks",
            len(queries), len(by_source), min(jobs, len(tasks)), len(tasks),
        )
        combined_labels: Dict[str, str] = {}
        combined_audit: List[Tuple[object, ...]] = []
        progress_step = max(1, len(tasks) // 20)
        for completed, (task_labels, task_audit) in enumerate(
            _bounded_process_map(
                _annotation_group_worker, tasks, min(jobs, len(tasks)),
            ),
            1,
        ):
            combined_labels.update(task_labels)
            combined_audit.extend(task_audit)
            if completed == len(tasks) or completed % progress_step == 0:
                LOG.info(
                    "Residual annotation progress: %d/%d tasks; %d/%d "
                    "intervals complete",
                    completed, len(tasks), len(combined_labels), len(queries),
                )
        combined_audit.sort(key=lambda row: (
            str(row[1]), str(row[2]), int(row[3]), int(row[4]), str(row[0]),
        ))
        return combined_labels, combined_audit
    try:
        for source_fasta, members in sorted(by_source.items()):
            # Interval-only graph CIGAR does not need query bases. Keep the
            # source grouping for balanced per-local-graph scheduling without
            # opening tens of thousands of query FASTAs here.
            with nullcontext() as source:
                for query in members:
                    usable: List[PlacementEvidence] = []
                    for raw in by_record.get((query.partition, query.record_id), ()):
                        if raw.owner.local_end <= query.local_start or raw.owner.local_start >= query.local_end:
                            continue
                        reader = target_reader(raw.ref_haplotype)
                        if reader is None or raw.ref_contig not in names[raw.ref_haplotype]:
                            continue
                        paired = tuple(
                            clipped for block in raw.paired_blocks
                            for clipped in (clip_pair(
                                block, raw.strand,
                                query.local_start, query.local_end,
                            ),)
                            if clipped is not None
                        )
                        if not paired:
                            continue
                        qstart = max(
                            query.local_start,
                            min(raw.query_start, min(
                                q0 for q0, _q1, _r0, _r1 in paired
                            )),
                        )
                        qend = min(
                            query.local_end,
                            max(raw.query_end, max(
                                q1 for _q0, q1, _r0, _r1 in paired
                            )),
                        )
                        usable.append(PlacementEvidence(
                            query.query_id, raw.mapped_partition,
                            raw.ref_haplotype, raw.ref_contig, raw.strand,
                            raw.target_record, qstart, qend, paired,
                            raw.alignment_score,
                            sum(q1 - q0 for q0, q1, _r0, _r1 in paired),
                            raw.identity, raw.source, raw.ordinal,
                        ))
                    components = resolve_query_evidence(usable)
                    sequence = ""
                    if components:
                        references = {
                            (component.ref_haplotype, component.ref_contig):
                            readers[component.ref_haplotype]
                            for component in components
                        }
                        label = build_multilocus_graphic_path(
                            query_id=query.query_id,
                            query_start=query.local_start,
                            query_end=query.local_end,
                            query_sequence=sequence,
                            components=components,
                            references=references,
                            exact_bases=False,
                        )
                    else:
                        label = f"{query.local_end - query.local_start}I"
                    labels[query.query_id] = label
                    audit.append((
                        query.query_id, query.partition, query.record_id,
                        query.local_start, query.local_end, len(usable),
                        len(components), label,
                    ))
    finally:
        for reader in readers.values():
            reader.close()
    return labels, audit


def _annotation_group_worker(task):
    queries, evidence_rows, target_sources = task
    return annotation_labels(
        queries, evidence_rows, target_sources, jobs=1,
    )


def _score_query_bucket(task):
    groups, min_score = task
    retained = []
    rejected = []
    for source_fasta, queries in groups:
        with IndexedFasta(source_fasta) as source:
            for query in queries:
                sequence = source.fetch(
                    query.record_id, query.local_start, query.local_end,
                )
                if novelty_score(sequence, MASKED_BASE_WEIGHT) >= min_score:
                    retained.append(query)
                else:
                    rejected.append(query)
    return retained, rejected


def filter_residual_queries_parallel(
    queries: Sequence[NovelQuery], min_score: int, jobs: int,
) -> Tuple[List[NovelQuery], List[NovelQuery]]:
    if not queries:
        return [], []
    by_source: DefaultDict[str, List[NovelQuery]] = defaultdict(list)
    for query in queries:
        by_source[query.source_fasta].append(query)
    groups = [
        (
            sum(row.local_end - row.local_start for row in members),
            source_fasta,
            tuple(members),
        )
        for source_fasta, members in by_source.items()
    ]
    bucket_count = min(len(groups), max(1, jobs * 4))
    buckets = [(0, index, []) for index in range(bucket_count)]
    heapq.heapify(buckets)
    for cost, source_fasta, members in sorted(
        groups, key=lambda row: (-row[0], row[1]),
    ):
        load, bucket_index, bucket = heapq.heappop(buckets)
        bucket.append((source_fasta, members))
        heapq.heappush(buckets, (load + cost, bucket_index, bucket))
    tasks = [
        (tuple(bucket), min_score)
        for _load, _index, bucket in sorted(buckets, key=lambda row: row[1])
        if bucket
    ]
    LOG.info(
        "Checking the size/novelty cutoff for %d residual intervals using "
        "%d process workers across %d tasks",
        len(queries), min(jobs, len(tasks)), len(tasks),
    )
    retained: List[NovelQuery] = []
    rejected: List[NovelQuery] = []
    progress_step = max(1, len(tasks) // 20)
    for completed, (kept, dropped) in enumerate(
        _bounded_process_map(
            _score_query_bucket, tasks, min(jobs, len(tasks)),
        ),
        1,
    ):
        retained.extend(kept)
        rejected.extend(dropped)
        if completed == len(tasks) or completed % progress_step == 0:
            LOG.info(
                "Residual cutoff progress: %d/%d tasks; %d intervals "
                "retained and %d rejected",
                completed, len(tasks), len(retained), len(rejected),
            )
    retained.sort(key=lambda row: row.query_id)
    rejected.sort(key=lambda row: row.query_id)
    return retained, rejected


def filter_residual_queries_from_assemblies(
    queries: Sequence[NovelQuery],
    query_paths: str,
    min_score: int,
) -> Tuple[List[NovelQuery], List[NovelQuery]]:
    """Score residuals from original FAI sources, never `_samples.fasta`."""
    sources = read_assembly_sources(query_paths)
    by_haplotype: DefaultDict[str, List[NovelQuery]] = defaultdict(list)
    for query in queries:
        if query.haplotype not in sources:
            raise ValueError(
                f"{query.query_id}: haplotype {query.haplotype!r} is absent "
                f"from {query_paths}"
            )
        by_haplotype[query.haplotype].append(query)
    retained: List[NovelQuery] = []
    rejected: List[NovelQuery] = []
    for haplotype in sorted(by_haplotype):
        source = sources[haplotype]
        with IndexedFasta(source.fasta_path, source.fai_path) as reader:
            for query in by_haplotype[haplotype]:
                if query.source_contig not in reader.index:
                    raise KeyError(
                        f"{query.query_id}: {query.source_contig!r} is absent "
                        f"from {source.fasta_path}"
                    )
                sequence = reader.fetch(
                    query.source_contig, query.source_start, query.source_end,
                )
                destination = (
                    retained
                    if novelty_score(sequence, MASKED_BASE_WEIGHT) >= min_score
                    else rejected
                )
                destination.append(query)
    retained.sort(key=lambda row: row.query_id)
    rejected.sort(key=lambda row: row.query_id)
    return retained, rejected


def read_unified_inputs(
    unified_dir: str,
    partitions,
) -> Tuple[List[NovelQuery], Dict[Tuple[str, str], List[Interval]], List[UsedReferenceRow], List[Tuple[str, ...]]]:
    candidates = {
        row["candidate_id"]: row
        for row in read_rows(os.path.join(unified_dir, "all_novel_candidates.tsv"))
    }
    residual_queries: List[NovelQuery] = []
    for index, row in enumerate(read_rows(
        os.path.join(unified_dir, "residual_novel_intervals.tsv")
    ), 1):
        candidate = candidates[row["candidate_id"]]
        partition = row["partition"]
        residual_queries.append(NovelQuery(
            f"unified_residual_{index:012d}", "unified_residual",
            partition, row["candidate_id"],
            partition_sample_fasta(partitions[partition]),
            row["record_id"], row["haplotype"], row["source_contig"],
            int(row["source_start"]), int(row["source_end"]),
            int(row["local_start"]), int(row["local_end"]),
        ))
    masks: DefaultDict[Tuple[str, str], List[Interval]] = defaultdict(list)
    covered_audit = []
    for row in read_rows(os.path.join(unified_dir, "reference_covered_novel.tsv")):
        candidate = candidates[row["candidate_id"]]
        masks[(row["partition"], row["record_id"])].append((
            int(row["local_start"]), int(row["local_end"]),
        ))
        covered_audit.append(tuple(row.values()))
    refined_rows = [
        UsedReferenceRow(
            row["partition"], row["haplotype"], row["contig"],
            int(row["start"]), int(row["end"]),
            (
                f"unified_refined_{row['partition']}_{row['contig']}_"
                f"{row['start']}_{row['end']}"
            ),
        )
        for row in read_rows(os.path.join(
            unified_dir, "refined_template_intervals.tsv"
        ))
    ]
    return residual_queries, dict(masks), refined_rows, covered_audit


def fingerprint(path: str):
    stat = os.stat(path)
    return os.path.abspath(path), stat.st_size, stat.st_mtime_ns


def output_fingerprints(output_dir, partitions, annotations_only):
    """Validate every promised table, including the per-partition manifests."""
    names = ["residual_novel_annotations.tsv"]
    suffix = "annotations" if annotations_only else "adjusted"
    names.extend(f"{name}/{name}.{suffix}.tsv" for name in sorted(partitions))
    if not annotations_only:
        names.extend((
            "all_adjusted_paths.tsv", "partition_sources.tsv",
            "reference_covered_mapped_novel.tsv",
        ))
    result = []
    for name in names:
        path = os.path.join(output_dir, name)
        if not os.path.isfile(path) or os.path.getsize(path) == 0:
            raise ValueError(f"missing or empty refinement output: {path}")
        result.append(list(fingerprint(path)))
    return result


def parse_args(argv: Optional[Sequence[str]] = None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("-a", "--analysis-dir", required=True)
    parser.add_argument("-u", "--unified-refinement-dir", required=True)
    parser.add_argument("-q", "--query-paths", required=True)
    parser.add_argument("-r", "--main-reference", required=True)
    parser.add_argument("-o", "--output-dir", required=True)
    parser.add_argument("--local-candidates")
    parser.add_argument("--global-queries")
    parser.add_argument(
        "--alignment-file", action="append", default=[], metavar="TSV",
        help=(
            "old alignment evidence table; repeat for any number of local, "
            "global, large-novel, or small-novel tables"
        ),
    )
    parser.add_argument(
        "--target-fasta", action="append", default=[],
        metavar="HAPLOTYPE=FASTA",
        help=(
            "target sequence source for non-assembly evidence, e.g. "
            "large_novel=novel_loci.fa or small_novel=smallnovel.fa"
        ),
    )
    parser.add_argument("--min-score", type=int, default=100)
    parser.add_argument(
        "--annotations-only", action="store_true",
        help=(
            "write immutable graph sidecar annotations only; do not create "
            "or rewrite adjusted graph-coordinate manifests"
        ),
    )
    parser.add_argument("-j", "--jobs", type=int, default=16)
    parser.add_argument("--resume", action="store_true")
    parser.add_argument(
        "--log-level", choices=("DEBUG", "INFO", "WARNING", "ERROR"),
        default="INFO",
    )
    args = parser.parse_args(argv)
    if args.jobs < 1 or args.min_score < 1:
        parser.error("--jobs and --min-score must be positive")
    return args


def run(args: argparse.Namespace) -> None:
    analysis_dir = os.path.abspath(args.analysis_dir)
    unified_dir = os.path.abspath(args.unified_refinement_dir)
    query_paths = os.path.abspath(args.query_paths)
    main_reference = os.path.abspath(args.main_reference)
    output_dir = os.path.abspath(args.output_dir)
    local_candidates_path = os.path.abspath(
        args.local_candidates or os.path.join(
            analysis_dir, "local_novel_candidates.tsv",
        )
    )
    global_queries_path = (
        os.path.abspath(args.global_queries) if args.global_queries else None
    )
    alignment_files = [os.path.abspath(path) for path in args.alignment_file]
    target_files = []
    for value in args.target_fasta:
        if "=" not in value:
            raise ValueError("--target-fasta must be HAPLOTYPE=FASTA")
        target_files.append(os.path.abspath(value.split("=", 1)[1]))
    required = [
        os.path.join(analysis_dir, "partition_inputs.tsv"),
        os.path.join(analysis_dir, "_COORDINATE_ANALYSIS_SUCCESS.json"),
        os.path.join(analysis_dir, "all_unique_regions.bed"),
        os.path.join(analysis_dir, "all_used_novel_loci.bed"),
        os.path.join(unified_dir, "_UNIFIED_TEMPLATE_REFINEMENT_SUCCESS.json"),
        os.path.join(unified_dir, "all_novel_candidates.tsv"),
        os.path.join(unified_dir, "residual_novel_intervals.tsv"),
        os.path.join(unified_dir, "reference_covered_novel.tsv"),
        os.path.join(unified_dir, "refined_template_intervals.tsv"),
        query_paths, main_reference, *alignment_files, *target_files,
    ]
    if os.path.isfile(local_candidates_path):
        required.append(local_candidates_path)
    if global_queries_path:
        required.append(global_queries_path)
    for path in required:
        if not os.path.isfile(path):
            raise FileNotFoundError(path)
    partitions = read_partition_inputs(os.path.join(
        analysis_dir, "partition_inputs.tsv",
    ))
    from fixed_alternatives import fixed_partitions
    fixed = fixed_partitions(partitions)
    signature = hashlib.sha256(json.dumps({
        "version": VERSION,
        "annotation_reference_policy": ANNOTATION_POLICY,
        "files": [fingerprint(path) for path in required],
        "fixed_templates": fixed,
        "target_fasta": sorted(args.target_fasta),
        "min_score": args.min_score,
    }, sort_keys=True).encode()).hexdigest()
    marker = os.path.join(
        output_dir,
        "_ANNOTATIONS_SUCCESS.json"
        if args.annotations_only else "_ADJUSTED_PATHS_SUCCESS.json",
    )
    if args.resume and os.path.isfile(marker):
        try:
            with open(marker) as handle:
                saved = json.load(handle)
            if (saved.get("signature") == signature
                    and saved.get("outputs") == output_fingerprints(
                        output_dir, partitions, args.annotations_only)):
                LOG.info("RESUME: verified unified outputs in %s", output_dir)
                return
        except (OSError, ValueError, AttributeError):
            LOG.info("Incomplete refinement outputs or marker; rebuilding %s", output_dir)
    # A failed rebuild must not leave an earlier success marker behind.
    Path(marker).unlink(missing_ok=True)
    Path(output_dir).mkdir(parents=True, exist_ok=True)

    unique_rows = read_unique_rows(analysis_dir, partitions)
    original_used = read_used_reference_rows(analysis_dir)
    reference_records = parse_main_reference(main_reference)
    main_haplotypes = {row.haplotype for row in reference_records.values()}
    with open(os.path.join(
        analysis_dir, "_COORDINATE_ANALYSIS_SUCCESS.json",
    )) as handle:
        analysis_metadata = json.load(handle)
    reference_haplotypes = set(
        analysis_metadata.get("reference_haplotypes", ())
    )
    primary_haplotypes = main_haplotypes & reference_haplotypes
    if len(primary_haplotypes) != 1:
        raise ValueError(
            "main reference catalog must contain exactly one primary source "
            "haplotype named by coordinate-analysis metadata; catalog="
            f"{sorted(main_haplotypes)}, primary={sorted(reference_haplotypes)}"
        )
    main_haplotype = next(iter(primary_haplotypes))
    residuals, masks, refined_main, _covered_audit = read_unified_inputs(
        unified_dir, partitions,
    )
    # A source sample label never authorizes replacing an imported record.
    used_rows, refined_partitions = select_used_templates(
        original_used, refined_main, fixed, main_haplotype,
    )
    used_partitions = {row.partition for row in used_rows}
    missing_used = set(partitions) - used_partitions
    if missing_used:
        raise ValueError(
            "partitions lack both an original and refined block-reference "
            "interval: " + ",".join(sorted(missing_used)[:20])
        )
    LOG.info(
        "Template routing: %d partitions use unified refinements; %d retain "
        "their unchanged original templates",
        len(refined_partitions), len(partitions) - len(refined_partitions),
    )
    residuals, rejected_residuals = filter_residual_queries_from_assemblies(
        residuals, query_paths, args.min_score,
    )
    for query in rejected_residuals:
        masks.setdefault((query.partition, query.record_id), []).append((
            query.local_start, query.local_end,
        ))
    LOG.info(
        "Residual cutoff complete: %d retained for annotation/path output; "
        "%d below-cutoff intervals suppressed",
        len(residuals), len(rejected_residuals),
    )

    catalog = load_old_query_catalog(
        local_candidates_path, global_queries_path,
    )
    LOG.info(
        "Loaded %d old local/global query identifiers for annotation",
        len(catalog),
    )
    evidence = read_old_evidence(alignment_files, catalog)
    target_sources, assembly_sources = read_target_sources(
        query_paths, args.target_fasta,
    )
    target_sources[ANNOTATION_HAPLOTYPE] = (main_reference, main_reference + ".fai")
    labels, annotation_audit = annotation_labels(
        residuals, evidence, target_sources, args.jobs,
    )
    annotation_fields = (
        "query_id", "partition", "record_id", "local_start", "local_end",
        "input_evidence", "resolved_components", "graphic_cigar",
    )
    write_tsv(
        os.path.join(output_dir, "residual_novel_annotations.tsv"),
        annotation_fields,
        annotation_audit,
    )
    if args.annotations_only:
        by_partition: DefaultDict[str, List[Tuple[object, ...]]] = defaultdict(list)
        for row in annotation_audit:
            by_partition[str(row[1])].append(row)
        for partition in sorted(partitions):
            rows_for_partition = by_partition.get(partition, [])
            directory = os.path.join(output_dir, partition)
            Path(directory).mkdir(parents=True, exist_ok=True)
            write_tsv(
                os.path.join(directory, f"{partition}.annotations.tsv"),
                annotation_fields,
                rows_for_partition,
            )
        temporary = marker + ".tmp"
        with open(temporary, "wt") as out:
            json.dump({
                "version": VERSION + "-annotations-only",
                "signature": signature,
                "outputs": output_fingerprints(output_dir, partitions, True),
                "partitions": len(partitions),
                "annotations": len(annotation_audit),
                "below_cutoff_residuals": len(rejected_residuals),
                "immutable_graph_coordinates": True,
                "graph_files_modified": False,
            }, out, sort_keys=True)
            out.write("\n")
        os.replace(temporary, marker)
        LOG.info(
            "Wrote %d immutable sidecar annotations for %d partitions; "
            "no graph/alignment/linear files were modified",
            len(annotation_audit), len(partitions),
        )
        return
    rows = build_adjusted_paths(
        partitions=partitions,
        unique_rows=unique_rows,
        used_rows=used_rows,
        local_candidates=(),
        local_best={},
        local_masks={},
        global_queries=residuals,
        global_best={},
        reference_records=reference_records,
        main_reference=main_reference,
        merge_gap=0,
        min_score=args.min_score,
        assembly_sources=assembly_sources,
        jobs=args.jobs,
        small_novel_annotations=labels,
        precovered_masks=masks,
        force_reference_rows=True,
        include_reference_unique_spans=True,
        suppress_unique_haplotypes=reference_haplotypes,
        prevalidated_query_ids={row.query_id for row in residuals},
    )
    write_adjusted_path_outputs(
        rows, output_dir, args.jobs,
        {
            partition: expected_partition_sample_fasta(files)
            for partition, files in partitions.items()
        },
    )
    # This exact exclusion audit is what build_local_graphs.py records beside
    # the generated graph FASTAs.
    shutil.copyfile(
        os.path.join(unified_dir, "reference_covered_novel.tsv"),
        os.path.join(output_dir, "reference_covered_mapped_novel.tsv"),
    )
    temporary = marker + ".tmp"
    with open(temporary, "wt") as out:
        json.dump({
            "version": VERSION, "signature": signature,
            "outputs": output_fingerprints(output_dir, partitions, False),
            "partitions": len(partitions), "paths": len(rows),
            "residual_queries": len(residuals),
            "below_cutoff_residuals": len(rejected_residuals),
            "old_alignment_files": alignment_files,
            "old_evidence_rows": len(evidence),
            "policy": "old_alignments_annotation_only",
            "template_merge_gap": 0,
            "reference_unique_expansion": True,
        }, out, sort_keys=True)
        out.write("\n")
    os.replace(temporary, marker)
    LOG.info(
        "Wrote %d adjusted paths; annotated %d residual novel intervals",
        len(rows), len(residuals),
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
    raise SystemExit(main())
