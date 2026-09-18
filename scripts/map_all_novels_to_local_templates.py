#!/usr/bin/env python3
"""Refine local-graph templates with unsplit ``.unique.bed`` novels.

This stage deliberately runs before local/global novelty classification.  One
candidate is created for every original non-reference row in every partition's
``.unique.bed`` catalog.  Candidates are never split according to
``interval_status.tsv``. Within each partition they are aligned to the main-
reference template windows. Each local graph gets one query FASTA, one target
FASTA, and one Minimap2 invocation.  The shared target uses the largest
``max(5000, 2 * novel_length)`` flank required by any query in that graph, so
every query receives at least its requested search region. Every qualifying placement
contributes: its aligned query blocks are excluded from the novel and its
target blocks expand the template.  Downstream alignment results are not read
and therefore cannot influence refinement.
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
from pathlib import Path
from typing import DefaultDict, Dict, Iterable, Iterator, List, Mapping, Optional, Sequence, Tuple

from map_partition_local_novel import (
    AssemblyFastaSource,
    MINIMAP_PARAMS,
    read_assembly_sources,
)
from minsetref_align import AlignmentHit, ensure_executable, run_minimap2
from minsetref_core import IndexedFasta, count_unmasked, mp_context, wrap_fasta
from refine_partition_paths import (
    UniqueRow,
    UsedReferenceRow,
    merge_intervals,
    parse_main_reference,
    read_partition_inputs,
    read_unique_rows,
    read_used_reference_rows,
    subtract_intervals,
)
from summarize_partition_novelty import (
    HeaderRecord,
    PartitionFiles,
    parse_header_records,
    write_tsv,
)


LOG = logging.getLogger("map_all_novels_to_local_templates")
Interval = Tuple[int, int]
VERSION = "unified-template-refinement-v5-fixed-imported-alternatives"


@dataclasses.dataclass(frozen=True)
class UnifiedNovel:
    candidate_id: str
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

    @property
    def length(self) -> int:
        return self.local_end - self.local_start


@dataclasses.dataclass(frozen=True)
class TemplateWindow:
    template_id: str
    partition: str
    haplotype: str
    contig: str
    template_start: int
    template_end: int
    window_start: int
    window_end: int


@dataclasses.dataclass(frozen=True)
class AcceptedPlacement:
    placement_id: str
    candidate_id: str
    partition: str
    template_id: str
    ref_haplotype: str
    ref_contig: str
    ref_start: int
    ref_end: int
    strand: str
    aligner: str
    alignment_score: float
    identity: float
    unmasked_aligned_bases: int
    query_blocks: Tuple[Interval, ...]
    target_blocks: Tuple[Interval, ...]
    paired_blocks: Tuple[Tuple[int, int, int, int], ...]


@dataclasses.dataclass(frozen=True)
class PartitionResult:
    partition: str
    candidates: Tuple[UnifiedNovel, ...]
    placements: Tuple[AcceptedPlacement, ...]
    templates: Tuple[TemplateWindow, ...]
    residuals: Tuple[Tuple[str, int, int], ...]


def interval_distance(left: Interval, right: Interval) -> int:
    if left[1] < right[0]:
        return right[0] - left[1]
    if right[1] < left[0]:
        return left[0] - right[1]
    return 0


def unified_candidates(
    rows: Sequence[UniqueRow],
    reference_haplotypes: Iterable[str],
) -> List[UnifiedNovel]:
    references = set(reference_haplotypes)
    output: List[UnifiedNovel] = []
    for row in rows:
        if row.region_class == "reference" or row.haplotype in references:
            continue
        output.append(UnifiedNovel(
            f"unified_novel_{len(output) + 1:012d}",
            row.partition, row.record_id, row.haplotype,
            row.source_contig, row.source_start, row.source_end,
            row.local_start, row.local_end, row.region_name,
            row.region_class,
        ))
    return output


def merge_template_rows(
    rows: Sequence[UsedReferenceRow],
    main_reference: str,
) -> Dict[str, List[Tuple[str, int, int]]]:
    grouped: DefaultDict[Tuple[str, str], List[Interval]] = defaultdict(list)
    for row in rows:
        if row.haplotype == main_reference:
            grouped[(row.partition, row.source_contig)].append((
                row.source_start, row.source_end,
            ))
    output: DefaultDict[str, List[Tuple[str, int, int]]] = defaultdict(list)
    for (partition, contig), intervals in sorted(grouped.items()):
        for start, end in merge_intervals(intervals):
            output[partition].append((contig, start, end))
    return dict(output)


def build_template_windows(
    partition: str,
    templates: Sequence[Tuple[str, int, int]],
    maximum_novel_length: int,
    contig_lengths: Mapping[str, int],
    minimum_flank: int = 5000,
) -> List[TemplateWindow]:
    flank = max(minimum_flank, 2 * maximum_novel_length)
    output = []
    for index, (contig, start, end) in enumerate(templates, 1):
        if contig not in contig_lengths:
            raise ValueError(
                f"{partition}: reference lacks template contig {contig!r}"
            )
        contig_length = contig_lengths[contig]
        if start < 0 or end <= start or end > contig_length:
            raise ValueError(
                f"{partition}: invalid template interval "
                f"{contig}:{start}-{end} for contig length "
                f"{contig_length}"
            )
        output.append(TemplateWindow(
            f"template_{index:06d}", partition, "", contig, start, end,
            max(0, start - flank), min(contig_length, end + flank),
        ))
    return output


def write_partition_fastas(
    candidates: Sequence[UnifiedNovel],
    templates: Sequence[TemplateWindow],
    header_records: Mapping[str, HeaderRecord],
    assembly_sources: Mapping[str, AssemblyFastaSource],
    reference_sequences: Mapping[str, str],
    query_path: str,
    target_path: str,
    intervals: Optional[Mapping[str, Sequence[Interval]]] = None,
) -> Dict[str, Tuple[str, int, int]]:
    query_meta: Dict[str, Tuple[str, int, int]] = {}
    requests: DefaultDict[
        str, List[Tuple[str, str, str, int, int]]
    ] = defaultdict(list)
    ordinal = 0
    for candidate in candidates:
        header = header_records.get(candidate.record_id)
        if header is None:
            raise ValueError(
                f"{candidate.partition}.header lacks {candidate.record_id!r}"
            )
        if (
            header.haplotype != candidate.haplotype
            or header.source_contig != candidate.source_contig
            or header.source_start + candidate.local_start
            != candidate.source_start
            or header.source_start + candidate.local_end
            != candidate.source_end
        ):
            raise ValueError(
                f"{candidate.candidate_id}: .header coordinates disagree "
                "with all_unique_regions.bed"
            )
        if candidate.haplotype not in assembly_sources:
            raise ValueError(
                f"--query-paths lacks candidate haplotype "
                f"{candidate.haplotype!r}"
            )
        pieces = (
            list(intervals.get(candidate.candidate_id, ()))
            if intervals is not None else
            [(candidate.local_start, candidate.local_end)]
        )
        for start, end in pieces:
            ordinal += 1
            query_id = f"query_{ordinal:08d}"
            requests[candidate.haplotype].append((
                query_id, candidate.source_contig,
                candidate.candidate_id,
                header.source_start + start,
                header.source_start + end,
            ))
            query_meta[query_id] = (candidate.candidate_id, start, end)

    # A partition may contain records from hundreds of assemblies. Open one
    # assembly at a time so each worker holds only one FASTA descriptor and
    # transiently reads only one .fai, regardless of partition size.
    with open(query_path, "wt") as out:
        for haplotype in sorted(requests):
            source = assembly_sources[haplotype]
            with IndexedFasta(source.fasta_path, source.fai_path) as reader:
                for query_id, contig, _candidate_id, start, end in requests[
                    haplotype
                ]:
                    sequence = reader.fetch(contig, start, end)
                    out.write(f">{query_id}\n{wrap_fasta(sequence)}\n")
    with open(target_path, "wt") as out:
        for template in templates:
            sequence = reference_sequences[template.contig][
                template.window_start:template.window_end
            ]
            out.write(f">{template.template_id}\n{wrap_fasta(sequence)}\n")
    return query_meta


def clipped_hit_pairs(
    hit: AlignmentHit,
    piece_length: int,
) -> List[Tuple[int, int, int, int]]:
    output = []
    for q0, q1, t0, t1 in hit.aligned_pairs:
        left = max(0, int(q0))
        right = min(piece_length, int(q1))
        if right <= left:
            continue
        trim_left = left - int(q0)
        trim_right = int(q1) - right
        if hit.strand == "+":
            nt0, nt1 = int(t0) + trim_left, int(t1) - trim_right
        else:
            nt0, nt1 = int(t0) + trim_right, int(t1) - trim_left
        if nt1 > nt0 and nt1 - nt0 == right - left:
            output.append((left, right, nt0, nt1))
    return output


def qualifying_hit(
    hit: AlignmentHit,
    query_meta: Mapping[str, Tuple[str, int, int, str]],
    targets: Mapping[str, TemplateWindow],
    maximum_template_gap: int,
    minimum_unmasked: int,
) -> Optional[Tuple[str, Tuple[Interval, ...], Tuple[Interval, ...], Tuple[Tuple[int, int, int, int], ...], int]]:
    meta = query_meta.get(hit.query_id)
    target = targets.get(hit.target_id)
    if meta is None or target is None:
        return None
    candidate_id, piece_start, piece_end, sequence = meta
    pairs = clipped_hit_pairs(hit, piece_end - piece_start)
    if not pairs:
        return None
    query_blocks = tuple(merge_intervals([
        (piece_start + q0, piece_start + q1) for q0, q1, _t0, _t1 in pairs
    ]))
    target_blocks = tuple(merge_intervals([
        (target.window_start + t0, target.window_start + t1)
        for _q0, _q1, t0, t1 in pairs
    ]))
    if min(
        interval_distance(block, (target.template_start, target.template_end))
        for block in target_blocks
    ) >= maximum_template_gap:
        return None
    unmasked = 0
    for start, end in merge_intervals([(q0, q1) for q0, q1, _t0, _t1 in pairs]):
        unmasked += count_unmasked(sequence[start:end])
    if unmasked <= minimum_unmasked:
        return None
    paired = tuple(
        (piece_start + q0, piece_start + q1,
         target.window_start + t0, target.window_start + t1)
        for q0, q1, t0, t1 in pairs
    )
    return candidate_id, query_blocks, target_blocks, paired, unmasked


_COORDINATE_QUERY_META: Mapping[str, Tuple[str, int, int]] = {}
_COORDINATE_TARGETS: Mapping[str, TemplateWindow] = {}
_COORDINATE_REFERENCE_SEQUENCES: Mapping[str, str] = {}
_COORDINATE_UNMASKED_RUNS: Mapping[
    str, Tuple[Tuple[int, ...], Tuple[int, ...]]
] = {}
_COORDINATE_PARTITION = ""
_COORDINATE_REFERENCE_HAPLOTYPE = ""
_COORDINATE_MAXIMUM_GAP = 0
_COORDINATE_MINIMUM_UNMASKED = 0
_COORDINATE_MINIMUM_IDENTITY = 0.0


def _initialize_coordinate_paf_worker(
    query_meta: Mapping[str, Tuple[str, int, int]],
    targets: Mapping[str, TemplateWindow],
    reference_sequences: Mapping[str, str],
    partition: str,
    reference_haplotype: str,
    maximum_gap: int,
    minimum_unmasked: int,
    minimum_identity: float,
) -> None:
    global _COORDINATE_QUERY_META
    global _COORDINATE_TARGETS
    global _COORDINATE_REFERENCE_SEQUENCES
    global _COORDINATE_UNMASKED_RUNS
    global _COORDINATE_PARTITION
    global _COORDINATE_REFERENCE_HAPLOTYPE
    global _COORDINATE_MAXIMUM_GAP
    global _COORDINATE_MINIMUM_UNMASKED
    global _COORDINATE_MINIMUM_IDENTITY
    _COORDINATE_QUERY_META = query_meta
    _COORDINATE_TARGETS = targets
    _COORDINATE_REFERENCE_SEQUENCES = reference_sequences
    unmasked_runs: Dict[str, Tuple[Tuple[int, ...], Tuple[int, ...]]] = {}
    for target_id, target in targets.items():
        sequence = reference_sequences[target.contig][
            target.window_start:target.window_end
        ]
        runs = tuple(
            (match.start(), match.end())
            for match in re.finditer(r"[ACGT]+", sequence)
        )
        unmasked_runs[target_id] = (
            tuple(start for start, _end in runs),
            tuple(end for _start, end in runs),
        )
    _COORDINATE_UNMASKED_RUNS = unmasked_runs
    _COORDINATE_PARTITION = partition
    _COORDINATE_REFERENCE_HAPLOTYPE = reference_haplotype
    _COORDINATE_MAXIMUM_GAP = maximum_gap
    _COORDINATE_MINIMUM_UNMASKED = minimum_unmasked
    _COORDINATE_MINIMUM_IDENTITY = minimum_identity


def _coordinate_placement_from_paf_line(
    raw: str,
) -> Optional[AcceptedPlacement]:
    """Accept one PAF row using spans only; no CIGAR is generated or parsed."""
    fields = raw.rstrip("\n").split("\t")
    if len(fields) < 12:
        return None
    query_id, target_id = fields[0], fields[5]
    meta = _COORDINATE_QUERY_META.get(query_id)
    target = _COORDINATE_TARGETS.get(target_id)
    if meta is None or target is None:
        return None
    try:
        query_start, query_end = int(fields[2]), int(fields[3])
        target_start, target_end = int(fields[7]), int(fields[8])
        matches, block_length = int(fields[9]), int(fields[10])
    except ValueError:
        return None
    strand = fields[4]
    if (
        strand not in {"+", "-"}
        or query_end <= query_start
        or target_end <= target_start
        or block_length <= 0
    ):
        return None
    identity = 100.0 * matches / block_length
    if identity < _COORDINATE_MINIMUM_IDENTITY:
        return None
    candidate_id, piece_start, piece_end = meta
    piece_length = piece_end - piece_start
    query_start = max(0, min(query_start, piece_length))
    query_end = max(query_start, min(query_end, piece_length))
    if query_end <= query_start:
        return None
    absolute_target = (
        target.window_start + target_start,
        target.window_start + target_end,
    )
    if interval_distance(
        absolute_target, (target.template_start, target.template_end),
    ) >= _COORDINATE_MAXIMUM_GAP:
        return None
    runs = _COORDINATE_UNMASKED_RUNS.get(target_id)
    if runs is None:
        return None
    # Masking is determined from the aligned reference span. Uppercase runs
    # are indexed once per target, so thousands of secondary hits do not
    # repeatedly scan or copy the same reference bases.
    run_starts, run_ends = runs
    run_index = bisect.bisect_right(run_ends, target_start)
    unmasked = 0
    while run_index < len(run_starts) and run_starts[run_index] < target_end:
        unmasked += max(
            0,
            min(target_end, run_ends[run_index])
            - max(target_start, run_starts[run_index]),
        )
        run_index += 1
    if unmasked <= _COORDINATE_MINIMUM_UNMASKED:
        return None
    absolute_query = (
        piece_start + query_start, piece_start + query_end,
    )
    score = float(matches)
    for tag in fields[12:]:
        if tag.startswith("AS:i:"):
            try:
                score = float(tag[5:])
            except ValueError:
                pass
            break
    paired = ((
        absolute_query[0], absolute_query[1],
        absolute_target[0], absolute_target[1],
    ),)
    return AcceptedPlacement(
        "", candidate_id, _COORDINATE_PARTITION, target_id,
        _COORDINATE_REFERENCE_HAPLOTYPE, target.contig,
        absolute_target[0], absolute_target[1], strand, "minimap2",
        score, identity, unmasked, (absolute_query,), (absolute_target,),
        paired,
    )


def _parse_coordinate_paf_range(
    task: Tuple[str, int, int],
) -> Tuple[int, List[AcceptedPlacement]]:
    paf_path, range_start, range_end = task
    output: List[AcceptedPlacement] = []
    lines = 0
    with open(paf_path, "rb") as handle:
        if range_start:
            handle.seek(range_start - 1)
            if handle.read(1) != b"\n":
                handle.readline()
        else:
            handle.seek(0)
        while handle.tell() < range_end:
            raw = handle.readline()
            if not raw:
                break
            lines += 1
            if not raw.strip() or raw.startswith(b"#"):
                continue
            placement = _coordinate_placement_from_paf_line(
                raw.decode("utf-8"),
            )
            if placement is not None:
                output.append(placement)
    return lines, output


def parse_coordinate_paf(
    paf_path: str,
    query_meta: Mapping[str, Tuple[str, int, int]],
    targets: Mapping[str, TemplateWindow],
    reference_sequences: Mapping[str, str],
    partition: str,
    reference_haplotype: str,
    maximum_gap: int,
    minimum_unmasked: int,
    minimum_identity: float,
) -> List[AcceptedPlacement]:
    """Parse and qualify one partition PAF in its owning worker process."""
    size = os.path.getsize(paf_path)
    if size == 0:
        return []
    _initialize_coordinate_paf_worker(
        query_meta, targets, reference_sequences, partition,
        reference_haplotype, maximum_gap, minimum_unmasked,
        minimum_identity,
    )
    processed, accepted = _parse_coordinate_paf_range((paf_path, 0, size))
    LOG.info(
        "%s: coordinate-only PAF processing scanned %d lines; "
        "%d qualifying placements",
        partition, processed, len(accepted),
    )
    return accepted


def fingerprint(path: str) -> Tuple[str, int, int]:
    stat = os.stat(path)
    return os.path.abspath(path), stat.st_size, stat.st_mtime_ns


def read_partition_list(path: str) -> List[str]:
    """Read exact partition names, accepting either names or folder paths."""
    names: List[str] = []
    seen = set()
    with open(path, "rt") as handle:
        for line_number, raw in enumerate(handle, 1):
            stripped = raw.strip()
            if not stripped or stripped.startswith("#"):
                continue
            value = stripped.rstrip("/")
            name = os.path.basename(value)
            if not name or name in {".", ".."} or os.path.basename(name) != name:
                raise ValueError(
                    f"{path}:{line_number}: unsafe partition entry {raw.strip()!r}"
                )
            if name in seen:
                raise ValueError(
                    f"{path}:{line_number}: duplicate partition {name!r}"
                )
            seen.add(name)
            names.append(name)
    if not names:
        raise ValueError(f"{path}: no partition names")
    return names


def partition_result_signature(
    signature: str,
    partition: str,
    candidates: Sequence[UnifiedNovel],
    template_coords: Sequence[Tuple[str, int, int]],
) -> str:
    return hashlib.sha256(json.dumps({
        "signature": signature,
        "partition": partition,
        "candidates": [dataclasses.asdict(row) for row in candidates],
        "templates": list(template_coords),
    }, sort_keys=True).encode()).hexdigest()


def run_partition(
    partition: str,
    candidates: Sequence[UnifiedNovel],
    template_coords: Sequence[Tuple[str, int, int]],
    partition_files: PartitionFiles,
    assembly_sources: Mapping[str, AssemblyFastaSource],
    reference_source: AssemblyFastaSource,
    reference_sequences: Mapping[str, str],
    work_root: str,
    args: argparse.Namespace,
    signature: str,
    return_result: bool = True,
) -> Optional[PartitionResult]:
    workdir = os.path.join(work_root, "partitions", partition)
    result_path = os.path.join(workdir, "accepted_placements.tsv")
    residual_path = os.path.join(workdir, "residuals.tsv")
    template_path = os.path.join(workdir, "template_windows.tsv")
    marker_path = os.path.join(workdir, "_SUCCESS.json")
    partition_signature = partition_result_signature(
        signature, partition, candidates, template_coords,
    )
    if args.resume and os.path.isfile(marker_path):
        try:
            saved = json.load(open(marker_path))
            if saved.get("signature") == partition_signature and all(
                os.path.isfile(path) for path in
                (result_path, residual_path, template_path)
            ):
                if return_result:
                    return read_partition_result(
                        partition, candidates, result_path, residual_path,
                        template_path,
                    )
                return None
        except (OSError, ValueError, AttributeError):
            pass
    if os.path.isdir(workdir):
        # On NFS/Lustre, recently closed files can transiently remain as
        # hidden directory entries.  Stale scratch is safe because every
        # candidate FASTA/PAF is overwritten below.
        shutil.rmtree(workdir, ignore_errors=True)
    Path(workdir).mkdir(parents=True, exist_ok=True)
    lengths = {
        name: len(sequence) for name, sequence in reference_sequences.items()
    }
    try:
        needed_records = {row.record_id for row in candidates}
        header_records = {
            row.record_id: row
            for row in parse_header_records(partition_files, needed_records)
        }
        missing_records = needed_records - set(header_records)
        if missing_records:
            raise ValueError(
                f"{partition}.header lacks {len(missing_records)} candidate "
                "records; examples: " + ",".join(sorted(missing_records)[:10])
            )
        windows = build_template_windows(
            partition, template_coords,
            max(row.length for row in candidates), lengths,
            args.minimum_flank,
        )
        windows = [
            dataclasses.replace(row, haplotype=reference_source.haplotype)
            for row in windows
        ]
        target_by_id = {row.template_id: row for row in windows}
        stage = os.path.join(workdir, "combined_alignment")
        Path(stage).mkdir(parents=True, exist_ok=True)
        query_fasta = os.path.join(stage, "all_novels.fa")
        target_fasta = os.path.join(stage, "expanded_templates.fa")
        query_meta = write_partition_fastas(
            candidates, windows, header_records, assembly_sources,
            reference_sequences,
            query_fasta, target_fasta,
        )
        paf_path = os.path.join(stage, "nearby_template.minimap2.paf")
        run_minimap2(
            query_fasta, target_fasta, paf_path, args.cores, LOG,
            params={**MINIMAP_PARAMS, "N": args.max_secondary},
            emit_cigar=False,
        )
        LOG.info(
            "%s: parsing coordinate-only PAF in its partition worker "
            "(%.2f GiB)",
            partition,
            os.path.getsize(paf_path) / (1024 ** 3),
        )
        accepted = parse_coordinate_paf(
            paf_path, query_meta, target_by_id, reference_sequences, partition,
            reference_source.haplotype, args.maximum_template_gap,
            args.minimum_unmasked, args.min_identity,
        )
        new_masks: DefaultDict[str, List[Interval]] = defaultdict(list)
        seen = set()
        unique_accepted: List[AcceptedPlacement] = []
        accepted.sort(key=lambda row: (
            row.candidate_id, row.ref_contig, row.ref_start, row.ref_end,
            row.query_blocks, row.strand, row.template_id,
            -row.alignment_score, -row.identity,
        ))
        for placement in accepted:
            identity = (
                placement.candidate_id, placement.template_id,
                placement.strand, placement.paired_blocks,
            )
            if identity in seen:
                continue
            seen.add(identity)
            new_masks[placement.candidate_id].extend(placement.query_blocks)
            unique_accepted.append(placement)
        accepted = unique_accepted
        residuals: Dict[str, List[Interval]] = {
            candidate.candidate_id: subtract_intervals(
                [(candidate.local_start, candidate.local_end)],
                merge_intervals(new_masks.get(candidate.candidate_id, ())),
            )
            for candidate in candidates
        }
        shutil.rmtree(stage, ignore_errors=True)
    finally:
        # Successful partitions only need the compact result tables written
        # below. Best-effort cleanup must never turn completed alignment work
        # into an error on a network filesystem.
        shutil.rmtree(
            os.path.join(workdir, "combined_alignment"),
            ignore_errors=True,
        )
    accepted = [
        dataclasses.replace(
            row,
            placement_id=f"{partition}.placement_{index:012d}",
        )
        for index, row in enumerate(accepted, 1)
    ]
    residual_rows = tuple(
        (candidate_id, start, end)
        for candidate_id in sorted(residuals)
        for start, end in residuals[candidate_id]
    )
    write_partition_outputs(
        result_path, residual_path, template_path,
        accepted, residual_rows, windows,
    )
    with open(marker_path + ".tmp", "wt") as out:
        json.dump({
            "signature": partition_signature,
            "candidates": len(candidates),
            "placements": len(accepted),
            "residuals": len(residual_rows),
        }, out, sort_keys=True)
        out.write("\n")
    os.replace(marker_path + ".tmp", marker_path)
    if return_result:
        return PartitionResult(
            partition, tuple(candidates), tuple(accepted), tuple(windows),
            residual_rows,
        )
    return None


_PARTITION_FILES: Mapping[str, PartitionFiles] = {}
_PARTITION_CANDIDATES: Mapping[str, Sequence[UnifiedNovel]] = {}
_PARTITION_TEMPLATES: Mapping[
    str, Sequence[Tuple[str, int, int]]
] = {}
_PARTITION_ASSEMBLIES: Mapping[str, AssemblyFastaSource] = {}
_PARTITION_REFERENCE_SOURCE: Optional[AssemblyFastaSource] = None
_PARTITION_REFERENCE_SEQUENCES: Mapping[str, str] = {}
_PARTITION_WORK_ROOT = ""
_PARTITION_ARGS: Optional[argparse.Namespace] = None
_PARTITION_SIGNATURE = ""


def _initialize_partition_process(
    partition_files: Mapping[str, PartitionFiles],
    candidates: Mapping[str, Sequence[UnifiedNovel]],
    templates: Mapping[str, Sequence[Tuple[str, int, int]]],
    assemblies: Mapping[str, AssemblyFastaSource],
    reference_source: AssemblyFastaSource,
    reference_sequences: Mapping[str, str],
    work_root: str,
    args: argparse.Namespace,
    signature: str,
) -> None:
    global _PARTITION_FILES
    global _PARTITION_CANDIDATES
    global _PARTITION_TEMPLATES
    global _PARTITION_ASSEMBLIES
    global _PARTITION_REFERENCE_SOURCE
    global _PARTITION_REFERENCE_SEQUENCES
    global _PARTITION_WORK_ROOT
    global _PARTITION_ARGS
    global _PARTITION_SIGNATURE
    _PARTITION_FILES = partition_files
    _PARTITION_CANDIDATES = candidates
    _PARTITION_TEMPLATES = templates
    _PARTITION_ASSEMBLIES = assemblies
    _PARTITION_REFERENCE_SOURCE = reference_source
    _PARTITION_REFERENCE_SEQUENCES = reference_sequences
    _PARTITION_WORK_ROOT = work_root
    _PARTITION_ARGS = args
    _PARTITION_SIGNATURE = signature


def _run_partition_process(partition: str) -> str:
    if _PARTITION_REFERENCE_SOURCE is None or _PARTITION_ARGS is None:
        raise RuntimeError("partition process was not initialized")
    run_partition(
        partition,
        _PARTITION_CANDIDATES[partition],
        _PARTITION_TEMPLATES[partition],
        _PARTITION_FILES[partition],
        _PARTITION_ASSEMBLIES,
        _PARTITION_REFERENCE_SOURCE,
        _PARTITION_REFERENCE_SEQUENCES,
        _PARTITION_WORK_ROOT,
        _PARTITION_ARGS,
        _PARTITION_SIGNATURE,
        return_result=False,
    )
    return partition


def run_partitions_multiprocess(
    active: Sequence[str],
    partitions: Mapping[str, PartitionFiles],
    candidates: Mapping[str, Sequence[UnifiedNovel]],
    templates: Mapping[str, Sequence[Tuple[str, int, int]]],
    assemblies: Mapping[str, AssemblyFastaSource],
    reference_source: AssemblyFastaSource,
    reference_sequences: Mapping[str, str],
    work_root: str,
    args: argparse.Namespace,
    signature: str,
) -> List[str]:
    if not active:
        return []
    worker_count = min(args.jobs, len(active))
    context = mp_context()
    LOG.info(
        "Processing %d partitions with %d persistent Python processes; "
        "each process runs Minimap2 with %d threads",
        len(active), worker_count, args.cores,
    )
    completed: List[str] = []
    progress_step = max(1, len(active) // 20)
    with context.Pool(
        processes=worker_count,
        initializer=_initialize_partition_process,
        initargs=(
            partitions, candidates, templates, assemblies, reference_source,
            reference_sequences, work_root, args, signature,
        ),
    ) as pool:
        for partition in pool.imap_unordered(
            _run_partition_process, active, chunksize=1,
        ):
            completed.append(partition)
            if len(completed) == len(active) or len(completed) % progress_step == 0:
                LOG.info(
                    "Partition progress: %d/%d complete",
                    len(completed), len(active),
                )
    return completed


PLACEMENT_FIELDS = (
    "placement_id", "candidate_id", "partition", "template_id",
    "ref_haplotype", "ref_contig", "ref_start", "ref_end", "strand",
    "aligner", "alignment_score", "identity", "unmasked_aligned_bases",
    "query_blocks", "target_blocks", "paired_blocks",
)
TEMPLATE_FIELDS = (
    "template_id", "partition", "haplotype", "contig",
    "template_start", "template_end", "window_start", "window_end",
)


def placement_row(row: AcceptedPlacement) -> Tuple[object, ...]:
    return (
        row.placement_id, row.candidate_id, row.partition, row.template_id,
        row.ref_haplotype, row.ref_contig, row.ref_start, row.ref_end,
        row.strand, row.aligner, row.alignment_score, row.identity,
        row.unmasked_aligned_bases,
        json.dumps(row.query_blocks, separators=(",", ":")),
        json.dumps(row.target_blocks, separators=(",", ":")),
        json.dumps(row.paired_blocks, separators=(",", ":")),
    )


def write_partition_outputs(
    placement_path: str,
    residual_path: str,
    template_path: str,
    placements: Sequence[AcceptedPlacement],
    residuals: Sequence[Tuple[str, int, int]],
    templates: Sequence[TemplateWindow],
) -> None:
    write_tsv(placement_path, PLACEMENT_FIELDS, map(placement_row, placements))
    write_tsv(
        residual_path, ("candidate_id", "local_start", "local_end"),
        residuals,
    )
    write_tsv(
        template_path, TEMPLATE_FIELDS,
        (dataclasses.astuple(row) for row in templates),
    )


def read_partition_result(
    partition: str,
    candidates: Sequence[UnifiedNovel],
    placement_path: str,
    residual_path: str,
    template_path: str,
) -> PartitionResult:
    placements = []
    with open(placement_path, newline="") as handle:
        for row in csv.DictReader(handle, delimiter="\t"):
            placements.append(AcceptedPlacement(
                row["placement_id"], row["candidate_id"], row["partition"],
                row["template_id"], row["ref_haplotype"], row["ref_contig"],
                int(row["ref_start"]), int(row["ref_end"]), row["strand"],
                row["aligner"], float(row["alignment_score"]),
                float(row["identity"]), int(row["unmasked_aligned_bases"]),
                tuple(map(tuple, json.loads(row["query_blocks"]))),
                tuple(map(tuple, json.loads(row["target_blocks"]))),
                tuple(map(tuple, json.loads(row["paired_blocks"]))),
            ))
    with open(residual_path, newline="") as handle:
        residuals = tuple(
            (row["candidate_id"], int(row["local_start"]), int(row["local_end"]))
            for row in csv.DictReader(handle, delimiter="\t")
        )
    with open(template_path, newline="") as handle:
        templates = tuple(
            TemplateWindow(
                row["template_id"], row["partition"], row["haplotype"],
                row["contig"], int(row["template_start"]),
                int(row["template_end"]), int(row["window_start"]),
                int(row["window_end"]),
            )
            for row in csv.DictReader(handle, delimiter="\t")
        )
    return PartitionResult(
        partition, tuple(candidates), tuple(placements), templates, residuals,
    )


def read_completed_partition_result(
    partition: str,
    candidates: Sequence[UnifiedNovel],
    template_coords: Sequence[Tuple[str, int, int]],
    work_root: str,
    signature: str,
) -> Optional[PartitionResult]:
    """Read one shard result only when its exact success signature matches."""
    workdir = os.path.join(work_root, "partitions", partition)
    result_path = os.path.join(workdir, "accepted_placements.tsv")
    residual_path = os.path.join(workdir, "residuals.tsv")
    template_path = os.path.join(workdir, "template_windows.tsv")
    marker_path = os.path.join(workdir, "_SUCCESS.json")
    if not all(os.path.isfile(path) for path in (
        result_path, residual_path, template_path, marker_path,
    )):
        return None
    expected = partition_result_signature(
        signature, partition, candidates, template_coords,
    )
    try:
        with open(marker_path, "rt") as handle:
            marker = json.load(handle)
        if marker.get("signature") != expected:
            return None
        return read_partition_result(
            partition, candidates, result_path, residual_path, template_path,
        )
    except (OSError, ValueError, KeyError, TypeError, json.JSONDecodeError):
        return None


def write_public_outputs(
    output_dir: str,
    results: Sequence[PartitionResult],
) -> None:
    candidates = [row for result in results for row in result.candidates]
    placements = [row for result in results for row in result.placements]
    residuals = [row for result in results for row in result.residuals]
    templates = [row for result in results for row in result.templates]
    candidate_by_id = {row.candidate_id: row for row in candidates}

    write_tsv(
        os.path.join(output_dir, "all_novel_candidates.tsv"),
        tuple(field.name for field in dataclasses.fields(UnifiedNovel)),
        (dataclasses.astuple(row) for row in candidates),
    )
    write_tsv(
        os.path.join(output_dir, "all_novel_template_alignments.tsv"),
        PLACEMENT_FIELDS, map(placement_row, placements),
    )
    write_tsv(
        os.path.join(output_dir, "residual_novel_intervals.tsv"),
        ("candidate_id", "partition", "record_id", "haplotype",
         "source_contig", "source_start", "source_end", "local_start",
         "local_end"),
        (
            (
                candidate_id, candidate_by_id[candidate_id].partition,
                candidate_by_id[candidate_id].record_id,
                candidate_by_id[candidate_id].haplotype,
                candidate_by_id[candidate_id].source_contig,
                candidate_by_id[candidate_id].source_start
                + start - candidate_by_id[candidate_id].local_start,
                candidate_by_id[candidate_id].source_start
                + end - candidate_by_id[candidate_id].local_start,
                start, end,
            )
            for candidate_id, start, end in residuals
        ),
    )
    masks: DefaultDict[str, List[Interval]] = defaultdict(list)
    placements_by_candidate: DefaultDict[
        str, List[AcceptedPlacement]
    ] = defaultdict(list)
    for placement in placements:
        masks[placement.candidate_id].extend(placement.query_blocks)
        placements_by_candidate[placement.candidate_id].append(placement)
    covered_rows = []
    for candidate_id, intervals in masks.items():
        candidate = candidate_by_id[candidate_id]
        for start, end in merge_intervals(intervals):
            ids = sorted({
                placement.placement_id
                for placement in placements_by_candidate[candidate_id]
                if any(
                    a < end and start < b
                    for a, b in placement.query_blocks
                )
            })
            covered_rows.append((
                candidate.partition, candidate_id, candidate.record_id,
                candidate.haplotype, candidate.source_contig,
                candidate.source_start + start - candidate.local_start,
                candidate.source_start + end - candidate.local_start,
                start, end, ",".join(ids),
            ))
    write_tsv(
        os.path.join(output_dir, "reference_covered_novel.tsv"),
        ("partition", "candidate_id", "record_id", "haplotype",
         "source_contig", "source_start", "source_end", "local_start",
         "local_end", "placement_ids"),
        covered_rows,
    )
    placement_by_template: DefaultDict[Tuple[str, str], List[AcceptedPlacement]] = defaultdict(list)
    for placement in placements:
        placement_by_template[(placement.partition, placement.template_id)].append(
            placement
        )
    provisional_refined: DefaultDict[
        Tuple[str, str, str], List[Tuple[int, int, str]]
    ] = defaultdict(list)
    for template in templates:
        members = placement_by_template.get(
            (template.partition, template.template_id), ()
        )
        starts = [template.template_start]
        ends = [template.template_end]
        for placement in members:
            starts.extend(start for start, _end in placement.target_blocks)
            ends.extend(end for _start, end in placement.target_blocks)
        provisional_refined[(
            template.partition, template.haplotype, template.contig,
        )].append((
            min(starts), max(ends), "expanded" if members else "unchanged",
        ))
    refined_rows = []
    for (partition, haplotype, contig), spans in sorted(
        provisional_refined.items()
    ):
        merged = merge_intervals((start, end) for start, end, _status in spans)
        for index, (start, end) in enumerate(merged, 1):
            status = "expanded" if any(
                state == "expanded" and a < end and start < b
                for a, b, state in spans
            ) else "unchanged"
            refined_rows.append((
                partition, index, haplotype, contig, start, end, status,
            ))
    write_tsv(
        os.path.join(output_dir, "refined_template_intervals.tsv"),
        ("partition", "template_index", "haplotype", "contig", "start",
         "end", "status"),
        refined_rows,
    )


def write_success_marker(
    marker: str,
    signature: str,
    results: Sequence[PartitionResult],
    query_paths: str,
    main_reference_fasta: str,
) -> None:
    temporary = marker + f".tmp.{os.getpid()}"
    try:
        with open(temporary, "wt") as out:
            json.dump({
                "version": VERSION,
                "signature": signature,
                "partitions": len(results),
                "candidates": sum(len(row.candidates) for row in results),
                "placements": sum(len(row.placements) for row in results),
                "residuals": sum(len(row.residuals) for row in results),
                "query_paths": query_paths,
                "main_reference": main_reference_fasta,
            }, out, sort_keys=True)
            out.write("\n")
            out.flush()
            os.fsync(out.fileno())
        os.replace(temporary, marker)
    finally:
        try:
            os.remove(temporary)
        except FileNotFoundError:
            pass


def parse_args(argv: Optional[Sequence[str]] = None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("-a", "--analysis-dir", required=True)
    parser.add_argument("-q", "--query-paths", required=True)
    parser.add_argument("-r", "--main-reference", required=True)
    parser.add_argument("-o", "--output-dir", required=True)
    parser.add_argument(
        "-l", "--partition-list",
        help=(
            "run only listed partition names/folder paths and write only "
            "per-partition checkpoint results; use --finalize-only after all "
            "disjoint shard jobs finish"
        ),
    )
    parser.add_argument(
        "--finalize-only", action="store_true",
        help=(
            "run no alignments; validate all per-partition shard checkpoints "
            "and create the aggregate unified-refinement TSVs"
        ),
    )
    parser.add_argument("-j", "--jobs", type=int, default=16)
    parser.add_argument("-c", "--cores", type=int, default=4)
    parser.add_argument("--minimum-flank", type=int, default=5000)
    parser.add_argument("--maximum-template-gap", type=int, default=100)
    parser.add_argument("--minimum-unmasked", type=int, default=100)
    parser.add_argument(
        "--min-identity", type=float, default=0.0,
        help=(
            "optional extra identity filter; default 0 applies only the "
            "requested gap and unmasked-base acceptance rules"
        ),
    )
    parser.add_argument(
        "--max-secondary", type=int, default=1000,
        help=(
            "Minimap2 -N value; deliberately high so multiple qualifying "
            "placements can all refine a template (default: 1000)"
        ),
    )
    parser.add_argument("--resume", action="store_true")
    parser.add_argument(
        "--log-level", choices=("DEBUG", "INFO", "WARNING", "ERROR"),
        default="INFO",
    )
    args = parser.parse_args(argv)
    if args.partition_list and args.finalize_only:
        parser.error("--partition-list and --finalize-only are mutually exclusive")
    if args.jobs < 1 or args.cores < 1 or args.max_secondary < 1:
        parser.error("--jobs, --cores, and --max-secondary must be positive")
    if args.minimum_flank < 0 or args.maximum_template_gap < 1:
        parser.error("flank/gap parameters are invalid")
    if args.minimum_unmasked < 0:
        parser.error("--minimum-unmasked cannot be negative")
    if not 0 <= args.min_identity <= 100:
        parser.error("--min-identity must be between 0 and 100")
    return args


def run(args: argparse.Namespace) -> None:
    analysis_dir = os.path.abspath(args.analysis_dir)
    output_dir = os.path.abspath(args.output_dir)
    query_paths = os.path.abspath(args.query_paths)
    main_reference_fasta = os.path.abspath(args.main_reference)
    Path(output_dir).mkdir(parents=True, exist_ok=True)
    if not args.finalize_only:
        ensure_executable("minimap2")
    partitions = read_partition_inputs(
        os.path.join(analysis_dir, "partition_inputs.tsv")
    )
    unique = read_unique_rows(analysis_dir, partitions)
    used = read_used_reference_rows(analysis_dir)
    reference_records = parse_main_reference(main_reference_fasta)
    main_haplotypes = {row.haplotype for row in reference_records.values()}
    analysis_marker = os.path.join(
        analysis_dir, "_COORDINATE_ANALYSIS_SUCCESS.json"
    )
    metadata = json.load(open(analysis_marker))
    reference_haplotypes = set(metadata.get("reference_haplotypes", ()))
    primary_haplotypes = main_haplotypes & reference_haplotypes
    if len(primary_haplotypes) != 1:
        raise ValueError(
            "--main-reference must contain exactly one primary source "
            "haplotype named by coordinate-analysis metadata; found catalog="
            f"{sorted(main_haplotypes)}, primary={sorted(reference_haplotypes)}"
        )
    reference_haplotype = next(iter(primary_haplotypes))
    sources = read_assembly_sources(query_paths)
    if reference_haplotype not in sources:
        raise ValueError(
            f"--query-paths lacks main reference {reference_haplotype!r}"
        )
    reference_haplotypes = metadata.get(
        "reference_haplotypes", [reference_haplotype]
    )
    candidates = unified_candidates(unique, reference_haplotypes)
    by_partition: DefaultDict[str, List[UnifiedNovel]] = defaultdict(list)
    for candidate in candidates:
        by_partition[candidate.partition].append(candidate)
    from fixed_alternatives import fixed_partitions
    fixed = fixed_partitions(partitions)
    # Imported alternatives can be alignment targets; their source genome
    # is not a template that may be expanded by refinement.
    templates = merge_template_rows(
        [row for row in used if row.encoded_name not in fixed.get(row.partition, {})],
        reference_haplotype)
    active_all = sorted(
        partition for partition, rows in by_partition.items()
        if rows and templates.get(partition)
    )
    listed_names: Optional[List[str]] = None
    if args.partition_list:
        listed_names = read_partition_list(os.path.abspath(args.partition_list))
        known = set(partitions)
        unknown = [name for name in listed_names if name not in known]
        if unknown:
            LOG.warning(
                "Ignoring %d listed folders absent from partition_inputs.tsv "
                "(examples: %s)",
                len(unknown), ",".join(unknown[:5]),
            )
        selected = set(listed_names) & known
        active = [name for name in active_all if name in selected]
    else:
        active = list(active_all)
    required_reference_contigs = sorted({
        contig
        for partition in active
        for contig, _start, _end in templates[partition]
    })
    required = [
        analysis_marker,
        os.path.join(analysis_dir, "all_unique_regions.bed"),
        os.path.join(analysis_dir, "all_used_novel_loci.bed"),
        query_paths, main_reference_fasta,
        sources[reference_haplotype].fasta_path,
        sources[reference_haplotype].fai_path,
    ]
    signature = hashlib.sha256(json.dumps({
        "version": VERSION,
        "files": [fingerprint(path) for path in required],
        "fixed_templates": fixed,
        "minimum_flank": args.minimum_flank,
        "maximum_template_gap": args.maximum_template_gap,
        "minimum_unmasked": args.minimum_unmasked,
        "min_identity": args.min_identity,
        "max_secondary": args.max_secondary,
        "aligner": "minimap2",
    }, sort_keys=True).encode()).hexdigest()
    marker = os.path.join(
        output_dir, "_UNIFIED_TEMPLATE_REFINEMENT_SUCCESS.json",
    )
    public_outputs = (
        "all_novel_candidates.tsv",
        "all_novel_template_alignments.tsv",
        "residual_novel_intervals.tsv",
        "reference_covered_novel.tsv",
        "refined_template_intervals.tsv",
    )
    if (
        not args.partition_list and not args.finalize_only
        and args.resume and os.path.isfile(marker)
    ):
        try:
            saved = json.load(open(marker))
            if saved.get("signature") == signature and all(
                os.path.isfile(os.path.join(output_dir, name))
                for name in public_outputs
            ):
                LOG.info(
                    "RESUME: unified template refinement is complete in %s",
                    output_dir,
                )
                return
        except (OSError, ValueError, AttributeError):
            pass

    work_root = os.path.join(output_dir, "alignment_work")

    # A shard or finalization attempt invalidates an older aggregate success
    # marker. Public TSVs are rebuilt only after all exact checkpoints exist.
    if args.partition_list or args.finalize_only:
        try:
            os.remove(marker)
        except FileNotFoundError:
            pass

    if args.finalize_only:
        results: List[PartitionResult] = []
        incomplete: List[str] = []
        for partition in active_all:
            result = read_completed_partition_result(
                partition, by_partition[partition], templates[partition],
                work_root, signature,
            )
            if result is None:
                incomplete.append(partition)
            else:
                results.append(result)
        if incomplete:
            raise RuntimeError(
                f"cannot finalize: {len(incomplete)} active partitions lack "
                "a current successful shard checkpoint (examples: "
                + ",".join(incomplete[:10])
                + "); rerun their shard lists with --resume"
            )
        # Partitions without a main-reference template need no alignment; all
        # candidates in them remain residual.
        for partition in sorted(set(by_partition) - set(active_all)):
            rows = tuple(by_partition[partition])
            results.append(PartitionResult(
                partition, rows, (), (), tuple(
                    (row.candidate_id, row.local_start, row.local_end)
                    for row in rows
                ),
            ))
        results.sort(key=lambda row: row.partition)
        write_public_outputs(output_dir, results)
        write_success_marker(
            marker, signature, results, query_paths, main_reference_fasta,
        )
        LOG.info(
            "FINALIZE: wrote aggregate outputs for %d partitions", len(results),
        )
        return

    reference_sequences: Dict[str, str] = {}
    with IndexedFasta(
        sources[reference_haplotype].fasta_path,
        sources[reference_haplotype].fai_path,
    ) as reference:
        available = set(reference.names())
        missing = set(required_reference_contigs) - available
        if missing:
            raise ValueError(
                "main-reference assembly lacks template contigs: "
                + ",".join(sorted(missing))
            )
        LOG.info(
            "Loading %d required %s reference contigs into shared RAM",
            len(required_reference_contigs), reference_haplotype,
        )
        loaded_bases = 0
        for index, contig in enumerate(required_reference_contigs, 1):
            sequence = reference.sequence(contig)
            reference_sequences[contig] = sequence
            loaded_bases += len(sequence)
            if index == len(required_reference_contigs) or index % 5 == 0:
                LOG.info(
                    "Reference RAM load: %d/%d contigs (%.2f GiB)",
                    index, len(required_reference_contigs),
                    loaded_bases / (1024 ** 3),
                )
    completed = run_partitions_multiprocess(
        active, partitions, by_partition, templates, sources,
        sources[reference_haplotype], reference_sequences,
        work_root, args, signature,
    )
    if args.partition_list:
        assert listed_names is not None
        LOG.info(
            "SHARD COMPLETE: %d listed folders, %d known partitions, %d "
            "alignment checkpoints completed/reused. Run --finalize-only "
            "after every shard job succeeds.",
            len(listed_names),
            len(set(listed_names) & set(partitions)),
            len(completed),
        )
        return

    results: List[PartitionResult] = []
    for partition in active:
        result = read_completed_partition_result(
            partition, by_partition[partition], templates[partition],
            work_root, signature,
        )
        if result is None:
            raise RuntimeError(
                f"{partition}: completed worker produced no valid checkpoint"
            )
        results.append(result)

    # Partitions without a main-reference template retain all candidates.
    for partition in sorted(set(by_partition) - set(active)):
        rows = tuple(by_partition[partition])
        results.append(PartitionResult(
            partition, rows, (), (), tuple(
                (row.candidate_id, row.local_start, row.local_end)
                for row in rows
            ),
        ))
    results.sort(key=lambda row: row.partition)
    write_public_outputs(output_dir, results)
    write_success_marker(
        marker, signature, results, query_paths, main_reference_fasta,
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
