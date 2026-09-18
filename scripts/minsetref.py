#!/usr/bin/env python3
"""Min-set reference builder v3 orchestrator.

Later assemblies are first processed independently against one immutable main
reference. Candidates passing only the base min-segment rule are compared
genome-wide, weighted-gap-remerged on their source contigs, then lifted with
10-kb source anchors and size-filtered only at the end before being finalized as
a separately annotated ``novel_loci.fa``. Every original occurrence—including failed or redundant
promotion candidates—is then independently aligned to those novel loci.
Remaining sequence is cleaned only with insertions sharing an overlapping lift
window on either a main chromosome or novel-locus contig.
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
import subprocess
import sys
import threading
from collections import deque
from concurrent.futures import ThreadPoolExecutor, as_completed
from pathlib import Path
from typing import Callable, Dict, List, Optional, Sequence, Set, Tuple

from minsetref_core import (
    IndexedFasta,
    Piece,
    iter_n_free_segments,
    mp_context,
    read_primary_fasta_names,
    sanitize_id,
    source_interval_id,
    wrap_fasta,
    write_fasta_pieces,
)
from minsetref_align import (
    build_minimap2_index,
    build_winnowmap_rep_kmers,
    collect_alignment_hits,
    coverage_by_query_id,
    ensure_executable,
)
from minsetref_segments import (
    count_masked,
    count_unmasked,
    find_unaligned_segments,
    merge_intervals,
    subtract_intervals,
)
from minsetref_selfclean import self_clean_converge
from minsetref_localize import (
    AnchorBlock,
    LocalizeResult,
    Placement,
    anchor_blocks_from_hits,
    localize_insertion,
    localize_windows_from_blocks,
)

WINNOW_PARAMS = dict(k=19, w=10, m=300, p=0.001, N=100)
MINIMAP_PARAMS = dict(
    preset="asm5", p=0.001, N=20, f=0.001, K="100M",
    retry_on_sigkill=True, retry_threads=32, retry_K="100M",
)
SELF_MINIMAP_PARAMS = dict(
    MINIMAP_PARAMS,
    p=0.001,
    K="100M",
    retry_threads=16,
    retry_K="100M",
)
CHECKPOINT_VERSION = 3
LIFT_BLOCKS_NAME = "lift_blocks.tsv"
STAGE_LOGIC_VERSIONS = {
    # Version 1 remains compatible with checkpoints written before explicit
    # per-stage versioning was introduced.
    # Version 2 aligns every query assembly to the complete first assembly,
    # rather than its self-cleaned/truncated main_chroms.fa representation.
    "main_discovery": 2,
    # Version 2 removes inheritance of a stale pre-trimming main placement.
    # Version 9 consumes coordinate-named, freshly lifted novel loci.
    # Version 10 supports separating unprotected unlifted representatives from
    # the alignment target without changing completed alignment checkpoints.
    "novel_alignment": 10,
}
NOVEL_CLEANUP_LOGIC_VERSION = 12
LOCAL_CLEANUP_LOGIC_VERSION = 10
NOVEL_INPUT_SHARD_VERSION = 2
DEFAULT_PRIOR = "CHM13_h1,HG38_h1,HG002"
DEFAULT_KEEP_UNMAPPED_SAMPLES = ("HG38_h1",)
PRIOR_SCORE_STEP = 10_000_000_000
NOVEL_ANCHOR_BP = 10_000
NOVEL_MERGE_WEIGHT_THRESHOLD = 500.0
NOVEL_MERGE_MASKED_WEIGHT = 0.1
FILTERED_REFERENCE_FASTA = "filtered_sequences.fa"
UNMAPPED_NOVEL_FASTA = "unmapped_novel.fa"
SOURCE_COORD_RE = re.compile(r"^(.*):([^:]+):(\d+)-(\d+)([+-])$")
CANDIDATE_TABLE_FIELDS = [
    "candidate_id", "sample", "contig", "start", "end", "strand",
    "length", "unmasked", "final_id", "parent_id", "note", "origin_id",
    "lift_id", "anchor_status", "anchors_json", "anchor_note",
]


@dataclasses.dataclass(frozen=True)
class BlacklistContig:
    """One contig's merged BED intervals plus arrays for logarithmic lookup."""

    intervals: Tuple[Tuple[int, int], ...]
    starts: Tuple[int, ...]
    ends: Tuple[int, ...]


BlacklistIndex = Dict[str, BlacklistContig]


def file_sha256(path: str) -> str:
    digest = hashlib.sha256()
    with open(path, "rb") as handle:
        while chunk := handle.read(1024 * 1024):
            digest.update(chunk)
    return digest.hexdigest()


def read_blacklist_bed(path: str) -> BlacklistIndex:
    """Read and merge a standard 0-based, half-open BED3+ blacklist."""
    by_contig: Dict[str, List[Tuple[int, int]]] = {}
    with open(path) as handle:
        for line_number, raw in enumerate(handle, 1):
            line = raw.strip()
            if not line or line.startswith(("#", "track ", "browser ")):
                continue
            fields = line.split("\t")
            if len(fields) < 3:
                fields = line.split()
            if len(fields) < 3:
                raise ValueError(
                    f"blacklist BED line {line_number}: expected at least 3 columns"
                )
            contig = fields[0]
            try:
                start, end = int(fields[1]), int(fields[2])
            except ValueError as error:
                raise ValueError(
                    f"blacklist BED line {line_number}: start/end must be integers"
                ) from error
            if start < 0 or end <= start:
                raise ValueError(
                    f"blacklist BED line {line_number}: require 0 <= start < end, "
                    f"got {start}, {end}"
                )
            by_contig.setdefault(contig, []).append((start, end))

    output: BlacklistIndex = {}
    for contig, values in by_contig.items():
        merged = tuple(merge_intervals(values))
        output[contig] = BlacklistContig(
            merged,
            tuple(start for start, _end in merged),
            tuple(end for _start, end in merged),
        )
    return output


def configure_blacklist_args(args) -> None:
    """Validate/load --BlacklistRegion and attach its content signature."""
    path = getattr(args, "blacklist_region", None)
    if not path:
        args.blacklist_region = None
        args.blacklist_intervals = {}
        args.blacklist_region_signature = None
        return
    real_path = os.path.realpath(path)
    if not os.path.isfile(real_path):
        raise FileNotFoundError(f"Blacklist BED not found: {path}")
    args.blacklist_region = real_path
    args.blacklist_intervals = read_blacklist_bed(real_path)
    args.blacklist_region_signature = file_sha256(real_path)


def blacklist_overlaps(
    blacklist: BlacklistIndex, contig: str, start: int, end: int,
    sample: Optional[str] = None,
) -> Tuple[Tuple[int, int], ...]:
    """Return intervals for either CONTIG or SAMPLE#CONTIG BED naming."""
    if end <= start:
        return ()
    names = [contig]
    if sample and not contig.startswith(sample + "#"):
        names.append(f"{sample}#{contig}")
    overlaps: List[Tuple[int, int]] = []
    for name in names:
        indexed = blacklist.get(name)
        if indexed is None:
            continue
        first = bisect.bisect_right(indexed.ends, start)
        last = bisect.bisect_left(indexed.starts, end, lo=first)
        overlaps.extend(indexed.intervals[first:last])
    return tuple(merge_intervals(overlaps))


def _keep_unmapped_samples(args) -> Tuple[str, ...]:
    """Return a stable, de-duplicated sample exemption policy."""
    raw = getattr(args, "keep_unmapped_samples", DEFAULT_KEEP_UNMAPPED_SAMPLES)
    values = raw.split(",") if isinstance(raw, str) else raw
    return tuple(dict.fromkeys(
        str(sample).strip() for sample in values if str(sample).strip()
    ))


def _restart_signature(args) -> str:
    """Fingerprint options that can change biological results."""
    names = (
        "window_size", "min_segment", "min_identity", "promote_size", "flank",
        "local_extension", "novel_locus_extension", "cohort_window_size", "unmapped_main_size",
        "min_unmasked_fraction", "max_cycles", "cycle_minimap2",
        "blast_word_size", "blast_evalue", "blast_max_target_seqs",
        "no_repkmers", "skip_blastn",
    )
    reference_signature = getattr(
        args, "main_discovery_reference_signature", None,
    )
    if reference_signature is None:
        assembly_list = getattr(args, "assembly_list", None)
        if assembly_list:
            references = read_assembly_list(assembly_list)
            reference_signature = main_discovery_reference_signature(
                references[0],
            )
    payload = {name: getattr(args, name, None) for name in names}
    payload.update({
        "checkpoint_version": CHECKPOINT_VERSION,
        "main_discovery_target": "complete_first_assembly_v2",
        "main_discovery_reference_signature": reference_signature,
        "minimap_params": MINIMAP_PARAMS,
        "self_minimap_params": SELF_MINIMAP_PARAMS,
        "winnow_params": WINNOW_PARAMS,
    })
    if getattr(args, "full_main_discovery_cleaning", False):
        # Preserve legacy/default checkpoint hashes; only the opt-in graph
        # discovery path needs a distinct biological-result signature.
        payload["full_main_discovery_cleaning"] = True
    blacklist_signature = getattr(args, "blacklist_region_signature", None)
    if blacklist_signature is not None:
        # Omitting the key when no BED is supplied preserves checkpoint hashes
        # from runs predating blacklist support.
        payload["blacklist_region_signature"] = blacklist_signature
    encoded = json.dumps(payload, sort_keys=True, separators=(",", ":")).encode()
    return hashlib.sha256(encoded).hexdigest()


def _legacy_restart_signature_v1(args) -> str:
    """Signature used before complete-reference main discovery (read-only)."""
    names = (
        "window_size", "min_segment", "min_identity", "promote_size", "flank",
        "local_extension", "novel_locus_extension", "cohort_window_size",
        "unmapped_main_size", "min_unmasked_fraction", "max_cycles",
        "cycle_minimap2", "blast_word_size", "blast_evalue",
        "blast_max_target_seqs", "no_repkmers", "skip_blastn",
    )
    payload = {name: getattr(args, name, None) for name in names}
    payload.update({
        "checkpoint_version": CHECKPOINT_VERSION,
        "minimap_params": MINIMAP_PARAMS,
        "self_minimap_params": SELF_MINIMAP_PARAMS,
        "winnow_params": WINNOW_PARAMS,
    })
    return hashlib.sha256(
        json.dumps(payload, sort_keys=True, separators=(",", ":")).encode()
    ).hexdigest()


def main_discovery_reference_signature(assembly: AssemblyEntry) -> str:
    """Fingerprint the complete first assembly used as discovery template."""
    stat = os.stat(assembly.path)
    payload = {
        "sample": assembly.sample,
        "path": os.path.realpath(assembly.path),
        "size": stat.st_size,
        "mtime_ns": stat.st_mtime_ns,
    }
    return hashlib.sha256(
        json.dumps(payload, sort_keys=True, separators=(",", ":")).encode()
    ).hexdigest()


def _downstream_signature(args) -> str:
    """Fingerprint discovery settings plus the representative-priority policy."""
    payload = {
        "restart_signature": _restart_signature(args),
        "prior": getattr(args, "prior", DEFAULT_PRIOR),
        "priority_policy": "prior_bonus_10G_plus_length_v1",
    }
    encoded = json.dumps(payload, sort_keys=True, separators=(",", ":")).encode()
    return hashlib.sha256(encoded).hexdigest()


def _local_jobs(args) -> int:
    return getattr(args, "local_group_jobs", None) or getattr(args, "jobs", 1)


def _local_threads(args) -> int:
    return getattr(args, "local_group_threads", None) or args.threads


def _global_threads(args) -> int:
    """Configured total CPU budget: jobs times threads per job."""
    return max(1, int(getattr(args, "jobs", 1)) * int(args.threads))


def _parallel_task_threads(args, task_count: int) -> int:
    """Redistribute the CPU budget when fewer independent tasks than jobs exist."""
    override = getattr(args, "local_group_threads", None)
    if override is not None:
        return int(override)
    workers = max(1, min(_local_jobs(args), max(1, int(task_count))))
    return max(1, _global_threads(args) // workers)


def _merge_localize_results(results: Sequence[LocalizeResult]) -> LocalizeResult:
    """Combine placement evidence without duplicating an identical sequence."""
    placements: List[Placement] = []
    seen = set()
    notes: List[str] = []
    for result in results:
        if result.note and result.note not in notes:
            notes.append(result.note)
        for placement in result.placements:
            key = (placement.main, placement.start, placement.end, placement.strand)
            if key in seen:
                continue
            seen.add(key)
            placements.append(placement)
    if not placements:
        return LocalizeResult("unmapped", [], note=";".join(notes) or "no flank anchor")
    if len(placements) > 1:
        return LocalizeResult("both", placements, note=";".join(notes) or "multiple valid lifts")
    status = "mapped"
    if any(result.status == "one_sided" for result in results):
        status = "one_sided"
    return LocalizeResult(status, placements, note=";".join(notes))


def consolidate_identical_records(
    records: Sequence[Tuple[Piece, LocalizeResult]],
) -> List[Tuple[Piece, LocalizeResult]]:
    """Merge placement clones of the same exact source sequence.

    This is intentionally narrower than biological redundancy cleaning: records
    from different source occurrences remain independent even when their sequence
    is identical. Only clones sharing origin/source coordinates and sequence are
    folded into one sequence record carrying every valid placement.
    """
    grouped: Dict[Tuple[str, str, str, int, int, str], List[Tuple[Piece, LocalizeResult]]] = {}
    order: List[Tuple[str, str, str, int, int, str]] = []
    for piece, result in records:
        key = (piece.origin_id, piece.sample, piece.contig, piece.start, piece.end, piece.seq)
        if key not in grouped:
            grouped[key] = []
            order.append(key)
        grouped[key].append((piece, result))
    output: List[Tuple[Piece, LocalizeResult]] = []
    for key in order:
        members = grouped[key]
        representative = dataclasses.replace(members[0][0], final_id="")
        output.append((representative, _merge_localize_results([r for _p, r in members])))
    return output


def _relevant_lift_blocks(
    pieces: Sequence[Piece], hits: Sequence[object], flank: int,
) -> Dict[str, List[AnchorBlock]]:
    """Keep accepted alignment blocks near candidate source intervals only."""
    regions: Dict[str, List[Tuple[int, int]]] = {}
    for piece in pieces:
        regions.setdefault(piece.contig, []).append((
            max(0, piece.start - flank), piece.end + flank,
        ))
    region_index: Dict[str, Tuple[List[int], List[int]]] = {}
    for contig, vals in regions.items():
        merged = merge_intervals(vals)
        region_index[contig] = (
            [start for start, _end in merged],
            [end for _start, end in merged],
        )
    output: Dict[str, List[AnchorBlock]] = {}
    for hit in hits:
        query_id = getattr(hit, "query_id", "")
        indexed_regions = region_index.get(query_id)
        if not indexed_regions:
            continue
        starts, ends = indexed_regions
        score = float(getattr(hit, "identity", 0.0)) * float(
            getattr(hit, "aligned_bases", 0),
        ) / 100.0
        for q0, q1, t0, t1 in getattr(hit, "aligned_pairs", []):
            # Merged intervals are ordered and disjoint. Locate the first region
            # ending after q0 in O(log n), rather than scanning every residual
            # interval for every CIGAR block (which became quadratic on genomes
            # containing many thousands of discovered gaps).
            index = bisect.bisect_right(ends, q0)
            if index >= len(starts) or starts[index] >= q1:
                continue
            output.setdefault(query_id, []).append(AnchorBlock(
                int(q0), int(q1), int(t0), int(t1),
                str(getattr(hit, "target_id", "")), str(getattr(hit, "strand", "+")),
                float(getattr(hit, "identity", 0.0)), score,
            ))
    return output


def write_lift_blocks(
    path: str, pieces: Sequence[Piece], hits: Sequence[object], flank: int,
) -> Dict[str, List[AnchorBlock]]:
    blocks = _relevant_lift_blocks(pieces, hits, flank)
    Path(path).parent.mkdir(parents=True, exist_ok=True)
    with open(path, "w") as out:
        out.write("query_contig\tq0\tq1\ttarget\tt0\tt1\tstrand\tidentity\talignment_score\n")
        for contig in sorted(blocks):
            for block in blocks[contig]:
                out.write("\t".join([
                    contig, str(block.q0), str(block.q1), block.main,
                    str(block.t0), str(block.t1), block.strand,
                    repr(block.identity), repr(block.alignment_score),
                ]) + "\n")
    return blocks


def read_lift_blocks(
    path: str, contigs: Optional[Sequence[str]] = None,
) -> Dict[str, List[AnchorBlock]]:
    """Read compact lift evidence, optionally restricted while streaming."""
    if not os.path.isfile(path):
        raise RuntimeError(f"Completed stage is missing {LIFT_BLOCKS_NAME}: {path}")
    wanted = None if contigs is None else set(contigs)
    output: Dict[str, List[AnchorBlock]] = {}
    with open(path) as handle:
        next(handle, None)
        for line_no, raw in enumerate(handle, 2):
            if not raw.strip():
                continue
            cols = raw.rstrip("\n").split("\t")
            if len(cols) != 9:
                raise RuntimeError(f"Invalid lift block {path}:{line_no}")
            if wanted is not None and cols[0] not in wanted:
                continue
            try:
                block = AnchorBlock(
                    int(cols[1]), int(cols[2]), int(cols[4]), int(cols[5]),
                    cols[3], cols[6], float(cols[7]), float(cols[8]),
                )
            except ValueError as exc:
                raise RuntimeError(f"Invalid lift block {path}:{line_no}: {exc}") from exc
            output.setdefault(cols[0], []).append(block)
    return output


def localize_from_block_map(
    pieces: Sequence[Piece], blocks_by_contig: Dict[str, List[AnchorBlock]], flank: int,
) -> List[LocalizeResult]:
    """Lift pieces in input order with a streaming coordinate sweep."""
    return localize_windows_from_blocks(
        [(piece.contig, piece.start, piece.end) for piece in pieces],
        blocks_by_contig, flank,
    )


def prune_old_iterations(outdir: str, current_iteration: str, args,
                         log: logging.Logger) -> None:
    """Keep the latest successful iteration and remove older iteration work.

    Pruning happens only after ``current_iteration`` has completed successfully.
    A failed active iteration is therefore retained alongside the previous
    successful one. ``--keep-work`` disables all pruning.
    """
    if getattr(args, "keep_work", False):
        return
    root = os.path.join(outdir, "iterations")
    if not os.path.isdir(root):
        return
    current = os.path.realpath(current_iteration)
    for name in os.listdir(root):
        path = os.path.join(root, name)
        if not os.path.isdir(path) or os.path.realpath(path) == current:
            continue
        try:
            shutil.rmtree(path)
            log.debug("Pruned older iteration work: %s", path)
        except OSError as exc:
            log.warning("Could not prune older iteration %s: %s", path, exc)


# ---------------------------------------------------------------------------
# Assembly list + windowing
# ---------------------------------------------------------------------------

class AssemblyEntry:
    def __init__(self, sample: str, path: str):
        self.sample = sample
        self.path = path


def read_assembly_list(path: str) -> List[AssemblyEntry]:
    out: List[AssemblyEntry] = []
    seen_samples: Dict[str, int] = {}
    seen_folder_names: Dict[str, Tuple[str, int]] = {}
    with open(path) as fh:
        for line_no, raw in enumerate(fh, 1):
            line = raw.strip()
            if not line or line.startswith("#"):
                continue
            cols = line.split()
            if len(cols) != 2:
                raise ValueError(f"assembly list line {line_no}: expected exactly NAME FASTA; index must be FASTA.fai")
            cols[1] = os.path.abspath(os.path.expanduser(cols[1]))
            if not os.path.exists(cols[1]):
                raise FileNotFoundError(f"assembly list line {line_no}: {cols[1]} not found")
            sample = cols[0]
            if sample in seen_samples:
                raise ValueError(
                    f"assembly list line {line_no}: duplicate haplotype name {sample!r} "
                    f"(first seen on line {seen_samples[sample]})"
                )
            folder_name = sanitize_id(sample)
            if folder_name in seen_folder_names:
                previous_sample, previous_line = seen_folder_names[folder_name]
                raise ValueError(
                    f"assembly list line {line_no}: haplotype names {previous_sample!r} "
                    f"(line {previous_line}) and {sample!r} map to the same checkpoint "
                    f"folder {folder_name!r}"
                )
            seen_samples[sample] = line_no
            seen_folder_names[folder_name] = (sample, line_no)
            out.append(AssemblyEntry(sample, cols[1]))
    if not out:
        raise ValueError("assembly list is empty")
    return out


def build_sample_priorities(
    assemblies: Sequence[AssemblyEntry], prior_text: str,
) -> Tuple[Dict[str, int], List[str]]:
    """Expand ordered prior tiers into per-sample score bonuses."""
    tokens = [token.strip() for token in str(prior_text).split(",") if token.strip()]
    priorities: Dict[str, int] = {}
    unmatched: List[str] = []
    samples = [assembly.sample for assembly in assemblies]
    matched_groups: List[Tuple[str, List[str]]] = []
    for token in tokens:
        if "_" in token:
            matches = [sample for sample in samples if sample == token]
        else:
            matches = [
                sample for sample in samples
                if sample == token or sample.startswith(token + "_")
            ]
        if not matches:
            unmatched.append(token)
            continue
        matched_groups.append((token, matches))
    total_included = len(matched_groups)
    for index, (_token, matches) in enumerate(matched_groups):
        bonus = (total_included - index) * PRIOR_SCORE_STEP
        for sample in matches:
            # The earliest matching prior entry owns the stronger tier.
            priorities.setdefault(sample, bonus)
    return priorities, unmatched


def window_contig(reader: IndexedFasta, sample: str, source_path: str, contig: str,
                  window_size: int) -> List[Piece]:
    """Window one contig into ACGT/acgt pieces (split at N, all kept)."""
    pieces: List[Piece] = []
    safe_sample = sanitize_id(sample)
    safe_contig = sanitize_id(contig)
    length = reader.length(contig)
    seq = reader.sequence(contig)
    counter = 0
    for start in range(0, length, window_size):
        end = min(length, start + window_size)
        for la, lb, sub in iter_n_free_segments(seq[start:end]):
            counter += 1
            s0, s1 = start + la, start + lb
            pieces.append(Piece(f"win_{safe_sample}_{safe_contig}_{s0}_{s1}_{counter:08d}",
                                sample, contig, s0, s1, "+", sub, note="window"))
    return pieces


_MP_READER: Optional[IndexedFasta] = None


def _mp_init(reader: IndexedFasta) -> None:
    global _MP_READER
    _MP_READER = reader


def _mp_window_task(task):
    sample, source_path, contig, window_size = task
    return window_contig(_MP_READER, sample, source_path, contig, window_size)


def remove_blacklisted_piece_regions(
    pieces: Sequence[Piece], blacklist: BlacklistIndex,
) -> Tuple[List[Piece], int]:
    """Subtract source-coordinate BED intervals without scanning unrelated loci."""
    if not blacklist:
        return list(pieces), 0
    output: List[Piece] = []
    removed_bp = 0
    for piece in pieces:
        masks = blacklist_overlaps(
            blacklist, piece.contig, piece.start, piece.end, piece.sample,
        )
        if not masks:
            output.append(piece)
            continue
        retained = subtract_intervals([(piece.start, piece.end)], masks)
        removed_bp += piece.length - sum(end - start for start, end in retained)
        for index, (start, end) in enumerate(retained, 1):
            offset_start = start - piece.start
            offset_end = end - piece.start
            derived_id = f"{piece.temp_id}_blacklist{index:04d}"
            output.append(dataclasses.replace(
                piece,
                temp_id=derived_id,
                start=start,
                end=end,
                seq=piece.seq[offset_start:offset_end],
                origin_id=derived_id,
                lift_id=derived_id,
                note=(piece.note + ";blacklist_subtracted").lstrip(";"),
            ))
    return output, removed_bp


def window_assembly(assembly: AssemblyEntry, window_size: int, processes: int,
                    log: logging.Logger,
                    blacklist: Optional[BlacklistIndex] = None) -> List[Piece]:
    pieces: List[Piece] = []
    with IndexedFasta(assembly.path) as fa:
        contigs = fa.names()
        if processes and processes > 1 and len(contigs) > 1:
            tasks = [(assembly.sample, assembly.path, c, window_size) for c in contigs]
            with mp_context().Pool(processes=min(processes, len(contigs)),
                                   initializer=_mp_init, initargs=(fa,)) as pool:
                for sub in pool.map(_mp_window_task, tasks):
                    pieces.extend(sub)
        else:
            for c in contigs:
                pieces.extend(window_contig(fa, assembly.sample, assembly.path, c, window_size))
    before = len(pieces)
    pieces, removed_bp = remove_blacklisted_piece_regions(pieces, blacklist or {})
    if removed_bp:
        log.info(
            "%s: blacklist removed %d bp from first-reference windows (%d -> %d pieces)",
            assembly.sample, removed_bp, before, len(pieces),
        )
    log.info("Windowed %s into %d pieces", assembly.sample, len(pieces))
    return pieces


# ---------------------------------------------------------------------------
# Aligner closures
# ---------------------------------------------------------------------------

def make_self_align_fn(workdir: str, args, log: logging.Logger,
                       blast_mode: str = "independent",
                       broad_aligner: str = "winnowmap",
                       run_blastn: bool = True) -> Callable:
    """Return a self-align callable for self_clean_converge (winnowmap+blastn)."""
    minimap_params = (
        SELF_MINIMAP_PARAMS
        if broad_aligner == "minimap2"
        else MINIMAP_PARAMS
    )

    def align_fn(pieces: Sequence[Piece], cycle: int):
        cdir = os.path.join(workdir, f"cycle{cycle:02d}")
        Path(cdir).mkdir(parents=True, exist_ok=True)
        cand_fa = os.path.join(cdir, "cand.fa")
        write_fasta_pieces(pieces, cand_fa)
        write_piece_map(pieces, os.path.join(cdir, "cand.map.tsv"))
        cycle_rep_kmers = (
            getattr(args, "winnow_rep_kmers", None)
            if broad_aligner == "winnowmap" and not args.no_repkmers else None
        )
        return collect_alignment_hits(
            cand_fa, cand_fa, cdir, args.threads, args.min_identity, args.min_segment, log,
            winnow_params=WINNOW_PARAMS, rep_kmers=cycle_rep_kmers,
            blast_word_size=args.blast_word_size, blast_evalue=args.blast_evalue,
            blast_max_target_seqs=args.blast_max_target_seqs, prefix="self",
            skip_blastn=args.skip_blastn or not run_blastn,
            blast_mode=blast_mode,
            broad_aligner=broad_aligner, minimap_params=minimap_params,
        )

    # A completed cycle is reusable only under the same alignment policy.
    # In particular, changing Minimap2's secondary-chain policy must force
    # --resume to regenerate its PAF and kept sequence.
    align_fn.checkpoint_tag = json.dumps({  # type: ignore[attr-defined]
        "broad_aligner": broad_aligner,
        "blast_mode": blast_mode,
        "run_blastn": bool(run_blastn),
        "min_identity": args.min_identity,
        "min_segment": args.min_segment,
        "minimap_params": minimap_params,
        "winnow_params": WINNOW_PARAMS,
        "blast_word_size": args.blast_word_size,
        "blast_evalue": args.blast_evalue,
        "blast_max_target_seqs": args.blast_max_target_seqs,
        "skip_blastn": bool(args.skip_blastn),
    }, sort_keys=True, separators=(",", ":"))
    return align_fn


def self_clean_minimap_then_winnow(
    pieces: Sequence[Piece],
    workdir: str,
    args,
    mode: str,
    log: logging.Logger,
    label: str,
    unlimited: bool = False,
) -> List[Piece]:
    """Coarse minimap2 self-clean, then independent Winnowmap plus BLASTN."""
    if not pieces:
        return []
    global_args = argparse.Namespace(**vars(args))
    global_args.threads = _global_threads(args)
    log.info("%s: using %d threads for each cohort-wide aligner",
             label, global_args.threads)
    minimap_cycles = None if unlimited else 1 + args.cycle_minimap2
    winnow_cycles = None if unlimited else args.max_cycles
    if unlimited:
        log.info("%s: minimap2 and Winnowmap+BLASTN run until convergence (no cycle cap)",
                 label)
    else:
        log.info("%s: cycle limits are minimap2=%d and Winnowmap+BLASTN=%d",
                 label, minimap_cycles, winnow_cycles)
    minimap_dir = os.path.join(workdir, "minimap2")
    minimap_fn = make_self_align_fn(
        minimap_dir, global_args, log,
        blast_mode="residual", broad_aligner="minimap2", run_blastn=False,
    )
    coarse = self_clean_converge(
        pieces, minimap_fn, args.min_segment, mode, log,
        max_cycles=minimap_cycles,
        checkpoint_dir=minimap_dir,
        resume=bool(getattr(args, "resume", False)),
        sample_priorities=getattr(args, "sample_priorities", None),
    )
    log.info("%s: minimap2 self-clean retained %d/%d pieces",
             label, len(coarse), len(pieces))
    if not coarse:
        return []
    winnow_dir = os.path.join(workdir, "winnowmap_blastn")
    winnow_fn = make_self_align_fn(
        winnow_dir, global_args, log,
        blast_mode="independent", broad_aligner="winnowmap", run_blastn=True,
    )
    cleaned = self_clean_converge(
        coarse, winnow_fn, args.min_segment, mode, log,
        max_cycles=winnow_cycles,
        checkpoint_dir=winnow_dir,
        resume=bool(getattr(args, "resume", False)),
        sample_priorities=getattr(args, "sample_priorities", None),
    )
    log.info("%s: Winnowmap+BLASTN self-clean retained %d/%d minimap2 residual pieces",
             label, len(cleaned), len(coarse))
    return cleaned


def assign_final_ids(pieces: Sequence[Piece], prefix: str, counters: Dict[str, int]) -> None:
    for p in pieces:
        key = sanitize_id(prefix)
        counters[key] = counters.get(key, 0) + 1
        p.final_id = f"{key}_{counters[key]:06d}"


def write_piece_map(pieces: Sequence[Piece], path: str) -> None:
    """Persist the source-coordinate meaning of temporary PAF query names."""
    Path(path).parent.mkdir(parents=True, exist_ok=True)
    with open(path, "w") as out:
        out.write("query_id\tsample\tcontig\tstart\tend\tstrand\torigin_id\tlift_id\n")
        for piece in pieces:
            out.write("\t".join([
                piece.temp_id, piece.sample, piece.contig, str(piece.start), str(piece.end),
                piece.strand, piece.origin_id, piece.lift_id,
            ]) + "\n")


def merge_source_adjacent_pieces(pieces: Sequence[Piece]) -> List[Piece]:
    """Merge retained first-reference pieces that overlap or touch in source.

    Window boundaries are alignment implementation details, not final reference
    boundaries. Genuine removed source gaps remain separate records.
    """
    group_order: Dict[Tuple[str, str, str], int] = {}
    grouped: Dict[Tuple[str, str, str], List[Piece]] = {}
    for piece in pieces:
        key = (piece.sample, piece.contig, piece.strand)
        if key not in group_order:
            group_order[key] = len(group_order)
            grouped[key] = []
        grouped[key].append(piece)

    output: List[Piece] = []
    for key in sorted(grouped, key=lambda k: group_order[k]):
        members = sorted(grouped[key], key=lambda p: (p.start, p.end, p.temp_id))
        current: Optional[Piece] = None
        for piece in members:
            if current is None:
                current = dataclasses.replace(piece, final_id="")
                continue
            if piece.start > current.end:
                output.append(current)
                current = dataclasses.replace(piece, final_id="")
                continue
            overlap = current.end - piece.start
            if overlap >= piece.length:
                continue
            if overlap > 0 and current.seq[-overlap:] != piece.seq[:overlap]:
                raise ValueError(
                    f"Conflicting overlapping retained pieces at "
                    f"{piece.sample}:{piece.contig}:{piece.start}-{current.end}"
                )
            current = Piece(
                temp_id=current.temp_id,
                sample=current.sample,
                contig=current.contig,
                start=current.start,
                end=piece.end,
                strand=current.strand,
                seq=current.seq + piece.seq[overlap:],
                parent_id=current.parent_id,
                note="merged_source_adjacent",
                origin_id=current.origin_id,
            )
        if current is not None:
            output.append(current)
    return output


def read_cleaned_reference_intervals(
    cleaned_fa: str,
    expected_sample: Optional[str] = None,
) -> Tuple[str, Dict[str, List[Tuple[int, int]]]]:
    """Read original source intervals from a minsetref cleaned-main FASTA."""
    samples = set()
    intervals: Dict[str, List[Tuple[int, int]]] = {}
    with open(cleaned_fa) as handle:
        for line_number, raw in enumerate(handle, 1):
            if not raw.startswith(">"):
                continue
            fields = raw[1:].strip().split()
            if len(fields) < 2:
                raise ValueError(
                    f"Cleaned-reference header lacks source coordinates at "
                    f"{cleaned_fa}:{line_number}"
                )
            match = SOURCE_COORD_RE.match(fields[-1])
            if match is None:
                raise ValueError(
                    f"Cannot parse source coordinate {fields[-1]!r} at "
                    f"{cleaned_fa}:{line_number}"
                )
            sample, contig, start_text, end_text, strand = match.groups()
            start, end = int(start_text), int(end_text)
            if strand != "+":
                raise ValueError(
                    f"First-reference interval must use + strand at "
                    f"{cleaned_fa}:{line_number}"
                )
            if end <= start:
                raise ValueError(
                    f"Invalid source interval {start}-{end} at "
                    f"{cleaned_fa}:{line_number}"
                )
            if expected_sample is not None and sample != expected_sample:
                raise ValueError(
                    f"Cleaned-reference sample {sample!r} does not match "
                    f"expected sample {expected_sample!r}"
                )
            samples.add(sample)
            intervals.setdefault(contig, []).append((start, end))
    if not intervals:
        raise ValueError(f"No source-coordinate records found in {cleaned_fa}")
    if len(samples) != 1:
        raise ValueError(
            f"Expected one first-reference sample in {cleaned_fa}; found "
            f"{sorted(samples)!r}"
        )
    sample = expected_sample if expected_sample is not None else next(iter(samples))
    return sample, {
        contig: merge_intervals(contig_intervals)
        for contig, contig_intervals in intervals.items()
    }


def write_filtered_reference_fasta(
    source_fa: str,
    sample: str,
    retained_by_contig: Dict[str, Sequence[Tuple[int, int]]],
    output_fa: str,
    log: Optional[logging.Logger] = None,
) -> Tuple[int, int, int]:
    """Write N-free source sequence absent from the final cleaned reference.

    Returns ``(record_count, written_bp, omitted_ambiguous_bp)``. Coordinates
    are half-open source-assembly coordinates. Non-ACGT runs are omitted to
    mirror ``window_assembly``, where they were never alignment candidates.
    """
    output_real = os.path.realpath(output_fa)
    if output_real == os.path.realpath(source_fa):
        raise ValueError("Filtered-sequence output must differ from the source FASTA")
    Path(output_fa).parent.mkdir(parents=True, exist_ok=True)
    temporary = output_fa + ".tmp"
    record_count = 0
    written_bp = 0
    omitted_ambiguous_bp = 0
    safe_sample = sanitize_id(sample)
    try:
        with IndexedFasta(source_fa) as source, open(temporary, "w") as out:
            source_contigs = source.names()
            missing_contigs = sorted(set(retained_by_contig) - set(source_contigs))
            if missing_contigs:
                raise ValueError(
                    f"Retained cleaned-reference contigs are absent from {source_fa}: "
                    + ", ".join(missing_contigs[:10])
                )
            for contig in source_contigs:
                contig_length = source.length(contig)
                retained = []
                for start, end in merge_intervals(retained_by_contig.get(contig, ())):
                    if start < 0 or end > contig_length:
                        raise ValueError(
                            f"Retained interval {contig}:{start}-{end} lies outside "
                            f"source length {contig_length}"
                        )
                    retained.append((start, end))
                cursor = 0
                excluded_intervals: List[Tuple[int, int]] = []
                for start, end in retained:
                    if start > cursor:
                        excluded_intervals.append((cursor, start))
                    cursor = max(cursor, end)
                if cursor < contig_length:
                    excluded_intervals.append((cursor, contig_length))

                safe_contig = sanitize_id(contig)
                for start, end in excluded_intervals:
                    sequence = source.fetch(contig, start, end)
                    nfree_bp = 0
                    for local_start, local_end, subsequence in iter_n_free_segments(sequence):
                        source_start = start + local_start
                        source_end = start + local_end
                        record_count += 1
                        nfree_bp += len(subsequence)
                        written_bp += len(subsequence)
                        record_id = (
                            f"filtered_{safe_sample}_{safe_contig}_{record_count:08d}"
                        )
                        out.write(
                            f">{record_id} {sample}:{contig}:"
                            f"{source_start}-{source_end}+\n{wrap_fasta(subsequence)}\n"
                        )
                    omitted_ambiguous_bp += len(sequence) - nfree_bp
            out.flush()
            os.fsync(out.fileno())
        os.replace(temporary, output_fa)
        try:
            os.remove(output_fa + ".fai")
        except FileNotFoundError:
            pass
    finally:
        try:
            os.remove(temporary)
        except FileNotFoundError:
            pass
    if log is not None:
        log.info(
            "Filtered first-reference sequence: %d records, %d bp -> %s "
            "(%d ambiguous/N bp omitted)",
            record_count, written_bp, output_fa, omitted_ambiguous_bp,
        )
    return record_count, written_bp, omitted_ambiguous_bp


def _bidirectional_score_breaks(
    unique_intervals: Sequence[Tuple[int, int]],
) -> List[int]:
    """Break coordinates from forward/reverse +unique/-nonunique scoring."""
    if not unique_intervals:
        return []
    segments: List[Tuple[int, int, bool]] = []
    cursor = unique_intervals[0][0]
    for start, end in unique_intervals:
        if start > cursor:
            segments.append((cursor, start, False))
        segments.append((start, end, True))
        cursor = end

    breaks = set()
    score = 0
    for start, end, is_unique in segments:
        score += (end - start) if is_unique else -(end - start)
        if score < 0:
            breaks.add(end)
            score = 0

    score = 0
    for start, end, is_unique in reversed(segments):
        score += (end - start) if is_unique else -(end - start)
        if score < 0:
            breaks.add(start)
            score = 0
    return sorted(breaks)


def reassemble_unique_regions(
    pieces: Sequence[Piece],
    assembly_paths: Dict[str, str],
    max_gap: int = 5_000,
    log: Optional[logging.Logger] = None,
) -> List[Piece]:
    """Reconstruct scored same-contig spans from final unique source intervals.

    Gaps over ``max_gap`` and gaps containing N/n are unconditional breaks.
    Within each connected region, a nonunique gap is retained only when the
    bidirectional score rule leaves it inside an unbroken output span.
    """
    if not pieces:
        return []
    grouped: Dict[Tuple[str, str], List[Piece]] = {}
    group_order: List[Tuple[str, str]] = []
    for piece in pieces:
        key = (piece.sample, piece.contig)
        if key not in grouped:
            grouped[key] = []
            group_order.append(key)
        grouped[key].append(piece)

    output: List[Piece] = []
    output_counter = 0
    for sample, contig in group_order:
        fasta_path = assembly_paths.get(sample)
        if fasta_path is None:
            raise RuntimeError(f"No source assembly FASTA registered for {sample!r}")
        members = sorted(
            grouped[(sample, contig)], key=lambda p: (p.start, p.end, p.temp_id),
        )
        unique_intervals = merge_intervals([
            (piece.start, piece.end) for piece in members if piece.end > piece.start
        ])
        if not unique_intervals:
            continue
        with IndexedFasta(fasta_path) as source:
            if contig not in set(source.names()):
                raise RuntimeError(f"Source contig {contig!r} missing from {fasta_path}")
            clusters: List[List[Tuple[int, int]]] = []
            current = [unique_intervals[0]]
            for interval in unique_intervals[1:]:
                gap_start, gap_end = current[-1][1], interval[0]
                gap_size = max(0, gap_end - gap_start)
                hard_break = gap_size > max_gap
                if gap_size and not hard_break:
                    gap_sequence = source.fetch(contig, gap_start, gap_end)
                    hard_break = "N" in gap_sequence.upper()
                if hard_break:
                    clusters.append(current)
                    current = [interval]
                else:
                    current.append(interval)
            clusters.append(current)

            for cluster in clusters:
                cluster_start, cluster_end = cluster[0][0], cluster[-1][1]
                boundaries = [cluster_start]
                boundaries.extend(
                    point for point in _bidirectional_score_breaks(cluster)
                    if cluster_start < point < cluster_end
                )
                boundaries.append(cluster_end)
                boundaries = sorted(set(boundaries))
                for start, end in zip(boundaries, boundaries[1:]):
                    if end <= start:
                        continue
                    overlapping_unique = [
                        (a, b) for a, b in cluster if a < end and b > start
                    ]
                    if not overlapping_unique:
                        continue
                    source_members = [
                        piece for piece in members
                        if piece.start < end and piece.end > start
                    ]
                    if not source_members:
                        continue
                    representative = source_members[0]
                    output_counter += 1
                    unique_bp = sum(
                        min(end, b) - max(start, a) for a, b in overlapping_unique
                    )
                    output.append(Piece(
                        temp_id=(
                            f"reassembled_{sanitize_id(sample)}_{sanitize_id(contig)}_"
                            f"{start}_{end}_{output_counter:08d}"
                        ),
                        sample=sample,
                        contig=contig,
                        start=start,
                        end=end,
                        strand=representative.strand,
                        seq=source.fetch(contig, start, end),
                        parent_id=representative.temp_id,
                        note=f"score_reassembled;unique_bp={unique_bp}",
                        origin_id=representative.origin_id,
                        lift_id=representative.lift_id,
                        within_sample_rank=representative.within_sample_rank,
                    ))
    if log is not None:
        log.info(
            "Score reassembly: %d unique pieces -> %d reconstructed intervals "
            "(hard gap > %d bp)",
            len(pieces), len(output), max_gap,
        )
    return output


def merge_weighted_source_gaps(
    pieces: Sequence[Piece],
    assembly_paths: Dict[str, str],
    max_weighted_gap: float = NOVEL_MERGE_WEIGHT_THRESHOLD,
    masked_weight: float = NOVEL_MERGE_MASKED_WEIGHT,
    log: Optional[logging.Logger] = None,
) -> List[Piece]:
    """Merge residual source intervals before the final size filter.

    Adjacent intervals from the same sample and source contig are joined when
    the intervening source sequence has strict weighted distance

        unmasked_bp + masked_weight * masked_bp < max_weighted_gap.

    Lowercase A/C/G/T are masked and uppercase A/C/G/T are unmasked; other
    symbols contribute zero, matching the pipeline's existing novelty-score
    convention. The source FASTA supplies the entire merged sequence, so
    soft-masked gap sequence is retained rather than concatenating only the
    residual endpoints.
    """
    if not pieces:
        return []
    if max_weighted_gap <= 0:
        raise ValueError("max_weighted_gap must be positive")
    if masked_weight < 0:
        raise ValueError("masked_weight must be non-negative")

    grouped: Dict[Tuple[str, str], List[Piece]] = {}
    group_order: List[Tuple[str, str]] = []
    for piece in pieces:
        key = (piece.sample, piece.contig)
        if key not in grouped:
            grouped[key] = []
            group_order.append(key)
        grouped[key].append(piece)

    output: List[Piece] = []
    output_counter = 0
    merged_gap_count = 0
    merged_gap_bp = 0
    for sample, contig in group_order:
        fasta_path = assembly_paths.get(sample)
        if fasta_path is None:
            raise RuntimeError(f"No source assembly FASTA registered for {sample!r}")
        members = sorted(
            grouped[(sample, contig)], key=lambda p: (p.start, p.end, p.temp_id),
        )
        with IndexedFasta(fasta_path) as source:
            if contig not in set(source.names()):
                raise RuntimeError(f"Source contig {contig!r} missing from {fasta_path}")
            contig_length = source.length(contig)
            for piece in members:
                if piece.start < 0 or piece.end > contig_length or piece.end <= piece.start:
                    raise ValueError(
                        f"Invalid residual interval {sample}:{contig}:"
                        f"{piece.start}-{piece.end} for contig length {contig_length}"
                    )

            clusters: List[Tuple[int, int, Piece, int]] = []
            cluster_start = members[0].start
            cluster_end = members[0].end
            representative = members[0]
            source_piece_count = 1
            for piece in members[1:]:
                if piece.start <= cluster_end:
                    cluster_end = max(cluster_end, piece.end)
                    source_piece_count += 1
                    continue
                gap_start, gap_end = cluster_end, piece.start
                gap_sequence = source.fetch(contig, gap_start, gap_end)
                unmasked_bp = count_unmasked(gap_sequence)
                masked_bp = count_masked(gap_sequence)
                weighted_gap = unmasked_bp + masked_weight * masked_bp
                if weighted_gap < max_weighted_gap:
                    cluster_end = piece.end
                    source_piece_count += 1
                    merged_gap_count += 1
                    merged_gap_bp += gap_end - gap_start
                    continue
                clusters.append((
                    cluster_start, cluster_end, representative, source_piece_count,
                ))
                cluster_start = piece.start
                cluster_end = piece.end
                representative = piece
                source_piece_count = 1
            clusters.append((
                cluster_start, cluster_end, representative, source_piece_count,
            ))

            for start, end, representative, source_piece_count in clusters:
                output_counter += 1
                output.append(Piece(
                    temp_id=(
                        f"weighted_{sanitize_id(sample)}_{sanitize_id(contig)}_"
                        f"{start}_{end}_{output_counter:08d}"
                    ),
                    sample=sample,
                    contig=contig,
                    start=start,
                    end=end,
                    strand=representative.strand,
                    seq=source.fetch(contig, start, end),
                    parent_id=representative.temp_id,
                    note=(
                        "weighted_gap_remerged;"
                        f"source_pieces={source_piece_count};"
                        f"threshold={max_weighted_gap:g};"
                        f"masked_weight={masked_weight:g}"
                    ),
                    origin_id=representative.origin_id,
                    lift_id=representative.lift_id,
                    within_sample_rank=representative.within_sample_rank,
                ))
    if log is not None:
        log.info(
            "Weighted-gap remerge: %d residual pieces -> %d source intervals; "
            "bridged %d gaps (%d physical bp) with unmasked + %.3g*masked < %.3g",
            len(pieces), len(output), merged_gap_count, merged_gap_bp,
            masked_weight, max_weighted_gap,
        )
    return output


def write_anchored_reassembled_fasta(
    pieces: Sequence[Piece],
    assembly_paths: Dict[str, str],
    output_fa: str,
    anchor_bp: int = NOVEL_ANCHOR_BP,
) -> Dict[str, Tuple[int, int]]:
    """Write each reassembled core with up to ``anchor_bp`` source bases per side.

    Returns the core's half-open coordinates within each anchored FASTA record.
    The core itself is not modified; anchors exist only for the lift alignment.
    """
    Path(output_fa).parent.mkdir(parents=True, exist_ok=True)
    by_sample: Dict[str, List[Piece]] = {}
    sample_order: List[str] = []
    for piece in pieces:
        if piece.sample not in by_sample:
            by_sample[piece.sample] = []
            sample_order.append(piece.sample)
        by_sample[piece.sample].append(piece)

    core_coordinates: Dict[str, Tuple[int, int]] = {}
    with open(output_fa, "w") as out:
        for sample in sample_order:
            source_path = assembly_paths.get(sample)
            if source_path is None:
                raise RuntimeError(f"No source assembly FASTA registered for {sample!r}")
            with IndexedFasta(source_path) as source:
                available = set(source.names())
                for piece in by_sample[sample]:
                    if piece.contig not in available:
                        raise RuntimeError(
                            f"Source contig {piece.contig!r} missing from {source_path}"
                        )
                    anchor_start = max(0, piece.start - anchor_bp)
                    anchor_end = min(source.length(piece.contig), piece.end + anchor_bp)
                    core_start = piece.start - anchor_start
                    core_end = core_start + piece.length
                    sequence = source.fetch(piece.contig, anchor_start, anchor_end)
                    core_coordinates[piece.temp_id] = (core_start, core_end)
                    out.write(
                        f">{piece.temp_id} source={piece.sample}:{piece.contig}:"
                        f"{piece.start}-{piece.end} core={core_start}-{core_end} "
                        f"anchor_bp={anchor_bp}\n{wrap_fasta(sequence)}\n"
                    )
    return core_coordinates


def _lift_header(piece: Piece, result: LocalizeResult, use_final_id: bool) -> str:
    primary = piece.final_id if use_final_id and piece.final_id else piece.temp_id
    placements = ";".join(
        f"{placement.main}:{placement.start}-{placement.end}:{placement.strand}"
        for placement in result.placements
    ) or "."
    note = (result.note or ".").replace(" ", "_")
    return (
        f">{primary} source={piece.sample}:{piece.contig}:{piece.start}-{piece.end} "
        f"length={piece.length} lift_status={result.status} "
        f"reference={placements} lift_note={note}"
    )


def write_annotated_loci_fasta(
    records: Sequence[Tuple[Piece, LocalizeResult]],
    output_fa: str,
    use_final_id: bool = False,
) -> None:
    """Write core sequences with fresh main-lift annotations in their headers."""
    Path(output_fa).parent.mkdir(parents=True, exist_ok=True)
    with open(output_fa, "w") as out:
        for piece, result in records:
            out.write(f"{_lift_header(piece, result, use_final_id)}\n{wrap_fasta(piece.seq)}\n")


def lift_reassembled_loci(
    pieces: Sequence[Piece],
    assembly_paths: Dict[str, str],
    main_target: str,
    workdir: str,
    args,
    log: logging.Logger,
) -> List[Tuple[Piece, LocalizeResult]]:
    """Freshly localize weighted-gap-remerged cores using 10-kb source anchors."""
    if not pieces:
        return []
    Path(workdir).mkdir(parents=True, exist_ok=True)
    anchored_fa = os.path.join(
        os.path.dirname(workdir), "remerged.with_10kb_anchors.fa",
    )
    core_coordinates = write_anchored_reassembled_fasta(
        pieces, assembly_paths, anchored_fa, NOVEL_ANCHOR_BP,
    )
    lift_args = argparse.Namespace(**vars(args))
    lift_args.threads = _global_threads(args)
    log.info(
        "Reassembled-locus lift: aligning %d anchored records to main with %d threads",
        len(pieces), lift_args.threads,
    )
    hits = collect_alignment_hits(
        anchored_fa, main_target, workdir, lift_args.threads,
        args.min_identity, args.min_segment, log,
        winnow_params=WINNOW_PARAMS,
        blast_word_size=args.blast_word_size,
        blast_evalue=args.blast_evalue,
        blast_max_target_seqs=args.blast_max_target_seqs,
        prefix="remerged_to_main", skip_blastn=True,
        blast_mode="residual", broad_aligner="minimap2",
        minimap_params=MINIMAP_PARAMS,
    )
    hits_by_query: Dict[str, List[object]] = {}
    for hit in hits:
        hits_by_query.setdefault(hit.query_id, []).append(hit)

    records: List[Tuple[Piece, LocalizeResult]] = []
    for piece in pieces:
        core_start, core_end = core_coordinates[piece.temp_id]
        blocks = anchor_blocks_from_hits(
            hits_by_query.get(piece.temp_id, ()), piece.temp_id,
        )
        result = localize_insertion(
            core_start, core_end, blocks, flank=NOVEL_ANCHOR_BP,
        )
        records.append((piece, result))
    mapped = sum(bool(result.placements) for _piece, result in records)
    log.info(
        "Reassembled-locus lift complete: %d mapped/anchored, %d unmapped",
        mapped, len(records) - mapped,
    )
    return records


# ---------------------------------------------------------------------------
# Stage 1: build main chroms from the first genome
# ---------------------------------------------------------------------------

def build_main(first: AssemblyEntry, args, outdir: str, log: logging.Logger) -> Tuple[List[Piece], str, Optional[str]]:
    log.info("Building main chroms from first genome %s", first.sample)
    idir = os.path.join(outdir, "iterations", f"000_{sanitize_id(first.sample)}_main")
    Path(idir).mkdir(parents=True, exist_ok=True)
    candidates = window_assembly(
        first, args.window_size, args.processes, log,
        getattr(args, "blacklist_intervals", None),
    )
    global_args = argparse.Namespace(**vars(args))
    global_args.threads = _global_threads(args)
    log.info("First-reference self-clean: using %d cohort-wide aligner threads",
             global_args.threads)
    # Job 1, phase A: always run the initial minimap2 pass, followed by the
    # requested number of residual minimap2 cycles before Winnowmap begins.
    minimap_fn = make_self_align_fn(
        os.path.join(idir, "selfclean_minimap2"), global_args, log,
        blast_mode="residual", broad_aligner="minimap2", run_blastn=False,
    )
    coarse = self_clean_converge(
        candidates, minimap_fn, args.min_segment, "first_genome_minimap2", log,
        max_cycles=1 + args.cycle_minimap2,
        checkpoint_dir=os.path.join(idir, "selfclean_minimap2"),
        resume=bool(getattr(args, "resume", False)),
    )
    # Job 1, phase B: run repeat-aware Winnowmap plus an independent BLASTN
    # pass over every current candidate until that phase also converges.
    winnow_fn = make_self_align_fn(
        os.path.join(idir, "selfclean_winnowmap"), global_args, log,
        blast_mode="independent", broad_aligner="winnowmap", run_blastn=True,
    )
    kept = self_clean_converge(
        coarse, winnow_fn, args.min_segment, "first_genome_winnowmap", log,
        max_cycles=args.max_cycles,
        checkpoint_dir=os.path.join(idir, "selfclean_winnowmap"),
        resume=bool(getattr(args, "resume", False)),
    )
    if not kept:
        raise RuntimeError("first-genome self-clean kept nothing; check inputs/thresholds")
    before_merge = len(kept)
    kept = merge_source_adjacent_pieces(kept)
    log.info("Merged source-adjacent first-reference pieces: %d -> %d records",
             before_merge, len(kept))
    seen_main_ids = set()
    for piece in kept:
        piece.final_id = source_interval_id(
            piece.sample, piece.contig, piece.start, piece.end,
        )
        if piece.final_id in seen_main_ids:
            raise ValueError(f"duplicate coordinate-derived main id: {piece.final_id}")
        seen_main_ids.add(piece.final_id)
    main_fa = os.path.join(outdir, "main_chroms.fa")
    write_fasta_pieces(kept, main_fa, use_final_id=True)
    log.info("Main chroms: %d pieces -> %s", len(kept), main_fa)

    retained_by_contig: Dict[str, List[Tuple[int, int]]] = {}
    for piece in kept:
        retained_by_contig.setdefault(piece.contig, []).append((piece.start, piece.end))
    write_filtered_reference_fasta(
        first.path, first.sample, retained_by_contig,
        os.path.join(outdir, FILTERED_REFERENCE_FASTA), log,
    )

    removed_fastas, removed_bytes = compact_reference_iteration(idir, args.keep_work)
    if removed_fastas:
        log.info(
            "Reference cleanup removed %d intermediate FASTA files (%.2f GiB)",
            removed_fastas, removed_bytes / (1024 ** 3),
        )
    return kept, main_fa, None


def install_precleaned_main(source_fa: str, args, outdir: str,
                            log: logging.Logger) -> Tuple[List[Piece], str, Optional[str]]:
    """Install an already-cleaned main reference without self-aligning it."""
    source = os.path.realpath(source_fa)
    if not os.path.isfile(source):
        raise FileNotFoundError(f"Precleaned reference not found: {source_fa}")
    if source.endswith(".gz"):
        raise ValueError("--precleaned-reference must be an uncompressed FASTA")
    main_fa = os.path.join(outdir, "main_chroms.fa")
    if os.path.realpath(main_fa) == source:
        raise ValueError("Precleaned reference source and output main_chroms.fa are the same file")
    main_fai = main_fa + ".fai"
    source_fai = source + ".fai"
    for destination in (main_fa, main_fai):
        temporary = destination + f".tmp.{os.getpid()}"
        try:
            os.remove(temporary)
        except FileNotFoundError:
            pass
    try:
        fasta_temporary = main_fa + f".tmp.{os.getpid()}"
        shutil.copyfile(source, fasta_temporary)
        os.replace(fasta_temporary, main_fa)

        # The destination is byte-for-byte identical, so an adjacent source
        # .fai has valid offsets after the copy. This also lets explicitly
        # indexed FASTAs with nonuniform physical wrapping remain usable.
        try:
            os.remove(main_fai)
        except FileNotFoundError:
            pass
        if os.path.isfile(source_fai):
            fai_temporary = main_fai + f".tmp.{os.getpid()}"
            shutil.copyfile(source_fai, fai_temporary)
            os.replace(fai_temporary, main_fai)
    finally:
        for destination in (main_fa, main_fai):
            try:
                os.remove(destination + f".tmp.{os.getpid()}")
            except FileNotFoundError:
                pass
    # Validate that downstream random access and FASTA names are usable before
    # starting a long cohort run.
    with IndexedFasta(main_fa) as fa:
        names = fa.names()
        if not names:
            raise ValueError(f"Precleaned reference contains no FASTA records: {source}")
        total = sum(fa.length(name) for name in names)
    log.info("Installed precleaned main reference: %s (%d records, %d bp)",
             source, len(names), total)
    if os.path.isfile(source_fai):
        log.info("Installed trusted precleaned reference index: %s", main_fai)

    return [], main_fa, None


# ---------------------------------------------------------------------------
# Stage 2: insertions from a later assembly
# ---------------------------------------------------------------------------

def extract_insertions(
    assembly: AssemblyEntry, main_hits, min_segment: int, log: logging.Logger,
    blacklist: Optional[BlacklistIndex] = None,
) -> List[Piece]:
    """Novel unaligned segments (vs main) as insertion Pieces, per contig."""
    pieces: List[Piece] = []
    with IndexedFasta(assembly.path) as fa:
        clens = {c: fa.length(c) for c in fa.names()}
        cov = coverage_by_query_id(main_hits, clens)
        counter = 0
        blacklist_removed_bp = 0
        blacklist_split_filtered = 0
        for contig in fa.names():
            seq = fa.sequence(contig)
            segs = find_unaligned_segments(
                cov.get(contig, []), clens[contig], seq, min_segment,
            )
            if not segs:
                continue
            masks = blacklist_overlaps(
                blacklist or {}, contig, segs[0][0], segs[-1][1],
                assembly.sample,
            )
            if masks:
                before_bp = sum(end - start for start, end in segs)
                segs = subtract_intervals(segs, masks)
                blacklist_removed_bp += before_bp - sum(
                    end - start for start, end in segs
                )
            for a, b in segs:
                candidate_seq = seq[a:b]
                # A valid unaligned segment may split into sub-threshold pieces
                # after BED subtraction. Never bridge back across a blacklist.
                if count_unmasked(candidate_seq) < min_segment:
                    blacklist_split_filtered += 1
                    continue
                counter += 1
                pieces.append(Piece(
                    f"ins_{sanitize_id(assembly.sample)}_{sanitize_id(contig)}_{a}_{b}_{counter:08d}",
                    assembly.sample, contig, a, b, "+", candidate_seq,
                    note="insertion"))
    if blacklist_removed_bp or blacklist_split_filtered:
        log.info(
            "%s: blacklist removed %d unaligned bp; %d split fragments failed "
            "the %d-unmasked-bp segment filter",
            assembly.sample, blacklist_removed_bp, blacklist_split_filtered,
            min_segment,
        )
    log.info("%s: extracted %d insertion segments", assembly.sample, len(pieces))
    return pieces


def localize_all(insertions: Sequence[Piece], main_hits, flank: int) -> Dict[str, LocalizeResult]:
    """Localize final insertions after one sorted window/block overlap join."""
    by_contig: Dict[str, list] = {}
    for h in main_hits:
        by_contig.setdefault(h.query_id, []).append(h)
    candidate_contigs = {piece.contig for piece in insertions}
    blocks_by_contig = {
        contig: anchor_blocks_from_hits(by_contig.get(contig, []), contig)
        for contig in candidate_contigs
    }
    localized = localize_windows_from_blocks(
        [(piece.contig, piece.start, piece.end) for piece in insertions],
        blocks_by_contig, flank,
    )
    results: Dict[str, LocalizeResult] = {}
    for piece, result in zip(insertions, localized):
        results[piece.temp_id] = result
    return results


@dataclasses.dataclass
class LocalGroup:
    """New and prior insertions whose extended intervals overlap on one main contig."""
    main: str
    start: int
    end: int
    pieces: List[Piece]
    prior_pieces: List[Piece] = dataclasses.field(default_factory=list)


def split_lift_placements(
    pieces: Sequence[Piece], lifted_by_temp_id: Dict[str, LocalizeResult],
) -> Tuple[List[Piece], Dict[str, LocalizeResult]]:
    """Make one independently cleaned candidate for every valid lift placement."""
    expanded: List[Piece] = []
    by_lift: Dict[str, LocalizeResult] = {}
    for piece in pieces:
        result = lifted_by_temp_id[piece.temp_id]
        if len(result.placements) <= 1:
            expanded.append(piece)
            by_lift[piece.lift_id] = result
            continue
        total = len(result.placements)
        for index, placement in enumerate(result.placements, 1):
            lift_id = f"{piece.lift_id}_placement{index:02d}"
            clone = dataclasses.replace(
                piece,
                temp_id=f"{piece.temp_id}_placement{index:02d}",
                lift_id=lift_id,
                final_id="",
            )
            expanded.append(clone)
            note = ";".join(filter(None, [result.note, f"placement_{index}_of_{total}"]))
            by_lift[lift_id] = LocalizeResult(result.status, [placement], note=note)
    return expanded, by_lift


def _extended_placement(
    placement, main_lengths: Dict[str, int], extension: int,
) -> Optional[Tuple[str, int, int]]:
    length = main_lengths.get(placement.main)
    if length is None:
        return None
    start = max(0, min(placement.start, placement.end) - extension)
    end = min(length, max(placement.start, placement.end) + extension)
    if end <= start:
        return None
    return placement.main, start, end


def group_localized_insertions(
    pieces: Sequence[Piece], localized_by_lift: Dict[str, LocalizeResult],
    main_lengths: Dict[str, int], extension: int,
    placement_catalog: Sequence[Tuple[Piece, LocalizeResult]] = (),
    max_group_size: int = 1_000_000,
) -> Tuple[List[LocalGroup], List[Piece]]:
    """Merge extended new/prior placement intervals, chromosome-wise.

    A dictionary per main contig plus one coordinate sort makes this O(n log n),
    rather than repeatedly comparing every pair of insertion intervals.
    Every current placement-specific candidate is cleaned once. Prior accepted
    insertions participate in both interval expansion and local-reference content.
    """
    # Entry: start, end, current piece (or None), prior piece (or None).
    by_main: Dict[str, List[Tuple[int, int, Optional[Piece], Optional[Piece]]]] = {}
    unmapped: List[Piece] = []
    for piece in pieces:
        res = localized_by_lift.get(piece.lift_id)
        if res is None or not res.placements:
            unmapped.append(piece)
            continue
        extended = _extended_placement(res.placements[0], main_lengths, extension)
        if extended is None:
            unmapped.append(piece)
            continue
        main, start, end = extended
        by_main.setdefault(main, []).append((start, end, piece, None))

    for prior_piece, result in placement_catalog:
        for placement in result.placements:
            extended = _extended_placement(placement, main_lengths, extension)
            if extended is None:
                continue
            main, start, end = extended
            by_main.setdefault(main, []).append((start, end, None, prior_piece))

    groups: List[LocalGroup] = []
    for main in sorted(by_main):
        intervals = sorted(by_main[main], key=lambda x: (x[0], x[1]))
        current: Optional[LocalGroup] = None
        for start, end, piece, prior_piece in intervals:
            exceeds_cap = (current is not None
                           and max(current.end, end) - current.start > max_group_size)
            if current is None or start > current.end or exceeds_cap:
                if current is not None and current.pieces:
                    groups.append(current)
                current = LocalGroup(
                    main, start, end,
                    [piece] if piece is not None else [],
                    [prior_piece] if prior_piece is not None else [],
                )
            else:
                current.end = max(current.end, end)
                if piece is not None:
                    current.pieces.append(piece)
                if prior_piece is not None:
                    current.prior_pieces.append(prior_piece)
        if current is not None and current.pieces:
            groups.append(current)
    return groups, unmapped


def write_local_reference(
    main_fa: str, group: LocalGroup, output_fa: str,
) -> Tuple[List[Piece], Dict[str, Tuple[str, int]]]:
    """Write a main slice plus prior insertions and its coordinate translation."""
    Path(output_fa).parent.mkdir(parents=True, exist_ok=True)
    nearby_by_id: Dict[str, Piece] = {}
    for piece in group.prior_pieces:
        nearby_by_id.setdefault(piece.final_id or piece.temp_id, piece)
    nearby = list(nearby_by_id.values())
    with IndexedFasta(main_fa) as main, open(output_fa, "w") as out:
        sequence = main.fetch(group.main, group.start, group.end)
        target_id = sanitize_id(f"local_{group.main}_{group.start}_{group.end}")
        out.write(f">{target_id}\n{wrap_fasta(sequence)}\n")
        seen = {target_id, group.main}
        for piece in nearby:
            piece_id = piece.final_id or piece.temp_id
            if piece_id in seen:
                continue
            seen.add(piece_id)
            out.write(f">{piece_id}\n{wrap_fasta(piece.seq)}\n")
    return nearby, {target_id: (group.main, group.start)}


def write_multi_placement_reference(
    main_fa: str, result: LocalizeResult, extension: int, output_fa: str,
) -> Dict[str, Tuple[str, int]]:
    """Write merged main slices surrounding every valid placement of one locus."""
    Path(output_fa).parent.mkdir(parents=True, exist_ok=True)
    with IndexedFasta(main_fa) as main:
        lengths = {name: main.length(name) for name in main.names()}
        by_main: Dict[str, List[Tuple[int, int]]] = {}
        for placement in result.placements:
            extended = _extended_placement(placement, lengths, extension)
            if extended is None:
                continue
            name, start, end = extended
            by_main.setdefault(name, []).append((start, end))
        coordinate_map: Dict[str, Tuple[str, int]] = {}
        with open(output_fa, "w") as out:
            counter = 0
            for name in sorted(by_main):
                for start, end in merge_intervals(by_main[name]):
                    counter += 1
                    target_id = sanitize_id(f"local_{name}_{start}_{end}_{counter}")
                    out.write(f">{target_id}\n{wrap_fasta(main.fetch(name, start, end))}\n")
                    coordinate_map[target_id] = (name, start)
    if not coordinate_map:
        raise RuntimeError("Mapped novel-locus candidate has no valid main-reference slice")
    return coordinate_map


def qualifies_unmapped_main(piece: Piece, size_cutoff: int,
                            min_unmasked: int, min_fraction: float) -> bool:
    """Whether an unplaced sequence is substantial enough to become a main locus."""
    required_unmasked = max(float(min_unmasked), min_fraction * piece.length)
    return piece.length > size_cutoff and piece.unmasked >= required_unmasked


def qualifies_preclean_candidate(piece: Piece, min_segment: int) -> bool:
    """Retain a candidate before global cleaning using only min_segment."""
    return piece.unmasked >= min_segment


def qualifies_mapped_promotion(piece: Piece, promote_size: int,
                               min_unmasked: int, min_fraction: float) -> bool:
    """Whether a mapped sequence is large and sufficiently unmasked to promote."""
    required_unmasked = max(float(min_unmasked), min_fraction * piece.length)
    return piece.unmasked >= promote_size and piece.unmasked >= required_unmasked


def qualifies_promotion(piece: Piece, result: LocalizeResult, args) -> bool:
    if result.placements:
        return qualifies_mapped_promotion(
            piece, args.promote_size, args.min_segment, args.min_unmasked_fraction,
        )
    return qualifies_unmapped_main(
        piece, args.unmapped_main_size, args.min_segment, args.min_unmasked_fraction,
    )


def build_cohort_windows(
    records: Sequence[Tuple[Piece, LocalizeResult]],
    promoted_catalog: Sequence[Tuple[Piece, LocalizeResult]],
    main_lengths: Dict[str, int], extension: int = 3_000,
    preferred_size: int = 100_000,
) -> Tuple[List[LocalGroup], List[Piece]]:
    """Build dynamic cohort windows without splitting any extended insertion.

    Intervals are sorted per main contig. At an unspanned gap, start a new
    window when bridging to the next interval would make the current coordinate
    span exceed ``preferred_size``. Overlapping chains are never forcibly cut.
    """
    by_main: Dict[str, List[Tuple[int, int, Optional[Piece], Optional[Piece]]]] = {}
    unmapped: List[Piece] = []
    for piece, result in records:
        if not result.placements:
            unmapped.append(piece)
            continue
        # Discovery already split conflicting placements into separate records.
        extended = _extended_placement(result.placements[0], main_lengths, extension)
        if extended is None:
            unmapped.append(piece)
            continue
        main, start, end = extended
        by_main.setdefault(main, []).append((start, end, piece, None))
    for prior_piece, result in promoted_catalog:
        for placement in result.placements:
            extended = _extended_placement(placement, main_lengths, extension)
            if extended is None:
                continue
            main, start, end = extended
            by_main.setdefault(main, []).append((start, end, None, prior_piece))

    windows: List[LocalGroup] = []
    for main in sorted(by_main):
        intervals = sorted(by_main[main], key=lambda item: (item[0], item[1]))
        current: Optional[LocalGroup] = None
        for start, end, piece, prior_piece in intervals:
            gap = current is not None and start > current.end
            crosses_preferred = (current is not None
                                 and max(current.end, end) - current.start > preferred_size)
            if current is None or (gap and crosses_preferred):
                if current is not None and current.pieces:
                    windows.append(current)
                current = LocalGroup(
                    main, start, end,
                    [piece] if piece is not None else [],
                    [prior_piece] if prior_piece is not None else [],
                )
                continue
            current.end = max(current.end, end)
            if piece is not None:
                current.pieces.append(piece)
            if prior_piece is not None:
                current.prior_pieces.append(prior_piece)
        if current is not None and current.pieces:
            windows.append(current)
    return windows, unmapped


def classify_insertion(piece: Piece, promote_size: int) -> str:
    """Return the destination pool for a kept insertion."""
    return "main" if piece.unmasked >= promote_size else "small"


# ---------------------------------------------------------------------------
# Output table
# ---------------------------------------------------------------------------

TABLE_HEADER = ("insertion_id\tsample\tsource_contig\tsource_start\tsource_end\t"
                "seq_size\tunmasked\tpool\tstatus\tmain_chrom\tmain_start\tmain_end\tstrand\t"
                "representative_id\trepresentative_start\trepresentative_end\t"
                "representative_strand\tnote\n")


def placement_fields(res: Optional[LocalizeResult]) -> Tuple[str, str, str, str, str]:
    """Return (status, main_chrom, main_start, main_end, strand) for the primary placement."""
    if res is None or not res.placements:
        status = res.status if res else "unmapped"
        return status, ".", ".", ".", "."
    pl = res.placements[0]
    extra = ""
    if len(res.placements) > 1:
        extra = ";".join(f"{q.main}:{q.start}-{q.end}{q.strand}" for q in res.placements[1:])
    main = pl.main if not extra else f"{pl.main}|alt={extra}"
    return res.status, main, str(pl.start), str(pl.end), pl.strand


def table_row(piece: Piece, pool: str, res: Optional[LocalizeResult]) -> str:
    status, main, ms, me, strand = placement_fields(res)
    note = res.note if res else ""
    return "\t".join([
        piece.final_id or piece.temp_id, piece.sample, piece.contig,
        str(piece.start), str(piece.end), str(piece.length), str(piece.unmasked),
        pool, status, main, ms, me, strand, ".", ".", ".", ".", note or ".",
    ]) + "\n"


def redundant_row(piece: Piece, representative: Optional[Tuple[str, int, int, str]] = None) -> str:
    """A redundant insertion: record its new-assembly query coordinates."""
    if representative is None:
        rid, rs, re, rstrand = ".", ".", ".", "."
    else:
        rid, start, end, rstrand = representative
        rs, re = str(start), str(end)
    return "\t".join([
        piece.final_id or piece.temp_id, piece.sample, piece.contig,
        str(piece.start), str(piece.end), str(piece.length), str(piece.unmasked),
        ".", "redundant", ".", ".", ".", ".", rid, rs, re, rstrand,
        f"query:{piece.source_coord()}",
    ]) + "\n"


def _best_pool_representative(
    piece: Piece, hits,
    target_coordinate_map: Optional[Dict[str, Tuple[str, int]]] = None,
) -> Optional[Tuple[str, int, int, str]]:
    """Best existing representative for a query piece, using unique covered bp."""
    candidates = []
    for hit in hits:
        if hit.query_id != piece.temp_id:
            continue
        qcov = sum(b - a for a, b in merge_intervals(hit.query_intervals))
        if qcov <= 0 or not hit.target_intervals:
            continue
        ta = min(a for a, _ in hit.target_intervals)
        tb = max(b for _, b in hit.target_intervals)
        target = hit.target_id
        if target_coordinate_map and target in target_coordinate_map:
            target, offset = target_coordinate_map[target]
            ta += offset
            tb += offset
        candidates.append((qcov, hit.identity, target, ta, tb, hit.strand))
    if not candidates:
        return None
    _, _, target, a, b, strand = max(candidates)
    return target, a, b, strand


def exclude_against_reference_converge(
    pieces: Sequence[Piece], reference_fa: str, workdir: str, args,
    log: logging.Logger,
    representatives: Dict[str, Tuple[str, int, int, str]],
    blast_mode: str = "independent",
    alignment_sink: Optional[List] = None,
    winnow_rep_kmers: Optional[str] = None,
    broad_aligner: str = "winnowmap",
    run_blastn: bool = True,
    max_cycles: Optional[int] = None,
    target_coordinate_map: Optional[Dict[str, Tuple[str, int]]] = None,
    warn_on_limit: bool = True,
    unlimited: bool = False,
) -> List[Piece]:
    """Subtract represented intervals and repeat until nothing changes.

    If ``alignment_sink`` is supplied, hits are translated from piece-local
    coordinates back to the original assembly contig and appended. This lets
    newly exposed main alignments participate in later insertion localization.
    """
    current = list(pieces)
    if (broad_aligner == "winnowmap" and winnow_rep_kmers is None
            and not getattr(args, "no_repkmers", False)):
        winnow_rep_kmers = getattr(args, "winnow_rep_kmers", None)
    cycle_limit = args.max_cycles if max_cycles is None else max_cycles
    converged = False
    cycle = 0
    while unlimited or cycle < cycle_limit:
        if not current:
            converged = True
            break
        cdir = os.path.join(workdir, f"cycle{cycle:02d}")
        query_fa = os.path.join(cdir, "candidates.fa")
        write_fasta_pieces(current, query_fa)
        write_piece_map(current, os.path.join(cdir, "candidates.map.tsv"))
        hits = collect_alignment_hits(
            query_fa, reference_fa, cdir, args.threads, args.min_identity,
            args.min_segment, log, winnow_params=WINNOW_PARAMS, rep_kmers=winnow_rep_kmers,
            blast_word_size=args.blast_word_size, blast_evalue=args.blast_evalue,
            blast_max_target_seqs=args.blast_max_target_seqs, prefix="to_pool",
            skip_blastn=args.skip_blastn or not run_blastn, blast_mode=blast_mode,
            broad_aligner=broad_aligner, minimap_params=MINIMAP_PARAMS,
            blast_db_prefix=os.path.join(workdir, "shared_blastdb", "ref"))
        if alignment_sink is not None:
            by_id = {p.temp_id: p for p in current}
            for h in hits:
                p = by_id.get(h.query_id)
                if p is None:
                    continue
                alignment_sink.append(dataclasses.replace(
                    h,
                    query_id=p.contig,
                    query_intervals=[(p.start + a, p.start + b) for a, b in h.query_intervals],
                    query_start=p.start + h.query_start if h.query_start >= 0 else -1,
                    query_end=p.start + h.query_end if h.query_end >= 0 else -1,
                    aligned_pairs=[(p.start + q0, p.start + q1, t0, t1)
                                   for q0, q1, t0, t1 in h.aligned_pairs],
                ))
        cov = coverage_by_query_id(hits, {p.temp_id: p.length for p in current})
        hits_by_query: Dict[str, list] = {}
        for hit in hits:
            hits_by_query.setdefault(hit.query_id, []).append(hit)
        next_pieces: List[Piece] = []
        before = sum(p.length for p in current)
        for p in current:
            segs = find_unaligned_segments(cov.get(p.temp_id, []), p.length, p.seq,
                                           args.min_segment)
            if not segs:
                rep = _best_pool_representative(
                    p, hits_by_query.get(p.temp_id, ()), target_coordinate_map,
                )
                if rep is not None:
                    representatives[p.origin_id] = rep
            for idx, (a, b) in enumerate(segs, 1):
                next_pieces.append(Piece(
                    temp_id=f"{p.temp_id}_poolnr{idx:04d}", sample=p.sample,
                    contig=p.contig, start=p.start + a, end=p.start + b,
                    strand=p.strand, seq=p.seq[a:b], parent_id=p.temp_id,
                    note="nonredundant_vs_small_pool", origin_id=p.origin_id,
                    lift_id=p.lift_id,
                    within_sample_rank=p.within_sample_rank,
                ))
        removed = before - sum(p.length for p in next_pieces)
        log.info("Pool exclusion cycle %d: %d -> %d pieces, removed %d bp",
                 cycle, len(current), len(next_pieces), removed)
        current = next_pieces
        if removed == 0:
            converged = True
            break
        cycle += 1
    if current and not converged and warn_on_limit and not unlimited:
        log.warning("Reference exclusion reached max_cycles=%d before zero removal",
                    cycle_limit)
    return current


def clean_against_complete_reference(
    pieces: Sequence[Piece],
    minimap_target: str,
    reference_fasta: str,
    workdir: str,
    args,
    log: logging.Logger,
    label: str = "Post-self-clean full-reference exclusion",
) -> List[Piece]:
    """Remove all sequence represented by the complete first assembly.

    This accuracy pass is deliberately after cohort-wide candidate self-cleaning.
    It first runs residual Minimap2 exclusion until no additional coverage is
    found, then runs Winnowmap plus an independent BLASTN pass until the same
    convergence condition. The complete first-assembly FASTA—not its cleaned or
    windowed derivative—is the biological reference for both phases.
    """
    if not pieces:
        return []
    if not os.path.isfile(reference_fasta):
        raise FileNotFoundError(reference_fasta)
    global_args = argparse.Namespace(**vars(args))
    global_args.threads = _global_threads(args)
    representatives: Dict[str, Tuple[str, int, int, str]] = {}
    log.info(
        "%s: %d candidates; %d threads; both phases run until convergence",
        label, len(pieces), global_args.threads,
    )
    minimap_dir = os.path.join(workdir, "minimap2")
    minimap_residuals = exclude_against_reference_converge(
        pieces, minimap_target, minimap_dir, global_args, log,
        representatives, blast_mode="residual", broad_aligner="minimap2",
        run_blastn=False, warn_on_limit=False, unlimited=True,
    )
    log.info(
        "%s: Minimap2 retained %d/%d pieces",
        label, len(minimap_residuals), len(pieces),
    )
    if not minimap_residuals:
        return []
    winnow_dir = os.path.join(workdir, "winnowmap_blastn")
    cleaned = exclude_against_reference_converge(
        minimap_residuals, reference_fasta, winnow_dir, global_args, log,
        representatives, blast_mode="independent", broad_aligner="winnowmap",
        run_blastn=True,
        winnow_rep_kmers=getattr(global_args, "winnow_rep_kmers", None),
        warn_on_limit=False, unlimited=True,
    )
    log.info(
        "%s: Winnowmap+BLASTN retained %d/%d Minimap2 residual pieces",
        label, len(cleaned), len(minimap_residuals),
    )
    return cleaned


def exclude_against_pool_converge(pieces: Sequence[Piece], pool_fa: str, workdir: str,
                                  args, log: logging.Logger,
                                  representatives: Dict[str, Tuple[str, int, int, str]],
                                  rep_kmers: Optional[str] = None) -> List[Piece]:
    """Backward-compatible named wrapper for independent small-pool exclusion."""
    return exclude_against_reference_converge(
        pieces, pool_fa, workdir, args, log, representatives,
        blast_mode="independent", broad_aligner="winnowmap",
        winnow_rep_kmers=rep_kmers,
    )


# ---------------------------------------------------------------------------
# Stage 2 driver
# ---------------------------------------------------------------------------

def process_assembly(assembly: AssemblyEntry, iteration: int, main_fa: str, small_fa: str,
                     small_pieces: List[Piece],
                     placement_catalog: List[Tuple[Piece, LocalizeResult]], args, outdir: str,
                     counters: Dict[str, int], table_path: str, log: logging.Logger) -> Tuple[List[Piece], List[Piece]]:
    """Return (promoted_main_pieces, added_small_pieces)."""
    idir = os.path.join(outdir, "iterations", f"{iteration:03d}_{sanitize_id(assembly.sample)}")
    Path(idir).mkdir(parents=True, exist_ok=True)

    # Job 2, coarse pass: minimap2 removes the easy main-chromosome matches.
    # BLASTN is intentionally deferred until Winnowmap sees the reduced residual.
    minimap_hits = collect_alignment_hits(
        assembly.path, main_fa, os.path.join(idir, "to_main"), args.threads,
        args.min_identity, args.min_segment, log, winnow_params=WINNOW_PARAMS,
        blast_word_size=args.blast_word_size, blast_evalue=args.blast_evalue,
        blast_max_target_seqs=args.blast_max_target_seqs, prefix="to_main", skip_blastn=True,
        blast_mode="residual", broad_aligner="minimap2", minimap_params=MINIMAP_PARAMS)

    # Step b: novel segments.
    insertions = extract_insertions(
        assembly, minimap_hits, args.min_segment, log,
        getattr(args, "blacklist_intervals", None),
    )
    if not insertions:
        log.info("%s: no insertions", assembly.sample)
        prune_old_iterations(outdir, idir, args, log)
        return [], []

    # Job 2, phase A continuation: after the required whole-assembly minimap2
    # pass above, run exactly --cycle-minimap2 residual minimap2 cycles.
    main_representatives: Dict[str, Tuple[str, int, int, str]] = {}
    if args.cycle_minimap2:
        minimap_residuals = exclude_against_reference_converge(
            insertions, main_fa, os.path.join(idir, "to_main_minimap2_converge"), args, log,
            main_representatives, blast_mode="residual", alignment_sink=minimap_hits,
            broad_aligner="minimap2", run_blastn=False,
            max_cycles=args.cycle_minimap2,
        )
    else:
        minimap_residuals = insertions

    # Lift before expensive repeat-aware alignment, using only minimap2 evidence
    # accumulated from the full assembly and every requested residual cycle.
    for lift_index, piece in enumerate(minimap_residuals):
        piece.lift_id = f"lift{lift_index:08d}"
    lifted = localize_all(minimap_residuals, minimap_hits, args.flank)
    minimap_residuals, localized_by_lift = split_lift_placements(
        minimap_residuals, lifted,
    )
    with IndexedFasta(main_fa) as main_reader:
        main_lengths = {name: main_reader.length(name) for name in main_reader.names()}
    local_groups, unmapped = group_localized_insertions(
        minimap_residuals, localized_by_lift, main_lengths, args.local_extension,
        placement_catalog, args.max_local_group_size,
    )
    log.info("%s: %d mapped local groups, %d unmapped candidates",
             assembly.sample, len(local_groups), len(unmapped))

    # Winnowmap/BLASTN now operate only inside merged local regions. Each local
    # reference contains the main slice and previously accepted insertions at
    # that locus, combining main-reference and insertion-pool exclusion.
    pool_representatives: Dict[str, Tuple[str, int, int, str]] = {}
    self_provenance: Dict[str, set] = {}
    survivors: List[Piece] = []
    local_group_threads = _parallel_task_threads(args, len(local_groups))

    def clean_local_group(item):
        group_index, group = item
        local_args = argparse.Namespace(**vars(args))
        local_args.threads = local_group_threads
        local_representatives: Dict[str, Tuple[str, int, int, str]] = {}
        local_provenance: Dict[str, set] = {}
        group_dir = os.path.join(idir, "local_groups", f"group{group_index:06d}")
        local_fa = os.path.join(group_dir, "local_reference.fa")
        nearby, target_coordinate_map = write_local_reference(main_fa, group, local_fa)
        log.info("Local group %s:%d-%d: %d candidates, %d prior insertions",
                 group.main, group.start, group.end, len(group.pieces), len(nearby))
        local_rep_kmers = getattr(args, "winnow_rep_kmers", None)
        local_survivors = exclude_against_reference_converge(
            group.pieces, local_fa, os.path.join(group_dir, "to_local"), local_args, log,
            local_representatives, blast_mode="independent", broad_aligner="winnowmap",
            winnow_rep_kmers=local_rep_kmers,
            target_coordinate_map=target_coordinate_map,
        )
        self_fn = make_self_align_fn(
            os.path.join(group_dir, "selfclean"), local_args, log,
            blast_mode="independent", broad_aligner="winnowmap",
        )
        local_survivors = self_clean_converge(
            local_survivors, self_fn, local_args.min_segment, "insertion", log,
            max_cycles=local_args.max_cycles, provenance=local_provenance,
            sample_priorities=getattr(local_args, "sample_priorities", None),
        )
        return group_index, local_survivors, local_representatives, local_provenance

    if local_groups:
        workers = min(_local_jobs(args), len(local_groups))
        log.info("Cleaning %d local groups with %d parallel jobs x %d threads",
                 len(local_groups), workers, local_group_threads)
        with ThreadPoolExecutor(max_workers=workers, thread_name_prefix="local-group") as executor:
            group_results = list(executor.map(clean_local_group, enumerate(local_groups)))
        # executor.map preserves input order, keeping output IDs deterministic.
        for _index, local_survivors, local_representatives, local_provenance in group_results:
            survivors.extend(local_survivors)
            pool_representatives.update(local_representatives)
            for origin, links in local_provenance.items():
                self_provenance.setdefault(origin, set()).update(links)

    # Unmapped candidates are not compared or added to the small pool now.
    # Only substantial, sufficiently unmasked sequences become new main loci.
    unmapped_main = [piece for piece in unmapped if qualifies_unmapped_main(
        piece, args.unmapped_main_size, args.min_segment, args.min_unmasked_fraction,
    )]
    unmapped_main_objects = {id(piece) for piece in unmapped_main}
    ignored_origins = {
        piece.origin_id for piece in unmapped if id(piece) not in unmapped_main_objects
    }
    log.info("%s: promoting %d unmapped main loci; ignoring %d smaller/masked unmapped candidates",
             assembly.sample, len(unmapped_main), len(unmapped) - len(unmapped_main))
    survivors.extend(unmapped_main)
    unmapped_main_lifts = {piece.lift_id for piece in unmapped_main}

    promoted: List[Piece] = []
    added_small: List[Piece] = []
    kept_orig = {p.origin_id for p in survivors}
    # Assign all final IDs before resolving self-clean provenance so a redundant
    # original can point to the actual surviving representative record.
    for p in survivors:
        pool = "main" if p.lift_id in unmapped_main_lifts else classify_insertion(p, args.promote_size)
        if pool == "main":
            assign_final_ids([p], f"mainins_{assembly.sample}", counters)
            promoted.append(p)
        else:
            assign_final_ids([p], f"ins_{assembly.sample}", counters)
            added_small.append(p)
    survivors_by_origin: Dict[str, List[Piece]] = {}
    for p in survivors:
        survivors_by_origin.setdefault(p.origin_id, []).append(p)

    def resolve_self_representative(origin: str) -> Optional[Tuple[str, int, int, str]]:
        pending = deque(self_provenance.get(origin, set()))
        seen = {origin}
        while pending:
            item = pending.popleft()
            rep_origin, abs_start, abs_end, rep_strand = item
            if rep_origin in seen:
                continue
            seen.add(rep_origin)
            reps = survivors_by_origin.get(rep_origin, [])
            if reps:
                overlapping = [q for q in reps if q.start < abs_end and q.end > abs_start]
                rep = max(overlapping or reps, key=lambda q: q.length)
                a = max(rep.start, abs_start) - rep.start
                b = min(rep.end, abs_end) - rep.start
                if b <= a:
                    a, b = 0, rep.length
                return rep.final_id or rep.temp_id, a, b, rep_strand
            pending.extend(self_provenance.get(rep_origin, set()))
        return None

    with open(table_path, "a") as tbl:
        for p in survivors:
            res = localized_by_lift.get(p.lift_id, LocalizeResult("unmapped", []))
            pool = "main" if p.lift_id in unmapped_main_lifts else classify_insertion(p, args.promote_size)
            tbl.write(table_row(p, pool, res))
            placement_catalog.append((p, res))
        for ins in insertions:  # originals with no surviving sequence are redundant
            if ins.origin_id not in kept_orig and ins.origin_id not in ignored_origins:
                representative = main_representatives.get(ins.origin_id)
                if representative is None:
                    representative = pool_representatives.get(ins.origin_id)
                if representative is None:
                    representative = resolve_self_representative(ins.origin_id)
                tbl.write(redundant_row(ins, representative))

    # Grow the pools on disk.
    if added_small:
        _append_pieces(small_fa, added_small)
        small_pieces.extend(added_small)
    if promoted:
        _append_pieces(main_fa, promoted)
        log.info("%s: promoted %d insertions (>=%d unmasked) into main", assembly.sample,
                 len(promoted), args.promote_size)
    log.info("%s: +%d small, +%d main", assembly.sample, len(added_small), len(promoted))
    prune_old_iterations(outdir, idir, args, log)
    return promoted, added_small


def discover_assembly(
    assembly: AssemblyEntry, iteration: int, main_fa: str, args, outdir: str,
    counters: Dict[str, int], table_path: str, log: logging.Logger,
) -> Tuple[List[Tuple[Piece, LocalizeResult]], List[Tuple[Piece, LocalizeResult]]]:
    """Minimap2 discovery plus rare high-quality cleaning of promotion candidates.

    Returns ``(promoted_records, mapped_records_for_cohort_cleanup)``.
    """
    idir = os.path.join(outdir, "iterations", f"{iteration:03d}_{sanitize_id(assembly.sample)}")
    Path(idir).mkdir(parents=True, exist_ok=True)
    marker_payload = {
        "kind": "assembly", "iteration": iteration, "sample": assembly.sample,
        "path": os.path.realpath(assembly.path),
        "size": os.path.getsize(assembly.path),
        "mtime_ns": os.stat(assembly.path).st_mtime_ns,
        "config_signature": _restart_signature(args),
    }
    minimap_hits = collect_alignment_hits(
        assembly.path, main_fa, os.path.join(idir, "to_main"), args.threads,
        args.min_identity, args.min_segment, log, winnow_params=WINNOW_PARAMS,
        blast_word_size=args.blast_word_size, blast_evalue=args.blast_evalue,
        blast_max_target_seqs=args.blast_max_target_seqs, prefix="to_main", skip_blastn=True,
        blast_mode="residual", broad_aligner="minimap2", minimap_params=MINIMAP_PARAMS,
    )
    insertions = extract_insertions(
        assembly, minimap_hits, args.min_segment, log,
        getattr(args, "blacklist_intervals", None),
    )
    if not insertions:
        write_discovery_bundle(idir, [], [], checkpoint_records=[])
        _write_stage_results(idir, [])
        compact_discovery_iteration(idir, args.keep_work)
        _write_stage_marker(idir, "_READY", marker_payload)
        _commit_stage(
            idir, [], main_fa, os.path.join(outdir, "small_insertions.fa"),
            table_path, marker_payload,
        )
        return [], []

    minimap_residuals = insertions
    if args.cycle_minimap2:
        minimap_residuals = exclude_against_reference_converge(
            insertions, main_fa, os.path.join(idir, "to_main_minimap2_converge"), args, log,
            {}, blast_mode="residual", alignment_sink=minimap_hits,
            broad_aligner="minimap2", run_blastn=False,
            max_cycles=args.cycle_minimap2,
        )
    safe_sample = sanitize_id(assembly.sample)
    for lift_index, piece in enumerate(minimap_residuals):
        piece.lift_id = f"lift_{iteration:04d}_{safe_sample}_{lift_index:08d}"
    lifted = localize_all(minimap_residuals, minimap_hits, args.flank)
    discovered, localized_by_lift = split_lift_placements(minimap_residuals, lifted)
    all_discovered_records = [
        (dataclasses.replace(piece), localized_by_lift[piece.lift_id]) for piece in discovered
    ]

    promotion_candidates = [
        piece for piece in discovered
        if qualifies_promotion(piece, localized_by_lift[piece.lift_id], args)
    ]
    candidate_ids = {id(piece) for piece in promotion_candidates}
    cohort_candidates = [piece for piece in discovered if id(piece) not in candidate_ids]

    cleaned_large: List[Piece] = []
    if promotion_candidates:
        log.info("%s: high-quality Winnowmap+BLASTN cleaning for %d promotion candidates",
                 assembly.sample, len(promotion_candidates))
        main_rep_kmers = getattr(args, "winnow_rep_kmers", None)
        cleaned_large = exclude_against_reference_converge(
            promotion_candidates, main_fa, os.path.join(idir, "promotion_vs_main"), args, log,
            {}, blast_mode="residual", broad_aligner="winnowmap",
            winnow_rep_kmers=main_rep_kmers, run_blastn=True,
        )

    # Re-check after main cleaning; fragments that fall below promotion quality
    # join mapped cohort cleanup rather than being added prematurely.
    still_large: List[Piece] = []
    for piece in cleaned_large:
        result = localized_by_lift[piece.lift_id]
        if qualifies_promotion(piece, result, args):
            still_large.append(piece)
        elif result.placements:
            cohort_candidates.append(piece)

    # Remove redundancy among new large candidates only within the same mapped
    # locus. Unmapped large candidates are cleaned together because no locus is
    # available to distinguish occurrences.
    with IndexedFasta(main_fa) as main_reader:
        main_lengths = {name: main_reader.length(name) for name in main_reader.names()}
    large_groups, large_unmapped = group_localized_insertions(
        still_large, localized_by_lift, main_lengths, args.local_extension,
        max_group_size=sys.maxsize,
    )

    def selfclean_large_group(item):
        index, pieces = item
        local_args = argparse.Namespace(**vars(args))
        local_args.threads = promotion_group_threads
        align_fn = make_self_align_fn(
            os.path.join(idir, "promotion_selfclean", f"group{index:06d}"),
            local_args, log, blast_mode="independent", broad_aligner="winnowmap",
        )
        return self_clean_converge(
            pieces, align_fn, local_args.min_segment, "insertion", log,
            max_cycles=local_args.max_cycles,
            sample_priorities=getattr(local_args, "sample_priorities", None),
        )

    promotion_groups = [group.pieces for group in large_groups]
    if large_unmapped:
        promotion_groups.append(large_unmapped)
    promotion_group_threads = _parallel_task_threads(args, len(promotion_groups))
    cleaned_promotions: List[Piece] = []
    if promotion_groups:
        workers = min(_local_jobs(args), len(promotion_groups))
        with ThreadPoolExecutor(max_workers=workers, thread_name_prefix="promotion-clean") as executor:
            for pieces in executor.map(selfclean_large_group, enumerate(promotion_groups)):
                cleaned_promotions.extend(pieces)

    promoted_records: List[Tuple[Piece, LocalizeResult]] = []
    for piece in cleaned_promotions:
        result = localized_by_lift[piece.lift_id]
        if qualifies_promotion(piece, result, args):
            assign_final_ids([piece], f"mainins_{assembly.sample}", counters)
            promoted_records.append((piece, result))
        elif result.placements:
            cohort_candidates.append(piece)

    # Give cohort candidates stable globally unique FASTA names after any
    # iterative cleaner has temporarily re-based their ids.
    cohort_records: List[Tuple[Piece, LocalizeResult]] = []
    for index, piece in enumerate(cohort_candidates):
        piece.temp_id = f"disc_{iteration:04d}_{safe_sample}_{index:08d}"
        result = localized_by_lift[piece.lift_id]
        if result.placements:
            cohort_records.append((piece, result))

    promoted_pieces = [piece for piece, _ in promoted_records]
    checkpoint_records = (
        [("promoted", piece, result) for piece, result in promoted_records]
        + [("cohort", piece, result) for piece, result in cohort_records]
    )
    result_rows = [table_row(piece, "main", result)
                   for piece, result in promoted_records]
    # Persist every minimap-discovered insertion and its PAF evidence, plus the
    # exact post-cleaning records needed to continue without rerunning this
    # assembly.  Only after READY exists do we mutate cumulative outputs.
    write_discovery_bundle(
        idir, all_discovered_records, promoted_pieces,
        checkpoint_records=checkpoint_records,
    )
    _write_stage_results(idir, result_rows)
    compact_discovery_iteration(idir, args.keep_work)
    _write_stage_marker(idir, "_READY", marker_payload)
    _commit_stage(
        idir, checkpoint_records, main_fa, os.path.join(outdir, "small_insertions.fa"),
        table_path, marker_payload,
    )
    log.info("%s discovery: %d promoted, %d mapped candidates queued",
             assembly.sample, len(promoted_records), len(cohort_records))
    return promoted_records, cohort_records


def clean_cohort_insertions(
    records: Sequence[Tuple[Piece, LocalizeResult]],
    promoted_catalog: Sequence[Tuple[Piece, LocalizeResult]],
    main_fa: str, small_fa: str, args, outdir: str, counters: Dict[str, int],
    table_path: str, log: logging.Logger,
) -> Tuple[List[Piece], List[Piece]]:
    """Clean all mapped small-insertion candidates together in dynamic windows."""
    idir = os.path.join(outdir, "iterations", "cohort_cleanup")
    Path(idir).mkdir(parents=True, exist_ok=True)
    marker_payload = {
        "kind": "cohort_cleanup", "record_count": len(records),
        "config_signature": _restart_signature(args),
    }
    if not records:
        write_discovery_bundle(idir, [], [], checkpoint_records=[])
        _write_stage_results(idir, [])
        compact_discovery_iteration(idir, args.keep_work)
        _write_stage_marker(idir, "_READY", marker_payload)
        _commit_stage(idir, [], main_fa, small_fa, table_path, marker_payload)
        return [], []
    with IndexedFasta(main_fa) as main_reader:
        main_lengths = {name: main_reader.length(name) for name in main_reader.names()}
    windows, unmapped = build_cohort_windows(
        records, promoted_catalog, main_lengths,
        extension=args.local_extension,
        preferred_size=args.cohort_window_size,
    )
    if unmapped:
        log.info("Cohort cleanup ignores %d unplaced candidates for later handling", len(unmapped))
    log.info("Cohort cleanup: %d mapped candidates in %d dynamic windows",
             len(records) - len(unmapped), len(windows))
    cohort_window_threads = _parallel_task_threads(args, len(windows))

    def clean_window(item):
        window_index, window = item
        local_args = argparse.Namespace(**vars(args))
        local_args.threads = cohort_window_threads
        representatives: Dict[str, Tuple[str, int, int, str]] = {}
        provenance: Dict[str, set] = {}
        workdir = os.path.join(idir, "local_groups", f"window{window_index:06d}")
        local_fa = os.path.join(workdir, "local_reference.fa")
        nearby, coordinate_map = write_local_reference(main_fa, window, local_fa)
        log.info("Cohort window %s:%d-%d: %d candidates, %d promoted loci",
                 window.main, window.start, window.end, len(window.pieces), len(nearby))
        rep_kmers = getattr(local_args, "winnow_rep_kmers", None)
        survivors = exclude_against_reference_converge(
            window.pieces, local_fa, os.path.join(workdir, "to_local"), local_args, log,
            representatives, blast_mode="independent", broad_aligner="winnowmap",
            winnow_rep_kmers=rep_kmers, target_coordinate_map=coordinate_map,
            max_cycles=1, warn_on_limit=False,
        )
        self_fn = make_self_align_fn(
            os.path.join(workdir, "selfclean"), local_args, log,
            blast_mode="independent", broad_aligner="winnowmap",
        )
        survivors = self_clean_converge(
            survivors, self_fn, local_args.min_segment, "insertion", log,
            max_cycles=1, provenance=provenance, warn_on_limit=False,
            sample_priorities=getattr(local_args, "sample_priorities", None),
        )
        return window_index, survivors, representatives, provenance

    survivors: List[Piece] = []
    representatives: Dict[str, Tuple[str, int, int, str]] = {}
    if windows:
        workers = min(_local_jobs(args), len(windows))
        log.info("Cleaning cohort windows with %d parallel jobs x %d threads",
                 workers, cohort_window_threads)
        with ThreadPoolExecutor(max_workers=workers, thread_name_prefix="cohort-window") as executor:
            results = list(executor.map(clean_window, enumerate(windows)))
        for _index, kept, local_reps, _provenance in results:
            survivors.extend(kept)
            representatives.update(local_reps)

    result_by_lift = {piece.lift_id: result for piece, result in records}
    final_promoted: List[Piece] = []
    final_small: List[Piece] = []
    survivor_records: List[Tuple[Piece, LocalizeResult]] = []
    result_rows: List[str] = []
    checkpoint_records: List[Tuple[str, Piece, LocalizeResult]] = []
    for piece in survivors:
        result = result_by_lift[piece.lift_id]
        if qualifies_mapped_promotion(
            piece, args.promote_size, args.min_segment, args.min_unmasked_fraction,
        ):
            assign_final_ids([piece], f"mainins_{piece.sample}", counters)
            final_promoted.append(piece)
            pool = "main"
        else:
            assign_final_ids([piece], f"ins_{piece.sample}", counters)
            final_small.append(piece)
            pool = "small"
        survivor_records.append((piece, result))
        checkpoint_records.append((pool, piece, result))
        result_rows.append(table_row(piece, pool, result))

    kept_origins = {piece.origin_id for piece in survivors}
    written_redundant = set()
    for piece, _result in records:
        if piece.origin_id in kept_origins or piece.origin_id in written_redundant:
            continue
        written_redundant.add(piece.origin_id)
        result_rows.append(redundant_row(piece, representatives.get(piece.origin_id)))

    write_discovery_bundle(
        idir, survivor_records, final_promoted,
        checkpoint_records=checkpoint_records,
    )
    _write_stage_results(idir, result_rows)
    compact_discovery_iteration(idir, args.keep_work)
    _write_stage_marker(idir, "_READY", marker_payload)
    _commit_stage(
        idir, checkpoint_records, main_fa, small_fa, table_path, marker_payload,
    )
    log.info("Cohort cleanup complete: %d small, %d final promoted",
             len(final_small), len(final_promoted))
    return final_promoted, final_small


def _atomic_json(path: str, payload: dict) -> None:
    """Write a small checkpoint atomically on the same filesystem."""
    Path(path).parent.mkdir(parents=True, exist_ok=True)
    temporary = path + ".tmp"
    with open(temporary, "w") as out:
        json.dump(payload, out, sort_keys=True)
        out.write("\n")
        out.flush()
        os.fsync(out.fileno())
    os.replace(temporary, path)


def _stage_marker(iteration_dir: str, name: str) -> str:
    return os.path.join(iteration_dir, name)


def _write_stage_marker(iteration_dir: str, name: str, payload: dict) -> None:
    _atomic_json(_stage_marker(iteration_dir, name), payload)


def _read_stage_marker(iteration_dir: str, name: str) -> Optional[dict]:
    path = _stage_marker(iteration_dir, name)
    if not os.path.isfile(path):
        return None
    try:
        with open(path) as handle:
            value = json.load(handle)
    except (OSError, ValueError, TypeError):
        return None
    return value if isinstance(value, dict) else None


def _checkpoint_record(category: str, piece: Piece, result: LocalizeResult) -> dict:
    fields = dataclasses.asdict(piece)
    fields.pop("seq", None)
    return {
        "category": category,
        "piece": fields,
        "result": {
            "status": result.status,
            "placements": [dataclasses.asdict(placement) for placement in result.placements],
            "note": result.note,
        },
    }


def _write_resume_state(
    iteration_dir: str,
    records: Sequence[Tuple[str, Piece, LocalizeResult]],
) -> None:
    path = os.path.join(iteration_dir, "state.jsonl")
    with open(path, "w") as out:
        for category, piece, result in records:
            out.write(json.dumps(_checkpoint_record(category, piece, result), sort_keys=True) + "\n")


def _iter_resume_state(
    iteration_dir: str,
):
    """Stream restart records without retaining every sequence in memory."""
    metadata_path = os.path.join(iteration_dir, "state.jsonl")
    fasta_path = os.path.join(iteration_dir, "insertions.fa")
    if not os.path.isfile(metadata_path) or not os.path.isfile(fasta_path):
        raise RuntimeError(f"Incomplete restart state in {iteration_dir}")
    with IndexedFasta(fasta_path) as sequences, open(metadata_path) as metadata:
        for line_no, raw in enumerate(metadata, 1):
            if not raw.strip():
                continue
            try:
                row = json.loads(raw)
                fields = dict(row["piece"])
                # Compatibility with the short-lived rank-persisting format.
                # Sample priority is now always supplied by query_paths.txt.
                fields.pop("sample_rank", None)
                fields["seq"] = sequences.sequence(fields["temp_id"])
                piece = Piece(**fields)
                result_row = row["result"]
                result = LocalizeResult(
                    result_row["status"],
                    [Placement(**placement) for placement in result_row.get("placements", [])],
                    note=result_row.get("note", ""),
                )
                yield row["category"], piece, result
            except (KeyError, TypeError, ValueError) as exc:
                raise RuntimeError(
                    f"Invalid restart state {metadata_path}:{line_no}: {exc}"
                ) from exc


def _read_resume_state(
    iteration_dir: str,
) -> List[Tuple[str, Piece, LocalizeResult]]:
    """Materialize restart state for stages whose algorithms require a list."""
    return list(_iter_resume_state(iteration_dir))


def _candidate_anchor_text(result: LocalizeResult) -> str:
    return ";".join(
        f"{placement.main}:{placement.start}-{placement.end}:{placement.strand}"
        for placement in result.placements
    ) or "."


def _write_candidate_bundle(
    fasta_path: str,
    table_path: str,
    records,
) -> int:
    """Write restartable candidate sequence and metadata files atomically."""
    Path(fasta_path).parent.mkdir(parents=True, exist_ok=True)
    Path(table_path).parent.mkdir(parents=True, exist_ok=True)
    fasta_tmp = fasta_path + ".tmp"
    table_tmp = table_path + ".tmp"
    record_count = 0
    try:
        with open(fasta_tmp, "w") as fasta, open(table_tmp, "w", newline="") as table:
            writer = csv.DictWriter(
                table, fieldnames=CANDIDATE_TABLE_FIELDS,
                delimiter="\t", lineterminator="\n",
            )
            writer.writeheader()
            for piece, result in records:
                record_count += 1
                anchors = [dataclasses.asdict(placement) for placement in result.placements]
                anchor_text = _candidate_anchor_text(result)
                # Only the first token is interpreted as the FASTA identifier.
                # Source/anchor annotations remain attached for inspection and
                # survive the ordered merge into the cohort candidate FASTA.
                fasta.write(
                    f">{piece.temp_id} sample={sanitize_id(piece.sample)} "
                    f"source={sanitize_id(piece.contig)}:{piece.start}-{piece.end}{piece.strand} "
                    f"anchor_status={result.status} anchors={anchor_text}\n"
                    f"{wrap_fasta(piece.seq)}\n"
                )
                writer.writerow({
                    "candidate_id": piece.temp_id,
                    "sample": piece.sample,
                    "contig": piece.contig,
                    "start": piece.start,
                    "end": piece.end,
                    "strand": piece.strand,
                    "length": piece.length,
                    "unmasked": piece.unmasked,
                    "final_id": piece.final_id,
                    "parent_id": piece.parent_id,
                    "note": piece.note,
                    "origin_id": piece.origin_id,
                    "lift_id": piece.lift_id,
                    "anchor_status": result.status,
                    "anchors_json": json.dumps(anchors, separators=(",", ":")),
                    "anchor_note": result.note,
                })
        os.replace(fasta_tmp, fasta_path)
        os.replace(table_tmp, table_path)
    finally:
        for temporary in (fasta_tmp, table_tmp):
            try:
                os.remove(temporary)
            except FileNotFoundError:
                pass
    return record_count


def _read_candidate_bundle(
    fasta_path: str,
    table_path: str,
) -> List[Tuple[Piece, LocalizeResult]]:
    """Restore candidate objects from the merged disk bundle."""
    records: List[Tuple[Piece, LocalizeResult]] = []
    within_sample_counts: Dict[str, int] = {}
    with IndexedFasta(fasta_path) as sequences, open(table_path, newline="") as table:
        reader = csv.DictReader(table, delimiter="\t")
        missing = set(CANDIDATE_TABLE_FIELDS) - set(reader.fieldnames or ())
        if missing:
            raise RuntimeError(
                f"Candidate table {table_path} is missing columns: {', '.join(sorted(missing))}"
            )
        available = set(sequences.names())
        for line_number, row in enumerate(reader, 2):
            candidate_id = row["candidate_id"]
            if candidate_id not in available:
                raise RuntimeError(
                    f"Candidate FASTA {fasta_path} lacks {candidate_id!r} from "
                    f"{table_path}:{line_number}"
                )
            try:
                anchors = json.loads(row["anchors_json"] or "[]")
                sample = row["sample"]
                within_sample_rank = within_sample_counts.get(sample, 0)
                within_sample_counts[sample] = within_sample_rank + 1
                piece = Piece(
                    temp_id=candidate_id,
                    sample=sample,
                    contig=row["contig"],
                    start=int(row["start"]),
                    end=int(row["end"]),
                    strand=row["strand"],
                    seq=sequences.sequence(candidate_id),
                    final_id=row["final_id"],
                    parent_id=row["parent_id"],
                    note=row["note"],
                    origin_id=row["origin_id"],
                    lift_id=row["lift_id"],
                    within_sample_rank=within_sample_rank,
                )
                result = LocalizeResult(
                    row["anchor_status"],
                    [Placement(**placement) for placement in anchors],
                    note=row["anchor_note"],
                )
            except (TypeError, ValueError, KeyError, json.JSONDecodeError) as exc:
                raise RuntimeError(
                    f"Invalid candidate metadata {table_path}:{line_number}: {exc}"
                ) from exc
            if piece.length != int(row["length"]):
                raise RuntimeError(
                    f"Candidate length mismatch for {candidate_id!r} in {table_path}:{line_number}"
                )
            records.append((piece, result))
    return records


def _merge_candidate_bundles(
    shard_paths: Sequence[Tuple[str, str]],
    fasta_path: str,
    table_path: str,
    log: Optional[logging.Logger] = None,
) -> None:
    """Concatenate ordered per-assembly shards without materializing them."""
    Path(fasta_path).parent.mkdir(parents=True, exist_ok=True)
    fasta_tmp = fasta_path + ".tmp"
    table_tmp = table_path + ".tmp"
    try:
        with open(fasta_tmp, "w") as fasta_out, open(table_tmp, "w", newline="") as table_out:
            table_out.write("\t".join(CANDIDATE_TABLE_FIELDS) + "\n")
            total = len(shard_paths)
            progress_step = max(1, min(100, total // 20)) if total else 1
            for index, (shard_fasta, shard_table) in enumerate(shard_paths, 1):
                with open(shard_fasta) as fasta_in:
                    shutil.copyfileobj(fasta_in, fasta_out, length=1024 * 1024)
                with open(shard_table) as table_in:
                    next(table_in, None)
                    shutil.copyfileobj(table_in, table_out, length=1024 * 1024)
                if log is not None and (index == total or index % progress_step == 0):
                    log.info("Candidate merge progress: %d/%d shards", index, total)
        os.replace(fasta_tmp, fasta_path)
        os.replace(table_tmp, table_path)
    finally:
        for temporary in (fasta_tmp, table_tmp):
            try:
                os.remove(temporary)
            except FileNotFoundError:
                pass


def _write_stage_results(iteration_dir: str, rows: Sequence[str]) -> None:
    with open(os.path.join(iteration_dir, "results.tsv"), "w") as out:
        for row in rows:
            out.write(row if row.endswith("\n") else row + "\n")


def _append_table_results_unique(table_path: str, results_path: str) -> None:
    """Append a committed stage result once, keyed by the first TSV field."""
    with _TABLE_APPEND_LOCK:
        key = os.path.realpath(table_path)
        stat = os.stat(table_path) if os.path.isfile(table_path) else None
        cached = _TABLE_ID_CACHE.get(key)
        if (stat is not None and cached is not None
                and cached[:2] == (stat.st_size, stat.st_mtime_ns)):
            existing = cached[2]
        else:
            existing = set()
            if os.path.isfile(table_path):
                with open(table_path) as table:
                    next(table, None)
                    for raw in table:
                        if raw.strip():
                            existing.add(raw.split("\t", 1)[0])
        with open(table_path, "a") as table, open(results_path) as results:
            for raw in results:
                if not raw.strip():
                    continue
                record_id = raw.split("\t", 1)[0]
                if record_id in existing:
                    continue
                table.write(raw if raw.endswith("\n") else raw + "\n")
                existing.add(record_id)
        stat = os.stat(table_path)
        _TABLE_ID_CACHE[key] = (stat.st_size, stat.st_mtime_ns, existing)


_TABLE_APPEND_LOCK = threading.Lock()
_TABLE_ID_CACHE: Dict[str, Tuple[int, int, set]] = {}
_FASTA_ID_CACHE: Dict[str, Tuple[int, int, set]] = {}


def _known_fasta_ids(fa_path: str) -> set:
    """Read FASTA ids once per process, then track append-only growth in memory."""
    key = os.path.realpath(fa_path)
    if not os.path.isfile(fa_path):
        _FASTA_ID_CACHE[key] = (0, 0, set())
        return _FASTA_ID_CACHE[key][2]
    stat = os.stat(fa_path)
    cached = _FASTA_ID_CACHE.get(key)
    if cached is not None and cached[:2] == (stat.st_size, stat.st_mtime_ns):
        return cached[2]
    identifiers = set(read_primary_fasta_names(fa_path))
    _FASTA_ID_CACHE[key] = (stat.st_size, stat.st_mtime_ns, identifiers)
    return identifiers


def _append_pieces(fa_path: str, pieces: Sequence[Piece]) -> None:
    """Append FASTA records idempotently by their final/temporary identifier."""
    Path(fa_path).parent.mkdir(parents=True, exist_ok=True)
    existing = _known_fasta_ids(fa_path)
    with open(fa_path, "a") as out:
        for p in pieces:
            primary = p.final_id or p.temp_id
            if primary in existing:
                continue
            out.write(f">{primary} {p.source_coord()}\n{wrap_fasta(p.seq)}\n")
            existing.add(primary)
    stat = os.stat(fa_path)
    _FASTA_ID_CACHE[os.path.realpath(fa_path)] = (stat.st_size, stat.st_mtime_ns, existing)
    # A pre-existing samtools index becomes stale as soon as promoted records
    # are appended. IndexedFasta will safely rebuild its in-memory index later.
    try:
        os.remove(fa_path + ".fai")
    except FileNotFoundError:
        pass


def _append_annotated_novel_records(
    fa_path: str,
    records: Sequence[Tuple[Piece, LocalizeResult]],
) -> None:
    """Append final novel cores idempotently with lift annotations."""
    Path(fa_path).parent.mkdir(parents=True, exist_ok=True)
    existing = _known_fasta_ids(fa_path)
    with open(fa_path, "a") as out:
        for piece, result in records:
            primary = piece.final_id or piece.temp_id
            if primary in existing:
                continue
            out.write(f"{_lift_header(piece, result, True)}\n{wrap_fasta(piece.seq)}\n")
            existing.add(primary)
    stat = os.stat(fa_path)
    _FASTA_ID_CACHE[os.path.realpath(fa_path)] = (
        stat.st_size, stat.st_mtime_ns, existing,
    )
    try:
        os.remove(fa_path + ".fai")
    except FileNotFoundError:
        pass


def write_discovery_bundle(
    iteration_dir: str, records: Sequence[Tuple[Piece, LocalizeResult]],
    promoted: Sequence[Piece],
    checkpoint_records: Sequence[Tuple[str, Piece, LocalizeResult]] = (),
) -> None:
    """Persist the candidate sequences and placements needed by cohort cleanup."""
    Path(iteration_dir).mkdir(parents=True, exist_ok=True)
    # The FASTA is also the sequence store for restart metadata.  Include exact
    # post-cleaning checkpoint pieces in addition to the complete discovery
    # catalog, without duplicating a record that already has the same temp id.
    sequence_by_id: Dict[str, Piece] = {}
    for piece, _result in records:
        sequence_by_id.setdefault(piece.temp_id, piece)
    for _category, piece, _result in checkpoint_records:
        previous = sequence_by_id.get(piece.temp_id)
        if previous is not None and previous.seq != piece.seq:
            raise RuntimeError(f"Restart FASTA id collision: {piece.temp_id}")
        sequence_by_id[piece.temp_id] = piece
    write_fasta_pieces(list(sequence_by_id.values()),
                       os.path.join(iteration_dir, "insertions.fa"), use_final_id=False)
    write_fasta_pieces(promoted, os.path.join(iteration_dir, "promoted.fa"), use_final_id=True)
    with open(os.path.join(iteration_dir, "placements.tsv"), "w") as out:
        out.write("candidate_id\tsample\tcontig\tstart\tend\tlength\tunmasked\tstatus\tplacements\n")
        for piece, result in records:
            placements = ";".join(
                f"{p.main}:{p.start}-{p.end}:{p.strand}" for p in result.placements
            ) or "."
            out.write("\t".join([
                piece.final_id or piece.temp_id, piece.sample, piece.contig,
                str(piece.start), str(piece.end), str(piece.length), str(piece.unmasked),
                result.status, placements,
            ]) + "\n")
    _write_resume_state(iteration_dir, checkpoint_records)


def compact_discovery_iteration(iteration_dir: str, keep_work: bool = False) -> None:
    """Keep raw PAFs plus compact insertion/promotion records; remove other temp files."""
    if keep_work:
        return
    root_keep = {
        "insertions.fa", "promoted.fa", "placements.tsv", "state.jsonl",
        "results.tsv", LIFT_BLOCKS_NAME, "remerged.fa",
        "remerged.with_10kb_anchors.fa", "remerged.lifted.fa",
        "full_reference_cleaned.fa",
        "_READY", "_SUCCESS",
    }
    for root, dirs, files in os.walk(iteration_dir, topdown=False):
        for filename in files:
            path = os.path.join(root, filename)
            at_root = os.path.realpath(root) == os.path.realpath(iteration_dir)
            if (filename.endswith(".paf") or filename.endswith(".map.tsv")
                    or (at_root and filename in root_keep)):
                continue
            os.remove(path)
        for dirname in dirs:
            path = os.path.join(root, dirname)
            try:
                os.rmdir(path)
            except OSError:
                pass


def compact_reference_iteration(
    iteration_dir: str,
    keep_work: bool = False,
) -> Tuple[int, int]:
    """Compact completed first-reference work and report removed FASTA storage."""
    if keep_work or not os.path.isdir(iteration_dir):
        return 0, 0
    fasta_paths: List[Tuple[str, int]] = []
    for root, _dirs, files in os.walk(iteration_dir):
        for filename in files:
            lower = filename.lower()
            if lower.endswith((".fa", ".fasta", ".fna")):
                path = os.path.join(root, filename)
                try:
                    fasta_paths.append((path, os.path.getsize(path)))
                except FileNotFoundError:
                    pass
    compact_discovery_iteration(iteration_dir, keep_work=False)
    # Discovery compaction preserves a few root-level result FASTAs, but every
    # FASTA below the completed reference-clean work directory is intermediate.
    for path, _size in fasta_paths:
        for candidate in (path, path + ".fai"):
            try:
                os.remove(candidate)
            except FileNotFoundError:
                pass
    removed = [(path, size) for path, size in fasta_paths if not os.path.exists(path)]
    return len(removed), sum(size for _path, size in removed)


def _checkpoint_matches(
    marker: dict, assembly: AssemblyEntry, iteration: int,
    config_signature: Optional[str] = None,
    kind: str = "assembly",
    require_iteration: bool = True,
) -> bool:
    return (
        marker.get("kind") == kind
        and (not require_iteration or marker.get("iteration") == iteration)
        and marker.get("sample") == assembly.sample
        and os.path.realpath(str(marker.get("path", ""))) == os.path.realpath(assembly.path)
        and marker.get("size") == os.path.getsize(assembly.path)
        and marker.get("mtime_ns") == os.stat(assembly.path).st_mtime_ns
        and (config_signature is None
             or marker.get("config_signature") == config_signature)
    )


def _commit_stage(
    iteration_dir: str,
    state: Sequence[Tuple[str, Piece, LocalizeResult]],
    main_fa: str,
    small_fa: str,
    table_path: str,
    marker_payload: dict,
) -> None:
    """Idempotently install a READY stage into the cumulative outputs."""
    main_pieces = [piece for category, piece, _ in state
                   if category in {"promoted", "main"}]
    small_pieces = [piece for category, piece, _ in state if category == "small"]
    if main_pieces:
        _append_pieces(main_fa, main_pieces)
    if small_pieces:
        _append_pieces(small_fa, small_pieces)
    results_path = os.path.join(iteration_dir, "results.tsv")
    if not os.path.isfile(results_path):
        raise RuntimeError(f"Committed stage is missing results.tsv: {iteration_dir}")
    _append_table_results_unique(table_path, results_path)
    _write_stage_marker(iteration_dir, "_SUCCESS", marker_payload)


def _load_or_commit_finished_assembly(
    iteration_dir: str,
    assembly: AssemblyEntry,
    iteration: int,
    main_fa: str,
    small_fa: str,
    table_path: str,
    log: logging.Logger,
    config_signature: Optional[str] = None,
) -> Optional[Tuple[List[Tuple[Piece, LocalizeResult]],
                    List[Tuple[Piece, LocalizeResult]]]]:
    """Return checkpoint records when an assembly is READY or SUCCESS."""
    success = _read_stage_marker(iteration_dir, "_SUCCESS")
    ready = _read_stage_marker(iteration_dir, "_READY")
    marker = success or ready
    if marker is None:
        return None
    if not _checkpoint_matches(marker, assembly, iteration, config_signature):
        raise RuntimeError(
            f"Restart checkpoint does not match assembly {iteration} "
            f"({assembly.sample} {assembly.path}): {iteration_dir}"
        )
    state = _read_resume_state(iteration_dir)
    if success is None:
        log.info("Completing READY checkpoint for assembly %s", assembly.sample)
    # Recommit SUCCESS as well. Cumulative FASTA/table outputs are rebuilt on
    # every resume from these immutable stage states, so SUCCESS means the stage
    # computation is complete, not that a mutable aggregate file must exist.
    _commit_stage(iteration_dir, state, main_fa, small_fa, table_path, marker)
    promoted = [(piece, result) for category, piece, result in state
                if category == "promoted"]
    cohort = [(piece, result) for category, piece, result in state
              if category == "cohort"]
    log.info("RESUME: skipped finished assembly %s (%d promoted, %d cohort candidates)",
             assembly.sample, len(promoted), len(cohort))
    return promoted, cohort


def _load_or_commit_assembly_stage(
    iteration_dir: str,
    assembly: AssemblyEntry,
    iteration: int,
    kind: str,
    main_fa: str,
    small_fa: str,
    table_path: str,
    log: logging.Logger,
    config_signature: str,
) -> Optional[List[Tuple[str, Piece, LocalizeResult]]]:
    """Load a generic assembly-stage checkpoint, committing READY if needed."""
    success = _read_stage_marker(iteration_dir, "_SUCCESS")
    ready = _read_stage_marker(iteration_dir, "_READY")
    marker = success or ready
    if marker is None:
        return None
    expected_logic = STAGE_LOGIC_VERSIONS.get(kind, 1)
    observed_logic = int(marker.get("logic_version", 1))
    if observed_logic != expected_logic:
        log.info(
            "RESUME: invalidating %s checkpoint for %s (logic version %d -> %d)",
            kind, assembly.sample, observed_logic, expected_logic,
        )
        return None
    if not _checkpoint_matches(
        marker, assembly, iteration, None, kind=kind,
        require_iteration=(kind != "main_discovery"),
    ):
        raise RuntimeError(
            f"Restart checkpoint does not match {kind} assembly {iteration} "
            f"({assembly.sample} {assembly.path}): {iteration_dir}"
        )
    if marker.get("config_signature") != config_signature:
        log.info(
            "RESUME: invalidating %s checkpoint for %s because configuration changed",
            kind, assembly.sample,
        )
        return None
    state = _read_resume_state(iteration_dir)
    if success is None:
        log.info("Completing READY %s checkpoint for %s", kind, assembly.sample)
    _commit_stage(iteration_dir, state, main_fa, small_fa, table_path, marker)
    log.info("CHECKPOINT: loaded finished %s for %s", kind, assembly.sample)
    return state


def _load_main_discovery_checkpoint(
    iteration_dir: str,
    assembly: AssemblyEntry,
    iteration: int,
    main_fa: str,
    small_fa: str,
    table_path: str,
    log: logging.Logger,
    config_signature: str,
) -> bool:
    """Validate/commit a finished main stage without reading intermediate data."""
    success = _read_stage_marker(iteration_dir, "_SUCCESS")
    ready = _read_stage_marker(iteration_dir, "_READY")
    marker = success or ready
    if marker is None:
        return False
    expected_logic = STAGE_LOGIC_VERSIONS["main_discovery"]
    observed_logic = int(marker.get("logic_version", 1))
    if observed_logic != expected_logic:
        log.info(
            "RESUME: invalidating main_discovery checkpoint for %s "
            "(logic version %d -> %d)",
            assembly.sample, observed_logic, expected_logic,
        )
        return False
    if not _checkpoint_matches(
        marker, assembly, iteration, config_signature, kind="main_discovery",
        require_iteration=False,
    ):
        raise RuntimeError(
            f"Restart checkpoint does not match main_discovery assembly {iteration} "
            f"({assembly.sample} {assembly.path}): {iteration_dir}"
        )
    if success is None:
        log.info("Completing READY main_discovery checkpoint for %s", assembly.sample)
    # Main discovery has insertion-only state and therefore never appends a
    # cumulative FASTA. Its results.tsv is still recommitted idempotently.
    _commit_stage(iteration_dir, [], main_fa, small_fa, table_path, marker)
    log.info("CHECKPOINT: finished main_discovery for %s (marker only)", assembly.sample)
    return True


_COUNTER_RE = re.compile(r"^(.+)_([0-9]{6})$")


def _update_counters_from_pieces(counters: Dict[str, int], pieces: Sequence[Piece]) -> None:
    for piece in pieces:
        match = _COUNTER_RE.match(piece.final_id or piece.temp_id)
        if not match:
            continue
        prefix, number = match.group(1), int(match.group(2))
        counters[prefix] = max(counters.get(prefix, 0), number)


def _seed_counters_from_outputs(*fasta_paths: str) -> Dict[str, int]:
    counters: Dict[str, int] = {}
    for fasta_path in fasta_paths:
        if not os.path.isfile(fasta_path):
            continue
        for name in read_primary_fasta_names(fasta_path):
            match = _COUNTER_RE.match(name)
            if not match:
                continue
            prefix, number = match.group(1), int(match.group(2))
            counters[prefix] = max(counters.get(prefix, 0), number)
    return counters


def _looks_like_oom(exc: BaseException) -> bool:
    if isinstance(exc, subprocess.CalledProcessError):
        return exc.returncode in {-9, 9, 137}
    message = str(exc).lower()
    return "out of memory" in message or "oom" in message or "killed" in message


def _physical_memory_bytes() -> Optional[int]:
    try:
        pages = os.sysconf("SC_PHYS_PAGES")
        page_size = os.sysconf("SC_PAGE_SIZE")
        if pages > 0 and page_size > 0:
            return int(pages * page_size)
    except (AttributeError, OSError, ValueError):
        pass
    return None


def _read_memory_limit(path: str) -> Optional[int]:
    try:
        with open(path) as handle:
            text = handle.read().strip()
    except OSError:
        return None
    if not text or text.lower() == "max":
        return None
    try:
        value = int(text)
    except ValueError:
        return None
    # cgroup v1 uses values near LONG_MAX to mean unlimited.
    return value if 0 < value < (1 << 60) else None


def _cgroup_memory_limit_bytes() -> Optional[int]:
    """Return this process's cgroup memory ceiling on cgroup v1 or v2."""
    paths = {
        "/sys/fs/cgroup/memory.max",
        "/sys/fs/cgroup/memory/memory.limit_in_bytes",
    }
    try:
        with open("/proc/self/cgroup") as handle:
            for raw in handle:
                fields = raw.rstrip("\n").split(":", 2)
                if len(fields) != 3:
                    continue
                controllers, relative = fields[1], fields[2].lstrip("/")
                if not controllers:  # unified cgroup v2
                    paths.add(os.path.join("/sys/fs/cgroup", relative, "memory.max"))
                elif "memory" in controllers.split(","):  # cgroup v1
                    paths.add(os.path.join(
                        "/sys/fs/cgroup/memory", relative, "memory.limit_in_bytes",
                    ))
                    paths.add(os.path.join(
                        "/sys/fs/cgroup", relative, "memory.limit_in_bytes",
                    ))
    except OSError:
        pass
    limits = [value for path in paths if (value := _read_memory_limit(path)) is not None]
    return min(limits) if limits else None


def _slurm_memory_limit_bytes() -> Optional[int]:
    """Read Slurm's per-node allocation (documented in MiB)."""
    raw = os.environ.get("SLURM_MEM_PER_NODE", "").strip()
    if not raw:
        return None
    match = re.fullmatch(r"([0-9]+(?:\.[0-9]+)?)\s*([KMGT]?)B?", raw, re.I)
    if match is None:
        return None
    number = float(match.group(1))
    suffix = match.group(2).upper()
    # A suffix-free SLURM_MEM_PER_NODE value is always expressed in MiB.
    powers = {"": 2, "K": 1, "M": 2, "G": 3, "T": 4}
    return int(number * (1024 ** powers[suffix]))


def _visible_memory_info() -> Tuple[Optional[float], str]:
    candidates = [
        ("physical RAM", _physical_memory_bytes()),
        ("cgroup limit", _cgroup_memory_limit_bytes()),
        ("Slurm allocation", _slurm_memory_limit_bytes()),
    ]
    available = [(label, value) for label, value in candidates if value is not None]
    if not available:
        return None, "unknown"
    label, value = min(available, key=lambda item: item[1])
    return value / (1024.0 ** 3), label


def _visible_memory_gb() -> Optional[float]:
    """Compatibility wrapper returning the effective allocation, not host RAM."""
    return _visible_memory_info()[0]


def run_parallel_with_retries(
    items: Sequence[Tuple[int, AssemblyEntry]],
    worker: Callable[[Tuple[int, AssemblyEntry], int], object],
    jobs: int,
    retries: int,
    label: str,
    log: logging.Logger,
) -> Dict[int, object]:
    """Run assembly jobs concurrently and retry failures in smaller waves."""
    pending = list(items)
    results: Dict[int, object] = {}
    concurrency = max(1, min(jobs, len(pending))) if pending else 1
    for attempt in range(retries + 1):
        if not pending:
            break
        log.info("%s: attempt %d/%d with %d concurrent jobs",
                 label, attempt + 1, retries + 1, concurrency)
        failed: List[Tuple[int, AssemblyEntry]] = []
        failures: List[Tuple[Tuple[int, AssemblyEntry], BaseException]] = []
        with ThreadPoolExecutor(max_workers=concurrency,
                                thread_name_prefix=sanitize_id(label)) as executor:
            future_to_item = {
                executor.submit(worker, item, attempt): item for item in pending
            }
            completed = 0
            total_attempt = len(future_to_item)
            progress_step = max(1, min(100, total_attempt // 20))
            for future in as_completed(future_to_item):
                item = future_to_item[future]
                iteration, assembly = item
                try:
                    results[iteration] = future.result()
                except Exception as exc:  # retain every successful sibling checkpoint
                    failed.append(item)
                    failures.append((item, exc))
                    log.error("%s failed for %s (attempt %d): %s",
                              label, assembly.sample, attempt + 1, exc)
                finally:
                    completed += 1
                    if completed == total_attempt or completed % progress_step == 0:
                        log.info(
                            "%s progress: %d/%d jobs finished on attempt %d",
                            label, completed, total_attempt, attempt + 1,
                        )
        if not failed:
            break
        if attempt >= retries:
            details = "; ".join(
                f"{assembly.sample}: {exc}" for (_iteration, assembly), exc in failures
            )
            raise RuntimeError(
                f"{label}: {len(failed)} job(s) failed after {retries + 1} attempts: {details}"
            )
        if any(_looks_like_oom(exc) for _item, exc in failures):
            concurrency = max(1, concurrency // 2)
            log.warning("%s: possible OOM; reducing retry concurrency to %d",
                        label, concurrency)
        else:
            concurrency = max(1, min(concurrency, len(failed)))
        pending = sorted(failed, key=lambda item: item[0])
    return results


def run_group_jobs_with_retries(
    items: Sequence[Tuple[int, object]],
    worker: Callable[[Tuple[int, object], int], object],
    jobs: int,
    retries: int,
    label: str,
    log: logging.Logger,
) -> Dict[int, object]:
    """Retryable counterpart for independent lift-local cleaning groups."""
    pending = list(items)
    results: Dict[int, object] = {}
    concurrency = max(1, min(jobs, len(pending))) if pending else 1
    for attempt in range(retries + 1):
        if not pending:
            break
        failed: List[Tuple[int, object]] = []
        failures: List[Tuple[Tuple[int, object], BaseException]] = []
        log.info("%s: attempt %d/%d with %d parallel groups",
                 label, attempt + 1, retries + 1, concurrency)
        with ThreadPoolExecutor(max_workers=concurrency,
                                thread_name_prefix=sanitize_id(label)) as executor:
            future_to_item = {executor.submit(worker, item, attempt): item for item in pending}
            for future in as_completed(future_to_item):
                item = future_to_item[future]
                index = item[0]
                try:
                    results[index] = future.result()
                except Exception as exc:
                    failed.append(item)
                    failures.append((item, exc))
                    log.error("%s group %d failed (attempt %d): %s",
                              label, index, attempt + 1, exc)
        if not failed:
            break
        if attempt >= retries:
            details = "; ".join(f"group {item[0]}: {exc}" for item, exc in failures)
            raise RuntimeError(
                f"{label}: {len(failed)} group(s) failed after {retries + 1} attempts: {details}"
            )
        if any(_looks_like_oom(exc) for _item, exc in failures):
            concurrency = max(1, concurrency // 2)
            log.warning("%s: possible OOM; reducing retry concurrency to %d",
                        label, concurrency)
        else:
            concurrency = max(1, min(concurrency, len(failed)))
        pending = sorted(failed, key=lambda item: item[0])
    return results


# ---------------------------------------------------------------------------
# Cohort-v3 workflow: independent discovery -> novel loci -> local insertions
# ---------------------------------------------------------------------------

def _sample_stage_dir(outdir: str, kind: str, sample: str) -> str:
    """Stable per-haplotype checkpoint path, independent of query-list order."""
    return os.path.join(outdir, "iterations", kind, sanitize_id(sample))


def _set_current_assembly_rank(
    records: Sequence[Tuple[Piece, LocalizeResult]], iteration: int,
) -> List[Tuple[Piece, LocalizeResult]]:
    """Apply the current query-list priority to records loaded from a checkpoint."""
    return [
        (_piece_with_current_assembly_rank(piece, iteration), result)
        for piece, result in records
    ]


def _piece_with_current_assembly_rank(piece: Piece, iteration: int) -> Piece:
    prefix = f"mainlift_{iteration:04d}_"
    lift_id = re.sub(r"^mainlift_[0-9]{4}_", prefix, piece.lift_id)
    return dataclasses.replace(piece, lift_id=lift_id)

def _assembly_stage_marker(kind: str, assembly: AssemblyEntry, iteration: int,
                           args) -> dict:
    marker = {
        "kind": kind,
        "sample": assembly.sample,
        "path": os.path.realpath(assembly.path),
        "size": os.path.getsize(assembly.path),
        "mtime_ns": os.stat(assembly.path).st_mtime_ns,
        "config_signature": (
            _downstream_signature(args)
            if kind == "novel_alignment" else _restart_signature(args)
        ),
        "logic_version": STAGE_LOGIC_VERSIONS.get(kind, 1),
    }
    # Main discovery is haplotype-addressed and independent of list position.
    # Other legacy/per-order stages still use iteration in their directory key.
    if kind != "main_discovery":
        marker["iteration"] = iteration
    return marker


def _finish_assembly_stage(
    iteration_dir: str,
    marker: dict,
    all_records: Sequence[Tuple[Piece, LocalizeResult]],
    checkpoint_records: Sequence[Tuple[str, Piece, LocalizeResult]],
    result_rows: Sequence[str],
    main_fa: str,
    small_fa: str,
    table_path: str,
    args,
) -> None:
    write_discovery_bundle(
        iteration_dir, all_records, [], checkpoint_records=checkpoint_records,
    )
    _write_stage_results(iteration_dir, result_rows)
    compact_discovery_iteration(iteration_dir, args.keep_work)
    _write_stage_marker(iteration_dir, "_READY", marker)
    _commit_stage(
        iteration_dir, checkpoint_records, main_fa, small_fa, table_path, marker,
    )


def discover_main_novelty(
    assembly: AssemblyEntry,
    iteration: int,
    main_target: str,
    main_fa: str,
    args,
    outdir: str,
    table_path: str,
    log: logging.Logger,
    commit: bool = True,
) -> List[Tuple[Piece, LocalizeResult]]:
    """Independently discover every residual novel segment against fixed main."""
    idir = _sample_stage_dir(outdir, "main_discovery", assembly.sample)
    Path(idir).mkdir(parents=True, exist_ok=True)
    marker = _assembly_stage_marker("main_discovery", assembly, iteration, args)

    main_hits = collect_alignment_hits(
        assembly.path, main_target, os.path.join(idir, "to_main"), args.threads,
        args.min_identity, args.min_segment, log, winnow_params=WINNOW_PARAMS,
        blast_word_size=args.blast_word_size, blast_evalue=args.blast_evalue,
        blast_max_target_seqs=args.blast_max_target_seqs, prefix="to_main",
        skip_blastn=True, blast_mode="residual", broad_aligner="minimap2",
        minimap_params=MINIMAP_PARAMS,
    )
    candidates = extract_insertions(
        assembly, main_hits, args.min_segment, log,
        getattr(args, "blacklist_intervals", None),
    )
    if args.cycle_minimap2 and candidates:
        candidates = exclude_against_reference_converge(
            candidates, main_target, os.path.join(idir, "to_main_minimap2_cycles"),
            args, log, {}, blast_mode="residual", alignment_sink=main_hits,
            broad_aligner="minimap2", run_blastn=False,
            max_cycles=args.cycle_minimap2,
        )
    if getattr(args, "full_main_discovery_cleaning", False) and candidates:
        candidates = exclude_against_reference_converge(
            candidates,
            main_fa,
            os.path.join(idir, "to_main_winnowmap_blastn_cycles"),
            args,
            log,
            {},
            blast_mode="independent",
            alignment_sink=main_hits,
            broad_aligner="winnowmap",
            run_blastn=True,
            winnow_rep_kmers=getattr(args, "winnow_rep_kmers", None),
            warn_on_limit=False,
            unlimited=True,
        )

    # Persist accepted main-alignment blocks near the final minimap2 residuals.
    # They localize original occurrences in later novel/local alignment phases;
    # finalized remerged loci receive a separate fresh 10-kb-anchor lift.
    lift_blocks = write_lift_blocks(
        os.path.join(idir, LIFT_BLOCKS_NAME), candidates, main_hits, args.flank,
    )

    safe_sample = sanitize_id(assembly.sample)
    for index, piece in enumerate(candidates):
        piece.lift_id = f"mainlift_{safe_sample}_{index:08d}"
    lifted = {
        piece.temp_id: result
        for piece, result in zip(
            candidates, localize_from_block_map(candidates, lift_blocks, args.flank),
        )
    }
    discovered, by_lift = split_lift_placements(candidates, lifted)
    records: List[Tuple[Piece, LocalizeResult]] = []
    for index, piece in enumerate(discovered):
        piece.temp_id = f"mainnovel_{safe_sample}_{index:08d}"
        records.append((piece, by_lift[piece.lift_id]))

    checkpoint = [("insertion", piece, result) for piece, result in records]
    if commit:
        _finish_assembly_stage(
            idir, marker, records, checkpoint, [], main_fa,
            os.path.join(outdir, "small_insertions.fa"), table_path, args,
        )
    else:
        # Slurm-array workers own only their haplotype directory. The
        # controller promotes READY to SUCCESS and rebuilds shared aggregate
        # files after the array barrier, avoiding cross-job append races.
        write_discovery_bundle(
            idir, records, [], checkpoint_records=checkpoint,
        )
        _write_stage_results(idir, [])
        compact_discovery_iteration(idir, args.keep_work)
        _write_stage_marker(idir, "_READY", marker)
    log.info("%s main discovery complete: %d novel segments", assembly.sample, len(records))
    return records


def qualifies_local_insertion(piece: Piece, min_unmasked: int,
                              min_fraction: float) -> bool:
    """Filter masked/tiny residuals before lift-local cohort cleaning."""
    return piece.unmasked >= max(float(min_unmasked), min_fraction * piece.length)


def _commit_novel_loci_stage(
    iteration_dir: str,
    state: Sequence[Tuple[str, Piece, LocalizeResult]],
    novel_fa: str,
    unmapped_novel_fa: str,
    table_path: str,
    marker: dict,
) -> None:
    loci = [
        (piece, result) for category, piece, result in state
        if category == "novel_locus"
    ]
    if loci:
        _append_annotated_novel_records(novel_fa, loci)
    elif not os.path.exists(novel_fa):
        Path(novel_fa).touch()
    unmapped = [
        (piece, result) for category, piece, result in state
        if category == "unmapped_novel"
    ]
    if unmapped:
        _append_annotated_novel_records(unmapped_novel_fa, unmapped)
    elif not os.path.exists(unmapped_novel_fa):
        Path(unmapped_novel_fa).touch()
    results_path = os.path.join(iteration_dir, "results.tsv")
    if not os.path.isfile(results_path):
        raise RuntimeError(f"Committed novel-locus stage is missing results.tsv: {iteration_dir}")
    _append_table_results_unique(table_path, results_path)
    _write_stage_marker(iteration_dir, "_SUCCESS", marker)


def _upgrade_v8_novel_loci_split(
    iteration_dir: str,
    state: Sequence[Tuple[str, Piece, LocalizeResult]],
    marker: dict,
    keep_work: bool,
    keep_unmapped_samples: Sequence[str] = (),
    args=None,
) -> Tuple[List[Tuple[str, Piece, LocalizeResult]], dict]:
    """Reclassify/refilter a saved final state without repeating alignments."""
    records = [
        (piece, result) for category, piece, result in state
        if category in {"novel_locus", "unmapped_novel"}
    ]
    protected = set(keep_unmapped_samples)
    reference = [
        (piece, result) for piece, result in records
        if (
            (
                result.placements
                and (
                    args is None
                    or qualifies_mapped_promotion(
                        piece, args.promote_size, args.min_segment,
                        args.min_unmasked_fraction,
                    )
                )
            )
            or (not result.placements and piece.sample in protected)
        )
    ]
    unmapped = [
        (piece, result) for piece, result in records
        if not result.placements and piece.sample not in protected
    ]
    checkpoint = [
        ("novel_locus", piece, result) for piece, result in reference
    ] + [
        ("unmapped_novel", piece, result) for piece, result in unmapped
    ]
    write_discovery_bundle(
        iteration_dir, reference + unmapped,
        [piece for piece, _result in reference],
        checkpoint_records=checkpoint,
    )
    write_annotated_loci_fasta(
        reference, os.path.join(iteration_dir, "promoted.fa"), use_final_id=True,
    )
    _write_stage_results(
        iteration_dir,
        [table_row(piece, "novel_locus", result) for piece, result in reference]
        + [table_row(piece, "unmapped_novel", result)
           for piece, result in unmapped],
    )
    compact_discovery_iteration(iteration_dir, keep_work)
    upgraded_marker = dict(marker)
    upgraded_marker["logic_version"] = NOVEL_CLEANUP_LOGIC_VERSION
    upgraded_marker["keep_unmapped_samples"] = sorted(protected)
    _write_stage_marker(iteration_dir, "_READY", upgraded_marker)
    return checkpoint, upgraded_marker


def build_novel_loci(
    records: Sequence[Tuple[Piece, LocalizeResult]],
    main_fa: str,
    novel_fa: str,
    args,
    outdir: str,
    counters: Dict[str, int],
    table_path: str,
    log: logging.Logger,
    main_lift_blocks: Optional[Dict[Tuple[int, str], List[AnchorBlock]]] = None,
    reload_main_lift_blocks: bool = False,
    assembly_paths: Optional[Dict[str, str]] = None,
    complete_reference_target: Optional[str] = None,
    complete_reference_fasta: Optional[str] = None,
) -> List[Tuple[Piece, LocalizeResult]]:
    """Select one genome-wide representative per large novel sequence locus.

    Original assembly insertion records are never removed here.  This function
    operates on clones solely to decide which sequence representatives enter
    ``novel_loci.fa``. Unlifted representatives are normally written separately
    to ``unmapped_novel.fa`` and are not alignment targets. Samples named by
    ``--keep-unmapped-samples`` are exempt: their retained unlifted
    representatives remain in ``novel_loci.fa``. Every occurrence is
    reconsidered against the resulting locus FASTA in the next
    assembly-independent phase.
    """
    idir = os.path.join(outdir, "iterations", "novel_loci_cleanup")
    Path(idir).mkdir(parents=True, exist_ok=True)
    marker = {
        "kind": "novel_loci_cleanup",
        "record_count": len(records),
        "config_signature": _downstream_signature(args),
        "assembly_order_signature": getattr(args, "assembly_order_signature", None),
        "logic_version": NOVEL_CLEANUP_LOGIC_VERSION,
        "keep_unmapped_samples": sorted(_keep_unmapped_samples(args)),
    }
    _write_stage_marker(idir, "_IN_PROGRESS", {
        **marker,
        "kind": "novel_loci_cleanup_in_progress",
    })
    # Main-discovery may have cloned one exact sequence for two disagreeing
    # flanks. The global sequence set must contain that sequence only once while
    # retaining both possible placements as attributes.
    sequence_records = consolidate_identical_records(records)
    # Promotion-size and unmasked-fraction filtering belongs at the very end.
    # Cleaning may split one biological interval around a long soft-masked
    # region, and those residuals must remain available for weighted remerging.
    # Before cleaning, only the basic min_segment rule is allowed.
    candidates = [
        dataclasses.replace(piece, final_id="") for piece, _result in sequence_records
        if qualifies_preclean_candidate(piece, args.min_segment)
    ]
    log.info(
        "Novel-locus selection: %d/%d unique sequence records pass the "
        "pre-clean min_segment filter (no promotion-size filter yet)",
        len(candidates), len(sequence_records),
    )

    globally_cleaned: List[Piece] = []
    if candidates:
        globally_cleaned = self_clean_minimap_then_winnow(
            candidates, os.path.join(idir, "global_initial"),
            args, "insertion", log, "Initial novel-locus cleanup",
            unlimited=True,
        )
    reference_cleaned = globally_cleaned
    if globally_cleaned and complete_reference_fasta is not None:
        reference_cleaned = clean_against_complete_reference(
            globally_cleaned,
            (
                complete_reference_target
                if complete_reference_target is not None
                else complete_reference_fasta
            ),
            complete_reference_fasta,
            os.path.join(idir, "full_reference_after_selfclean"),
            args,
            log,
        )
    reference_cleaned_fa = os.path.join(idir, "full_reference_cleaned.fa")
    write_fasta_pieces(
        reference_cleaned, reference_cleaned_fa, use_final_id=False,
    )
    log.info(
        "Wrote %d post-self-clean full-reference residuals to %s",
        len(reference_cleaned), reference_cleaned_fa,
    )
    if assembly_paths is None:
        raise RuntimeError("Source assembly paths are required for weighted remerging and 10-kb lift")

    # Candidate self-cleaning and final full-reference exclusion both finish
    # before source-contig weighted-gap remerging.  Do not size-filter the
    # individual residuals: a long soft-masked bridge is deliberately cheap.
    remerged = merge_weighted_source_gaps(
        reference_cleaned,
        assembly_paths,
        max_weighted_gap=NOVEL_MERGE_WEIGHT_THRESHOLD,
        masked_weight=NOVEL_MERGE_MASKED_WEIGHT,
        log=log,
    )
    remerged_fa = os.path.join(idir, "remerged.fa")
    write_fasta_pieces(remerged, remerged_fa, use_final_id=False)
    log.info("Wrote %d weighted-gap-remerged cores to %s", len(remerged), remerged_fa)

    main_index = os.path.join(outdir, "indexes", "main_chroms.asm5.mmi")
    main_target = (
        complete_reference_target
        if complete_reference_target is not None
        else (main_index if os.path.isfile(main_index) else main_fa)
    )
    lifted_records = lift_reassembled_loci(
        remerged, assembly_paths, main_target,
        os.path.join(idir, "remerged_lift"), args, log,
    )
    lifted_fa = os.path.join(idir, "remerged.lifted.fa")
    write_annotated_loci_fasta(lifted_records, lifted_fa)
    log.info("Wrote freshly lifted reassembled cores to %s", lifted_fa)

    # Freshly lifted cores normally define reference loci. Retained unlifted
    # cores from explicitly protected samples also remain alignment targets.
    # Other unlifted cores are kept separately and never become targets.
    protected_samples = set(_keep_unmapped_samples(args))
    reference_records: List[Tuple[Piece, LocalizeResult]] = []
    unmapped_records: List[Tuple[Piece, LocalizeResult]] = []
    protected_unmapped_count = 0
    final_ids: Set[str] = set()
    failed_mapped_filter = 0
    failed_unmapped_filter = 0
    for piece, result in lifted_records:
        # Cleaning and weighted-gap remerging can materially change both total and
        # unmasked length. Reapply the biological promotion rule only after
        # the fresh lift, so novel_loci.fa never receives a stale pre-cleaning
        # classification. Unlifted records remain separately auditable in
        # unmapped_novel.fa (or protected by --keep-unmapped-samples).
        if not qualifies_promotion(piece, result, args):
            if result.placements:
                failed_mapped_filter += 1
            else:
                failed_unmapped_filter += 1
            continue
        piece.final_id = source_interval_id(
            piece.sample, piece.contig, piece.start, piece.end,
        )
        if piece.final_id in final_ids:
            raise ValueError(
                f"duplicate coordinate-derived novel-locus id: {piece.final_id}"
            )
        final_ids.add(piece.final_id)
        if result.placements or piece.sample in protected_samples:
            reference_records.append((piece, result))
            if not result.placements:
                protected_unmapped_count += 1
        else:
            unmapped_records.append((piece, result))
    log.info(
        "Final novel-locus split: %d reference loci retained in novel_loci.fa "
        "(%d protected unlifted from %s); %d other unlifted sequences retained "
        "in %s; final size/composition filter removed %d mapped and %d "
        "unmapped cores (mapped: unmasked >= %d and >= max(%d, %.3g * "
        "total); unmapped: total > %d and unmasked >= max(%d, %.3g * total))",
        len(reference_records), protected_unmapped_count,
        ",".join(sorted(protected_samples)) or "(none)",
        len(unmapped_records), UNMAPPED_NOVEL_FASTA,
        failed_mapped_filter, failed_unmapped_filter,
        args.promote_size, args.min_segment, args.min_unmasked_fraction,
        args.unmapped_main_size, args.min_segment, args.min_unmasked_fraction,
    )
    retained_records = reference_records + unmapped_records
    final_loci = [piece for piece, _result in reference_records]
    checkpoint = [
        ("novel_locus", piece, result) for piece, result in reference_records
    ] + [
        ("unmapped_novel", piece, result) for piece, result in unmapped_records
    ]
    result_rows = [
        table_row(piece, "novel_locus", result)
        for piece, result in reference_records
    ] + [
        table_row(piece, "unmapped_novel", result)
        for piece, result in unmapped_records
    ]
    write_discovery_bundle(
        idir, retained_records, final_loci, checkpoint_records=checkpoint,
    )
    write_annotated_loci_fasta(
        reference_records, os.path.join(idir, "promoted.fa"), use_final_id=True,
    )
    _write_stage_results(idir, result_rows)
    compact_discovery_iteration(idir, args.keep_work)
    try:
        os.remove(_stage_marker(idir, "_IN_PROGRESS"))
    except FileNotFoundError:
        pass
    _write_stage_marker(idir, "_READY", marker)
    _commit_novel_loci_stage(
        idir, checkpoint, novel_fa,
        os.path.join(outdir, UNMAPPED_NOVEL_FASTA), table_path, marker,
    )
    log.info(
        "Novel-locus selection complete: %d reference loci, %d separately "
        "retained unlifted records",
        len(reference_records), len(unmapped_records),
    )
    return reference_records


def align_insertions_to_novel_loci(
    assembly: AssemblyEntry,
    iteration: int,
    records: Sequence[Tuple[Piece, LocalizeResult]],
    novel_target: Optional[str],
    main_fa: str,
    args,
    outdir: str,
    table_path: str,
    log: logging.Logger,
) -> List[Tuple[Piece, LocalizeResult]]:
    """Subtract/lift every main residual against finalized novel loci.

    Large candidates that lost promotion remain in ``records``. They are only
    removed if this stage demonstrates representation by ``novel_loci.fa``;
    otherwise their residual sequence continues as a local insertion.
    """
    idir = os.path.join(
        outdir, "iterations", "novel_alignment",
        f"{iteration:04d}_{sanitize_id(assembly.sample)}",
    )
    Path(idir).mkdir(parents=True, exist_ok=True)
    marker = _assembly_stage_marker("novel_alignment", assembly, iteration, args)
    originals = [(dataclasses.replace(piece), result) for piece, result in records]
    candidates = [piece for piece, _result in originals]
    representatives: Dict[str, Tuple[str, int, int, str]] = {}
    novel_hits: List = []

    residual = candidates
    if novel_target and candidates:
        residual = exclude_against_reference_converge(
            candidates, novel_target, os.path.join(idir, "to_novel_loci"), args,
            log, representatives, blast_mode="residual", alignment_sink=novel_hits,
            broad_aligner="minimap2", run_blastn=False,
            max_cycles=1 + args.cycle_minimap2,
        )

    novel_blocks = write_lift_blocks(
        os.path.join(idir, LIFT_BLOCKS_NAME),
        [piece for piece, _result in originals], novel_hits, args.flank,
    )
    discovery_dir = _sample_stage_dir(outdir, "main_discovery", assembly.sample)
    main_blocks_path = os.path.join(discovery_dir, LIFT_BLOCKS_NAME)
    needed_contigs = {piece.contig for piece in residual}
    combined_blocks: Dict[str, List[AnchorBlock]] = {}
    if os.path.isfile(main_blocks_path):
        for contig, blocks in read_lift_blocks(main_blocks_path, needed_contigs).items():
            combined_blocks.setdefault(contig, []).extend(blocks)
    else:
        log.warning("Missing compact main lift evidence for %s; using checkpoint placements",
                    assembly.sample)
    for contig, blocks in novel_blocks.items():
        if contig not in needed_contigs:
            continue
        combined_blocks.setdefault(contig, []).extend(blocks)

    # Every split residual gets a new stable lift identity. This prevents two
    # fragments of one original insertion from overwriting each other's lift.
    safe_sample = sanitize_id(assembly.sample)
    for index, piece in enumerate(residual):
        new_lift = f"locallift_{iteration:04d}_{safe_sample}_{index:08d}"
        piece.lift_id = new_lift

    combined_results = localize_from_block_map(residual, combined_blocks, args.flank)
    selected: Dict[str, LocalizeResult] = {}
    for piece, result in zip(residual, combined_results):
        # A trimmed/split residual must earn its placement from blocks adjacent
        # to its final coordinates. Never inherit a stale pre-trimming lift.
        selected[piece.temp_id] = result
    expanded, result_by_lift = split_lift_placements(residual, selected)

    local_records: List[Tuple[Piece, LocalizeResult]] = []
    for index, piece in enumerate(expanded):
        if not qualifies_local_insertion(
            piece, args.min_segment, args.min_unmasked_fraction,
        ):
            continue
        piece.temp_id = f"localnovel_{iteration:04d}_{safe_sample}_{index:08d}"
        local_records.append((piece, result_by_lift[piece.lift_id]))

    survivor_origins = {piece.origin_id for piece, _result in local_records}
    result_rows: List[str] = []
    written_origins = set()
    for piece, _result in originals:
        if piece.origin_id in survivor_origins or piece.origin_id in written_origins:
            continue
        representative = representatives.get(piece.origin_id)
        if representative is None:
            continue
        written_origins.add(piece.origin_id)
        result_rows.append(redundant_row(piece, representative))

    checkpoint = [("local_candidate", piece, result)
                  for piece, result in local_records]
    _finish_assembly_stage(
        idir, marker, local_records, checkpoint, result_rows, main_fa,
        os.path.join(outdir, "small_insertions.fa"), table_path, args,
    )
    log.info("%s novel-locus alignment: %d -> %d local candidates",
             assembly.sample, len(records), len(local_records))
    return local_records


def run_novel_alignment_stage(
    assemblies: Sequence[AssemblyEntry],
    novel_target: Optional[str],
    main_fa: str,
    args,
    outdir: str,
    table_path: str,
    log: logging.Logger,
) -> List[Tuple[Piece, LocalizeResult]]:
    """Run the independently restartable residual-to-novel-loci stage.

    This is a top-level entry point so the expensive stage can be launched by
    ``align_remaining_to_novel_loci.py`` instead of being coupled to novel-locus
    construction. Its per-assembly checkpoints remain compatible with the
    original in-pipeline implementation.
    """
    main_config_signature = _restart_signature(args)
    config_signature = _downstream_signature(args)
    small_fa = os.path.join(outdir, "small_insertions.fa")
    assembly_items = list(enumerate(assemblies[1:], start=1))

    def worker(item: Tuple[int, AssemblyEntry], attempt: int):
        iteration, assembly = item
        idir = os.path.join(
            outdir, "iterations", "novel_alignment",
            f"{iteration:04d}_{sanitize_id(assembly.sample)}",
        )
        state = _load_or_commit_assembly_stage(
            idir, assembly, iteration, "novel_alignment", main_fa, small_fa,
            table_path, log, config_signature,
        )
        if state is not None:
            read_lift_blocks(os.path.join(idir, LIFT_BLOCKS_NAME), ())
            return [
                (piece, result) for category, piece, result in state
                if category == "local_candidate"
            ]
        if os.path.isdir(idir):
            shutil.rmtree(idir)
        discovery_dir = _sample_stage_dir(outdir, "main_discovery", assembly.sample)
        if not _load_main_discovery_checkpoint(
            discovery_dir, assembly, iteration, main_fa, small_fa,
            table_path, log, main_config_signature,
        ):
            raise RuntimeError(
                f"Missing completed main-discovery checkpoint for {assembly.sample}"
            )
        assembly_records = [
            (piece, result)
            for category, piece, result in _iter_resume_state(discovery_dir)
            if category == "insertion"
        ]
        assembly_records = _set_current_assembly_rank(assembly_records, iteration)
        return align_insertions_to_novel_loci(
            assembly, iteration, assembly_records, novel_target,
            main_fa, args, outdir, table_path, log,
        )

    local_results = run_parallel_with_retries(
        assembly_items, worker, args.jobs, args.retries,
        "novel-locus-alignment", log,
    )
    local_records = [
        record
        for iteration, _assembly in assembly_items
        for record in local_results.get(iteration, [])  # type: ignore[union-attr]
    ]
    log.info(
        "Novel-locus alignment complete: %d lift-local candidates",
        len(local_records),
    )
    return local_records


def clean_local_insertions_v3(
    records: Sequence[Tuple[Piece, LocalizeResult]],
    main_fa: str,
    novel_fa: str,
    small_fa: str,
    args,
    outdir: str,
    counters: Dict[str, int],
    table_path: str,
    log: logging.Logger,
) -> List[Piece]:
    """Clean small/local insertions only against their lifted reference locus."""
    idir = os.path.join(outdir, "iterations", "local_insertion_cleanup")
    Path(idir).mkdir(parents=True, exist_ok=True)
    marker = {
        "kind": "local_insertion_cleanup",
        "record_count": len(records),
        "config_signature": _downstream_signature(args),
        "assembly_order_signature": getattr(args, "assembly_order_signature", None),
        "logic_version": LOCAL_CLEANUP_LOGIC_VERSION,
    }

    target_path: Dict[str, str] = {}
    target_lengths: Dict[str, int] = {}
    for fasta_path in (main_fa, novel_fa):
        if not os.path.isfile(fasta_path) or os.path.getsize(fasta_path) == 0:
            continue
        with IndexedFasta(fasta_path) as reader:
            for name in reader.names():
                if name in target_path:
                    raise RuntimeError(f"Reference id occurs in main and novel loci: {name}")
                target_path[name] = fasta_path
                target_lengths[name] = reader.length(name)

    windows, unmapped = build_cohort_windows(
        records, [], target_lengths, extension=args.local_extension,
        preferred_size=args.cohort_window_size,
    )
    log.info("Local insertion cleanup: %d windows, %d unmapped records ignored",
             len(windows), len(unmapped))
    cleanup_threads = _parallel_task_threads(args, len(windows))

    def clean_window(item: Tuple[int, object], attempt: int):
        window_index, raw_window = item
        window = raw_window  # type: ignore[assignment]
        workdir = os.path.join(idir, "windows", f"window{window_index:06d}")
        if attempt:
            shutil.rmtree(workdir, ignore_errors=True)
        local_args = argparse.Namespace(**vars(args))
        local_args.threads = cleanup_threads
        fasta_path = target_path.get(window.main)
        if fasta_path is None:
            raise RuntimeError(f"No reference FASTA for lifted target {window.main}")
        local_fa = os.path.join(workdir, "local_reference.fa")
        _nearby, coordinate_map = write_local_reference(fasta_path, window, local_fa)
        rep_kmers = getattr(local_args, "winnow_rep_kmers", None)
        representatives: Dict[str, Tuple[str, int, int, str]] = {}
        survivors = exclude_against_reference_converge(
            window.pieces, local_fa, os.path.join(workdir, "to_local"),
            local_args, log, representatives, blast_mode="independent",
            broad_aligner="winnowmap", winnow_rep_kmers=rep_kmers,
            target_coordinate_map=coordinate_map, max_cycles=1,
            warn_on_limit=False,
        )
        self_fn = make_self_align_fn(
            os.path.join(workdir, "selfclean"), local_args, log,
            blast_mode="independent", broad_aligner="winnowmap",
        )
        survivors = self_clean_converge(
            survivors, self_fn, local_args.min_segment, "insertion", log,
            max_cycles=1, warn_on_limit=False,
            sample_priorities=getattr(local_args, "sample_priorities", None),
        )
        return survivors, representatives

    survivors: List[Piece] = []
    representatives: Dict[str, Tuple[str, int, int, str]] = {}
    if windows:
        cleanup_workers = min(_local_jobs(args), len(windows))
        log.info("Local insertion cleanup: %d parallel jobs x %d threads",
                 cleanup_workers, cleanup_threads)
        results = run_group_jobs_with_retries(
            list(enumerate(windows)), clean_window, _local_jobs(args),
            args.retries, "local-insertion-clean", log,
        )
        for index in sorted(results):
            kept, local_representatives = results[index]  # type: ignore[misc]
            survivors.extend(kept)
            representatives.update(local_representatives)

    # Re-lift only the final cleaned sequences. Combine the compact main blocks
    # saved during independent discovery with the novel-locus blocks saved after
    # the cohort-wide novel reference was finalized.
    prior_result_by_lift = {piece.lift_id: result for piece, result in records}
    final_results: List[Optional[LocalizeResult]] = [None] * len(survivors)
    by_assembly: Dict[Tuple[int, str], List[int]] = {}
    for index, piece in enumerate(survivors):
        match = re.match(r"^locallift_([0-9]{4})_", piece.lift_id)
        if match:
            by_assembly.setdefault((int(match.group(1)), piece.sample), []).append(index)
        else:
            final_results[index] = prior_result_by_lift.get(piece.lift_id)
    for (iteration, sample), indices in by_assembly.items():
        safe_sample = sanitize_id(sample)
        main_blocks_path = os.path.join(
            _sample_stage_dir(outdir, "main_discovery", sample),
            LIFT_BLOCKS_NAME,
        )
        novel_blocks_path = os.path.join(
            outdir, "iterations", "novel_alignment",
            f"{iteration:04d}_{safe_sample}", LIFT_BLOCKS_NAME,
        )
        assembly_pieces = [survivors[index] for index in indices]
        needed_contigs = {piece.contig for piece in assembly_pieces}
        combined: Dict[str, List[AnchorBlock]] = {}
        for path in (main_blocks_path, novel_blocks_path):
            for contig, blocks in read_lift_blocks(path, needed_contigs).items():
                combined.setdefault(contig, []).extend(blocks)
        for index, result in zip(
            indices, localize_from_block_map(assembly_pieces, combined, args.flank),
        ):
            piece = survivors[index]
            previous = prior_result_by_lift.get(piece.lift_id)
            # A placement clone was cleaned only in its own local window. Do not
            # reattach a sibling location whose clone was removed as redundant;
            # if both clones survive, consolidation below restores both.
            if "_placement" in piece.lift_id and previous and previous.placements:
                allowed = {
                    (placement.main, placement.strand) for placement in previous.placements
                }
                kept = [
                    placement for placement in result.placements
                    if (placement.main, placement.strand) in allowed
                ]
                result = LocalizeResult(
                    previous.status if kept else "unmapped", kept,
                    note=";".join(filter(None, [result.note, "placement_specific_clean"])),
                )
            final_results[index] = result

    lifted_survivors = consolidate_identical_records([
        (piece, result)
        for piece, result in zip(survivors, final_results)
        if result is not None and result.placements
        and qualifies_local_insertion(piece, args.min_segment, args.min_unmasked_fraction)
    ])
    final_records: List[Tuple[Piece, LocalizeResult]] = []
    result_rows: List[str] = []
    for piece, result in lifted_survivors:
        assign_final_ids([piece], f"ins_{piece.sample}", counters)
        final_records.append((piece, result))
        result_rows.append(table_row(piece, "small", result))

    kept_origins = {piece.origin_id for piece, _result in final_records}
    written_origins = set()
    for piece, _result in records:
        if piece.origin_id in kept_origins or piece.origin_id in written_origins:
            continue
        written_origins.add(piece.origin_id)
        result_rows.append(redundant_row(piece, representatives.get(piece.origin_id)))

    checkpoint = [("small", piece, result) for piece, result in final_records]
    write_discovery_bundle(
        idir, final_records, [], checkpoint_records=checkpoint,
    )
    _write_stage_results(idir, result_rows)
    compact_discovery_iteration(idir, args.keep_work)
    _write_stage_marker(idir, "_READY", marker)
    _commit_stage(idir, checkpoint, main_fa, small_fa, table_path, marker)
    log.info("Local insertion cleanup complete: %d retained", len(final_records))
    return [piece for piece, _result in final_records]


# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------

def setup_logger(outdir: str, verbose: bool, resume: bool = False) -> logging.Logger:
    Path(outdir, "logs").mkdir(parents=True, exist_ok=True)
    log = logging.getLogger("minsetref")
    log.setLevel(logging.DEBUG)
    for handler in log.handlers:
        handler.close()
    log.handlers.clear()
    fmt = logging.Formatter("[%(asctime)s] %(levelname)s: %(message)s", "%Y-%m-%d %H:%M:%S")
    sh = logging.StreamHandler(sys.stderr)
    sh.setLevel(logging.DEBUG if verbose else logging.INFO)
    sh.setFormatter(fmt)
    fh = logging.FileHandler(Path(outdir, "logs", "pipeline.log"),
                             mode="at" if resume else "wt")
    fh.setFormatter(fmt)
    log.addHandler(sh)
    log.addHandler(fh)
    return log


def parse_args(argv: Optional[Sequence[str]] = None):
    p = argparse.ArgumentParser(
        description="Min-set reference builder v3 (global novel loci + lift-local insertions).",
    )
    p.add_argument("-l", "--assembly-list", required=True)
    p.add_argument("-o", "--outdir", required=True)
    p.add_argument("--precleaned-reference",
                   help="use this cleaned FASTA as main, skip self-cleaning, and skip assembly-list row 1")
    p.add_argument(
        "--BlacklistRegion", "--blacklist-region",
        dest="blacklist_region", metavar="BED",
        help=(
            "BED3+ file of source-assembly regions to exclude (standard "
            "0-based, half-open coordinates); blacklisted sequence is "
            "subtracted before insertion filtering"
        ),
    )
    p.add_argument("--window-size", type=int, default=100_000)
    p.add_argument("-j", "--jobs", type=int, default=1,
                   help="independent assembly/local jobs to run concurrently")
    p.add_argument(
        "--candidate-shard-jobs", type=int,
        help=(
            "parallel workers used only while filtering per-assembly main-discovery "
            "checkpoints into candidate shards (default: -j/--jobs)"
        ),
    )
    p.add_argument("-t", "--threads", type=int, default=8,
                   help="threads used by each concurrent alignment job")
    p.add_argument("--retries", type=int, default=2,
                   help="automatic retries after a failed parallel job (default: 2)")
    p.add_argument("--processes", type=int, default=1, help="worker processes for windowing")
    p.add_argument("--min-segment", type=int, default=300, help="valid segment + island-merge cutoff (bp)")
    p.add_argument("--min-identity", type=float, default=95.0)
    p.add_argument("--promote-size", type=int, default=10_000,
                   help="unmasked bp required for a mapped novel-locus candidate")
    p.add_argument("--flank", type=int, default=5_000, help="localization flank (bp)")
    p.add_argument("--local-extension", type=int, default=3_000,
                   help="target-contig bases added on each side of a local insertion lift")
    p.add_argument("--novel-locus-extension", type=int, default=30_000,
                   help="deprecated compatibility option; no longer used")
    p.add_argument("--cohort-window-size", type=int, default=100_000,
                   help="preferred cohort window span; cuts occur only at unspanned gaps")
    p.add_argument("--local-group-jobs", type=int,
                   help="deprecated local-only override for -j/--jobs")
    p.add_argument("--local-group-threads", type=int,
                   help="deprecated local-only override for -t/--threads")
    p.add_argument("--unmapped-main-size", type=int, default=10_000,
                   help="unmapped sequences longer than this enter global cleaning; "
                        "final unlifted representatives normally go to "
                        "unmapped_novel.fa")
    p.add_argument(
        "--keep-unmapped-samples",
        default=",".join(DEFAULT_KEEP_UNMAPPED_SAMPLES),
        metavar="SAMPLE[,SAMPLE...]",
        help=(
            "comma-separated samples whose retained unlifted representatives "
            "remain in novel_loci.fa (default: %(default)s); pass an empty "
            "value to exempt no samples"
        ),
    )
    p.add_argument("--min-unmasked-fraction", type=float, default=0.2,
                   help="minimum unmasked fraction for any promoted locus")
    p.add_argument(
        "--prior", default=DEFAULT_PRIOR,
        help=(
            "ordered comma-separated priority tiers (default: %(default)s); "
            "an entry without '_' matches all ENTRY_* haplotypes"
        ),
    )
    p.add_argument("--max-cycles", type=int, default=10,
                   help="maximum cycles for capped cleaning phases; novel-locus "
                        "global_initial always runs to convergence")
    p.add_argument("--cycle-minimap2", type=int, default=1, metavar="N",
                   help="maximum additional residual minimap2 cycles in capped phases "
                        "(global_initial is uncapped)")
    p.add_argument("--blast-word-size", type=int, default=50)
    p.add_argument("--blast-evalue", default="1e-300")
    p.add_argument("--blast-max-target-seqs", type=int, default=100)
    p.add_argument("--no-repkmers", action="store_true", help="skip meryl/winnowmap -W repetitive db")
    p.add_argument("--skip-blastn", action="store_true", help="debug: skip BLASTN passes")
    p.add_argument(
        "--full-main-discovery-cleaning", action="store_true",
        help=(
            "after per-assembly Minimap2 exclusion, also run Winnowmap+BLASTN "
            "against the fixed complete reference to convergence before cohort self-cleaning"
        ),
    )
    p.add_argument("--keep-work", action="store_true",
                   help="retain cycle FASTAs, SAMs, and alignment databases instead of compacting")
    p.add_argument("--resume", action="store_true",
                   help="continue OUT, committing READY steps and skipping SUCCESS assemblies")
    p.add_argument(
        "--stop-after-novel-loci", action="store_true",
        help=(
            "deprecated compatibility flag; stopping after finalized novel "
            "loci is now the default"
        ),
    )
    p.add_argument(
        "--run-novel-alignment-in-pipeline", action="store_true",
        help=(
            "legacy opt-in: run residual-to-novel alignment and local cleanup "
            "inside minsetref.py; normally use "
            "align_remaining_to_novel_loci.py separately"
        ),
    )
    p.add_argument("--force", action="store_true")
    p.add_argument("--verbose", action="store_true")
    args = p.parse_args(argv)
    args.keep_unmapped_samples = _keep_unmapped_samples(args)
    for name in ("window_size", "jobs", "threads", "processes", "min_segment", "promote_size",
                 "flank", "local_extension", "novel_locus_extension", "unmapped_main_size",
                 "cohort_window_size", "max_cycles"):
        if getattr(args, name) <= 0:
            p.error(f"--{name.replace('_','-')} must be > 0")
    if args.cycle_minimap2 < 0:
        p.error("--cycle-minimap2 must be >= 0")
    if args.retries < 0:
        p.error("--retries must be >= 0")
    for name in ("candidate_shard_jobs", "local_group_jobs", "local_group_threads"):
        value = getattr(args, name)
        if value is not None and value <= 0:
            p.error(f"--{name.replace('_','-')} must be > 0")
    if not (0 <= args.min_unmasked_fraction <= 1):
        p.error("--min-unmasked-fraction must be in [0,1]")
    if not (0 < args.min_identity <= 100):
        p.error("--min-identity must be in (0,100]")
    if args.resume and args.force:
        p.error("--resume and --force are mutually exclusive")
    try:
        configure_blacklist_args(args)
    except (OSError, ValueError) as error:
        p.error(str(error))
    return args


def check_dependencies(args) -> None:
    ensure_executable("minimap2")
    # Winnowmap is used for first-reference cleaning, genome-wide/local
    # novel-locus cleaning, and final lift-local insertion cleaning.
    ensure_executable("winnowmap")
    if not args.no_repkmers:
        ensure_executable("meryl")
    if not args.skip_blastn:
        ensure_executable("blastn")
        ensure_executable("makeblastdb")


def main(argv: Optional[Sequence[str]] = None) -> int:
    args = parse_args(argv)
    if args.precleaned_reference:
        source = os.path.realpath(args.precleaned_reference)
        out_real = os.path.realpath(args.outdir)
        try:
            source_inside_output = os.path.commonpath([source, out_real]) == out_real
        except ValueError:
            source_inside_output = False
        if source_inside_output:
            raise ValueError("--precleaned-reference cannot be inside --outdir; --force could delete it")
    if args.resume and not os.path.isdir(args.outdir):
        raise FileNotFoundError(f"Cannot resume; output folder does not exist: {args.outdir}")
    if os.path.exists(args.outdir) and not args.force and not args.resume:
        raise FileExistsError(f"Output exists; use --force or --resume: {args.outdir}")
    if os.path.exists(args.outdir) and args.force:
        shutil.rmtree(args.outdir)
    Path(args.outdir).mkdir(parents=True, exist_ok=True)

    log = setup_logger(args.outdir, args.verbose, resume=args.resume)
    assemblies = read_assembly_list(args.assembly_list)
    args.main_discovery_reference_signature = (
        main_discovery_reference_signature(assemblies[0])
    )
    check_dependencies(args)
    log.info("Loaded %d assemblies", len(assemblies))
    if args.blacklist_region:
        log.info(
            "Loaded blacklist BED %s: %d contigs, %d merged intervals",
            args.blacklist_region,
            len(args.blacklist_intervals),
            sum(len(index.intervals) for index in args.blacklist_intervals.values()),
        )
    args.sample_priorities, unmatched_prior = build_sample_priorities(
        assemblies, args.prior,
    )
    if args.sample_priorities:
        log.info(
            "Prior-scored samples: %s",
            ", ".join(
                f"{sample}={bonus // PRIOR_SCORE_STEP}x10G"
                for sample, bonus in sorted(
                    args.sample_priorities.items(), key=lambda item: -item[1]
                )
            ),
        )
    if unmatched_prior:
        log.warning("--prior entries matched no samples: %s", ", ".join(unmatched_prior))
    # Priority no longer depends on ordinary query-list order.  This signature
    # tracks cohort membership/path identity only; --prior is in config_signature.
    order_payload = "\n".join(sorted(
        f"{assembly.sample}\0{os.path.realpath(assembly.path)}"
        for assembly in assemblies
    )).encode()
    args.assembly_order_signature = hashlib.sha256(order_payload).hexdigest()
    args.winnow_rep_kmers = None
    main_config_signature = _restart_signature(args)
    config_signature = _downstream_signature(args)
    format_path = os.path.join(args.outdir, "checkpoint_format.json")
    if args.resume:
        if not os.path.isfile(format_path):
            raise RuntimeError(
                "This output folder predates restart checkpoints and cannot be resumed "
                "safely. Start a new output folder (or rebuild it with --force)."
            )
        try:
            with open(format_path) as handle:
                checkpoint_format = json.load(handle)
        except (OSError, ValueError, TypeError) as exc:
            raise RuntimeError(f"Invalid checkpoint format marker: {format_path}") from exc
        if checkpoint_format.get("version") != CHECKPOINT_VERSION:
            raise RuntimeError(
                f"Unsupported checkpoint format {checkpoint_format.get('version')!r}"
            )
    else:
        _atomic_json(format_path, {"version": CHECKPOINT_VERSION})

    table_path = os.path.join(args.outdir, "insertions.tsv")
    main_fa = os.path.join(args.outdir, "main_chroms.fa")
    filtered_reference_fa = os.path.join(args.outdir, FILTERED_REFERENCE_FASTA)
    novel_fa = os.path.join(args.outdir, "novel_loci.fa")
    unmapped_novel_fa = os.path.join(args.outdir, UNMAPPED_NOVEL_FASTA)
    small_fa = os.path.join(args.outdir, "small_insertions.fa")

    reference_marker_path = os.path.join(args.outdir, "reference._SUCCESS")
    expected_reference = {
        "kind": "reference",
        "sample": assemblies[0].sample,
        "assembly_path": os.path.realpath(assemblies[0].path),
        "assembly_size": os.path.getsize(assemblies[0].path),
        "assembly_mtime_ns": os.stat(assemblies[0].path).st_mtime_ns,
        "precleaned_reference": (
            os.path.realpath(args.precleaned_reference) if args.precleaned_reference else None
        ),
        "precleaned_size": (
            os.path.getsize(args.precleaned_reference) if args.precleaned_reference else None
        ),
        "precleaned_mtime_ns": (
            os.stat(args.precleaned_reference).st_mtime_ns if args.precleaned_reference else None
        ),
        "config_signature": main_config_signature,
    }
    reference_marker = None
    if args.resume and os.path.isfile(reference_marker_path):
        try:
            with open(reference_marker_path) as handle:
                reference_marker = json.load(handle)
        except (OSError, ValueError, TypeError) as exc:
            raise RuntimeError(f"Invalid reference checkpoint: {reference_marker_path}") from exc
        if not os.path.isfile(main_fa):
            raise RuntimeError(f"Reference checkpoint exists but FASTA is missing: {main_fa}")
        recorded_base_size = reference_marker.get("base_size")
        installed_reference_verified = (
            isinstance(recorded_base_size, int)
            and os.path.getsize(main_fa) == recorded_base_size
        )
        if recorded_base_size is not None and not installed_reference_verified:
            raise RuntimeError(
                "Installed cleaned-reference size does not match its completed "
                f"checkpoint: checkpoint={recorded_base_size!r}, "
                f"current={os.path.getsize(main_fa)!r}"
            )
        legacy_reference_signature = _legacy_restart_signature_v1(args)
        ignored_precleaned_metadata = []
        for key, value in expected_reference.items():
            if (
                key == "config_signature"
                and reference_marker.get(key) == legacy_reference_signature
                and args.blacklist_region_signature is None
            ):
                # The first-reference self-cleaning algorithm did not change.
                # Reuse that expensive result while invalidating/recomputing
                # main-discovery checkpoints under logic version 2.
                continue
            if (
                key in {"precleaned_size", "precleaned_mtime_ns"}
                and installed_reference_verified
                and reference_marker.get("precleaned_reference")
                == expected_reference.get("precleaned_reference")
                and reference_marker.get(key) != value
            ):
                # Once the precleaned reference has been copied into outdir and
                # atomically checkpointed, that installed FASTA is authoritative.
                # The external source may later be replaced in place without
                # changing the completed run. Main discovery is independently
                # fingerprinted against the complete first assembly.
                ignored_precleaned_metadata.append(key)
                continue
            if reference_marker.get(key) != value:
                raise RuntimeError(
                    f"Resume reference mismatch for {key}: checkpoint={reference_marker.get(key)!r}, "
                    f"current={value!r}"
                )
        if ignored_precleaned_metadata:
            log.warning(
                "RESUME: external precleaned-reference metadata changed (%s); "
                "reusing verified installed cleaned reference %s",
                ", ".join(ignored_precleaned_metadata), main_fa,
            )
        if not args.no_repkmers:
            args.winnow_rep_kmers = build_winnowmap_rep_kmers(
                assemblies[0].path,
                os.path.join(args.outdir, "indexes", f"reference_repetitive_k{WINNOW_PARAMS['k']}.txt"),
                WINNOW_PARAMS["k"], log,
            )
        log.info("RESUME: using completed main reference %s", main_fa)
        if not args.precleaned_reference:
            if not os.path.isfile(filtered_reference_fa):
                sample, retained_by_contig = read_cleaned_reference_intervals(
                    main_fa, expected_sample=assemblies[0].sample,
                )
                write_filtered_reference_fasta(
                    assemblies[0].path, sample, retained_by_contig,
                    filtered_reference_fa, log,
                )
            first_dir = os.path.join(
                args.outdir, "iterations",
                f"000_{sanitize_id(assemblies[0].sample)}_main",
            )
            removed_fastas, removed_bytes = compact_reference_iteration(
                first_dir, args.keep_work,
            )
            if removed_fastas:
                log.info(
                    "Reference cleanup removed %d intermediate FASTA files (%.2f GiB)",
                    removed_fastas, removed_bytes / (1024 ** 3),
                )
    else:
        # A reference without its atomic success marker is incomplete.  It is
        # safe to rebuild because later assembly stages cannot become READY
        # before the reference checkpoint exists.
        if args.resume:
            for path in (
                main_fa, main_fa + ".fai", filtered_reference_fa,
                filtered_reference_fa + ".fai",
            ):
                try:
                    os.remove(path)
                except FileNotFoundError:
                    pass
            first_dir = os.path.join(
                args.outdir, "iterations", f"000_{sanitize_id(assemblies[0].sample)}_main",
            )
            shutil.rmtree(first_dir, ignore_errors=True)
            log.info("RESUME: incomplete reference step will be rebuilt")
        if not args.no_repkmers:
            args.winnow_rep_kmers = build_winnowmap_rep_kmers(
                assemblies[0].path,
                os.path.join(args.outdir, "indexes", f"reference_repetitive_k{WINNOW_PARAMS['k']}.txt"),
                WINNOW_PARAMS["k"], log,
            )
        if args.precleaned_reference:
            log.info("Skipping assembly-list row 1 (%s); using precleaned reference %s",
                     assemblies[0].sample, args.precleaned_reference)
            _main_pieces, main_fa, _rep_kmers = install_precleaned_main(
                args.precleaned_reference, args, args.outdir, log,
            )
        else:
            _main_pieces, main_fa, _rep_kmers = build_main(
                assemblies[0], args, args.outdir, log,
            )
        expected_reference["base_size"] = os.path.getsize(main_fa)
        _atomic_json(reference_marker_path, expected_reference)

    # These files are derived, cumulative views of immutable per-stage
    # checkpoints. Rebuild them on every resume after validating the reference,
    # then idempotently recommit every completed stage below. This prevents a
    # missing/truncated output from being mistaken for a valid completed run.
    if args.resume:
        log.info(
            "RESUME: rebuilding cumulative novel/unmapped/small FASTAs "
            "and insertion table"
        )
    with open(table_path, "w") as tbl:
        tbl.write(TABLE_HEADER)
    open(novel_fa, "w").close()
    open(unmapped_novel_fa, "w").close()
    open(small_fa, "w").close()
    _TABLE_ID_CACHE.pop(os.path.realpath(table_path), None)
    _FASTA_ID_CACHE.pop(os.path.realpath(novel_fa), None)
    _FASTA_ID_CACHE.pop(os.path.realpath(unmapped_novel_fa), None)
    _FASTA_ID_CACHE.pop(os.path.realpath(small_fa), None)

    total_worker_threads = args.jobs * args.threads
    available_cpus = os.cpu_count()
    log.info("Parallel policy: %d jobs x %d threads = %d worker threads; retries=%d",
             args.jobs, args.threads, total_worker_threads, args.retries)
    visible_memory, memory_source = _visible_memory_info()
    if visible_memory is not None:
        log.info("Effective memory limit: %.1f GiB from %s (%.1f GiB/job at %d jobs)",
                 visible_memory, memory_source, visible_memory / args.jobs, args.jobs)
    if available_cpus and total_worker_threads > available_cpus:
        log.warning("Requested %d worker threads but this process sees %d CPUs",
                    total_worker_threads, available_cpus)

    index_dir = os.path.join(args.outdir, "indexes")
    main_index = build_minimap2_index(
        main_fa, os.path.join(index_dir, "main_chroms.asm5.mmi"), log,
        preset=MINIMAP_PARAMS["preset"],
    )
    discovery_index = build_minimap2_index(
        assemblies[0].path,
        os.path.join(index_dir, "full_first_reference.asm5.mmi"),
        log,
        preset=MINIMAP_PARAMS["preset"],
    )
    counters = _seed_counters_from_outputs(main_fa, novel_fa, small_fa)
    assembly_items = list(enumerate(assemblies[1:], start=1))

    # Phase 1: every assembly sees the exact same complete first assembly.
    # The self-cleaned main_chroms.fa is deliberately not used here: sequence
    # removed or split during reference self-cleaning is still reference, not
    # biological novelty in a later assembly.
    def main_discovery_worker(item: Tuple[int, AssemblyEntry], attempt: int):
        iteration, assembly = item
        idir = _sample_stage_dir(args.outdir, "main_discovery", assembly.sample)
        finished = _load_main_discovery_checkpoint(
            idir, assembly, iteration, main_fa, small_fa,
            table_path, log, main_config_signature,
        )
        if finished:
            return None
        if os.path.isdir(idir):
            shutil.rmtree(idir)
        discover_main_novelty(
            assembly, iteration, discovery_index, main_fa, args, args.outdir,
            table_path, log,
        )
        return None

    run_parallel_with_retries(
        assembly_items, main_discovery_worker, args.jobs, args.retries,
        "main-discovery", log,
    )
    log.info("Main-discovery stage complete for all %d query assemblies", len(assembly_items))

    # Checkpoints written immediately before order fingerprints were introduced
    # can still be resumed safely when every stored numeric rank agrees with the
    # current list.  A reordered list deliberately fails this compatibility test.
    legacy_order_matches = True
    for iteration, assembly in assembly_items:
        stage_dir = _sample_stage_dir(args.outdir, "main_discovery", assembly.sample)
        stage_marker = (
            _read_stage_marker(stage_dir, "_SUCCESS")
            or _read_stage_marker(stage_dir, "_READY")
        )
        if stage_marker is None or stage_marker.get("iteration") != iteration:
            legacy_order_matches = False
            break

    def order_checkpoint_matches(marker: dict) -> bool:
        observed = marker.get("assembly_order_signature")
        if observed is not None:
            return observed == args.assembly_order_signature
        return legacy_order_matches

    # Phase 2 begins only after the Stage-1 barrier. Each assembly independently
    # selects its promotion candidates into a restartable FASTA+table shard.
    candidate_input_dir = os.path.join(args.outdir, "iterations", "novel_locus_input")

    def candidate_shard_paths(iteration: int, assembly: AssemblyEntry) -> Tuple[str, str, str]:
        shard_dir = os.path.join(
            candidate_input_dir, "shards", sanitize_id(assembly.sample),
        )
        return (
            shard_dir,
            os.path.join(shard_dir, "candidates.fa"),
            os.path.join(shard_dir, "candidates.tsv"),
        )

    def large_input_worker(item: Tuple[int, AssemblyEntry], attempt: int):
        iteration, assembly = item
        idir = _sample_stage_dir(args.outdir, "main_discovery", assembly.sample)
        shard_dir, shard_fasta, shard_table = candidate_shard_paths(iteration, assembly)
        marker = _read_stage_marker(shard_dir, "_SUCCESS") if args.resume else None
        if (marker is not None
                and marker.get("kind") == "novel_locus_input_shard"
                and marker.get("logic_version") == NOVEL_INPUT_SHARD_VERSION
                and marker.get("sample") == assembly.sample
                and marker.get("config_signature") == main_config_signature
                and os.path.isfile(shard_fasta)
                and os.path.isfile(shard_table)):
            return int(marker["candidate_count"]), int(marker["record_count"])
        log.info(
            "Building candidate shard for %s", assembly.sample,
        )
        record_count = 0

        def qualifying_records():
            nonlocal record_count
            for category, piece, result in _iter_resume_state(idir):
                if category != "insertion":
                    continue
                record_count += 1
                # Do not apply promotion-size or unmasked-fraction thresholds
                # before global cleaning/remerging.  A candidate only needs to
                # satisfy the base min_segment rule at this stage.
                if qualifies_preclean_candidate(piece, args.min_segment):
                    yield piece, result

        candidate_count = _write_candidate_bundle(
            shard_fasta, shard_table, qualifying_records(),
        )
        marker = {
            "kind": "novel_locus_input_shard",
            "logic_version": NOVEL_INPUT_SHARD_VERSION,
            "sample": assembly.sample,
            "record_count": record_count,
            "candidate_count": candidate_count,
            "config_signature": main_config_signature,
        }
        _write_stage_marker(shard_dir, "_SUCCESS", marker)
        log.info(
            "Finished candidate shard for %s: %d/%d records retained",
            assembly.sample, candidate_count, record_count,
        )
        return candidate_count, record_count

    candidate_shard_jobs = args.candidate_shard_jobs or args.jobs
    log.info(
        "Candidate-shard filtering: %d parallel workers (alignment -j remains %d)",
        candidate_shard_jobs, args.jobs,
    )
    shard_results = run_parallel_with_retries(
        assembly_items, large_input_worker, candidate_shard_jobs, args.retries,
        "novel-locus-input-shards", log,
    )
    total_discovered = sum(int(shard_results[i][1]) for i, _assembly in assembly_items)
    candidate_count = sum(int(shard_results[i][0]) for i, _assembly in assembly_items)

    # Merge shards in assembly-list order. The merge itself streams bytes and
    # never loads candidate sequences or anchor rows into a cohort-wide object.
    merged_candidate_fasta = os.path.join(candidate_input_dir, "all_candidates.fa")
    merged_candidate_table = os.path.join(candidate_input_dir, "all_candidates.tsv")
    ordered_shards = [
        candidate_shard_paths(iteration, assembly)[1:]
        for iteration, assembly in assembly_items
    ]
    log.info("Merging %d candidate shards into cohort-wide FASTA and table",
             len(ordered_shards))
    _merge_candidate_bundles(
        ordered_shards, merged_candidate_fasta, merged_candidate_table, log=log,
    )
    _write_stage_marker(candidate_input_dir, "_SUCCESS", {
        "kind": "novel_locus_input_merged",
        "logic_version": NOVEL_INPUT_SHARD_VERSION,
        "record_count": candidate_count,
        "source_record_count": total_discovered,
        "assembly_count": len(assembly_items),
        "config_signature": main_config_signature,
    })
    log.info("Independent main discovery complete: %d novel records; %d min-segment candidates",
             total_discovered, candidate_count)
    log.info("Merged candidate inputs: %s and %s",
             merged_candidate_table, merged_candidate_fasta)

    # Phase 2: choose globally non-redundant representatives, weighted-remerge
    # them, freshly lift 10-kb-anchored cores, and keep the occurrence catalog.
    novel_dir = os.path.join(args.outdir, "iterations", "novel_loci_cleanup")
    novel_success = _read_stage_marker(novel_dir, "_SUCCESS") if args.resume else None
    novel_ready = _read_stage_marker(novel_dir, "_READY") if args.resume else None
    novel_marker = novel_success or novel_ready
    keep_unmapped_samples = _keep_unmapped_samples(args)
    preserve_compatible_clean_cycles = False
    if (
        novel_marker is not None
        and NOVEL_CLEANUP_LOGIC_VERSION == 10
        and int(novel_marker.get("logic_version", 1)) in {8, 9}
        and novel_marker.get("kind") == "novel_loci_cleanup"
        and novel_marker.get("record_count") == candidate_count
        and novel_marker.get("config_signature") == config_signature
        and order_checkpoint_matches(novel_marker)
    ):
        old_state = _read_resume_state(novel_dir)
        prior_logic_version = int(novel_marker.get("logic_version", 1))
        upgraded_state, novel_marker = _upgrade_v8_novel_loci_split(
            novel_dir, old_state, novel_marker, args.keep_work,
            keep_unmapped_samples, args,
        )
        del old_state, upgraded_state
        log.info(
            "RESUME: upgraded novel-locus checkpoint %d -> %d by reapplying "
            "the final post-lift filter and mapped/unmapped split; global "
            "alignments were not rerun",
            prior_logic_version, NOVEL_CLEANUP_LOGIC_VERSION,
        )
    if (novel_marker is not None
            and int(novel_marker.get("logic_version", 1)) != NOVEL_CLEANUP_LOGIC_VERSION):
        preserve_compatible_clean_cycles = bool(
            int(novel_marker.get("logic_version", 1)) in {6, 8, 9, 10, 11}
            and novel_marker.get("kind") == "novel_loci_cleanup"
            and novel_marker.get("record_count") == candidate_count
            and novel_marker.get("config_signature") == config_signature
            and order_checkpoint_matches(novel_marker)
            and os.path.isdir(os.path.join(novel_dir, "global_initial"))
        )
        log.info(
            "RESUME: upgrading novel-locus cleanup checkpoint "
            "(logic version %d -> %d)%s",
            int(novel_marker.get("logic_version", 1)), NOVEL_CLEANUP_LOGIC_VERSION,
            "; preserving compatible global_initial cycles"
            if preserve_compatible_clean_cycles else "",
        )
        novel_marker = None
    if novel_marker is not None and not order_checkpoint_matches(novel_marker):
        log.info(
            "RESUME: invalidating novel-locus cleanup checkpoint because "
            "cohort membership changed",
        )
        novel_marker = None
    if (novel_marker is not None
            and novel_marker.get("config_signature") != config_signature):
        log.info(
            "RESUME: invalidating novel-locus cleanup checkpoint because "
            "configuration or --prior changed",
        )
        novel_marker = None
    if novel_marker is not None:
        if (novel_marker.get("kind") != "novel_loci_cleanup"
                or novel_marker.get("record_count") != candidate_count):
            raise RuntimeError("Novel-locus checkpoint does not match discovery records")
        saved_keep_unmapped = tuple(sorted(
            str(sample) for sample
            in novel_marker.get("keep_unmapped_samples", [])
        ))
        requested_keep_unmapped = tuple(sorted(keep_unmapped_samples))
        novel_state = _read_resume_state(novel_dir)
        if saved_keep_unmapped != requested_keep_unmapped:
            novel_state, novel_marker = _upgrade_v8_novel_loci_split(
                novel_dir, novel_state, novel_marker, args.keep_work,
                keep_unmapped_samples, args,
            )
            log.info(
                "RESUME: reclassified saved novel-locus records for "
                "--keep-unmapped-samples=%s; no alignments were rerun",
                ",".join(keep_unmapped_samples) or "(none)",
            )
        _commit_novel_loci_stage(
            novel_dir, novel_state, novel_fa, unmapped_novel_fa,
            table_path, novel_marker,
        )
        novel_records = [(piece, result) for category, piece, result in novel_state
                         if category == "novel_locus"]
        _update_counters_from_pieces(counters, [piece for piece, _ in novel_records])
        log.info("RESUME: skipped completed novel-locus construction")
    else:
        in_progress = (
            _read_stage_marker(novel_dir, "_IN_PROGRESS")
            if args.resume else None
        )
        in_progress_version = (
            int(in_progress.get("logic_version", 1))
            if in_progress is not None else None
        )
        compatible_prior_in_progress = bool(
            in_progress is not None
            and in_progress.get("kind") == "novel_loci_cleanup_in_progress"
            and in_progress_version in {6, 8, 9, 10, 11}
            and in_progress.get("record_count") == candidate_count
            and in_progress.get("config_signature") == config_signature
            and order_checkpoint_matches(in_progress)
            and os.path.isdir(os.path.join(novel_dir, "global_initial"))
        )
        preserve_compatible_clean_cycles = (
            preserve_compatible_clean_cycles or compatible_prior_in_progress
        )
        preserve_partial = bool(
            in_progress is not None
            and in_progress.get("kind") == "novel_loci_cleanup_in_progress"
            and in_progress_version == NOVEL_CLEANUP_LOGIC_VERSION
            and in_progress.get("record_count") == candidate_count
            and in_progress.get("config_signature") == config_signature
            and order_checkpoint_matches(in_progress)
        ) or preserve_compatible_clean_cycles
        if os.path.isdir(novel_dir) and not preserve_partial:
            shutil.rmtree(novel_dir)
        elif preserve_compatible_clean_cycles:
            # Compatible prior versions used the same independent full-candidate
            # self-clean. Preserve that expensive result, but rebuild the final
            # lift and mapped/unmapped output split.
            for name in os.listdir(novel_dir):
                if name == "global_initial":
                    continue
                path = os.path.join(novel_dir, name)
                if os.path.isdir(path):
                    shutil.rmtree(path)
                else:
                    os.remove(path)
            log.info(
                "RESUME: preserved compatible global_initial cycles; "
                "rebuilding weighted-gap remerge and fresh 10-kb lift",
            )
        elif preserve_partial:
            log.info("RESUME: preserving completed novel-locus cleanup cycles in %s",
                     novel_dir)
        open(novel_fa, "w").close()
        open(unmapped_novel_fa, "w").close()
        _FASTA_ID_CACHE.pop(os.path.realpath(novel_fa), None)
        _FASTA_ID_CACHE.pop(os.path.realpath(unmapped_novel_fa), None)
        large_candidates = _read_candidate_bundle(
            merged_candidate_fasta, merged_candidate_table,
        )
        novel_records = build_novel_loci(
            large_candidates, main_fa, novel_fa, args, args.outdir, counters,
            table_path, log, reload_main_lift_blocks=True,
            assembly_paths={assembly.sample: assembly.path for assembly in assemblies},
            complete_reference_target=discovery_index,
            complete_reference_fasta=assemblies[0].path,
        )
        del large_candidates

    if args.stop_after_novel_loci or not args.run_novel_alignment_in_pipeline:
        log.info(
            "Stopped after finalized novel loci. The optional "
            "residual-to-novel stage can be run with "
            "align_remaining_to_novel_loci.py. small_insertions.fa has not "
            "been finalized because local cleanup consumes that stage's output."
        )
        log.info(
            "Done through novel loci. main=%s novel_loci=%s "
            "unmapped_novel=%s",
            main_fa, novel_fa, unmapped_novel_fa,
        )
        return 0

    novel_target: Optional[str] = None
    if os.path.isfile(novel_fa) and os.path.getsize(novel_fa) > 0:
        novel_target = build_minimap2_index(
            novel_fa, os.path.join(index_dir, "novel_loci.asm5.mmi"), log,
            preset=MINIMAP_PARAMS["preset"],
        )

    # Phase 3: all original main residuals—including failed/duplicate promotion
    # candidates—independently align to the separate novel-locus reference.
    local_records = run_novel_alignment_stage(
        assemblies, novel_target, main_fa, args, args.outdir, table_path, log,
    )

    # Phase 4: local insertions are compared only with candidates sharing an
    # overlapping lift window on a main chromosome or a novel-locus contig.
    local_dir = os.path.join(args.outdir, "iterations", "local_insertion_cleanup")
    local_success = _read_stage_marker(local_dir, "_SUCCESS") if args.resume else None
    local_ready = _read_stage_marker(local_dir, "_READY") if args.resume else None
    local_marker = local_success or local_ready
    if (local_marker is not None
            and int(local_marker.get("logic_version", 1)) != LOCAL_CLEANUP_LOGIC_VERSION):
        log.info(
            "RESUME: invalidating local-cleanup checkpoint (logic version %d -> %d)",
            int(local_marker.get("logic_version", 1)), LOCAL_CLEANUP_LOGIC_VERSION,
        )
        local_marker = None
    if local_marker is not None and not order_checkpoint_matches(local_marker):
        log.info(
            "RESUME: invalidating local-cleanup checkpoint because "
            "cohort membership changed",
        )
        local_marker = None
    if (local_marker is not None
            and local_marker.get("config_signature") != config_signature):
        log.info(
            "RESUME: invalidating local-cleanup checkpoint because "
            "configuration or --prior changed",
        )
        local_marker = None
    if local_marker is not None:
        if (local_marker.get("kind") != "local_insertion_cleanup"
                or local_marker.get("record_count") != len(local_records)):
            raise RuntimeError("Local-cleanup checkpoint does not match local candidates")
        local_state = _read_resume_state(local_dir)
        _commit_stage(
            local_dir, local_state, main_fa, small_fa, table_path, local_marker,
        )
        log.info("RESUME: skipped completed local insertion cleanup")
    else:
        if os.path.isdir(local_dir):
            shutil.rmtree(local_dir)
        open(small_fa, "w").close()
        _FASTA_ID_CACHE.pop(os.path.realpath(small_fa), None)
        clean_local_insertions_v3(
            local_records, main_fa, novel_fa, small_fa, args, args.outdir,
            counters, table_path, log,
        )

    log.info("Done. main=%s novel_loci=%s unmapped_novel=%s small=%s table=%s",
             main_fa, novel_fa, unmapped_novel_fa, small_fa, table_path)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
