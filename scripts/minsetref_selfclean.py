#!/usr/bin/env python3
"""Self-clean with the simplified keep-earlier rule and a converge loop.

Redundancy rule (per confirmed spec): for a self-alignment where query interval
A aligns to target interval B, keep the EARLIER occurrence and remove the later
one over ``B - (A intersect B)``.  The intersection only matters for an internal
repeat (query and target are the same candidate); across two different
candidates the whole aligned target block is redundant.

Modes:
  * "insertion" (assemblies >= 2): keep the unaligned complement, bridging
    represented islands with < min_segment unmasked bases, then require
    unmasked >= min_segment.
  * "first_genome": redundant coverage is bridged across interior UNALIGNED
    islands carrying < min_segment MASKED bases (``find_redundant_and_kept``);
    kept segments are those with unmasked >= min_segment OR masked >= min_segment
    (so masked repeat blocks that block the bridge are preserved).

The self-align step is injected as a callable so the redundancy/keep logic can
be unit-tested with synthetic hits and reused with the real winnowmap+blastn
aligner in the orchestrator.
"""
from __future__ import annotations

import collections as cl
import dataclasses
import hashlib
import json
import logging
import os
import sys
import re
from pathlib import Path
from typing import Callable, Dict, List, Optional, Sequence, Tuple

from minsetref_core import (
    IndexedFasta,
    Piece,
    count_masked,
    count_unmasked,
    iter_n_free_segments,
    write_fasta_pieces,
)
from minsetref_segments import (
    find_redundant_and_kept,
    find_unaligned_segments,
    merge_intervals,
    novelty_score,
    subtract_intervals,
)

Interval = Tuple[int, int]

# A self-align callable takes the current pieces and a cycle index and returns
# alignment hits between candidate temp_ids (query_intervals/target_intervals are
# local to each candidate's sequence).
AlignFn = Callable[[Sequence[Piece], int], Sequence["object"]]
SELF_CLEAN_CONVERGENCE_BP = 500


def _pieces_signature(pieces: Sequence[Piece]) -> str:
    """Exact ordered signature for validating a completed cycle checkpoint."""
    digest = hashlib.sha256()
    for piece in pieces:
        fields = dataclasses.asdict(piece)
        sequence = fields.pop("seq")
        digest.update(json.dumps(fields, sort_keys=True, separators=(",", ":")).encode())
        digest.update(b"\0")
        digest.update(sequence.encode("ascii"))
        digest.update(b"\n")
    return digest.hexdigest()


def _checkpoint_input_signature(
    pieces: Sequence[Piece], align_fn: AlignFn,
) -> str:
    """Fingerprint both candidate sequence state and the alignment policy."""
    pieces_signature = _pieces_signature(pieces)
    checkpoint_tag = str(getattr(align_fn, "checkpoint_tag", ""))
    if not checkpoint_tag:
        return pieces_signature
    digest = hashlib.sha256()
    digest.update(checkpoint_tag.encode("utf-8"))
    digest.update(b"\0")
    digest.update(pieces_signature.encode("ascii"))
    return digest.hexdigest()


def _write_cycle_checkpoint(
    cycle_dir: str,
    pieces: Sequence[Piece],
    input_signature: str,
    removed_bp: int,
    converged: bool,
) -> None:
    Path(cycle_dir).mkdir(parents=True, exist_ok=True)
    fasta_path = os.path.join(cycle_dir, "kept.fa")
    state_path = os.path.join(cycle_dir, "kept.state.jsonl")
    marker_path = os.path.join(cycle_dir, "_SUCCESS.json")
    fasta_tmp = fasta_path + ".tmp"
    state_tmp = state_path + ".tmp"
    marker_tmp = marker_path + ".tmp"
    try:
        write_fasta_pieces(pieces, fasta_tmp)
        with open(state_tmp, "w") as state:
            for piece in pieces:
                fields = dataclasses.asdict(piece)
                fields.pop("seq", None)
                state.write(json.dumps(fields, sort_keys=True) + "\n")
        marker = {
            "input_signature": input_signature,
            "output_signature": _pieces_signature(pieces),
            "piece_count": len(pieces),
            "removed_bp": int(removed_bp),
            "converged": bool(converged),
        }
        with open(marker_tmp, "w") as out:
            json.dump(marker, out, sort_keys=True)
            out.write("\n")
            out.flush()
            os.fsync(out.fileno())
        os.replace(fasta_tmp, fasta_path)
        os.replace(state_tmp, state_path)
        os.replace(marker_tmp, marker_path)
    finally:
        for temporary in (fasta_tmp, state_tmp, marker_tmp):
            try:
                os.remove(temporary)
            except FileNotFoundError:
                pass


def _read_cycle_checkpoint(
    cycle_dir: str,
    input_signature: str,
) -> Tuple[List[Piece], int, bool] | None:
    fasta_path = os.path.join(cycle_dir, "kept.fa")
    state_path = os.path.join(cycle_dir, "kept.state.jsonl")
    marker_path = os.path.join(cycle_dir, "_SUCCESS.json")
    if not all(os.path.isfile(path) for path in (fasta_path, state_path, marker_path)):
        return None
    try:
        with open(marker_path) as handle:
            marker = json.load(handle)
        if marker.get("input_signature") != input_signature:
            return None
        pieces: List[Piece] = []
        with IndexedFasta(fasta_path) as sequences, open(state_path) as state:
            available = set(sequences.names())
            for raw in state:
                if not raw.strip():
                    continue
                fields = json.loads(raw)
                # Compatibility with the short-lived rank-persisting format.
                fields.pop("sample_rank", None)
                temp_id = fields["temp_id"]
                if temp_id not in available:
                    return None
                fields["seq"] = sequences.sequence(temp_id)
                pieces.append(Piece(**fields))
        if len(pieces) != int(marker.get("piece_count", -1)):
            return None
        if _pieces_signature(pieces) != marker.get("output_signature"):
            return None
        return pieces, int(marker.get("removed_bp", 0)), bool(marker.get("converged", False))
    except (OSError, ValueError, TypeError, KeyError, json.JSONDecodeError):
        return None


_ASSEMBLY_RANK_RE = re.compile(r"^(?:mainlift|locallift)_([0-9]+)_")


def _assembly_rank(piece: Piece) -> int | None:
    match = _ASSEMBLY_RANK_RE.match(piece.lift_id)
    return int(match.group(1)) if match is not None else None


def candidate_rank(
    pieces: Sequence[Piece], sample_priorities: Dict[str, int] | None = None,
) -> Dict[str, int]:
    """Deterministic keep priority: configured sample bonus plus sequence length."""
    contig_order: Dict[str, int] = {}
    for p in pieces:
        if p.contig not in contig_order:
            contig_order[p.contig] = len(contig_order)
    score_ranked = sample_priorities is not None
    ranked = any(_assembly_rank(piece) is not None for piece in pieces)
    if score_ranked:
        assert sample_priorities is not None
        order = sorted(pieces, key=lambda p: (
            -(sample_priorities.get(p.sample, 0) + p.length),
            p.within_sample_rank if p.within_sample_rank >= 0 else sys.maxsize,
            p.sample, contig_order[p.contig], p.start, p.end, p.temp_id,
        ))
    elif ranked:
        fallback = max(
            (_assembly_rank(piece) or 0 for piece in pieces),
            default=0,
        ) + 1
        order = sorted(pieces, key=lambda p: (
            _assembly_rank(p) if _assembly_rank(p) is not None else fallback,
            contig_order[p.contig], p.start, p.end, p.temp_id,
        ))
    else:
        # Preserve the historical behavior for first-reference/legacy pieces
        # that predate assembly-ranked lift identifiers.
        order = sorted(
            pieces, key=lambda p: (contig_order[p.contig], p.start, p.end, p.temp_id),
        )
    return {p.temp_id: i for i, p in enumerate(order)}


def compute_redundant_coverage(
    pieces: Sequence[Piece], hits: Sequence["object"],
    provenance: Dict[str, set] | None = None,
    sample_priorities: Dict[str, int] | None = None,
) -> Dict[str, List[Interval]]:
    """Redundant candidate-local intervals per candidate (keep-earlier).

    Aligner output direction is not stable: a homologous pair may be emitted as
    earlier-query -> later-target, later-query -> earlier-target, or both.  Always
    normalize the hit to the earlier/later occurrence before assigning redundant
    coverage so the biological result is independent of that reporting choice.
    """
    rank = candidate_rank(pieces, sample_priorities)
    by_id = {p.temp_id: p for p in pieces}
    cov: Dict[str, List[Interval]] = cl.defaultdict(list)
    for hit in hits:
        q, t = hit.query_id, hit.target_id
        if q not in by_id or t not in by_id:
            continue
        qblocks = merge_intervals(hit.query_intervals)
        tblocks = merge_intervals(hit.target_intervals)
        if not qblocks or not tblocks:
            continue
        if q == t:
            q0 = min(a for a, _ in qblocks)
            t0 = min(a for a, _ in tblocks)
            if q0 == t0:  # diagonal/self-equivalent interval: exclude nothing
                continue
            if q0 < t0:
                winner_id, loser_id = q, t
                winner_blocks, loser_blocks = qblocks, tblocks
            else:
                winner_id, loser_id = t, q
                winner_blocks, loser_blocks = tblocks, qblocks
            redundant = subtract_intervals(
                loser_blocks, winner_blocks,
            )  # later interval - (earlier intersect later)
        else:
            if rank[q] < rank[t]:
                winner_id, loser_id = q, t
                winner_blocks, loser_blocks = qblocks, tblocks
            else:
                winner_id, loser_id = t, q
                winner_blocks, loser_blocks = tblocks, qblocks
            # Different candidate coordinate systems do not overlap, so the
            # aligned interval on the later occurrence is wholly redundant.
            redundant = loser_blocks
        if redundant:
            cov[loser_id].extend(redundant)
            if provenance is not None:
                wa = min((a for a, _ in winner_blocks), default=0)
                wb = max((b for _, b in winner_blocks), default=by_id[winner_id].length)
                provenance.setdefault(by_id[loser_id].origin_id, set()).add((
                    by_id[winner_id].origin_id,
                    by_id[winner_id].start + wa,
                    by_id[winner_id].start + wb,
                    getattr(hit, "strand", "+"),
                ))
    return {tid: merge_intervals(v) for tid, v in cov.items()}


def _keep_ok(
    sub: str, min_segment: int, mode: str, masked_weight: float = 0.0,
) -> bool:
    if mode == "first_genome":
        return count_unmasked(sub) >= min_segment or count_masked(sub) >= min_segment
    return novelty_score(sub, masked_weight) >= min_segment


def keep_pieces(pieces: Sequence[Piece], coverage: Dict[str, List[Interval]],
                min_segment: int, mode: str,
                masked_weight: float = 0.0) -> Tuple[List[Piece], int]:
    """Return (kept_pieces, removed_bp_this_round)."""
    kept: List[Piece] = []
    removed_bp = 0
    for p in pieces:
        cov = merge_intervals(coverage.get(p.temp_id, []))
        if mode == "first_genome":
            redundant, keep_intervals = find_redundant_and_kept(cov, p.length, p.seq, min_segment)
        else:
            redundant = cov
            # Apply the same novelty rule used against the main/pool: preserve
            # unaligned sequence, bridge across represented islands containing
            # < min_segment unmasked bases, then retain valid residual segments.
            keep_intervals = find_unaligned_segments(
                cov, p.length, p.seq, min_segment, masked_weight,
            )
        subidx = 0
        for a, b in keep_intervals:
            for na, nb, nfree in iter_n_free_segments(p.seq[a:b]):
                if not _keep_ok(nfree, min_segment, mode, masked_weight):
                    continue
                subidx += 1
                s0, s1 = a + na, a + nb
                kept.append(Piece(
                    temp_id=f"{p.temp_id}_nr{subidx:04d}", sample=p.sample, contig=p.contig,
                    start=p.start + s0, end=p.start + s1, strand=p.strand, seq=nfree,
                    parent_id=p.temp_id, note="nonredundant_piece", origin_id=p.origin_id,
                    lift_id=p.lift_id,
                    within_sample_rank=p.within_sample_rank,
                ))
    removed_bp = sum(p.length for p in pieces) - sum(p.length for p in kept)
    return kept, max(0, removed_bp)


def self_clean_once(
    pieces: Sequence[Piece], hits: Sequence["object"], min_segment: int,
    mode: str, sample_priorities: Dict[str, int] | None = None,
    masked_weight: float = 0.0,
) -> Tuple[List[Piece], int]:
    coverage = compute_redundant_coverage(
        pieces, hits, sample_priorities=sample_priorities,
    )
    if masked_weight:
        return keep_pieces(pieces, coverage, min_segment, mode, masked_weight)
    return keep_pieces(pieces, coverage, min_segment, mode)


def self_clean_converge(pieces: Sequence[Piece], align_fn: AlignFn, min_segment: int,
                        mode: str, log: logging.Logger, max_cycles: Optional[int] = 10,
                        provenance: Dict[str, set] | None = None,
                        min_cycles: int = 1,
                        warn_on_limit: bool = True,
                        checkpoint_dir: str | None = None,
                        resume: bool = False,
                        sample_priorities: Dict[str, int] | None = None,
                        convergence_bp: int = SELF_CLEAN_CONVERGENCE_BP,
                        masked_weight: float = 0.0) -> List[Piece]:
    """Repeat self-align until convergence or the configured cycle cap.

    A cycle is also converged when it preserves the piece count and trims less
    than ``convergence_bp`` in total. Re-aligning after such a small edge trim
    cannot by itself contribute another valid output segment.
    """
    if convergence_bp <= 0:
        raise ValueError("convergence_bp must be positive")
    current: List[Piece] = list(pieces)
    converged = False
    cycle = 0
    while max_cycles is None or cycle < max_cycles:
        if not current:
            converged = True
            break
        # Re-base to short, unique ids each round so temp_ids don't accumulate
        # _nr suffixes across cycles and stay safe for winnowmap/BLAST headers.
        # align_fn writes the candidate FASTA from these ids, so it sees them.
        for i, p in enumerate(current):
            p.temp_id = f"cand{i:08d}"
        input_signature = _checkpoint_input_signature(current, align_fn)
        cycle_checkpoint_dir = (
            os.path.join(checkpoint_dir, f"cycle{cycle:02d}")
            if checkpoint_dir is not None else None
        )
        if resume and cycle_checkpoint_dir is not None:
            input_piece_count = len(current)
            checkpoint = _read_cycle_checkpoint(cycle_checkpoint_dir, input_signature)
            if checkpoint is not None:
                current, removed_bp, cycle_converged = checkpoint
                # Older checkpoints recorded convergence only at zero removal.
                # Re-evaluate the safe sub-threshold condition while loading so
                # an existing 604 -> 604, 47-bp cycle can stop immediately.
                cycle_converged = cycle_converged or (
                    len(current) == input_piece_count
                    and removed_bp < convergence_bp
                )
                log.info(
                    "RESUME self-clean [%s] cycle %d: loaded %d pieces "
                    "(removed %d bp)",
                    mode, cycle, len(current), removed_bp,
                )
                if cycle_converged and cycle + 1 >= min_cycles:
                    converged = True
                    break
                cycle += 1
                continue
        hits = align_fn(current, cycle)
        coverage = compute_redundant_coverage(
            current, hits, provenance, sample_priorities=sample_priorities,
        )
        if masked_weight:
            kept, removed_bp = keep_pieces(
                current, coverage, min_segment, mode, masked_weight,
            )
        else:
            kept, removed_bp = keep_pieces(current, coverage, min_segment, mode)
        log.info("Self-clean [%s] cycle %d: %d -> %d pieces, removed %d bp",
                 mode, cycle, len(current), len(kept), removed_bp)
        same_piece_count = len(kept) == len(current)
        current = kept
        cycle_converged = (
            (removed_bp == 0 or (
                same_piece_count and removed_bp < convergence_bp
            ))
            and cycle + 1 >= min_cycles
        )
        if cycle_checkpoint_dir is not None:
            _write_cycle_checkpoint(
                cycle_checkpoint_dir, current, input_signature,
                removed_bp, cycle_converged,
            )
        if cycle_converged:
            converged = True
            break
        cycle += 1
    if current and not converged and warn_on_limit and max_cycles is not None:
        log.warning("Self-clean [%s] reached max_cycles=%d before convergence",
                    mode, max_cycles)
    return current
