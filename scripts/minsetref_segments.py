#!/usr/bin/env python3
"""Segment/interval logic for the min-set reference builder (v2 plan).

This module isolates the new merge rule so it can be unit-tested without any
aligner.  Two contexts use the *same* bridging primitive with opposite polarity:

  * Assembly novelty: we keep UNALIGNED segments.  Two unaligned segments are
    merged across an interior ALIGNED island iff that island carries fewer than
    ``min_segment`` UNMASKED (uppercase A/C/G/T) bases.  A valid novel segment
    then needs >= ``min_segment`` unmasked bases.

  * First-genome self-redundancy: we keep the REDUNDANT (aligned A-(A&B))
    segments to remove.  Two redundant segments are merged across an interior
    UNALIGNED island iff that island carries fewer than ``min_segment`` MASKED
    (lowercase a/c/g/t) bases.  The kept sequence is the complement.

The asymmetry (unmasked for novelty islands, masked for redundancy islands) is
intentional and matches the confirmed spec.
"""
from __future__ import annotations

from typing import List, Sequence, Tuple

Interval = Tuple[int, int]


# ---------------------------------------------------------------------------
# C-level base counters
# ---------------------------------------------------------------------------

def count_unmasked(seq: str) -> int:
    """Uppercase A/C/G/T count (C-speed)."""
    return seq.count("A") + seq.count("C") + seq.count("G") + seq.count("T")


def count_masked(seq: str) -> int:
    """Lowercase a/c/g/t count (soft-masked bases), C-speed."""
    return seq.count("a") + seq.count("c") + seq.count("g") + seq.count("t")


def novelty_score(seq: str, masked_weight: float = 0.0) -> float:
    """Uppercase bases plus the configured fraction of lowercase bases."""
    return count_unmasked(seq) + float(masked_weight) * count_masked(seq)


def significant_bp(
    seq: str, a: int, b: int, kind: str, masked_weight: float = 0.0,
) -> float:
    """Count unmasked or masked bases in seq[a:b]."""
    sub = seq[a:b]
    if kind == "unmasked":
        return novelty_score(sub, masked_weight)
    if kind == "masked":
        return count_masked(sub)
    raise ValueError(f"kind must be 'unmasked' or 'masked', got {kind!r}")


# ---------------------------------------------------------------------------
# Interval primitives (0-based, half-open)
# ---------------------------------------------------------------------------

def merge_intervals(intervals: Sequence[Interval]) -> List[Interval]:
    vals = sorted((min(a, b), max(a, b)) for a, b in intervals if max(a, b) > min(a, b))
    if not vals:
        return []
    out = [vals[0]]
    for a, b in vals[1:]:
        la, lb = out[-1]
        if a <= lb:
            out[-1] = (la, max(lb, b))
        else:
            out.append((a, b))
    return out


def complement_intervals(intervals: Sequence[Interval], start: int, end: int) -> List[Interval]:
    merged = merge_intervals(intervals)
    out: List[Interval] = []
    cur = start
    for a, b in merged:
        a = max(start, min(a, end))
        b = max(start, min(b, end))
        if a > cur:
            out.append((cur, a))
        cur = max(cur, b)
    if cur < end:
        out.append((cur, end))
    return out


def subtract_intervals(intervals: Sequence[Interval], masks: Sequence[Interval]) -> List[Interval]:
    """Intervals minus masks by a two-pointer scan (half-open)."""
    source = merge_intervals(intervals)
    masks_merged = merge_intervals(masks)
    if not source or not masks_merged:
        return source
    out: List[Interval] = []
    mask_index = 0
    for a, b in source:
        while mask_index < len(masks_merged) and masks_merged[mask_index][1] <= a:
            mask_index += 1
        cursor = a
        scan = mask_index
        while scan < len(masks_merged) and masks_merged[scan][0] < b:
            ma, mb = masks_merged[scan]
            if ma > cursor:
                out.append((cursor, min(ma, b)))
            cursor = max(cursor, mb)
            if cursor >= b:
                break
            scan += 1
        if cursor < b:
            out.append((cursor, b))
        # Preserve a mask that extends beyond this source interval: it may also
        # overlap the next source interval. Fully consumed masks are skipped by
        # the leading while loop on the next iteration.
        mask_index = scan
    return out


# ---------------------------------------------------------------------------
# Bridging: merge same-class segments across small interior islands
# ---------------------------------------------------------------------------

def bridge_small_islands(
    target_intervals: Sequence[Interval],
    seq: str,
    min_content: int,
    content_kind: str,
    masked_weight: float = 0.0,
) -> List[Interval]:
    """Merge consecutive target segments across the interior gap between them
    when that gap ("island") carries fewer than ``min_content`` bases of the
    given kind.

    Only interior gaps (flanked by a target segment on both sides) are ever
    bridged; terminal caps before the first / after the last target segment are
    left untouched, so a small aligned cap is never reclassified as novel and a
    small unaligned cap is never reclassified as redundant.
    """
    target = merge_intervals(target_intervals)
    if len(target) <= 1:
        return target
    out: List[Interval] = [target[0]]
    for cur_a, cur_b in target[1:]:
        prev_a, prev_b = out[-1]
        # island is the interior gap between the previous kept segment and this one
        if significant_bp(
            seq, prev_b, cur_a, content_kind, masked_weight,
        ) < min_content:
            out[-1] = (prev_a, cur_b)  # bridge across the small island
        else:
            out.append((cur_a, cur_b))
    return out


# ---------------------------------------------------------------------------
# Context wrappers
# ---------------------------------------------------------------------------

def find_unaligned_segments(
    aligned_intervals: Sequence[Interval],
    contig_len: int,
    seq: str,
    min_segment: int,
    masked_weight: float = 0.0,
) -> List[Interval]:
    """Assembly-novelty segments.

    ``aligned_intervals`` is the union of valid alignment query-intervals.  The
    unaligned complement is bridged across interior aligned islands whose
    novelty score is < ``min_segment``, then only segments with novelty score
    >= ``min_segment`` are returned. Novelty score is uppercase A/C/G/T plus
    ``masked_weight`` times lowercase a/c/g/t.
    """
    aligned = merge_intervals(aligned_intervals)
    unaligned = complement_intervals(aligned, 0, contig_len)
    bridged = bridge_small_islands(
        unaligned, seq, min_segment, "unmasked", masked_weight,
    )
    return [
        (a, b)
        for a, b in bridged
        if novelty_score(seq[a:b], masked_weight) >= min_segment
    ]


def find_redundant_and_kept(
    redundant_intervals: Sequence[Interval],
    contig_len: int,
    seq: str,
    min_segment: int,
) -> Tuple[List[Interval], List[Interval]]:
    """First-genome self-redundancy.

    ``redundant_intervals`` is the A-(A&B) redundant coverage (already restricted
    to alignments with >= min_segment unmasked upstream).  Redundant segments are
    bridged across interior unaligned islands whose *masked* content is
    < ``min_segment``.  Returns (merged_redundant, kept) where kept is the
    complement over the contig.
    """
    redundant = merge_intervals(redundant_intervals)
    bridged = bridge_small_islands(redundant, seq, min_segment, "masked")
    kept = complement_intervals(bridged, 0, contig_len)
    return bridged, kept


if __name__ == "__main__":  # lightweight self-test / demonstration
    import sys

    def show(title, result):
        print(f"{title}: {result}")

    # ---- Novelty demo ---------------------------------------------------
    # contig: [0,300) unaligned(unmasked) | [300,350) aligned small(unmasked)
    #         [350,700) unaligned(unmasked) | [700,1100) aligned big(unmasked)
    #         [1100,1450) unaligned
    seq = ("A" * 300) + ("A" * 50) + ("A" * 350) + ("A" * 400) + ("A" * 350)
    aligned = [(300, 350), (700, 1100)]
    novel = find_unaligned_segments(aligned, len(seq), seq, 300)
    # 50bp aligned island (<300 unmasked) bridges 0-300 with 350-700 -> 0-700;
    # 400bp aligned island (>=300) stays; 1100-1450 separate.
    assert novel == [(0, 700), (1100, 1450)], novel
    show("novelty bridged (small aligned island absorbed)", novel)

    # small aligned CAP at the end must NOT be bridged into novelty
    seq2 = ("A" * 400) + ("A" * 50)  # 400 unaligned + 50 aligned cap
    novel2 = find_unaligned_segments([(400, 450)], len(seq2), seq2, 300)
    assert novel2 == [(0, 400)], novel2
    show("novelty terminal aligned cap left alone", novel2)

    # a novel segment below 300 unmasked is dropped
    seq3 = ("A" * 100) + ("C" * 800)  # 100 unaligned(<300) then aligned
    novel3 = find_unaligned_segments([(100, 900)], len(seq3), seq3, 300)
    assert novel3 == [], novel3
    show("novelty sub-threshold segment dropped", novel3)

    # ---- Redundancy demo (first genome; island polarity = MASKED) -------
    # redundant [0,400) | unaligned island [400,450) with few masked -> bridge
    # redundant [450,900) | unaligned island [900,1300) masked-rich -> keep
    seqR = ("A" * 400) + ("A" * 50) + ("A" * 450) + ("a" * 400) + ("A" * 300)
    redundant_in = [(0, 400), (450, 900)]
    red, kept = find_redundant_and_kept(redundant_in, len(seqR), seqR, 300)
    # island [400,450) has 0 masked (<300) -> bridge to [0,900)
    assert red == [(0, 900)], red
    assert kept == [(900, 1600)], kept
    show("redundancy bridged across low-masked island", red)
    show("redundancy kept complement", kept)

    # unaligned island that is masked-rich (>=300 masked) blocks the bridge
    seqR2 = ("A" * 400) + ("a" * 350) + ("A" * 450)
    red2, kept2 = find_redundant_and_kept([(0, 400), (750, 1200)], len(seqR2), seqR2, 300)
    assert red2 == [(0, 400), (750, 1200)], red2
    assert kept2 == [(400, 750)], kept2
    show("redundancy masked-rich island preserved as kept", kept2)

    print("all segment self-tests passed", file=sys.stderr)
