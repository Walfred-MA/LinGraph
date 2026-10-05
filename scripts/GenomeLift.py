#!/usr/bin/env python3
import shutil
import argparse
import bisect
import collections as cl
import gzip
import math
import os
import re
import sys
import tempfile
from dataclasses import dataclass
from multiprocessing import Pool
from typing import Dict, List, Tuple

from assembly_contigs import HG38_MAIN_CONTIGS, accepted_contig

CANONICAL_HG38 = HG38_MAIN_CONTIGS

GLOBAL_REF_ALLELES = None
GLOBAL_PRIORITY_BY_PREFIX = None
GLOBAL_REF_ANCHOR_INDEX = None
GLOBAL_ASSIGNABLE_REF_ALLELES = None
GLOBAL_REF_TAG_OPTIONS = None
GLOBAL_NON_REFERENCE_TAGS = frozenset()
GLOBAL_ALIGNMENT_SCORES = None
CONTIG_MERGE_FAN_IN = 32
DATACLASS_OPTIONS = {"slots": True} if sys.version_info >= (3, 10) else {}


def init_worker(
    ref_alleles,
    priority_by_prefix=None,
    non_reference_tags=frozenset(),
    alignment_scores=None,
):
    global GLOBAL_REF_ALLELES, GLOBAL_PRIORITY_BY_PREFIX
    global GLOBAL_REF_ANCHOR_INDEX, GLOBAL_ASSIGNABLE_REF_ALLELES
    global GLOBAL_REF_TAG_OPTIONS, GLOBAL_NON_REFERENCE_TAGS
    global GLOBAL_ALIGNMENT_SCORES
    GLOBAL_REF_ALLELES = ref_alleles
    GLOBAL_PRIORITY_BY_PREFIX = priority_by_prefix or {}
    GLOBAL_NON_REFERENCE_TAGS = non_reference_tags
    # Keep ``None`` distinct from an explicitly loaded (possibly empty)
    # alignment file.  The former preserves the historical ownership order;
    # the latter enables the same-graph score/span tie-breaker.
    GLOBAL_ALIGNMENT_SCORES = alignment_scores
    GLOBAL_REF_ANCHOR_INDEX = build_reference_anchor_index(ref_alleles)
    GLOBAL_ASSIGNABLE_REF_ALLELES = eligible_ref_alleles(ref_alleles)
    GLOBAL_REF_TAG_OPTIONS = {
        False: reference_tag_pair_options(ref_alleles, ignore_distance=False),
        True: reference_tag_pair_options(ref_alleles, ignore_distance=True),
    }


def open_maybe_gzip(path):
    path = str(path)
    if path.endswith(".gz"):
        return gzip.open(path, "rt")
    return open(path, "r")


_ALIGNMENT_OPERATION_RE = re.compile(r"(\d+)([=MXIDHS])([A-Za-z]*)")
_ALIGNMENT_HEADER = (
    "hotspot_index", "contig", "start", "end", "strand", "queryname",
    "graph_path", "graphcigar", "refpositions", "qpositions",
)
_ALIGNMENT_COMPONENT_RE = re.compile(r"([<>])([^<>]+)")
_ALIGNMENT_SEGMENT_RE = re.compile(r"([<>])([^:<>]+):([^<>]*)")
# Final lift stages (see _elect_final_stages).
FINAL_STAGE_HIGH_COVERAGE = 0.9
FINAL_STAGE_HIGH_KMER_SIMILARITY = 0.9
FINAL_STAGE_MIN_COVERAGE = 0.5
FINAL_STAGE_MIN_KMER_SIMILARITY = 0.5
FINAL_STAGE_ANCHOR_BONUS = 0.2
FINAL_STAGE_PROMOTION_MARGIN = 1.5
FINAL_STAGE_MULTI_MIN_SHARE = 0.1
FINAL_STAGE_MULTI_MIN_REFERENCE_COVERAGE = 0.5
FINAL_STAGE_MULTI_MAX_SHARED = 0.5
FINAL_STAGE_MIN_PATH_COVERAGE = 2_000
FINAL_STAGE_COLUMNS = (
    "assignment_stage", "assignment_tier", "location_mapping",
    "sequence_mapping",
)


def graph_alignment_affine_score(graph_cigar: str) -> int:
    """Score every aligned graph segment in one ``_align.txt`` row.

    Unlike a linear reference CIGAR, every named segment here belongs to the
    same graph traversal; later named segments are not secondary insertions.
    The score is therefore accumulated over all segments while allowing a gap
    operation to continue across a segment boundary:

        matches - 4*mismatches - 4*gap_openings - gap_extensions
    """
    matches = 0
    mismatches = 0
    gap_openings = 0
    gap_extensions = 0
    active_gap = ""
    found_segment = False
    for segment_match in re.finditer(r"[<>][^:<>]+:([^<>]*)", graph_cigar or ""):
        found_segment = True
        body = segment_match.group(1)
        for operation_match in _ALIGNMENT_OPERATION_RE.finditer(body):
            length = int(operation_match.group(1))
            operation = operation_match.group(2)
            if length <= 0:
                continue
            if operation in {"I", "D"}:
                if active_gap == operation:
                    gap_extensions += length
                else:
                    gap_openings += 1
                    gap_extensions += max(0, length - 1)
                    active_gap = operation
                continue
            if operation == "H":
                continue
            active_gap = ""
            if operation in {"=", "M"}:
                matches += length
            elif operation == "X":
                mismatches += length
    if not found_segment:
        return 0
    return (
        matches
        - 4 * mismatches
        - 4 * gap_openings
        - gap_extensions
    )


def canonical_alignment_allele_name(name: str) -> str:
    """Normalize only zero padding in the terminal PA numeric index."""
    try:
        group, sample, haplotype, index, _genome = parse_allele_name(name)
    except (TypeError, ValueError):
        return name
    match = re.fullmatch(r"0*(\d+)(.*)", index)
    if match is None:
        return name
    return (
        f"{group}_{sample}_{haplotype}_{int(match.group(1))}"
        f"{match.group(2)}"
    )


def graph_alignment_affine_score_in_query(
    graph_cigar: str, query_positions, query_start: int, query_end: int,
) -> int:
    """Score only the part of an ``_align.txt`` row inside a query window.

    ``query_positions`` is the row's column 10: the oriented query interval
    of every graph segment, in segment order.  Query-consuming operations are
    clipped to ``[query_start, query_end)``; a deletion counts when it lies
    strictly inside the window.  For the complete query this equals
    ``graph_alignment_affine_score``.
    """
    matches = mismatches = gap_openings = gap_extensions = 0
    active_gap = ""
    bodies = [
        match.group(1)
        for match in re.finditer(r"[<>][^:<>]+:([^<>]*)", graph_cigar or "")
    ]
    if not bodies:
        return 0
    if len(bodies) != len(query_positions):
        raise ValueError(
            "graph CIGAR segment and query-position counts differ: "
            f"{len(bodies)}/{len(query_positions)}"
        )
    for body, (segment_query_start, _segment_query_end) in zip(
        bodies, query_positions,
    ):
        cursor = segment_query_start
        for operation_match in _ALIGNMENT_OPERATION_RE.finditer(body):
            length = int(operation_match.group(1))
            operation = operation_match.group(2)
            if length <= 0 or operation == "H":
                continue
            if operation == "D":
                if query_start < cursor < query_end:
                    if active_gap == "D":
                        gap_extensions += length
                    else:
                        gap_openings += 1
                        gap_extensions += max(0, length - 1)
                        active_gap = "D"
                continue
            if operation not in {"=", "M", "X", "I"}:
                continue
            inside = (
                min(cursor + length, query_end) - max(cursor, query_start)
            )
            cursor += length
            if inside <= 0:
                continue
            if operation == "I":
                if active_gap == "I":
                    gap_extensions += inside
                else:
                    gap_openings += 1
                    gap_extensions += max(0, inside - 1)
                    active_gap = "I"
                continue
            active_gap = ""
            if operation == "X":
                mismatches += inside
            else:
                matches += inside
    return (
        matches
        - 4 * mismatches
        - 4 * gap_openings
        - gap_extensions
    )


def read_pa_loci(path: str) -> dict:
    """Index RefMatch/GenomeLift-input PA loci by (graph, assembly contig)."""
    loci = cl.defaultdict(list)
    for parts in iter_rows(path):
        if parts[0] == "allelename":
            continue
        try:
            contig, start, end, strand = parse_loc(parts[2])
        except ValueError:
            continue
        loci[(allele_matrix_name(parts[0]), contig)].append(
            (start, end, strand or "+", parts[0])
        )
    return {
        key: (tuple(value[0] for value in values), tuple(values))
        for key, values in (
            (key, sorted(values)) for key, values in loci.items()
        )
    }


def read_alignment_scores(path: str, pa_loci) -> Dict[str, int]:
    """Score every PA from the ``_align.txt`` row that contains its locus.

    `_align.txt` query names use hotspot-local indices, not GenomeLift PA
    names, so each PA is matched by locus exactly as graphcigartoref_persample
    ``locate_alignment_row`` does: same graph and contig, an alignment
    interval containing the PA, preferring the same strand, then the smallest
    span.  The PA is scored on its own query sub-interval of that row, so PAs
    sharing one hotspot alignment are ranked by their own alignment quality.
    Unaligned (``*``) rows give no score.
    """
    if not path or not pa_loci:
        return {}
    best = {}
    alignment_header = "\t".join(_ALIGNMENT_HEADER)
    with open(path, "rt") as handle:
        for line_number, raw in enumerate(handle, 1):
            if not raw.strip() or raw.startswith("#"):
                continue
            fields = raw.rstrip("\r\n").split("\t")
            if tuple(fields) == _ALIGNMENT_HEADER:
                continue
            if len(fields) != len(_ALIGNMENT_HEADER):
                raise ValueError(
                    f"{path}:{line_number}: expected ten alignment columns, "
                    f"found {len(fields)}"
                )
            indexed = pa_loci.get((allele_matrix_name(fields[5].strip()),
                                   fields[1].strip()))
            if indexed is None:
                continue
            try:
                row_start, row_end = int(fields[2]), int(fields[3])
            except ValueError as error:
                raise ValueError(
                    f"{path}:{line_number}: invalid alignment interval"
                ) from error
            row_strand = fields[4].strip() or "+"
            aligned = fields[6].strip() not in {"", "*"} and fields[7].strip() not in {"", "*"}
            starts, loci = indexed
            query_positions = None
            index = bisect.bisect_left(starts, row_start)
            while index < len(loci) and loci[index][0] < row_end:
                start, end, strand, name = loci[index]
                index += 1
                if end > row_end:
                    continue
                key = (
                    0 if row_strand == strand else 1,
                    row_end - row_start,
                    row_start,
                    row_end,
                    0 if aligned else 1,
                    line_number,
                )
                previous = best.get(name)
                if previous is not None and previous[0] <= key:
                    continue
                score = None
                if aligned:
                    if query_positions is None:
                        query_positions = _parse_alignment_position_list(
                            fields[9].strip(), f"{path}:{line_number}",
                        )
                    if row_strand == "-":
                        local = (row_end - end, row_end - start)
                    else:
                        local = (start - row_start, end - row_start)
                    score = graph_alignment_affine_score_in_query(
                        fields[7].strip(), query_positions, *local,
                    )
                best[name] = (key, score)
    return {
        name: score
        for name, (_key, score) in best.items()
        if score is not None
    }


def alignment_score_for_allele(
    allelename: str, alignment_scores=None,
):
    if alignment_scores is None:
        alignment_scores = GLOBAL_ALIGNMENT_SCORES
    if alignment_scores is None:
        return None
    if allelename in alignment_scores:
        return alignment_scores[allelename]
    return alignment_scores.get(canonical_alignment_allele_name(allelename))


def parse_non_reference_tags(value: str):
    """Parse exact non-reference allele names, or the special value ``all``."""
    tokens = tuple(
        token.strip()
        for token in (value or "").split(",")
        if token.strip()
    )
    if not tokens:
        return frozenset()
    lowered = {token.lower() for token in tokens}
    if "all" in lowered:
        if len(tokens) != 1:
            raise ValueError("--non-reference-tags all cannot be combined with allele names")
        return "all"
    return frozenset(tokens)


def non_reference_tag_enabled(allele: str, ref_alleles, selection) -> bool:
    """Return whether an otherwise unfixed non-reference allele may self-seed."""
    if allele in ref_alleles:
        return False
    return selection == "all" or allele in selection


def overlaploci(regions):
    if not regions:
        return []

    events = []
    sizes = []
    for idx, (start, end, *_) in enumerate(regions):
        if end < start:
            start, end = end, start
        events.append((start, 0, idx))
        events.append((end, 1, idx))
        sizes.append(end - start)

    events.sort(key=lambda e: (e[0], e[1]))

    active = set()
    best_idx = None
    result = []

    for coord, etype, idx in events:
        if etype == 0:
            active.add(idx)
            if best_idx is None or sizes[idx] > sizes[best_idx]:
                if best_idx is not None and result:
                    result[-1][1] = coord
                result.append([coord, regions[idx][1], idx])
                best_idx = idx
        else:
            if idx in active:
                active.remove(idx)
            if best_idx == idx:
                if active:
                    new_best = max(active, key=lambda j: sizes[j])
                    if new_best != best_idx:
                        if result:
                            result[-1][1] = coord
                        result.append([coord, regions[new_best][1], new_best])
                        best_idx = new_best
                    else:
                        if result:
                            result[-1][1] = coord
                else:
                    if result:
                        result[-1][1] = coord
                    best_idx = None
    return result


def priority_score_for_allele(allelename: str, priority_by_prefix=None):
    """Return the compact integer priority rank for an allele prefix.

    Higher values have higher priority. Each matrix prefix is ranked first by
    the sum of its reference-supported sequence scores:

        reference size / max(RefMatch distance, 1e-6)

    A matrix without reference support uses its total size across samples.
    Total size and the prefix are deterministic tie-breakers. Unknown prefixes
    get priority 0.
    """
    if priority_by_prefix is None:
        priority_by_prefix = GLOBAL_PRIORITY_BY_PREFIX or {}
    return priority_by_prefix.get(allele_matrix_name(allelename), 0)


def collapse_unique_regions(
    regions, priority_by_prefix=None, priority_by_allele=None,
):
    """
    Compute the longest priority-owned segment for each allele.

    Overlaps are assigned first by matrix-prefix priority. Within one known
    matrix/graph prefix only, optional complete-graph alignment scores rank
    competing PAs. Thus alignment quality cannot flip ownership between
    different graphs. A missing or tied alignment score falls back to the
    larger PA span and then the historical stable row order.

    Regions are stored compactly as (start, end, score, original_index) for the
    sweep, keeping memory use low. Coordinates are 0-based, right-open.
    """
    if not regions:
        return {}

    compact = []
    allcoords = []
    norm_regions = []

    if priority_by_prefix is None:
        priority_by_prefix = GLOBAL_PRIORITY_BY_PREFIX or {}
    alignment_ranking_enabled = priority_by_allele is not None
    priority_by_allele = priority_by_allele or {}
    prefix_counts = cl.Counter(
        allele_matrix_name(name) for _start, _end, name in regions
    )
    scored_prefixes = {
        prefix
        for prefix, count in prefix_counts.items()
        if (
            alignment_ranking_enabled
            and prefix in priority_by_prefix
            and count > 1
        )
    }

    for idx, (start, end, name) in enumerate(regions):
        if end < start:
            start, end = end, start
        prefix = allele_matrix_name(name)
        prefix_score = priority_score_for_allele(name, priority_by_prefix)
        alignment_score = alignment_score_for_allele(
            name, priority_by_allele,
        )
        if prefix in scored_prefixes:
            score = (
                prefix_score,
                1,
                int(alignment_score is not None),
                int(alignment_score or 0),
                end - start,
            )
        else:
            # This is exactly the old prefix-only ordering. ``region_i`` below
            # remains the final stable tie-breaker.
            score = (prefix_score, 0, 0, 0, 0)
        compact.append((start, end, score, idx))
        norm_regions.append((start, end, name))
        allcoords.extend((start, end))

    coord_order = sorted(range(len(allcoords)), key=lambda x: (allcoords[x], x % 2, x))

    active = set()  # (score, compact_index); max(active) is the current owner
    highest = ((-math.inf, 0, 0, 0, 0), -1)
    curr_start = 0
    result = []

    for coord_index in coord_order:
        coord = allcoords[coord_index]
        region_i = coord_index // 2
        start, end, score, original_idx = compact[region_i]
        active_key = (score, region_i)

        if coord_index % 2 == 0:  # start event
            if not active:
                curr_start = coord

            active.add(active_key)

            if active_key > highest:
                if highest[1] >= 0:
                    _, high_region_i = highest
                    _, _, _, high_original_idx = compact[high_region_i]
                    result.append((curr_start, coord, high_original_idx))
                highest = active_key
                curr_start = coord

        else:  # end event
            if highest[1] >= 0:
                _, high_region_i = highest
                _, _, _, high_original_idx = compact[high_region_i]
                result.append((curr_start, coord, high_original_idx))

            active.discard(active_key)
            highest = (
                max(active)
                if active else ((-math.inf, 0, 0, 0, 0), -1)
            )
            curr_start = coord

    result = [x for x in result if x[1] - x[0] > 0]

    seg_by_idx = [[] for _ in regions]
    for sub_start, sub_end, idx in result:
        seg_by_idx[idx].append((sub_start, sub_end))

    # The sweep emits a boundary whenever *any* active interval starts or
    # ends.  When the winning allele does not change, this can split one
    # continuous owned interval into several touching pieces.  Coalesce those
    # pieces before applying the one-component output rule; otherwise choosing
    # only the largest piece invents gaps and corrupts columns 13/14.
    for idx, pieces in enumerate(seg_by_idx):
        merged = []
        for sub_start, sub_end in sorted(pieces):
            if merged and sub_start <= merged[-1][1]:
                merged[-1] = (
                    merged[-1][0], max(merged[-1][1], sub_end),
                )
            else:
                merged.append((sub_start, sub_end))
        seg_by_idx[idx] = merged

    out = {}
    for i, (start, end, name) in enumerate(norm_regions):
        if not seg_by_idx[i]:
            out[name] = (-1, -1, 0, -1, -1)
        else:
            # The output model stores one assigned interval. Keep the largest
            # owned fragment if a lower-priority allele is split by overlaps.
            best = max(seg_by_idx[i], key=lambda x: (x[1] - x[0], -x[0], -x[1]))
            sub_start, sub_end = best
            out[name] = (
                sub_start,
                sub_end,
                sub_end - sub_start,
                sub_start - start,
                end - sub_end,
            )
    return out


NEG_INF_STR = "-inf"



def compute_priority_display_offsets_for_regions(regions, unique_map):
    """Compute neighbor-aware ownership tokens for output columns 13/14.

    unique_map is produced by collapse_unique_regions() and records the segment
    owned by each allele after prefix-priority overlap assignment. Its left/right
    offsets are truncation amounts relative to the original allele interval.

    Every side with an adjacent positive-ownership allele is emitted as
    ``exact_neighbor_allele:value``.  ``value`` retains the ownership sign:

    * ``+N``: this allele owns the gap on that side;
    * ``N``: the adjacent allele has first priority over the gap;
    * ``-N``: N overlapping bases are excluded from this allele;
    * ``-inf``: this allele has no usable owned interval.

    Adjacency is determined from the priority-owned intervals, not the original
    overlapping interval starts. This makes the two sides reciprocal even when
    a high-priority interval truncates the left or right side of a larger one.

    For non-overlapping neighboring valid alleles, keep the true 0-based,
    right-open gap distance:
        gap = curr_start - prev_end
    so [16834473,16942797) followed by [16942798,17007856) reports 1, not 0.

    For touching neighbors, report signed zero as before.
    For overlaps, emit a negative truncation amount when an allele side was
    cut by a higher-priority neighbor; use signed zero on the owning side.
    """
    if not regions:
        return {}

    norm_all = []
    for input_index, (start, end, name) in enumerate(regions):
        if end < start:
            start, end = end, start
        norm_all.append((start, end, name, input_index))

    offsets = {
        name: [NEG_INF_STR, NEG_INF_STR]
        for _, _, name, _ in norm_all
    }
    neighbors = {
        name: ["", ""]
        for _, _, name, _ in norm_all
    }
    owned = []

    for start, end, name, input_index in norm_all:
        sub_start, sub_end, unique_len, left_trim, right_trim = unique_map.get(
            name, (-1, -1, 0, -1, -1)
        )
        if unique_len > 0:
            offsets[name] = [
                f"-{left_trim}" if left_trim > 0 else "0",
                f"-{right_trim}" if right_trim > 0 else "0",
            ]
            owned.append((sub_start, sub_end, name, input_index))

    owned.sort(key=lambda value: (
        value[0], value[1], value[2], value[3],
    ))

    def is_zero_offset(x):
        return str(x) in {"0", "+0", "-0"}

    for i in range(1, len(owned)):
        _, prev_end, prev_name, _ = owned[i - 1]
        curr_start, _, curr_name, _ = owned[i]
        delta = curr_start - prev_end

        # Exact query allele names are deliberately recorded on both sides.
        # Downstream cross-window validation requires reciprocal A*B / B*A
        # evidence and cannot recover it reliably from distance alone.
        neighbors[prev_name][1] = curr_name
        neighbors[curr_name][0] = prev_name

        if delta > 0:
            # True gap in 0-based, right-open coordinates.
            if is_zero_offset(offsets[prev_name][1]):
                offsets[prev_name][1] = str(delta)
            if is_zero_offset(offsets[curr_name][0]):
                offsets[curr_name][0] = f"+{delta}"
        elif delta == 0:
            # Directly touching half-open intervals.
            if is_zero_offset(offsets[prev_name][1]):
                offsets[prev_name][1] = "+0"
            if is_zero_offset(offsets[curr_name][0]):
                offsets[curr_name][0] = "-0"
        else:
            # Priority-owned intervals should not overlap. Retain deterministic
            # signed-zero ownership if legacy/corrupt input nevertheless does.
            if is_zero_offset(offsets[prev_name][1]):
                offsets[prev_name][1] = "+0"
            if is_zero_offset(offsets[curr_name][0]):
                offsets[curr_name][0] = "-0"

    def neighbor_token(neighbor, value):
        if not neighbor or value == NEG_INF_STR:
            return value
        return f"{neighbor}:{value}"

    return {
        name: (
            neighbor_token(neighbors[name][0], values[0]),
            neighbor_token(neighbors[name][1], values[1]),
        )
        for name, values in offsets.items()
    }


def compute_priority_offsets_for_regions(regions, unique_map):
    """
    Report display offsets using local pairwise priority.

    Coordinates are 0-based, right-open.

    This function should report the raw local gap/overlap relationship only.
    Any 20 kb capping or lower-priority leftover splitting is handled later by
    build_massive_bed_updated.py.

    For each adjacent pair in sorted order:
    - positive gap d: previous allele gets right=d, current allele gets left=+d
    - touching: previous allele gets right=+0, current allele gets left=-0
    - overlap d<0: previous allele gets right=d, current allele gets left=-0

    For zero offsets, the sign records local priority:
    - +0 means the allele end has higher priority over its adjacent end
    - -0 means the allele end has lower priority than its adjacent end

    Invalid alleles (no truly unique segment) remain -inf / -inf and must not
    participate in deciding which neighboring valid allele owns a gap/overlap.
    """
    if not regions:
        return {}

    norm_all = []
    for start, end, name in regions:
        if end < start:
            start, end = end, start
        norm_all.append((start, end, name))
    norm_all.sort(key=lambda x: (x[0], x[1], x[2]))

    offsets = {name: [NEG_INF_STR, NEG_INF_STR] for _, _, name in norm_all}

    # Only alleles with a positive unique segment are allowed to participate in
    # gap/overlap ownership. Alleles truncated to zero still stay in the output
    # map as -inf/-inf placeholders, but they are invisible to the adjacency
    # calculation below.
    norm_valid = []
    for start, end, name in norm_all:
        if unique_map.get(name, (-1, -1, 0, -1, -1))[2] > 0:
            offsets[name] = ["0", "0"]
            norm_valid.append((start, end, name))

    for i in range(1, len(norm_valid)):
        _, prev_end, prev_name = norm_valid[i - 1]
        curr_start, _, curr_name = norm_valid[i]
        delta = curr_start - prev_end

        if delta > 0:
            offsets[prev_name][1] = str(delta)
            offsets[curr_name][0] = f"+{delta}"
        elif delta == 0:
            offsets[prev_name][1] = "+0"
            offsets[curr_name][0] = "-0"
        else:
            offsets[prev_name][1] = str(delta)
            offsets[curr_name][0] = "-0"

    return {name: (vals[0], vals[1]) for name, vals in offsets.items()}



def parse_loc(s: str):
    s = s.strip()
    match = re.match(r"^(.+):(-?\d+)-(-?\d+)([+-])?$", s)
    if not match:
        raise ValueError(f"Malformed location: {s}")
    contig, start_s, end_s, strand = match.groups()
    start, end = int(start_s), int(end_s)
    strand = strand if strand else "+"
    if end < start:
        start, end = end, start
    return contig, start, end, strand


def parse_allele_name(allelename: str):
    parts = allelename.split("_")
    if len(parts) < 4:
        raise ValueError(f"Unexpected allelename format: {allelename}")

    group_name = "_".join(parts[:-3])
    sample_name = parts[-3]
    hap = parts[-2]
    idx = parts[-1]
    genome = f"{sample_name}_{hap}"
    return group_name, sample_name, hap, idx, genome


def allele_belongs_to_genome(allelename: str, genome: str):
    """Return whether an allele's encoded sample/haplotype is ``genome``."""
    try:
        return parse_allele_name(allelename)[4] == genome
    except (TypeError, ValueError):
        return False


def allele_matrix_name(allelename: str):
    """Return the allele matrix/group prefix before sample, hap, and index."""
    try:
        group_name, _, _, _, _ = parse_allele_name(allelename)
        return group_name
    except Exception:
        return allelename.rsplit("_", 3)[0]


def same_allele_matrix(a: str, b: str):
    return allele_matrix_name(a) == allele_matrix_name(b)


def determine_assembly_genome(allelename: str, asm_contig: str):
    _, sample_name, hap, _, genome = parse_allele_name(allelename)
    if genome == "HG38_h1" and not accepted_contig(genome, asm_contig):
        return "HG38_h1"
    return genome


def is_hg38_alt_assembly(allelename: str, asm_contig: str) -> bool:
    try:
        _, _, _, _, genome = parse_allele_name(allelename)
    except Exception:
        return False
    return not accepted_contig(genome, asm_contig)


def split_similar_refs(s: str):
    s = s.strip()
    if not s:
        return []
    toks = re.split(r"[;,]", s)
    out = []
    for tok in toks:
        tok = tok.strip()
        if not tok or ":" not in tok:
            continue
        name, val = tok.rsplit(":", 1)
        try:
            out.append((name, float(val)))
        except ValueError:
            continue
    return out


class DSU:
    def __init__(self):
        self.parent = {}

    def add(self, x):
        if x not in self.parent:
            self.parent[x] = x

    def find(self, x):
        self.add(x)
        while self.parent[x] != x:
            self.parent[x] = self.parent[self.parent[x]]
            x = self.parent[x]
        return x

    def union(self, a, b):
        ra = self.find(a)
        rb = self.find(b)
        if ra != rb:
            if ra < rb:
                self.parent[rb] = ra
            else:
                self.parent[ra] = rb


@dataclass(**DATACLASS_OPTIONS)
class RefAllele:
    name: str
    genome: str
    contig: str
    start: int
    end: int
    strand: str
    gene_state: str = ""
    delegate_only: bool = False
    unique_len: int = 0
    sub_start: int = -1
    sub_end: int = -1
    left_offset: int = -1
    right_offset: int = -1
    # Initiated/normal tags. For seeded unique reference alleles these are
    # their direct identity tags; for unseeded reference alleles these may be
    # filled by normal reference propagation.
    left_tag: str = ""
    right_tag: str = ""
    tag_fixed: bool = False

    # Reference-only propagated/control tags. These record what each reference
    # allele end receives from propagation even when the allele also has an
    # initiated identity tag. They are used for comparison only, not seeding.
    propagated_left_tag: str = ""
    propagated_right_tag: str = ""
    # Comparison-only pairs derived independently from the immediately
    # adjacent reference block on each physical side. The upstream and
    # downstream neighbors are normally different blocks, so keep both
    # coherent alternatives instead of combining their endpoint labels into
    # one artificial pair.
    adjacent_tag_pairs: tuple = ()


@dataclass(**DATACLASS_OPTIONS)
class EndRec:
    allele: str
    genome: str
    contig: str
    pos: int
    side: str
    strand: str
    gene_state: str
    tag: str = ""

    def key(self):
        return (self.allele, self.side)


@dataclass(**DATACLASS_OPTIONS)
class AnchorHit:
    refname: str
    contig: str
    pos: int
    strand: str
    side: str


@dataclass(**DATACLASS_OPTIONS)
class AlleleRecord:
    allelename: str
    best_ref_allele: str
    assembly_loc: str
    ref_loc: str
    intron_or_exon: str
    similarity: float
    similarity_raw: str
    class_type: str
    similar_refs_raw: str
    genome: str
    asm_contig: str
    asm_start: int
    asm_end: int
    asm_strand: str
    gene_state: str = ""
    matched_ref_allele: str = ""
    matched_ref_coords: str = ""
    parsed_similar_refs: tuple = ()
    allele_unique_start: int = -1
    allele_unique_end: int = -1
    allele_unique_len: int = 0
    allele_left_offset: object = NEG_INF_STR
    allele_right_offset: object = NEG_INF_STR
    assigned_ref_alleles_bylocation: str = ""
    assigned_ref_alleles_coordinates: str = ""
    assigned_ref_difference: float = math.inf
    grouped_ref_names: tuple = ()
    grouped_ref_original_locs: tuple = ()
    grouped_part_index: int = 0
    grouped_part_count: int = 0
    # A grouped positional assignment is serialized on every part row.  Keep
    # the original query-side intervals and ownership extensions in parallel
    # lists so later stages can compare the union without losing the exact
    # interval belonging to each individual query allele.
    grouped_query_names: tuple = ()
    grouped_query_original_locs: tuple = ()
    grouped_query_left_extensions: tuple = ()
    grouped_query_right_extensions: tuple = ()
    # Assignment pass used by the exclusive pseudo-linear owner election.
    # 1: initial mutual sequence match; 2: coordinate-aware tag pass;
    # 4: distance-normalized tag pass; 3: final positional fallback.
    assignment_stage: int = 0


@dataclass(**DATACLASS_OPTIONS)
class FinalLiftRow:
    fields: tuple
    allele: str
    genome: str
    contig: str
    start: int
    end: int
    strand: str
    assigned_refs: tuple
    class_type: str
    part_index: int
    left_extension: str
    right_extension: str
    assignment_stage: int = 0
    final_stage_present: bool = False
    # Effective column-16 tier and every reference named by the column-17/18
    # location and sequence mappings.
    assignment_tier: str = ""
    mapped_refs: tuple = ()


def final_stage_and_tier(fields):
    """Effective (stage, tier, present) of a GenomeLift row.

    Columns 15/16 join the location and sequence mappings with ";" (e.g.
    ``2;3`` / ``primary;secondary``); the primary one is effective, else the
    first. ``present`` is False for rows without a numeric stage.
    """
    stages = (fields[14] if len(fields) > 14 else "").strip().split(";")
    tiers = (fields[15] if len(fields) > 15 else "").strip().lower().split(";")
    index = tiers.index("primary") if "primary" in tiers else 0
    stage_text = stages[index] if index < len(stages) else stages[0]
    if not stage_text.isdigit():
        return 0, (tiers[0] if tiers else ""), False
    return int(stage_text), (tiers[index] if index < len(tiers) else ""), True


def final_mapping_refs(fields):
    """References named by the column-17/18 mappings (tier|refs|interval|score)."""
    references = []
    for value in fields[16:18] if len(fields) > 17 else ():
        parts = value.strip().split("|")
        if len(parts) >= 2 and parts[1] not in {"", "."}:
            references.extend(
                name.strip() for name in parts[1].split(";") if name.strip()
            )
    return tuple(dict.fromkeys(references))


def final_source_class(fields):
    """Source RefMatch class: column 17 in 17-column files, else column 7."""
    if len(fields) == 17 and fields[16].strip() and "|" not in fields[16]:
        return fields[16].strip()
    return fields[6].strip()


@dataclass
class FinalStageAssignmentRow:
    """One materialized GenomeLift row during the final stage election."""

    fields: list
    allele: str
    genome: str
    matrix: str
    asm_contig: str
    query_tags: tuple
    assigned_refs: tuple
    coordinate_refs: tuple
    kmer_differences: dict
    stage2_refs: tuple = ()
    final_stage: int = 0
    # Parts of one grouped placement share a claim on their reference run.
    group_key: str = ""
    # ``-inf``/``-inf`` rows own no assembly sequence and never claim.
    placeholder: bool = False
    # Reference-haplotype rows describe the reference catalog itself.
    reference_row: bool = False
    # Assembly locus of this PA (or grouped part), used to find the
    # containing ``_align.txt`` row: alignment rows are named by hotspot-local
    # indices that differ from GenomeLift PA names.
    asm_start: int = 0
    asm_end: int = 0
    asm_strand: str = "+"
    # Source RefMatch class (diagnostic only).
    source_class: str = ""
    # part_N of a grouped placement (0 when ungrouped).
    part_index: int = 0
    # Columns 17/18: the location (stage 2) and sequence (stage 3/4)
    # mappings, dict(stage, tier, refs, interval, score); stage 1 sets both.
    location: dict = None
    sequence: dict = None


def _final_stage_part_value(value: str, part_index: int) -> str:
    values = value.split(";")
    if part_index and len(values) >= part_index:
        return values[part_index - 1].strip()
    return value.strip()


def _final_stage_kmer_differences(fields, part_index: int) -> dict:
    """Return the minimum RefMatch difference for every named reference."""
    result = {}
    if len(fields) > 7:
        for token in re.split(r"[;,]", fields[7]):
            token = token.strip()
            if not token or ":" not in token:
                continue
            reference, raw_difference = token.rsplit(":", 1)
            reference = reference.strip().split("[", 1)[0]
            try:
                difference = float(raw_difference)
            except ValueError:
                continue
            if not reference or not math.isfinite(difference) or difference < 0:
                continue
            result[reference] = min(
                difference, result.get(reference, math.inf),
            )

    best = _final_stage_part_value(fields[1], part_index)
    if best.endswith("]") and "[" in best:
        best = best.split("[", 1)[0]
    try:
        best_difference = float(fields[5])
    except (ValueError, IndexError):
        best_difference = math.inf
    if (
        best not in {"", ".", "NA", "na", "None", "none"}
        and math.isfinite(best_difference)
        and best_difference >= 0
    ):
        result[best] = min(
            best_difference, result.get(best, math.inf),
        )
    return result


def _final_stage_anchor_rank(query_tags, reference_options):
    """Use GenomeLift's exact-then-distance-normalized flank comparison."""
    for normalized_pass in (0, 1):
        compared_query = (
            tuple(normalize_tag(tag, ignore_distance=True) for tag in query_tags)
            if normalized_pass else query_tags
        )
        best = None
        for mode, raw_pair in reference_options:
            compared_pair = (
                tuple(normalize_tag(tag, ignore_distance=True) for tag in raw_pair)
                if normalized_pass else raw_pair
            )
            matches = sum(
                bool(query and reference and query == reference)
                for query, reference in zip(compared_query, compared_pair)
            )
            if not matches:
                continue
            score = (
                normalized_pass,
                0 if matches == 2 else 1,
                0 if str(mode) == "initiated" else 1,
            )
            if best is None or score < best:
                best = score
        if best is not None:
            return best
    return None


def _reference_interval_index(ref_alleles):
    """Build sorted per-graph/per-contig indexes for column-10 lookup."""
    result = {}
    for name, reference in ref_alleles.items():
        key = (allele_matrix_name(name), reference.contig)
        result.setdefault(key, []).append((
            reference.start, reference.end, name,
        ))
    for key, values in result.items():
        values.sort(key=lambda value: (value[0], value[1], value[2]))
        prefix_max_end = []
        maximum = -1
        for _start, end, _name in values:
            maximum = max(maximum, end)
            prefix_max_end.append(maximum)
        result[key] = (
            tuple(value[0] for value in values), tuple(values),
            tuple(prefix_max_end),
        )
    return result


def _coordinate_reference_candidates(
    coordinate: str, matrix: str, interval_index,
) -> tuple:
    """Resolve column 10 by overlap without scanning every reference PA."""
    if not coordinate:
        return ()
    try:
        contig, start, end, _strand = parse_loc(coordinate)
    except ValueError:
        return ()
    starts, intervals, prefix_max_end = interval_index.get(
        (matrix, contig), ((), (), ()),
    )
    stop = bisect.bisect_left(starts, end)
    candidates = []
    index = stop - 1
    while index >= 0 and prefix_max_end[index] > start:
        ref_start, ref_end, name = intervals[index]
        overlap = min(end, ref_end) - max(start, ref_start)
        if overlap > 0:
            candidates.append((
                -overlap, ref_end - ref_start, ref_start, ref_end, name,
            ))
        index -= 1
    if not candidates:
        return ()
    return (min(candidates)[-1],)


def _parse_final_stage_rows(path: str, ref_alleles, targets):
    """Read the 14 public columns and retain comments for atomic rewriting."""
    with open(path, "rt") as handle:
        raw_lines = list(handle)
    interval_index = _reference_interval_index(targets)
    parsed = []
    row_by_line = {}
    for line_index, raw in enumerate(raw_lines):
        if not raw.strip() or raw.startswith("#"):
            continue
        fields = raw.rstrip("\r\n").split("\t")
        if fields[0] == "allelename" or fields[0] in {"DEL", "NA"}:
            continue
        if len(fields) < 14:
            raise ValueError(
                f"{path}:{line_index + 1}: expected at least 14 GenomeLift "
                f"columns, found {len(fields)}"
            )
        # Source RefMatch class (grouped rows show part_N in column 7).
        # Re-running over staged output replaces columns 15-17.
        source_class = final_source_class(fields)
        fields = fields[:14]
        part_match = re.fullmatch(r"part_(\d+)", fields[6].strip())
        part_index = int(part_match.group(1)) if part_match else 0
        try:
            asm_contig, asm_start, asm_end, asm_strand = parse_loc(
                _final_stage_part_value(fields[2], part_index)
            )
            genome = parse_allele_name(fields[0])[4]
        except ValueError:
            continue
        matrix = allele_matrix_name(fields[0])
        assigned = tuple(dict.fromkeys(
            reference.strip()
            for reference in fields[8].split(";")
            if reference.strip() in targets
            and allele_matrix_name(reference.strip()) == matrix
        ))
        coordinate = _final_stage_part_value(fields[9], part_index)
        coordinate_refs = _coordinate_overlapping_references(
            coordinate, matrix, interval_index,
        )
        # Every part of a grouped placement repeats the complete query-locus
        # list and merged reference interval, so together they identify the
        # group (the same token graphcigartoref_persample uses).
        group_key = (
            "\x1f".join((
                matrix, asm_contig, fields[9].strip(),
                fields[2].strip() if ";" in fields[2] else "",
            ))
            if part_index else fields[0]
        )
        placeholder = (
            _final_stage_part_value(fields[12], part_index) == "-inf"
            and _final_stage_part_value(fields[13], part_index) == "-inf"
        )
        row = FinalStageAssignmentRow(
            fields=fields,
            allele=fields[0],
            genome=genome,
            matrix=matrix,
            asm_contig=asm_contig,
            query_tags=(fields[10].strip(), fields[11].strip()),
            assigned_refs=assigned,
            coordinate_refs=coordinate_refs,
            kmer_differences=_final_stage_kmer_differences(
                fields, part_index,
            ),
            group_key=group_key,
            placeholder=placeholder,
            reference_row=fields[0] in ref_alleles,
            asm_start=asm_start,
            asm_end=asm_end,
            asm_strand=asm_strand or "+",
            source_class=source_class,
            part_index=part_index,
        )
        parsed.append(row)
        row_by_line[line_index] = row
    return raw_lines, parsed, row_by_line


def _stage_assignment_coordinates(references, ref_alleles):
    resolved = tuple(dict.fromkeys(
        location_refname_for_output(reference, ref_alleles)
        for reference in references
        if location_refname_for_output(reference, ref_alleles)
    ))
    if not resolved:
        return (), ""
    if len(resolved) == 1:
        return resolved, location_coords_for_output(resolved[0], ref_alleles)
    coordinates = [ref_alleles[name] for name in resolved if name in ref_alleles]
    if not coordinates or len({item.contig for item in coordinates}) != 1:
        return resolved, ""
    strands = {item.strand for item in coordinates}
    strand = next(iter(strands)) if len(strands) == 1 else "+"
    return resolved, format_ref_anchor_interval(
        coordinates[0].contig,
        min(item.start for item in coordinates),
        max(item.end for item in coordinates),
        strand,
    )


def _merge_intervals(intervals) -> tuple:
    """Return the sorted union of half-open intervals."""
    if not intervals:
        return ()
    ordered = sorted(intervals)
    merged = []
    current_start, current_end = ordered[0]
    for start, end in ordered[1:]:
        if start <= current_end:
            current_end = max(current_end, end)
        else:
            merged.append((current_start, current_end))
            current_start, current_end = start, end
    merged.append((current_start, current_end))
    return tuple(merged)


def _parse_alignment_position_list(text: str, context: str) -> list:
    intervals = []
    for token in text.split(";"):
        values = token.split("_")
        if len(values) != 2:
            raise ValueError(f"{context}: invalid path interval {token!r}")
        try:
            start, end = map(int, values)
        except ValueError as error:
            raise ValueError(
                f"{context}: invalid path interval {token!r}"
            ) from error
        if start < 0 or end < start:
            raise ValueError(f"{context}: invalid path interval {token!r}")
        intervals.append((start, end))
    return intervals


def read_alignment_projection_rows(path: str, keep_keys=None) -> dict:
    """Index `_align.txt` rows by (graph, contig) for locus-based lookup.

    `_align.txt` query names use hotspot-local indices that differ from
    GenomeLift PA names, so a PA is matched to the alignment row whose
    interval contains its locus, exactly as graphcigartoref_persample does.
    Only columns 2-7, 9 and 10 are sliced; the graph CIGAR in column 8 is
    never materialized.  Each retained row is
    ``(start, end, strand, components)`` where every component is
    ``(path, orientation, path_start, path_end, query_start, query_end)``
    in the row's oriented query coordinates.
    """
    if not path:
        return {}
    keep = set(keep_keys) if keep_keys is not None else None
    result = {}
    alignment_header = "\t".join(_ALIGNMENT_HEADER)
    with open(path, "rt") as handle:
        for line_number, raw in enumerate(handle, 1):
            if not raw.strip() or raw.startswith("#"):
                continue
            line = raw.rstrip("\r\n")
            if line == alignment_header:
                continue
            tabs = [index for index, value in enumerate(line) if value == "\t"]
            if len(tabs) != len(_ALIGNMENT_HEADER) - 1:
                raise ValueError(
                    f"{path}:{line_number}: expected ten alignment columns, "
                    f"found {len(tabs) + 1}"
                )
            query_name = line[tabs[4] + 1:tabs[5]].strip()
            matrix = allele_matrix_name(query_name)
            contig = line[tabs[0] + 1:tabs[1]].strip()
            key = (matrix, contig)
            if keep is not None and key not in keep:
                continue
            context = f"{path}:{line_number}"
            try:
                start = int(line[tabs[1] + 1:tabs[2]])
                end = int(line[tabs[2] + 1:tabs[3]])
            except ValueError as error:
                raise ValueError(f"{context}: invalid alignment interval") from error
            strand = line[tabs[3] + 1:tabs[4]].strip() or "+"
            graph_path = line[tabs[5] + 1:tabs[6]].strip()
            ref_positions = line[tabs[7] + 1:tabs[8]].strip()
            query_positions = line[tabs[8] + 1:].strip()
            components = ()
            if (
                graph_path not in {"", "*"}
                and ref_positions not in {"", "."}
                and query_positions not in {"", "."}
            ):
                paths = _ALIGNMENT_COMPONENT_RE.findall(graph_path)
                path_intervals = _parse_alignment_position_list(
                    ref_positions, context,
                )
                query_intervals = _parse_alignment_position_list(
                    query_positions, context,
                )
                if not (len(paths) == len(path_intervals) == len(query_intervals)):
                    raise ValueError(
                        f"{context}: graph path, reference-position, and "
                        "query-position column counts do not match"
                    )
                components = tuple(
                    (name, orientation, path_start, path_end,
                     query_start, query_end)
                    for (orientation, name), (path_start, path_end),
                    (query_start, query_end)
                    in zip(paths, path_intervals, query_intervals)
                )
            result.setdefault(key, []).append((start, end, strand, components))
    return result


def _select_containing_alignment(rows, start: int, end: int, strand: str):
    """Mirror graphcigartoref_persample.locate_alignment_row's interval rule."""
    candidates = [
        (
            0 if row[2] == strand else 1,
            row[1] - row[0],
            row[0],
            row[1],
            0 if row[3] else 1,
            order,
            row,
        )
        for order, row in enumerate(rows)
        if row[0] <= start and row[1] >= end
    ]
    return min(candidates)[-1] if candidates else None


def project_locus_path_intervals(
    alignment_rows, matrix: str, contig: str, start: int, end: int,
    strand: str,
) -> dict:
    """Project one PA locus onto graph-path intervals of its alignment row.

    The PA's oriented query sub-interval is clipped against every component's
    query range (column 10) and mapped linearly onto that component's path
    range (column 9), respecting ``<`` traversal.  A PA that is only part of a
    longer hotspot alignment is therefore not credited with its neighbours'
    path intervals.
    """
    if end <= start:
        return {}
    row = _select_containing_alignment(
        alignment_rows.get((matrix, contig), ()), start, end, strand,
    )
    if row is None or not row[3]:
        return {}
    row_start, row_end, row_strand, components = row
    if row_strand == "-":
        local_start, local_end = row_end - end, row_end - start
    else:
        local_start, local_end = start - row_start, end - row_start
    intervals_by_path = {}
    for name, orientation, path_start, path_end, query_start, query_end in components:
        clipped_start = max(query_start, local_start)
        clipped_end = min(query_end, local_end)
        if clipped_end <= clipped_start or query_end <= query_start:
            continue
        scale = (path_end - path_start) / (query_end - query_start)
        offset_start = (clipped_start - query_start) * scale
        offset_end = (clipped_end - query_start) * scale
        if orientation == "<":
            projected = (path_end - offset_end, path_end - offset_start)
        else:
            projected = (path_start + offset_start, path_start + offset_end)
        projected_start = max(path_start, int(round(projected[0])))
        projected_end = min(path_end, int(round(projected[1])))
        if projected_end > projected_start:
            intervals_by_path.setdefault(name, []).append(
                (projected_start, projected_end)
            )
    return {
        name: _merge_intervals(intervals)
        for name, intervals in intervals_by_path.items()
    }


def _reference_path_interval_index(reference_intervals) -> dict:
    """Index reference-PA graph-path intervals by path name.

    ``reference_intervals`` maps each reference PA to the graph-path
    intervals projected from the reference haplotype's ``_align.txt``.  Each
    path maps to sorted ``(start, end, reference)`` tuples plus a prefix
    maximum of ends, so a query interval finds every overlapping reference PA
    by binary search.
    """
    by_path = cl.defaultdict(list)
    for reference, intervals_by_path in reference_intervals.items():
        for path_name, intervals in intervals_by_path.items():
            for start, end in intervals:
                by_path[path_name].append((start, end, reference))
    result = {}
    for path_name, values in by_path.items():
        values.sort(key=lambda value: (value[0], value[1], value[2]))
        prefix_max_end = []
        maximum = -1
        for _start, end, _reference in values:
            maximum = max(maximum, end)
            prefix_max_end.append(maximum)
        result[path_name] = (
            tuple(value[0] for value in values),
            tuple(values),
            tuple(prefix_max_end),
        )
    return result


def reference_pa_path_intervals(reference_alignments: str, eligible) -> dict:
    """Project every eligible reference PA onto the reference alignments."""
    keys = {
        (allele_matrix_name(name), reference.contig)
        for name, reference in eligible.items()
    }
    alignment_rows = read_alignment_projection_rows(reference_alignments, keys)
    result = {}
    for name, reference in eligible.items():
        intervals = project_locus_path_intervals(
            alignment_rows, allele_matrix_name(name), reference.contig,
            reference.start, reference.end, reference.strand or "+",
        )
        if intervals:
            result[name] = intervals
    return result


def _stage4_reference_coverages(
    row, path_intervals, reference_path_index,
) -> dict:
    """Sum query/reference-PA overlap on the graph paths both aligned to.

    Query and reference intervals are in the same path-local coordinates, so
    overlap is measured directly on each shared graph path, whichever
    haplotype's path it is.  Only same-graph reference PAs are scored.
    """
    coverages = cl.defaultdict(int)
    for path_name, query_intervals in path_intervals.items():
        indexed = reference_path_index.get(path_name)
        if indexed is None:
            continue
        reference_starts, references, reference_prefix_max_end = indexed
        for query_start, query_end in query_intervals:
            stop = bisect.bisect_left(reference_starts, query_end)
            index = stop - 1
            while (
                index >= 0
                and reference_prefix_max_end[index] > query_start
            ):
                ref_start, ref_end, reference = references[index]
                overlap = min(query_end, ref_end) - max(
                    query_start, ref_start,
                )
                if (
                    overlap > 0
                    and allele_matrix_name(reference) == row.matrix
                ):
                    coverages[reference] += overlap
                index -= 1
    return dict(coverages)


def _target_references(ref_alleles) -> dict:
    """Every reference PA with its own locus, including -inf placeholders."""
    return {
        name: reference
        for name, reference in ref_alleles.items()
        if allele_belongs_to_genome(name, reference.genome)
        and not reference.delegate_only
    }


def _coordinate_overlapping_references(
    coordinate: str, matrix: str, interval_index,
) -> tuple:
    """Return every same-graph reference PA overlapping a column-10 interval."""
    if not coordinate:
        return ()
    try:
        contig, start, end, _strand = parse_loc(coordinate)
    except ValueError:
        return ()
    starts, intervals, prefix_max_end = interval_index.get(
        (matrix, contig), ((), (), ()),
    )
    index = bisect.bisect_left(starts, end) - 1
    names = []
    while index >= 0 and prefix_max_end[index] > start:
        ref_start, ref_end, name = intervals[index]
        if min(end, ref_end) - max(start, ref_start) > 0:
            names.append(name)
        index -= 1
    return tuple(sorted(names))


def _slice_alignment_row(fields, window_start: int, window_end: int) -> dict:
    """Exact path intervals of one query window of an ``_align.txt`` row.

    The graph CIGAR (column 8) is walked component by component from each
    component's query start (column 10) and path range (column 9; ``<``
    components run backwards). Only =/X/M bases inside the window are kept;
    I/S (including flanks) and D contribute no shared base.
    """
    segments = _ALIGNMENT_SEGMENT_RE.findall(fields[7])
    try:
        path_ranges = _parse_alignment_position_list(fields[8], "refpositions")
        query_ranges = _parse_alignment_position_list(fields[9], "qpositions")
    except ValueError:
        return {}
    if not (len(segments) == len(path_ranges) == len(query_ranges)):
        return {}
    intervals = cl.defaultdict(list)
    for (orientation, path, body), (path_start, path_end), (query_start, _query_end) in zip(
        segments, path_ranges, query_ranges,
    ):
        query_cursor = query_start
        path_cursor = path_start if orientation == ">" else path_end
        for length_text, operation, _payload in _ALIGNMENT_OPERATION_RE.findall(body):
            length = int(length_text)
            if operation in {"=", "X", "M"}:
                clip_start = max(query_cursor, window_start)
                clip_end = min(query_cursor + length, window_end)
                if clip_end > clip_start:
                    if orientation == ">":
                        span = (path_cursor + clip_start - query_cursor,
                                path_cursor + clip_end - query_cursor)
                    else:
                        span = (path_cursor - (clip_end - query_cursor),
                                path_cursor - (clip_start - query_cursor))
                    intervals[path].append(span)
                query_cursor += length
                path_cursor += length if orientation == ">" else -length
            elif operation in {"I", "S"}:
                query_cursor += length
            elif operation == "D":
                path_cursor += length if orientation == ">" else -length
    return {path: _merge_intervals(values) for path, values in intervals.items()}


def read_exact_pa_path_intervals(path: str, loci) -> dict:
    """Slice every PA from the ``_align.txt`` row containing its locus.

    ``loci`` maps a PA name to ``(graph, contig, start, end, strand)``.
    Alignment rows are matched by locus, never by name (their query names use
    hotspot-local indices): same graph and contig, the row interval contains
    the PA, preferring the same strand, then the smallest span, as
    graphcigartoref_persample.locate_alignment_row. The file is streamed; only
    the selected PA slices are kept.
    """
    if not path or not loci:
        return {}
    by_key = cl.defaultdict(list)
    for name, (matrix, contig, start, end, strand) in loci.items():
        by_key[(matrix, contig)].append((start, end, strand or "+", name))
    indexed = {}
    for key, values in by_key.items():
        values.sort()
        indexed[key] = ([value[0] for value in values], values)
    best = {}
    with open(path, "rt") as handle:
        for line_number, raw in enumerate(handle, 1):
            if not raw.strip() or raw.startswith("#"):
                continue
            fields = raw.rstrip("\r\n").split("\t")
            if len(fields) != len(_ALIGNMENT_HEADER) or tuple(fields) == _ALIGNMENT_HEADER:
                continue
            entry = indexed.get((allele_matrix_name(fields[5].strip()), fields[1].strip()))
            if entry is None:
                continue
            try:
                row_start, row_end = int(fields[2]), int(fields[3])
            except ValueError:
                continue
            row_strand = fields[4].strip() or "+"
            aligned = (
                fields[6].strip() not in {"", "*"}
                and fields[7].strip() not in {"", "*"}
            )
            starts, values = entry
            index = bisect.bisect_left(starts, row_start)
            while index < len(values) and values[index][0] < row_end:
                start, end, strand, name = values[index]
                index += 1
                if end > row_end:
                    continue
                rank = (
                    0 if row_strand == strand else 1, row_end - row_start,
                    row_start, row_end, 0 if aligned else 1, line_number,
                )
                if name in best and best[name][0] <= rank:
                    continue
                sliced = {}
                if aligned:
                    window = (
                        (row_end - end, row_end - start) if row_strand == "-"
                        else (start - row_start, end - row_start)
                    )
                    sliced = _slice_alignment_row(fields, *window)
                best[name] = (rank, sliced)
    return {name: sliced for name, (_rank, sliced) in best.items()}


def _reference_overlap_intervals(query_intervals, reference_index, matrix) -> dict:
    """Map each same-graph reference PA to its overlap with the query slice."""
    overlaps = cl.defaultdict(list)
    for path_name, intervals in query_intervals.items():
        indexed = reference_index.get(path_name)
        if indexed is None:
            continue
        starts, references, prefix_max_end = indexed
        for query_start, query_end in intervals:
            index = bisect.bisect_left(starts, query_end) - 1
            while index >= 0 and prefix_max_end[index] > query_start:
                ref_start, ref_end, reference = references[index]
                start = max(query_start, ref_start)
                end = min(query_end, ref_end)
                if end > start and allele_matrix_name(reference) == matrix:
                    overlaps[reference].append((path_name, start, end))
                index -= 1
    return overlaps


def _union_overlap(overlaps, references) -> int:
    """Query bases overlapping any reference PA of a (multi-matching) set."""
    by_path = cl.defaultdict(list)
    for reference in references:
        for path_name, start, end in overlaps.get(reference, ()):
            by_path[path_name].append((start, end))
    return sum(
        end - start
        for intervals in by_path.values()
        for start, end in _merge_intervals(intervals)
    )


def _fill_reference_run(members, targets, interval_index, available):
    """Extend a reference set to its whole run, filling gaps.

    Same-graph reference PAs lying mostly (>= half of their length) inside
    the members' span on the same contig join the run, provided they have
    graph-path intervals in the reference ``_align.txt`` (``available``).
    Returns None when the members are not on one graph and contig.
    """
    if any(name not in targets for name in members):
        return None
    references = [targets[name] for name in members]
    matrices = {allele_matrix_name(name) for name in members}
    contigs = {reference.contig for reference in references}
    if len(matrices) != 1 or len(contigs) != 1:
        return None
    matrix, contig = next(iter(matrices)), next(iter(contigs))
    start = min(reference.start for reference in references)
    end = max(reference.end for reference in references)
    filled = set(members)
    for name in _coordinate_overlapping_references(
        f"{contig}:{start}-{end}", matrix, interval_index,
    ):
        if name in filled or name not in targets or name not in available:
            continue
        other = targets[name]
        inside = min(end, other.end) - max(start, other.start)
        if inside * 2 >= max(1, other.end - other.start):
            filled.add(name)
    return tuple(sorted(filled))


def _elect_final_stages(
    rows, ref_alleles, targets, overlaps_by_row, merged_spans=None,
    available=None,
):
    """Four-stage lift election per genome, on query units.

    A unit is either one query PA, or a merged query: all parts of one
    GenomeLift grouped placement (part_N rows) or one query placed on several
    references (a multi-reference column 9). A merged query is scored as one
    unit against its merged reference. Gaps are part of both intervals: the
    merged reference is the whole run of its references
    (_fill_reference_run), and the merged query is the span from its first to
    its last part, sliced from the alignment row containing it
    (``merged_spans``; the parts' own slices when no row contains the span).
    Coverage is the overlap of the merged query with the merged reference
    divided by the merged query length; the unit is anchored by its column-9
    placement and has no k-mer similarity. Every part receives the unit's
    stage, tier and references.

    A single query PA is scored against single reference PAs and one
    coverage-based multi-matching set: the references each covering >=
    FINAL_STAGE_MULTI_MIN_SHARE of the query and >=
    FINAL_STAGE_MULTI_MIN_REFERENCE_COVERAGE of their own length (paralogs
    covering the same query bases excluded), filled to their whole run
    (_fill_reference_run). Similarity is 1 - the RefMatch k-mer difference
    (best member of a set). A candidate is anchored when a member is named in
    column 9, overlaps the column-10 interval, or matches the query's
    propagated anchor tags.

    Lifting has two sources that may disagree: location (anchors) and
    sequence (k-mer similarity, coverage). Each row gets a location mapping
    and/or a sequence mapping:

      1  location and sequence agree: coverage >= 0.9 and similarity >= 0.9
         after a +0.2 anchor bonus to both, the unit's highest-similarity
         candidate. Fills both mappings.
      2  location: the unit's anchors, with no coverage or k-mer condition.
         GenomeLift's positional column 9 (a merged query's gap-filled run),
         else its column-10 interval (the references inside it, possibly
         none), else the best propagated-tag match. Units competing for a
         location are ordered by coverage on it; the first is primary.
      3  sequence. Beside a stage-2 location (each query PA, merged-query
         members included): its best k-mer reference at similarity >= 0.9
         when that disagrees with the location. It is a duplicate unless
         the location is a duplicate, the best distance is 1.5x clear of the
         second and the reference is unclaimed; then it is primary and the
         location is dropped (translocation). Without a location: coverage
         >= 0.5 or similarity > 0.5, promoted to primary when the reference
         is unclaimed and the best overlap exceeds 1.5x the second (or 1.5x
         the best distance is below the second); otherwise a duplicate.
      4  sequence, without stages 1-3: >= 2,000 overlapping bases, duplicate.

    Each reference PA is claimed by one primary mapping. A set claims only
    the references it matched, not the gap references filled between them
    (those count for coverage and column 9 only). Units owning no sequence
    (all parts -inf/-inf) may hold a stage-2 location but never a primary
    sequence mapping.
    """
    reference_tag_pairs = reference_tag_pair_options(
        ref_alleles, ignore_distance=False,
    )
    interval_index = _reference_interval_index(targets)
    merged_spans = merged_spans or {}
    available = set(targets) if available is None else available
    rows_by_genome = cl.defaultdict(list)
    for row in rows:
        rows_by_genome[row.genome].append(row)

    def single_unit(row):
        length = max(1, row.asm_end - row.asm_start)
        overlaps = overlaps_by_row.get(row.allele, {})
        overlap_bp = {
            reference: _union_overlap(overlaps, (reference,))
            for reference in overlaps
        }
        similarity = {
            reference: 1.0 - difference
            for reference, difference in row.kmer_differences.items()
            if reference in targets
            and allele_matrix_name(reference) == row.matrix
        }
        singles = set(overlap_bp) | set(similarity)
        anchored = set(row.assigned_refs) | set(row.coordinate_refs)
        for reference in singles:
            if _final_stage_anchor_rank(
                row.query_tags, reference_tag_pairs.get(reference, ()),
            ) is not None:
                anchored.add(reference)
        candidate_sets = {(reference,) for reference in singles}
        claims = {}
        eligible = sorted((
            reference for reference, value in overlap_bp.items()
            if value >= FINAL_STAGE_MULTI_MIN_SHARE * length
            and value >= FINAL_STAGE_MULTI_MIN_REFERENCE_COVERAGE * max(
                1, targets[reference].end - targets[reference].start,
            )
        ), key=lambda reference: (-overlap_bp[reference], reference))
        # A single-to-many neighbour covers other query bases; a paralog
        # (sharing the graph path) covers the same ones. From the largest
        # overlap down, a reference joins only when at most
        # FINAL_STAGE_MULTI_MAX_SHARED of its overlap is already covered by
        # the members chosen so far.
        chosen = []
        for reference in eligible:
            added = (
                _union_overlap(overlaps, chosen + [reference])
                - _union_overlap(overlaps, chosen)
            )
            if (
                not chosen
                or overlap_bp[reference] - added
                <= FINAL_STAGE_MULTI_MAX_SHARED * overlap_bp[reference]
            ):
                chosen.append(reference)
        shared = tuple(sorted(chosen))
        if len(shared) > 1:
            filled = _fill_reference_run(
                shared, targets, interval_index, available,
            )
            if filled:
                candidate_sets.add(filled)
                # Gap references count for coverage only; the unit claims
                # just the references it matched.
                claims[filled] = shared
        candidates = []
        for members in candidate_sets:
            coverage = _union_overlap(overlaps, members) / length
            best_similarity = max(
                (similarity.get(reference, 0.0) for reference in members),
                default=0.0,
            )
            bonus = (
                FINAL_STAGE_ANCHOR_BONUS
                if anchored.intersection(members) else 0.0
            )
            candidates.append(dict(
                members=members, coverage=coverage,
                similarity=best_similarity, anchored=bool(bonus),
                effective_coverage=coverage + bonus,
                effective_similarity=best_similarity + bonus,
            ))
        units = [
            ((reference,), value) for reference, value in overlap_bp.items()
        ] + [
            (members, _union_overlap(overlaps, members)) for members in claims
        ]
        units.sort(key=lambda item: (-item[1], item[0]))
        # Location: GenomeLift's positional placement (column 9, else its
        # column-10 interval), else the best propagated-tag match.
        coordinate = _final_stage_part_value(row.fields[9], row.part_index)
        location = None
        if row.assigned_refs:
            location = dict(refs=tuple(row.assigned_refs),
                            interval=coordinate, source="col9")
        elif coordinate:
            location = dict(refs=tuple(row.coordinate_refs),
                            interval=coordinate, source="col10")
        else:
            ranks = {}
            for reference in singles:
                rank = _final_stage_anchor_rank(
                    row.query_tags, reference_tag_pairs.get(reference, ()),
                )
                if rank is not None:
                    ranks[reference] = rank
            if ranks:
                best_rank = min(ranks.values())
                location = dict(refs=tuple(sorted(
                    reference for reference, rank in ranks.items()
                    if rank == best_rank
                )), interval="", source="tags")
        return dict(
            rows=[row], key=row.group_key, name=row.allele,
            placeholder=row.placeholder, length=length,
            similarity=similarity, candidates=candidates, units=units,
            claims=claims, merged=False, location=location,
            overlap=lambda references: _union_overlap(overlaps, references),
        )

    def merged_unit(unit_rows):
        named = tuple(sorted({
            reference for row in unit_rows for reference in row.assigned_refs
        }))
        references = _fill_reference_run(
            named, targets, interval_index, available,
        ) or named
        span = merged_spans.get(unit_rows[0].group_key)
        outside = cl.Counter()
        if span is not None:
            span_key, length = span
            overlaps = overlaps_by_row.get(span_key, {})

            def overlap_of(members):
                return _union_overlap(overlaps, members)

            for reference in overlaps:
                if reference not in references:
                    outside[reference] = _union_overlap(overlaps, (reference,))
        else:
            length = sum(max(1, row.asm_end - row.asm_start) for row in unit_rows)

            def overlap_of(members):
                return sum(
                    _union_overlap(overlaps_by_row.get(row.allele, {}), members)
                    for row in unit_rows
                )

            for row in unit_rows:
                overlaps = overlaps_by_row.get(row.allele, {})
                for reference in overlaps:
                    if reference not in references:
                        outside[reference] += _union_overlap(
                            overlaps, (reference,),
                        )
        overlap = overlap_of(references)
        length = max(1, length)
        # Location: the grouped placement's column 9 (gap-filled) and its
        # merged column-10 interval; without a same-graph column 9, the
        # references inside that interval.
        coordinate = _final_stage_part_value(
            unit_rows[0].fields[9], unit_rows[0].part_index,
        )
        location_refs = references or tuple(sorted({
            reference for row in unit_rows for reference in row.coordinate_refs
        }))
        location = (
            dict(refs=location_refs, interval=coordinate, source="col9")
            if location_refs or coordinate else None
        )
        coverage = overlap / length
        units = [(references, overlap)] + [
            ((reference,), value) for reference, value in outside.items()
        ]
        units.sort(key=lambda item: (-item[1], item[0]))
        return dict(
            rows=list(unit_rows), key=unit_rows[0].group_key,
            name=min(row.allele for row in unit_rows),
            placeholder=all(row.placeholder for row in unit_rows),
            length=length, similarity={},
            candidates=[dict(
                members=references, coverage=coverage, similarity=0.0,
                anchored=True,
                effective_coverage=coverage + FINAL_STAGE_ANCHOR_BONUS,
                effective_similarity=FINAL_STAGE_ANCHOR_BONUS,
            )] if references else [],
            # Gap references count for coverage only; the unit claims just
            # its column-9 references.
            units=units, claims={references: named}, merged=True,
            location=location, overlap=overlap_of,
        )

    for genome_rows in rows_by_genome.values():
        units = []
        grouped = cl.defaultdict(list)
        for row in genome_rows:
            if row.reference_row:
                if row.allele in row.assigned_refs:
                    mapping = dict(
                        stage=1, tier="primary", refs=(row.allele,),
                        interval=row.fields[9].strip(), score="reference",
                    )
                    row.location, row.sequence = mapping, dict(mapping)
            elif row.part_index or len(row.assigned_refs) > 1:
                grouped[row.group_key].append(row)
            else:
                units.append(single_unit(row))
        units.extend(merged_unit(unit_rows) for unit_rows in grouped.values())

        # Reference PA -> claiming key (a unit's key, or a PA name for a
        # member's own stage-3 sequence mapping).
        claimed = {}

        def claim_set(unit, members):
            return unit["claims"].get(tuple(members), tuple(members))

        def free(key, references):
            return all(
                claimed.get(reference, key) == key for reference in references
            )

        def claim(key, references):
            for reference in references:
                claimed[reference] = key

        def prefix(unit, rule):
            return ("merged_" if unit["merged"] else "") + rule

        def mapped(unit):
            row = unit["rows"][0]
            return row.location is not None or row.sequence is not None

        def set_sequence(unit, members, stage, tier, score):
            for row in unit["rows"]:
                row.sequence = dict(
                    stage=stage, tier=tier, refs=tuple(members), interval="",
                    score=prefix(unit, score),
                )
            if tier == "primary":
                claim(unit["key"], claim_set(unit, members))

        def candidate_score(candidate):
            return (
                f"cov={candidate['coverage']:.3f},"
                f"sim={candidate['similarity']:.3f},"
                f"anchored={int(candidate['anchored'])}"
            )

        # Stage 1: sequence and location agree (anchor bonus); the unit's
        # highest-similarity candidate. Fills both mappings.
        stage1 = []
        for unit in units:
            if unit["placeholder"]:
                continue
            passing = [
                candidate for candidate in unit["candidates"]
                if candidate["effective_coverage"] >= FINAL_STAGE_HIGH_COVERAGE
                and candidate["effective_similarity"]
                >= FINAL_STAGE_HIGH_KMER_SIMILARITY
            ]
            if passing:
                best = max(passing, key=lambda candidate: (
                    candidate["effective_similarity"],
                    candidate["effective_coverage"],
                    -len(candidate["members"]),
                ))
                stage1.append((
                    -best["effective_similarity"], -best["effective_coverage"],
                    unit["name"], unit, best,
                ))
        for *_key, unit, best in sorted(stage1, key=lambda item: item[:3]):
            if free(unit["key"], claim_set(unit, best["members"])):
                score = prefix(unit, "high_confidence,") + candidate_score(best)
                for row in unit["rows"]:
                    row.location = dict(
                        stage=1, tier="primary", refs=tuple(best["members"]),
                        interval="", score=score,
                    )
                    row.sequence = dict(row.location)
                claim(unit["key"], claim_set(unit, best["members"]))

        # Stage 2: location by anchors, no coverage or k-mer condition.
        # Units competing for a location are ordered by coverage on it.
        stage2 = []
        for unit in units:
            if mapped(unit) or unit["location"] is None:
                continue
            location = unit["location"]
            coverage = unit["overlap"](location["refs"]) / unit["length"]
            stage2.append((-coverage, unit["name"], unit, location, coverage))
        for _coverage, _name, unit, location, coverage in sorted(
            stage2, key=lambda item: item[:2],
        ):
            references = claim_set(unit, location["refs"])
            tier = "primary" if free(unit["key"], references) else "secondary"
            if tier == "primary":
                claim(unit["key"], references)
            for row in unit["rows"]:
                row.location = dict(
                    stage=2, tier=tier, refs=tuple(location["refs"]),
                    interval=location["interval"],
                    score=prefix(unit, f"anchor,source={location['source']},"
                                       f"cov={coverage:.3f}"),
                )

        # Stage 3 beside a location: each query PA's own k-mer best
        # reference (>= 0.9) where it disagrees with the location. Primary
        # only when the location is a duplicate, the best distance is 1.5x
        # clear of the second and the reference is unclaimed; the location
        # is then dropped (translocation).
        hard = []
        for unit in units:
            for row in unit["rows"]:
                if row.location is None or row.location["stage"] != 2:
                    continue
                by_distance = sorted(
                    (difference, reference)
                    for reference, difference in row.kmer_differences.items()
                    if reference in targets
                    and allele_matrix_name(reference) == row.matrix
                )
                if (
                    not by_distance
                    or 1.0 - by_distance[0][0] < FINAL_STAGE_HIGH_KMER_SIMILARITY
                    or by_distance[0][1] in row.location["refs"]
                ):
                    continue
                second = by_distance[1][0] if len(by_distance) > 1 else math.inf
                hard.append((by_distance[0][0], row.allele, unit, row,
                             by_distance[0][1], second))
        for distance, _name, unit, row, reference, second in sorted(
            hard, key=lambda item: item[:2],
        ):
            clear = FINAL_STAGE_PROMOTION_MARGIN * distance < second
            tier = (
                "primary"
                if row.location["tier"] == "secondary" and clear
                and not row.placeholder and free(row.allele, (reference,))
                else "secondary"
            )
            row.sequence = dict(
                stage=3, tier=tier, refs=(reference,), interval="",
                score=f"kmer_high,sim={1.0 - distance:.3f},"
                      f"dist={distance:.4f},second_dist={second:.4f}",
            )
            if tier == "primary":
                claim(row.allele, (reference,))
                row.location = None

        # Stage 3 without a location: coverage >= 0.5 or similarity > 0.5;
        # promote a clear unclaimed best reference, otherwise a duplicate.
        stage3 = []
        for unit in units:
            if mapped(unit):
                continue
            unit_overlaps = unit["units"]
            by_distance = sorted(
                ((1.0 - value, reference)
                 for reference, value in unit["similarity"].items()),
            )
            coverage_unit = (
                unit_overlaps[0][0]
                if unit_overlaps and unit_overlaps[0][1] / unit["length"]
                >= FINAL_STAGE_MIN_COVERAGE else None
            )
            kmer_unit = (
                (by_distance[0][1],)
                if by_distance
                and 1.0 - by_distance[0][0] > FINAL_STAGE_MIN_KMER_SIMILARITY
                else None
            )
            if coverage_unit is None and kmer_unit is None:
                continue
            promotable = promoted_route = None
            routes = {}
            if coverage_unit is not None:
                # A set's own members do not compete with it.
                second = max((
                    value for members, value in unit_overlaps[1:]
                    if not set(members) <= set(coverage_unit)
                ), default=0)
                routes["coverage"] = (
                    f"cov={unit_overlaps[0][1] / unit['length']:.3f},"
                    f"second_cov={second / unit['length']:.3f}"
                )
                if unit_overlaps[0][1] > FINAL_STAGE_PROMOTION_MARGIN * second:
                    promotable, promoted_route = coverage_unit, "coverage"
            if kmer_unit is not None:
                second = by_distance[1][0] if len(by_distance) > 1 else math.inf
                routes["kmer"] = (
                    f"dist={by_distance[0][0]:.4f},second_dist={second:.4f}"
                )
                if (
                    promotable is None
                    and FINAL_STAGE_PROMOTION_MARGIN * by_distance[0][0] < second
                ):
                    promotable, promoted_route = kmer_unit, "kmer"
            selected = coverage_unit if coverage_unit is not None else kmer_unit
            chosen = promotable or selected
            stage3.append((
                -dict(unit_overlaps).get(chosen, 0),
                1.0 - max(
                    (unit["similarity"].get(reference, 0.0)
                     for reference in chosen), default=0.0,
                ),
                unit["name"], unit, promotable, selected, routes,
                promoted_route, "coverage" if coverage_unit else "kmer",
            ))
        for *_key, unit, promotable, selected, routes, promoted_route, \
                selected_route in sorted(
            stage3, key=lambda item: item[:3],
        ):
            if (
                promotable is not None and not unit["placeholder"]
                and free(unit["key"], claim_set(unit, promotable))
            ):
                set_sequence(unit, promotable, 3, "primary",
                             f"{promoted_route},{routes[promoted_route]}")
            else:
                reason = (
                    "placeholder" if unit["placeholder"] and promotable
                    else "claimed" if promotable else "unclear"
                )
                set_sequence(
                    unit, selected, 3, "secondary",
                    f"{selected_route},{routes[selected_route]},reason={reason}",
                )

        # Stage 4: >= 2,000 overlapping bases, always a duplicate (sequence).
        for unit in units:
            if mapped(unit):
                continue
            unit_overlaps = unit["units"]
            if (
                unit_overlaps
                and unit_overlaps[0][1] >= FINAL_STAGE_MIN_PATH_COVERAGE
            ):
                set_sequence(
                    unit, unit_overlaps[0][0], 4, "secondary",
                    f"overlap,overlap_bp={unit_overlaps[0][1]},"
                    f"cov={unit_overlaps[0][1] / unit['length']:.3f}",
                )


def finalize_lift_stages(
    path: str,
    ref_alleles,
    alignments: str = "",
    reference_alignments: str = "",
):
    """Elect final mappings, rewrite columns 9/10 and 15-18, return counts.

    Coverage needs both ``alignments`` (the sample ``_align.txt``) and
    ``reference_alignments`` (the reference haplotype ``_align.txt``); PAs are
    matched to alignment rows by locus, never by name.

      9/10   the primary mapping's references and interval (the location
             when neither mapping is primary)
      15/16  stage and tier of the location and sequence mappings, joined
             with ";" in that order (``2;3`` / ``primary;secondary``);
             stage 1 is ``1`` / ``primary``; ``0`` / "" when unmapped
      17     location mapping (stage 1/2): ``tier|references|interval|score``
      18     sequence mapping (stage 1/3/4), same format
             (``.`` when absent; references ;-joined; score starts with the
             rule: high_confidence, anchor, kmer_high, coverage, kmer,
             overlap, reference; merged_ for merged queries)

    Counts are by the stage shown in columns 9/10.
    """
    targets = _target_references(ref_alleles)
    raw_lines, rows, row_by_line = _parse_final_stage_rows(
        path, ref_alleles, targets,
    )
    query_loci = {
        row.allele: (row.matrix, row.asm_contig, row.asm_start, row.asm_end,
                     row.asm_strand)
        for row in rows if not row.reference_row
    }
    # Merged queries (grouped parts, multi-reference column 9) are scored on
    # their whole span, gaps between parts included.
    merged_rows = cl.defaultdict(list)
    for row in rows:
        if not row.reference_row and (
            row.part_index or len(row.assigned_refs) > 1
        ):
            merged_rows[row.group_key].append(row)
    span_loci = {}
    for key, group_rows in merged_rows.items():
        if len({row.asm_contig for row in group_rows}) != 1:
            continue
        first = min(group_rows, key=lambda row: (row.asm_start, row.asm_end))
        span_loci["\x00merged\x1f" + key] = (key, (
            first.matrix, first.asm_contig,
            min(row.asm_start for row in group_rows),
            max(row.asm_end for row in group_rows),
            first.asm_strand,
        ))
    query_loci.update({
        span_key: locus for span_key, (_key, locus) in span_loci.items()
    })
    reference_loci = {
        name: (allele_matrix_name(name), reference.contig, reference.start,
               reference.end, reference.strand or "+")
        for name, reference in targets.items()
    }
    query_slices = read_exact_pa_path_intervals(alignments, query_loci)
    reference_slices = {
        name: intervals
        for name, intervals in read_exact_pa_path_intervals(
            reference_alignments, reference_loci,
        ).items()
        if intervals
    } if query_slices else {}
    reference_index = _reference_path_interval_index(reference_slices)
    overlaps_by_row = {
        row.allele: _reference_overlap_intervals(
            query_slices.get(row.allele, {}), reference_index, row.matrix,
        )
        for row in rows if not row.reference_row
    }
    merged_spans = {}
    for span_key, (key, locus) in span_loci.items():
        if not query_slices.get(span_key):
            # No single alignment row contains the whole span: the parts'
            # own slices are used instead.
            continue
        overlaps_by_row[span_key] = _reference_overlap_intervals(
            query_slices[span_key], reference_index, locus[0],
        )
        merged_spans[key] = (span_key, locus[3] - locus[2])
    _elect_final_stages(
        rows, ref_alleles, targets, overlaps_by_row,
        merged_spans=merged_spans, available=set(reference_slices),
    )

    counts = cl.defaultdict(int)
    for row in rows:
        fields = row.fields
        mappings = []
        for mapping in (row.location, row.sequence):
            if mapping is None:
                continue
            references, coordinates = _stage_assignment_coordinates(
                mapping["refs"], ref_alleles,
            )
            mapping["refs"] = references
            # A location keeps GenomeLift's lifted interval when it has one.
            mapping["interval"] = mapping["interval"] or coordinates
            mappings.append(mapping)
        if mappings:
            # Columns 9/10: the primary mapping (location first).
            shown = next(
                (mapping for mapping in mappings
                 if mapping["tier"] == "primary"),
                mappings[0],
            )
            fields[8] = ";".join(shown["refs"])
            fields[9] = shown["interval"]
            row.final_stage = shown["stage"]
        elif not row.reference_row:
            fields[8] = ""
        if row.location is not None and row.location["stage"] == 1:
            stage_text, tier_text = "1", "primary"
        else:
            stage_text = ";".join(str(m["stage"]) for m in mappings) or "0"
            tier_text = ";".join(m["tier"] for m in mappings)

        def mapping_text(mapping):
            if mapping is None:
                return "."
            return "|".join((
                mapping["tier"], ";".join(mapping["refs"]) or ".",
                mapping["interval"] or ".", mapping["score"],
            ))

        fields.extend((
            stage_text, tier_text,
            mapping_text(row.location), mapping_text(row.sequence),
        ))
        counts[row.final_stage] += 1

    temporary = path + f".final-stages.tmp.{os.getpid()}"
    try:
        with open(temporary, "wt") as output:
            for line_index, raw in enumerate(raw_lines):
                row = row_by_line.get(line_index)
                if row is None:
                    if raw.rstrip("\r\n").split("\t", 1)[0] == "allelename":
                        header = raw.rstrip("\r\n").split("\t")[:14]
                        output.write("\t".join(
                            header + list(FINAL_STAGE_COLUMNS)
                        ) + "\n")
                    else:
                        output.write(raw)
                else:
                    output.write("\t".join(row.fields) + "\n")
        os.replace(temporary, path)
    finally:
        try:
            os.remove(temporary)
        except FileNotFoundError:
            pass
    return dict(counts)


def _final_extension_value(token: str) -> str:
    token = (token or "").strip()
    if ":" in token:
        token = token.rsplit(":", 1)[1].strip()
    return token


def final_lift_effective_interval(
    row: FinalLiftRow,
    max_extension: int,
):
    """Apply output columns 13/14, returning None for unowned rows."""
    values = []
    for token in (row.left_extension, row.right_extension):
        token = _final_extension_value(token)
        if token in {"-inf", "+inf", "inf"}:
            return None
        if token in {"", ".", "NA", "None", "none", "null"}:
            values.append(("", 0, 0))
            continue
        sign = token[0] if token[:1] in {"+", "-"} else ""
        body = token[1:] if sign else token
        if not body.isdigit():
            raise ValueError(
                f"{row.allele}: malformed GenomeLift extension {token!r}"
            )
        raw = int(body)
        if sign == "-":
            effective = raw
        elif sign == "+":
            effective = min(raw, max_extension)
        else:
            effective = min(
                max(0, raw - min(raw, max_extension)), max_extension,
            )
        values.append((sign, raw, effective))

    (left_sign, left_raw, left_effective), (
        right_sign, right_raw, right_effective,
    ) = values
    start = (
        row.start + left_raw
        if left_sign == "-"
        else max(0, row.start - left_effective)
    )
    end = (
        row.end - right_raw
        if right_sign == "-"
        else row.end + right_effective
    )
    if end <= start:
        return None
    return start, end


def read_final_lift_rows(path: str) -> List[FinalLiftRow]:
    rows = []
    with open(path, "rt") as handle:
        for raw in handle:
            if not raw.strip() or raw.startswith("#"):
                continue
            fields = raw.rstrip("\r\n").split("\t")
            if fields[0] == "allelename" or len(fields) < 14:
                continue
            match = re.fullmatch(r"part_(\d+)", fields[6].strip())
            part_index = int(match.group(1)) if match else 0

            def part_value(value: str) -> str:
                values = value.split(";")
                if part_index and len(values) >= part_index:
                    return values[part_index - 1].strip()
                return value.strip()

            try:
                contig, start, end, strand = parse_loc(part_value(fields[2]))
                genome = parse_allele_name(fields[0])[4]
            except ValueError:
                continue
            assigned_refs = tuple(
                value.strip() for value in fields[8].split(";")
                if value.strip()
            )
            stage, tier, stage_present = final_stage_and_tier(fields)
            rows.append(FinalLiftRow(
                fields=tuple(fields),
                allele=fields[0],
                genome=genome,
                contig=contig,
                start=start,
                end=end,
                strand=strand,
                assigned_refs=assigned_refs,
                class_type=final_source_class(fields),
                part_index=part_index,
                left_extension=part_value(fields[12]),
                right_extension=part_value(fields[13]),
                assignment_stage=stage,
                final_stage_present=stage_present,
                assignment_tier=tier,
                mapped_refs=final_mapping_refs(fields),
            ))
    return rows


def _positive_interval_overlap(left, right):
    return min(left[1], right[1]) > max(left[0], right[0])


def _choose_final_deletion_anchors(
    left_rows, right_rows, max_extension=10000,
):
    """Choose flanking query anchors whose extended regions overlap.

    The two sorted sweeps avoid the former all-pairs comparison.  A missing
    reference block is a deletion only when both neighboring reference blocks
    and both query anchors have overlapping effective extension intervals.
    """
    left_by_contig = cl.defaultdict(list)
    right_by_contig = cl.defaultdict(list)
    for row in left_rows:
        effective = final_lift_effective_interval(row, max_extension)
        if effective is not None:
            left_by_contig[row.contig].append((effective, row))
    for row in right_rows:
        effective = final_lift_effective_interval(row, max_extension)
        if effective is not None:
            right_by_contig[row.contig].append((effective, row))

    candidates = []
    for contig in sorted(set(left_by_contig).intersection(right_by_contig)):
        left_local = sorted(
            left_by_contig[contig],
            key=lambda item: (item[0][0], item[0][1], item[1].allele),
        )
        right_local = sorted(
            right_by_contig[contig],
            key=lambda item: (item[0][0], item[0][1], item[1].allele),
        )
        left_index = right_index = 0
        while left_index < len(left_local) and right_index < len(right_local):
            left_effective, left = left_local[left_index]
            right_effective, right = right_local[right_index]
            if _positive_interval_overlap(left_effective, right_effective):
                if left.start + left.end <= right.start + right.end:
                    strand = "+"
                    left_facing, right_facing = left.end, right.start
                else:
                    strand = "-"
                    left_facing, right_facing = left.start, right.end
                candidates.append((
                    abs(right_facing - left_facing),
                    left.contig,
                    min(left.start, right.start),
                    max(left.end, right.end),
                    left.allele,
                    right.allele,
                    left,
                    right,
                    max(0, (left_facing + right_facing) // 2),
                    strand,
                ))
            if left_effective[1] <= right_effective[1]:
                left_index += 1
            else:
                right_index += 1
    if not candidates:
        return None
    candidates.sort(key=lambda value: value[:6])
    value = candidates[0]
    return value[6], value[7], value[8], value[9]


def final_row_reference_presence(
    row: FinalLiftRow,
    refhaplo: str,
):
    """Return reference blocks whose presence is supported by this query row.

    Assigned column-9 references are the strongest evidence. When column 9 is
    blank, retain the RefMatch best hit from column 2 solely as
    deletion-suppression evidence. This includes both unowned ``-inf`` rows and
    finite rows whose propagated tags infer only a boundary (possibly a
    zero-width column-10 interval). A column-2 fallback is never promoted into
    an anchor or graph-CIGAR assignment.

    Final-stage tables follow the same rule: column 9 holds the elected
    reference of every stage-1 to stage-4 row, and a stage-0 row (blank
    column 9) falls back to column 2.  A stage-3/4 copy therefore suppresses a
    DEL at its source; only primary rows (stages 1/2 and promoted 3) may act
    as deletion anchors.
    """
    if row.assigned_refs or row.mapped_refs:
        # Columns 17/18 may name a second (location or sequence) target.
        return tuple(dict.fromkeys(row.assigned_refs + row.mapped_refs))
    if len(row.fields) < 2:
        return ()

    references = []
    for token in row.fields[1].split(";"):
        token = token.strip()
        # GenomeLift can display a gene-state annotation as REF[STATE].
        token = re.sub(r"\[[^\[\]]*\]$", "", token)
        if token and allele_belongs_to_genome(token, refhaplo):
            references.append(token)
    return tuple(references)


def build_final_unmatched_rows(
    rows: List[FinalLiftRow],
    refhaplo: str,
    max_extension: int = 10000,
):
    """Create DEL/NA status rows for unmatched reference runs."""
    references_by_chrom = cl.defaultdict(list)
    for row in rows:
        if row.genome != refhaplo or row.allele not in row.assigned_refs:
            continue
        effective = final_lift_effective_interval(row, max_extension)
        if effective is None:
            continue
        references_by_chrom[row.contig].append((row, effective))

    query_genomes = sorted({
        row.genome for row in rows
        if row.genome != refhaplo and row.class_type.lower() != "del"
    })
    output = []
    for query_genome in query_genomes:
        matched_refs = set()
        anchors = cl.defaultdict(list)
        for row in rows:
            if row.genome != query_genome:
                continue
            effective = final_lift_effective_interval(row, max_extension)
            presence_refs = final_row_reference_presence(row, refhaplo)
            matched_refs.update(presence_refs)
            if effective is None:
                continue
            if (
                row.part_index == 0
                and len(row.assigned_refs) == 1
                and (
                    not row.final_stage_present
                    # Only primary mappings anchor: stage 1, a primary
                    # stage-2 location, or a primary stage-3 sequence.
                    or (
                        row.assignment_tier == "primary"
                        and row.assignment_stage in {1, 2, 3}
                    )
                )
            ):
                anchors[row.assigned_refs[0]].append(row)

        for chrom, entries in sorted(references_by_chrom.items()):
            entries = sorted(entries, key=lambda value: (
                value[0].start, value[0].end, value[0].allele,
            ))
            index = 0
            while index < len(entries):
                if entries[index][0].allele in matched_refs:
                    index += 1
                    continue
                run_start = index
                while (
                    index < len(entries)
                    and entries[index][0].allele not in matched_refs
                ):
                    index += 1
                run_end = index
                missing = entries[run_start:run_end]
                ref_start = min(effective[0] for _row, effective in missing)
                ref_end = max(effective[1] for _row, effective in missing)
                if ref_end <= ref_start:
                    continue

                chosen = None
                if run_start > 0 and run_end < len(entries):
                    left_ref = entries[run_start - 1][0]
                    right_ref = entries[run_end][0]
                    left_reference_effective = entries[run_start - 1][1]
                    right_reference_effective = entries[run_end][1]
                    if _positive_interval_overlap(
                        left_reference_effective,
                        right_reference_effective,
                    ):
                        chosen = _choose_final_deletion_anchors(
                            anchors.get(left_ref.allele, ()),
                            anchors.get(right_ref.allele, ()),
                            max_extension,
                        )
                status = "DEL" if chosen is not None else "NA"
                missing_names = tuple(row.allele for row, _ in missing)
                missing_original_loci = tuple(row.fields[2] for row, _ in missing)
                joined_names = ";".join(missing_names)
                reference_loc = f"{chrom}:{ref_start}-{ref_end}+"
                left_anchor = chosen[0].allele if chosen is not None else ""
                right_anchor = chosen[1].allele if chosen is not None else ""
                output.append("\t".join((
                    status,
                    joined_names,
                    ".",
                    ";".join(missing_original_loci),
                    query_genome,
                    ".",
                    status,
                    joined_names,
                    joined_names,
                    reference_loc,
                    left_anchor,
                    right_anchor,
                    "0",
                    "0",
                    "0",
                )))
    return output


def append_final_unmatched_rows(
    path: str,
    refhaplo: str,
    max_extension: int = 10000,
) -> Tuple[int, int]:
    unmatched_rows = build_final_unmatched_rows(
        read_final_lift_rows(path), refhaplo, max_extension,
    )
    if unmatched_rows:
        with open(path, "a") as output:
            for row in unmatched_rows:
                output.write(row + "\n")
    deletion_count = sum(row.startswith("DEL\t") for row in unmatched_rows)
    unknown_count = sum(row.startswith("NA\t") for row in unmatched_rows)
    return deletion_count, unknown_count


def end_sign(side: str, strand: str):
    if strand == "+":
        return "-" if side == "left" else "+"
    else:
        return "+" if side == "left" else "-"


def normalize_tag(tag: str, ignore_distance=False):
    if not tag:
        return tag
    if ignore_distance:
        return re.sub(r"(?:_distance\d+)+$", "", tag)
    return tag


def extract_signed_chunks_from_tag(tag: str):
    tag = normalize_tag(tag, ignore_distance=True).strip()
    out = []
    for m in re.finditer(r"((?:[^_]+_){3}[^_]+)([+\-I])", tag):
        out.append((m.group(1), m.group(2)))
    return out


def best_location_ref_for_delegate(refa: RefAllele, ref_alleles: Dict[str, RefAllele]):
    if refa.gene_state and refa.gene_state in ref_alleles:
        gs = ref_alleles[refa.gene_state]
        if not gs.delegate_only:
            return refa.gene_state

    candidates = []
    for tag in (refa.left_tag, refa.right_tag):
        if not tag:
            continue
        raw = normalize_tag(tag, ignore_distance=False)
        candidates.append(strip_terminal_gene_state_suffix(raw, refa.gene_state))
        candidates.append(normalize_tag(raw, ignore_distance=True))

    seen = set()
    for cand in candidates:
        cand = cand.strip()
        if not cand or cand in seen:
            continue
        seen.add(cand)

        exact = re.sub(r"[+\-I]$", "", cand)
        if exact in ref_alleles and not ref_alleles[exact].delegate_only:
            return exact

        for name, sign in extract_signed_chunks_from_tag(cand):
            if name in ref_alleles and not ref_alleles[name].delegate_only:
                return name

    return None


def location_refname_for_output(refname: str, ref_alleles: Dict[str, RefAllele]):
    refa = ref_alleles.get(refname)
    if refa is None:
        return refname

    if not refa.delegate_only:
        return refname

    resolved = best_location_ref_for_delegate(refa, ref_alleles)
    return resolved if resolved is not None else ""


def location_coords_for_output(refname: str, ref_alleles: Dict[str, RefAllele]):
    location_refname = location_refname_for_output(refname, ref_alleles)
    if not location_refname:
        return ""

    ra = ref_alleles.get(location_refname)
    if ra is None:
        return ""

    return fixed_ref_coords(ra) or original_ref_coords(ra)


def fixed_ref_coords(ra: RefAllele):
    if ra.sub_start >= 0 and ra.sub_end >= 0:
        return f"{ra.contig}:{ra.sub_start}-{ra.sub_end}"
    return ""


def original_ref_coords(ra: RefAllele):
    return f"{ra.contig}:{ra.start}-{ra.end}"


def format_ref_anchor_interval(contig: str, start: int, end: int, strand: str = ""):
    if strand in {"+", "-"}:
        return f"{contig}:{start}-{end}{strand}"
    return f"{contig}:{start}-{end}"


def build_reference_anchor_index(ref_alleles: Dict[str, RefAllele]):
    anchor_by_seed = {}
    for refname, refa in ref_alleles.items():
        left_seed, right_seed = seed_tags_for_reference(refname, refa)
        if left_seed:
            anchor_by_seed[left_seed] = AnchorHit(
                refname=refname,
                contig=refa.contig,
                pos=refa.start,
                strand=refa.strand,
                side="left",
            )
        if right_seed:
            anchor_by_seed[right_seed] = AnchorHit(
                refname=refname,
                contig=refa.contig,
                pos=refa.end,
                strand=refa.strand,
                side="right",
            )
    ordered_seeds = sorted(anchor_by_seed.keys(), key=len, reverse=True)
    return anchor_by_seed, ordered_seeds


def resolve_reference_anchor_from_tag(tag: str, anchor_by_seed, ordered_seeds):
    raw = normalize_tag(tag, ignore_distance=True).strip()
    if not raw:
        return None

    hit = anchor_by_seed.get(raw)
    if hit is not None:
        return hit

    for seed in ordered_seeds:
        if raw == seed or raw.startswith(seed + "_"):
            return anchor_by_seed[seed]

    return None


def infer_ref_coords_from_asm_tags(left_tag: str, right_tag: str, anchor_by_seed, ordered_seeds):
    left_hit = resolve_reference_anchor_from_tag(left_tag, anchor_by_seed, ordered_seeds)
    right_hit = resolve_reference_anchor_from_tag(right_tag, anchor_by_seed, ordered_seeds)

    if left_hit is None or right_hit is None:
        return ""
    if left_hit.contig != right_hit.contig:
        return ""

    start = min(left_hit.pos, right_hit.pos)
    end = max(left_hit.pos, right_hit.pos)
    strand = left_hit.strand if left_hit.strand == right_hit.strand else ""
    return format_ref_anchor_interval(left_hit.contig, start, end, strand)


def contig_sort_key(contig: str):
    """
    Natural-ish contig ordering for assembly-coordinate output sorting.

    Canonical chr1..chr22, chrX, chrY, chrM are ordered biologically.
    Other contigs are ordered after canonical contigs using a tokenized
    alphanumeric key, so contig2 sorts before contig10.
    """
    canonical_order = {f"chr{i}": i for i in range(1, 23)}
    canonical_order.update({"chrX": 23, "chrY": 24, "chrM": 25, "chrMT": 25})

    if contig in canonical_order:
        return (0, canonical_order[contig])

    m = re.match(r"^chr(\d+)$", contig)
    if m:
        return (0, int(m.group(1)))

    tokens = []
    for tok in re.split(r"(\d+)", contig):
        if not tok:
            continue
        if tok.isdigit():
            tokens.append((0, int(tok)))
        else:
            tokens.append((1, tok))

    return (1, tokens)


def allele_record_assembly_sort_key(r):
    """Sort rows within each genome/haplotype by assembly coordinates."""
    return (
        r.genome,
        contig_sort_key(r.asm_contig),
        r.asm_start,
        r.asm_end,
        r.allelename,
    )


def direct_ref_seed_tags(refname: str, refa: RefAllele):
    return (
        f"{refname}{end_sign('left', refa.strand)}",
        f"{refname}{end_sign('right', refa.strand)}",
    )


def strip_terminal_gene_state_suffix(tag: str, gene_state: str):
    if not tag or not gene_state:
        return tag
    return re.sub(rf"_{re.escape(gene_state)}[+-]$", "", tag)


def synthesize_delegate_seed_tags(refa: RefAllele):
    candidates = []
    for tag in (refa.left_tag, refa.right_tag):
        tag = normalize_tag(tag, ignore_distance=False)
        if not tag:
            continue
        candidates.append(strip_terminal_gene_state_suffix(tag, refa.gene_state))

    if not candidates:
        return ("", "")

    base = min(candidates, key=lambda x: (len(x), x))
    return (
        f"{base}_{refa.gene_state}{end_sign('left', refa.strand)}",
        f"{base}_{refa.gene_state}{end_sign('right', refa.strand)}",
    )


def seed_tags_for_reference(refname: str, refa: RefAllele):
    if refa.delegate_only:
        return synthesize_delegate_seed_tags(refa)
    return direct_ref_seed_tags(refname, refa)


def finalize_delegate_reference_tags(ref_alleles: Dict[str, RefAllele]):
    for ra in ref_alleles.values():
        if not ra.delegate_only:
            continue
        left, right = synthesize_delegate_seed_tags(ra)
        if left and right:
            ra.left_tag = left
            ra.right_tag = right


def eligible_ref_alleles(ref_alleles: Dict[str, RefAllele]):
    return {
        name: ra
        for name, ra in ref_alleles.items()
        if ra.unique_len > 0 and allele_belongs_to_genome(name, ra.genome)
    }


def iter_rows(path):
    with open_maybe_gzip(path) as f:
        for line in f:
            line = line.rstrip("\n")
            if not line or line.startswith("#"):
                continue
            parts = line.split("\t")
            if len(parts) < 8:
                continue
            yield parts


def build_prefix_priority_index(genome_files, refhaplo: str, filter_hg38_alts: bool = True):
    """Build the historical sample/RefMatch priority ranks.

    This helper is retained for compatibility with older callers and tests.
    The production ``main`` path now uses
    :func:`build_reference_size_priority_index`, so assembly lengths and
    RefMatch distances no longer affect current GenomeLift ownership.

    For every row matched to an allele from refhaplo, compute:

      reference_score = reference_size / max(RefMatch_distance, 1e-6)

    For each allele-matrix prefix, retain:
      total_reference_score = sum of reference_score among its supported rows
      total_allele_size = total size of that prefix across all split genomes

    A prefix without a reference-supported row uses total_allele_size as its
    primary score. Prefixes are sorted by (primary_score, total_allele_size,
    prefix), then converted to compact integer ranks. Larger rank means higher
    overlap priority.
    """
    total_reference_score = cl.defaultdict(float)
    prefixes_with_reference = set()
    total_allele_size = cl.Counter()

    for _, genome_file in genome_files.items():
        for parts in iter_rows(genome_file):
            allelename = parts[0]
            best_ref_allele = parts[1]
            assembly_loc = parts[2]
            reference_loc = parts[3]
            distance_raw = parts[5]

            try:
                asm_contig, asm_start, asm_end, _ = parse_loc(assembly_loc)
                if filter_hg38_alts and is_hg38_alt_assembly(allelename, asm_contig):
                    continue
                genome = determine_assembly_genome(allelename, asm_contig)
                prefix = allele_matrix_name(allelename)
            except Exception:
                continue

            size = max(0, asm_end - asm_start)
            total_allele_size[prefix] += size

            try:
                _, _, _, _, matched_genome = parse_allele_name(best_ref_allele)
                _, ref_start, ref_end, _ = parse_loc(reference_loc)
                distance = float(distance_raw)
            except (TypeError, ValueError):
                continue

            if matched_genome != refhaplo:
                continue

            if not math.isfinite(distance):
                distance = math.inf

            reference_size = max(0, ref_end - ref_start)
            total_reference_score[prefix] += (
                reference_size / max(distance, 1e-6)
            )
            prefixes_with_reference.add(prefix)

    def primary_score(prefix):
        if prefix in prefixes_with_reference:
            return total_reference_score[prefix]
        return float(total_allele_size[prefix])

    ordered = sorted(
        total_allele_size.keys(),
        key=lambda prefix: (
            primary_score(prefix),
            total_allele_size[prefix],
            prefix,
        ),
    )
    return {prefix: rank + 1 for rank, prefix in enumerate(ordered)}


def build_reference_size_priority_index(ref_alleles):
    """Build stable matrix priorities from the reference haplotype only.

    Assembly intervals and KmerMatch distances can vary between samples even
    when the same graph block is being lifted.  They therefore must not decide
    which overlapping assembly interval owns a region.  The reference
    haplotype is the shared coordinate system: rank each matrix by the total
    size of its CHM13/reference blocks, with the matrix name as the
    deterministic tie-breaker.  Larger ranks retain the existing
    ``collapse_unique_regions`` convention of higher priority.

    ``ref_alleles`` should be the catalog collected from the requested
    reference haplotype *before* promoted alleles from other genomes are
    added.  Alleles with no positive reference span are omitted and therefore
    receive the normal unknown-prefix priority of zero.
    """
    total_reference_size = cl.Counter()
    for name, ref in (ref_alleles or {}).items():
        try:
            size = max(0, int(ref.end) - int(ref.start))
        except (AttributeError, TypeError, ValueError):
            continue
        if size <= 0:
            continue
        total_reference_size[allele_matrix_name(name)] += size

    ordered = sorted(
        total_reference_size,
        key=lambda prefix: (total_reference_size[prefix], prefix),
    )
    return {prefix: rank + 1 for rank, prefix in enumerate(ordered)}


def assign_reference_unique_regions(ref_alleles: Dict[str, RefAllele], priority_by_prefix=None):
    by_gc = cl.defaultdict(list)
    for ra in ref_alleles.values():
        by_gc[(ra.genome, ra.contig)].append((ra.start, ra.end, ra.name))

    for _, regions in by_gc.items():
        regions.sort(key=lambda x: (x[0], x[1], x[2]))
        unique_map = collapse_unique_regions(regions, priority_by_prefix=priority_by_prefix)
        for name, vals in unique_map.items():
            sub_start, sub_end, unique_len, left_off, right_off = vals
            ra = ref_alleles[name]
            ra.sub_start = sub_start
            ra.sub_end = sub_end
            ra.unique_len = unique_len
            ra.left_offset = left_off
            ra.right_offset = right_off


def build_gene_states_from_edges(ref_alleles: Dict[str, RefAllele], ref_edges):
    dsu = DSU()
    ref_names = set(ref_alleles.keys())
    for name in ref_names:
        dsu.add(name)

    for a, others in ref_edges.items():
        if a not in ref_names:
            continue
        for b in others:
            if b in ref_names:
                dsu.union(a, b)

    groups = cl.defaultdict(list)
    for name in ref_names:
        groups[dsu.find(name)].append(name)

    for _, members in groups.items():
        members.sort()
        label = members[0]
        for m in members:
            ref_alleles[m].gene_state = label

    return groups


def initialize_reference_end_tags(ref_alleles: Dict[str, RefAllele]):
    gene_state_count = cl.Counter(ra.gene_state for ra in ref_alleles.values())
    ends = {}

    for ra in ref_alleles.values():
        left = EndRec(
            allele=ra.name, genome=ra.genome, contig=ra.contig,
            pos=ra.start, side="left", strand=ra.strand,
            gene_state=ra.gene_state, tag=""
        )
        right = EndRec(
            allele=ra.name, genome=ra.genome, contig=ra.contig,
            pos=ra.end, side="right", strand=ra.strand,
            gene_state=ra.gene_state, tag=""
        )

        if ra.unique_len > 0 and gene_state_count[ra.gene_state] == 1 and not ra.delegate_only:
            left.tag = f"{ra.name}{end_sign('left', ra.strand)}"
            right.tag = f"{ra.name}{end_sign('right', ra.strand)}"
            ra.tag_fixed = True
        else:
            ra.tag_fixed = False

        ends[left.key()] = left
        ends[right.key()] = right

    return ends


def propagate_reference_tags(ref_alleles: Dict[str, RefAllele], ref_ends: Dict[Tuple[str, str], EndRec], max_near=1000):
    import heapq

    by_gc = cl.defaultdict(list)
    for end in ref_ends.values():
        by_gc[(end.genome, end.contig)].append(end)

    partner = {}
    for ra in ref_alleles.values():
        partner[(ra.name, "left")] = (ra.name, "right")
        partner[(ra.name, "right")] = (ra.name, "left")

    for _, arr in by_gc.items():
        arr.sort(key=lambda e: (e.pos, e.allele, e.side))
        n = len(arr)
        idx_of = {e.key(): i for i, e in enumerate(arr)}
        best_dist = {e.key(): math.inf for e in arr}
        heap = []
        push_id = 0

        def push(end_key, dist, tag):
            nonlocal push_id
            if dist < best_dist[end_key]:
                best_dist[end_key] = dist
                heapq.heappush(heap, (dist, push_id, end_key, tag))
                push_id += 1

        for e in arr:
            if e.tag:
                push(e.key(), 0, e.tag)

        if not heap:
            continue

        while heap:
            dist, _, ekey, tag = heapq.heappop(heap)
            if dist != best_dist[ekey]:
                continue

            e = ref_ends[ekey]
            if not e.tag:
                e.tag = tag

            pkey = partner[ekey]
            pe = ref_ends[pkey]
            if not pe.tag:
                ra = ref_alleles[e.allele]
                # An unfixed reference block has no tag of its own. It is a
                # transparent interval: preserve the incoming label verbatim
                # and charge the physical block length while crossing it.
                push(pkey, dist + max(0, ra.end - ra.start), tag)

            i = idx_of[ekey]

            if i > 0:
                nb = arr[i - 1]
                if not nb.tag:
                    nd = dist + (e.pos - nb.pos)
                    if nd <= max_near:
                        ntag = tag
                    else:
                        ntag = f"{normalize_tag(tag, ignore_distance=False)}_distance{nd}"
                    push(nb.key(), nd, ntag)

            if i + 1 < n:
                nb = arr[i + 1]
                if not nb.tag:
                    nd = dist + (nb.pos - e.pos)
                    if nd <= max_near:
                        ntag = tag
                    else:
                        ntag = f"{normalize_tag(tag, ignore_distance=False)}_distance{nd}"
                    push(nb.key(), nd, ntag)

    for ra in ref_alleles.values():
        ra.left_tag = ref_ends[(ra.name, "left")].tag
        ra.right_tag = ref_ends[(ra.name, "right")].tag


def assign_reference_propagated_tags(
    ref_alleles: Dict[str, RefAllele],
    ref_ends: Dict[Tuple[str, str], EndRec],
    max_near=1000,
):
    """
    Build a second, reference-only propagated tag pair for every reference allele.

    left_tag/right_tag remain the initiated/normal reference tags. For unique
    reference anchors, those are identity tags and they still seed propagation.

    propagated_left_tag/propagated_right_tag record what each reference endpoint
    receives from propagation from other seeded reference alleles. Initiated tags
    do not block this propagated container. These propagated tags are for
    comparison only and are never used to seed assembly propagation.
    """
    import heapq

    by_gc = cl.defaultdict(list)
    for end in ref_ends.values():
        by_gc[(end.genome, end.contig)].append(end)

    idx_lookup = {}
    for arr in by_gc.values():
        arr.sort(key=lambda e: (e.pos, e.allele, e.side))
        for i, end in enumerate(arr):
            idx_lookup[end.key()] = (arr, i)

    partner = {}
    for ra in ref_alleles.values():
        partner[(ra.name, "left")] = (ra.name, "right")
        partner[(ra.name, "right")] = (ra.name, "left")

    def make_partner_tag(tag, ra, side):
        return f"{normalize_tag(tag, ignore_distance=False)}_{ra.gene_state}{end_sign(side, ra.strand)}"

    # For each endpoint, keep the two best propagated states from distinct
    # source alleles. That is enough to recover the best non-self source for
    # the endpoint even when the closest source is the endpoint's own allele.
    best_states = cl.defaultdict(list)
    heap = []
    push_id = 0

    def state_sort_key(state):
        dist, source_allele, tag = state
        return (dist, tag, source_allele)

    def push(end_key, dist, tag, source_allele):
        nonlocal push_id
        states = best_states[end_key]

        replaced = False
        for i, (old_dist, old_source, old_tag) in enumerate(states):
            if old_source == source_allele:
                if (dist, tag) < (old_dist, old_tag):
                    states[i] = (dist, source_allele, tag)
                    replaced = True
                break
        else:
            if len(states) < 2:
                states.append((dist, source_allele, tag))
                replaced = True
            else:
                worst_i = max(range(len(states)), key=lambda j: state_sort_key(states[j]))
                if state_sort_key((dist, source_allele, tag)) < state_sort_key(states[worst_i]):
                    states[worst_i] = (dist, source_allele, tag)
                    replaced = True

        if not replaced:
            return

        states.sort(key=state_sort_key)
        heapq.heappush(heap, (dist, push_id, end_key, tag, source_allele))
        push_id += 1

    for e in ref_ends.values():
        # Only genuinely fixed reference alleles emit propagation seeds.
        # Labels inherited by an unresolved reference block must not be
        # reclassified as that block's own source in this secondary index.
        if ref_alleles[e.allele].tag_fixed and e.tag:
            push(e.key(), 0, e.tag, e.allele)

    while heap:
        dist, _, ekey, tag, source_allele = heapq.heappop(heap)
        states = best_states.get(ekey, [])
        if not any(
            d == dist and src == source_allele and t == tag
            for d, src, t in states
        ):
            continue

        e = ref_ends[ekey]

        pkey = partner[ekey]
        pe = ref_ends[pkey]
        ra = ref_alleles[e.allele]
        if ra.tag_fixed:
            ptag = make_partner_tag(tag, ra, pe.side)
            partner_dist = dist
        else:
            ptag = tag
            partner_dist = dist + max(0, ra.end - ra.start)
        push(pkey, partner_dist, ptag, source_allele)

        arr, i = idx_lookup[ekey]

        if i > 0:
            nb = arr[i - 1]
            nd = dist + (e.pos - nb.pos)
            if nd <= max_near:
                ntag = tag
            else:
                ntag = f"{normalize_tag(tag, ignore_distance=False)}_distance{nd}"
            push(nb.key(), nd, ntag, source_allele)

        if i + 1 < len(arr):
            nb = arr[i + 1]
            nd = dist + (nb.pos - e.pos)
            if nd <= max_near:
                ntag = tag
            else:
                ntag = f"{normalize_tag(tag, ignore_distance=False)}_distance{nd}"
            push(nb.key(), nd, ntag, source_allele)

    for ra in ref_alleles.values():
        left_key = (ra.name, "left")
        right_key = (ra.name, "right")

        left_external = next(
            (tag for _, source_allele, tag in best_states.get(left_key, [])
             if source_allele != ra.name),
            "",
        )
        right_external = next(
            (tag for _, source_allele, tag in best_states.get(right_key, [])
             if source_allele != ra.name),
            "",
        )

        ra.propagated_left_tag = left_external
        ra.propagated_right_tag = right_external

    # Also retain one coherent pair propagated from each directly adjacent
    # reference block. Normal propagation above keeps only the globally best
    # external state and can therefore discard the label from one physical
    # side, or choose a longer same-source composite over the direct edge.
    # These pairs are comparison-only; official tags and propagation seeds are
    # unchanged.
    ordered_by_gc = cl.defaultdict(list)
    for ra in ref_alleles.values():
        if ra.unique_len > 0:
            ordered_by_gc[(ra.genome, ra.contig)].append(ra)

    def with_distance(tag, distance):
        if not tag:
            return ""
        if distance <= max_near:
            return tag
        return f"{normalize_tag(tag, ignore_distance=False)}_distance{distance}"

    for arr in ordered_by_gc.values():
        arr.sort(key=lambda item: (item.start, item.end, item.name))
        for index, ra in enumerate(arr):
            pairs = []

            if index > 0:
                upstream = arr[index - 1]
                upstream_seed = with_distance(
                    upstream.right_tag,
                    max(0, ra.start - upstream.end),
                )
                if upstream_seed:
                    upstream_partner = (
                        make_partner_tag(upstream_seed, ra, "right")
                        if ra.tag_fixed else upstream_seed
                    )
                    pairs.append((
                        upstream_seed,
                        upstream_partner,
                    ))

            if index + 1 < len(arr):
                downstream = arr[index + 1]
                downstream_seed = with_distance(
                    downstream.left_tag,
                    max(0, downstream.start - ra.end),
                )
                if downstream_seed:
                    downstream_partner = (
                        make_partner_tag(downstream_seed, ra, "left")
                        if ra.tag_fixed else downstream_seed
                    )
                    pairs.append((
                        downstream_partner,
                        downstream_seed,
                    ))

            # Preserve order (upstream first, downstream second) while
            # removing an exact duplicate at coincident/equivalent edges.
            ra.adjacent_tag_pairs = tuple(dict.fromkeys(pairs))


def split_input_by_genome(input_path: str, outdir: str, filter_hg38_alts: bool = True):
    os.makedirs(outdir, exist_ok=True)
    handles = {}
    genome_files = {}

    try:
        for parts in iter_rows(input_path):
            allelename = parts[0]
            assembly_loc = parts[2]

            try:
                asm_contig, _, _, _ = parse_loc(assembly_loc)

                if filter_hg38_alts and is_hg38_alt_assembly(allelename, asm_contig):
                    continue

                genome = determine_assembly_genome(allelename, asm_contig)
            except Exception:
                continue

            if genome not in handles:
                path = os.path.join(outdir, f"{genome}.tsv")
                handles[genome] = open(path, "w")
                genome_files[genome] = path

            handles[genome].write("\t".join(parts[:8]) + "\n")
    finally:
        for h in handles.values():
            h.close()

    return genome_files


def collect_reference_from_refgenome_file(ref_file: str, refhaplo: str, gene_state_cutoff: float):
    """
    Build reference alleles and the reference-only similarity graph.

    Every block in the reference-haplotype file defines itself from its
    assembly coordinates.  Its RefMatch-selected best hit is additional
    similarity information, not the authority for whether that reference
    block exists.  This distinction matters for masked/no-kmer blocks and for
    tied paralogs whose selected best hit is a different reference block.

    gene_state edges are defined ONLY by the reference-haplotype self rows,
    and ONLY when reported difference < gene_state_cutoff.
    """
    ref_alleles = {}
    ref_edges = cl.defaultdict(set)

    for parts in iter_rows(ref_file):
        allelename = parts[0]
        best_ref_allele = parts[1]
        assembly_loc = parts[2]
        ref_loc = parts[3]
        class_type = parts[6].strip().lower()
        similar_refs_raw = parts[7]

        # The row itself is a reference allele even if KmerMatch reported no
        # usable hit (best_ref_allele is NA) or selected a tied paralog.  Use
        # the reference assembly interval carried by this row as its canonical
        # coordinates.
        if allele_belongs_to_genome(allelename, refhaplo):
            try:
                self_contig, self_start, self_end, self_strand = parse_loc(
                    assembly_loc,
                )
            except Exception:
                pass
            else:
                ref_alleles[allelename] = RefAllele(
                    name=allelename,
                    genome=refhaplo,
                    contig=self_contig,
                    start=self_start,
                    end=self_end,
                    strand=self_strand,
                    delegate_only=False,
                )

        # Column 2 must identify an allele from the requested reference
        # haplotype.  A malformed/mixed RefMatch file must not seed query
        # alleles into the reference catalog.
        if not allele_belongs_to_genome(best_ref_allele, refhaplo):
            continue

        try:
            ref_contig, ref_start, ref_end, ref_strand = parse_loc(ref_loc)
        except Exception:
            continue

        if best_ref_allele not in ref_alleles:
            ref_alleles[best_ref_allele] = RefAllele(
                name=best_ref_allele,
                genome=refhaplo,
                contig=ref_contig,
                start=ref_start,
                end=ref_end,
                strand=ref_strand,
                delegate_only=(class_type == "novel"),
            )
        elif class_type == "novel":
            ref_alleles[best_ref_allele].delegate_only = True

        for other, diff in split_similar_refs(similar_refs_raw):
            if (
                diff < gene_state_cutoff
                and allele_belongs_to_genome(other, refhaplo)
            ):
                ref_edges[best_ref_allele].add(other)

    return ref_alleles, ref_edges


def supplement_reference_alleles_from_all_genomes(genome_files, ref_alleles, refhaplo):
    """
    Promote best_ref_allele values observed outside the refhaplo file into the
    reference-haplotype allele set, as long as they have a usable ref_loc.

    IMPORTANT:
    This function does NOT add any similarity edges to ref_edges.
    Non-reference genomes may help supply coordinates / delegate_only flags,
    but they cannot connect reference alleles into the same gene_state.
    """
    coord_support = cl.defaultdict(cl.Counter)
    coord_best_similarity = {}
    delegate_only_names = set()
    promoted = 0

    for _, genome_file in genome_files.items():
        for parts in iter_rows(genome_file):
            best_ref_allele = parts[1].strip()
            ref_loc = parts[3].strip()
            similarity_raw = parts[5].strip()
            class_type = parts[6].strip().lower()

            if not allele_belongs_to_genome(best_ref_allele, refhaplo):
                continue

            if class_type == "novel":
                delegate_only_names.add(best_ref_allele)

            try:
                ref_contig, ref_start, ref_end, ref_strand = parse_loc(ref_loc)
            except Exception:
                continue

            coord_key = (ref_contig, ref_start, ref_end, ref_strand)
            coord_support[best_ref_allele][coord_key] += 1

            try:
                similarity = float(similarity_raw)
            except ValueError:
                similarity = math.inf

            score_key = (best_ref_allele, coord_key)
            if similarity < coord_best_similarity.get(score_key, math.inf):
                coord_best_similarity[score_key] = similarity

    for allele_name in delegate_only_names:
        if allele_name in ref_alleles:
            ref_alleles[allele_name].delegate_only = True

    for allele_name, counts in coord_support.items():
        if allele_name in ref_alleles or not counts:
            continue

        best_coord = min(
            counts.items(),
            key=lambda kv: (
                -kv[1],
                coord_best_similarity.get((allele_name, kv[0]), math.inf),
                kv[0][0], kv[0][1], kv[0][2], kv[0][3],
            ),
        )[0]
        ref_contig, ref_start, ref_end, ref_strand = best_coord

        ref_alleles[allele_name] = RefAllele(
            name=allele_name,
            genome=refhaplo,
            contig=ref_contig,
            start=ref_start,
            end=ref_end,
            strand=ref_strand,
            delegate_only=(allele_name in delegate_only_names),
        )
        promoted += 1

    return promoted


def load_genome_records(genome_file: str, ref_alleles: Dict[str, RefAllele], filter_hg38_alts: bool = True):
    records = []
    for parts in iter_rows(genome_file):
        allelename = parts[0]
        best_ref_allele = parts[1]
        assembly_loc = parts[2]
        ref_loc = parts[3]
        intron_or_exon = parts[4]
        similarity_raw = parts[5]
        try:
            similarity = float(similarity_raw)
        except ValueError:
            similarity = math.inf
        class_type = parts[6]
        similar_refs_raw = parts[7]

        try:
            asm_contig, asm_start, asm_end, asm_strand = parse_loc(assembly_loc)

            if filter_hg38_alts and is_hg38_alt_assembly(allelename, asm_contig):
                continue

            genome = determine_assembly_genome(allelename, asm_contig)
        except Exception:
            continue

        rec = AlleleRecord(
            allelename=allelename,
            best_ref_allele=best_ref_allele,
            assembly_loc=assembly_loc,
            ref_loc=ref_loc,
            intron_or_exon=intron_or_exon,
            similarity=similarity,
            similarity_raw=similarity_raw,
            class_type=class_type,
            similar_refs_raw=similar_refs_raw,
            genome=genome,
            asm_contig=asm_contig,
            asm_start=asm_start,
            asm_end=asm_end,
            asm_strand=asm_strand,
            parsed_similar_refs=tuple(split_similar_refs(similar_refs_raw)),
        )

        ra = ref_alleles.get(best_ref_allele)
        rec.gene_state = ra.gene_state if ra is not None else best_ref_allele
        records.append(rec)

    return records


def assign_assembly_unique_regions(
    rows: List[AlleleRecord], alignment_scores=None,
):
    by_gc = cl.defaultdict(list)
    row_map = {}

    for r in rows:
        by_gc[(r.genome, r.asm_contig)].append((r.asm_start, r.asm_end, r.allelename))
        row_map[r.allelename] = r

    priority_by_prefix = GLOBAL_PRIORITY_BY_PREFIX or {}
    if alignment_scores is None:
        alignment_scores = GLOBAL_ALIGNMENT_SCORES

    for _, regions in by_gc.items():
        regions.sort(key=lambda x: (x[0], x[1], x[2]))
        unique_map = collapse_unique_regions(
            regions,
            priority_by_prefix=priority_by_prefix,
            priority_by_allele=alignment_scores,
        )
        display_offsets = compute_priority_display_offsets_for_regions(regions, unique_map)

        for allelename, vals in unique_map.items():
            sub_start, sub_end, unique_len, left_off, right_off = vals
            r = row_map[allelename]
            r.allele_unique_start = sub_start
            r.allele_unique_end = sub_end
            r.allele_unique_len = unique_len
            if unique_len > 0:
                r.allele_left_offset, r.allele_right_offset = display_offsets.get(
                    allelename,
                    (str(left_off), str(right_off)),
                )
            else:
                r.allele_left_offset = NEG_INF_STR
                r.allele_right_offset = NEG_INF_STR


def confident_one_to_one(rows: List[AlleleRecord], ref_alleles: Dict[str, RefAllele], sim_cutoff=0.1):
    all_ref_alleles = ref_alleles
    if ref_alleles is GLOBAL_REF_ALLELES and GLOBAL_ASSIGNABLE_REF_ALLELES is not None:
        assignable_ref_alleles = GLOBAL_ASSIGNABLE_REF_ALLELES
    else:
        assignable_ref_alleles = eligible_ref_alleles(ref_alleles)

    anchors = {}
    used_ref_keys = set()

    # Coordinate-identical blocks from the same matrix are stronger evidence
    # than KmerMatch ties.  This is especially important for reference
    # bootstraps: several paralogous blocks can all have distance 0, while the
    # query/reference block at the identical source interval is unambiguous.
    # Keep the same mutual-uniqueness rule used by the normal candidate pass so
    # duplicated query intervals cannot both claim one reference block.
    exact_refs_by_key = cl.defaultdict(list)
    for refname, refa in assignable_ref_alleles.items():
        exact_key = (
            allele_matrix_name(refname),
            refa.contig,
            refa.start,
            refa.end,
            refa.strand,
        )
        exact_refs_by_key[exact_key].append(refname)

    exact_claims = cl.defaultdict(list)
    row_to_exact_ref = {}
    for r in rows:
        exact_key = (
            allele_matrix_name(r.allelename),
            r.asm_contig,
            r.asm_start,
            r.asm_end,
            r.asm_strand,
        )
        exact_refs = sorted(set(exact_refs_by_key.get(exact_key, ())))
        if len(exact_refs) != 1:
            continue
        refname = exact_refs[0]
        refkey = (r.asm_contig, refname)
        exact_claims[refkey].append(r.allelename)
        row_to_exact_ref[r.allelename] = refname

    for r in rows:
        refname = row_to_exact_ref.get(r.allelename)
        if refname is None:
            continue
        refkey = (r.asm_contig, refname)
        if len(exact_claims[refkey]) == 1 and refkey not in used_ref_keys:
            anchors[r.allelename] = refname
            used_ref_keys.add(refkey)

    candidate_claims = cl.defaultdict(list)
    row_to_unique_ref = {}

    for r in rows:
        if r.allelename in anchors:
            continue
        if r.class_type.lower().endswith("kmerrescue"):
            # match_partition_blocks emits this explicit marker only when
            # ordinary RefMatch retained no hit and one reference had >5000
            # weighted common k-mers and >10x the runner-up. Keep the genuine
            # normalized distance but allow that separately qualified evidence
            # to enter the same mutual one-to-one claim resolution.
            refs_filt = [
                name for name, _val in r.parsed_similar_refs
                if name in all_ref_alleles
            ]
        else:
            refs_filt = [
                name for name, val in r.parsed_similar_refs
                if val < sim_cutoff and name in all_ref_alleles
            ]
        refs_filt = sorted(set(refs_filt))
        if len(refs_filt) == 1:
            refname = refs_filt[0]
            refkey = (r.asm_contig, refname)
            candidate_claims[refkey].append(r.allelename)
            row_to_unique_ref[r.allelename] = refname

    for r in rows:
        if r.allelename in anchors:
            continue
        refname = row_to_unique_ref.get(r.allelename)
        if refname is None:
            continue
        if refname not in assignable_ref_alleles:
            continue
        refkey = (r.asm_contig, refname)
        if len(candidate_claims[refkey]) == 1 and refkey not in used_ref_keys:
            anchors[r.allelename] = refname
            used_ref_keys.add(refkey)

    return anchors


def initialize_assembly_tags(
    rows,
    ref_alleles,
    anchors,
    non_reference_tags=frozenset(),
):
    ends = {}
    for r in rows:
        left = EndRec(
            allele=r.allelename, genome=r.genome, contig=r.asm_contig,
            pos=r.asm_start, side="left", strand=r.asm_strand,
            gene_state=r.gene_state, tag=""
        )
        right = EndRec(
            allele=r.allelename, genome=r.genome, contig=r.asm_contig,
            pos=r.asm_end, side="right", strand=r.asm_strand,
            gene_state=r.gene_state, tag=""
        )

        if r.allelename in anchors:
            refname = anchors[r.allelename]
            refa = ref_alleles[refname]
            seed_left, seed_right = seed_tags_for_reference(refname, refa)
            same_orientation = (r.asm_strand == refa.strand)
            if same_orientation:
                left.tag = seed_left
                right.tag = seed_right
            else:
                left.tag = seed_right
                right.tag = seed_left
        elif non_reference_tag_enabled(
            r.allelename, ref_alleles, non_reference_tags,
        ):
            # Explicitly enabled non-reference alleles are fixed synthetic
            # anchors and emit their exact own tags. They never inherit or
            # rewrite a neighboring reference label.
            left.tag = f"{r.allelename}{end_sign('left', r.asm_strand)}"
            right.tag = f"{r.allelename}{end_sign('right', r.asm_strand)}"

        ends[(r.allelename, "left")] = left
        ends[(r.allelename, "right")] = right

    return ends


def propagate_assembly_tags(rows, asm_ends, max_near=1000):
    import heapq

    by_contig = cl.defaultdict(list)
    row_by_name = {r.allelename: r for r in rows}

    for end in asm_ends.values():
        by_contig[end.contig].append(end)
    for arr in by_contig.values():
        arr.sort(key=lambda e: (e.pos, e.allele, e.side))

    partner = {}
    for r in rows:
        partner[(r.allelename, "left")] = (r.allelename, "right")
        partner[(r.allelename, "right")] = (r.allelename, "left")

    for _, arr in by_contig.items():
        idx_of = {e.key(): i for i, e in enumerate(arr)}
        best_dist = {e.key(): math.inf for e in arr}
        heap = []
        push_id = 0

        def push(end_key, dist, tag):
            nonlocal push_id
            if dist < best_dist[end_key]:
                best_dist[end_key] = dist
                heapq.heappush(heap, (dist, push_id, end_key, tag))
                push_id += 1

        for e in arr:
            if e.tag:
                push(e.key(), 0, e.tag)

        if not heap:
            continue

        while heap:
            dist, _, ekey, tag = heapq.heappop(heap)
            if dist != best_dist[ekey]:
                continue

            e = asm_ends[ekey]
            if not e.tag:
                e.tag = tag

            pkey = partner[ekey]
            pe = asm_ends[pkey]
            if not pe.tag:
                # Unfixed query blocks are transparent. Do not invent a tag
                # from their allele name, and do not erase the incoming tag.
                # Charging the interval length allows opposite fixed flanks to
                # control opposite ends of a consecutive unresolved run.
                row = row_by_name[e.allele]
                push(
                    pkey,
                    dist + max(0, row.asm_end - row.asm_start),
                    tag,
                )

            i = idx_of[ekey]

            if i > 0:
                nb = arr[i - 1]
                if not nb.tag:
                    nd = dist + (e.pos - nb.pos)
                    if nd <= max_near:
                        ntag = tag
                    else:
                        ntag = f"{normalize_tag(tag, ignore_distance=False)}_distance{nd}"
                    push(nb.key(), nd, ntag)

            if i + 1 < len(arr):
                nb = arr[i + 1]
                if not nb.tag:
                    nd = dist + (nb.pos - e.pos)
                    if nd <= max_near:
                        ntag = tag
                    else:
                        ntag = f"{normalize_tag(tag, ignore_distance=False)}_distance{nd}"
                    push(nb.key(), nd, ntag)


def tags_of_record(r: AlleleRecord, asm_ends):
    lrec = asm_ends.get((r.allelename, "left"))
    rrec = asm_ends.get((r.allelename, "right"))
    l = lrec.tag if lrec is not None else ""
    rr = rrec.tag if rrec is not None else ""
    return l, rr


def reference_tag_pair_options(ref_alleles: Dict[str, RefAllele], ignore_distance=False):
    out = {}

    def add_option(options, seen, mode, pair):
        left, right = pair
        if not (left or right):
            return
        key = (left, right)
        if key not in seen:
            seen.add(key)
            options.append((mode, key))

        # Opposite-strand assembly alleles can carry reference-boundary tags on
        # opposite physical ends. Add a swapped comparison-only option.
        if left and right and left != right:
            swapped = (right, left)
            if swapped not in seen:
                seen.add(swapped)
                options.append((mode + ":swapped", swapped))

    for name, ra in ref_alleles.items():
        options = []
        seen = set()

        initiated_pair = (
            normalize_tag(ra.left_tag, ignore_distance=ignore_distance),
            normalize_tag(ra.right_tag, ignore_distance=ignore_distance),
        )
        add_option(options, seen, "initiated", initiated_pair)

        propagated_pair = (
            normalize_tag(ra.propagated_left_tag, ignore_distance=ignore_distance),
            normalize_tag(ra.propagated_right_tag, ignore_distance=ignore_distance),
        )
        add_option(options, seen, "propagated", propagated_pair)

        for index, pair in enumerate(ra.adjacent_tag_pairs):
            adjacent_pair = tuple(
                normalize_tag(tag, ignore_distance=ignore_distance)
                for tag in pair
            )
            add_option(
                options, seen, f"adjacent:{index}", adjacent_pair,
            )

        out[name] = options
    return out

def tag_pair_bydifference(asm_tags, ref_tags, diffval):
    matches = 0
    if asm_tags[0] and ref_tags[0] and asm_tags[0] == ref_tags[0]:
        matches += 1
    if asm_tags[1] and ref_tags[1] and asm_tags[1] == ref_tags[1]:
        matches += 1

    if matches == 0:
        return None

    return (0 if matches == 2 else 1, diffval)


def tag_pair_options_bydifference(asm_tags, ref_tag_options, diffval):
    """Return best score if assembly tags match any reference tag version.

    Match quality is still primarily two-end vs one-end. For equal quality and
    same diff, initiated and propagated matches are both accepted, but initiated
    is sorted first for deterministic behavior.
    """
    best = None
    for mode, ref_tags in ref_tag_options:
        base_score = tag_pair_bydifference(asm_tags, ref_tags, diffval)
        if base_score is None:
            continue

        mode_penalty = 0 if mode == "initiated" else 1
        full_score = (base_score[0], diffval, mode_penalty)
        if best is None or full_score < best:
            best = full_score

    return best


def rebuild_assembly_tags_from_matches(
    rows,
    ref_alleles,
    asm_ends,
    matched_to_ref,
    max_near=1000,
    non_reference_tags=frozenset(),
):
    new_ends = {}

    for r in rows:
        left = EndRec(
            allele=r.allelename, genome=r.genome, contig=r.asm_contig,
            pos=r.asm_start, side="left", strand=r.asm_strand,
            gene_state=r.gene_state, tag=""
        )
        right = EndRec(
            allele=r.allelename, genome=r.genome, contig=r.asm_contig,
            pos=r.asm_end, side="right", strand=r.asm_strand,
            gene_state=r.gene_state, tag=""
        )

        refname = matched_to_ref.get(r.allelename)
        if refname is not None and refname in ref_alleles:
            refa = ref_alleles[refname]
            seed_left, seed_right = seed_tags_for_reference(refname, refa)
            same_orientation = (r.asm_strand == refa.strand)
            if same_orientation:
                left.tag = seed_left
                right.tag = seed_right
            else:
                left.tag = seed_right
                right.tag = seed_left
        elif non_reference_tag_enabled(
            r.allelename, ref_alleles, non_reference_tags,
        ):
            left.tag = f"{r.allelename}{end_sign('left', r.asm_strand)}"
            right.tag = f"{r.allelename}{end_sign('right', r.asm_strand)}"

        new_ends[(r.allelename, "left")] = left
        new_ends[(r.allelename, "right")] = right

    asm_ends.clear()
    asm_ends.update(new_ends)

    propagate_assembly_tags(rows, asm_ends, max_near=max_near)


def _group_has_exact_outer_tag_match(
    query_rows: List[AlleleRecord],
    asm_ends,
    upstream_facing_tag: str,
    downstream_facing_tag: str,
):
    """Require the unmatched run to carry both facing anchor-end tags."""
    if not query_rows:
        return False
    query_left_raw = tags_of_record(query_rows[0], asm_ends)[0]
    query_right_raw = tags_of_record(query_rows[-1], asm_ends)[1]
    if not (
        query_left_raw
        and query_right_raw
        and upstream_facing_tag
        and downstream_facing_tag
    ):
        return False

    for ignore_distance in (False, True):
        query_left = normalize_tag(
            query_left_raw, ignore_distance=ignore_distance,
        )
        query_right = normalize_tag(
            query_right_raw, ignore_distance=ignore_distance,
        )
        if (
            query_left == normalize_tag(
                upstream_facing_tag, ignore_distance=ignore_distance,
            )
            and query_right == normalize_tag(
                downstream_facing_tag, ignore_distance=ignore_distance,
            )
        ):
            return True
    return False


def _merged_reference_location(
    reference_names: List[str],
    ref_alleles: Dict[str, RefAllele],
):
    refs = [ref_alleles[name] for name in reference_names]
    if not refs or len({ra.contig for ra in refs}) != 1:
        return ""
    start = min(ra.start for ra in refs)
    end = max(ra.end for ra in refs)
    strands = {ra.strand for ra in refs}
    if len(strands) != 1:
        return ""
    strand = refs[0].strand
    return format_ref_anchor_interval(refs[0].contig, start, end, strand)


def assign_grouped_anchor_fallback(
    rows: List[AlleleRecord],
    ref_alleles: Dict[str, RefAllele],
    asm_ends,
    matched_asm,
    matched_ref_keys,
    matched_to_ref,
    allowed_shape_priorities=None,
    sim_cutoff=0.1,
    assignment_stage=3,
):
    """Assign position-bracketed blocks between established fixed anchors.

    KmerMatch candidates and distances are intentionally not consulted. The
    fallback is allowed only when both outer propagated tags agree exactly.
    Shape priority is one-reference/many-query, one-query/many-reference, then
    many-query/many-reference, followed by one-query/one-reference after the
    normal similarity passes are exhausted.

    A position-only assignment is sufficient evidence that the reference is
    present, but it is not automatically a new propagation seed. An assigned
    query emits the reference's identity tags only when that exact reference
    also occurs among its KmerMatch candidates below ``sim_cutoff``.
    """
    assignable_refs = eligible_ref_alleles(ref_alleles)
    refs_by_contig = cl.defaultdict(list)
    for name, ra in assignable_refs.items():
        refs_by_contig[ra.contig].append(name)
    for names in refs_by_contig.values():
        names.sort(key=lambda name: (
            ref_alleles[name].start,
            ref_alleles[name].end,
            name,
        ))

    candidates = []
    rows_by_contig = cl.defaultdict(list)
    for row in rows:
        rows_by_contig[row.asm_contig].append(row)

    for asm_contig, contig_rows in rows_by_contig.items():
        ordered = sorted(contig_rows, key=lambda row: (
            row.asm_start, row.asm_end, row.allelename,
        ))
        anchor_indexes = [
            index for index, row in enumerate(ordered)
            if row.allelename in matched_to_ref
            and ";" not in matched_to_ref[row.allelename]
        ]
        for left_index, right_index in zip(anchor_indexes, anchor_indexes[1:]):
            query_group = [
                row for row in ordered[left_index + 1:right_index]
                if row.allelename not in matched_asm
            ]
            if not query_group:
                continue
            # Do not jump over another already assigned row.
            if len(query_group) != right_index - left_index - 1:
                continue

            left_anchor_name = matched_to_ref[ordered[left_index].allelename]
            right_anchor_name = matched_to_ref[ordered[right_index].allelename]
            left_anchor = ref_alleles.get(left_anchor_name)
            right_anchor = ref_alleles.get(right_anchor_name)
            if (
                left_anchor is None
                or right_anchor is None
                or left_anchor.contig != right_anchor.contig
            ):
                continue

            reference_forward = (
                (left_anchor.start, left_anchor.end)
                <= (right_anchor.start, right_anchor.end)
            )
            if reference_forward:
                interval_start = left_anchor.end
                interval_end = right_anchor.start
            else:
                interval_start = right_anchor.end
                interval_end = left_anchor.start
            if interval_end < interval_start:
                continue

            reference_group = []
            for name in refs_by_contig.get(left_anchor.contig, ()):
                ra = ref_alleles[name]
                if name in {left_anchor_name, right_anchor_name}:
                    continue
                if (asm_contig, name) in matched_ref_keys:
                    continue
                if ra.start >= interval_start and ra.end <= interval_end:
                    reference_group.append(name)
            if not reference_forward:
                reference_group.reverse()
            if not reference_group:
                continue

            query_count = len(query_group)
            reference_count = len(reference_group)
            if reference_count == 1 and query_count > 1:
                shape_priority = 0
            elif query_count == 1 and reference_count > 1:
                shape_priority = 1
            elif query_count > 1 and reference_count > 1:
                shape_priority = 2
            elif query_count == 1 and reference_count == 1:
                shape_priority = 3
            else:
                continue
            if (
                allowed_shape_priorities is not None
                and shape_priority not in allowed_shape_priorities
            ):
                continue
            if not _group_has_exact_outer_tag_match(
                query_group,
                asm_ends,
                tags_of_record(ordered[left_index], asm_ends)[1],
                tags_of_record(ordered[right_index], asm_ends)[0],
            ):
                continue
            candidates.append((
                shape_priority,
                asm_contig,
                query_group[0].asm_start,
                tuple(row.allelename for row in query_group),
                tuple(reference_group),
            ))

    row_by_name = {row.allelename: row for row in rows}
    assigned_groups = 0
    for _, asm_contig, _, query_names, reference_names in sorted(candidates):
        if any(name in matched_asm for name in query_names):
            continue
        if any((asm_contig, name) in matched_ref_keys for name in reference_names):
            continue
        merged_location = _merged_reference_location(
            list(reference_names), ref_alleles,
        )
        if not merged_location:
            continue
        original_locations = tuple(
            format_ref_anchor_interval(
                ref_alleles[name].contig,
                ref_alleles[name].start,
                ref_alleles[name].end,
                ref_alleles[name].strand,
            )
            for name in reference_names
        )
        joined_names = ";".join(reference_names)
        part_count = len(query_names)
        query_rows = tuple(row_by_name[name] for name in query_names)
        query_original_locations = tuple(row.assembly_loc for row in query_rows)
        query_left_extensions = tuple(
            str(row.allele_left_offset) for row in query_rows
        )
        query_right_extensions = tuple(
            str(row.allele_right_offset) for row in query_rows
        )
        for part_index, query_name in enumerate(query_names, 1):
            row = row_by_name[query_name]
            matched_asm.add(query_name)
            # Assignment and tag emission are deliberately separate. The
            # positional bracket is enough to map this row and suppress a
            # false DEL, whereas only corroborating sequence similarity may
            # promote a single-reference assignment into a fixed tag source.
            if len(reference_names) == 1 and any(
                name == reference_names[0] and difference < sim_cutoff
                for name, difference in row.parsed_similar_refs
            ):
                matched_to_ref[query_name] = reference_names[0]
            row.matched_ref_allele = joined_names
            row.assigned_ref_difference = math.inf
            row.assigned_ref_alleles_bylocation = joined_names
            row.assigned_ref_alleles_coordinates = merged_location
            row.matched_ref_coords = merged_location
            row.grouped_ref_names = tuple(reference_names)
            row.grouped_ref_original_locs = original_locations
            row.grouped_part_index = part_index if part_count > 1 else 0
            row.grouped_part_count = part_count
            row.grouped_query_names = tuple(query_names)
            row.grouped_query_original_locs = query_original_locations
            row.grouped_query_left_extensions = query_left_extensions
            row.grouped_query_right_extensions = query_right_extensions
            row.assignment_stage = assignment_stage
        for reference_name in reference_names:
            matched_ref_keys.add((asm_contig, reference_name))
        assigned_groups += 1
    return assigned_groups


def assign_matches(
    rows,
    ref_alleles,
    asm_ends,
    anchors,
    max_near=1000,
    sim_cutoff=0.1,
    non_reference_tags=frozenset(),
):
    all_ref_alleles = ref_alleles
    if ref_alleles is GLOBAL_REF_ALLELES and GLOBAL_ASSIGNABLE_REF_ALLELES is not None:
        assignable_ref_alleles = GLOBAL_ASSIGNABLE_REF_ALLELES
    else:
        assignable_ref_alleles = eligible_ref_alleles(ref_alleles)

    matched_ref_keys = set()
    matched_asm = set()
    matched_to_ref = {}

    def commit_match(r, refname, diffval, assignment_stage):
        refkey = (r.asm_contig, refname)
        matched_asm.add(r.allelename)
        matched_ref_keys.add(refkey)
        matched_to_ref[r.allelename] = refname

        r.matched_ref_allele = refname
        r.assigned_ref_difference = diffval

        ra = assignable_ref_alleles[refname]
        r.assigned_ref_alleles_bylocation = location_refname_for_output(refname, all_ref_alleles)
        r.assigned_ref_alleles_coordinates = location_coords_for_output(refname, all_ref_alleles)
        r.matched_ref_coords = original_ref_coords(ra)
        r.assignment_stage = assignment_stage

    for r in rows:
        if r.allelename in anchors:
            refname = anchors[r.allelename]
            if refname not in assignable_ref_alleles:
                continue

            diffval = math.inf
            for nm, dv in r.parsed_similar_refs:
                if nm == refname:
                    diffval = dv
                    break

            commit_match(r, refname, diffval, 1)

    rebuild_assembly_tags_from_matches(
        rows,
        all_ref_alleles,
        asm_ends,
        matched_to_ref,
        max_near=max_near,
        non_reference_tags=non_reference_tags,
    )

    row_map = {r.allelename: r for r in rows}

    for ignore_distance in (False, True):
        if all_ref_alleles is GLOBAL_REF_ALLELES and GLOBAL_REF_TAG_OPTIONS is not None:
            ref_tags = GLOBAL_REF_TAG_OPTIONS[ignore_distance]
        else:
            ref_tags = reference_tag_pair_options(
                all_ref_alleles,
                ignore_distance=ignore_distance,
            )

        progress = True
        while progress:
            progress = False

            # A reference block that was split into multiple consecutive query
            # blocks must be claimed as one group before a KmerMatch-bearing
            # member can consume it alone. Both outer fixed anchors are
            # required, so this remains position-constrained even when some
            # group members are Novel and have no reference candidate.
            if assign_grouped_anchor_fallback(
                rows,
                all_ref_alleles,
                asm_ends,
                matched_asm,
                matched_ref_keys,
                matched_to_ref,
                allowed_shape_priorities={0},
                sim_cutoff=sim_cutoff,
                assignment_stage=2 if not ignore_distance else 4,
            ):
                rebuild_assembly_tags_from_matches(
                    rows,
                    all_ref_alleles,
                    asm_ends,
                    matched_to_ref,
                    max_near=max_near,
                    non_reference_tags=non_reference_tags,
                )
                progress = True
                continue

            exact_new = []
            for r in rows:
                if r.allelename in matched_asm:
                    continue

                asm_tags = tuple(
                    normalize_tag(x, ignore_distance=ignore_distance)
                    for x in tags_of_record(r, asm_ends)
                )

                candidates = []
                for refname, diffval in r.parsed_similar_refs:
                    if diffval >= sim_cutoff:
                        continue
                    refkey = (r.asm_contig, refname)
                    if refname not in all_ref_alleles or refkey in matched_ref_keys:
                        continue
                    score = tag_pair_options_bydifference(asm_tags, ref_tags[refname], diffval)
                    if score is not None and score[0] == 0:
                        candidates.append((score, refname, diffval))

                if len(candidates) == 1:
                    _, refname, diffval = candidates[0]
                    if refname in assignable_ref_alleles:
                        exact_new.append((r.allelename, refname, diffval))

            if exact_new:
                for allelename, refname, diffval in exact_new:
                    if allelename in matched_asm:
                        continue
                    refkey = (row_map[allelename].asm_contig, refname)
                    if refkey in matched_ref_keys:
                        continue
                    commit_match(
                        row_map[allelename], refname, diffval,
                        2 if not ignore_distance else 4,
                    )

                rebuild_assembly_tags_from_matches(
                    rows,
                    all_ref_alleles,
                    asm_ends,
                    matched_to_ref,
                    max_near=max_near,
                    non_reference_tags=non_reference_tags,
                )
                progress = True
                continue

            # Exact two-end one-to-one matching is exhausted. Recover a
            # split/merged run before the older one-end greedy pass can consume
            # one member of a legitimate one-to-many group.
            if assign_grouped_anchor_fallback(
                rows,
                all_ref_alleles,
                asm_ends,
                matched_asm,
                matched_ref_keys,
                matched_to_ref,
                allowed_shape_priorities={0, 1, 2},
                sim_cutoff=sim_cutoff,
                assignment_stage=2 if not ignore_distance else 4,
            ):
                rebuild_assembly_tags_from_matches(
                    rows,
                    all_ref_alleles,
                    asm_ends,
                    matched_to_ref,
                    max_near=max_near,
                    non_reference_tags=non_reference_tags,
                )
                progress = True
                continue

            row_best = []

            for r in rows:
                if r.allelename in matched_asm:
                    continue

                asm_tags = tuple(
                    normalize_tag(x, ignore_distance=ignore_distance)
                    for x in tags_of_record(r, asm_ends)
                )

                row_candidates = []
                for refname, diffval in r.parsed_similar_refs:
                    if diffval >= sim_cutoff:
                        continue
                    refkey = (r.asm_contig, refname)
                    if refname not in all_ref_alleles or refkey in matched_ref_keys:
                        continue

                    score = tag_pair_options_bydifference(asm_tags, ref_tags[refname], diffval)
                    if score is not None:
                        row_candidates.append((score, refname, diffval))

                if not row_candidates:
                    continue

                row_candidates.sort()
                best_score, best_refname, best_diffval = row_candidates[0]

                if best_refname not in assignable_ref_alleles:
                    continue

                row_best.append((best_score, r.allelename, r.asm_contig, best_refname, best_diffval))

            row_best.sort()

            accepted = False
            for score, allelename, asm_contig, refname, diffval in row_best:
                refkey = (asm_contig, refname)
                if allelename in matched_asm or refkey in matched_ref_keys:
                    continue

                commit_match(
                    row_map[allelename], refname, diffval,
                    2 if not ignore_distance else 4,
                )
                rebuild_assembly_tags_from_matches(
                    rows,
                    all_ref_alleles,
                    asm_ends,
                    matched_to_ref,
                    max_near=max_near,
                    non_reference_tags=non_reference_tags,
                )
                progress = True
                accepted = True
                break

            if not accepted:
                break

    if assign_grouped_anchor_fallback(
        rows,
        all_ref_alleles,
        asm_ends,
        matched_asm,
        matched_ref_keys,
        matched_to_ref,
        sim_cutoff=sim_cutoff,
        assignment_stage=3,
    ):
        rebuild_assembly_tags_from_matches(
            rows,
            all_ref_alleles,
            asm_ends,
            matched_to_ref,
            max_near=max_near,
            non_reference_tags=non_reference_tags,
        )


def write_genome_output(rows, asm_ends, ref_alleles, outpath, write_header=False, emit_nonunique_placeholder=False):
    header = [
        "allelename",
        "best_ref_allele",
        "assembly_loc",
        "ref_loc",
        "IntronorExon",
        "similarity",
        "classified_type",
        "all_similar_ref_alleles",
        "assigned_ref_alleles_bylocation",
        "assigned_ref_alleles_coordinates",
        "asm_left_tag",
        "asm_right_tag",
        "allele_left_offset",
        "allele_right_offset",
        *FINAL_STAGE_COLUMNS,
    ]

    if ref_alleles is GLOBAL_REF_ALLELES and GLOBAL_REF_ANCHOR_INDEX is not None:
        anchor_by_seed, ordered_seeds = GLOBAL_REF_ANCHOR_INDEX
    else:
        anchor_by_seed, ordered_seeds = build_reference_anchor_index(ref_alleles)

    with open(outpath, "w") as w:
        if write_header:
            w.write("\t".join(header) + "\n")

        for r in sorted(rows, key=allele_record_assembly_sort_key):
            if emit_nonunique_placeholder and r.allele_unique_len <= 0:
                row = [
                    r.allelename,
                    r.best_ref_allele,
                    r.assembly_loc,
                    r.ref_loc,
                    r.intron_or_exon,
                    r.similarity_raw,
                    r.class_type,
                    "",
                    "",
                    "",
                    "",
                    "",
                    NEG_INF_STR,
                    NEG_INF_STR,
                    # Columns 15-18; finalize_lift_stages replaces them.
                    "0",
                    "",
                    ".",
                    ".",
                ]
                w.write("\t".join(row) + "\n")
                continue

            ltag, rtag = tags_of_record(r, asm_ends)

            best_ref_display = r.best_ref_allele
            if best_ref_display and r.gene_state and r.gene_state != r.best_ref_allele:
                best_ref_display = f"{r.best_ref_allele}[{r.gene_state}]"

            ref_loc_display = r.ref_loc
            assembly_loc_display = r.assembly_loc
            similarity_display = r.similarity_raw
            class_display = r.class_type
            similar_refs_display = r.similar_refs_raw
            left_extension_display = str(r.allele_left_offset)
            right_extension_display = str(r.allele_right_offset)
            if r.grouped_ref_names:
                joined_names = ";".join(r.grouped_ref_names)
                best_ref_display = joined_names
                ref_loc_display = ";".join(r.grouped_ref_original_locs)
                similarity_display = "."
                similar_refs_display = joined_names
                if r.grouped_part_index:
                    class_display = f"part_{r.grouped_part_index}"
                    # Repeat the complete, index-aligned query metadata on
                    # every part row.  The row's part_N value selects its own
                    # entry, while downstream grouped comparisons retain all
                    # original block boundaries.
                    assembly_loc_display = ";".join(
                        r.grouped_query_original_locs
                    )
                    left_extension_display = ";".join(
                        r.grouped_query_left_extensions
                    )
                    right_extension_display = ";".join(
                        r.grouped_query_right_extensions
                    )

            assigned_ref_coords = r.assigned_ref_alleles_coordinates
            if (not r.assigned_ref_alleles_bylocation) and (not assigned_ref_coords):
                assigned_ref_coords = infer_ref_coords_from_asm_tags(ltag, rtag, anchor_by_seed, ordered_seeds)

            row = [
                r.allelename,
                best_ref_display,
                assembly_loc_display,
                ref_loc_display,
                r.intron_or_exon,
                similarity_display,
                class_display,
                similar_refs_display,
                r.assigned_ref_alleles_bylocation,
                assigned_ref_coords,
                ltag,
                rtag,
                left_extension_display,
                right_extension_display,
                # Columns 15-18; finalize_lift_stages replaces them
                # (17/18: location and sequence mappings).
                str(r.assignment_stage),
                "",
                ".",
                ".",
            ]
            w.write("\t".join(row) + "\n")


def process_loaded_rows(
    all_rows,
    ref_alleles,
    sim_cutoff,
    near_dist,
    outpath,
    non_reference_tags=None,
    assembly_only=False,
):
    """Process one genome or one contig and write its sorted output rows."""
    if non_reference_tags is None:
        non_reference_tags = GLOBAL_NON_REFERENCE_TAGS
    assign_assembly_unique_regions(all_rows)

    if assembly_only:
        # Reference-free mode deliberately stops after assembly-coordinate
        # ownership has been resolved.  Empty end/reference maps leave the
        # reference-assignment and tag columns blank while preserving the
        # assembly gap/overlap results in columns 13 and 14.
        write_genome_output(
            all_rows,
            {},
            {},
            outpath,
            write_header=False,
            emit_nonunique_placeholder=True,
        )
        return outpath

    valid_rows = [r for r in all_rows if r.allele_unique_len > 0]

    if not valid_rows:
        write_genome_output(
            all_rows,
            {},
            ref_alleles,
            outpath,
            write_header=False,
            emit_nonunique_placeholder=True,
        )
        return outpath

    anchors = confident_one_to_one(valid_rows, ref_alleles, sim_cutoff=sim_cutoff)
    asm_ends = initialize_assembly_tags(
        valid_rows,
        ref_alleles,
        anchors,
        non_reference_tags=non_reference_tags,
    )
    propagate_assembly_tags(valid_rows, asm_ends, max_near=near_dist)
    assign_matches(
        valid_rows,
        ref_alleles,
        asm_ends,
        anchors,
        max_near=near_dist,
        sim_cutoff=sim_cutoff,
        non_reference_tags=non_reference_tags,
    )

    write_genome_output(
        all_rows,
        asm_ends,
        ref_alleles,
        outpath,
        write_header=False,
        emit_nonunique_placeholder=True,
    )
    return outpath


def process_one_genome(task):
    (
        genome,
        genome_file,
        sim_cutoff,
        near_dist,
        outdir,
        filter_hg38_alts,
        refhaplo,
        assembly_only,
    ) = task
    ref_alleles = GLOBAL_REF_ALLELES

    print("Processing", genome, flush=True)
    all_rows = load_genome_records(
        genome_file,
        ref_alleles,
        filter_hg38_alts=filter_hg38_alts,
    )
    outpath = os.path.join(outdir, f"{genome}.out.tsv")
    return process_loaded_rows(
        all_rows,
        ref_alleles,
        sim_cutoff,
        near_dist,
        outpath,
        assembly_only=assembly_only,
    )


def process_one_contig(task):
    """Worker entry point for contig-level multiprocessing."""
    task_index, rows, sim_cutoff, near_dist, part_path, assembly_only = task
    process_loaded_rows(
        rows,
        GLOBAL_REF_ALLELES,
        sim_cutoff,
        near_dist,
        part_path,
        assembly_only=assembly_only,
    )
    return task_index, part_path


def output_line_assembly_sort_key(line):
    """Recreate allele_record_assembly_sort_key from one output line."""
    parts = line.rstrip("\n").split("\t")
    allelename = parts[0] if parts else ""
    assembly_loc = parts[2] if len(parts) > 2 else ""
    try:
        contig, start, end, _ = parse_loc(assembly_loc)
        return (contig_sort_key(contig), start, end, allelename)
    except Exception:
        return ((2, assembly_loc), 0, 0, allelename)


def _merge_sorted_output_files(part_paths, outpath):
    """Merge one descriptor-bounded batch of sorted output files."""
    import heapq

    handles = [open(path, "r") for path in part_paths]
    try:
        with open(outpath, "w") as output:
            for line in heapq.merge(
                *handles,
                key=output_line_assembly_sort_key,
            ):
                output.write(line)
    finally:
        for handle in handles:
            handle.close()


def merge_contig_outputs(
    part_paths, outpath, fan_in=CONTIG_MERGE_FAN_IN,
):
    """Merge sorted contig outputs with a bounded number of open files."""
    part_paths = list(part_paths)
    fan_in = int(fan_in)
    if fan_in < 2:
        raise ValueError("contig merge fan-in must be at least 2")

    output_dir = os.path.dirname(os.path.abspath(outpath)) or "."
    merge_dir = tempfile.mkdtemp(
        prefix=f".{os.path.basename(outpath)}.merge.",
        dir=output_dir,
    )
    try:
        current_paths = part_paths
        round_index = 0
        while len(current_paths) > fan_in:
            next_paths = []
            for batch_index, start in enumerate(
                range(0, len(current_paths), fan_in)
            ):
                merged_path = os.path.join(
                    merge_dir,
                    f"round_{round_index:04d}_{batch_index:06d}.tsv",
                )
                _merge_sorted_output_files(
                    current_paths[start:start + fan_in], merged_path,
                )
                next_paths.append(merged_path)
            current_paths = next_paths
            round_index += 1

        final_path = os.path.join(merge_dir, "final.tsv")
        _merge_sorted_output_files(current_paths, final_path)
        os.replace(final_path, outpath)
    finally:
        shutil.rmtree(merge_dir, ignore_errors=True)


def process_one_genome_by_contig(
    task,
    ref_alleles,
    priority_by_prefix,
    threads,
    worker_pool=None,
    non_reference_tags=None,
    alignment_scores=None,
):
    """Use multiple processes across independent contigs of one genome."""
    (
        genome,
        genome_file,
        sim_cutoff,
        near_dist,
        outdir,
        filter_hg38_alts,
        refhaplo,
        assembly_only,
    ) = task
    print("Processing", genome, "by contig", flush=True)

    all_rows = load_genome_records(
        genome_file,
        ref_alleles,
        filter_hg38_alts=filter_hg38_alts,
    )
    rows_by_contig = cl.defaultdict(list)
    for row in all_rows:
        rows_by_contig[row.asm_contig].append(row)

    outpath = os.path.join(outdir, f"{genome}.out.tsv")
    if not rows_by_contig:
        with open(outpath, "w"):
            pass
        return outpath

    ordered_contigs = sorted(
        rows_by_contig,
        key=lambda contig: (contig_sort_key(contig), contig),
    )
    parts_dir = os.path.join(outdir, f"{genome}.parts")
    os.makedirs(parts_dir, exist_ok=True)
    contig_tasks = [
        (
            index,
            rows_by_contig[contig],
            sim_cutoff,
            near_dist,
            os.path.join(parts_dir, f"{index:06d}.tsv"),
            assembly_only,
        )
        for index, contig in enumerate(ordered_contigs)
    ]

    workers = min(threads, len(contig_tasks))
    if non_reference_tags is None:
        non_reference_tags = GLOBAL_NON_REFERENCE_TAGS
    if worker_pool is not None:
        completed = list(worker_pool.imap_unordered(
            process_one_contig,
            contig_tasks,
            chunksize=1,
        ))
    elif workers == 1:
        init_worker(
            ref_alleles,
            priority_by_prefix,
            non_reference_tags,
            alignment_scores,
        )
        completed = [process_one_contig(item) for item in contig_tasks]
    else:
        with Pool(
            processes=workers,
            initializer=init_worker,
            initargs=(
                ref_alleles,
                priority_by_prefix,
                non_reference_tags,
                alignment_scores,
            ),
        ) as pool:
            completed = list(pool.imap_unordered(
                process_one_contig,
                contig_tasks,
                chunksize=1,
            ))

    completed.sort()
    merge_contig_outputs([path for _, path in completed], outpath)
    return outpath


def count_genome_contigs(genome_file, filter_hg38_alts=True):
    """Count usable assembly contigs for multiprocessing strategy selection."""
    contigs = set()
    for parts in iter_rows(genome_file):
        try:
            allelename = parts[0]
            contig, _, _, _ = parse_loc(parts[2])
            if filter_hg38_alts and is_hg38_alt_assembly(allelename, contig):
                continue
            contigs.add(contig)
        except Exception:
            continue
    return len(contigs)


def main():
    parser = argparse.ArgumentParser(
        description="Split-first streaming/multiprocess allele annotation."
    )
    parser.add_argument("-i", "--input", required=True, help="Input TSV (.gz supported)")
    parser.add_argument("-o", "--output", required=True, help="Output TSV")
    parser.add_argument("--refhaplo", default="CHM13_h1")
    parser.add_argument(
        "--no-reference",
        action="store_true",
        help=(
            "run without a reference haplotype; compute only assembly-contig "
            "gap/overlap ownership (output columns 13 and 14)"
        ),
    )
    parser.add_argument(
        "--sim-cutoff",
        type=float,
        default=0.1,
        help="Strict diff cutoff for confident unique mapping anchors (default: 0.1 = similarity > 0.9)",
    )
    parser.add_argument(
        "--gene-state-cutoff",
        type=float,
        default=0.1,
        help="Maximum reported difference allowed to connect reference alleles into the same gene_state (default: 0.1)",
    )
    parser.add_argument(
        "--non-reference-tags",
        default="",
        help=(
            "comma-separated exact non-reference allele names allowed to "
            "seed their own tags, or 'all' to enable every non-reference "
            "allele (default: none)"
        ),
    )
    parser.add_argument(
        "--alignments",
        default="",
        help=(
            "QUERY_align.txt[,REFERENCE_align.txt]: the sample _align.txt "
            "ranks overlapping PAs and, with the reference haplotype "
            "_align.txt, gives the query/reference coverage used by the "
            "final lift stages. Without the second path the reference file "
            "is taken from REFHAPLO/REFHAPLO.align.txt beside the sample "
            "folder"
        ),
    )
    parser.add_argument("--near-dist", type=int, default=1000)
    parser.add_argument(
        "--max-extension",
        type=int,
        default=10000,
        help=(
            "maximum GenomeLift column-13/14 context included in appended "
            "DEL/NA unmatched-reference intervals (default: 10000)"
        ),
    )
    parser.add_argument("-t","--threads", type=int, default=8)
    parser.add_argument(
        "--parallel-unit",
        choices=("auto", "genome", "contig"),
        default="auto",
        help=(
            "multiprocessing unit (default: auto; uses contigs when too few "
            "genomes are available to occupy the requested workers)"
        ),
    )
    parser.add_argument("--tmpdir", default="", help="Temporary directory or existing debug dir")
    parser.add_argument("--reuse-split", action="store_true", help="Reuse existing split files in tmpdir/genomes")
    parser.add_argument(
        "--filter-hg38-alts",
        dest="filter_hg38_alts",
        action="store_true",
        help="Filter HG38_h1 alternative/non-canonical contigs (default: on)",
    )
    parser.add_argument(
        "--keep-hg38-alts",
        dest="filter_hg38_alts",
        action="store_false",
        default=True,
        help="Keep HG38_h1 alternative/non-canonical contigs",
    )
    args = parser.parse_args()
    if args.threads < 1:
        parser.error("--threads must be positive")
    if args.max_extension < 0:
        parser.error("--max-extension must be non-negative")
    try:
        non_reference_tags = parse_non_reference_tags(
            args.non_reference_tags,
        )
    except ValueError as exc:
        parser.error(str(exc))
    alignment_paths = [
        value.strip() for value in args.alignments.split(",")
        if value.strip()
    ]
    if len(alignment_paths) > 2:
        parser.error(
            "--alignments takes QUERY_align.txt[,REFERENCE_align.txt]"
        )
    for alignment_path in alignment_paths:
        if not os.path.isfile(alignment_path):
            parser.error(f"alignment file not found: {alignment_path}")
    args.alignments = alignment_paths[0] if alignment_paths else ""
    args.reference_alignments = (
        alignment_paths[1] if len(alignment_paths) > 1 else ""
    )
    if args.alignments and not args.no_reference:
        if not args.reference_alignments:
            # SAMPLES/SAMPLE/SAMPLE.align.txt -> SAMPLES/REFHAPLO/REFHAPLO.align.txt
            args.reference_alignments = os.path.join(
                os.path.dirname(os.path.dirname(
                    os.path.abspath(args.alignments)
                )),
                args.refhaplo, f"{args.refhaplo}.align.txt",
            )
            if not os.path.isfile(args.reference_alignments):
                parser.error(
                    "final lift stages need the reference haplotype "
                    "_align.txt for coverage; pass --alignments "
                    "SAMPLE_align.txt,REFERENCE_align.txt (not found: "
                    f"{args.reference_alignments})"
                )
        print(
            "Final-stage coverage alignments: "
            f"{args.alignments},{args.reference_alignments}",
            flush=True,
        )

    alignment_scores = None
    if args.alignments:
        try:
            alignment_scores = read_alignment_scores(
                args.alignments, read_pa_loci(args.input),
            )
        except (OSError, ValueError) as exc:
            parser.error(str(exc))
        print(
            f"Loaded locus-matched alignment scores for "
            f"{len(alignment_scores)} PAs",
            flush=True,
        )

    tmpdir = args.tmpdir if args.tmpdir else "GenomeLiftTemp"
    remove_tmpdir_after_finish = (not args.tmpdir)

    try:
        genome_split_dir = os.path.join(tmpdir, "genomes")
        genome_out_dir = os.path.join(tmpdir, "out")
        os.makedirs(genome_out_dir, exist_ok=True)

        if args.reuse_split:
            genome_files = {
                fn[:-4]: os.path.join(genome_split_dir, fn)
                for fn in os.listdir(genome_split_dir)
                if fn.endswith(".tsv")
            }
        else:
            genome_files = split_input_by_genome(
                args.input,
                genome_split_dir,
                filter_hg38_alts=args.filter_hg38_alts,
            )

        if not genome_files:
            raise ValueError(f"No genome split files found in {genome_split_dir}")

        if args.no_reference:
            print(
                "No-reference mode: computing assembly-contig gaps/overlaps only",
                flush=True,
            )
            ref_alleles = {}
            priority_by_prefix = {}
        else:
            ref_file = genome_files.get(args.refhaplo)
            if ref_file is None:
                raise ValueError(
                    f"Reference genome file not found: {args.refhaplo}"
                )

            print("Analyzing reference alleles", flush=True)
            ref_alleles, ref_edges = collect_reference_from_refgenome_file(
                ref_file,
                args.refhaplo,
                args.gene_state_cutoff,
            )

            # Priority must be shared across every sample.  Capture it from
            # the original reference catalog before promoted alleles are
            # added; sample assembly lengths and RefMatch distances are
            # deliberately excluded.
            print("Building reference-size priority index", flush=True)
            priority_by_prefix = build_reference_size_priority_index(
                ref_alleles,
            )

            print("Promoting novel reference alleles", flush=True)
            promoted = supplement_reference_alleles_from_all_genomes(
                genome_files,
                ref_alleles,
                args.refhaplo,
            )
            print(f"Added {promoted} promoted reference alleles", flush=True)

            print("Finding unique regions", flush=True)
            assign_reference_unique_regions(
                ref_alleles,
                priority_by_prefix=priority_by_prefix,
            )

            print("Building gene states", flush=True)
            build_gene_states_from_edges(ref_alleles, ref_edges)

            print("Annotating reference ends", flush=True)
            ref_ends = initialize_reference_end_tags(ref_alleles)

            print("Propagating reference tags", flush=True)
            propagate_reference_tags(
                ref_alleles,
                ref_ends,
                max_near=args.near_dist,
            )
            finalize_delegate_reference_tags(ref_alleles)

            print("Building reference propagated/control tags", flush=True)
            assign_reference_propagated_tags(
                ref_alleles,
                ref_ends,
                max_near=args.near_dist,
            )

            print("Finished annotating reference", flush=True)

        tasks = [
            (
                genome,
                gfile,
                args.sim_cutoff,
                args.near_dist,
                genome_out_dir,
                args.filter_hg38_alts,
                args.refhaplo,
                args.no_reference,
            )
            for genome, gfile in genome_files.items()
        ]

        parallel_unit = args.parallel_unit
        contig_counts = None
        if parallel_unit == "auto":
            contig_counts = [
                count_genome_contigs(
                    task[1],
                    filter_hg38_alts=args.filter_hg38_alts,
                )
                for task in tasks
            ]
            if (
                args.threads > len(tasks)
                and max(contig_counts, default=0) > 1
            ):
                parallel_unit = "contig"
            else:
                parallel_unit = "genome"

        print(f"Multiprocessing by {parallel_unit}", flush=True)

        if args.threads == 1:
            init_worker(
                ref_alleles,
                priority_by_prefix,
                non_reference_tags,
                alignment_scores,
            )
            outfiles = [process_one_genome(t) for t in tasks]
        elif parallel_unit == "contig":
            # Reuse one pool for every genome so the annotated reference index
            # is initialized only once in each worker.
            if contig_counts is None:
                contig_counts = [
                    count_genome_contigs(
                        task[1],
                        filter_hg38_alts=args.filter_hg38_alts,
                    )
                    for task in tasks
                ]
            contig_workers = min(
                args.threads,
                max(contig_counts, default=1),
            )
            with Pool(
                processes=contig_workers,
                initializer=init_worker,
                initargs=(
                    ref_alleles,
                    priority_by_prefix,
                    non_reference_tags,
                    alignment_scores,
                ),
            ) as pool:
                outfiles = [
                    process_one_genome_by_contig(
                        task,
                        ref_alleles,
                        priority_by_prefix,
                        args.threads,
                        worker_pool=pool,
                        non_reference_tags=non_reference_tags,
                        alignment_scores=alignment_scores,
                    )
                    for task in tasks
                ]
        else:
            with Pool(
                processes=min(args.threads, len(tasks)),
                initializer=init_worker,
                initargs=(
                    ref_alleles,
                    priority_by_prefix,
                    non_reference_tags,
                    alignment_scores,
                ),
            ) as pool:
                outfiles = pool.map(process_one_genome, tasks)

        with open(args.output, "w") as w:
            for outfile in sorted(outfiles):
                with open(outfile, "r") as f:
                    for line in f:
                        w.write(line)

        if not args.no_reference:
            stage_counts = finalize_lift_stages(
                args.output, ref_alleles, args.alignments,
                args.reference_alignments,
            )
            print(
                "Final lift stages: "
                + ", ".join(
                    f"stage {stage}={stage_counts.get(stage, 0)}"
                    for stage in range(5)
                ),
                flush=True,
            )
            deletion_count, unknown_count = append_final_unmatched_rows(
                args.output,
                args.refhaplo,
                args.max_extension,
            )
            print(
                f"Appended {deletion_count} bracketed DEL row(s) and "
                f"{unknown_count} unanchored NA row(s) for unmatched "
                "reference runs",
                flush=True,
            )

    finally:
        if remove_tmpdir_after_finish and os.path.isdir(tmpdir):
            # shutil.rmtree(tmpdir)
            pass


if __name__ == "__main__":
    main()
