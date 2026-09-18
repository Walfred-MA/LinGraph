#!/usr/bin/env python3
import shutil
import argparse
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


def read_alignment_scores(path: str) -> Dict[str, int]:
    """Read best complete-graph alignment score for every PA name."""
    if not path:
        return {}
    scores: Dict[str, int] = {}
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
            allele = fields[5].strip()
            graph_cigar = fields[7].strip()
            if not allele or not graph_cigar or graph_cigar == "*":
                continue
            score = graph_alignment_affine_score(graph_cigar)
            for key in (allele, canonical_alignment_allele_name(allele)):
                scores[key] = max(score, scores.get(key, score))
    return scores


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
            rows.append(FinalLiftRow(
                fields=tuple(fields),
                allele=fields[0],
                genome=genome,
                contig=contig,
                start=start,
                end=end,
                strand=strand,
                assigned_refs=assigned_refs,
                class_type=fields[6].strip(),
                part_index=part_index,
                left_extension=part_value(fields[12]),
                right_extension=part_value(fields[13]),
            ))
    return rows


def _choose_final_deletion_anchors(left_rows, right_rows):
    candidates = []
    for left in left_rows:
        for right in right_rows:
            if left.contig != right.contig:
                continue
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
    """
    if row.assigned_refs:
        return row.assigned_refs
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
            if row.part_index == 0 and len(row.assigned_refs) == 1:
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
                    chosen = _choose_final_deletion_anchors(
                        anchors.get(left_ref.allele, ()),
                        anchors.get(right_ref.allele, ()),
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
        return re.sub(r"_distance\d+$", "", tag)
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

    def commit_match(r, refname, diffval):
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

            commit_match(r, refname, diffval)

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
                    commit_match(row_map[allelename], refname, diffval)

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

                commit_match(row_map[allelename], refname, diffval)
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
            "optional sample _align.txt used only to rank overlapping PAs "
            "from the same graph; omitted preserves legacy ownership"
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
    if args.alignments and not os.path.isfile(args.alignments):
        parser.error(f"alignment file not found: {args.alignments}")

    alignment_scores = None
    if args.alignments:
        try:
            alignment_scores = read_alignment_scores(args.alignments)
        except (OSError, ValueError) as exc:
            parser.error(str(exc))
        print(
            f"Loaded alignment scores for {len(alignment_scores)} PA "
            "name/alias entries",
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
