#!/usr/bin/env python3

"""

QuickTriocheck.py

Check trio consistency in a multi-sample VCF or separate haplotype VCFs.

Multiple child VCFs produce combined and per-haplotype summaries, in input order.

This version supports two sample encodings:

1. Standard VCF GT fields, for example:

     FORMAT=GT

     sample=0/1, 1/1, 1|0

2. HSV fields like the uploaded graphcigartoref_clean.smoke.sv.vcf, for example:

     FORMAT=HSV

     sample=1:DEL:-376:<3H3D27=346D:...

     sample=0|1:DEL:-376:...

     sample=0, 0|0, or . for non-carriers

Default matching:

- exact sample matching by default

- use --prefix-match to allow HG00514 to match HG00514_h1 and HG00514_h2

- child/mother/father selected columns must not overlap

- SV matching is fuzzy by merged-VCF breakpoint and size ratio

- --original-distance switches SV matching to reconstructed per-haplotype original REFERENCE coordinates

  recovered from row POS plus the artificial leading H padding in each HSV EXTENDGRAPHCIGAR

  * --original-distance 0: original reference intervals must overlap

  * --original-distance N: original reference intervals may be separated by at most N bp

  * --original-distance exact: original reference start and end coordinates must be identical

- SNP/INDEL matching is exact

- child VarINS / HSV INS_ alleles are excluded by default; use --VarINS to include them

- in separate-haplotype mode, parental missingness is determined from ##referenceCoverage header intervals

- child-carried alleles are excluded if either parent lacks reference coverage at the child reference anchor

- separate-haplotype mode trims child calls within 10000 bp of either end of

  the child assembly query contig *and* the corresponding raw referenceCoverage

  intervals by default. Query-contig lengths come from the FASTA indexes in

  ``--query-paths``; reference-side boundaries come from ``##referenceCoverage``.

  ``--trim-edges [BP]`` overrides both trims; ``--trim-edges 0`` disables them.

- child-only PA-edge QC reconstructs every PA from QUERYCOORD + LABEL_H and

  excludes a child record when any PA is <1000 bp from either contig edge by

  default. ``--pa-edge-distance 0`` disables it. Parents and parental

  referenceCoverage is affected only by ``--trim-edges``; PA-edge filtering remains

  child-query-side only.

"""

from __future__ import annotations

import argparse

import gzip

import os

import re

import shlex

import sys

from bisect import bisect_right

from collections import Counter, defaultdict

from dataclasses import dataclass, replace

from typing import Dict, List, Optional, Sequence, Tuple, Union


DEFAULT_BREAKPOINT_SLOP = 500

DEFAULT_MIN_SIZE = 50

SIZE_RATIO_MIN = 0.7

DEFAULT_PA_EDGE_DISTANCE = 1000

SCRIPT_VERSION = "decomposed-format-v6.6-reference-edge-and-mixed-sv-components-2026-09-26"


@dataclass(frozen=True)

class OriginalLocus:

    sample: str

    chrom: str

    start: int

    end: int


@dataclass(frozen=True)

class QueryLocus:

    contig: str

    start: int

    end: int

    strand: str


@dataclass(frozen=True)

class VariantRecord:

    chrom: str

    pos: int

    end: int

    svtype: str

    size: int

    vid: str

    ref: str

    alt: str

    allele_index: int

    variant_class: str

    source: str

    line_no: int

    original_loci: Tuple[OriginalLocus, ...] = ()

    query_loci: Tuple[QueryLocus, ...] = ()

    pa_loci: Tuple[QueryLocus, ...] = ()


@dataclass(frozen=True)

class SVComponent:

    svtype: str

    size: int

    original_loci: Tuple[OriginalLocus, ...] = ()

    query_loci: Tuple[QueryLocus, ...] = ()

    pa_loci: Tuple[QueryLocus, ...] = ()


INFO_RE = re.compile(r"([^=;]+)(?:=([^;]*))?")

SVTYPE_FROM_ALT = re.compile(r"^<([^>]+)>$")

INT_RE = re.compile(r"^-?\d+$")

GT_SEPS = re.compile(r"[|/]")

REFERENCE_COVERAGE_RE = re.compile(

    r'^##referenceCoverage=<Sample="([^"]+)",Chrom="([^"]+)",Start=(\d+),End=(\d+)>$'

)



def merge_reference_coverage_intervals(

    intervals_by_chrom: Dict[str, List[Tuple[int, int]]],

    edge_trim: int = 0,

) -> Dict[str, List[Tuple[int, int]]]:

    """Trim, then merge 0-based half-open referenceCoverage intervals.

    ``edge_trim`` bp are removed independently from both ends of every raw

    coverage interval *before* merging.  Trimming before merging is important:

    it preserves an unbenchmarkable flank at each alignment/coverage boundary

    instead of accidentally erasing internal contig/alignment edges.

    """

    trim = max(0, int(edge_trim))

    merged: Dict[str, List[Tuple[int, int]]] = {}

    for chrom, intervals in intervals_by_chrom.items():

        clean: List[Tuple[int, int]] = []

        for s, e in intervals:

            start = int(s) + trim

            end = int(e) - trim

            if end > start:

                clean.append((start, end))

        clean.sort()

        out: List[Tuple[int, int]] = []

        for start, end in clean:

            if not out or start > out[-1][1]:

                out.append((start, end))

            else:

                out[-1] = (out[-1][0], max(out[-1][1], end))

        merged[chrom] = out

    return merged


def reference_coverage_contains(

    coverage_by_chrom: Dict[str, List[Tuple[int, int]]],

    chrom: str,

    coord0: int,

) -> bool:

    """Return True when 0-based coord0 is inside a half-open coverage interval."""

    intervals = coverage_by_chrom.get(chrom, [])

    if not intervals:

        return False

    i = bisect_right(intervals, (int(coord0), 10**30)) - 1

    if i < 0:

        return False

    start, end = intervals[i]

    return start <= int(coord0) < end


def child_reference_anchors0(rec: VariantRecord) -> List[int]:

    """Return child reference anchor(s) in 0-based coordinates for coverage lookup.

    With --original-distance, use reconstructed per-haplotype original reference

    starts. Otherwise use the VCF POS converted from 1-based to 0-based.

    """

    if rec.original_loci:

        anchors = [max(0, int(loc.start) - 1) for loc in rec.original_loci]

        if anchors:

            return anchors

    return [max(0, int(rec.pos) - 1)]


def filter_child_records_to_query_core(

    records: Sequence[VariantRecord],

    query_bounds_by_source: Dict[str, Dict[str, Tuple[int, int]]],

    edge_trim: int,

    stats_by_source: Optional[Dict[str, Dict[str, int]]] = None,

) -> Tuple[List[VariantRecord], int, int]:

    """Remove child calls within ``edge_trim`` bp of query-contig ends.

    Trimming uses assembly query coordinates. Reference coverage remains an

    independent, untrimmed parental callability test at the child reference

    anchor.

    """

    trim = max(0, int(edge_trim))

    if trim == 0:

        return list(records), 0, 0

    kept: List[VariantRecord] = []

    excluded = 0

    missing = 0

    for rec in records:

        bounds = query_bounds_by_source.get(rec.source, {})

        if not rec.query_loci or any(

            locus.contig not in bounds for locus in rec.query_loci

        ):

            missing += 1

            if stats_by_source is not None:

                source_stats = stats_by_source.setdefault(rec.source, {})

                source_stats['records_excluded_missing_query_coordinates'] = (

                    source_stats.get(

                        'records_excluded_missing_query_coordinates', 0

                    ) + 1

                )

            continue

        inside = all(

            locus.start >= bounds[locus.contig][0] + trim

            and locus.end <= bounds[locus.contig][1] - trim

            for locus in rec.query_loci

        )

        if inside:

            kept.append(rec)

        else:

            excluded += 1

            if stats_by_source is not None:

                source_stats = stats_by_source.setdefault(rec.source, {})

                source_stats['records_excluded_edge_trim'] = (

                    source_stats.get('records_excluded_edge_trim', 0) + 1

                )

    return kept, excluded, missing


def filter_child_records_to_reference_core(

    records: Sequence[VariantRecord],

    coverage_by_source: Dict[str, Dict[str, List[Tuple[int, int]]]],

    stats_by_source: Optional[Dict[str, Dict[str, int]]] = None,

) -> Tuple[List[VariantRecord], int, int]:

    """Keep child calls whose reference anchor remains in trimmed referenceCoverage.

    ``coverage_by_source`` is expected to have already been trimmed by
    ``merge_reference_coverage_intervals(..., edge_trim=--trim-edges)``.
    This makes --trim-edges symmetric: child QUERYCOORD is checked against the
    assembly-contig core and the child reference anchor is checked against the
    reference-side alignment/coverage core.
    """

    kept: List[VariantRecord] = []

    excluded = 0

    missing = 0

    for rec in records:

        source_cov = coverage_by_source.get(rec.source, {})

        chrom_cov = source_cov.get(rec.chrom, [])

        if not chrom_cov:

            missing += 1

            if stats_by_source is not None:

                source_stats = stats_by_source.setdefault(rec.source, {})

                source_stats['records_excluded_missing_reference_coverage'] = (

                    source_stats.get('records_excluded_missing_reference_coverage', 0) + 1

                )

            continue

        anchors0 = child_reference_anchors0(rec)

        inside = bool(anchors0) and all(

            reference_coverage_contains(source_cov, rec.chrom, anchor0)

            for anchor0 in anchors0

        )

        if inside:

            kept.append(rec)

            continue

        excluded += 1

        if stats_by_source is not None:

            source_stats = stats_by_source.setdefault(rec.source, {})

            source_stats['records_excluded_reference_edge_trim'] = (

                source_stats.get('records_excluded_reference_edge_trim', 0) + 1

            )

            source_stats['records_excluded_edge_trim'] = (

                source_stats.get('records_excluded_edge_trim', 0) + 1

            )

    return kept, excluded, missing


def filter_child_records_by_pa_edge(

    records: Sequence[VariantRecord],

    query_bounds_by_source: Dict[str, Dict[str, Tuple[int, int]]],

    edge_distance: int,

    stats_by_source: Optional[Dict[str, Dict[str, int]]] = None,

) -> Tuple[List[VariantRecord], int, int]:

    """Exclude CHILD calls when any associated PA approaches a contig edge.

    PA intervals are reconstructed from child QUERYCOORD + LABEL_H.  For a

    ``+`` query interval, ``leftH_rightH`` gives

        PA = [QUERYCOORD.start - leftH, QUERYCOORD.end + rightH]

    and for ``-`` the H sides are swapped:

        PA = [QUERYCOORD.start - rightH, QUERYCOORD.end + leftH].

    A record is excluded when *any* reconstructed child PA has strictly less

    than ``edge_distance`` bp between either PA boundary and the corresponding

    child assembly contig edge.  Parents are never passed to this function.

    Records without usable PA metadata are retained, but counted separately.

    """

    distance = max(0, int(edge_distance))

    if distance == 0:

        return list(records), 0, 0

    kept: List[VariantRecord] = []

    excluded = 0

    missing = 0

    for rec in records:

        bounds = query_bounds_by_source.get(rec.source, {})

        usable = [

            locus for locus in rec.pa_loci

            if locus.contig in bounds

        ]

        if not usable or len(usable) != len(rec.pa_loci):

            missing += 1

            if stats_by_source is not None:

                source_stats = stats_by_source.setdefault(rec.source, {})

                source_stats['records_missing_pa_coordinates'] = (

                    source_stats.get('records_missing_pa_coordinates', 0) + 1

                )

            # PA-edge filtering is only defined for reconstructable PAs; do not

            # silently discard a call solely because optional PA metadata is absent.

            kept.append(rec)

            continue

        reaches_edge = any(

            locus.start < bounds[locus.contig][0] + distance

            or locus.end > bounds[locus.contig][1] - distance

            for locus in usable

        )

        if reaches_edge:

            excluded += 1

            if stats_by_source is not None:

                source_stats = stats_by_source.setdefault(rec.source, {})

                source_stats['records_excluded_pa_edge'] = (

                    source_stats.get('records_excluded_pa_edge', 0) + 1

                )

            continue

        kept.append(rec)

    return kept, excluded, missing


def open_text(path: str):

    if path.endswith(".gz"):

        return gzip.open(path, "rt")

    return open(path, "r")


def open_text_write(path: str):

    if path.endswith(".gz"):

        return gzip.open(path, "wt")

    return open(path, "w")


def read_query_path_contig_lengths(

    query_paths: str,

) -> Dict[str, Dict[str, int]]:

    """Read NAME FASTA [FAI] and return authoritative contig lengths."""

    listing = os.path.abspath(os.path.expanduser(query_paths))

    base = os.path.dirname(listing)

    by_sample: Dict[str, Dict[str, int]] = {}

    with open(listing, "rt", encoding="utf-8") as handle:

        for line_number, raw in enumerate(handle, 1):

            if not raw.strip() or raw.lstrip().startswith("#"):

                continue

            fields = shlex.split(raw, comments=True)

            if len(fields) not in {2, 3}:

                raise ValueError(

                    f"{listing}:{line_number}: expected NAME FASTA [FAI]"

                )

            sample, fasta_text = fields[:2]

            fasta = os.path.expanduser(fasta_text)

            if not os.path.isabs(fasta):

                fasta = os.path.join(base, fasta)

            fai_text = fields[2] if len(fields) == 3 else fasta + ".fai"

            fai = os.path.expanduser(fai_text)

            if not os.path.isabs(fai):

                fai = os.path.join(base, fai)

            if sample in by_sample:

                raise ValueError(

                    f"{listing}:{line_number}: duplicate sample {sample!r}"

                )

            if not os.path.isfile(fai):

                raise FileNotFoundError(

                    f"{listing}:{line_number}: FASTA index not found: {fai}"

                )

            lengths: Dict[str, int] = {}

            with open(fai, "rt", encoding="utf-8") as fai_handle:

                for fai_line_number, fai_raw in enumerate(fai_handle, 1):

                    columns = fai_raw.split()

                    if len(columns) < 2:

                        continue

                    try:

                        length = int(columns[1])

                    except ValueError:

                        raise ValueError(

                            f"{fai}:{fai_line_number}: invalid contig length"

                        ) from None

                    contig = columns[0]

                    if length < 0 or contig in lengths:

                        raise ValueError(

                            f"{fai}:{fai_line_number}: negative length or "

                            f"duplicate contig {contig!r}"

                        )

                    lengths[contig] = length

            if not lengths:

                raise ValueError(f"{fai}: no contig lengths found")

            by_sample[sample] = lengths

    if not by_sample:

        raise ValueError(f"{listing}: no assembly rows found")

    return by_sample


def query_bounds_for_vcf_samples(

    sample_names: Sequence[str],

    lengths_by_sample: Dict[str, Dict[str, int]],

    vcf_path: str,

) -> Dict[str, Tuple[int, int]]:

    """Resolve the query-contig bounds for one child haplotype VCF."""

    missing = [name for name in sample_names if name not in lengths_by_sample]

    if missing:

        raise ValueError(

            f"{vcf_path}: sample(s) absent from --query-paths: "

            + ",".join(missing)

        )

    bounds: Dict[str, Tuple[int, int]] = {}

    for sample in sample_names:

        for contig, length in lengths_by_sample[sample].items():

            value = (0, int(length))

            previous = bounds.get(contig)

            if previous is not None and previous != value:

                raise ValueError(

                    f"{vcf_path}: conflicting lengths for query contig "

                    f"{contig!r} across selected samples"

                )

            bounds[contig] = value

    return bounds


def parse_info(info: str) -> Dict[str, str]:

    out: Dict[str, str] = {}

    if not info or info == ".":

        return out

    for part in info.split(";"):

        if not part:

            continue

        m = INFO_RE.fullmatch(part)

        if not m:

            continue

        key = m.group(1)

        val = m.group(2)

        out[key] = "" if val is None else val

    return out


def parse_int(text: Optional[str]) -> Optional[int]:

    if text is None or text == "":

        return None

    if not INT_RE.match(text):

        return None

    return int(text)


def sample_prefix(name: str) -> str:

    return name.split("_", 1)[0]


def resolve_sample_indices(

    vcf_samples: Sequence[str],

    selector: str,

    prefix_match: bool,

) -> List[int]:

    if prefix_match:

        return [

            i

            for i, sample in enumerate(vcf_samples)

            if sample == selector or sample_prefix(sample) == selector

        ]

    return [i for i, sample in enumerate(vcf_samples) if sample == selector]


def check_no_role_overlap(

    child_idx: Sequence[int],

    mother_idx: Sequence[int],

    father_idx: Sequence[int],

    matched_samples: Dict[str, List[str]],

) -> None:

    if (set(child_idx) & set(mother_idx)) or (set(child_idx) & set(father_idx)) or (set(mother_idx) & set(father_idx)):

        raise SystemExit(

            "ERROR: child/mother/father sample columns overlap.\n"

            "This would make trio matching invalid and can cause artificial 100% match.\n"

            f"child_columns:  {','.join(matched_samples.get('child', []))}\n"

            f"mother_columns: {','.join(matched_samples.get('mother', []))}\n"

            f"father_columns: {','.join(matched_samples.get('father', []))}"

        )


def split_format_sample(sample_field: str, fmt_keys: Sequence[str]) -> List[str]:

    """Split a VCF sample field according to FORMAT keys.

    New graphreftovcf output uses decomposed FORMAT fields such as

    GT:TYPE:SIZE:EXTENDGRAPHCIGAR:ASSEMBLYCONTIG:QUERYCOORD:ALLELENAME:LABEL_H.

    EXTENDGRAPHCIGAR may itself contain ':' in named graph chunks.  When that

    happens, preserve all excess ':' pieces inside EXTENDGRAPHCIGAR by anchoring

    the fields before and after it.

    """

    if not fmt_keys:

        return []

    parts = sample_field.split(":")

    if len(parts) == len(fmt_keys):

        return parts

    cigar_key = "EXTENDGRAPHCIGAR" if "EXTENDGRAPHCIGAR" in fmt_keys else None

    if cigar_key is not None and len(parts) > len(fmt_keys):

        ci = fmt_keys.index(cigar_key)

        right_n = len(fmt_keys) - ci - 1

        if right_n == 0:

            return parts[:ci] + [":".join(parts[ci:])]

        if len(parts) >= ci + 1 + right_n:

            return parts[:ci] + [":".join(parts[ci:len(parts)-right_n])] + parts[len(parts)-right_n:]

    # Ordinary VCF fields should have exactly one value per FORMAT key.  For

    # malformed short rows, pad with missing values so callers fail softly.

    if len(parts) < len(fmt_keys):

        parts = parts + ["."] * (len(fmt_keys) - len(parts))

    return parts[:len(fmt_keys)]


def get_format_value(sample_field: str, fmt_keys: Sequence[str], key: str) -> Optional[str]:

    if key not in fmt_keys:

        return None

    idx = fmt_keys.index(key)

    parts = split_format_sample(sample_field, fmt_keys)

    if idx >= len(parts):

        return None

    return parts[idx]


def observation_values(value: Optional[str]) -> List[str]:

    """Return all per-observation values from one decomposed FORMAT field."""

    if value is None or value == "":

        return [""]

    return str(value).split(",")


def aligned_observation_values(*values: Optional[str]) -> List[Tuple[str, ...]]:

    """Transpose decomposed FORMAT columns without confusing ALT and calls.

    TYPE/SIZE/EXTENDGRAPHCIGAR are Number=. per-observation columns. Their comma

    index is an observation index, not the VCF ALT allele index. Scalar columns

    are broadcast; inconsistent non-scalar columns are rejected softly.

    """

    columns = [observation_values(value) for value in values]

    count = max((len(column) for column in columns), default=0)

    if count == 0 or any(len(column) not in {1, count} for column in columns):

        return []

    return [

        tuple(column[0] if len(column) == 1 else column[index] for column in columns)

        for index in range(count)

    ]


def get_hsv_value(sample_field: str, fmt_keys: Sequence[str]) -> Optional[str]:

    """Return the HSV value for graphcigartoref-style sample fields.

    These files commonly declare FORMAT=HSV while the HSV value itself contains

    colons, e.g. 1:DEL:-376:... . In that common case, the whole sample field

    is the HSV value and must not be truncated at the first colon.

    """

    if "HSV" not in fmt_keys:

        return None

    if len(fmt_keys) == 1 and fmt_keys[0] == "HSV":

        return sample_field

    hsv_idx = fmt_keys.index("HSV")

    if hsv_idx == len(fmt_keys) - 1:

        # HSV is colon-rich. If it is the last FORMAT item, split only the

        # preceding FORMAT fields and preserve the rest verbatim as HSV.

        parts = sample_field.split(":", hsv_idx)

        if len(parts) > hsv_idx:

            return parts[hsv_idx]

    return get_format_value(sample_field, fmt_keys, "HSV")


def hsv_has_allele(hsv_value: str, allele_index: int) -> bool:

    """

    Parse HSV values used by graphcigartoref-style VCFs.

    Examples:

      0                         absent

      .                         missing

      0|0                       absent

      1:DEL:-376:<3H3D...       ALT allele 1 present

      0|1:DEL:-376:...          ALT allele 1 present

      2:INS:...                 ALT allele 2 present

    """

    if not hsv_value or hsv_value == ".":

        return False

    target = str(allele_index)

    for token in hsv_value.split("|"):

        token = token.strip()

        if not token or token in {"0", "."}:

            continue

        allele = token.split(":", 1)[0]

        if allele == target:

            return True

    return False


def hsv_allele_type_contains(hsv_value: str, allele_index: int, text: str) -> bool:

    """Return True if the carried HSV allele has a type field containing text.

    For an HSV token like 1:DEL:-376:..., the type field is DEL.

    For 1:INS_xxx:..., this can be used to detect VarINS-like calls.

    """

    if not hsv_value or hsv_value == ".":

        return False

    target = str(allele_index)

    for token in hsv_value.split("|"):

        token = token.strip()

        if not token or token in {"0", "."}:

            continue

        parts = token.split(":")

        if len(parts) < 2:

            continue

        allele = parts[0]

        event_type = parts[1]

        if allele == target and text in event_type:

            return True

    return False


def vcf_unescape(text: str) -> str:

    """Undo the escaping used inside graphcigartoref HSV fields."""

    if text is None or text == ".":

        return ""

    return (

        str(text)

        .replace("@0A", "\n")

        .replace("@09", "\t")

        .replace("@2C", ",")

        .replace("@3B", ";")

        .replace("@3A", ":")

        .replace("@40", "@")

    )


def _looks_like_label_h_coord(text: str) -> bool:

    # One LABEL_H component is either one H coordinate (SNP-style) or a

    # left/right H pair (SV-style). Cross-graph provenance may join several

    # components with '&'.

    component = r"(?:[+-]?\d+H)(?:_?(?:[+-]?\d+H))?"

    return bool(re.fullmatch(rf"(?:{component})(?:&{component})*|\.", text or ""))


def parse_hsv_allele_fields(allele: str) -> Optional[List[str]]:

    """Parse GT:TYPE:SIZE:CIGAR:assemblycontig:query_range:allelename:label_H.

    The CIGAR/payload field may itself contain colons, so parse the first three

    fields from the left and the final metadata fields from the right. Both the

    current 8-field form and the legacy 7-field form are accepted. If LABEL_H

    is absent, empty, or '.', it is normalized to 0H_0H (zero left/right label-H distances).

    """

    first = allele.split(":", 3)

    if len(first) != 4:

        return None

    rest = first[3]

    # Current 8-field form. Also accept an explicitly empty/missing LABEL_H.

    last_new = rest.rsplit(":", 4)

    if len(last_new) == 5:

        label_h = last_new[-1].strip()

        if label_h in {"", "."} or _looks_like_label_h_coord(label_h):

            last_new[-1] = "0H_0H" if label_h in {"", "."} else label_h

            return first[:3] + last_new

    # Legacy 7-field form: LABEL_H is completely absent. Treat it as 0H_0H.

    last_old = rest.rsplit(":", 3)

    if len(last_old) == 4:

        return first[:3] + last_old + ["0H_0H"]

    return None


def parse_query_range_point(qrange: str) -> Optional[Tuple[int, int, str]]:

    """Parse an HSV query_range such as 14500-14876+ or 14500-."""

    text = vcf_unescape(qrange or "")

    m = re.match(r"^(\d+)(?:-(\d+))?([+-]?)$", text)

    if not m:

        return None

    a = int(m.group(1))

    b = int(m.group(2)) if m.group(2) is not None else a

    strand = m.group(3) or "+"

    if b < a:

        a, b = b, a

    return a, b, strand


# graphreftovcf writes each sample event from its own original reference start.

# If that start differs from the merged row POS, it prefixes an artificial

# standalone >NH or <NH chunk before the real sample CIGAR.  Therefore the

# pre-clustering reference start is recoverable from row POS + this prefix.

_CIGAR_OP_RE = re.compile(r"(\d+)([=MXIDHS])([A-Za-z]*)")

REF_CONSUMING_OPS = {"M", "=", "X", "D"}


def iter_top_level_cigar_chunks(encoded: str):

    """Yield (direction, chunk_text) for top-level >... / <... chunks."""

    if not encoded:

        return

    i = 0

    n = len(encoded)

    while i < n:

        if encoded[i] not in "><":

            i += 1

            continue

        direction = encoded[i]

        j = i + 1

        while j < n and encoded[j] not in "><":

            j += 1

        yield direction, encoded[i + 1:j]

        i = j


def split_row_position_padding(cigar: str) -> Tuple[str, str]:

    """Split graphreftovcf's artificial leading row-POS H padding.

    The writer only treats an H-only leading chunk as row-position padding when

    another top-level chunk follows immediately, e.g. >20H>335I... or

    <7H>426D....  This mirrors the writer's own split_row_position_padding().

    """

    prefix = ""

    body = cigar or ""

    while True:

        m = re.match(r"^[<>]\d+H(?=[<>])", body)

        if not m:

            break

        prefix += m.group(0)

        body = body[m.end():]

    return prefix, body


def row_padding_shift(prefix: str) -> int:

    """Convert artificial >NH/<NH row padding to signed reference shift."""

    shift = 0

    for direction, n_s in re.findall(r"([<>])(\d+)H", prefix or ""):

        n = int(n_s)

        shift += n if direction == ">" else -n

    return shift


def pairwise_ref_span(body: str) -> int:

    """Reference-consuming span in one pathless/main CIGAR body."""

    total = 0

    pos = 0

    while pos < len(body):

        if body[pos] in "><":

            break

        m = _CIGAR_OP_RE.match(body, pos)

        if not m:

            pos += 1

            continue

        n = int(m.group(1))

        op = m.group(2)

        if op in REF_CONSUMING_OPS:

            total += n

        pos = m.end()

    return total


def main_reference_span_from_cigar(cigar_body: str) -> int:

    """Recover the event span on the main/reference chromosome.

    Named chunks such as >NC_060932.1:... are encoded insertion mappings onto

    another reference/graph segment. They describe inserted query sequence and

    do not consume the main reference chromosome at the insertion anchor, so

    they are intentionally excluded. Pathless chunks are the main-reference

    stream and M/=X/D consume reference coordinates; H/I/S do not.

    """

    body = cigar_body or ""

    chunks = list(iter_top_level_cigar_chunks(body))

    if not chunks:

        return pairwise_ref_span(body)

    total = 0

    for _direction, chunk in chunks:

        if ":" in chunk:

            continue

        total += pairwise_ref_span(chunk)

    return total


def recover_original_reference_locus(

    row_chrom: str,

    row_pos: int,

    sample_name: str,

    cigar: str,

    event_type: str,

    size_text: str,

) -> Optional[OriginalLocus]:

    """Recover one sample/haplotype's pre-clustering reference interval.

    row_pos is the displayed merged VCF POS (the cluster's rounded mean start).

    The leading standalone H chunk stores sample_ref_start - row_pos.  The real

    pathless CIGAR then gives the main-reference span. For a pure insertion this

    span is zero, producing a point interval. A DEL-size fallback is used for

    legacy/simple records whose CIGAR lacks a parseable reference-consuming op.

    """

    if row_pos <= 0:

        return None

    prefix, body = split_row_position_padding(cigar or "")

    original_start = int(row_pos) + row_padding_shift(prefix)

    ref_span = main_reference_span_from_cigar(body)

    typ = (event_type or "").upper()

    if ref_span == 0 and (typ == "DEL" or typ.endswith("_DEL")):

        try:

            ref_span = abs(int(float(size_text)))

        except Exception:

            pass

    original_end = original_start + max(0, int(ref_span))

    return OriginalLocus(

        sample=sample_name,

        chrom=row_chrom,

        start=original_start,

        end=original_end,

    )


def hsv_original_loci(

    hsv_value: str,

    allele_index: int,

    sample_name: str,

    row_chrom: str,

    row_pos: int,

) -> List[OriginalLocus]:

    """Return original REFERENCE loci for one carried ALT allele.

    The original start is encoded by an optional leading row-position H padding

    in EXTENDGRAPHCIGAR. If that H padding is absent, its shift is zero, so the

    original start is simply the VCF row POS. LABEL_H is *not* required for

    this reconstruction. Legacy HSV alleles that omit LABEL_H are therefore

    benchmarkable and are treated as having zero LABEL_H.

    If the full legacy/new HSV metadata parser cannot split a token, fall back

    to the mandatory GT:TYPE:SIZE prefix and use row POS with zero H shift.

    This prevents a missing LABEL_H from discarding an otherwise valid allele.

    """

    if not hsv_value or hsv_value == "." or allele_index <= 0:

        return []

    target = str(allele_index)

    out: List[OriginalLocus] = []

    for token in hsv_value.split("|"):

        token = token.strip()

        if not token or token in {"0", "."}:

            continue

        # The allele number is always the first field and can be checked even

        # when trailing legacy metadata (including LABEL_H) is absent.

        prefix = token.split(":", 3)

        if len(prefix) < 3 or prefix[0] != target:

            continue

        vals = parse_hsv_allele_fields(token)

        if vals is not None and len(vals) >= 8:

            cigar = vals[3]

            event_type = vals[1]

            size_text = vals[2]

        else:

            # Missing/unparseable trailing HSV metadata must not make original

            # coordinates unavailable. The absence of row-position H padding

            # means H shift = 0; DEL span can still be recovered from SIZE.

            cigar = ""

            event_type = prefix[1]

            size_text = prefix[2]

        locus = recover_original_reference_locus(

            row_chrom=row_chrom,

            row_pos=row_pos,

            sample_name=sample_name,

            cigar=cigar,

            event_type=event_type,

            size_text=size_text,

        )

        if locus is not None:

            out.append(locus)

    return out


def sample_group_original_loci(

    sample_fields: Sequence[str],

    sample_names: Sequence[str],

    fmt_keys: Sequence[str],

    allele_index: int,

    row_chrom: str,

    row_pos: int,

) -> List[OriginalLocus]:

    """Collect reconstructed original reference loci from selected samples.

    Supports both legacy FORMAT=HSV payloads and the newer decomposed FORMAT:

    GT:TYPE:SIZE:EXTENDGRAPHCIGAR:ASSEMBLYCONTIG:QUERYCOORD:ALLELENAME:LABEL_H.

    The PA fields may both be absent; neither is needed for reconstruction.

    A missing leading H chunk in EXTENDGRAPHCIGAR means a row-position shift of 0.

    """

    out: List[OriginalLocus] = []

    if "HSV" in fmt_keys:

        for sample_field, sample_name in zip(sample_fields, sample_names):

            hsv = get_hsv_value(sample_field, fmt_keys)

            out.extend(

                hsv_original_loci(

                    hsv or "",

                    allele_index,

                    sample_name,

                    row_chrom=row_chrom,

                    row_pos=row_pos,

                )

            )

        return out

    # New decomposed graphreftovcf FORMAT.  The old code returned [] here for

    # every record because it required an HSV key even though all information

    # needed to reconstruct the reference locus is present explicitly.

    if "GT" not in fmt_keys or "EXTENDGRAPHCIGAR" not in fmt_keys:

        return out

    for sample_field, sample_name in zip(sample_fields, sample_names):

        gt = get_format_value(sample_field, fmt_keys, "GT")

        if not gt_has_allele(gt or "", allele_index):

            continue

        observations = aligned_observation_values(

            get_format_value(sample_field, fmt_keys, "EXTENDGRAPHCIGAR"),

            get_format_value(sample_field, fmt_keys, "TYPE"),

            get_format_value(sample_field, fmt_keys, "SIZE"),

        )

        for cigar, event_type, size_text in observations:

            # '.' / empty CIGAR means no explicit row-position H shift. For

            # deletions, recover_original_reference_locus() can use SIZE.

            if cigar in {"", "."}:

                cigar = ""

            locus = recover_original_reference_locus(

                row_chrom=row_chrom,

                row_pos=row_pos,

                sample_name=sample_name,

                cigar=vcf_unescape(cigar),

                event_type=event_type,

                size_text=size_text,

            )

            if locus is not None:

                out.append(locus)

    return out


def sample_group_query_loci(

    sample_fields: Sequence[str],

    fmt_keys: Sequence[str],

    allele_index: int,

) -> List[QueryLocus]:

    """Return carried observations' assembly query intervals.

    Current VCFs store the contig and interval in decomposed

    ``ASSEMBLYCONTIG`` and ``QUERYCOORD`` fields. Legacy HSV payloads carry

    the same values in fields five and six.

    """

    out: List[QueryLocus] = []

    if "HSV" in fmt_keys:

        target = str(allele_index)

        for sample_field in sample_fields:

            hsv = get_hsv_value(sample_field, fmt_keys) or ""

            for token in hsv.split("|"):

                values = parse_hsv_allele_fields(token.strip())

                if values is None or values[0] != target:

                    continue

                parsed = parse_query_range_point(values[5])

                contig = vcf_unescape(values[4])

                if parsed is not None and contig not in {"", "."}:

                    out.append(QueryLocus(contig, *parsed))

    elif (

        "GT" in fmt_keys

        and "ASSEMBLYCONTIG" in fmt_keys

        and "QUERYCOORD" in fmt_keys

    ):

        for sample_field in sample_fields:

            gt = get_format_value(sample_field, fmt_keys, "GT")

            if not gt_has_allele(gt or "", allele_index):

                continue

            observations = aligned_observation_values(

                get_format_value(sample_field, fmt_keys, "ASSEMBLYCONTIG"),

                get_format_value(sample_field, fmt_keys, "QUERYCOORD"),

            )

            for contig_text, coordinate_text in observations:

                contig = vcf_unescape(contig_text)

                parsed = parse_query_range_point(coordinate_text)

                if parsed is not None and contig not in {"", "."}:

                    out.append(QueryLocus(contig, *parsed))

    # Cross-graph observations can repeat the same assembly interval.

    return list(dict.fromkeys(out))

def parse_label_h_pairs(label_h: str) -> List[Tuple[int, int]]:

    """Parse SV LABEL_H components into ``(left_H, right_H)`` distances.

    Cross-graph provenance may contain multiple PA components joined by ``&``.

    Single-H SNP-style LABEL_H values do not define both PA boundaries and are

    therefore not used for PA-edge reconstruction.

    """

    text = vcf_unescape(label_h or "").strip()

    if not text or text == ".":

        return []

    out: List[Tuple[int, int]] = []

    for component in text.split("&"):

        component = component.strip()

        m = re.fullmatch(r"([+-]?\d+)H_?([+-]?\d+)H", component)

        if not m:

            continue

        out.append((int(m.group(1)), int(m.group(2))))

    return out


def reconstruct_pa_loci(

    contig: str,

    query_locus: QueryLocus,

    label_h: str,

) -> List[QueryLocus]:

    """Reconstruct PA query intervals from QUERYCOORD and LABEL_H."""

    out: List[QueryLocus] = []

    for left_h, right_h in parse_label_h_pairs(label_h):

        if query_locus.strand == "-":

            start = query_locus.start - right_h

            end = query_locus.end + left_h

        else:

            start = query_locus.start - left_h

            end = query_locus.end + right_h

        if end < start:

            start, end = end, start

        out.append(QueryLocus(contig, start, end, query_locus.strand))

    return out


def sample_group_pa_loci(

    sample_fields: Sequence[str],

    fmt_keys: Sequence[str],

    allele_index: int,

) -> List[QueryLocus]:

    """Return all reconstructable PAs for carried observations.

    This uses only child-side assembly metadata when called from separate-file

    mode.  Legacy HSV stores contig/query/LABEL_H in fields 5/6/8; decomposed

    FORMAT stores them as ASSEMBLYCONTIG, QUERYCOORD and LABEL_H.

    """

    out: List[QueryLocus] = []

    if "HSV" in fmt_keys:

        target = str(allele_index)

        for sample_field in sample_fields:

            hsv = get_hsv_value(sample_field, fmt_keys) or ""

            for token in hsv.split("|"):

                raw_token = token.strip()

                values = parse_hsv_allele_fields(raw_token)

                if values is None or values[0] != target:

                    continue

                # parse_hsv_allele_fields() intentionally normalizes an absent

                # legacy LABEL_H to 0H_0H for other coordinate logic.  For PA

                # reconstruction, however, 0H_0H is usable only when LABEL_H

                # was actually present in the input; do not invent a PA from

                # missing legacy metadata.

                first = raw_token.split(":", 3)

                if len(first) != 4:

                    continue

                last_new = first[3].rsplit(":", 4)

                if len(last_new) != 5:

                    continue

                raw_label_h = last_new[-1].strip()

                if raw_label_h in {"", "."} or not _looks_like_label_h_coord(raw_label_h):

                    continue

                parsed = parse_query_range_point(values[5])

                contig = vcf_unescape(values[4])

                if parsed is None or contig in {"", "."}:

                    continue

                qloc = QueryLocus(contig, *parsed)

                out.extend(reconstruct_pa_loci(contig, qloc, raw_label_h))

    elif (

        "GT" in fmt_keys

        and "ASSEMBLYCONTIG" in fmt_keys

        and "QUERYCOORD" in fmt_keys

        and "LABEL_H" in fmt_keys

    ):

        for sample_field in sample_fields:

            gt = get_format_value(sample_field, fmt_keys, "GT")

            if not gt_has_allele(gt or "", allele_index):

                continue

            observations = aligned_observation_values(

                get_format_value(sample_field, fmt_keys, "ASSEMBLYCONTIG"),

                get_format_value(sample_field, fmt_keys, "QUERYCOORD"),

                get_format_value(sample_field, fmt_keys, "LABEL_H"),

            )

            for contig_text, coordinate_text, label_h in observations:

                contig = vcf_unescape(contig_text)

                parsed = parse_query_range_point(coordinate_text)

                if parsed is None or contig in {"", "."}:

                    continue

                qloc = QueryLocus(contig, *parsed)

                out.extend(reconstruct_pa_loci(contig, qloc, label_h))

    return list(dict.fromkeys(out))


def _normalize_component_svtype(event_type: str) -> str:

    """Normalize decomposed observation TYPE values for DEL/INS matching."""

    typ = (event_type or "").strip().upper()

    if typ == "DEL" or typ.startswith("DEL_") or typ.endswith("_DEL"):

        return "DEL"

    if typ == "INS" or typ.startswith("INS_") or typ.endswith("_INS"):

        return "INS"

    return typ


def _parse_component_size(size_text: str) -> Optional[int]:

    try:

        return abs(int(float(str(size_text).strip())))

    except (TypeError, ValueError):

        return None


def _legacy_hsv_raw_label_h(raw_token: str) -> str:

    """Return an explicitly present legacy LABEL_H, otherwise an empty string."""

    first = raw_token.split(":", 3)

    if len(first) != 4:

        return ""

    last_new = first[3].rsplit(":", 4)

    if len(last_new) != 5:

        return ""

    raw_label_h = last_new[-1].strip()

    if raw_label_h in {"", "."} or not _looks_like_label_h_coord(raw_label_h):

        return ""

    return raw_label_h


def sample_group_sv_components(

    sample_fields: Sequence[str],

    sample_names: Sequence[str],

    fmt_keys: Sequence[str],

    allele_index: int,

    row_chrom: str,

    row_pos: int,

    use_original_coordinates: bool = False,

) -> List[SVComponent]:

    """Return carried DEL/INS observation components for one ALT allele.

    Decomposed graphreftovcf FORMAT columns are observation-indexed. A replacement
    allele can therefore contain (for example) TYPE=DEL,INS and SIZE=-120,85.
    Components are grouped by unique (TYPE, absolute SIZE), with corresponding
    reference/query/PA metadata aggregated onto each component. The caller decides
    whether row-level or component-level representation should be used.
    """

    observations = []

    target = str(allele_index)

    if "HSV" in fmt_keys:

        for sample_field, sample_name in zip(sample_fields, sample_names):

            hsv = get_hsv_value(sample_field, fmt_keys) or ""

            for token in hsv.split("|"):

                raw_token = token.strip()

                if not raw_token or raw_token in {"0", "."}:

                    continue

                prefix = raw_token.split(":", 3)

                if len(prefix) < 3 or prefix[0] != target:

                    continue

                values = parse_hsv_allele_fields(raw_token)

                event_type = values[1] if values is not None else prefix[1]

                size_text = values[2] if values is not None else prefix[2]

                typ = _normalize_component_svtype(event_type)

                size = _parse_component_size(size_text)

                if typ not in {"DEL", "INS"} or size is None:

                    continue

                cigar = values[3] if values is not None else ""

                original_locus = None

                if use_original_coordinates:

                    original_locus = recover_original_reference_locus(

                        row_chrom=row_chrom,

                        row_pos=row_pos,

                        sample_name=sample_name,

                        cigar=vcf_unescape(cigar),

                        event_type=event_type,

                        size_text=size_text,

                    )

                query_locus = None

                pa_loci: Tuple[QueryLocus, ...] = ()

                if values is not None and len(values) >= 8:

                    contig = vcf_unescape(values[4])

                    parsed = parse_query_range_point(values[5])

                    if parsed is not None and contig not in {"", "."}:

                        query_locus = QueryLocus(contig, *parsed)

                        raw_label_h = _legacy_hsv_raw_label_h(raw_token)

                        if raw_label_h:

                            pa_loci = tuple(reconstruct_pa_loci(contig, query_locus, raw_label_h))

                observations.append((typ, size, original_locus, query_locus, pa_loci))

    elif "GT" in fmt_keys and "TYPE" in fmt_keys and "SIZE" in fmt_keys:

        for sample_field, sample_name in zip(sample_fields, sample_names):

            gt = get_format_value(sample_field, fmt_keys, "GT")

            if not gt_has_allele(gt or "", allele_index):

                continue

            aligned = aligned_observation_values(

                get_format_value(sample_field, fmt_keys, "TYPE"),

                get_format_value(sample_field, fmt_keys, "SIZE"),

                get_format_value(sample_field, fmt_keys, "EXTENDGRAPHCIGAR"),

                get_format_value(sample_field, fmt_keys, "ASSEMBLYCONTIG"),

                get_format_value(sample_field, fmt_keys, "QUERYCOORD"),

                get_format_value(sample_field, fmt_keys, "LABEL_H"),

            )

            for event_type, size_text, cigar, contig_text, coordinate_text, label_h in aligned:

                typ = _normalize_component_svtype(event_type)

                size = _parse_component_size(size_text)

                if typ not in {"DEL", "INS"} or size is None:

                    continue

                original_locus = None

                if use_original_coordinates:

                    original_locus = recover_original_reference_locus(

                        row_chrom=row_chrom,

                        row_pos=row_pos,

                        sample_name=sample_name,

                        cigar="" if cigar in {"", "."} else vcf_unescape(cigar),

                        event_type=event_type,

                        size_text=size_text,

                    )

                query_locus = None

                pa_loci: Tuple[QueryLocus, ...] = ()

                contig = vcf_unescape(contig_text)

                parsed = parse_query_range_point(coordinate_text)

                if parsed is not None and contig not in {"", "."}:

                    query_locus = QueryLocus(contig, *parsed)

                    pa_loci = tuple(reconstruct_pa_loci(contig, query_locus, label_h))

                observations.append((typ, size, original_locus, query_locus, pa_loci))

    grouped = {}

    order = []

    for typ, size, original_locus, query_locus, pa_loci in observations:

        key = (typ, size)

        if key not in grouped:

            grouped[key] = {"original": [], "query": [], "pa": []}

            order.append(key)

        bucket = grouped[key]

        if original_locus is not None and original_locus not in bucket["original"]:

            bucket["original"].append(original_locus)

        if query_locus is not None and query_locus not in bucket["query"]:

            bucket["query"].append(query_locus)

        for pa_locus in pa_loci:

            if pa_locus not in bucket["pa"]:

                bucket["pa"].append(pa_locus)

    return [

        SVComponent(

            svtype=typ,

            size=size,

            original_loci=tuple(grouped[(typ, size)]["original"]),

            query_loci=tuple(grouped[(typ, size)]["query"]),

            pa_loci=tuple(grouped[(typ, size)]["pa"]),

        )

        for typ, size in order

    ]


def build_group_variant_records(

    *,

    chrom: str,

    pos: int,

    vid: str,

    ref: str,

    alt_allele: str,

    info: Dict[str, str],

    allele_index: int,

    source: str,

    line_no: int,

    sample_fields: Sequence[str],

    sample_names: Sequence[str],

    fmt_keys: Sequence[str],

    min_size: int,

    include_all: bool,

    use_original_coordinates: bool,

    include_query_metadata: bool = False,

) -> Tuple[List[VariantRecord], int]:

    """Build carried records for one role/sample group and one ALT allele.

    Mixed DEL+INS alleles are atomized before the size cutoff, so a replacement
    can contribute two SVs when both components independently exceed --min-size.
    The second return value counts otherwise benchmarkable components/records
    dropped only because original reference coordinates could not be recovered.
    """

    if not sample_group_has_allele(sample_fields, fmt_keys, allele_index):

        return [], 0

    row_svtype = infer_svtype(ref, alt_allele, info)

    components = sample_group_sv_components(

        sample_fields,

        sample_names,

        fmt_keys,

        allele_index,

        chrom,

        pos,

        use_original_coordinates=use_original_coordinates,

    )

    component_types = {component.svtype for component in components}

    use_components = bool(components) and (

        len(component_types) > 1

        or row_svtype not in {"DEL", "INS"}

        or row_svtype not in component_types

    )

    if use_components:

        out: List[VariantRecord] = []

        missing_original = 0

        for component_number, component in enumerate(components, 1):

            variant_class = classify_variant(

                ref, alt_allele, component.svtype, component.size, min_size

            )

            if not include_all and (

                component.size <= min_size or variant_class != "SV"

            ):

                continue

            if use_original_coordinates and not component.original_loci:

                missing_original += 1

                continue

            component_end = pos + max(component.size, 1) - 1

            base_vid = vid if vid and vid != "." else f"{chrom}:{pos}:{allele_index}"

            out.append(VariantRecord(

                chrom=chrom,

                pos=pos,

                end=component_end,

                svtype=component.svtype,

                size=component.size,

                vid=f"{base_vid}:component{component_number}:{component.svtype}:{component.size}",

                ref=ref,

                alt=alt_allele,

                allele_index=allele_index,

                variant_class=variant_class,

                source=source,

                line_no=line_no,

                original_loci=component.original_loci,

                query_loci=component.query_loci if include_query_metadata else (),

                pa_loci=component.pa_loci if include_query_metadata else (),

            ))

        return out, missing_original

    svtype = row_svtype

    end, size = infer_end_and_size_for_allele(

        pos, ref, alt_allele, info, allele_index

    )

    variant_class = classify_variant(ref, alt_allele, svtype, size, min_size)

    if not include_all:

        if size is None or end is None or size <= min_size or variant_class != "SV":

            return [], 0

    else:

        if end is None:

            if size is None:

                end = pos + max(len(ref), len(alt_allele)) - 1

            else:

                end = pos + max(abs(size), 1) - 1

        if size is None:

            size = max(len(ref), len(alt_allele))

    if end is None or size is None:

        return [], 0

    original_loci = []

    if use_original_coordinates:

        original_loci = sample_group_original_loci(

            sample_fields, sample_names, fmt_keys, allele_index, chrom, pos

        )

        if not original_loci:

            return [], 1

    query_loci = (

        sample_group_query_loci(sample_fields, fmt_keys, allele_index)

        if include_query_metadata else []

    )

    pa_loci = (

        sample_group_pa_loci(sample_fields, fmt_keys, allele_index)

        if include_query_metadata else []

    )

    return [VariantRecord(

        chrom=chrom,

        pos=pos,

        end=end,

        svtype=svtype,

        size=size,

        vid=vid if vid and vid != "." else f"{chrom}:{pos}:{svtype}:{size}:{allele_index}",

        ref=ref,

        alt=alt_allele,

        allele_index=allele_index,

        variant_class=variant_class,

        source=source,

        line_no=line_no,

        original_loci=tuple(original_loci),

        query_loci=tuple(query_loci),

        pa_loci=tuple(pa_loci),

    )], 0


def gt_has_allele(gt_value: str, allele_index: int) -> bool:

    if not gt_value or gt_value in {".", "./.", ".|."}:

        return False

    target = str(allele_index)

    alleles = [a for a in GT_SEPS.split(gt_value) if a not in {"", "."}]

    return any(a == target for a in alleles)


def sample_group_has_allele(

    sample_fields: Sequence[str],

    fmt_keys: Sequence[str],

    allele_index: int,

) -> bool:

    """

    Return True if at least one selected sample column carries allele_index.

    Supports FORMAT containing GT or HSV. If neither is declared, it falls back

    to HSV-like parsing because this is safer for graphcigartoref output than

    treating missing GT as positive.

    """

    if not sample_fields or allele_index <= 0:

        return False

    has_gt = "GT" in fmt_keys

    has_hsv = "HSV" in fmt_keys

    for sample in sample_fields:

        if not sample or sample == ".":

            continue

        if has_gt:

            gt = get_format_value(sample, fmt_keys, "GT")

            if gt_has_allele(gt or "", allele_index):

                return True

        elif has_hsv:

            # In FORMAT=HSV, the whole sample field is the HSV value.

            # If FORMAT=HSV:OTHER, take only the HSV component.

            hsv = get_hsv_value(sample, fmt_keys)

            if hsv_has_allele(hsv or "", allele_index):

                return True

        else:

            # Fallback for non-standard files with no GT and no declared HSV.

            # Only explicit allele tokens count as positive.

            if hsv_has_allele(sample, allele_index):

                return True

    return False


def sample_value_is_missing(sample_field: str, fmt_keys: Sequence[str]) -> bool:

    """Return True when a selected sample field is missing/no-call.

    For HSV-style files, a literal ``.`` is treated as missing rather than as

    a confirmed non-carrier. For GT-style files, full or partial missing GTs

    such as ``.``, ``./.``, ``.|.``, and ``0/.`` are also treated as missing.

    """

    if not sample_field or sample_field == ".":

        return True

    has_gt = "GT" in fmt_keys

    has_hsv = "HSV" in fmt_keys

    if has_gt:

        gt = get_format_value(sample_field, fmt_keys, "GT")

        if gt is None or gt == "":

            return True

        alleles = GT_SEPS.split(gt)

        return any(a in {"", "."} for a in alleles)

    if has_hsv:

        hsv = get_hsv_value(sample_field, fmt_keys)

        if not hsv or hsv == ".":

            return True

        for token in hsv.split("|"):

            token = token.strip()

            if not token:

                return True

            allele = token.split(":", 1)[0]

            if allele == ".":

                return True

        return False

    # Fallback for non-standard files with no GT and no declared HSV key.

    if sample_field == ".":

        return True

    for token in sample_field.split("|"):

        token = token.strip()

        if not token:

            return True

        allele = token.split(":", 1)[0]

        if allele == ".":

            return True

    return False


def sample_group_has_missing(

    sample_fields: Sequence[str],

    fmt_keys: Sequence[str],

) -> bool:

    """Return True if any selected sample column has a missing/no-call value."""

    if not sample_fields:

        return True

    return any(sample_value_is_missing(sample, fmt_keys) for sample in sample_fields)


def sample_group_hsv_type_contains(

    sample_fields: Sequence[str],

    fmt_keys: Sequence[str],

    allele_index: int,

    text: str,

) -> bool:

    """Return True if a carried allele's TYPE contains text.

    Supports both legacy HSV payloads and the new decomposed GT:TYPE:... FORMAT.

    """

    if not sample_fields or allele_index <= 0:

        return False

    has_gt = "GT" in fmt_keys

    has_hsv = "HSV" in fmt_keys

    has_type = "TYPE" in fmt_keys

    for sample in sample_fields:

        if not sample or sample == ".":

            continue

        if has_hsv:

            hsv = get_hsv_value(sample, fmt_keys)

            if hsv_allele_type_contains(hsv or "", allele_index, text):

                return True

        elif has_gt and has_type:

            gt = get_format_value(sample, fmt_keys, "GT")

            if not gt_has_allele(gt or "", allele_index):

                continue

            event_types = observation_values(

                get_format_value(sample, fmt_keys, "TYPE")

            )

            if any(text in event_type for event_type in event_types):

                return True

        elif not has_gt:

            # Fallback for non-standard HSV-like files with no declared HSV key.

            if hsv_allele_type_contains(sample, allele_index, text):

                return True

    return False


def infer_svtype(ref: str, alt: str, info: Dict[str, str]) -> str:

    if info.get("SVTYPE"):

        return info["SVTYPE"].upper()

    m = SVTYPE_FROM_ALT.match(alt)

    if m:

        return m.group(1).upper()

    if len(ref) > len(alt):

        return "DEL"

    if len(alt) > len(ref):

        return "INS"

    return "UNK"


def pick_info_value(info: Dict[str, str], key: str, allele_index: int) -> Optional[str]:

    val = info.get(key)

    if val is None or val == "":

        return None

    parts = val.split(",")

    if len(parts) == 1:

        return parts[0]

    idx = allele_index - 1

    if 0 <= idx < len(parts):

        return parts[idx]

    return parts[0]


def infer_end_and_size(

    pos: int,

    ref: str,

    alt: str,

    info: Dict[str, str],

) -> Tuple[Optional[int], Optional[int]]:

    end = parse_int(info.get("END"))

    svlen = parse_int(info.get("SVLEN"))

    if svlen is None and "SVLEN" in info and "," in info["SVLEN"]:

        svlen = parse_int(info["SVLEN"].split(",", 1)[0])

    if end is None and svlen is not None:

        end = pos + max(abs(svlen), 1) - 1

    if svlen is None and end is not None:

        svlen = end - pos + 1

    if svlen is None:

        if alt.startswith("<") and alt.endswith(">"):

            svlen = None

        else:

            svlen = len(alt) - len(ref)

            if svlen == 0:

                svlen = None

    if end is None and svlen is not None:

        end = pos + max(abs(svlen), 1) - 1

    size = abs(svlen) if svlen is not None else None

    return end, size


def infer_end_and_size_for_allele(

    pos: int,

    ref: str,

    alt: str,

    info: Dict[str, str],

    allele_index: int,

) -> Tuple[Optional[int], Optional[int]]:

    allele_info = dict(info)

    for key in ("END", "SVLEN"):

        val = pick_info_value(info, key, allele_index)

        if val is not None:

            allele_info[key] = val

    return infer_end_and_size(pos, ref, alt, allele_info)


def classify_variant(ref: str, alt: str, svtype: str, size: Optional[int], min_size: int) -> str:

    if len(ref) == 1 and len(alt) == 1 and not alt.startswith("<") and ref != alt:

        return "SNP"

    if size is not None and size > min_size:

        return "SV"

    if not alt.startswith("<"):

        return "INDEL"

    if svtype != "UNK":

        return "SV"

    return "OTHER"


def load_trio_records(

    vcf_path: str,

    child_selector: str,

    mother_selector: str,

    father_selector: str,

    min_size: int,

    include_all: bool,

    prefix_match: bool,

    exclude_child_varins: bool = True,

    use_original_coordinates: bool = False,

) -> Tuple[

    List[VariantRecord],

    List[VariantRecord],

    List[VariantRecord],

    Dict[str, List[str]],

    int,

    List[str],

    Dict[int, str],

    Dict[str, int],

]:

    child_records: List[VariantRecord] = []

    mother_records: List[VariantRecord] = []

    father_records: List[VariantRecord] = []

    matched_samples: Dict[str, List[str]] = {}

    child_idx_vcf: List[int] = []

    mother_idx_vcf: List[int] = []

    father_idx_vcf: List[int] = []

    selected_col_indices: List[int] = []

    max_needed_col = 8

    saw_header = False

    total_records = 0

    header_lines: List[str] = []

    raw_by_line: Dict[int, str] = {}

    parent_missing_stats: Dict[str, int] = {

        "excluded_child_records": 0,

        "mother_missing": 0,

        "father_missing": 0,

        "both_missing": 0,

        "child_missing_original": 0,

        "mother_missing_original": 0,

        "father_missing_original": 0,

    }

    with open_text(vcf_path) as fh:

        for line_no, raw in enumerate(fh, 1):

            if raw.startswith("#"):

                header_lines.append(raw.rstrip("\n"))

            if raw.startswith("##"):

                continue

            if raw.startswith("#CHROM"):

                saw_header = True

                header_fields = raw.rstrip("\n").split("\t")

                if len(header_fields) < 10:

                    raise SystemExit("ERROR: VCF has no sample columns.")

                samples = header_fields[9:]

                child_idx_sample = resolve_sample_indices(samples, child_selector, prefix_match)

                mother_idx_sample = resolve_sample_indices(samples, mother_selector, prefix_match)

                father_idx_sample = resolve_sample_indices(samples, father_selector, prefix_match)

                matched_samples["child"] = [samples[i] for i in child_idx_sample]

                matched_samples["mother"] = [samples[i] for i in mother_idx_sample]

                matched_samples["father"] = [samples[i] for i in father_idx_sample]

                missing = [

                    role

                    for role, idxs in [

                        ("child", child_idx_sample),

                        ("mother", mother_idx_sample),

                        ("father", father_idx_sample),

                    ]

                    if not idxs

                ]

                if missing:

                    raise SystemExit(

                        "ERROR: No VCF sample columns matched: "

                        + ", ".join(missing)

                        + "\nDefault matching is exact. Use --prefix-match for h1/h2 names."

                    )

                check_no_role_overlap(

                    child_idx_sample,

                    mother_idx_sample,

                    father_idx_sample,

                    matched_samples,

                )

                child_idx_vcf = [9 + i for i in child_idx_sample]

                mother_idx_vcf = [9 + i for i in mother_idx_sample]

                father_idx_vcf = [9 + i for i in father_idx_sample]

                selected_col_indices = child_idx_vcf + mother_idx_vcf + father_idx_vcf

                max_needed_col = max([8] + selected_col_indices)

                continue

            if raw.startswith("#"):

                continue

            if not saw_header:

                raise SystemExit("ERROR: VCF header line #CHROM was not found before records.")

            total_records += 1

            # Split only up to the last sample column needed. This matters for huge VCFs.

            fields = raw.rstrip("\n").split("\t", max_needed_col + 1)

            if len(fields) <= max_needed_col or len(fields) < 9:

                continue

            chrom, pos_s, vid, ref, alt, _qual, _flt, info_s = fields[:8]

            try:

                pos = int(pos_s)

            except ValueError:

                continue

            info = parse_info(info_s)

            fmt_keys = fields[8].split(":") if fields[8] else []

            child_group = [fields[i] for i in child_idx_vcf]

            mother_group = [fields[i] for i in mother_idx_vcf]

            father_group = [fields[i] for i in father_idx_vcf]

            for allele_index, alt_allele in enumerate(alt.split(","), 1):

                child_has = sample_group_has_allele(child_group, fmt_keys, allele_index)

                mother_has = sample_group_has_allele(mother_group, fmt_keys, allele_index)

                father_has = sample_group_has_allele(father_group, fmt_keys, allele_index)

                child_excluded_varins = (

                    child_has

                    and exclude_child_varins

                    and sample_group_hsv_type_contains(child_group, fmt_keys, allele_index, "INS_")

                )

                child_role_records: List[VariantRecord] = []

                child_missing_original_count = 0

                if child_has and not child_excluded_varins:

                    child_role_records, child_missing_original_count = build_group_variant_records(

                        chrom=chrom,

                        pos=pos,

                        vid=vid,

                        ref=ref,

                        alt_allele=alt_allele,

                        info=info,

                        allele_index=allele_index,

                        source=os.path.basename(vcf_path),

                        line_no=line_no,

                        sample_fields=child_group,

                        sample_names=matched_samples["child"],

                        fmt_keys=fmt_keys,

                        min_size=min_size,

                        include_all=include_all,

                        use_original_coordinates=use_original_coordinates,

                    )

                mother_role_records, mother_missing_original_count = build_group_variant_records(

                    chrom=chrom,

                    pos=pos,

                    vid=vid,

                    ref=ref,

                    alt_allele=alt_allele,

                    info=info,

                    allele_index=allele_index,

                    source=os.path.basename(vcf_path),

                    line_no=line_no,

                    sample_fields=mother_group,

                    sample_names=matched_samples["mother"],

                    fmt_keys=fmt_keys,

                    min_size=min_size,

                    include_all=include_all,

                    use_original_coordinates=use_original_coordinates,

                )

                father_role_records, father_missing_original_count = build_group_variant_records(

                    chrom=chrom,

                    pos=pos,

                    vid=vid,

                    ref=ref,

                    alt_allele=alt_allele,

                    info=info,

                    allele_index=allele_index,

                    source=os.path.basename(vcf_path),

                    line_no=line_no,

                    sample_fields=father_group,

                    sample_names=matched_samples["father"],

                    fmt_keys=fmt_keys,

                    min_size=min_size,

                    include_all=include_all,

                    use_original_coordinates=use_original_coordinates,

                )

                if child_missing_original_count:

                    parent_missing_stats["child_missing_original"] += child_missing_original_count

                if mother_missing_original_count:

                    parent_missing_stats["mother_missing_original"] += mother_missing_original_count

                if father_missing_original_count:

                    parent_missing_stats["father_missing_original"] += father_missing_original_count

                mother_missing = sample_group_has_missing(mother_group, fmt_keys)

                father_missing = sample_group_has_missing(father_group, fmt_keys)

                # Preserve the historical conservative rule: if a carried parent
                # component should be benchmarkable but lacks original coordinates,
                # treat that parent as missing for child benchmarking.

                mother_original_missing = (

                    use_original_coordinates and mother_has and mother_missing_original_count > 0

                )

                father_original_missing = (

                    use_original_coordinates and father_has and father_missing_original_count > 0

                )

                child_excluded_parent_missing = bool(child_role_records) and (

                    mother_missing or father_missing or mother_original_missing or father_original_missing

                )

                if child_excluded_parent_missing:

                    excluded_n = len(child_role_records)

                    parent_missing_stats["excluded_child_records"] += excluded_n

                    if mother_missing or mother_original_missing:

                        parent_missing_stats["mother_missing"] += excluded_n

                    if father_missing or father_original_missing:

                        parent_missing_stats["father_missing"] += excluded_n

                    if (mother_missing or mother_original_missing) and (father_missing or father_original_missing):

                        parent_missing_stats["both_missing"] += excluded_n

                else:

                    child_records.extend(child_role_records)

                    if child_role_records:

                        # Keep only rows actually carried by the child and included in
                        # the benchmark, so --fp/--tp follows the same denominator.

                        raw_by_line.setdefault(line_no, raw.rstrip("\n"))

                mother_records.extend(mother_role_records)

                father_records.extend(father_role_records)

    if not saw_header:

        raise SystemExit("ERROR: VCF header line #CHROM was not found.")

    return child_records, mother_records, father_records, matched_samples, total_records, header_lines, raw_by_line, parent_missing_stats


def size_ratio_ok(size_a: int, size_b: int, min_ratio: float) -> bool:

    if size_a <= 0 or size_b <= 0:

        return False

    return min(size_a, size_b) / max(size_a, size_b) >= min_ratio


def endpoints_close(a: VariantRecord, b: VariantRecord, slop: int) -> bool:

    return abs(a.pos - b.pos) <= slop and abs(a.end - b.end) <= slop


OriginalDistance = Union[int, str]


def match_sv(child: VariantRecord, parent: VariantRecord, slop: int, min_ratio: float) -> bool:

    if child.chrom != parent.chrom:

        return False

    if child.svtype != "UNK" and parent.svtype != "UNK" and child.svtype != parent.svtype:

        return False

    if not endpoints_close(child, parent, slop):

        return False

    if not size_ratio_ok(child.size, parent.size, min_ratio):

        return False

    return True


def original_interval_match(a: OriginalLocus, b: OriginalLocus, distance: OriginalDistance) -> bool:

    """Compare reconstructed pre-clustering reference intervals."""

    if a.chrom != b.chrom:

        return False

    if distance == "exact":

        return a.start == b.start and a.end == b.end

    max_gap = int(distance)

    if a.end >= b.start and b.end >= a.start:

        return True

    if a.end < b.start:

        gap = b.start - a.end

    else:

        gap = a.start - b.end

    return gap <= max_gap


def match_sv_original(

    child: VariantRecord,

    parent: VariantRecord,

    original_distance: OriginalDistance,

    min_ratio: float,

) -> bool:

    if child.chrom != parent.chrom:

        return False

    if child.svtype != "UNK" and parent.svtype != "UNK" and child.svtype != parent.svtype:

        return False

    if not size_ratio_ok(child.size, parent.size, min_ratio):

        return False

    if not child.original_loci or not parent.original_loci:

        return False

    return any(

        original_interval_match(child_locus, parent_locus, original_distance)

        for child_locus in child.original_loci

        for parent_locus in parent.original_loci

    )


def find_parent_matches(

    child_records: List[VariantRecord],

    parent_records: List[VariantRecord],

    slop: int,

    min_ratio: float,

    original_distance: Optional[OriginalDistance] = None,

) -> List[bool]:

    exact_set = set()

    if original_distance is not None:

        # Original coordinates are reconstructed on the reference chromosome

        # from per-sample HSV CIGAR row-position padding. Index parental SV

        # intervals into 100-kb bins so large trio VCFs do not require an

        # O(child_SVs * parent_SVs_per_chromosome) scan. The final acceptance

        # test remains match_sv_original(), so matching semantics are unchanged.

        original_bin_size = 100_000

        parent_bins: Dict[Tuple[str, int], List[VariantRecord]] = defaultdict(list)

        for rec in parent_records:

            if rec.variant_class in {"SNP", "INDEL"}:

                exact_set.add((rec.chrom, rec.pos, rec.ref, rec.alt, rec.variant_class))

                continue

            for locus in rec.original_loci:

                b0 = max(0, locus.start) // original_bin_size

                b1 = max(0, locus.end) // original_bin_size

                for b in range(b0, b1 + 1):

                    parent_bins[(locus.chrom, b)].append(rec)

        out: List[bool] = []

        gap = 0 if original_distance == "exact" else int(original_distance)

        for child in child_records:

            if child.variant_class in {"SNP", "INDEL"}:

                key = (child.chrom, child.pos, child.ref, child.alt, child.variant_class)

                out.append(key in exact_set)

                continue

            candidates: List[VariantRecord] = []

            seen_ids = set()

            for locus in child.original_loci:

                qstart = max(0, locus.start - gap)

                qend = max(qstart, locus.end + gap)

                b0 = qstart // original_bin_size

                b1 = qend // original_bin_size

                for b in range(b0, b1 + 1):

                    for parent in parent_bins.get((locus.chrom, b), []):

                        parent_id = id(parent)

                        if parent_id not in seen_ids:

                            seen_ids.add(parent_id)

                            candidates.append(parent)

            matched = any(

                match_sv_original(child, parent, original_distance, min_ratio)

                for parent in candidates

            )

            out.append(matched)

        return out

    bin_size = max(1, slop)

    sv_bins: Dict[Tuple[str, int], List[VariantRecord]] = defaultdict(list)

    for rec in parent_records:

        if rec.variant_class in {"SNP", "INDEL"}:

            exact_set.add((rec.chrom, rec.pos, rec.ref, rec.alt, rec.variant_class))

        else:

            sv_bins[(rec.chrom, rec.pos // bin_size)].append(rec)

    out: List[bool] = []

    for child in child_records:

        if child.variant_class in {"SNP", "INDEL"}:

            key = (child.chrom, child.pos, child.ref, child.alt, child.variant_class)

            out.append(key in exact_set)

            continue

        center_bin = child.pos // bin_size

        matched = False

        for b in range(center_bin - 1, center_bin + 2):

            for parent in sv_bins.get((child.chrom, b), []):

                if match_sv(child, parent, slop, min_ratio):

                    matched = True

                    break

            if matched:

                break

        out.append(matched)

    return out


def _unique_join(values: Sequence[str]) -> str:

    seen = set()

    out = []

    for value in values:

        value = str(value)

        if value not in seen:

            seen.add(value)

            out.append(value)

    return ",".join(out)


def add_info_annotations(raw_vcf_row: str, annotations: Sequence[str]) -> str:

    fields = raw_vcf_row.split("\t")

    if len(fields) < 8:

        return raw_vcf_row

    extra = ";".join(a for a in annotations if a)

    if not extra:

        return raw_vcf_row

    if not fields[7] or fields[7] == ".":

        fields[7] = extra

    else:

        fields[7] = fields[7] + ";" + extra

    return "\t".join(fields)


def write_fp_vcf(

    outpath: str,

    header_lines: Sequence[str],

    raw_by_line: Dict[int, str],

    child_records: Sequence[VariantRecord],

    mom_hits: Sequence[bool],

    dad_hits: Sequence[bool],

) -> Tuple[int, int]:

    """

    Write child-only variants to a VCF.

    Definition of FP here: a child-carried allele that is not matched in mother

    and not matched in father under the same trio-check matching rules.

    For multi-allelic rows, this writes the original VCF row once and annotates

    which child ALT allele index or indices are child-only. This avoids changing

    ALT/sample fields and preserves VCF validity for HSV fields.

    """

    fp_by_line: Dict[int, List[VariantRecord]] = defaultdict(list)

    for rec, mom_hit, dad_hit in zip(child_records, mom_hits, dad_hits):

        if not mom_hit and not dad_hit:

            fp_by_line[rec.line_no].append(rec)

    fp_row_count = 0

    fp_allele_count = sum(len(v) for v in fp_by_line.values())

    with open_text_write(outpath) as out:

        inserted_info_headers = False

        for header in header_lines:

            if header.startswith("#CHROM") and not inserted_info_headers:

                out.write('##INFO=<ID=TRIOCHECK_FP,Number=0,Type=Flag,Description="Child-carried allele not matched in either parent by Triocheck">\n')

                out.write('##INFO=<ID=TRIOCHECK_FP_ALLELES,Number=.,Type=Integer,Description="ALT allele indices that are child-only under Triocheck">\n')

                out.write('##INFO=<ID=TRIOCHECK_FP_TYPES,Number=.,Type=String,Description="SVTYPE values for child-only alleles">\n')

                out.write('##INFO=<ID=TRIOCHECK_FP_CLASSES,Number=.,Type=String,Description="Variant classes for child-only alleles">\n')

                out.write('##INFO=<ID=TRIOCHECK_FP_LINES,Number=1,Type=Integer,Description="Input VCF line number reported by Triocheck">\n')

                inserted_info_headers = True

            out.write(header + "\n")

        if not inserted_info_headers:

            out.write('##INFO=<ID=TRIOCHECK_FP,Number=0,Type=Flag,Description="Child-carried allele not matched in either parent by Triocheck">\n')

            out.write('##INFO=<ID=TRIOCHECK_FP_ALLELES,Number=.,Type=Integer,Description="ALT allele indices that are child-only under Triocheck">\n')

            out.write('##INFO=<ID=TRIOCHECK_FP_TYPES,Number=.,Type=String,Description="SVTYPE values for child-only alleles">\n')

            out.write('##INFO=<ID=TRIOCHECK_FP_CLASSES,Number=.,Type=String,Description="Variant classes for child-only alleles">\n')

            out.write('##INFO=<ID=TRIOCHECK_FP_LINES,Number=1,Type=Integer,Description="Input VCF line number reported by Triocheck">\n')

        for line_no in sorted(fp_by_line):

            raw = raw_by_line.get(line_no)

            if raw is None:

                continue

            recs = fp_by_line[line_no]

            annotations = [

                "TRIOCHECK_FP",

                "TRIOCHECK_FP_ALLELES=" + _unique_join(str(r.allele_index) for r in recs),

                "TRIOCHECK_FP_TYPES=" + _unique_join(r.svtype for r in recs),

                "TRIOCHECK_FP_CLASSES=" + _unique_join(r.variant_class for r in recs),

                f"TRIOCHECK_FP_LINES={line_no}",

            ]

            out.write(add_info_annotations(raw, annotations) + "\n")

            fp_row_count += 1

    return fp_row_count, fp_allele_count



def write_tp_vcf(

    outpath: str,

    header_lines: Sequence[str],

    raw_by_line: Dict[int, str],

    child_records: Sequence[VariantRecord],

    mom_hits: Sequence[bool],

    dad_hits: Sequence[bool],

) -> Tuple[int, int]:

    """

    Write true-positive/inherited child variants to a VCF.

    Definition of TP here: a child-carried allele that is matched in mother

    and/or father under the same trio-check matching rules used in the summary.

    For multi-allelic rows, this writes the original VCF row once and annotates

    which child ALT allele index or indices are true positives. It also records

    which TP alleles matched the mother and which matched the father. This avoids

    changing ALT/sample fields and preserves VCF validity for HSV fields.

    """

    tp_by_line: Dict[int, List[Tuple[VariantRecord, bool, bool]]] = defaultdict(list)

    for rec, mom_hit, dad_hit in zip(child_records, mom_hits, dad_hits):

        if mom_hit or dad_hit:

            tp_by_line[rec.line_no].append((rec, mom_hit, dad_hit))

    tp_row_count = 0

    tp_allele_count = sum(len(v) for v in tp_by_line.values())

    tp_headers = [

        '##INFO=<ID=TRIOCHECK_TP,Number=0,Type=Flag,Description="Child-carried allele matched in at least one parent by Triocheck">',

        '##INFO=<ID=TRIOCHECK_TP_ALLELES,Number=.,Type=Integer,Description="ALT allele indices matched in at least one parent under Triocheck">',

        '##INFO=<ID=TRIOCHECK_TP_MOTHER_ALLELES,Number=.,Type=Integer,Description="TP ALT allele indices matched in mother under Triocheck">',

        '##INFO=<ID=TRIOCHECK_TP_FATHER_ALLELES,Number=.,Type=Integer,Description="TP ALT allele indices matched in father under Triocheck">',

        '##INFO=<ID=TRIOCHECK_TP_TYPES,Number=.,Type=String,Description="SVTYPE values for true-positive child alleles">',

        '##INFO=<ID=TRIOCHECK_TP_CLASSES,Number=.,Type=String,Description="Variant classes for true-positive child alleles">',

        '##INFO=<ID=TRIOCHECK_TP_LINES,Number=1,Type=Integer,Description="Input VCF line number reported by Triocheck">',

    ]

    with open_text_write(outpath) as out:

        inserted_info_headers = False

        for header in header_lines:

            if header.startswith("#CHROM") and not inserted_info_headers:

                for h in tp_headers:

                    out.write(h + "\n")

                inserted_info_headers = True

            out.write(header + "\n")

        if not inserted_info_headers:

            for h in tp_headers:

                out.write(h + "\n")

        for line_no in sorted(tp_by_line):

            raw = raw_by_line.get(line_no)

            if raw is None:

                continue

            hits = tp_by_line[line_no]

            recs = [x[0] for x in hits]

            mother_recs = [rec for rec, mom_hit, _dad_hit in hits if mom_hit]

            father_recs = [rec for rec, _mom_hit, dad_hit in hits if dad_hit]

            annotations = [

                "TRIOCHECK_TP",

                "TRIOCHECK_TP_ALLELES=" + _unique_join(str(r.allele_index) for r in recs),

                "TRIOCHECK_TP_TYPES=" + _unique_join(r.svtype for r in recs),

                "TRIOCHECK_TP_CLASSES=" + _unique_join(r.variant_class for r in recs),

                f"TRIOCHECK_TP_LINES={line_no}",

            ]

            if mother_recs:

                annotations.append(

                    "TRIOCHECK_TP_MOTHER_ALLELES="

                    + _unique_join(str(r.allele_index) for r in mother_recs)

                )

            if father_recs:

                annotations.append(

                    "TRIOCHECK_TP_FATHER_ALLELES="

                    + _unique_join(str(r.allele_index) for r in father_recs)

                )

            out.write(add_info_annotations(raw, annotations) + "\n")

            tp_row_count += 1

    return tp_row_count, tp_allele_count


def summarize(

    vcf_path: str,

    child_selector: str,

    mother_selector: str,

    father_selector: str,

    child_records: List[VariantRecord],

    mother_records: List[VariantRecord],

    father_records: List[VariantRecord],

    matched_samples: Dict[str, List[str]],

    total_records: int,

    slop: int,

    min_ratio: float,

    min_size: int,

    include_all: bool,

    prefix_match: bool,

    exclude_child_varins: bool = True,

    original_distance: Optional[OriginalDistance] = None,

    parent_missing_stats: Optional[Dict[str, int]] = None,

    mom_hits: Optional[List[bool]] = None,

    dad_hits: Optional[List[bool]] = None,

) -> str:

    if mom_hits is None:

        mom_hits = find_parent_matches(child_records, mother_records, slop, min_ratio, original_distance)

    if dad_hits is None:

        dad_hits = find_parent_matches(child_records, father_records, slop, min_ratio, original_distance)

    if parent_missing_stats is None:

        parent_missing_stats = {

            "excluded_child_records": 0,

            "mother_missing": 0,

            "father_missing": 0,

            "both_missing": 0,

            "child_missing_original": 0,

            "mother_missing_original": 0,

            "father_missing_original": 0,

        }

    total = len(child_records)

    found_in_mother = sum(mom_hits)

    found_in_father = sum(dad_hits)

    found_in_either = sum(1 for m, d in zip(mom_hits, dad_hits) if m or d)

    found_in_both = sum(1 for m, d in zip(mom_hits, dad_hits) if m and d)

    child_only = total - found_in_either

    type_counter: Dict[str, List[int]] = {}

    for rec, m, d in zip(child_records, mom_hits, dad_hits):

        key = rec.variant_class if include_all else rec.svtype

        vals = type_counter.setdefault(key, [0, 0])

        vals[0] += 1

        if m or d:

            vals[1] += 1

    lines: List[str] = []

    lines.append(f"vcf\t{vcf_path}")

    lines.append(f"total_vcf_records_scanned\t{total_records}")

    lines.append(f"child_selector\t{child_selector}")

    lines.append(f"child_columns\t{','.join(matched_samples['child'])}")

    lines.append(f"mother_selector\t{mother_selector}")

    lines.append(f"mother_columns\t{','.join(matched_samples['mother'])}")

    lines.append(f"father_selector\t{father_selector}")

    lines.append(f"father_columns\t{','.join(matched_samples['father'])}")

    lines.append(f"sample_match_mode\t{'prefix' if prefix_match else 'exact'}")

    lines.append(f"min_size_bp\t{min_size}")

    if original_distance is None:

        lines.append("coordinate_match_mode\tmerged_vcf")

        lines.append(f"breakpoint_slop_bp\t{slop}")

        lines.append("breakpoint_slop_applied\ttrue")

        lines.append("original_distance\tNA")

    else:

        lines.append("coordinate_match_mode\toriginal_reference_from_hsv_cigar")

        lines.append(f"original_distance\t{original_distance}")

        lines.append("original_distance_definition\texact=reference_start_and_end_identical;0=reference_interval_overlap;N=reference_interval_gap_le_N_bp")

        lines.append("original_coordinate_source\tVCF_CHROM_and_POS_plus_per_sample_EXTENDGRAPHCIGAR_row_padding_and_main_reference_span")

        lines.append(f"breakpoint_slop_bp\t{slop}")

        lines.append("breakpoint_slop_applied\tfalse")

    lines.append(f"size_ratio_min\t{min_ratio}")

    lines.append(f"mode\t{'all' if include_all else 'sv_only'}")

    lines.append(f"exclude_child_varins\t{str(exclude_child_varins).lower()}")

    lines.append(f"include_child_varins\t{str(not exclude_child_varins).lower()}")

    lines.append("exclude_parent_missing\ttrue")

    lines.append("parent_missing_definition\tany selected mother/father sample column with '.', missing GT, or partial missing GT")

    lines.append(f"child_records_excluded_parent_missing\t{parent_missing_stats.get('excluded_child_records', 0)}")

    lines.append(f"child_records_excluded_mother_missing\t{parent_missing_stats.get('mother_missing', 0)}")

    lines.append(f"child_records_excluded_father_missing\t{parent_missing_stats.get('father_missing', 0)}")

    lines.append(f"child_records_excluded_both_parents_missing\t{parent_missing_stats.get('both_missing', 0)}")

    lines.append(f"child_records_excluded_missing_original_coordinates\t{parent_missing_stats.get('child_missing_original', 0)}")

    lines.append(f"mother_carried_records_missing_original_coordinates\t{parent_missing_stats.get('mother_missing_original', 0)}")

    lines.append(f"father_carried_records_missing_original_coordinates\t{parent_missing_stats.get('father_missing_original', 0)}")

    lines.append(f"child_records_with_original_coordinates\t{sum(1 for r in child_records if r.original_loci)}")

    lines.append(f"mother_records_with_original_coordinates\t{sum(1 for r in mother_records if r.original_loci)}")

    lines.append(f"father_records_with_original_coordinates\t{sum(1 for r in father_records if r.original_loci)}")

    lines.append(f"child_records_loaded\t{len(child_records)}")

    lines.append(f"mother_records_loaded\t{len(mother_records)}")

    lines.append(f"father_records_loaded\t{len(father_records)}")

    lines.append(f"{'child_total_variants' if include_all else 'child_total_sv'}\t{total}")

    lines.append(f"child_found_in_mother\t{found_in_mother}")

    lines.append(f"child_found_in_father\t{found_in_father}")

    lines.append(f"child_found_in_either_parent\t{found_in_either}")

    lines.append(f"child_found_in_both_parents\t{found_in_both}")

    lines.append(f"child_not_found_in_parents\t{child_only}")

    if total > 0:

        lines.append(f"child_found_in_either_parent_fraction\t{found_in_either / total:.6f}")

        lines.append(f"child_not_found_in_parents_fraction\t{child_only / total:.6f}")

    else:

        lines.append("child_found_in_either_parent_fraction\tNA")

        lines.append("child_not_found_in_parents_fraction\tNA")

    lines.append("")

    lines.append(f"{'type' if include_all else 'svtype'}\ttotal\tfound_in_either_parent\tfraction_found")

    sort_order = {"SNP": 0, "INDEL": 1, "SV": 2, "OTHER": 3}

    for key in sorted(type_counter, key=lambda x: (sort_order.get(x, 99), x)):

        total_t, found_t = type_counter[key]

        frac = found_t / total_t if total_t else 0.0

        lines.append(f"{key}\t{total_t}\t{found_t}\t{frac:.6f}")

    return "\n".join(lines)



def load_haplotype_vcf(

    vcf_path: str,

    min_size: int,

    include_all: bool,

    exclude_varins: bool = False,

    use_original_coordinates: bool = False,

    reference_edge_trim: int = 0,

) -> Tuple[

    List[VariantRecord],

    List[str],

    int,

    List[str],

    Dict[int, str],

    int,

    Dict[str, List[Tuple[int, int]]],

    int,

]:

    """Load one independently generated haplotype VCF.

    In separate-file mode, parental coverage is taken from the VCF header's

    ``##referenceCoverage=<...>`` records. Those intervals are declared by

    graphcigartoref to be 0-based half-open and are independent of variant-row

    presence/absence.

    ``reference_edge_trim`` removes that many bases from both ends of every raw

    referenceCoverage interval before intervals are merged.

    Returns carried alleles, merged reference-coverage intervals, and the raw

    number of referenceCoverage header records parsed from this VCF.

    """

    records: List[VariantRecord] = []

    header_lines: List[str] = []

    raw_by_line: Dict[int, str] = {}

    sample_names: List[str] = []

    total_records = 0

    missing_original = 0

    saw_header = False

    reference_coverage: Dict[str, List[Tuple[int, int]]] = defaultdict(list)

    reference_coverage_count = 0

    with open_text(vcf_path) as fh:

        for line_no, raw in enumerate(fh, 1):

            if raw.startswith('#'):

                header_line = raw.rstrip('\n')

                header_lines.append(header_line)

                m_cov = REFERENCE_COVERAGE_RE.match(header_line)

                if m_cov:

                    _sample, cov_chrom, cov_start_s, cov_end_s = m_cov.groups()

                    cov_start = int(cov_start_s)

                    cov_end = int(cov_end_s)

                    if cov_end > cov_start:

                        reference_coverage[cov_chrom].append((cov_start, cov_end))

                    reference_coverage_count += 1

            if raw.startswith('##'):

                continue

            if raw.startswith('#CHROM'):

                saw_header = True

                header_fields = raw.rstrip('\n').split('\t')

                if len(header_fields) < 10:

                    raise SystemExit(f'ERROR: {vcf_path} has no sample columns.')

                sample_names = header_fields[9:]

                continue

            if raw.startswith('#'):

                continue

            if not saw_header:

                raise SystemExit(f'ERROR: #CHROM header not found before records in {vcf_path}')

            total_records += 1

            fields = raw.rstrip('\n').split('\t')

            if len(fields) < 10:

                continue

            chrom, pos_s, vid, ref, alt, _qual, _flt, info_s = fields[:8]

            try:

                pos = int(pos_s)

            except ValueError:

                continue

            info = parse_info(info_s)

            fmt_keys = fields[8].split(':') if fields[8] else []

            sample_fields = fields[9:]

            for allele_index, alt_allele in enumerate(alt.split(','), 1):

                if not sample_group_has_allele(sample_fields, fmt_keys, allele_index):

                    continue

                if exclude_varins and sample_group_hsv_type_contains(

                    sample_fields, fmt_keys, allele_index, 'INS_'

                ):

                    continue

                allele_records, allele_missing_original = build_group_variant_records(

                    chrom=chrom,

                    pos=pos,

                    vid=vid,

                    ref=ref,

                    alt_allele=alt_allele,

                    info=info,

                    allele_index=allele_index,

                    source=vcf_path,

                    line_no=line_no,

                    sample_fields=sample_fields,

                    sample_names=sample_names,

                    fmt_keys=fmt_keys,

                    min_size=min_size,

                    include_all=include_all,

                    use_original_coordinates=use_original_coordinates,

                    include_query_metadata=True,

                )

                missing_original += allele_missing_original

                if not allele_records:

                    continue

                records.extend(allele_records)

                raw_by_line.setdefault(line_no, raw.rstrip('\n'))

    if not saw_header:

        raise SystemExit(f'ERROR: #CHROM header line was not found in {vcf_path}')

    reference_coverage = merge_reference_coverage_intervals(

        reference_coverage, edge_trim=reference_edge_trim

    )

    return (

        records,

        sample_names,

        total_records,

        header_lines,

        raw_by_line,

        missing_original,

        reference_coverage,

        reference_coverage_count,

    )


def exclude_child_parent_missing_separate(

    child_records: Sequence[VariantRecord],

    mother_paths: Sequence[str],

    father_paths: Sequence[str],

    coverage_by_source: Dict[str, Dict[str, List[Tuple[int, int]]]],

    stats_by_source: Optional[Dict[str, Dict[str, int]]] = None,

) -> Tuple[List[VariantRecord], Dict[str, int], Dict[str, int]]:

    """Exclude child alleles when a parental haplotype lacks reference coverage.

    Coverage comes only from ``##referenceCoverage`` header intervals. These are

    0-based half-open. For each child allele, the child reference anchor is the

    reconstructed original reference start when available; otherwise VCF POS-1.

    To mirror the original multi-sample rule where *any* selected parental

    haplotype missing makes that parent missing, mother_missing is True when any

    maternal haplotype VCF does not cover the child anchor; likewise for father.

    """

    kept: List[VariantRecord] = []

    stats = {

        'excluded_child_records': 0,

        'mother_missing': 0,

        'father_missing': 0,

        'both_missing': 0,

    }

    qc = {

        'mother_child_records_covered_all_haplotypes': 0,

        'father_child_records_covered_all_haplotypes': 0,

        'mother_child_records_not_covered_one_or_more_haplotypes': 0,

        'father_child_records_not_covered_one_or_more_haplotypes': 0,

    }

    for rec in child_records:

        anchors0 = child_reference_anchors0(rec)

        def path_covers(path: str) -> bool:

            cov = coverage_by_source.get(path, {})

            # A record can theoretically expose more than one original locus.

            # Require every child anchor represented by this record to be covered.

            return bool(anchors0) and all(

                reference_coverage_contains(cov, rec.chrom, anchor0)

                for anchor0 in anchors0

            )

        mother_covered = [path_covers(path) for path in mother_paths]

        father_covered = [path_covers(path) for path in father_paths]

        mother_missing = (not mother_covered) or (not all(mother_covered))

        father_missing = (not father_covered) or (not all(father_covered))

        if stats_by_source is not None:

            source_stats = stats_by_source.setdefault(rec.source, {})

            for key, excluded in (

                ('records_excluded_parent_missing', mother_missing or father_missing),

                ('records_excluded_mother_missing', mother_missing),

                ('records_excluded_father_missing', father_missing),

                ('records_excluded_both_parents_missing', mother_missing and father_missing),

            ):

                source_stats[key] = source_stats.get(key, 0) + int(excluded)

        if mother_missing:

            qc['mother_child_records_not_covered_one_or_more_haplotypes'] += 1

        else:

            qc['mother_child_records_covered_all_haplotypes'] += 1

        if father_missing:

            qc['father_child_records_not_covered_one_or_more_haplotypes'] += 1

        else:

            qc['father_child_records_covered_all_haplotypes'] += 1

        if mother_missing or father_missing:

            stats['excluded_child_records'] += 1

            if mother_missing:

                stats['mother_missing'] += 1

            if father_missing:

                stats['father_missing'] += 1

            if mother_missing and father_missing:

                stats['both_missing'] += 1

            continue

        kept.append(rec)

    return kept, stats, qc

def _output_path_for_source(basepath: str, source: str, multiple: bool) -> str:

    if not multiple:

        return basepath

    base = basepath

    if base.endswith('.vcf.gz'):

        root, ext = base[:-7], '.vcf.gz'

    elif base.endswith('.vcf'):

        root, ext = base[:-4], '.vcf'

    elif base.endswith('.gz'):

        root, ext = base[:-3], '.gz'

    else:

        root, ext = base, '.vcf'

    tag = os.path.basename(source)

    if tag.endswith('.vcf.gz'):

        tag = tag[:-7]

    elif tag.endswith('.vcf'):

        tag = tag[:-4]

    tag = re.sub(r'[^A-Za-z0-9_.-]+', '_', tag)

    return f'{root}.{tag}{ext}'


def write_split_child_outputs(

    basepath: str,

    mode: str,

    child_sources: Sequence[str],

    headers_by_source: Dict[str, List[str]],

    raw_by_source: Dict[str, Dict[int, str]],

    child_records: Sequence[VariantRecord],

    mom_hits: Sequence[bool],

    dad_hits: Sequence[bool],

) -> List[Tuple[str, int, int]]:

    """Write TP/FP VCFs separately for each child haplotype source."""

    out_stats: List[Tuple[str, int, int]] = []

    multiple = len(child_sources) > 1

    for source in child_sources:

        idxs = [i for i, rec in enumerate(child_records) if rec.source == source]

        recs = [child_records[i] for i in idxs]

        mh = [mom_hits[i] for i in idxs]

        dh = [dad_hits[i] for i in idxs]

        outpath = _output_path_for_source(basepath, source, multiple)

        if mode == 'fp':

            rows, alleles = write_fp_vcf(

                outpath, headers_by_source[source], raw_by_source[source], recs, mh, dh

            )

        else:

            rows, alleles = write_tp_vcf(

                outpath, headers_by_source[source], raw_by_source[source], recs, mh, dh

            )

        out_stats.append((outpath, rows, alleles))

    return out_stats


def summarize_child_haplotypes(

    child_paths: Sequence[str],

    child_records: Sequence[VariantRecord],

    mom_hits: Sequence[bool],

    dad_hits: Sequence[bool],

    include_all: bool,

    stats_by_source: Optional[Dict[str, Dict[str, int]]] = None,

) -> str:

    """Report retained calls by source without repeating parental matching."""

    counts = defaultdict(Counter)

    types = defaultdict(lambda: defaultdict(Counter))

    for rec, mother, father in zip(child_records, mom_hits, dad_hits):

        found = mother or father

        counts[rec.source].update({

            'total': 1,

            'found_in_mother': int(mother),

            'found_in_father': int(father),

            'found_in_either_parent': int(found),

            'found_in_both_parents': int(mother and father),

            'not_found_in_parents': int(not found),

        })

        kind = rec.variant_class if include_all else rec.svtype

        types[rec.source][kind].update(total=1, found=int(found))

    lines = []

    stats_by_source = stats_by_source or {}

    for number, path in enumerate(child_paths, 1):

        prefix = f'child_h{number}'

        values = counts[path]

        total = values['total']

        lines.append(f'{prefix}_vcf\t{path}')

        for key in (

            'vcf_records_scanned',

            'records_excluded_edge_trim',

            'records_excluded_reference_edge_trim',

            'records_excluded_missing_reference_coverage',

            'records_excluded_pa_edge',

            'records_missing_pa_coordinates',

            'records_excluded_missing_query_coordinates',

            'records_excluded_parent_missing',

            'records_excluded_mother_missing',

            'records_excluded_father_missing',

            'records_excluded_both_parents_missing',

            'records_excluded_missing_original_coordinates',

        ):

            lines.append(f'{prefix}_{key}\t{stats_by_source.get(path, {}).get(key, 0)}')

        lines.append(f'{prefix}_records_loaded\t{total}')

        total_key = 'total_variants' if include_all else 'total_sv'

        lines.append(f'{prefix}_{total_key}\t{total}')

        for key in ('found_in_mother', 'found_in_father', 'found_in_either_parent',

                    'found_in_both_parents', 'not_found_in_parents'):

            lines.append(f'{prefix}_{key}\t{values[key]}')

        for key in ('found_in_either_parent', 'not_found_in_parents'):

            fraction = f'{values[key] / total:.6f}' if total else 'NA'

            lines.append(f'{prefix}_{key}_fraction\t{fraction}')

        lines.append('')

    kind_label = 'type' if include_all else 'svtype'

    lines.append(f'child_haplotype\t{kind_label}\ttotal\tfound_in_either_parent\tfraction_found')

    sort_order = {'SNP': 0, 'INDEL': 1, 'SV': 2, 'OTHER': 3}

    for number, path in enumerate(child_paths, 1):

        for kind in sorted(types[path], key=lambda value: (sort_order.get(value, 99), value)):

            values = types[path][kind]

            total, found = values['total'], values['found']

            lines.append(f'h{number}\t{kind}\t{total}\t{found}\t{found / total:.6f}')

    return '\n'.join(lines)


def summarize_separate_files(

    child_paths: Sequence[str],

    mother_paths: Sequence[str],

    father_paths: Sequence[str],

    child_records: List[VariantRecord],

    mother_records: List[VariantRecord],

    father_records: List[VariantRecord],

    scanned_by_role: Dict[str, int],

    slop: int,

    min_ratio: float,

    min_size: int,

    include_all: bool,

    exclude_child_varins: bool,

    original_distance: Optional[OriginalDistance],

    missing_original_by_role: Dict[str, int],

    parent_missing_stats: Dict[str, int],

    parent_missing_lookup_qc: Dict[str, int],

    edge_trim_bp: int,

    edge_trim_excluded_by_role: Dict[str, int],

    reference_edge_excluded_child: int,

    child_missing_reference_coverage: int,

    pa_edge_distance_bp: int,

    pa_edge_excluded_child: int,

    pa_edge_missing_child: int,

    mom_hits: Sequence[bool],

    dad_hits: Sequence[bool],

    child_stats_by_source: Optional[Dict[str, Dict[str, int]]] = None,

    query_paths: Optional[str] = None,

) -> str:

    total = len(child_records)

    found_in_mother = sum(mom_hits)

    found_in_father = sum(dad_hits)

    found_in_either = sum(1 for m, d in zip(mom_hits, dad_hits) if m or d)

    found_in_both = sum(1 for m, d in zip(mom_hits, dad_hits) if m and d)

    child_only = total - found_in_either

    type_counter: Dict[str, List[int]] = {}

    for rec, m, d in zip(child_records, mom_hits, dad_hits):

        key = rec.variant_class if include_all else rec.svtype

        vals = type_counter.setdefault(key, [0, 0])

        vals[0] += 1

        if m or d:

            vals[1] += 1

    lines: List[str] = []

    lines.append(f'script_version\t{SCRIPT_VERSION}')

    lines.append('input_mode\tseparate_haplotype_vcfs')

    lines.append('child_vcfs\t' + ','.join(child_paths))

    lines.append('mother_vcfs\t' + ','.join(mother_paths))

    lines.append('father_vcfs\t' + ','.join(father_paths))

    lines.append('query_paths\t' + (query_paths or 'NA'))

    lines.append('query_contig_length_source\tfai' if query_paths else 'query_contig_length_source\tNA')

    lines.append('separate_file_parent_missing_semantics\treferenceCoverage_header_0_based_half_open_at_child_reference_anchor')

    lines.append(f'edge_trim_bp\t{edge_trim_bp}')

    lines.append('edge_trim_definition\ttrim_child_QUERYCOORD_and_every_raw_referenceCoverage_interval_before_merge')

    lines.append(f'child_records_excluded_edge_trim\t{edge_trim_excluded_by_role.get("child", 0)}')

    lines.append(

        f'child_records_excluded_query_edge_trim\t'

        f'{max(0, edge_trim_excluded_by_role.get("child", 0) - reference_edge_excluded_child)}'

    )

    lines.append(f'child_records_excluded_reference_edge_trim\t{reference_edge_excluded_child}')

    lines.append(f'child_records_excluded_missing_reference_coverage\t{child_missing_reference_coverage}')

    lines.append(f'parent_referenceCoverage_edge_trim_bp\t{edge_trim_bp}')

    lines.append(f'pa_edge_distance_bp\t{pa_edge_distance_bp}')

    lines.append('pa_edge_definition\tchild_only; reconstruct_each_PA_from_QUERYCOORD_plus_LABEL_H; exclude_record_if_any_PA_edge_distance_lt_threshold')

    lines.append(f'child_records_excluded_pa_edge\t{pa_edge_excluded_child}')

    lines.append(f'child_records_without_usable_pa_coordinates\t{pa_edge_missing_child}')

    lines.append(

        'child_records_excluded_missing_query_coordinates\t'

        + str(sum(

            values.get('records_excluded_missing_query_coordinates', 0)

            for values in (child_stats_by_source or {}).values()

        ))

    )

    lines.append(f'child_vcf_records_scanned\t{scanned_by_role.get("child", 0)}')

    lines.append(f'mother_vcf_records_scanned\t{scanned_by_role.get("mother", 0)}')

    lines.append(f'father_vcf_records_scanned\t{scanned_by_role.get("father", 0)}')

    lines.append(f'min_size_bp\t{min_size}')

    if original_distance is None:

        lines.append('coordinate_match_mode\tvcf_coordinates')

        lines.append(f'breakpoint_slop_bp\t{slop}')

        lines.append('breakpoint_slop_applied\ttrue')

        lines.append('original_distance\tNA')

    else:

        lines.append('coordinate_match_mode\toriginal_reference_from_hsv_cigar')

        lines.append(f'original_distance\t{original_distance}')

        lines.append('original_distance_definition\texact=reference_start_and_end_identical;0=reference_interval_overlap;N=reference_interval_gap_le_N_bp')

        lines.append('breakpoint_slop_applied\tfalse')

    lines.append(f'size_ratio_min\t{min_ratio}')

    lines.append(f'mode\t{"all" if include_all else "sv_only"}')

    lines.append(f'exclude_child_varins\t{str(exclude_child_varins).lower()}')

    lines.append('exclude_parent_missing\ttrue')

    lines.append('parent_missing_definition\tany parental haplotype whose edge-trimmed ##referenceCoverage intervals do not contain the child reference anchor')

    lines.append(f"child_records_excluded_parent_missing\t{parent_missing_stats.get('excluded_child_records', 0)}")

    lines.append(f"child_records_excluded_mother_missing\t{parent_missing_stats.get('mother_missing', 0)}")

    lines.append(f"child_records_excluded_father_missing\t{parent_missing_stats.get('father_missing', 0)}")

    lines.append(f"child_records_excluded_both_parents_missing\t{parent_missing_stats.get('both_missing', 0)}")

    lines.append(f"mother_child_records_covered_all_haplotypes\t{parent_missing_lookup_qc.get('mother_child_records_covered_all_haplotypes', 0)}")

    lines.append(f"father_child_records_covered_all_haplotypes\t{parent_missing_lookup_qc.get('father_child_records_covered_all_haplotypes', 0)}")

    lines.append(f"mother_child_records_not_covered_one_or_more_haplotypes\t{parent_missing_lookup_qc.get('mother_child_records_not_covered_one_or_more_haplotypes', 0)}")

    lines.append(f"father_child_records_not_covered_one_or_more_haplotypes\t{parent_missing_lookup_qc.get('father_child_records_not_covered_one_or_more_haplotypes', 0)}")

    lines.append(f'child_records_excluded_missing_original_coordinates\t{missing_original_by_role.get("child", 0)}')

    lines.append(f'mother_records_excluded_missing_original_coordinates\t{missing_original_by_role.get("mother", 0)}')

    lines.append(f'father_records_excluded_missing_original_coordinates\t{missing_original_by_role.get("father", 0)}')

    lines.append(f'child_records_loaded\t{len(child_records)}')

    lines.append(f'mother_records_loaded\t{len(mother_records)}')

    lines.append(f'father_records_loaded\t{len(father_records)}')

    lines.append(f'{"child_total_variants" if include_all else "child_total_sv"}\t{total}')

    lines.append(f'child_found_in_mother\t{found_in_mother}')

    lines.append(f'child_found_in_father\t{found_in_father}')

    lines.append(f'child_found_in_either_parent\t{found_in_either}')

    lines.append(f'child_found_in_both_parents\t{found_in_both}')

    lines.append(f'child_not_found_in_parents\t{child_only}')

    if total:

        lines.append(f'child_found_in_either_parent_fraction\t{found_in_either / total:.6f}')

        lines.append(f'child_not_found_in_parents_fraction\t{child_only / total:.6f}')

    else:

        lines.append('child_found_in_either_parent_fraction\tNA')

        lines.append('child_not_found_in_parents_fraction\tNA')

    lines.append('')

    lines.append(f'{"type" if include_all else "svtype"}\ttotal\tfound_in_either_parent\tfraction_found')

    sort_order = {'SNP': 0, 'INDEL': 1, 'SV': 2, 'OTHER': 3}

    for key in sorted(type_counter, key=lambda x: (sort_order.get(x, 99), x)):

        total_t, found_t = type_counter[key]

        frac = found_t / total_t if total_t else 0.0

        lines.append(f'{key}\t{total_t}\t{found_t}\t{frac:.6f}')

    if len(child_paths) > 1:

        lines.append('')

        lines.append(summarize_child_haplotypes(

            child_paths, child_records, mom_hits, dad_hits, include_all,

            child_stats_by_source,

        ))

    return '\n'.join(lines)

def main() -> int:

    parser = argparse.ArgumentParser(

        description=(

            "Check trio consistency in either one multi-sample VCF or independent "

            "child/mother/father haplotype VCFs with GT or HSV fields. "

            "Multiple child VCFs are reported both together and separately."

        )

    )

    parser.add_argument(

        "-i", "--input", default=None,

        help="Legacy mode: one multi-sample VCF, optionally .gz. Omit for separate-haplotype-VCF mode.",

    )

    parser.add_argument(

        "--child", nargs="+", required=True,

        help="Legacy mode: one child sample selector. Separate mode: one or more child haplotype VCF paths.",

    )

    parser.add_argument(

        "--mother", nargs="+", required=True,

        help="Legacy mode: one mother sample selector. Separate mode: one or more mother haplotype VCF paths.",

    )

    parser.add_argument(

        "--father", nargs="+", required=True,

        help="Legacy mode: one father sample selector. Separate mode: one or more father haplotype VCF paths.",

    )

    parser.add_argument("--min-size", type=int, default=DEFAULT_MIN_SIZE, help="Minimum SV size in bp. Default: 50")

    parser.add_argument("--slop", type=int, default=DEFAULT_BREAKPOINT_SLOP, help="Breakpoint tolerance in bp for ordinary VCF-coordinate SV matching. Default: 500")

    parser.add_argument(

        "--original-distance", default=None, metavar="BP|exact",

        help=(

            "Use reconstructed original per-haplotype REFERENCE coordinates for SV matching. "

            "0 requires interval overlap; a positive integer permits that many bp between intervals; "

            "'exact' requires identical start and end. When set, --slop is ignored for SV matching."

        ),

    )

    parser.add_argument(

        "--trim-edges", "--edge-trim", dest="edge_trim", nargs="?", const=10000,

        default=None, type=int, metavar="BP",

        help=(

            "Separate-haplotype mode: trim BP from both ends of child assembly "

            "query contigs and from every raw ##referenceCoverage interval before "

            "benchmarking. QUERYCOORD uses contig lengths from --query-paths; "

            "reference-side trimming applies to child and parental coverage. Default: "

            "10000 bp in separate-haplotype mode; "

            "use --trim-edges 0 to disable. In legacy multi-sample mode the default is 0."

        ),

    )

    parser.add_argument(

        "--pa-edge-distance", "--pa-edge-trim", dest="pa_edge_distance",

        default=None, type=int, metavar="BP",

        help=(

            "Separate-haplotype mode, CHILD ONLY: reconstruct every PA from "

            "QUERYCOORD + LABEL_H and exclude a child record when any PA is "

            "strictly closer than BP to either assembly-contig edge. Default: "

            "1000 bp. Use --pa-edge-distance 0 to disable. Parents and parental "

            "referenceCoverage are never filtered by this option."

        ),

    )

    parser.add_argument(

        "-q", "--query-paths", default=None, metavar="FILE",

        help=(

            "NAME FASTA [FAI] assembly list. Required in separate-haplotype "

            "mode when --trim-edges or --pa-edge-distance is greater than zero; "

            "child query-contig lengths are read from the listed .fai files."

        ),

    )

    parser.add_argument("--size-ratio", type=float, default=SIZE_RATIO_MIN, help="Minimum size ratio min(size)/max(size). Default: 0.7")

    parser.add_argument("--all", action="store_true", help="Check all variant classes instead of only SV > min-size")

    parser.add_argument("--prefix-match", action="store_true", help="Legacy single-VCF mode: allow selectors such as HG00514 to match HG00514_h1/HG00514_h2")

    parser.add_argument(

        "--VarINS", "--varins", dest="varins", action="store_true",

        help="Include child-carried HSV alleles whose child genotype type field contains INS_. Default: exclude them.",

    )

    parser.add_argument(

        "--fp", "--fp-vcf", dest="fp_vcf", nargs="?", const="fp.vcf", default=None,

        help=(

            "Write child-only variants. In separate-file mode with multiple child VCFs, "

            "one VCF is written per child haplotype by inserting the child filename before the extension."

        ),

    )

    parser.add_argument(

        "--tp", "--tp-vcf", dest="tp_vcf", nargs="?", const="tp.vcf", default=None,

        help=(

            "Write inherited/TP variants. In separate-file mode with multiple child VCFs, "

            "one VCF is written per child haplotype by inserting the child filename before the extension."

        ),

    )

    args = parser.parse_args()

    print(f"Triocheck_version\t{SCRIPT_VERSION}", file=sys.stderr)

    # Edge trimming is enabled by default for separate-haplotype benchmarking.

    # Legacy multi-sample VCF mode has no per-source ##referenceCoverage edges,

    # so preserve its historical behavior unless the user explicitly requests

    # an unsupported positive trim (which is rejected below).

    if args.edge_trim is None:

        args.edge_trim = 0 if args.input is not None else 10000

    if args.pa_edge_distance is None:

        args.pa_edge_distance = 0 if args.input is not None else DEFAULT_PA_EDGE_DISTANCE

    if args.min_size < 0:

        parser.error("--min-size must be >= 0")

    if args.slop < 0:

        parser.error("--slop must be >= 0")

    if args.edge_trim < 0:

        parser.error("--trim-edges/--edge-trim must be >= 0")

    if args.pa_edge_distance < 0:

        parser.error("--pa-edge-distance/--pa-edge-trim must be >= 0")

    if not (0 < args.size_ratio <= 1.0):

        parser.error("--size-ratio must be in (0, 1]")

    original_distance: Optional[OriginalDistance] = None

    if args.original_distance is not None:

        raw_original_distance = str(args.original_distance).strip().lower()

        if raw_original_distance == "exact":

            original_distance = "exact"

        else:

            try:

                original_distance = int(raw_original_distance)

            except ValueError:

                parser.error("--original-distance must be a non-negative integer or 'exact'")

            if original_distance < 0:

                parser.error("--original-distance must be >= 0 or 'exact'")

    # Legacy mode: one multi-sample VCF and one selector per role.

    if args.input is not None:

        if args.pa_edge_distance > 0:

            parser.error(

                "--pa-edge-distance is child-assembly metadata filtering and currently "

                "requires separate-haplotype VCF mode"

            )

        if args.edge_trim > 0:

            parser.error(

                "--trim-edges currently requires separate-haplotype VCF mode because "

                "the filter is defined from each VCF's ##referenceCoverage intervals"

            )

        if len(args.child) != 1 or len(args.mother) != 1 or len(args.father) != 1:

            parser.error("With --input, --child/--mother/--father must each contain exactly one sample selector")

        child_selector, mother_selector, father_selector = args.child[0], args.mother[0], args.father[0]

        (

            child_records, mother_records, father_records, matched_samples,

            total_records, header_lines, raw_by_line, parent_missing_stats,

        ) = load_trio_records(

            vcf_path=args.input,

            child_selector=child_selector,

            mother_selector=mother_selector,

            father_selector=father_selector,

            min_size=args.min_size,

            include_all=args.all,

            prefix_match=args.prefix_match,

            exclude_child_varins=not args.varins,

            use_original_coordinates=original_distance is not None,

        )

        if original_distance is not None and child_records and not any(r.original_loci for r in child_records):

            parser.error("--original-distance was requested, but no benchmarkable child original reference coordinates could be reconstructed from HSV CIGARs")

        mom_hits = find_parent_matches(child_records, mother_records, args.slop, args.size_ratio, original_distance)

        dad_hits = find_parent_matches(child_records, father_records, args.slop, args.size_ratio, original_distance)

        print(summarize(

            vcf_path=args.input,

            child_selector=child_selector,

            mother_selector=mother_selector,

            father_selector=father_selector,

            child_records=child_records,

            mother_records=mother_records,

            father_records=father_records,

            matched_samples=matched_samples,

            total_records=total_records,

            slop=args.slop,

            min_ratio=args.size_ratio,

            min_size=args.min_size,

            include_all=args.all,

            prefix_match=args.prefix_match,

            exclude_child_varins=not args.varins,

            original_distance=original_distance,

            parent_missing_stats=parent_missing_stats,

            mom_hits=mom_hits,

            dad_hits=dad_hits,

        ))

        if args.fp_vcf:

            fp_rows, fp_alleles = write_fp_vcf(args.fp_vcf, header_lines, raw_by_line, child_records, mom_hits, dad_hits)

            print(f"fp_vcf\t{args.fp_vcf}")

            print(f"fp_vcf_rows\t{fp_rows}")

            print(f"fp_vcf_child_alleles\t{fp_alleles}")

        if args.tp_vcf:

            tp_rows, tp_alleles = write_tp_vcf(args.tp_vcf, header_lines, raw_by_line, child_records, mom_hits, dad_hits)

            print(f"tp_vcf\t{args.tp_vcf}")

            print(f"tp_vcf_rows\t{tp_rows}")

            print(f"tp_vcf_child_alleles\t{tp_alleles}")

        return 0

    # Separate-file mode: each role is the union of independently generated haplotype VCFs.

    for role, paths in (("child", args.child), ("mother", args.mother), ("father", args.father)):

        for path in paths:

            if not os.path.isfile(path):

                parser.error(f"{role} VCF does not exist: {path}")

    if (args.edge_trim > 0 or args.pa_edge_distance > 0) and not args.query_paths:

        parser.error(

            "--query-paths is required for child query/PA edge filtering; "

            "supply query_paths.txt or disable both --trim-edges and --pa-edge-distance"

        )

    try:

        query_lengths_by_sample = (

            read_query_path_contig_lengths(args.query_paths)

            if (args.edge_trim > 0 or args.pa_edge_distance > 0) else {}

        )

    except (OSError, ValueError) as error:

        parser.error(str(error))

    role_records: Dict[str, List[VariantRecord]] = {"child": [], "mother": [], "father": []}

    scanned_by_role = {"child": 0, "mother": 0, "father": 0}

    missing_original_by_role = {"child": 0, "mother": 0, "father": 0}

    child_headers: Dict[str, List[str]] = {}

    child_raw: Dict[str, Dict[int, str]] = {}

    child_stats_by_source: Dict[str, Dict[str, int]] = {}

    coverage_by_source: Dict[str, Dict[str, List[Tuple[int, int]]]] = {}

    coverage_header_count_by_source: Dict[str, int] = {}

    query_bounds_by_source: Dict[str, Dict[str, Tuple[int, int]]] = {}

    for role, paths in (("child", args.child), ("mother", args.mother), ("father", args.father)):

        for path in paths:

            recs, sample_names, scanned, headers, raw_rows, missing_original, reference_coverage, reference_coverage_count = load_haplotype_vcf(

                vcf_path=path,

                min_size=args.min_size,

                include_all=args.all,

                exclude_varins=(role == "child" and not args.varins),

                use_original_coordinates=original_distance is not None,

                reference_edge_trim=args.edge_trim,

            )

            role_records[role].extend(recs)

            scanned_by_role[role] += scanned

            missing_original_by_role[role] += missing_original

            coverage_by_source[path] = reference_coverage

            coverage_header_count_by_source[path] = reference_coverage_count

            if role == "child":

                if args.edge_trim > 0 or args.pa_edge_distance > 0:

                    try:

                        query_bounds_by_source[path] = (

                            query_bounds_for_vcf_samples(

                                sample_names, query_lengths_by_sample, path

                            )

                        )

                    except ValueError as error:

                        parser.error(str(error))

                child_headers[path] = headers

                child_raw[path] = raw_rows

                child_stats_by_source[path] = {

                    'vcf_records_scanned': scanned,

                    'records_excluded_missing_original_coordinates': missing_original,

                    'records_excluded_missing_query_coordinates': 0,

                    'records_excluded_edge_trim': 0,

                    'records_excluded_reference_edge_trim': 0,

                    'records_excluded_missing_reference_coverage': 0,

                    'records_excluded_pa_edge': 0,

                    'records_missing_pa_coordinates': 0,

                }

            print(f"loaded_{role}_vcf\t{path}\tsamples={','.join(sample_names)}\tcarried_records={len(recs)}\treferenceCoverage_intervals={reference_coverage_count}", file=sys.stderr)

    if args.edge_trim > 0:

        missing_child_coverage_headers = [

            path for path in args.child

            if coverage_header_count_by_source.get(path, 0) == 0

        ]

        if missing_child_coverage_headers:

            parser.error(

                "--trim-edges reference-side filtering requires "

                "##referenceCoverage header intervals in every child VCF; "

                "missing in: " + ", ".join(missing_child_coverage_headers)

            )

    # --trim-edges is symmetric in separate-file mode: child calls must be in
    # both the assembly-query core and the child's trimmed referenceCoverage core.
    # Parental referenceCoverage was already trimmed while loading and therefore
    # also affects the later parent-callability test.
    edge_trim_excluded_by_role = {"child": 0}

    reference_edge_excluded_child = 0

    child_missing_query_coordinates = 0

    child_missing_reference_coverage = 0

    if args.edge_trim > 0:

        (

            role_records["child"],

            edge_trim_excluded_by_role["child"],

            child_missing_query_coordinates,

        ) = filter_child_records_to_query_core(

            role_records["child"], query_bounds_by_source, args.edge_trim,

            stats_by_source=child_stats_by_source,

        )

        if child_missing_query_coordinates:

            print(

                "WARNING: excluded "

                f"{child_missing_query_coordinates} child record(s) without "

                "usable ASSEMBLYCONTIG/QUERYCOORD for edge trimming",

                file=sys.stderr,

            )

    if args.edge_trim > 0:

        (

            role_records["child"],

            reference_edge_excluded_child,

            child_missing_reference_coverage,

        ) = filter_child_records_to_reference_core(

            role_records["child"],

            coverage_by_source,

            stats_by_source=child_stats_by_source,

        )

        edge_trim_excluded_by_role["child"] += reference_edge_excluded_child

        if child_missing_reference_coverage:

            print(

                "WARNING: excluded "

                f"{child_missing_reference_coverage} child record(s) without "

                "usable child ##referenceCoverage for reference-side edge trimming",

                file=sys.stderr,

            )

    # PA edge filtering remains CHILD ONLY. Mother/father records are not removed

    # by PA metadata; their referenceCoverage has already been trimmed by

    # --trim-edges during loading.

    pa_edge_excluded_child = 0

    pa_edge_missing_child = 0

    if args.pa_edge_distance > 0:

        (

            role_records["child"],

            pa_edge_excluded_child,

            pa_edge_missing_child,

        ) = filter_child_records_by_pa_edge(

            role_records["child"],

            query_bounds_by_source,

            args.pa_edge_distance,

            stats_by_source=child_stats_by_source,

        )

        if pa_edge_missing_child:

            print(

                "WARNING: retained "

                f"{pa_edge_missing_child} child record(s) without fully usable "

                "PA coordinates/LABEL_H for PA-edge filtering",

                file=sys.stderr,

            )

    child_records_all = role_records["child"]

    mother_records = role_records["mother"]

    father_records = role_records["father"]

    required_coverage_paths = list(args.mother) + list(args.father)

    missing_coverage_headers = [

        path for path in required_coverage_paths

        if coverage_header_count_by_source.get(path, 0) == 0

    ]

    if missing_coverage_headers:

        parser.error(

            "Separate-file parent-missing filtering requires "

            "##referenceCoverage header intervals in every parental VCF; "

            "missing in: " + ", ".join(missing_coverage_headers)

        )

    if original_distance is not None and not child_records_all:

        parser.error(

            "No benchmarkable child records were loaded before parental coverage filtering. "

            f"child_records_missing_original_coordinates={missing_original_by_role.get('child', 0)}. "

            "LABEL_H/ALLELENAME PA metadata is not used for original-reference "

            "coordinate reconstruction."

        )

    child_records, parent_missing_stats, parent_missing_lookup_qc = exclude_child_parent_missing_separate(

        child_records=child_records_all,

        mother_paths=args.mother,

        father_paths=args.father,

        coverage_by_source=coverage_by_source,

        stats_by_source=child_stats_by_source,

    )

    if original_distance is not None and child_records_all and not child_records:

        parser.error(

            "All child records with reconstructed original coordinates were excluded by parental "

            "##referenceCoverage filtering; this is not a LABEL_H reconstruction failure. "

            f"loaded_before_coverage={len(child_records_all)}, "

            f"excluded_parent_missing={parent_missing_stats.get('excluded_child_records', 0)}, "

            f"mother_missing={parent_missing_stats.get('mother_missing', 0)}, "

            f"father_missing={parent_missing_stats.get('father_missing', 0)}"

        )

    mom_hits = find_parent_matches(child_records, mother_records, args.slop, args.size_ratio, original_distance)

    dad_hits = find_parent_matches(child_records, father_records, args.slop, args.size_ratio, original_distance)

    print(summarize_separate_files(

        child_paths=args.child,

        mother_paths=args.mother,

        father_paths=args.father,

        child_records=child_records,

        mother_records=mother_records,

        father_records=father_records,

        scanned_by_role=scanned_by_role,

        slop=args.slop,

        min_ratio=args.size_ratio,

        min_size=args.min_size,

        include_all=args.all,

        exclude_child_varins=not args.varins,

        original_distance=original_distance,

        missing_original_by_role=missing_original_by_role,

        parent_missing_stats=parent_missing_stats,

        parent_missing_lookup_qc=parent_missing_lookup_qc,

        edge_trim_bp=args.edge_trim,

        edge_trim_excluded_by_role=edge_trim_excluded_by_role,

        reference_edge_excluded_child=reference_edge_excluded_child,

        child_missing_reference_coverage=child_missing_reference_coverage,

        pa_edge_distance_bp=args.pa_edge_distance,

        pa_edge_excluded_child=pa_edge_excluded_child,

        pa_edge_missing_child=pa_edge_missing_child,

        mom_hits=mom_hits,

        dad_hits=dad_hits,

        child_stats_by_source=child_stats_by_source,

        query_paths=args.query_paths,

    ))

    if args.fp_vcf:

        for outpath, rows, alleles in write_split_child_outputs(

            args.fp_vcf, "fp", args.child, child_headers, child_raw,

            child_records, mom_hits, dad_hits,

        ):

            print(f"fp_vcf\t{outpath}")

            print(f"fp_vcf_rows\t{rows}")

            print(f"fp_vcf_child_alleles\t{alleles}")

    if args.tp_vcf:

        for outpath, rows, alleles in write_split_child_outputs(

            args.tp_vcf, "tp", args.child, child_headers, child_raw,

            child_records, mom_hits, dad_hits,

        ):

            print(f"tp_vcf\t{outpath}")

            print(f"tp_vcf_rows\t{rows}")

            print(f"tp_vcf_child_alleles\t{alleles}")

    return 0


if __name__ == "__main__":

    raise SystemExit(main())
