#!/usr/bin/env python3
"""Match blocked partition-graph BED regions with reference 31-mers.

The two BED inputs are the BED6 outputs from ``block_partition_alignments.py``.
Rows are compared only when BED column 5 (the 1-based hotspot/graph index) is
the same.  BED column 4 is used as the graph name and is checked for agreement
between the samples.

For every graph, this script:

* extracts every BED region from its sample assembly;
* calls the compiled KmerMatch executable beside this script;
* compares query regions only against reference regions with KmerMatch;
* converts the rectangular report to the matrix consumed by RefMatch.py; and
* runs the pasted RefMatch selection/tagging logic in the worker process.

Query regions are never compared with other query regions.  Consequently, two
similar queries without a reference match remain separate novel alleles; one is
not selected as a representative for the other.

The reference sample supplies the reference alleles. Results from every graph are merged
into one deterministic, RefMatch-compatible eight-column table.  With
``--norm-folder``, the per-graph FASTA and upper-triangular ``_norm.txt.gz``
files are also retained and can be passed directly to RefMatch.py.

The cohort mode accepts a partition batch list, one partition-local
``_align.txt`` containing every sample, and the normalized query-path list.
It reproduces the already-committed blocked BED from the graph cache so every
blocked interval retains its source sample, then compares all non-reference
samples with the reference blocks in one KmerMatch call. Reference rows are
included in the resulting ``PARTITION_refmatch.tsv`` table.
"""

from __future__ import annotations

import argparse
import csv
import collections as cl
import dataclasses
import gc
import gzip
import logging
import math
import multiprocessing
import os
import re
import shlex
import struct
import subprocess
import tempfile
from concurrent.futures import ProcessPoolExecutor, as_completed
from pathlib import Path
from typing import Dict, Iterable, List, Optional, Sequence, Tuple

from align_partition_hotspots import graph_name_from_list_entry
from block_partition_alignments import (
    block_one_alignment,
    cached_cohort_assemblysmall_config,
    cached_cohort_block_paths,
    cached_cohort_merge_links,
)
from minsetref_core import IndexedFasta, mp_context, sanitize_id, wrap_fasta
from summarize_partition_hotspot_segments import (
    query_name_matches_haplotype,
    read_alignment_output,
)
from uniform_graph_blocks import load_graphfixbreaks


LOG = logging.getLogger("match_partition_blocks")
REFMATCH_INITIAL_DISTANCE = 0.3
REFMATCH_DISTANCE_MARGIN = 0.2
KMER_RESCUE_MIN_COMMON = 5000.0
KMER_RESCUE_MIN_RATIO = 10.0


def _float32(value: float) -> float:
    """Round a number with the same precision as NumPy's float32 scalar."""
    return struct.unpack("=f", struct.pack("=f", float(value)))[0]


@dataclasses.dataclass(frozen=True)
class BedRegion:
    sample_index: int
    input_index: int
    hotspot_index: int
    graph_name: str
    contig: str
    start: int
    end: int
    strand: str


@dataclasses.dataclass(frozen=True)
class GraphTask:
    hotspot_index: int
    graph_name: str
    regions: Tuple[BedRegion, ...]


@dataclasses.dataclass(frozen=True)
class Allele:
    header: str
    locus: str
    info: str
    is_reference: bool


@dataclasses.dataclass(frozen=True)
class GraphResult:
    hotspot_index: int
    graph_name: str
    lines: Tuple[str, ...]
    region_count: int
    fasta_path: str = ""
    norm_path: str = ""


@dataclasses.dataclass(frozen=True)
class CohortTask:
    partition: str
    directory: str
    alignment: str
    blocking: str
    cache: str
    output: str
    resume: bool


@dataclasses.dataclass(frozen=True)
class CohortResult:
    partition: str
    status: str
    alignment_rows: int = 0
    blocked_regions: int = 0
    output_rows: int = 0
    output: str = ""
    message: str = ""


_FASTAS: Tuple[Optional[IndexedFasta], Optional[IndexedFasta]] = (None, None)
_SAMPLE_NAMES: Tuple[str, str] = ("", "")
_MASK_LOWERCASE = True
_WHITELIST: Tuple[str, ...] = ()
_NORM_FOLDER = ""
_KMERMATCH = ""
_COHORT_SAMPLE_FASTAS: Dict[str, str] = {}
_COHORT_SAMPLE_ORDER: Dict[str, int] = {}
_COHORT_REFERENCE = "CHM13_h1"
_COHORT_QUERY_PATHS = ""
_COHORT_SELECTED_SAMPLES: Tuple[str, ...] = ()
_COHORT_SAMPLE_LIST = ""
_COHORT_GFIXBREAKS = None
_COHORT_MERGE_SMALL = 20_000


def sample_name_from_fasta(path: str) -> str:
    name = os.path.basename(path)
    name = re.sub(r"(?i)\.(?:fasta|fna|fa)$", "", name)
    return sanitize_id(name)


def parse_whitelist(value: str) -> List[str]:
    """Accept a comma-separated value or a file containing such values."""
    if not value:
        return []
    if os.path.exists(value):
        opener = gzip.open if value.endswith(".gz") else open
        with opener(value, "rt") as handle:
            text = handle.read()
    else:
        text = value
    text = text.replace("\r", ",").replace("\n", ",")
    return [part.strip() for part in text.split(",") if part.strip()]


def header_whitelist_keys(header: str) -> List[str]:
    parts = header.split("_")
    keys = [header]
    if len(parts) > 1:
        keys.append("_".join(parts[:-1]))
    if len(parts) >= 3:
        keys.append("_".join(parts[-3:-1]))
    return keys


def whitelist_rank(header: str, whitelist_items: Sequence[str]) -> int:
    if not whitelist_items:
        return 10**9
    keys = set(header_whitelist_keys(header))
    for rank, token in enumerate(whitelist_items):
        if token in keys:
            return rank
    for rank, token in enumerate(whitelist_items):
        if token in header:
            return rank
    return len(whitelist_items) + 10**6


def read_block_bed(path: str, sample_index: int) -> List[BedRegion]:
    rows: List[BedRegion] = []
    graph_names: Dict[int, str] = {}
    with open(path, "rt") as handle:
        for line_number, raw in enumerate(handle, 1):
            stripped = raw.strip()
            if not stripped or stripped.startswith("#"):
                continue
            fields = stripped.split()
            if fields[:3] == ["chrom", "start", "end"]:
                continue
            if len(fields) < 6:
                raise ValueError(
                    f"{path}:{line_number}: expected BED6 from "
                    "block_partition_alignments.py"
                )
            try:
                start = int(fields[1])
                end = int(fields[2])
                hotspot_index = int(fields[4])
            except ValueError as error:
                raise ValueError(
                    f"{path}:{line_number}: BED start, end, and score/graph "
                    "index must be integers"
                ) from error
            if start < 0 or end <= start:
                raise ValueError(
                    f"{path}:{line_number}: invalid interval {start}-{end}"
                )
            if hotspot_index < 1:
                raise ValueError(
                    f"{path}:{line_number}: graph index must be positive"
                )
            graph_name = fields[3]
            previous_name = graph_names.setdefault(hotspot_index, graph_name)
            if previous_name != graph_name:
                raise ValueError(
                    f"{path}:{line_number}: graph index {hotspot_index} has "
                    f"both {previous_name!r} and {graph_name!r}"
                )
            strand = fields[5]
            if strand not in {"+", "-"}:
                raise ValueError(
                    f"{path}:{line_number}: BED strand must be + or -"
                )
            rows.append(BedRegion(
                sample_index=sample_index,
                input_index=len(rows),
                hotspot_index=hotspot_index,
                graph_name=graph_name,
                contig=fields[0],
                start=start,
                end=end,
                strand=strand,
            ))
    return rows


def build_graph_tasks(
    bed_ref: str,
    bed_query: str,
) -> Tuple[List[GraphTask], Tuple[int, int]]:
    rows_by_sample = (
        read_block_bed(bed_ref, 0),
        read_block_bed(bed_query, 1),
    )
    grouped: Dict[int, List[BedRegion]] = cl.defaultdict(list)
    names: Dict[int, str] = {}
    for rows in rows_by_sample:
        for row in rows:
            previous = names.setdefault(row.hotspot_index, row.graph_name)
            if previous != row.graph_name:
                raise ValueError(
                    f"graph index {row.hotspot_index} is named {previous!r} in "
                    f"one BED and {row.graph_name!r} in the other"
                )
            grouped[row.hotspot_index].append(row)
    tasks = [
        GraphTask(index, names[index], tuple(sorted(
            grouped[index],
            key=lambda region: (region.sample_index, region.input_index),
        )))
        for index in sorted(grouped)
    ]
    return tasks, (len(rows_by_sample[0]), len(rows_by_sample[1]))


def kmermatch_input_sequence(sequence: str, mask_lowercase: bool) -> str:
    """Preserve FASTA case for KmerMatch's symmetric lowercase weighting."""
    # The argument remains part of this helper's interface for callers and
    # tests; KmerMatch itself applies either 0.1 or 1.0 weight according to -m.
    del mask_lowercase
    return sequence


def load_kmermatch_matrix(
    report_path: str,
    alleles: Sequence[Allele],
    sparse: bool = False,
):
    headers = [allele.header for allele in alleles]
    header_to_index = {header: index for index, header in enumerate(headers)}
    query_indexes = {
        index for index, allele in enumerate(alleles)
        if not allele.is_reference
    }
    reference_indexes = {
        index for index, allele in enumerate(alleles)
        if allele.is_reference
    }
    size = len(headers)
    # Cohort partitions can contain tens of thousands of blocked regions but
    # only a small number of reference regions. KmerMatch reports a rectangle,
    # so retaining an N x N Python matrix would be unnecessarily quadratic.
    matrix = (
        [cl.defaultdict(float) for _ in range(size)]
        if sparse else [[0.0] * size for _ in range(size)]
    )
    seen = set()
    query_variances: Dict[int, float] = {}
    reference_variances: Dict[int, float] = {}
    required = {
        "seqAname", "seqBname", "common_kmers",
        "kmers_only_in_A", "kmers_only_in_B",
        "variance_A", "variance_B",
    }
    with open(report_path, newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        if reader.fieldnames is None or not required.issubset(reader.fieldnames):
            raise ValueError(
                f"{report_path}: unrecognized KmerMatch header "
                f"{reader.fieldnames!r}"
            )
        for line_number, row in enumerate(reader, 2):
            left_header = row["seqAname"]
            right_header = row["seqBname"]
            if left_header not in header_to_index or right_header not in header_to_index:
                raise ValueError(
                    f"{report_path}:{line_number}: KmerMatch returned an "
                    "unknown sequence name"
                )
            left = header_to_index[left_header]
            right = header_to_index[right_header]
            if left not in query_indexes or right not in reference_indexes:
                raise ValueError(
                    f"{report_path}:{line_number}: expected a query-reference "
                    f"pair, found {left_header!r}, {right_header!r}"
                )
            pair = (left, right)
            if pair in seen:
                raise ValueError(
                    f"{report_path}:{line_number}: duplicate KmerMatch pair "
                    f"{left_header!r}, {right_header!r}"
                )
            seen.add(pair)
            try:
                common = float(row["common_kmers"])
                only_left = float(row["kmers_only_in_A"])
                only_right = float(row["kmers_only_in_B"])
                variance_left = float(row["variance_A"])
                variance_right = float(row["variance_B"])
            except (TypeError, ValueError) as error:
                raise ValueError(
                    f"{report_path}:{line_number}: non-numeric KmerMatch value"
                ) from error
            if not all(math.isfinite(value) for value in (
                common, only_left, only_right,
                variance_left, variance_right,
            )):
                raise ValueError(
                    f"{report_path}:{line_number}: non-finite KmerMatch value"
                )
            if min(
                common, only_left, only_right,
                variance_left, variance_right,
            ) < 0.0:
                raise ValueError(
                    f"{report_path}:{line_number}: negative KmerMatch value"
                )
            previous = query_variances.setdefault(left, variance_left)
            if not math.isclose(
                previous, variance_left, rel_tol=1e-12, abs_tol=1e-12,
            ):
                raise ValueError(
                    f"{report_path}:{line_number}: inconsistent variance for "
                    f"query {left_header!r}"
                )
            previous = reference_variances.setdefault(right, variance_right)
            if not math.isclose(
                previous, variance_right, rel_tol=1e-12, abs_tol=1e-12,
            ):
                raise ValueError(
                    f"{report_path}:{line_number}: inconsistent variance for "
                    f"reference {right_header!r}"
                )
            matrix[left][right] = common
            matrix[right][left] = common

    expected = len(query_indexes) * len(reference_indexes)
    if len(seen) != expected:
        raise ValueError(
            f"{report_path}: expected {expected} query-reference rows, "
            f"found {len(seen)}"
        )
    for index, variance in query_variances.items():
        matrix[index][index] = variance
    for index, variance in reference_variances.items():
        matrix[index][index] = variance
    return matrix


def empty_kmermatch_matrix(allele_count: int, sparse: bool = False):
    matrix = (
        [cl.defaultdict(float) for _ in range(allele_count)]
        if sparse else [[0.0] * allele_count for _ in range(allele_count)]
    )
    for index in range(allele_count):
        matrix[index][index] = 1.0
    return matrix


def run_kmermatch_files(
    alleles: Sequence[Allele], query_path: str, reference_path: str,
    report_path: str, sparse: bool = False,
):
    if not _KMERMATCH:
        raise RuntimeError("KmerMatch executable was not initialized")
    completed = subprocess.run(
        [
            _KMERMATCH,
            "-i", query_path,
            "-r", reference_path,
            "-o", report_path,
            "-m", "1" if _MASK_LOWERCASE else "0",
        ],
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
        text=True,
        check=False,
    )
    if completed.returncode != 0:
        detail = completed.stderr.strip() or completed.stdout.strip()
        raise RuntimeError(
            f"KmerMatch exited with status {completed.returncode}"
            + (f": {detail}" if detail else "")
        )
    if not os.path.isfile(report_path):
        raise RuntimeError("KmerMatch did not create its report")
    return load_kmermatch_matrix(report_path, alleles, sparse=sparse)


def run_kmermatch(
    task: GraphTask,
    alleles: Sequence[Allele],
    sequences: Sequence[str],
    sparse: bool = False,
):
    query_records = [
        (allele, sequence)
        for allele, sequence in zip(alleles, sequences)
        if not allele.is_reference
    ]
    reference_records = [
        (allele, sequence)
        for allele, sequence in zip(alleles, sequences)
        if allele.is_reference
    ]
    if not query_records or not reference_records:
        # Preserve valid self norms without inventing cross-sample matches.
        return empty_kmermatch_matrix(len(alleles), sparse)
    with tempfile.TemporaryDirectory(
        prefix=f"kmermatch_graph_{task.hotspot_index}_",
    ) as directory:
        query_path = os.path.join(directory, "query.fasta")
        reference_path = os.path.join(directory, "reference.fasta")
        report_path = os.path.join(directory, "kmermatch.tsv")
        with open(query_path, "wt") as output:
            for allele, sequence in query_records:
                prepared = kmermatch_input_sequence(sequence, _MASK_LOWERCASE)
                output.write(f">{allele.header}\n{wrap_fasta(prepared)}\n")
        with open(reference_path, "wt") as output:
            for allele, sequence in reference_records:
                prepared = kmermatch_input_sequence(sequence, _MASK_LOWERCASE)
                output.write(f">{allele.header}\n{wrap_fasta(prepared)}\n")
        return run_kmermatch_files(
            alleles, query_path, reference_path, report_path, sparse,
        )


def _is_kmermatch_sigkill(error: BaseException) -> bool:
    return "KmerMatch exited with status -9" in str(error)


def kmer_distance(
    matrix: Sequence[Sequence[float]], left: int, right: int,
) -> float:
    # RefMatch.py stores the norm matrix as np.float32.  Preserve that
    # arithmetic here so comparisons exactly on its 0.3/0.2 retention
    # boundary receive the same classification without requiring NumPy.
    left_variance = max(1.0, _float32(matrix[left][left]))
    right_variance = max(1.0, _float32(matrix[right][right]))
    variance_product = _float32(left_variance * right_variance)
    denominator = _float32(math.sqrt(variance_product))
    ratio = _float32(
        _float32(matrix[left][right]) / denominator
    )
    return _float32(1.0 - ratio)


def _haplotype(header: str) -> str:
    return header.rsplit("_", 1)[0] if "_" in header else header


def refmatch_lines(
    alleles: Sequence[Allele],
    matrix: Sequence[Sequence[float]],
    whitelist_items: Sequence[str] = (),
) -> List[str]:
    """Apply the pasted RefMatch.py matching and tagging behavior."""
    reference_indexes = [
        index for index, allele in enumerate(alleles) if allele.is_reference
    ]
    reference_set = set(reference_indexes)
    query_priority = [
        whitelist_rank(allele.header, whitelist_items) for allele in alleles
    ]

    pairs = [
        (query_priority[left], kmer_distance(matrix, left, right), left, right)
        for left in range(len(alleles))
        if left not in reference_set
        for right in reference_indexes
        if matrix[left][right] > 0
    ]
    pairs.sort(key=lambda value: (value[0], value[1], value[2], value[3]))

    matched_query: Dict[int, float] = {}
    matched_ref: Dict[Tuple[str, str], float] = {}
    all_matches: Dict[int, List[Tuple[float, int, int]]] = cl.defaultdict(list)
    for _, distance, left, right in pairs:
        haplotype = _haplotype(alleles[left].header)
        ref_header = alleles[right].header
        pair = (min(haplotype, ref_header), max(haplotype, ref_header))
        best_query = matched_query.get(left, REFMATCH_INITIAL_DISTANCE)
        best_ref = matched_ref.get(pair, REFMATCH_INITIAL_DISTANCE)
        if (
            distance >= best_query + REFMATCH_DISTANCE_MARGIN
            and distance >= best_ref + REFMATCH_DISTANCE_MARGIN
        ):
            continue
        duplicate = 0
        if pair in matched_ref:
            duplicate = 1
        elif left in all_matches and all_matches[left][0][2] == 1:
            first = all_matches[left][0]
            all_matches[left][0] = (first[0], first[1], 2)
            duplicate = 2
        if left not in matched_query:
            matched_query[left] = distance
        if pair not in matched_ref:
            matched_ref[pair] = distance
        all_matches[left].append((distance, right, duplicate))

    # A merged query can contain a strong reference block plus substantial
    # flanking sequence.  Its normalized distance can consequently fail the
    # ordinary RefMatch cutoff even when one reference has overwhelming raw
    # weighted-k-mer support.  Rescue only queries for which RefMatch retained
    # no hit at all, and only when the best common-kmer value is both strictly
    # above the absolute threshold and strictly more than the requested ratio
    # over the runner-up.  Keep the real normalized distance in the output.
    unmatched_queries = [
        index for index in range(len(alleles))
        if index not in all_matches and index not in reference_set
    ]
    rescued_queries = set()
    for left in unmatched_queries:
        ranked = sorted(
            (
                (matrix[left][right], right)
                for right in reference_indexes
            ),
            key=lambda value: (-value[0], value[1]),
        )
        if not ranked:
            continue
        best_common, right = ranked[0]
        second_common = ranked[1][0] if len(ranked) > 1 else 0.0
        if not (
            best_common > KMER_RESCUE_MIN_COMMON
            and best_common > KMER_RESCUE_MIN_RATIO * second_common
        ):
            continue

        distance = kmer_distance(matrix, left, right)
        haplotype = _haplotype(alleles[left].header)
        ref_header = alleles[right].header
        pair = (min(haplotype, ref_header), max(haplotype, ref_header))
        duplicate = 1 if pair in matched_ref else 0
        matched_query[left] = distance
        if pair not in matched_ref:
            matched_ref[pair] = distance
        all_matches[left].append((distance, right, duplicate))
        rescued_queries.add(left)

    novel_matches: Dict[int, List[Tuple[float, int, int]]] = cl.defaultdict(list)
    novel_tags: Dict[int, str] = {}
    novel_candidates = [
        index for index in range(len(alleles))
        if index not in all_matches and index not in reference_set
    ]
    novel_candidates.sort(key=lambda index: (query_priority[index], index))
    for index in novel_candidates:
        # An unmatched query is always its own novel allele.  In particular,
        # similarity to another unmatched query must not make the latter a
        # novel representative for this record.
        novel_matches[index].append((0.0, index, 0))
        novel_tags[index] = "Novel"

    tags = ("Pri", "Dup", "Cov")
    output = []
    for index, allele in enumerate(alleles):
        if index in all_matches:
            distance, match_index, duplicate = all_matches[index][0]
            matches = ";".join(
                f"{alleles[target].header}:{round(value, 6)}"
                for value, target, _ in all_matches[index]
                if matrix[index][target] > 0
            )
            tag = tags[duplicate]
            if index in rescued_queries:
                # Preserve the normal Pri/Dup/Cov prefix while explicitly
                # carrying the evidence type into GenomeLift. The normalized
                # distance remains unchanged and is never forged as <0.1.
                tag += "KmerRescue"
        elif index in novel_matches:
            distance, match_index, _ = novel_matches[index][0]
            matches = ";".join(
                f"{alleles[target].header}:{round(value, 6)}"
                for value, target, _ in novel_matches[index]
                # Keep the explicit self representative even when a short or
                # otherwise unusable novel allele has no usable 31-mers.
                if target == index or matrix[index][target] > 0
            )
            tag = novel_tags[index]
        else:
            distance = 0.0
            match_index = -1
            matches = ""
            tag = "Ref" if allele.is_reference else "Novel"

        match_header = alleles[match_index].header if match_index >= 0 else "NA"
        match_locus = alleles[match_index].locus if match_index >= 0 else "NA"
        output.append("\t".join((
            allele.header,
            match_header,
            allele.locus,
            match_locus,
            str(int("ENSE" in allele.info)),
            str(distance),
            tag,
            matches,
        )) + "\n")
    return output


def _atomic_text(path: str, text: str) -> None:
    temporary = path + f".tmp.{os.getpid()}"
    try:
        with open(temporary, "wt") as output:
            output.write(text)
        os.replace(temporary, path)
    finally:
        try:
            os.remove(temporary)
        except FileNotFoundError:
            pass


def _write_norm_artifacts(
    task: GraphTask,
    alleles: Sequence[Allele],
    sequences: Sequence[str],
    matrix: Sequence[Sequence[float]],
) -> Tuple[str, str]:
    if not _NORM_FOLDER:
        return "", ""
    prefix = f"{task.hotspot_index:06d}_{sanitize_id(task.graph_name)}"
    fasta_path = os.path.join(_NORM_FOLDER, prefix + ".fasta")
    norm_path = fasta_path + "_norm.txt.gz"
    fasta_text = "".join(
        f">{allele.header} {allele.locus} {allele.info}\n{wrap_fasta(sequence)}\n"
        for allele, sequence in zip(alleles, sequences)
    )
    _atomic_text(fasta_path, fasta_text)

    temporary = norm_path + f".tmp.{os.getpid()}"
    try:
        with gzip.open(temporary, "wt") as output:
            for index, row in enumerate(matrix):
                output.write(",".join(
                    format(value, ".17g") for value in row[index:]
                ) + "\n")
        os.replace(temporary, norm_path)
    finally:
        try:
            os.remove(temporary)
        except FileNotFoundError:
            pass
    return fasta_path, norm_path


def _close_worker_fastas() -> None:
    global _FASTAS
    for fasta in _FASTAS:
        if fasta is not None:
            fasta.close()
    _FASTAS = (None, None)


def _init_worker(
    fasta_ref: str,
    fasta_query: str,
    sample_names: Tuple[str, str],
    mask_lowercase: bool,
    whitelist: Tuple[str, ...],
    norm_folder: str,
    kmermatch: str,
) -> None:
    global _FASTAS, _SAMPLE_NAMES, _MASK_LOWERCASE
    global _WHITELIST, _NORM_FOLDER, _KMERMATCH
    gc.disable()
    _close_worker_fastas()
    _FASTAS = (IndexedFasta(fasta_ref), IndexedFasta(fasta_query))
    _SAMPLE_NAMES = sample_names
    _MASK_LOWERCASE = mask_lowercase
    _WHITELIST = whitelist
    _NORM_FOLDER = norm_folder
    _KMERMATCH = kmermatch


def _process_graph(task: GraphTask) -> GraphResult:
    sequences: List[str] = []
    alleles: List[Allele] = []
    graph_prefix = sanitize_id(task.graph_name)
    for region in task.regions:
        fasta = _FASTAS[region.sample_index]
        if fasta is None:
            raise RuntimeError("worker FASTA readers were not initialized")
        if region.contig not in fasta.index:
            raise ValueError(
                f"sample {_SAMPLE_NAMES[region.sample_index]} contig "
                f"{region.contig!r} is absent from its FASTA"
            )
        contig_length = fasta.length(region.contig)
        if region.end > contig_length:
            raise ValueError(
                f"sample {_SAMPLE_NAMES[region.sample_index]} region "
                f"{region.contig}:{region.start}-{region.end} exceeds FASTA "
                f"length {contig_length}"
            )
        sequence = fasta.fetch(region.contig, region.start, region.end)
        sample_name = _SAMPLE_NAMES[region.sample_index]
        header = f"{graph_prefix}_{sample_name}_{region.input_index + 1}"
        # GenomeLift derives sample/haplotype from the allele name and expects
        # this field to contain only CONTIG:START-ENDSTRAND. Prefixing the
        # contig with the sample would make HG38 canonical contigs look like
        # alternative contigs to GenomeLift's default filter.
        locus = f"{region.contig}:{region.start}-{region.end}{region.strand}"
        sequences.append(sequence)
        alleles.append(Allele(
            header=header,
            locus=locus,
            info=".",
            is_reference=region.sample_index == 0,
        ))

    matrix = run_kmermatch(task, alleles, sequences)
    lines = refmatch_lines(alleles, matrix, _WHITELIST)
    fasta_path, norm_path = _write_norm_artifacts(
        task, alleles, sequences, matrix,
    )
    return GraphResult(
        hotspot_index=task.hotspot_index,
        graph_name=task.graph_name,
        lines=tuple(lines),
        region_count=len(alleles),
        fasta_path=fasta_path,
        norm_path=norm_path,
    )


def write_merged_output(
    path: str,
    results: Iterable[GraphResult],
    include_header: bool,
) -> None:
    absolute = os.path.abspath(path)
    Path(os.path.dirname(absolute) or ".").mkdir(parents=True, exist_ok=True)
    temporary = absolute + f".tmp.{os.getpid()}"
    try:
        with open(temporary, "wt") as output:
            if include_header:
                output.write(
                    "query\tmatch\tquery_locus\tmatch_locus\t"
                    "ensembl\tdistance\ttag\tmatches\n"
                )
            for result in sorted(results, key=lambda value: value.hotspot_index):
                output.writelines(result.lines)
        os.replace(temporary, absolute)
    finally:
        try:
            os.remove(temporary)
        except FileNotFoundError:
            pass


def write_norm_manifest(folder: str, results: Sequence[GraphResult]) -> None:
    if not folder:
        return
    lines = ["hotspot_index\tgraph\tfasta\tnorm\tregions\n"]
    for result in sorted(results, key=lambda value: value.hotspot_index):
        lines.append("\t".join((
            str(result.hotspot_index),
            result.graph_name,
            result.fasta_path,
            result.norm_path,
            str(result.region_count),
        )) + "\n")
    _atomic_text(os.path.join(folder, "manifest.tsv"), "".join(lines))


def read_cohort_query_paths(path: str) -> Dict[str, str]:
    """Read NAME FASTA [FAI] while preserving normalized cohort order."""
    source = Path(path).expanduser().resolve()
    output: Dict[str, str] = {}
    with source.open("rt", encoding="utf-8") as handle:
        for line_number, raw in enumerate(handle, 1):
            if not raw.strip() or raw.lstrip().startswith("#"):
                continue
            fields = shlex.split(raw, comments=True)
            if len(fields) not in {2, 3}:
                raise ValueError(
                    f"{source}:{line_number}: expected NAME FASTA [FAI]"
                )
            sample = fields[0]
            fasta = Path(fields[1]).expanduser()
            if not fasta.is_absolute():
                fasta = source.parent / fasta
            fasta = fasta.resolve()
            if sample in output:
                raise ValueError(
                    f"{source}:{line_number}: duplicate sample {sample!r}"
                )
            if not fasta.is_file():
                raise FileNotFoundError(
                    f"{source}:{line_number}: FASTA not found: {fasta}"
                )
            output[sample] = str(fasta)
    if not output:
        raise ValueError(f"{source}: no cohort samples")
    return output


def cohort_sample_for_query(query_name: str, samples: Iterable[str]) -> str:
    matches = [
        sample for sample in samples
        if query_name_matches_haplotype(query_name, (sample,))
    ]
    if not matches:
        raise ValueError(
            f"cannot identify a query-path sample from {query_name!r}"
        )
    longest = max(map(len, matches))
    selected = [sample for sample in matches if len(sample) == longest]
    if len(selected) != 1:
        raise ValueError(
            f"ambiguous query-path samples for {query_name!r}: "
            + ",".join(sorted(selected))
        )
    return selected[0]


def read_cohort_partition_names(path: str) -> List[str]:
    names = []
    with open(path, "rt", encoding="utf-8") as handle:
        for line_number, raw in enumerate(handle, 1):
            if not raw.strip() or raw.lstrip().startswith("#"):
                continue
            token = raw.strip().split()[0]
            names.append(graph_name_from_list_entry(
                token, f"{path}:{line_number}",
            ))
    if not names:
        raise ValueError(f"partition list is empty: {path}")
    duplicate = next(
        (name for name, count in cl.Counter(names).items() if count > 1),
        None,
    )
    if duplicate is not None:
        raise ValueError(f"partition list contains duplicate {duplicate!r}")
    return names


def build_cohort_tasks(args: argparse.Namespace) -> List[CohortTask]:
    root = Path(args.graph_folder).expanduser().resolve()
    if not root.is_dir():
        raise NotADirectoryError(root)
    output_root = Path(
        args.output_graph_folder or args.graph_folder
    ).expanduser().resolve()
    output_root.mkdir(parents=True, exist_ok=True)
    tasks = []
    for partition in read_cohort_partition_names(args.partition_list):
        directory = root / partition
        output_directory = output_root / partition
        tasks.append(CohortTask(
            partition=partition,
            directory=str(directory),
            alignment=str(directory / f"{partition}_align.txt"),
            blocking=str(
                directory / f"{partition}_breakpoint_consistency.bed"
            ),
            cache=str(directory / f"{partition}cache.json"),
            output=str(output_directory / f"{partition}_refmatch.tsv"),
            resume=bool(args.resume),
        ))
    return tasks


def _initialize_cohort_worker(
    sample_fastas: Dict[str, str],
    reference: str,
    query_paths: str,
    graphfixbreaks: str,
    merge_small: int,
    kmermatch: str,
    masking: bool,
    whitelist: Tuple[str, ...],
    selected_samples: Tuple[str, ...],
    sample_list: str,
) -> None:
    global _COHORT_SAMPLE_FASTAS, _COHORT_SAMPLE_ORDER
    global _COHORT_REFERENCE, _COHORT_QUERY_PATHS, _COHORT_GFIXBREAKS
    global _COHORT_SELECTED_SAMPLES, _COHORT_SAMPLE_LIST
    global _COHORT_MERGE_SMALL, _KMERMATCH, _MASK_LOWERCASE, _WHITELIST
    gc.disable()
    _COHORT_SAMPLE_FASTAS = dict(sample_fastas)
    _COHORT_SAMPLE_ORDER = {
        sample: index for index, sample in enumerate(sample_fastas)
    }
    _COHORT_REFERENCE = reference
    _COHORT_QUERY_PATHS = os.path.abspath(query_paths)
    _COHORT_SELECTED_SAMPLES = tuple(selected_samples)
    _COHORT_SAMPLE_LIST = sample_list
    _COHORT_GFIXBREAKS = load_graphfixbreaks(graphfixbreaks)
    _COHORT_MERGE_SMALL = int(merge_small)
    _KMERMATCH = os.path.abspath(kmermatch)
    _MASK_LOWERCASE = bool(masking)
    _WHITELIST = whitelist


def _cohort_output_is_current(task: CohortTask) -> bool:
    if not task.resume:
        return False
    try:
        # Script timestamps do not invalidate a completed RefMatch table, but
        # its three semantic graph inputs must. Breakpoint consistency can be
        # regenerated without rebuilding the local graph; reusing an older
        # RefMatch table after that update would mix two blocking protocols.
        return (
            os.path.isfile(task.output)
            and os.path.getsize(task.output) > 0
            and os.stat(task.output).st_mtime_ns >= max(
                os.stat(path).st_mtime_ns
                for path in (task.alignment, task.blocking, task.cache)
            )
        )
    except OSError:
        return False


def _read_committed_blocking(path: str) -> List[str]:
    with open(path, "rt", encoding="utf-8") as handle:
        return [raw.rstrip("\r\n") for raw in handle if raw.strip()]


def _process_cohort_partition(task: CohortTask) -> CohortResult:
    missing = [
        path for path in (task.alignment, task.blocking, task.cache)
        if not os.path.isfile(path)
    ]
    if missing:
        return CohortResult(
            task.partition, "skipped_missing_input", message=missing[0],
        )
    if not os.path.getsize(task.alignment) or not os.path.getsize(task.cache):
        return CohortResult(
            task.partition, "skipped_empty_input",
            message=task.alignment if not os.path.getsize(task.alignment)
            else task.cache,
        )
    if _cohort_output_is_current(task):
        return CohortResult(
            task.partition, "reused", output=task.output,
        )
    rows = read_alignment_output(task.alignment)
    if any(row.hotspot_index != 1 for row in rows):
        bad = next(row.hotspot_index for row in rows if row.hotspot_index != 1)
        raise ValueError(
            f"{task.partition}: alignment hotspot index is {bad}, expected 1"
        )
    row_samples = [
        cohort_sample_for_query(row.query_name, _COHORT_SAMPLE_FASTAS)
        for row in rows
    ]
    if _COHORT_GFIXBREAKS is None:
        raise RuntimeError("cohort block-matching worker was not initialized")

    graphdb = _COHORT_GFIXBREAKS.graphDB.load_json(task.cache)
    control = cached_cohort_assemblysmall_config(graphdb)
    if control is None:
        raise ValueError(
            f"{task.partition}: cache has no frozen cohort assemblysmall control"
        )
    if control.get("scope") != "all_samples":
        raise ValueError(
            f"{task.partition}: cohort assemblysmall scope is "
            f"{control.get('scope')!r}, expected 'all_samples'"
        )
    if int(control.get("cutoff", -1)) != _COHORT_MERGE_SMALL:
        raise ValueError(
            f"{task.partition}: cohort assemblysmall cutoff is "
            f"{control.get('cutoff')!r}, expected {_COHORT_MERGE_SMALL}"
        )
    cohort_links = cached_cohort_merge_links(graphdb)
    cohort_paths = cached_cohort_block_paths(graphdb)
    if cohort_links is None or cohort_paths is None:
        raise ValueError(
            f"{task.partition}: incomplete cohort assemblysmall control"
        )
    blocked = []
    for input_index, row in enumerate(rows):
        row_blocks, _links = block_one_alignment(
            _COHORT_GFIXBREAKS, graphdb, task.partition, input_index, row,
            _COHORT_MERGE_SMALL,
            reference_merge_links=cohort_links,
            reference_block_paths=cohort_paths,
            fixed_merge_control=True,
        )
        blocked.extend(row_blocks)
    blocked.sort(key=lambda value: (value.input_index, value.region_index))
    generated = [row.line().rstrip("\n") for row in blocked]
    committed = _read_committed_blocking(task.blocking)
    if generated != committed:
        difference = next((
            index for index, (left, right) in enumerate(
                zip(generated, committed), 1,
            ) if left != right
        ), min(len(generated), len(committed)) + 1)
        generated_row = (
            generated[difference - 1]
            if difference <= len(generated) else "<missing>"
        )
        committed_row = (
            committed[difference - 1]
            if difference <= len(committed) else "<missing>"
        )
        raise ValueError(
            f"{task.partition}: reproduced blocking output differs from "
            f"{task.blocking} at row {difference} "
            f"({len(generated)} generated, {len(committed)} committed); "
            f"generated={generated_row!r}; committed={committed_row!r}"
        )

    ordered_blocks = sorted(blocked, key=lambda value: (
        row_samples[value.input_index] != _COHORT_REFERENCE,
        _COHORT_SAMPLE_ORDER[row_samples[value.input_index]],
        value.input_index, value.region_index,
    ))
    selected_samples = set(_COHORT_SELECTED_SAMPLES)
    ordered_blocks = [
        block for block in ordered_blocks
        if row_samples[block.input_index] in selected_samples
    ]
    sample_region_counts: Dict[str, int] = cl.defaultdict(int)
    alleles: List[Allele] = []
    block_indexes_by_sample: Dict[str, List[int]] = cl.defaultdict(list)
    for block in ordered_blocks:
        sample = row_samples[block.input_index]
        sample_region_counts[sample] += 1
        header = (
            f"{sanitize_id(task.partition)}_{sanitize_id(sample)}_"
            f"{sample_region_counts[sample]}"
        )
        locus = f"{block.contig}:{block.start}-{block.end}{block.strand}"
        block_indexes_by_sample[sample].append(len(alleles))
        alleles.append(Allele(
            header, locus, ".", sample == _COHORT_REFERENCE,
        ))

    has_reference = any(allele.is_reference for allele in alleles)
    has_query = any(not allele.is_reference for allele in alleles)
    if not has_reference or not has_query:
        matrix = empty_kmermatch_matrix(len(alleles), sparse=True)
    else:
        # Stream one sample at a time. A large cohort partition can represent
        # gigabases of blocked sequence; retaining all extracted strings in a
        # Python list would duplicate KmerMatch's own input memory.
        with tempfile.TemporaryDirectory(
            prefix=f"cohort_kmermatch_{sanitize_id(task.partition)}_",
        ) as temporary:
            query_path = os.path.join(temporary, "query.fasta")
            reference_path = os.path.join(temporary, "reference.fasta")
            report_path = os.path.join(temporary, "kmermatch.tsv")
            with open(query_path, "wt") as query_output, open(
                reference_path, "wt"
            ) as reference_output:
                for sample, allele_indexes in block_indexes_by_sample.items():
                    destination = (
                        reference_output
                        if sample == _COHORT_REFERENCE else query_output
                    )
                    reader = IndexedFasta(_COHORT_SAMPLE_FASTAS[sample])
                    try:
                        for allele_index in allele_indexes:
                            block = ordered_blocks[allele_index]
                            if block.contig not in reader.index:
                                raise KeyError(
                                    f"{task.partition}: {sample} contig "
                                    f"{block.contig!r} is absent from "
                                    f"{_COHORT_SAMPLE_FASTAS[sample]}"
                                )
                            if block.end > reader.length(block.contig):
                                raise ValueError(
                                    f"{task.partition}: {sample} interval "
                                    f"{block.contig}:{block.start}-{block.end} "
                                    f"exceeds contig length "
                                    f"{reader.length(block.contig)}"
                                )
                            sequence = reader.fetch(
                                block.contig, block.start, block.end,
                            )
                            destination.write(
                                f">{alleles[allele_index].header}\n"
                            )
                            destination.write(
                                kmermatch_input_sequence(
                                    sequence, _MASK_LOWERCASE,
                                ) + "\n"
                            )
                    finally:
                        reader.close()
            matrix = run_kmermatch_files(
                alleles, query_path, reference_path, report_path, sparse=True,
            )
    lines = refmatch_lines(alleles, matrix, _WHITELIST)
    write_merged_output(
        task.output,
        [GraphResult(1, task.partition, tuple(lines), len(alleles))],
        True,
    )
    return CohortResult(
        task.partition, "built", len(rows), len(blocked), len(lines),
        task.output,
    )


def write_cohort_manifest(
    path: str, results: Sequence[CohortResult],
) -> None:
    lines = [
        "partition\tstatus\talignment_rows\tblocked_regions\t"
        "output_rows\toutput\tmessage\n"
    ]
    for result in results:
        message = result.message.replace("\t", " ").replace("\n", " ")
        lines.append("\t".join(map(str, (
            result.partition, result.status, result.alignment_rows,
            result.blocked_regions, result.output_rows, result.output,
            message,
        ))) + "\n")
    absolute = os.path.abspath(path)
    Path(os.path.dirname(absolute) or ".").mkdir(parents=True, exist_ok=True)
    _atomic_text(absolute, "".join(lines))


def run_cohort(args: argparse.Namespace) -> int:
    sample_fastas = read_cohort_query_paths(args.query_paths)
    if args.reference_haplotype not in sample_fastas:
        raise ValueError(
            f"reference {args.reference_haplotype!r} is absent from "
            f"{args.query_paths}"
        )
    selected_samples = tuple(sample_fastas)
    sample_list = ""
    if args.sample_list:
        selected_samples = tuple(parse_whitelist(args.sample_list))
        if os.path.isfile(args.sample_list):
            sample_list = os.path.abspath(os.path.expanduser(args.sample_list))
        if not selected_samples:
            raise ValueError("--sample-list selects no samples")
        duplicate = next((
            sample for sample, count in cl.Counter(selected_samples).items()
            if count > 1
        ), None)
        if duplicate is not None:
            raise ValueError(
                f"--sample-list contains duplicate {duplicate!r}"
            )
        unknown = [
            sample for sample in selected_samples if sample not in sample_fastas
        ]
        if unknown:
            raise ValueError(
                "--sample-list contains sample(s) absent from --query-paths: "
                + ", ".join(unknown)
            )
        if args.reference_haplotype not in selected_samples:
            selected_samples = (args.reference_haplotype, *selected_samples)
    kmermatch = os.path.abspath(os.path.expanduser(args.kmermatch))
    if not os.path.isfile(kmermatch) or not os.access(kmermatch, os.X_OK):
        raise FileNotFoundError(
            f"KmerMatch executable is missing or not executable: {kmermatch}"
        )
    graphfixbreaks = os.path.abspath(os.path.expanduser(args.graphfixbreaks))
    if not os.path.isfile(graphfixbreaks):
        raise FileNotFoundError(graphfixbreaks)
    tasks = build_cohort_tasks(args)
    workers = min(args.jobs, len(tasks))
    initializer_args = (
        sample_fastas, args.reference_haplotype, args.query_paths,
        graphfixbreaks, args.merge_small, kmermatch, bool(args.masking),
        tuple(parse_whitelist(args.whitelist)), selected_samples, sample_list,
    )
    LOG.info(
        "Matching all-sample blocked regions for %d partitions with %d "
        "spawned worker process%s",
        len(tasks), workers, "" if workers == 1 else "es",
    )
    results: List[CohortResult] = []
    if workers == 1:
        _initialize_cohort_worker(*initializer_args)
        for task in tasks:
            try:
                results.append(_process_cohort_partition(task))
            except Exception as error:
                LOG.error("%s: %s", task.partition, error)
                results.append(CohortResult(
                    task.partition, "failed", message=str(error),
                ))
    else:
        context = multiprocessing.get_context("spawn")
        with ProcessPoolExecutor(
            max_workers=workers,
            mp_context=context,
            initializer=_initialize_cohort_worker,
            initargs=initializer_args,
        ) as executor:
            pending = {
                executor.submit(_process_cohort_partition, task): task
                for task in tasks[:workers * 2]
            }
            next_index = len(pending)
            while pending:
                for future in as_completed(tuple(pending)):
                    task = pending.pop(future)
                    try:
                        results.append(future.result())
                    except Exception as error:
                        LOG.error("%s: %s", task.partition, error)
                        results.append(CohortResult(
                            task.partition, "failed", message=str(error),
                        ))
                    if next_index < len(tasks):
                        next_task = tasks[next_index]
                        next_index += 1
                        pending[executor.submit(
                            _process_cohort_partition, next_task,
                        )] = next_task
                    break
                if len(results) % max(1, len(tasks) // 20) == 0:
                    LOG.info(
                        "Cohort block-match progress: %d/%d partitions",
                        len(results), len(tasks),
                    )
    order = {task.partition: index for index, task in enumerate(tasks)}
    results.sort(key=lambda result: order[result.partition])
    write_cohort_manifest(args.manifest, results)
    counts = cl.Counter(result.status for result in results)
    LOG.info(
        "Cohort block matching complete: %d built, %d reused, %d skipped, "
        "%d failed",
        counts["built"], counts["reused"],
        sum(count for status, count in counts.items() if status.startswith("skipped")),
        counts["failed"],
    )
    return 1 if args.strict and counts["failed"] else 0


def parse_args(argv: Optional[Sequence[str]] = None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--bed-ref", help="reference blocked BED6")
    parser.add_argument("--fasta-ref", help="reference assembly FASTA")
    parser.add_argument("--sample-ref", default="", help="reference sample label")
    parser.add_argument("--bed-query", help="query blocked BED6")
    parser.add_argument("--fasta-query", help="query assembly FASTA")
    parser.add_argument("--sample-query", default="", help="query sample label")
    parser.add_argument("-o", "--output")
    parser.add_argument(
        "--partition-list",
        help=(
            "cohort mode: one active partition or batch-list entry per line"
        ),
    )
    parser.add_argument(
        "-G", "--graph-folder",
        help="cohort mode: root containing partition graph folders",
    )
    parser.add_argument(
        "--output-graph-folder", default="",
        help=(
            "cohort mode: write PARTITION_refmatch.tsv under this separate "
            "partition root; default: --graph-folder"
        ),
    )
    parser.add_argument(
        "-q", "--query-paths",
        help="cohort mode: normalized NAME FASTA [FAI] assembly list",
    )
    parser.add_argument(
        "--manifest",
        help="cohort mode: batch status manifest written atomically",
    )
    parser.add_argument(
        "--reference-haplotype", default="CHM13_h1",
        help="cohort reference sample/haplotype (default: CHM13_h1)",
    )
    parser.add_argument(
        "--sample-list", default="",
        help=(
            "cohort mode: optional comma-separated or file-based subset of "
            "sample/haplotype labels to match; the reference is always added"
        ),
    )
    parser.add_argument(
        "--graphfixbreaks",
        default=os.path.join(
            os.path.dirname(os.path.abspath(__file__)), "gfixbreaks.py",
        ),
        help="cohort mode: bundled gfixbreaks.py implementation",
    )
    parser.add_argument(
        "--merge-small", type=int, default=20_000,
        help="cohort mode: gfixbreaks assemblysmall cutoff (default: 20000)",
    )
    parser.add_argument(
        "--resume", action="store_true",
        help="cohort mode: reuse current per-partition RefMatch tables",
    )
    parser.add_argument(
        "--strict", action="store_true",
        help="cohort mode: fail the batch when any partition fails",
    )
    parser.add_argument("-j", "--jobs", type=int, default=16)
    parser.add_argument(
        "--kmermatch",
        default=os.path.join(os.path.dirname(os.path.abspath(__file__)), "KmerMatch"),
        help="compiled KmerMatch executable (default: beside this script)",
    )
    parser.add_argument(
        "-m", "--masking", choices=(0, 1), type=int, default=1,
        help=(
            "give k-mers overlapping lowercase bases 0.1 variance weight; "
            "use 0 for full weight (default: 1)"
        ),
    )
    parser.add_argument(
        "-w", "--whitelist", default="",
        help="comma-separated priority sample/haplotype labels, or a file",
    )
    parser.add_argument(
        "--norm-folder", default="",
        help="optionally retain each graph's FASTA and RefMatch norm matrix",
    )
    parser.add_argument("--header", action="store_true")
    parser.add_argument("-v", "--verbose", action="store_true")
    args = parser.parse_args(argv)
    if args.jobs < 1:
        parser.error("--jobs must be positive")
    if args.merge_small < 0:
        parser.error("--merge-small must be nonnegative")
    if args.partition_list:
        missing = [
            option for option, value in (
                ("--graph-folder", args.graph_folder),
                ("--query-paths", args.query_paths),
                ("--manifest", args.manifest),
            ) if not value
        ]
        if missing:
            parser.error("cohort mode requires " + ", ".join(missing))
    else:
        missing = [
            option for option, value in (
                ("--bed-ref", args.bed_ref),
                ("--fasta-ref", args.fasta_ref),
                ("--bed-query", args.bed_query),
                ("--fasta-query", args.fasta_query),
                ("--output", args.output),
            ) if not value
        ]
        if missing:
            parser.error("two-sample mode requires " + ", ".join(missing))
    return args


def run(args: argparse.Namespace) -> int:
    fasta_paths = (
        os.path.abspath(os.path.expanduser(args.fasta_ref)),
        os.path.abspath(os.path.expanduser(args.fasta_query)),
    )
    for fasta in fasta_paths:
        if not os.path.isfile(fasta):
            raise FileNotFoundError(f"assembly FASTA not found: {fasta}")
    kmermatch = os.path.abspath(os.path.expanduser(args.kmermatch))
    if not os.path.isfile(kmermatch):
        raise FileNotFoundError(
            f"KmerMatch executable not found: {kmermatch}; build it with "
            f"`make -C {os.path.dirname(os.path.abspath(__file__))} KmerMatch`"
        )
    if not os.access(kmermatch, os.X_OK):
        raise PermissionError(f"KmerMatch is not executable: {kmermatch}")
    sample_names = (
        sanitize_id(args.sample_ref)
        if args.sample_ref else sample_name_from_fasta(fasta_paths[0]),
        sanitize_id(args.sample_query)
        if args.sample_query else sample_name_from_fasta(fasta_paths[1]),
    )
    if sample_names[0] == sample_names[1]:
        raise ValueError(
            "sample labels must differ; pass --sample-ref and --sample-query"
        )

    tasks, row_counts = build_graph_tasks(args.bed_ref, args.bed_query)
    norm_folder = ""
    if args.norm_folder:
        norm_folder = os.path.abspath(os.path.expanduser(args.norm_folder))
        Path(norm_folder).mkdir(parents=True, exist_ok=True)
    whitelist = tuple(parse_whitelist(args.whitelist))
    if not tasks:
        write_merged_output(args.output, (), args.header)
        write_norm_manifest(norm_folder, ())
        LOG.info("Both BED files are empty; wrote an empty merged output")
        return 0

    workers = min(args.jobs, len(tasks))
    initializer_args = (
        fasta_paths[0], fasta_paths[1], sample_names,
        bool(args.masking), whitelist, norm_folder, kmermatch,
    )
    LOG.info(
        "Matching %d reference + %d query BED regions across %d graphs with "
        "%d worker process%s",
        row_counts[0], row_counts[1], len(tasks), workers,
        "" if workers == 1 else "es",
    )
    results: List[GraphResult] = []
    if workers == 1:
        _init_worker(*initializer_args)
        try:
            for task in tasks:
                results.append(_process_graph(task))
        finally:
            _close_worker_fastas()
    else:
        gc.collect()
        sigkill_tasks: List[GraphTask] = []
        with ProcessPoolExecutor(
            max_workers=workers,
            mp_context=mp_context(),
            initializer=_init_worker,
            initargs=initializer_args,
        ) as executor:
            future_to_task = {
                executor.submit(_process_graph, task): task for task in tasks
            }
            report_every = max(1, len(tasks) // 20)
            for completed, future in enumerate(as_completed(future_to_task), 1):
                task = future_to_task[future]
                try:
                    results.append(future.result())
                except Exception as error:
                    if _is_kmermatch_sigkill(error):
                        sigkill_tasks.append(task)
                        LOG.warning(
                            "graph %d (%s): KmerMatch received SIGKILL; "
                            "deferring one low-concurrency retry",
                            task.hotspot_index, task.graph_name,
                        )
                    else:
                        raise RuntimeError(
                            f"graph {task.hotspot_index} "
                            f"({task.graph_name}): {error}"
                        ) from error
                if completed % report_every == 0 or completed == len(tasks):
                    LOG.info("K-mer matching progress: %d/%d graphs", completed, len(tasks))

        if sigkill_tasks:
            LOG.warning(
                "Retrying %d SIGKILL-terminated KmerMatch graph%s serially "
                "after all parallel workers have exited",
                len(sigkill_tasks), "" if len(sigkill_tasks) == 1 else "s",
            )
            _init_worker(*initializer_args)
            try:
                for retry_index, task in enumerate(
                    sorted(sigkill_tasks, key=lambda value: value.hotspot_index),
                    1,
                ):
                    try:
                        results.append(_process_graph(task))
                    except Exception as error:
                        raise RuntimeError(
                            f"graph {task.hotspot_index} ({task.graph_name}) "
                            f"failed again during serial SIGKILL retry: {error}"
                        ) from error
                    LOG.info(
                        "Serial KmerMatch retry progress: %d/%d graphs",
                        retry_index, len(sigkill_tasks),
                    )
            finally:
                _close_worker_fastas()

    write_merged_output(args.output, results, args.header)
    write_norm_manifest(norm_folder, results)
    LOG.info(
        "Wrote %d RefMatch rows from %d graphs to %s",
        sum(len(result.lines) for result in results), len(results), args.output,
    )
    return 0


def main(argv: Optional[Sequence[str]] = None) -> int:
    args = parse_args(argv)
    logging.basicConfig(
        level=logging.DEBUG if args.verbose else logging.INFO,
        format="[%(asctime)s] %(levelname)s: %(message)s",
    )
    try:
        return run_cohort(args) if args.partition_list else run(args)
    except Exception as error:
        LOG.error("%s", error)
        if args.verbose:
            LOG.exception("partition block matching failed")
        return 1


if __name__ == "__main__":
    raise SystemExit(main())
