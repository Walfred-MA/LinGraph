#!/usr/bin/env python3
"""Run graphcigartoref directly from GenomeLift and whole graph alignments.

This per-sample pipeline does not materialize graph-CIGAR, pair, record-FASTA,
or combined-graph intermediate files:

* GenomeLift query/reference pairs are resolved in memory;
* whole query/reference alignment graph CIGARs are sliced to allele intervals;
* multi-match query/reference sets are composed and aligned all-to-all once,
  retaining one complete contextual CIGAR plus per-block ownership metadata;
* the sample assembly is loaded into RAM once and query allele bases are sliced
  from its contigs;
* required graph path records are found in the same ``-G``/``-L`` graph layout
  used by ``align_partition_hotspots.py``.

The conversion stages are provided by graphcigartoref.py.  This per-sample
wrapper serializes their result in a seven-column layout that also carries the
complete query/reference whole-alignment intervals as downstream metadata.
"""

from __future__ import annotations

import argparse
import csv
import dataclasses
from contextlib import ExitStack, contextmanager
import gc
import json
import multiprocessing as mp
import os
import re
import sys
import time
import traceback
from typing import Dict, Iterable, Iterator, List, Mapping, Optional, Sequence, Set, Tuple

import graphcigartoref as core
from alternative_intervals import read_tag as read_alternative_tag
from graph_cigar_payloads import query_only_graph_cigar
from align_partition_hotspots import graph_name_from_list_entry
from block_partition_alignments import graph_prefix
from GenomeLift import allele_matrix_name, parse_allele_name
from local_graph_whole_cigar import WholeGraphCigar, complete_query_coverage
from local_reference_templates import LocalTemplate, read_templates
from minsetref_core import IndexedFasta, sanitize_id
from minsetref_segments import count_unmasked
from summarize_partition_hotspot_segments import read_alignment_output


GENOMELIFT_COLUMN_COUNT = 14
QUERY_CONSUMING_OPS = {"=", "M", "X", "I"}
GRAPH_SEGMENT_RE = re.compile(r"([<>])([^:<>]+):")
DEFAULT_REFERENCE_EXTENSION = 15_000
DEFAULT_EDGE_ALIGNMENT_SCORE = 100
DEFAULT_REFERENCE_FLANK_RESCUE = 10_000
MIN_REFERENCE_CANDIDATE_QUERY_COVERAGE = 0.5
HIGH_SIMILARITY_QUERY_COVERAGE = 0.9
MIN_LAST_RESORT_SIMILARITY_ALIGNED_BASES = 2_000


@dataclasses.dataclass(frozen=True)
class GenomeLiftRow:
    allele: str
    locus: core.Coord
    locus_text: str
    assigned_ref: str
    best_ref: str = ""
    class_type: str = ""
    assignment_source: str = ""
    assigned_refs: Tuple[str, ...] = ()
    merged_ref_locus: Optional[core.Coord] = None
    merged_ref_locus_text: str = ""
    part_index: int = 0
    group_token: str = ""
    left_extension: str = "0"
    right_extension: str = "0"


@dataclasses.dataclass(frozen=True)
class GraphEntry:
    hotspot_index: int
    graph_name: str
    graph_prefix: str
    fasta_path: str


@dataclasses.dataclass(frozen=True)
class PreparedInputRecord:
    allele: str
    graph_cigar: str
    hotspot_index: int
    # The interval actually represented by graph_cigar.  This can be wider
    # (or narrower) than GenomeLift column 3 after columns 13/14 are applied.
    locus: Optional[core.Coord] = None
    # Complete whole-alignment interval from --align-query/--align-ref.  This
    # is emitted separately from ``locus`` without changing the graph CIGAR.
    alignment_locus: Optional[core.Coord] = None
    # The ownership-trimmed but unextended block.  Stage1 aligns this region
    # before either flank, preventing repetitive extension paths from changing
    # the core placement.
    core_graph_cigar: Optional[str] = None
    core_locus: Optional[core.Coord] = None


class EmptySyntheticDeletionError(ValueError):
    """A synthetic deletion has no reference bases after safe clipping."""


@dataclasses.dataclass(frozen=True)
class SyntheticDeletionRow:
    """A reference-only run bracketed by two query-side anchors."""

    query_name: str
    query_coord: core.Coord
    ref_name: str
    ref_coord: core.Coord
    label: str = "DEL"

    def finished_key(self) -> Tuple[str, str, str]:
        # Per-sample-v2 column 5 is '.' for explicit DEL rows.
        return (self.query_name, self.ref_name, ".")

    def to_tsv(self, chrom_length: int) -> str:
        start = max(0, min(self.ref_coord.start, chrom_length))
        end = max(start, min(self.ref_coord.end, chrom_length))
        if end <= start:
            raise EmptySyntheticDeletionError(
                f"{self.ref_name}: empty synthetic deletion interval after "
                f"clipping to {self.ref_coord.chrom} length {chrom_length}"
            )
        operations = []
        if start:
            operations.append(f"{start}H")
        operations.append(f"{end - start}D")
        if end < chrom_length:
            operations.append(f"{chrom_length - end}H")
        cigar = f">{self.ref_coord.chrom}:" + "".join(operations)
        query_coord_text = coord_text(self.query_coord)
        ref_coord_text = coord_text(core.Coord(
            self.ref_coord.chrom, start, end, "+",
        ))
        return "\t".join((
            self.query_name,
            self.ref_name,
            query_coord_text,
            ref_coord_text,
            ".",
            ref_coord_text,
            cigar,
        ))


def render_synthetic_deletion_or_warn(
    row: SyntheticDeletionRow,
    chrom_length: int,
    source: str,
) -> Optional[str]:
    """Render a deletion, dropping a chromosome-clipped empty interval."""
    try:
        return row.to_tsv(chrom_length)
    except EmptySyntheticDeletionError as error:
        sys.stderr.write(
            f"[{source}] warning: {error}; dropping synthetic deletion "
            f"row {row.query_name!r}\n"
        )
        return None


@dataclasses.dataclass(frozen=True)
class PreparationTask:
    row: GenomeLiftRow
    is_reference: bool
    # Column-10 reference seeds are tied to one exact _align.txt row. Keeping
    # that row on the task prevents a second containment lookup from choosing
    # a different overlapping hotspot alignment. The seed selects/anchors the
    # candidate, while the complete row remains available for extension.
    alignment: Optional[object] = None


_PREP_QUERY_ALIGNMENT_INDEX = None
_PREP_REFERENCE_ALIGNMENT_INDEX = None
_PREP_MAX_EXTENSION = 10000
_PREP_CROSS_VALIDATION_EXTENSION = 0
_PREP_REFERENCE_EXTENSION = DEFAULT_REFERENCE_EXTENSION


def read_graph_entries(graph_list: str, graph_folder: str) -> List[GraphEntry]:
    graph_root = os.path.abspath(os.path.expanduser(graph_folder))
    if not os.path.isdir(graph_root):
        raise NotADirectoryError(f"graph folder not found: {graph_root}")
    entries: List[GraphEntry] = []
    with open(graph_list, "rt") as handle:
        for line_number, raw in enumerate(handle, 1):
            content = raw.strip()
            if not content or content.startswith("#"):
                continue
            graph_name = graph_name_from_list_entry(
                content.split()[0], f"{graph_list}:{line_number}",
            )
            entries.append(GraphEntry(
                hotspot_index=len(entries) + 1,
                graph_name=graph_name,
                graph_prefix=sanitize_id(graph_prefix(graph_name)),
                fasta_path=os.path.join(
                    graph_root, graph_name, graph_name + ".FA",
                ),
            ))
    if not entries:
        raise ValueError(f"graph list is empty: {graph_list}")
    counts: Dict[str, int] = {}
    for entry in entries:
        counts[entry.graph_prefix] = counts.get(entry.graph_prefix, 0) + 1
    duplicate = next((name for name, count in counts.items() if count > 1), None)
    if duplicate is not None:
        raise ValueError(
            f"graph list has multiple names that normalize to {duplicate!r}"
        )
    missing = [
        entry.fasta_path
        for entry in entries
        if not os.path.isfile(entry.fasta_path)
    ]
    if missing:
        examples = ", ".join(missing[:5])
        suffix = f" (and {len(missing) - 5} more)" if len(missing) > 5 else ""
        sys.stderr.write(
            "[graphcigartoref_persample] warning: "
            f"{len(missing)} graph path(s) selected by -L do not exist under "
            f"-G; affected candidate groups will be skipped: "
            f"{examples}{suffix}\n"
        )
    return entries


def read_genomelift(path: str) -> List[GenomeLiftRow]:
    rows: List[GenomeLiftRow] = []
    with open(path, "rt") as handle:
        for line_number, raw in enumerate(handle, 1):
            if not raw.strip() or raw.startswith("#"):
                continue
            fields = raw.rstrip("\r\n").split("\t")
            if fields[0] == "allelename":
                continue
            if len(fields) < GENOMELIFT_COLUMN_COUNT:
                raise ValueError(
                    f"{path}:{line_number}: expected at least "
                    f"{GENOMELIFT_COLUMN_COUNT} GenomeLift columns, found "
                    f"{len(fields)}"
                )
            # GenomeLift appends reference-status records with no assembly
            # locus. They are consumed separately after normal rows are read.
            if fields[0] in {"DEL", "NA"}:
                continue
            part_match = re.fullmatch(r"part_(\d+)", fields[6].strip())
            part_index = int(part_match.group(1)) if part_match else 0

            def part_value(value: str) -> str:
                values = value.split(";")
                if part_index and len(values) >= part_index:
                    return values[part_index - 1].strip()
                return value.strip()

            locus_text = part_value(fields[2])
            locus = core.parse_coord(locus_text)
            if locus is None:
                raise ValueError(
                    f"{path}:{line_number}: malformed assembly locus "
                    f"{locus_text!r}"
                )
            location_ref = fields[8].strip()
            # Grouped ``part_N`` rows can carry one semicolon-delimited value
            # per part in the display/reference columns.  Resolve column 2 at
            # the same time as the assembly locus so a similarity fallback is
            # tied to this part rather than to the whole display list.
            best_ref = part_value(fields[1])
            # GenomeLift can append a gene-state display value in brackets to
            # column 2.  The graph alignment records use the undecorated
            # reference allele name.
            if best_ref.endswith("]") and "[" in best_ref:
                best_ref = best_ref.split("[", 1)[0]
            if best_ref in {"", ".", "NA", "na", "None", "none"}:
                best_ref = ""
            class_type = fields[6].strip()
            # Retain column 9 for diagnostics. Column 10 remains the primary
            # graph-CIGAR seed; sanitized column 2 can be used later only as a
            # lower-priority same-graph fallback.
            if location_ref:
                assigned_ref = location_ref
                assignment_source = "location"
            else:
                assigned_ref = ""
                assignment_source = ""
            assigned_refs = tuple(
                value.strip() for value in location_ref.split(";")
                if value.strip()
            )
            merged_ref_locus_text = fields[9].strip()
            merged_ref_locus = core.parse_coord(merged_ref_locus_text)
            if merged_ref_locus is None and merged_ref_locus_text:
                # Direct GenomeLift placements historically omit the strand
                # in column 10. Coordinates are forward-reference intervals;
                # use '+' solely as their traversal default.
                merged_ref_locus = core.parse_coord(
                    merged_ref_locus_text + "+"
                )
            # GenomeLift is globally sorted by assembly position. Rows from a
            # different overlapping graph can therefore occur between part_1
            # and later parts of the same group. Identify a group from its
            # biological assignment instead of relying on line adjacency.
            group_token = ""
            if part_index:
                group_token = "\x1f".join((
                    allele_matrix_name(fields[0]),
                    locus.chrom,
                    merged_ref_locus_text,
                    fields[2].strip() if ";" in fields[2] else "",
                ))
            rows.append(GenomeLiftRow(
                allele=fields[0],
                locus=locus,
                locus_text=locus_text,
                assigned_ref=assigned_ref,
                best_ref=best_ref,
                class_type=class_type,
                assignment_source=assignment_source,
                assigned_refs=assigned_refs,
                merged_ref_locus=merged_ref_locus,
                merged_ref_locus_text=merged_ref_locus_text,
                part_index=part_index,
                group_token=group_token,
                left_extension=part_value(fields[12]),
                right_extension=part_value(fields[13]),
            ))
    return rows


def select_effective_query_rows(
    lift_rows: Sequence[GenomeLiftRow],
    query_genome: str,
    refhaplo: str,
) -> Tuple[List[GenomeLiftRow], int, Dict[str, int]]:
    """Select rows participating in the three-tier reference cascade.

    Priority 0 is a high-coverage same-graph column-2 match, priority 1 is one
    unmerged reference allele in column 9, priority 2 retries the column-2
    result at the ordinary threshold, and priority 3 is the coordinate in
    column 10. Candidate construction later verifies the named reference rows
    and their actual ``_align.txt`` intervals.
    """
    selected: List[GenomeLiftRow] = []
    query_row_count = 0
    assignment_counts = {
        "column9_single": 0,
        "bestref_similarity": 0,
        "column10": 0,
    }

    for row in lift_rows:
        try:
            genome = parse_allele_name(row.allele)[4]
        except ValueError:
            continue
        if genome != query_genome:
            continue
        if row.class_type.lower() == "del":
            continue
        query_row_count += 1
        matrix = allele_matrix_name(row.allele)

        def valid_named_reference(value: str) -> bool:
            value = value.strip()
            if not value or allele_matrix_name(value) != matrix:
                return False
            try:
                return parse_allele_name(value)[4] == refhaplo
            except ValueError:
                return False

        single_column9 = (
            len(row.assigned_refs) == 1
            and valid_named_reference(row.assigned_refs[0])
        )
        has_bestref = any(
            valid_named_reference(value)
            for value in row.best_ref.split(";")
        )
        if has_bestref:
            assignment_source = "bestref_similarity"
        elif single_column9:
            assignment_source = "column9_single"
        elif row.merged_ref_locus is not None:
            assignment_source = "column10"
        else:
            continue

        group_token = row.group_token
        if row.part_index:
            grouped_query_loci = (
                row.group_token.rsplit("\x1f", 1)[-1]
                if row.group_token else ""
            )
            group_token = "\x1f".join((
                allele_matrix_name(row.allele),
                row.locus.chrom,
                row.merged_ref_locus_text,
                grouped_query_loci,
            ))
        selected_row = dataclasses.replace(
            row,
            assigned_ref="",
            assigned_refs=(),
            assignment_source=assignment_source,
            group_token=group_token,
        )
        selected.append(selected_row)
        assignment_counts[assignment_source] += 1
    return selected, query_row_count, assignment_counts


def read_explicit_deletion_rows(
    path: str,
    lift_rows: Sequence[GenomeLiftRow],
    query_genome: str,
) -> List[SyntheticDeletionRow]:
    """Read explicit ``DEL REF . REFLOC ...`` GenomeLift status rows."""
    rows_by_name = {row.allele: row for row in lift_rows}
    deletions = []
    with open(path, "rt") as handle:
        for line_number, raw in enumerate(handle, 1):
            if not raw.strip() or raw.startswith("#"):
                continue
            fields = raw.rstrip("\r\n").split("\t")
            if not fields or fields[0] != "DEL":
                continue
            if len(fields) < GENOMELIFT_COLUMN_COUNT:
                raise ValueError(
                    f"{path}:{line_number}: DEL row has {len(fields)} columns; "
                    f"expected at least {GENOMELIFT_COLUMN_COUNT}"
                )
            row_genome = fields[4].strip()
            if row_genome and row_genome != query_genome:
                continue
            reference_names = tuple(
                value.strip() for value in fields[8].split(";")
                if value.strip()
            )
            ref_coord = core.parse_coord(fields[9].strip())
            left_anchor = rows_by_name.get(fields[10].strip())
            right_anchor = rows_by_name.get(fields[11].strip())
            if (
                not reference_names
                or ref_coord is None
                or left_anchor is None
                or right_anchor is None
                or left_anchor.locus.chrom != right_anchor.locus.chrom
            ):
                raise ValueError(
                    f"{path}:{line_number}: DEL row lacks its reference "
                    "interval or two same-contig query anchors"
                )
            if (
                left_anchor.locus.start + left_anchor.locus.end
                <= right_anchor.locus.start + right_anchor.locus.end
            ):
                query_strand = "+"
                left_facing = left_anchor.locus.end
                right_facing = right_anchor.locus.start
            else:
                query_strand = "-"
                left_facing = left_anchor.locus.start
                right_facing = right_anchor.locus.end
            breakpoint = max(0, (left_facing + right_facing) // 2)
            safe_chrom = sanitize_id(ref_coord.chrom)
            query_name = (
                f"DEL{safe_chrom}p{ref_coord.start}p{ref_coord.end}_"
                f"{query_genome}_{len(deletions) + 1}"
            )
            deletions.append(SyntheticDeletionRow(
                query_name=query_name,
                query_coord=core.Coord(
                    left_anchor.locus.chrom,
                    breakpoint,
                    breakpoint,
                    query_strand,
                ),
                ref_name=";".join(reference_names),
                ref_coord=core.Coord(
                    ref_coord.chrom, ref_coord.start, ref_coord.end, "+",
                ),
            ))
    return deletions


def infer_query_genome(query_alignment: str) -> str:
    name = os.path.basename(query_alignment)
    for suffix in (".align.txt", ".align.tsv", ".txt", ".tsv"):
        if name.endswith(suffix):
            name = name[:-len(suffix)]
            break
    if not name:
        raise ValueError("could not infer query genome from alignment filename")
    return name


def infer_reference_alignment(query_alignment: str, refhaplo: str) -> str:
    query_dir = os.path.dirname(os.path.abspath(query_alignment))
    alignment_root = os.path.dirname(query_dir)
    return os.path.join(alignment_root, refhaplo, f"{refhaplo}.align.txt")


def graph_cigar_query_span(graph_cigar: str) -> int:
    return sum(
        operation.n
        for segment in core.parse_graphic_segments(graph_cigar, "span")
        for operation in segment.ops
        if operation.op in QUERY_CONSUMING_OPS
    )


def slice_graph_cigar_by_query(
    graph_cigar: str, qstart: int, qend: int,
) -> str:
    qstart = int(qstart)
    qend = int(qend)
    if qstart < 0 or qend <= qstart:
        raise ValueError(f"invalid graph-CIGAR query slice {qstart}-{qend}")
    pieces: List[str] = []
    cursor = 0
    for segment in core.parse_graphic_segments(graph_cigar, "slice"):
        segment_span = sum(
            operation.n
            for operation in segment.ops
            if operation.op in QUERY_CONSUMING_OPS
        )
        overlap_start = max(qstart, cursor)
        overlap_end = min(qend, cursor + segment_span)
        if overlap_end > overlap_start:
            piece = core.slice_query_segment_by_query(
                segment,
                overlap_start - cursor,
                overlap_end - cursor,
            )
            if piece:
                pieces.append(piece)
        cursor += segment_span
    if qend > cursor:
        raise ValueError(
            f"graph-CIGAR query slice {qstart}-{qend} exceeds query span {cursor}"
        )
    if not pieces:
        raise ValueError(
            f"graph-CIGAR query slice {qstart}-{qend} contains no usable path"
        )
    sliced = "".join(pieces)
    observed = graph_cigar_query_span(sliced)
    expected = qend - qstart
    if observed != expected:
        raise ValueError(
            f"sliced graph-CIGAR query span is {observed}, expected {expected}"
        )
    return sliced


def build_alignment_index(
    alignments,
    graph_entries: Sequence[GraphEntry],
) -> Dict[Tuple[str, str], Tuple[object, ...]]:
    """Index alignments by graph prefix, including fresh unprefixed rows.

    Distributed/legacy alignments may already encode a stable graph prefix in
    ``query_name``; retain that identity even if the current list was reordered.
    Fresh ``align_partition_hotspots.py`` rows are named only
    ``SAMPLE_HAPLOTYPE_INDEX``.  For those rows, recover the missing prefix from
    alignment column 1 (``hotspot_index``) and the same ordered ``-L`` list.
    """
    known_prefixes = {entry.graph_prefix for entry in graph_entries}
    prefix_by_hotspot = {
        entry.hotspot_index: entry.graph_prefix for entry in graph_entries
    }
    grouped: Dict[Tuple[str, str], List[object]] = {}
    for row in alignments:
        query_name = str(getattr(row, "query_name", "") or "").strip()
        if not query_name:
            raise ValueError("alignment row has no query name")
        encoded_matrix = allele_matrix_name(query_name)
        if encoded_matrix in known_prefixes:
            matrix = encoded_matrix
        else:
            hotspot_index = int(getattr(row, "hotspot_index", 0) or 0)
            matrix = prefix_by_hotspot.get(hotspot_index, "")
            if not matrix:
                raise ValueError(
                    f"alignment row {query_name!r} has hotspot index "
                    f"{hotspot_index}, which is absent from the graph list"
                )
        grouped.setdefault((matrix, row.contig), []).append(row)
    return {key: tuple(rows) for key, rows in grouped.items()}


def _alignment_row_locus_operations(row, locus: core.Coord):
    if not row.graph_path or row.graph_path == "*":
        return []
    required = (
        "query_name", "graph_cigar", "ref_positions", "query_positions",
    )
    if any(not hasattr(row, field) for field in required):
        return []
    try:
        qstart, qend = oriented_query_interval(row, locus)
        completed = complete_query_coverage(
            WholeGraphCigar(
                row.query_name,
                row.graph_path,
                row.graph_cigar,
                row.ref_positions,
                row.query_positions,
            ),
            row.end - row.start,
        )
        sliced = slice_graph_cigar_by_query(
            completed.graph_cigar, qstart, qend,
        )
        return [
            operation
            for segment in core.parse_graphic_segments(
                sliced, "alignment-row selection",
            )
            for operation in segment.ops
        ]
    except Exception:
        return []


def alignment_row_locus_support(row, locus: core.Coord) -> Tuple[int, int]:
    """Return aligned and match-like query bases within one requested locus."""
    operations = _alignment_row_locus_operations(row, locus)
    aligned = sum(
        operation.n for operation in operations
        if operation.op in {"=", "M", "X"}
    )
    match_like = sum(
        operation.n for operation in operations
        if operation.op in {"=", "M"}
    )
    return aligned, match_like


def alignment_row_locus_score(row, locus: core.Coord) -> int:
    """Score one row over LOCUS using the downstream alignment policy."""
    score = 0
    for operation in _alignment_row_locus_operations(row, locus):
        if operation.op in {"=", "M"}:
            score += operation.n
        elif operation.op in {"X", "I", "D"}:
            score -= operation.n + core.EDGE_ALIGNMENT_GAP_PENALTY
    return score


def locate_alignment_row(
    allele: str,
    locus: core.Coord,
    alignment_index: Mapping[Tuple[str, str], Sequence[object]],
):
    matrix = allele_matrix_name(allele)
    candidates = []
    for row_order, row in enumerate(
        alignment_index.get((matrix, locus.chrom), ())
    ):
        if row.start > locus.start or row.end < locus.end:
            continue
        alignment_score = alignment_row_locus_score(row, locus)
        candidates.append((
            0 if row.strand == locus.strand else 1,
            row.end - row.start,
            row.start,
            row.end,
            0 if row.graph_path and row.graph_path != "*" else 1,
            -alignment_score,
            row_order,
            row,
        ))
    if not candidates:
        raise ValueError(
            f"{allele}: no {matrix!r} alignment contains "
            f"{locus.chrom}:{locus.start}-{locus.end}{locus.strand}"
        )
    candidates.sort(key=lambda value: value[:7])
    selected = candidates[0]
    tied_interval = [
        candidate for candidate in candidates
        if candidate[:5] == selected[:5]
    ]
    if len(tied_interval) > 1:
        score_text = ", ".join(
            f"{getattr(candidate[7], 'query_name', '?')}="
            f"{-candidate[5]}"
            for candidate in tied_interval
        )
        chosen = selected[7]
        print(
            f"WARNING: {allele}: {len(tied_interval)} containing alignment "
            f"rows tie on interval {locus.chrom}:"
            f"{chosen.start}-{chosen.end}{chosen.strand}; alignment scores "
            f"[{score_text}]; selected first highest-scoring row "
            f"{getattr(chosen, 'query_name', '?')}",
            file=sys.stderr,
            flush=True,
        )
    return selected[7]


def oriented_query_interval(row, locus: core.Coord) -> Tuple[int, int]:
    if row.strand == "+":
        return locus.start - row.start, locus.end - row.start
    return row.end - locus.end, row.end - locus.start


def _extension_value(token: str) -> str:
    """Return the numeric part of an optional ``neighbor:value`` token."""
    token = (token or "").strip()
    if ":" in token:
        token = token.rsplit(":", 1)[1].strip()
    return token


def parse_extension_token(token: str, max_extension: int) -> Tuple[str, int, int]:
    """Return ``(sign, raw_size, effective_size)`` for GenomeLift offsets.

    ``+N`` owns the adjacent gap and receives the first capped extension.
    Unsigned ``N`` is the lower-priority side and receives only the capped
    remainder after the owner.  ``-N`` excludes an overlap from this interval.
    """
    token = _extension_value(token)
    if token in {"", ".", "NA", "None", "none", "null"}:
        return "", 0, 0
    if token in {"-inf", "+inf", "inf"}:
        return "-", 10**18, 10**18
    sign = token[0] if token[:1] in {"+", "-"} else ""
    body = token[1:] if sign else token
    if not body.isdigit():
        raise ValueError(f"malformed GenomeLift extension token {token!r}")
    raw = int(body)
    if sign == "-":
        effective = raw
    elif sign == "+":
        effective = min(raw, max_extension)
    else:
        effective = min(max(0, raw - min(raw, max_extension)), max_extension)
    return sign, raw, effective


def apply_locus_extensions(
    row: GenomeLiftRow, max_extension: int,
) -> core.Coord:
    """Apply GenomeLift columns 13/14 in absolute coordinate order."""
    locus = row.locus
    left_sign, left_raw, left_effective = parse_extension_token(
        row.left_extension, max_extension,
    )
    right_sign, right_raw, right_effective = parse_extension_token(
        row.right_extension, max_extension,
    )
    start = (
        locus.start + left_raw
        if left_sign == "-"
        else max(0, locus.start - left_effective)
    )
    end = (
        locus.end - right_raw
        if right_sign == "-"
        else locus.end + right_effective
    )
    if end <= start:
        raise ValueError(
            f"{row.allele}: GenomeLift extensions "
            f"{row.left_extension!r}/{row.right_extension!r} exclude the "
            f"entire interval {row.locus_text}"
        )
    return core.Coord(locus.chrom, start, end, locus.strand)


def apply_alignment_locus_extensions(
    row: GenomeLiftRow, max_extension: int,
) -> core.Coord:
    """Return the query interval made available to graph alignment.

    Negative GenomeLift offsets describe *ownership*, not sequence that must
    be discarded before alignment.  Keep that sequence here so a lower-
    priority interval can fill bases that the higher-priority neighbor did not
    align.  Positive/unsigned offsets still add their normal bounded context.
    Ownership is resolved later, per aligned interval, by graphreftovcf.py.
    """
    locus = row.locus
    left_sign, _left_raw, left_effective = parse_extension_token(
        row.left_extension, max_extension,
    )
    right_sign, _right_raw, right_effective = parse_extension_token(
        row.right_extension, max_extension,
    )
    start = (
        locus.start
        if left_sign == "-"
        else max(0, locus.start - left_effective)
    )
    end = (
        locus.end
        if right_sign == "-"
        else locus.end + right_effective
    )
    return core.Coord(locus.chrom, start, end, locus.strand)


def prepare_input_record(
    row: GenomeLiftRow,
    alignment_index,
    max_extension: int = 10000,
    cross_validation_extension: int = 0,
    symmetric_extension: Optional[int] = None,
) -> PreparedInputRecord:
    # Align the complete raw query block.  Negative GenomeLift values encode
    # lower-priority overlap ownership and must not pre-trim sequence here.
    # graphreftovcf.py resolves the overlap per aligned base after every row is
    # available.
    core_locus = row.locus
    # Preserve the complete graph CIGAR represented by the selected interval.
    if symmetric_extension is None:
        cigar_locus = apply_alignment_locus_extensions(row, max_extension)
    else:
        symmetric_extension = int(symmetric_extension)
        if symmetric_extension < 0:
            raise ValueError("symmetric extension must be >= 0")
        cigar_locus = core.Coord(
            row.locus.chrom,
            max(0, row.locus.start - symmetric_extension),
            row.locus.end + symmetric_extension,
            row.locus.strand,
        )
    try:
        alignment = locate_alignment_row(
            row.allele, cigar_locus, alignment_index,
        )
    except ValueError as extended_error:
        # Extensions are optional context.  A block can be valid while the
        # requested context extends beyond the whole alignment that generated
        # it.  In that case retain every available extended base by clipping
        # to the alignment containing the original, unextended block.
        try:
            alignment = locate_alignment_row(
                row.allele, row.locus, alignment_index,
            )
        except ValueError:
            raise extended_error
        clipped_start = max(cigar_locus.start, alignment.start)
        clipped_end = min(cigar_locus.end, alignment.end)
        if clipped_end <= clipped_start:
            raise extended_error
        cigar_locus = core.Coord(
            cigar_locus.chrom,
            clipped_start,
            clipped_end,
            cigar_locus.strand,
        )
    whole_alignment = core.Coord(
        alignment.contig, alignment.start, alignment.end, alignment.strand,
    )
    if cross_validation_extension > 0:
        # Keep this context extension separate from GenomeLift ownership.  It
        # exposes nearby sequence for cross-graph validation; the containing
        # whole alignment remains authoritative.
        cigar_locus = core.extend_coord_for_cross_validation(
            cigar_locus, cross_validation_extension,
        )
        cigar_locus = core.Coord(
            cigar_locus.chrom,
            max(cigar_locus.start, whole_alignment.start),
            min(cigar_locus.end, whole_alignment.end),
            cigar_locus.strand,
        )
        if cigar_locus.end <= cigar_locus.start:
            raise ValueError(
                f"{row.allele}: cross-validation context does not intersect "
                "the graph alignment"
            )
    qstart, qend = oriented_query_interval(alignment, cigar_locus)
    core_qstart, core_qend = oriented_query_interval(alignment, core_locus)
    try:
        completed = complete_query_coverage(
            WholeGraphCigar(
                alignment.query_name,
                alignment.graph_path,
                alignment.graph_cigar,
                alignment.ref_positions,
                alignment.query_positions,
            ),
            alignment.end - alignment.start,
        )
        graph_cigar = slice_graph_cigar_by_query(
            completed.graph_cigar, qstart, qend,
        )
        core_graph_cigar = slice_graph_cigar_by_query(
            completed.graph_cigar, core_qstart, core_qend,
        )
    except ValueError as error:
        raise ValueError(
            f"{row.allele}: failed to slice alignment hotspot "
            f"{alignment.hotspot_index} {alignment.contig}:"
            f"{alignment.start}-{alignment.end}{alignment.strand}: {error}"
        ) from error
    expected = cigar_locus.end - cigar_locus.start
    observed = graph_cigar_query_span(graph_cigar)
    if observed != expected:
        raise ValueError(
            f"{row.allele}: locus span {expected} disagrees with sliced "
            f"graph-CIGAR query span {observed}"
        )
    core_expected = core_locus.end - core_locus.start
    core_observed = graph_cigar_query_span(core_graph_cigar)
    if core_observed != core_expected:
        raise ValueError(
            f"{row.allele}: core locus span {core_expected} disagrees with "
            f"sliced graph-CIGAR query span {core_observed}"
        )
    return PreparedInputRecord(
        row.allele,
        graph_cigar,
        alignment.hotspot_index,
        cigar_locus,
        whole_alignment,
        core_graph_cigar,
        core_locus,
    )


def merged_locus(
    rows: Sequence[GenomeLiftRow],
    context: str,
    outer_strand: Optional[str] = None,
) -> core.Coord:
    if not rows:
        raise ValueError(f"{context}: cannot merge an empty row set")
    chroms = {row.locus.chrom for row in rows}
    strands = {row.locus.strand for row in rows}
    if len(chroms) != 1 or (outer_strand is None and len(strands) != 1):
        raise ValueError(
            f"{context}: grouped intervals must share one contig and strand"
        )
    if outer_strand is not None and outer_strand not in {"+", "-"}:
        raise ValueError(f"{context}: invalid outer strand {outer_strand!r}")
    first = rows[0].locus
    return core.Coord(
        first.chrom,
        min(row.locus.start for row in rows),
        max(row.locus.end for row in rows),
        outer_strand or first.strand,
    )


def coord_text(coord: core.Coord) -> str:
    return f"{coord.chrom}:{coord.start}-{coord.end}{coord.strand}"


def oriented_subinterval(outer: core.Coord, inner: core.Coord) -> Tuple[int, int]:
    if inner.chrom != outer.chrom:
        raise ValueError("grouped interval contig mismatch")
    if inner.start < outer.start or inner.end > outer.end:
        raise ValueError("grouped interval lies outside its merged interval")
    if outer.strand == "+":
        return inner.start - outer.start, inner.end - outer.start
    return outer.end - inner.end, outer.end - inner.start


def add_query_gap_to_graph_cigar(
    graph_cigar: str, size: int, *, at_start: bool,
) -> str:
    """Put a query-only gap inside a terminal graph-segment CIGAR.

    A bare ``10I`` cannot be placed between two ``>path:CIGAR`` records: if the
    preceding record ends in ``900H``, concatenating ``900H10I`` makes that H
    internal, so the core parser incorrectly treats the clipped 900 graph bases
    as part of the represented interval.  Keep terminal H operations terminal
    by inserting a prefix gap after the first leading H or a suffix gap before
    the final trailing H.
    """
    size = int(size)
    if size <= 0:
        return graph_cigar
    segment_matches = list(GRAPH_SEGMENT_RE.finditer(graph_cigar))
    if not segment_matches:
        raise ValueError("cannot attach a query gap to an empty graph CIGAR")
    insertion = f"{size}I"
    if at_start:
        body_start = segment_matches[0].end()
        leading_h = re.match(r"\d+H", graph_cigar[body_start:])
        insert_at = body_start + (leading_h.end() if leading_h else 0)
    else:
        trailing_h = re.search(r"\d+H$", graph_cigar)
        insert_at = trailing_h.start() if trailing_h else len(graph_cigar)
    return graph_cigar[:insert_at] + insertion + graph_cigar[insert_at:]


def compose_prepared_graph_cigar(
    rows: Sequence[GenomeLiftRow],
    outer: core.Coord,
    prepared: Mapping[str, PreparedInputRecord],
    context: str,
    *,
    core_only: bool = False,
) -> str:
    """Cover one merged interval with its prepared part graph CIGARs.

    Parts can originate from different graph partitions. Query-coordinate gaps
    become plain insertions, while overlaps are trimmed from the later part.
    """
    pieces = []
    for row in rows:
        record = prepared.get(row.allele)
        if record is None:
            raise KeyError(f"{context}: no prepared graph CIGAR for {row.allele!r}")
        record_locus = (
            record.core_locus if core_only else record.locus
        ) or row.locus
        start, end = oriented_subinterval(outer, record_locus)
        graph_cigar = (
            record.core_graph_cigar if core_only else record.graph_cigar
        ) or record.graph_cigar
        if record_locus.strand != outer.strand:
            graph_cigar = core.reverse_mapping_gcigar_for_reverse_path(
                graph_cigar,
            )
        pieces.append((start, end, row.allele, graph_cigar))
    pieces.sort(key=lambda value: (value[0], value[1], value[2]))

    output: List[str] = []
    leading_gap = 0
    cursor = 0
    for start, end, allele, graph_cigar in pieces:
        if end <= cursor:
            continue
        if start > cursor:
            gap_size = start - cursor
            if output:
                output[-1] = add_query_gap_to_graph_cigar(
                    output[-1], gap_size, at_start=False,
                )
            else:
                leading_gap += gap_size
            cursor = start
        trim_start = max(0, cursor - start)
        trim_end = end - start
        if trim_start or trim_end != graph_cigar_query_span(graph_cigar):
            graph_cigar = slice_graph_cigar_by_query(
                graph_cigar, trim_start, trim_end,
            )
        if leading_gap:
            graph_cigar = add_query_gap_to_graph_cigar(
                graph_cigar, leading_gap, at_start=True,
            )
            leading_gap = 0
        output.append(graph_cigar)
        cursor = end
    total = outer.end - outer.start
    if cursor < total:
        if not output:
            raise ValueError(f"{context}: merged interval contains no graph CIGAR")
        output[-1] = add_query_gap_to_graph_cigar(
            output[-1], total - cursor, at_start=False,
        )
    composed = "".join(output)
    observed = graph_cigar_query_span(composed)
    if observed != total:
        raise ValueError(
            f"{context}: composed graph-CIGAR span {observed}, expected {total}"
        )
    return composed


def merged_prepared_locus(
    rows: Sequence[GenomeLiftRow],
    prepared: Mapping[str, PreparedInputRecord],
    context: str,
    outer_strand: Optional[str] = None,
    *,
    core_only: bool = False,
) -> core.Coord:
    effective_rows: List[GenomeLiftRow] = []
    for row in rows:
        record = prepared.get(row.allele)
        if record is None:
            raise KeyError(f"{context}: no prepared graph CIGAR for {row.allele!r}")
        effective_rows.append(dataclasses.replace(
            row,
            locus=(
                record.core_locus if core_only else record.locus
            ) or row.locus,
        ))
    return merged_locus(effective_rows, context, outer_strand=outer_strand)


def composed_part_slices(
    rows: Sequence[GenomeLiftRow],
    outer: core.Coord,
    prepared: Mapping[str, PreparedInputRecord],
) -> Dict[str, Tuple[int, int, core.Coord]]:
    """Return each part's full interval as offsets on the merged query.

    The merged composition may have overlapping parts.  Retain each full part
    interval as logical ownership metadata; the VCF stage uses it to project
    calls from the one complete merged CIGAR back to the individual blocks.
    """
    result: Dict[str, Tuple[int, int, core.Coord]] = {}
    for row in rows:
        record = prepared.get(row.allele)
        if record is None:
            continue
        locus = record.locus or row.locus
        represented_start, represented_end = oriented_subinterval(outer, locus)
        if outer.strand == "+":
            coord = core.Coord(
                outer.chrom,
                outer.start + represented_start,
                outer.start + represented_end,
                outer.strand,
            )
        else:
            coord = core.Coord(
                outer.chrom,
                outer.end - represented_end,
                outer.end - represented_start,
                outer.strand,
            )
        result[row.allele] = (represented_start, represented_end, coord)
    return result


def _init_prepare_worker(
    query_index, reference_index, max_extension: int,
    cross_validation_extension: int, reference_extension: int,
) -> None:
    global _PREP_QUERY_ALIGNMENT_INDEX
    global _PREP_REFERENCE_ALIGNMENT_INDEX
    global _PREP_MAX_EXTENSION
    global _PREP_CROSS_VALIDATION_EXTENSION
    global _PREP_REFERENCE_EXTENSION
    _PREP_QUERY_ALIGNMENT_INDEX = query_index
    _PREP_REFERENCE_ALIGNMENT_INDEX = reference_index
    _PREP_MAX_EXTENSION = max_extension
    _PREP_CROSS_VALIDATION_EXTENSION = cross_validation_extension
    _PREP_REFERENCE_EXTENSION = reference_extension


def _prepare_input_worker(task: PreparationTask) -> PreparedInputRecord:
    alignment_index = (
        _PREP_REFERENCE_ALIGNMENT_INDEX
        if task.is_reference else _PREP_QUERY_ALIGNMENT_INDEX
    )
    if task.is_reference:
        if task.alignment is not None:
            return prepare_seeded_reference_record(
                task.row, task.alignment,
            )
        return prepare_input_record(
            task.row,
            alignment_index,
            cross_validation_extension=0,
            symmetric_extension=_PREP_REFERENCE_EXTENSION,
        )
    return prepare_input_record(
        task.row,
        alignment_index,
        _PREP_MAX_EXTENSION,
        _PREP_CROSS_VALIDATION_EXTENSION,
    )


def prepare_seeded_reference_record(
    seed_row: GenomeLiftRow,
    alignment,
) -> PreparedInputRecord:
    """Prepare the complete graph row selected by a column-10 seed.

    ``seed_row`` is deliberately retained as ``core_*`` metadata so callers
    can report and rank the exact intersection that selected this candidate.
    It is not an alignment boundary: the primary graph comparison uses the
    complete ``graph_cigar`` and ``locus`` returned here.
    """
    alignment_locus = core.Coord(
        alignment.contig, alignment.start, alignment.end, alignment.strand,
    )
    exact_index = {
        (allele_matrix_name(seed_row.allele), alignment.contig): (alignment,)
    }
    seed_record = prepare_input_record(
        seed_row,
        exact_index,
        cross_validation_extension=0,
        symmetric_extension=0,
    )
    full_row = dataclasses.replace(
        seed_row,
        locus=alignment_locus,
        locus_text=coord_text(alignment_locus),
    )
    full_record = prepare_input_record(
        full_row,
        exact_index,
        cross_validation_extension=0,
        symmetric_extension=0,
    )
    return dataclasses.replace(
        full_record,
        core_graph_cigar=seed_record.graph_cigar,
        core_locus=seed_record.locus,
    )


def _prepare_input_records(
    tasks: Sequence[PreparationTask],
    query_index,
    reference_index,
    processes: int,
    chunksize: int,
    start_method: str,
    max_extension: int,
    cross_validation_extension: int = 0,
    reference_extension: int = DEFAULT_REFERENCE_EXTENSION,
) -> Dict[str, PreparedInputRecord]:
    if processes <= 1 or len(tasks) <= 1:
        records = {}
        for task in tasks:
            index = reference_index if task.is_reference else query_index
            if task.is_reference:
                if task.alignment is not None:
                    record = prepare_seeded_reference_record(
                        task.row, task.alignment,
                    )
                else:
                    record = prepare_input_record(
                        task.row, index,
                        cross_validation_extension=0,
                        symmetric_extension=reference_extension,
                    )
            else:
                record = prepare_input_record(
                    task.row, index, max_extension,
                    cross_validation_extension,
                )
            records[record.allele] = record
        return records

    worker_count = min(processes, len(tasks))
    actual_chunksize = (
        chunksize if chunksize > 0
        else max(1, min(64, (len(tasks) + worker_count * 8 - 1) // (worker_count * 8)))
    )
    context = mp.get_context(start_method)
    global _PREP_QUERY_ALIGNMENT_INDEX
    global _PREP_REFERENCE_ALIGNMENT_INDEX
    global _PREP_MAX_EXTENSION
    global _PREP_CROSS_VALIDATION_EXTENSION
    global _PREP_REFERENCE_EXTENSION
    _PREP_QUERY_ALIGNMENT_INDEX = query_index
    _PREP_REFERENCE_ALIGNMENT_INDEX = reference_index
    _PREP_MAX_EXTENSION = max_extension
    _PREP_CROSS_VALIDATION_EXTENSION = cross_validation_extension
    _PREP_REFERENCE_EXTENSION = reference_extension
    initializer = None
    initargs = ()
    if start_method != "fork":
        initializer = _init_prepare_worker
        initargs = (
            query_index, reference_index, max_extension,
            cross_validation_extension, reference_extension,
        )
    records: Dict[str, PreparedInputRecord] = {}
    with context.Pool(
        processes=worker_count,
        initializer=initializer,
        initargs=initargs,
    ) as pool:
        for record in pool.imap(_prepare_input_worker, tasks, actual_chunksize):
            records[record.allele] = record
    return records


def complete_alignment_record(alignment) -> PreparedInputRecord:
    """Represent one complete ``_align.txt`` row as an input record."""
    completed = complete_query_coverage(
        WholeGraphCigar(
            alignment.query_name,
            alignment.graph_path,
            alignment.graph_cigar,
            alignment.ref_positions,
            alignment.query_positions,
        ),
        alignment.end - alignment.start,
    )
    locus = core.Coord(
        alignment.contig, alignment.start, alignment.end, alignment.strand,
    )
    return PreparedInputRecord(
        alignment.query_name,
        completed.graph_cigar,
        alignment.hotspot_index,
        locus,
        locus,
        completed.graph_cigar,
        locus,
    )


def local_template_path_overlap(alignment, template: LocalTemplate) -> int:
    """Return graph-template bases traversed by one complete alignment."""
    try:
        record = complete_alignment_record(alignment)
        template_key = core.path_key(template.graph_path)
        template_length = template.end - template.start
        return sum(
            max(0, min(segment.end, template_length) - max(segment.start, 0))
            for segment in core.parse_graphic_segments(
                record.graph_cigar, alignment.query_name,
            )
            if core.path_key(segment.path) == template_key
        )
    except Exception:
        return 0


def add_local_template_groups(
    args: argparse.Namespace,
    graph_entries: Sequence[GraphEntry],
    query_alignments: Sequence[object],
    templates: Sequence[LocalTemplate],
    group_specs: List[tuple],
    prepared: Dict[str, PreparedInputRecord],
    reference_output_name_by_allele: Dict[str, str],
    reference_candidate_priority_by_allele: Dict[str, int],
    reference_candidate_source_by_allele: Dict[str, str],
) -> int:
    """Add one best-overlap query comparison for each local template.

    A partition can contain more than one reference-free template, and two
    distinct alignment rows can carry the same historical PA name.  Preserve
    the first occurrence for compatibility, but give every later emitted PA
    the next unused numeric suffix within its graph/sample/haplotype.  The
    alias is local to these generated comparisons; the source alignment row is
    otherwise left unchanged.
    """
    entry_by_prefix = {entry.graph_prefix: entry for entry in graph_entries}
    alignments_by_graph: Dict[str, List[object]] = {}
    for alignment in query_alignments:
        query_name = str(getattr(alignment, "query_name", "") or "").strip()
        if not query_name:
            raise ValueError("alignment row has no query name")
        graph_prefix = query_name.split("_", 1)[0]
        alignments_by_graph.setdefault(graph_prefix, []).append(alignment)

    reserved_query_names = {
        str(getattr(alignment, "query_name", "") or "").strip()
        for alignment in query_alignments
    }
    used_output_query_names = {
        row.allele
        for specification in group_specs
        for row in specification[0]
    }
    next_index_by_identity: Dict[Tuple[str, str, str], int] = {}
    for name in reserved_query_names | used_output_query_names:
        try:
            group_name, sample_name, haplotype, index, _genome = (
                parse_allele_name(name)
            )
        except (TypeError, ValueError):
            continue
        if not index.isdigit():
            continue
        identity = (group_name, sample_name, haplotype)
        next_index_by_identity[identity] = max(
            next_index_by_identity.get(identity, 1), int(index) + 1,
        )

    def unique_output_query_name(source_name: str) -> str:
        if source_name not in used_output_query_names:
            used_output_query_names.add(source_name)
            return source_name
        group_name, sample_name, haplotype, _index, _genome = (
            parse_allele_name(source_name)
        )
        identity = (group_name, sample_name, haplotype)
        index = next_index_by_identity.get(identity, 1)
        while True:
            candidate = (
                f"{group_name}_{sample_name}_{haplotype}_{index}"
            )
            index += 1
            if (
                candidate not in reserved_query_names
                and candidate not in used_output_query_names
            ):
                next_index_by_identity[identity] = index
                reserved_query_names.add(candidate)
                used_output_query_names.add(candidate)
                return candidate

    added = 0
    for index, template in enumerate(templates, 1):
        entry = entry_by_prefix.get(template.graph_prefix)
        if entry is None:
            raise ValueError(
                f"local template {template.name!r} refers to unknown graph "
                f"prefix {template.graph_prefix!r}"
            )
        candidates = []
        for alignment in alignments_by_graph.get(template.graph_prefix, ()):
            overlap = local_template_path_overlap(alignment, template)
            if overlap <= 0:
                continue
            candidates.append((
                -overlap,
                -(alignment.end - alignment.start),
                alignment.query_name,
                alignment.contig,
                alignment.start,
                alignment.end,
                alignment,
            ))
        if not candidates:
            continue
        alignment = min(candidates)[-1]
        source_query_record = complete_alignment_record(alignment)
        query_locus = source_query_record.locus
        if query_locus is None:
            raise ValueError(f"{alignment.query_name}: missing alignment locus")
        output_query_name = unique_output_query_name(alignment.query_name)
        query_record = dataclasses.replace(
            source_query_record, allele=output_query_name,
        )
        query_row = GenomeLiftRow(
            allele=output_query_name,
            locus=query_locus,
            locus_text=coord_text(query_locus),
            assigned_ref="",
            class_type="LocalTemplate",
            assignment_source="local_template",
        )

        reference_name = (
            f"{template.graph_prefix}_{args.refhaplo}_"
            f"localtemplate{index:09d}"
        )
        reference_locus = core.Coord(
            template.output_contig, template.start, template.end, "+",
        )
        reference_gcigar = (
            f">{template.graph_path}:"
            f"{template.end - template.start}="
        )
        reference_row = GenomeLiftRow(
            allele=reference_name,
            locus=reference_locus,
            locus_text=coord_text(reference_locus),
            assigned_ref=reference_name,
            best_ref=reference_name,
            class_type="LocalTemplate",
            assignment_source="local_template",
            assigned_refs=(reference_name,),
            merged_ref_locus=reference_locus,
            merged_ref_locus_text=coord_text(reference_locus),
        )
        reference_record = PreparedInputRecord(
            reference_name,
            reference_gcigar,
            entry.hotspot_index,
            reference_locus,
            reference_locus,
            reference_gcigar,
            reference_locus,
        )
        old_query = prepared.get(query_record.allele)
        if old_query is not None and old_query != query_record:
            raise ValueError(
                f"query alignment name {query_record.allele!r} resolves to "
                "conflicting complete rows"
            )
        prepared[query_record.allele] = query_record
        prepared[reference_name] = reference_record
        reference_output_name_by_allele[reference_name] = template.name
        # Tier 1 accepts a usable partial comparison, which is important when
        # the sample carries a large insertion/deletion relative to template.
        reference_candidate_priority_by_allele[reference_name] = 1
        reference_candidate_source_by_allele[reference_name] = "local_template"
        group_specs.append((
            (query_row,), (reference_row,), reference_locus, False,
            ((query_row.allele, template.name),),
        ))
        added += 1
    return added


def build_direct_inputs(
    args: argparse.Namespace,
    graph_entries: Sequence[GraphEntry],
) -> Tuple[List[core.PairRow], Dict[str, str]]:
    query_genome = args.query_genome or infer_query_genome(args.align_query)
    reference_alignment = args.align_ref or infer_reference_alignment(
        args.align_query, args.refhaplo,
    )
    for path in (
        args.genomelift, args.align_query, reference_alignment, args.graph_list,
    ):
        if not os.path.isfile(path):
            raise FileNotFoundError(path)

    local_template_path = getattr(args, "local_reference_templates", "")
    if local_template_path and not os.path.isfile(local_template_path):
        raise FileNotFoundError(local_template_path)
    local_templates = (
        read_templates(local_template_path)
        if local_template_path else []
    )

    lift_rows = read_genomelift(args.genomelift)
    (
        selected_query_rows,
        query_row_count,
        assignment_counts,
    ) = select_effective_query_rows(lift_rows, query_genome, args.refhaplo)
    missing_reference_rows = query_row_count - len(selected_query_rows)
    if not selected_query_rows and not local_templates:
        raise ValueError(
            f"no assigned GenomeLift rows found for query genome {query_genome!r}"
        )

    query_alignments = read_alignment_output(args.align_query)
    reference_alignments = read_alignment_output(reference_alignment)
    query_index = build_alignment_index(query_alignments, graph_entries)
    reference_index = build_alignment_index(reference_alignments, graph_entries)
    usable_query_rows: List[GenomeLiftRow] = []
    for row in selected_query_rows:
        locate_alignment_row(row.allele, row.locus, query_index)
        usable_query_rows.append(row)
    selected_query_rows = usable_query_rows
    if not selected_query_rows and not local_templates:
        raise ValueError(
            f"no convertible GenomeLift rows found for query genome "
            f"{query_genome!r}"
        )

    grouped_queries: List[List[GenomeLiftRow]] = []
    group_by_token: Dict[str, List[GenomeLiftRow]] = {}
    for query_row in selected_query_rows:
        if query_row.group_token:
            group = group_by_token.get(query_row.group_token)
            if group is None:
                group = []
                group_by_token[query_row.group_token] = group
                grouped_queries.append(group)
            group.append(query_row)
        else:
            grouped_queries.append([query_row])

    # Column 10 supplies the primary placement, while the reference GenomeLift
    # rows supply stable allele labels.  Column 2 can also locate a lower-
    # priority, same-graph reference candidate when the genomic placement is
    # absent or yields no usable alignment.
    reference_lift_rows: Dict[Tuple[str, str], List[GenomeLiftRow]] = {}
    reference_lift_by_allele: Dict[str, GenomeLiftRow] = {}
    for row in lift_rows:
        try:
            row_genome = parse_allele_name(row.allele)[4]
        except ValueError:
            continue
        if row_genome != args.refhaplo:
            continue
        reference_lift_by_allele.setdefault(row.allele, row)
        reference_lift_rows.setdefault(
            (allele_matrix_name(row.allele), row.locus.chrom), [],
        ).append(row)
    for rows in reference_lift_rows.values():
        rows.sort(key=lambda row: (
            row.locus.start, row.locus.end, row.locus.strand, row.allele,
        ))

    group_specs = []
    reference_alignment_by_allele = {}
    reference_output_name_by_allele = {}
    reference_candidate_priority_by_allele = {}
    reference_candidate_retry_priority_by_allele = {}
    reference_candidate_last_retry_priority_by_allele = {}
    reference_candidate_source_by_allele = {}
    reference_slice_number = 0
    groups_without_reference_candidate = 0
    bestref_similarity_candidates = 0
    preprepared: Dict[str, PreparedInputRecord] = {}
    for query_rows_unsorted in grouped_queries:
        query_rows = (
            sorted(query_rows_unsorted, key=lambda row: row.part_index)
            if query_rows_unsorted[0].part_index
            else query_rows_unsorted
        )
        if query_rows[0].part_index:
            observed_parts = [row.part_index for row in query_rows]
            if len(set(observed_parts)) != len(observed_parts):
                raise ValueError(
                    f"{query_rows[0].allele}: duplicate grouped part indexes "
                    f"within graph {allele_matrix_name(query_rows[0].allele)}: "
                    f"{observed_parts}"
                )
            expected_parts = list(range(1, len(query_rows) + 1))
            if observed_parts != expected_parts:
                # GenomeLift part numbers describe the complete positional
                # group, which can contain alleles from several local graphs.
                # This converter deliberately separates those graphs before
                # reference intersection, so a graph-local subset can legally
                # contain original parts such as [2, 4]. Re-number only that
                # local subset while preserving its original order.
                query_rows = [
                    dataclasses.replace(row, part_index=index)
                    for index, row in enumerate(query_rows, 1)
                ]
        placements = {
            row.merged_ref_locus for row in query_rows
            if row.merged_ref_locus is not None
        }
        if len(placements) > 1:
            raise ValueError(
                f"{query_rows[0].allele}: grouped query parts disagree on "
                "their column-10 reference interval"
            )
        matrices = {allele_matrix_name(row.allele) for row in query_rows}
        if len(matrices) != 1:
            raise ValueError(
                f"{query_rows[0].allele}: grouped query parts cross graphs"
            )
        matrix = next(iter(matrices))
        candidate_group = (tuple(row.allele for row in query_rows),)
        candidate_specs = []
        seen_alignment_candidates = set()

        def add_seed_candidates(
            seed: core.Coord,
            priority: int,
            source: str,
            preferred_output_name: Optional[str] = None,
        ) -> None:
            alignments = sorted(
                reference_index.get((matrix, seed.chrom), ()),
                key=lambda row: (
                    row.start, row.end, row.strand, row.hotspot_index,
                    row.query_name, row.graph_path, row.graph_cigar,
                ),
            )
            for alignment in alignments:
                start = max(seed.start, alignment.start)
                end = min(seed.end, alignment.end)
                if end <= start:
                    continue
                alignment_key = (
                    alignment.contig, alignment.start, alignment.end,
                    alignment.strand, alignment.hotspot_index,
                    alignment.query_name, alignment.graph_path,
                    alignment.graph_cigar,
                )
                # Deduplicate within one decision tier. The same physical row
                # must remain available in a later tier because tier 2 has a
                # >50% success rule while tiers 1 and 3 intentionally do not.
                tiered_alignment_key = (priority, alignment_key)
                if tiered_alignment_key in seen_alignment_candidates:
                    continue
                seen_alignment_candidates.add(tiered_alignment_key)
                output_name = preferred_output_name
                if output_name is None:
                    labelled_reference_rows = []
                    for labelled_row in reference_lift_rows.get(
                        (matrix, seed.chrom), (),
                    ):
                        overlap = min(end, labelled_row.locus.end) - max(
                            start, labelled_row.locus.start,
                        )
                        if overlap > 0:
                            labelled_reference_rows.append((
                                -overlap,
                                labelled_row.locus.end - labelled_row.locus.start,
                                labelled_row.locus.start,
                                labelled_row.locus.end,
                                labelled_row.allele,
                            ))
                    if labelled_reference_rows:
                        output_name = min(labelled_reference_rows)[-1]
                candidate_specs.append((
                    priority, source, seed, alignment, start, end, output_name,
                ))

        # Priorities 0 and 2 share one column-2 alignment. It is accepted first
        # above 90%; if not, the cached result is reconsidered above 50% only
        # after the unmerged column-9 tier has had its chance.
        for query_row in query_rows:
            for best_ref in query_row.best_ref.split(";"):
                best_ref = best_ref.strip()
                if not best_ref or allele_matrix_name(best_ref) != matrix:
                    continue
                reference_row = reference_lift_by_allele.get(best_ref)
                if reference_row is None:
                    continue
                before = len(candidate_specs)
                add_seed_candidates(
                    reference_row.locus,
                    0,
                    "bestref_similarity",
                    preferred_output_name=best_ref,
                )
                bestref_similarity_candidates += len(candidate_specs) - before

        # Priority 1: exactly one query row and one unmerged column-9
        # reference. The named reference row supplies the seed and label.
        if len(query_rows) == 1 and len(query_rows[0].assigned_refs) == 1:
            assigned_ref = query_rows[0].assigned_refs[0]
            if allele_matrix_name(assigned_ref) == matrix:
                reference_row = reference_lift_by_allele.get(assigned_ref)
                if reference_row is not None:
                    add_seed_candidates(
                        reference_row.locus,
                        1,
                        "column9_single",
                        preferred_output_name=assigned_ref,
                    )

        # Priority 3: the current column-10 coordinate-intersection behavior,
        # including merged and unnamed placements. It runs only if tiers 1 and
        # 2 fail their respective success rules.
        placement = next(iter(placements)) if placements else None
        if placement is not None:
            add_seed_candidates(placement, 3, "column10")

        if not candidate_specs:
            groups_without_reference_candidate += 1
            continue

        # Each intersecting _align.txt interval is a separate candidate. Its
        # seed is not a hard boundary: alignment uses the complete reference
        # graph row. Candidate priority implements the explicit roll-down.
        for (
            candidate_priority, candidate_source, seed, alignment, start, end,
            preferred_output_name,
        ) in candidate_specs:
            reference_slice_number += 1
            source_tag = {
                0: "bestref",
                1: "col9",
                3: "col10",
            }.get(candidate_priority, "candidate")
            reference_name = (
                f"{matrix}_{args.refhaplo}_{source_tag}"
                f"{reference_slice_number:09d}"
            )
            reference_locus = core.Coord(
                seed.chrom, start, end, alignment.strand,
            )
            reference_alignment_locus = core.Coord(
                alignment.contig,
                alignment.start,
                alignment.end,
                alignment.strand,
            )
            reference_output_name = (
                preferred_output_name or reference_name
            )
            reference_row = GenomeLiftRow(
                allele=reference_name,
                locus=reference_locus,
                locus_text=coord_text(reference_locus),
                assigned_ref=reference_name,
                best_ref=reference_name,
                class_type="Ref",
                assignment_source=candidate_source,
                assigned_refs=(reference_name,),
                merged_ref_locus=reference_locus,
                merged_ref_locus_text=coord_text(reference_locus),
            )
            reference_alignment_by_allele[reference_name] = alignment
            reference_output_name_by_allele[reference_name] = (
                reference_output_name
            )
            reference_candidate_priority_by_allele[reference_name] = (
                candidate_priority
            )
            if candidate_source == "bestref_similarity":
                reference_candidate_retry_priority_by_allele[reference_name] = 2
                # Reuse the same expensive alignment once more, after all
                # positional candidates fail. A substantial absolute match is
                # useful for otherwise missing long query alleles even when it
                # covers less than half of the query interval.
                reference_candidate_last_retry_priority_by_allele[
                    reference_name
                ] = 4
            reference_candidate_source_by_allele[reference_name] = (
                candidate_source
            )
            group_specs.append((
                tuple(query_rows), (reference_row,), reference_alignment_locus,
                len(query_rows) > 1, candidate_group,
            ))

    local_group_count = add_local_template_groups(
        args,
        graph_entries,
        query_alignments,
        local_templates,
        group_specs,
        preprepared,
        reference_output_name_by_allele,
        reference_candidate_priority_by_allele,
        reference_candidate_source_by_allele,
    )
    if local_templates:
        sys.stderr.write(
            "[graphcigartoref_persample] local-template fallback: "
            f"selected {local_group_count}/{len(local_templates)} "
            "reference-free template comparison(s) by maximum graph-path "
            "interval overlap\n"
        )

    if groups_without_reference_candidate:
        sys.stderr.write(
            "[graphcigartoref_persample] skipped "
            f"{groups_without_reference_candidate} query group(s) with "
            "no column-9, column-2, or column-10 same-graph reference "
            "_align.txt candidate\n"
        )
    if bestref_similarity_candidates:
        sys.stderr.write(
            "[graphcigartoref_persample] registered "
            f"{bestref_similarity_candidates} priority-0/2/4 column-2 same-graph "
            "similarity candidate(s) (one cached alignment each)\n"
        )
    task_by_allele: Dict[str, PreparationTask] = {}
    for query_rows, reference_rows, _, _, _ in group_specs:
        for query_row in query_rows:
            if query_row.allele in preprepared:
                continue
            task_by_allele.setdefault(
                query_row.allele, PreparationTask(query_row, False),
            )
        for ref_row in reference_rows:
            if ref_row.allele in preprepared:
                continue
            task_by_allele.setdefault(
                ref_row.allele,
                PreparationTask(
                    ref_row, True,
                    reference_alignment_by_allele[ref_row.allele],
                ),
            )
    tasks = list(task_by_allele.values())
    workers = min(args.processes, len(tasks))
    sharing = (
        "serial" if workers <= 1 else
        "shared fork indexes" if args.mp_start_method == "fork" else
        f"{args.mp_start_method} copies"
    )
    sys.stderr.write(
        f"[graphcigartoref_persample] preparing {len(tasks)} unique allele "
        f"records with {workers} process{'es' if workers != 1 else ''} "
        f"({sharing})\n"
    )
    prepared = _prepare_input_records(
        tasks, query_index, reference_index,
        args.processes, args.chunksize, args.mp_start_method,
        getattr(args, "max_extension", 10000),
        getattr(args, "cross_validation_extension", 4000),
        getattr(args, "reference_extension", DEFAULT_REFERENCE_EXTENSION),
    )
    prepared.update(preprepared)
    sys.stderr.write(
        "[graphcigartoref_persample] reference cascade: column-2 same-graph "
        "similarity (>90% aligned query coverage), one unmerged column-9 "
        "allele (any usable alignment), cached column-2 similarity (>50%), "
        "column-10 coordinate candidates (any usable alignment), then cached "
        "column-2 similarity (>=2000 aligned unmasked query bases) as a last "
        "resort. Seeds "
        "select complete graph rows and never hard-crop alignment to the seed\n"
    )
    cross_validation_extension = getattr(args, "cross_validation_extension", 4000)
    if cross_validation_extension:
        sys.stderr.write(
            "[graphcigartoref_persample] added "
            f"{cross_validation_extension}-bp symmetric query context for "
            "cross-graph validation\n"
        )
    gcigars = {name: record.graph_cigar for name, record in prepared.items()}
    pairs: List[core.PairRow] = []
    line_number = 0
    grouped_pair_count = 0
    for (
        query_rows_all, reference_rows_all, reference_outer_hint, grouped,
        candidate_group,
    ) in group_specs:
        query_rows = tuple(
            row for row in query_rows_all if row.allele in prepared
        )
        reference_rows = tuple(
            row for row in reference_rows_all if row.allele in prepared
        )
        if not query_rows or not reference_rows:
            continue
        reference_outer = merged_prepared_locus(
            reference_rows,
            prepared,
            f"reference group for {query_rows[0].allele}",
            outer_strand=reference_outer_hint.strand,
        )
        reference_core_outer = merged_prepared_locus(
            reference_rows,
            prepared,
            f"unextended reference core for {query_rows[0].allele}",
            outer_strand=reference_outer_hint.strand,
            core_only=True,
        )
        ref_name = ";".join(row.allele for row in reference_rows)
        output_ref_name = ";".join(
            reference_output_name_by_allele.get(row.allele, row.allele)
            for row in reference_rows
        )
        if grouped:
            gcigars[ref_name] = compose_prepared_graph_cigar(
                reference_rows,
                reference_outer,
                prepared,
                f"reference group {ref_name}",
            )

        # Build one merged comparison per query orientation.  For a multi-match
        # group the complete CIGAR is written once, together with one ownership
        # interval per original query part.  graphreftovcf expands those labels
        # logically and assigns events back to the individual blocks without
        # discarding the flanking anchor context used by this alignment.
        query_rows_by_strand: Dict[str, List[GenomeLiftRow]] = {}
        for query_row in query_rows:
            query_rows_by_strand.setdefault(query_row.locus.strand, []).append(
                query_row
            )
        for query_strand, strand_rows_list in query_rows_by_strand.items():
            strand_rows = tuple(strand_rows_list)
            query_base_name = ";".join(row.allele for row in strand_rows)
            if grouped:
                query_outer = merged_prepared_locus(
                    strand_rows,
                    prepared,
                    f"query group {query_base_name}",
                    outer_strand=query_strand,
                )
                query_name = (
                    query_base_name
                    if len(query_rows_by_strand) == 1
                    else f"{query_base_name}|outer_strand={query_strand}"
                )
                gcigars[query_name] = compose_prepared_graph_cigar(
                    strand_rows,
                    query_outer,
                    prepared,
                    f"query group {query_base_name} ({query_strand})",
                )
                part_slices = composed_part_slices(
                    strand_rows, query_outer, prepared,
                )
                grouped_pair_count += len(strand_rows)
            else:
                query_name = strand_rows[0].allele
                query_outer = (
                    prepared[query_name].locus or strand_rows[0].locus
                )
                part_slices = {
                    query_name: (0, query_outer.end - query_outer.start, query_outer),
                }
            query_core_outer = merged_prepared_locus(
                strand_rows,
                prepared,
                f"unextended query core {query_base_name}",
                outer_strand=query_strand,
                core_only=True,
            )
            query_core_slice = oriented_subinterval(
                query_outer, query_core_outer,
            )
            reference_core_slice = oriented_subinterval(
                reference_outer, reference_core_outer,
            )

            reference_coord_text = ";".join(
                coord_text(
                    prepared[row.allele].locus
                    or row.locus
                )
                for row in reference_rows
            )
            reference_alignment_regions = ";".join(
                coord_text(
                    prepared[row.allele].alignment_locus
                    or prepared[row.allele].locus
                    or row.locus
                )
                for row in reference_rows
            )
            for query_row in strand_rows:
                part = part_slices.get(query_row.allele)
                if part is None:
                    sys.stderr.write(
                        "[graphcigartoref_persample] warning: skipping "
                        f"fully overlapped grouped part {query_row.allele}\n"
                    )
                    continue
                qstart, qend, qcoord = part
                qrecord = prepared[query_row.allele]
                line_number += 1
                pair = core.PairRow(
                    line_no=line_number,
                    raw_line="",
                    label=query_row.allele,
                    ref_label=ref_name,
                    query_coord_text=coord_text(query_outer),
                    ref_coord_text=coord_text(reference_outer),
                    label_coord=coord_text(qcoord),
                    query_name=query_name,
                    query_coord=query_outer,
                    ref_name=ref_name,
                    ref_coord=reference_outer,
                )
                # These dynamic fields intentionally keep PairRow compatible
                # with graphcigartoref.py while carrying the per-part output
                # view used after the merged comparison is built.
                pair.output_query_name = query_row.allele
                pair.output_ref_name = output_ref_name
                pair.output_query_coord = qcoord
                pair.output_query_core_coord = (
                    qrecord.core_locus or query_row.locus
                )
                pair.output_ref_coord_text = reference_coord_text
                pair.output_query_alignment = coord_text(
                    qrecord.alignment_locus or qrecord.locus or query_row.locus
                )
                pair.output_ref_alignment = reference_alignment_regions
                pair.output_slice = (qstart, qend)
                # The selected tier supplies a seed for one complete same-graph
                # reference candidate. Never route it through the isolated-core
                # aligner: the seed is metadata, not an alignment boundary.
                pair.query_core_slice = None
                pair.ref_core_slice = None
                pair.column10_reference_seed_slice = reference_core_slice
                pair.column10_query_owned_slice = query_core_slice
                pair.column10_candidate_key = (
                    candidate_group, query_strand, query_name,
                    coord_text(query_outer),
                )
                pair.reference_candidate_priority = (
                    reference_candidate_priority_by_allele.get(ref_name, 0)
                )
                pair.reference_candidate_retry_priority = (
                    reference_candidate_retry_priority_by_allele.get(ref_name)
                )
                pair.reference_candidate_last_retry_priority = (
                    reference_candidate_last_retry_priority_by_allele.get(
                        ref_name
                    )
                )
                pair.reference_candidate_source = (
                    reference_candidate_source_by_allele.get(
                        ref_name, "column10",
                    )
                )
                pair.column10_intersection_span = (
                    reference_core_outer.end - reference_core_outer.start
                )
                pair.column10_reference_sort = (
                    reference_outer.chrom,
                    reference_outer.start,
                    reference_outer.end,
                    reference_outer.strand,
                    ref_name,
                )
                pairs.append(pair)
    if missing_reference_rows:
        local_note = (
            f"; {local_group_count} reference-free template comparison(s) "
            "were retained separately"
            if local_group_count else ""
        )
        sys.stderr.write(
            f"[graphcigartoref_persample] skipped "
            f"{missing_reference_rows} query rows at reference selection "
            "because columns 9, 2, and 10 supplied no usable reference "
            "candidate; "
            "after later ownership/alignment filters retained "
            f"{assignment_counts['column9_single']} single-column-9 row(s), "
            f"{assignment_counts['bestref_similarity']} column-2 row(s), and "
            f"{assignment_counts['column10']} column-10 row(s)"
            f"{local_note}\n"
        )
    if grouped_pair_count:
        sys.stderr.write(
            f"[graphcigartoref_persample] prepared {grouped_pair_count} grouped "
            "query part(s) through merged query/reference graph intervals\n"
        )
    return pairs, gcigars


_WORKER_GCIGARS = None
_WORKER_GRAPH_SEQUENCES = None
_WORKER_GRAPH_PATH_COORDS = None
_WORKER_CHROM_LENGTHS = None
_WORKER_GRAPH_MAPPINGS = None
_WORKER_COORD_BY_NAME = None
_WORKER_ALLOWCHROMS = None
_WORKER_ASSEMBLY_SEQUENCES = None
_WORKER_REF_READER = None
_SHARED_REFERENCE_READER = None
_PARENT_REFERENCE_READER = None
_WORKER_QUERY_CACHE = None

# Query slices are copied out of the in-memory assembly.  Keeping one slice for
# every block makes a long run grow without bound (and is especially costly when
# the optional cross-validation context widens every slice).  A bounded cache is
# sufficient because a comparison normally reuses only the few records in its
# current group; cache misses simply take another cheap string slice.
_QUERY_CACHE_MAX_ENTRIES = 256


def _iter_fasta_records(
    path: str, required: Optional[Set[str]] = None,
) -> Iterator[Tuple[str, Tuple[str, ...], str]]:
    """Yield requested ``(name, header_fields, sequence)`` records.

    Supplying normalized path keys avoids joining sequence strings for
    unrelated records while streaming a consolidated cohort FASTA.
    """
    fields: Optional[Tuple[str, ...]] = None
    chunks: List[str] = []
    keep = False
    with open(path, "rt") as handle:
        for raw in handle:
            if raw.startswith(">"):
                if fields is not None and keep:
                    yield fields[0], fields, "".join(chunks)
                parsed = tuple(raw[1:].strip().split())
                if not parsed:
                    raise ValueError(f"{path}: empty FASTA header")
                fields = parsed
                keep = required is None or core.path_key(fields[0]) in required
                chunks = []
            elif raw.strip():
                if fields is None:
                    raise ValueError(f"{path}: sequence occurs before first FASTA header")
                if keep:
                    chunks.append(raw.strip())
    if fields is not None and keep:
        yield fields[0], fields, "".join(chunks)


def _store_alias(mapping: Dict, key: str, value, source: str) -> None:
    # dict.get bypasses ReferenceSliceSequences materialization so stored
    # coordinate markers compare against markers, not fetched sequences.
    previous = dict.get(mapping, key)
    if previous is not None and previous != value:
        raise ValueError(f"conflicting graph record alias {key!r} while reading {source}")
    mapping[key] = value


def _reference_reader_for_sequences():
    reader = (
        _WORKER_REF_READER
        or _SHARED_REFERENCE_READER
        or _PARENT_REFERENCE_READER
    )
    if reader is None:
        raise RuntimeError(
            "no reference reader is active for reference-slice graph "
            "sequences"
        )
    return reader


class _GraphReferenceReader(IndexedFasta):
    """Read graph-build reference slices without a shared seek position."""

    def fetch(self, name, start, end, strand="+"):
        length, offset, line_bases, line_width = self.index[name]
        if not 0 <= start <= end <= length:
            raise ValueError(f"{self.path}: invalid interval {name}:{start}-{end}")
        if start == end:
            return ""
        self.reopen()
        first = offset + (start // line_bases) * line_width + start % line_bases
        last = offset + ((end - 1) // line_bases) * line_width + (end - 1) % line_bases + 1
        sequence = self._pread(last - first, first).replace(b"\n", b"").replace(b"\r", b"").decode("ascii")
        if len(sequence) != end - start:
            raise ValueError(f"{self.path}: FASTA/index mismatch for {name}:{start}-{end}")
        return core.revcomp(sequence) if strand == "-" else sequence


@dataclasses.dataclass(frozen=True)
class _GraphReferenceSlice:
    reader: _GraphReferenceReader
    contig: str
    start: int
    end: int
    strand: str

    def sequence(self):
        return self.reader.fetch(self.contig, self.start, self.end, self.strand)


@contextmanager
def _graph_summary_reference(summary_dir, calling_reader, calling_haplotype):
    """Resolve the package's compaction reference independently of calling.

    Graph builds record the selected reference in inputs/assemblies.json.
    Older builds used the first normalized assembly. No chromosome alias or
    sequence similarity is used to choose a different assembly.
    """
    inputs = os.path.join(os.path.dirname(os.path.abspath(summary_dir)), "inputs")
    metadata_path = os.path.join(inputs, "assemblies.json")
    fasta = fai = haplotype = None
    if os.path.isfile(metadata_path):
        with open(metadata_path) as handle:
            metadata = json.load(handle)
        haplotype = metadata.get("reference_name")
        fasta = metadata.get("reference_fasta")
        if haplotype:
            row = next((row for row in metadata.get("assemblies", ())
                        if row["name"] == haplotype), {})
            fasta = fasta or row.get("fasta")
            fai = row.get("fai")
            if not fasta:
                raise ValueError(f"{metadata_path}: reference {haplotype!r} lacks a FASTA")
    else:
        query_path = os.path.join(inputs, "query_paths.normalized.txt")
        if os.path.isfile(query_path):
            with open(query_path) as handle:
                for raw in handle:
                    if not raw.strip() or raw.lstrip().startswith("#"):
                        continue
                    fields = raw.split()
                    if len(fields) < 2:
                        raise ValueError(f"{query_path}: expected NAME FASTA [FAI]")
                    haplotype, fasta = fields[:2]
                    fai = fields[2] if len(fields) > 2 else None
                    break
    if not fasta:
        yield calling_reader, calling_haplotype
        return

    def input_path(value):
        value = os.path.expanduser(value)
        return os.path.abspath(value if os.path.isabs(value) else os.path.join(inputs, value))

    fasta = input_path(fasta)
    fai = input_path(fai) if fai and fai not in {".", "null", "None"} else fasta + ".fai"
    if os.path.realpath(fasta) == os.path.realpath(getattr(calling_reader, "fasta_path", "")):
        yield calling_reader, haplotype
        return
    for path in (fasta, fai):
        if not os.path.isfile(path):
            raise FileNotFoundError(f"graph-build reference required to restore compacted paths: {path}")
    reader = _GraphReferenceReader(fasta, fai)
    try:
        sys.stderr.write(
            f"[graphcigartoref_persample] restoring graph paths from build reference "
            f"{haplotype or '.'}: {fasta}; calling reference remains {calling_haplotype}\n"
        )
        yield reader, haplotype
    finally:
        # Slice markers retain the index and reopen on demand. Their os.pread
        # reads are safe across fork workers; spawn reopens via IndexedFasta.
        reader.close()


class ReferenceSliceSequences(dict):
    """Graph-path sequences with reference-slice virtualization.

    A value is a sequence string, a calling-reference tuple, or an indexed
    graph-build-reference slice. Compacted paths use the reference selected
    during graph building, while stored FASTA records are virtualized only
    when byte-identical to the calling reference. The parent retains only
    coordinates for either kind of marker.
    """

    __slots__ = ()

    def __getitem__(self, key):
        value = dict.__getitem__(self, key)
        if isinstance(value, _GraphReferenceSlice):
            return value.sequence()
        if type(value) is tuple:
            contig, start, end, strand, needs_upper = value
            sequence = _reference_reader_for_sequences().fetch(
                contig, start, end, strand,
            )
            return sequence.upper() if needs_upper else sequence
        return value

    def get(self, key, default=None):
        try:
            return self[key]
        except KeyError:
            return default

    def close(self):
        readers = {id(value.reader): value.reader for value in dict.values(self)
                   if isinstance(value, _GraphReferenceSlice)}
        for reader in readers.values():
            reader.close()


GRAPH_SUMMARY_FIELDS = {
    "partition", "path_order", "type", "source_haplotype",
    "source_contig", "source_start", "source_end", "strand", "segments",
}
GRAPH_REFERENCE_SLICE_RE = re.compile(r"^(.+):(\d+)-(\d+)([+-])$")


def _summary_public_path_id(row: Mapping[str, str]) -> str:
    """Return the public path ID used by summary/alternatives.fasta."""
    partition = row["partition"]
    source_contig = row["source_contig"]
    start = int(row["source_start"])
    end = int(row["source_end"])
    if start < 0 or end <= start:
        raise ValueError(
            f"{partition}: invalid summary graph interval {start}-{end}"
        )
    partition_fields = partition.split("_")
    partition_id = partition_fields[0]
    if partition_id == "merged" and len(partition_fields) > 1:
        partition_id = partition_fields[1]
    if (
        not partition_id
        or not source_contig
        or any(character.isspace() for character in source_contig)
        or "/" in source_contig
    ):
        raise ValueError(
            f"{partition}: cannot construct summary graph path ID from "
            f"{source_contig!r}"
        )
    return f"{partition_id}_{source_contig}_{start}_{end}"


def _summary_reference_contig(
    source_contig: str, ref_reader: core.FastaRegionReader,
) -> Optional[str]:
    """Resolve the NC/chr aliases used by graph packages and references."""
    candidates = [source_contig]
    nc_to_chr = {
        f"NC_{accession:06d}.1": f"chr{chromosome}"
        for accession, chromosome in zip(range(60925, 60947), range(1, 23))
    }
    nc_to_chr.update({"NC_060947.1": "chrX", "NC_060948.1": "chrY"})
    chr_to_nc = {chromosome: accession for accession, chromosome in nc_to_chr.items()}
    alias = nc_to_chr.get(source_contig) or chr_to_nc.get(source_contig)
    if alias:
        candidates.append(alias)
    return next((name for name in candidates if name in ref_reader.index), None)


def _summary_reference_slice(
    row: Mapping[str, str], ref_reader: core.FastaRegionReader,
) -> Optional[Tuple[str, int, int, str, str]]:
    """Restore a graph-build reference slice after checking its bounds."""
    encoded = row.get("reference_slice", ".")
    if encoded in {None, "", "."}:
        return None
    match = GRAPH_REFERENCE_SLICE_RE.fullmatch(encoded)
    if match is None:
        raise ValueError(f"invalid graph-summary reference_slice {encoded!r}")
    table_contig, start_text, end_text, strand = match.groups()
    start, end = int(start_text), int(end_text)
    expected = int(row["source_end"]) - int(row["source_start"])
    contig = _summary_reference_contig(table_contig, ref_reader)
    if (
        contig is None
        or start < 0
        or end <= start
        or end - start != expected
        or end > ref_reader.index[contig][0]
    ):
        raise ValueError(
            f"invalid graph-summary reference_slice {encoded!r} for "
            f"{expected}-base path"
        )
    sequence = ref_reader.fetch(contig, start, end, strand)
    return contig, start, end, strand, sequence


def _summary_header_mapping(
    fields: Sequence[str], ref_reader: core.FastaRegionReader,
) -> str:
    """Read either compact legacy or six-column packaged graph headers."""
    # Current alternatives.fasta headers put the authoritative graphic CIGAR
    # in column 4. Older graph FASTAs can put it in column 3, so scan both.
    for token in fields[2:]:
        if token.startswith((">", "<")):
            return token
    # Prefer the mapped coordinate (column 3) before the source-haplotype
    # coordinate (column 2) when a package has no explicit mapping CIGAR.
    for token in fields[2:3]:
        if token and token != "." and ";" not in token:
            mapping = core.synthesize_reference_mapping_from_coord(
                token, ref_reader,
            )
            if mapping:
                return mapping
    return (
        core.synthesize_reference_mapping_from_coord(fields[1], ref_reader)
        if len(fields) >= 2 else ""
    )


def _store_graph_record(
    *,
    name: str,
    fields: Sequence[str],
    sequence: str,
    graph_prefix_value: str,
    source: str,
    sequences: Dict[str, str],
    coordinates: Dict[str, core.Coord],
    mappings: Dict[str, str],
    ref_reader: core.FastaRegionReader,
    template_by_graph_path: Mapping[Tuple[str, str], LocalTemplate],
    mapping_override: str = "",
    coord_override: Optional[core.Coord] = None,
    sequence_marker: Optional[object] = None,
) -> str:
    """Store one graph record and return its normalized path key."""
    base = core.path_key(name)
    aliases = tuple(dict.fromkeys((name, core.strip_match_name(name), base)))
    stored_value: object = sequence if sequence_marker is None else sequence_marker
    if sequence_marker is None and sequence and len(fields) >= 2:
        # Reference-derived paths need no stored copy: when the reference
        # slice at the header coordinate is byte-identical to the sequence,
        # keep only the coordinates and materialize on demand.  The
        # verification fetch here is transient.
        primary_names = getattr(ref_reader, "_primary_names", None)
        if primary_names is not None:
            for token in fields[1:7]:
                if not token or token in {".", "reference", "original"}:
                    continue
                slice_coord = core.parse_coord(token)
                if slice_coord is None:
                    continue
                chrom = slice_coord.chrom
                if chrom not in primary_names and ":" in chrom:
                    # Header coordinates may carry a haplotype prefix, e.g.
                    # ``CHM13_h1:NC_060929.1:...``.
                    chrom = chrom.split(":", 1)[1]
                if (
                    chrom not in primary_names
                    or slice_coord.end - slice_coord.start != len(sequence)
                    or slice_coord.start < 0
                    or slice_coord.end > ref_reader.index[chrom][0]
                ):
                    continue
                fetched = ref_reader.fetch(
                    chrom,
                    slice_coord.start,
                    slice_coord.end,
                    slice_coord.strand,
                )
                if fetched == sequence:
                    needs_upper = False
                elif fetched.upper() == sequence:
                    # Soft-masked reference vs uppercased package bytes:
                    # materialize with .upper() so consumers see exactly the
                    # packaged sequence.
                    needs_upper = True
                else:
                    continue
                stored_value = (
                    chrom,
                    slice_coord.start,
                    slice_coord.end,
                    slice_coord.strand,
                    needs_upper,
                )
                break
    for alias in aliases:
        _store_alias(sequences, alias, stored_value, source)

    coord = coord_override
    if coord is None and len(fields) >= 2:
        coord = core.parse_coord(fields[1])
    if coord is not None:
        for alias in aliases:
            _store_alias(coordinates, alias, coord, source)

    mapping = mapping_override or _summary_header_mapping(fields, ref_reader)
    if mapping:
        for alias in aliases:
            _store_alias(mappings, alias, mapping, source)

    template = template_by_graph_path.get((graph_prefix_value, base))
    if template is not None:
        template_coord = core.Coord(
            template.output_contig,
            template.start,
            template.end,
            "+",
        )
        template_mapping = core.synthesize_reference_mapping_from_coord(
            coord_text(template_coord), ref_reader,
        )
        for alias in aliases:
            coordinates[alias] = template_coord
            mappings[alias] = template_mapping
    return base


def _read_required_graph_summary_data(
    summary_dir: str,
    graph_entries: Sequence[GraphEntry],
    required: Set[str],
    ref_reader: core.FastaRegionReader,
    reference_haplotype: str,
    local_templates: Sequence[LocalTemplate],
    sequences: Dict[str, str],
    coordinates: Dict[str, core.Coord],
    mappings: Dict[str, str],
) -> Set[str]:
    with ExitStack() as stack:
        source = []

        def graph_source():
            if not source:
                source.append(stack.enter_context(_graph_summary_reference(
                    summary_dir, ref_reader, reference_haplotype,
                )))
            return source[0]

        return _read_required_graph_summary_data_from_sources(
            summary_dir, graph_entries, required, ref_reader, reference_haplotype,
            local_templates, sequences, coordinates, mappings, graph_source,
        )


def _read_required_graph_summary_data_from_sources(
    summary_dir, graph_entries, required, ref_reader, reference_haplotype,
    local_templates, sequences, coordinates, mappings, graph_source,
) -> Set[str]:
    """Resolve required paths from the consolidated graph-build package.

    Paths are copied from alternatives.fasta unless their summary row carries
    a ``reference_slice``. Such paths, plus an unstored main-reference
    ``original`` row, are reconstructed directly from the reference
    used to build the package. The calling reference may be a different
    assembly. Any unresolved path is left to the legacy per-partition fallback in
    ``read_required_graph_data``.
    """
    summary_path = os.path.join(summary_dir, "local_graphs.tsv")
    alternatives_path = os.path.join(summary_dir, "alternatives.fasta")
    for path in (summary_path, alternatives_path):
        if not os.path.isfile(path):
            raise FileNotFoundError(path)

    entry_by_partition = {entry.graph_name: entry for entry in graph_entries}
    rows_by_partition: Dict[str, List[Dict[str, str]]] = {}
    with open(summary_path, "rt", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        if reader.fieldnames is None or not GRAPH_SUMMARY_FIELDS.issubset(
            reader.fieldnames
        ):
            raise ValueError(
                f"{summary_path}: expected columns "
                f"{sorted(GRAPH_SUMMARY_FIELDS)}"
            )
        for row in reader:
            partition = row["partition"]
            if partition in entry_by_partition:
                rows_by_partition.setdefault(partition, []).append(row)

    row_by_record: Dict[str, Tuple[GraphEntry, Dict[str, str]]] = {}
    for partition, rows in rows_by_partition.items():
        orders = [int(row["path_order"]) for row in rows]
        if len(orders) != len(set(orders)):
            raise ValueError(
                f"{summary_path}: {partition} has duplicate path_order"
            )
        has_refined_reference = any(row["type"] == "reference" for row in rows)
        graph_rows = [
            row for row in rows
            if row["type"] != "original" or not has_refined_reference
        ]
        entry = entry_by_partition[partition]
        for row in graph_rows:
            record_id = _summary_public_path_id(row)
            previous = row_by_record.get(record_id)
            if previous is not None:
                raise ValueError(
                    f"{summary_path}: duplicate materialized graph path "
                    f"{record_id!r} in {previous[0].graph_name!r} and "
                    f"{partition!r}"
                )
            row_by_record[record_id] = (entry, row)

    template_by_graph_path = {
        (template.graph_prefix, core.path_key(template.graph_path)): template
        for template in local_templates
    }
    found: Set[str] = set()
    for name, fields, sequence in _iter_fasta_records(
        alternatives_path, required,
    ):
        base = core.path_key(name)
        if base not in required:
            continue
        entry_and_row = row_by_record.get(base)
        if entry_and_row is None:
            # The package can include selected templates outside the active
            # calling partition list. They are irrelevant to this sample.
            continue
        entry, _row = entry_and_row
        found.add(_store_graph_record(
            name=name,
            fields=fields,
            sequence=sequence,
            graph_prefix_value=entry.graph_prefix,
            source=alternatives_path,
            sequences=sequences,
            coordinates=coordinates,
            mappings=mappings,
            ref_reader=ref_reader,
            template_by_graph_path=template_by_graph_path,
        ))

    # Recreate graph-build reference-backed paths and legacy unstored
    # main-reference originals without reopening individual graph FASTAs.
    for record_id in sorted(required.difference(found)):
        entry_and_row = row_by_record.get(record_id)
        if entry_and_row is None:
            continue
        entry, row = entry_and_row
        if row.get("reference_slice", ".") in {None, "", "."} and row["type"] != "original":
            continue
        graph_reader, graph_haplotype = graph_source()
        compacted = _summary_reference_slice(row, graph_reader)
        if compacted is not None:
            contig, start, end, strand, sequence = compacted
            source_coord = core.Coord(
                f"{row['source_haplotype']}:{row['source_contig']}",
                int(row["source_start"]), int(row["source_end"]),
                row["strand"],
            )
            stored_header = row.get("fasta_header", ".")
            if stored_header in {None, "", "."}:
                raise ValueError(
                    f"{summary_path}: compacted path {record_id} lacks "
                    "fasta_header"
                )
            fields = tuple(stored_header.split())
            if not fields or fields[0] != record_id:
                raise ValueError(
                    f"{summary_path}: compacted fasta_header ID does not "
                    f"match {record_id}"
                )
            mapping = _summary_header_mapping(fields, graph_reader) or (
                core.synthesize_reference_mapping_from_coord(
                    f"{contig}:{start}-{end}{strand}", graph_reader,
                )
            )
            found.add(_store_graph_record(
                name=record_id,
                fields=fields,
                sequence=sequence,
                graph_prefix_value=entry.graph_prefix,
                source=summary_path,
                sequences=sequences,
                coordinates=coordinates,
                mappings=mappings,
                ref_reader=ref_reader,
                template_by_graph_path=template_by_graph_path,
                mapping_override=mapping,
                coord_override=source_coord,
                # The source is the graph-build reference, which may differ
                # from the calling reference. Keep only its coordinates.
                sequence_marker=(
                    (contig, start, end, strand, False)
                    if graph_reader is ref_reader else
                    _GraphReferenceSlice(graph_reader, contig, start, end, strand)
                ),
            ))
            continue
        if (
            row["type"] != "original"
            or row["source_haplotype"] not in {reference_haplotype, graph_haplotype}
        ):
            continue
        original_reader = graph_reader if row["source_haplotype"] == graph_haplotype else ref_reader
        contig = _summary_reference_contig(row["source_contig"], original_reader)
        if contig is None:
            continue
        start = int(row["source_start"])
        end = int(row["source_end"])
        strand = row["strand"]
        if (
            start < 0
            or end <= start
            or end > original_reader.index[contig][0]
            or strand not in {"+", "-"}
        ):
            raise ValueError(
                f"{summary_path}: invalid main-reference interval "
                f"{contig}:{start}-{end}{strand}"
            )
        source_coord = core.Coord(
            f"{row['source_haplotype']}:{row['source_contig']}",
            start,
            end,
            strand,
        )
        mapping = core.synthesize_reference_mapping_from_coord(
            f"{contig}:{start}-{end}{strand}", original_reader,
        )
        found.add(_store_graph_record(
            name=record_id,
            fields=(record_id,),
            # A main-reference original IS the reference interval by
            # construction; store only its coordinates.
            sequence="",
            sequence_marker=(
                (contig, start, end, strand, False)
                if original_reader is ref_reader else
                _GraphReferenceSlice(original_reader, contig, start, end, strand)
            ),
            graph_prefix_value=entry.graph_prefix,
            source=summary_path,
            sequences=sequences,
            coordinates=coordinates,
            mappings=mappings,
            ref_reader=ref_reader,
            template_by_graph_path=template_by_graph_path,
            mapping_override=mapping,
            coord_override=source_coord,
        ))
    return found


def _candidate_graph_entries(
    graph_entries: Sequence[GraphEntry], pairs: Sequence[core.PairRow]
) -> List[GraphEntry]:
    """Prefer the pair matrices, falling back to the complete ordered list."""
    by_prefix = {entry.graph_prefix: entry for entry in graph_entries}
    labels = pair_matrix_labels(pairs)
    if labels and labels.issubset(by_prefix):
        wanted = labels
        return [entry for entry in graph_entries if entry.graph_prefix in wanted]
    unknown = sorted(labels.difference(by_prefix))
    if unknown:
        examples = ", ".join(repr(value) for value in unknown[:5])
        suffix = f" (and {len(unknown) - 5} more)" if len(unknown) > 5 else ""
        sys.stderr.write(
            "[graphcigartoref_persample] warning: pair matrix label(s) not found "
            f"in -L: {examples}{suffix}; searching the full graph list\n"
        )
    return list(graph_entries)


def pair_matrix_labels(pairs: Sequence[core.PairRow]) -> Set[str]:
    """Return every original matrix represented by simple/composite pairs."""
    labels: Set[str] = set()
    for pair in pairs:
        # query_name/ref_name are the actual prepared graph-CIGAR record names.
        # label/ref_label may instead be display/ownership fields, so consult
        # them only for legacy PairRows that have not yet resolved record names.
        for encoded in (
            pair.query_name or pair.label,
            pair.ref_name or pair.ref_label,
        ):
            if not encoded:
                continue
            # Mixed-strand composite query names carry this internal suffix.
            encoded = encoded.split("|outer_strand=", 1)[0]
            for allele in encoded.split(";"):
                allele = allele.strip()
                if allele:
                    labels.add(allele_matrix_name(allele))
    return labels


def required_graph_paths(
    pairs: Sequence[core.PairRow], gcigars: Dict[str, str]
) -> Set[str]:
    required: Set[str] = set()
    for pair in pairs:
        for row_name in (pair.query_name, pair.ref_name):
            if row_name is None:
                continue
            graph_cigar = gcigars.get(row_name)
            if graph_cigar is None:
                raise KeyError(f"no prepared graph CIGAR for {row_name!r}")
            for segment in core.parse_graphic_segments(graph_cigar, row_name):
                required.add(core.path_key(segment.path))
    return required


def read_required_graph_data(
    graph_entries: Sequence[GraphEntry],
    pairs: Sequence[core.PairRow],
    gcigars: Dict[str, str],
    ref_reader: core.FastaRegionReader,
    local_templates: Sequence[LocalTemplate] = (),
    summary_dir: str = "",
    reference_haplotype: str = "",
) -> Tuple[Dict[str, str], Dict[str, core.Coord], Dict[str, str], int]:
    """Read only graph records referenced by the current graph CIGARs.

    Prefer a graph-build ``summary/`` package when supplied. Any path not
    represented there retains the legacy per-partition FASTA fallback.
    """
    required = required_graph_paths(pairs, gcigars)
    if not required:
        return {}, {}, {}, 0

    sequences: Dict[str, str] = ReferenceSliceSequences()
    coordinates: Dict[str, core.Coord] = {}
    mappings: Dict[str, str] = {}
    found: Set[str] = set()
    files_read = 0
    pair_labels = pair_matrix_labels(pairs)
    template_by_graph_path = {
        (template.graph_prefix, core.path_key(template.graph_path)): template
        for template in local_templates
    }

    if summary_dir:
        summary_found = _read_required_graph_summary_data(
            summary_dir,
            graph_entries,
            required,
            ref_reader,
            reference_haplotype,
            local_templates,
            sequences,
            coordinates,
            mappings,
        )
        found.update(summary_found)
        sys.stderr.write(
            "[graphcigartoref_persample] consolidated graph summary resolved "
            f"{len(summary_found)}/{len(required)} required graph path(s); "
            f"{len(required.difference(found))} require partition-FASTA "
            "fallback\n"
        )

    for entry in _candidate_graph_entries(graph_entries, pairs):
        if required.issubset(found):
            break
        if not os.path.isfile(entry.fasta_path):
            # A filtered graph list can retain paths removed upstream. Missing
            # required records are reported below and their complete candidate
            # groups are omitted before conversion.
            continue
        files_read += 1
        for name, fields, sequence in _iter_fasta_records(
            entry.fasta_path, required,
        ):
            base = core.path_key(name)
            if base not in required or base in found:
                continue
            found.add(_store_graph_record(
                name=name,
                fields=fields,
                sequence=sequence,
                graph_prefix_value=entry.graph_prefix,
                source=entry.fasta_path,
                sequences=sequences,
                coordinates=coordinates,
                mappings=mappings,
                ref_reader=ref_reader,
                template_by_graph_path=template_by_graph_path,
            ))

    missing = sorted(required.difference(found))
    if missing:
        examples = ", ".join(repr(value) for value in missing[:10])
        suffix = f" (and {len(missing) - 10} more)" if len(missing) > 10 else ""
        sys.stderr.write(
            "[graphcigartoref_persample] warning: "
            f"{len(missing)} graph path(s) used by the alignments were not "
            f"found via -G/-L; affected candidate groups will be skipped: "
            f"{examples}{suffix}\n"
        )
    return sequences, coordinates, mappings, files_read


def skip_unavailable_graph_candidates(
    pairs: Sequence[core.PairRow],
    gcigars: Mapping[str, str],
    graph_sequences: Mapping[str, object],
) -> List[core.PairRow]:
    """Drop complete candidate groups that need an unavailable graph path."""
    available = {core.path_key(name) for name in graph_sequences}
    kept: List[core.PairRow] = []
    skipped: List[Tuple[str, Tuple[str, ...]]] = []
    for candidate_group in group_column10_candidates(pairs):
        group_pairs = [
            pair
            for _base_pair, output_pairs in candidate_group
            for pair in output_pairs
        ]
        missing = tuple(sorted(
            required_graph_paths(group_pairs, dict(gcigars)).difference(available)
        ))
        if missing:
            label = str(candidate_group[0][0].label)
            skipped.append((label, missing))
            continue
        kept.extend(group_pairs)
    if skipped:
        examples = ", ".join(
            f"{label} ({'|'.join(paths[:3])})"
            for label, paths in skipped[:5]
        )
        suffix = f" (and {len(skipped) - 5} more)" if len(skipped) > 5 else ""
        sys.stderr.write(
            "[graphcigartoref_persample] warning: skipped "
            f"{len(skipped)} candidate group(s) because graph paths listed "
            f"through -L are unavailable under -G: {examples}{suffix}\n"
        )
    return kept


def _graph_cigar_query_span(graph_cigar: str, row_name: str) -> int:
    return sum(
        core.query_consume(operation.op, operation.n)
        for segment in core.parse_graphic_segments(graph_cigar, row_name)
        for operation in segment.ops
    )


def _query_record_sequences(
    pair: core.PairRow,
    gcigars: Dict[str, str],
    assembly_sequences: Dict[str, str],
    cache: Dict[Tuple[str, int, int, str], str],
) -> Dict[str, str]:
    if pair.query_name is None or pair.query_coord is None:
        raise ValueError(
            f"pair {pair.line_no} lacks a resolved query name/coordinate"
        )
    coord = pair.query_coord
    key = (coord.chrom, coord.start, coord.end, coord.strand)
    sequence = cache.get(key)
    if sequence is None:
        contig_sequence = assembly_sequences.get(coord.chrom)
        if contig_sequence is None:
            contig_sequence = assembly_sequences.get(
                core.strip_match_name(coord.chrom)
            )
        if contig_sequence is None:
            raise KeyError(
                f"contig {coord.chrom!r} is missing from --fasta-query"
            )
        if coord.start < 0 or coord.end < coord.start or coord.end > len(contig_sequence):
            raise ValueError(
                f"{pair.query_name}: query interval {pair.query_coord_text} is "
                f"outside contig {coord.chrom!r} length {len(contig_sequence)}"
            )
        sequence = contig_sequence[coord.start:coord.end]
        if coord.strand == "-":
            sequence = core.revcomp(sequence)
        if len(cache) >= _QUERY_CACHE_MAX_ENTRIES:
            # Do not retain thousands of large cross-validation slices in each
            # worker.  This is deliberately a simple bounded cache rather than
            # an LRU: grouped comparisons are adjacent, so clearing at the
            # boundary preserves the useful locality without bookkeeping or
            # changing any alignment result.
            cache.clear()
        cache[key] = sequence
    expected = _graph_cigar_query_span(gcigars[pair.query_name], pair.query_name)
    if len(sequence) != expected:
        raise ValueError(
            f"{pair.query_name}: fetched {len(sequence)} bases from "
            f"{pair.query_coord_text}, but its graph CIGAR consumes "
            f"{expected} query bases"
        )
    return {
        pair.query_name: sequence,
        core.strip_match_name(pair.query_name): sequence,
    }


def truncate_grouped_provisional_row(row_text: str, pair: core.PairRow) -> str:
    """Legacy helper that trims a merged comparison to one ownership part.

    The production path now performs this reduction through
    :func:`format_provisional_row`; retaining this helper preserves the older
    operation for callers/tests that explicitly request it.
    """
    if pair.query_name == pair.label or not pair.label_coord:
        return row_text
    if pair.query_coord is None:
        raise ValueError(f"grouped pair {pair.line_no} has no merged query coordinate")
    slice_coord_text = getattr(pair, "slice_coord_text", pair.label_coord)
    part_coord = core.parse_coord(slice_coord_text)
    if part_coord is None:
        raise ValueError(
            f"grouped pair {pair.line_no} has malformed part coordinate "
            f"{slice_coord_text!r}"
        )
    qstart, qend = oriented_subinterval(pair.query_coord, part_coord)
    fields = row_text.split("\t")
    if len(fields) < 7:
        raise ValueError(
            f"grouped pair {pair.line_no} produced fewer than seven columns"
        )
    fields[6] = slice_graph_cigar_by_query(fields[6], qstart, qend)
    observed = graph_cigar_query_span(fields[6])
    expected = part_coord.end - part_coord.start
    if observed != expected:
        raise ValueError(
            f"grouped pair {pair.line_no} truncated query span {observed}, "
            f"expected {expected}"
        )
    fields[0] = pair.label
    fields[1] = slice_coord_text
    fields[4] = pair.label
    fields[5] = pair.label_coord
    return "\t".join(fields)


def alignment_score_query_bounds(
    graph_cigar: str,
    minimum_score: int = DEFAULT_EDGE_ALIGNMENT_SCORE,
) -> Optional[Tuple[int, int]]:
    """Return query bounds retained by bidirectional alignment scoring."""
    operations = [
        operation
        for segment in core.parse_graphic_segments(
            graph_cigar, "edge alignment scoring",
        )
        for operation in segment.ops
    ]
    bounds = core.alignment_quality_operation_bounds(
        operations, minimum_score=minimum_score,
    )
    if bounds is None:
        return None
    left_index, right_index = bounds
    qstart = sum(
        core.query_consume(operation.op, operation.n)
        for operation in operations[:left_index]
    )
    qend = sum(
        core.query_consume(operation.op, operation.n)
        for operation in operations[:right_index]
    )
    if qend <= qstart:
        return None
    return qstart, qend


def slice_query_coord(
    coord: core.Coord, qstart: int, qend: int, total: int,
) -> core.Coord:
    """Project a query-axis slice back to strand-aware assembly coordinates."""
    qstart = int(qstart)
    qend = int(qend)
    total = int(total)
    if coord.end - coord.start != total:
        raise ValueError(
            f"query coordinate {coord_text(coord)} spans "
            f"{coord.end - coord.start}, but graph CIGAR spans {total}"
        )
    if qstart < 0 or qend <= qstart or qend > total:
        raise ValueError(f"invalid query-coordinate slice {qstart}-{qend}/{total}")
    if coord.strand == "+":
        start = coord.start + qstart
        end = coord.start + qend
    else:
        start = coord.end - qend
        end = coord.end - qstart
    return core.Coord(coord.chrom, start, end, coord.strand)


def _trim_reference_operation_edges(graph_cigar: str, left: int, right: int) -> str:
    """Apply reference-only edge cuts without rewriting retained operations."""
    chunks = []
    contexts = {}
    operation_index = 0
    for index, match in enumerate(re.finditer(r"([<>])([^<>]*)", graph_cigar)):
        direction, value = match.groups()
        name, body = value.split(":", 1) if ":" in value else (None, value)
        context = 0 if index == 0 or name is None else index
        chunk = {"prefix": direction + (name + ":" if name is not None else ""),
                 "head": "", "tail": "", "tokens": [], "drop": False}
        chunks.append(chunk)
        members, tokens = contexts.setdefault(context, ([], []))
        members.append(chunk)
        for token in core._CIGAR_RE.finditer(body):
            size, operation = int(token.group(1)), token.group(2)
            kept = left <= operation_index < right
            item = [token.group(), size, operation, operation_index, kept]
            chunk["tokens"].append(item)
            tokens.append(item)
            operation_index += 1

    for context, (members, tokens) in contexts.items():
        removed = [token for token in tokens if token[2] != "H" and not token[4]]
        if not removed:
            continue
        if any(token[2] != "D" for token in removed):
            raise ValueError("reference-only edge trimming would remove query bases")
        kept = [token for token in tokens if token[2] != "H" and token[4]]
        head = tokens[0] if tokens and tokens[0][2] == "H" else None
        tail = tokens[-1] if tokens and tokens[-1][2] == "H" and tokens[-1] is not head else None
        left_ref = sum(token[1] for token in removed if token[2] in "=MXD" and token[3] < left)
        right_ref = sum(token[1] for token in removed if token[2] in "=MXD" and token[3] >= right)
        if not kept:
            for chunk in members:
                chunk["drop"] = True
            if context == 0:
                # Keep the original main declaration ahead of retained alternatives.
                first = members[0]
                first.update(drop=False, tokens=[], head=f"{(head[1] if head else 0) + left_ref}H")
            continue
        for token in removed:
            token[0] = ""
        if left_ref:
            if head is not None:
                head[0] = f"{head[1] + left_ref}H"
            else:
                members[0]["head"] = f"{left_ref}H"
        if right_ref:
            if tail is not None:
                tail[0] = f"{tail[1] + right_ref}H"
            else:
                members[-1]["tail"] = f"{right_ref}H"

    output = []
    for index, chunk in enumerate(chunks):
        if chunk["drop"]:
            continue
        body = chunk["head"] + "".join(token[0] for token in chunk["tokens"]) + chunk["tail"]
        if not body and index == 0:
            body = "0H"
        if body:
            output.append(chunk["prefix"] + body)
    return "".join(output)


def trim_query_edges_by_alignment_score(
    graph_cigar: str,
    query_coord: core.Coord,
    minimum_score: int = DEFAULT_EDGE_ALIGNMENT_SCORE,
) -> Optional[Tuple[str, core.Coord, int, int]]:
    """Discard low-quality query edges using bidirectional CIGAR scoring.

    The return offsets are in the input query-CIGAR frame.  ``None`` means the
    comparison never accumulates sufficient confidence from both edges.
    Rows with at most one non-H operation bypass score trimming.
    """
    total = graph_cigar_query_span(graph_cigar)
    # Count across the complete row; H and continuation markers do not count.
    operations = [
        operation
        for segment in core.parse_graphic_segments(
            graph_cigar, "single-operation trimming exemption",
        )
        for operation in segment.ops
    ]
    non_h_count = sum(operation.op != "H" for operation in operations)
    if non_h_count <= 1:
        if query_coord.end - query_coord.start != total:
            raise ValueError("query coordinate span does not match CIGAR query span")
        return graph_cigar, query_coord, 0, total
    bounds = core.alignment_quality_operation_bounds(operations, minimum_score=minimum_score)
    if bounds is None:
        return None
    left, right = bounds
    qstart = sum(core.query_consume(op.op, op.n) for op in operations[:left])
    qend = sum(core.query_consume(op.op, op.n) for op in operations[:right])
    trimmed_coord = slice_query_coord(query_coord, qstart, qend, total)
    if qstart == 0 and qend == total:
        trimmed_cigar = (
            graph_cigar
            if all(op.op == "H" for op in operations[:left] + operations[right:])
            else _trim_reference_operation_edges(graph_cigar, left, right)
        )
    else:
        # Keep the existing query-clipping/extension path unchanged.
        trimmed_cigar = core.slice_graph_cigar_by_query(
            graph_cigar, qstart, qend, "alignment-score edge trim",
            include_left_boundary_deletions=False,
        )
    if graph_cigar_query_span(trimmed_cigar) != qend - qstart:
        raise ValueError("operation-edge trimming changed the retained query span")
    return trimmed_cigar, trimmed_coord, qstart, qend


def _align_reference_flank_ops(
    reference_sequence: str,
    query_sequence: str,
    *,
    reverse_from_core: bool,
) -> List[core.CigarOp]:
    """Align one raw-reference/query flank and return traversal-oriented ops."""
    if not reference_sequence and not query_sequence:
        return []
    if not reference_sequence:
        return [core.CigarOp(len(query_sequence), "I", query_sequence)]
    if not query_sequence:
        return [core.CigarOp(len(reference_sequence), "D", "")]

    reference = reference_sequence
    query = query_sequence
    if reverse_from_core:
        # Upstream is aligned core-first. Restore the original far-to-core
        # traversal after alignment so it can be prepended to the core CIGAR.
        reference = core.revcomp(reference)
        query = core.revcomp(query)
    if reference.upper() == query.upper():
        cigar = f"{len(reference)}="
    else:
        if max(len(reference), len(query)) < 1_000:
            aligned, _r0, _q0, _r1, _q1 = core._tm_global_insert_align_dp(
                reference, query,
            )
        else:
            # Pairs of >=1000 bp always try minimap2 first: SSW's complete
            # n*m matrix is intractable for the near-megabase query gaps
            # this stage can see.  When minimap2 finds no usable hit, pairs
            # under 10 kb on both sides still get an exact SSW matrix (at
            # most ~1e8 cells, well under a second); anything larger keeps
            # the conservative unaligned D/I representation, because a big
            # pair minimap2 cannot anchor has no homology worth minutes of
            # quadratic search.
            try:
                aligned = core.minimap2_payload_ops(reference, query)
            except RuntimeError:
                if max(len(reference), len(query)) < 10_000:
                    aligned = core.swspy(reference, query)
                else:
                    aligned = []
                    core._tm_add_op(aligned, len(reference), 'D', '')
                    core._tm_add_op(aligned, len(query), 'I', query)
        cigar = core._tm_format_pairwise_ops(aligned)
    if reverse_from_core:
        cigar = core._reverse_template_cigar(cigar)
    operations = core.parse_cigar_ops(cigar)
    observed_reference = sum(
        operation.n for operation in operations
        if operation.op in {"=", "X", "D"}
    )
    observed_query = sum(
        operation.n for operation in operations
        if operation.op in {"=", "X", "I"}
    )
    if observed_reference != len(reference_sequence):
        raise ValueError(
            "reference-flank alignment changed the reference span: "
            f"{observed_reference}/{len(reference_sequence)}"
        )
    if observed_query != len(query_sequence):
        raise ValueError(
            "reference-flank alignment changed the query span: "
            f"{observed_query}/{len(query_sequence)}"
        )
    return operations


def _reference_flank_gcigar(
    boundary: core.GraphicSegment,
    query_sequence: str,
    ref_reader: core.FastaRegionReader,
    *,
    upstream: bool,
    maximum_reference_extension: int,
) -> str:
    """Align a clipped query flank to raw reference outside one CIGAR edge."""
    query_size = len(query_sequence)
    if query_size <= 0 or maximum_reference_extension <= 0:
        return ""
    reference_size = (
        min(maximum_reference_extension, 3 * query_size)
        if query_size < 100 else maximum_reference_extension
    )
    reference_name = boundary.path
    if reference_name not in ref_reader.index:
        reference_name = core.strip_match_name(reference_name)
    if reference_name not in ref_reader.index:
        return ""
    contig_size = ref_reader.index[reference_name][0]
    direction = boundary.direction
    if upstream:
        if direction == ">":
            start, end = max(0, boundary.start - reference_size), boundary.start
        else:
            start, end = boundary.end, min(contig_size, boundary.end + reference_size)
    else:
        if direction == ">":
            start, end = boundary.end, min(contig_size, boundary.end + reference_size)
        else:
            start, end = max(0, boundary.start - reference_size), boundary.start
    if end <= start:
        return ""
    strand = "+" if direction == ">" else "-"
    reference_sequence = ref_reader.fetch(reference_name, start, end, strand)
    try:
        operations = _align_reference_flank_ops(
            reference_sequence,
            query_sequence,
            reverse_from_core=upstream,
        )
    except (FileNotFoundError, RuntimeError, ValueError):
        # A rescue is optional evidence. Preserve the already score-clipped
        # core when neither in-process SSW nor minimap2 can align this flank.
        return ""
    return core._format_interval_gcigar(
        direction,
        reference_name,
        contig_size,
        start,
        end,
        operations,
    )


def _reference_alignment_regions_from_gcigar(graph_cigar: str) -> str:
    """Return compact genomic regions actually present in a reference CIGAR."""
    by_key: Dict[Tuple[str, str], List[Tuple[int, int]]] = {}
    for segment in core.parse_graphic_segments(
        graph_cigar, "rescued reference alignment regions",
    ):
        if segment.end <= segment.start:
            continue
        strand = "+" if segment.direction == ">" else "-"
        by_key.setdefault((segment.path, strand), []).append(
            (segment.start, segment.end)
        )
    output = []
    for (path, strand), intervals in sorted(by_key.items()):
        intervals.sort()
        start, end = intervals[0]
        for next_start, next_end in intervals[1:]:
            if next_start <= end:
                end = max(end, next_end)
            else:
                output.append(f"{path}:{start}-{end}{strand}")
                start, end = next_start, next_end
        output.append(f"{path}:{start}-{end}{strand}")
    return ";".join(output)


def _trimmed_reference_regions(graph_cigar: str) -> str:
    """Account for main continuations in a coordinate-only view of the CIGAR."""
    main, alternatives = [], []
    for index, match in enumerate(re.finditer(r"([<>])([^<>]*)", graph_cigar)):
        direction, value = match.groups()
        if index == 0:
            main.append(direction + value)
        elif ":" not in value:
            main.append(value)
        else:
            alternatives.append(direction + value)
    return _reference_alignment_regions_from_gcigar("".join(main + alternatives))


def rescue_score_clipped_query_edges(
    graph_cigar: str,
    query_coord: core.Coord,
    query_sequence: Optional[str],
    ref_reader: Optional[core.FastaRegionReader],
    maximum_reference_extension: int = DEFAULT_REFERENCE_FLANK_RESCUE,
) -> Tuple[Optional[Tuple[str, core.Coord, int, int]], Optional[str]]:
    """Realign score-clipped flanks to raw reference and then clip again.

    The first score pass supplies the trusted core. Upstream is aligned in
    reverse-complement orientation so both flank alignments grow away from that
    core. The returned reference-region string also reflects reference bases
    removed by trimming, even when no raw-reference flank was added.
    """
    initial = trim_query_edges_by_alignment_score(graph_cigar, query_coord)
    if initial is None:
        return None, None
    cigar, _trimmed_coord, qstart, qend = initial
    trimmed_regions = None
    if (qstart, qend) == (0, graph_cigar_query_span(graph_cigar)) and cigar != graph_cigar:
        before_regions = _trimmed_reference_regions(graph_cigar)
        after_regions = _trimmed_reference_regions(cigar)
        if before_regions != after_regions:
            trimmed_regions = after_regions
    total = graph_cigar_query_span(graph_cigar)
    if (
        ref_reader is None
        or query_sequence is None
        or len(query_sequence) != total
        or (qstart == 0 and qend == total)
        or maximum_reference_extension <= 0
    ):
        return initial, trimmed_regions
    segments = core.parse_graphic_segments(cigar, "score-clipped core")
    if not segments:
        return initial, trimmed_regions
    upstream = _reference_flank_gcigar(
        segments[0], query_sequence[:qstart], ref_reader,
        upstream=True,
        maximum_reference_extension=maximum_reference_extension,
    ) if qstart else ""
    downstream = _reference_flank_gcigar(
        segments[-1], query_sequence[qend:], ref_reader,
        upstream=False,
        maximum_reference_extension=maximum_reference_extension,
    ) if qend < total else ""
    if not upstream and not downstream:
        return initial, trimmed_regions
    combined = upstream + cigar + downstream
    if graph_cigar_query_span(combined) != total:
        # Raw-reference flank rescue is optional evidence layered on top of
        # an already score-validated graph alignment.  A rare payload or
        # segment-boundary representation can fail to round-trip through the
        # combined graph-CIGAR parser even though the retained core itself is
        # valid.  Never let that optional extension invalidate every reference
        # candidate for the query: discard the inconsistent rescue and keep
        # the original score-clipped core.
        return initial, trimmed_regions
    polished = trim_query_edges_by_alignment_score(combined, query_coord)
    if polished is None:
        return initial, trimmed_regions
    regions = (
        _trimmed_reference_regions(polished[0])
        if polished[2:] == (0, total) and polished[0] != combined
        else _reference_alignment_regions_from_gcigar(polished[0])
    )
    return polished, regions


def format_provisional_row(
    row_text: str,
    pair: core.PairRow,
    query_sequence: Optional[str] = None,
    ref_reader: Optional[core.FastaRegionReader] = None,
    reference_flank_rescue: int = DEFAULT_REFERENCE_FLANK_RESCUE,
) -> Optional[str]:
    """Format one individual output view from a merged provisional row."""
    fields = row_text.split("\t")
    if len(fields) < 7:
        raise ValueError(
            f"pair {pair.line_no} produced fewer than seven output columns"
        )
    # A core-seeded comparison deliberately emits an empty provisional CIGAR
    # when its isolated core never reaches the alignment-confidence score. Such a
    # comparison is unusable by definition; skip it before attempting an
    # ownership/query slice whose requested span cannot exist.
    if not fields[6] or graph_cigar_query_span(fields[6]) == 0:
        return None
    # The core operates on the merged query/reference records.  Reduce that
    # result to the current individual query part before writing it.  The
    # reference name/coordinate may remain a semicolon-separated multi-match;
    # graphreftovcf uses the complete graph CIGAR and preserves that provenance.
    output_name = getattr(pair, "output_query_name", None)
    if output_name:
        slice_start, slice_end = getattr(
            pair, "output_slice", (0, graph_cigar_query_span(fields[6])),
        )
        cigar = slice_graph_cigar_by_query(fields[6], slice_start, slice_end)
        output_coord = getattr(pair, "output_query_coord", None)
        if output_coord is None:
            raise ValueError(f"pair {pair.line_no}: missing output query coordinate")
        expected = output_coord.end - output_coord.start
        observed = graph_cigar_query_span(cigar)
        if observed != expected:
            raise ValueError(
                f"pair {pair.line_no}: part graph-CIGAR span {observed}, "
                f"expected {expected} for {output_name}"
            )
        part_query_sequence = (
            query_sequence[slice_start:slice_end]
            if query_sequence is not None else None
        )
        anchored, rescued_reference_regions = rescue_score_clipped_query_edges(
            cigar,
            output_coord,
            part_query_sequence,
            ref_reader,
            reference_flank_rescue,
        )
        if anchored is None:
            return None
        cigar, output_coord, _trim_start, _trim_end = anchored
        # A full-span ownership slice can already remove terminal D. Carry
        # that reference change through even if score trimming has no work.
        if rescued_reference_regions is None and (slice_start, slice_end) == (0, graph_cigar_query_span(fields[6])):
            before_regions = _trimmed_reference_regions(fields[6])
            after_regions = _trimmed_reference_regions(cigar)
            if before_regions != after_regions:
                rescued_reference_regions = after_regions
        cigar = query_only_graph_cigar(cigar)
        return core.append_tag("\t".join((
            output_name,
            getattr(pair, "output_ref_name", pair.ref_name or fields[2]),
            coord_text(output_coord),
            getattr(pair, "output_ref_coord_text", fields[3]),
            getattr(pair, "output_query_alignment", pair.query_coord_text),
            rescued_reference_regions or getattr(
                pair, "output_ref_alignment", pair.ref_coord_text,
            ),
            cigar,
        )), read_alternative_tag(fields))

    # Compatibility path for callers constructing a direct PairRow themselves.
    # Keep the requested seven-column per-sample layout even in that case.
    query_coord = core.parse_coord(fields[1]) or pair.query_coord
    if query_coord is None:
        raise ValueError(f"pair {pair.line_no}: missing query coordinate")
    anchored, rescued_reference_regions = rescue_score_clipped_query_edges(
        fields[6], query_coord, query_sequence, ref_reader,
        reference_flank_rescue,
    )
    if anchored is None:
        return None
    cigar, query_coord, _trim_start, _trim_end = anchored
    cigar = query_only_graph_cigar(cigar)
    return core.append_tag("\t".join((
        pair.query_name or fields[0],
        pair.ref_name or fields[2],
        coord_text(query_coord),
        fields[3],
        getattr(pair, "query_alignment_regions_text", pair.query_coord_text),
        rescued_reference_regions or getattr(
            pair, "reference_alignment_regions_text", pair.ref_coord_text,
        ),
        cigar,
    )), read_alternative_tag(fields))


def format_merged_comparison_row(
    row_text: str,
    base_pair: core.PairRow,
    output_pairs: Sequence[core.PairRow],
    query_sequence: Optional[str] = None,
    ref_reader: Optional[core.FastaRegionReader] = None,
    reference_flank_rescue: int = DEFAULT_REFERENCE_FLANK_RESCUE,
) -> Optional[str]:
    """Keep one complete multi-match CIGAR plus per-part ownership metadata.

    ``graphreftovcf.py``'s per-sample-v2 reader is deliberately able to turn
    this single physical row into one logical ownership view per query label.
    The CIGAR must remain the complete merged comparison: slicing it here makes
    a short repetitive part lose the flanking anchors that disambiguated its
    placement in the first place.
    """
    fields = row_text.split("\t")
    if len(fields) < 7:
        raise ValueError(
            f"pair {base_pair.line_no} produced fewer than seven output columns"
        )
    if not fields[6] or graph_cigar_query_span(fields[6]) == 0:
        return None
    if (
        base_pair.query_name is None
        or base_pair.ref_name is None
        or base_pair.query_coord is None
    ):
        raise ValueError(
            f"pair {base_pair.line_no}: incomplete merged comparison metadata"
        )

    # build_provisional_row may already have narrowed the comparison relative
    # to the PairRow's prepared context, so the coordinate written beside its
    # CIGAR is authoritative here.
    query_coord = core.parse_coord(fields[1]) or base_pair.query_coord
    if query_coord is None:
        raise ValueError(
            f"pair {base_pair.line_no}: missing merged query coordinate"
        )
    anchored, rescued_reference_regions = rescue_score_clipped_query_edges(
        fields[6], query_coord, query_sequence, ref_reader,
        reference_flank_rescue,
    )
    if anchored is None:
        return None
    cigar, query_coord, _trim_start, _trim_end = anchored
    query_alignment_regions = ";".join(
        getattr(pair, "output_query_alignment", base_pair.query_coord_text)
        for pair in output_pairs
    )
    reference_coordinates = getattr(
        base_pair, "output_ref_coord_text", base_pair.ref_coord_text,
    )
    reference_alignment_regions = (
        rescued_reference_regions
        or getattr(
            base_pair, "output_ref_alignment", base_pair.ref_coord_text,
        )
    )
    cigar = query_only_graph_cigar(cigar)
    return core.append_tag("\t".join((
        base_pair.query_name,
        getattr(base_pair, "output_ref_name", base_pair.ref_name),
        coord_text(query_coord),
        reference_coordinates,
        query_alignment_regions,
        reference_alignment_regions,
        cigar,
    )), read_alternative_tag(fields))


def _build_row(
    pair: core.PairRow,
    gcigars,
    graph_sequences,
    graph_path_coords,
    chrom_lengths,
    graph_mappings,
    ref_reader,
    coord_by_name,
    allowchroms,
    assembly_sequences,
    query_cache,
    reference_flank_rescue=DEFAULT_REFERENCE_FLANK_RESCUE,
) -> Optional[str]:
    record_sequences = _query_record_sequences(
        pair, gcigars, assembly_sequences, query_cache
    )
    row_text = core.build_provisional_row(
        pair,
        gcigars,
        graph_sequences,
        record_sequences,
        graph_path_coords,
        chrom_lengths,
        graph_mappings,
        ref_reader,
        coord_by_name,
        allowchroms=allowchroms,
    )
    return format_provisional_row(
        row_text,
        pair,
        record_sequences.get(pair.query_name or ""),
        ref_reader,
        reference_flank_rescue,
    )


def _build_comparison_rows(
    comparison: Tuple[core.PairRow, Sequence[core.PairRow]],
    gcigars,
    graph_sequences,
    graph_path_coords,
    chrom_lengths,
    graph_mappings,
    ref_reader,
    coord_by_name,
    allowchroms,
    assembly_sequences,
    query_cache,
    reference_flank_rescue=DEFAULT_REFERENCE_FLANK_RESCUE,
) -> List[str]:
    """Build one comparison and retain a grouped CIGAR as one physical row."""
    base_pair, output_pairs = comparison
    record_sequences = _query_record_sequences(
        base_pair, gcigars, assembly_sequences, query_cache,
    )
    row_text = core.build_provisional_row(
        base_pair,
        gcigars,
        graph_sequences,
        record_sequences,
        graph_path_coords,
        chrom_lengths,
        graph_mappings,
        ref_reader,
        coord_by_name,
        allowchroms=allowchroms,
    )
    if is_grouped_pair(base_pair) or len(output_pairs) > 1:
        formatted = format_merged_comparison_row(
            row_text,
            base_pair,
            output_pairs,
            record_sequences.get(base_pair.query_name or ""),
            ref_reader,
            reference_flank_rescue,
        )
    else:
        formatted = format_provisional_row(
            row_text,
            output_pairs[0],
            record_sequences.get(base_pair.query_name or ""),
            ref_reader,
            reference_flank_rescue,
        )
    if formatted is not None:
        formatted = core.append_tag(
            formatted, getattr(base_pair, 'alternative_intervals', None),
        )
    return [formatted] if formatted is not None else []


def _init_worker_state(
    gcigars,
    graph_sequences,
    graph_path_coords,
    chrom_lengths,
    graph_mappings,
    coord_by_name,
    assembly_sequences,
    allowchroms,
    reference_path: str,
    local_reference_templates: str,
    realignment_enabled: bool,
) -> None:
    global _WORKER_GCIGARS
    global _WORKER_GRAPH_SEQUENCES
    global _WORKER_GRAPH_PATH_COORDS
    global _WORKER_CHROM_LENGTHS
    global _WORKER_GRAPH_MAPPINGS
    global _WORKER_COORD_BY_NAME
    global _WORKER_ALLOWCHROMS
    global _WORKER_ASSEMBLY_SEQUENCES
    global _WORKER_REF_READER
    global _WORKER_QUERY_CACHE
    core.set_sv_realignment_enabled(realignment_enabled)
    _WORKER_GCIGARS = gcigars
    _WORKER_GRAPH_SEQUENCES = graph_sequences
    _WORKER_GRAPH_PATH_COORDS = graph_path_coords
    _WORKER_CHROM_LENGTHS = chrom_lengths
    _WORKER_GRAPH_MAPPINGS = graph_mappings
    _WORKER_COORD_BY_NAME = coord_by_name
    _WORKER_ASSEMBLY_SEQUENCES = assembly_sequences
    _WORKER_ALLOWCHROMS = set(allowchroms or [])
    if _SHARED_REFERENCE_READER is not None:
        # Fork start method: adopt the parent's in-RAM reference (sequence
        # pages stay shared copy-on-write).  The template reader keeps a
        # seekable file handle, which forked siblings must not share, so
        # reopen it in this worker only.
        _WORKER_REF_READER = _SHARED_REFERENCE_READER
        if _WORKER_REF_READER._template_reader is not None:
            _WORKER_REF_READER._template_reader = core.FastaRegionReader(
                local_reference_templates,
            )
    else:
        _WORKER_REF_READER = core.FastaRegionReader(
            reference_path,
            template_fasta=(local_reference_templates or None),
        )
    _WORKER_QUERY_CACHE = {}


def _worker_build_row(
    task: Tuple[int, core.PairRow]
) -> Tuple[int, Optional[str], Optional[str]]:
    index, pair = task
    try:
        row = _build_row(
            pair,
            _WORKER_GCIGARS,
            _WORKER_GRAPH_SEQUENCES,
            _WORKER_GRAPH_PATH_COORDS,
            _WORKER_CHROM_LENGTHS,
            _WORKER_GRAPH_MAPPINGS,
            _WORKER_REF_READER,
            _WORKER_COORD_BY_NAME,
            _WORKER_ALLOWCHROMS,
            _WORKER_ASSEMBLY_SEQUENCES,
            _WORKER_QUERY_CACHE,
        )
        return index, row, None
    except Exception:
        label = pair.query_name or pair.label or f"line {pair.line_no}"
        return (
            index,
            None,
            f"failed on pair index={index} line={pair.line_no} "
            f"label={label}:\n{traceback.format_exc()}",
        )


def _worker_build_comparison(
    task: Tuple[int, Tuple[core.PairRow, Sequence[core.PairRow]]]
) -> Tuple[int, Optional[List[str]], Optional[str]]:
    index, comparison = task
    base_pair, _output_pairs = comparison
    try:
        rows = _build_comparison_rows(
            comparison,
            _WORKER_GCIGARS,
            _WORKER_GRAPH_SEQUENCES,
            _WORKER_GRAPH_PATH_COORDS,
            _WORKER_CHROM_LENGTHS,
            _WORKER_GRAPH_MAPPINGS,
            _WORKER_REF_READER,
            _WORKER_COORD_BY_NAME,
            _WORKER_ALLOWCHROMS,
            _WORKER_ASSEMBLY_SEQUENCES,
            _WORKER_QUERY_CACHE,
        )
        return index, rows, None
    except Exception:
        label = base_pair.query_name or base_pair.label or f"line {base_pair.line_no}"
        return (
            index,
            None,
            f"failed on comparison index={index} line={base_pair.line_no} "
            f"label={label}:\n{traceback.format_exc()}",
        )


def group_comparison_pairs(
    pairs: Sequence[core.PairRow],
) -> List[Tuple[core.PairRow, Tuple[core.PairRow, ...]]]:
    """Collapse repeated per-part views into one expensive comparison task."""
    grouped: Dict[Tuple[object, ...], List[core.PairRow]] = {}
    order: List[Tuple[object, ...]] = []
    for pair in pairs:
        key = (
            pair.query_name,
            pair.ref_name,
            pair.query_coord,
            pair.ref_coord,
        )
        if key not in grouped:
            grouped[key] = []
            order.append(key)
        grouped[key].append(pair)
    return [
        (grouped[key][0], tuple(grouped[key]))
        for key in order
    ]


def group_column10_candidates(
    pairs: Sequence[core.PairRow],
) -> List[Tuple[Tuple[core.PairRow, Tuple[core.PairRow, ...]], ...]]:
    """Group alternative column-10 reference intersections per query.

    Ordinary callers without column-10 candidate metadata retain one
    comparison per group, preserving the legacy behavior.
    """
    comparisons = group_comparison_pairs(pairs)
    grouped: Dict[Tuple[object, ...], List[Tuple[core.PairRow, Tuple[core.PairRow, ...]]]] = {}
    order: List[Tuple[object, ...]] = []
    for index, comparison in enumerate(comparisons):
        base_pair = comparison[0]
        key = getattr(
            base_pair,
            "column10_candidate_key",
            ("legacy_comparison", index),
        )
        if key not in grouped:
            grouped[key] = []
            order.append(key)
        grouped[key].append(comparison)
    return [tuple(grouped[key]) for key in order]


def graph_cigar_aligned_overlap(row_text: str) -> int:
    """Count bases consumed on both axes in a seven-column output row."""
    fields = row_text.split("\t")
    if len(fields) < 7:
        raise ValueError("formatted graph-CIGAR row has fewer than seven columns")
    return sum(
        operation.n
        for segment in core.parse_graphic_segments(
            fields[6], "column-10 candidate overlap",
        )
        for operation in segment.ops
        if operation.op in {"=", "M", "X"}
    )


def graph_cigar_aligned_unmasked_query_bases(
    row_text: str,
    assembly_sequences: Mapping[str, str],
) -> int:
    """Count uppercase A/C/G/T query bases in aligned CIGAR operations.

    Only query bases paired to a graph/reference base (``=``, ``M``, or ``X``)
    contribute. Insertions and soft/hard-clipped query sequence do not satisfy
    the 2,000-bp last-resort placement threshold.
    """
    fields = row_text.split("\t")
    if len(fields) < 7:
        raise ValueError("formatted graph-CIGAR row has fewer than seven columns")
    query_coord = core.parse_coord(fields[2])
    if query_coord is None:
        raise ValueError(
            f"formatted graph-CIGAR row has malformed query coordinate "
            f"{fields[2]!r}"
        )
    contig_sequence = assembly_sequences.get(query_coord.chrom)
    if contig_sequence is None:
        contig_sequence = assembly_sequences.get(
            core.strip_match_name(query_coord.chrom)
        )
    if contig_sequence is None:
        raise KeyError(
            f"query contig {query_coord.chrom!r} is absent from the loaded "
            "assembly while counting aligned unmasked bases"
        )
    if (
        query_coord.start < 0
        or query_coord.end < query_coord.start
        or query_coord.end > len(contig_sequence)
    ):
        raise ValueError(
            f"query coordinate {fields[2]!r} is outside contig "
            f"{query_coord.chrom!r} length {len(contig_sequence)}"
        )
    query_sequence = contig_sequence[query_coord.start:query_coord.end]
    if query_coord.strand == "-":
        query_sequence = core.revcomp(query_sequence)

    query_position = 0
    unmasked = 0
    for segment in core.parse_graphic_segments(
        fields[6], "last-resort unmasked aligned query bases",
    ):
        for operation in segment.ops:
            consumed = core.query_consume(operation.op, operation.n)
            if operation.op in {"=", "M", "X"}:
                unmasked += count_unmasked(
                    query_sequence[
                        query_position:query_position + operation.n
                    ]
                )
            query_position += consumed
    if query_position != len(query_sequence):
        raise ValueError(
            "formatted graph-CIGAR query span disagrees with its assembly "
            f"coordinate: {query_position}/{len(query_sequence)}"
        )
    return unmasked


def _build_column10_candidate_rows(
    candidate_group: Sequence[
        Tuple[core.PairRow, Sequence[core.PairRow]]
    ],
    gcigars,
    graph_sequences,
    graph_path_coords,
    chrom_lengths,
    graph_mappings,
    ref_reader,
    coord_by_name,
    allowchroms,
    assembly_sequences,
    query_cache,
) -> List[str]:
    """Align every reference intersection and retain the best overlap."""
    rows, _stats = _evaluate_column10_candidate_rows(
        candidate_group,
        gcigars,
        graph_sequences,
        graph_path_coords,
        chrom_lengths,
        graph_mappings,
        ref_reader,
        coord_by_name,
        allowchroms,
        assembly_sequences,
        query_cache,
    )
    return rows


def _evaluate_column10_candidate_rows(
    candidate_group: Sequence[
        Tuple[core.PairRow, Sequence[core.PairRow]]
    ],
    gcigars,
    graph_sequences,
    graph_path_coords,
    chrom_lengths,
    graph_mappings,
    ref_reader,
    coord_by_name,
    allowchroms,
    assembly_sequences,
    query_cache,
) -> Tuple[List[str], Dict[str, int]]:
    """Try candidate tiers in order and retain the best usable overlap.

    Priority 0 accepts column-2 similarity above 90%. Priority 1 accepts any
    usable alignment to one unmerged column-9 reference. Priority 2 reuses the
    cached column-2 result above 50%. Priority 3 accepts any usable general
    column-10 coordinate result. Priority 4 is reached only when all of those
    fail and reuses column 2 if it has at least 2,000 aligned unmasked query
    bases.
    """
    selected = None
    errors = []
    unusable = 0
    usable = 0
    low_coverage = 0
    fallback_used = 0
    attempted = 0
    candidates_by_priority: Dict[
        int, List[Tuple[core.PairRow, Sequence[core.PairRow]]]
    ] = {}
    for comparison in candidate_group:
        base_pair = comparison[0]
        priority = int(getattr(base_pair, "reference_candidate_priority", 0))
        candidates_by_priority.setdefault(priority, []).append(comparison)
        retry_priority = getattr(
            base_pair, "reference_candidate_retry_priority", None,
        )
        if retry_priority is not None:
            candidates_by_priority.setdefault(
                int(retry_priority), [],
            ).append(comparison)
        last_retry_priority = getattr(
            base_pair, "reference_candidate_last_retry_priority", None,
        )
        if last_retry_priority is not None:
            candidates_by_priority.setdefault(
                int(last_retry_priority), [],
            ).append(comparison)

    evaluation_cache = {}
    for priority in sorted(candidates_by_priority):
        tier_successful = []
        for comparison in candidates_by_priority[priority]:
            base_pair = comparison[0]
            cache_key = id(base_pair)
            cached = evaluation_cache.get(cache_key)
            new_evaluation = cached is None
            if new_evaluation:
                attempted += 1
                try:
                    rows = _build_comparison_rows(
                        comparison,
                        gcigars,
                        graph_sequences,
                        graph_path_coords,
                        chrom_lengths,
                        graph_mappings,
                        ref_reader,
                        coord_by_name,
                        allowchroms,
                        assembly_sequences,
                        query_cache,
                    )
                    cached = (rows, None)
                except Exception:
                    error_text = traceback.format_exc()
                    errors.append(error_text)
                    cached = (None, error_text)
                evaluation_cache[cache_key] = cached
            rows, error_text = cached
            if error_text is not None:
                continue
            if not rows:
                # Formatting returns no row when the comparison has no usable
                # sufficient bidirectional alignment-confidence score.
                if new_evaluation:
                    unusable += 1
                continue
            span = int(getattr(base_pair, "column10_intersection_span", 0))
            reference_sort = tuple(getattr(
                base_pair,
                "column10_reference_sort",
                (base_pair.ref_coord_text, base_pair.ref_name or ""),
            ))
            for row in rows:
                aligned_overlap = graph_cigar_aligned_overlap(row)
                aligned_unmasked = (
                    graph_cigar_aligned_unmasked_query_bases(
                        row, assembly_sequences,
                    )
                    if priority == 4 else -1
                )
                query_span = (
                    base_pair.query_coord.end - base_pair.query_coord.start
                    if base_pair.query_coord is not None else 0
                )
                aligned_fraction = (
                    aligned_overlap / query_span if query_span > 0 else 0.0
                )
                tier_successful.append((
                    -aligned_overlap,
                    -span,
                    reference_sort,
                    row,
                    aligned_fraction,
                    aligned_unmasked,
                ))
            if new_evaluation:
                usable += len(rows)
        if not tier_successful:
            continue
        tier_successful.sort(key=lambda value: value[:3])

        if priority == 0:
            passing_similarity = [
                candidate for candidate in tier_successful
                if candidate[4] > HIGH_SIMILARITY_QUERY_COVERAGE
            ]
            if passing_similarity:
                passing_similarity.sort(key=lambda value: value[:3])
                selected = passing_similarity[0]
                fallback_used = 1
                break
            continue

        if priority in {1, 3}:
            # Single-column-9 and column-10 tiers deliberately accept partial
            # alignments so deletions and substitutions remain visible.
            selected = tier_successful[0]
            break

        # Priority 2 reuses the already computed column-2 alignment.
        if priority == 4:
            passing_similarity = [
                candidate for candidate in tier_successful
                if candidate[5]
                >= MIN_LAST_RESORT_SIMILARITY_ALIGNED_BASES
            ]
            if passing_similarity:
                passing_similarity.sort(key=lambda value: value[:3])
                selected = passing_similarity[0]
                fallback_used = 1
                break
            continue

        passing_similarity = [
            candidate for candidate in tier_successful
            if candidate[4] > MIN_REFERENCE_CANDIDATE_QUERY_COVERAGE
        ]
        low_coverage += len(tier_successful) - len(passing_similarity)
        if passing_similarity:
            passing_similarity.sort(key=lambda value: value[:3])
            selected = passing_similarity[0]
            fallback_used = 1
            break

    if selected is None:
        if errors and len(errors) == attempted:
            raise RuntimeError(
                "all genomic and same-graph similarity candidates failed:\n"
                + "\n".join(errors)
            )
        return [], {
            "candidate_comparisons": attempted,
            "usable_candidates": usable,
            "unusable_candidates": unusable,
            "low_coverage_candidates": low_coverage,
            "bestref_similarity_groups": 0,
            "candidate_errors": len(errors),
            "emitted_groups": 0,
            "empty_groups": 1,
        }
    return [selected[3]], {
        "candidate_comparisons": attempted,
        "usable_candidates": usable,
        "unusable_candidates": unusable,
        "low_coverage_candidates": low_coverage,
        "bestref_similarity_groups": fallback_used,
        "candidate_errors": len(errors),
        "emitted_groups": 1,
        "empty_groups": 0,
    }


def _worker_build_column10_candidates(
    task: Tuple[
        int,
        Sequence[Tuple[core.PairRow, Sequence[core.PairRow]]],
    ],
) -> Tuple[int, Optional[List[str]], Optional[str]]:
    index, candidate_group = task
    base_pair = candidate_group[0][0]
    try:
        rows, stats = _evaluate_column10_candidate_rows(
            candidate_group,
            _WORKER_GCIGARS,
            _WORKER_GRAPH_SEQUENCES,
            _WORKER_GRAPH_PATH_COORDS,
            _WORKER_CHROM_LENGTHS,
            _WORKER_GRAPH_MAPPINGS,
            _WORKER_REF_READER,
            _WORKER_COORD_BY_NAME,
            _WORKER_ALLOWCHROMS,
            _WORKER_ASSEMBLY_SEQUENCES,
            _WORKER_QUERY_CACHE,
        )
        return index, rows, stats, None
    except Exception:
        label = base_pair.query_name or base_pair.label or f"line {base_pair.line_no}"
        return (
            index,
            None,
            None,
            f"failed on column-10 candidate group index={index} "
            f"line={base_pair.line_no} label={label}:\n{traceback.format_exc()}",
        )


def is_grouped_pair(pair: core.PairRow) -> bool:
    return (
        bool(pair.query_name and ";" in pair.query_name)
        or bool(pair.ref_name and ";" in pair.ref_name)
    )


def output_pair_finished_key(pair: core.PairRow) -> Tuple[str, str, str]:
    """Resume key matching the individual seven-column output row."""
    output_name = getattr(pair, "output_query_name", None)
    if output_name:
        return (
            output_name,
            getattr(pair, "output_ref_name", pair.ref_name or ""),
            getattr(pair, "output_query_alignment", pair.label),
        )
    return core.finished_pair_key(pair)


def comparison_finished_key(
    comparison: Tuple[core.PairRow, Sequence[core.PairRow]],
) -> Tuple[str, str, str]:
    """Return the exact key written for one physical comparison row.

    Simple comparisons still write an individual row.  A grouped comparison
    writes one composite query/reference row whose fifth column contains the
    ordered ownership-alignment interval for every query label.  Resume must
    therefore keep or recompute that physical row atomically; filtering its
    component PairRows independently can silently recreate a partial group.
    """
    base_pair, output_pairs = comparison
    if is_grouped_pair(base_pair) or len(output_pairs) > 1:
        return (
            base_pair.query_name or "",
            getattr(base_pair, "output_ref_name", base_pair.ref_name or ""),
            ";".join(
                getattr(
                    pair,
                    "output_query_alignment",
                    base_pair.query_coord_text,
                )
                for pair in output_pairs
            ),
        )
    return output_pair_finished_key(output_pairs[0])


def _coord_by_record_name(pairs: Iterable[core.PairRow]) -> Dict[str, core.Coord]:
    result: Dict[str, core.Coord] = {}
    for pair in pairs:
        for name, coord in (
            (pair.query_name, pair.query_coord),
            (pair.ref_name, pair.ref_coord),
        ):
            if name is None or coord is None:
                continue
            for alias in dict.fromkeys((name, core.strip_match_name(name))):
                old = result.get(alias)
                if old is not None and old != coord:
                    raise ValueError(
                        f"record {alias!r} has conflicting GenomeLift coordinates"
                    )
                result[alias] = coord
    return result


def _validate_direct_pairs(
    pairs: Sequence[core.PairRow], gcigars: Dict[str, str]
) -> None:
    for pair in pairs:
        if (
            pair.query_name is None
            or pair.ref_name is None
            or pair.query_coord is None
            or pair.ref_coord is None
        ):
            raise ValueError(
                f"internal pair {pair.line_no} is not resolved"
            )
        missing = [
            name
            for name in (pair.query_name, pair.ref_name)
            if name not in gcigars
        ]
        if missing:
            raise KeyError(
                f"internal pair {pair.line_no}: missing prepared graph CIGAR(s): "
                + ", ".join(repr(value) for value in missing)
            )


class _CompletedRowBuffer:
    """Batch completed, potentially out-of-order worker rows for output.

    A zero buffer size is the diagnostic mode: every completed row is written
    and flushed immediately.  This leaves an exact on-disk record of completed
    comparisons when another worker is stuck on a pathological pair.
    """

    def __init__(
        self,
        output,
        buffer_size: int,
        max_bytes: int = 64 * 1024 * 1024,
    ) -> None:
        self.output = output
        self.buffer_size = int(buffer_size)
        self.max_bytes = int(max_bytes)
        self.rows: List[str] = []
        self.bytes = 0

    def write(self, row: str) -> None:
        if self.buffer_size == 0:
            self.output.write(row + "\n")
            self.output.flush()
            return
        self.rows.append(row)
        self.bytes += len(row) + 1
        if (
            len(self.rows) >= self.buffer_size
            or (self.max_bytes > 0 and self.bytes >= self.max_bytes)
        ):
            self.flush()

    def flush(self) -> None:
        if self.rows:
            self.output.write("\n".join(self.rows) + "\n")
            self.rows.clear()
            self.bytes = 0
        self.output.flush()


CONVERSION_STAT_KEYS = (
    "candidate_comparisons",
    "usable_candidates",
    "unusable_candidates",
    "low_coverage_candidates",
    "bestref_similarity_groups",
    "candidate_errors",
    "emitted_groups",
    "empty_groups",
)


def _empty_conversion_stats() -> Dict[str, int]:
    return {key: 0 for key in CONVERSION_STAT_KEYS}


def _add_conversion_stats(total: Dict[str, int], values) -> None:
    if not values:
        return
    for key in CONVERSION_STAT_KEYS:
        total[key] += int(values.get(key, 0))


def _row_query_interval(row_text: str) -> Optional[core.Coord]:
    fields = row_text.split("\t", 3)
    if len(fields) < 3:
        return None
    return core.parse_coord(fields[2])


def _merged_interval_map(
    coordinates: Iterable[core.Coord],
) -> Dict[str, List[Tuple[int, int]]]:
    by_contig: Dict[str, List[Tuple[int, int]]] = {}
    seen = set()
    for coord in coordinates:
        if coord is None or coord.end <= coord.start:
            continue
        key = (coord.chrom, coord.start, coord.end)
        if key in seen:
            continue
        seen.add(key)
        by_contig.setdefault(coord.chrom, []).append((coord.start, coord.end))
    merged: Dict[str, List[Tuple[int, int]]] = {}
    for contig, intervals in by_contig.items():
        intervals.sort()
        output = []
        start, end = intervals[0]
        for next_start, next_end in intervals[1:]:
            if next_start <= end:
                end = max(end, next_end)
            else:
                output.append((start, end))
                start, end = next_start, next_end
        output.append((start, end))
        merged[contig] = output
    return merged


def _interval_map_size(intervals: Mapping[str, Sequence[Tuple[int, int]]]) -> int:
    return sum(
        end - start
        for values in intervals.values()
        for start, end in values
    )


def _interval_map_intersection_size(
    left: Mapping[str, Sequence[Tuple[int, int]]],
    right: Mapping[str, Sequence[Tuple[int, int]]],
) -> int:
    total = 0
    for contig, left_values in left.items():
        right_values = right.get(contig, ())
        left_index = 0
        right_index = 0
        while left_index < len(left_values) and right_index < len(right_values):
            left_start, left_end = left_values[left_index]
            right_start, right_end = right_values[right_index]
            total += max(
                0, min(left_end, right_end) - max(left_start, right_start),
            )
            if left_end <= right_end:
                left_index += 1
            else:
                right_index += 1
    return total


def _read_output_query_intervals(path: str) -> List[core.Coord]:
    intervals = []
    with open(path, "rt") as handle:
        for raw in handle:
            if not raw.strip() or raw.startswith("#"):
                continue
            coord = _row_query_interval(raw.rstrip("\r\n"))
            if coord is not None and coord.end > coord.start:
                intervals.append(coord)
    return intervals


def _log_conversion_coverage(
    eligible_coordinates: Iterable[core.Coord],
    emitted_coordinates: Iterable[core.Coord],
    stats: Mapping[str, int],
) -> None:
    eligible = _merged_interval_map(eligible_coordinates)
    emitted = _merged_interval_map(emitted_coordinates)
    eligible_size = _interval_map_size(eligible)
    emitted_size = _interval_map_size(emitted)
    shared_size = _interval_map_intersection_size(eligible, emitted)
    uncovered = eligible_size - shared_size
    outside = emitted_size - shared_size
    percentage = (100.0 * uncovered / eligible_size) if eligible_size else 0.0
    sys.stderr.write(
        "[graphcigartoref_persample] coverage summary: "
        f"eligible core union={eligible_size} bp "
        f"({eligible_size / 1_000_000_000:.3f} Gb); "
        f"emitted query union={emitted_size} bp "
        f"({emitted_size / 1_000_000_000:.3f} Gb); "
        f"eligible bases not emitted={uncovered} bp ({percentage:.2f}%); "
        f"emitted context outside eligible cores={outside} bp\n"
    )
    sys.stderr.write(
        "[graphcigartoref_persample] candidate summary for comparisons run: "
        f"groups emitted={stats.get('emitted_groups', 0)}, "
        f"groups with no usable row={stats.get('empty_groups', 0)}, "
        f"reference candidates attempted={stats.get('candidate_comparisons', 0)}, "
        f"usable={stats.get('usable_candidates', 0)}, "
        f"rejected without a usable score-clipped row="
        f"{stats.get('unusable_candidates', 0)}, "
        f"usable but below the 50% decision threshold="
        f"{stats.get('low_coverage_candidates', 0)}, "
        f"groups selected by column-2 similarity="
        f"{stats.get('bestref_similarity_groups', 0)}, "
        f"candidate exceptions={stats.get('candidate_errors', 0)}\n"
    )


def _write_parallel(
    pairs: Sequence[core.PairRow],
    output,
    args: argparse.Namespace,
    shared_values: Tuple,
) -> Tuple[Dict[str, int], List[core.Coord]]:
    context = mp.get_context(args.mp_start_method)
    buffer_size = int(getattr(args, "buffer_size", 1000))
    buffer_bytes = int(getattr(args, "buffer_bytes", 64 * 1024 * 1024))
    progress_every = int(getattr(args, "progress_every", 1000))
    progress_seconds = float(getattr(args, "progress_seconds", 60.0))
    maxtasksperchild = int(getattr(args, "maxtasksperchild", 512))
    candidate_groups = group_column10_candidates(pairs)
    comparison_count = sum(len(group) for group in candidate_groups)
    conversion_stats = _empty_conversion_stats()
    emitted_coordinates: List[core.Coord] = []
    sys.stderr.write(
        "[graphcigartoref_persample] converting "
        f"{comparison_count} reference-intersection comparison(s) in "
        f"{len(candidate_groups)} query candidate group(s) with "
        f"{args.processes} worker(s); "
        f"row buffer={buffer_size}, byte cap={buffer_bytes or 'off'}, "
        f"worker recycle={maxtasksperchild or 'off'}\n"
    )

    def run_pool(values: Tuple) -> None:
        initargs = values + (
            assembly_sequences,
            tuple(sorted(allowchroms)),
            args.reference,
            getattr(args, "local_reference_templates", ""),
            bool(getattr(args, "realignment", False)),
        )
        writer = _CompletedRowBuffer(output, buffer_size, buffer_bytes)
        completed = 0
        emitted = 0
        started = time.monotonic()
        # Move the fully built shared structures (assembly, graph sequences,
        # CIGAR tables) into the GC permanent generation before forking.
        # Cyclic-GC passes write to object headers, which would otherwise
        # unshare copy-on-write pages in every worker; freezing keeps those
        # pages shared.  Collection behavior for new objects is unchanged.
        gc.freeze()
        try:
            with context.Pool(
                processes=args.processes,
                initializer=_init_worker_state,
                initargs=initargs,
                maxtasksperchild=(maxtasksperchild or None),
            ) as pool:
                # Conversion costs vary by orders of magnitude. Dispatch one
                # merged comparison per pool task and consume results as they
                # finish. A multi-match produces one physical row containing
                # the complete CIGAR plus all logical ownership intervals.
                result_iter = pool.imap_unordered(
                    _worker_build_column10_candidates,
                    enumerate(candidate_groups),
                    1,
                )
                while completed < len(candidate_groups):
                    try:
                        _index, rows, stats, error = result_iter.next(
                            timeout=(progress_seconds if progress_seconds > 0 else None)
                        )
                    except mp.TimeoutError:
                        elapsed = max(1e-6, time.monotonic() - started)
                        sys.stderr.write(
                            "[graphcigartoref_persample] still waiting: "
                            f"{completed}/{len(candidate_groups)} candidate "
                            "groups completed "
                            f"after {elapsed:.0f}s; the next result may be a "
                            "pathological long comparison\n"
                        )
                        continue
                    completed += 1
                    if error is not None:
                        base_pair = candidate_groups[_index][0][0]
                        if is_grouped_pair(base_pair):
                            _add_conversion_stats(conversion_stats, {
                                "candidate_comparisons": len(
                                    candidate_groups[_index]
                                ),
                                "candidate_errors": len(
                                    candidate_groups[_index]
                                ),
                                "empty_groups": 1,
                            })
                            sys.stderr.write(
                                "[graphcigartoref_persample] warning: skipping "
                                f"unconvertible grouped comparison:\n{error}\n"
                            )
                            continue
                        pool.terminate()
                        raise RuntimeError(error)
                    _add_conversion_stats(conversion_stats, stats)
                    for row in rows or ():
                        writer.write(row)
                        coord = _row_query_interval(row)
                        if coord is not None:
                            emitted_coordinates.append(coord)
                        emitted += 1
                    if (
                        progress_every > 0
                        and (
                            completed % progress_every == 0
                            or completed == len(candidate_groups)
                        )
                    ):
                        elapsed = max(1e-6, time.monotonic() - started)
                        rate = completed / elapsed
                        last_pair = candidate_groups[_index][0][0]
                        last_label = last_pair.query_name or last_pair.label
                        if len(last_label) > 240:
                            last_label = last_label[:237] + "..."
                        sys.stderr.write(
                            "[graphcigartoref_persample] completed "
                            f"{completed}/{len(candidate_groups)} candidate groups, "
                            f"{emitted} output rows ({rate:.2f}/s); "
                            f"last={last_label}\n"
                        )
        finally:
            # Preserve every result already returned to the parent on normal
            # completion, Python exceptions, or an interactive interruption.
            writer.flush()

    values = shared_values[:6]
    assembly_sequences = shared_values[6]
    allowchroms = shared_values[7]
    if args.share_mode == "manager":
        with context.Manager() as manager:
            managed = tuple(manager.dict(value) for value in values)
            run_pool(managed)
    else:
        run_pool(values)
    return conversion_stats, emitted_coordinates


def parse_args(argv: Optional[Sequence[str]] = None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=(
            "Direct per-sample graphcigartoref from GenomeLift, whole query/reference "
            "alignments, a faidxed sample assembly, and the local graph list"
        )
    )
    parser.add_argument(
        "--genomelift", required=True,
        help="GenomeLift TSV for the sample and reference alleles",
    )
    parser.add_argument(
        "--align-query", required=True,
        help="sample align_partition_hotspots alignment file",
    )
    parser.add_argument(
        "--align-ref", default="",
        help="reference alignment file; inferred from --align-query when omitted",
    )
    parser.add_argument(
        "-q",
        "--fasta-query",
        required=True,
        help="sample assembly FASTA loaded entirely into RAM; FASTA.fai must exist",
    )
    parser.add_argument(
        "-G",
        "--graph-folder",
        required=True,
        help="root containing GRAPH_NAME/GRAPH_NAME.FA local graphs",
    )
    parser.add_argument(
        "--graph-summary",
        default="",
        help=(
            "optional graph-build summary directory containing "
            "local_graphs.tsv and alternatives.fasta; required paths are "
            "loaded from this consolidated package before falling back to "
            "individual GRAPH_NAME.FA files"
        ),
    )
    parser.add_argument(
        "-L",
        "--graph-list",
        required=True,
        help="ordered graph-name/path list used by align_partition_hotspots.py",
    )
    parser.add_argument(
        "-r", "--reference", required=True, help="faidxed reference FASTA"
    )
    parser.add_argument(
        "--local-reference-templates",
        default="",
        help=(
            "optional sparse local-template FASTA produced by "
            "local_reference_templates.py for graphs lacking the requested "
            "reference haplotype"
        ),
    )
    parser.add_argument(
        "--query-genome", default="",
        help="query sample/haplotype name; inferred from --align-query when omitted",
    )
    parser.add_argument("--refhaplo", default="CHM13_h1")
    parser.add_argument("-o", "--output", default="")
    parser.add_argument(
        "-u",
        "--unreplace",
        action="store_true",
        help="clean/resume an existing --output and append unfinished pairs",
    )
    parser.add_argument(
        "-t",
        "--processes",
        type=int,
        default=1,
        help="number of worker processes; 1 disables multiprocessing",
    )
    parser.add_argument(
        "--chunksize",
        type=int,
        default=0,
        help=(
            "record-preparation Pool.imap chunksize; 0 chooses automatically. "
            "CIGAR conversion always dispatches one pair per worker task"
        ),
    )
    parser.add_argument(
        "--buffer-size",
        type=int,
        default=1000,
        metavar="ROWS",
        help=(
            "number of completed parallel result rows buffered in RAM before "
            "writing (default: 1000); 0 writes and flushes every completed "
            "row immediately for debugging"
        ),
    )
    parser.add_argument(
        "--buffer-bytes",
        type=int,
        default=64 * 1024 * 1024,
        metavar="BYTES",
        help=(
            "hard upper bound for completed output rows held in RAM before a "
            "flush (default: 67108864; 0 disables the byte bound)"
        ),
    )
    parser.add_argument(
        "--progress-every",
        type=int,
        default=1000,
        metavar="N",
        help=(
            "report completed comparison progress every N results; 0 disables "
            "progress messages (default: 1000)"
        ),
    )
    parser.add_argument(
        "--progress-seconds",
        type=float,
        default=60.0,
        metavar="SECONDS",
        help=(
            "emit a waiting diagnostic when no comparison finishes within this "
            "many seconds; 0 disables the watchdog (default: 60)"
        ),
    )
    parser.add_argument(
        "--maxtasksperchild",
        type=int,
        default=512,
        metavar="N",
        help=(
            "recycle each worker after N comparison tasks to bound allocator/cache "
            "growth; 0 keeps workers for the whole run (default: 512)"
        ),
    )
    parser.add_argument(
        "--mp-start-method",
        choices=("fork", "spawn", "forkserver"),
        default="fork",
    )
    parser.add_argument(
        "--share-mode",
        choices=("fork", "manager"),
        default="fork",
        help="share graph dictionaries copy-on-write with fork, or through a Manager",
    )
    parser.add_argument(
        "--add-alt",
        action="store_true",
        help="keep alignments to alternate/non-reference chromosomes",
    )
    realignment_group = parser.add_mutually_exclusive_group()
    realignment_group.add_argument(
        "--realignment",
        dest="realignment",
        action="store_true",
        default=True,
        help=(
            "enable optional local and score-linked SV-region realignment; "
            "enabled by default"
        ),
    )
    realignment_group.add_argument(
        "--no-realignment",
        dest="realignment",
        action="store_false",
        help="disable local and score-linked SV-region realignment",
    )
    parser.add_argument(
        "--max-extension",
        type=int,
        default=10000,
        metavar="BP",
        help=(
            "maximum query sequence context taken from either side according "
            "to GenomeLift columns 13/14 (default: 10000)"
        ),
    )
    parser.add_argument(
        "--reference-extension",
        type=int,
        default=DEFAULT_REFERENCE_EXTENSION,
        metavar="BP",
        help=(
            "legacy compatibility option; reference tiers now seed a complete "
            "same-graph reference _align.txt row, so this fixed-distance "
            "extension is not used for those candidates"
        ),
    )
    parser.add_argument(
        "--cross-validation-extension",
        "--cross-graph-extension",
        type=int,
        default=4000,
        metavar="BP",
        help=(
            "additional symmetric context around each GenomeLift block for "
            "cross-graph validation after ownership extensions "
            "(default: 4000; use 0 to disable)"
        ),
    )
    parser.add_argument(
        '--alternative', default='',
        help='four-column original-template BED; carry projected calling intervals to graphreftovcf',
    )
    return parser.parse_args(argv)


def run(args: argparse.Namespace) -> int:
    alternative_path = getattr(args, 'alternative', '')
    alternative_policy = core.AlternativeIntervals(
        alternative_path,
        read_templates(args.local_reference_templates) if args.local_reference_templates else (),
    ) if alternative_path else None
    if alternative_policy and args.unreplace and args.output and os.path.exists(args.output):
        raise ValueError('--alternative requires a fresh output; omit --unreplace to avoid retaining rows from an older interval policy')
    core.set_sv_realignment_enabled(bool(getattr(args, "realignment", False)))
    if args.processes < 1:
        raise ValueError("--processes must be >= 1")
    if args.chunksize < 0:
        raise ValueError("--chunksize must be >= 0")
    if args.buffer_size < 0:
        raise ValueError("--buffer-size must be >= 0")
    if getattr(args, "buffer_bytes", 64 * 1024 * 1024) < 0:
        raise ValueError("--buffer-bytes must be >= 0")
    if getattr(args, "progress_every", 1000) < 0:
        raise ValueError("--progress-every must be >= 0")
    if getattr(args, "progress_seconds", 60.0) < 0:
        raise ValueError("--progress-seconds must be >= 0")
    if getattr(args, "maxtasksperchild", 512) < 0:
        raise ValueError("--maxtasksperchild must be >= 0")
    if args.max_extension < 0:
        raise ValueError("--max-extension must be >= 0")
    if getattr(args, "reference_extension", DEFAULT_REFERENCE_EXTENSION) < 0:
        raise ValueError("--reference-extension must be >= 0")
    if getattr(args, "cross_validation_extension", 4000) < 0:
        raise ValueError("--cross-validation-extension must be >= 0")
    if (
        args.processes > 1
        and args.mp_start_method != "fork"
    ):
        sys.stderr.write(
            "[graphcigartoref_persample] warning: a non-fork start method copies "
            "the in-memory sample assembly into each worker; use "
            "--mp-start-method fork on Linux/macOS to share it copy-on-write\n"
        )

    graph_entries = read_graph_entries(args.graph_list, args.graph_folder)
    pairs, gcigars = build_direct_inputs(args, graph_entries)
    _validate_direct_pairs(pairs, gcigars)
    eligible_core_coordinates = list({
        (
            getattr(pair, "output_query_core_coord", pair.query_coord).chrom,
            getattr(pair, "output_query_core_coord", pair.query_coord).start,
            getattr(pair, "output_query_core_coord", pair.query_coord).end,
            getattr(pair, "output_query_core_coord", pair.query_coord).strand,
        )
        for pair in pairs
        if getattr(pair, "output_query_core_coord", pair.query_coord) is not None
    })
    eligible_core_coordinates = [
        core.Coord(chrom, start, end, strand)
        for chrom, start, end, strand in eligible_core_coordinates
    ]
    query_genome = args.query_genome or infer_query_genome(args.align_query)
    lift_rows = read_genomelift(args.genomelift)
    synthetic_deletions = read_explicit_deletion_rows(
        args.genomelift, lift_rows, query_genome,
    )
    if synthetic_deletions:
        sys.stderr.write(
            "[graphcigartoref_persample] read "
            f"{len(synthetic_deletions)} explicit GenomeLift DEL row(s) for "
            "whole-deletion CIGAR output\n"
        )

    if args.unreplace:
        if not args.output:
            raise ValueError("--unreplace requires --output")
        candidate_groups = group_column10_candidates(pairs)
        before_pairs = len(pairs)
        before_comparisons = sum(len(group) for group in candidate_groups)
        before_candidate_groups = len(candidate_groups)
        wanted = {
            comparison_finished_key(comparison)
            for group in candidate_groups
            for comparison in group
        }
        wanted.update(row.finished_key() for row in synthetic_deletions)
        finished, old, kept, removed, malformed = (
            core.prune_existing_output_for_pairs(args.output, wanted)
        )
        pending_candidate_groups = [
            group for group in candidate_groups
            if not any(
                comparison_finished_key(comparison) in finished
                for comparison in group
            )
        ]
        pairs = [
            pair for group in pending_candidate_groups
            for _base_pair, output_pairs in group for pair in output_pairs
        ]
        synthetic_deletions = [
            row for row in synthetic_deletions
            if row.finished_key() not in finished
        ]
        sys.stderr.write(
            f"[graphcigartoref_persample] --unreplace: cleaned {args.output}; "
            f"kept {kept}/{old} existing rows; removed {removed}"
            f"{f' ({malformed} malformed/blank)' if malformed else ''}; "
            f"skipped {before_candidate_groups - len(pending_candidate_groups)}/"
            f"{before_candidate_groups} current query candidate group(s) "
            f"from {before_comparisons} reference-intersection comparison(s); "
            f"({before_pairs - len(pairs)}/{before_pairs} ownership part(s))\n"
        )

    try:
        query_index_reader = core.FastaRegionReader(args.fasta_query)
    except FileNotFoundError as error:
        raise FileNotFoundError(
            f"--fasta-query must be faidxed: {error}"
        ) from error
    query_index = dict(query_index_reader.index)
    query_index_reader.handle.close()
    assembly_sequences = core.read_fasta(args.fasta_query)
    missing_contigs = sorted(
        contig for contig in query_index if contig not in assembly_sequences
    )
    if missing_contigs:
        raise ValueError(
            "--fasta-query FAI contains contigs absent from the FASTA: "
            + ", ".join(repr(value) for value in missing_contigs[:10])
        )
    wrong_lengths = [
        (contig, len(assembly_sequences[contig]), values[0])
        for contig, values in query_index.items()
        if len(assembly_sequences[contig]) != values[0]
    ]
    if wrong_lengths:
        contig, observed, indexed = wrong_lengths[0]
        raise ValueError(
            f"--fasta-query and its FAI disagree for {contig!r}: FASTA has "
            f"{observed} bases, FAI reports {indexed}"
        )
    sys.stderr.write(
        f"[graphcigartoref_persample] loaded sample assembly into RAM: "
        f"{len(query_index)} contig(s), "
        f"{sum(values[0] for values in query_index.values())} bases\n"
    )
    try:
        # The reference lives in RAM once in this parent; fork workers share
        # the pages copy-on-write and skip per-fetch file seeks entirely.
        ref_reader = core.InMemoryFastaRegionReader(
            args.reference,
            template_fasta=(args.local_reference_templates or None),
        )
    except FileNotFoundError as error:
        raise FileNotFoundError(f"--reference must be faidxed: {error}") from error
    global _SHARED_REFERENCE_READER, _PARENT_REFERENCE_READER
    _PARENT_REFERENCE_READER = ref_reader
    if args.mp_start_method == "fork":
        _SHARED_REFERENCE_READER = ref_reader
    sys.stderr.write(
        "[graphcigartoref_persample] loaded reference into RAM: "
        f"{len(ref_reader._sequences)} contig(s), "
        f"{sum(map(len, ref_reader._sequences.values()))} bases\n"
    )

    output = None
    graph_sequences = None
    conversion_stats = _empty_conversion_stats()
    emitted_coordinates: List[core.Coord] = []
    try:
        graph_sequences, graph_path_coords, graph_mappings, files_read = (
            read_required_graph_data(
                graph_entries,
                pairs,
                gcigars,
                ref_reader,
                read_templates(args.local_reference_templates)
                if args.local_reference_templates else (),
                summary_dir=getattr(args, "graph_summary", ""),
                reference_haplotype=args.refhaplo,
            )
        )
        pairs = skip_unavailable_graph_candidates(
            pairs, gcigars, graph_sequences,
        )
        if alternative_policy is not None:
            restricted = alternative_policy.annotate_pairs(pairs, gcigars, graph_path_coords)
            sys.stderr.write(
                f'[graphcigartoref_persample] restricted {restricted}/{len(pairs)} '
                'comparisons to original alternative/novel intervals\n'
            )
        graph_source = (
            f"consolidated summary plus {files_read} fallback graph FASTA(s)"
            if getattr(args, "graph_summary", "")
            else f"{files_read} graph FASTA(s)"
        )
        virtual_slices = [
            value for value in dict.values(graph_sequences)
            if type(value) is tuple
        ]
        sys.stderr.write(
            "[graphcigartoref_persample] virtualized "
            f"{len(virtual_slices)} reference-slice graph record(s) "
            f"({sum(end - start for _c, start, end, _s, _u in virtual_slices)} "
            "bases served from the in-RAM reference instead of stored "
            "copies)\n"
        )
        graph_slices = [value for value in dict.values(graph_sequences)
                        if isinstance(value, _GraphReferenceSlice)]
        if graph_slices:
            sys.stderr.write(
                f"[graphcigartoref_persample] virtualized {len(graph_slices)} "
                f"graph-build-reference records ({sum(value.end - value.start for value in graph_slices)} "
                "bases served from indexed graph-build FASTA)\n"
            )
        sys.stderr.write(
            f"[graphcigartoref_persample] loaded "
            f"{len({core.path_key(key) for key in graph_sequences})} required "
            f"graph path(s) from {graph_source}\n"
        )
        chrom_lengths = {
            name: values[0] for name, values in ref_reader.index.items()
        }
        coord_by_name = _coord_by_record_name(pairs)
        allowchroms = set() if args.add_alt else core.default_allowed_reference_chroms(pairs)
        if args.local_reference_templates:
            allowchroms.update(
                template.output_contig
                for template in read_templates(args.local_reference_templates)
            )
        mode = "a" if args.unreplace else "w"
        output = open(args.output, mode) if args.output else sys.stdout
        if args.processes == 1:
            cache: Dict[Tuple[str, int, int, str], str] = {}
            writer = _CompletedRowBuffer(
                output,
                args.buffer_size,
                getattr(args, "buffer_bytes", 64 * 1024 * 1024),
            )
            try:
                for candidate_group in group_column10_candidates(pairs):
                    base_pair = candidate_group[0][0]
                    try:
                        row_texts, candidate_stats = (
                            _evaluate_column10_candidate_rows(
                            candidate_group,
                            gcigars,
                            graph_sequences,
                            graph_path_coords,
                            chrom_lengths,
                            graph_mappings,
                            ref_reader,
                            coord_by_name,
                            allowchroms,
                            assembly_sequences,
                            cache,
                            )
                        )
                    except Exception as error:
                        if not is_grouped_pair(base_pair):
                            raise
                        _add_conversion_stats(conversion_stats, {
                            "candidate_comparisons": len(candidate_group),
                            "candidate_errors": len(candidate_group),
                            "empty_groups": 1,
                        })
                        sys.stderr.write(
                            "[graphcigartoref_persample] warning: skipping "
                            f"unconvertible grouped candidate set "
                            f"{base_pair.label!r}: {error}\n"
                        )
                        continue
                    _add_conversion_stats(conversion_stats, candidate_stats)
                    for row_text in row_texts:
                        writer.write(row_text)
                        coord = _row_query_interval(row_text)
                        if coord is not None:
                            emitted_coordinates.append(coord)
            finally:
                writer.flush()
        elif pairs:
            shared_values = (
                gcigars,
                graph_sequences,
                graph_path_coords,
                chrom_lengths,
                graph_mappings,
                coord_by_name,
                assembly_sequences,
                allowchroms,
            )
            conversion_stats, emitted_coordinates = _write_parallel(
                pairs, output, args, shared_values,
            )
        for deletion in synthetic_deletions:
            chrom_length = chrom_lengths.get(deletion.ref_coord.chrom)
            if chrom_length is None:
                raise KeyError(
                    f"synthetic deletion reference chromosome "
                    f"{deletion.ref_coord.chrom!r} is absent from "
                    f"{args.reference}.fai"
                )
            rendered_deletion = render_synthetic_deletion_or_warn(
                deletion,
                chrom_length,
                "graphcigartoref_persample",
            )
            if rendered_deletion is not None:
                if alternative_policy is not None:
                    mask = alternative_policy.direct_tag(deletion.ref_coord.chrom)
                    if mask is None and any(
                        name.split('_', 1)[0] in alternative_policy.by_prefix
                        for name in deletion.ref_name.split(';')
                    ):
                        # No graph-CIGAR anchor exists for this synthetic row.
                        mask = alternative_policy.tag([])
                    rendered_deletion = core.append_tag(rendered_deletion, mask)
                output.write(rendered_deletion + "\n")
        output.flush()
        # A normal run already collected every emitted coordinate in the
        # parent and should not rescan a potentially multi-gigabyte CIGAR file.
        # Resume mode must include rows retained from the previous output.
        if args.output and args.unreplace:
            emitted_coordinates = _read_output_query_intervals(args.output)
        _log_conversion_coverage(
            eligible_core_coordinates,
            emitted_coordinates,
            conversion_stats,
        )
    finally:
        if output is not None and output is not sys.stdout:
            output.close()
        if isinstance(graph_sequences, ReferenceSliceSequences):
            graph_sequences.close()
        ref_reader.close()
    return 0


def main(argv: Optional[Sequence[str]] = None) -> int:
    try:
        return run(parse_args(argv))
    except Exception as error:
        print(f"graphcigartoref_persample.py: error: {error}", file=sys.stderr)
        return 1


# Worker-side stage profiling: pool workers exit through os._exit, so the
# core module's atexit dump never fires in them.  When profiling is enabled,
# wrap the per-group evaluators here and dump the accumulated core stage
# times periodically from inside the workers instead.
if getattr(core, "_PROFILE_ENABLED", False):
    def _persample_profile_wrap(fn, dump_every: int = 0):
        name = "persample." + fn.__name__

        def wrapped(*args, **kwargs):
            begin = core.time.perf_counter()
            try:
                return fn(*args, **kwargs)
            finally:
                core._STAGE_TIMES[name] = (
                    core._STAGE_TIMES.get(name, 0.0)
                    + core.time.perf_counter() - begin
                )
                count = core._STAGE_COUNTS.get(name, 0) + 1
                core._STAGE_COUNTS[name] = count
                if dump_every and count % dump_every == 0:
                    core._dump_stage_profile()
        wrapped.__name__ = fn.__name__
        return wrapped

    for _profiled_name, _dump_every in (
        ("_build_row", 512),
        ("_evaluate_column10_candidate_rows", 512),
        ("_query_record_sequences", 0),
        ("_worker_build_comparison", 0),
    ):
        _profiled_fn = globals().get(_profiled_name)
        if callable(_profiled_fn):
            globals()[_profiled_name] = _persample_profile_wrap(
                _profiled_fn, _dump_every,
            )


if __name__ == "__main__":
    raise SystemExit(main())
