#!/usr/bin/env python3
"""Determine conservative initial blocks for one partition graph.

This is the single-graph counterpart of
``summarize_partition_hotspot_segments.py``.  It intentionally reuses the same
projection, repeat-redundancy stabilization, validity-mask, and neighbor-join
functions so automatic ``minsetgraphmake --fixbreaks`` initialization has the
same semantics as the cohort-wide command.
"""
from __future__ import annotations

import argparse
import dataclasses
import logging
import os
from pathlib import Path
from typing import Iterable, List, Optional, Sequence

from align_partition_hotspots import OutputRow
from summarize_partition_hotspot_segments import (
    DEFAULT_MINIMUM_SEGMENT_SIZE,
    DEFAULT_SEGMENT_CUTOFF,
    DEFAULT_VALID_ALIGNMENT_SIZE,
    DEFAULT_VALID_REGION_MERGE_GAP,
    FinalSegment,
    OriginalBedInterval,
    build_graph_valid_regions,
    graph_name_from_fasta,
    intersect_segments_with_valid_regions,
    merge_projected_intervals,
    project_alignment_row,
    read_alignment_output,
    read_graph_path_records,
    select_alignment_rows,
    stabilize_segment_redundancy,
    validate_nonoverlapping_segments,
    write_final_segments,
)


LOG = logging.getLogger("partition_hotspot_segments")


def read_reference_bed(
    path: str, reference_haplotypes: Iterable[str],
) -> List[OriginalBedInterval]:
    accepted = set(reference_haplotypes)
    output: List[OriginalBedInterval] = []
    with open(path, "rt") as handle:
        for line_number, raw in enumerate(handle, 1):
            if not raw.strip() or raw.startswith("#"):
                continue
            fields = raw.rstrip("\r\n").split("\t")
            if len(fields) < 7:
                fields = raw.split()
            if len(fields) < 7:
                raise ValueError(f"{path}:{line_number}: expected BED7+")
            try:
                start, end = int(fields[1]), int(fields[2])
            except ValueError as error:
                raise ValueError(
                    f"{path}:{line_number}: invalid BED coordinates"
                ) from error
            strand, haplotype = fields[5], fields[6]
            if start < 0 or end <= start or strand not in {"+", "-"}:
                raise ValueError(f"{path}:{line_number}: invalid BED interval")
            if not accepted or haplotype in accepted:
                output.append(OriginalBedInterval(
                    haplotype, fields[0], start, end, strand,
                ))
    return output


def graph_path_reference_bed(
    graph_fasta: str, reference_haplotypes: Iterable[str],
) -> List[OriginalBedInterval]:
    """Use complete reference graph paths when no explicit BED is available."""
    paths = list(read_graph_path_records(graph_fasta).values())
    if not paths:
        raise ValueError(
            f"{graph_fasta}: graph paths lack HAPLOTYPE:CONTIG:START-END headers"
        )
    accepted = set(reference_haplotypes)
    output = [
        OriginalBedInterval(
            path.haplotype, path.contig, path.start, path.end, path.strand,
        )
        for path in paths if not accepted or path.haplotype in accepted
    ]
    if not output:
        raise ValueError(
            f"{graph_fasta}: no graph path matches reference haplotype(s) "
            + ",".join(sorted(accepted))
        )
    return output


def build_partition_segments(
    rows: Sequence[OutputRow],
    graph_fasta: str,
    reference_bed: Optional[str] = None,
    reference_haplotypes: Iterable[str] = (),
    scutoff: int = DEFAULT_SEGMENT_CUTOFF,
    dcutoff: int = DEFAULT_SEGMENT_CUTOFF,
    minimum_segment_size: int = DEFAULT_MINIMUM_SEGMENT_SIZE,
    valid_minimum_alignment: int = DEFAULT_VALID_ALIGNMENT_SIZE,
    valid_merge_gap: int = DEFAULT_VALID_REGION_MERGE_GAP,
    only_reference: bool = False,
) -> List[FinalSegment]:
    """Build conservative, non-overlapping initial blocks for one graph."""
    graph_fasta = os.path.abspath(graph_fasta)
    graph_name = graph_name_from_fasta(graph_fasta)
    accepted = tuple(reference_haplotypes)
    filtered = select_alignment_rows(
        [row for row in rows if row.hotspot_index == 1],
        accepted, only_reference,
    )
    if not filtered:
        raise ValueError("single-partition alignment contains no hotspot-index 1 rows")
    paths = read_graph_path_records(graph_fasta)
    bed_haplotypes = accepted if only_reference else ()
    if reference_bed:
        beds = read_reference_bed(reference_bed, bed_haplotypes)
    else:
        adjacent = os.path.join(os.path.dirname(graph_fasta), graph_name + ".bed")
        beds = (
            read_reference_bed(adjacent, bed_haplotypes)
            if os.path.isfile(adjacent)
            else graph_path_reference_bed(graph_fasta, bed_haplotypes)
        )
    if not beds:
        raise ValueError("reference BED/path selection produced no initial intervals")

    valid_regions = build_graph_valid_regions(
        filtered, paths, beds, valid_minimum_alignment, valid_merge_gap,
    )
    projected = [
        interval
        for row in filtered
        for interval in project_alignment_row(row, graph_name, paths, beds)
    ]
    merged = merge_projected_intervals(projected, scutoff, dcutoff)
    stabilized = stabilize_segment_redundancy(merged, minimum_segment_size)
    clipped = intersect_segments_with_valid_regions(
        stabilized, valid_regions, minimum_segment_size,
    )
    output = [dataclasses.replace(segment, strand="+") for segment in clipped]
    output.sort(key=lambda segment: (
        segment.query_contig, segment.query_start, segment.query_end,
        segment.graph_start, segment.graph_end,
    ))
    validate_nonoverlapping_segments(output)
    return output


def parse_args(argv: Optional[Sequence[str]] = None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("-i", "--input", required=True, help="ten-column alignment TSV")
    parser.add_argument("-g", "--graph", required=True, help="single graph .FA")
    parser.add_argument("-b", "--reference-bed")
    parser.add_argument("-o", "--output", required=True)
    parser.add_argument("--reference-haplotype", action="append", default=[])
    parser.add_argument(
        "--onlyref", "--only-ref", action="store_true",
        help="restrict rows and source intervals to reference haplotype(s)",
    )
    parser.add_argument("--scutoff", type=int, default=DEFAULT_SEGMENT_CUTOFF)
    parser.add_argument("--dcutoff", type=int, default=DEFAULT_SEGMENT_CUTOFF)
    parser.add_argument(
        "--min-segment-size", type=int, default=DEFAULT_MINIMUM_SEGMENT_SIZE,
    )
    parser.add_argument(
        "--valid-min-alignment", type=int, default=DEFAULT_VALID_ALIGNMENT_SIZE,
    )
    parser.add_argument(
        "--valid-merge-gap", type=int, default=DEFAULT_VALID_REGION_MERGE_GAP,
    )
    parser.add_argument("-v", "--verbose", action="store_true")
    args = parser.parse_args(argv)
    if args.scutoff < 0 or args.dcutoff < 1:
        parser.error("--scutoff must be nonnegative and --dcutoff positive")
    if args.min_segment_size < 0:
        parser.error("--min-segment-size must be nonnegative")
    return args


def main(argv: Optional[Sequence[str]] = None) -> int:
    args = parse_args(argv)
    logging.basicConfig(
        level=logging.DEBUG if args.verbose else logging.INFO,
        format="[%(asctime)s] %(levelname)s: %(message)s",
    )
    try:
        segments = build_partition_segments(
            read_alignment_output(args.input), args.graph,
            args.reference_bed,
            args.reference_haplotype or (["CHM13_h1"] if args.onlyref else []),
            args.scutoff, args.dcutoff, args.min_segment_size,
            args.valid_min_alignment, args.valid_merge_gap,
            args.onlyref,
        )
        write_final_segments(args.output, segments)
        LOG.info("Wrote %d initial partition blocks: %s", len(segments), args.output)
        return 0
    except Exception as error:
        LOG.error("%s", error)
        if args.verbose:
            LOG.exception("single-partition segment determination failed")
        return 1


if __name__ == "__main__":
    raise SystemExit(main())
