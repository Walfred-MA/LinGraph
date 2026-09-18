#!/usr/bin/env python3
"""Organize small-novel mappings and annotate their FASTA headers.

This utility is also the repair tool for a completed
``map_small_novel_minset.py`` run.  It does not run an aligner.  Instead, it
combines the accepted reference, large-novel, and final cleaned-minset
placements already written by that pipeline.
"""

from __future__ import annotations

import argparse
import bisect
import csv
import glob
import hashlib
import json
import logging
import os
import re
from collections import defaultdict
from itertools import chain
from pathlib import Path
from typing import Dict, Iterable, Iterator, List, Mapping, Optional, Sequence, Tuple

from minsetref_align import (
    PAF_INDEX_THRESHOLD_BYTES,
    filter_hits_by_unmasked_query,
    load_filtered_paf_hits_indexed,
    parse_blast17_sam,
    parse_winnowmap_paf,
)
from minsetref_core import fasta_lengths, merge_intervals, wrap_fasta
from minsetref_light import MASKED_BASE_WEIGHT
from summarize_partition_novelty import write_tsv


LOG = logging.getLogger("organize_small_novel_alignments")
ASSIGNMENT_ORDER = {
    "reference": 0,
    "large_novel": 1,
    "cleaned_small_novel": 2,
}
PIECE_SOURCE_RE = re.compile(r"^([^:]+):(.+):(\d+)-(\d+)([+-])$")
CIGAR_RE = re.compile(r"(\d+)([MIDNSHP=X])")


def read_dict_rows(path: str) -> Iterator[Dict[str, str]]:
    with open(path, "rt", newline="") as handle:
        first = handle.readline()
        if not first:
            return
        if first.startswith("#"):
            first = first[1:]
        fields = next(csv.reader([first], delimiter="\t"))
        yield from csv.DictReader(handle, delimiter="\t", fieldnames=fields)


def root_query_id(query_id: str) -> str:
    return query_id.split(".residual.", 1)[0]


def graphic_query_bases(path: str) -> int:
    return sum(
        int(length)
        for length, operation in CIGAR_RE.findall(path)
        if operation in {"M", "I", "=", "X"}
    )


def complete_graphic_path(path: str, query_length: int) -> str:
    """Make a saved graph path consume every emitted query base."""
    if not path or path == ".":
        return f"{query_length}I"
    if not (path.startswith((">", "<")) or re.fullmatch(r"\d+I", path)):
        return f"{query_length}I"
    missing = query_length - graphic_query_bases(path)
    if missing <= 0:
        return path
    match = re.search(r"(\d+H)$", path)
    if match:
        return path[:match.start()] + f"{missing}I" + match.group(1)
    return path + f"{missing}I"


def parse_json_intervals(value: str) -> List[Tuple[int, int]]:
    if not value or value == ".":
        return []
    decoded = json.loads(value)
    output: List[Tuple[int, int]] = []
    for item in decoded:
        if not isinstance(item, list) or len(item) != 2:
            raise ValueError(f"invalid interval JSON: {value!r}")
        start, end = int(item[0]), int(item[1])
        if end <= start:
            raise ValueError(f"invalid interval {start}-{end}")
        output.append((start, end))
    return output


def normalized_json(value: str, empty: str = "[]") -> str:
    if not value or value == ".":
        return empty
    return json.dumps(json.loads(value), separators=(",", ":"))


def target_coordinates(
    target_kind: str,
    target_haplotype: str,
    target_contig: str,
    target_record: str,
    start: str,
    end: str,
    strand: str,
) -> str:
    coordinate_name = target_contig if target_contig not in {"", "."} else target_record
    prefix = (
        f"{target_haplotype}:" if target_haplotype not in {"", "."}
        else ""
    )
    return f"{target_kind}:{prefix}{coordinate_name}:{start}-{end}{strand}"


def placement_rows(path: str, assignment: str) -> Iterator[Dict[str, str]]:
    if not os.path.isfile(path):
        raise FileNotFoundError(path)
    for row in read_dict_rows(path):
        query_id = row.get("query_id") or row.get("candidate_id")
        if not query_id:
            raise ValueError(f"{path}: placement row lacks query/candidate ID")
        target_record = row["target_record"]
        target_haplotype = row.get("ref_haplotype") or "."
        target_contig = row.get("ref_contig") or target_record
        target_start = row.get("ref_start") or "."
        target_end = row.get("ref_end") or "."
        strand = row.get("strand") or "."
        mapped = target_coordinates(
            assignment, target_haplotype, target_contig, target_record,
            target_start, target_end, strand,
        )
        liftover = row.get("liftover_path") or "."
        label = liftover if liftover != "." else mapped
        yield {
            "original_query_id": root_query_id(query_id),
            "query_id": query_id,
            "partition": (
                row.get("partition") or row.get("mapped_partition") or "."
            ),
            "mapped_type": assignment,
            "target_kind": assignment,
            "target_record": target_record,
            "target_haplotype": target_haplotype,
            "target_contig": target_contig,
            "target_start": target_start,
            "target_end": target_end,
            "strand": strand,
            "alignment_score": row.get("alignment_score") or "0",
            "aligned_score": (
                row.get("aligned_score") or row.get("alignment_score") or "0"
            ),
            "aligned_bases": row.get("aligned_bases") or "0",
            "identity": row.get("identity") or "0",
            "source": row.get("source") or ".",
            "query_blocks": normalized_json(row.get("query_blocks") or "[]"),
            "target_blocks": normalized_json(row.get("target_blocks") or "[]"),
            "paired_blocks": normalized_json(row.get("paired_blocks") or "[]"),
            "liftover_path": liftover,
            "mapped_coordinates": mapped,
            "label": label,
        }


def multilocus_reporting_rows(path: str) -> Iterator[Dict[str, str]]:
    """Read the corrected all-stage table, preserving every component."""
    for row in read_dict_rows(path):
        haplotype = row.get("ref_haplotype") or "."
        if haplotype == "large_novel":
            assignment = "large_novel"
        elif haplotype == "small_novel":
            assignment = "cleaned_small_novel"
        else:
            assignment = "reference"
        # Reuse the normal placement conversion without reopening the file.
        query_id = row.get("query_id") or row.get("candidate_id")
        if not query_id:
            raise ValueError(f"{path}: reporting row lacks candidate_id")
        target_record = row["target_record"]
        target_contig = row.get("ref_contig") or target_record
        start = row.get("ref_start") or "."
        end = row.get("ref_end") or "."
        strand = row.get("strand") or "."
        mapped = target_coordinates(
            assignment, haplotype, target_contig, target_record,
            start, end, strand,
        )
        liftover = row.get("liftover_path") or "."
        yield {
            "original_query_id": root_query_id(query_id),
            "query_id": query_id,
            "partition": (
                row.get("partition") or row.get("mapped_partition") or "."
            ),
            "mapped_type": assignment,
            "target_kind": assignment,
            "target_record": target_record,
            "target_haplotype": haplotype,
            "target_contig": target_contig,
            "target_start": start,
            "target_end": end,
            "strand": strand,
            "alignment_score": row.get("alignment_score") or "0",
            "aligned_score": row.get("alignment_score") or "0",
            "aligned_bases": row.get("aligned_bases") or "0",
            "identity": row.get("identity") or "0",
            "source": row.get("source") or ".",
            "query_blocks": normalized_json(row.get("query_blocks") or "[]"),
            "target_blocks": normalized_json(row.get("target_blocks") or "[]"),
            "paired_blocks": normalized_json(row.get("paired_blocks") or "[]"),
            "liftover_path": liftover,
            "mapped_coordinates": mapped,
            "label": liftover if liftover != "." else mapped,
        }


def final_annotation_rows(path: str) -> Iterator[Dict[str, str]]:
    if not os.path.isfile(path):
        raise FileNotFoundError(path)
    for row in read_dict_rows(path):
        if row.get("status") != "mapped":
            continue
        raw_kind = row.get("target_kind") or "."
        if raw_kind == "small_novel":
            assignment = "cleaned_small_novel"
        elif raw_kind == "large_novel":
            assignment = "large_novel"
        elif raw_kind in {"reference", "cleaned_reference"}:
            assignment = "reference"
        else:
            continue
        query_id = row["query_id"]
        target_record = row.get("target_record") or "."
        target_haplotype = row.get("target_kind") or "."
        target_contig = target_record
        target_start = row.get("target_start") or "."
        target_end = row.get("target_end") or "."
        strand = row.get("strand") or "."
        mapped = target_coordinates(
            assignment, target_haplotype, target_contig, target_record,
            target_start, target_end, strand,
        )
        liftover = row.get("liftover_path") or "."
        label = row.get("label") or "."
        if label == ".":
            label = liftover if liftover != "." else mapped
        yield {
            "original_query_id": root_query_id(query_id),
            "query_id": query_id,
            "partition": row.get("partition") or ".",
            "mapped_type": assignment,
            "target_kind": raw_kind,
            "target_record": target_record,
            "target_haplotype": target_haplotype,
            "target_contig": target_contig,
            "target_start": target_start,
            "target_end": target_end,
            "strand": strand,
            "alignment_score": row.get("alignment_score") or "0",
            "aligned_score": row.get("aligned_score") or "0",
            "aligned_bases": row.get("aligned_bases") or "0",
            "identity": row.get("identity") or "0",
            "source": row.get("source") or ".",
            "query_blocks": normalized_json(row.get("query_blocks") or "[]"),
            "target_blocks": normalized_json(row.get("target_blocks") or "[]"),
            "paired_blocks": normalized_json(row.get("paired_blocks") or "[]"),
            "liftover_path": liftover,
            "mapped_coordinates": mapped,
            "label": label,
        }


def direct_cleaned_rows(catalog_path: str) -> Iterator[Dict[str, str]]:
    if not os.path.isfile(catalog_path):
        raise FileNotFoundError(catalog_path)
    for row in read_dict_rows(catalog_path):
        query_id = row["origin_query_id"]
        record_id = row["smallnovel_id"]
        local_start = int(row["local_start"])
        local_end = int(row["local_end"])
        length = int(row["length"])
        mapped = f"cleaned_small_novel:{record_id}:0-{length}+"
        yield {
            "original_query_id": root_query_id(query_id),
            "query_id": query_id,
            "partition": row["partition"],
            "mapped_type": "cleaned_small_novel",
            "target_kind": "small_novel",
            "target_record": record_id,
            "target_haplotype": "small_novel",
            "target_contig": record_id,
            "target_start": "0",
            "target_end": str(length),
            "strand": "+",
            "alignment_score": str(row.get("novelty_score") or length),
            "aligned_score": str(row.get("novelty_score") or length),
            "aligned_bases": str(length),
            "identity": "100",
            "source": "retained_self_identity",
            "query_blocks": json.dumps([[local_start, local_end]], separators=(",", ":")),
            "target_blocks": json.dumps([[0, length]], separators=(",", ":")),
            "paired_blocks": json.dumps(
                [[local_start, local_end, 0, length]], separators=(",", ":"),
            ),
            "liftover_path": f">{record_id}:{length}=",
            "mapped_coordinates": mapped,
            "label": f"small_novel=>{record_id}:{length}=",
        }


def recovery_piece_sources(
    candidate_fasta: str,
) -> Dict[str, Tuple[str, str, int, int, str]]:
    output: Dict[str, Tuple[str, str, int, int, str]] = {}
    for record_id, description, _sequence in iter_fasta(candidate_fasta):
        tokens = description.split()
        if len(tokens) < 2:
            raise ValueError(
                f"{candidate_fasta}: {record_id!r} lacks a source coordinate"
            )
        match = PIECE_SOURCE_RE.fullmatch(tokens[1])
        if match is None:
            raise ValueError(
                f"{candidate_fasta}: invalid source coordinate {tokens[1]!r}"
            )
        sample, contig, start, end, strand = match.groups()
        output[record_id] = (sample, contig, int(start), int(end), strand)
    return output


def query_for_recovery_piece(
    sample: str,
    contig: str,
    start: int,
    end: int,
    query_index: Mapping[
        Tuple[str, str],
        Tuple[Sequence[int], Sequence[Tuple[str, Mapping[str, str]]]],
    ],
) -> Tuple[str, Mapping[str, str]]:
    starts, members = query_index.get((sample, contig), ((), ()))
    # Only intervals starting no later than this recovery piece can contain
    # it.  Grouping by assembly record avoids the former all-query scan for
    # every PAF/SAM hit (O(hits * all_queries)).
    stop = bisect.bisect_right(starts, start)
    matches = []
    for index in range(stop):
        query_id, query = members[index]
        if end <= int(query["local_end"]):
            matches.append((query_id, query))
    if not matches:
        raise ValueError(
            f"cannot assign recovery alignment piece "
            f"{sample}:{contig}:{start}-{end} to an original global query"
        )
    minimum = min(
        int(query["local_end"]) - int(query["local_start"])
        for _query_id, query in matches
    )
    best = [
        (query_id, query) for query_id, query in matches
        if int(query["local_end"]) - int(query["local_start"]) == minimum
    ]
    roots = {root_query_id(query_id) for query_id, _query in best}
    if len(roots) != 1:
        raise ValueError(
            f"ambiguous recovery provenance for "
            f"{sample}:{contig}:{start}-{end}: {','.join(sorted(roots))}"
        )
    return best[0]


def build_recovery_query_index(
    queries: Mapping[str, Mapping[str, str]],
) -> Dict[
    Tuple[str, str],
    Tuple[List[int], List[Tuple[str, Mapping[str, str]]]],
]:
    """Index global queries by assembly record and local interval start."""
    grouped: Dict[
        Tuple[str, str], List[Tuple[str, Mapping[str, str]]]
    ] = defaultdict(list)
    for query_id, query in queries.items():
        grouped[(query["haplotype"], query["record_id"])].append(
            (query_id, query)
        )
    output = {}
    for key, members in grouped.items():
        members.sort(key=lambda item: (
            int(item[1]["local_start"]),
            int(item[1]["local_end"]),
            item[0],
        ))
        output[key] = (
            [int(query["local_start"]) for _query_id, query in members],
            members,
        )
    return output


def find_blast_seqfile(cycle_dir: str) -> Optional[str]:
    candidates = (
        os.path.join(cycle_dir, "align.blastdb", "ref.seq"),
        os.path.join(os.path.dirname(cycle_dir), "reference_blastdb", "ref.seq"),
    )
    for path in candidates:
        if os.path.isfile(path):
            return path
    return None


def accepted_recovery_hits(
    cycle_dir: str,
    candidate_fasta: str,
    min_identity: float,
    min_score: int,
):
    hits = []
    for source in ("minimap2", "winnowmap"):
        paf = os.path.join(cycle_dir, f"align.{source}.paf")
        if not os.path.isfile(paf):
            continue
        if os.path.getsize(paf) > PAF_INDEX_THRESHOLD_BYTES:
            hits.extend(load_filtered_paf_hits_indexed(
                paf, candidate_fasta, min_identity, min_score, min_score,
                MASKED_BASE_WEIGHT, None, LOG, source=source,
                retain_aligned_pairs=True, aggregate_query_coverage=False,
            ))
        else:
            parsed = list(parse_winnowmap_paf(
                paf, min_identity, min_score, source=source,
                retain_aligned_pairs=True,
            ))
            hits.extend(filter_hits_by_unmasked_query(
                parsed, candidate_fasta, min_score, MASKED_BASE_WEIGHT,
            ))
    sam = os.path.join(cycle_dir, "align.blastn.sam")
    seqfile = find_blast_seqfile(cycle_dir)
    if os.path.isfile(sam) and seqfile is not None:
        parsed_blast = list(parse_blast17_sam(
            sam, seqfile, min_identity, min_score,
            query_lengths=fasta_lengths(candidate_fasta),
        ))
        hits.extend(filter_hits_by_unmasked_query(
            parsed_blast, candidate_fasta, min_score, MASKED_BASE_WEIGHT,
        ))
    return hits


def recovery_temp_alignment_rows(
    mapping_dir: str,
    queries: Mapping[str, Mapping[str, str]],
    min_identity: float,
    min_score: int,
) -> Iterator[Dict[str, str]]:
    """Recover exact cleaned-target mappings from preserved cycle outputs."""
    root = os.path.join(
        mapping_dir, "alignment_work", "self_clean", "coverage_recovery",
        "original_to_cleaned",
    )
    cycle_dirs = sorted(set(
        os.path.dirname(path) for path in glob.glob(
            os.path.join(root, "*", "cycle*", "candidates.fa"),
        )
    ))
    if not cycle_dirs:
        return
    LOG.info(
        "Recovering cleaned-small alignment coordinates from %d preserved "
        "alignment cycles under %s",
        len(cycle_dirs), root,
    )
    query_index = build_recovery_query_index(queries)
    provenance_cache: Dict[
        Tuple[str, str, int, int], Tuple[str, Mapping[str, str]]
    ] = {}
    for cycle_dir in cycle_dirs:
        candidate_fasta = os.path.join(cycle_dir, "candidates.fa")
        sources = recovery_piece_sources(candidate_fasta)
        for hit in accepted_recovery_hits(
            cycle_dir, candidate_fasta, min_identity, min_score,
        ):
            source = sources.get(hit.query_id)
            if source is None or hit.target_id in {"", "*"}:
                continue
            sample, contig, piece_start, piece_end, _piece_strand = source
            provenance_key = (sample, contig, piece_start, piece_end)
            provenance = provenance_cache.get(provenance_key)
            if provenance is None:
                provenance = query_for_recovery_piece(
                    sample, contig, piece_start, piece_end, query_index,
                )
                provenance_cache[provenance_key] = provenance
            original_query_id, query = provenance
            query_blocks = [
                (piece_start + start, piece_start + end)
                for start, end in hit.query_intervals if end > start
            ]
            target_blocks = [
                (start, end) for start, end in hit.target_intervals
                if end > start
            ]
            paired_blocks = [
                (piece_start + q0, piece_start + q1, t0, t1)
                for q0, q1, t0, t1 in hit.aligned_pairs
            ]
            if not query_blocks:
                continue
            target_start = min(
                (start for start, _end in target_blocks),
                default=hit.target_start,
            )
            target_end = max(
                (end for _start, end in target_blocks),
                default=hit.target_end,
            )
            mapped = (
                f"cleaned_small_novel:{hit.target_id}:"
                f"{target_start}-{target_end}{hit.strand}"
            )
            yield {
                "original_query_id": root_query_id(original_query_id),
                "query_id": (
                    f"{root_query_id(original_query_id)}.residual."
                    f"{piece_start}-{piece_end}"
                ),
                "partition": query["partition"],
                "mapped_type": "cleaned_small_novel",
                "target_kind": "small_novel",
                "target_record": hit.target_id,
                "target_haplotype": "small_novel",
                "target_contig": hit.target_id,
                "target_start": str(target_start),
                "target_end": str(target_end),
                "strand": hit.strand,
                "alignment_score": f"{hit.alignment_score:.6f}",
                "aligned_score": f"{hit.alignment_score:.6f}",
                "aligned_bases": str(hit.aligned_bases),
                "identity": f"{hit.identity:.6f}",
                "source": f"recovered_{hit.source}",
                "query_blocks": json.dumps(query_blocks, separators=(",", ":")),
                "target_blocks": json.dumps(target_blocks, separators=(",", ":")),
                "paired_blocks": json.dumps(paired_blocks, separators=(",", ":")),
                "liftover_path": ".",
                "mapped_coordinates": mapped,
                "label": mapped,
            }


ALIGNMENT_FIELDS = (
    "alignment_id", "original_query_id", "query_id", "partition",
    "mapped_type", "target_kind", "target_record", "target_haplotype",
    "target_contig", "target_start", "target_end", "strand",
    "alignment_score", "aligned_score", "aligned_bases", "identity",
    "source", "query_blocks", "target_blocks", "paired_blocks",
    "liftover_path", "mapped_coordinates", "label",
)


def deduplicate_alignments(rows: Iterable[Dict[str, str]]) -> List[Dict[str, str]]:
    selected: Dict[Tuple[str, ...], Dict[str, str]] = {}
    for row in rows:
        key = (
            row["original_query_id"], row["query_id"], row["mapped_type"],
            row["target_record"], row["strand"], row["query_blocks"],
            row["target_blocks"],
        )
        old = selected.get(key)
        if (
            old is None
            or float(row["alignment_score"]) > float(old["alignment_score"])
            or (
                float(row["alignment_score"]) == float(old["alignment_score"])
                and old["liftover_path"] == "."
                and row["liftover_path"] != "."
            )
        ):
            selected[key] = row
    output = sorted(selected.values(), key=lambda row: (
        row["original_query_id"],
        ASSIGNMENT_ORDER[row["mapped_type"]],
        row["query_id"], row["target_record"], row["target_start"],
    ))
    for index, row in enumerate(output, 1):
        row["alignment_id"] = f"smallnovel_alignment_{index:012d}"
    return output


def read_query_catalog(path: str) -> Dict[str, Dict[str, str]]:
    output: Dict[str, Dict[str, str]] = {}
    for row in read_dict_rows(path):
        query_id = row.get("original_query_id") or root_query_id(row["query_id"])
        if query_id in output:
            raise ValueError(f"{path}: duplicate original query {query_id!r}")
        output[query_id] = row
    return output


def build_summary(
    queries: Mapping[str, Mapping[str, str]],
    alignments: Sequence[Mapping[str, str]],
    terminal_statuses: Optional[Mapping[str, str]] = None,
) -> Tuple[List[Tuple[object, ...]], Dict[str, Dict[str, str]]]:
    terminal = terminal_statuses or {}
    by_query: Dict[str, List[Mapping[str, str]]] = defaultdict(list)
    for row in alignments:
        by_query[row["original_query_id"]].append(row)
    rows: List[Tuple[object, ...]] = []
    lookup: Dict[str, Dict[str, str]] = {}
    for query_id, query in queries.items():
        members = by_query.get(query_id, [])
        mapped_types = sorted(
            {row["mapped_type"] for row in members},
            key=ASSIGNMENT_ORDER.__getitem__,
        )
        terminal_reason = terminal.get(query_id)
        primary = (
            mapped_types[0] if mapped_types
            else terminal_reason or "unmapped"
        )
        coordinates = list(dict.fromkeys(
            row["mapped_coordinates"] for row in members
        ))
        combined_paths = list(dict.fromkeys(
            row["liftover_path"] for row in members
            if row["liftover_path"] != "."
        ))
        labels = combined_paths or list(dict.fromkeys(
            row["label"] for row in members if row["label"] != "."
        ))
        intervals: List[Tuple[int, int]] = []
        for member in members:
            intervals.extend(parse_json_intervals(member["query_blocks"]))
        covered = sum(end - start for start, end in merge_intervals(intervals))
        query_length = int(query["local_end"]) - int(query["local_start"])
        if covered >= query_length and query_length > 0:
            status = "mapped"
        elif covered > 0:
            status = "partial"
        elif terminal_reason is not None:
            status = "filtered"
        else:
            status = "unmapped"
        mapped_types_text = (
            ",".join(mapped_types) if mapped_types
            else terminal_reason or "."
        )
        coordinates_text = ";".join(coordinates) if coordinates else "."
        labels_text = ";".join(
            complete_graphic_path(label, query_length)
            for label in labels
        ) if labels else f"{query_length}I"
        assembly_coordinates = (
            f"{query['haplotype']}:{query['source_contig']}:"
            f"{query['source_start']}-{query['source_end']}+"
        )
        rows.append((
            query_id, query["partition"], query["haplotype"],
            query["source_contig"], query["source_start"],
            query["source_end"], assembly_coordinates, status,
            mapped_types_text, primary, coordinates_text, labels_text,
            len(members), covered, query_length,
            f"{covered / query_length:.6f}" if query_length else "0.000000",
        ))
        lookup[query_id] = {
            "status": status,
            "mapped_types": mapped_types_text,
            "primary_mapped_type": primary,
            "mapped_coordinates": coordinates_text,
            "labels": labels_text,
            "assembly_coordinates": assembly_coordinates,
        }
    return rows, lookup


def query_roots_in_tsv(path: str) -> set[str]:
    """Read provenance roots from either current or recovery audit TSVs."""
    if not os.path.isfile(path):
        return set()
    output: set[str] = set()
    for row in read_dict_rows(path):
        query_id = (
            row.get("original_query_id")
            or row.get("origin_query_id")
            or row.get("query_id")
        )
        if query_id:
            output.add(root_query_id(query_id))
    return output


def terminal_query_statuses(
    mapping_dir: str,
    layout_name: str,
    queries: Mapping[str, Mapping[str, str]],
) -> Dict[str, str]:
    """Identify query roots legitimately removed by the novelty threshold.

    Current runs record exact stage-04/stage-07 filtered intervals.  A
    completed legacy recovery predates those bundles, but its 08+09 union is
    the complete self-clean input universe.  Original roots outside that
    universe reached no score-passing self-clean segment and are therefore a
    terminal pre-self-clean filtered class, not a missing target assignment.
    """
    audit = os.path.join(mapping_dir, "intermediates")
    output: Dict[str, str] = {}
    for name in (
        "04_reference_filtered_small.tsv",
        "07_large_novel_filtered_small.tsv",
        "self_clean_recovery_filtered_below_score.tsv",
    ):
        for query_id in query_roots_in_tsv(os.path.join(audit, name)):
            output[query_id] = "filtered_below_score"

    self_clean_universe = set()
    for name in (
        "08_self_clean_retained.tsv",
        "09_self_clean_removed.tsv",
    ):
        self_clean_universe.update(
            query_roots_in_tsv(os.path.join(audit, name))
        )
    # One original query can yield both a below-threshold fragment and a
    # score-passing self-clean fragment.  The latter must still have an exact
    # target assignment; never let the small fragment hide missing evidence.
    for query_id in self_clean_universe:
        output.pop(query_id, None)

    if layout_name == "recovered_legacy":
        recovery_marker = os.path.join(
            audit, "_SELF_CLEAN_RECOVERY_SUCCESS.json",
        )
        if os.path.isfile(recovery_marker):
            for query_id in set(queries) - self_clean_universe:
                output.setdefault(query_id, "filtered_before_self_clean")
    return output


def write_unassigned_diagnostics(
    mapping_dir: str,
    unassigned: Sequence[str],
) -> str:
    """Write stage-membership evidence for genuinely unresolved queries."""
    audit = os.path.join(mapping_dir, "intermediates")
    stages = (
        "04_reference_filtered_small",
        "07_large_novel_filtered_small",
        "08_self_clean_retained",
        "09_self_clean_removed",
        "self_clean_recovery_input",
        "self_clean_recovery_filtered_below_score",
        "self_clean_recovery_represented",
        "self_clean_recovery_recovered",
    )
    memberships = {
        stage: query_roots_in_tsv(os.path.join(audit, stage + ".tsv"))
        for stage in stages
    }
    path = os.path.join(audit, "unassigned_small_novel_queries.tsv")
    write_tsv(
        path,
        ("original_query_id", *stages),
        (
            (
                query_id,
                *("1" if query_id in memberships[stage] else "0"
                  for stage in stages),
            )
            for query_id in unassigned
        ),
    )
    return path


SUMMARY_FIELDS = (
    "original_query_id", "partition", "haplotype", "source_contig",
    "source_start", "source_end", "assembly_coordinates", "status",
    "mapped_types", "primary_mapped_type", "mapped_coordinates", "labels",
    "alignment_count", "covered_bases", "query_length", "coverage_fraction",
)


def iter_fasta(path: str) -> Iterator[Tuple[str, str, str]]:
    name: Optional[str] = None
    description = ""
    chunks: List[str] = []
    with open(path, "rt") as handle:
        for raw in handle:
            if raw.startswith(">"):
                if name is not None:
                    yield name, description, "".join(chunks)
                description = raw[1:].strip()
                name = description.split()[0]
                chunks = []
            else:
                chunks.append(raw.strip())
    if name is not None:
        yield name, description, "".join(chunks)


def replace_header_fields(description: str, updates: Mapping[str, str]) -> str:
    parts = description.split()
    record_id = parts[0]
    retained = [
        token for token in parts[1:]
        if token.split("=", 1)[0] not in updates
    ]
    retained.extend(f"{key}={value}" for key, value in updates.items())
    return " ".join([record_id, *retained])


def smallnovel_header_updates(
    row: Mapping[str, str],
    summary: Mapping[str, Mapping[str, str]],
) -> Dict[str, str]:
    origin = root_query_id(row["origin_query_id"])
    origin_summary = summary.get(origin, {})
    length = int(row["length"])
    assembly = (
        f"{row['haplotype']}:{row['source_contig']}:"
        f"{row['source_start']}-{row['source_end']}+"
    )
    return {
        "assembly_coordinates": assembly,
        "mapped_type": "cleaned_small_novel",
        "mapped_coordinates": f"{row['smallnovel_id']}:0-{length}+",
        "graphic_cigar": origin_summary.get("labels", "."),
        "origin_mapped_types": origin_summary.get("mapped_types", "."),
    }


def annotate_smallnovel_fasta(
    fasta_path: str,
    catalog_path: str,
    summary: Mapping[str, Mapping[str, str]],
) -> None:
    catalog = {row["smallnovel_id"]: row for row in read_dict_rows(catalog_path)}
    temporary = fasta_path + f".tmp.{os.getpid()}"
    try:
        with open(temporary, "wt") as output:
            for record_id, description, sequence in iter_fasta(fasta_path):
                row = catalog.get(record_id)
                if row is None:
                    raise ValueError(f"{fasta_path}: {record_id!r} absent from {catalog_path}")
                header = replace_header_fields(
                    description, smallnovel_header_updates(row, summary),
                )
                output.write(f">{header}\n{wrap_fasta(sequence)}\n")
        os.replace(temporary, fasta_path)
        try:
            os.unlink(fasta_path + ".fai")
        except FileNotFoundError:
            pass
    finally:
        try:
            os.unlink(temporary)
        except FileNotFoundError:
            pass


def annotate_combined_minset_fasta(
    fasta_path: str,
    catalog_path: str,
    summary: Mapping[str, Mapping[str, str]],
) -> None:
    """Apply the same annotations to small records in novel_minset.fa."""
    catalog = {row["smallnovel_id"]: row for row in read_dict_rows(catalog_path)}
    temporary = fasta_path + f".tmp.{os.getpid()}"
    try:
        with open(temporary, "wt") as output:
            for record_id, description, sequence in iter_fasta(fasta_path):
                row = catalog.get(record_id)
                header = (
                    replace_header_fields(
                        description, smallnovel_header_updates(row, summary),
                    )
                    if row is not None else description
                )
                output.write(f">{header}\n{wrap_fasta(sequence)}\n")
        os.replace(temporary, fasta_path)
        try:
            os.unlink(fasta_path + ".fai")
        except FileNotFoundError:
            pass
    finally:
        try:
            os.unlink(temporary)
        except FileNotFoundError:
            pass


def annotate_all_query_fasta(
    mapping_dir: str,
    query_catalog_path: str,
    queries: Mapping[str, Mapping[str, str]],
    summary: Mapping[str, Mapping[str, str]],
) -> str:
    audit = os.path.join(mapping_dir, "intermediates")
    source_fasta = os.path.join(audit, "00_all_global_novel_queries.fa")
    audit_catalog = os.path.join(audit, "00_all_global_novel_queries.tsv")
    target = os.path.join(mapping_dir, "all_small_novel_queries.annotated.fa")
    temporary = target + f".tmp.{os.getpid()}"
    try:
        with open(temporary, "wt") as output:
            if (
                os.path.isfile(source_fasta)
                and os.path.isfile(audit_catalog)
                and os.path.samefile(query_catalog_path, audit_catalog)
            ):
                artifact_to_query = {
                    row["artifact_id"]: row["original_query_id"]
                    for row in read_dict_rows(query_catalog_path)
                }
                for artifact_id, description, sequence in iter_fasta(source_fasta):
                    query_id = artifact_to_query[artifact_id]
                    row = summary[query_id]
                    header = replace_header_fields(description, {
                        "assembly_coordinates": row["assembly_coordinates"],
                        "mapped_type": row["mapped_types"],
                        "mapped_coordinates": row["mapped_coordinates"],
                        "graphic_cigar": row["labels"],
                        "mapping_status": row["status"],
                    })
                    output.write(f">{header}\n{wrap_fasta(sequence)}\n")
            else:
                cache_fasta, cache_index = find_query_cache(mapping_dir)
                write_cached_query_annotations(
                    output, cache_fasta, cache_index, queries, summary,
                )
        os.replace(temporary, target)
    finally:
        try:
            os.unlink(temporary)
        except FileNotFoundError:
            pass
    return target


def find_query_cache(mapping_dir: str) -> Tuple[str, str]:
    cohort_dir = os.path.dirname(mapping_dir)
    candidates = (
        os.path.join(mapping_dir, "intermediates", "query_sequence_cache"),
        os.path.join(mapping_dir, "query_extraction.work"),
        os.path.join(
            cohort_dir, "global_novel_mapping", "global_query_extraction.work",
        ),
    )
    for directory in candidates:
        fasta = os.path.join(directory, "candidate_sequences.fa")
        index = os.path.join(directory, "candidate_sequences.index.tsv")
        if os.path.isfile(fasta) and os.path.isfile(index):
            return fasta, index
    raise FileNotFoundError(
        "no candidate_sequences.fa/index.tsv query cache found under "
        f"{mapping_dir} or its sibling global_novel_mapping"
    )


def write_cached_query_annotations(
    output,
    cache_fasta: str,
    cache_index: str,
    queries: Mapping[str, Mapping[str, str]],
    summary: Mapping[str, Mapping[str, str]],
) -> None:
    index: Dict[str, Tuple[int, int, int, int]] = {}
    for row in read_dict_rows(cache_index):
        index[row["candidate_id"]] = (
            int(row["local_start"]), int(row["local_end"]),
            int(row["sequence_offset"]), int(row["sequence_length"]),
        )
    descriptor = os.open(cache_fasta, os.O_RDONLY)
    try:
        for query_id, query in queries.items():
            cache_row = index.get(query_id)
            if cache_row is None:
                raise ValueError(
                    f"{cache_index}: missing cached sequence for {query_id!r}"
                )
            cache_start, cache_end, offset, sequence_length = cache_row
            encoded = os.pread(descriptor, sequence_length, offset)
            if len(encoded) != sequence_length:
                raise ValueError(
                    f"{cache_fasta}: incomplete cached sequence for {query_id!r}"
                )
            query_start, query_end = int(query["local_start"]), int(query["local_end"])
            if query_start < cache_start or query_end > cache_end:
                raise ValueError(
                    f"{query_id}: query {query_start}-{query_end} is outside "
                    f"cached interval {cache_start}-{cache_end}"
                )
            sequence = encoded[
                query_start - cache_start:query_end - cache_start
            ].decode("ascii")
            info = summary[query_id]
            output.write(
                f">{query_id} source={query['haplotype']}:"
                f"{query['source_contig']}:{query['source_start']}-"
                f"{query['source_end']}+ partition={query['partition']} "
                f"assembly_coordinates={info['assembly_coordinates']} "
                f"mapped_type={info['mapped_types']} "
                f"mapped_coordinates={info['mapped_coordinates']} "
                f"graphic_cigar={info['labels']} "
                f"mapping_status={info['status']}\n"
                f"{wrap_fasta(sequence)}\n"
            )
    finally:
        os.close(descriptor)


def fingerprint(path: str) -> Tuple[str, int, int]:
    stat = os.stat(path)
    return os.path.abspath(path), stat.st_size, stat.st_mtime_ns


def existing_paths(paths: Iterable[str]) -> List[str]:
    return [path for path in paths if os.path.isfile(path)]


def recursive_placement_tables(workdir: str) -> List[str]:
    """Find consolidated placement tables in old and new work layouts."""
    if not os.path.isdir(workdir):
        return []
    return sorted(set(glob.glob(
        os.path.join(workdir, "**", "placements.tsv"), recursive=True,
    )))


def tsv_has_fields(path: str, required: Sequence[str]) -> bool:
    if not os.path.isfile(path):
        return False
    with open(path, "rt", newline="") as handle:
        first = handle.readline()
    if first.startswith("#"):
        first = first[1:]
    fields = set(next(csv.reader([first], delimiter="\t"), ()))
    return set(required).issubset(fields)


def select_existing_layout(mapping_dir: str) -> Dict[str, object]:
    """Select either the current audit layout or the recovered legacy layout."""
    audit = os.path.join(mapping_dir, "intermediates")
    cohort_dir = os.path.dirname(mapping_dir)
    global_dir = os.path.join(cohort_dir, "global_novel_mapping")
    audit_queries = os.path.join(audit, "00_all_global_novel_queries.tsv")
    query_catalog = (
        audit_queries if os.path.isfile(audit_queries)
        else os.path.join(global_dir, "all_global_queries.tsv")
    )
    audit_reference = os.path.join(
        audit, "02_reference_aligned.placements.tsv",
    )
    reference_paths = (
        [audit_reference]
        if os.path.isfile(audit_reference)
        else [
            *existing_paths((
                os.path.join(
                    global_dir, "all_global_reference_placements.tsv",
                ),
            )),
            *recursive_placement_tables(os.path.join(
                mapping_dir, "alignment_work", "reference_completion",
            )),
        ]
    )
    audit_large = os.path.join(audit, "05_large_novel_aligned.placements.tsv")
    large_paths = (
        [audit_large]
        if os.path.isfile(audit_large)
        else recursive_placement_tables(os.path.join(
            mapping_dir, "alignment_work", "large_exclusion",
        ))
    )
    cleaned_paths = existing_paths((
        os.path.join(audit, "08_self_clean_mappings.placements.tsv"),
    ))
    recovered_catalog = os.path.join(audit, "08_self_clean_retained.tsv")
    recovered_fasta = os.path.join(audit, "08_self_clean_retained.fa")
    recovered_minset = os.path.join(mapping_dir, "novel_minset.recovered.fa")
    if (
        os.path.isfile(recovered_fasta)
        and tsv_has_fields(
            recovered_catalog, ("smallnovel_id", "origin_query_id"),
        )
    ):
        small_catalog = recovered_catalog
        small_fasta = recovered_fasta
        combined_fasta = (
            recovered_minset if os.path.isfile(recovered_minset)
            else os.path.join(mapping_dir, "novel_minset.fa")
        )
        layout = "recovered_legacy"
    else:
        small_catalog = os.path.join(mapping_dir, "smallnovel.tsv")
        small_fasta = os.path.join(mapping_dir, "smallnovel.fa")
        combined_fasta = os.path.join(mapping_dir, "novel_minset.fa")
        layout = "current"
    return {
        "layout": layout,
        "query_catalog": query_catalog,
        "reference_paths": reference_paths,
        "large_paths": large_paths,
        "cleaned_paths": cleaned_paths,
        "small_catalog": small_catalog,
        "small_fasta": small_fasta,
        "combined_fasta": combined_fasta,
        "final_annotations": os.path.join(
            mapping_dir, "all_small_novel_annotations.tsv",
        ),
        "cleaned_annotations": os.path.join(
            mapping_dir, "cleaned_small_novel_annotations.tsv",
        ),
        "multilocus_reporting": os.path.join(
            mapping_dir, "all_small_novel_reporting_alignments.tsv",
        ),
    }


def organize_mapping(
    mapping_dir: str,
    require_all: bool = True,
    min_identity: float = 95.0,
    min_score: int = 100,
) -> Dict[str, int]:
    mapping_dir = os.path.abspath(mapping_dir)
    audit = os.path.join(mapping_dir, "intermediates")
    layout = select_existing_layout(mapping_dir)
    query_catalog_path = str(layout["query_catalog"])
    reference_paths = list(layout["reference_paths"])
    large_paths = list(layout["large_paths"])
    cleaned_paths = list(layout["cleaned_paths"])
    final_path = str(layout["final_annotations"])
    cleaned_annotations_path = str(layout["cleaned_annotations"])
    multilocus_reporting_path = str(layout["multilocus_reporting"])
    small_catalog = str(layout["small_catalog"])
    small_fasta = str(layout["small_fasta"])
    combined_fasta = str(layout["combined_fasta"])
    required = (
        query_catalog_path, small_catalog, small_fasta,
        combined_fasta, *reference_paths, *large_paths,
        *cleaned_paths,
        *((final_path,) if os.path.isfile(final_path) else ()),
        *((cleaned_annotations_path,)
          if os.path.isfile(cleaned_annotations_path) else ()),
        *((multilocus_reporting_path,)
          if os.path.isfile(multilocus_reporting_path) else ()),
    )
    for path in required:
        if not os.path.isfile(path):
            raise FileNotFoundError(path)
    if not reference_paths and not os.path.isfile(multilocus_reporting_path):
        raise FileNotFoundError(
            "no reference placement table found in intermediates, "
            "global_novel_mapping, or alignment_work/reference_completion"
        )
    LOG.info(
        "Reference placement tables: %s",
        ",".join(reference_paths),
    )
    if not large_paths:
        LOG.warning(
            "No accepted large-novel placement table was found; continuing "
            "because this is valid when that stage accepted zero mappings"
        )
    else:
        LOG.info(
            "Large-novel placement tables: %s",
            ",".join(large_paths),
        )

    queries = read_query_catalog(query_catalog_path)
    has_multilocus = os.path.isfile(multilocus_reporting_path)
    recovery_rows = (
        recovery_temp_alignment_rows(
            mapping_dir, queries, min_identity, min_score,
        )
        if layout["layout"] == "recovered_legacy" and not has_multilocus
        else ()
    )
    final_rows = (
        final_annotation_rows(final_path)
        if os.path.isfile(final_path) and not has_multilocus else ()
    )
    cleaned_annotation_rows = (
        final_annotation_rows(cleaned_annotations_path)
        if os.path.isfile(cleaned_annotations_path) and not has_multilocus
        else ()
    )
    if not os.path.isfile(final_path) and not has_multilocus:
        LOG.info(
            "Final all-query annotation table is absent; reconstructing "
            "cleaned mappings from preserved recovery alignment cycles"
        )
    if has_multilocus:
        LOG.info(
            "Using corrected multi-locus all-stage report: %s",
            multilocus_reporting_path,
        )
        alignment_sources = (
            multilocus_reporting_rows(multilocus_reporting_path),
        )
    else:
        alignment_sources = (
            *(placement_rows(path, "reference") for path in reference_paths),
            *(placement_rows(path, "large_novel") for path in large_paths),
            *(placement_rows(path, "cleaned_small_novel") for path in cleaned_paths),
            final_rows,
            cleaned_annotation_rows,
            recovery_rows,
            direct_cleaned_rows(small_catalog),
        )
    alignments = deduplicate_alignments(chain(*alignment_sources))
    unknown = {row["original_query_id"] for row in alignments} - set(queries)
    if unknown:
        raise ValueError(
            "alignment summary contains unknown original queries: "
            + ",".join(sorted(unknown)[:20])
        )
    detail_path = os.path.join(mapping_dir, "all_small_novel_alignments.tsv")
    write_tsv(
        detail_path, ALIGNMENT_FIELDS,
        (tuple(row[field] for field in ALIGNMENT_FIELDS) for row in alignments),
    )
    terminal_status = terminal_query_statuses(
        mapping_dir, str(layout["layout"]), queries,
    )
    summary_rows, summary_lookup = build_summary(
        queries, alignments, terminal_status,
    )
    summary_path = os.path.join(mapping_dir, "small_novel_alignment_summary.tsv")
    write_tsv(summary_path, SUMMARY_FIELDS, summary_rows)
    unassigned = [
        query_id for query_id, row in summary_lookup.items()
        if row["status"] == "unmapped"
    ]
    if require_all and unassigned:
        diagnostic_path = write_unassigned_diagnostics(
            mapping_dir, unassigned,
        )
        raise ValueError(
            f"{len(unassigned)} small-novel queries have no accepted "
            "reference, large-novel, or cleaned-small alignment; first: "
            + ",".join(unassigned[:20])
            + f"; stage diagnostics: {diagnostic_path}"
        )
    annotate_smallnovel_fasta(small_fasta, small_catalog, summary_lookup)
    annotate_combined_minset_fasta(
        combined_fasta, small_catalog, summary_lookup,
    )
    annotate_all_query_fasta(
        mapping_dir, query_catalog_path, queries, summary_lookup,
    )

    marker = os.path.join(audit, "_ALIGNMENT_SUMMARY_SUCCESS.json")
    payload = {
        "version": "small-novel-alignment-summary-v1",
        "layout": layout["layout"],
        "min_identity": min_identity,
        "min_score": min_score,
        "inputs": [fingerprint(path) for path in required],
        "queries": len(queries),
        "alignments": len(alignments),
        "unassigned_queries": len(unassigned),
        "annotated_small_fasta": small_fasta,
        "annotated_combined_fasta": combined_fasta,
        "detail": detail_path,
        "summary": summary_path,
    }
    temporary = marker + f".tmp.{os.getpid()}"
    with open(temporary, "wt") as output:
        json.dump(payload, output, sort_keys=True)
        output.write("\n")
    os.replace(temporary, marker)
    LOG.info(
        "Organized %d alignments for %d original small-novel queries; "
        "%d remain unassigned",
        len(alignments), len(queries), len(unassigned),
    )
    return {
        "queries": len(queries),
        "alignments": len(alignments),
        "unassigned": len(unassigned),
    }


def parse_args(argv: Optional[Sequence[str]] = None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "-m", "--mapping-dir", required=True,
        help="completed output directory of map_small_novel_minset.py",
    )
    parser.add_argument(
        "--allow-unmapped", action="store_true",
        help="write unmapped summary rows instead of failing",
    )
    parser.add_argument("--min-identity", type=float, default=95.0)
    parser.add_argument("--min-score", type=int, default=100)
    parser.add_argument(
        "--log-level", choices=("DEBUG", "INFO", "WARNING", "ERROR"),
        default="INFO",
    )
    args = parser.parse_args(argv)
    if not 0 < args.min_identity <= 100:
        parser.error("--min-identity must be in (0,100]")
    if args.min_score < 1:
        parser.error("--min-score must be positive")
    return args


def main(argv: Optional[Sequence[str]] = None) -> int:
    args = parse_args(argv)
    logging.basicConfig(
        level=getattr(logging, args.log_level),
        format="[%(asctime)s] %(levelname)s: %(message)s",
        datefmt="%Y-%m-%d %H:%M:%S",
    )
    try:
        organize_mapping(
            args.mapping_dir,
            require_all=not args.allow_unmapped,
            min_identity=args.min_identity,
            min_score=args.min_score,
        )
    except Exception as error:
        LOG.error("%s", error)
        if args.log_level == "DEBUG":
            raise
        return 1
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
