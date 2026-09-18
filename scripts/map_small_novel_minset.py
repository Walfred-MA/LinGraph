#!/usr/bin/env python3
"""Build and annotate a nonredundant small-novel minset.

This stage runs after ``map_global_novel.py``. Reference-uncovered global-query
pieces are first cleaned against the established large ``novel_loci.fa``.
Surviving pieces are self-cleaned with Minimap2 followed by
Winnowmap+BLASTN. The fixed large loci plus the selected small representatives
form one novel minset, and every original global query is finally mapped to
that combined target. All mappings are annotations: original query sequence is
never removed or replaced.
"""

from __future__ import annotations

import argparse
import csv
import dataclasses
import hashlib
import json
import logging
import os
import re
import shutil
import sys
from collections import defaultdict
from collections.abc import Mapping as MappingABC
from pathlib import Path
from types import SimpleNamespace
from typing import Dict, Iterator, List, Mapping, Optional, Sequence, Tuple

from map_global_novel import (
    annotate_global_reporting_components,
    build_disk_query_cache,
    build_graphic_liftover_path,
)
from map_partition_local_novel import (
    REPORTING_RESOLUTION_VERSION,
    CachedCandidateSequence,
    load_candidate_header_records,
    read_assembly_sources,
    read_partition_inputs,
    read_reporting_components,
    reporting_components_to_evidence,
    write_reporting_components,
)
from minsetref_core import (
    IndexedFasta,
    Piece,
    count_masked,
    count_unmasked,
    merge_intervals,
    subtract_intervals,
    wrap_fasta,
    write_fasta_pieces,
)
from minsetref_light import (
    InputRecord,
    MASKED_BASE_WEIGHT,
    build_haplotype_priorities,
    self_clean_two_stage,
)
from minsetref_segments import novelty_score
from minsetref_graphic import build_multilocus_graphic_path
from minsetref_multilocus import (
    PlacementEvidence,
    ResolvedPlacementComponent,
    resolve_all_evidence,
)
from organize_small_novel_alignments import organize_mapping
from refine_partition_paths import (
    GlobalPlacement,
    NovelQuery,
    ReferenceRecord,
    align_global_queries,
    global_mapping_candidates,
    global_residual_queries,
    read_global_placements,
    read_novel_queries,
)
from summarize_partition_novelty import write_tsv
from fixed_alternatives import (
    ANNOTATION_POLICY, annotation_reference_records, annotation_reference_sources,
)


LOG = logging.getLogger("map_small_novel_minset")
VERSION = "small-novel-minset-v3-audit-intermediates"


@dataclasses.dataclass(frozen=True)
class SmallNovelRecord:
    smallnovel_id: str
    haplotype: str
    source_contig: str
    source_start: int
    source_end: int
    partition: str
    record_id: str
    local_start: int
    local_end: int
    origin_query_id: str
    novelty_score: float
    length: int
    sequence: str = dataclasses.field(repr=False, compare=False)


class ResidualQueryCache(MappingABC[str, CachedCandidateSequence]):
    """Resolve residual IDs through the disk cache of their original query."""

    retain_cycle_sequences = False

    def __init__(
        self,
        original: Mapping[str, CachedCandidateSequence],
        query_ids: Sequence[str],
    ) -> None:
        self.original = original
        self.query_ids = tuple(query_ids)

    @staticmethod
    def original_id(query_id: str) -> str:
        return query_id.split(".residual.", 1)[0]

    def __getitem__(self, query_id: str) -> CachedCandidateSequence:
        return self.original[self.original_id(query_id)]

    def __iter__(self) -> Iterator[str]:
        return iter(self.query_ids)

    def __len__(self) -> int:
        return len(self.query_ids)


def fingerprint(path: str) -> Tuple[str, int, int]:
    stat = os.stat(path)
    return os.path.abspath(path), stat.st_size, stat.st_mtime_ns


def fasta_target_records(
    fasta_path: str,
    small_ids: Sequence[str] = (),
) -> Dict[str, ReferenceRecord]:
    small = set(small_ids)
    output: Dict[str, ReferenceRecord] = {}
    with IndexedFasta(fasta_path) as fasta:
        for record_id in fasta.names():
            if record_id in output:
                raise ValueError(f"{fasta_path}: duplicate FASTA ID {record_id!r}")
            length = fasta.length(record_id)
            output[record_id] = ReferenceRecord(
                record_id,
                "small_novel" if record_id in small else "large_novel",
                record_id,
                0,
                length,
                "+",
                length,
            )
    return output


def concatenate_fastas(paths: Sequence[str], output_path: str) -> None:
    temporary = output_path + f".tmp.{os.getpid()}"
    try:
        with open(temporary, "wb") as output:
            for path in paths:
                if not os.path.isfile(path) or os.path.getsize(path) == 0:
                    continue
                with open(path, "rb") as source:
                    shutil.copyfileobj(source, output, length=16 * 1024 * 1024)
                    source.seek(-1, os.SEEK_END)
                    if source.read(1) not in {b"\n", b"\r"}:
                        output.write(b"\n")
        os.replace(temporary, output_path)
        try:
            os.unlink(output_path + ".fai")
        except FileNotFoundError:
            pass
    finally:
        try:
            os.unlink(temporary)
        except FileNotFoundError:
            pass


def query_sequence(
    query: NovelQuery,
    cache: Mapping[str, CachedCandidateSequence],
) -> str:
    cached = cache[ResidualQueryCache.original_id(query.query_id)]
    start = query.local_start - cached.local_start
    end = query.local_end - cached.local_start
    sequence = cached.sequence[start:end]
    if len(sequence) != query.local_end - query.local_start:
        raise ValueError(
            f"{query.query_id}: incomplete sequence cache for "
            f"{query.local_start}-{query.local_end}"
        )
    return sequence


def annotate_fixed_target_reporting_components(
    queries: Sequence[NovelQuery],
    components: Sequence[ResolvedPlacementComponent],
    cached: Mapping[str, CachedCandidateSequence],
    target_fasta: str,
) -> List[ResolvedPlacementComponent]:
    """Attach ordered paths for large-novel or cleaned-small targets."""
    query_by_id = {query.query_id: query for query in queries}
    grouped: Dict[str, List[ResolvedPlacementComponent]] = defaultdict(list)
    for component in components:
        grouped[component.query_id].append(component)
    output: List[ResolvedPlacementComponent] = []
    with IndexedFasta(target_fasta) as target:
        for query_id in sorted(grouped):
            query = query_by_id.get(query_id)
            if query is None:
                raise ValueError(
                    f"reporting component references unknown query {query_id!r}"
                )
            members = sorted(grouped[query_id], key=lambda row: (
                row.query_start, row.query_end, row.component_index,
            ))
            references = {}
            for member in members:
                if member.ref_contig not in target.index:
                    raise ValueError(
                        f"{query_id}: target record {member.ref_contig!r} "
                        f"is absent from {target_fasta}"
                    )
                references[(
                    member.ref_haplotype, member.ref_contig,
                )] = target
            sequence = query_sequence(query, cached)
            path = build_multilocus_graphic_path(
                query_id=query_id,
                query_start=query.local_start,
                query_end=query.local_end,
                query_sequence=sequence,
                components=members,
                references=references,
            )
            output.extend(
                dataclasses.replace(member, liftover_path=path)
                for member in members
            )
    return output


def reporting_components_from_placements(
    placements: Mapping[str, GlobalPlacement],
) -> List[ResolvedPlacementComponent]:
    """Compatibility conversion for synthetic/legacy placement tables."""
    output: List[ResolvedPlacementComponent] = []
    for placement in placements.values():
        if not placement.query_blocks or not placement.paired_blocks:
            continue
        output.append(ResolvedPlacementComponent(
            placement.query_id, 1, placement.partition,
            placement.ref_haplotype, placement.ref_contig,
            placement.ref_start, placement.ref_end, placement.strand,
            placement.target_record,
            min(start for start, _end in placement.query_blocks),
            max(end for _start, end in placement.query_blocks),
            tuple(placement.query_blocks), tuple(placement.target_blocks),
            tuple(placement.paired_blocks), placement.alignment_score,
            placement.aligned_bases, placement.identity, placement.source,
            0, placement.liftover_path,
        ))
    return output


def read_stage_reporting_or_placements(
    path: str,
    placements: Mapping[str, GlobalPlacement],
) -> List[ResolvedPlacementComponent]:
    if os.path.isfile(path):
        return read_reporting_components((path,))
    return reporting_components_from_placements(placements)


def combine_small_novel_reporting_components(
    queries: Sequence[NovelQuery],
    stage_components: Sequence[
        Tuple[str, Sequence[ResolvedPlacementComponent]]
    ],
    cached: Mapping[str, CachedCandidateSequence],
    assembly_sources: Mapping[str, object],
    large_fasta: str,
    small_fasta: str,
) -> List[ResolvedPlacementComponent]:
    """Build one full reference/large/small path per original query.

    Stage order mirrors cleaning and supplies deterministic tie order.
    Accepted cores from successive stages are normally disjoint because each
    later stage receives only the previous stage's residual query intervals.
    """
    query_by_id = {query.query_id: query for query in queries}
    evidence: List[PlacementEvidence] = []
    valid_stages = {
        "reference", "large_novel", "cleaned_small_novel",
    }
    ordinal = 0
    for stage, components in stage_components:
        if stage not in valid_stages:
            raise ValueError(f"unknown reporting stage {stage!r}")
        for component in components:
            root = query_root_id(component.query_id)
            if root not in query_by_id:
                raise ValueError(
                    f"{stage}: unknown original global query {root!r}"
                )
            for row in reporting_components_to_evidence((component,)):
                evidence.append(dataclasses.replace(
                    row,
                    query_id=root,
                    alignment_score=row.alignment_score,
                    source=f"{stage}:{row.source}",
                    ordinal=ordinal,
                ))
                ordinal += 1
    resolved = resolve_all_evidence(evidence)
    readers: Dict[str, IndexedFasta] = {}
    output: List[ResolvedPlacementComponent] = []
    try:
        large_reader = IndexedFasta(large_fasta)
        small_reader = IndexedFasta(small_fasta)
        readers["large_novel"] = large_reader
        readers["small_novel"] = small_reader
        for query_id in sorted(resolved):
            query = query_by_id[query_id]
            members = resolved[query_id]
            references = {}
            for member in members:
                if member.ref_haplotype == "large_novel":
                    reader = large_reader
                elif member.ref_haplotype == "small_novel":
                    reader = small_reader
                else:
                    source = assembly_sources.get(member.ref_haplotype)
                    if source is None:
                        raise ValueError(
                            f"{query_id}: no reference sequence source for "
                            f"{member.ref_haplotype!r}"
                        )
                    reader = readers.get(member.ref_haplotype)
                    if reader is None:
                        reader = IndexedFasta(
                            source.fasta_path, source.fai_path,
                        )
                        readers[member.ref_haplotype] = reader
                references[(
                    member.ref_haplotype, member.ref_contig,
                )] = reader
            path = build_multilocus_graphic_path(
                query_id=query_id,
                query_start=query.local_start,
                query_end=query.local_end,
                query_sequence=query_sequence(query, cached),
                components=members,
                references=references,
            )
            output.extend(
                dataclasses.replace(member, liftover_path=path)
                for member in members
            )
    finally:
        for reader in readers.values():
            reader.close()
    return output


def query_interval(
    query: NovelQuery,
    local_start: int,
    local_end: int,
) -> NovelQuery:
    """Return one source-coordinate-preserving subinterval of a query."""
    if (
        local_start < query.local_start
        or local_end <= local_start
        or local_end > query.local_end
    ):
        raise ValueError(
            f"{query.query_id}: invalid audit interval "
            f"{local_start}-{local_end} outside "
            f"{query.local_start}-{query.local_end}"
        )
    source_start = (
        query.source_start + local_start - query.local_start
    )
    return dataclasses.replace(
        query,
        source_start=source_start,
        source_end=source_start + local_end - local_start,
        local_start=local_start,
        local_end=local_end,
    )


def query_root_id(query_id: str) -> str:
    return query_id.split(".residual.", 1)[0]


def normalized_query_segments(
    base_queries: Sequence[NovelQuery],
    segments: Sequence[NovelQuery],
) -> List[NovelQuery]:
    """Union segment coordinates and express them on the base query rows."""
    intervals: Dict[str, List[Tuple[int, int]]] = defaultdict(list)
    for segment in segments:
        intervals[query_root_id(segment.query_id)].append((
            segment.local_start, segment.local_end,
        ))
    output: List[NovelQuery] = []
    for query in base_queries:
        for start, end in merge_intervals(
            intervals.get(query_root_id(query.query_id), ())
        ):
            start = max(start, query.local_start)
            end = min(end, query.local_end)
            if end > start:
                output.append(query_interval(query, start, end))
    return output


def placement_query_segments(
    queries: Sequence[NovelQuery],
    placements: Mapping[str, GlobalPlacement],
) -> List[NovelQuery]:
    """Materialize the accepted query blocks represented by placements."""
    output: List[NovelQuery] = []
    for query in queries:
        placement = placements.get(query.query_id)
        if placement is None:
            continue
        for start, end in merge_intervals(placement.query_blocks):
            start = max(start, query.local_start)
            end = min(end, query.local_end)
            if end > start:
                output.append(query_interval(query, start, end))
    return output


def subtract_query_segments(
    queries: Sequence[NovelQuery],
    masks: Sequence[NovelQuery],
) -> List[NovelQuery]:
    """Subtract query-local mask intervals without an all-pairs scan."""
    mask_by_query: Dict[str, List[Tuple[int, int]]] = defaultdict(list)
    for mask in masks:
        mask_by_query[query_root_id(mask.query_id)].append((
            mask.local_start, mask.local_end,
        ))
    for query_id in tuple(mask_by_query):
        mask_by_query[query_id] = merge_intervals(mask_by_query[query_id])
    output: List[NovelQuery] = []
    for query in queries:
        intervals = subtract_intervals(
            [(query.local_start, query.local_end)],
            mask_by_query.get(query_root_id(query.query_id), ()),
        )
        output.extend(
            query_interval(query, start, end)
            for start, end in intervals
        )
    return output


def retained_piece_queries(
    pieces: Sequence[Piece],
    queries: Sequence[NovelQuery],
) -> List[NovelQuery]:
    query_by_id = {query.query_id: query for query in queries}
    output: List[NovelQuery] = []
    for piece in pieces:
        query = query_by_id.get(piece.origin_id)
        if query is None:
            raise KeyError(
                f"self-clean output has unknown origin {piece.origin_id!r}"
            )
        output.append(query_interval(query, piece.start, piece.end))
    return normalized_query_segments(queries, output)


def write_query_audit_bundle(
    prefix: str,
    label: str,
    queries: Sequence[NovelQuery],
    cached: Mapping[str, CachedCandidateSequence],
) -> Tuple[str, str, str]:
    """Write one exact query interval set as FASTA, BED, and TSV."""
    fasta_path = prefix + ".fa"
    bed_path = prefix + ".bed"
    tsv_path = prefix + ".tsv"
    Path(os.path.dirname(prefix)).mkdir(parents=True, exist_ok=True)
    fasta_tmp = fasta_path + f".tmp.{os.getpid()}"
    bed_tmp = bed_path + f".tmp.{os.getpid()}"
    tsv_tmp = tsv_path + f".tmp.{os.getpid()}"
    fields = (
        "artifact_id", "stage", "query_id", "original_query_id",
        "origin", "partition", "interval_id", "haplotype",
        "source_contig", "source_start", "source_end", "record_id",
        "local_start", "local_end", "length", "unmasked_bases",
        "masked_bases", "novelty_score",
    )
    try:
        with open(fasta_tmp, "wt") as fasta_out, open(
            bed_tmp, "wt"
        ) as bed_out, open(tsv_tmp, "wt", newline="") as tsv_out:
            writer = csv.writer(
                tsv_out, delimiter="\t", lineterminator="\n",
            )
            writer.writerow(fields)
            for index, query in enumerate(queries, 1):
                artifact_id = f"{label}_{index:012d}"
                sequence = query_sequence(query, cached)
                unmasked = count_unmasked(sequence)
                masked = count_masked(sequence)
                score = novelty_score(sequence, MASKED_BASE_WEIGHT)
                fasta_out.write(
                    f">{artifact_id} source={query.haplotype}:"
                    f"{query.source_contig}:{query.source_start}-"
                    f"{query.source_end} query={query.query_id} "
                    f"partition={query.partition} record={query.record_id} "
                    f"local={query.local_start}-{query.local_end} "
                    f"score={score:.1f}\n{wrap_fasta(sequence)}\n"
                )
                bed_out.write(
                    f"{query.source_contig}\t{query.source_start}\t"
                    f"{query.source_end}\t{artifact_id}\t0\t+\t"
                    f"{query.haplotype}\t{query.partition}\t"
                    f"{query.record_id}\t{query_root_id(query.query_id)}\t"
                    f"{query.local_start}\t{query.local_end}\t"
                    f"{score:.1f}\n"
                )
                writer.writerow((
                    artifact_id, label, query.query_id,
                    query_root_id(query.query_id), query.origin,
                    query.partition, query.interval_id, query.haplotype,
                    query.source_contig, query.source_start,
                    query.source_end, query.record_id, query.local_start,
                    query.local_end, len(sequence), unmasked, masked,
                    f"{score:.1f}",
                ))
        os.replace(fasta_tmp, fasta_path)
        os.replace(bed_tmp, bed_path)
        os.replace(tsv_tmp, tsv_path)
    finally:
        for path in (fasta_tmp, bed_tmp, tsv_tmp):
            try:
                os.unlink(path)
            except FileNotFoundError:
                pass
    return fasta_path, bed_path, tsv_path


def write_placement_audit(
    path: str,
    placements: Sequence[GlobalPlacement],
) -> None:
    fields = [field.name for field in dataclasses.fields(GlobalPlacement)]
    write_tsv(
        path,
        fields,
        (
            (
                *dataclasses.astuple(row)[:-4],
                json.dumps(row.query_blocks, separators=(",", ":")),
                json.dumps(row.target_blocks, separators=(",", ":")),
                json.dumps(row.paired_blocks, separators=(",", ":")),
                row.liftover_path,
            )
            for row in sorted(
                placements,
                key=lambda item: (item.query_id, item.target_record),
            )
        ),
    )


def atomic_copy(source: str, target: str) -> None:
    temporary = target + f".tmp.{os.getpid()}"
    try:
        shutil.copyfile(source, temporary)
        os.replace(temporary, target)
    finally:
        try:
            os.unlink(temporary)
        except FileNotFoundError:
            pass


def priority_input_records(
    queries: Sequence[NovelQuery],
) -> List[InputRecord]:
    records: List[InputRecord] = []
    for order, query in enumerate(queries):
        group = (
            query.haplotype.rsplit("_", 1)[0]
            if "_" in query.haplotype else query.haplotype
        )
        records.append(InputRecord(
            query.query_id,
            query.query_id,
            group,
            query.haplotype,
            order,
            query.local_end - query.local_start,
            query.source_contig,
            query.source_start,
            query.source_end,
        ))
    return records


def queries_to_pieces(
    queries: Sequence[NovelQuery],
    cache: Mapping[str, CachedCandidateSequence],
) -> List[Piece]:
    within_haplotype: Dict[str, int] = defaultdict(int)
    output: List[Piece] = []
    for index, query in enumerate(queries):
        rank = within_haplotype[query.haplotype]
        within_haplotype[query.haplotype] += 1
        output.append(Piece(
            temp_id=f"small_input_{index:012d}",
            sample=query.haplotype,
            contig=query.record_id,
            start=query.local_start,
            end=query.local_end,
            strand="+",
            seq=query_sequence(query, cache),
            origin_id=query.query_id,
            lift_id=query.query_id,
            within_sample_rank=rank,
        ))
    return output


def materialize_small_minset(
    pieces: Sequence[Piece],
    queries: Sequence[NovelQuery],
    priorities: Mapping[str, int],
    output_dir: str,
) -> List[SmallNovelRecord]:
    query_by_id = {query.query_id: query for query in queries}
    ordered = sorted(pieces, key=lambda piece: (
        -priorities.get(piece.sample, 0),
        -piece.length,
        piece.sample,
        piece.contig,
        piece.start,
        piece.end,
        piece.origin_id,
    ))
    records: List[SmallNovelRecord] = []
    for index, piece in enumerate(ordered, 1):
        query = query_by_id.get(piece.origin_id)
        if query is None:
            raise KeyError(
                f"self-clean piece has unknown origin {piece.origin_id!r}"
            )
        if (
            piece.contig != query.record_id
            or piece.start < query.local_start
            or piece.end > query.local_end
        ):
            raise ValueError(
                f"{piece.origin_id}: invalid retained piece "
                f"{piece.contig}:{piece.start}-{piece.end}"
            )
        source_start = query.source_start + piece.start - query.local_start
        source_end = source_start + piece.length
        records.append(SmallNovelRecord(
            f"smallnovel_{index:012d}",
            piece.sample,
            query.source_contig,
            source_start,
            source_end,
            query.partition,
            query.record_id,
            piece.start,
            piece.end,
            piece.origin_id,
            novelty_score(piece.seq, MASKED_BASE_WEIGHT),
            piece.length,
            piece.seq,
        ))

    fasta_path = os.path.join(output_dir, "smallnovel.fa")
    bed_path = os.path.join(output_dir, "smallnovel.bed")
    catalog_path = os.path.join(output_dir, "smallnovel.tsv")
    fasta_tmp = fasta_path + f".tmp.{os.getpid()}"
    bed_tmp = bed_path + f".tmp.{os.getpid()}"
    try:
        with open(fasta_tmp, "wt") as fasta, open(bed_tmp, "wt") as bed:
            bed.write(
                "#source_contig\tsource_start\tsource_end\t"
                "smallnovel_id\tnovelty_score\tstrand\thaplotype\t"
                "partition\torigin_query_id\trecord_id\tlocal_start\t"
                "local_end\n"
            )
            for record in records:
                fasta.write(
                    f">{record.smallnovel_id} "
                    f"source={record.haplotype}:{record.source_contig}:"
                    f"{record.source_start}-{record.source_end}+ "
                    f"assembly_coordinates={record.haplotype}:"
                    f"{record.source_contig}:{record.source_start}-"
                    f"{record.source_end}+ "
                    f"partition={record.partition} "
                    f"origin={record.origin_query_id} "
                    f"mapped_type=cleaned_small_novel "
                    f"mapped_coordinates={record.smallnovel_id}:0-"
                    f"{record.length}+\n"
                    f"{wrap_fasta(record.sequence)}\n"
                )
                bed.write(
                    f"{record.source_contig}\t{record.source_start}\t"
                    f"{record.source_end}\t{record.smallnovel_id}\t"
                    f"{record.novelty_score:.1f}\t+\t{record.haplotype}\t"
                    f"{record.partition}\t{record.origin_query_id}\t"
                    f"{record.record_id}\t{record.local_start}\t"
                    f"{record.local_end}\n"
                )
        os.replace(fasta_tmp, fasta_path)
        os.replace(bed_tmp, bed_path)
        try:
            os.unlink(fasta_path + ".fai")
        except FileNotFoundError:
            pass
    finally:
        for path in (fasta_tmp, bed_tmp):
            try:
                os.unlink(path)
            except FileNotFoundError:
                pass

    fields = [
        field.name for field in dataclasses.fields(SmallNovelRecord)
        if field.name != "sequence"
    ]
    write_tsv(
        catalog_path,
        fields,
        (
            tuple(
                getattr(record, field)
                for field in fields
            )
            for record in records
        ),
    )
    return records


def annotate_minset_placements(
    queries: Sequence[NovelQuery],
    placements: Dict[str, GlobalPlacement],
    cache: Mapping[str, CachedCandidateSequence],
    minset_fasta: str,
    output_path: str,
) -> Dict[str, str]:
    labels: Dict[str, str] = {}
    with IndexedFasta(minset_fasta) as target:
        for query in queries:
            placement = placements.get(query.query_id)
            if placement is None:
                continue
            path = build_graphic_liftover_path(
                query, placement, cache[query.query_id], target,
            )
            if path == ".":
                continue
            labels[query.query_id] = f"{placement.ref_haplotype}={path}"

    fields = (
        "query_id", "partition", "status", "target_kind",
        "target_record", "target_start", "target_end", "strand",
        "alignment_score", "aligned_score", "aligned_bases", "identity",
        "source", "query_blocks", "target_blocks", "paired_blocks",
        "liftover_path", "label",
    )
    rows = []
    for query in queries:
        placement = placements.get(query.query_id)
        if placement is None:
            rows.append((
                query.query_id, query.partition, "unmapped",
                ".", ".", ".", ".", ".", 0, 0, 0, 0,
                ".", "[]", "[]", "[]", ".", ".",
            ))
            continue
        rows.append((
            query.query_id,
            query.partition,
            "mapped",
            placement.ref_haplotype,
            placement.target_record,
            placement.ref_start,
            placement.ref_end,
            placement.strand,
            placement.alignment_score,
            placement.aligned_score,
            placement.aligned_bases,
            placement.identity,
            placement.source,
            json.dumps(placement.query_blocks, separators=(",", ":")),
            json.dumps(placement.target_blocks, separators=(",", ":")),
            json.dumps(placement.paired_blocks, separators=(",", ":")),
            labels.get(query.query_id, ".").split("=", 1)[-1],
            labels.get(query.query_id, "."),
        ))
    write_tsv(output_path, fields, rows)
    return labels


def parse_args(argv: Optional[Sequence[str]] = None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=(
            "Clean reference-uncovered global novelty against large novel "
            "loci and a priority-aware small self-minset"
        ),
    )
    parser.add_argument("-a", "--analysis-dir", required=True)
    parser.add_argument("-g", "--global-mapping-dir", required=True)
    parser.add_argument(
        "-l", "--local-mapping-dir",
        help=(
            "output of map_partition_local_novel.py; defaults to the "
            "partition_novelty directory beside --global-mapping-dir"
        ),
    )
    parser.add_argument("-n", "--large-novel-fasta", required=True)
    parser.add_argument("-q", "--query-paths", required=True)
    parser.add_argument("-o", "--output-dir", required=True)
    parser.add_argument("-j", "--jobs", type=int, default=16)
    parser.add_argument("-c", "--cores", type=int, default=16)
    parser.add_argument("-p", "--priority", default="CHM13_h1,HG38_h1,HG002")
    parser.add_argument("--alignment-anchor", type=int, default=50)
    parser.add_argument("--min-score", type=int, default=100)
    parser.add_argument("--min-identity", type=float, default=95.0)
    parser.add_argument("--minimap-query-batch", default="50M")
    parser.add_argument("--max-cycles", type=int, default=0)
    parser.add_argument("--self-align-min", type=int, default=50)
    parser.add_argument("--self-anchor", type=int, default=50)
    parser.add_argument("--self-blast-word-size", type=int, default=28)
    parser.add_argument("--blast-word-size", type=int, default=50)
    parser.add_argument("--blast-evalue", default="1e-300")
    parser.add_argument("--blast-max-target-seqs", type=int, default=100)
    parser.add_argument("--skip-blastn", action="store_true")
    parser.add_argument("--resume", action="store_true")
    parser.add_argument(
        "--log-level", choices=("DEBUG", "INFO", "WARNING", "ERROR"),
        default="INFO",
    )
    args = parser.parse_args(argv)
    for name in (
        "jobs", "cores", "min_score", "self_align_min",
        "self_blast_word_size", "blast_word_size", "blast_max_target_seqs",
    ):
        if getattr(args, name) < 1:
            parser.error(f"--{name.replace('_', '-')} must be positive")
    if args.alignment_anchor < 0 or args.self_anchor < 0 or args.max_cycles < 0:
        parser.error("anchor sizes and --max-cycles cannot be negative")
    if not 0 < args.min_identity <= 100:
        parser.error("--min-identity must be in (0,100]")
    if re.fullmatch(
        r"[1-9][0-9]*[kKmMgG]?", args.minimap_query_batch,
    ) is None:
        parser.error(
            "--minimap-query-batch must be a positive integer with an "
            "optional K/M/G suffix"
        )
    return args


def run(args: argparse.Namespace) -> None:
    analysis_dir = os.path.abspath(args.analysis_dir)
    global_dir = os.path.abspath(args.global_mapping_dir)
    local_dir = os.path.abspath(
        getattr(args, "local_mapping_dir", None)
        if getattr(args, "local_mapping_dir", None)
        else os.path.join(os.path.dirname(global_dir), "partition_novelty")
    )
    large_fasta = os.path.abspath(args.large_novel_fasta)
    query_paths = os.path.abspath(args.query_paths)
    output_dir = os.path.abspath(args.output_dir)
    Path(output_dir).mkdir(parents=True, exist_ok=True)

    partition_inputs = os.path.join(analysis_dir, "partition_inputs.tsv")
    global_queries_path = os.path.join(global_dir, "all_global_queries.tsv")
    global_placements_path = os.path.join(
        global_dir, "all_global_reference_placements.tsv",
    )
    global_reporting_path = os.path.join(
        global_dir, "all_global_reference_reporting_alignments.tsv",
    )
    global_marker = os.path.join(
        global_dir, "_GLOBAL_NOVEL_MAPPING_SUCCESS.json",
    )
    local_candidates_path = os.path.join(
        analysis_dir, "local_novel_candidates.tsv",
    )
    local_alignments_path = os.path.join(
        local_dir, "all_local_novel_alignments.tsv",
    )
    required = (
        partition_inputs,
        local_candidates_path,
        local_alignments_path,
        global_queries_path,
        global_placements_path,
        global_marker,
        large_fasta,
        query_paths,
    )
    for path in required:
        if not os.path.isfile(path):
            raise FileNotFoundError(path)
    if os.path.isfile(global_reporting_path):
        required = (*required, global_reporting_path)
    with open(global_marker, "rt") as handle:
        global_metadata = json.load(handle)
    if global_metadata.get("version") != (
        "global-novel-mapping-v5-graphic-annotation"
    ):
        raise ValueError(
            "global mapping predates annotation-only graphic liftover; "
            "rerun map_global_novel.py"
        )
    if float(global_metadata.get("min_score", -1)) != args.min_score:
        raise ValueError(
            f"--min-score {args.min_score} does not match global mapping "
            f"min_score {global_metadata.get('min_score')!r}"
        )
    main_reference_value = global_metadata.get("main_reference")
    if not main_reference_value:
        raise ValueError(
            "global mapping marker does not record main_reference"
        )
    main_reference = os.path.abspath(main_reference_value)
    if not os.path.isfile(main_reference):
        raise FileNotFoundError(main_reference)
    main_reference_records = annotation_reference_records(main_reference)

    assembly_sources = read_assembly_sources(query_paths)
    reference_sources = annotation_reference_sources(assembly_sources, main_reference)
    signature_files = [
        *required,
        main_reference,
        *(
            path
            for source in sorted(
                assembly_sources.values(), key=lambda row: row.haplotype,
            )
            for path in (source.fasta_path, source.fai_path)
        ),
    ]
    settings = {
        "annotation_reference_policy": ANNOTATION_POLICY,
        "priority": args.priority,
        "alignment_anchor": args.alignment_anchor,
        "min_score": args.min_score,
        "min_identity": args.min_identity,
        "minimap_query_batch": args.minimap_query_batch,
        "max_cycles": args.max_cycles,
        "self_align_min": args.self_align_min,
        "self_anchor": args.self_anchor,
        "self_blast_word_size": args.self_blast_word_size,
        "blast_word_size": args.blast_word_size,
        "blast_evalue": args.blast_evalue,
        "blast_max_target_seqs": args.blast_max_target_seqs,
        "skip_blastn": args.skip_blastn,
    }
    signature_file_fingerprints = [
        fingerprint(path) for path in signature_files
    ]
    signature = hashlib.sha256(json.dumps({
        "version": VERSION,
        "files": signature_file_fingerprints,
        "settings": settings,
    }, sort_keys=True).encode()).hexdigest()
    upstream_alignment_signature = signature
    marker = os.path.join(output_dir, "_SMALL_NOVEL_MINSET_SUCCESS.json")
    audit_dir = os.path.join(output_dir, "intermediates")
    audit_stems = (
        "00_all_global_novel_queries",
        "01_reference_exclusion_input",
        "02_reference_aligned",
        "03_reference_remaining",
        "04_reference_filtered_small",
        "05_large_novel_aligned",
        "06_large_novel_remaining",
        "07_large_novel_filtered_small",
        "08_self_clean_retained",
        "09_self_clean_removed",
    )
    expected = tuple(os.path.join(output_dir, name) for name in (
        "smallnovel.fa", "smallnovel.bed", "smallnovel.tsv",
        "novel_minset.fa", "all_small_novel_annotations.tsv",
        "cleaned_small_novel_annotations.tsv",
        "all_small_novel_alignments.tsv",
        "small_novel_alignment_summary.tsv",
        "all_small_novel_queries.annotated.fa",
        "all_small_novel_reporting_alignments.tsv",
    )) + tuple(
        os.path.join(audit_dir, stem + extension)
        for stem in audit_stems
        for extension in (".fa", ".bed", ".tsv")
    ) + (
        os.path.join(audit_dir, "02_reference_aligned.placements.tsv"),
        os.path.join(audit_dir, "02_reference_reporting_alignments.tsv"),
        os.path.join(audit_dir, "05_large_novel_aligned.placements.tsv"),
        os.path.join(audit_dir, "05_large_novel_reporting_alignments.tsv"),
        os.path.join(audit_dir, "08_self_clean_mappings.placements.tsv"),
        os.path.join(audit_dir, "08_cleaned_small_reporting_alignments.tsv"),
        os.path.join(audit_dir, "local_novel_candidates.tsv"),
        os.path.join(audit_dir, "local_novel_alignments.tsv"),
        os.path.join(audit_dir, "manifest.tsv"),
        os.path.join(
            audit_dir, "query_sequence_cache", "candidate_sequences.fa",
        ),
        os.path.join(
            audit_dir, "query_sequence_cache",
            "candidate_sequences.index.tsv",
        ),
        os.path.join(
            audit_dir, "query_sequence_cache", "_SUCCESS.json",
        ),
        os.path.join(audit_dir, "_ALIGNMENT_SUMMARY_SUCCESS.json"),
    )
    if args.resume and os.path.isfile(marker):
        try:
            with open(marker, "rt") as handle:
                saved = json.load(handle)
            if (
                saved.get("version") == VERSION
                and saved.get("signature") == signature
                and saved.get("reporting_resolution")
                == REPORTING_RESOLUTION_VERSION
                and all(os.path.isfile(path) for path in expected)
            ):
                LOG.info("RESUME: small-novel minset is complete in %s", output_dir)
                return
        except (OSError, ValueError, AttributeError):
            pass

    queries = read_novel_queries(global_queries_path)
    reference_placements = read_global_placements(global_placements_path)
    if os.path.isfile(global_reporting_path):
        initial_reference_reporting = read_reporting_components((
            global_reporting_path,
        ))
    else:
        LOG.warning(
            "Global multi-locus report is absent; using legacy selected "
            "placements. Replay map_global_novel.py "
            "--reuse-alignments-only before production reporting."
        )
        initial_reference_reporting = reporting_components_from_placements(
            reference_placements,
        )
    unknown = set(reference_placements) - {query.query_id for query in queries}
    if unknown:
        raise ValueError(
            "reference placements contain unknown global queries: "
            + ",".join(sorted(unknown)[:20])
        )
    partitions = read_partition_inputs(partition_inputs)
    candidates = global_mapping_candidates(queries)
    source_records = load_candidate_header_records(
        candidates, partitions, args.jobs,
    )
    cache_digest = hashlib.sha256()
    cache_digest.update(b"small-novel-query-cache-v2-independent\n")
    cache_digest.update(json.dumps({
        "files": [
            fingerprint(path)
            for path in (
                global_queries_path, partition_inputs, query_paths,
                *(
                    path
                    for source in sorted(
                        assembly_sources.values(),
                        key=lambda row: row.haplotype,
                    )
                    for path in (source.fasta_path, source.fai_path)
                ),
            )
        ],
        "alignment_anchor": args.alignment_anchor,
    }, sort_keys=True).encode())
    for candidate in candidates:
        record = source_records[(
            candidate.novel_partition, candidate.record_id,
        )]
        cache_digest.update((
            f"\n{candidate.candidate_id}\t{candidate.haplotype}\t"
            f"{candidate.source_contig}\t{candidate.source_start}\t"
            f"{candidate.source_end}\t{candidate.novel_partition}\t"
            f"{candidate.record_id}\t{candidate.local_start}\t"
            f"{candidate.local_end}\t{record.record_length}"
        ).encode())
    cache_signature = cache_digest.hexdigest()
    cache_dir = os.path.join(audit_dir, "query_sequence_cache")
    cached = build_disk_query_cache(
        candidates,
        source_records,
        assembly_sources,
        args.alignment_anchor,
        args.jobs,
        cache_dir,
        cache_signature,
        args.resume,
    )
    del source_records
    succeeded = False
    try:
        reference_input = global_residual_queries(
            queries, reference_placements, args.min_score, cached,
        )
        initial_reference_residual_count = len(reference_input)
        reference_residuals = list(reference_input)
        completion_placements: Dict[str, GlobalPlacement] = {}
        completion_reporting: List[ResolvedPlacementComponent] = []
        if reference_residuals and main_reference_records:
            completion_cache = ResidualQueryCache(
                cached, [query.query_id for query in reference_residuals],
            )
            completion_placements = align_global_queries(
                reference_residuals,
                global_mapping_candidates(reference_residuals),
                completion_cache,
                main_reference,
                main_reference_records,
                os.path.join(
                    output_dir, "alignment_work", "reference_completion",
                ),
                args,
                upstream_alignment_signature + ":reference-completion",
            )
            completion_reporting = annotate_global_reporting_components(
                reference_residuals,
                read_stage_reporting_or_placements(os.path.join(
                    output_dir, "alignment_work", "reference_completion",
                    "partitions", "global_reference",
                    "reporting_placements.tsv",
                ), completion_placements),
                completion_cache,
                reference_sources,
            )
            reference_residuals = global_residual_queries(
                reference_residuals,
                completion_placements,
                args.min_score,
                cached,
            )
        LOG.info(
            "Reference completion: %d residual query pieces from conservative "
            "global liftover -> %d score-passing reference-uncovered pieces",
            initial_reference_residual_count, len(reference_residuals),
        )
        large_records = fasta_target_records(large_fasta)
        large_placements: Dict[str, GlobalPlacement] = {}
        large_reporting: List[ResolvedPlacementComponent] = []
        residual_cache = ResidualQueryCache(
            cached, [query.query_id for query in reference_residuals],
        )
        if reference_residuals and large_records:
            large_placements = align_global_queries(
                reference_residuals,
                global_mapping_candidates(reference_residuals),
                residual_cache,
                large_fasta,
                large_records,
                os.path.join(output_dir, "alignment_work", "large_exclusion"),
                args,
                upstream_alignment_signature + ":large-exclusion",
            )
            large_reporting = annotate_fixed_target_reporting_components(
                reference_residuals,
                read_stage_reporting_or_placements(os.path.join(
                    output_dir, "alignment_work", "large_exclusion",
                    "partitions", "global_reference",
                    "reporting_placements.tsv",
                ), large_placements),
                residual_cache,
                large_fasta,
            )
        small_candidates = global_residual_queries(
            reference_residuals,
            large_placements,
            args.min_score,
            cached,
        )
        LOG.info(
            "Small-novel construction: %d global queries -> %d "
            "reference residuals -> %d large-locus residuals",
            len(queries), len(reference_residuals), len(small_candidates),
        )

        priority_records = priority_input_records(small_candidates)
        priorities, unmatched = build_haplotype_priorities(
            priority_records, args.priority,
        )
        if unmatched:
            LOG.warning(
                "Priority selectors matched no small-novel haplotype: %s",
                ",".join(unmatched),
            )
        pieces = queries_to_pieces(small_candidates, cached)
        self_dir = os.path.join(output_dir, "alignment_work", "self_clean")
        if os.path.isdir(self_dir) and not args.resume:
            shutil.rmtree(self_dir)
        self_args = SimpleNamespace(
            cores=args.cores,
            min_identity=args.min_identity,
            min_unmasked=args.min_score,
            self_align_min=args.self_align_min,
            self_anchor=args.self_anchor,
            self_blast_word_size=args.self_blast_word_size,
            blast_word_size=args.blast_word_size,
            blast_evalue=args.blast_evalue,
            blast_max_target_seqs=args.blast_max_target_seqs,
            skip_blastn=args.skip_blastn,
            max_cycles=args.max_cycles,
            resume=args.resume,
        )
        retained = self_clean_two_stage(
            pieces,
            self_dir,
            self_args,
            priorities,
            "Small-novel minset",
        )

        # Stable audit products. Unlike the aligner cycle FASTAs, these use
        # absolute assembly coordinates, preserve provenance, exclude
        # temporary anchors, and explicitly materialize retained/removed sets.
        Path(audit_dir).mkdir(parents=True, exist_ok=True)
        self_retained = retained_piece_queries(retained, small_candidates)
        self_removed = subtract_query_segments(
            small_candidates, self_retained,
        )
        global_reference_aligned = placement_query_segments(
            queries, reference_placements,
        )
        completion_reference_aligned = placement_query_segments(
            reference_input, completion_placements,
        )
        reference_aligned = normalized_query_segments(
            queries,
            [*global_reference_aligned, *completion_reference_aligned],
        )
        reference_remaining = normalized_query_segments(
            queries, reference_residuals,
        )
        reference_unaligned = subtract_query_segments(
            queries, reference_aligned,
        )
        reference_filtered = subtract_query_segments(
            reference_unaligned, reference_remaining,
        )
        large_aligned = normalized_query_segments(
            reference_residuals,
            placement_query_segments(reference_residuals, large_placements),
        )
        large_remaining = normalized_query_segments(
            reference_residuals, small_candidates,
        )
        large_unaligned = subtract_query_segments(
            reference_residuals, large_aligned,
        )
        large_filtered = subtract_query_segments(
            large_unaligned, large_remaining,
        )
        audit_sets = (
            ("00_all_global_novel_queries", queries),
            ("01_reference_exclusion_input", reference_input),
            ("02_reference_aligned", reference_aligned),
            ("03_reference_remaining", reference_remaining),
            ("04_reference_filtered_small", reference_filtered),
            ("05_large_novel_aligned", large_aligned),
            ("06_large_novel_remaining", large_remaining),
            ("07_large_novel_filtered_small", large_filtered),
            ("08_self_clean_retained", self_retained),
            ("09_self_clean_removed", self_removed),
        )
        for label, members in audit_sets:
            write_query_audit_bundle(
                os.path.join(audit_dir, label), label, members, cached,
            )
        write_placement_audit(
            os.path.join(
                audit_dir, "02_reference_aligned.placements.tsv",
            ),
            [*reference_placements.values(), *completion_placements.values()],
        )
        write_reporting_components(
            os.path.join(
                audit_dir, "02_reference_reporting_alignments.tsv",
            ),
            [*initial_reference_reporting, *completion_reporting],
        )
        write_placement_audit(
            os.path.join(
                audit_dir, "05_large_novel_aligned.placements.tsv",
            ),
            list(large_placements.values()),
        )
        write_reporting_components(
            os.path.join(
                audit_dir, "05_large_novel_reporting_alignments.tsv",
            ),
            large_reporting,
        )
        atomic_copy(
            local_candidates_path,
            os.path.join(audit_dir, "local_novel_candidates.tsv"),
        )
        atomic_copy(
            local_alignments_path,
            os.path.join(audit_dir, "local_novel_alignments.tsv"),
        )
        manifest_rows = (
            ("00", "all_global_novel_queries", "all score-passing intact global queries"),
            ("01", "reference_exclusion_input", "residual/skipped pool entering reference completion"),
            ("02", "reference_aligned", "query intervals represented by the main reference"),
            ("03", "reference_remaining", "score-passing query intervals remaining after reference exclusion"),
            ("04", "reference_filtered_small", "unaligned reference residuals removed by score/N filtering"),
            ("05", "large_novel_aligned", "query intervals represented by established large novel loci"),
            ("06", "large_novel_remaining", "score-passing intervals entering self-clean"),
            ("07", "large_novel_filtered_small", "unaligned large-locus residuals removed by score/N filtering"),
            ("08", "self_clean_retained", "representatives retained by priority-aware self-clean"),
            ("09", "self_clean_removed", "query intervals represented by the retained self-clean set"),
        )
        manifest_output = [
            (
                stage,
                name,
                f"{stage}_{name}.fa",
                f"{stage}_{name}.bed",
                f"{stage}_{name}.tsv",
                description,
            )
            for stage, name, description in manifest_rows
        ]
        manifest_output.append((
            "cache",
            "query_sequence_cache",
            "query_sequence_cache/candidate_sequences.fa",
            ".",
            "query_sequence_cache/candidate_sequences.index.tsv",
            "persistent disk-backed source-window cache reused on resume",
        ))
        manifest_output.extend((
            (
                "alignment",
                "all_small_novel_alignments",
                ".",
                ".",
                "../all_small_novel_alignments.tsv",
                "accepted reference, large-novel, and cleaned-small "
                "alignments with exact query/target blocks",
            ),
            (
                "alignment",
                "small_novel_alignment_summary",
                "../all_small_novel_queries.annotated.fa",
                ".",
                "../small_novel_alignment_summary.tsv",
                "one organized assignment summary per original query",
            ),
            (
                "alignment",
                "cleaned_small_novel_annotations",
                ".",
                ".",
                "../cleaned_small_novel_annotations.tsv",
                "fixed-target mapping of every stage-06 residual to its "
                "final cleaned-small representative",
            ),
        ))
        write_tsv(
            os.path.join(audit_dir, "manifest.tsv"),
            ("stage", "name", "fasta", "bed", "tsv", "description"),
            manifest_output,
        )

        small_records = materialize_small_minset(
            retained, small_candidates, priorities, output_dir,
        )
        small_fasta = os.path.join(output_dir, "smallnovel.fa")
        cleaned_target_records = fasta_target_records(
            small_fasta,
            [record.smallnovel_id for record in small_records],
        )
        cleaned_placements: Dict[str, GlobalPlacement] = {}
        cleaned_reporting: List[ResolvedPlacementComponent] = []
        cleaned_cache = ResidualQueryCache(
            cached, [query.query_id for query in small_candidates],
        )
        if small_candidates and cleaned_target_records:
            cleaned_placements = align_global_queries(
                small_candidates,
                global_mapping_candidates(small_candidates),
                cleaned_cache,
                small_fasta,
                cleaned_target_records,
                os.path.join(
                    output_dir, "alignment_work", "cleaned_small_novel",
                ),
                args,
                signature + ":cleaned-small-novel",
            )
            cleaned_reporting = annotate_fixed_target_reporting_components(
                small_candidates,
                read_stage_reporting_or_placements(os.path.join(
                    output_dir, "alignment_work", "cleaned_small_novel",
                    "partitions", "global_reference",
                    "reporting_placements.tsv",
                ), cleaned_placements),
                cleaned_cache,
                small_fasta,
            )
        missing_cleaned_assignments = {
            query.query_id for query in small_candidates
        } - set(cleaned_placements)
        if missing_cleaned_assignments:
            raise ValueError(
                f"{len(missing_cleaned_assignments)} score-passing "
                "self-clean inputs have no final cleaned-small placement; "
                "first: "
                + ",".join(sorted(missing_cleaned_assignments)[:20])
            )
        write_placement_audit(
            os.path.join(
                audit_dir, "08_self_clean_mappings.placements.tsv",
            ),
            list(cleaned_placements.values()),
        )
        write_reporting_components(
            os.path.join(
                audit_dir, "08_cleaned_small_reporting_alignments.tsv",
            ),
            cleaned_reporting,
        )
        annotate_minset_placements(
            small_candidates,
            cleaned_placements,
            cleaned_cache,
            small_fasta,
            os.path.join(
                output_dir, "cleaned_small_novel_annotations.tsv",
            ),
        )
        combined_reporting = combine_small_novel_reporting_components(
            queries,
            (
                (
                    "reference",
                    [*initial_reference_reporting, *completion_reporting],
                ),
                ("large_novel", large_reporting),
                ("cleaned_small_novel", cleaned_reporting),
            ),
            cached,
            reference_sources,
            large_fasta,
            small_fasta,
        )
        write_reporting_components(
            os.path.join(
                output_dir, "all_small_novel_reporting_alignments.tsv",
            ),
            combined_reporting,
        )
        large_ids = set(large_records)
        duplicate_ids = large_ids & {
            record.smallnovel_id for record in small_records
        }
        if duplicate_ids:
            raise ValueError(
                "large/small novel FASTAs have duplicate IDs: "
                + ",".join(sorted(duplicate_ids)[:20])
            )
        minset_fasta = os.path.join(output_dir, "novel_minset.fa")
        concatenate_fastas(
            (large_fasta, small_fasta),
            minset_fasta,
        )
        minset_records = fasta_target_records(
            minset_fasta,
            [record.smallnovel_id for record in small_records],
        )
        final_placements: Dict[str, GlobalPlacement] = {}
        if queries and minset_records:
            final_placements = align_global_queries(
                queries,
                candidates,
                cached,
                minset_fasta,
                minset_records,
                os.path.join(output_dir, "alignment_work", "final_minset"),
                args,
                signature + ":final-minset",
            )
        annotations_path = os.path.join(
            output_dir, "all_small_novel_annotations.tsv",
        )
        annotate_minset_placements(
            queries, final_placements, cached, minset_fasta, annotations_path,
        )
        organize_mapping(output_dir, require_all=True)

        temporary = marker + f".tmp.{os.getpid()}"
        with open(temporary, "wt") as output:
            json.dump({
                "version": VERSION,
                "signature": signature,
                "global_queries": len(queries),
                "initial_reference_residuals": (
                    initial_reference_residual_count
                ),
                "reference_completion_placements": len(
                    completion_placements
                ),
                "reference_residuals": len(reference_residuals),
                "large_locus_residuals": len(small_candidates),
                "large_novel_records": len(large_records),
                "small_novel_records": len(small_records),
                "self_clean_records": len(retained),
                "cleaned_small_placements": len(cleaned_placements),
                "reporting_resolution": REPORTING_RESOLUTION_VERSION,
                "reporting_components": len(combined_reporting),
                "annotated_queries": len(final_placements),
                "large_novel_fasta": large_fasta,
                "query_paths": query_paths,
                "min_score": args.min_score,
            }, output, sort_keys=True)
            output.write("\n")
        os.replace(temporary, marker)
        succeeded = True
        LOG.info(
            "Small-novel minset complete: %d fixed large loci + %d selected "
            "small representatives; annotated %d/%d original queries",
            len(large_records), len(small_records),
            len(final_placements), len(queries),
        )
    finally:
        cached.close()


def main(argv: Optional[Sequence[str]] = None) -> int:
    args = parse_args(argv)
    logging.basicConfig(
        level=getattr(logging, args.log_level),
        format="[%(asctime)s] %(levelname)s: %(message)s",
        datefmt="%Y-%m-%d %H:%M:%S",
    )
    try:
        run(args)
    except Exception as error:
        LOG.error("%s", error)
        if args.log_level == "DEBUG":
            raise
        return 1
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
