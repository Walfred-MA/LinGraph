#!/usr/bin/env python3
"""Clean redundant sequence from one small homologous-sequence FASTA.

All records from the highest-priority selected haplotype form one fixed
reference set. Every other record is cleaned against that combined reference,
and the resulting residual pieces are self-cleaned together. Each phase runs
minimap2 to convergence, followed by Winnowmap plus an independent full-query
BLASTN pass to convergence.

Output coordinates are 0-based, half-open, and local to the input FASTA
records.  The final retained core intervals are expanded by the requested
anchor size and clipped to their original record lengths.
"""
from __future__ import annotations

import argparse
import dataclasses
import filecmp
import hashlib
import json
import logging
import os
import re
import shutil
import subprocess
import sys
import tempfile
from collections import defaultdict
from pathlib import Path
from typing import Dict, List, Optional, Sequence, Set, Tuple

from minsetref_align import (
    AlignmentHit,
    build_minimap2_index,
    collect_alignment_hits,
    coverage_by_query_id,
    ensure_executable,
)
from minsetref_core import (
    IndexedFasta,
    Piece,
    iter_n_free_segments,
    wrap_fasta,
    write_fasta_pieces,
)
from minsetref_segments import find_unaligned_segments, merge_intervals, novelty_score
from minsetref_selfclean import (
    _pieces_signature,
    _read_cycle_checkpoint,
    _write_cycle_checkpoint,
    self_clean_converge,
)


LOG = logging.getLogger("minsetref_light")
PRIOR_SCORE_STEP = 10_000_000_000
MASKED_BASE_WEIGHT = 0.3
MINIMAP_PARAMS = dict(
    preset="asm5", p=0.001, N=20, f=0.001, K="100M",
    retry_on_sigkill=True, retry_threads=32, retry_K="100M",
)
LIGHT_WINNOW_PARAMS = dict(k=19, w=10, m=100, p=0.001, N=100)
GENOMIC_COORD_RE = re.compile(r"^(.+):(\d+)-(\d+)(?:[+-])?$")
COORDINATE_REFERENCE_NAME_RE = re.compile(r"^(.+)_(\d+)_(\d+)$")
AUTO_FAIDX_MIN_BYTES = 100 * 1024 * 1024
QUERY_REFERENCE_MIN_SCORE = 500
PIECE_SOURCE_RE = re.compile(r"^([^:]+):(.*):(\d+)-(\d+)([+-])$")
KMER_FRAGMENT_RE = re.compile(r"^(.*):(\d+)-(\d+)$")


@dataclasses.dataclass(frozen=True)
class InputRecord:
    record_id: str
    description: str
    group: str
    haplotype: str
    order: int
    length: int
    genomic_contig: Optional[str] = None
    genomic_start: Optional[int] = None
    genomic_end: Optional[int] = None


@dataclasses.dataclass(frozen=True)
class AnchoredCore:
    """Core coordinates within one temporary anchored self-alignment record."""
    core_start: int
    core_end: int
    query_length: int


def _parse_genomic_token(token: str) -> Optional[Tuple[str, int, int]]:
    match = GENOMIC_COORD_RE.match(token)
    if match is None:
        return None
    contig, start_text, end_text = match.groups()
    start, end = int(start_text), int(end_text)
    if end <= start:
        raise ValueError(f"invalid genomic interval {token!r}: end must exceed start")
    return contig, start, end


def read_input_records(path: str) -> List[InputRecord]:
    """Read record identity/header metadata and indexed sequence lengths."""
    if str(path).endswith(".gz"):
        raise ValueError(
            "minsetref_light requires an uncompressed FASTA; decompress the input first"
        )
    descriptions: List[Tuple[str, str, str, str, Optional[Tuple[str, int, int]]]] = []
    seen: Dict[str, int] = {}
    with open(path, "rt") as handle:
        for line_number, raw in enumerate(handle, 1):
            if not raw.startswith(">"):
                continue
            description = raw[1:].strip()
            fields = description.split()
            if not fields:
                raise ValueError(f"empty FASTA header at {path}:{line_number}")
            record_id = fields[0]
            if record_id in seen:
                raise ValueError(
                    f"duplicate FASTA record id {record_id!r} at {path}:{line_number}; "
                    f"first seen at header {seen[record_id]}"
                )
            seen[record_id] = line_number
            id_fields = record_id.split("_")
            if len(id_fields) < 3:
                raise ValueError(
                    f"record id {record_id!r} must contain at least three '_' fields "
                    "(for example g29647_CHM13_h1_1)"
                )
            group = id_fields[1]
            haplotype = "_".join(id_fields[1:3])
            genomic = None
            for token in fields[1:]:
                genomic = _parse_genomic_token(token)
                if genomic is not None:
                    break
            descriptions.append((record_id, description, group, haplotype, genomic))
    if not descriptions:
        raise ValueError(f"input FASTA contains no records: {path}")

    with IndexedFasta(path) as fasta:
        indexed_names = fasta.names()
        expected_names = [fields[0] for fields in descriptions]
        if indexed_names != expected_names:
            raise ValueError(
                "FASTA headers and random-access index disagree; remove a stale .fai "
                "or run samtools faidx again"
            )
        records: List[InputRecord] = []
        for order, (record_id, description, group, haplotype, genomic) in enumerate(
            descriptions
        ):
            kwargs = {}
            if genomic is not None:
                kwargs = {
                    "genomic_contig": genomic[0],
                    "genomic_start": genomic[1],
                    "genomic_end": genomic[2],
                }
            records.append(
                InputRecord(
                    record_id=record_id,
                    description=description,
                    group=group,
                    haplotype=haplotype,
                    order=order,
                    length=fasta.length(record_id),
                    **kwargs,
                )
            )
    return records


def ensure_large_input_fai(
    fasta_path: str,
    minimum_bytes: int = AUTO_FAIDX_MIN_BYTES,
) -> bool:
    """Create or refresh ``FASTA.fai`` for inputs of at least 100 MiB.

    Returns true when samtools was run. A current nonempty index is reused.
    """
    fasta_size = os.path.getsize(fasta_path)
    if fasta_size <= minimum_bytes:
        return False
    fai_path = fasta_path + ".fai"
    index_current = (
        os.path.isfile(fai_path)
        and os.path.getsize(fai_path) > 0
        and os.stat(fai_path).st_mtime_ns >= os.stat(fasta_path).st_mtime_ns
    )
    if index_current:
        LOG.info(
            "Reusing current FASTA index for %.2f GiB input: %s",
            fasta_size / (1024 ** 3), fai_path,
        )
        return False
    ensure_executable("samtools")
    LOG.info(
        "Input is %.2f GiB; building FASTA index with samtools faidx: %s",
        fasta_size / (1024 ** 3), fasta_path,
    )
    subprocess.run(["samtools", "faidx", fasta_path], check=True)
    if not os.path.isfile(fai_path) or os.path.getsize(fai_path) == 0:
        raise RuntimeError(f"samtools faidx did not create a nonempty {fai_path}")
    return True


def _parse_coordinate_reference_name(
    name: str, context: Optional[str] = None,
) -> Optional[Tuple[str, int, int]]:
    """Parse HAPLOTYPE_CONTIG_SOURCESTART_SOURCEEND, with strict spans."""
    match = COORDINATE_REFERENCE_NAME_RE.fullmatch(name)
    if match is None:
        return None
    source_name, start_text, end_text = match.groups()
    source_start, source_end = int(start_text), int(end_text)
    if source_name.count("_") < 2 or source_end <= source_start:
        location = f"{context}: " if context else ""
        raise ValueError(
            f"{location}invalid coordinate reference name {name!r}; expected "
            "HAPLOTYPE_CONTIG_SOURCESTART_SOURCEEND with SOURCEEND > SOURCESTART"
        )
    return source_name, source_start, source_end


def _translate_coordinate_reference_interval(
    name: str,
    start: int,
    end: int,
    context: Optional[str] = None,
    coordinate_space: Optional[str] = None,
) -> Optional[Tuple[str, int, int, int, int]]:
    """Translate a coordinate-named BED interval to source and local space.

    New local-block BEDs use absolute source coordinates in columns 2-3.  The
    legacy representation used 0-based coordinates local to the encoded
    interval.  Exact source bounds identify the new representation; otherwise
    the legacy local form is preferred where the two numeric spaces overlap.
    """
    parsed = _parse_coordinate_reference_name(name, context)
    if parsed is None:
        return None
    source_name, source_start, source_end = parsed
    span = source_end - source_start
    location = f"{context}: " if context else ""
    if coordinate_space not in {None, "source"}:
        raise ValueError(
            f"{location}unsupported coordinate space {coordinate_space!r}"
        )
    if coordinate_space == "source":
        if source_start <= start < end <= source_end:
            return (
                source_name,
                start,
                end,
                start - source_start,
                end - source_start,
            )
        raise ValueError(
            f"{location}source-coordinate BED interval {start}-{end} is "
            f"outside {source_start}-{source_end} for {name!r}"
        )
    if start == source_start and end == source_end:
        return source_name, start, end, 0, span
    if 0 <= start < end <= span:
        return (
            source_name,
            source_start + start,
            source_start + end,
            start,
            end,
        )
    if source_start <= start < end <= source_end:
        return (
            source_name,
            start,
            end,
            start - source_start,
            end - source_start,
        )
    raise ValueError(
        f"{location}BED interval {start}-{end} is outside coordinate-named "
        f"reference length {span} and source interval "
        f"{source_start}-{source_end} for {name!r}"
    )


def _selectors_from_file(path: str) -> List[str]:
    selectors: List[str] = []
    with open(path, "rt") as handle:
        for line_number, raw in enumerate(handle, 1):
            content = raw.strip()
            if not content or content.startswith("#"):
                continue
            fields = content.split()
            if len(fields) >= 3:
                try:
                    start, end = int(fields[1]), int(fields[2])
                except ValueError:
                    pass
                else:
                    # Current local-block BED6+ metadata:
                    # CONTIG CORE_START CORE_END BLOCK SCORE STRAND SAMPLE
                    # FASTA_ID FASTA_START FASTA_END source_coordinates
                    if (
                        len(fields) >= 11
                        and fields[10] == "source_coordinates"
                    ):
                        try:
                            fasta_start, fasta_end = (
                                int(fields[8]), int(fields[9])
                            )
                        except ValueError as error:
                            raise ValueError(
                                f"{path}:{line_number}: invalid anchored FASTA "
                                "coordinates"
                            ) from error
                        if not (
                            fasta_start <= start < end <= fasta_end
                            and fields[0]
                            and fields[6]
                            and fields[7]
                        ):
                            raise ValueError(
                                f"{path}:{line_number}: core {start}-{end} is "
                                f"outside anchored FASTA interval "
                                f"{fasta_start}-{fasta_end}"
                            )
                        selectors.append(fields[7])
                        continue
                    numbered_block = False
                    if len(fields) >= 4:
                        try:
                            numbered_block = int(fields[3]) > 0
                        except ValueError:
                            pass
                    coordinate_space = (
                        "source"
                        if len(fields) >= 7
                        and fields[6] == "source_coordinates"
                        else None
                    )
                    if numbered_block or coordinate_space == "source":
                        translated = _translate_coordinate_reference_interval(
                            fields[0], start, end, f"{path}:{line_number}",
                            coordinate_space,
                        )
                        if translated is None:
                            bed_kind = (
                                "numbered block BED"
                                if numbered_block
                                else "source-coordinate block BED"
                            )
                            raise ValueError(
                                f"{path}:{line_number}: {bed_kind} "
                                f"contig {fields[0]!r} is not coordinate-named; "
                                "expected "
                                "HAPLOTYPE_CONTIG_SOURCESTART_SOURCEEND. "
                                "Run temp_fixes/"
                                "fix_existing_block_bed_names.py first."
                            )
                        if coordinate_space == "source":
                            assert translated is not None
                            start, end = translated[3], translated[4]
                    selectors.append(f"{fields[0]}:{start}-{end}")
                    continue
            for field in fields:
                selectors.extend(token.strip() for token in field.split(",") if token.strip())
    return selectors


def expand_reference_selectors(values: Sequence[str]) -> List[str]:
    """Expand repeatable -r values, including selector and BED-like files."""
    selectors: List[str] = []
    for value in values:
        if os.path.isfile(value):
            selectors.extend(_selectors_from_file(value))
        else:
            selectors.extend(token.strip() for token in value.split(",") if token.strip())
    return selectors


def _name_embeds_haplotype(name: str, haplotype: str) -> bool:
    """Whether a synthetic contig name contains a complete haplotype token."""
    return (
        name == haplotype
        or name.startswith(haplotype + "_")
        or name.endswith("_" + haplotype)
        or f"_{haplotype}_" in name
    )


def _coordinate_named_reference_overlap(
    record: InputRecord, name: str, start: int, end: int,
) -> Optional[bool]:
    """Match a BED interval on SAMPLE_CONTIG_SOURCESTART_SOURCEEND."""
    if (
        record.genomic_contig is None
        or record.genomic_start is None
        or record.genomic_end is None
    ):
        return None
    translated = _translate_coordinate_reference_interval(name, start, end)
    if translated is None:
        return None
    source_name, translated_start, translated_end, _local_start, _local_end = (
        translated
    )
    if source_name != f"{record.haplotype}_{record.genomic_contig}":
        return None
    return (
        record.genomic_start < translated_end
        and translated_start < record.genomic_end
    )


def _selector_matches(record: InputRecord, selector: str) -> bool:
    genomic = _parse_genomic_token(selector)
    if genomic is not None:
        contig, start, end = genomic
        coordinate_named_match = _coordinate_named_reference_overlap(
            record, contig, start, end,
        )
        if coordinate_named_match is not None:
            return coordinate_named_match
        exact_genomic_match = (
            record.genomic_contig == contig
            and record.genomic_start is not None
            and record.genomic_end is not None
            and record.genomic_start < end
            and start < record.genomic_end
        )
        if exact_genomic_match:
            return True
        # Local-block reference BEDs can use synthetic names such as
        # main_CHM13_h1_NC_060927.1_merged_000206. Identify the corresponding
        # haplotype and source contig embedded in that name, then require the
        # same half-open overlap used for an ordinary genomic BED contig.
        if _name_embeds_haplotype(contig, record.haplotype):
            return (
                record.genomic_contig is not None
                and record.genomic_contig in contig
                and record.genomic_start is not None
                and record.genomic_end is not None
                and record.genomic_start < end
                and start < record.genomic_end
            )
        return False
    if selector == record.record_id:
        return True
    if _name_embeds_haplotype(selector, record.haplotype):
        return True
    if "_" in selector:
        return record.haplotype == selector
    return record.group == selector


def select_reference_ids(
    records: Sequence[InputRecord],
    selectors: Sequence[str],
    priorities: Optional[Dict[str, int]] = None,
) -> Tuple[Set[str], List[str], bool]:
    """Select all records from one highest-priority reference haplotype.

    Matching ``-r`` selectors restrict which haplotypes may win, but do not
    restrict the winning haplotype to only the overlapping record: every one of
    its records is required because a partition can contain multiple sequences
    from the same sample. If no selector matches, every input haplotype is
    eligible and the priority/length fallback chooses one.
    """
    priorities = priorities or {}
    selector_matches: Set[str] = set()
    unmatched: List[str] = []
    for selector in selectors:
        matches = [record.record_id for record in records if _selector_matches(record, selector)]
        if matches:
            selector_matches.update(matches)
        else:
            unmatched.append(selector)
    fallback = not selector_matches
    eligible_haplotypes = (
        {record.haplotype for record in records if record.record_id in selector_matches}
        if selector_matches
        else {record.haplotype for record in records}
    )
    by_haplotype: Dict[str, List[InputRecord]] = defaultdict(list)
    for record in records:
        if record.haplotype in eligible_haplotypes:
            by_haplotype[record.haplotype].append(record)
    winning_haplotype = max(
        by_haplotype,
        key=lambda haplotype: (
            priorities.get(haplotype, 0)
            + max(record.length for record in by_haplotype[haplotype]),
            sum(record.length for record in by_haplotype[haplotype]),
            -min(record.order for record in by_haplotype[haplotype]),
        ),
    )
    selected = {
        record.record_id for record in records
        if record.haplotype == winning_haplotype
    }
    return selected, unmatched, fallback


def build_haplotype_priorities(
    records: Sequence[InputRecord], priority_text: str,
) -> Tuple[Dict[str, int], List[str]]:
    """Expand ordered priority tiers to haplotype score bonuses."""
    tokens = [token.strip() for token in str(priority_text).split(",") if token.strip()]
    unique_records: Dict[str, InputRecord] = {}
    for record in records:
        unique_records.setdefault(record.haplotype, record)
    matched_groups: List[List[str]] = []
    unmatched: List[str] = []
    for token in tokens:
        if "_" in token:
            matches = [
                haplotype
                for haplotype in unique_records
                if haplotype == token
            ]
        else:
            matches = [
                haplotype
                for haplotype, record in unique_records.items()
                if record.group == token
            ]
        if matches:
            matched_groups.append(matches)
        else:
            unmatched.append(token)
    priorities: Dict[str, int] = {}
    total_included = len(matched_groups)
    for index, matches in enumerate(matched_groups):
        bonus = (total_included - index) * PRIOR_SCORE_STEP
        for haplotype in matches:
            priorities.setdefault(haplotype, bonus)
    return priorities, unmatched


def load_initial_pieces(
    fasta_path: str, records: Sequence[InputRecord],
) -> Tuple[List[Piece], Dict[str, int]]:
    pieces: List[Piece] = []
    within_haplotype: Dict[str, int] = {}
    lengths = {record.record_id: record.length for record in records}
    with IndexedFasta(fasta_path) as fasta:
        for record in records:
            rank = within_haplotype.get(record.haplotype, 0)
            within_haplotype[record.haplotype] = rank + 1
            pieces.append(
                Piece(
                    temp_id=f"input{record.order:08d}",
                    sample=record.haplotype,
                    contig=record.record_id,
                    start=0,
                    end=record.length,
                    strand="+",
                    seq=fasta.sequence(record.record_id),
                    origin_id=record.record_id,
                    lift_id=f"light_{record.order:08d}",
                    within_sample_rank=rank,
                )
            )
    return pieces, lengths


def ensure_fasta_fai(fasta_path: str) -> str:
    """Ensure the exact on-disk index required by KmerRefCleaner exists."""
    fai_path = fasta_path + ".fai"
    current = (
        os.path.isfile(fai_path)
        and os.path.getsize(fai_path) > 0
        and os.stat(fai_path).st_mtime_ns >= os.stat(fasta_path).st_mtime_ns
    )
    if not current:
        ensure_executable("samtools")
        subprocess.run(["samtools", "faidx", fasta_path], check=True)
    if not os.path.isfile(fai_path) or os.path.getsize(fai_path) == 0:
        raise RuntimeError(f"samtools faidx did not create a nonempty {fai_path}")
    return fai_path


def _read_plain_fasta(path: str) -> List[Tuple[str, str]]:
    """Read complete FASTA descriptions and sequences in file order."""
    output: List[Tuple[str, str]] = []
    description: Optional[str] = None
    sequence: List[str] = []
    with open(path, "rt") as handle:
        for line_number, raw in enumerate(handle, 1):
            if raw.startswith(">"):
                if description is not None:
                    output.append((description, "".join(sequence)))
                description = raw[1:].strip()
                sequence = []
                if not description:
                    raise ValueError(f"{path}:{line_number}: empty FASTA header")
            elif raw.strip():
                if description is None:
                    raise ValueError(f"{path}:{line_number}: sequence before header")
                sequence.append(raw.strip())
    if description is not None:
        output.append((description, "".join(sequence)))
    return output


def load_kmer_precleaned_pieces(
    input_fasta: str,
    records: Sequence[InputRecord],
    selected_ids: Set[str],
    output_fasta: str,
    executable: str,
    threads: int,
) -> Tuple[List[Piece], Dict[str, int]]:
    """Run KmerRefCleaner and lift its fragments back to original records.

    KmerRefCleaner appends ``:START-END`` to a split query ID.  This adapter
    deliberately keeps those implementation IDs out of the Python pipeline:
    every returned Piece uses the original record ID and source coordinate.
    """
    if threads < 1:
        raise ValueError("KmerRefCleaner thread count must be positive")
    if os.path.sep in executable or (os.path.altsep and os.path.altsep in executable):
        resolved = os.path.abspath(os.path.expanduser(executable))
        if not os.path.isfile(resolved) or not os.access(resolved, os.X_OK):
            raise FileNotFoundError(f"KmerRefCleaner executable is not runnable: {resolved}")
    else:
        resolved = shutil.which(executable) or ""
        if not resolved:
            raise FileNotFoundError(f"KmerRefCleaner executable not found on PATH: {executable}")

    ensure_fasta_fai(input_fasta)
    command = [
        resolved, "-i", input_fasta, "-o", output_fasta,
        "-t", str(threads),
    ]
    if selected_ids:
        command.extend(["-r", ",".join(
            record.record_id for record in records
            if record.record_id in selected_ids
        )])
    LOG.info("Running exact 31-mer pre-clean with %s", resolved)
    subprocess.run(command, check=True)
    if not os.path.isfile(output_fasta):
        raise RuntimeError(f"KmerRefCleaner did not create {output_fasta}")

    by_id = {record.record_id: record for record in records}
    original_ids = sorted(by_id, key=len, reverse=True)
    within_haplotype: Dict[str, int] = {}
    pieces: List[Piece] = []
    for output_order, (description, sequence) in enumerate(
        _read_plain_fasta(output_fasta)
    ):
        output_id = description.split()[0]
        original_id: Optional[str] = None
        start = 0
        end = 0
        if output_id in by_id:
            original_id = output_id
            end = by_id[output_id].length
        else:
            for candidate in original_ids:
                prefix = candidate + ":"
                if not output_id.startswith(prefix):
                    continue
                match = KMER_FRAGMENT_RE.fullmatch(output_id)
                if match is None or match.group(1) != candidate:
                    continue
                original_id = candidate
                start, end = int(match.group(2)), int(match.group(3))
                break
        if original_id is None:
            raise ValueError(
                f"{output_fasta}: KmerRefCleaner output ID {output_id!r} "
                "does not map to an input record"
            )
        record = by_id[original_id]
        if not (0 <= start < end <= record.length):
            raise ValueError(
                f"{output_fasta}: invalid KmerRefCleaner fragment "
                f"{output_id!r} for input length {record.length}"
            )
        if len(sequence) != end - start:
            raise ValueError(
                f"{output_fasta}: fragment {output_id!r} contains "
                f"{len(sequence)} bases but encodes {end - start}"
            )
        rank = within_haplotype.get(record.haplotype, 0)
        within_haplotype[record.haplotype] = rank + 1
        pieces.append(Piece(
            temp_id=f"kmer{output_order:08d}",
            sample=record.haplotype,
            contig=original_id,
            start=start,
            end=end,
            strand="+",
            seq=sequence,
            origin_id=original_id,
            lift_id=f"light_{record.order:08d}_kmer_{start}_{end}",
            within_sample_rank=rank,
        ))
    retained = sum(piece.length for piece in pieces)
    original = sum(record.length for record in records)
    LOG.info(
        "KmerRefCleaner retained %d fragments and %d/%d bp",
        len(pieces), retained, original,
    )
    return pieces, {record.record_id: record.length for record in records}


def _max_cycles(args: argparse.Namespace) -> Optional[int]:
    return None if args.max_cycles == 0 else args.max_cycles


def write_anchored_candidate_fasta(
    pieces: Sequence[Piece],
    source_fasta: str,
    output_fasta: str,
    anchor: int,
    old_cutoff: int,
) -> Tuple[Dict[str, AnchoredCore], Dict[str, float]]:
    """Write insertion cores with source flanks used only for self-alignment."""
    mapping: Dict[str, AnchoredCore] = {}
    thresholds: Dict[str, float] = {}
    Path(output_fasta).parent.mkdir(parents=True, exist_ok=True)
    with IndexedFasta(source_fasta) as source, open(output_fasta, "wt") as out:
        available = set(source.names())
        for piece in pieces:
            if piece.contig not in available:
                raise ValueError(
                    f"insertion source record {piece.contig!r} is absent from "
                    f"{source_fasta}"
                )
            source_length = source.length(piece.contig)
            if (
                piece.start < 0
                or piece.end > source_length
                or piece.end - piece.start != piece.length
            ):
                raise ValueError(
                    f"invalid insertion core {piece.contig}:{piece.start}-{piece.end} "
                    f"for source length {source_length}"
                )
            anchored_start = max(0, piece.start - anchor)
            anchored_end = min(source_length, piece.end + anchor)
            sequence = source.fetch(piece.contig, anchored_start, anchored_end)
            core_start = piece.start - anchored_start
            core_end = core_start + piece.length
            mapping[piece.temp_id] = AnchoredCore(
                core_start=core_start,
                core_end=core_end,
                query_length=len(sequence),
            )
            thresholds[piece.temp_id] = min(
                float(old_cutoff), 0.5 * len(sequence),
            )
            out.write(
                f">{piece.temp_id} {piece.sample}:{piece.contig}:"
                f"{anchored_start}-{anchored_end}+ "
                f"core={core_start}-{core_end}\n{wrap_fasta(sequence)}\n"
            )
    return mapping, thresholds


def _clip_intervals_to_core(
    intervals: Sequence[Tuple[int, int]], core: AnchoredCore,
) -> List[Tuple[int, int]]:
    return merge_intervals([
        (
            max(left, core.core_start) - core.core_start,
            min(right, core.core_end) - core.core_start,
        )
        for left, right in intervals
        if min(right, core.core_end) > max(left, core.core_start)
    ])


def _clip_paired_block_to_cores(
    block: Tuple[int, int, int, int],
    strand: str,
    query_core: AnchoredCore,
    target_core: AnchoredCore,
) -> Optional[Tuple[int, int, int, int]]:
    """Clip a paired alignment block so both axes lie inside their cores."""
    q0, q1, t0, t1 = block
    block_length = min(q1 - q0, t1 - t0)
    if block_length <= 0:
        return None
    offset_start = 0
    offset_end = block_length
    if strand == "-":
        offset_start = max(
            offset_start,
            q1 - query_core.core_end,
            target_core.core_start - t0,
        )
        offset_end = min(
            offset_end,
            q1 - query_core.core_start,
            target_core.core_end - t0,
        )
        if offset_end <= offset_start:
            return None
        clipped_q0 = q1 - offset_end
        clipped_q1 = q1 - offset_start
    else:
        offset_start = max(
            offset_start,
            query_core.core_start - q0,
            target_core.core_start - t0,
        )
        offset_end = min(
            offset_end,
            query_core.core_end - q0,
            target_core.core_end - t0,
        )
        if offset_end <= offset_start:
            return None
        clipped_q0 = q0 + offset_start
        clipped_q1 = q0 + offset_end
    clipped_t0 = t0 + offset_start
    clipped_t1 = t0 + offset_end
    return (
        clipped_q0 - query_core.core_start,
        clipped_q1 - query_core.core_start,
        clipped_t0 - target_core.core_start,
        clipped_t1 - target_core.core_start,
    )


def remap_anchored_self_hits(
    hits: Sequence[AlignmentHit],
    mapping: Dict[str, AnchoredCore],
) -> List[AlignmentHit]:
    """Clip anchor-space hits to paired core-space coverage."""
    remapped: List[AlignmentHit] = []
    for hit in hits:
        query_core = mapping.get(hit.query_id)
        target_core = mapping.get(hit.target_id)
        if query_core is None or target_core is None:
            continue
        pairs = [
            clipped
            for block in hit.aligned_pairs
            if (
                clipped := _clip_paired_block_to_cores(
                    block, hit.strand, query_core, target_core,
                )
            ) is not None
        ]
        if hit.aligned_pairs:
            if not pairs:
                continue
            query_intervals = merge_intervals([
                (q0, q1) for q0, q1, _t0, _t1 in pairs
            ])
            target_intervals = merge_intervals([
                (t0, t1) for _q0, _q1, t0, t1 in pairs
            ])
        else:
            # All real PAF/SAM hits carry paired blocks. Keep this fallback for
            # synthetic/API callers, but still require core coverage on both axes.
            query_intervals = _clip_intervals_to_core(
                hit.query_intervals, query_core,
            )
            target_intervals = _clip_intervals_to_core(
                hit.target_intervals, target_core,
            )
        if not query_intervals or not target_intervals:
            continue
        aligned_bases = sum(right - left for left, right in query_intervals)
        remapped.append(dataclasses.replace(
            hit,
            aligned_bases=aligned_bases,
            query_intervals=query_intervals,
            target_intervals=target_intervals,
            query_start=min(left for left, _right in query_intervals),
            query_end=max(right for _left, right in query_intervals),
            target_start=min(left for left, _right in target_intervals),
            target_end=max(right for _left, right in target_intervals),
            aligned_pairs=pairs,
        ))
    return remapped


def _collect_hits(
    query_fasta: str,
    target_fasta: str,
    cycle_dir: str,
    args: argparse.Namespace,
    broad_aligner: str,
    run_blastn: bool,
    blast_db_prefix: Optional[str] = None,
    short_self_alignment: bool = False,
    min_score_by_query: Optional[Dict[str, float]] = None,
    aggregate_query_coverage: bool = False,
    reuse_existing: bool = False,
    query_to_reference: bool = False,
):
    if query_to_reference:
        alignment_minimum = max(
            int(args.min_unmasked), QUERY_REFERENCE_MIN_SCORE,
        )
        parse_minimum = alignment_minimum
    else:
        alignment_minimum = int(args.min_unmasked)
        parse_minimum = (
            args.self_align_min if short_self_alignment else alignment_minimum
        )
    return collect_alignment_hits(
        query_fasta=query_fasta,
        target_fasta=target_fasta,
        workdir=cycle_dir,
        threads=args.cores,
        min_identity=args.min_identity,
        min_segment=alignment_minimum,
        log=LOG,
        winnow_params={
            **LIGHT_WINNOW_PARAMS,
            "m": parse_minimum,
        },
        blast_word_size=(
            args.self_blast_word_size
            if short_self_alignment else args.blast_word_size
        ),
        blast_evalue=args.blast_evalue,
        blast_max_target_seqs=args.blast_max_target_seqs,
        prefix="align",
        skip_blastn=args.skip_blastn or not run_blastn,
        blast_mode="independent",
        broad_aligner=broad_aligner,
        minimap_params=MINIMAP_PARAMS,
        blast_db_prefix=blast_db_prefix,
        masked_weight=MASKED_BASE_WEIGHT,
        min_score_by_query=min_score_by_query,
        parse_min_segment=parse_minimum,
        retain_aligned_pairs=short_self_alignment,
        aggregate_query_coverage=aggregate_query_coverage,
        query_workers=args.cores,
        reuse_existing=reuse_existing,
    )


def _self_clean_phase(
    pieces: Sequence[Piece],
    phase_dir: str,
    args: argparse.Namespace,
    broad_aligner: str,
    run_blastn: bool,
    priorities: Dict[str, int],
    anchored_source_fasta: Optional[str] = None,
) -> List[Piece]:
    Path(phase_dir).mkdir(parents=True, exist_ok=True)

    def align_fn(current: Sequence[Piece], cycle: int):
        cycle_dir = os.path.join(phase_dir, f"cycle{cycle:02d}")
        Path(cycle_dir).mkdir(parents=True, exist_ok=True)
        candidate_fasta = os.path.join(cycle_dir, "candidates.fa")
        candidate_tmp = candidate_fasta + f".tmp.{os.getpid()}"
        anchor_mapping = None
        query_thresholds = None
        if anchored_source_fasta is not None:
            anchor_mapping, query_thresholds = write_anchored_candidate_fasta(
                current,
                anchored_source_fasta,
                candidate_tmp,
                args.self_anchor,
                args.min_unmasked,
            )
        else:
            write_fasta_pieces(current, candidate_tmp)
        reuse_alignment = (
            bool(getattr(args, "resume", False))
            and os.path.isfile(candidate_fasta)
            and filecmp.cmp(candidate_tmp, candidate_fasta, shallow=False)
        )
        if reuse_alignment:
            os.unlink(candidate_tmp)
            LOG.info("REUSE: self-clean cycle candidate FASTA %s", candidate_fasta)
        else:
            os.replace(candidate_tmp, candidate_fasta)
        hits = _collect_hits(
            candidate_fasta,
            candidate_fasta,
            cycle_dir,
            args,
            broad_aligner=broad_aligner,
            run_blastn=run_blastn,
            short_self_alignment=anchored_source_fasta is not None,
            min_score_by_query=query_thresholds,
            reuse_existing=reuse_alignment,
        )
        if anchor_mapping is not None:
            hits = remap_anchored_self_hits(hits, anchor_mapping)
        return hits

    align_fn.checkpoint_tag = repr((
        "light-self-clean-v2",
        broad_aligner,
        bool(run_blastn),
        float(args.min_identity),
        int(args.min_unmasked),
        int(args.self_align_min),
        int(args.self_anchor),
        int(args.self_blast_word_size),
        int(args.blast_word_size),
        str(args.blast_evalue),
        int(args.blast_max_target_seqs),
        bool(args.skip_blastn),
        tuple(sorted(priorities.items())),
    ))

    return self_clean_converge(
        pieces,
        align_fn,
        args.min_unmasked,
        "insertion",
        LOG,
        max_cycles=_max_cycles(args),
        sample_priorities=priorities,
        convergence_bp=args.min_unmasked,
        masked_weight=MASKED_BASE_WEIGHT,
        checkpoint_dir=os.path.join(phase_dir, "checkpoints"),
        resume=bool(getattr(args, "resume", False)),
    )


def self_clean_two_stage(
    pieces: Sequence[Piece],
    workdir: str,
    args: argparse.Namespace,
    priorities: Dict[str, int],
    label: str,
    anchored_source_fasta: Optional[str] = None,
) -> List[Piece]:
    if not pieces:
        return []
    before = len(pieces)
    coarse = _self_clean_phase(
        pieces,
        os.path.join(workdir, "minimap2"),
        args,
        broad_aligner="minimap2",
        run_blastn=False,
        priorities=priorities,
        anchored_source_fasta=anchored_source_fasta,
    )
    LOG.info("%s minimap2 self-clean retained %d/%d pieces", label, len(coarse), before)
    if not coarse:
        return []
    cleaned = _self_clean_phase(
        coarse,
        os.path.join(workdir, "winnowmap_blastn"),
        args,
        broad_aligner="winnowmap",
        run_blastn=True,
        priorities=priorities,
        anchored_source_fasta=anchored_source_fasta,
    )
    LOG.info(
        "%s Winnowmap+BLASTN self-clean retained %d/%d minimap2 residual pieces",
        label,
        len(cleaned),
        len(coarse),
    )
    return cleaned


def residual_pieces(
    pieces: Sequence[Piece], coverage: Dict[str, List[Tuple[int, int]]],
    min_unmasked: int, masked_weight: float = MASKED_BASE_WEIGHT,
) -> List[Piece]:
    """Subtract represented query coverage and retain valid N-free residuals."""
    residuals: List[Piece] = []
    for piece in pieces:
        intervals = find_unaligned_segments(
            coverage.get(piece.temp_id, ()),
            piece.length,
            piece.seq,
            min_unmasked,
            masked_weight,
        )
        child_index = 0
        for left, right in intervals:
            for nleft, nright, sequence in iter_n_free_segments(piece.seq[left:right]):
                if novelty_score(sequence, masked_weight) < min_unmasked:
                    continue
                child_index += 1
                source_left = left + nleft
                source_right = left + nright
                residuals.append(
                    Piece(
                        temp_id=f"{piece.temp_id}_res{child_index:04d}",
                        sample=piece.sample,
                        contig=piece.contig,
                        start=piece.start + source_left,
                        end=piece.start + source_right,
                        strand=piece.strand,
                        seq=sequence,
                        parent_id=piece.temp_id,
                        note="reference_residual",
                        origin_id=piece.origin_id,
                        lift_id=piece.lift_id,
                        within_sample_rank=piece.within_sample_rank,
                    )
                )
    return residuals


def _fixed_phase_policy_signature(
    current: Sequence[Piece],
    reference_fasta: str,
    args: argparse.Namespace,
    broad_aligner: str,
    run_blastn: bool,
) -> str:
    """Fingerprint candidate state, reference file, and exclusion policy."""
    reference_stat = os.stat(reference_fasta)
    policy = {
        "format": "fixed-reference-exclusion-v2",
        "broad_aligner": broad_aligner,
        "run_blastn": bool(run_blastn),
        "reference_path": os.path.realpath(reference_fasta),
        "reference_size": reference_stat.st_size,
        "reference_mtime_ns": reference_stat.st_mtime_ns,
        "min_identity": float(args.min_identity),
        "min_unmasked": int(args.min_unmasked),
        "query_reference_min_score": QUERY_REFERENCE_MIN_SCORE,
        "blast_word_size": int(args.blast_word_size),
        "blast_evalue": str(args.blast_evalue),
        "blast_max_target_seqs": int(args.blast_max_target_seqs),
        "skip_blastn": bool(args.skip_blastn),
        "pieces": _pieces_signature(current),
    }
    return hashlib.sha256(json.dumps(
        policy, sort_keys=True, separators=(",", ":"),
    ).encode()).hexdigest()


def _alignment_target_signature_from_pieces(pieces: Sequence[Piece]) -> str:
    """Hash the ordered IDs and bases visible to an aligner."""
    digest = hashlib.sha256()
    for piece in pieces:
        digest.update(piece.temp_id.encode("utf-8"))
        digest.update(b"\0")
        digest.update(piece.seq.encode("ascii"))
        digest.update(b"\n")
    return digest.hexdigest()


def _alignment_target_signature_from_fasta(path: str) -> str:
    ensure_large_input_fai(path)
    digest = hashlib.sha256()
    with IndexedFasta(path) as fasta:
        for name in fasta.names():
            digest.update(name.encode("utf-8"))
            digest.update(b"\0")
            digest.update(fasta.sequence(name).encode("ascii"))
            digest.update(b"\n")
    return digest.hexdigest()


def _load_fixed_candidate_fasta(
    path: str,
    original_pieces: Sequence[Piece],
) -> List[Piece]:
    """Restore cycle candidates, borrowing stable provenance from originals."""
    descriptions: List[Tuple[str, str]] = []
    with open(path, "rt") as handle:
        for line_number, raw in enumerate(handle, 1):
            if not raw.startswith(">"):
                continue
            fields = raw[1:].strip().split()
            if len(fields) < 2:
                raise ValueError(
                    f"{path}:{line_number}: candidate header lacks source coordinates"
                )
            descriptions.append((fields[0], fields[1]))
    if not descriptions:
        raise ValueError(f"{path}: no candidate FASTA records")
    provenance: Dict[str, List[Piece]] = defaultdict(list)
    for piece in original_pieces:
        provenance[piece.contig].append(piece)
    restored: List[Piece] = []
    ensure_large_input_fai(path)
    with IndexedFasta(path) as sequences:
        available = set(sequences.names())
        for temp_id, source in descriptions:
            match = PIECE_SOURCE_RE.fullmatch(source)
            if match is None:
                raise ValueError(
                    f"{path}: invalid candidate source coordinate {source!r}"
                )
            sample, contig, start_text, end_text, strand = match.groups()
            start, end = int(start_text), int(end_text)
            possible = [
                piece
                for piece in provenance.get(contig, [])
                if (
                    piece.sample == sample
                    and piece.strand == strand
                    and piece.start <= start
                    and end <= piece.end
                )
            ]
            if not possible:
                raise ValueError(
                    f"{path}: candidate source {sample}:{contig}:"
                    f"{start}-{end}{strand} is not contained in any current "
                    "input piece"
                )
            # Residual pieces from one original record can share a contig.
            # The narrowest containing parent is the immediate provenance.
            original = min(possible, key=lambda piece: (piece.length, piece.start))
            if temp_id not in available:
                raise ValueError(
                    f"{path}: candidate record {temp_id!r} is absent from its index"
                )
            stored_length = sequences.length(temp_id)
            source_left = start - original.start
            source_right = end - original.start
            sequence = original.seq[source_left:source_right]
            if end - start != stored_length or stored_length != len(sequence):
                raise ValueError(
                    f"{path}: {temp_id!r} coordinates span {end - start} bp "
                    f"but its stored/current-source lengths are "
                    f"{stored_length}/{len(sequence)} bp"
                )
            restored.append(Piece(
                temp_id=temp_id,
                sample=sample,
                contig=contig,
                start=start,
                end=end,
                strand=strand,
                seq=sequence,
                parent_id=original.parent_id,
                note=original.note,
                origin_id=original.origin_id,
                lift_id=original.lift_id,
                within_sample_rank=original.within_sample_rank,
            ))
    return restored


def _load_reusable_next_phase_input(
    candidate_fasta: str,
    original_pieces: Sequence[Piece],
    source_fasta: str,
) -> Optional[List[Piece]]:
    """Load the serialized output of a completed preceding phase.

    Cycle 0 of the next phase is materialized from the preceding phase's
    returned pieces before the next aligner starts.  It is therefore a later,
    stronger resume point than any PAF inside the preceding phase.
    """
    if (
        not os.path.isfile(candidate_fasta)
        or os.path.getsize(candidate_fasta) == 0
        or os.stat(candidate_fasta).st_mtime_ns
            < os.stat(source_fasta).st_mtime_ns
    ):
        return None
    try:
        return _load_fixed_candidate_fasta(candidate_fasta, original_pieces)
    except (OSError, ValueError) as error:
        LOG.warning(
            "Cannot use next-phase resume FASTA %s: %s",
            candidate_fasta,
            error,
        )
        return None


def _latest_reusable_fixed_cycle(
    phase_dir: str,
    broad_aligner: str,
    original_pieces: Sequence[Piece],
    source_fasta: str,
) -> Optional[Tuple[int, List[Piece], str]]:
    """Find a current legacy/incomplete cycle without regenerating its FASTA."""
    # The cleaned reference is a derived file and may be regenerated with
    # identical contents during resume. Its newer mtime must not invalidate
    # every existing cycle. The original input is the stable source dependency;
    # within each cycle the PAF must be newer than its own candidates.fa.
    source_mtime_ns = os.stat(source_fasta).st_mtime_ns
    suffix = f"align.{broad_aligner}.paf"
    candidates: Dict[int, Tuple[str, str]] = {}
    if not os.path.isdir(phase_dir):
        return None
    for entry in os.scandir(phase_dir):
        match = re.fullmatch(r"cycle(\d+)", entry.name)
        if match is None or not entry.is_dir(follow_symlinks=False):
            continue
        candidate_fasta = os.path.join(entry.path, "candidates.fa")
        paf = os.path.join(entry.path, suffix)
        if not os.path.isfile(candidate_fasta) or not os.path.isfile(paf):
            continue
        if (
            os.path.getsize(candidate_fasta) == 0
            or os.path.getsize(paf) == 0
            or os.stat(candidate_fasta).st_mtime_ns < source_mtime_ns
            or os.stat(paf).st_mtime_ns < os.stat(candidate_fasta).st_mtime_ns
        ):
            continue
        candidates[int(match.group(1))] = (candidate_fasta, paf)
    if not candidates:
        return None
    # A high-numbered directory alone is not sufficient evidence: it could be
    # stale debris from a different run. Require the legacy cycle chain to be
    # contiguous, and require each next candidate state to have been written
    # after the preceding cycle's complete PAF.
    latest = -1
    previous_paf: Optional[str] = None
    for cycle in range(max(candidates) + 1):
        if cycle not in candidates:
            break
        candidate_fasta, paf = candidates[cycle]
        if (
            previous_paf is not None
            and os.stat(candidate_fasta).st_mtime_ns
                < os.stat(previous_paf).st_mtime_ns
        ):
            break
        latest = cycle
        previous_paf = paf
    if latest < 0:
        return None
    cycle = latest
    candidate_fasta, _paf = candidates[cycle]
    restored = _load_fixed_candidate_fasta(candidate_fasta, original_pieces)
    return cycle, restored, candidate_fasta


def exclude_against_reference_phase(
    pieces: Sequence[Piece],
    reference_fasta: str,
    phase_dir: str,
    args: argparse.Namespace,
    broad_aligner: str,
    run_blastn: bool,
    source_fasta: Optional[str] = None,
) -> List[Piece]:
    """Clean all query pieces together against one fixed reference FASTA."""
    current = list(pieces)
    if not current:
        return []
    Path(phase_dir).mkdir(parents=True, exist_ok=True)
    target = reference_fasta
    blast_db_prefix = None
    if broad_aligner == "minimap2":
        target = build_minimap2_index(
            reference_fasta,
            os.path.join(phase_dir, "reference.asm5.mmi"),
            LOG,
            preset=str(MINIMAP_PARAMS["preset"]),
        )
    elif run_blastn and not args.skip_blastn:
        blast_db_prefix = os.path.join(phase_dir, "reference_blastdb", "ref")

    converged = False
    cycle = 0
    limit = _max_cycles(args)
    checkpoint_root = os.path.join(phase_dir, "checkpoints")
    resume = bool(getattr(args, "resume", False))
    stable_source_fasta = source_fasta or getattr(args, "input", None)
    restored_candidate_fasta: Optional[str] = None

    # First replay completed compact checkpoints. This reads only saved state;
    # it does not regenerate any earlier candidate FASTA.
    if resume:
        if os.path.isdir(checkpoint_root):
            while limit is None or cycle < limit:
                # Candidate IDs are deterministically re-based before every
                # alignment cycle. Checkpoint validation must fingerprint that
                # same state, not the caller's pre-cycle temporary IDs.
                for index, piece in enumerate(current):
                    piece.temp_id = f"query{index:08d}"
                signature = _fixed_phase_policy_signature(
                    current, reference_fasta, args, broad_aligner, run_blastn,
                )
                checkpoint = _read_cycle_checkpoint(
                    os.path.join(checkpoint_root, f"cycle{cycle:02d}"), signature,
                )
                if checkpoint is None:
                    break
                current, removed_bp, cycle_converged = checkpoint
                LOG.info(
                    "RESUME fixed-reference [%s] cycle %d: loaded %d residual "
                    "pieces (removed %d bp)",
                    broad_aligner, cycle, len(current), removed_bp,
                )
                if cycle_converged:
                    return current
                cycle += 1

        # Backward-compatible recovery for a run created before fixed-target
        # state checkpoints existed. The latest current candidates.fa is the
        # input state of that unfinished cycle and its PAF can be reused.
        if source_fasta is not None:
            legacy = _latest_reusable_fixed_cycle(
                phase_dir, broad_aligner, pieces, source_fasta,
            )
            if legacy is not None and legacy[0] >= cycle:
                cycle, current, candidate_fasta = legacy
                restored_candidate_fasta = os.path.realpath(candidate_fasta)
                LOG.info(
                    "RESUME fixed-reference [%s] directly at cycle %d from "
                    "existing candidate FASTA %s",
                    broad_aligner, cycle, candidate_fasta,
                )

    while limit is None or cycle < limit:
        if not current:
            converged = True
            break
        for index, piece in enumerate(current):
            piece.temp_id = f"query{index:08d}"
        cycle_dir = os.path.join(phase_dir, f"cycle{cycle:02d}")
        Path(cycle_dir).mkdir(parents=True, exist_ok=True)
        candidate_fasta = os.path.join(cycle_dir, "candidates.fa")
        expected_paf = os.path.join(
            cycle_dir, f"align.{broad_aligner}.paf",
        )
        input_signature = _fixed_phase_policy_signature(
            current, reference_fasta, args, broad_aligner, run_blastn,
        )
        candidate_signature_path = os.path.join(
            cycle_dir, "candidates.input_signature.txt",
        )
        saved_candidate_signature = None
        try:
            with open(candidate_signature_path, "rt") as signature_handle:
                saved_candidate_signature = signature_handle.read().strip()
        except OSError:
            pass
        candidate_state_matches = (
            saved_candidate_signature == input_signature
            or (
                restored_candidate_fasta is not None
                and os.path.realpath(candidate_fasta)
                    == restored_candidate_fasta
            )
        )
        reuse_alignment = (
            resume
            and candidate_state_matches
            and os.path.isfile(candidate_fasta)
            and os.path.getsize(candidate_fasta) > 0
            and os.path.isfile(expected_paf)
            and os.path.getsize(expected_paf) > 0
            and (
                stable_source_fasta is None
                or os.stat(candidate_fasta).st_mtime_ns
                    >= os.stat(stable_source_fasta).st_mtime_ns
            )
            and os.stat(expected_paf).st_mtime_ns
                >= os.stat(candidate_fasta).st_mtime_ns
        )
        if reuse_alignment:
            LOG.info(
                "REUSE: fixed-target %s cycle candidate FASTA and PAF %s",
                broad_aligner,
                candidate_fasta,
            )
        else:
            candidate_tmp = candidate_fasta + f".tmp.{os.getpid()}"
            signature_tmp = candidate_signature_path + f".tmp.{os.getpid()}"
            write_fasta_pieces(current, candidate_tmp)
            try:
                with open(signature_tmp, "wt") as signature_handle:
                    signature_handle.write(input_signature + "\n")
                    signature_handle.flush()
                    os.fsync(signature_handle.fileno())
                os.replace(candidate_tmp, candidate_fasta)
                os.replace(signature_tmp, candidate_signature_path)
            finally:
                for temporary in (candidate_tmp, signature_tmp):
                    try:
                        os.remove(temporary)
                    except FileNotFoundError:
                        pass
        # Legacy work directories often contain a very large candidates.fa
        # but no .fai. Without this check IndexedFasta silently scans the whole
        # file before the PAF parser emits its next progress message.
        ensure_large_input_fai(candidate_fasta)
        hits = _collect_hits(
            candidate_fasta,
            target,
            cycle_dir,
            args,
            broad_aligner=broad_aligner,
            run_blastn=run_blastn,
            blast_db_prefix=blast_db_prefix,
            aggregate_query_coverage=True,
            reuse_existing=reuse_alignment,
            query_to_reference=True,
        )
        coverage = coverage_by_query_id(
            hits, {piece.temp_id: piece.length for piece in current},
        )
        kept = residual_pieces(current, coverage, args.min_unmasked)
        removed_bp = max(
            0,
            sum(piece.length for piece in current)
            - sum(piece.length for piece in kept),
        )
        same_piece_count = len(kept) == len(current)
        LOG.info(
            "Reference exclusion [%s] cycle %d: %d -> %d pieces, removed %d bp",
            broad_aligner,
            cycle,
            len(current),
            len(kept),
            removed_bp,
        )
        current = kept
        cycle_converged = removed_bp == 0 or (
            same_piece_count and removed_bp < args.min_unmasked
        )
        _write_cycle_checkpoint(
            os.path.join(checkpoint_root, f"cycle{cycle:02d}"),
            current,
            input_signature,
            removed_bp,
            cycle_converged,
        )
        if cycle_converged:
            converged = True
            break
        cycle += 1
    if current and not converged and limit is not None:
        LOG.warning(
            "Reference exclusion [%s] reached max_cycles=%d before convergence",
            broad_aligner,
            limit,
        )
    return current


def exclude_against_reference_two_stage(
    pieces: Sequence[Piece],
    reference_fasta: str,
    workdir: str,
    args: argparse.Namespace,
) -> List[Piece]:
    winnowmap_dir = os.path.join(workdir, "winnowmap_blastn")
    coarse: Optional[List[Piece]] = None
    if bool(getattr(args, "resume", False)):
        next_phase_input = os.path.join(
            winnowmap_dir, "cycle00", "candidates.fa",
        )
        coarse = _load_reusable_next_phase_input(
            next_phase_input, pieces, args.input,
        )
        if coarse is not None:
            LOG.info(
                "RESUME: using %s as the completed minimap2-phase output "
                "(%d pieces); minimap2 cycles will not be reparsed",
                next_phase_input,
                len(coarse),
            )
    if coarse is None:
        coarse = exclude_against_reference_phase(
            pieces,
            reference_fasta,
            os.path.join(workdir, "minimap2"),
            args,
            broad_aligner="minimap2",
            run_blastn=False,
            source_fasta=args.input,
        )
    LOG.info(
        "Query-vs-reference minimap2 retained %d/%d pieces", len(coarse), len(pieces),
    )
    if not coarse:
        return []
    cleaned = exclude_against_reference_phase(
        coarse,
        reference_fasta,
        winnowmap_dir,
        args,
        broad_aligner="winnowmap",
        run_blastn=True,
        source_fasta=args.input,
    )
    LOG.info(
        "Query-vs-reference Winnowmap+BLASTN retained %d/%d minimap2 residual pieces",
        len(cleaned),
        len(coarse),
    )
    return cleaned


def write_bed(
    output_path: str,
    reference_pieces: Sequence[Piece],
    query_pieces: Sequence[Piece],
    original_lengths: Dict[str, int],
    record_order: Dict[str, int],
    anchor: int,
    record_coordinates: Dict[str, Tuple[str, int, int]],
) -> int:
    """Write BED6+haplotype+class+core-novelty-score atomically."""
    rows: List[Tuple[int, int, int, int, Piece, str]] = []
    serial = 0
    for class_rank, (kind, pieces) in enumerate(
        (("reference", reference_pieces), ("unique_query", query_pieces))
    ):
        for piece in pieces:
            serial += 1
            rows.append(
                (
                    record_order[piece.contig],
                    piece.start,
                    piece.end,
                    class_rank,
                    piece,
                    kind,
                )
            )
    rows.sort(key=lambda row: row[:4])
    Path(output_path).parent.mkdir(parents=True, exist_ok=True)
    temporary = output_path + ".tmp"
    try:
        with open(temporary, "wt") as out:
            for index, (_order, _start, _end, _class_rank, piece, kind) in enumerate(
                rows, 1
            ):
                start = max(0, piece.start - anchor)
                end = min(original_lengths[piece.contig], piece.end + anchor)
                try:
                    source_contig, source_start, source_end = (
                        record_coordinates[piece.contig]
                    )
                except KeyError as error:
                    raise ValueError(
                        f"input record {piece.contig!r} lacks absolute assembly "
                        "coordinates in its FASTA header"
                    ) from error
                if source_end - source_start != original_lengths[piece.contig]:
                    raise ValueError(
                        f"input record {piece.contig!r} has FASTA length "
                        f"{original_lengths[piece.contig]} but header interval "
                        f"{source_contig}:{source_start}-{source_end} has length "
                        f"{source_end - source_start}"
                    )
                absolute_start = source_start + start
                absolute_end = source_start + end
                out.write(
                    "\t".join(
                        [
                            piece.contig,
                            str(start),
                            str(end),
                            f"unique_region_{index:06d}_{source_contig}_"
                            f"{absolute_start}_{absolute_end}",
                            "0",
                            piece.strand,
                            piece.sample,
                            kind,
                            f"{novelty_score(piece.seq, MASKED_BASE_WEIGHT):.1f}",
                        ]
                    )
                    + "\n"
                )
        os.replace(temporary, output_path)
    finally:
        try:
            os.remove(temporary)
        except FileNotFoundError:
            pass
    return len(rows)


def write_fasta(
    output_path: str,
    reference_pieces: Sequence[Piece],
    query_pieces: Sequence[Piece],
    original_lengths: Dict[str, int],
    record_order: Dict[str, int],
    anchor: int,
    record_coordinates: Dict[str, Tuple[str, int, int]],
    source_fasta: str,
) -> int:
    """Write retained anchored intervals as an ordinary multi-FASTA.

    IDs end in source ``START_END`` coordinates so the result can be consumed
    directly by graphmake-style tools and by the local graph linearizer.
    """
    rows: List[Tuple[int, int, int, int, Piece, str]] = []
    for class_rank, (kind, pieces) in enumerate(
        (("reference", reference_pieces), ("unique_query", query_pieces))
    ):
        for piece in pieces:
            rows.append((
                record_order[piece.contig], piece.start, piece.end,
                class_rank, piece, kind,
            ))
    rows.sort(key=lambda row: row[:4])
    Path(output_path).parent.mkdir(parents=True, exist_ok=True)
    temporary = output_path + ".tmp"
    try:
        with IndexedFasta(source_fasta) as source, open(temporary, "wt") as out:
            for index, (_order, _start, _end, _rank, piece, kind) in enumerate(
                rows, 1
            ):
                start = max(0, piece.start - anchor)
                end = min(original_lengths[piece.contig], piece.end + anchor)
                sequence = source.fetch(piece.contig, start, end)
                source_contig, source_start, source_end = record_coordinates.get(
                    piece.contig, (piece.contig, 0, original_lengths[piece.contig]),
                )
                if source_end - source_start != original_lengths[piece.contig]:
                    raise ValueError(
                        f"input record {piece.contig!r} has FASTA length "
                        f"{original_lengths[piece.contig]} but header interval "
                        f"{source_contig}:{source_start}-{source_end} has length "
                        f"{source_end - source_start}"
                    )
                absolute_start = source_start + start
                absolute_end = source_start + end
                safe_contig = re.sub(r"[^A-Za-z0-9.#:-]+", "_", source_contig)
                record_id = (
                    f"minset{index:08d}_{piece.sample}_{safe_contig}_"
                    f"{absolute_start}_{absolute_end}"
                )
                out.write(
                    f">{record_id} {piece.sample}:{source_contig}:{absolute_start}-"
                    f"{absolute_end}{piece.strand} source={piece.contig}:"
                    f"{start}-{end}{piece.strand} class={kind} "
                    f"novelty_score={novelty_score(piece.seq, MASKED_BASE_WEIGHT):.1f}\n"
                    f"{wrap_fasta(sequence)}\n"
                )
        os.replace(temporary, output_path)
    finally:
        try:
            os.remove(temporary)
        except FileNotFoundError:
            pass
    return len(rows)


def run(args: argparse.Namespace) -> int:
    for executable in ("minimap2", "winnowmap"):
        ensure_executable(executable)
    if not args.skip_blastn:
        for executable in ("makeblastdb", "blastn"):
            ensure_executable(executable)

    ensure_large_input_fai(args.input)
    records = read_input_records(args.input)
    priorities, unmatched_priorities = build_haplotype_priorities(
        records, args.priority,
    )
    for token in unmatched_priorities:
        LOG.warning("Priority selector matched no input haplotype: %s", token)
    if priorities:
        LOG.info("Applied %d ordered haplotype priority assignment(s)", len(priorities))
    else:
        LOG.info("No priority tier applied; longer current pieces win")

    selectors = expand_reference_selectors(args.reference)
    selected_ids, unmatched_selectors, used_fallback = select_reference_ids(
        records, selectors, priorities,
    )
    for selector in unmatched_selectors:
        LOG.warning("Reference selector matched no input record: %s", selector)
    if used_fallback:
        selected_records = [
            record for record in records if record.record_id in selected_ids
        ]
        selected_haplotype = selected_records[0].haplotype
        selection_rule = (
            "highest-priority fallback haplotype"
            if priorities
            else "largest-sequence fallback haplotype"
        )
        if selectors:
            LOG.warning(
                "No -r selector matched; using %s %s (%d records, %d bp)",
                selection_rule,
                selected_haplotype,
                len(selected_records),
                sum(record.length for record in selected_records),
            )
        else:
            LOG.info(
                "No -r selector supplied; using %s %s (%d records, %d bp)",
                selection_rule,
                selected_haplotype,
                len(selected_records),
                sum(record.length for record in selected_records),
            )
    LOG.info(
        "Loaded %d records; selected %d complete record(s) as the reference set",
        len(records),
        len(selected_ids),
    )
    for record in records:
        if record.record_id in selected_ids:
            LOG.info(
                "Selected reference record: %s (haplotype=%s, length=%d)",
                record.record_id,
                record.haplotype,
                record.length,
            )

    if getattr(args, "kmer_preclean", False):
        initial, original_lengths = load_kmer_precleaned_pieces(
            args.input,
            records,
            selected_ids,
            os.path.join(args.active_workdir, "kmer_precleaned.fa"),
            args.kmer_ref_cleaner,
            args.cores,
        )
    else:
        initial, original_lengths = load_initial_pieces(args.input, records)
    reference_input = [piece for piece in initial if piece.contig in selected_ids]
    query_input = [piece for piece in initial if piece.contig not in selected_ids]

    # First remove redundancy within the selected reference haplotype. Only
    # this cleaned reference set may act as the fixed target for the remaining
    # samples. The retained reference pieces join the query residuals once more
    # in the final all-vs-all self-clean, so the final BED contains only pieces
    # that survived the same cohort-wide redundancy decision.
    fixed_reference = self_clean_two_stage(
        reference_input,
        os.path.join(args.active_workdir, "reference_selfclean"),
        args,
        priorities,
        "Selected-reference",
    )
    if not fixed_reference:
        raise RuntimeError(
            "selected-reference self-clean retained no valid reference piece"
        )
    reference_fasta = os.path.join(args.active_workdir, "cleaned_reference.fa")
    reuse_reference_fasta = (
        bool(getattr(args, "resume", False))
        and os.path.isfile(reference_fasta)
        and os.path.getsize(reference_fasta) > 0
    )
    if reuse_reference_fasta:
        reuse_reference_fasta = (
            _alignment_target_signature_from_fasta(reference_fasta)
            == _alignment_target_signature_from_pieces(fixed_reference)
        )
    if reuse_reference_fasta:
        LOG.info(
            "RESUME: exact content match; retaining existing cleaned "
            "reference FASTA without regeneration: %s",
            reference_fasta,
        )
    else:
        write_fasta_pieces(fixed_reference, reference_fasta)
        # Old fixed-target PAFs were produced against a different reference.
        args.resume = False

    query_residuals = exclude_against_reference_two_stage(
        query_input,
        reference_fasta,
        os.path.join(args.active_workdir, "query_vs_reference"),
        args,
    )
    jointly_cleaned = self_clean_two_stage(
        [*fixed_reference, *query_residuals],
        os.path.join(args.active_workdir, "reference_query_selfclean"),
        args,
        priorities,
        "Reference-plus-query residual",
        anchored_source_fasta=args.input,
    )
    cleaned_reference = [
        piece for piece in jointly_cleaned if piece.contig in selected_ids
    ]
    query_cleaned = [
        piece for piece in jointly_cleaned if piece.contig not in selected_ids
    ]
    if not cleaned_reference:
        raise RuntimeError(
            "joint self-clean removed every selected reference piece; make "
            "the intended reference haplotype earlier in --priority"
        )

    record_coordinates = {
        record.record_id: (
            record.genomic_contig,
            record.genomic_start,
            record.genomic_end,
        )
        for record in records
        if (
            record.genomic_contig is not None
            and record.genomic_start is not None
            and record.genomic_end is not None
        )
    }
    output_format = args.output_format
    if output_format == "auto":
        output_format = (
            "bed" if str(args.output).lower().endswith((".bed", ".bed.gz"))
            else "fasta"
        )
    if output_format == "bed":
        count = write_bed(
            args.output,
            cleaned_reference,
            query_cleaned,
            original_lengths,
            {record.record_id: record.order for record in records},
            args.anchor,
            record_coordinates,
        )
    else:
        count = write_fasta(
            args.output,
            cleaned_reference,
            query_cleaned,
            original_lengths,
            {record.record_id: record.order for record in records},
            args.anchor,
            record_coordinates,
            args.input,
        )
    LOG.info(
        "Wrote %d retained local %s records (%d reference, %d unique query) to %s",
        count,
        output_format,
        len(cleaned_reference),
        len(query_cleaned),
        args.output,
    )
    return count


def parse_args(argv: Optional[Sequence[str]] = None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=(
            "Remove redundancy from one homologous-sequence FASTA and emit "
            "retained intervals as BED or retained sequences as FASTA."
        )
    )
    parser.add_argument("-i", "--input", required=True, help="input multi-FASTA")
    parser.add_argument("-o", "--output", required=True, help="output BED or FASTA path")
    parser.add_argument(
        "--output-format", choices=("auto", "bed", "fasta"), default="auto",
        help="output format; auto uses BED for .bed/.bed.gz and FASTA otherwise",
    )
    parser.add_argument(
        "-r",
        "--reference",
        action="append",
        default=[],
        help=(
            "repeatable reference selector or selector-file path; selectors are "
            "contig:start-end, an exact fields-2+3 haplotype, or an exact field-2 group"
        ),
    )
    parser.add_argument(
        "-p",
        "--priority",
        default="",
        help=(
            "ordered comma-separated haplotype/group priorities; each earlier "
            "matched tier adds 10G to the current-piece length score"
        ),
    )
    parser.add_argument(
        "-c", "--cores", type=int, default=8, help="threads for every aligner (default: 8)"
    )
    parser.add_argument(
        "--kmer-preclean",
        action="store_true",
        help="run KmerRefCleaner as a direct executable before alignment cleaning",
    )
    parser.add_argument(
        "--kmer-ref-cleaner",
        default=str(Path(__file__).resolve().with_name("KmerRefCleaner")),
        help="KmerRefCleaner executable path or command name",
    )
    parser.add_argument(
        "--min-score",
        "--min-unmasked",
        dest="min_unmasked",
        type=int,
        default=100,
        help=(
            "minimum uppercase + 0.3*lowercase A/C/G/T score in a valid "
            "retained/aligned piece (default: 100)"
        ),
    )
    parser.add_argument(
        "--anchor",
        type=int,
        default=50,
        help="expand final retained BED intervals on each side (default: 50)",
    )
    parser.add_argument(
        "--min-identity",
        type=float,
        default=95.0,
        help="minimum alignment identity percentage (default: 95)",
    )
    parser.add_argument(
        "--max-cycles",
        type=int,
        default=0,
        help="cycle cap per phase; 0 means run to convergence (default: 0)",
    )
    parser.add_argument("--blast-word-size", type=int, default=50)
    parser.add_argument(
        "--self-align-min",
        type=int,
        default=50,
        help=(
            "minimum raw alignment length parsed during anchored insertion "
            "self-alignment (default: 50)"
        ),
    )
    parser.add_argument(
        "--self-blast-word-size",
        type=int,
        default=28,
        help="BLASTN word size for anchored insertion self-alignment (default: 28)",
    )
    parser.add_argument(
        "--self-anchor",
        type=int,
        default=50,
        help="source context added to each insertion self-alignment side (default: 50)",
    )
    parser.add_argument("--blast-evalue", default="1e-300")
    parser.add_argument("--blast-max-target-seqs", type=int, default=100)
    parser.add_argument(
        "--skip-blastn",
        action="store_true",
        help="run only minimap2 and Winnowmap (not recommended for maximum sensitivity)",
    )
    parser.add_argument(
        "--workdir",
        help="explicit persistent work directory (always retained)",
    )
    parser.add_argument(
        "--keep-work",
        action="store_true",
        help="retain an OUTPUT.work directory when --workdir is not supplied",
    )
    parser.add_argument(
        "--resume",
        action="store_true",
        help=(
            "reuse validated cycle checkpoints and completed alignments from "
            "--workdir or the stable --keep-work OUTPUT.work directory"
        ),
    )
    parser.add_argument(
        "--log-level",
        choices=("DEBUG", "INFO", "WARNING", "ERROR"),
        default="INFO",
    )
    args = parser.parse_args(argv)
    if args.cores < 1:
        parser.error("--cores must be at least 1")
    if args.min_unmasked < 1:
        parser.error("--min-score must be at least 1")
    if args.anchor < 0:
        parser.error("--anchor cannot be negative")
    if args.max_cycles < 0:
        parser.error("--max-cycles cannot be negative")
    if args.self_align_min < 1:
        parser.error("--self-align-min must be at least 1")
    if args.self_blast_word_size < 4:
        parser.error("--self-blast-word-size must be at least 4")
    if args.self_anchor < 0:
        parser.error("--self-anchor cannot be negative")
    if not (0.0 < args.min_identity <= 100.0):
        parser.error("--min-identity must be in (0, 100]")
    if os.path.realpath(args.input) == os.path.realpath(args.output):
        parser.error("--input and --output must be different paths")
    if args.resume and not (args.workdir or args.keep_work):
        parser.error("--resume requires --workdir or --keep-work")
    return args


def main(argv: Optional[Sequence[str]] = None) -> int:
    args = parse_args(argv)
    logging.basicConfig(
        level=getattr(logging, args.log_level),
        format="[%(asctime)s] %(levelname)s: %(message)s",
        datefmt="%Y-%m-%d %H:%M:%S",
    )
    generated_workdir = False
    if args.workdir:
        args.active_workdir = os.path.abspath(args.workdir)
        Path(args.active_workdir).mkdir(parents=True, exist_ok=True)
        LOG.info("Persistent work directory: %s", args.active_workdir)
    elif args.keep_work:
        args.active_workdir = os.path.abspath(args.output + ".work")
        Path(args.active_workdir).mkdir(parents=True, exist_ok=True)
        LOG.info("Persistent work directory: %s", args.active_workdir)
    else:
        output_path = Path(args.output).resolve()
        output_path.parent.mkdir(parents=True, exist_ok=True)
        args.active_workdir = tempfile.mkdtemp(
            prefix=f".{output_path.name}.work.",
            dir=str(output_path.parent),
        )
        generated_workdir = True
        LOG.info("Temporary work directory: %s", args.active_workdir)
    succeeded = False
    try:
        run(args)
        succeeded = True
    except subprocess.CalledProcessError as error:
        if error.stderr:
            LOG.error("%s", error.stderr.rstrip())
        LOG.error("External command failed with exit code %d: %s", error.returncode, error.cmd)
    except (OSError, RuntimeError, ValueError) as error:
        LOG.error("%s", error)
    finally:
        if generated_workdir and succeeded:
            shutil.rmtree(args.active_workdir)
            LOG.info("Removed successful temporary work directory: %s", args.active_workdir)
        elif generated_workdir:
            LOG.error("Retained failed-job work directory: %s", args.active_workdir)
    return 0 if succeeded else 1


if __name__ == "__main__":
    sys.exit(main())
