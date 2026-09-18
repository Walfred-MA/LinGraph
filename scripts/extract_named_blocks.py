#!/usr/bin/env python3
"""Extract BED blocks and name them in original source coordinates."""

from __future__ import annotations

import argparse
from contextlib import ExitStack
import gzip
import os
import re
import sys
from collections import defaultdict
from dataclasses import dataclass
from pathlib import Path
from typing import Dict, Iterable, Iterator, List, TextIO, Tuple


SOURCE_PATTERN = re.compile(
    r"^(?P<sample>[^:]+):(?P<contig>.+):"
    r"(?P<start>[0-9]+)-(?P<end>[0-9]+)(?P<strand>[+-])?$"
)


@dataclass(frozen=True)
class BedBlock:
    record: str
    start: int
    end: int
    group: str
    line_number: int


@dataclass(frozen=True)
class SourceMapping:
    sample: str
    contig: str
    start: int
    end: int
    strand: str


def open_text(path: Path) -> TextIO:
    if path.name.lower().endswith((".gz", ".bgz", ".bgzf")):
        return gzip.open(path, "rt", encoding="ascii")
    return path.open("rt", encoding="ascii")


def normalize_group(value: str, line_number: int) -> str:
    value = value.replace("_", "").replace("-", "")
    if re.fullmatch(r"g[0-9]+", value):
        return value
    if value.isdigit():
        return f"g{value}"
    if value in {"", "."}:
        return f"g{line_number}"
    if re.fullmatch(r"[A-Za-z0-9.]+", value):
        return value
    raise ValueError(
        f"BED line {line_number}: column 4 must contain only letters, "
        f"numbers, or dots (underscores and hyphens are removed), not "
        f"{value!r}"
    )


def read_bed(path: Path) -> Tuple[
    Dict[str, List[BedBlock]], int
]:
    by_record: Dict[str, List[BedBlock]] = defaultdict(list)
    count = 0
    with path.open("rt", encoding="utf-8") as handle:
        for line_number, line in enumerate(handle, 1):
            if not line.strip() or line.startswith("#"):
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 4:
                fields = line.split()
            if len(fields) < 4:
                raise ValueError(
                    f"{path}:{line_number}: expected at least four BED fields"
                )
            start = int(fields[1])
            end = int(fields[2])
            if start < 0 or end <= start:
                raise ValueError(
                    f"{path}:{line_number}: invalid interval {start}-{end}"
                )
            group = normalize_group(fields[3], line_number)
            by_record[fields[0]].append(
                BedBlock(fields[0], start, end, group, line_number)
            )
            count += 1
    if count == 0:
        raise ValueError(f"{path}: no BED blocks")
    return by_record, count


def fasta_records(path: Path) -> Iterator[Tuple[str, str]]:
    header = None
    sequence_parts: List[str] = []
    with open_text(path) as handle:
        for line_number, line in enumerate(handle, 1):
            if line.startswith(">"):
                if header is not None:
                    yield header, "".join(sequence_parts)
                header = line[1:].rstrip("\r\n")
                if not header:
                    raise ValueError(f"{path}:{line_number}: empty FASTA header")
                sequence_parts = []
            else:
                if header is None:
                    if line.strip():
                        raise ValueError(
                            f"{path}:{line_number}: sequence before header"
                        )
                    continue
                # Preserve uppercase/lowercase exactly. Calling .upper() here
                # destroys soft masking and makes repeats look unmasked.
                sequence_parts.append("".join(line.split()))
    if header is not None:
        yield header, "".join(sequence_parts)


def parse_source_mapping(
    header: str,
    fallback_sample: str = "",
    record_length: int = 0,
) -> SourceMapping:
    fields = header.split()
    record = fields[0]
    candidates = []
    for field in fields[1:]:
        if field.startswith("source="):
            candidates.insert(0, field[len("source="):])
        elif "=" not in field:
            candidates.append(field)

    for candidate in candidates:
        match = SOURCE_PATTERN.fullmatch(candidate)
        if match is None:
            continue
        start = int(match.group("start"))
        end = int(match.group("end"))
        if end <= start:
            raise ValueError(
                f"FASTA record {record}: invalid source mapping {candidate!r}"
            )
        return SourceMapping(
            sample=match.group("sample"),
            contig=match.group("contig"),
            start=start,
            end=end,
            strand=match.group("strand") or "+",
        )
    if fallback_sample:
        if record_length <= 0:
            raise ValueError(
                f"FASTA record {record}: cannot use --sample for an empty record"
            )
        return SourceMapping(
            sample=fallback_sample,
            contig=record,
            start=0,
            end=record_length,
            strand="+",
        )
    raise ValueError(
        f"FASTA record {record}: no SAMPLE:CONTIG:START-END[+-] "
        "or source=SAMPLE:CONTIG:START-END mapping; pass --sample for an "
        "ordinary assembly FASTA"
    )


def safe_name_component(value: str) -> str:
    return re.sub(r"[^A-Za-z0-9._#-]+", "_", value)


def source_coordinates(
    mapping: SourceMapping, local_start: int, local_end: int
) -> Tuple[int, int]:
    source_length = mapping.end - mapping.start
    if local_end > source_length:
        raise ValueError(
            f"local interval {local_start}-{local_end} exceeds mapped "
            f"source length {source_length}"
        )
    if mapping.strand == "-":
        return mapping.end - local_end, mapping.end - local_start
    return mapping.start + local_start, mapping.start + local_end


def wrap_sequence(sequence: str, width: int) -> Iterable[str]:
    for start in range(0, len(sequence), width):
        yield sequence[start:start + width]


def parse_arguments() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=(
            "Extract BED intervals from merged/novel FASTA records and name "
            "them with original sample, contig, and zero-based coordinates."
        )
    )
    parser.add_argument("-f", "--fasta", required=True, type=Path)
    parser.add_argument("-b", "--bed", required=True, type=Path)
    parser.add_argument("-o", "--output", required=True, type=Path)
    parser.add_argument(
        "--mapping-bed", type=Path,
        help=(
            "also write an exact 16-column reference-to-block mapping BED; "
            "used to bypass alignment for a supplied pre-grouped BED"
        ),
    )
    parser.add_argument(
        "--sample",
        default="",
        help=(
            "source haplotype for ordinary assembly FASTAs whose headers do "
            "not include SAMPLE:CONTIG:START-END metadata"
        ),
    )
    parser.add_argument(
        "--line-width", type=int, default=80,
        help="output FASTA sequence width (default: 80)",
    )
    return parser.parse_args()


def main() -> int:
    arguments = parse_arguments()
    if arguments.line_width <= 0:
        raise ValueError("--line-width must be positive")

    by_record, expected_blocks = read_bed(arguments.bed)
    arguments.output.parent.mkdir(parents=True, exist_ok=True)
    temporary = Path(str(arguments.output) + ".tmp")
    mapping_temporary = None
    if arguments.mapping_bed is not None:
        if arguments.mapping_bed.resolve() == arguments.output.resolve():
            raise ValueError("--mapping-bed must differ from --output")
        arguments.mapping_bed.parent.mkdir(parents=True, exist_ok=True)
        mapping_temporary = Path(str(arguments.mapping_bed) + ".tmp")
    written = 0
    seen_names = set()

    try:
        with ExitStack() as stack:
            output = stack.enter_context(temporary.open("wt", encoding="ascii"))
            mapping_output = (
                stack.enter_context(mapping_temporary.open("wt", encoding="utf-8"))
                if mapping_temporary is not None else None
            )
            for full_header, sequence in fasta_records(arguments.fasta):
                record = full_header.split()[0]
                blocks = by_record.pop(record, None)
                if not blocks:
                    continue
                mapping = parse_source_mapping(
                    full_header, arguments.sample, len(sequence)
                )
                mapping_length = mapping.end - mapping.start
                if len(sequence) != mapping_length:
                    raise ValueError(
                        f"FASTA record {record}: sequence length "
                        f"{len(sequence)} differs from source mapping length "
                        f"{mapping_length}"
                    )

                for block in blocks:
                    if block.end > len(sequence):
                        raise ValueError(
                            f"{arguments.bed}:{block.line_number}: interval "
                            f"exceeds FASTA record {record} length "
                            f"{len(sequence)}"
                        )
                    source_start, source_end = source_coordinates(
                        mapping, block.start, block.end
                    )
                    name = (
                        f"{block.group}_"
                        f"{safe_name_component(mapping.sample)}_"
                        f"{safe_name_component(mapping.contig)}_"
                        f"{source_start}_{source_end}"
                    )
                    if name in seen_names:
                        raise ValueError(f"duplicate output name {name}")
                    seen_names.add(name)
                    output.write(
                        f">{name} source={mapping.sample}:{mapping.contig}:"
                        f"{source_start}-{source_end}{mapping.strand}"
                        + (" sequence_role=imported_alternative"
                           if "sequence_role=imported_alternative" in full_header.split() else "")
                        + "\n"
                    )
                    subsequence = sequence[block.start:block.end]
                    for sequence_line in wrap_sequence(
                        subsequence, arguments.line_width
                    ):
                        output.write(sequence_line + "\n")
                    if mapping_output is not None:
                        unmasked = sum(base in "ACGT" for base in subsequence)
                        required = min(1000.0, 0.5 * unmasked)
                        query_sample = mapping.sample
                        mapping_output.write("\t".join((
                            block.record,
                            str(block.start),
                            str(block.end),
                            name,
                            "1000",
                            "+",
                            query_sample,
                            "0",
                            str(len(subsequence)),
                            str(unmasked),
                            str(unmasked),
                            f"{required:.1f}",
                            "100.000000",
                            "1",
                            "0",
                            "0",
                        )) + "\n")
                    written += 1

        if by_record:
            examples = ", ".join(sorted(by_record)[:5])
            raise ValueError(
                f"{len(by_record)} BED record names were absent from FASTA; "
                f"examples: {examples}"
            )
        if written != expected_blocks:
            raise RuntimeError(
                f"wrote {written} of {expected_blocks} BED blocks"
            )
        os.replace(temporary, arguments.output)
        if mapping_temporary is not None and arguments.mapping_bed is not None:
            os.replace(mapping_temporary, arguments.mapping_bed)
    except Exception:
        for path in (temporary, mapping_temporary):
            if path is None:
                continue
            try:
                path.unlink()
            except FileNotFoundError:
                pass
        raise

    print(
        f"records={written} output={arguments.output}"
        + (
            f" mapping={arguments.mapping_bed}"
            if arguments.mapping_bed is not None else ""
        ),
        file=sys.stderr,
    )
    return 0


if __name__ == "__main__":
    try:
        raise SystemExit(main())
    except (OSError, ValueError, RuntimeError) as error:
        print(f"extract_named_blocks.py: error: {error}", file=sys.stderr)
        raise SystemExit(1)
