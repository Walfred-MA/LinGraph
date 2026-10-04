#!/usr/bin/env python3
"""Add one sample/haplotype prefix to FASTA contig identifiers."""

from __future__ import annotations

import argparse
import gzip
import os
from pathlib import Path
import re
import sys
import tempfile
from typing import Optional, Sequence, Tuple


LIST_NAME = re.compile(r"^(?P<sample>[A-Za-z0-9.-]+)_h(?P<hap>[1-9][0-9]*)$")
FASTA_PREFIX = re.compile(r"^(?P<sample>[A-Za-z0-9.-]+)#(?P<hap>[1-9][0-9]*)$")
HAPLOTYPE_PREFIX = re.compile(r"^([A-Za-z0-9.-]+#[1-9][0-9]*)#")


def canonical_prefix(name: str) -> str:
    """Accept either an assembly-list label or its FASTA prefix."""
    value = name.strip()
    match = LIST_NAME.fullmatch(value)
    if match:
        return f"{match.group('sample')}#{match.group('hap')}"
    if FASTA_PREFIX.fullmatch(value):
        return value
    raise ValueError(
        f"invalid sample/haplotype name {name!r}; use SAMPLE_hN or SAMPLE#N, "
        "for example HG002_h1 or HG002#1"
    )


def normalize_record_name(record_name: str, prefix: str) -> str:
    """Return a record name with exactly one leading prefix."""
    prefix_parts = prefix.split("#")
    record_parts = record_name.split("#")
    size = len(prefix_parts)
    copies = 0
    while record_parts[:size] == prefix_parts:
        record_parts = record_parts[size:]
        copies += 1
    if copies == 1:
        return record_name
    return "#".join((*prefix_parts, *record_parts))


def check_haplotype_collision(record_name: str, prefix: str, source: Path, line: int) -> None:
    """Reject a different embedded haplotype prefix instead of nesting it."""
    remainder = record_name[len(prefix) + 1:] if record_name.startswith(prefix + "#") else record_name
    match = HAPLOTYPE_PREFIX.match(remainder)
    if match and match.group(1) != prefix:
        raise ValueError(
            f"{source}:{line}: contig {record_name!r} already contains the "
            f"different haplotype prefix {match.group(1)!r}; expected {prefix!r}. "
            "Correct the FASTA header or assembly-list name first; adding another "
            "prefix would create an ambiguous contig name."
        )


def is_gzip(path: Path) -> bool:
    with path.open("rb") as handle:
        return handle.read(2) == b"\x1f\x8b"


def default_output_path(source: Path) -> Path:
    """Return a sibling FASTA path that is safe to use as a new output."""
    name = source.name
    for suffix in (".bgzf", ".bgz", ".gz"):
        if name.lower().endswith(suffix):
            name = name[:-len(suffix)]
            break
    stem, extension = os.path.splitext(name)
    if not extension:
        extension = ".fa"
    return source.with_name(f"{stem}.namefixed{extension}")


def normalize_fasta(source: Path, destination: Path, prefix: str) -> Tuple[int, int]:
    """Write a prefixed FASTA atomically; return (renamed, record count)."""
    source = source.expanduser().resolve()
    destination = destination.expanduser().resolve()
    if source == destination:
        raise ValueError("input and output must be different files; the input is never overwritten")
    if not source.is_file():
        raise FileNotFoundError(f"input FASTA not found: {source}")

    destination.parent.mkdir(parents=True, exist_ok=True)
    descriptor, temporary_name = tempfile.mkstemp(
        prefix=f".{destination.name}.", suffix=".tmp", dir=destination.parent,
    )
    os.close(descriptor)
    temporary = Path(temporary_name)
    input_open = gzip.open if is_gzip(source) else open
    output_open = gzip.open if destination.name.lower().endswith((".gz", ".bgz", ".bgzf")) else open
    seen = set()
    renamed = 0
    records = 0
    try:
        with input_open(source, "rt", encoding="utf-8", newline="") as inp, \
                output_open(temporary, "wt", encoding="utf-8", newline="") as out:
            for line_number, line in enumerate(inp, 1):
                if not line.startswith(">"):
                    out.write(line)
                    continue

                header = line[1:].rstrip("\r\n")
                fields = header.split(None, 1)
                if not fields:
                    raise ValueError(f"{source}:{line_number}: empty FASTA header")
                original = fields[0]
                updated = normalize_record_name(original, prefix)
                check_haplotype_collision(updated, prefix, source, line_number)
                if updated in seen:
                    raise ValueError(
                        f"{source}:{line_number}: name fixing would create duplicate "
                        f"FASTA identifier {updated!r}; correct the input headers first"
                    )
                seen.add(updated)
                records += 1
                renamed += updated != original
                description = header[len(original):]
                out.write(f">{updated}{description}\n")

        if records == 0:
            raise ValueError(f"no FASTA records found: {source}")
        os.replace(temporary, destination)
        stale_index = Path(str(destination) + ".fai")
        if stale_index.is_file():
            stale_index.unlink()
        return renamed, records
    except BaseException:
        if temporary.exists():
            temporary.unlink()
        raise


def parse_args(argv: Optional[Sequence[str]] = None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Prefix FASTA contig names with a sample and haplotype.",
        epilog=(
            "Example: python3 tools/namecontigsfix.py -i raw.fa -n HG002_h1 "
            "writes raw.namefixed.fa. Then run samtools faidx on that output "
            "and use it in the assembly list."
        ),
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    parser.add_argument("-i", "--input", required=True, help="input FASTA (plain or gzip-compressed)")
    parser.add_argument(
        "-n", "--name", "--sample", dest="name", required=True,
        help="assembly name SAMPLE_hN or FASTA prefix SAMPLE#N, e.g. HG002_h1",
    )
    parser.add_argument(
        "-o", "--output",
        help="output FASTA (default: input name with .namefixed before its extension)",
    )
    return parser.parse_args(argv)


def main(argv: Optional[Sequence[str]] = None) -> int:
    args = parse_args(argv)
    source = Path(args.input).expanduser()
    prefix = canonical_prefix(args.name)
    destination = Path(args.output).expanduser() if args.output else default_output_path(source)
    renamed, records = normalize_fasta(source, destination, prefix)
    print(
        f"[namecontigsfix] wrote {destination.resolve()}: "
        f"{renamed}/{records} contig name(s) changed; prefix {prefix!r}"
    )
    print(f"[namecontigsfix] index with: samtools faidx {destination.resolve()}")
    return 0


if __name__ == "__main__":
    try:
        raise SystemExit(main())
    except Exception as error:
        print(f"namecontigsfix.py: error: {error}", file=sys.stderr)
        raise SystemExit(1)
