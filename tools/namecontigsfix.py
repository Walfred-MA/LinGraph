#!/usr/bin/env python3
"""Add one sample/haplotype prefix to FASTA record names, idempotently."""

from __future__ import annotations

import argparse
import gzip
import os
from pathlib import Path
import sys
import tempfile
from typing import Optional, Sequence


def normalize_record_name(record_name: str, prefix: str) -> str:
    """Return ``record_name`` with exactly one leading ``prefix``."""
    prefix_parts = [part for part in prefix.split("#") if part]
    if not prefix_parts:
        raise ValueError("--name must contain a nonempty FASTA prefix")

    record_parts = [part for part in record_name.split("#") if part]
    prefix_size = len(prefix_parts)
    leading_copies = 0
    while record_parts[:prefix_size] == prefix_parts:
        record_parts = record_parts[prefix_size:]
        leading_copies += 1

    # An already correct header is returned byte-for-byte at the identifier
    # level.  Older files containing a duplicated leading prefix are repaired.
    if leading_copies == 1:
        return record_name
    return "#".join((*prefix_parts, *record_parts))


def normalize_fasta(input_path: Path, output_path: Path, prefix: str) -> bool:
    """Normalize FASTA headers and return whether any header changed."""
    input_path = input_path.resolve()
    output_path = output_path.resolve()
    output_path.parent.mkdir(parents=True, exist_ok=True)

    descriptor, temporary_name = tempfile.mkstemp(
        prefix=f".{output_path.name}.",
        suffix=".tmp",
        dir=output_path.parent,
        text=True,
    )
    os.close(descriptor)
    temporary = Path(temporary_name)

    changed = False
    saw_header = False
    source_open = gzip.open if input_path.name.lower().endswith(".gz") else open
    output_open = gzip.open if output_path.name.lower().endswith(".gz") else open
    try:
        with source_open(input_path, "rt") as source, output_open(temporary, "wt") as output:
            for line_number, line in enumerate(source, 1):
                if not line.startswith(">"):
                    output.write(line)
                    continue

                saw_header = True
                header = line[1:].rstrip("\r\n")
                fields = header.split()
                if not fields:
                    raise ValueError(
                        f"{input_path}:{line_number}: empty FASTA header"
                    )
                normalized = normalize_record_name(fields[0], prefix)
                changed = changed or normalized != fields[0]
                if len(fields) == 1:
                    output.write(f">{normalized}\n")
                else:
                    description = "\t".join(fields[1:])
                    output.write(f">{normalized}\t{description}\n")

        if not saw_header:
            raise ValueError(f"no FASTA records found: {input_path}")

        output_replaced = changed or input_path != output_path
        if output_replaced:
            os.replace(temporary, output_path)
        else:
            temporary.unlink()
        if output_replaced:
            stale_index = Path(str(output_path) + ".fai")
            if stale_index.is_file():
                stale_index.unlink()
        return changed
    except BaseException:
        if temporary.exists():
            temporary.unlink()
        raise


def parse_args(argv: Optional[Sequence[str]] = None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Idempotently prefix FASTA contig names for PATs",
    )
    parser.add_argument("-i", "--input", required=True)
    parser.add_argument("-n", "--name", required=True)
    parser.add_argument("-o", "--output", required=True)
    return parser.parse_args(argv)


def main(argv: Optional[Sequence[str]] = None) -> int:
    args = parse_args(argv)
    changed = normalize_fasta(Path(args.input), Path(args.output), args.name)
    state = "updated" if changed else "already normalized"
    print(f"[namecontigsfix] {state}: {Path(args.output).resolve()}")
    return 0


if __name__ == "__main__":
    try:
        raise SystemExit(main())
    except Exception as error:
        print(f"namecontigsfix.py: error: {error}", file=sys.stderr)
        raise SystemExit(1)
