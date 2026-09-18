#!/usr/bin/env python3
"""Extract template-defined grouped BED blocks from cohort assemblies.

The input BED is the five-column output of ``normalize-template-bed``:
CONTIG, START, END, GROUP, TEMPLATE_HAPLOTYPE. Each interval is fetched from
its named assembly, so graph construction needs no global reference sample.
"""

from __future__ import annotations

import argparse
import os
from pathlib import Path
import sys
from typing import Optional, Sequence

from graph_build_snakemake.workflow.scripts.pipeline_inputs import (
    read_assembly_list,
)
from minsetref_core import IndexedFasta, sanitize_id, wrap_fasta


def parse_args(argv: Optional[Sequence[str]] = None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--query-paths", required=True)
    parser.add_argument("--bed", required=True)
    parser.add_argument("--output", required=True)
    parser.add_argument("--mapping-bed", required=True)
    return parser.parse_args(argv)


def run(args: argparse.Namespace) -> None:
    assemblies = {
        row.name: row for row in read_assembly_list(args.query_paths)
    }
    output = Path(args.output).expanduser().resolve()
    mapping = Path(args.mapping_bed).expanduser().resolve()
    output.parent.mkdir(parents=True, exist_ok=True)
    mapping.parent.mkdir(parents=True, exist_ok=True)
    output_tmp = Path(str(output) + f".tmp.{os.getpid()}")
    mapping_tmp = Path(str(mapping) + f".tmp.{os.getpid()}")
    readers = {}
    seen = set()
    written = 0
    try:
        with output_tmp.open("wt", encoding="ascii") as fasta, \
                mapping_tmp.open("wt", encoding="utf-8") as bed, \
                open(args.bed, "rt", encoding="utf-8") as source:
            for line_number, raw in enumerate(source, 1):
                if not raw.strip() or raw.lstrip().startswith("#"):
                    continue
                fields = raw.rstrip("\r\n").split("\t")
                if len(fields) != 5:
                    raise ValueError(
                        f"{args.bed}:{line_number}: expected normalized BED5"
                    )
                contig, start_text, end_text, group, sample = fields
                start, end = int(start_text), int(end_text)
                assembly = assemblies.get(sample)
                if assembly is None:
                    raise ValueError(
                        f"{args.bed}:{line_number}: unknown template {sample!r}"
                    )
                reader = readers.get(sample)
                if reader is None:
                    reader = IndexedFasta(assembly.fasta, assembly.fai)
                    readers[sample] = reader
                if contig not in reader.index:
                    raise KeyError(
                        f"{args.bed}:{line_number}: {contig!r} is absent from "
                        f"{sample}"
                    )
                sequence = reader.fetch(contig, start, end)
                record = (
                    f"{group}_{sanitize_id(sample)}_{sanitize_id(contig)}_"
                    f"{start}_{end}"
                )
                if record in seen:
                    raise ValueError(
                        f"{args.bed}:{line_number}: duplicate template block "
                        f"{record!r}"
                    )
                seen.add(record)
                unmasked = sum(base in "ACGT" for base in sequence)
                required = min(1000.0, 0.5 * unmasked)
                fasta.write(
                    f">{record} source={sample}:{contig}:{start}-{end}+\n"
                    f"{wrap_fasta(sequence)}\n"
                )
                bed.write("\t".join((
                    contig, str(start), str(end), record, "1000", "+",
                    sample, "0", str(len(sequence)), str(unmasked),
                    str(unmasked), f"{required:.1f}", "100.000000", "1",
                    "0", "0",
                )) + "\n")
                written += 1
        if written == 0:
            raise ValueError(f"{args.bed}: no template blocks")
        os.replace(output_tmp, output)
        os.replace(mapping_tmp, mapping)
    finally:
        for reader in readers.values():
            reader.close()
        for temporary in (output_tmp, mapping_tmp):
            try:
                temporary.unlink()
            except FileNotFoundError:
                pass


def main(argv: Optional[Sequence[str]] = None) -> int:
    args = parse_args(argv)
    run(args)
    return 0


if __name__ == "__main__":
    try:
        raise SystemExit(main())
    except (OSError, ValueError, RuntimeError) as error:
        print(
            f"extract_grouped_template_blocks.py: error: {error}",
            file=sys.stderr,
        )
        raise SystemExit(1)
