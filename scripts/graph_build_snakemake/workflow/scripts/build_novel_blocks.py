#!/usr/bin/env python3
"""Make whole-record blocks for small novel loci and k-mer blocks for large loci."""
from __future__ import annotations

import argparse
import hashlib
import os
import subprocess
import sys
from pathlib import Path
from typing import Dict, Iterable, List, Optional, Sequence, Tuple

from pipeline_inputs import atomic_write


def read_fai_allow_empty(path: str) -> List[Tuple[str, int]]:
    records: List[Tuple[str, int]] = []
    seen = set()
    with open(path, "rt", encoding="utf-8") as handle:
        for line_number, raw in enumerate(handle, 1):
            if not raw.strip():
                continue
            fields = raw.rstrip("\n").split("\t")
            if len(fields) < 5:
                raise ValueError(f"{path}:{line_number}: expected FAI5+")
            name = fields[0]
            try:
                length = int(fields[1])
                tuple(int(value) for value in fields[2:5])
            except ValueError as error:
                raise ValueError(f"{path}:{line_number}: invalid FAI field") from error
            if not name or length <= 0 or name in seen:
                raise ValueError(f"{path}:{line_number}: invalid or duplicate record {name!r}")
            seen.add(name)
            records.append((name, length))
    return records


def extract_records(source: str, wanted: Iterable[str], output: str) -> None:
    wanted_set = set(wanted)
    observed = set()
    destination = Path(output)
    destination.parent.mkdir(parents=True, exist_ok=True)
    temporary = destination.with_name(destination.name + f".tmp.{os.getpid()}")
    keep = False
    try:
        with open(source, "rt", encoding="ascii") as handle, temporary.open(
            "wt", encoding="ascii"
        ) as out:
            for line_number, raw in enumerate(handle, 1):
                if raw.startswith(">"):
                    fields = raw[1:].strip().split()
                    if not fields:
                        raise ValueError(f"{source}:{line_number}: empty FASTA header")
                    name = fields[0]
                    keep = name in wanted_set
                    if keep:
                        observed.add(name)
                if keep:
                    out.write(raw if raw.endswith("\n") else raw + "\n")
        missing = wanted_set - observed
        if missing:
            raise ValueError(
                f"{source}: requested large novel records were absent: "
                + ", ".join(sorted(missing)[:5])
            )
        os.replace(temporary, destination)
    finally:
        try:
            temporary.unlink()
        except FileNotFoundError:
            pass


def read_auto_bed(path: str, lengths: Dict[str, int]) -> Dict[str, List[Tuple[int, int]]]:
    rows: Dict[str, List[Tuple[int, int]]] = {}
    with open(path, "rt", encoding="utf-8") as handle:
        for line_number, raw in enumerate(handle, 1):
            if not raw.strip() or raw.lstrip().startswith("#"):
                continue
            fields = raw.split()
            if len(fields) < 3:
                raise ValueError(f"{path}:{line_number}: expected BED3+")
            contig = fields[0]
            try:
                start, end = int(fields[1]), int(fields[2])
            except ValueError as error:
                raise ValueError(f"{path}:{line_number}: invalid coordinates") from error
            if contig not in lengths or start < 0 or end <= start or end > lengths[contig]:
                raise ValueError(f"{path}:{line_number}: invalid novel block interval")
            rows.setdefault(contig, []).append((start, end))
    return rows


def block_name(contig: str, index: int) -> str:
    digest = hashlib.sha256(contig.encode("utf-8")).hexdigest()[:12]
    return f"novel-{digest}-{index:05d}"


def run(command: Sequence[str]) -> None:
    print("[NOVEL-BLOCKS] " + " ".join(command), flush=True)
    subprocess.run(list(command), check=True)


def parse_args(argv: Optional[Sequence[str]] = None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--fasta", required=True)
    parser.add_argument("--fai", required=True)
    parser.add_argument("--output", required=True)
    parser.add_argument("--workdir", required=True)
    parser.add_argument("--samtools", default="samtools")
    parser.add_argument("--kmerindex", required=True)
    parser.add_argument("--kmerblocking", required=True)
    parser.add_argument("--threads", type=int, required=True)
    parser.add_argument("--large-threshold", type=int, default=1_000_000)
    args = parser.parse_args(argv)
    if args.threads < 1 or args.large_threshold < 1:
        parser.error("--threads and --large-threshold must be positive")
    return args


def main(argv: Optional[Sequence[str]] = None) -> int:
    args = parse_args(argv)
    records = read_fai_allow_empty(args.fai)
    if not records:
        atomic_write(args.output, "")
        return 0
    lengths = dict(records)
    large = [name for name, length in records if length > args.large_threshold]
    auto: Dict[str, List[Tuple[int, int]]] = {}
    if large:
        workdir = Path(args.workdir).resolve()
        workdir.mkdir(parents=True, exist_ok=True)
        large_fasta = workdir / "novel_loci.gt1Mb.fa"
        large_index = workdir / "novel_loci.gt1Mb.bin"
        raw_bed = workdir / "novel_loci.gt1Mb.raw.bed"
        extract_records(args.fasta, large, str(large_fasta))
        run([args.samtools, "faidx", str(large_fasta)])
        run([
            args.kmerindex, "-f", str(large_fasta), "-o", str(large_index),
            "-t", str(args.threads),
        ])
        run([
            args.kmerblocking, "-i", str(large_index), "-f", str(large_fasta) + ".fai",
            "-o", str(raw_bed), "-t", str(args.threads),
        ])
        auto = read_auto_bed(str(raw_bed), lengths)

    lines = []
    for contig, length in records:
        intervals = auto.get(contig) if length > args.large_threshold else None
        if not intervals:
            intervals = [(0, length)]
        for index, (start, end) in enumerate(intervals, 1):
            lines.append(
                f"{contig}\t{start}\t{end}\t{block_name(contig, index)}\n"
            )
    atomic_write(args.output, "".join(lines))
    return 0


if __name__ == "__main__":
    try:
        raise SystemExit(main())
    except (OSError, ValueError, RuntimeError, subprocess.CalledProcessError) as error:
        print(f"build_novel_blocks.py: error: {error}", file=sys.stderr)
        raise SystemExit(1)
