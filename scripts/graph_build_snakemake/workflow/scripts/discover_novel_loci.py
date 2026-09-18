#!/usr/bin/env python3
"""K-mer prefilter assemblies, then run the full minsetref novel-locus cleanup."""
from __future__ import annotations

import argparse
import hashlib
import json
import os
import re
import shutil
import subprocess
import sys
from pathlib import Path
from typing import Dict, List, Optional, Sequence

from pipeline_inputs import Assembly, atomic_write, read_assembly_list, read_fai


def file_fingerprint(path: str) -> Dict[str, object]:
    stat = os.stat(path)
    return {
        "path": os.path.realpath(path),
        "size": stat.st_size,
        "mtime_ns": stat.st_mtime_ns,
    }


def atomic_copy(source: str, destination: str) -> None:
    target = Path(destination)
    target.parent.mkdir(parents=True, exist_ok=True)
    temporary = target.with_name(target.name + f".tmp.{os.getpid()}")
    try:
        shutil.copyfile(source, temporary)
        os.replace(temporary, target)
    finally:
        try:
            temporary.unlink()
        except FileNotFoundError:
            pass


def run(command: Sequence[str]) -> None:
    print("[NOVEL] " + " ".join(command), flush=True)
    subprocess.run(list(command), check=True)


def residual_source_intervals(original, residual) -> Dict[str, List[tuple]]:
    """Recover original half-open coordinates encoded by KmerRefCleaner."""
    original_lengths = {name: original.length(name) for name in original.names()}
    intervals: Dict[str, List[tuple]] = {}
    for residual_name in residual.names():
        residual_length = residual.length(residual_name)
        if residual_name in original_lengths and residual_length == original_lengths[residual_name]:
            contig, start, end = residual_name, 0, residual_length
        else:
            match = re.fullmatch(r"(.+):([0-9]+)-([0-9]+)", residual_name)
            if match is None or match.group(1) not in original_lengths:
                raise ValueError(
                    f"cannot map KmerRefCleaner residual {residual_name!r} "
                    "back to its original FASTA record"
                )
            contig = match.group(1)
            start, end = int(match.group(2)), int(match.group(3))
            if start < 0 or end <= start or end > original_lengths[contig]:
                raise ValueError(f"invalid KmerRefCleaner residual interval {residual_name!r}")
            if end - start != residual_length:
                raise ValueError(
                    f"KmerRefCleaner residual {residual_name!r} has {residual_length} bp, "
                    f"expected {end - start}"
                )
        intervals.setdefault(contig, []).append((start, end))
    return intervals


def merge_context_intervals(
    intervals: Sequence[tuple], length: int, flank: int,
) -> List[tuple]:
    expanded = sorted(
        (max(0, start - flank), min(length, end + flank))
        for start, end in intervals
    )
    merged: List[List[int]] = []
    for start, end in expanded:
        if not merged or start > merged[-1][1]:
            merged.append([start, end])
        else:
            merged[-1][1] = max(merged[-1][1], end)
    return [(start, end) for start, end in merged]


def write_context_fasta(
    assembly: Assembly,
    residual_path: Path,
    output_path: Path,
    flank: int,
    indexed_fasta_class,
    wrap_fasta,
) -> None:
    """Write only residual neighborhoods, retaining enough source flank to lift."""
    temporary = output_path.with_name(output_path.name + f".tmp.{os.getpid()}")
    output_path.parent.mkdir(parents=True, exist_ok=True)
    try:
        with indexed_fasta_class(assembly.fasta, assembly.fai) as original, \
                indexed_fasta_class(str(residual_path)) as residual:
            intervals = residual_source_intervals(original, residual)
            with temporary.open("wt", encoding="ascii") as output:
                for contig in original.names():
                    if contig not in intervals:
                        continue
                    contexts = merge_context_intervals(
                        intervals[contig], original.length(contig), flank,
                    )
                    for index, (start, end) in enumerate(contexts, 1):
                        record = (
                            contig if start == 0 and end == original.length(contig)
                            else f"{contig}:context-{start}-{end}-{index:04d}"
                        )
                        output.write(
                            f">{record} source={assembly.name}:{contig}:{start}-{end}+\n"
                        )
                        for chunk_start in range(start, end, 8_000_000):
                            sequence = original.fetch(
                                contig, chunk_start, min(end, chunk_start + 8_000_000)
                            )
                            output.write(f"{wrap_fasta(sequence)}\n")
        os.replace(temporary, output_path)
    finally:
        try:
            temporary.unlink()
        except FileNotFoundError:
            pass


def parse_args(argv: Optional[Sequence[str]] = None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--fixed-reference", required=True)
    parser.add_argument("--assemblies", required=True)
    parser.add_argument("--working-directory", help="launch directory for relative assembly paths")
    parser.add_argument("--output", required=True)
    parser.add_argument(
        "--covered-bed",
        help=(
            "BED3 audit file for final novel 31-mer spans covered by the fixed "
            "reference/alternative catalog; default: beside --output"
        ),
    )
    parser.add_argument("--workdir", required=True)
    parser.add_argument("--kmer-ref-cleaner", required=True)
    parser.add_argument("--minsetref", required=True)
    parser.add_argument("--python", default=sys.executable)
    parser.add_argument("--samtools", default="samtools")
    parser.add_argument("--threads", type=int, required=True)
    parser.add_argument("--jobs", type=int, default=1)
    parser.add_argument(
        "--context-flank", type=int, default=30_000,
        help="original source flank retained around each k-mer residual (default: 30000)",
    )
    parser.add_argument("--exclude-name", action="append", default=[])
    parser.add_argument("--exclude-fasta", action="append", default=[])
    parser.add_argument("--skip-blastn", action="store_true")
    parser.add_argument("--keep-work", action="store_true")
    parser.add_argument(
        "--minimum-unmasked-bases", type=int, default=1000,
        help=(
            "discard final novel records with fewer remaining uppercase A/C/G/T "
            "bases after reference-kmer masking (default: 1000)"
        ),
    )
    args = parser.parse_args(argv)
    if (args.threads < 1 or args.jobs < 1 or args.context_flank < 10_000 or
            args.minimum_unmasked_bases < 0):
        parser.error(
            "--threads/--jobs must be positive, --context-flank must be >= 10000, "
            "and --minimum-unmasked-bases must be nonnegative"
        )
    return args


def validate_fixed(path: str) -> str:
    fasta = os.path.realpath(path)
    if not os.path.isfile(fasta):
        raise FileNotFoundError(fasta)
    if fasta.lower().endswith((".gz", ".bgz", ".bgzf")):
        raise ValueError(f"fixed reference must be uncompressed: {fasta}")
    fai = fasta + ".fai"
    if not os.path.isfile(fai):
        raise FileNotFoundError(
            f"fixed-reference index {fai} not found; run `samtools faidx {fasta}`"
        )
    read_fai(fai)
    return fasta


def selected_queries(args: argparse.Namespace) -> List[Assembly]:
    excluded_names = set(args.exclude_name)
    excluded_paths = {os.path.realpath(path) for path in args.exclude_fasta}
    selected = []
    for assembly in read_assembly_list(args.assemblies, args.working_directory):
        if assembly.name in excluded_names:
            continue
        if os.path.realpath(assembly.fasta) in excluded_paths:
            continue
        selected.append(assembly)
    return selected


def main(argv: Optional[Sequence[str]] = None) -> int:
    args = parse_args(argv)
    fixed = validate_fixed(args.fixed_reference)
    queries = selected_queries(args)
    workdir = Path(args.workdir).resolve()
    output = Path(args.output).resolve()
    covered_bed = (
        Path(args.covered_bed).resolve()
        if args.covered_bed
        else output.with_name(output.stem + ".reference_kmer_covered.bed")
    )
    if covered_bed in (output, Path(str(output) + ".fai")):
        raise ValueError("--covered-bed must differ from --output and its .fai")
    signature_payload = {
        "version": 2,
        "fixed": file_fingerprint(fixed),
        "fixed_fai": file_fingerprint(fixed + ".fai"),
        "assemblies": [
            {
                "name": row.name,
                "fasta": file_fingerprint(row.fasta),
                "fai": file_fingerprint(row.fai),
            }
            for row in queries
        ],
        "exclude_all": True,
        "skip_blastn": bool(args.skip_blastn),
        "context_flank": args.context_flank,
    }
    signature = hashlib.sha256(
        json.dumps(signature_payload, sort_keys=True, separators=(",", ":")).encode()
    ).hexdigest()
    signature_path = workdir / "input_signature.json"
    prior_signature = None
    if signature_path.is_file():
        try:
            prior_signature = json.loads(signature_path.read_text()).get("signature")
        except (OSError, ValueError, TypeError):
            prior_signature = None
    if workdir.exists() and prior_signature != signature:
        shutil.rmtree(workdir)
    workdir.mkdir(parents=True, exist_ok=True)
    atomic_write(
        str(signature_path),
        json.dumps(
            {"signature": signature, "inputs": signature_payload},
            indent=2,
            sort_keys=True,
        ) + "\n",
    )

    residual_dir = workdir / "kmer_residuals"
    residual_dir.mkdir(parents=True, exist_ok=True)
    input_list = workdir / "kmer_inputs.list"
    output_list = workdir / "kmer_outputs.list"
    residual_by_sample = {
        row.name: residual_dir / f"{row.name}.fa" for row in queries
    }
    atomic_write(
        str(input_list),
        "".join(f"{row.name}\t{row.fasta}\n" for row in queries),
    )
    atomic_write(
        str(output_list),
        "".join(
            f"{row.name}\t{residual_by_sample[row.name]}\n" for row in queries
        ),
    )
    kmer_marker = workdir / "kmer_prefilter.complete"
    residuals_complete = kmer_marker.is_file() and all(
        path.is_file() for path in residual_by_sample.values()
    )
    kmer_ran = False
    if queries and not residuals_complete:
        run([
            args.kmer_ref_cleaner,
            "-I", str(input_list),
            "-O", str(output_list),
            "-r", fixed,
            "--exclude-all",
            "-t", str(args.threads),
        ])
        atomic_write(str(kmer_marker), signature + "\n")
        kmer_ran = True

    nonempty = []
    source_root = Path(args.minsetref).resolve().parent
    if str(source_root) not in sys.path:
        sys.path.insert(0, str(source_root))
    from minsetref_core import IndexedFasta, wrap_fasta
    for row in queries:
        residual = residual_by_sample[row.name]
        if not residual.is_file():
            raise RuntimeError(f"KmerRefCleaner did not create {residual}")
        if residual.stat().st_size == 0:
            continue
        fai = Path(str(residual) + ".fai")
        if not fai.is_file() or fai.stat().st_mtime_ns < residual.stat().st_mtime_ns:
            run([args.samtools, "faidx", str(residual)])
        context = residual_dir / f"{row.name}.context.fa"
        context_fai = Path(str(context) + ".fai")
        if kmer_ran or not context.is_file() or not context_fai.is_file():
            write_context_fasta(
                row, residual, context, args.context_flank, IndexedFasta, wrap_fasta,
            )
            run([args.samtools, "faidx", str(context)])
        nonempty.append((row.name, context, context_fai))

    if not nonempty:
        output.parent.mkdir(parents=True, exist_ok=True)
        atomic_write(str(output), "")
        atomic_write(str(output) + ".fai", "")
        atomic_write(str(covered_bed), "")
        print("[NOVEL] K-mer prefilter removed every query sequence; novel_loci.fa is empty")
        return 0

    cleaned_list = workdir / "query_paths.kmer_cleaned.txt"
    cleaned_rows = [f"FIXEDCATALOG_h1\t{fixed}\n"]
    cleaned_rows.extend(
        f"{sample}\t{residual}\n"
        for sample, residual, fai in nonempty
    )
    atomic_write(str(cleaned_list), "".join(cleaned_rows))

    minset_out = workdir / "minsetref"
    command = [
        args.python, args.minsetref,
        "-l", str(cleaned_list),
        "-o", str(minset_out),
        "--precleaned-reference", fixed,
        "--stop-after-novel-loci",
        "--full-main-discovery-cleaning",
        "-j", str(min(args.jobs, len(nonempty))),
        "-t", str(max(1, args.threads // min(args.jobs, len(nonempty)))),
    ]
    if args.skip_blastn:
        command.append("--skip-blastn")
    if args.keep_work:
        command.append("--keep-work")
    if (minset_out / "checkpoint_format.json").is_file():
        command.append("--resume")
    elif minset_out.exists():
        shutil.rmtree(minset_out)
    run(command)

    discovered = minset_out / "novel_loci.fa"
    if not discovered.is_file():
        raise RuntimeError(f"minsetref did not create {discovered}")
    if discovered.stat().st_size == 0:
        atomic_write(str(output), "")
        atomic_write(str(output) + ".fai", "")
        atomic_write(str(covered_bed), "")
        print("[NOVEL] Alignment cleaning produced no novel loci")
        return 0

    discovered_fai = Path(str(discovered) + ".fai")
    if (not discovered_fai.is_file() or
            discovered_fai.stat().st_mtime_ns < discovered.stat().st_mtime_ns):
        run([args.samtools, "faidx", str(discovered)])

    masked = workdir / "novel_loci.reference_kmer_masked.fa"
    work_bed = workdir / "novel_loci.reference_kmer_covered.bed"
    run([
        args.kmer_ref_cleaner,
        "-i", str(discovered),
        "-o", str(masked),
        "-r", fixed,
        "--exclude-all",
        "--mask-covered",
        "--min-unmasked-bases", str(args.minimum_unmasked_bases),
        "--covered-bed", str(work_bed),
        "-t", str(args.threads),
    ])
    atomic_copy(str(masked), str(output))
    atomic_copy(str(work_bed), str(covered_bed))
    if output.stat().st_size:
        run([args.samtools, "faidx", str(output)])
    else:
        atomic_write(str(output) + ".fai", "")
    return 0


if __name__ == "__main__":
    try:
        raise SystemExit(main())
    except (OSError, ValueError, RuntimeError, subprocess.CalledProcessError) as error:
        print(f"discover_novel_loci.py: error: {error}", file=sys.stderr)
        raise SystemExit(1)
