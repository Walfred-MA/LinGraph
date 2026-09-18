#!/usr/bin/env python3
"""Prepare PATs assemblies with WindowMasker and optional contig-name repair."""

from __future__ import annotations

import argparse
from concurrent.futures import ThreadPoolExecutor, as_completed
from dataclasses import dataclass
import gzip
import os
from pathlib import Path
import re
import shlex
import shutil
import subprocess
import sys
import tempfile
from typing import Dict, Iterable, List, Optional, Sequence, Tuple


SAMPLE_RE = re.compile(
    r"^(?P<sample>[A-Za-z0-9.-]+)_h(?P<haplotype>[1-9][0-9]*)$"
)
PREFIX_RE = re.compile(
    r"^(?P<sample>[A-Za-z0-9.-]+)#(?P<haplotype>[1-9][0-9]*)$"
)
LOWERCASE_DNA = frozenset(b"acgt")


@dataclass(frozen=True)
class Assembly:
    name: str
    fasta: Path
    output: Path


@dataclass(frozen=True)
class FastaInspection:
    records: int
    bases: int
    lowercase_bases: int
    names_correct: bool

    @property
    def is_masked(self) -> bool:
        return self.lowercase_bases > 0


@dataclass(frozen=True)
class PreparationResult:
    assembly: Assembly
    was_masked: bool
    masking_ran: bool
    names_were_correct: bool
    names_fixed: bool


def is_gzip(path: Path) -> bool:
    with path.open("rb") as handle:
        return handle.read(2) == b"\x1f\x8b"


def open_fasta(path: Path):
    if is_gzip(path):
        return gzip.open(path, "rb")
    return path.open("rb")


def expected_prefix(assembly_name: str) -> str:
    match = SAMPLE_RE.fullmatch(assembly_name)
    if match is None:
        raise ValueError(
            f"invalid assembly name {assembly_name!r}; expected "
            "SAMPLE_hHAPLOTYPE, for example HG002_h1"
        )
    return f"{match.group('sample')}#{match.group('haplotype')}"


def canonical_name(value: str) -> Tuple[str, str]:
    """Return ``(assembly-list name, FASTA prefix)`` for --name."""
    match = SAMPLE_RE.fullmatch(value)
    if match is None:
        match = PREFIX_RE.fullmatch(value)
    if match is None:
        raise ValueError(
            f"invalid --name {value!r}; expected SAMPLE_hHAPLOTYPE or "
            "SAMPLE#HAPLOTYPE, for example HG002_h1 or HG002#1"
        )
    sample = match.group("sample")
    haplotype = match.group("haplotype")
    return f"{sample}_h{haplotype}", f"{sample}#{haplotype}"


def inspect_fasta(path: Path, prefix: str) -> FastaInspection:
    """Validate a FASTA and detect lowercase masking and header conformance."""
    records = 0
    bases = 0
    lowercase_bases = 0
    names_correct = True
    seen_names = set()
    saw_sequence = False

    with open_fasta(path) as handle:
        for line_number, raw in enumerate(handle, 1):
            if raw.startswith(b">"):
                header = raw[1:].strip()
                if not header:
                    raise ValueError(f"{path}:{line_number}: empty FASTA header")
                identifier_bytes = header.split(None, 1)[0]
                try:
                    identifier = identifier_bytes.decode("utf-8")
                except UnicodeDecodeError as error:
                    raise ValueError(
                        f"{path}:{line_number}: FASTA identifier is not UTF-8"
                    ) from error
                if identifier in seen_names:
                    raise ValueError(
                        f"{path}:{line_number}: duplicate FASTA identifier "
                        f"{identifier!r}"
                    )
                seen_names.add(identifier)
                records += 1
                names_correct = names_correct and (
                    identifier == prefix or identifier.startswith(prefix + "#")
                )
                saw_sequence = True
                continue

            sequence = b"".join(raw.split())
            if not sequence:
                continue
            if not saw_sequence:
                raise ValueError(
                    f"{path}:{line_number}: sequence data occurs before a FASTA header"
                )
            bases += len(sequence)
            # bytes.count runs in C; a Python-level loop over a whole human
            # assembly makes the detection step unnecessarily expensive.
            lowercase_bases += sum(sequence.count(base) for base in LOWERCASE_DNA)

    if records == 0:
        raise ValueError(f"no FASTA records found: {path}")
    if bases == 0:
        raise ValueError(f"no FASTA sequence data found: {path}")
    return FastaInspection(records, bases, lowercase_bases, names_correct)


def normalize_record_name(record_name: str, prefix: str) -> str:
    """Return a record name with exactly one leading sample/haplotype prefix."""
    prefix_parts = prefix.split("#")
    record_parts = record_name.split("#")
    prefix_size = len(prefix_parts)
    leading_copies = 0
    while record_parts[:prefix_size] == prefix_parts:
        record_parts = record_parts[prefix_size:]
        leading_copies += 1
    if leading_copies == 1:
        return record_name
    return "#".join((*prefix_parts, *record_parts))


def fix_contig_names(source: Path, destination: Path, prefix: str) -> None:
    """Write an uncompressed FASTA with normalized, collision-free IDs."""
    seen_names = set()
    saw_header = False
    with open_fasta(source) as input_handle, destination.open("wb") as output_handle:
        for line_number, raw in enumerate(input_handle, 1):
            if not raw.startswith(b">"):
                output_handle.write(raw)
                continue

            saw_header = True
            header = raw[1:].rstrip(b"\r\n")
            fields = header.split(None, 1)
            if not fields:
                raise ValueError(f"{source}:{line_number}: empty FASTA header")
            try:
                identifier = fields[0].decode("utf-8")
            except UnicodeDecodeError as error:
                raise ValueError(
                    f"{source}:{line_number}: FASTA identifier is not UTF-8"
                ) from error
            normalized = normalize_record_name(identifier, prefix)
            if normalized in seen_names:
                raise ValueError(
                    f"{source}:{line_number}: contig-name repair would create "
                    f"duplicate identifier {normalized!r}"
                )
            seen_names.add(normalized)

            suffix = header[len(fields[0]) :]
            output_handle.write(b">" + normalized.encode("utf-8") + suffix + b"\n")

    if not saw_header:
        raise ValueError(f"no FASTA records found: {source}")


def copy_uncompressed(source: Path, destination: Path) -> None:
    with open_fasta(source) as input_handle, destination.open("wb") as output_handle:
        shutil.copyfileobj(input_handle, output_handle, length=1024 * 1024)


def run_command(command: Sequence[str]) -> None:
    print("[prepare] run:", shlex.join(command), flush=True)
    completed = subprocess.run(list(command))
    if completed.returncode != 0:
        raise RuntimeError(
            f"command exited with status {completed.returncode}: {shlex.join(command)}"
        )


def find_executable(value: str) -> str:
    if os.path.sep in value:
        candidate = Path(value).expanduser().resolve()
        if not candidate.is_file() or not os.access(candidate, os.X_OK):
            raise FileNotFoundError(f"executable not found or not executable: {candidate}")
        return str(candidate)
    resolved = shutil.which(value)
    if resolved is None:
        raise FileNotFoundError(f"required executable not found on PATH: {value}")
    return resolved


def run_windowmasker(
    windowmasker: str,
    source: Path,
    counts: Path,
    output: Path,
) -> None:
    run_command(
        [
            windowmasker,
            "-mk_counts",
            "-infmt",
            "fasta",
            "-in",
            str(source),
            "-out",
            str(counts),
        ]
    )
    run_command(
        [
            windowmasker,
            "-ustat",
            str(counts),
            "-dust",
            "yes",
            "-outfmt",
            "fasta",
            "-infmt",
            "fasta",
            "-in",
            str(source),
            "-out",
            str(output),
        ]
    )
    if not output.is_file() or output.stat().st_size == 0:
        raise RuntimeError(f"WindowMasker did not create a nonempty FASTA: {output}")


def output_basename(path: Path) -> str:
    name = path.name
    if name.lower().endswith(".gz"):
        name = name[:-3]
    if not name:
        raise ValueError(f"cannot derive an output filename from {path}")
    return name


def read_assembly_list(path: Path, output_folder: Path) -> List[Assembly]:
    assemblies: List[Assembly] = []
    seen_names: Dict[str, int] = {}
    seen_outputs: Dict[Path, Tuple[str, int]] = {}
    with path.open("rt", encoding="utf-8") as handle:
        for line_number, raw in enumerate(handle, 1):
            if not raw.strip() or raw.lstrip().startswith("#"):
                continue
            fields = shlex.split(raw, comments=True)
            if len(fields) != 2:
                raise ValueError(
                    f"{path}:{line_number}: expected exactly NAME FASTA; index must be FASTA.fai"
                )
            name = fields[0]
            expected_prefix(name)
            if name in seen_names:
                raise ValueError(
                    f"{path}:{line_number}: duplicate assembly name {name!r}; "
                    f"first used on line {seen_names[name]}"
                )
            fasta = Path(fields[1]).expanduser()
            fasta = fasta.resolve()
            if not fasta.is_file():
                raise FileNotFoundError(f"{path}:{line_number}: FASTA not found: {fasta}")

            output = (output_folder / output_basename(fasta)).resolve()
            if output in seen_outputs:
                previous_name, previous_line = seen_outputs[output]
                raise ValueError(
                    f"{path}:{line_number}: output filename collision between "
                    f"{previous_name!r} (line {previous_line}) and {name!r}: {output}"
                )
            if output == fasta:
                raise ValueError(
                    f"refusing to overwrite input FASTA {fasta}; choose a different "
                    "--output-folder"
                )

            seen_names[name] = line_number
            seen_outputs[output] = (name, line_number)
            assemblies.append(Assembly(name, fasta, output))

    if not assemblies:
        raise ValueError(f"{path}: no assembly records")
    return assemblies


def prepare_one(
    assembly: Assembly,
    inspection: FastaInspection,
    remask: bool,
    contignamefix: bool,
    windowmasker: Optional[str],
    samtools: str,
) -> PreparationResult:
    prefix = expected_prefix(assembly.name)
    masking_ran = remask or not inspection.is_masked
    names_fixed = contignamefix and not inspection.names_correct

    assembly.output.parent.mkdir(parents=True, exist_ok=True)
    with tempfile.TemporaryDirectory(
        prefix=f".{assembly.output.name}.{assembly.name}.",
        dir=str(assembly.output.parent),
    ) as temporary_name:
        temporary = Path(temporary_name)
        windowmasker_input = assembly.fasta
        if is_gzip(assembly.fasta):
            windowmasker_input = temporary / "input.fa"
            copy_uncompressed(assembly.fasta, windowmasker_input)

        if masking_ran:
            if windowmasker is None:
                raise RuntimeError("internal error: WindowMasker was not resolved")
            masked = temporary / "masked.fa"
            run_windowmasker(
                windowmasker,
                windowmasker_input,
                temporary / "windowmasker.counts",
                masked,
            )
            staged = masked
        else:
            staged = temporary / "copied.fa"
            copy_uncompressed(assembly.fasta, staged)

        # WindowMasker preserves FASTA IDs, so the input inspection determines
        # whether repair is needed. The resulting FASTA is inspected again below.
        if names_fixed:
            fixed = temporary / "namefixed.fa"
            fix_contig_names(staged, fixed, prefix)
            staged = fixed

        final_inspection = inspect_fasta(staged, prefix)
        if contignamefix and not final_inspection.names_correct:
            raise RuntimeError(
                f"contig-name repair failed for {assembly.name}: {assembly.fasta}"
            )

        os.replace(staged, assembly.output)

    index = Path(str(assembly.output) + ".fai")
    try:
        index.unlink()
    except FileNotFoundError:
        pass
    run_command([samtools, "faidx", str(assembly.output)])
    if not index.is_file() or index.stat().st_size == 0:
        raise RuntimeError(f"samtools did not create a nonempty index: {index}")

    return PreparationResult(
        assembly=assembly,
        was_masked=inspection.is_masked,
        masking_ran=masking_ran,
        names_were_correct=inspection.names_correct,
        names_fixed=names_fixed,
    )


def atomic_write_lines(path: Path, lines: Iterable[str]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    descriptor, temporary_name = tempfile.mkstemp(
        prefix=f".{path.name}.", suffix=".tmp", dir=str(path.parent), text=True
    )
    try:
        with os.fdopen(descriptor, "wt", encoding="utf-8") as handle:
            for line in lines:
                handle.write(line)
        os.replace(temporary_name, path)
    finally:
        try:
            os.unlink(temporary_name)
        except FileNotFoundError:
            pass


def parse_args(argv: Optional[Sequence[str]] = None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=(
            "Prepare one FASTA or every FASTA in a PATs assembly list with "
            "WindowMasker, optional contig-name repair, and samtools faidx."
        )
    )
    inputs = parser.add_mutually_exclusive_group(required=True)
    inputs.add_argument(
        "-q", "--assemblies",
        help="two-column SAMPLE_hHAPLOTYPE FASTA assembly list",
    )
    inputs.add_argument(
        "-i", "--input",
        help="one assembly FASTA; requires --name",
    )
    parser.add_argument(
        "-O", "--output-folder", required=True,
        help="folder for prepared FASTAs, indexes, and the rewritten list",
    )
    parser.add_argument(
        "-j", "--jobs", type=int, default=1,
        help="assemblies to prepare concurrently (default: 1)",
    )
    parser.add_argument(
        "--remask", action="store_true",
        help="run WindowMasker even when lowercase soft masking is detected",
    )
    parser.add_argument(
        "--contignamefix", action="store_true",
        help=(
            "detect headers that do not start SAMPLE#HAPLOTYPE and prefix them "
            "only when needed"
        ),
    )
    parser.add_argument(
        "--name",
        help=(
            "sample/haplotype name required by -i, or override for a one-entry "
            "assembly list; accepts HG002_h1 or HG002#1"
        ),
    )
    parser.add_argument(
        "--prepared-list",
        help="output assembly-list path (default: OUTPUT/query_paths.prepared.txt)",
    )
    parser.add_argument("--windowmasker", default="windowmasker")
    parser.add_argument("--samtools", default="samtools")
    args = parser.parse_args(argv)
    if args.jobs < 1:
        parser.error("--jobs must be positive")
    if args.input and not args.name:
        parser.error("--name is required with --input")
    return args


def main(argv: Optional[Sequence[str]] = None) -> int:
    args = parse_args(argv)
    output_folder = Path(args.output_folder).expanduser().resolve()
    output_folder.mkdir(parents=True, exist_ok=True)
    prepared_list = (
        Path(args.prepared_list).expanduser().resolve()
        if args.prepared_list
        else output_folder / "query_paths.prepared.txt"
    )

    if args.input:
        input_fasta = Path(args.input).expanduser().resolve()
        if not input_fasta.is_file():
            raise FileNotFoundError(f"assembly FASTA not found: {input_fasta}")
        name, _prefix = canonical_name(args.name)
        output = (output_folder / output_basename(input_fasta)).resolve()
        if output == input_fasta:
            raise ValueError(
                f"refusing to overwrite input FASTA {input_fasta}; choose a "
                "different --output-folder"
            )
        assemblies = [Assembly(name, input_fasta, output)]
    else:
        query_path = Path(args.assemblies).expanduser().resolve()
        if not query_path.is_file():
            raise FileNotFoundError(f"assembly list not found: {query_path}")
        assemblies = read_assembly_list(query_path, output_folder)
        if args.name:
            if len(assemblies) != 1:
                raise ValueError(
                    "--name can only be used when the assembly list has exactly "
                    "one entry"
                )
            name, _prefix = canonical_name(args.name)
            original = assemblies[0]
            assemblies = [Assembly(name, original.fasta, original.output)]
    inspections: Dict[str, FastaInspection] = {}
    inspection_errors: List[str] = []
    with ThreadPoolExecutor(max_workers=args.jobs) as executor:
        futures = {
            executor.submit(
                inspect_fasta, assembly.fasta, expected_prefix(assembly.name)
            ): assembly
            for assembly in assemblies
        }
        for future in as_completed(futures):
            assembly = futures[future]
            try:
                inspection = future.result()
            except Exception as error:
                inspection_errors.append(f"{assembly.name}: {error}")
                continue
            inspections[assembly.name] = inspection
            mask_state = "masked" if inspection.is_masked else "unmasked"
            name_state = "correct" if inspection.names_correct else "needs repair"
            print(
                f"[prepare] {assembly.name}: {mask_state}; contig names {name_state}; "
                f"{inspection.records} records",
                flush=True,
            )
    if inspection_errors:
        raise RuntimeError(
            "assembly inspection failed:\n  " + "\n  ".join(inspection_errors)
        )

    samtools = find_executable(args.samtools)
    needs_windowmasker = args.remask or any(
        not inspection.is_masked for inspection in inspections.values()
    )
    windowmasker = (
        find_executable(args.windowmasker) if needs_windowmasker else None
    )

    results: Dict[str, PreparationResult] = {}
    errors: List[str] = []
    with ThreadPoolExecutor(max_workers=args.jobs) as executor:
        futures = {
            executor.submit(
                prepare_one,
                assembly,
                inspections[assembly.name],
                args.remask,
                args.contignamefix,
                windowmasker,
                samtools,
            ): assembly
            for assembly in assemblies
        }
        for future in as_completed(futures):
            assembly = futures[future]
            try:
                result = future.result()
            except Exception as error:
                errors.append(f"{assembly.name}: {error}")
                continue
            results[assembly.name] = result
            mask_action = "WindowMasker run" if result.masking_ran else "mask retained"
            if result.names_fixed:
                name_action = "contig names repaired"
            elif args.contignamefix:
                name_action = "contig names already correct"
            else:
                name_action = "contig-name repair disabled"
            print(
                f"[prepare] complete {assembly.name}: {mask_action}; {name_action}; "
                f"{assembly.output}",
                flush=True,
            )

    if errors:
        raise RuntimeError("assembly preparation failed:\n  " + "\n  ".join(errors))

    atomic_write_lines(
        prepared_list,
        (
            f"{assembly.name}\t{assembly.output}\n"
            for assembly in assemblies
        ),
    )
    print(f"[prepare] prepared assembly list: {prepared_list}", flush=True)
    return 0


if __name__ == "__main__":
    try:
        raise SystemExit(main())
    except (OSError, ValueError, RuntimeError, subprocess.SubprocessError) as error:
        print(f"prepare_assemblies.py: error: {error}", file=sys.stderr)
        raise SystemExit(1)
