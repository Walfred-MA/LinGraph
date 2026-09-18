#!/usr/bin/env python3
"""Input validation and small orchestration helpers for the graph workflow."""
from __future__ import annotations

import argparse
from collections import Counter
from concurrent.futures import ThreadPoolExecutor, as_completed
import fcntl
import hashlib
import json
import os
import re
import shlex
import shutil
import subprocess
import sys
import tempfile
import threading
import time
from dataclasses import dataclass
from pathlib import Path
from typing import Dict, Iterable, List, Optional, Sequence, Tuple


SAMPLE_RE = re.compile(r"^(?P<sample>[A-Za-z0-9.-]+)_h(?P<haplotype>[1-9][0-9]*)$")
SAFE_BLOCK_RE = re.compile(r"^[A-Za-z0-9.]+$")
UNSAFE_BLOCK_RE = re.compile(r"[^A-Za-z0-9.]+")
FASTA_SUFFIX_RE = re.compile(r"(?:\.fasta|\.fna|\.fa)(?:\.gz)?$", re.IGNORECASE)
SOURCE_COORD_TOKEN_RE = re.compile(
    r"^source=[^:]+:.+:(?P<start>[0-9]+)-(?P<end>[0-9]+)"
    r"(?P<strand>[+-]?)$"
)
MAX_BLOCK_NAME_LENGTH = 96
MINSET_KMER_PRECLEAN_PROTOCOL = "kmer-preclean-v2-blast17-identifiers\n"
MINSET_MAX_RIGOROUS_PROTOCOL = (
    "max-rigorous-minimap2-first-v2-blast17-identifiers\n"
)
INTEGRATED_PARTITION_PROTOCOL = "integrated-partition-v1\n"
PARTITION_TIMEOUT_EXIT_CODE = 124
DEFAULT_ALIGNMENT_TIMEOUT = 21_600
_FAILURE_LIST_THREAD_LOCK = threading.Lock()


def integrated_partition_protocol(args: argparse.Namespace) -> str:
    if getattr(args, "static_block", False):
        dynamic_args = argparse.Namespace(**{**vars(args), "static_block": False})
        return "static-block-v1\n" + integrated_partition_protocol(dynamic_args)
    fast_mode = bool(getattr(args, "fast_mode", False))
    slow_rigorous = bool(getattr(args, "slow_rigorous", False))
    max_rigorous = bool(getattr(args, "max_rigorous", False))
    linear = bool(getattr(args, "linear", True))
    template_defined = bool(getattr(args, "template_defined_reference", False))
    if (
        not fast_mode and not slow_rigorous and not max_rigorous
        and linear and not template_defined
    ):
        # Preserve resume compatibility with every pre-fast-mode partition.
        return INTEGRATED_PARTITION_PROTOCOL
    if fast_mode:
        # Preserve the exact marker written by historical --fast-mode runs.
        return (
            "integrated-partition-v3;"
            f"fast=1;linear={int(linear)};"
            f"template-defined={int(template_defined)};"
            "max-unmapped-gap=50\n"
        )
    alignment_mode = (
        "max-rigorous" if max_rigorous else
        "slow-rigorous" if slow_rigorous else
        "legacy"
    )
    return (
        "integrated-partition-v4;"
        f"alignment-mode={alignment_mode};linear={int(linear)};"
        f"template-defined={int(template_defined)};"
        f"max-unmapped-gap={50 if fast_mode else 'not-applicable'}\n"
    )


def minset_completion_marker(
    unique_bed: Path, args: argparse.Namespace,
) -> Tuple[Path, str]:
    """Return the mode-specific minset completion marker and protocol."""
    if bool(getattr(args, "max_rigorous", False)):
        return (
            Path(str(unique_bed) + ".max-rigorous-v2.complete"),
            MINSET_MAX_RIGOROUS_PROTOCOL,
        )
    return (
        Path(str(unique_bed) + ".kmer-preclean-v2.complete"),
        MINSET_KMER_PRECLEAN_PROTOCOL,
    )


@dataclass(frozen=True)
class Assembly:
    name: str
    fasta: str
    fai: str


def atomic_write(path: str, text: str) -> None:
    destination = Path(path)
    destination.parent.mkdir(parents=True, exist_ok=True)
    temporary = destination.with_name(destination.name + f".tmp.{os.getpid()}")
    try:
        temporary.write_text(text, encoding="utf-8")
        os.replace(temporary, destination)
    finally:
        try:
            temporary.unlink()
        except FileNotFoundError:
            pass


def block_policy_inputs(graph: Path, static_block: bool) -> List[str]:
    """Track mode changes without invalidating pre-option dynamic builds."""
    path = graph / "inputs" / "block_policy.txt"
    if not static_block and not path.exists():
        return []
    value = "static-block-v1\n" if static_block else "dynamic-block-v1\n"
    if not path.is_file() or path.read_text() != value:
        atomic_write(str(path), value)
    return [str(path)]


def partition_failure_list(args: argparse.Namespace) -> Path:
    configured = str(getattr(args, "fail_list", "") or "").strip()
    if configured:
        return Path(configured).expanduser().resolve()
    return Path(args.output_dir).expanduser().resolve() / "fail.list"


def update_partition_failure(
    path: Path, partition: str, reason: Optional[str], diagnostic: str = "",
) -> None:
    """Upsert or clear one partition in a process-safe cohort failure list."""
    if reason is None and not path.exists():
        return
    path.parent.mkdir(parents=True, exist_ok=True)
    lock_path = Path(str(path) + ".lock")
    with _FAILURE_LIST_THREAD_LOCK:
        with lock_path.open("a+") as lock_handle:
            fcntl.flock(lock_handle.fileno(), fcntl.LOCK_EX)
            try:
                failures: Dict[str, Tuple[str, str]] = {}
                try:
                    with path.open("rt", encoding="utf-8") as input_handle:
                        for raw in input_handle:
                            fields = raw.rstrip("\r\n").split("\t")
                            if not fields or fields[0] in {"", "partition"}:
                                continue
                            failures[fields[0]] = (
                                fields[1] if len(fields) > 1 else "failed",
                                fields[2] if len(fields) > 2 else "",
                            )
                except FileNotFoundError:
                    pass
                if reason is None:
                    failures.pop(partition, None)
                else:
                    clean_reason = reason.replace("\t", " ").replace("\n", " ")
                    clean_diagnostic = (
                        diagnostic.replace("\t", " ").replace("\n", " ")
                    )
                    failures[partition] = (clean_reason, clean_diagnostic)
                rows = ["partition\treason\tcheckpoint"]
                rows.extend(
                    f"{name}\t{values[0]}\t{values[1]}"
                    for name, values in sorted(failures.items())
                )
                atomic_write(str(path), "\n".join(rows) + "\n")
            finally:
                fcntl.flock(lock_handle.fileno(), fcntl.LOCK_UN)


def run_partition_graph_command(
    args: argparse.Namespace, partition: str, graph_dir: Path,
    command: Sequence[str],
) -> bool:
    """Run one graph command, converting only status 124 into cohort skip."""
    completed = subprocess.run(list(command), check=False)
    if completed.returncode == 0:
        update_partition_failure(
            partition_failure_list(args), partition, None,
        )
        return True
    if completed.returncode == PARTITION_TIMEOUT_EXIT_CODE:
        timeout = int(getattr(
            args, "alignment_timeout", DEFAULT_ALIGNMENT_TIMEOUT,
        ))
        checkpoint = graph_dir / f".{partition}_align.txt.direct.work"
        failures = partition_failure_list(args)
        update_partition_failure(
            failures, partition, f"alignment_timeout_{timeout}s",
            str(checkpoint),
        )
        print(
            f"[WARNING] {partition}: alignment exceeded {timeout} seconds; "
            f"skipping this partition, retaining resumable checkpoints at "
            f"{checkpoint}, and recording it in {failures}",
            file=sys.stderr, flush=True,
        )
        return False
    raise subprocess.CalledProcessError(completed.returncode, list(command))


def preparation_message(folder: str) -> str:
    command = os.path.join(os.path.abspath(folder), "prepare_assemblies.py")
    return (
        "Prepare or repair the assemblies with "
        f"`python3 {command} -q ASSEMBLIES -O PREPARED --contignamefix`, "
        "then rerun this workflow with PREPARED/query_paths.prepared.txt."
    )


def read_assembly_list(path: str, working_directory=None) -> List[Assembly]:
    source = Path(path).expanduser().resolve()
    base = Path(working_directory) if working_directory else Path.cwd()
    assemblies: List[Assembly] = []
    seen: Dict[str, int] = {}
    with source.open("rt", encoding="utf-8") as handle:
        for line_number, raw in enumerate(handle, 1):
            if not raw.strip() or raw.lstrip().startswith("#"):
                continue
            fields = shlex.split(raw, comments=True)
            if len(fields) != 2:
                raise ValueError(
                    f"{source}:{line_number}: expected exactly NAME FASTA; index must be FASTA.fai"
                )
            name = fields[0]
            if name in seen:
                raise ValueError(
                    f"{source}:{line_number}: duplicate assembly name {name!r}; "
                    f"first used on line {seen[name]}"
                )
            fasta = Path(fields[1]).expanduser()
            if not fasta.is_absolute():
                fasta = base / fasta
            fasta = fasta.resolve()
            if not fasta.is_file():
                raise FileNotFoundError(f"{source}:{line_number}: {fasta} not found")
            fai = Path(str(fasta) + ".fai")
            if not fai.is_file():
                raise FileNotFoundError(
                    f"{source}:{line_number}: FASTA index {fai} not found; "
                    f"run `samtools faidx {fasta}` before starting the workflow"
                )
            seen[name] = line_number
            assemblies.append(Assembly(name, str(fasta), str(fai)))
    if not assemblies:
        raise ValueError(f"{source}: no assembly records")
    return assemblies


def first_fasta_name(path: str) -> str:
    if path.lower().endswith((".gz", ".bgz", ".bgzf")):
        raise ValueError(f"graph construction requires an uncompressed FASTA: {path}")
    with open(path, "rt", encoding="ascii") as handle:
        for raw in handle:
            if raw.startswith(">"):
                name = raw[1:].strip().split()[0]
                if name:
                    return name
                break
    raise ValueError(f"{path}: no FASTA records")


def derive_reference_name(path: str, existing: Iterable[str]) -> str:
    first = first_fasta_name(path)
    fields = first.split("#")
    if len(fields) >= 2 and fields[0] and fields[1].isdigit() and int(fields[1]) > 0:
        candidate = f"{fields[0]}_h{int(fields[1])}"
    else:
        stem = FASTA_SUFFIX_RE.sub("", Path(path).name)
        stem = re.sub(r"[^A-Za-z0-9.-]+", "-", stem).strip(".-") or "reference"
        candidate = stem if SAMPLE_RE.fullmatch(stem) else f"{stem}_h1"
    used = set(existing)
    if candidate not in used:
        return candidate
    digest = hashlib.sha1(os.path.realpath(path).encode()).hexdigest()[:8]
    base = SAMPLE_RE.fullmatch(candidate).group("sample")  # type: ignore[union-attr]
    hap = SAMPLE_RE.fullmatch(candidate).group("haplotype")  # type: ignore[union-attr]
    return f"{base}-{digest}_h{hap}"


def resolve_reference_and_assemblies(query_paths: str, reference: str, working_directory=None) -> List[Assembly]:
    assemblies = read_assembly_list(query_paths, working_directory)
    by_name = {row.name: row for row in assemblies}
    reference_text = str(reference).strip()
    if not reference_text:
        raise ValueError("reference cannot be empty")
    if reference_text in by_name:
        selected = by_name[reference_text]
    else:
        candidate = Path(reference_text).expanduser()
        if not candidate.is_absolute():
            candidate = Path(working_directory or Path.cwd()) / candidate
        candidate = candidate.resolve()
        if not candidate.is_file():
            raise ValueError(
                f"reference {reference_text!r} is neither an assembly-list name "
                "nor an existing FASTA"
            )
        matches = [row for row in assemblies if os.path.samefile(row.fasta, candidate)]
        if matches:
            selected = matches[0]
        else:
            fai = Path(str(candidate) + ".fai")
            if not fai.is_file():
                raise FileNotFoundError(
                    f"reference FASTA index {fai} not found; run "
                    f"`samtools faidx {candidate}` before starting the workflow"
                )
            selected = Assembly(
                derive_reference_name(str(candidate), by_name),
                str(candidate),
                str(fai.resolve()),
            )
            assemblies.insert(0, selected)
    return [selected] + [row for row in assemblies if row.name != selected.name]


def read_fai(path: str) -> List[Tuple[str, int]]:
    records: List[Tuple[str, int]] = []
    seen = set()
    with open(path, "rt", encoding="utf-8") as handle:
        for line_number, raw in enumerate(handle, 1):
            if not raw.strip():
                continue
            fields = raw.rstrip("\n").split("\t")
            if len(fields) < 5:
                raise ValueError(f"{path}:{line_number}: expected FAI5+")
            try:
                length = int(fields[1])
                tuple(int(value) for value in fields[2:5])
            except ValueError as error:
                raise ValueError(f"{path}:{line_number}: invalid FAI numeric field") from error
            name = fields[0]
            if not name or length <= 0:
                raise ValueError(f"{path}:{line_number}: empty contig name/sequence")
            if name in seen:
                raise ValueError(f"{path}:{line_number}: duplicate contig {name!r}")
            seen.add(name)
            records.append((name, length))
    if not records:
        raise ValueError(f"{path}: empty FASTA index")
    return records


def parse_labeled_fasta(value: str) -> Tuple[str, str]:
    """Parse SAMPLE=FASTA used by the catalog-building helper."""
    if "=" not in value:
        raise ValueError(
            f"invalid labeled FASTA {value!r}; expected SAMPLE=FASTA (use .=FASTA "
            "to preserve existing source annotations)"
        )
    sample, path = value.split("=", 1)
    sample = sample.strip()
    fasta = os.path.abspath(os.path.expanduser(path))
    if not sample or not path:
        raise ValueError(f"invalid labeled FASTA {value!r}; expected SAMPLE=FASTA")
    if ":" in sample or any(character.isspace() for character in sample):
        raise ValueError(f"catalog sample label contains an unsafe character: {sample!r}")
    if not os.path.isfile(fasta):
        raise FileNotFoundError(fasta)
    if fasta.lower().endswith((".gz", ".bgz", ".bgzf")):
        raise ValueError(f"catalog FASTA must be uncompressed: {fasta}")
    if not os.path.isfile(fasta + ".fai"):
        raise FileNotFoundError(
            f"FASTA index {fasta}.fai not found; run `samtools faidx {fasta}` first"
        )
    return sample, fasta


def combine_fastas(args: argparse.Namespace) -> None:
    """Stream labeled FASTAs into one source-annotated, duplicate-free catalog."""
    sources = [parse_labeled_fasta(value) for value in args.input]
    output = Path(args.output).resolve()
    output.parent.mkdir(parents=True, exist_ok=True)
    temporary = output.with_name(output.name + f".tmp.{os.getpid()}")
    seen = set()
    record_count = 0
    try:
        with temporary.open("wt", encoding="ascii") as destination:
            for sample, fasta in sources:
                if os.path.getsize(fasta) == 0 and os.path.getsize(fasta + ".fai") == 0:
                    continue
                lengths = dict(read_fai(fasta + ".fai"))
                observed = set()
                with open(fasta, "rt", encoding="ascii") as source:
                    for line_number, raw in enumerate(source, 1):
                        if not raw.startswith(">"):
                            destination.write(raw)
                            continue
                        fields = raw[1:].strip().split()
                        if not fields:
                            raise ValueError(f"{fasta}:{line_number}: empty FASTA header")
                        name = fields[0]
                        if name not in lengths:
                            raise ValueError(
                                f"{fasta}:{line_number}: record {name!r} is absent from its .fai"
                            )
                        if name in observed:
                            raise ValueError(f"{fasta}:{line_number}: duplicate record {name!r}")
                        if name in seen:
                            raise ValueError(
                                f"duplicate FASTA record {name!r} while combining fixed/novel catalogs"
                            )
                        observed.add(name)
                        seen.add(name)
                        record_count += 1
                        if sample == "alternative":
                            # The source sample/coordinate is provenance, not
                            # permission to refine this imported sequence.
                            fields = [field for field in fields
                                      if not field.startswith("sequence_role=")]
                            fields.append("sequence_role=imported_alternative")
                            raw = ">" + " ".join(fields) + "\n"
                        source_mapping = None
                        source_mapping_index = None
                        for field_index, field in enumerate(fields[1:], 1):
                            match = SOURCE_COORD_TOKEN_RE.fullmatch(field)
                            if match is not None:
                                source_mapping = match
                                source_mapping_index = field_index
                                break
                        if source_mapping is not None and (
                            int(source_mapping.group("end"))
                            - int(source_mapping.group("start")) != lengths[name]
                        ):
                            raise ValueError(
                                f"{fasta}:{line_number}: record {name!r} has a "
                                "source= interval whose length differs from its .fai"
                            )
                        if source_mapping is not None:
                            if not source_mapping.group("strand"):
                                assert source_mapping_index is not None
                                fields[source_mapping_index] += "+"
                                destination.write(">" + " ".join(fields) + "\n")
                            else:
                                destination.write(
                                    raw if raw.endswith("\n") else raw + "\n"
                                )
                        elif sample == ".":
                            destination.write(
                                raw if raw.endswith("\n") else raw + "\n"
                            )
                        else:
                            # Keep any non-coordinate description fields while
                            # adding an addressable source for ordinary FASTAs.
                            description = raw[1:].strip()
                            destination.write(
                                f">{description} "
                                f"source={sample}:{name}:0-{lengths[name]}+\n"
                            )
                if observed != set(lengths):
                    missing = sorted(set(lengths) - observed)
                    raise ValueError(
                        f"{fasta}: .fai contains records absent from FASTA: "
                        + ", ".join(missing[:5])
                    )
                destination.write("\n")
        if not record_count:
            raise ValueError("cannot create an empty fixed/catalog FASTA")
        os.replace(temporary, output)
    finally:
        try:
            temporary.unlink()
        except FileNotFoundError:
            pass


def augment_bed(args: argparse.Namespace) -> None:
    """Add one whole-record block for every alternative record absent from BED."""
    rows = read_bed_rows(args.input, automatic=False)
    mentioned = {row[0] for row in rows}
    used_groups = {row[3] for row in rows}
    additions: List[List[str]] = []
    for fai_path in args.alternative_fai:
        for record, length in read_fai(fai_path):
            if record in mentioned:
                continue
            digest = hashlib.sha256(record.encode("utf-8")).hexdigest()[:12]
            group = f"alternative{digest}"
            if group in used_groups:
                raise ValueError(f"generated alternative block group collision: {group}")
            used_groups.add(group)
            additions.append([record, "0", str(length), group])
    atomic_write(
        args.output,
        "".join("\t".join(row) + "\n" for row in rows + additions),
    )


def combine_beds(args: argparse.Namespace) -> None:
    chunks: List[str] = []
    records = 0
    group_owner: Dict[str, int] = {}
    for source_index, path in enumerate(args.input):
        with open(path, "rt", encoding="utf-8") as handle:
            for raw in handle:
                if not raw.strip() or raw.lstrip().startswith("#"):
                    continue
                fields = raw.split()
                if len(fields) >= 4:
                    group = fields[3]
                    owner = group_owner.setdefault(group, source_index)
                    if owner != source_index:
                        raise ValueError(
                            f"generated block group {group!r} collides between "
                            f"{args.input[owner]} and {path}"
                        )
                chunks.append(raw.rstrip("\r\n") + "\n")
                records += 1
    if not records:
        raise ValueError("combined block BED would be empty")
    atomic_write(args.output, "".join(chunks))


def validate_sample_name(name: str) -> re.Match[str]:
    match = SAMPLE_RE.fullmatch(name)
    if match is None:
        raise ValueError(
            f"invalid assembly name {name!r}; expected SAMPLE_hHAPLOTYPE, "
            "for example HG002_h1"
        )
    return match


def conforms_to_header_prefix(name: str, records: Sequence[Tuple[str, int]]) -> bool:
    match = validate_sample_name(name)
    prefix = f"{match.group('sample')}#{int(match.group('haplotype'))}"
    if not records:
        return False
    first_contig = records[0][0]
    return first_contig == prefix or first_contig.startswith(prefix + "#")


def prepare_inputs(args: argparse.Namespace) -> None:
    template_defined = bool(getattr(args, "template_defined", False))
    working_directory = getattr(args, "working_directory", None)
    if template_defined:
        assemblies = read_assembly_list(args.query_paths, working_directory)
        reference = None
    else:
        assemblies = resolve_reference_and_assemblies(
            args.query_paths, args.reference, working_directory,
        )
        reference = assemblies[0]
    catalogs: Dict[str, List[Tuple[str, int]]] = {}
    errors: List[str] = []
    for assembly in assemblies:
        try:
            validate_sample_name(assembly.name)
            if any(character.isspace() for character in assembly.fasta):
                raise ValueError(
                    "FASTA paths containing whitespace are not supported by "
                    "the downstream query_paths readers"
                )
            if assembly.fasta.lower().endswith((".gz", ".bgz", ".bgzf")):
                raise ValueError("FASTA is compressed; random-access stages require plain FASTA")
            catalogs[assembly.name] = read_fai(assembly.fai)
        except (OSError, ValueError) as error:
            errors.append(f"{assembly.name} ({assembly.fasta}): {error}")

    if errors:
        raise ValueError("\n".join(errors) + "\n" + preparation_message(args.preparation_folder))

    lines = [
        f"{row.name}\t{row.fasta}"
        for row in assemblies
    ]
    atomic_write(args.output, "\n".join(lines) + "\n")
    metadata = {
        "reference_mode": (
            "per_block_template" if template_defined else "cohort_reference"
        ),
        "reference_name": reference.name if reference is not None else None,
        "reference_fasta": reference.fasta if reference is not None else None,
        "assembly_count": len(assemblies),
        "assemblies": [
            {
                "name": row.name,
                "fasta": row.fasta,
                "fai": row.fai,
                "contigs": len(catalogs[row.name]),
                "header_prefix_conforming": conforms_to_header_prefix(
                    row.name, catalogs[row.name]
                ),
                "is_reference": (
                    reference is not None and row.name == reference.name
                ),
            }
            for row in assemblies
        ],
    }
    atomic_write(args.metadata, json.dumps(metadata, indent=2, sort_keys=True) + "\n")


def run_faidx(args: argparse.Namespace) -> None:
    expected_output = os.path.abspath(args.fasta) + ".fai"
    requested_output = os.path.abspath(args.output)
    if requested_output != expected_output:
        raise ValueError(
            f"faidx output must be adjacent to its FASTA: expected "
            f"{expected_output}, got {requested_output}"
        )
    Path(requested_output).parent.mkdir(parents=True, exist_ok=True)
    completed = subprocess.run(
        [args.samtools, "faidx", os.path.abspath(args.fasta)],
        text=True,
    )
    if completed.returncode != 0:
        raise RuntimeError(f"samtools faidx failed for {args.fasta}")
    read_fai(requested_output)


def normalize_block_name(value: str, line_number: int) -> str:
    value = value.replace("_", "").replace("-", "")
    if value in {"", "."}:
        return f"g{line_number}"
    if value.isdigit():
        return f"g{value}"
    if SAFE_BLOCK_RE.fullmatch(value) and len(value) <= MAX_BLOCK_NAME_LENGTH:
        return value

    # Block names become FASTA IDs, directory names, and filename components.
    # Keep a readable prefix while using a digest to distinguish punctuation
    # variants and to bound very long annotation-derived names.  The mapping is
    # based only on the original value, so repeated names continue to identify
    # one pre-grouped block.
    digest = hashlib.sha256(value.encode("utf-8")).hexdigest()[:12]
    slug = UNSAFE_BLOCK_RE.sub("", value).strip(".") or "block"
    prefix_length = MAX_BLOCK_NAME_LENGTH - len(digest)
    slug = slug[:prefix_length].rstrip(".") or "block"
    return f"{slug}{digest}"


def read_bed_rows(path: str, automatic: bool = False) -> List[List[str]]:
    rows: List[List[str]] = []
    with open(path, "rt", encoding="utf-8") as handle:
        for line_number, raw in enumerate(handle, 1):
            if not raw.strip() or raw.lstrip().startswith("#"):
                continue
            fields = raw.rstrip("\r\n").split("\t")
            if len(fields) < 3:
                fields = raw.split()
            if len(fields) < 3 or (not automatic and len(fields) < 4):
                expected = "BED3+" if automatic else "BED4+"
                raise ValueError(f"{path}:{line_number}: expected {expected}")
            try:
                start, end = int(fields[1]), int(fields[2])
            except ValueError as error:
                raise ValueError(f"{path}:{line_number}: invalid coordinates") from error
            if not fields[0] or start < 0 or end <= start:
                raise ValueError(f"{path}:{line_number}: invalid BED interval")
            name = f"g{len(rows) + 1}" if automatic else normalize_block_name(
                fields[3], len(rows) + 1
            )
            row = [fields[0], str(start), str(end), name]
            if not automatic and len(fields) > 4:
                row.append(fields[-1])
            rows.append(row)
    if not rows:
        raise ValueError(f"{path}: no BED intervals")
    return rows


def bed_has_groups(path: str) -> bool:
    rows = read_bed_rows(path, automatic=False)
    names = [row[3] for row in rows]
    return len(names) != len(set(names))


def normalize_bed(args: argparse.Namespace) -> None:
    rows = read_bed_rows(args.input, args.automatic)
    fai = dict(read_fai(args.reference_fai))
    for line_number, row in enumerate(rows, 1):
        contig, start_text, end_text = row[:3]
        if contig not in fai:
            raise ValueError(
                f"{args.input}:{line_number}: contig {contig!r} is absent from "
                f"the reference index {args.reference_fai}"
            )
        if int(end_text) > fai[contig]:
            raise ValueError(
                f"{args.input}:{line_number}: {contig}:{start_text}-{end_text} "
                f"exceeds reference length {fai[contig]}"
            )
    atomic_write(args.output, "".join("\t".join(row) + "\n" for row in rows))
    names = [row[3] for row in rows]
    counts = dict(sorted(Counter(names).items()))
    metadata = {
        "input": os.path.abspath(args.input),
        "records": len(rows),
        "grouped": any(count > 1 for count in counts.values()),
        "groups": counts,
        "automatic": bool(args.automatic),
    }
    atomic_write(args.metadata, json.dumps(metadata, indent=2, sort_keys=True) + "\n")


def normalize_template_bed(args: argparse.Namespace) -> None:
    """Validate grouped BED rows against their template assembly origins.

    BED column 5 may name the template haplotype. When omitted, the contig
    must occur in exactly one assembly FAI. The normalized output always has
    five columns and therefore carries an explicit template for every block.
    """
    assemblies = read_assembly_list(args.query_paths)
    by_name = {row.name: row for row in assemblies}
    contig_owners: Dict[str, List[Tuple[str, int]]] = {}
    for assembly in assemblies:
        for contig, length in read_fai(assembly.fai):
            contig_owners.setdefault(contig, []).append((assembly.name, length))

    rows = []
    group_counts: Counter[str] = Counter()
    template_counts: Counter[str] = Counter()
    with open(args.input, "rt", encoding="utf-8") as handle:
        for line_number, raw in enumerate(handle, 1):
            if not raw.strip() or raw.lstrip().startswith("#"):
                continue
            fields = raw.rstrip("\r\n").split("\t")
            if len(fields) < 4:
                fields = raw.split()
            if len(fields) < 4:
                raise ValueError(
                    f"{args.input}:{line_number}: expected BED4+"
                )
            if len(fields) > 5:
                raise ValueError(
                    f"{args.input}:{line_number}: template-defined input must "
                    "be BED4 or BED5; BED column 5 is the template assembly"
                )
            contig = fields[0]
            try:
                start, end = int(fields[1]), int(fields[2])
            except ValueError as error:
                raise ValueError(
                    f"{args.input}:{line_number}: invalid coordinates"
                ) from error
            group = normalize_block_name(fields[3], len(rows) + 1)
            explicit = fields[4] if len(fields) == 5 else ""
            if explicit and explicit not in by_name:
                raise ValueError(
                    f"{args.input}:{line_number}: unknown template assembly "
                    f"{explicit!r} in BED column 5"
                )
            owners = contig_owners.get(contig, [])
            if explicit:
                length = next(
                    (size for sample, size in owners if sample == explicit),
                    None,
                )
                if length is None:
                    raise ValueError(
                        f"{args.input}:{line_number}: template {explicit!r} "
                        f"does not contain contig {contig!r}"
                    )
                template = explicit
            elif len(owners) == 1:
                template, length = owners[0]
            elif not owners:
                raise ValueError(
                    f"{args.input}:{line_number}: contig {contig!r} is absent "
                    "from every query assembly"
                )
            else:
                raise ValueError(
                    f"{args.input}:{line_number}: contig {contig!r} occurs in "
                    f"multiple assemblies ({','.join(sample for sample, _ in owners)}); "
                    "put the template assembly name in BED column 5"
                )
            if start < 0 or end <= start or end > int(length):
                raise ValueError(
                    f"{args.input}:{line_number}: {contig}:{start}-{end} "
                    f"exceeds {template} contig length {length}"
                )
            rows.append((contig, start, end, group, template))
            group_counts[group] += 1
            template_counts[template] += 1
    if not rows:
        raise ValueError(f"{args.input}: no BED intervals")
    atomic_write(
        args.output,
        "".join("\t".join(map(str, row)) + "\n" for row in rows),
    )
    metadata = {
        "input": os.path.abspath(args.input),
        "records": len(rows),
        "grouped": any(count > 1 for count in group_counts.values()),
        "groups": dict(sorted(group_counts.items())),
        "template_defined": True,
        "templates": dict(sorted(template_counts.items())),
    }
    atomic_write(
        args.metadata, json.dumps(metadata, indent=2, sort_keys=True) + "\n",
    )


def grouped_merge(args: argparse.Namespace) -> None:
    rows = read_bed_rows(args.bed, automatic=False)
    target_records = read_fai(args.target_fai)
    known_groups = {row[3] for row in rows}
    mappings: List[Tuple[str, str]] = []
    for target, _length in target_records:
        group = target.split("_", 1)[0]
        if group not in known_groups:
            raise ValueError(
                f"target FASTA record {target!r} does not begin with a known "
                f"input BED group ({', '.join(sorted(known_groups)[:10])})"
            )
        mappings.append((target, group))
    atomic_write(
        args.output,
        "".join("\t".join(row[:4]) + "\n" for row in rows),
    )
    atomic_write(
        args.groups,
        "partition\tmerged_partition\n"
        + "".join(f"{target}\t{group}\n" for target, group in mappings),
    )


def natural_key(value: str) -> Tuple[object, ...]:
    return tuple(
        int(token) if token.isdigit() else token
        for token in re.split(r"([0-9]+)", value)
    )


def list_partitions(args: argparse.Namespace) -> None:
    root = Path(args.localblocks).resolve()
    completed_breakpoints = None
    if args.breakpoint_consistent:
        if not args.breakpoint_manifest:
            raise ValueError(
                "--breakpoint-consistent requires --breakpoint-manifest"
            )
        completed_breakpoints = set()
        with open(args.breakpoint_manifest, "rt", encoding="utf-8") as handle:
            header = handle.readline().rstrip("\r\n").split("\t")
            try:
                partition_column = header.index("partition")
                status_column = header.index("status")
            except ValueError as error:
                raise ValueError(
                    f"{args.breakpoint_manifest}: invalid breakpoint manifest header"
                ) from error
            for line_number, raw in enumerate(handle, 2):
                if not raw.strip():
                    continue
                fields = raw.rstrip("\r\n").split("\t")
                if max(partition_column, status_column) >= len(fields):
                    raise ValueError(
                        f"{args.breakpoint_manifest}:{line_number}: short row"
                    )
                if fields[status_column] in {"built", "reused"}:
                    completed_breakpoints.add(fields[partition_column])
    paths = []
    excluded = []
    for directory in root.iterdir():
        if not directory.is_dir():
            continue
        partition = directory.name
        if args.breakpoint_consistent:
            if partition not in completed_breakpoints:
                excluded.append((
                    partition, str(directory), "BREAKPOINT_STATUS",
                ))
                continue
            graph = directory / f"{partition}.FA"
            bed = directory / f"{partition}.bed"
            alignment = directory / f"{partition}_align.txt"
            initial = directory / f"{partition}_breakpoint_consistency.initial.tsv"
            blocks = directory / f"{partition}_breakpoint_consistency.bed"
            cache = directory / f"{partition}cache.json"
            refmatch = directory / f"{partition}_refmatch.tsv"
            missing = []
            for label, path, require_nonempty in (
                ("FA", graph, True),
                ("BED", bed, True),
                ("ALIGNMENT", alignment, True),
                ("INITIAL", initial, True),
                ("BLOCKS", blocks, False),
                ("CACHE", cache, True),
            ):
                if not path.is_file() or (
                    require_nonempty and path.stat().st_size == 0
                ):
                    missing.append(label)
            if not missing:
                try:
                    if cache.stat().st_mtime < max(
                        graph.stat().st_mtime,
                        bed.stat().st_mtime,
                        alignment.stat().st_mtime,
                        initial.stat().st_mtime,
                        blocks.stat().st_mtime,
                    ):
                        missing.append("CURRENT_CACHE")
                except OSError:
                    missing.append("CURRENT_CACHE")
            if getattr(args, "refmatch_consistent", False):
                if not refmatch.is_file() or refmatch.stat().st_size == 0:
                    missing.append("REFMATCH")
                elif not missing:
                    try:
                        if refmatch.stat().st_mtime < max(
                            alignment.stat().st_mtime,
                            blocks.stat().st_mtime,
                            cache.stat().st_mtime,
                        ):
                            missing.append("CURRENT_REFMATCH")
                    except OSError:
                        missing.append("CURRENT_REFMATCH")
            if missing:
                excluded.append((partition, str(directory), ",".join(missing)))
                continue
            paths.append(graph.resolve())
            continue
        suffix = "_samples.fasta" if args.samples else ".fasta"
        path = directory / f"{partition}{suffix}"
        if path.is_file() and path.stat().st_size > 0:
            paths.append(path.resolve())
    paths.sort(key=lambda value: natural_key(value.parent.name))
    if args.excluded:
        atomic_write(
            args.excluded,
            "partition\tdirectory\treason\n"
            + "".join(
                f"{partition}\t{directory}\t{reason}\n"
                for partition, directory, reason in sorted(
                    excluded, key=lambda row: natural_key(row[0]),
                )
            ),
        )
    if not paths:
        kind = (
            "breakpoint-consistent final graphs"
            if args.breakpoint_consistent else
            ("sample FASTAs" if args.samples else "partition FASTAs")
        )
        raise ValueError(f"{root}: no nonempty {kind}")
    atomic_write(args.output, "".join(str(path) + "\n" for path in paths))


def make_batches(args: argparse.Namespace) -> None:
    with open(args.input, "rt", encoding="utf-8") as handle:
        targets = [line.strip() for line in handle if line.strip()]
    if not targets:
        raise ValueError(f"{args.input}: no partition targets")
    duplicate_targets = sorted(
        target for target, copies in Counter(targets).items() if copies > 1
    )
    if duplicate_targets:
        preview = ", ".join(duplicate_targets[:5])
        raise ValueError(
            f"{args.input}: duplicate partition target(s): {preview}"
        )
    count = min(max(1, args.count), len(targets))
    output = Path(args.output)
    output.mkdir(parents=True, exist_ok=True)
    batches: List[List[str]] = [[] for _ in range(count)]
    for index, target in enumerate(targets):
        batches[index % count].append(target)
    ids = []
    for index, batch in enumerate(batches, 1):
        batch_id = f"batch{index:03d}"
        ids.append(batch_id)
        atomic_write(
            str(output / f"{batch_id}.list"),
            "".join(target + "\n" for target in batch),
        )
    atomic_write(args.ids, "".join(batch_id + "\n" for batch_id in ids))


def write_dummy(args: argparse.Namespace) -> None:
    atomic_write(args.fasta, ">index_dummy\nNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNN\n")
    atomic_write(args.list, f"index_dummy\t{Path(args.fasta).resolve()}\n")


def preflight(args: argparse.Namespace) -> None:
    missing = []
    for executable in args.executable:
        if os.path.sep in executable:
            found = os.path.isfile(executable) and os.access(executable, os.X_OK)
        else:
            found = shutil.which(executable) is not None
        if not found:
            missing.append(executable)
    for path in args.file:
        if not os.path.isfile(path):
            missing.append(path)
    if missing:
        raise FileNotFoundError(
            "required graph-pipeline tools/files are unavailable: "
            + ", ".join(missing)
            + ". "
            + preparation_message(args.preparation_folder)
        )
    atomic_write(args.output, "ok\n")


def wait_nonempty(args: argparse.Namespace) -> None:
    """Wait for completed-job outputs to become visible and nonempty.

    A Slurm compute node can close a file before its metadata is immediately
    visible on the submission node.  An instantaneous ``test -s`` therefore
    creates false Snakemake failures on shared filesystems.
    """
    timeout = float(args.timeout)
    interval = float(args.interval)
    if timeout < 0:
        raise ValueError("--timeout must be non-negative")
    if interval <= 0:
        raise ValueError("--interval must be positive")
    paths = [Path(value) for value in args.paths]
    if not paths:
        raise ValueError("wait-nonempty requires at least one path")

    deadline = time.monotonic() + timeout
    while True:
        unavailable = []
        for path in paths:
            try:
                size = path.stat().st_size
            except FileNotFoundError:
                unavailable.append(f"{path} (missing)")
            else:
                if size == 0:
                    unavailable.append(f"{path} (empty)")
        if not unavailable:
            return
        remaining = deadline - time.monotonic()
        if remaining <= 0:
            raise FileNotFoundError(
                "required output(s) did not become nonempty within "
                f"{timeout:g} seconds: " + ", ".join(unavailable)
            )
        time.sleep(min(interval, remaining))


def read_manifest(path: str) -> List[str]:
    with open(path, "rt", encoding="utf-8") as handle:
        values = [line.strip() for line in handle if line.strip()]
    if not values:
        raise ValueError(f"{path}: empty batch")
    return values


def resolve_partition_header_manifest_entry(path: str) -> Path:
    """Resolve current and stale batch entries to a partition header file."""
    entry = Path(path).resolve()
    partition = entry.parent.name
    canonical = entry.parent / f"{partition}_samples.header"
    legacy = entry.parent / f"{partition}.header"
    recognized = {
        canonical.name,
        legacy.name,
        f"{partition}_samples.fasta",
    }
    if entry.name not in recognized:
        raise ValueError(f"invalid partition header manifest: {entry}")
    for candidate in (canonical, legacy):
        if candidate.is_file() and candidate.stat().st_size > 0:
            return candidate
    raise ValueError(
        f"invalid partition header manifest: {entry}; expected a nonempty "
        f"{canonical}"
    )


def run_slurm(args: argparse.Namespace) -> None:
    if args.cpus < 1:
        raise ValueError("--cpus must be positive")
    command = list(args.slurm_command)
    if command and command[0] == "--":
        command.pop(0)
    if not command:
        raise ValueError("run-slurm requires a command after --")
    job_name = re.sub(r"[^A-Za-z0-9_.-]+", "-", args.job_name).strip(".-")
    if not job_name:
        raise ValueError("--job-name contains no safe characters")
    log_dir = Path(args.log_dir).resolve()
    log_dir.mkdir(parents=True, exist_ok=True)
    sbatch_args = list(getattr(args, "sbatch_arg", []) or [])
    if any(
        value == "--" or value == "--wrap" or value.startswith("--wrap=")
        for value in sbatch_args
    ):
        raise ValueError("--sbatch-arg cannot replace the managed --wrap command")

    def has_option(long_name: str, short_name: str = "") -> bool:
        for value in sbatch_args:
            if value == long_name or value.startswith(long_name + "="):
                return True
            if short_name and (value == short_name or (
                value.startswith(short_name) and not value.startswith("--")
            )):
                return True
        return False

    sbatch_command = [
        args.sbatch,
        "--parsable",
        "--wait",
        "--export=ALL",
        f"--job-name={job_name}",
    ]
    if not has_option("--cpus-per-task", "-c"):
        sbatch_command.append(f"--cpus-per-task={args.cpus}")
    if not has_option("--mem") and not has_option("--mem-per-cpu"):
        sbatch_command.append(f"--mem={args.memory}")
    if args.account and not has_option("--account", "-A"):
        sbatch_command.append(f"--account={args.account}")
    if args.time and not has_option("--time", "-t"):
        sbatch_command.append(f"--time={args.time}")
    if args.partition and not has_option("--partition", "-p"):
        sbatch_command.append(f"--partition={args.partition}")
    sbatch_command.extend([
        f"--output={log_dir / (job_name + '-%j.out')}",
        f"--error={log_dir / (job_name + '-%j.err')}",
    ])
    sbatch_command.extend(sbatch_args)
    sbatch_command.extend(["--wrap", shlex.join(command)])
    print(f"[SLURM] {shlex.join(sbatch_command)}", flush=True)
    completed = subprocess.run(sbatch_command, check=False)
    if completed.returncode != 0:
        raise RuntimeError(
            f"Slurm job {job_name!r} failed with status {completed.returncode}"
        )


def _run_locked_batch_partition(
    args: argparse.Namespace, header: Path,
) -> List[str]:
    """Process one partition while its interprocess lock is held."""
    partition = header.parent.name
    unique_bed = header.parent / f"{partition}.unique.bed"
    completed_workdir: Path
    completion_marker: Optional[Path] = None
    print(f"[BATCH] Starting partition: {partition}", flush=True)
    if args.mode == "integrated":
        return _run_integrated_partition(args, header)
    if args.mode == "minset":
        completed_workdir = Path(str(unique_bed) + ".work")
        completion_marker, completion_protocol = minset_completion_marker(
            unique_bed, args,
        )
        current = (
            unique_bed.is_file()
            and unique_bed.stat().st_size > 0
            and unique_bed.stat().st_mtime_ns >= header.stat().st_mtime_ns
            and completion_marker.is_file()
            and completion_marker.read_text(encoding="utf-8")
            == completion_protocol
            and completion_marker.stat().st_mtime_ns
            >= unique_bed.stat().st_mtime_ns
        )
        if not current:
            scratch_parent = os.environ.get("SLURM_TMPDIR") or None
            with tempfile.TemporaryDirectory(
                prefix=f"minsetref_{partition}_", dir=scratch_parent,
            ) as scratch:
                samples = Path(scratch) / f"{partition}_samples.fasta"
                subprocess.run([
                    args.sample_fasta_builder,
                    "--materialize-header", str(header),
                    "-q", args.query_paths,
                    "--output-fasta", str(samples),
                    "-j", str(args.threads),
                ], check=True)
                command = [
                    args.python, args.script,
                    "-i", str(samples),
                    "-o", str(unique_bed),
                    "-r", args.reference_name,
                    "-p", args.reference_name,
                    "-c", str(args.threads),
                    "--workdir", str(completed_workdir),
                    "--resume",
                ]
                if not getattr(args, "max_rigorous", False):
                    command.extend((
                        "--kmer-preclean",
                        "--kmer-ref-cleaner", args.kmer_ref_cleaner,
                    ))
                if args.skip_blastn:
                    command.append("--skip-blastn")
                subprocess.run(command, check=True)
        outputs = [unique_bed]
    else:
        if not unique_bed.is_file() or unique_bed.stat().st_size == 0:
            raise FileNotFoundError(unique_bed)
        adjusted = (
            Path(args.adjusted_dir).resolve()
            / partition / f"{partition}.adjusted.tsv"
        )
        if not adjusted.is_file() or adjusted.stat().st_size == 0:
            raise FileNotFoundError(adjusted)
        partition_dir = Path(args.output_dir).resolve() / partition
        partition_dir.mkdir(parents=True, exist_ok=True)
        completed_workdir = partition_dir / f"{partition}.work"
        graph = partition_dir / f"{partition}.FA"
        align = partition_dir / f"{partition}_align.txt"
        linear = partition_dir / f"{partition}_linear.txt"
        command = [
            args.python, args.script,
            "-i", str(adjusted),
            "--input-format", "adjusted",
            "--hotspot-header", str(header),
            "--sample-fasta-builder", args.sample_fasta_builder,
            "-q", args.query_paths,
            "-o", str(graph),
            "--align", str(align),
            "-r", args.reference_name,
            "-p", args.reference_name,
            "-t", str(args.threads),
            "--timeout", str(getattr(
                args, "alignment_timeout", DEFAULT_ALIGNMENT_TIMEOUT,
            )),
            "--workdir", str(completed_workdir),
            "--resume",
        ]
        if getattr(args, "fast_mode", False):
            command.append("--fast-mode")
        elif (
            getattr(args, "slow_rigorous", False)
            or getattr(args, "max_rigorous", False)
        ):
            command.append("--slow-rigorous")
        outputs = [graph, align]
        if getattr(args, "linear", True):
            command.extend(["--linear", str(linear)])
            outputs.append(linear)
        if args.skip_blastn:
            command.append("--skip-blastn")
        command.append(
            "--alignment-payload"
            if args.alignment_payload
            else "--no-alignment-payload"
        )
        if not run_partition_graph_command(
            args, partition, partition_dir, command,
        ):
            return []
    for output in outputs:
        if not output.is_file() or output.stat().st_size == 0:
            raise RuntimeError(
                f"{partition}: expected nonempty output missing: {output}"
            )
    if completion_marker is not None:
        atomic_write(str(completion_marker), completion_protocol)
    if completed_workdir.exists():
        shutil.rmtree(completed_workdir)
        print(
            f"[CLEANUP] Removed completed partition work directory: "
            f"{completed_workdir}",
            flush=True,
        )
    print(f"[BATCH] Finished partition: {partition}", flush=True)
    return [
        f"{partition}\t{output.resolve()}\t{output.stat().st_size}"
        for output in outputs
    ]


def _copy_partition_adjusted_outputs(source: Path, destination: Path) -> None:
    """Atomically install one partition's compact local-refinement outputs."""
    destination.mkdir(parents=True, exist_ok=True)
    for source_path in source.iterdir():
        if not source_path.is_file():
            continue
        temporary = destination / (
            source_path.name + f".tmp.{os.getpid()}.{threading.get_ident()}"
        )
        try:
            shutil.copy2(source_path, temporary)
            os.replace(temporary, destination / source_path.name)
        finally:
            try:
                temporary.unlink()
            except FileNotFoundError:
                pass


def partition_template_haplotypes(header: Path) -> List[str]:
    """Read authoritative template haplotypes from PARTITION.bed."""
    partition = header.parent.name
    bed = header.parent / f"{partition}.bed"
    output: List[str] = []
    seen = set()
    with bed.open("rt", encoding="utf-8") as handle:
        for line_number, raw in enumerate(handle, 1):
            if not raw.strip() or raw.startswith("#"):
                continue
            fields = raw.rstrip("\r\n").split("\t")
            if len(fields) < 7:
                fields = raw.split()
            if len(fields) < 7 or not fields[6]:
                raise ValueError(
                    f"{bed}:{line_number}: expected template haplotype in "
                    "BED column 7"
                )
            if fields[6] not in seen:
                output.append(fields[6])
                seen.add(fields[6])
    if not output:
        raise ValueError(f"{bed}: no template haplotypes")
    return output


def selected_minset_reference_haplotype(path: Path) -> str:
    """Return the single haplotype labeled reference by minsetref_light."""
    selected = set()
    with path.open("rt", encoding="utf-8") as handle:
        for line_number, raw in enumerate(handle, 1):
            if not raw.strip() or raw.startswith("#"):
                continue
            fields = raw.rstrip("\r\n").split("\t")
            if len(fields) < 9:
                fields = raw.split()
            if len(fields) < 9:
                raise ValueError(
                    f"{path}:{line_number}: expected minset BED9+"
                )
            if fields[7] == "reference":
                selected.add(fields[6])
    if len(selected) != 1:
        raise ValueError(
            f"{path}: expected exactly one selected template haplotype, "
            f"found {sorted(selected)}"
        )
    return next(iter(selected))


def _run_integrated_partition(
    args: argparse.Namespace, header: Path,
) -> List[str]:
    """Materialize once, then clean, refine, graph, align, and linearize."""
    partition = header.parent.name
    unique_bed = header.parent / f"{partition}.unique.bed"
    unique_work = Path(str(unique_bed) + ".work")
    cleaning_marker, cleaning_protocol = minset_completion_marker(
        unique_bed, args,
    )
    graph_dir = Path(args.output_dir).resolve() / partition
    graph_dir.mkdir(parents=True, exist_ok=True)
    graph = graph_dir / f"{partition}.FA"
    align = graph_dir / f"{partition}_align.txt"
    linear = graph_dir / f"{partition}_linear.txt"
    graph_work = graph_dir / f"{partition}.work"
    lifecycle_marker = graph_dir / f"{partition}.integrated-v1.complete"
    persistent_adjusted = (
        Path(args.adjusted_dir).resolve()
        / partition / f"{partition}.adjusted.tsv"
    )
    static_block = bool(getattr(args, "static_block", False))
    linear_enabled = bool(getattr(args, "linear", True))
    protocol = integrated_partition_protocol(args)
    sys.path.insert(0, str(Path(__file__).resolve().parents[3]))
    from fixed_alternatives import has_fixed_templates, validate_fixed_paths
    original_templates = header.parent / f'{partition}.fasta'
    fixed = has_fixed_templates(original_templates)
    # An updated policy must not reuse graphs built from expanded alternatives.
    if fixed:
        protocol += 'fixed-imported-alternatives-v1\n'
    template_defined = bool(getattr(args, "template_defined_reference", False))
    template_haplotypes = (
        partition_template_haplotypes(header) if template_defined else []
    )
    final_outputs = ([] if static_block else [unique_bed]) + [persistent_adjusted, graph, align]
    if linear_enabled:
        final_outputs.append(linear)
    lifecycle_current = (
        lifecycle_marker.is_file()
        and lifecycle_marker.read_text(encoding="utf-8")
        == protocol
        and all(path.is_file() and path.stat().st_size > 0 for path in final_outputs)
        and lifecycle_marker.stat().st_mtime_ns
        >= max(path.stat().st_mtime_ns for path in [header, header.parent / f'{partition}.fasta', *final_outputs])
        and (static_block or (
            cleaning_marker.is_file()
            and cleaning_marker.read_text(encoding="utf-8") == cleaning_protocol
        ))
    )
    if lifecycle_current and static_block:
        from static_blocks import validate_static_paths
        try:
            validate_static_paths(original_templates, persistent_adjusted, graph)
        except (ValueError, KeyError, OSError) as error:
            lifecycle_current = False
            print(f'[LIFECYCLE] {partition}: invalid static templates; rebuilding: {error}', flush=True)
    if lifecycle_current and fixed:
        try:
            validate_fixed_paths(original_templates, persistent_adjusted, graph)
        except (ValueError, KeyError, OSError) as error:
            lifecycle_current = False
            print(f'[LIFECYCLE] {partition}: invalid imported template; rebuilding: {error}', flush=True)
            # An invalid completed graph must not seed --resume caches.
            if graph_work.exists():
                shutil.rmtree(graph_work)
    if lifecycle_current:
        update_partition_failure(
            partition_failure_list(args), partition, None,
        )
        for workdir in (unique_work, graph_work):
            if workdir.exists():
                shutil.rmtree(workdir)
        print(
            f"[BATCH] Reusing completed integrated partition: {partition}",
            flush=True,
        )
        return [
            f"{partition}\t{output.resolve()}\t{output.stat().st_size}"
            for output in final_outputs
        ]

    if graph_work.exists() and (
        static_block or (lifecycle_marker.is_file()
                         and lifecycle_marker.read_text().startswith("static-block-v1\n"))
    ):
        shutil.rmtree(graph_work)
    scratch_parent = os.environ.get("SLURM_TMPDIR") or None
    with tempfile.TemporaryDirectory(
        prefix=f"minsetref_integrated_{partition}_", dir=scratch_parent,
    ) as scratch_text:
        scratch = Path(scratch_text)
        samples = scratch / f"{partition}_samples.fasta"
        subprocess.run([
            args.sample_fasta_builder,
            "--materialize-header", str(header),
            "-q", args.query_paths,
            "--output-fasta", str(samples),
            "-j", str(args.threads),
        ], check=True)
        if not samples.is_file() or samples.stat().st_size == 0:
            raise RuntimeError(
                f"{partition}: materialized sample FASTA is missing or empty"
            )
        print(
            f"[LIFECYCLE] {partition}: materialized {samples}", flush=True,
        )

        if static_block:
            from static_blocks import write_static_paths
            adjusted = scratch / "adjusted"
            write_static_paths(original_templates, adjusted)
            partition_reference = (
                template_haplotypes[0] if template_defined else args.reference_name
            )
        else:
            cleaning_current = (
                unique_bed.is_file()
                and unique_bed.stat().st_size > 0
                and unique_bed.stat().st_mtime_ns >= header.stat().st_mtime_ns
                and cleaning_marker.is_file()
                and cleaning_marker.read_text(encoding="utf-8")
                == cleaning_protocol
                and cleaning_marker.stat().st_mtime_ns
                >= unique_bed.stat().st_mtime_ns
            )
            if not cleaning_current:
                command = [
                    args.python, args.script,
                    "-i", str(samples),
                    "-o", str(unique_bed),
                    "-c", str(args.threads),
                    "--workdir", str(unique_work),
                    "--resume",
                ]
                if not getattr(args, "max_rigorous", False):
                    command.extend((
                        "--kmer-preclean",
                        "--kmer-ref-cleaner", args.kmer_ref_cleaner,
                    ))
                if template_defined:
                    for haplotype in template_haplotypes:
                        command.extend(("-r", haplotype))
                    command.extend(("-p", ",".join(template_haplotypes)))
                else:
                    command.extend((
                        "-r", args.reference_name,
                        "-p", args.reference_name,
                    ))
                if args.skip_blastn:
                    command.append("--skip-blastn")
                subprocess.run(command, check=True)
                if not unique_bed.is_file() or unique_bed.stat().st_size == 0:
                    raise RuntimeError(
                        f"{partition}: minset cleaning produced no unique BED"
                    )
                atomic_write(str(cleaning_marker), cleaning_protocol)
            print(
                f"[LIFECYCLE] {partition}: unique BED ready: {unique_bed}",
                flush=True,
            )

            if template_defined:
                partition_reference = selected_minset_reference_haplotype(
                    unique_bed,
                )
                if partition_reference not in template_haplotypes:
                    raise ValueError(
                        f"{partition}: selected reference {partition_reference!r} "
                        f"is not one of its defined templates {template_haplotypes}"
                    )
                assembly_by_name = {
                    row.name: row for row in read_assembly_list(args.query_paths)
                }
                try:
                    partition_main_reference = assembly_by_name[
                        partition_reference
                    ].fasta
                except KeyError as error:
                    raise ValueError(
                        f"{partition}: template reference {partition_reference!r} "
                        "is absent from query paths"
                    ) from error
            else:
                partition_reference = args.reference_name
                partition_main_reference = args.main_reference

            localblocks = scratch / "localblocks"
            localblocks.mkdir()
            os.symlink(header.parent, localblocks / partition, target_is_directory=True)
            analysis = scratch / "coordinate_analysis"
            refinement = scratch / "local_refinement"
            adjusted = scratch / "adjusted"
            coordinate_command = [
                args.python, args.coordinate_script,
                "-i", str(localblocks),
                "-o", str(analysis),
                "-j", "1",
                "--sweep-jobs", "1",
                "--input-anchor", "50",
                "--reference-haplotypes", partition_reference,
                "--allow-missing-reference-haplotypes",
                "--ignore-contigs", "chrM",
            ]
            if getattr(args, "allow_missing_bed_haplotypes", False):
                coordinate_command.append("--allow-missing-bed-haplotypes")
            subprocess.run(coordinate_command, check=True)
            subprocess.run([
                args.python, args.refinement_script,
                "-a", str(analysis),
                "-q", args.query_paths,
                "-r", partition_main_reference,
                "-o", str(refinement),
                "-j", "1",
                "-c", str(args.threads),
                "--minimum-flank", "5000",
                "--maximum-template-gap", "100",
                "--minimum-unmasked", "100",
            ], check=True)
            subprocess.run([
                args.python, args.apply_script,
                "-a", str(analysis),
                "-u", str(refinement),
                "-q", args.query_paths,
                "-r", partition_main_reference,
                "-o", str(adjusted),
                "--min-score", "100",
                "-j", str(args.threads),
            ], check=True)
        adjusted_partition = adjusted / partition
        adjusted_manifest = adjusted_partition / f"{partition}.adjusted.tsv"
        if not adjusted_manifest.is_file() or adjusted_manifest.stat().st_size == 0:
            raise RuntimeError(
                f"{partition}: local refinement produced no adjusted manifest"
            )
        if fixed:
            validate_fixed_paths(original_templates, adjusted_manifest)
        print(
            f"[LIFECYCLE] {partition}: " + (
                "static block templates ready" if static_block else "local refinement complete"
            ), flush=True,
        )

        graph_command = [
            args.python, args.graph_script,
            "-i", str(adjusted_manifest),
            "--input-format", "adjusted",
            "--hotspot-header", str(header),
            "--materialized-samples", str(samples),
            "-q", args.query_paths,
            "-o", str(graph),
            "--align", str(align),
            "-r", partition_reference,
            "-p", partition_reference,
            "-t", str(args.threads),
            "--timeout", str(getattr(
                args, "alignment_timeout", DEFAULT_ALIGNMENT_TIMEOUT,
            )),
            "--workdir", str(graph_work),
            "--resume",
            (
                "--alignment-payload"
                if args.alignment_payload else "--no-alignment-payload"
            ),
        ]
        if getattr(args, "fast_mode", False):
            graph_command.append("--fast-mode")
        elif (
            getattr(args, "slow_rigorous", False)
            or getattr(args, "max_rigorous", False)
        ):
            graph_command.append("--slow-rigorous")
        graph_outputs = [graph, align]
        if linear_enabled:
            graph_command.extend(["--linear", str(linear)])
            graph_outputs.append(linear)
        if args.skip_blastn:
            graph_command.append("--skip-blastn")
        if not run_partition_graph_command(
            args, partition, graph_dir, graph_command,
        ):
            return []
        for output in graph_outputs:
            if not output.is_file() or output.stat().st_size == 0:
                raise RuntimeError(
                    f"{partition}: expected graph output missing: {output}"
                )
        if fixed:
            validate_fixed_paths(original_templates, adjusted_manifest, graph)
        if static_block:
            from static_blocks import validate_static_paths
            validate_static_paths(original_templates, adjusted_manifest, graph)
        _copy_partition_adjusted_outputs(
            adjusted_partition, persistent_adjusted.parent,
        )
        if not persistent_adjusted.is_file() or persistent_adjusted.stat().st_size == 0:
            raise RuntimeError(
                f"{partition}: persistent adjusted manifest is missing"
            )
        atomic_write(str(lifecycle_marker), protocol)
        print(
            f"[LIFECYCLE] {partition}: graph and alignment ready"
            + ("; linear output ready" if linear_enabled else ""),
            flush=True,
        )

    for workdir in (unique_work, graph_work):
        if workdir.exists():
            shutil.rmtree(workdir)
            print(
                f"[CLEANUP] Removed completed cycle work directory: {workdir}",
                flush=True,
            )
    print(
        f"[CLEANUP] Deleted materialized sample FASTA after all outputs: "
        f"{partition}",
        flush=True,
    )
    return [
        f"{partition}\t{output.resolve()}\t{output.stat().st_size}"
        for output in final_outputs
    ]


def run_batch_partition(args: argparse.Namespace, header_path: str) -> List[str]:
    """Serialize the same partition across concurrent batch processes."""
    header = resolve_partition_header_manifest_entry(header_path)
    partition = header.parent.name
    lock_path = header.parent / f".{partition}.pipeline.lock"
    with lock_path.open("a+") as lock_handle:
        print(f"[LOCK] Waiting for partition lock: {partition}", flush=True)
        fcntl.flock(lock_handle.fileno(), fcntl.LOCK_EX)
        try:
            print(f"[LOCK] Acquired partition lock: {partition}", flush=True)
            return _run_locked_batch_partition(args, header)
        finally:
            fcntl.flock(lock_handle.fileno(), fcntl.LOCK_UN)


def run_batch(args: argparse.Namespace) -> None:
    if args.threads < 1:
        raise ValueError("--threads must be positive")
    if args.jobs < 1:
        raise ValueError("--jobs must be positive")
    if getattr(args, "alignment_timeout", DEFAULT_ALIGNMENT_TIMEOUT) < 1:
        raise ValueError("--alignment-timeout must be positive")
    if not args.query_paths:
        raise ValueError("run-batch requires --query-paths")
    if not args.sample_fasta_builder:
        raise ValueError("run-batch requires --sample-fasta-builder")
    if (
        not getattr(args, "template_defined_reference", False)
        and hasattr(args, "reference_name")
        and not args.reference_name
    ):
        raise ValueError("run-batch requires --reference-name")
    static_block = bool(getattr(args, "static_block", False))
    if static_block and args.mode != "integrated":
        raise ValueError("--static-block requires integrated mode")
    if args.mode in {"minset", "integrated"} and not static_block and not args.kmer_ref_cleaner:
        raise ValueError(
            f"{args.mode} run-batch requires --kmer-ref-cleaner"
        )
    if args.mode == "graph" and not args.adjusted_dir:
        raise ValueError("graph run-batch requires --adjusted-dir")
    if args.mode == "integrated":
        required = {
            "--graph-script": args.graph_script,
            "--coordinate-script": args.coordinate_script,
            "--refinement-script": args.refinement_script,
            "--apply-script": args.apply_script,
            "--adjusted-dir": args.adjusted_dir,
        }
        if static_block:
            for option in ("--coordinate-script", "--refinement-script", "--apply-script"):
                required.pop(option)
        elif not getattr(args, "template_defined_reference", False):
            required["--main-reference"] = args.main_reference
        missing = [option for option, value in required.items() if not value]
        if missing:
            raise ValueError(
                "integrated run-batch requires " + ", ".join(missing)
            )
    targets = read_manifest(args.manifest)
    canonical_targets = [
        str(resolve_partition_header_manifest_entry(target)) for target in targets
    ]
    duplicate_targets = sorted(
        target
        for target, copies in Counter(canonical_targets).items()
        if copies > 1
    )
    if duplicate_targets:
        preview = ", ".join(duplicate_targets[:5])
        raise ValueError(
            f"{args.manifest}: duplicate partition assignment(s): {preview}"
        )
    targets = canonical_targets
    Path(args.output_dir).mkdir(parents=True, exist_ok=True)
    completed_by_index: List[Optional[List[str]]] = [None] * len(targets)
    if args.jobs == 1:
        for index, header_path in enumerate(targets):
            completed_by_index[index] = run_batch_partition(args, header_path)
    else:
        with ThreadPoolExecutor(max_workers=args.jobs) as executor:
            future_indexes = {
                executor.submit(run_batch_partition, args, header_path): index
                for index, header_path in enumerate(targets)
            }
            try:
                for future in as_completed(future_indexes):
                    completed_by_index[future_indexes[future]] = future.result()
            except BaseException:
                for future in future_indexes:
                    future.cancel()
                raise
    completed = [
        row
        for rows in completed_by_index
        if rows is not None
        for row in rows
    ]
    atomic_write(args.done, "partition\toutput\tbytes\n" + "\n".join(completed) + "\n")


def parse_args(argv: Optional[Sequence[str]] = None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    subparsers = parser.add_subparsers(dest="command", required=True)

    faidx = subparsers.add_parser("faidx")
    faidx.add_argument("--samtools", default="samtools")
    faidx.add_argument("--fasta", required=True)
    faidx.add_argument("--output", required=True)

    prepare = subparsers.add_parser("prepare-inputs")
    prepare.add_argument("--working-directory", help="launch directory for relative assembly paths")
    prepare.add_argument("--query-paths", required=True)
    prepare.add_argument("--reference", required=True)
    prepare.add_argument("--output", required=True)
    prepare.add_argument("--metadata", required=True)
    prepare.add_argument("--preparation-folder", required=True)
    prepare.add_argument(
        "--template-defined", action="store_true",
        help="preserve input order and record references as per-block templates",
    )

    combine_fasta = subparsers.add_parser("combine-fastas")
    combine_fasta.add_argument(
        "--input", action="append", required=True, metavar="SAMPLE=FASTA",
    )
    combine_fasta.add_argument("--output", required=True)

    augment = subparsers.add_parser("augment-bed")
    augment.add_argument("--input", required=True)
    augment.add_argument("--alternative-fai", action="append", default=[])
    augment.add_argument("--output", required=True)

    combine_bed = subparsers.add_parser("combine-beds")
    combine_bed.add_argument("--input", action="append", required=True)
    combine_bed.add_argument("--output", required=True)

    normalize = subparsers.add_parser("normalize-bed")
    normalize.add_argument("--input", required=True)
    normalize.add_argument("--reference-fai", required=True)
    normalize.add_argument("--output", required=True)
    normalize.add_argument("--metadata", required=True)
    normalize.add_argument("--automatic", action="store_true")

    normalize_template = subparsers.add_parser("normalize-template-bed")
    normalize_template.add_argument("--input", required=True)
    normalize_template.add_argument("--query-paths", required=True)
    normalize_template.add_argument("--output", required=True)
    normalize_template.add_argument("--metadata", required=True)

    merge = subparsers.add_parser("grouped-merge")
    merge.add_argument("--bed", required=True)
    merge.add_argument("--target-fai", required=True)
    merge.add_argument("--output", required=True)
    merge.add_argument("--groups", required=True)

    listing = subparsers.add_parser("list-partitions")
    listing.add_argument("--localblocks", required=True)
    listing.add_argument("--output", required=True)
    listing.add_argument("--samples", action="store_true")
    listing.add_argument(
        "--breakpoint-consistent", action="store_true",
        help=(
            "write final .FA paths only when alignment, breakpoint blocks, "
            "and the all-sample cache are complete/current"
        ),
    )
    listing.add_argument(
        "--excluded",
        help="optional TSV report of partitions omitted from the final list",
    )
    listing.add_argument(
        "--breakpoint-manifest",
        help="breakpoint-consistency status TSV used for final inclusion",
    )
    listing.add_argument(
        "--refmatch-consistent", action="store_true",
        help=(
            "also require a current nonempty PARTITION_refmatch.tsv written "
            "after the cache and blocking output"
        ),
    )

    batches = subparsers.add_parser("make-batches")
    batches.add_argument("--input", required=True)
    batches.add_argument("--output", required=True)
    batches.add_argument("--ids", required=True)
    batches.add_argument("--count", type=int, default=100)

    dummy = subparsers.add_parser("write-dummy")
    dummy.add_argument("--fasta", required=True)
    dummy.add_argument("--list", required=True)

    check = subparsers.add_parser("preflight")
    check.add_argument("--executable", action="append", default=[])
    check.add_argument("--file", action="append", default=[])
    check.add_argument("--output", required=True)
    check.add_argument("--preparation-folder", required=True)

    wait = subparsers.add_parser("wait-nonempty")
    wait.add_argument("--timeout", type=float, default=60.0)
    wait.add_argument("--interval", type=float, default=1.0)
    wait.add_argument("paths", nargs="+")

    slurm = subparsers.add_parser("run-slurm")
    slurm.add_argument("--sbatch", default="sbatch")
    slurm.add_argument("--cpus", type=int, required=True)
    slurm.add_argument("--memory", required=True)
    slurm.add_argument("--account", default="")
    slurm.add_argument("--time", default="")
    slurm.add_argument("--partition", default="")
    slurm.add_argument(
        "--sbatch-arg", action="append", default=[],
        help="one extra sbatch argument; repeat as needed",
    )
    slurm.add_argument("--job-name", required=True)
    slurm.add_argument("--log-dir", required=True)
    # Do not call this positional ``command``: the subparser destination uses
    # that name to select the action in main().  Reusing it replaces
    # ``run-slurm`` with a list and makes actions[args.command] fail.
    slurm.add_argument("slurm_command", nargs=argparse.REMAINDER)

    batch = subparsers.add_parser("run-batch")
    batch.add_argument(
        "--mode", choices=("minset", "graph", "integrated"), required=True,
    )
    batch.add_argument("--manifest", required=True)
    batch.add_argument("--script", required=True)
    batch.add_argument("--graph-script", default="")
    batch.add_argument("--coordinate-script", default="")
    batch.add_argument("--refinement-script", default="")
    batch.add_argument("--apply-script", default="")
    batch.add_argument("--python", default=sys.executable)
    batch.add_argument("--reference-name", default="")
    batch.add_argument(
        "--template-defined-reference", action="store_true",
        help=(
            "select each partition reference from its local block templates "
            "instead of imposing --reference-name cohort-wide"
        ),
    )
    batch.add_argument(
        "--allow-missing-bed-haplotypes", action="store_true",
        help=(
            "allow PARTITION.bed templates from an external alternative "
            "FASTA even when their haplotypes are absent from query paths"
        ),
    )
    batch.add_argument("--main-reference", default="")
    batch.add_argument("--query-paths", default="")
    batch.add_argument("--sample-fasta-builder", default="")
    batch.add_argument("--kmer-ref-cleaner", default="")
    batch.add_argument("--adjusted-dir", default="")
    batch.add_argument("--output-dir", required=True)
    batch.add_argument("--threads", type=int, required=True)
    batch.add_argument(
        "--alignment-timeout", type=int, default=DEFAULT_ALIGNMENT_TIMEOUT,
        help="seconds allowed for each graph alignment command (default: 21600)",
    )
    batch.add_argument(
        "--fail-list", default="",
        help="timed-out partition report (default: OUTPUT_DIR/fail.list)",
    )
    batch.add_argument(
        "--jobs", type=int, default=1,
        help="partitions processed concurrently within this batch (default: 1)",
    )
    batch.add_argument("--done", required=True)
    batch.add_argument("--skip-blastn", action="store_true")
    batch.add_argument("--static-block", action="store_true")
    alignment_modes = batch.add_mutually_exclusive_group()
    alignment_modes.add_argument(
        "--fast", "--fast-mode", dest="fast_mode", action="store_true",
    )
    alignment_modes.add_argument("--slow-rigorous", action="store_true")
    alignment_modes.add_argument("--max-rigorous", action="store_true")
    linear_group = batch.add_mutually_exclusive_group()
    linear_group.add_argument("--linear", dest="linear", action="store_true")
    linear_group.add_argument("--no-linear", dest="linear", action="store_false")
    batch.set_defaults(linear=None)
    payload_group = batch.add_mutually_exclusive_group()
    payload_group.add_argument(
        "--alignment-payload", dest="alignment_payload", action="store_true",
        default=True,
    )
    payload_group.add_argument(
        "--no-alignment-payload", dest="alignment_payload", action="store_false",
    )
    args = parser.parse_args(argv)
    if args.command == "run-batch" and args.linear is None:
        args.linear = not args.fast_mode
    return args


def main(argv: Optional[Sequence[str]] = None) -> int:
    args = parse_args(argv)
    actions = {
        "faidx": run_faidx,
        "prepare-inputs": prepare_inputs,
        "combine-fastas": combine_fastas,
        "augment-bed": augment_bed,
        "combine-beds": combine_beds,
        "normalize-bed": normalize_bed,
        "normalize-template-bed": normalize_template_bed,
        "grouped-merge": grouped_merge,
        "list-partitions": list_partitions,
        "make-batches": make_batches,
        "write-dummy": write_dummy,
        "preflight": preflight,
        "wait-nonempty": wait_nonempty,
        "run-slurm": run_slurm,
        "run-batch": run_batch,
    }
    actions[args.command](args)
    return 0


if __name__ == "__main__":
    try:
        raise SystemExit(main())
    except (OSError, ValueError, RuntimeError, subprocess.CalledProcessError) as error:
        print(f"graph pipeline: error: {error}", file=sys.stderr)
        raise SystemExit(1)
