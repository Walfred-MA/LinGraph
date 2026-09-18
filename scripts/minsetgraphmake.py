#!/usr/bin/env python3
"""Build, align, block, and linearize one minimum-set local graph.

Accepted inputs are a minsetref-style multi-FASTA, BED regions plus an indexed
``-q/--query-paths`` table, or (with ``--resume``) an existing graph ``.FA``.
The normal build path is::

    BED -> _samples.fasta -> KmerRefCleaner -> minsetref_light -> graph .FA

or the corresponding FASTA path without the first extraction step.  Optional
stages are enabled by output-path arguments: ``--align``, ``--fixbreaks``, and
``--linear``.
"""
from __future__ import annotations

import argparse
import logging
import os
import re
import shutil
import subprocess
import sys
import tempfile
from pathlib import Path
from typing import Dict, Iterable, List, Optional, Sequence, Tuple

from align_partition_hotspots import (
    DirectSequence,
    read_direct_bed,
)
from block_partition_alignments import block_one_alignment
from build_local_graphs import (
    build_one as build_one_adjusted_partition,
    public_path_id,
    read_manifest as read_adjusted_manifest,
)
from build_folder_uniformbreaks import uniformbreaks_from_blocks
from localgraphtolinear import graphtolinear as local_graph_to_linear
from minsetref_core import IndexedFasta, wrap_fasta
from partition_hotspot_segments import build_partition_segments
from summarize_partition_hotspot_segments import (
    read_alignment_output,
    select_alignment_rows,
    write_final_segments,
)
from uniform_graph_blocks import load_graphfixbreaks, read_initial_blocks


LOG = logging.getLogger("minsetgraphmake")
HERE = Path(__file__).resolve().parent
SOURCE_COORDINATE_RE = re.compile(r"^(.+):([^:]+):(\d+)-(\d+)([+-])$")
UNIQUE_REGION_RE = re.compile(r"^unique_region_\d+_(.+)_(\d+)_(\d+)$")
MINSET_BED_CLASSES = {"reference", "unique_query"}


def output_is_current(output: str, inputs: Iterable[Optional[str]]) -> bool:
    if not os.path.isfile(output) or os.path.getsize(output) == 0:
        return False
    output_time = os.stat(output).st_mtime_ns
    existing = [path for path in inputs if path and os.path.isfile(path)]
    return bool(existing) and output_time >= max(
        os.stat(path).st_mtime_ns for path in existing
    )


def atomic_copy(source: str, destination: str) -> None:
    source = os.path.abspath(source)
    destination = os.path.abspath(destination)
    if source == destination:
        return
    Path(destination).parent.mkdir(parents=True, exist_ok=True)
    temporary = destination + f".tmp.{os.getpid()}"
    try:
        shutil.copy2(source, temporary)
        os.replace(temporary, destination)
    finally:
        try:
            os.remove(temporary)
        except FileNotFoundError:
            pass


def is_minset_bed(path: str) -> bool:
    with open(path, "rt") as handle:
        for raw in handle:
            if not raw.strip() or raw.startswith("#"):
                continue
            fields = raw.rstrip("\r\n").split("\t")
            if len(fields) < 9:
                fields = raw.split()
            return len(fields) >= 9 and fields[7] in MINSET_BED_CLASSES
    return False


def write_direct_records_fasta(
    records: Sequence[DirectSequence], path: str, prefix: str = "g1",
    include_sample_in_coordinate: bool = False,
) -> str:
    """Write extracted BED records with IDs accepted by minsetref_light."""
    Path(path).parent.mkdir(parents=True, exist_ok=True)
    temporary = path + f".tmp.{os.getpid()}"
    try:
        with open(temporary, "wt") as output:
            for index, record in enumerate(records, 1):
                sample_fields = record.sample.split("_", 1)
                sample = (
                    record.sample if len(sample_fields) == 2
                    else record.sample + "_h1"
                )
                safe_contig = re.sub(r"[^A-Za-z0-9.#:-]+", "_", record.contig)
                record_id = (
                    f"{prefix}{index}_{sample}_{safe_contig}_"
                    f"{record.start}_{record.end}"
                )
                coordinate_contig = (
                    f"{sample}:{record.contig}"
                    if include_sample_in_coordinate else record.contig
                )
                class_token = (
                    f" class={record.record_class}" if record.record_class else ""
                )
                output.write(
                    f">{record_id} {coordinate_contig}:{record.start}-"
                    f"{record.end} block={prefix} region={index} "
                    f"source_name={record.query_name} strand={record.strand}"
                    f"{class_token}\n"
                    f"{wrap_fasta(record.sequence)}\n"
                )
        os.replace(temporary, path)
    finally:
        try:
            os.remove(temporary)
        except FileNotFoundError:
            pass
    return path


def materialize_bed_input(
    path: str, query_paths: str, output: str,
    include_sample_in_coordinate: bool = False,
) -> str:
    return write_direct_records_fasta(
        read_direct_bed(path, query_paths), output,
        include_sample_in_coordinate=include_sample_in_coordinate,
    )


def materialize_minset_bed(
    bed_path: str, source_fasta: str, output_fasta: str,
) -> str:
    """Extract local-coordinate minset BED rows from their source multi-FASTA."""
    rows: List[Tuple[str, int, int, str, str, str, str]] = []
    with open(bed_path, "rt") as handle:
        for line_number, raw in enumerate(handle, 1):
            if not raw.strip() or raw.startswith("#"):
                continue
            fields = raw.rstrip("\r\n").split("\t")
            if len(fields) < 9:
                fields = raw.split()
            if len(fields) < 9:
                raise ValueError(f"{bed_path}:{line_number}: expected minset BED9")
            rows.append((
                fields[0], int(fields[1]), int(fields[2]),
                fields[5], fields[6], fields[7], fields[3],
            ))
    with IndexedFasta(source_fasta) as source:
        records = []
        for index, (contig, start, end, strand, sample, kind, name) in enumerate(rows, 1):
            if contig not in source.index:
                raise KeyError(f"{bed_path}: {contig!r} is absent from {source_fasta}")
            sequence = source.fetch(contig, start, end)
            source_contig, source_start, source_end = contig, start, end
            match = UNIQUE_REGION_RE.fullmatch(name)
            if match is not None:
                source_contig, source_start, source_end = (
                    match.group(1), int(match.group(2)), int(match.group(3))
                )
                if source_end - source_start != len(sequence):
                    raise ValueError(
                        f"{bed_path}: absolute interval in {name!r} does not "
                        "match the retained sequence length"
                    )
            records.append(DirectSequence(
                sample, source_contig, source_start, source_end, strand,
                f"minset_{kind}_{index:08d}", sequence, kind,
            ))
    return write_direct_records_fasta(
        records, output_fasta, include_sample_in_coordinate=True,
    )


def run_checked(command: Sequence[str], label: str) -> None:
    LOG.info("%s", label)
    LOG.debug("Command: %s", " ".join(map(str, command)))
    subprocess.run(list(command), check=True)


def run_minsetref_light(
    args: argparse.Namespace, source_fasta: str, work_root: str,
) -> Tuple[str, str]:
    suffix = ".minset.bed" if args.minset_format == "bed" else ".minset.fasta"
    minset_output = os.path.abspath(args.output + suffix)
    minset_work = minset_output + ".work"
    if args.resume and output_is_current(minset_output, [source_fasta]):
        LOG.info("RESUME: reusing minset output: %s", minset_output)
    else:
        command = [
            sys.executable, str(HERE / "minsetref_light.py"),
            "-i", source_fasta,
            "-o", minset_output,
            "--output-format", args.minset_format,
            "-c", str(args.threads),
            "--kmer-preclean",
            "--kmer-ref-cleaner", args.kmer_ref_cleaner,
            "--workdir", minset_work,
            "--min-score", str(args.min_score),
            "--anchor", str(args.anchor),
            "--min-identity", str(args.min_identity),
            "--max-cycles", str(args.max_cycles),
        ]
        for selector in args.reference:
            command.extend(["-r", selector])
        if args.priority:
            command.extend(["-p", args.priority])
        if args.skip_blastn:
            command.append("--skip-blastn")
        if args.resume:
            command.append("--resume")
        run_checked(command, "Running minsetref_light with KmerRefCleaner pre-clean")

    if args.minset_format == "fasta":
        return minset_output, minset_output
    cleaned_fasta = os.path.join(work_root, "minset_retained_samples.fasta")
    if not (args.resume and output_is_current(cleaned_fasta, [minset_output, source_fasta])):
        materialize_minset_bed(minset_output, source_fasta, cleaned_fasta)
    return cleaned_fasta, minset_output


def graphmake_candidates(configured: str = "") -> Iterable[str]:
    if configured:
        yield configured
    environment = os.environ.get("MINSETREF_GRAPHMAKE")
    if environment:
        yield environment
    yield str(HERE / "graphmake.py")
    yield str(HERE.parent / "Graph" / "graphmake.py")
    yield str(HERE.parent / "graphmake.py")
    yield str(HERE.parent / "Ctyper2" / "graphmake.py")


def resolve_graphmake(configured: str = "") -> str:
    path = next(
        (os.path.abspath(os.path.expanduser(value)) for value in graphmake_candidates(configured)
         if value and os.path.isfile(os.path.abspath(os.path.expanduser(value)))),
        None,
    )
    if path is None:
        raise FileNotFoundError(
            "graphmake.py was not found; use --graphmake or MINSETREF_GRAPHMAKE"
        )
    return path


def run_graphmake(
    args: argparse.Namespace, cleaned_fasta: str, work_root: str,
) -> str:
    legacy = resolve_graphmake(args.graphmake)
    generated = cleaned_fasta + "_graph.FA"
    if args.resume and output_is_current(args.output, [cleaned_fasta]):
        LOG.info("RESUME: reusing graph: %s", args.output)
        return os.path.abspath(args.output)
    graph_work = os.path.join(work_root, "graphmake") + os.sep
    Path(graph_work).mkdir(parents=True, exist_ok=True)
    command = [
        sys.executable, legacy,
        "-i", cleaned_fasta,
        "-c", "1", "-f", "1", "-a", "0", "-l", "0", "-s", "1",
        "-t", str(args.threads),
        "-d", graph_work,
        "-u", "1" if args.resume else "0",
        "-q", str(args.highqual),
    ]
    if args.priority:
        command.extend(["-p", args.priority])
    run_checked(command, f"Building graph with {legacy}")
    if not os.path.isfile(generated) or os.path.getsize(generated) == 0:
        raise RuntimeError(f"graphmake did not create expected graph: {generated}")
    atomic_copy(generated, args.output)
    return os.path.abspath(args.output)


def build_graph(
    args: argparse.Namespace, cleaned_fasta: str, work_root: str,
) -> str:
    """Create the graph .FA directly or through the legacy fine-graph builder."""
    if args.graph_builder == "legacy":
        return run_graphmake(args, cleaned_fasta, work_root)
    if not (args.resume and output_is_current(args.output, [cleaned_fasta])):
        atomic_copy(cleaned_fasta, args.output)
    else:
        LOG.info("RESUME: reusing graph: %s", args.output)
    LOG.info(
        "Built minimum-set graph directly from retained FASTA paths: %s",
        args.output,
    )
    return os.path.abspath(args.output)


def write_adjacent_reference_bed(graph_fasta: str) -> Optional[str]:
    """Derive the reference seed BED from graph path source metadata."""
    graph_name = os.path.basename(graph_fasta)
    if graph_name.lower().endswith(".fa"):
        graph_name = graph_name[:-3]
    bed_path = os.path.join(os.path.dirname(graph_fasta), graph_name + ".bed")
    if os.path.isfile(bed_path) and os.path.getsize(bed_path) > 0:
        return bed_path
    rows = []
    with open(graph_fasta, "rt") as handle:
        for line_number, raw in enumerate(handle, 1):
            if not raw.startswith(">"):
                continue
            fields = raw[1:].strip().split()
            if len(fields) < 2:
                continue
            source = SOURCE_COORDINATE_RE.fullmatch(fields[1])
            if source is None:
                continue
            haplotype, contig, start_text, end_text, strand = source.groups()
            if "class=reference" not in fields[2:]:
                continue
            start, end = int(start_text), int(end_text)
            rows.append((contig, start, end, graph_name, 1000, strand, haplotype))
    if not rows:
        return None
    temporary = bed_path + f".tmp.{os.getpid()}"
    with open(temporary, "wt") as output:
        for row in rows:
            output.write("\t".join(map(str, row)) + "\n")
    os.replace(temporary, bed_path)
    return bed_path


def run_alignment(
    args: argparse.Namespace, graph_fasta: str, original_input: str,
) -> str:
    alignment = os.path.abspath(args.align or args.alignment)
    if args.align:
        fast_mode = bool(getattr(args, "fast_mode", False))
        slow_rigorous = bool(getattr(args, "slow_rigorous", False))
        dependencies = [graph_fasta, original_input, args.query_paths]
        mode_marker = alignment + ".mode"
        if slow_rigorous:
            mode_protocol = "slow-rigorous-v1;blastn=all;winnowmap=all\n"
        elif fast_mode:
            mode_protocol = "fast-mode-v2;max-unmapped-gap=50\n"
        else:
            mode_protocol = "normal-mode-v1\n"
        mode_protocol = mode_protocol.rstrip('\n') + f';payload={int(bool(args.alignment_payload))}\n'
        marker_current = False
        try:
            with open(mode_marker, "rt") as marker_handle:
                marker_current = marker_handle.read() == mode_protocol
        except FileNotFoundError:
            # Both alignment mode and saved payload policy must be known.
            marker_current = False
        if (
            args.resume and marker_current
            and output_is_current(alignment, dependencies)
        ):
            LOG.info("RESUME: reusing alignment: %s", alignment)
            return alignment
        input_format = args.alignment_source_format
        command = [
            sys.executable, str(HERE / "align_partition_hotspots.py"),
            "-i", original_input,
            "--input-format", input_format,
            "--graph", graph_fasta,
            "-o", alignment,
            "-t", str(args.threads),
            "--threads-per-job", str(args.threads_per_job),
            "--timeout", str(args.timeout),
            "--header",
        ]
        if fast_mode:
            command.append("--fast-mode")
        elif slow_rigorous:
            command.append("--slow-rigorous")
        else:
            command.append("--cycle")
        if args.jobs is not None:
            command.extend(["-j", str(args.jobs)])
        if args.query_paths:
            command.extend(["-q", args.query_paths])
        if args.sample:
            command.extend(["-s", args.sample])
        if args.keep_work:
            command.append("--keep-work")
        command.append(
            "--alignment-payload"
            if args.alignment_payload
            else "--no-alignment-payload"
        )
        if args.verbose:
            command.append("-v")
        label = (
            "Fast-aligning all original sample records to the graph"
            if fast_mode else
            "Slow-rigorous alignment of every original sample record with "
            "both BLASTN and Winnowmap"
            if slow_rigorous else
            "Aligning all original sample records to the graph until convergence"
        )
        run_checked(command, label)
        temporary_marker = mode_marker + f".tmp.{os.getpid()}"
        with open(temporary_marker, "wt") as marker_handle:
            marker_handle.write(mode_protocol)
        os.replace(temporary_marker, mode_marker)
    elif not os.path.isfile(alignment):
        raise FileNotFoundError(alignment)
    return alignment


def save_graphdb_atomic(graphdb, path: str) -> None:
    absolute = os.path.abspath(path)
    Path(os.path.dirname(absolute) or ".").mkdir(parents=True, exist_ok=True)
    temporary = absolute + f".tmp.{os.getpid()}"
    graphdb.save_json(temporary)
    os.replace(temporary, absolute)


def build_or_load_graphdb(
    args: argparse.Namespace,
    graph_fasta: str,
    alignment_path: str,
    rows,
):
    module = load_graphfixbreaks(args.graphfixbreaks)
    if args.cache_load:
        graphdb = module.graphDB.load_json(args.cache_load)
        LOG.info("Loaded uniform-break cache: %s", args.cache_load)
    else:
        initial_path = args.fixbreaks_init
        generated_initial = False
        if not initial_path:
            initial_path = args.fixbreaks + ".initial.tsv"
            segments = build_partition_segments(
                rows, graph_fasta, args.reference_bed,
                args.reference_haplotype,
                args.scutoff, args.dcutoff, args.min_segment_size,
                args.valid_min_alignment, args.valid_merge_gap,
                args.onlyref,
            )
            write_final_segments(initial_path, segments)
            generated_initial = True
            LOG.info(
                "Automatically determined %d conservative initial blocks: %s",
                len(segments), initial_path,
            )
        blocks = read_initial_blocks(initial_path)
        if not blocks:
            raise ValueError(f"initial block file contains no blocks: {initial_path}")
        graphdb, path_lengths, usable, skipped, uniform_count = uniformbreaks_from_blocks(
            module, rows, blocks, graph_fasta,
        )
        graphdb.cache_metadata.update({
            "format": "minsetgraphmake-gfixbreaks-cache-v1",
            "coordinate_system": "0-based-half-open",
            "validregions_semantics": (
                "reference-block-graph-closure-v1"
                if args.onlyref else
                "cohort-source-block-graph-closure-v1"
            ),
            "graph_fasta": os.path.abspath(graph_fasta),
            "reference_alignment": os.path.abspath(alignment_path),
            "reference_blocks": os.path.abspath(initial_path),
            "graph_path_lengths": path_lengths,
            "usable_reference_block_count": usable,
            "skipped_reference_block_count": skipped,
            "uniform_reference_region_count": uniform_count,
            "automatic_initial_blocks": generated_initial,
            "alignment_scope": (
                "only_reference" if args.onlyref else "all_samples"
            ),
        })
    if args.cache_save:
        save_graphdb_atomic(graphdb, args.cache_save)
        LOG.info("Saved uniform-break cache: %s", args.cache_save)
    return module, graphdb


def write_fixed_blocks(
    args: argparse.Namespace,
    graph_fasta: str,
    alignment_path: str,
) -> str:
    rows = select_alignment_rows(
        read_alignment_output(alignment_path),
        args.reference_haplotype or ("CHM13_h1",),
        args.onlyref,
    )
    if args.onlyref and not rows:
        raise ValueError("--onlyref selected no alignment rows")
    module, graphdb = build_or_load_graphdb(
        args, graph_fasta, alignment_path, rows,
    )
    graph_name = os.path.basename(graph_fasta)
    if graph_name.lower().endswith(".fa"):
        graph_name = graph_name[:-3]
    blocked = []
    for input_index, row in enumerate(rows):
        row_blocks, _links = block_one_alignment(
            module, graphdb, graph_name, input_index, row,
            args.merge_small,
        )
        blocked.extend(row_blocks)
    blocked.sort(key=lambda value: (
        value.contig, value.start, value.end, value.input_index,
        value.region_index,
    ))
    output_path = os.path.abspath(args.fixbreaks)
    Path(os.path.dirname(output_path) or ".").mkdir(parents=True, exist_ok=True)
    temporary = output_path + f".tmp.{os.getpid()}"
    try:
        with open(temporary, "wt") as output:
            for row in blocked:
                output.write(row.line())
        os.replace(temporary, output_path)
    finally:
        try:
            os.remove(temporary)
        except FileNotFoundError:
            pass
    LOG.info("Wrote %d uniform graph blocks: %s", len(blocked), output_path)
    return output_path


def linear_query_fasta(
    args: argparse.Namespace, original_input: str, work_root: str,
) -> str:
    if args.alignment_source_format == "fasta":
        return original_input
    records = read_direct_bed(original_input, args.query_paths)
    path = os.path.join(work_root, "linear_samples.fasta")
    Path(path).parent.mkdir(parents=True, exist_ok=True)
    with open(path, "wt") as output:
        for record in records:
            output.write(
                f">{record.query_name} {record.contig}:{record.start}-"
                f"{record.end}{record.strand}\n{wrap_fasta(record.sequence)}\n"
            )
    return path


def determine_source_format(args: argparse.Namespace) -> str:
    if args.input_format != "auto":
        return args.input_format
    if args.input.endswith(".FA"):
        return "graph"
    lowered = args.input.lower()
    if lowered.endswith(".adjusted.tsv"):
        return "adjusted"
    if lowered.endswith(".bed"):
        return "bed"
    if lowered.endswith((".fa", ".fasta", ".fna", ".fas")):
        return "fasta"
    raise ValueError(
        f"cannot infer input type from {args.input!r}; use --input-format"
    )


def determine_alignment_source_format(args: argparse.Namespace) -> str:
    """Infer the direct-alignment input independently of the graph input."""
    if args.align_input_format != "auto":
        return args.align_input_format
    if not args.align_input and getattr(args, "source_format", "") in {
        "fasta", "bed", "graph",
    }:
        return "bed" if args.source_format == "bed" else "fasta"
    path = args.align_input or args.input
    if getattr(args, "source_format", "") == "adjusted" and not args.align_input:
        return "fasta"
    if path.lower().endswith(".bed"):
        return "bed"
    if path.lower().endswith((".fa", ".fasta", ".fna", ".fas")):
        return "fasta"
    raise ValueError(
        f"cannot infer alignment input type from {path!r}; "
        "use --align-input-format"
    )


def materialize_hotspot_header(
    args: argparse.Namespace, work_root: str,
) -> Tuple[str, str]:
    """Stream one header manifest into a short-lived partition FASTA."""
    scratch_parent = (
        work_root if args.keep_work else os.environ.get("SLURM_TMPDIR", work_root)
    )
    Path(scratch_parent).mkdir(parents=True, exist_ok=True)
    scratch = tempfile.mkdtemp(prefix="minsetref_partition_", dir=scratch_parent)
    samples_fasta = os.path.join(scratch, "input_samples.fasta")
    command = [
        args.sample_fasta_builder,
        "--materialize-header", args.hotspot_header,
        "-q", args.query_paths,
        "--output-fasta", samples_fasta,
        "-j", str(max(1, args.threads)),
    ]
    run_checked(command, "Materializing one partition from its hotspot header")
    if not os.path.isfile(samples_fasta) or os.path.getsize(samples_fasta) == 0:
        raise RuntimeError(
            f"sample FASTA builder produced no records: {samples_fasta}"
        )
    return samples_fasta, scratch


def read_block_templates(path: str):
    """Read one partition's small block-template FASTA into a source index."""
    output = {}
    description: Optional[str] = None
    pieces: List[str] = []

    def finish() -> None:
        if description is None:
            return
        fields = description.split()
        attributes = {
            key: value
            for field in fields[1:] if "=" in field
            for key, value in (field.split("=", 1),)
        }
        coordinate = next((
            field for field in fields[1:]
            if re.fullmatch(r"[^:]+:\d+-\d+", field)
        ), None)
        if coordinate is None or "sample" not in attributes:
            raise ValueError(
                f"{path}: block-template header lacks source coordinates: "
                f"{description!r}"
            )
        contig, interval = coordinate.rsplit(":", 1)
        start_text, end_text = interval.split("-", 1)
        key = (
            attributes["sample"], contig, int(start_text), int(end_text),
        )
        if key in output:
            raise ValueError(f"{path}: duplicate block-template source {key}")
        output[key] = (
            "".join(pieces), attributes.get("strand", "+"),
        )

    with open(path, "rt") as source:
        for raw in source:
            if raw.startswith(">"):
                finish()
                description = raw[1:].strip()
                pieces = []
            else:
                pieces.append(raw.strip())
    finish()
    if not output:
        raise ValueError(f"{path}: no block-template records")
    return output


def write_adjusted_graph_fasta(
    args: argparse.Namespace,
    partition: str,
    rows,
    generated: str,
) -> None:
    """Install refined paths, retaining original templates when unrefined."""
    has_refined_reference = any(row["role"] == "reference" for row in rows)
    temporary = args.output + f".tmp.{os.getpid()}"
    try:
        with open(temporary, "wt") as destination:
            if not has_refined_reference:
                partition_template = os.path.join(
                    os.path.dirname(args.hotspot_header),
                    f"{partition}.fasta",
                )
                legacy_template = (
                    os.path.splitext(args.hotspot_header)[0] + ".fasta"
                )
                block_fasta = (
                    partition_template
                    if os.path.isfile(partition_template)
                    else legacy_template
                )
                if not os.path.isfile(block_fasta):
                    raise FileNotFoundError(
                        f"{partition}: unrefined graph requires its block "
                        "template FASTA; checked "
                        f"{partition_template} and {legacy_template}"
                    )
                templates = read_block_templates(block_fasta)
                for row in rows:
                    if row["role"] != "original":
                        continue
                    key = (
                        row["source_haplotype"], row["source_contig"],
                        int(row["source_start"]), int(row["source_end"]),
                    )
                    if key not in templates:
                        raise KeyError(
                            f"{partition}: original template {key} is absent "
                            f"from {block_fasta}"
                        )
                    sequence, strand = templates[key]
                    source = (
                        f"{row['source_haplotype']}:{row['source_contig']}:"
                        f"{row['source_start']}-{row['source_end']}{strand}"
                    )
                    local = (
                        f"{row['source_contig']}:{row['source_start']}-"
                        f"{row['source_end']}{strand}"
                    )
                    destination.write(
                        f">{public_path_id(row)} {source} {local} . {local} "
                        f"reference\n{wrap_fasta(sequence)}\n"
                    )
            if os.path.isfile(generated) and os.path.getsize(generated):
                with open(generated, "rt") as source:
                    shutil.copyfileobj(source, destination, 1024 * 1024)
        if os.path.getsize(temporary) == 0:
            raise RuntimeError(
                f"locally adjusted manifest produced no graph paths: {args.input}"
            )
        os.replace(temporary, args.output)
    finally:
        try:
            os.remove(temporary)
        except FileNotFoundError:
            pass


def build_adjusted_graph(
    args: argparse.Namespace,
    samples_fasta: str,
    work_root: str,
) -> str:
    partition = os.path.basename(os.path.dirname(args.input))
    expected = f"{partition}.adjusted.tsv"
    if os.path.basename(args.input) != expected:
        raise ValueError(
            f"adjusted manifest must be PARTITION/{expected}: {args.input}"
        )
    build_root = os.path.join(work_root, "locally_refined_graph")
    rows = read_adjusted_manifest(args.input, partition)
    from fixed_alternatives import has_fixed_templates, validate_fixed_paths
    template_fasta = Path(args.hotspot_header).parent / (partition + '.fasta')
    fixed = template_fasta.is_file() and has_fixed_templates(template_fasta)
    if fixed:
        validate_fixed_paths(template_fasta, args.input)
    _partition, records, _bases = build_one_adjusted_partition(
        (partition, args.input), build_root, args.resume, samples_fasta,
    )
    generated = os.path.join(
        build_root, partition, f"{partition}_samples.fasta",
    )
    if not (args.resume and output_is_current(args.output, [args.input, generated])):
        write_adjusted_graph_fasta(args, partition, rows, generated)
    if fixed:
        validate_fixed_paths(template_fasta, args.input, args.output)
    if not args.keep_work:
        os.remove(generated)
        try:
            os.rmdir(os.path.dirname(generated))
        except OSError:
            pass
    LOG.info(
        "Built immutable locally refined graph%s: %s",
        " from resumed adjusted paths" if records < 0 else "", args.output,
    )
    return os.path.abspath(args.output)


def run(args: argparse.Namespace) -> int:
    args.input = os.path.abspath(args.input)
    args.output = os.path.abspath(args.output)
    args.query_paths = (
        os.path.abspath(args.query_paths) if args.query_paths else None
    )
    args.align_input = (
        os.path.abspath(args.align_input) if args.align_input else None
    )
    args.materialized_samples = (
        os.path.abspath(args.materialized_samples)
        if args.materialized_samples else None
    )
    args.source_format = determine_source_format(args)
    args.alignment_source_format = determine_alignment_source_format(args)
    if args.source_format == "graph" and not args.resume:
        raise ValueError("existing graph .FA input requires --resume")
    if (
        args.source_format == "adjusted"
        and not args.hotspot_header
        and not args.materialized_samples
    ):
        raise ValueError(
            "adjusted manifest input requires --hotspot-header or "
            "--materialized-samples"
        )
    if args.hotspot_header and not args.query_paths:
        raise ValueError("--hotspot-header requires -q/--query-paths")
    if args.source_format == "bed" and not args.query_paths:
        raise ValueError("BED input requires -q/--query-paths")
    if args.alignment_source_format == "bed" and not args.query_paths:
        raise ValueError("BED alignment input requires -q/--query-paths")

    work_root = os.path.abspath(args.workdir or (args.output + ".work"))
    Path(work_root).mkdir(parents=True, exist_ok=True)
    materialized_scratch: Optional[str] = None
    materialized_samples: Optional[str] = args.materialized_samples
    if materialized_samples:
        if not os.path.isfile(materialized_samples) or not os.path.getsize(
            materialized_samples
        ):
            raise FileNotFoundError(
                f"materialized sample FASTA is missing or empty: "
                f"{materialized_samples}"
            )
        LOG.info(
            "Using caller-owned materialized partition FASTA: %s",
            materialized_samples,
        )
        if args.hotspot_header:
            args.hotspot_header = os.path.abspath(args.hotspot_header)
    elif args.hotspot_header:
        args.hotspot_header = os.path.abspath(args.hotspot_header)
        materialized_samples, materialized_scratch = materialize_hotspot_header(
            args, work_root,
        )
    original_input = args.align_input or materialized_samples or args.input

    succeeded = False
    try:
        if args.source_format == "graph":
            atomic_copy(args.input, args.output)
            graph_fasta = args.output
            cleaned_fasta = args.input
        elif args.source_format == "adjusted":
            assert materialized_samples is not None
            cleaned_fasta = materialized_samples
            graph_fasta = build_adjusted_graph(
                args, materialized_samples, work_root,
            )
        else:
            if args.source_format == "bed":
                precleaned_bed = is_minset_bed(args.input)
                samples_fasta = os.path.join(work_root, "input_samples.fasta")
                if not (args.resume and output_is_current(samples_fasta, [args.input, args.query_paths])):
                    materialize_bed_input(
                        args.input, args.query_paths, samples_fasta,
                        include_sample_in_coordinate=precleaned_bed,
                    )
                source_fasta = samples_fasta
            else:
                precleaned_bed = False
                source_fasta = args.input

            if args.source_format == "bed" and precleaned_bed:
                LOG.info("RESUME: input BED is already a minsetref_light BED")
                cleaned_fasta = source_fasta
            else:
                cleaned_fasta, _minset_output = run_minsetref_light(
                    args, source_fasta, work_root,
                )
            graph_fasta = build_graph(args, cleaned_fasta, work_root)

        adjacent_bed = write_adjacent_reference_bed(graph_fasta)
        if adjacent_bed:
            LOG.info("Reference seed BED: %s", adjacent_bed)

        alignment_path: Optional[str] = None
        if args.align or args.alignment or args.fixbreaks or args.linear:
            if not (args.align or args.alignment):
                raise ValueError(
                    "--fixbreaks/--linear requires --align OUTPUT or --alignment EXISTING"
                )
            alignment_path = run_alignment(
                args, graph_fasta, original_input,
            )

        if args.fixbreaks:
            assert alignment_path is not None
            write_fixed_blocks(args, graph_fasta, alignment_path)

        if args.linear:
            assert alignment_path is not None
            query_fasta = linear_query_fasta(args, original_input, work_root)
            local_graph_to_linear(
                alignment_path, query_fasta, graph_fasta,
                os.path.abspath(args.linear), args.threads,
                os.path.join(work_root, "stretcherouts") + os.sep,
                args.polish,
            )
            LOG.info("Wrote MSA-linearized graph: %s", os.path.abspath(args.linear))
        succeeded = True
    finally:
        if materialized_scratch and succeeded and not args.keep_work:
            shutil.rmtree(materialized_scratch)
            LOG.info("Deleted short-lived partition sample FASTA: %s", materialized_samples)
        elif materialized_scratch and not succeeded:
            LOG.error(
                "Retained failed partition materialization directory: %s",
                materialized_scratch,
            )
    return 0


def parse_args(argv: Optional[Sequence[str]] = None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("-i", "--input", required=True)
    parser.add_argument("-o", "--output", required=True, help="output graph .FA")
    parser.add_argument(
        "--input-format", choices=("auto", "fasta", "bed", "graph", "adjusted"),
        default="auto",
    )
    parser.add_argument(
        "-q", "--query-paths",
        help="SAMPLE FASTA (with adjacent FASTA.fai) table required for BED input",
    )
    parser.add_argument("-r", "--reference", action="append", default=[])
    parser.add_argument("-p", "--priority", default="")
    parser.add_argument("-s", "--sample", help="fallback sample name for FASTA")
    parser.add_argument("-t", "--threads", type=int, default=8)
    parser.add_argument("-j", "--jobs", type=int)
    parser.add_argument("--threads-per-job", type=int, default=2)
    parser.add_argument(
        "--timeout", type=int, default=21_600,
        help="timeout per external alignment command (default: 21600, six hours)",
    )
    parser.add_argument(
        "--fast", "--fast-mode", dest="fast_mode", action="store_true",
        help=(
            "use bulk BLASTN plus 256-MiB Winnowmap query batches for "
            "--align; --fast-mode is a compatibility alias"
        ),
    )
    parser.add_argument(
        "--slow-rigorous", action="store_true",
        help=(
            "use bulk BLASTN and Winnowmap for every original sample record "
            "during --align"
        ),
    )
    parser.add_argument("--resume", action="store_true")
    parser.add_argument("--workdir")
    parser.add_argument("--keep-work", action="store_true")
    parser.add_argument(
        "--hotspot-header",
        help="PARTITION_samples.header used to materialize one short-lived sample FASTA",
    )
    parser.add_argument(
        "--materialized-samples",
        help=(
            "caller-owned partition FASTA already materialized from the "
            "header; use it without creating or deleting another copy"
        ),
    )
    parser.add_argument(
        "--sample-fasta-builder",
        default=str(HERE / "build_local_graphs_parallel_hotspots"),
        help="C++ builder supporting --materialize-header",
    )
    payload_group = parser.add_mutually_exclusive_group()
    payload_group.add_argument(
        "--alignment-payload", dest="alignment_payload", action="store_true",
        default=False,
        help=(
            "include query sequence payloads on I/X operations in --align "
            "output; default is compact CIGARs with bases restored during calling"
        ),
    )
    payload_group.add_argument(
        "--no-alignment-payload", dest="alignment_payload", action="store_false",
        help="omit query sequence payloads from I/X operations in --align output",
    )

    parser.add_argument("--minset-format", choices=("fasta", "bed"), default="fasta")
    parser.add_argument(
        "--kmer-ref-cleaner", default=str(HERE / "KmerRefCleaner"),
    )
    parser.add_argument("--min-score", type=int, default=100)
    parser.add_argument("--anchor", type=int, default=50)
    parser.add_argument("--min-identity", type=float, default=95.0)
    parser.add_argument("--max-cycles", type=int, default=0)
    parser.add_argument("--skip-blastn", action="store_true")

    parser.add_argument("--graphmake", help="legacy graphmake.py path")
    parser.add_argument(
        "--graph-builder", choices=("direct", "legacy"), default="direct",
        help=(
            "direct treats minsetref_light's retained FASTA paths as the graph; "
            "legacy runs graphmake.py's additional fine-graph pass"
        ),
    )
    parser.add_argument("--highqual", type=int, choices=(0, 1), default=1)

    parser.add_argument("--align", metavar="TSV", help="align original input to graph")
    parser.add_argument(
        "--align-input", metavar="FASTA_OR_BED",
        help=(
            "original FASTA/BED to align when -i is a resumed graph; "
            "defaults to -i"
        ),
    )
    parser.add_argument(
        "--align-input-format", choices=("auto", "fasta", "bed"),
        default="auto",
    )
    parser.add_argument(
        "--alignment", metavar="TSV",
        help="reuse an existing ten-column alignment instead of running --align",
    )
    parser.add_argument("--fixbreaks", metavar="BED", help="uniform blocked BED output")
    parser.add_argument("--fixbreaks-init", metavar="TSV")
    parser.add_argument("--cache-load", metavar="JSON")
    parser.add_argument("--cache-save", metavar="JSON")
    parser.add_argument("--reference-bed")
    parser.add_argument("--reference-haplotype", action="append", default=[])
    parser.add_argument(
        "--onlyref", "--only-ref", action="store_true",
        help=(
            "restrict breakpoint-cache construction and blocking output to "
            "reference rows; default: use all samples"
        ),
    )
    parser.add_argument("--graphfixbreaks", default=str(HERE / "gfixbreaks.py"))
    parser.add_argument("--merge-small", type=int, default=20_000)
    parser.add_argument("--scutoff", type=int, default=15_000)
    parser.add_argument("--dcutoff", type=int, default=15_000)
    parser.add_argument("--min-segment-size", type=int, default=1_000)
    parser.add_argument("--valid-min-alignment", type=int, default=100)
    parser.add_argument("--valid-merge-gap", type=int, default=1_000)

    parser.add_argument("--linear", metavar="GAF", help="MSA-linearized graph output")
    parser.add_argument("--polish", type=int, choices=(0, 1), default=0)
    parser.add_argument("-v", "--verbose", action="store_true")
    args = parser.parse_args(argv)
    if args.threads < 1 or args.threads_per_job < 1:
        parser.error("thread counts must be positive")
    if args.jobs is not None and args.jobs < 1:
        parser.error("--jobs must be positive")
    if args.timeout < 1:
        parser.error("--timeout must be positive")
    if args.fast_mode and args.slow_rigorous:
        parser.error("--fast and --slow-rigorous are mutually exclusive")
    if args.min_score < 1 or args.anchor < 0 or args.max_cycles < 0:
        parser.error("invalid minset cleaning size/cycle option")
    if not 0 < args.min_identity <= 100:
        parser.error("--min-identity must be in (0,100]")
    if args.merge_small < 0:
        parser.error("--merge-small must be nonnegative")
    if args.cache_load and not args.fixbreaks:
        parser.error("--cache-load requires --fixbreaks")
    if args.cache_save and not args.fixbreaks:
        parser.error("--cache-save requires --fixbreaks")
    if args.fixbreaks_init and not args.fixbreaks:
        parser.error("--fixbreaks-init requires --fixbreaks")
    if args.align and args.alignment:
        parser.error("use only one of --align and --alignment")
    if os.path.abspath(args.input) == os.path.abspath(args.output) and not args.resume:
        parser.error("--input and --output must differ unless resuming a graph")
    return args


def main(argv: Optional[Sequence[str]] = None) -> int:
    args = parse_args(argv)
    logging.basicConfig(
        level=logging.DEBUG if args.verbose else logging.INFO,
        format="[%(asctime)s] %(levelname)s: %(message)s",
    )
    try:
        return run(args)
    except subprocess.CalledProcessError as error:
        if error.returncode == 124:
            LOG.warning(
                "Partition alignment timed out; returning status 124 so the "
                "cohort runner can record and skip this partition"
            )
            return 124
        LOG.error("External command failed (%d): %s", error.returncode, error.cmd)
        return 1
    except Exception as error:
        LOG.error("%s", error)
        if args.verbose:
            LOG.exception("minset graph pipeline failed")
        return 1


if __name__ == "__main__":
    raise SystemExit(main())
