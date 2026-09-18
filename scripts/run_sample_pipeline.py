#!/usr/bin/env python3
"""Run one assembly from partition hotspots through a single-sample VCF.

The driver keeps all filenames and sample labels consistent across:

  align_partition_hotspots.py
  summarize_partition_hotspot_segments.py
  block_partition_alignments.py
  match_partition_blocks.py
  GenomeLift.py
  local_reference_templates.py
  graphcigartoref_persample.py
  fill_graphcigartoref_gaps.py
  graphreftovcf_persample.py

The reference alignment and blocked BED are precomputed inputs shared by all
query samples. By default, each query hotspot is aligned with both BLASTN and
Winnowmap (``slow-rigorous`` mode). Existing nonempty query-stage outputs are skipped unless
``--rerun`` is supplied. With ``--resume``, only the consecutive completed
prefix is skipped; the first incomplete stage and every downstream stage are
run so newly generated intermediates cannot be combined with stale outputs.

Before alignment, the query assembly is inspected for lowercase repeat
masking and the expected sample/haplotype contig prefix.  Missing preparation
is performed with WindowMasker and namecontigsfix.py, and missing query or
reference FASTA indexes are created with samtools faidx.  When no hotspot
table is supplied, KmerSearcher generates one from the prepared query FASTA
and the same partition target list used by the alignment stages.

By default, complete local graphs with no requested-reference alignment are
called against their locus-template path. Use
``--no-local-template-fallback`` to disable this behavior.
"""

from __future__ import annotations

import argparse
import dataclasses
import gzip
import os
from pathlib import Path
import shlex
import shutil
import subprocess
import sys
import tempfile
from typing import Iterable, List, Optional, Sequence, Tuple

from graph_list_cache import (
    ensure_graph_list_cache, ensure_template_list_cache, graph_list_path,
    template_list_path,
)


HIDDEN_LINE_SEPARATORS = "\r\n\u2028\u2029"


@dataclasses.dataclass(frozen=True)
class Stage:
    name: str
    command: Tuple[str, ...]
    inputs: Tuple[Path, ...]
    outputs: Tuple[Path, ...]
    output_text: Optional[str] = None
    readiness_marker: Optional[Path] = None
    readiness_protocol: Optional[str] = None
    allow_missing_readiness_marker: bool = False
    write_readiness_marker: bool = False


@dataclasses.dataclass(frozen=True)
class PipelinePaths:
    reference_align: Path
    reference_blocks: Path
    query_align: Path
    query_summary: Path
    query_blocks: Path
    refmatch: Path
    genomelift: Path
    genomeliftfix: Path
    local_templates: Path
    graphcigartoref: Path
    graphcigartoreffix: Path
    vcf: Path
    genomelift_tmp: Path
    norm_folder: Path
    log_folder: Path


def clean_path(value: str, option: str) -> Path:
    cleaned = str(value).strip(HIDDEN_LINE_SEPARATORS)
    if any(character in cleaned for character in HIDDEN_LINE_SEPARATORS):
        raise ValueError(
            f"{option} contains an embedded newline or Unicode line separator"
        )
    if not cleaned:
        raise ValueError(f"{option} is empty")
    if cleaned != value:
        print(
            f"[partition_pipeline] warning: removed a hidden line separator "
            f"from {option}: {cleaned!r}",
            file=sys.stderr,
        )
    return Path(os.path.abspath(os.path.expanduser(cleaned)))


def output_paths(
    output_root: Path,
    sample: str,
    reference_align: Path,
    reference_blocks: Path,
) -> PipelinePaths:
    query_dir = output_root / sample
    query_align = query_dir / f"{sample}.align.txt"
    return PipelinePaths(
        reference_align=reference_align,
        reference_blocks=reference_blocks,
        query_align=query_align,
        query_summary=query_dir / f"{sample}.partition_segments.tsv",
        query_blocks=query_dir / f"{sample}.align.txt_blocks.bed",
        refmatch=query_dir / f"{sample}.refmatch.tsv",
        genomelift=query_dir / f"{sample}.genomelift.tsv",
        genomeliftfix=query_dir / f"{sample}.genomeliftfix.tsv",
        local_templates=query_dir / "local_reference_templates.fa",
        graphcigartoref=query_dir / f"{sample}.graphcigartoref.tsv",
        graphcigartoreffix=query_dir / f"{sample}.graphcigartoreffix.tsv",
        vcf=query_dir / f"{sample}.vcf",
        genomelift_tmp=query_dir / f"{sample}.GenomeLiftTemp",
        norm_folder=query_dir / f"{sample}.per_graph_norms",
        log_folder=query_dir / "pipeline_logs",
    )


def pipeline_script_path(script_dir: Path, script: str) -> Path:
    if script == "namecontigsfix.py":
        return script_dir.parent / "tools" / script
    return script_dir / script


def script_command(script_dir: Path, script: str, *arguments: object) -> Tuple[str, ...]:
    return (
        sys.executable,
        str(pipeline_script_path(script_dir, script)),
        *(str(argument) for argument in arguments),
    )


def prefer_local_executable(executable: str, script_dir: Path) -> str:
    """Prefer an executable beside this driver before searching ``PATH``."""
    if os.path.sep in executable:
        return executable
    candidates = (
        ("KmerSearcher", "kmer_searcher")
        if executable in {"KmerSearcher", "kmer_searcher"}
        else (executable,)
    )
    for candidate in candidates:
        local = script_dir / candidate
        if local.is_file() and os.access(local, os.X_OK):
            return str(local)
    if executable in {"KmerSearcher", "kmer_searcher"}:
        # Never silently use a stale user-wide KmerSearcher. Returning the
        # expected local path makes preflight report the missing installation.
        return str(script_dir / candidates[0])
    for candidate in candidates:
        found = shutil.which(candidate)
        if found is not None:
            return found
    return executable


def fasta_stem(path: Path) -> str:
    """Strip common FASTA/compression suffixes from one filename."""
    name = path.name
    if name.lower().endswith(".gz"):
        name = name[:-3]
    for suffix in (".fasta", ".fna", ".fa"):
        if name.lower().endswith(suffix):
            return name[:-len(suffix)]
    return name


def contig_prefix_for_sample(sample: str) -> Optional[str]:
    """Return PAT-style prefix, preserving coordinate-named references."""
    if sample.startswith(("CHM13_", "HG38_")):
        return None
    return sample.replace("_h", "#")


def inspect_query_fasta(path: Path, prefix: Optional[str]) -> Tuple[bool, bool]:
    """Return ``(has_lowercase_masking, headers_are_normalized)``."""
    opener = gzip.open if path.name.lower().endswith(".gz") else open
    has_lowercase_masking = False
    headers_are_normalized = True
    saw_header = False

    with opener(path, "rt") as handle:
        for line_number, line in enumerate(handle, 1):
            if line.startswith(">"):
                saw_header = True
                fields = line[1:].split()
                if not fields:
                    raise ValueError(f"{path}:{line_number}: empty FASTA header")
                if prefix is not None and not (
                    fields[0] == prefix or fields[0].startswith(prefix + "#")
                ):
                    headers_are_normalized = False
            elif any(base in line for base in "acgt"):
                has_lowercase_masking = True

    if not saw_header:
        raise ValueError(f"no FASTA records found: {path}")
    return has_lowercase_masking, headers_are_normalized


def build_fasta_preparation_stages(
    args: argparse.Namespace,
    script_dir: Path,
) -> Tuple[Path, List[Stage]]:
    """Plan masking, idempotent contig naming, and FASTA indexing."""
    source = clean_path(args.fasta_query, "--fasta-query")
    if not source.is_file():
        raise FileNotFoundError(source)

    output_root = clean_path(args.output_root, "--output-root")
    sample_dir = output_root / args.sample
    stem = fasta_stem(source)
    final_local = sample_dir / f"{stem}_masked.fa"
    maskinfo = sample_dir / f"{stem}_maskinfo"
    prefix = contig_prefix_for_sample(args.sample)
    is_masked, headers_ready = inspect_query_fasta(source, prefix)
    stages: List[Stage] = []
    naming_input = source

    if not is_masked:
        mask_output = (
            sample_dir / f"{stem}_masked.unnamed.fa"
            if not headers_ready and prefix is not None
            else final_local
        )
        stages.extend((
            Stage(
                name="query_mask_counts",
                command=(
                    args.windowmasker,
                    "-mk_counts", "-in", str(source), "-out", str(maskinfo),
                ),
                inputs=(source,),
                outputs=(maskinfo,),
            ),
            Stage(
                name="query_mask",
                command=(
                    args.windowmasker,
                    "-ustat", str(maskinfo),
                    "-dust", "yes",
                    "-outfmt", "fasta",
                    "-in", str(source),
                    "-out", str(mask_output),
                ),
                inputs=(source, maskinfo),
                outputs=(mask_output,),
            ),
        ))
        naming_input = mask_output

    if prefix is not None and not headers_ready:
        stages.append(Stage(
            name="query_contig_names",
            command=script_command(
                script_dir,
                "namecontigsfix.py",
                "-i", naming_input,
                "-n", prefix,
                "-o", final_local,
            ),
            inputs=(naming_input,),
            outputs=(final_local,),
        ))
        prepared = final_local
    else:
        prepared = naming_input

    prepared_index = Path(str(prepared) + ".fai")
    stages.append(Stage(
        name="query_faidx",
        command=(args.samtools, "faidx", str(prepared)),
        inputs=(prepared,),
        outputs=(prepared_index,),
    ))
    return prepared, stages


def build_reference_index_stage(args: argparse.Namespace) -> Stage:
    reference = clean_path(args.fasta_reference, "--fasta-reference")
    return Stage(
        name="reference_faidx",
        command=(args.samtools, "faidx", str(reference)),
        inputs=(reference,),
        outputs=(Path(str(reference) + ".fai"),),
    )


def build_hotspot_stages(
    args: argparse.Namespace,
    prepared_query: Path,
) -> Tuple[Path, List[Stage]]:
    """Use a supplied hotspot table or plan one KmerSearcher invocation."""
    if args.hotspot_query:
        return clean_path(args.hotspot_query, "--hotspot-query"), []

    output_root = clean_path(args.output_root, "--output-root")
    graph_list = clean_path(args.graph_list, "--graph-list")
    template_list = (
        clean_path(args.template_list, "--template-list")
        if getattr(args, "template_list", None) else template_list_path(graph_list)
    )
    template_cache = Path(str(template_list) + ".bin")
    sample_dir = output_root / args.sample
    query_list = sample_dir / f"{args.sample}.kmer_searcher.query.list"
    hotspot = sample_dir / f"{args.sample}_hotspot.txt"
    query_list_text = f"{args.sample}\t{prepared_query}\n"

    return hotspot, [
        Stage(
            name="hotspot_query_list",
            command=("<write-query-list>",),
            inputs=(prepared_query,),
            outputs=(query_list,),
            output_text=query_list_text,
        ),
        Stage(
            name="hotspot_search",
            command=(
                args.kmer_searcher,
                "-T", str(template_list),
                "-I", str(query_list),
                "-o", str(sample_dir) + os.sep,
                "-c", str(args.hotspot_cutoff),
                "-w", str(args.hotspot_window),
                "-N", str(args.threads),
            ),
            inputs=(
                prepared_query, query_list, template_list, template_cache,
            ),
            outputs=(hotspot,),
            readiness_marker=Path(str(hotspot) + '.targets.mode'),
            readiness_protocol=(
                f'template-hotspots-v1\n{template_list}\n'
                f'cutoff={args.hotspot_cutoff};window={args.hotspot_window}\n'
            ),
            write_readiness_marker=True,
        ),
    ]


def build_stages(args: argparse.Namespace, script_dir: Path) -> Tuple[PipelinePaths, List[Stage]]:
    output_root = clean_path(args.output_root, "--output-root")
    graph_folder = clean_path(args.graph_folder, "--graph-folder")
    graph_list = clean_path(args.graph_list, "--graph-list")
    fasta_query = clean_path(args.fasta_query, "--fasta-query")
    fasta_reference = clean_path(args.fasta_reference, "--fasta-reference")
    fasta_query_index = Path(str(fasta_query) + ".fai")
    fasta_reference_index = Path(str(fasta_reference) + ".fai")
    hotspot_query = clean_path(args.hotspot_query, "--hotspot-query")
    reference_align = clean_path(args.reference_align, "--reference-align")
    reference_blocks = clean_path(args.reference_blocks, "--reference-blocks")
    paths = output_paths(
        output_root, args.sample, reference_align, reference_blocks,
    )

    align_threads = args.align_threads or args.threads
    summary_jobs = args.summary_jobs or min(args.threads, 16)
    block_jobs = args.block_jobs or args.threads
    match_jobs = args.match_jobs or args.threads

    def align_stage(name: str, hotspot: Path, fasta: Path, sample: str, output: Path) -> Stage:
        align_arguments: List[object] = [
            "-i", hotspot,
            "-q", fasta,
            "-L", graph_list,
            "-G", graph_folder,
            "-s", sample,
            "-o", output,
            "-t", align_threads,
            "--threads-per-job", args.threads_per_job,
            "--extension", args.extension,
            "--merge-gap", args.merge_gap,
            "--timeout", args.timeout,
            "--header",
        ]
        align_arguments.append("--fast" if args.fast else "--slow-rigorous")
        alignment_payload = bool(getattr(args, "alignment_payload", False))
        align_arguments.append(
            "--alignment-payload" if alignment_payload else "--no-alignment-payload"
        )
        alignment_protocol = (
            "fast-mode-v2;max-unmapped-gap=50\n" if args.fast else
            "slow-rigorous-v1;blastn=all;winnowmap=all\n"
        )
        alignment_protocol = alignment_protocol.rstrip('\n') + f';payload={int(alignment_payload)}\n'
        return Stage(
            name=name,
            command=script_command(
                script_dir,
                "align_partition_hotspots.py",
                *align_arguments,
            ),
            inputs=(
                hotspot, fasta, Path(str(fasta) + ".fai"),
                graph_list, graph_folder,
            ),
            outputs=(output,),
            readiness_marker=Path(str(output) + ".mode"),
            readiness_protocol=alignment_protocol,
            # Neither supported mode may reuse an unmarked alignment because
            # its protocol cannot be verified.
            allow_missing_readiness_marker=False,
        )

    def summary_stage(name: str, alignment: Path, output: Path) -> Stage:
        return Stage(
            name=name,
            command=script_command(
                script_dir,
                "summarize_partition_hotspot_segments.py",
                "-i", alignment,
                "-G", graph_folder,
                "-L", graph_list,
                "-o", output,
                "-j", summary_jobs,
                "--reference-haplotype", args.reference_sample,
            ),
            inputs=(alignment, graph_list, graph_folder),
            outputs=(output,),
        )

    def block_stage(name: str, alignment: Path, output: Path) -> Stage:
        return Stage(
            name=name,
            command=script_command(
                script_dir,
                "block_partition_alignments.py",
                "-i", alignment,
                "-G", graph_folder,
                "-L", graph_list,
                "-o", output,
                "-j", block_jobs,
            ),
            inputs=(alignment, graph_list, graph_folder),
            outputs=(output,),
        )

    stages = [
        align_stage(
            "query_align", hotspot_query, fasta_query,
            args.sample, paths.query_align,
        ),
        summary_stage("query_summary", paths.query_align, paths.query_summary),
        block_stage("query_blocks", paths.query_align, paths.query_blocks),
        Stage(
            name="local_reference_templates",
            command=script_command(
                script_dir,
                "local_reference_templates.py",
                "--graph-folder", graph_folder,
                "--graph-list", graph_list,
                "--align-ref", reference_align,
                "--reference-haplotype", args.reference_sample,
                "--output", paths.local_templates.with_suffix('.unlifted.fa'),
                *(("--disabled",) if not args.local_template_fallback else ()),
            ),
            inputs=(reference_align, graph_list, graph_folder,
                    script_dir / 'local_reference_templates.py'),
            outputs=(
                paths.local_templates.with_suffix('.unlifted.fa'),
                Path(str(paths.local_templates.with_suffix('.unlifted.fa')) + ".fai"),
            ),
            readiness_marker=Path(str(paths.local_templates.with_suffix('.unlifted.fa'))+'.selection.complete'),
            readiness_protocol=f'local-template-selection-v4-fixed-alternatives:{args.reference_sample}:{reference_align}\n',
        ),
        Stage(
            name="lift_local_templates",
            command=script_command(
                script_dir, "lift_local_templates.py",
                "--input", paths.local_templates.with_suffix('.unlifted.fa'),
                "--reference-fasta", fasta_reference,
                "--reference-haplotype", args.reference_sample,
                "--source", args.sample, fasta_query,
                "--output", paths.local_templates, "--threads", args.threads,
                *(("--assemblies", args.template_assemblies) if getattr(args, 'template_assemblies', '') else ()),
            ),
            inputs=(paths.local_templates.with_suffix('.unlifted.fa'), fasta_reference,
                    fasta_query, script_dir / 'lift_local_templates.py',
                    script_dir / 'minsetref_localize.py', script_dir / 'assembly_contigs.py',
                    *((Path(args.template_assemblies),) if getattr(args, 'template_assemblies', '') else ())),
            outputs=(paths.local_templates, Path(str(paths.local_templates)+'.fai')),
            readiness_marker=Path(str(paths.local_templates)+'.lift.complete'),
            readiness_protocol=f'local-template-lift-v2-fixed-alternatives:{args.reference_sample}:{fasta_reference}\n',
        ),
    ]

    match_arguments: List[object] = [
        "--bed-query", paths.query_blocks,
        "--fasta-query", fasta_query,
        "--sample-query", args.sample,
        "--bed-ref", paths.reference_blocks,
        "--fasta-ref", fasta_reference,
        "--sample-ref", args.reference_sample,
        "--output", paths.refmatch,
        "--jobs", match_jobs,
        "--header",
    ]
    match_inputs = [
        paths.query_blocks, fasta_query, fasta_query_index,
        paths.reference_blocks, fasta_reference, fasta_reference_index,
    ]
    if args.kmermatch:
        kmermatch = clean_path(args.kmermatch, "--kmermatch")
        match_arguments.extend(("--kmermatch", kmermatch))
        match_inputs.append(kmermatch)
    if args.keep_norms:
        match_arguments.extend(("--norm-folder", paths.norm_folder))

    genomelift_arguments: List[object] = [
        "--input", paths.refmatch,
        "--alignments", paths.query_align,
        "--output", paths.genomelift,
        "--refhaplo", args.reference_sample,
        "--threads", args.threads,
        "--parallel-unit", "contig",
        "--tmpdir", paths.genomelift_tmp,
        "--max-extension", args.max_extension,
    ]
    if args.non_reference_tags:
        genomelift_arguments.extend((
            "--non-reference-tags", args.non_reference_tags,
        ))
    vcf_arguments: List[object] = [
        "-i", paths.graphcigartoreffix,
        "-r", fasta_reference,
        "--local-reference-templates", paths.local_templates,
        "--fasta-query", fasta_query,
        "-m", f"{paths.genomelift},{paths.genomeliftfix}",
        "-o", paths.vcf,
        "--columns", args.sample,
        "-t", args.threads,
        "--format-processes", args.format_processes,
        "--edge-blackregion", args.edge_blackregion,
        "--max-impute", args.max_extension,
    ]
    if args.seqcompress:
        vcf_arguments.append("--seqcompress")
    if args.var_in_insert > 0:
        vcf_arguments.extend(("--var-in-insert", args.var_in_insert))
    if args.svcutoff is not None:
        vcf_arguments.extend(("--svcutoff", args.svcutoff))
    vcf_arguments.append("--PAtag" if args.pa_tag else "--noPAtag")

    stages.extend((
        Stage(
            name="refmatch",
            command=script_command(
                script_dir, "match_partition_blocks.py", *match_arguments,
            ),
            inputs=tuple(match_inputs),
            outputs=(paths.refmatch,),
        ),
        Stage(
            name="genomelift",
            command=script_command(
                script_dir,
                "GenomeLift.py",
                *genomelift_arguments,
            ),
            inputs=(paths.refmatch, paths.query_align),
            outputs=(paths.genomelift,),
        ),
        Stage(
            name="graphcigartoref",
            command=script_command(
                script_dir,
                "graphcigartoref_persample.py",
                "--genomelift", paths.genomelift,
                "--align-query", paths.query_align,
                "--align-ref", paths.reference_align,
                "--query-genome", args.sample,
                "--refhaplo", args.reference_sample,
                "--fasta-query", fasta_query,
                "-G", graph_folder,
                "-L", graph_list,
                "-r", fasta_reference,
                "--local-reference-templates", paths.local_templates,
                "-o", paths.graphcigartoref,
                "-t", args.threads,
                "--max-extension", args.max_extension,
                "--buffer-size", args.buffer_size,
                "--buffer-bytes", args.buffer_bytes,
                "--progress-every", args.progress_every,
                "--progress-seconds", args.progress_seconds,
                "--maxtasksperchild", args.maxtasksperchild,
                "--cross-validation-extension", args.cross_validation_extension,
            ),
            inputs=(
                paths.genomelift, paths.query_align, paths.reference_align,
                fasta_query, fasta_query_index,
                fasta_reference, fasta_reference_index,
                graph_list, graph_folder,
                paths.local_templates,
                Path(str(paths.local_templates) + ".fai"),
            ),
            outputs=(paths.graphcigartoref,),
        ),
        Stage(
            name="graphcigartoref_gapfill",
            command=script_command(
                script_dir,
                "fill_graphcigartoref_gaps.py",
                "-i", paths.graphcigartoref,
                "--genomelift", paths.genomelift,
                "--genomelift-fix-output", paths.genomeliftfix,
                "--query-genome", args.sample,
                "--fasta-query", fasta_query,
                "--reference", fasta_reference,
                "--local-reference-templates", paths.local_templates,
                "--output", paths.graphcigartoreffix,
                "--processes", args.threads,
            ),
            inputs=(
                paths.graphcigartoref, paths.genomelift,
                fasta_query, fasta_query_index,
                fasta_reference, fasta_reference_index,
                paths.local_templates,
                Path(str(paths.local_templates) + ".fai"),
            ),
            outputs=(paths.graphcigartoreffix, paths.genomeliftfix),
        ),
        Stage(
            name="vcf",
            command=script_command(
                script_dir,
                "graphreftovcf_persample.py",
                *vcf_arguments,
            ),
            inputs=(
                paths.graphcigartoreffix, fasta_reference, fasta_query,
                script_dir / "graphreftovcf.py",
                script_dir / "graphreftovcf_persample.py",
                fasta_reference_index, fasta_query_index, paths.genomelift,
                paths.genomeliftfix,
                paths.local_templates,
                Path(str(paths.local_templates) + ".fai"),
            ),
            outputs=(paths.vcf,),
            readiness_marker=Path(str(paths.vcf) + ".complete"),
            readiness_protocol=f"insertion-snps-v1 svcutoff={args.svcutoff}\n",
        ),
    ))
    return paths, stages


def output_ready(path: Path) -> bool:
    return path.is_file() and path.stat().st_size > 0


def stage_ready(stage: Stage) -> bool:
    if not all(output_ready(path) for path in stage.outputs):
        return False
    if stage.name in ('hotspot_search', 'local_reference_templates', 'lift_local_templates', 'vcf'):
        oldest = min(path.stat().st_mtime_ns for path in stage.outputs)
        if any(not path.exists() or path.stat().st_mtime_ns > oldest for path in stage.inputs):
            return False
    if stage.readiness_marker is not None:
        try:
            marker_text = stage.readiness_marker.read_text()
        except FileNotFoundError:
            if not stage.allow_missing_readiness_marker:
                return False
        else:
            if marker_text != stage.readiness_protocol:
                return False
    if stage.output_text is None:
        return True
    if len(stage.outputs) != 1:
        return False
    try:
        return stage.outputs[0].read_text() == stage.output_text
    except OSError:
        return False


def missing_inputs(paths: Iterable[Path]) -> List[Path]:
    return [path for path in paths if not path.exists()]


def stage_run_plan(
    stages: Sequence[Stage],
    *,
    resume: bool,
    rerun: bool,
) -> List[bool]:
    """Return whether each ordered pipeline stage should be run.

    Normal mode preserves the historical behavior of independently reusing
    every existing nonempty output. Resume mode reuses only a completed prefix:
    once an incomplete stage is found, all later stages are rebuilt even if a
    later output happens to exist. This keeps downstream outputs consistent
    with the newly generated intermediate files.
    """
    if rerun:
        return [True] * len(stages)
    if not resume:
        plan = [
            not stage_ready(stage)
            for stage in stages
        ]
        # Changing hotspot targets invalidates every alignment/call derived
        # from them, including migration from the old graph-sequence index.
        for index, (stage, should_run) in enumerate(zip(stages, plan)):
            if stage.name == 'hotspot_search' and should_run:
                plan[index:] = [True] * (len(plan) - index)
                break
        return plan

    run_later = False
    plan: List[bool] = []
    for stage in stages:
        if not run_later:
            run_later = not stage_ready(stage)
        plan.append(run_later)
    return plan


def write_command_log(log, command: Sequence[str]) -> None:
    line = shlex.join(list(command))
    log.write(f"COMMAND: {line}\n")
    log.flush()
    print(f"[partition_pipeline] command: {line}", file=sys.stderr)


def run_stage(stage: Stage, log_folder: Path, rerun: bool, dry_run: bool) -> None:
    ready = stage_ready(stage)
    if ready and not rerun:
        print(
            f"[partition_pipeline] skip {stage.name}: existing output(s) "
            + ", ".join(str(path) for path in stage.outputs),
            file=sys.stderr,
        )
        return

    absent = missing_inputs(stage.inputs)
    if absent and not dry_run:
        raise FileNotFoundError(
            f"stage {stage.name} is missing input(s): "
            + ", ".join(str(path) for path in absent)
        )
    for output in stage.outputs:
        output.parent.mkdir(parents=True, exist_ok=True)

    print(f"[partition_pipeline] start {stage.name}", file=sys.stderr)
    if dry_run:
        if stage.output_text is not None:
            print(
                f"[partition_pipeline] write: {stage.outputs[0]}",
                file=sys.stderr,
            )
        else:
            print(
                f"[partition_pipeline] command: {shlex.join(list(stage.command))}",
                file=sys.stderr,
            )
        return

    log_folder.mkdir(parents=True, exist_ok=True)
    log_path = log_folder / f"{stage.name}.log"
    if stage.output_text is not None:
        if len(stage.outputs) != 1:
            raise ValueError(
                f"internal text stage {stage.name} must have exactly one output"
            )
        output = stage.outputs[0]
        descriptor, temporary_name = tempfile.mkstemp(
            prefix=f".{output.name}.",
            suffix=".tmp",
            dir=output.parent,
            text=True,
        )
        try:
            with os.fdopen(descriptor, "wt") as temporary:
                temporary.write(stage.output_text)
            os.replace(temporary_name, output)
        except BaseException:
            try:
                os.unlink(temporary_name)
            except FileNotFoundError:
                pass
            raise
        with open(log_path, "wt") as log:
            log.write(f"WROTE: {output}\n")
        print(f"[partition_pipeline] finished {stage.name}", file=sys.stderr)
        return

    if stage.write_readiness_marker and stage.readiness_marker is not None:
        stage.readiness_marker.unlink(missing_ok=True)
    with open(log_path, "wt") as log:
        write_command_log(log, stage.command)
        process = subprocess.Popen(
            stage.command,
            stdout=subprocess.PIPE,
            stderr=subprocess.STDOUT,
            text=True,
            bufsize=1,
        )
        try:
            assert process.stdout is not None
            for line in process.stdout:
                sys.stderr.write(line)
                log.write(line)
            return_code = process.wait()
        except BaseException:
            process.terminate()
            process.wait()
            raise
        finally:
            if process.stdout is not None:
                process.stdout.close()
    if return_code:
        # A streaming stage can leave a nonempty partial TSV after SIGTERM or
        # a worker failure.  Do not let that file satisfy ``--resume`` and feed
        # incomplete data to a downstream stage.
        removed = []
        for output in stage.outputs:
            try:
                output.unlink()
                removed.append(str(output))
            except FileNotFoundError:
                pass
        if stage.readiness_marker is not None:
            try:
                stage.readiness_marker.unlink()
                removed.append(str(stage.readiness_marker))
            except FileNotFoundError:
                pass
        if removed:
            print(
                f"[partition_pipeline] removed partial output(s) after failed "
                f"stage {stage.name}: {', '.join(removed)}",
                file=sys.stderr,
            )
        raise RuntimeError(
            f"stage {stage.name} failed with exit code {return_code}; "
            f"see {log_path}"
        )
    not_created = [path for path in stage.outputs if not output_ready(path)]
    if not not_created and stage.write_readiness_marker and stage.readiness_marker is not None:
        stage.readiness_marker.write_text(stage.readiness_protocol or '')
    if not not_created and stage.readiness_marker is not None:
        try:
            marker_text = stage.readiness_marker.read_text()
        except FileNotFoundError:
            marker_text = None
        if marker_text != stage.readiness_protocol:
            raise RuntimeError(
                f"stage {stage.name} completed but did not write expected "
                f"mode marker {stage.readiness_marker}"
            )
    if not_created:
        raise RuntimeError(
            f"stage {stage.name} completed but did not create nonempty output(s): "
            + ", ".join(str(path) for path in not_created)
        )
    print(f"[partition_pipeline] finished {stage.name}", file=sys.stderr)


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description=(
            "Run one sample through partition alignment, block matching, "
            "GenomeLift, graph-CIGAR conversion, and VCF generation."
        )
    )
    parser.add_argument("--sample", required=True, help="query sample/haplotype, e.g. HG00513_h1")
    parser.add_argument(
        "--fasta-query",
        required=True,
        help="query FASTA; masking, contig naming, and faidx are prepared as needed",
    )
    parser.add_argument(
        "--hotspot-query",
        help=(
            "existing query sample _hotspot.txt; when omitted, KmerSearcher "
            "writes <output-root>/<sample>/<sample>_hotspot.txt"
        ),
    )
    parser.add_argument("--reference-sample", default="CHM13_h1")
    parser.add_argument(
        "--no-local-template-fallback",
        dest="local_template_fallback",
        action="store_false",
        help=(
            "disable the default behavior that calls reference-free local "
            "graphs against their locus template"
        ),
    )
    parser.set_defaults(local_template_fallback=True)
    parser.add_argument(
        "--fasta-reference",
        required=True,
        help="reference FASTA; a missing .fai is created automatically",
    )
    parser.add_argument(
        "--reference-align",
        required=True,
        help="precomputed reference align_partition_hotspots.py output",
    )
    parser.add_argument(
        "--reference-blocks",
        required=True,
        help="precomputed reference block_partition_alignments.py BED6 output",
    )
    parser.add_argument("-G", "--graph-folder", required=True)
    parser.add_argument(
        "-L", "--graph-list",
        help="ordered graph alignment list; default: <graph-folder>.list",
    )
    parser.add_argument("-O", "--output-root", default="LocalGraphsAligns")
    parser.add_argument("-t", "--threads", type=int, default=16)
    parser.add_argument(
        "--align-threads",
        type=int,
        default=0,
        help="alignment CPU budget; default: use the full -t/--threads value",
    )
    parser.add_argument("--threads-per-job", type=int, default=4)
    alignment_mode = parser.add_mutually_exclusive_group()
    alignment_mode.add_argument(
        "--fast",
        action="store_true",
        help=(
            "use cohort-style fast local-graph alignment: bulk BLASTN for "
            "each graph and Winnowmap only when the BLASTN result is not "
            "sufficiently complete"
        ),
    )
    alignment_mode.add_argument(
        "--rigorous", "--accurate",
        dest="fast",
        action="store_false",
        help=(
            "run BLASTN and Winnowmap for every query hotspot (default)"
        ),
    )
    parser.set_defaults(fast=False)
    parser.add_argument("--summary-jobs", type=int, default=0, help="default: min(-t, 16)")
    parser.add_argument("--block-jobs", type=int, default=0, help="default: -t")
    parser.add_argument("--match-jobs", type=int, default=0, help="default: -t")
    parser.add_argument("--format-processes", type=int, default=16)
    parser.add_argument(
        "--svonly", "--svcutoff", dest="svcutoff", nargs="?", const=20,
        default=None, type=int, metavar="BP",
        help="exclude SNPs (including inside insertions) and retain SVs larger than BP [20]",
    )
    parser.add_argument(
        "--seqcompress",
        action="store_true",
        help=(
            "store query sequence in EXTENDGRAPHCIGAR payloads instead of "
            "plain INFO/SEQ"
        ),
    )
    parser.add_argument(
        "--var-in-insert",
        nargs="?",
        const=100,
        default=0,
        type=int,
        metavar="BP",
        help=(
            "recursively compare merged insertion observations with minimap2 "
            "when the largest has more than BP uppercase bases; a bare option "
            "uses 100 (default: disabled)"
        ),
    )
    pa_group = parser.add_mutually_exclusive_group()
    pa_group.add_argument(
        "--PAtag", "--pa-tag",
        dest="pa_tag",
        action="store_true",
        default=True,
        help=(
            "retain final VCF ALLELENAME provenance (default); LABEL_H is "
            "always retained"
        ),
    )
    pa_group.add_argument(
        "--noPAtag", "--no-pa-tag",
        dest="pa_tag",
        action="store_false",
        help=(
            "remove final VCF ALLELENAME while retaining LABEL_H"
        ),
    )
    parser.add_argument(
        "--buffer-size",
        type=int,
        default=1000,
        metavar="ROWS",
        help=(
            "completed graphcigartoref rows buffered before writing "
            "(default: 1000); use 0 to flush every completed parallel job "
            "for diagnosis"
        ),
    )
    parser.add_argument(
        "--buffer-bytes",
        type=int,
        default=64 * 1024 * 1024,
        metavar="BYTES",
        help=(
            "hard graphcigartoref result-buffer byte cap (default: 67108864; "
            "0 disables the cap)"
        ),
    )
    parser.add_argument(
        "--progress-every",
        type=int,
        default=1000,
        metavar="N",
        help="graphcigartoref progress interval; 0 disables (default: 1000)",
    )
    parser.add_argument(
        "--progress-seconds",
        type=float,
        default=60.0,
        metavar="SECONDS",
        help="graphcigartoref waiting diagnostic interval; 0 disables (default: 60)",
    )
    parser.add_argument(
        "--maxtasksperchild",
        type=int,
        default=512,
        metavar="N",
        help="recycle graphcigartoref workers after N comparisons (default: 512; 0 disables)",
    )
    parser.add_argument(
        "--edge-blackregion",
        type=int,
        default=0,
        help=(
            "exclude VCF events within this many bp of a reference whole-"
            "alignment edge (default: 0, disabled; set e.g. 500 to enable)"
        ),
    )
    parser.add_argument("--extension", type=int, default=30000)
    payload_group = parser.add_mutually_exclusive_group()
    payload_group.add_argument(
        "--alignment-payload", dest="alignment_payload", action="store_true",
        default=False, help="store I/X sequence payloads in intermediate alignment files",
    )
    payload_group.add_argument(
        "--no-alignment-payload", dest="alignment_payload", action="store_false",
        help="save compact alignments and restore payloads per comparison (default)",
    )
    parser.add_argument('--template-assemblies', default='', help='NAME FASTA [FAI] sources for lifting fallback templates')
    parser.add_argument(
        "--max-extension",
        type=int,
        default=10000,
        help=(
            "maximum GenomeLift ownership extension used for graph-CIGAR "
            "slices and VCF labeling (default: 10000)"
        ),
    )
    parser.add_argument(
        "--cross-validation-extension",
        type=int,
        default=4000,
        metavar="BP",
        help=(
            "symmetric graphcigartoref context for cross-graph validation "
            "(default: 4000; use 0 for a diagnostic compatibility run)"
        ),
    )
    parser.add_argument(
        "--non-reference-tags",
        default="",
        help=(
            "comma-separated exact non-reference allele names allowed to "
            "seed GenomeLift tags, or 'all' (default: none)"
        ),
    )
    parser.add_argument("--merge-gap", type=int, default=30000)
    parser.add_argument("--timeout", type=int, default=3600)
    parser.add_argument(
        "--windowmasker",
        default="windowmasker",
        help="WindowMasker executable used when the query FASTA has no lowercase a/c/g/t",
    )
    parser.add_argument(
        "--samtools",
        default="samtools",
        help="samtools executable used to create missing FASTA indexes",
    )
    parser.add_argument(
        "--kmer-searcher",
        default="KmerSearcher",
        help="KmerSearcher executable used when --hotspot-query is omitted",
    )
    parser.add_argument(
        "--hotspot-cutoff",
        type=int,
        default=1000,
        help="KmerSearcher -c cutoff used for generated hotspot tables",
    )
    parser.add_argument(
        "--hotspot-window",
        type=int,
        default=2000,
        help="KmerSearcher -w window used for generated hotspot tables",
    )
    parser.add_argument("--kmermatch", help="compiled KmerMatch executable; default: beside match_partition_blocks.py")
    parser.add_argument("--keep-norms", action="store_true", help="retain per-graph KmerMatch/RefMatch matrices")
    run_mode = parser.add_mutually_exclusive_group()
    run_mode.add_argument(
        "--resume",
        action="store_true",
        help=(
            "skip the consecutive completed stages at the start, then rerun "
            "the first incomplete stage and all downstream stages"
        ),
    )
    run_mode.add_argument(
        "--rerun",
        action="store_true",
        help="rerun every stage even when its output already exists",
    )
    parser.add_argument("--dry-run", action="store_true", help="print pending commands without running them")
    return parser


def parse_args(argv: Optional[Sequence[str]] = None) -> argparse.Namespace:
    parser = build_parser()
    args = parser.parse_args(argv)
    for name in (
        "threads", "threads_per_job", "format_processes",
    ):
        if getattr(args, name) < 1:
            parser.error(f"--{name.replace('_', '-')} must be at least 1")
    for name in ("align_threads", "summary_jobs", "block_jobs", "match_jobs"):
        if getattr(args, name) < 0:
            parser.error(f"--{name.replace('_', '-')} cannot be negative")
    if args.edge_blackregion < 0:
        parser.error("--edge-blackregion cannot be negative")
    if args.var_in_insert < 0:
        parser.error("--var-in-insert cannot be negative")
    if args.svcutoff is not None and args.svcutoff < 0:
        parser.error("--svcutoff cannot be negative")
    if args.max_extension < 0:
        parser.error("--max-extension cannot be negative")
    if args.buffer_size < 0:
        parser.error("--buffer-size cannot be negative")
    if args.buffer_bytes < 0:
        parser.error("--buffer-bytes cannot be negative")
    if args.progress_every < 0:
        parser.error("--progress-every cannot be negative")
    if args.progress_seconds < 0:
        parser.error("--progress-seconds cannot be negative")
    if args.maxtasksperchild < 0:
        parser.error("--maxtasksperchild cannot be negative")
    if args.cross_validation_extension < 0:
        parser.error("--cross-validation-extension cannot be negative")
    if args.hotspot_cutoff < 1:
        parser.error("--hotspot-cutoff must be at least 1")
    if args.hotspot_window < 1:
        parser.error("--hotspot-window must be at least 1")
    return args


def preflight(
    args: argparse.Namespace,
    script_dir: Path,
    stages: Sequence[Stage],
    run_plan: Optional[Sequence[bool]] = None,
) -> None:
    missing_scripts = [
        pipeline_script_path(script_dir, name)
        for name in (
            "align_partition_hotspots.py",
            "summarize_partition_hotspot_segments.py",
            "block_partition_alignments.py",
            "match_partition_blocks.py",
            "GenomeLift.py",
            "local_reference_templates.py",
            "lift_local_templates.py",
            "assembly_contigs.py",
            "graphcigartoref_persample.py",
            "fill_graphcigartoref_gaps.py",
            "graphreftovcf_persample.py",
            "namecontigsfix.py",
        )
        if not pipeline_script_path(script_dir, name).is_file()
    ]
    if missing_scripts:
        raise FileNotFoundError(
            "missing pipeline script(s): "
            + ", ".join(str(path) for path in missing_scripts)
        )

    if not args.kmermatch:
        executable = script_dir / "KmerMatch"
        if not executable.is_file() or not os.access(executable, os.X_OK):
            raise FileNotFoundError(
                f"compiled KmerMatch not found/executable at {executable}; "
                f"run: make -C {shlex.quote(str(script_dir))} KmerMatch"
            )

    # Check only roots required by the first pending stage; later generated
    # inputs are checked immediately before their consuming stage.
    if run_plan is None:
        run_plan = stage_run_plan(
            stages,
            resume=args.resume,
            rerun=args.rerun,
        )
    generated_outputs = set()
    for stage, should_run in zip(stages, run_plan):
        if should_run:
            if stage.output_text is None:
                executable = stage.command[0]
                if os.path.sep in executable:
                    executable_found = (
                        Path(executable).is_file()
                        and os.access(executable, os.X_OK)
                    )
                else:
                    executable_found = shutil.which(executable) is not None
                if not executable_found:
                    raise FileNotFoundError(
                        f"stage {stage.name} requires executable {executable!r}"
                    )
            absent = missing_inputs(stage.inputs)
            absent = [path for path in absent if path not in generated_outputs]
            if absent and not args.dry_run:
                raise FileNotFoundError(
                    f"stage {stage.name} is missing source input(s): "
                    + ", ".join(str(path) for path in absent)
                )
        generated_outputs.update(stage.outputs)


def main(argv: Optional[Sequence[str]] = None) -> int:
    args = parse_args(argv)
    script_dir = Path(__file__).resolve().parent
    args.kmer_searcher = prefer_local_executable(
        args.kmer_searcher, script_dir,
    )
    graph_folder = clean_path(args.graph_folder, "--graph-folder")
    canonical_graph_list = graph_list_path(graph_folder)
    requested_graph_list = (
        clean_path(args.graph_list, "--graph-list")
        if args.graph_list else canonical_graph_list
    )
    if requested_graph_list == canonical_graph_list:
        if args.dry_run:
            args.graph_list = str(canonical_graph_list)
            print(
                f"[partition_pipeline] would ensure graph list/template cache: "
                f"{canonical_graph_list}, {template_list_path(canonical_graph_list)}.bin",
                file=sys.stderr,
            )
        else:
            listing, cache, rebuilt = ensure_graph_list_cache(
                graph_folder,
                jobs=args.threads,
                script_dir=script_dir,
            )
            args.graph_list = str(listing)
            args.template_list = str(cache.with_suffix(''))
            if rebuilt:
                print(
                    f"[partition_pipeline] regenerated graph list/template cache: "
                    f"{listing}, {cache}",
                    file=sys.stderr,
                )
    else:
        args.graph_list = str(requested_graph_list)
        if not args.dry_run and not args.hotspot_query:
            templates, _cache, _rebuilt = ensure_template_list_cache(
                graph_folder, requested_graph_list, jobs=args.threads,
                script_dir=script_dir,
            )
            args.template_list = str(templates)
    prepared_query, preparation_stages = build_fasta_preparation_stages(
        args, script_dir,
    )
    args.fasta_query = str(prepared_query)
    hotspot_query, hotspot_stages = build_hotspot_stages(
        args, prepared_query,
    )
    args.hotspot_query = str(hotspot_query)
    paths, sample_stages = build_stages(args, script_dir)
    stages = [
        *preparation_stages,
        *hotspot_stages,
        build_reference_index_stage(args),
        *sample_stages,
    ]
    run_plan = stage_run_plan(
        stages,
        resume=args.resume,
        rerun=args.rerun,
    )
    if args.resume:
        first_pending = next(
            (stage.name for stage, should_run in zip(stages, run_plan) if should_run),
            None,
        )
        if first_pending is None:
            print(
                "[partition_pipeline] resume: all stages are already complete",
                file=sys.stderr,
            )
        else:
            print(
                f"[partition_pipeline] resume from {first_pending}; "
                "all downstream stages will be rebuilt",
                file=sys.stderr,
            )
    preflight(args, script_dir, stages, run_plan)
    paths.log_folder.mkdir(parents=True, exist_ok=True)
    for stage, should_run in zip(stages, run_plan):
        run_stage(stage, paths.log_folder, should_run, args.dry_run)
    print(
        f"[partition_pipeline] complete: {paths.vcf}",
        file=sys.stderr,
    )
    return 0


if __name__ == "__main__":
    try:
        raise SystemExit(main())
    except KeyboardInterrupt:
        print("[partition_pipeline] interrupted", file=sys.stderr)
        raise SystemExit(130)
    except Exception as error:
        print(f"run_sample_pipeline.py: error: {error}", file=sys.stderr)
        raise SystemExit(1)
