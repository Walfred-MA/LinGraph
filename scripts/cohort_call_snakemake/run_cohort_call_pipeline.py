#!/usr/bin/env python3
"""Launch the checkpointed cohort graph-to-VCF workflow."""

from __future__ import annotations

import argparse
from collections import Counter
import hashlib
import json
import os
from pathlib import Path
import shlex
import shutil
import subprocess
import sys
from typing import List, Optional, Sequence


REQUIRED_SNAKEMAKE = "6.15.1"
HERE = Path(__file__).resolve().parent
REPOSITORY = HERE.parent
GRAPH_INPUT_HELPERS = REPOSITORY / "graph_build_snakemake" / "workflow" / "scripts"
sys.path.insert(0, str(GRAPH_INPUT_HELPERS))
sys.path.insert(0, str(REPOSITORY))

from align_partition_hotspots import graph_name_from_list_entry
from graph_list_cache import graph_list_path
from pipeline_inputs import read_assembly_list
import graphvcfmerge_checkpoints as merge_checkpoints
import cohort_vcf_merge as cohort_merge


def absolute(value: str, base: Optional[Path] = None) -> str:
    path = Path(value).expanduser()
    if not path.is_absolute():
        path = (base or Path.cwd()) / path
    return str(path.resolve())


def read_sample_selection(value: str) -> List[str]:
    if not value:
        return []
    candidate = Path(value).expanduser()
    if candidate.is_file():
        output = []
        with candidate.open("rt", encoding="utf-8") as handle:
            for raw in handle:
                fields = shlex.split(raw, comments=True)
                if fields:
                    # Accept either a names-only file or the same
                    # NAME FASTA format used by query_paths.txt.
                    output.extend(
                        token.strip() for token in fields[0].split(",")
                        if token.strip()
                    )
        return output
    return [token.strip() for token in value.split(",") if token.strip()]


def read_partition_names(path: str) -> List[str]:
    names: List[str] = []
    seen = set()
    with open(path, "rt", encoding="utf-8") as handle:
        for line_number, raw in enumerate(handle, 1):
            if not raw.strip() or raw.lstrip().startswith("#"):
                continue
            name = graph_name_from_list_entry(
                raw.strip().split()[0], f"{path}:{line_number}",
            )
            if name in seen:
                raise ValueError(
                    f"{path}:{line_number}: duplicate partition {name!r}"
                )
            seen.add(name)
            names.append(name)
    if not names:
        raise ValueError(f"partition list is empty: {path}")
    return names


def write_text_if_changed(path: Path, text: str) -> bool:
    try:
        if path.read_text(encoding="utf-8") == text:
            return False
    except FileNotFoundError:
        pass
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_name(path.name + f".tmp.{os.getpid()}")
    try:
        temporary.write_text(text, encoding="utf-8")
        os.replace(temporary, path)
    finally:
        try:
            temporary.unlink()
        except FileNotFoundError:
            pass
    return True


def write_json_if_changed(path: Path, value: dict) -> bool:
    text = json.dumps(value, indent=2, sort_keys=True) + "\n"
    return write_text_if_changed(path, text)


def merge_needs_resume(root: Path, context: dict) -> bool:
    """Do not let an old aggregate marker hide an unfinished chromosome."""
    chrom_list = root / "chroms.txt"
    if not chrom_list.is_file():
        return False
    with chrom_list.open() as handle:
        chroms = [line.split("\t", 1)[0].strip() for line in handle if line.strip()]
    return any(not merge_checkpoints.is_done(root, chrom, context, adopt_manual=False)
               for chrom in chroms)


def write_stage_config(
    path: Path,
    value: dict,
    migration_markers: Optional[Sequence[Path]] = None,
) -> None:
    """Write a semantic stage config without invalidating legacy results.

    A newly introduced config can be older than existing stage outputs:
    this migrates an already prepared run without forcing a one-time rerun.
    Subsequent content changes retain their real timestamp and invalidate only
    the stage that declares this file as an input.
    """
    existed = path.is_file()
    changed = write_json_if_changed(path, value)
    existing_markers = [
        marker for marker in (migration_markers or ()) if marker.is_file()
    ]
    if changed and not existed and existing_markers:
        # Leave a full second of separation for shared filesystems whose
        # timestamp resolution is coarser than nanoseconds.
        inherited = max(
            1,
            min(marker.stat().st_mtime_ns for marker in existing_markers)
            - 1_000_000_000,
        )
        os.utime(path, ns=(inherited, inherited))


def check_snakemake(command: str, skip: bool) -> None:
    executable = shutil.which(command) if os.path.sep not in command else command
    if not executable or not os.path.isfile(executable):
        raise FileNotFoundError(
            f"Snakemake {REQUIRED_SNAKEMAKE} is required but {command!r} "
            "was not found"
        )
    if skip:
        return
    version = subprocess.check_output([command, "--version"], text=True).strip()
    if version != REQUIRED_SNAKEMAKE:
        raise RuntimeError(
            f"this workflow requires Snakemake {REQUIRED_SNAKEMAKE}; "
            f"{command!r} reports {version!r}"
        )


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description=(
            "Extract selected graph-build cohort alignments, lift every "
            "haplotype to an included reference, write one VCF per "
            "haplotype, and finish with a disk-sharded graphvcfmerge.py "
            "multi-sample SV callset."
        )
    )
    parser.add_argument("-r", "--reference", default="",
                        help="reference sample/haplotype in the graph cohort (default: first assembly)")
    parser.add_argument("-G", "--graph-folder", default="",
                        help="read-only completed graph-build root")
    parser.add_argument('--alternative', default='',
                        help='original-template BED (default: GRAPH/summary/alternative_intervals.bed when present)')
    parser.add_argument("-O", "--output-folder", default="",
                        help="separate cohort-calling output root")
    parser.add_argument(
        "--unlock", action="store_true",
        help=(
            "unlock the cohort Snakemake working directory using the existing "
            "-O/--output-folder run configuration, then exit; only -O is "
            "required in this mode"
        ),
    )
    parser.add_argument(
        "-q", "--assemblies", default="",
        help=(
            "full graph-cohort NAME FASTA list; default: "
            "GRAPH/inputs/query_paths.normalized.txt"
        ),
    )
    parser.add_argument(
        "-L", "--partition-list", default="",
        help="active graph partition list; default: GRAPH/Graphs.list",
    )
    parser.add_argument(
        "--samples", default="",
        help=(
            "optional comma-separated names, names-only file, or NAME "
            "FASTA query-paths file to call; the reference is always "
            "included"
        ),
    )
    parser.add_argument("-j", "--cores", type=int, default=64)
    parser.add_argument("--sample-threads", type=int, default=16)
    parser.add_argument(
        "--match-batches", type=int, default=20,
        help="number of partition-match batch jobs in the DAG (default: 20)",
    )
    parser.add_argument(
        "--match-threads", type=int, default=4,
        help="parallel local graphs processed inside each match batch (default: 4)",
    )
    parser.add_argument(
        "--io-jobs", type=int, default=16,
        help=(
            "worker processes inside the single extraction/concatenation "
            "job (default: 16)"
        ),
    )
    parser.add_argument(
        "--extract-batch-size", type=int, default=10_000,
        help=(
            "maximum partition files read per internal extraction batch "
            "(default: 10000)"
        ),
    )
    parser.add_argument("--format-processes", type=int, default=16)
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
        default=100,
        type=int,
        metavar="BP",
        help=(
            "recursively promote supported insertions within merged insertion "
            "templates at least BP long; a bare option and the default use "
            "100, while 0 disables recursion"
        ),
    )
    pa_group = parser.add_mutually_exclusive_group()
    pa_group.add_argument(
        "--PAtag", "--pa-tag",
        dest="pa_tag",
        action="store_true",
        default=True,
        help=(
            "retain final VCF ALLELENAME and LABEL_H provenance (default)"
        ),
    )
    pa_group.add_argument(
        "--noPAtag", "--no-pa-tag",
        dest="pa_tag",
        action="store_false",
        help="remove final VCF ALLELENAME and LABEL_H provenance",
    )
    conflict_group = parser.add_mutually_exclusive_group()
    conflict_group.add_argument(
        "--resolve-conflicts-by-alignment-score",
        dest="resolve_conflicts_by_alignment_score",
        action="store_true",
        default=True,
        help=(
            "retain an otherwise rejected cross-graph SV only when its PA "
            "has the uniquely best main-path affine alignment score (default)"
        ),
    )
    conflict_group.add_argument(
        "--no-resolve-conflicts-by-alignment-score",
        dest="resolve_conflicts_by_alignment_score",
        action="store_false",
        help="use the legacy conflict-removal rule",
    )
    parser.add_argument(
        "--separate-adjacent-indels",
        "--separate-adjacent-indel-components",
        dest="separate_adjacent_indels",
        action="store_true",
        help=(
            "emit directly adjacent insertion and deletion components as "
            "separate per-sample VCF rows"
        ),
    )
    mode_group = cohort_merge.add_modes(parser)
    parser.add_argument('--mc-graph', '--MC-graph', dest='mc_graph', action='store_true',
                        help='select --all by default for subsequent GFA export')
    stages = parser.add_mutually_exclusive_group()
    stages.add_argument('--merge-only', action='store_true', help='merge existing sample VCFs without upstream calling')
    stages.add_argument('--recall-only', action='store_true', help='recall VCFs and merge from existing saved inputs; do not schedule upstream jobs')
    mode_group.add_argument(
        "--SVonly", "--sv-only",
        dest="sv_only_size",
        nargs="?",
        const=20,
        default=None,
        type=int,
        metavar="SIZE",
        help=(
            "write SV/indel-only VCFs; SIZE defaults to 20 bp when omitted, "
            "and clusters with MAXSIZE > 0.6*SIZE are retained"
        ),
    )
    parser.add_argument(
        "--svcutoff", "--minsvsize", dest="minsvsize", type=int, default=20, metavar="BP",
        help=(
            "minimum representative size in the final merged cohort SV VCF "
            "when --SVonly is absent (default: 20); smaller indels go to "
            "cohort.indel.vcf with --svindel; --SVonly SIZE overrides "
            "this final cutoff"
        ),
    )
    parser.add_argument(
        "--merge-distance", type=int, default=500, metavar="BP",
        help="maximum reference distance for final cohort SV merging (default: 500)",
    )
    parser.add_argument(
        "--size-similarity", type=float, default=0.7, metavar="FRACTION",
        help="minimum min(size)/max(size) for final SV merging (default: 0.7)",
    )
    parser.add_argument(
        "--sequence-similarity", type=float, default=0.7,
        metavar="FRACTION",
        help="compatibility option; final insertion merging uses fixed KmerMatch cosine >=0.5",
    )
    parser.add_argument(
        "--merge-processes", type=int, default=None, metavar="INT",
        help=(
            "parallel per-sample VCF readers in graphvcfmerge.py "
            "(default: min(16, --cores))"
        ),
    )
    parser.add_argument(
        "--merge-locus-processes", type=int, default=None, metavar="INT",
        help=(
            "deprecated compatibility option; the current graphvcfmerge.py "
            "uses --merge-processes for all worker pools"
        ),
    )
    parser.add_argument(
        "--merge-chr1-resource-multiplier", type=int, default=2, metavar="INT",
        help=(
            "Slurm chr1 merge CPU/total-memory multiplier in a dedicated job "
            "(default: 2); 1 disables; includes CHM13/GRCh38 RefSeq chr1 names"
        ),
    )
    parser.add_argument(
        "--merge-tmpdir", default="", metavar="DIR",
        help=(
            "graphvcfmerge.py shard directory (default: "
            "OUTPUT/tmp/graphvcfmerge)"
        ),
    )
    parser.add_argument(
        "--keep-merge-tmpdir", action="store_true",
        help=(
            "retain merge intermediates after successful VCF publication "
            "(local and Slurm); failed merges retain their scratch directories"
        ),
    )
    parser.add_argument("--max-extension", type=int, default=10_000)
    parser.add_argument("--cross-validation-extension", type=int, default=4_000)
    realignment_group = parser.add_mutually_exclusive_group()
    realignment_group.add_argument(
        "--realignment",
        dest="realignment",
        action="store_true",
        default=True,
        help=(
            "enable optional local and score-linked SV-region realignment "
            "during graph-CIGAR and VCF generation (default)"
        ),
    )
    realignment_group.add_argument(
        "--no-realignment",
        dest="realignment",
        action="store_false",
        help="disable graph-CIGAR and VCF realignment",
    )
    parser.add_argument("--edge-blackregion", type=int, default=0)
    parser.add_argument("--buffer-size", type=int, default=1000)
    parser.add_argument("--buffer-bytes", type=int, default=64 * 1024 * 1024)
    parser.add_argument("--progress-every", type=int, default=1000)
    parser.add_argument("--progress-seconds", type=float, default=60.0)
    parser.add_argument("--maxtasksperchild", type=int, default=512)
    parser.add_argument("--non-reference-tags", default="")
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
    parser.add_argument("--script-folder", default=str(REPOSITORY))
    parser.add_argument("--kmermatch", default=str(REPOSITORY / "KmerMatch"))
    parser.add_argument("--snakemake", default="snakemake")
    parser.add_argument("--skip-version-check", action="store_true")
    parser.add_argument("--dry-run", "-n", action="store_true")
    parser.add_argument("--slurm", action="store_true")
    parser.add_argument("--slurm-jobs", type=int, default=20)
    parser.add_argument("--slurm-account", default="")
    parser.add_argument("--slurm-time", default="200:00:00")
    parser.add_argument("--slurm-partition", default="")
    parser.add_argument("--slurm-args", default="")
    parser.add_argument("--match-memory", default="32G")
    parser.add_argument("--extraction-memory", default="64G")
    parser.add_argument("--genomelift-memory", default="32G")
    parser.add_argument("--graphcigar-memory", default="32G")
    parser.add_argument("--gapfill-memory", default="32G")
    parser.add_argument("--vcf-memory", default="32G")
    parser.add_argument(
        "--merge-memory", default="", metavar="MEMORY",
        help=(
            "base memory for merge-stage Slurm jobs; by default each stage "
            "requests 2G per allocated CPU; chr1 applies its resource multiplier"
        ),
    )
    parser.add_argument(
        "snakemake_args", nargs=argparse.REMAINDER,
        help="extra Snakemake arguments after --",
    )
    return parser


def parse_args(argv: Optional[Sequence[str]] = None) -> argparse.Namespace:
    parser = build_parser()
    args = parser.parse_args(argv)
    if args.unlock:
        if not args.output_folder:
            parser.error("--unlock requires -O/--output-folder")
        return args
    missing = [
        option for option, value in (
            ("-G/--graph-folder", args.graph_folder),
            ("-O/--output-folder", args.output_folder),
        ) if not value
    ]
    if missing:
        parser.error("the following arguments are required: " + ", ".join(missing))
    if args.merge_processes is None:
        args.merge_processes = min(16, args.cores)
    if args.merge_locus_processes is None:
        args.merge_locus_processes = min(4, args.cores)
    for name in (
        "cores", "sample_threads", "match_batches", "match_threads",
        "io_jobs", "extract_batch_size", "format_processes", "slurm_jobs",
        "merge_processes", "merge_locus_processes", "merge_chr1_resource_multiplier",
    ):
        if getattr(args, name) < 1:
            parser.error(f"--{name.replace('_', '-')} must be positive")
    for name in (
        "max_extension", "cross_validation_extension", "edge_blackregion",
        "var_in_insert",
        "buffer_size", "buffer_bytes", "progress_every", "progress_seconds",
        "maxtasksperchild",
    ):
        if getattr(args, name) < 0:
            parser.error(f"--{name.replace('_', '-')} cannot be negative")
    if args.sv_only_size is not None and args.sv_only_size < 1:
        parser.error("--SVonly SIZE must be positive")
    if args.sv_only_size is not None:
        args.merge_mode = 'svonly'
    args.merge_mode = cohort_merge.resolve_mode(args)
    if args.minsvsize < 1:
        parser.error("--minsvsize must be positive")
    if args.merge_distance < 0:
        parser.error("--merge-distance cannot be negative")
    if not 0 <= args.size_similarity <= 1:
        parser.error("--size-similarity must be between 0 and 1")
    if not 0 <= args.sequence_similarity <= 1:
        parser.error("--sequence-similarity must be between 0 and 1")
    if args.seqcompress:
        parser.error(
            "--seqcompress is incompatible with cohort merging because the "
            "merger requires a plain INFO/SEQ allele for every insertion"
        )
    try:
        args.slurm_args = shlex.join(shlex.split(args.slurm_args))
    except ValueError as error:
        parser.error(f"invalid --slurm-args: {error}")
    if args.slurm and not args.slurm_args and not args.slurm_account:
        parser.error("--slurm-account is required with --slurm unless --slurm-args supplies it")
    if args.slurm and not args.slurm_args and not args.slurm_partition:
        parser.error("--slurm-partition is required with --slurm unless --slurm-args supplies it")
    if args.slurm_args and not args.slurm:
        parser.error("--slurm-args requires --slurm")
    if not args.slurm and not args.recall_only and args.sample_threads > args.cores:
        parser.error("--sample-threads cannot exceed --cores without --slurm")
    if not args.slurm and not args.recall_only and args.match_threads > args.cores:
        parser.error("--match-threads cannot exceed --cores without --slurm")
    if not args.slurm and args.merge_processes > args.cores:
        parser.error("--merge-processes cannot exceed --cores without --slurm")
    if not args.slurm and args.merge_locus_processes > args.cores:
        parser.error("--merge-locus-processes cannot exceed --cores without --slurm")
    return args


def main(argv: Optional[Sequence[str]] = None) -> int:
    args = parse_args(argv)
    if args.recall_only and not args.unlock:
        if args.slurm:
            raise ValueError('--recall-only uses the current allocation; use LinGraph --slurm or an allocated shell')
        if args.assemblies or args.snakemake_args:
            raise ValueError('--recall-only uses saved assembly paths and no Snakemake arguments; use --samples for selection')
        from cohort_vcf_recall import run as recall
        recall(absolute(args.output_folder), args.merge_mode, graph=absolute(args.graph_folder),
               reference=args.reference, samples=read_sample_selection(args.samples) if args.samples else None,
               processes=args.cores, cutoff=args.sv_only_size or args.minsvsize, dry_run=args.dry_run)
        return 0
    check_snakemake(args.snakemake, args.skip_version_check)

    if args.unlock:
        output = Path(absolute(args.output_folder))
        run_config = output / "inputs" / "cohort_call.run.json"
        if not run_config.is_file():
            raise FileNotFoundError(
                f"cannot unlock from {output}: existing cohort run "
                f"configuration is missing: {run_config}"
            )
        command = [
            args.snakemake,
            "--snakefile", str(HERE / "Snakefile"),
            "--unlock",
            "--config", f"run_config={run_config}",
        ]
        return subprocess.run(command, cwd=str(HERE)).returncode

    if args.merge_only:
        if args.slurm:
            raise ValueError('--merge-only runs inside the current allocation; use LinGraph --slurm or an allocated shell')
        cohort_merge.run(absolute(args.output_folder), args.merge_mode, processes=args.merge_processes,
                         cutoff=args.sv_only_size or args.minsvsize, merge_distance=args.merge_distance,
                         size_similarity=args.size_similarity, sequence_similarity=args.sequence_similarity,
                         var_in_insert=args.var_in_insert, kmermatch=args.kmermatch, dry_run=args.dry_run,
                         keep_merge_tmpdir=args.keep_merge_tmpdir)
        return 0

    graph = Path(absolute(args.graph_folder))
    alternative_candidate = Path(absolute(args.alternative)) if args.alternative else graph / 'summary' / 'alternative_intervals.bed'
    if args.alternative and not alternative_candidate.is_file():
        raise FileNotFoundError(alternative_candidate)
    alternative_bed = str(alternative_candidate) if alternative_candidate.is_file() else ''
    if alternative_bed:
        from alternative_intervals import AlternativeIntervals
        alternative_digest = AlternativeIntervals(alternative_bed).digest
    else:
        alternative_digest = ''
    output = Path(absolute(args.output_folder))
    if graph == output:
        raise ValueError("--output-folder must differ from --graph-folder")
    graph_partitions = graph / "Graphs"
    if not graph_partitions.is_dir():
        raise NotADirectoryError(f"completed graph partition root not found: {graph_partitions}")
    script_folder = Path(absolute(args.script_folder))

    query_paths = absolute(
        args.assemblies or str(graph / "inputs" / "query_paths.normalized.txt")
    )
    if args.partition_list:
        partition_list = absolute(args.partition_list)
    else:
        canonical_list = graph_list_path(graph_partitions)
        canonical_ready = (
            canonical_list.is_file() and canonical_list.stat().st_size > 0
        )
        legacy = next((
            candidate for candidate in (
                graph / "summary" / "Graphs.list",
                graph / "partition_targets.list",
                graph / "summary" / "partition_targets.list",
            )
            if candidate.is_file() and candidate.stat().st_size > 0
        ), None)
        if args.dry_run:
            partition_list = str(
                canonical_list if canonical_ready
                else (legacy or canonical_list)
            )
        else:
            if not canonical_ready and legacy is not None:
                names = read_partition_names(str(legacy))
                write_text_if_changed(
                    canonical_list,
                    "".join(
                        str(graph_partitions / name / f"{name}.FA") + "\n"
                        for name in names
                    ),
                )
            partition_list = str(canonical_list)
    assemblies = read_assembly_list(query_paths)
    assembly_by_name = {row.name: row for row in assemblies}
    if not args.reference:
        if not assemblies:
            raise ValueError('empty graph cohort assembly list')
        args.reference = assemblies[0].name
    if args.reference not in assembly_by_name:
        raise ValueError(
            f"reference {args.reference!r} is absent from {query_paths}"
        )

    requested = read_sample_selection(args.samples)
    if requested:
        duplicate = next((
            name for name, count in Counter(requested).items()
            if count > 1
        ), None)
        if duplicate is not None:
            raise ValueError(f"sample selection contains duplicate {duplicate!r}")
        unknown = [name for name in requested if name not in assembly_by_name]
        if unknown:
            raise ValueError(
                "sample selection is absent from the graph cohort: "
                + ", ".join(unknown)
            )
        selected_set = set(requested)
        selected = [row.name for row in assemblies if row.name in selected_set]
    else:
        selected = [row.name for row in assemblies]
    selected = [args.reference] + [name for name in selected if name != args.reference]

    partitions = read_partition_names(partition_list)
    required_scripts = (
        "match_partition_blocks.py", "GenomeLift.py",
        "graphcigartoref_persample.py", "fill_graphcigartoref_gaps.py",
        "alternative_intervals.py",
        "graphreftovcf_persample.py", "local_reference_templates.py", "lift_local_templates.py",
        "assembly_contigs.py", "fixed_alternatives.py", "minsetref_localize.py", "minsetref_align.py",
        "gfixbreaks.py", "graphvcfmerge.py", "graphvcfmerge_snp.py", "graphvcfmerge_snp_compact.py", "graphvcfmerge_snp_runs.py", "cohort_vcf_merge.py", "graphvcfmerge_kmer.py",
        "graphvcfmerge_nested.py", "graphvcfmerge_scheduler.py", "graphvcfmerge_checkpoints.py",
    )
    missing_scripts = [
        script_folder / name for name in required_scripts
        if not (script_folder / name).is_file()
    ]
    if missing_scripts:
        raise FileNotFoundError(
            "missing calling script(s): " + ", ".join(map(str, missing_scripts))
        )
    from graphvcfmerge_kmer import resolve_kmermatch
    kmermatch = Path(resolve_kmermatch(absolute(args.kmermatch)))

    assemblies_payload = [
        {"name": row.name, "fasta": row.fasta, "fai": row.fai}
        for row in assemblies
    ]
    inputs_root = output / "inputs"
    checkpoints_root = output / "checkpoints"
    layout_config_path = inputs_root / "cohort_call.layout.json"
    liftover_config_path = inputs_root / "cohort_call.liftover.json"
    graphcigar_config_path = inputs_root / "cohort_call.graphcigar.json"
    local_template_config_path = inputs_root / "cohort_call.local_templates.json"
    vcf_config_path = inputs_root / "cohort_call.vcf.json"
    merge_config_path = inputs_root / "cohort_call.merge.json"
    previous_run_config = {}
    previous_run_path = inputs_root / "cohort_call.run.json"
    try:
        previous_run_config = json.loads(
            previous_run_path.read_text(encoding="utf-8")
        )
    except (FileNotFoundError, json.JSONDecodeError):
        pass
    layout_migration_markers = [
        checkpoints_root / "layout.json",
        checkpoints_root / "Graphs.list",
        inputs_root / "query_paths.normalized.txt",
        inputs_root / "selected_query_paths.txt",
        inputs_root / "selected_samples.list",
        *sorted((inputs_root / "match_batches").glob("batch*.list")),
        *sorted((inputs_root / "extraction_batches").glob("batch*.list")),
    ]

    def previous_path_is(key: str, current: str) -> bool:
        previous = previous_run_config.get(key)
        return bool(previous) and absolute(str(previous)) == current

    layout_unchanged = bool(previous_run_config) and all((
        previous_path_is("graph_folder", str(graph)),
        previous_path_is("output_folder", str(output)),
        previous_path_is("query_paths", query_paths),
        previous_path_is("partition_list", partition_list),
        previous_run_config.get("partitions") == partitions,
        previous_run_config.get("reference") == args.reference,
        previous_run_config.get("selected_samples") == selected,
        previous_run_config.get("extract_batch_size") == args.extract_batch_size,
        previous_run_config.get("match_batches") == args.match_batches,
    ))

    # Only semantic inputs belong in stage configs. Runtime capacity and Slurm
    # settings stay in the main run config so changing them can resume existing
    # outputs without refreshing unrelated checkpoints.
    write_stage_config(
        layout_config_path,
        {
            "protocol": "cohort-call-v1",
            "layout_protocol": "cohort-call-layout-v2",
            "graph_folder": str(graph),
            "output_folder": str(output),
            "partitions": partitions,
            "reference": args.reference,
            "selected_samples": selected,
            "assemblies": assemblies_payload,
            "extract_batch_size": args.extract_batch_size,
            "match_batches": args.match_batches,
        },
        migration_markers=(layout_migration_markers if layout_unchanged else None),
    )
    write_stage_config(
        liftover_config_path,
        {
            "protocol": "cohort-call-liftover-v2-same-graph-alignment-priority",
            "max_extension": args.max_extension,
            "non_reference_tags": args.non_reference_tags,
        },
    )
    write_stage_config(
        graphcigar_config_path,
        {
            "protocol": "cohort-call-graphcigar-v1",
            "alternative": alternative_bed,
            "alternative_sha256": alternative_digest,
            "max_extension": args.max_extension,
            "cross_validation_extension": args.cross_validation_extension,
            "realignment": args.realignment,
        },
    )
    write_stage_config(
        local_template_config_path,
        {
            "protocol": "cohort-call-local-templates-v3-fixed-alternatives",
            "local_template_fallback": args.local_template_fallback,
            "reference": args.reference,
            "reference_fasta": assembly_by_name[args.reference].fasta,
            "flank": 10000,
        },
    )
    write_stage_config(
        vcf_config_path,
        {
            "protocol": "cohort-call-vcf-v5-default-insertion-snps",
            "sv_only_size": args.sv_only_size,
            "sv_only_retention_cutoff": (
                (3 * (args.sv_only_size or args.minsvsize)) // 5
                if args.merge_mode == 'svonly' else None
            ),
            "seqcompress": args.seqcompress,
            "pa_tag": args.pa_tag,
            "resolve_conflicts_by_alignment_score": (
                args.resolve_conflicts_by_alignment_score
            ),
            "separate_adjacent_indels": args.separate_adjacent_indels,
            "realignment": args.realignment,
            "edge_blackregion": args.edge_blackregion,
            "max_impute": args.max_extension,
        },
    )
    final_minsvsize = (
        args.sv_only_size
        if args.sv_only_size is not None else args.minsvsize
    )
    cohort_vcf = output / "tmp" / "cohort_merge" / "cohort.sv.vcf"
    merge_settings = {
        "minsvsize": final_minsvsize, "merge_distance": args.merge_distance,
        "size_similarity": args.size_similarity, "sequence_similarity": args.sequence_similarity,
        "var_in_insert": args.var_in_insert, "emit_small": args.merge_mode in ('all', 'svindel'),
    }
    merge_signature = hashlib.sha256(json.dumps([merge_settings, args.merge_mode == "all", "worker-snp-v3"], sort_keys=True).encode()).hexdigest()[:16]
    merge_tmpdir = Path(absolute(
        args.merge_tmpdir or str(output / "tmp" / "graphvcfmerge" / merge_signature)
    ))
    write_stage_config(
        merge_config_path,
        {
            "protocol": "cohort-call-merge-v4-worker-snps",
            "merge_mode": args.merge_mode,
            **merge_settings,
        },
        migration_markers=(cohort_vcf,),
    )

    run_config = {
        "protocol": "cohort-call-v1",
        "graph_folder": str(graph),
        "output_folder": str(output),
        "query_paths": query_paths,
        "partition_list": partition_list,
        "partitions": partitions,
        "reference": args.reference,
        "reference_fasta": assembly_by_name[args.reference].fasta,
        "reference_fai": assembly_by_name[args.reference].fai,
        "selected_samples": selected,
        "assemblies": assemblies_payload,
        "layout_config": str(layout_config_path),
        "liftover_config": str(liftover_config_path),
        "graphcigar_config": str(graphcigar_config_path),
        "local_template_config": str(local_template_config_path),
        "vcf_config": str(vcf_config_path),
        "merge_config": str(merge_config_path),
        "cohort_vcf": str(cohort_vcf),
        "script_folder": str(script_folder),
        "kmermatch": str(kmermatch),
        "cores": args.cores,
        "extract_batch_size": args.extract_batch_size,
        "match_batches": args.match_batches,
        "match_threads": args.match_threads,
        "sample_threads": args.sample_threads,
        "io_jobs": args.io_jobs,
        "format_processes": args.format_processes,
        "merge_processes": args.merge_processes,
        "merge_chr1_resource_multiplier": args.merge_chr1_resource_multiplier,
        "merge_locus_processes": args.merge_locus_processes,
        "merge_tmpdir": str(merge_tmpdir),
        "keep_merge_tmpdir": args.keep_merge_tmpdir,
        # VCF-only switches are supplied to Snakemake separately below. Keep
        # these compatibility keys stable so changing output presentation does
        # not invalidate layout/matching/extraction checkpoints.
        "sv_only_size": None,
        "sv_only_retention_cutoff": None,
        "max_extension": args.max_extension,
        "cross_validation_extension": args.cross_validation_extension,
        "realignment": args.realignment,
        "edge_blackregion": args.edge_blackregion,
        "buffer_size": args.buffer_size,
        "buffer_bytes": args.buffer_bytes,
        "progress_every": args.progress_every,
        "progress_seconds": args.progress_seconds,
        "maxtasksperchild": args.maxtasksperchild,
        "non_reference_tags": args.non_reference_tags,
        "local_template_fallback": args.local_template_fallback,
        "slurm": args.slurm,
        "slurm_jobs": args.slurm_jobs,
        "slurm_account": args.slurm_account,
        "slurm_time": args.slurm_time,
        "slurm_partition": args.slurm_partition,
        "slurm_args": args.slurm_args,
        "match_memory": args.match_memory,
        "extraction_memory": args.extraction_memory,
        "genomelift_memory": args.genomelift_memory,
        "graphcigar_memory": args.graphcigar_memory,
        "gapfill_memory": args.gapfill_memory,
        "vcf_memory": args.vcf_memory,
        "merge_memory": args.merge_memory,
    }
    config_path = output / "inputs" / "cohort_call.run.json"
    write_json_if_changed(config_path, run_config)

    command = [
        args.snakemake,
        "--snakefile", str(HERE / "Snakefile"),
        "--cores", str(max(args.cores, args.slurm_jobs) if args.slurm else args.cores),
        "--rerun-incomplete",
        "--printshellcmds",
    ]
    resources = [f"io_jobs={args.io_jobs}"]
    if args.slurm:
        resources.append(f"slurm_jobs={args.slurm_jobs}")
    command.extend(["--resources", *resources])
    command.extend(["--config", f"run_config={config_path}"])
    if args.dry_run:
        command.append("--dry-run")
    if args.slurm and (merge_tmpdir / "chroms.complete.tsv").is_file() and merge_needs_resume(
        merge_tmpdir, {
            "minsvsize": final_minsvsize, "merge_distance": args.merge_distance,
            "size_similarity": args.size_similarity, "sequence_similarity": args.sequence_similarity,
            "var_in_insert": args.var_in_insert,
            "emit_small": args.merge_mode in ('all', 'svindel'),
        },
    ):
        # The helper will skip checked chromosomes and submit only unfinished
        # ones, even if a previous launch left a whole-stage marker behind.
        command.append("--forcerun=graphvcfmerge_chroms")
    extras = list(args.snakemake_args)
    if extras and extras[0] == "--":
        extras.pop(0)
    command.extend(extras)
    return subprocess.run(command, cwd=str(HERE)).returncode


if __name__ == "__main__":
    try:
        raise SystemExit(main())
    except (OSError, ValueError, RuntimeError, subprocess.CalledProcessError) as error:
        print(f"run_cohort_call_pipeline.py: error: {error}", file=sys.stderr)
        raise SystemExit(1)
