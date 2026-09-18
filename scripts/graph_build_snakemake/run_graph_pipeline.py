#!/usr/bin/env python3
"""Launch the graph workflow with the Snakemake 6.15.1 command line."""
from __future__ import annotations

import argparse
import json
import os
import shlex
import shutil
import subprocess
import sys
from pathlib import Path
from typing import Optional, Sequence


REQUIRED_SNAKEMAKE = "6.15.1"
REQUIRED_PARTITION_DRIVER_OPTIONS = (
    "--reference",
    "--reference-haplotype",
    "--annotations-only",
)
HERE = Path(__file__).resolve().parent
REPOSITORY = HERE.parent
sys.path.insert(0, str(HERE / "workflow" / "scripts"))
from pipeline_inputs import read_assembly_list


def absolute(value: str, base: Optional[Path] = None) -> str:
    path = Path(value).expanduser()
    if not path.is_absolute():
        path = (base or Path.cwd()) / path
    return str(path.resolve())


def annotation_outputs_complete(marker: Path) -> bool:
    """Check the saved output inventory before Snakemake trusts its marker."""
    try:
        saved = json.loads(marker.read_text())
        outputs = saved.get("outputs")
        if not isinstance(outputs, list) or not outputs:
            return False
        for path, size, mtime in outputs:
            artifact = Path(path)
            stat = artifact.stat()
            if (not artifact.is_file() or size <= 0
                    or stat.st_size != size or stat.st_mtime_ns != mtime):
                return False
        return True
    except (OSError, ValueError, TypeError, AttributeError):
        return False


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description=(
            "Build local minimum-set graphs from an assembly list and either "
            "a cohort reference or template-defined grouped BED blocks."
        )
    )
    parser.add_argument(
        "-r", "--reference", default="",
        help=(
            "reference FASTA or NAME in the assembly list; optional with "
            "--bed-grouped, where each supplied block defines its template"
        ),
    )
    parser.add_argument("-q", "--assemblies", required=True,
                        help="NAME FASTA assembly list")
    parser.add_argument("-b", "--bed", help="optional BED4+ block file")
    parser.add_argument(
        "--alternative", metavar="FASTA",
        help=(
            "optional multi-record alternative-locus FASTA; it is appended to "
            "the reference and used by k-mer and alignment cleaning"
        ),
    )
    parser.add_argument(
        "--find-novel-loci", nargs="?", const="__USE_QUERY_LIST__",
        metavar="ASSEMBLY_LIST",
        help=(
            "discover and self-clean novel loci before building blocks; with no "
            "argument, reuse -q/--assemblies"
        ),
    )
    parser.add_argument(
        "--static-block", action="store_true",
        help="keep initial block templates fixed; skip novel discovery and local template refinement",
    )
    parser.add_argument(
        "--novel-locus-jobs", type=int, default=4, metavar="INT",
        help="parallel minsetref discovery jobs within the total -j CPU budget (default: 4)",
    )
    parser.add_argument(
        "--keep-novel-work", action="store_true",
        help="retain minsetref alignment intermediates under GRAPH/novelloci",
    )
    parser.add_argument(
        "--bed-grouped", nargs="?", const=True, default=False, metavar="BED",
        help=(
            "the BED already enumerates repeat-grouped loci; optionally give "
            "the BED path here instead of using -b/--bed"
        ),
    )
    parser.add_argument("-G", "--graph-folder", required=True,
                        help="intermediate/checkpoint graph folder")
    parser.add_argument(
        "-O", "--output-folder", default="", metavar="DIR",
        help=(
            "deprecated compatibility option; graph outputs now remain "
            "under -G/--graph-folder"
        ),
    )
    parser.add_argument(
        "--continue-cohort-call", metavar="DIR",
        help=(
            "after graph construction succeeds, continue with the cohort "
            "variant-calling workflow and write it under DIR"
        ),
    )
    parser.add_argument(
        "--cohort-call-args", default="", metavar="TEXT",
        help=(
            "additional shell-quoted options for the continued cohort caller, "
            "for example '--SVonly 50 --merge-locus-processes 4'"
        ),
    )
    parser.add_argument("-j", "--cores", type=int, default=64)
    parser.add_argument("--script-folder", default=str(REPOSITORY))
    parser.add_argument("--preparation-folder",
                        default=str(REPOSITORY.parent / "tools"))
    parser.add_argument("--large-novels",
                        help="optional established novel_loci.fa")
    parser.add_argument("--partition-batches", type=int, default=100)
    parser.add_argument(
        "--slurm", action="store_true",
        help=(
            "submit novel-locus discovery, per-assembly k-mer searches, and "
            "per-batch minset jobs through sbatch; other rules remain local"
        ),
    )
    parser.add_argument(
        "--slurm-jobs", type=int, default=20, metavar="INT",
        help="maximum simultaneous sbatch wrapper jobs (default: 20)",
    )
    parser.add_argument(
        "--slurm-account", metavar="NAME",
        help=(
            "Slurm account passed to sbatch (required with --slurm unless "
            "--slurm-args is supplied)"
        ),
    )
    parser.add_argument(
        "--slurm-time", default="24:00:00", metavar="TIME",
        help="Slurm time limit passed to sbatch (default: 24:00:00)",
    )
    parser.add_argument(
        "--slurm-partition", metavar="NAME",
        help=(
            "Slurm partition passed to sbatch (required with --slurm unless "
            "--slurm-args is supplied)"
        ),
    )
    parser.add_argument(
        "--slurm-args", default="", metavar="TEXT",
        help=(
            "additional shell-quoted sbatch options, for example "
            "'--time=200:00:00 --account=mchaisso_100 --partition=qcb "+
            "--mem=64G --cpus-per-task=32'; supplied resource-selection "+
            "options override the workflow defaults"
        ),
    )
    parser.add_argument(
        "--novel-slurm-cpus", type=int, default=32, metavar="INT",
        help="CPUs for the optional novel-locus discovery Slurm job (default: 32)",
    )
    parser.add_argument(
        "--novel-slurm-memory", default="128G", metavar="MEMORY",
        help="memory for the optional novel-locus discovery Slurm job (default: 128G)",
    )
    parser.add_argument(
        "--minset-partition-jobs", type=int, default=None, metavar="INT",
        help=(
            "partitions run concurrently inside each minset batch job; "
            "defaults to 1 in fast mode and 4 in rigorous modes. Threads per "
            "partition are calculated automatically from the job CPU allocation"
        ),
    )
    parser.add_argument(
        "--compressgraph", action="store_true",
        help=(
            "deprecated compatibility flag; the graph summary package is now "
            "always written under GRAPH/summary"
        ),
    )
    parser.add_argument(
        "--annotate-novels", action="store_true",
        help=(
            "run the optional cohort local/global/small novel annotation "
            "pass; disabled by default"
        ),
    )
    parser.add_argument("--skip-blastn", action="store_true",
                        help="skip BLAST stages in minsetref_light/minsetgraphmake")
    alignment_modes = parser.add_mutually_exclusive_group()
    alignment_modes.add_argument(
        "--fast", "--fast-mode", dest="alignment_mode",
        action="store_const", const="fast",
        help=(
            "use bulk BLASTN and send only incompletely aligned repeat-heavy "
            "queries to Winnowmap (default); --fast-mode is retained as a "
            "compatibility alias"
        ),
    )
    alignment_modes.add_argument(
        "--rigorous", "--slow-rigorous", dest="alignment_mode",
        action="store_const", const="slow-rigorous",
        help=(
            "run both bulk BLASTN and Winnowmap for every sequence in each "
            "local-graph batch"
        ),
    )
    alignment_modes.add_argument(
        "--max-rigorous", dest="alignment_mode",
        action="store_const", const="max-rigorous",
        help=(
            "disable the k-mer redundancy prefilter, begin redundancy "
            "cleaning with Minimap2, and use --slow-rigorous graph alignment"
        ),
    )
    parser.set_defaults(alignment_mode="fast")
    parser.add_argument(
        "--alignment-timeout", type=int, default=21_600, metavar="SECONDS",
        help=(
            "timeout for each partition alignment command in cohort mode "
            "(default: 21600, six hours); timed-out partitions are recorded "
            "in GRAPH/Graphs/fail.list and skipped"
        ),
    )
    parser.add_argument(
        "--linear", action="store_true",
        help=(
            "generate each _linear.txt output in fast mode; rigorous modes "
            "continue to generate linear output by default"
        ),
    )
    payload_group = parser.add_mutually_exclusive_group()
    payload_group.add_argument(
        "--alignment-payload", dest="alignment_payload", action="store_true",
        default=True,
        help=(
            "include oriented query sequence payloads on I/X operations in "
            "_align.txt (default)"
        ),
    )
    payload_group.add_argument(
        "--no-alignment-payload", dest="alignment_payload", action="store_false",
        help="omit query sequence payloads from I/X operations in _align.txt",
    )
    parser.add_argument("--dry-run", "-n", action="store_true")
    parser.add_argument(
        "--unlock", action="store_true",
        help=(
            "unlock the graph workflow and, with --continue-cohort-call, the "
            "cohort workflow, then exit without running jobs; also accepted "
            "after --. Use only when no other workflow instance is running"
        ),
    )
    parser.add_argument("--snakemake", default="snakemake")
    parser.add_argument("--skip-version-check", action="store_true")
    parser.add_argument(
        "snakemake_args", nargs=argparse.REMAINDER,
        help="extra Snakemake arguments after --",
    )
    return parser


def parse_args(argv: Optional[Sequence[str]] = None) -> argparse.Namespace:
    parser = build_parser()
    args = parser.parse_args(argv)
    # The legacy Snakemake spelling must have the same unlock-only behavior
    # as the launcher option, including the continued cohort workflow.
    args.unlock = args.unlock or "--unlock" in args.snakemake_args
    args.fast_mode = args.alignment_mode == "fast"
    args.slow_rigorous = args.alignment_mode == "slow-rigorous"
    args.max_rigorous = args.alignment_mode == "max-rigorous"
    if isinstance(args.bed_grouped, str):
        grouped_bed = args.bed_grouped
        if args.bed:
            left = Path(args.bed).expanduser().resolve()
            right = Path(grouped_bed).expanduser().resolve()
            if left != right:
                parser.error(
                    "-b/--bed and the BED passed to --bed-grouped disagree"
                )
        args.bed = grouped_bed
        args.bed_grouped = True
    try:
        slurm_args = shlex.split(args.slurm_args)
    except ValueError as error:
        parser.error(f"invalid --slurm-args quoting: {error}")
    args.slurm_args = shlex.join(slurm_args) if slurm_args else ""
    try:
        args.cohort_call_args = shlex.split(args.cohort_call_args)
    except ValueError as error:
        parser.error(f"invalid --cohort-call-args quoting: {error}")
    if args.cohort_call_args and not args.continue_cohort_call:
        parser.error("--cohort-call-args requires --continue-cohort-call DIR")
    if args.cores < 4:
        parser.error("--cores must be at least 4")
    if args.partition_batches < 1:
        parser.error("--partition-batches must be positive")
    if args.slurm_jobs < 1:
        parser.error("--slurm-jobs must be positive")
    if args.novel_locus_jobs < 1:
        parser.error("--novel-locus-jobs must be positive")
    if args.novel_slurm_cpus < 1:
        parser.error("--novel-slurm-cpus must be positive")
    if args.minset_partition_jobs is None:
        args.minset_partition_jobs = 1 if args.fast_mode else 4
    if args.minset_partition_jobs < 1:
        parser.error("--minset-partition-jobs must be positive")
    if args.alignment_timeout < 1:
        parser.error("--alignment-timeout must be positive")
    if not args.novel_slurm_memory:
        parser.error("--novel-slurm-memory cannot be empty")
    if args.slurm and not args.slurm_args and not args.slurm_account:
        parser.error("--slurm-account is required with --slurm")
    if args.slurm and not args.slurm_args and not args.slurm_partition:
        parser.error("--slurm-partition is required with --slurm")
    if args.slurm and not args.slurm_time:
        parser.error("--slurm-time cannot be empty with --slurm")
    if args.slurm_args and not args.slurm:
        parser.error("--slurm-args requires --slurm")
    if args.bed_grouped and not args.bed:
        parser.error("--bed-grouped requires -b/--bed")
    if not args.reference and not args.bed_grouped:
        parser.error(
            "-r/--reference is required unless --bed-grouped supplies "
            "template-defined blocks"
        )
    if not args.reference and (
        args.alternative
        or args.find_novel_loci is not None
        or args.large_novels
        or args.annotate_novels
    ):
        parser.error(
            "template-defined reference-free mode does not accept "
            "--alternative, --find-novel-loci, --large-novels, or "
            "--annotate-novels because those stages require "
            "one cohort-wide reference"
        )
    if args.find_novel_loci is not None and args.large_novels:
        parser.error("--find-novel-loci and --large-novels are mutually exclusive")
    if args.static_block and (args.find_novel_loci is not None or args.annotate_novels):
        parser.error("--static-block cannot be combined with --find-novel-loci or --annotate-novels")
    return args


def cohort_reference_name(reference: str, assemblies) -> str:
    """Resolve the graph reference to the assembly name needed by calling."""
    names = {row.name for row in assemblies}
    if reference in names:
        return reference
    reference_path = absolute(reference)
    matches = [row.name for row in assemblies if row.fasta == reference_path]
    if len(matches) == 1:
        return matches[0]
    raise ValueError(
        "--continue-cohort-call requires -r/--reference to name an assembly "
        "or resolve to exactly one FASTA in -q/--assemblies"
    )


def cohort_call_command(
    args: argparse.Namespace,
    *,
    query_list: str,
    graph_folder: str,
    reference_name: str,
) -> list[str]:
    """Build the post-graph cohort-caller command with shared resources."""
    output = absolute(args.continue_cohort_call)
    if output == graph_folder:
        raise ValueError(
            "--continue-cohort-call DIR must differ from --graph-folder"
        )
    local_limit = args.cores if not args.slurm else 16
    command = [
        sys.executable,
        str(REPOSITORY / "cohort_call_snakemake" /
            "run_cohort_call_pipeline.py"),
        "--reference", reference_name,
        "--graph-folder", graph_folder,
        "--output-folder", output,
        "--assemblies", query_list,
        "--cores", str(args.cores),
        "--sample-threads", str(min(16, local_limit)),
        "--match-threads", str(min(4, local_limit)),
        "--io-jobs", str(min(16, local_limit)),
        "--format-processes", str(min(16, local_limit)),
        "--script-folder", absolute(args.script_folder),
        "--kmermatch", str(Path(absolute(args.script_folder)) / "KmerMatch"),
        "--snakemake", args.snakemake,
        "--skip-version-check",
    ]
    if args.slurm:
        command.extend([
            "--slurm",
            "--slurm-jobs", str(args.slurm_jobs),
            "--slurm-time", args.slurm_time,
        ])
        if args.slurm_account:
            command.extend(["--slurm-account", args.slurm_account])
        if args.slurm_partition:
            command.extend(["--slurm-partition", args.slurm_partition])
        if args.slurm_args:
            command.extend(["--slurm-args", args.slurm_args])
    command.extend(args.cohort_call_args)
    return command


def check_snakemake(command: str, skip: bool) -> None:
    executable = shutil.which(command) if os.path.sep not in command else command
    if not executable or not os.path.isfile(executable):
        raise FileNotFoundError(
            f"Snakemake {REQUIRED_SNAKEMAKE} is required but {command!r} "
            "was not found. Create the environment described in "
            f"{HERE / 'envs' / 'snakemake-6.15.1.yaml'}."
        )
    if skip:
        return
    version = subprocess.check_output([command, "--version"], text=True).strip()
    if version != REQUIRED_SNAKEMAKE:
        raise RuntimeError(
            f"this workflow requires Snakemake {REQUIRED_SNAKEMAKE}; "
            f"{command!r} reports {version!r}. Use --skip-version-check only "
            "after independently confirming compatibility."
        )


def check_partition_driver(script_folder: str) -> None:
    """Reject a partially deployed driver before launching a long workflow."""
    driver = Path(script_folder) / "run_partition_novelty_pipeline.sh"
    if not driver.is_file():
        raise FileNotFoundError(f"partition novelty driver not found: {driver}")
    completed = subprocess.run(
        ["bash", str(driver), "--help"],
        text=True,
        stdout=subprocess.PIPE,
        stderr=subprocess.STDOUT,
    )
    if completed.returncode != 0:
        raise RuntimeError(
            f"cannot inspect partition novelty driver {driver}: "
            f"--help exited with {completed.returncode}"
        )
    advertised = set(completed.stdout.split())
    missing = [
        option for option in REQUIRED_PARTITION_DRIVER_OPTIONS
        if option not in advertised
    ]
    if missing:
        raise RuntimeError(
            f"incompatible partition novelty driver {driver}; missing "
            + ", ".join(missing)
            + ". Synchronize run_partition_novelty_pipeline.sh with this "
            "graph_build_snakemake workflow before starting."
        )


def main(argv: Optional[Sequence[str]] = None) -> int:
    args = parse_args(argv)
    check_snakemake(args.snakemake, args.skip_version_check)
    if args.annotate_novels and not args.unlock:
        check_partition_driver(absolute(args.script_folder))
    query_list = absolute(args.assemblies)
    alternative = absolute(args.alternative) if args.alternative else ""
    if alternative and not os.path.isfile(alternative):
        raise FileNotFoundError(f"alternative FASTA not found: {alternative}")
    if alternative and not os.path.isfile(alternative + ".fai"):
        raise FileNotFoundError(
            f"alternative FASTA index {alternative}.fai not found; run "
            f"`samtools faidx {alternative}` before starting the workflow"
        )
    if args.find_novel_loci == "__USE_QUERY_LIST__":
        novel_query_list = query_list
    elif args.find_novel_loci:
        novel_query_list = absolute(args.find_novel_loci)
        # Validate the independent discovery list, including every adjacent FAI.
        read_assembly_list(novel_query_list)
    else:
        novel_query_list = ""
    input_assemblies = read_assembly_list(query_list)
    graph_folder = absolute(args.graph_folder)
    continued_reference = (
        cohort_reference_name(args.reference, input_assemblies)
        if args.continue_cohort_call and not args.unlock else ""
    )
    reference_defined = bool(args.reference)
    reference = args.reference
    assembly_names = {row.name for row in input_assemblies}
    if reference not in assembly_names and reference not in {"", "."}:
        candidate = Path(reference).expanduser()
        if candidate.exists() or os.path.sep in reference:
            reference = absolute(reference)
    values = {
        "reference": reference,
        "reference_defined": "true" if reference_defined else "false",
        "query_paths": query_list,
        "input_working_directory": str(Path.cwd()),
        "bed": absolute(args.bed) if args.bed else "",
        "bed_grouped": "true" if args.bed_grouped else "false",
        "alternative_fasta": alternative,
        "find_novel_loci": "true" if novel_query_list else "false",
        "static_block": "true" if args.static_block else "false",
        "novel_query_paths": novel_query_list,
        "novel_locus_jobs": str(args.novel_locus_jobs),
        "keep_novel_work": "true" if args.keep_novel_work else "false",
        "graph_folder": graph_folder,
        "script_folder": absolute(args.script_folder),
        "preparation_folder": absolute(args.preparation_folder),
        "large_novel_fasta": absolute(args.large_novels) if args.large_novels else "",
        "threads": str(args.cores),
        "partition_batches": str(args.partition_batches),
        "slurm": "true" if args.slurm else "false",
        "slurm_jobs": str(args.slurm_jobs),
        "slurm_account": args.slurm_account or "",
        "slurm_time": args.slurm_time,
        "slurm_partition": args.slurm_partition or "",
        "slurm_args": args.slurm_args,
        "novel_slurm_cpus": str(args.novel_slurm_cpus),
        "novel_slurm_memory": args.novel_slurm_memory,
        "minset_partition_jobs": str(args.minset_partition_jobs),
        "compressgraph": "true" if args.compressgraph else "false",
        "annotate_novels": "true" if args.annotate_novels else "false",
        "skip_blastn": "true" if args.skip_blastn else "false",
        "alignment_mode": args.alignment_mode,
        "fast_mode": "true" if args.fast_mode else "false",
        "slow_rigorous": "true" if args.slow_rigorous else "false",
        "max_rigorous": "true" if args.max_rigorous else "false",
        "alignment_timeout": str(args.alignment_timeout),
        "linear": "true" if (not args.fast_mode or args.linear) else "false",
        "alignment_payload": "true" if args.alignment_payload else "false",
    }
    command = [
        args.snakemake,
        "--snakefile", str(HERE / "Snakefile"),
        "--cores", str(max(args.cores, args.slurm_jobs) if args.slurm else args.cores),
        "--rerun-incomplete",
        "--printshellcmds",
    ]
    if args.slurm:
        command.extend(["--resources", f"slurm_jobs={args.slurm_jobs}"])
    annotation_marker = Path(graph_folder) / "cohort_annotations" / "_ANNOTATIONS_SUCCESS.json"
    if (args.annotate_novels and not args.unlock and annotation_marker.is_file()
            and not annotation_outputs_complete(annotation_marker)):
        print("Annotation outputs are incomplete; scheduling their reconstruction.", flush=True)
        command.extend(["--forcerun", "refine_and_annotate_novel_sequences"])
    command.append("--config")
    command.extend(
        f"{key}={value}" for key, value in values.items() if value != ""
    )
    if args.dry_run and not args.unlock:
        command.append("--dry-run")
    extras = list(args.snakemake_args)
    if extras and extras[0] == "--":
        extras.pop(0)
    command.extend(option for option in extras if option != "--unlock")
    if args.unlock:
        command.append("--unlock")
        print("Unlocking graph-building workflow.", flush=True)
    graph_result = subprocess.run(command, cwd=str(HERE))
    if graph_result.returncode != 0 or not args.continue_cohort_call:
        return graph_result.returncode
    if args.unlock:
        # Unlock is a maintenance action, not a successful graph build. Do
        # not require graphs.complete or forward normal calling/Slurm options.
        print("Unlocking cohort-calling workflow.", flush=True)
        return subprocess.run(
            [
                sys.executable,
                str(REPOSITORY / "cohort_call_snakemake" /
                    "run_cohort_call_pipeline.py"),
                "--unlock", "--output-folder", absolute(args.continue_cohort_call),
                "--snakemake", args.snakemake, "--skip-version-check",
            ],
            cwd=str(REPOSITORY),
        ).returncode
    cohort_command = cohort_call_command(
        args,
        query_list=str(Path(graph_folder) / "inputs" / "query_paths.normalized.txt"),
        graph_folder=graph_folder,
        reference_name=continued_reference,
    )
    if args.dry_run:
        print(
            "Cohort calling will run after graph construction:\n"
            + shlex.join(cohort_command)
        )
        return 0
    graph_complete = Path(graph_folder) / "Graphs" / "graphs.complete"
    if not graph_complete.is_file():
        raise RuntimeError(
            "graph Snakemake exited successfully but the graph completion "
            f"marker is absent: {graph_complete}; cohort calling was not "
            "started"
        )
    return subprocess.run(cohort_command, cwd=str(REPOSITORY)).returncode


if __name__ == "__main__":
    try:
        raise SystemExit(main())
    except (OSError, ValueError, RuntimeError, subprocess.CalledProcessError) as error:
        print(f"run_graph_pipeline.py: error: {error}", file=sys.stderr)
        raise SystemExit(1)
