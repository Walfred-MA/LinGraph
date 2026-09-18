#!/usr/bin/env python3
"""LinGraph: a small, resumable interface to graph building and assembly calling."""
from __future__ import annotations

import argparse
from contextlib import nullcontext
import csv
import fcntl
import hashlib
import json
import os
from pathlib import Path
import re
import shlex
import shutil
import subprocess
import sys

import run_sample_pipeline as sample_pipeline
import cohort_vcf_merge as cohort_merge
import lingraph_options as cli_options
from graph_list_cache import graph_list_path, template_list_path

ROOT = Path(__file__).resolve().parent
TOOLS = ROOT.parent / "tools"
SAMPLE = re.compile(r"^[A-Za-z0-9.-]+_h[1-9][0-9]*$")


def say(message):
    print(f"[LinGraph] {message}", flush=True)


def absolute(value, base=None):
    path = Path(value).expanduser()
    return ((base or Path.cwd()) / path).resolve() if not path.is_absolute() else path.resolve()


def parser(show_advanced=False):
    p = argparse.ArgumentParser(
        description="LinGraph — build graphs and call assembly variants.",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""Common run options (place after graph or singular):
  -G/--graph DIR          Graph save/resume directory; existing graph for singular
  -O/--output DIR         Calling output directory
  -I/--input-list FILE    Assembly list: NAME FASTA per row
  -r/--reference NAME     Reference assembly name or FASTA
  -t/--threads INT        CPU budget (default: up to 16)
  --rigorous             Run both aligners; singular always uses this mode
  --all / --svonly / --svindel / --snp
                         Variant selection (default: --svonly, or --all
                         with --mc-graph); --svcutoff defaults to 20 bp
  --merge-only           Merge existing per-sample VCFs
  --slurm                Submit work through SLURM
  --slurm-args TEXT      Quoted sbatch options
  --dry-run              Show commands without running or writing files

Graph-only options:
  --mc-graph             Also export cohort.gfa
  --gfa-only             Export existing merged VCFs
  --fast / --fast-mode    Selective alignment (default for graph mode)
  --static-block         Keep initial templates; skip additional novel loci
  --recall-only          Recall and merge from saved calling alignments
  --slurm-jobs INT       Maximum simultaneous SLURM jobs (default: 20)
  --dag-dry-run          Preview actual Snakemake jobs for a completed graph;
                         refresh saved configurations without running jobs

Full options:
  python LinGraph.py graph --help-all
  python LinGraph.py singular --help-all
  python LinGraph.py prepare --help

Examples:
  python LinGraph.py graph -I cohort.list -G savegraph -O calls --all
  python LinGraph.py singular -I new_samples.list -G graph -O calls
  python LinGraph.py --unlock
""")
    p.add_argument("--unlock", action="store_true",
                   help="clear graph-building and cohort-calling Snakemake locks for this repository, then exit")
    p.add_argument("--dry-run", action="store_true", help="with --unlock, show lock directories without removing them")
    modes = p.add_subparsers(dest="mode", title="commands")
    prepare = modes.add_parser('prepare', add_help=False, help='explicitly prepare assemblies; use prepare --help for all options')
    prepare.add_argument('arguments', nargs=argparse.REMAINDER)
    for mode in ("graph", "singular"):
        q = modes.add_parser(mode, formatter_class=argparse.RawDescriptionHelpFormatter,
            help=("build/resume a graph, call cohort samples, and merge VCFs (default: --fast)"
                  if mode == "graph" else
                  "call one or many samples independently against an existing graph (rigorous alignment only)"),
            description=("Build/resume a graph, call every cohort sample, and merge VCFs." if mode == "graph"
                         else "Call one or many sample/haplotype assemblies independently against an existing graph; write one VCF and coverage report per input."),
            epilog=("Examples:\n  python LinGraph.py graph -I cohort.list -G savegraph -O calls\n"
                    "  python LinGraph.py graph -I cohort.list -r CHM13_h1 -G savegraph -O calls --MC-graph"
                    if mode == "graph" else
                    "Examples:\n  python LinGraph.py singular -i query.fa --sample HG002_h1 -G graph -O calls\n"
                    "  python LinGraph.py singular -I samples.list -G graph -r reference.fa -O calls\n"
                    "\nEach input is called independently. Add --merge to also merge the resulting VCFs.")
                    + "\n\nLists: NAME FASTA, one sample/haplotype per row; index: FASTA.fai.\n"
                    "Prepare assemblies separately with LinGraph.py prepare.\n"
                    "LinGraph checks the first sequence and adjacent index; it never prepares inputs.\n"
                    "Relative FASTA paths are relative to the list. Omit -L to use all partitions.\n"
                    "Run locally by default. Repeat the command to resume. --dry-run only prints the plan.")
        q.add_argument("-G", "--graph", "--graph-folder", required=True,
                       help="directory to save or resume graph construction" if mode == "graph" else
                            "existing compact summary/graph directory")
        q.add_argument("-O", "--output", "--output-folder", required=True, help="output directory for VCFs, logs, and run metadata")
        inputs = q.add_mutually_exclusive_group()
        inputs.add_argument("-I", "--input-list", help=(
            "NAME FASTA list for one or many samples; each input is called independently"
            if mode == "singular" else "NAME FASTA list; graph defaults to its saved cohort list"))
        if mode == "singular":
            q.set_defaults(mc_graph=False, gfa_only=False, insertion_only=None, alternative=[])
            inputs.add_argument("-i", "--input", help="one query FASTA (requires --sample)")
            q.add_argument("--sample", help="query name, e.g. HG002_h1, for -i")
            q.add_argument("--reference-caches", metavar="DIR",
                           help="use supplied reference alignment and blocks without cache validation; default: GRAPH/references/NAME_rig")
            q.add_argument("--merge", action="store_true", help="also merge independent calls from multiple inputs; explicit variant-selection options also enable merging")
        q.add_argument("-r", "--reference", help="reference NAME or FASTA; default: first saved cohort assembly")
        q.add_argument("-L", "--graph-list", help="partition names or paths, one per row; default: all")
        q.add_argument("-t", "-j", "--threads", "--cores", type=int,
                       default=max(4 if mode == "graph" else 1, min(16, os.cpu_count() or 1)),
                       help="local CPU budget / calling worker CPUs with SLURM; graph requires >=4 (default: %(default)s)")
        if mode == "graph":
            q.add_argument("--MC-graph", "--mc-graph", dest="mc_graph", action="store_true", help="convert merged VCF to cohort.gfa; defaults to --all, otherwise --svonly")
        alignment = q.add_mutually_exclusive_group()
        if mode == "graph":
            alignment.add_argument("--fast", "--fast-mode", dest="alignment_mode", action="store_const", const="fast",
                                   help="BLASTN first; Winnowmap for qualifying unfinished queries (graph default)")
        alignment.add_argument("--rigorous", "--slow-rigorous", dest="alignment_mode", action="store_const", const="rigorous",
                               help="run BLASTN and Winnowmap for every query (always used by singular)")
        cohort_merge.add_modes(q)
        stages = q.add_mutually_exclusive_group()
        stages.add_argument("--merge-only", action="store_true", help="merge existing per-sample VCFs" +
                            ("; add --mc-graph to also export GFA" if mode == "graph" else ""))
        if mode == "graph":
            stages.add_argument("--gfa-only", action="store_true", help="export existing merged VCFs to GFA, without calling or merging")
        q.add_argument("--vcf-list", help="input VCF list for --merge-only")
        if mode == 'graph':
            q.add_argument("--static-block", action="store_true",
                           help="keep initial block templates fixed; skip additional novel-locus discovery")
            stages.add_argument("--recall-only", action="store_true", help="recall VCFs from saved gap-filled alignments and merge; require existing upstream inputs")
            q.add_argument("--dag-dry-run", action="store_true",
                           help="show actual Snakemake jobs for an existing graph; refreshes saved run configurations but runs no jobs")
            q.add_argument("--snakemake-args", default="", metavar="TEXT",
                           help="quoted extra Snakemake options for cohort calling, e.g. '--forcerun graphreftovcf_sample'")
        q.add_argument("--svcutoff", type=int, default=20, help="SV size boundary; --svindel writes smaller SVs to .indel.vcf (default: 20)")
        if mode == "graph":
            q.add_argument("--insertion-only", nargs="?", const=50, type=int, metavar="SIZE",
                           help="export only insertions in --MC-graph, minimum size [50]")
        q.add_argument("-slurm", "--slurm", action="store_true",
                       help="submit graph workflow stages as SLURM jobs; singular/partial runs use one allocation")
        q.add_argument("--slurm-args", default="", metavar="TEXT", help="quoted sbatch options, passed to each submitted job")
        if mode == "graph":
            q.add_argument("--slurm-jobs", type=int, default=None,
                           help="maximum simultaneous workflow SLURM jobs (default: 20)")
        q.add_argument("--dry-run", action="store_true", default=argparse.SUPPRESS,
                       help="validate inputs and show commands without running or writing files")
        q.add_argument("--help-all", action="help", help="show help including additional settings")
        a = q.add_argument_group("additional settings")
        a.add_argument("--reference-name", help="sample/haplotype name for an external reference FASTA")
        if mode == 'graph':
            a.add_argument("--alternative", action="append", default=[], metavar="FASTA",
                           help="imported alternatives for graph construction")
            a.add_argument("--snakemake", default="snakemake", help="Snakemake executable (requires 6.15.1)")
        a.add_argument("--slurm-account", help="SLURM account, if required by your cluster")
        a.add_argument("--slurm-partition", help="SLURM partition, otherwise cluster default")
        a.add_argument("--slurm-memory", help="override memory per SLURM job; default: stage-specific, or 64G for one allocation")
        a.add_argument("--slurm-time", default="200:00:00", help="allocation time limit (default: %(default)s)")
        if mode == "graph":
            a.add_argument("-b", "--bed", help="optional graph construction BED")
            a.add_argument("--bed-grouped", action="store_true", help="BED already defines repeat groups (requires -b)")
            a.add_argument("--find-novel-loci", nargs="?", const=True,
                           help="discover novel loci; optional assembly list, otherwise use -I")
        else:
            a.add_argument("--satellite-bed", help="optional annotation for coverage categories")
            a.add_argument("--graph-assemblies", help="optional NAME FASTA cohort metadata; summary reconstruction does not read cohort assemblies")
        if not show_advanced:
            for option in a._group_actions:
                option.help = argparse.SUPPRESS
        if mode == "graph":
            cli_options.add_options(q, 'call', 'calling and merge settings', show_advanced)
            cli_options.add_options(q, 'build', 'graph construction settings', show_advanced)
        else:
            cli_options.add_options(q, 'sample', 'sample calling settings', show_advanced)
            cli_options.add_options(q, 'merge', 'merge settings', show_advanced)
        if mode == "graph":
            cli_options.add_options(q, 'gfa', 'GFA export settings (require --mc-graph)', show_advanced)
    return p


def unlock_workflows(dry_run=False):
    # Both launchers run Snakemake with cwd set to their workflow directory.
    # Snakemake 6.15.1 cleanup_locks() removes exactly .snakemake/locks. Do the
    # same without loading a Snakefile, which may require unavailable run data.
    for workflow in ("graph_build_snakemake", "cohort_call_snakemake"):
        locks = ROOT / workflow / ".snakemake" / "locks"
        if not locks.exists():
            say(f"No locks: {locks}")
            continue
        if dry_run:
            say(f"Would unlock: {locks}")
            continue
        shutil.rmtree(locks)
        say(f"Unlocked: {locks}")
    return 0


def check_fasta(path):
    if not path.is_file() or not path.stat().st_size:
        raise ValueError(f"FASTA missing or empty: {path}")
    if any(c.isspace() for c in str(path)):
        raise ValueError(f"FASTA paths cannot contain whitespace: {path}")
    if path.suffix.lower() in {".gz", ".bgz", ".bgzf"}:
        raise ValueError(f"Use an uncompressed assembly FASTA: {path}; run tools/prepare_assemblies.py separately")


def check_prepared(name, fasta):
    """Inspect only the first FASTA record and first FAI row; never modify inputs."""
    native_reference = name in {"CHM13_h1", "HG38_h1"}

    def fail(reason):
        if native_reference:
            raise ValueError(
                f"Assembly {name} is not ready: {reason}.\n"
                "Repair the FASTA/index while preserving its contig names. "
                "Regenerate the adjacent index with: "
                + shlex.join(["samtools", "faidx", str(fasta)])
            )
        command = [sys.executable, str(TOOLS / "prepare_assemblies.py"),
                   "-i", str(fasta), "--name", name, "-O", "prepared_assemblies",
                   "--remask", "--contignamefix"]
        raise ValueError(f"Assembly {name} is not ready: {reason}.\n"
                         "Please mask/format/index it separately with tools/prepare_assemblies.py, for example:\n"
                         + shlex.join(command) + "\nThen use prepared_assemblies/query_paths.prepared.txt with -I. "
                         "LinGraph has not modified this assembly.")

    index = Path(str(fasta) + ".fai")
    if not index.is_file():
        fail(f"missing adjacent index {index}")
    with fasta.open("rb") as handle:
        raw = handle.readline()
        if not raw.startswith(b">"):
            fail("the first line must be a FASTA header")
        fields = raw[1:].split()
        if not fields:
            fail("the first FASTA header is empty")
        contig = fields[0].decode("ascii", errors="replace")
        prefix = name.rsplit("_h", 1)[0] + "#" + name.rsplit("_h", 1)[1]
        if not native_reference and contig != prefix and not contig.startswith(prefix + "#"):
            fail(f"first sequence {contig!r} must be {prefix!r} or start with {prefix + '#'!r}")
        offset = handle.tell()
        with index.open() as fai:
            record = fai.readline().split()
        try:
            valid_index = (len(record) >= 5 and record[0] == contig and int(record[1]) > 0
                           and int(record[2]) == offset and int(record[3]) > 0
                           and int(record[4]) >= int(record[3]))
        except ValueError:
            valid_index = False
        if not valid_index:
            fail(f"the first entry of {index} is invalid or does not match the FASTA")
        if native_reference:
            return
        masked = False
        for line in handle:
            if line.startswith(b">"):
                break
            if any(base in line for base in (b"a", b"c", b"g", b"t")):
                masked = True
                break
        if not masked:
            fail("no lowercase a/c/g/t masking was found in the first sequence")


def read_samples(path, legacy=False, validate=True):
    rows = []
    seen = set()
    for number, line in enumerate(path.read_text().splitlines(), 1):
        if not line.strip() or line.lstrip().startswith("#"):
            continue
        fields = line.split()
        if len(fields) != 2 and not (legacy and len(fields) == 3):
            raise ValueError(f"{path}:{number}: expected NAME FASTA; the index is always FASTA.fai")
        name, filename = fields[:2]
        if not SAMPLE.fullmatch(name) or name in seen:
            raise ValueError(f"{path}:{number}: duplicate or invalid name {name!r}; use SAMPLE_h1, SAMPLE_h2, …")
        fasta = absolute(filename, path.parent)
        if validate:
            check_fasta(fasta)
        if validate and len(fields) == 3 and absolute(fields[2], path.parent) != Path(str(fasta) + ".fai"):
            raise ValueError(f"{path}:{number}: move the index to {fasta}.fai")
        rows.append((name, fasta))
        seen.add(name)
    if not rows:
        raise ValueError(f"Empty sample list: {path}")
    return rows


def write_text(path, content, dry=False):
    if dry:
        return
    path.parent.mkdir(parents=True, exist_ok=True)
    if path.is_file() and path.read_text() == content:
        return
    temporary = path.with_name(path.name + f".tmp.{os.getpid()}")
    temporary.write_text(content)
    temporary.replace(path)


def stamps(paths):
    result = []
    for path in paths:
        path = Path(path)
        stat = path.stat() if path.exists() else None
        result.append([str(path), stat.st_size if stat else None, stat.st_mtime_ns if stat else None])
    return result


def partition_dependency_paths(graphroot, names):
    """Legacy checkpoint paths, constructed without inspecting partition files."""
    return [str(graphroot / name / (name + suffix)) for name in names
            for suffix in ('.FA', '.fasta', 'cache.json')]


class Runner:
    def __init__(self, args):
        self.args = args
        self.work = absolute(args.output) / "lingraph"

    def run(self, name, command, inputs=(), outputs=(), force=False, before=None,
            graph_inputs=(), legacy_graph_inputs=()):
        command = list(map(str, command))
        marker = self.work / "checkpoints" / f"{name}.json"
        stage_inputs = tuple(inputs)
        inputs = (*stage_inputs, *graph_inputs)
        signature = {"command": command, "inputs": stamps(inputs)}
        signature["scripts"] = stamps(Path(c) for c in command if c.endswith(".py") and Path(c).is_file())
        try:
            previous = json.loads(marker.read_text())
        except (FileNotFoundError, ValueError):
            previous = None
        if not force and outputs and all(Path(x).is_file() for x in outputs):
            expected = dict(signature, outputs=stamps(outputs))
            if previous == expected:
                say(f"Reuse {name}")
                return
            if legacy_graph_inputs and isinstance(previous, dict):
                # Replace only the old graph-dependency tail for comparison;
                # commands, scripts, stage inputs and outputs must still match.
                old_inputs = previous.get('inputs', [])
                count = len(stage_inputs)
                completed_at = marker.stat().st_mtime_ns
                if (isinstance(old_inputs, list)
                        and [row[0] for row in old_inputs[count:]
                             if isinstance(row, list) and len(row) == 3] == list(legacy_graph_inputs)
                        and len(old_inputs) == count + len(legacy_graph_inputs)
                        and dict(previous, inputs=old_inputs[:count] + signature['inputs'][count:]) == expected
                        and all(row[2] is None or row[2] <= completed_at
                                for row in signature['inputs'][count:])):
                    write_text(marker, json.dumps(expected, indent=2) + '\n', self.args.dry_run)
                    say(f"Reuse {name} (upgraded graph checkpoint metadata)")
                    return
        say(f"{name}: {shlex.join(command)}")
        if getattr(self.args, 'dag_dry_run', False) and name in {'build-graph', 'call-cohort'}:
            code = subprocess.call(command)
            if code:
                raise RuntimeError(f"{name} DAG preview failed (exit {code})")
            return
        if self.args.dry_run:
            return
        marker.unlink(missing_ok=True)
        for output in outputs:
            Path(output).parent.mkdir(parents=True, exist_ok=True)
        logs = self.work / "logs"
        logs.mkdir(parents=True, exist_ok=True)
        if before:
            before()
        with (logs / f"{name}.log").open("w") as log:
            log.write(shlex.join(command) + "\n")
            log.flush()
            with subprocess.Popen(command, stdout=subprocess.PIPE, stderr=subprocess.STDOUT, text=True) as process:
                try:
                    for line in process.stdout:
                        print(line, end="", flush=True)
                        log.write(line)
                        log.flush()
                    code = process.wait()
                except BaseException:
                    process.terminate()
                    process.wait()
                    raise
            if code:
                raise RuntimeError(f"{name} failed (exit {code}). See {logs / (name + '.log')}; repeat the command to resume.")
        missing = [str(x) for x in outputs if not Path(x).is_file()]
        if missing:
            raise RuntimeError(f"{name} finished without expected output(s): {', '.join(missing)}")
        # Capture post-run inputs: faidx and graph caches may be created by a stage.
        signature["inputs"] = stamps(inputs)
        write_text(marker, json.dumps(dict(signature, outputs=stamps(outputs)), indent=2) + "\n")

    def script(self, name, script, arguments, **kwargs):
        return self.run(name, [sys.executable, ROOT / script, *arguments], **kwargs)

def reference(args, cohort, samples):
    known = dict(cohort)
    known.update(dict(samples))
    if args.reference in known:
        return args.reference, known[args.reference]
    if args.reference:
        fasta = absolute(args.reference)
        check_fasta(fasta)
        matching = [name for name, path in known.items() if path == fasta]
        with fasta.open() as handle:
            first_fields = handle.readline().lstrip(">").split()
        if not first_fields:
            raise ValueError(f"Empty reference FASTA header: {fasta}; run tools/prepare_assemblies.py separately")
        first = first_fields[0]
        prefix = first.split("#")
        inferred = prefix[0] + "_h" + prefix[1] if len(prefix) >= 2 else sample_pipeline.fasta_stem(fasta)
        name = args.reference_name or (matching[0] if matching else inferred)
        if not SAMPLE.fullmatch(name):
            name += "_h1"
        if not SAMPLE.fullmatch(name):
            raise ValueError("Supply --reference-name SAMPLE_h1 for this reference FASTA")
        if (name in known and known[name] != fasta
                and (args.mode != 'singular' or name in dict(samples))):
            raise ValueError(f"Reference name {name} already identifies another FASTA")
        return name, fasta
    defaults = cohort or (samples if args.mode == "graph" else [])
    if defaults:
        return defaults[0]
    raise ValueError("Supply -r reference.fa (or a reference NAME from the saved graph cohort)")


def partition_names(path):
    from align_partition_hotspots import graph_name_from_list_entry
    names = []
    seen = set()
    for line in path.read_text().splitlines():
        if line.strip() and not line.lstrip().startswith("#"):
            name = graph_name_from_list_entry(line.split()[0], str(path))
            if name in {".", ".."} or not re.fullmatch(r"[A-Za-z0-9_.#-]+", name):
                raise ValueError(f"Unsafe graph partition name: {name}")
            if name in seen:
                raise ValueError(f"Duplicate partition {name} in {path}")
            names.append(name)
            seen.add(name)
    if not names:
        raise ValueError(f"Empty partition list: {path}")
    return names


def graph_complete(graph):
    listing = graph_list_path(graph / "Graphs")
    if (
        not (graph / "Graphs" / "graphs.complete").is_file()
        or not listing.is_file() or not listing.stat().st_size
    ):
        return False
    return all((graph / "Graphs" / name / (name + suffix)).is_file()
               and (graph / "Graphs" / name / (name + suffix)).stat().st_size > 0
               for name in partition_names(listing)
               for suffix in (".FA", "_align.txt", "_breakpoint_consistency.bed", "cache.json"))


def graph_mode(args, runner, graph, output, samples, ref):
    name, fasta = ref
    rows = [(name, fasta)] + [(n, f) for n, f in samples if n != name]
    normalized = runner.work / "cohort.list"
    write_text(normalized, "".join(f"{n}\t{f}\n" for n, f in rows),
               args.dry_run and not args.dag_dry_run)
    common = ["-r", name, "-G", graph, "-O", output, "-q", normalized,
              "-j", args.threads, "--snakemake", args.snakemake,
              "--sample-threads", min(16, args.threads), "--match-threads", min(4, args.threads),
              "--io-jobs", min(16, args.threads), "--format-processes", min(16, args.threads),
              "--svcutoff", args.svcutoff, "--" + cohort_merge.resolve_mode(args)]
    common += workflow_slurm_options(args)
    common += cli_options.forward(args, 'call')
    if args.dag_dry_run:
        common.append('--dry-run')
    from cohort_call_snakemake.run_cohort_call_pipeline import parse_args as parse_call
    from graph_build_snakemake.run_graph_pipeline import parse_args as parse_build
    parse_call(list(map(str, common)))
    if args.graph_list:
        common += ["-L", absolute(args.graph_list)]
    if args.snakemake_args:
        common += ['--', *shlex.split(args.snakemake_args)]
    # With explicit construction inputs, let Snakemake check its DAG even for a completed graph.
    call_tuning = {action.dest for action in cli_options.tuning_actions('call')}
    build_tuning = set(cli_options.values(args, 'build')) - call_tuning
    build = bool(args.input_list or args.bed or args.find_novel_loci or args.alternative or build_tuning
                 or args.alignment_mode is not None or not graph_complete(graph))
    if build:
        say("Build/resume graph before calling all cohort assemblies")
        command = ["-r", name, "-q", normalized, "-G", graph, "-j", args.threads,
                   "--snakemake", args.snakemake,
                   "--rigorous" if args.alignment_mode == "rigorous" else "--fast",
                   *workflow_slurm_options(args)]
        command += cli_options.forward(args, 'build')
        if args.static_block:
            command.append("--static-block")
        if args.bed:
            command += ["-b", absolute(args.bed)]
        if args.bed_grouped:
            command += ["--bed-grouped"]
        if args.find_novel_loci:
            command += ["--find-novel-loci"]
            if args.find_novel_loci is not True:
                command.append(absolute(args.find_novel_loci))
        if args.alternative:
            command += ["--alternative", absolute(args.alternative[0])]
        if args.dag_dry_run:
            command.append('--dry-run')
        parse_build(list(map(str, command)))
        runner.script("build-graph", "graph_build_snakemake/run_graph_pipeline.py", command, force=True)
        if not args.dry_run and not graph_complete(graph):
            raise RuntimeError("Graph build did not produce complete calling inputs; inspect the build log")
    # Keep the compact package manifest synchronized with the canonical
    # <GraphFolder>.list emitted by graph construction.
    targets = graph_list_path(graph / "Graphs")
    summary_targets = graph / "summary" / "Graphs.list"
    say(f"Copy valid partition list: {targets} -> {summary_targets}")
    if not args.dry_run:
        write_text(summary_targets, targets.read_text())
    runner.script("call-cohort", "cohort_call_snakemake/run_cohort_call_pipeline.py", common,
                  outputs=cohort_merge.output_paths(output, cohort_merge.resolve_mode(args)), force=True)
    if getattr(args, 'samples', ''):
        from cohort_call_snakemake.run_cohort_call_pipeline import read_sample_selection
        selected = {name, *read_sample_selection(args.samples)}
        rows = [(n, f) for n, f in rows if n in selected]
    return [output / "samples" / n / f"{n}.vcf" for n, _ in rows]


def reference_lengths_hash(fai):
    """Hash sorted contig-name/length pairs, ignoring FASTA bytes and FAI offsets."""
    lengths = {}
    with Path(fai).open() as handle:
        for number, line in enumerate(handle, 1):
            fields = line.split()
            if not fields:
                continue
            try:
                name, length = fields[0], int(fields[1])
            except (IndexError, ValueError):
                raise ValueError(f"{fai}:{number}: invalid chromosome length") from None
            if length < 0 or name in lengths:
                raise ValueError(f"{fai}:{number}: negative length or duplicate chromosome {name}")
            lengths[name] = length
    if not lengths:
        raise ValueError(f"Empty reference index: {fai}")
    content = ''.join(f'{name}\t{length}\n' for name, length in sorted(lengths.items()))
    return hashlib.sha256(content.encode()).hexdigest()


def singular_reference(args, runner, graph, ref, graphroot, listing, index,
                       graph_dependencies, ordered_names, searcher):
    """Reuse completed reference files, or resume an unfinished preparation."""
    refname, reffa = ref
    template_listing = Path(index).with_suffix('')
    supplied = getattr(args, 'reference_caches', None)
    mode = 'fast' if args.alignment_mode == 'fast' else 'rigorous'
    graph_base = graph.parent if graph.name in {'summary', 'Graphs'} else graph
    refdir = absolute(supplied) if supplied else (
        graph_base / 'references' / (refname + ('_fast' if mode == 'fast' else '_rig'))
    )
    align = refdir / f'{refname}.align.txt'
    blocks = refdir / f'{refname}.align.txt_blocks.bed'
    if supplied:
        say(f'Use supplied reference results without validation: {refdir}')
        return align, blocks

    # The directory is selected before reading any identity metadata. A lock
    # lets simultaneous callers share one completed reference preparation.
    if not args.dry_run:
        refdir.mkdir(parents=True, exist_ok=True)
    lock_context = nullcontext() if args.dry_run else (refdir / '.cache.lock').open('a+')
    with lock_context as lock:
        if lock is not None:
            fcntl.flock(lock.fileno(), fcntl.LOCK_EX)
        have_align = align.is_file() and align.stat().st_size > 0
        have_blocks = blocks.is_file() and blocks.stat().st_size > 0
        lengths_hash = reference_lengths_hash(Path(str(reffa) + '.fai'))
        metadata = refdir / 'reference_cache.json'
        try:
            state = json.loads(metadata.read_text())
        except (FileNotFoundError, ValueError):
            state = {}
        old_signature = state.get('signature', {}) if isinstance(state, dict) else {}
        recorded_lengths = old_signature.get('fai_lengths_sha256') if isinstance(old_signature, dict) else None
        if (have_align or have_blocks) and recorded_lengths != lengths_hash:
            raise ValueError(f'Reference cache {refdir} has existing outputs but its FAI-length signature '
                             f'is missing or does not match {reffa}.fai. Use a separate reference name '
                             'for a different reference version, or restore the matching reference_cache.json.')
        if have_align and have_blocks:
            say(f'Reuse reference results: {refdir} (alignment, blocks and FAI lengths match)')
            return align, blocks

        alignment_options = {key: getattr(args, key) for key in ('extension', 'merge_gap', 'timeout')
                             if hasattr(args, key)}
        alignment_payload = bool(getattr(args, 'alignment_payload', False))

        def signature():
            return dict(protocol='reference-results-v3-graph-manifests', reference=refname, mode=mode,
                        fai_lengths_sha256=lengths_hash,
                        partitions_sha256=hashlib.sha256('\n'.join(ordered_names).encode()).hexdigest(),
                        graph_inputs_sha256=hashlib.sha256(json.dumps(stamps(graph_dependencies)).encode()).hexdigest(),
                        template_inputs_sha256=hashlib.sha256(json.dumps(stamps(
                            [template_listing, index]
                        )).encode()).hexdigest(),
                        alignment_options=alignment_options, alignment_payload=alignment_payload)

        expected = signature()
        valid_state = isinstance(state, dict) and isinstance(state.get('completed'), list)
        old_signature = state.get('signature', {}) if valid_state else {}
        if isinstance(old_signature, dict) and old_signature.get('protocol') == 'reference-results-v2-template-hotspots':
            # v2 and v3 use identical alignment/hotspot semantics. Only the
            # graph dependency bookkeeping changed; do not discard completed
            # scientific results just to upgrade that bookkeeping.
            bookkeeping = {'protocol', 'graph_inputs_sha256', 'template_inputs_sha256'}
            same_identity = ({key: value for key, value in old_signature.items() if key not in bookkeeping}
                             == {key: value for key, value in expected.items() if key not in bookkeeping})
            recorded_at = metadata.stat().st_mtime_ns
            if same_identity and all(row[2] is None or row[2] <= recorded_at
                                     for row in stamps(graph_dependencies)):
                state['signature'] = expected
                write_text(metadata, json.dumps(state, indent=2) + '\n', args.dry_run)
                say('Upgraded reference cache metadata; preserving completed stages')
        if (not isinstance(state, dict) or state.get('signature') != expected
                or not isinstance(state.get('completed'), list)):
            state = dict(signature=expected, completed=[])
        completed = set(state.get('completed', []))
        have_align = 'reference-align' in completed and have_align
        have_blocks = 'reference-blocks' in completed and have_blocks

        reflist = refdir / 'reference.list'
        hotspot = refdir / f'{refname}_hotspot.txt'
        stage_order = ('reference-hotspots', 'reference-align', 'reference-blocks')

        def build(stage, action):
            # Invalidate this stage and its dependants before execution; a
            # partial file from a failed job must never become a cache hit.
            completed.difference_update(stage_order[stage_order.index(stage):])
            state.update(signature=expected, completed=sorted(completed))
            write_text(metadata, json.dumps(state, indent=2) + '\n', args.dry_run)
            action()
            completed.add(stage)
            if stage in {'reference-hotspots', 'reference-blocks'}:
                expected.update(signature())
            state.update(signature=expected, completed=sorted(completed))
            write_text(metadata, json.dumps(state, indent=2) + '\n', args.dry_run)

        if not have_align:
            if 'reference-hotspots' not in completed or not hotspot.is_file():
                write_text(reflist, f'{refname}\t{reffa}\n', args.dry_run)
                build('reference-hotspots', lambda: runner.run(
                    'reference-hotspots', [searcher, '-T', template_listing, '-I', reflist,
                     '-o', str(refdir) + os.sep, '-c', 1000, '-w', 2000, '-N', args.threads],
                    inputs=[template_listing, index, reflist, reffa], outputs=[hotspot], force=True))
            options = [item for key, value in alignment_options.items()
                       for item in ('--' + key.replace('_', '-'), value)]
            build('reference-align', lambda: runner.script(
                'reference-align', 'align_partition_hotspots.py',
                ['-i', hotspot, '-q', reffa, '-L', listing, '-G', graphroot, '-s', refname,
                 '-o', align, '-t', getattr(args, 'align_threads', 0) or args.threads,
                 '--threads-per-job', getattr(args, 'threads_per_job', min(4, args.threads)),
                 '--' + mode, '--header',
                 '--alignment-payload' if alignment_payload else '--no-alignment-payload',
                 *options],
                inputs=[hotspot, reffa, listing, *graph_dependencies], outputs=[align], force=True))
        if not have_align or not have_blocks:
            build('reference-blocks', lambda: runner.script(
                'reference-blocks', 'block_partition_alignments.py',
                ['-i', align, '-G', graphroot, '-L', listing, '-o', blocks, '-j', args.threads],
                inputs=[align, listing, *graph_dependencies], outputs=[blocks], force=True))
        return align, blocks


def singular_mode(args, runner, graph, output, samples, ref, cohort):
    refname, reffa = ref
    package = graph if (graph / "local_graphs.tsv").is_file() else graph / "summary"
    graphroot = graph / "Graphs" if (graph / "Graphs").is_dir() else graph
    canonical_listing = graph_list_path(graphroot)
    canonical_names = None
    # A summary-only package still needs reconstruction. An existing graph is
    # used as supplied: do not probe or repair individual partition files.
    if (not canonical_listing.is_file() and (package / "local_graphs.tsv").is_file()
            and (graphroot == graph or not (package / "Graphs.list").is_file())):
        reconstructed = runner.work / "reconstructed"
        graphroot = reconstructed / "Graphs"
        canonical_listing = graph_list_path(graphroot)
        say(f"Reconstruct compact graph package into {reconstructed}")
        command = ["-i", package, "-r", reffa, "--reference-haplotype", refname,
                   "-o", reconstructed, "-j", args.threads, "--resume"]
        runner.script("reconstruct", "reconstruct_local_graph_folders.py", command, force=True)
        if args.dry_run:
            if (package / "Graphs.list").is_file():
                canonical_names = partition_names(package / "Graphs.list")
            else:
                with (package / "local_graphs.tsv").open() as handle:
                    canonical_names = list(dict.fromkeys(
                        row["partition"] for row in csv.DictReader(handle, delimiter="\t")))
    elif not canonical_listing.is_file() and (package / "Graphs.list").is_file():
        canonical_listing = package / "Graphs.list"

    if canonical_names is None and canonical_listing.is_file():
        canonical_names = partition_names(canonical_listing)
    listing = absolute(args.graph_list) if args.graph_list else canonical_listing
    names = partition_names(listing) if args.graph_list else canonical_names
    if not names:
        raise ValueError(f"No graph partition list found in {graph}; supply an existing graph list with -L")
    for name in names:
        if not re.fullmatch(r"[A-Za-z0-9_.#-]+", name) or name in {".", ".."}:
            raise ValueError(f"Invalid partition name: {name}")
    if canonical_names and set(names) == set(canonical_names):
        listing, names = canonical_listing, canonical_names

    # The template list and its prebuilt index belong to the same ordered
    # partition list. Selecting a subset requires its own companion files.
    template_listing = template_list_path(
        graph_list_path(graphroot) if listing == canonical_listing else listing)
    if template_listing.is_file() and partition_names(template_listing) != names:
        template_listing = template_listing.with_name(
            template_listing.name.replace('.template.list', '.selected.template.list'))
        if template_listing.is_file() and partition_names(template_listing) != names:
            raise ValueError(f"Template list order does not match {listing}: {template_listing}")
    index = Path(str(template_listing) + ".bin")
    # Track package-level metadata for resume, not hundreds of thousands of
    # partition paths on the shared filesystem.
    graph_dependencies = [listing, template_listing, index, graphroot / "graphs.complete"]
    legacy_graph_dependencies = partition_dependency_paths(graphroot, names)
    searcher = sample_pipeline.prefer_local_executable("KmerSearcher", ROOT)
    align, blocks = singular_reference(
        args, runner, graph, ref, graphroot, listing, index, graph_dependencies,
        names, searcher,
    )
    vcfs = []
    template_sources = {name: fasta for name, fasta in cohort if fasta.is_file()}
    template_sources.update(samples)
    template_sources[refname] = reffa
    template_list = runner.work / 'template-assemblies.list'
    write_text(template_list, ''.join(f'{n}\t{f}\n' for n, f in template_sources.items()), args.dry_run)
    for number, (name, fasta) in enumerate(samples, 1):
        say(f"Sample {number}/{len(samples)}: {name}")
        options = sample_pipeline.parse_args([
            "--sample", name, "--fasta-query", str(fasta), "--fasta-reference", str(reffa),
            "--reference-sample", refname, "--reference-align", str(align), "--reference-blocks", str(blocks),
            "-G", str(graphroot), "-L", str(listing), "-O", str(output / "samples"),
            "-t", str(args.threads), "--threads-per-job", str(min(4, args.threads)),
            "--format-processes", str(min(16, args.threads)),
            *cli_options.forward(args, 'sample')])
        options.kmer_searcher = searcher
        options.template_list = str(template_listing)
        options.fast = args.alignment_mode == "fast"
        options.template_assemblies = str(template_list)
        options.svcutoff = (
            (3 * args.svcutoff) // 5
            if cohort_merge.resolve_mode(args) == 'svonly' else None
        )
        # Preparation is deliberately separate. Use the validated input unchanged.
        prepared = fasta
        options.fasta_query = str(prepared)
        queryhotspot, hotspots = sample_pipeline.build_hotspot_stages(options, prepared)
        options.hotspot_query = str(queryhotspot)
        paths, stages = sample_pipeline.build_stages(options, ROOT)
        for stage in [*hotspots, *stages]:
            if stage.output_text is not None:
                write_text(stage.outputs[0], stage.output_text, args.dry_run)
            else:
                runner.run(name + "-" + stage.name, stage.command,
                           inputs=stage.inputs, outputs=stage.outputs,
                           graph_inputs=graph_dependencies, legacy_graph_inputs=legacy_graph_dependencies)
                if stage.readiness_marker and not args.dry_run:
                    write_text(stage.readiness_marker, stage.readiness_protocol or "complete\n")
        prefix = output / "samples" / name / f"{name}.coverage"
        coverage = ["--vcf", paths.vcf, "--graphcigar", paths.graphcigartoreffix,
                    "--reference-fai", str(reffa) + ".fai", "--query-fai", str(prepared) + ".fai",
                    "--reference-fasta", reffa, "--query-fasta", prepared, "--sample", name,
                    "-o", prefix, "--threads", args.threads]
        if args.satellite_bed:
            coverage += ["--satellite-bed", absolute(args.satellite_bed)]
        runner.script(name + "-coverage", "compare_vcf_coverage.py", coverage,
                      inputs=[paths.vcf, paths.graphcigartoreffix, reffa, prepared,
                              *([absolute(args.satellite_bed)] if args.satellite_bed else [])],
                      outputs=[Path(str(prefix) + suffix) for suffix in (".summary.tsv", ".missing_reference.bed", ".missing_query.bed")])
        vcfs.append(paths.vcf)
    if len(vcfs) > 1 and (args.merge or args.merge_mode is not None):
        vcflist = runner.work / "vcfs.list"
        write_text(vcflist, "".join(str(p) + "\n" for p in vcfs), args.dry_run)
        mode = cohort_merge.resolve_mode(args)
        runner.script("merge", "cohort_vcf_merge.py", ["run", "--" + mode, "-I", vcflist,
                      "-O", output, "--svcutoff", args.svcutoff, "-t", args.threads,
                      *cli_options.forward(args, 'merge')],
                      inputs=[vcflist, *vcfs], outputs=cohort_merge.output_paths(output, mode))
    return vcfs


def mc_graph(args, runner, graph, output, cohort, samples, ref):
    """Use the calling backbone and its freshly lifted fallback templates."""
    graph_root = graph.parent if graph.name == "summary" else graph
    fixed = ref[1] if args.reference else graph_root / "inputs" / "reference_alternatives_novels.fa"
    templates = (output / 'checkpoints/local_reference_templates.fa' if args.mode == 'graph'
                 else output / 'samples' / samples[0][0] / 'local_reference_templates.fa')
    for fasta in (fixed, templates):
        if not args.dry_run and not fasta.is_file():
            raise ValueError(f"MC-graph requires the backbone/catalog and lifted calling templates: {fasta}")
    catalog = graph_root / "summary" / "alternatives.fasta"
    if not catalog.is_file() and (graph / "alternatives.fasta").is_file():
        catalog = graph / "alternatives.fasta"
    query_sources = dict(cohort)
    query_sources.update(samples)
    query_sources[ref[0]] = ref[1]
    query_list = runner.work / "mc-query_paths.list"
    write_text(query_list, "".join(f"{name}\t{path}\n" for name, path in query_sources.items()), args.dry_run)
    merge_mode = cohort_merge.resolve_mode(args)
    vcf_inputs = cohort_merge.output_paths(output, merge_mode)
    command = ["-v", *vcf_inputs, "-q", query_list,
               "--graph-folder", graph_root, "--reference-haplotype", ref[0],
               "--local-reference-templates", templates, "-o", output / "cohort.gfa",
               "--gfa-mode", "rgfa", "--processes", args.threads,
               "--svonly" if merge_mode == 'svonly' else "--all", "--svcutoff", args.svcutoff]
    if args.reference:
        command += ['-r', fixed]
    if args.insertion_only is not None:
        command += ['--insertion-only', args.insertion_only]
    command += cli_options.forward(args, 'gfa')
    inputs = [*vcf_inputs, query_list, fixed, templates, *query_sources.values()]
    if args.mode == 'singular':
        for value in args.alternative:
            fasta = absolute(value)
            command += ['-a', fasta]
            inputs += [fasta, Path(str(fasta)+'.fai')]
    if catalog.is_file() or args.mode == "graph":
        command += ["--local-path-fasta", catalog]
        inputs.append(catalog)
    inputs += [Path(str(path) + ".fai") for path in (fixed, templates, *query_sources.values(), catalog)]
    inputs += [ROOT / name for name in ("gfa_interval_pipeline.py", "gfa_interval_metadata.py",
               "gfa_query_anchors.py", "gfa_stable_coords.py", "gfa_source_catalog.py",
               "minsetref_core.py", "minsetref_segments.py", "assembly_contigs.py")]
    if args.slurm:
        wrapper = [sys.executable, ROOT / 'graph_build_snakemake/workflow/scripts/pipeline_inputs.py',
                   'run-slurm', '--cpus', args.threads, '--memory', args.slurm_memory or '64G',
                   '--job-name', 'LinGraph-GFA', '--log-dir', runner.work / 'slurm_logs',
                   '--time', args.slurm_time]
        for option, value in (('--account', args.slurm_account), ('--partition', args.slurm_partition)):
            if value:
                wrapper += [option, value]
        wrapper += ['--sbatch-arg=' + value for value in shlex.split(args.slurm_args)]
        runner.run('MC-graph', [*wrapper, '--', sys.executable, ROOT / 'merged_vcf_to_gfa.py', *command],
                   inputs=inputs, outputs=[output / 'cohort.gfa'])
    else:
        runner.script("MC-graph", "merged_vcf_to_gfa.py", command,
                      inputs=inputs, outputs=[output / "cohort.gfa"])


RECALL_SETTINGS = {'format_processes', 'edge_blackregion', 'max_extension', 'realignment',
                   'resolve_conflicts_by_alignment_score', 'pa_tag', 'separate_adjacent_indels',
                   'seqcompress'}


def recall_settings(args):
    return {key: value for key, value in cli_options.values(args, 'call').items()
            if key in RECALL_SETTINGS}


def workflow_slurm_options(args):
    if not args.slurm:
        return []
    options = ["--slurm", "--slurm-jobs", str(args.slurm_jobs or 20),
               "--slurm-time", args.slurm_time]
    extra = (["--mem=" + args.slurm_memory] if args.slurm_memory else []) + shlex.split(args.slurm_args)
    for flag, value in (("--slurm-account", args.slurm_account),
                        ("--slurm-partition", args.slurm_partition),
                        ("--slurm-args", shlex.join(extra))):
        if value:
            options.append(flag + "=" + value)
    return options


def submit(args, argv):
    output = absolute(args.output)
    work = output / "lingraph"
    # Exact option tokens are removed; values and paths are never shell-expanded.
    child = []
    tokens = iter(argv)
    for token in tokens:
        if token in {"--slurm", "-slurm"}:
            continue
        if token in {"--slurm-args", "--slurm-jobs"}:
            next(tokens)
            continue
        if token.startswith(("--slurm-args=", "--slurm-jobs=")):
            continue
        child.append(token)
    command = [sys.executable, str(ROOT / "LinGraph.py"), *child]
    script = work / "run.slurm.sh"
    content = "#!/bin/bash\nset -euo pipefail\ncd " + shlex.quote(str(Path.cwd())) + "\nexec " + shlex.join(command) + "\n"
    write_text(script, content, args.dry_run)
    submit_command = ["sbatch", "--wait", "--parsable", "--job-name=LinGraph", "--nodes=1", "--ntasks=1",
                      f"--cpus-per-task={args.threads}", f"--mem={args.slurm_memory or '64G'}", f"--time={args.slurm_time}",
                      "--output=" + str(work / "slurm-%j.log")]
    if args.slurm_account:
        submit_command += ["--account=" + args.slurm_account]
    if args.slurm_partition:
        submit_command += ["--partition=" + args.slurm_partition]
    submit_command += shlex.split(args.slurm_args)
    submit_command.append(str(script))
    say("SLURM allocation: " + shlex.join(submit_command))
    if args.dry_run:
        print(content)
        return 0
    say(f"Waiting for the job; live pipeline output will be in {work}/slurm-JOBID.log")
    return subprocess.call(submit_command)


def main(argv=None):
    argv = list(sys.argv[1:] if argv is None else argv)
    if argv and argv[0] == 'prepare':
        return subprocess.call([sys.executable, str(TOOLS / 'prepare_assemblies.py'), *argv[1:]])
    p = parser(show_advanced="--help-all" in argv)
    args = p.parse_args(argv)
    if args.unlock:
        if args.mode is not None:
            p.error("--unlock is a standalone command; omit graph/singular and run options")
        return unlock_workflows(args.dry_run)
    if args.mode is None:
        p.error("choose graph or singular, or use --unlock")
    if args.mode == 'prepare':
        p.error('put preparation options after prepare: LinGraph.py prepare --help')
    if getattr(args, 'dag_dry_run', False):
        if args.merge_only or args.recall_only or args.gfa_only:
            p.error("partial runs do not use Snakemake; use --dry-run")
        args.dry_run = True
        if not graph_complete(absolute(args.graph)):
            p.error("--dag-dry-run requires an existing completed graph; use --dry-run to preview a new build")
    try:
        extra_snakemake = shlex.split(getattr(args, 'snakemake_args', ''))
    except ValueError as error:
        p.error(f"invalid --snakemake-args: {error}")
    if any(option.split('=', 1)[0] in {'--unlock', '--dry-run', '-n'} for option in extra_snakemake):
        p.error("use LinGraph --unlock or --dag-dry-run instead of passing them in --snakemake-args")
    try:
        args.slurm_args = shlex.join(shlex.split(args.slurm_args))
    except ValueError as error:
        p.error(f"invalid --slurm-args: {error}")
    if (args.slurm_args or getattr(args, 'slurm_jobs', None) is not None) and not args.slurm:
        p.error("--slurm-args and --slurm-jobs require --slurm")
    if getattr(args, 'slurm_jobs', None) is not None and args.slurm_jobs < 1:
        p.error("--slurm-jobs must be positive")
    partial = args.merge_only or getattr(args, 'recall_only', False) or args.gfa_only
    if partial and getattr(args, 'reference_caches', None):
        p.error('--reference-caches applies to singular calling, not merge/export-only runs')
    if getattr(args, 'static_block', False) and (
        args.find_novel_loci is not None or getattr(args, 'annotate_novels', False)
    ):
        p.error("--static-block cannot be combined with --find-novel-loci or --annotate-novels")
    if getattr(args, 'slurm_jobs', None) is not None and partial:
        p.error("--slurm-jobs controls the full graph workflow; partial runs use one allocation")
    if partial and (args.alignment_mode is not None or extra_snakemake):
        p.error("alignment modes and --snakemake-args apply to full runs, not partial stages")
    if partial and any(getattr(args, name, None) for name in ('bed', 'bed_grouped', 'find_novel_loci', 'graph_list')):
        p.error("construction inputs and partition selection apply to full runs, not partial stages")
    if partial:
        settings = cli_options.values(args, 'call' if args.mode == 'graph' else 'sample')
        if args.mode == 'graph':
            settings.update(cli_options.values(args, 'build'))
        allowed = set(cli_options.values(args, 'merge')) if not args.gfa_only else set()
        if getattr(args, 'recall_only', False):
            allowed |= RECALL_SETTINGS | {'samples'}
        unused = sorted(set(settings) - allowed)
        if unused:
            p.error("these settings require a different stage: " + ', '.join('--' + name.replace('_', '-') for name in unused))
    if args.vcf_list and not args.merge_only:
        p.error("--vcf-list requires --merge-only")
    if args.gfa_only:
        args.mc_graph = True
    if args.alignment_mode == 'rigorous' and getattr(args, 'skip_blastn', False):
        p.error("--rigorous requires both aligners and cannot be combined with --skip-blastn")
    if (cli_options.values(args, 'gfa') or args.insertion_only is not None) and not args.mc_graph:
        p.error("GFA export settings require --mc-graph")
    cli_options.merge_settings(args)
    if args.threads < 1:
        p.error("--threads must be positive")
    if args.svcutoff < 1 or (args.insertion_only is not None and args.insertion_only < 0):
        p.error("--svcutoff must be positive and --insertion-only nonnegative")
    if args.mode == "graph" and args.threads < 4 and not partial:
        p.error("graph mode requires at least 4 threads (existing graph-builder requirement)")
    graph, output = absolute(args.graph), absolute(args.output)
    if graph == output or graph in output.parents or output in graph.parents:
        p.error("-G and -O must be separate, non-nested directories")
    if getattr(args, 'recall_only', False):
        if args.slurm:
            return submit(args, argv)
        from cohort_vcf_recall import run as recall
        from cohort_call_snakemake.run_cohort_call_pipeline import read_sample_selection
        recall(output, cohort_merge.resolve_mode(args), graph=graph, reference=args.reference,
               assembly_list=args.input_list, processes=args.threads, cutoff=args.svcutoff,
               samples=read_sample_selection(args.samples) if getattr(args, 'samples', '') else None,
               caller_settings=recall_settings(args), merge_settings=cli_options.merge_settings(args),
               dry_run=args.dry_run)
        if not args.mc_graph:
            return 0
    if args.merge_only:
        if args.slurm:
            return submit(args, argv)
        cohort_merge.run(output, cohort_merge.resolve_mode(args),
                         listing=absolute(args.vcf_list) if args.vcf_list else None,
                         cutoff=args.svcutoff, dry_run=args.dry_run, **cli_options.merge_settings(args))
        if not args.mc_graph:
            return 0
    if args.gfa_only and args.slurm:
        return submit(args, argv)
    if args.mode == "singular":
        if not (args.input or args.input_list) and not partial:
            p.error("singular calling requires -i/--input or -I/--input-list")
        if not graph.is_dir():
            p.error("singular mode requires an existing graph directory or compact summary directory")
        if bool(args.input) != bool(args.sample):
            p.error("-i and --sample must be used together")
        if args.sample and not SAMPLE.fullmatch(args.sample):
            p.error("--sample must have the form SAMPLE_h1")
    else:
        if args.bed_grouped and not args.bed:
            p.error("--bed-grouped requires -b BED")
        if len(args.alternative) > 1:
            p.error("graph construction accepts one --alternative FASTA")
    for filename in [args.graph_list, *args.alternative, getattr(args, "bed", None), getattr(args, "satellite_bed", None)]:
        if filename and not absolute(filename).is_file():
            p.error(f"Input file not found: {filename}")
    saved = graph / "inputs" / "query_paths.normalized.txt"
    if graph.name == "summary" and not saved.is_file():
        saved = graph.parent / "inputs" / "query_paths.normalized.txt"
    if getattr(args, "graph_assemblies", None):
        saved = absolute(args.graph_assemblies)
        if not saved.is_file():
            p.error(f"Graph assembly list not found: {saved}")
    # Saved build-cohort paths are provenance in singular mode. A distributed
    # summary must work after those source assemblies have been removed.
    cohort = read_samples(saved, legacy=True, validate=args.mode != 'singular') if saved.is_file() else []
    if args.input_list:
        samples = read_samples(absolute(args.input_list))
    elif args.mode == "singular" and args.input:
        fasta = absolute(args.input)
        check_fasta(fasta)
        samples = [(args.sample, fasta)]
    else:
        samples = cohort
    if not samples:
        p.error("Supply -I cohort.list to build or resume this graph")
    ref = reference(args, cohort, samples)
    checked_samples = samples
    if args.mode == 'graph' and getattr(args, 'samples', ''):
        from cohort_call_snakemake.run_cohort_call_pipeline import read_sample_selection
        selected = read_sample_selection(args.samples)
        if not selected:
            p.error('--samples must select at least one sample')
        if len(set(selected)) != len(selected):
            p.error('--samples contains duplicate names')
        unknown = set(selected) - {ref[0], *(name for name, _fasta in samples)}
        if unknown:
            p.error('--samples contains names absent from the cohort: ' + ', '.join(sorted(unknown)))
        checked_samples = [(name, fasta) for name, fasta in samples if name in selected]
    checked = set()
    # Saved graph sources may use legacy contig names. Only assemblies being
    # called and the chosen reference must satisfy the new-input convention.
    for name, fasta in [*checked_samples, ref]:
        if (name, fasta) not in checked:
            check_prepared(name, fasta)
            checked.add((name, fasta))
    for alternative in args.alternative:
        fasta = absolute(alternative)
        check_fasta(fasta)
        if not Path(str(fasta) + ".fai").is_file():
            raise ValueError(f"Missing adjacent index {fasta}.fai; index the alternative FASTA separately before running LinGraph")
    say(f"Preparation checks passed for {len(checked)} assemblies (first sequence and adjacent index only)")
    if args.slurm and args.mode == "singular":
        return submit(args, argv)
    runner = Runner(args)
    say(f"Mode: {args.mode}; samples: {len(samples)}; reference: {ref[0]}; CPUs: {args.threads}; execution: {'SLURM workflow jobs' if args.slurm else 'local'}")
    say(f"Graph: {graph}\nOutput: {output}")
    write_text(runner.work / "run.json", json.dumps({"arguments": vars(args), "reference": [ref[0], str(ref[1])],
               "samples": [[n, str(f)] for n, f in samples]}, indent=2) + "\n", args.dry_run)
    if partial:
        vcfs = []
    elif args.mode == "graph":
        vcfs = graph_mode(args, runner, graph, output, samples, ref)
    else:
        vcfs = singular_mode(args, runner, graph, output, samples, ref, cohort)
    if args.mc_graph:
        mc_graph(args, runner, graph, output, cohort, samples, ref)
    if getattr(args, 'dag_dry_run', False):
        say('DAG preview complete; no workflow jobs were executed.')
    else:
        say("Plan complete; no commands were executed." if args.dry_run else "Complete.")
    for path in vcfs:
        say(f"Sample VCF: {path}")
    if (args.mode == "graph" or partial or len(vcfs) > 1) and (
            args.mode == "graph" or args.merge or args.mc_graph or args.merge_mode is not None):
        for path in cohort_merge.output_paths(output, cohort_merge.resolve_mode(args)):
            say(f"Merged VCF: {path}")
    if args.mc_graph:
        say(f"GFA: {output / 'cohort.gfa'}")
    if args.mode == "singular":
        say(f"Coverage: {output}/samples/NAME/NAME.coverage.summary.tsv (plus missing-region BEDs)")
    return 0


if __name__ == "__main__":
    try:
        raise SystemExit(main())
    except (ValueError, OSError, RuntimeError) as error:
        print(f"[LinGraph] ERROR: {error}", file=sys.stderr)
        raise SystemExit(1)
    except KeyboardInterrupt:
        print("[LinGraph] Interrupted; repeat the command to resume.", file=sys.stderr)
        raise SystemExit(130)
