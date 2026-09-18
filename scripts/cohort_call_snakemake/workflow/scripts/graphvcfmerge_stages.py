#!/usr/bin/env python3
"""Run graphvcfmerge chromosome shards locally or in balanced Slurm jobs.

The scan stage writes ``chroms.txt`` only after its manifest and per-locus
metadata are complete.  This helper ranks the shard folders by reference or
alternative-locus length, reserves dedicated boosted jobs for chromosome 1,
and distributes the other loci round-robin. It waits for every job and writes
a marker only after all jobs succeed. The job limit bounds concurrent jobs.
Ranking is O(N log N).
"""

from __future__ import annotations

import argparse
from concurrent.futures import ThreadPoolExecutor, as_completed
import os
from pathlib import Path
import re
import shlex
import subprocess
import sys
from typing import List, Optional, Sequence, Tuple

sys.path.insert(0, str(Path(__file__).resolve().parents[3]))
import graphvcfmerge_checkpoints as checkpoints


ChromEntry = Tuple[str, str, int]


def is_chromosome_one(name: str) -> bool:
    """Match primary chr1 names, including GRCh38 and T2T-CHM13 RefSeq IDs."""
    return bool(re.fullmatch(
        r"(?:chr1|1|NC_000001(?:\.\d+)?|NC_060925(?:\.\d+)?)",
        name, flags=re.IGNORECASE,
    ))


def take_option(values, long_name, short_name=""):
    """Remove value-bearing options, returning the last supplied value."""
    rest = []
    value = None
    index = 0
    while index < len(values):
        token = values[index]
        if token == long_name or (short_name and token == short_name):
            if index + 1 >= len(values):
                raise ValueError(f"missing value for {token}")
            index += 1
            value = values[index]
        elif token.startswith(long_name + "="):
            value = token[len(long_name) + 1:]
        elif short_name and token.startswith(short_name) and not token.startswith("--"):
            value = token[len(short_name):]
        else:
            rest.append(token)
        index += 1
    return rest, value


def scaled_memory(value: str, multiplier: int) -> str:
    """Scale a Slurm memory quantity in its original unit (bare values are MB)."""
    match = re.fullmatch(r"(\d+)([KMGT]?)", str(value).strip(), flags=re.IGNORECASE)
    if match is None:
        raise ValueError(f"invalid Slurm memory {value!r}; use e.g. 64G or 64000")
    # Slurm's special --mem=0 requests all node memory and remains zero.
    return f"{int(match[1]) * multiplier}{match[2].upper()}"


def job_resources(sbatch_args, cpus, memory, command, multiplier):
    """Scale effective allocation overrides and give boosted jobs matching workers."""
    options, override_cpus = take_option(sbatch_args, "--cpus-per-task", "-c")
    options, override_memory = take_option(options, "--mem")
    options, per_cpu_memory = take_option(options, "--mem-per-cpu")
    if override_memory is not None and per_cpu_memory is not None:
        raise ValueError("choose either --mem or --mem-per-cpu in --sbatch-arg")
    allocated_cpus = int(override_cpus if override_cpus is not None else cpus) * multiplier
    if allocated_cpus < 1:
        raise ValueError("--cpus-per-task must be positive")
    if per_cpu_memory is not None:
        # Keeping memory per CPU constant doubles total memory with the CPUs.
        memory_option = "--mem-per-cpu=" + scaled_memory(per_cpu_memory, 1)
    else:
        memory_option = "--mem=" + scaled_memory(
            override_memory if override_memory is not None else memory, multiplier,
        )
    options.extend([f"--cpus-per-task={allocated_cpus}", memory_option])
    worker_command = list(command)
    if multiplier != 1:
        worker_command, _old_processes = take_option(worker_command, "--processes", "-t")
        worker_command.extend(["--processes", str(allocated_cpus)])
    return options, worker_command, allocated_cpus, memory_option


def read_chrom_entries(path: str) -> List[ChromEntry]:
    """Return ``(name, shard-directory, locus-length)`` scan entries."""
    entries: List[ChromEntry] = []
    seen = set()
    with open(path, "rt", encoding="utf-8") as handle:
        for line_number, raw in enumerate(handle, 1):
            line = raw.rstrip("\r\n")
            if not line:
                continue
            fields = line.split("\t")
            name = fields[0].strip()
            target = fields[2].strip() if len(fields) >= 3 else name
            if not name or not target:
                raise ValueError(
                    f"{path}:{line_number}: incomplete chromosome entry"
                )
            if name in seen:
                raise ValueError(
                    f"{path}:{line_number}: duplicate chromosome {name!r}"
                )
            seen.add(name)
            size_field = fields[3] if len(fields) >= 4 else fields[1]
            try:
                locus_length = int(size_field)
            except (ValueError, IndexError) as error:
                raise ValueError(
                    f"{path}:{line_number}: invalid locus length"
                ) from error
            if locus_length < 0:
                raise ValueError(
                    f"{path}:{line_number}: negative locus length"
                )
            entries.append((name, target, locus_length))
    return entries


def read_folder_list(path: str) -> List[str]:
    """Read one balanced job's stage-1 shard-folder paths."""
    folders: List[str] = []
    with open(path, "rt", encoding="utf-8") as handle:
        for line_number, raw in enumerate(handle, 1):
            folder = raw.rstrip("\r\n")
            if not folder:
                continue
            if "," in folder:
                raise ValueError(
                    f"{path}:{line_number}: shard folder contains a comma"
                )
            if not os.path.isfile(os.path.join(folder, "chrom.txt")):
                raise ValueError(
                    f"{path}:{line_number}: not a stage-1 shard folder: "
                    f"{folder}"
                )
            folders.append(folder)
    if not folders:
        raise ValueError(f"empty shard-folder list: {path}")
    return folders


def command_after_separator(values: Sequence[str]) -> List[str]:
    command = list(values)
    if command and command[0] == "--":
        command.pop(0)
    if not command:
        raise ValueError("a graphvcfmerge chromosome command is required")
    if "--chrom" in command or any(
        value.startswith("--chrom=") for value in command
    ):
        raise ValueError("the managed chromosome command must omit --chrom")
    return command


def safe_job_name(value: str) -> str:
    name = re.sub(r"[^A-Za-z0-9_.-]+", "-", value).strip(".-")
    if not name:
        raise ValueError("--job-name contains no safe characters")
    return name


def has_option(values: Sequence[str], long_name: str, short_name: str = "") -> bool:
    for value in values:
        if value == long_name or value.startswith(long_name + "="):
            return True
        if short_name and value.startswith(short_name) and not value.startswith("--"):
            return True
    return False


def atomic_marker(path: str, count: int) -> None:
    destination = Path(path).resolve()
    destination.parent.mkdir(parents=True, exist_ok=True)
    temporary = destination.with_name(
        destination.name + f".tmp.{os.getpid()}"
    )
    try:
        temporary.write_text(
            f"completed\t{count}\n", encoding="utf-8",
        )
        os.replace(str(temporary), str(destination))
    finally:
        try:
            temporary.unlink()
        except FileNotFoundError:
            pass


def atomic_lines(path: Path, lines: Sequence[str]) -> None:
    """Atomically replace a small text manifest."""
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_name(path.name + f".tmp.{os.getpid()}")
    try:
        temporary.write_text("".join(lines), encoding="utf-8")
        os.replace(str(temporary), str(path))
    finally:
        try:
            temporary.unlink()
        except FileNotFoundError:
            pass


def run_group(args: argparse.Namespace) -> None:
    """Process one balanced list of shard folders inside a Slurm job."""
    folders = read_folder_list(args.folder_list)
    command = command_after_separator(args.merge_command)
    _rest, root = take_option(command, "--shards-dir")
    if root:
        context = completion_context(command)
        folders = [folder for folder in folders if not checkpoints.is_done(
            root, (Path(folder) / "chrom.txt").read_text().strip(), context,
        )]
        if not folders:
            print("[merge:skip] every chromosome in this job is checked complete", flush=True)
            return
    command.extend(["--chrom", ",".join(folders)])
    print(
        f"[merge:chrom-job] {len(folders)} shard folder(s): "
        f"{shlex.join(command)}",
        flush=True,
    )
    completed = subprocess.run(command, check=False)
    if completed.returncode != 0:
        raise RuntimeError(
            f"chromosome group failed with status {completed.returncode}"
        )


def completion_context(command):
    context = dict(checkpoints.DEFAULT_SETTINGS)
    for name, cast in (("minsvsize", int), ("merge_distance", int),
                       ("size_similarity", float), ("sequence_similarity", float),
                       ("var_in_insert", int)):
        _remaining, value = take_option(command, "--" + name.replace("_", "-"))
        if value is not None:
            context[name] = cast(value)
    # The merger also accepts a bare --svcutoff, whose argparse const is 20.
    for index, token in enumerate(command):
        if token.startswith("--svcutoff="):
            context["minsvsize"] = int(token.split("=", 1)[1])
        elif token == "--svcutoff":
            following = command[index + 1] if index + 1 < len(command) else ""
            context["minsvsize"] = int(following) if re.fullmatch(r"[+-]?\d+", following) else 20
    context["emit_small"] = "--output-small" in command
    _remaining, insertion_snps = take_option(command, "--insertion-snps")
    if insertion_snps:
        context["insertion_snps"] = insertion_snps
    return context


def run_chroms(args: argparse.Namespace) -> None:
    """Run every saved chromosome and publish completion only on success."""
    entries = read_chrom_entries(args.chrom_list)
    command = command_after_separator(args.merge_command)
    _remaining, command_root = take_option(command, "--shards-dir")
    root = Path(args.shards_dir or command_root or Path(args.chrom_list).resolve().parent)
    context = completion_context(command)
    all_entries = entries
    entries = []
    for entry in all_entries:
        state = checkpoints.status(root, entry[0], context)
        print(f"[merge:status] {entry[0]}: {state}", flush=True)
        if state != "done":
            entries.append(entry)
        else:
            checkpoints.is_done(root, entry[0], context)  # Bind an explicit manual marker to this run.
    marker = Path(args.marker).resolve()
    try:
        marker.unlink()
    except FileNotFoundError:
        pass

    if not entries:
        atomic_marker(str(marker), len(all_entries))
        return

    def publish_complete():
        missing = [name for name, _target, _length in all_entries
                   if not checkpoints.is_done(root, name, context)]
        if missing:
            raise RuntimeError(f"chromosome jobs returned without checked completion: {missing[:6]}")
        atomic_marker(str(marker), len(all_entries))

    if not args.slurm:
        for index, (name, target, _locus_length) in enumerate(entries, 1):
            chrom_command = [*command, "--chrom", target]
            print(
                f"[merge:chrom] {index}/{len(entries)}: {name}: "
                f"{shlex.join(chrom_command)}",
                flush=True,
            )
            completed = subprocess.run(chrom_command, check=False)
            if completed.returncode != 0:
                raise RuntimeError(
                    f"chromosome task {name!r} failed with status "
                    f"{completed.returncode}"
                )
        publish_complete()
        return

    if args.cpus < 1:
        raise ValueError("--cpus must be positive")
    if args.jobs < 1:
        raise ValueError("--jobs must be positive")
    if args.chr1_resource_multiplier < 1:
        raise ValueError("--chr1-resource-multiplier must be positive")
    base_job_name = safe_job_name(args.job_name)
    log_dir = Path(args.log_dir).resolve()
    log_dir.mkdir(parents=True, exist_ok=True)
    sbatch_args = list(args.sbatch_arg or ())
    forbidden = ("--wrap", "--array", "--wait")
    if any(
        value == "--" or any(
            value == option or value.startswith(option + "=")
            for option in forbidden
        )
        for value in sbatch_args
    ):
        raise ValueError(
            "--sbatch-arg cannot set --array, --wait, or --wrap"
        )

    ranked = sorted(
        (
            (locus_length, name, str(Path(target).resolve()))
            for name, target, locus_length in entries
        ),
        key=lambda item: (-item[0], item[1], item[2]),
    )
    boosted = [item for item in ranked
               if args.chr1_resource_multiplier != 1 and is_chromosome_one(item[1])]
    ordinary = [item for item in ranked if item not in boosted]
    groups = [[item] for item in boosted]
    multipliers = [args.chr1_resource_multiplier] * len(groups)
    if ordinary:
        ordinary_count = min(len(ordinary), max(1, args.jobs - len(groups)))
        ordinary_groups = [[] for _ in range(ordinary_count)]
        for rank, item in enumerate(ordinary):
            ordinary_groups[rank % ordinary_count].append(item)
        groups.extend(ordinary_groups)
        multipliers.extend([1] * ordinary_count)
    job_count = len(groups)

    job_list_dir = Path(
        args.job_list_dir
        or (Path(args.chrom_list).resolve().parent / "chrom_jobs")
    ).resolve()
    job_list_dir.mkdir(parents=True, exist_ok=True)
    ordinal_width = len(str(job_count))
    group_specs = []
    assignment_lines = ["job\trank\tlocus_length\tchrom\tfolder\tcpus\tmemory_request\n"]
    rank_by_folder = {
        folder: rank for rank, (_size, _name, folder) in enumerate(ranked, 1)
    }
    for group_index, (group, multiplier) in enumerate(zip(groups, multipliers), 1):
        ordinal = str(group_index).zfill(ordinal_width)
        folder_list = job_list_dir / f"job_{ordinal}.folders.txt"
        atomic_lines(
            folder_list,
            [folder + "\n" for _size, _name, folder in group],
        )
        group_bases = sum(size for size, _name, _folder in group)
        options, worker_command, cpus, memory_option = job_resources(
            sbatch_args, args.cpus, args.memory, command, multiplier,
        )
        group_specs.append((group_index, folder_list, len(group), group_bases,
                            options, worker_command))
        for size, name, folder in group:
            assignment_lines.append(
                f"{group_index}\t{rank_by_folder[folder]}\t{size}\t"
                f"{name}\t{folder}\t{cpus}\t{memory_option}\n"
            )
    atomic_lines(job_list_dir / "assignments.tsv", assignment_lines)

    def submit_group(
        index: int, folder_list: Path, folder_count: int, group_bases: int,
        options: List[str], worker_command: List[str],
    ) -> Tuple[int, int]:
        ordinal = str(index).zfill(ordinal_width)
        job_name = safe_job_name(f"{base_job_name}-{ordinal}")[:128]
        worker = [
            sys.executable,
            str(Path(__file__).resolve()),
            "run-group",
            "--folder-list", str(folder_list),
            "--",
            *worker_command,
        ]
        sbatch_command = [
            args.sbatch,
            "--parsable",
            "--wait",
            "--export=ALL",
            f"--job-name={job_name}",
        ]
        if args.account and not has_option(sbatch_args, "--account", "-A"):
            sbatch_command.append(f"--account={args.account}")
        if args.time and not has_option(sbatch_args, "--time", "-t"):
            sbatch_command.append(f"--time={args.time}")
        if args.partition and not has_option(sbatch_args, "--partition", "-p"):
            sbatch_command.append(f"--partition={args.partition}")
        sbatch_command.extend([
            f"--output={log_dir / (job_name + '-%j.out')}",
            f"--error={log_dir / (job_name + '-%j.err')}",
            *options,
            "--wrap", shlex.join(worker),
        ])
        print(
            f"[SLURM CHROM] job {index}/{job_count}: "
            f"{folder_count} folder(s), {group_bases} locus bases: "
            f"{shlex.join(sbatch_command)}",
            flush=True,
        )
        completed = subprocess.run(sbatch_command, check=False)
        return index, completed.returncode

    failures = []
    with ThreadPoolExecutor(max_workers=min(args.jobs, job_count)) as executor:
        futures = {
            executor.submit(submit_group, *spec): spec[0]
            for spec in group_specs
        }
        for future in as_completed(futures):
            index = futures[future]
            try:
                _completed_index, returncode = future.result()
            except Exception as error:
                failures.append(f"job {index}: {error}")
            else:
                if returncode != 0:
                    failures.append(f"job {index}: status {returncode}")
    if failures:
        details = "; ".join(failures)
        raise RuntimeError(f"Slurm chromosome job(s) failed: {details}")
    publish_complete()


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser()
    subparsers = parser.add_subparsers(dest="action", required=True)

    group = subparsers.add_parser("run-group")
    group.add_argument("--folder-list", required=True)
    group.add_argument("merge_command", nargs=argparse.REMAINDER)

    run = subparsers.add_parser("run-chroms")
    run.add_argument("--chrom-list", required=True)
    run.add_argument("--shards-dir", help="persistent merge folder; inferred from the merger command by default")
    run.add_argument("--marker", required=True)
    run.add_argument("--slurm", action="store_true")
    run.add_argument("--sbatch", default="sbatch")
    run.add_argument("--cpus", type=int, default=1)
    run.add_argument("--memory", default="2G")
    run.add_argument(
        "--chr1-resource-multiplier", type=int, default=2, metavar="INT",
        help="Slurm chr1 CPU/total-memory multiplier, in a dedicated job [2]; 1 disables",
    )
    run.add_argument("--account", default="")
    run.add_argument("--time", default="")
    run.add_argument("--partition", default="")
    run.add_argument("--sbatch-arg", action="append", default=[])
    run.add_argument(
        "--jobs", type=int, default=25,
        help="maximum concurrent Slurm chromosome jobs [25]",
    )
    run.add_argument("--job-list-dir", default="")
    run.add_argument("--job-name", default="merge-cohort-chrom")
    run.add_argument("--log-dir", required=True)
    run.add_argument("merge_command", nargs=argparse.REMAINDER)
    return parser


def main(argv: Optional[Sequence[str]] = None) -> int:
    args = build_parser().parse_args(argv)
    if args.action == "run-group":
        run_group(args)
    else:
        run_chroms(args)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
