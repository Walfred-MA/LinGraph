#!/usr/bin/env python3
"""Build and install every C++ tool, then provision runtime dependencies.

Successful C++ builds are copied from ``scripts/src/<tool>/`` into ``scripts/``.
The normal invocation creates or updates a dedicated Conda environment with
the external tools and pinned Python stack. ``--check-only`` and
``--build-only`` provide non-installing modes. The script exits nonzero when
any requested action fails.
"""

from __future__ import annotations

import argparse
import dataclasses
import json
import os
from pathlib import Path
import re
import shlex
import shutil
import subprocess
import sys
import tempfile
from typing import List, Optional, Sequence, Tuple


ROOT = Path(__file__).resolve().parent
DEFAULT_SRC = ROOT / "src"
MINIMUM_PYTHON = (3, 9)
MINIMUM_HTSLIB = (1, 19)
DEFAULT_CONDA_ENV = "minsetref-snakemake-6.15.1"
CONDA_PACKAGES: Tuple[str, ...] = (
    "python>=3.9",
    "snakemake-minimal=6.15.1",
    "tabulate=0.8.10",
    "numpy",
    "polars>=1.0",
    "mappy",
    "minimap2",
    "samtools",
    "htslib>=1.19",
    "winnowmap",
    "blast>=2.14.0",
    "bcftools",
    "meryl",
)
PIP_ONLY_PACKAGES: Tuple[str, ...] = ("ssw-py", "parasail", "edlib")
ROOT_BINARY_ALIASES = {
    # KmerSearcher builds the lower-case workflow executable.  Install the
    # historical mixed-case name as well so both run_sample_pipeline.py's
    # default and graph_build_snakemake's command name work after one build.
    "KmerSearcher": ("KmerSearcher", "kmer_searcher"),
    "KmerStrd": ("KmerStrd",),
    # Keep the workflow-facing name explicit.  The directory Makefile builds
    # this executable and compile_all() copies it into scripts/.
    "build_local_graphs_parallel_hotspots": (
        "build_local_graphs_parallel_hotspots",
    ),
}


@dataclasses.dataclass(frozen=True)
class ExternalDependency:
    label: str
    command: str
    provided_by: str
    version_args: Tuple[str, ...] = ()
    required_version: Optional[str] = None
    minimum_version: Optional[Tuple[int, ...]] = None


@dataclasses.dataclass(frozen=True)
class PythonDependency:
    module: str
    distribution: str
    purpose: str
    required_version: Optional[str] = None


@dataclasses.dataclass(frozen=True)
class CheckResult:
    label: str
    ok: bool
    detail: str


@dataclasses.dataclass(frozen=True)
class BuildResult:
    directory: Path
    ok: bool
    detail: str
    output: Optional[Path] = None


EXTERNAL_DEPENDENCIES: Tuple[ExternalDependency, ...] = (
    ExternalDependency("minimap2", "minimap2", "minimap2"),
    ExternalDependency("samtools", "samtools", "samtools"),
    ExternalDependency(
        "Winnowmap2", "winnowmap", "winnowmap (Winnowmap2 package)",
    ),
    ExternalDependency(
        "Nucleotide-Nucleotide BLAST", "blastn", "NCBI BLAST+ >=2.14.0",
        version_args=("-version",), minimum_version=(2, 14, 0),
    ),
    ExternalDependency("makeblastdb", "makeblastdb", "NCBI BLAST+"),
    ExternalDependency("bcftools", "bcftools", "bcftools"),
    ExternalDependency("meryl", "meryl", "meryl"),
    ExternalDependency(
        "Snakemake CLI", "snakemake", "snakemake 6.15.1 Python package",
        ("--version",), "6.15.1",
    ),
)


PYTHON_DEPENDENCIES: Tuple[PythonDependency, ...] = (
    PythonDependency(
        "polars", "polars", "parallel in-memory sorting of compact SNP chromosomes",
    ),
    PythonDependency(
        "numpy", "numpy", "chromosome arrays and clustering used by graphvcfmerge.py",
    ),
    PythonDependency(
        "mappy", "mappy", "minimap2 Python bindings used by graphreftovcf.py",
    ),
    PythonDependency(
        "parasail", "parasail", "accelerated alignment in graphcigartoref.py and graphvcfmerge.py",
    ),
    PythonDependency(
        "edlib", "edlib", "exact alignment acceleration in graphvcfmerge.py",
    ),
    PythonDependency(
        "ssw", "ssw-py", "Smith-Waterman payload realignment",
    ),
    PythonDependency(
        "snakemake", "snakemake", "graph-building workflow engine", "6.15.1",
    ),
    PythonDependency(
        "tabulate", "tabulate", "Snakemake 6.15.1 DAG-statistics compatibility",
        "0.8.10",
    ),
)


def command_text(command: Sequence[str]) -> str:
    return shlex.join(str(value) for value in command)


def environment_python(prefix: str) -> Path:
    root = Path(prefix)
    if os.name == "nt":
        candidates = (root / "python.exe", root / "Scripts" / "python.exe")
    else:
        binary_dir = root / "bin"
        candidates = (binary_dir / "python", binary_dir / "python3")
        versioned = sorted(
            path for path in binary_dir.glob("python[0-9]*")
            if re.fullmatch(r"python\d+(?:\.\d+)*", path.name)
        )
        candidates = (*candidates, *versioned)
    for candidate in candidates:
        if candidate.is_file() and os.access(candidate, os.X_OK):
            return candidate
    return candidates[0]


def environment_bin(prefix: str) -> Path:
    return Path(prefix) / ("Scripts" if os.name == "nt" else "bin")


def environment_search_path(prefix: str) -> str:
    current = os.environ.get("PATH", "")
    return os.pathsep.join(value for value in (str(environment_bin(prefix)), current) if value)


def numeric_version(value: str) -> Tuple[int, ...]:
    match = re.search(r"(\d+(?:\.\d+)*)", value.strip())
    if match is None:
        return ()
    return tuple(int(part) for part in match.group(1).split("."))


def check_conda_htslib(prefix: Path) -> CheckResult:
    """Check the installed HTSlib package version from Conda metadata."""
    records = sorted((prefix / "conda-meta").glob("htslib-*.json"))
    versions: List[str] = []
    for record in records:
        try:
            metadata = json.loads(record.read_text())
        except (OSError, ValueError):
            continue
        if metadata.get("name") == "htslib" and metadata.get("version"):
            versions.append(str(metadata["version"]))
    if not versions:
        return CheckResult(
            "HTSlib", False,
            f"no HTSlib package metadata found under {prefix}",
        )
    version_text = max(versions, key=numeric_version)
    version = numeric_version(version_text)
    required = ".".join(map(str, MINIMUM_HTSLIB))
    if version >= MINIMUM_HTSLIB:
        return CheckResult(
            "HTSlib", True,
            f"{prefix} ({version_text}); minimum is {required}",
        )
    return CheckResult(
        "HTSlib", False,
        f"{prefix} has {version_text}; minimum is {required}",
    )


def command_tokens(value: str, option: str) -> List[str]:
    try:
        tokens = shlex.split(value)
    except ValueError as error:
        raise ValueError(f"{option}: cannot parse command {value!r}: {error}") from error
    if not tokens:
        raise ValueError(f"{option}: command cannot be empty")
    resolved = shutil.which(tokens[0])
    if resolved is None:
        raise FileNotFoundError(
            f"{option}: executable not found on PATH: {tokens[0]}"
        )
    tokens[0] = resolved
    return tokens


def source_directories(src_dir: Path) -> List[Path]:
    if not src_dir.is_dir():
        raise FileNotFoundError(f"source directory not found: {src_dir}")
    directories = sorted(
        path for path in src_dir.iterdir()
        if path.is_dir() and not path.name.startswith(".")
    )
    if not directories:
        raise ValueError(f"source directory contains no build folders: {src_dir}")
    return directories


def standalone_source(directory: Path) -> Path:
    suffixes = {".cpp", ".cc", ".cxx"}
    sources = sorted(
        path for path in directory.iterdir()
        if path.is_file() and path.suffix.lower() in suffixes
    )
    preferred = [
        path for path in sources
        if path.stem.lower() == directory.name.lower()
    ]
    if len(preferred) == 1:
        return preferred[0]
    if len(sources) == 1:
        return sources[0]
    if not sources:
        raise ValueError(
            f"{directory}: no Makefile or top-level C++ source was found"
        )
    raise ValueError(
        f"{directory}: multiple standalone C++ sources require a Makefile: "
        + ", ".join(path.name for path in sources)
    )


def built_executable(directory: Path) -> Path:
    candidates = sorted(
        path for path in directory.iterdir()
        if path.is_file() and os.access(path, os.X_OK)
    )
    preferred_names = (directory.name, directory.name.lower())
    for name in preferred_names:
        preferred = directory / name
        if preferred in candidates:
            return preferred
    if len(candidates) == 1:
        return candidates[0]
    if not candidates:
        raise ValueError(f"build produced no top-level executable in {directory}")
    raise ValueError(
        f"build produced multiple executables in {directory}; cannot select one: "
        + ", ".join(path.name for path in candidates)
    )


def install_compiled_binary(result: BuildResult, install_dir: Path) -> BuildResult:
    if not result.ok or result.output is None:
        return result
    source = result.output
    names = (source.name, *ROOT_BINARY_ALIASES.get(result.directory.name, ()))
    installed: List[Path] = []
    try:
        install_dir.mkdir(parents=True, exist_ok=True)
        for name in dict.fromkeys(names):
            destination = install_dir / name
            shutil.copy2(source, destination)
            if not destination.is_file() or not os.access(destination, os.X_OK):
                raise OSError(f"installed file is not executable: {destination}")
            installed.append(destination)
    except OSError as error:
        return BuildResult(
            result.directory, False,
            f"compiled successfully but could not install executable: {error}",
            result.output,
        )
    destinations = ", ".join(str(path) for path in installed)
    return BuildResult(
        result.directory, True,
        f"compiled successfully; installed {destinations}",
        installed[0],
    )


def run_build(command: Sequence[str]) -> Tuple[bool, str]:
    print(f"    $ {command_text(command)}", flush=True)
    try:
        completed = subprocess.run(command, check=False)
    except OSError as error:
        return False, str(error)
    if completed.returncode != 0:
        return False, f"command exited with status {completed.returncode}"
    return True, "compiled successfully"


def standalone_dependency_flags(
    directory: Path,
    dependency_prefix: Optional[Path],
) -> Tuple[List[str], List[str]]:
    """Return compile/link flags for standalone tools with native libraries."""
    if directory.name != "KmerRegionHotspot":
        return [], []
    if dependency_prefix is None:
        raise ValueError(
            "KmerRegionHotspot requires HTSlib development files; activate the "
            "minsetref Conda environment or pass --dependency-prefix PREFIX"
        )
    prefix = dependency_prefix.expanduser().resolve()
    header = prefix / "include" / "htslib" / "faidx.h"
    library_dir = prefix / "lib"
    libraries = tuple(
        library_dir / name for name in (
            "libhts.so", "libhts.dylib", "libhts.a",
        )
    )
    if not header.is_file():
        raise ValueError(
            f"KmerRegionHotspot requires {header}; install "
            f"htslib>={'.'.join(map(str, MINIMUM_HTSLIB))}"
        )
    if not any(path.is_file() for path in libraries):
        raise ValueError(
            f"KmerRegionHotspot requires libhts in {library_dir}; "
            f"install htslib>={'.'.join(map(str, MINIMUM_HTSLIB))}"
        )
    compile_flags = [f"-I{prefix / 'include'}"]
    link_flags = [
        f"-L{library_dir}", f"-Wl,-rpath,{library_dir}", "-lhts",
    ]
    return compile_flags, link_flags


def make_dependency_arguments(
    directory: Path,
    dependency_prefix: Optional[Path],
) -> List[str]:
    """Return Make variable overrides for tools with native dependencies."""
    if directory.name != "KmerSearcher":
        return []
    if dependency_prefix is None:
        raise ValueError(
            "KmerSearcher requires HTSlib development files; activate the "
            "minsetref Conda environment or pass --dependency-prefix PREFIX"
        )
    prefix = dependency_prefix.expanduser().resolve()
    header = prefix / "include" / "htslib" / "faidx.h"
    library_dir = prefix / "lib"
    libraries = tuple(
        library_dir / name for name in (
            "libhts.so", "libhts.dylib", "libhts.a",
        )
    )
    if not header.is_file():
        raise ValueError(
            f"KmerSearcher requires {header}; install "
            f"htslib>={'.'.join(map(str, MINIMUM_HTSLIB))}"
        )
    if not any(path.is_file() for path in libraries):
        raise ValueError(
            f"KmerSearcher requires libhts in {library_dir}; "
            f"install htslib>={'.'.join(map(str, MINIMUM_HTSLIB))}"
        )
    include = shlex.quote(str(prefix / "include"))
    library = shlex.quote(str(library_dir))
    return [
        "CXXFLAGS=-Wall -std=c++17 -O3 -flto "
        f"-DCOUNTINT=uint8_t -I{include}",
        f"LDFLAGS=-L{library} -Wl,-rpath,{library}",
        "LDLIBS=-lhts -lcurl -lz -lpthread",
    ]


def compile_directory(
    directory: Path,
    *,
    jobs: int,
    compiler: Sequence[str],
    make_command: Sequence[str],
    force: bool,
    dependency_prefix: Optional[Path] = None,
) -> BuildResult:
    print(f"[BUILD] {directory.name}", flush=True)
    makefile = directory / "Makefile"
    if makefile.is_file():
        try:
            dependency_arguments = make_dependency_arguments(
                directory, dependency_prefix,
            )
        except (OSError, ValueError) as error:
            return BuildResult(directory, False, str(error))
        command = list(make_command) + ["-C", str(directory)]
        if force:
            command.append("-B")
        command.extend([
            "-j", str(jobs), f"CXX={command_text(compiler)}",
            *dependency_arguments,
        ])
        ok, detail = run_build(command)
        if not ok:
            return BuildResult(directory, False, detail)
        try:
            output = built_executable(directory)
        except (OSError, ValueError) as error:
            return BuildResult(directory, False, str(error))
        return BuildResult(directory, True, detail, output)

    try:
        source = standalone_source(directory)
    except (OSError, ValueError) as error:
        return BuildResult(directory, False, str(error))
    output = directory / source.stem
    try:
        compile_flags, link_flags = standalone_dependency_flags(
            directory, dependency_prefix,
        )
    except (OSError, ValueError) as error:
        return BuildResult(directory, False, str(error))
    command = list(compiler) + [
        "-O3", "-std=c++17", "-Wall", "-Wextra", "-Wpedantic",
        "-pthread", *compile_flags, "-o", str(output), str(source),
        *link_flags,
    ]
    ok, detail = run_build(command)
    if ok and (not output.is_file() or not os.access(output, os.X_OK)):
        ok = False
        detail = f"compiler returned success but output is not executable: {output}"
    return BuildResult(directory, ok, detail, output)


def compile_all(
    src_dir: Path,
    *,
    jobs: int,
    cxx: str,
    make: str,
    force: bool,
    install_dir: Optional[Path] = None,
    dependency_prefix: Optional[Path] = None,
) -> List[BuildResult]:
    directories = source_directories(src_dir)
    try:
        compiler = command_tokens(cxx, "--cxx")
    except (FileNotFoundError, ValueError) as error:
        return [BuildResult(directory, False, str(error)) for directory in directories]
    try:
        make_command = command_tokens(make, "--make")
        make_error = None
    except (FileNotFoundError, ValueError) as error:
        make_command = []
        make_error = str(error)

    results: List[BuildResult] = []
    for directory in directories:
        if (directory / "Makefile").is_file() and make_error is not None:
            print(f"[BUILD] {directory.name}", flush=True)
            results.append(BuildResult(directory, False, make_error))
            continue
        result = compile_directory(
            directory,
            jobs=jobs,
            compiler=compiler,
            make_command=make_command,
            force=force,
            dependency_prefix=dependency_prefix,
        )
        if install_dir is not None:
            result = install_compiled_binary(result, install_dir)
        results.append(result)
    return results


def check_external_dependencies(search_path: Optional[str] = None) -> List[CheckResult]:
    results: List[CheckResult] = []
    for dependency in EXTERNAL_DEPENDENCIES:
        path = shutil.which(dependency.command, path=search_path)
        if path is None:
            results.append(CheckResult(
                dependency.label,
                False,
                f"missing command {dependency.command!r} "
                f"(provided by {dependency.provided_by})",
            ))
        else:
            if (
                dependency.required_version is None
                and dependency.minimum_version is None
            ):
                results.append(CheckResult(
                    dependency.label, True, f"{dependency.command} -> {path}",
                ))
                continue
            try:
                completed = subprocess.run(
                    [path, *dependency.version_args],
                    check=False,
                    stdout=subprocess.PIPE,
                    stderr=subprocess.PIPE,
                    text=True,
                )
            except OSError as error:
                results.append(CheckResult(
                    dependency.label, False,
                    f"cannot run {path} for its version: {error}",
                ))
                continue
            version_lines = completed.stdout.strip().splitlines()
            version_text = version_lines[0].strip() if version_lines else ""
            if dependency.required_version is not None:
                version_ok = version_text == dependency.required_version
            else:
                version_ok = (
                    numeric_version(version_text)
                    >= (dependency.minimum_version or ())
                )
            ok = completed.returncode == 0 and version_ok
            if ok:
                detail = f"{path} (version {version_text})"
            elif completed.returncode != 0:
                detail = (
                    f"{path} version check exited with status "
                    f"{completed.returncode}"
                )
            elif dependency.required_version is not None:
                detail = (
                    f"{path} has version {version_text or 'unknown'}, required "
                    f"exactly {dependency.required_version}"
                )
            else:
                required = ".".join(map(str, dependency.minimum_version or ()))
                detail = (
                    f"{path} has version {version_text or 'unknown'}, required "
                    f"at least {required}"
                )
            results.append(CheckResult(dependency.label, ok, detail))
    return results


def check_blast_database_health(search_path: Optional[str] = None) -> CheckResult:
    """Build and query a tiny v4 database with the selected BLAST tools."""
    label = "BLAST v4 database smoke test"
    makeblastdb = shutil.which("makeblastdb", path=search_path)
    blastn = shutil.which("blastn", path=search_path)
    if makeblastdb is None or blastn is None:
        missing = [
            name for name, path in (
                ("makeblastdb", makeblastdb), ("blastn", blastn),
            )
            if path is None
        ]
        return CheckResult(
            label, False,
            "cannot run smoke test; missing " + ", ".join(missing),
        )

    wrapper = ROOT / "runmakeblastdb"
    if not wrapper.is_file():
        return CheckResult(label, False, f"missing database wrapper: {wrapper}")

    environment = os.environ.copy()
    if search_path is not None:
        environment["PATH"] = search_path
    # Deterministic, non-low-complexity DNA so BLAST's default DUST masking
    # cannot turn this into a false-negative smoke test.
    state = 17
    subject_bases: List[str] = []
    for _ in range(320):
        state = (1103515245 * state + 12345) & 0x7FFFFFFF
        subject_bases.append("ACGT"[(state >> 16) & 3])
    subject_sequence = "".join(subject_bases)
    query_sequence = subject_sequence[40:200]
    try:
        with tempfile.TemporaryDirectory(prefix="minsetref-blast-check-") as temporary:
            directory = Path(temporary)
            subject = directory / "subject.fa"
            query = directory / "query.fa"
            prefix = directory / "subject_db"
            subject.write_text(f">subject\n{subject_sequence}\n")
            query.write_text(f">query\n{query_sequence}\n")
            build = subprocess.run(
                [
                    "bash", str(wrapper),
                    "-in", str(subject), "-out", str(prefix),
                    "-blastdb_version", "4",
                ],
                check=False,
                stdout=subprocess.PIPE,
                stderr=subprocess.PIPE,
                text=True,
                env=environment,
            )
            if build.returncode != 0:
                detail_lines = build.stderr.strip().splitlines()
                detail = detail_lines[-1] if detail_lines else "no diagnostic"
                return CheckResult(
                    label, False,
                    f"v4 database build failed with status "
                    f"{build.returncode}: {detail}",
                )

            required = [
                Path(str(prefix) + suffix)
                for suffix in (".nhr", ".nin", ".nsq", ".seq")
            ]
            missing_or_empty = [
                path.name for path in required
                if not path.is_file() or path.stat().st_size == 0
            ]
            if missing_or_empty:
                return CheckResult(
                    label, False,
                    "v4 build left missing/empty component(s): "
                    + ", ".join(missing_or_empty),
                )
            if Path(str(prefix) + ".ndb").exists():
                return CheckResult(
                    label, False,
                    "database wrapper produced an LMDB/v5 .ndb component",
                )

            search = subprocess.run(
                [
                    blastn, "-query", str(query), "-db", str(prefix),
                    "-task", "blastn", "-word_size", "7", "-outfmt", "6",
                ],
                check=False,
                stdout=subprocess.PIPE,
                stderr=subprocess.PIPE,
                text=True,
                env=environment,
            )
            if search.returncode != 0 or not search.stdout.strip():
                detail_lines = search.stderr.strip().splitlines()
                detail = detail_lines[-1] if detail_lines else "no alignment returned"
                return CheckResult(
                    label, False,
                    f"blastn could not query the generated v4 database: {detail}",
                )
    except OSError as error:
        return CheckResult(label, False, f"smoke test could not run: {error}")

    return CheckResult(
        label, True,
        f"{makeblastdb} created a complete v4 database and {blastn} queried it",
    )


def check_python_dependencies(
    python_executable: Optional[str] = None,
) -> List[CheckResult]:
    executable = python_executable or sys.executable
    version_script = (
        "import sys; print('.'.join(str(value) for value in sys.version_info[:3]))"
    )
    try:
        version_completed = subprocess.run(
            [executable, "-c", version_script],
            check=False,
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
            text=True,
        )
    except OSError as error:
        unavailable = CheckResult(
            "Python", False, f"cannot run Python interpreter {executable}: {error}",
        )
        return [unavailable, *[
            CheckResult(
                f"Python module {dependency.module}", False,
                f"not checked because Python interpreter {executable} is unavailable",
            )
            for dependency in PYTHON_DEPENDENCIES
        ]]
    version_text = version_completed.stdout.strip()
    try:
        version = tuple(int(value) for value in version_text.split("."))
    except ValueError:
        version = ()
    python_ok = version_completed.returncode == 0 and version >= MINIMUM_PYTHON
    results = [CheckResult(
        "Python",
        python_ok,
        f"{executable} ({version_text or 'version unknown'}); "
        f"minimum is {'.'.join(map(str, MINIMUM_PYTHON))}",
    )]
    script = """
import importlib
import importlib.metadata
import sys

module = importlib.import_module(sys.argv[1])
try:
    version = importlib.metadata.version(sys.argv[2])
except importlib.metadata.PackageNotFoundError:
    version = getattr(module, "__version__", "version unknown")
print(version)
"""
    for dependency in PYTHON_DEPENDENCIES:
        try:
            completed = subprocess.run(
                [
                    executable, "-c", script,
                    dependency.module, dependency.distribution,
                ],
                check=False,
                stdout=subprocess.PIPE,
                stderr=subprocess.PIPE,
                text=True,
            )
        except OSError as error:
            results.append(CheckResult(
                f"Python module {dependency.module}", False,
                f"cannot run {executable}: {error}",
            ))
            continue
        if completed.returncode == 0:
            output = completed.stdout.strip().splitlines()
            version_text = output[-1] if output else "version unknown"
            version_ok = (
                dependency.required_version is None
                or version_text == dependency.required_version
            )
            if version_ok:
                detail = (
                    f"imported successfully ({version_text}); "
                    f"{dependency.purpose}"
                )
            else:
                detail = (
                    f"found version {version_text}, required exactly "
                    f"{dependency.required_version}; {dependency.purpose}"
                )
            results.append(CheckResult(
                f"Python module {dependency.module}", version_ok, detail,
            ))
        else:
            error_lines = completed.stderr.strip().splitlines()
            error_text = error_lines[-1] if error_lines else "import failed"
            results.append(CheckResult(
                f"Python module {dependency.module}",
                False,
                f"{error_text}; install distribution {dependency.distribution!r}",
            ))
    return results


def conda_context() -> Optional[List[str]]:
    conda_value = os.environ.get("CONDA_EXE", "conda")
    try:
        conda_command = command_tokens(conda_value, "Conda")
    except (FileNotFoundError, ValueError):
        print(
            "[FAILED] Conda is required but was not found. Install Miniconda "
            "or Miniforge, initialize it for this shell, and rerun install.py.",
            file=sys.stderr,
        )
        return None
    return conda_command


def conda_environment_prefix(
    conda_command: Sequence[str], environment_name: str,
) -> Optional[str]:
    try:
        completed = subprocess.run(
            [*conda_command, "info", "--base"],
            check=False,
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
            text=True,
        )
    except OSError as error:
        print(
            f"[FAILED] Could not query the Conda installation: {error}",
            file=sys.stderr,
        )
        return None
    base_lines = completed.stdout.strip().splitlines()
    if completed.returncode != 0 or not base_lines:
        error_lines = completed.stderr.strip().splitlines()
        detail = error_lines[-1] if error_lines else "unknown Conda error"
        print(f"[FAILED] Could not locate the Conda base: {detail}", file=sys.stderr)
        return None
    conda_base = Path(base_lines[-1]).expanduser().resolve()
    return str(conda_base / "envs" / environment_name)


def selected_conda_environment_prefix(
    conda_command: Sequence[str], environment_name: Optional[str],
) -> Optional[str]:
    """Use an explicitly named environment, otherwise the active non-base env."""
    if environment_name is not None:
        return conda_environment_prefix(conda_command, environment_name)
    active_prefix = os.environ.get("CONDA_PREFIX")
    active_name = os.environ.get("CONDA_DEFAULT_ENV")
    if active_prefix and active_name != "base":
        return str(Path(active_prefix).expanduser().resolve())
    return conda_environment_prefix(conda_command, DEFAULT_CONDA_ENV)


def select_fast_solver(
    conda_command: Sequence[str], conda_prefix: str,
) -> Tuple[List[str], List[str], str]:
    conda_base = Path(conda_prefix).parent.parent
    mamba_name = "mamba.exe" if os.name == "nt" else "mamba"
    micromamba_name = "micromamba.exe" if os.name == "nt" else "micromamba"
    candidates = (
        os.environ.get("MAMBA_EXE"),
        shutil.which("mamba"),
        str(environment_bin(str(conda_base)) / mamba_name),
        shutil.which("micromamba"),
        str(environment_bin(str(conda_base)) / micromamba_name),
    )
    for candidate in candidates:
        if not candidate:
            continue
        path = Path(candidate).expanduser()
        if path.is_file() and os.access(path, os.X_OK):
            return [str(path.resolve())], [], path.name
    # Never allow an implicit fallback to the classic solver. If the libmamba
    # plugin is unavailable or broken, Conda fails quickly with an actionable
    # error instead of spending hours on the full historical index.
    return list(conda_command), ["--solver", "libmamba"], "Conda libmamba"


def install_conda_dependencies(
    environment_name: Optional[str] = None,
) -> Optional[str]:
    conda_command = conda_context()
    if conda_command is None:
        return None
    conda_prefix = selected_conda_environment_prefix(
        conda_command, environment_name,
    )
    if conda_prefix is None:
        return None
    environment_exists = (Path(conda_prefix) / "conda-meta").is_dir()
    if environment_exists:
        executable_results = check_external_dependencies(
            environment_search_path(conda_prefix)
        )
        executable_results.append(check_blast_database_health(
            environment_search_path(conda_prefix)
        ))
        python_results = check_python_dependencies(
            str(environment_python(conda_prefix))
        )
        htslib_result = check_conda_htslib(Path(conda_prefix))
        if all(result.ok for result in [
            *executable_results, *python_results, htslib_result,
        ]):
            print(
                f"\n[OK] Conda environment is already complete: {conda_prefix}",
                flush=True,
            )
            return conda_prefix
    conda_verb = "install" if environment_exists else "create"
    solver_command, solver_options, solver_label = select_fast_solver(
        conda_command, conda_prefix
    )
    print(f"\n[INFO] Using {solver_label} for dependency resolution.", flush=True)
    commands = [
        (
            "Conda runtime and Python dependencies",
            [
                *solver_command, conda_verb, *solver_options, "--yes",
                "--prefix", conda_prefix, "--override-channels",
                "--strict-channel-priority", "-c", "conda-forge", "-c",
                "bioconda", *CONDA_PACKAGES,
            ],
        ),
        (
            "pip-only dependency inside the Conda environment",
            [
                str(environment_python(conda_prefix)), "-m", "pip", "install",
                *PIP_ONLY_PACKAGES,
            ],
        ),
    ]

    for label, command in commands:
        print(f"\n[INSTALL] {label}", flush=True)
        print(f"    $ {command_text(command)}", flush=True)
        try:
            completed = subprocess.run(command, check=False)
        except OSError as error:
            print(
                f"[FAILED] Python dependency installation: {error}",
                file=sys.stderr,
            )
            return None
        if completed.returncode != 0:
            print(
                "[FAILED] Dependency installation exited with status "
                f"{completed.returncode}. Install or repair mamba/libmamba "
                "in the Conda base environment; the installer will not fall "
                "back to the classic solver.",
                file=sys.stderr,
            )
            return None
    print(
        "\nConda environment ready. Activate it with:\n"
        f"  conda activate {conda_prefix}",
        flush=True,
    )
    return conda_prefix


def print_build_summary(results: Sequence[BuildResult]) -> int:
    failures = 0
    print("\nBuild summary")
    print("-------------")
    for result in results:
        status = "OK" if result.ok else "FAILED"
        print(f"[{status}] {result.directory.name}: {result.detail}")
        failures += not result.ok
    return failures


def print_check_summary(
    executable_results: Sequence[CheckResult],
    python_results: Sequence[CheckResult],
) -> int:
    failures = 0
    print("\nEnvironment summary")
    print("-------------------")
    for result in [*executable_results, *python_results]:
        status = "OK" if result.ok else "MISSING"
        print(f"[{status}] {result.label}: {result.detail}")
        failures += not result.ok
    return failures


def print_install_guidance(
    executable_results: Sequence[CheckResult],
    python_results: Sequence[CheckResult],
) -> None:
    missing_executables = [result for result in executable_results if not result.ok]
    missing_python = [result for result in python_results if not result.ok]
    if not missing_executables and not missing_python:
        return
    print("\nRequired dependencies are still missing.")
    print(
        "After installing Conda, rerun install.py. A new environment can also "
        "be created manually with:\n"
        f"  conda create -n {DEFAULT_CONDA_ENV} -c conda-forge -c bioconda "
        + " ".join(CONDA_PACKAGES)
        + f"\n  conda activate {DEFAULT_CONDA_ENV}"
        + "\n  python -m pip install "
        + " ".join(PIP_ONLY_PACKAGES)
    )


def parse_args(argv: Optional[Sequence[str]] = None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=(
            "Force-compile every immediate scripts/src/ tool directory, copy its "
            "executable into scripts/, and install runtime "
            "dependencies into a dedicated Conda environment."
        )
    )
    parser.add_argument(
        "--src", type=Path, default=DEFAULT_SRC,
        help=f"source directory containing tool folders (default: {DEFAULT_SRC})",
    )
    parser.add_argument(
        "-j", "--jobs", type=int, default=max(1, os.cpu_count() or 1),
        help="parallel jobs passed to each Makefile build",
    )
    parser.add_argument(
        "--cxx", default=os.environ.get("CXX", "c++"),
        help="C++ compiler command (default: CXX or c++)",
    )
    parser.add_argument(
        "--make", default=os.environ.get("MAKE", "make"),
        help="Make command (default: MAKE or make)",
    )
    parser.add_argument(
        "--dependency-prefix", type=Path,
        help=(
            "prefix containing native dependency headers/libraries; defaults "
            "to the provisioned or active Conda environment"
        ),
    )
    parser.add_argument(
        "--incremental", action="store_true",
        help="allow Makefiles to reuse existing objects instead of forcing rebuilds",
    )
    parser.add_argument(
        "--no-install-deps", action="store_true",
        help="build and check only; do not create or modify a Conda environment",
    )
    parser.add_argument(
        "--conda-env", default=None, metavar="NAME",
        help=(
            "Conda environment name; by default use the active non-base "
            f"environment, or {DEFAULT_CONDA_ENV!r} when none is active"
        ),
    )
    actions = parser.add_mutually_exclusive_group()
    actions.add_argument(
        "--check-only", action="store_true",
        help="check dependencies without compiling src/",
    )
    actions.add_argument(
        "--build-only", action="store_true",
        help="compile src/ without checking runtime dependencies",
    )
    args = parser.parse_args(argv)
    if args.jobs < 1:
        parser.error("--jobs must be at least 1")
    if (
        args.conda_env is not None
        and (
            not args.conda_env
            or args.conda_env in {".", ".."}
            or Path(args.conda_env).name != args.conda_env
        )
    ):
        parser.error("--conda-env must be a simple environment name, not a path")
    return args


def main(argv: Optional[Sequence[str]] = None) -> int:
    args = parse_args(argv)
    failures = 0
    conda_prefix: Optional[str] = None
    dependency_prefix = (
        args.dependency_prefix.expanduser().resolve()
        if args.dependency_prefix is not None else None
    )
    if dependency_prefix is None and os.environ.get("CONDA_PREFIX"):
        dependency_prefix = Path(os.environ["CONDA_PREFIX"]).expanduser().resolve()

    should_install = (
        not args.check_only and not args.build_only and not args.no_install_deps
    )
    if should_install:
        conda_prefix = install_conda_dependencies(args.conda_env)
        if conda_prefix is None:
            print(
                "[FAILED] Dependency provisioning did not complete; "
                "source compilation was not started.",
                file=sys.stderr,
            )
            return 1
        elif args.dependency_prefix is None:
            dependency_prefix = Path(conda_prefix)

    if not args.check_only:
        try:
            build_results = compile_all(
                args.src.resolve(),
                jobs=args.jobs,
                cxx=args.cxx,
                make=args.make,
                force=not args.incremental,
                install_dir=ROOT,
                dependency_prefix=dependency_prefix,
            )
        except (OSError, ValueError) as error:
            print(f"[FAILED] build discovery: {error}", file=sys.stderr)
            failures += 1
        else:
            failures += print_build_summary(build_results)

    if not args.build_only:
        if conda_prefix is not None:
            selected_path = environment_search_path(conda_prefix)
            executable_results = check_external_dependencies(
                selected_path
            )
            executable_results.append(check_blast_database_health(selected_path))
            python_results = check_python_dependencies(
                str(environment_python(conda_prefix))
            )
        else:
            executable_results = check_external_dependencies()
            executable_results.append(check_blast_database_health())
            python_results = check_python_dependencies()
        htslib_prefix = (
            Path(conda_prefix)
            if conda_prefix is not None
            else dependency_prefix
        )
        if (
            htslib_prefix is not None
            and (htslib_prefix / "conda-meta").is_dir()
        ):
            executable_results.append(check_conda_htslib(htslib_prefix))
        failures += print_check_summary(executable_results, python_results)
        print_install_guidance(executable_results, python_results)

    if failures:
        print(f"\nInstallation check failed with {failures} problem(s).")
        return 1
    if args.check_only:
        print("\nAll environment checks passed.")
    elif args.build_only:
        print("\nAll source folders compiled and installed into the project root.")
    else:
        print("\nAll binaries and Conda environment checks passed.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
