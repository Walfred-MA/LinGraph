#!/usr/bin/env python3
"""Maintain graph alignment lists and separate template-only hotspot indexes."""

from __future__ import annotations

import argparse
import fcntl
import os
from pathlib import Path
import re
import subprocess
import tempfile
from typing import Iterable, List, Optional, Sequence, Tuple


def graph_list_path(graph_folder: Path) -> Path:
    """Return the sibling ``<GraphFolder>.list`` path."""
    folder = graph_folder.expanduser().resolve()
    return Path(str(folder) + ".list")


def template_list_path(graph_listing: Path) -> Path:
    """Name the template list paired with an ordered graph alignment list."""
    listing = Path(graph_listing)
    stem = str(listing)[:-5] if str(listing).endswith(".list") else str(listing)
    return Path(stem + ".template.list")


def _natural_key(path: Path) -> Tuple[object, ...]:
    return tuple(
        int(piece) if piece.isdigit() else piece.lower()
        for piece in re.split(r"(\d+)", path.parent.name)
    )


def discover_valid_graph_fastas(graph_folder: Path) -> List[Path]:
    """Return naturally ordered, nonempty canonical graph FASTAs."""
    folder = graph_folder.expanduser().resolve()
    if not folder.is_dir():
        raise NotADirectoryError(f"graph folder not found: {folder}")
    graphs = []
    for directory in folder.iterdir():
        if not directory.is_dir() or directory.name.startswith("."):
            continue
        graph = directory / f"{directory.name}.FA"
        if graph.is_file() and graph.stat().st_size > 0:
            graphs.append(graph.resolve())
    graphs.sort(key=_natural_key)
    if not graphs:
        raise ValueError(f"no valid local graph FASTAs found under {folder}")
    return graphs


def _validated_graph_fastas(
    graph_folder: Path, graph_fastas: Optional[Iterable[Path]],
) -> List[Path]:
    folder = graph_folder.expanduser().resolve()
    if graph_fastas is None:
        return discover_valid_graph_fastas(folder)
    graphs: List[Path] = []
    seen = set()
    for value in graph_fastas:
        graph = Path(value).expanduser().resolve()
        if graph in seen:
            continue
        expected = folder / graph.parent.name / f"{graph.parent.name}.FA"
        if graph != expected or not graph.is_file() or graph.stat().st_size == 0:
            raise ValueError(f"invalid reconstructed graph FASTA: {graph}")
        seen.add(graph)
        graphs.append(graph)
    graphs.sort(key=_natural_key)
    if not graphs:
        raise ValueError(f"no valid local graph FASTAs found under {folder}")
    return graphs


def _read_valid_existing_list(listing: Path, graph_folder: Path) -> Optional[List[Path]]:
    """Read a canonical list without widening its saved valid partition set."""
    try:
        lines = listing.read_text().splitlines()
    except (OSError, UnicodeError):
        return None
    folder = graph_folder.expanduser().resolve()
    graphs: List[Path] = []
    seen = set()
    for raw in lines:
        value = raw.strip()
        if not value or value.startswith("#") or len(value.split()) != 1:
            return None
        graph = Path(value).expanduser().resolve()
        expected = folder / graph.parent.name / f"{graph.parent.name}.FA"
        if (
            graph in seen or graph != expected
            or not graph.is_file() or graph.stat().st_size == 0
        ):
            return None
        seen.add(graph)
        graphs.append(graph)
    return graphs or None


def _local_kmer_searcher(script_dir: Path) -> Path:
    for name in ("KmerSearcher", "kmer_searcher"):
        candidate = script_dir / name
        if candidate.is_file() and os.access(candidate, os.X_OK):
            return candidate.resolve()
    raise FileNotFoundError(
        f"KmerSearcher is missing from script folder {script_dir}; "
        "run install.py from the active Conda environment"
    )


def _atomic_write(path: Path, content: str) -> None:
    temporary = path.with_name(f".{path.name}.tmp.{os.getpid()}")
    try:
        temporary.write_text(content)
        os.replace(temporary, path)
    finally:
        try:
            temporary.unlink()
        except FileNotFoundError:
            pass


def _cache_is_current(cache: Path, listing: Path, graphs: Iterable[Path]) -> bool:
    if not cache.is_file() or cache.stat().st_size == 0:
        return False
    cache_mtime = cache.stat().st_mtime_ns
    return cache_mtime >= listing.stat().st_mtime_ns and all(
        cache_mtime >= graph.stat().st_mtime_ns for graph in graphs
    )


def ensure_graph_list_cache(
    graph_folder: Path,
    *,
    jobs: int,
    graph_fastas: Optional[Iterable[Path]] = None,
    force_rebuild: bool = False,
    script_dir: Optional[Path] = None,
) -> Tuple[Path, Path, bool]:
    """Return the graph list, template-only index, and whether it was rebuilt.

    Graph alignment uses PARTITION.FA; hotspot discovery exclusively uses
    PARTITION.fasta through a separate .template.list and .template.list.bin.
    """
    if jobs < 1:
        raise ValueError("KmerSearcher jobs must be positive")
    folder = graph_folder.expanduser().resolve()
    listing = graph_list_path(folder)
    lock_path = Path(str(listing) + ".lock")
    lock_path.parent.mkdir(parents=True, exist_ok=True)

    with lock_path.open("a+") as lock:
        fcntl.flock(lock.fileno(), fcntl.LOCK_EX)
        saved_graphs = (
            _read_valid_existing_list(listing, folder)
            if graph_fastas is None else None
        )
        graphs = (
            saved_graphs
            if saved_graphs is not None
            else _validated_graph_fastas(folder, graph_fastas)
        )
        content = "".join(f"{graph}\n" for graph in graphs)
        try:
            existing = listing.read_text()
        except FileNotFoundError:
            existing = None
        regenerated = existing != content
        if regenerated:
            _atomic_write(listing, content)

        _templates, cache, rebuilt = ensure_template_list_cache(
            folder, listing, jobs=jobs, force_rebuild=force_rebuild or regenerated,
            script_dir=script_dir,
        )
        return listing, cache, rebuilt


def ensure_template_list_cache(
    graph_folder: Path, graph_listing: Path, *, jobs: int,
    force_rebuild: bool = False, script_dir: Optional[Path] = None,
) -> Tuple[Path, Path, bool]:
    """Index original templates in exactly the graph list's partition order.

    Preserve an existing build-time template list with a different partition
    set/order. A separate selected template list serves that calling subset.
    Never substitute graph paths when original template FASTAs are missing.
    """
    from align_partition_hotspots import graph_name_from_list_entry

    if jobs < 1:
        raise ValueError("KmerSearcher jobs must be positive")
    folder = Path(graph_folder).expanduser().resolve()
    graph_listing = Path(graph_listing).expanduser().resolve()
    names = [
        graph_name_from_list_entry(line.split()[0], str(graph_listing))
        for line in graph_listing.read_text().splitlines()
        if line.strip() and not line.lstrip().startswith("#")
    ]
    if not names or len(names) != len(set(names)):
        raise ValueError(f"empty or duplicate graph partitions in {graph_listing}")
    if any(not re.fullmatch(r"[A-Za-z0-9_.#-]+", name) or name in {".", ".."}
           for name in names):
        raise ValueError(f"invalid graph partition name in {graph_listing}")
    templates = [folder / name / f"{name}.fasta" for name in names]
    missing = [path for path in templates if not path.is_file() or not path.stat().st_size]
    if missing:
        raise FileNotFoundError(
            "Hotspot discovery requires original PARTITION.fasta templates; "
            "reconstruct/resume the graph templates before calling. Missing: "
            + ", ".join(map(str, missing[:5]))
        )
    listing = template_list_path(graph_listing)
    content = "".join(f"{path}\n" for path in templates)
    if listing.is_file():
        old_paths = [Path(line) for line in listing.read_text().splitlines() if line.strip()]
        if (old_paths and old_paths != templates
                and all(path.is_absolute() and path.is_file() for path in old_paths)):
            listing = listing.with_name(listing.name.replace(".template.list", ".selected.template.list"))
    cache = Path(str(listing) + ".bin")
    marker = Path(str(cache) + ".fasta_target_unique_v2")
    incomplete = Path(str(cache) + ".tmp")
    lock_path = Path(str(listing) + ".lock")
    lock_path.parent.mkdir(parents=True, exist_ok=True)
    with lock_path.open("a+") as lock:
        fcntl.flock(lock.fileno(), fcntl.LOCK_EX)
        regenerated = not listing.is_file() or listing.read_text() != content
        if regenerated:
            _atomic_write(listing, content)
        rebuild = force_rebuild or regenerated or not _cache_is_current(
            cache, listing, templates,
        )
        if not rebuild:
            return listing, cache, False

        for path in (cache, marker, incomplete):
            try:
                path.unlink()
            except FileNotFoundError:
                pass

        tools = (script_dir or Path(__file__).resolve().parent).resolve()
        searcher = _local_kmer_searcher(tools)
        with tempfile.TemporaryDirectory(
            prefix=f".{listing.name}.kmersearcher.", dir=str(listing.parent),
        ) as temporary_name:
            temporary = Path(temporary_name)
            dummy_fasta = temporary / "index_dummy.fa"
            dummy_list = temporary / "index_dummy.list"
            dummy_fasta.write_text(">index_dummy\n" + "N" * 32 + "\n")
            dummy_list.write_text(f"index_dummy\t{dummy_fasta}\n")
            completed = subprocess.run(
                [
                    str(searcher), "-T", str(listing),
                    "-I", str(dummy_list), "-o", str(temporary) + os.sep,
                    "-c", "1000", "-w", "2000", "-n", str(jobs),
                ],
                check=False,
                stdout=subprocess.PIPE,
                stderr=subprocess.STDOUT,
                text=True,
            )
        if completed.returncode != 0:
            detail = "\n".join(completed.stdout.strip().splitlines()[-20:])
            message = (
                f"KmerSearcher could not build {cache}; status "
                f"{completed.returncode}"
            )
            if detail:
                message += f":\n{detail}"
            raise RuntimeError(message)
        if not _cache_is_current(cache, listing, templates):
            raise RuntimeError(
                f"KmerSearcher did not create a current nonempty cache: {cache}"
            )
        return listing, cache, True


def main(argv: Optional[Sequence[str]] = None) -> None:
    parser = argparse.ArgumentParser(
        description=(
            "create or validate <GraphFolder>.list and the separate "
            "<GraphFolder>.template.list[.bin] hotspot targets"
        )
    )
    parser.add_argument("-G", "--graph-folder", required=True)
    parser.add_argument("-j", "--jobs", type=int, default=1)
    parser.add_argument(
        "--script-folder", default=str(Path(__file__).resolve().parent),
        help="folder containing the repository KmerSearcher executable",
    )
    parser.add_argument("--force", action="store_true")
    args = parser.parse_args(argv)
    listing, cache, rebuilt = ensure_graph_list_cache(
        Path(args.graph_folder),
        jobs=args.jobs,
        force_rebuild=args.force,
        script_dir=Path(args.script_folder),
    )
    action = "rebuilt" if rebuilt else "current"
    print(f"[graph_list_cache] {action}: {listing}")
    print(f"[graph_list_cache] {action}: {cache}")


if __name__ == "__main__":
    main()
