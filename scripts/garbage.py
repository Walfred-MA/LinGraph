#!/usr/bin/env python3
"""Move finished temporary trees aside and delete them in the background.

Deleting a large work tree file by file on a network filesystem can stall a
stage for many minutes.  ``discard(path)`` replaces ``shutil.rmtree(path)``:
when ``LINGRAPH_GARBAGE`` names a garbage folder (``GRAPH/garbage`` or
``OUTPUT/garbage``, set by the launchers), the path is renamed into it, which
is one metadata operation, and a detached low-priority process deletes the
garbage folder's contents in parallel.  The original path is gone on return,
exactly as after ``rmtree``.

Without ``LINGRAPH_GARBAGE``, or when a rename is not possible (another
filesystem, a symlink, or a garbage folder inside the path), ``discard`` falls
back to ``shutil.rmtree`` with the same arguments.  Anything a killed job
leaves in a garbage folder is removed by the next ``sweep`` or deleter.

Command line (used for the detached deleter): ``garbage.py purge GARBAGE_DIR``.
"""
from __future__ import annotations

from concurrent.futures import ThreadPoolExecutor
import fcntl
import os
from pathlib import Path
import shutil
import subprocess
import sys
from typing import Optional, Union
import uuid

ENVIRONMENT_VARIABLE = "LINGRAPH_GARBAGE"
LOCK_NAME = ".purge.lock"
PURGE_THREADS = 8

PathLike = Union[str, "os.PathLike[str]"]


def garbage_root() -> Optional[Path]:
    value = os.environ.get(ENVIRONMENT_VARIABLE, "").strip()
    return Path(value).expanduser().absolute() if value else None


def use_garbage(root: PathLike) -> Path:
    """Send this process's (and its children's) discards to ``root``."""
    path = Path(root).expanduser().absolute()
    os.environ[ENVIRONMENT_VARIABLE] = str(path)
    return path


def discard(path: PathLike, ignore_errors: bool = False) -> None:
    """Remove ``path`` like ``shutil.rmtree``, deferring the deletion."""
    root = garbage_root()
    source = Path(path)
    if root is None or source.is_symlink() or not source.is_dir():
        # Missing paths, files and symlinks keep rmtree's exact behavior.
        shutil.rmtree(path, ignore_errors=ignore_errors)
        return
    absolute = source.absolute()
    try:
        inside = root == absolute or root.is_relative_to(absolute)
    except AttributeError:  # Python < 3.9 has no is_relative_to
        inside = str(root).startswith(str(absolute) + os.sep) or root == absolute
    if inside:
        shutil.rmtree(path, ignore_errors=ignore_errors)
        return
    try:
        root.mkdir(parents=True, exist_ok=True)
        target = root / f"{source.name}.{os.getpid()}.{uuid.uuid4().hex[:12]}"
        os.rename(source, target)
    except OSError:
        # Typically EXDEV: the garbage folder is on another filesystem.
        shutil.rmtree(path, ignore_errors=ignore_errors)
        return
    start_purger(root)


def _purger_running(root: Path) -> bool:
    try:
        with open(root / LOCK_NAME, "a") as handle:
            try:
                fcntl.flock(handle, fcntl.LOCK_EX | fcntl.LOCK_NB)
            except BlockingIOError:
                return True
            fcntl.flock(handle, fcntl.LOCK_UN)
            return False
    except OSError:
        # Locking unsupported here: start a deleter; extra ones are harmless.
        return False


def start_purger(root: PathLike) -> None:
    """Start a detached deleter for ``root`` unless one is already running."""
    root = Path(root)
    if _purger_running(root):
        return
    try:
        subprocess.Popen(
            [sys.executable, str(Path(__file__).resolve()), "purge", str(root)],
            stdin=subprocess.DEVNULL,
            stdout=subprocess.DEVNULL,
            stderr=subprocess.DEVNULL,
            close_fds=True,
            start_new_session=True,
        )
    except OSError:
        pass  # The next discard or sweep retries.


def sweep(root: PathLike) -> None:
    """Delete leftovers from earlier runs in the background, if any."""
    root = Path(root)
    if root.is_dir() and _entries(root):
        start_purger(root)


def _entries(root: Path):
    with os.scandir(root) as entries:
        return [entry for entry in entries if entry.name != LOCK_NAME]


def _remove(path: str) -> None:
    if os.path.isdir(path) and not os.path.islink(path):
        shutil.rmtree(path, ignore_errors=True)
    else:
        try:
            os.unlink(path)
        except OSError:
            pass


def _remove_tree(path: str, pool: ThreadPoolExecutor) -> None:
    """Delete one discarded item, spreading its top-level entries over threads."""
    if os.path.isdir(path) and not os.path.islink(path):
        try:
            children = [entry.path for entry in os.scandir(path)]
        except OSError:
            children = []
        list(pool.map(_remove, children))
    _remove(path)


def purge(root: PathLike) -> int:
    """Empty ``root`` while holding its lock; exit when nothing is left."""
    root = Path(root)
    try:
        os.nice(10)
    except (AttributeError, OSError):
        pass
    root.mkdir(parents=True, exist_ok=True)
    with open(root / LOCK_NAME, "a") as handle:
        try:
            fcntl.flock(handle, fcntl.LOCK_EX | fcntl.LOCK_NB)
        except BlockingIOError:
            return 0  # Another deleter owns this folder.
        except OSError:
            pass  # No locking on this filesystem; delete anyway.
        with ThreadPoolExecutor(max_workers=PURGE_THREADS) as pool:
            previous = None
            while True:
                names = {entry.name for entry in _entries(root)}
                # Stop when nothing is left or a pass removed nothing (e.g. a
                # file still held open on NFS); a later sweep retries those.
                if not names or names == previous:
                    break
                previous = names
                for name in names:
                    _remove_tree(str(root / name), pool)
    # A discard that saw our lock just before we released it did not start a
    # deleter; pick up anything that arrived meanwhile.
    if {entry.name for entry in _entries(root)} - (previous or set()):
        start_purger(root)
    return 0


def main(argv=None) -> int:
    argv = sys.argv[1:] if argv is None else list(argv)
    if len(argv) == 2 and argv[0] == "purge":
        return purge(argv[1])
    print("usage: garbage.py purge GARBAGE_DIR", file=sys.stderr)
    return 2


if __name__ == "__main__":
    raise SystemExit(main())
