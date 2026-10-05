"""Keep Snakemake run state (``.snakemake/``) with each run, not in the code.

Both workflow launchers pass ``--directory RUN_DIR`` (the graph folder or the
cohort output folder), so locks, incomplete-job markers, metadata and logs live
in ``RUN_DIR/.snakemake``.  Separate runs from one checkout then lock and resume
independently, and unlocking one run cannot release another.
"""
from __future__ import annotations

from pathlib import Path
import shutil
from typing import List


def state_dir(run_dir: Path) -> Path:
    return Path(run_dir) / ".snakemake"


def check_no_legacy_incomplete(workflow_dir: Path) -> None:
    """Refuse a new run while the old in-tree state has unfinished jobs.

    Before run state moved to the run folder, ``WORKFLOW/.snakemake/incomplete``
    recorded jobs that had started but not finished.  The new location does not
    see those markers, so resuming such a run here could accept half-written
    outputs as complete.
    """
    incomplete = state_dir(workflow_dir) / "incomplete"
    if incomplete.is_dir() and any(incomplete.iterdir()):
        raise RuntimeError(
            f"{incomplete} holds unfinished jobs from a run started before "
            "Snakemake state moved to the run folder. Finish or clean that run "
            "with the scripts that started it; if it is abandoned, remove "
            f"{state_dir(workflow_dir)} and start again."
        )


def unlock(run_dir: Path) -> List[str]:
    """Remove only this run's Snakemake locks; return what was removed."""
    locks = state_dir(run_dir) / "locks"
    if not locks.exists():
        return []
    shutil.rmtree(locks)
    return [str(locks)]
