#!/usr/bin/env python3
"""Nested alignment refinement with a killable two-hour cluster deadline."""
from __future__ import annotations

import atexit
import json
import math
from multiprocessing.util import Finalize
import os
from pathlib import Path
import pickle
import select
import signal
import subprocess
import sys
import tempfile
import time
import traceback

NESTED_CLUSTER_SECONDS = 2 * 60 * 60
NESTED_PROGRESS_SECONDS = 30


class AlignmentSession:
    """Reusable interpreter in its own process group, including native children."""
    def __init__(self, command=None, environment=None, directory=None):
        self.temp = tempfile.TemporaryDirectory(prefix="merge-nested-", dir=directory)
        self.root = Path(self.temp.name)
        self.closed = False
        try:
            self.process = subprocess.Popen(
                command or [sys.executable, str(Path(__file__).resolve()), "--worker", str(self.root)],
                stdin=subprocess.PIPE, stdout=subprocess.PIPE, start_new_session=True,
                env=environment,
            )
        except BaseException:
            self.temp.cleanup()
            raise

    def close(self):
        if self.closed:
            return
        self.closed = True
        # Kill the group even if its leader has exited: an external aligner
        # or inner pool worker may still be running in that group.
        try:
            os.killpg(self.process.pid, signal.SIGKILL)
        except ProcessLookupError:
            pass
        self.process.wait()
        self.process.stdin.close()
        self.process.stdout.close()
        self.temp.cleanup()

    def run(self, payload, seconds, check_cancelled=None):
        started = time.monotonic()
        report = started + NESTED_PROGRESS_SECONDS
        if check_cancelled is not None:
            check_cancelled()
        path_id = (
            payload[3].get("nested_path_id", "unknown")
            if isinstance(payload, tuple) and len(payload) == 5 else "unknown"
        )
        with (self.root / "request.pkl").open("wb") as out:
            pickle.dump(payload, out, protocol=pickle.HIGHEST_PROTOCOL)
        message = ('run\t' + json.dumps(str(self.root)) + '\n').encode() if check_cancelled else b'run\n'
        self.process.stdin.write(message)
        self.process.stdin.flush()
        acknowledgement = b""
        while acknowledgement != b"done\n":
            if check_cancelled is not None:
                check_cancelled()
            ready, _, _ = select.select(
                [self.process.stdout], [], [],
                min(.1 if check_cancelled else NESTED_PROGRESS_SECONDS,
                    max(0, seconds - (time.monotonic() - started))),
            )
            if not ready:
                elapsed = time.monotonic() - started
                if elapsed >= seconds:
                    self.close()
                    raise TimeoutError(f"nested cluster exceeded {seconds:g} seconds")
                if check_cancelled is not None or time.monotonic() < report:
                    continue
                print(
                    f"[merge:nested:waiting] path={path_id} alignment_worker={self.process.pid} "
                    f"elapsed={elapsed:.0f}s limit={seconds:g}s",
                    file=sys.stderr, flush=True,
                )
                report = time.monotonic() + NESTED_PROGRESS_SECONDS
                continue
            chunk = os.read(self.process.stdout.fileno(), 5 - len(acknowledgement))
            if not chunk:
                raise RuntimeError("nested alignment worker exited without a result")
            acknowledgement += chunk
            if not b"done\n".startswith(acknowledgement):
                raise RuntimeError("invalid nested alignment worker acknowledgement")
        with (self.root / "result.pkl").open("rb") as source:
            success, result = pickle.load(source)
        (self.root / "result.pkl").unlink()
        if not success:
            raise RuntimeError(f"nested alignment failed:\n{result}")
        return result


_session = None
_owner = None
_finalizer = None


def close_session():
    global _session, _owner
    if _session is not None and _owner == os.getpid():
        _session.close()
    _session, _owner = None, None


atexit.register(close_session)


def refine_nested(members, sequences, row_ids, context, pooled_map=None, workers=1,
                  alignment_scheduler=None):
    """Try alignment first; only expiry, not no-hit/error, triggers KmerMatch."""
    global _session, _owner, _finalizer
    if len(members) == 1:
        return [(members[0], [members[0] + (None,)])]
    if alignment_scheduler is not None:
        return _refine_shared(members, sequences, row_ids, context, alignment_scheduler)
    seconds = float(context.get("nested_cluster_timeout", NESTED_CLUSTER_SECONDS))
    if not math.isfinite(seconds) or seconds <= 0:
        raise ValueError("nested cluster timeout must be finite and positive")
    if _session is None or _owner != os.getpid():
        # A forked process must not signal a session belonging to its parent.
        _session = AlignmentSession()
        _owner = os.getpid()
        _finalizer = Finalize(None, close_session, exitpriority=10)
    # Do not serialize the cohort manifest for each small derived cluster.
    alignment_context = {key: context[key] for key in
                         ("size_similarity", "sequence_similarity", "winnowmap", "meryl", "kmermatch", "nested_path_id")
                         if key in context}
    alignment_context["legacy_fast"] = True
    alignment_context["records_dir"] = str(_session.root)
    loaded_sequences = list(sequences)
    if len(loaded_sequences) != len(members):
        raise ValueError("nested sequence/member count mismatch")
    try:
        return _session.run(
            (members, loaded_sequences, row_ids, alignment_context, workers),
            seconds,
        )
    except TimeoutError:
        close_session()
        print(f"[merge:nested-fallback] path={context.get('nested_path_id', 'unknown')} "
              f"cluster at {members[0][0]} "
              f"({len(members)} members) exceeded {seconds:g}s; alignment process group "
              f"killed; restarting this cluster with KmerMatch >=0.5 and mappy ({workers} worker(s))",
              file=sys.stderr, flush=True)
        from graphvcfmerge_kmer import refine_kmer_first
        return refine_kmer_first(
            members, loaded_sequences, row_ids, context, pooled_map, workers,
        )
    except BaseException:
        close_session()
        raise


def _bodies_task(payload):
    import graphvcfmerge as merger
    import graphreftovcf as core
    group_number, representative_sequence, candidates, context = payload
    sequences = [representative_sequence] + [seq for _, seq in candidates]
    bodies = merger._final_member_bodies(list(range(len(sequences))), 0, sequences, core, context)
    return group_number, [(index, bodies[slot + 1]) for slot, (index, _) in enumerate(candidates)]


def _refine_shared(members, sequences, row_ids, context, scheduler):
    """Ordered grouping with alignments/fallbacks using the same global slots."""
    import graphvcfmerge as merger
    import graphreftovcf as core
    from graphvcfmerge_scheduler import ClusterTimeout
    seconds = float(context.get('nested_cluster_timeout', NESTED_CLUSTER_SECONDS))
    path = context.get('nested_path_id', 'unknown')
    sequences = list(sequences)
    if len(sequences) != len(members):
        raise ValueError('nested sequence/member count mismatch')
    alignment_context = {key: context[key] for key in
                         ('size_similarity', 'sequence_similarity', 'winnowmap', 'meryl',
                          'kmermatch', 'nested_path_id') if key in context}
    alignment_context.update(legacy_fast=True, records_dir=str(scheduler.root))
    try:
        with scheduler.group(path, seconds) as service:
            groups = merger._nearest_online_subclusters(
                members, sequences, row_ids, core, alignment_context,
                alignment_service=service,
            )
            representatives = [merger._current_representative_index(group, members, row_ids)
                               for group in groups]
            bodies = [{rep: None} for rep in representatives]
            total = sum(len(group) - 1 for group in groups)
            chunk = max(1, min(64, -(-total // (scheduler.workers * 4))))

            def tasks():
                for number, (group, rep) in enumerate(zip(groups, representatives)):
                    batch, bases = [], len(sequences[rep])
                    for index in group:
                        if index != rep:
                            sequence = sequences[index]
                            if batch and (len(batch) >= chunk or bases + len(sequence) > 1024 * 1024):
                                yield number, sequences[rep], batch, alignment_context
                                batch, bases = [], len(sequences[rep])
                            batch.append((index, sequence))
                            bases += len(sequence)
                    if batch:
                        yield number, sequences[rep], batch, alignment_context

            for number, results in service.map(_bodies_task, tasks()):
                bodies[number].update(results)
            service.check()
            return [(members[rep], sorted(
                (members[i] + (bodies[number][i],) for i in group),
                key=lambda member: (member[0], member[1], member[2], member[4]),
            )) for number, (group, rep) in enumerate(zip(groups, representatives))]
    except ClusterTimeout:
        # Exiting the group has already killed/reaped only its active tasks.
        print(f'[merge:nested-fallback] path={path} cluster at {members[0][0]} '
              f'({len(members)} members) exceeded {seconds:g}s; cancelled its alignments; '
              f'restarting with KmerMatch on {scheduler.workers} shared worker(s)',
              file=sys.stderr, flush=True)
        from graphvcfmerge_kmer import refine_kmer_first
        with scheduler.group(path + ':fallback') as service:
            return refine_kmer_first(members, sequences, row_ids, context,
                                     service.map, scheduler.workers)


def _aligned_groups(payload):
    import graphvcfmerge as merger
    import graphreftovcf as core
    members, sequences, row_ids, context, workers = payload
    groups = merger._nearest_online_subclusters(members, sequences, row_ids, core, context)
    representatives = [merger._current_representative_index(group, members, row_ids) for group in groups]
    bodies = [{rep: None} for rep in representatives]
    total = sum(len(group) - 1 for group in groups)
    chunk = max(1, min(64, -(-total // (max(1, workers) * 4))))
    tasks = []
    for number, (group, rep) in enumerate(zip(groups, representatives)):
        pending = [index for index in group if index != rep]
        for start in range(0, len(pending), chunk):
            tasks.append((number, sequences[rep], [(i, sequences[i]) for i in pending[start:start + chunk]], context))
    if workers > 1 and len(tasks) > 1:
        # Legacy standalone callers may request an inner pool. Main multi-core
        # nested rounds use _refine_shared and never enter this branch.
        with merger.mp_context().Pool(min(workers, len(tasks))) as pool:
            for number, results in pool.imap_unordered(_bodies_task, tasks, chunksize=1):
                bodies[number].update(results)
    else:
        for task in tasks:
            number, results = _bodies_task(task)
            bodies[number].update(results)
    refined = []
    for number, (group, rep) in enumerate(zip(groups, representatives)):
        chosen = sorted((members[i] + (bodies[number][i],) for i in group),
                        key=lambda member: (member[0], member[1], member[2], member[4]))
        refined.append((members[rep], chosen))
    merger._truvari_seqsim_cached.cache_clear()
    return refined


def worker(root):
    root = Path(root)
    acknowledgements = sys.stdout
    sys.stdout = sys.stderr
    for command in sys.stdin.buffer:
        if command != b"run\n":
            raise ValueError("invalid nested worker command")
        try:
            with (root / "request.pkl").open("rb") as source:
                payload = pickle.load(source)
            result = (True, _aligned_groups(payload))
        except Exception:
            result = (False, traceback.format_exc())
        with (root / "result.pkl").open("wb") as out:
            pickle.dump(result, out, protocol=pickle.HIGHEST_PROTOCOL)
        del result
        if "payload" in locals():
            del payload
        acknowledgements.write("done\n")
        acknowledgements.flush()


if __name__ == "__main__":
    if len(sys.argv) != 3 or sys.argv[1] != "--worker":
        raise SystemExit("Internal worker; run graphvcfmerge.py or recover_graphvcfmerge.py")
    worker(sys.argv[2])
