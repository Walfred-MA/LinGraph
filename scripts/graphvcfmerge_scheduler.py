"""A bounded, killable alignment pool shared by nested-path coordinators.

Path coordinators are threads outside this pool. Each pool slot owns a reusable
Python subprocess/process group, so a cluster timeout can kill its running
native alignments without interrupting other clusters. No executors or native
thread pools are created inside an alignment worker.
"""
from __future__ import annotations

from collections import deque
from concurrent.futures import CancelledError, Future, TimeoutError as FutureTimeout
from contextlib import contextmanager
import importlib
import json
import math
import os
from pathlib import Path
import pickle
import queue
import signal
import sys
import tempfile
import threading
import time
import traceback


class ClusterTimeout(TimeoutError):
    """The cluster's alignment-first wall-clock budget expired."""


def payload_bytes(value):
    """Conservative input charge; references to disk-spooled data stay small."""
    if isinstance(value, (str, bytes)):
        return len(value) + 64
    if isinstance(value, (tuple, list)):
        return 64 + sum(payload_bytes(item) for item in value)
    if isinstance(value, dict):
        return 64 + sum(payload_bytes(k) + payload_bytes(v) for k, v in value.items())
    return 32


class AlignmentGroup:
    def __init__(self, pool, path, seconds):
        if seconds is not None and (not math.isfinite(seconds) or seconds <= 0):
            raise ValueError('cluster timeout must be finite and positive')
        self.pool = pool
        self.path = path
        self.started = time.monotonic()
        self.deadline = None if seconds is None else self.started + seconds
        self.seconds = seconds
        self.cancelled = False
        self.active = 0
        self.futures = set()

    def check(self):
        if self.pool.closed or self.cancelled:
            raise CancelledError('alignment group cancelled')
        if self.deadline is not None and time.monotonic() >= self.deadline:
            raise ClusterTimeout(f'nested cluster exceeded {self.seconds:g} seconds')

    def result(self, future):
        report = time.monotonic() + 30
        while True:
            self.check()
            try:
                return future.result(timeout=.1)
            except FutureTimeout:
                # A completed task may itself have raised TimeoutError.
                if future.done():
                    return future.result()
                if time.monotonic() >= report:
                    stats = self.pool.stats()
                    print(f'[merge:nested:waiting] path={self.path} '
                          f'elapsed={time.monotonic() - self.started:.0f}s '
                          f'limit={self.seconds if self.seconds is not None else "fallback"} '
                          f'shared_active={stats["active"]}/{self.pool.workers} '
                          f'shared_queued={stats["queued"]}', file=sys.stderr, flush=True)
                    report = time.monotonic() + 30

    def map(self, function, tasks):
        """Consume lazily, retaining at most 2N futures for this caller."""
        try:
            return self._map(function, tasks)
        except BaseException:
            # Callers may remove sequence/score files in their own finally
            # blocks. Stop readers before unwinding into that cleanup.
            self.close()
            raise

    def _map(self, function, tasks):
        pending = deque()
        items = iter(tasks)
        exhausted = False
        results = []
        while pending or not exhausted:
            while not exhausted and len(pending) < self.pool.max_pending:
                try:
                    task = next(items)
                except StopIteration:
                    exhausted = True
                    break
                pending.append(self.pool.submit(self, function, task))
            if pending:
                results.append(self.result(pending.popleft()))
        return results

    def align(self, query, reference):
        return self.result(self.pool.submit(self, pair_alignment_task, (query, reference)))

    def align_pairs(self, pairs):
        return self.map(pair_alignment_task, pairs)

    def close(self):
        # Queued tasks can be cancelled without touching their input files.
        # Wait only for running tasks; unrelated long tasks must not block us.
        with self.pool.condition:
            self.cancelled = True
            for future in tuple(self.futures):
                future.cancel()
            self.pool.condition.notify_all()
            while self.active:
                self.pool.condition.wait(.1)


class AlignmentScheduler:
    def __init__(self, workers, directory=None, max_pending_bytes=64 * 1024 * 1024):
        self.workers = int(workers)
        if self.workers < 1 or max_pending_bytes < 1:
            raise ValueError('alignment workers and queue byte budget must be positive')
        self.max_pending = 2 * self.workers
        self.max_pending_bytes = int(max_pending_bytes)
        self.condition = threading.Condition()
        self.queue = queue.Queue()
        self.closed = False
        self.pending = self.pending_bytes = self.active = self.completed = 0
        self.peak_active = self.peak_pending = self.peak_pending_bytes = 0
        self.pids = set()
        self.temporary = tempfile.TemporaryDirectory(prefix='shared-alignments-', dir=directory)
        self.root = Path(self.temporary.name)
        self.threads = [threading.Thread(target=self._slot, name=f'alignment-slot-{i}',
                                         daemon=True) for i in range(self.workers)]
        self.previous_sigterm = None
        for thread in self.threads:
            thread.start()

    def __enter__(self):
        if threading.current_thread() is threading.main_thread():
            self.previous_sigterm = signal.getsignal(signal.SIGTERM)
            signal.signal(signal.SIGTERM, self._terminated)
        return self

    @staticmethod
    def _terminated(signum, _frame):
        raise SystemExit(128 + signum)

    def __exit__(self, *_exc):
        try:
            self.close()
        finally:
            if self.previous_sigterm is not None:
                signal.signal(signal.SIGTERM, self.previous_sigterm)

    @contextmanager
    def group(self, path, seconds=None):
        group = AlignmentGroup(self, path, seconds)
        try:
            yield group
        finally:
            group.close()

    def submit(self, group, function, payload):
        module, name = function.__module__, function.__name__
        if module not in ('graphvcfmerge_scheduler', 'graphvcfmerge_nested',
                          'graphvcfmerge_kmer', 'graphvcfmerge'):
            raise ValueError(f'unsupported alignment task {module}.{name}')
        weight = payload_bytes(payload)
        with self.condition:
            while True:
                group.check()
                # A single unusually large pair can make progress alone.
                if self.pending < self.max_pending and (
                    not self.pending or self.pending_bytes + weight <= self.max_pending_bytes
                ):
                    break
                self.condition.wait(.1)
            future = Future()
            group.futures.add(future)
            self.pending += 1
            self.pending_bytes += weight
            self.peak_pending = max(self.peak_pending, self.pending)
            self.peak_pending_bytes = max(self.peak_pending_bytes, self.pending_bytes)
            self.queue.put((group, future, (module, name, payload), weight))
            return future

    def stats(self):
        with self.condition:
            return dict(active=self.active, queued=self.pending - self.active,
                        completed=self.completed, peak_active=self.peak_active,
                        peak_pending=self.peak_pending, peak_pending_bytes=self.peak_pending_bytes,
                        worker_pids=sorted(self.pids))

    def _slot(self):
        from graphvcfmerge_nested import AlignmentSession
        session = None
        try:
            while True:
                entry = self.queue.get()
                if entry is None:
                    return
                group, future, payload, weight = entry
                running = False
                try:
                    with self.condition:
                        if not future.set_running_or_notify_cancel():
                            continue
                        group.check()
                        self.active += 1
                        group.active += 1
                        running = True
                        self.peak_active = max(self.peak_active, self.active)
                    if session is None:
                        environment = dict(os.environ)
                        for key in ('OMP_NUM_THREADS', 'OPENBLAS_NUM_THREADS', 'MKL_NUM_THREADS',
                                    'NUMEXPR_NUM_THREADS', 'VECLIB_MAXIMUM_THREADS'):
                            environment[key] = '1'
                        session = AlignmentSession(
                            [sys.executable, str(Path(__file__).resolve()), '--worker'],
                            environment=environment, directory=str(self.root),
                        )
                        with self.condition:
                            self.pids.add(session.process.pid)
                    result = session.run(payload, float('inf'), check_cancelled=group.check)
                    future.set_result(result)
                except BaseException as error:
                    if session is not None:
                        session.close()
                        session = None
                    if not future.done():
                        future.set_exception(error)
                finally:
                    with self.condition:
                        if running:
                            self.active -= 1
                            group.active -= 1
                        self.pending -= 1
                        self.pending_bytes -= weight
                        self.completed += 1
                        group.futures.discard(future)
                        self.condition.notify_all()
                    # Do not retain sequences/results while waiting for more work.
                    entry = payload = future = group = None
                    if 'result' in locals():
                        del result
        finally:
            if session is not None:
                session.close()

    def close(self):
        with self.condition:
            if not self.closed:
                self.closed = True
                self.condition.notify_all()
                for _ in self.threads:
                    self.queue.put(None)
        for thread in self.threads:
            thread.join()
        self.temporary.cleanup()


_pair_cache = {}


def pair_alignment_task(payload):
    import graphvcfmerge as merger
    import graphreftovcf as core
    query, reference = payload
    # Reuse a small number of mappy indexes without retaining every template
    # seen during a chromosome. Cache keys are sequence contents, never IDs.
    if reference not in _pair_cache and (
        len(_pair_cache) >= 8 or sum(map(len, _pair_cache)) + len(reference) > 64 * 1024 * 1024
    ):
        _pair_cache.clear()
    return merger._fast_pair_alignment(query, reference, _pair_cache, core)


def worker():
    """AlignmentSession wire protocol; import original functions in a clean process."""
    acknowledgements = sys.stdout
    sys.stdout = sys.stderr
    # Import once before the first request; native thread limits are set at exec.
    import graphvcfmerge as merger
    import graphreftovcf
    del graphreftovcf
    for command in sys.stdin.buffer:
        fields = command.decode().rstrip('\n').split('\t')
        if len(fields) != 2 or fields[0] != 'run':
            raise ValueError('invalid shared alignment worker request')
        root = Path(json.loads(fields[1]))
        try:
            with (root / 'request.pkl').open('rb') as source:
                module, name, payload = pickle.load(source)
            result = (True, getattr(importlib.import_module(module), name)(payload))
        except Exception:
            result = (False, traceback.format_exc())
        with (root / 'result.pkl').open('wb') as out:
            pickle.dump(result, out, protocol=pickle.HIGHEST_PROTOCOL)
        del result
        if 'payload' in locals():
            del payload
        merger._truvari_seqsim_cached.cache_clear()
        acknowledgements.write('done\n')
        acknowledgements.flush()


if __name__ == '__main__':
    if sys.argv[1:] != ['--worker']:
        raise SystemExit('Internal worker; run graphvcfmerge.py')
    worker()
