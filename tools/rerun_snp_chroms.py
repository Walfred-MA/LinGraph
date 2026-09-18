#!/usr/bin/env python3
"""Submit every saved SNP chromosome again as a separate Slurm job.

Existing parts never cause a job to be skipped. Shards and existing complete
parts are retained; each merger atomically replaces its part after success.
Stop previous controllers/jobs first. Active chromosome locks are respected.
This launcher submits chromosome jobs only, without SV or insertion dependencies.
"""
from __future__ import annotations

import argparse
import fcntl
import hashlib
import importlib.util
import json
from pathlib import Path
import re
import shlex
import shutil
import subprocess
import sys
import tempfile


ROOT = Path(__file__).resolve().parents[1] / "scripts"


def preflight(work):
    shards = work / 'merge/tmp/snp_shards'
    manifest = json.loads((shards / 'manifest.json').read_text())
    if manifest.get('protocol') not in (2, 3, 4):
        raise ValueError('expected a compact SNP manifest (protocol 2, 3 or 4)')
    keys = list(manifest.get('chrom_order', manifest['chroms']))
    if len(keys) != len(set(keys)) or set(keys) != set(manifest['chroms']):
        raise ValueError('chrom_order does not contain each chromosome exactly once')
    if any(not re.fullmatch(r'[0-9a-f]{64}', key) for key in keys):
        raise ValueError('invalid chromosome key in SNP manifest')
    directory = Path(manifest['directory'])
    if not directory.is_absolute() or not directory.is_dir():
        raise ValueError(f'saved SNP generation directory is unavailable: {directory}')
    # A stale lock file is harmless. Only a lock held by a live process prevents
    # submission; never unlink it and allow two writers into the same part.
    for key in keys:
        locus = 'locus_' + hashlib.sha256(key.encode()).hexdigest()[:24]
        lock = directory / 'chroms' / locus / 'merge.lock'
        if lock.exists():
            with lock.open('r') as handle:
                try:
                    fcntl.flock(handle, fcntl.LOCK_EX | fcntl.LOCK_NB)
                except BlockingIOError as error:
                    raise RuntimeError(f'{manifest["chroms"][key]} still has an active merge worker; '
                                       f'stop its old job before rerunning ({lock})') from error
                fcntl.flock(handle, fcntl.LOCK_UN)
    return shards, manifest, keys


def submit(args):
    if args.jobs < 1 or args.processes < 1:
        raise ValueError('--jobs and --processes must be positive')
    work = args.work.expanduser().resolve()
    shards, manifest, keys = preflight(work)
    if not keys:
        print('[snp-rerun] No SNP chromosomes in the manifest.', flush=True)
        return []
    script = ROOT / 'graphvcfmerge_snp.py'
    if not script.is_file() or not (ROOT / 'graphvcfmerge_snp_memory.py').is_file():
        raise ValueError('sync the current graphvcfmerge_snp*.py files before submitting')
    for module in ('numpy', 'polars'):
        if importlib.util.find_spec(module) is None:
            raise RuntimeError(f'{module} is missing from {sys.executable}; install with: '
                               f'{shlex.join([sys.executable, "-m", "pip", "install", module])}')
    sbatch = shutil.which(args.sbatch)
    if not sbatch:
        raise FileNotFoundError(f'sbatch command not found: {args.sbatch}')
    logs = Path(tempfile.mkdtemp(prefix='snp-rerun-', dir=work))
    print(f'[snp-rerun] Rerunning all {len(keys)} chromosomes, including existing parts; '
          f'{args.processes} CPUs, {args.memory} per job, at most {args.jobs} concurrent jobs.', flush=True)
    print(f'[snp-rerun] Logs and submitted job IDs: {logs}', flush=True)
    ids = []
    with (logs / 'jobs.jsonl').open('w') as ledger:
        for index, key in enumerate(keys):
            worker = [sys.executable, '-u', str(script), '--stage', 'chrom',
                      '--shards-dir', str(shards), '--chrom', key, '--processes', str(args.processes)]
            command = [sbatch, '--parsable', '--export=ALL', '--ntasks=1',
                       f'--job-name=merge-snp-{key[:12]}', f'--cpus-per-task={args.processes}',
                       f'--mem={args.memory}', f'--time={args.time}',
                       f'--output={logs}/{key}-%j.out', f'--error={logs}/{key}-%j.err']
            for option in ('account', 'partition'):
                if getattr(args, option):
                    command.append(f'--{option}={getattr(args, option)}')
            # Separate sbatch jobs, arranged into independent chains to bound
            # concurrency without needing a persistent launcher or job arrays.
            # Failure of one chromosome must not prevent the others from running.
            if index >= args.jobs:
                command.append(f'--dependency=afterany:{ids[index - args.jobs]}')
            command.extend(['--wrap', shlex.join(worker)])
            result = subprocess.run(command, capture_output=True, text=True)
            if result.returncode:
                raise RuntimeError(f'sbatch failed for {manifest["chroms"][key]}: '
                                   f'{result.stderr.strip() or result.stdout.strip()}; '
                                   f'{len(ids)} jobs already submitted (see {logs}/jobs.jsonl)')
            response = result.stdout.strip()
            match = re.fullmatch(r'(\d+)(?:;[^\s;]+)?', response)
            if not match:
                raise RuntimeError(f'uncertain sbatch response {response!r}; inspect squeue before retrying; '
                                   f'see {logs}/jobs.jsonl for earlier submissions')
            job_id = match.group(1)
            ids.append(job_id)
            ledger.write(json.dumps(dict(job_id=job_id, key=key, chrom=manifest['chroms'][key],
                                         command=command)) + '\n')
            ledger.flush()
            print(f'[snp-rerun] submitted {manifest["chroms"][key]}: job {job_id}', flush=True)
            if result.stderr.strip():
                print(result.stderr.strip(), file=sys.stderr, flush=True)
    print('[snp-rerun] All chromosome jobs submitted. Run SNP concat after they all succeed.', flush=True)
    return ids


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--work', type=Path, required=True, help='existing unfinished-graphcigar directory')
    parser.add_argument('-t', '--processes', type=int, default=32)
    parser.add_argument('--memory', default='64G')
    parser.add_argument('--jobs', type=int, default=40, help='maximum concurrent chromosome jobs [40]')
    parser.add_argument('--time', default='5:00:00')
    parser.add_argument('--account', default='')
    parser.add_argument('--partition', default='')
    parser.add_argument('--sbatch', default='sbatch')
    args = parser.parse_args(argv)
    try:
        submit(args)
    except (OSError, ValueError, KeyError, RuntimeError) as error:
        parser.exit(1, f'[snp-rerun] ERROR: {error}\n')
    return 0


if __name__ == '__main__':
    raise SystemExit(main())
