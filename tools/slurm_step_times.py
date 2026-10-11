#!/usr/bin/env python
"""Run time of each step of a pipeline run, from its Slurm logs and sacct.

Job IDs come from the log names in RUN/slurm_logs (NAME-JOBID.out/.err, as
the pipelines' run-slurm writes them); sacct gives each job's state, run time,
CPUs and queue wait. One unit of work is one job name (one sample, batch or
once-step); a merge chromosome job (merge-cohort-chrom-N, merge-sv-chrom-N:
N is renumbered on every resume) is identified by the chromosomes its log
reports done instead. Of a unit run more than once, only its latest
COMPLETED job counts (failed, cancelled, timed-out and out-of-memory attempts
are skipped); a unit that never completed shows its latest job, flagged.

Merge steps count only the final merge run: merge-* jobs submitted at or
after the last COMPLETED merge-cohort-scan (merge-sv-scan for merge_grvcfs.py)
start; earlier merge runs (other code, other chromosome grouping) are
skipped. A COMPLETED chromosome job whose log reports no chromosome done only
confirmed finished work on a resume, and is skipped too.

Steps (job name without its sample, batch or chromosome) are listed in
pipeline order:

  step  units  completed  skipped_attempts  job_h  core_h  queue_h  longest_h  first_start  last_end  span_h

job_h = sum of the kept jobs' run times (no queue waiting), core_h = run time
x allocated CPUs, longest_h = the slowest unit, span_h = first start to last
end of the step (wall clock). One row per kept job goes to OUTPUT .jobs.tsv.

Steps are grouped into stages, each with a subtotal (change with
--group STEP=STAGE; unknown steps go to 'other'):
  1-preprocessing  prepare (LinGraph-prepare-*), match, extract-samples,
                   local-templates, lift-local-templates, alternative-loci
  2-alignment      kmersearch, minset, discover-novel-loci (graph build)
  3-linearization  lift, gcigar, gapfill, vcf
  4-merge          merge-* (final merge run only)
Steps the pipelines run on the launching node (not as Slurm jobs) have no
log here and are not counted.

Usage:
  python tools/slurm_step_times.py graph cohort_calls [-o pipeline_steps.tsv]
  (each RUN's slurm_logs/ and lingraph/slurm_logs/ are read)
"""
import argparse
from pathlib import Path
import re
import sys

sys.path.insert(0, str(Path(__file__).resolve().parent))
from merge_timecost import sacct_rows, seconds_between  # noqa: E402

LOG_NAME = re.compile(r'(.+)-(\d+)\.(?:out|err)')
CHROM_DONE = re.compile(r'chrom stage done: (\S+)')
CHROM_JOB = re.compile(r'merge-(?:cohort|sv|snp)-chrom(?:-\d+)?')
# Exact once-step names first: 'lift-local-templates' is not a 'lift-' sample.
ONCE = ('extract-samples', 'local-templates', 'lift-local-templates', 'alternative-loci',
        'merge-cohort-scan', 'merge-small-scan', 'merge-cohort-concat', 'merge-small-concat',
        'merge-cohort-publish', 'merge-sv-scan', 'merge-snp-prepare', 'merge-sv-concat',
        'merge-snp-concat', 'merge-publish')
PER_UNIT = (('LinGraph-prepare-', 'prepare'), ('kmersearch-', 'kmersearch'), ('minset-', 'minset'),
            ('match-', 'match'), ('lift-', 'lift'), ('gcigar-', 'gcigar'), ('gapfill-', 'gapfill'),
            ('vcf-', 'vcf'), ('merge-small-', 'merge-small-chrom'))
STAGES = {'1-preprocessing': ('prepare', 'match', 'extract-samples', 'local-templates',
                              'lift-local-templates', 'alternative-loci'),
          '2-alignment': ('kmersearch', 'minset', 'discover-novel-loci'),
          '3-linearization': ('lift', 'gcigar', 'gapfill', 'vcf')}
ORDER = ['prepare', 'kmersearch', 'discover-novel-loci', 'minset',
         'match', 'extract-samples', 'local-templates', 'lift-local-templates', 'alternative-loci',
         'lift', 'gcigar', 'gapfill', 'vcf', 'merge-cohort-scan', 'merge-sv-scan', 'merge-small-scan',
         'merge-snp-prepare', 'merge-cohort-chrom', 'merge-sv-chrom', 'merge-small-chrom',
         'merge-snp-chrom', 'merge-cohort-concat', 'merge-sv-concat', 'merge-small-concat',
         'merge-snp-concat', 'merge-cohort-publish', 'merge-publish']


def step_of(name):
    if name in ONCE:
        return name
    if CHROM_JOB.fullmatch(name):
        return re.sub(r'-\d+$', '', name)
    for prefix, step in PER_UNIT:
        if name.startswith(prefix):
            return step
    return re.sub(r'-\d+$', '', name)


def stage_of(step, groups):
    if step in groups:
        return groups[step]
    for stage, steps in STAGES.items():
        if step in steps:
            return stage
    return '4-merge' if step.startswith('merge-') else 'other'


def log_tail(logs, name, job_id, size=1 << 16):
    """End of a job's .err and .out (the chromosome lines come last)."""
    text = ''
    for suffix in ('err', 'out'):
        path = Path(logs[job_id]) / f'{name}-{job_id}.{suffix}'
        try:
            with open(path, 'rb') as handle:
                handle.seek(max(0, path.stat().st_size - size))
                text += handle.read().decode(errors='replace')
        except OSError:
            pass
    return text


def unit_of(row, logs):
    """Job name; a merge chromosome job: its chromosomes reported done."""
    name = row['JobName']
    if CHROM_JOB.fullmatch(name) and row['State'].startswith('COMPLETED'):
        chroms = sorted(set(CHROM_DONE.findall(log_tail(logs, name, row['JobID']))))
        if chroms:
            return f'{step_of(name)}:{",".join(chroms)}'
    return name


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('runs', nargs='+', metavar='RUN',
                        help='pipeline output folders (e.g. graph cohort_calls); their slurm_logs/ and '
                             'lingraph/slurm_logs/ are read, or give slurm_logs folders directly')
    parser.add_argument('--group', action='append', default=[], metavar='STEP=STAGE',
                        help='put STEP into STAGE (repeatable), e.g. --group discover-novel-loci=1-preprocessing')
    parser.add_argument('-o', '--output', help='also write the table here (kept jobs to OUTPUT .jobs.tsv)')
    args = parser.parse_args(argv)
    groups = dict(item.split('=', 1) for item in args.group)
    folders = []
    for run in map(Path, args.runs):
        found = [folder for folder in (run / 'slurm_logs', run / 'lingraph' / 'slurm_logs') if folder.is_dir()]
        folders += found or ([run] if run.is_dir() else [])
    names, logs = {}, {}
    for folder in folders:
        count = 0
        for path in folder.iterdir():
            match = LOG_NAME.fullmatch(path.name)
            if match:
                names[match.group(2)] = match.group(1)
                logs[match.group(2)] = folder
                count += 1
        print(f'[steps] {folder}: {count} log file(s)', file=sys.stderr, flush=True)
    if not names:
        raise ValueError(f'no NAME-JOBID.out/.err logs in {", ".join(map(str, folders)) or args.runs}')
    print(f'[steps] {len(names)} job(s); asking sacct', file=sys.stderr, flush=True)
    rows = sacct_rows(sorted(names, key=int))
    missing = len(names) - len(rows)
    if missing:
        print(f'[steps] sacct has no record of {missing} job(s) (too old or purged); they are not counted',
              file=sys.stderr, flush=True)

    kept, steps = [], {}
    new_entry = lambda: dict(units=0, completed=0, skipped=0, job_s=0.0, core_s=0.0,
                             queue_s=0.0, longest=0.0, first=None, last=None)
    # Only the final merge run: from the last completed merge scan on.
    scans = [row['Submit'] for row in rows if row['JobName'] in ('merge-cohort-scan', 'merge-sv-scan')
             and row['State'].startswith('COMPLETED')]
    final_merge = max(scans) if scans else None
    earlier = 0
    if final_merge:
        current = []
        for row in rows:
            if row['JobName'].startswith('merge-') and row['Submit'] < final_merge:
                steps.setdefault(step_of(row['JobName']), new_entry())['skipped'] += 1
                earlier += 1
            else:
                current.append(row)
        rows = current
        print(f'[steps] merge: final run from {final_merge}; {earlier} earlier merge job(s) skipped',
              file=sys.stderr, flush=True)

    # Per unit: the latest COMPLETED job, else the latest job.
    units = {}
    noop = 0
    for row in sorted(rows, key=lambda row: int(re.sub(r'\D.*', '', row['JobID']) or 0)):
        unit = unit_of(row, logs)
        if CHROM_JOB.fullmatch(row['JobName']) and row['State'].startswith('COMPLETED') and ':' not in unit:
            steps.setdefault(step_of(row['JobName']), new_entry())['skipped'] += 1
            noop += 1                        # a resume's "already complete" job
            continue
        units.setdefault(unit, []).append(row)
    if noop:
        print(f'[steps] {noop} completed chromosome job(s) reported no chromosome done; skipped',
              file=sys.stderr, flush=True)
    # A failed chromosome job's log names no chromosome: when its step has
    # completed jobs, it counts as a skipped attempt, not as its own unit.
    chrom_steps_done = {unit.split(':')[0] for unit in units if ':' in unit}
    for unit, attempts in units.items():
        done = [row for row in attempts if row['State'].startswith('COMPLETED')]
        step = step_of(attempts[-1]['JobName'])
        if not done and ':' not in unit and step in chrom_steps_done:
            steps.setdefault(step, new_entry())['skipped'] += len(attempts)
            continue
        row = (done or attempts)[-1]
        entry = steps.setdefault(step, new_entry())
        run = float(row['ElapsedRaw'] or 0)
        cpus = int(row['AllocCPUS'] or 0)
        queue = seconds_between(row['Submit'], row['Start'])
        entry['units'] += 1
        entry['completed'] += bool(done)
        entry['skipped'] += len(attempts) - 1
        entry['job_s'] += run
        entry['core_s'] += run * cpus
        entry['queue_s'] += queue
        entry['longest'] = max(entry['longest'], run)
        if row['Start'][:1].isdigit():
            entry['first'] = min(entry['first'] or row['Start'], row['Start'])
        if row['End'][:1].isdigit():
            entry['last'] = max(entry['last'] or row['End'], row['End'])
        kept.append((step, unit, row, run, cpus, queue, bool(done)))

    hours = lambda seconds: f'{seconds / 3600:.3f}'
    order = {name: index for index, name in enumerate(ORDER)}
    fields = ('units', 'completed', 'skipped', 'job_s', 'core_s', 'queue_s')
    lines = ['stage\tstep\tunits\tcompleted\tskipped_attempts\tjob_h\tcore_h\tqueue_h\tlongest_h'
             '\tfirst_start\tlast_end\tspan_h']

    def line(stage, step, entry, first=None, last=None, longest=None):
        span = seconds_between(first, last) if first and last else None
        return '\t'.join(map(str, [
            stage, step, entry['units'], entry['completed'], entry['skipped'], hours(entry['job_s']),
            hours(entry['core_s']), hours(entry['queue_s']), '' if longest is None else hours(longest),
            first or '', last or '', '' if span is None else hours(span)]))

    total = dict.fromkeys(fields, 0)
    by_stage = {}
    for step, entry in steps.items():
        if entry['units']:                   # else only skipped jobs (e.g. an old merge job name)
            by_stage.setdefault(stage_of(step, groups), []).append(step)
    for stage in sorted(by_stage):
        subtotal = dict.fromkeys(fields, 0)
        first = last = None
        for step in sorted(by_stage[stage], key=lambda name: (order.get(name, len(ORDER)), name)):
            entry = steps[step]
            for key in fields:
                subtotal[key] += entry[key]
            first = min(filter(None, [first, entry['first']]), default=None)
            last = max(filter(None, [last, entry['last']]), default=None)
            lines.append(line(stage, step, entry, entry['first'], entry['last'], entry['longest']))
        for key in fields:
            total[key] += subtotal[key]
        lines.append(line(stage, 'TOTAL', subtotal, first, last))
    lines.append(line('all', 'TOTAL', total))
    incomplete = total['units'] - total['completed']
    if incomplete:
        lines.append(f'# {incomplete} unit(s) never completed: their latest job is counted '
                     '(state in the jobs table)')
    if final_merge:
        lines.append(f'# merge steps: final merge run only (from {final_merge}); '
                     f'{earlier} earlier merge job(s) skipped')
    hidden = [step for step, entry in steps.items() if not entry['units']]
    other = by_stage.get('other')
    if other:
        lines.append(f'# steps in no stage (\'other\'; --group STEP=STAGE): {", ".join(sorted(other))}')
    if hidden:
        lines.append(f'# steps with only skipped jobs, not shown: {", ".join(sorted(hidden))}')
    text = '\n'.join(lines) + '\n'
    jobs_tsv = Path(args.output).with_suffix('.jobs.tsv') if args.output else folders[0] / 'steps.jobs.tsv'
    with open(jobs_tsv, 'w') as handle:
        handle.write('stage\tstep\tunit\tjob_id\tjob_name\tstate\tcpus\trun_s\tcore_s\tqueue_s\tsubmit\tstart\tend\n')
        for step, unit, row, run, cpus, queue, _done in sorted(
                kept, key=lambda item: (order.get(item[0], len(ORDER)), item[2]['Start'])):
            handle.write('\t'.join(map(str, [stage_of(step, groups), step, unit, row['JobID'], row['JobName'], row['State'], cpus,
                                              f'{run:.0f}', f'{run * cpus:.0f}', f'{queue:.0f}',
                                              row['Submit'], row['Start'], row['End']])) + '\n')
    if args.output:
        Path(args.output).write_text(text)
    sys.stdout.write(text)
    print(f'[steps] one row per kept job: {jobs_tsv}', file=sys.stderr)
    return 0


if __name__ == '__main__':
    try:
        raise SystemExit(main())
    except (ValueError, OSError) as error:
        print(f'[steps] ERROR: {error}', file=sys.stderr)
        raise SystemExit(1)
