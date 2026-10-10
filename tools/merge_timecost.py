#!/usr/bin/env python
"""Time cost of cohort merges at growing cohort sizes, on Slurm.

Each experiment picks the per-sample grVCFs matching a regular expression
(as grep -E) in their path or in their samples' lines of the cohort's
query_paths.normalized.txt (name and assembly FASTA), and merges them with
tools/merge_grvcfs.py --slurm (its default job resources). Experiments run one
after another; this driver waits for each merge and records:

  TIMECOST/NAME/mergevcf.list    the selected grVCFs
  TIMECOST/NAME/query_paths.txt  their assemblies (--mode exact)
  TIMECOST/NAME/merge.log        merge messages, each prefixed with elapsed time
  TIMECOST/NAME/phases.tsv       wall time of each merge phase
  TIMECOST/NAME/jobs.tsv         each Slurm job's run time, CPUs and queue wait (sacct)
  TIMECOST/NAME/timecost.json    totals
  TIMECOST/NAME/merge/           merge output (Slurm logs in merge/slurm_logs)
  TIMECOST/summary.tsv           one line per finished experiment

Wall times run from start to end of the merge, Slurm queue waiting included.
Added-up times leave the waiting out: job_h = sum of the jobs' run times
(sacct Elapsed), core_h = sum of run time x allocated CPUs, queue_h = sum of
the jobs' waits (Start - Submit). Phases: plan, sv-scan, snp-prepare, chrom
(all SV and SNP chromosome jobs), sv-concat, snp-concat, publish (until the
merge exits).

Default experiments (name, pattern, expected haplotypes):
  trio NA192 6 | panarabic KSA|Japan 28 | hgsvc3 hgsvc3 122 | hprc2 HPRC2 464 | all . 1152

Examples:
  python tools/merge_timecost.py start -C cohort_calls -O timecost -t 32 \\
      --slurm-account mchaisso_100 --slurm-partition qcb
  python tools/merge_timecost.py start ... --only trio panarabic
  python tools/merge_timecost.py summary -O timecost [--refresh]
"""
import argparse
from datetime import datetime
import json
from pathlib import Path
import re
import shlex
import shutil
import subprocess
import sys
import time

ROOT = Path(__file__).resolve().parents[1]
SCRIPTS = ROOT / 'scripts' if (ROOT / 'scripts' / 'cohort_vcf_merge.py').is_file() else ROOT
MERGE_TOOL = Path(__file__).resolve().with_name('merge_grvcfs.py')
EXPERIMENTS = (('trio', 'NA192', 6), ('panarabic', 'KSA|Japan', 28), ('hgsvc3', 'hgsvc3', 122),
               ('hprc2', 'HPRC2', 464), ('all', '.', 1152))
# merge_grvcfs.py prints "[merge_grvcfs] <stage>[ n/N chrom]: <command>" as it
# starts each stage (or submits its Slurm job).
STAGE_LINE = re.compile(r'^\[merge_grvcfs\] ([a-z-]+)(?: \d+/\d+ .*?)?: \S*python\S* ')
PHASE = {'sv-chrom': 'chrom', 'snp-chrom': 'chrom'}
# sbatch --parsable prints "JOBID" or "JOBID;CLUSTER" on its own line.
JOB_ID_LINE = re.compile(r'^\[\s*[\d.]+s\] (\d+)(?:;\S+)?$')


# ---------------------------------------------------------------- selection

def cohort_vcfs(cohort, listing):
    sys.path.insert(0, str(SCRIPTS))
    import graphvcfmerge
    listing = Path(listing) if listing else Path(cohort) / 'inputs' / 'mergevcf.list'
    if listing.is_file():
        return [str(Path(path).resolve()) for path in graphvcfmerge.read_vcf_input_lists([str(listing)])], listing
    paths = sorted(str(path.resolve()) for path in (Path(cohort) / 'samples').glob('*/*.vcf'))
    if not paths:
        raise ValueError(f'no {listing} and no {cohort}/samples/*/*.vcf')
    return paths, Path(cohort) / 'samples'


def vcf_samples(path):
    """Sample names of a grVCF's #CHROM line (the names --exact looks up)."""
    import gzip
    opener = gzip.open if str(path).endswith('.gz') else open
    with opener(path, 'rt') as handle:
        for line in handle:
            if line.startswith('#CHROM'):
                return line.rstrip('\n').split('\t')[9:]
            if not line.startswith('#'):
                break
    raise ValueError(f'{path}: no #CHROM header line')


def subset_query_paths(cohort, query_paths, names):
    """Lines (NAME FASTA [FAI], paths made absolute) of the selected samples
    only, so the merge does not see assemblies outside the experiment."""
    query_paths = Path(query_paths) if query_paths else Path(cohort) / 'inputs' / 'query_paths.normalized.txt'
    if not query_paths.is_file():
        raise ValueError(f'--mode exact needs {query_paths} (or --query-paths)')
    base = query_paths.resolve().parent
    wanted, lines = set(names), []
    for raw in query_paths.read_text().splitlines():
        fields = raw.split()
        if len(fields) >= 2 and not fields[0].startswith('#') and fields[0] in wanted:
            lines.append('\t'.join([fields[0], *(str(base / value) for value in fields[1:3])]))
            wanted.discard(fields[0])
    if wanted:
        raise ValueError(f'{len(wanted)} selected sample(s) not in {query_paths}: {sorted(wanted)[:5]}')
    return lines


def match_texts(args, vcfs):
    """Per grVCF, the text the patterns search: its path plus its samples'
    query-path lines (KSA, HPRC2, ... are often only in assembly paths)."""
    query_paths = Path(args.query_paths) if args.query_paths else (
        Path(args.cohort) / 'inputs' / 'query_paths.normalized.txt')
    lines = {}
    if query_paths.is_file():
        for raw in query_paths.read_text().splitlines():
            fields = raw.split()
            if fields and not fields[0].startswith('#'):
                lines[fields[0]] = raw
        print(f'[timecost] patterns search grVCF paths and {query_paths}', file=sys.stderr, flush=True)
    else:
        print(f'[timecost] no {query_paths}: patterns search grVCF paths only', file=sys.stderr, flush=True)
    return {path: '\n'.join([path, *(lines.get(sample, sample) for sample in vcf_samples(path))])
            for path in vcfs}


def parse_experiment(text):
    try:
        name, pattern, count = text.rsplit(':', 2)
        return name, pattern, int(count)
    except ValueError:
        raise argparse.ArgumentTypeError(f'expected NAME:REGEX:HAPLOTYPES, got {text!r}')


def exit_code(directory):
    """The recorded merge exit code; None when the run never finished (it
    may still be running)."""
    try:
        return json.loads((directory / 'timecost.json').read_text()).get('exit_code')
    except (OSError, ValueError):
        return None


def select(args, name, pattern, expected, vcfs, texts):
    """The experiment's grVCFs (and assemblies); None if it finished before."""
    directory = Path(args.output).resolve() / name
    if directory.exists() and any(directory.iterdir()) and not args.force:
        code = exit_code(directory)
        if code == 0:
            print(f'[timecost] {name}: finished before, skipped (--force reruns it)', file=sys.stderr, flush=True)
            return None
        if code is None:
            print(f'[timecost] {name}: earlier run never finished (interrupted), rerunning from scratch; '
                  'make sure no earlier driver or merge job is still running', file=sys.stderr, flush=True)
        else:
            print(f'[timecost] {name}: failed before (exit {code}), rerunning from scratch',
                  file=sys.stderr, flush=True)
    regex = re.compile(pattern)
    selected = [path for path in vcfs if regex.search(texts[path])]
    print(f'[timecost] {name}: pattern {pattern!r} selects {len(selected)} grVCF(s) '
          f'(expected {expected})', file=sys.stderr, flush=True)
    if len(selected) != expected and not args.allow_count_mismatch:
        raise ValueError(f'{name}: {len(selected)} grVCFs match {pattern!r}, expected {expected} '
                         '(check the pattern, or pass --allow-count-mismatch)')
    if not selected:
        raise ValueError(f'{name}: no grVCF matches {pattern!r}')
    lines = []
    if args.mode == 'exact':
        names = [sample for path in selected for sample in vcf_samples(path)]
        lines = subset_query_paths(args.cohort, args.query_paths, names)
    return dict(directory=directory, name=name, pattern=pattern, selected=selected, lines=lines)


def write_experiment(args, item, references):
    """A fresh experiment folder with its lists and merge command."""
    directory = item['directory']
    if directory.exists():
        shutil.rmtree(directory)
    directory.mkdir(parents=True)
    (directory / 'mergevcf.list').write_text(''.join(path + '\n' for path in item['selected']))
    command = [sys.executable, str(MERGE_TOOL), '-I', str(directory / 'mergevcf.list'),
               '-O', str(directory / 'merge'), '-t', str(args.threads), '--slurm',
               '--slurm-jobs', str(args.slurm_jobs)]
    if args.slurm_account:
        command += ['--slurm-account', args.slurm_account]
    if args.slurm_partition:
        command += ['--slurm-partition', args.slurm_partition]
    if args.mode == 'exact':
        (directory / 'query_paths.txt').write_text(''.join(line + '\n' for line in item['lines']))
        command += ['--exact', str(directory / 'query_paths.txt')]
        for reference in references:
            command += ['--reference-fasta', reference]
    command += shlex.split(args.merge_args)
    config = dict(name=item['name'], pattern=item['pattern'], haplotypes=len(item['selected']),
                  mode=args.mode, threads=args.threads, command=command)
    (directory / 'experiment.json').write_text(json.dumps(config, indent=1) + '\n')
    return directory


def exact_references(args):
    if args.mode != 'exact':
        return []
    if args.reference_fasta:
        return [str(Path(path).resolve()) for path in args.reference_fasta]
    sys.path.insert(0, str(SCRIPTS))
    import cohort_vcf_merge
    # Same defaults as the cohort merge: run reference + alternative loci.
    # Absolute, as the cohort folder may be given relative to here.
    return [str(Path(path).resolve()) for path in cohort_vcf_merge.exact_inputs(args.cohort, 'auto')[1]]


def start(args):
    chosen = list(args.experiment or EXPERIMENTS)
    if args.only:
        unknown = set(args.only) - {name for name, _, _ in chosen}
        if unknown:
            raise ValueError(f'--only: unknown experiment(s) {sorted(unknown)}')
        chosen = [item for item in chosen if item[0] in args.only]
    vcfs, source = cohort_vcfs(args.cohort, args.vcf_list)
    print(f'[timecost] {len(vcfs)} grVCF(s) from {source}', file=sys.stderr, flush=True)
    references = exact_references(args)
    texts = match_texts(args, vcfs)
    # Select everything first, so a wrong pattern stops before any merge starts.
    items = [item for item in (select(args, *experiment, vcfs, texts) for experiment in chosen) if item]
    failed = []
    for item in items:
        if run_experiment(write_experiment(args, item, references)):
            failed.append(item['name'])
    if failed:
        print(f'[timecost] failed: {", ".join(failed)}; repeating the command reruns them', file=sys.stderr)
    return 1 if failed else 0


# ---------------------------------------------------------------- timing

def run_experiment(directory):
    """Run the merge, timestamp its log and time its phases; its exit code."""
    config = json.loads((directory / 'experiment.json').read_text())
    started_at = time.time()
    started = time.monotonic()
    print(f'[timecost] {config["name"]}: {config["haplotypes"]} haplotypes started '
          f'{time.strftime("%Y-%m-%d %H:%M:%S")}', file=sys.stderr, flush=True)
    child = subprocess.Popen(config['command'], stdout=subprocess.PIPE, stderr=subprocess.STDOUT)
    # The once stages wait for their job, so a phase ends where the next starts.
    phases = [['plan', 0.0, None]]
    with open(directory / 'merge.log', 'w') as log:
        log.write('# ' + shlex.join(config['command']) + '\n')
        for raw in child.stdout:
            line = raw.decode(errors='replace')
            now = time.monotonic() - started
            match = STAGE_LINE.match(line)
            if match:
                phase = PHASE.get(match.group(1), match.group(1))
                if phase != phases[-1][0]:
                    phases[-1][2] = now
                    phases.append([phase, now, None])
            log.write(f'[{now:10.1f}s] {line}')
            log.flush()
    child.wait()
    wall = time.monotonic() - started
    phases[-1][2] = wall
    totals = {}
    with open(directory / 'phases.tsv', 'w') as handle:
        handle.write('phase\tstart_s\tend_s\twall_s\n')
        for name, begin, end in phases:
            handle.write(f'{name}\t{begin:.1f}\t{end:.1f}\t{end - begin:.1f}\n')
            totals[name] = round(totals.get(name, 0) + end - begin, 1)
    result = dict(
        experiment=config['name'], pattern=config['pattern'], haplotypes=config['haplotypes'],
        mode=config['mode'], threads=config['threads'], exit_code=child.returncode,
        wall_s=round(wall, 1), phases=totals,
        started=time.strftime('%Y-%m-%d %H:%M:%S', time.localtime(started_at)),
        finished=time.strftime('%Y-%m-%d %H:%M:%S'),
    )
    result['slurm'] = slurm_usage(directory)
    (directory / 'timecost.json').write_text(json.dumps(result, indent=1) + '\n')
    print(f'[timecost] {config["name"]}: exit {child.returncode}, wall {wall / 3600:.2f} h',
          file=sys.stderr, flush=True)
    if child.returncode:
        report_failure(directory)
    summary(directory.parent)
    return child.returncode


def sacct_rows(job_ids):
    """sacct allocation rows (no steps) of the jobs, JobID order kept."""
    fields = ['JobID', 'JobName', 'State', 'ElapsedRaw', 'AllocCPUS', 'Submit', 'Start', 'End']
    rows = []
    for begin in range(0, len(job_ids), 200):
        output = subprocess.run(
            ['sacct', '-X', '-n', '-P', '-j', ','.join(job_ids[begin:begin + 200]),
             '--format=' + ','.join(fields)],
            check=True, capture_output=True, text=True).stdout
        rows += [dict(zip(fields, line.split('|'))) for line in output.splitlines() if line.strip()]
    return rows


def job_phase(name):
    """merge-sv-chrom-3 -> chrom, merge-sv-scan -> sv-scan."""
    stage = re.sub(r'-\d+$', '', name[len('merge-'):] if name.startswith('merge-') else name)
    return PHASE.get(stage, stage)


def seconds_between(first, second):
    try:
        return max(0.0, (datetime.fromisoformat(second) - datetime.fromisoformat(first)).total_seconds())
    except ValueError:
        return 0.0


def slurm_usage(directory, attempts=6):
    """Added-up run, core and queue seconds of the merge's Slurm jobs (job
    IDs from merge.log), overall and per phase; also writes jobs.tsv. None
    when sacct is unavailable."""
    job_ids = []
    for line in (directory / 'merge.log').read_text(errors='replace').splitlines():
        match = JOB_ID_LINE.match(line)
        if match and match.group(1) not in job_ids:
            job_ids.append(match.group(1))
    if not job_ids:
        return None
    for attempt in range(attempts):
        try:
            rows = sacct_rows(job_ids)
        except (OSError, subprocess.CalledProcessError) as error:
            print(f'[timecost] sacct failed ({error}); `summary --refresh` retries later',
                  file=sys.stderr, flush=True)
            return None
        # Accounting can lag a few seconds behind the job's end.
        if len(rows) >= len(job_ids) or attempt == attempts - 1:
            break
        time.sleep(10)
    total = dict(jobs=0, job_s=0.0, core_s=0.0, queue_s=0.0)
    phases = {}
    with open(directory / 'jobs.tsv', 'w') as handle:
        handle.write('job_id	job_name	phase	state	cpus	run_s	core_s	queue_s	submit	start	end\n')
        for row in rows:
            run = float(row['ElapsedRaw'] or 0)
            cpus = int(row['AllocCPUS'] or 0)
            queue = seconds_between(row['Submit'], row['Start'])
            phase = job_phase(row['JobName'])
            handle.write('\t'.join(map(str, [row['JobID'], row['JobName'], phase, row['State'], cpus,
                                              f'{run:.0f}', f'{run * cpus:.0f}', f'{queue:.0f}',
                                              row['Submit'], row['Start'], row['End']])) + '\n')
            for bucket in (total, phases.setdefault(phase, dict(jobs=0, job_s=0.0, core_s=0.0, queue_s=0.0))):
                bucket['jobs'] += 1
                bucket['job_s'] += run
                bucket['core_s'] += run * cpus
                bucket['queue_s'] += queue
    missing = len(job_ids) - len(rows)
    if missing:
        print(f'[timecost] sacct has no record of {missing} of {len(job_ids)} job(s); '
              '`summary --refresh` retries', file=sys.stderr, flush=True)
    return dict(total, missing_jobs=missing, phases=phases)


def refresh(output):
    """Query sacct again for every recorded experiment."""
    for path in sorted(Path(output).glob('*/timecost.json')):
        result = json.loads(path.read_text())
        result['slurm'] = slurm_usage(path.parent)
        path.write_text(json.dumps(result, indent=1) + '\n')


def report_failure(directory, lines=15):
    """End of the newest non-empty Slurm error log (else of merge.log)."""
    logs = sorted((path for path in (directory / 'merge' / 'slurm_logs').glob('*.err')
                   if path.stat().st_size), key=lambda path: path.stat().st_mtime)
    source = logs[-1] if logs else directory / 'merge.log'
    try:
        tail = source.read_text(errors='replace').splitlines()[-lines:]
    except OSError:
        return
    print(f'[timecost] end of {source}:', file=sys.stderr)
    for line in tail:
        print(f'    {line}', file=sys.stderr)
    sys.stderr.flush()


def summary(output):
    output = Path(output)
    results = []
    for path in sorted(output.glob('*/timecost.json')):
        try:
            results.append(json.loads(path.read_text()))
        except (OSError, ValueError):
            continue
    results.sort(key=lambda item: item['haplotypes'])
    names = []
    for item in results:
        names += [name for name in item['phases'] if name not in names]
    hours = lambda seconds: '' if seconds is None else f'{seconds / 3600:.3f}'
    columns = ['experiment', 'haplotypes', 'mode', 'threads', 'exit_code', 'wall_h',
               'job_h', 'core_h', 'queue_h', 'jobs', 'started', 'finished',
               *(f'{name}_wall_h' for name in names), *(f'{name}_job_h' for name in names),
               *(f'{name}_core_h' for name in names)]
    lines = ['\t'.join(columns)]
    for item in results:
        slurm = item.get('slurm') or {}
        row = dict(item, wall_h=hours(item['wall_s']), job_h=hours(slurm.get('job_s')),
                   core_h=hours(slurm.get('core_s')), queue_h=hours(slurm.get('queue_s')),
                   jobs=slurm.get('jobs', ''))
        for name in names:
            row[f'{name}_wall_h'] = hours(item['phases'].get(name))
            phase = (slurm.get('phases') or {}).get(name, {})
            row[f'{name}_job_h'] = hours(phase.get('job_s'))
            row[f'{name}_core_h'] = hours(phase.get('core_s'))
        lines.append('\t'.join(str(row.get(column, '')) for column in columns))
    (output / 'summary.tsv').write_text('\n'.join(lines) + '\n')
    return output / 'summary.tsv'


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    commands = parser.add_subparsers(dest='command', required=True)
    begin = commands.add_parser('start', help='select the grVCFs of each experiment and merge them one by one')
    begin.add_argument('-C', '--cohort', required=True, help='cohort output folder (cohort_calls)')
    begin.add_argument('-O', '--output', required=True, help='experiments folder (timecost)')
    begin.add_argument('-t', '--threads', type=int, default=32,
                       help='merge_grvcfs.py -t, same for every experiment (default 32)')
    begin.add_argument('--mode', choices=('exact', 'cigar'), default='exact',
                       help='exact (default, as the cohort pipeline merge) or cigar (merge_grvcfs default)')
    begin.add_argument('-I', '--vcf-list', help='default: COHORT/inputs/mergevcf.list')
    begin.add_argument('--query-paths', help='--mode exact: default COHORT/inputs/query_paths.normalized.txt')
    begin.add_argument('--reference-fasta', action='append', default=[],
                       help='--mode exact: default the cohort run reference + checkpoints/alternative_loci.fa')
    begin.add_argument('--experiment', action='append', type=parse_experiment, metavar='NAME:REGEX:HAPLOTYPES',
                       help='replaces the default experiments (repeatable)')
    begin.add_argument('--only', nargs='+', metavar='NAME', help='run only these experiments')
    begin.add_argument('--slurm-account', default='')
    begin.add_argument('--slurm-partition', default='')
    begin.add_argument('--slurm-jobs', type=int, default=20, help='merge_grvcfs.py --slurm-jobs (default 20)')
    begin.add_argument('--merge-args', default='', help='quoted extra merge_grvcfs.py options')
    begin.add_argument('--allow-count-mismatch', action='store_true')
    begin.add_argument('--force', action='store_true', help='also rerun experiments that finished')
    collect = commands.add_parser('summary', help='rewrite OUTPUT/summary.tsv')
    collect.add_argument('-O', '--output', required=True)
    collect.add_argument('--refresh', action='store_true', help='query sacct again for every experiment')
    args = parser.parse_args(argv)
    if args.command == 'start':
        return start(args)
    if args.refresh:
        refresh(args.output)
    print(summary(args.output).read_text(), end='')
    return 0


if __name__ == '__main__':
    try:
        raise SystemExit(main())
    except (ValueError, OSError) as error:
        print(f'[timecost] ERROR: {error}', file=sys.stderr)
        raise SystemExit(1)
