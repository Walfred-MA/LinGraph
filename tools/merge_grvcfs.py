#!/usr/bin/env python3
"""Merge individual grVCFs into cohort grVCFs: LinGraph's merge step only.

Runs the pipeline's merge (scripts/cohort_vcf_merge.py) on the given files
and writes OUTPUT/cohort.sv.vcf, cohort.indel.vcf and cohort.snp.vcf.

  CIGAR version (default): each sample keeps its own breakpoint, size and
      alignment to its row's representative.
  --exact QUERY_PATHS: members are realigned against their assemblies and
      moved onto the representative's breakpoint (needs the assemblies,
      NAME FASTA [FAI] per line, and --reference-fasta).

The merge runs as stages, in order: SV scan, SV merge per chromosome, SV
concat, SNP prepare, SNP merge per chromosome, SNP concat, publish (stages a
mode does not need are skipped). Each stage is its own process:
  default   locally, one after another (chromosomes one at a time), -t workers
  --slurm   one Slurm job per stage; run-once stages get -t CPUs and
            --slurm-memory (default 64G below 100 input grVCFs, 128G from 100
            on), chromosome stages run as up to --slurm-jobs jobs at once
            (SV: min(16, -t) CPUs and 2G per CPU, or 64G from 100 input
            grVCFs on, chr1 doubled; SNP: 32 CPUs, 64G)
Finished stages are recorded under OUTPUT/tmp/merge_only/; repeating the
command resumes after the last finished stage or chromosome.

Examples:
  python tools/merge_grvcfs.py -i a.vcf b.vcf c.vcf -O merged -t 16
  python tools/merge_grvcfs.py -I vcfs.list -O merged --exact query_paths.txt \\
      --reference-fasta chm13.fa -t 16
  python tools/merge_grvcfs.py -I vcfs.list -O merged -t 32 --slurm --slurm-jobs 20 \\
      --slurm-account ACCOUNT --slurm-partition PARTITION
"""
import argparse
from concurrent.futures import ThreadPoolExecutor
import hashlib
import json
from pathlib import Path
import re
import shlex
import subprocess
import sys

# Pipeline scripts: ROOT/scripts (repository layout) or ROOT (flat copy).
ROOT = Path(__file__).resolve().parents[1]
SCRIPTS = ROOT / 'scripts' if (ROOT / 'scripts' / 'cohort_vcf_merge.py').is_file() else ROOT
sys.path.insert(0, str(SCRIPTS))

import cohort_vcf_merge  # noqa: E402
import graphvcfmerge  # noqa: E402

SLURM_HELPER = SCRIPTS / 'graph_build_snakemake' / 'workflow' / 'scripts' / 'pipeline_inputs.py'
STAGES = ('sv-scan', 'sv-chrom', 'sv-concat', 'snp-prepare', 'snp-chrom', 'snp-concat', 'publish')
# Primary chr1 names (as graphvcfmerge_stages.is_chromosome_one): its SV job gets
# twice the CPUs and memory.
CHR1 = re.compile(r'(?:chr1|1|NC_000001(?:\.\d+)?|NC_060925(?:\.\d+)?)', re.IGNORECASE)


def input_vcfs(args):
    paths = list(args.input or [])
    if args.vcf_list:
        paths.extend(graphvcfmerge.read_vcf_input_lists([args.vcf_list]))
    paths = graphvcfmerge.expand_vcf_inputs(paths)
    if not paths:
        raise ValueError('no input grVCFs: give -i FILE ... and/or -I LIST')
    missing = [path for path in paths if not Path(path).is_file()]
    if missing:
        raise FileNotFoundError(f'missing input grVCF(s): {missing[:5]}')
    return paths


def parse_args(argv=None):
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('-i', '--input', nargs='+', metavar='VCF', help='individual grVCFs (globs allowed)')
    parser.add_argument('-I', '--vcf-list', metavar='FILE', help='one grVCF path per line (relative to the list)')
    parser.add_argument('-O', '--output', metavar='DIR', help='output folder')
    parser.add_argument('-t', '--threads', type=int, default=1,
                        help='worker processes; with --slurm also the CPUs of run-once jobs (default: 1)')
    parser.add_argument('--exact', metavar='QUERY_PATHS',
                        help='realign against the assemblies listed here (exact version)')
    parser.add_argument('--reference-fasta', action='append', default=[], metavar='FASTA',
                        help='--exact: indexed reference FASTA; repeat for local reference templates')
    kinds = parser.add_mutually_exclusive_group()
    kinds.add_argument('--svonly', dest='mode', action='store_const', const='svonly',
                       help='only cohort.sv.vcf (SVs at or above --svcutoff)')
    kinds.add_argument('--svindel', dest='mode', action='store_const', const='svindel',
                       help='cohort.sv.vcf and cohort.indel.vcf, no SNPs')
    kinds.add_argument('--snp', dest='mode', action='store_const', const='snp',
                       help='only cohort.snp.vcf')
    parser.add_argument('--svcutoff', type=int, default=20, help='SV size cutoff (default: 20)')
    parser.add_argument('--merge-distance', type=int, default=500)
    parser.add_argument('--size-similarity', type=float, default=.7)
    parser.add_argument('--sequence-similarity', type=float, default=.7)
    parser.add_argument('--var-in-insert', type=int, default=100,
                        help='call variants inside merged insertions (0: off)')
    parser.add_argument('--keep-full-locus-dup-insertions', action='store_true',
                        help='retain full-locus duplication parent insertions')
    parser.add_argument('--keep-merge-tmpdir', action='store_true',
                        help='keep OUTPUT/tmp/merge_only after the merged files are written')
    slurm = parser.add_argument_group('Slurm (each stage as a job; default: local, one stage at a time)')
    slurm.add_argument('--slurm', action='store_true', help='submit the stages as Slurm jobs')
    slurm.add_argument('--slurm-jobs', type=int, default=20,
                       help='chromosome jobs running at once (default: 20)')
    slurm.add_argument('--slurm-account', default='')
    slurm.add_argument('--slurm-partition', default='')
    slurm.add_argument('--slurm-time', default='200:00:00', help='per job (default: 200:00:00)')
    slurm.add_argument('--slurm-memory', default='',
                       help='memory of the run-once jobs (default: 64G below 100 input grVCFs, else 128G)')
    slurm.add_argument('--slurm-args', default='', metavar='TEXT', help='quoted extra sbatch options for every job')
    # One stage of a plan written by the driver (how the driver runs each stage).
    parser.add_argument('--stage', choices=STAGES, help=argparse.SUPPRESS)
    parser.add_argument('--plan', help=argparse.SUPPRESS)
    parser.add_argument('--chrom', help=argparse.SUPPRESS)
    args = parser.parse_args(argv)
    if args.threads < 1:
        parser.error('--threads must be positive')
    if args.stage:
        if not args.plan:
            parser.error('--stage needs --plan')
        return args
    if not args.output:
        parser.error('-O/--output is required')
    if args.exact and not args.reference_fasta:
        parser.error('--exact needs --reference-fasta')
    if args.exact and args.mode not in (None, 'all'):
        parser.error('--exact writes all three files; drop --svonly/--svindel/--snp')
    if args.reference_fasta and not args.exact:
        parser.error('--reference-fasta is only used with --exact')
    if args.slurm_jobs < 1:
        parser.error('--slurm-jobs must be positive')
    if not args.slurm and (args.slurm_account or args.slurm_partition or args.slurm_memory or args.slurm_args):
        parser.error('Slurm settings require --slurm')
    try:
        args.slurm_args = shlex.split(args.slurm_args)
    except ValueError as error:
        parser.error(f'invalid --slurm-args: {error}')
    return args


def run_stage(args):
    """One stage of the plan, in this process."""
    plan = json.loads(Path(args.plan).read_text())
    threads = args.threads
    if args.stage in ('sv-chrom', 'snp-chrom') and not args.chrom:
        raise ValueError(f'--stage {args.stage} needs --chrom')
    {
        'sv-scan': lambda: cohort_vcf_merge.stage_sv_scan(plan, threads),
        'sv-chrom': lambda: cohort_vcf_merge.stage_sv_chrom(plan, [args.chrom], threads),
        'sv-concat': lambda: cohort_vcf_merge.stage_sv_concat(plan, threads),
        'snp-prepare': lambda: cohort_vcf_merge.stage_snp_prepare(plan, threads),
        'snp-chrom': lambda: cohort_vcf_merge.stage_snp_chrom(plan, args.chrom, threads),
        'snp-concat': lambda: cohort_vcf_merge.stage_snp_concat(plan),
        'publish': lambda: cohort_vcf_merge.stage_publish(plan),
    }[args.stage]()
    return 0


class Stages:
    """Runs the plan's stages locally or as Slurm jobs; skips recorded ones."""

    def __init__(self, args, plan, plan_path):
        self.args, self.plan, self.plan_path = args, plan, plan_path
        self.done = Path(plan['root']) / 'stages'
        self.done.mkdir(parents=True, exist_ok=True)
        self.output = Path(plan['output'])
        self.once_memory = args.slurm_memory or ('64G' if len(plan['paths']) < 100 else '128G')

    def _marker(self, stage, chrom=None):
        name = stage if chrom is None else f"{stage}.{hashlib.sha1(chrom.encode()).hexdigest()[:16]}"
        return self.done / f'{name}.done'

    def _command(self, stage, threads, chrom=None):
        return [sys.executable, str(Path(__file__).resolve()), '--stage', stage,
                '--plan', str(self.plan_path), '-t', str(threads),
                *(['--chrom', chrom] if chrom is not None else [])]

    def _slurm(self, command, cpus, memory, job):
        args = self.args
        wrapper = [sys.executable, str(SLURM_HELPER), 'run-slurm', '--cpus', str(cpus),
                   '--memory', memory, '--job-name', job,
                   '--log-dir', str(self.output / 'slurm_logs'), '--time', args.slurm_time]
        if args.slurm_account:
            wrapper += ['--account', args.slurm_account]
        if args.slurm_partition:
            wrapper += ['--partition', args.slurm_partition]
        wrapper += ['--sbatch-arg=' + value for value in args.slurm_args]
        return [*wrapper, '--', *command]

    def _run(self, command, what):
        print(f'[merge_grvcfs] {what}: {shlex.join(command)}', file=sys.stderr, flush=True)
        return subprocess.run(command, check=False).returncode

    def once(self, stage):
        marker = self._marker(stage)
        if marker.is_file():
            print(f'[merge_grvcfs] {stage}: done before', file=sys.stderr, flush=True)
            return
        command = self._command(stage, self.args.threads)
        if self.args.slurm:
            command = self._slurm(command, self.args.threads, self.once_memory, f'merge-{stage}')
        if self._run(command, stage):
            raise RuntimeError(f'stage {stage} failed; see the messages above'
                               + (f' and {self.output / "slurm_logs"}' if self.args.slurm else ''))
        marker.touch()

    def _resources(self, stage, chrom):
        if stage == 'snp-chrom':
            return 32, '64G'
        multiplier = 2 if CHR1.fullmatch(chrom) else 1
        cpus = min(16, self.args.threads) * multiplier
        # As the pipeline: 2G per CPU, or 64G per job from 100 input grVCFs on.
        return cpus, (f'{64 * multiplier}G' if len(self.plan['paths']) >= 100 else f'{2 * cpus}G')

    def each(self, stage, chroms, names=None):
        """Run a chromosome stage; ``names`` labels manifest keys in messages."""
        label = (names or {}).get
        pending = [chrom for chrom in chroms if not self._marker(stage, chrom).is_file()]
        if len(pending) < len(chroms):
            print(f'[merge_grvcfs] {stage}: {len(chroms) - len(pending)} of {len(chroms)} '
                  'chromosome(s) done before', file=sys.stderr, flush=True)
        if not self.args.slurm:
            for number, chrom in enumerate(pending, 1):
                what = f'{stage} {number}/{len(pending)} {label(chrom, chrom)}'
                if self._run(self._command(stage, self.args.threads, chrom), what):
                    raise RuntimeError(f'{stage} failed for {label(chrom, chrom)}; see the messages above')
                self._marker(stage, chrom).touch()
            return

        def submit(item):
            number, chrom = item
            cpus, memory = self._resources(stage, label(chrom, chrom))
            command = self._slurm(self._command(stage, cpus, chrom), cpus, memory,
                                  f'merge-{stage}-{number}')
            code = self._run(command, f'{stage} {number}/{len(pending)} {label(chrom, chrom)}')
            if not code:
                self._marker(stage, chrom).touch()
            return label(chrom, chrom), code

        if not pending:
            return
        with ThreadPoolExecutor(min(self.args.slurm_jobs, len(pending))) as pool:
            failed = [chrom for chrom, code in pool.map(submit, enumerate(pending, 1)) if code]
        if failed:
            raise RuntimeError(f'{stage} failed for {len(failed)} chromosome(s): {", ".join(failed[:8])}; '
                               f'see {self.output / "slurm_logs"}; repeat the command to rerun only these')


def main(argv=None):
    args = parse_args(argv)
    if args.stage:
        return run_stage(args)
    paths = input_vcfs(args)
    output = Path(args.output).resolve()
    output.mkdir(parents=True, exist_ok=True)
    mode = 'all' if args.exact else (args.mode or 'all')
    print(f'[merge_grvcfs] {len(paths)} grVCF(s) -> {output} '
          f'({"exact" if args.exact else "CIGAR"} version; '
          f'{"Slurm jobs" if args.slurm else "local, one stage at a time"})', file=sys.stderr, flush=True)
    plan = cohort_vcf_merge.merge_plan(
        output, mode, paths, cutoff=args.svcutoff, merge_distance=args.merge_distance,
        size_similarity=args.size_similarity, sequence_similarity=args.sequence_similarity,
        var_in_insert=args.var_in_insert,
        ignore_full_locus_dup_insertions=not args.keep_full_locus_dup_insertions,
        exact=str(Path(args.exact).resolve()) if args.exact else None,
        reference_fastas=[str(Path(path).resolve()) for path in args.reference_fasta])
    plan_path = Path(plan['root']) / 'plan.json'
    plan_path.write_text(json.dumps(plan, indent=1) + '\n')
    stages = Stages(args, plan, plan_path)
    if mode != 'snp':
        stages.once('sv-scan')
        stages.each('sv-chrom', cohort_vcf_merge.sv_chroms(plan))
        stages.once('sv-concat')
    if mode in ('all', 'snp'):
        stages.once('snp-prepare')
        manifest = Path(plan['root']) / 'snp_shards' / 'manifest.json'
        names = {key: str(value.get('chrom', key)) if isinstance(value, dict) else str(value)
                 for key, value in json.loads(manifest.read_text())['chroms'].items()}
        stages.each('snp-chrom', cohort_vcf_merge.snp_chroms(plan), names)
        stages.once('snp-concat')
    stages.once('publish')
    targets = cohort_vcf_merge.output_paths(output, mode)
    if not args.keep_merge_tmpdir:
        cohort_vcf_merge.cleanup_plan(plan, targets)
    for leftover in (output / 'tmp' / 'merge_only', output / 'tmp'):
        try:
            leftover.rmdir()          # only when the merge left it empty
        except OSError:
            pass
    for target in targets:
        print(target)
    return 0


if __name__ == '__main__':
    try:
        raise SystemExit(main())
    except (ValueError, OSError, RuntimeError) as error:
        print(f'[merge_grvcfs] ERROR: {error}', file=sys.stderr)
        raise SystemExit(1)
