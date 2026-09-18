"""Recall individual VCFs from saved cohort alignments, then merge them.

This route does not construct a Snakemake DAG or update upstream stage configs.
All required alignment inputs must already exist before any VCF is replaced.
"""
from __future__ import annotations

import json
from pathlib import Path
import re
import shlex
import subprocess
import sys
import tempfile

import cohort_vcf_merge as merge


ROOT = Path(__file__).resolve().parent


def _read_json(path):
    try:
        with Path(path).open() as handle:
            return json.load(handle)
    except FileNotFoundError as error:
        raise FileNotFoundError(f'--recall-only requires saved calling configuration: {path}') from error


def _stage_config(run, key):
    return _read_json(run[key]) if run.get(key) else run


def _stamps(paths):
    result = []
    for path in paths:
        path = Path(path)
        info = path.stat()
        result.append([str(path), info.st_size, info.st_mtime_ns])
    return result


def _write_json(path, value):
    with merge.snp.atomic_output(path) as handle:
        json.dump(value, handle, indent=2)
        handle.write('\n')


def _selected(run, assembly_list=None, samples=None):
    available = run['selected_samples']
    assemblies = {row['name']: row for row in run['assemblies']}
    requested = list(samples) if samples is not None else list(available)
    if assembly_list:
        listing = Path(assembly_list).expanduser().resolve()
        requested = []
        for number, raw in enumerate(listing.read_text().splitlines(), 1):
            if not raw.strip() or raw.lstrip().startswith('#'):
                continue
            fields = shlex.split(raw, comments=True)
            if len(fields) != 2:
                raise ValueError(f'{listing}:{number}: expected NAME FASTA')
            name, filename = fields
            path = Path(filename).expanduser()
            path = path if path.is_absolute() else listing.parent / path
            if name not in assemblies or path.resolve() != Path(assemblies[name]['fasta']).resolve():
                raise ValueError(f'{name}: --recall-only requires the same assembly path used for the saved alignments')
            requested.append(name)
    if not requested or len(requested) != len(set(requested)):
        raise ValueError('--recall-only requires a nonempty selection without duplicate samples')
    for name in [run['reference'], *requested]:
        if not re.fullmatch(r'[A-Za-z0-9.-]+_h[1-9][0-9]*', name):
            raise ValueError(f'invalid saved sample name: {name!r}')
        if name not in available or name not in assemblies:
            raise ValueError(f'{name}: no saved individual calling inputs; --recall-only cannot call new samples')
    return [run['reference'], *[name for name in requested if name != run['reference']]], assemblies


def _run_caller(command, destination, log_path):
    """Keep the previous VCF intact until its replacement succeeds."""
    log_path.parent.mkdir(parents=True, exist_ok=True)
    with tempfile.TemporaryDirectory(prefix='.recall-', dir=destination.parent) as temporary:
        temporary_vcf = Path(temporary) / destination.name
        command = list(command)
        command[command.index('--output') + 1] = str(temporary_vcf)
        with log_path.open('w') as log:
            log.write(shlex.join(command) + '\n')
            log.flush()
            with subprocess.Popen(command, stdout=subprocess.PIPE, stderr=subprocess.STDOUT, text=True) as process:
                try:
                    for line in process.stdout:
                        print(line, end='', flush=True)
                        log.write(line)
                        log.flush()
                    code = process.wait()
                except BaseException:
                    process.terminate()
                    process.wait()
                    raise
            if code:
                raise RuntimeError(f'VCF recall failed for {destination.parent.name}; previous VCF preserved. See {log_path}')
        if not temporary_vcf.is_file() or not temporary_vcf.stat().st_size:
            raise RuntimeError(f'VCF recall produced no VCF: {destination}')
        temporary_vcf.replace(destination)
        report = temporary_vcf.with_suffix('.inconsistencies.tsv')
        if report.is_file():
            report.replace(destination.with_suffix('.inconsistencies.tsv'))


def run(output, mode='svonly', *, graph=None, reference=None, assembly_list=None,
        samples=None, processes=1, cutoff=20, dry_run=False,
        caller_settings=None, merge_settings=None):
    output = Path(output).expanduser().resolve()
    saved = _read_json(output / 'inputs/cohort_call.run.json')
    if saved.get('protocol') != 'cohort-call-v1':
        raise ValueError('--recall-only requires a cohort-call-v1 run configuration')
    if Path(saved['output_folder']).resolve() != output:
        raise ValueError('--recall-only requires the original calling output directory')
    if graph and Path(graph).expanduser().resolve() != Path(saved['graph_folder']).resolve():
        raise ValueError('--recall-only graph differs from the graph used for the saved alignments')
    if reference and reference != saved['reference'] and Path(reference).expanduser().resolve() != Path(saved['reference_fasta']).resolve():
        raise ValueError(f"--recall-only cannot change reference: saved alignments use {saved['reference']}")
    if processes < 1 or cutoff < 1 or mode not in ('all', 'svonly', 'snp', 'svindel'):
        raise ValueError('invalid recall mode, process count, or cutoff')
    selected, assemblies = _selected(saved, assembly_list, samples)
    vcf = _stage_config(saved, 'vcf_config')
    merge_config = _stage_config(saved, 'merge_config')
    caller_settings = dict(caller_settings or {})
    if 'max_extension' in caller_settings:
        caller_settings['max_impute'] = caller_settings.pop('max_extension')
    format_processes = int(caller_settings.pop('format_processes', saved.get('format_processes', 16)))
    if format_processes < 1:
        raise ValueError('--format-processes must be positive')
    vcf.update(caller_settings)
    merge_settings = dict(merge_settings or {})
    merge_processes = merge_settings.pop('processes', processes)
    kmermatch_override = merge_settings.pop('kmermatch', None)
    merge_config.update(merge_settings)
    if vcf.get('seqcompress', saved.get('seqcompress', False)):
        raise ValueError('--recall-only merging requires plain INFO/SEQ; saved seqcompress is enabled')

    # Reuse saved calling settings unless explicitly overridden. No template,
    # alignment, gap-fill, or layout config is written.
    options = [
        '--processes', str(processes), '--format-processes', str(min(processes, format_processes)),
        '--edge-blackregion', str(vcf.get('edge_blackregion', saved.get('edge_blackregion', 0))),
        '--max-impute', str(vcf.get('max_impute', saved.get('max_extension', 10000))),
        '--realignment' if vcf.get('realignment', saved.get('realignment', True)) else '--no-realignment',
        '--resolve-conflicts-by-alignment-score' if vcf.get('resolve_conflicts_by_alignment_score', True) else '--no-resolve-conflicts-by-alignment-score',
        '--PAtag' if vcf.get('pa_tag', saved.get('pa_tag', True)) else '--noPAtag',
    ]
    if vcf.get('separate_adjacent_indels', False):
        options.append('--separate-adjacent-indels')
    if mode == 'svonly':
        options.extend(['--svcutoff', str(3 * cutoff // 5)])
    template = output / 'checkpoints/local_reference_templates.fa'
    shared = [Path(saved['reference_fasta']), Path(saved['reference_fai']),
              template, Path(str(template) + '.fai')]
    scripts = [ROOT / name for name in ('graphreftovcf_persample.py', 'graphreftovcf.py',
               'graph_cigar_payloads.py', 'alternative_intervals.py', 'local_reference_templates.py')]
    jobs = []
    missing = set()
    for name in selected:
        sample_root = output / 'samples' / name
        graphcigar = sample_root / f'{name}.graphcigartoreffix.tsv'
        lift = sample_root / f'{name}.genomelift.tsv'
        liftfix = sample_root / f'{name}.genomeliftfix.tsv'
        assembly = assemblies[name]
        inputs = [*shared, *scripts, graphcigar, lift, liftfix,
                  Path(assembly['fasta']), Path(assembly['fai'])]
        missing.update(str(path) for path in inputs if not path.is_file())
        destination = sample_root / f'{name}.vcf'
        command = [sys.executable, str(scripts[0]), '--input', str(graphcigar),
                   '--ref', saved['reference_fasta'], '--local-reference-templates', str(template),
                   '--fasta-query', assembly['fasta'], '--coord-map', f'{lift},{liftfix}',
                   '--output', str(destination), '--columns', name, *options]
        jobs.append((name, inputs, command, destination))
    if missing:
        raise FileNotFoundError('--recall-only requires existing inputs; it will not rebuild them:\n' + '\n'.join(sorted(missing)))
    settings = {key: merge_config.get(key, default) for key, default in (
        ('merge_distance', 500), ('size_similarity', .7), ('sequence_similarity', .7), ('var_in_insert', 100))}
    kmermatch = kmermatch_override or saved.get('kmermatch', merge.sv.DEFAULT_KMERMATCH)
    if mode != 'snp' and not dry_run:
        from graphvcfmerge_kmer import resolve_kmermatch
        kmermatch = resolve_kmermatch(kmermatch)
    print(f'[cohort-recall] {len(jobs)} samples; reference={saved["reference"]}; mode={mode}; {processes} workers', flush=True)
    for number, (name, inputs, command, destination) in enumerate(jobs, 1):
        print(f'[cohort-recall] {number}/{len(jobs)}: {name}', flush=True)
        semantic_command = list(command)
        for capacity in ('--processes', '--format-processes'):
            index = semantic_command.index(capacity)
            del semantic_command[index:index + 2]
        signature = {'protocol': 'cohort-vcf-recall-v1', 'command': semantic_command, 'inputs': _stamps(inputs)}
        marker = output / 'checkpoints/vcf_recall' / f'{name}.json'
        if marker.is_file() and destination.is_file():
            previous = _read_json(marker)
            if previous == dict(signature, outputs=_stamps([destination])):
                print(f'[cohort-recall] reuse completed VCF: {destination}', flush=True)
                continue
        print('[cohort-recall] ' + shlex.join(command), flush=True)
        if not dry_run:
            log_path = output / 'logs/vcf_recall' / f'{name}.log'
            _run_caller(command, destination, log_path)
            _write_json(marker, dict(signature, outputs=_stamps([destination])))
    return merge.run(output, mode, paths=[job[3] for job in jobs], processes=merge_processes,
                     cutoff=cutoff, kmermatch=kmermatch, dry_run=dry_run,
                     keep_merge_tmpdir=merge_config.get('keep_merge_tmpdir', saved.get('keep_merge_tmpdir', False)),
                     **settings)
