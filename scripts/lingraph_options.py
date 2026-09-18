"""Expose backend tuning options through LinGraph without duplicating defaults."""
import argparse
import copy
from functools import lru_cache


# LinGraph owns run selection, paths, stage control, and public mode names.
# Backend-only maintenance commands are intentionally not copied into a run.
COMMON = {'help', 'dry_run', 'unlock', 'snakemake', 'snakemake_args', 'slurm',
          'slurm_jobs', 'slurm_account', 'slurm_time', 'slurm_partition', 'slurm_args'}
MANAGED = {
    'build': COMMON | {'reference', 'assemblies', 'bed', 'alternative', 'find_novel_loci',
        'bed_grouped', 'graph_folder', 'output_folder', 'continue_cohort_call',
        'cohort_call_args', 'cores', 'alignment_mode', 'compressgraph', 'static_block'},
    'call': COMMON | {'reference', 'graph_folder', 'output_folder', 'assemblies',
        'partition_list', 'cores', 'sv_only_size', 'minsvsize', 'merge_mode',
        'mc_graph', 'merge_only', 'recall_only'},
    'sample': COMMON | {'sample', 'fasta_query', 'hotspot_query', 'reference_sample',
        'fasta_reference', 'reference_align', 'reference_blocks', 'graph_folder',
        'graph_list', 'output_root', 'threads', 'fast', 'svcutoff',
        'template_assemblies', 'kmer_searcher', 'resume', 'rerun', 'windowmasker', 'samtools'},
    'gfa': COMMON | {'vcf', 'variant_mode', 'svcutoff', 'insertion_only',
        'query_fasta_list', 'graph_folder', 'reference_fasta', 'reference_haplotype',
        'reference_fai', 'local_reference_templates', 'alternatives_fasta',
        'local_path_fasta', 'output', 'processes', 'bed_only', 'validate_only'},
}


@lru_cache(None)
def backend_parser(stage):
    if stage == 'build':
        from graph_build_snakemake.run_graph_pipeline import build_parser
    elif stage in {'call', 'merge'}:
        from cohort_call_snakemake.run_cohort_call_pipeline import build_parser
    elif stage == 'sample':
        from run_sample_pipeline import build_parser
    elif stage == 'gfa':
        from merged_vcf_to_gfa import build_parser
    else:
        raise ValueError(stage)
    return build_parser()


def frontend_dest(stage, action):
    # Graph construction's alternatives are FASTA; calling's override is BED.
    return 'alternative_bed' if stage == 'call' and action.dest == 'alternative' else action.dest


def tuning_actions(stage):
    if stage == 'merge':
        return [action for action in tuning_actions('call') if action.dest in {
            'merge_distance', 'size_similarity', 'sequence_similarity', 'merge_processes',
            'var_in_insert', 'kmermatch', 'keep_merge_tmpdir'}]
    return [action for action in backend_parser(stage)._actions
            if action.option_strings and action.dest not in MANAGED[stage]]


def add_options(parser, stage, title, show_advanced):
    group = parser.add_argument_group(title)
    mutexes = {}
    for action in tuning_actions(stage):
        names = (['--alternative-bed'] if stage == 'call' and action.dest == 'alternative'
                 else [name for name in action.option_strings if name.startswith('--')])
        if any(name in parser._option_string_actions for name in names):
            continue
        cloned = copy.copy(action)
        cloned.option_strings = names
        cloned.dest = frontend_dest(stage, action)
        # Only explicitly supplied tuning values are forwarded. Backends keep
        # their defaults, and a value shared by two stages reaches both.
        cloned.default = argparse.SUPPRESS
        cloned.required = False
        if not show_advanced:
            cloned.help = argparse.SUPPRESS
        container = group
        for source in backend_parser(stage)._mutually_exclusive_groups:
            if action in source._group_actions:
                if id(source) not in mutexes:
                    mutexes[id(source)] = group.add_mutually_exclusive_group()
                container = mutexes[id(source)]
                break
        container._add_action(cloned)


def values(args, stage):
    """Return explicit values keyed by their backend destination name."""
    return {action.dest: getattr(args, frontend_dest(stage, action))
            for action in tuning_actions(stage)
            if hasattr(args, frontend_dest(stage, action))}


def forward(args, stage):
    result = []
    emitted = set()
    for action in tuning_actions(stage):
        dest = frontend_dest(stage, action)
        if dest in emitted or not hasattr(args, dest):
            continue
        value = getattr(args, dest)
        flag = next(name for name in action.option_strings if name.startswith('--'))
        if stage == 'merge' and action.dest == 'merge_processes':
            flag = '--processes'
        if isinstance(action, argparse._StoreConstAction):
            if value != action.const:
                continue
            result.append(flag)
        elif isinstance(action, argparse._AppendAction):
            for item in value:
                result.extend([flag, *map(str, item if isinstance(item, (list, tuple)) else [item])])
        elif isinstance(value, (list, tuple)):
            result.extend([flag, *map(str, value)])
        else:
            result.append(flag + '=' + str(value))
        emitted.add(dest)
    return result


def merge_settings(args):
    settings = values(args, 'merge')
    settings['processes'] = settings.pop('merge_processes', args.threads)
    if settings['processes'] < 1:
        raise ValueError('--merge-processes must be positive')
    for name in ('merge_distance', 'var_in_insert'):
        if settings.get(name, 0) < 0:
            raise ValueError('--' + name.replace('_', '-') + ' cannot be negative')
    for name in ('size_similarity', 'sequence_similarity'):
        if not 0 <= settings.get(name, .7) <= 1:
            raise ValueError('--' + name.replace('_', '-') + ' must be between 0 and 1')
    return settings
