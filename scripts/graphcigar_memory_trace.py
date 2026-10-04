"""Opt-in OOM diagnostics; GRAPH_CIGARTOREF_MEMORY_TRACE=<query substring>.

Use 1 to trace every candidate group. Stage entry is flushed before work so
SIGKILL leaves a useful last event. Sequence contents are never logged. Linux
PSS/private values distinguish inherited shared pages from worker allocations.
"""
from functools import wraps
import json
import os
import resource
import sys
import time


def memory_mib():
    peak = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss
    values = {'peak_rss_mib': peak / (1024 ** 2 if sys.platform == 'darwin' else 1024)}
    try:
        with open('/proc/self/smaps_rollup') as handle:
            fields = {parts[0].rstrip(':'): int(parts[1])
                      for line in handle if len(parts := line.split()) == 3
                      and parts[2] == 'kB'}
        values.update(rss_mib=fields['Rss'] / 1024,
                      pss_mib=fields['Pss'] / 1024,
                      private_mib=(fields.get('Private_Clean', 0)
                                   + fields.get('Private_Dirty', 0)) / 1024)
    except (OSError, ValueError, KeyError):
        pass
    return {key: round(value, 2) for key, value in values.items()}


def emit(event, stage, **details):
    record = dict(pid=os.getpid(), time=round(time.time(), 3),
                  event=event, stage=stage, **memory_mib(), **details)
    sys.stderr.write('[graphcigar memory] ' + json.dumps(record, sort_keys=True) + '\n')
    sys.stderr.flush()


def install(core, namespace, query_filter):
    active = False

    def wrap(fn, name):
        @wraps(fn)
        def traced(*args, **kwargs):
            if not active:
                return fn(*args, **kwargs)
            sizes = {f'arg{i}_len': len(value) for i, value in enumerate(args)
                     if isinstance(value, (str, bytes, list, tuple))}
            sizes.update({f'{key}_len': len(value) for key, value in kwargs.items()
                          if isinstance(value, (str, bytes, list, tuple))})
            emit('enter', name, **sizes)
            try:
                result = fn(*args, **kwargs)
            except BaseException:
                emit('error', name)
                raise
            emit('exit', name)
            return result
        return traced

    for name in (
        'build_provisional_row', 'build_stage1', 'weighted_lcs',
        '_segment_local_cuts', 'build_stage2', 'build_stage2_core_seeded_extensions',
        'build_stage3_polished', 'build_stage4_linear',
        '_realignment_hardness_allows', '_realign_tandem_stage2',
        '_minimap2_tandem_hits', '_mappy_payload_ops', 'minimap2_payload_ops',
        '_global_affine_payload_align', '_affine_global_fallback',
        '_tm_global_insert_align_dp', 'swspy',
        '_realign_masked_repeat_windows', '_tm_realign_secondary_sv_window',
        '_tm_consolidate_sv_ops', '_tm_realign_score_linked_indel_windows',
        'build_query_coverages', 'encode_large_insertions_in_piece',
        '_canonicalize_stage3_against_sequences', 'serialize_stage3_reflalign',
    ):
        fn = getattr(core, name, None)
        if callable(fn):
            setattr(core, name, wrap(fn, name))
    if hasattr(core, 'MaskedSVScorer'):
        cls = core.MaskedSVScorer
        cls.__init__ = wrap(cls.__init__, 'MaskedSVScorer.__init__')
    for name in ('_build_comparison_rows', 'format_merged_comparison_row',
                 'format_provisional_row', '_align_reference_flank_ops'):
        if callable(namespace.get(name)):
            namespace[name] = wrap(namespace[name], name)

    evaluate = namespace['_evaluate_column10_candidate_rows']

    @wraps(evaluate)
    def trace_group(candidate_group, *args, **kwargs):
        nonlocal active
        names = [base.query_name or base.label or '' for base, _members in candidate_group]
        if query_filter != '1' and not any(query_filter in name for name in names):
            return evaluate(candidate_group, *args, **kwargs)
        previous, active = active, True
        details = [dict(query=base.query_name, reference=base.ref_name,
                        query_coord=base.query_coord_text, reference_coord=base.ref_coord_text)
                   for base, _members in candidate_group]
        try:
            emit('enter', 'candidate_group', candidates=details)
            result = evaluate(candidate_group, *args, **kwargs)
            emit('exit', 'candidate_group')
            return result
        finally:
            active = previous

    namespace['_evaluate_column10_candidate_rows'] = trace_group
