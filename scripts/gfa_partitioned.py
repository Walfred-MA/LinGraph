"""Partitioned rGFA export: bounded memory for very large cohorts.

``--partition-records N`` processes the variants in partitions of about N VCF
records instead of all at once, trading extra passes for memory. The GFA and
its sidecars are identical to a one-pass run (only the unresolved report is
ordered by partition).

A partition is a set of whole variant trees: a variant, everything nested in
it or copying from it, and its parents. Trees only share FASTA roots
(reference, alternative, template, catalog and duplication paths), so:

1. scan: one light pass over the VCFs (no sample columns) records each data
   line's byte offset and its tree links, the target names that reachable
   variants need, and the full-locus duplication alignments.
2. global roots: templates, catalog paths and duplication clusters are
   resolved once, from those names, exactly as the one-pass run does.
3. partition pass: each partition's records are read by offset and run
   through the usual extraction, verification and topology. Its leaves,
   cut points, hits and sidecar parts are kept on disk.
4. numbering: root segments are numbered first, as in one pass. Event ranks
   are the merge of the partitions' dependency orders by VCF line (the
   one-pass heap order), and event segment IDs follow from them.
5. write pass: each partition writes its S and P blocks and its links, then
   the blocks are copied in global rank order, links are sorted globally and
   the sidecars merged.
"""
from array import array
from collections import OrderedDict
import gc
import heapq
import os
from pathlib import Path
import pickle
import shutil
import tempfile
import time

import numpy as np

import gfa_topology as topology
from gfa_interval_pipeline import (
    CHUNK, Event, Hit, PathSpec, QueryCandidates, Root, Run, _copy_parts,
    _drop_unverified, _extract_queries, _flatten, _index_roots, _link_leaves,
    _path_leaves, _path_and_link_leaves, _prepare_intervals, _reachable_events, _report_unresolved,
    _valid_name, _write_links, _write_metadata, _write_parts, _write_segments,
    _write_source_sidecars, _write_variant_index)
from gfa_interval_metadata import TEMPLATE_PATH
from gfa_query_anchors import read_query_list


class _Last(dict):
    """Records only the most recent assignment (rows() side channels)."""
    key = last = None

    def __setitem__(self, key, value):
        self.key, self.last = key, value

    def get_for(self, key, default=None):
        """The value last set for ``key``; ``default`` if it was set for another."""
        return self.last if self.key == key else default


def _digits(names):
    return {int(name) for name in names
            if name.isascii() and name.isdigit() and name == str(int(name))}


def _skip(ordinals, forbidden):
    """Map 1-based segment ordinals to IDs that skip numeric path names."""
    numbers = np.asarray(ordinals, np.int64).copy()
    for value in sorted(forbidden):
        numbers[numbers >= value] += 1
    return numbers


# ---------------------------------------------------------------- scan

def _scan(args, roots, log):
    from gfa_catalog_paths import read_alternative_loci
    from gfa_interval_metadata import RecordIndex, rows, vcf_paths
    index = RecordIndex()
    kinds, lines = _Last(), _Last()
    dup_alignments = {}
    yielded = array('q')
    retained = array('b')
    parents = {}            # non-SNP variant ID -> data line (possible parents)
    targets = {}            # data line -> non-root names it references
    hashes = array('q')
    # Breakpoint keys (parent name hash mixed with the start or end), so rows
    # that abut on one parent (gfa_junctions) are put in one tree.
    start_keys, end_keys = array('q'), array('q')
    loci = set(read_alternative_loci(list(vcf_paths(args.vcf))))
    collisions = set()
    started = time.monotonic()
    for record in rows(args.vcf, frozenset(), args.size_cutoff,
                       insertion_only=getattr(args, 'insertion_only', None) is not None,
                       event_kinds=kinds, record_indexes=lines,
                       dup_alignments=dup_alignments, index=index, template_sites=True):
        identifier, chrom, pos, ref_end, _length, values, _query, _reason, keep = record
        # Only SNP/DEL/SUB rows record a kind; an insertion never matches.
        kind = kinds.get_for(identifier, 'insertion')
        line = lines.get_for(identifier)
        if getattr(args, 'variant_mode', 'all') == 'svonly' and kind == 'snp':
            keep = False
        yielded.append(line)
        retained.append(bool(keep))
        hashes.append(hash(identifier))
        # Top-level breakpoints as StableResolver.breakpoints; nested rows
        # only ever join rows of their own parent, already their tree.
        start = pos - 1 if kind == 'snp' else pos
        deleted = sum(value[5] - value[4] for value in values if value[2] == 'D')
        end = (pos if kind == 'snp' else start + deleted if deleted else
               start if ref_end == pos else ref_end)
        parent = hash(chrom)
        start_keys.append(_mix(parent, start))
        end_keys.append(_mix(parent, end))
        if kind != 'snp':
            parents[identifier] = line
        names = {chrom, *(value[3] for value in values if value[3])}
        names.difference_update(roots)
        if names:
            targets[line] = tuple(names)
        if identifier.startswith(('alt_', 'dup_')) or identifier in loci:
            collisions.add(identifier)
        if len(yielded) % 1_000_000 == 0:
            log(f'Scan: {len(yielded)} variants')
    hashes = np.frombuffer(hashes, np.int64)
    ordered = np.sort(hashes)
    if len(ordered) > 1 and (ordered[1:] == ordered[:-1]).any():
        repeated = set(ordered[1:][ordered[1:] == ordered[:-1]].tolist())
        lines_at = np.frombuffer(yielded, np.int64)[np.isin(hashes, list(repeated))]
        raise ValueError('duplicate insertion ID or variant ID at VCF data lines '
                         + ', '.join(str(value + 1) for value in lines_at[:10].tolist()))
    log(f'Scan: {len(yielded)} variants, {len(targets)} with non-reference targets, '
        f'{len(dup_alignments)} duplication alignments in {time.monotonic() - started:.1f}s')
    return dict(index=index, yielded=np.frombuffer(yielded, np.int64),
                retained=np.frombuffer(retained, np.int8).astype(bool),
                parents=parents, targets=targets, dup_alignments=dup_alignments,
                collisions=collisions, abutting=_abutting(
                    np.frombuffer(yielded, np.int64), np.frombuffer(start_keys, np.int64),
                    np.frombuffer(end_keys, np.int64)))


def _mix(parent, point):
    """One int64 key for (parent name hash, point); collisions only join trees."""
    return ((parent * 1_000_003 + point) & 0xFFFFFFFFFFFFFFFF) - (1 << 63)


def _abutting(lines, start_keys, end_keys):
    """Line pairs whose union puts rows ending where others start in one
    tree: each such row with the first row starting there, and the rows
    starting there with each other (linear in the rows at one point)."""
    order = np.argsort(start_keys, kind='stable')
    ordered = start_keys[order]
    low = np.searchsorted(ordered, end_keys, side='left')
    high = np.searchsorted(ordered, end_keys, side='right')
    pairs, groups = [], set()
    for row in np.flatnonzero(high > low).tolist():
        first = int(order[low[row]])
        if first != row or high[row] - low[row] > 1:
            pairs.append((int(lines[row]), int(lines[first])))
            groups.add(int(low[row]))
    for start in sorted(groups):
        stop = int(np.searchsorted(ordered, ordered[start], side='right'))
        pairs.extend((int(lines[order[start]]), int(lines[other]))
                     for other in order[start + 1:stop].tolist())
    return pairs


def _trees_and_needed(scan):
    """Union variant trees; collect names reachable variants need as roots."""
    parents, targets = scan['parents'], scan['targets']
    parent_of = {}

    def find(line):
        root = line
        while parent_of.get(root, root) != root:
            root = parent_of[root]
        while parent_of.get(line, line) != root:
            parent_of[line], line = root, parent_of[line]
        return root

    def union(line, other):
        a, b = find(line), find(other)
        if a != b:
            parent_of[max(a, b)] = min(a, b)

    templates = {}
    for line, names in targets.items():
        for name in names:
            other = parents.get(name)
            if other is not None:
                union(line, other)
            elif TEMPLATE_PATH.fullmatch(name):
                # A --exact duplication site and the rows on its template.
                union(line, templates.setdefault(name, line))
    # Rows that follow each other in a carrier's walk (gfa_junctions).
    for line, other in scan['abutting']:
        union(line, other)
    # Reachable = retained variants and everything they are built on.
    yielded, retained = scan['yielded'], scan['retained']
    needed, seen = set(), set()
    pending = [line for line in targets
               if retained[np.searchsorted(yielded, line)]]
    while pending:
        line = pending.pop()
        if line in seen:
            continue
        seen.add(line)
        for name in targets.get(line, ()):
            if name in parents:
                pending.append(parents[name])
            else:
                needed.add(name)
    return find, parent_of, needed


def _partitions(scan, find, parent_of, size):
    """Pack whole trees, in VCF order, into partitions of about ``size`` records."""
    yielded = scan['yielded']
    part = np.empty(len(yielded), np.int32)
    tree_part = {}
    current, count = 0, 0
    for number, line in enumerate(yielded.tolist()):
        root = find(line) if line in parent_of else line
        if root == line:
            if count >= size:
                current, count = current + 1, 0
            tree_part[line] = current
        part[number] = tree_part[root]
        count += 1
    return part, current + 1 if len(yielded) else 0


def _in_child(function, *arguments):
    """Run one partition step in a forked child, so all its memory is freed.

    Returns the function's result; re-raises its error. A child killed by the
    system (e.g. out of memory) is reported with its signal.
    """
    import multiprocessing as mp
    import traceback
    if 'fork' not in mp.get_all_start_methods():
        return function(*arguments)
    context = mp.get_context('fork')
    receive, send = context.Pipe(duplex=False)

    def target():
        try:
            send.send(('ok', function(*arguments)))
        except BaseException as error:  # reported to the parent
            try:
                send.send(('error', error, traceback.format_exc()))
            except Exception:
                send.send(('error', RuntimeError(str(error)), traceback.format_exc()))
        finally:
            send.close()

    child = context.Process(target=target)
    child.start()
    send.close()
    try:
        message = receive.recv()
    except EOFError:
        child.join()
        raise RuntimeError(f'partition worker exited with code {child.exitcode} '
                           '(a negative code is the signal; -9 is usually out of memory)')
    child.join()
    if message[0] == 'error':
        error = message[1]
        if isinstance(error, ValueError):
            raise ValueError(f'{error}') from None
        raise RuntimeError(f'partition worker failed:\n{message[2]}')
    return message[1]


# ---------------------------------------------------------------- partition pass

def _partition_pass(args, number, select, context, work, log):
    """Build one partition as in the one-pass run and save its state."""
    from gfa_duplications import apply_duplication_paths
    from gfa_stable_coords import StableResolver
    from gfa_interval_metadata import rows
    sources, roots, source_aliases = context['sources'], context['roots'], context['aliases']
    header_ordinals, header_lengths = context['header_ordinals'], context['header_lengths']
    folder = work / f'part{number:05d}'
    folder.mkdir()
    with tempfile.TemporaryFile(dir=work) as spool:
        candidates = QueryCandidates(spool)
        events = OrderedDict()
        kinds, record_indexes, sequence_checks = {}, {}, {}
        repeated = set()
        for record in rows(args.vcf, sources, args.size_cutoff, sequence_checks=sequence_checks,
                           candidate_sink=candidates.add,
                           insertion_only=getattr(args, 'insertion_only', None) is not None,
                           event_kinds=kinds, record_indexes=record_indexes,
                           select=select, index=context['index'], repeated=repeated):
            identifier, chrom, pos, ref_end, length, values, query, reason, keep = record
            if identifier in events:
                raise ValueError(f'duplicate insertion ID or variant ID {identifier!r}')
            kind = kinds.get(identifier, 'insertion')
            if getattr(args, 'variant_mode', 'all') == 'svonly' and kind == 'snp':
                keep = False
            # The global data line keeps VCF order across partitions.
            events[identifier] = Event(record_indexes[identifier], identifier, chrom, pos,
                                       ref_end, length, [Run(*run) for run in values],
                                       query, reason, keep, kind)
        kinds = None
        aliases = {name: f'ins_{header_ordinals[name]}' for name in events if name in header_ordinals}
        reachable = _reachable_events(events, include_parents=True)
        before = dict(context['dup_plan']['stats'])
        apply_duplication_paths(events, context['dup_plan'], Run)
        dup_delta = {key: value - before.get(key, 0)
                     for key, value in context['dup_plan']['stats'].items()}
        hits = {}
        resolver = StableResolver(events, roots, header_lengths, hits, reachable,
                                  args.nested_pos_base, source_aliases=source_aliases)
        resolver.index_variants()
        init_order = list(resolver.order)
        # Ranks are assigned in init order before any drop removes some.
        base_rank = resolver.ranks[init_order[0]] if init_order else 0
        stable_names = [aliases.get(name, name) for name in init_order]
        if len(stable_names) != len(set(stable_names)) or any(name in roots for name in stable_names):
            raise ValueError('insertion aliases collide with another insertion or FASTA root')
        for name in stable_names:
            _valid_name(name)
        intervals, missing = _prepare_intervals(events, reachable, sources, args.anchor,
                                                query_flanks=False)
        unresolved = folder / 'unresolved.tsv'
        drop = bool(getattr(args, 'drop_unverified', False))
        unrecoverable = [item for item in missing if item[0] not in candidates.positions]
        dropped = [] if drop else None
        if drop:
            dropped.extend(unrecoverable)
        else:
            _report_unresolved(unresolved, unrecoverable, aliases, log)
        intervals, extracted = _extract_queries(
            args, events, sources, intervals, missing, sequence_checks, candidates,
            folder, unresolved, aliases, log, dropped=dropped)
        hits.update(extracted)
        sequence_checks = candidates = None
        if drop:
            _drop_unverified(events, resolver, intervals, dropped, unresolved, aliases, log)
        unique_minimum = max(args.size_cutoff, args.unique_size_cutoff)
        _write_metadata(folder / 'anchors.bed', events, intervals, aliases, resolver,
                        unique_minimum, args.anchor,
                        mapping_prefix=lambda name: f'{record_indexes[name]}\t', sidecars=False)
        _write_variant_index(folder / 'variants.tsv', events, record_indexes, aliases, resolver)

        position = {name: offset for offset, name in enumerate(init_order)}
        specs, spec_event = [], array('q')
        for name in resolver.order:
            event = events[name]
            alias = aliases.get(name, name)
            specs.append(PathSpec(alias, event.kind, name, 0, event.length))
            spec_event.append(position[name])
            for run in event.runs:
                if (event.retained and run.operation in ('I', 'S') and run.unique and
                        run.qend - run.qstart >= unique_minimum):
                    specs.append(PathSpec(f'{alias}#unique{run.unique}', 'unique',
                                          name, run.qstart, run.qend))
                    spec_event.append(position[name])
        # Link-only walks write no P line, so they have no spec_event.
        from gfa_junctions import junction_specs
        specs.extend(junction_specs(args.vcf, events, resolver, aliases, PathSpec, log,
                                    select=select, index=context['index'], repeated=repeated))
        names = set()
        for spec in specs:
            _valid_name(spec.name)
            if spec.name in names or spec.name in roots:
                raise ValueError(f'duplicate GFA path name {spec.name!r}')
            names.add(spec.name)

        root_count = len(context['root_ids'])
        ids = dict(context['root_ids'])
        ids.update((name, root_count + offset) for offset, name in enumerate(init_order))
        leaves = topology.collect_leaves(specs, roots, resolver, args.anchor, ids,
                                         _path_leaves, _link_leaves, args.processes, log,
                                         paired_leaves=_path_and_link_leaves)
        ids = None
        keys = topology.boundaries(leaves, True, args.max_node_length)
        split = np.searchsorted(keys, root_count << topology.POINT_BITS)
        root_keys, event_keys = keys[:split], keys[split:]
        # Event pairs only ever come from this partition's own leaves.
        rows_on_events = leaves.source >= root_count
        first, last = topology.pair_ranges(event_keys, leaves.source[rows_on_events],
                                           leaves.start[rows_on_events], leaves.end[rows_on_events])
        link = leaves.link_row[rows_on_events] & (last > first)
        marks = (np.bincount(first[link], minlength=len(event_keys) + 1)
                 - np.bincount(last[link], minlength=len(event_keys) + 1))
        used_events = np.cumsum(marks)[:len(event_keys)] > 0
        counts = np.bincount((event_keys[used_events] >> topology.POINT_BITS) - root_count,
                             minlength=len(init_order))
        hit_rows = [(offset, hits[name]) for offset, name in enumerate(init_order) if name in hits]
        paths = sorted({hit.path for _offset, hit in hit_rows})
        path_number = {path: value for value, path in enumerate(paths)}
        hit_table = np.full((len(init_order), 5), -1, np.int64)
        for offset, hit in hit_rows:
            hit_table[offset] = (path_number[hit.path], hit.offset, hit.total, hit.left, hit.right)
        np.savez(folder / 'state.npz', source=leaves.source, start=leaves.start, end=leaves.end,
                 orient=leaves.orient, link_counts=leaves.link_counts,
                 path_counts=leaves.path_counts, ranks=leaves.ranks, root_keys=root_keys,
                 event_keys=event_keys, used_events=used_events, counts=counts,
                 order_lines=np.array([events[name].order for name in init_order], np.int64),
                 spec_event=np.frombuffer(spec_event, np.int64) if spec_event else np.zeros(0, np.int64),
                 hits=hit_table)
        with open(folder / 'state.pickle', 'wb') as handle:
            pickle.dump(dict(
                init_order=init_order, base_rank=base_rank,
                specs=[(spec.name, spec.kind, spec.source, spec.start, spec.end) for spec in specs],
                aliases={name: aliases[name] for name in init_order if name in aliases},
                hit_paths=paths, digits=_digits(names | set(stable_names)),
                name_hashes=[hash(name) for name in names | set(stable_names)]), handle,
                protocol=pickle.HIGHEST_PROTOCOL)
    log(f'Partition {number + 1}: {len(events)} variants, {len(specs)} paths, '
        f'{int(used_events.sum())} event segments')
    return folder, dup_delta


# ---------------------------------------------------------------- write pass

def _write_partition(args, folder, context, numbering, log):
    """Write this partition's S and P blocks and links with global numbers."""
    state = np.load(folder / 'state.npz')
    with open(folder / 'state.pickle', 'rb') as handle:
        meta = pickle.load(handle)
    roots = context['roots']
    root_count = len(context['root_ids'])
    init_order = meta['init_order']
    global_rank = numbering['ranks'][folder.name]            # per init position
    bases = numbering['bases'][folder.name]                   # first ordinal per event
    event_keys, used_events = state['event_keys'], state['used_events']
    event_of_pair = (event_keys >> topology.POINT_BITS) - root_count
    running = np.cumsum(used_events) - used_events            # used pairs before each pair
    starts = np.zeros(len(init_order) + 1, np.int64)
    np.cumsum(state['counts'], out=starts[1:])
    ordinals = bases[event_of_pair] + running - starts[event_of_pair]
    event_numbers = _skip(ordinals, numbering['forbidden'])

    keys = np.concatenate((numbering['root_keys'], event_keys))
    used = np.concatenate((np.zeros(len(numbering['root_keys']), bool), used_events))
    numbers = np.concatenate((numbering['root_numbers'], event_numbers))
    ranks = state['ranks'].copy()
    local = ranks >= meta['base_rank']
    ranks[local] = global_rank[ranks[local] - meta['base_rank']]
    leaves = topology.Leaves(state['source'], state['start'], state['end'], state['orient'],
                             state['link_counts'], state['path_counts'], ranks)
    specs = [PathSpec(*values) for values in meta['specs']]
    first, last = topology.pair_ranges(keys, leaves.source, leaves.start, leaves.end)
    covered = np.concatenate((numbering['root_used'], used_events))
    topology.check_path_coverage(covered, first, last, leaves, specs)

    names = context['root_names'] + init_order
    is_query = np.concatenate((np.zeros(root_count, bool), np.ones(len(init_order), bool)))
    segments = topology.Segments(keys, used, numbers, names, is_query)
    hit_table = state['hits']
    hits = {name: Hit(meta['hit_paths'][int(row[0])], *map(int, row[1:]))
            for name, row in zip(init_order, hit_table) if row[0] >= 0}
    rank_of = dict(zip(init_order, global_rank.tolist()))
    s_parts, p_parts, (left, right, link_ranks) = _write_parts(
        folder, specs, roots, leaves, first, last, numbers, segments, '', hits,
        meta['aliases'], rank_of, True, args.processes,
        preload_reference=getattr(args, 'preload_reference', True), log=log)
    for name, values in (('left', left), ('right', right), ('ranks', link_ranks)):
        np.save(folder / f'links.{name}.npy', values)
    # Per-event byte sizes let the blocks be copied in global rank order.
    s_bytes = _block_bytes(s_parts, folder / 'segments.gfa', state['counts'])
    p_count = np.bincount(state['spec_event'], minlength=len(init_order))
    p_bytes = _block_bytes(p_parts, folder / 'paths.gfa', p_count)
    np.savez(folder / 'blocks.npz', s_bytes=s_bytes, p_bytes=p_bytes)
    log(f'{folder.name}: wrote {len(segments)} segments, {len(specs)} paths, {len(left)} links')


def _block_bytes(parts, target, lines_per_block):
    """Concatenate parts into ``target``; return each block's byte size."""
    sizes = np.zeros(len(lines_per_block), np.int64)
    block, remaining = 0, None
    with open(target, 'wb') as output:
        for path in parts:
            with open(path, 'rb') as handle:
                for line in handle:
                    while remaining == 0 or remaining is None:
                        if remaining == 0:
                            block += 1
                        remaining = int(lines_per_block[block])
                    sizes[block] += len(line)
                    remaining -= 1
                    output.write(line)
            os.unlink(path)
    return sizes


# ---------------------------------------------------------------- driver

def run_partitioned(args, bed, log):
    from gfa_catalog_paths import resolve_catalog_paths
    from gfa_duplications import log_duplication_plan, plan_duplication_paths
    from gfa_interval_metadata import read_header_contigs, vcf_paths
    from gfa_source_catalog import resolve_local_sources
    from gfa_stable_coords import StableResolver

    if args.gfa_mode != 'rgfa' or args.bed_only or args.validate_only:
        raise ValueError('--partition-records needs --gfa-mode rgfa and a GFA output '
                         '(not --bed-only or --validate-only)')
    sources = read_query_list(args.query_fasta_list)
    header_ordinals, header_lengths = read_header_contigs(args.vcf)
    roots = _index_roots(args)
    log('Scanning VCF records for partitioning')
    scan = _scan(args, roots, log)
    find, parent_of, needed = _trees_and_needed(scan)
    part, count = _partitions(scan, find, parent_of, args.partition_records)
    log(f'Partitioned {len(part)} variants into {count} partition(s) of about '
        f'{args.partition_records} records (whole variant trees)')

    source_aliases = resolve_local_sources(
        {}, roots, (), sources, list(_flatten(args.local_path_fasta)), log,
        template_catalogs=args.local_reference_templates,
        backbone=getattr(args, 'reference_haplotype', None), needed=needed)
    from gfa_duplications import add_template_roots
    add_template_roots(needed, roots, source_aliases, Root, strict=False)
    catalog = getattr(args, 'alternative_catalog', None)
    if catalog:
        # DUP_ templates are not catalog records; their loci may be.
        source_aliases.update(resolve_catalog_paths(
            scan['collisions'], roots, (), catalog, list(vcf_paths(args.vcf)),
            getattr(args, 'reference_haplotype', None), args.reference_fasta, args.output or bed,
            Root, log, threads=args.processes,
            needed={name for name in needed if not TEMPLATE_PATH.fullmatch(name)}))
    add_template_roots(needed, roots, source_aliases, Root)
    plan = plan_duplication_paths(scan['dup_alignments'], roots, source_aliases, Root,
                                  None, scan['collisions'])
    dup_count = len(scan['dup_alignments'])
    root_resolver = StableResolver({}, roots, header_lengths, {}, set(),
                                   args.nested_pos_base, source_aliases=source_aliases)
    root_names = sorted(roots, key=lambda name: (root_resolver.ranks[name], name))
    context = dict(sources=sources, roots=roots, aliases=source_aliases,
                   header_ordinals=header_ordinals, header_lengths=header_lengths,
                   index=scan['index'], dup_plan=plan,
                   root_ids={name: number for number, name in enumerate(root_names)},
                   root_names=root_names)
    yielded = scan['yielded']
    del scan
    gc.collect()

    output = Path(args.output).resolve()
    output.parent.mkdir(parents=True, exist_ok=True)
    unresolved = Path(str(bed) + '.unresolved.tsv')
    with tempfile.TemporaryDirectory(prefix='gfa-partitions-', dir=bed.parent) as directory:
        work = Path(directory)
        folders = []
        for number in range(count):
            started = time.monotonic()
            select = np.sort(yielded[part == number])
            try:
                folder, dup_delta = _in_child(_partition_pass, args, number, select,
                                              context, work, log)
            except ValueError:
                failed = work / f'part{number:05d}' / 'unresolved.tsv'
                if failed.is_file():
                    shutil.copyfile(failed, unresolved)
                log(f'Stopped in partition {number + 1} of {count}; later partitions were '
                    'not checked')
                raise
            folders.append(folder)
            for key, value in dup_delta.items():
                plan['stats'][key] += value
            log(f'Partition {number + 1}/{count} done in {time.monotonic() - started:.1f}s')
        if dup_count:
            log_duplication_plan(plan, dup_count, log)
        numbering, root_part = _number(args, folders, context, root_resolver, log)
        for folder in folders:
            _in_child(_write_partition, args, folder, context, numbering, log)
        _assemble(args, output, folders, context, numbering, root_part, work, log)
        _merge_sidecars(bed, output, folders, context, log)
    return 0


def _data_lines(path):
    with open(path) as handle:
        next(handle, None)
        yield from handle


def _merge_sidecars(bed, output, folders, context, log):
    """Merge partition sidecars into the one-pass files and order."""
    def bed_key(line):
        f = line.rstrip('\n').split('\t')
        return f[6], f[0], int(f[1]), int(f[2]), f[3], f[5], int(f[7]), int(f[8])
    with open(bed, 'w', buffering=CHUNK) as out:
        out.write('#contig\tstart\tend\tpath\tscore\tstrand\tsample\t'
                  'insertion_start\tinsertion_end\n')
        out.writelines(heapq.merge(*[_data_lines(folder / 'anchors.bed') for folder in folders],
                                   key=bed_key))
    with open(str(bed) + '.mapping.tsv', 'w', buffering=CHUNK) as out:
        out.write('piece\tevent\tpiece_offset\tpiece_size\tsample\tcontig\tquery_start\t'
                  'query_end\tstrand\tref_contig\tref_start\tref_end\toperation\ttarget\t'
                  'target_start\ttarget_end\ttarget_strand\tref_strand\n')
        for line in heapq.merge(*[_data_lines(folder / 'anchors.bed.mapping.tsv')
                                  for folder in folders],
                                key=lambda line: int(line[:line.index('\t')])):
            out.write(line[line.index('\t') + 1:])
    with open(str(output) + '.variants.tsv', 'w', buffering=CHUNK) as out:
        out.write('#vcf_index\tvariant_id\tpath\tkind\tparent\tparent_start\t'
                  'parent_end\tparent_strand\n')
        out.writelines(heapq.merge(*[_data_lines(folder / 'variants.tsv') for folder in folders],
                                   key=lambda line: int(line[:line.index('\t')])))
    reported = 0
    with open(str(bed) + '.unresolved.tsv', 'w') as out:
        out.write('variant_id\tgfa_name\treason\n')
        for folder in folders:
            if (folder / 'unresolved.tsv').is_file():
                for line in _data_lines(folder / 'unresolved.tsv'):
                    out.write(line)
                    reported += 1
    _write_source_sidecars(bed, context['aliases'], context['roots'])
    if reported:
        log(f'{reported} variant(s) not exported are listed in {bed}.unresolved.tsv')


def _number(args, folders, context, root_resolver, log):
    """Number root segments, merge event ranks, and place event segment blocks."""
    roots = context['roots']
    root_ids = context['root_ids']
    specs = [PathSpec(root.name, root.kind, root.name, 0, root.length)
             for root in sorted((root for root in roots.values() if root.emit),
                                key=lambda value: value.order)]
    leaves = topology.collect_leaves(specs, roots, root_resolver, args.anchor, root_ids,
                                     _path_leaves, _link_leaves, args.processes, log,
                                     paired_leaves=_path_and_link_leaves)
    key_parts = [topology.boundaries(leaves, True, args.max_node_length)]
    forbidden = _digits(roots)
    hashes = [hash(name) for name in roots]
    sequences, metas = [], []
    for folder in folders:
        state = np.load(folder / 'state.npz')
        key_parts.append(state['root_keys'])
        with open(folder / 'state.pickle', 'rb') as handle:
            meta = pickle.load(handle)
        forbidden |= meta['digits']
        hashes.extend(meta['name_hashes'])
        sequences.append(state['order_lines'])
        metas.append((meta['base_rank'], state['counts']))
    hashes = np.sort(np.array(hashes, np.int64))
    if len(hashes) > 1 and (hashes[1:] == hashes[:-1]).any():
        raise ValueError('duplicate GFA path name across partitions; rerun without '
                         '--partition-records to see the name')
    root_keys = np.unique(np.concatenate(key_parts))
    # Root pairs used by any link leaf: root paths and every partition.
    marks = np.zeros(len(root_keys) + 1, np.int64)

    def mark(source, start, end, link_row):
        rows = link_row & (source < len(root_ids))
        first, last = topology.pair_ranges(root_keys, source[rows], start[rows], end[rows])
        valid = last > first
        marks[:] += (np.bincount(first[valid], minlength=len(marks))
                     - np.bincount(last[valid], minlength=len(marks)))
    mark(leaves.source, leaves.start, leaves.end, leaves.link_row)
    for folder in folders:
        state = np.load(folder / 'state.npz')
        part_leaves = topology.Leaves(state['source'], state['start'], state['end'],
                                      state['orient'], state['link_counts'],
                                      state['path_counts'], state['ranks'])
        mark(part_leaves.source, part_leaves.start, part_leaves.end, part_leaves.link_row)
    root_used = np.cumsum(marks)[:len(root_keys)] > 0
    root_numbers = _skip(np.cumsum(root_used), forbidden)
    total_root = int(root_used.sum())
    # Global dependency order: the one-pass heap always takes the smallest
    # VCF line among the partitions' next variants, i.e. heapq.merge.
    base_rank = max(root_resolver.ranks.values(), default=0) + 1
    def tagged(number, sequence):
        for offset, line in enumerate(sequence.tolist()):
            yield line, number, offset
    merged = heapq.merge(*[tagged(number, sequence) for number, sequence in enumerate(sequences)])
    ranks = [np.zeros(len(sequence), np.int64) for sequence in sequences]
    bases = [np.zeros(len(sequence), np.int64) for sequence in sequences]
    order = array('i')
    ordinal = total_root + 1
    for position, (_line, number, offset) in enumerate(merged):
        ranks[number][offset] = base_rank + position
        bases[number][offset] = ordinal
        ordinal += int(metas[number][1][offset])
        order.append(number)
    log(f'Numbered {total_root} root segments and {ordinal - total_root - 1} variant segments')
    numbering = dict(root_keys=root_keys, root_used=root_used, root_numbers=root_numbers,
                     forbidden=forbidden,
                     ranks={folder.name: ranks[number] for number, folder in enumerate(folders)},
                     bases={folder.name: bases[number] for number, folder in enumerate(folders)},
                     order=np.frombuffer(order, np.int32) if order else np.zeros(0, np.int32))
    root_part = (specs, leaves)
    return numbering, root_part


def _assemble(args, output, folders, context, numbering, root_part, work, log):
    """Header, root S, variant S blocks, links, root P, variant P blocks."""
    specs, leaves = root_part
    roots = context['roots']
    root_keys, root_numbers = numbering['root_keys'], numbering['root_numbers']
    first, last = topology.pair_ranges(root_keys, leaves.source, leaves.start, leaves.end)
    topology.check_path_coverage(numbering['root_used'], first, last, leaves, specs)
    segments = topology.Segments(root_keys, numbering['root_used'], root_numbers,
                                 context['root_names'],
                                 np.zeros(len(context['root_names']), bool))
    ranks = {name: 0 if roots[name].kind == 'reference' else 1 for name in roots}
    s_parts, p_parts, root_links = _write_parts(
        work, specs, roots, leaves, first, last, root_numbers, segments, '', {}, {},
        ranks, True, args.processes, tag='roots.',
        preload_reference=getattr(args, 'preload_reference', True), log=log)
    link_parts = [root_links] + [
        tuple(np.load(folder / f'links.{name}.npy', mmap_mode='r')
              for name in ('left', 'right', 'ranks')) for folder in folders]
    blocks = [np.load(folder / 'blocks.npz') for folder in folders]
    order = numbering['order']
    temporary = str(output) + f'.tmp.{os.getpid()}'
    try:
        import gzip
        opener = gzip.open if str(output).endswith('.gz') else open
        with opener(temporary, 'wb') as out:
            out.write(b'H\tVN:Z:1.0\tTS:Z:merged_vcf_to_gfa.py\n')
            _copy_parts(s_parts, out)
            _copy_blocks(out, folders, [block['s_bytes'] for block in blocks], order, 'segments.gfa')
            written = _merge_links(out, link_parts)
            _copy_parts(p_parts, out)
            _copy_blocks(out, folders, [block['p_bytes'] for block in blocks], order, 'paths.gfa')
        os.replace(temporary, output)
    finally:
        if os.path.exists(temporary):
            os.unlink(temporary)
    log(f'Wrote {written} links and the partitioned graph to {output}')


def _merge_links(out, parts, window=20_000_000):
    """Write the union of sorted link lists, ~``window`` links in RAM at a time.

    Each part is sorted by (left, right); links are processed in ranges of the
    left end, so one range from every part is reduced and written at a time.
    """
    total = sum(len(part[0]) for part in parts)
    if not total:
        return 0
    lows = [int(part[0][0]) for part in parts if len(part[0])]
    highs = [int(part[0][-1]) for part in parts if len(part[0])]
    low, high = min(lows), max(highs) + 1
    step = max(1, (high - low) * window // total)
    written = 0
    cursors = [0] * len(parts)
    for boundary in list(range(low + step, high, step)) + [high]:
        chunk = []
        for number, (left, right, ranks) in enumerate(parts):
            start = cursors[number]
            stop = int(np.searchsorted(left, boundary, side='left'))
            if stop > start:
                chunk.append((np.asarray(left[start:stop]), np.asarray(right[start:stop]),
                              np.asarray(ranks[start:stop])))
            cursors[number] = stop
        if chunk:
            left, right, ranks = topology._reduce_links(chunk)
            _write_links(out, left, right, ranks, '', True)
            written += len(left)
    return written


def _copy_blocks(out, folders, sizes, order, name):
    """Copy each partition's per-variant blocks in global rank order."""
    handles = [open(folder / name, 'rb') for folder in folders]
    cursor = [0] * len(folders)
    try:
        if not len(order):
            return
        change = np.flatnonzero(np.diff(order)) + 1
        starts = np.concatenate(([0], change))
        ends = np.concatenate((change, [len(order)]))
        for low, high in zip(starts.tolist(), ends.tolist()):
            number = int(order[low])
            count = high - low
            size = int(sizes[number][cursor[number]:cursor[number] + count].sum())
            cursor[number] += count
            while size:
                block = handles[number].read(min(CHUNK, size))
                if not block:
                    raise ValueError(f'{folders[number].name}/{name} is truncated')
                out.write(block)
                size -= len(block)
    finally:
        for handle in handles:
            handle.close()
