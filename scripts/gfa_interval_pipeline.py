"""Simple, query-backed GFA construction without a database."""
from bisect import bisect_left, bisect_right
from collections import Counter, defaultdict, OrderedDict
from dataclasses import dataclass
import gzip
import hashlib
import json
import os
from pathlib import Path
import re
import shutil
import tempfile

from gfa_query_anchors import interval_chunks, read_query_list, scaffold_pool
from minsetref_core import IndexedFasta

CHUNK = 1024 * 1024
VALID_SEQUENCE = re.compile(r'[A-Za-z=.]+')
VG_SEQUENCE = re.compile(r'[ACGTNacgtn]+')


@dataclass
class Run:
    qstart: int
    qend: int
    operation: str
    target: object
    rstart: int
    rend: int
    orientation: str
    unique: int


@dataclass
class Event:
    order: int
    identifier: str
    chrom: str
    pos: int
    ref_end: int
    length: int
    runs: list
    query: object
    query_reason: str
    retained: bool
    kind: str = 'insertion'


@dataclass
class Root:
    name: str
    path: str
    fai: object
    length: int
    emit: bool
    kind: str
    order: int
    lift: object = None


@dataclass
class Interval:
    event: str
    sample: str
    contig: str
    low: int
    high: int
    start: int
    end: int
    strand: str
    left: int
    right: int


@dataclass
class Hit:
    path: str
    offset: int
    total: int
    left: int
    right: int


@dataclass
class PathSpec:
    name: str
    kind: str
    source: str
    start: int
    end: int


class QueryCandidates:
    """Disk-backed alternatives; retain only a file offset and cursor per event."""
    def __init__(self, spool):
        self.spool = spool
        self.positions = {}

    def add(self, identifier, candidates):
        if len(candidates) > 1:
            self.positions[identifier] = (self.spool.tell(), 0)
            self.spool.write((json.dumps(candidates[1:], separators=(',', ':')) + '\n').encode())

    def next(self, identifier):
        position = self.positions.get(identifier)
        if position is None:
            return None
        offset, cursor = position
        self.spool.seek(offset)
        candidates = json.loads(self.spool.readline())
        value = tuple(candidates[cursor])
        if cursor + 1 == len(candidates):
            self.discard(identifier)
        else:
            self.positions[identifier] = (offset, cursor + 1)
        return value

    def discard(self, identifier):
        self.positions.pop(identifier, None)


def _valid_name(value):
    if (not value or value[0] in '*=' or ',' in value or
            any(not 33 <= ord(character) <= 126 for character in value)):
        raise ValueError(f'invalid GFA name: {value!r}')


def _tag(value):
    return str(value).replace('\t', ' ').replace('\r', ' ').replace('\n', ' ')


def _flatten(values):
    for group in values or ():
        for value in group:
            yield value


def _index_roots(args):
    from assembly_contigs import accepted_contig
    from gfa_source_catalog import indexed_header, source_interval
    stable = args.gfa_mode == 'rgfa'
    inputs = []
    if args.reference_fasta:
        inputs.append((args.reference_fasta, args.reference_fai, 'reference'))
    inputs.extend((value, None, 'alternative') for value in _flatten(args.alternatives_fasta))
    roots = {}
    emitted = 0
    seen_inputs = set()
    for value, fai, input_kind in inputs:
        path = str(Path(value).resolve())
        key = path, str(fai or ''), input_kind
        if key in seen_inputs:
            continue
        seen_inputs.add(key)
        reader = IndexedFasta(path, fai)
        try:
            for name, index in reader.index.items():
                length = index[0]
                fields = indexed_header(reader, name) if stable and length else []
                fixed = 'sequence_role=imported_alternative' in fields
                source = source_interval(fields, length) if fields else None
                if input_kind == 'reference' and not fixed and source and not getattr(args, 'reference_haplotype', None):
                    args.reference_haplotype = source[0]
                sample = source[0] if source else (getattr(args, 'reference_haplotype', None) if input_kind == 'reference' else None)
                contig = source[1] if source else name
                if not accepted_contig(sample, contig):
                    continue
                emit = bool(length and (stable or (input_kind != 'reference' and
                            name.startswith(('alternative', 'novel')) and
                            length >= args.size_cutoff)))
                kind = ('novel' if input_kind == 'alternative' and name.startswith('novel')
                        else input_kind)
                if fixed or (input_kind == 'reference' and source and source[0] != getattr(args, 'reference_haplotype', None)):
                    kind = 'novel' if name.startswith('novel') else 'alternative'
                previous = roots.get(name)
                if previous and (previous.path != path or previous.length != length):
                    raise ValueError(f'FASTA root {name!r} is defined more than once')
                if previous:
                    if emit and not previous.emit:
                        previous.emit = True
                        previous.kind = kind
                        previous.order = emitted
                        emitted += 1
                    continue
                roots[name] = Root(name, path, fai, length, emit, kind,
                                   emitted if emit else -1)
                if emit:
                    emitted += 1
        finally:
            reader.close()
    return roots


def _extract_batch(task):
    source, contig, requests, output = task
    reader = IndexedFasta(*source)
    hits = []
    failures = []
    try:
        with open(output, 'w+b', buffering=CHUNK) as sequences:
            for request in requests:
                event, low, high, start, end, strand = request[:6]
                checks = request[6] if len(request) > 6 else ()
                offset = sequences.tell()
                written = 0
                for sequence in interval_chunks(reader, contig, low, high, strand):
                    if not VALID_SEQUENCE.fullmatch(sequence):
                        raise ValueError(f'{contig}: sequence contains characters invalid in GFA')
                    data = sequence.encode('ascii')
                    sequences.write(data)
                    written += len(data)
                if written != high-low:
                    raise ValueError(f'{contig}:{low}-{high}: incomplete FASTA interval')
                left = start-low if strand == '+' else high-end
                mismatch = None
                for qstart, qend, expected in checks:
                    sequences.seek(offset+left+qstart)
                    remaining = qend-qstart
                    digest = hashlib.sha256()
                    while remaining:
                        block = sequences.read(min(CHUNK, remaining))
                        if not block:
                            raise ValueError(f'{event}: incomplete query sequence for validation')
                        digest.update(block.upper())
                        remaining -= len(block)
                    if digest.digest() != expected:
                        mismatch = (f'sequence_mismatch: selected query {contig}:{start}-{end}{strand} '
                                    f'does not match INFO/SEQ at insertion offsets '
                                    f'{qstart}-{qend}; SIZE alone does not identify '
                                    'the representative allele')
                        break
                if mismatch is not None:
                    # A failed candidate must not become a graph source. Reuse
                    # its temporary space for the next request in this batch.
                    sequences.seek(offset)
                    sequences.truncate()
                    failures.append((event, mismatch))
                    continue
                sequences.seek(offset+written)
                hits.append((event, output, offset, written, left,
                             written-left-(end-start)))
    finally:
        reader.close()
    return hits, failures


class Resolver:
    def __init__(self, events, roots, header_lengths, hits):
        self.events = events
        self.roots = roots
        self.header_lengths = header_lengths
        self.hits = hits

    def target_length(self, name, required=True):
        if name in self.events:
            return self.events[name].length
        if name in self.roots:
            return self.roots[name].length
        length = self.header_lengths.get(name)
        if length is None and required:
            raise ValueError(f'cannot determine coordinate length of graph target {name!r}')
        return length

    def target_interval(self, name, start, end, orientation):
        length = self.target_length(name, required=orientation == '-')
        low, high = ((start, end) if orientation == '+' else
                     (length-end, length-start))
        if low < 0 or high < low or (length is not None and high > length):
            raise ValueError(f'{name}: mapped interval {low}-{high} is out of bounds')
        return low, high

    @staticmethod
    def reverse(leaves):
        return [(kind, source, start, end, '-' if orientation == '+' else '+')
                for kind, source, start, end, orientation in reversed(leaves)]

    def expand(self, name, start, end, orientation='+', stack=()):
        event = self.events.get(name)
        if event is None:
            leaves = [('root', name, start, end, '+')]
            return leaves if orientation == '+' else self.reverse(leaves)
        if name in stack:
            raise ValueError('cyclic insertion graph: ' + ' -> '.join((*stack, name)))
        if not 0 <= start <= end <= event.length:
            raise ValueError(f'{name}: interval {start}-{end} is outside 0-{event.length}')
        leaves = []
        covered = start
        for run in event.runs:
            if run.qend <= start or run.qstart >= end:
                continue
            low, high = max(start, run.qstart), min(end, run.qend)
            if low != covered:
                raise ValueError(f'{name}: graph definition has a gap at {covered}')
            if run.operation == '=':
                target_low, target_high = self.target_interval(
                    run.target, run.rstart, run.rend, run.orientation)
                if run.orientation == '+':
                    a = target_low + low-run.qstart
                    b = target_low + high-run.qstart
                else:
                    a = target_high - (high-run.qstart)
                    b = target_high - (low-run.qstart)
                leaves.extend(self.expand(run.target, a, b, run.orientation,
                                          (*stack, name)))
            else:
                hit = self.hits.get(name)
                if hit is None:
                    raise ValueError(f'{name}: query sequence is required for {run.operation}')
                leaves.append(('query', name, hit.left+low, hit.left+high, '+'))
            covered = high
        if covered != end:
            raise ValueError(f'{name}: graph definition ends at {covered}, expected {end}')
        return leaves if orientation == '+' else self.reverse(leaves)

    def local_path(self, name, start, end, flank):
        event = self.events[name]
        core_start, core_end = max(0, start-flank), min(event.length, end+flank)
        leaves = self.expand(name, core_start, core_end)
        hit = self.hits.get(name)
        if hit:
            left = min(hit.left, max(0, flank-start))
            right = min(hit.right, max(0, flank-(event.length-end)))
            if left:
                leaves.insert(0, ('query', name, hit.left-left, hit.left, '+'))
            if right:
                boundary = hit.left + event.length
                leaves.append(('query', name, boundary, boundary+right, '+'))
        return leaves


def _reachable_events(events, include_parents=False):
    reachable = set()
    pending = [event.identifier for event in events.values() if event.retained]
    while pending:
        identifier = pending.pop()
        if identifier in reachable:
            continue
        reachable.add(identifier)
        event = events[identifier]
        if include_parents and event.chrom in events:
            pending.append(event.chrom)
        for run in event.runs:
            if (run.operation == '=' or include_parents) and run.target in events:
                pending.append(run.target)
    return reachable


def _query_coordinate(event, run):
    if event.query is None:
        return None
    sample, contig, start, end, strand = event.query
    if strand == '+':
        low, high = start+run.qstart, start+run.qend
    else:
        low, high = end-run.qend, end-run.qstart
    return sample, contig, low, high, strand


def _prepare_intervals(events, reachable, sources, anchor, query_flanks=True):
    needed = []
    literal = set()
    failures = {}
    for identifier in reachable:
        event = events[identifier]
        has_literal = any(run.operation in ('X', 'I', 'S') for run in event.runs)
        if has_literal:
            literal.add(identifier)
            if event.query is None:
                failures[identifier] = event.query_reason
        if event.query and (has_literal or (query_flanks and event.retained and anchor)):
            needed.append(event)
    needed.sort(key=lambda event: (event.query[0], event.query[1],
                                   event.query[2], event.query[3], event.order))
    intervals = {}
    current_sample = None
    reader = None
    try:
        for event in needed:
            sample, contig, start, end, strand = event.query
            if sample != current_sample:
                if reader:
                    reader.close()
                reader = IndexedFasta(*sources[sample])
                current_sample = sample
            if contig not in reader.index:
                failures[event.identifier] = (
                    f'scaffold_missing: query scaffold {contig!r} is absent from FASTA '
                    f'for sample {sample!r}')
                continue
            contig_length = reader.index[contig][0]
            if not 0 <= start < end <= contig_length:
                failures[event.identifier] = (
                    f'coordinate_out_of_bounds: QUERYCOORD {contig}:{start}-{end} is outside '
                    f'scaffold length {contig_length}')
                continue
            flank = anchor if event.retained else 0
            low, high = max(0, start-flank), min(reader.index[contig][0], end+flank)
            left = start-low if strand == '+' else high-end
            right = high-end if strand == '+' else start-low
            intervals[event.identifier] = Interval(event.identifier, sample, contig,
                                                   low, high, start, end, strand,
                                                   left, right)
    finally:
        if reader:
            reader.close()
    missing = [(value, failures.get(value, 'query interval was not prepared'))
               for value in sorted(literal-set(intervals),
                                   key=lambda item: events[item].order)]
    return intervals, missing


def _report_unresolved(path, failures, aliases, log):
    with path.open('w') as output:
        output.write('variant_id\tgfa_name\treason\n')
        for identifier, reason in failures:
            output.write(f'{identifier}\t{aliases.get(identifier, identifier)}\t{reason}\n')
    if failures:
        for identifier, reason in failures[:20]:
            log(f'Unresolved variant {identifier}: {reason}')
        if len(failures) > 20:
            log(f'{len(failures)-20} additional unresolved variants are listed in {path}')
        categories = Counter(reason.split(':', 1)[0] for _identifier, reason in failures)
        for reason, count in categories.most_common():
            log(f'Unresolved query intervals: {reason}={count}')
        first, reason = failures[0]
        raise ValueError(f'{len(failures)} reachable insertions lack usable query intervals; '
                         f'first unresolved variant is {first!r}: {reason}; see {path}')


def _extract_queries(args, events, sources, intervals, missing, sequence_checks,
                     candidates, temporary_root, unresolved, aliases, log):
    """Retry failed candidates in scaffold batches, preserving successful reads."""
    verified_intervals, hits = {}, {}
    attempts = Counter()
    round_number = 0
    _report_unresolved(unresolved, [], aliases, log)
    while intervals or missing:
        round_number += 1
        for identifier in set(intervals) | {name for name, _reason in missing}:
            attempts[identifier] += int(events[identifier].query is not None)
        sorted_intervals = sorted(intervals.values(), key=lambda value: (
            value.sample, value.contig, value.low, value.high, events[value.event].order))
        tasks, batch = [], []
        group = None
        def submit():
            tasks.append((sources[group[0]], group[1], list(batch),
                          str(temporary_root / f'{round_number}-{len(tasks)}.sequence')))
        for interval in sorted_intervals:
            key = interval.sample, interval.contig
            if key != group or len(batch) >= 2048:
                if batch:
                    submit()
                group, batch = key, []
            batch.append((interval.event, interval.low, interval.high,
                          interval.start, interval.end, interval.strand,
                          sequence_checks.get(interval.event, ()) if sequence_checks is not None else ()))
        if batch:
            submit()

        failures = dict(missing)
        if tasks:
            log(f'Extracting {len(intervals)} sorted query intervals with {args.processes} workers '
                f'(candidate round {round_number})')
            with scaffold_pool(min(args.processes, len(tasks))) as pool:
                for completed, (result, rejected) in enumerate(pool.imap_unordered(
                        _extract_batch, tasks, chunksize=1), 1):
                    for identifier, path, offset, total, left, right in result:
                        hits[identifier] = Hit(path, offset, total, left, right)
                        verified_intervals[identifier] = intervals[identifier]
                        if candidates is not None:
                            candidates.discard(identifier)
                    failures.update(rejected)
                    if completed % 20 == 0 or completed == len(tasks):
                        log(f'Extracted {completed}/{len(tasks)} scaffold batches')
        if set(intervals) - set(hits) - set(failures):
            raise ValueError('one or more query intervals were not extracted')
        retry, exhausted = set(), []
        for identifier in sorted(failures, key=lambda name: events[name].order):
            candidate = candidates.next(identifier) if candidates is not None else None
            if candidate is None:
                reason = failures[identifier]
                if candidates is not None:
                    reason += f'; exhausted {attempts[identifier]} eligible query candidate(s)'
                exhausted.append((identifier, reason))
            else:
                events[identifier].query = candidate
                retry.add(identifier)
        _report_unresolved(unresolved, exhausted, aliases, log)
        if not retry:
            break
        log(f'Retrying {len(retry)} insertions with later SIZE/span-matched observations')
        intervals, missing = _prepare_intervals(events, retry, sources, args.anchor,
                                               query_flanks=args.gfa_mode != 'rgfa')
    return verified_intervals, hits


def _write_metadata(bed, events, intervals, aliases, resolver, unique_minimum, anchor):
    bed_rows = []
    for event in events.values():
        interval = intervals.get(event.identifier)
        if not event.retained or interval is None:
            continue
        alias = aliases.get(event.identifier, event.identifier)
        for run in event.runs:
            if run.operation not in ('I', 'S'):
                continue
            if run.unique and run.qend-run.qstart < unique_minimum:
                continue
            coordinate = _query_coordinate(event, run)
            sample, contig, low, high, strand = coordinate
            name = alias if run.unique == 0 else f'{alias}#unique{run.unique}'
            bed_rows.append((sample, contig, max(interval.low, low-anchor),
                             min(interval.high, high+anchor), name, strand, low, high))
    bed_rows.sort()
    with open(bed, 'w', buffering=CHUNK) as out:
        out.write('#contig\tstart\tend\tpath\tscore\tstrand\tsample\t'
                  'insertion_start\tinsertion_end\n')
        for sample, contig, start, end, name, strand, low, high in bed_rows:
            out.write(f'{contig}\t{start}\t{end}\t{name}\t0\t{strand}\t{sample}\t'
                      f'{low}\t{high}\n')

    with open(str(bed)+'.mapping.tsv', 'w', buffering=CHUNK) as out:
        out.write('piece\tevent\tpiece_offset\tpiece_size\tsample\tcontig\tquery_start\t'
                  'query_end\tstrand\tref_contig\tref_start\tref_end\toperation\ttarget\t'
                  'target_start\ttarget_end\ttarget_strand\tref_strand\n')
        for event in events.values():
            if hasattr(resolver, 'ranks') and event.identifier not in resolver.ranks:
                continue
            for number, run in enumerate(event.runs):
                coordinate = _query_coordinate(event, run)
                if coordinate:
                    sample, contig, qstart, qend, strand = coordinate
                else:
                    sample = contig = qstart = qend = strand = '.'
                if run.target:
                    tstart, tend = resolver.target_interval(
                        run.target, run.rstart, run.rend, run.orientation)
                    target, tstrand = run.target, run.orientation
                    if hasattr(resolver, 'canonical_interval'):
                        target, tstart, tend, tstrand = resolver.canonical_interval(
                            target, tstart, tend, tstrand)
                    target = aliases.get(target, target)
                else:
                    target = tstart = tend = '.'
                    tstrand = run.orientation
                event_name = aliases.get(event.identifier, event.identifier)
                piece = (f'{event_name}#unique{run.unique}' if run.unique else event_name)
                ref_start, ref_end = (resolver.breakpoints(event.identifier)
                                      if hasattr(resolver, 'ranks') else
                                      (event.pos, event.ref_end))
                ref_contig, ref_strand = event.chrom, '+'
                if hasattr(resolver, 'canonical_interval'):
                    ref_contig, ref_start, ref_end, ref_strand = resolver.canonical_interval(
                        ref_contig, ref_start, ref_end)
                out.write(f'{piece}\t{event_name}\t{run.qstart}\t{run.qend-run.qstart}\t'
                          f'{sample}\t{contig}\t{qstart}\t{qend}\t{strand}\t'
                          f'{aliases.get(ref_contig,ref_contig)}\t{ref_start}\t{ref_end}\t'
                          f'{run.operation}\t{target}\t{tstart}\t{tend}\t{tstrand}\t{ref_strand}\n')
    if hasattr(resolver, 'source_aliases'):
        with open(str(bed)+'.local-paths.tsv', 'w') as out:
            out.write('local_path\ttarget\tstart\tend\tstrand\n')
            for name, alias in sorted(resolver.source_aliases.items()):
                out.write(f'{name}\t{alias.target}\t{alias.start}\t{alias.end}\t{alias.strand}\n')
        with open(str(bed)+'.template-lifts.tsv', 'w') as out:
            out.write('template\tbackbone\tstatus\tplacements\tnote\n')
            for name, root in sorted(resolver.roots.items()):
                if root.lift:
                    lift = root.lift
                    out.write(f'{name}\t{lift["backbone"]}\t{lift["status"]}\t{lift["placements"]}\t{lift["note"]}\n')


def _path_leaves(spec, roots, resolver, anchor):
    if spec.kind in ('reference', 'alternative', 'novel'):
        return resolver.expand(spec.source, spec.start, spec.end)
    if hasattr(resolver, 'ranks') and spec.kind == 'deletion':
        return resolver.local_path(spec.source, spec.start, spec.end, max(1, anchor))
    if hasattr(resolver, 'ranks') and spec.kind in ('insertion', 'snp', 'substitution'):
        # Keep the named insertion path at its original 0-based core origin.
        # Parent flanks are used to form graph edges, not prepended to this P.
        return resolver.expand(spec.source, spec.start, spec.end)
    return resolver.local_path(spec.source, spec.start, spec.end, anchor)


def _link_leaves(spec, roots, resolver, anchor):
    if hasattr(resolver, 'ranks') and spec.source in roots and roots[spec.source].lift:
        return resolver.local_path(spec.source, spec.start, spec.end, max(1, anchor))
    if hasattr(resolver, 'ranks') and spec.kind in ('insertion', 'snp', 'substitution', 'deletion'):
        # Even --anchor 0 must attach the allele to its parent graph.
        return resolver.local_path(spec.source, spec.start, spec.end, max(1, anchor))
    return _path_leaves(spec, roots, resolver, anchor)


def _atomic_keys(leaf, boundaries):
    kind, source, start, end, orientation = leaf
    points = boundaries[kind, source]
    first, last = bisect_left(points, start), bisect_right(points, end)
    selected = points[first:last]
    pairs = list(zip(selected, selected[1:]))
    if orientation == '-':
        pairs.reverse()
    for low, high in pairs:
        yield (kind, source, low, high), orientation


class _OpenFiles:
    def __init__(self, maximum=32):
        self.maximum = maximum
        self.files = OrderedDict()

    def get(self, path):
        handle = self.files.pop(path, None)
        if handle is None:
            handle = open(path, 'rb')
        self.files[path] = handle
        while len(self.files) > self.maximum:
            self.files.popitem(last=False)[1].close()
        return handle

    def close(self):
        for handle in self.files.values():
            handle.close()
        self.files.clear()


def _copy(source, offset, size, output, uppercase=False):
    source.seek(offset)
    while size:
        block = source.read(min(CHUNK, size))
        if not block:
            raise ValueError('truncated temporary query sequence')
        if uppercase and not VG_SEQUENCE.fullmatch(block.decode('ascii')):
            raise ValueError('query sequence contains bases outside ACGTN; VG would '
                             'change these bases on import')
        output.write(block.upper() if uppercase else block)
        size -= len(block)


def _write_segments(output, segment_ids, hits, roots, aliases, prefix, ranks=None):
    query_files = _OpenFiles()
    root_reader = None
    root_path = None
    try:
        for (kind, source, start, end), number in segment_ids.items():
            output.write(f'S\t{prefix}{number}\t'.encode('ascii'))
            if kind == 'query':
                hit = hits[source]
                if ranks is not None and not hit.left <= start < end <= hit.total-hit.right:
                    raise ValueError(f'{source}: query flank cannot have an insertion coordinate')
                _copy(query_files.get(hit.path), hit.offset+start, end-start, output,
                      uppercase=ranks is not None)
                tags = (f'LN:i:{end-start}\tSN:Z:{_tag(aliases.get(source,source))}\t'
                        f'SO:i:{start-hit.left}\tTP:Z:query')
            else:
                root = roots.get(source)
                if root is None:
                    output.write(b'*')
                    root_kind = 'unknown'
                else:
                    if root.path != root_path:
                        if root_reader:
                            root_reader.close()
                        root_reader = IndexedFasta(root.path, root.fai)
                        root_path = root.path
                    written = 0
                    for sequence in interval_chunks(root_reader, source, start, end, '+'):
                        if not VALID_SEQUENCE.fullmatch(sequence):
                            raise ValueError(f'{source}: sequence contains characters invalid in GFA')
                        if ranks is not None and not VG_SEQUENCE.fullmatch(sequence):
                            raise ValueError(f'{source}: sequence contains bases outside ACGTN; '
                                             'VG would change these bases on import')
                        output.write((sequence.upper() if ranks is not None else sequence).encode('ascii'))
                        written += len(sequence)
                    if written != end-start:
                        raise ValueError(f'{source}:{start}-{end}: incomplete FASTA interval')
                    root_kind = root.kind
                tags = (f'LN:i:{end-start}\tSN:Z:{_tag(source)}\tSO:i:{start}\t'
                        f'TP:Z:{root_kind}')
            if ranks is not None:
                tags += f'\tSR:i:{ranks[source]}'
            output.write(('\t' + tags + '\n').encode('ascii'))
    finally:
        query_files.close()
        if root_reader:
            root_reader.close()


def _write_graph(output_path, temporary_root, specs, roots, resolver, anchor,
                 boundaries, segment_ids, prefix, hits, aliases):
    path_file = temporary_root / 'paths.gfa'
    links = {}
    ranks = getattr(resolver, 'ranks', None)
    with open(path_file, 'wb', buffering=CHUNK) as output:
        for spec in specs:
            previous = None
            rank = ranks[spec.source] if ranks is not None else 0
            for leaf in _link_leaves(spec, roots, resolver, anchor):
                for key, orientation in _atomic_keys(leaf, boundaries):
                    token = segment_ids[key], orientation
                    if previous:
                        edge = previous[0], previous[1], token[0], token[1]
                        if ranks is not None:
                            reverse_edge = (edge[2], '-' if edge[3] == '+' else '+',
                                            edge[0], '-' if edge[1] == '+' else '+')
                            edge = min(edge, reverse_edge)
                        links[edge] = min(links.get(edge, rank), rank)
                    previous = token
            output.write(f'P\t{spec.name}\t'.encode('ascii'))
            first = True
            for leaf in _path_leaves(spec, roots, resolver, anchor):
                for key, orientation in _atomic_keys(leaf, boundaries):
                    token = segment_ids[key], orientation
                    if not first:
                        output.write(b',')
                    output.write(f'{prefix}{token[0]}{token[1]}'.encode('ascii'))
                    first = False
            if first:
                raise ValueError(f'{spec.name}: empty GFA path')
            lift = roots[spec.source].lift if spec.source in roots else None
            lift_tag = f'\tLS:Z:{lift["status"]}' if lift else ''
            if lift and lift['status'] not in ('mapped', 'one_sided'):
                lift_tag += '\tUP:Z:unplaced'
            output.write(f'\t*\tTP:Z:{spec.kind}{lift_tag}\n'.encode('ascii'))

    temporary = str(output_path) + f'.tmp.{os.getpid()}'
    try:
        opener = gzip.open if str(output_path).endswith('.gz') else open
        with opener(temporary, 'wb') as output:
            output.write(b'H\tVN:Z:1.0\tTS:Z:merged_vcf_to_gfa.py\n')
            _write_segments(output, segment_ids, hits, roots, aliases, prefix, ranks)
            for edge in sorted(links):
                left, left_orientation, right, right_orientation = edge
                tags = f'\tSR:i:{links[edge]}' if ranks is not None else ''
                output.write(f'L\t{prefix}{left}\t{left_orientation}\t{prefix}{right}\t'
                             f'{right_orientation}\t0M{tags}\n'.encode('ascii'))
            with open(path_file, 'rb') as paths:
                shutil.copyfileobj(paths, output, length=CHUNK)
        os.replace(temporary, output_path)
    finally:
        if os.path.exists(temporary):
            os.unlink(temporary)
    return len(links)


def run(args, _read_vcf=None, _chunks=None, _operations=None, log=print):
    from gfa_interval_metadata import vcf_paths
    bed = Path(args.anchor_bed or (str(args.output or next(vcf_paths(args.vcf))) + '.anchors.bed')).resolve()
    bed.parent.mkdir(parents=True, exist_ok=True)
    with tempfile.TemporaryFile(dir=bed.parent) as spool:
        return _run(args, QueryCandidates(spool), bed, log)


def _run(args, candidates, bed, log):
    from gfa_interval_metadata import read_header_contigs, rows

    sources = read_query_list(args.query_fasta_list)
    header_ordinals, header_lengths = read_header_contigs(args.vcf)
    roots = _index_roots(args)
    sequence_checks = {} if args.gfa_mode == 'rgfa' else None
    log('Streaming VCF metadata and SIZE/span-matched query intervals')
    events = OrderedDict()
    event_kinds = {}
    for order, record in enumerate(rows(args.vcf, sources, args.size_cutoff,
                                       sequence_checks=sequence_checks,
                                       candidate_sink=candidates.add if sequence_checks is not None else None,
                                       insertion_only=getattr(args, 'insertion_only', None) is not None,
                                       event_kinds=event_kinds)):
        identifier, chrom, pos, ref_end, length, values, query, query_reason, retained = record
        if identifier in events:
            raise ValueError(f'duplicate insertion ID or variant ID {identifier!r}')
        kind = event_kinds.get(identifier, 'insertion')
        if getattr(args, 'variant_mode', 'all') == 'svonly' and kind == 'snp':
            retained = False
        if kind != 'insertion' and args.gfa_mode == 'query':
            raise ValueError('SNP/deletion/substitution branches require --gfa-mode rgfa; '
                             'use --insertion-only for the legacy query mode')
        events[identifier] = Event(order, identifier, chrom, pos, ref_end, length,
                                   [Run(*run) for run in values], query, query_reason, retained, kind)
        if (order+1) % 10000 == 0:
            log(f'Indexed {order+1} insertion definitions')

    aliases = {name: f'ins_{header_ordinals[name]}'
               for name in events if name in header_ordinals}
    stable = args.gfa_mode == 'rgfa'
    reachable = _reachable_events(events, include_parents=stable)
    hits = {}
    if stable:
        from gfa_source_catalog import resolve_local_sources
        from gfa_stable_coords import StableResolver
        source_aliases = resolve_local_sources(events, roots, reachable, sources,
                                              list(_flatten(args.local_path_fasta)), log,
                                              template_catalogs=args.local_reference_templates,
                                              backbone=getattr(args, 'reference_haplotype', None))
        resolver = StableResolver(events, roots, header_lengths, hits, reachable,
                                  args.nested_pos_base, source_aliases=source_aliases)
        stable_names = [aliases.get(name, name) for name in resolver.order] + list(roots)
        if len(stable_names) != len(set(stable_names)):
            raise ValueError('insertion aliases collide with another insertion or FASTA root')
        for name in stable_names:
            _valid_name(name)
        stable_names = set(stable_names)
    else:
        resolver = Resolver(events, roots, header_lengths, hits)
    intervals, missing = _prepare_intervals(events, reachable, sources, args.anchor,
                                          query_flanks=not stable)
    unresolved = Path(str(bed)+'.unresolved.tsv')
    unique_minimum = max(args.size_cutoff, args.unique_size_cutoff)
    # Reject known-unrecoverable metadata before reading any assembly sequence.
    _report_unresolved(unresolved, [item for item in missing
                       if not stable or item[0] not in candidates.positions], aliases, log)
    if args.bed_only:
        _report_unresolved(unresolved, missing, aliases, log)
        _write_metadata(bed, events, intervals, aliases, resolver,
                        unique_minimum, args.anchor)
        log(f'Wrote {bed} and {bed}.mapping.tsv (query sequences not verified)')
        return 0

    with tempfile.TemporaryDirectory(prefix='gfa-build-', dir=bed.parent) as directory:
        temporary_root = Path(directory)
        intervals, extracted = _extract_queries(
            args, events, sources, intervals, missing, sequence_checks,
            candidates if stable else None, temporary_root, unresolved, aliases, log)
        hits.update(extracted)
        _write_metadata(bed, events, intervals, aliases, resolver,
                        unique_minimum, args.anchor)

        specs = []
        for root in sorted((root for root in roots.values() if root.emit),
                           key=lambda value: value.order):
            specs.append(PathSpec(root.name, root.kind, root.name, 0, root.length))
        ordered_events = ([events[name] for name in resolver.order] if stable
                          else events.values())
        for event in ordered_events:
            if not stable and not event.retained:
                continue
            alias = aliases.get(event.identifier, event.identifier)
            specs.append(PathSpec(alias, event.kind, event.identifier, 0, event.length))
            for run in event.runs:
                if (event.retained and run.operation in ('I', 'S') and run.unique and
                        run.qend-run.qstart >= unique_minimum):
                    specs.append(PathSpec(f'{alias}#unique{run.unique}', 'unique',
                                          event.identifier, run.qstart, run.qend))

        names = set()
        for spec in specs:
            _valid_name(spec.name)
            if spec.name in names:
                raise ValueError(f'duplicate GFA path name {spec.name!r}')
            names.add(spec.name)
        prefix = '__gfa_segment_'
        while any(name.startswith(prefix) for name in names):
            prefix = '_' + prefix

        log(f'Building shared topology for {len(specs)} GFA paths')
        cuts = defaultdict(set)
        for spec in specs:
            getters = ((_path_leaves, _link_leaves) if stable and (spec.kind in ('insertion', 'snp', 'substitution', 'deletion') or
                        (spec.source in roots and roots[spec.source].lift))
                       else (_path_leaves,))
            for get_leaves in getters:
                for kind, source, start, end, _orientation in get_leaves(
                        spec, roots, resolver, args.anchor):
                    cuts[kind, source].update((start, end))
                    if stable:
                        # Chop only covered intervals, never large unused gaps
                        # between disjoint pieces of a query source.
                        cuts[kind, source].update(range(start+args.max_node_length, end,
                                                       args.max_node_length))
        boundaries = {source: sorted(points) for source, points in cuts.items()}
        segment_keys = set()
        for spec in specs:
            for leaf in _link_leaves(spec, roots, resolver, args.anchor):
                segment_keys.update(key for key, _orientation in
                                    _atomic_keys(leaf, boundaries))
        if stable:
            prefix = ''  # Preserve numeric node IDs through VG import/export.
            segment_ids = {}
            number = 0
            for key in sorted(segment_keys, key=lambda k: (resolver.ranks[k[1]], k)):
                number += 1
                while str(number) in names or str(number) in stable_names:
                    number += 1
                segment_ids[key] = number
        else:
            segment_ids = {key: number for number, key in enumerate(sorted(segment_keys), 1)}
        if args.validate_only:
            if stable:
                class Discard:
                    def write(self, _data):
                        pass
                _write_segments(Discard(), segment_ids, hits, roots, aliases, prefix,
                                resolver.ranks)
            log(f'Validated {len(segment_ids)} segments and {len(specs)} paths')
            return 0

        output = Path(args.output).resolve()
        output.parent.mkdir(parents=True, exist_ok=True)
        log('Writing GFA segments, links, and paths')
        links = _write_graph(output, temporary_root, specs, roots, resolver, args.anchor,
                             boundaries, segment_ids, prefix, hits, aliases)
        log(f'Wrote {len(segment_ids)} segments, {links} links, and '
            f'{len(specs)} paths to {output}')
        return 0
