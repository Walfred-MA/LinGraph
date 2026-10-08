"""Simple, query-backed GFA construction without a database."""
from collections import Counter, defaultdict, OrderedDict
from contextlib import ExitStack
from dataclasses import dataclass
import gzip
import hashlib
import json
import os
from pathlib import Path
import re
import shutil
import tempfile
import time

from gfa_query_anchors import interval_chunks, read_query_list, reverse_complement, scaffold_pool
from minsetref_core import IndexedFasta

CHUNK = 1024 * 1024
ROOT_WINDOW = 4 * CHUNK
QUERY_WINDOW = 64 * 1024
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
    # A duplication path copies [base, base+length) of FASTA record ``record``.
    record: str = None
    base: int = 0


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
                if not fixed and not accepted_contig(sample, contig):
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
        window_start, window = 0, ''
        # Requests are sorted and bounded to one scaffold. Limit read-ahead to
        # this batch, so sparse requests do not read past its final interval.
        window_limit = max((request[2] for request in requests), default=0)
        with open(output, 'w+b', buffering=CHUNK) as sequences:
            for request in requests:
                event, low, high, start, end, strand = request[:6]
                checks = request[6] if len(request) > 6 else ()
                offset = sequences.tell()
                if high - low <= QUERY_WINDOW:
                    if not (window_start <= low and high <= window_start + len(window)):
                        window_start = low
                        window = reader.fetch(contig, low, min(window_limit, low + QUERY_WINDOW))
                    sequence = window[low-window_start:high-window_start]
                    if strand != '+':
                        sequence = reverse_complement(sequence)
                    if not VALID_SEQUENCE.fullmatch(sequence):
                        raise ValueError(f'{contig}: sequence contains characters invalid in GFA')
                    data = sequence.encode('ascii')
                    if len(data) != high-low:
                        raise ValueError(f'{contig}:{low}-{high}: incomplete FASTA interval')
                    left = start-low if strand == '+' else high-end
                    # Validate short alleles directly before buffering their
                    # bytes. Seeking back into w+b after every SNP flushes the
                    # write buffer and turns millions of SNPs into tiny I/Os.
                    mismatch = None
                    for qstart, qend, expected in checks:
                        block = data[left+qstart:left+qend]
                        if len(block) != qend-qstart:
                            raise ValueError(f'{event}: incomplete query sequence for validation')
                        if hashlib.sha256(block.upper()).digest() != expected:
                            mismatch = _sequence_mismatch(contig, start, end, strand, qstart, qend)
                            break
                    if mismatch is not None:
                        failures.append((event, mismatch))
                        continue
                    sequences.write(data)
                    hits.append((event, output, offset, len(data), left,
                                 len(data)-left-(end-start)))
                    continue
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
                        mismatch = _sequence_mismatch(contig, start, end, strand, qstart, qend)
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


def _sequence_mismatch(contig, start, end, strand, qstart, qend):
    return (f'sequence_mismatch: selected query {contig}:{start}-{end}{strand} '
            f'does not match INFO/SEQ at insertion offsets '
            f'{qstart}-{qend}; SIZE alone does not identify '
            'the representative allele')


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


def _needed_targets(events, reachable):
    """Names reachable variants are placed on that are not variants."""
    needed = set()
    for name in reachable:
        event = events[name]
        needed.add(event.chrom)
        needed.update(run.target for run in event.runs if run.target)
    return needed - set(events)


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


def _report_unresolved(path, failures, aliases, log, fatal=True):
    with path.open('w') as output:
        output.write('variant_id\tgfa_name\treason\n')
        for identifier, reason in failures:
            output.write(f'{identifier}\t{aliases.get(identifier, identifier)}\t{reason}\n')
    if failures and not fatal:
        categories = Counter(reason.split(':', 1)[0] for _identifier, reason in failures)
        log(f'Dropped {len(failures)} variant(s) whose sequence could not be verified '
            '(--drop-unverified); ' + ', '.join(f'{reason}={count}' for reason, count
                                               in categories.most_common()) + f'; see {path}')
        return
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
                     candidates, temporary_root, unresolved, aliases, log, dropped=None):
    """Retry failed candidates in scaffold batches, preserving successful reads.

    With ``dropped`` (a list), exhausted variants are collected there instead
    of stopping the run; the caller reports and removes them.
    """
    verified_intervals, hits = {}, {}
    attempts = Counter()
    round_number = 0
    if dropped is None:
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
        if dropped is None:
            _report_unresolved(unresolved, exhausted, aliases, log)
        else:
            dropped.extend(exhausted)
        if not retry:
            break
        log(f'Retrying {len(retry)} insertions with later SIZE/span-matched observations')
        intervals, missing = _prepare_intervals(events, retry, sources, args.anchor,
                                               query_flanks=args.gfa_mode != 'rgfa')
    return verified_intervals, hits


def _drop_unverified(events, resolver, intervals, dropped, unresolved, aliases, log):
    """Remove unverifiable variants and everything built on them from the graph.

    A dropped variant is never exported with unverified bases. Variants nested
    in it or aligned onto it are dropped too; the sidecar then marks all of
    them as not exported (path '.').
    """
    reasons = OrderedDict()
    for identifier, reason in dropped:
        reasons.setdefault(identifier, reason)
    order = getattr(resolver, 'order', None)
    if order is not None:
        for name in order:  # parents precede children
            if name in reasons:
                continue
            event = events[name]
            parent = next((target for target in (event.chrom, *(run.target for run in event.runs))
                           if target in reasons), None)
            if parent is not None:
                reasons[name] = f'depends_on_dropped: {aliases.get(parent, parent)}'
        resolver.order = [name for name in order if name not in reasons]
        resolver.variant_intervals = None
        resolver.source_intervals = None
        for name in reasons:
            resolver.ranks.pop(name, None)
    else:
        for name in reasons:
            events[name].retained = False
    for name in reasons:
        intervals.pop(name, None)
    _report_unresolved(unresolved, list(reasons.items()), aliases, log, fatal=False)


def _write_metadata(bed, events, intervals, aliases, resolver, unique_minimum, anchor,
                    mapping_prefix=None, sidecars=True):
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
                if mapping_prefix is not None:
                    out.write(mapping_prefix(event.identifier))
                out.write(f'{piece}\t{event_name}\t{run.qstart}\t{run.qend-run.qstart}\t'
                          f'{sample}\t{contig}\t{qstart}\t{qend}\t{strand}\t'
                          f'{aliases.get(ref_contig,ref_contig)}\t{ref_start}\t{ref_end}\t'
                          f'{run.operation}\t{target}\t{tstart}\t{tend}\t{tstrand}\t{ref_strand}\n')
    if sidecars and hasattr(resolver, 'source_aliases'):
        _write_source_sidecars(bed, resolver.source_aliases, resolver.roots)


def _write_source_sidecars(bed, source_aliases, roots):
    with open(str(bed)+'.local-paths.tsv', 'w') as out:
        out.write('local_path\ttarget\tstart\tend\tstrand\n')
        for name, alias in sorted(source_aliases.items()):
            out.write(f'{name}\t{alias.target}\t{alias.start}\t{alias.end}\t{alias.strand}\n')
    with open(str(bed)+'.template-lifts.tsv', 'w') as out:
        out.write('template\tbackbone\tstatus\tplacements\tnote\n')
        for name, root in sorted(roots.items()):
            if root.lift:
                lift = root.lift
                out.write('\t'.join(str(value) for value in (
                    name, lift.get('backbone', '.'), lift.get('status', '.'),
                    lift.get('placements', '.'), lift.get('note', '.'))) + '\n')


def _path_leaves(spec, roots, resolver, anchor):
    if spec.kind == 'junction':
        return resolver.junction_leaves(spec.start, max(1, anchor))
    if spec.kind in ('reference', 'alternative', 'novel', 'duplication'):
        return resolver.expand(spec.source, spec.start, spec.end)
    if hasattr(resolver, 'ranks') and spec.kind == 'deletion':
        return resolver.local_path(spec.source, spec.start, spec.end, max(1, anchor))
    if hasattr(resolver, 'ranks') and spec.kind in ('insertion', 'snp', 'substitution'):
        # Keep the named insertion path at its original 0-based core origin.
        # Parent flanks are used to form graph edges, not prepended to this P.
        return resolver.expand(spec.source, spec.start, spec.end)
    return resolver.local_path(spec.source, spec.start, spec.end, anchor)


def _link_leaves(spec, roots, resolver, anchor):
    if spec.kind == 'junction':
        return resolver.junction_leaves(spec.start, max(1, anchor))
    if hasattr(resolver, 'ranks') and spec.source in roots and roots[spec.source].lift:
        return resolver.local_path(spec.source, spec.start, spec.end, max(1, anchor))
    if hasattr(resolver, 'ranks') and spec.kind in ('insertion', 'snp', 'substitution', 'deletion'):
        # Even --anchor 0 must attach the allele to its parent graph.
        return resolver.local_path(spec.source, spec.start, spec.end, max(1, anchor))
    return _path_leaves(spec, roots, resolver, anchor)


def _path_and_link_leaves(spec, roots, resolver, anchor):
    """Resolve the core once and reuse it when attaching the parent flanks."""
    path = _path_leaves(spec, roots, resolver, anchor)
    if not hasattr(resolver, 'ranks') or spec.kind in ('junction', 'deletion'):
        return path, path
    attached = (spec.kind in ('insertion', 'snp', 'substitution') or
                (spec.source in roots and roots[spec.source].lift))
    if not attached:
        return path, path
    if spec.kind not in ('reference', 'alternative', 'novel', 'duplication',
                         'insertion', 'snp', 'substitution'):
        return path, _link_leaves(spec, roots, resolver, anchor)
    flank = max(1, anchor)
    link = (resolver._flank(spec.source, spec.start, flank, 'left') + path +
            resolver._flank(spec.source, spec.end, flank, 'right'))
    return path, link


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


def _write_segments(output, segment_ids, hits, roots, aliases, prefix, ranks=None,
                    root_sequences=None):
    query_files = _OpenFiles()
    root_reader = None
    root_path = None
    window = None
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
                    written = 0
                    record = root.record or source
                    low, high = root.base+start, root.base+end
                    if root_sequences is not None:
                        root_sequence = root_sequences[root.path, record]
                        pieces = (root_sequence[pos:min(high, pos+CHUNK)]
                                  for pos in range(low, high, CHUNK))
                    else:
                        if root.path != root_path:
                            if root_reader:
                                root_reader.close()
                            root_reader = IndexedFasta(root.path, root.fai)
                            root_path = root.path
                            window = None
                    if root_sequences is None and high - low <= ROOT_WINDOW:
                        # Segments of one root arrive in coordinate order: read
                        # the FASTA in large windows, not one call per segment.
                        if not (window and window[0] == record and
                                window[1] <= low and high <= window[1] + len(window[2])):
                            limit = root_reader.index[record][0] if record in root_reader.index else high
                            window = (record, low, ''.join(interval_chunks(
                                root_reader, record, low, max(high, min(limit, low + ROOT_WINDOW)), '+')))
                        pieces = (window[2][low-window[1]:high-window[1]],)
                    elif root_sequences is None:
                        pieces = interval_chunks(root_reader, record, low, high, '+')
                    for sequence in pieces:
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


_WRITE_STATE = None


def _preload_root_sequences(segment_ids, roots, log):
    """Load only FASTA records supplying S lines, before the writer fork.

    Duplication slices share their backing record. Lookup-only catalog names
    and query assemblies are excluded. Immutable strings share their large
    data buffers through copy-on-write, even across many writer processes.
    """
    import numpy as np
    from gfa_topology import POINT_BITS

    sequences = {}
    with ExitStack() as stack:
        readers, wanted = {}, {}
        for source in np.flatnonzero(~segment_ids.is_query).tolist():
            root = roots.get(segment_ids.names[source])
            if root is None:
                continue
            low, high = np.searchsorted(segment_ids.keys,
                                        [source << POINT_BITS, (source+1) << POINT_BITS])
            if not segment_ids.used[low:high].any():
                continue
            record = root.record or root.name
            key = root.path, record
            if key in wanted:
                continue
            if root.path not in readers:
                readers[root.path] = stack.enter_context(IndexedFasta(root.path, root.fai))
            wanted[key] = readers[root.path].index[record][0]
        if wanted:
            log(f'Preloading {len(wanted)} reference/alternative FASTA records '
                f'({sum(wanted.values()) / (1024**3):.2f} GiB of bases) for shared GFA writers')
        started = time.monotonic()
        for (path, record), length in wanted.items():
            sequence = readers[path].sequence(record)
            if len(sequence) != length:
                raise ValueError(f'{record}: incomplete FASTA sequence')
            sequences[path, record] = sequence
        if wanted:
            log(f'Preloaded reference/alternative sequences in {time.monotonic() - started:.1f}s')
    return sequences


def _write_part(task):
    """Write one ordered block of S or P lines to its own temporary file."""
    kind, low, high, path = task
    specs, roots, leaves, first, last, numbers, segments, prefix, hits, aliases, ranks, root_sequences = _WRITE_STATE
    import gfa_topology as topology
    with open(path, 'wb', buffering=CHUNK) as output:
        if kind == 'S':
            _write_segments(output, segments.part(low, high), hits, roots, aliases, prefix,
                            ranks, root_sequences)
            return path
        for index in range(low, high):
            spec = specs[index]
            if spec.kind == 'junction':
                continue  # link-only walk
            text = topology.path_text(index, leaves, first, last, numbers, prefix)
            if not text:
                raise ValueError(f'{spec.name}: empty GFA path')
            lift = roots[spec.source].lift if spec.source in roots else None
            lift_tag = f'\tLS:Z:{lift["status"]}' if lift else ''
            if lift and lift['status'] not in ('mapped', 'one_sided'):
                lift_tag += '\tUP:Z:unplaced'
            output.write(f'P\t{spec.name}\t{text}\t*\tTP:Z:{spec.kind}{lift_tag}\n'.encode('ascii'))
    return path


def _write_parts(temporary_root, specs, roots, leaves, first, last, numbers, segment_ids,
                 prefix, hits, aliases, ranks, stable, processes=1, tag='',
                 preload_reference=True, log=print):
    """Write ordered S and P blocks (forked workers) and compute the links.

    Returns (S part files, P part files, (left, right, rank) link arrays).
    """
    global _WRITE_STATE
    import gc
    import multiprocessing as mp
    import gfa_topology as topology
    blocks = max(1, processes) * 4
    tasks = []
    for kind, total in (('S', len(segment_ids)), ('P', len(specs))):
        size = max(1, -(-total // blocks))
        for number, low in enumerate(range(0, total, size)):
            tasks.append((kind, low, min(total, low + size),
                          temporary_root / f'part.{tag}{kind}.{number:06d}.gfa'))
    root_sequences = _preload_root_sequences(segment_ids, roots, log) if preload_reference else None
    _WRITE_STATE = (specs, roots, leaves, first, last, numbers, segment_ids, prefix,
                    hits, aliases, ranks, root_sequences)
    try:
        if processes > 1 and len(tasks) > 1 and 'fork' in mp.get_all_start_methods():
            gc.freeze()
            try:
                with mp.get_context('fork').Pool(processes, maxtasksperchild=1) as pool:
                    pending = pool.map_async(_write_part, tasks, 1)
                    # Links are computed here while the workers write S and P.
                    links = topology.links(leaves, first, last, numbers, stable)
                    pending.get()
            finally:
                gc.unfreeze()
        else:
            for task in tasks:
                _write_part(task)
            links = topology.links(leaves, first, last, numbers, stable)
    finally:
        _WRITE_STATE = None
    return ([task[3] for task in tasks if task[0] == 'S'],
            [task[3] for task in tasks if task[0] == 'P'], links)


def _copy_parts(parts, output, remove=True):
    for path in parts:
        with open(path, 'rb') as part:
            shutil.copyfileobj(part, output, length=CHUNK)
        if remove:
            os.unlink(path)


def _write_links(output, left, right, link_ranks, prefix, ranked):
    step = 100_000
    for low in range(0, len(left), step):
        part = slice(low, low + step)
        rows = zip((left[part] >> 1).tolist(), (left[part] & 1).tolist(),
                   (right[part] >> 1).tolist(), (right[part] & 1).tolist(),
                   link_ranks[part].tolist())
        output.write(''.join(
            f'L\t{prefix}{a}\t{"+-"[ao]}\t{prefix}{b}\t{"+-"[bo]}\t0M'
            + (f'\tSR:i:{rank}\n' if ranked else '\n')
            for a, ao, b, bo, rank in rows).encode('ascii'))


def _write_graph(output_path, temporary_root, specs, roots, resolver, leaves,
                 first, last, numbers, segment_ids, prefix, hits, aliases, stable,
                 processes=1, preload_reference=True, log=print):
    """Write H, S, L, P lines; S and P blocks are written by forked workers."""
    ranks = getattr(resolver, 'ranks', None)
    s_parts, p_parts, (left, right, link_ranks) = _write_parts(
        temporary_root, specs, roots, leaves, first, last, numbers, segment_ids, prefix,
        hits, aliases, ranks, stable, processes, preload_reference=preload_reference, log=log)
    temporary = str(output_path) + f'.tmp.{os.getpid()}'
    try:
        opener = gzip.open if str(output_path).endswith('.gz') else open
        with opener(temporary, 'wb') as output:
            output.write(b'H\tVN:Z:1.0\tTS:Z:merged_vcf_to_gfa.py\n')
            _copy_parts(s_parts, output)
            _write_links(output, left, right, link_ranks, prefix, ranks is not None)
            _copy_parts(p_parts, output)
        os.replace(temporary, output_path)
    finally:
        if os.path.exists(temporary):
            os.unlink(temporary)
    return len(left)


def run(args, _read_vcf=None, _chunks=None, _operations=None, log=print):
    from gfa_interval_metadata import vcf_paths
    bed = Path(args.anchor_bed or (str(args.output or next(vcf_paths(args.vcf))) + '.anchors.bed')).resolve()
    bed.parent.mkdir(parents=True, exist_ok=True)
    if getattr(args, 'partition_records', 0):
        from gfa_partitioned import run_partitioned
        return run_partitioned(args, bed, log)
    with tempfile.TemporaryFile(dir=bed.parent) as spool:
        return _run(args, QueryCandidates(spool), bed, log)


def _run(args, candidates, bed, log):
    from gfa_interval_metadata import read_header_contigs, rows

    log('Indexing query list, VCF headers, and reference FASTAs')
    sources = read_query_list(args.query_fasta_list)
    header_ordinals, header_lengths = read_header_contigs(args.vcf)
    roots = _index_roots(args)
    sequence_checks = {} if args.gfa_mode == 'rgfa' else None
    log('Streaming VCF metadata and SIZE/span-matched query intervals')
    started = time.monotonic()
    events = OrderedDict()
    event_kinds = {}
    record_indexes = {}
    dup_alignments = {} if args.gfa_mode == 'rgfa' else None
    repeated = set()
    for order, record in enumerate(rows(args.vcf, sources, args.size_cutoff,
                                       sequence_checks=sequence_checks,
                                       candidate_sink=candidates.add if sequence_checks is not None else None,
                                       insertion_only=getattr(args, 'insertion_only', None) is not None,
                                       event_kinds=event_kinds, record_indexes=record_indexes,
                                       dup_alignments=dup_alignments, repeated=repeated)):
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
    log(f'Indexed {len(events)} variant definitions in {time.monotonic() - started:.1f}s')

    aliases = {name: f'ins_{header_ordinals[name]}'
               for name in events if name in header_ordinals}
    stable = args.gfa_mode == 'rgfa'
    reachable = _reachable_events(events, include_parents=stable)
    hits = {}
    if stable:
        log('Resolving local source catalogs and variant dependencies')
        started = time.monotonic()
        from gfa_source_catalog import resolve_local_sources
        from gfa_stable_coords import StableResolver
        source_aliases = resolve_local_sources(events, roots, reachable, sources,
                                              list(_flatten(args.local_path_fasta)), log,
                                              template_catalogs=args.local_reference_templates,
                                              backbone=getattr(args, 'reference_haplotype', None))
        from gfa_duplications import add_template_roots
        from gfa_interval_metadata import TEMPLATE_PATH
        needed = _needed_targets(events, reachable)
        add_template_roots(needed, roots, source_aliases, Root, strict=False)
        catalog = getattr(args, 'alternative_catalog', None)
        if catalog:
            from gfa_catalog_paths import resolve_catalog_paths
            from gfa_interval_metadata import vcf_paths
            # DUP_ templates are not catalog records; their loci may be.
            source_aliases.update(resolve_catalog_paths(
                events, roots, reachable, catalog, list(vcf_paths(args.vcf)),
                getattr(args, 'reference_haplotype', None), args.reference_fasta, args.output or bed,
                Root, log, threads=args.processes,
                needed={name for name in needed if not TEMPLATE_PATH.fullmatch(name)}))
        add_template_roots(needed, roots, source_aliases, Root)
        from gfa_duplications import add_duplication_paths
        add_duplication_paths(events, roots, dup_alignments, source_aliases, Run, Root, log)
        resolver = StableResolver(events, roots, header_lengths, hits, reachable,
                                  args.nested_pos_base, source_aliases=source_aliases)
        coordinates = resolver.index_variants()
        log(f'Indexed variant coordinates on {len(coordinates.by_parent)} chromosomes/insertion parents')
        del coordinates
        stable_names = [aliases.get(name, name) for name in resolver.order] + list(roots)
        if len(stable_names) != len(set(stable_names)):
            raise ValueError('insertion aliases collide with another insertion or FASTA root')
        for name in stable_names:
            _valid_name(name)
        stable_names = set(stable_names)
        log(f'Resolved source catalogs and dependencies in {time.monotonic() - started:.1f}s')
    else:
        resolver = Resolver(events, roots, header_lengths, hits)
    log('Preparing query intervals')
    intervals, missing = _prepare_intervals(events, reachable, sources, args.anchor,
                                          query_flanks=not stable)
    unresolved = Path(str(bed)+'.unresolved.tsv')
    unique_minimum = max(args.size_cutoff, args.unique_size_cutoff)
    # Reject known-unrecoverable metadata before reading any assembly sequence.
    drop = bool(getattr(args, 'drop_unverified', False))
    unrecoverable = [item for item in missing
                     if not stable or item[0] not in candidates.positions]
    dropped = [] if drop else None
    if drop:
        # Collected here and reported once, after extraction.
        dropped.extend(unrecoverable)
    else:
        _report_unresolved(unresolved, unrecoverable, aliases, log)
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
            candidates if stable else None, temporary_root, unresolved, aliases, log,
            dropped=dropped)
        hits.update(extracted)
        if drop:
            _drop_unverified(events, resolver, intervals, dropped, unresolved, aliases, log)
        log('Writing anchor and mapping metadata')
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
        if stable:
            from gfa_junctions import junction_specs
            log('Resolving carrier junctions')
            specs.extend(junction_specs(args.vcf, events, resolver, aliases, PathSpec, log,
                                        repeated=repeated))

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
        import gfa_topology as topology
        started = time.monotonic()
        ids, source_names, is_query = topology.source_order(
            roots, events, resolver, header_lengths, stable)
        leaves = topology.collect_leaves(specs, roots, resolver, args.anchor, ids,
                                         _path_leaves, _link_leaves, args.processes, log,
                                         paired_leaves=_path_and_link_leaves)
        del ids
        log(f'Topology: {len(leaves.source)} leaves for {len(specs)} paths in '
            f'{time.monotonic() - started:.1f}s')
        started = time.monotonic()
        keys = topology.boundaries(leaves, stable, args.max_node_length)
        first, last = topology.pair_ranges(keys, leaves.source, leaves.start, leaves.end)
        forbidden = set()
        if stable:
            prefix = ''  # Preserve numeric node IDs through VG import/export.
            forbidden = {int(name) for name in names | stable_names
                         if name.isascii() and name.isdigit() and name == str(int(name))}
        used, numbers = topology.number_segments(keys, leaves, first, last, stable, forbidden)
        topology.check_path_coverage(used, first, last, leaves, specs)
        segment_ids = topology.Segments(keys, used, numbers, source_names, is_query)
        log(f'Topology: numbered {len(segment_ids)} segments from {len(keys)} cut points in '
            f'{time.monotonic() - started:.1f}s')
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
        started = time.monotonic()
        links = _write_graph(output, temporary_root, specs, roots, resolver, leaves,
                             first, last, numbers, segment_ids, prefix, hits, aliases, stable,
                             args.processes, getattr(args, 'preload_reference', True), log)
        log(f'Wrote {len(segment_ids)} segments, {links} links, and '
            f'{len(specs)} paths to {output} in {time.monotonic() - started:.1f}s')
        if stable:
            index = Path(str(output) + '.variants.tsv')
            count = _write_variant_index(index, events, record_indexes, aliases, resolver)
            log(f'Wrote {count} variant-to-path rows to {index}')
        return 0


def _write_variant_index(path, events, record_indexes, aliases, resolver):
    from gfa_interval_metadata import TEMPLATE_PATH
    """Map each VCF data line to its GFA allele path and parent breakpoints.

    ``vcf_index`` numbers VCF data lines from 0 across the -v inputs in order.
    ``path`` is the allele's P line ('.' when not exported). Parent intervals
    are in the parent P path's own coordinates; ``parent_strand`` is '-' when
    a local parent name maps in reverse onto its included path. Deletion paths
    carry flanks, so a GAF walk uses the parent interval, not their P line.
    """
    rows = []
    for event in events.values():
        exported = event.identifier in resolver.ranks
        parent = start = end = strand = '.'
        allele_path = aliases.get(event.identifier, event.identifier) if exported else '.'
        if exported:
            start, end = resolver.breakpoints(event.identifier)
            parent, strand = event.chrom, '+'
            if parent in events and events[parent].kind == 'deletion':
                # Rows nested in a deletion keep their offset (allele order).
                start = end = event.pos - resolver.nested_pos_base
            if parent not in events:
                parent, start, end, strand = resolver.canonical_interval(parent, start, end)
            parent = aliases.get(parent, parent)
            runs = event.runs
            if (len(runs) == 1 and runs[0].operation == '=' and runs[0].target in resolver.roots
                    and TEMPLATE_PATH.fullmatch(runs[0].target)
                    and resolver.roots[runs[0].target].kind == 'duplication'
                    and runs[0].rend - runs[0].rstart == resolver.roots[runs[0].target].length):
                # A duplication site's allele is its template path, which
                # carries the copy's nested rows.
                allele_path = runs[0].target
        rows.append((record_indexes.get(event.identifier, -1), event.identifier, allele_path,
                     event.kind, parent, start, end, strand))
    rows.sort(key=lambda row: row[0])
    with open(path, 'w', buffering=CHUNK) as out:
        out.write('#vcf_index\tvariant_id\tpath\tkind\tparent\tparent_start\t'
                  'parent_end\tparent_strand\n')
        for row in rows:
            out.write('\t'.join(map(str, row)) + '\n')
    return len(rows)
