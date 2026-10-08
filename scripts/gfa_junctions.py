"""Link walks for variants that follow each other in a carrier's allele.

Every variant path is attached to its parent by flank walks (a base on each
side), which give the links for entering and leaving one variant. When a
carrier's allele has two variants of one parent with no parent base between
them (adjacent SNPs, an insertion next to a deletion, the kept-base
insertions nested in a --exact deletion, ...), its walk goes from one
variant's allele straight into the next one's. Those links come from here.

Rows that abut a sibling on their parent (one's end is the other's start),
insertions with two alleles of one sample, and everything nested in them are read again from the VCFs for their
carriers' QUERYCOORD. A carrier's siblings that abut on the parent and touch
on the query form a chain; its link walk is the parent's flank, each
variant's allele as that carrier walks it (its own nested rows, recursively,
in query order), then the parent's flank. Identical walks are kept once and
written as link-only walks (no P line). A --exact duplication template
(DUP_ path) has no flank into its sites' parents, so each carrier's site
allele with nested rows gets a walk too.
"""
from array import array
from collections import defaultdict
from itertools import groupby
from operator import itemgetter
import re
import time

import numpy as np

from allele_edges import edge_owned
from gfa_interval_metadata import TEMPLATE_PATH, _file_lines, _selected_lines, vcf_paths

_COORD = re.compile(r'(\d+)(?:-(\d+))?([+-]?)')


def _allele_parent(event, roots):
    """The name the rows nested in this allele are placed on."""
    runs = event.runs
    if (len(runs) == 1 and runs[0].operation == '=' and runs[0].target in roots
            and TEMPLATE_PATH.fullmatch(runs[0].target)
            and roots[runs[0].target].kind == 'duplication'):
        return runs[0].target
    return event.identifier


def _involved(events, resolver, repeated=()):
    """Exported rows in a junction, and their breakpoints on the parent.

    ``repeated``: insertions with two alleles of one sample, which can
    follow themselves at their breakpoint."""
    coordinates = resolver.index_variants()
    spans, children = coordinates.spans, coordinates.by_parent
    seeds = set()
    for rows in children.values():
        # Merge sorted boundaries, rather than allocating two dictionaries of
        # lists keyed by every chromosome/position in a SNP cohort.
        starts = groupby(rows, key=itemgetter(0))
        ends = groupby(sorted(rows, key=itemgetter(1)), key=itemgetter(1))
        first, last = next(starts, None), next(ends, None)
        while first is not None and last is not None:
            if first[0] < last[0]:
                first = next(starts, None)
            elif last[0] < first[0]:
                last = next(ends, None)
            else:
                starters = [row[2] for row in first[1]]
                enders = [row[2] for row in last[1]]
                # A lone zero-length insertion only abuts itself.
                if len(enders) > 1 or len(starters) > 1 or enders[0] != starters[0]:
                    seeds.update(enders)
                    seeds.update(starters)
                first, last = next(starts, None), next(ends, None)
    for name, (start, end) in spans.items():
        if name in repeated and start == end:
            seeds.add(name)
        parent = _allele_parent(events[name], resolver.roots)
        if parent != name and children.get(parent):
            seeds.add(name)
    involved, pending = set(), list(seeds)
    while pending:
        name = pending.pop()
        if name not in involved:
            involved.add(name)
            pending.extend(row[2] for row in children.get(
                _allele_parent(events[name], resolver.roots), ()))
    return involved, spans


def _carriers(vcfs, wanted, events, select=None, index=None):
    """Carriers' observations of the wanted rows, parsed as
    gfa_sample_gaf.build_shards does, as compact columns: (row names,
    {sample, contig, row, low, high, reverse} numpy arrays)."""
    rows, row_ids, sample_ids, contig_ids = [], {}, {}, {}
    columns = {name: array(code) for name, code in (
        ('sample', 'i'), ('contig', 'i'), ('row', 'i'), ('low', 'q'), ('high', 'q'),
        ('reverse', 'b'))}
    for number, path in enumerate(vcf_paths(vcfs)):
        lines = (_selected_lines(path, select, index, number) if select is not None
                 else _file_lines(path, [0]))
        samples = ()
        for _label, line, data in lines:
            if line.startswith('#CHROM\t'):
                samples = [sample_ids.setdefault(name, len(sample_ids))
                           for name in line.rstrip('\r\n').split('\t')[9:]]
                continue
            if data is None:
                continue
            head = line.split('\t', 3)
            if len(head) < 4 or head[2] not in wanted:
                continue
            fields = line.rstrip('\r\n').split('\t')
            if len(fields) < 10:
                continue
            keys = fields[8].split(':')
            try:
                gt, ctg, qc = keys.index('GT'), keys.index('ASSEMBLYCONTIG'), keys.index('QUERYCOORD')
            except ValueError:
                continue
            snp = getattr(events[fields[2]], 'kind', 'insertion') == 'snp'
            row = None
            for sample, text in zip(samples, fields[9:]):
                if not text.startswith('1'):
                    continue
                values = text.split(':')
                if len(values) != len(keys) or values[gt] != '1':
                    continue
                names, coords = values[ctg].split(','), values[qc].split(',')
                for position, coord in enumerate(coords):
                    match = _COORD.fullmatch(coord)
                    if match is None:
                        continue
                    low = int(match[1])
                    if snp and not match[2] and match[3] == '-':
                        low -= 1  # '-' points are traversal boundaries
                    high = int(match[2]) if match[2] else low + (1 if snp else 0)
                    if row is None:
                        row = row_ids.setdefault(fields[2], len(row_ids))
                        if row == len(rows):
                            rows.append(fields[2])
                    contig = names[position if len(names) > 1 else 0]
                    for name, value in (('sample', sample),
                                        ('contig', contig_ids.setdefault(contig, len(contig_ids))),
                                        ('row', row), ('low', min(low, high)),
                                        ('high', max(low, high)), ('reverse', match[3] == '-')):
                        columns[name].append(value)
    arrays = {name: np.frombuffer(values, dtype=np.dtype(values.typecode)).astype(
        np.int64 if name in ('low', 'high') else np.int32) if values else
        np.zeros(0, np.int64) for name, values in columns.items()}
    return rows, arrays


def junction_specs(vcfs, events, resolver, aliases, spec_type, log, select=None, index=None,
                   repeated=()):
    """Set resolver.junctions and return their link-only path specs."""
    started = time.monotonic()
    involved, spans = _involved(events, resolver, repeated)
    resolver.junctions = []
    if not involved:
        return []
    rows, columns = _carriers(vcfs, involved, events, select, index)
    # Observations sorted by carrier contig, then query; blocks per contig.
    order = np.lexsort((columns['high'], columns['low'], columns['contig'], columns['sample']))
    sample, contig = columns['sample'][order], columns['contig'][order]
    lows, highs = columns['low'][order], columns['high'][order]
    names = [rows[row] for row in columns['row'][order].tolist()]
    strands = ['-' if value else '+' for value in columns['reverse'][order].tolist()]
    edges = np.flatnonzero((sample[1:] != sample[:-1]) | (contig[1:] != contig[:-1])) + 1
    blocks = list(zip([0, *edges.tolist()], [*edges.tolist(), len(order)]))
    lows_list, highs_list = lows.tolist(), highs.tolist()
    roots = resolver.roots

    def nested(row, block, low, high, strand):
        """The carrier's rows nested in this allele, in its walk order."""
        parent = _allele_parent(events[row], roots)
        deletion = getattr(events[row], 'kind', 'insertion') == 'deletion'
        first, last = block
        found, starts, ends, kept_spans = [], [], [], []
        current = events[row].length
        for position in range(first + int(np.searchsorted(lows[first:last], low)), last):
            a, b, child = lows_list[position], highs_list[position], names[position]
            if a > high:
                break
            if b <= high and child != row and events[child].chrom == parent:
                item = (a, b, child, strands[position])
                if a == b and low < high and a in (low, high) and not deletion:
                    edge = starts if a == (high if strand == '-' else low) else ends
                    edge.append((*spans[child], item))
                else:
                    found.append(item)
                    kept_spans.append(spans[child])
                    current += (b - a) - (spans[child][1] - spans[child][0])
        if not deletion:
            # Edge points shared with an abutting allele (allele_edges.py).
            found.extend(edge_owned(starts, ends, kept_spans, current, high - low,
                                    events[row].length))

        def along(item):
            walk = item[0] if strand == '+' else -item[1]
            # A deletion's nested rows keep their offsets (all at breakpoint 0).
            return ((events[item[2]].pos, walk) if deletion else
                    (*spans[item[2]], walk))
        found.sort(key=along)
        return tuple((child, nested(child, block, a, b, child_strand))
                     for a, b, child, child_strand in found)

    walks = set()
    for block in blocks:
        groups = defaultdict(list)
        for position in range(*block):
            groups[events[names[position]].chrom, strands[position]].append(
                (lows_list[position], highs_list[position], names[position]))
        for (parent, strand), members in groups.items():
            # Walk order along the parent: query order, and parent order for
            # rows at one query point (e.g. abutting deletions).
            members.sort(key=(lambda m: (m[0], m[1], *spans[m[2]])) if strand == '+' else
                         (lambda m: (-m[1], -m[0], *spans[m[2]])))
            chain = []
            for member in members + [None]:
                if chain and member is not None:
                    low0, high0, previous = chain[-1]
                    low, high, row = member
                    touching = high0 == low if strand == '+' else low0 == high
                    if touching and spans[previous][1] == spans[row][0]:
                        chain.append(member)
                        continue
                if chain:
                    items = tuple((row, nested(row, block, low, high, strand))
                                  for low, high, row in chain)
                    single = items[0]
                    if len(items) > 1 or (
                            _allele_parent(events[single[0]], roots) != single[0] and single[1]):
                        walks.add((parent, spans[chain[0][2]][0], spans[chain[-1][2]][1], items))
                chain = [member] if member is not None else []

    def exported(item):
        return item[0] in resolver.ranks and all(exported(child) for child in item[1])

    resolver.junctions = sorted(walk for walk in walks
                                if all(exported(item) for item in walk[3]))
    specs = []
    for number, (_parent, _start, _end, items) in enumerate(resolver.junctions):
        first = items[0][0]
        specs.append(spec_type(f'{aliases.get(first, first)}#junction{number + 1}',
                               'junction', first, number, number))
    log(f'Junctions: {len(specs)} link walk(s) from {len(rows)} abutting row(s) in '
        f'{time.monotonic() - started:.1f}s')
    return specs
