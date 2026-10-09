#!/usr/bin/env python3
"""Check a GrvcfGraph segment GFA against the Python resolver and the inputs.

0. Paths: GrvcfGraph's --roots-table equals the Python resolution
   (_index_roots, resolve_local_sources, add_template_roots): kind, emitted,
   output order, FASTA record/offset/length, alias and lift of every name.
1. Placements: every row GrvcfGraph kept is recomputed with
   gfa_stable_coords.StableResolver (the Python converter's coordinates) and
   mapped to the sequence that carries its bases (alias -> included path,
   duplication site -> DUP_ template, row nested in a deletion -> that
   deletion's interval). Both tools must keep the same rows and agree on
   sequence, start, end and strand.
2. Segments: S lines of one sequence (SN) tile it from offset 0 without gaps,
   none is longer than --max-node-length, IDs run 1..N in file order, and the
   bases equal the FASTA record (paths) or the VCF allele (rows).
3. Cuts: every placement's start and end on its sequence is a segment boundary.
4. Links: the L lines are exactly the links of the Python converter's link
   walks (gfa_interval_pipeline._link_leaves with a 1-base flank, i.e. each
   row's allele entered from the base before it and left to the base after
   it, through parents, deletions, duplication sites and template lifts; each
   emitted path), mapped onto these segments. Carrier junctions are not
   expected (they come from the GAF walks).
5. Paths: every P line walks exactly the segments of the Python path walk.

    python tools/check_grvcf_nodes.py -g cohort.gfa -p placements.tsv \
        --roots-table roots.tsv -v VCF... -r REF.fa [-a ALT.fa] [--local-reference-templates T.fa] \
        [--local-path-fasta CATALOG] [--reference-haplotype CHM13_h1]
"""
import argparse
from bisect import bisect_right
from collections import defaultdict
from pathlib import Path
import re
import sys

# Repository layout (scripts/) or the flat cluster copy (tools/ next to the scripts).
for _folder in (Path(__file__).resolve().parents[1] / 'scripts', Path(__file__).resolve().parents[1]):
    if (_folder / 'gfa_stable_coords.py').is_file():
        sys.path.insert(0, str(_folder))
        break

_CHUNK = re.compile(r'[<>]([^<>]*)')
_OP = re.compile(r'(\d+)([=MXIDHS])([A-Za-z]*)')


def handle_of(segment, orientation):
    return 2 * int(segment) + (orientation == '-')


def canonical(a, b):
    """A link and its reverse complement are one link."""
    return min((a, b), (b ^ 1, a ^ 1))


def literal(text, size):
    """INFO/SEQ bases (plain, or graph-encoded with literal payloads only)."""
    from gfa_query_anchors import _unescape
    text = _unescape(text)
    if text[:1] not in ('<', '>'):
        return text
    bases = []
    for chunk in _CHUNK.finditer(text):
        body = chunk[1].split(':', 1)[-1]
        for length, operation, payload in _OP.findall(body):
            if operation in 'IXS' and int(length):
                bases.append(payload)
    sequence = ''.join(bases)
    if len(sequence) != size:
        raise ValueError('decoded SEQ length differs from SVLEN')
    return sequence


def row_alleles(paths, wanted):
    """ID -> allele bases for rows with their own sequence, from columns 1-8."""
    from gfa_interval_metadata import _file_lines, vcf_paths
    alleles = {}
    for path in vcf_paths(paths):
        for _lineno, line, data in _file_lines(path, [0]):
            if data is None:
                continue
            fields = line.rstrip('\r\n').split('\t', 8)
            if fields[2] not in wanted:
                continue
            info = dict(item.split('=', 1) for item in fields[7].split(';') if '=' in item)
            svtype = info.get('SVTYPE', '')
            if svtype == 'INS':
                alleles[fields[2]] = literal(info['SEQ'], abs(int(info['SVLEN'])))
            elif svtype == 'SUB':
                from gfa_query_anchors import _unescape
                alleles[fields[2]] = _unescape(info['SEQ'])
            else:
                alleles[fields[2]] = fields[4]
    return alleles


def resolve_roots(targets, reference, alternatives=(), templates=(), catalogs=(), backbone=None):
    """The Python converter's FASTA paths and aliases for these target names."""
    from types import SimpleNamespace
    from gfa_duplications import add_template_roots
    from gfa_interval_pipeline import Root, _index_roots
    from gfa_source_catalog import resolve_local_sources

    args = SimpleNamespace(gfa_mode='rgfa', reference_fasta=reference, reference_fai=None,
                           alternatives_fasta=[list(alternatives)],
                           reference_haplotype=backbone, size_cutoff=0)
    roots = _index_roots(args)
    aliases = resolve_local_sources(None, roots, None, {}, list(catalogs), lambda message: None,
                                    template_catalogs=list(templates),
                                    backbone=args.reference_haplotype, needed=set(targets))
    add_template_roots(set(targets), roots, aliases, Root, strict=True)
    return roots, aliases


def roots_table(roots, aliases):
    """name -> comparable fields, as GrvcfGraph --roots-table writes them."""
    table = {}
    for root in roots.values():
        alias = aliases.get(root.name)
        lift = root.lift or {}

        def side(name):
            value = lift.get(name)
            return f'{value[0]}:{value[1]}:{value[2]}' if value else '.'
        table[root.name] = (root.kind, '1' if root.emit else '0', str(root.order),
                            root.record or root.name, str(root.base), str(root.length),
                            *((alias.target, str(alias.start), str(alias.end), alias.strand) if alias
                              else ('.', '.', '.', '.')),
                            lift.get('status', '.'), side('left'), side('right'))
    return table


def python_placements(args):
    """ID -> (sequence, start, end, strand, inside_deletion) from StableResolver."""
    from gfa_interval_metadata import read_header_contigs, rows, TEMPLATE_PATH
    from gfa_interval_pipeline import Event, Hit, Run, _needed_targets, _reachable_events
    from gfa_stable_coords import StableResolver

    _ordinals, header_lengths = read_header_contigs(args.vcf)
    events, kinds = {}, {}
    for order, record in enumerate(rows(args.vcf, {}, 0, sequence_checks={}, event_kinds=kinds)):
        identifier, chrom, pos, ref_end, length, values, query, reason, retained = record
        events[identifier] = Event(order, identifier, chrom, pos, ref_end, length,
                                   [Run(*run) for run in values], query, reason, retained,
                                   kinds.get(identifier, 'insertion'))
    reachable = _reachable_events(events, include_parents=True)
    needed = _needed_targets(events, reachable)
    roots, aliases = resolve_roots(needed, args.reference_fasta, args.alternatives_fasta,
                                   args.local_reference_templates, args.local_path_fasta,
                                   args.reference_haplotype)

    def site_template(name):
        runs = events[name].runs
        if (len(runs) == 1 and runs[0].operation == '=' and runs[0].target in roots
                and TEMPLATE_PATH.fullmatch(runs[0].target)):
            return runs[0].target
        return None

    # Literal alleles as 'query' leaves in their own coordinates (no flanks).
    hits = {name: Hit('', 0, events[name].length, 0, 0) for name in reachable
            if events[name].kind != 'deletion' and not site_template(name)}
    resolver = StableResolver(events, roots, header_lengths, hits, reachable, 1,
                              source_aliases=aliases)

    placements = {}

    def place(name):
        if name in placements:
            return placements[name]
        event = events[name]
        parent = events.get(event.chrom)
        if parent is not None and parent.kind == 'deletion':
            value = place(event.chrom)[:4] + (1,)
        else:
            start, end = resolver.breakpoints(name)
            if parent is not None:
                value = (site_template(event.chrom) or event.chrom, start, end, '+', 0)
            elif event.chrom in aliases:
                value = aliases[event.chrom].interval(start, end) + (0,)
            else:
                value = (event.chrom, start, end, '+', 0)
        placements[name] = value
        return value

    for name in reachable:
        place(name)
    sequences = {}
    for root in roots.values():
        if root.emit and root.name not in aliases and root.length:
            sequences[root.name] = ('root', root)
    for name in reachable:
        event = events[name]
        if event.kind != 'deletion' and not site_template(name) and event.length:
            sequences[name] = ('row', event.length)
    # Link walks and path walks, as gfa_interval_pipeline builds them.
    link_walks, path_walks = [], {}
    for name in reachable:
        link_walks.append(resolver.local_path(name, 0, events[name].length, 1))
    for root in roots.values():
        if root.emit and root.length:
            path_walks[root.name] = resolver.expand(root.name, 0, root.length)
            link_walks.append(resolver.local_path(root.name, 0, root.length, 1) if root.lift
                              else path_walks[root.name])
    return placements, sequences, roots_table(roots, aliases), link_walks, path_walks


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__.split('\n\n')[0])
    parser.add_argument('-g', '--gfa', required=True)
    parser.add_argument('-p', '--placements', required=True, help='GrvcfGraph --placements TSV')
    parser.add_argument('--roots-table', help='GrvcfGraph --roots-table TSV (compared when given)')
    parser.add_argument('-v', '--vcf', nargs='+', required=True, action='extend')
    parser.add_argument('-r', '--reference-fasta', required=True)
    parser.add_argument('-a', '--alternatives-fasta', nargs='+', action='extend', default=[])
    parser.add_argument('--local-reference-templates', nargs='+', action='extend', default=[])
    parser.add_argument('--local-path-fasta', nargs='+', action='extend', default=[])
    parser.add_argument('--reference-haplotype')
    parser.add_argument('--max-node-length', type=int, default=1024)
    parser.add_argument('--show', type=int, default=10, help='problems printed per check')
    args = parser.parse_args(argv)

    problems = defaultdict(list)

    def problem(kind, text):
        problems[kind].append(text)

    expected, sequences, python_roots, link_walks, path_walks = python_placements(args)
    if args.roots_table:
        observed_roots = {}
        with open(args.roots_table) as handle:
            for line in handle:
                if line.startswith('#'):
                    continue
                fields = line.rstrip('\n').split('\t')
                # Skip the FASTA path column: the Python side keeps it as given.
                observed_roots[fields[0]] = (fields[1], fields[2], fields[3], *fields[5:])
        for name in sorted(set(python_roots) | set(observed_roots)):
            if python_roots.get(name) != observed_roots.get(name):
                problem('roots', f'{name}: python {python_roots.get(name)} grvcfgraph {observed_roots.get(name)}')
        print(f'paths: {len(observed_roots)} GrvcfGraph, {len(python_roots)} Python, '
              f'{len(problems.get("roots", []))} differ')
    observed = {}
    with open(args.placements) as handle:
        next(handle)
        for line in handle:
            fields = line.rstrip('\n').split('\t')
            observed[fields[0]] = (fields[5], int(fields[6]), int(fields[7]), fields[8], int(fields[9]))
    for name in sorted(set(expected) - set(observed)):
        problem('kept_only_by_python', name)
    for name in sorted(set(observed) - set(expected)):
        problem('kept_only_by_grvcfgraph', name)
    for name in sorted(set(expected) & set(observed)):
        if expected[name] != observed[name]:
            problem('placement', f'{name}: python {expected[name]} grvcfgraph {observed[name]}')
    print(f'placements: {len(observed)} GrvcfGraph rows, {len(expected)} Python rows, '
          f'{len(problems.get("placement", []))} differ')

    from minsetref_core import IndexedFasta
    readers = {}
    alleles = row_alleles(args.vcf, {name for name, (kind, _value) in sequences.items() if kind == 'row'})
    cuts = defaultdict(set)
    for _name, (sequence, start, end, _strand, inside) in observed.items():
        if not inside:
            cuts[sequence].update((start, end))
    seen = set()
    counts = defaultdict(int)

    def check(name, pieces):
        seen.add(name)
        if name not in sequences:
            problem('unexpected_sequence', name)
            return
        kind, value = sequences[name]
        offsets = [offset for offset, _bases in pieces]
        cursor = 0
        for offset, bases in pieces:
            if offset != cursor:
                problem('tiling', f'{name}: segment at {offset}, expected {cursor}')
                return
            if len(bases) > args.max_node_length:
                problem('node_length', f'{name}:{offset} has {len(bases)} bases')
            cursor += len(bases)
        spelled = ''.join(bases for _offset, bases in pieces)
        if kind == 'root':
            root = value
            if root.path not in readers:
                readers[root.path] = IndexedFasta(root.path, root.fai)
            source = readers[root.path].fetch(root.record or root.name, root.base,
                                              root.base + root.length)
        else:
            source = alleles.get(name, '')
        if spelled != source.upper():
            problem('bases', f'{name}: {len(spelled)} spelled bases differ from the source '
                             f'({len(source)} bases)')
        boundaries = set(offsets) | {cursor}
        missing = sorted(cuts.get(name, set()) - boundaries)
        if missing:
            problem('cuts', f'{name}: breakpoints {missing[:5]} are not segment boundaries')

    segments = defaultdict(list)    # SN -> [(SO, segment ID)] in file order
    gfa_links, link_lines, gfa_paths, walk_links = set(), 0, {}, 0
    current, pieces, expected_id = None, [], 1
    with open(args.gfa) as handle:
        for line in handle:
            if line.startswith('L\t'):
                fields = line.rstrip('\n').split('\t')
                if not any(field.startswith('SR:i:') for field in fields[6:]):
                    walk_links += 1      # added by the GAF walks (no rank): not a link walk
                    continue
                gfa_links.add(canonical(handle_of(fields[1], fields[2]), handle_of(fields[3], fields[4])))
                link_lines += 1
                continue
            if line.startswith('P\t'):
                fields = line.split('\t')
                gfa_paths[fields[1]] = [handle_of(step[:-1], step[-1]) for step in fields[2].split(',')]
                continue
            if not line.startswith('S\t'):
                continue
            fields = line.rstrip('\n').split('\t')
            if int(fields[1]) != expected_id:
                problem('ids', f'segment {fields[1]} where {expected_id} was expected')
                expected_id = int(fields[1])
            expected_id += 1
            tags = dict(field.split(':', 2)[0::2] for field in fields[3:])
            name = tags['SN']
            counts['segments'] += 1
            if name != current:
                if current is not None:
                    if current in seen:
                        problem('split_sequence', f'{current}: segments are not contiguous in the file')
                    check(current, pieces)
                current, pieces = name, []
            pieces.append((int(tags['SO']), fields[2]))
            segments[name].append((int(tags['SO']), int(fields[1])))
    if current is not None:
        check(current, pieces)
    for name in sorted(set(sequences) - seen):
        problem('missing_sequence', name)
    print(f'segments: {counts["segments"]} in {len(seen)} sequences ({len(sequences)} expected)')

    starts = {name: [offset for offset, _id in items] for name, items in segments.items()}

    def walk_handles(leaves):
        """Python leaves -> oriented segments (contiguous leaves of one source merged)."""
        merged = []
        for _kind, source, start, end, strand in leaves:
            if start == end:
                continue
            last = merged[-1] if merged else None
            if last and last[0] == source and last[3] == strand and (
                    (strand == '+' and last[2] == start) or (strand == '-' and last[1] == end)):
                last[1], last[2] = min(last[1], start), max(last[2], end)
            else:
                merged.append([source, start, end, strand])
        handles = []
        for source, start, end, strand in merged:
            if source not in segments:
                raise KeyError(f'walk uses {source}, which has no segments')
            first = bisect_right(starts[source], start) - 1
            last = bisect_right(starts[source], end - 1) - 1
            ids = [segments[source][index][1] for index in range(first, last + 1)]
            handles.extend((2 * value + (strand == '-')) for value in (ids if strand == '+' else ids[::-1]))
        return handles

    expected_links = set()
    for leaves in link_walks:
        steps = walk_handles(leaves)
        expected_links.update(canonical(a, b) for a, b in zip(steps, steps[1:]))
    describe = lambda link: f'{link[0] >> 1}{"+-"[link[0] & 1]} -> {link[1] >> 1}{"+-"[link[1] & 1]}'
    for link in sorted(expected_links - gfa_links):
        problem('link_missing', describe(link))
    for link in sorted(gfa_links - expected_links):
        problem('link_extra', describe(link))
    if link_lines != len(gfa_links):
        problem('link_duplicates', f'{link_lines - len(gfa_links)} repeated L lines')
    print(f'links: {len(gfa_links)} in the GFA, {len(expected_links)} from the Python link walks'
          f' (+{walk_links} added by GAF walks)')
    for name in sorted(set(path_walks) | set(gfa_paths)):
        if name not in gfa_paths or name not in path_walks:
            problem('path_missing' if name not in gfa_paths else 'path_extra', name)
        elif walk_handles(path_walks[name]) != gfa_paths[name]:
            problem('path_walk', f'{name}: P line differs from the Python path walk')
    print(f'paths: {len(gfa_paths)} P lines, {len(path_walks)} expected')
    if problems:
        for kind, items in sorted(problems.items()):
            print(f'{kind}: {len(items)}')
            for item in items[:args.show]:
                print(f'  {item}')
        print('result: FAILED')
        return 1
    print('result: ok')
    return 0


if __name__ == '__main__':
    raise SystemExit(main())
