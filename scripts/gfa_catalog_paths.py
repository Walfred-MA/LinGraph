"""Whole alternative paths from the graph catalog, placed on the backbone.

VCF targets outside the backbone (local templates, alternatives, novels) are
found in the graph folder's ``summary/alternatives.fasta`` by source interval,
never by name. A target's source interval comes from the VCF header's
``##alternativeLocus`` line, or from its own catalog header.

Graph construction can cut one alternative sequence into several local loci.
On each source contig, all non-backbone catalog records whose source intervals
overlap or touch are stitched back into one forward path covering their union;
overlapping bases must agree, or the pieces stay separate. Every VCF target
maps onto its merged path at its source offset (reversed for a '-' source).

Placement: a path made of exactly one record uses that record's catalog
placement (header fields 3-4, a graph CIGAR on a backbone contig). Any other
path is aligned once, as a whole, with lift_local_templates.lift_templates.
Only ends aligned exactly at the path boundary attach to the backbone; other
paths stay in the graph unplaced. Backbone records map onto the backbone root.
"""
from collections import defaultdict
from pathlib import Path
import re

from gfa_query_anchors import reverse_complement
from gfa_source_catalog import SOURCE, SourceAlias, indexed_header
from minsetref_core import IndexedFasta

_META = re.compile(r'(\w+)=("(?:[^"\\]|\\.)*"|[^,>]*)')
_SEGMENT = re.compile(r'([<>])([^:<>]+):([^<>]+)')
_OP = re.compile(r'(\d+)([=XMIDHSN])((?<=[IX])[A-Za-z]*)?')


def read_alternative_loci(vcf_paths):
    """ID -> (sample, contig, start, end, strand) from ##alternativeLocus."""
    import gzip
    loci = {}
    for path in vcf_paths:
        opener = gzip.open if str(path).endswith('.gz') else open
        with opener(path, 'rt') as handle:
            for line in handle:
                if not line.startswith('##'):
                    break
                if not line.startswith('##alternativeLocus=<'):
                    continue
                values = {}
                for key, value in _META.findall(line[len('##alternativeLocus=<'):]):
                    if value.startswith('"'):
                        value = value[1:-1].replace('\\"', '"').replace('\\\\', '\\')
                    values[key] = value
                try:
                    loci[values['ID']] = (values['SourceHaplotype'], values['SourceContig'],
                                          int(values['SourceStart']), int(values['SourceEnd']),
                                          values.get('SourceStrand', '+') or '+')
                except (KeyError, ValueError):
                    continue
    return loci


def placement_lift(cigar_field, roots, backbone, reverse, length):
    """Lift dict from a single-segment catalog graph CIGAR, or None to relift.

    The CIGAR must parse completely and consume exactly ``length`` query
    bases; otherwise its end positions would not be the record's ends.
    """
    match = _SEGMENT.fullmatch(cigar_field or '')
    if match is None:
        return None
    marker, target, body = match.groups()
    root = roots.get(target)
    if root is None or root.kind != 'reference':
        return None  # placed on another genome: not usable for this backbone
    tokens = list(_OP.finditer(body))
    if not tokens or ''.join(token[0] for token in tokens) != body:
        return None
    ops = [(int(token[1]), token[2]) for token in tokens]
    lead = ops[0][0] if ops[0][1] == 'H' else 0
    body_ops = [(n, op) for n, op in ops if op != 'H']
    consumed = sum(n for n, op in body_ops if op in '=XMDN')
    query_ops = [op for n, op in body_ops if op in '=XMIS']
    if (not query_ops or lead + consumed > root.length
            or sum(n for n, op in body_ops if op in '=XMIS') != length):
        return None
    left_ok, right_ok = query_ops[0] in '=XM', query_ops[-1] in '=XM'
    strand = '+' if marker == '>' else '-'
    low, high = ((lead, lead + consumed) if strand == '+' else
                 (root.length - lead - consumed, root.length - lead))
    left = (target, low if strand == '+' else high, strand)
    right = (target, high if strand == '+' else low, strand)
    if reverse:
        # The merged path is the reverse complement of this record.
        flip = '-' if strand == '+' else '+'
        left, right = (target, right[1], flip), (target, left[1], flip)
        left_ok, right_ok = right_ok, left_ok
    status = 'mapped' if left_ok and right_ok else 'one_sided' if left_ok or right_ok else 'unmapped'
    return dict(status=status, backbone=backbone, left=left if left_ok else None,
                right=right if right_ok else None, note='catalog placement',
                placements=f'{target}:{low}-{high}:{strand}')


def _catalog_records(reader, backbone, stats):
    """(sample, contig) -> sorted records; name -> record; backbone excluded."""
    by_contig, by_name = defaultdict(list), {}
    for name, index in reader.index.items():
        fields = indexed_header(reader, name)
        match = SOURCE.fullmatch(fields[1]) if len(fields) > 1 else None
        if match is None:
            stats['records_without_source'] += 1
            continue
        sample, contig, start, end, strand = match.groups()
        record = dict(name=name, sample=sample, contig=contig, start=int(start),
                      end=int(end), strand=strand or '+', length=index[0],
                      placement=fields[3] if len(fields) > 3 else '.')
        if record['end'] - record['start'] != record['length']:
            stats['records_length_mismatch'] += 1
            continue
        by_name[name] = record
        if sample != backbone:
            by_contig[sample, contig].append(record)
    for records in by_contig.values():
        records.sort(key=lambda record: (record['start'], record['end'], record['name']))
    return by_contig, by_name


def _components(records, targets):
    """Records grouped by overlapping/touching source intervals (coordinates
    only), keeping groups that contain at least one target interval."""
    groups, current, end = [], [], -1
    for record in records:
        if current and record['start'] > end:
            groups.append(current)
            current = []
        current.append(record)
        end = max(end, record['end']) if len(current) > 1 else record['end']
    if current:
        groups.append(current)
    kept = []
    for group in groups:
        low, high = group[0]['start'], max(record['end'] for record in group)
        if any(low <= a and b <= high for a, b, _name in targets):
            kept.append(group)
    return kept


def _stitch(records, reader, stats):
    """Merge overlapping/touching records into [start, end, sequence, members].

    Overlapping bases must agree; a disagreeing record starts a new piece.
    """
    merged = []
    for record in records:
        sequence = ''.join(reader.fetch(record['name'], 0, record['length']))
        if record['strand'] == '-':
            sequence = reverse_complement(sequence)
        current = merged[-1] if merged else None
        if current and record['start'] <= current[1]:
            overlap = current[1] - record['start']
            shared = min(overlap, record['length'])
            existing = _tail(current, record['start'], shared)
            if existing.upper() == sequence[:shared].upper():
                if record['end'] > current[1]:
                    current[2].append(sequence[shared:])
                    current[1] = record['end']
                current[3].append(record)
                continue
            stats['overlap_mismatch'] += 1
        merged.append([record['start'], record['end'], [sequence], [record], 0])
    return [[start, end, ''.join(chunks), members] for start, end, chunks, members, _ in merged]


def _tail(piece, start, size):
    """Bases [start, start+size) of a piece held as a list of chunks."""
    chunks = piece[2]
    if len(chunks) > 64:
        chunks[:] = [''.join(chunks)]
    position = piece[1]
    wanted_end = start + size
    parts = []
    for chunk in reversed(chunks):
        chunk_start = position - len(chunk)
        if chunk_start < wanted_end and position > start:
            parts.append(chunk[max(0, start - chunk_start):wanted_end - chunk_start])
        if chunk_start <= start:
            break
        position = chunk_start
    return ''.join(reversed(parts))


def resolve_catalog_paths(events, roots, reachable, catalog, vcf_paths, backbone,
                          backbone_fasta, output_prefix, root_type, log, *,
                          threads=1, relift=None, needed=None):
    """Add merged catalog paths as roots; return {VCF name: SourceAlias}.

    With a precomputed ``needed`` set, ``events`` only needs to contain the
    variant IDs that could collide with catalog names (``alt_*`` and
    ``##alternativeLocus`` IDs), not every variant.
    """
    if needed is None:
        needed = set()
        for name in reachable:
            event = events[name]
            needed.add(event.chrom)
            needed.update(run.target for run in event.runs if run.target)
        needed.difference_update(events)
    needed = set(needed) - set(roots)
    loci = read_alternative_loci(vcf_paths)
    # Templates the samples were called against, even without variants, so
    # GAF walks can pass through them. Best effort: skipped if not covered.
    optional = set(loci) - needed - set(events) - set(roots)
    if not needed and not optional:
        return {}
    if not backbone:
        raise ValueError('--alternative-catalog needs --reference-haplotype to tell '
                         'backbone records from alternatives')
    if not Path(str(catalog) + '.fai').is_file():
        raise ValueError(f'alternative catalog requires an index: {catalog}.fai; run samtools faidx')
    stats = defaultdict(int)
    aliases, unresolved = {}, []
    with IndexedFasta(catalog) as reader:
        by_contig, by_name = _catalog_records(reader, backbone, stats)
        wanted = {}
        for name in sorted(needed | optional):
            if name in loci:
                wanted[name] = loci[name]
            elif name in by_name:
                record = by_name[name]
                wanted[name] = (record['sample'], record['contig'], record['start'],
                                record['end'], record['strand'])
            elif name in needed:
                unresolved.append((name, 'no ##alternativeLocus line or catalog record'))
        for name, (sample, contig, start, end, strand) in wanted.items():
            if sample == backbone:
                if contig in roots and roots[contig].kind == 'reference' and end <= roots[contig].length:
                    aliases[name] = SourceAlias(contig, start, end, strand)
                    base = roots[contig]
                    roots[name] = root_type(name, base.path, base.fai, end - start, False,
                                            base.kind, -1)
                    stats['backbone_names'] += 1
                elif name in needed:
                    unresolved.append((name, f'backbone interval {contig}:{start}-{end} not in -r'))
        by_key = defaultdict(list)
        for name, (sample, contig, start, end, _strand) in wanted.items():
            if name not in aliases and sample != backbone:
                by_key[sample, contig].append((start, end, name))
        paths = []
        for key in sorted(by_key):
            for group in _components(by_contig.get(key, []), by_key[key]):
                paths.extend((key, piece) for piece in _stitch(group, reader, stats))


        fasta = Path(str(output_prefix) + '.alternative_paths.fa')
        fasta.parent.mkdir(parents=True, exist_ok=True)
        kept, index_rows = [], []
        with fasta.open('w') as out:
            for (sample, contig), (start, end, sequence, members) in paths:
                covered = [name for a, b, name in by_key[sample, contig]
                           if start <= a <= b <= end]
                if not covered:
                    continue
                path_name = f'alt_{sample}_{contig}_{start}_{end}'
                if path_name in roots or path_name in events:
                    raise ValueError(f'alternative path name {path_name!r} collides with another path')
                out.write(f'>{path_name}\n')
                offset = out.tell()
                out.write(sequence + '\n')
                index_rows.append(f'{path_name}\t{len(sequence)}\t{offset}\t'
                                  f'{max(1, len(sequence))}\t{len(sequence) + 1}\n')
                kept.append((path_name, sample, contig, start, end, sequence, members, covered))
        Path(str(fasta) + '.fai').write_text(''.join(index_rows))

    next_order = max((root.order for root in roots.values()), default=-1) + 1
    lifts = {}
    to_lift = []
    for path_name, _sample, _contig, _start, _end, sequence, members, _covered in kept:
        lift = None
        if len(members) == 1:
            lift = placement_lift(members[0]['placement'], roots, backbone,
                                  reverse=members[0]['strand'] == '-',
                                  length=members[0]['length'])
        if lift is None:
            to_lift.append((path_name, _sample, _contig, _start, _end, sequence))
        else:
            lifts[path_name] = lift
            stats['catalog_placements'] += 1
    if to_lift:
        log(f'Aligning {len(to_lift)} merged alternative path(s) without a backbone placement')
        lifted = (relift or _relift)(to_lift, roots, backbone, backbone_fasta,
                                     Path(str(output_prefix) + '.alternative_lift'), threads)
        lifts.update(lifted)
        stats['relifted'] += len(to_lift)

    audit = Path(str(output_prefix) + '.alternative_paths.tsv')
    with audit.open('w') as out:
        out.write('path\tsample\tcontig\tstart\tend\trecords\tlift_status\tplacement\tnote\tvcf_names\n')
        for path_name, sample, contig, start, end, sequence, members, covered in kept:
            lift = lifts.get(path_name) or dict(status='unmapped', backbone=backbone,
                                                left=None, right=None,
                                                note='not lifted', placements='.')
            roots[path_name] = root_type(path_name, str(fasta), str(fasta) + '.fai', len(sequence),
                                         True, 'alternative', next_order, lift)
            next_order += 1
            for name in covered:
                s, c, a, b, strand = wanted[name]
                aliases[name] = SourceAlias(path_name, a - start, b - start, strand)
                # Alias roots only record a length; the resolver translates
                # them onto the merged path, which carries the lift.
                roots[name] = root_type(name, str(fasta), str(fasta) + '.fai', b - a,
                                        False, 'alternative', -1)
            stats[f'lift_{lift["status"]}'] += 1
            out.write(f'{path_name}\t{sample}\t{contig}\t{start}\t{end}\t'
                      f'{",".join(m["name"] for m in members)}\t{lift["status"]}\t'
                      f'{lift.get("placements", ".")}\t{lift.get("note", ".")}\t{",".join(covered)}\n')
    for name in wanted:
        if name not in aliases:
            if name in needed:
                unresolved.append((name, 'source interval not covered by merged catalog records'))
            else:
                stats['header_only_not_in_catalog'] += 1
    stats['header_only_included'] = sum(name in aliases for name in optional)
    if unresolved:
        preview = '; '.join(f'{name}: {reason}' for name, reason in sorted(unresolved)[:10])
        raise ValueError(f'{len(unresolved)} VCF target(s) not found in {catalog}: {preview}')
    log(f'Alternative catalog: {len(kept)} merged path(s) for {len(aliases)} VCF name(s); '
        + ', '.join(f'{key}={value}' for key, value in sorted(stats.items())))
    return aliases


def _relift(paths, roots, backbone, backbone_fasta, workdir, threads):
    """Align whole merged paths as fixed sequences; ends attach only if unique."""
    from gfa_source_catalog import template_lift
    from lift_local_templates import lift_templates
    from local_reference_templates import LocalTemplate
    templates = [LocalTemplate(
        name=name, graph_name=name, graph_prefix=name, graph_path=name,
        source_haplotype=sample, source_contig=contig, output_contig=name,
        start=0, end=len(sequence), source_strand='+', source_start=start,
        source_end=end, sequence=sequence, fixed=True)
        for name, sample, contig, start, end, sequence in paths]
    lifts = {}
    for template in lift_templates(templates, {}, backbone_fasta, backbone, workdir, threads):
        fields = [template.name, f'backbone={template.backbone}',
                  f'lift_status={template.lift_status}', f'reference={template.reference}',
                  f'lift_note={template.lift_note}']
        lifts[template.name] = template_lift(fields, roots, backbone)
    return lifts
