"""Translate local graph coordinates without exporting duplicate local paths."""
from collections import defaultdict
from contextlib import ExitStack
from dataclasses import dataclass
import hashlib
import os
from pathlib import Path
import re

from gfa_query_anchors import interval_chunks
from minsetref_core import IndexedFasta
from assembly_contigs import accepted_contig

SOURCE = re.compile(r'([^:]+):(.+):(\d+)-(\d+)([+-]?)$')


@dataclass(frozen=True)
class SourceAlias:
    target: str
    start: int
    end: int
    strand: str

    def interval(self, start, end, orientation='+'):
        low, high = ((self.start + start, self.start + end) if self.strand == '+'
                     else (self.end - end, self.end - start))
        strand = '+' if orientation == self.strand else '-'
        return self.target, low, high, strand


def indexed_header(reader, name):
    """Read only the header preceding the indexed sequence, never scan the catalog."""
    offset = reader.index[name][1]
    width = 4096
    while True:
        start = max(0, offset-width)
        block = os.pread(reader._fd, offset-start, start).rstrip(b'\r\n')
        marker = block.rfind(b'\n>')
        header = block[marker+2:] if marker >= 0 else (block[1:] if block.startswith(b'>') else b'')
        fields = header.decode('utf-8').split()
        if fields and fields[0] == name and b'\n' not in header:
            return fields
        if start == 0 or width >= 1024 * 1024:
            raise ValueError(f'{reader.path}: cannot read indexed header for {name!r}; check .fai')
        width *= 2


def source_interval(fields, length):
    values = [token[7:] for token in fields[1:] if token.startswith('source=')]
    attrs = dict(token.split('=', 1) for token in fields[1:] if '=' in token)
    if not values and all(key in attrs for key in ('source_haplotype', 'source_contig', 'source_start', 'source_end', 'source_strand')):
        values = [f'{attrs["source_haplotype"]}:{attrs["source_contig"]}:{attrs["source_start"]}-{attrs["source_end"]}{attrs["source_strand"]}']
    if not values and len(fields) > 1:
        values = [fields[1]]
    if not values:
        return None
    match = SOURCE.fullmatch(values[0])
    if not match:
        return None
    sample, contig, start, end, strand = match.groups()
    start, end = int(start), int(end)
    if end-start != length:
        raise ValueError(f'{fields[0]}: source interval length {end-start} differs from FASTA length {length}')
    return sample, contig, start, end, strand or '+'


def sequence_digest(reader, name, start, end, strand='+'):
    digest = hashlib.sha256()
    observed = 0
    for chunk in interval_chunks(reader, name, start, end, strand):
        observed += len(chunk)
        digest.update(chunk.upper().encode('ascii'))
    if observed != end-start:
        raise ValueError(f'{name}: truncated FASTA interval {start}-{end}')
    return digest.digest()


def template_lift(fields, roots, backbone):
    attrs = dict(token.split('=', 1) for token in fields[1:] if '=' in token)
    status = attrs.get('lift_status')
    if not attrs.get('backbone') or status not in ('mapped', 'one_sided', 'both', 'unmapped'):
        raise ValueError(f'{fields[0]}: fallback template must be lifted before export; run lift_local_templates.py')
    if backbone and attrs['backbone'] != backbone:
        raise ValueError(f'{fields[0]}: stale template lift to {attrs["backbone"]}, expected {backbone}')
    lift = dict(status=status, backbone=attrs['backbone'], left=None, right=None,
                note=attrs.get('lift_note', '.'), placements=attrs.get('reference', '.'))
    if status not in ('mapped', 'one_sided'):
        return lift  # Ambiguous/disagreeing/unmapped placements remain unplaced.
    match = re.fullmatch(r'(.+):(\d+)-(\d+):([+-])', attrs.get('reference', ''))
    if not match:
        raise ValueError(f'{fields[0]}: malformed lifted reference placement')
    target, start, end, strand = match.groups()
    start, end = int(start), int(end)
    if (target not in roots or roots[target].kind != 'reference' or
            not 0 <= start <= end <= roots[target].length):
        raise ValueError(f'{fields[0]}: lift target {target}:{start}-{end} is outside the selected backbone')
    left = (target, start if strand == '+' else end, strand)
    right = (target, end if strand == '+' else start, strand)
    if status == 'mapped':
        lift.update(left=left, right=right)
    elif start != end:
        raise ValueError(f'{fields[0]}: one-sided lift must identify one breakpoint')
    elif attrs.get('lift_note') == 'one_side_up':
        lift['left'] = left
    elif attrs.get('lift_note') == 'one_side_down':
        lift['right'] = right
    else:
        raise ValueError(f'{fields[0]}: one-sided lift lacks side metadata')
    return lift


def resolve_local_sources(events, roots, reachable, sources, catalogs, log,
                          template_catalogs=(), backbone=None):
    """Map needed catalog records to contained, sequence-identical included roots."""
    from gfa_interval_pipeline import Root

    counts = defaultdict(lambda: [0, 0])
    for root in roots.values():
        if root.emit:
            counts[root.kind][0] += 1
            counts[root.kind][1] += root.length
    log('Included FASTA paths: ' + '; '.join(
        f'{kind}={count} paths/{bases} bases' for kind, (count, bases) in sorted(counts.items())))
    needed = set()
    for name in reachable:
        event = events[name]
        needed.add(event.chrom)
        needed.update(run.target for run in event.runs if run.target)
    needed.difference_update(events)
    needed.difference_update(roots)
    if (not needed or not catalogs) and not template_catalogs:
        return {}

    aliases = {}
    with ExitStack() as stack:
        readers = {}

        def reader_for(path, fai=None):
            key = str(Path(path).resolve()), fai
            if key not in readers:
                readers[key] = stack.enter_context(IndexedFasta(path, fai))
            return readers[key]

        included = defaultdict(list)
        for name, root in roots.items():
            reader = reader_for(root.path, root.fai)
            source = source_interval(indexed_header(reader, name), root.length)
            if source is not None:
                included[source[:2]].append((name, source))
            else:
                # A plain reference supplied directly as a query assembly has
                # an unambiguous whole-contig source even without source= tags.
                for sample, (path, _fai) in sources.items():
                    if Path(path).resolve() == Path(root.path).resolve():
                        included[sample, name].append((name, (sample, name, 0, root.length, '+')))
                if backbone and root.kind == 'reference' and not included[backbone, name]:
                    included[backbone, name].append((name, (backbone, name, 0, root.length, '+')))

        catalog_roles = [(str(path), True) for path in dict.fromkeys(template_catalogs)]
        catalog_roles.extend((str(path), False) for path in dict.fromkeys(catalogs))
        next_order = max((root.order for root in roots.values()), default=-1) + 1
        added = reused = 0
        for catalog, is_template in catalog_roles:
            # Loading the FAI is bounded by record count; require it for this
            # potentially multi-GB lookup input rather than scanning sequences.
            if not Path(catalog + '.fai').is_file():
                raise ValueError(f'local path catalog requires an adjacent index: {catalog}.fai; run samtools faidx')
            reader = reader_for(catalog)
            for name in sorted(reader.index if is_template else needed.intersection(reader.index)):
                if name in roots:
                    if is_template:
                        template_lift(indexed_header(reader, name), roots, backbone)
                    prior = aliases.get(name, SourceAlias(name, 0, roots[name].length, '+'))
                    root = roots[prior.target]
                    if (reader.index[name][0] != prior.end-prior.start or
                            sequence_digest(reader, name, 0, reader.index[name][0]) !=
                            sequence_digest(reader_for(root.path, root.fai), prior.target,
                                            prior.start, prior.end, prior.strand)):
                        raise ValueError(f'{name}: conflicting sequence definitions in FASTA catalogs')
                    continue
                length = reader.index[name][0]
                fields = indexed_header(reader, name)
                source = source_interval(fields, length)
                if source is None:
                    raise ValueError(f'{name}: local path lacks source SAMPLE:CONTIG:START-END metadata')
                sample, contig, low, high, strand = source
                if not accepted_contig(sample, contig):
                    if is_template:
                        continue
                    raise ValueError(f'{name}: excluded HG38 non-primary source contig {contig}')
                lift = template_lift(fields, roots, backbone) if is_template else None
                matching = []
                local_digest = None
                for target, origin in included[sample, contig]:
                    _sample, _contig, start, end, orientation = origin
                    if not start <= low <= high <= end:
                        continue
                    a, b = ((low-start, high-start) if orientation == '+' else (end-high, end-low))
                    alias = SourceAlias(target, a, b, '+' if strand == orientation else '-')
                    if local_digest is None:
                        local_digest = sequence_digest(reader, name, 0, length)
                    root = roots[target]
                    if local_digest == sequence_digest(reader_for(root.path, root.fai), target, a, b, alias.strand):
                        matching.append(alias)
                matching = set(matching)
                if not matching:
                    if is_template:
                        roots[name] = Root(name, catalog, None, length, True, 'alternative', next_order, lift)
                        included[sample, contig].append((name, source))
                        next_order += 1
                        added += 1
                        continue
                    hint = ('no --local-reference-templates supplied; lift the selected fallback '
                            'catalog with lift_local_templates.py and pass it with '
                            '--local-reference-templates' if not template_catalogs else
                            'check that --local-reference-templates contains the selected fallback '
                            'for this source interval and the same backbone used for cohort calling')
                    raise ValueError(f'{name}: no sequence-identical covering source in included -r/-a FASTAs '
                                     f'or selected templates for {sample}:{contig}:{low}-{high}{strand}; '
                                     f'the local catalog is lookup-only and cannot add this path; {hint}')
                if len(matching) != 1:
                    raise ValueError(f'{name}: ambiguous source mapping to multiple included paths')
                alias = matching.pop()
                if name in aliases and aliases[name] != alias:
                    raise ValueError(f'{name}: conflicting mappings in local path catalogs')
                aliases[name] = alias
                root = roots[alias.target]
                roots[name] = Root(name, root.path, root.fai, length, is_template,
                                   'alternative' if is_template else root.kind,
                                   next_order if is_template else -1, lift)
                if is_template:
                    next_order += 1
                    reused += 1
        if template_catalogs:
            log(f'Fallback templates: {added} additional sequence paths, {reused} paths reusing included sequence')
    log(f'Resolved {len(aliases)} local path names onto included paths; no local catalog paths exported')
    return aliases
