#!/usr/bin/env python3
"""Alternative loci of one calling run: graph blocks the linear reference lacks.

Which graphs lack the linear reference is known only once the reference is
aligned (local_reference_templates.py). Their original blocks (the partition
intervals of original_intervals.bed, or of local_graphs.tsv) are the
alternative loci of this run, each one an independent contig from 0 to its
end, like a chromosome:

  - a block on its own keeps its block name (its graph path ID, e.g.
    g60825_chr2_34470763_34511500);
  - blocks of one source contig and strand that touch end to end are one
    locus named CONTIG#STRAND#START#END (e.g. chrEBV#+#35361#171823), unless
    the original fixed alternatives (--novel-loci, e.g. summary/novel_loci.fa)
    hold a record with exactly that source= interval and the same bases: the
    locus then keeps that record's name (e.g. HG38_h1_chrEBV_35361_171823).

A locus's bases are the bases of the templates the calls were made on, cut at
the block boundaries; template flanks outside the blocks are not part of any
locus. A template that holds none of its graph's blocks (its source is another
contig, e.g. a sample contig in a graph built on a backbone window) is kept
whole as a locus under its own template name.

Outputs (PREFIX defaults to alternative_loci next to the templates):
  PREFIX.tsv         name, length, kind, source interval, blocks of every locus
  PREFIX.fa(.fai)    locus sequences
  PREFIX.pieces.tsv  template interval -> locus interval, for VCF conversion
"""
from __future__ import annotations

import argparse
from collections import defaultdict
import os
from pathlib import Path
import re
import sys

from alternative_intervals import partition_prefix, read_original_bed, original_rows
from local_reference_templates import read_templates

SOURCE = re.compile(r'([^:]+):(.+):(\d+)-(\d+)([+-]?)')

LOCI_PROTOCOL = "alternative-loci-v1"
LOCI_COLUMNS = ("name", "length", "kind", "haplotype", "contig", "start", "end",
                "strand", "blocks")
PIECE_COLUMNS = ("template", "template_start", "template_end", "locus", "locus_start")


def merged_name(contig, strand, start, end):
    return f"{contig}#{strand}#{start}#{end}"


class OriginalAlternatives:
    """Records of the original fixed alternatives FASTA, by source= interval."""

    def __init__(self, path):
        from minsetref_core import IndexedFasta
        self.reader = IndexedFasta(path)
        self.by_source = {}
        with open(path) as handle:
            for raw in handle:
                if not raw.startswith('>'):
                    continue
                fields = raw[1:].split()
                for field in fields[1:]:
                    match = SOURCE.fullmatch(field[7:]) if field.startswith('source=') else None
                    if match:
                        haplotype, contig, start, end, strand = match.groups()
                        key = (haplotype, contig, strand or '+', int(start), int(end))
                        self.by_source.setdefault(key, fields[0])
                        break

    def name(self, haplotype, contig, strand, start, end, sequence):
        """(record name or None, 'exact' / 'different_bases' / 'absent')."""
        record = self.by_source.get((haplotype, contig, strand, start, end))
        if record is None:
            return None, 'absent'
        if self.reader.fetch(record, 0, self.reader.length(record)).upper() != sequence.upper():
            return None, 'different_bases'
        return record, 'exact'


def build_loci(templates, blocks, log=print, originals=None):
    """(loci, pieces). loci: dicts with LOCI_COLUMNS plus 'sequence';
    pieces: (template, template_start, template_end, locus, locus_start).
    originals: OriginalAlternatives naming merged loci that equal a record."""
    by_partition = defaultdict(list)
    by_prefix = defaultdict(list)
    for row in blocks:
        contig, start, end, partition, _score, strand, haplotype, path_id = row[:8]
        block = (haplotype, contig, strand, int(start), int(end), path_id)
        by_partition[partition].append(block)
        by_prefix[partition_prefix(partition)].append(block)
    owner = {}          # block -> (template, template_start, sequence)
    by_name = {template.name: template for template in templates}
    whole = []          # templates holding none of their blocks
    counts = defaultdict(int)
    imported = set()    # blocks that are an imported record's own block
    # An imported alternative is its own sequence: it owns only the block that
    # is exactly its record, and no input template is matched to that block
    # by provenance. Imported templates are therefore placed first.
    for template in sorted(templates, key=lambda t: (not t.fixed, t.name)):
        candidates = by_partition.get(template.graph_name)
        if candidates is None:
            candidates = by_prefix.get(partition_prefix(template.graph_name), [])
        length = template.end - template.start
        if len(template.sequence) != length:
            raise ValueError(f"{template.name}: template sequence is not loaded")
        held = 0
        for block in candidates:
            haplotype, contig, strand, start, end, _path_id = block
            if (haplotype, contig, strand) != (template.source_haplotype, template.source_contig,
                                               template.source_strand):
                continue
            if template.fixed and (start, end) != (template.source_start, template.source_end):
                continue
            if not template.fixed and block in imported:
                continue
            if not (template.source_start <= start and end <= template.source_end):
                if template.source_start < end and start < template.source_end:
                    counts['blocks_partly_in_template'] += 1
                continue
            offset = (start - template.source_start if strand == '+'
                      else template.source_end - end)
            if block in owner:
                counts['blocks_in_two_templates'] += 1
                held += 1
                continue
            owner[block] = (template.name, offset, template.sequence[offset:offset + end - start])
            if template.fixed:
                imported.add(block)
            held += 1
        if not held:
            whole.append(template)
    # Blocks touching end to end on one source contig and strand are one locus.
    groups = defaultdict(list)
    for block in owner:
        # An imported record is never merged with another block.
        groups[('imported', block) if block in imported else block[:3]].append(block)
    loci, pieces = [], []
    for key, members in sorted(groups.items(), key=lambda item: (item[0][0] == 'imported', str(item[0]))):
        haplotype, contig, strand = members[0][:3]
        members.sort(key=lambda block: (block[3], block[4]))
        runs = []
        for block in members:
            if runs and runs[-1][-1][4] == block[3]:
                runs[-1].append(block)
            elif runs and block[3] < runs[-1][-1][4]:
                raise ValueError(f'overlapping original blocks {runs[-1][-1][5]} and {block[5]}')
            else:
                runs.append([block])
        for run in runs:
            start, end = run[0][3], run[-1][4]
            ordered = run if strand == '+' else list(reversed(run))
            sequence = ''.join(owner[block][2] for block in ordered)
            name = run[0][5]
            if len(run) > 1:
                name = merged_name(contig, strand, start, end)
                if originals is not None:
                    record, status = originals.name(haplotype, contig, strand, start, end, sequence)
                    counts['merged_original_' + status] += 1
                    if status == 'different_bases':
                        log(f'[alternative-loci] warning: {name} has the source interval of an '
                            'original alternative but different bases; keeping the merged name')
                    name = record or name
            first = by_name[owner[run[0]][0]]
            loci.append(dict(name=name, length=end - start, kind='block' if len(run) == 1 else 'merged',
                             haplotype=haplotype, contig=contig, start=start, end=end,
                             strand=strand, blocks=','.join(block[5] for block in run),
                             sequence=sequence, graph=first.graph_name,
                             graph_prefix=first.graph_prefix,
                             fixed=all(by_name[owner[block][0]].fixed for block in run)))
            for block in run:
                template, offset, bases = owner[block]
                locus_start = block[3] - start if strand == '+' else end - block[4]
                pieces.append((template, offset, offset + len(bases), name, locus_start))
    for template in whole:
        length = template.end - template.start
        loci.append(dict(name=template.output_contig, length=length, kind='template',
                         haplotype=template.source_haplotype, contig=template.source_contig,
                         start=template.source_start, end=template.source_end,
                         strand=template.source_strand, blocks='.', sequence=template.sequence,
                         graph=template.graph_name, graph_prefix=template.graph_prefix,
                         fixed=template.fixed))
        pieces.append((template.name, 0, length, template.output_contig, 0))
    names = defaultdict(list)
    for locus in loci:
        names[locus['name']].append(locus)
    clashes = sorted(name for name, values in names.items() if len(values) > 1)
    if clashes:
        raise ValueError(f'alternative locus names collide: {clashes[:10]}')
    counts['loci'] = len(loci)
    counts['block_loci'] = sum(locus['kind'] == 'block' for locus in loci)
    counts['merged_loci'] = sum(locus['kind'] == 'merged' for locus in loci)
    counts['whole_template_loci'] = len(whole)
    counts['blocks'] = len(owner)
    log('[alternative-loci] ' + ', '.join(f'{key}={value}' for key, value in sorted(counts.items())))
    return loci, pieces


def _replace(path, write):
    path = Path(path)
    temporary = path.with_name(path.name + f'.tmp.{os.getpid()}')
    try:
        with temporary.open('w') as handle:
            write(handle)
        os.replace(temporary, path)
    finally:
        if temporary.exists():
            temporary.unlink()


def locus_records(loci, templates):
    """LocalTemplate records of the loci that are not already a template.

    A locus with a template's name and bases (a template that is exactly its
    block, e.g. an imported alternative) stays in local_reference_templates.fa
    only, so no name is in both files; the same name with other bases is an
    error.
    """
    from local_reference_templates import LocalTemplate
    sequences = {template.name: template.sequence for template in templates}
    records = []
    for locus in loci:
        known = sequences.get(locus['name'])
        if known is not None:
            if known.upper() != locus['sequence'].upper():
                raise ValueError(f"alternative locus {locus['name']} collides with a local "
                                 "template of the same name and different bases")
            continue
        records.append(LocalTemplate(
            name=locus['name'], graph_name=locus['graph'], graph_prefix=locus['graph_prefix'],
            graph_path=locus['name'], source_haplotype=locus['haplotype'],
            source_contig=locus['contig'], output_contig=locus['name'], start=0,
            end=locus['length'], source_strand=locus['strand'], source_start=locus['start'],
            source_end=locus['end'], sequence=locus['sequence'], fixed=locus['fixed']))
    return records


def lift_records(records, backbone_fasta, backbone, threads=1):
    """Place each locus on the backbone as a whole, by its own sequence.

    lift_local_templates leaves out non-imported records on excluded HG38
    contigs; they are kept here, unplaced."""
    import tempfile
    from dataclasses import replace
    from lift_local_templates import lift_templates
    with tempfile.TemporaryDirectory(prefix='.alternative_loci_lift.') as workdir:
        lifted = {record.name: record
                  for record in lift_templates(records, {}, backbone_fasta, backbone, workdir, threads)}
    return [lifted.get(record.name) or replace(
                record, backbone=backbone, lift_status='unmapped', reference='.',
                lift_note='not_lifted') for record in records]


def write_loci(prefix, loci, pieces, records=()):
    """PREFIX.tsv and PREFIX.pieces.tsv for every locus; PREFIX.fa(.fai) in
    the local_reference_templates.fa format for the records given."""
    from local_reference_templates import write_templates
    prefix = str(prefix)
    Path(prefix).parent.mkdir(parents=True, exist_ok=True)

    def table(handle):
        handle.write(f'#{LOCI_PROTOCOL}\n' + '\t'.join(LOCI_COLUMNS) + '\n')
        for locus in loci:
            handle.write('\t'.join(str(locus[column]) for column in LOCI_COLUMNS) + '\n')

    def piece_table(handle):
        handle.write(f'#{LOCI_PROTOCOL}\n' + '\t'.join(PIECE_COLUMNS) + '\n')
        for piece in sorted(pieces):
            handle.write('\t'.join(map(str, piece)) + '\n')

    write_templates(prefix + '.fa', list(records))
    _replace(prefix + '.tsv', table)
    _replace(prefix + '.pieces.tsv', piece_table)


def read_pieces(path):
    """template -> [(template_start, template_end, locus, locus_start)], sorted."""
    pieces = defaultdict(list)
    with open(path) as handle:
        for raw in handle:
            if raw.startswith('#') or raw.startswith(PIECE_COLUMNS[0] + '\t') or not raw.strip():
                continue
            template, t0, t1, locus, l0 = raw.rstrip('\n').split('\t')
            pieces[template].append((int(t0), int(t1), locus, int(l0)))
    for values in pieces.values():
        values.sort()
    return dict(pieces)


def read_loci(path):
    """name -> dict of LOCI_COLUMNS (ints for length/start/end)."""
    loci = {}
    with open(path) as handle:
        for raw in handle:
            if raw.startswith('#') or raw.startswith(LOCI_COLUMNS[0] + '\t') or not raw.strip():
                continue
            values = dict(zip(LOCI_COLUMNS, raw.rstrip('\n').split('\t')))
            for key in ('length', 'start', 'end'):
                values[key] = int(values[key])
            loci[values['name']] = values
    return loci


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('-t', '--templates', required=True,
                        help="this run's local_reference_templates.fa")
    source = parser.add_mutually_exclusive_group(required=True)
    source.add_argument('-b', '--original-intervals', help='GRAPH/summary/original_intervals.bed')
    source.add_argument('-s', '--summary', help='GRAPH/summary/local_graphs.tsv (when the BED is absent)')
    parser.add_argument('-n', '--novel-loci', default='',
                        help='original fixed alternatives FASTA with source= headers '
                             '(GRAPH/summary/novel_loci.fa); merged loci equal to a '
                             'record take its name')
    parser.add_argument('-o', '--prefix', default='',
                        help='output prefix (default: alternative_loci next to --templates)')
    parser.add_argument('-r', '--lift-reference', default='',
                        help='backbone FASTA to lift the loci onto (as lift_local_templates.py)')
    parser.add_argument('--reference-haplotype', default='', help='backbone name for --lift-reference')
    parser.add_argument('-j', '--threads', type=int, default=1)
    args = parser.parse_args(argv)
    if bool(args.lift_reference) != bool(args.reference_haplotype):
        parser.error('--lift-reference and --reference-haplotype go together')
    templates = read_templates(args.templates, load_sequences=True)
    blocks = (read_original_bed(args.original_intervals) if args.original_intervals
              else original_rows(args.summary))
    originals = OriginalAlternatives(args.novel_loci) if args.novel_loci else None
    loci, pieces = build_loci(templates, blocks,
                              log=lambda text: print(text, file=sys.stderr),
                              originals=originals)
    prefix = args.prefix or str(Path(args.templates).parent / 'alternative_loci')
    records = locus_records(loci, templates)
    if args.lift_reference and records:
        records = lift_records(records, args.lift_reference, args.reference_haplotype, args.threads)
    print(f'[alternative-loci] {len(records)} loci written to {prefix}.fa; '
          f'{len(loci) - len(records)} are local templates already', file=sys.stderr)
    write_loci(prefix, loci, pieces, records)
    return 0


if __name__ == '__main__':
    raise SystemExit(main())
