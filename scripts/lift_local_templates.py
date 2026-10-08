#!/usr/bin/env python3
"""Place each local template on the backbone by its own stored sequence."""
from collections import defaultdict
from dataclasses import replace
import argparse
import logging
from pathlib import Path
import tempfile

from assembly_contigs import accepted_contig
from gfa_query_anchors import reverse_complement
from local_reference_templates import read_templates, write_templates
from minsetref_core import IndexedFasta, wrap_fasta
from minsetref_localize import LocalizeResult, Placement, anchor_blocks_from_hits


def annotate(template, result, backbone, source_forward=True):
    """Retain public orientation for source-forward or fixed-sequence alignments."""
    placements = result.placements
    note = result.note or '.'
    orientation = template.source_strand if source_forward else '+'
    if orientation == '-':
        placements = list(reversed(placements))
        note = {'one_side_up': 'one_side_down', 'one_side_down': 'one_side_up'}.get(note, note)
    coordinates = ';'.join(
        f'{p.main}:{p.start}-{p.end}:{"+" if p.strand == orientation else "-"}'
        for p in placements) or '.'
    return replace(template, name=template.output_contig, storage_strand='+',
                   backbone=backbone, lift_status=result.status,
                   reference=coordinates, lift_note=note)


def localize_fixed_template(length, blocks):
    """Attach only uniquely aligned endpoints of the independent sequence."""
    left = {(b.main, b.strand, b.map_qpos(0)) for b in blocks if b.q0 == 0}
    right = {(b.main, b.strand, b.map_qpos(length)) for b in blocks if b.q1 == length}
    if len(left) > 1 or len(right) > 1:
        return LocalizeResult('unmapped', [], 'ambiguous template endpoints')
    if not left and not right:
        return LocalizeResult('unmapped', [], 'no aligned template endpoint')
    if left and right:
        a, b = next(iter(left)), next(iter(right))
        if a[:2] != b[:2] or (a[2] > b[2] if a[1] == '+' else a[2] < b[2]):
            return LocalizeResult('unmapped', [], 'noncollinear template endpoints')
        return LocalizeResult('mapped', [Placement(a[0], min(a[2], b[2]), max(a[2], b[2]), a[1])])
    side = next(iter(left or right))
    return LocalizeResult('one_sided', [Placement(side[0], side[2], side[2], side[1])],
                          'one_side_up' if left else 'one_side_down')


def lift_templates(templates, sources, backbone_fasta, backbone, workdir,
                   threads=1, flank=10000, min_identity=95.0, min_segment=300,
                   collect=None):
    """Place every template by its own stored sequence.

    A template is an independent piece: input assemblies are never opened, so
    a different assembly version under the same sample name cannot change or
    block its lift.  ``sources`` and ``flank`` are kept for call compatibility.
    """
    templates = [t for t in templates
                 if t.fixed or accepted_contig(t.source_haplotype, t.source_contig)]
    results, queries = {}, {}
    workdir = Path(workdir)
    workdir.mkdir(parents=True, exist_ok=True)
    anchored = workdir / 'template_sequences.fa'
    with IndexedFasta(backbone_fasta) as reference, anchored.open('w') as out:
        for index, template in enumerate(templates):
            sequence = template.sequence.upper()
            contig, start, end = template.source_contig, template.source_start, template.source_end
            # A backbone-sourced template is placed directly only when its own
            # bases equal the backbone at the recorded interval.
            if (not template.fixed and template.source_haplotype == backbone
                    and contig in reference.index and 0 <= start < end <= reference.index[contig][0]):
                core = reference.fetch(contig, start, end).upper()
                if sequence in (core, reverse_complement(core)):
                    strand = '+' if sequence == core else '-'
                    results[index] = LocalizeResult('mapped', [Placement(contig, start, end, strand)], 'direct source')
                    continue
            query = f'template_{index}'
            queries[query] = index, 0, len(template.sequence)
            out.write(f'>{query}\n{wrap_fasta(template.sequence)}\n')
        if queries:
            # The aligner opens this file separately; flush the final buffered
            # records even when the complete query FASTA is smaller than 8 KiB.
            out.flush()
            if collect is None:
                from minsetref_align import collect_alignment_hits
                collect = collect_alignment_hits
            hits = collect(str(anchored), backbone_fasta, str(workdir / 'align'), threads,
                           min_identity, min_segment, logging.getLogger('lift-local-templates'),
                           skip_blastn=True, broad_aligner='minimap2',
                           minimap_params=dict(preset=None, p=0.001, N=100, f=0.001, K='100M'))
            grouped = defaultdict(list)
            for hit in hits:
                grouped[hit.query_id].append(hit)
            for query, (index, start, end) in queries.items():
                blocks = [b for b in anchor_blocks_from_hits(grouped[query], query)
                          if b.main in reference.index and accepted_contig(backbone, b.main)]
                results[index] = localize_fixed_template(end, blocks)
    return [annotate(template, results[index], backbone, source_forward=False)
            for index, template in enumerate(templates)]


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--input', required=True)
    parser.add_argument('--reference-fasta', required=True)
    parser.add_argument('--reference-haplotype', required=True)
    # Accepted for existing pipeline commands; templates never read input assemblies.
    parser.add_argument('--assemblies', help='ignored: templates are lifted from their own sequence')
    parser.add_argument('--source', action='append', nargs=2, default=[], metavar=('SAMPLE', 'FASTA'),
                        help='ignored: templates are lifted from their own sequence')
    parser.add_argument('--output', required=True)
    parser.add_argument('--threads', type=int, default=1)
    parser.add_argument('--flank', type=int, default=10000)
    args = parser.parse_args(argv)
    if args.threads < 1 or args.flank < 1:
        parser.error('threads and flank must be positive')
    output = Path(args.output).resolve()
    output.parent.mkdir(parents=True, exist_ok=True)
    logging.basicConfig(level=logging.INFO, format='[lift-local-templates] %(message)s')
    with tempfile.TemporaryDirectory(prefix='template-lift-', dir=output.parent) as work:
        templates = lift_templates(read_templates(args.input, load_sequences=True), {},
                                   args.reference_fasta, args.reference_haplotype,
                                   work, args.threads, args.flank)
        write_templates(str(output), templates)
    Path(str(output)+'.lift.complete').write_text(
        f'local-template-lift-v2-fixed-alternatives:{args.reference_haplotype}:{args.reference_fasta}\n')
    counts = defaultdict(int)
    for template in templates:
        counts[template.lift_status] += 1
    print(f'[lift-local-templates] {len(templates)} templates: {dict(counts)}; backbone={args.reference_haplotype}')
    return 0


if __name__ == '__main__':
    raise SystemExit(main())
