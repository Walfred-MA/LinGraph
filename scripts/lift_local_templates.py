#!/usr/bin/env python3
"""Place fixed template sequences or refinable templates with source flanks."""
from collections import defaultdict
from dataclasses import replace
import argparse
import logging
from pathlib import Path
import tempfile

from assembly_contigs import accepted_contig
from gfa_query_anchors import read_query_list, reverse_complement
from local_reference_templates import read_templates, write_templates
from minsetref_core import IndexedFasta, wrap_fasta
from minsetref_localize import LocalizeResult, Placement, anchor_blocks_from_hits, localize_insertion


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
    templates = [t for t in templates if accepted_contig(t.source_haplotype, t.source_contig)]
    results, queries = {}, {}
    by_sample = defaultdict(list)
    for index, template in enumerate(templates):
        by_sample[template.source_haplotype].append((index, template))
    workdir = Path(workdir)
    workdir.mkdir(parents=True, exist_ok=True)
    anchored = workdir / 'source_flanks.fa'
    with IndexedFasta(backbone_fasta) as reference, anchored.open('w') as out:
        for sample, group in sorted(by_sample.items()):
            for index, template in group:
                if template.fixed:
                    query = f'template_{index}'
                    queries[query] = index, 0, len(template.sequence)
                    out.write(f'>{query}\n{wrap_fasta(template.sequence)}\n')
            group = [(index, template) for index, template in group if not template.fixed]
            if not group:
                continue
            source = sources.get(sample)
            if source is None:
                for index, _template in group:
                    results[index] = LocalizeResult('unmapped', [], 'source assembly unavailable')
                continue
            with IndexedFasta(*source) as reader:
                for index, template in group:
                    contig, start, end = template.source_contig, template.source_start, template.source_end
                    if contig not in reader.index or not 0 <= start < end <= reader.index[contig][0]:
                        raise ValueError(f'{template.name}: source interval is outside its assembly')
                    core = reader.fetch(contig, start, end)
                    expected = core if template.source_strand == '+' else reverse_complement(core)
                    if expected.upper() != template.sequence.upper():
                        raise ValueError(f'{template.name}: template sequence differs from source provenance')
                    if sample == backbone and contig in reference.index:
                        if end <= reference.index[contig][0] and reference.fetch(contig, start, end).upper() == core.upper():
                            results[index] = LocalizeResult('mapped', [Placement(contig, start, end, '+')], 'direct source')
                            continue
                    low, high = max(0, start-flank), min(reader.index[contig][0], end+flank)
                    query = f'template_{index}'
                    queries[query] = index, start-low, end-low
                    out.write(f'>{query}\n{wrap_fasta(reader.fetch(contig, low, high))}\n')
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
                           minimap_params=dict(preset='asm5', p=0.001, N=100, f=0.001, K='100M'))
            grouped = defaultdict(list)
            for hit in hits:
                grouped[hit.query_id].append(hit)
            for query, (index, start, end) in queries.items():
                blocks = [b for b in anchor_blocks_from_hits(grouped[query], query)
                          if b.main in reference.index and accepted_contig(backbone, b.main)]
                results[index] = (localize_fixed_template(end, blocks) if templates[index].fixed else
                                  localize_insertion(start, end, blocks, flank, require_unique=True))
    return [annotate(template, results[index], backbone, source_forward=not template.fixed)
            for index, template in enumerate(templates)]


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--input', required=True)
    parser.add_argument('--reference-fasta', required=True)
    parser.add_argument('--reference-haplotype', required=True)
    parser.add_argument('--assemblies', help='NAME FASTA [FAI] source assembly list')
    parser.add_argument('--source', action='append', nargs=2, default=[], metavar=('SAMPLE', 'FASTA'))
    parser.add_argument('--output', required=True)
    parser.add_argument('--threads', type=int, default=1)
    parser.add_argument('--flank', type=int, default=10000)
    args = parser.parse_args(argv)
    if args.threads < 1 or args.flank < 1:
        parser.error('threads and flank must be positive')
    sources = read_query_list(args.assemblies) if args.assemblies else {}
    sources.update({sample: (fasta, None) for sample, fasta in args.source})
    sources[args.reference_haplotype] = args.reference_fasta, None
    output = Path(args.output).resolve()
    output.parent.mkdir(parents=True, exist_ok=True)
    logging.basicConfig(level=logging.INFO, format='[lift-local-templates] %(message)s')
    with tempfile.TemporaryDirectory(prefix='template-lift-', dir=output.parent) as work:
        templates = lift_templates(read_templates(args.input, load_sequences=True), sources,
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
