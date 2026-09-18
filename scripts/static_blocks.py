"""Keep each partition's original template sequences as its complete path set."""
from dataclasses import replace
import json
from pathlib import Path

from fixed_alternatives import (
    FIXED_ROLE, FIXED_TAG, STATIC_ROLE, fasta_records, template_source,
)


def write_static_paths(template_fasta, output_dir):
    """Write the ordinary adjusted-manifest format without discovering sequences."""
    from build_local_graphs import public_path_id_from_values
    from refine_partition_paths import PathRow, write_adjusted_path_outputs

    template_fasta = Path(template_fasta).resolve()
    partition = template_fasta.parent.name
    originals, references, names, records = [], [], set(), set()
    for fields, length, _digest in fasta_records(template_fasta):
        sample, contig, start, end, strand = template_source(fields)
        if start < 0 or end - start != length or length == 0 or strand not in {'+', '-'}:
            raise ValueError(f'{template_fasta}: invalid template interval for {fields[0]}')
        name = public_path_id_from_values(partition, contig, start, end)
        if name in names or fields[0] in records:
            raise ValueError(f'{template_fasta}: duplicate template path {name}')
        names.add(name)
        records.add(fields[0])
        kind = FIXED_ROLE if FIXED_TAG in fields else STATIC_ROLE
        segment = dict(fasta=str(template_fasta), record=fields[0], start=0,
                       end=length, strand='+', source_haplotype=sample,
                       source_contig=contig, source_start=start, source_end=end,
                       source_strand=strand)
        originals.append(PathRow(partition, 0, '', 'original', 'original',
            'original_template', sample, contig, start, end, strand, '.', fields[0]))
        references.append(PathRow(partition, 0, '', 'reference', 'reference',
            kind, sample, contig, start, end, strand,
            json.dumps([segment], separators=(',', ':')), fields[0]))
    if not originals:
        raise ValueError(f'{template_fasta}: no block templates')
    rows = [replace(row, path_order=order, path_id=f'{partition}_{row.role}_{order:06d}')
            for order, row in enumerate(originals + references, 1)]
    write_adjusted_path_outputs(rows, str(output_dir), jobs=1)
    return Path(output_dir) / partition / f'{partition}.adjusted.tsv'


def validate_static_paths(template_fasta, manifest, graph_fasta=None):
    """Reject added, expanded, or changed paths when resuming a static partition."""
    from build_local_graphs import public_path_id_from_values, read_manifest, public_path_id

    template_fasta = Path(template_fasta).resolve()
    partition = template_fasta.parent.name
    expected = {}
    for fields, length, digest in fasta_records(template_fasta):
        source = template_source(fields)
        name = public_path_id_from_values(partition, source[1], source[2], source[3])
        if name in expected:
            raise ValueError(f'{template_fasta}: duplicate template path {name}')
        expected[name] = (source, length, digest, fields[0])
    observed = set()
    for row in read_manifest(str(manifest), partition):
        if row['role'] == 'original':
            continue
        name = public_path_id(row)
        if name not in expected or name in observed or row['role'] != 'reference':
            raise ValueError(f'{manifest}: unexpected static path {name}')
        source, length, _digest, record = expected[name]
        actual = (row['source_haplotype'], row['source_contig'],
                  int(row['source_start']), int(row['source_end']), row['strand'])
        segments = json.loads(row['segments'])
        if (actual != source or row['source_kind'] not in {STATIC_ROLE, FIXED_ROLE}
                or len(segments) != 1 or segments[0]['record'] != record
                or Path(segments[0]['fasta']).resolve() != template_fasta
                or (segments[0]['start'], segments[0]['end'], segments[0]['strand']) != (0, length, '+')):
            raise ValueError(f'{manifest}: changed static template {name}')
        observed.add(name)
    if not expected or observed != set(expected):
        raise ValueError(f'{manifest}: missing static templates')
    if graph_fasta is not None:
        observed = set()
        for fields, length, digest in fasta_records(graph_fasta):
            name = fields[0]
            if name not in expected or name in observed:
                raise ValueError(f'{graph_fasta}: unexpected static graph path {name}')
            if (template_source(fields), length, digest) != expected[name][:3]:
                raise ValueError(f'{graph_fasta}: changed static graph sequence {name}')
            observed.add(name)
        if observed != set(expected):
            raise ValueError(f'{graph_fasta}: missing static graph templates')
