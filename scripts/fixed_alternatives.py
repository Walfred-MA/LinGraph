"""Imported alternatives are immutable templates; source labels are provenance."""
import json
from pathlib import Path
import re
import csv
import hashlib

FIXED_ROLE = 'imported_alternative'
STATIC_ROLE = 'static_template'
STORED_TEMPLATE_ROLES = frozenset({FIXED_ROLE, STATIC_ROLE})
FIXED_TAG = 'sequence_role=' + FIXED_ROLE
POLICY = 'fixed-imported-alternatives-v1'
VALIDATION_POLICY = 'fixed-imported-alternatives-validated-v2'
ANNOTATION_POLICY = 'fixed-imported-record-coordinates-v1'
ANNOTATION_HAPLOTYPE = 'alternative'


def annotation_reference_records(path):
    """Address imported sequences by their catalog record, not source provenance.

    An imported record can differ from its source assembly or have no source
    assembly available. Its stored sequence and forward record coordinates
    remain authoritative for every annotation stage.
    """
    from dataclasses import replace
    from refine_partition_paths import parse_main_reference
    from summarize_partition_novelty import fasta_header_descriptions
    records = parse_main_reference(path)
    for _, description in fasta_header_descriptions(path):
        name = description.split()[0]
        record = records[name]
        if is_fixed(description) or record.haplotype == ANNOTATION_HAPLOTYPE:
            records[name] = replace(
                record, haplotype=ANNOTATION_HAPLOTYPE,
                source_contig=name, source_start=0,
                source_end=record.length, strand='+')
    return records


def annotation_reference_sources(assembly_sources, catalog):
    """Add the fixed catalog as a sequence source without adding a query sample."""
    from map_partition_local_novel import AssemblyFastaSource
    path = str(Path(catalog).resolve())
    sources = dict(assembly_sources)
    if ANNOTATION_HAPLOTYPE in sources:
        raise ValueError(f'{ANNOTATION_HAPLOTYPE!r} is reserved for imported annotation targets')
    sources[ANNOTATION_HAPLOTYPE] = AssemblyFastaSource(
        ANNOTATION_HAPLOTYPE, path, path + '.fai')
    return sources


def attributes(description):
    return dict(token.split('=', 1) for token in description.split()[1:] if '=' in token)


def is_fixed(description):
    return attributes(description).get('sequence_role') == FIXED_ROLE


def has_fixed_templates(path):
    with open(path) as source:
        return any(line.startswith('>') and is_fixed(line[1:]) for line in source)


def template_source(fields):
    attrs = attributes(' '.join(fields))
    match = re.fullmatch(r'(.+):(\d+)-(\d+)', fields[1]) if len(fields) > 1 else None
    if match is not None and 'sample' in attrs:
        contig, start, end = match.groups()
        return attrs['sample'], contig, int(start), int(end), attrs.get('strand', '+')
    for field in fields[1:]:
        match = re.fullmatch(r'([^:]+):(.+):(\d+)-(\d+)([+-]?)', field.removeprefix('source='))
        if match:
            sample, contig, start, end, strand = match.groups()
            return sample, contig, int(start), int(end), strand or '+'
    raise ValueError(f'{fields[0]}: malformed fixed-template provenance')


def partition_templates(partition, files):
    """Index explicitly fixed records in the persistent block-template FASTA."""
    path = Path(files.reference_bed).parent / (partition + '.fasta')
    result = {}
    if not path.is_file():
        return result
    with path.open() as source:
        for line in source:
            if not line.startswith('>') or not is_fixed(line[1:]):
                continue
            fields = line[1:].split()
            sample, contig, start, end, strand = template_source(fields)
            result[fields[0]] = dict(
                fasta=str(path.resolve()), record=fields[0], start=0,
                end=int(end)-int(start), strand='+',
                source_haplotype=sample, source_contig=contig,
                source_start=int(start), source_end=int(end),
                source_strand=strand)
    return result


def fixed_partitions(partitions):
    return {name: records for name, files in partitions.items()
            if (records := partition_templates(name, files))}


def unbacked_template_haplotypes(partitions, used_rows, observed_haplotypes):
    """Accept absent source assemblies only for verified imported BED records."""
    fixed = fixed_partitions(partitions)
    missing = set()
    for row in used_rows:
        if row.haplotype in observed_haplotypes:
            continue
        original = fixed.get(row.partition, {}).get(row.encoded_name)
        if (original is None
                or original['source_haplotype'] != row.haplotype
                or original['source_contig'] != row.source_contig
                or not original['source_start'] <= row.source_start < row.source_end <= original['source_end']):
            missing.add(row.haplotype)
    return missing


def fixed_segment(row):
    segments = json.loads(row['segments'])
    if len(segments) != 1 or not isinstance(segments[0], dict):
        raise ValueError('fixed alternative requires one original FASTA segment')
    return segments[0]


def materialize_fixed(row, partition, reader=None):
    """Read an imported or static block record from its stored template FASTA."""
    from minsetref_core import IndexedFasta
    segment = fixed_segment(row)
    if reader is None:
        with IndexedFasta(segment['fasta']) as indexed:
            return materialize_fixed(row, partition, indexed)
    sequence = reader.fetch(segment['record'], segment['start'], segment['end'])
    if len(sequence) != int(row['source_end'])-int(row['source_start']):
        raise ValueError(f'{partition}: fixed alternative sequence length changed')
    from build_local_graphs import public_path_id
    coord = f"{row['source_contig']}:{row['source_start']}-{row['source_end']}{row['strand']}"
    tag = 'sequence_role=' + row['source_kind']
    header = (f">{public_path_id(row)} {row['source_haplotype']}:{coord} "
              f"{coord} . {coord} {tag} reference")
    return header, sequence


def fasta_records(path):
    """Stream header, length and case-sensitive SHA256, without loading bases."""
    fields = None
    length = 0
    digest = None
    with open(path, 'rb') as source:
        for raw in source:
            if raw.startswith(b'>'):
                if fields is not None:
                    yield fields, length, digest.hexdigest()
                fields = raw[1:].decode().split()
                if not fields:
                    raise ValueError(f'{path}: empty FASTA header')
                length, digest = 0, hashlib.sha256()
            elif raw.strip() and not raw.startswith(b'#'):
                if fields is None:
                    raise ValueError(f'{path}: sequence before FASTA header')
                sequence = raw.rstrip(b'\r\n')
                length += len(sequence)
                digest.update(sequence)
        if fields is not None:
            yield fields, length, digest.hexdigest()


def validate_fixed_paths(template_fasta, manifest, graph_fasta=None, record_ids=None):
    """Verify imported path identity, bounds, orientation and exact sequence.

    record_ids allows a migration audit of originals whose headers predate the
    explicit fixed tag. Callers must obtain those IDs by checking the import
    catalog, not by trusting the graph completion marker.
    """
    from build_local_graphs import public_path_id
    from minsetref_core import IndexedFasta
    selected = None if record_ids is None else set(record_ids)
    expected = {}
    for fields, length, digest in fasta_records(template_fasta):
        included = fields[0] in selected if selected is not None else FIXED_TAG in fields
        if not included:
            continue
        if fields[0] in expected:
            raise ValueError(f'{template_fasta}: duplicate imported record {fields[0]}')
        source = template_source(fields)
        if source[3]-source[2] != length or source[4] not in {'+','-'}:
            raise ValueError(f'{fields[0]}: original template bounds/orientation do not match FASTA')
        expected[fields[0]] = (source, length, digest)
    if selected is not None and set(expected) != selected:
        raise ValueError(f'{template_fasta}: imported records missing from original FASTA')
    found, graph_expected = set(), {}
    with open(manifest, newline='') as source:
        for row in csv.DictReader(source, delimiter='\t'):
            if row.get('source_kind') != FIXED_ROLE:
                continue
            segment = fixed_segment(row)
            record = segment['record']
            if record not in expected or record in found:
                raise ValueError(f'{manifest}: unexpected/duplicate imported template {record}')
            found.add(record)
            provenance, length, digest = expected[record]
            actual = (row['source_haplotype'], row['source_contig'],
                      int(row['source_start']), int(row['source_end']), row['strand'])
            segment_source = (segment['source_haplotype'],segment['source_contig'],
                              int(segment['source_start']),int(segment['source_end']),segment['source_strand'])
            if (row['role'] != 'reference' or actual != provenance or segment_source != provenance
                    or int(segment['start']) != 0 or int(segment['end']) != length
                    or segment['strand'] != '+'):
                raise ValueError(f'{record}: imported template interval/orientation changed')
            # A staged repair can point back to the original live FASTA. Verify
            # the bytes rather than requiring identical filesystem paths.
            observed = hashlib.sha256()
            with IndexedFasta(segment['fasta']) as reader:
                if reader.length(record) != length:
                    raise ValueError(f'{record}: imported source FASTA length changed')
                for start in range(0,length,1024*1024):
                    observed.update(reader.fetch(record,start,min(start+1024*1024,length)).encode())
            if observed.hexdigest() != digest:
                raise ValueError(f'{record}: imported source FASTA sequence changed')
            name = public_path_id(row)
            if name in graph_expected:
                raise ValueError(f'{manifest}: duplicate imported public path {name}')
            graph_expected[name] = (provenance,length,digest)
    if found != set(expected):
        raise ValueError(f'{manifest}: missing fixed imported templates: {sorted(set(expected)-found)}')
    if graph_fasta is not None:
        seen = set()
        for fields,length,digest in fasta_records(graph_fasta):
            name = fields[0]
            if FIXED_TAG not in fields and name not in graph_expected:
                continue
            if name in seen or name not in graph_expected:
                raise ValueError(f'{graph_fasta}: unexpected/duplicate imported graph path {name}')
            seen.add(name)
            if (FIXED_TAG not in fields or fields[-1] != 'reference'
                    or (template_source(fields),length,digest) != graph_expected[name]):
                raise ValueError(f'{name}: imported graph path identity, coordinates or sequence changed')
        if seen != set(graph_expected):
            raise ValueError(f'{graph_fasta}: missing imported graph paths: {sorted(set(graph_expected)-seen)}')
    return {'templates':len(expected),'bases':sum(row[1] for row in expected.values()),
            'policy':VALIDATION_POLICY}
