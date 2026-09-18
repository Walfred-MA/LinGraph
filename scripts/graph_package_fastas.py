"""Export final partition FASTAs without reopening their source assemblies."""
from __future__ import annotations

import csv
import hashlib
import os
from pathlib import Path
import re
import shutil
import tempfile


def select_manifests(manifests, listing):
    """Preserve the finalized list's order and never widen its partition set."""
    available = dict(manifests)
    selected = []
    seen = set()
    with open(listing) as handle:
        for number, raw in enumerate(handle, 1):
            if not raw.strip() or raw.lstrip().startswith('#'):
                continue
            fields = raw.split()
            if len(fields) != 1:
                raise ValueError(f'{listing}:{number}: expected one partition name or FASTA path')
            name = fields[0].rstrip('/').rsplit('/', 1)[-1]
            for suffix in ('.FA', '.fasta', '.fa'):
                if name.endswith(suffix):
                    name = name[:-len(suffix)]
                    break
            if not re.fullmatch(r'[A-Za-z0-9_.#-]+', name) or name in {'.', '..'}:
                raise ValueError(f'{listing}:{number}: invalid partition {name!r}')
            if name in seen or name not in available:
                raise ValueError(f'{listing}:{number}: duplicate partition or missing adjusted manifest: {name}')
            seen.add(name)
            selected.append((name, available[name]))
    if not selected:
        raise ValueError(f'{listing}: empty partition list')
    return selected


def _records(path):
    """Read one record at a time; preserve graph headers and sequence wrapping."""
    header = None
    lines = []
    with open(path) as handle:
        for number, raw in enumerate(handle, 1):
            if raw.startswith('>'):
                if header is not None:
                    yield header, ''.join(line.strip() for line in lines), ''.join(lines)
                header = raw[1:].rstrip('\r\n')
                if not header.split():
                    raise ValueError(f'{path}:{number}: empty FASTA header')
                lines = []
            elif raw.strip():
                if header is None:
                    raise ValueError(f'{path}:{number}: sequence before FASTA header')
                lines.append(raw)
        if header is not None:
            yield header, ''.join(line.strip() for line in lines), ''.join(lines)


def _source_key(pack, source):
    haplotype, contig, start, end, strand = source
    return haplotype, pack.CHR_TO_NC.get(contig, contig), int(start), int(end), strand


def _row_source(pack, row):
    return _source_key(pack, (row['source_haplotype'], row['source_contig'],
                             row['source_start'], row['source_end'], row['strand']))


def _graph_rows(rows):
    refined = any(row['role'] == 'reference' for row in rows)
    return [row for row in rows if row['role'] != 'original' or not refined]


def _copy_partition(pack, partition, manifest, folder, destination, global_ids):
    from fixed_alternatives import template_source
    rows = pack.read_manifest(manifest, partition)
    summaries = {row['path_order']: pack._package_summary_row(partition, row) for row in rows}
    expected = {}
    for row in _graph_rows(rows):
        record_id = pack.public_path_id(row)
        if record_id in expected:
            raise ValueError(f'{partition}: duplicate graph path ID {record_id}')
        expected[record_id] = row
    originals = {pack.public_path_id(row): row for row in rows if row['role'] == 'original'}
    graph = folder / partition / (partition + '.FA')
    stored = {}
    count = bases = 0

    def original_header(row, flags=()):
        provenance = (f"{row['source_haplotype']}:{row['source_contig']}:"
                      f"{row['source_start']}-{row['source_end']}{row['strand']}")
        return ' '.join([pack.public_path_id(row), provenance, '.', '.', '.', 'reference', *flags])

    def append(row, header, sequence, body):
        nonlocal count, bases
        record_id = pack.public_path_id(row)
        if record_id in global_ids:
            raise ValueError(f'{graph}: duplicate FASTA ID {record_id}')
        if not sequence or re.search(r'\s', sequence):
            raise ValueError(f'{graph}: empty sequence or whitespace inside {record_id}')
        global_ids.add(record_id)
        destination.write('>' + header + '\n' + body)
        if not body.endswith('\n'):
            destination.write('\n')
        count += 1
        bases += len(sequence)

    for header, sequence, body in _records(graph):
        record_id = header.split()[0]
        if record_id not in expected:
            raise ValueError(f'{graph}: path absent from refined manifest: {record_id}')
        row = expected[record_id]
        source = _source_key(pack, template_source(header.split()))
        if source != _row_source(pack, row):
            raise ValueError(f'{graph}: source coordinates/strand disagree with manifest for {record_id}')
        if len(sequence) != int(row['source_end']) - int(row['source_start']):
            raise ValueError(f'{graph}: sequence length disagrees with manifest for {record_id}')
        append(row, header, sequence, body)
        if record_id in originals:
            original = originals[record_id]
            original_source = _row_source(pack, original)
            if source[:4] != original_source[:4]:
                raise ValueError(f'{graph}: original and graph path ID {record_id} refer to different sources')
            oriented = sequence if source[4] == original_source[4] else pack.reverse_complement(sequence)
            stored[record_id] = hashlib.sha256(oriented.encode('ascii')).digest()
        else:
            stored[record_id] = None
    missing = expected.keys() - stored.keys()
    if missing:
        raise ValueError(f'{graph}: missing refined paths: {", ".join(sorted(missing)[:10])}')

    # Original block FASTA IDs are not necessarily the public graph IDs.
    # Match provenance and rewrite only additional original-template headers.
    templates = folder / partition / (partition + '.fasta')
    by_source = {_row_source(pack, row)[:4]: (record_id, row)
                 for record_id, row in originals.items()}
    found = set()
    if templates.is_file():
        for header, sequence, _body in _records(templates):
            source = _source_key(pack, template_source(header.split()))
            if source[:4] not in by_source:
                raise ValueError(f'{templates}: template absent from manifest: {header}')
            record_id, row = by_source[source[:4]]
            if record_id in found:
                raise ValueError(f'{templates}: duplicate original template {record_id}')
            found.add(record_id)
            if len(sequence) != source[3] - source[2]:
                raise ValueError(f'{templates}: sequence length disagrees with coordinates for {record_id}')
            if source[4] != row['strand']:
                sequence = pack.reverse_complement(sequence)
            if record_id in stored:
                if hashlib.sha256(sequence.encode('ascii')).digest() != stored[record_id]:
                    raise ValueError(f'{templates}: conflicting original/graph sequences for shared ID {record_id}')
                continue
            flags = [field for field in header.split()[1:] if field.startswith('sequence_role=')]
            description = original_header(row, flags)
            append(row, description, sequence, pack.wrap_fasta(sequence) + '\n')
            stored[record_id] = hashlib.sha256(sequence.encode('ascii')).digest()
    missing = originals.keys() - stored.keys()
    if missing:
        raise ValueError(f'{templates}: missing original hotspot templates: {", ".join(sorted(missing)[:10])}')
    return [summaries[row['path_order']] for row in rows], count, bases


def write_graph_package(manifests, partition_sources, output_dir, jobs,
                        query_paths=None, reference_haplotype=None,
                        embed_all_templates=False, *, graph_folder):
    """Use retained final FASTAs; extract only partitions whose files are gone."""
    import build_local_graphs as pack
    folder = Path(graph_folder)
    ready = set()
    fallback = []
    for partition, manifest in manifests:
        graph = folder / partition / (partition + '.FA')
        template = folder / partition / (partition + '.fasta')
        rows = pack.read_manifest(manifest, partition)
        graph_ids = {pack.public_path_id(row) for row in _graph_rows(rows)}
        extra_originals = any(row['role'] == 'original' and pack.public_path_id(row) not in graph_ids
                              for row in rows)
        if (graph.is_file() and graph.stat().st_size
                and (not extra_originals or (template.is_file() and template.stat().st_size))):
            ready.add(partition)
        else:
            fallback.append((partition, manifest))
    if not ready:
        pack.LOG.info('No complete retained partition FASTAs; using existing sequence extraction')
        return pack.write_package_indexes(manifests, partition_sources, output_dir, jobs,
                                          query_paths, reference_haplotype, embed_all_templates)
    pack.LOG.info('Exporting %d partitions directly from final FASTAs; %d require sequence extraction',
                  len(ready), len(fallback))
    workdir = tempfile.mkdtemp(prefix='.build_local_graphs.package.work.', dir=output_dir)
    success = False
    try:
        fallback_rows = {}
        fallback_ranges = {}
        if fallback:
            fallback_dir = Path(workdir) / 'extracted'
            fallback_dir.mkdir()
            pack.write_package_indexes(fallback,
                                       {name: partition_sources[name] for name, _ in fallback},
                                       str(fallback_dir), jobs, query_paths,
                                       reference_haplotype, embed_all_templates)
            with (fallback_dir / 'local_graphs.tsv').open() as handle:
                for row in csv.DictReader(handle, delimiter='\t'):
                    fallback_rows.setdefault(row['partition'], []).append(row)
            fallback_ranges = pack.index_fasta_byte_ranges(str(fallback_dir / 'alternatives.fasta'))
        summary_path = Path(workdir) / 'local_graphs.tsv'
        sequences_path = Path(workdir) / 'alternatives.fasta'
        path_count = stored_count = bases = slice_count = 0
        global_ids = set()
        with summary_path.open('w', newline='') as summary, sequences_path.open('w') as sequences:
            writer = csv.DictWriter(summary, fieldnames=pack.SUMMARY_FIELDS, delimiter='\t', lineterminator='\n')
            writer.writeheader()
            for index, (partition, manifest) in enumerate(manifests, 1):
                if partition in ready:
                    rows, count, size = _copy_partition(pack, partition, manifest, folder, sequences, global_ids)
                    stored_count += count
                    bases += size
                else:
                    rows = fallback_rows[partition]
                    copied = set()
                    with (fallback_dir / 'alternatives.fasta').open('rb') as source:
                        for row in rows:
                            record_id = pack.public_path_id(row)
                            if record_id not in fallback_ranges or record_id in copied:
                                continue
                            if record_id in global_ids:
                                raise ValueError(f'duplicate FASTA ID across partitions: {record_id}')
                            global_ids.add(record_id)
                            copied.add(record_id)
                            start, end = fallback_ranges[record_id]
                            source.seek(start)
                            record = source.read(end - start).decode('utf-8')
                            sequences.write(record)
                            if not record.endswith('\n'):
                                sequences.write('\n')
                            stored_count += 1
                            bases += int(row['source_end']) - int(row['source_start'])
                writer.writerows(rows)
                path_count += len(rows)
                slice_count += sum(row['reference_slice'] != '.' for row in rows)
                if index % max(1, len(manifests) // 20) == 0 or index == len(manifests):
                    pack.LOG.info('Packaged %d/%d partitions (%d stored paths)', index, len(manifests), stored_count)
        pack.validate_embedded_templates(workdir, reference_haplotype)
        os.replace(summary_path, Path(output_dir) / 'local_graphs.tsv')
        os.replace(sequences_path, Path(output_dir) / 'alternatives.fasta')
        success = True
        return path_count, stored_count + slice_count, stored_count, bases, slice_count
    finally:
        if success:
            shutil.rmtree(workdir)
        else:
            pack.LOG.error('Retained failed package work directory: %s', workdir)
