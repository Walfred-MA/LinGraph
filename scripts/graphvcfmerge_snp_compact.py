"""Numeric SNP spools: one observation is 21 bytes, with interned provenance.

CHROM is stored in section metadata outside the records. Bases use two IUPAC
nibbles. Scalar LABEL_H uses signed zigzag encoding; joined cross-graph values
use an interned dictionary index in the same uint32 field. Its encoding mode,
strand and SNP type live in the query dictionary. No Python SNP objects or VCF
strings are retained during chromosome sorting.
"""
from __future__ import annotations

from bisect import bisect_right
from collections import defaultdict
import hashlib
import json
import multiprocessing
import os
from pathlib import Path
import re
import struct
import sys
import tempfile
import time

import numpy as np
import graphvcfmerge as vcf

DTYPE = np.dtype([('pos', '<u4'), ('bases', 'u1'), ('query_pos', '<u4'),
                  ('query', '<u4'), ('allele', '<u4'), ('label_h', '<u4')], align=False)
PACK = struct.Struct('<IBIIII')
BASES = 'ACGTRYSWKMBDHVN'
BASE_INDEX = {base: index for index, base in enumerate(BASES)}
PROTOCOL = 4
LABEL_GROUP = 2  # Query metadata: False = missing, True = scalar, 2 = joined.


def chrom_key(chrom):
    return hashlib.sha256(chrom.encode()).hexdigest()


def u32(value):
    value = int(value)
    if not 0 <= value < 2**32:
        raise ValueError(f'SNP coordinate/index outside uint32: {value}')
    return value


def pack_bases(ref, alt):
    try:
        return (BASE_INDEX[ref.upper()] << 4) | BASE_INDEX[alt.upper()]
    except KeyError as error:
        raise ValueError(f'SNP requires single IUPAC REF/ALT bases: {ref!r}/{alt!r}') from error


def label_number(text):
    if text == '.':
        return 0, False
    if not re.fullmatch(r'[+-]?\d+H', text):
        raise ValueError(f'invalid SNP LABEL_H: {text!r}')
    value = int(text[:-1])
    if not -(2**31) <= value < 2**31:
        raise ValueError(f'SNP LABEL_H outside signed 32-bit encoding: {text}')
    return (value << 1) ^ (value >> 31), True


def encode_label(text, groups):
    if '&' not in text:
        return label_number(text)
    if text not in groups:
        # Contributors are independent of ALLELENAME cardinality. Preserve the
        # entire joined value, including order and explicit missing components.
        for component in text.split('&'):
            label_number(component)
        groups[text] = u32(len(groups))
    return groups[text], LABEL_GROUP


def label_text(value, mode=True, groups=()):
    if mode == LABEL_GROUP:
        return groups[int(value)]
    if not mode:
        return '.'
    return f'{(int(value) >> 1) ^ -(int(value) & 1)}H'


def snp_observations(field, fmt):
    """Decode split FORMAT columns without an ambiguous colon-joined payload."""
    if field.split(':', 1)[0] in ('', '.', './.', '0', '0/0'):
        return []
    if fmt == 'HSNP':
        result = []
        for allele in vcf.split_sample_alleles(field, fmt):
            values = vcf.parse_hsv_allele(allele)
            if values is None:
                raise ValueError('malformed SNP observation')
            result.append([vcf.vcf_unescape(value) if value != '.' else '.' for value in values])
        return result
    columns = field.split(':')
    names = fmt.split(':')
    if len(columns) != len(names):
        raise ValueError('malformed allele: SNP FORMAT column count')
    vectors = [column.split(',') for column in columns[1:]]
    if not vectors or any(len(vector) != len(vectors[0]) for vector in vectors):
        raise ValueError('malformed allele: SNP FORMAT observation count')
    result = []
    for index in range(len(vectors[0])):
        values = dict(zip(names[1:], (vector[index] for vector in vectors)))
        values['GT'] = columns[0]
        result.append([vcf.vcf_unescape(values.get(name, '.')) if values.get(name, '.') != '.' else '.'
                       for name in vcf.SNP_EVENT_FIELDS])
    return result


class Writer:
    """One binary file per scan worker; chromosome spans live in metadata."""
    def __init__(self, directory):
        self.directory = Path(directory)
        self.directory.mkdir(parents=True, exist_ok=True)
        self.queries, self.alleles, self.label_groups = {}, {}, {}
        self.metadata = {}
        self.chroms, self.coverage, self.lengths = {}, defaultdict(list), {}
        self.handle = (self.directory / "records.bin").open("wb")
        self.sections = defaultdict(list)
        self.counts = defaultdict(int)
        self.last_key = None
        self.offset = 0

    def add(self, chrom, pos, ref, alt, sample, values=None, state='1', filt='PASS'):
        if int(pos) < 1:
            raise ValueError('SNP VCF position must be positive')
        if values is None:
            query, qpos, strand, kind, allele, label, label_mode = '.', 0, '+', 'SNP', '.', 0, False
        else:
            if len(values) < 8 or values[0] != '1' or 'SNP' not in values[1] or values[2] != '1':
                raise ValueError('malformed SNP observation')
            if values[3].upper() != alt.upper():
                raise ValueError('SNP observation/base disagreement')
            match = re.fullmatch(r'(\d+)([+-])', values[5])
            if not match:
                raise ValueError(f'invalid SNP QUERYCOORD: {values[5]!r}')
            qpos, strand = int(match[1]), match[2]
            query, kind, allele = values[4], values[1], values[6]
            label, label_mode = encode_label(values[7], self.label_groups)
        query_key = (sample, query, strand, kind, state, label_mode, filt)
        qi = self.queries.setdefault(query_key, len(self.queries))
        ai = self.alleles.setdefault(allele, len(self.alleles))
        packed = PACK.pack(u32(pos), pack_bases(ref, alt), u32(qpos), u32(qi), u32(ai), u32(label))
        key = chrom_key(chrom)
        self.chroms[key] = chrom
        if key != self.last_key:
            self.sections[key].append([self.offset, 0])
            self.last_key = key
        self.sections[key][-1][1] += 1
        self.counts[key] += 1
        self.offset += PACK.size
        self.handle.write(packed)

    def consume(self, row, samples):
        if len(row.samples) != len(samples):
            raise ValueError('sample column count differs from header')
        if not vcf.is_snp_format(row.fmt):
            return
        for sample, field in zip(samples, row.samples):
            values = snp_observations(field, row.fmt)
            state = field.split(':', 1)[0]
            if not values:
                if state not in ('', '.', './.', '0', '0/0'):
                    raise ValueError('malformed allele')
                self.add(row.chrom, row.pos, row.ref, row.alt, sample,
                         state='0' if state in ('0', '0/0') else '.', filt=row.filt)
            for parsed in values:
                self.add(row.chrom, row.pos, row.ref, row.alt, sample, parsed, filt=row.filt)

    def close(self):
        self.handle.close()

    def finish(self):
        self.close()
        data = dict(protocol=PROTOCOL, directory=str(self.directory.resolve()),
                    queries=list(self.queries), alleles=list(self.alleles),
                    label_groups=list(self.label_groups), chroms=self.chroms,
                    coverage=dict(self.coverage), lengths=self.lengths,
                    counts=dict(self.counts), sections=dict(self.sections), metadata=list(self.metadata))
        vcf.checkpoints.write_json(self.directory / 'source.json', data)
        return data


def _scan_one(task):
    index, paths, directory = task
    writer = Writer(Path(directory) / str(index))
    all_samples, metadata = [], []
    try:
        for path in paths:
            samples = []
            with vcf.open_text(str(path)) as source:
                for number, raw in enumerate(source, 1):
                    try:
                        if raw.startswith('##referenceCoverage=<'):
                            data = vcf._parse_structured_meta(raw.strip(), 'referenceCoverage')
                            if data.get('Sample'):
                                chrom = data['Chrom']
                                writer.chroms[chrom_key(chrom)] = chrom
                                writer.coverage[chrom].append([data['Sample'], int(data['Start']), int(data['End'])])
                        elif raw.startswith('##'):
                            metadata.append(raw.strip())
                        elif raw.startswith('#CHROM\t'):
                            samples = raw.rstrip('\r\n').split('\t')[9:]
                            if not samples or len(samples) != len(set(samples)):
                                raise ValueError('expected unique sample names in the VCF header')
                        elif raw.strip() and not raw.startswith('#'):
                            if not samples:
                                raise ValueError('missing #CHROM sample header')
                            row = vcf.parse_vcf_record(raw, number)
                            writer.consume(row, samples)
                    except (ValueError, KeyError) as error:
                        raise ValueError(f'{path}:{number}: {error}') from error
            if not samples:
                raise ValueError(f'{path}: missing #CHROM sample header')
            all_samples.extend(sample for sample in samples if sample not in all_samples)
        source = writer.finish()
        return dict(source=str(writer.directory / 'source.json'), samples=all_samples,
                    metadata=metadata, chroms=source['chroms'])
    finally:
        writer.close()


def scan(paths, root, processes=1, manifest_output=None, insertion_snps=None):
    if not paths or processes < 1:
        raise ValueError('SNP scan needs inputs and positive --processes')
    root = Path(os.path.abspath(root))
    root.mkdir(parents=True, exist_ok=True)
    directory = tempfile.mkdtemp(prefix='scan-', dir=root)
    tasks = [(i, group, directory) for i, group in enumerate(worker_groups(paths, processes))]
    if processes > 1 and len(tasks) > 1:
        with multiprocessing.get_context('spawn').Pool(min(processes, len(tasks))) as pool:
            results = pool.map(_scan_one, tasks)
    else:
        results = [_scan_one(task) for task in tasks]
    return manifest_from_sources(results, root, directory, manifest_output, insertion_snps)


def worker_groups(paths, processes):
    # Contiguous batches preserve input/sample order independently of worker count.
    count = min(max(1, processes), len(paths))
    return [[str(p) for p in paths[len(paths)*i//count:len(paths)*(i+1)//count]] for i in range(count)]


def prepare(root, raw_manifest, manifest_output=None, insertion_snps=None):
    raw = json.loads(Path(raw_manifest).read_text())
    result = dict(source=None, sources=raw['sources'], samples=raw['samples'],
                  metadata=raw['metadata'], chroms=raw['loci'])
    return manifest_from_sources([result], Path(root), raw['directory'], manifest_output,
                                 insertion_snps or raw.get('insertion_snps'))


def manifest_from_sources(results, root, directory, manifest_output=None, insertion_snps=None):
    root = Path(root)
    root.mkdir(parents=True, exist_ok=True)
    samples, metadata, chroms, sources, definitions = [], [], {}, [], {}
    for result in results:
        samples.extend(sample for sample in result['samples'] if sample not in samples)
        chroms.update(result['chroms'])
        sources.extend(result['sources'] if 'sources' in result else [result['source']])
        for line in result['metadata']:
            if line.startswith('##contig=<'):
                data = vcf._parse_structured_meta(line, 'contig')
                name, length = data.get('ID'), data.get('length')
                if name in definitions:
                    if definitions[name] != length:
                        raise ValueError(f'conflicting contig lengths for {name!r}')
                    continue
                definitions[name] = length
            if line not in metadata:
                metadata.append(line)
    order = vcf._contig_order(metadata)
    chroms = dict(sorted(chroms.items(), key=lambda item: (order.get(item[1], len(order)), item[1])))
    by_chrom = defaultdict(list)
    for source in sources:
        for key in json.loads(Path(source).read_text())['chroms']:
            by_chrom[key].append(source)
    # Newly discovered insertion loci are small and independent. Keep them out
    # of chromosome jobs; concat opens their manifest only after SVs finish.
    jobs = dict(chroms)
    job_loci = {key: [key] for key in jobs}
    locus_order = list(jobs)
    # Declare contigs in the same order as the concatenated chromosome parts.
    contigs = {vcf._parse_structured_meta(line, 'contig')['ID']: line
               for line in metadata if line.startswith('##contig=<')}
    metadata = [line for line in metadata if not line.startswith('##contig=<')]
    for key in locus_order:
        line = contigs.pop(chroms[key], None)
        if line:
            metadata.append(line)
    metadata.extend(contigs.values())
    manifest = dict(protocol=PROTOCOL, directory=directory, samples=samples,
                    metadata=metadata, chroms=jobs, chrom_order=list(jobs), loci=chroms,
                    job_loci=dict(job_loci), locus_order=locus_order,
                    sources=sources, sources_by_chrom=dict(by_chrom),
                    insertion_snps=os.path.abspath(insertion_snps) if insertion_snps else None)
    vcf.checkpoints.write_json(root / 'manifest.json', manifest)
    if manifest_output:
        vcf.checkpoints.write_json(manifest_output, manifest)
    print(f'[merge:snp] saved 21-byte records: {len(sources)} shards, {len(chroms)} loci in {len(jobs)} chromosome jobs', file=sys.stderr)
    return manifest


def source_sections(data, key):
    """Validate physical spool and return byte offset/count spans (v2 readable)."""
    directory = Path(data['directory'])
    if 'sections' in data:
        binary = directory / 'records.bin'
        expected = sum(data['counts'].values()) * DTYPE.itemsize
        spans = data['sections'].get(key, [])
    else:
        binary = directory / (key + '.bin')
        expected = data['counts'].get(key, 0) * DTYPE.itemsize
        spans = [(0, expected // DTYPE.itemsize)] if expected else []
    if not binary.is_file() and expected == 0:
        return binary, spans
    if not binary.is_file() or binary.stat().st_size != expected:
        raise ValueError(f'missing or truncated SNP shard: {binary}')
    if sum(count for _, count in spans) != data['counts'].get(key, 0):
        raise ValueError(f'invalid SNP section counts: {binary}')
    for offset, count in spans:
        if offset < 0 or count < 0 or offset % DTYPE.itemsize or offset + count * DTYPE.itemsize > expected:
            raise ValueError(f'invalid SNP section bounds: {binary}')
    return binary, spans


def _load_source(task):
    # Insertion loci fit in memory. Original chromosome sorting uses bounded spans.
    path, key, locus = task
    data = json.loads(Path(path).read_text())
    binary, spans = source_sections(data, key)
    chunks = [np.fromfile(binary, dtype=DTYPE, offset=offset, count=count) for offset, count in spans]
    array = np.concatenate(chunks) if chunks else np.empty(0, dtype=DTYPE)
    return array, data, key, locus


def remap_records(data, query_map, allele_map, group_map, grouped_queries, source):
    """Translate local dictionary IDs without changing scalar LABEL_H values."""
    if len(data) and (data['query'].max() >= len(query_map) or data['allele'].max() >= len(allele_map)):
        raise ValueError(f'invalid SNP dictionary index: {source}')
    grouped = grouped_queries[data['query']]
    group_ids = data['label_h'][grouped]
    if len(group_ids) and group_ids.max() >= len(group_map):
        raise ValueError(f'invalid SNP LABEL_H dictionary index: {source}')
    data['label_h'][grouped] = group_map[group_ids]
    data['query'] = query_map[data['query']]
    data['allele'] = allele_map[data['allele']]


def concat_insertions(manifest, insertion_snps=None):
    """Plan insertion append and complete the header, without loading SNP arrays."""
    root = insertion_snps or manifest.get('insertion_snps')
    if not root:
        return manifest['metadata'], []
    index = json.loads((Path(root) / 'manifest.json').read_text())
    # Old manifests included insertion spools in chromosome jobs. Their rows
    # already exist in those parts, so explicitly passing the same root is safe.
    embedded = {str(Path(path).resolve()) for path in manifest.get('sources', [])}
    originals = set(manifest.get('loci', manifest['chroms']).values())
    samples = set(manifest['samples'])
    definitions = {}
    for line in manifest['metadata']:
        if line.startswith('##contig=<'):
            data = vcf._parse_structured_meta(line, 'contig')
            definitions[data['ID']] = data.get('length')
    by_locus, names, lengths = defaultdict(list), {}, {}
    paths = dict.fromkeys(str(Path(path).resolve()) for path in index['sources'])
    for path in paths:
        if path in embedded:
            continue
        data = json.loads(Path(path).read_text())
        if any(query[0] not in samples for query in data['queries']):
            raise ValueError(f'insertion SNP source contains an unknown sample: {path}')
        for key, name in data['chroms'].items():
            if name in originals:
                raise ValueError(f'insertion contig already occurs in chromosome parts: {name!r}')
            length = str(data['lengths'][name])
            if (name in definitions and definitions[name] != length or
                    name in lengths and lengths[name] != length):
                raise ValueError(f'conflicting insertion length for {name!r}')
            lengths[name] = length
            names[key] = name
            by_locus[key].append(path)
    plan = [(key, names[key], paths) for key, paths in sorted(by_locus.items(), key=lambda item: names[item[0]])]
    # Insertions follow every original chromosome in both the body and header.
    metadata = [line for line in manifest['metadata'] if not (
        line.startswith('##contig=<') and vcf._parse_structured_meta(line, 'contig')['ID'] in lengths)]
    metadata.extend(f'##contig=<ID={name},length={lengths[name]}>' for _, name, _ in plan)
    return metadata, plan


def append_insertions(handle, manifest, plan):
    """Load, sort and emit one insertion locus at a time; never create sort runs."""
    total = 0
    for key, name, paths in plan:
        queries, alleles, groups = {}, {}, {}
        coverage, arrays = defaultdict(list), []
        for path in paths:
            array, data, _, _ = _load_source((path, key, 0))
            query_map = np.array([queries.setdefault(tuple(q), len(queries)) for q in data['queries']], dtype='<u4')
            allele_map = np.array([alleles.setdefault(a, len(alleles)) for a in data['alleles']], dtype='<u4')
            group_map = np.array([groups.setdefault(g, len(groups)) for g in data.get('label_groups', [])], dtype='<u4')
            grouped_queries = np.array([q[5] == LABEL_GROUP for q in data['queries']], dtype=bool)
            remap_records(array, query_map, allele_map, group_map, grouped_queries, path)
            arrays.append(array)
            for sample, start, end in data['coverage'].get(name, []):
                coverage[sample].append((start, end))
        ordered = np.concatenate(arrays)
        del arrays, array
        ordered.sort(kind='mergesort', order=DTYPE.names)
        total += write_locus(handle, name, ordered, list(queries), list(alleles), list(groups),
                             coverage, manifest['samples'])
        del ordered
    if plan:
        print(f'[merge:snp] concat appended {total} sites from {len(plan)} insertion loci', file=sys.stderr)
    return total


def merge_chrom(root, chrom, processes=1):
    manifest = json.loads((Path(root) / 'manifest.json').read_text())
    key = chrom if chrom in manifest['chroms'] else chrom_key(chrom)
    if key not in manifest['chroms']:
        raise ValueError(f'unknown SNP chromosome: {chrom}')
    with vcf.checkpoints.chrom_lock(manifest['directory'], key):
        return _emit_chrom(manifest, key, processes)


def _emit_chrom(manifest, key, processes):
    from graphvcfmerge_snp import atomic_output, part_path
    from graphvcfmerge_snp_memory import sorted_locus
    chrom = manifest['chroms'][key]
    loci = manifest['job_loci'][key]
    parts = Path(manifest['directory']) / 'parts'
    parts.mkdir(exist_ok=True)
    count = 0
    with atomic_output(part_path(manifest, key, 'snp')) as handle:
        for locus in loci:
            ordered, queries, alleles, groups, coverage = sorted_locus(manifest, locus, max(1, processes))
            print(f'[merge:snp] {manifest["loci"][locus]}: writing VCF from {len(ordered):,} observations',
                  file=sys.stderr, flush=True)
            count += write_locus(handle, manifest['loci'][locus], ordered,
                                 queries, alleles, groups, coverage, manifest['samples'], progress=True)
            del ordered, queries, alleles, groups, coverage
    with atomic_output(part_path(manifest, key, 'indel')):
        pass
    print(f'[merge:snp] {chrom}: {count} sites; in-memory merge complete', file=sys.stderr, flush=True)
    return key


def write_locus(handle, chrom, ordered, query_list, allele_list, label_groups, coverage, samples, *, progress=False):
    """Emit sorted SNP observations with identical semantics for either merge stage."""
    coverage = {sample: vcf._merge_intervals(values) for sample, values in coverage.items()}

    def covered(sample, pos):
        intervals = coverage.get(sample, ())
        start = max(0, pos - 1)
        index = bisect_right(intervals, (start, float('inf'))) - 1
        return index >= 0 and intervals[index][0] <= start and pos <= intervals[index][1]

    missing_field = vcf._state_sample_field('.', vcf.SNP_FORMAT)
    reference_field = vcf._state_sample_field('0', vcf.SNP_FORMAT)
    chrom_token = vcf._safe_variant_token(chrom)
    last_report = time.monotonic()
    left = count = 0
    while left < len(ordered):
        pos, bases = int(ordered[left]['pos']), int(ordered[left]['bases'])
        right = left + 1
        while right < len(ordered) and ordered[right]['pos'] == pos and ordered[right]['bases'] == bases:
            right += 1
        ref, alt = BASES[bases >> 4], BASES[bases & 15]
        observations, states, filters = defaultdict(set), defaultdict(set), set()
        for record in ordered[left:right]:
            sample, query, strand, kind, state, label_mode, filt = query_list[int(record['query'])]
            filters.update(filt.split(';'))
            if state != '1':
                states[sample].add(state)
                continue
            observations[sample].add((
                '1', kind, '1', alt, query, f'{int(record["query_pos"])}{strand}',
                allele_list[int(record['allele'])], label_text(record['label_h'], label_mode, label_groups),
            ))
        if observations:
            fields = []
            for sample in samples:
                values = observations.get(sample)
                if values:
                    values = sorted(values)
                    fields.append('1:' + ':'.join(','.join(vcf.vcf_escape(value[index]) for value in values)
                                                for index in range(1, len(vcf.SNP_EVENT_FIELDS))))
                else:
                    sample_states = states.get(sample, ())
                    fields.append(missing_field if '.' in sample_states else reference_field
                                  if '0' in sample_states or covered(sample, pos) else missing_field)
            digest = hashlib.sha256(json.dumps([chrom, [pos, ref, alt]]).encode()).hexdigest()[:20]
            filt = 'PASS' if 'PASS' in filters else ';'.join(sorted(filters - {'.', ''})) or '.'
            handle.write('\t'.join([chrom, str(pos), f'SNP_{chrom_token}_{pos}_{digest}',
                                    ref, alt, '.', filt, f'NSUP={len(observations)}', vcf.SNP_FORMAT, *fields]) + '\n')
            count += 1
        left = right
        if progress and count % 10000 == 0 and time.monotonic() - last_report >= 10:
            print(f'[merge:snp] {chrom}: wrote {count:,} sites; {left:,}/{len(ordered):,} observations',
                  file=sys.stderr, flush=True)
            last_report = time.monotonic()
    return count


def insertion_owner(root, path_id):
    pointer = Path(root) / (chrom_key(path_id) + '.json')
    return json.loads(pointer.read_text()).get('owner', path_id) if pointer.exists() else path_id


def save_insertion(root, path_id, sequences, bodies, sources, rep_slot, sample_names, owner=None):
    """Call SNPs on the final insertion representative using existing CIGARs."""
    target = Path(root) / chrom_key(path_id)
    target.parent.mkdir(parents=True, exist_ok=True)
    # Each final insertion has one owner. A retry atomically replaces its data.
    directory = Path(tempfile.mkdtemp(prefix=target.name + '-', dir=target.parent))
    writer = Writer(directory)
    template = sequences[rep_slot].upper()
    writer.lengths[path_id] = len(template)
    writer.chroms[chrom_key(path_id)] = path_id
    try:
        for slot, (sequence, body, source) in enumerate(zip(sequences, bodies, sources)):
            sample_index, query, label, label_h, qstart, qend, strand = source
            sample = sample_names[sample_index]
            body = f'{len(template)}=' if slot == rep_slot else body
            if not body:
                raise ValueError(f'missing member alignment for insertion SNPs: {path_id}')
            vcf._fast_validate_alignment_body(body, len(sequence), len(template), context=f'insertion SNP {path_id}')
            _shift, body = vcf._fast_split_row_position_shift(body)
            # The individual caller's LABEL_H distances use genomic left/right,
            # independently of traversal strand.
            pair = vcf.parse_label_h_pair(label_h)
            label_start = min(qstart, qend) - pair[0] if pair is not None else None
            tpos = qpos = 0
            for count_text, op in vcf._FAST_CHUNK_OPS.findall(body):
                count = int(count_text)
                if op in '=MX':
                    writer.coverage[path_id].append([sample, tpos, tpos + count])
                    if op != '=':
                        for i in range(count):
                            ref, alt = template[tpos + i], sequence[qpos + i].upper()
                            if ref == alt or ref not in 'ACGT' or alt not in 'ACGT':
                                continue
                            query_pos = qend - (qpos + i) if strand == '-' else qstart + qpos + i
                            label_point = f'{query_pos - label_start}H' if label_start is not None else '.'
                            writer.add(path_id, tpos + i + 1, ref, alt, sample, [
                                '1', 'INS_SNP', '1', alt, query, f'{query_pos}{strand}', label, label_point,
                            ])
                    tpos += count
                    qpos += count
                elif op in 'DH':
                    tpos += count
                elif op in 'IS':
                    qpos += count
        data = writer.finish()
        data['owner'] = owner or path_id
        # A stable pointer is published last; incomplete/replaced generations
        # are never discovered by SNP jobs.
        vcf.checkpoints.write_json(Path(str(target) + '.json'), data)
    finally:
        writer.close()


def finalize_insertions(root, vcf_paths=None, *, insertion_ids=None):
    """Publish only spools for emitted insertions, optionally collected at copy."""
    root = Path(root)
    root.mkdir(parents=True, exist_ok=True)
    if insertion_ids is not None:
        if vcf_paths is not None:
            raise ValueError('provide VCF paths or insertion IDs, not both')
        pointers = {root / (chrom_key(name) + '.json') for name in insertion_ids}
    elif vcf_paths is None:
        pointers = {path for path in root.glob('*.json') if path.name != 'manifest.json'}
    else:
        pointers = set()
        for path in vcf_paths:
            with vcf.open_text(str(path)) as handle:
                for raw in handle:
                    if raw.startswith('#'):
                        continue
                    # Do not split the potentially huge cohort sample columns.
                    end = -1
                    for _ in range(8):
                        end = raw.find('\t', end + 1)
                        if end < 0:
                            raise ValueError('invalid VCF row while finalizing insertion SNPs')
                    fields = raw[:end].split('\t')
                    if vcf.parse_info_field(fields[7]).get('SVTYPE') == 'INS':
                        pointers.add(root / (chrom_key(fields[2]) + '.json'))
    if any(not path.is_file() for path in pointers):
        raise ValueError('SV output is missing an insertion SNP spool')
    sources = [str(path.resolve()) for path in sorted(pointers)]
    vcf.checkpoints.write_json(root / 'manifest.json', dict(protocol=PROTOCOL, sources=sources))
