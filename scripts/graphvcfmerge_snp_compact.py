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
from contextlib import nullcontext
import hashlib
import heapq
import io
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
import graphvcfmerge_insertion_store as insertion_store

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
    def __init__(self, directory, durable=True, *, memory=False):
        # Insertion records stay in RAM up to 1 MiB, then spill to scratch.
        # They are published with their metadata in a shared bundle at finish.
        self.directory = Path(directory)
        self.durable = durable
        if durable and not memory:
            self.directory.mkdir(parents=True, exist_ok=True)
        self.queries, self.alleles, self.label_groups = {}, {}, {}
        self.metadata = {}
        self.chroms, self.coverage, self.lengths = {}, defaultdict(list), {}
        self.handle = (tempfile.SpooledTemporaryFile(max_size=1024 * 1024, mode='w+b')
                       if memory else (self.directory / "records.bin").open("wb"))
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

    def data(self):
        return dict(protocol=PROTOCOL,
                    directory=str(self.directory.resolve() if self.durable
                                  else os.path.abspath(self.directory)),
                    queries=list(self.queries), alleles=list(self.alleles),
                    label_groups=list(self.label_groups), chroms=self.chroms,
                    coverage=dict(self.coverage), lengths=self.lengths,
                    counts=dict(self.counts), sections=dict(self.sections), metadata=list(self.metadata))

    def finish(self):
        self.close()
        data = self.data()
        vcf.checkpoints.write_json(self.directory / 'source.json', data, self.durable)
        return data


def _scan_one(task):
    index, paths, directory = task
    writer = Writer(Path(directory) / str(index))
    all_samples, metadata = [], []
    try:
        for path in paths:
            samples, alternative_entries = [], set()
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
                            if raw.startswith('##pseudoLinearMapping=<'):
                                entry = vcf.alternative_path_entry(raw.strip())
                                if entry is not None:
                                    alternative_entries.add(entry)
                            metadata.append(raw.strip())
                        elif raw.startswith('#CHROM\t'):
                            samples = raw.rstrip('\r\n').split('\t')[9:]
                            if not samples or len(samples) != len(set(samples)):
                                raise ValueError('expected unique sample names in the VCF header')
                        elif raw.strip() and not raw.startswith('#'):
                            if not samples:
                                raise ValueError('missing #CHROM sample header')
                            row = vcf.parse_vcf_record(raw, number)
                            # A copy's differences against its source (Path=alt
                            # entries only) are annotations, never merged.
                            if (('PAMAP=' in raw or 'PAPATH=alt' in raw)
                                    and vcf.alternative_path_record(row.info, alternative_entries)):
                                continue
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


def _load_source(task, loaded=None):
    # Insertion loci fit in memory. Original chromosome sorting uses bounded spans.
    # ``loaded``: this bundled source's (metadata, records), read in a batch.
    path, key, locus = task
    if isinstance(path, dict):
        data, records = loaded if loaded is not None else insertion_store.read(path)
        if len(records) != sum(data['counts'].values()) * DTYPE.itemsize:
            raise ValueError(f'invalid bundled SNP record count: {path}')
        spans = data['sections'].get(key, [])
        if sum(count for _, count in spans) != data['counts'].get(key, 0):
            raise ValueError(f'invalid bundled SNP section counts: {path}')
        chunks = []
        for offset, count in spans:
            if offset < 0 or count < 0 or offset % DTYPE.itemsize or offset + count * DTYPE.itemsize > len(records):
                raise ValueError(f'invalid bundled SNP section bounds: {path}')
            chunks.append(np.frombuffer(records, dtype=DTYPE, offset=offset, count=count))
        array = np.concatenate(chunks) if chunks else np.empty(0, dtype=DTYPE)
        return array, data, key, locus
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
    if index.get('layout') == MANIFEST_LAYOUT:
        return _concat_finished(manifest, root, index['chromosomes'])
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
    paths = {}
    for source in index['sources']:
        identity = (source['bundle'], source['key'], source['offset'], source['digest']) if isinstance(source, dict) else str(Path(source).resolve())
        paths[identity] = source if isinstance(source, dict) else identity
    finished, taken = _usable_finished(root, manifest['samples'], paths, embedded)
    for key, (name, length) in taken.items():
        if name in originals:
            raise ValueError(f'insertion contig already occurs in chromosome parts: {name!r}')
        if (name in definitions and definitions[name] != length or
                name in lengths and lengths[name] != length):
            raise ValueError(f'conflicting insertion length for {name!r}')
        lengths[name] = length
    # Read source metadata a batch at a time and keep only contig names and
    # lengths: holding every source's metadata at once grows with the cohort.
    pending = [path for identity, path in paths.items() if identity not in embedded
               and not (isinstance(path, dict) and path['key'] in taken)]
    for start in range(0, len(pending), METADATA_BATCH_SOURCES):
        batch = pending[start:start + METADATA_BATCH_SOURCES]
        bundled = [path for path in batch if isinstance(path, dict)]
        metadata_of = {id(path): entry[0] for path, entry in zip(
            bundled, insertion_store.read_many(bundled, with_records=False))}
        for path in batch:
            data = (metadata_of[id(path)] if isinstance(path, dict)
                    else json.loads(Path(path).read_text()))
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
        del metadata_of, bundled, batch
    plan = InsertionPlan(
        [(key, names[key], paths) for key, paths in sorted(by_locus.items(), key=lambda item: names[item[0]])],
        finished)
    # Insertions follow every original chromosome in both the body and header.
    metadata = [line for line in manifest['metadata'] if not (
        line.startswith('##contig=<') and vcf._parse_structured_meta(line, 'contig')['ID'] in lengths)]
    ordered_names = sorted([name for _, name, _ in plan] + [name for name, _ in taken.values()])
    metadata.extend(f'##contig=<ID={name},length={lengths[name]}>' for name in ordered_names)
    return metadata, plan


APPEND_BATCH_LOCI = 512
# Insertion SNPs finished by the SV-merge chromosome jobs (finish_chrom_insertions).
FINISHED_DIRECTORY = 'finished'
FINISHED_PROTOCOL = 2             # LABEL.part rows, LABEL.contigs, LABEL.json
LEGACY_FINISHED_PROTOCOL = 1      # LABEL.json listing every locus with its digest
# manifest.json of a store whose every chromosome finished its insertions:
# only the chromosome labels, nothing per insertion.
MANIFEST_LAYOUT = 'chromosome-finished'


class InsertionPlan(list):
    """Insertion loci still to merge at concat, plus chromosome-finished parts:
    (part, names to take or None for all, sample columns to reorder or None)."""
    def __init__(self, loci=(), finished=(), finished_loci=None):
        super().__init__(loci)
        self.finished = list(finished)
        self.finished_loci = finished_loci


def finished_paths(root, label):
    directory = Path(root) / FINISHED_DIRECTORY
    return tuple(directory / (label + suffix) for suffix in ('.json', '.part', '.contigs'))


def finished_info(root, label):
    """LABEL.json of a complete current-format finished chromosome, else None."""
    info_path, part, contigs = finished_paths(root, label)
    try:
        info = json.loads(info_path.read_text())
    except (OSError, ValueError):
        return None
    if info.get('protocol') != FINISHED_PROTOCOL or not part.is_file() or not contigs.is_file():
        return None
    return info


def discard_finished(root, label):
    """A chromosome about to save its insertions again drops its old part first."""
    for path in finished_paths(root, label):    # .json first: the part is then unusable
        path.unlink(missing_ok=True)


def write_finished_manifest(root, labels):
    """manifest.json naming the finished chromosomes; nothing per insertion."""
    missing = [label for label in labels if finished_info(root, label) is None]
    if missing:
        raise ValueError(f'{len(missing)} chromosome(s) have not finished their insertion SNPs: '
                         f'{missing[:6]}')
    vcf.checkpoints.write_json(Path(root) / 'manifest.json',
                               dict(protocol=PROTOCOL, layout=MANIFEST_LAYOUT, chromosomes=list(labels)))


class _InsertionHeader:
    """Header lines: the small-variant metadata, then one ##contig per insertion
    in name order, streamed from the finished chromosomes' .contigs files."""

    def __init__(self, metadata, contigs, replaced):
        self.metadata, self.contigs, self.replaced = metadata, contigs, replaced

    def __iter__(self):
        for line in self.metadata:
            if not (line.startswith('##contig=<') and
                    vcf._parse_structured_meta(line, 'contig')['ID'] in self.replaced):
                yield line
        for name, length in _contig_rows(self.contigs):
            yield f'##contig=<ID={name},length={length}>'


def _contig_rows(paths):
    def rows(path):
        with open(path) as handle:
            for raw in handle:
                name, length = raw.rstrip('\n').split('\t')
                yield name, length
    return heapq.merge(*(rows(path) for path in paths))


def _concat_finished(manifest, root, labels):
    """Plan a store whose chromosomes all finished: their parts are streamed,
    nothing is read per insertion."""
    originals = set(manifest.get('loci', manifest['chroms']).values())
    definitions = {}
    for line in manifest['metadata']:
        if line.startswith('##contig=<'):
            data = vcf._parse_structured_meta(line, 'contig')
            definitions[data['ID']] = data.get('length')
    finished, contigs = [], []
    for label in labels:
        info = finished_info(root, label)
        if info is None:
            raise ValueError(f'{root}: chromosome {label} has not finished its insertion SNPs; '
                             'run graphvcfmerge.py --sv --stage finish-insertions')
        _info, part, contig_path = finished_paths(root, label)
        order = None
        if info['samples'] != list(manifest['samples']):
            if sorted(info['samples']) != sorted(manifest['samples']):
                raise ValueError(f'{part}: samples differ from the small-variant merge')
            position = {sample: index for index, sample in enumerate(info['samples'])}
            order = [position[sample] for sample in manifest['samples']]
        finished.append((str(part), None, order))
        contigs.append(contig_path)
    # One pass checks names against the chromosomes and earlier definitions.
    replaced, previous, loci = set(), None, 0
    for name, length in _contig_rows(contigs):
        if name == previous:
            # The store holds one entry per name: the other chromosome's was overwritten.
            raise ValueError(f'insertion name {name!r} is used by two chromosomes')
        if name in originals:
            raise ValueError(f'insertion contig already occurs in chromosome parts: {name!r}')
        if name in definitions:
            if definitions[name] != length:
                raise ValueError(f'conflicting insertion length for {name!r}')
            replaced.add(name)
        previous, loci = name, loci + 1
    return (_InsertionHeader(manifest['metadata'], contigs, replaced),
            InsertionPlan([], finished, finished_loci=loci))
# Insertion sources whose metadata concat_insertions holds at once.
METADATA_BATCH_SOURCES = 4096


def append_insertions(handle, manifest, plan, root=None):
    """Append insertion SNP sites after the original chromosomes, by insertion name.

    Loci finished by the SV-merge chromosome jobs are copied from their parts;
    the rest are merged here. Both streams are already in name order, so they
    are interleaved without loading either."""
    finished = getattr(plan, 'finished', [])
    if not finished:
        total = 0
        for count, text in _merge_planned_loci(manifest, plan, root):
            handle.write(text)
            total += count
        if plan:
            print(f'[merge:snp] concat appended {total} sites from {len(plan)} insertion loci', file=sys.stderr)
        return total
    streams = [_finished_rows(*item) for item in finished]
    streams.append(line for _count, text in _merge_planned_loci(manifest, plan, root)
                   for line in text.splitlines(keepends=True))
    total = 0
    for line in heapq.merge(*streams, key=_row_chrom):
        handle.write(line)
        total += 1
    finished_loci = plan.finished_loci
    if finished_loci is None:
        finished_loci = sum(len(names) for _part, names, _order in finished)
    loci = len(plan) + finished_loci
    print(f'[merge:snp] concat appended {total} sites from {loci} insertion loci '
          f'({loci - len(plan)} finished by chromosome jobs)', file=sys.stderr)
    return total


def _merge_planned_loci(manifest, plan, root):
    """Yield (sites, text) per planned locus, reading APPEND_BATCH_LOCI loci at a time.

    Bundled sources and realignment records are read with one open per bundle
    per batch, not one per locus."""
    realign_root = Path(root) / 'realign' if root else None
    # Older runs wrote realignment records as one JSON file per path.
    legacy = ({path.stem for path in realign_root.glob('*.json')}
              if realign_root is not None and realign_root.is_dir() else set())
    for start in range(0, len(plan), APPEND_BATCH_LOCI):
        chunk = plan[start:start + APPEND_BATCH_LOCI]
        bundled = [path for _key, _name, paths in chunk for path in paths if isinstance(path, dict)]
        loaded = {id(path): entry for path, entry in zip(bundled, insertion_store.read_many(bundled))}
        realignments = {}
        if realign_root is not None:
            found = insertion_store.find_many(realign_root, [key for key, _name, _paths in chunk])
            keys = list(found)
            realignments = dict(zip(keys, (entry[0] for entry in insertion_store.read_many(
                [found[key] for key in keys], with_records=False))))
        for key, name, paths in chunk:
            text = io.StringIO()
            count = _emit_insertion_locus(text, manifest['samples'], key, name, paths, loaded,
                                          realignments.get(key), realign_root, legacy)
            yield count, text.getvalue()
        del loaded, realignments


def _row_chrom(line):
    return line[:line.index('\t')]


def _finished_rows(part, names, order=None):
    """A finished part's rows: those of NAMES (None: all), sample columns in ORDER."""
    with open(part) as handle:
        for line in handle:
            if names is not None and _row_chrom(line) not in names:
                continue
            if order is not None:
                fields = line.rstrip('\n').split('\t')
                line = '\t'.join(fields[:9] + [fields[9 + index] for index in order]) + '\n'
            yield line


def _usable_finished(root, samples, paths, embedded):
    """Chromosome-finished parts whose loci still match this store.

    A locus is taken from a part only when its single bundled source has the
    digest the chromosome job merged, and the part used this sample order;
    everything else is merged at concat as before."""
    directory = Path(root) / FINISHED_DIRECTORY
    if not directory.is_dir():
        return [], {}
    sources = defaultdict(list)
    for identity, path in paths.items():
        if identity in embedded:
            continue
        if not isinstance(path, dict):
            return [], {}  # Legacy per-locus spools: keys unknown without reading them.
        sources[path['key']].append(path)
    finished, taken = [], {}
    for index in sorted(directory.glob('*.json')):
        try:
            info = json.loads(index.read_text())
        except (OSError, ValueError):
            continue
        part = index.with_suffix('.part')
        if info.get('protocol') != LEGACY_FINISHED_PROTOCOL or not part.is_file():
            continue
        if info.get('samples') != list(samples):
            print(f'[merge:snp] {index.name}: sample order differs; merging its loci at concat',
                  file=sys.stderr)
            continue
        names = set()
        for key, name, length, digest in info.get('loci', []):
            found = sources.get(key, [])
            if key not in taken and len(found) == 1 and found[0]['digest'] == digest:
                taken[key] = (name, str(length))
                names.add(name)
        if names:
            finished.append((str(part), names, None))
    return finished, taken


DIRECT_READ_BATCH = 1    # one insertion in memory per worker: a large one can take GBs


def backfill_index(root):
    """Entry offsets for a backfill, whose SV merge (the only writer) is done:
    manifest.json's sources, plus one walk of realign/. None when the store
    has no manifest or holds legacy per-insertion spools."""
    root = Path(root)
    path = root / 'manifest.json'
    if not path.is_file():
        return None
    sources = json.loads(path.read_text()).get('sources')
    if sources is None or any(not isinstance(source, dict) for source in sources):
        return None
    realign_root = root / 'realign'
    realign, legacy = {}, set()
    if realign_root.is_dir():
        legacy = {path.stem for path in realign_root.glob('*.json')}
        realign = {source['key']: source for source in insertion_store.sources(realign_root)}
    return {source['key']: source for source in sources}, realign, legacy


def finish_chrom_insertions(root, part_paths, samples, label, index=None):
    """Merge one chromosome's insertion SNPs right after its SV merge.

    The SV merge has just saved every insertion of this chromosome (its INS
    rows' owner paths), with nested realignments. Each is merged alone, as at
    concat, into ROOT/finished/LABEL.part in insertion-name order, with one
    "name<TAB>length" line per insertion in LABEL.contigs (the ##contig header
    lines). LABEL.json, written last, records only the sample order and counts.

    INDEX (backfill_index) is for a backfill only: entries are then read at
    their indexed offsets through open, unlocked bundles, a batch at a time in
    file order, instead of locating and locking per insertion."""
    from graphvcfmerge_snp import atomic_output
    root = Path(root)
    owned = {}
    for part in part_paths:
        if not part or not Path(part).is_file():
            continue
        with open(part) as handle:
            for raw in handle:
                fields = raw.split('\t', 8)
                if len(fields) >= 8 and 'SVTYPE=INS' in fields[7].split(';'):
                    name = vcf.checkpoints.insertion_snp_owner(fields[2], fields[7])
                    owned[chrom_key(name)] = name
    realign_root = root / 'realign'
    if index is None:
        found = insertion_store.find_many(root, owned)
        legacy, realign_found = set(), {}
        if realign_root.is_dir():
            legacy = {path.stem for path in realign_root.glob('*.json')}
            realign_found = insertion_store.find_many(realign_root, found)
    else:
        sources, realign_sources, legacy = index
        found = {key: sources[key] for key in owned if key in sources}
        realign_found = {key: realign_sources[key] for key in found if key in realign_sources}
    missing = sorted(owned[key] for key in owned.keys() - found.keys())
    if missing:
        raise ValueError(f'{label}: {len(missing)} insertion(s) have no saved SNP entry in {root}: '
                         f'{missing[:6]}')
    info_path, part_path, contigs_path = finished_paths(root, label)
    info_path.unlink(missing_ok=True)
    samples = list(samples)
    allowed = set(samples)
    loci = total = 0
    items = sorted((owned[key], key) for key in found)
    batch_size = 1 if index is None else DIRECT_READ_BATCH
    with atomic_output(part_path) as handle, atomic_output(contigs_path) as contigs, \
            (insertion_store.DirectReader() if index is not None else nullcontext()) as reader:
        read_many = insertion_store.read_many if reader is None else reader.read_many
        for start in range(0, len(items), batch_size):
            batch = items[start:start + batch_size]
            entries = read_many([found[key] for _name, key in batch])
            realigned = [key for _name, key in batch if key in realign_found]
            realigned = dict(zip(realigned, read_many([realign_found[key] for key in realigned],
                                                      with_records=False)))
            for (name, key), entry in zip(batch, entries):
                source = found[key]
                if any(query[0] not in allowed for query in entry[0]['queries']):
                    raise ValueError(f'insertion SNP source contains an unknown sample: {source}')
                realignment = realigned[key][0] if key in realigned else None
                total += _emit_insertion_locus(handle, samples, key, name, [source], {id(source): entry},
                                               realignment, realign_root, legacy)
                contigs.write(f"{name}\t{entry[0]['lengths'][name]}\n")
                loci += 1
            del entries, realigned
    vcf.checkpoints.write_json(info_path, dict(protocol=FINISHED_PROTOCOL, samples=samples,
                                               insertions=loci, sites=total))
    print(f'[merge:snp] {label}: finished {total} insertion SNP sites from {loci} insertions',
          file=sys.stderr, flush=True)
    return total


def _emit_insertion_locus(handle, samples, key, name, paths, loaded, realignment, realign_root, legacy):
    """Merge one insertion's SNP observations from all its sources and write its sites."""
    queries, alleles, groups = {}, {}, {}
    coverage, arrays = defaultdict(list), []
    for path in paths:
        array, data, _, _ = _load_source((path, key, 0), loaded.get(id(path)))
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
    query_list, allele_list = list(queries), list(alleles)
    if realignment is not None:
        # --exact nested realignment records for this insertion path.
        ordered = apply_realignment(realign_root / (key + '.json'), ordered, query_list,
                                    allele_list, kind='INS_SNP', data=realignment)
    elif key in legacy:
        ordered = apply_realignment(realign_root / (key + '.json'), ordered, query_list,
                                    allele_list, kind='INS_SNP')
    return write_locus(handle, name, ordered, query_list, allele_list, list(groups),
                       coverage, samples)


def merge_chrom(root, chrom, processes=1):
    manifest = json.loads((Path(root) / 'manifest.json').read_text())
    key = chrom if chrom in manifest['chroms'] else chrom_key(chrom)
    if key not in manifest['chroms']:
        raise ValueError(f'unknown SNP chromosome: {chrom}')
    with vcf.checkpoints.chrom_lock(manifest['directory'], key):
        return _emit_chrom(manifest, key, processes, root=root)


def apply_realignment(path, ordered, queries, alleles, kind='SNP', data=None):
    """Apply graphvcfmerge --exact realignment records to one sorted locus.

    Drops remove the named SNP observations (superseded by a moved SV's
    flank representation); adds are that representation's new SNPs. A drop
    that matches nothing is an error: the SV merge replaced it. kind is SNP
    for reference loci and INS_SNP for insertion paths. ``data``: the
    records, already read (batched callers).
    """
    path = Path(path)
    source = (insertion_store.find(path.parent, path.stem)
              if data is None and re.fullmatch(r'[0-9a-f]{64}', path.stem) else None)
    if data is not None:
        pass
    elif source is not None:
        data, _ = insertion_store.read(source)
    elif path.is_file():
        data = json.loads(path.read_text())
    else:
        return ordered
    if not data['drops'] and not data['adds']:
        return ordered
    by_query = defaultdict(set)
    for number, query in enumerate(queries):
        if query[3] == kind and query[4] == '1':
            by_query[(query[0], query[1], query[2])].add(number)
    keep = np.ones(len(ordered), dtype=bool)
    for sample, pos, alt, contig, qpos, strand in data['drops']:
        ids = by_query.get((sample, contig, strand), set())
        first = np.searchsorted(ordered['pos'], pos, side='left')
        last = np.searchsorted(ordered['pos'], pos, side='right')
        rows = [row for row in range(first, last)
                if keep[row] and int(ordered['query'][row]) in ids
                and int(ordered['query_pos'][row]) == qpos
                and int(ordered['bases'][row]) & 15 == BASE_INDEX[alt.upper()]]
        if not rows:
            raise ValueError(f'realignment drop not found: {sample} {data["chrom"]}:{pos} {alt} {contig}:{qpos}{strand}')
        keep[rows] = False
    added = np.empty(len(data['adds']), dtype=DTYPE)
    query_index = {tuple(query): number for number, query in enumerate(queries)}
    allele = alleles.index('.') if '.' in alleles else None
    if allele is None:
        alleles.append('.')
        allele = len(alleles) - 1
    for row, (sample, pos, ref, alt, contig, qpos, strand) in enumerate(data['adds']):
        query = (sample, contig, strand, kind, '1', False, 'PASS')
        if query not in query_index:
            queries.append(query)
            query_index[query] = len(queries) - 1
        added[row] = (u32(pos), pack_bases(ref, alt), u32(qpos), query_index[query], allele, 0)
    result = np.concatenate([ordered[keep], added])
    result.sort(kind='mergesort', order=DTYPE.names)
    print(f'[merge:snp] {data["chrom"]}: realignment dropped {int((~keep).sum())} and added '
          f'{len(added)} SNP observation(s)', file=sys.stderr, flush=True)
    return result


def _emit_chrom(manifest, key, processes, root=None):
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
            if root is not None:
                ordered = apply_realignment(Path(root) / 'realign' / (locus + '.json'),
                                            ordered, queries, alleles)
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


_OWNERS = {}    # (root, insertion key) -> owner, from this process's puts or a prefetch


def prefetch_owners(root, path_ids):
    """Owners of many insertion paths with one lookup per bundle (nested rounds)."""
    keys = {chrom_key(path_id) for path_id in path_ids} - {
        key for (cached_root, key) in _OWNERS if cached_root == str(root)}
    found = insertion_store.find_many(root, keys)
    for key, (data, _) in zip(found, insertion_store.read_many(list(found.values()),
                                                               with_records=False)):
        _OWNERS[str(root), key] = data.get('owner')
    for key in keys - set(found):
        _OWNERS[str(root), key] = None      # no spool (e.g. a deletion parent)


def insertion_owner(root, path_id):
    cached = _OWNERS.get((str(root), chrom_key(path_id)), False)
    if cached is not False:
        return cached or path_id
    source = insertion_store.find(root, chrom_key(path_id))
    if source is not None:
        return insertion_store.read(source, with_records=False)[0].get('owner', path_id)
    pointer = Path(root) / (chrom_key(path_id) + '.json')
    return json.loads(pointer.read_text()).get('owner', path_id) if pointer.exists() else path_id


def save_insertion(root, path_id, sequences, bodies, sources, rep_slot, sample_names, owner=None):
    """Call SNPs on the final insertion representative using existing CIGARs."""
    # Dictionary IDs and record offsets remain local to this insertion. Bundle
    # append publishes its complete metadata and records together, even when empty.
    writer = Writer(root, durable=False, memory=True)
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
        data = writer.data()
        data['owner'] = owner or path_id
        insertion_store.put(root, chrom_key(path_id), data, writer.handle)
        _OWNERS[str(root), chrom_key(path_id)] = data['owner']
    finally:
        writer.close()


def finalize_insertions(root, vcf_paths=None, *, insertion_ids=None):
    """Publish only spools for emitted insertions, optionally collected at copy."""
    root = Path(root)
    root.mkdir(parents=True, exist_ok=True)
    if insertion_ids is not None:
        if vcf_paths is not None:
            raise ValueError('provide VCF paths or insertion IDs, not both')
        keys = {chrom_key(name) for name in insertion_ids}
    elif vcf_paths is None:
        keys = {path.stem for path in root.glob('*.json') if path.name != 'manifest.json'}
        keys.update(source['key'] for source in insertion_store.sources(root))
    else:
        keys = set()
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
                        keys.add(chrom_key(vcf.checkpoints.insertion_snp_owner(fields[2], fields[7])))
    bundled = insertion_store.find_many(root, keys)
    sources = [bundled[key] if key in bundled else legacy_insertion_source(root, key)
               for key in sorted(keys)]
    vcf.checkpoints.write_json(root / 'manifest.json', dict(protocol=PROTOCOL, sources=sources))


def insertion_source(root, key):
    """Prefer the bundled generation; existing per-insertion spools stay readable."""
    source = insertion_store.find(root, key)
    if source is not None:
        return source
    return legacy_insertion_source(root, key)


def legacy_insertion_source(root, key):
    """A per-insertion spool of an older run (key known not to be bundled)."""
    pointer = Path(root) / (key + '.json')
    if not pointer.is_file():
        raise ValueError(f'SV output is missing an insertion SNP spool: {pointer}')
    return str(pointer.resolve())


def save_insertion_realignment(root, path_id, records):
    insertion_store.put(Path(root) / 'realign', chrom_key(path_id), records)
