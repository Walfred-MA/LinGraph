#!/usr/bin/env python3
"""Recover insertion SNPs from merged grVCFs and indexed sample assemblies.

Reuse saved member alignments and compose partial insertion-to-insertion
coordinates before merging SNP calls. No SV clustering or alignment is
performed. Bare X/M operations require assembly bases.
"""
from __future__ import annotations

import argparse
from collections import OrderedDict, defaultdict
from contextlib import ExitStack
import gzip
import itertools
import os
from pathlib import Path
import re
import shutil
import sqlite3
import sys
import tempfile
import time
from types import SimpleNamespace

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
sys.path.insert(0, str(ROOT / "scripts"))
import graphvcfmerge as vcf
import graphvcfmerge_snp_compact as compact
from graph_cigar_payloads import query_only_graph_cigar
from tools.grvcf_to_vcf import Catalog, Sequences, OP, open_input, parse_info
from tools.vcf_utils import IndexedFasta
from tools.insertion_coordinates import Coordinates, cigar_chunks

DNA = re.compile(r'[ACGTRYSWKMBDHVN]+', re.I)
COMPLEMENT = str.maketrans('ACGTRYSWKMBDHVN', 'TGCAYRSWMKVHDBN')


def assemblies(path):
    """Two tab-separated columns: exact VCF sample name, indexed FASTA path."""
    result = {}
    if path is None:
        return result
    path = Path(path).resolve()
    for number, raw in enumerate(path.read_text().splitlines(), 1):
        if not raw.strip() or raw.startswith('#'):
            continue
        fields = raw.split('\t')
        if len(fields) != 2 or not all(fields):
            raise ValueError(f'{path}:{number}: expected SAMPLE<TAB>FASTA')
        sample, fasta = fields
        if sample in result:
            raise ValueError(f'{path}:{number}: duplicate sample {sample!r}')
        fasta = Path(fasta).expanduser()
        result[sample] = str((fasta if fasta.is_absolute() else path.parent / fasta).resolve())
    return result


class AssemblyReaders:
    def __init__(self, mapping):
        self.mapping, self.opened = mapping, OrderedDict()

    def fetch(self, sample, contig, start, end, strand):
        if sample not in self.mapping:
            raise ValueError(f'{sample}: alternate bases are absent from the VCF; '
                             'supply --assemblies with this sample and its FASTA')
        path = self.mapping[sample]
        if path not in self.opened:
            if path.endswith(('.gz', '.bgz')):
                raise ValueError(f'{sample}: use an uncompressed assembly FASTA with its .fai index')
            if len(self.opened) >= 16:
                self.opened.popitem(last=False)[1].close()
            self.opened[path] = IndexedFasta(path)
        self.opened.move_to_end(path)
        sequence = self.opened[path].fetch(contig, start, end)
        if len(sequence) != end - start or (sequence and not DNA.fullmatch(sequence)):
            raise ValueError(f'{sample}:{contig}:{start}-{end}: truncated or invalid assembly sequence')
        return sequence.translate(COMPLEMENT)[::-1] if strand == '-' else sequence

    def close(self):
        for reader in self.opened.values():
            reader.close()
        self.opened.clear()


def observations(field, fmt):
    """Decode parallel FORMAT vectors without splitting escaped contig names."""
    if fmt == 'HSV':
        if field in ('.', '0', '0/0', './.'):
            return []
        values = [vcf.parse_hsv_allele(x) for x in vcf.split_sample_alleles(field, fmt)]
        if not values or any(x is None for x in values):
            raise ValueError('malformed legacy HSV observation')
        return [dict(zip(vcf.SV_EVENT_FIELDS, x)) for x in values]
    names = fmt.split(':')
    required = {'GT', 'TYPE', 'EXTENDGRAPHCIGAR', 'ASSEMBLYCONTIG', 'QUERYCOORD'}
    if len(names) != len(set(names)) or not required <= set(names):
        raise ValueError('expected merged SV FORMAT with CIGAR and query coordinates')
    if field in ('.', '0', '0/0', './.'):
        return []
    columns = field.split(':')
    if len(columns) != len(names):
        raise ValueError('sample field count differs from FORMAT')
    values = dict(zip(names, columns))
    if values['GT'] in ('.', './.', '0', '0/0'):
        return []
    if values['GT'] != '1':
        raise ValueError(f"expected merged haploid GT, got {values['GT']!r}")
    vectors = {name: text.split(',') for name, text in values.items() if name != 'GT'}
    sizes = {len(vector) for vector in vectors.values()}
    if len(sizes) != 1:
        raise ValueError('inconsistent observation counts across FORMAT columns')
    return [dict(GT='1', **{name: vcf.vcf_unescape(vector[i]) if vector[i] != '.' else '.'
                           for name, vector in vectors.items()}) for i in range(sizes.pop())]


def operations(cigar, identifier):
    text = query_only_graph_cigar(cigar)
    shift, text = vcf._fast_split_row_position_shift(text)
    if shift or not text.startswith('>') or any(x in text[1:] for x in '<>'):
        raise ValueError('expected one forward alignment against the merged insertion representative')
    body = text[1:]
    if ':' in body:
        name, body = body.split(':', 1)
        if name != identifier:
            raise ValueError(f'alignment target {name!r} differs from insertion {identifier!r}')
    result, end = [], 0
    for match in OP.finditer(body):
        n, op, payload = int(match[1]), match[2], match[3].upper()
        if match.start() != end or n < 1:
            raise ValueError('malformed member alignment')
        if payload and (op not in 'XIS' or len(payload) != n or not DNA.fullmatch(payload)):
            raise ValueError(f'invalid {op} payload')
        result.append((n, op, payload))
        end = match.end()
    if not result or end != len(body):
        raise ValueError('malformed member alignment suffix')
    return result


def member_operations(cigar, identifier, coordinates):
    # Old row-position padding describes the genomic breakpoint, not the
    # insertion's local reference origin. Named graph chunks have their own H.
    _shift, text = vcf._fast_split_row_position_shift(cigar)
    chunks = cigar_chunks(text)
    result, qpos, main_pos = [], 0, 0
    for direction, name, ops in chunks:
        target = name or identifier
        tpos = 0 if name else main_pos
        for n, op, payload in ops:
            if payload and (op not in 'XIS' or len(payload) != n or not DNA.fullmatch(payload)):
                raise ValueError(f'invalid {op} payload')
            result.append((n, op, payload, target, direction == '<', tpos, qpos))
            if op in '=MXDH':
                tpos += n
            if op in '=MXIS':
                qpos += n
        length = coordinates.length(target)
        if tpos > length or (len(chunks) == 1 and name is None and tpos != length):
            raise ValueError(f'alignment spans template={tpos}; expected template length={length}')
        if name is None:
            main_pos = tpos
    return result, qpos


def extract_row(fields, samples, readers, writer, coordinates):
    identifier = fields[2]
    if len(fields) != 9 + len(samples):
        raise ValueError('sample column count differs from header')
    count = 0
    for sample, field in zip(samples, fields[9:]):
        for value in observations(field, fields[8]):
            if value['TYPE'] != 'INS' or value['GT'] != '1':
                raise ValueError('expected insertion observation in merged INS row')
            qstart, qend, strand = vcf.parse_query_range_required(value['QUERYCOORD'])
            ops, qspan = member_operations(value['EXTENDGRAPHCIGAR'], identifier, coordinates)
            if qspan != qend - qstart:
                raise ValueError(f'{sample}: alignment spans query={qspan}; expected QUERYCOORD length={qend-qstart}')
            needs_bases = any(op == 'M' or (op == 'X' and not payload) for _, op, payload, *_ in ops)
            sequence = (readers.fetch(sample, value['ASSEMBLYCONTIG'], qstart, qend, strand)
                        if needs_bases or sample in readers.mapping else None)
            pair = vcf.parse_label_h_pair(value.get('LABEL_H', '.'))
            label_start = qstart - pair[0] if pair is not None else None
            for n, op, payload, target, reverse, tpos, qpos in ops:
                if sequence is not None and payload and sequence[qpos:qpos+n] != payload:
                    raise ValueError(f'{sample}: assembly bases disagree with embedded CIGAR payload')
                if op in '=MX':
                    length = coordinates.length(target)
                    lo, hi = (length-tpos-n, length-tpos) if reverse else (tpos, tpos+n)
                    template_segment = coordinates.reference(target, lo, hi)
                    if reverse:
                        template_segment = template_segment.translate(COMPLEMENT)[::-1]
                    if op == '=' and sequence is not None:
                        if template_segment != sequence[qpos:qpos+n] and any(a != b and a in 'ACGT' and b in 'ACGT'
                               for a, b in zip(template_segment, sequence[qpos:qpos+n])):
                            raise ValueError(f'{sample}: assembly disagrees with an = alignment; '
                                             'check the FASTA, contig and query coordinates')
                    # Also examine = blocks: the source representative itself
                    # can differ from the insertion to which it was mapped.
                    for offset, width, destination, start, step in coordinates.oriented(target, tpos, n, reverse):
                        low = min(start, start+(width-1)*step)
                        ref_segment = coordinates.reference(destination, low, low+width)
                        query_offset = qpos + offset
                        if sequence is not None:
                            segment = sequence[query_offset:query_offset+width]
                        elif op == '=':
                            segment = template_segment[offset:offset+width]
                        else:
                            segment = payload[offset:offset+width]
                        if step < 0:
                            segment = segment.translate(COMPLEMENT)
                            ref_segment = ref_segment[::-1]
                        writer.chroms[compact.chrom_key(destination)] = destination
                        writer.lengths[destination] = coordinates.length(destination)
                        writer.coverage[destination].append([sample, low, low+width])
                        if segment == ref_segment:
                            continue
                        for i, (ref, alt) in enumerate(zip(ref_segment, segment)):
                            if ref == alt or ref not in 'ACGT' or alt not in 'ACGT':
                                continue
                            point = qend - query_offset - i if strand == '-' else qstart + query_offset + i
                            label = f'{point-label_start}H' if label_start is not None else '.'
                            writer.add(destination, start+i*step+1, ref, alt, sample, [
                                '1', 'INS_SNP', '1', alt, value['ASSEMBLYCONTIG'],
                                f'{point}{strand}', value.get('ALLELENAME', '.'), label,
                            ])
                            count += 1
    return count


def index_inputs(paths, db, work):
    db.executescript('CREATE TABLE alleles (id TEXT PRIMARY KEY, info TEXT); '
                     'CREATE TABLE resolved (id TEXT PRIMARY KEY, sequence TEXT); '
                     'CREATE TABLE records (id TEXT PRIMARY KEY, path TEXT, offset INTEGER, '
                     'bytes INTEGER, line INTEGER);')
    samples = None
    indexed, progress = 0, time.monotonic()
    for index, original in enumerate(paths):
        path = original
        if str(path).endswith(('.gz', '.bgz')):
            path = work / f'input-{index}.vcf'
            with gzip.open(original, 'rb') as source, path.open('wb') as out:
                shutil.copyfileobj(source, out)
        header = None
        with open(path, 'rb') as source:
            number = 0
            while True:
                offset = source.tell()
                raw = source.readline()
                if not raw:
                    break
                number += 1
                if raw.startswith(b'#CHROM\t'):
                    if header is not None:
                        raise ValueError(f'{original}: duplicate #CHROM header')
                    header = raw.rstrip(b'\r\n').decode().split('\t')[9:]
                    if not header or len(set(header)) != len(header):
                        raise ValueError(f'{original}: expected unique cohort sample names')
                    if samples is not None and header != samples:
                        raise ValueError(f'{original}: input VCF sample names/order differ')
                    samples = header
                elif raw.startswith(b'#') or not raw.strip():
                    continue
                else:
                    if header is None:
                        raise ValueError(f'{original}:{number}: missing #CHROM header')
                    fields = raw.split(b'\t', 8)
                    if len(fields) != 9:
                        raise ValueError(f'{original}:{number}: invalid VCF row')
                    info_text, identifier = fields[7].decode(), fields[2].decode()
                    kind = parse_info(info_text).get('SVTYPE')
                    if kind not in ('INS', 'SUB'):
                        continue
                    if identifier == '.' or not re.fullmatch(r'[^\s,:<>]+', identifier):
                        raise ValueError(f'{original}:{number}: insertion requires a usable unique ID')
                    try:
                        db.execute('INSERT INTO alleles VALUES (?,?)', (identifier, info_text))
                        if kind == 'INS':
                            db.execute('INSERT INTO records VALUES (?,?,?,?,?)',
                                       (identifier, str(path), offset, len(raw), number))
                            indexed += 1
                    except sqlite3.IntegrityError as error:
                        raise ValueError(f'{original}:{number}: duplicate insertion ID {identifier!r}') from error
                    if time.monotonic() - progress >= 15:
                        print(f'[insertion-snps] indexed {indexed} insertions; {original}:{number}',
                              file=sys.stderr, flush=True)
                        progress = time.monotonic()
        if header is None:
            raise ValueError(f'{original}: missing #CHROM header')
    db.commit()
    return samples


_WORKER = None


def init_worker(database, mapping, samples, work, aux_fasta):
    global _WORKER
    db = sqlite3.connect(Path(database).as_uri() + '?mode=ro', uri=True)
    stack = ExitStack()
    catalog = Catalog(aux_fasta[0], aux_fasta[1:], stack) if aux_fasta else SimpleNamespace(readers={})
    _WORKER = (db, AssemblyReaders(mapping), samples, work, OrderedDict(), Coordinates(db, catalog), stack)


def close_worker():
    db, readers, _samples, _work, inputs, _coordinates, stack = _WORKER
    db.close()
    readers.close()
    for source in inputs.values():
        source.close()
    stack.close()


def process_batch(task):
    number, rows = task
    db, readers, samples, work, inputs, coordinates, _stack = _WORKER
    directory, count = Path(work) / f'batch-{number}', 0
    writer = compact.Writer(directory)
    try:
        for identifier, path, offset, size, line in rows:
            try:
                if path not in inputs:
                    if len(inputs) >= 8:
                        inputs.popitem(last=False)[1].close()
                    inputs[path] = open(path, 'rb')
                inputs.move_to_end(path)
                source = inputs[path]
                source.seek(offset)
                raw = source.read(size)
                fields = raw.decode().rstrip('\r\n').split('\t')
                if fields[2] != identifier:
                    raise ValueError('input changed after indexing')
                count += extract_row(fields, samples, readers, writer, coordinates)
            except (ValueError, KeyError, OSError) as error:
                raise ValueError(f'{path}:{line}: insertion {identifier}: {error}') from error
        for name, intervals in writer.coverage.items():
            by_sample = defaultdict(list)
            for sample, start, end in intervals:
                by_sample[sample].append((start, end))
            writer.coverage[name] = [[sample, start, end] for sample, values in by_sample.items()
                                     for start, end in vcf._merge_intervals(values)]
        data = writer.finish()
    finally:
        writer.close()
    return str(directory / 'source.json'), data['chroms'], count, len(rows)


def batches(db):
    batch, size, number = [], 0, 0
    for row in db.execute('SELECT id,path,offset,bytes,line FROM records ORDER BY id'):
        batch.append(row)
        size += row[3]
        if size >= 8 * 1024 * 1024 or len(batch) >= 64:
            yield number, batch
            batch, size, number = [], 0, number + 1
    if batch:
        yield number, batch


def read_header(path):
    metadata = []
    with open_input(path) as source:
        for raw in source:
            if raw.startswith('##'):
                metadata.append(raw.rstrip('\r\n'))
            elif raw.startswith('#CHROM\t'):
                return metadata, raw.rstrip('\r\n').split('\t')[9:]
            else:
                raise ValueError(f'{path}: malformed VCF header')
    raise ValueError(f'{path}: missing #CHROM header')


def run(args):
    output = Path(args.output).resolve()
    paths = [Path(path).resolve() for path in args.input]
    mapping = assemblies(args.assemblies)
    protected = paths + [Path(p) for p in mapping.values()] + [Path(p+'.fai') for p in mapping.values()]
    protected += [Path(p).resolve() for p in [args.assemblies, args.append_to, *args.aux_fasta] if p]
    if output in protected:
        raise ValueError('output must not replace an input VCF, assembly, or sample map')
    if output.exists() and not args.force:
        raise ValueError(f'output already exists: {output}; use --force to replace it atomically')
    if args.processes < 1:
        raise ValueError('--processes must be positive')
    output.parent.mkdir(parents=True, exist_ok=True)
    stamps = {path: vcf.checkpoints.file_stamp(path) for path in paths}
    with tempfile.TemporaryDirectory(prefix='insertion-snps-', dir=args.tmpdir or output.parent) as temporary, ExitStack() as stack:
        work = Path(temporary)
        database = work / 'index.sqlite'
        db = sqlite3.connect(database)
        stack.callback(db.close)
        print('[insertion-snps] indexing insertion IDs and saved alignments', file=sys.stderr, flush=True)
        samples = index_inputs(paths, db, work)
        catalog = Catalog(args.aux_fasta[0], args.aux_fasta[1:], stack) if args.aux_fasta else SimpleNamespace(readers={})
        resolver = Sequences(db, catalog)
        count = db.execute('SELECT COUNT(*) FROM records').fetchone()[0]
        progress = time.monotonic()
        for number, (identifier,) in enumerate(db.execute('SELECT id FROM records ORDER BY id'), 1):
            resolver.target(identifier)
            if time.monotonic() - progress >= 15:
                print(f'[insertion-snps] resolved {number}/{count} insertion templates', file=sys.stderr, flush=True)
                progress = time.monotonic()
        db.commit()
        coordinates = Coordinates(db, catalog, build=True)
        progress = time.monotonic()
        for number, (identifier,) in enumerate(db.execute('SELECT id FROM records ORDER BY id'), 1):
            coordinates.segments(identifier)
            if time.monotonic() - progress >= 15:
                print(f'[insertion-snps] mapped {number}/{count} insertion coordinate systems', file=sys.stderr, flush=True)
                progress = time.monotonic()
        db.commit()
        workers = min(args.processes, max(1, count))
        print(f'[insertion-snps] extracting {count} insertion(s) with {workers} worker(s)', file=sys.stderr, flush=True)
        done = observations_count = 0
        sources = defaultdict(list)
        pool = stack.enter_context(vcf.mp_context().Pool(workers, initializer=init_worker,
                                   initargs=(str(database), mapping, samples, str(work), args.aux_fasta))) if workers > 1 else None
        if pool is None:
            init_worker(str(database), mapping, samples, str(work), args.aux_fasta)
            stack.callback(close_worker)
        tasks = iter(batches(db))
        while window := list(itertools.islice(tasks, workers * 2)):
            results = pool.imap(process_batch, window) if pool else map(process_batch, window)
            for path, chroms, written, processed in results:
                for name in chroms.values():
                    sources[name].append(path)
                observations_count += written
                done += processed
                print(f'[insertion-snps] {done}/{count} insertions; {observations_count} SNP observations', file=sys.stderr, flush=True)
        if pool is not None:
            pool.close()
            pool.join()
        plan = [(compact.chrom_key(name), name, paths) for name, paths in sorted(sources.items())]
        for name in sources:
            db.execute('INSERT OR IGNORE INTO loci VALUES (?,?)', (name, coordinates.length(name)))
        db.commit()
        metadata = ['##fileformat=VCFv4.2']
        if args.append_to:
            metadata, previous_samples = read_header(args.append_to)
            if previous_samples != samples:
                raise ValueError('--append-to SNP VCF sample names/order differ from the SV VCF')
        declared = {}
        for line in metadata:
            if line.startswith('##contig=<'):
                values = vcf._parse_structured_meta(line, 'contig')
                declared[values['ID']] = values.get('length')
        for identifier, length in db.execute('SELECT id,length FROM loci ORDER BY id'):
            if identifier in declared and declared[identifier] != str(length):
                raise ValueError(f'conflicting contig length in --append-to: {identifier}')
            if identifier not in declared:
                metadata.append(f'##contig=<ID={identifier},length={length}>')
        metadata.append('##source=extract_insertion_snps')
        with vcf.open_output_text(str(output)) as out:
            vcf.write_vcf_header(out, metadata, samples, 'small')
            if args.append_to:
                seen = {}
                with open_input(args.append_to) as source:
                    for raw in source:
                        if raw.startswith('#'):
                            continue
                        chrom = raw.split('\t', 1)[0]
                        if chrom not in seen:
                            seen[chrom] = chrom in sources
                        if seen[chrom]:
                            raise ValueError(f'--append-to already contains insertion contig {chrom!r}; refusing duplicate calls')
                        out.write(raw if raw.endswith('\n') else raw + '\n')
            sites = compact.append_insertions(out, {'samples': samples}, plan)
            if any(vcf.checkpoints.file_stamp(path) != stamp for path, stamp in stamps.items()):
                raise ValueError('input VCF changed during extraction')
        print(f'[insertion-snps] wrote {sites} recovered SNP rows to {output}', file=sys.stderr)
    return sites


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('-i', '--input', nargs='+', required=True, help='merged SV VCF, optionally followed by the indel VCF')
    parser.add_argument('-o', '--output', required=True, help='output SNP VCF (plain or gzip)')
    parser.add_argument('--assemblies', help='tab-separated SAMPLE and FASTA path; paths relative to this map')
    parser.add_argument('-a', '--aux-fasta', action='append', default=[], help='indexed graph FASTA for encoded INFO/SEQ dependencies')
    parser.add_argument('-t', '--processes', type=int, default=1)
    parser.add_argument('--append-to', help='copy this existing SNP VCF first and merge insertion contig headers; input is retained')
    parser.add_argument('--tmpdir', help='parent directory for temporary sequence index and worker output')
    parser.add_argument('--force', action='store_true', help='atomically replace an existing output after successful extraction')
    args = parser.parse_args(argv)
    try:
        run(args)
    except (OSError, ValueError, KeyError, RuntimeError, sqlite3.Error) as error:
        print(f'ERROR: {error}', file=sys.stderr)
        return 1
    return 0


if __name__ == '__main__':
    raise SystemExit(main())
