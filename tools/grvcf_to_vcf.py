#!/usr/bin/env python3
"""Convert LinGraph grVCF representative alleles into reference-based VCF."""
from __future__ import annotations

import argparse
from collections import Counter, OrderedDict
from contextlib import ExitStack
import gzip
import os
from pathlib import Path
import re
import shutil
import sqlite3
import subprocess
import sys
import tempfile

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
sys.path.insert(0, str(ROOT / "scripts"))
from tools.vcf_utils import IndexedFasta, parse_info, format_info, normalize_alleles
from graph_cigar_payloads import query_only_graph_cigar

DNA = re.compile(r'[ACGTNacgtn]+')
OP = re.compile(r'(\d+)([=MXIDHS])([A-Za-z]*)')
GT = re.compile(r'(?:\.|\d+)(?:[|/](?:\.|\d+))*')


def unescape(text):
    # Replace @40 last so a literal escape-looking string stays literal.
    for encoded, plain in (('@0A', '\n'), ('@09', '\t'), ('@2C', ','),
                           ('@3B', ';'), ('@3A', ':'), ('@7C', '|'), ('@40', '@')):
        text = text.replace(encoded, plain)
    return text


def revcomp(sequence):
    return sequence.translate(str.maketrans('ACGTNacgtn', 'TGCANtgcan'))[::-1]


def bases(sequence, label):
    if not DNA.fullmatch(sequence):
        raise ValueError(f'{label}: expected A/C/G/T/N sequence')
    return sequence.upper()


def integer(info, key, default=None):
    value = info.get(key)
    if value in (None, '', '.'):
        return default
    try:
        return int(value)
    except ValueError as error:
        raise ValueError(f'INFO/{key} is not a single integer: {value!r}') from error


def open_input(path):
    return gzip.open(path, 'rt') if str(path).endswith(('.gz', '.bgz')) else open(path)


class Catalog:
    """Indexed FASTAs; the first file defines the output reference contigs."""
    def __init__(self, reference, auxiliary, stack):
        self.readers = {}
        self.contigs = OrderedDict()
        for number, path in enumerate([reference, *auxiliary]):
            reader = IndexedFasta(str(path))
            stack.callback(reader.close)
            for name, entry in reader.index.items():
                if not re.fullmatch(r'[^\s,:<>\[\]]+', name):
                    raise ValueError(f'{path}: FASTA name {name!r} cannot be used as a VCF contig')
                if name in self.readers:
                    raise ValueError(f'duplicate FASTA contig {name!r}; provide unambiguous catalogs')
                if entry[0] < 1 or entry[2] < 1 or entry[3] < entry[2]:
                    raise ValueError(f'{path}.fai: invalid index for {name}')
                self.readers[name] = reader
                if number == 0:
                    self.contigs[name] = entry[0]

    def length(self, name):
        return self.readers[name].index[name][0]

    def fetch(self, name, start, end):
        sequence = self.readers[name].fetch(name, start, end)
        if len(sequence) != end-start:
            raise ValueError(f'{name}:{start}-{end}: truncated FASTA or stale index')
        return bases(sequence, f'{name}:{start}-{end}') if sequence else ''


class Sequences:
    """Resolve encoded alleles, including forward references to other VCF IDs."""
    def __init__(self, database, catalog):
        self.db, self.catalog = database, catalog
        self.active = set()

    def target(self, name):
        cached = self.db.execute('SELECT sequence FROM resolved WHERE id=?', (name,)).fetchone()
        if cached:
            return cached[0]
        row = self.db.execute('SELECT info FROM alleles WHERE id=?', (name,)).fetchone()
        if row is None:
            raise ValueError(f'graph-CIGAR target {name!r} is absent; supply its FASTA with -a')
        if name in self.active:
            raise ValueError(f'circular sequence dependency at {name!r}')
        self.active.add(name)
        try:
            sequence = self.allele(parse_info(row[0]))
        finally:
            self.active.remove(name)
        self.db.execute('INSERT INTO resolved VALUES (?,?)', (name, sequence))
        return sequence

    def interval(self, name, start, end, reverse=False):
        if name in self.catalog.readers:
            if self.db.execute('SELECT 1 FROM alleles WHERE id=?', (name,)).fetchone():
                raise ValueError(f'ambiguous target {name!r}: both FASTA contig and graph allele ID')
            length = self.catalog.length(name)
            if reverse:
                start, end = length-end, length-start
            sequence = self.catalog.fetch(name, start, end)
        else:
            full = self.target(name)
            if reverse:
                start, end = len(full)-end, len(full)-start
            if not 0 <= start <= end <= len(full):
                raise ValueError(f'{name}: graph-CIGAR interval is outside the target')
            sequence = full[start:end]
        return revcomp(sequence) if reverse else sequence

    def decode(self, encoded):
        text = query_only_graph_cigar(unescape(encoded))
        if not text or text[0] not in '<>':
            raise ValueError('encoded sequence must start with > or <')
        output, previous = [], 0
        for chunk in re.finditer(r'([<>])([^<>]+)', text):
            if chunk.start() != previous:
                raise ValueError('malformed graph-CIGAR chunk')
            previous = chunk.end()
            direction, body = chunk.groups()
            target, separator, rest = body.partition(':')
            if separator:
                body = rest
            else:
                target = None
            cursor = position = 0
            operations = list(OP.finditer(body))
            for number, match in enumerate(operations):
                if match.start() != position:
                    raise ValueError('malformed graph-CIGAR operation')
                position = match.end()
                count, operation, payload = int(match[1]), match[2], match[3]
                if payload and (operation not in 'XIS' or len(payload) != count):
                    raise ValueError(f'invalid {operation} payload length')
                if operation == 'H':
                    if number == 0:
                        cursor = count
                    elif number != len(operations)-1:
                        raise ValueError('internal hard clip in graph CIGAR')
                elif operation in '=M':
                    if operation == 'M':
                        raise ValueError('ambiguous M operation; exact =/X or INFO/SEQ is required')
                    if not target:
                        raise ValueError('unnamed match has no sequence target; INFO/SEQ is required')
                    output.append(self.interval(target, cursor, cursor+count, direction == '<'))
                    cursor += count
                elif operation in 'XIS':
                    if count and not payload:
                        raise ValueError(f'{operation} lacks query bases; INFO/SEQ is required')
                    # Inline bases are already in query order, even for < chunks.
                    if payload:
                        output.append(bases(payload, 'graph-CIGAR payload'))
                    if operation == 'X':
                        cursor += count
                elif operation == 'D':
                    cursor += count
            if position != len(body):
                raise ValueError('malformed graph-CIGAR suffix')
        if previous != len(text) or not output:
            raise ValueError('graph CIGAR contains no complete query sequence')
        return ''.join(output)

    def allele(self, info):
        for key in ('SEQ', 'SVINSSEQ', 'EXTENDGRAPHCIGAR', 'CIGAR'):
            text = info.get(key)
            if text in (None, '', '.'):
                continue
            if text.startswith(('<', '>')):
                sequence = self.decode(text)
            elif key in {'SEQ', 'SVINSSEQ'}:
                sequence = bases(text, 'INFO/' + key)
            else:
                raise ValueError(f'INFO/{key} is not a graph CIGAR')
            if info.get('SVTYPE') == 'INS' and integer(info, 'SVLEN') not in (None, len(sequence)):
                raise ValueError('INFO/SVLEN differs from the inserted sequence length')
            return sequence
        raise ValueError('insertion/replacement has no sequence in INFO/SEQ or encoded CIGAR')


def genotypes(fields, sample_count):
    if not sample_count:
        return []
    names = fields[8].split(':')
    if names in (['HSV'], ['HSNP']):
        output = []
        for text in fields[9:]:
            if GT.fullmatch(text):
                output.append(text)
                continue
            # Pipes separate legacy observations, not chromosomes in a genotype.
            values = {unescape(event).split(':', 1)[0] for event in text.split('|')}
            if len(values) != 1 or not GT.fullmatch(next(iter(values))):
                raise ValueError('legacy sample has inconsistent or malformed genotypes')
            output.append(values.pop())
    else:
        if len(names) != len(set(names)) or 'GT' not in names:
            raise ValueError('sample FORMAT must contain one GT field (or legacy HSV/HSNP)')
        index = names.index('GT')
        output = []
        for text in fields[9:]:
            values = text.split(':')
            if text == '.':
                output.append('.')
            elif len(values) <= index or len(values) > len(names):
                raise ValueError('sample column does not match FORMAT')
            else:
                output.append(values[index])
    for value in output:
        if not GT.fullmatch(value):
            raise ValueError(f'invalid GT {value!r}')
    return output


def counts(genotype, alternate_count):
    ac, an, ns = [0]*alternate_count, 0, 0
    for value in genotype:
        called = False
        for allele in re.split(r'[|/]', value):
            if allele == '.':
                continue
            index = int(allele)
            if index > alternate_count:
                raise ValueError(f'GT {value!r} refers to an absent ALT allele')
            an += 1
            called = True
            if index:
                ac[index-1] += 1
        ns += called
    return OrderedDict(AC=','.join(map(str, ac)), AN=str(an),
                       AF=','.join(format(value/an, '.8g') if an else '.' for value in ac), NS=str(ns))


def convert_record(fields, catalog, resolver, sample_count):
    chrom, text_pos, identifier, ref, alt = fields[:5]
    pos = int(text_pos)
    if pos < 0:
        raise ValueError('POS cannot be negative')
    info = parse_info(fields[7])
    genotype = genotypes(fields, sample_count)
    symbolic = alt.startswith('<') and alt.endswith('>')
    if symbolic:
        kind = alt[1:-1].upper()
        declared = info.get('SVTYPE', kind)
        if declared != kind:
            raise ValueError('ALT and INFO/SVTYPE disagree')
        if kind not in {'INS', 'DEL', 'SUB'}:
            raise ValueError(f'unsupported symbolic allele {alt!r}')
        if pos > catalog.contigs[chrom]:
            raise ValueError('POS is beyond the reference contig')
        if pos:
            anchor = catalog.fetch(chrom, pos-1, pos)
            if ref.upper() != anchor and ref.upper() != 'N':
                raise ValueError(f'REF mismatch at {chrom}:{pos}: {ref} != {anchor}')
        end = integer(info, 'END', pos)
        if end < pos:
            raise ValueError('END is before POS')
        if kind == 'INS':
            if end != pos:
                raise ValueError('INS has a nonzero reference span; expected END=POS')
            deletion = ''
        else:
            if end == pos and kind == 'DEL':
                end += abs(integer(info, 'SVLEN', 0))
            if end == pos:
                raise ValueError(f'{kind} has no deleted/replaced interval')
            deletion = catalog.fetch(chrom, pos, end)
        insertion = resolver.allele(info) if kind in {'INS', 'SUB'} else ''
        if kind == 'INS' and integer(info, 'SVLEN') not in (None, len(insertion)):
            raise ValueError('INFO/SVLEN differs from the inserted sequence length')
        if pos:
            ref, alt = anchor + deletion, anchor + insertion
        else:
            # grVCF can use boundary zero. VCF requires a right padding base.
            anchor = catalog.fetch(chrom, end, end+1)
            pos, ref, alt = 1, deletion+anchor, insertion+anchor
        pos, ref, alt = normalize_alleles(pos, ref, alt)
    else:
        if pos < 1:
            raise ValueError('sequence-resolved VCF POS must be at least 1')
        ref = bases(ref, 'REF')
        alternatives = [bases(value, 'ALT') for value in alt.split(',')]
        if catalog.fetch(chrom, pos-1, pos-1+len(ref)) != ref:
            raise ValueError(f'REF does not match the reference at {chrom}:{pos}')
        alt = ','.join(alternatives)
        if len(alternatives) == 1:
            pos, ref, alt = normalize_alleles(pos, ref, alt)
    if ref in alt.split(','):
        raise ValueError('ALT equals REF; refusing to emit a monomorphic variant')
    output_info = counts(genotype, len(alt.split(','))) if sample_count else OrderedDict()
    # Sequence alleles fully describe replacements; do not split a SUB into
    # independent deletion and insertion records or label it by net size.
    row = [chrom, str(pos), identifier, ref, alt, fields[5], fields[6], format_info(output_info)]
    if sample_count:
        row += ['GT', *genotype]
    return row


def load_input(path, database):
    metadata, header = [], None
    with open_input(path) as handle:
        for number, raw in enumerate(handle, 1):
            if raw.startswith('##'):
                if header is not None:
                    raise ValueError(f'{path}:{number}: metadata after #CHROM')
                metadata.append(raw.rstrip('\r\n'))
                continue
            if raw.startswith('#CHROM\t'):
                if header is not None:
                    raise ValueError('duplicate #CHROM header')
                header = raw.rstrip('\r\n').split('\t')
                if header[:8] != ['#CHROM', 'POS', 'ID', 'REF', 'ALT', 'QUAL', 'FILTER', 'INFO']:
                    raise ValueError('invalid VCF column names')
                if len(header) not in (8,) and len(header) < 10:
                    raise ValueError('expected 8 VCF columns or FORMAT plus sample columns')
                if len(header) > 8 and header[8] != 'FORMAT':
                    raise ValueError('ninth VCF column must be FORMAT')
                if len(set(header[9:])) != len(header[9:]):
                    raise ValueError('duplicate VCF sample names')
                continue
            if not raw.strip():
                continue
            if header is None or raw.startswith('#'):
                raise ValueError(f'{path}:{number}: missing or invalid #CHROM header')
            fields = raw.rstrip('\r\n').split('\t')
            if len(fields) != len(header):
                raise ValueError(f'{path}:{number}: expected {len(header)} columns, found {len(fields)}')
            database.execute('INSERT INTO records VALUES (?,?)', (number, raw.rstrip('\r\n')))
            if fields[2] != '.' and fields[4] in {'<INS>', '<SUB>'}:
                try:
                    database.execute('INSERT INTO alleles VALUES (?,?)', (fields[2], fields[7]))
                except sqlite3.IntegrityError as error:
                    raise ValueError(f'{path}:{number}: duplicate graph allele ID {fields[2]!r}') from error
    if header is None:
        raise ValueError('input is missing #CHROM header')
    database.commit()
    return metadata, header


def write_header(handle, metadata, header, catalog, reference, filters):
    handle.write('##fileformat=VCFv4.3\n##source=LinGraph_grvcf_to_vcf\n')
    handle.write('##reference=' + Path(reference).resolve().as_uri() + '\n')
    handle.write('##LinGraphConversion="Representative alleles; graph-only records are in the unplaced grVCF"\n')
    for name, length in catalog.contigs.items():
        handle.write(f'##contig=<ID={name},length={length}>\n')
    declared = set()
    for line in metadata:
        if line.startswith('##FILTER=<ID='):
            declared.add(line.split('ID=', 1)[1].split(',', 1)[0])
            handle.write(line + '\n')
    for name in sorted(filters - declared - {'PASS', '.'}):
        if not re.fullmatch(r'[A-Za-z_][A-Za-z0-9_.-]*', name):
            raise ValueError(f'invalid FILTER ID {name!r}')
        handle.write(f'##FILTER=<ID={name},Description="Filter retained from input grVCF">\n')
    if len(header) > 8:
        handle.write('##INFO=<ID=AC,Number=A,Type=Integer,Description="Alternate allele counts in called genotypes">\n'
                     '##INFO=<ID=AN,Number=1,Type=Integer,Description="Total called allele count">\n'
                     '##INFO=<ID=AF,Number=A,Type=Float,Description="Alternate allele frequencies among called alleles">\n'
                     '##INFO=<ID=NS,Number=1,Type=Integer,Description="Samples with at least one called allele">\n'
                     '##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">\n')
    handle.write('\t'.join(header) + '\n')


def run(args):
    source, output = Path(args.input).resolve(), Path(args.output).resolve()
    unplaced = Path(args.unplaced).resolve() if args.unplaced else output.with_name(
        re.sub(r'\.vcf(?:\.gz)?$', '', output.name) + '.unplaced.gr.vcf')
    inputs = [source, Path(args.reference).resolve(), *map(lambda value: Path(value).resolve(), args.aux_fasta)]
    if output == unplaced or any(path in inputs for path in (output, unplaced)):
        raise ValueError('input and output paths must be different')
    if not output.name.endswith(('.vcf', '.vcf.gz')):
        raise ValueError('--output must end with .vcf or .vcf.gz')
    for path in (output, unplaced, Path(str(output)+'.tbi'), Path(str(output)+'.csi')):
        if path.exists() and not args.force:
            raise FileExistsError(f'output exists: {path}; use --force to replace it')
    bcftools = None
    if args.normalize or str(output).endswith('.gz'):
        bcftools = shutil.which(args.bcftools)
        if not bcftools:
            raise FileNotFoundError('bcftools is required for --normalize or .vcf.gz output')
    output.parent.mkdir(parents=True, exist_ok=True)
    unplaced.parent.mkdir(parents=True, exist_ok=True)
    stats = Counter()
    with ExitStack() as stack, tempfile.TemporaryDirectory(prefix='.grvcf-', dir=output.parent) as tmp:
        work = Path(tmp)
        catalog = Catalog(args.reference, args.aux_fasta, stack)
        database = sqlite3.connect(work / 'records.sqlite')
        stack.callback(database.close)
        database.executescript('''PRAGMA journal_mode=OFF; PRAGMA temp_store=FILE;
            CREATE TABLE records (ordinal INTEGER PRIMARY KEY, text TEXT);
            CREATE TABLE alleles (id TEXT PRIMARY KEY, info TEXT);
            CREATE TABLE resolved (id TEXT PRIMARY KEY, sequence TEXT);
            CREATE TABLE converted (chrom INTEGER, pos INTEGER, ordinal INTEGER, text TEXT);
        ''')
        metadata, header = load_input(source, database)
        resolver = Sequences(database, catalog)
        rank = {name: index for index, name in enumerate(catalog.contigs)}
        sample_count = max(0, len(header)-9)
        filters = set()
        staged_unplaced = work / 'unplaced.gr.vcf'
        with staged_unplaced.open('w') as rejected:
            rejected.write('\n'.join(metadata) + '\n' + '\t'.join(header) + '\n')
            for number, raw in database.execute('SELECT ordinal,text FROM records ORDER BY ordinal'):
                stats['read'] += 1
                fields = raw.split('\t')
                if fields[0] not in catalog.contigs:
                    rejected.write(raw + '\n')
                    stats['unplaced'] += 1
                    continue
                try:
                    row = convert_record(fields, catalog, resolver, sample_count)
                except (ValueError, KeyError) as error:
                    raise ValueError(f'{source}:{number} ({fields[2]}): {error}') from error
                database.execute('INSERT INTO converted VALUES (?,?,?,?)',
                    (rank[row[0]], int(row[1]), number, '\t'.join(row)))
                filters.update(row[6].split(';'))
                stats['converted'] += 1
        database.commit()
        staged = work / 'converted.vcf'
        with staged.open('w') as handle:
            write_header(handle, metadata, header, catalog, args.reference, filters)
            for row, in database.execute('SELECT text FROM converted ORDER BY chrom,pos,ordinal'):
                handle.write(row + '\n')
        if bcftools:
            def command(argv):
                result = subprocess.run([bcftools, *map(str, argv)], capture_output=True, text=True)
                if result.returncode:
                    raise RuntimeError(result.stderr.strip() or 'bcftools failed')
            if args.normalize:
                normalized = work / 'normalized.vcf'
                command(['norm', '-f', args.reference, '-c', 'e', '-Ov', '-o', normalized, staged])
                staged = normalized
            sorted_output = work / ('output.vcf.gz' if str(output).endswith('.gz') else 'output.vcf')
            command(['sort', '-Oz' if str(output).endswith('.gz') else '-Ov', '-o', sorted_output, staged])
            staged = sorted_output
            if str(output).endswith('.gz'):
                command(['index', '--csi', staged])
        # Publish only after every record and optional normalization succeeds.
        with tempfile.NamedTemporaryFile(prefix='.grvcf-unplaced-', dir=unplaced.parent, delete=False) as handle:
            side_tmp = Path(handle.name)
        try:
            shutil.copyfile(staged_unplaced, side_tmp)
            os.replace(side_tmp, unplaced)
        finally:
            side_tmp.unlink(missing_ok=True)
        os.replace(staged, output)
        for suffix in ('.tbi', '.csi'):
            Path(str(output)+suffix).unlink(missing_ok=True)
        if str(output).endswith('.gz'):
            os.replace(str(staged)+'.csi', str(output)+'.csi')
    print(f'[grvcf_to_vcf] read={stats["read"]} converted={stats["converted"]} '
          f'unplaced={stats["unplaced"]}\nVCF: {output}\nUnplaced grVCF: {unplaced}', file=sys.stderr)
    return 0


def build_parser():
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog='Example:\n  python tools/grvcf_to_vcf.py -i calls/cohort.sv.vcf \\\n    -r reference.fa -o calls/cohort.sv.standard.vcf\n\n'
        'Keeps the merged representative alleles and original genotype ploidy.\n'
        'Graph-specific INFO/FORMAT fields are replaced by GT and allele counts.\n'
        'Rows outside the reference FASTA are preserved in a separate grVCF.\n'
        'Use --normalize for bcftools left alignment; .vcf.gz also writes a CSI index.')
    parser.add_argument('-i', '--input', required=True, help='individual or merged grVCF (.vcf or .vcf.gz)')
    parser.add_argument('-r', '--reference', required=True, help='calling reference FASTA, with adjacent .fai')
    parser.add_argument('-o', '--output', required=True, help='standard output .vcf or BGZF .vcf.gz')
    parser.add_argument('-a', '--aux-fasta', action='append', default=[], metavar='FASTA',
                        help='repeatable indexed graph/alternative FASTA for decoding encoded INFO sequences')
    parser.add_argument('--unplaced', help='graph-only grVCF output (default: OUTPUT.unplaced.gr.vcf)')
    parser.add_argument('--normalize', action='store_true', help='left-align with bcftools norm and re-sort')
    parser.add_argument('--bcftools', default='bcftools', help='bcftools executable')
    parser.add_argument('--force', action='store_true', help='replace existing output files after successful conversion')
    return parser


def main(argv=None):
    args = build_parser().parse_args(argv)
    try:
        return run(args)
    except (OSError, ValueError, KeyError, RuntimeError, sqlite3.Error) as error:
        print(f'[grvcf_to_vcf] ERROR: {error}', file=sys.stderr)
        return 1


if __name__ == '__main__':
    raise SystemExit(main())
