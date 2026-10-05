#!/usr/bin/env python3
"""Convert a LinGraph grVCF into a plain VCF, in one streaming pass.

Per row: nested rows (CHROM = a parent row ID or DUP_ template: INS_, DEL_,
SUB_, DUP_) are dropped; an <INS> with plain INFO/SEQ becomes REF/ALT bases
(ALT = REF + SEQ); other symbolic alleles stay symbolic (END/SVLEN); INFO
keeps SVTYPE, END, SVLEN and NSUP; FORMAT keeps GT only. No reference is
read: each row is the merged representative allele, not a base-for-base
haplotype.

The sequence helpers below (Catalog, Sequences, ...) decode graph-encoded
alleles for tools/insertion_coordinates.py.
"""
from __future__ import annotations

import argparse
from collections import Counter
import gzip
import os
from pathlib import Path
import re
import shutil
import subprocess
import sys

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
sys.path.insert(0, str(ROOT / "scripts"))
from tools.vcf_utils import IndexedFasta, parse_info
from graph_cigar_payloads import query_only_graph_cigar

DNA = re.compile(r'[ACGTNacgtn]+')
OP = re.compile(r'(\d+)([=MXIDHS])([A-Za-z]*)')
# CHROM of a nested row: its parent row's ID (graphvcfmerge INS_/DEL_/SUB_
# rows) or a shared full-locus-dup template path (DUP_).
NESTED_CHROM = re.compile(r'(?:INS|DEL|SUB|DUP)_')


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


KEEP_INFO = ('SVTYPE', 'END', 'SVLEN', 'NSUP')
INFO_HEADER = (
    '##INFO=<ID=SVTYPE,Number=1,Type=String,Description="Structural variant type">',
    '##INFO=<ID=END,Number=1,Type=Integer,Description="End position of a symbolic allele">',
    '##INFO=<ID=SVLEN,Number=1,Type=Integer,Description="Insertion positive, deletion negative">',
    '##INFO=<ID=NSUP,Number=1,Type=Integer,Description="Haplotypes carrying the variant">',
    '##ALT=<ID=INS,Description="Insertion (sequence not given)">',
    '##ALT=<ID=DEL,Description="Deletion">',
    '##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">',
)


def contig_id(line):
    match = re.match(r'##contig=<ID=([^,>]+)', line)
    return match.group(1) if match else None


def convert(source, out, stats):
    """Stream one grVCF into a plain VCF (see the module docstring)."""
    contigs, filters, samples = [], [], None
    for raw in source:
        line = raw.rstrip('\r\n')
        if line.startswith('##'):
            name = contig_id(line)
            if name is not None:
                if not NESTED_CHROM.match(name):
                    contigs.append(line)
            elif line.startswith('##FILTER=<'):
                filters.append(line)
            continue
        if line.startswith('#CHROM'):
            columns = line.split('\t')
            samples = len(columns) > 9
            out.write('##fileformat=VCFv4.2\n##source=LinGraph_grvcf_to_vcf\n'
                      '##LinGraphConversion="Representative alleles; nested rows removed"\n')
            for item in (*contigs, *filters, *INFO_HEADER):
                if samples or not item.startswith('##FORMAT'):
                    out.write(item + '\n')
            out.write('\t'.join(columns[:9 if samples else 8] + columns[9:]) + '\n')
            continue
        if not line.strip():
            continue
        if samples is None:
            raise ValueError('missing #CHROM header')
        fields = line.split('\t')
        stats['read'] += 1
        if NESTED_CHROM.match(fields[0]):
            stats['nested_removed'] += 1
            continue
        if fields[1] == '0':
            # Boundary zero has no anchor base; VCF positions start at 1.
            stats['position_zero_removed'] += 1
            continue
        info = parse_info(fields[7])
        if fields[4] == '<INS>' and DNA.fullmatch(info.get('SEQ') or ''):
            fields[4] = (fields[3] + info['SEQ']).upper()
            info.pop('END', None)
            stats['insertions_with_bases'] += 1
        fields[7] = ';'.join(f'{key}={info[key]}' for key in KEEP_INFO
                             if info.get(key) not in (None, '', '.')) or '.'
        if samples:
            keys = fields[8].split(':')
            index = keys.index('GT') if 'GT' in keys else None
            fields[8:] = ['GT'] + [
                (value.split(':')[index] if index is not None and len(value.split(':')) > index
                 else '.') for value in fields[9:]]
        else:
            fields = fields[:8]
        out.write('\t'.join(fields) + '\n')
        stats['written'] += 1


def run(args):
    source, output = Path(args.input).resolve(), Path(args.output).resolve()
    if source == output:
        raise ValueError('input and output paths must be different')
    if output.exists() and not args.force:
        raise FileExistsError(f'output exists: {output}; use --force to replace it')
    compress = output.name.endswith('.gz')
    bgzip = shutil.which('bgzip') if compress else None
    if compress and not bgzip:
        raise FileNotFoundError('bgzip is required for .vcf.gz output')
    stats = Counter()
    staged = output.with_name('.' + output.name + '.tmp')
    with open_input(source) as handle:
        if compress:
            with open(staged, 'wb') as target:
                process = subprocess.Popen([bgzip, '-c'], stdin=subprocess.PIPE, stdout=target,
                                           text=True)
                try:
                    convert(handle, process.stdin, stats)
                finally:
                    process.stdin.close()
                    if process.wait():
                        raise RuntimeError('bgzip failed')
        else:
            with open(staged, 'w') as target:
                convert(handle, target, stats)
    os.replace(staged, output)
    print(f'[grvcf_to_vcf] read={stats["read"]} written={stats["written"]} '
          f'nested_removed={stats["nested_removed"]} '
          f'insertions_with_bases={stats["insertions_with_bases"]} '
          f'position_zero_removed={stats["position_zero_removed"]}\nVCF: {output}', file=sys.stderr)
    return 0


def build_parser():
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog='Example:\n  python3 tools/grvcf_to_vcf.py -i calls/cohort.sv.vcf '
               '-o calls/cohort.sv.standard.vcf.gz')
    parser.add_argument('-i', '--input', required=True, help='individual or merged grVCF (.vcf or .vcf.gz)')
    parser.add_argument('-o', '--output', required=True, help='plain .vcf, or .vcf.gz (bgzip)')
    parser.add_argument('-r', '--reference', help='accepted for older commands; not read')
    parser.add_argument('--force', action='store_true', help='replace an existing output')
    return parser


def main(argv=None):
    args = build_parser().parse_args(argv)
    try:
        return run(args)
    except (OSError, ValueError, RuntimeError) as error:
        print(f'[grvcf_to_vcf] ERROR: {error}', file=sys.stderr)
        return 1


if __name__ == '__main__':
    raise SystemExit(main())
