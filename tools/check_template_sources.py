#!/usr/bin/env python3
"""Check that local reference templates are their original sequences.

local_reference_templates.fa (OUTPUT/checkpoints/ in graph mode) holds one
record per alternative-locus graph path. Its header names the path's origin,
source=HAPLOTYPE:CONTIG:START-END[+-]. This tool fetches that interval from
the original files and compares it with the record, case-insensitively
(reverse-complemented for '-'):

  - HAPLOTYPE in the assembly list (-q, NAME FASTA per line): that assembly;
  - otherwise a record of a --source-fasta whose own header has the same
    source= (e.g. LinGraph/data/alternatives.fa for imported alternatives)
    or a covering one;
  - otherwise a --source-fasta record named CONTIG (e.g. novel_loci.fa).

When every record is exact, checks against the templates (e.g.
check_merge_lossless.py -r local_reference_templates.fa) are checks against
the original assemblies, alternatives and novel sequences. Read-only.

Example:
  python tools/check_template_sources.py \\
    -t trio_calls/checkpoints/local_reference_templates.fa -q cohort.list \\
    --source-fasta LinGraph/data/alternatives.fa -o template_sources.tsv
"""
import argparse
from collections import Counter, defaultdict
import os
import re
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from check_vcf_lossless import FastaIndex, revcomp  # noqa: E402

SOURCE = re.compile(r'^([^:]+):(.+):(\d+)-(\d+)([+-]?)$')


def records(path):
    """(header fields, sequence) of a FASTA, skipping '#' lines."""
    fields, pieces = None, []
    with open(path) as handle:
        for raw in handle:
            if raw.startswith('#'):
                continue
            if raw.startswith('>'):
                if fields is not None:
                    yield fields, ''.join(pieces)
                fields, pieces = raw[1:].split(), []
            elif raw.strip():
                pieces.append(raw.strip())
    if fields is not None:
        yield fields, ''.join(pieces)


def source_of(fields):
    """(haplotype, contig, start, end, strand) from source=..., else None."""
    for field in fields[1:]:
        if field.startswith('source='):
            match = SOURCE.fullmatch(field[7:])
            if match:
                hap, contig, start, end, strand = match.groups()
                return hap, contig, int(start), int(end), strand or '+'
    return None


def read_assemblies(path):
    base = os.path.dirname(os.path.abspath(path))
    output = {}
    with open(path) as handle:
        for raw in handle:
            fields = raw.split()
            if len(fields) >= 2 and not fields[0].startswith('#'):
                output[fields[0]] = os.path.join(base, fields[1])
    return output


class Sources:
    """Original sequences: assemblies by haplotype, --source-fasta records by
    their own source= interval and by name."""

    def __init__(self, assemblies, fastas):
        self.assemblies = assemblies
        self.readers = {}
        self.by_source = defaultdict(list)      # (hap, contig) -> [(start, end, strand, fasta, name)]
        self.by_name = {}
        for fasta in fastas:
            reader = FastaIndex([fasta])
            self.readers[fasta] = reader
            for fields, _sequence in records(fasta):
                self.by_name.setdefault(fields[0], fasta)
                source = source_of(fields)
                if source:
                    hap, contig, start, end, strand = source
                    self.by_source[(hap, contig)].append((start, end, strand, fasta, fields[0]))

    def reader(self, fasta):
        if fasta not in self.readers:
            self.readers[fasta] = FastaIndex([fasta])
        return self.readers[fasta]

    def fetch(self, hap, contig, start, end, strand):
        """(kind, sequence in the template's orientation) or (kind, None)."""
        if hap in self.assemblies:
            reader = self.reader(self.assemblies[hap])
            if contig in reader and end <= reader.length(contig):
                bases = reader.fetch(contig, start, end)
                return 'assembly', revcomp(bases) if strand == '-' else bases
            return 'assembly', None
        for r0, r1, r_strand, fasta, name in self.by_source.get((hap, contig), ()):
            if r0 <= start and end <= r1:
                reader = self.reader(fasta)
                # The record holds [r0, r1) in its own source orientation.
                if r_strand == '-':
                    bases = revcomp(reader.fetch(name, r1 - end, r1 - start))
                else:
                    bases = reader.fetch(name, start - r0, end - r0)
                return 'source_fasta', revcomp(bases) if strand == '-' else bases
        fasta = self.by_name.get(contig)
        if fasta is not None:
            reader = self.reader(fasta)
            if end <= reader.length(contig):
                bases = reader.fetch(contig, start, end)
                return 'named_record', revcomp(bases) if strand == '-' else bases
        return 'not_found', None


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('-t', '--templates', required=True,
                        help='local_reference_templates.fa')
    parser.add_argument('-q', '--query-paths', required=True,
                        help='NAME FASTA per line (e.g. the cohort list)')
    parser.add_argument('--source-fasta', action='append', default=[], metavar='FASTA',
                        help='indexed FASTA of other originals; repeat (alternatives.fa, '
                             'novel_loci.fa, ...)')
    parser.add_argument('-o', '--report', required=True, help='per-template TSV')
    args = parser.parse_args(argv)
    sources = Sources(read_assemblies(args.query_paths), args.source_fasta)
    counts = Counter()
    with open(args.report, 'w') as out:
        out.write('template\tsource\tkind\tstatus\tlength\tfirst_mismatch\n')
        for fields, sequence in records(args.templates):
            source = source_of(fields)
            if source is None:
                kind, status, first = 'no_source', 'no_source', '.'
            else:
                kind, original = sources.fetch(*source)
                if original is None:
                    status, first = 'source_missing', '.'
                elif original.upper() == sequence.upper():
                    status, first = 'exact', '.'
                else:
                    status = 'mismatch'
                    first = next((str(i) for i, (a, b) in enumerate(zip(original.upper(), sequence.upper()))
                                  if a != b), str(min(len(original), len(sequence))))
            counts[(kind, status)] += 1
            label = '{}:{}:{}-{}{}'.format(*source) if source else '.'
            out.write(f'{fields[0]}\t{label}\t{kind}\t{status}\t{len(sequence)}\t{first}\n')
    print('kind\tstatus\ttemplates')
    for (kind, status), number in sorted(counts.items()):
        print(f'{kind}\t{status}\t{number}')
    exact = sum(n for (_k, status), n in counts.items() if status == 'exact')
    total = sum(counts.values())
    print(f'all_exact\t{"yes" if exact == total else "no"}\t{exact}/{total}')
    return 0 if exact == total else 1


if __name__ == '__main__':
    raise SystemExit(main())
