#!/usr/bin/env python3
"""Rename reference chromosome names in a VCF, including LinGraph grVCFs.

Renamed: the CHROM column and, for names declared by ##contig, the
reference-contig fields of ##contig, ##alternativeLocus (ID),
##referenceCoverage (Chrom) and ##pseudoLinearMapping (Reference). Nested
grVCF rows (CHROM = parent ID INS_/DEL_/SUB_/DUP_...) and record IDs are kept,
so nested rows still point at their parents.

Exactly one renaming mode is required:

  --CHM13fix        NC_060925.1 .. NC_060948.1 -> chr1 .. chr22, chrX, chrY
  --noprefix        drop everything up to the last '#' (HG19#1#chr1 -> chr1)
  --fixtable FILE   two-column table: replace column 1 with column 2

Names that the chosen mode does not match are left unchanged.
"""

from __future__ import annotations

import argparse
import gzip
import re
import sys
from typing import Callable, Dict, Optional, TextIO


CHM13_NAMES: Dict[str, str] = {
    f"NC_{60925 + i:06d}.1": name
    for i, name in enumerate([f"chr{n}" for n in range(1, 23)] + ["chrX", "chrY"])
}

# Reference-contig fields in LinGraph grVCF headers. Only names declared by a
# ##contig line are renamed there, so sample/source contig paths (alt-layer
# ##pseudoLinearMapping Reference, ##alternativeLocus SourceContig, Query)
# keep their names.
HEADER_FIELDS = (
    ("##contig=<", re.compile(r'((?:<|,)ID=)([^,>]+)')),
    ("##alternativeLocus=<", re.compile(r'((?:<|,)ID=)([^,>]+)')),
    ("##referenceCoverage=<", re.compile(r'(,Chrom=")([^"]*)(?=")')),
    ("##pseudoLinearMapping=<", re.compile(r'(,Reference=")([^"]*)(?=:\d+-\d+[+-]")')),
)
# CHROM of a nested grVCF row is its parent row's ID, not a contig.
NESTED_CHROM = re.compile(r'(?:INS|DEL|SUB|DUP)_')


def open_text(path: str, mode: str) -> TextIO:
    if path == "-":
        return sys.stdin if "r" in mode else sys.stdout
    if path.endswith(".gz"):
        return gzip.open(path, mode + "t")
    return open(path, mode + "t")


def load_fixtable(path: str) -> Dict[str, str]:
    """Read OLD NEW pairs (tab or whitespace separated); skip blank and '#' lines."""
    table: Dict[str, str] = {}
    with open_text(path, "r") as handle:
        for line_number, raw in enumerate(handle, 1):
            line = raw.strip()
            if not line or line.startswith("#"):
                continue
            fields = line.split()
            if len(fields) < 2:
                raise ValueError(f"{path}:{line_number}: expected two columns")
            old, new = fields[0], fields[1]
            if old in table and table[old] != new:
                raise ValueError(
                    f"{path}:{line_number}: {old!r} mapped to both "
                    f"{table[old]!r} and {new!r}"
                )
            table[old] = new
    return table


def build_renamer(args: argparse.Namespace) -> Callable[[str], str]:
    if args.CHM13fix:
        return lambda name: CHM13_NAMES.get(name, name)
    if args.noprefix:
        return lambda name: name.rsplit("#", 1)[-1]
    table = load_fixtable(args.fixtable)
    return lambda name: table.get(name, name)


def fix_vcf(source: TextIO, sink: TextIO, rename: Callable[[str], str]) -> Dict[str, str]:
    """Stream the VCF, renaming contigs; return the renames that were applied."""
    applied: Dict[str, str] = {}

    def mapped(name: str) -> str:
        new = rename(name)
        if new != name:
            applied[name] = new
        return new

    # The header is buffered so ##contig IDs are known before any other
    # header line that refers to them.
    header = []
    line = ""
    for line in source:
        if not line.startswith("##"):
            break
        header.append(line)
    declared = {
        match.group(2)
        for raw in header if raw.startswith("##contig=<")
        for match in [HEADER_FIELDS[0][1].search(raw)] if match
    }
    seen_contigs = set()
    for raw in header:
        for prefix, pattern in HEADER_FIELDS:
            if raw.startswith(prefix):
                raw = pattern.sub(
                    lambda m: m.group(1) + (mapped(m.group(2)) if m.group(2) in declared
                                            else m.group(2)),
                    raw, count=1)
                if prefix == "##contig=<":
                    contig = pattern.search(raw).group(2)
                    if contig in seen_contigs:
                        print(f"warning: dropped duplicate ##contig header for {contig}",
                              file=sys.stderr)
                        raw = ""
                    seen_contigs.add(contig)
                break
        sink.write(raw)

    def records():
        if line and not line.startswith("##"):
            yield line
        yield from source

    for record in records():
        if not record.startswith("#"):
            chrom, sep, rest = record.partition("\t")
            if not NESTED_CHROM.match(chrom):
                record = mapped(chrom) + sep + rest
        sink.write(record)
    return applied


def parse_args(argv: Optional[list] = None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("-i", "--input", default="-",
                        help="input VCF (plain or .gz; default stdin)")
    parser.add_argument("-o", "--output", default="-",
                        help="output VCF (.gz suffix writes gzip; default stdout)")
    mode = parser.add_mutually_exclusive_group(required=True)
    mode.add_argument("--CHM13fix", action="store_true",
                      help="rename CHM13 RefSeq accessions NC_060925.1-NC_060948.1 "
                           "to chr1-chr22, chrX, chrY")
    mode.add_argument("--noprefix", action="store_true",
                      help="remove the prefix up to the last '#' (HG19#1#chr1 -> chr1)")
    mode.add_argument("--fixtable", metavar="FILE",
                      help="two-column table; names in column 1 become column 2")
    return parser.parse_args(argv)


def main(argv: Optional[list] = None) -> int:
    args = parse_args(argv)
    try:
        rename = build_renamer(args)
        with open_text(args.input, "r") as source, open_text(args.output, "w") as sink:
            applied = fix_vcf(source, sink, rename)
    except (OSError, ValueError) as error:
        print(f"error: {error}", file=sys.stderr)
        return 1
    print(f"renamed {len(applied)} chromosome name(s)", file=sys.stderr)
    return 0


if __name__ == "__main__":
    sys.exit(main())
