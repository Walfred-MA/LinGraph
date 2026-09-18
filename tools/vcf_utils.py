"""Indexed FASTA access and allele helpers for the grVCF converter.

These helpers preserve the implementation previously imported from the retired
vcf_to_truvari.py script, so the public converter has no archive dependency.
"""
from __future__ import annotations

import os
from collections import OrderedDict
from typing import Optional, Tuple


class IndexedFasta:
    """Minimal random-access FASTA reader backed by a samtools .fai index."""

    def __init__(self, path: str):
        self.path = os.path.abspath(path)
        self.index_path = self.path + ".fai"
        if not os.path.isfile(self.path):
            raise FileNotFoundError(self.path)
        if not os.path.isfile(self.index_path):
            raise FileNotFoundError(
                f"missing FASTA index {self.index_path}; run: samtools faidx {self.path}"
            )
        self.index: "OrderedDict[str, Tuple[int, int, int, int]]" = OrderedDict()
        with open(self.index_path, "rt") as handle:
            for line_number, raw in enumerate(handle, 1):
                if not raw.strip():
                    continue
                fields = raw.rstrip("\n").split("\t")
                if len(fields) < 5:
                    raise ValueError(
                        f"{self.index_path}:{line_number}: expected at least 5 columns"
                    )
                self.index[fields[0]] = tuple(map(int, fields[1:5]))
        if not self.index:
            raise ValueError(f"empty FASTA index: {self.index_path}")
        self.fd = os.open(self.path, os.O_RDONLY)

    def close(self) -> None:
        if self.fd >= 0:
            os.close(self.fd)
            self.fd = -1

    def fetch(self, contig: str, start: int, end: int) -> str:
        if contig not in self.index:
            raise KeyError(f"contig {contig!r} is absent from {self.index_path}")
        length, offset, line_bases, line_width = self.index[contig]
        start = int(start)
        end = int(end)
        if start < 0 or end < start or end > length:
            raise ValueError(
                f"reference interval outside {contig}:0-{length}: {start}-{end}"
            )
        output = bytearray()
        position = start
        while position < end:
            line_index = position // line_bases
            within_line = position % line_bases
            take = min(end - position, line_bases - within_line)
            byte_offset = offset + line_index * line_width + within_line
            output.extend(os.pread(self.fd, take, byte_offset))
            position += take
        return output.decode("ascii").upper()


def parse_info(text: str) -> "OrderedDict[str, Optional[str]]":
    output: "OrderedDict[str, Optional[str]]" = OrderedDict()
    if not text or text == ".":
        return output
    for item in text.split(";"):
        if not item:
            continue
        if "=" in item:
            key, value = item.split("=", 1)
            output[key] = value
        else:
            output[item] = None
    return output


def format_info(info: "OrderedDict[str, Optional[str]]") -> str:
    if not info:
        return "."
    return ";".join(
        key if value is None else f"{key}={value}"
        for key, value in info.items()
    )


def normalize_alleles(pos: int, ref: str, alt: str) -> Tuple[int, str, str]:
    """Minimize a biallelic record without allowing an empty VCF allele."""
    ref = ref.upper()
    alt = alt.upper()
    while len(ref) > 1 and len(alt) > 1 and ref[-1] == alt[-1]:
        ref = ref[:-1]
        alt = alt[:-1]
    while len(ref) > 1 and len(alt) > 1 and ref[0] == alt[0]:
        ref = ref[1:]
        alt = alt[1:]
        pos += 1
    return pos, ref, alt
