#!/usr/bin/env python3
"""Rewrite of graphcigartoref/refexalign VCF conversion.

Design goals
------------
1. Parse each input alignment row independently.
2. Within a row/sample/assembly-contig, group nearby SV tokens on the assembly
   CIGAR stream using an index/score system.  Background tokens (=, M, X) are
   gaps; D/I/encoded insertions are hits.
3. Across samples, cluster raw SV units on reference coordinates.  Insertion
   clusters use the member with the most uppercase (unmasked) inserted bases
   as their representative; optional recursive minimap2 encoding rewrites the
   other member CIGARs against that representative.
4. SNPs are clustered only by exact CHROM/POS/REF/ALT.
5. Keep memory bounded: stream rows, store compact raw units, cluster SVs in
   chromosome blocks, write a temporary VCF, optionally perform per-line
   insertion encoding, then replace the final output.

This intentionally does not patch the older script.  It reuses compatible input
conventions and conservative helper logic.
"""

from __future__ import annotations

import argparse
import gc
import gzip
import heapq
import math
import multiprocessing as mp
import os
import pickle
import re
import shutil
import sys
import tempfile
import time
import urllib.parse
import bisect
from array import array
from bisect import bisect_left
from collections import defaultdict
from dataclasses import dataclass, field, replace
from statistics import median
from typing import Callable, Dict, Iterable, Iterator, List, Mapping, Optional, Sequence, Tuple

from local_reference_templates import LocalTemplate, read_templates
from graph_cigar_payloads import query_only_graph_cigar

try:
    from ssw import AlignmentMgr  # type: ignore
    try:
        from ssw.alignmentmgr import BitwiseAlignmentFlag  # type: ignore
    except Exception:  # pragma: no cover
        BitwiseAlignmentFlag = None
except Exception:  # pragma: no cover
    AlignmentMgr = None
    BitwiseAlignmentFlag = None

try:
    import mappy  # type: ignore
except Exception:  # pragma: no cover
    mappy = None

QUERY_OPS = {"M", "=", "X", "I", "S"}
REF_OPS = {"M", "=", "X", "D"}
SINGLE_HAPLOTYPE_SAMPLES = {"HG38", "CHM13"}
MAX_INFO_DELSEQ_BP = 1000
INFO_RAW_DELSEQ_MAX_BP = 500
SV_EVENT_FIELDS_LEGACY = (
    "GT", "TYPE", "SIZE", "EXTENDGRAPHCIGAR", "ASSEMBLYCONTIG",
    "QUERYCOORD", "ALLELENAME", "LABEL_H",
)
# Keep the internal SV payload in its historical order so all existing
# coordinate/source indexes remain stable.  TEMPLATEOFFSET is the signed
# displacement from the VCF row POS to this observation's actual reference
# breakpoint.  It is emitted as its own FORMAT field rather than hidden in a
# leading H chunk of EXTENDGRAPHCIGAR.
SV_EVENT_FIELDS = SV_EVENT_FIELDS_LEGACY + ("TEMPLATEOFFSET",)
SNP_EVENT_FIELDS = (
    "GT", "TYPE", "SIZE", "BASE", "ASSEMBLYCONTIG", "QUERYCOORD",
    "ALLELENAME", "LABEL_H",
)
# Keep PA provenance as the final two external FORMAT fields so --noPAtag can
# remove that suffix while retaining TEMPLATEOFFSET.
SV_FORMAT_FIELDS = (
    *SV_EVENT_FIELDS_LEGACY[:6], "TEMPLATEOFFSET",
    *SV_EVENT_FIELDS_LEGACY[6:],
)
SNP_FORMAT_FIELDS = SNP_EVENT_FIELDS
SV_FORMAT = ":".join(SV_FORMAT_FIELDS)
SNP_FORMAT = ":".join(SNP_FORMAT_FIELDS)
SV_FORMAT_BASE_FIELDS = SV_FORMAT_FIELDS[:7]
SNP_FORMAT_BASE_FIELDS = SNP_FORMAT_FIELDS[:-2]
SV_FORMAT_BASE = ":".join(SV_FORMAT_BASE_FIELDS)
SNP_FORMAT_BASE = ":".join(SNP_FORMAT_BASE_FIELDS)
SV_FORMAT_LEGACY = ":".join(SV_EVENT_FIELDS_LEGACY)
SNP_FORMAT_LEGACY = ":".join(SNP_EVENT_FIELDS)
SV_FORMAT_LEGACY_BASE = ":".join(SV_EVENT_FIELDS_LEGACY[:-2])
SNP_FORMAT_LEGACY_BASE = ":".join(SNP_EVENT_FIELDS[:-2])
DEFAULT_BLOCK_GAP = 1_000_000
DEFAULT_FORMAT_PROCESSES = 16
HIDDEN_LINE_SEPARATORS = "\r\n\u2028\u2029"




def mp_context():
    import platform
    if platform.system() == "Darwin":
        # fork is unsafe on macOS (ObjC runtime); use spawn to avoid silent worker crashes
        return mp.get_context("spawn")
    try:
        return mp.get_context("fork")
    except (ValueError, RuntimeError):
        return mp.get_context()


def _freeze_shared_heap():
    """Freeze the parent heap before forking workers.

    Cyclic-GC passes write into object headers, unsharing
    copy-on-write pages in every forked worker.  Freezing the
    already-built shared structures keeps those pages shared;
    collection of objects created afterwards is unchanged.
    """
    gc.freeze()


def short_log_text(text: str, limit: int = 240) -> str:
    text = str(text)
    if len(text) <= limit:
        return text
    half = max(1, limit // 2)
    return text[:half] + f"...[{len(text) - 2 * half} omitted]..." + text[-half:]


def normalize_cli_path(path: Optional[str], option: str) -> Optional[str]:
    """Remove pasted line separators at path edges and reject internal ones."""
    if path is None:
        return None
    cleaned = path.strip(HIDDEN_LINE_SEPARATORS)
    if any(character in cleaned for character in HIDDEN_LINE_SEPARATORS):
        raise ValueError(
            f"{option} contains an embedded newline or Unicode line separator"
        )
    if not cleaned:
        raise ValueError(f"{option} is empty after removing line separators")
    if cleaned != path:
        print(
            f"[graphreftovcf] warning: removed a hidden line separator "
            f"from {option}: {cleaned!r}",
            file=sys.stderr,
        )
    return cleaned

# ---------------------------------------------------------------------------
# Basic utilities
# ---------------------------------------------------------------------------

def revcomp(seq: str) -> str:
    table = str.maketrans("ACGTNacgtn", "TGCANtgcan")
    return seq.translate(table)[::-1]


def strip_match_name(name: str) -> str:
    return re.sub(r"\[.*$", "", name or "")


def open_text(path: str):
    if path.endswith(".gz"):
        return gzip.open(path, "rt")
    return open(path, "r")


def vcf_escape(x) -> str:
    if x is None or x == "":
        return "."
    s = str(x)
    return (
        s.replace("@", "@40")   # @ must be escaped first
         .replace(":", "@3A")
         .replace(";", "@3B")
         .replace(",", "@2C")
         .replace("\t", "@09")
         .replace("\n", "@0A")
    )


def vcf_unescape(x: str) -> str:
    if x is None or x == ".":
        return ""
    s = str(x)
    return (
        s.replace("@0A", "\n")
         .replace("@09", "\t")
         .replace("@2C", ",")
         .replace("@3B", ";")
         .replace("@3A", ":")
         .replace("@40", "@")   # @ must be unescaped last
    )



def get_haplotype(allelename: str) -> str:
    allelename = strip_match_name((allelename or "").strip())
    if not allelename:
        return allelename
    parts = allelename.split("_")
    for i in range(len(parts) - 1, 0, -1):
        if re.fullmatch(r"h[12]", parts[i]):
            return f"{parts[i - 1]}_{parts[i]}"
    return allelename


def get_event_haplotype(query_name: str, label_name: str) -> str:
    """Choose the biological sample for an event, independent of its matrix.

    Current per-sample query IDs have the form
    ``graph_sample_haplotype_lineindex`` while the ownership label can be only
    the graph/matrix name. Prefer a recognizable haplotype in the query ID and
    retain the label-based behavior as a compatibility fallback for older
    graphcigartoref files.
    """
    query_haplotype = get_haplotype(query_name)
    if re.search(r"_h[12]$", query_haplotype or ""):
        return query_haplotype
    label_haplotype = get_haplotype(label_name)
    if label_haplotype:
        return label_haplotype
    return query_haplotype


def first_label_token(label_text: str) -> str:
    label_text = (label_text or "").strip()
    if not label_text or label_text == ".":
        return ""
    return label_text.split(";", 1)[0].strip()


def normalize_requested_haplotypes(columns_list: Sequence[str]) -> List[str]:
    seen = set()
    out: List[str] = []
    for token in columns_list:
        token = token.strip()
        if not token:
            continue
        if re.search(r"_h[12]$", token):
            expanded = [token]
        elif token in SINGLE_HAPLOTYPE_SAMPLES:
            expanded = [f"{token}_h1"]
        else:
            expanded = [f"{token}_h1", f"{token}_h2"]
        for hap in expanded:
            if hap not in seen:
                seen.add(hap)
                out.append(hap)
    return out


# ---------------------------------------------------------------------------
# FASTA and coordinate metadata
# ---------------------------------------------------------------------------

class IndexedFastaReader:
    """Lazy FASTA reader backed by .fai/.faidx."""

    def __init__(self, fasta_path: str, index_path: Optional[str] = None):
        self.fasta_path = fasta_path
        candidates = []
        if index_path:
            candidates.append(index_path)
        candidates.extend([fasta_path + ".fai", fasta_path + ".faidx"])
        root, _ext = os.path.splitext(fasta_path)
        candidates.append(root + ".faidx")
        self.index_path = next((x for x in candidates if x and os.path.exists(x)), None)
        self.index: Dict[str, Tuple[int, int, int, int]] = {}
        if self.index_path is None:
            self._build_index_from_fasta()
        else:
            with open(self.index_path) as handle:
                for raw in handle:
                    if not raw.strip():
                        continue
                    cols = raw.rstrip("\n").split("\t")
                    if len(cols) < 5:
                        continue
                    name = cols[0]
                    length, offset, line_bases, line_width = map(int, cols[1:5])
                    self.index[name] = (length, offset, line_bases, line_width)
        self.fd = os.open(self.fasta_path, os.O_RDONLY)

    def _build_index_from_fasta(self):
        """Build a transient faidx-like index when no .fai/.faidx is present."""
        name: Optional[str] = None
        seq_offset = 0
        length = 0
        line_bases = 0
        line_width = 0

        def finish_record():
            if name is not None:
                self.index[name] = (length, seq_offset, max(1, line_bases), max(1, line_width))

        with open(self.fasta_path, "rb") as handle:
            while True:
                line_start = handle.tell()
                raw = handle.readline()
                if not raw:
                    break
                if raw.startswith(b">"):
                    finish_record()
                    header = raw[1:].strip().split()
                    name = header[0].decode("utf-8") if header else ""
                    seq_offset = handle.tell()
                    length = 0
                    line_bases = 0
                    line_width = 0
                    continue
                if name is None:
                    continue
                bases = raw.rstrip(b"\r\n")
                if not bases:
                    continue
                if line_bases == 0:
                    line_bases = len(bases)
                    line_width = len(raw)
                    seq_offset = line_start
                length += len(bases)
        finish_record()

    def __getstate__(self):
        state = self.__dict__.copy()
        state["fd"] = -1  # file descriptors cannot be pickled; reopen in __setstate__
        return state

    def __setstate__(self, state):
        self.__dict__.update(state)
        self.fd = os.open(self.fasta_path, os.O_RDONLY)

    def __contains__(self, name: str) -> bool:
        return name in self.index

    def length(self, name: str) -> int:
        return self.index[name][0]

    def fetch(self, name: str, start: int, end: int, strand: str = "+") -> str:
        if name not in self.index:
            raise KeyError(f"{name!r} not found in FASTA index {self.index_path}")
        length, offset, line_bases, line_width = self.index[name]
        start = max(0, min(int(start), length))
        end = max(start, min(int(end), length))
        out = bytearray()
        pos = start
        while pos < end:
            line_idx = pos // line_bases
            within = pos % line_bases
            n = min(end - pos, line_bases - within)
            byte_offset = offset + line_idx * line_width + within
            out.extend(os.pread(self.fd, n, byte_offset))
            pos += n
        seq = out.decode("ascii")
        return revcomp(seq) if strand == "-" else seq


class SparseLocalTemplateReader:
    """Expose each interval FASTA record as an independent local contig."""

    def __init__(self, path: str, templates: Sequence[LocalTemplate]):
        self.path = path
        self.reader = IndexedFastaReader(path)
        self.intervals: Dict[str, List[LocalTemplate]] = defaultdict(list)
        self.index: Dict[str, Tuple[int, int, int, int]] = {}
        for template in templates:
            self.intervals[template.output_contig].append(template)
            old_length = self.index.get(
                template.output_contig, (0, 0, 1, 1),
            )[0]
            self.index[template.output_contig] = (
                max(old_length, template.end), 0, 1, 1,
            )
        for rows in self.intervals.values():
            rows.sort(key=lambda row: (
                row.start, row.end, row.graph_name, row.graph_path,
            ))

    def __contains__(self, name: str) -> bool:
        return name in self.index

    def length(self, name: str) -> int:
        return self.index[name][0]

    def fetch(self, name: str, start: int, end: int, strand: str = "+") -> str:
        if name not in self.index:
            raise KeyError(f"{name!r} is absent from local templates")
        length = self.index[name][0]
        start = max(0, min(int(start), length))
        end = max(start, min(int(end), length))
        cursor = start
        pieces = []
        while cursor < end:
            candidates = [
                row for row in self.intervals[name]
                if row.start <= cursor < row.end
            ]
            if not candidates:
                next_start = min(
                    (
                        row.start for row in self.intervals[name]
                        if row.start > cursor
                    ),
                    default=end,
                )
                gap_end = min(end, next_start)
                pieces.append("N" * (gap_end - cursor))
                cursor = gap_end
                continue
            row = min(candidates, key=lambda value: (
                value.end, value.graph_name, value.graph_path,
            ))
            piece_end = min(end, row.end)
            local_start = cursor - row.start
            local_end = piece_end - row.start
            if row.storage_strand == "+":
                piece = self.reader.fetch(
                    row.name, local_start, local_end,
                )
            else:
                length = row.end - row.start
                piece = self.reader.fetch(
                    row.name,
                    length - local_end,
                    length - local_start,
                    "-",
                )
            pieces.append(piece)
            cursor = piece_end
        sequence = "".join(pieces)
        return revcomp(sequence) if strand == "-" else sequence


@dataclass
class SeqRecord:
    name: str
    allelename: str
    contig: str
    start: int
    end: int
    strand: str = "+"
    seqlen: int = 0
    left_ext: str = "0"
    right_ext: str = "0"
    left_neighbor: str = ""
    right_neighbor: str = ""
    reader: Optional[IndexedFastaReader] = None
    fasta_name: Optional[str] = None
    source: str = ""
    orig_contig: Optional[str] = None
    orig_start: Optional[int] = None
    orig_end: Optional[int] = None
    orig_strand: Optional[str] = None
    fasta_offset: int = 0

    def __post_init__(self):
        if self.seqlen <= 0:
            self.seqlen = max(0, int(self.end) - int(self.start))
        if self.fasta_name is None:
            self.fasta_name = self.name
        if self.orig_contig is None:
            self.orig_contig = self.contig
        if self.orig_start is None:
            self.orig_start = int(self.start)
        if self.orig_end is None:
            self.orig_end = int(self.end)
        if self.orig_strand is None:
            self.orig_strand = self.strand

    def fetch_local(self, start: int, end: int, strand: str = "+") -> str:
        start = max(0, min(int(start), self.seqlen))
        end = max(start, min(int(end), self.seqlen))
        if self.reader is None or self.fasta_name not in self.reader.index:
            # Missing FASTA / alias failure means sequence is unavailable, not
            # biologically N.  Return empty so callers do not create fake N gaps.
            return ""
        # FASTA record is in local query orientation. If the metadata strand is
        # '-', the FASTA bases are still row-local orientation for query records.
        seq = self.reader.fetch(self.fasta_name, self.fasta_offset + start, self.fasta_offset + end, "+")
        return revcomp(seq) if strand == "-" else seq

    def genomic_interval_from_local(self, off0: int, off1: int) -> Tuple[int, int]:
        off0 = max(0, min(int(off0), self.seqlen))
        off1 = max(0, min(int(off1), self.seqlen))
        if off1 < off0:
            off0, off1 = off1, off0
        if self.strand == "+":
            return self.start + off0, self.start + off1
        return self.end - off1, self.end - off0

    def genomic_point_from_local(self, off: int) -> int:
        off = max(0, min(int(off), self.seqlen))
        return self.start + off if self.strand == "+" else self.end - off


def parse_coord_token(token: str) -> Optional[Tuple[str, int, int, str]]:
    token = strip_match_name((token or "").strip())
    m = re.match(r"^[<>]?(.+):(\d+)-(\d+)([+-]?)$", token)
    if not m:
        return None
    chrom = m.group(1)
    start = int(m.group(2))
    end = int(m.group(3))
    strand = m.group(4) or "+"
    if end < start:
        start, end = end, start
    return chrom, start, end, strand


def split_extension_neighbor_token(tok: str) -> Tuple[str, str]:
    """Return (neighbor_name, extension_value).

    The -m extension columns can be either the old numeric form ("+500",
    "500", "-500") or the new neighbor-aware form
    "neighborAllele:+500".  Only the last ':' separates the neighbor name
    from the numeric extension value.
    """
    tok = (tok or "").strip()
    if ":" not in tok:
        return "", tok
    neighbor, value = tok.rsplit(":", 1)
    return neighbor.strip(), value.strip()


def extension_token_value(tok: str) -> str:
    _neighbor, value = split_extension_neighbor_token(tok)
    return value


def is_extension_token(tok: str) -> bool:
    tok = extension_token_value(tok)
    return bool(re.fullmatch(r"(?:[+-]?\d+|-inf|\+inf|inf)", tok or ""))


def parse_extend_token(tok: str, max_impute: int) -> Tuple[str, int, int]:
    """Return (sign, raw, effective).  Same convention as the old script.

    -N  shrinks/excludes N bp from that side.
    +N  high-priority extension by min(N, max_impute).
    N   lower-priority gap: min(max(0, N-max_impute), max_impute).
    inf is treated as an excluded side.
    """
    tok = extension_token_value(tok)
    tok = (tok or "").strip()
    if tok in {"", ".", "NA", "None", "none", "null"}:
        return "", 0, 0
    if tok in {"-inf", "+inf", "inf"}:
        return "-", 10**18, 10**18
    sign = ""
    body = tok
    if tok and tok[0] in "+-":
        sign = tok[0]
        body = tok[1:]
    if not body.isdigit():
        return sign, 0, 0
    raw = int(body)
    if sign == "-":
        effective = raw
    elif sign == "+":
        effective = min(raw, max_impute)
    else:
        effective = min(max(0, raw - max_impute), max_impute)
    return sign, raw, effective


def apply_interval_extensions_to_coord(start: int, end: int, strand: str, left_ext="0", right_ext="0", max_impute=10000) -> Tuple[int, int]:
    """Apply -m extensions in absolute genomic-coordinate order.

    The -m extension columns are not orientation-aware:
      * left_ext  always applies to the lower-coordinate/start side.
      * right_ext always applies to the higher-coordinate/end side.

    ``strand`` is intentionally ignored and is kept only for API compatibility
    with older callers.
    """
    del strand
    start = int(start)
    end = int(end)
    lsign, lraw, leff = parse_extend_token(left_ext, max_impute)
    rsign, rraw, reff = parse_extend_token(right_ext, max_impute)

    # Absolute coordinate-left side.
    if lsign == "-":
        start += lraw
    else:
        start -= leff

    # Absolute coordinate-right side.
    if rsign == "-":
        end -= rraw
    else:
        end += reff

    if start < 0:
        start = 0
    if end < start:
        end = start
    return start, end


def add_alias(records: Dict[str, SeqRecord], rec: SeqRecord):
    for name in {rec.name, strip_match_name(rec.name), rec.allelename, strip_match_name(rec.allelename)}:
        if name:
            records[name] = rec
    m_short = re.search(r"([A-Za-z0-9]+_h[12](?:_\d+)?)$", rec.name)
    if m_short:
        records.setdefault(m_short.group(1), rec)
        records.setdefault(strip_match_name(m_short.group(1)), rec)


def read_reference_records(path: str) -> Tuple[Dict[str, SeqRecord], List[str], IndexedFastaReader]:
    reader = IndexedFastaReader(path)
    records: Dict[str, SeqRecord] = {}
    order: List[str] = []
    # The per-sample wrapper's in-memory reader can expose virtual allele
    # views through .index.  When the reference haplotype's assembly IS the
    # reference FASTA, those views must not become reference contigs (and
    # from there ##contig header lines): the reference is only the FASTA's
    # real records. Reference-free loci keep their template sequences and are
    # declared through both ##alternativeLocus and the required ##contig line.
    view_names = getattr(reader, "views", None) or {}
    for name, (length, _offset, _lb, _lw) in reader.index.items():
        if name in view_names:
            continue
        rec = SeqRecord(name=name, allelename=name, contig=name, start=0, end=length, strand="+", seqlen=length, reader=reader, fasta_name=name, source="ref")
        add_alias(records, rec)
        order.append(name)
    return records, order, reader


def add_local_template_reference_records(
    records: Dict[str, SeqRecord],
    order: List[str],
    path: str,
) -> List[LocalTemplate]:
    """Add sparse template intervals as reference-coordinate records."""
    templates = read_templates(path)
    if not templates:
        return []
    local_reader = SparseLocalTemplateReader(path, templates)
    for contig, (length, _offset, _line_bases, _line_width) in (
        local_reader.index.items()
    ):
        if contig in records:
            raise ValueError(
                "local-template output contig collides with requested "
                f"reference FASTA: {contig!r}"
            )
        record = SeqRecord(
            name=contig,
            allelename=contig,
            contig=contig,
            start=0,
            end=length,
            strand="+",
            seqlen=length,
            reader=local_reader,
            fasta_name=contig,
            source="local_template",
        )
        add_alias(records, record)
        order.append(contig)
    return templates


def read_fasta_metadata(path: str, require_coord: bool, max_impute: int, attach_reader: bool = True) -> Tuple[Dict[str, SeqRecord], List[str], Optional[IndexedFastaReader]]:
    reader: Optional[IndexedFastaReader] = None
    if attach_reader:
        try:
            reader = IndexedFastaReader(path)
        except Exception:
            reader = None
    records: Dict[str, SeqRecord] = {}
    order: List[str] = []
    with open_text(path) as fh:
        for raw in fh:
            if not raw.startswith(">"):
                continue
            parts = raw[1:].strip().split()
            if not parts:
                continue
            name = parts[0]
            coord = parse_coord_token(parts[1]) if len(parts) >= 2 else None
            allele = parts[2] if len(parts) >= 3 else name
            if coord is None:
                if require_coord:
                    raise ValueError(f"Header for {name} missing coordinate field")
                length = reader.length(name) if reader is not None and name in reader.index else 0
                coord = (name, 0, length, "+")
            contig, start, end, strand = coord
            seqlen = reader.length(name) if reader is not None and name in reader.index else max(0, end - start)
            rec = SeqRecord(name=name, allelename=allele, contig=contig, start=start, end=end, strand=strand, seqlen=seqlen, reader=reader, fasta_name=name, source="query")
            add_alias(records, rec)
            order.append(name)
    return records, order, reader


def read_liftover_coord_map(path: str, max_impute: int) -> Tuple[Dict[str, SeqRecord], List[str]]:
    records: Dict[str, SeqRecord] = {}
    order: List[str] = []
    with open_text(path) as fh:
        for line_no, raw in enumerate(fh, start=1):
            line = raw.strip()
            if not line or line.startswith("#"):
                continue
            parts = line.split()
            if not parts:
                continue
            name = parts[0]
            part_match = (
                re.fullmatch(r"part_(\d+)", parts[6])
                if len(parts) > 6 else None
            )
            part_index = int(part_match.group(1)) if part_match else 0

            def part_value(value: str) -> str:
                values = value.split(";")
                if part_index and len(values) >= part_index:
                    return values[part_index - 1].strip()
                return value.strip()

            # The original label/query coordinate is now taken explicitly from
            # column 3 of the -m/--coord-map file.  The token format is the
            # same old coordinate format accepted by parse_coord_token(), e.g.
            # contig:start-end+ or contig:start-end-.
            coord = (
                parse_coord_token(part_value(parts[2]))
                if len(parts) >= 3 else None
            )
            if coord is None:
                continue
            left_ext = "0"
            right_ext = "0"
            left_neighbor = ""
            right_neighbor = ""
            if len(parts) >= 3:
                selected_left = part_value(parts[-2])
                selected_right = part_value(parts[-1])
                if (
                    is_extension_token(selected_left)
                    and is_extension_token(selected_right)
                ):
                    left_neighbor, left_ext = split_extension_neighbor_token(
                        selected_left
                    )
                    right_neighbor, right_ext = split_extension_neighbor_token(
                        selected_right
                    )
            contig, start, end, strand = coord
            ext_start, ext_end = apply_interval_extensions_to_coord(start, end, strand, left_ext, right_ext, max_impute)
            rec = SeqRecord(
                name=name,
                allelename=name,
                contig=contig,
                start=ext_start,
                end=ext_end,
                strand=strand,
                seqlen=max(0, ext_end - ext_start),
                left_ext=left_ext,
                right_ext=right_ext,
                left_neighbor=left_neighbor,
                right_neighbor=right_neighbor,
                source="coord_map",
                orig_contig=contig,
                orig_start=start,
                orig_end=end,
                orig_strand=strand,
            )
            add_alias(records, rec)
            order.append(name)

    # Current GenomeLift files store numeric ownership offsets but do not
    # repeat the adjacent block name in those tokens.  Recover the same
    # source/confirmation relationship used by graphvcfmerge from the ordered
    # positive-ownership rows.  Explicit ``neighbor:value`` metadata, when
    # present in older files, always wins.
    by_sample_contig: Dict[Tuple[str, str], List[SeqRecord]] = defaultdict(list)
    seen_records = set()
    for name in order:
        rec = records.get(name)
        if rec is None or id(rec) in seen_records:
            continue
        seen_records.add(id(rec))
        if extension_token_value(rec.left_ext) in {"-inf", "+inf", "inf"}:
            continue
        if extension_token_value(rec.right_ext) in {"-inf", "+inf", "inf"}:
            continue
        by_sample_contig[(get_haplotype(name), rec.orig_contig or rec.contig)].append(rec)
    for neighboring in by_sample_contig.values():
        neighboring.sort(key=lambda rec: (
            int(rec.orig_start if rec.orig_start is not None else rec.start),
            int(rec.orig_end if rec.orig_end is not None else rec.end),
            rec.name,
        ))
        for index in range(1, len(neighboring)):
            previous = neighboring[index - 1]
            current = neighboring[index]
            if not previous.right_neighbor:
                previous.right_neighbor = current.name
            if not current.left_neighbor:
                current.left_neighbor = previous.name
    return records, order


def split_coord_map_paths(value: Optional[str]) -> List[str]:
    """Return comma-ordered coordinate maps, rejecting empty elements."""
    if value is None:
        return []
    raw_paths = value.split(",")
    if any(not path.strip() for path in raw_paths):
        raise ValueError(
            "--coord-map contains an empty comma-separated path"
        )
    return [
        normalize_cli_path(path.strip(), "--coord-map")
        for path in raw_paths
    ]


def read_liftover_coord_maps(
    paths: Sequence[str], max_impute: int,
) -> Tuple[Dict[str, SeqRecord], List[str]]:
    """Overlay coordinate maps from left to right; later PA names win."""
    records: Dict[str, SeqRecord] = {}
    order: List[str] = []
    seen_order = set()
    for path in paths:
        current, current_order = read_liftover_coord_map(
            path, max_impute=max_impute,
        )
        records.update(current)
        for name in current_order:
            if name not in seen_order:
                order.append(name)
                seen_order.add(name)
    return records, order


def read_query_metadata(path: str, require_coord: bool, max_impute: int) -> Tuple[Dict[str, SeqRecord], List[str], Optional[IndexedFastaReader]]:
    with open_text(path) as fh:
        for raw in fh:
            if not raw.strip():
                continue
            if raw.startswith(">"):
                return read_fasta_metadata(path, require_coord=require_coord, max_impute=max_impute, attach_reader=True)
            records, order = read_liftover_coord_map(path, max_impute=max_impute)
            return records, order, None
    return {}, [], None


def merge_coord_metadata(query_records: Dict[str, SeqRecord], coord_records: Dict[str, SeqRecord]) -> Dict[str, SeqRecord]:
    """Return lookup where -m coordinate records override coordinates/extensions.

    Actual sequence readers from query_records are preserved when names match.
    """
    if not coord_records:
        return query_records
    out = dict(query_records)
    for key, crec in coord_records.items():
        qrec = query_records.get(key)
        if qrec is not None and qrec.reader is not None:
            merged = SeqRecord(
                name=qrec.name,
                allelename=crec.allelename or qrec.allelename,
                contig=crec.contig,
                start=crec.start,
                end=crec.end,
                strand=crec.strand,
                seqlen=qrec.seqlen,
                left_ext=crec.left_ext,
                right_ext=crec.right_ext,
                left_neighbor=crec.left_neighbor,
                right_neighbor=crec.right_neighbor,
                reader=qrec.reader,
                fasta_name=qrec.fasta_name,
                source="query+coord_map",
            )
        else:
            merged = crec
        add_alias(out, merged)
    return out


# ---------------------------------------------------------------------------
# CIGAR parsing and stream handling
# ---------------------------------------------------------------------------

@dataclass
class CigarOp:
    op: str
    n: int
    payload: str = ""

_CIGAR_RE = re.compile(r"(\d+)([=MXIDHS])([A-Za-z]*)")

INFINITE_OWNERSHIP_TOKENS = frozenset({"-inf", "+inf", "inf"})
MIN_INFINITE_FALLBACK_ALIGNMENT_SCORE = 1000
VALID_OWNERSHIP_PRIORITY_SCORE = 1_000_000_000_000


def parse_pairwise_ops(body: str) -> List[CigarOp]:
    ops: List[CigarOp] = []
    pos = 0
    n = len(body)
    while pos < n:
        if body[pos] in "><":
            # A named encoded chunk has leaked into a body.  Stop here; caller
            # will handle it as a top-level encoded insertion.
            break
        m = _CIGAR_RE.match(body, pos)
        if not m:
            # Ignore separators or legacy punctuation rather than crashing the row.
            pos += 1
            continue
        length = int(m.group(1))
        op = m.group(2)
        payload = m.group(3) or ""
        if op == "M":
            op = "="
        ops.append(CigarOp(op=op, n=length, payload=payload))
        pos = m.end()
    return ops


def validate_graph_cigar_syntax(
    encoded: str, *, context: str = "graph CIGAR",
) -> None:
    """Require complete, unambiguous token coverage of a graph CIGAR."""
    text = encoded or ""
    if not text or text[0] not in "><":
        raise ValueError(f"{context}: CIGAR must begin with > or <: {text!r}")
    cursor = 0
    for chunk in re.finditer(r"([<>])([^<>]*)", text):
        if chunk.start() != cursor:
            raise ValueError(f"{context}: malformed graph CIGAR {text!r}")
        value = chunk.group(2)
        if ":" in value:
            name, body = value.split(":", 1)
            if not name:
                raise ValueError(f"{context}: empty graph target in {text!r}")
        else:
            body = value
        if not body:
            raise ValueError(f"{context}: empty alignment chunk in {text!r}")
        body_cursor = 0
        while body_cursor < len(body):
            match = _CIGAR_RE.match(body, body_cursor)
            if match is None:
                raise ValueError(
                    f"{context}: invalid alignment text at offset "
                    f"{chunk.start(2) + body_cursor}: {text!r}"
                )
            length = int(match.group(1))
            operation = match.group(2)
            payload = match.group(3) or ""
            # A named 0H-only chunk is an explicit coordinate-zero anchor
            # before an encoded insertion. It consumes no query/reference.
            zero_origin_anchor = operation == 'H' and body == '0H' and ':' in value
            if length <= 0 and not zero_origin_anchor:
                raise ValueError(f"{context}: zero-length CIGAR operation")
            if payload and operation not in {"I", "X"}:
                raise ValueError(
                    f"{context}: unexpected payload on {operation} operation"
                )
            if payload and len(payload) != length:
                raise ValueError(
                    f"{context}: {operation} payload has {len(payload)} bases, "
                    f"expected {length}"
                )
            body_cursor = match.end()
        cursor = chunk.end()
    if cursor != len(text):
        raise ValueError(f"{context}: malformed graph CIGAR {text!r}")


def drop_empty_encoded_segments(encoded: str) -> str:
    """Discard later named H-only chunks without changing the main declaration.

    Older saved alignments can contain an empty encoded segment such as
    ``>template:100H0H``. It has no query bases or aligned reference bases.
    Keep every other chunk verbatim so malformed real operations still fail
    validation and pathless main continuations retain their original meaning.
    """
    chunks = list(re.finditer(r"([<>])([^<>]*)", encoded or ""))
    if not chunks or chunks[0].start() != 0 or ":" not in chunks[0].group(2):
        return encoded
    return "".join(
        chunk.group(0) for index, chunk in enumerate(chunks)
        if index == 0 or re.fullmatch(r"[^:<>]+:(?:[0-9]+H)+", chunk.group(2)) is None
    )


def format_op(op: CigarOp, include_payload: bool = True) -> str:
    if op.op in {"=", "D", "H", "S"} or not include_payload:
        return f"{op.n}{op.op}"
    return f"{op.n}{op.op}{op.payload}"


def iter_top_level_chunks(encoded: str) -> Iterator[Tuple[str, str]]:
    if not encoded:
        return
    i = 0
    n = len(encoded)
    while i < n:
        if encoded[i] not in "><":
            i += 1
            continue
        direction = encoded[i]
        j = i + 1
        while j < n and encoded[j] not in "><":
            j += 1
        yield direction, encoded[i + 1:j]
        i = j


def main_path_affine_score_components(
    chunks: Sequence[Tuple[str, str]],
) -> Tuple[int, int, int, int]:
    """Return M, X, gap-open, and gap-extension counts for the main path.

    The first named graph-CIGAR chunk and subsequent pathless chunks form the
    main alignment.  Every later named chunk is a secondary/encoded traversal:
    its complete query span counts as insertion sequence on the main axis and
    its internal matches cannot improve the main-path score.  Hard clips are
    coordinate padding and do not interrupt an adjacent affine gap.
    """
    scoring_ops: List[CigarOp] = []
    main_seen = False
    for _direction, text in chunks:
        if ":" in text:
            _name, body = text.split(":", 1)
            if not main_seen:
                main_seen = True
                scoring_ops.extend(parse_pairwise_ops(body))
            else:
                query_span = sum(
                    operation.n
                    for operation in parse_pairwise_ops(body)
                    if operation.op in QUERY_OPS
                )
                if query_span > 0:
                    scoring_ops.append(CigarOp(
                        op="I", n=query_span, payload="",
                    ))
        elif main_seen:
            scoring_ops.extend(parse_pairwise_ops(text))

    matches = 0
    mismatches = 0
    gap_opens = 0
    gap_extensions = 0
    active_gap = ""
    for operation in scoring_ops:
        if operation.n <= 0:
            continue
        if operation.op in {"I", "D"}:
            if active_gap == operation.op:
                gap_extensions += operation.n
            else:
                gap_opens += 1
                gap_extensions += max(0, operation.n - 1)
                active_gap = operation.op
            continue
        if operation.op == "H":
            continue
        active_gap = ""
        if operation.op in {"=", "M"}:
            matches += operation.n
        elif operation.op == "X":
            mismatches += operation.n
    return matches, mismatches, gap_opens, gap_extensions


def infinite_fallback_alignment_score(encoded: str) -> int:
    """Score a complete graph alignment for ``-inf`` ownership fallback.

    The score requested for choosing between overlapping unowned blocks is:

        matched bases - 4 * number of insertion/deletion operations

    Mismatch bases are neither matches nor insertion/deletion events. Indel
    length is deliberately not charged a second time: one I or D CIGAR token
    is one event regardless of the number of bases that token consumes.
    """
    score = 0
    for _direction, text in iter_top_level_chunks(encoded or ""):
        body = text.split(":", 1)[1] if ":" in text else text
        for operation in parse_pairwise_ops(body):
            if operation.op == "=":
                score += operation.n
            elif operation.op in {"I", "D"}:
                score -= 4
    return score


def graph_cigar_operation_count(chunks: Sequence[Tuple[str, str]]) -> int:
    """Count positive non-H operations across main and encoded segments."""
    return sum(
        op.n > 0 and op.op != 'H'
        for _direction, text in chunks
        for op in parse_pairwise_ops(text.split(':', 1)[-1])
    )


def terminal_indel_indexes(
    ops: Sequence[CigarOp], *, row_operation_count: Optional[int] = None,
) -> set:
    """Skip edge I/D only when the whole row has more than one non-H op."""
    if row_operation_count is None:
        row_operation_count = sum(op.n > 0 and op.op != 'H' for op in ops)
    if row_operation_count <= 1:
        return set()
    aligned = [index for index, op in enumerate(ops) if op.n > 0 and op.op in {'=', 'M', 'X'}]
    if not aligned:
        return {index for index, op in enumerate(ops) if op.op in {'I', 'D'}}
    return {index for index, op in enumerate(ops)
            if op.op in {'I', 'D'} and (index < aligned[0] or index > aligned[-1])}


def terminal_main_indel_indexes(
    chunks: Sequence[Tuple[str, str]], *, row_operation_count: Optional[int] = None,
) -> set:
    if row_operation_count is None:
        row_operation_count = graph_cigar_operation_count(chunks)
    ops = []
    for index, (_direction, text) in enumerate(chunks):
        if index == 0 or ':' not in text:
            body = text.split(':', 1)[-1]
            ops.extend(op for op in parse_pairwise_ops(body) if op.n > 0 and op.op != 'H')
    return terminal_indel_indexes(ops, row_operation_count=row_operation_count)


def terminal_main_insertion_indexes(chunks: Sequence[Tuple[str, str]]) -> set:
    # Compatibility name; callers now skip both types of terminal indel.
    return terminal_main_indel_indexes(chunks)


def encoded_query_len(text: str) -> int:
    total = 0
    for n_s, op in re.findall(r"(\d+)([=MXIDHS])", text or ""):
        if op in QUERY_OPS:
            total += int(n_s)
    return total


def is_named_encoded_segment(text: str) -> bool:
    return bool(text and text[0] in "><" and ":" in text)


def split_first_graph_segment(encoded: str) -> Tuple[str, str, str]:
    for direction, text in iter_top_level_chunks(encoded):
        if ":" in text:
            name, body = text.split(":", 1)
            return direction, name, body
    return ">", "", encoded


def segment_clip_lengths(ops: Sequence[CigarOp]) -> Tuple[int, int]:
    lead_h = ops[0].n if ops and ops[0].op == "H" else 0
    trail_h = ops[-1].n if ops and ops[-1].op == "H" else 0
    return lead_h, trail_h


def ref_len_for_rec(rec: Optional[SeqRecord]) -> int:
    return rec.seqlen if rec is not None else 0


def segment_interval_to_forward(ref_rec: SeqRecord, direction: str, lead_h: int, trail_h: int, off0: int, off1: int) -> Tuple[int, int]:
    off0 = int(off0)
    off1 = int(off1)
    if off1 < off0:
        off0, off1 = off1, off0
    L = ref_len_for_rec(ref_rec)
    if direction == "<":
        end = L - lead_h
        return end - off1, end - off0
    start = lead_h
    return start + off0, start + off1


def segment_point_to_forward(ref_rec: SeqRecord, direction: str, lead_h: int, trail_h: int, off: int) -> int:
    L = ref_len_for_rec(ref_rec)
    if direction == "<":
        return L - lead_h - int(off)
    return lead_h + int(off)


def local_ref_to_genomic(ref_rec: SeqRecord, f0: int, f1: int) -> Tuple[int, int]:
    # Reference/path records may themselves have a genomic strand.  Convert
    # forward local offsets on that record to genomic coordinates.
    return ref_rec.genomic_interval_from_local(f0, f1)


def local_ref_point_to_genomic(ref_rec: SeqRecord, off: int) -> int:
    return ref_rec.genomic_point_from_local(off)


def resolve_record(name: str, records: Dict[str, SeqRecord]) -> Optional[SeqRecord]:
    if not name:
        return None
    return records.get(name) or records.get(strip_match_name(name))


def map_query_interval_to_label_h(query_rec: SeqRecord, label_rec: SeqRecord, q0: int, q1: int) -> Optional[Tuple[int, int]]:
    """Map a query/object-local interval onto the label/ownership local axis.

    Variant parsing walks the -i query/object CIGAR, so q0/q1 are local to
    query_rec.  Ownership labels, extensions, and uncertain regions belong to
    the -i label coordinate plus -m neighbor/offsite metadata.  This helper
    converts through assembly coordinates and clips to the label interval.
    """
    if query_rec.contig != label_rec.contig:
        return None
    q0 = int(q0)
    q1 = int(q1)
    if q1 < q0:
        q0, q1 = q1, q0
    if q0 == q1:
        g = query_rec.genomic_point_from_local(q0)
        if g < label_rec.start or g > label_rec.end:
            return None
        # Label-local coordinates are absolute genomic offsets from the
        # lower-coordinate/start side.  Do not flip by label_rec.strand.
        local = g - label_rec.start
        local = max(0, min(label_rec.seqlen, int(local)))
        return local, local

    g0, g1 = query_rec.genomic_interval_from_local(q0, q1)
    if g1 < g0:
        g0, g1 = g1, g0
    ov0 = max(g0, label_rec.start)
    ov1 = min(g1, label_rec.end)
    if ov1 <= ov0:
        return None
    # Label-local coordinates are absolute genomic offsets from the
    # lower-coordinate/start side.  Do not flip by label_rec.strand.
    l0 = ov0 - label_rec.start
    l1 = ov1 - label_rec.start
    l0 = max(0, min(label_rec.seqlen, int(l0)))
    l1 = max(0, min(label_rec.seqlen, int(l1)))
    if l1 < l0:
        l0, l1 = l1, l0
    return l0, l1


def object_interval_to_genomic_span(query_rec: SeqRecord, q0: int, q1: int) -> Tuple[int, int]:
    """Return the un-clipped genomic span of a query/object-local interval."""
    q0 = int(q0)
    q1 = int(q1)
    if q1 < q0:
        q0, q1 = q1, q0
    if q1 == q0:
        g = query_rec.genomic_point_from_local(q0)
        return g, g
    g0, g1 = query_rec.genomic_interval_from_local(q0, q1)
    if g1 < g0:
        g0, g1 = g1, g0
    return g0, g1


def format_label_h_value(value: int) -> str:
    """Format one label-boundary distance as an H token.

    Values are deliberately not clipped: negative and zero distances are valid.
    """
    return f"{int(value)}H"


def format_label_h_point(label_rec: SeqRecord, assembly_pos: int) -> str:
    """Return SNP label-boundary distance as {leftsize}H.

    Breakpoint/label-H distances are absolute genomic-coordinate distances.
    Strand is intentionally ignored.
    """
    assembly_pos = int(assembly_pos)
    leftsize = assembly_pos - int(label_rec.start)
    return format_label_h_value(leftsize)


def format_label_h_span(query_rec: SeqRecord, label_rec: SeqRecord, q0: int, q1: int) -> str:
    """Return SV/indel label-boundary distances as {leftsize}H_{rightsize}H.

    Breakpoint/label-H distances are absolute genomic-coordinate distances:
      leftsize  = SVstart - regionstart
      rightsize = regionend - SVend

    Strand is intentionally ignored.  Coordinates are computed from the
    un-clipped assembly endpoints, so values can be negative or zero without
    changing label-assignment eligibility.
    """
    if query_rec.contig != label_rec.contig:
        return "."
    sv_start, sv_end = object_interval_to_genomic_span(query_rec, q0, q1)
    region_start = int(label_rec.start)
    region_end = int(label_rec.end)
    leftsize = int(sv_start) - region_start
    rightsize = region_end - int(sv_end)
    return f"{format_label_h_value(leftsize)}_{format_label_h_value(rightsize)}"

_ACTIVE_INTERVAL_OWNERSHIP: Dict[
    str, Tuple[str, Tuple[Tuple[int, int], ...]]
] = {}


def owned_object_local_intervals(
    query_rec: SeqRecord,
    base_label: str,
    q0: int,
    q1: int,
) -> Optional[List[Tuple[int, int]]]:
    """Return owned portions of a query-local interval.

    ``None`` means no interval-sweep metadata exists for this label and asks
    callers to preserve legacy ownership behavior.  An empty list means the
    interval is explicitly owned by a higher-priority aligned row.
    """
    ownership = (
        _ACTIVE_INTERVAL_OWNERSHIP.get(base_label)
        or _ACTIVE_INTERVAL_OWNERSHIP.get(strip_match_name(base_label))
    )
    if ownership is None:
        return None
    contig, genomic_intervals = ownership
    if query_rec.contig != contig:
        return []
    q0 = int(q0)
    q1 = int(q1)
    if q1 < q0:
        q0, q1 = q1, q0
    is_point = q1 == q0
    if is_point:
        genomic_point = query_rec.genomic_point_from_local(q0)
        for start, end in genomic_intervals:
            # A deletion consumes no query base.  Keep a boundary deletion on
            # either touching owned interval so split deletions from adjacent
            # rows can still be cross-row validated and merged.
            if start <= genomic_point <= end:
                return [(q0, q1)]
        return []
    genomic_start, genomic_end = query_rec.genomic_interval_from_local(q0, q1)
    if genomic_end < genomic_start:
        genomic_start, genomic_end = genomic_end, genomic_start
    accepted_genomic = intersect_intervals(
        [(genomic_start, genomic_end)], genomic_intervals,
    )
    accepted_local: List[Tuple[int, int]] = []
    for start, end in accepted_genomic:
        if query_rec.strand == "+":
            local_start = start - query_rec.start
            local_end = end - query_rec.start
        else:
            local_start = query_rec.end - end
            local_end = query_rec.end - start
        local_start = max(q0, min(q1, local_start))
        local_end = max(q0, min(q1, local_end))
        if local_end > local_start:
            accepted_local.append((local_start, local_end))
    return merge_intervals(accepted_local)


def object_interval_is_owned(
    query_rec: SeqRecord,
    base_label: str,
    q0: int,
    q1: int,
) -> Optional[bool]:
    intervals = owned_object_local_intervals(query_rec, base_label, q0, q1)
    return None if intervals is None else bool(intervals)


def label_for_object_interval(query_rec: SeqRecord, label_rec: SeqRecord, base_label: str, q0: int, q1: int, is_snp: bool, max_impute: int, bp_uncertain: int) -> str:
    owned = object_interval_is_owned(query_rec, base_label, q0, q1)
    if owned is False:
        return ""
    mapped = map_query_interval_to_label_h(query_rec, label_rec, q0, q1)
    if mapped is None:
        return ""
    l0, l1 = mapped
    label = label_for_query_interval(
        label_rec, base_label, l0, l1, is_snp, max_impute, bp_uncertain,
    )
    # A negative-overlap fallback is deliberately outside the legacy
    # core/extension label regions.  The global interval sweep has already
    # proved that no aligned primary owner covers it, so admit that interval as
    # this label's ordinary ownership instead of rejecting the whole block.
    if not label and owned is True:
        return base_label
    return label


# ---------------------------------------------------------------------------
# Score merge by index
# ---------------------------------------------------------------------------

def score_merge_by_index(
    items: Sequence[dict],
    isgap: Callable[[dict], bool],
    gap_break: int = 100,
    gap_penalty: int = 4,
    hit_score: Optional[Callable[[dict], int]] = None,
    gap_size: Optional[Callable[[dict], int]] = None,
) -> List[List[dict]]:
    """Return groups of non-gap items using two-pass AND score logic.

    `isgap(x)` identifies background items.  Hits are variant tokens.  A merge
    boundary survives only if both left-to-right and right-to-left passes agree.
    """
    if not items:
        return []
    if hit_score is None:
        hit_score = lambda x: 1
    if gap_size is None:
        gap_size = lambda x: 1

    def one_pass(seq: Sequence[dict]) -> List[bool]:
        flags = [False] * len(seq)
        score = 0
        consecutive_gap = 0
        active = False
        for i, x in enumerate(seq):
            if isgap(x):
                g = max(1, int(gap_size(x)))
                consecutive_gap += g
                score -= gap_penalty * g
                if score < 0 or consecutive_gap > gap_break:
                    score = 0
                    consecutive_gap = 0
                    active = False
                elif active and i > 0:
                    flags[i] = True
            else:
                if active and i > 0:
                    flags[i] = True
                score += max(1, int(hit_score(x)))
                consecutive_gap = 0
                active = True
        return flags

    n = len(items)
    fwd = one_pass(items)
    rev = one_pass(list(reversed(items)))
    bwd = [False] * n
    for j in range(1, n):
        if rev[j]:
            bwd[n - j] = True

    merge = [False] * n
    for i in range(1, n):
        merge[i] = fwd[i] and bwd[i]

    groups: List[List[dict]] = []
    cur: List[dict] = []
    for i, x in enumerate(items):
        if i > 0 and not merge[i]:
            if cur:
                groups.append(cur)
            cur = []
        if not isgap(x):
            cur.append(x)
    if cur:
        groups.append(cur)
    return groups


def mark_scored_interval_replacements(
    groups: Sequence[Sequence[dict]],
) -> List[List[dict]]:
    """Mark only groups accepted by the original within-path score merger.

    ``score_merge_by_index`` enforces both the score test and the hard
    ``within_gap_break`` limit (100 bp by default).  A multi-operation group
    accepted there is emitted as one full interval replacement.  This helper
    must never perform additional distance-based linking: doing so would
    bypass the hard alignment-gap limit.
    """
    output: List[List[dict]] = []
    for group in groups:
        ordered = sorted(group, key=lambda item: int(item["idx"]))
        if len(ordered) > 1:
            ordered = [
                dict(item, _interval_replacement=True) for item in ordered
            ]
        output.append(ordered)
    return output


def split_colocated_opposite_indels(
    groups: Sequence[Sequence[dict]],
) -> List[List[dict]]:
    """Keep directly adjacent insertion/deletion components separate.

    The within-line score merger normally joins a zero-background ``D/I`` or
    ``I/D`` pair into one interval replacement.  When requested, split only
    that exact opposite-axis boundary so the two components become separate
    VCF rows.  Same-axis components and operations separated by alignment
    background retain the existing merger behavior.
    """
    output: List[List[dict]] = []

    def axis(item: dict) -> str:
        operation = item.get("op")
        if operation == "D":
            return "D"
        if operation in {"I", "ENCODED_INS"}:
            return "I"
        return ""

    for group in groups:
        current: List[dict] = []
        previous = None
        for item in sorted(group, key=lambda value: int(value["idx"])):
            if previous is not None:
                previous_axis = axis(previous)
                item_axis = axis(item)
                if (
                    int(item["idx"]) == int(previous["idx"]) + 1
                    and previous_axis
                    and item_axis
                    and previous_axis != item_axis
                ):
                    if current:
                        output.append(current)
                    current = []
            current.append(item)
            previous = item
        if current:
            output.append(current)
    return output


# ---------------------------------------------------------------------------
# Raw event containers
# ---------------------------------------------------------------------------

@dataclass
class SVUnit:
    id: int
    query: str
    sample: str
    allelename: str
    ref_contig: str
    ref_start: int
    ref_end: int
    qry_contig: str
    qry_start: int
    qry_end: int
    qry_strand: str
    label: str
    ins_size: int
    del_size: int
    label_h: str = "."
    tokens: List[dict] = field(default_factory=list)
    cigar: str = ""
    raw_ins_seq: str = ""
    encoded: bool = False
    type_prefix: str = ""
    line_no: int = 0
    qry_local_start: int = -1
    qry_local_end: int = -1
    # A pure deletion at a row edge is held until all rows are available.  It
    # is rescued only if it joins a deletion from another block/line.
    edge_candidate: bool = False
    cross_row_merged: bool = False
    main_match_bases: int = 0
    complete_query_length: int = 0
    main_mismatch_bases: int = 0
    main_gap_opens: int = 0
    main_gap_extensions: int = 0
    alternative_regions: Optional[Tuple[Tuple[str, int, int], ...]] = None

    @property
    def max_size(self) -> int:
        return max(self.ins_size, self.del_size)

    @property
    def svtype(self) -> str:
        if self.ins_size > 0 and self.del_size > 0 and self.ins_size == self.del_size:
            return "SUB"
        if self.ins_size > self.del_size:
            return "INS"
        return "DEL"

    @property
    def signed_size(self) -> int:
        if self.svtype == "INS":
            return self.ins_size
        if self.svtype == "DEL":
            return -self.del_size
        return self.ins_size - self.del_size

    @property
    def main_alignment_score(self) -> int:
        return (
            self.main_match_bases
            - 4 * self.main_mismatch_bases
            - 4 * self.main_gap_opens
            - self.main_gap_extensions
        )


@dataclass
class SNPUnit:
    query: str
    sample: str
    allelename: str
    chrom: str
    pos0: int
    ref_base: str
    alt_base: str
    qry_contig: str
    qry_pos: int
    qry_strand: str
    label: str
    label_h: str = "."
    type_prefix: str = ""
    qry_local_start: int = -1
    qry_local_end: int = -1


@dataclass(frozen=True)
class MainPathAlignmentEvidence:
    query: str
    sample: str
    allele_name: str
    qry_contig: str
    main_match_bases: int
    complete_query_length: int
    line_no: int
    main_mismatch_bases: int = 0
    main_gap_opens: int = 0
    main_gap_extensions: int = 0

    @property
    def alignment_score(self) -> int:
        return (
            self.main_match_bases
            - 4 * self.main_mismatch_bases
            - 4 * self.main_gap_opens
            - self.main_gap_extensions
        )


@dataclass
class RowParseResult:
    svs: List[SVUnit]
    snps: List[SNPUnit]
    coverage_spans: List[Tuple[str, str, int, int]]
    filtered_spans: List[Tuple[str, str, int, int]]
    seen_queries: List[str]
    row_filter_reason: str = ""
    edge_filtered_events: int = 0
    alternative_filtered_events: int = 0
    main_path_evidence: List[MainPathAlignmentEvidence] = field(
        default_factory=list,
    )


# ---------------------------------------------------------------------------
# Labeling regions
# ---------------------------------------------------------------------------

def label_regions(rec: SeqRecord, max_impute: int, bp_uncertain: int) -> Tuple[
    List[Tuple[int, int]],
    List[Tuple[int, int]],
    List[Tuple[int, int]],
    List[Tuple[int, int]],
    List[Tuple[int, int]],
    List[Tuple[int, int]],
    List[Tuple[int, int]],
    List[Tuple[int, int]],
    List[Tuple[int, int]],
]:
    """Return local intervals for label assignment.

    Returns:
      core_clean, extension_clean, neighbor_overlap, left_neighbor_only,
      right_neighbor_only, self_left_uncertain, self_right_uncertain,
      core_original, extension_original

    Core/extension use the old ownership rules.  Neighbor-side uncertain
    windows are removed from core/extension and reported with leading-star
    neighbor labels.  Self-side uncertain windows are returned separately so
    they can be labeled as selflabel*neighbor... before the plain core/ext
    checks run.
    """
    seqlen = rec.seqlen
    lsign, lraw, leff = parse_extend_token(rec.left_ext, max_impute)
    rsign, rraw, reff = parse_extend_token(rec.right_ext, max_impute)

    left_boundary = 0
    right_boundary = seqlen
    left_ext_interval: Optional[Tuple[int, int]] = None
    right_ext_interval: Optional[Tuple[int, int]] = None

    if lsign == "-":
        left_boundary = min(seqlen, lraw)
    else:
        le = min(seqlen, leff)
        left_boundary = le
        if le > 0:
            left_ext_interval = (0, le)

    if rsign == "-":
        right_boundary = max(0, seqlen - rraw)
    else:
        reff2 = min(seqlen, reff)
        right_boundary = max(0, seqlen - reff2)
        if reff2 > 0:
            right_ext_interval = (right_boundary, seqlen)

    if right_boundary < left_boundary:
        right_boundary = left_boundary

    core = [(left_boundary, right_boundary)] if right_boundary > left_boundary else []
    ext: List[Tuple[int, int]] = []
    if left_ext_interval and left_ext_interval[1] > left_ext_interval[0]:
        ext.append(left_ext_interval)
    if right_ext_interval and right_ext_interval[1] > right_ext_interval[0]:
        ext.append(right_ext_interval)

    # One-sided breakpoint uncertainty model.
    #
    # The self-side interval is the portion of the uncertain region that
    # overlaps this label/self.  The neighbor-side interval is the portion that
    # does not overlap self and therefore belongs to the neighboring label.
    #
    #   self_left_uncertain      = [left_break, left_break + bp_uncertain]
    #   neighbor_left_uncertain  = [left_break - 2*bp_uncertain, left_break]
    #   self_right_uncertain     = [right_break - bp_uncertain, right_break]
    #   neighbor_right_uncertain = [right_break, right_break + 2*bp_uncertain]
    left_self_uncertain: List[Tuple[int, int]] = []
    right_self_uncertain: List[Tuple[int, int]] = []
    left_neighbor_uncertain: List[Tuple[int, int]] = []
    right_neighbor_uncertain: List[Tuple[int, int]] = []
    if bp_uncertain > 0:
        # These endpoints already are the effective ownership boundaries after
        # applying each positive extension or negative overlap above.  Do not
        # apply lraw/rraw a second time here.  The old double application moved
        # a large-overlap boundary deep into the block and falsely marked
        # ordinary internal variants as cross-graph events.
        left_endpoint = left_boundary
        right_endpoint = right_boundary

        if 0 <= lraw < 10**12:
            left_break = left_endpoint

            a = max(0, left_break)
            b = min(seqlen, left_break + bp_uncertain)
            if b > a:
                left_self_uncertain.append((a, b))

            a = max(0, left_break - 2 * bp_uncertain)
            b = min(seqlen, left_break)
            if b > a:
                left_neighbor_uncertain.append((a, b))

        if 0 <= rraw < 10**12:
            right_break = right_endpoint

            a = max(0, right_break - bp_uncertain)
            b = min(seqlen, right_break)
            if b > a:
                right_self_uncertain.append((a, b))

            a = max(0, right_break)
            b = min(seqlen, right_break + 2 * bp_uncertain)
            if b > a:
                right_neighbor_uncertain.append((a, b))

    # Explicit self windows are checked before neighbor windows so self has
    # priority whenever an uncertain region touches both sides.
    left_self_uncertain = merge_intervals(left_self_uncertain)
    right_self_uncertain = merge_intervals(right_self_uncertain)

    left_neighbor_uncertain = merge_intervals(left_neighbor_uncertain)
    right_neighbor_uncertain = merge_intervals(right_neighbor_uncertain)
    neighbor_all = merge_intervals(left_neighbor_uncertain + right_neighbor_uncertain)
    neighbor_overlap = intersect_intervals(left_neighbor_uncertain, right_neighbor_uncertain)
    left_only = subtract_intervals(left_neighbor_uncertain, neighbor_overlap)
    right_only = subtract_intervals(right_neighbor_uncertain, neighbor_overlap)

    return (
        subtract_intervals(core, neighbor_all),
        subtract_intervals(ext, neighbor_all),
        neighbor_overlap,
        left_only,
        right_only,
        left_self_uncertain,
        right_self_uncertain,
        core,
        ext,
    )


def old_region_label(base_label: str, q0: int, q1: int, core: Sequence[Tuple[int, int]], ext: Sequence[Tuple[int, int]]) -> str:
    """Return the old self label for a region: self first, then self+."""
    if overlaps_any(q0, q1, core):
        return base_label
    if overlaps_any(q0, q1, ext):
        return base_label + "+"
    return base_label


def _uncertain_neighbor_suffix(rec: SeqRecord, left_hit: bool, right_hit: bool) -> List[str]:
    """Return neighbor labels in left/right order for uncertain windows."""
    neighbors: List[str] = []
    if left_hit and rec.left_neighbor:
        neighbors.append(rec.left_neighbor)
    if right_hit and rec.right_neighbor:
        neighbors.append(rec.right_neighbor)
    return neighbors


def uncertain_self_label(base_label: str, rec: SeqRecord, left_hit: bool, right_hit: bool, self_label: Optional[str] = None) -> str:
    """Label for the self-overlapping side of an uncertain window.

    Self-side uncertainty keeps the self label first, then appends any neighbor
    labels after '*', for example: labelA*neighborL or
    labelA*neighborL*neighborR.  If neighbor metadata is unavailable, keep a
    trailing '*' so the call is still marked as breakpoint-uncertain.
    """
    # Star-separated components are cross-window identity keys and must be
    # exact allele names. ``self_label`` may end in '+' to describe ordinary
    # extension ownership; that marker is not part of allele identity.
    label = base_label
    neighbors = _uncertain_neighbor_suffix(rec, left_hit, right_hit)
    if not neighbors:
        return label + "*"
    return label + "*" + "*".join(neighbors)


def uncertain_neighbor_label(base_label: str, rec: SeqRecord, left_hit: bool, right_hit: bool, self_label: Optional[str] = None) -> str:
    """Label for the neighbor-only side of an uncertain window.

    Neighbor-only uncertainty gets a leading '*', followed by the self label
    and then any neighboring labels: *labelA*neighborL.
    """
    label = base_label
    neighbors = _uncertain_neighbor_suffix(rec, left_hit, right_hit)
    if not neighbors:
        return "*" + label
    return "*" + label + "*" + "*".join(neighbors)

def merge_intervals(intervals: Sequence[Tuple[int, int]]) -> List[Tuple[int, int]]:
    vals = sorted((min(a, b), max(a, b)) for a, b in intervals if b > a)
    if not vals:
        return []
    out = [vals[0]]
    for a, b in vals[1:]:
        la, lb = out[-1]
        if a <= lb:
            out[-1] = (la, max(lb, b))
        else:
            out.append((a, b))
    return out


def subtract_intervals(intervals: Sequence[Tuple[int, int]], masks: Sequence[Tuple[int, int]]) -> List[Tuple[int, int]]:
    intervals = merge_intervals(intervals)
    masks = merge_intervals(masks)
    out: List[Tuple[int, int]] = []
    mask_index = 0
    for start, end in intervals:
        cursor = start
        while (
            mask_index < len(masks)
            and masks[mask_index][1] <= cursor
        ):
            mask_index += 1
        scan_index = mask_index
        while scan_index < len(masks) and masks[scan_index][0] < end:
            mask_start, mask_end = masks[scan_index]
            if cursor < mask_start:
                out.append((cursor, min(mask_start, end)))
            cursor = max(cursor, mask_end)
            if cursor >= end:
                break
            scan_index += 1
        if cursor < end:
            out.append((cursor, end))
        # Input intervals are merged and disjoint, so masks passed completely
        # by this interval can never affect a later interval. Keep a spanning
        # mask at scan_index for the next interval.
        mask_index = scan_index
        while (
            mask_index < len(masks)
            and masks[mask_index][1] <= end
        ):
            mask_index += 1
    return out


def intersect_intervals(a_intervals: Sequence[Tuple[int, int]], b_intervals: Sequence[Tuple[int, int]]) -> List[Tuple[int, int]]:
    a_vals = merge_intervals(a_intervals)
    b_vals = merge_intervals(b_intervals)
    out: List[Tuple[int, int]] = []
    i = 0
    j = 0
    while i < len(a_vals) and j < len(b_vals):
        a0, a1 = a_vals[i]
        b0, b1 = b_vals[j]
        x0 = max(a0, b0)
        x1 = min(a1, b1)
        if x1 > x0:
            out.append((x0, x1))
        if a1 < b1:
            i += 1
        else:
            j += 1
    return out


def interval_overlap(a0: int, a1: int, b0: int, b1: int) -> int:
    """Return positive overlap length for half-open intervals, else 0.

    Coordinates are normalized so reversed intervals do not accidentally become
    truthy/falsey for the wrong reason.
    """
    if a1 < a0:
        a0, a1 = a1, a0
    if b1 < b0:
        b0, b1 = b1, b0
    return max(0, min(a1, b1) - max(a0, b0))


def iter_interval_pairs(intervals):
    """Yield (start, end) pairs from flat or accidentally nested interval lists.

    Accepts both [(0, 10), (20, 30)] and [[(0, 10)], [(20, 30)]].
    This avoids false failures when callers pass grouped interval lists.
    """
    if intervals is None:
        return
    for item in intervals:
        if item is None:
            continue
        if (
            isinstance(item, (tuple, list))
            and len(item) == 2
            and isinstance(item[0], (int, float))
            and isinstance(item[1], (int, float))
        ):
            yield int(item[0]), int(item[1])
        elif isinstance(item, (tuple, list)):
            yield from iter_interval_pairs(item)


def overlaps_any(q0: int, q1: int, intervals: Sequence[Tuple[int, int]]) -> bool:
    if q1 <= q0:
        q1 = q0 + 1
    return any(interval_overlap(q0, q1, a, b) > 0 for a, b in iter_interval_pairs(intervals))


def label_for_query_interval(rec: SeqRecord, base_label: str, q0: int, q1: int, is_snp: bool, max_impute: int, bp_uncertain: int) -> str:
    if q1 <= q0:
        q0 = max(0, q0 - 1)
        q1 = min(rec.seqlen, q1 + 1)
    # SNP behavior stays old: self / self+ only, no breakpoint-uncertain '*'.
    core, ext, neighbor_overlap, left_neighbor, right_neighbor, left_self, right_self, core_old, ext_old = label_regions(
        rec,
        max_impute=max_impute,
        bp_uncertain=0 if is_snp else bp_uncertain,
    )
    # Priority:
    #   1. self-side uncertain windows: selflabel*neighbor...
    #   2. old core/extension labels outside uncertain windows
    #   3. neighbor-only uncertain windows: *selflabel*neighbor...
    if not is_snp:
        self_label = old_region_label(base_label, q0, q1, core_old, ext_old)
        left_self_hit = overlaps_any(q0, q1, left_self)
        right_self_hit = overlaps_any(q0, q1, right_self)
        if left_self_hit or right_self_hit:
            return uncertain_self_label(base_label, rec, left_self_hit, right_self_hit, self_label)

    if overlaps_any(q0, q1, core):
        return base_label
    if overlaps_any(q0, q1, ext):
        return base_label + "+"

    if not is_snp:
        self_label = old_region_label(base_label, q0, q1, core_old, ext_old)
        if overlaps_any(q0, q1, neighbor_overlap):
            return uncertain_neighbor_label(base_label, rec, True, True, self_label)
        if overlaps_any(q0, q1, left_neighbor):
            return uncertain_neighbor_label(base_label, rec, True, False, self_label)
        if overlaps_any(q0, q1, right_neighbor):
            return uncertain_neighbor_label(base_label, rec, False, True, self_label)
    return ""


# ---------------------------------------------------------------------------
# Per-line parser
# ---------------------------------------------------------------------------

def parse_graphcigartoref_row(line: str, line_no: int) -> Optional[dict]:
    line = line.strip()
    if not line or line.startswith("#"):
        return None
    parts = line.split()
    if len(parts) < 4:
        raise ValueError(f"line {line_no}: expected >=4 columns, got {len(parts)}: {line[:200]!r}")
    row = {
        "line_no": line_no,
        "query": parts[0],
        "query_coord": parts[1],
        "reference": ".",
        "reference_coord": ".",
        "labels": ".",
        "label_coords": ".",
        "encoded": "",
        "raw": line,
    }
    if (
        len(parts) >= 7
        and parse_coord_token(parts[1]) is None
        and parse_coord_token(parts[2].split(";", 1)[0]) is not None
    ):
        # graphcigartoref_persample v2 keeps the complete comparison and its
        # full whole-alignment boundaries.  Individual ownership names and
        # coordinates are resolved from -m/--coord-map below.
        row.update({
            "schema": "persample_v2",
            "reference": parts[1],
            "query_coord": parts[2],
            "reference_coord": parts[3],
            "labels": parts[0].split("|outer_strand=", 1)[0],
            "label_coords": ".",
            "query_alignment_regions": parts[4],
            "reference_alignment_regions": parts[5],
            "encoded": normalize_encoded_cigar_text(parts[6]),
        })
    elif len(parts) >= 7:
        row["schema"] = "legacy"
        row["reference"] = parts[2]
        row["reference_coord"] = parts[3]
        row["labels"] = parts[4]
        row["label_coords"] = parts[5]
        row["encoded"] = normalize_encoded_cigar_text(parts[6])
    elif len(parts) >= 6:
        row["reference"] = parts[2]
        row["reference_coord"] = parts[3]
        row["labels"] = parts[4]
        row["encoded"] = normalize_encoded_cigar_text(parts[5])
    else:
        row["labels"] = parts[2]
        row["encoded"] = normalize_encoded_cigar_text(parts[3])
    # Accept older gap/refinement output at the input boundary, then keep one
    # query-only convention for validation and variant calling. No alignment
    # or genomic coordinates change when dropping the reference half of X.
    try:
        row["encoded"] = query_only_graph_cigar(row["encoded"])
    except ValueError as error:
        raise ValueError(
            f"line {line_no} query {row['query']!r}: {error}"
        ) from error
    from alternative_intervals import read_tag
    row['alternative_intervals'] = read_tag(parts)
    return row


def expand_graphcigartoref_ownership_rows(row: dict) -> List[dict]:
    """Expand one v2 merged comparison into individual ownership views.

    The expensive graph-to-reference CIGAR is intentionally retained intact.
    Each logical view changes only the ownership allele and selects the
    corresponding full query-alignment interval, allowing VCF calls to be
    clipped against the original GenomeLift block independently.
    """
    if row.get("schema") != "persample_v2":
        return [row]
    if str(row.get("query", "")).startswith("DEL"):
        deletion = dict(row)
        deletion["labels"] = "DEL"
        return [deletion]

    labels = [
        value.strip()
        for value in str(row.get("labels", "")).split(";")
        if value.strip()
    ]
    if not labels:
        return [row]
    query_regions = [
        value.strip()
        for value in str(row.get("query_alignment_regions", "")).split(";")
        if value.strip()
    ]
    expanded = []
    for index, label in enumerate(labels):
        logical = dict(row)
        logical["labels"] = label
        if len(query_regions) == len(labels):
            logical["query_alignment_region"] = query_regions[index]
        else:
            logical["query_alignment_region"] = str(
                row.get("query_alignment_regions", "")
            )
        expanded.append(logical)
    return expanded


def _ownership_coordinate_intervals(
    rec: SeqRecord,
    max_impute: int,
) -> Tuple[str, Tuple[int, int], List[Tuple[int, int]]]:
    """Return ``(contig, primary, fallback)`` ownership coordinates.

    ``primary`` is the ordinary core plus positive/unsigned extension after
    applying the established ownership cap.  ``fallback`` is sequence removed
    only by a negative overlap token.  The latter remains alignable, but is
    accepted only where no primary interval actually aligned.
    """
    contig = rec.orig_contig or rec.contig
    start = int(rec.orig_start if rec.orig_start is not None else rec.start)
    end = int(rec.orig_end if rec.orig_end is not None else rec.end)
    strand = rec.orig_strand or rec.strand
    left_sign, _left_raw, left_effective = parse_extend_token(
        rec.left_ext, max_impute,
    )
    right_sign, _right_raw, right_effective = parse_extend_token(
        rec.right_ext, max_impute,
    )
    full_start = (
        start if left_sign == "-" else max(0, start - left_effective)
    )
    full_end = end if right_sign == "-" else end + right_effective
    primary_start, primary_end = apply_interval_extensions_to_coord(
        start, end, strand, rec.left_ext, rec.right_ext, max_impute,
    )
    primary_start = max(full_start, primary_start)
    primary_end = min(full_end, primary_end)
    if primary_end < primary_start:
        primary_end = primary_start
    fallback = subtract_intervals(
        [(full_start, full_end)],
        [(primary_start, primary_end)] if primary_end > primary_start else [],
    )
    return contig, (primary_start, primary_end), fallback


def _has_infinite_ownership_token(rec: SeqRecord) -> bool:
    """Return true when GenomeLift says the row has no owned interval."""
    return any(
        extension_token_value(token) in INFINITE_OWNERSHIP_TOKENS
        for token in (rec.left_ext, rec.right_ext)
    )


def _original_record_interval(rec: SeqRecord) -> Tuple[str, int, int]:
    contig = rec.orig_contig or rec.contig
    start = int(rec.orig_start if rec.orig_start is not None else rec.start)
    end = int(rec.orig_end if rec.orig_end is not None else rec.end)
    if end < start:
        start, end = end, start
    return contig, start, end


def _valid_region_union_by_sample_contig(
    label_lookup: Mapping[str, SeqRecord],
    max_impute: int,
) -> Dict[Tuple[str, str], List[Tuple[int, int]]]:
    """Return the union of every finite GenomeLift ownership region.

    This deliberately uses GenomeLift ownership coordinates whether or not a
    graph alignment retained that entire interval. A valid region therefore
    always excludes an overlapping ``-inf`` fallback candidate.
    """
    regions: Dict[Tuple[str, str], List[Tuple[int, int]]] = defaultdict(list)
    seen_records = set()
    for rec in label_lookup.values():
        record_id = id(rec)
        if record_id in seen_records:
            continue
        seen_records.add(record_id)
        if _has_infinite_ownership_token(rec):
            continue
        contig, primary, _fallback = _ownership_coordinate_intervals(
            rec, max_impute,
        )
        start, end = primary
        if end <= start:
            continue
        regions[(get_haplotype(rec.name), contig)].append((start, end))
    return {
        sample_contig: merge_intervals(intervals)
        for sample_contig, intervals in regions.items()
    }


def resolve_infinite_fallback_candidates(
    candidates_by_sample_contig: Mapping[
        Tuple[str, str], Sequence[Tuple[int, int, int, str, int]]
    ],
    valid_regions_by_sample_contig: Mapping[
        Tuple[str, str], Sequence[Tuple[int, int]]
    ],
) -> Tuple[Dict[str, List[Tuple[int, int]]], int, int]:
    """Assign each uncovered ``-inf`` base to exactly one qualified label.

    Candidate tuples are ``(start, end, total_score, label, line_number)``.
    Valid GenomeLift intervals participate in the same coordinate sweep with
    score 1,000,000,000,000, so they always outrank ``-inf`` candidates. Only
    winning segments with a smaller score are recorded. Label and input line
    are deterministic tie-breakers for equal-scoring ``-inf`` rows. Runtime is
    O((N + M) log(N + M)) for N candidates and M valid-region intervals.

    Returns ``(intervals_by_label, selected_bases, competing_bases)``.
    """
    selected_by_label: Dict[str, List[Tuple[int, int]]] = defaultdict(list)
    selected_bases = 0
    competing_bases = 0

    for sample_contig, raw_candidates in candidates_by_sample_contig.items():
        candidates: List[Tuple[int, int, int, str, int, bool]] = []
        for start, end, score, label, line_no in raw_candidates:
            if (
                score >= MIN_INFINITE_FALLBACK_ALIGNMENT_SCORE
                and score < VALID_OWNERSHIP_PRIORITY_SCORE
                and end > start
            ):
                candidates.append((
                    start, end, score, label, line_no, False,
                ))
        if not candidates:  # No recordable -inf candidate in this bucket.
            continue

        for valid_index, (start, end) in enumerate(
            valid_regions_by_sample_contig.get(sample_contig, ())
        ):
            if end > start:
                candidates.append((
                    start,
                    end,
                    VALID_OWNERSHIP_PRIORITY_SCORE,
                    "",
                    valid_index,
                    True,
                ))

        # Events are grouped by coordinate during the scan. The active state
        # before processing events at X describes the half-open segment ending
        # at X; after processing it describes the segment starting at X.
        CANDIDATE_END = 0
        CANDIDATE_START = 1
        events: List[Tuple[int, int, int]] = []
        for candidate_id, candidate in enumerate(candidates):
            start, end = candidate[:2]
            events.append((start, CANDIDATE_START, candidate_id))
            events.append((end, CANDIDATE_END, candidate_id))
        events.sort(key=lambda value: value[:2])

        active_ids = set()
        active_label_counts: Dict[str, int] = defaultdict(int)
        score_heap: List[Tuple[int, int, str, int, int]] = []
        previous_coordinate: Optional[int] = None
        event_index = 0

        while event_index < len(events):
            coordinate = events[event_index][0]
            if (
                previous_coordinate is not None
                and coordinate > previous_coordinate
            ):
                while score_heap and score_heap[0][4] not in active_ids:
                    heapq.heappop(score_heap)
                if score_heap:
                    winner_id = score_heap[0][4]
                    winner = candidates[winner_id]
                    winner_score = winner[2]
                    winner_label = winner[3]
                    winner_is_valid = winner[5]
                    # Valid regions participate in ranking but are masks only;
                    # record exclusively the lower-scoring -inf winner.
                    if (
                        not winner_is_valid
                        and winner_score < VALID_OWNERSHIP_PRIORITY_SCORE
                    ):
                        intervals = selected_by_label[winner_label]
                        if (
                            intervals
                            and intervals[-1][1] == previous_coordinate
                        ):
                            intervals[-1] = (intervals[-1][0], coordinate)
                        else:
                            intervals.append((
                                previous_coordinate, coordinate,
                            ))
                        segment_size = coordinate - previous_coordinate
                        selected_bases += segment_size
                        if len(active_label_counts) > 1:
                            competing_bases += segment_size

            while (
                event_index < len(events)
                and events[event_index][0] == coordinate
            ):
                _coordinate, event_type, candidate_id = events[event_index]
                if event_type == CANDIDATE_START:
                    (
                        _start, _end, score, label, line_no, is_valid,
                    ) = candidates[candidate_id]
                    active_ids.add(candidate_id)
                    if not is_valid:
                        active_label_counts[label] += 1
                    heapq.heappush(
                        score_heap,
                        (
                            -score,
                            0 if is_valid else 1,
                            label,
                            line_no,
                            candidate_id,
                        ),
                    )
                else:
                    active_ids.discard(candidate_id)
                    label = candidates[candidate_id][3]
                    is_valid = candidates[candidate_id][5]
                    if not is_valid:
                        remaining = active_label_counts[label] - 1
                        if remaining:
                            active_label_counts[label] = remaining
                        else:
                            del active_label_counts[label]
                event_index += 1
            previous_coordinate = coordinate

    return dict(selected_by_label), selected_bases, competing_bases


def build_effective_interval_ownership(
    graphcigar_path: str,
    label_lookup: Dict[str, SeqRecord],
    max_impute: int,
) -> Dict[str, Tuple[str, Tuple[Tuple[int, int], ...]]]:
    """Resolve finite and ``-inf`` ownership from emitted alignments.

    Finite calls retain their score-clipped graph-alignment ownership. For
    deciding whether ``-inf`` may be used, however, every valid GenomeLift
    region participates in one combined sweep at priority 1,000,000,000,000.
    An ``-inf`` row requires total alignment score >= 1000 and is recorded only
    where it wins below that sentinel priority.
    """
    primary_by_label: Dict[str, List[Tuple[int, int]]] = defaultdict(list)
    fallback_by_label: Dict[str, List[Tuple[int, int]]] = defaultdict(list)
    infinite_by_sample_contig: Dict[
        Tuple[str, str], List[Tuple[int, int, int, str, int]]
    ] = defaultdict(list)
    label_contig: Dict[str, str] = {}
    label_sample: Dict[str, str] = {}
    infinite_below_score = 0
    valid_regions_by_sample_contig = _valid_region_union_by_sample_contig(
        label_lookup, max_impute,
    )
    with open_text(graphcigar_path) as handle:
        for line_no, raw in enumerate(handle, start=1):
            row = parse_graphcigartoref_row(raw, line_no)
            if row is None:
                continue
            query_coord = parse_coord_token(row.get("query_coord", ""))
            if query_coord is None:
                continue
            query_contig, query_start, query_end, _query_strand = query_coord
            if query_end <= query_start:
                continue
            for ownership_row in expand_graphcigartoref_ownership_rows(row):
                label = first_label_token(ownership_row.get("labels", ""))
                if not label or label == "DEL":
                    continue
                rec = resolve_record(label, label_lookup)
                if rec is None:
                    continue
                if _has_infinite_ownership_token(rec):
                    contig, original_start, original_end = (
                        _original_record_interval(rec)
                    )
                    if contig != query_contig or original_end <= original_start:
                        continue
                    aligned = intersect_intervals(
                        [(query_start, query_end)],
                        [(original_start, original_end)],
                    )
                    if not aligned:
                        continue
                    score = infinite_fallback_alignment_score(
                        ownership_row.get("encoded", "")
                    )
                    if score < MIN_INFINITE_FALLBACK_ALIGNMENT_SCORE:
                        infinite_below_score += 1
                        continue
                    sample = get_haplotype(label)
                    for aligned_start, aligned_end in aligned:
                        infinite_by_sample_contig[(sample, contig)].append((
                            aligned_start,
                            aligned_end,
                            score,
                            label,
                            int(row.get("line_no", line_no)),
                        ))
                    label_contig[label] = contig
                    label_sample[label] = sample
                    continue
                contig, primary, fallback = _ownership_coordinate_intervals(
                    rec, max_impute,
                )
                if contig != query_contig:
                    continue
                full_parts = (
                    ([primary] if primary[1] > primary[0] else []) + fallback
                )
                aligned = intersect_intervals(
                    [(query_start, query_end)], full_parts,
                )
                if not aligned:
                    continue
                if label in label_contig and label_contig[label] != contig:
                    raise ValueError(
                        f"ownership label {label!r} occurs on multiple contigs"
                    )
                label_contig[label] = contig
                label_sample[label] = get_haplotype(label)
                if primary[1] > primary[0]:
                    primary_by_label[label].extend(
                        intersect_intervals(aligned, [primary])
                    )
                if fallback:
                    fallback_by_label[label].extend(
                        intersect_intervals(aligned, fallback)
                    )

    primary_by_sample_contig: Dict[
        Tuple[str, str], List[Tuple[int, int]]
    ] = defaultdict(list)
    for label, intervals in primary_by_label.items():
        primary_by_sample_contig[(
            label_sample[label], label_contig[label],
        )].extend(intervals)
    for key, intervals in list(primary_by_sample_contig.items()):
        primary_by_sample_contig[key] = merge_intervals(intervals)

    effective: Dict[str, Tuple[str, Tuple[Tuple[int, int], ...]]] = {}
    finite_by_label: Dict[str, List[Tuple[int, int]]] = {}
    labels = set(primary_by_label) | set(fallback_by_label)
    primary_bases = 0
    fallback_bases = 0
    for label in labels:
        primary = merge_intervals(primary_by_label.get(label, ()))
        fallback = subtract_intervals(
            merge_intervals(fallback_by_label.get(label, ())),
            primary_by_sample_contig.get((
                label_sample[label], label_contig[label],
            ), ()),
        )
        intervals = merge_intervals(primary + fallback)
        if intervals:
            finite_by_label[label] = intervals
            effective[label] = (label_contig[label], tuple(intervals))
        primary_bases += sum(end - start for start, end in primary)
        fallback_bases += sum(end - start for start, end in fallback)

    infinite_by_label, infinite_bases, competing_infinite_bases = (
        resolve_infinite_fallback_candidates(
            infinite_by_sample_contig,
            valid_regions_by_sample_contig,
        )
    )
    for label, intervals in infinite_by_label.items():
        effective[label] = (label_contig[label], tuple(intervals))
    print(
        "[graphreftovcf] interval ownership: "
        f"{len(effective)} aligned label(s); "
        f"primary={primary_bases} bp; uncovered negative fallback="
        f"{fallback_bases} bp; selected -inf fallback={infinite_bases} bp "
        f"(resolved {competing_infinite_bases} competing bp, rejected "
        f"{infinite_below_score} row(s) below score "
        f"{MIN_INFINITE_FALLBACK_ALIGNMENT_SCORE}; valid-region priority="
        f"{VALID_OWNERSHIP_PRIORITY_SCORE})",
        file=sys.stderr,
    )
    return effective


def normalize_encoded_cigar_text(text: str) -> str:
    s = (text or "").strip()
    if not s:
        return s
    # Old @3A hex-encoding: decode all @XX tokens back to their characters.
    if any(tok in s.lower() for tok in ("@3a", "@3b", "@2c", "@09", "@0a", "@40")):
        s = re.sub(r'@([0-9A-Fa-f]{2})', lambda m: chr(int(m.group(1), 16)), s)
    # Legacy %-encoded graph CIGAR payloads.
    elif any(tok in s.lower() for tok in ("%3a", "%3d", "%3e", "%3c", "%2d", "%25")):
        s = urllib.parse.unquote(s)
    # Legacy form sometimes used bare '%' as ':' in graph CIGAR chunks.
    if "%" in s and ":" not in s and any(ch in s for ch in (">", "<", "*")):
        s = s.replace("%", ":")
    return s


def x_payload_query(payload: str, n: int) -> str:
    if not payload:
        return ""
    if len(payload) == n:
        return payload
    if len(payload) == 2 * n:
        return payload[n:]
    return payload[:n]


def ref_base_from_rec(ref_rec: SeqRecord, local_off: int) -> str:
    seq = ref_rec.fetch_local(local_off, local_off + 1, "+")
    return (seq[:1] or "").upper()


def has_real_n(seq: str) -> bool:
    """True only for observed N bases. Empty means unavailable, not N."""
    return bool(seq) and "N" in seq.upper()


_NON_ACGT_RE = re.compile(r"[^ACGT]", re.IGNORECASE)


def has_non_acgt(seq: str) -> bool:
    """True for an observed sequence containing any non-A/C/G/T base.

    Empty means that sequence is unavailable, not that it is ambiguous.
    Lowercase A/C/G/T (for example from a repeat-masked FASTA) is valid.
    """
    return bool(seq) and bool(_NON_ACGT_RE.search(seq))


def cigar_has_non_acgt_in_iops(cigar: str) -> bool:
    """True if an insertion payload carries a non-A/C/G/T base."""
    for match in _CIGAR_RE.finditer(cigar or ""):
        if match.group(2) != "I":
            continue
        length = int(match.group(1))
        payload = (match.group(3) or "")[:length]
        if payload and has_non_acgt(payload):
            return True
    return False


def valid_called_base(base: str) -> bool:
    return len(base or "") == 1 and base.upper() in {"A", "C", "G", "T"}


def query_edge_distance(q0: int, q1: int, query_length: int) -> int:
    """Minimum distance from a query-local event interval to either edge."""
    q0 = int(q0)
    q1 = int(q1)
    if q1 < q0:
        q0, q1 = q1, q0
    query_length = max(0, int(query_length))
    q0 = max(0, min(q0, query_length))
    q1 = max(q0, min(q1, query_length))
    return min(q0, query_length - q1)


def parse_coord_list(text: str) -> List[Tuple[str, int, int, str]]:
    return [
        coord
        for coord in (
            parse_coord_token(token)
            for token in (text or "").split(";")
        )
        if coord is not None
    ]


def event_inside_trimmed_alignment(
    contig: str,
    start: int,
    end: int,
    alignments: Sequence[Tuple[str, int, int, str]],
    trim: int,
) -> bool:
    """Return true when an event lies wholly inside one trimmed alignment.

    Insertions/deletion query breakpoints are zero-width points and use closed
    boundary checks. Nonzero spans must be wholly contained in the half-open
    reliable interior. An alignment shorter than twice ``trim`` has no
    callable interior.
    """
    start, end = int(start), int(end)
    if end < start:
        start, end = end, start
    trim = max(0, int(trim))
    for achrom, astart, aend, _strand in alignments:
        if achrom != contig:
            continue
        safe_start = int(astart) + trim
        safe_end = int(aend) - trim
        if safe_end < safe_start:
            continue
        if start == end:
            if safe_start <= start <= safe_end:
                return True
        elif safe_start <= start and end <= safe_end:
            return True
    return False


def query_seq(rec: SeqRecord, q0: int, q1: int) -> str:
    return rec.fetch_local(q0, q1, "+")


def fetch_allele_query_sequence_with_record(vals: List[str], lookup: Dict[str, SeqRecord]) -> Tuple[str, Optional[SeqRecord], int, int, str]:
    """Slice the full insertion query sequence from the FASTA using sample field coordinates.

    vals[4] = vcf-escaped query contig name
    vals[5] = genomic range  e.g. "1000-2000+" (0-based coords from header, strand-suffixed)
    Coordinates are 0-based.

    Finds the FASTA segment (SeqRecord) whose genomic interval fully covers
    0-based half-open [g0, g1), converts directly to local FASTA coordinates,
    and returns the sequence in the requested strand orientation (revcomp when
    rec.strand != requested strand).
    """
    if len(vals) < 6 or not lookup:
        return "", None, 0, 0, "+"
    qry_contig = vcf_unescape(vals[4])
    m = re.match(r'^(\d+)(?:-(\d+))?([+-]?)$', vals[5])
    if not m:
        return "", None, 0, 0, "+"
    g0 = int(m.group(1))
    g1 = int(m.group(2)) if m.group(2) else g0
    strd_l = m.group(3) or "+"
    if g1 <= g0:
        return "", None, 0, 0, strd_l
    rec = _find_covering_record(qry_contig, g0, g1, lookup)
    if rec is None:
        return "", None, 0, 0, strd_l
    # Convert genomic 0-based half-open [g0, g1) to local FASTA coordinates.
    # No +/-1 adjustment belongs here.
    if rec.strand == "+":
        local_start = g0 - rec.start
        local_end = g1 - rec.start
    else:
        local_start = rec.end - g1
        local_end = rec.end - g0
    # If record strand and requested strand differ, revcomp via fetch_local("-").
    fetch_strand = "+" if rec.strand == strd_l else "-"
    if local_start < 0 or local_end > rec.seqlen or local_end <= local_start:
        return "", rec, local_start, local_end, strd_l
    return rec.fetch_local(local_start, local_end, fetch_strand), rec, local_start, local_end, strd_l


def fetch_allele_query_sequence(vals: List[str], lookup: Dict[str, SeqRecord]) -> str:
    seq, _rec, _local_start, _local_end, _strd_l = fetch_allele_query_sequence_with_record(vals, lookup)
    return seq


def query_range_len(qrange: str) -> int:
    m = re.match(r'^(\d+)(?:-(\d+))?([+-]?)$', qrange or "")
    if not m:
        return 0
    g0 = int(m.group(1))
    g1 = int(m.group(2)) if m.group(2) else g0
    return max(0, g1 - g0)


def _looks_like_label_h_coord(text: str) -> bool:
    component = r"(?:[+-]?\d+H)(?:_?(?:[+-]?\d+H))?"
    atom = rf"(?:{component}|\.)"
    return bool(re.fullmatch(rf"{atom}(?:&{atom})*", text or ""))


def parse_hsv_allele(allele: str) -> Optional[List[str]]:
    """Parse one FORMAT/HSV allele while preserving ':' inside EXTENDGRAPHCIGAR.

    Internal SV format is
    GT:TYPE:SIZE:CIGAR:assemblycontig:query_range:allelename:label_H:offset.
    Legacy eight- and seven-field alleles are accepted and receive a missing
    TEMPLATEOFFSET (and, for seven fields, a missing LABEL_H).
    Encoded CIGAR chunks contain reference names like ``NC_060932.1:57561709...``,
    so a plain split(":") corrupts the trailing fields.

    Current VCF output declares and transposes these components into separate
    FORMAT fields. This parser handles the reconstructed internal payload as
    well as escaped and unescaped legacy HSV/HSNP envelope values.
    """
    if ":" not in allele:
        allele = vcf_unescape(allele.replace("@7C", "|"))
    first = allele.split(":", 3)
    if len(first) != 4:
        return None
    rest = first[3]
    last_offset = rest.rsplit(":", 5)
    if (
        len(last_offset) == 6
        and _looks_like_label_h_coord(last_offset[-2])
        and re.fullmatch(r"[+-]?\d+|\.", last_offset[-1] or "")
    ):
        return first[:3] + last_offset
    last_new = rest.rsplit(":", 4)
    if len(last_new) == 5 and _looks_like_label_h_coord(last_new[-1]):
        return first[:3] + last_new + ["."]
    last_old = rest.rsplit(":", 3)
    if len(last_old) == 4:
        return first[:3] + last_old + [".", "."]
    return None


def format_hsv_allele(vals: Sequence[str]) -> str:
    """Return the internal colon-delimited representation of one event."""
    return ":".join(vals)


def is_sv_format(format_text: str) -> bool:
    return format_text in {
        "HSV", SV_FORMAT, SV_FORMAT_BASE,
        SV_FORMAT_LEGACY, SV_FORMAT_LEGACY_BASE,
    }


def is_snp_format(format_text: str) -> bool:
    return format_text in {
        "HSNP", SNP_FORMAT, SNP_FORMAT_BASE,
        SNP_FORMAT_LEGACY, SNP_FORMAT_LEGACY_BASE,
    }


def _format_schema(format_text: str):
    if format_text in {
        SV_FORMAT, SV_FORMAT_BASE,
        SV_FORMAT_LEGACY, SV_FORMAT_LEGACY_BASE,
    }:
        return SV_EVENT_FIELDS, tuple(format_text.split(":"))
    if format_text in {
        SNP_FORMAT, SNP_FORMAT_BASE,
        SNP_FORMAT_LEGACY, SNP_FORMAT_LEGACY_BASE,
    }:
        return SNP_EVENT_FIELDS, tuple(format_text.split(":"))
    return (), ()


def split_sample_alleles(
    field_text: str, format_text: str,
) -> List[str]:
    """Decode one sample column into backward-compatible event payloads."""
    text = (field_text or ".").strip()
    if text in {"", ".", "0", "0/0", "./."}:
        return []
    if format_text in {"HSV", "HSNP"}:
        return [
            vcf_unescape(value.replace("@7C", "|"))
            if ":" not in value else value
            for value in text.split("|")
            if value and value not in {".", "0", "0/0", "./."}
        ]
    canonical, present = _format_schema(format_text)
    if not present:
        return []
    components = text.split(":")
    if len(components) != len(present):
        return []
    genotype = components[0]
    if genotype in {"", ".", "0", "0/0", "./."}:
        return []
    value_columns = [component.split(",") for component in components[1:]]
    event_count = max((len(values) for values in value_columns), default=0)
    if event_count == 0:
        return []
    # Every Number=. FORMAT column describes the same ordered observation
    # vector.  Broadcasting a truncated singleton column across multiple
    # events can silently assign the wrong contig, query coordinate, or
    # template offset, so reject inconsistent cardinalities.
    if any(len(values) != event_count for values in value_columns):
        return []
    output: List[str] = []
    for event_index in range(event_count):
        by_name = {"GT": genotype}
        for name, column in zip(present[1:], value_columns):
            encoded = column[event_index] if len(column) > 1 else column[0]
            # Keep VCF missing values as explicit placeholders while rebuilding
            # the internal colon-delimited allele. An empty placeholder makes
            # the right-anchored CIGAR parser shift all following fields.
            by_name[name] = "." if encoded == "." else vcf_unescape(encoded)
        values = [by_name.get(name, ".") for name in canonical]
        output.append(format_hsv_allele(values))
    return output


def format_sample_alleles(
    alleles: Sequence[str], format_text: str,
) -> str:
    """Transpose event payloads into declared VCF FORMAT fields."""
    values = []
    for event_index, allele in enumerate(alleles, 1):
        if not allele:
            continue
        value = parse_hsv_allele(allele)
        if value is None:
            raise ValueError(
                f"cannot format malformed event {event_index} for VCF "
                f"FORMAT {format_text!r}: {allele!r}"
            )
        values.append(value)
    if not values:
        return (
            "." if format_text in {"HSV", "HSNP"}
            else ":".join(["."] * len(format_text.split(":")))
        )
    if format_text in {"HSV", "HSNP"}:
        return "|".join(
            vcf_escape(format_hsv_allele(value[:8])).replace("|", "@7C")
            for value in values
        )
    canonical, present = _format_schema(format_text)
    if not present:
        raise ValueError(f"unsupported VCF FORMAT schema: {format_text!r}")
    canonical_indexes = {name: index for index, name in enumerate(canonical)}
    columns = [values[0][canonical_indexes["GT"]]]
    for name in present[1:]:
        field_index = canonical_indexes[name]
        columns.append(",".join(
            vcf_escape(value[field_index]) for value in values
        ))
    return ":".join(columns)


def format_missing_sample(symbol: str, format_text: str) -> str:
    if format_text in {"HSV", "HSNP"}:
        return symbol
    return ":".join(
        [symbol] + ["."] * (len(format_text.split(":")) - 1)
    )


def _add_coverage_span(spans: List[Tuple[str, str, int, int]], query_name: str, ref_rec: SeqRecord, f0: int, f1: int):
    g0, g1 = local_ref_to_genomic(ref_rec, f0, f1)
    if g1 < g0:
        g0, g1 = g1, g0
    if g1 == g0:
        g1 = g0 + 1
    spans.append((query_name, ref_rec.contig, g0, g1))


def query_seq_oriented(rec: SeqRecord, q0: int, q1: int, direction: str) -> str:
    """Return query bases in the orientation of an encoded insertion segment."""
    seq = rec.fetch_local(q0, q1, "+")
    return revcomp(seq) if direction == "<" else seq


def _plain_cigar_body_from_items(items: Sequence[dict], direction: str) -> str:
    body = "".join(str(x.get("text", "")) for x in items)
    if body and not body.startswith((">", "<")):
        return direction + body
    return body


def _interval_replacement_cigar(
    direction: str,
    ref_size: int,
    query_size: int,
    query_sequence: str,
) -> str:
    """Encode one merged interval as deletion plus full query replacement."""
    parts: List[str] = []
    if ref_size > 0:
        parts.append(f"{int(ref_size)}D")
    if query_size > 0:
        parts.append(f"{int(query_size)}I{query_sequence or ''}")
    return (direction or ">") + "".join(parts)


def _event_query_interval_for_overlap(a: int, b: int) -> Tuple[int, int]:
    if b < a:
        a, b = b, a
    if b == a:
        return max(0, a - 1), a + 1
    return a, b


def nested_sv_redundant(candidate: dict, main_svs: Sequence[SVUnit], accepted_nested: Sequence[SVUnit], line_no: object = ".") -> bool:
    """Drop a nested SV only when it is truly equivalent to a parent/nearby call.

    Equivalence is the same two-component size rule used for cross-sample SV
    clustering:

        max(abs(i1-i2), abs(d1-d2)) < 0.3 * max(i1,i2,d1,d2)

    plus same query/assembly contig and overlapping query interval.  This means
    a nested 515D is *not* redundant with a parent 130253 bp insertion:
    max(abs(0-130253), abs(515-0)) is far larger than 0.3 * 130253.
    """
    ci = int(candidate.get("ins_size", 0))
    cd = int(candidate.get("del_size", 0))
    cq0, cq1 = _event_query_interval_for_overlap(int(candidate.get("qry_start", 0)), int(candidate.get("qry_end", 0)))
    for sv in list(main_svs) + list(accepted_nested):
        if sv.query != candidate.get("query") or sv.qry_contig != candidate.get("qry_contig"):
            continue
        m = max(ci, sv.ins_size, cd, sv.del_size)
        diff = max(abs(ci - sv.ins_size), abs(cd - sv.del_size)) if m > 0 else 0
        threshold = 0.3 * m
        compatible = m > 0 and diff < threshold
        sq0, sq1 = _event_query_interval_for_overlap(sv.qry_start, sv.qry_end)
        ov = interval_overlap(cq0, cq1, sq0, sq1)
        if compatible and ov > 0:
            return True
    return False


def call_variants_inside_encoded_insertion(
    encoded_text: str,
    parent_query_start: int,
    query_rec: SeqRecord,
    label_rec: SeqRecord,
    source_rec: SeqRecord,
    base_label: str,
    sample: str,
    ref_lookup: Dict[str, SeqRecord],
    max_impute: int,
    bp_uncertain: int,
    min_ref_span: int,
    gap_break: int,
    gap_penalty: int,
    line_no: object = ".",
    separate_adjacent_indels: bool = False,
    snps_only: bool = False,
    row_operation_count: Optional[int] = None,
) -> Tuple[List[dict], List[SNPUnit], List[Tuple[str, str, int, int]]]:
    """Call additional variants inside one encoded insertion segment.

    The returned variants are normal raw calls for later clustering.  The only
    special marker is ``type_prefix='INS_'`` in the genotype TYPE field.  No
    recursive encoding is attempted inside an encoded insertion. With
    ``snps_only``, gaps only advance coordinates and no nested SVs are created.
    """
    nested_sv_candidates: List[dict] = []
    nested_snps: List[SNPUnit] = []
    filtered_spans: List[Tuple[str, str, int, int]] = []
    if not encoded_text or encoded_text[0] not in "><" or ":" not in encoded_text:
        return nested_sv_candidates, nested_snps, filtered_spans

    direction = encoded_text[0]
    ref_name, body = encoded_text[1:].split(":", 1)
    nested_ref_rec = resolve_record(ref_name, ref_lookup)
    if nested_ref_rec is None:
        return nested_sv_candidates, nested_snps, filtered_spans
    ops = parse_pairwise_ops(body)
    if not ops:
        return nested_sv_candidates, nested_snps, filtered_spans

    lead_h, trail_h = segment_clip_lengths(ops)
    ref_span = sum(op.n for op in ops if op.op in {"=", "X", "D"})
    if min_ref_span > 0 and ref_span < min_ref_span:
        return nested_sv_candidates, nested_snps, filtered_spans

    stream: List[dict] = []
    qcur = 0
    rcur = 0

    def add_gap(n: int, text: str):
        if n > 0 and not snps_only:
            stream.append({"kind": "gap", "text": text, "gap_len": n})

    def add_variant(tok: dict):
        tok["idx"] = len(stream)
        stream.append(tok)

    ignored_edges = terminal_indel_indexes(ops, row_operation_count=row_operation_count)
    for index, op in enumerate(ops):
        if op.n <= 0 or op.op == 'H':
            continue
        if index in ignored_edges:
            if op.op == 'I':
                qcur += op.n
            elif op.op == 'D':
                rcur += op.n
            continue
        if op.op in {"=", "M"}:
            add_gap(op.n, f"{op.n}=")
            rcur += op.n
            qcur += op.n
            continue
        if op.op == "X":
            q_payload = x_payload_query(op.payload, op.n)
            for i in range(op.n):
                f0, f1 = segment_interval_to_forward(nested_ref_rec, direction, lead_h, trail_h, rcur + i, rcur + i + 1)
                lf0, lf1 = min(f0, f1), max(f0, f1)
                rb = ref_base_from_rec(nested_ref_rec, lf0)
                if q_payload and i < len(q_payload):
                    qb = q_payload[i:i + 1].upper()
                else:
                    qb = source_rec.fetch_local(parent_query_start + qcur + i, parent_query_start + qcur + i + 1).upper()
                # Query bases follow traversal order. Reverse target
                # coordinates already reverse their positions; complement
                # each base without reversing the order a second time.
                if direction == "<":
                    qb = revcomp(qb)
                if not valid_called_base(rb) or not valid_called_base(qb):
                    g0, g1 = local_ref_to_genomic(nested_ref_rec, lf0, lf1)
                    if g1 < g0:
                        g0, g1 = g1, g0
                    filtered_spans.append((query_rec.name, nested_ref_rec.contig, g0, g1 if g1 > g0 else g0 + 1))
                    continue
                if rb == qb:
                    continue
                g0, g1 = local_ref_to_genomic(nested_ref_rec, lf0, lf1)
                if g1 < g0:
                    g0, g1 = g1, g0
                q_local = parent_query_start + qcur + i
                qg = query_rec.genomic_point_from_local(q_local)
                label = label_for_object_interval(query_rec, label_rec, base_label, q_local, q_local + 1, True, max_impute, bp_uncertain)
                if label:
                    nested_snps.append(SNPUnit(
                        query=query_rec.name,
                        sample=sample,
                        allelename=base_label,
                        chrom=nested_ref_rec.contig,
                        pos0=min(g0, g1),
                        ref_base=rb,
                        alt_base=qb,
                        qry_contig=query_rec.contig,
                        qry_pos=qg,
                        qry_strand=query_rec.strand,
                        label=label,
                        label_h=format_label_h_point(label_rec, qg),
                        type_prefix="INS_",
                        qry_local_start=q_local,
                        qry_local_end=q_local + 1,
                    ))
            add_gap(op.n, f"{op.n}X")
            rcur += op.n
            qcur += op.n
            continue
        if op.op == "D":
            if snps_only:
                rcur += op.n
                continue
            f0, f1 = segment_interval_to_forward(nested_ref_rec, direction, lead_h, trail_h, rcur, rcur + op.n)
            lf0, lf1 = min(f0, f1), max(f0, f1)
            g0, g1 = local_ref_to_genomic(nested_ref_rec, lf0, lf1)
            if g1 < g0:
                g0, g1 = g1, g0
            ref_deleted_seq = nested_ref_rec.fetch_local(lf0, lf1, "+")
            if has_non_acgt(ref_deleted_seq):
                filtered_spans.append((query_rec.name, nested_ref_rec.contig, g0, g1 if g1 > g0 else g0 + 1))
                rcur += op.n
                continue
            q_local = parent_query_start + qcur
            qg0, qg1 = query_rec.genomic_interval_from_local(q_local, q_local)
            add_variant({
                "kind": "sv",
                "op": "D",
                "text": f"{op.n}D",
                "ins_size": 0,
                "del_size": op.n,
                "ref_contig": nested_ref_rec.contig,
                "ref_start": g0,
                "ref_end": g1,
                "qry_start_local": q_local,
                "qry_end_local": q_local,
                "qry_start": qg0,
                "qry_end": qg1,
                "raw_seq": "",
                "encoded": False,
                "score": op.n,
            })
            rcur += op.n
            continue
        if op.op == "I":
            if snps_only:
                qcur += op.n
                continue
            if op.payload and len(op.payload) >= op.n:
                seq = op.payload[:op.n]
                if direction == "<":
                    seq = revcomp(seq)
            else:
                seq = query_seq_oriented(source_rec, parent_query_start + qcur, parent_query_start + qcur + op.n, direction)
            seq = (seq or "")[:op.n]
            anchor_local = segment_point_to_forward(nested_ref_rec, direction, lead_h, trail_h, rcur)
            anchor = local_ref_point_to_genomic(nested_ref_rec, anchor_local)
            if has_non_acgt(seq):
                filtered_spans.append((query_rec.name, nested_ref_rec.contig, anchor, anchor + 1))
                qcur += op.n
                continue
            q_local0 = parent_query_start + qcur
            q_local1 = parent_query_start + qcur + op.n
            qg0, qg1 = query_rec.genomic_interval_from_local(q_local0, q_local1)
            add_variant({
                "kind": "sv",
                "op": "I",
                "text": f"{op.n}I{seq}",
                "ins_size": op.n,
                "del_size": 0,
                "ref_contig": nested_ref_rec.contig,
                "ref_start": anchor,
                "ref_end": anchor,
                "qry_start_local": q_local0,
                "qry_end_local": q_local1,
                "qry_start": qg0,
                "qry_end": qg1,
                "raw_seq": seq,
                "encoded": False,
                "score": op.n,
            })
            qcur += op.n
            continue
        if op.op == "S":
            add_gap(op.n, f"{op.n}S")
            qcur += op.n

    if snps_only:
        return nested_sv_candidates, nested_snps, filtered_spans

    groups = score_merge_by_index(
        stream,
        isgap=lambda x: x.get("kind") == "gap",
        gap_break=gap_break,
        gap_penalty=gap_penalty,
        hit_score=lambda x: int(x.get("score", 1)),
        gap_size=lambda x: int(x.get("gap_len", 1)),
    )
    if separate_adjacent_indels:
        groups = split_colocated_opposite_indels(groups)
    groups = mark_scored_interval_replacements(groups)

    for group in groups:
        if not group:
            continue
        first_idx = min(g["idx"] for g in group)
        last_idx = max(g["idx"] for g in group)
        span_items = stream[first_idx:last_idx + 1]
        qloc0 = min(int(g["qry_start_local"]) for g in group)
        qloc1 = max(int(g["qry_end_local"]) for g in group)
        ref_start = min(int(g["ref_start"]) for g in group)
        ref_end = max(int(g["ref_end"]) for g in group)
        interval_replacement = any(
            bool(g.get("_interval_replacement")) for g in group
        )
        if interval_replacement:
            ins_size = max(0, qloc1 - qloc0)
            del_size = max(0, ref_end - ref_start)
            raw_ins_seq = query_seq_oriented(
                source_rec, qloc0, qloc1, direction,
            )
            event_cigar = _interval_replacement_cigar(
                direction, del_size, ins_size, raw_ins_seq,
            )
        else:
            ins_size = sum(int(g.get("ins_size", 0)) for g in group)
            del_size = sum(int(g.get("del_size", 0)) for g in group)
            raw_ins_seq = "".join(
                str(g.get("raw_seq", ""))
                for g in group if g.get("op") == "I"
            )
            event_cigar = _plain_cigar_body_from_items(span_items, direction)
        if ins_size <= 0 and del_size <= 0:
            continue
        label = label_for_object_interval(query_rec, label_rec, base_label, qloc0, qloc1, False, max_impute, bp_uncertain)
        if not label:
            continue
        label_h = format_label_h_span(query_rec, label_rec, qloc0, qloc1)
        qg0, qg1 = query_rec.genomic_interval_from_local(qloc0, qloc1)
        if qg1 < qg0:
            qg0, qg1 = qg1, qg0
        nested_sv_candidates.append({
            "query": query_rec.name,
            "sample": sample,
            "allelename": base_label,
            "ref_contig": group[0].get("ref_contig", nested_ref_rec.contig),
            "ref_start": ref_start,
            "ref_end": ref_end,
            "qry_contig": query_rec.contig,
            "qry_start": qg0,
            "qry_end": qg1,
            "qry_strand": query_rec.strand,
            "qry_local_start": qloc0,
            "qry_local_end": qloc1,
            "label": label,
            "label_h": label_h,
            "ins_size": ins_size,
            "del_size": del_size,
            "tokens": list(group),
            "cigar": event_cigar,
            "raw_ins_seq": raw_ins_seq,
            "encoded": False,
            "type_prefix": "INS_",
        })
    return nested_sv_candidates, nested_snps, filtered_spans


def parse_row_events(
    row: dict,
    ref_lookup: Dict[str, SeqRecord],
    query_lookup: Dict[str, SeqRecord],
    label_lookup: Dict[str, SeqRecord],
    sequence_lookup: Dict[str, SeqRecord],
    max_impute: int,
    bp_uncertain: int,
    gap_break: int,
    gap_penalty: int,
    var_in_insert_min_bp: int,
    next_id_start: int,
    strict: bool = True,
    edge_blackregion: int = 0,
    filter_non_acgt_sequences: bool = True,
    score_main_alignment: bool = False,
    separate_adjacent_indels: bool = False,
    call_insertion_snps: bool = True,
) -> RowParseResult:
    query_name = row["query"]
    source_rec = resolve_record(query_name, sequence_lookup) or resolve_record(query_name, query_lookup)
    base_label = first_label_token(row.get("labels", "")) or query_name
    synthetic_reference_deletion = (
        query_name.startswith("DEL") and base_label == "DEL"
    )
    alignment_edge_schema = row.get("schema") == "persample_v2"
    query_alignment_regions = parse_coord_list(
        row.get("query_alignment_region", "")
        or row.get("query_alignment_regions", "")
    )
    reference_alignment_regions = parse_coord_list(
        row.get("reference_alignment_regions", "")
    )
    # CIGAR-local q coordinates are local to the row query/object sequence
    # (column 1/2), whose bases come from -s.  Column 6 is the ownership label
    # interval and is handled below as label_rec; mixing its start with the -s
    # sequence length shifts large insertions outside the actual FASTA record.
    query_coord = parse_coord_token(row.get("query_coord", "")) or parse_coord_token(row.get("label_coords", ""))
    if query_coord is None:
        msg = f"line {row.get('line_no', '.')}: bad query coordinate for {query_name!r}: {row.get('query_coord', '')!r}"
        if strict:
            raise ValueError(msg)
        return RowParseResult([], [], [], [], [])
    q_contig, q_start, q_end, q_strand = query_coord
    if (source_rec is not None and source_rec.contig == q_contig
            and source_rec.strand == q_strand
            and source_rec.seqlen == source_rec.end - source_rec.start
            and source_rec.start <= q_start <= q_end <= source_rec.end):
        offset = q_start - source_rec.start if q_strand == '+' else source_rec.end - q_end
        source_rec = replace(source_rec, start=q_start, end=q_end, seqlen=q_end - q_start,
                             fasta_offset=source_rec.fasta_offset + offset)
    query_rec = SeqRecord(
        name=query_name,
        allelename=base_label,
        contig=q_contig,
        start=q_start,
        end=q_end,
        strand=q_strand,
        seqlen=(source_rec.seqlen if source_rec is not None and source_rec.seqlen > 0 else max(0, q_end - q_start)),
        reader=(source_rec.reader if source_rec is not None else None),
        fasta_name=(source_rec.fasta_name if source_rec is not None else query_name),
        source="input_row",
        fasta_offset=(source_rec.fasta_offset if source_rec is not None else 0),
    )
    if source_rec is None:
        source_rec = query_rec
    elif sequence_lookup is not query_lookup and resolve_record(query_name, sequence_lookup) is None:
        pass
    if filter_non_acgt_sequences:
        observed_query_sequence = query_rec.fetch_local(0, query_rec.seqlen, "+")
        if has_non_acgt(observed_query_sequence):
            return RowParseResult(
                [], [], [], [], [query_rec.name],
                row_filter_reason="non_acgt_query_sequence",
            )
    label_meta = resolve_record(base_label, label_lookup)
    if label_meta is not None:
        # Original label coordinates now come from column 3 of -m/--coord-map.
        # label_meta.start/end may already include extension for other metadata
        # uses, so use the preserved orig_* values here and apply extension
        # below exactly once for label assignment and H-distance output.
        l_contig = label_meta.orig_contig or label_meta.contig
        l_start = int(label_meta.orig_start if label_meta.orig_start is not None else label_meta.start)
        l_end = int(label_meta.orig_end if label_meta.orig_end is not None else label_meta.end)
        l_strand = label_meta.orig_strand or label_meta.strand
    else:
        # Compatibility fallback for older inputs without a matching -m row:
        # use the label coordinate embedded in the -i row and no extensions.
        # Per-sample v2 rows carry the individual query coordinate in column 3
        # but intentionally leave the legacy label-coordinate column empty.
        # Use that coordinate as the self-label fallback when no GenomeLift
        # coord-map was supplied; an explicit coord-map still takes priority
        # above and preserves its original/extended interval metadata.
        label_coord = (
            parse_coord_token(row.get("label_coords", ""))
            or parse_coord_token(row.get("query_coord", ""))
        )
        if label_coord is None:
            msg = f"line {row.get('line_no', '.')}: no -m coordinate for label {base_label!r} and bad input label coordinate: {row.get('label_coords', '')!r}"
            if strict:
                raise ValueError(msg)
            return RowParseResult([], [], [], [], [])
        l_contig, l_start, l_end, l_strand = label_coord
        label_meta = SeqRecord(
            name=base_label,
            allelename=base_label,
            contig=l_contig,
            start=l_start,
            end=l_end,
            strand=l_strand,
            seqlen=max(0, l_end - l_start),
        )

    # Build the genomic interval used for label assignment/H-distance output.
    # Positive extensions enlarge the absolute genomic start/end sides here.
    # Negative tokens are exclusions, not extensions; label_regions() applies
    # those exclusions on the local ownership axis below, so do not shrink the
    # genomic interval here or the excluded side would be counted twice.
    lsign, _lraw, leff = parse_extend_token(label_meta.left_ext, max_impute)
    rsign, _rraw, reff = parse_extend_token(label_meta.right_ext, max_impute)
    l_ext_start, l_ext_end = l_start, l_end
    if lsign != "-":
        l_ext_start = max(0, l_ext_start - leff)
    if rsign != "-":
        l_ext_end += reff
    if l_ext_end < l_ext_start:
        l_ext_end = l_ext_start
    label_rec = SeqRecord(
        name=base_label,
        allelename=base_label,
        contig=l_contig,
        start=l_ext_start,
        end=l_ext_end,
        strand=l_strand,
        seqlen=max(0, l_ext_end - l_ext_start),
        left_ext=label_meta.left_ext,
        right_ext=label_meta.right_ext,
        left_neighbor=label_meta.left_neighbor,
        right_neighbor=label_meta.right_neighbor,
        source="input_label+coord_map",
    )
    sample = get_event_haplotype(query_name, base_label)

    encoded_row = drop_empty_encoded_segments(row.get("encoded", ""))
    if strict:
        validate_graph_cigar_syntax(
            encoded_row,
            context=f"line {row.get('line_no', '.')} query {query_name!r}",
        )
    direction, ref_name, first_body = split_first_graph_segment(encoded_row)
    ref_rec = resolve_record(ref_name, ref_lookup)
    if ref_rec is None:
        msg = f"line {row.get('line_no', '.')}: reference/graph path {ref_name!r} not found for query {query_name!r}"
        if strict:
            raise KeyError(msg)
        return RowParseResult([], [], [], [], [query_rec.name])

    # Build one main alignment token stream: background ops and SV ops in row order.
    stream: List[dict] = []
    snps: List[SNPUnit] = []
    nested_sv_candidates: List[dict] = []
    coverage_spans: List[Tuple[str, str, int, int]] = []
    filtered_spans: List[Tuple[str, str, int, int]] = []

    qpos = 0
    rpos = 0
    main_ops_for_clip: List[CigarOp] = []
    chunks = list(iter_top_level_chunks(encoded_row))
    row_operation_count = graph_cigar_operation_count(chunks)
    ignore_terminal_indels = terminal_main_indel_indexes(
        chunks, row_operation_count=row_operation_count,
    )
    main_op_index = 0
    first_named_seen = False
    # We need the full main op list only for lead/trail H.  It is safe to collect
    # pathless chunks but not named encoded insertions.
    for d, text in chunks:
        if ":" in text:
            if not first_named_seen:
                _name, body = text.split(":", 1)
                main_ops_for_clip.extend(parse_pairwise_ops(body))
                first_named_seen = True
            else:
                continue
        elif first_named_seen:
            main_ops_for_clip.extend(parse_pairwise_ops(text))
    lead_h, trail_h = segment_clip_lengths(main_ops_for_clip)

    def add_owned_match_coverage(
        query_start_local: int,
        query_end_local: int,
        reference_start_local: int,
    ) -> None:
        """Add only reference spans whose equal-length query bases are owned."""
        owned = owned_object_local_intervals(
            query_rec, base_label, query_start_local, query_end_local,
        )
        if owned is None:
            owned = [(query_start_local, query_end_local)]
        for owned_start, owned_end in owned:
            if owned_end <= owned_start:
                continue
            delta_start = owned_start - query_start_local
            delta_end = owned_end - query_start_local
            forward_start, forward_end = segment_interval_to_forward(
                ref_rec,
                direction,
                lead_h,
                trail_h,
                reference_start_local + delta_start,
                reference_start_local + delta_end,
            )
            _add_coverage_span(
                coverage_spans,
                query_rec.name,
                ref_rec,
                forward_start,
                forward_end,
            )

    def add_gap(n: int, text: str):
        if n <= 0:
            return
        stream.append({"kind": "gap", "text": text, "gap_len": n})

    def add_variant(tok: dict):
        tok["idx"] = len(stream)
        stream.append(tok)

    def handle_ops(ops: Sequence[CigarOp], chunk_direction: str):
        nonlocal qpos, rpos
        del chunk_direction
        nonlocal main_op_index
        for op in ops:
            if op.n <= 0:
                continue
            if op.op == "H":
                continue
            current_main_op_index = main_op_index
            main_op_index += 1
            if current_main_op_index in ignore_terminal_indels:
                if op.op == 'I':
                    qpos += op.n
                elif op.op == 'D':
                    rpos += op.n
                continue
            if op.op in {"=", "M"}:
                add_owned_match_coverage(qpos, qpos + op.n, rpos)
                add_gap(op.n, f"{op.n}=")
                rpos += op.n
                qpos += op.n
            elif op.op == "X":
                add_owned_match_coverage(qpos, qpos + op.n, rpos)
                q_payload = x_payload_query(op.payload, op.n)
                for i in range(op.n):
                    fbase0, fbase1 = segment_interval_to_forward(ref_rec, direction, lead_h, trail_h, rpos + i, rpos + i + 1)
                    rb = ref_base_from_rec(ref_rec, min(fbase0, fbase1))
                    qb = (q_payload[i:i+1] if q_payload and i < len(q_payload) else source_rec.fetch_local(qpos + i, qpos + i + 1)).upper()
                    if direction == "<":
                        qb = revcomp(qb)
                    if not valid_called_base(rb) or not valid_called_base(qb):
                        g0, g1 = local_ref_to_genomic(ref_rec, min(fbase0, fbase1), max(fbase0, fbase1))
                        if g1 < g0:
                            g0, g1 = g1, g0
                        filtered_spans.append((query_rec.name, ref_rec.contig, g0, g1 if g1 > g0 else g0 + 1))
                        continue
                    if rb != qb:
                        g0, g1 = local_ref_to_genomic(ref_rec, min(fbase0, fbase1), max(fbase0, fbase1))
                        qg = query_rec.genomic_point_from_local(qpos + i)
                        label = label_for_object_interval(query_rec, label_rec, base_label, qpos + i, qpos + i + 1, True, max_impute, bp_uncertain)
                        if label:
                            snps.append(SNPUnit(
                                query=query_rec.name,
                                sample=sample,
                                allelename=base_label,
                                chrom=ref_rec.contig,
                                pos0=min(g0, g1),
                                ref_base=rb,
                                alt_base=qb,
                                qry_contig=query_rec.contig,
                                qry_pos=qg,
                                qry_strand=query_rec.strand,
                                label=label,
                                label_h=format_label_h_point(label_rec, qg),
                                qry_local_start=qpos + i,
                                qry_local_end=qpos + i + 1,
                            ))
                add_gap(op.n, f"{op.n}X")
                rpos += op.n
                qpos += op.n
            elif op.op == "D":
                f0, f1 = segment_interval_to_forward(ref_rec, direction, lead_h, trail_h, rpos, rpos + op.n)
                # A bracketed reference-only row emitted by
                # graphcigartoref_persample is an explicit absence call, not
                # evidence that the sample covers the deleted reference span.
                deletion_owned = object_interval_is_owned(
                    query_rec, base_label, qpos, qpos,
                )
                if (
                    not synthetic_reference_deletion
                    and deletion_owned is not False
                ):
                    _add_coverage_span(coverage_spans, query_rec.name, ref_rec, f0, f1)
                lf0, lf1 = min(f0, f1), max(f0, f1)
                g0, g1 = local_ref_to_genomic(ref_rec, lf0, lf1)
                if g1 < g0:
                    g0, g1 = g1, g0
                ref_deleted_seq = ref_rec.fetch_local(lf0, lf1, "+")
                if has_non_acgt(ref_deleted_seq):
                    filtered_spans.append((query_rec.name, ref_rec.contig, g0, g1 if g1 > g0 else g0 + 1))
                    rpos += op.n
                    continue
                qg0, qg1 = query_rec.genomic_interval_from_local(qpos, qpos)
                add_variant({
                    "kind": "sv",
                    "op": "D",
                    "text": f"{op.n}D",
                    "ins_size": 0,
                    "del_size": op.n,
                    "ref_start": g0,
                    "ref_end": g1,
                    "qry_start_local": qpos,
                    "qry_end_local": qpos,
                    "qry_start": qg0,
                    "qry_end": qg1,
                    "raw_seq": "",
                    "encoded": False,
                    "score": op.n,
                })
                rpos += op.n
            elif op.op == "I":
                seq = op.payload if op.payload and len(op.payload) >= op.n else source_rec.fetch_local(qpos, qpos + op.n)
                seq = (seq or "")[:op.n]
                if len(seq) != op.n:
                    message = (
                        f"line {row.get('line_no', '.')}: insertion sequence "
                        f"fetch returned {len(seq)} of {op.n} bases for "
                        f"{query_rec.name!r} at row-local query interval "
                        f"{qpos}-{qpos + op.n}"
                    )
                    if strict:
                        raise ValueError(message)
                    qpos += op.n
                    continue
                anchor_local = segment_point_to_forward(ref_rec, direction, lead_h, trail_h, rpos)
                anchor = local_ref_point_to_genomic(ref_rec, anchor_local)
                if has_non_acgt(seq):
                    filtered_spans.append((query_rec.name, ref_rec.contig, anchor, anchor + 1))
                    qpos += op.n
                    continue
                qg0, qg1 = query_rec.genomic_interval_from_local(qpos, qpos + op.n)
                add_variant({
                    "kind": "sv",
                    "op": "I",
                    "text": f"{op.n}I{seq}",
                    "ins_size": op.n,
                    "del_size": 0,
                    "ref_start": anchor,
                    "ref_end": anchor,
                    "qry_start_local": qpos,
                    "qry_end_local": qpos + op.n,
                    "qry_start": qg0,
                    "qry_end": qg1,
                    "raw_seq": seq,
                    "encoded": False,
                    "score": op.n,
                })
                qpos += op.n
            elif op.op == "S":
                # Soft clips are row-local sequence context, not SV calls.
                add_gap(op.n, f"{op.n}S")
                qpos += op.n

    first_named_seen = False
    for d, text in chunks:
        if ":" in text:
            name, body = text.split(":", 1)
            if not first_named_seen:
                # The first named chunk is the main segment.  Use only its body.
                first_named_seen = True
                handle_ops(parse_pairwise_ops(body), d)
            else:
                # Later named chunks are encoded insertions at the current main cursor.
                full = d + text
                _encoded_name, encoded_body = text.split(":", 1)
                qlen = sum(
                    operation.n
                    for operation in parse_pairwise_ops(encoded_body)
                    if operation.op in QUERY_OPS
                )
                if qlen <= 0:
                    continue
                raw = source_rec.fetch_local(qpos, qpos + qlen)
                if len(raw) != qlen:
                    message = (
                        f"line {row.get('line_no', '.')}: encoded insertion "
                        f"sequence fetch returned {len(raw)} of {qlen} bases "
                        f"for {query_rec.name!r} at row-local query interval "
                        f"{qpos}-{qpos + qlen}"
                    )
                    if strict:
                        raise ValueError(message)
                    qpos += qlen
                    continue
                anchor_local = segment_point_to_forward(ref_rec, direction, lead_h, trail_h, rpos)
                anchor = local_ref_point_to_genomic(ref_rec, anchor_local)
                # The encoded graph path describes the alignment, but the
                # biological insertion sequence always comes from the sample
                # assembly.  Preserve it in full for INFO/SEQ.
                raw_for_info = raw
                if call_insertion_snps or var_in_insert_min_bp > 0:
                    n_svs, n_snps, n_filtered = call_variants_inside_encoded_insertion(
                        full,
                        parent_query_start=qpos,
                        query_rec=query_rec,
                        label_rec=label_rec,
                        source_rec=source_rec,
                        base_label=base_label,
                        sample=sample,
                        ref_lookup=ref_lookup,
                        max_impute=max_impute,
                        bp_uncertain=bp_uncertain,
                        min_ref_span=var_in_insert_min_bp,
                        gap_break=gap_break,
                        gap_penalty=gap_penalty,
                        line_no=row.get("line_no", "."),
                        separate_adjacent_indels=separate_adjacent_indels,
                        snps_only=var_in_insert_min_bp <= 0,
                        row_operation_count=row_operation_count,
                    )
                    nested_sv_candidates.extend(n_svs)
                    if call_insertion_snps:
                        snps.extend(n_snps)
                    filtered_spans.extend(n_filtered)
                qg0, qg1 = query_rec.genomic_interval_from_local(qpos, qpos + qlen)
                add_variant({
                    "kind": "sv",
                    "op": "ENCODED_INS",
                    "text": full,
                    "ins_size": qlen,
                    "del_size": 0,
                    "ref_start": anchor,
                    "ref_end": anchor,
                    "qry_start_local": qpos,
                    "qry_end_local": qpos + qlen,
                    "qry_start": qg0,
                    "qry_end": qg1,
                    "raw_seq": raw_for_info,
                    "encoded": True,
                    "score": qlen,
                })
                qpos += qlen
        elif first_named_seen:
            # Pathless scattered/main continuation.  It is part of the main linear
            # alignment, not an encoded insertion.
            handle_ops(parse_pairwise_ops(text), d)

    declared_query_span = max(0, int(q_end) - int(q_start))
    if qpos != declared_query_span:
        message = (
            f"line {row.get('line_no', '.')}: graph CIGAR consumes {qpos} "
            f"query bases but QUERYCOORD spans {declared_query_span} for "
            f"{query_name!r} ({q_contig}:{q_start}-{q_end}{q_strand})"
        )
        if strict:
            raise ValueError(message)
        return RowParseResult([], [], [], [], [query_rec.name])

    def isgap(item: dict) -> bool:
        return item.get("kind") == "gap"

    groups = score_merge_by_index(
        stream,
        isgap=isgap,
        gap_break=gap_break,
        gap_penalty=gap_penalty,
        hit_score=lambda x: int(x.get("score", 1)),
        gap_size=lambda x: int(x.get("gap_len", 1)),
    )
    if separate_adjacent_indels:
        groups = split_colocated_opposite_indels(groups)
    groups = mark_scored_interval_replacements(groups)

    svs: List[SVUnit] = []
    next_id = next_id_start
    for group in groups:
        if not group:
            continue
        first_idx = min(g["idx"] for g in group)
        last_idx = max(g["idx"] for g in group)
        span_items = stream[first_idx:last_idx + 1]
        ref_start = min(int(g["ref_start"]) for g in group)
        ref_end = max(int(g["ref_end"]) for g in group)
        # For pure insertions, keep a point interval.  For pure deletions, query
        # can be a point; for complex groups, span all local anchors/ranges.
        qloc0 = min(int(g["qry_start_local"]) for g in group)
        qloc1 = max(int(g["qry_end_local"]) for g in group)
        interval_replacement = any(
            bool(g.get("_interval_replacement")) for g in group
        )
        if interval_replacement:
            ins_size = max(0, qloc1 - qloc0)
            del_size = max(0, ref_end - ref_start)
        else:
            ins_size = sum(int(g.get("ins_size", 0)) for g in group)
            del_size = sum(int(g.get("del_size", 0)) for g in group)
        if ins_size <= 0 and del_size <= 0:
            continue
        qg0, qg1 = query_rec.genomic_interval_from_local(qloc0, qloc1)
        if qg1 < qg0:
            qg0, qg1 = qg1, qg0
        label = (
            base_label
            if synthetic_reference_deletion
            else label_for_object_interval(
                query_rec, label_rec, base_label, qloc0, qloc1,
                False, max_impute, bp_uncertain,
            )
        )
        if not label:
            continue
        label_h = (
            "0H_0H"
            if synthetic_reference_deletion
            else format_label_h_span(query_rec, label_rec, qloc0, qloc1)
        )
        if interval_replacement:
            raw_ins_seq = source_rec.fetch_local(qloc0, qloc1, "+")
            fmt_cigar = _interval_replacement_cigar(
                direction, del_size, ins_size, raw_ins_seq,
            )
            has_encoded = False
        else:
            body_parts: List[str] = []
            raw_seq_parts: List[str] = []
            has_encoded = False
            for item in span_items:
                if item.get("kind") == "gap":
                    # Preserve short internal background between grouped variants.
                    body_parts.append(item.get("text", ""))
                elif item.get("op") == "ENCODED_INS":
                    body_parts.append(item.get("text", ""))
                    raw_seq_parts.append(item.get("raw_seq", ""))
                    has_encoded = True
                else:
                    body_parts.append(item.get("text", ""))
                    if item.get("op") == "I":
                        raw_seq_parts.append(item.get("raw_seq", ""))
            # If the first part is not already a named encoded segment, prefix
            # the row/main direction.  This makes pure deletions like >426D.
            body = "".join(body_parts)
            if body and not body.startswith((">", "<")):
                fmt_cigar = direction + body
            else:
                fmt_cigar = body
            raw_ins_seq = "".join(raw_seq_parts)
        sv_unit = SVUnit(
            id=next_id,
            query=query_rec.name,
            sample=sample,
            allelename=base_label,
            ref_contig=ref_rec.contig,
            ref_start=ref_start,
            ref_end=ref_end,
            qry_contig=query_rec.contig,
            qry_start=qg0,
            qry_end=qg1,
            qry_strand=query_rec.strand,
            label=label,
            ins_size=ins_size,
            label_h=label_h,
            del_size=del_size,
            tokens=list(group),
            cigar=fmt_cigar,
            raw_ins_seq=raw_ins_seq,
            encoded=has_encoded,
            type_prefix="",
            line_no=row.get("line_no", 0),
            qry_local_start=qloc0,
            qry_local_end=qloc1,
        )
        if has_non_acgt(sv_unit.raw_ins_seq) or cigar_has_non_acgt_in_iops(sv_unit.cigar):
            filtered_spans.append((
                query_rec.name,
                sv_unit.ref_contig,
                sv_unit.ref_start,
                max(sv_unit.ref_end, sv_unit.ref_start + 1),
            ))
            next_id += 1
            continue
        svs.append(sv_unit)
        next_id += 1

    accepted_nested: List[SVUnit] = []
    for cand in nested_sv_candidates:
        if nested_sv_redundant(cand, svs, accepted_nested, line_no=row.get("line_no", ".")):
            continue
        nested_unit = SVUnit(
            id=next_id,
            query=cand["query"],
            sample=cand["sample"],
            allelename=cand["allelename"],
            ref_contig=cand["ref_contig"],
            ref_start=int(cand["ref_start"]),
            ref_end=int(cand["ref_end"]),
            qry_contig=cand["qry_contig"],
            qry_start=int(cand["qry_start"]),
            qry_end=int(cand["qry_end"]),
            qry_strand=cand["qry_strand"],
            label=cand["label"],
            ins_size=int(cand["ins_size"]),
            label_h=cand.get("label_h", "."),
            del_size=int(cand["del_size"]),
            tokens=list(cand.get("tokens", [])),
            cigar=cand.get("cigar", ""),
            raw_ins_seq=cand.get("raw_ins_seq", ""),
            encoded=bool(cand.get("encoded", False)),
            type_prefix=cand.get("type_prefix", ""),
            line_no=row.get("line_no", 0),
            qry_local_start=int(cand.get("qry_local_start", -1)),
            qry_local_end=int(cand.get("qry_local_end", -1)),
        )
        if has_non_acgt(nested_unit.raw_ins_seq) or cigar_has_non_acgt_in_iops(nested_unit.cigar):
            filtered_spans.append((
                query_rec.name,
                nested_unit.ref_contig,
                nested_unit.ref_start,
                max(nested_unit.ref_end, nested_unit.ref_start + 1),
            ))
            next_id += 1
            continue
        accepted_nested.append(nested_unit)
        next_id += 1
    svs.extend(accepted_nested)

    edge_filtered_events = 0
    edge_blackregion = max(0, int(edge_blackregion))
    if edge_blackregion:
        retained_svs: List[SVUnit] = []
        for sv in svs:
            if alignment_edge_schema:
                edge_hit = (
                    not synthetic_reference_deletion
                    and not event_inside_trimmed_alignment(
                        sv.ref_contig,
                        sv.ref_start,
                        sv.ref_end,
                        reference_alignment_regions,
                        edge_blackregion,
                    )
                )
            else:
                edge_hit = (
                    not synthetic_reference_deletion
                    and sv.qry_local_start >= 0
                    and query_edge_distance(
                        sv.qry_local_start,
                        sv.qry_local_end,
                        query_rec.seqlen,
                    ) < edge_blackregion
                )
            if edge_hit:
                if not alignment_edge_schema and "*" in (sv.label or ""):
                    # Defer uncertain cross-graph events of every SV type. The
                    # global validation pass either confirms/resolves the
                    # redundant group or removes it as a singleton/inconsistent
                    # graphvcfmerge-style group.
                    sv.edge_candidate = True
                    retained_svs.append(sv)
                    continue
                filtered_spans.append((
                    query_rec.name,
                    sv.ref_contig,
                    sv.ref_start,
                    max(sv.ref_end, sv.ref_start + 1),
                ))
                edge_filtered_events += 1
            else:
                retained_svs.append(sv)
        svs = retained_svs

        retained_snps: List[SNPUnit] = []
        for snp in snps:
            if alignment_edge_schema:
                edge_hit = not event_inside_trimmed_alignment(
                    snp.chrom,
                    snp.pos0,
                    snp.pos0 + 1,
                    reference_alignment_regions,
                    edge_blackregion,
                )
            else:
                edge_hit = (
                    snp.qry_local_start >= 0
                    and query_edge_distance(
                        snp.qry_local_start,
                        snp.qry_local_end,
                        query_rec.seqlen,
                    ) < edge_blackregion
                )
            if edge_hit:
                filtered_spans.append((
                    query_rec.name, snp.chrom, snp.pos0, snp.pos0 + 1,
                ))
                edge_filtered_events += 1
            else:
                retained_snps.append(snp)
        snps = retained_snps

    main_path_evidence: List[MainPathAlignmentEvidence] = []
    if score_main_alignment:
        matches, mismatches, gap_opens, gap_extensions = (
            main_path_affine_score_components(chunks)
        )
        evidence = MainPathAlignmentEvidence(
            query=query_rec.name,
            sample=sample,
            allele_name=base_label,
            qry_contig=query_rec.contig,
            main_match_bases=matches,
            complete_query_length=max(0, int(query_rec.seqlen)),
            line_no=int(row.get("line_no", 0)),
            main_mismatch_bases=mismatches,
            main_gap_opens=gap_opens,
            main_gap_extensions=gap_extensions,
        )
        main_path_evidence.append(evidence)
        for sv in svs:
            sv.main_match_bases = evidence.main_match_bases
            sv.complete_query_length = evidence.complete_query_length
            sv.main_mismatch_bases = evidence.main_mismatch_bases
            sv.main_gap_opens = evidence.main_gap_opens
            sv.main_gap_extensions = evidence.main_gap_extensions

    from alternative_intervals import filter_result
    return filter_result(RowParseResult(
        svs,
        snps,
        coverage_spans,
        filtered_spans,
        [query_rec.name],
        edge_filtered_events=edge_filtered_events,
        main_path_evidence=main_path_evidence,
    ), row.get('alternative_intervals'))


# ---------------------------------------------------------------------------
# Coverage maps and genotyping helpers
# ---------------------------------------------------------------------------

def build_sample_interval_map(spans: Iterable[Tuple[str, str, int, int]], haplotype_mode: bool, allowed_samples: Optional[set] = None) -> Dict[str, Dict[str, List[Tuple[int, int]]]]:
    tmp: Dict[str, Dict[str, List[Tuple[int, int]]]] = defaultdict(lambda: defaultdict(list))
    for query, chrom, a, b in spans:
        if b < a:
            a, b = b, a
        if b == a:
            b = a + 1
        sample = get_haplotype(query) if haplotype_mode else query
        if allowed_samples is not None and sample not in allowed_samples:
            continue
        tmp[sample][chrom].append((a, b))
    out: Dict[str, Dict[str, List[Tuple[int, int]]]] = {}
    for sample, by_chrom in tmp.items():
        out[sample] = {chrom: merge_intervals(vals) for chrom, vals in by_chrom.items()}
    return out


def sample_has_overlap(interval_map: Dict[str, Dict[str, List[Tuple[int, int]]]], sample: str, chrom: str, a: int, b: int) -> bool:
    by_chrom = interval_map.get(sample)
    if not by_chrom:
        return False
    intervals = by_chrom.get(chrom)
    if not intervals:
        return False
    if b < a:
        a, b = b, a
    if b == a:
        b = a + 1
    i = bisect_left(intervals, (b, -10**30))
    if i <= 0:
        return False
    s, e = intervals[i - 1]
    return s < b and a < e


def missing_symbol(sample: str, chrom: str, a: int, b: int, coverage_map, filtered_map) -> str:
    if sample_has_overlap(filtered_map, sample, chrom, a, b):
        return "."
    if sample_has_overlap(coverage_map, sample, chrom, a, b):
        return "0"
    return "."


# ---------------------------------------------------------------------------
# Across-sample clustering
# ---------------------------------------------------------------------------

def _interval_gap(a0: int, a1: int, b0: int, b1: int) -> int:
    """Half-open interval gap, also handling zero-length breakpoint spans."""
    if b0 > a1:
        return b0 - a1
    if a0 > b1:
        return a0 - b1
    return 0


def _join_unique_contributors(values: Iterable[str]) -> str:
    result: List[str] = []
    seen = set()
    for value in values:
        for part in (value or ".").split("&"):
            if part not in seen:
                seen.add(part)
                result.append(part)
    return "&".join(result) if result else "."


@dataclass(frozen=True)
class CrossGraphSVEvent:
    sv: SVUnit
    first_label: str
    confirm_labels: Tuple[str, ...]
    leading_star: bool


@dataclass(frozen=True)
class CrossGraphInconsistency:
    group_id: int
    reason: str
    members: Tuple[CrossGraphSVEvent, ...]
    observed_first_labels: Tuple[str, ...]
    observed_confirmation_labels: Tuple[str, ...]
    selected_event_ids: Tuple[int, ...] = tuple()
    rejected_event_ids: Tuple[int, ...] = tuple()
    selected_variant_names: Tuple[str, ...] = tuple()
    rejected_variant_names: Tuple[str, ...] = tuple()
    selected_alignment_scores: Tuple[str, ...] = tuple()
    rejected_alignment_scores: Tuple[str, ...] = tuple()


def is_main_path_sv_observation(sv: SVUnit) -> bool:
    """Return whether an SV was called directly from the first/main path."""
    return not sv.encoded and not (sv.type_prefix or "").startswith("INS_")


def split_cross_graph_label(label: str) -> Tuple[str, Tuple[str, ...], bool, bool]:
    """Parse the graphvcfmerge ``A*B`` / ``*A*B`` label convention."""
    raw = label or ""
    if not raw or raw == ".":
        return "", tuple(), False, False
    leading = raw.startswith("*")
    # A legacy uncertain label can carry the extension-ownership suffix on
    # its self component (``A+*B``), whereas neighbor metadata contains the
    # exact allele name ``A``. Canonicalize that terminal ownership marker so
    # old and new files use the same reciprocal A*B / B*A identities.
    parts = [
        value.strip()[:-1]
        if value.strip().endswith("+") else value.strip()
        for value in raw.split("*") if value.strip()
    ]
    if not parts:
        return "", tuple(), "*" in raw, leading
    return parts[0], tuple(parts[1:]), "*" in raw, leading


def _canonical_cross_graph_label(label: str) -> str:
    value = (label or "").strip()
    return value[:-1] if value.endswith("+") else value


def _event_main_path_evidence(
    event: CrossGraphSVEvent,
) -> MainPathAlignmentEvidence:
    sv = event.sv
    return MainPathAlignmentEvidence(
        query=sv.query,
        sample=sv.sample,
        allele_name=sv.allelename or event.first_label,
        qry_contig=sv.qry_contig,
        main_match_bases=max(0, int(sv.main_match_bases)),
        complete_query_length=max(0, int(sv.complete_query_length)),
        line_no=int(sv.line_no),
        main_mismatch_bases=max(0, int(sv.main_mismatch_bases)),
        main_gap_opens=max(0, int(sv.main_gap_opens)),
        main_gap_extensions=max(0, int(sv.main_gap_extensions)),
    )


def _compare_main_path_evidence(
    left: MainPathAlignmentEvidence,
    right: MainPathAlignmentEvidence,
) -> int:
    """Compare affine score, matched bases, then earlier input row."""
    if left.alignment_score != right.alignment_score:
        return 1 if left.alignment_score > right.alignment_score else -1
    if left.main_match_bases != right.main_match_bases:
        return 1 if left.main_match_bases > right.main_match_bases else -1
    if left.line_no != right.line_no:
        return 1 if left.line_no < right.line_no else -1
    return 0


def _alignment_score_detail(evidence: MainPathAlignmentEvidence) -> str:
    return (
        f"{evidence.allele_name}={evidence.alignment_score}"
        f"(M={evidence.main_match_bases},"
        f"X={evidence.main_mismatch_bases},"
        f"GO={evidence.main_gap_opens},"
        f"GE={evidence.main_gap_extensions})"
    )


def _best_main_path_evidence_by_label(
    svs: Sequence[SVUnit],
    evidence_rows: Optional[Sequence[MainPathAlignmentEvidence]],
) -> Dict[Tuple[str, str, str], MainPathAlignmentEvidence]:
    best: Dict[Tuple[str, str, str], MainPathAlignmentEvidence] = {}

    def admit(evidence: MainPathAlignmentEvidence) -> None:
        label = _canonical_cross_graph_label(evidence.allele_name)
        if not label:
            return
        key = (evidence.sample, evidence.qry_contig, label)
        previous = best.get(key)
        if (
            previous is None
            or _compare_main_path_evidence(evidence, previous) > 0
        ):
            best[key] = evidence

    for evidence in evidence_rows or ():
        admit(evidence)
    for sv in svs:
        first, _confirms, _has_star, _leading = split_cross_graph_label(
            sv.label,
        )
        admit(_event_main_path_evidence(CrossGraphSVEvent(
            sv=sv,
            first_label=first or sv.allelename,
            confirm_labels=tuple(),
            leading_star=False,
        )))
    return best


def _cross_graph_size_ratio_ok(
    left: CrossGraphSVEvent,
    right: CrossGraphSVEvent,
    ratio_min: float,
    ratio_max: float,
) -> bool:
    left_size = max(1, left.sv.max_size)
    right_size = max(1, right.sv.max_size)
    ratio = min(left_size, right_size) / float(max(left_size, right_size))
    return ratio_min <= ratio <= ratio_max


def _pure_deletion_pair_can_join(
    left: SVUnit, right: SVUnit, query_slop: int = 1,
) -> bool:
    if not (
        left.ins_size == 0 and left.del_size > 0
        and right.ins_size == 0 and right.del_size > 0
        and left.ref_contig == right.ref_contig
        and left.qry_strand == right.qry_strand
    ):
        return False
    left_direction = left.cigar[:1] if left.cigar[:1] in {">", "<"} else ""
    right_direction = right.cigar[:1] if right.cigar[:1] in {">", "<"} else ""
    if left_direction != right_direction:
        return False
    if _interval_gap(
        left.qry_start, left.qry_end, right.qry_start, right.qry_end,
    ) > max(0, int(query_slop)):
        return False
    return _interval_gap(
        left.ref_start, left.ref_end, right.ref_start, right.ref_end,
    ) == 0


def _cross_graph_confirm_set(members: Sequence[CrossGraphSVEvent]) -> set:
    confirms = set()
    for member in members:
        confirms.update(value for value in member.confirm_labels if value)
        if (
            member.leading_star
            and member.first_label
            and member.confirm_labels
        ):
            confirms.add(member.first_label)
    return confirms


def _resolve_leading_cross_graph_label(
    members: Sequence[CrossGraphSVEvent],
) -> Tuple[str, Optional[CrossGraphSVEvent]]:
    ordered = sorted(members, key=lambda event: (
        event.sv.qry_start, event.sv.qry_end,
        event.sv.line_no, event.sv.id,
    ))
    if not ordered:
        return "", None
    by_source: Dict[str, CrossGraphSVEvent] = {}
    for event in ordered:
        if event.first_label and event.first_label not in by_source:
            by_source[event.first_label] = event
    first = ordered[0]
    if not first.confirm_labels:
        return first.first_label, by_source.get(first.first_label, first)
    seen_sources = {first.first_label}
    target = first.confirm_labels[0]
    while True:
        source_event = by_source.get(target)
        if source_event is None:
            return target, None
        source = source_event.first_label
        next_target = (
            source_event.confirm_labels[0]
            if source_event.confirm_labels else ""
        )
        if not next_target or next_target in seen_sources:
            return source, source_event
        seen_sources.add(source)
        target = next_target


def _cross_graph_candidate(
    current: CrossGraphSVEvent,
    candidate: CrossGraphSVEvent,
    ratio_min: float,
    ratio_max: float,
    max_distance: int,
) -> bool:
    if current.sv.svtype != candidate.sv.svtype:
        return False
    if _pure_deletion_pair_can_join(current.sv, candidate.sv):
        return True
    if not _cross_graph_size_ratio_ok(
        current, candidate, ratio_min, ratio_max,
    ):
        return False
    current_size = max(1, current.sv.max_size)
    candidate_size = max(1, candidate.sv.max_size)
    distance_limit = min(
        max(0, int(max_distance)),
        0.1 * min(current_size, candidate_size),
    )
    return candidate.sv.qry_start - current.sv.qry_start <= distance_limit


def _combined_cross_graph_provenance(
    keep: CrossGraphSVEvent,
    members: Sequence[CrossGraphSVEvent],
) -> Tuple[str, str]:
    ordered = [keep] + [member for member in members if member is not keep]
    pairs: List[Tuple[str, str]] = []
    seen = set()
    for member in ordered:
        label = member.first_label
        if not label or label in seen:
            continue
        seen.add(label)
        pairs.append((label, member.sv.label_h or "."))
    return (
        "&".join(label for label, _label_h in pairs) or ".",
        "&".join(label_h for _label, label_h in pairs) or ".",
    )


def _assemble_validated_deletion(
    keep: CrossGraphSVEvent,
    members: Sequence[CrossGraphSVEvent],
) -> Tuple[SVUnit, bool]:
    """Assemble connected deletion pieces; preserve a complete call if present."""
    if not members or not all(
        member.sv.ins_size == 0 and member.sv.del_size > 0
        for member in members
    ):
        return keep.sv, False
    ordered = sorted(members, key=lambda member: (
        member.sv.ref_start, member.sv.ref_end, member.sv.id,
    ))
    ref_start = ordered[0].sv.ref_start
    ref_end = ordered[0].sv.ref_end
    for member in ordered[1:]:
        if member.sv.ref_start > ref_end:
            return keep.sv, False
        ref_end = max(ref_end, member.sv.ref_end)
    complete = next((
        member for member in members
        if member.sv.ref_start == ref_start and member.sv.ref_end == ref_end
    ), None)
    if complete is not None:
        return complete.sv, False

    result = keep.sv
    result.ref_start = ref_start
    result.ref_end = ref_end
    result.ins_size = 0
    result.del_size = max(0, ref_end - ref_start)
    result.qry_start = min(member.sv.qry_start for member in members)
    result.qry_end = max(member.sv.qry_end for member in members)
    result.qry_local_start = -1
    result.qry_local_end = -1
    direction = result.cigar[:1] if result.cigar[:1] in {">", "<"} else ">"
    result.cigar = f"{direction}{result.del_size}D"
    result.tokens = [
        token for member in members for token in member.sv.tokens
    ]
    return result, True


def validate_cross_graph_svs(
    svs: Sequence[SVUnit],
    ratio_min: float = 0.7,
    ratio_max: float = 1.0,
    max_distance: int = 500,
    query_lookup: Optional[Mapping[str, SeqRecord]] = None,
    inconsistency_groups: Optional[List[CrossGraphInconsistency]] = None,
    main_path_evidence: Optional[
        Sequence[MainPathAlignmentEvidence]
    ] = None,
    resolve_conflicts_by_alignment_score: bool = False,
) -> Tuple[List[SVUnit], int, int, int]:
    """Apply graphvcfmerge-style confirmation to every star-labelled SV type.

    Consistent duplicate INS/SUB/DEL groups retain one full representative.
    A consistent connected deletion with no complete representative is
    assembled across its rows.  Singleton or inconsistent ``*`` groups are
    removed exactly as in graphvcfmerge.py unless score resolution is
    explicitly enabled.  In score mode, one conflicting PA survives only if
    it is strictly better than every opposing PA after affine score, matched
    bases, and input-order tie breaking.
    """
    output = [sv for sv in svs if "*" not in (sv.label or "")]
    if resolve_conflicts_by_alignment_score:
        # Cross-validation is only between main-path observations.  Preserve
        # secondary-path variants, but strip their star syntax so they cannot
        # be treated as reciprocal main-path support later.
        for sv in svs:
            if "*" not in (sv.label or "") or is_main_path_sv_observation(sv):
                continue
            first, _confirms, _has_star, _leading = split_cross_graph_label(
                sv.label,
            )
            sv.label = first or sv.allelename or "."
            output.append(sv)
    evidence_by_label = _best_main_path_evidence_by_label(
        svs, main_path_evidence,
    )
    grouped: Dict[Tuple[str, str], List[CrossGraphSVEvent]] = defaultdict(list)
    for sv in svs:
        first, confirms, has_star, leading = split_cross_graph_label(sv.label)
        if not has_star or (
            resolve_conflicts_by_alignment_score
            and not is_main_path_sv_observation(sv)
        ):
            continue
        grouped[(sv.sample, sv.qry_contig)].append(CrossGraphSVEvent(
            sv=sv,
            first_label=first,
            confirm_labels=confirms,
            leading_star=leading,
        ))

    validated_groups = 0
    assembled_deletions = 0
    dropped_events = 0
    next_inconsistency_group = 1
    for events in grouped.values():
        values = sorted(events, key=lambda event: (
            event.sv.qry_start, event.sv.qry_end,
            event.sv.line_no, event.sv.id,
        ))
        occupied = set()
        for index, current in enumerate(values):
            if current.sv.id in occupied:
                continue
            members = [current]
            current_size = max(1, current.sv.max_size)
            for candidate in values[index + 1:]:
                if candidate.sv.id in occupied:
                    continue
                distance = candidate.sv.qry_start - current.sv.qry_start
                if distance > max(1, int(0.1 * current_size)):
                    break
                if candidate.first_label == current.first_label:
                    break
                if _cross_graph_candidate(
                    current, candidate,
                    ratio_min, ratio_max, max_distance,
                ):
                    members.append(candidate)

            first_labels = {
                member.first_label for member in members
                if member.first_label
            }
            confirms = _cross_graph_confirm_set(members)
            consistent = (
                len(members) >= 2
                and len(first_labels) >= 2
                and first_labels == confirms
            )
            occupied.update(member.sv.id for member in members)
            if not consistent:
                if len(members) < 2:
                    base_reason = "singleton"
                elif len(first_labels) < 2:
                    base_reason = "non_distinct_source_labels"
                else:
                    base_reason = "confirmation_label_mismatch"

                best_member = members[0]
                best_evidence = _event_main_path_evidence(best_member)
                for member in members[1:]:
                    evidence = _event_main_path_evidence(member)
                    if _compare_main_path_evidence(
                        evidence, best_evidence,
                    ) > 0:
                        best_member = member
                        best_evidence = evidence

                opponent_evidence: List[MainPathAlignmentEvidence] = []
                seen_opponents = set()

                def add_opponent(
                    evidence: Optional[MainPathAlignmentEvidence],
                ) -> None:
                    if evidence is None:
                        return
                    identity = (
                        evidence.query,
                        evidence.allele_name,
                        evidence.line_no,
                        evidence.main_match_bases,
                        evidence.main_mismatch_bases,
                        evidence.main_gap_opens,
                        evidence.main_gap_extensions,
                    )
                    if identity in seen_opponents:
                        return
                    seen_opponents.add(identity)
                    opponent_evidence.append(evidence)

                for member in members:
                    if member is not best_member:
                        add_opponent(_event_main_path_evidence(member))
                opposing_labels = set(confirms)
                opposing_labels.update(first_labels)
                opposing_labels.discard(best_member.first_label)
                for label in sorted(opposing_labels):
                    add_opponent(evidence_by_label.get((
                        best_member.sv.sample,
                        best_member.sv.qry_contig,
                        _canonical_cross_graph_label(label),
                    )))

                score_winner = (
                    resolve_conflicts_by_alignment_score
                    and bool(opponent_evidence)
                    and all(
                        _compare_main_path_evidence(
                            best_evidence, opponent,
                        ) > 0
                        for opponent in opponent_evidence
                    )
                )
                if score_winner:
                    reason = (
                        f"{base_reason}_resolved_by_main_alignment_score"
                    )
                    best_member.sv.label = (
                        best_member.first_label
                        or best_member.sv.allelename
                        or "."
                    )
                    best_member.sv.cross_row_merged = True
                    output.append(best_member.sv)
                    validated_groups += 1
                    rejected_members = [
                        member for member in members
                        if member is not best_member
                    ]
                    dropped_events += len(rejected_members)
                    selected_evidence = [best_evidence]
                    rejected_evidence = opponent_evidence
                else:
                    reason = (
                        f"{base_reason}_rejected_by_main_alignment_score"
                        if resolve_conflicts_by_alignment_score
                        and opponent_evidence
                        else f"{base_reason}_removed_as_conflict"
                    )
                    rejected_members = list(members)
                    dropped_events += len(rejected_members)
                    selected_evidence = (
                        [
                            opponent for opponent in opponent_evidence
                            if _compare_main_path_evidence(
                                opponent, best_evidence,
                            ) >= 0
                        ]
                        if resolve_conflicts_by_alignment_score else []
                    )
                    rejected_evidence = [
                        _event_main_path_evidence(member)
                        for member in rejected_members
                    ]

                def unique_names(
                    values: Sequence[MainPathAlignmentEvidence],
                ) -> Tuple[str, ...]:
                    return tuple(dict.fromkeys(
                        value.allele_name or value.query for value in values
                    ))

                if inconsistency_groups is not None:
                    inconsistency_groups.append(CrossGraphInconsistency(
                        group_id=next_inconsistency_group,
                        reason=reason,
                        members=tuple(members),
                        observed_first_labels=tuple(sorted(first_labels)),
                        observed_confirmation_labels=tuple(sorted(confirms)),
                        selected_event_ids=(
                            (best_member.sv.id,) if score_winner else tuple()
                        ),
                        rejected_event_ids=tuple(
                            member.sv.id for member in rejected_members
                        ),
                        selected_variant_names=unique_names(
                            selected_evidence,
                        ),
                        rejected_variant_names=unique_names(
                            rejected_evidence,
                        ),
                        selected_alignment_scores=tuple(
                            _alignment_score_detail(value)
                            for value in selected_evidence
                        ),
                        rejected_alignment_scores=tuple(
                            _alignment_score_detail(value)
                            for value in rejected_evidence
                        ),
                    ))
                next_inconsistency_group += 1
                continue

            insertion_members = [
                member for member in members if member.sv.ins_size > 0
            ]
            nonleading: List[CrossGraphSVEvent] = []
            if insertion_members:
                insertion_targets = insertion_sequence_targets(
                    [member.sv for member in insertion_members],
                    query_lookup or {},
                )
                # A validated duplicate insertion group must retain its
                # largest observation, irrespective of input row order or
                # whether that observation carried the leading-star spelling.
                keep_event = max(insertion_members, key=lambda member: (
                    insertion_targets[member.sv.id].unmasked_bases,
                    len(insertion_targets[member.sv.id].sequence),
                    member.sv.ins_size,
                    -member.sv.line_no,
                    -member.sv.id,
                ))
            else:
                nonleading = [
                    member for member in members if not member.leading_star
                ]
            if not insertion_members and nonleading:
                keep_event = min(nonleading, key=lambda member: (
                    member.sv.line_no, member.sv.id,
                ))
            elif not insertion_members:
                resolved_label, resolved_source = (
                    _resolve_leading_cross_graph_label(members)
                )
                keep_event = resolved_source or min(members, key=lambda member: (
                    member.sv.line_no, member.sv.id,
                ))
                if resolved_label and resolved_label != keep_event.first_label:
                    keep_event = min(members, key=lambda member: (
                        0 if member.first_label == resolved_label else 1,
                        member.sv.line_no, member.sv.id,
                    ))

            kept_sv, assembled = _assemble_validated_deletion(
                keep_event, members,
            )
            actual_keep = next(
                (member for member in members if member.sv is kept_sv),
                keep_event,
            )
            combined_label, combined_label_h = (
                _combined_cross_graph_provenance(actual_keep, members)
            )
            kept_sv.label = combined_label
            kept_sv.label_h = combined_label_h
            kept_sv.allelename = _join_unique_contributors(
                member.first_label for member in members
            )
            kept_sv.edge_candidate = any(
                member.sv.edge_candidate for member in members
            )
            kept_sv.cross_row_merged = True
            output.append(kept_sv)
            validated_groups += 1
            assembled_deletions += int(assembled)

    output.sort(key=lambda sv: (
        sv.ref_contig, sv.ref_start, sv.ref_end, sv.id,
    ))
    return output, validated_groups, assembled_deletions, dropped_events


def default_inconsistency_report_path(output_path: str) -> str:
    if output_path.lower().endswith(".vcf"):
        return output_path[:-4] + ".inconsistencies.tsv"
    return output_path + ".inconsistencies.tsv"


def _report_field(value) -> str:
    return str(value).replace("\t", " ").replace("\r", " ").replace("\n", " ")


def _report_coordinate(sv: SVUnit, query: bool = False) -> str:
    if query:
        return f"{sv.qry_contig}:{sv.qry_start}-{sv.qry_end}{sv.qry_strand}"
    return f"{sv.ref_contig}:{sv.ref_start}-{sv.ref_end}"


def write_cross_graph_inconsistency_report(
    output_path: str,
    groups: Sequence[CrossGraphInconsistency],
) -> None:
    """Atomically write one diagnostic row per conflicting event."""
    header = (
        "collision_group", "reason", "decision", "sample", "event_id",
        "source_line", "query_record", "variant_name", "allele_name",
        "main_path_matched_bases", "complete_row_query_length",
        "main_path_mismatch_bases", "main_path_gap_opens",
        "main_path_gap_extensions", "main_path_alignment_score",
        "selected_variant_names", "selected_alignment_scores",
        "rejected_variant_names", "rejected_alignment_scores",
        "variant_type", "variant_size", "insertion_size", "deletion_size",
        "reference_contig", "reference_start_0based",
        "reference_end_0based_exclusive", "query_contig",
        "query_start_0based", "query_end_0based_exclusive", "query_strand",
        "label", "first_label", "confirmation_labels",
        "observed_first_labels", "observed_confirmation_labels",
        "collision_member_count", "peer_event_ids",
        "peer_reference_coordinates", "peer_query_coordinates",
    )
    temporary = output_path + f".tmp.{os.getpid()}"
    try:
        with open(temporary, "wt", encoding="utf-8") as handle:
            handle.write("\t".join(header) + "\n")
            for group in groups:
                for member in group.members:
                    sv = member.sv
                    decision = (
                        "selected"
                        if sv.id in group.selected_event_ids else "rejected"
                    )
                    peers = [
                        peer for peer in group.members if peer is not member
                    ]
                    row = (
                        group.group_id, group.reason, decision, sv.sample,
                        sv.id, sv.line_no, sv.query, sv.allelename,
                        sv.allelename, sv.main_match_bases,
                        sv.complete_query_length, sv.main_mismatch_bases,
                        sv.main_gap_opens, sv.main_gap_extensions,
                        sv.main_alignment_score,
                        "&".join(group.selected_variant_names) or ".",
                        ";".join(group.selected_alignment_scores) or ".",
                        "&".join(group.rejected_variant_names) or ".",
                        ";".join(group.rejected_alignment_scores) or ".",
                        sv.svtype, sv.max_size, sv.ins_size, sv.del_size,
                        sv.ref_contig, sv.ref_start, sv.ref_end,
                        sv.qry_contig, sv.qry_start, sv.qry_end,
                        sv.qry_strand, sv.label or ".",
                        member.first_label or ".",
                        "&".join(member.confirm_labels) or ".",
                        "&".join(group.observed_first_labels) or ".",
                        "&".join(group.observed_confirmation_labels) or ".",
                        len(group.members),
                        ",".join(str(peer.sv.id) for peer in peers) or ".",
                        ",".join(
                            _report_coordinate(peer.sv) for peer in peers
                        ) or ".",
                        ",".join(
                            _report_coordinate(peer.sv, query=True)
                            for peer in peers
                        ) or ".",
                    )
                    handle.write("\t".join(map(_report_field, row)) + "\n")
        os.replace(temporary, output_path)
    finally:
        if os.path.exists(temporary):
            os.remove(temporary)


def deduplicate_cross_graph_snps(
    snps: Sequence[SNPUnit], assembly_window: int = 10,
) -> Tuple[List[SNPUnit], int, int]:
    """Deduplicate SNP observations from overlapping graph paths.

    A position can legitimately have more than one alternate observation when
    independent local graph paths cover the same reference base (for example,
    ``A>C`` on one path and ``A>G`` on another).  The old implementation
    treated that as a fatal disagreement.  That is too strict for the
    per-block graph-CIGAR output: the formatter already groups observations by
    ``(chrom, pos, REF, ALT)`` and can emit the distinct alternate alleles as
    separate HSNP rows.  We therefore deduplicate only observations with the
    same REF/ALT and retain distinct alleles instead of aborting conversion.

    ``assembly_window`` still collapses repeated observations of the same
    allele on the same assembly contig, preserving all contributing labels.

    Returns ``(deduplicated, merged_count, multiallelic_site_count)``.  The
    validation and grouping run as sorted-run sweeps with O(1) extra memory
    instead of per-site dictionaries of sets/lists (multiple transient GB at
    cohort scale); emitted rows and counts are unchanged.
    """
    snps = list(snps)
    snps.sort(key=lambda snp: (
        snp.chrom, snp.pos0, snp.ref_base, snp.alt_base,
    ))
    total = len(snps)
    multiallelic_sites = 0
    index = 0
    while index < total:
        chrom = snps[index].chrom
        pos0 = snps[index].pos0
        refs = set()
        alleles = 0
        previous_allele = None
        site_end = index
        while (
            site_end < total
            and snps[site_end].chrom == chrom
            and snps[site_end].pos0 == pos0
        ):
            refs.add(snps[site_end].ref_base)
            allele = (snps[site_end].ref_base, snps[site_end].alt_base)
            if allele != previous_allele:
                alleles += 1
                previous_allele = allele
            site_end += 1
        if len(refs) > 1:
            raise ValueError(
                "cross-graph SNP reference disagreement at "
                f"{chrom}:{pos0 + 1}: {', '.join(sorted(refs))}"
            )
        if alleles > 1:
            multiallelic_sites += 1
        index = site_end

    snps.sort(key=lambda snp: (
        snp.sample, snp.chrom, snp.pos0, snp.ref_base, snp.alt_base,
        snp.qry_contig, snp.qry_pos, snp.query, snp.label,
    ))
    output: List[SNPUnit] = []
    merged = 0
    window = max(0, int(assembly_window))
    index = 0
    while index < total:
        run_key = (
            snps[index].sample, snps[index].chrom, snps[index].pos0,
            snps[index].ref_base, snps[index].alt_base,
        )
        run_end = index
        while run_end < total and (
            snps[run_end].sample, snps[run_end].chrom, snps[run_end].pos0,
            snps[run_end].ref_base, snps[run_end].alt_base,
        ) == run_key:
            run_end += 1
        members = snps[index:run_end]
        index = run_end
        kept_groups: List[List[SNPUnit]] = []
        for snp in members:
            duplicate_group = next((
                group for group in kept_groups
                if group[0].qry_contig == snp.qry_contig
                and abs(group[0].qry_pos - snp.qry_pos) < window
            ), None)
            if duplicate_group is None:
                kept_groups.append([snp])
            else:
                duplicate_group.append(snp)
        for group in kept_groups:
            keep = group[0]
            if len(group) > 1:
                keep.label = _join_unique_contributors(
                    snp.label for snp in group
                )
                keep.label_h = _join_unique_contributors(
                    snp.label_h for snp in group
                )
                keep.allelename = _join_unique_contributors(
                    snp.allelename for snp in group
                )
                merged += len(group) - 1
            output.append(keep)
    output.sort(key=lambda snp: (
        snp.chrom, snp.pos0, snp.alt_base,
        snp.sample, snp.qry_contig, snp.qry_pos,
    ))
    return output, merged, multiallelic_sites


def merge_cross_row_deletions(
    svs: Sequence[SVUnit],
    query_breakpoint_slop: int = 1,
) -> Tuple[List[SVUnit], int, int]:
    """Join one deletion split across adjacent graphcigartoref rows.

    Only pure, top-level deletions from different input lines are candidates.
    They must use the same sample, assembly contig/strand and CIGAR direction;
    their reference intervals must overlap or touch and their assembly
    breakpoints must be the same (allowing one base of coordinate rounding).
    This deliberately narrow rule avoids combining merely nearby deletions.

    The returned FORMAT labels retain every contributing local block in order,
    joined by ``&`` in both the block-name and block-coordinate fields.
    """
    candidates: Dict[Tuple[str, str, str, str, str], List[SVUnit]] = defaultdict(list)
    passthrough: List[SVUnit] = []
    for sv in svs:
        direction = sv.cigar[:1] if sv.cigar[:1] in {">", "<"} else ""
        if (
            sv.ins_size == 0
            and sv.del_size > 0
            and not sv.encoded
            and not sv.type_prefix
        ):
            candidates[(
                sv.sample, sv.ref_contig, sv.qry_contig,
                sv.qry_strand, direction,
            )].append(sv)
        else:
            passthrough.append(sv)

    merged_groups = 0
    absorbed_units = 0
    output = list(passthrough)
    slop = max(0, int(query_breakpoint_slop))
    for group in candidates.values():
        group = sorted(group, key=lambda value: (
            value.ref_start, value.ref_end, value.qry_start,
            value.qry_end, value.line_no, value.id,
        ))
        parent = list(range(len(group)))

        def root(index: int) -> int:
            while parent[index] != index:
                parent[index] = parent[parent[index]]
                index = parent[index]
            return index

        def union(left: int, right: int) -> None:
            left_root = root(left)
            right_root = root(right)
            if left_root != right_root:
                parent[max(left_root, right_root)] = min(left_root, right_root)

        active: List[int] = []
        for index, current in enumerate(group):
            active = [
                old_index for old_index in active
                if group[old_index].ref_end >= current.ref_start
            ]
            for old_index in active:
                previous = group[old_index]
                if previous.line_no == current.line_no:
                    continue
                if _interval_gap(
                    previous.qry_start, previous.qry_end,
                    current.qry_start, current.qry_end,
                ) > slop:
                    continue
                union(old_index, index)
            active.append(index)

        components: Dict[int, List[SVUnit]] = defaultdict(list)
        for index, sv in enumerate(group):
            components[root(index)].append(sv)
        for members in components.values():
            members.sort(key=lambda value: (
                value.ref_start, value.ref_end, value.line_no, value.id,
            ))
            if len(members) == 1:
                output.append(members[0])
                continue
            ref_start = min(value.ref_start for value in members)
            ref_end = max(value.ref_end for value in members)
            first = members[0]
            direction = first.cigar[:1] if first.cigar[:1] in {">", "<"} else ">"
            first.ref_start = ref_start
            first.ref_end = ref_end
            first.del_size = max(0, ref_end - ref_start)
            first.qry_start = min(value.qry_start for value in members)
            first.qry_end = max(value.qry_end for value in members)
            first.qry_local_start = -1
            first.qry_local_end = -1
            first.query = _join_unique_contributors(
                value.query for value in members
            )
            first.allelename = _join_unique_contributors(
                value.allelename for value in members
            )
            first.label = _join_unique_contributors(
                value.label for value in members
            )
            first.label_h = _join_unique_contributors(
                value.label_h for value in members
            )
            first.tokens = [
                token for value in members for token in value.tokens
            ]
            first.cigar = f"{direction}{first.del_size}D"
            first.line_no = min(value.line_no for value in members)
            first.edge_candidate = any(value.edge_candidate for value in members)
            first.cross_row_merged = True
            output.append(first)
            merged_groups += 1
            absorbed_units += len(members) - 1

    output.sort(key=lambda value: (
        value.ref_contig, value.ref_start, value.ref_end, value.id,
    ))
    return output, merged_groups, absorbed_units

def ref_interval_distance(a0: int, a1: int, b0: int, b1: int) -> int:
    if b0 > a1:
        return b0 - a1
    if a0 > b1:
        return a0 - b1
    return 0


def sv_size_compatible(a: Tuple[int, int], b: Tuple[int, int]) -> bool:
    i1, d1 = a
    i2, d2 = b
    m = max(i1, i2, d1, d2)
    if m <= 0:
        return False
    return max(abs(i1 - i2), abs(d1 - d2)) < 0.3 * m


def sv_can_merge_tuple(x: Tuple, y: Tuple) -> bool:
    # tuple: (local_idx, global_id, start, end, ins, del, sample)
    _li1, _gid1, s1, e1, i1, d1, _sample1 = x
    _li2, _gid2, s2, e2, i2, d2, _sample2 = y
    if not sv_size_compatible((i1, d1), (i2, d2)):
        return False
    r1 = max(i1, d1, 1)
    r2 = max(i2, d2, 1)
    return ref_interval_distance(s1, e1, s2, e2) <= min(r1, r2)


def cluster_sv_block(
    block: List[Tuple[int, int, int, int, int, str]],
) -> List[List[int]]:
    """Cluster a single chromosome block.

    Input tuples are (global_id, start, end, ins_size, del_size, sample) and
    already all come from the same chromosome/block.  A cluster may contain at
    most one event from each haplotype/sample; same-path calls separated by the
    within-path hard gap must remain distinct VCF rows.  Return clusters as
    global ids.
    """
    if not block:
        return []
    items = [
        (i, gid, s, e, ins, dele, sample)
        for i, (gid, s, e, ins, dele, sample) in enumerate(block)
    ]
    items.sort(key=lambda x: (x[2], x[3], x[1]))
    parent = list(range(len(items)))
    component_samples = [{str(item[6])} for item in items]

    def root(i: int) -> int:
        while parent[i] != i:
            parent[i] = parent[parent[i]]
            i = parent[i]
        return i

    def union(a: int, b: int):
        ra = root(a)
        rb = root(b)
        if ra == rb or component_samples[ra].intersection(component_samples[rb]):
            return False
        if ra < rb:
            parent[rb] = ra
            component_samples[ra].update(component_samples[rb])
        else:
            parent[ra] = rb
            component_samples[rb].update(component_samples[ra])
        return True

    active: List[int] = []
    for idx, item in enumerate(items):
        _li, _gid, s, e, ins, dele, _sample = item
        new_active: List[int] = []
        for j in active:
            prev = items[j]
            _pli, _pgid, ps, pe, pins, pdel, _psample = prev
            # If gap exceeds previous radius, previous can never merge with this
            # or later items because items are sorted by start.
            if s - pe <= max(pins, pdel, 1):
                new_active.append(j)
        active = new_active
        for j in active:
            if sv_can_merge_tuple(item, items[j]):
                union(idx, j)
        active.append(idx)

    grouped: Dict[int, List[int]] = defaultdict(list)
    for i, item in enumerate(items):
        grouped[root(i)].append(item[1])
    return [sorted(vals) for vals in grouped.values()]


def make_sv_blocks(
    svs: Sequence[SVUnit],
    block_gap: int = DEFAULT_BLOCK_GAP,
) -> List[List[Tuple[int, int, int, int, int, str]]]:
    vals = sorted([
        (
            sv.ref_contig, sv.ref_start, sv.ref_end, sv.id,
            sv.ins_size, sv.del_size, sv.sample,
        )
        for sv in svs
    ])
    blocks: List[List[Tuple[int, int, int, int, int, str]]] = []
    cur: List[Tuple[int, int, int, int, int, str]] = []
    cur_chrom = None
    cur_end = None
    for chrom, start, end, gid, ins, dele, sample in vals:
        if cur and (chrom != cur_chrom or (cur_end is not None and start - cur_end > block_gap)):
            blocks.append(cur)
            cur = []
        cur.append((gid, start, end, ins, dele, sample))
        cur_chrom = chrom
        cur_end = end if cur_end is None else max(cur_end, end)
    if cur:
        blocks.append(cur)
    return blocks


def cluster_svs_parallel(svs: Sequence[SVUnit], processes: int, block_gap: int) -> List[List[int]]:
    blocks = make_sv_blocks(svs, block_gap=block_gap)
    if processes <= 1 or len(blocks) <= 1:
        out: List[List[int]] = []
        for b in blocks:
            out.extend(cluster_sv_block(b))
        return out
    _freeze_shared_heap()
    with mp_context().Pool(processes=processes) as pool:
        parts = pool.map(cluster_sv_block, blocks)
    out: List[List[int]] = []
    for part in parts:
        out.extend(part)
    return out


# ---------------------------------------------------------------------------
# Formatting
# ---------------------------------------------------------------------------

def choose_representative_svs(members: Sequence[SVUnit]) -> SVUnit:
    return max(members, key=lambda x: (x.max_size, x.ins_size + x.del_size, -x.ref_start, x.query))


def count_unmasked_bases(sequence: str) -> int:
    """Count uppercase A/C/G/T bases; lowercase sequence is masked."""
    return sum(base in "ACGT" for base in (sequence or ""))


@dataclass
class InsertionSequenceTarget:
    """One insertion sequence plus its address in the source query record."""

    sv: SVUnit
    sequence: str
    path_name: str
    record_length: int
    local_start: int
    local_end: int
    direction: str = ">"

    @property
    def unmasked_bases(self) -> int:
        return count_unmasked_bases(self.sequence)


@dataclass
class RecursiveInsertionPiece:
    """A residual I segment that may be replaced by a nested alignment."""

    sequence: str
    source: InsertionSequenceTarget
    anchor: int
    ordinal: int
    replacement: Optional[List[object]] = None

    @property
    def unmasked_bases(self) -> int:
        return count_unmasked_bases(self.sequence)


def insertion_sequence_target(
    sv: SVUnit,
    query_lookup: Mapping[str, SeqRecord],
) -> InsertionSequenceTarget:
    """Fetch the complete query side of an insertion observation.

    Local row coordinates are authoritative and preserve the CIGAR/query
    orientation.  The genomic sample-field coordinates are retained as a
    fallback for older inputs whose SVUnit lacks local offsets.
    """
    rec = resolve_record(sv.query, query_lookup)
    local_start = int(sv.qry_local_start)
    local_end = int(sv.qry_local_end)
    direction = ">"
    sequence = ""
    if (
        rec is not None
        and local_start >= 0
        and local_end > local_start
        and local_end <= rec.seqlen
    ):
        sequence = rec.fetch_local(local_start, local_end, "+")
    else:
        sequence, rec, local_start, local_end, requested_strand = (
            fetch_allele_query_sequence_with_record(
                [
                    "", "", "", "", vcf_escape(sv.qry_contig),
                    format_query_range(sv),
                ],
                query_lookup,
            )
        )
        if rec is not None and rec.strand != requested_strand:
            direction = "<"
    if not sequence:
        sequence = sv.raw_ins_seq or extract_raw_inserted_sequence_from_cigar(
            sv.cigar,
        )
        rec = None
        local_start = 0
        local_end = len(sequence)
        direction = ">"
    return InsertionSequenceTarget(
        sv=sv,
        sequence=sequence,
        path_name=(rec.name if rec is not None else ""),
        record_length=(rec.seqlen if rec is not None else len(sequence)),
        local_start=max(0, int(local_start)),
        local_end=max(0, int(local_end)),
        direction=direction,
    )


def insertion_sequence_targets(
    members: Sequence[SVUnit],
    query_lookup: Mapping[str, SeqRecord],
) -> Dict[int, InsertionSequenceTarget]:
    return {
        member.id: insertion_sequence_target(member, query_lookup)
        for member in members
    }


def largest_insertion_target(
    members: Sequence[SVUnit],
    targets: Mapping[int, InsertionSequenceTarget],
) -> InsertionSequenceTarget:
    """Choose the largest grouped insertion by uppercase sequence content."""
    return max(
        (targets[member.id] for member in members),
        key=lambda target: (
            target.unmasked_bases,
            len(target.sequence),
            target.sv.ins_size,
            -target.sv.line_no,
            -target.sv.id,
        ),
    )


def choose_info_representative_svs(
    members: Sequence[SVUnit],
    svtype: str,
    query_lookup: Optional[Mapping[str, SeqRecord]] = None,
) -> SVUnit:
    """Choose the public INFO representative.

    Every insertion row uses the member with the most uppercase (unmasked)
    inserted bases.  This replaces the old encoded-first selection rule.
    """
    if svtype == "INS":
        targets = insertion_sequence_targets(members, query_lookup or {})
        return largest_insertion_target(members, targets).sv
    return choose_representative_svs(members)


def strip_cigar_payloads(text: str) -> str:
    """Return an extended graph CIGAR without inline sequence payloads.

    X/I/S consume query sequence and may carry bases immediately after the
    operation.  Named graph paths and all operation lengths are retained.
    """
    if not text:
        return ""
    return re.sub(r"(\d+)([IXS])([A-Za-z]+)", r"\1\2", str(text))


def strip_d_payloads_keep_ix(text: str) -> str:
    if not text:
        return ""
    return re.sub(r"(\d+)D[A-Za-z]+", r"\1D", str(text))


def cluster_svtype(ins: int, dele: int) -> str:
    if ins > 0 and dele > 0 and ins == dele:
        return "SUB"
    if ins > dele:
        return "INS"
    return "DEL"


def fetch_ref_base(ref_records: Dict[str, SeqRecord], chrom: str, pos0: int) -> str:
    rec = resolve_record(chrom, ref_records)
    if rec is None:
        return "N"
    return (rec.fetch_local(pos0, pos0 + 1, "+")[:1] or "N").upper()


def fetch_ref_seq(ref_records: Dict[str, SeqRecord], chrom: str, a: int, b: int) -> str:
    rec = resolve_record(chrom, ref_records)
    if rec is None:
        return ""
    if b < a:
        a, b = b, a
    return rec.fetch_local(a, b, "+")


def format_query_range(sv: SVUnit) -> str:
    if sv.qry_end == sv.qry_start:
        return f"{sv.qry_start}{sv.qry_strand}"
    return f"{sv.qry_start}-{sv.qry_end}{sv.qry_strand}"


def _ensure_cigar_chunk(text: str) -> str:
    """Return a top-level EXTENDCIGAR chunk for plain sequence/body text."""
    text = text or ""
    if not text:
        return text
    if text[0] in "><":
        return text
    if re.fullmatch(r"[A-Za-z]+", text):
        return f">{len(text)}I{text}"
    return ">" + text


def add_row_position_padding(cigar: str, sample_ref_start: int, row_ref_start: int) -> str:
    """Prefix H padding so every sample genotype is expressed from row POS.

    X = row_ref_start, Y = sample_ref_start.  If Y>X, prefix >{Y-X}H.
    If Y<X, prefix <{X-Y}H.  H is always added at the beginning, as requested.
    """
    cigar = cigar or ""
    shift = int(sample_ref_start) - int(row_ref_start)
    if shift == 0:
        return cigar
    body = _ensure_cigar_chunk(cigar)
    if shift > 0:
        return f">{shift}H" + body
    return f"<{abs(shift)}H" + body


def split_row_position_padding(cigar: str) -> Tuple[str, str]:
    """Split artificial leading row-POS H padding from the real sample CIGAR.

    Padding inserted by add_row_position_padding() is a standalone leading
    top-level chunk followed immediately by another top-level chunk, e.g.
    >20H>335IACGT or <7H>path:... .  Ordinary encoded/raw CIGARs such as
    >20H335= are left untouched.
    """
    prefix = ""
    body = cigar or ""
    while True:
        m = re.match(r"^[<>]\d+H(?=[<>])", body)
        if not m:
            break
        prefix += m.group(0)
        body = body[m.end():]
    return prefix, body


def format_sv_sample_event(
    sv: SVUnit,
    row_ref_start: Optional[int] = None,
    cigar_override: Optional[str] = None,
) -> str:
    type_label = (sv.type_prefix or "") + sv.svtype
    cigar = sv.cigar if cigar_override is None else cigar_override
    template_offset = (
        0 if row_ref_start is None
        else int(sv.ref_start) - int(row_ref_start)
    )
    return format_hsv_allele([
        "1",
        type_label,
        str(sv.signed_size),
        cigar or ".",
        vcf_escape(sv.qry_contig),
        format_query_range(sv),
        vcf_escape(sv.label),
        vcf_escape(sv.label_h),
        str(template_offset),
    ])


def format_snp_sample_event(snp: SNPUnit) -> str:
    type_label = (snp.type_prefix or "") + "SNP"
    return format_hsv_allele([
        "1",
        type_label,
        "1",
        vcf_escape(snp.alt_base),
        vcf_escape(snp.qry_contig),
        f"{snp.qry_pos}{snp.qry_strand}",
        vcf_escape(snp.label),
        vcf_escape(snp.label_h),
    ])


def validate_sv_query_invariants(sv: SVUnit, *, context: str) -> None:
    """Fail before a nonzero insertion can lose its assembly interval."""
    if int(sv.ins_size) <= 0:
        return
    query_span = int(sv.qry_end) - int(sv.qry_start)
    sequence_length = len(sv.raw_ins_seq or "")
    cigar_query_span = encoded_query_len(sv.cigar or "")
    if (
        query_span != int(sv.ins_size)
        or sequence_length != int(sv.ins_size)
        or cigar_query_span != int(sv.ins_size)
    ):
        raise ValueError(
            f"{context}: insertion invariant failed for event {sv.id}: "
            f"SIZE={int(sv.ins_size)}, QUERYCOORD span={query_span}, "
            f"sequence length={sequence_length}, CIGAR query span="
            f"{cigar_query_span}"
        )


def _vcf_meta_quote(value: str) -> str:
    return '"' + str(value).replace("\\", "\\\\").replace('"', '\\"') + '"'


def write_reference_coverage_header(handle, sample_names, reference_coverage_map) -> None:
    """Write merged reference-consuming main graph-CIGAR spans as metadata."""
    handle.write("##referenceCoverageSource=main_graph_cigar_aligned_reference_bases\n")
    handle.write("##referenceCoverageCoordinateSystem=0-based-half-open\n")
    interval_count = sum(
        len(intervals)
        for sample in sample_names
        for intervals in reference_coverage_map.get(sample, {}).values()
    )
    handle.write(f"##referenceCoverageIntervalCount={interval_count}\n")
    for sample in sample_names:
        by_chrom = reference_coverage_map.get(sample, {})
        for chrom in sorted(by_chrom):
            for start, end in by_chrom[chrom]:
                handle.write(
                    "##referenceCoverage=<"
                    f"Sample={_vcf_meta_quote(sample)},"
                    f"Chrom={_vcf_meta_quote(chrom)},"
                    f"Start={int(start)},End={int(end)}>\n"
                )


def write_vcf_header(
    handle,
    sample_names: Sequence[str],
    kind: str,
    reference_coverage_map=None,
    edge_blackregion: int = 0,
    filter_non_acgt_sequences: bool = True,
    reference_contigs: Sequence[Tuple[str, int]] = (),
    local_template_coordinates: bool = False,
    local_templates: Sequence[LocalTemplate] = (),
):
    handle.write("##fileformat=VCFv4.2\n")
    handle.write("##source=refexalign_vcf_rewrite.py\n")
    handle.write(
        "##querySequenceAlphabetFilter="
        + ("ACGT_only" if filter_non_acgt_sequences else "disabled")
        + "\n"
    )
    handle.write(f"##queryEdgeBlackregion={max(0, int(edge_blackregion))}\n")
    handle.write(
        f"##graphAlignmentEdgeBlackregion={max(0, int(edge_blackregion))}\n"
    )
    handle.write("##perSampleV2QueryAlignmentEdgeBlackregion=0\n")
    handle.write(
        "##perSampleV2ReferenceAlignmentEdgeBlackregion="
        f"{max(0, int(edge_blackregion))}\n"
    )
    handle.write("##localBlockContributorSeparator=&\n")
    local_template_contigs: List[Tuple[str, int]] = []
    if local_template_coordinates or local_templates:
        handle.write("##referenceFreeLocalTemplateFallback=true\n")
        handle.write(
            "##referenceFreeLocusCoordinates=<CHROM=AlternativeLocusFastaID,"
            "POS=OneBasedAlternativeLocusPosition,"
            "Description=\"Reference-free local graphs use the locus "
            "FASTA sequence as an independent REF contig; original assembly "
            "coordinates are provenance only\">\n"
        )
        seen_loci = set()
        for template in local_templates:
            if template.output_contig in seen_loci:
                continue
            seen_loci.add(template.output_contig)
            local_template_contigs.append((
                template.output_contig,
                template.end - template.start,
            ))
            handle.write(
                "##alternativeLocus=<"
                f"ID={template.output_contig},"
                f"Length={template.end - template.start},"
                f"SourceHaplotype={_vcf_meta_quote(template.source_haplotype)},"
                f"SourceContig={_vcf_meta_quote(template.source_contig)},"
                f"SourceStart={template.source_start},"
                f"SourceEnd={template.source_end},"
                f"SourceStrand={template.source_strand},"
                "SourceCoordinates=ZeroBasedHalfOpen>\n"
            )
    seen_contigs = set()
    # Every value emitted in the VCF CHROM column needs a standard ##contig
    # declaration. Reference-free loci retain their richer
    # ##alternativeLocus metadata above and are also declared here for HTSlib
    # compatibility.
    for contig, length in (*reference_contigs, *local_template_contigs):
        if not contig or contig in seen_contigs:
            continue
        seen_contigs.add(contig)
        handle.write(f"##contig=<ID={contig},length={int(length)}>\n")
    handle.write('##INFO=<ID=SVTYPE,Number=1,Type=String,Description="Structural variant type">\n')
    handle.write('##INFO=<ID=END,Number=1,Type=Integer,Description="End coordinate, 1-based inclusive for symbolic alleles">\n')
    handle.write('##INFO=<ID=NSUP,Number=1,Type=Integer,Description="Number of samples/haplotypes with a valid variant label">\n')
    handle.write('##INFO=<ID=SVLEN,Number=1,Type=Integer,Description="Representative signed length: insertion positive, deletion negative">\n')
    handle.write('##INFO=<ID=MAXSIZE,Number=1,Type=Integer,Description="Maximum insertion/deletion component size among cluster members">\n')
    handle.write('##INFO=<ID=CIGAR,Number=1,Type=String,Description="Representative EXTENDGRAPHCIGAR with payloads; large alleles prefer encoded representatives when available">\n')
    handle.write('##INFO=<ID=SVINSSEQ,Number=1,Type=String,Description="Complete representative inserted sequence; required for every allele with an inserted component">\n')
    handle.write('##INFO=<ID=SVDELSEQ,Number=1,Type=String,Description="Representative deleted reference sequence for alleles with deleted component <=500 bp, otherwise .">\n')
    handle.write('##INFO=<ID=QUERYSEQ,Number=1,Type=String,Description="VCF-escaped representative query contig and 0-based query coordinate or range with strand, separated by |">\n')
    handle.write('##FORMAT=<ID=GT,Number=1,Type=String,Description="Haploid genotype: 1 variant, 0 callable reference, . no call">\n')
    handle.write('##FORMAT=<ID=TYPE,Number=.,Type=String,Description="Per-observation variant type">\n')
    handle.write('##FORMAT=<ID=SIZE,Number=.,Type=Integer,Description="Per-observation signed variant size">\n')
    handle.write('##FORMAT=<ID=EXTENDGRAPHCIGAR,Number=.,Type=String,Description="Per-observation extended graph CIGAR">\n')
    handle.write('##FORMAT=<ID=BASE,Number=.,Type=String,Description="Per-observation SNP alternate base">\n')
    handle.write('##FORMAT=<ID=ASSEMBLYCONTIG,Number=.,Type=String,Description="Per-observation assembly contig">\n')
    handle.write('##FORMAT=<ID=QUERYCOORD,Number=.,Type=String,Description="Per-observation 0-based query coordinate or range with strand">\n')
    handle.write('##FORMAT=<ID=TEMPLATEOFFSET,Number=.,Type=Integer,Description="Signed displacement from the VCF row POS to the observation breakpoint">\n')
    handle.write('##FORMAT=<ID=ALLELENAME,Number=.,Type=String,Description="Per-observation allele name; cross-graph contributors are joined by &">\n')
    handle.write('##FORMAT=<ID=LABEL_H,Number=.,Type=String,Description="Per-observation PA-relative label coordinate; cross-graph contributors are joined by &">\n')
    handle.write('##FORMAT=<ID=HSV,Number=1,Type=String,Description="Legacy escaped structural-variant sample payload accepted by downstream readers">\n')
    handle.write('##FORMAT=<ID=HSNP,Number=1,Type=String,Description="Legacy escaped SNP sample payload accepted by downstream readers">\n')
    write_reference_coverage_header(
        handle, sample_names, reference_coverage_map or {},
    )
    cols = ["#CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER", "INFO", "FORMAT"] + list(sample_names)
    handle.write("\t".join(cols) + "\n")


def sv_cluster_to_row(
    cluster_id: int,
    members: List[SVUnit],
    sample_names: Sequence[str],
    coverage_map,
    filtered_map,
    ref_records,
    query_lookup: Optional[Mapping[str, SeqRecord]] = None,
    var_in_insert_min_bp: int = 0,
) -> Optional[str]:
    if not members:
        return None
    members = sorted(members, key=lambda m: (m.ref_start, m.ref_end, m.query, m.id))
    for member in members:
        validate_sv_query_invariants(
            member,
            context=(
                f"VCF cluster {cluster_id} on {member.ref_contig}:"
                f"{member.ref_start} sample {member.sample!r}"
            ),
        )
    chrom = members[0].ref_contig
    # Full reference span of the cluster, used for coverage/no-call and
    # deletion sequence extraction.  This remains min/max so no member is lost.
    ref_start = min(m.ref_start for m in members)
    ref_end = max(m.ref_end for m in members)
    # Representative row POS: use rounded mean of member reference starts.
    # This avoids always anchoring adjacent tiny indel clusters to the leftmost
    # member and then padding nearly every genotype by 1 bp.
    row_pos = int((sum(m.ref_start for m in members) / float(len(members))) + 0.5)
    # For pure insertion clusters, keep END at POS/span start.
    if ref_end < ref_start:
        ref_end = ref_start
    ins_med = int(round(median([m.ins_size for m in members])))
    del_med = int(round(median([m.del_size for m in members])))
    svtype = cluster_svtype(ins_med, del_med)
    maxsize = max(max(m.ins_size, m.del_size) for m in members)
    insertion_targets: Dict[int, InsertionSequenceTarget] = {}
    if svtype == "INS":
        insertion_targets = insertion_sequence_targets(
            members, query_lookup or {},
        )
        representative_target = largest_insertion_target(
            members, insertion_targets,
        )
        rep = representative_target.sv
        # The merged insertion event is represented by its largest unmasked
        # member, so INFO/SVLEN and INFO/SEQ/CIGAR describe the same allele.
        svlen = rep.ins_size
    else:
        representative_target = None
        rep = choose_info_representative_svs(
            members, svtype, query_lookup,
        )
        svlen = -del_med if svtype == "DEL" else ins_med - del_med
    pos0 = max(0, row_pos - 1)
    pos1 = max(ref_end, row_pos)
    ref_base = fetch_ref_base(ref_records, chrom, pos0)
    alt = f"<{svtype}>"
    event_by_sample: Dict[str, List[SVUnit]] = defaultdict(list)
    for m in members:
        if m.label:
            event_by_sample[m.sample].append(m)
    nsup = len(event_by_sample)
    if nsup == 0:
        return None
    # Deletion sequence remains a bounded legacy staging field. Insertion
    # sequence is never size-limited: the final uncompressed VCF promises a
    # complete INFO/SEQ value for every INS/SUB observation, and SVINSSEQ is
    # the lossless fallback if QUERYSEQ coordinate recovery is unavailable.
    svdelseq = "."
    del_component = max((m.del_size for m in members), default=0)
    if del_component > 0 and del_component <= INFO_RAW_DELSEQ_MAX_BP and ref_end > ref_start:
        svdelseq = fetch_ref_seq(ref_records, chrom, ref_start, ref_end)
        if len(svdelseq) > INFO_RAW_DELSEQ_MAX_BP:
            svdelseq = "."

    svinsseq = "."
    ins_component = max((m.ins_size for m in members), default=0)
    if ins_component > 0:
        raw_candidates = [m.raw_ins_seq for m in members if m.raw_ins_seq]
        insertion_target = (
            representative_target
            if representative_target is not None
            else insertion_sequence_target(rep, query_lookup or {})
        )
        representative_sequence = insertion_target.sequence
        if representative_sequence:
            svinsseq = representative_sequence
        elif rep.raw_ins_seq:
            svinsseq = rep.raw_ins_seq
        elif raw_candidates:
            svinsseq = max(raw_candidates, key=len)
        if svinsseq == ".":
            raise ValueError(
                f"cannot emit {svtype} record on {chrom}:{row_pos}: "
                "the representative insertion has no sample-assembly "
                "sequence"
            )

    # INFO/CIGAR always comes from the same largest insertion representative
    # used for INFO/SEQ and INFO/SVLEN.
    info_cigar = rep.cigar
    # Encode the representative allele's query coordinates so the reencode pass
    # can fetch the full insertion sequence without needing CIGAR reconstruction.
    # Format: "<vcf_escaped_contig>|<g0>-<g1><strand>" — pipe avoids VCF metacharacters.
    queryseq_coord = f"{vcf_escape(rep.qry_contig)}|{format_query_range(rep)}" if rep.qry_contig else ""
    row_id = f"{svtype}_{chrom}_{row_pos}_{cluster_id}"
    cigar_overrides: Dict[int, str] = {}
    if representative_target is not None:
        cigar_overrides = recursive_insertion_cigar_overrides(
            members,
            insertion_targets,
            representative_target,
            int(var_in_insert_min_bp),
            row_id,
        )
    info_parts = [
        f"SVTYPE={svtype}",
        f"END={max(row_pos, ref_end)}",
        f"NSUP={nsup}",
        f"SVLEN={svlen}",
        f"MAXSIZE={maxsize}",
        f"CIGAR={info_cigar}",
        f"SVINSSEQ={vcf_escape(svinsseq)}",
        f"SVDELSEQ={vcf_escape(svdelseq)}",
    ]
    if queryseq_coord:
        info_parts.append(f"QUERYSEQ={queryseq_coord}")
    sample_fields: List[str] = []
    for sample in sample_names:
        evs = event_by_sample.get(sample, [])
        if evs:
            # Record each observation's exact breakpoint displacement from
            # the displayed row POS in TEMPLATEOFFSET.  The alignment CIGAR
            # stays solely on its query/template axes.
            row_offset_origin = row_pos
            sample_fields.append(format_sample_alleles([
                format_sv_sample_event(
                    e,
                    row_ref_start=row_offset_origin,
                    cigar_override=cigar_overrides.get(e.id),
                )
                for e in sorted(
                    evs, key=lambda x: (x.qry_contig, x.qry_start, x.id),
                )
            ], SV_FORMAT))
        else:
            symbol = missing_symbol(
                sample, chrom, ref_start, max(ref_end, ref_start + 1),
                coverage_map, filtered_map,
            )
            sample_fields.append(format_missing_sample(symbol, SV_FORMAT))
    row = [chrom, str(row_pos), row_id, ref_base, alt, ".", "PASS", ";".join(info_parts), SV_FORMAT] + sample_fields
    return "\t".join(row)


def snp_clusters(snps: Sequence[SNPUnit]) -> Dict[Tuple[str, int, str, str], List[SNPUnit]]:
    out: Dict[Tuple[str, int, str, str], List[SNPUnit]] = defaultdict(list)
    for s in snps:
        out[(s.chrom, s.pos0, s.ref_base, s.alt_base)].append(s)
    return out


def snp_cluster_to_row(key: Tuple[str, int, str, str], members: List[SNPUnit], sample_names: Sequence[str], coverage_map, filtered_map) -> Optional[str]:
    chrom, pos0, ref, alt = key
    event_by_sample: Dict[str, List[SNPUnit]] = defaultdict(list)
    for m in members:
        if m.label:
            event_by_sample[m.sample].append(m)
    nsup = len(event_by_sample)
    if nsup == 0:
        return None
    info = f"NSUP={nsup}"
    sample_fields = []
    for sample in sample_names:
        evs = event_by_sample.get(sample, [])
        if evs:
            sample_fields.append(format_sample_alleles([
                format_snp_sample_event(e)
                for e in sorted(
                    evs, key=lambda x: (x.qry_contig, x.qry_pos, x.query),
                )
            ], SNP_FORMAT))
        else:
            symbol = missing_symbol(
                sample, chrom, pos0, pos0 + 1, coverage_map, filtered_map,
            )
            sample_fields.append(format_missing_sample(symbol, SNP_FORMAT))
    row = [chrom, str(pos0 + 1), f"SNP_{chrom}_{pos0 + 1}_{alt}", ref, alt, ".", "PASS", info, SNP_FORMAT] + sample_fields
    return "\t".join(row)




def parse_info_field(info: str) -> Dict[str, str]:
    out: Dict[str, str] = {}
    for part in info.split(";"):
        if not part:
            continue
        if "=" in part:
            k, v = part.split("=", 1)
            out[k] = v
        else:
            out[part] = ""
    return out


def format_info_field(info: Dict[str, str], original_order: Optional[List[str]] = None) -> str:
    if original_order is None:
        original_order = list(info.keys())
    seen = set()
    parts = []
    for k in original_order:
        if k in info and k not in seen:
            seen.add(k)
            parts.append(f"{k}={info[k]}" if info[k] != "" else k)
    for k, v in info.items():
        if k not in seen:
            parts.append(f"{k}={v}" if v != "" else k)
    return ";".join(parts)


# ---------------------------------------------------------------------------
# Compact sample-level insertion encoding pipeline (SSW local alignment)
# ---------------------------------------------------------------------------

REALIGN_MIN_ALIGNED_BP = 200
REALIGN_MIN_SIMILARITY = 0.90
REALIGN_MAX_SSW_SEQ_LEN = 50_000  # skip SSW for sequences larger than this to avoid O(m*n) hangs


def _add_extend_op(ops: List[Tuple[str, int, str]], op: str, n: int, payload: str = "") -> None:
    if n <= 0:
        return
    if op == "M":
        op = "="
    if op in {"=", "D", "H"}:
        payload = ""
    if ops and ops[-1][0] == op:
        old_op, old_n, old_payload = ops[-1]
        ops[-1] = (old_op, old_n + n, old_payload + payload)
    else:
        ops.append((op, n, payload))


def _format_extend_ops(ops: Sequence[Tuple[str, int, str]]) -> str:
    out: List[str] = []
    for op, n, payload in ops:
        if n <= 0:
            continue
        out.append(f"{n}{op}{payload}" if payload else f"{n}{op}")
    return "".join(out)


def cigarextend(cigars: str, refseq: str, queryseq: str) -> str:
    """Fill and polish an SSW CIGAR against actual target/query sequence.

    Current convention:
      =/M/X consume target and query; I consumes query only; D/H consume target only.
      M/= are recomputed into exact = and X runs. X/I/S payloads are query bases.
    """
    toks = re.findall(r"(\d+)([M=XHDIS])([A-Za-z]*)", cigars or "")
    refposi = 0
    qposi = 0
    out: List[str] = []
    for n_str, op, _payload in toks:
        n = int(n_str)
        if n <= 0:
            continue
        if op in ("=", "M"):
            ref_block = refseq[refposi:refposi + n].ljust(n, "N")
            qry_block = queryseq[qposi:qposi + n].ljust(n, "N")
            last_i = 0
            in_mismatch = False
            for i, (r, q) in enumerate(zip(ref_block, qry_block)):
                same = (r.upper() == q.upper())
                if (not in_mismatch) and (not same):
                    if i - last_i > 0:
                        out.append(f"{i - last_i}=")
                    last_i = i
                    in_mismatch = True
                elif in_mismatch and same:
                    if i - last_i > 0:
                        out.append(f"{i - last_i}X{qry_block[last_i:i]}")
                    last_i = i
                    in_mismatch = False
            if n - last_i > 0:
                if in_mismatch:
                    out.append(f"{n - last_i}X{qry_block[last_i:n]}")
                else:
                    out.append(f"{n - last_i}=")
            refposi += n
            qposi += n
        elif op in ("X", "S"):
            payload = queryseq[qposi:qposi + n].ljust(n, "N")
            out.append(f"{n}{op}{payload}")
            refposi += n
            qposi += n
        elif op == "I":
            payload = queryseq[qposi:qposi + n].ljust(n, "N")
            out.append(f"{n}I{payload}")
            qposi += n
        elif op in ("D", "H"):
            out.append(f"{n}{op}")
            refposi += n
    return "".join(out)


def summarize_extend_body(body: str) -> Tuple[int, float]:
    aligned = 0
    matches = 0
    for n_s, op in re.findall(r"(\d+)([=MX])", body or ""):
        n = int(n_s)
        aligned += n
        if op in {"=", "M"}:
            matches += n
    return aligned, (matches / aligned if aligned else 0.0)


def cigar_findqrange(cigars: str, qstart: int, qend: int, posi: int = 0) -> List[object]:
    """Cut a standard EXTEND/SSW CIGAR to query interval [qstart, qend).

    This is the user's helper with current convention:
      =/M/X/S consume query + reference; D/H consume reference only; I consumes query only.
    Returns [refstart, refend, cigar] where ref offsets are relative to the
    start of the original alignment CIGAR's target/reference interval.
    """
    toks = re.findall(r"\d+[a-zA-Z=]", cigars or "")
    refposi = 0
    newcigars: List[str] = []
    refstart = -1
    refend = -1
    for cigar in toks:
        lastposi = posi
        thesize = int(cigar[:-1])
        thetype = cigar[-1]
        if thetype in {"=", "M", "X", "S"}:
            posi += thesize
            refposi += thesize
        elif thetype in {"D", "H"}:
            refposi += thesize
        elif thetype == "I":
            posi += thesize

        if posi <= qstart:
            continue

        if lastposi <= qstart:
            if thetype in {"M", "=", "X", "S", "D", "H"}:
                refstart = refposi - thesize + qstart - lastposi
            else:
                refstart = refposi
            corrsize = min(qend, posi) - qstart
            if corrsize > 0:
                newcigars.append(f"{corrsize}{thetype}")

        if posi >= qend:
            if thetype in {"M", "=", "X", "S", "D", "H"}:
                refend = refposi - posi + qend
            else:
                refend = refposi
            if lastposi > qstart:
                corrsize = qend - max(qstart, lastposi)
                if corrsize > 0:
                    newcigars.append(f"{corrsize}{thetype}")
            break

        if lastposi > qstart:
            newcigars.append(cigar)

    return [refstart, refend, "".join(newcigars)]


def graphicalign(allaligns: Sequence[Sequence[object]], allow_gap: int = REALIGN_MIN_ALIGNED_BP) -> List[List[int]]:
    """Assign overlapping query intervals to the highest-scoring alignment.

    Input rows are [qstart, qend, score, alignment_index].  Output rows are
    [alignment_index, assigned_qstart, assigned_qend].  Segments <= allow_gap
    are returned to raw/unassigned sequence.
    """
    if not allaligns:
        return []
    all_coordinates = [int(x) for y in allaligns for x in y[:2]]
    scores = [float(y[2]) for y in allaligns]
    lineindex = [int(y[3]) for y in allaligns]
    sort_index = sorted(range(len(all_coordinates)), key=lambda x: (all_coordinates[x], x % 2))
    current_alignscores: List[List[float]] = []
    allgroups: List[List[float]] = [[0, 0, 0, 0]]

    for index in sort_index:
        coordinate = all_coordinates[index]
        score = scores[index // 2]
        if index % 2 == 0:
            if len(current_alignscores) == 0 or score > current_alignscores[-1][0]:
                if len(current_alignscores):
                    allgroups[-1][-1] = coordinate
                allgroups.append([score, index // 2, coordinate, coordinate])
            bisect.insort(current_alignscores, [score, -(index // 2)])
        else:
            allgroups[-1][-1] = coordinate
            if not current_alignscores:
                continue
            highest_index = -int(current_alignscores[-1][1])
            current_alignscores.remove([scores[index // 2], -(index // 2)])
            if index // 2 == highest_index and len(current_alignscores):
                allgroups.append([current_alignscores[-1][0], -current_alignscores[-1][1], coordinate, coordinate])

    out = [[lineindex[int(x[1])], int(x[2]), int(x[3])] for x in allgroups[1:]]
    merged: List[List[int]] = []
    last_index = -1
    for index, start, end in out:
        if end <= start:
            continue
        if index == last_index and merged:
            merged[-1][-1] = end
        else:
            merged.append([index, start, end])
        last_index = index
    return [x for x in merged if x[2] - x[1] > allow_gap]


def extend_query_span_from_cigar(cigar: str) -> int:
    total = 0
    for n_s, op in re.findall(r"(\d+)([=MXIDHS])", cigar or ""):
        if op in QUERY_OPS:
            total += int(n_s)
    return total


def sample_cigar_has_named_segment(cigar: str) -> bool:
    return bool(re.search(r"[<>][^:<>]+:", cigar or ""))


def info_cigar_has_encoded_target(info_cigar: str) -> bool:
    if not info_cigar:
        return False
    for _direction, text in iter_top_level_chunks(info_cigar):
        if ":" not in text:
            continue
        _name, body = text.split(":", 1)
        ops = parse_pairwise_ops(body)
        ref_span = sum(op.n for op in ops if op.op in REF_OPS)
        qry_span = sum(op.n for op in ops if op.op in QUERY_OPS)
        if ref_span > 0 and qry_span > 0:
            return True
    return False


def extract_raw_inserted_sequence_from_cigar(cigar: str) -> str:
    """Return raw inserted bases from a raw sample/INFO insertion CIGAR."""
    if not cigar or sample_cigar_has_named_segment(cigar):
        return ""
    seqs: List[str] = []
    saw_i = False
    for _direction, text in iter_top_level_chunks(cigar):
        if ":" in text:
            return ""
        for op in parse_pairwise_ops(text):
            if op.op == "I":
                saw_i = True
                seqs.append((op.payload or "")[:op.n])
    if saw_i:
        return "".join(seqs)
    if re.fullmatch(r"[A-Za-z]+", cigar or ""):
        return cigar
    return ""


def extract_insertion_payloads_from_encoded_cigar(cigar: str) -> str:
    """Return all embedded I-op payload sequences from an encoded (named-segment) CIGAR.

    Works for both raw CIGARs and encoded CIGARs with named segments.  For named
    segments the I-op payload is the unaligned insertion sequence embedded inline.
    For raw top-level chunks the I-op payload is the entire raw insertion sequence.
    """
    if not cigar:
        return ""
    seqs: List[str] = []
    for _direction, text in iter_top_level_chunks(cigar):
        body = text.split(":", 1)[1] if ":" in text else text
        for op in parse_pairwise_ops(body):
            if op.op == "I" and op.payload:
                seqs.append(op.payload[:op.n])
    return "".join(seqs)


def reconstruct_full_query_from_encoded_cigar(cigar: str, ref_lookup: Dict) -> str:
    """Reconstruct the full variant range (complete insertion sequence) from an encoded CIGAR.

    The full variant range = all query-consuming bases in order:
      - Named segment M/X/= ops  → fetch from reference via ref_lookup
      - Named segment I-ops       → embedded payload (insertion of the insertion)
      - Raw top-level chunks      → embedded I-op sequence
    """
    if not cigar:
        return ""
    parts: List[str] = []
    for direction, text in iter_top_level_chunks(cigar):
        if ":" not in text:
            for op in parse_pairwise_ops(text):
                if op.op == "I" and op.payload:
                    parts.append(op.payload[:op.n])
        else:
            name, body = text.split(":", 1)
            rec = resolve_record(name, ref_lookup) if ref_lookup else None
            h_acc = 0
            for op in parse_pairwise_ops(body):
                if op.op == "H":
                    h_acc += op.n
                elif op.op in ("=", "M", "X"):
                    if rec:
                        if direction == ">":
                            seq = rec.fetch_local(h_acc, h_acc + op.n, "+")
                        else:
                            local_end = rec.seqlen - h_acc
                            local_start = local_end - op.n
                            seq = rec.fetch_local(max(0, local_start), max(0, local_end), "+")
                            seq = revcomp(seq)
                        parts.append(seq)
                    h_acc += op.n
                elif op.op == "I":
                    if op.payload:
                        parts.append(op.payload[:op.n])
                elif op.op == "D":
                    h_acc += op.n
    return "".join(parts)


def decode_graph_cigar_to_seq(
    cigar: str,
    ref_lookup: Dict,
    mainrefname: str = "",
    mainrefpos: int = 0,
    mainstrd: str = "+",
) -> str:
    """Reconstruct the full query (insertion) sequence from a graph-encoded CIGAR.

    Main-path segments (no ':') all share one continuous main reference at
    mainrefname:mainrefpos on strand mainstrd.  Alt-path segments (with ':')
    each carry their own independent reference template.

    Strategy: pre-fetch each reference region once (trimming H-clips), then
    walk ops in forward order with a per-segment cursor.
    """
    if not cigar:
        return ""

    chunks = list(iter_top_level_chunks(cigar))

    # --- Pre-fetch main reference ---
    main_ref_total = sum(
        op.n
        for _dir, text in chunks
        if ":" not in text
        for op in parse_pairwise_ops(text)
        if op.op in REF_OPS
    )
    main_ref_seq = ""
    if mainrefname and ref_lookup and main_ref_total > 0:
        main_rec = resolve_record(mainrefname, ref_lookup)
        if main_rec:
            if mainstrd == "+":
                s, e = mainrefpos, mainrefpos + main_ref_total
            else:
                e = main_rec.seqlen - mainrefpos
                s = max(0, e - main_ref_total)
            main_ref_seq = main_rec.fetch_local(s, e, "+")
            if mainstrd == "-":
                main_ref_seq = revcomp(main_ref_seq)

    parts: List[str] = []
    main_cursor = 0

    for direction, text in chunks:
        if ":" not in text:
            # ── Main-path segment ──────────────────────────────────────────
            for op in parse_pairwise_ops(text):
                if op.op in ("=", "X"):
                    if main_ref_seq:
                        parts.append(main_ref_seq[main_cursor:main_cursor + op.n])
                    main_cursor += op.n
                elif op.op == "I":
                    if op.payload:
                        parts.append(op.payload[:op.n])
                elif op.op == "D":
                    main_cursor += op.n
                # H / S: skip

        else:
            # ── Alt-path segment (independent reference) ───────────────────
            name, body = text.split(":", 1)
            rec = resolve_record(name, ref_lookup) if ref_lookup else None
            ops = parse_pairwise_ops(body)
            if not ops:
                continue

            lead_h = ops[0].n if ops[0].op == "H" else 0
            # Sum ref-consuming ops only (excludes H since H ∉ REF_OPS)
            ref_span = sum(op.n for op in ops if op.op in REF_OPS)

            # Pre-fetch alt reference region (between H-clips)
            alt_ref_seq = ""
            if rec and ref_span > 0:
                if direction == ">":
                    s, e = lead_h, lead_h + ref_span
                else:
                    trail_h = ops[-1].n if ops[-1].op == "H" else 0
                    e = rec.seqlen - trail_h
                    s = max(0, e - ref_span)
                alt_ref_seq = rec.fetch_local(s, e, "+")
                if direction == "<":
                    alt_ref_seq = revcomp(alt_ref_seq)

            # For reverse segments: invert op order and revcomp I payloads so
            # we can walk forward with a single cursor (same as main path).
            if direction == "<":
                ops = [
                    CigarOp(o.op, o.n, revcomp(o.payload) if o.payload else "")
                    for o in reversed(ops)
                ]

            alt_cursor = 0
            for op in ops:
                if op.op == "H":
                    continue  # already accounted for in pre-fetch bounds
                elif op.op in ("=", "X"):
                    if alt_ref_seq:
                        parts.append(alt_ref_seq[alt_cursor:alt_cursor + op.n])
                    alt_cursor += op.n
                elif op.op == "I":
                    if op.payload:
                        parts.append(op.payload[:op.n])
                elif op.op == "D":
                    alt_cursor += op.n
                # S: skip

    return "".join(parts)


def collect_encoded_reference_locations(encoded_cigar: str, ref_lookup: Dict[str, SeqRecord]) -> List[dict]:
    """Collect encoded INFO target locations in insertion-query order."""
    locs: List[dict] = []
    order = 0
    q_cursor = 0
    for direction, text in iter_top_level_chunks(encoded_cigar):
        if ":" not in text:
            q_cursor += encoded_query_len(text)
            continue
        name, body = text.split(":", 1)
        rec = resolve_record(name, ref_lookup)
        if rec is None:
            q_cursor += encoded_query_len(body)
            continue
        ops = parse_pairwise_ops(body)
        if not ops:
            continue
        lead_h, trail_h = segment_clip_lengths(ops)
        ref_span = sum(op.n for op in ops if op.op in REF_OPS)
        q_span = sum(op.n for op in ops if op.op in QUERY_OPS)
        if ref_span <= 0 or q_span <= 0:
            q_cursor += max(0, q_span)
            continue
        f0, f1 = segment_interval_to_forward(rec, direction, lead_h, trail_h, 0, ref_span)
        if f1 < f0:
            f0, f1 = f1, f0
        f0 = max(0, min(rec.seqlen, f0))
        f1 = max(f0, min(rec.seqlen, f1))
        if f1 <= f0:
            q_cursor += q_span
            continue
        seq = rec.fetch_local(f0, f1, "+")
        if has_non_acgt(seq):
            q_cursor += q_span
            continue
        if direction == "<":
            seq = revcomp(seq)
        locs.append({
            "kind": "encoded",
            "direction": direction,
            "name": name,
            "rec": rec,
            "start": f0,
            "end": f1,
            "seq": seq,
            "q_start": q_cursor,
            "q_end": q_cursor + q_span,
            "q_span": q_span,
            "match_bases": sum(op.n for op in ops if op.op in {"=", "M"}),
            "order": order,
        })
        q_cursor += q_span
        order += 1
    locs.sort(key=lambda x: (-int(x.get("match_bases", 0)), -int(x.get("q_span", 0)), int(x.get("order", 0))))
    return locs


def _alignment_value(aln, *names, default=None):
    if isinstance(aln, dict):
        for name in names:
            if name in aln:
                return aln[name]
    for name in names:
        if hasattr(aln, name):
            return getattr(aln, name)
    return default


def _alignment_cigar(aln) -> str:
    cigar = _alignment_value(aln, "CIGAR", "cigar", "cigar_string", default="")
    if cigar is None:
        return ""
    if isinstance(cigar, str):
        return cigar
    if isinstance(cigar, bytes):
        return cigar.decode()
    if isinstance(cigar, (list, tuple)):
        parts: List[str] = []
        for item in cigar:
            if isinstance(item, (list, tuple)) and len(item) >= 2:
                # ssw-py cigar_pair_list may be (op, n) or (n, op); accept both.
                a, b = item[0], item[1]
                if isinstance(a, str):
                    parts.append(f"{int(b)}{a}")
                else:
                    parts.append(f"{int(a)}{b}")
        return "".join(parts)
    return str(cigar)


def ssw_local_align(query_seq: str, target_seq: str) -> Optional[dict]:
    """Local-to-local align query to target with ssw-py.

    Returns query and target half-open coordinates plus a polished payloaded body.
    """
    if AlignmentMgr is None or not query_seq or not target_seq:
        return None
    mgr = AlignmentMgr(match_score=2, mismatch_penalty=2)
    mgr.set_read(query_seq.upper())
    mgr.set_reference(target_seq.upper())
    kwargs = {"gap_open": 3, "gap_extension": 1}
    if BitwiseAlignmentFlag is not None:
        kwargs["bitwise_flag"] = int(getattr(BitwiseAlignmentFlag, "best_idxs", 1))
    aln = mgr.align(**kwargs)

    cigar = _alignment_cigar(aln)
    score = _alignment_value(aln, "optimal_score", "score", default=0) or 0
    q0 = _alignment_value(aln, "read_start", "query_begin", "query_start", default=0)
    q1 = _alignment_value(aln, "read_end", "query_end", default=-1)
    t0 = _alignment_value(aln, "reference_start", "ref_begin", "target_start", default=0)
    t1 = _alignment_value(aln, "reference_end", "ref_end", "target_end", default=-1)
    q0 = int(q0)
    q1 = int(q1) + 1
    t0 = int(t0)
    t1 = int(t1) + 1
    score = float(score)
    if not cigar or score <= 0 or q1 <= q0 or t1 <= t0:
        return None
    q0 = max(0, min(len(query_seq), q0))
    q1 = max(q0, min(len(query_seq), q1))
    t0 = max(0, min(len(target_seq), t0))
    t1 = max(t0, min(len(target_seq), t1))
    if q1 <= q0 or t1 <= t0:
        return None
    q_sub = query_seq[q0:q1]
    t_sub = target_seq[t0:t1]
    body = cigarextend(cigar, t_sub, q_sub)
    aligned, sim = summarize_extend_body(body)
    return {
        "q0": q0,
        "q1": q1,
        "target_t0": t0,
        "target_t1": t1,
        "cigar": cigar,
        "body": body,
        "score": score,
        "aligned": aligned,
        "similarity": sim,
    }


def minimap2_local_align(query_seq: str, target_seq: str) -> Optional[dict]:
    """Local align query to target using mappy (minimap2 Python bindings).

    Used for sequences > 1000 bp where SSW O(m*n) DP is too slow.
    Returns the same dict format as ssw_local_align.
    """
    if mappy is None or not query_seq or not target_seq:
        return None
    try:
        aligner = mappy.Aligner(seq=target_seq, preset="asm10", best_n=1)
    except Exception as exc:
        return None
    hit = None
    for h in aligner.map(query_seq):
        hit = h
        break
    if hit is None:
        return None
    q0 = int(hit.q_st)
    q1 = int(hit.q_en)
    t0 = int(hit.r_st)
    t1 = int(hit.r_en)
    cigar = hit.cigar_str

    if not cigar or q1 <= q0 or t1 <= t0:
        return None
    q_sub = query_seq[q0:q1]
    t_sub = target_seq[t0:t1]
    
    body = cigarextend(cigar, t_sub, q_sub)
    aligned, sim = summarize_extend_body(body)
    score = float(hit.mlen)
    return {
        "q0": q0,
        "q1": q1,
        "target_t0": t0,
        "target_t1": t1,
        "cigar": cigar,
        "body": body,
        "score": score,
        "aligned": aligned,
        "similarity": sim,
    }


def _minimap2_insertion_alignment(
    query_seq: str,
    target_seq: str,
    aligner=None,
) -> Optional[dict]:
    """Return the best forward minimap2 alignment for two insertion alleles.

    Insertions in one reference-coordinate cluster are already oriented to the
    same reference strand.  Reverse hits are therefore not valid encodings for
    the sample graph CIGAR and are ignored.
    """
    if mappy is None or not query_seq or not target_seq:
        return None
    if aligner is None:
        try:
            aligner = _make_insertion_minimap2_aligner(target_seq)
        except Exception:
            return None
    candidates = []
    try:
        for hit in aligner.map(query_seq):
            if int(getattr(hit, "strand", 1)) != 1:
                continue
            q0 = int(hit.q_st)
            q1 = int(hit.q_en)
            t0 = int(hit.r_st)
            t1 = int(hit.r_en)
            cigar = hit.cigar_str
            if not cigar or q1 <= q0 or t1 <= t0:
                continue
            q_sub = query_seq[q0:q1]
            t_sub = target_seq[t0:t1]
            body = cigarextend(cigar, t_sub, q_sub)
            aligned, similarity = summarize_extend_body(body)
            if aligned <= 0:
                continue
            candidates.append({
                "q0": q0,
                "q1": q1,
                "target_t0": t0,
                "target_t1": t1,
                "cigar": cigar,
                "body": body,
                "score": float(getattr(hit, "mlen", aligned)),
                "aligned": aligned,
                "similarity": similarity,
                "mapq": int(getattr(hit, "mapq", 0)),
            })
    except Exception:
        return None
    if not candidates:
        return None
    return max(
        candidates,
        key=lambda value: (
            int(value["aligned"]),
            float(value["similarity"]),
            int(value["mapq"]),
            float(value["score"]),
        ),
    )


def _make_insertion_minimap2_aligner(target_sequence: str):
    """Build a minimap2 index suitable for the current recursion level.

    The assembly preset intentionally suppresses very short mappings, whereas
    --var-in-insert accepts a bare threshold of 100 bp.  Minimap2's short-read
    preset keeps those nested 100--1000 bp comparisons callable; longer
    insertion targets retain the assembly preset.
    """
    if mappy is None:
        raise RuntimeError("mappy is unavailable")
    preset = "sr" if len(target_sequence) < 1000 else "asm10"
    return mappy.Aligner(seq=target_sequence, preset=preset, best_n=5)


def _slice_insertion_target(
    target: InsertionSequenceTarget,
    start: int,
    end: int,
) -> InsertionSequenceTarget:
    """Address an oriented subsequence in the same source query path."""
    start = max(0, min(len(target.sequence), int(start)))
    end = max(start, min(len(target.sequence), int(end)))
    if target.direction == "<":
        local_start = target.local_end - end
        local_end = target.local_end - start
    else:
        local_start = target.local_start + start
        local_end = target.local_start + end
    return InsertionSequenceTarget(
        sv=target.sv,
        sequence=target.sequence[start:end],
        path_name=target.path_name,
        record_length=target.record_length,
        local_start=local_start,
        local_end=local_end,
        direction=target.direction,
    )


def _target_chunk(
    target: InsertionSequenceTarget,
    start: int,
    end: int,
    body: str,
) -> str:
    """Render one named query-path chunk for an oriented target subinterval."""
    if not target.path_name or not body or end <= start:
        return ""
    start = max(0, min(len(target.sequence), int(start)))
    end = max(start, min(len(target.sequence), int(end)))
    if target.direction == "<":
        forward_start = target.local_end - end
        forward_end = target.local_end - start
        leading_h = max(0, target.record_length - forward_end)
        trailing_h = max(0, forward_start)
    else:
        forward_start = target.local_start + start
        forward_end = target.local_start + end
        leading_h = max(0, forward_start)
        trailing_h = max(0, target.record_length - forward_end)
    return (
        f"{target.direction}{target.path_name}:"
        f"{_hclip_text(leading_h)}{body}{_hclip_text(trailing_h)}"
    )


def _render_recursive_insertion_parts(parts: Sequence[object]) -> str:
    output: List[str] = []
    for part in parts:
        if isinstance(part, str):
            output.append(part)
        elif isinstance(part, RecursiveInsertionPiece):
            if part.replacement is None:
                output.append(encode_raw_piece_as_insertion(part.sequence))
            else:
                output.append(_render_recursive_insertion_parts(part.replacement))
    return "".join(output)


def _alignment_parts_against_insertion_target(
    query: InsertionSequenceTarget,
    target: InsertionSequenceTarget,
    ordinal: List[int],
    aligner=None,
) -> Tuple[Optional[List[object]], List[RecursiveInsertionPiece]]:
    """Build graph-CIGAR pieces and expose residual I segments for recursion."""
    alignment = _minimap2_insertion_alignment(
        query.sequence,
        target.sequence,
        aligner=aligner,
    )
    if alignment is None or not target.path_name:
        return None, []

    q0 = int(alignment["q0"])
    q1 = int(alignment["q1"])
    t0 = int(alignment["target_t0"])
    t1 = int(alignment["target_t1"])
    target_length = len(target.sequence)
    parts: List[object] = []
    residuals: List[RecursiveInsertionPiece] = []

    def add_residual(start: int, end: int, anchor: int) -> None:
        if end <= start:
            return
        ordinal[0] += 1
        source = _slice_insertion_target(query, start, end)
        piece = RecursiveInsertionPiece(
            sequence=source.sequence,
            source=source,
            anchor=max(0, int(anchor)),
            ordinal=ordinal[0],
        )
        parts.append(piece)
        residuals.append(piece)

    add_residual(0, q0, t0)
    qpos = q0
    tpos = t0
    block_start = 0
    block_ops: List[str] = []
    if t0 > 0:
        block_ops.append(f"{t0}D")

    def flush_block() -> None:
        nonlocal block_start, block_ops
        if not block_ops or tpos <= block_start:
            block_start = tpos
            block_ops = []
            return
        chunk = _target_chunk(
            target,
            block_start,
            tpos,
            "".join(block_ops),
        )
        if chunk:
            parts.append(chunk)
        block_start = tpos
        block_ops = []

    for op in parse_pairwise_ops(str(alignment["body"])):
        if op.n <= 0:
            continue
        if op.op == "I":
            flush_block()
            add_residual(qpos, qpos + op.n, tpos)
            qpos += op.n
            block_start = tpos
            continue
        payload = op.payload[:op.n] if op.payload else ""
        block_ops.append(
            f"{op.n}{op.op}{payload}" if payload else f"{op.n}{op.op}"
        )
        if op.op in QUERY_OPS:
            qpos += op.n
        if op.op in REF_OPS or op.op == "H":
            tpos += op.n

    # Mappy's local target end is authoritative.  Add any target suffix as a
    # deletion within this insertion comparison rather than silently clipping
    # it away.
    tpos = max(tpos, t1)
    if target_length > tpos:
        block_ops.append(f"{target_length - tpos}D")
        tpos = target_length
    flush_block()
    add_residual(q1, len(query.sequence), t1)
    if not any(isinstance(part, str) for part in parts):
        return None, []
    return parts, residuals


def _group_recursive_insertion_pieces(
    pieces: Sequence[RecursiveInsertionPiece],
) -> List[List[RecursiveInsertionPiece]]:
    """Group residual insertions at the same or minimally shifted anchor."""
    ordered = sorted(pieces, key=lambda piece: (
        piece.anchor, -piece.unmasked_bases, piece.ordinal,
    ))
    groups: List[List[RecursiveInsertionPiece]] = []
    for piece in ordered:
        if not groups:
            groups.append([piece])
            continue
        previous = groups[-1][-1]
        distance_limit = max(
            1,
            min(
                500,
                int(0.1 * max(
                    1,
                    min(previous.unmasked_bases, piece.unmasked_bases),
                )),
            ),
        )
        if piece.anchor - previous.anchor <= distance_limit:
            groups[-1].append(piece)
        else:
            groups.append([piece])
    return groups


def _recursively_encode_insertion_pieces(
    pieces: Sequence[RecursiveInsertionPiece],
    min_unmasked_bases: int,
    ordinal: List[int],
    depth: int = 0,
) -> None:
    """Recursively align supported residual insertions to their largest member."""
    if depth >= 32:
        return
    for group in _group_recursive_insertion_pieces(pieces):
        if len(group) < 2:
            continue
        representative = max(
            group,
            key=lambda piece: (
                piece.unmasked_bases,
                len(piece.sequence),
                -piece.ordinal,
            ),
        )
        if representative.unmasked_bases <= int(min_unmasked_bases):
            continue
        target = representative.source
        if not target.path_name:
            continue
        try:
            aligner = _make_insertion_minimap2_aligner(target.sequence)
        except Exception:
            continue
        nested: List[RecursiveInsertionPiece] = []
        for piece in group:
            if piece is representative:
                continue
            replacement, residuals = _alignment_parts_against_insertion_target(
                piece.source,
                target,
                ordinal,
                aligner=aligner,
            )
            if replacement is None:
                continue
            piece.replacement = replacement
            nested.extend(residuals)
        if nested:
            _recursively_encode_insertion_pieces(
                nested,
                min_unmasked_bases,
                ordinal,
                depth + 1,
            )


def recursive_insertion_cigar_overrides(
    members: Sequence[SVUnit],
    targets: Mapping[int, InsertionSequenceTarget],
    representative: InsertionSequenceTarget,
    min_unmasked_bases: int,
    fallback_target_name: str,
) -> Dict[int, str]:
    """Align all other observations to the largest insertion, recursively."""
    if (
        int(min_unmasked_bases) <= 0
        or len(members) < 2
        or representative.unmasked_bases <= int(min_unmasked_bases)
        or mappy is None
    ):
        return {}
    root_target = representative
    if not root_target.path_name:
        root_target = InsertionSequenceTarget(
            sv=representative.sv,
            sequence=representative.sequence,
            path_name=fallback_target_name,
            record_length=len(representative.sequence),
            local_start=0,
            local_end=len(representative.sequence),
            direction=">",
        )
    try:
        aligner = _make_insertion_minimap2_aligner(root_target.sequence)
    except Exception:
        return {}

    ordinal = [0]
    parts_by_id: Dict[int, List[object]] = {}
    residuals: List[RecursiveInsertionPiece] = []
    for member in members:
        if member is representative.sv:
            continue
        query = targets[member.id]
        if not query.sequence or any(
            base.upper() not in "ACGT" for base in query.sequence
        ):
            continue
        parts, member_residuals = _alignment_parts_against_insertion_target(
            query,
            root_target,
            ordinal,
            aligner=aligner,
        )
        if parts is None:
            continue
        parts_by_id[member.id] = parts
        residuals.extend(member_residuals)
    if residuals:
        _recursively_encode_insertion_pieces(
            residuals,
            int(min_unmasked_bases),
            ordinal,
        )
    return {
        member_id: _render_recursive_insertion_parts(parts)
        for member_id, parts in parts_by_id.items()
    }


def _hclip_text(n: int) -> str:
    return f"{n}H" if n > 0 else ""


def _dclip_text(n: int) -> str:
    return f"{n}D" if n > 0 else ""


def addH_raw(target_len: int, refstart: int, refend: int, cigar_body: str) -> str:
    refstart = max(0, int(refstart))
    refend = max(refstart, int(refend))
    left = refstart
    right = max(0, int(target_len) - refend)
    return ">" + _hclip_text(left) + cigar_body + _hclip_text(right)


def encoded_forward_interval(loc: dict, oriented_start: int, oriented_end: int) -> Tuple[int, int]:
    a = int(oriented_start)
    b = int(oriented_end)
    if b < a:
        a, b = b, a
    if loc.get("direction") == "<":
        f0 = int(loc["end"]) - b
        f1 = int(loc["end"]) - a
    else:
        f0 = int(loc["start"]) + a
        f1 = int(loc["start"]) + b
    if f1 < f0:
        f0, f1 = f1, f0
    return f0, f1


def addH_encoded(loc: dict, refstart: int, refend: int, cigar_body: str) -> str:
    rec: SeqRecord = loc["rec"]
    f0, f1 = encoded_forward_interval(loc, refstart, refend)
    f0 = max(0, min(rec.seqlen, f0))
    f1 = max(f0, min(rec.seqlen, f1))
    direction = loc.get("direction", ">")
    if direction == "<":
        left = max(0, rec.seqlen - f1)
        right = max(0, f0)
    else:
        left = max(0, f0)
        right = max(0, rec.seqlen - f1)
    return f"{direction}{loc.get('name', '')}:" + _hclip_text(left) + cigar_body + _hclip_text(right)


def encode_raw_piece_as_insertion(seq: str) -> str:
    return f">{len(seq)}I{seq}" if seq else ""


def count_raw_unaligned_bp(encoded: str) -> int:
    """Count query bases that are insertions (I-ops) not placed against any reference.

    Both raw chunks (no ':') and named-segment chunks (have ':') are scanned for
    I-ops only.  '=' and 'X' ops always represent bases aligned to a reference —
    either the named segment's reference or the main alignment reference — and are
    NOT counted as unaligned.
    """
    total = 0
    for _direction, text in iter_top_level_chunks(encoded):
        body = text.split(":", 1)[1] if ":" in text else text
        for n_s, op in re.findall(r"(\d+)([MIDHS=X])", body):
            if op == "I":
                total += int(n_s)
    return total


def subtract_assigned_from_remaining(remaining: Sequence[Tuple[int, int]], assigned: Sequence[Tuple[int, int]]) -> List[Tuple[int, int]]:
    return subtract_intervals(list(remaining), list(assigned))


def build_alignment_candidates_for_target(raw_seq: str, remaining: Sequence[Tuple[int, int]], target: dict,
                                          min_align: int, sim_min: float) -> List[dict]:
    candidates: List[dict] = []
    target_seq = target.get("seq", "")
    if not target_seq:
        return candidates
    for a, b in remaining:
        if b - a <= min_align:
            continue
        qpiece = raw_seq[a:b]
        use_minimap2 = len(qpiece) > 1000 or len(target_seq) > 1000
        if use_minimap2 and mappy is None:
            if len(qpiece) > REALIGN_MAX_SSW_SEQ_LEN or len(target_seq) > REALIGN_MAX_SSW_SEQ_LEN:
                continue
        aln = minimap2_local_align(qpiece, target_seq) if use_minimap2 else ssw_local_align(qpiece, target_seq)
        if aln is None:
            continue
        full_q0 = a + int(aln["q0"])
        full_q1 = a + int(aln["q1"])
        if full_q1 - full_q0 <= min_align:
            continue
        if int(aln["aligned"]) <= min_align or float(aln["similarity"]) < sim_min:
            continue
        cand = dict(aln)
        cand.update({
            "kind": target.get("kind"),
            "target": target,
            "q0": full_q0,
            "q1": full_q1,
        })
        candidates.append(cand)
    return candidates


def render_assigned_alignment(cand: dict, qstart: int, qend: int, raw_seq: str) -> Tuple[int, int, str]:
    refstart_rel, refend_rel, cut_cigar = cigar_findqrange(str(cand["cigar"]), qstart, qend, int(cand["q0"]))
    if not cut_cigar or refstart_rel < 0 or refend_rel < 0:
        return qstart, qend, encode_raw_piece_as_insertion(raw_seq[qstart:qend])
    abs_refstart = int(cand["target_t0"]) + int(refstart_rel)
    abs_refend = int(cand["target_t0"]) + int(refend_rel)
    target_seq = cand["target"].get("seq", "")
    query_sub = raw_seq[qstart:qend]
    ref_sub = target_seq[abs_refstart:abs_refend]
    polished = cigarextend(cut_cigar, ref_sub, query_sub)
    if cand.get("kind") == "encoded":
        encoded = addH_encoded(cand["target"], abs_refstart, abs_refend, polished)
        f0, f1 = encoded_forward_interval(cand["target"], abs_refstart, abs_refend)
    else:
        encoded = addH_raw(len(target_seq), abs_refstart, abs_refend, polished)
    return qstart, qend, encoded


def render_alignment_with_raw_flanks(cand: dict, raw_seq: str) -> str:
    """Render one local INFO alignment plus raw prefix/suffix outside its query span."""
    q0 = max(0, min(len(raw_seq), int(cand.get("q0", 0))))
    q1 = max(q0, min(len(raw_seq), int(cand.get("q1", q0))))
    out: List[str] = []
    if q0 > 0:
        out.append(encode_raw_piece_as_insertion(raw_seq[:q0]))
    if q1 > q0:
        target = cand.get("target", {})
        target_seq = target.get("seq", "")
        body = str(cand.get("body", ""))
        t0 = max(0, min(len(target_seq), int(cand.get("target_t0", 0))))
        t1 = max(t0, min(len(target_seq), int(cand.get("target_t1", t0))))
        if body and target.get("kind") == "encoded":
            direction = target.get("direction", ">")
            left = max(0, t0)
            right = max(0, len(target_seq) - t1)
            out.append(f"{direction}{target.get('name', '')}:" + _dclip_text(left) + body + _dclip_text(right))
        elif body:
            left = max(0, t0)
            right = max(0, len(target_seq) - t1)
            out.append(">" + _dclip_text(left) + body + _dclip_text(right))
        else:
            _a, _b, enc = render_assigned_alignment(cand, q0, q1, raw_seq)
            out.append(enc)
    if q1 < len(raw_seq):
        out.append(encode_raw_piece_as_insertion(raw_seq[q1:]))
    return "".join(x for x in out if x)


def segmented_encode_raw_insertion(raw_seq: str, encoded_locs: Sequence[dict], raw_target_seq: str,
                                   cutoff: int = REALIGN_MIN_ALIGNED_BP,
                                   sim_min: float = REALIGN_MIN_SIMILARITY,
                                   min_align: int = REALIGN_MIN_ALIGNED_BP) -> str:
    """Encode raw insertion sequence against INFO target(s), preserving full query coverage.

    Aligned pieces use the INFO target locations and H clips.  Unaligned pieces
    are emitted as >{size}I{rawseq}.  If no alignment is accepted, return empty
    so the caller keeps the original plain sample CIGAR.
    """
    if not raw_seq:
        return ""
    if AlignmentMgr is None and mappy is None:
        return ""

    targets: List[dict] = []
    for loc in encoded_locs:
        if loc.get("seq"):
            targets.append(loc)
    if not targets and raw_target_seq:
        targets.append({"kind": "raw", "seq": raw_target_seq, "score": len(raw_target_seq)})

    if not targets:
        return ""

    remaining: List[Tuple[int, int]] = [(0, len(raw_seq))]
    rendered: List[Tuple[int, int, str]] = []

    for target in targets:
        while True:
            candidates = build_alignment_candidates_for_target(raw_seq, remaining, target, min_align, sim_min)
            if not candidates:
                break
            alignments = [[c["q0"], c["q1"], c["score"], i] for i, c in enumerate(candidates)]
            assignments = graphicalign(alignments, allow_gap=cutoff)
            if not assignments:
                break
            assigned_intervals: List[Tuple[int, int]] = []
            for cand_i, q0, q1 in assignments:
                cand = candidates[int(cand_i)]
                if q1 <= q0:
                    continue
                rendered.append(render_assigned_alignment(cand, int(q0), int(q1), raw_seq))
                assigned_intervals.append((int(q0), int(q1)))
            if not assigned_intervals:
                break
            new_remaining = subtract_assigned_from_remaining(remaining, assigned_intervals)
            if new_remaining == remaining:
                break
            remaining = new_remaining

    if not rendered:
        return ""

    for a, b in remaining:
        if b > a:
            raw_piece = encode_raw_piece_as_insertion(raw_seq[a:b])
            rendered.append((a, b, raw_piece))

    rendered.sort(key=lambda x: (x[0], x[1]))
    # Sanity: encoded/raw pieces must partition the entire query interval.
    cursor = 0
    out: List[str] = []
    for a, b, txt in rendered:
        if a > cursor:
            raw_piece = encode_raw_piece_as_insertion(raw_seq[cursor:a])
            out.append(raw_piece)
        elif a < cursor:
            # This should not happen after graphicalign/subtraction; skip the overlapped prefix safely.
            trim = cursor - a
            if b <= cursor:
                continue
            # If an overlap appears, fall back to raw for the unconsumed suffix.
            txt = encode_raw_piece_as_insertion(raw_seq[cursor:b])
            a = cursor
        out.append(txt)
        cursor = max(cursor, b)
    if cursor < len(raw_seq):
        raw_piece = encode_raw_piece_as_insertion(raw_seq[cursor:])
        out.append(raw_piece)
    encoded_out = "".join(x for x in out if x)
    return encoded_out

def _best_encoded_for_allele(
    raw_seq: str,
    encoded_locs: List[dict],
    info_cigar: str,
    info_raw_seq: str,
    current_unaligned: int,
    target_name: str = "INFO",
) -> Tuple[str, bool]:
    """Return (best_encoding, skip_row_prefix) or ('', False) if no improvement."""
    if not raw_seq:
        return "", False
    best, best_u = "", current_unaligned
    best_is_info_cigar = False
    
    if info_cigar and encoded_query_len(info_cigar) == len(raw_seq):
        u_info = count_raw_unaligned_bp(info_cigar)
        if u_info < best_u:
            best, best_u = info_cigar, u_info
            best_is_info_cigar = True

    info_raw_exact = bool(info_raw_seq) and raw_seq.upper() == info_raw_seq.upper()
    if info_raw_exact:
        if (
            encoded_locs and info_cigar
            and encoded_query_len(info_cigar) == len(raw_seq)
        ):
            u_info = count_raw_unaligned_bp(info_cigar)
            if u_info <= best_u:
                return info_cigar, True
        exact = f">{target_name}:{len(raw_seq)}="
        if count_raw_unaligned_bp(exact) < best_u:
            return exact, False
    
    if AlignmentMgr is None and mappy is None:
        return best, best_is_info_cigar
            
    if encoded_locs:
        # Thresholds: cutoff=200, sim_min=0.0 (per request), min_align=200
        enc = segmented_encode_raw_insertion(raw_seq, encoded_locs, "", cutoff=200, sim_min=0.0, min_align=200)
        u = count_raw_unaligned_bp(enc) if enc else -1
        if enc and u < best_u:
            best, best_u = enc, u
            best_is_info_cigar = False
            
    if info_raw_seq:
        min_align = max(1, int(0.50 * len(raw_seq)))
        dummy_rec = SeqRecord(name=target_name, allelename=target_name, contig=target_name, start=0, end=len(info_raw_seq), seqlen=len(info_raw_seq))
        target = {"kind": "encoded", "direction": ">", "name": target_name, "rec": dummy_rec, "start": 0, "end": len(info_raw_seq), "seq": info_raw_seq}
        
        cands = build_alignment_candidates_for_target(raw_seq, [(0, len(raw_seq))], target, min_align=min_align, sim_min=0.0)
        if cands:
            best_cand = max(cands, key=lambda c: c["aligned"])
            enc_raw = render_alignment_with_raw_flanks(best_cand, raw_seq)
            u_raw = count_raw_unaligned_bp(enc_raw)
            if enc_raw and u_raw < best_u:
                best, best_u = enc_raw, u_raw
                best_is_info_cigar = False
                    
    if best and encoded_query_len(best) != len(raw_seq):
        raise ValueError(
            "candidate insertion encoding does not consume the complete "
            f"query allele: {encoded_query_len(best)} versus {len(raw_seq)}"
        )
    return best, best_is_info_cigar

# ---------------------------------------------------------------------------
# Optional per-line encoding pass
# ---------------------------------------------------------------------------


# Interval index: contig/allelename → all SeqRecords for that contig, for
# covering-range lookup when multiple FASTA segments share the same contig name.
_QUERY_INTERVAL_INDEX: Dict[str, List[SeqRecord]] = {}
# The lookup dict _QUERY_INTERVAL_INDEX was built from.  When a caller passes
# this exact dict, the index provably contains every record the covering-record
# acceptance check could match, so the O(n) fallback scan can be skipped.
_QUERY_INTERVAL_INDEX_SOURCE: Optional[Dict[str, SeqRecord]] = None
_REENCODE_WORKER_CONTEXT = {}


def build_contig_interval_index(lookup: Dict[str, SeqRecord]) -> Dict[str, List[SeqRecord]]:
    """Build contig→[SeqRecord, ...] index over all records in lookup.

    Used to find a specific FASTA segment that fully covers a requested
    genomic range when multiple segments share the same contig/allelename.
    """
    index: Dict[str, List[SeqRecord]] = {}
    seen: set = set()
    for rec in lookup.values():
        rid = id(rec)
        if rid in seen:
            continue
        seen.add(rid)
        for key in {rec.name, strip_match_name(rec.name), rec.contig, strip_match_name(rec.contig), rec.allelename, strip_match_name(rec.allelename)}:
            if key:
                index.setdefault(key, []).append(rec)
    return index


def _init_reencode_worker(context: dict) -> None:
    """Install read-only sequence state once in each realignment worker."""
    global _REENCODE_WORKER_CONTEXT, _QUERY_INTERVAL_INDEX
    global _QUERY_INTERVAL_INDEX_SOURCE
    _REENCODE_WORKER_CONTEXT = context
    _QUERY_INTERVAL_INDEX = context["query_interval_index"]
    _QUERY_INTERVAL_INDEX_SOURCE = context["query_lookup"]


def _reencode_vcf_line_worker(line: str) -> str:
    context = _REENCODE_WORKER_CONTEXT
    if not context:
        raise RuntimeError("realignment worker context is not initialized")
    return reencode_vcf_line(
        line,
        context["ref_lookup"],
        context["force_realign"],
        context["query_lookup"],
    )


def _find_covering_record(contig: str, g0: int, g1: int, lookup: Dict[str, SeqRecord]) -> Optional[SeqRecord]:
    """Return the SeqRecord for contig whose genomic interval covers [g0, g1]."""
    candidates: List[SeqRecord] = []
    seen: set = set()

    def add_candidate(r: Optional[SeqRecord]) -> None:
        if r is None or id(r) in seen:
            return
        if r.start <= g0 and r.end >= g1 and (
            r.name == contig
            or strip_match_name(r.name) == contig
            or r.contig == contig
            or strip_match_name(r.contig) == contig
            or r.allelename == contig
            or strip_match_name(r.allelename) == contig
        ):
            seen.add(id(r))
            candidates.append(r)

    # Header-index path: contig/name/allelename -> all candidate FASTA intervals.
    for key in {contig, strip_match_name(contig)}:
        for r in _QUERY_INTERVAL_INDEX.get(key, []):
            add_candidate(r)
    add_candidate(resolve_record(contig, lookup))
    # O(n) last-resort fallback, skipped when the active index was built from
    # this exact lookup: the index is keyed by all six identity fields the
    # acceptance check above compares, so it already yielded every possible
    # candidate and the scan cannot add more.
    if lookup is not _QUERY_INTERVAL_INDEX_SOURCE:
        for r in lookup.values():
            add_candidate(r)
    if not candidates:
        return None
    return min(candidates, key=lambda r: (r.end - r.start, r.name))


def reencode_vcf_line(line: str, ref_lookup: Optional[Dict[str, SeqRecord]] = None, force_realign: int = 200, query_lookup: Optional[Dict[str, SeqRecord]] = None) -> str:
    if ref_lookup is None: ref_lookup = {}
    if query_lookup is None: query_lookup = ref_lookup
    if line.startswith("#"): return line
    cols = line.rstrip("\n").split("\t")
    if len(cols) < 10 or not is_sv_format(cols[8]): return line
    format_text = cols[8]
    
    info_order = [p.split("=", 1)[0] for p in cols[7].split(";") if p]
    info = parse_info_field(cols[7])
    
    # QUERYSEQ is an internal coordinate carrier used by the final VCF schema
    # pass to materialize INFO/SEQ.  Keep it until that pass has fetched the
    # authoritative query bases.
    queryseq_coord = info.get("QUERYSEQ", "") or ""
    has_queryseq = bool(queryseq_coord)
    
    info_cigar = info.get("CIGAR", "")
    # Check for total insertion size in CIGAR; if <= 50, stop work but return the line
    if count_raw_unaligned_bp(info_cigar) <= 50:
        if has_queryseq:
            cols[7] = format_info_field(info, info_order)
            return "\t".join(cols) + "\n"
        return line
    
    encoded_locs = collect_encoded_reference_locations(info_cigar, ref_lookup) if info_cigar_has_encoded_target(info_cigar) else []
    
    info_raw_from_fasta = ""
    if queryseq_coord and "|" in queryseq_coord:
        qname_esc, qrange = queryseq_coord.split("|", 1)
        info_raw_from_fasta, _, _, _, _ = fetch_allele_query_sequence_with_record(["", "", "", "", qname_esc, qrange], query_lookup)
    # Fallback: when QUERYSEQ is absent, the INFO CIGAR itself may carry the raw
    # insertion sequence inline (e.g. >10764I{GGCAG...}).  Extract it so that
    # PRIORITY 2 in _best_encoded_for_allele can align sample alleles against it.
    if not info_raw_from_fasta and info_cigar:
        info_raw_from_fasta = extract_raw_inserted_sequence_from_cigar(info_cigar)

    # Pre-scan allele work
    has_work = False
    for field in cols[9:]:
        for allele in split_sample_alleles(field, format_text):
            vals = parse_hsv_allele(allele)
            if vals is not None and vals[0] == "1" and "I" in vals[3]:
                _pad, core = split_row_position_padding(vals[3])
                is_named = sample_cigar_has_named_segment(core)
                if (is_named and count_raw_unaligned_bp(core) >= force_realign) or (not is_named):
                    has_work = True; break
        if has_work: break
        
    if not has_work or (AlignmentMgr is None and mappy is None) or not info_cigar:
        cols[7] = format_info_field(info, info_order)
        return "\t".join(cols) + "\n"
    
    for ci in range(9, len(cols)):
        field = cols[ci]
        alleles = split_sample_alleles(field, format_text)
        if not alleles:
            continue
        parts_out = []
        for ai, allele in enumerate(alleles):
            vals = parse_hsv_allele(allele)
            if vals is None or vals[0] != "1" or "I" not in vals[3]:
                parts_out.append(allele); continue
            
            pad_prefix, core_cigar = split_row_position_padding(vals[3])
            is_named = sample_cigar_has_named_segment(core_cigar)
            if is_named and count_raw_unaligned_bp(core_cigar) < force_realign:
                parts_out.append(allele); continue
                
            raw_seq = fetch_allele_query_sequence(vals, query_lookup)
            if not raw_seq:
                raw_seq = extract_raw_inserted_sequence_from_cigar(core_cigar)
            if not raw_seq or "N" in raw_seq.upper():
                parts_out.append(allele); continue

            query_span = query_range_len(vals[5])
            if query_span != len(raw_seq):
                raise ValueError(
                    f"{cols[0]}:{cols[1]}:{cols[2]} sample column {ci - 8} "
                    f"allele {ai + 1}: QUERYCOORD span={query_span} but "
                    f"sequence length={len(raw_seq)}"
                )
                
            current_unaligned = count_raw_unaligned_bp(core_cigar) if is_named else len(raw_seq)
            old_query_len = encoded_query_len(core_cigar)
            if old_query_len != len(raw_seq):
                raise ValueError(
                    f"{cols[0]}:{cols[1]}:{cols[2]} sample column {ci - 8} "
                    f"allele {ai + 1}: EXTENDGRAPHCIGAR consumes "
                    f"{old_query_len} query bases but sequence length is "
                    f"{len(raw_seq)}"
                )
            
            encoded, skip_prefix = _best_encoded_for_allele(
                            raw_seq, encoded_locs, info_cigar, info_raw_from_fasta, 
                            current_unaligned, target_name=cols[2]
                        )
            
            if encoded:
                if encoded_query_len(encoded) != len(raw_seq):
                    raise ValueError(
                        f"{cols[0]}:{cols[1]}:{cols[2]} sample column "
                        f"{ci - 8} allele {ai + 1}: re-encoded CIGAR consumes "
                        f"{encoded_query_len(encoded)} query bases, expected "
                        f"{len(raw_seq)}"
                    )
                vals[3] = encoded if skip_prefix else pad_prefix + encoded
            parts_out.append(format_hsv_allele(vals))
        cols[ci] = format_sample_alleles(parts_out, format_text)
        
    cols[7] = format_info_field(info, info_order)
    return "\t".join(cols) + "\n"



def line_level_reencode_vcf(tmp_path: str, out_path: str, realignment: bool, ref_lookup: Optional[Dict[str, SeqRecord]] = None, query_lookup: Optional[Dict[str, SeqRecord]] = None, processes: int = 1, force_realign: int = 200):
    """Stream final VCF and optionally re-encode insertion CIGARs in parallel."""
    out_dir = os.path.dirname(os.path.abspath(out_path)) or "."
    fd, stage_path = tempfile.mkstemp(
        prefix=os.path.basename(out_path) + ".",
        suffix=".reencode.tmp",
        dir=out_dir,
    )
    os.close(fd)

    try:
        if not realignment:
            with open(tmp_path) as src, open(stage_path, "w") as dst:
                shutil.copyfileobj(src, dst, length=1024 * 1024)
            os.replace(stage_path, out_path)
            return

        # Without an aligner the re-encoding operation is a no-op. Copying in
        # large blocks avoids dispatching hundreds of thousands of unchanged
        # VCF lines through Python workers.
        if AlignmentMgr is None and mappy is None:
            warned_no_aligner = False
            with open(tmp_path) as src, open(stage_path, "w") as dst:
                for line in src:
                    dst.write(line)
                    if not warned_no_aligner and not line.startswith("#"):
                        print(
                            "[WARN] --realignment requested but no supported "
                            "aligner is importable; keeping raw sample insertion cigars",
                            file=sys.stderr,
                        )
                        warned_no_aligner = True
            os.replace(stage_path, out_path)
            return

        ref_lookup = ref_lookup or {}
        query_lookup = query_lookup or {}
        context = {
            "ref_lookup": ref_lookup,
            "query_lookup": query_lookup,
            "query_interval_index": build_contig_interval_index(query_lookup),
            "force_realign": force_realign,
        }
        global _REENCODE_WORKER_CONTEXT, _QUERY_INTERVAL_INDEX
        global _QUERY_INTERVAL_INDEX_SOURCE
        _REENCODE_WORKER_CONTEXT = context
        _QUERY_INTERVAL_INDEX = context["query_interval_index"]
        _QUERY_INTERVAL_INDEX_SOURCE = context["query_lookup"]

        completed = 0
        next_report = 100_000

        def consume(results, destination) -> None:
            nonlocal completed, next_report
            for result in results:
                destination.write(result)
                completed += 1
                if completed >= next_report:
                    print(
                        f"[graphreftovcf] realigned {completed} VCF rows",
                        file=sys.stderr,
                    )
                    while next_report <= completed:
                        next_report += 100_000

        with open(tmp_path) as src, open(stage_path, "w") as dst:
            # Headers need no re-encoding. Write them directly before the pool
            # begins and submit only variant records. imap preserves input
            # order, so parallel work cannot reorder the final VCF.
            first_data_line = None
            for line in src:
                if line.startswith("#"):
                    dst.write(line)
                else:
                    first_data_line = line
                    break

            def data_lines():
                if first_data_line is not None:
                    yield first_data_line
                yield from src

            worker_count = max(1, int(processes))
            print(
                f"[graphreftovcf] realigning VCF rows with "
                f"{worker_count} process(es) ...",
                file=sys.stderr,
            )
            if worker_count > 1:
                _freeze_shared_heap()
                with mp_context().Pool(
                    processes=worker_count,
                    initializer=_init_reencode_worker,
                    initargs=(context,),
                    maxtasksperchild=64,
                ) as pool:
                    consume(
                        pool.imap(
                            _reencode_vcf_line_worker,
                            data_lines(),
                            chunksize=32,
                        ),
                        dst,
                    )
            else:
                consume(
                    (_reencode_vcf_line_worker(line) for line in data_lines()),
                    dst,
                )
        print(
            f"[graphreftovcf] finished realigning {completed} VCF rows",
            file=sys.stderr,
        )
        os.replace(stage_path, out_path)
    except Exception:
        try:
            os.remove(stage_path)
        except OSError:
            pass
        raise


_INTERNAL_SEQUENCE_INFO_IDS = {
    "CIGAR", "SVINSSEQ", "SVDELSEQ", "QUERYSEQ",
    "SEQ", "EXTENDGRAPHCIGAR",
}
_OPTIONAL_PA_FORMAT_IDS = {"ALLELENAME", "LABEL_H"}


def _plain_info_sequence(
    info: Mapping[str, str],
    query_lookup: Mapping[str, SeqRecord],
    cigar: str,
) -> str:
    """Return the representative query allele sequence for final INFO/SEQ."""
    existing = vcf_unescape(info.get("SEQ", ""))
    if existing:
        return existing

    coordinate = info.get("QUERYSEQ", "")
    if coordinate and coordinate != "." and "|" in coordinate:
        query_name, query_range = coordinate.split("|", 1)
        sequence, _record, _start, _end, _strand = (
            fetch_allele_query_sequence_with_record(
                ["", "", "", "", query_name, query_range],
                query_lookup,
            )
        )
        if sequence:
            return sequence

    legacy = vcf_unescape(info.get("SVINSSEQ", ""))
    if legacy:
        return legacy
    return extract_raw_inserted_sequence_from_cigar(cigar)


def _rewrite_final_sample_columns(
    columns: List[str], *, seqcompress: bool, pa_tag: bool,
) -> None:
    """Apply final payload and provenance policy to FORMAT/sample columns."""
    if len(columns) < 9:
        return
    format_names = columns[8].split(":")
    keep_indexes = [
        index for index, name in enumerate(format_names)
        if pa_tag or name not in _OPTIONAL_PA_FORMAT_IDS
    ]
    cigar_index = (
        format_names.index("EXTENDGRAPHCIGAR")
        if "EXTENDGRAPHCIGAR" in format_names else None
    )
    for column_index in range(9, len(columns)):
        values = columns[column_index].split(":")
        if len(values) != len(format_names):
            raise ValueError(
                "VCF sample column does not match FORMAT while applying "
                "the final sequence/provenance schema"
            )
        if cigar_index is not None and not seqcompress:
            cigar_values = values[cigar_index].split(",")
            values[cigar_index] = ",".join(
                vcf_escape(strip_cigar_payloads(vcf_unescape(value)))
                if value not in {"", "."} else "."
                for value in cigar_values
            )
        columns[column_index] = ":".join(values[index] for index in keep_indexes)
    columns[8] = ":".join(format_names[index] for index in keep_indexes)


def _final_sequence_header_lines(seqcompress: bool, pa_tag: bool) -> List[str]:
    payload_mode = "inline" if seqcompress else "absent"
    sequence_mode = "omitted" if seqcompress else "plain"
    cigar_description = (
        "Representative extended graph alignment CIGAR with inline query "
        "sequence payloads"
        if seqcompress else
        "Representative extended graph alignment CIGAR without inline "
        "sequence payloads; the plain query allele is in INFO/SEQ"
    )
    return [
        "##minsetrefSequenceEncoding=<"
        f"SEQ={sequence_mode},EXTENDGRAPHCIGARPayload={payload_mode}>\n",
        '##INFO=<ID=SEQ,Number=1,Type=String,Description="Representative '
        'plain query allele sequence; . when sequence compression is enabled '
        'or the event has no query sequence">\n',
        '##INFO=<ID=EXTENDGRAPHCIGAR,Number=1,Type=String,Description="'
        + cigar_description + '">\n',
        "##crossLocalRegionProvenanceTags="
        + (
            "retained"
            if pa_tag else
            "ALLELENAME_and_LABEL_H_removed"
        ) + "\n",
    ]


def write_final_vcf_schema(
    input_path: str,
    output_path: str,
    query_lookup: Optional[Dict[str, SeqRecord]] = None,
    *,
    seqcompress: bool = False,
    pa_tag: bool = True,
) -> None:
    """Write the public VCF schema from the richer internal staging VCF.

    Internal formatting retains QUERYSEQ and legacy sequence fields long
    enough for optional realignment.  The public file contains exactly the
    requested INFO/SEQ and INFO/EXTENDGRAPHCIGAR representation. ALLELENAME
    and its PA-relative LABEL_H coordinate are retained together by default
    and both are removed by --noPAtag.
    """
    query_lookup = query_lookup or {}
    global _QUERY_INTERVAL_INDEX, _QUERY_INTERVAL_INDEX_SOURCE
    _QUERY_INTERVAL_INDEX = build_contig_interval_index(query_lookup)
    _QUERY_INTERVAL_INDEX_SOURCE = query_lookup
    output_directory = os.path.dirname(os.path.abspath(output_path)) or "."
    descriptor, temporary = tempfile.mkstemp(
        prefix=os.path.basename(output_path) + ".",
        suffix=".schema.tmp",
        dir=output_directory,
    )
    os.close(descriptor)
    wrote_encoding_headers = False

    def write_encoding_headers(destination) -> None:
        nonlocal wrote_encoding_headers
        if wrote_encoding_headers:
            return
        destination.writelines(_final_sequence_header_lines(seqcompress, pa_tag))
        wrote_encoding_headers = True

    row_count = 0
    fast_rows = 0
    sv_rows = 0
    next_report = 500_000
    started = time.monotonic()
    print(
        f"[graphreftovcf] writing final VCF schema for {output_path} ...",
        file=sys.stderr,
    )
    # Rows whose FORMAT is exactly the full SNP schema only need the final two
    # PA-provenance fields removed; precompute that rewrite.
    snp_strippable_formats = {} if pa_tag else {
        SNP_FORMAT: (SNP_FORMAT_BASE, len(SNP_FORMAT_FIELDS) - 1, 2),
    }
    try:
        with open(input_path) as source, open(temporary, "w") as destination:
            for raw in source:
                if raw.startswith("##INFO=<ID="):
                    identifier = raw.split("##INFO=<ID=", 1)[1].split(",", 1)[0]
                    if identifier in _INTERNAL_SEQUENCE_INFO_IDS:
                        if not wrote_encoding_headers:
                            write_encoding_headers(destination)
                        continue
                if raw.startswith("##minsetrefSequenceEncoding=") or raw.startswith(
                    "##crossLocalRegionProvenanceTags="
                ):
                    continue
                if raw.startswith("##localBlockContributorSeparator=") and not pa_tag:
                    continue
                if raw.startswith("##FORMAT=<ID=EXTENDGRAPHCIGAR,"):
                    payload_text = (
                        "with inline query sequence payloads"
                        if seqcompress else "without inline sequence payloads"
                    )
                    destination.write(
                        '##FORMAT=<ID=EXTENDGRAPHCIGAR,Number=.,Type=String,'
        f'Description="Per-observation extended graph CIGAR '
        f'{payload_text}">\n'
                    )
                    continue
                if not pa_tag and (
                    raw.startswith("##FORMAT=<ID=ALLELENAME,")
                    or raw.startswith("##FORMAT=<ID=LABEL_H,")
                ):
                    continue
                if raw.startswith("#CHROM"):
                    write_encoding_headers(destination)
                    destination.write(raw)
                    continue
                if raw.startswith("#"):
                    destination.write(raw)
                    continue

                row_count += 1
                if row_count >= next_report:
                    elapsed = max(1e-6, time.monotonic() - started)
                    print(
                        f"[graphreftovcf] final schema: {row_count} rows "
                        f"({row_count / elapsed:.0f}/s)",
                        file=sys.stderr,
                    )
                    while next_report <= row_count:
                        next_report += 500_000
                head = raw.split("\t", 9)
                if len(head) < 9:
                    raise ValueError("malformed VCF record in final schema pass")
                format_text = head[8] if len(head) > 9 else head[8].rstrip("\r\n")
                if not is_sv_format(format_text):
                    stripped = snp_strippable_formats.get(format_text)
                    if stripped is not None:
                        base_format, format_colons, remove_count = stripped
                        columns = raw.rstrip("\r\n").split("\t")
                        if all(
                            columns[index].count(":") == format_colons
                            for index in range(9, len(columns))
                        ):
                            columns[8] = base_format
                            for index in range(9, len(columns)):
                                columns[index] = columns[index].rsplit(
                                    ":", remove_count,
                                )[0]
                            destination.write("\t".join(columns) + "\n")
                            fast_rows += 1
                            continue
                        # A field-count mismatch falls through to the generic
                        # path so it raises the same schema error as before.
                    else:
                        format_names = format_text.split(":")
                        if (
                            raw.endswith("\n")
                            and not raw.endswith("\r\n")
                            and (
                                pa_tag
                                or not _OPTIONAL_PA_FORMAT_IDS.intersection(
                                    format_names
                                )
                            )
                            and (
                                seqcompress
                                or "EXTENDGRAPHCIGAR" not in format_names
                            )
                        ):
                            destination.write(raw)
                            fast_rows += 1
                            continue
                columns = raw.rstrip("\r\n").split("\t")
                if len(columns) < 9:
                    raise ValueError("malformed VCF record in final schema pass")
                if is_sv_format(columns[8]):
                    sv_rows += 1
                    info_order = [
                        part.split("=", 1)[0]
                        for part in columns[7].split(";") if part
                    ]
                    info = parse_info_field(columns[7])
                    cigar = info.get("EXTENDGRAPHCIGAR", "") or info.get(
                        "CIGAR", ""
                    )
                    sequence = (
                        "" if seqcompress else
                        _plain_info_sequence(info, query_lookup, cigar)
                    )
                    svtype = str(info.get("SVTYPE", "")).upper()
                    if (
                        not seqcompress
                        and ("INS" in svtype or "SUB" in svtype)
                        and not sequence
                    ):
                        raise ValueError(
                            "cannot write an uncompressed insertion with "
                            "INFO/SEQ=.: "
                            f"{columns[0]}:{columns[1]}:{columns[2]} "
                            f"(QUERYSEQ={info.get('QUERYSEQ', '.')})"
                        )
                    final_cigar = (
                        cigar if seqcompress else strip_cigar_payloads(cigar)
                    )
                    final_info = {
                        key: value for key, value in info.items()
                        if key not in _INTERNAL_SEQUENCE_INFO_IDS
                    }
                    final_info["SEQ"] = vcf_escape(sequence) if sequence else "."
                    final_info["EXTENDGRAPHCIGAR"] = final_cigar or "."
                    final_order: List[str] = []
                    inserted = False
                    for key in info_order:
                        if key in _INTERNAL_SEQUENCE_INFO_IDS:
                            if not inserted:
                                final_order.extend(("SEQ", "EXTENDGRAPHCIGAR"))
                                inserted = True
                            continue
                        if key not in final_order:
                            final_order.append(key)
                    if not inserted:
                        final_order.extend(("SEQ", "EXTENDGRAPHCIGAR"))
                    columns[7] = format_info_field(final_info, final_order)
                _rewrite_final_sample_columns(
                    columns, seqcompress=seqcompress, pa_tag=pa_tag,
                )
                destination.write("\t".join(columns) + "\n")
        os.replace(temporary, output_path)
        elapsed = max(1e-6, time.monotonic() - started)
        print(
            f"[graphreftovcf] final schema complete for {output_path}: "
            f"{row_count} rows ({sv_rows} SV, {fast_rows} fast-path) in "
            f"{elapsed:.0f}s ({row_count / elapsed:.0f}/s)",
            file=sys.stderr,
        )
    except Exception:
        try:
            os.remove(temporary)
        except OSError:
            pass
        raise

# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def parse_columns_arg(columns: Optional[str]) -> Optional[List[str]]:
    if not columns:
        return None
    if os.path.exists(columns):
        vals: List[str] = []
        with open(columns) as fh:
            for raw in fh:
                raw = raw.strip()
                if raw:
                    vals.append(raw.split()[0])
        return vals
    return [x.strip() for x in columns.split(",") if x.strip()]


_ROW_WORKER_CONTEXT = {}
_FORMAT_WORKER_CONTEXT = {}


def _init_row_worker(ctx: dict) -> None:
    global _ROW_WORKER_CONTEXT
    global _ACTIVE_INTERVAL_OWNERSHIP
    _ROW_WORKER_CONTEXT = ctx
    _ACTIVE_INTERVAL_OWNERSHIP = ctx.get("interval_ownership", {})


def _row_worker_parse(item):
    line_no, raw = item
    row = parse_graphcigartoref_row(raw, line_no)
    if row is None:
        return line_no, None
    ctx = _ROW_WORKER_CONTEXT
    if not ctx:
        raise RuntimeError("row worker context is not initialized")
    results = []
    for ownership_row in expand_graphcigartoref_ownership_rows(row):
        results.append(parse_row_events(
            ownership_row,
            ref_lookup=ctx["ref_lookup"],
            query_lookup=ctx["query_records"],
            label_lookup=ctx["label_records"],
            sequence_lookup=ctx["sequence_records"],
            max_impute=ctx["max_impute"],
            bp_uncertain=ctx["breakpoint_uncertain"],
            gap_break=ctx["within_gap_break"],
            gap_penalty=ctx["within_gap_penalty"],
            var_in_insert_min_bp=ctx["var_in_insert_min_bp"],
            next_id_start=1,
            strict=ctx["strict"],
            edge_blackregion=ctx["edge_blackregion"],
            filter_non_acgt_sequences=ctx["filter_non_acgt_sequences"],
            score_main_alignment=ctx.get("score_main_alignment", False),
            separate_adjacent_indels=ctx.get(
                "separate_adjacent_indels", False,
            ),
            call_insertion_snps=ctx.get("call_insertion_snps", True),
        ))
    if len(results) == 1:
        return line_no, results[0]
    combined = RowParseResult([], [], [], [], [])
    all_non_acgt = bool(results)
    for result in results:
        combined.svs.extend(result.svs)
        combined.snps.extend(result.snps)
        combined.coverage_spans.extend(result.coverage_spans)
        combined.filtered_spans.extend(result.filtered_spans)
        combined.seen_queries.extend(result.seen_queries)
        combined.main_path_evidence.extend(result.main_path_evidence)
        combined.edge_filtered_events += result.edge_filtered_events
        combined.alternative_filtered_events += result.alternative_filtered_events
        all_non_acgt = (
            all_non_acgt
            and result.row_filter_reason == "non_acgt_query_sequence"
        )
    if all_non_acgt:
        combined.row_filter_reason = "non_acgt_query_sequence"
    return line_no, combined


def _init_format_worker(ctx: dict) -> None:
    global _FORMAT_WORKER_CONTEXT
    _FORMAT_WORKER_CONTEXT = ctx


def _format_snp_batch(batch):
    ctx = _FORMAT_WORKER_CONTEXT
    snps = ctx["snps"]
    rows: List[str] = []
    for key, start, end in batch:
        row = snp_cluster_to_row(
            key,
            snps[start:end],
            ctx["sample_names"],
            ctx["coverage_map"],
            ctx["filtered_map"],
        )
        if row:
            rows.append(row)
    return rows


def _format_sv_batch(batch):
    ctx = _FORMAT_WORKER_CONTEXT
    id_to_sv = ctx["id_to_sv"]
    cutoff = ctx["svcutoff"]
    rows: List[str] = []
    for cluster_id, member_ids in batch:
        members = [id_to_sv[value] for value in member_ids if value in id_to_sv]
        if not members:
            continue
        if cutoff is not None and max(sv.max_size for sv in members) <= cutoff:
            continue
        row = sv_cluster_to_row(
            cluster_id,
            members,
            ctx["sample_names"],
            ctx["coverage_map"],
            ctx["filtered_map"],
            ctx["ref_records"],
            ctx["query_lookup"],
            ctx["var_in_insert_min_bp"],
        )
        if row:
            rows.append(row)
    return rows


def _iter_snp_batch_specs(
    snps: Sequence[SNPUnit], batch_sites: int = 1024,
):
    batch = []
    index = 0
    while index < len(snps):
        snp = snps[index]
        key = (snp.chrom, snp.pos0, snp.ref_base, snp.alt_base)
        end = index + 1
        while end < len(snps):
            other = snps[end]
            other_key = (other.chrom, other.pos0, other.ref_base, other.alt_base)
            if other_key != key:
                break
            end += 1
        batch.append((key, index, end))
        if len(batch) >= batch_sites:
            yield batch
            batch = []
        index = end
    if batch:
        yield batch


def _iter_sv_batch_specs(
    cluster_ids: Sequence[Sequence[int]], batch_clusters: int = 256,
):
    batch = []
    for cluster_id, member_ids in enumerate(cluster_ids, 1):
        batch.append((cluster_id, tuple(member_ids)))
        if len(batch) >= batch_clusters:
            yield batch
            batch = []
    if batch:
        yield batch


def _format_rows_to_body(
    kind: str,
    specs,
    body_path: str,
    context: dict,
    processes: int,
) -> int:
    formatter = _format_snp_batch if kind == "SNP" else _format_sv_batch
    completed = 0
    next_report = 100_000
    global _FORMAT_WORKER_CONTEXT
    _FORMAT_WORKER_CONTEXT = context

    def consume(results, output):
        nonlocal completed, next_report
        for rows in results:
            for row in rows:
                output.write(row + "\n")
            completed += len(rows)
            if completed >= next_report:
                print(
                    f"[graphreftovcf] formatted {completed} {kind} VCF rows",
                    file=sys.stderr,
                )
                while next_report <= completed:
                    next_report += 100_000

    with open(body_path, "w") as output:
        if processes > 1:
            _freeze_shared_heap()
            # maxtasksperchild caps copy-on-write growth: SV batches index the
            # shared id_to_sv arena at random, and every touched object's
            # refcount write privatizes its page.  Recycling a worker after a
            # few batches resets its dirtied set, bounding per-worker RSS.
            with mp_context().Pool(
                processes=processes,
                initializer=_init_format_worker,
                initargs=(context,),
                maxtasksperchild=16,
            ) as pool:
                consume(pool.imap(formatter, specs, chunksize=1), output)
        else:
            consume((formatter(batch) for batch in specs), output)
    print(
        f"[graphreftovcf] finished formatting {completed} {kind} VCF rows",
        file=sys.stderr,
    )
    return completed


def _temporary_body_path(output_path: str, kind: str) -> str:
    directory = os.path.dirname(os.path.abspath(output_path)) or "."
    fd, path = tempfile.mkstemp(
        prefix=os.path.basename(output_path) + f".{kind.lower()}.",
        suffix=".body.tmp",
        dir=directory,
    )
    os.close(fd)
    return path


def _write_header_and_bodies(
    output_path: str,
    sample_names: Sequence[str],
    kind: str,
    body_paths: Sequence[str],
    reference_coverage_map=None,
    edge_blackregion: int = 0,
    filter_non_acgt_sequences: bool = True,
    reference_contigs: Sequence[Tuple[str, int]] = (),
    local_template_coordinates: bool = False,
    local_templates: Sequence[LocalTemplate] = (),
) -> None:
    started = time.monotonic()
    print(
        f"[graphreftovcf] assembling header and body into {output_path} ...",
        file=sys.stderr,
    )
    with open(output_path, "w") as output:
        write_vcf_header(
            output,
            sample_names,
            kind=kind,
            reference_coverage_map=reference_coverage_map,
            edge_blackregion=edge_blackregion,
            filter_non_acgt_sequences=filter_non_acgt_sequences,
            reference_contigs=reference_contigs,
            local_template_coordinates=local_template_coordinates,
            local_templates=local_templates,
        )
        for body_path in body_paths:
            with open(body_path) as body:
                shutil.copyfileobj(body, output, length=1024 * 1024)
    print(
        f"[graphreftovcf] assembled {output_path}: "
        f"{os.path.getsize(output_path)} bytes in "
        f"{time.monotonic() - started:.0f}s",
        file=sys.stderr,
    )


def main():
    parser = argparse.ArgumentParser(description="Rewrite graphcigartoref/refexalign output to VCF with safe within-line and across-sample clustering.")
    parser.add_argument("-i", "--input", required=True, help="graphcigartoref/refexalign output")
    parser.add_argument("--ref", "-r", required=True, help="reference FASTA with .fai")
    parser.add_argument(
        "--local-reference-templates",
        default=None,
        help="optional sparse local-template FASTA for reference-free graphs",
    )
    parser.add_argument("--query", "-q", default=None, help="legacy query metadata FASTA or liftover/coord TSV")
    parser.add_argument("-s", "--seq", "--fasta-query", dest="seq", default=None, help="FASTA containing actual query sequences; optional if --query is FASTA")
    parser.add_argument(
        "-m", "--coord-map", default=None,
        help=(
            "optional comma-ordered coordinate metadata TSVs; column 3 "
            "supplies label/query coordinates and trailing columns supply "
            "ownership metadata; later files replace matching PA names"
        ),
    )
    parser.add_argument("-o", "--output", required=True, help="output VCF, or prefix with --all")
    parser.add_argument("--all", action="store_true", help="write output.snp.vcf and output.indelsv.vcf")
    parser.add_argument("--svcutoff", "--svonly", nargs="?", const=20, default=None, type=int,
                        help="SV-only output mode: when provided, output only non-SNP SV/indel clusters with row MAXSIZE > cutoff; if provided without a value, use 20")
    parser.add_argument(
        "--seqcompress",
        action="store_true",
        help=(
            "store representative query sequence inline in "
            "INFO/EXTENDGRAPHCIGAR payloads and emit INFO/SEQ=.; by default "
            "INFO/SEQ is plain sequence and all EXTENDGRAPHCIGAR payloads "
            "are removed"
        ),
    )
    pa_group = parser.add_mutually_exclusive_group()
    pa_group.add_argument(
        "--PAtag", "--pa-tag",
        dest="pa_tag",
        action="store_true",
        default=True,
        help=(
            "retain PA provenance FORMAT tags ALLELENAME and LABEL_H "
            "(default)"
        ),
    )
    pa_group.add_argument(
        "--noPAtag", "--no-pa-tag",
        dest="pa_tag",
        action="store_false",
        help=(
            "remove final VCF ALLELENAME and PA-relative LABEL_H"
        ),
    )
    parser.add_argument("--lenient", action="store_true", help="keep old skip/filter behavior for unresolved rows or N-containing plain insertions instead of raising")
    parser.add_argument("--max-impute", type=int, default=10000, help="imputation size/cap used for extension ownership and breakpoint-uncertain neighbor-break calculation")
    parser.add_argument("--breakpoint-uncertain", type=int, default=200, help="breakpoint uncertainty for SV labels: self-side uses this bp size and keeps old label; neighbor-side uses 2x this size and emits *selflabel*neighbor")
    parser.add_argument("--within-gap-break", type=int, default=100, help="hard background bp break for within-line score merge")
    parser.add_argument("--within-gap-penalty", type=int, default=4, help="score penalty per background bp for within-line merge")
    parser.add_argument(
        "--separate-adjacent-indels",
        "--separate-adjacent-indel-components",
        dest="separate_adjacent_indels",
        action="store_true",
        help=(
            "emit directly adjacent insertion and deletion components as "
            "separate VCF rows instead of one interval replacement"
        ),
    )
    parser.add_argument(
        "--var-in-insert",
        nargs="?",
        const=100,
        default=0,
        type=int,
        metavar="BP",
        help=(
            "for merged insertions with more than BP uppercase/unmasked bases, "
            "align every other observation to the largest insertion with "
            "minimap2 and recursively encode supported residual insertions; "
            "if provided without a value, use 100 (default: disabled)"
        ),
    )
    parser.add_argument(
        "--edge-blackregion",
        type=int,
        default=0,
        metavar="BP",
        help=(
            "filter calls touching the first/last BP of the full reference "
            "graph alignment (per-sample v2); legacy rows retain their "
            "query-slice behavior (default: 0, disabled; set e.g. 500 to enable)"
        ),
    )
    parser.add_argument(
        "--allow-non-acgt-sequences",
        action="store_true",
        help=(
            "do not exclude complete query/allele rows containing bases other "
            "than A/C/G/T (default: exclude them)"
        ),
    )
    parser.add_argument("--cluster-block-gap", type=int, default=DEFAULT_BLOCK_GAP, help="split SV clustering blocks at reference gaps larger than this")
    parser.add_argument(
        "--validate-cross-graph-svs",
        "--merge-cross-row-deletions",
        dest="validate_cross_graph_svs",
        action="store_true",
        help=(
            "resolve all star-confirmed cross-graph SV types and assemble "
            "validated split deletions; enabled automatically by "
            "graphreftovcf_persample.py"
        ),
    )
    conflict_group = parser.add_mutually_exclusive_group()
    conflict_group.add_argument(
        "--resolve-conflicts-by-alignment-score",
        dest="resolve_conflicts_by_alignment_score",
        action="store_true",
        default=True,
        help=(
            "retain a singleton or inconsistent cross-graph SV only when "
            "its PA has the uniquely best main-path affine alignment score; "
            "enabled by default"
        ),
    )
    conflict_group.add_argument(
        "--no-resolve-conflicts-by-alignment-score",
        dest="resolve_conflicts_by_alignment_score",
        action="store_false",
        help="use the September 3 conflict-removal rule",
    )
    parser.add_argument(
        "--cross-graph-ratio-min", type=float, default=0.7,
        help="minimum complete-SV size ratio for cross-graph confirmation [0.7]",
    )
    parser.add_argument(
        "--cross-graph-ratio-max", type=float, default=1.0,
        help="maximum complete-SV size ratio for cross-graph confirmation [1.0]",
    )
    parser.add_argument(
        "--cross-graph-max-distance", type=int, default=500,
        help="absolute maximum assembly-start distance for confirmation [500]",
    )
    parser.add_argument(
        "--cross-graph-snp-window", type=int, default=10,
        help="same-assembly SNP duplicate distance window [10]",
    )
    parser.add_argument(
        "--inconsistency-report",
        default=None,
        metavar="TSV",
        help=(
            "write selected/rejected cross-graph conflict details; score "
            "resolution writes OUTPUT.inconsistencies.tsv by default"
        ),
    )
    parser.add_argument("-c", "--columns", help="comma-separated sample/haplotype names or file")
    parser.add_argument("--per-query-columns", action="store_true", help="use one VCF sample column per query record")
    realignment_group = parser.add_mutually_exclusive_group()
    realignment_group.add_argument("--realignment", dest="realignment", action="store_true", default=True, help="after writing the raw temp VCF, re-encode raw insertion sample CIGARs against the INFO target (default)")
    realignment_group.add_argument("--no-realignment", dest="realignment", action="store_false", help="disable final VCF insertion-CIGAR realignment")
    parser.add_argument("--force-realign", type=int, default=200, metavar="SIZE", help="if a sample CIGAR still has >= SIZE bp of unaligned sequence after encoding, retry encoding against the raw INFO insertion sequence and keep whichever has less unaligned; set to 0 to always retry (default: 200)")
    parser.add_argument("-t", "--processes", type=int, default=1, help="worker processes for input-row parsing, SV clustering, VCF formatting, and optional realignment")
    parser.add_argument("--format-processes", type=int, default=DEFAULT_FORMAT_PROCESSES, metavar="N", help=f"maximum workers for memory-heavy VCF formatting and optional realignment (default: {DEFAULT_FORMAT_PROCESSES}; parsing and clustering still use -t)")
    args = parser.parse_args()

    if args.processes < 1:
        parser.error("-t/--processes must be at least 1")
    if args.format_processes < 1:
        parser.error("--format-processes must be at least 1")
    if args.edge_blackregion < 0:
        parser.error("--edge-blackregion must be at least 0")
    if args.var_in_insert < 0:
        parser.error("--var-in-insert must be at least 0")
    if args.var_in_insert > 0 and mappy is None:
        parser.error(
            "--var-in-insert requires the mappy minimap2 Python bindings"
        )
    if not (
        0 < args.cross_graph_ratio_min
        <= args.cross_graph_ratio_max
        <= 1.0
    ):
        parser.error(
            "cross-graph ratios must satisfy 0 < MIN <= MAX <= 1"
        )
    if args.cross_graph_max_distance < 0:
        parser.error("--cross-graph-max-distance must be at least 0")
    if args.cross_graph_snp_window < 0:
        parser.error("--cross-graph-snp-window must be at least 0")
    if args.resolve_conflicts_by_alignment_score:
        args.validate_cross_graph_svs = True
    if args.inconsistency_report:
        args.validate_cross_graph_svs = True
    for attribute, option in (
        ("input", "--input"),
        ("ref", "--ref"),
        ("local_reference_templates", "--local-reference-templates"),
        ("query", "--query"),
        ("seq", "--fasta-query"),
        ("coord_map", "--coord-map"),
        ("output", "--output"),
        ("inconsistency_report", "--inconsistency-report"),
    ):
        try:
            setattr(
                args,
                attribute,
                normalize_cli_path(getattr(args, attribute), option),
            )
        except ValueError as error:
            parser.error(str(error))
    inconsistency_report_path = args.inconsistency_report
    if (
        inconsistency_report_path is None
        and args.resolve_conflicts_by_alignment_score
    ):
        inconsistency_report_path = default_inconsistency_report_path(
            args.output,
        )
    format_processes = min(args.processes, args.format_processes)
    if format_processes < args.processes:
        print(
            f"[graphreftovcf] using {format_processes} formatting worker(s); "
            f"-t {args.processes} still applies to parsing and clustering "
            f"(override with --format-processes)",
            file=sys.stderr,
        )

    ref_records, ref_order, _ref_reader = read_reference_records(args.ref)
    # Primary reference records are collected here. Local graph templates are
    # added separately by write_vcf_header so every possible CHROM receives
    # both its standard declaration and richer alternative-locus metadata.
    primary_reference_order = tuple(ref_order)
    local_templates = (
        add_local_template_reference_records(
            ref_records, ref_order, args.local_reference_templates,
        )
        if args.local_reference_templates else []
    )
    reference_contigs = [
        (name, ref_records[name].seqlen)
        for name in primary_reference_order
        if name in ref_records
    ]
    # Restrict the primary-reference portion of the ##contig header to records
    # named in the reference FASTA's own .fai. Local-template contigs are
    # appended separately by write_vcf_header; virtual per-block views remain
    # excluded (the reference-haplotype sample shares one FASTA between --ref
    # and --fasta-query, which used to leak those view names here).
    reference_fai_path = args.ref + ".fai"
    if os.path.isfile(reference_fai_path):
        fai_names = set()
        with open(reference_fai_path) as fai_handle:
            for fai_line in fai_handle:
                fields = fai_line.rstrip("\n").split("\t")
                if fields and fields[0]:
                    fai_names.add(fields[0])
        dropped_contigs = [
            name for name, _length in reference_contigs
            if name not in fai_names
        ]
        if dropped_contigs:
            print(
                f"[graphreftovcf] excluded {len(dropped_contigs)} "
                "non-reference record(s) from the ##contig header "
                f"(e.g. {dropped_contigs[:3]!r})",
                file=sys.stderr,
            )
        reference_contigs = [
            (name, length) for name, length in reference_contigs
            if name in fai_names
        ]
    if args.query:
        query_records, query_order, _query_reader = read_query_metadata(args.query, require_coord=False, max_impute=args.max_impute)
    else:
        query_records, query_order, _query_reader = {}, [], None
    sequence_records: Dict[str, SeqRecord] = {}
    if args.seq:
        sequence_records, _seq_order, _seq_reader = read_fasta_metadata(args.seq, require_coord=False, max_impute=args.max_impute, attach_reader=True)
    else:
        sequence_records = query_records
    coord_records: Dict[str, SeqRecord] = {}
    coord_map_paths = split_coord_map_paths(args.coord_map)
    if coord_map_paths:
        for coord_map_path in coord_map_paths:
            if not os.path.isfile(coord_map_path):
                parser.error(
                    f"--coord-map file does not exist: {coord_map_path}"
                )
        coord_records, coord_order = read_liftover_coord_maps(
            coord_map_paths, max_impute=args.max_impute,
        )
        query_records = merge_coord_metadata(query_records, coord_records)
        if coord_order:
            # Preserve original q order first, then metadata-only records.
            seen = set(query_order)
            for x in coord_order:
                if x not in seen:
                    query_order.append(x)
                    seen.add(x)

    # GenomeLift (optionally overlaid by genomeliftfix.tsv) is the sole PA
    # ownership authority.  Do not perform a second graph-CIGAR-dependent
    # election while formatting VCF records: doing so would make one-map and
    # two-map calls follow different ownership rules and could override the
    # decision already finalized by fill_graphcigartoref_gaps.py.
    interval_ownership = {}
    if coord_map_paths:
        print(
            "[graphreftovcf] loaded "
            f"{len(coord_map_paths)} comma-ordered coordinate map(s); later "
            "PA rows override earlier rows and coordinate-map ownership is "
            "authoritative",
            file=sys.stderr,
        )

    ref_lookup = dict(ref_records)
    # Reference-like graph paths can be resolved from query/sequence metadata too.
    ref_lookup.update(query_records)
    ref_lookup.update(sequence_records)
    query_sequence_lookup = dict(query_records)
    query_sequence_lookup.update(sequence_records)

    svs: List[SVUnit] = []
    snps: List[SNPUnit] = []
    coverage_spans_set = set()
    filtered_spans_set = set()
    seen_queries: List[str] = []
    seen_q_set = set()
    main_path_evidence: Dict[
        Tuple[str, str, str], MainPathAlignmentEvidence
    ] = {}
    next_id = 1
    non_acgt_filtered_rows = 0
    edge_filtered_events = 0

    alternative_filtered_events = 0

    # SNP observations are spilled to a temp pickle stream as parse results
    # arrive and reloaded only after the SV structures are freed, so the
    # parse-phase peak never holds SV and SNP objects together.  Sample
    # counts collected here keep every existing filter/print exact.
    spill_snps = os.environ.get("GRAPH_REFTOVCF_SPILL_SNPS", "1") not in {"0", ""}
    raw_columns = parse_columns_arg(args.columns)
    haplotype_mode = not args.per_query_columns
    snp_sample_counts: Dict[str, int] = {}
    spilled_snp_count = 0
    snp_spill_buffer: List[SNPUnit] = []
    snp_spill_handle = None
    snp_spill_path = None
    if spill_snps:
        _fd, snp_spill_path = tempfile.mkstemp(
            prefix=os.path.basename(args.output or "graphreftovcf") + ".",
            suffix=".snpspill.tmp",
            dir=os.path.dirname(os.path.abspath(args.output)) or "."
            if args.output else None,
        )
        snp_spill_handle = os.fdopen(_fd, "wb")

    def _flush_snp_spill():
        nonlocal snp_spill_buffer
        if snp_spill_buffer:
            pickle.dump(snp_spill_buffer, snp_spill_handle,
                        protocol=pickle.HIGHEST_PROTOCOL)
            snp_spill_buffer = []

    def admit_main_path_evidence(
        evidence: MainPathAlignmentEvidence,
    ) -> None:
        if not haplotype_mode:
            evidence = MainPathAlignmentEvidence(
                query=evidence.query,
                sample=evidence.query,
                allele_name=evidence.allele_name,
                qry_contig=evidence.qry_contig,
                main_match_bases=evidence.main_match_bases,
                complete_query_length=evidence.complete_query_length,
                line_no=evidence.line_no,
                main_mismatch_bases=evidence.main_mismatch_bases,
                main_gap_opens=evidence.main_gap_opens,
                main_gap_extensions=evidence.main_gap_extensions,
            )
        key = (
            evidence.sample,
            evidence.qry_contig,
            _canonical_cross_graph_label(evidence.allele_name),
        )
        previous = main_path_evidence.get(key)
        if (
            previous is None
            or _compare_main_path_evidence(evidence, previous) > 0
        ):
            main_path_evidence[key] = evidence

    def consume_row_result(line_no: int, result: Optional[RowParseResult]):
        nonlocal next_id, non_acgt_filtered_rows, edge_filtered_events
        nonlocal alternative_filtered_events
        nonlocal spilled_snp_count
        if result is None:
            return
        if result.row_filter_reason == "non_acgt_query_sequence":
            non_acgt_filtered_rows += 1
        edge_filtered_events += result.edge_filtered_events
        alternative_filtered_events += result.alternative_filtered_events
        for sv in result.svs:
            sv.id = next_id
            next_id += 1
        svs.extend(result.svs)
        if spill_snps:
            for snp in result.snps:
                if not haplotype_mode:
                    snp.sample = snp.query
                snp_sample_counts[snp.sample] = (
                    snp_sample_counts.get(snp.sample, 0) + 1
                )
            spilled_snp_count += len(result.snps)
            snp_spill_buffer.extend(result.snps)
            if len(snp_spill_buffer) >= 100_000:
                _flush_snp_spill()
        else:
            snps.extend(result.snps)
        coverage_spans_set.update(result.coverage_spans)
        filtered_spans_set.update(result.filtered_spans)
        for evidence in result.main_path_evidence:
            admit_main_path_evidence(evidence)
        for q in result.seen_queries:
            if q and q not in seen_q_set:
                seen_q_set.add(q)
                seen_queries.append(q)

    global _ROW_WORKER_CONTEXT
    global _ACTIVE_INTERVAL_OWNERSHIP
    _ACTIVE_INTERVAL_OWNERSHIP = interval_ownership
    _ROW_WORKER_CONTEXT = {
        "ref_lookup": ref_lookup,
        "query_records": query_records,
        "label_records": coord_records,
        "sequence_records": sequence_records,
        "max_impute": args.max_impute,
        "breakpoint_uncertain": args.breakpoint_uncertain,
        "within_gap_break": args.within_gap_break,
        "within_gap_penalty": args.within_gap_penalty,
        # SNPs come directly from the existing insertion-template alignments,
        # independently of cluster-level recursive SV encoding and its size
        # threshold. Keep the old row-local nested SV caller disabled.
        "var_in_insert_min_bp": 0,
        "call_insertion_snps": args.svcutoff is None,
        "edge_blackregion": args.edge_blackregion,
        "filter_non_acgt_sequences": not args.allow_non_acgt_sequences,
        "strict": not args.lenient,
        "interval_ownership": interval_ownership,
        "score_main_alignment": args.resolve_conflicts_by_alignment_score,
        "separate_adjacent_indels": args.separate_adjacent_indels,
    }

    #print(f"[graphreftovcf] parsing input rows (processes={args.processes}) ...", file=sys.stderr)
    with open_text(args.input) as fh:
        row_iter = ((line_no, raw) for line_no, raw in enumerate(fh, start=1))
        if args.processes > 1:
            _freeze_shared_heap()
            # Recycle parse workers periodically so refcount writes into the
            # shared record/sequence structures cannot accumulate into large
            # private copies over a worker's lifetime (64 chunks of 16 rows).
            with mp_context().Pool(processes=args.processes,
                                   initializer=_init_row_worker,
                                   initargs=(_ROW_WORKER_CONTEXT,),
                                   maxtasksperchild=64) as pool:
                for line_no, result in pool.imap(_row_worker_parse, row_iter, chunksize=16):
                    consume_row_result(line_no, result)
        else:
            for item in row_iter:
                line_no, result = _row_worker_parse(item)
                consume_row_result(line_no, result)
    _flush_snp_spill()
    parsed_snp_total = spilled_snp_count if spill_snps else len(snps)
    print(f"[graphreftovcf] parsed: {len(svs)} SVs, {parsed_snp_total} SNPs from {len(seen_queries)} queries", file=sys.stderr)

    if not haplotype_mode:
        # In per-query mode, sample columns are exact query records, not haplotype labels.
        for sv in svs:
            sv.sample = sv.query
        for snp in snps:
            snp.sample = snp.query
    if raw_columns is not None:
        sample_names = normalize_requested_haplotypes(raw_columns) if haplotype_mode else raw_columns
    elif args.per_query_columns:
        sample_names = query_order or seen_queries
    else:
        sample_names = []
        seen_samples = set()
        for q in query_order + seen_queries:
            sample = q if not haplotype_mode else get_haplotype(q)
            if sample and sample not in seen_samples:
                sample_names.append(sample)
                seen_samples.add(sample)
        # Some graphcigartoref rows carry a fuller overlap label than the FASTA
        # record name.  Include observed event samples so valid calls are not
        # silently dropped when --columns is absent.
        observed_snp_samples = (
            list(snp_sample_counts) if spill_snps
            else [x.sample for x in snps]
        )
        for sample in [x.sample for x in svs] + observed_snp_samples:
            if sample and sample not in seen_samples:
                sample_names.append(sample)
                seen_samples.add(sample)
    allowed = set(sample_names)

    print(
        f"[graphreftovcf] VCF sample columns ({len(sample_names)}): "
        f"{sample_names[:20]}{'...' if len(sample_names) > 20 else ''}",
        file=sys.stderr,
    )
    # Drop events outside requested columns early.  Spilled SNPs are filtered
    # at reload; their per-sample counts keep these checks and prints exact.
    observed_event_samples = sorted(
        {value.sample for value in svs}
        | (set(snp_sample_counts) if spill_snps else {value.sample for value in snps})
    )
    event_count_before_filter = len(svs) + parsed_snp_total
    if sample_names:
        svs = [sv for sv in svs if sv.sample in allowed]
        main_path_evidence = {
            key: evidence
            for key, evidence in main_path_evidence.items()
            if evidence.sample in allowed
        }
        if not spill_snps:
            snps = [s for s in snps if s.sample in allowed]
    filtered_snp_total = (
        (sum(count for sample, count in snp_sample_counts.items()
             if sample in allowed) if sample_names else spilled_snp_count)
        if spill_snps else len(snps)
    )
    if raw_columns is not None and event_count_before_filter and not (
        svs or filtered_snp_total
    ):
        raise ValueError(
            "--columns matched no parsed event samples; requested "
            f"{sample_names!r}, observed {observed_event_samples[:20]!r}"
        )
    print(f"[graphreftovcf] after sample filter: {len(svs)} SVs, {filtered_snp_total} SNPs", file=sys.stderr)
    if args.validate_cross_graph_svs:
        inconsistency_groups: List[CrossGraphInconsistency] = []
        validated_cross_graph_groups = 0
        assembled_cross_graph_deletions = 0
        dropped_cross_graph_events = 0
        if svs:
            (
                svs,
                validated_cross_graph_groups,
                assembled_cross_graph_deletions,
                dropped_cross_graph_events,
            ) = validate_cross_graph_svs(
                svs,
                ratio_min=args.cross_graph_ratio_min,
                ratio_max=args.cross_graph_ratio_max,
                max_distance=args.cross_graph_max_distance,
                query_lookup=query_sequence_lookup,
                inconsistency_groups=inconsistency_groups,
                main_path_evidence=tuple(main_path_evidence.values()),
                resolve_conflicts_by_alignment_score=(
                    args.resolve_conflicts_by_alignment_score
                ),
            )
            for event_id, sv in enumerate(svs, 1):
                sv.id = event_id
        if inconsistency_report_path:
            write_cross_graph_inconsistency_report(
                inconsistency_report_path, inconsistency_groups,
            )
        print(
            "[graphreftovcf] cross-graph SV validation: resolved "
            f"{validated_cross_graph_groups} group(s), assembled "
            f"{assembled_cross_graph_deletions} split deletion(s), removed "
            f"{dropped_cross_graph_events} inconsistent/singleton event(s)"
            + (
                f"; report: {inconsistency_report_path}"
                if inconsistency_report_path else ""
            ),
            file=sys.stderr,
        )
    if args.validate_cross_graph_svs and snps:
        # Distinct alternates at one reference base are retained as a valid
        # multi-allelic observation; the dedup pass counts them while it
        # sweeps instead of building a per-site dictionary here.
        (
            snps,
            merged_cross_graph_snps,
            multiallelic_sites,
        ) = deduplicate_cross_graph_snps(
            snps, assembly_window=args.cross_graph_snp_window,
        )
        print(
            "[graphreftovcf] cross-graph SNP validation: merged "
            f"{merged_cross_graph_snps} duplicate observation(s); retained "
            f"{multiallelic_sites} multi-allelic site(s)",
            file=sys.stderr,
        )
    deferred_edge_events = [
        sv for sv in svs
        if sv.edge_candidate and not sv.cross_row_merged
    ]
    if deferred_edge_events:
        for sv in deferred_edge_events:
            filtered_spans_set.add((
                sv.query,
                sv.ref_contig,
                sv.ref_start,
                max(sv.ref_end, sv.ref_start + 1),
            ))
        deferred_ids = {id(sv) for sv in deferred_edge_events}
        svs = [sv for sv in svs if id(sv) not in deferred_ids]
        edge_filtered_events += len(deferred_edge_events)
    print(
        "[graphreftovcf] sequence/edge filters: "
        f"{non_acgt_filtered_rows} query row(s) excluded for non-ACGT bases; "
        f"{edge_filtered_events} event(s) excluded within "
        f"{args.edge_blackregion} bp of a reference graph-alignment edge",
        file=sys.stderr,
    )
    coverage_map = build_sample_interval_map(coverage_spans_set, haplotype_mode=haplotype_mode, allowed_samples=allowed if sample_names else None)
    filtered_map = build_sample_interval_map(filtered_spans_set, haplotype_mode=haplotype_mode, allowed_samples=allowed if sample_names else None)
    # Header coverage must describe reference bases actually traversed by the
    # accepted main graph CIGAR (=, X, or D), not the coarse/assigned reference
    # coordinate carried in graphcigartoref column 4.
    reference_coverage_map = coverage_map
    coverage_spans_set.clear()
    filtered_spans_set.clear()

    print(f"[graphreftovcf] clustering {len(svs)} SVs ...", file=sys.stderr)
    from alternative_intervals import contains
    # Cross-row deletion assembly can enlarge an event after row parsing.
    svs = [sv for sv in svs if sv.alternative_regions is None or contains(
        sv.alternative_regions, sv.ref_contig, sv.ref_start, sv.ref_end)]
    print(f'[graphreftovcf] original alternative intervals: excluded {alternative_filtered_events} event(s)', file=sys.stderr)
    id_to_sv = {sv.id: sv for sv in svs}
    sv_cluster_ids = cluster_svs_parallel(svs, processes=max(1, args.processes), block_gap=args.cluster_block_gap) if svs else []
    # Do not let clustering bridge two separately projected calling windows.
    split_clusters = []
    for cluster in sv_cluster_ids:
        by_window = {}
        for event_id in cluster:
            event = id_to_sv[event_id]
            window = next((region for region in event.alternative_regions or ()
                           if contains((region,), event.ref_contig, event.ref_start, event.ref_end)), None)
            by_window.setdefault(window, []).append(event_id)
        split_clusters.extend(by_window.values())
    sv_cluster_ids = split_clusters
    print(
        f"[graphreftovcf] {len(sv_cluster_ids)} SV clusters, "
        f"{filtered_snp_total if spill_snps else len(snps)} SNP sites",
        file=sys.stderr,
    )
    def sv_cluster_sort_key(member_ids):
        first_contig = ""
        min_start = None
        min_id = None
        for value in member_ids:
            member = id_to_sv.get(value)
            if member is None:
                continue
            if min_id is None:
                first_contig = member.ref_contig
                min_start = member.ref_start
                min_id = member.id
            else:
                min_start = min(min_start, member.ref_start)
                min_id = min(min_id, member.id)
        return first_contig, min_start or 0, min_id or 0

    print("[graphreftovcf] ordering SV clusters ...", file=sys.stderr)
    sv_cluster_ids.sort(key=sv_cluster_sort_key)
    svcutoff_mode = args.svcutoff is not None

    def split_output_name(path: str, suffix: str) -> str:
        if path.endswith(".vcf"):
            return path[:-4] + suffix + ".vcf"
        return path + suffix + ".vcf"

    sv_body = _temporary_body_path(args.output, "sv")
    snp_body = None
    try:
        sv_context = {
            "id_to_sv": id_to_sv,
            "sample_names": tuple(sample_names),
            "coverage_map": coverage_map,
            "filtered_map": filtered_map,
            "ref_records": ref_records,
            "query_lookup": query_sequence_lookup,
            "var_in_insert_min_bp": max(0, int(args.var_in_insert or 0)),
            "svcutoff": args.svcutoff,
        }
        print(
            f"[graphreftovcf] formatting {len(sv_cluster_ids)} SV clusters "
            f"with {format_processes} process(es) ...",
            file=sys.stderr,
        )
        _format_rows_to_body(
            "SV",
            _iter_sv_batch_specs(sv_cluster_ids),
            sv_body,
            sv_context,
            format_processes,
        )
        # The SV body is self-contained on disk; free every SV structure
        # before the spilled SNP observations come back, so SV and SNP
        # objects never coexist at full size.
        sv_context.clear()
        sv_cluster_ids.clear()
        id_to_sv.clear()
        svs.clear()

        if not svcutoff_mode:
            if spill_snps:
                snp_spill_handle.close()
                reload_started = time.monotonic()
                allowed_samples = allowed if sample_names else None
                reloaded: List[SNPUnit] = []
                with open(snp_spill_path, "rb") as spill_in:
                    while True:
                        try:
                            chunk = pickle.load(spill_in)
                        except EOFError:
                            break
                        if allowed_samples is None:
                            reloaded.extend(chunk)
                        else:
                            reloaded.extend(
                                snp for snp in chunk
                                if snp.sample in allowed_samples
                            )
                snps = reloaded
                os.remove(snp_spill_path)
                snp_spill_path = None
                print(
                    f"[graphreftovcf] reloaded {len(snps)} spilled SNP "
                    f"observation(s) in "
                    f"{time.monotonic() - reload_started:.1f}s",
                    file=sys.stderr,
                )
                if args.validate_cross_graph_svs and snps:
                    (
                        snps,
                        merged_cross_graph_snps,
                        multiallelic_sites,
                    ) = deduplicate_cross_graph_snps(
                        snps, assembly_window=args.cross_graph_snp_window,
                    )
                    print(
                        "[graphreftovcf] cross-graph SNP validation: merged "
                        f"{merged_cross_graph_snps} duplicate observation(s);"
                        f" retained {multiallelic_sites} multi-allelic "
                        "site(s)",
                        file=sys.stderr,
                    )
            print("[graphreftovcf] ordering SNP observations ...", file=sys.stderr)
            snps.sort(key=lambda value: (
                value.chrom, value.pos0, value.ref_base, value.alt_base,
            ))
            snp_body = _temporary_body_path(args.output, "snp")
            snp_context = {
                "snps": snps,
                "sample_names": tuple(sample_names),
                "coverage_map": coverage_map,
                "filtered_map": filtered_map,
            }
            print(
                f"[graphreftovcf] formatting {len(snps)} SNP observations "
                f"with {format_processes} process(es) ...",
                file=sys.stderr,
            )
            _format_rows_to_body(
                "SNP",
                _iter_snp_batch_specs(snps),
                snp_body,
                snp_context,
                format_processes,
            )
            snp_context.clear()
            snps.clear()

        def finalize_vcf_temp(
            tmp_path: str, output_path: str, *, realign: bool,
        ) -> None:
            schema_input = tmp_path
            realigned_path = output_path + f".realigned.tmp.{os.getpid()}"
            if realign:
                line_level_reencode_vcf(
                    tmp_path,
                    realigned_path,
                    True,
                    ref_lookup,
                    query_sequence_lookup,
                    processes=format_processes,
                    force_realign=args.force_realign,
                )
                schema_input = realigned_path
                os.remove(tmp_path)
            try:
                write_final_vcf_schema(
                    schema_input,
                    output_path,
                    query_sequence_lookup,
                    seqcompress=args.seqcompress,
                    pa_tag=args.pa_tag,
                )
            finally:
                if os.path.exists(schema_input):
                    os.remove(schema_input)

        if svcutoff_mode:
            # --svcutoff is an SV-only output mode; SNP formatting is skipped.
            indelsv_out = (
                split_output_name(args.output, ".indelsv")
                if args.all else args.output
            )
            indelsv_tmp = indelsv_out + ".tmp"
            _write_header_and_bodies(
                indelsv_tmp, sample_names, "indelsv", [sv_body],
                reference_coverage_map,
                edge_blackregion=args.edge_blackregion,
                filter_non_acgt_sequences=not args.allow_non_acgt_sequences,
                reference_contigs=reference_contigs,
                local_template_coordinates=bool(local_templates),
                local_templates=local_templates,
            )
            finalize_vcf_temp(
                indelsv_tmp, indelsv_out, realign=args.realignment,
            )
        elif args.all:
            snp_out = split_output_name(args.output, ".snp")
            indelsv_out = split_output_name(args.output, ".indelsv")
            snp_tmp = snp_out + ".tmp"
            indelsv_tmp = indelsv_out + ".tmp"
            _write_header_and_bodies(
                snp_tmp, sample_names, "snp", [snp_body],
                reference_coverage_map,
                edge_blackregion=args.edge_blackregion,
                filter_non_acgt_sequences=not args.allow_non_acgt_sequences,
                reference_contigs=reference_contigs,
                local_template_coordinates=bool(local_templates),
                local_templates=local_templates,
            )
            _write_header_and_bodies(
                indelsv_tmp, sample_names, "indelsv", [sv_body],
                reference_coverage_map,
                edge_blackregion=args.edge_blackregion,
                filter_non_acgt_sequences=not args.allow_non_acgt_sequences,
                reference_contigs=reference_contigs,
                local_template_coordinates=bool(local_templates),
                local_templates=local_templates,
            )
            finalize_vcf_temp(snp_tmp, snp_out, realign=False)
            finalize_vcf_temp(
                indelsv_tmp, indelsv_out, realign=args.realignment,
            )
        else:
            tmp = args.output + ".tmp"
            _write_header_and_bodies(
                tmp, sample_names, "all", [sv_body, snp_body],
                reference_coverage_map,
                edge_blackregion=args.edge_blackregion,
                filter_non_acgt_sequences=not args.allow_non_acgt_sequences,
                reference_contigs=reference_contigs,
                local_template_coordinates=bool(local_templates),
                local_templates=local_templates,
            )
            finalize_vcf_temp(tmp, args.output, realign=args.realignment)
    finally:
        for body_path in (sv_body, snp_body):
            if body_path and os.path.exists(body_path):
                os.remove(body_path)
        if snp_spill_handle is not None:
            try:
                snp_spill_handle.close()
            except Exception:
                pass
        if snp_spill_path is not None and os.path.exists(snp_spill_path):
            os.remove(snp_spill_path)


if __name__ == "__main__":
    main()
