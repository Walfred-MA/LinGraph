#!/usr/bin/env python3
"""Graph CIGAR to reference liftover/encoding pipeline.

Production-clean version with debug/contribution tracing removed.
Behavior is intended to match graphcigartoref_lessbuggy_formatfix_realign_ssw_orientfix_maincont.py.
"""
from __future__ import annotations
import argparse
import bisect
import multiprocessing as mp
import os
import re
import shutil
import subprocess
import sys
import tempfile
import time
import traceback
from dataclasses import dataclass
from typing import Dict, List, Optional, Sequence, Set, Tuple
try:
    import numpy as _np
except ImportError:  # pure-python fallback keeps the pipeline runnable
    _np = None
try:
    import parasail as _parasail
except ImportError:  # numpy/pure-python DP fallback remains available
    _parasail = None
_PARASAIL_MATRIX_CACHE: dict = {}
_SV_REALIGNMENT_ENABLED = True
try:
    import mappy as _mappy
except ImportError:  # subprocess minimap2 fallback remains available
    _mappy = None
from local_reference_templates import read_templates
from graph_cigar_payloads import query_only_graph_cigar
from local_graph_whole_cigar import add_graph_cigar_payloads
from alternative_intervals import AlternativeIntervals, append_tag

_SEG_RE = re.compile('([<>])([^:<>]+):')
_COORD_RE = re.compile('^(.+):(\\d+)-(\\d+)([+-])$')
_CIGAR_RE = re.compile('(\\d+)([A-Z=])([A-Za-z]*)')


def set_sv_realignment_enabled(enabled: bool) -> None:
    """Enable or disable optional local and score-linked SV realignment."""
    global _SV_REALIGNMENT_ENABLED
    _SV_REALIGNMENT_ENABLED = bool(enabled)

def strip_match_name(name: str) -> str:
    if not name:
        return name
    if '[' in name and name.endswith(']'):
        return name.split('[', 1)[0]
    return name

def revcomp(seq: str) -> str:
    tbl = str.maketrans('ACGTacgtNn', 'TGCAtgcaNn')
    return seq.translate(tbl)[::-1]

class FastaRegionReader:

    def __init__(
        self,
        fasta_path: str,
        index_path: Optional[str]=None,
        template_fasta: Optional[str]=None,
    ):
        self.fasta_path = fasta_path
        candidates = []
        if index_path:
            candidates.append(index_path)
        candidates.extend([fasta_path + '.fai', fasta_path + '.faidx'])
        root, _ext = os.path.splitext(fasta_path)
        candidates.append(root + '.faidx')
        candidates = list(dict.fromkeys(candidates))
        self.index_path = next((x for x in candidates if x and os.path.exists(x)), None)
        if self.index_path is None:
            raise FileNotFoundError('Reference FASTA index not found: ' + ', '.join(candidates))
        self.index: Dict[str, Tuple[int, int, int, int]] = {}
        with open(self.index_path) as handle:
            for raw in handle:
                if not raw.strip():
                    continue
                cols = raw.rstrip('\n').split('\t')
                if len(cols) < 5:
                    continue
                name = cols[0]
                length, offset, line_bases, line_width = map(int, cols[1:5])
                self.index[name] = (length, offset, line_bases, line_width)
        self.handle = open(self.fasta_path, 'rb')
        self._primary_names = set(self.index)
        self._template_reader = None
        self._template_intervals = {}
        if template_fasta:
            self._template_reader = FastaRegionReader(template_fasta)
            for template in read_templates(template_fasta):
                if template.output_contig in self._primary_names:
                    raise ValueError(
                        'local-template output contig collides with requested '
                        f'reference FASTA: {template.output_contig!r}'
                    )
                self._template_intervals.setdefault(
                    template.output_contig, [],
                ).append(template)
                old_length = self.index.get(
                    template.output_contig, (0, 0, 1, 1),
                )[0]
                self.index[template.output_contig] = (
                    max(old_length, template.end), 0, 1, 1,
                )
            for rows in self._template_intervals.values():
                rows.sort(key=lambda row: (
                    row.start, row.end, row.graph_name, row.graph_path,
                ))

    def __contains__(self, name: str) -> bool:
        return name in self.index

    def fetch(self, name: str, start: int, end: int, strand: str='+') -> str:
        if name not in self.index:
            raise KeyError(f'{name!r} not found in FASTA index {self.index_path}')
        length, offset, line_bases, line_width = self.index[name]
        start = max(0, min(int(start), length))
        end = max(start, min(int(end), length))
        if end <= start:
            return ''
        if name not in self._primary_names:
            seq = self._fetch_template_interval(name, start, end)
            return revcomp(seq) if strand == '-' else seq
        start_line = start // line_bases
        end_line = end // line_bases
        start_byte = offset + start_line * line_width + start % line_bases
        end_byte_remainder = end % line_bases
        if end_byte_remainder == 0 and end_line > 0:
            end_byte = offset + (end_line - 1) * line_width + line_bases
        else:
            end_byte = offset + end_line * line_width + end_byte_remainder
        self.handle.seek(start_byte)
        raw_bytes = self.handle.read(end_byte - start_byte)
        seq = raw_bytes.replace(b'\n', b'').replace(b'\r', b'').decode('ascii')
        if strand == '-':
            seq = revcomp(seq)
        return seq

    def _fetch_template_interval(self, name: str, start: int, end: int) -> str:
        if self._template_reader is None:
            raise KeyError(f'no local-template reader configured for {name!r}')
        rows = self._template_intervals.get(name, ())
        cursor = start
        pieces = []
        while cursor < end:
            candidates = [
                row for row in rows
                if row.start <= cursor < row.end
            ]
            if not candidates:
                next_start = min(
                    (row.start for row in rows if row.start > cursor),
                    default=end,
                )
                gap_end = min(end, next_start)
                pieces.append('N' * (gap_end - cursor))
                cursor = gap_end
                continue
            row = min(candidates, key=lambda value: (
                value.end, value.graph_name, value.graph_path,
            ))
            piece_end = min(end, row.end)
            local_start = cursor - row.start
            local_end = piece_end - row.start
            if row.storage_strand == "+":
                piece = self._template_reader.fetch(
                    row.name, local_start, local_end,
                )
            else:
                length = row.end - row.start
                piece = self._template_reader.fetch(
                    row.name,
                    length - local_end,
                    length - local_start,
                    "-",
                )
            pieces.append(piece)
            cursor = piece_end
        return ''.join(pieces)

    def close(self) -> None:
        self.handle.close()
        if self._template_reader is not None:
            self._template_reader.close()


class InMemoryFastaRegionReader(FastaRegionReader):
    """FastaRegionReader serving primary contigs from RAM.

    Loading happens once in the parent; forked workers share the sequence
    pages copy-on-write, and primary fetches become in-memory slices with
    the same clamping and strand semantics as the file-backed reader.
    Template-backed virtual contigs keep the inherited on-disk path.
    """

    def __init__(
        self,
        fasta_path: str,
        index_path: Optional[str]=None,
        template_fasta: Optional[str]=None,
    ):
        super().__init__(
            fasta_path,
            index_path=index_path,
            template_fasta=template_fasta,
        )
        self._sequences: Dict[str, str] = {}
        for name in self._primary_names:
            length, offset, line_bases, line_width = self.index[name]
            if length <= 0:
                self._sequences[name] = ''
                continue
            full_lines = length // line_bases
            remainder = length % line_bases
            span = full_lines * line_width + remainder
            self.handle.seek(offset)
            raw = self.handle.read(span)
            sequence = raw.replace(b'\n', b'').replace(b'\r', b'').decode('ascii')
            if len(sequence) < length:
                # The final line may be shorter than line_width; read the tail.
                extra = self.handle.read(line_width)
                sequence += (
                    extra.replace(b'\n', b'').replace(b'\r', b'')
                    .decode('ascii')[: length - len(sequence)]
                )
            if len(sequence) != length:
                raise ValueError(
                    f'{fasta_path}: contig {name!r} has {len(sequence)} bases '
                    f'but the FAI reports {length}'
                )
            self._sequences[name] = sequence

    def fetch(self, name: str, start: int, end: int, strand: str='+') -> str:
        sequence = self._sequences.get(name)
        if sequence is None:
            return super().fetch(name, start, end, strand)
        length = self.index[name][0]
        start = max(0, min(int(start), length))
        end = max(start, min(int(end), length))
        if end <= start:
            return ''
        piece = sequence[start:end]
        return revcomp(piece) if strand == '-' else piece


def path_key(path: str) -> str:
    base = strip_match_name(path)
    if '~' in base:
        base = base.split('~', 1)[0]
    return base

def path_instance_name(path: str, segment_index: int, child_index: int=0) -> str:
    return f'{path}~{segment_index},{child_index}'

@dataclass(frozen=True)
class Coord:
    chrom: str
    start: int
    end: int
    strand: str


def extend_coord_for_cross_validation(coord: Coord, extension: int = 0) -> Coord:
    """Add symmetric context around a coordinate for cross-graph validation.

    This is deliberately separate from GenomeLift's ownership extension.  The
    latter decides which neighboring gap/overlap belongs to a block; this
    helper only requests additional comparison context.  Callers must still
    intersect the result with the complete graph alignment and its reliable
    interior before slicing a graph CIGAR.
    """
    extension = int(extension)
    if extension < 0:
        raise ValueError("cross-validation extension must be >= 0")
    if extension == 0:
        return coord
    return Coord(
        coord.chrom,
        max(0, coord.start - extension),
        coord.end + extension,
        coord.strand,
    )

@dataclass
class PairRow:
    line_no: int
    raw_line: str
    label: str
    ref_label: str
    query_coord_text: str
    ref_coord_text: str
    label_coord: str = ''
    query_name: Optional[str] = None
    query_coord: Optional[Coord] = None
    ref_name: Optional[str] = None
    ref_coord: Optional[Coord] = None

@dataclass(frozen=True)
class CigarOp:
    n: int
    op: str
    payload: str

@dataclass
class GraphicSegment:
    row_name: str
    seg_index: int
    direction: str
    path: str
    path_len: int
    start: int
    end: int
    cigar: str
    ops: List[CigarOp]

@dataclass
class EndpointRef:
    side: str
    segment_index: int
    is_end: int

@dataclass
class UnifiedBreakpoint:
    path: str
    group_index: int
    coord: int
    members: List[int]

@dataclass
class UnifiedSubsegment:
    side: str
    original_index: int
    path: str
    strand: str
    group_start: int
    group_end: int
    coord_start: int
    coord_end: int
    token_id: int = 0
    stage2_start: Optional[int] = None
    stage2_end: Optional[int] = None

    @property
    def span(self) -> int:
        return max(0, self.coord_end - self.coord_start)

@dataclass
class Stage1Result:
    ref_segments: List[GraphicSegment]
    qry_segments: List[GraphicSegment]
    unified_breaks: Dict[str, List[UnifiedBreakpoint]]
    ref_unified: List[UnifiedSubsegment]
    qry_unified: List[UnifiedSubsegment]
    ref_tokens: List[int]
    qry_tokens: List[int]
    lcs_pairs: List[Tuple[int, int]]
    ref_ordered: List[UnifiedSubsegment]
    qry_ordered: List[UnifiedSubsegment]
    reverse_ref_view: bool = False

@dataclass
class MatchedRun:
    path: str
    strand: str
    ref_original_index: int
    qry_original_index: int
    group_start: int
    group_end: int
    coord_start: int
    coord_end: int

@dataclass
class Stage2Result:
    matched_runs: List[MatchedRun]
    pieces: List['PairwisePiece']
    pairwise_cigar: str

@dataclass
class TemplateMergeSideInterval:
    r0: int
    r1: int
    op: str
    seq: Optional[str]
    q0: int
    q1: int

@dataclass
class TemplateMergeSideInsertion:
    seq: str
    q0: int
    q1: int
    has_payload: bool

@dataclass
class TemplateMergeSideTrack:
    intervals: List[TemplateMergeSideInterval]
    insertions: Dict[int, TemplateMergeSideInsertion]
    r_end: int
    q_end: int

@dataclass
class PairwisePiece:
    kind: str
    path: str
    ref_original_index: Optional[int]
    qry_original_index: Optional[int]
    coord_start: int
    coord_end: int
    ref_path_start: int
    ref_path_end: int
    cigar: str

@dataclass
class SegmentCutInfo:
    cuts: List[int]
    block_cigars: List[str]
    qposis: List[int]
    block_intervals: List[Tuple[int, int]]
    block_by_interval: Dict[Tuple[int, int], int]

@dataclass
class RawAssemblyUnit:
    ref_cigar: str
    qry_cigar: str
    template_seq: str
    ref_side_seq: Optional[str] = None
    qry_side_seq: Optional[str] = None
    break_side: Optional[str] = None
    break_len: int = 0

@dataclass
class OrderedRegion:
    side: str
    kind: str
    original_index: int
    path: str
    block_index: int
    cigar: str
    template_seq: str
    side_seq: Optional[str]
    seq_span: int
    anchor_index: int
    ordered_index: Optional[int]
    uses_template: bool = True
    coord_start: int = 0
    coord_end: int = 0

@dataclass
class Stage3PolishResult:
    pieces: List[PairwisePiece]
    pairwise_cigar: str

@dataclass
class LinearPiece:
    piece: PairwisePiece
    chrom: str
    strand: str
    genome_start: int
    genome_end: int
    cigar: str
    piece_index: int = -1
    piece_size: int = 0
    main_split_parent: int = -1

@dataclass
class Stage3Result:
    linear_pieces: List[LinearPiece]

@dataclass
class LocalDuplicateEncodedPiece:
    q0: int
    q1: int
    linear_piece: LinearPiece
    match_bases: int
_SECONDARY_LCS_MIN_MATCH_BASES = 200
_STAGE1_TERMINAL_LIFTOVER_MAX_SPAN = 200
_STAGE2_SIDE_ONLY_VALID_MIN_SPAN = 200

@dataclass
class QueryCoverage:
    segment_index: int
    q0: int
    q1: int
    path_q0: int
    path_q1: int
    gcigar: str
    mapping_gcigar: str
    path_seq_override: str = ''

@dataclass
class AnnealTemplate:
    q0: int
    q1: int
    chrom: str
    strand: str
    ref_start: int
    ref_end: int
    qry_gcigar: str
    ref_gcigar: str
    path_seq: str

def parse_coord(text: str) -> Optional[Coord]:
    if not text:
        return None
    m = _COORD_RE.match(text.strip())
    if not m:
        return None
    chrom, s, e, strand = m.groups()
    return Coord(chrom=chrom, start=int(s), end=int(e), strand=strand)

def coord_overlap(a: Coord, b: Coord) -> int:
    if a.chrom != b.chrom:
        return 0
    return max(0, min(a.end, b.end) - max(a.start, b.start))

def read_pairs(path: str) -> List[PairRow]:
    out: List[PairRow] = []
    with open(path) as fh:
        for i, line in enumerate(fh, 1):
            line = line.rstrip('\n')
            if not line:
                continue
            parts = line.split('\t')
            if len(parts) >= 5 and parse_coord(parts[1]) is not None and (parse_coord(parts[3]) is not None):
                label_coord = parts[5] if len(parts) >= 7 and parts[5] else parts[3]
                out.append(PairRow(line_no=i, raw_line=line, label=parts[4], label_coord=label_coord, ref_label=parts[2], query_coord_text=parts[1], ref_coord_text=parts[3], query_name=parts[0], ref_name=parts[2], query_coord=parse_coord(parts[1]), ref_coord=parse_coord(parts[3])))
                continue
            if len(parts) < 4:
                raise ValueError(f'{path}:{i}: expected >=4 columns')
            out.append(PairRow(line_no=i, raw_line=line, label=parts[0], label_coord=parts[2], ref_label=parts[1], query_coord_text=parts[2], ref_coord_text=parts[3]))
    return out

def read_fasta(path: str) -> Dict[str, str]:
    seqs: Dict[str, str] = {}
    name: Optional[str] = None
    chunks: List[str] = []
    with open(path) as fh:
        for line in fh:
            line = line.rstrip('\n')
            if not line:
                continue
            if line.startswith('>'):
                if name is not None:
                    seq = ''.join(chunks)
                    seqs[name] = seq
                    seqs[strip_match_name(name)] = seq
                name = line[1:].split()[0]
                chunks = []
            else:
                chunks.append(line.strip())
    if name is not None:
        seq = ''.join(chunks)
        seqs[name] = seq
        seqs[strip_match_name(name)] = seq
    return seqs

def read_fasta_header_aliases(path: str, alias_column: int=2) -> Tuple[Dict[str, str], List[Tuple[str, Coord, str]], Dict[str, Coord]]:
    aliases: Dict[str, str] = {}
    coord_index: List[Tuple[str, Coord, str]] = []
    coord_by_name: Dict[str, Coord] = {}
    with open(path) as fh:
        for line in fh:
            if not line.startswith('>'):
                continue
            cols = line[1:].rstrip('\n').split()
            if not cols:
                continue
            rec_name = cols[0]
            aliases[rec_name] = rec_name
            aliases[strip_match_name(rec_name)] = rec_name
            if alias_column > 0 and len(cols) >= alias_column:
                alias = cols[alias_column - 1]
                aliases[alias] = rec_name
                aliases[strip_match_name(alias)] = rec_name
                coord = parse_coord(alias)
                if coord is not None:
                    coord_index.append((rec_name, coord, alias))
                    coord_by_name[rec_name] = coord
    return (aliases, coord_index, coord_by_name)

def read_graphic_cigars(path: str) -> Dict[str, str]:
    out: Dict[str, str] = {}
    with open(path) as fh:
        for i, line in enumerate(fh, 1):
            line = line.rstrip('\n')
            if not line:
                continue
            parts = line.split('\t')
            if len(parts) < 3:
                raise ValueError(f'{path}:{i}: expected >=3 columns')
            out[parts[0]] = parts[2]
    return out

def read_graph_sequences(path: str) -> Dict[str, str]:
    return read_fasta(path)

def read_graph_path_coords(path: str, coord_column: int=2) -> Dict[str, Coord]:
    out: Dict[str, Coord] = {}
    with open(path) as fh:
        for line in fh:
            if not line.startswith('>'):
                continue
            cols = line[1:].rstrip('\n').split()
            if len(cols) < coord_column:
                continue
            name = cols[0]
            coord = parse_coord(cols[coord_column - 1])
            if coord is None:
                continue
            out[name] = coord
            out[path_key(name)] = coord
    return out

def read_fai_lengths(path: str) -> Dict[str, int]:
    fai = path + '.fai'
    out: Dict[str, int] = {}
    with open(fai) as fh:
        for line in fh:
            cols = line.rstrip('\n').split('\t')
            if len(cols) >= 2:
                out[cols[0]] = int(cols[1])
    return out

def synthesize_reference_mapping_from_coord(coord_text: str, ref_reader: Optional[FastaRegionReader]) -> str:
    coord = parse_coord(strip_match_name(coord_text))
    if coord is None:
        return ''
    chrom = coord.chrom
    start = coord.start
    end = coord.end
    span = end - start
    if span <= 0:
        return ''
    strand = coord.strand
    if ref_reader is not None and chrom in ref_reader.index:
        chrom_len = ref_reader.index[chrom][0]
        if strand == '-':
            return f'<{chrom}:{max(0, chrom_len - end)}H{span}={max(0, start)}H'
        return f'>{chrom}:{max(0, start)}H{span}={max(0, chrom_len - end)}H'
    if strand == '-':
        return f'<{chrom}:0H{span}={start}H'
    return f'>{chrom}:{start}H{span}='

def read_graph_header_mappings(graph_fasta: str, ref_reader: Optional[FastaRegionReader], column: int=3, coord_column: int=2) -> Dict[str, str]:
    out: Dict[str, str] = {}
    with open(graph_fasta) as fh:
        for raw in fh:
            if not raw.startswith('>'):
                continue
            fields = raw[1:].strip().split()
            if not fields:
                continue
            name = fields[0]
            chosen = ''
            if len(fields) >= column and fields[column - 1].startswith(('>', '<')):
                chosen = fields[column - 1]
            elif len(fields) >= coord_column:
                chosen = synthesize_reference_mapping_from_coord(fields[coord_column - 1], ref_reader)
            if chosen:
                out[name] = chosen
                out[strip_match_name(name)] = chosen
    return out

def slice_graphic_mapping_by_query(mapping_gcigar: str, qstart: int, qend: int, ref_reader: Optional[FastaRegionReader]) -> str:
    qstart = int(qstart)
    qend = int(qend)
    if qend <= qstart:
        return ''
    out: List[str] = []
    qpos_global = 0
    for seg in parse_graphic_segments(mapping_gcigar, 'mapping'):
        body = _strip_terminal_h(seg.ops)
        qpos = qpos_global
        rpos = seg.start if seg.direction != '<' else seg.end
        ref_start: Optional[int] = None
        ref_end: Optional[int] = None
        sub: List[CigarOp] = []

        def mark_ref(a: int, b: int) -> None:
            nonlocal ref_start, ref_end
            if b < a:
                a, b = b, a
            if b <= a:
                return
            ref_start = a if ref_start is None else min(ref_start, a)
            ref_end = b if ref_end is None else max(ref_end, b)

        for tok in body:
            q_consume = query_consume(tok.op, tok.n)
            r_consume = ref_consume_pair(tok.op, tok.n)
            if q_consume > 0:
                ov0 = max(qstart, qpos)
                ov1 = min(qend, qpos + q_consume)
                if ov1 > ov0:
                    take = ov1 - ov0
                    off = ov0 - qpos
                    if tok.op in {'=', 'X'}:
                        if seg.direction == '<':
                            local_r1 = rpos - off
                            local_r0 = local_r1 - take
                        else:
                            local_r0 = rpos + off
                            local_r1 = local_r0 + take
                    else:
                        local_r0 = rpos
                        local_r1 = rpos
                    mark_ref(local_r0, local_r1)
                    payload = ''
                    if tok.payload and len(tok.payload) == tok.n:
                        payload = tok.payload[off:off + take]
                    sub.append(CigarOp(n=take, op=tok.op, payload=payload))
            elif tok.op == 'D' and r_consume > 0:
                if qstart <= qpos < qend:
                    if seg.direction == '<':
                        del_r0 = rpos - r_consume
                        del_r1 = rpos
                    else:
                        del_r0 = rpos
                        del_r1 = rpos + r_consume
                    mark_ref(del_r0, del_r1)
                    sub.append(tok)
            qpos += q_consume
            if seg.direction == '<':
                rpos -= r_consume
            else:
                rpos += r_consume
        qpos_global = qpos
        if not sub:
            continue
        if ref_start is None:
            ref_start = seg.start
        if ref_end is None:
            ref_end = ref_start
        out.append(_format_interval_gcigar(seg.direction, seg.path, seg.path_len, ref_start, ref_end, sub))
    return ''.join(out)

def invert_path_to_ref_slice_for_compare(sliced_gcigar: str, path_name: str, path_start: int, path_end: int, path_len: int, path_direction: str, ref_reader: Optional[FastaRegionReader]) -> Tuple[str, Optional[Tuple[str, str, int, int, Optional[int]]]]:
    """Invert a graph-path-to-reference slice for local comparison.

                ``sliced_gcigar`` is a mapping where the graph path is the query and the
                linear reference chromosome is the path.  For comparing the sample query
                against the lifted reference sequence on the graph path, we need the inverse
                representation: the linear reference sequence as the query and the graph
                path as the path.

                CIGAR operation conversion is therefore:
                        graph->ref "="  -> ref->graph "="
                        graph->ref "X"  -> ref->graph "X" with the *reference* base payload
                        graph->ref "I"  -> ref->graph "D" (graph has extra bases)
                        graph->ref "D"  -> ref->graph "I" with the reference inserted bases

                Coordinates returned in ref_meta are always forward chromosome coordinates
                [ref_start, ref_end), independent of reference strand.
                """
    segs = parse_graphic_segments(sliced_gcigar, 'mapping_slice')
    if len(segs) != 1:
        return ('', None)
    seg = segs[0]
    ref_direction = seg.direction
    ref_chrom = seg.path
    body = _strip_terminal_h(seg.ops)
    chrom_len = ref_reader.index[ref_chrom][0] if ref_reader is not None and ref_chrom in ref_reader.index else None
    rpos = seg.start if ref_direction == '>' else seg.end
    ref_min: Optional[int] = None
    ref_max: Optional[int] = None

    def mark_ref(a: int, b: int) -> None:
        nonlocal ref_min, ref_max
        if b < a:
            a, b = (b, a)
        if b <= a:
            return
        ref_min = a if ref_min is None else min(ref_min, a)
        ref_max = b if ref_max is None else max(ref_max, b)

    def fetch_ref_interval(a: int, b: int, strand: str) -> str:
        if b < a:
            a, b = (b, a)
        n = max(0, b - a)
        if n == 0:
            return ''
        if ref_reader is None or ref_chrom not in ref_reader:
            return 'N' * n
        seq = ref_reader.fetch(ref_chrom, a, b, strand)
        if len(seq) < n:
            seq += 'N' * (n - len(seq))
        return seq[:n]

    def consume_reference_bases(n: int) -> str:
        nonlocal rpos
        if n <= 0:
            return ''
        if ref_direction == '>':
            a, b = (rpos, rpos + n)
            rpos = b
            mark_ref(a, b)
            return fetch_ref_interval(a, b, '+')
        a, b = (rpos - n, rpos)
        rpos = a
        mark_ref(a, b)
        return fetch_ref_interval(a, b, '-')
    inv_body: List[CigarOp] = []
    for tok in body:
        if tok.op == '=':
            consume_reference_bases(tok.n)
            inv_body.append(CigarOp(n=tok.n, op='=', payload=''))
        elif tok.op == 'X':
            consume_reference_bases(tok.n)
            inv_body.append(CigarOp(n=tok.n, op='X', payload='N' * tok.n))
        elif tok.op == 'D':
            inv_body.append(CigarOp(n=tok.n, op='I', payload=consume_reference_bases(tok.n)))
        elif tok.op == 'I':
            inv_body.append(CigarOp(n=tok.n, op='D', payload=''))
        else:
            raise ValueError(f'unsupported op {tok.op!r} in graph-to-reference mapping')
    if not inv_body:
        return ('', None)
    if ref_min is None or ref_max is None:
        ref_min = ref_max = rpos
    direction = path_direction
    body_ops = inv_body
    if path_direction == '<':
        flipped: List[CigarOp] = []
        for tok in reversed(inv_body):
            payload = tok.payload
            if payload and tok.op in {'I', 'X'}:
                payload = revcomp(payload)
            flipped.append(CigarOp(tok.n, tok.op, payload))
        body_ops = flipped
    gcigar = _format_interval_gcigar(direction, path_name, path_len, path_start, path_end, body_ops)
    return (gcigar, (ref_direction, ref_chrom, ref_min, ref_max, chrom_len))

def resolve_pairs_against_q_records(pairs: Sequence[PairRow], gcigars: Dict[str, str], header_aliases: Dict[str, str], coord_index: Sequence[Tuple[str, Coord, str]]) -> List[PairRow]:
    out: List[PairRow] = []
    for pair in pairs:
        qcoord = parse_coord(pair.query_coord_text)
        rcoord = parse_coord(pair.ref_coord_text)
        pair.query_coord = qcoord
        pair.ref_coord = rcoord
        qname = None
        for cand in (pair.label, strip_match_name(pair.label)):
            if cand in gcigars:
                qname = cand
                break
        if qname is None:
            qname = header_aliases.get(pair.label) or header_aliases.get(strip_match_name(pair.label))
        if qname is None and qcoord is not None:
            best = max(coord_index, key=lambda row: coord_overlap(qcoord, row[1]), default=None)
            if best is not None and coord_overlap(qcoord, best[1]) > 0:
                qname = best[0]
        rname = None
        for cand in (pair.ref_label, strip_match_name(pair.ref_label)):
            if cand in gcigars:
                rname = cand
                break
        if rname is None:
            rname = header_aliases.get(pair.ref_label) or header_aliases.get(strip_match_name(pair.ref_label))
        if rname is None and rcoord is not None:
            best = max(coord_index, key=lambda row: coord_overlap(rcoord, row[1]), default=None)
            if best is not None and coord_overlap(rcoord, best[1]) > 0:
                rname = best[0]
        if qname not in gcigars or rname not in gcigars:
            continue
        pair.query_name = qname
        pair.ref_name = rname
        out.append(pair)
    return out

def parse_cigar_ops(text: str) -> List[CigarOp]:
    ops: List[CigarOp] = []
    for m in _CIGAR_RE.finditer(text):
        op = m.group(2)
        if op == 'M':
            op = '='
        ops.append(CigarOp(n=int(m.group(1)), op=op, payload=m.group(3) or ''))
    return ops

def ref_consume(op: str, n: int) -> int:
    return n if op in {'=', 'X', 'D', 'H'} else 0

def parse_graphic_segments(gcigar: str, row_name: str) -> List[GraphicSegment]:
    out: List[GraphicSegment] = []
    matches = list(_SEG_RE.finditer(gcigar))
    for seg_index, m in enumerate(matches):
        direction = m.group(1)
        path = m.group(2)
        body_start = m.end()
        body_end = matches[seg_index + 1].start() if seg_index + 1 < len(matches) else len(gcigar)
        cigar = gcigar[body_start:body_end]
        ops = parse_cigar_ops(cigar)
        left_h = ops[0].n if ops and ops[0].op == 'H' else 0
        right_h = ops[-1].n if ops and ops[-1].op == 'H' else 0
        body = list(ops)
        if body and body[0].op == 'H':
            body = body[1:]
        if body and body[-1].op == 'H':
            body = body[:-1]
        span = sum((ref_consume(tok.op, tok.n) for tok in body))
        path_len = left_h + span + right_h
        if direction == '>':
            start = left_h
            end = path_len - right_h
        else:
            start = right_h
            end = path_len - left_h
        out.append(GraphicSegment(row_name=row_name, seg_index=seg_index, direction=direction, path=path, path_len=path_len, start=start, end=end, cigar=cigar, ops=ops))
    return out

def graphic_segment_to_gcigar(segment: GraphicSegment) -> str:
    return f"{segment.direction}{segment.path}:{''.join((_op_to_cigar(tok) for tok in segment.ops))}"

def query_consume(op: str, n: int) -> int:
    return n if op in {'=', 'X', 'I'} else 0

def ref_consume_pair(op: str, n: int) -> int:
    return n if op in {'=', 'X', 'D', 'H'} else 0

def _format_interval_gcigar(direction: str, path: str, path_len: int, start: int, end: int, body_ops: Sequence[CigarOp]) -> str:
    span = max(0, end - start)
    if direction == '<':
        left_h = max(0, path_len - end)
        right_h = max(0, start)
    else:
        left_h = max(0, start)
        right_h = max(0, path_len - end)
    ops: List[str] = []
    if left_h > 0:
        ops.append(f'{left_h}H')
    if span > 0 or body_ops:
        for tok in body_ops:
            ops.append(_op_to_cigar(tok))
    if right_h > 0:
        ops.append(f'{right_h}H')
    return f"{direction}{path}:{''.join(ops)}"

def _x_payload_kind(payload: str, n: int) -> str:
    """Return mismatch payload encoding: query-only (n) or block pairwise (2n).

                The pairwise block format is ref_bases + query_bases.  For example,
                4XNNNNGAGG means ref=NNNN and query=GAGG.
                """
    if not payload:
        return 'empty'
    if len(payload) == n:
        return 'query'
    if len(payload) == 2 * n:
        return 'block'
    raise ValueError(f'X payload length mismatch for {n}X: payload length={len(payload)}')

def _slice_x_payload(payload: str, full_n: int, offset: int, take: int) -> str:
    """Slice an X payload while preserving its encoding."""
    if not payload:
        return ''
    kind = _x_payload_kind(payload, full_n)
    if kind == 'query':
        return payload[offset:offset + take]
    if kind == 'block':
        ref_payload = payload[:full_n]
        qry_payload = payload[full_n:]
        return ref_payload[offset:offset + take] + qry_payload[offset:offset + take]
    return ''

def _x_payload_query(payload: str, n: int) -> str:
    """Return the query-side bases from an X payload."""
    if not payload:
        return 'N' * n
    kind = _x_payload_kind(payload, n)
    if kind == 'query':
        return payload
    if kind == 'block':
        return payload[n:]
    return 'N' * n

def _merge_payload_for_same_op(op: str, old_n: int, old_payload: str, n: int, payload: str) -> str:
    """Merge payloads when adjacent CIGAR ops of the same type are coalesced."""
    if op == 'X':
        old_kind = _x_payload_kind(old_payload, old_n) if old_payload else 'empty'
        new_kind = _x_payload_kind(payload, n) if payload else 'empty'
        if old_kind == 'block' and new_kind == 'block':
            old_ref, old_qry = (old_payload[:old_n], old_payload[old_n:])
            new_ref, new_qry = (payload[:n], payload[n:])
            return old_ref + new_ref + old_qry + new_qry
        if old_kind == 'query' and new_kind == 'query':
            return old_payload + payload
        if old_kind == 'empty':
            return payload
        if new_kind == 'empty':
            return old_payload
        old_ref = old_payload[:old_n] if old_kind == 'block' else 'N' * old_n
        old_qry = old_payload[old_n:] if old_kind == 'block' else old_payload
        new_ref = payload[:n] if new_kind == 'block' else 'N' * n
        new_qry = payload[n:] if new_kind == 'block' else payload
        return old_ref + new_ref + old_qry + new_qry
    return old_payload + payload

def _reverse_mapping_payload(tok: CigarOp) -> str:
    payload = tok.payload
    if not payload:
        return payload
    if tok.op == 'I':
        return revcomp(payload)
    if tok.op == 'X':
        if len(payload) == tok.n:
            return revcomp(payload)
        if len(payload) == 2 * tok.n:
            return _revcomp_x_payload(payload, tok.n)
    return payload

def _reverse_graphic_segment_for_reverse_path(seg: GraphicSegment) -> str:
    body = _strip_terminal_h(seg.ops)
    rev_body: List[CigarOp] = []
    for tok in reversed(body):
        rev_body.append(CigarOp(tok.n, tok.op, _reverse_mapping_payload(tok)))
    new_direction = '<' if seg.direction == '>' else '>'
    return _format_interval_gcigar(new_direction, seg.path, seg.path_len, seg.start, seg.end, rev_body)

def reverse_mapping_gcigar_for_reverse_path(mapping_gcigar: str) -> str:
    segs = parse_graphic_segments(mapping_gcigar, 'mapping_reverse')
    return ''.join((_reverse_graphic_segment_for_reverse_path(seg) for seg in reversed(segs)))
_PAYLOAD_RE = re.compile('(\\d+)([=XDIHM])([A-Za-z]*)')

def parse_cigar_with_query_payload(cigar: str) -> List[Tuple[int, str, str]]:
    ops: List[Tuple[int, str, str]] = []
    for m in _PAYLOAD_RE.finditer(cigar):
        size = int(m.group(1))
        op = m.group(2)
        payload = m.group(3) or ''
        if payload and op not in {'I', 'X'}:
            raise ValueError(f'Unexpected payload on op {op}: {size}{op}{payload}')
        if op == 'I' and payload and (len(payload) != size):
            raise ValueError(f'Payload length mismatch for {size}{op}')
        ops.append((size, op, payload))
    return ops

def format_cigar_ops(ops: Sequence[Tuple[int, str, str]]) -> str:
    out: List[str] = []
    for size, op, payload in ops:
        if size <= 0:
            continue
        if payload:
            out.append(f'{size}{op}{payload}')
        else:
            out.append(f'{size}{op}')
    return ''.join(out)

def split_cigar_by_rranges(cigar: str, rranges: Sequence[Tuple[int, int]]) -> List[Tuple[str, Optional[int], Optional[int]]]:
    if not rranges:
        return []
    for i in range(1, len(rranges)):
        if rranges[i][0] != rranges[i - 1][1]:
            raise ValueError('rranges must be consecutive')
    segments = parse_cigar_with_query_payload(cigar)
    cutted: List[List[Tuple[int, str, str]]] = [[] for _ in rranges]
    q_starts: List[Optional[int]] = [None for _ in rranges]
    q_ends: List[Optional[int]] = [None for _ in rranges]
    rpos = 0
    qpos = 0
    last_i = len(rranges) - 1

    def add_piece(i: int, size: int, op: str, payload: str='', q_s: Optional[int]=None, q_e: Optional[int]=None) -> None:
        if size <= 0:
            return
        if payload:
            if op == 'X':
                _x_payload_kind(payload, size)
            elif len(payload) != size:
                raise ValueError(f'Payload length mismatch while adding {size}{op}: payload length={len(payload)}')
        if cutted[i] and cutted[i][-1][1] == op:
            old_size, old_op, old_payload = cutted[i][-1]
            if op in {'I', 'X'}:
                merged = _merge_payload_for_same_op(op, old_size, old_payload, size, payload)
                cutted[i][-1] = (old_size + size, old_op, merged)
            else:
                cutted[i][-1] = (old_size + size, old_op, '')
        else:
            cutted[i].append((size, op, payload if op in {'I', 'X'} else ''))
        if q_s is not None:
            if q_starts[i] is None:
                q_starts[i] = q_s
            q_ends[i] = q_e
    for size, op, payload in segments:
        if op in {'=', 'M', 'X', 'D', 'H'}:
            r0 = rpos
            r1 = rpos + size
            for i, (s, e) in enumerate(rranges):
                ov_s = max(s, r0)
                ov_e = min(e, r1)
                if ov_e <= ov_s:
                    continue
                take = ov_e - ov_s
                offset = ov_s - r0
                if op in {'D', 'H'}:
                    add_piece(i, take, op, payload='', q_s=qpos, q_e=qpos)
                elif op == 'X':
                    sub_payload = _slice_x_payload(payload, size, offset, take) if payload else ''
                    add_piece(i, take, 'X', payload=sub_payload, q_s=qpos + offset, q_e=qpos + offset + take)
                else:
                    add_piece(i, take, op, payload='', q_s=qpos + offset, q_e=qpos + offset + take)
            rpos += size
            if op not in {'D', 'H'}:
                qpos += size
            continue
        if op == 'I':
            assigned = False
            for i in range(1, len(rranges)):
                prev_s, prev_e = rranges[i - 1]
                curr_s, _curr_e = rranges[i]
                if prev_e == curr_s == rpos:
                    add_piece(i - 1, size, 'I', payload=payload, q_s=qpos, q_e=qpos + size)
                    assigned = True
                    break
            if not assigned:
                for i, (s, e) in enumerate(rranges):
                    is_last = i == last_i
                    if s <= rpos < e or (is_last and rpos == e):
                        add_piece(i, size, 'I', payload=payload, q_s=qpos, q_e=qpos + size)
                        assigned = True
                        break
            qpos += size
    results: List[Tuple[str, Optional[int], Optional[int]]] = []
    for cut, qs, qe in zip(cutted, q_starts, q_ends):
        if qs is None:
            results.append(('', None, None))
        else:
            results.append((format_cigar_ops(cut), qs, qe))
    return results

def split_cigar_by_rposi(cigar: str, rposis: List[int], reference_seq: Optional[str]=None, reference_start: int=0) -> Tuple[List[str], List[int]]:
    if any((rposis[i] >= rposis[i + 1] for i in range(len(rposis) - 1))):
        raise ValueError('rposis must be strictly increasing')
    cuts = [int(x) for x in rposis]
    n_blocks = len(cuts) + 1
    out: List[List[Tuple[int, str, str]]] = [[] for _ in range(n_blocks)]
    cut_to_q: Dict[int, Optional[int]] = {cut: None for cut in cuts}

    def ref_slice(r0: int, r1: int) -> str:
        n = max(0, r1 - r0)
        if n == 0:
            return ''
        if reference_seq is None:
            return 'N' * n
        s = r0 - reference_start
        e = r1 - reference_start
        left_pad = max(0, -s)
        right_pad = max(0, e - len(reference_seq))
        s = max(0, s)
        e = min(len(reference_seq), e)
        seq = 'N' * left_pad + reference_seq[s:e] + 'N' * right_pad
        if len(seq) < n:
            seq += 'N' * (n - len(seq))
        return seq[:n]

    def parse_cigar_with_payload(text: str) -> List[Tuple[int, str, str]]:
        ops = []
        pos = 0
        pat = re.compile('(\\d+)([=MXIDH])([A-Za-z]*)')
        while pos < len(text):
            m = pat.match(text, pos)
            if not m:
                raise ValueError(f'Cannot parse CIGAR near: {text[pos:pos + 80]}')
            n = int(m.group(1))
            op = m.group(2)
            payload = m.group(3) or ''
            if op == '=':
                op = 'M'
            if payload and op not in {'I', 'D', 'X'}:
                raise ValueError(f'Unexpected payload on {n}{op}: {payload}')
            if op in {'I', 'D'} and payload and (len(payload) != n):
                raise ValueError(f'Payload length mismatch for {n}{op}: payload length={len(payload)}')
            if op == 'X' and payload and (len(payload) not in {n, 2 * n}):
                raise ValueError(f'X payload length mismatch for {n}X: payload length={len(payload)}')
            ops.append((n, op, payload))
            pos = m.end()
        return ops

    def add_piece(block_i: int, n: int, op: str, payload: str='') -> None:
        if n <= 0:
            return
        if op == 'M':
            payload = ''
        if op in {'I', 'D', 'X'} and (not payload):
            payload = 'N' * n
        if payload and op == 'X':
            _x_payload_kind(payload, n)
        if out[block_i] and out[block_i][-1][1] == op:
            old_n, old_op, old_payload = out[block_i][-1]
            merged = _merge_payload_for_same_op(op, old_n, old_payload, n, payload)
            out[block_i][-1] = (old_n + n, old_op, merged)
        else:
            out[block_i].append((n, op, payload))

    def format_ops(ops: List[Tuple[int, str, str]]) -> str:
        parts: List[str] = []
        for n, op, payload in ops:
            if n <= 0:
                continue
            emit_op = '=' if op == 'M' else op
            if payload:
                parts.append(f'{n}{emit_op}{payload}')
            else:
                parts.append(f'{n}{emit_op}')
        return ''.join(parts)

    def block_for_insertion_anchor(rpos: int) -> int:
        """
                                I at rpos is before reference base rpos.
                                Assign it to the block containing rpos - 1.
                                """
        if not cuts:
            return 0
        anchor = rpos - 1
        if anchor < cuts[0]:
            return 0
        for i in range(1, len(cuts)):
            if cuts[i - 1] <= anchor < cuts[i]:
                return i
        return len(cuts)

    def x_query_payload(payload: str, full_n: int, offset: int, take: int) -> str:
        if not payload:
            return 'N' * take
        if len(payload) == 2 * full_n:
            return payload[full_n + offset:full_n + offset + take]
        return payload[offset:offset + take]

    def normal_payload(payload: str, full_n: int, op: str, offset: int, take: int, r_abs_start: int) -> str:
        if op == 'M':
            return ''
        if payload:
            if op == 'X' and len(payload) == 2 * full_n:
                return _slice_x_payload(payload, full_n, offset, take)
            if op == 'X' and len(payload) == full_n:
                qry_payload = payload[offset:offset + take]
                ref_payload = ref_slice(r_abs_start, r_abs_start + take)
                return ref_payload + qry_payload
            return payload[offset:offset + take]
        if op == 'D':
            return ref_slice(r_abs_start, r_abs_start + take)
        if op in {'I', 'X'}:
            return 'N' * take
        return ''

    def record_qpos_for_ref_interval(cut: int, r0: int, q0: int, size: int, op: str) -> None:
        """
                                Record qpos immediately before reference base `cut`.

                                M/X:
                                                q advances with r.

                                D/H:
                                                q does not advance.

                                I at the same cut can overwrite this later because I at rpos is
                                before reference base rpos.
                                """
        if op in {'M', 'X'}:
            if r0 <= cut <= r0 + size:
                if cut_to_q[cut] is None:
                    cut_to_q[cut] = q0 + (cut - r0)
        elif op in {'D', 'H'}:
            if r0 <= cut <= r0 + size:
                if cut_to_q[cut] is None:
                    cut_to_q[cut] = q0
    ops = parse_cigar_with_payload(cigar)
    rpos = 0
    qpos = 0
    neg_inf = -10 ** 30
    pos_inf = 10 ** 30
    boundaries = [neg_inf] + cuts + [pos_inf]
    for size, op, payload in ops:
        for cut in cuts:
            record_qpos_for_ref_interval(cut, rpos, qpos, size, op)
        if op == 'H':
            rpos += size
            continue
        if op == 'I':
            block_i = block_for_insertion_anchor(rpos)
            add_piece(block_i, size, 'I', payload if payload else 'N' * size)
            if rpos in cut_to_q:
                cut_to_q[rpos] = qpos + size
            qpos += size
            continue
        if op not in {'M', 'X', 'D'}:
            raise ValueError(f'Unsupported op: {op}')
        r0 = rpos
        r1 = rpos + size
        for block_i in range(n_blocks):
            b0 = boundaries[block_i]
            b1 = boundaries[block_i + 1]
            ov_s = max(r0, b0)
            ov_e = min(r1, b1)
            if ov_e <= ov_s:
                continue
            take = ov_e - ov_s
            offset = ov_s - r0
            is_outer_block = block_i == 0 or block_i == n_blocks - 1
            if is_outer_block:
                if op in {'M', 'X'}:
                    if op == 'X':
                        i_payload = x_query_payload(payload, size, offset, take)
                    else:
                        i_payload = ref_slice(ov_s, ov_e)
                    add_piece(block_i, take, 'I', i_payload)
                elif op == 'D':
                    pass
            else:
                sub_payload = normal_payload(payload, size, op, offset, take, ov_s)
                add_piece(block_i, take, op, sub_payload)
        rpos = r1
        if op in {'M', 'X'}:
            qpos += size
    for cut in cuts:
        if cut_to_q[cut] is None:
            cut_to_q[cut] = qpos
    cutted_cigars = [format_ops(block_ops) for block_ops in out]
    qposis = [int(cut_to_q[cut]) for cut in cuts]
    qposis.append(qpos)
    return (cutted_cigars, qposis)

def _mapping_piece_query_ranges(mapping_gcigar: str) -> List[Tuple[GraphicSegment, Tuple[int, int]]]:
    out: List[Tuple[GraphicSegment, Tuple[int, int]]] = []
    qpos = 0
    for seg in parse_graphic_segments(mapping_gcigar, 'mapping_pieces'):
        qlen = 0
        for tok in _strip_terminal_h(seg.ops):
            qlen += query_consume(tok.op, tok.n)
        out.append((seg, (qpos, qpos + qlen)))
        qpos += qlen
    return out

def _aligned_size_for_piece(cigar: str) -> int:
    return sum((int(tok[:-1]) for tok in re.findall('\\d+[=XD]', cigar)))

def build_query_coverages(qry_segments: Sequence[GraphicSegment], graph_sequences: Optional[Dict[str, str]]=None, graph_mappings: Optional[Dict[str, str]]=None) -> List[QueryCoverage]:
    out: List[QueryCoverage] = []
    qpos = 0
    for seg in qry_segments:
        qlen = 0
        for tok in _strip_terminal_h(seg.ops):
            qlen += query_consume(tok.op, tok.n)
        seg_q0 = qpos
        seg_q1 = qpos + qlen
        mapping = None
        if graph_mappings is not None:
            mapping = graph_mappings.get(seg.path) or graph_mappings.get(strip_match_name(seg.path))
        if mapping:
            mapping_for_split = mapping
            mapping_parts = _mapping_piece_query_ranges(mapping_for_split)
            if len(mapping_parts) > 1:
                split_seg = seg
                if seg.direction == '<':
                    split_seg = parse_graphic_segments(_reverse_graphic_segment_for_reverse_path(seg), 'qry_cov_split_forward')[0]
                clipped_parts: List[Tuple[GraphicSegment, int, int, int, int]] = []
                for mseg, (rr0, rr1) in mapping_parts:
                    clip0 = max(split_seg.start, rr0)
                    clip1 = min(split_seg.end, rr1)
                    if clip1 <= clip0:
                        continue
                    clipped_parts.append((mseg, rr0, rr1, clip0, clip1))
                if not clipped_parts:
                    qpos = seg_q1
                    continue
                child_pieces = split_cigar_by_rranges(split_seg.cigar, [(clip0, clip1) for _mseg, _rr0, _rr1, clip0, clip1 in clipped_parts])
                for part_i, ((mseg, rr0, _rr1, clip0, clip1), (cut_cigar, local_q0, local_q1)) in enumerate(zip(clipped_parts, child_pieces), start=1):
                    if not cut_cigar or local_q0 is None or local_q1 is None:
                        continue
                    if _aligned_size_for_piece(cut_cigar) < 200:
                        continue
                    map_local0 = clip0 - rr0
                    map_local1 = clip1 - rr0
                    child_mapping_gcigar = slice_graphic_mapping_by_query(graphic_segment_to_gcigar(mseg), map_local0, map_local1, None)
                    if not child_mapping_gcigar:
                        continue
                    child_path = path_instance_name(split_seg.path, seg.seg_index, part_i)
                    child_gcigar = f'{split_seg.direction}{child_path}:{cut_cigar}'
                    if seg.direction == '<':
                        # child_mapping_gcigar must stay in the forward graph-path coordinate
                        # frame.  It is later sliced by qry_slice_seg.start/end, which are
                        # forward coordinates on the synthetic child path; invert_path_to_ref_slice_for_compare()
                        # applies qry_slice_seg.direction to orient the comparison.
                        child_gcigar = reverse_mapping_gcigar_for_reverse_path(child_gcigar)
                        q_lo = qlen - local_q1
                        q_hi = qlen - local_q0
                    else:
                        q_lo = local_q0
                        q_hi = local_q1
                    child_path_seq = ''
                    if graph_sequences is not None:
                        parent_seq = graph_sequences.get(split_seg.path) or graph_sequences.get(path_key(split_seg.path))
                        if parent_seq is not None:
                            if clip0 < 0 or clip1 > len(parent_seq):
                                raise ValueError(f'{seg.row_name} segment {seg.seg_index}: child path slice {clip0}-{clip1} out of bounds for {split_seg.path}')
                            child_path_seq = parent_seq[clip0:clip1]
                    out.append(QueryCoverage(segment_index=seg.seg_index, q0=seg_q0 + q_lo, q1=seg_q0 + q_hi, path_q0=0, path_q1=clip1 - clip0, gcigar=child_gcigar, mapping_gcigar=child_mapping_gcigar, path_seq_override=child_path_seq))
                qpos = seg_q1
                continue
        base_gcigar = graphic_segment_to_gcigar(seg)
        out.append(QueryCoverage(segment_index=seg.seg_index, q0=seg_q0, q1=seg_q1, path_q0=0, path_q1=seg.path_len, gcigar=base_gcigar, mapping_gcigar=mapping or '', path_seq_override=''))
        qpos += qlen
    return out

def slice_query_segment_by_query(
    segment: GraphicSegment,
    rel_qstart: int,
    rel_qend: int,
    *,
    include_left_boundary_deletions: bool = True,
    include_right_boundary_deletions: bool = False,
) -> str:
    rel_qstart = int(rel_qstart)
    rel_qend = int(rel_qend)
    if rel_qend <= rel_qstart:
        return ''
    body = _strip_terminal_h(segment.ops)
    qpos = 0
    rpos = segment.start if segment.direction != '<' else segment.end
    ref_min: Optional[int] = None
    ref_max: Optional[int] = None
    anchor_pos: Optional[int] = None
    sub: List[CigarOp] = []

    def mark_ref(a: int, b: int) -> None:
        nonlocal ref_min, ref_max
        if b < a:
            a, b = (b, a)
        if b <= a:
            return
        ref_min = a if ref_min is None else min(ref_min, a)
        ref_max = b if ref_max is None else max(ref_max, b)
    for tok in body:
        q_consume = query_consume(tok.op, tok.n)
        r_consume = ref_consume_pair(tok.op, tok.n)
        if q_consume > 0:
            ov0 = max(rel_qstart, qpos)
            ov1 = min(rel_qend, qpos + q_consume)
            if ov1 > ov0:
                take = ov1 - ov0
                off = ov0 - qpos
                if tok.op in {'=', 'X'}:
                    if segment.direction == '<':
                        local_r1 = rpos - off
                        local_r0 = local_r1 - take
                    else:
                        local_r0 = rpos + off
                        local_r1 = local_r0 + take
                    mark_ref(local_r0, local_r1)
                    anchor_pos = local_r0
                elif tok.op == 'I':
                    local_r0 = rpos
                    local_r1 = rpos
                    if anchor_pos is None:
                        anchor_pos = rpos
                else:
                    local_r0 = rpos
                    local_r1 = rpos
                    if anchor_pos is None:
                        anchor_pos = rpos
                payload = ''
                if tok.op == 'I' and tok.payload and (len(tok.payload) == tok.n):
                    payload = tok.payload[off:off + take]
                elif tok.op == 'X' and tok.payload:
                    payload = _slice_x_payload(tok.payload, tok.n, off, take)
                sub.append(CigarOp(n=take, op=tok.op, payload=payload))
        elif tok.op == 'D' and r_consume > 0:
            if (
                rel_qstart < qpos < rel_qend
                or (
                    include_left_boundary_deletions
                    and rel_qstart == qpos < rel_qend
                )
                or (
                    include_right_boundary_deletions
                    and rel_qstart < qpos == rel_qend
                )
            ):
                if segment.direction == '<':
                    del_r0 = rpos - r_consume
                    del_r1 = rpos
                else:
                    del_r0 = rpos
                    del_r1 = rpos + r_consume
                mark_ref(del_r0, del_r1)
                anchor_pos = del_r0
                sub.append(CigarOp(tok.n, tok.op, tok.payload))
        qpos += q_consume
        if segment.direction == '<':
            rpos -= r_consume
        else:
            rpos += r_consume
    if not sub:
        return ''
    if ref_min is None or ref_max is None:
        anchor = anchor_pos if anchor_pos is not None else segment.start
        ref_min = anchor
        ref_max = anchor
    return _format_interval_gcigar(segment.direction, segment.path, segment.path_len, ref_min, ref_max, sub)


def graph_cigar_query_span(graph_cigar: str, row_name: str='graph-cigar') -> int:
    """Return the number of query bases represented by a graph CIGAR."""
    return sum(
        query_consume(operation.op, operation.n)
        for segment in parse_graphic_segments(graph_cigar, row_name)
        for operation in segment.ops
    )


def slice_graph_cigar_by_query(
    graph_cigar: str,
    qstart: int,
    qend: int,
    row_name: str='graph-cigar',
    *,
    include_left_boundary_deletions: bool = True,
    include_right_boundary_deletions: bool = False,
) -> str:
    """Slice a complete graph CIGAR in its query-coordinate frame.

    This belongs in the core because both an input graph alignment and a final
    reference CIGAR must be clipped with exactly the same payload-preserving
    operation.  In particular, X/I payload bases are retained by
    :func:`slice_query_segment_by_query`.
    """
    qstart = int(qstart)
    qend = int(qend)
    if qstart < 0 or qend <= qstart:
        raise ValueError(f'invalid graph-CIGAR query slice {qstart}-{qend}')
    pieces: List[str] = []
    cursor = 0
    first_output_piece = True
    for segment in parse_graphic_segments(graph_cigar, row_name):
        segment_span = sum(
            query_consume(operation.op, operation.n)
            for operation in segment.ops
        )
        overlap_start = max(qstart, cursor)
        overlap_end = min(qend, cursor + segment_span)
        if overlap_end > overlap_start:
            piece = slice_query_segment_by_query(
                segment,
                overlap_start - cursor,
                overlap_end - cursor,
                include_left_boundary_deletions=(
                    include_left_boundary_deletions or not first_output_piece
                ),
                include_right_boundary_deletions=(
                    include_right_boundary_deletions
                    and overlap_end == qend
                ),
            )
            if piece:
                pieces.append(piece)
                first_output_piece = False
        cursor += segment_span
    if qend > cursor:
        raise ValueError(
            f'graph-CIGAR query slice {qstart}-{qend} exceeds query span {cursor}'
        )
    if not pieces:
        raise ValueError(
            f'graph-CIGAR query slice {qstart}-{qend} contains no usable path'
        )
    sliced = ''.join(pieces)
    observed = graph_cigar_query_span(sliced, row_name)
    expected = qend - qstart
    if observed != expected:
        raise ValueError(
            f'sliced graph-CIGAR query span is {observed}, expected {expected}'
        )
    return sliced


def reference_coord_from_graph_cigar(
    graph_cigar: str,
    reference_chrom: str,
    fallback: Optional[Coord]=None,
    row_name: str='graph-cigar',
) -> Coord:
    """Recover the reference interval actually retained in a sliced CIGAR."""
    wanted = path_key(reference_chrom)
    matches = [
        segment
        for segment in parse_graphic_segments(graph_cigar, row_name)
        if path_key(segment.path) == wanted
    ]
    if not matches:
        if fallback is None:
            raise ValueError(
                f'no {reference_chrom!r} segment remains in sliced graph CIGAR'
            )
        return fallback
    start = min(segment.start for segment in matches)
    end = max(segment.end for segment in matches)
    directions = {segment.direction for segment in matches}
    if len(directions) == 1:
        strand = '-' if next(iter(directions)) == '<' else '+'
    elif fallback is not None:
        strand = fallback.strand
    else:
        strand = '+'
    return Coord(reference_chrom, start, end, strand)

def _segment_forward_for_template(segment: GraphicSegment) -> GraphicSegment:
    if segment.direction == '>':
        return segment
    return parse_graphic_segments(_reverse_graphic_segment_for_reverse_path(segment), f'{segment.row_name}_forward')[0]

def _segment_forward_to_local_coord(segment: GraphicSegment, coord: int) -> int:
    """Convert forward path coordinate to this segment's CIGAR/template coordinate."""
    coord = int(coord)
    if segment.direction == '<':
        return segment.path_len - coord
    return coord

def _segment_forward_interval_to_local(segment: GraphicSegment, start: int, end: int) -> Tuple[int, int]:
    """Convert a forward path interval to local CIGAR-order template coordinates."""
    start = int(start)
    end = int(end)
    if end < start:
        start, end = (end, start)
    if segment.direction == '<':
        return (segment.path_len - end, segment.path_len - start)
    return (start, end)

def _segment_local_interval_to_forward(segment: GraphicSegment, start: int, end: int) -> Tuple[int, int]:
    """Convert a local CIGAR-order template interval back to forward path coordinates."""
    start = int(start)
    end = int(end)
    if end < start:
        start, end = (end, start)
    if segment.direction == '<':
        return (segment.path_len - end, segment.path_len - start)
    return (start, end)

def _oriented_path_seq(path_seq: str, start: int, end: int, direction: str) -> str:
    """Return path bases for forward interval [start,end) in the requested local direction."""
    start = max(0, int(start))
    end = max(start, min(int(end), len(path_seq)))
    seq = path_seq[start:end]
    return revcomp(seq) if direction == '<' else seq

def _segment_local_template_seq(segment: GraphicSegment, local_start: int, local_end: int, path_sequences: Dict[str, str]) -> str:
    path_seq = path_sequences.get(segment.path) or path_sequences.get(path_key(segment.path))
    if path_seq is None:
        raise KeyError(f'missing graph/path sequence for {segment.path}')
    f0, f1 = _segment_local_interval_to_forward(segment, local_start, local_end)
    return _oriented_path_seq(path_seq, f0, f1, segment.direction)

def _segment_forward_template_seq(segment: GraphicSegment, start: int, end: int, path_sequences: Dict[str, str]) -> str:
    path_seq = path_sequences.get(segment.path) or path_sequences.get(path_key(segment.path))
    if path_seq is None:
        raise KeyError(f'missing graph/path sequence for {segment.path}')
    return _oriented_path_seq(path_seq, start, end, segment.direction)

def _reverse_template_cigar(cigar: str) -> str:
    """Reverse-complement an ordinary same-template CIGAR body into opposite template order."""
    ops: List[CigarOp] = []
    for tok in reversed(parse_cigar_ops(cigar)):
        payload = tok.payload
        if payload and tok.op in {'I', 'X'}:
            payload = _reverse_mapping_payload(tok)
        ops.append(CigarOp(tok.n, tok.op, payload))
    return ''.join((_op_to_cigar(tok) for tok in ops))

def _segment_query_length(segment: GraphicSegment) -> int:
    return sum((query_consume(tok.op, tok.n) for tok in _strip_terminal_h(segment.ops)))

def ordered_valid_unified_segments(unified_segments: Sequence[UnifiedSubsegment], shared_ids: Set[int], reverse_view: bool=False) -> List[UnifiedSubsegment]:
    """Return predefined valid chunks in original segment order.

                Do not drop unshared chunks here. Failed chunks must remain visible to
                LCS/Stage2 so they can become D/I boundaries with their attached leftovers.
                Do not reverse the ref order before LCS; final strand handling belongs to
                projection/serialization, not token ordering.
                """
    del shared_ids, reverse_view
    by_original: Dict[int, List[UnifiedSubsegment]] = {}
    order: List[int] = []
    for seg in unified_segments:
        if seg.original_index not in by_original:
            order.append(seg.original_index)
            by_original[seg.original_index] = []
        by_original[seg.original_index].append(seg)
    flat: List[UnifiedSubsegment] = []
    for orig_idx in order:
        parts = by_original[orig_idx]
        parts.sort(key=lambda x: (x.group_start, x.group_end))
        if parts and parts[0].strand == '<':
            parts = list(reversed(parts))
        flat.extend(parts)
    return flat

def weighted_lcs(ref_ids: Sequence[int], qry_ids: Sequence[int], weight_by_id: Dict[int, int]) -> List[Tuple[int, int]]:
    n = len(ref_ids)
    m = len(qry_ids)
    dp = [[0] * (m + 1) for _ in range(n + 1)]
    for i in range(n - 1, -1, -1):
        for j in range(m - 1, -1, -1):
            best = dp[i + 1][j] if dp[i + 1][j] >= dp[i][j + 1] else dp[i][j + 1]
            if ref_ids[i] == qry_ids[j]:
                w = weight_by_id.get(ref_ids[i], 1)
                cand = dp[i + 1][j + 1] + w
                if cand > best:
                    best = cand
            dp[i][j] = best
    out: List[Tuple[int, int]] = []
    i = j = 0
    while i < n and j < m:
        if ref_ids[i] == qry_ids[j]:
            w = weight_by_id.get(ref_ids[i], 1)
            if dp[i][j] == dp[i + 1][j + 1] + w:
                out.append((i, j))
                i += 1
                j += 1
                continue
        if dp[i + 1][j] >= dp[i][j + 1]:
            i += 1
        else:
            j += 1
    return out

def _same_path_for_stage1(a: GraphicSegment, b: GraphicSegment) -> bool:
    return path_key(a.path) == path_key(b.path) and a.direction == b.direction

def _segment_local_cuts(ref_segments: Sequence[GraphicSegment], qry_segments: Sequence[GraphicSegment]) -> Tuple[List[List[int]], List[List[int]]]:
    """
        Build Stage1 cuts using the original pairwise ref-vs-query overlap rule,
        then share only reference-visible breakpoints across other occurrences of
        the same graph path.

        Why this is different from the previous graph-wide version:
          * Query-only endpoints do not create LCS chunks.  They stay as local
            leftovers/insertions and later go through duplicate-LCS or fallback.
          * A breakpoint becomes shared across duplicated query occurrences only
            after it is visible on a reference-overlapping chunk of the same graph.

        Example:
          ref P:342-131669
          qry P:342-125152 and duplicate P:97138-125152

        The duplicate creates a reference-visible cut at 97138 on P, and that cut is
        then also applied to qry P:342-125152, restoring the old primary LCS:
        342-97138 followed by 97138-125152.  But a query-only path with coordinates
        59-65-6521 and no reference P match will not be split into a tiny 59-65
        failed LCS chunk.
        """
    ref_cuts: List[List[int]] = [[seg.start, seg.end] for seg in ref_segments]
    qry_cuts: List[List[int]] = [[seg.start, seg.end] for seg in qry_segments]
    for ri, rseg in enumerate(ref_segments):
        for qi, qseg in enumerate(qry_segments):
            if not _same_path_for_stage1(rseg, qseg):
                continue
            ov0 = max(rseg.start, qseg.start)
            ov1 = min(rseg.end, qseg.end)
            if ov1 <= ov0:
                continue
            for x in (ov0, ov1):
                if rseg.start <= x <= rseg.end:
                    ref_cuts[ri].append(x)
                if qseg.start <= x <= qseg.end:
                    qry_cuts[qi].append(x)
    ref_visible_points: Dict[Tuple[str, str], Set[int]] = {}
    for seg, cuts in zip(ref_segments, ref_cuts):
        key = (path_key(seg.path), seg.direction)
        ref_visible_points.setdefault(key, set()).update(cuts)

    def add_ref_visible_points(cuts_by_segment: List[List[int]], segments: Sequence[GraphicSegment]) -> None:
        for i, seg in enumerate(segments):
            key = (path_key(seg.path), seg.direction)
            for x in ref_visible_points.get(key, set()):
                if seg.start <= x <= seg.end:
                    cuts_by_segment[i].append(x)
    add_ref_visible_points(ref_cuts, ref_segments)
    add_ref_visible_points(qry_cuts, qry_segments)
    ref_cuts = [sorted(set(xs)) for xs in ref_cuts]
    qry_cuts = [sorted(set(xs)) for xs in qry_cuts]
    return (ref_cuts, qry_cuts)

def _truncate_segments_from_local_cuts(side: str, segments: Sequence[GraphicSegment], cuts_by_segment: Sequence[Sequence[int]]) -> List[UnifiedSubsegment]:
    """
        Create UnifiedSubsegment objects from segment-local cuts.

        group_start/group_end are now local coordinate cuts, not global path-name groups.
        """
    out: List[UnifiedSubsegment] = []
    for seg_i, (seg, cuts) in enumerate(zip(segments, cuts_by_segment)):
        if len(cuts) < 2:
            continue
        for a, b in zip(cuts, cuts[1:]):
            if b <= a:
                continue
            out.append(UnifiedSubsegment(side=side, original_index=seg_i, path=seg.path, strand=seg.direction, group_start=a, group_end=b, coord_start=a, coord_end=b))
    return out

def _stage2_coord_start(seg: UnifiedSubsegment) -> int:
    return seg.coord_start if seg.stage2_start is None else seg.stage2_start

def _stage2_coord_end(seg: UnifiedSubsegment) -> int:
    return seg.coord_end if seg.stage2_end is None else seg.stage2_end

def _init_stage2_interval(seg: UnifiedSubsegment) -> None:
    if seg.stage2_start is None:
        seg.stage2_start = seg.coord_start
    if seg.stage2_end is None:
        seg.stage2_end = seg.coord_end

def _suppress_terminal_liftovers_for_lcs(chunks: Sequence[UnifiedSubsegment], min_span: int=_STAGE1_TERMINAL_LIFTOVER_MAX_SPAN) -> List[UnifiedSubsegment]:
    """Remove tiny terminal chunks from LCS and attach them to a neighbor.

                The core coord_start/coord_end of retained chunks stays unchanged, so token
                matching/LCS still sees only the reliable interval. The stage2_start/end of
                the retained chunk expands over suppressed terminal fragments, so Stage2 later
                prepends/appends the original sliced CIGAR to that nearest chunk before the
                same-template comparison.
                """
    for chunk in chunks:
        _init_stage2_interval(chunk)
    by_original: Dict[int, List[UnifiedSubsegment]] = {}
    for chunk in chunks:
        by_original.setdefault(chunk.original_index, []).append(chunk)
    drop_ids: Set[int] = set()
    attachments: List[Tuple[UnifiedSubsegment, UnifiedSubsegment, str]] = []
    for _orig, parts in by_original.items():
        parts.sort(key=lambda x: (x.coord_start, x.coord_end))
        if len(parts) <= 1:
            continue
        left = 0
        right = len(parts) - 1
        while left < right and parts[left].span < min_span:
            left += 1
        while right > left and parts[right].span < min_span:
            right -= 1
        if left > 0:
            keeper = parts[left]
            _init_stage2_interval(keeper)
            keeper.stage2_start = min(_stage2_coord_start(keeper), parts[0].coord_start)
            for dropped in parts[:left]:
                drop_ids.add(id(dropped))
                attachments.append((dropped, keeper, 'prefix'))
        if right < len(parts) - 1:
            keeper = parts[right]
            _init_stage2_interval(keeper)
            keeper.stage2_end = max(_stage2_coord_end(keeper), parts[-1].coord_end)
            for dropped in parts[right + 1:]:
                drop_ids.add(id(dropped))
                attachments.append((dropped, keeper, 'suffix'))
    if not drop_ids:
        return list(chunks)
    return [chunk for chunk in chunks if id(chunk) not in drop_ids]

def _assign_token_ids_by_path_interval(ref_unified: List[UnifiedSubsegment], qry_unified: List[UnifiedSubsegment]) -> None:
    """
        Assign same token to the same path occurrence interval.

        This still allows ref/qry chunks covering the same path interval to match,
        but does not let one query occurrence inherit another occurrence's cuts.
        """
    key_to_id: Dict[Tuple[str, str, int, int], int] = {}
    next_id = 1
    for seg in list(ref_unified) + list(qry_unified):
        key = (path_key(seg.path), seg.strand, seg.coord_start, seg.coord_end)
        if key not in key_to_id:
            key_to_id[key] = next_id
            next_id += 1
        seg.token_id = key_to_id[key]

def build_stage1(ref_gcigar: str, qry_gcigar: str, ref_name: str, qry_name: str, max_merge_gap: int=100, reverse_ref_view: bool=False) -> Stage1Result:
    del max_merge_gap
    ref_segments = parse_graphic_segments(ref_gcigar, ref_name)
    qry_segments = parse_graphic_segments(qry_gcigar, qry_name)
    ref_cuts, qry_cuts = _segment_local_cuts(ref_segments, qry_segments)
    ref_unified = _truncate_segments_from_local_cuts('ref', ref_segments, ref_cuts)
    qry_unified = _truncate_segments_from_local_cuts('qry', qry_segments, qry_cuts)
    ref_unified = _suppress_terminal_liftovers_for_lcs(ref_unified)
    qry_unified = _suppress_terminal_liftovers_for_lcs(qry_unified)
    _assign_token_ids_by_path_interval(ref_unified, qry_unified)
    ref_shared_ids = {seg.token_id for seg in ref_unified}
    qry_shared_ids = {seg.token_id for seg in qry_unified}
    shared_ids = ref_shared_ids & qry_shared_ids
    ref_ordered = ordered_valid_unified_segments(ref_unified, shared_ids, reverse_view=False)
    qry_ordered = ordered_valid_unified_segments(qry_unified, shared_ids, reverse_view=False)
    ref_ids = [seg.token_id for seg in ref_ordered]
    qry_ids = [seg.token_id for seg in qry_ordered]
    weight_by_id: Dict[int, int] = {}
    for seg in ref_unified:
        weight_by_id[seg.token_id] = seg.span
    for seg in qry_unified:
        weight_by_id.setdefault(seg.token_id, seg.span)
    lcs_pairs = weighted_lcs(ref_ids, qry_ids, weight_by_id)
    return Stage1Result(ref_segments=ref_segments, qry_segments=qry_segments, unified_breaks={}, ref_unified=ref_unified, qry_unified=qry_unified, ref_tokens=ref_ids, qry_tokens=qry_ids, lcs_pairs=lcs_pairs, ref_ordered=ref_ordered, qry_ordered=qry_ordered, reverse_ref_view=reverse_ref_view)

def _strip_terminal_h(ops: Sequence[CigarOp]) -> List[CigarOp]:
    body = list(ops)
    if body and body[0].op == 'H':
        body = body[1:]
    if body and body[-1].op == 'H':
        body = body[:-1]
    return body

def _canonical_body_ops(segment: GraphicSegment) -> List[CigarOp]:
    body = _strip_terminal_h(segment.ops)
    if segment.direction != '<':
        return body
    out: List[CigarOp] = []
    for tok in reversed(body):
        payload = tok.payload
        if payload and tok.op == 'I':
            payload = revcomp(payload)
        elif payload and tok.op == 'X':
            payload = _revcomp_x_payload(payload, tok.n)
        out.append(CigarOp(tok.n, tok.op, payload))
    return out

def _revcomp_x_payload(payload: str, n: int) -> str:
    if n <= 0 or not payload:
        return payload
    kind = _x_payload_kind(payload, n)
    if kind == 'query':
        return revcomp(payload)
    if kind == 'block':
        ref_payload = payload[:n]
        qry_payload = payload[n:]
        return revcomp(ref_payload) + revcomp(qry_payload)
    return payload

def build_matched_runs(stage1: Stage1Result) -> List[MatchedRun]:
    runs: List[MatchedRun] = []
    for ref_idx, qry_idx in stage1.lcs_pairs:
        ref_seg = stage1.ref_ordered[ref_idx]
        qry_seg = stage1.qry_ordered[qry_idx]
        if not runs:
            runs.append(MatchedRun(path=ref_seg.path, strand=ref_seg.strand, ref_original_index=ref_seg.original_index, qry_original_index=qry_seg.original_index, group_start=ref_seg.group_start, group_end=ref_seg.group_end, coord_start=ref_seg.coord_start, coord_end=ref_seg.coord_end))
            continue
        last = runs[-1]
        if last.path == ref_seg.path and last.strand == ref_seg.strand and (last.ref_original_index == ref_seg.original_index) and (last.qry_original_index == qry_seg.original_index) and (last.group_end == ref_seg.group_start) and (last.group_end == qry_seg.group_start):
            last.group_end = ref_seg.group_end
            last.coord_end = ref_seg.coord_end
        else:
            runs.append(MatchedRun(path=ref_seg.path, strand=ref_seg.strand, ref_original_index=ref_seg.original_index, qry_original_index=qry_seg.original_index, group_start=ref_seg.group_start, group_end=ref_seg.group_end, coord_start=ref_seg.coord_start, coord_end=ref_seg.coord_end))
    return runs

def _tm_parse_template_cigar(cigar: str) -> List[Tuple[int, str, str]]:
    out: List[Tuple[int, str, str]] = []
    pos = 0
    pat = re.compile('(\\d+)([MXID=H])([A-Za-z]*)')
    while pos < len(cigar):
        m = pat.match(cigar, pos)
        if not m:
            raise ValueError(f'Cannot parse CIGAR near: {cigar[pos:pos + 80]}')
        n = int(m.group(1))
        op = m.group(2)
        payload = m.group(3) or ''
        if op == '=':
            op = 'M'
        if payload and op not in {'X', 'I', 'D'}:
            raise ValueError(f'Unexpected payload on {n}{op}: {payload}')
        if op in {'I', 'D'} and payload and (len(payload) != n):
            pass
        if op == 'X' and payload and (len(payload) not in {n, 2 * n}):
            raise ValueError(f'Payload length mismatch for {n}{op}: payload length={len(payload)}')
        out.append((n, op, payload))
        pos = m.end()
    return out

def _tm_template_slice(template_seq: Optional[str], r0: int, r1: int, template_offset: int=0) -> str:
    n = max(0, r1 - r0)
    if n == 0:
        return ''
    if template_seq is None:
        return 'N' * n
    s = r0 - template_offset
    e = r1 - template_offset
    left_pad = max(0, -s)
    right_pad = max(0, e - len(template_seq))
    s = max(0, s)
    e = min(len(template_seq), e)
    return 'N' * left_pad + template_seq[s:e] + 'N' * right_pad

def _tm_side_x_payload(n: int, payload: str, side_role: str, template_seq: Optional[str], rpos: int, template_offset: int=0) -> str:
    """Return the side/sample bases represented by an X op.

                Important: merge_cigars_on_same_template() compares two side tracks
                against a shared template.  Its arguments are not the final ref/query
                CIGAR.  Therefore an X payload must be interpreted as the bases on that
                side, independent of whether the side is called "ref" or "qry" in the
                merge.  For pairwise/block X payloads (ref_bases + query_bases), the
                side/sample bases are the query half.  Using the first half for the
                "ref" track makes identical side CIGARs compare as mismatches.
                """
    del side_role, template_seq, rpos, template_offset
    if not payload:
        return 'N' * n
    if len(payload) == 2 * n:
        return payload[n:]
    if len(payload) == n:
        return payload
    raise ValueError(f'X payload length mismatch: {n}X has payload length {len(payload)}')

def _tm_side_i_payload(n: int, payload: str) -> str:
    if not payload:
        return 'N' * n
    if len(payload) == n:
        return payload
    raise ValueError(f'I payload length mismatch: {n}I has payload length {len(payload)}')

def _tm_build_side_track(cigar: str, template_seq: Optional[str]=None, side_seq: Optional[str]=None, template_offset: int=0, side_role: str='qry') -> TemplateMergeSideTrack:
    intervals: List[TemplateMergeSideInterval] = []
    insertions: Dict[int, TemplateMergeSideInsertion] = {}
    rpos = 0
    qpos = 0
    for n, op, payload in _tm_parse_template_cigar(cigar):
        if op == 'M':
            if side_seq is not None:
                seq = side_seq[qpos:qpos + n]
                if len(seq) != n:
                    seq = seq.ljust(n, 'N')
            else:
                seq = _tm_template_slice(template_seq, rpos, rpos + n, template_offset)
            intervals.append(TemplateMergeSideInterval(rpos, rpos + n, 'M', seq, qpos, qpos + n))
            rpos += n
            qpos += n
        elif op == 'X':
            if payload:
                seq = _tm_side_x_payload(n, payload, side_role=side_role, template_seq=template_seq, rpos=rpos, template_offset=template_offset)
            elif side_seq is not None:
                seq = side_seq[qpos:qpos + n]
                if len(seq) != n:
                    seq = seq.ljust(n, 'N')
            else:
                seq = _tm_side_x_payload(n, payload, side_role=side_role, template_seq=template_seq, rpos=rpos, template_offset=template_offset)
            intervals.append(TemplateMergeSideInterval(rpos, rpos + n, 'X', seq, qpos, qpos + n))
            rpos += n
            qpos += n
        elif op in {'D', 'H'}:
            intervals.append(TemplateMergeSideInterval(rpos, rpos + n, 'D', None, qpos, qpos))
            rpos += n
        elif op == 'I':
            if payload and payload != 'n':
                seq = _tm_side_i_payload(n, payload)
                has_payload = True
            elif side_seq is not None:
                seq = side_seq[qpos:qpos + n]
                if len(seq) != n:
                    seq = seq.ljust(n, 'N')
                has_payload = True
            elif payload == 'n':
                seq = 'N' * n
                has_payload = False
            else:
                seq = 'N' * n
                has_payload = False
            if rpos in insertions:
                old = insertions[rpos]
                insertions[rpos] = TemplateMergeSideInsertion(seq=old.seq + seq, q0=old.q0, q1=qpos + n, has_payload=old.has_payload and has_payload)
            else:
                insertions[rpos] = TemplateMergeSideInsertion(seq=seq, q0=qpos, q1=qpos + n, has_payload=has_payload)
            qpos += n
    return TemplateMergeSideTrack(intervals=intervals, insertions=insertions, r_end=rpos, q_end=qpos)

def _tm_find_interval(intervals: Sequence[TemplateMergeSideInterval], a: int, b: int) -> Optional[TemplateMergeSideInterval]:
    for iv in intervals:
        if iv.r0 <= a and b <= iv.r1:
            return iv
    return None

def _tm_interval_seq(iv: Optional[TemplateMergeSideInterval], a: int, b: int) -> str:
    if iv is None or iv.op == 'D':
        return ''
    if iv.seq is None:
        return 'N' * (b - a)
    off = a - iv.r0
    return iv.seq[off:off + (b - a)]

def _tm_add_op(out: List[Tuple[int, str, str]], n: int, op: str, payload: str='') -> None:
    if n <= 0:
        return
    if op == 'M':
        payload = ''
    if out and out[-1][1] == op:
        old_n, old_op, old_payload = out[-1]
        if op == 'M':
            out[-1] = (old_n + n, old_op, '')
        else:
            merged = _merge_payload_for_same_op(op, old_n, old_payload, n, payload)
            out[-1] = (old_n + n, old_op, merged)
    else:
        out.append((n, op, payload))

def _match_mismatch_runs(ref_upper: str, qry_upper: str) -> List[Tuple[int, int, bool]]:
    """Split two equal-length uppercased sequences into (start, end, is_match) runs.

    Emitting one op per run through the existing merge helpers produces output
    byte-identical to the historical per-base loops: M/= payloads are empty, and
    both X payload encodings (block ref+qry and query-only) merge adjacent runs
    by concatenating per side, which is associative.
    """
    m = len(ref_upper)
    runs: List[Tuple[int, int, bool]] = []
    i = 0
    while i < m:
        is_match = ref_upper[i] == qry_upper[i]
        j = i + 1
        if is_match:
            while j < m:
                k = min(j + 512, m)
                if ref_upper[j:k] == qry_upper[j:k]:
                    j = k
                    continue
                while j < k and ref_upper[j] == qry_upper[j]:
                    j += 1
                break
        else:
            while j < m and ref_upper[j] != qry_upper[j]:
                j += 1
        runs.append((i, j, is_match))
        i = j
    return runs

def _tm_compare_equal_order(ref_seq: str, qry_seq: str, out: List[Tuple[int, str, str]]) -> None:
    m = min(len(ref_seq), len(qry_seq))
    ref_upper = ref_seq[:m].upper()
    qry_upper = qry_seq[:m].upper()
    if ref_upper == qry_upper:
        _tm_add_op(out, m, 'M')
    else:
        for i, j, is_match in _match_mismatch_runs(ref_upper, qry_upper):
            if is_match:
                _tm_add_op(out, j - i, 'M')
            else:
                _tm_add_op(out, j - i, 'X', ref_seq[i:j] + qry_seq[i:j])
    if len(ref_seq) > m:
        _tm_add_op(out, len(ref_seq) - m, 'D', ref_seq[m:])
    if len(qry_seq) > m:
        _tm_add_op(out, len(qry_seq) - m, 'I', qry_seq[m:])
# Payload pairs whose longer side is below this use the in-process exact
# global DP (parasail-accelerated); larger pairs go to minimap2.  The old
# 100 bp threshold predates the fast DP: a subprocess spawn costs ~6 ms,
# while the DP now aligns a sub-1kb pair in well under a millisecond, and
# minimap2 seeding frequently finds no hit below ~1 kb anyway.
_PAYLOAD_MINIMAP2_THRESHOLD = 1000

def _tm_ops_from_aligned_strings(aln_ref: str, aln_qry: str) -> List[Tuple[int, str, str]]:
    out: List[Tuple[int, str, str]] = []
    if len(aln_ref) != len(aln_qry):
        raise ValueError('swspy aligned strings have different lengths')
    ref_upper = aln_ref.upper()
    qry_upper = aln_qry.upper()
    m = len(aln_ref)
    i = 0
    while i < m:
        ref_dash = aln_ref[i] == '-'
        qry_dash = aln_qry[i] == '-'
        if ref_dash and qry_dash:
            i += 1
            continue
        j = i + 1
        if ref_dash:
            while j < m and aln_ref[j] == '-' and aln_qry[j] != '-':
                j += 1
            _tm_add_op(out, j - i, 'I', aln_qry[i:j])
        elif qry_dash:
            while j < m and aln_qry[j] == '-' and aln_ref[j] != '-':
                j += 1
            _tm_add_op(out, j - i, 'D', '')
        elif ref_upper[i] == qry_upper[i]:
            while (
                j < m
                and aln_ref[j] != '-'
                and aln_qry[j] != '-'
                and ref_upper[j] == qry_upper[j]
            ):
                j += 1
            _tm_add_op(out, j - i, '=', '')
        else:
            while (
                j < m
                and aln_ref[j] != '-'
                and aln_qry[j] != '-'
                and ref_upper[j] != qry_upper[j]
            ):
                j += 1
            _tm_add_op(out, j - i, 'X', aln_ref[i:j] + aln_qry[i:j])
        i = j
    return out

def _tm_global_insert_align_dp(ref_seq: str, qry_seq: str) -> Tuple[List[Tuple[int, str, str]], int, int, int, int]:
    """Global DP alignment of insertion payloads in forward orientation.

    Return ops plus aligned ref/query start/end coordinates. D has no payload;
    I carries query payload; X carries ref+query payload.

    The matrix fill is vectorized with numpy when available, replicating the
    per-cell recurrence exactly: match +2 / mismatch -3 / gap -4 with strict
    ">" tie-breaking that prefers diagonal, then delete, then insert.  The
    insert chain curr[j] = max(cand[j], curr[j-1] - 4) is computed as a
    running maximum of cand[j] + 4j (a prefix max with linear decay).
    """
    n = len(ref_seq)
    m = len(qry_seq)
    match_score = 2
    mismatch_score = -3
    gap_score = -4
    if n == 0 or m == 0:
        trivial: List[Tuple[int, str, str]] = []
        if n:
            _tm_add_op(trivial, n, 'D', '')
        if m:
            _tm_add_op(trivial, m, 'I', qry_seq)
        return (trivial, 0, 0, n, m)
    if _parasail is not None:
        # parasail computes only the score TABLE (SIMD, same recurrence:
        # +2/-3 matrix with open=extend=4 equals this DP's linear gap model);
        # the traceback below recomputes this function's own tie-breaking
        # (diagonal > delete > insert, strict ">") from the table values, so
        # the emitted ops are byte-identical to the numpy/pure-python paths.
        ref_upper_text = ref_seq.upper()
        qry_upper_text = qry_seq.upper()
        alphabet = ''.join(sorted(set(ref_upper_text) | set(qry_upper_text)))
        matrix = _PARASAIL_MATRIX_CACHE.get(alphabet)
        if matrix is None:
            matrix = _parasail.matrix_create(
                alphabet, match_score, mismatch_score,
            )
            _PARASAIL_MATRIX_CACHE[alphabet] = matrix
        # The score_table is a view into the result's C buffer; the result
        # object must stay referenced for as long as the table is read.
        parasail_result = _parasail.nw_table(
            qry_upper_text, ref_upper_text, -gap_score, -gap_score, matrix,
        )
        table = parasail_result.score_table

        def cell(i: int, j: int) -> int:
            # Full-matrix accessor over parasail's boundary-free table:
            # table[query-1, ref-1]; row/column zero follow the gap ramp.
            if i == 0:
                return j * gap_score
            if j == 0:
                return i * gap_score
            return int(table[j - 1, i - 1])

        rev_steps: List[Tuple[int, str, str]] = []
        i = n
        j = m
        while i > 0 or j > 0:
            if i > 0 and j > 0:
                here = cell(i, j)
                diag = cell(i - 1, j - 1) + (
                    match_score
                    if ref_upper_text[i - 1] == qry_upper_text[j - 1]
                    else mismatch_score
                )
                delete = cell(i - 1, j) + gap_score
                if here > diag and here > delete:
                    bt = 2
                elif here > diag:
                    bt = 1
                else:
                    bt = 0
            else:
                bt = 1 if j == 0 else 2
            if i > 0 and j > 0 and bt == 0:
                rb = ref_seq[i - 1]
                qb = qry_seq[j - 1]
                if rb.upper() == qb.upper():
                    rev_steps.append((1, '=', ''))
                else:
                    rev_steps.append((1, 'X', rb + qb))
                i -= 1
                j -= 1
            elif i > 0 and (j == 0 or bt == 1):
                rev_steps.append((1, 'D', ''))
                i -= 1
            elif j > 0:
                rev_steps.append((1, 'I', qry_seq[j - 1]))
                j -= 1
            else:
                break
        ops = []
        for n0, op, payload in reversed(rev_steps):
            _tm_add_op(ops, n0, op, payload)
        return (ops, 0, 0, n, m)
    if _np is not None and n > 0 and m > 0:
        ref_upper = _np.frombuffer(
            ref_seq.upper().encode('latin-1'), dtype=_np.uint8,
        )
        qry_upper = _np.frombuffer(
            qry_seq.upper().encode('latin-1'), dtype=_np.uint8,
        )
        trace = _np.zeros((n + 1, m + 1), dtype=_np.uint8)
        trace[0, 1:] = 2
        trace[1:, 0] = 1
        prev = gap_score * _np.arange(m + 1, dtype=_np.int64)
        decay = -gap_score * _np.arange(m + 1, dtype=_np.int64)
        chain = _np.empty(m + 1, dtype=_np.int64)
        for i in range(1, n + 1):
            sub_row = _np.where(
                qry_upper == ref_upper[i - 1], match_score, mismatch_score,
            )
            diag = prev[:-1] + sub_row
            delete = prev[1:] + gap_score
            cand = _np.maximum(diag, delete)
            row_trace = trace[i]
            row_trace[1:][delete > diag] = 1
            chain[0] = i * gap_score
            _np.add(cand, decay[1:], out=chain[1:])
            _np.maximum.accumulate(chain, out=chain)
            curr = chain - decay
            row_trace[1:][curr[1:] > cand] = 2
            prev = curr
    else:
        trace = [bytearray(m + 1) for _ in range(n + 1)]
        prev = [j * gap_score for j in range(m + 1)]
        for j in range(1, m + 1):
            trace[0][j] = 2
        for i in range(1, n + 1):
            curr = [i * gap_score] + [0] * m
            trace[i][0] = 1
            rb = ref_seq[i - 1]
            for j in range(1, m + 1):
                qb = qry_seq[j - 1]
                diag = prev[j - 1] + (match_score if rb.upper() == qb.upper() else mismatch_score)
                delete = prev[j] + gap_score
                insert = curr[j - 1] + gap_score
                best = diag
                bt = 0
                if delete > best:
                    best = delete
                    bt = 1
                if insert > best:
                    best = insert
                    bt = 2
                curr[j] = best
                trace[i][j] = bt
            prev = curr
    rev_steps: List[Tuple[int, str, str]] = []
    i = n
    j = m
    while i > 0 or j > 0:
        bt = trace[i][j]
        if i > 0 and j > 0 and (bt == 0):
            rb = ref_seq[i - 1]
            qb = qry_seq[j - 1]
            if rb.upper() == qb.upper():
                rev_steps.append((1, '=', ''))
            else:
                rev_steps.append((1, 'X', rb + qb))
            i -= 1
            j -= 1
        elif i > 0 and (j == 0 or bt == 1):
            rev_steps.append((1, 'D', ''))
            i -= 1
        elif j > 0:
            rev_steps.append((1, 'I', qry_seq[j - 1]))
            j -= 1
        else:
            break
    ops: List[Tuple[int, str, str]] = []
    for n0, op, payload in reversed(rev_steps):
        _tm_add_op(ops, n0, op, payload)
    return (ops, 0, 0, len(ref_seq), len(qry_seq))

def _tm_ops_from_cigar_on_sequences(ref_seq: str, qry_seq: str, cigar: str, ref_start: int=0, qry_start: int=0) -> List[Tuple[int, str, str]]:
    """Convert an aligner CIGAR into pairwise ops with payloads.

    The aligner CIGAR is interpreted on forward-strand ref/query strings. If the
    aligner starts internally, leading ref-only/query-only sequence is emitted as
    D/I so the returned CIGAR is global over the full payloads.
    """
    out: List[Tuple[int, str, str]] = []
    r = max(0, int(ref_start))
    q = max(0, int(qry_start))
    if r > 0:
        _tm_add_op(out, r, 'D', '')
    if q > 0:
        _tm_add_op(out, q, 'I', qry_seq[:q])
    for n_text, op, payload in _PAYLOAD_RE.findall(cigar):
        n = int(n_text)
        if n <= 0:
            continue
        if op == 'M':
            op = '='
        if op in {'=', 'X'}:
            ref_chunk = ref_seq[r:r + n]
            if len(ref_chunk) < n:
                ref_chunk += 'N' * (n - len(ref_chunk))
            qry_chunk = qry_seq[q:q + n]
            if len(qry_chunk) < n:
                qry_chunk += 'N' * (n - len(qry_chunk))
            ref_upper = ref_chunk.upper()
            qry_upper = qry_chunk.upper()
            if ref_upper == qry_upper:
                _tm_add_op(out, n, '=', '')
            else:
                for i, j, is_match in _match_mismatch_runs(
                    ref_upper, qry_upper,
                ):
                    if is_match:
                        _tm_add_op(out, j - i, '=', '')
                    else:
                        _tm_add_op(
                            out, j - i, 'X',
                            ref_chunk[i:j] + qry_chunk[i:j],
                        )
            r += n
            q += n
        elif op == 'D':
            _tm_add_op(out, n, 'D', '')
            r += n
        elif op == 'I':
            seq = payload if payload and len(payload) == n else qry_seq[q:q + n]
            if len(seq) < n:
                seq += 'N' * (n - len(seq))
            _tm_add_op(out, n, 'I', seq[:n])
            q += n
        elif op == 'S':
            seq = qry_seq[q:q + n]
            if len(seq) < n:
                seq += 'N' * (n - len(seq))
            _tm_add_op(out, n, 'I', seq[:n])
            q += n
        elif op == 'H':
            continue
        else:
            raise ValueError(f'unsupported SSW CIGAR op {op!r}')
    if r < len(ref_seq):
        _tm_add_op(out, len(ref_seq) - r, 'D', '')
    if q < len(qry_seq):
        _tm_add_op(out, len(qry_seq) - q, 'I', qry_seq[q:])
    return out

def _result_get(result, *names, default=None):
    if isinstance(result, dict):
        for name in names:
            if name in result:
                return result[name]
    for name in names:
        if hasattr(result, name):
            return getattr(result, name)
    return default

def _normalize_swspy_result(ref_seq: str, qry_seq: str, result) -> List[Tuple[int, str, str]]:
    if isinstance(result, str):
        return _tm_ops_from_cigar_on_sequences(ref_seq, qry_seq, result, 0, 0)
    if isinstance(result, list):
        out: List[Tuple[int, str, str]] = []
        for n0, op, payload in result:
            if op == '=':
                op = 'M'
            if op == 'D':
                payload = ''
            _tm_add_op(out, int(n0), op, payload)
        return out
    if isinstance(result, tuple):
        if len(result) == 3 and isinstance(result[2], str):
            return _tm_ops_from_cigar_on_sequences(ref_seq, qry_seq, result[2], int(result[0]), int(result[1]))
        if len(result) >= 5 and isinstance(result[-1], str):
            return _tm_ops_from_cigar_on_sequences(ref_seq, qry_seq, result[-1], int(result[-5]), int(result[-3]))
        if len(result) == 2 and isinstance(result[0], str) and isinstance(result[1], str):
            return _tm_ops_from_aligned_strings(result[0], result[1])
    aln_ref = _result_get(result, 'aligned_ref', 'ref_aln', 'target_aln')
    aln_qry = _result_get(result, 'aligned_query', 'query_aln', 'qry_aln')
    if aln_ref is not None and aln_qry is not None:
        return _tm_ops_from_aligned_strings(str(aln_ref), str(aln_qry))
    cigar = _result_get(result, 'cigar', 'cigar_string')
    if cigar is not None:
        ref_start = _result_get(result, 'ref_start', 'rstart', 'target_start', 'start_ref', default=0)
        qry_start = _result_get(result, 'query_start', 'qry_start', 'qstart', 'start_query', default=0)
        return _tm_ops_from_cigar_on_sequences(ref_seq, qry_seq, str(cigar), int(ref_start), int(qry_start))
    raise RuntimeError('swspy() must return either normalized ops, a CIGAR string, (ref_start, query_start, cigar), aligned strings, or an object/dict with cigar plus ref_start/query_start')

def _tm_ssw_cigar_text(aln) -> str:
    """Return a CIGAR string from an ssw-py alignment object."""
    cigar = _result_get(aln, 'CIGAR', 'cigar', 'cigar_string', 'cigar_pair_list', default='')
    if cigar is None:
        return ''
    if isinstance(cigar, str):
        return cigar
    if isinstance(cigar, bytes):
        return cigar.decode()
    if isinstance(cigar, (list, tuple)):
        parts: List[str] = []
        op_map = {0: 'M', 1: 'I', 2: 'D', 4: 'S', 5: 'H', 7: '=', 8: 'X'}
        for item in cigar:
            if isinstance(item, (list, tuple)) and len(item) >= 2:
                a, b = (item[0], item[1])
                if isinstance(a, str):
                    parts.append(f'{int(b)}{a}')
                elif isinstance(b, str):
                    parts.append(f'{int(a)}{b}')
                elif int(a) in op_map:
                    parts.append(f'{int(b)}{op_map[int(a)]}')
                else:
                    parts.append(f"{int(a)}{op_map.get(int(b), 'M')}")
        return ''.join(parts)
    return str(cigar)

def swspy(ref_seq: str, qry_seq: str) -> List[Tuple[int, str, str]]:
    """Long insertion-payload aligner adapter.

    The function name is intentionally ``swspy`` to keep the call site isolated,
    but the implementation uses the installed ``ssw`` / ssw-py package, matching
    the FYI script.  It aligns ``qry_seq`` to ``ref_seq`` in forward orientation
    and converts the local SSW result into a global pairwise CIGAR over the full
    insertion payloads by adding leading/trailing D/I as needed.
    """
    try:
        from ssw import AlignmentMgr
        try:
            from ssw.alignmentmgr import BitwiseAlignmentFlag
        except Exception:
            BitwiseAlignmentFlag = None
    except ImportError as exc:
        raise RuntimeError('Insertion-vs-insertion realignment for payloads >=100 bp requires the ssw / ssw-py package') from exc
    if not ref_seq and (not qry_seq):
        return []
    if not ref_seq:
        return [(len(qry_seq), 'I', qry_seq)] if qry_seq else []
    if not qry_seq:
        return [(len(ref_seq), 'D', '')] if ref_seq else []
    mgr = AlignmentMgr(match_score=2, mismatch_penalty=2)
    mgr.set_read(qry_seq.upper())
    mgr.set_reference(ref_seq.upper())
    kwargs = {'gap_open': 3, 'gap_extension': 1}
    if BitwiseAlignmentFlag is not None:
        kwargs['bitwise_flag'] = int(getattr(BitwiseAlignmentFlag, 'best_idxs', 1))
    aln = mgr.align(**kwargs)
    cigar = _tm_ssw_cigar_text(aln)
    score = _result_get(aln, 'optimal_score', 'score', default=0) or 0
    ref_start = _result_get(aln, 'reference_start', 'ref_begin', 'target_start', default=0)
    qry_start = _result_get(aln, 'read_start', 'query_begin', 'query_start', default=0)
    if not cigar or float(score) <= 0:
        out: List[Tuple[int, str, str]] = []
        _tm_add_op(out, len(ref_seq), 'D', '')
        _tm_add_op(out, len(qry_seq), 'I', qry_seq)
        return out
    return _tm_ops_from_cigar_on_sequences(ref_seq, qry_seq, cigar, int(ref_start), int(qry_start))


def _minimap2_executable() -> Optional[str]:
    """Return the configured minimap2 executable, if it is available."""
    configured = os.environ.get('GRAPH_CIGARTOREF_MINIMAP2', 'minimap2')
    if os.path.isabs(configured) and os.access(configured, os.X_OK):
        return configured
    return shutil.which(configured)


def _minimap2_threads() -> int:
    """Return threads for one payload-alignment subprocess.

    graphcigartoref_persample already parallelizes independent comparisons.
    Keeping each nested minimap2 invocation single-threaded prevents, for
    example, 32 Python workers from silently creating roughly 96 minimap2
    threads.  An advanced caller can still override this for a serial run.
    """
    raw = os.environ.get('GRAPH_CIGARTOREF_MINIMAP2_THREADS', '1')
    try:
        value = int(raw)
    except ValueError as error:
        raise ValueError(
            'GRAPH_CIGARTOREF_MINIMAP2_THREADS must be a positive integer'
        ) from error
    if value <= 0:
        raise ValueError(
            'GRAPH_CIGARTOREF_MINIMAP2_THREADS must be a positive integer'
        )
    return value


def _payload_ops_from_hit(
    ref_seq: str, qry_seq: str, cigar: str, ref_start: int, qry_start: int,
) -> List[Tuple[int, str, str]]:
    ops = _tm_ops_from_cigar_on_sequences(
        ref_seq, qry_seq, cigar, ref_start, qry_start,
    )
    ref_span = sum((ref_consume_pair(op, n) for n, op, _payload in ops))
    qry_span = sum((query_consume(op, n) for n, op, _payload in ops))
    if ref_span != len(ref_seq) or qry_span != len(qry_seq):
        raise RuntimeError(
            'minimap2 payload CIGAR did not cover complete sequences '
            f'(reference {ref_span}/{len(ref_seq)}, query {qry_span}/{len(qry_seq)})'
        )
    return ops


def _mappy_payload_ops(ref_seq: str, qry_seq: str) -> List[Tuple[int, str, str]]:
    """In-process minimap2 payload alignment via mappy.

    Mirrors the subprocess contract: asm5 preset, forward-only primary hits,
    best hit by (matches, block length), and the same no-usable-hit error so
    caller fallbacks behave identically.
    """
    try:
        aligner = _mappy.Aligner(seq=ref_seq, preset='asm5', n_threads=1)
    except Exception as error:
        raise RuntimeError(
            f'minimap2 payload alignment failed: {error}'
        ) from error
    best: Optional[Tuple[int, int, str, int, int]] = None
    for hit in aligner.map(qry_seq):
        if hit.strand != 1 or not hit.is_primary:
            continue
        if hit.q_en <= hit.q_st or hit.r_en <= hit.r_st:
            continue
        cigar = hit.cigar_str
        if not cigar:
            continue
        candidate = (hit.mlen, hit.blen, cigar, hit.r_st, hit.q_st)
        if best is None or candidate[:2] > best[:2]:
            best = candidate
    if best is None:
        raise RuntimeError('minimap2 produced no usable forward payload alignment')
    _matches, _block_len, cigar, ref_start, qry_start = best
    return _payload_ops_from_hit(ref_seq, qry_seq, cigar, ref_start, qry_start)


def minimap2_payload_ops(ref_seq: str, qry_seq: str) -> List[Tuple[int, str, str]]:
    """Align one long insertion payload pair with minimap2.

    The graph-CIGAR merge requires a *global* payload CIGAR, while minimap2
    reports a local PAF alignment.  The PAF ``cg`` tag and its target/query
    starts are therefore expanded with :func:`_tm_ops_from_cigar_on_sequences`
    so unaligned prefixes and suffixes remain explicit D/I operations.  The
    forward-only mode matches the graph-template merge orientation contract.
    Callers conservatively retain the original D/I representation when
    minimap2 cannot produce a usable forward hit.
    """
    if not ref_seq and not qry_seq:
        return []
    if not ref_seq:
        return [(len(qry_seq), 'I', qry_seq)] if qry_seq else []
    if not qry_seq:
        return [(len(ref_seq), 'D', '')] if ref_seq else []
    if _mappy is not None and os.environ.get(
        'GRAPH_CIGARTOREF_NO_MAPPY', '',
    ) in {'', '0'}:
        # In-process minimap2 (same library the CLI wraps): no fork/exec, no
        # temp FASTA files, ~6 ms of fixed cost removed per call.  Set
        # GRAPH_CIGARTOREF_NO_MAPPY=1 to force the subprocess path.
        return _mappy_payload_ops(ref_seq, qry_seq)
    executable = _minimap2_executable()
    if executable is None:
        raise FileNotFoundError(
            'minimap2 is not available; set GRAPH_CIGARTOREF_MINIMAP2 '
            'to its executable path'
        )

    timeout_text = os.environ.get('GRAPH_CIGARTOREF_MINIMAP2_TIMEOUT', '300')
    try:
        timeout = float(timeout_text)
    except ValueError:
        timeout = 300.0
    timeout_value = timeout if timeout > 0 else None

    with tempfile.TemporaryDirectory(prefix='graphcigartoref-mm2-') as directory:
        target_path = os.path.join(directory, 'target.fa')
        query_path = os.path.join(directory, 'query.fa')
        with open(target_path, 'wt') as handle:
            handle.write('>graph_ref\n')
            handle.write(ref_seq)
            handle.write('\n')
        with open(query_path, 'wt') as handle:
            handle.write('>graph_query\n')
            handle.write(qry_seq)
            handle.write('\n')
        command = [
            executable,
            '-x', 'asm5',
            '-t', str(_minimap2_threads()),
            '-c',
            '--eqx',
            '--secondary=no',
            '--for-only',
            target_path,
            query_path,
        ]
        try:
            completed = subprocess.run(
                command,
                check=False,
                stdout=subprocess.PIPE,
                stderr=subprocess.PIPE,
                text=True,
                timeout=timeout_value,
            )
        except (OSError, subprocess.TimeoutExpired) as error:
            raise RuntimeError(f'minimap2 payload alignment failed: {error}') from error

    if completed.returncode != 0:
        message = completed.stderr.strip().splitlines()
        detail = message[-1] if message else f'exit code {completed.returncode}'
        raise RuntimeError(f'minimap2 payload alignment failed: {detail}')

    best: Optional[Tuple[int, int, str, int, int]] = None
    for raw in completed.stdout.splitlines():
        fields = raw.split('\t')
        if len(fields) < 12 or fields[0] != 'graph_query' or fields[5] != 'graph_ref':
            continue
        # The forward-only request should make this true.  Keep the explicit
        # check because a reverse hit cannot be safely projected into the
        # forward payload CIGAR contract without reversing the reference track.
        if fields[4] != '+':
            continue
        try:
            qstart = int(fields[2])
            qend = int(fields[3])
            tstart = int(fields[7])
            tend = int(fields[8])
            matches = int(fields[9])
            block_len = int(fields[10])
        except ValueError:
            continue
        tags: Dict[str, str] = {}
        for tag in fields[12:]:
            parts = tag.split(':', 2)
            if len(parts) == 3:
                tags[parts[0]] = parts[2]
        cigar = tags.get('cg')
        if not cigar or qend <= qstart or tend <= tstart:
            continue
        candidate = (matches, block_len, cigar, tstart, qstart)
        if best is None or candidate[:2] > best[:2]:
            best = candidate
    if best is None:
        raise RuntimeError('minimap2 produced no usable forward payload alignment')

    _matches, _block_len, cigar, ref_start, qry_start = best
    ops = _tm_ops_from_cigar_on_sequences(
        ref_seq, qry_seq, cigar, ref_start, qry_start,
    )
    ref_span = sum((ref_consume_pair(op, n) for n, op, _payload in ops))
    qry_span = sum((query_consume(op, n) for n, op, _payload in ops))
    if ref_span != len(ref_seq) or qry_span != len(qry_seq):
        raise RuntimeError(
            'minimap2 payload CIGAR did not cover complete sequences '
            f'(reference {ref_span}/{len(ref_seq)}, query {qry_span}/{len(qry_seq)})'
        )
    return ops

def _tm_realign_insertion_payloads(ref_seq: str, qry_seq: str, out: List[Tuple[int, str, str]]) -> None:
    """Realign insertion payloads and append pairwise ops to ``out``.

    This is the early, graph-template-aware rescue inherited from the original
    implementation.  When both graph sides carry sequence at the same empty-
    template insertion anchor, align those complete payloads before graph-path
    assembly can fragment them into a later D/.../I pattern.

    Payloads shorter than 1000 bp use the in-process global DP implementation.
    Payloads of at least 1000 bp use minimap2.  Minimap2's local PAF alignment is
    expanded to a global pairwise CIGAR by retaining unaligned prefixes and
    suffixes as D/I.  If minimap2 is unavailable or produces no usable forward
    hit, retain the conservative D/I representation.  The later FASTA-backed
    adjacent and anchored-window rescues remain independent safety nets.
    """
    if not ref_seq and (not qry_seq):
        return
    # Exact long replacements do not need the local aligner.  They commonly
    # arise when two different graph paths carry the same sequence.
    if len(ref_seq) == len(qry_seq) and ref_seq.upper() == qry_seq.upper():
        _tm_add_op(out, len(ref_seq), '=', '')
        return
    max_payload = max(len(ref_seq), len(qry_seq))
    try:
        if max_payload < _PAYLOAD_MINIMAP2_THRESHOLD:
            ops, _ref_start, _qry_start, _ref_end, _qry_end = _tm_global_insert_align_dp(ref_seq, qry_seq)
        else:
            ops = minimap2_payload_ops(ref_seq, qry_seq)
    except (FileNotFoundError, ImportError, RuntimeError):
        # Keep the conservative pre-realignment representation when neither
        # required external aligner is available.  This preserves the old D/I
        # result instead of aborting an otherwise valid graph conversion.
        _tm_add_op(out, len(ref_seq), 'D', '')
        _tm_add_op(out, len(qry_seq), 'I', qry_seq)
        return
    for n0, op, payload in ops:
        if op == 'D':
            payload = ''
        _tm_add_op(out, int(n0), op, payload)

def _tm_compare_insertions(ref_ins: Optional[TemplateMergeSideInsertion], qry_ins: Optional[TemplateMergeSideInsertion], out: List[Tuple[int, str, str]]) -> None:
    ref_seq = ref_ins.seq if ref_ins is not None else ''
    qry_seq = qry_ins.seq if qry_ins is not None else ''
    if not ref_seq and (not qry_seq):
        return
    if ref_seq and (not qry_seq):
        _tm_add_op(out, len(ref_seq), 'D', '')
        return
    if qry_seq and (not ref_seq):
        _tm_add_op(out, len(qry_seq), 'I', qry_seq)
        return
    if ref_ins is not None and qry_ins is not None and (not ref_ins.has_payload or not qry_ins.has_payload):
        _tm_add_op(out, len(ref_seq), 'D', '')
        _tm_add_op(out, len(qry_seq), 'I', qry_seq)
        return
    _tm_realign_insertion_payloads(ref_seq, qry_seq, out)

def _tm_format_pairwise_ops(ops: Sequence[Tuple[int, str, str]]) -> str:
    parts: List[str] = []
    for n, op, payload in ops:
        if n <= 0:
            continue
        emit_op = '=' if op == 'M' else op
        if payload:
            parts.append(f'{n}{emit_op}{payload}')
        else:
            parts.append(f'{n}{emit_op}')
    return ''.join(parts)

_SV_CONSOLIDATE_MAX_GAP = 100
_SV_CONSOLIDATE_FLANK = 50
_SV_CONSOLIDATE_MAX_WINDOW = 5000


_SV_WINDOW_GAP_OPEN = 12
_SV_WINDOW_GAP_EXTEND = 1

# The second consolidation pass has no biological window-size cutoff.  This
# limit only selects the alignment backend: quadratic affine DP is preferred
# for ordinary local windows, while minimap2 is used for very large ones.
_SV_SECONDARY_AFFINE_MAX_CELLS = 100_000_000
_SV_SECONDARY_AFFINE_MAX_AXIS = 50_000


def _sv_window_affine_align(ref_win: str, qry_win: str) -> List[Tuple[int, str, str]]:
    """Affine-gap global alignment for SV-consolidation windows.

    A linear gap model is indifferent between one long gap and many
    fragments of the same total length, which is exactly the tandem-
    repeat shattering this pass exists to undo; a high open / cheap
    extend cost makes the single consolidated indel decisively optimal.
    Falls back to the linear DP when parasail is unavailable.
    """
    if not ref_win:
        out: List[Tuple[int, str, str]] = []
        if qry_win:
            _tm_add_op(out, len(qry_win), 'I', qry_win)
        return out
    if not qry_win:
        out = []
        _tm_add_op(out, len(ref_win), 'D', ref_win)
        return out
    if _parasail is None:
        ops, _r0, _q0, _r1, _q1 = _tm_global_insert_align_dp(ref_win, qry_win)
        return [(n, 'M' if op == '=' else op, p) for n, op, p in ops]
    ref_upper = ref_win.upper()
    qry_upper = qry_win.upper()
    alphabet = ''.join(sorted(set(ref_upper) | set(qry_upper)))
    matrix = _PARASAIL_MATRIX_CACHE.get(alphabet)
    if matrix is None:
        matrix = _parasail.matrix_create(alphabet, 2, -3)
        _PARASAIL_MATRIX_CACHE[alphabet] = matrix
    result = _parasail.nw_trace_striped_16(
        qry_upper, ref_upper, _SV_WINDOW_GAP_OPEN, _SV_WINDOW_GAP_EXTEND,
        matrix,
    )
    if result.saturated:
        result = _parasail.nw_trace_striped_32(
            qry_upper, ref_upper, _SV_WINDOW_GAP_OPEN,
            _SV_WINDOW_GAP_EXTEND, matrix,
        )
    body = result.cigar.decode
    if isinstance(body, bytes):
        body = body.decode('ascii')
    out = []
    ref_pos = 0
    qry_pos = 0
    for m in re.finditer(r'(\d+)([=XIDM])', body):
        n = int(m.group(1))
        op = m.group(2)
        if op in '=M':
            _tm_add_op(out, n, 'M')
            ref_pos += n
            qry_pos += n
        elif op == 'X':
            _tm_add_op(
                out, n, 'X',
                ref_win[ref_pos:ref_pos + n]
                + qry_win[qry_pos:qry_pos + n],
            )
            ref_pos += n
            qry_pos += n
        elif op == 'I':
            _tm_add_op(out, n, 'I', qry_win[qry_pos:qry_pos + n])
            qry_pos += n
        elif op == 'D':
            _tm_add_op(out, n, 'D', ref_win[ref_pos:ref_pos + n])
            ref_pos += n
    if ref_pos != len(ref_win) or qry_pos != len(qry_win):
        ops, _r0, _q0, _r1, _q1 = _tm_global_insert_align_dp(ref_win, qry_win)
        return [(n, 'M' if op == '=' else op, p) for n, op, p in ops]
    return out


def _split_anchor_op(n: int, op: str, payload: str, left_cols: int) -> Tuple[Tuple[int, str, str], Tuple[int, str, str]]:
    """Split an aligned M/X op into (first left_cols, remainder)."""
    right_cols = n - left_cols
    if op != 'X' or not payload:
        return (left_cols, op, ''), (right_cols, op, '')
    if len(payload) == 2 * n:
        return (
            (left_cols, op, payload[:left_cols] + payload[n:n + left_cols]),
            (right_cols, op, payload[left_cols:n] + payload[n + left_cols:]),
        )
    return (
        (left_cols, op, payload[:left_cols]),
        (right_cols, op, payload[left_cols:]),
    )


def _tm_score_linked_indel_groups(
    ops: Sequence[Tuple[int, str, str]],
) -> List[Tuple[int, int]]:
    """Select aggressive candidate windows for sequence realignment.

    I, D, and X lengths are positive hit scores.  Only intervening exact-match
    M bases cost one point; a mismatch is itself variant evidence and must not
    be charged as background.  There is deliberately no hard distance cutoff.
    A boundary is joined only when both directional sweeps retain a
    non-negative score, matching
    :func:`graphreftovcf.score_merge_by_index`'s two-pass AND rule.  Returned
    pairs are inclusive positions in the variant-op index list.  These pairs
    never define VCF event grouping; graphreftovcf applies its legacy event
    merger after consuming the resulting realigned CIGAR.
    """
    variants = [i for i, (_n, op, _payload) in enumerate(ops) if op in 'IDX']
    if len(variants) < 2:
        return []
    gaps = [
        sum(n for n, op, _payload in ops[a + 1:b] if op == 'M')
        for a, b in zip(variants, variants[1:])
    ]

    def one_pass(order: Sequence[int]) -> Set[int]:
        joined: Set[int] = set()
        score = 0
        active = False
        previous: Optional[int] = None
        for position in order:
            if previous is not None:
                boundary = min(previous, position)
                score -= gaps[boundary]
                if score < 0:
                    score = 0
                    active = False
                if active:
                    joined.add(boundary)
            score += max(1, int(ops[variants[position]][0]))
            active = True
            previous = position
        return joined

    forward = one_pass(list(range(len(variants))))
    reverse = one_pass(list(range(len(variants) - 1, -1, -1)))
    joined = forward & reverse
    groups: List[Tuple[int, int]] = []
    start = 0
    for boundary in range(len(variants) - 1):
        if boundary not in joined:
            if boundary > start:
                groups.append((start, boundary))
            start = boundary + 1
    if len(variants) - 1 > start:
        groups.append((start, len(variants) - 1))
    # Pure mismatch groups do not need SV realignment.  A mismatch can bridge
    # and strengthen any group containing an insertion or deletion.
    return [
        (start, end)
        for start, end in groups
        if any(ops[variants[k]][1] in 'ID' for k in range(start, end + 1))
    ]


def _tm_alignment_offsets(
    ops: Sequence[Tuple[int, str, str]],
) -> Tuple[List[int], List[int], List[int]]:
    """Return alignment-column, reference, and query prefix offsets."""
    columns = [0]
    reference = [0]
    query = [0]
    for n, op, _payload in ops:
        columns.append(columns[-1] + n)
        reference.append(reference[-1] + (n if op in 'MXD' else 0))
        query.append(query[-1] + (n if op in 'MXI' else 0))
    return columns, reference, query


def _tm_axis_offset_at_column(
    ops: Sequence[Tuple[int, str, str]],
    columns: Sequence[int],
    reference: Sequence[int],
    query: Sequence[int],
    column: int,
) -> Tuple[int, int]:
    """Map an alignment-column boundary to reference/query offsets."""
    if column <= 0:
        return 0, 0
    if column >= columns[-1]:
        return reference[-1], query[-1]
    op_index = bisect.bisect_right(columns, column) - 1
    inside = column - columns[op_index]
    _n, op, _payload = ops[op_index]
    return (
        reference[op_index] + (inside if op in 'MXD' else 0),
        query[op_index] + (inside if op in 'MXI' else 0),
    )


def _tm_slice_ops_by_columns(
    ops: Sequence[Tuple[int, str, str]],
    columns: Sequence[int],
    start: int,
    end: int,
) -> List[Tuple[int, str, str]]:
    """Slice pairwise ops by alignment columns without losing payloads."""
    out: List[Tuple[int, str, str]] = []
    if end <= start:
        return out
    op_index = bisect.bisect_right(columns, start) - 1
    while op_index < len(ops) and columns[op_index] < end:
        n, op, payload = ops[op_index]
        piece_start = max(start, columns[op_index]) - columns[op_index]
        piece_end = min(end, columns[op_index + 1]) - columns[op_index]
        take = piece_end - piece_start
        if take > 0:
            if op == 'X':
                piece_payload = _slice_x_payload(
                    payload, n, piece_start, take,
                )
            elif op in 'ID' and payload:
                piece_payload = payload[piece_start:piece_end]
            else:
                piece_payload = ''
            _tm_add_op(out, take, op, piece_payload)
        op_index += 1
    return out


def _tm_column_for_axis_pair(
    ops: Sequence[Tuple[int, str, str]],
    ref_target: int,
    qry_target: int,
) -> Optional[int]:
    """Return the alignment-column boundary at one exact axis-coordinate pair.

    A realignment may move a gap across a requested splice boundary.  In that
    case the old reference/query coordinate pair is absent from the new path
    and there is no safe way to retain the original anchors; return ``None``
    so the caller keeps the old core alignment.
    """
    ref_pos = 0
    qry_pos = 0
    column = 0
    for n, op, _payload in ops:
        if (ref_pos, qry_pos) == (ref_target, qry_target):
            return column
        if op in 'MX':
            ref_delta = ref_target - ref_pos
            qry_delta = qry_target - qry_pos
            if ref_delta == qry_delta and 0 <= ref_delta <= n:
                return column + ref_delta
            ref_pos += n
            qry_pos += n
        elif op == 'D':
            ref_delta = ref_target - ref_pos
            if qry_target == qry_pos and 0 <= ref_delta <= n:
                return column + ref_delta
            ref_pos += n
        elif op == 'I':
            qry_delta = qry_target - qry_pos
            if ref_target == ref_pos and 0 <= qry_delta <= n:
                return column + qry_delta
            qry_pos += n
        column += n
    if (ref_pos, qry_pos) == (ref_target, qry_target):
        return column
    return None


def _tm_project_realigned_core(
    replacement: Sequence[Tuple[int, str, str]],
    core_ref_start: int,
    core_qry_start: int,
    core_ref_end: int,
    core_qry_end: int,
) -> Optional[List[Tuple[int, str, str]]]:
    """Extract only the original non-anchor core from an anchored alignment."""
    replacement = list(replacement)
    replacement_columns, _reference, _query = _tm_alignment_offsets(
        replacement,
    )
    left = _tm_column_for_axis_pair(
        replacement, core_ref_start, core_qry_start,
    )
    right = _tm_column_for_axis_pair(
        replacement, core_ref_end, core_qry_end,
    )
    if left is None or right is None or right < left:
        return None
    core = _tm_slice_ops_by_columns(
        replacement, replacement_columns, left, right,
    )
    ref_span = sum(n for n, op, _payload in core if op in 'MXD')
    qry_span = sum(n for n, op, _payload in core if op in 'MXI')
    if (
        ref_span != core_ref_end - core_ref_start
        or qry_span != core_qry_end - core_qry_start
    ):
        return None
    return core


def _tm_realignment_score(
    ops: Sequence[Tuple[int, str, str]],
) -> int:
    """Score one candidate core before deciding whether to splice it.

    Score = matched bases - mismatched bases - inserted/deleted bases
    - 4 * gap openings.  Consecutive operations on the same gap axis count as
    one opening; switching directly between I and D starts a new gap.
    """
    score = 0
    open_gap = ""
    for n, op, _payload in ops:
        n = int(n)
        if op in 'M=':
            score += n
            open_gap = ""
        elif op == 'X':
            score -= n
            open_gap = ""
        elif op in 'ID':
            score -= n
            if open_gap != op:
                score -= 4
            open_gap = op
        else:
            open_gap = ""
    return score


def _tm_expand_indel_group_window(
    ops: Sequence[Tuple[int, str, str]],
    columns: Sequence[int],
    first_op: int,
    last_op: int,
    anchor: int,
) -> Tuple[int, int]:
    """Expand an indel core by ``anchor`` bases on both sequence axes."""
    left = columns[first_op]
    ref_need = qry_need = anchor
    op_index = first_op - 1
    while op_index >= 0 and (ref_need > 0 or qry_need > 0):
        n, op, _payload = ops[op_index]
        if op in 'MX':
            take = min(n, max(ref_need, qry_need))
        elif op == 'D':
            # Cross the complete single-axis op when the other axis still
            # needs anchor sequence; otherwise a partial boundary is valid.
            take = n if qry_need > 0 else min(n, ref_need)
        elif op == 'I':
            take = n if ref_need > 0 else min(n, qry_need)
        else:
            take = n
        left = columns[op_index + 1] - take
        if op in 'MXD':
            ref_need = max(0, ref_need - take)
        if op in 'MXI':
            qry_need = max(0, qry_need - take)
        if take < n:
            break
        op_index -= 1

    right = columns[last_op + 1]
    ref_need = qry_need = anchor
    op_index = last_op + 1
    while op_index < len(ops) and (ref_need > 0 or qry_need > 0):
        n, op, _payload = ops[op_index]
        if op in 'MX':
            take = min(n, max(ref_need, qry_need))
        elif op == 'D':
            take = n if qry_need > 0 else min(n, ref_need)
        elif op == 'I':
            take = n if ref_need > 0 else min(n, qry_need)
        else:
            take = n
        right = columns[op_index] + take
        if op in 'MXD':
            ref_need = max(0, ref_need - take)
        if op in 'MXI':
            qry_need = max(0, qry_need - take)
        if take < n:
            break
        op_index += 1
    return left, right


def _tm_realign_secondary_sv_window(
    ref_win: str, qry_win: str,
) -> Optional[List[Tuple[int, str, str]]]:
    """Realign a score-linked window with a scalable global backend."""
    def normalize(
        aligned: Sequence[Tuple[int, str, str]],
    ) -> List[Tuple[int, str, str]]:
        return [
            (n, 'M', '') if op in '=M' else (n, op, payload)
            for n, op, payload in aligned
        ]

    if not ref_win:
        return [(len(qry_win), 'I', qry_win)] if qry_win else []
    if not qry_win:
        return [(len(ref_win), 'D', '')] if ref_win else []
    use_affine = (
        max(len(ref_win), len(qry_win)) <= _SV_SECONDARY_AFFINE_MAX_AXIS
        and len(ref_win) * len(qry_win) <= _SV_SECONDARY_AFFINE_MAX_CELLS
    )
    if use_affine and _parasail is not None:
        return normalize(_sv_window_affine_align(ref_win, qry_win))
    try:
        return normalize(minimap2_payload_ops(ref_win, qry_win))
    except (FileNotFoundError, RuntimeError, subprocess.SubprocessError):
        # A tiny window can still be handled safely by the legacy in-process
        # fallback.  For a large repetitive window, preserving the original
        # decomposition is safer than forcing a quadratic/linear-DP result.
        if use_affine and max(len(ref_win), len(qry_win)) <= 1000:
            return normalize(_sv_window_affine_align(ref_win, qry_win))
        return None


def _tm_realign_score_linked_indel_windows(
    ops: List[Tuple[int, str, str]], ref_full: str, qry_full: str,
) -> List[Tuple[int, str, str]]:
    """Second-pass realignment of score-linked insertion regions.

    Unlike the original event merger, this pass has no 100-bp hard break.
    Its groups select realignment regions only.  Each core is expanded on both
    the reference and query by ``min(1000, 10% * max(core axis sizes))``;
    overlapping expanded regions are unioned and aligned exactly once.  The
    anchors remain original; only a strictly better-scoring projected core is
    spliced back.
    """
    variants = [i for i, (_n, op, _payload) in enumerate(ops) if op in 'IDX']
    groups = _tm_score_linked_indel_groups(ops)
    if not groups:
        return ops
    columns, reference, query = _tm_alignment_offsets(ops)
    windows: List[Tuple[int, int, int, int]] = []
    for start, end in groups:
        first_op = variants[start]
        last_op = variants[end]
        core_left = columns[first_op]
        core_right = columns[last_op + 1]
        core_ref = reference[last_op + 1] - reference[first_op]
        core_qry = query[last_op + 1] - query[first_op]
        anchor = int(min(1000, 0.1 * max(core_ref, core_qry)))
        if max(core_ref, core_qry) > 0:
            anchor = max(1, anchor)
        left, right = _tm_expand_indel_group_window(
            ops, columns, first_op, last_op, anchor,
        )
        windows.append((left, right, core_left, core_right))

    # Dynamic anchors from neighboring groups can overlap.  Union them before
    # alignment so no sequence is aligned twice.  Keep each original core
    # separately: the anchors and sequence between cores are never replaced.
    windows.sort()
    unioned: List[List[object]] = []
    for left, right, core_left, core_right in windows:
        if unioned and left <= unioned[-1][1]:
            unioned[-1][1] = max(unioned[-1][1], right)
            unioned[-1][2].append((core_left, core_right))
        else:
            unioned.append([left, right, [(core_left, core_right)]])

    rebuilt: List[Tuple[int, str, str]] = []
    cursor = 0
    for left, right, cores in unioned:
        left = int(left)
        right = int(right)
        r0, q0 = _tm_axis_offset_at_column(
            ops, columns, reference, query, left,
        )
        r1, q1 = _tm_axis_offset_at_column(
            ops, columns, reference, query, right,
        )
        replacement = _tm_realign_secondary_sv_window(
            ref_full[r0:r1], qry_full[q0:q1],
        )
        replacement_ref = sum(
            n for n, op, _payload in (replacement or ()) if op in 'MXD'
        )
        replacement_qry = sum(
            n for n, op, _payload in (replacement or ()) if op in 'MXI'
        )
        replacement_valid = not (
            replacement is None
            or replacement_ref != r1 - r0
            or replacement_qry != q1 - q0
        )
        for core_left, core_right in sorted(cores):
            original_core = _tm_slice_ops_by_columns(
                ops, columns, core_left, core_right,
            )
            for piece in _tm_slice_ops_by_columns(
                ops, columns, cursor, core_left,
            ):
                _tm_add_op(rebuilt, *piece)
            projected = None
            if replacement_valid:
                core_r0, core_q0 = _tm_axis_offset_at_column(
                    ops, columns, reference, query, core_left,
                )
                core_r1, core_q1 = _tm_axis_offset_at_column(
                    ops, columns, reference, query, core_right,
                )
                projected = _tm_project_realigned_core(
                    replacement,
                    core_r0 - r0,
                    core_q0 - q0,
                    core_r1 - r0,
                    core_q1 - q0,
                )
            if (
                projected is None
                or _tm_realignment_score(projected)
                <= _tm_realignment_score(original_core)
            ):
                projected = original_core
            for piece in projected:
                _tm_add_op(rebuilt, *piece)
            cursor = core_right
    for piece in _tm_slice_ops_by_columns(
        ops, columns, cursor, columns[-1],
    ):
        _tm_add_op(rebuilt, *piece)

    in_ref = reference[-1]
    in_qry = query[-1]
    out_ref = sum(n for n, op, _payload in rebuilt if op in 'MXD')
    out_qry = sum(n for n, op, _payload in rebuilt if op in 'MXI')
    if (in_ref, in_qry) != (out_ref, out_qry):
        return ops
    return rebuilt


def _tm_consolidate_sv_ops(ops: List[Tuple[int, str, str]], ref_full: str, qry_full: str) -> List[Tuple[int, str, str]]:
    """Select local and aggressive windows, then realign their sequences.

    Sweep the block's I/D/X variant runs first-to-last with
    score += variant_length - exact_match_gap.  A negative score or an exact
    match gap over 100 columns breaks the chain; the sweep runs in both
    directions and the breaks are unioned.  Each surviving multi-indel
    group's region - extended exactly 50 columns into the flanking
    anchors, splitting anchor ops at the cut - is realigned with affine
    gaps so micro-anchored tandem-repeat indels collapse into one
    canonical decomposition.  Anchors remain unchanged and the core is
    replaced only on a strict alignment-score improvement.  A second pass uses
    unit gap penalty without a hard distance cutoff solely to select broader
    sequence-realignment regions.
    """
    if not _SV_REALIGNMENT_ENABLED:
        return ops
    # Payload realignment helpers emit exact matches as ``=`` while this
    # realignment-window stage historically uses ``M`` internally.  Canonicalize
    # them before all axis and gap accounting; both serialize back to ``=``.
    ops = [
        (n, 'M', '') if op == '=' else (n, op, payload)
        for n, op, payload in ops
    ]
    variants = [i for i, (_n, op, _p) in enumerate(ops) if op in 'IDX']
    if (
        len(variants) < 2
        or not any(ops[index][1] in 'ID' for index in variants)
    ):
        return ops
    # The realignment slices ref_full/qry_full by op-derived offsets, so
    # the captured sequences must consume exactly as the ops do.  If a
    # rare walk desync makes them disagree, skip consolidation rather
    # than realign against shifted bases.
    if (
        len(ref_full) != sum(n for n, op, _p in ops if op in 'MXD')
        or len(qry_full) != sum(n for n, op, _p in ops if op in 'MXI')
    ):
        return ops
    gaps: List[int] = []
    for a, b in zip(variants, variants[1:]):
        gaps.append(sum(n for n, op, _p in ops[a + 1:b] if op == 'M'))
    breaks = set()

    def sweep(order):
        score = 0
        previous = None
        for pos in order:
            length = ops[variants[pos]][0]
            if previous is None:
                score = length
            else:
                gap = gaps[min(previous, pos)]
                if gap > _SV_CONSOLIDATE_MAX_GAP:
                    breaks.add(min(previous, pos))
                    score = length
                else:
                    score += length - gap
                    if score < 0:
                        breaks.add(min(previous, pos))
                        score = length
            previous = pos
    sweep(range(len(variants)))
    sweep(range(len(variants) - 1, -1, -1))
    groups: List[Tuple[int, int]] = []
    start = 0
    for k in range(len(variants) - 1):
        if k in breaks:
            groups.append((start, k))
            start = k + 1
    groups.append((start, len(variants) - 1))
    groups = [
        (s, e) for s, e in groups
        if e > s
        and any(ops[variants[k]][1] in 'ID' for k in range(s, e + 1))
    ]
    if not groups:
        return _tm_realign_score_linked_indel_windows(
            ops, ref_full, qry_full,
        )
    # Pass 1: flank cut points in (op index, column) space.
    cuts_by_op: Dict[int, set] = {}

    def add_cut(op_index: int, col: int) -> None:
        n = ops[op_index][0]
        if 0 < col < n:
            cuts_by_op.setdefault(op_index, set()).add(col)

    windows: List[
        Tuple[
            Tuple[int, int], Tuple[int, int],
            Tuple[int, int], Tuple[int, int],
        ]
    ] = []
    previous_right: Tuple[int, int] = (0, 0)
    for group_index, (s, e) in enumerate(groups):
        first_op = variants[s]
        last_op = variants[e]
        cap_op = (
            variants[groups[group_index + 1][0]]
            if group_index + 1 < len(groups) else len(ops)
        )
        # walk left up to 50 columns over M/X anchors, never crossing
        # the previous group's window
        need = _SV_CONSOLIDATE_FLANK
        left_boundary = (first_op, 0)
        j = first_op - 1
        while need > 0 and j >= previous_right[0] and ops[j][1] in 'MX':
            reserved = previous_right[1] if j == previous_right[0] else 0
            available = ops[j][0] - reserved
            if available <= 0:
                break
            if available <= need and reserved == 0:
                need -= available
                left_boundary = (j, 0)
                j -= 1
            else:
                take = min(need, available)
                left_boundary = (j, ops[j][0] - take)
                add_cut(j, ops[j][0] - take)
                need = 0
        need = _SV_CONSOLIDATE_FLANK
        right_boundary = (last_op + 1, 0)
        j = last_op + 1
        while need > 0 and j < cap_op and ops[j][1] in 'MX':
            if ops[j][0] <= need:
                need -= ops[j][0]
                right_boundary = (j + 1, 0)
                j += 1
            else:
                right_boundary = (j, need)
                add_cut(j, need)
                need = 0
        windows.append((
            left_boundary, right_boundary,
            (first_op, 0), (last_op + 1, 0),
        ))
        previous_right = right_boundary
    # Pass 2: materialize the split op list and remap boundaries.
    split_ops: List[Tuple[int, str, str]] = []
    boundary_map: Dict[Tuple[int, int], int] = {}
    for i, (n, op, payload) in enumerate(ops):
        boundary_map[(i, 0)] = len(split_ops)
        cut_cols = sorted(cuts_by_op.get(i, ()))
        previous = 0
        remainder = (n, op, payload)
        for col in cut_cols:
            first, remainder = _split_anchor_op(
                remainder[0], remainder[1], remainder[2], col - previous,
            )
            split_ops.append(first)
            boundary_map[(i, col)] = len(split_ops)
            previous = col
        split_ops.append(remainder)
    boundary_map[(len(ops), 0)] = len(split_ops)
    ref_off = [0] * (len(split_ops) + 1)
    qry_off = [0] * (len(split_ops) + 1)
    for i, (n, op, _p) in enumerate(split_ops):
        ref_off[i + 1] = ref_off[i] + (n if op in 'MXD' else 0)
        qry_off[i + 1] = qry_off[i] + (n if op in 'MXI' else 0)
    new_ops: List[Tuple[int, str, str]] = []
    cursor = 0
    for (
        left_boundary, right_boundary, core_left_boundary,
        core_right_boundary,
    ) in windows:
        left = boundary_map[left_boundary]
        right = boundary_map[right_boundary]
        core_left = boundary_map[core_left_boundary]
        core_right = boundary_map[core_right_boundary]
        ref_win = ref_full[ref_off[left]:ref_off[right]]
        qry_win = qry_full[qry_off[left]:qry_off[right]]
        if max(len(ref_win), len(qry_win)) > _SV_CONSOLIDATE_MAX_WINDOW:
            continue
        for i in range(cursor, core_left):
            _tm_add_op(new_ops, *split_ops[i])
        replacement = _sv_window_affine_align(ref_win, qry_win)
        original_core = split_ops[core_left:core_right]
        projected = _tm_project_realigned_core(
            replacement,
            ref_off[core_left] - ref_off[left],
            qry_off[core_left] - qry_off[left],
            ref_off[core_right] - ref_off[left],
            qry_off[core_right] - qry_off[left],
        )
        if (
            projected is None
            or _tm_realignment_score(projected)
            <= _tm_realignment_score(original_core)
        ):
            projected = original_core
        for piece in projected:
            _tm_add_op(new_ops, *piece)
        cursor = core_right
    for i in range(cursor, len(split_ops)):
        _tm_add_op(new_ops, *split_ops[i])
    # Safety net: consolidation must be alignment-preserving; if the
    # rebuilt ops consume a different ref/query span than the input,
    # keep the original decomposition rather than corrupt the row.
    in_ref = sum(n for n, op, _p in ops if op in 'MXD')
    in_qry = sum(n for n, op, _p in ops if op in 'MXI')
    out_ref = sum(n for n, op, _p in new_ops if op in 'MXD')
    out_qry = sum(n for n, op, _p in new_ops if op in 'MXI')
    if (in_ref, in_qry) != (out_ref, out_qry):
        return ops
    return _tm_realign_score_linked_indel_windows(
        new_ops, ref_full, qry_full,
    )



def merge_cigars_on_same_template(ref_cigar: str, qry_cigar: str, template_seq: Optional[str]=None, refseq: Optional[str]=None, queryseq: Optional[str]=None, template_offset: int=0) -> str:
    ref = _tm_build_side_track(ref_cigar, template_seq=template_seq, side_seq=refseq, template_offset=template_offset, side_role='ref')
    qry = _tm_build_side_track(qry_cigar, template_seq=template_seq, side_seq=queryseq, template_offset=template_offset, side_role='qry')
    breakpoints = {0, ref.r_end, qry.r_end}
    for iv in ref.intervals:
        breakpoints.add(iv.r0)
        breakpoints.add(iv.r1)
    for iv in qry.intervals:
        breakpoints.add(iv.r0)
        breakpoints.add(iv.r1)
    for p in ref.insertions:
        breakpoints.add(p)
    for p in qry.insertions:
        breakpoints.add(p)
    points = sorted(breakpoints)
    out: List[Tuple[int, str, str]] = []
    ref_parts: List[str] = []
    qry_parts: List[str] = []
    for idx, p in enumerate(points):
        point_ref_ins = ref.insertions.get(p)
        point_qry_ins = qry.insertions.get(p)
        if (
            _SV_REALIGNMENT_ENABLED
            and point_ref_ins is not None
            and point_ref_ins.seq
        ):
            ref_parts.append(point_ref_ins.seq)
        if (
            _SV_REALIGNMENT_ENABLED
            and point_qry_ins is not None
            and point_qry_ins.seq
        ):
            qry_parts.append(point_qry_ins.seq)
        _tm_compare_insertions(point_ref_ins, point_qry_ins, out)
        if idx + 1 >= len(points):
            continue
        a = p
        b = points[idx + 1]
        if b <= a:
            continue
        ref_iv = _tm_find_interval(ref.intervals, a, b)
        qry_iv = _tm_find_interval(qry.intervals, a, b)
        ref_present = ref_iv is not None and ref_iv.op in {'M', 'X'}
        qry_present = qry_iv is not None and qry_iv.op in {'M', 'X'}
        if not ref_present and (not qry_present):
            continue
        if ref_present and (not qry_present):
            ref_seq = _tm_interval_seq(ref_iv, a, b)
            _tm_add_op(out, len(ref_seq), 'D', ref_seq)
            if _SV_REALIGNMENT_ENABLED:
                ref_parts.append(ref_seq)
            continue
        if qry_present and (not ref_present):
            qry_seq = _tm_interval_seq(qry_iv, a, b)
            _tm_add_op(out, len(qry_seq), 'I', qry_seq)
            if _SV_REALIGNMENT_ENABLED:
                qry_parts.append(qry_seq)
            continue
        ref_seq = _tm_interval_seq(ref_iv, a, b)
        qry_seq = _tm_interval_seq(qry_iv, a, b)
        _tm_compare_equal_order(ref_seq, qry_seq, out)
        if _SV_REALIGNMENT_ENABLED:
            ref_parts.append(ref_seq)
            qry_parts.append(qry_seq)
    if _SV_REALIGNMENT_ENABLED:
        out = _tm_consolidate_sv_ops(
            out, ''.join(ref_parts), ''.join(qry_parts),
        )
    return _tm_format_pairwise_ops(out)

def _build_segment_cut_infos(segments: Sequence[GraphicSegment], ordered: Sequence[UnifiedSubsegment], path_sequences: Dict[str, str]) -> Dict[int, SegmentCutInfo]:
    by_original: Dict[int, List[UnifiedSubsegment]] = {}
    for seg in ordered:
        by_original.setdefault(seg.original_index, []).append(seg)
    out: Dict[int, SegmentCutInfo] = {}
    for orig_i, segment in enumerate(segments):
        parts = by_original.get(orig_i, [])
        forward_cuts = sorted({p.coord_start for p in parts} | {p.coord_end for p in parts} | {_stage2_coord_start(p) for p in parts} | {_stage2_coord_end(p) for p in parts})
        cuts = sorted({_segment_forward_to_local_coord(segment, c) for c in forward_cuts})
        local_path_seq = _segment_local_template_seq(segment, 0, segment.path_len, path_sequences)
        block_cigars, qposis = split_cigar_by_rposi(segment.cigar, cuts, reference_seq=local_path_seq, reference_start=0)
        local_bounds = [0] + cuts + [segment.path_len]
        block_intervals: List[Tuple[int, int]] = []
        for a, b in zip(local_bounds, local_bounds[1:]):
            block_intervals.append(_segment_local_interval_to_forward(segment, a, b))
        block_by_interval: Dict[Tuple[int, int], int] = {}
        for idx, (f0, f1) in enumerate(block_intervals):
            if f1 > f0:
                block_by_interval[f0, f1] = idx
        out[orig_i] = SegmentCutInfo(cuts=cuts, block_cigars=block_cigars, qposis=qposis, block_intervals=block_intervals, block_by_interval=block_by_interval)
    return out

def _segment_block_local_interval(segment: GraphicSegment, cut_info: SegmentCutInfo, block_index: int) -> Tuple[int, int]:
    if not cut_info.cuts:
        return (0, segment.path_len)
    if block_index == 0:
        return (0, cut_info.cuts[0])
    if block_index == len(cut_info.block_cigars) - 1:
        return (cut_info.cuts[-1], segment.path_len)
    return (cut_info.cuts[block_index - 1], cut_info.cuts[block_index])

def _segment_block_template_seq(segment: GraphicSegment, cut_info: SegmentCutInfo, block_index: int, path_sequences: Dict[str, str]) -> str:
    a, b = _segment_block_local_interval(segment, cut_info, block_index)
    if b <= a:
        return ''
    return _segment_local_template_seq(segment, a, b, path_sequences)

def _build_ordered_regions(side: str, segments: Sequence[GraphicSegment], ordered: Sequence[UnifiedSubsegment], cut_infos: Dict[int, SegmentCutInfo], path_sequences: Dict[str, str]) -> List[OrderedRegion]:
    """Build the predefined Stage2 chunk stream.

                For each original segment, split_cigar_by_rposi() defines n+1 blocks:
                        block 0      = left leftover
                        block 1..n-1 = valid chunks between unified breakpoint cuts
                        block n      = right leftover

                Standalone leftover/gap regions are not emitted. Instead, block 0 is
                converted to an insertion-style gap and prepended to the first valid block;
                block n is converted and appended to the last valid block. If the owner
                valid block fails LCS, the attached leftovers convert with it into the same
                D/I boundary. If it aligns, the attached leftovers stay in the merge run.
                """
    out: List[OrderedRegion] = []
    by_original: Dict[int, List[Tuple[int, UnifiedSubsegment]]] = {}
    for ordered_idx, seg in enumerate(ordered):
        by_original.setdefault(seg.original_index, []).append((ordered_idx, seg))

    def block_query_len_from_cigar(cigar: str) -> int:
        if not cigar:
            return 0
        return sum((query_consume(tok.op, tok.n) for tok in parse_cigar_ops(cigar)))

    def leftover_as_i_cigar(cigar: str, template_seq: str='') -> str:
        """Normalize an outer leftover block to pure insertion sequence.

                                Outer leftovers are unanchored sequence carried by that side. H and D are
                                dropped. M/= uses template bases, X uses query-side payload, and I keeps
                                its payload. The result is one I operation or an empty string.
                                """
        if not cigar:
            return ''
        seq_chunks: List[str] = []
        tpos = 0
        for tok in parse_cigar_ops(cigar):
            if tok.op in {'H', 'D'}:
                if tok.op == 'H':
                    tpos += tok.n
                continue
            if tok.op == 'I':
                seq_chunks.append(tok.payload if tok.payload and tok.payload != 'n' else 'N' * tok.n)
                continue
            if tok.op == 'X':
                seq_chunks.append(_x_payload_query(tok.payload, tok.n))
                tpos += tok.n
                continue
            if tok.op == '=':
                chunk = template_seq[tpos:tpos + tok.n] if template_seq else ''
                if len(chunk) != tok.n:
                    chunk = chunk.ljust(tok.n, 'N')
                seq_chunks.append(chunk)
                tpos += tok.n
                continue
            raise ValueError(f'unsupported leftover op {tok.op!r}')
        seq = ''.join(seq_chunks)
        return f'{len(seq)}I{seq}' if seq else ''

    def interval_block_indexes(segment: GraphicSegment, cut_info: SegmentCutInfo, start: int, end: int) -> List[int]:
        idxs: List[int] = []
        if end <= start:
            return idxs
        local_start, local_end = _segment_forward_interval_to_local(segment, start, end)
        for idx in range(len(cut_info.cuts) - 1):
            a = cut_info.cuts[idx]
            b = cut_info.cuts[idx + 1]
            if local_start <= a and b <= local_end:
                idxs.append(idx + 1)
        return idxs

    def interval_template_seq(segment: GraphicSegment, start: int, end: int) -> str:
        if end <= start:
            return ''
        return _segment_forward_template_seq(segment, start, end, path_sequences)
    for orig_i, segment in enumerate(segments):
        cut_info = cut_infos[orig_i]
        parts = by_original.get(orig_i, [])
        parts.sort(key=lambda row: (*_segment_forward_interval_to_local(segment, row[1].coord_start, row[1].coord_end), row[0]))
        if not parts:
            whole_chunks: List[str] = []
            whole_len = 0
            for block_i, block in enumerate(cut_info.block_cigars):
                if not block:
                    continue
                tmpl = _segment_block_template_seq(segment, cut_info, block_i, path_sequences)
                if block_i == 0 or block_i == len(cut_info.block_cigars) - 1:
                    block = leftover_as_i_cigar(block, tmpl)
                whole_chunks.append(block)
                whole_len += block_query_len_from_cigar(block)
            whole_cigar = ''.join(whole_chunks)
            if whole_cigar:
                out.append(OrderedRegion(side=side, kind='side_only', original_index=orig_i, path=segment.path, block_index=-1, cigar=whole_cigar, template_seq='', side_seq=None, seq_span=whole_len, anchor_index=-1, ordered_index=None, uses_template=False))
            continue
        valid_rows: List[Tuple[int, UnifiedSubsegment, int, int, int, str, str]] = []
        for ordered_idx, part in parts:
            core_block_i = cut_info.block_by_interval.get((part.coord_start, part.coord_end))
            if core_block_i is None:
                continue
            stage2_start = _stage2_coord_start(part)
            stage2_end = _stage2_coord_end(part)
            block_ids = interval_block_indexes(segment, cut_info, stage2_start, stage2_end)
            if not block_ids:
                continue
            cigar = ''.join((cut_info.block_cigars[i] for i in block_ids if cut_info.block_cigars[i]))
            if not cigar:
                continue
            template_seq = interval_template_seq(segment, stage2_start, stage2_end)
            valid_rows.append((ordered_idx, part, core_block_i, block_ids[0], block_ids[-1], cigar, template_seq))
        if not valid_rows:
            continue
        first_valid_i = valid_rows[0][3]
        last_valid_i = valid_rows[-1][4]

        def attach_blocks_as_i(block_ids: Sequence[int]) -> str:
            chunks: List[str] = []
            for block_i in block_ids:
                if 0 <= block_i < len(cut_info.block_cigars):
                    chunks.append(leftover_as_i_cigar(cut_info.block_cigars[block_i], _segment_block_template_seq(segment, cut_info, block_i, path_sequences)))
            return ''.join(chunks)
        left_attach = attach_blocks_as_i(range(0, first_valid_i))
        right_attach = attach_blocks_as_i(range(last_valid_i + 1, len(cut_info.block_cigars)))
        for ordered_idx, part, core_block_i, first_block_i, last_block_i, region_cigar, template_seq in valid_rows:
            cigar = region_cigar
            if first_block_i == first_valid_i and left_attach:
                cigar = left_attach + cigar
            if last_block_i == last_valid_i and right_attach:
                cigar = cigar + right_attach
            out.append(OrderedRegion(side=side, kind='valid', original_index=orig_i, path=segment.path, block_index=core_block_i, cigar=cigar, template_seq=template_seq, side_seq=None, seq_span=block_query_len_from_cigar(cigar), anchor_index=ordered_idx, ordered_index=ordered_idx, uses_template=True, coord_start=_stage2_coord_start(part), coord_end=_stage2_coord_end(part)))
    return out

def build_stage2(stage1: Stage1Result, path_sequences: Dict[str, str]) -> Stage2Result:
    runs = build_matched_runs(stage1)
    ref_cut_infos = _build_segment_cut_infos(stage1.ref_segments, stage1.ref_ordered, path_sequences)
    qry_cut_infos = _build_segment_cut_infos(stage1.qry_segments, stage1.qry_ordered, path_sequences)
    ref_regions = _build_ordered_regions('ref', stage1.ref_segments, stage1.ref_ordered, ref_cut_infos, path_sequences)
    qry_regions = _build_ordered_regions('qry', stage1.qry_segments, stage1.qry_ordered, qry_cut_infos, path_sequences)
    ref_aligned_to_qry = {ref_idx: qry_idx for ref_idx, qry_idx in stage1.lcs_pairs}
    qry_region_pos_by_valid: Dict[int, int] = {region.ordered_index: idx for idx, region in enumerate(qry_regions) if region.kind == 'valid' and region.ordered_index is not None}
    units: List[RawAssemblyUnit] = []

    def _template_query_len(cigar: str) -> int:
        return sum((n if op in {'M', 'X', 'I'} else 0 for n, op, _payload in _tm_parse_template_cigar(cigar)))

    def _sequence_from_region_cigar(region: OrderedRegion) -> str:
        """Return the side/sample sequence carried by a Stage-2 region."""
        out: List[str] = []
        qpos = 0
        tpos = 0
        for tok in parse_cigar_ops(region.cigar):
            if tok.op == 'I':
                if tok.payload and tok.payload != 'n':
                    out.append(tok.payload)
                elif region.side_seq is not None:
                    out.append(region.side_seq[qpos:qpos + tok.n])
                else:
                    out.append('N' * tok.n)
                qpos += tok.n
            elif tok.op == 'X':
                if tok.payload:
                    seq = _x_payload_query(tok.payload, tok.n)
                elif region.side_seq is not None:
                    seq = region.side_seq[qpos:qpos + tok.n]
                else:
                    seq = 'N' * tok.n
                out.append(seq if len(seq) == tok.n else seq.ljust(tok.n, 'N'))
                qpos += tok.n
                tpos += tok.n
            elif tok.op == '=':
                if region.side_seq is not None:
                    seq = region.side_seq[qpos:qpos + tok.n]
                else:
                    seq = region.template_seq[tpos:tpos + tok.n]
                out.append(seq if len(seq) == tok.n else seq.ljust(tok.n, 'N'))
                qpos += tok.n
                tpos += tok.n
            elif tok.op == 'D':
                tpos += tok.n
            elif tok.op == 'H':
                continue
            else:
                raise ValueError(f'unsupported region op {tok.op!r}')
        return ''.join(out)

    def _seq_as_template_insertion(seq: str) -> str:
        return f'{len(seq)}I{seq}' if seq else ''

    def _region_unit_cigar(region: OrderedRegion) -> Tuple[str, str, Optional[str]]:
        """Return (cigar, template_seq, side_seq) for normal merge-run units.

                                After LCS, path/template identity should not break the merge run. Valid
                                regions keep their template CIGAR/template sequence. Gap leftovers from
                                split_cigar_by_rposi() are not true ref-only/query-only evidence; they
                                are represented as insertion tracks on an empty local template so the
                                ref/qry leftovers can compare with each other inside the same run.
                                """
        if region.kind == 'gap' or not region.uses_template:
            return (_seq_as_template_insertion(_sequence_from_region_cigar(region)), '', None)
        return (region.cigar, region.template_seq, region.side_seq)

    def _region_side_length(region: OrderedRegion) -> int:
        if region.kind == 'gap' or not region.uses_template:
            return len(_sequence_from_region_cigar(region))
        return _template_query_len(region.cigar)

    def add_run_region(region: OrderedRegion) -> None:
        cigar, template_seq, side_seq = _region_unit_cigar(region)
        if not cigar:
            return
        label = f'{region.side}.{region.kind}.orig{region.original_index}.block{region.block_index}'
        if region.side == 'ref':
            units.append(RawAssemblyUnit(cigar, '', template_seq, ref_side_seq=side_seq))
        else:
            units.append(RawAssemblyUnit('', cigar, template_seq, qry_side_seq=side_seq))

    def add_run_pair(ref_region: OrderedRegion, qry_region: OrderedRegion) -> None:
        ref_cigar, ref_template, ref_side_seq = _region_unit_cigar(ref_region)
        qry_cigar, qry_template, qry_side_seq = _region_unit_cigar(qry_region)
        if not ref_cigar and (not qry_cigar):
            return

        def region_direction(region: OrderedRegion) -> str:
            segments = stage1.ref_segments if region.side == 'ref' else stage1.qry_segments
            if 0 <= region.original_index < len(segments):
                return segments[region.original_index].direction
            return '>'

        def pad_to_union(region: OrderedRegion, cigar: str, union_start: int, union_end: int) -> str:
            if not cigar or not region.uses_template:
                return cigar
            if region_direction(region) == '<':
                prefix = max(0, union_end - region.coord_end)
                suffix = max(0, region.coord_start - union_start)
            else:
                prefix = max(0, region.coord_start - union_start)
                suffix = max(0, union_end - region.coord_end)
            if prefix:
                cigar = f'{prefix}D' + cigar
            if suffix:
                cigar = cigar + f'{suffix}D'
            return cigar
        ref_dir = region_direction(ref_region)
        qry_dir = region_direction(qry_region)
        template_seq = qry_template if qry_template else ref_template
        if ref_region.uses_template and qry_region.uses_template and (path_key(ref_region.path) == path_key(qry_region.path)) and (ref_region.coord_end > ref_region.coord_start) and (qry_region.coord_end > qry_region.coord_start):
            union_start = min(ref_region.coord_start, qry_region.coord_start)
            union_end = max(ref_region.coord_end, qry_region.coord_end)
            path_seq = path_sequences.get(ref_region.path) or path_sequences.get(path_key(ref_region.path))
            if path_seq is not None and union_start >= 0 and (union_end <= len(path_seq)):
                template_seq = _oriented_path_seq(path_seq, union_start, union_end, qry_dir)
                if union_end > union_start and (union_start != ref_region.coord_start or union_end != ref_region.coord_end or union_start != qry_region.coord_start or (union_end != qry_region.coord_end)):
                    ref_cigar = pad_to_union(ref_region, ref_cigar, union_start, union_end)
                    qry_cigar = pad_to_union(qry_region, qry_cigar, union_start, union_end)
                if ref_dir != qry_dir:
                    ref_cigar = _reverse_template_cigar(ref_cigar)
        elif ref_region.uses_template and qry_region.uses_template and (ref_dir != qry_dir):
            ref_cigar = _reverse_template_cigar(ref_cigar)
            if ref_template and (not qry_template):
                template_seq = revcomp(ref_template)
        units.append(RawAssemblyUnit(ref_cigar, qry_cigar, template_seq, ref_side_seq=ref_side_seq, qry_side_seq=qry_side_seq))

    def add_side_only_boundary(region: OrderedRegion, extra_regions: Sequence[OrderedRegion]=()) -> None:
        n = _region_side_length(region) + sum((_region_side_length(extra) for extra in extra_regions))
        if n <= 0:
            return
        units.append(RawAssemblyUnit('', '', '', break_side=region.side, break_len=n))

    def emit_unaligned_regions(regions: Sequence[OrderedRegion]) -> None:
        i = 0
        while i < len(regions):
            group = [regions[i]]
            j = i + 1
            while j < len(regions) and regions[j].side == regions[i].side and (regions[j].original_index == regions[i].original_index):
                group.append(regions[j])
                j += 1
            boundary_regions = [r for r in group if r.kind in {'side_only', 'valid'}]
            if boundary_regions:
                representative = boundary_regions[0]
                add_side_only_boundary(representative, [r for r in group if r is not representative])
            else:
                for region in group:
                    add_run_region(region)
            i = j

    def known_unaligned_sequence(regions: Sequence[OrderedRegion]) -> Optional[str]:
        """Return a complete A/C/G/T sequence for one unmatched side."""
        chunks: List[str] = []
        expected = 0
        for region in regions:
            expected += _region_side_length(region)
            chunks.append(_sequence_from_region_cigar(region))
        sequence = ''.join(chunks)
        if not sequence or len(sequence) != expected:
            return None
        if any(base.upper() not in {'A', 'C', 'G', 'T'} for base in sequence):
            return None
        return sequence

    def emit_unaligned_window(
        ref_window: Sequence[OrderedRegion],
        qry_window: Sequence[OrderedRegion],
    ) -> None:
        """Sequence-align a known two-sided window between LCS anchors.

        Previously the two sides were emitted independently as ``Dn`` and
        ``In`` boundaries.  That bypassed the existing insertion-payload
        realigner, so different graph paths carrying identical bases became a
        false equal-length substitution.  Encoding both sequences as
        insertions on one empty local template routes them through the existing
        insertion-vs-insertion alignment.  One-sided or unknown-base windows
        retain the original boundary behavior.
        """
        if ref_window and qry_window:
            ref_sequence = known_unaligned_sequence(ref_window)
            qry_sequence = known_unaligned_sequence(qry_window)
            if ref_sequence is not None and qry_sequence is not None:
                units.append(RawAssemblyUnit(
                    _seq_as_template_insertion(ref_sequence),
                    _seq_as_template_insertion(qry_sequence),
                    '',
                ))
                return
        emit_unaligned_regions(ref_window)
        emit_unaligned_regions(qry_window)

    ref_pos = 0
    qry_pos = 0
    while True:
        next_ref_aligned_pos: Optional[int] = None
        for idx in range(ref_pos, len(ref_regions)):
            region = ref_regions[idx]
            if region.kind == 'valid' and region.ordered_index in ref_aligned_to_qry:
                next_ref_aligned_pos = idx
                break
        if next_ref_aligned_pos is None:
            emit_unaligned_window(
                ref_regions[ref_pos:], qry_regions[qry_pos:],
            )
            break
        ref_region = ref_regions[next_ref_aligned_pos]
        assert ref_region.ordered_index is not None
        qry_valid_idx = ref_aligned_to_qry[ref_region.ordered_index]
        qry_aligned_pos = qry_region_pos_by_valid[qry_valid_idx]
        emit_unaligned_window(
            ref_regions[ref_pos:next_ref_aligned_pos],
            qry_regions[qry_pos:qry_aligned_pos],
        )
        add_run_pair(ref_region, qry_regions[qry_aligned_pos])
        ref_pos = next_ref_aligned_pos + 1
        qry_pos = qry_aligned_pos + 1

    def _unit_is_boundary(unit: RawAssemblyUnit) -> bool:
        return unit.break_side in {'ref', 'qry'}

    def _concat_side_seq(run: Sequence[RawAssemblyUnit], side: str) -> Optional[str]:
        chunks: List[str] = []
        for unit in run:
            cigar = unit.ref_cigar if side == 'ref' else unit.qry_cigar
            seq = unit.ref_side_seq if side == 'ref' else unit.qry_side_seq
            qlen = _template_query_len(cigar)
            if qlen == 0:
                continue
            if seq is None or len(seq) != qlen:
                return None
            chunks.append(seq)
        return ''.join(chunks) if chunks else None
    merge_run_index = 0

    def _merge_unit_run(run: Sequence[RawAssemblyUnit]) -> str:
        nonlocal merge_run_index
        if not run:
            return ''
        ref_chain = ''.join((unit.ref_cigar for unit in run))
        qry_chain = ''.join((unit.qry_cigar for unit in run))
        template_chain = ''.join((unit.template_seq for unit in run))
        refseq = _concat_side_seq(run, 'ref')
        queryseq = _concat_side_seq(run, 'qry')
        merged = merge_cigars_on_same_template(ref_chain, qry_chain, template_seq=template_chain, refseq=refseq, queryseq=queryseq)
        merge_run_index += 1
        return merged
    pairwise_parts: List[str] = []
    run_units: List[RawAssemblyUnit] = []

    def flush_run() -> None:
        if not run_units:
            return
        pairwise_parts.append(_merge_unit_run(run_units))
        run_units.clear()
    for unit in units:
        if _unit_is_boundary(unit):
            flush_run()
            if unit.break_side == 'ref':
                boundary_cigar = f'{unit.break_len}Dn' if unit.break_len > 0 else ''
            else:
                boundary_cigar = f'{unit.break_len}In' if unit.break_len > 0 else ''
            pairwise_parts.append(boundary_cigar)
            continue
        run_units.append(unit)
    flush_run()
    pairwise_cigar = ''.join((part for part in pairwise_parts if part))
    main_piece = PairwisePiece(kind='main', path=stage1.ref_segments[0].path if stage1.ref_segments else '', ref_original_index=0 if stage1.ref_segments else None, qry_original_index=0 if stage1.qry_segments else None, coord_start=0, coord_end=pairwise_ref_span(pairwise_cigar), ref_path_start=0, ref_path_end=pairwise_ref_span(pairwise_cigar), cigar=pairwise_cigar)
    return Stage2Result(matched_runs=runs, pieces=[main_piece], pairwise_cigar=pairwise_cigar)

def build_stage3_polished(stage2: Stage2Result) -> Stage3PolishResult:
    kept: List[PairwisePiece] = []
    for piece in stage2.pieces:
        if piece.kind in {'main', 'match'}:
            kept.append(polish_match_piece(piece))
        else:
            kept.append(piece)
    return Stage3PolishResult(pieces=kept, pairwise_cigar=''.join((piece.cigar for piece in kept)))

def project_linear_piece(piece: PairwisePiece, path_coords: Dict[str, Coord]) -> Optional[LinearPiece]:
    coord = path_coords.get(piece.path) or path_coords.get(path_key(piece.path))
    if coord is None:
        return None
    span = coord.end - coord.start
    if span <= 0:
        return None
    if piece.ref_path_start < 0:
        return None
    ref_span = pairwise_ref_span(piece.cigar)
    if ref_span <= 0:
        return None
    if coord.strand == '+':
        g0 = coord.start + piece.ref_path_start
        g1 = g0 + ref_span
    else:
        g1 = coord.end - piece.ref_path_start
        g0 = g1 - ref_span
    if g0 < 0:
        return None
    return LinearPiece(piece=piece, chrom=coord.chrom, strand=coord.strand, genome_start=g0, genome_end=g1, cigar=piece.cigar)

def project_linear_piece_on_reference_backbone(piece: PairwisePiece, backbone: Coord, ref_offset_start: int) -> Optional[LinearPiece]:
    ref_span = pairwise_ref_span(piece.cigar)
    if ref_span <= 0:
        return None
    if backbone.strand == '+':
        g0 = backbone.start + ref_offset_start
        g1 = g0 + ref_span
    else:
        g1 = backbone.end - ref_offset_start
        g0 = g1 - ref_span
    if g0 < 0:
        return None
    return LinearPiece(piece=piece, chrom=backbone.chrom, strand=backbone.strand, genome_start=g0, genome_end=g1, cigar=piece.cigar)

def _zero_ref_linear_piece(cigar: str) -> LinearPiece:
    dummy_piece = PairwisePiece(kind='plain_ins', path='', ref_original_index=None, qry_original_index=None, coord_start=0, coord_end=0, ref_path_start=0, ref_path_end=0, cigar=cigar)
    return LinearPiece(piece=dummy_piece, chrom='', strand='+', genome_start=0, genome_end=0, cigar=cigar)

def _insertion_cigar_from_payload(n: int, payload: Optional[str]) -> str:
    if payload is not None and len(payload) == n:
        return f'{n}I{payload}'
    return f'{n}I'

def _payload_fragment(payload: Optional[str], insert_q0: int, frag_q0: int, frag_q1: int) -> Optional[str]:
    if payload is None:
        return None
    rel0 = frag_q0 - insert_q0
    rel1 = frag_q1 - insert_q0
    if rel0 < 0 or rel1 < rel0 or rel1 > len(payload):
        return None
    return payload[rel0:rel1]

def _query_sequence_fragment(query_sequence: Optional[str], q0: int, q1: int) -> Optional[str]:
    """Return query FASTA bases for a row-local query interval, if valid.

                Large-insertion encoding should use the real query sequence from -q as
                the source of inserted bases.  This avoids recovering bases from graph
                CIGAR coordinates plus -g path sequence.
                """
    if query_sequence is None:
        return None
    q0 = int(q0)
    q1 = int(q1)
    if q0 < 0 or q1 < q0 or q1 > len(query_sequence):
        return None
    return query_sequence[q0:q1]

def _lookup_query_sequence(query_name: Optional[str], record_sequences: Optional[Dict[str, str]]) -> Optional[str]:
    if not query_name or record_sequences is None:
        return None
    return record_sequences.get(query_name) or record_sequences.get(strip_match_name(query_name))

def _pad_or_trim(seq: str, n: int) -> str:
    if len(seq) < n:
        return seq + 'N' * (n - len(seq))
    return seq[:n]

def _template_path_sequence_from_query_gcigar(qry_gcigar: str, query_fragment: Optional[str]) -> Optional[str]:
    """Build a synthetic local template from -q query bases.

                The older code recovered this template by slicing graph_sequences (-g)
                using graph-CIGAR projected coordinates.  That is fragile for sliced
                graph paths and can request impossible ranges.  This function uses the
                actual query FASTA fragment from -q instead.
                """
    if query_fragment is None:
        return None
    segs = parse_graphic_segments(qry_gcigar, 'query_template_from_q')
    if not segs:
        return None
    qpos = 0
    out: List[str] = []
    for tok in _strip_terminal_h(segs[0].ops):
        if tok.op == '=':
            out.append(_pad_or_trim(query_fragment[qpos:qpos + tok.n], tok.n))
            qpos += tok.n
        elif tok.op == 'X':
            if tok.payload:
                out.append(_pad_or_trim(_x_payload_query(tok.payload, tok.n), tok.n))
            else:
                out.append(_pad_or_trim(query_fragment[qpos:qpos + tok.n], tok.n))
            qpos += tok.n
        elif tok.op == 'I':
            qpos += tok.n
        elif tok.op == 'D':
            out.append('N' * tok.n)
    return ''.join(out)

def _pairwise_ref_consume_op(tok: CigarOp) -> int:
    return tok.n if tok.op in {'=', 'X', 'D'} else 0

def _linear_subpiece_from_ref_offsets(linear: LinearPiece, cigar: str, ref_offset_start: int, ref_offset_end: int) -> Optional[LinearPiece]:
    if not cigar:
        return None
    ref_span = pairwise_ref_span(cigar)
    if ref_span == 0:
        return _zero_ref_linear_piece(cigar)
    if ref_offset_end < ref_offset_start:
        return None
    if linear.strand == '+':
        genome_start = linear.genome_start + ref_offset_start
        genome_end = linear.genome_start + ref_offset_end
    else:
        genome_start = linear.genome_end - ref_offset_end
        genome_end = linear.genome_end - ref_offset_start
    src = linear.piece
    new_piece = PairwisePiece(kind=src.kind, path=src.path, ref_original_index=src.ref_original_index, qry_original_index=src.qry_original_index, coord_start=src.coord_start, coord_end=src.coord_end, ref_path_start=src.ref_path_start, ref_path_end=src.ref_path_end, cigar=cigar)
    return LinearPiece(piece=new_piece, chrom=linear.chrom, strand=linear.strand, genome_start=genome_start, genome_end=genome_end, cigar=cigar)

def encode_query_interval_by_reference(insert_q0: int, insert_q1: int, insert_payload: Optional[str], query_sequence: Optional[str], stage1: Stage1Result, graph_sequences: Dict[str, str], graph_mappings: Dict[str, str], ref_reader: Optional[FastaRegionReader], query_coverages: Optional[Sequence[QueryCoverage]]=None) -> List[LinearPiece]:
    """Replace one large query insertion interval by lifted reference pieces.

                insert_q0/insert_q1 are query-local coordinates on the original graph CIGAR
                row, not coordinates local to a single path segment.  This is the direct
                implementation of the AllinsertsQcoordi x AllpathsQcoordi overlap step.
                """
    insert_q0 = int(insert_q0)
    insert_q1 = int(insert_q1)
    if insert_q1 <= insert_q0:
        return []

    def fallback_piece(a: int, b: int) -> LinearPiece:
        frag = _payload_fragment(insert_payload, insert_q0, a, b)
        if frag is None:
            frag = _query_sequence_fragment(query_sequence, a, b)
        return _zero_ref_linear_piece(_insertion_cigar_from_payload(b - a, frag))
    coverages = query_coverages if query_coverages is not None else build_query_coverages(stage1.qry_segments, graph_sequences, graph_mappings)
    overlaps: List[Tuple[int, int, QueryCoverage]] = []
    for cov in coverages:
        ov0 = max(insert_q0, cov.q0)
        ov1 = min(insert_q1, cov.q1)
        if ov1 > ov0:
            overlaps.append((ov0, ov1, cov))
    overlaps.sort(key=lambda row: (row[0], row[1], row[2].segment_index))
    if not overlaps:
        return [fallback_piece(insert_q0, insert_q1)]
    work_items: List[Tuple[str, object]] = []
    cursor = insert_q0
    for ov0, ov1, cov in overlaps:
        if ov0 > cursor:
            work_items.append(('plain', fallback_piece(cursor, ov0)))
        rel_q0 = ov0 - cov.q0
        rel_q1 = ov1 - cov.q0
        cov_seg = parse_graphic_segments(cov.gcigar, 'qry_cov_slice')[0]
        qry_slice_gcigar = slice_query_segment_by_query(cov_seg, rel_q0, rel_q1)
        if not qry_slice_gcigar:
            work_items.append(('plain', fallback_piece(ov0, ov1)))
            cursor = ov1
            continue
        qry_slice_seg = parse_graphic_segments(qry_slice_gcigar, 'qry_slice')[0]
        mapping = cov.mapping_gcigar
        if not mapping:
            work_items.append(('plain', fallback_piece(ov0, ov1)))
            cursor = ov1
            continue
        mapped_q0 = qry_slice_seg.start - cov.path_q0
        mapped_q1 = qry_slice_seg.end - cov.path_q0
        mapped_slice = slice_graphic_mapping_by_query(mapping, mapped_q0, mapped_q1, ref_reader)
        if not mapped_slice:
            work_items.append(('plain', fallback_piece(ov0, ov1)))
            cursor = ov1
            continue
        inv_gcigar, ref_meta = invert_path_to_ref_slice_for_compare(mapped_slice, qry_slice_seg.path, qry_slice_seg.start, qry_slice_seg.end, qry_slice_seg.path_len, qry_slice_seg.direction, ref_reader)
        if not inv_gcigar or ref_meta is None:
            work_items.append(('plain', fallback_piece(ov0, ov1)))
            cursor = ov1
            continue
        ref_direction, ref_chrom, ref_start, ref_end, _chrom_len = ref_meta
        strand = '+' if qry_slice_seg.direction == ref_direction else '-'
        norm_q_gcigar, _norm_q_span = _normalize_gcigar_for_anneal(qry_slice_gcigar, qry_slice_seg.path)
        norm_r_gcigar, _norm_r_span = _normalize_gcigar_for_anneal(inv_gcigar, qry_slice_seg.path)
        query_fragment = _query_sequence_fragment(query_sequence, ov0, ov1)
        if query_fragment is None:
            query_fragment = _payload_fragment(insert_payload, insert_q0, ov0, ov1)
        if query_fragment is None:
            work_items.append(('plain', fallback_piece(ov0, ov1)))
            cursor = ov1
            continue
        template_query_fragment = revcomp(query_fragment) if qry_slice_seg.direction == '<' else query_fragment
        template_path_seq = _template_path_sequence_from_query_gcigar(norm_q_gcigar, template_query_fragment)
        if template_path_seq is None:
            work_items.append(('plain', fallback_piece(ov0, ov1)))
            cursor = ov1
            continue
        work_items.append(('template', AnnealTemplate(q0=ov0, q1=ov1, chrom=ref_chrom, strand=strand, ref_start=ref_start, ref_end=ref_end, qry_gcigar=norm_q_gcigar, ref_gcigar=norm_r_gcigar, path_seq=template_path_seq)))
        cursor = ov1
    if cursor < insert_q1:
        work_items.append(('plain', fallback_piece(cursor, insert_q1)))
    out: List[LinearPiece] = []
    pending_templates: List[AnnealTemplate] = []

    def flush_templates() -> None:
        nonlocal pending_templates
        if not pending_templates:
            return
        for template in pending_templates:
            q_seg = parse_graphic_segments(template.qry_gcigar, 'anneal_compare_q')[0]
            r_seg = parse_graphic_segments(template.ref_gcigar, 'anneal_compare_r')[0]
            path_sequences = {q_seg.path: template.path_seq}
            enc_cigar, _ref_path_start, _ref_path_end = compare_sliced_gcigars(template.ref_gcigar, template.qry_gcigar, path_sequences)
            anchored = _anchor_annealed_linear_piece(template, enc_cigar)
            if anchored is not None:
                out.append(anchored)
            else:
                out.append(fallback_piece(template.q0, template.q1))
        pending_templates = []
    for kind, value in work_items:
        if kind == 'template':
            pending_templates.append(value)
            continue
        flush_templates()
        out.append(value)
    flush_templates()
    return out

def _promote_same_backbone_encoded_pieces_to_main(pieces: List[LinearPiece], host_linear_piece: Optional[LinearPiece]) -> None:
    """Promote only an exact, non-overlapping continuation of the main path.

    A recovered query piece normally stays ``encoded_qry``.  The exception is
    a piece immediately adjacent to a main piece on the query axis whose
    reference interval touches that main piece on the same chromosome and
    strand.  Such a piece is not an insertion: it is a split serialization of
    one continuous reference walk.  Leaving it encoded emits a second named
    graph-CIGAR segment, which downstream code correctly interprets as an
    alternative traversal and can therefore turn into a huge false insertion.

    Reference overlap is deliberately forbidden.  This keeps secondary
    tandem/duplicate matches from consuming reference bases already owned by
    the primary alignment.
    """
    if host_linear_piece is None or not pieces:
        return

    def overlaps_host(candidate: LinearPiece) -> bool:
        if candidate.chrom != host_linear_piece.chrom:
            return False
        return (
            max(candidate.genome_start, host_linear_piece.genome_start)
            < min(candidate.genome_end, host_linear_piece.genome_end)
        )

    def is_supported(candidate: LinearPiece) -> bool:
        return (
            candidate.piece.kind == 'encoded_qry'
            and candidate.chrom == host_linear_piece.chrom
            and candidate.strand == host_linear_piece.strand
            and pairwise_ref_span(candidate.cigar) > 0
            and _piece_match_size_from_cigar(candidate.cigar) > 0
            and not overlaps_host(candidate)
        )

    def follows(left: LinearPiece, right: LinearPiece) -> bool:
        if (
            left.piece.kind != 'main'
            or left.chrom != right.chrom
            or left.strand != right.strand
        ):
            return False
        if left.strand == '+':
            return left.genome_end == right.genome_start
        return left.genome_start == right.genome_end

    # Repeat so a run of recovered boundary pieces can extend the main path
    # one exact interval at a time, without jumping over a gap or insertion.
    changed = True
    while changed:
        changed = False
        for index, candidate in enumerate(pieces):
            if not is_supported(candidate):
                continue
            left = pieces[index - 1] if index > 0 else None
            right = pieces[index + 1] if index + 1 < len(pieces) else None
            parent = None
            if left is not None and follows(left, candidate):
                parent = left
            elif right is not None and follows(candidate, right):
                parent = right
            if parent is None:
                continue
            candidate.piece.kind = 'main'
            candidate.main_split_parent = parent.main_split_parent
            changed = True

def encode_large_insertions_in_piece(piece: PairwisePiece, piece_q0: int, linear_piece: Optional[LinearPiece], query_sequence: Optional[str], stage1: Stage1Result, graph_sequences: Dict[str, str], graph_mappings: Dict[str, str], ref_reader: Optional[FastaRegionReader], min_encode_query_span: int, local_duplicate_pieces: Optional[Sequence[LocalDuplicateEncodedPiece]]=None, query_coverages: Optional[Sequence[QueryCoverage]]=None) -> List[LinearPiece]:
    """Split one pairwise piece and encode every I op longer than threshold."""
    ops = parse_cigar_ops(piece.cigar)
    if not ops:
        return []
    has_large_insert = any((tok.op == 'I' and tok.n > min_encode_query_span for tok in ops))
    if not has_large_insert:
        if linear_piece is not None:
            return [linear_piece]
        return [_zero_ref_linear_piece(piece.cigar)]
    out: List[LinearPiece] = []
    current_ops: List[CigarOp] = []
    current_ref_start = 0
    qpos = int(piece_q0)
    ref_offset = 0
    main_split_parent = -1
    if linear_piece is not None and piece.kind == 'main':
        main_split_parent = linear_piece.main_split_parent if linear_piece.main_split_parent >= 0 else id(linear_piece)

    def flush_current() -> None:
        nonlocal current_ops, current_ref_start
        if not current_ops:
            return
        cigar = ''.join((_op_to_cigar(tok) for tok in current_ops))
        if linear_piece is None:
            out.append(_zero_ref_linear_piece(cigar))
        else:
            sub = _linear_subpiece_from_ref_offsets(linear_piece, cigar, current_ref_start, ref_offset)
            if sub is not None:
                if sub.piece.kind == 'main' and main_split_parent >= 0:
                    sub.main_split_parent = main_split_parent
                out.append(sub)
        current_ops = []
        current_ref_start = ref_offset
    for tok in ops:
        if tok.op == 'I' and tok.n > min_encode_query_span:
            flush_current()
            payload = tok.payload if len(tok.payload) == tok.n else None
            out.extend(encode_query_interval_with_local_duplicates(qpos, qpos + tok.n, payload, query_sequence, stage1, graph_sequences, graph_mappings, ref_reader, local_duplicate_pieces, query_coverages))
            qpos += tok.n
            current_ref_start = ref_offset
            continue
        if not current_ops:
            current_ref_start = ref_offset
        current_ops.append(tok)
        qpos += query_consume(tok.op, tok.n)
        ref_offset += _pairwise_ref_consume_op(tok)
    flush_current()
    local_duplicate_ids = {id(item.linear_piece) for item in local_duplicate_pieces or []}
    _promote_same_backbone_encoded_pieces_to_main(out, linear_piece)
    for encoded_piece in out:
        if (
            id(encoded_piece) in local_duplicate_ids
            and encoded_piece.piece.kind != 'main'
        ):
            encoded_piece.piece.kind = 'encoded_qry'
    return out

def compare_sliced_gcigars(ref_gcigar: str, qry_gcigar: str, path_sequences: Dict[str, str]) -> Tuple[str, int, int]:
    local_stage1 = build_stage1(ref_gcigar, qry_gcigar, 'ref_slice', 'qry_slice', max_merge_gap=1)
    local_stage2 = build_stage2(local_stage1, path_sequences)
    cigar = ''.join((piece.cigar for piece in local_stage2.pieces))
    ref_spans = [(piece.ref_path_start, piece.ref_path_end) for piece in local_stage2.pieces if piece.ref_path_end > piece.ref_path_start]
    if ref_spans:
        ref_path_start = min((a for a, _b in ref_spans))
        ref_path_end = max((b for _a, b in ref_spans))
    else:
        ref_seg = parse_graphic_segments(ref_gcigar, 'ref_slice')[0]
        ref_path_start = ref_seg.start
        ref_path_end = ref_seg.start
    return (cigar, ref_path_start, ref_path_end)

def _normalize_gcigar_for_anneal(gcigar: str, path_name: str) -> Tuple[str, int]:
    seg = parse_graphic_segments(gcigar, 'anneal_norm')[0]
    body_ops = _canonical_body_ops(seg)
    span = sum((ref_consume_pair(tok.op, tok.n) for tok in body_ops))
    return (_format_interval_gcigar('>', path_name, span, 0, span, body_ops), span)

def _anchor_annealed_linear_piece(template: AnnealTemplate, cigar: str) -> Optional[LinearPiece]:
    ref_span = pairwise_ref_span(cigar)
    if ref_span <= 0:
        return None
    piece = PairwisePiece(kind='encoded_qry', path=parse_graphic_segments(template.qry_gcigar, 'anneal_final_q')[0].path, ref_original_index=None, qry_original_index=None, coord_start=0, coord_end=0, ref_path_start=0, ref_path_end=0, cigar=cigar)
    if template.strand == '+':
        genome_start = template.ref_start
        genome_end = genome_start + ref_span
    else:
        genome_end = template.ref_end
        genome_start = genome_end - ref_span
    if genome_start < 0:
        return None
    return LinearPiece(piece=piece, chrom=template.chrom, strand=template.strand, genome_start=genome_start, genome_end=genome_end, cigar=cigar)

def build_stage4_linear(stage1: Stage1Result, stage3: Stage3PolishResult, graph_sequences: Dict[str, str], path_coords: Dict[str, Coord], chrom_lengths: Dict[str, int], graph_mappings: Dict[str, str], ref_reader: Optional[FastaRegionReader], ref_backbone_coord: Optional[Coord]=None, min_encode_query_span: int=50, record_sequences: Optional[Dict[str, str]]=None, query_name: Optional[str]=None) -> Stage3Result:
    query_sequence = _lookup_query_sequence(query_name or (stage1.qry_segments[0].row_name if stage1.qry_segments else None), record_sequences)
    local_duplicate_pieces = build_secondary_duplicate_lcs_pieces(stage1, graph_sequences, ref_backbone_coord)
    has_large_insertions = any(tok.op == 'I' and tok.n > min_encode_query_span for piece in stage3.pieces for tok in parse_cigar_ops(piece.cigar))
    query_coverages = build_query_coverages(stage1.qry_segments, graph_sequences, graph_mappings) if has_large_insertions else None
    if ref_backbone_coord is not None:
        main_cigar = ''.join((piece.cigar for piece in stage3.pieces))
        if not main_cigar:
            return Stage3Result([])
        ref_pieces = [piece for piece in stage3.pieces if pairwise_ref_span(piece.cigar) > 0]
        if ref_pieces:
            first_piece = ref_pieces[0]
            last_piece = ref_pieces[-1]
            main_ref_original_index = first_piece.ref_original_index
            main_qry_original_index = first_piece.qry_original_index
            main_coord_start = first_piece.coord_start
            main_coord_end = last_piece.coord_end
            main_ref_path_start = first_piece.ref_path_start
            main_ref_path_end = last_piece.ref_path_end
        else:
            main_ref_original_index = stage3.pieces[0].ref_original_index if stage3.pieces else None
            main_qry_original_index = stage3.pieces[0].qry_original_index if stage3.pieces else None
            main_coord_start = 0
            main_coord_end = 0
            main_ref_path_start = 0
            main_ref_path_end = 0
        main_piece = PairwisePiece(kind='main', path=stage1.ref_segments[0].path if stage1.ref_segments else '', ref_original_index=main_ref_original_index, qry_original_index=main_qry_original_index, coord_start=main_coord_start, coord_end=main_coord_end, ref_path_start=main_ref_path_start, ref_path_end=main_ref_path_end, cigar=main_cigar)
        main_linear = project_linear_piece_on_reference_backbone(main_piece, ref_backbone_coord, main_ref_path_start)
        if main_linear is None:
            return Stage3Result([])
        main_piece, main_linear = _rescue_internal_di_windows_multiscale(
            main_piece, main_linear, 0, query_sequence, ref_reader)
        main_piece, main_linear = _rescue_long_di_pairs_in_linear_piece(
            main_piece, main_linear, 0, query_sequence, ref_reader)
        encoded_linear = encode_large_insertions_in_piece(main_piece, 0, main_linear, query_sequence, stage1, graph_sequences, graph_mappings, ref_reader, min_encode_query_span, local_duplicate_pieces, query_coverages)
        assign_global_piece_metadata(encoded_linear)
        return Stage3Result(encoded_linear)
    linear_pieces: List[LinearPiece] = []
    qcursor = 0
    ref_segment_offsets: List[int] = []
    running_ref_offset = 0
    for seg in stage1.ref_segments:
        ref_segment_offsets.append(running_ref_offset)
        running_ref_offset += _segment_query_length(seg)
    for piece in stage3.pieces:
        piece_q0 = qcursor
        qcursor += pairwise_query_span(piece.cigar)
        piece_ref_span = pairwise_ref_span(piece.cigar)
        if piece_ref_span == 0:
            linear_pieces.extend(encode_large_insertions_in_piece(piece, piece_q0, None, query_sequence, stage1, graph_sequences, graph_mappings, ref_reader, min_encode_query_span, local_duplicate_pieces, query_coverages))
            continue
        if ref_backbone_coord is not None:
            if piece.ref_original_index is None:
                linear = None
            else:
                ref_seg = stage1.ref_segments[piece.ref_original_index]
                local_ref_offset = piece.ref_path_start - ref_seg.start
                if local_ref_offset < 0:
                    raise ValueError(f'piece local ref offset became negative: seg={piece.ref_original_index} seg_start={ref_seg.start} piece_ref_path_start={piece.ref_path_start}')
                backbone_offset = ref_segment_offsets[piece.ref_original_index] + local_ref_offset
                linear = project_linear_piece_on_reference_backbone(piece, ref_backbone_coord, backbone_offset)
        else:
            linear = project_linear_piece(piece, path_coords)
        if linear is None:
            continue
        piece, linear = _rescue_internal_di_windows_multiscale(
            piece, linear, piece_q0, query_sequence, ref_reader)
        piece, linear = _rescue_long_di_pairs_in_linear_piece(
            piece, linear, piece_q0, query_sequence, ref_reader)
        linear_pieces.extend(encode_large_insertions_in_piece(piece, piece_q0, linear, query_sequence, stage1, graph_sequences, graph_mappings, ref_reader, min_encode_query_span, local_duplicate_pieces, query_coverages))
    assign_global_piece_metadata(linear_pieces)
    return Stage3Result(linear_pieces=linear_pieces)

def pairwise_ref_span(cigar: str) -> int:
    span = 0
    for tok in parse_cigar_ops(cigar):
        if tok.op in {'=', 'X', 'D'}:
            span += tok.n
    return span

def pairwise_query_span(cigar: str) -> int:
    span = 0
    for tok in parse_cigar_ops(cigar):
        if tok.op in {'=', 'X', 'I'}:
            span += tok.n
    return span


def default_allowed_reference_chroms(pairs: Sequence[PairRow]) -> Set[str]:
    """Return the primary/reference chromosomes seen on the reference side.

    This is intentionally based on the pair/reference coordinates, not on every
    sequence present in the reference FASTA.  For example, if the pair file uses
    chr1..chr22/chrX as reference backbones, chr2 can still be emitted inside a
    chr1 row, but chr2_alt will be de-aligned unless it appears as a reference
    backbone in the pair file or --add-alt is used.
    """
    out: Set[str] = set()
    for pair in pairs:
        coord = pair.ref_coord or parse_coord(pair.ref_coord_text)
        if coord is not None and coord.chrom:
            out.add(coord.chrom)
    return out


def _query_payload_from_record_sequence(query_sequence: Optional[str], q0: int, q1: int, context: str) -> str:
    """Return exact query bases for [q0,q1), or raise a hard error.

    De-aligning an alternate-chromosome alignment to an insertion must preserve
    the inserted payload.  The query FASTA is mandatory for this pipeline, so a
    missing/short sequence is a real bug and must not silently become Ns.
    """
    q0 = int(q0)
    q1 = int(q1)
    if q1 < q0:
        raise ValueError(f'{context}: invalid query interval {q0}-{q1}')
    if query_sequence is None:
        raise ValueError(f'{context}: missing query sequence while de-aligning non-reference chromosome segment')
    if q0 < 0 or q1 > len(query_sequence):
        raise ValueError(f'{context}: query interval {q0}-{q1} out of bounds for query sequence length {len(query_sequence)}')
    payload = query_sequence[q0:q1]
    if len(payload) != q1 - q0:
        raise ValueError(f'{context}: failed to recover exact query payload for interval {q0}-{q1}')
    return payload


def _dealign_disallowed_linear_pieces(pieces: Sequence[LinearPiece], allowchroms: Optional[Set[str]], query_sequence: Optional[str], query_name: str='') -> List[LinearPiece]:
    """Turn aligned pieces on non-allowed chromosomes into plain insertions.

    The pieces are in query order.  We therefore recover the exact insertion
    payload from the row query FASTA using a running query cursor, instead of
    trying to synthesize bases from the piece CIGAR.  If that sequence is not
    available, this raises a hard error by design.
    """
    if not allowchroms:
        return list(pieces)
    out: List[LinearPiece] = []
    qpos = 0
    label = query_name or 'query'
    for piece in pieces:
        qspan = pairwise_query_span(piece.cigar)
        ref_span = pairwise_ref_span(piece.cigar)
        if ref_span > 0 and piece.chrom not in allowchroms:
            payload = _query_payload_from_record_sequence(query_sequence, qpos, qpos + qspan, f'{label} {piece.chrom}:{piece.genome_start}-{piece.genome_end}')
            if qspan > 0:
                cigar = f'{qspan}I{payload}'
                new_pairwise = PairwisePiece(kind='plain_ins', path='', ref_original_index=None, qry_original_index=piece.piece.qry_original_index, coord_start=0, coord_end=0, ref_path_start=0, ref_path_end=0, cigar=cigar)
                out.append(LinearPiece(piece=new_pairwise, chrom='', strand='+', genome_start=0, genome_end=0, cigar=cigar))
        else:
            out.append(piece)
        qpos += qspan
    return out

def _insertion_anchor_belongs_to_path_range(segment: GraphicSegment, anchor: int, path_start: int, path_end: int) -> bool:
    if path_end < path_start:
        path_start, path_end = (path_end, path_start)
    return path_start < anchor <= path_end or anchor == path_start == segment.start

def _slice_segment_ops_by_path_range(segment: GraphicSegment, path_start: int, path_end: int) -> List[CigarOp]:
    path_start = int(path_start)
    path_end = int(path_end)
    if path_end < path_start:
        path_start, path_end = (path_end, path_start)
    out: List[CigarOp] = []
    rpos = segment.start if segment.direction != '<' else segment.end

    def add(tok: CigarOp) -> None:
        if tok.n <= 0:
            return
        if out and out[-1].op == tok.op:
            old = out[-1]
            payload = _merge_payload_for_same_op(tok.op, old.n, old.payload, tok.n, tok.payload)
            out[-1] = CigarOp(old.n + tok.n, old.op, payload if tok.op in {'I', 'X', 'D'} else '')
        else:
            out.append(tok)
    for tok in _strip_terminal_h(segment.ops):
        q_consume = query_consume(tok.op, tok.n)
        r_consume = ref_consume_pair(tok.op, tok.n)
        if tok.op == 'I' and q_consume > 0 and (r_consume == 0):
            if _insertion_anchor_belongs_to_path_range(segment, rpos, path_start, path_end):
                add(CigarOp(tok.n, tok.op, tok.payload))
            continue
        if r_consume <= 0:
            continue
        if segment.direction == '<':
            op_start = rpos - r_consume
            op_end = rpos
        else:
            op_start = rpos
            op_end = rpos + r_consume
        ov0 = max(path_start, op_start)
        ov1 = min(path_end, op_end)
        if ov1 > ov0:
            take = ov1 - ov0
            if segment.direction == '<':
                offset = op_end - ov1
            else:
                offset = ov0 - op_start
            payload = ''
            if tok.op == 'X' and tok.payload:
                payload = _slice_x_payload(tok.payload, tok.n, offset, take)
            elif tok.op in {'D', 'I'} and tok.payload and (len(tok.payload) == tok.n):
                payload = tok.payload[offset:offset + take]
            add(CigarOp(take, tok.op, payload))
        if segment.direction == '<':
            rpos -= r_consume
        else:
            rpos += r_consume
    return out

def _slice_segment_body_by_path_range(segment: GraphicSegment, path_start: int, path_end: int) -> str:
    return ''.join((_op_to_cigar(tok) for tok in _slice_segment_ops_by_path_range(segment, path_start, path_end)))

def _slice_graphic_segment_by_path_range(segment: GraphicSegment, path_start: int, path_end: int) -> str:
    body = _slice_segment_ops_by_path_range(segment, path_start, path_end)
    if not body:
        return ''
    path_start = int(path_start)
    path_end = int(path_end)
    if path_end < path_start:
        path_start, path_end = (path_end, path_start)
    return _format_interval_gcigar(segment.direction, segment.path, segment.path_len, path_start, path_end, body)

def _segment_query_offsets(segments: Sequence[GraphicSegment]) -> List[int]:
    offsets: List[int] = []
    qpos = 0
    for segment in segments:
        offsets.append(qpos)
        qpos += _segment_query_length(segment)
    return offsets

def _segment_path_range_query_bounds(segment: GraphicSegment, path_start: int, path_end: int) -> Tuple[int, int]:
    path_start = int(path_start)
    path_end = int(path_end)
    if path_end < path_start:
        path_start, path_end = (path_end, path_start)
    qpos = 0
    rpos = segment.start if segment.direction != '<' else segment.end
    points: List[int] = []
    for tok in _strip_terminal_h(segment.ops):
        q_consume = query_consume(tok.op, tok.n)
        r_consume = ref_consume_pair(tok.op, tok.n)
        if tok.op == 'I' and q_consume > 0 and (r_consume == 0):
            if _insertion_anchor_belongs_to_path_range(segment, rpos, path_start, path_end):
                points.extend([qpos, qpos + tok.n])
            qpos += q_consume
            continue
        if r_consume > 0:
            if segment.direction == '<':
                op_start = rpos - r_consume
                op_end = rpos
            else:
                op_start = rpos
                op_end = rpos + r_consume
            ov0 = max(path_start, op_start)
            ov1 = min(path_end, op_end)
            if ov1 > ov0:
                take = ov1 - ov0
                if q_consume > 0:
                    if segment.direction == '<':
                        offset = op_end - ov1
                    else:
                        offset = ov0 - op_start
                    points.extend([qpos + offset, qpos + offset + take])
                else:
                    points.append(qpos)
            if segment.direction == '<':
                rpos -= r_consume
            else:
                rpos += r_consume
        qpos += q_consume
    if not points:
        return (qpos, qpos)
    return (min(points), max(points))

def _ordered_series_query_bounds(stage1: Stage1Result, side: str, ordered_indices: Sequence[int]) -> Tuple[int, int]:
    ordered = stage1.ref_ordered if side == 'ref' else stage1.qry_ordered
    segments = stage1.ref_segments if side == 'ref' else stage1.qry_segments
    offsets = _segment_query_offsets(segments)
    points: List[int] = []
    for ordered_idx in ordered_indices:
        if ordered_idx < 0 or ordered_idx >= len(ordered):
            continue
        chunk = ordered[ordered_idx]
        if chunk.original_index < 0 or chunk.original_index >= len(segments):
            continue
        local_q0, local_q1 = _segment_path_range_query_bounds(segments[chunk.original_index], _stage2_coord_start(chunk), _stage2_coord_end(chunk))
        base = offsets[chunk.original_index]
        points.extend([base + local_q0, base + local_q1])
    if not points:
        return (0, 0)
    return (min(points), max(points))


def _split_graph_cigar_at_core(
    graph_cigar: str,
    core_bounds: Tuple[int, int],
    row_name: str,
) -> str:
    """Make the two core boundaries explicit graph-segment boundaries."""
    total = graph_cigar_query_span(graph_cigar, row_name)
    core_start, core_end = map(int, core_bounds)
    if core_start < 0 or core_end <= core_start or core_end > total:
        raise ValueError(
            f"{row_name}: invalid core query interval "
            f"{core_start}-{core_end}/{total}"
        )
    pieces: List[str] = []
    for start, end in (
        (0, core_start),
        (core_start, core_end),
        (core_end, total),
    ):
        if end > start:
            pieces.append(slice_graph_cigar_by_query(
                graph_cigar, start, end, row_name,
            ))
    return "".join(pieces)


def _ordered_chunk_query_bounds(
    stage1: Stage1Result, side: str,
) -> List[Tuple[int, int]]:
    """Return every ordered LCS chunk's interval on that side's query axis."""
    if side == "ref":
        ordered = stage1.ref_ordered
        segments = stage1.ref_segments
    elif side == "qry":
        ordered = stage1.qry_ordered
        segments = stage1.qry_segments
    else:
        raise ValueError(f"unknown Stage1 side {side!r}")
    offsets = _segment_query_offsets(segments)
    output: List[Tuple[int, int]] = []
    for chunk in ordered:
        local_start, local_end = _segment_path_range_query_bounds(
            segments[chunk.original_index],
            _stage2_coord_start(chunk),
            _stage2_coord_end(chunk),
        )
        base = offsets[chunk.original_index]
        output.append((base + local_start, base + local_end))
    return output


def _core_partitioned_indices(
    bounds: Sequence[Tuple[int, int]],
    core_bounds: Tuple[int, int],
) -> Tuple[List[int], List[int], List[int]]:
    """Partition ordered chunks into upstream/core/downstream query regions."""
    core_start, core_end = core_bounds
    upstream: List[int] = []
    core_indices: List[int] = []
    downstream: List[int] = []
    for index, (start, end) in enumerate(bounds):
        if end <= core_start:
            upstream.append(index)
        elif start >= core_end:
            downstream.append(index)
        else:
            core_indices.append(index)
    return upstream, core_indices, downstream


def build_stage1_core_first(
    ref_gcigar: str,
    qry_gcigar: str,
    ref_name: str,
    qry_name: str,
    ref_core_bounds: Tuple[int, int],
    qry_core_bounds: Tuple[int, int],
    max_merge_gap: int = 100,
    reverse_ref_view: bool = False,
) -> Stage1Result:
    """Build Stage1 with the unextended core insulated from both flanks.

    The core is aligned first.  The downstream token streams are aligned in
    normal order, while the upstream streams are reversed before LCS so their
    core-facing (right) edge wins ambiguous ties.  Flank repeats therefore
    cannot replace a better core match.
    """
    split_ref = _split_graph_cigar_at_core(
        ref_gcigar, ref_core_bounds, f"{ref_name} reference core",
    )
    split_qry = _split_graph_cigar_at_core(
        qry_gcigar, qry_core_bounds, f"{qry_name} query core",
    )
    stage1 = build_stage1(
        split_ref,
        split_qry,
        ref_name,
        qry_name,
        max_merge_gap=max_merge_gap,
        reverse_ref_view=reverse_ref_view,
    )

    ref_regions = _core_partitioned_indices(
        _ordered_chunk_query_bounds(stage1, "ref"), ref_core_bounds,
    )
    qry_regions = _core_partitioned_indices(
        _ordered_chunk_query_bounds(stage1, "qry"), qry_core_bounds,
    )
    weight_by_id: Dict[int, int] = {}
    for chunk in stage1.ref_unified:
        weight_by_id[chunk.token_id] = max(
            weight_by_id.get(chunk.token_id, 0), chunk.span,
        )
    for chunk in stage1.qry_unified:
        weight_by_id[chunk.token_id] = max(
            weight_by_id.get(chunk.token_id, 0), chunk.span,
        )

    def align_region(
        ref_indices: Sequence[int],
        qry_indices: Sequence[int],
        *, reverse: bool,
    ) -> List[Tuple[int, int]]:
        ref_order = list(ref_indices)
        qry_order = list(qry_indices)
        if reverse:
            ref_order.reverse()
            qry_order.reverse()
        local_pairs = weighted_lcs(
            [stage1.ref_tokens[index] for index in ref_order],
            [stage1.qry_tokens[index] for index in qry_order],
            weight_by_id,
        )
        return [
            (ref_order[ref_index], qry_order[qry_index])
            for ref_index, qry_index in local_pairs
        ]

    upstream = align_region(ref_regions[0], qry_regions[0], reverse=True)
    core_pairs = align_region(ref_regions[1], qry_regions[1], reverse=False)
    downstream = align_region(ref_regions[2], qry_regions[2], reverse=False)
    combined = upstream + core_pairs + downstream
    combined.sort(key=lambda pair: (pair[1], pair[0]))
    stage1.lcs_pairs = combined
    return stage1


EDGE_ALIGNMENT_CONFIDENCE_SCORE = 100
EDGE_ALIGNMENT_GAP_PENALTY = 4


def alignment_quality_operation_bounds(
    ops: Sequence[CigarOp],
    minimum_score: int = EDGE_ALIGNMENT_CONFIDENCE_SCORE,
    gap_penalty: int = EDGE_ALIGNMENT_GAP_PENALTY,
) -> Optional[Tuple[int, int]]:
    """Return the operation range retained by bidirectional edge scoring.

    Scanning inward from either edge, ``=``/``M`` contributes its length and
    each mismatch, insertion, or deletion contributes ``-(length + 4)``.  A
    negative score clips through that operation and resets the score. Once the
    score is greater than 100, that edge is locked and all remaining internal
    operations are retained. If a short alignment never reaches 100, the
    positive range surviving the complete scan is retained. Operations outside
    the first/last match are never eligible.
    """
    minimum_score = int(minimum_score)
    gap_penalty = int(gap_penalty)
    if minimum_score < 0:
        raise ValueError('minimum alignment confidence score must be >= 0')
    if gap_penalty < 0:
        raise ValueError('alignment gap penalty must be >= 0')
    match_indices = [
        index for index, operation in enumerate(ops)
        if operation.op in {'=', 'M'} and operation.n > 0
    ]
    if not match_indices:
        return None
    first_match = match_indices[0]
    last_match = match_indices[-1]

    def contribution(operation: CigarOp) -> int:
        if operation.op in {'=', 'M'}:
            return operation.n
        if operation.op in {'X', 'I', 'D'}:
            return -(operation.n + gap_penalty)
        return 0

    left = first_match
    score = 0
    for index in range(first_match, last_match + 1):
        score += contribution(ops[index])
        if score < 0:
            left = index + 1
            score = 0
            continue
        if score > minimum_score:
            break
    right = last_match + 1
    score = 0
    for index in range(last_match, first_match - 1, -1):
        score += contribution(ops[index])
        if score < 0:
            right = index
            score = 0
            continue
        if score > minimum_score:
            break
    if right <= left:
        return None
    return left, right


def _core_pairwise_anchor_window(
    pairwise_cigar: str,
    minimum_score: int = EDGE_ALIGNMENT_CONFIDENCE_SCORE,
) -> Optional[Tuple[str, int, int]]:
    """Return a score-clipped core CIGAR and retained reference interval.

    Weak query sequence removed from either edge remains represented as an
    insertion, while weak reference sequence is reassigned to the corresponding
    reference extension. Coordinates are local to the unpolished core.
    """
    ops = parse_cigar_ops(pairwise_cigar)
    bounds = alignment_quality_operation_bounds(
        ops, minimum_score=minimum_score,
    )
    if bounds is None:
        return None
    first_index, right_index = bounds
    last_index = right_index - 1
    prefix = ops[:first_index]
    middle = ops[first_index:last_index + 1]
    suffix = ops[last_index + 1:]
    ref_start = sum(
        tok.n for tok in prefix if tok.op in {'=', 'X', 'D'}
    )
    ref_end = ref_start + sum(
        tok.n for tok in middle if tok.op in {'=', 'X', 'D'}
    )
    prefix_query = sum(
        tok.n for tok in prefix if tok.op in {'=', 'X', 'I'}
    )
    suffix_query = sum(
        tok.n for tok in suffix if tok.op in {'=', 'X', 'I'}
    )
    anchored: List[CigarOp] = []
    if prefix_query:
        _append_coalesced_op(anchored, CigarOp(prefix_query, 'I', ''))
    for tok in middle:
        _append_coalesced_op(anchored, tok)
    if suffix_query:
        _append_coalesced_op(anchored, CigarOp(suffix_query, 'I', ''))
    return (
        ''.join(_op_to_cigar(tok) for tok in anchored),
        ref_start,
        ref_end,
    )


def _align_graph_cigar_extension(
    ref_gcigar: str,
    qry_gcigar: str,
    ref_name: str,
    qry_name: str,
    path_sequences: Dict[str, str],
    *,
    outward_from_right_edge: bool,
    reverse_ref_view: bool = False,
) -> str:
    """Run graph-path then base alignment for one complete flank.

    Upstream flanks are reversed on both sides before alignment, making the
    core-facing right edge the first edge considered by weighted LCS and by the
    following base-alignment stage. The pairwise result is restored to the
    original left-to-right orientation before it is returned.
    """
    ref_span = graph_cigar_query_span(ref_gcigar, ref_name) if ref_gcigar else 0
    qry_span = graph_cigar_query_span(qry_gcigar, qry_name) if qry_gcigar else 0
    if ref_span == 0 and qry_span == 0:
        return ''
    if ref_span == 0:
        return f'{qry_span}I' if qry_span else ''
    if qry_span == 0:
        return f'{ref_span}D' if ref_span else ''

    aligned_ref = ref_gcigar
    aligned_qry = qry_gcigar
    if outward_from_right_edge:
        aligned_ref = reverse_mapping_gcigar_for_reverse_path(ref_gcigar)
        aligned_qry = reverse_mapping_gcigar_for_reverse_path(qry_gcigar)
    stage1 = build_stage1(
        aligned_ref,
        aligned_qry,
        ref_name,
        qry_name,
        reverse_ref_view=reverse_ref_view,
    )
    cigar = build_stage2(stage1, path_sequences).pairwise_cigar
    if outward_from_right_edge:
        cigar = _reverse_template_cigar(cigar)
    observed_ref = pairwise_ref_span(cigar)
    observed_qry = pairwise_query_span(cigar)
    if observed_ref != ref_span or observed_qry != qry_span:
        raise ValueError(
            f'{qry_name}: extension alignment span drifted '
            f'(reference {observed_ref}/{ref_span}, query {observed_qry}/{qry_span})'
        )
    return cigar


def build_stage2_core_seeded_extensions(
    ref_gcigar: str,
    qry_gcigar: str,
    ref_name: str,
    qry_name: str,
    ref_core_bounds: Tuple[int, int],
    qry_core_bounds: Tuple[int, int],
    path_sequences: Dict[str, str],
    reverse_ref_view: bool = False,
    minimum_score: int = EDGE_ALIGNMENT_CONFIDENCE_SCORE,
) -> Tuple[Stage1Result, Optional[Stage2Result]]:
    """Align the core completely, then grow both flanks from its endpoints.

    The core first undergoes graph-path matching and Stage-2 base alignment in
    isolation. Bidirectional alignment-quality clipping determines the
    reference splice coordinates. The upstream and downstream extensions then
    independently repeat graph-path matching followed by base alignment,
    starting at those resolved core-facing coordinates. No flank can change
    the core placement.
    """
    full_ref_span = graph_cigar_query_span(ref_gcigar, ref_name)
    full_qry_span = graph_cigar_query_span(qry_gcigar, qry_name)
    ref_core_start, ref_core_end = map(int, ref_core_bounds)
    qry_core_start, qry_core_end = map(int, qry_core_bounds)
    if not (0 <= ref_core_start < ref_core_end <= full_ref_span):
        raise ValueError(
            f'{ref_name}: invalid reference core interval '
            f'{ref_core_start}-{ref_core_end}/{full_ref_span}'
        )
    if not (0 <= qry_core_start < qry_core_end <= full_qry_span):
        raise ValueError(
            f'{qry_name}: invalid query core interval '
            f'{qry_core_start}-{qry_core_end}/{full_qry_span}'
        )

    # Keep a full Stage1 object for downstream duplicate/liftover annotation;
    # its Stage2 output is deliberately not used for the primary alignment.
    full_stage1 = build_stage1_core_first(
        ref_gcigar,
        qry_gcigar,
        ref_name,
        qry_name,
        ref_core_bounds,
        qry_core_bounds,
        reverse_ref_view=reverse_ref_view,
    )

    core_ref_gcigar = slice_graph_cigar_by_query(
        ref_gcigar, ref_core_start, ref_core_end,
        f'{ref_name} isolated core',
    )
    core_qry_gcigar = slice_graph_cigar_by_query(
        qry_gcigar, qry_core_start, qry_core_end,
        f'{qry_name} isolated core',
    )
    core_stage1 = build_stage1(
        core_ref_gcigar,
        core_qry_gcigar,
        f'{ref_name} core',
        f'{qry_name} core',
        reverse_ref_view=reverse_ref_view,
    )
    core_pairwise = build_stage2(core_stage1, path_sequences).pairwise_cigar
    anchored_core = _core_pairwise_anchor_window(
        core_pairwise, minimum_score=minimum_score,
    )
    if anchored_core is None:
        # Extensions are evidence only when they grow from a usable core.
        return full_stage1, None
    core_cigar, core_ref_local_start, core_ref_local_end = anchored_core

    upstream_ref_end = ref_core_start + core_ref_local_start
    downstream_ref_start = ref_core_start + core_ref_local_end
    upstream_ref_gcigar = (
        slice_graph_cigar_by_query(
            ref_gcigar, 0, upstream_ref_end,
            f'{ref_name} upstream extension',
        )
        if upstream_ref_end > 0 else ''
    )
    upstream_qry_gcigar = (
        slice_graph_cigar_by_query(
            qry_gcigar, 0, qry_core_start,
            f'{qry_name} upstream extension',
        )
        if qry_core_start > 0 else ''
    )
    downstream_ref_gcigar = (
        slice_graph_cigar_by_query(
            ref_gcigar, downstream_ref_start, full_ref_span,
            f'{ref_name} downstream extension',
        )
        if downstream_ref_start < full_ref_span else ''
    )
    downstream_qry_gcigar = (
        slice_graph_cigar_by_query(
            qry_gcigar, qry_core_end, full_qry_span,
            f'{qry_name} downstream extension',
        )
        if qry_core_end < full_qry_span else ''
    )

    upstream_cigar = _align_graph_cigar_extension(
        upstream_ref_gcigar,
        upstream_qry_gcigar,
        f'{ref_name} upstream',
        f'{qry_name} upstream',
        path_sequences,
        outward_from_right_edge=True,
        reverse_ref_view=reverse_ref_view,
    )
    downstream_cigar = _align_graph_cigar_extension(
        downstream_ref_gcigar,
        downstream_qry_gcigar,
        f'{ref_name} downstream',
        f'{qry_name} downstream',
        path_sequences,
        outward_from_right_edge=False,
        reverse_ref_view=reverse_ref_view,
    )
    pairwise_cigar = _concat_cigar_bodies(
        upstream_cigar, core_cigar, downstream_cigar,
    )
    observed_ref = pairwise_ref_span(pairwise_cigar)
    observed_qry = pairwise_query_span(pairwise_cigar)
    if observed_ref != full_ref_span or observed_qry != full_qry_span:
        raise ValueError(
            f'{qry_name}: core-seeded alignment span drifted '
            f'(reference {observed_ref}/{full_ref_span}, '
            f'query {observed_qry}/{full_qry_span})'
        )
    main_piece = PairwisePiece(
        kind='main',
        path=full_stage1.ref_segments[0].path if full_stage1.ref_segments else '',
        ref_original_index=0 if full_stage1.ref_segments else None,
        qry_original_index=0 if full_stage1.qry_segments else None,
        coord_start=0,
        coord_end=full_ref_span,
        ref_path_start=0,
        ref_path_end=full_ref_span,
        cigar=pairwise_cigar,
    )
    return full_stage1, Stage2Result(
        matched_runs=[], pieces=[main_piece], pairwise_cigar=pairwise_cigar,
    )

def _intervals_touch_or_overlap(a0: int, a1: int, b0: int, b1: int) -> bool:
    if a1 < a0:
        a0, a1 = (a1, a0)
    if b1 < b0:
        b0, b1 = (b1, b0)
    return max(a0, b0) <= min(a1, b1)

def _gcigar_from_ordered_indices(stage1: Stage1Result, side: str, ordered_indices: Sequence[int]) -> str:
    ordered = stage1.ref_ordered if side == 'ref' else stage1.qry_ordered
    segments = stage1.ref_segments if side == 'ref' else stage1.qry_segments
    if not ordered_indices:
        return ''
    chunks = [ordered[i] for i in ordered_indices if 0 <= i < len(ordered)]
    if not chunks:
        return ''
    pieces: List[str] = []
    i = 0
    while i < len(chunks):
        chunk = chunks[i]
        orig = chunk.original_index
        run_start = _stage2_coord_start(chunk)
        run_end = _stage2_coord_end(chunk)
        j = i + 1
        while j < len(chunks):
            nxt = chunks[j]
            if nxt.original_index != orig:
                break
            nxt_start = _stage2_coord_start(nxt)
            nxt_end = _stage2_coord_end(nxt)
            if not _intervals_touch_or_overlap(run_start, run_end, nxt_start, nxt_end):
                break
            run_start = min(run_start, nxt_start)
            run_end = max(run_end, nxt_end)
            j += 1
        if 0 <= orig < len(segments):
            gcigar = _slice_graphic_segment_by_path_range(segments[orig], run_start, run_end)
            if gcigar:
                pieces.append(gcigar)
        i = j
    return ''.join(pieces)

def _stage1_weight_by_token(stage1: Stage1Result) -> Dict[int, int]:
    weight_by_id: Dict[int, int] = {}
    for chunk in stage1.ref_unified:
        weight_by_id[chunk.token_id] = max(weight_by_id.get(chunk.token_id, 0), chunk.span)
    for chunk in stage1.qry_unified:
        weight_by_id[chunk.token_id] = max(weight_by_id.get(chunk.token_id, 0), chunk.span)
    return weight_by_id

def _secondary_lcs_match_bases(pairs: Sequence[Tuple[int, int]], ref_ids: Sequence[int], weight_by_id: Dict[int, int]) -> int:
    return sum((weight_by_id.get(ref_ids[ref_i], 1) for ref_i, _qry_i in pairs))

def _first_secondary_duplicate_group(stage1: Stage1Result, remaining_qry: Set[int], ref_token_ids: Set[int], min_match_bases: int) -> Optional[List[int]]:
    i = 0
    while i < len(stage1.qry_ordered):
        if i not in remaining_qry:
            i += 1
            continue
        group: List[int] = []
        shared_match_bases = 0
        while i < len(stage1.qry_ordered) and i in remaining_qry:
            group.append(i)
            chunk = stage1.qry_ordered[i]
            if chunk.token_id in ref_token_ids:
                shared_match_bases += chunk.span
            i += 1
        if group and shared_match_bases > min_match_bases:
            return group
    return None

def _build_one_secondary_duplicate_piece(stage1: Stage1Result, ref_ordered_indices: Sequence[int], qry_ordered_indices: Sequence[int], path_sequences: Dict[str, str], ref_backbone_coord: Coord, match_bases: int) -> Optional[LocalDuplicateEncodedPiece]:
    ref_gcigar = _gcigar_from_ordered_indices(stage1, 'ref', ref_ordered_indices)
    qry_gcigar = _gcigar_from_ordered_indices(stage1, 'qry', qry_ordered_indices)
    if not ref_gcigar or not qry_gcigar:
        return None
    local_stage1 = build_stage1(ref_gcigar, qry_gcigar, 'secondary_dup_ref', 'secondary_dup_qry', max_merge_gap=1)
    local_stage2 = build_stage2(local_stage1, path_sequences)
    raw_cigar = local_stage2.pairwise_cigar
    if not raw_cigar or pairwise_ref_span(raw_cigar) <= 0:
        return None
    q0, q1 = _ordered_series_query_bounds(stage1, 'qry', qry_ordered_indices)
    ref_q0, _ref_q1 = _ordered_series_query_bounds(stage1, 'ref', ref_ordered_indices)
    if q1 <= q0:
        return None
    piece = PairwisePiece(kind='encoded_qry', path=stage1.ref_segments[0].path if stage1.ref_segments else '', ref_original_index=None, qry_original_index=None, coord_start=0, coord_end=pairwise_ref_span(raw_cigar), ref_path_start=0, ref_path_end=pairwise_ref_span(raw_cigar), cigar=raw_cigar)
    piece = polish_match_piece(piece)
    ref_span = pairwise_ref_span(piece.cigar)
    if ref_span <= 0:
        return None
    if pairwise_query_span(piece.cigar) != q1 - q0:
        return None
    ref_offset = ref_q0 + piece.ref_path_start
    if ref_backbone_coord.strand == '+':
        genome_start = ref_backbone_coord.start + ref_offset
        genome_end = genome_start + ref_span
    else:
        genome_end = ref_backbone_coord.end - ref_offset
        genome_start = genome_end - ref_span
    if genome_start < 0:
        return None
    linear = LinearPiece(piece=piece, chrom=ref_backbone_coord.chrom, strand=ref_backbone_coord.strand, genome_start=genome_start, genome_end=genome_end, cigar=piece.cigar)
    return LocalDuplicateEncodedPiece(q0=q0, q1=q1, linear_piece=linear, match_bases=match_bases)

def build_secondary_duplicate_lcs_pieces(stage1: Stage1Result, path_sequences: Dict[str, str], ref_backbone_coord: Optional[Coord], min_match_bases: int=_SECONDARY_LCS_MIN_MATCH_BASES) -> List[LocalDuplicateEncodedPiece]:
    """Encode duplicate/translocation query chunks missed by the primary LCS.

        The primary Stage1 LCS is left untouched.  This pass repeatedly scans the
        remaining query chunks that failed that LCS, finds the first consecutive
        failed block with >min_match_bases of graph-template token support in the
        reference, re-aligns that block to the full reference chunk stream with
        the same weighted LCS, and turns the resulting local LCS/window into one
        encoded insertion piece anchored on the reference backbone.  Query chunks
        not consumed by these secondary pieces are left for the existing global
        liftover insertion encoder.
        """
    if ref_backbone_coord is None:
        return []
    if not stage1.ref_ordered or not stage1.qry_ordered:
        return []
    matched_qry = {qry_i for _ref_i, qry_i in stage1.lcs_pairs}
    remaining_qry: Set[int] = {i for i in range(len(stage1.qry_ordered)) if i not in matched_qry}
    ref_token_ids = set(stage1.ref_tokens)
    weight_by_id = _stage1_weight_by_token(stage1)
    out: List[LocalDuplicateEncodedPiece] = []
    blocked: Set[int] = set()
    iteration = 0
    while True:
        searchable = remaining_qry - blocked
        group = _first_secondary_duplicate_group(stage1, searchable, ref_token_ids, min_match_bases)
        if group is None:
            break
        iteration += 1
        qry_ids = [stage1.qry_ordered[i].token_id for i in group]
        lcs_pairs = weighted_lcs(stage1.ref_tokens, qry_ids, weight_by_id)
        match_bases = _secondary_lcs_match_bases(lcs_pairs, stage1.ref_tokens, weight_by_id)
        if not lcs_pairs or match_bases <= min_match_bases:
            blocked.update(group)
            continue
        local_q_indexes = [qry_i for _ref_i, qry_i in lcs_pairs]
        q_window_start = min(local_q_indexes)
        q_window_end = max(local_q_indexes) + 1
        qry_window = group[q_window_start:q_window_end]
        ref_indexes = [ref_i for ref_i, _qry_i in lcs_pairs]
        ref_window = list(range(min(ref_indexes), max(ref_indexes) + 1))
        try:
            encoded = _build_one_secondary_duplicate_piece(stage1, ref_window, qry_window, path_sequences, ref_backbone_coord, match_bases)
        except Exception as exc:
            encoded = None
        if encoded is None:
            blocked.update(group)
            continue
        out.append(encoded)
        remaining_qry.difference_update(qry_window)
        blocked.difference_update(qry_window)
    out.sort(key=lambda item: (item.q0, item.q1, -item.match_bases))
    non_overlapping: List[LocalDuplicateEncodedPiece] = []
    cursor = -1
    for item in out:
        if item.q0 < cursor:
            continue
        non_overlapping.append(item)
        cursor = item.q1
    return non_overlapping

def encode_query_interval_with_local_duplicates(insert_q0: int, insert_q1: int, insert_payload: Optional[str], query_sequence: Optional[str], stage1: Stage1Result, graph_sequences: Dict[str, str], graph_mappings: Dict[str, str], ref_reader: Optional[FastaRegionReader], local_duplicate_pieces: Optional[Sequence[LocalDuplicateEncodedPiece]], query_coverages: Optional[Sequence[QueryCoverage]]=None) -> List[LinearPiece]:
    if not local_duplicate_pieces:
        return encode_query_interval_by_reference(insert_q0, insert_q1, insert_payload, query_sequence, stage1, graph_sequences, graph_mappings, ref_reader, query_coverages)
    usable = [item for item in local_duplicate_pieces if insert_q0 <= item.q0 and item.q1 <= insert_q1 and (item.q1 > item.q0)]
    usable.sort(key=lambda item: (item.q0, item.q1))
    if not usable:
        return encode_query_interval_by_reference(insert_q0, insert_q1, insert_payload, query_sequence, stage1, graph_sequences, graph_mappings, ref_reader, query_coverages)
    out: List[LinearPiece] = []
    cursor = insert_q0
    for item in usable:
        if item.q0 < cursor:
            continue
        if item.q0 > cursor:
            frag_payload = _payload_fragment(insert_payload, insert_q0, cursor, item.q0)
            out.extend(encode_query_interval_by_reference(cursor, item.q0, frag_payload, query_sequence, stage1, graph_sequences, graph_mappings, ref_reader, query_coverages))
        out.append(item.linear_piece)
        cursor = item.q1
    if cursor < insert_q1:
        frag_payload = _payload_fragment(insert_payload, insert_q0, cursor, insert_q1)
        out.extend(encode_query_interval_by_reference(cursor, insert_q1, frag_payload, query_sequence, stage1, graph_sequences, graph_mappings, ref_reader, query_coverages))
    return out


_LONG_DI_REALIGN_THRESHOLD = 1000
_LONG_DI_MIN_MATCHES = 1000
_SHORT_DI_REALIGN_MIN_SPAN = 31
_SMALL_DI_REALIGN_MAX_SPAN = 30
_SMALL_DI_MATCH_SCORE = 1
_SMALL_DI_MISMATCH_SCORE = -1
_SMALL_DI_GAP_OPEN_SCORE = -4
_SMALL_DI_GAP_EXTEND_SCORE = -1
_COMPOSITE_DI_BROAD_EXACT_ANCHOR = 1000
_COMPOSITE_DI_LOCAL_EXACT_ANCHOR = 31


def _realign_small_di_sequences(
    ref_seq: str,
    qry_seq: str,
) -> Optional[List[Tuple[int, str, str]]]:
    """Globally align one small adjacent D/I pair with affine gap costs.

    Only pairs whose reference and query sides are both at most 30 bp are
    eligible.  A one-base gap costs -4 and each additional base in that same gap
    costs -1.  Mismatch and match scores are -1 and +1, respectively.  X and I
    payloads contain query bases, matching the final graph-CIGAR convention.
    """
    if (
        not ref_seq
        or not qry_seq
        or max(len(ref_seq), len(qry_seq)) > _SMALL_DI_REALIGN_MAX_SPAN
    ):
        return None

    match_state = 0
    delete_state = 1
    insert_state = 2
    unreachable = -10**9
    n = len(ref_seq)
    m = len(qry_seq)

    match = [[unreachable] * (m + 1) for _ in range(n + 1)]
    delete = [[unreachable] * (m + 1) for _ in range(n + 1)]
    insert = [[unreachable] * (m + 1) for _ in range(n + 1)]
    trace_match = [[-1] * (m + 1) for _ in range(n + 1)]
    trace_delete = [[-1] * (m + 1) for _ in range(n + 1)]
    trace_insert = [[-1] * (m + 1) for _ in range(n + 1)]
    match[0][0] = 0

    for i in range(1, n + 1):
        if i == 1:
            delete[i][0] = _SMALL_DI_GAP_OPEN_SCORE
            trace_delete[i][0] = match_state
        else:
            delete[i][0] = delete[i - 1][0] + _SMALL_DI_GAP_EXTEND_SCORE
            trace_delete[i][0] = delete_state
    for j in range(1, m + 1):
        if j == 1:
            insert[0][j] = _SMALL_DI_GAP_OPEN_SCORE
            trace_insert[0][j] = match_state
        else:
            insert[0][j] = insert[0][j - 1] + _SMALL_DI_GAP_EXTEND_SCORE
            trace_insert[0][j] = insert_state

    def choose(candidates: Sequence[Tuple[int, int]]) -> Tuple[int, int]:
        # Candidate order is the deterministic tie breaker.
        best_score, best_state = candidates[0]
        for score, state in candidates[1:]:
            if score > best_score:
                best_score, best_state = score, state
        return best_score, best_state

    for i in range(1, n + 1):
        rb = ref_seq[i - 1]
        for j in range(1, m + 1):
            qb = qry_seq[j - 1]
            substitution = (
                _SMALL_DI_MATCH_SCORE
                if rb.upper() == qb.upper()
                else _SMALL_DI_MISMATCH_SCORE
            )
            previous, trace_match[i][j] = choose((
                (match[i - 1][j - 1], match_state),
                (delete[i - 1][j - 1], delete_state),
                (insert[i - 1][j - 1], insert_state),
            ))
            match[i][j] = previous + substitution

            delete[i][j], trace_delete[i][j] = choose((
                (delete[i - 1][j] + _SMALL_DI_GAP_EXTEND_SCORE, delete_state),
                (match[i - 1][j] + _SMALL_DI_GAP_OPEN_SCORE, match_state),
                (insert[i - 1][j] + _SMALL_DI_GAP_OPEN_SCORE, insert_state),
            ))
            insert[i][j], trace_insert[i][j] = choose((
                (insert[i][j - 1] + _SMALL_DI_GAP_EXTEND_SCORE, insert_state),
                (match[i][j - 1] + _SMALL_DI_GAP_OPEN_SCORE, match_state),
                (delete[i][j - 1] + _SMALL_DI_GAP_OPEN_SCORE, delete_state),
            ))

    _score, state = choose((
        (match[n][m], match_state),
        (delete[n][m], delete_state),
        (insert[n][m], insert_state),
    ))
    reverse_ops: List[Tuple[int, str, str]] = []
    i = n
    j = m
    while i > 0 or j > 0:
        if state == match_state:
            if i <= 0 or j <= 0:
                raise ValueError('invalid small D/I alignment traceback through match')
            rb = ref_seq[i - 1]
            qb = qry_seq[j - 1]
            previous_state = trace_match[i][j]
            reverse_ops.append(
                (1, '=', '') if rb.upper() == qb.upper() else (1, 'X', qb)
            )
            i -= 1
            j -= 1
        elif state == delete_state:
            if i <= 0:
                raise ValueError('invalid small D/I alignment traceback through deletion')
            previous_state = trace_delete[i][j]
            reverse_ops.append((1, 'D', ''))
            i -= 1
        elif state == insert_state:
            if j <= 0:
                raise ValueError('invalid small D/I alignment traceback through insertion')
            previous_state = trace_insert[i][j]
            reverse_ops.append((1, 'I', qry_seq[j - 1]))
            j -= 1
        else:
            raise ValueError('invalid small D/I alignment traceback state')
        state = previous_state

    out: List[Tuple[int, str, str]] = []
    for size, op, payload in reversed(reverse_ops):
        _tm_add_op(out, size, op, payload)
    return out


def _realign_long_di_sequences(ref_seq: str, qry_seq: str) -> Optional[List[Tuple[int, str, str]]]:
    """Return a supported global alignment for one adjacent D/I sequence pair.

    Graph-CIGAR input often omits the payload after an ``I`` operation.  Stage 2
    must therefore represent a homologous reference/query pair as ``D`` + ``I``
    temporarily.  Once the piece has been projected to real reference
    coordinates, the two sequences are available and can be compared directly.
    Short pairs are aligned with the in-process global DP implementation; long
    pairs continue to use minimap2.  A rescue is accepted only when there is
    substantial sequence support, so a genuine unrelated deletion/insertion
    remains D/I.
    """
    if not ref_seq or not qry_seq:
        return None
    min_span = min(len(ref_seq), len(qry_seq))
    if min_span < _SHORT_DI_REALIGN_MIN_SPAN:
        return None

    # Different graph paths frequently carry an identical replacement.  Avoid
    # invoking either aligner for this common and unambiguous case.
    if len(ref_seq) == len(qry_seq) and ref_seq.upper() == qry_seq.upper():
        return [(len(ref_seq), '=', '')]

    max_span = max(len(ref_seq), len(qry_seq))
    if max_span < _LONG_DI_REALIGN_THRESHOLD:
        ops, _r0, _q0, _r1, _q1 = _tm_global_insert_align_dp(
            ref_seq, qry_seq,
        )
    else:
        try:
            ops = minimap2_payload_ops(ref_seq, qry_seq)
        except (FileNotFoundError, RuntimeError, ValueError):
            # Never send a >=1 kb rescue to quadratic DP/SSW.  The original
            # D/I representation is conservative when minimap2 is unavailable.
            return None
    matches = sum(n for n, op, _payload in ops if op == '=')
    aligned = sum(n for n, op, _payload in ops if op in {'=', 'X'})
    required_matches = min(
        _LONG_DI_MIN_MATCHES,
        max(_SHORT_DI_REALIGN_MIN_SPAN, (min_span + 1) // 2),
    )
    if matches < required_matches:
        return None
    if aligned < required_matches:
        return None
    if matches < 0.5 * min_span:
        return None
    return ops


def _rescue_internal_di_windows_in_linear_piece(
    piece: PairwisePiece,
    linear: LinearPiece,
    piece_q0: int,
    query_sequence: Optional[str],
    ref_reader: Optional[FastaRegionReader],
    strong_exact_min: int = 31,
) -> Tuple[PairwisePiece, LinearPiece]:
    """Realign a composite D/I window bounded by strong exact anchors.

    A graph traversal can describe one homologous replacement as many deletion
    operations separated by tiny exact matches, followed later by an insertion.
    The adjacent-D/I rescue cannot see that these operations belong to one
    replacement.  This pass globally realigns the complete internal window
    against the authoritative query and reference FASTAs.  Unanchored leading
    and trailing windows remain untouched.
    """
    if ref_reader is None or query_sequence is None:
        return (piece, linear)
    ops = parse_cigar_ops(piece.cigar)
    strong_indices = [
        i for i, tok in enumerate(ops)
        if tok.op == '=' and tok.n >= strong_exact_min
    ]
    if len(strong_indices) < 2:
        return (piece, linear)

    next_strong = {
        left: right for left, right in zip(strong_indices, strong_indices[1:])
    }
    out: List[CigarOp] = []
    rpos = linear.genome_start if linear.strand == '+' else linear.genome_end
    qpos = int(piece_q0)
    changed = False

    def append_and_advance(tok: CigarOp) -> None:
        nonlocal rpos, qpos
        _append_coalesced_op(out, tok)
        if tok.op in {'=', 'X', 'D'}:
            rpos = rpos + tok.n if linear.strand == '+' else rpos - tok.n
        if tok.op in {'=', 'X', 'I'}:
            qpos += tok.n

    i = 0
    while i < len(ops):
        tok = ops[i]
        right = next_strong.get(i)
        if right is None:
            append_and_advance(tok)
            i += 1
            continue

        # Keep the left anchor. rpos/qpos then identify the uncertain window.
        append_and_advance(tok)
        window = ops[i + 1:right]
        ref_len = sum(
            item.n for item in window if item.op in {'=', 'X', 'D'}
        )
        qry_len = sum(
            item.n for item in window if item.op in {'=', 'X', 'I'}
        )
        has_deletion = any(item.op == 'D' for item in window)
        has_insertion = any(item.op == 'I' for item in window)
        rescued_ops: Optional[List[Tuple[int, str, str]]] = None
        if (
            has_deletion
            and has_insertion
            and min(ref_len, qry_len) >= _SHORT_DI_REALIGN_MIN_SPAN
            and rpos >= 0
            and qpos >= 0
            and qpos + qry_len <= len(query_sequence)
        ):
            try:
                ref_seq = _fetch_linear_piece_ref_chunk(
                    ref_reader, linear, rpos, ref_len,
                )[0]
            except Exception:
                ref_seq = ''
            qry_seq = query_sequence[qpos:qpos + qry_len]
            if len(ref_seq) == ref_len and len(qry_seq) == qry_len:
                candidate = _realign_long_di_sequences(ref_seq, qry_seq)
                if candidate is not None:
                    candidate_ref_span = sum(
                        n for n, op, _payload in candidate
                        if op in {'=', 'X', 'D'}
                    )
                    candidate_qry_span = sum(
                        n for n, op, _payload in candidate
                        if op in {'=', 'X', 'I'}
                    )
                    if (
                        candidate_ref_span == ref_len
                        and candidate_qry_span == qry_len
                    ):
                        rescued_ops = candidate

        if rescued_ops is None:
            for item in window:
                append_and_advance(item)
        else:
            for n, op, payload in rescued_ops:
                append_and_advance(CigarOp(
                    int(n), op, payload if op in {'I', 'X'} else '',
                ))
            changed = True
        i = right

    if not changed:
        return (piece, linear)
    new_cigar = ''.join(_op_to_cigar(tok) for tok in out)
    src = piece
    new_piece = PairwisePiece(
        kind=src.kind,
        path=src.path,
        ref_original_index=src.ref_original_index,
        qry_original_index=src.qry_original_index,
        coord_start=src.coord_start,
        coord_end=src.coord_end,
        ref_path_start=src.ref_path_start,
        ref_path_end=src.ref_path_end,
        cigar=new_cigar,
    )
    new_linear = LinearPiece(
        piece=new_piece,
        chrom=linear.chrom,
        strand=linear.strand,
        genome_start=linear.genome_start,
        genome_end=linear.genome_end,
        cigar=new_cigar,
        piece_index=linear.piece_index,
        piece_size=linear.piece_size,
        main_split_parent=linear.main_split_parent,
    )
    return (new_piece, new_linear)


def _rescue_internal_di_windows_multiscale(
    piece: PairwisePiece,
    linear: LinearPiece,
    piece_q0: int,
    query_sequence: Optional[str],
    ref_reader: Optional[FastaRegionReader],
) -> Tuple[PairwisePiece, LinearPiece]:
    """Rescue broad repeat replacements before considering local windows.

    Exact runs of only a few dozen or hundred bases are common inside tandem
    repeats and are not reliable boundaries for a large replacement.  A first
    pass bounded by >=1 kb exact anchors lets the authoritative FASTAs resolve
    the entire repeat block.  The existing 31-bp pass then handles genuinely
    local D/I mixtures that the broad pass cannot bracket.
    """
    for minimum_exact in (
        _COMPOSITE_DI_BROAD_EXACT_ANCHOR,
        _COMPOSITE_DI_LOCAL_EXACT_ANCHOR,
    ):
        piece, linear = _rescue_internal_di_windows_in_linear_piece(
            piece,
            linear,
            piece_q0,
            query_sequence,
            ref_reader,
            strong_exact_min=minimum_exact,
        )
    return (piece, linear)


def _rescue_long_di_pairs_in_linear_piece(
    piece: PairwisePiece,
    linear: LinearPiece,
    piece_q0: int,
    query_sequence: Optional[str],
    ref_reader: Optional[FastaRegionReader],
) -> Tuple[PairwisePiece, LinearPiece]:
    """Replace an eligible adjacent D/I pair with a sequence alignment.

    This is deliberately placed after reference projection and before large-I
    encoding.  At that point ``linear`` supplies the exact reference interval,
    while ``query_sequence`` supplies the omitted graph-CIGAR I payload.  The
    original stage1/stage2 matching is unchanged.  Pairs with both sides at most
    30 bp use affine global alignment; pairs with both sides at least 31 bp keep
    the existing sequence-supported long-pair rescue.
    """
    if ref_reader is None or query_sequence is None:
        return (piece, linear)
    ops = parse_cigar_ops(piece.cigar)
    if len(ops) < 2:
        return (piece, linear)

    out: List[CigarOp] = []
    rpos = linear.genome_start if linear.strand == '+' else linear.genome_end
    qpos = int(piece_q0)
    changed = False
    i = 0

    def append_and_advance(tok: CigarOp) -> None:
        nonlocal rpos, qpos
        _append_coalesced_op(out, tok)
        if tok.op in {'=', 'X', 'D'}:
            rpos = rpos + tok.n if linear.strand == '+' else rpos - tok.n
        if tok.op in {'=', 'X', 'I'}:
            qpos += tok.n

    while i < len(ops):
        first = ops[i]
        if first.op not in {'D', 'I'}:
            append_and_advance(first)
            i += 1
            continue

        first_op = first.op
        second_op = 'I' if first_op == 'D' else 'D'
        j = i
        first_len = 0
        while j < len(ops) and ops[j].op == first_op:
            first_len += ops[j].n
            j += 1
        k = j
        second_len = 0
        while k < len(ops) and ops[k].op == second_op:
            second_len += ops[k].n
            k += 1
        if second_len <= 0:
            append_and_advance(first)
            i += 1
            continue

        ref_len = first_len if first_op == 'D' else second_len
        qry_len = second_len if first_op == 'D' else first_len
        is_small_pair = max(ref_len, qry_len) <= _SMALL_DI_REALIGN_MAX_SPAN
        is_supported_long_pair = (
            min(ref_len, qry_len) >= _SHORT_DI_REALIGN_MIN_SPAN
        )
        if not is_small_pair and not is_supported_long_pair:
            append_and_advance(first)
            i += 1
            continue
        if rpos < 0 or qpos < 0 or qpos + qry_len > len(query_sequence):
            append_and_advance(first)
            i += 1
            continue
        try:
            ref_seq = _fetch_linear_piece_ref_chunk(ref_reader, linear, rpos, ref_len)[0]
        except Exception:
            append_and_advance(first)
            i += 1
            continue
        qry_seq = query_sequence[qpos:qpos + qry_len]
        if len(ref_seq) != ref_len or len(qry_seq) != qry_len:
            append_and_advance(first)
            i += 1
            continue
        if is_small_pair:
            rescued = _realign_small_di_sequences(ref_seq, qry_seq)
        else:
            rescued = _realign_long_di_sequences(ref_seq, qry_seq)
        if rescued is None:
            append_and_advance(first)
            i += 1
            continue
        for n, op, payload in rescued:
            append_and_advance(CigarOp(int(n), op, payload if op in {'I', 'X'} else ''))
        changed = True
        i = k

    if not changed:
        return (piece, linear)
    new_cigar = ''.join((_op_to_cigar(tok) for tok in out))
    src = piece
    new_piece = PairwisePiece(
        kind=src.kind,
        path=src.path,
        ref_original_index=src.ref_original_index,
        qry_original_index=src.qry_original_index,
        coord_start=src.coord_start,
        coord_end=src.coord_end,
        ref_path_start=src.ref_path_start,
        ref_path_end=src.ref_path_end,
        cigar=new_cigar,
    )
    new_linear = LinearPiece(
        piece=new_piece,
        chrom=linear.chrom,
        strand=linear.strand,
        genome_start=linear.genome_start,
        genome_end=linear.genome_end,
        cigar=new_cigar,
        piece_index=linear.piece_index,
        piece_size=linear.piece_size,
        main_split_parent=linear.main_split_parent,
    )
    return (new_piece, new_linear)

def _piece_match_size_from_cigar(cigar: str) -> int:
    """Return the raw piece size used by main annealing.

                This intentionally counts only strong match bases (=/M). It is computed
                before emission-time gap filling or overlap trimming, so jump decisions do
                not change when serialization edits a CIGAR boundary.
                """
    return sum((tok.n for tok in parse_cigar_ops(cigar) if tok.op == '='))

def assign_global_piece_metadata(linear_pieces: Sequence[LinearPiece]) -> None:
    """Assign stable post-encoding indexes and raw match sizes to all pieces."""
    for i, piece in enumerate(linear_pieces):
        piece.piece_index = i
        piece.piece_size = _piece_match_size_from_cigar(piece.cigar)

@dataclass
class PieceEmitPlan:
    linked_left: bool = False
    linked_right: bool = False
    prefix: str = ''
    suffix: str = ''
    bridge_suffix: str = ''
    trim_prefix_ref: int = 0
    trim_suffix_ref: int = 0
    claimed_by_main: bool = False
    claimed_by_insertion: bool = False
    linked_left_from_index: Optional[int] = None
    linked_right_to_index: Optional[int] = None

def _linear_piece_interval(piece: LinearPiece) -> Tuple[int, int]:
    a = min(piece.genome_start, piece.genome_end)
    b = max(piece.genome_start, piece.genome_end)
    return (a, b)

def _interval_contains_or_same(left: LinearPiece, right: LinearPiece) -> bool:
    a0, a1 = _linear_piece_interval(left)
    b0, b1 = _linear_piece_interval(right)
    len_a = a1 - a0
    len_b = b1 - b0
    if len_a <= 0 or len_b <= 0:
        return False
    outer_span = max(a1, b1) - min(a0, b0)
    return outer_span == max(len_a, len_b)

def _linear_piece_gap_or_overlap(left: LinearPiece, right: LinearPiece) -> Tuple[str, int, int, int]:
    """Return relation, size, gap_start, gap_end for two genomic intervals.

                The returned relation is "gap", "overlap", or "touch".  gap_start/end are
                forward genomic coordinates for the bridge interval when relation == "gap".
                """
    a0, a1 = _linear_piece_interval(left)
    b0, b1 = _linear_piece_interval(right)
    if a1 < b0:
        return ('gap', b0 - a1, a1, b0)
    if b1 < a0:
        return ('gap', a0 - b1, b1, a0)
    overlap = min(a1, b1) - max(a0, b0)
    if overlap > 0:
        return ('overlap', overlap, 0, 0)
    return ('touch', 0, 0, 0)

def _linear_pieces_can_anneal(left: LinearPiece, right: LinearPiece, max_abs_gap: int=200) -> bool:
    if left.chrom != right.chrom or left.strand != right.strand:
        return False
    if pairwise_ref_span(left.cigar) <= 0 or pairwise_ref_span(right.cigar) <= 0:
        return False
    if left.piece_size <= 0 or right.piece_size <= 0:
        return False
    if _interval_contains_or_same(left, right):
        return False
    relation, size, _gap0, _gap1 = _linear_piece_gap_or_overlap(left, right)
    if relation == 'touch':
        return True
    return size < max_abs_gap

def _fetch_linear_gap_by_coords(chrom: str, strand: str, gap_start: int, gap_end: int, ref_reader: Optional[FastaRegionReader]) -> str:
    if gap_end <= gap_start:
        return ''
    if ref_reader is None or chrom not in ref_reader:
        return 'N' * (gap_end - gap_start)
    return ref_reader.fetch(chrom, gap_start, gap_end, strand)

def _append_coalesced_op(out: List[CigarOp], tok: CigarOp) -> None:
    if tok.n <= 0:
        return
    if out and out[-1].op == tok.op:
        old = out[-1]
        payload = _merge_payload_for_same_op(tok.op, old.n, old.payload, tok.n, tok.payload)
        out[-1] = CigarOp(old.n + tok.n, tok.op, payload if tok.op in {'I', 'X', 'D'} else '')
    else:
        out.append(tok)

def _fetch_piece_ref_bases_by_offset(piece: LinearPiece, ref_offset: int, n: int, ref_reader: Optional[FastaRegionReader]) -> str:
    if n <= 0:
        return ''
    if ref_reader is None or not piece.chrom or piece.chrom not in ref_reader:
        return 'N' * n
    if piece.strand == '+':
        a = piece.genome_start + ref_offset
        b = a + n
        return ref_reader.fetch(piece.chrom, a, b, '+')
    b = piece.genome_end - ref_offset
    a = b - n
    return ref_reader.fetch(piece.chrom, a, b, '-')

def _removed_query_bases(piece: LinearPiece, tok: CigarOp, op_offset: int, take: int, ref_offset: int, ref_reader: Optional[FastaRegionReader]) -> str:
    if take <= 0:
        return ''
    if tok.op == '=':
        return _fetch_piece_ref_bases_by_offset(piece, ref_offset, take, ref_reader)
    if tok.op == 'X':
        payload = _slice_x_payload(tok.payload, tok.n, op_offset, take) if tok.payload else ''
        return _x_payload_query(payload, take)
    if tok.op == 'I':
        if tok.payload and tok.payload != 'n':
            return tok.payload[op_offset:op_offset + take]
        return 'N' * take
    return ''

def _trim_cigar_prefix_by_ref(piece: LinearPiece, cigar: str, trim_ref: int, ref_reader: Optional[FastaRegionReader]) -> Tuple[str, str]:
    if trim_ref <= 0:
        return (cigar, '')
    ops = parse_cigar_ops(cigar)
    remaining = int(trim_ref)
    ref_offset = 0
    out: List[CigarOp] = []
    removed: List[str] = []
    for tok in ops:
        ref_n = _pairwise_ref_consume_op(tok)
        if remaining > 0:
            if tok.op == 'I':
                removed.append(_removed_query_bases(piece, tok, 0, tok.n, ref_offset, ref_reader))
                continue
            if ref_n > 0:
                take = min(remaining, ref_n)
                removed.append(_removed_query_bases(piece, tok, 0, take, ref_offset, ref_reader))
                keep = tok.n - take
                ref_offset += take
                remaining -= take
                if keep > 0:
                    payload = ''
                    if tok.op == 'X' and tok.payload:
                        payload = _slice_x_payload(tok.payload, tok.n, take, keep)
                    elif tok.op == 'D' and tok.payload:
                        payload = tok.payload[take:take + keep]
                    _append_coalesced_op(out, CigarOp(keep, tok.op, payload))
                    ref_offset += keep
                continue
        _append_coalesced_op(out, tok)
        ref_offset += ref_n
    return (''.join((_op_to_cigar(tok) for tok in out)), ''.join(removed))

def _trim_cigar_suffix_by_ref(piece: LinearPiece, cigar: str, trim_ref: int, ref_reader: Optional[FastaRegionReader]) -> Tuple[str, str]:
    if trim_ref <= 0:
        return (cigar, '')
    ops = parse_cigar_ops(cigar)
    total_ref = pairwise_ref_span(cigar)
    keep_ref_end = max(0, total_ref - int(trim_ref))
    ref_offset = 0
    out: List[CigarOp] = []
    removed: List[str] = []
    for tok in ops:
        ref_n = _pairwise_ref_consume_op(tok)
        ref0 = ref_offset
        ref1 = ref_offset + ref_n
        if tok.op == 'I':
            if ref0 >= keep_ref_end:
                removed.append(_removed_query_bases(piece, tok, 0, tok.n, ref0, ref_reader))
            else:
                _append_coalesced_op(out, tok)
            continue
        if ref_n == 0:
            _append_coalesced_op(out, tok)
            continue
        if ref1 <= keep_ref_end:
            _append_coalesced_op(out, tok)
        elif ref0 >= keep_ref_end:
            removed.append(_removed_query_bases(piece, tok, 0, tok.n, ref0, ref_reader))
        else:
            keep = keep_ref_end - ref0
            remove = ref1 - keep_ref_end
            if keep > 0:
                payload = ''
                if tok.op == 'X' and tok.payload:
                    payload = _slice_x_payload(tok.payload, tok.n, 0, keep)
                elif tok.op == 'D' and tok.payload:
                    payload = tok.payload[:keep]
                _append_coalesced_op(out, CigarOp(keep, tok.op, payload))
            removed.append(_removed_query_bases(piece, tok, keep, remove, keep_ref_end, ref_reader))
        ref_offset = ref1
    return (''.join((_op_to_cigar(tok) for tok in out)), ''.join(removed))

def _main_side_chain(pieces: Sequence[LinearPiece], anchor_index: int, forward: bool, max_abs_gap: int=200) -> List[int]:
    """Find main-annealed pieces on one side of the anchor.

                The endpoint stack implements the requested jump rule: a later candidate may
                link to an earlier endpoint only when its raw match size is larger than the
                sum of the skipped endpoint sizes.  The returned indexes exclude anchor_index.
                """
    anchor = pieces[anchor_index]
    endpoints: List[int] = [anchor_index]
    valid_size = 1
    iterator = range(anchor_index + 1, len(pieces)) if forward else range(anchor_index - 1, -1, -1)
    for idx in iterator:
        piece = pieces[idx]
        if piece.piece.kind != 'main':
            continue
        if pairwise_ref_span(piece.cigar) <= 0 or piece.piece_size <= 0:
            continue
        piece_size_sum = 0
        matched = False
        for pos in range(valid_size - 1, -1, -1):
            ep_idx = endpoints[pos]
            left_idx, right_idx = (ep_idx, idx) if forward else (idx, ep_idx)
            if piece.piece_size > piece_size_sum and _linear_pieces_can_anneal(pieces[left_idx], pieces[right_idx], max_abs_gap):
                new_valid_size = pos + 2
                if len(endpoints) < new_valid_size:
                    endpoints.append(idx)
                else:
                    endpoints[pos + 1] = idx
                valid_size = new_valid_size
                del endpoints[valid_size:]
                matched = True
                break
            if piece_size_sum > piece.piece_size:
                break
            piece_size_sum += pieces[ep_idx].piece_size
        if matched:
            continue
    return endpoints[1:valid_size]

def _main_anneal_group(pieces: Sequence[LinearPiece], max_abs_gap: int=200) -> List[int]:
    main_indexes = [i for i, piece in enumerate(pieces) if piece.piece.kind == 'main' and pairwise_ref_span(piece.cigar) > 0]
    if not main_indexes:
        return []
    anchor_index = main_indexes[0]
    rev_side = _main_side_chain(pieces, anchor_index, forward=False, max_abs_gap=max_abs_gap)
    fwd_side = _main_side_chain(pieces, anchor_index, forward=True, max_abs_gap=max_abs_gap)
    return list(reversed(rev_side)) + [anchor_index] + fwd_side

def _add_link_adjustment(left: LinearPiece, right: LinearPiece, plans: Dict[int, PieceEmitPlan], ref_reader: Optional[FastaRegionReader]) -> None:
    left_plan = plans[left.piece_index]
    right_plan = plans[right.piece_index]
    left_plan.linked_right = True
    right_plan.linked_left = True
    left_plan.linked_right_to_index = right.piece_index
    right_plan.linked_left_from_index = left.piece_index
    relation, size, gap0, gap1 = _linear_piece_gap_or_overlap(left, right)
    if relation == 'gap' and size > 0:
        gap_seq = _fetch_linear_gap_by_coords(left.chrom, left.strand, gap0, gap1, ref_reader)
        left_plan.bridge_suffix += f'{size}D{gap_seq}'
    elif relation == 'overlap' and size > 0:
        if left.strand == '+':
            right_plan.trim_prefix_ref += size
        else:
            left_plan.trim_suffix_ref += size

def _build_emit_plans(pieces: Sequence[LinearPiece], ref_reader: Optional[FastaRegionReader], max_abs_gap: int=200) -> Dict[int, PieceEmitPlan]:
    plans: Dict[int, PieceEmitPlan] = {piece.piece_index: PieceEmitPlan() for piece in pieces}
    main_group = _main_anneal_group(pieces, max_abs_gap=max_abs_gap)
    claimed_by_main = set(main_group)
    for idx in claimed_by_main:
        plans[idx].claimed_by_main = True
    for left_idx, right_idx in zip(main_group, main_group[1:]):
        _add_link_adjustment(pieces[left_idx], pieces[right_idx], plans, ref_reader)

    # Only main-chain pieces are annealed/merged across segment boundaries.
    # Encoded insertion alignments are deliberately left as independent anchored
    # segments, even when two encoded_qry pieces are adjacent or separated only
    # by zero-reference insertion payloads.  This keeps encoded insertion calls
    # stable and avoids silently changing their segment structure.
    return plans

def _adjust_body_for_emit(piece: LinearPiece, plan: PieceEmitPlan, ref_reader: Optional[FastaRegionReader]) -> str:
    body = piece.cigar
    prefix_insert = ''
    suffix_insert = ''
    if plan.trim_prefix_ref > 0:
        body, removed = _trim_cigar_prefix_by_ref(piece, body, plan.trim_prefix_ref, ref_reader)
        if removed:
            prefix_insert = _concat_cigar_bodies(prefix_insert, f'{len(removed)}I{removed}')
    if plan.trim_suffix_ref > 0:
        body, removed = _trim_cigar_suffix_by_ref(piece, body, plan.trim_suffix_ref, ref_reader)
        if removed:
            suffix_insert = _concat_cigar_bodies(suffix_insert, f'{len(removed)}I{removed}')
    return _concat_cigar_bodies(plan.prefix, prefix_insert, body, suffix_insert, plan.suffix, plan.bridge_suffix)

def _query_payload_from_x(tok: CigarOp) -> str:
    if not tok.payload:
        return 'N' * tok.n
    return _x_payload_query(tok.payload, tok.n)

def _trimmed_query_payload(ops: Sequence[CigarOp]) -> str:
    out: List[str] = []
    for tok in ops:
        if tok.op == 'I':
            out.append(tok.payload)
        elif tok.op == 'X':
            out.append(_query_payload_from_x(tok))
        elif tok.op == '=':
            raise ValueError('cannot reconstruct exact payload from = during polish')
        elif tok.op in {'D', 'H'}:
            continue
        else:
            raise ValueError(f'cannot collapse weak flank containing op {tok.op}')
    return ''.join(out)

def _trimmed_ref_span(ops: Sequence[CigarOp]) -> int:
    span = 0
    for tok in ops:
        if tok.op in {'D', 'X', '='}:
            span += tok.n
        elif tok.op in {'I', 'H'}:
            continue
        else:
            raise ValueError(f'cannot trim weak flank containing op {tok.op}')
    return span

def _trimmed_query_size(ops: Sequence[CigarOp]) -> int:
    size = 0
    for tok in ops:
        if tok.op in {'X', '=', 'I'}:
            size += tok.n
        elif tok.op in {'D', 'H'}:
            continue
        else:
            raise ValueError(f'cannot collapse weak flank containing op {tok.op}')
    return size

def _op_to_cigar(tok: CigarOp) -> str:
    return f'{tok.n}{tok.op}{tok.payload}'

def _concat_cigar_bodies(*parts: str) -> str:
    """Concatenate CIGAR bodies while coalescing adjacent same ops.

    This is used when zero-reference insertion pieces are attached to the
    previous or next real segment.  In particular, adjacent insertions must
    become one payload-carrying I op, e.g. 3Iabc + 2Izz -> 5Iabczz.
    """
    out: List[CigarOp] = []
    for part in parts:
        if not part:
            continue
        for tok in parse_cigar_ops(part):
            _append_coalesced_op(out, tok)
    return ''.join(_op_to_cigar(tok) for tok in out)

def _count_only_segments(cigar: str) -> List[str]:
    return re.findall('\\d+[=XDIH]', cigar)

def _segment_size(seg: str) -> int:
    return int(seg[:-1])

def _segment_op(seg: str) -> str:
    return seg[-1]

def _find_trim_side(segments: Sequence[str], strong_match_min: int=20) -> Tuple[int, int, int]:
    length = len(segments)
    truncate_qsize = 0
    truncate_rsize = 0
    score = 0
    scorestart = 0
    truncate_qsize_start = 0
    truncate_rsize_start = 0
    i = length
    for i, seg in enumerate(segments):
        size = _segment_size(seg)
        op = _segment_op(seg)
        if op == '=' and size > strong_match_min:
            break
        if op == '=':
            score += size
        else:
            score -= 4 * size
        if score > 20:
            i = scorestart
            truncate_qsize = truncate_qsize_start
            truncate_rsize = truncate_rsize_start
            break
        if score >= 0:
            scorestart = i
            truncate_qsize_start = truncate_qsize
            truncate_rsize_start = truncate_rsize
        else:
            scorestart = length
        if op in {'X', '=', 'I'}:
            truncate_qsize += size
        if op in {'D', 'X', '='}:
            truncate_rsize += size
    return (i, truncate_qsize, truncate_rsize)

def polish_match_piece(piece: PairwisePiece, strong_match_min: int=20) -> PairwisePiece:
    ops = parse_cigar_ops(piece.cigar)
    segments = _count_only_segments(piece.cigar)
    if len(segments) != len(ops):
        raise ValueError('count-only cigar segmentation drifted from parsed ops')
    if not ops:
        return piece
    if not any((tok.op == '=' and tok.n > strong_match_min for tok in ops)):
        return piece
    left_boundary, left_qry_size, left_ref_trim = _find_trim_side(segments, strong_match_min)
    right_boundary, right_qry_size, right_ref_trim = _find_trim_side(list(reversed(segments)), strong_match_min)
    if left_boundary + right_boundary > len(ops):
        return piece
    prefix = ops[:left_boundary]
    middle = ops[left_boundary:len(ops) - right_boundary]
    suffix = ops[len(ops) - right_boundary:]
    if left_ref_trim != _trimmed_ref_span(prefix) or left_qry_size != _trimmed_query_size(prefix):
        raise ValueError('left polish trim accounting drifted from parsed ops')
    if right_ref_trim != _trimmed_ref_span(suffix) or right_qry_size != _trimmed_query_size(suffix):
        raise ValueError('right polish trim accounting drifted from parsed ops')
    rebuilt: List[str] = []
    if left_qry_size:
        try:
            left_payload = _trimmed_query_payload(prefix)
        except ValueError:
            left_payload = ''
        rebuilt.append(f'{left_qry_size}I{left_payload}' if left_payload else f'{left_qry_size}I')
    for tok in middle:
        rebuilt.append(_op_to_cigar(tok))
    if right_qry_size:
        try:
            right_payload = _trimmed_query_payload(suffix)
        except ValueError:
            right_payload = ''
        rebuilt.append(f'{right_qry_size}I{right_payload}' if right_payload else f'{right_qry_size}I')
    return PairwisePiece(kind=piece.kind, path=piece.path, ref_original_index=piece.ref_original_index, qry_original_index=piece.qry_original_index, coord_start=piece.coord_start + left_ref_trim, coord_end=piece.coord_end - right_ref_trim, ref_path_start=piece.ref_path_start + left_ref_trim, ref_path_end=piece.ref_path_end - right_ref_trim, cigar=''.join(rebuilt))

def anchor_linear_piece(piece: LinearPiece, chrom_lengths: Dict[str, int]) -> Tuple[str, str, int, str, int]:
    chrom_len = chrom_lengths.get(piece.chrom, 2000000000)
    if chrom_len is None:
        raise KeyError(f'missing chromosome length for {piece.chrom}')
    ref_span = pairwise_ref_span(piece.cigar)
    if ref_span == 0:
        raise ValueError('cannot anchor a zero-reference-span piece directly')
    if ref_span != piece.genome_end - piece.genome_start:
        raise ValueError(f'linear piece span mismatch for {piece.chrom}:{piece.genome_start}-{piece.genome_end}: ref_span={ref_span}')
    if piece.strand == '+':
        direction = '>'
        left_h = piece.genome_start
        right_h = chrom_len - piece.genome_end
    else:
        direction = '<'
        left_h = chrom_len - piece.genome_end
        right_h = piece.genome_start
    return (direction, piece.chrom, left_h, piece.cigar, right_h)

def main_anchor_text_from_coord(coord: Coord, chrom_lengths: Dict[str, int]) -> str:
    """Return the dummy first segment that declares the main/backbone alignment.

                The first graphic-CIGAR segment is the main alignment indicator.  If a
                leading encoded insertion/reference-consuming piece appears before the first
                real main body, we still need to emit this zero-body anchor first.  It must
                come from the full reference backbone coordinate, not from an internal
                LinearPiece selected by LCS/liftover.
                """
    chrom_len = chrom_lengths.get(coord.chrom)
    if chrom_len is None:
        raise KeyError(f'missing chromosome length for {coord.chrom}')
    if coord.strand == '-':
        return f'<{coord.chrom}:{max(0, chrom_len - coord.end)}H'
    return f'>{coord.chrom}:{max(0, coord.start)}H'

def serialize_stage3_reflalign(stage3: Stage3Result, chrom_lengths: Dict[str, int], ref_reader: Optional[FastaRegionReader]=None, main_anchor_coord: Optional[Coord]=None, allowchroms: Optional[Set[str]]=None, query_sequence: Optional[str]=None, query_name: str='') -> str:
    pieces = stage3.linear_pieces
    if not pieces:
        return ''
    pieces = _dealign_disallowed_linear_pieces(pieces, allowchroms, query_sequence, query_name=query_name)
    if not pieces:
        return ''
    assign_global_piece_metadata(pieces)
    plans = _build_emit_plans(pieces, ref_reader)
    piece_by_index = {piece.piece_index: piece for piece in pieces}
    first_ref_piece_index: Optional[int] = None
    first_main_piece_index: Optional[int] = None
    first_main_anchor_text = ''
    for anchor_piece in pieces:
        if pairwise_ref_span(anchor_piece.cigar) <= 0:
            continue
        if first_ref_piece_index is None:
            first_ref_piece_index = anchor_piece.piece_index
        if anchor_piece.piece.kind == 'main':
            first_main_piece_index = anchor_piece.piece_index
            if main_anchor_coord is not None:
                first_main_anchor_text = main_anchor_text_from_coord(main_anchor_coord, chrom_lengths)
            else:
                direction0, chrom0, left_h0, _body0, _right_h0 = anchor_linear_piece(anchor_piece, chrom_lengths)
                first_main_anchor_text = f'{direction0}{chrom0}:{left_h0}H'
            break
    emit_front_main_anchor = first_ref_piece_index is not None and first_main_piece_index is not None and (first_ref_piece_index != first_main_piece_index)
    pending_prefix = ''
    last_ref_piece_index: Optional[int] = None
    for piece in pieces:
        if pairwise_ref_span(piece.cigar) == 0:
            if last_ref_piece_index is None:
                pending_prefix = _concat_cigar_bodies(pending_prefix, piece.cigar)
            else:
                plans[last_ref_piece_index].suffix = _concat_cigar_bodies(plans[last_ref_piece_index].suffix, piece.cigar)
            continue
        if pending_prefix:
            plans[piece.piece_index].prefix = _concat_cigar_bodies(pending_prefix, plans[piece.piece_index].prefix)
            pending_prefix = ''
        last_ref_piece_index = piece.piece_index
    if pending_prefix:
        if main_anchor_coord is not None:
            anchor = main_anchor_text_from_coord(main_anchor_coord, chrom_lengths)
            if anchor.endswith(':0H'):
                anchor = anchor[:-2]
            return anchor + pending_prefix
        return pending_prefix
    rendered: List[str] = []
    if emit_front_main_anchor:
        rendered.append(first_main_anchor_text)
    for piece in pieces:
        if pairwise_ref_span(piece.cigar) == 0:
            continue
        plan = plans[piece.piece_index]
        direction, chrom, left_h, _body, right_h = anchor_linear_piece(piece, chrom_lengths)
        body = _adjust_body_for_emit(piece, plan, ref_reader)
        left_h_text = '' if plan.linked_left or left_h == 0 else f'{left_h}H'
        right_h_text = '' if plan.linked_right or right_h == 0 else f'{right_h}H'
        main_anchor_already_emitted = emit_front_main_anchor and piece.piece_index == first_main_piece_index
        if main_anchor_already_emitted:
            rendered.append(f'{direction}{body}{right_h_text}')
        elif plan.linked_left:
            adjacent_main_main = plan.claimed_by_main and piece.piece.kind == 'main' and (plan.linked_left_from_index == piece.piece_index - 1)
            left_link_piece = piece_by_index.get(plan.linked_left_from_index) if plan.linked_left_from_index is not None else None
            intervening_ref_piece = False
            if plan.linked_left_from_index is not None:
                lo = min(plan.linked_left_from_index, piece.piece_index) + 1
                hi = max(plan.linked_left_from_index, piece.piece_index)
                for mid_i in range(lo, hi):
                    mid_piece = piece_by_index.get(mid_i)
                    if mid_piece is not None and pairwise_ref_span(mid_piece.cigar) > 0:
                        intervening_ref_piece = True
                        break
            same_main_split = plan.claimed_by_main and piece.piece.kind == 'main' and (piece.main_split_parent >= 0) and (left_link_piece is not None) and (left_link_piece.piece.kind == 'main') and (left_link_piece.main_split_parent == piece.main_split_parent)
            if adjacent_main_main:
                rendered.append(f'{body}{right_h_text}')
            elif same_main_split:
                if intervening_ref_piece:
                    rendered.append(f'{direction}{body}{right_h_text}')
                else:
                    rendered.append(f'{body}{right_h_text}')
            elif plan.claimed_by_main and piece.piece.kind == 'main' and not intervening_ref_piece:
                # A zero-reference-span insertion can sit between two main pieces
                # after non-reference chromosomes are de-aligned.  That insertion
                # is attached to the left piece's suffix above, so it must not
                # force a fresh anchor on the right main piece.
                rendered.append(f'{body}{right_h_text}')
            elif plan.claimed_by_main and piece.piece.kind == 'main':
                clip = f'{left_h}H' if left_h else ''
                rendered.append(f'{direction}{chrom}:{clip}{body}{right_h_text}')
            else:
                rendered.append(f'{body}{right_h_text}')
        else:
            rendered.append(f'{direction}{chrom}:{left_h_text}{body}{right_h_text}')
    return ''.join(rendered)

def _fetch_output_ref_chunk(ref_reader: FastaRegionReader, segment: GraphicSegment, rpos: int, n: int) -> Tuple[str, int]:
    if n <= 0:
        return ('', rpos)
    if segment.direction == '>':
        return (ref_reader.fetch(segment.path, rpos, rpos + n, '+'), rpos + n)
    return (ref_reader.fetch(segment.path, rpos - n, rpos, '-'), rpos - n)

def _append_recanonicalized_op(out: List[CigarOp], op: str, n: int, payload: str='') -> None:
    if n <= 0:
        return
    if op in {'=', 'D', 'H'}:
        payload = ''
    elif op == 'X' and payload:
        payload = _x_payload_query(payload, n)
    if out and out[-1].op == op:
        last = out[-1]
        if op in {'=', 'D', 'H'}:
            merged_payload = ''
        else:
            merged_payload = _merge_payload_for_same_op(op, last.n, last.payload, n, payload)
        out[-1] = CigarOp(last.n + n, op, merged_payload)
    else:
        out.append(CigarOp(n, op, payload))

def _append_ref_query_comparison(out: List[CigarOp], ref_seq: str, query_seq: str) -> None:
    m = min(len(ref_seq), len(query_seq))
    ref_upper = ref_seq[:m].upper()
    qry_upper = query_seq[:m].upper()
    if ref_upper == qry_upper:
        _append_recanonicalized_op(out, '=', m, '')
        return
    for i, j, is_match in _match_mismatch_runs(ref_upper, qry_upper):
        if is_match:
            _append_recanonicalized_op(out, '=', j - i, '')
        else:
            _append_recanonicalized_op(out, 'X', j - i, query_seq[i:j])

def _fetch_linear_piece_ref_chunk(ref_reader: FastaRegionReader, piece: LinearPiece, rpos: int, n: int) -> Tuple[str, int]:
    if n <= 0:
        return ('', rpos)
    if not piece.chrom:
        return ('N' * n, rpos)
    if piece.strand == '+':
        return (ref_reader.fetch(piece.chrom, rpos, rpos + n, '+'), rpos + n)
    return (ref_reader.fetch(piece.chrom, rpos - n, rpos, '-'), rpos - n)

def _canonicalize_stage3_against_sequences(pair: PairRow, stage3: Stage3Result, record_sequences: Dict[str, str], ref_reader: Optional[FastaRegionReader]) -> Stage3Result:
    """Canonicalize LinearPiece CIGAR bodies before serialization.

                This is intentionally done on Stage3Result/LinearPiece objects instead of on
                the serialized graphic CIGAR string. The serialized output may contain
                pathless main-continuation chunks (for example, ``<1376=...H``), which the
                normal named graphic-CIGAR parser cannot parse. Rebuilding here preserves the
                output grammar while still forcing every =/X/I operation to agree with the
                final -r and -q sequences. Final payload convention is query-only
                payloads for X/I and no D payload.
                """
    if ref_reader is None:
        return stage3
    query_name = pair.query_name or pair.label
    actual_query = record_sequences.get(query_name) or record_sequences.get(strip_match_name(query_name)) or ''
    if not actual_query:
        return stage3
    qpos = 0
    new_pieces: List[LinearPiece] = []
    for linear in stage3.linear_pieces:
        rpos = linear.genome_start if linear.strand == '+' else linear.genome_end
        body: List[CigarOp] = []
        for tok in parse_cigar_ops(linear.cigar):
            if tok.op in {'=', 'X'}:
                try:
                    ref_chunk, rpos = _fetch_linear_piece_ref_chunk(ref_reader, linear, rpos, tok.n)
                except Exception:
                    ref_chunk = 'N' * tok.n
                query_chunk = actual_query[qpos:qpos + tok.n]
                if len(query_chunk) != tok.n:
                    query_chunk = _x_payload_query(tok.payload, tok.n) if tok.op == 'X' else ref_chunk
                qpos += tok.n
                _append_ref_query_comparison(body, ref_chunk, query_chunk)
            elif tok.op == 'D':
                try:
                    ref_chunk, rpos = _fetch_linear_piece_ref_chunk(ref_reader, linear, rpos, tok.n)
                except Exception:
                    ref_chunk = tok.payload if len(tok.payload) == tok.n else 'N' * tok.n
                _append_recanonicalized_op(body, 'D', tok.n, '')
            elif tok.op == 'I':
                query_chunk = actual_query[qpos:qpos + tok.n]
                if len(query_chunk) != tok.n:
                    query_chunk = tok.payload if len(tok.payload) == tok.n else 'N' * tok.n
                qpos += tok.n
                _append_recanonicalized_op(body, 'I', tok.n, query_chunk)
            elif tok.op == 'H':
                _append_recanonicalized_op(body, 'H', tok.n, tok.payload)
            else:
                _append_recanonicalized_op(body, tok.op, tok.n, tok.payload)
        new_cigar = ''.join((_op_to_cigar(tok) for tok in body))
        if new_cigar == linear.cigar:
            new_pieces.append(linear)
            continue
        src = linear.piece
        new_pairwise = PairwisePiece(kind=src.kind, path=src.path, ref_original_index=src.ref_original_index, qry_original_index=src.qry_original_index, coord_start=src.coord_start, coord_end=src.coord_end, ref_path_start=src.ref_path_start, ref_path_end=src.ref_path_end, cigar=new_cigar)
        new_pieces.append(LinearPiece(piece=new_pairwise, chrom=linear.chrom, strand=linear.strand, genome_start=linear.genome_start, genome_end=linear.genome_end, cigar=new_cigar))
    return Stage3Result(new_pieces)


def _rescue_serialized_di_pairs(
    graph_cigar: str,
    query_sequence: Optional[str],
    ref_reader: Optional[FastaRegionReader],
) -> str:
    """Rescue D/I pairs that become adjacent only during serialization.

    Large insertions are temporarily split into zero-reference LinearPieces.
    Consequently a reference ``D`` and its homologous query ``I`` can live in
    separate Stage-4 pieces and evade the earlier within-piece rescue.  The
    serializer attaches such a zero-reference piece back to its neighboring
    reference segment.  At that point the pair is adjacent again and can be
    checked against the authoritative query/reference FASTAs without changing
    either sequence span or the selected reference placement.
    """
    if not graph_cigar or query_sequence is None or ref_reader is None:
        return graph_cigar
    segments = parse_graphic_segments(graph_cigar, 'serialized D/I rescue')
    if not segments:
        return graph_cigar
    rebuilt: List[str] = []
    query_offset = 0
    changed = False
    for segment in segments:
        body_ops = _strip_terminal_h(segment.ops)
        body_cigar = ''.join(_op_to_cigar(tok) for tok in body_ops)
        body_ref_span = pairwise_ref_span(body_cigar)
        body_query_span = pairwise_query_span(body_cigar)
        if not body_cigar or not segment.path:
            rebuilt.append(graphic_segment_to_gcigar(segment))
            query_offset += body_query_span
            continue
        pairwise = PairwisePiece(
            kind='main',
            path=segment.path,
            ref_original_index=None,
            qry_original_index=None,
            coord_start=segment.start,
            coord_end=segment.end,
            ref_path_start=segment.start,
            ref_path_end=segment.end,
            cigar=body_cigar,
        )
        linear = LinearPiece(
            piece=pairwise,
            chrom=segment.path,
            strand='-' if segment.direction == '<' else '+',
            genome_start=segment.start,
            genome_end=segment.end,
            cigar=body_cigar,
        )
        rescued, rescued_linear = _rescue_internal_di_windows_multiscale(
            pairwise,
            linear,
            query_offset,
            query_sequence,
            ref_reader,
        )
        rescued, _rescued_linear = _rescue_long_di_pairs_in_linear_piece(
            rescued,
            rescued_linear,
            query_offset,
            query_sequence,
            ref_reader,
        )
        if _SV_REALIGNMENT_ENABLED:
            try:
                segment_ref_sequence = ref_reader.fetch(
                    segment.path,
                    segment.start,
                    segment.end,
                    '-' if segment.direction == '<' else '+',
                )
            except Exception:
                segment_ref_sequence = ''
            segment_query_sequence = query_sequence[
                query_offset:query_offset + body_query_span
            ]
            serialized_ops = [
                tok for tok in parse_cigar_ops(rescued.cigar)
                if tok.op != 'H'
            ]
            # This late pass is only for insertions that were scattered across
            # Stage-4 piece boundaries and first meet in the serialized row.
            if (
                sum(tok.op == 'I' for tok in serialized_ops) >= 2
                and len(segment_ref_sequence) == body_ref_span
                and len(segment_query_sequence) == body_query_span
            ):
                consolidated_ops = _tm_consolidate_sv_ops(
                    [
                        (
                            tok.n,
                            'M' if tok.op in {'=', 'M'} else tok.op,
                            tok.payload,
                        )
                        for tok in serialized_ops
                    ],
                    segment_ref_sequence,
                    segment_query_sequence,
                )
                final_ops: List[CigarOp] = []
                for n, op, payload in consolidated_ops:
                    emit_op = '=' if op in {'M', '='} else op
                    if emit_op == 'X' and payload:
                        payload = _x_payload_query(payload, n)
                    _append_recanonicalized_op(
                        final_ops, emit_op, n, payload,
                    )
                consolidated_cigar = ''.join(
                    _op_to_cigar(tok) for tok in final_ops
                )
                if (
                    pairwise_ref_span(consolidated_cigar) == body_ref_span
                    and pairwise_query_span(consolidated_cigar)
                    == body_query_span
                ):
                    rescued.cigar = consolidated_cigar
        if rescued.cigar != body_cigar:
            changed = True
        rescued_ops = parse_cigar_ops(rescued.cigar)
        rebuilt.append(_format_interval_gcigar(
            segment.direction,
            segment.path,
            segment.path_len,
            segment.start,
            segment.end,
            rescued_ops,
        ))
        query_offset += body_query_span
    if not changed:
        return graph_cigar
    output = ''.join(rebuilt)
    if graph_cigar_query_span(output, 'serialized D/I rescue output') != len(query_sequence):
        raise ValueError(
            'serialized D/I rescue changed the complete query span'
        )
    return output
_WORKER_GCIGARS = None
_WORKER_GRAPH_SEQUENCES = None
_WORKER_RECORD_SEQUENCES = None
_WORKER_GRAPH_PATH_COORDS = None
_WORKER_CHROM_LENGTHS = None
_WORKER_GRAPH_MAPPINGS = None
_WORKER_COORD_BY_NAME = None
_WORKER_ALLOWCHROMS = None
_WORKER_REF_READER = None

def _init_worker_state(gcigars, graph_sequences, record_sequences, graph_path_coords, chrom_lengths, graph_mappings, coord_by_name, allowchroms, reference_path: str, realignment_enabled: bool=True) -> None:
    global _WORKER_GCIGARS
    global _WORKER_GRAPH_SEQUENCES
    global _WORKER_RECORD_SEQUENCES
    global _WORKER_GRAPH_PATH_COORDS
    global _WORKER_CHROM_LENGTHS
    global _WORKER_GRAPH_MAPPINGS
    global _WORKER_COORD_BY_NAME
    global _WORKER_ALLOWCHROMS
    global _WORKER_REF_READER
    set_sv_realignment_enabled(realignment_enabled)
    _WORKER_GCIGARS = gcigars
    _WORKER_GRAPH_SEQUENCES = graph_sequences
    _WORKER_RECORD_SEQUENCES = record_sequences
    _WORKER_GRAPH_PATH_COORDS = graph_path_coords
    _WORKER_CHROM_LENGTHS = chrom_lengths
    _WORKER_GRAPH_MAPPINGS = graph_mappings
    _WORKER_COORD_BY_NAME = coord_by_name
    _WORKER_ALLOWCHROMS = set(allowchroms or [])
    _WORKER_REF_READER = FastaRegionReader(reference_path) if reference_path else None

def _worker_build_provisional_row(task: Tuple[int, PairRow]) -> Tuple[int, Optional[str], Optional[str]]:
    idx, pair = task
    try:
        row = build_provisional_row(pair, _WORKER_GCIGARS, _WORKER_GRAPH_SEQUENCES, _WORKER_RECORD_SEQUENCES, _WORKER_GRAPH_PATH_COORDS, _WORKER_CHROM_LENGTHS, _WORKER_GRAPH_MAPPINGS, _WORKER_REF_READER, _WORKER_COORD_BY_NAME, allowchroms=_WORKER_ALLOWCHROMS)
        return (idx, row, None)
    except Exception:
        label = pair.query_name or pair.label or f'line {pair.line_no}'
        msg = f'failed on pair index={idx} line={pair.line_no} label={label}:\n{traceback.format_exc()}'
        return (idx, None, msg)


def finished_pair_key_from_values(query_name: Optional[str], ref_name: Optional[str], label: str) -> Tuple[str, str, str]:
    return (query_name or '', ref_name or '', label or '')

def finished_pair_key(pair: PairRow) -> Tuple[str, str, str]:
    return finished_pair_key_from_values(
        pair.query_name,
        pair.ref_name,
        getattr(pair, 'query_alignment_regions_text', pair.label),
    )

def finished_output_row_key(cols: Sequence[str]) -> Tuple[str, str, str]:
    """Return the resume key for legacy or per-sample-v2 seven-column rows."""
    if (
        len(cols) >= 7
        and parse_coord(cols[1]) is None
        and parse_coord(cols[2].split(';', 1)[0]) is not None
    ):
        return finished_pair_key_from_values(cols[0], cols[1], cols[4])
    return finished_pair_key_from_values(cols[0], cols[2], cols[4])

def read_finished_pair_keys(save_path: str) -> Set[Tuple[str, str, str]]:
    """Read already completed legacy or per-sample-v2 output rows.

    Legacy rows key on query/reference/label (columns 1/3/5); v2 rows key on
    query/reference/query-alignment-regions (columns 1/2/5).
    """
    done: Set[Tuple[str, str, str]] = set()
    if not save_path or not os.path.exists(save_path):
        return done
    with open(save_path) as fh:
        for line_no, raw in enumerate(fh, 1):
            line = raw.rstrip('\n')
            if not line:
                continue
            cols = line.split('\t')
            if len(cols) < 5:
                sys.stderr.write(f'[graphcigartoref] warning: ignoring malformed existing output line {line_no} in {save_path}: expected >=5 columns\n')
                continue
            done.add(finished_output_row_key(cols))
    return done

def prune_existing_output_for_pairs(save_path: str, wanted_keys: Set[Tuple[str, str, str]]) -> Tuple[Set[Tuple[str, str, str]], int, int, int, int]:
    """Keep only existing output rows still requested by the current match-pair file.

    --unreplace uses the schema-specific completed-pair key described by
    ``read_finished_pair_keys``.
    This cleanup removes stale rows whose key is no longer in the current pair list,
    malformed/blank rows, and duplicate completed rows for the same key.
    It returns the keys that remain in the output after cleanup.
    """
    finished_keys: Set[Tuple[str, str, str]] = set()
    if not save_path or not os.path.exists(save_path):
        return (finished_keys, 0, 0, 0, 0)

    tmp_path = f'{save_path}.tmp.{os.getpid()}'
    total_rows = 0
    kept_rows = 0
    removed_rows = 0
    malformed_rows = 0
    try:
        with open(save_path) as in_fh, open(tmp_path, 'w') as out_fh:
            for line_no, raw in enumerate(in_fh, 1):
                total_rows += 1
                line = raw.rstrip('\n')
                if not line:
                    removed_rows += 1
                    malformed_rows += 1
                    continue
                cols = line.split('\t')
                if any(value.startswith('ALTERNATIVE_INTERVALS:Z:') for value in cols[7:]):
                    raise ValueError('cannot --unreplace interval-restricted output; rerun into a fresh output so every row uses the requested BED policy')
                if len(cols) < 5:
                    sys.stderr.write(f'[graphcigartoref] warning: removing malformed existing output line {line_no} in {save_path}: expected >=5 columns\n')
                    removed_rows += 1
                    malformed_rows += 1
                    continue
                key = finished_output_row_key(cols)
                if key not in wanted_keys:
                    removed_rows += 1
                    continue
                if key in finished_keys:
                    removed_rows += 1
                    continue
                finished_keys.add(key)
                kept_rows += 1
                out_fh.write(raw if raw.endswith('\n') else raw + '\n')
        os.replace(tmp_path, save_path)
    except Exception:
        try:
            if os.path.exists(tmp_path):
                os.unlink(tmp_path)
        finally:
            raise
    return (finished_keys, total_rows, kept_rows, removed_rows, malformed_rows)

def exclude_finished_pairs(pairs: Sequence[PairRow], finished_keys: Set[Tuple[str, str, str]]) -> List[PairRow]:
    if not finished_keys:
        return list(pairs)
    return [pair for pair in pairs if finished_pair_key(pair) not in finished_keys]


def _exact_same_coordinate_gcigar(
    pair: PairRow,
    query_sequence: Optional[str],
    ref_reader: Optional[FastaRegionReader],
) -> Optional[str]:
    """Return an exact whole-record mapping when coordinates and bases prove it.

    A graph traversal through a tandem repeat can choose different copies on
    its two sides and manufacture an internal deletion even when a query and
    reference record describe the same genomic interval, or the query is an
    exact coordinate-defined subinterval of the reference record.  This
    shortcut is deliberately strict: containment, orientation, and every
    authoritative FASTA base must agree.  It therefore cannot erase a true
    variant merely because two records have similar graph paths.
    """
    if ref_reader is None or not query_sequence:
        return None
    query_coord = pair.query_coord or parse_coord(pair.query_coord_text)
    ref_coord = pair.ref_coord or parse_coord(pair.ref_coord_text)
    if query_coord is None or ref_coord is None:
        return None
    if (
        query_coord.chrom != ref_coord.chrom
        or query_coord.strand != ref_coord.strand
        or query_coord.start < ref_coord.start
        or query_coord.end > ref_coord.end
    ):
        return None
    span = query_coord.end - query_coord.start
    if span <= 0 or len(query_sequence) != span:
        return None
    if query_coord.chrom not in ref_reader.index:
        return None
    try:
        reference_sequence = ref_reader.fetch(
            query_coord.chrom,
            query_coord.start,
            query_coord.end,
            query_coord.strand,
        )
    except Exception:
        return None
    if query_sequence.upper() != reference_sequence.upper():
        return None
    path_len = ref_reader.index[query_coord.chrom][0]
    direction = '<' if query_coord.strand == '-' else '>'
    return _format_interval_gcigar(
        direction,
        query_coord.chrom,
        path_len,
        query_coord.start,
        query_coord.end,
        [CigarOp(span, '=', '')],
    )


def restore_record_payloads(
    graph_cigar: str,
    name: str,
    record_sequences: Dict[str, str],
    coord: Optional[Coord] = None,
    reader: Optional[FastaRegionReader] = None,
) -> str:
    """Hydrate a compact row only for the duration of its comparison.

    Existing complete payloads remain usable without a sequence lookup. Never
    mutate the shared CIGAR dictionary: retaining hydrated rows there would
    duplicate large FASTA slices in every worker after fork.
    """
    missing = any(
        token.n > 0 and token.op in {'I', 'X'} and not token.payload
        for segment in parse_graphic_segments(graph_cigar, name)
        for token in segment.ops
    )
    if not missing:
        return graph_cigar
    sequence = _lookup_query_sequence(name, record_sequences)
    if sequence is None and reader is not None and coord is not None:
        if coord.chrom not in reader.index:
            raise KeyError(f'{name}: reference contig {coord.chrom!r} is missing')
        if not 0 <= coord.start <= coord.end <= reader.index[coord.chrom][0]:
            raise ValueError(f'{name}: reference interval is outside its contig: {coord}')
        sequence = reader.fetch(coord.chrom, coord.start, coord.end, coord.strand)
    if sequence is None:
        raise ValueError(
            f'{name}: compact graph CIGAR requires its oriented local FASTA '
            'sequence to restore I/X payloads'
        )
    return add_graph_cigar_payloads(graph_cigar, sequence, name)


def build_provisional_row(pair: PairRow, gcigars: Dict[str, str], graph_sequences: Dict[str, str], record_sequences: Dict[str, str], graph_path_coords: Dict[str, Coord], chrom_lengths: Dict[str, int], graph_mappings: Dict[str, str], ref_reader: Optional[FastaRegionReader], coord_by_name: Dict[str, Coord], canonicalize_output: bool=False, allowchroms: Optional[Set[str]]=None) -> str:
    if pair.query_name is None or pair.ref_name is None:
        raise ValueError(f'unresolved pair: {pair.label}')
    ref_backbone = coord_by_name.get(pair.ref_name) or pair.ref_coord
    query_sequence = _lookup_query_sequence(pair.query_name, record_sequences)
    cigar = _exact_same_coordinate_gcigar(
        pair, query_sequence, ref_reader,
    )
    reverse_ref_view = False
    if pair.query_coord is not None and pair.ref_coord is not None:
        reverse_ref_view = pair.query_coord.strand != pair.ref_coord.strand
    if cigar is None:
        # The per-sample caller has fetched the oriented query interval; the
        # standalone caller has -q records. The reference side uses its own
        # -q record when available, otherwise its oriented reference interval.
        query_gcigar = restore_record_payloads(
            gcigars[pair.query_name], pair.query_name, record_sequences,
        )
        reference_gcigar = restore_record_payloads(
            gcigars[pair.ref_name], pair.ref_name, record_sequences,
            ref_backbone, ref_reader,
        )
        query_core_bounds = getattr(pair, 'query_core_slice', None)
        reference_core_bounds = getattr(pair, 'ref_core_slice', None)
        if query_core_bounds is not None and reference_core_bounds is not None:
            result, stage2 = build_stage2_core_seeded_extensions(
                reference_gcigar,
                query_gcigar,
                pair.ref_name,
                pair.query_name,
                reference_core_bounds,
                query_core_bounds,
                graph_sequences,
                reverse_ref_view=reverse_ref_view,
            )
        else:
            result = build_stage1(reference_gcigar, query_gcigar, pair.ref_name, pair.query_name, reverse_ref_view=reverse_ref_view)
            stage2 = build_stage2(result, graph_sequences)
        if stage2 is None:
            stage4 = Stage3Result([])
        else:
            stage3 = build_stage3_polished(stage2)
            stage4 = build_stage4_linear(result, stage3, graph_sequences, graph_path_coords, chrom_lengths, graph_mappings, ref_reader, ref_backbone_coord=ref_backbone, record_sequences=record_sequences, query_name=pair.query_name)
        stage4 = _canonicalize_stage3_against_sequences(pair, stage4, record_sequences, ref_reader)
        cigar = serialize_stage3_reflalign(stage4, chrom_lengths, ref_reader=ref_reader, main_anchor_coord=ref_backbone, allowchroms=allowchroms, query_sequence=query_sequence, query_name=pair.query_name or pair.label)
        cigar = _rescue_serialized_di_pairs(cigar, query_sequence, ref_reader)
    qcoord_rec = coord_by_name.get(pair.query_name)
    rcoord_rec = coord_by_name.get(pair.ref_name)
    qcoord = f'{qcoord_rec.chrom}:{qcoord_rec.start}-{qcoord_rec.end}{qcoord_rec.strand}' if qcoord_rec is not None else pair.query_coord_text
    rcoord = f'{rcoord_rec.chrom}:{rcoord_rec.start}-{rcoord_rec.end}{rcoord_rec.strand}' if rcoord_rec is not None else pair.ref_coord_text
    cigar = query_only_graph_cigar(cigar)
    row = '\t'.join([pair.query_name, qcoord, pair.ref_name, rcoord, pair.label, pair.label_coord, cigar])
    return append_tag(row, getattr(pair, 'alternative_intervals', None))

def _default_chunksize(n_items: int, n_processes: int) -> int:
    if n_processes <= 1:
        return 1
    return max(1, min(32, (n_items + n_processes * 8 - 1) // (n_processes * 8)))

def _write_rows_serial(pairs: Sequence[PairRow], out_fh, gcigars: Dict[str, str], graph_sequences: Dict[str, str], record_sequences: Dict[str, str], graph_path_coords: Dict[str, Coord], chrom_lengths: Dict[str, int], graph_mappings: Dict[str, str], ref_reader: Optional[FastaRegionReader], coord_by_name: Dict[str, Coord], allowchroms: Optional[Set[str]]) -> None:
    for pair in pairs:
        out_fh.write(build_provisional_row(pair, gcigars, graph_sequences, record_sequences, graph_path_coords, chrom_lengths, graph_mappings, ref_reader, coord_by_name, allowchroms=allowchroms) + '\n')

def _write_rows_parallel(pairs: Sequence[PairRow], out_fh, processes: int, chunksize: int, start_method: str, share_mode: str, gcigars: Dict[str, str], graph_sequences: Dict[str, str], record_sequences: Dict[str, str], graph_path_coords: Dict[str, Coord], chrom_lengths: Dict[str, int], graph_mappings: Dict[str, str], coord_by_name: Dict[str, Coord], allowchroms: Optional[Set[str]], reference_path: str, realignment_enabled: bool=True) -> None:
    ctx = mp.get_context(start_method) if start_method else mp.get_context()
    chunksize = chunksize if chunksize > 0 else _default_chunksize(len(pairs), processes)

    def run_pool(shared_state) -> None:
        with ctx.Pool(processes=processes, initializer=_init_worker_state, initargs=shared_state) as pool:
            for _idx, row, err in pool.imap(_worker_build_provisional_row, enumerate(pairs), chunksize):
                if err is not None:
                    pool.terminate()
                    raise RuntimeError(err)
                out_fh.write(row + '\n')
    if share_mode == 'manager':
        with ctx.Manager() as manager:
            shared_state = (manager.dict(gcigars), manager.dict(graph_sequences), manager.dict(record_sequences), manager.dict(graph_path_coords), manager.dict(chrom_lengths), manager.dict(graph_mappings), manager.dict(coord_by_name), tuple(sorted(allowchroms or [])), reference_path, realignment_enabled)
            run_pool(shared_state)
    else:
        shared_state = (gcigars, graph_sequences, record_sequences, graph_path_coords, chrom_lengths, graph_mappings, coord_by_name, tuple(sorted(allowchroms or [])), reference_path, realignment_enabled)
        run_pool(shared_state)

def main() -> None:
    ap = argparse.ArgumentParser(description='Production graphcigartoref with encoded-insertion liftover.')
    ap.add_argument('-i', '--input', required=True)
    ap.add_argument('-m', '--matchpair', required=True)
    ap.add_argument('-q', '--query', required=True)
    ap.add_argument('-g', '--graph', required=True)
    ap.add_argument('-r', '--reference', required=True)
    ap.add_argument('-o', '--output', default='')
    ap.add_argument('-u', '--unreplace', action='store_true', help='resume safely: read existing --output and skip pairs whose output columns 1, 3, and 5 are already present; append remaining rows instead of replacing the file')
    ap.add_argument('--query-header-alias-column', type=int, default=2)
    ap.add_argument('-t', '--processes', type=int, default=1, help='number of worker processes; 1 disables multiprocessing')
    ap.add_argument('--chunksize', type=int, default=0, help='Pool.imap chunksize; 0 chooses a balanced default')
    ap.add_argument('--mp-start-method', choices=('fork', 'spawn', 'forkserver'), default='fork', help='multiprocessing start method; default fork shares -q/-g dictionaries copy-on-write')
    ap.add_argument('--share-mode', choices=('fork', 'manager'), default='fork', help='manager shares large dictionaries through multiprocessing.Manager; fork uses ordinary globals/copy-on-write and is faster on Linux/macOS fork')
    ap.add_argument('--add-alt', action='store_true', help='keep alignments to non-reference/alternate chromosomes; by default these pieces are de-aligned to plain insertions with exact query payload')
    realignment_group = ap.add_mutually_exclusive_group()
    realignment_group.add_argument('--realignment', dest='realignment', action='store_true', default=True, help='enable local and score-linked SV-region realignment (default)')
    realignment_group.add_argument('--no-realignment', dest='realignment', action='store_false', help='disable local and score-linked SV-region realignment')
    ap.add_argument('--alternative', default='', help='four-column original-template BED; restrict VCF calls to projected intervals while preserving alignment context')
    args = ap.parse_args()
    alternative_policy = AlternativeIntervals(args.alternative) if args.alternative else None
    if alternative_policy and args.unreplace and args.output and os.path.exists(args.output):
        raise SystemExit('--alternative requires a fresh output; omit --unreplace to avoid reusing calls from a different interval policy')
    set_sv_realignment_enabled(args.realignment)
    if args.processes < 1:
        raise SystemExit('--processes must be >= 1')
    if args.chunksize < 0:
        raise SystemExit('--chunksize must be >= 0')
    if args.processes > 1 and args.share_mode == 'fork' and (args.mp_start_method != 'fork'):
        sys.stderr.write('[graphcigartoref] warning: --share-mode fork with a non-fork start method may copy the large -q sequence dictionary per worker; use --mp-start-method fork on Linux.\n')
    pairs = read_pairs(args.matchpair)
    gcigars = read_graphic_cigars(args.input)
    graph_sequences = read_graph_sequences(args.graph)
    record_sequences = read_fasta(args.query)
    graph_path_coords = read_graph_path_coords(args.graph)
    chrom_lengths = read_fai_lengths(args.reference)
    ref_reader = FastaRegionReader(args.reference)
    graph_mappings = read_graph_header_mappings(args.graph, ref_reader)
    aliases, coord_index, coord_by_name = read_fasta_header_aliases(args.query, args.query_header_alias_column)
    if any((pair.query_name is None or pair.ref_name is None for pair in pairs)):
        pairs = resolve_pairs_against_q_records(pairs, gcigars, aliases, coord_index)
    if alternative_policy is not None:
        alternative_policy.annotate_pairs(pairs, gcigars, graph_path_coords, coord_by_name)
    if args.unreplace:
        if not args.output:
            raise SystemExit('--unreplace requires --output so the existing save file can be scanned and appended')
        before_count = len(pairs)
        wanted_keys = {finished_pair_key(pair) for pair in pairs}
        finished_keys, old_rows, kept_rows, removed_rows, malformed_rows = prune_existing_output_for_pairs(args.output, wanted_keys)
        pairs = exclude_finished_pairs(pairs, finished_keys)
        skipped_count = before_count - len(pairs)
        sys.stderr.write(
            f'[graphcigartoref] --unreplace: cleaned {args.output}; '
            f'kept {kept_rows}/{old_rows} existing rows; removed {removed_rows}'
            f'{f" ({malformed_rows} malformed/blank)" if malformed_rows else ""}; '
            f'skipped {skipped_count}/{before_count} current pairs\n'
        )
    allowchroms = set() if args.add_alt else default_allowed_reference_chroms(pairs)
    output_mode = 'a' if args.unreplace else 'w'
    out_fh = open(args.output, output_mode) if args.output else sys.stdout
    if args.processes == 1:
        _write_rows_serial(pairs, out_fh, gcigars, graph_sequences, record_sequences, graph_path_coords, chrom_lengths, graph_mappings, ref_reader, coord_by_name, allowchroms)
    else:
        _write_rows_parallel(pairs, out_fh, args.processes, args.chunksize, args.mp_start_method, args.share_mode, gcigars, graph_sequences, record_sequences, graph_path_coords, chrom_lengths, graph_mappings, coord_by_name, allowchroms, args.reference, args.realignment)
    if out_fh is not sys.stdout:
        out_fh.close()
# --------------------------------------------------------------------------
# Opt-in stage profiler: set GRAPH_CIGARTOREF_PROFILE=1 to accumulate wall
# time per pipeline stage in every process (workers included) and print one
# summary line per process periodically and at exit.  Off by default; wraps
# call sites only when enabled, so normal runs are untouched.  Nested stages
# double-count deliberately (a child's time also appears in its parent).
_PROFILE_STAGE_NAMES = (
    'build_stage1',
    'build_stage1_core_first',
    'build_stage2',
    'build_stage2_core_seeded_extensions',
    'build_stage3_polished',
    'build_stage4_linear',
    '_align_graph_cigar_extension',
    'invert_path_to_ref_slice_for_compare',
    'slice_graphic_mapping_by_query',
    'compare_sliced_gcigars',
    'weighted_lcs',
    '_segment_local_cuts',
    'build_query_coverages',
    'encode_large_insertions_in_piece',
    '_rescue_internal_di_windows_multiscale',
    '_rescue_long_di_pairs_in_linear_piece',
    '_canonicalize_stage3_against_sequences',
    'slice_graph_cigar_by_query',
    'split_cigar_by_rposi',
    'parse_graphic_segments',
    'swspy',
    'minimap2_payload_ops',
    '_tm_global_insert_align_dp',
    'align_payload_ops',
)
_STAGE_TIMES: Dict[str, float] = {}
_STAGE_COUNTS: Dict[str, int] = {}
_PROFILE_ENABLED = (
    os.environ.get('GRAPH_CIGARTOREF_PROFILE', '') not in {'', '0'}
)


def _dump_stage_profile() -> None:
    total = sum(_STAGE_TIMES.values())
    if total <= 0:
        return
    ranked = sorted(_STAGE_TIMES.items(), key=lambda item: -item[1])
    sys.stderr.write(
        f'[graphcigartoref profile pid={os.getpid()}] '
        + '; '.join(
            f'{name}={seconds:.1f}s/{_STAGE_COUNTS[name]}'
            for name, seconds in ranked[:10]
        )
        + '\n'
    )


if os.environ.get('GRAPH_CIGARTOREF_PROFILE', '') not in {'', '0'}:
    import atexit as _atexit

    def _profile_wrap(fn, dump_every: int = 0):
        name = fn.__name__

        def wrapped(*args, **kwargs):
            begin = time.perf_counter()
            try:
                return fn(*args, **kwargs)
            finally:
                _STAGE_TIMES[name] = (
                    _STAGE_TIMES.get(name, 0.0) + time.perf_counter() - begin
                )
                count = _STAGE_COUNTS.get(name, 0) + 1
                _STAGE_COUNTS[name] = count
                if dump_every and count % dump_every == 0:
                    _dump_stage_profile()
        wrapped.__name__ = name
        return wrapped

    for _stage_name in _PROFILE_STAGE_NAMES:
        _stage_fn = globals().get(_stage_name)
        if callable(_stage_fn):
            globals()[_stage_name] = _profile_wrap(
                _stage_fn,
                dump_every=1024 if _stage_name == 'build_stage1' else 0,
            )
    _atexit.register(_dump_stage_profile)


if __name__ == '__main__':
    main()
