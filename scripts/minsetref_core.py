#!/usr/bin/env python3
"""Shared primitives for the min-set reference builder (v2).

Random-access FASTA, N-free segmentation, FASTA IO helpers, and a
multiprocessing context.  Interval/counter primitives live in
``minsetref_segments`` and are re-exported here for convenience.
"""
from __future__ import annotations

import dataclasses
import gzip
import multiprocessing as mp
import os
import re
import sys
from pathlib import Path
from typing import Dict, Iterator, List, Optional, Tuple

from minsetref_segments import (  # re-export
    complement_intervals,
    count_masked,
    count_unmasked,
    merge_intervals,
    subtract_intervals,
)

DNA_WRAP = 80
SAFE_ID_RE = re.compile(r"[^A-Za-z0-9_.:-]+")
SOURCE_ID_UNSAFE_RE = re.compile(r"[^A-Za-z0-9._#:-]+")
NFREE_RE = re.compile(r"[ACGTacgt]+")


def sanitize_id(text: str) -> str:
    text = SAFE_ID_RE.sub("_", str(text).strip()).strip("_")
    return text or "seq"


def source_interval_id(
    sample: str, contig: str, start: int, end: int,
) -> str:
    """Return a stable FASTA id that preserves its source coordinate interval."""
    if start < 0 or end <= start:
        raise ValueError(f"invalid source interval {contig}:{start}-{end}")
    safe_sample = SOURCE_ID_UNSAFE_RE.sub("_", str(sample).strip()).strip("_")
    safe_contig = SOURCE_ID_UNSAFE_RE.sub("_", str(contig).strip()).strip("_")
    if not safe_sample or not safe_contig:
        raise ValueError("sample and source contig must have non-empty FASTA ids")
    return (
        f"{safe_sample}_{safe_contig}_"
        f"{int(start)}_{int(end)}"
    )


def wrap_fasta(seq: str, width: int = DNA_WRAP) -> str:
    return "\n".join(seq[i:i + width] for i in range(0, len(seq), width))


def open_maybe_gzip(path: str, mode: str = "rt"):
    return gzip.open(path, mode) if str(path).endswith(".gz") else open(path, mode)


def read_primary_fasta_names(path: str) -> List[str]:
    names: List[str] = []
    with open_maybe_gzip(path, "rt") as fh:
        for raw in fh:
            if raw.startswith(">"):
                names.append(raw[1:].strip().split()[0])
    return names


def iter_n_free_segments(seq: str) -> Iterator[Tuple[int, int, str]]:
    """Yield maximal (start, end, subseq) runs of A/C/G/T/a/c/g/t (C-speed)."""
    for m in NFREE_RE.finditer(seq):
        yield m.start(), m.end(), m.group()


def mp_context():
    """spawn on macOS (fork is unsafe with the ObjC runtime), else fork."""
    if sys.platform.startswith("darwin"):
        return mp.get_context("spawn")
    try:
        return mp.get_context("fork")
    except (ValueError, RuntimeError):
        return mp.get_context()


class IndexedFasta:
    """Picklable random-access FASTA reader (os.pread), uncompressed only."""

    def __init__(
        self, fasta_path: str, fai_path: Optional[str] = None,
    ):
        self.path = str(fasta_path)
        if self.path.endswith(".gz"):
            raise ValueError(f"Random access requires an uncompressed FASTA: {self.path}")
        self.index: Dict[str, Tuple[int, int, int, int]] = {}
        fai = str(fai_path) if fai_path is not None else self.path + ".fai"
        if os.path.exists(fai):
            self._read_fai(fai)
        else:
            self._build_index()
        self._fd = os.open(self.path, os.O_RDONLY)

    def _read_fai(self, fai: str) -> None:
        with open(fai, "rt") as fh:
            for raw in fh:
                if not raw.strip():
                    continue
                cols = raw.rstrip("\n").split("\t")
                if len(cols) >= 5:
                    self.index[cols[0]] = tuple(map(int, cols[1:5]))  # type: ignore[assignment]

    def _build_index(self) -> None:
        name: Optional[str] = None
        length = 0
        seq_offset = 0
        line_bases = 0
        line_width = 0
        prev_short = False

        def finish() -> None:
            if name is not None:
                self.index[name] = (length, seq_offset, max(1, line_bases), max(1, line_width))

        with open(self.path, "rb") as fh:
            while True:
                line_start = fh.tell()
                raw = fh.readline()
                if not raw:
                    break
                if raw.startswith(b">"):
                    finish()
                    header = raw[1:].strip().split()
                    name = header[0].decode("utf-8") if header else ""
                    length = 0
                    seq_offset = fh.tell()
                    line_bases = 0
                    line_width = 0
                    prev_short = False
                    continue
                if name is None:
                    continue
                bases = raw.rstrip(b"\r\n")
                if line_bases == 0:
                    if not bases:
                        continue
                    line_bases = len(bases)
                    line_width = len(raw)
                    seq_offset = line_start
                    length += len(bases)
                    continue
                if prev_short or len(bases) > line_bases:
                    raise ValueError(
                        f"FASTA {self.path!r} is not uniformly wrapped at contig "
                        f"{name!r}: sequence line widths vary before the final line. "
                        f"Run `samtools faidx {self.path}` first (its .fai is trusted), "
                        f"or rewrap the FASTA to a fixed line width."
                    )
                if len(bases) < line_bases:
                    prev_short = True
                length += len(bases)
        finish()

    def __enter__(self):
        return self

    def __exit__(self, exc_type, exc, tb):
        self.close()

    def __getstate__(self):
        state = self.__dict__.copy()
        state["_fd"] = -1
        return state

    def __setstate__(self, state):
        self.__dict__.update(state)
        self._fd = os.open(self.path, os.O_RDONLY)

    def close(self) -> None:
        fd = getattr(self, "_fd", -1)
        if fd is not None and fd >= 0:
            os.close(fd)
            self._fd = -1

    def reopen(self) -> None:
        """Reopen a closed FASTA while retaining already loaded index metadata."""
        if getattr(self, "_fd", -1) < 0:
            self._fd = os.open(self.path, os.O_RDONLY)

    def _pread(self, n: int, offset: int) -> bytes:
        out = bytearray()
        while len(out) < n:
            chunk = os.pread(self._fd, n - len(out), offset + len(out))
            if not chunk:
                break
            out.extend(chunk)
        return bytes(out)

    def names(self) -> List[str]:
        return list(self.index.keys())

    def length(self, name: str) -> int:
        return self.index[name][0]

    def fetch(self, name: str, start: int, end: int) -> str:
        if name not in self.index:
            raise KeyError(f"{name!r} not found in FASTA {self.path}")
        length, offset, line_bases, line_width = self.index[name]
        start = max(0, min(int(start), length))
        end = max(start, min(int(end), length))
        out = bytearray()
        pos = start
        while pos < end:
            line_idx = pos // line_bases
            within = pos % line_bases
            n = min(end - pos, line_bases - within)
            out.extend(self._pread(n, offset + line_idx * line_width + within))
            pos += n
        return out.decode("ascii")

    def sequence(self, name: str) -> str:
        if name not in self.index:
            raise KeyError(f"{name!r} not found in FASTA {self.path}")
        length, offset, line_bases, line_width = self.index[name]
        if length == 0:
            return ""
        n_full = (length - 1) // line_bases
        last = length - n_full * line_bases
        total = n_full * line_width + last
        raw = self._pread(total, offset)
        if line_width != line_bases:
            raw = raw.replace(b"\r", b"").replace(b"\n", b"")
        return raw.decode("ascii")


def fasta_lengths(path: str) -> Dict[str, int]:
    with IndexedFasta(path) as fa:
        return {name: fa.length(name) for name in fa.names()}


@dataclasses.dataclass
class Piece:
    """One candidate/kept sequence with its source-assembly coordinates."""
    temp_id: str
    sample: str
    contig: str
    start: int
    end: int
    strand: str
    seq: str
    final_id: str = ""
    parent_id: str = ""
    note: str = ""
    # Stable identifier for the original candidate.  temp_id is deliberately
    # re-based during iterative self-cleaning, so it cannot be used for
    # provenance across cycles.
    origin_id: str = ""
    # Stable identity of one post-minimap residual. One original insertion can
    # split into residuals that lift to different local regions.
    lift_id: str = ""
    # Stable record order inside one haplotype. The sample's rank is deliberately
    # not stored; it is supplied dynamically from the current query list.
    within_sample_rank: int = -1

    def __post_init__(self) -> None:
        if not self.origin_id:
            self.origin_id = self.temp_id
        if not self.lift_id:
            self.lift_id = self.temp_id

    @property
    def length(self) -> int:
        return len(self.seq)

    @property
    def unmasked(self) -> int:
        return count_unmasked(self.seq)

    @property
    def masked(self) -> int:
        return count_masked(self.seq)

    def source_coord(self) -> str:
        return f"{self.sample}:{self.contig}:{self.start}-{self.end}{self.strand}"


def write_fasta_pieces(pieces, path: str, use_final_id: bool = False) -> None:
    Path(path).parent.mkdir(parents=True, exist_ok=True)
    with open(path, "wt") as out:
        for p in pieces:
            primary = p.final_id if use_final_id and p.final_id else p.temp_id
            out.write(f">{primary} {p.source_coord()}\n{wrap_fasta(p.seq)}\n")
