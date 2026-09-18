#!/usr/bin/env python3
"""Aligner wrappers and alignment parsing for the min-set reference builder (v2).

Stage 1 uses minimap2 for large-main and first-genome alignment. Stage 2 uses
Winnowmap for insertion self/pool alignment. Both emit ``-c --eqx`` PAF and
share the same parser. BLASTN is either a residual rescue or an independent
second pass, depending on the stage.
"""
from __future__ import annotations

import dataclasses
import heapq
import json
import logging
import os
import re
import shutil
import subprocess
import threading
from pathlib import Path
from typing import Dict, Iterator, List, Optional, Sequence, Tuple

from minsetref_core import (
    IndexedFasta,
    fasta_lengths,
    iter_n_free_segments,
    read_primary_fasta_names,
    mp_context,
    wrap_fasta,
)
from minsetref_segments import (
    complement_intervals,
    count_unmasked,
    merge_intervals,
    novelty_score,
)

QUERY_COVER_OPS = {"M", "=", "X"}
REF_CONSUME_OPS = {"M", "=", "X", "D", "N"}
CIGAR_RE = re.compile(r"(\d+)([MIDNSHP=X])")

# Default winnowmap parameters (repetitive-aware, fast); overridable.
WINNOWMAP_DEFAULTS = dict(k=19, w=10, m=300, p=0.001, N=100)
MINIMAP2_DEFAULTS = dict(preset="asm5", p=0.001, N=20, f=0.001, K="100M")
PAF_INDEX_THRESHOLD_BYTES = 100_000_000


# ---------------------------------------------------------------------------
# CIGAR helpers
# ---------------------------------------------------------------------------

def parse_cigar(cigar: str) -> List[Tuple[int, str]]:
    if not cigar or cigar == "*":
        return []
    return [(int(n), op) for n, op in CIGAR_RE.findall(cigar)]


def cigar_ref_span(cigar: str) -> int:
    return sum(n for n, op in parse_cigar(cigar) if op in REF_CONSUME_OPS)


def paired_coverage_and_blocks_from_paf_cigar(
    cigar: str, qstart: int, qend: int, tstart: int, strand: str,
    retain_blocks: bool = True,
) -> Tuple[
    List[Tuple[int, int]], List[Tuple[int, int]], int, int,
    List[Tuple[int, int, int, int]],
]:
    """Coverage statistics and paired blocks from one PAF CIGAR traversal.

    M/=/X consume and represent both sequences; query I/S and target D/N are
    gaps.  Reverse-strand query coordinates are mapped back to the forward query
    axis; target stays forward.
    """
    q_intervals: List[Tuple[int, int]] = []
    t_intervals: List[Tuple[int, int]] = []
    blocks: List[Tuple[int, int, int, int]] = []
    qpos = int(qstart) if strand != "-" else int(qend)
    tpos = int(tstart)
    aligned = 0
    matches = 0
    for n, op in parse_cigar(cigar):
        if op in QUERY_COVER_OPS:
            if strand == "-":
                q0, q1 = qpos - n, qpos
                qpos -= n
            else:
                q0, q1 = qpos, qpos + n
                qpos += n
            tq0, tq1 = min(q0, q1), max(q0, q1)
            target_block = (tpos, tpos + n)
            t_intervals.append(target_block)
            if retain_blocks:
                block = (tq0, tq1, target_block[0], target_block[1])
                if blocks:
                    pq0, pq1, pt0, pt1 = blocks[-1]
                    if strand == "-" and tq1 == pq0 and target_block[0] == pt1:
                        blocks[-1] = (tq0, pq1, pt0, target_block[1])
                    elif strand != "-" and pq1 == tq0 and pt1 == target_block[0]:
                        blocks[-1] = (pq0, tq1, pt0, target_block[1])
                    else:
                        blocks.append(block)
                else:
                    blocks.append(block)
            tpos += n
            q_intervals.append((tq0, tq1))
            aligned += n
            if op in {"=", "M"}:
                matches += n
        elif op in {"I", "S"}:
            qpos = qpos - n if strand == "-" else qpos + n
        elif op in {"D", "N"}:
            tpos += n
    return (
        merge_intervals(q_intervals),
        merge_intervals(t_intervals),
        aligned,
        matches,
        blocks,
    )


def paired_coverage_stats_from_paf_cigar(
    cigar: str, qstart: int, qend: int, tstart: int, strand: str
) -> Tuple[List[Tuple[int, int]], List[Tuple[int, int]], int, int]:
    q_intervals, t_intervals, aligned, matches, _blocks = (
        paired_coverage_and_blocks_from_paf_cigar(
            cigar, qstart, qend, tstart, strand,
        )
    )
    return q_intervals, t_intervals, aligned, matches


def paired_blocks_from_paf_cigar(cigar: str, qstart: int, qend: int, tstart: int,
                                 strand: str) -> List[Tuple[int, int, int, int]]:
    """Forward-axis M/=/X blocks (q0,q1,t0,t1) from a PAF cg CIGAR."""
    _query, _target, _aligned, _matches, blocks = (
        paired_coverage_and_blocks_from_paf_cigar(
            cigar, qstart, qend, tstart, strand,
        )
    )
    return blocks


def paired_blocks_from_sam_cigar(cigar: str, flag: int, qlen: Optional[int],
                                 target_start: int) -> List[Tuple[int, int, int, int]]:
    """Forward-axis M/=/X blocks (q0,q1,t0,t1) from a BLASTN outfmt-17 CIGAR."""
    blocks: List[Tuple[int, int, int, int]] = []
    qpos = 0
    tpos = max(0, int(target_start))
    for n, op in parse_cigar(cigar):
        if op in QUERY_COVER_OPS:
            blocks.append((qpos, qpos + n, tpos, tpos + n))
            qpos += n
            tpos += n
        elif op in {"I", "S", "H"}:
            qpos += n
        elif op in {"D", "N"}:
            tpos += n
    if (flag & 16) and qlen is not None:
        blocks = [(qlen - q1, qlen - q0, t0, t1) for q0, q1, t0, t1 in blocks]
    return blocks


def paired_coverage_from_sam_cigar(
    cigar: str, flag: int, qlen: Optional[int], target_start: int
) -> Tuple[List[Tuple[int, int]], List[Tuple[int, int]], int]:
    """Represented query/target blocks from BLASTN outfmt-17 SAM."""
    q_axis: List[Tuple[int, int]] = []
    t_intervals: List[Tuple[int, int]] = []
    qpos = 0
    tpos = max(0, int(target_start))
    aligned = 0
    for n, op in parse_cigar(cigar):
        if op in QUERY_COVER_OPS:
            q_axis.append((qpos, qpos + n))
            t_intervals.append((tpos, tpos + n))
            qpos += n
            tpos += n
            aligned += n
        elif op in {"I", "S", "H"}:
            qpos += n
        elif op in {"D", "N"}:
            tpos += n
    if qlen is None:
        qlen = qpos
    if flag & 16:
        q_intervals = [(max(0, qlen - b), max(0, qlen - a)) for a, b in q_axis]
    else:
        q_intervals = q_axis
    return merge_intervals(q_intervals), merge_intervals(t_intervals), aligned


def query_alignment_span_from_sam_cigar(cigar: str, flag: int = 0, qlen: Optional[int] = None) -> Tuple[int, int]:
    ops = parse_cigar(cigar)
    if not ops:
        return -1, -1
    consuming = {"M", "=", "X", "I", "S", "H"}
    total = sum(n for n, op in ops if op in consuming)
    if qlen is None:
        qlen = total
    left = 0
    for n, op in ops:
        if op in {"H", "S"}:
            left += n
        else:
            break
    right = 0
    for n, op in reversed(ops):
        if op in {"H", "S"}:
            right += n
        else:
            break
    a = max(0, left)
    b = max(a, total - right)
    if flag & 16:
        return max(0, qlen - b), max(0, qlen - a)
    return a, b


# ---------------------------------------------------------------------------
# Alignment record
# ---------------------------------------------------------------------------

@dataclasses.dataclass
class AlignmentHit:
    query_id: str
    target_id: str
    identity: float
    aligned_bases: int
    query_intervals: List[Tuple[int, int]]
    source: str
    target_start: int = -1
    target_end: int = -1
    query_start: int = -1
    query_end: int = -1
    strand: str = "+"
    target_intervals: List[Tuple[int, int]] = dataclasses.field(default_factory=list)
    # Paired M/=/X blocks as (q0, q1, t0, t1) on forward query and forward target
    # axes. For strand '-', within a block query q1<->target t0 and q0<->target t1.
    aligned_pairs: List[Tuple[int, int, int, int]] = dataclasses.field(default_factory=list)
    alignment_score: float = 0.0


def tags_from_paf(cols: Sequence[str]) -> Dict[str, str]:
    tags: Dict[str, str] = {}
    for tok in cols[12:]:
        parts = tok.split(":", 2)
        if len(parts) == 3:
            tags[parts[0]] = parts[2]
    return tags


def _parse_winnowmap_paf_line(
    raw: str,
    min_identity: float,
    min_align_size: int,
    source: str = "winnowmap",
    retain_aligned_pairs: bool = True,
) -> Optional[AlignmentHit]:
    if not raw.strip() or raw.startswith("#"):
        return None
    cols = raw.rstrip("\n").split("\t")
    if len(cols) < 12:
        return None
    qname = cols[0]
    qstart = int(cols[2])
    qend = int(cols[3])
    strand = cols[4] if len(cols) > 4 else "+"
    tname = cols[5]
    tstart = int(cols[7])
    tend = int(cols[8])
    nmatch = int(cols[9])
    block_len = int(cols[10])
    # The represented M/=/X length cannot exceed the PAF block length. Reject
    # a short row before decoding its potentially large cg:Z CIGAR.
    if block_len < max(1, min_align_size):
        return None
    paf_tags = tags_from_paf(cols)
    cigar = paf_tags.get("cg", "")
    if cigar:
        q_intervals, t_intervals, aligned, matches, pairs = (
            paired_coverage_and_blocks_from_paf_cigar(
                cigar, qstart=qstart, qend=qend, tstart=tstart,
                strand=strand, retain_blocks=retain_aligned_pairs,
            )
        )
        identity = 100.0 * matches / aligned if aligned else 0.0
    else:
        identity = 100.0 * nmatch / block_len
        matches = nmatch
        q_intervals = [(qstart, qend)]
        t_intervals = [(tstart, tend)]
        pairs = (
            [(qstart, qend, tstart, tend)]
            if retain_aligned_pairs else []
        )
        aligned = min(max(0, qend - qstart), max(0, tend - tstart))
    try:
        alignment_score = float(paf_tags.get("AS", matches))
    except ValueError:
        alignment_score = float(matches)
    if identity < min_identity or aligned < min_align_size:
        return None
    return AlignmentHit(
        query_id=qname, target_id=tname, identity=identity,
        aligned_bases=aligned, query_intervals=q_intervals, source=source,
        target_start=tstart, target_end=tend, query_start=qstart,
        query_end=qend, strand=strand, target_intervals=t_intervals,
        aligned_pairs=pairs, alignment_score=alignment_score,
    )


def parse_winnowmap_paf(
    path: str,
    min_identity: float,
    min_align_size: int,
    source: str = "winnowmap",
    retain_aligned_pairs: bool = True,
) -> Iterator[AlignmentHit]:
    """Parse winnowmap/minimap2 PAF (``-c --eqx``).

    Identity is computed on represented M/=/X blocks only, so query I/S and
    target D/N gaps remain candidate sequence.
    """
    if not os.path.exists(path):
        return
    with open(path, "rt") as fh:
        for raw in fh:
            hit = _parse_winnowmap_paf_line(
                raw, min_identity, min_align_size, source,
                retain_aligned_pairs,
            )
            if hit is not None:
                yield hit


# ---------------------------------------------------------------------------
# BLASTN outfmt-17 SAM
# ---------------------------------------------------------------------------

def load_blast_seq_names(seqfile: str) -> Dict[str, str]:
    mapping: Dict[str, str] = {}
    if not seqfile or not os.path.exists(seqfile):
        return mapping
    with open(seqfile, "rt") as fh:
        for i, raw in enumerate(fh):
            name = raw.strip().split()[0] if raw.strip() else ""
            if name:
                mapping[f"BL_ORD_ID:{i}"] = name
    return mapping


def parse_sam_tags(cols: Sequence[str]) -> Dict[str, str]:
    tags: Dict[str, str] = {}
    for tok in cols[11:]:
        parts = tok.split(":", 2)
        if len(parts) == 3:
            tags[parts[0]] = parts[2]
    return tags


def replace_blast_ord_name(text: str, ord_map: Dict[str, str]) -> str:
    return re.sub(r"BL_ORD_ID:\d+", lambda m: ord_map.get(m.group(0), m.group(0)), text)


def load_blast_query_names(
    query_lengths: Dict[str, int],
) -> Dict[str, str]:
    """Map BLAST's one-based Query_N aliases in FASTA record order."""
    return {
        f"Query_{index}": name
        for index, name in enumerate(query_lengths, 1)
    }


def parse_blast17_sam(path: str, seqfile: str, min_identity: float, min_align_size: int,
                      query_lengths: Optional[Dict[str, int]] = None) -> Iterator[AlignmentHit]:
    """Parse BLAST outfmt 17 while restoring its inverted SAM identifiers.

    BLAST writes database subjects in SAM column 1 (often ``BL_ORD_ID:N``)
    and input queries as the SAM reference in column 3 (often ``Query_N``).
    The public ``AlignmentHit`` contract is the opposite: ``query_id`` and
    query coordinates describe the input query FASTA, while ``target_id`` and
    target coordinates describe the database sequence.
    """
    ord_map = load_blast_seq_names(seqfile)
    query_lengths = query_lengths or {}
    query_name_map = load_blast_query_names(query_lengths)
    if not os.path.exists(path):
        return
    with open(path, "rt", errors="replace") as fh:
        for line_number, raw in enumerate(fh, 1):
            if not raw.strip() or raw.startswith("@") or raw.startswith("#"):
                continue
            cols = raw.rstrip("\n").split("\t")
            if len(cols) < 11:
                cols = raw.strip().split()
                if len(cols) < 11:
                    continue
            target_name = replace_blast_ord_name(cols[0], ord_map)
            query_name = query_name_map.get(cols[2], cols[2])
            if re.fullmatch(r"BL_ORD_ID:\d+", target_name):
                raise ValueError(
                    f"{path}:{line_number}: cannot restore database subject "
                    f"{target_name!r}; the BLAST database .seq map is missing "
                    "or incomplete"
                )
            if re.fullmatch(r"Query_\d+", query_name) and query_name not in query_lengths:
                raise ValueError(
                    f"{path}:{line_number}: cannot restore BLAST query alias "
                    f"{query_name!r}; original query FASTA names are unavailable"
                )
            try:
                flag = int(cols[1])
            except Exception:
                flag = 0
            if flag & 4:
                continue
            if query_name == "*" or target_name == "*":
                continue
            try:
                query_reference_start = max(0, int(cols[3]) - 1)
            except Exception:
                query_reference_start = -1
            cigar = cols[5]
            if cigar == "*":
                continue
            if query_reference_start < 0:
                continue
            target_intervals, query_intervals, aligned = paired_coverage_from_sam_cigar(
                cigar, flag=flag, qlen=None,
                target_start=query_reference_start,
            )
            if aligned < min_align_size:
                continue
            tags = parse_sam_tags(cols)
            identity: Optional[float] = None
            if "PI" in tags:
                try:
                    identity = float(tags["PI"])
                except ValueError:
                    identity = None
            if identity is None and "NM" in tags:
                try:
                    nm = int(float(tags["NM"]))
                    identity = 100.0 * max(0, aligned - nm) / aligned if aligned else 0.0
                except ValueError:
                    identity = None
            if identity is None:
                # No usable PI/NM: identity cannot be validated -> ignore hit.
                continue
            if identity < min_identity:
                continue
            target_start, target_end = query_alignment_span_from_sam_cigar(
                cigar, flag=flag, qlen=None,
            )
            query_end = query_reference_start + cigar_ref_span(cigar)
            sam_pairs = paired_blocks_from_sam_cigar(
                cigar, flag, None, query_reference_start,
            )
            pairs = [
                (query_start, query_stop, target_start_block, target_stop)
                for target_start_block, target_stop, query_start, query_stop
                in sam_pairs
            ]
            try:
                alignment_score = float(tags.get("AS", aligned * identity / 100.0))
            except ValueError:
                alignment_score = float(aligned) * identity / 100.0
            yield AlignmentHit(
                query_id=query_name, target_id=target_name,
                identity=identity, aligned_bases=aligned,
                query_intervals=query_intervals, source="blastn",
                target_start=target_start, target_end=target_end,
                query_start=query_reference_start, query_end=query_end,
                strand="-" if flag & 16 else "+",
                target_intervals=target_intervals,
                aligned_pairs=pairs,
                alignment_score=alignment_score,
            )


# ---------------------------------------------------------------------------
# External command wrappers
# ---------------------------------------------------------------------------

def run_command(cmd: Sequence[str], log: logging.Logger, stdout_path: Optional[str] = None) -> None:
    log.info("RUN: %s", " ".join(map(str, cmd)))
    if stdout_path:
        temporary = f"{stdout_path}.tmp.{os.getpid()}.{threading.get_ident()}"
        try:
            with open(temporary, "wt") as out:
                proc = subprocess.run(
                    cmd, stdout=out, stderr=subprocess.PIPE, text=True,
                )
            if proc.returncode == 0:
                os.replace(temporary, stdout_path)
        finally:
            try:
                os.remove(temporary)
            except FileNotFoundError:
                pass
    else:
        proc = subprocess.run(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)
    if proc.stdout:
        log.debug(proc.stdout.rstrip())
    if proc.stderr:
        log.debug(proc.stderr.rstrip())
    if proc.returncode != 0:
        raise subprocess.CalledProcessError(proc.returncode, cmd, output=proc.stdout, stderr=proc.stderr)


def ensure_executable(name: str) -> None:
    if shutil.which(name) is None:
        raise RuntimeError(f"Required executable not found on PATH: {name}")


def build_blast_db(ref_fasta: str, db_prefix: str, log: logging.Logger) -> str:
    dbp = Path(db_prefix)
    dbp.parent.mkdir(parents=True, exist_ok=True)
    run_command(["makeblastdb", "-in", ref_fasta, "-dbtype", "nucl",
                 "-blastdb_version", "4", "-out", str(dbp)], log)
    seq_path = str(dbp) + ".seq"
    with open(seq_path, "wt") as out:
        for name in read_primary_fasta_names(ref_fasta):
            out.write(name + "\n")
    return str(dbp)


def build_winnowmap_rep_kmers(ref_fasta: str, out_txt: str, k: int, log: logging.Logger,
                              distinct: float = 0.9998) -> Optional[str]:
    """Build one reusable Winnowmap repetitive-k-mer file via Meryl.

    Returns the path, or None when meryl is unavailable (winnowmap then runs
    without -W, i.e. as a repetitive-unaware minimap2).
    """
    Path(out_txt).parent.mkdir(parents=True, exist_ok=True)
    if (os.path.isfile(out_txt) and os.path.getsize(out_txt) > 0
            and os.path.getmtime(out_txt) >= os.path.getmtime(ref_fasta)):
        log.info("REUSE: repetitive-k-mer file %s", out_txt)
        return out_txt
    if shutil.which("meryl") is None:
        log.warning("meryl not found; winnowmap will run without a repetitive-k-mer db (-W)")
        return None
    db_dir = out_txt + ".meryl"
    if os.path.exists(db_dir):
        shutil.rmtree(db_dir, ignore_errors=True)
    try:
        run_command(["meryl", "count", f"k={k}", "output", db_dir, ref_fasta], log)
        run_command(["meryl", "print", "greater-than", f"distinct={distinct}", db_dir], log,
                    stdout_path=out_txt)
    finally:
        # Winnowmap only consumes the text list. Keeping the much larger Meryl
        # database would add unnecessary cohort disk usage.
        shutil.rmtree(db_dir, ignore_errors=True)
    return out_txt


def run_winnowmap(query_fasta: str, target_fasta: str, out_paf: str, threads: int,
                  log: logging.Logger, params: Optional[dict] = None,
                  rep_kmers: Optional[str] = None, reuse_existing: bool = False) -> None:
    if reuse_existing and os.path.exists(out_paf):
        log.info("REUSE: existing winnowmap PAF %s", out_paf)
        return
    p = dict(WINNOWMAP_DEFAULTS)
    if params:
        p.update(params)
    cmd = ["winnowmap", "-c", "--eqx",
           "-k", str(p["k"]), "-w", str(p["w"]), "-m", str(p["m"]), "-p", str(p["p"]),
           "-N", str(p["N"]), "-t", str(threads)]
    if rep_kmers and os.path.exists(rep_kmers):
        cmd += ["-W", rep_kmers]
    cmd += [target_fasta, query_fasta]
    run_command(cmd, log, stdout_path=out_paf)


def run_minimap2(query_fasta: str, target_fasta: str, out_paf: str, threads: int,
                 log: logging.Logger, params: Optional[dict] = None,
                 reuse_existing: bool = False,
                 emit_cigar: bool = True) -> None:
    """Run the fast stage-1 assembly aligner using the older pipeline policy."""
    if reuse_existing and os.path.exists(out_paf):
        log.info("REUSE: existing minimap2 PAF %s", out_paf)
        return
    p = dict(MINIMAP2_DEFAULTS)
    if params:
        p.update(params)

    def command(worker_threads: int, query_batch: Optional[str]) -> List[str]:
        cmd = ["minimap2", "-x", str(p["preset"])]
        if emit_cigar:
            cmd.extend(["-c", "--eqx"])
        cmd.extend([
            "-N", str(p["N"]), "-p", str(p["p"]), "-f", str(p["f"]),
        ])
        if query_batch:
            cmd.extend(["-K", str(query_batch)])
        cmd.extend([
            "-t", str(worker_threads),
            target_fasta, query_fasta,
        ])
        return cmd

    try:
        run_command(
            command(threads, p.get("K")),
            log,
            stdout_path=out_paf,
        )
    except subprocess.CalledProcessError as error:
        if (
            error.returncode not in {-9, 137}
            or not bool(p.get("retry_on_sigkill"))
        ):
            raise
        retry_threads = max(
            1, min(int(p.get("retry_threads", 32)), threads),
        )
        retry_batch = str(p.get("retry_K", "100M"))
        log.warning(
            "Minimap2 was killed; retrying once with %d threads and "
            "-K %s to reduce peak RAM",
            retry_threads,
            retry_batch,
        )
        try:
            os.remove(out_paf)
        except FileNotFoundError:
            pass
        run_command(
            command(retry_threads, retry_batch),
            log,
            stdout_path=out_paf,
        )


def build_minimap2_index(ref_fasta: str, out_index: str, log: logging.Logger,
                         preset: str = "asm5") -> str:
    """Build one reusable minimap2 index for parallel assembly jobs."""
    Path(out_index).parent.mkdir(parents=True, exist_ok=True)
    if os.path.isfile(out_index) and os.path.getsize(out_index) > 0:
        if os.path.getmtime(out_index) >= os.path.getmtime(ref_fasta):
            log.info("REUSE: minimap2 index %s", out_index)
            return out_index
    temporary = out_index + ".tmp"
    try:
        os.remove(temporary)
    except FileNotFoundError:
        pass
    run_command(["minimap2", "-x", str(preset), "-d", temporary, ref_fasta], log)
    if not os.path.isfile(temporary) or os.path.getsize(temporary) == 0:
        raise RuntimeError(f"minimap2 produced no index: {temporary}")
    os.replace(temporary, out_index)
    return out_index


def run_blastn17(query_fasta: str, db_prefix: str, out_sam: str, threads: int, min_identity: float,
                 word_size: int, evalue: str, max_target_seqs: int, log: logging.Logger,
                 reuse_existing: bool = False) -> None:
    if not os.path.exists(query_fasta) or os.path.getsize(query_fasta) == 0:
        raise RuntimeError(f"BLASTN residual query FASTA missing/empty: {query_fasta}")
    if reuse_existing and os.path.exists(out_sam):
        log.info("REUSE: completed BLASTN SAM %s", out_sam)
        return
    run_command([
        "blastn", "-task", "megablast", "-query", query_fasta, "-db", db_prefix,
        "-outfmt", "17", "-word_size", str(word_size), "-num_threads", str(threads),
        "-evalue", str(evalue), "-dust", "yes", "-lcase_masking",
        "-perc_identity", str(min_identity), "-qcov_hsp_perc", "0.1",
        "-max_target_seqs", str(max_target_seqs), "-parse_deflines",
        "-out", out_sam,
    ], log)


# ---------------------------------------------------------------------------
# Coverage + residual selection
# ---------------------------------------------------------------------------

def coverage_by_query_id(hits: Sequence[AlignmentHit], query_lengths: Dict[str, int]) -> Dict[str, List[Tuple[int, int]]]:
    """Union of represented query intervals per query id (clamped to length)."""
    cov: Dict[str, List[Tuple[int, int]]] = {}
    for hit in hits:
        qlen = query_lengths.get(hit.query_id)
        if qlen is None:
            continue
        acc = cov.setdefault(hit.query_id, [])
        for a, b in hit.query_intervals:
            a = max(0, min(int(a), qlen))
            b = max(0, min(int(b), qlen))
            if b > a:
                acc.append((a, b))
    return {qid: merge_intervals(vals) for qid, vals in cov.items()}


_FILTER_HITS: Sequence[AlignmentHit] = ()
_FILTER_HIT_INDICES: Dict[str, List[int]] = {}
_FILTER_SEQUENCES: Dict[str, str] = {}
_FILTER_MIN_UNMASKED = 0
_FILTER_MASKED_WEIGHT = 0.0
_FILTER_THRESHOLDS: Dict[str, float] = {}


def _filter_one_query_hit_indices(qid: str) -> List[int]:
    return _filter_one_query_hit_indices_from(
        qid,
        _FILTER_HITS,
        _FILTER_HIT_INDICES,
        _FILTER_SEQUENCES,
        _FILTER_MIN_UNMASKED,
        _FILTER_MASKED_WEIGHT,
        _FILTER_THRESHOLDS,
    )


def _filter_one_query_hit_indices_from(
    qid: str,
    hits: Sequence[AlignmentHit],
    hit_indices: Dict[str, List[int]],
    sequences: Dict[str, str],
    min_unmasked: int,
    masked_weight: float,
    thresholds: Dict[str, float],
) -> List[int]:
    indices = hit_indices.get(qid, ())
    seq = sequences.get(qid)
    if seq is None:
        return []
    qlen = len(seq)
    query_threshold = float(
        thresholds.get(qid, min_unmasked)
    )
    all_intervals: List[Tuple[int, int]] = []
    endpoints = set()
    for index in indices:
        for raw_a, raw_b in hits[index].query_intervals:
            a = max(0, min(int(raw_a), qlen))
            b = max(0, min(int(raw_b), qlen))
            a, b = min(a, b), max(a, b)
            if b <= a:
                continue
            all_intervals.append((a, b))
            endpoints.add(a)
            endpoints.add(b)
    union = merge_intervals(all_intervals)
    prefix: Dict[int, float] = {}
    total = 0.0
    ordered_endpoints = sorted(endpoints)
    if ordered_endpoints:
        prefix[ordered_endpoints[0]] = 0.0
    union_index = 0
    for left, right in zip(ordered_endpoints, ordered_endpoints[1:]):
        while union_index < len(union) and union[union_index][1] <= left:
            union_index += 1
        if (
            union_index < len(union)
            and union[union_index][0] <= left
            and right <= union[union_index][1]
        ):
            segment = seq[left:right]
            total += (
                novelty_score(segment, masked_weight)
                if masked_weight
                else count_unmasked(segment)
            )
        prefix[right] = total
    kept: List[int] = []
    for index in indices:
        hit = hits[index]
        if hit.aligned_bases < query_threshold:
            continue
        intervals = merge_intervals([
            (max(0, min(int(a), qlen)), max(0, min(int(b), qlen)))
            for a, b in hit.query_intervals
        ])
        if sum(prefix[b] - prefix[a] for a, b in intervals) >= query_threshold:
            kept.append(index)
    return kept


def _filter_query_chunk(qids: Sequence[str]) -> List[int]:
    kept: List[int] = []
    for qid in qids:
        kept.extend(_filter_one_query_hit_indices(qid))
    return kept


def filter_hits_by_unmasked_query(hits: Sequence[AlignmentHit], query_fasta: str,
                                  min_unmasked: int,
                                  masked_weight: float = 0.0,
                                  min_score_by_query: Optional[Dict[str, float]] = None,
                                  workers: int = 1,
                                  ) -> List[AlignmentHit]:
    """Keep hits whose aligned query sequence reaches the novelty threshold.

    Score is uppercase A/C/G/T plus ``masked_weight`` times lowercase a/c/g/t.
    The default weight of zero preserves the original unmasked-only rule.
    Intervals are unioned before scoring so adjacent CIGAR blocks cannot double
    count query sequence.
    """
    if min_unmasked <= 0:
        return list(hits)
    by_query: Dict[str, List[int]] = {}
    for index, hit in enumerate(hits):
        by_query.setdefault(hit.query_id, []).append(index)
    with IndexedFasta(query_fasta) as fa:
        available = set(fa.names())
        sequences = {
            qid: fa.sequence(qid)
            for qid in by_query
            if qid in available
        }

    qids = list(sequences)
    worker_count = min(max(1, int(workers)), max(1, len(qids)))
    context = mp_context()
    use_processes = (
        worker_count > 1
        and context.get_start_method() == "fork"
        and len(qids) >= worker_count * 4
        and threading.current_thread() is threading.main_thread()
    )
    if use_processes:
        global _FILTER_HITS
        global _FILTER_HIT_INDICES
        global _FILTER_SEQUENCES
        global _FILTER_MIN_UNMASKED
        global _FILTER_MASKED_WEIGHT
        global _FILTER_THRESHOLDS
        _FILTER_HITS = hits
        _FILTER_HIT_INDICES = by_query
        _FILTER_SEQUENCES = sequences
        _FILTER_MIN_UNMASKED = min_unmasked
        _FILTER_MASKED_WEIGHT = masked_weight
        _FILTER_THRESHOLDS = min_score_by_query or {}
        try:
            # Greedily balance by hit count. Workers inherit the read-only hit
            # and sequence tables through fork; only integer indices return.
            buckets: List[Tuple[int, int, List[str]]] = [
                (0, index, []) for index in range(worker_count)
            ]
            heapq.heapify(buckets)
            for qid in sorted(
                qids, key=lambda name: -len(by_query[name]),
            ):
                load, bucket_index, members = heapq.heappop(buckets)
                members.append(qid)
                heapq.heappush(buckets, (
                    load + len(by_query[qid]), bucket_index, members,
                ))
            chunks = [
                members for _load, _index, members in sorted(
                    buckets, key=lambda row: row[1],
                ) if members
            ]
            with context.Pool(processes=len(chunks)) as pool:
                selected_chunks = pool.map(_filter_query_chunk, chunks)
            selected_indices = sorted(
                index
                for chunk in selected_chunks
                for index in chunk
            )
        finally:
            _FILTER_HITS = ()
            _FILTER_HIT_INDICES = {}
            _FILTER_SEQUENCES = {}
            _FILTER_THRESHOLDS = {}
    else:
        thresholds = min_score_by_query or {}
        selected_indices = sorted(
            index
            for qid in qids
            for index in _filter_one_query_hit_indices_from(
                qid,
                hits,
                by_query,
                sequences,
                min_unmasked,
                masked_weight,
                thresholds,
            )
        )
    return [hits[index] for index in selected_indices]


PafQueryIndex = List[Tuple[str, List[Tuple[int, int]]]]


# Per-process state for indexed large-PAF filtering.  A process pool is used
# instead of Python threads because CIGAR parsing and interval filtering are
# CPU-heavy Python work.  Each worker owns its PAF descriptor and FASTA reader,
# avoiding shared seek state while retaining copy-on-write configuration data.
_PAF_WORKER_HANDLE = None
_PAF_WORKER_SEQUENCES: Optional[IndexedFasta] = None
_PAF_WORKER_AVAILABLE: set[str] = set()
_PAF_WORKER_MIN_IDENTITY = 0.0
_PAF_WORKER_PARSER_MINIMUM = 0
_PAF_WORKER_MIN_SEGMENT = 0
_PAF_WORKER_MASKED_WEIGHT = 0.0
_PAF_WORKER_THRESHOLDS: Dict[str, float] = {}
_PAF_WORKER_SOURCE = "winnowmap"
_PAF_WORKER_RETAIN_PAIRS = True
_PAF_WORKER_AGGREGATE = False
_PAF_WORKER_INTERLEAVED = False


def _initialize_paf_process_worker(
    paf_path: str,
    query_fasta: str,
    min_identity: float,
    parser_minimum: int,
    min_segment: int,
    masked_weight: float,
    thresholds: Dict[str, float],
    source: str,
    retain_aligned_pairs: bool,
    aggregate_query_coverage: bool,
    interleaved: bool,
) -> None:
    global _PAF_WORKER_HANDLE
    global _PAF_WORKER_SEQUENCES
    global _PAF_WORKER_AVAILABLE
    global _PAF_WORKER_MIN_IDENTITY
    global _PAF_WORKER_PARSER_MINIMUM
    global _PAF_WORKER_MIN_SEGMENT
    global _PAF_WORKER_MASKED_WEIGHT
    global _PAF_WORKER_THRESHOLDS
    global _PAF_WORKER_SOURCE
    global _PAF_WORKER_RETAIN_PAIRS
    global _PAF_WORKER_AGGREGATE
    global _PAF_WORKER_INTERLEAVED
    _PAF_WORKER_HANDLE = open(paf_path, "rb")
    # Linux/fork workers inherit the already parsed FASTA index copy-on-write.
    # Spawn workers (macOS) construct it here instead.
    if (
        _PAF_WORKER_SEQUENCES is None
        or _PAF_WORKER_SEQUENCES.path != query_fasta
    ):
        _PAF_WORKER_SEQUENCES = IndexedFasta(query_fasta)
        _PAF_WORKER_AVAILABLE = set(_PAF_WORKER_SEQUENCES.names())
    else:
        _PAF_WORKER_SEQUENCES.reopen()
    _PAF_WORKER_MIN_IDENTITY = min_identity
    _PAF_WORKER_PARSER_MINIMUM = parser_minimum
    _PAF_WORKER_MIN_SEGMENT = min_segment
    _PAF_WORKER_MASKED_WEIGHT = masked_weight
    _PAF_WORKER_THRESHOLDS = thresholds
    _PAF_WORKER_SOURCE = source
    _PAF_WORKER_RETAIN_PAIRS = retain_aligned_pairs
    _PAF_WORKER_AGGREGATE = aggregate_query_coverage
    _PAF_WORKER_INTERLEAVED = interleaved


def _filter_indexed_paf_process_chunk(
    task: Tuple[int, Sequence[Tuple[str, List[Tuple[int, int]]]]],
) -> Tuple[int, int, List[AlignmentHit], List[Tuple[int, AlignmentHit]]]:
    """Filter one byte-local query chunk inside a persistent worker process."""
    chunk_index, entries = task
    paf = _PAF_WORKER_HANDLE
    sequences = _PAF_WORKER_SEQUENCES
    if paf is None or sequences is None:
        raise RuntimeError("indexed PAF process worker was not initialized")
    retained: List[AlignmentHit] = []
    retained_with_offsets: List[Tuple[int, AlignmentHit]] = []
    for qid, ranges in entries:
        if qid not in _PAF_WORKER_AVAILABLE:
            continue
        query_hits: List[AlignmentHit] = []
        query_offsets: List[int] = []
        for range_start, range_end in ranges:
            paf.seek(range_start)
            while paf.tell() < range_end:
                offset = paf.tell()
                raw = paf.readline()
                if not raw:
                    raise IOError(
                        f"unexpected EOF in indexed PAF range "
                        f"{range_start}-{range_end}"
                    )
                hit = _parse_winnowmap_paf_line(
                    raw.decode("utf-8"),
                    _PAF_WORKER_MIN_IDENTITY,
                    _PAF_WORKER_PARSER_MINIMUM,
                    _PAF_WORKER_SOURCE,
                    _PAF_WORKER_RETAIN_PAIRS,
                )
                if hit is None:
                    continue
                if hit.query_id != qid:
                    raise ValueError(
                        f"PAF query index expected {qid!r} at byte {offset}, "
                        f"found {hit.query_id!r}"
                    )
                query_offsets.append(offset)
                query_hits.append(hit)
        if not query_hits:
            continue
        if _PAF_WORKER_MIN_SEGMENT <= 0:
            selected_indices = range(len(query_hits))
        else:
            selected_indices = _filter_one_query_hit_indices_from(
                qid,
                query_hits,
                {qid: list(range(len(query_hits)))},
                {qid: sequences.sequence(qid)},
                _PAF_WORKER_MIN_SEGMENT,
                _PAF_WORKER_MASKED_WEIGHT,
                _PAF_WORKER_THRESHOLDS,
            )
        if _PAF_WORKER_AGGREGATE:
            intervals = merge_intervals([
                interval
                for index in selected_indices
                for interval in query_hits[index].query_intervals
            ])
            if intervals:
                retained.append(AlignmentHit(
                    query_id=qid,
                    target_id="*",
                    identity=100.0,
                    aligned_bases=sum(end - start for start, end in intervals),
                    query_intervals=intervals,
                    source=f"{_PAF_WORKER_SOURCE}_query_coverage",
                ))
        elif _PAF_WORKER_INTERLEAVED:
            retained_with_offsets.extend(
                (query_offsets[index], query_hits[index])
                for index in selected_indices
            )
        else:
            retained.extend(query_hits[index] for index in selected_indices)
    return chunk_index, len(entries), retained, retained_with_offsets


def _read_paf_query_index(path: str) -> PafQueryIndex:
    by_query: Dict[str, List[Tuple[int, int]]] = {}
    order: List[str] = []
    with open(path, "rt") as handle:
        for line_number, raw in enumerate(handle, 1):
            fields = raw.rstrip("\n").split("\t")
            if len(fields) != 3:
                raise ValueError(f"{path}:{line_number}: expected QUERY START END")
            qid, start_text, end_text = fields
            start, end = int(start_text), int(end_text)
            if not qid or start < 0 or end <= start:
                raise ValueError(f"{path}:{line_number}: invalid PAF byte range")
            if qid not in by_query:
                by_query[qid] = []
                order.append(qid)
            by_query[qid].append((start, end))
    return [(qid, by_query[qid]) for qid in order]


def index_paf_by_query(
    paf_path: str,
    log: logging.Logger,
) -> PafQueryIndex:
    """Persist contiguous byte ranges for each PAF query and return the index."""
    if not os.path.isfile(paf_path):
        return []
    index_path = paf_path + ".query_ranges.tsv"
    metadata_path = index_path + ".meta.json"
    paf_stat = os.stat(paf_path)
    signature = {
        "size": paf_stat.st_size,
        "mtime_ns": paf_stat.st_mtime_ns,
    }
    try:
        with open(metadata_path, "rt") as handle:
            metadata = json.load(handle)
        if metadata == signature and os.path.isfile(index_path):
            indexed = _read_paf_query_index(index_path)
            log.info(
                "REUSE: PAF query-range index %s (%d queries)",
                index_path, len(indexed),
            )
            return indexed
    except (OSError, ValueError, TypeError, json.JSONDecodeError):
        pass

    log.info("Indexing PAF byte ranges by query: %s", paf_path)
    by_query: Dict[str, List[Tuple[int, int]]] = {}
    order: List[str] = []
    current_query: Optional[str] = None
    current_start = 0
    with open(paf_path, "rb") as paf:
        while True:
            line_start = paf.tell()
            raw = paf.readline()
            if not raw:
                if current_query is not None:
                    by_query[current_query].append((current_start, line_start))
                break
            stripped = raw.strip()
            tab = raw.find(b"\t")
            query = None
            if stripped and not raw.startswith(b"#") and tab > 0:
                query = raw[:tab].decode("utf-8")
            if query != current_query:
                if current_query is not None:
                    by_query[current_query].append((current_start, line_start))
                current_query = query
                if query is not None:
                    current_start = line_start
                    if query not in by_query:
                        by_query[query] = []
                        order.append(query)

    temporary_index = f"{index_path}.tmp.{os.getpid()}"
    temporary_metadata = f"{metadata_path}.tmp.{os.getpid()}"
    try:
        with open(temporary_index, "wt") as output:
            for qid in order:
                for start, end in by_query[qid]:
                    output.write(f"{qid}\t{start}\t{end}\n")
        with open(temporary_metadata, "wt") as output:
            json.dump(signature, output, sort_keys=True)
            output.write("\n")
        os.replace(temporary_index, index_path)
        os.replace(temporary_metadata, metadata_path)
    finally:
        for temporary in (temporary_index, temporary_metadata):
            try:
                os.remove(temporary)
            except FileNotFoundError:
                pass
    indexed = [(qid, by_query[qid]) for qid in order]
    log.info(
        "Indexed %d PAF queries in %d contiguous byte ranges",
        len(indexed), sum(len(ranges) for _qid, ranges in indexed),
    )
    return indexed


def load_filtered_paf_hits_indexed(
    paf_path: str,
    query_fasta: str,
    min_identity: float,
    parser_minimum: int,
    min_segment: int,
    masked_weight: float,
    min_score_by_query: Optional[Dict[str, float]],
    log: logging.Logger,
    *,
    source: str,
    retain_aligned_pairs: bool,
    aggregate_query_coverage: bool = False,
    workers: int = 1,
) -> List[AlignmentHit]:
    """Parse/filter one indexed PAF query at a time with bounded raw-hit RAM."""
    query_index = index_paf_by_query(paf_path, log)
    if not query_index:
        return []
    interleaved = (
        not aggregate_query_coverage
        and any(len(ranges) > 1 for _qid, ranges in query_index)
    )
    thresholds = min_score_by_query or {}
    progress_step = max(1, len(query_index) // 20)
    sequences = IndexedFasta(query_fasta)
    available = set(sequences.names())

    def load_chunk(
        entries: Sequence[Tuple[str, List[Tuple[int, int]]]],
    ) -> Tuple[List[AlignmentHit], List[Tuple[int, AlignmentHit]]]:
        chunk_retained: List[AlignmentHit] = []
        chunk_offsets: List[Tuple[int, AlignmentHit]] = []
        with open(paf_path, "rb") as paf:
            for qid, ranges in entries:
                if qid not in available:
                    continue
                query_hits: List[AlignmentHit] = []
                query_offsets: List[int] = []
                for range_start, range_end in ranges:
                    paf.seek(range_start)
                    while paf.tell() < range_end:
                        offset = paf.tell()
                        raw = paf.readline()
                        if not raw:
                            raise IOError(
                                f"unexpected EOF in indexed PAF range "
                                f"{paf_path}:{range_start}-{range_end}"
                            )
                        hit = _parse_winnowmap_paf_line(
                            raw.decode("utf-8"), min_identity, parser_minimum,
                            source, retain_aligned_pairs,
                        )
                        if hit is None:
                            continue
                        if hit.query_id != qid:
                            raise ValueError(
                                f"PAF query index expected {qid!r} at byte "
                                f"{offset}, found {hit.query_id!r}"
                            )
                        query_offsets.append(offset)
                        query_hits.append(hit)
                if not query_hits:
                    continue
                if min_segment <= 0:
                    selected_indices = range(len(query_hits))
                else:
                    selected_indices = _filter_one_query_hit_indices_from(
                        qid,
                        query_hits,
                        {qid: list(range(len(query_hits)))},
                        {qid: sequences.sequence(qid)},
                        min_segment,
                        masked_weight,
                        thresholds,
                    )
                if aggregate_query_coverage:
                    intervals = merge_intervals([
                        interval
                        for index in selected_indices
                        for interval in query_hits[index].query_intervals
                    ])
                    if intervals:
                        chunk_retained.append(AlignmentHit(
                            query_id=qid,
                            target_id="*",
                            identity=100.0,
                            aligned_bases=sum(
                                end - start for start, end in intervals
                            ),
                            query_intervals=intervals,
                            source=f"{source}_query_coverage",
                        ))
                elif interleaved:
                    chunk_offsets.extend(
                        (query_offsets[index], query_hits[index])
                        for index in selected_indices
                    )
                else:
                    chunk_retained.extend(
                        query_hits[index] for index in selected_indices
                    )
        return chunk_retained, chunk_offsets

    worker_count = min(max(1, int(workers)), len(query_index))
    try:
        if worker_count > 1:
            # Query-index order follows first appearance in the PAF, so these
            # chunks are mostly byte-contiguous. Separate processes provide
            # real CPU parallelism for Python CIGAR parsing and filtering.
            chunk_size = max(1, (
                len(query_index) + worker_count * 8 - 1
            ) // (worker_count * 8))
            chunks = [
                query_index[start:start + chunk_size]
                for start in range(0, len(query_index), chunk_size)
            ]
            log.info(
                "Filtering indexed PAF with %d processes in %d byte-local "
                "query chunks",
                worker_count, len(chunks),
            )
            results: Dict[
                int, Tuple[List[AlignmentHit], List[Tuple[int, AlignmentHit]]]
            ] = {}
            completed_queries = 0
            next_progress = progress_step
            context = mp_context()
            if context.get_start_method() == "fork":
                # Share the large parsed .fai dictionary by copy-on-write.
                global _PAF_WORKER_SEQUENCES
                global _PAF_WORKER_AVAILABLE
                _PAF_WORKER_SEQUENCES = sequences
                _PAF_WORKER_AVAILABLE = available
            with context.Pool(
                processes=worker_count,
                initializer=_initialize_paf_process_worker,
                initargs=(
                    paf_path,
                    query_fasta,
                    min_identity,
                    parser_minimum,
                    min_segment,
                    masked_weight,
                    thresholds,
                    source,
                    retain_aligned_pairs,
                    aggregate_query_coverage,
                    interleaved,
                ),
            ) as pool:
                tasks = enumerate(chunks)
                for (
                    chunk_index,
                    query_count,
                    rows,
                    offsets,
                ) in pool.imap_unordered(
                    _filter_indexed_paf_process_chunk,
                    tasks,
                    chunksize=1,
                ):
                    results[chunk_index] = (rows, offsets)
                    completed_queries += query_count
                    if (
                        completed_queries >= next_progress
                        or completed_queries == len(query_index)
                    ):
                        retained_count = sum(
                            len(rows) + len(offsets)
                            for rows, offsets in results.values()
                        )
                        log.info(
                            "PAF parallel filter progress: %d/%d queries; "
                            "%d hits retained",
                            completed_queries,
                            len(query_index),
                            retained_count,
                        )
                        next_progress = completed_queries + progress_step
            retained = []
            retained_with_offsets = []
            for chunk_index in range(len(chunks)):
                rows, offsets = results[chunk_index]
                retained.extend(rows)
                retained_with_offsets.extend(offsets)
        else:
            retained, retained_with_offsets = load_chunk(query_index)
            log.info(
                "PAF streaming/filter progress: %d/%d queries; %d hits retained",
                len(query_index), len(query_index),
                len(retained) + len(retained_with_offsets),
            )
    finally:
        sequences.close()

    if interleaved:
        retained_with_offsets.sort(key=lambda row: row[0])
        retained = [hit for _offset, hit in retained_with_offsets]
    return retained


def aggregate_hits_by_query_coverage(
    hits: Sequence[AlignmentHit], source: str,
) -> List[AlignmentHit]:
    """Collapse accepted hits to unioned query coverage."""
    by_query: Dict[str, List[Tuple[int, int]]] = {}
    for hit in hits:
        by_query.setdefault(hit.query_id, []).extend(hit.query_intervals)
    return [
        AlignmentHit(
            query_id=qid,
            target_id="*",
            identity=100.0,
            aligned_bases=sum(end - start for start, end in intervals),
            query_intervals=intervals,
            source=source,
        )
        for qid, raw_intervals in by_query.items()
        if (intervals := merge_intervals(raw_intervals))
    ]


def write_residual_blast_queries(query_fasta: str, broad_hits: Sequence[AlignmentHit],
                                 out_fasta: str, min_segment: int,
                                 log: logging.Logger,
                                 masked_weight: float = 0.0) -> Dict[str, Tuple[str, int]]:
    """Write per-query gaps not covered by the broad pass, as BLAST queries.

    A gap is emitted when its novelty score reaches ``min_segment``. Returns a
    map from residual id to (original_query_id, offset).
    """
    mapping: Dict[str, Tuple[str, int]] = {}
    Path(out_fasta).parent.mkdir(parents=True, exist_ok=True)
    with IndexedFasta(query_fasta) as fa, open(out_fasta, "wt") as out:
        qlens = {name: fa.length(name) for name in fa.names()}
        cov = coverage_by_query_id(broad_hits, qlens)
        for qid in fa.names():
            qlen = qlens[qid]
            gaps = complement_intervals(cov.get(qid, []), 0, qlen)
            if not gaps:
                continue
            seq = fa.sequence(qid)
            subidx = 0
            for a, b in gaps:
                sub = seq[a:b]
                for na, nb, nfree in iter_n_free_segments(sub):
                    score = (
                        novelty_score(nfree, masked_weight)
                        if masked_weight
                        else count_unmasked(nfree)
                    )
                    if score < min_segment:
                        continue
                    subidx += 1
                    s0 = a + na
                    rid = f"{qid}_resgap{subidx:04d}_{s0}_{a + nb}"
                    mapping[rid] = (qid, s0)
                    out.write(f">{rid} {qid}:{s0}-{a + nb}+\n{wrap_fasta(nfree)}\n")
    log.info("Wrote %d residual BLAST query gaps", len(mapping))
    return mapping


def remap_residual_blast_hits(hits: Sequence[AlignmentHit],
                              residual_map: Dict[str, Tuple[str, int]]) -> List[AlignmentHit]:
    out: List[AlignmentHit] = []
    for hit in hits:
        if hit.query_id not in residual_map:
            continue
        parent, offset = residual_map[hit.query_id]
        q_intervals = [(offset + a, offset + b) for a, b in hit.query_intervals if b > a]
        if not q_intervals:
            continue
        pairs = [(offset + q0, offset + q1, t0, t1) for q0, q1, t0, t1 in hit.aligned_pairs]
        out.append(dataclasses.replace(hit, query_id=parent, query_intervals=q_intervals,
                                       aligned_pairs=pairs, source="blastn_residual"))
    return out


def collect_alignment_hits(query_fasta: str, target_fasta: str, workdir: str, threads: int,
                           min_identity: float, min_segment: int, log: logging.Logger,
                           winnow_params: Optional[dict] = None, rep_kmers: Optional[str] = None,
                           blast_word_size: int = 50, blast_evalue: str = "1e-300",
                           blast_max_target_seqs: int = 100, prefix: str = "align",
                           skip_blastn: bool = False, reuse_existing: bool = False,
                           blast_mode: str = "residual", broad_aligner: str = "winnowmap",
                           minimap_params: Optional[dict] = None,
                           blast_db_prefix: Optional[str] = None,
                           masked_weight: float = 0.0,
                           min_score_by_query: Optional[Dict[str, float]] = None,
                           parse_min_segment: Optional[int] = None,
                           query_workers: int = 1,
                           retain_aligned_pairs: bool = True,
                           aggregate_query_coverage: bool = False) -> List[AlignmentHit]:
    """Collect valid Winnowmap and BLASTN hits.

    ``blast_mode='residual'`` BLASTs only broad-aligner gaps (the efficient mode
    for large-main alignment). ``blast_mode='independent'`` BLASTs the full
    query independently (used for insertion pools and insertion self-alignment).
    Every returned hit reaches ``min_segment`` query novelty score. With the
    default ``masked_weight=0``, this is the original unmasked-base rule.
    """
    if blast_mode not in {"residual", "independent"}:
        raise ValueError("blast_mode must be 'residual' or 'independent'")
    if broad_aligner not in {"minimap2", "winnowmap"}:
        raise ValueError("broad_aligner must be 'minimap2' or 'winnowmap'")
    parser_minimum = (
        min_segment if parse_min_segment is None else int(parse_min_segment)
    )
    if parser_minimum < 1:
        raise ValueError("parse_min_segment must be positive")
    Path(workdir).mkdir(parents=True, exist_ok=True)
    paf = os.path.join(workdir, f"{prefix}.{broad_aligner}.paf")
    if broad_aligner == "minimap2":
        run_minimap2(query_fasta, target_fasta, paf, threads, log,
                     params=minimap_params, reuse_existing=reuse_existing)
    else:
        run_winnowmap(query_fasta, target_fasta, paf, threads, log,
                      params=winnow_params, rep_kmers=rep_kmers,
                      reuse_existing=reuse_existing)
    paf_size = os.path.getsize(paf) if os.path.isfile(paf) else 0
    log.info(
        "Parsing existing %s PAF once (%.2f GiB)",
        broad_aligner, paf_size / (1024 ** 3),
    )
    if paf_size > PAF_INDEX_THRESHOLD_BYTES:
        hits = load_filtered_paf_hits_indexed(
            paf,
            query_fasta,
            min_identity,
            parser_minimum,
            min_segment,
            masked_weight,
            min_score_by_query,
            log,
            source=broad_aligner,
            retain_aligned_pairs=retain_aligned_pairs,
            aggregate_query_coverage=aggregate_query_coverage,
            workers=query_workers,
        )
        if query_workers > 1:
            log.info(
                "Parsed %d valid %s hits with indexed multiprocessing "
                "(%d requested workers)",
                len(hits), broad_aligner, query_workers,
            )
        else:
            log.info(
                "Parsed %d valid %s hits with indexed single-process "
                "query streaming",
                len(hits), broad_aligner,
            )
    else:
        hits = filter_hits_by_unmasked_query(
            list(parse_winnowmap_paf(
                paf,
                min_identity,
                parser_minimum,
                source=broad_aligner,
                retain_aligned_pairs=retain_aligned_pairs,
            )),
            query_fasta,
            min_segment,
            masked_weight,
            min_score_by_query,
            query_workers,
        )
        if aggregate_query_coverage:
            hits = aggregate_hits_by_query_coverage(
                hits, f"{broad_aligner}_query_coverage",
            )
        log.info(
            "Parsed %d valid %s hits in memory (PAF <= 100 MB)",
            len(hits), broad_aligner,
        )
    if skip_blastn:
        return hits

    residual_map: Optional[Dict[str, Tuple[str, int]]]
    if blast_mode == "residual":
        blast_query = os.path.join(workdir, f"{prefix}.blast_residual.fa")
        residual_map = write_residual_blast_queries(
            query_fasta, hits, blast_query, min_segment, log, masked_weight,
        )
        if not residual_map:
            log.info(
                "BLASTN rescue skipped for %s: no residual gaps reach score %d",
                prefix, min_segment,
            )
            return hits
    else:
        blast_query = query_fasta
        residual_map = None

    db_prefix = blast_db_prefix or os.path.join(workdir, f"{prefix}.blastdb", "ref")
    db_sequence_names = db_prefix + ".seq"
    database_current = (
        os.path.isfile(db_sequence_names)
        and os.path.getsize(db_sequence_names) > 0
        and os.stat(db_sequence_names).st_mtime_ns
            >= os.stat(target_fasta).st_mtime_ns
    )
    if database_current:
        log.info("REUSE: fixed-reference BLAST database %s", db_prefix)
    else:
        build_blast_db(target_fasta, db_prefix, log)
    sam = os.path.join(workdir, f"{prefix}.blastn.sam")
    run_blastn17(blast_query, db_prefix, sam, threads, min_identity,
                 blast_word_size, blast_evalue, blast_max_target_seqs, log)
    bl_raw = filter_hits_by_unmasked_query(
        list(parse_blast17_sam(sam, db_prefix + ".seq", min_identity, parser_minimum,
                               query_lengths=fasta_lengths(blast_query))),
        blast_query, min_segment, masked_weight, min_score_by_query,
        query_workers,
    )
    bl_hits = remap_residual_blast_hits(bl_raw, residual_map) if residual_map is not None else bl_raw
    log.info("Parsed %d %s BLASTN hits; retained %d", len(bl_raw), blast_mode, len(bl_hits))
    hits.extend(bl_hits)
    if aggregate_query_coverage:
        hits = aggregate_hits_by_query_coverage(
            hits, "combined_query_coverage",
        )
    return hits


def collect_blastn_hits(
    query_fasta: str,
    target_fasta: str,
    workdir: str,
    threads: int,
    min_identity: float,
    min_segment: int,
    log: logging.Logger,
    *,
    blast_word_size: int = 50,
    blast_evalue: str = "1e-300",
    blast_max_target_seqs: int = 100,
    prefix: str = "align",
    blast_db_prefix: Optional[str] = None,
    masked_weight: float = 0.0,
    min_score_by_query: Optional[Dict[str, float]] = None,
    parse_min_segment: Optional[int] = None,
    reuse_existing: bool = False,
    query_workers: int = 1,
) -> List[AlignmentHit]:
    """Run only BLASTN and return query-score-filtered alignment hits.

    This is used by staged cleaners that must inspect the residual after a
    broad aligner before deciding whether BLASTN is necessary.  Unlike
    ``collect_alignment_hits``, it does not rerun Minimap2 or Winnowmap.
    """
    parser_minimum = (
        min_segment if parse_min_segment is None else int(parse_min_segment)
    )
    if parser_minimum < 1:
        raise ValueError("parse_min_segment must be positive")
    Path(workdir).mkdir(parents=True, exist_ok=True)
    db_prefix = blast_db_prefix or os.path.join(
        workdir, f"{prefix}.blastdb", "ref",
    )
    db_sequence_names = db_prefix + ".seq"
    database_current = (
        os.path.isfile(db_sequence_names)
        and os.path.getsize(db_sequence_names) > 0
        and os.stat(db_sequence_names).st_mtime_ns
            >= os.stat(target_fasta).st_mtime_ns
    )
    if database_current:
        log.info("REUSE: fixed-reference BLAST database %s", db_prefix)
    else:
        build_blast_db(target_fasta, db_prefix, log)
    sam = os.path.join(workdir, f"{prefix}.blastn.sam")
    run_blastn17(
        query_fasta,
        db_prefix,
        sam,
        threads,
        min_identity,
        blast_word_size,
        blast_evalue,
        blast_max_target_seqs,
        log,
        reuse_existing=reuse_existing,
    )
    raw_hits = list(parse_blast17_sam(
        sam,
        db_prefix + ".seq",
        min_identity,
        parser_minimum,
        query_lengths=fasta_lengths(query_fasta),
    ))
    hits = filter_hits_by_unmasked_query(
        raw_hits,
        query_fasta,
        min_segment,
        masked_weight,
        min_score_by_query,
        query_workers,
    )
    log.info("Parsed %d valid staged BLASTN hits", len(hits))
    return hits
