#!/usr/bin/env python3
"""Align one query sequence to one reference sequence through a large tandem
repeat array (e.g. a centromeric HOR array) and write one graph CIGAR.

Steps (the k-mer-block method of the remaining-gap fill):
  1. k-mer blocks (kmer_blocks.py): the reference is cut into repeat-unit copies
     with its own canonical 31-mers seen < --max-occurrence times, the query with
     the reference's usable 31-mers. The two sequences are cut at the same time.
  2. BlockMatchDP (threaded k-mer matching): --pair-score log, --gap-factor 0.1;
     the two ordered copy lists are aligned end to end.
  3. block_guided_align.py with --anchored-ends: every paired copy is aligned at
     base level with +-15 kb flanks (in parallel), pieces are trimmed and merged
     by similarity, and every query gap is mapped to the whole reference with
     minimap2 (threaded); a clear best hit is a relocation, the rest insertion.
     --no-relocation inserts every gap instead. --replace-uncertain re-places
     expansions whose copies cannot be told apart in a fixed way (off by default).

Output (-o): one line, a graph CIGAR covering the whole query and the whole
reference. The main path is >REF:OPS; a relocated piece is its own segment
>REF:STARTH...ENDH, or <REF:...H...H on the reverse strand (H counts the
reference outside it, in its own direction), placed at its query position; a
following >OPS resumes the main path. Ops are =, X (payload: reference bases,
then query bases), I (payload: the inserted query bases) and D. Reference the
main path skips is D, query bases nothing aligned are I.

Intermediate files are written to OUTPUT.horalign/ and removed after success
unless --keep-all is given; after a failure they are kept.

    python LinGraph.py horalign -r reference.fa -q query.fa -o result.gcigar -t 16
"""
import argparse
import gzip
import multiprocessing as mp
import os
import shutil
import subprocess
import sys
import time
from pathlib import Path

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
import block_guided_align  # noqa: E402
import kmer_blocks  # noqa: E402

MARKER = ".horalign"


def single_record(path):
    """(name, sequence) of a FASTA (plain or gzip) that must hold one record."""
    opener = gzip.open if str(path).endswith(".gz") else open
    found, parts = [], []
    with opener(path, "rt") as handle:
        for line in handle:
            if line.startswith(">"):
                found.append(line[1:].split()[0] if line[1:].split() else "")
                if len(found) > 1:
                    raise SystemExit(f"{path}: horalign takes one sequence per FASTA, found more")
            elif found:
                parts.append(line.strip())
    if not found or not parts:
        raise SystemExit(f"{path}: no sequence")
    return found[0], "".join(parts)


def cut_cells(task):
    """k-mer blocks of one sequence; the query uses the reference's usable k-mers."""
    label, seq, reference, max_occurrence, path = task
    valid = None
    if reference is not None:
        valid = kmer_blocks.valid_kmer_array(reference, max_occurrence)
    blocks = kmer_blocks.segment(seq, valid, max_occurrence)
    kmer_blocks.write_cells(path, label, seq, blocks)
    return label, len(blocks)


def find_blockmatchdp(given):
    for candidate in (given, HERE / "BlockMatchDP", shutil.which("BlockMatchDP")):
        if candidate and Path(candidate).is_file() and os.access(candidate, os.X_OK):
            return str(candidate)
    raise SystemExit("BlockMatchDP not found: build it with `python scripts/install.py` "
                     "or pass --blockmatchdp")


def log(message):
    print(f"[horalign] {message}", file=sys.stderr, flush=True)


def main(argv=None):
    parser = argparse.ArgumentParser(prog="LinGraph.py horalign", description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("-r", "--reference", required=True, help="reference FASTA, one sequence")
    parser.add_argument("-q", "--query", required=True, help="query FASTA, one sequence (same strand)")
    parser.add_argument("-o", "--output", required=True, help="output graph CIGAR file")
    parser.add_argument("-t", "--threads", type=int, default=8)
    parser.add_argument("--keep-all", action="store_true",
                        help="keep the intermediate files in OUTPUT.horalign/")
    parser.add_argument("--replace-uncertain", action="store_true",
                        help="re-place expansions whose copies cannot be told apart: pair them in order "
                             "from the left (off by default)")
    parser.add_argument("--no-relocation", action="store_true",
                        help="do not map query gaps to the reference; every gap is inserted sequence")
    advanced = parser.add_argument_group("advanced")
    advanced.add_argument("--max-occurrence", type=int, default=255,
                          help="use 31-mers seen fewer times than this (default: 255)")
    advanced.add_argument("--top-cells", type=int, default=15_000,
                          help="BlockMatchDP --top-cells (default: 15000)")
    advanced.add_argument("--flank", type=int, default=15_000,
                          help="bp added on both sides of each paired copy (default: 15000)")
    advanced.add_argument("--clear-ratio", type=float, default=0.9,
                          help="a relocation hit is clear when every overlapping rival scores "
                               "below this fraction of it (default: 0.9)")
    advanced.add_argument("--blockmatchdp", help="BlockMatchDP executable (default: beside this script)")
    advanced.add_argument("--minimap2", default="minimap2")
    args = parser.parse_args(argv)
    if args.threads < 1:
        parser.error("-t must be at least 1")

    started = time.time()
    blockmatchdp = find_blockmatchdp(args.blockmatchdp)
    ref_name, ref_seq = single_record(args.reference)
    query_name, query_seq = single_record(args.query)
    output = Path(args.output)
    work = Path(str(output) + ".horalign")
    if work.exists() and not (work / MARKER).exists():
        raise SystemExit(f"{work} exists and was not made by horalign; remove it or choose another -o")
    work.mkdir(parents=True, exist_ok=True)
    (work / MARKER).write_text("horalign work folder\n")
    for name, (record, seq) in (("reference.fa", (ref_name, ref_seq)), ("query.fa", (query_name, query_seq))):
        with open(work / name, "w") as out:
            out.write(f">{record}\n{seq}\n")
    log(f"reference {ref_name} {len(ref_seq)} bp; query {query_name} {len(query_seq)} bp; "
        f"{args.threads} threads")

    # 1. k-mer blocks, both sequences at once.
    tasks = [("r", ref_seq, None, args.max_occurrence, str(work / "r.kcells.fa")),
             ("q", query_seq, ref_seq, args.max_occurrence, str(work / "q.kcells.fa"))]
    if args.threads > 1:
        with mp.get_context("fork").Pool(2) as pool:
            cells = dict(pool.map(cut_cells, tasks))
    else:
        cells = dict(map(cut_cells, tasks))
    log(f"k-mer blocks: reference {cells['r']}, query {cells['q']}")

    # 2. BlockMatchDP.
    with open(work / "BlockMatchDP.log", "w") as handle:
        result = subprocess.run([blockmatchdp, "-q", str(work / "q.kcells.fa"),
                                 "-t", str(work / "r.kcells.fa"), "-o", str(work / "dp"),
                                 "--max-occurrence", str(args.max_occurrence),
                                 "--top-cells", str(args.top_cells), "--gap-factor", "0.1",
                                 "--pair-score", "log", "-T", str(args.threads)],
                                stdout=handle, stderr=subprocess.STDOUT)
    if result.returncode:
        raise SystemExit(f"BlockMatchDP failed; see {work / 'BlockMatchDP.log'}")
    log("BlockMatchDP done")

    # 3. Block-guided base alignment and the graph CIGAR.
    tmp = Path(str(output) + ".tmp")
    block_guided_align.main(
        ["-q", str(work / "query.fa"), "-t", str(work / "reference.fa"),
         "--dp", str(work / "dp.dp.tsv"), "-o", str(work / "bga"), "-T", str(args.threads),
         "--flank", str(args.flank), "--clear-ratio", str(args.clear_ratio),
         "--minimap2", args.minimap2, "--anchored-ends",
         "--graph-cigar", str(tmp), "--reference-name", ref_name]
        + (["--no-relocation"] if args.no_relocation else [])
        + (["--replace-uncertain"] if args.replace_uncertain else []))
    os.replace(tmp, output)

    summary = dict(line.rstrip("\n").split("\t", 1) for line in open(work / "bga.summary"))
    log(f"query bp aligned {int(summary['query_bp_block']) + int(summary['query_bp_gapfill'])}, "
        f"relocated {summary['query_bp_relocated']}, inserted {summary['query_bp_ins']}; "
        f"relocation segments {summary['graph_cigar_relocations']}; "
        f"{time.time() - started:.0f} s")
    if args.keep_all:
        log(f"wrote {output}; intermediate files in {work}")
    else:
        shutil.rmtree(work)
        log(f"wrote {output}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
