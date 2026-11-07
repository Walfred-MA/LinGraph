#!/usr/bin/env python3
"""Align the gaps the normal gap fill left, with k-mer blocks.

Input gaps are tools/find_pseudolinear_gaps.py's PREFIX.query.bed: UNMAPPED
query sequence between two collinear anchors (same reference contig and
strand), with the reference gap between the anchors. Run the finder with
--min-size 0 to take every remaining gap.

For every gap:
  1. regions: the query gap and the reference gap, each with --flank bp (30,000)
     of anchoring sequence on both sides (clipped at contig ends); a '-' gap's
     query region is reverse-complemented into reference orientation;
  2. k-mer blocks (tools/kmer_blocks.py): the reference region cut with its own
     31-mers occurring < 255 times, the query region with that same k-mer set;
  3. BlockMatchDP: --pair-score log (similarity capped at 0.999), --gap-factor
     0.1, --max-occurrence 255, --top-cells 15000;
  4. tools/block_guided_align.py (±15 kb block alignments, trimming, merge by
     similarity; query gaps relocated to a clear best reference hit inside the
     region, else inserted);
  5. events in genome coordinates: DEL overlapping the reference gap, INS/REL
     whose query interval overlaps the query gap (flank events are the anchors'
     business and are dropped).

Work files go to PREFIX.work/GAP/ (a finished gap is skipped on rerun; --only
runs a subset, e.g. one SLURM array task). The combined outputs cover every
finished gap of the BED:
  PREFIX.events.tsv   gap, type, reference interval, size, query interval, source
  PREFIX.summary.tsv  per gap: sizes, blocks, DP pairs, aligned / relocated /
                      inserted query bp, events

    python tools/fill_remaining_gaps.py --gaps SAMPLE.query.bed -q SAMPLE.fa \
        -r CHM13.fa -o SAMPLE.remaining --blockmatchdp minsetref/BlockMatchDP -T 16
"""
import argparse
import csv
import gzip
import os
import subprocess
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
import block_guided_align  # noqa: E402
import kmer_blocks  # noqa: E402

COMPLEMENT = str.maketrans("ACGTacgtNn", "TGCAtgcaNn")


def read_gaps(path):
    gaps = []
    with open(path) as handle:
        header = handle.readline().lstrip("#").rstrip("\n").split("\t")
        for line in handle:
            row = dict(zip(header, line.rstrip("\n").split("\t")))
            gaps.append(dict(name=row["name"], query_contig=row["query_contig"],
                             query_start=int(row["query_start"]), query_end=int(row["query_end"]),
                             strand=row["strand"], reference_contig=row["reference_contig"],
                             reference_start=int(row["reference_start"]),
                             reference_end=int(row["reference_end"])))
    return gaps


def load_contigs(path, wanted):
    """The wanted records of a FASTA (plain or gzip), uppercase/lowercase kept."""
    opener = gzip.open if path.endswith(".gz") else open
    contigs, name, parts = {}, None, []
    with opener(path, "rt") as handle:
        for line in handle:
            if line.startswith(">"):
                if name in wanted:
                    contigs[name] = "".join(parts)
                name, parts = line[1:].split()[0], []
            elif name in wanted:
                parts.append(line.strip())
    if name in wanted:
        contigs[name] = "".join(parts)
    missing = wanted - set(contigs)
    if missing:
        raise SystemExit(f"{path}: missing {', '.join(sorted(missing)[:5])}")
    return contigs


def run_gap(gap, query, reference, args, work):
    q_seq, r_seq = query[gap["query_contig"]], reference[gap["reference_contig"]]
    qa, qb = max(0, gap["query_start"] - args.flank), min(len(q_seq), gap["query_end"] + args.flank)
    ra = max(0, gap["reference_start"] - args.flank)
    rb = min(len(r_seq), gap["reference_end"] + args.flank)
    region_q = q_seq[qa:qb]
    if gap["strand"] == "-":
        region_q = region_q.translate(COMPLEMENT)[::-1]
    region_r = r_seq[ra:rb]
    work.mkdir(parents=True, exist_ok=True)
    for name, seq in (("query.fa", region_q), ("reference.fa", region_r)):
        with open(work / name, "w") as out:
            out.write(f">{name[:-3]}\n{seq}\n")

    valid = kmer_blocks.valid_kmers(region_r, args.max_occurrence)
    r_blocks = kmer_blocks.segment(region_r, None, args.max_occurrence)
    q_blocks = kmer_blocks.segment(region_q, valid, args.max_occurrence)
    kmer_blocks.write_cells(str(work / "r.kcells.fa"), "r", region_r, r_blocks)
    kmer_blocks.write_cells(str(work / "q.kcells.fa"), "q", region_q, q_blocks)
    result = subprocess.run([args.blockmatchdp, "-q", str(work / "q.kcells.fa"),
                             "-t", str(work / "r.kcells.fa"), "-o", str(work / "dp"),
                             "--max-occurrence", str(args.max_occurrence),
                             "--top-cells", str(args.top_cells), "--gap-factor", "0.1",
                             "--pair-score", "log", "-T", str(args.threads)]
                            + (["--occurrence-source", args.occurrence_source]
                               if args.occurrence_source != "both" else [])
                            + (["--gap-open-mean", str(args.gap_open_mean)] if args.gap_open_mean else [])
                            + (["--gap-open-factor", str(args.gap_open_factor)]
                               if args.gap_open_factor is not None else []),
                            capture_output=True, text=True)
    if result.returncode:
        sys.stderr.write(result.stderr)
        raise SystemExit(f"{gap['name']}: BlockMatchDP failed")
    block_guided_align.main(["-q", str(work / "query.fa"), "-t", str(work / "reference.fa"),
                             "--dp", str(work / "dp.dp.tsv"), "-o", str(work / "bga"),
                             "-T", str(args.threads), "--minimap2", args.minimap2,
                             # the region ends are the flanks beside the gap's anchors
                             "--anchored-ends"]
                            + (["--replace-uncertain"] if args.replace_uncertain else []))

    # Events in genome coordinates, limited to the gap itself.
    def query_interval(s, e):
        return (qa + s, qa + e) if gap["strand"] == "+" else (qb - e, qb - s)

    events = []
    for row in csv.DictReader(open(work / "bga.events.tsv"), delimiter="\t"):
        rs, re_ = ra + int(row["ref_start"]), ra + int(row["ref_end"])
        qs, qe = query_interval(int(row["query_start"]), int(row["query_end"]))
        if row["type"] == "DEL":
            keep = rs < gap["reference_end"] and re_ > gap["reference_start"]
        else:
            keep = qs < gap["query_end"] and qe > gap["query_start"]
        if keep:
            events.append((gap["name"], row["type"], gap["reference_contig"], rs, re_, row["size"],
                           gap["query_contig"], qs, qe, row["source"]))
    with open(work / "events.genome.tsv", "w") as out:
        for e in events:
            out.write("\t".join(map(str, e)) + "\n")
    bga = dict(line.rstrip("\n").split("\t") for line in open(work / "bga.summary"))
    summary = [gap["name"], gap["query_end"] - gap["query_start"],
               gap["reference_end"] - gap["reference_start"], len(q_blocks), len(r_blocks),
               bga["paired_blocks"], bga["query_bp_block"], bga["query_bp_relocated"],
               bga["query_bp_ins"], sum(e[1] == "INS" for e in events),
               sum(e[1] == "DEL" for e in events), sum(e[1] == "REL" for e in events)]
    with open(work / "gap.summary.tsv", "w") as out:
        out.write("\t".join(map(str, summary)) + "\n")
    (work / "done").write_text("ok\n")


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--gaps", required=True, help="find_pseudolinear_gaps.py PREFIX.query.bed")
    parser.add_argument("-q", "--query", required=True, help="the sample's assembly FASTA")
    parser.add_argument("-r", "--reference", required=True, help="reference FASTA")
    parser.add_argument("-o", "--output", required=True, help="output prefix")
    parser.add_argument("-T", "--threads", type=int, default=8)
    parser.add_argument("--flank", type=int, default=30_000)
    parser.add_argument("--max-occurrence", type=int, default=255)
    parser.add_argument("--top-cells", type=int, default=15_000)
    parser.add_argument("--blockmatchdp", default="BlockMatchDP")
    parser.add_argument("--occurrence-source", choices=("both", "target"), default="both",
                        help="BlockMatchDP --occurrence-source: count k-mer occurrences for "
                             "--max-occurrence in both sequences (default) or the reference only")
    parser.add_argument("--replace-uncertain", dest="replace_uncertain", action="store_true", default=False,
                        help="block_guided_align.py --replace-uncertain: re-place expansions whose copies "
                             "cannot be told apart (off by default)")
    parser.add_argument("--no-replace-uncertain", dest="replace_uncertain", action="store_false",
                        help="keep BlockMatchDP's pairs (default)")
    parser.add_argument("--gap-open-factor", type=float,
                        help="BlockMatchDP --gap-open-factor: cost of the first unpaired block of a run, "
                             "times its self-match (default: not passed, i.e. equal to the 0.1 gap factor)")
    parser.add_argument("--gap-open-mean", type=float, default=0.0,
                        help="BlockMatchDP --gap-open-mean: extra cost per run of unpaired blocks, in "
                             "mean block self-matches (default: 0, off)")
    parser.add_argument("--minimap2", default="minimap2")
    parser.add_argument("--only", help="comma-separated gap names to run")
    args = parser.parse_args(argv)

    gaps = read_gaps(args.gaps)
    todo_names = set(args.only.split(",")) if args.only else None
    work_root = Path(args.output + ".work")
    todo = [g for g in gaps if (todo_names is None or g["name"] in todo_names)
            and not (work_root / g["name"] / "done").is_file()]
    if todo:
        query = load_contigs(args.query, {g["query_contig"] for g in todo})
        reference = load_contigs(args.reference, {g["reference_contig"] for g in todo})
        for number, gap in enumerate(todo, 1):
            print(f"[remaining] {number}/{len(todo)} {gap['name']}: query "
                  f"{gap['query_end'] - gap['query_start']} bp, reference "
                  f"{gap['reference_end'] - gap['reference_start']} bp", file=sys.stderr, flush=True)
            run_gap(gap, query, reference, args, work_root / gap["name"])

    finished = [g for g in gaps if (work_root / g["name"] / "done").is_file()]
    with open(args.output + ".events.tsv", "w") as out:
        out.write("gap\ttype\tref_contig\tref_start\tref_end\tsize\tquery_contig\tquery_start\t"
                  "query_end\tsource\n")
        for g in finished:
            out.write(open(work_root / g["name"] / "events.genome.tsv").read())
    with open(args.output + ".summary.tsv", "w") as out:
        out.write("gap\tquery_gap_bp\treference_gap_bp\tquery_blocks\treference_blocks\tdp_pairs\t"
                  "query_bp_block\tquery_bp_relocated\tquery_bp_inserted\tins\tdel\trel\n")
        for g in finished:
            out.write(open(work_root / g["name"] / "gap.summary.tsv").read())
    print(f"[remaining] {len(finished)}/{len(gaps)} gaps finished -> {args.output}.events.tsv",
          file=sys.stderr)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
