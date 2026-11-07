#!/usr/bin/env python3
"""Find large unfilled gaps in pseudo-linear assignments that two anchors bound.

Along each query contig, consecutive PRIMARY rows of the location layer
(``pa_path`` = ``pri``) are anchors. An anchor's direction combines the strands
of its query and reference intervals ('+' when they agree, '-' otherwise), so
a contig placed in reverse runs backwards along the reference. When UNMAPPED query sequence lies between
two consecutive anchors and both anchors map to the same reference contig, on
the same strand, in collinear order, the gap has the same two ends on the
query and on the reference and can still be filled. It is reported as

    query gap      = left anchor query end .. right anchor query start
    reference gap  = the reference span between the two anchors

For a ``-`` strand pair the reference gap is right anchor end .. left anchor
start. When the anchors overlap on the reference, the reference gap is the
empty interval at the left anchor's reference end and ``reference_overlap``
gives the overlap. A gap is kept when the query or the reference gap is at
least ``--min-size`` bp (default 10,000).

Gaps are numbered first; optional filters then drop gaps without changing
the names of the others:
  --vote-window BP       safety check (one input; --genomelift needed): the
                         placed blocks within BP upstream and within BP
                         downstream each vote, weighted by length; each side's
                         winner (reference contig + direction) must be the
                         gap's, hold > 50% of its side and lie mostly on that
                         side of the reference gap; the contig must extend at
                         least --min-contig-flank bp (1 Mb) beyond the gap on
                         both sides (--fasta-query needed); the GenomeLift
                         placements
                         of the blocks in the same two windows must each put
                         > 50% of their bp on the gap's reference contig, on
                         that side of the reference gap; and GenomeLift
                         placements inside the query gap, when there are any,
                         must (> 50% of their bp) fall within --liftover-slack
                         of the reference gap
  --min-anchor-score S   both anchors' alignment_score (the complete-row score
                         fill_graphcigartoref_gaps.py gates at 1000) must be >= S
  --fasta-query/--fasta-reference   (one input only; .fai required) gaps whose
                         query or reference gap contains an N are dropped

Outputs, for ``-o PREFIX`` (all coordinates zero-based, half-open):

    PREFIX.query.bed      query gap BED6 + reference gap and details
    PREFIX.reference.bed  reference gap BED6 + query gap and details
    PREFIX.dropped.bed    the gaps the filters dropped, with the reason: anchor_score,
                          N (an N in the query or reference gap) or safety
    PREFIX.summary.tsv    counts per sample

Both BEDs carry the same gap names (``SAMPLE_gapN``) so a gap can be paired
across them. Rows are sorted by coordinate.
"""
from __future__ import annotations

import argparse
from collections import defaultdict
import os
from pathlib import Path
import re
import sys

COORD_RE = re.compile(r"^(.+):(\d+)-(\d+)([+-])$")
DETAIL_FIELDS = (
    "query_gap_bp", "reference_gap_bp", "reference_overlap", "unmapped_bp",
    "unmapped_rows", "insertion_bp", "left_anchor_bp", "right_anchor_bp",
    "left_anchor", "right_anchor", "left_anchor_score", "right_anchor_score",
    "query_n_bp", "reference_n_bp", "upstream_vote", "downstream_vote", "liftover_vote",
)


class Fai:
    """Random access to an uncompressed, samtools-indexed FASTA."""

    def __init__(self, path):
        if not os.path.isfile(path + ".fai"):
            raise ValueError(f"{path}: needs a samtools faidx index ({path}.fai)")
        self.path, self.index = path, {}
        for line in open(path + ".fai"):
            name, length, offset, bases, width = line.split("\t")[:5]
            self.index[name] = (int(length), int(offset), int(bases), int(width))

    def count_n(self, name, start, end):
        if name not in self.index:
            raise ValueError(f"{self.path}: no record {name}")
        length, offset, bases, width = self.index[name]
        start, end = max(0, start), min(length, end)
        if end <= start:
            return 0
        first = offset + start // bases * width + start % bases
        last = offset + (end - 1) // bases * width + (end - 1) % bases
        with open(self.path, "rb") as handle:
            handle.seek(first)
            chunk = handle.read(last - first + 1)
        return chunk.count(b"N") + chunk.count(b"n")


def parse_coord(text):
    match = COORD_RE.match(text)
    if match is None:
        return None
    contig, start, end, strand = match.groups()
    return contig, int(start), int(end), strand


def sample_name(path):
    name = os.path.basename(path)
    for suffix in (".pseudolinear.tsv", ".tsv"):
        if name.endswith(suffix):
            return name[:-len(suffix)]
    return name


def read_rows(path):
    """Return per-query-contig anchors and gap evidence, streaming the file."""
    anchors = defaultdict(list)    # contig -> (qs, qe, ref, rs, re, strand, score)
    unmapped = defaultdict(list)   # contig -> (qs, qe)
    inserted = defaultdict(list)   # contig -> (qs, qe)
    with open(path) as handle:
        header = handle.readline().rstrip("\n").split("\t")
        try:
            kind_at = header.index("type")
            query_at = header.index("query_interval")
            reference_at = header.index("reference_interval")
        except ValueError:
            raise ValueError(f"{path}: not a readable pseudolinear TSV "
                             "(needs type, query_interval, reference_interval)") from None
        path_at = header.index("pa_path") if "pa_path" in header else None
        score_at = header.index("alignment_score") if "alignment_score" in header else None
        width = max(kind_at, query_at, reference_at, path_at or 0, score_at or 0) + 1
        for line_number, line in enumerate(handle, 2):
            fields = line.rstrip("\n").split("\t", width)
            if len(fields) < width:
                continue
            if path_at is not None and fields[path_at] not in ("pri", "", "."):
                continue
            kind = fields[kind_at]
            query = parse_coord(fields[query_at])
            if query is None or query[2] <= query[1]:
                continue
            contig, qs, qe, query_strand = query
            if kind == "PRIMARY":
                reference = parse_coord(fields[reference_at])
                if reference is None:
                    continue
                # The placement direction combines both strands: a '-' query
                # interval on a '+' reference runs backwards along the reference.
                ref_contig, rs, re_, ref_strand = reference
                direction = "+" if query_strand == ref_strand else "-"
                try:
                    score = float(fields[score_at]) if score_at is not None else None
                except ValueError:
                    score = None
                anchors[contig].append((qs, qe, ref_contig, rs, re_, direction, score))
            elif kind == "UNMAPPED":
                unmapped[contig].append((qs, qe))
            elif kind == "INSERTION":
                inserted[contig].append((qs, qe))
    return anchors, unmapped, inserted


def overlap_bp(intervals, start, end):
    return sum(max(0, min(e, end) - max(s, start)) for s, e in intervals)


def find_gaps(sample, anchors, unmapped, inserted, min_size):
    gaps = []
    for contig in sorted(anchors):
        rows = sorted(anchors[contig], key=lambda anchor: anchor[:6])
        holes = sorted(unmapped.get(contig, ()))
        insertions = inserted.get(contig, ())
        cursor = 0
        for left, right in zip(rows, rows[1:]):
            qstart, qend = left[1], right[0]
            # Only UNMAPPED query sequence makes an unfilled gap.
            while cursor < len(holes) and holes[cursor][1] <= qstart:
                cursor += 1
            inside = []
            index = cursor
            while index < len(holes) and holes[index][0] < qend:
                inside.append(holes[index])
                index += 1
            unmapped_bp = overlap_bp(inside, qstart, qend)
            if unmapped_bp == 0:
                continue
            _, _, lref, lrs, lre, lstrand, lscore = left
            _, _, rref, rrs, rre, rstrand, rscore = right
            if lref != rref or lstrand != rstrand:
                continue
            if lstrand == "+":
                if rrs < lrs:
                    continue
                gap_start, gap_end = lre, rrs
            else:
                if rre > lre:
                    continue
                gap_start, gap_end = rre, lrs
            overlap = max(0, gap_start - gap_end)
            if overlap:
                gap_start = gap_end = lre if lstrand == "+" else lrs
            query_bp, reference_bp = qend - qstart, gap_end - gap_start
            if max(query_bp, reference_bp) < min_size:
                continue
            gaps.append(dict(
                sample=sample, query_contig=contig, query_start=qstart,
                query_end=qend, reference_contig=lref, reference_start=gap_start,
                reference_end=gap_end, strand=lstrand,
                query_gap_bp=query_bp, reference_gap_bp=reference_bp,
                reference_overlap=overlap, unmapped_bp=unmapped_bp,
                unmapped_rows=len(inside),
                insertion_bp=overlap_bp(insertions, qstart, qend),
                left_anchor_bp=left[1] - left[0], right_anchor_bp=right[1] - right[0],
                left_anchor=f"{lref}:{lrs}-{lre}{lstrand}",
                right_anchor=f"{rref}:{rrs}-{rre}{rstrand}",
                left_anchor_score="." if lscore is None else f"{lscore:g}",
                right_anchor_score="." if rscore is None else f"{rscore:g}",
                query_n_bp=".", reference_n_bp=".",
                upstream_vote=".", downstream_vote=".", liftover_vote=".",
            ))
    for number, gap in enumerate(sorted(gaps, key=lambda g: (
            g["query_contig"], g["query_start"])), 1):
        gap["name"] = f"{sample}_gap{number}"
    return gaps


def read_liftover(paths):
    """GenomeLift placements: query contig -> [(qs, qe, ref contig, rs, re)].

    The final location mapping (column 17, tier|references|interval|score) is
    used, column 10 when it is empty. Later files (genomeliftfix) override rows
    of earlier ones with the same block name and assembly locus."""
    rows = {}
    for path in paths:
        with open(path) as handle:
            for line in handle:
                fields = line.rstrip("\n").split("\t")
                if len(fields) < 10 or fields[0] in ("allelename", "DEL", "NA") or line.startswith("#"):
                    continue
                interval = ""
                if len(fields) >= 17 and fields[16].strip() not in ("", "."):
                    parts = fields[16].split("|", 3)
                    interval = parts[2].strip() if len(parts) == 4 else ""
                if interval in ("", "."):
                    interval = fields[9].strip()
                ref = parse_coord(interval) or parse_coord(interval + "+") if interval else None
                for locus_text in fields[2].split(";"):
                    locus = parse_coord(locus_text.strip())
                    if locus is not None:
                        rows[(fields[0], locus_text.strip())] = (locus, ref)
    lifts = defaultdict(list)
    for locus, ref in rows.values():
        if ref is not None:
            lifts[locus[0]].append((locus[1], locus[2], ref[0], ref[1], ref[2]))
    for v in lifts.values():
        v.sort()
    return lifts


def side_vote(rows, start, end):
    """Length-weighted vote of anchors clipped to [start, end):
    ((contig, direction), share, {key: [(weight, reference midpoint)]})."""
    weights = defaultdict(list)
    for qs, qe, ref, rs, re_, direction, _ in rows:
        w = min(qe, end) - max(qs, start)
        if w > 0:
            weights[(ref, direction)].append((w, (rs + re_) / 2))
    total = sum(w for v in weights.values() for w, _ in v)
    if not total:
        return None, 0.0, weights
    key = max(weights, key=lambda k: sum(w for w, _ in weights[k]))
    return key, sum(w for w, _ in weights[key]) / total, weights


def safety_check(gap, rows, lifts, window, slack, min_flank, contig_length):
    """The contig must extend ``min_flank`` bp beyond the gap on both sides;
    majority votes of the placed blocks within ``window`` upstream and
    downstream, of the GenomeLift placements of the blocks in the same
    windows, and of the GenomeLift placements inside the query gap if any."""
    notes, ok = {}, True
    flanks = (gap["query_start"], contig_length - gap["query_end"])
    notes["contig_flanks"] = f"{flanks[0] / 1e6:.2f}/{flanks[1] / 1e6:.2f}Mb"
    if min(flanks) < min_flank:
        ok = False
    key = (gap["reference_contig"], gap["strand"])
    middle = (gap["reference_start"] + gap["reference_end"]) / 2
    for side, start, end, before in (("upstream", gap["query_start"] - window, gap["query_start"], True),
                                     ("downstream", gap["query_end"], gap["query_end"] + window, False)):
        winner, share, weights = side_vote(rows, start, end)
        # '+' gaps: upstream blocks lie before the reference gap, downstream after; '-' the reverse.
        expect_before = before if gap["strand"] == "+" else not before
        placed = weights.get(key, [])
        total = sum(w for w, _ in placed)
        ordered = sum(w for w, m in placed if (m <= middle) == expect_before) / total if total else 0.0
        notes[side] = (f"{winner[0]}{winner[1]}" if winner else ".") + \
            f":{share:.2f}:order{ordered:.2f}:support{total / 1e6:.2f}Mb"
        if winner != key or share <= 0.5 or ordered <= 0.5:
            ok = False
    lifted = lifts.get(gap["query_contig"], ())
    # GenomeLift of the flanking blocks: same chromosome, on the right side
    # (GenomeLift strands do not give a direction, so only order is checked).
    for side, start, end, before in (("upstream", gap["query_start"] - window, gap["query_start"], True),
                                     ("downstream", gap["query_end"], gap["query_end"] + window, False)):
        weights = defaultdict(list)
        for qs, qe, ref, rs, re_ in lifted:
            w = min(qe, end) - max(qs, start)
            if w > 0:
                weights[ref].append((w, (rs + re_) / 2))
        total = sum(w for v in weights.values() for w, _ in v)
        placed = weights.get(gap["reference_contig"], [])
        support = sum(w for w, _ in placed)
        expect_before = before if gap["strand"] == "+" else not before
        ordered = sum(w for w, m in placed if (m <= middle) == expect_before) / support if support else 0.0
        share = support / total if total else 0.0
        notes["lift_" + side] = f"{share:.2f}:order{ordered:.2f}:support{support / 1e6:.2f}Mb"
        if share <= 0.5 or ordered <= 0.5:
            ok = False
    # GenomeLift inside the query gap, when there is any, must agree too.
    agree = votes = 0
    for qs, qe, ref, rs, re_ in lifted:
        w = min(qe, gap["query_end"]) - max(qs, gap["query_start"])
        if w <= 0:
            continue
        votes += w
        if ref == gap["reference_contig"] and rs < gap["reference_end"] + slack and re_ > gap["reference_start"] - slack:
            agree += w
    notes["liftover"] = (f"up {notes['lift_upstream']} down {notes['lift_downstream']} "
                         f"inside {agree}/{votes}; contig beyond gap {notes['contig_flanks']}")
    if votes and agree / votes <= 0.5:
        ok = False
    return ok, notes


def write_outputs(prefix, gaps, summaries):
    Path(prefix).parent.mkdir(parents=True, exist_ok=True)
    details = lambda g: [str(g[k]) for k in DETAIL_FIELDS]  # noqa: E731
    score = lambda g: str(min(1000, max(g["query_gap_bp"], g["reference_gap_bp"]) // 1000))  # noqa: E731
    query_rows = sorted(gaps, key=lambda g: (g["sample"], g["query_contig"], g["query_start"]))
    with open(prefix + ".query.bed", "w") as out:
        out.write("#query_contig\tquery_start\tquery_end\tname\tscore\tstrand\t"
                  "reference_contig\treference_start\treference_end\t"
                  + "\t".join(DETAIL_FIELDS) + "\n")
        for g in query_rows:
            out.write("\t".join([
                g["query_contig"], str(g["query_start"]), str(g["query_end"]),
                g["name"], score(g), g["strand"], g["reference_contig"],
                str(g["reference_start"]), str(g["reference_end"]), *details(g)]) + "\n")
    reference_rows = sorted(gaps, key=lambda g: (
        g["reference_contig"], g["reference_start"], g["reference_end"], g["name"]))
    with open(prefix + ".reference.bed", "w") as out:
        out.write("#reference_contig\treference_start\treference_end\tname\tscore\tstrand\t"
                  "query_contig\tquery_start\tquery_end\t"
                  + "\t".join(DETAIL_FIELDS) + "\n")
        for g in reference_rows:
            out.write("\t".join([
                g["reference_contig"], str(g["reference_start"]), str(g["reference_end"]),
                g["name"], score(g), g["strand"], g["query_contig"],
                str(g["query_start"]), str(g["query_end"]), *details(g)]) + "\n")
    with open(prefix + ".summary.tsv", "w") as out:
        out.write("sample\tanchors\tunmapped_rows\tgaps\tquery_gap_bp\treference_gap_bp\t"
                  "unmapped_bp_in_gaps\n")
        for sample, (anchors, holes) in summaries.items():
            mine = [g for g in gaps if g["sample"] == sample]
            out.write("\t".join(map(str, (
                sample, anchors, holes, len(mine),
                sum(g["query_gap_bp"] for g in mine),
                sum(g["reference_gap_bp"] for g in mine),
                sum(g["unmapped_bp"] for g in mine)))) + "\n")


def write_dropped(path, dropped):
    """The gaps the filters removed, in query coordinates, with the reason."""
    with open(path, "w") as out:
        out.write("#query_contig\tquery_start\tquery_end\tname\treason\tstrand\t"
                  "reference_contig\treference_start\treference_end\tquery_gap_bp\treference_gap_bp\n")
        for g, reason in sorted(dropped, key=lambda x: (x[0]["sample"], x[0]["query_contig"],
                                                        x[0]["query_start"])):
            out.write("\t".join(map(str, (
                g["query_contig"], g["query_start"], g["query_end"], g["name"], reason, g["strand"],
                g["reference_contig"], g["reference_start"], g["reference_end"],
                g["query_gap_bp"], g["reference_gap_bp"]))) + "\n")


def main(argv=None):
    parser = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("inputs", nargs="+", help="SAMPLE.pseudolinear.tsv file(s)")
    parser.add_argument("-o", "--output", required=True, help="output prefix")
    parser.add_argument("--min-size", type=int, default=10_000,
                        help="keep gaps whose query or reference span is at least "
                             "this many bp (default: 10000)")
    parser.add_argument("--min-anchor-score", type=float, default=None,
                        help="drop gaps unless both anchors' alignment_score is at least "
                             "this (fill_graphcigartoref_gaps.py uses 1000; default: off)")
    parser.add_argument("--vote-window", type=int, default=0, metavar="BP",
                        help="safety check: placed blocks within BP upstream and downstream must "
                             "each vote (by length, > 50%%) for the gap's reference contig and "
                             "direction, on the right side of the reference gap; GenomeLift "
                             "placements in the same windows, and inside the gap if any, must "
                             "agree (needs --genomelift; default: 0, off)")
    parser.add_argument("--genomelift", help="comma-separated GenomeLift tables (genomelift.tsv,"
                                             "genomeliftfix.tsv) for --vote-window")
    parser.add_argument("--min-contig-flank", type=int, default=1_000_000,
                        help="with --vote-window, the contig must extend at least this many bp "
                             "beyond the gap on both sides (default: 1000000)")
    parser.add_argument("--liftover-slack", type=int, default=1_000_000,
                        help="a GenomeLift placement agrees when it lies within this many bp "
                             "of the reference gap (default: 1000000)")
    parser.add_argument("--fasta-query", help="the sample's FASTA: drop gaps with N in the query gap")
    parser.add_argument("--fasta-reference", help="reference FASTA: drop gaps with N in the "
                                                  "reference gap")
    args = parser.parse_args(argv)
    if args.min_size < 0:
        parser.error("--min-size must be nonnegative")
    if args.fasta_query and len(args.inputs) != 1:
        parser.error("--fasta-query needs exactly one input (one sample)")
    if args.vote_window and (not args.genomelift or not args.fasta_query or len(args.inputs) != 1):
        parser.error("--vote-window needs --genomelift, --fasta-query and exactly one input (one sample)")
    lifts = read_liftover(args.genomelift.split(",")) if args.vote_window else None
    query_fai = Fai(args.fasta_query) if args.fasta_query else None
    reference_fai = Fai(args.fasta_reference) if args.fasta_reference else None
    gaps, summaries, dropped = [], {}, []
    for path in args.inputs:
        sample = sample_name(path)
        if sample in summaries:
            parser.error(f"duplicate sample name {sample!r} from {path}")
        anchors, unmapped, inserted = read_rows(path)
        found = find_gaps(sample, anchors, unmapped, inserted, args.min_size)
        summaries[sample] = (sum(map(len, anchors.values())),
                             sum(map(len, unmapped.values())))
        total = len(found)
        before = found
        if args.min_anchor_score is not None:
            if any(a[6] is None for rows in anchors.values() for a in rows):
                raise ValueError(f"{path}: anchors without a numeric alignment_score")
            found = [g for g in found if min(float(g["left_anchor_score"]),
                                             float(g["right_anchor_score"])) >= args.min_anchor_score]
        by_score = total - len(found)
        dropped.extend((g, "anchor_score") for g in before if g not in found)
        if query_fai or reference_fai:
            for g in found:
                if query_fai:
                    g["query_n_bp"] = query_fai.count_n(g["query_contig"], g["query_start"], g["query_end"])
                if reference_fai:
                    g["reference_n_bp"] = reference_fai.count_n(
                        g["reference_contig"], g["reference_start"], g["reference_end"])
            kept = [g for g in found if g["query_n_bp"] in (".", 0) and g["reference_n_bp"] in (".", 0)]
            dropped.extend((g, "N") for g in found if g not in kept)
            by_n, found = len(found) - len(kept), kept
        else:
            by_n = 0
        by_vote = 0
        if args.vote_window:
            kept = []
            for g in found:
                ok, notes = safety_check(g, anchors.get(g["query_contig"], ()), lifts,
                                         args.vote_window, args.liftover_slack,
                                         args.min_contig_flank,
                                         query_fai.index[g["query_contig"]][0])
                g["upstream_vote"], g["downstream_vote"] = notes["upstream"], notes["downstream"]
                g["liftover_vote"] = notes["liftover"]
                if ok:
                    kept.append(g)
                else:
                    dropped.append((g, "safety"))
                    print(f"[gaps] {sample}: safety check rejects {g['name']} "
                          f"({g['query_gap_bp']} / {g['reference_gap_bp']} bp): up {notes['upstream']}, "
                          f"down {notes['downstream']}, liftover {notes['liftover']}", flush=True)
            by_vote, found = len(found) - len(kept), kept
        if args.min_anchor_score is not None or query_fai or reference_fai or args.vote_window:
            print(f"[gaps] {sample}: {total} gap(s) found; dropped {by_score} for anchor score < "
                  f"{args.min_anchor_score}, {by_n} with N, {by_vote} by the safety check", flush=True)
        print(f"[gaps] {sample}: {len(found)} gap(s) >= {args.min_size} bp, "
              f"{sum(g['query_gap_bp'] for g in found)} query bp, "
              f"{sum(g['reference_gap_bp'] for g in found)} reference bp", flush=True)
        gaps.extend(found)
    write_outputs(args.output, gaps, summaries)
    write_dropped(args.output + ".dropped.bed", dropped)
    print(f"[gaps] Wrote {args.output}.query.bed, {args.output}.reference.bed, "
          f"{args.output}.summary.tsv, {args.output}.dropped.bed", flush=True)
    return 0


if __name__ == "__main__":
    try:
        raise SystemExit(main())
    except (OSError, ValueError) as error:
        print(f"[gaps] ERROR: {error}", file=sys.stderr)
        raise SystemExit(1)
