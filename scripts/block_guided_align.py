#!/usr/bin/env python3
"""Block-guided base-level alignment of a haplotype region to a reference region.

Inputs are the two sequences whose k-mer blocks BlockMatchDP aligned (-q, -t;
records are concatenated, as for the blocks) and its .dp.tsv. Block names carry
their offsets on those sequences: NAME:START-END|...

1. Every paired block (DP op "pair") is aligned with --flank bp (default 15000)
   added on both sides of the query block and of the reference block, using the
   normal-region payload aligner (graphcigartoref.minimap2_payload_ops: minimap2
   -c --eqx --for-only, best single hit or combined tandem hits, expanded to a
   global CIGAR). The alignment is trimmed with the gap-fill edge-trim policy
   (+1 per matched base, -(n + 4) per X/I/D operation; an edge is kept once the
   score exceeds 100), so flank bases stay only where they align well.
2. Merge: pieces are taken by decreasing similarity, = / (= + X + I + D). Each
   keeps only the query bases no more similar piece has claimed, and is then
   edge-trimmed again.
3. Query-adjacent pieces that move forward on the reference but overlap there
   (tandem copies): the overlap is cut from the less similar piece; its cut
   query bases become inserted sequence.
4. Query gaps between query-adjacent pieces (default): every gap is mapped to
   the whole reference with minimap2 (--relocate-options). When its best hit
   (by AS) is clear, i.e. no other hit covering at least half of the best hit's
   query span scores >= --clear-ratio (0.9) times the best, that hit becomes a
   "relocated" piece (either strand) and the rest of the gap is inserted
   sequence; otherwise the whole gap is inserted sequence. Insertions sit at the
   left neighbour's reference end (a remainder right of a relocated hit at the
   right neighbour's reference start).
   With --gap-fill, gaps are instead treated as in fill_graphcigartoref_gaps.py:
   touching on the reference: a pure insertion; both gaps positive and both <
   --max-gap (1 Mb): aligned with --anchor (100) bp anchors on both sides
   (graphcigartoref_persample._align_reference_flank_ops), anchors removed
   afterwards; otherwise (a backward jump, or a gap >= --max-gap) unresolved.
   In both modes, touching on the query with a reference gap gives a deletion of
   the part of the reference gap that no piece covers. With --anchored-ends the
   region's corners count as fixed points, so query before the first and after
   the last piece is treated as a gap too (otherwise it is left out).
5. Events: insertions and deletions of at least --min-size bp (default 50) from
   the piece CIGARs and the gaps, in reference coordinates; each relocated piece
   is also one REL event (its reference interval; size = its query span).

Outputs (PREFIX):
  .pieces.tsv   final alignment pieces in query order (kind block/relocated/ins/
                del, or gapfill/unresolved with --gap-fill), query and reference
                intervals, similarity, strand, CIGAR (a '-' piece's CIGAR walks the
                reverse complement of its query interval)
  .events.tsv   INS/DEL/REL events >= --min-size bp
  .gaps.fa, .gaps.paf  the query gaps and their minimap2 hits (default mode)
  .summary      counts
  --graph-cigar FILE    the whole alignment as one graph CIGAR (horalign.py): the
                main path is the heaviest reference-forward chain of aligned pieces;
                relocated gaps and other aligned pieces are relocation segments
                (>REF:STARTH...ENDH or <REF:...); checked to spell the query.

--no-relocation inserts every query gap instead of mapping it to the reference.
--replace-uncertain (opt-in): before step 1, expansions whose copies cannot be told
apart are re-placed (see replace_uncertain); by default the DP pairs are kept.

    python tools/block_guided_align.py -q HAP_query.fasta -t CHM13_query.fasta \
        --dp kdp_HAP_chm.dp.tsv -o HAP.bga -T 16
"""
import argparse
import bisect
import multiprocessing as mp
import os
import re
import subprocess
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
import graphcigartoref as core  # noqa: E402
import graphcigartoref_persample as persample  # noqa: E402

QUERY_OPS, REF_OPS, MATCH_OPS = set("=XMI"), set("=XMD"), set("=M")
COMPLEMENT = str.maketrans("ACGTNacgtn", "TGCANtgcan")
Q = T = ""          # sequences, shared with forked workers
ARGS = None


def records(path):
    name, parts = None, []
    with open(path) as handle:
        for line in handle:
            if line.startswith(">"):
                if name is not None:
                    yield name, "".join(parts)
                name, parts = line[1:].split()[0], []
            else:
                parts.append(line.strip())
    if name is not None:
        yield name, "".join(parts)


def span(name):
    a, b = re.search(r":(\d+)-(\d+)\|", name).groups()
    return int(a), int(b)


class Piece:
    __slots__ = ("kind", "qs", "ts", "ops", "qe", "te", "matches", "columns", "label", "strand")

    def __init__(self, kind, qs, ts, ops, label="", strand="+"):
        self.kind, self.qs, self.ts, self.ops, self.label = kind, qs, ts, ops, label
        self.strand = strand
        self.qe = qs + sum(n for n, op in ops if op in QUERY_OPS)
        self.te = ts + sum(n for n, op in ops if op in REF_OPS)
        self.matches = sum(n for n, op in ops if op in MATCH_OPS)
        self.columns = sum(n for n, op in ops if op in "=XMID")

    @property
    def similarity(self):
        return self.matches / self.columns if self.columns else 0.0


def edge_bounds(ops, minimum_score=100, gap_penalty=4):
    """fill_graphcigartoref_gaps._gapfill_edge_operation_bounds, on (n, op) pairs."""
    matches = [i for i, (n, op) in enumerate(ops) if op in MATCH_OPS and n > 0]
    if not matches:
        return None
    first, last = matches[0], matches[-1]

    def contribution(n, op):
        if op in MATCH_OPS:
            return n
        if op in "XID":
            return -(n + gap_penalty)
        return 0

    left, score = first, 0
    for i in range(first, last + 1):
        score += contribution(*ops[i])
        if score < 0:
            left, score = i + 1, 0
        elif score > minimum_score:
            break
    right, score = last + 1, 0
    for i in range(last, first - 1, -1):
        score += contribution(*ops[i])
        if score < 0:
            right, score = i, 0
        elif score > minimum_score:
            break
    return (left, right) if right > left else None


def trim(piece):
    """Edge-trim a piece; unaligned I/D operations left at either edge are dropped."""
    if piece is None:
        return None
    bounds = edge_bounds(piece.ops)
    if bounds is None:
        return None
    left, right = bounds
    while left < right and piece.ops[left][1] in "ID":
        left += 1
    while right > left and piece.ops[right - 1][1] in "ID":
        right -= 1
    if right <= left:
        return None
    qs = piece.qs + sum(n for n, op in piece.ops[:left] if op in QUERY_OPS)
    ts = piece.ts + sum(n for n, op in piece.ops[:left] if op in REF_OPS)
    return Piece(piece.kind, qs, ts, piece.ops[left:right], piece.label)


def slice_piece(piece, a, b, on_query=True):
    """Keep the bases of one sequence in [a, b); the other sequence's operations
    (D when slicing the query, I when slicing the reference) are kept only
    strictly inside the interval."""
    own, other = (QUERY_OPS, REF_OPS) if on_query else (REF_OPS, QUERY_OPS)
    other_only = "D" if on_query else "I"
    q, r = piece.qs, piece.ts
    out, start = [], None
    for n, op in piece.ops:
        pos = q if on_query else r
        if op in own:
            lo, hi = max(a, pos), min(b, pos + n)
            if hi > lo:
                if start is None:
                    shift = lo - pos
                    start = (lo, r + shift if op in other else r) if on_query else \
                            (q + shift if op in other else q, lo)
                out.append((hi - lo, op))
        elif op == other_only and a < pos < b and start is not None:
            out.append((n, op))
        if op in QUERY_OPS:
            q += n
        if op in REF_OPS:
            r += n
    if start is None:
        return None
    while out and out[-1][1] == other_only:
        out.pop()
    return Piece(piece.kind, start[0], start[1], out, piece.label)


def align_pair(task):
    index, qa, qb, ta, tb = task
    f = ARGS.flank
    q0, q1 = max(0, qa - f), min(len(Q), qb + f)
    t0, t1 = max(0, ta - f), min(len(T), tb + f)
    try:
        ops = core.minimap2_payload_ops(T[t0:t1], Q[q0:q1])
    except (RuntimeError, ValueError):
        return index, None, "no_hit"
    piece = trim(Piece("block", q0, t0, [(n, op) for n, op, _ in ops], label=f"pair{index}"))
    return index, piece, "aligned" if piece else "trimmed_away"


def align_gap(task):
    """Align one query/reference gap with anchors on both sides (normal-region gap fill)."""
    index, aq, at, la, bq, bt, ra = task
    reference, query = T[at - la:bt + ra], Q[aq - la:bq + ra]
    try:
        ops = [(op.n, op.op) for op in persample._align_reference_flank_ops(
            reference, query, reverse_from_core=False)]
    except (RuntimeError, ValueError):  # no usable alignment: keep the unaligned D/I representation
        ops = None
    sub = None
    if ops is not None:
        sub = slice_piece(Piece("gapfill", aq - la, at - la, ops), aq, bq, on_query=True)
    if sub is None or sub.ts < at or sub.te > bt:
        # anchors did not land on themselves: align the exact gaps without anchors
        try:
            ops = [(op.n, op.op) for op in persample._align_reference_flank_ops(
                T[at:bt], Q[aq:bq], reverse_from_core=False)]
        except (RuntimeError, ValueError):
            ops = [(bt - at, "D"), (bq - aq, "I")]
        return index, Piece("gapfill", aq, at, ops)
    ops = ([(sub.ts - at, "D")] if sub.ts > at else []) + sub.ops + \
          ([(bt - sub.te, "D")] if bt > sub.te else [])
    return index, Piece("gapfill", aq, at, ops)


def relocate_gaps(kept, args):
    """Map every query gap to the whole reference; a clear best hit becomes a
    relocated piece, the rest of the gap inserted sequence."""
    gaps = [(a, b) for a, b in zip(kept, kept[1:]) if b.qs > a.qe]
    if not gaps:
        return [], {"relocated": 0, "ambiguous": 0, "no_hit": 0}
    with open(args.output + ".gaps.fa", "w") as out:
        for k, (a, b) in enumerate(gaps):
            out.write(f">g{k}\n{Q[a.qe:b.qs]}\n")
    target = args.target
    if sum(1 for line in open(args.target) if line.startswith(">")) != 1:
        target = args.output + ".target.fa"
        with open(target, "w") as out:
            out.write(f">target\n{T}\n")
    with open(args.output + ".gaps.paf", "w") as paf:
        result = subprocess.run([args.minimap2, *args.relocate_options.split(), "-t", str(args.threads),
                                 target, args.output + ".gaps.fa"], stdout=paf, stderr=subprocess.PIPE,
                                text=True)
    if result.returncode:
        sys.stderr.write(result.stderr)
        raise SystemExit("minimap2 failed on the query gaps")
    if target != args.target:
        os.remove(target)
    hits = {}
    for line in open(args.output + ".gaps.paf"):
        f = line.rstrip("\n").split("\t")
        tags = dict((t[:2], t[5:]) for t in f[12:])
        if "cg" in tags:
            hits.setdefault(int(f[0][1:]), []).append(
                (int(tags.get("AS", 0)), int(f[2]), int(f[3]), f[4], int(f[7]), tags["cg"]))
    pieces, stats = [], {"relocated": 0, "ambiguous": 0, "no_hit": 0}
    for k, (a, b) in enumerate(gaps):
        h = sorted(hits.get(k, []), key=lambda x: (-x[0], x[1], x[4]))
        best = h[0] if h else None
        clear = False
        if best is not None:
            rivals = [x for x in h[1:] if min(best[2], x[2]) - max(best[1], x[1]) >= 0.5 * (best[2] - best[1])]
            clear = not rivals or rivals[0][0] < args.clear_ratio * best[0]
        if not clear:
            stats["ambiguous" if best else "no_hit"] += 1
            pieces.append(Piece("ins", a.qe, a.te, [(b.qs - a.qe, "I")], label="ambiguous" if best else "no_hit"))
            continue
        stats["relocated"] += 1
        _, qs, qe, strand, ts, cigar = best
        ops = [(int(n), op) for n, op in re.findall(r"(\d+)([=XMID])", cigar)]
        pieces.append(Piece("relocated", a.qe + qs, ts, ops, label=f"g{k}", strand=strand))
        if qs > 0:
            pieces.append(Piece("ins", a.qe, a.te, [(qs, "I")], label="relocation_flank"))
        if a.qe + qe < b.qs:
            pieces.append(Piece("ins", a.qe + qe, b.ts, [(b.qs - a.qe - qe, "I")], label="relocation_flank"))
    return pieces, stats


def subtract(a, b, covered):
    """Parts of [a, b) not inside the sorted, merged intervals ``covered``."""
    out, cursor = [], a
    k = max(0, bisect.bisect_right(covered, (a, float("inf"))) - 1)
    while k < len(covered) and covered[k][0] < b:
        s, e = covered[k]
        if e > cursor:
            if s > cursor:
                out.append((cursor, min(s, b)))
            cursor = max(cursor, e)
        k += 1
    if cursor < b:
        out.append((cursor, b))
    return out


def merge_intervals(intervals):
    merged = []
    for s, e in sorted(intervals):
        if merged and s <= merged[-1][1]:
            merged[-1][1] = max(merged[-1][1], e)
        else:
            merged.append([s, e])
    return [tuple(x) for x in merged]


def replace_uncertain(rows, best, anchor_sim=0.3, anchor_margin=0.1, min_extra=20):
    """DP pairs, as (query cell, target cell), with uncertain stretches re-placed.

    An anchor is a DP pair whose target cell is the query cell's clear best match
    anywhere (.best.tsv: similarity >= anchor_sim, runner-up at least anchor_margin
    below). Between two anchors, a stretch holding pairs and at least min_extra more
    query cells than target cells is an expansion whose copies cannot be told apart:
    which ones the DP left unpaired is arbitrary. Its cells are paired in order from
    the left and the extra query cells stay unpaired at the right end, so identical
    arrays are aligned identically. Returns (pairs, number of stretches re-placed)."""
    out, seg, replaced = [], [], 0

    def flush():
        nonlocal replaced
        qs = [r["query_cell"] for r in seg if r["op"] in ("pair", "query_gap")]
        ts = [r["target_cell"] for r in seg if r["op"] in ("pair", "target_gap")]
        if any(r["op"] == "pair" for r in seg) and len(qs) - len(ts) >= min_extra:
            out.extend(zip(qs, ts))
            replaced += 1
        else:
            out.extend((r["query_cell"], r["target_cell"]) for r in seg if r["op"] == "pair")
        seg.clear()

    for r in rows:
        if r["op"] == "pair":
            b = best.get(r["query_cell"])
            if b and b[0] == r["target_cell"] and b[1] >= anchor_sim and b[1] - b[2] >= anchor_margin:
                flush()
                out.append((r["query_cell"], r["target_cell"]))
                continue
        seg.append(r)
    flush()
    return out, replaced


def main_chain(pieces):
    """The pieces (in query order) with the most matched bases that also move
    forward on the reference without overlapping there."""
    order = sorted(range(len(pieces)), key=lambda i: (pieces[i].qs, pieces[i].ts))
    ends = sorted({p.te for p in pieces})
    tree = [(0, -1)] * (len(ends) + 1)  # Fenwick tree: best (score, piece) by reference end
    best, back = {}, {}
    for i in order:
        p, top = pieces[i], (0, -1)
        k = bisect.bisect_right(ends, p.ts)
        while k > 0:
            top, k = max(top, tree[k]), k - (k & -k)
        best[i], back[i] = top[0] + p.matches, top[1]
        k = bisect.bisect_left(ends, p.te) + 1
        while k <= len(ends):
            tree[k], k = max(tree[k], (best[i], i)), k + (k & -k)
    chain, i = [], max(((best[i], i) for i in order), default=(0, -1))[1]
    while i >= 0:
        chain.append(pieces[i])
        i = back[i]
    return chain[::-1]


def add_op(ops, n, op, payload=""):
    """Append (n, op, payload), merging with an equal last op. An X payload is a
    (reference, query) pair; an I payload the query bases."""
    if n <= 0:
        return
    if op == "X" and not payload:
        payload = ("", "")
    if ops and ops[-1][1] == op:
        last = ops[-1][2]
        ops[-1] = (ops[-1][0] + n, op, (last[0] + payload[0], last[1] + payload[1])
                   if op == "X" else last + payload)
    else:
        ops.append((n, op, payload))


def piece_payload_ops(p):
    """(n, op, payload) of a piece along the query. A '-' piece walks the
    reverse-complemented reference."""
    out = []
    ops = p.ops if p.strand == "+" else p.ops[::-1]
    q, r = p.qs, (p.ts if p.strand == "+" else p.te)
    for n, op in ops:
        if op in "=XM":
            reference = T[r:r + n] if p.strand == "+" else T[r - n:r].translate(COMPLEMENT)[::-1]
            query = Q[q:q + n]
            if op == "=":
                add_op(out, n, "=")
            elif op == "X":
                add_op(out, n, "X", (reference, query))
            else:
                for a, b in zip(reference, query):
                    add_op(out, 1, "=") if a.upper() == b.upper() else add_op(out, 1, "X", (a, b))
        elif op == "I":
            add_op(out, n, "I", Q[q:q + n])
        elif op == "D":
            add_op(out, n, "D")
        q += n if op in QUERY_OPS else 0
        r += (n if p.strand == "+" else -n) if op in REF_OPS else 0
    return out


def format_ops(ops):
    """X payloads are written as the reference bases, then the query bases."""
    return "".join(f"{n}{op}{''.join(payload) if op == 'X' else payload}" for n, op, payload in ops)


def graph_cigar(chain, relocations, name):
    """One graph CIGAR of the whole query on the whole reference: the main path
    (>NAME:...) walks the chain, with D for reference it skips and I for query
    nothing aligned; each relocation is its own segment >NAME:STARTH...ENDH (or
    <NAME:... on the reverse strand) at its query position, after which >... resumes
    the main path."""
    length, chunks, main = len(T), [], []

    def add(n, op, payload=""):
        add_op(main, n, op, payload)

    def flush(head=False):
        if main or (head and not chunks):
            prefix = ">" if chunks else f">{name}:" + ("0H" if not main else "")
            chunks.append(prefix + format_ops(main))
            main.clear()

    relocations = sorted(relocations, key=lambda p: p.qs)
    q = r = k = 0
    for c in chain + [Piece("end", len(Q), length, [])]:
        add(c.ts - r, "D")
        while k < len(relocations) and relocations[k].qs < c.qs:
            p = relocations[k]
            if p.qs < q or p.qe > c.qs:
                raise ValueError(f"relocated piece {p.qs}-{p.qe} overlaps the main path on the query")
            add(p.qs - q, "I", Q[q:p.qs])
            flush(head=True)
            body = format_ops(piece_payload_ops(p))
            chunks.append(f">{name}:{p.ts}H{body}{length - p.te}H" if p.strand == "+" else
                          f"<{name}:{length - p.te}H{body}{p.ts}H")
            q, k = p.qe, k + 1
        add(c.qs - q, "I", Q[q:c.qs])
        for op in piece_payload_ops(c):
            add(*op)
        q, r = c.qe, c.te
    flush()
    return "".join(chunks)


def decode_graph_cigar(text):
    """(query, main-path reference length) spelled by a graph CIGAR, for checking."""
    query, main_length, reverse = [], 0, T.translate(COMPLEMENT)[::-1]
    for index, (direction, body) in enumerate(re.findall(r"([<>])([^<>]*)", text)):
        named = ":" in body
        body = body.split(":", 1)[1] if named else body
        is_main = index == 0 or not named
        sequence = T if direction == ">" else reverse
        pos = main_length if is_main else 0
        for n, op, payload in re.findall(r"(\d+)([=XIDH])([A-Za-z]*)", body):
            n = int(n)
            if op == "=":
                query.append(sequence[pos:pos + n])
            elif op == "X":
                if payload[:n].upper() != sequence[pos:pos + n].upper():
                    raise ValueError("X payload does not match the reference")
                query.append(payload[n:])
            elif op == "I":
                query.append(payload)
            pos += n if op in "=XDH" else 0
        if is_main:
            main_length = pos
    return "".join(query), main_length


def main(argv=None):
    global Q, T, ARGS
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("-q", "--query", required=True, help="haplotype FASTA (the DP query)")
    parser.add_argument("-t", "--target", required=True, help="reference FASTA (the DP target)")
    parser.add_argument("--dp", required=True, help="BlockMatchDP .dp.tsv")
    parser.add_argument("-o", "--output", required=True, help="output prefix")
    parser.add_argument("-T", "--threads", type=int, default=8)
    parser.add_argument("--flank", type=int, default=15_000)
    parser.add_argument("--anchor", type=int, default=100)
    parser.add_argument("--max-gap", type=int, default=1_000_000)
    parser.add_argument("--min-size", type=int, default=50)
    parser.add_argument("--gap-fill", action="store_true",
                        help="fill gaps with anchored alignments instead of relocating them")
    parser.add_argument("--minimap2", default="minimap2")
    parser.add_argument("--relocate-options", default="-cx asm20 --eqx --secondary=yes -N 50 -p 0.5",
                        help="minimap2 options for mapping query gaps to the whole reference")
    parser.add_argument("--anchored-ends", action="store_true",
                        help="the region's two ends are anchored (e.g. flanks beside known anchors): "
                             "query before the first and after the last piece is a gap too")
    parser.add_argument("--clear-ratio", type=float, default=0.9,
                        help="a best hit is clear when every overlapping rival scores below this "
                             "fraction of it (default: 0.9)")
    parser.add_argument("--no-relocation", action="store_true",
                        help="do not map query gaps to the reference; every gap is inserted sequence")
    parser.add_argument("--graph-cigar",
                        help="also write the whole alignment as one graph CIGAR to this file")
    parser.add_argument("--reference-name",
                        help="reference name in the graph CIGAR (default: the first -t record)")
    parser.add_argument("--replace-uncertain", dest="replace_uncertain", action="store_true", default=False,
                        help="re-place expansions whose copies cannot be told apart: between clear-best "
                             "anchor pairs, pair in order from the left (off by default: it lowered trio "
                             "consistency on 4 trios)")
    parser.add_argument("--no-replace-uncertain", dest="replace_uncertain", action="store_false",
                        help="use BlockMatchDP's pairs as they are (default)")
    ARGS = args = parser.parse_args(argv)

    Q = "".join(s for _, s in records(args.query))
    T = "".join(s for _, s in records(args.target))
    pairs, rows, names = [], [], {}
    with open(args.dp) as handle:
        header = handle.readline().rstrip("\n").split("\t")
        for line in handle:
            row = dict(zip(header, line.rstrip("\n").split("\t")))
            if row["op"] == "pair":
                pairs.append((len(pairs), *span(row["query_name"]), *span(row["target_name"])))
            for side in ("query", "target"):
                if row[f"{side}_name"] not in ("", ".", "-"):
                    row[f"{side}_cell"] = int(row[f"{side}_cell"])
                    names[row[f"{side}_cell"]] = row[f"{side}_name"]
            rows.append(row)
    replaced = 0
    if args.replace_uncertain:
        best_path = re.sub(r"\.dp\.tsv$", "", args.dp) + ".best.tsv"
        if not os.path.isfile(best_path):
            raise SystemExit(f"{best_path} not found (BlockMatchDP output); use --no-replace-uncertain")
        best = {}
        with open(best_path) as handle:
            header = handle.readline().rstrip("\n").split("\t")
            for line in handle:
                row = dict(zip(header, line.rstrip("\n").split("\t")))
                if row["best_target_cell"] != ".":
                    best[int(row["query_cell"])] = (int(row["best_target_cell"]), float(row["similarity"]),
                                                    float(row["second_similarity"]))
        cells, replaced = replace_uncertain(rows, best)
        if replaced:
            pairs = [(i, *span(names[q]), *span(names[t])) for i, (q, t) in enumerate(cells)]
    print(f"[bga] {len(pairs)} paired blocks" +
          (f" ({replaced} uncertain stretches re-placed)" if args.replace_uncertain else "") +
          f"; query {len(Q)} bp, reference {len(T)} bp",
          file=sys.stderr, flush=True)

    # 1. Paired blocks with flanks.
    status = {"aligned": 0, "no_hit": 0, "trimmed_away": 0}
    raw = []
    context = mp.get_context("fork")
    with context.Pool(args.threads) as pool:
        for index, piece, state in pool.imap_unordered(align_pair, pairs, chunksize=4):
            status[state] += 1
            if piece is not None:
                raw.append((index, piece))
            done = sum(status.values())
            if done % 200 == 0:
                print(f"[bga] aligned {done}/{len(pairs)}", file=sys.stderr, flush=True)
    raw.sort()

    # 2. Merge by similarity: each piece keeps query bases no more similar piece claimed.
    claimed, kept = [], []
    for index, piece in sorted(raw, key=lambda x: (-x[1].similarity, -x[1].matches, x[0])):
        for a, b in subtract(piece.qs, piece.qe, claimed):
            sub = trim(slice_piece(piece, a, b, on_query=True))
            if sub is not None:
                kept.append(sub)
                bisect.insort(claimed, (sub.qs, sub.qe))
    kept.sort(key=lambda p: (p.qs, p.ts))

    # 3. Forward reference overlaps between query neighbours: cut the less similar piece.
    cuts, i = 0, 0
    while i + 1 < len(kept):
        a, b = kept[i], kept[i + 1]
        if b.ts < a.te and b.ts >= a.ts:
            cuts += 1
            if a.similarity < b.similarity:
                new = trim(slice_piece(a, a.ts, b.ts, on_query=False))
                if new is None:
                    del kept[i]
                    i = max(0, i - 1)
                    continue
                kept[i] = new
            else:
                new = trim(slice_piece(b, a.te, b.te, on_query=False))
                if new is None:
                    del kept[i + 1]
                    continue
                kept[i + 1] = new
            continue
        i += 1

    # 4. Gaps between query neighbours.
    extra, tasks, unresolved = [], [], []
    relocation = {"relocated": 0, "ambiguous": 0, "no_hit": 0}
    # With --anchored-ends, zero-length end points at the region's corners
    # bound the gaps before the first and after the last piece.
    bounded = ([Piece("end", 0, 0, [])] + kept + [Piece("end", len(Q), len(T), [])]
               if args.anchored_ends else kept)
    if args.no_relocation and not args.gap_fill:
        extra = [Piece("ins", a.qe, a.te, [(b.qs - a.qe, "I")], label="no_relocation")
                 for a, b in zip(bounded, bounded[1:]) if b.qs > a.qe]
    elif not args.gap_fill:
        extra, relocation = relocate_gaps(bounded, args)
    covered = merge_intervals([(p.ts, p.te) for p in kept] +
                              [(p.ts, p.te) for p in extra if p.kind == "relocated"])
    for k, (a, b) in enumerate(zip(bounded, bounded[1:])):
        qgap, tgap = b.qs - a.qe, b.ts - a.te
        if qgap == 0 and tgap > 0:
            for s, e in subtract(a.te, b.ts, covered):
                extra.append(Piece("del", a.qe, s, [(e - s, "D")]))
        elif not args.gap_fill:
            continue
        elif tgap < 0 or (qgap > 0 and tgap > 0 and max(qgap, tgap) >= args.max_gap):
            if qgap > 0:
                unresolved.append(Piece("unresolved", a.qe, a.te, [(qgap, "I")]))
        elif tgap == 0 and qgap > 0:
            extra.append(Piece("ins", a.qe, a.te, [(qgap, "I")]))
        elif qgap > 0 and tgap > 0:
            la = min(args.anchor, a.qe - a.qs, a.te - a.ts)
            ra = min(args.anchor, b.qe - b.qs, b.te - b.ts)
            tasks.append((k, a.qe, a.te, la, b.qs, b.ts, ra))
    if tasks:
        with context.Pool(args.threads) as pool:
            for _, piece in pool.imap_unordered(align_gap, tasks, chunksize=4):
                extra.append(piece)
    final = sorted(kept + extra + unresolved, key=lambda p: (p.qs, p.ts, p.kind))

    # 5. Events.
    events = []
    for p in final:
        if p.kind == "unresolved":
            continue
        if p.kind == "relocated" and p.qe - p.qs >= args.min_size:
            events.append(("REL", p.ts, p.te, p.qe - p.qs, p.qs, p.qe, p.kind))
        # a '-' piece's CIGAR walks its query interval from the end
        q, r, step = (p.qs, p.ts, 1) if p.strand == "+" else (p.qe, p.ts, -1)
        for n, op in p.ops:
            if op == "I" and n >= args.min_size:
                events.append(("INS", r, r, n, *sorted((q, q + step * n)), p.kind))
            elif op == "D" and n >= args.min_size:
                events.append(("DEL", r, r + n, n, q, q, p.kind))
            q += step * n if op in QUERY_OPS else 0
            r += n if op in REF_OPS else 0
    events.sort(key=lambda e: (e[1], e[2], e[0]))

    with open(args.output + ".pieces.tsv", "w") as out:
        out.write("kind\tquery_start\tquery_end\tref_start\tref_end\tsimilarity\tmatches\tsource\t"
                  "strand\tcigar\n")
        for p in final:
            cigar = "".join(f"{n}{op}" for n, op in p.ops)
            out.write(f"{p.kind}\t{p.qs}\t{p.qe}\t{p.ts}\t{p.te}\t{p.similarity:.4f}\t{p.matches}\t"
                      f"{p.label or '.'}\t{p.strand}\t{cigar}\n")
    with open(args.output + ".events.tsv", "w") as out:
        out.write("type\tref_start\tref_end\tsize\tquery_start\tquery_end\tsource\n")
        for e in events:
            out.write("\t".join(map(str, e)) + "\n")
    bp = lambda kind: sum(p.qe - p.qs for p in final if p.kind == kind)  # noqa: E731
    first, last = (final[0].qs, final[-1].qe) if final else (0, 0)
    summary = [
        ("paired_blocks", len(pairs)), *status.items(),
        *([("uncertain_stretches_replaced", replaced)] if args.replace_uncertain else []),
        ("block_pieces", sum(p.kind == "block" for p in final)),
        ("reference_overlap_cuts", cuts),
        ("gap_alignments", len(tasks)),
        ("gaps_relocated", relocation["relocated"]), ("gaps_ambiguous", relocation["ambiguous"]),
        ("gaps_no_hit", relocation["no_hit"]),
        ("query_span", f"{first}-{last}"),
        ("query_bp_block", bp("block")), ("query_bp_relocated", bp("relocated")),
        ("query_bp_gapfill", bp("gapfill")),
        ("query_bp_ins", bp("ins")), ("query_bp_unresolved", bp("unresolved")),
        ("block_matches", sum(p.matches for p in final if p.kind == "block")),
        ("relocated_matches", sum(p.matches for p in final if p.kind == "relocated")),
        ("gapfill_matches", sum(p.matches for p in final if p.kind == "gapfill")),
        ("events_ins", sum(e[0] == "INS" for e in events)),
        ("events_del", sum(e[0] == "DEL" for e in events)),
        ("events_rel", sum(e[0] == "REL" for e in events)),
        ("events_ins_bp", sum(e[3] for e in events if e[0] == "INS")),
        ("events_del_bp", sum(e[3] for e in events if e[0] == "DEL")),
    ]
    if args.graph_cigar:
        # Main path: the heaviest reference-forward chain of aligned pieces; the
        # other aligned pieces and the relocated gaps are relocation segments.
        aligned = [p for p in final if p.kind in ("block", "gapfill")]
        chain = main_chain(aligned)
        on_chain = set(map(id, chain))
        relocations = [p for p in final if p.kind == "relocated" or
                       (p.kind in ("block", "gapfill") and id(p) not in on_chain)]
        name = args.reference_name or next(records(args.target))[0]
        text = graph_cigar(chain, relocations, name)
        query, main_length = decode_graph_cigar(text)
        if query.upper() != Q.upper() or main_length != len(T):
            raise SystemExit("graph CIGAR check failed: it does not spell the query on the reference")
        with open(args.graph_cigar, "w") as out:
            out.write(text + "\n")
        summary += [("graph_cigar_relocations", len(relocations)),
                    ("graph_cigar_offchain_blocks", len(relocations) - relocation["relocated"])]
    with open(args.output + ".summary", "w") as out:
        out.write("".join(f"{k}\t{v}\n" for k, v in summary))
    print("[bga] " + "; ".join(f"{k} {v}" for k, v in summary), file=sys.stderr)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
