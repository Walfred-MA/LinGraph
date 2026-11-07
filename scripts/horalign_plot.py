#!/usr/bin/env python3
"""Plot horalign graph CIGARs: the reference on top, each query below it.

Every query is reduced to aligned blocks, insertions (query only) and deletions
(reference only) of at least --min-size bp; smaller indels stay inside their
aligned block. Relocation segments are drawn as insertions, merged with any
inserted sequence next to them. All CIGARs must be against the same reference.

Columns follow the reference. A reference interval is one column, drawn as an
aligned block where a query has it and as a dotted gap where the query deletes
it; the reference track is coloured by the first query (green where that query
deletes it). Insertions at the same reference position (within --tolerance bp,
one per query) share a column as wide as the longest of them; queries without
one there show a dotted gap. Blocks are house-shaped pentagons with length
proportional to size.

Coordinates are local (0-based) unless --ref-start / --query-start give the
genome position of the first base; with --query-strand -, the query aligned is
the reverse complement of its contig interval. Rows after the first report how
many of the first query's events they share exactly (type, reference position,
size).

    python scripts/horalign_plot.py -g HG002_h2.gcigar -g HG004_h1.gcigar \\
        --label HG002_h2 --label "HG004_h1 (mother)" --ref-label "CHM13 chr6" \\
        --ref-start 58270000 --query-start 58460000 --query-start 110590000 \\
        --query-strand + --query-strand - -o HG002_h2.horalign
"""
import argparse
import bisect
import math
import re
import sys

COLORS = {"aligned": "#0072B2", "ins": "#D55E00", "del": "#009E73"}  # Wong (2011)
GAP, MARK, TEXT, AXIS = "#52514e", "#b5b4ae", "#0b0b0b", "#52514e"
HALF_HEIGHT = 0.16


def read_gcigar(path, min_size):
    """(reference name, reference length, events, query length) of one graph CIGAR.
    Events are dicts kind/r0/r1/q0/q1 in query order."""
    text = open(path).read().strip()
    chunks = re.findall(r"([<>])([^<>]*)", text)
    if not chunks or ":" not in chunks[0][1]:
        raise SystemExit(f"{path}: not a graph CIGAR (expected >NAME:...)")
    name = chunks[0][1].split(":", 1)[0]
    prims, r, q = [], 0, 0  # (kind, r, q, ref length, query length)

    def add(kind, rl, ql):
        if prims and prims[-1][0] == kind and kind in "ID":
            k, r0, q0, a, b = prims[-1]
            prims[-1] = (k, r0, q0, a + rl, b + ql)
        else:
            prims.append((kind, r, q, rl, ql))

    for index, (direction, body) in enumerate(chunks):
        named = ":" in body
        body = body.split(":", 1)[1] if named else body
        ops = [(int(n), op) for n, op in re.findall(r"(\d+)([=XMIDH])[A-Za-z]*", body)]
        if index and named:  # a relocation: inserted sequence at this place
            span = sum(n for n, op in ops if op in "=XMI")
            add("I", 0, span)
            q += span
            continue
        if direction != ">":
            raise SystemExit(f"{path}: the main path must be forward")
        for n, op in ops:
            if op == "H":
                add("H", n, 0)
            elif op in "=XM":
                add("M", n, n)
            elif op == "I":
                add("I", 0, n)
            elif op == "D":
                add("D", n, 0)
            q += n if op in "=XMI" else 0
            r += n if op in "=XMDH" else 0

    events, block = [], None
    for kind, r0, q0, rl, ql in prims:
        if kind == "M" or (kind in "ID" and max(rl, ql) < min_size):
            if block is None:
                block = dict(kind="aligned", r0=r0, q0=q0, r1=r0, q1=q0)
            block["r1"], block["q1"] = r0 + rl, q0 + ql
            continue
        if block is not None:
            events.append(block)
            block = None
        if kind != "H":
            events.append(dict(kind="ins" if kind == "I" else "del", r0=r0, q0=q0, r1=r0 + rl, q1=q0 + ql))
    if block is not None:
        events.append(block)
    for e in events:  # an "aligned" run of small indels only
        if e["kind"] == "aligned" and e["r1"] == e["r0"]:
            e["kind"] = "ins"
        elif e["kind"] == "aligned" and e["q1"] == e["q0"]:
            e["kind"] = "del"
    return name, r, events, q


def layout(tracks, ref_length, tolerance):
    """Columns in drawing order: ('ref', r0, r1) or ('ins', r, {track: event})."""
    cuts = {0, ref_length}
    insertions = []
    for t, track in enumerate(tracks):
        for e in track["events"]:
            if e["kind"] == "ins":
                insertions.append((e["r0"], t, e))
            else:
                cuts.update((e["r0"], e["r1"]))
    groups = []
    for r, t, e in sorted(insertions, key=lambda x: (x[0], x[1])):
        g = groups[-1] if groups else None
        if g is None or r - g["r"] > tolerance or t in g["members"]:
            groups.append(dict(r=r, members={t: e}))
        else:
            g["members"][t] = e
    cuts.update(g["r"] for g in groups)
    cuts = sorted(cuts)
    by_position = {}
    for g in groups:
        by_position.setdefault(g["r"], []).append(g)
    columns = []
    for a, b in zip(cuts, cuts[1:] + [None]):
        for g in by_position.get(a, []):
            width = max(e["q1"] - e["q0"] for e in g["members"].values())
            columns.append(dict(type="ins", r0=a, r1=a, width=width, members=g["members"]))
        if b is not None and b > a:
            columns.append(dict(type="ref", r0=a, r1=b, width=b - a))
    x = 0
    for c in columns:
        c["x0"], x = x, x + c["width"]
    return columns, x


def state_at(track, r0):
    """Index of the aligned/deletion event of a track covering reference r0, or None."""
    starts = track["ref_starts"]
    i = bisect.bisect_right(starts, r0) - 1
    if i >= 0:
        e = track["events"][track["ref_events"][i]]
        if e["r0"] <= r0 < e["r1"]:
            return track["ref_events"][i]
    return None


def spans(columns, track, colour_from=None):
    """Drawing spans of one row: (x0, x1, kind or None for a gap, event). Columns
    of one event merge into one span; another query's insertion splits it."""
    out = []
    for c in columns:
        if c["type"] == "ins":
            e = c["members"].get(track["index"]) if colour_from is None else None
            if e is not None:
                out.append([c["x0"], c["x0"] + e["q1"] - e["q0"], "ins", e])
                if e["q1"] - e["q0"] < c["width"]:
                    out.append([c["x0"] + e["q1"] - e["q0"], c["x0"] + c["width"], None, None])
            else:
                out.append([c["x0"], c["x0"] + c["width"], None, None])
            continue
        source = track if colour_from is None else colour_from
        i = state_at(source, c["r0"])
        e = None if i is None else source["events"][i]
        if colour_from is not None:  # the reference row
            kind = "aligned" if e is None or e["kind"] == "aligned" else "del"
        else:
            kind = None if e is None or e["kind"] == "del" else "aligned"
        last = out[-1] if out else None
        if last and last[2] == kind and last[3] is e and last[1] == c["x0"]:
            last[1] = c["x0"] + c["width"]
        else:
            out.append([c["x0"], c["x0"] + c["width"], kind, e])
    return out


def query_segments(columns, track):
    """(q0, q1, x0, x1) pieces mapping the query of a track onto the drawing."""
    out = []
    for c in columns:
        if c["type"] == "ins":
            e = c["members"].get(track["index"])
            if e is not None:
                out.append((e["q0"], e["q1"], c["x0"], c["x0"] + e["q1"] - e["q0"]))
            continue
        i = state_at(track, c["r0"])
        e = None if i is None else track["events"][i]
        if e is None or e["kind"] != "aligned":
            continue
        scale = (e["q1"] - e["q0"]) / (e["r1"] - e["r0"])
        out.append((e["q0"] + (c["r0"] - e["r0"]) * scale, e["q0"] + (c["r1"] - e["r0"]) * scale,
                    c["x0"], c["x0"] + c["width"]))
    return out


def nice_ticks(lo, hi, target=8):
    """Round tick values in [lo, hi] and their step."""
    span = max(hi - lo, 1)
    raw = span / target
    power = 10 ** math.floor(math.log10(raw))
    step = min((m * power for m in (1, 2, 2.5, 5, 10)), key=lambda s: abs(span / s - target))
    return [v * step for v in range(math.ceil(lo / step), math.floor(hi / step) + 1)], step


def tick_label(value, step):
    """Mb with as many decimals as the step needs."""
    digits = 0
    while abs(round(step / 1e6, digits) - step / 1e6) > 1e-9 and digits < 6:
        digits += 1
    return f"{value / 1e6:.{digits}f}"


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("-g", "--gcigar", action="append", required=True,
                        help="horalign graph CIGAR; repeat for more queries (first = focal)")
    parser.add_argument("-o", "--output", required=True, help="output prefix (.pdf, .png)")
    parser.add_argument("--label", action="append", default=[], help="query label, one per -g")
    parser.add_argument("--query-start", action="append", type=int, default=[],
                        help="genome position of the query's first base (contig start), one per -g")
    parser.add_argument("--query-strand", action="append", default=[], choices="+-",
                        help="strand of the query on its contig, one per -g")
    parser.add_argument("--query-axis", help="bottom axis title (first query; default: its label)")
    parser.add_argument("--ref-label", help="reference label (default: the CIGAR's reference name)")
    parser.add_argument("--ref-start", type=int, default=0, help="genome position of the reference's first base")
    parser.add_argument("--marks", help="comma-separated reference genome positions drawn as dashed lines")
    parser.add_argument("--marks-label", default="boundary", help="legend text of --marks")
    parser.add_argument("--min-size", type=int, default=50,
                        help="smallest insertion/deletion drawn as its own block (default: 50)")
    parser.add_argument("--tolerance", type=int, default=100,
                        help="insertions of different queries this close share a column (default: 100)")
    parser.add_argument("--width", type=float, default=18, help="figure width in inches (default: 18)")
    parser.add_argument("--formats", default="pdf,png", help="comma-separated formats (default: pdf,png)")
    parser.add_argument("--dpi", type=int, default=300)
    args = parser.parse_args(argv)
    n = len(args.gcigar)
    for option in ("label", "query_start", "query_strand"):
        values = getattr(args, option)
        if values and len(values) != n:
            parser.error(f"--{option.replace('_', '-')} must be given once per -g")

    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from matplotlib.legend_handler import HandlerPatch
    from matplotlib.lines import Line2D
    from matplotlib.patches import Patch, Polygon
    plt.rcParams.update({"font.family": "sans-serif",
                         "font.sans-serif": ["Helvetica", "Arial", "Liberation Sans", "DejaVu Sans"],
                         "pdf.fonttype": 42, "ps.fonttype": 42, "svg.fonttype": "none"})

    tracks, ref_name, ref_length = [], None, None
    for t, path in enumerate(args.gcigar):
        name, length, events, qlen = read_gcigar(path, args.min_size)
        if ref_name is None:
            ref_name, ref_length = name, length
        elif (name, length) != (ref_name, ref_length):
            raise SystemExit(f"{path}: reference {name} ({length} bp) differs from {ref_name} ({ref_length} bp)")
        ref_events = [i for i, e in enumerate(events) if e["kind"] != "ins"]
        ref_events.sort(key=lambda i: events[i]["r0"])
        tracks.append(dict(index=t, path=path, events=events, qlen=qlen, ref_events=ref_events,
                           ref_starts=[events[i]["r0"] for i in ref_events],
                           label=args.label[t] if args.label else re.sub(r"\.gcigar$", "", path.rsplit("/", 1)[-1]),
                           start=args.query_start[t] if args.query_start else 0,
                           strand=args.query_strand[t] if args.query_strand else None))
    columns, total = layout(tracks, ref_length, args.tolerance)

    def sv(track):
        return {(e["kind"], e["r0"], max(e["r1"] - e["r0"], e["q1"] - e["q0"]))
                for e in track["events"] if e["kind"] != "aligned"}

    focal = sv(tracks[0])
    for track in tracks:
        kinds = [e["kind"] for e in track["events"]]
        track["shared"] = len(sv(track) & focal)
        print(f"{track['label']}: {kinds.count('ins')} insertions, {kinds.count('del')} deletions"
              + (f"; {track['shared']}/{len(focal)} of the first query's events shared"
                 if track["index"] else ""))

    rows = n + 1
    fig, ax = plt.subplots(figsize=(args.width, 1.8 + 1.0 * rows))
    fig.subplots_adjust(left=0.14, right=0.99, top=1 - 0.95 / (1.8 + rows), bottom=0.75 / (1.8 + rows))
    ax.set_xlim(0, total)
    ax.set_ylim(-0.45, n + 0.45)
    ax.set_yticks([])
    for side in ("left", "right", "top"):
        ax.spines[side].set_visible(False)
    tip_max = 0.004 * total

    def pentagon(x0, x1, y, colour):
        length = x1 - x0
        tip = min(0.4 * length, tip_max)
        ax.add_patch(Polygon([(x0, y - HALF_HEIGHT), (x1 - tip, y - HALF_HEIGHT), (x1, y),
                              (x1 - tip, y + HALF_HEIGHT), (x0, y + HALF_HEIGHT)],
                             closed=True, facecolor=colour, edgecolor="white", linewidth=0.4))

    def draw(row_spans, y):
        for x0, x1, kind, _ in row_spans:
            if x1 <= x0:
                continue
            if kind is None:
                ax.plot([x0, x1], [y, y], color=GAP, linestyle=(0, (1, 2)), linewidth=0.9)
            else:
                pentagon(x0, x1, y, COLORS[kind])

    y_ref = n
    draw(spans(columns, tracks[0], colour_from=tracks[0]), y_ref)
    for track in tracks:
        draw(spans(columns, track), n - 1 - track["index"])

    def ref_x(position):
        local = position - args.ref_start
        for c in columns:
            if c["type"] == "ref" and c["r0"] <= local <= c["r1"]:
                return c["x0"] + local - c["r0"]
        return None

    if args.marks:
        for m in args.marks.split(","):
            x = ref_x(int(m))
            if x is not None:
                ax.plot([x, x], [-0.35, y_ref + 0.35], color=MARK, linestyle=(0, (4, 3)), linewidth=0.8, zorder=0)

    def axis_ticks(lo, hi, to_x):
        values, step = nice_ticks(lo, hi)
        found = [(to_x(v), v) for v in values]
        return [x for x, _ in found if x is not None], [tick_label(v, step) for x, v in found if x is not None]

    top = ax.twiny()
    top.set_xlim(ax.get_xlim())
    xs, labels = axis_ticks(args.ref_start, args.ref_start + ref_length, ref_x)
    top.set_xticks(xs, labels)
    top.set_xlabel(f"{args.ref_label or ref_name} (Mb)", color=AXIS, labelpad=6)
    for side in ("left", "right", "bottom"):
        top.spines[side].set_visible(False)
    if xs:
        top.spines["top"].set_bounds(min(xs), max(xs))

    first = tracks[0]
    segments = sorted(query_segments(columns, first))
    seg_starts = [s[0] for s in segments]

    def query_local(track, position):
        return position - track["start"] if track["strand"] != "-" else track["start"] + track["qlen"] - position

    def query_x(position):
        local = query_local(first, position)
        i = bisect.bisect_right(seg_starts, local) - 1
        for q0, q1, x0, x1 in segments[max(0, i - 1):i + 2]:
            if q0 <= local <= q1 and q1 > q0:
                return x0 + (local - q0) / (q1 - q0) * (x1 - x0)
        return None

    xs, labels = axis_ticks(first["start"], first["start"] + first["qlen"], query_x)
    order = sorted(range(len(xs)), key=lambda i: xs[i])
    ax.set_xticks([xs[i] for i in order], [labels[i] for i in order])
    if xs:
        ax.spines["bottom"].set_bounds(min(xs), max(xs))
    ax.set_xlabel(f"{args.query_axis or first['label']} (Mb)", color=AXIS, labelpad=6)
    for axis in (ax, top):
        axis.tick_params(colors=AXIS, labelsize=9, length=3)
        for spine in axis.spines.values():
            spine.set_color(AXIS)

    def side_label(y, title, start, length, strand=None):
        lo, hi = start, start + length
        text = (f"{title}{f' ({strand})' if strand else ''}\n"
                f"{lo / 1e6:.2f}-{hi / 1e6:.2f} Mb ({length / 1e6:.2f} Mb)")
        ax.text(-0.006 * total, y, text, ha="right", va="center", fontsize=9.5, color=TEXT, clip_on=False)

    side_label(y_ref, args.ref_label or ref_name, args.ref_start, ref_length)
    for track in tracks:
        y = n - 1 - track["index"]
        side_label(y, track["label"], track["start"], track["qlen"], track["strand"])
        if track["index"]:
            ax.text(-0.006 * total, y - 0.42, f"{track['shared']}/{len(focal)} events shared",
                    ha="right", va="center", fontsize=9.5, color=TEXT, clip_on=False)

    class PentagonHandler(HandlerPatch):
        def create_artists(self, legend, orig, xdescent, ydescent, width, height, fontsize, trans):
            x0, y0, tip = -xdescent, -ydescent, 0.3 * width
            return [Polygon([(x0, y0), (x0 + width - tip, y0), (x0 + width, y0 + height / 2),
                             (x0 + width - tip, y0 + height), (x0, y0 + height)], closed=True,
                            facecolor=orig.get_facecolor(), edgecolor="none", transform=trans)]

    shown = {e["kind"] for t in tracks for e in t["events"]} | {"aligned"}
    names = {"aligned": "aligned block", "ins": "insertion (query only)", "del": "deletion (reference only)"}
    handles = [Patch(facecolor=COLORS[k], label=names[k]) for k in ("aligned", "ins", "del") if k in shown]
    handles.append(Line2D([], [], color=GAP, linestyle=(0, (1, 2)), linewidth=1.4, label="gap in this sequence"))
    if args.marks:
        handles.append(Line2D([], [], color=MARK, linestyle=(0, (4, 3)), linewidth=1.4, label=args.marks_label))
    fig.legend(handles=handles, loc="upper center", ncol=len(handles), frameon=False, fontsize=11,
               handler_map={Patch: PentagonHandler()}, handlelength=2.2, columnspacing=2.2,
               bbox_to_anchor=(0.5 + 0.07, 1.0))

    written = []
    for fmt in args.formats.split(","):
        path = f"{args.output}.{fmt.strip()}"
        fig.savefig(path, dpi=args.dpi)
        written.append(path)
    print("wrote", " ".join(written), file=sys.stderr)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
