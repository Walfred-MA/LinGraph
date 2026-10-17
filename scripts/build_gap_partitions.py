#!/usr/bin/env python3
"""Turn the large gaps of a block BED profile into CHM13-only partitions.

``prepare`` takes the reference contigs from the FASTA index, subtracts the
union of the block BED intervals, and keeps every remaining gap longer than
``--min-gap`` (default 10,000 bp). Each kept gap gets one partition folder::

    OUT/Graphs/NAME/NAME.fasta            gap template (localblocks style)
    OUT/Graphs/NAME/NAME.bed              its source-coordinate BED row
    OUT/Graphs/NAME/NAME_samples.fasta    retained gap sequence in --window records
    OUT/Graphs/NAME/NAME_samples.fasta.fai
    OUT/Graphs/NAME/NAME_samples.fasta.json  sample layout and filter settings
    OUT/Graphs/NAME/NAME_samples.selfkmer.bed  spans removed by the k-mer filter

``NAME`` is ``gapI_SAMPLE_CONTIG_START_END`` with zero-based, half-open gap
coordinates. Every partition also carries ``--flank`` bp (default 30,000) of
anchoring reference sequence on each side, clipped at the contig ends, so the
template, the sample FASTA, and the ``align`` query span the anchored interval
(``anchor_start``-``anchor_end`` in ``partitions.tsv``). The flanks are always
retained whole. The sample FASTA splits the gap into near-equal records of at
most ``--window`` bp (default 50,000, the graph block size), named
``gapI_SAMPLE_K`` with their reference coordinates. Aligner memory grows with
query length times satellite copy number, so one multi-megabase satellite
record is avoided; all windows of a gap are still self-cleaned together.

Before windowing, an exact self k-mer filter (``--self-kmer``, default 63; 0
disables it) walks the gap in genome order and removes every base covered by a
canonical k-mer that already occurred earlier in the same gap, the
keep-earlier rule of the alignment self-clean. Removed spans separated by less
than ``--kmer-merge-gap`` bp are joined, and retained fragments shorter than
``--min-fragment`` bp are dropped. ``OUT/gaps.bed`` lists every gap,
``OUT/gaps.selected.bed`` the kept ones, and ``OUT/partitions.tsv`` the kept
partitions in genome order with their record count and retained bases.

``build`` makes the condensed graph ``NAME.FA`` for each partition from the
reference sequence alone: ``minsetref_light.py`` removes the self-redundant
copies by alignment (``--min-identity``, default 95). Satellites are soft-masked,
so ``--min-score`` (uppercase + 0.3 x lowercase bases) defaults to 30 rather
than the pipeline's 100: a 100 bp fully masked fragment stays. KmerRefCleaner is not
used: it removes only non-reference sequence, and every record is reference.
The retained FASTA is the graph, as with
``minsetgraphmake.py --graph-builder direct``.

``align`` maps the original, unfiltered gap onto its condensed ``NAME.FA``
and writes the graph-CIGAR table ``NAME_align.txt``, the same format as the
pipeline's ``_align.txt``. The gap is split at the same near-equal
``--window`` grid into query records ``gapI_SAMPLE_wK`` (``NAME_query.fasta``).
``align_partition_hotspots.py --slow-rigorous --primary-aligner minimap2
--no-winnowmap`` aligns every query with minimap2 only; ``--with-winnowmap``
also runs Winnowmap on every query and pools both aligners' hits before the
graph-CIGAR conversion.
"""
from __future__ import annotations

import argparse
import collections
import csv
from concurrent.futures import ThreadPoolExecutor, as_completed
import json
import os
from pathlib import Path
import re
import shutil
import subprocess
import sys

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))

from minsetref_core import IndexedFasta  # noqa: E402

LINE_WIDTH = 80
# Byte -> 2-bit base code; 4 marks every non-A/C/G/T byte.
KMER_CODES = bytes("ACGTacgt".find(chr(i)) % 4 if chr(i) in "ACGTacgt" else 4
                   for i in range(256))
BUILD_PROTOCOL = "gap-condensed-v1"
ALIGN_PROTOCOL = "gap-align-v1"


def atomic_write(path, content):
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_name(path.name + f".tmp.{os.getpid()}")
    temporary.write_text(content)
    temporary.replace(path)


def read_fai(path):
    lengths = {}
    with open(path) as handle:
        for line in handle:
            fields = line.split("\t")
            if len(fields) >= 2 and fields[0]:
                lengths[fields[0]] = int(fields[1])
    if not lengths:
        raise ValueError(f"empty FASTA index: {path}")
    return lengths


def read_bed_intervals(path, lengths):
    intervals = {name: [] for name in lengths}
    unknown = set()
    with open(path) as handle:
        for line_number, line in enumerate(handle, 1):
            if not line.strip() or line.startswith(("#", "track", "browser")):
                continue
            fields = line.split()
            if len(fields) < 3:
                raise ValueError(f"{path}:{line_number}: fewer than three BED columns")
            contig, start, end = fields[0], int(fields[1]), int(fields[2])
            if contig not in lengths:
                unknown.add(contig)
                continue
            if not 0 <= start <= end <= lengths[contig]:
                raise ValueError(f"{path}:{line_number}: interval outside {contig}")
            intervals[contig].append((start, end))
    return intervals, unknown


def find_gaps(lengths, intervals):
    gaps = []
    for contig, length in lengths.items():
        cursor = 0
        for start, end in sorted(intervals[contig]):
            if start > cursor:
                gaps.append((contig, cursor, start))
            cursor = max(cursor, end)
        if cursor < length:
            gaps.append((contig, cursor, length))
    return gaps


def write_fasta_with_fai(path, records):
    """Write wrapped (id, description, sequence) records and their exact index."""
    path = Path(path)
    temporary = path.with_name(path.name + f".tmp.{os.getpid()}")
    index = []
    offset = 0
    with temporary.open("w") as out:
        for record_id, description, sequence in records:
            header = f">{record_id}\t{description}\n" if description else f">{record_id}\n"
            body = "".join(sequence[i:i + LINE_WIDTH] + "\n"
                           for i in range(0, len(sequence), LINE_WIDTH))
            out.write(header)
            out.write(body)
            offset += len(header.encode())
            index.append(f"{record_id}\t{len(sequence)}\t{offset}\t{LINE_WIDTH}\t{LINE_WIDTH + 1}\n")
            offset += len(body)
    temporary.replace(path)
    atomic_write(str(path) + ".fai", "".join(index))


def split_windows(start, end, window):
    """Near-equal consecutive pieces of at most ``window`` bp (0: one piece)."""
    count = 1 if window == 0 else -(-(end - start) // window)
    return [(start + (end - start) * i // count, start + (end - start) * (i + 1) // count)
            for i in range(count)]


def self_kmer_covered(sequence, k):
    """Merged spans of k-mers whose canonical form occurred earlier.

    The k-mers are exact 2-bit codes of uppercase-equivalent A/C/G/T; any
    other base breaks the k-mer run.
    """
    codes = sequence.encode().translate(KMER_CODES)
    mask = (1 << (2 * k)) - 1
    shift = 2 * (k - 1)
    seen = set()
    covered = []
    forward = reverse = 0
    run = 0
    for end, code in enumerate(codes, 1):
        if code > 3:
            run = forward = reverse = 0
            continue
        forward = ((forward << 2) | code) & mask
        reverse = (reverse >> 2) | ((3 - code) << shift)
        run += 1
        if run < k:
            continue
        kmer = forward if forward < reverse else reverse
        if kmer in seen:
            start = end - k
            if covered and covered[-1][1] >= start:
                covered[-1][1] = end
            else:
                covered.append([start, end])
        else:
            seen.add(kmer)
    return covered


def retained_intervals(sequence, k, merge_gap, min_fragment):
    """Gap-relative spans kept by the self k-mer filter (k=0: everything)."""
    length = len(sequence)
    if k == 0:
        return [(0, length)]
    removed = []
    for start, end in self_kmer_covered(sequence, k):
        if removed and start - removed[-1][1] < merge_gap:
            removed[-1][1] = end
        else:
            removed.append([start, end])
    kept, cursor = [], 0
    for start, end in removed + [[length, length]]:
        if start - cursor >= min_fragment:
            kept.append((cursor, start))
        cursor = max(cursor, end)
    return kept


def merge_spans(spans):
    merged = []
    for start, end in sorted(spans):
        if merged and start <= merged[-1][1]:
            merged[-1] = (merged[-1][0], max(merged[-1][1], end))
        else:
            merged.append((start, end))
    return merged


def window_pieces(kept, length, window):
    """Cut retained spans at the near-equal --window grid of the whole gap."""
    grid = split_windows(0, length, window)
    pieces = []
    for start, end in kept:
        for left, right in grid:
            if left < end and start < right:
                pieces.append((max(start, left), min(end, right)))
    return pieces


def partition_paths(graphs, name):
    folder = graphs / name
    return dict(
        folder=folder,
        template=folder / f"{name}.fasta",
        bed=folder / f"{name}.bed",
        samples=folder / f"{name}_samples.fasta",
        graph=folder / f"{name}.FA",
        done=folder / f"{name}.FA.done",
        work=folder / f"{name}.FA.work",
        query=folder / f"{name}_query.fasta",
        align=folder / f"{name}_align.txt",
        align_done=folder / f"{name}_align.txt.done",
    )


def read_partitions(output):
    table = Path(output) / "partitions.tsv"
    if not table.is_file():
        raise FileNotFoundError(f"{table} is missing; run prepare first")
    rows = []
    with table.open() as handle:
        header = handle.readline().rstrip("\n").split("\t")
        for line in handle:
            if line.strip():
                rows.append(dict(zip(header, line.rstrip("\n").split("\t"))))
    return rows


def prepare(args):
    reference = Path(args.reference).resolve()
    fai = Path(str(reference) + ".fai")
    if not fai.is_file():
        raise FileNotFoundError(f"{fai} is missing; run samtools faidx {reference}")
    lengths = read_fai(fai)
    for contig in args.ignore_contig:
        lengths.pop(contig, None)
    intervals, unknown = read_bed_intervals(args.bed, lengths)
    if unknown:
        print(f"[gaps] Ignored {len(unknown)} BED contig(s) absent from the reference "
              f"(or ignored): {', '.join(sorted(unknown)[:10])}", flush=True)
    empty = [contig for contig, rows in intervals.items() if not rows]
    if empty:
        print(f"[gaps] Contig(s) with no BED interval become whole-contig gaps: "
              f"{', '.join(empty)}", flush=True)

    gaps = find_gaps(lengths, intervals)
    selected = [gap for gap in gaps if gap[2] - gap[1] > args.min_gap]
    output = Path(args.output).resolve()
    graphs = output / "Graphs"
    state = dict(reference=str(reference), bed=str(Path(args.bed).resolve()),
                 sample=args.sample, min_gap=args.min_gap,
                 ignore_contig=sorted(args.ignore_contig), gaps=selected)
    config = output / "gap_partitions.json"
    if config.exists():
        if json.loads(config.read_text()) != json.loads(json.dumps(state)):
            raise ValueError(f"{output} was prepared from different inputs; use a new output directory")
    elif output.exists() and any(output.iterdir()):
        raise ValueError(f"{output} is not empty and is not a gap-partition run")

    rows = []
    for index, (contig, start, end) in enumerate(selected, 1):
        rows.append((f"gap{index}_{args.sample}_{contig}_{start}_{end}",
                     f"gap{index}_{args.sample}", contig, start, end))
    output.mkdir(parents=True, exist_ok=True)
    atomic_write(output / "gaps.bed", "".join(
        f"{c}\t{s}\t{e}\t{e - s}\n" for c, s, e in gaps))
    atomic_write(output / "gaps.selected.bed", "".join(
        f"{c}\t{s}\t{e}\t{name}\n" for name, _, c, s, e in rows))

    built = resampled = 0
    layout = dict(flank=args.flank, window=args.window, self_kmer=args.self_kmer,
                  kmer_merge_gap=args.kmer_merge_gap, min_fragment=args.min_fragment)
    summary = {}
    with IndexedFasta(str(reference)) as fasta:
        for name, prefix, contig, start, end in rows:
            paths = partition_paths(graphs, name)
            marker = Path(str(paths["samples"]) + ".json")
            try:
                saved = json.loads(marker.read_text())
            except (OSError, ValueError):
                saved = {}
            if (paths["template"].is_file() and paths["bed"].is_file()
                    and paths["samples"].is_file()
                    and Path(str(paths["samples"]) + ".fai").is_file()
                    and saved.get("layout") == layout):
                summary[name] = (saved["records"], saved["retained_bp"],
                                 saved["anchor_start"], saved["anchor_end"])
                continue
            if paths["samples"].is_file():
                # A different layout invalidates any earlier build.
                paths["done"].unlink(missing_ok=True)
                if paths["work"].exists():
                    shutil.rmtree(paths["work"])
                resampled += 1
            else:
                built += 1
            paths["folder"].mkdir(parents=True, exist_ok=True)
            # Anchoring sequence: up to --flank bp on each side, clipped at the
            # contig ends. The flanks are always retained whole.
            astart = max(0, start - args.flank)
            aend = min(lengths[contig], end + args.flank)
            sequence = fasta.fetch(contig, astart, aend)
            if len(sequence) != aend - astart:
                raise ValueError(f"{name}: fetched {len(sequence)} bp for {aend - astart} bp")
            region, core = f"{contig}:{astart}-{aend}", f"{contig}:{start}-{end}"
            write_fasta_with_fai(paths["template"], [(
                name, f"{region} block={name} core={core} role=reference "
                f"sample={args.sample} strand=+ origin=bed_gap", sequence)])
            atomic_write(paths["bed"], "\t".join(map(str, (
                contig, astart, aend, name, 1000, "+", args.sample, name,
                astart, aend, "source_coordinates"))) + "\n")
            kept = retained_intervals(sequence, args.self_kmer,
                                      args.kmer_merge_gap, args.min_fragment)
            kept = merge_spans(kept + [span for span in (
                (0, start - astart), (end - astart, aend - astart)) if span[1] > span[0]])
            pieces = window_pieces(kept, len(sequence), args.window)
            write_fasta_with_fai(paths["samples"], [
                (f"{prefix}_{k}", f"{contig}:{astart + s}-{astart + e}\tpartition={name}",
                 sequence[s:e])
                for k, (s, e) in enumerate(pieces, 1)])
            removed, cursor = [], 0
            for s, e in kept + [(len(sequence), len(sequence))]:
                if s > cursor:
                    removed.append((astart + cursor, astart + s))
                cursor = e
            atomic_write(str(paths["samples"])[:-len(".fasta")] + ".selfkmer.bed", "".join(
                f"{contig}\t{s}\t{e}\t{name}\n" for s, e in removed))
            summary[name] = (len(pieces), sum(e - s for s, e in pieces), astart, aend)
            atomic_write(marker, json.dumps(dict(
                layout=layout, records=summary[name][0], retained_bp=summary[name][1],
                anchor_start=astart, anchor_end=aend)) + "\n")
            print(f"[gaps] {name}: kept {summary[name][1]}/{aend - astart} bp "
                  f"(gap {end - start} bp + flanks) in {summary[name][0]} record(s)",
                  flush=True)

    atomic_write(output / "partitions.tsv",
                 "partition\trecords\tcontig\tstart\tend\tlength\tretained\t"
                 "anchor_start\tanchor_end\tfolder\n"
                 + "".join(
                     f"{name}\t{summary[name][0]}\t{c}\t{s}\t{e}\t{e - s}\t"
                     f"{summary[name][1]}\t{summary[name][2]}\t{summary[name][3]}\t"
                     f"{graphs / name}\n"
                     for name, _, c, s, e in rows))
    if not config.exists():
        atomic_write(config, json.dumps(state, indent=2) + "\n")
    print(f"[gaps] {len(gaps)} gaps, {sum(e - s for _, s, e in gaps)} bp; "
          f"{len(rows)} longer than {args.min_gap} bp, "
          f"{sum(e - s for _, _, _, s, e in rows)} bp, "
          f"{sum(item[1] for item in summary.values())} bp retained with flanks; "
          f"{built} partition folder(s) written, {resampled} rewritten",
          flush=True)
    print(f"[gaps] Partitions: {output / 'partitions.tsv'}", flush=True)
    return 0


def build_command(args, paths):
    command = [
        sys.executable, str(HERE / "minsetref_light.py"),
        "-i", str(paths["samples"]),
        "-o", str(paths["graph"]),
        "--output-format", "fasta",
        "-c", str(args.threads),
        "-r", args.sample,
        "-p", args.sample,
        "--min-identity", str(args.min_identity),
        "--min-score", str(args.min_score),
        "--workdir", str(paths["work"]),
        "--resume",
    ]
    if args.skip_blastn:
        command.append("--skip-blastn")
    return command


def build_one(args, row):
    name = row["partition"]
    if row.get("records") == "0":
        return name, "skipped (the k-mer filter retained no sequence)"
    paths = partition_paths(Path(row["folder"]).parent, name)
    command = build_command(args, paths)
    # Thread count does not change the result, so it is not part of the key.
    threads = command.index("-c")
    signature = json.dumps(dict(protocol=BUILD_PROTOCOL,
                                command=command[2:threads] + command[threads + 2:],
                                samples_size=paths["samples"].stat().st_size,
                                samples_mtime_ns=paths["samples"].stat().st_mtime_ns))
    if (paths["done"].is_file() and paths["done"].read_text() == signature
            and paths["graph"].is_file() and paths["graph"].stat().st_size > 0):
        return name, "reused"
    paths["done"].unlink(missing_ok=True)
    log = paths["folder"] / f"{name}.FA.log"
    print(f"[build] {name} ({row['length']} bp): log {log}", flush=True)
    with log.open("w") as handle:
        result = subprocess.run(command, stdout=handle, stderr=subprocess.STDOUT)
    if result.returncode != 0:
        return name, f"failed (exit {result.returncode}); see {log}"
    if not paths["graph"].is_file() or paths["graph"].stat().st_size == 0:
        return name, f"failed (no graph); see {log}"
    atomic_write(paths["done"], signature)
    if not args.keep_work and paths["work"].exists():
        shutil.rmtree(paths["work"])
    return name, "built"


def align_command(args, paths):
    return [
        sys.executable, str(HERE / "align_partition_hotspots.py"),
        "-i", str(paths["query"]),
        "--input-format", "fasta",
        "--graph", str(paths["graph"]),
        "-o", str(paths["align"]),
        "-s", args.sample,
        "-t", str(args.threads),
        "--timeout", str(args.timeout),
        "--slow-rigorous",
        "--primary-aligner", "minimap2",
        *([] if args.with_winnowmap else ["--no-winnowmap"]),
        "--header",
        "--alignment-payload" if args.alignment_payload else "--no-alignment-payload",
    ]


def align_one(args, row):
    name = row["partition"]
    paths = partition_paths(Path(row["folder"]).parent, name)
    if row.get("records") == "0":
        return name, "skipped (no condensed graph: the k-mer filter retained no sequence)"
    if not (paths["done"].is_file() and paths["graph"].is_file()):
        return name, "failed (NAME.FA is not built; run build first)"
    contig = row["contig"]
    start = int(row.get("anchor_start", row["start"]))
    end = int(row.get("anchor_end", row["end"]))
    prefix = name.split("_")[0] + "_" + args.sample
    pieces = split_windows(start, end, args.window)
    names = [f"{prefix}_w{k}" for k in range(1, len(pieces) + 1)]
    layout = [(n, str(e - s)) for n, (s, e) in zip(names, pieces)]
    try:
        with open(str(paths["query"]) + ".fai") as handle:
            current = [tuple(line.split("\t")[:2]) for line in handle if line.strip()]
    except FileNotFoundError:
        current = None
    if current != layout:
        with IndexedFasta(str(paths["template"])) as template:
            sequence = template.fetch(name, 0, end - start)
        if len(sequence) != end - start:
            raise ValueError(f"{name}: template has {len(sequence)} bp, expected {end - start}")
        write_fasta_with_fai(paths["query"], [
            (n, f"{args.sample}:{contig}:{s}-{e}+", sequence[s - start:e - start])
            for n, (s, e) in zip(names, pieces)])
    command = align_command(args, paths)
    threads = command.index("-t")
    stamp = lambda path: [path.stat().st_size, path.stat().st_mtime_ns]  # noqa: E731
    signature = json.dumps(dict(protocol=ALIGN_PROTOCOL,
                                command=command[2:threads] + command[threads + 2:],
                                graph=stamp(paths["graph"]), query=stamp(paths["query"])))
    if (paths["align_done"].is_file() and paths["align_done"].read_text() == signature
            and paths["align"].is_file() and paths["align"].stat().st_size > 0):
        return name, "reused"
    paths["align_done"].unlink(missing_ok=True)
    log = paths["folder"] / f"{name}_align.txt.log"
    print(f"[align] {name} ({len(pieces)} queries): log {log}", flush=True)
    with log.open("w") as handle:
        result = subprocess.run(command, stdout=handle, stderr=subprocess.STDOUT)
    if result.returncode != 0:
        return name, f"failed (exit {result.returncode}); see {log}"
    if not paths["align"].is_file() or paths["align"].stat().st_size == 0:
        return name, f"failed (no alignment table); see {log}"
    atomic_write(paths["align_done"], signature)
    return name, "aligned"


def build(args, worker=build_one, label="build"):
    rows = read_partitions(args.output)
    if args.only:
        keys = {row["partition"]: row for row in rows}
        keys.update((str(index), row) for index, row in enumerate(rows, 1))
        missing = [key for key in args.only if key not in keys]
        if missing:
            raise ValueError(f"unknown partition name/index: {', '.join(missing)}")
        rows = list({keys[key]["partition"]: keys[key] for key in args.only}.values())
    # Start the largest gaps first so they do not finish last on their own.
    rows = sorted(rows, key=lambda row: -int(row.get("retained", row["length"])))
    failures = []
    with ThreadPoolExecutor(max_workers=args.jobs) as pool:
        futures = [pool.submit(worker, args, row) for row in rows]
        for future in as_completed(futures):
            name, status = future.result()
            print(f"[{label}] {name}: {status}", flush=True)
            if status.startswith("failed"):
                failures.append(name)
    if failures:
        print(f"[{label}] {len(failures)} partition(s) failed: {', '.join(failures)}",
              file=sys.stderr, flush=True)
        return 1
    if label == "build" and not args.only:
        listing = Path(args.output).resolve() / "Graphs.list"
        atomic_write(listing, "".join(
            f"{partition_paths(Path(row['folder']).parent, row['partition'])['graph']}\n"
            for row in read_partitions(args.output) if row.get("records") != "0"))
        print(f"[build] Graph list: {listing}", flush=True)
    return 0


def find_tool(name):
    """A sibling script or binary in either the repo (scripts/, tools/) or flat layout."""
    for candidate in (HERE / name, HERE / "tools" / name, HERE.parent / "tools" / name):
        if candidate.is_file():
            return candidate
    raise FileNotFoundError(f"{name} not found beside {HERE}")


def anchor_hits(paf, min_cover=0.5):
    """Primary and secondary hits covering at least min_cover of each anchor."""
    hits = {"left": [], "right": []}
    with open(paf) as handle:
        for line in handle:
            f = line.rstrip("\n").split("\t")
            if len(f) < 12 or f[0] not in hits:
                continue
            hit = dict(qlen=int(f[1]), qstart=int(f[2]), qend=int(f[3]), strand=f[4],
                       contig=f[5], tlen=int(f[6]), tstart=int(f[7]), tend=int(f[8]),
                       matches=int(f[9]), mapq=int(f[11]))
            if hit["qend"] - hit["qstart"] >= min_cover * hit["qlen"]:
                hits[f[0]].append(hit)
    return hits


def anchor_interval(left, right):
    """Interval between a left and a right anchor hit, or None if not collinear."""
    if left["contig"] != right["contig"] or left["strand"] != right["strand"]:
        return None
    if left["strand"] == "+":
        start = left["tend"] + (left["qlen"] - left["qend"])
        end = right["tstart"] - right["qstart"]
    else:
        start = right["tend"] + right["qstart"]
        end = left["tstart"] - (left["qlen"] - left["qend"])
    start, end = max(0, start), min(left["tlen"], end)
    return (start, end) if end > start else None


def locate_region(hits, min_mapq):
    """Best collinear left/right pair (same contig and strand, in order) with at
    least one uniquely placed anchor; a duplicated anchor is resolved by the
    other one. Returns contig, start, end, strand, left mapq, right mapq."""
    best = None
    for left in hits["left"]:
        for right in hits["right"]:
            if max(left["mapq"], right["mapq"]) < min_mapq:
                continue
            interval = anchor_interval(left, right)
            if interval is None:
                continue
            key = (left["matches"] + right["matches"], -(interval[1] - interval[0]))
            if best is None or key > best[0]:
                best = (key, left, right, interval)
    if best is None:
        found = {side: sorted({f"{h['contig']}{h['strand']}(mapq {h['mapq']})" for h in found})[:3]
                 for side, found in hits.items()}
        raise ValueError(f"no collinear anchor pair with mapq >= {min_mapq} on either side; "
                         f"left hits {found['left'] or 'none'}, right hits {found['right'] or 'none'}")
    _, left, right, (start, end) = best
    return left["contig"], start, end, left["strand"], left["mapq"], right["mapq"]


def run_logged(command, log):
    with open(log, "w") as handle:
        handle.write(" ".join(map(str, command)) + "\n")
        handle.flush()
        code = subprocess.call(list(map(str, command)), stdout=handle, stderr=subprocess.STDOUT)
    if code:
        raise RuntimeError(f"command failed (exit {code}); see {log}")


COMPLEMENT = str.maketrans("ACGTacgtNn", "TGCAtgcaNn")


def fasta_records(path):
    name, header, parts = None, None, []
    with open(path) as handle:
        for line in handle:
            if line.startswith(">"):
                if name is not None:
                    yield name, header, "".join(parts)
                header = line[1:].strip()
                name, parts = header.split()[0], []
            else:
                parts.append(line.strip())
    if name is not None:
        yield name, header, "".join(parts)


def cell_intervals(cells_fasta, windows_fasta):
    """Genome interval of every cell, in file order.

    Windows carry SAMPLE:CONTIG:START-END(+|-) headers (the strand of the
    extraction). Cells are WINDOW:A-B|STEP slices of the oriented window; the
    orientation is read from the sequences themselves."""
    windows = {}
    for name, header, sequence in fasta_records(windows_fasta):
        token = header.split()[1] if len(header.split()) > 1 else ""
        match = re.fullmatch(r"(?:[^:]+:)?(.+):(\d+)-(\d+)([+-])", token)
        if match is None:
            raise ValueError(f"{windows_fasta}: no coordinates in header {header!r}")
        windows[name] = (sequence, match.group(1), int(match.group(2)), int(match.group(3)),
                         match.group(4))
    flipped, result = {}, []
    for name, _, sequence in fasta_records(cells_fasta):
        match = re.fullmatch(r"(.+):(\d+)-(\d+)\|([<>].+)", name)
        if match is None:
            raise ValueError(f"{cells_fasta}: unexpected cell name {name!r}")
        window, a, b = match.group(1), int(match.group(2)), int(match.group(3))
        wseq, contig, cs, ce, strand = windows[window]
        if window not in flipped:
            if wseq[a:b] == sequence:
                flipped[window] = False
            elif wseq.translate(COMPLEMENT)[::-1][a:b] == sequence:
                flipped[window] = True
            else:
                raise ValueError(f"{name}: sequence matches neither strand of {window}")
        length = len(wseq)
        wa, wb = (length - b, length - a) if flipped[window] else (a, b)
        start, end = (cs + wa, cs + wb) if strand == "+" else (ce - wb, ce - wa)
        result.append(dict(name=name, contig=contig, start=start, end=end, length=b - a))
    return result


def chain_cells(best_tsv, query_cells, target_cells, prefix, min_sim, max_step, max_unplaced):
    """Place each query cell at its best target cell and chain collinear runs.

    A run continues while the best target advances by 1..max_step cells in one
    direction, allowing up to max_unplaced unplaced query cells in between.
    Events between runs: DEL (target cells skipped), DUP (target cells revisited),
    INV (direction change) and INS (unplaced query cells)."""
    rows = list(csv.DictReader(open(best_tsv), delimiter="\t"))
    if len(rows) != len(query_cells):
        raise ValueError(f"{best_tsv}: {len(rows)} rows for {len(query_cells)} query cells")
    placed = [int(r["target_rank"]) if r["target_rank"] != "." and float(r["similarity"]) >= min_sim
              else None for r in rows]
    sims = [float(r["similarity"]) for r in rows]
    segments, current, unplaced = [], None, 0
    for i, rank in enumerate(placed):
        if rank is None:
            unplaced += 1
            if current is not None and unplaced > max_unplaced:
                segments.append(current)
                current = None
            continue
        unplaced = 0
        if current is not None:
            step = rank - current["t_last"]
            direction = current["dir"] or (1 if step > 0 else -1)
            if 1 <= direction * step <= max_step:
                current.update(q_last=i, t_last=rank, dir=direction, cells=current["cells"] + 1,
                               sim=current["sim"] + sims[i])
                continue
            segments.append(current)
        current = dict(q_first=i, q_last=i, t_first=rank, t_last=rank, dir=0, cells=1, sim=sims[i])
    if current is not None:
        segments.append(current)

    def span(cells, first, last):
        lo, hi = min(first, last), max(first, last)
        part = cells[lo:hi + 1]
        return part[0]["contig"], min(c["start"] for c in part), max(c["end"] for c in part), \
            sum(c["length"] for c in part)

    with open(prefix + ".chain.segments.tsv", "w") as out:
        out.write("segment\tquery_cells\tquery_contig\tquery_start\tquery_end\tquery_bp\t"
                  "target_cells\tref_contig\tref_start\tref_end\tref_bp\tdirection\tcells\tmean_similarity\n")
        for n, seg in enumerate(segments, 1):
            qc, qs, qe, qbp = span(query_cells, seg["q_first"], seg["q_last"])
            tc, ts, te, tbp = span(target_cells, seg["t_first"], seg["t_last"])
            out.write(f"{n}\t{seg['q_first']}-{seg['q_last']}\t{qc}\t{qs}\t{qe}\t{qbp}\t"
                      f"{seg['t_first']}-{seg['t_last']}\t{tc}\t{ts}\t{te}\t{tbp}\t"
                      f"{'+-'[seg['dir'] < 0] if seg['dir'] else '.'}\t{seg['cells']}\t"
                      f"{seg['sim'] / seg['cells']:.3f}\n")

    events = []
    for a, b in zip(segments, segments[1:]):
        boundary = target_cells[a["t_last"]]
        ref_point = boundary["end"] if (a["dir"] or 1) > 0 else boundary["start"]
        query_gap = [c for c in range(a["q_last"] + 1, b["q_first"])]
        if query_gap:
            qc, qs, qe, qbp = span(query_cells, query_gap[0], query_gap[-1])
            events.append(("INS", boundary["contig"], ref_point, ref_point, qbp, qc, qs, qe))
        if a["dir"] and b["dir"] and a["dir"] != b["dir"]:
            q = query_cells[b["q_first"]]
            events.append(("INV", boundary["contig"], ref_point, ref_point, 0, q["contig"], q["start"], q["start"]))
            continue
        direction = a["dir"] or b["dir"] or 1
        jump = direction * (b["t_first"] - a["t_last"])
        q = query_cells[b["q_first"]]
        if jump > 1:
            tc, ts, te, tbp = span(target_cells, a["t_last"] + direction, b["t_first"] - direction)
            events.append(("DEL", tc, ts, te, tbp, q["contig"], q["start"], q["start"]))
        elif jump <= 0:
            tc, ts, te, tbp = span(target_cells, b["t_first"], a["t_last"])
            events.append(("DUP", tc, ts, te, tbp, q["contig"], q["start"], q["start"]))
    with open(prefix + ".chain.events.tsv", "w") as out:
        out.write("type\tref_contig\tref_start\tref_end\tsize\tquery_contig\tquery_start\tquery_end\n")
        for e in events:
            out.write("\t".join(map(str, e)) + "\n")
    placed_n = sum(p is not None for p in placed)
    lines = [f"query_cells\t{len(placed)}", f"placed_cells\t{placed_n}", f"segments\t{len(segments)}",
             f"segments_ge_10_cells\t{sum(s['cells'] >= 10 for s in segments)}",
             f"cells_in_segments_ge_10\t{sum(s['cells'] for s in segments if s['cells'] >= 10)}"]
    for kind in ("DEL", "DUP", "INS", "INV"):
        sizes = [e[4] for e in events if e[0] == kind]
        lines.append(f"{kind}\t{len(sizes)}\tge_10kb\t{sum(x >= 10_000 for x in sizes)}"
                     f"\tbp_ge_10kb\t{sum(x for x in sizes if x >= 10_000)}")
    atomic_write(prefix + ".chain.summary", "\n".join(lines) + "\n")
    return placed_n, len(segments), events


def dp_events(dp_tsv, query_cells, target_cells, prefix):
    """SV calls from the glocal NW: runs of unpaired query cells are insertions
    (placed at the end of the last paired reference cell), runs of skipped
    reference cells inside the aligned span are deletions."""
    qindex = {c["name"]: c for c in query_cells}
    tindex = {c["name"]: c for c in target_cells}
    events, run, last_target, last_query = [], [], None, None

    def close():
        if not run:
            return
        kind = run[0][0]
        cells = [qindex[name] if kind == "query_gap" else tindex[name] for _, name in run]
        size = sum(c["length"] for c in cells)
        if kind == "query_gap":
            point = last_target["end"] if last_target else target_cells[0]["start"]
            contig = last_target["contig"] if last_target else target_cells[0]["contig"]
            events.append(("INS", contig, point, point, size, cells[0]["contig"],
                           min(c["start"] for c in cells), max(c["end"] for c in cells), len(cells)))
        else:
            qpoint = last_query["end"] if last_query else query_cells[0]["start"]
            qcontig = last_query["contig"] if last_query else query_cells[0]["contig"]
            events.append(("DEL", cells[0]["contig"], min(c["start"] for c in cells),
                           max(c["end"] for c in cells), size, qcontig, qpoint, qpoint, len(cells)))
        run.clear()

    for row in csv.DictReader(open(dp_tsv), delimiter="\t"):
        op = row["op"]
        if op == "pair":
            close()
            last_target, last_query = tindex[row["target_name"]], qindex[row["query_name"]]
            continue
        name = row["query_name"] if op == "query_gap" else row["target_name"]
        if run and run[0][0] != op:
            close()
        run.append((op, name))
    close()
    with open(prefix + ".dp.events.tsv", "w") as out:
        out.write("type\tref_contig\tref_start\tref_end\tsize\tquery_contig\tquery_start\tquery_end\tcells\n")
        for e in events:
            out.write("\t".join(map(str, e)) + "\n")
    return events


def haplotype_one(args, row, reference):
    """Block-level DP of one haplotype's copy of a gap against the reference cells."""
    name, sample = row["partition"], args.sample
    paths = partition_paths(Path(row["folder"]).parent, name)
    if not (paths["align"].is_file() and paths["query"].is_file() and paths["graph"].is_file()):
        raise ValueError(f"{name}: run build and align for the reference first")
    folder = paths["folder"] / "haplotypes" / sample
    folder.mkdir(parents=True, exist_ok=True)
    prefix = name.split("_")[0]
    out = lambda suffix: folder / f"{sample}{suffix}"  # noqa: E731
    contig = row["contig"]
    astart = int(row.get("anchor_start", row["start"]))
    aend = int(row.get("anchor_end", row["end"]))

    # 1. Locate the haplotype's copy with outer anchors beside the anchored interval.
    region_file = out(".region.tsv")
    if not region_file.is_file():
        with IndexedFasta(reference) as fasta:
            length = fasta.length(contig)
            spans = {"left": (max(0, astart - args.anchor_length), astart),
                     "right": (aend, min(length, aend + args.anchor_length))}
            if any(e - s < 1000 for s, e in spans.values()):
                raise ValueError(f"{name}: an outer anchor is shorter than 1 kb (contig end)")
            write_fasta_with_fai(out(".anchors.fa"), [
                (side, f"{contig}:{s}-{e}", fasta.fetch(contig, s, e))
                for side, (s, e) in spans.items()])
        paf = out(".anchors.paf")
        with open(paf, "w") as handle:
            code = subprocess.call(["minimap2", "-x", "asm5", "--secondary=yes", "-N", "10",
                                    "-t", str(args.threads),
                                    args.fasta, str(out(".anchors.fa"))], stdout=handle,
                                   stderr=subprocess.DEVNULL)
        if code:
            raise RuntimeError(f"{name} {sample}: minimap2 failed (exit {code})")
        hcontig, hstart, hend, strand, lmapq, rmapq = locate_region(anchor_hits(paf), args.min_mapq)
        atomic_write(region_file, "contig\tstart\tend\tstrand\treference\tleft_mapq\tright_mapq\n"
                     f"{hcontig}\t{hstart}\t{hend}\t{strand}\t{contig}:{astart}-{aend}\t{lmapq}\t{rmapq}\n")
    with open(region_file) as handle:
        handle.readline()
        hcontig, hstart, hend, strand = handle.readline().split("\t")[:4]
    hstart, hend = int(hstart), int(hend)

    # 2. The region in reference orientation (one record unless --window splits it).
    query = out("_query.fasta")
    layout_file = Path(str(query) + ".json")
    layout = {"window": args.window, "region": f"{hcontig}:{hstart}-{hend}{strand}"}
    try:
        same_layout = json.loads(layout_file.read_text()) == layout
    except (OSError, ValueError):
        same_layout = False
    if not (query.is_file() and same_layout):
        with IndexedFasta(args.fasta) as fasta:
            sequence = fasta.fetch(hcontig, hstart, hend)
        if strand == "-":
            sequence = sequence.translate(str.maketrans("ACGTacgtNn", "TGCAtgcaNn"))[::-1]
        records = []
        for k, (s, e) in enumerate(split_windows(0, len(sequence), args.window), 1):
            cs, ce = (hstart + s, hstart + e) if strand == "+" else (hend - e, hend - s)
            records.append((f"{prefix}_{sample}_w{k}", f"{sample}:{hcontig}:{cs}-{ce}{strand}",
                            sequence[s:e]))
        write_fasta_with_fai(query, records)
        atomic_write(layout_file, json.dumps(layout) + "\n")

    def stale(output, *inputs):
        output = Path(output)
        return not output.is_file() or any(
            Path(i).stat().st_mtime_ns > output.stat().st_mtime_ns for i in inputs)

    # 3. minimap2-only alignment to the partition graph.
    align = out("_align.txt")
    if stale(align, query):
        run_logged([sys.executable, find_tool("align_partition_hotspots.py"),
                    "-i", query, "--input-format", "fasta", "--graph", paths["graph"],
                    "-o", align, "-s", sample, "-t", args.threads, "--timeout", args.timeout,
                    "--slow-rigorous", "--primary-aligner", "minimap2", "--no-winnowmap",
                    "--header", "--no-alignment-payload"], out("_align.txt.log"))

    # 4. Cells in genome order, for the haplotype and (once) for the reference.
    cells_tool = find_tool("alignment_cells_fasta.py")
    reference_cells = paths["folder"] / f"{prefix}.cells_genome.fa"
    for alignment, sequences, cells in ((paths["align"], paths["query"], reference_cells),
                                        (align, query, out(".cells_genome.fa"))):
        if stale(cells, alignment):
            run_logged([sys.executable, cells_tool, "-a", alignment, "-q", sequences,
                        "-o", cells, "--genome-order"], str(cells) + ".log")

    # 5. KmerNorm uint16 matches within graph blocks and glocal NW.
    result = out(f".vs{args.reference_sample}")
    summary_file = Path(str(result) + ".dp.summary")
    previous = (dict(line.rstrip("\n").split("\t") for line in open(summary_file))
                if summary_file.is_file() else {})
    open_factor = args.gap_factor if args.gap_open_factor is None else args.gap_open_factor
    current = (previous.get("gap_factor") is not None
               and float(previous["gap_factor"]) == args.gap_factor
               and float(previous.get("gap_open_factor", previous["gap_factor"])) == open_factor
               and Path(str(result) + ".best.tsv").is_file()
               and not stale(summary_file, out(".cells_genome.fa"), reference_cells))
    if not current:
        run_logged([args.blockmatchdp, "-q", out(".cells_genome.fa"), "-t", reference_cells,
                    "-g", paths["graph"], "-o", result, "--max-occurrence", args.max_occurrence,
                    "--top-cells", args.top_cells, "--gap-factor", args.gap_factor,
                    "--gap-open-factor", open_factor, "-T", args.threads], str(result) + ".log")
    summary = dict(line.rstrip("\n").split("\t") for line in open(str(result) + ".dp.summary"))

    # 6. Best match per cell, chained into collinear segments and events.
    placed, segments, events = chain_cells(
        str(result) + ".best.tsv",
        cell_intervals(out(".cells_genome.fa"), query),
        cell_intervals(reference_cells, paths["query"]),
        str(result), args.chain_min_similarity, args.chain_max_step, args.chain_max_unplaced)
    big = collections.Counter(e[0] for e in events if e[4] >= 10_000)
    calls = dp_events(str(result) + ".dp.tsv", cell_intervals(out(".cells_genome.fa"), query),
                      cell_intervals(reference_cells, paths["query"]), str(result))
    big_calls = collections.Counter(e[0] for e in calls if e[4] >= 10_000)
    return name, (f"{sample} {hcontig}:{hstart}-{hend}{strand}; score {summary['score']}; "
                  f"paired {summary['paired']}/{summary['query_cells']} "
                  f"(same block {summary['paired_same_group']}); "
                  f"query unpaired {summary['query_unpaired']}; chain: {placed} cells placed, "
                  f"{segments} segments, events >= 10 kb {dict(big)}; "
                  f"DP calls {len(calls)}, >= 10 kb {dict(big_calls)}")


def haplotype(args):
    config = json.loads((Path(args.output) / "gap_partitions.json").read_text())
    rows = read_partitions(args.output)
    keys = {row["partition"]: row for row in rows}
    keys.update((str(index), row) for index, row in enumerate(rows, 1))
    missing = [key for key in args.only if key not in keys]
    if missing:
        raise ValueError(f"unknown partition name/index: {', '.join(missing)}")
    failures = 0
    for key in args.only:
        try:
            name, status = haplotype_one(args, keys[key], config["reference"])
            print(f"[haplotype] {name}: {status}", flush=True)
        except (ValueError, RuntimeError) as error:
            failures += 1
            print(f"[haplotype] {keys[key]['partition']} {args.sample}: failed: {error}",
                  file=sys.stderr, flush=True)
    return 1 if failures else 0


def build_parser():
    parser = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    commands = parser.add_subparsers(dest="command", required=True)

    p = commands.add_parser("prepare", help="find gaps and write partition folders")
    p.add_argument("-b", "--bed", required=True, help="block BED profile (BED3+)")
    p.add_argument("-r", "--reference", required=True,
                   help="uncompressed reference FASTA with adjacent .fai (CHM13)")
    p.add_argument("-o", "--output", required=True, help="dedicated output directory")
    p.add_argument("-s", "--sample", default="CHM13_h1",
                   help="reference haplotype name (default: CHM13_h1)")
    p.add_argument("--min-gap", type=int, default=10_000,
                   help="keep gaps strictly longer than this many bp (default: 10000)")
    p.add_argument("--flank", type=int, default=30_000,
                   help="anchoring reference sequence added on each side of a gap, "
                        "clipped at contig ends and always retained (default: 30000)")
    p.add_argument("--window", type=int, default=50_000,
                   help="maximum sample record length; 0 keeps each gap as one record "
                        "(default: 50000)")
    p.add_argument("--self-kmer", type=int, default=63,
                   help="remove spans covered by canonical k-mers seen earlier in the gap; "
                        "0 disables (default: 63)")
    p.add_argument("--kmer-merge-gap", type=int, default=100,
                   help="join removed spans closer than this (default: 100)")
    p.add_argument("--min-fragment", type=int, default=100,
                   help="drop retained fragments shorter than this (default: 100)")
    p.add_argument("--ignore-contig", action="append", default=[],
                   help="reference contig to skip, e.g. chrM (repeatable)")

    b = commands.add_parser("build", help="build the CHM13-only condensed .FA per partition")
    b.add_argument("-o", "--output", required=True, help="directory written by prepare")
    b.add_argument("-s", "--sample", default="CHM13_h1",
                   help="reference haplotype name used by prepare (default: CHM13_h1)")
    b.add_argument("-t", "--threads", type=int, default=8, help="threads per partition")
    b.add_argument("-j", "--jobs", type=int, default=1, help="partitions run at once")
    b.add_argument("--min-identity", type=float, default=95.0,
                   help="alignment identity %% at which a copy is redundant (default: 95)")
    b.add_argument("--min-score", type=int, default=30,
                   help="minsetref_light piece score (uppercase + 0.3 x lowercase bases); "
                        "30 keeps a 100 bp fully soft-masked fragment (default: 30)")
    b.add_argument("--only", action="append", default=[],
                   help="partition name or 1-based partitions.tsv row (repeatable), "
                        "e.g. --only $SLURM_ARRAY_TASK_ID")
    b.add_argument("--skip-blastn", action="store_true")
    b.add_argument("--keep-work", action="store_true",
                   help="keep NAME.FA.work after a successful build")

    a = commands.add_parser("align", help="map the original gap onto its condensed .FA "
                                          "with minimap2 and Winnowmap")
    a.add_argument("-o", "--output", required=True, help="directory written by prepare")
    a.add_argument("-s", "--sample", default="CHM13_h1",
                   help="reference haplotype name used by prepare (default: CHM13_h1)")
    a.add_argument("-t", "--threads", type=int, default=8, help="threads per partition")
    a.add_argument("-j", "--jobs", type=int, default=1, help="partitions run at once")
    a.add_argument("--window", type=int, default=0,
                   help="split the query into records of at most this length; "
                        "0 aligns the whole sequence as one record (default: 0)")
    a.add_argument("--timeout", type=int, default=21_600,
                   help="seconds per aligner command (default: 21600)")
    a.add_argument("--with-winnowmap", action="store_true",
                   help="also run Winnowmap on every query and pool its hits "
                        "(slow on satellites; default: minimap2 only)")
    a.add_argument("--alignment-payload", action="store_true",
                   help="attach query bases to I/X operations (default: compact CIGARs)")
    a.add_argument("--only", action="append", default=[],
                   help="partition name or 1-based partitions.tsv row (repeatable)")

    h = commands.add_parser("haplotype", help="locate one haplotype's copy of a gap, cut it "
                                              "into graph blocks and align them to the reference cells")
    h.add_argument("-o", "--output", required=True, help="directory written by prepare")
    h.add_argument("-s", "--sample", required=True, help="haplotype name, e.g. NA19240_h1")
    h.add_argument("-f", "--fasta", required=True, help="haplotype assembly FASTA (with .fai)")
    h.add_argument("--only", action="append", required=True,
                   help="partition name or 1-based partitions.tsv row (repeatable)")
    h.add_argument("--reference-sample", default="CHM13_h1")
    h.add_argument("-t", "--threads", type=int, default=8)
    h.add_argument("--anchor-length", type=int, default=50_000,
                   help="reference bp mapped on each side, outside the anchored interval "
                        "(default: 50000)")
    h.add_argument("--min-mapq", type=int, default=20)
    h.add_argument("--window", type=int, default=0,
                   help="split the region into records of at most this length; "
                        "0 aligns the whole region as one record (default: 0)")
    h.add_argument("--timeout", type=int, default=21_600)
    h.add_argument("--max-occurrence", type=int, default=255)
    h.add_argument("--top-cells", type=int, default=15_000)
    h.add_argument("--gap-factor", type=float, default=0.1,
                   help="DP cost of an unpaired cell, as a fraction of its self-match (default: 0.1)")
    h.add_argument("--gap-open-factor", type=float, default=None,
                   help="DP cost of the first unpaired cell of a run (default: same as --gap-factor)")
    h.add_argument("--chain-min-similarity", type=float, default=0.5,
                   help="minimum normalized match to place a cell (default: 0.5)")
    h.add_argument("--chain-max-step", type=int, default=5,
                   help="largest target-cell step that continues a segment (default: 5)")
    h.add_argument("--chain-max-unplaced", type=int, default=2,
                   help="unplaced query cells allowed inside a segment (default: 2)")
    h.add_argument("--blockmatchdp", default=str(HERE / "BlockMatchDP"),
                   help="BlockMatchDP executable (default: beside this script)")
    return parser


def main(argv=None):
    parser = build_parser()
    args = parser.parse_args(argv)
    if not re.fullmatch(r"[A-Za-z0-9.-]+_h[1-9][0-9]*", args.sample):
        parser.error("--sample must have the form NAME_h1")
    if args.command == "haplotype":
        return haplotype(args)
    if args.command == "prepare":
        if min(args.min_gap, args.window, args.kmer_merge_gap, args.min_fragment,
               args.flank) < 0:
            parser.error("--min-gap, --flank, --window, --kmer-merge-gap and "
                         "--min-fragment must be nonnegative")
        if not 0 <= args.self_kmer <= 64 or 0 < args.self_kmer < 15:
            parser.error("--self-kmer must be 0 or between 15 and 64")
        return prepare(args)
    if args.threads < 1 or args.jobs < 1:
        parser.error("--threads and --jobs must be positive")
    if args.command == "align":
        if args.window < 0 or args.timeout < 1:
            parser.error("--window must be >= 0 and --timeout positive")
        return build(args, align_one, "align")
    if not 0 < args.min_identity <= 100:
        parser.error("--min-identity must be in (0, 100]")
    if args.min_score < 1:
        parser.error("--min-score must be positive")
    return build(args)


if __name__ == "__main__":
    try:
        raise SystemExit(main())
    except (OSError, ValueError) as error:
        print(f"[gaps] ERROR: {error}", file=sys.stderr)
        raise SystemExit(1)
