#!/usr/bin/env python3
"""Reference-only coordinate hardness: local complexity, repeats and anchors.

Default seeding is minimap2's ordinary (non-HPC) k=15, w=10. No mapping
preset is applied. The score is a heuristic, not MAPQ or an error probability.
See docs/reference_hardness.md for definitions and limitations.
"""
from collections import Counter
import argparse
import gzip
import hashlib
import json
import math
from pathlib import Path
import re
from typing import Iterator, NamedTuple


_DNA_RUN = re.compile(r"[ACGTU]+", re.IGNORECASE)
_INVALID_BASE = re.compile(r"[^ACGTURYSWKMBDHVN]", re.IGNORECASE)
_NT = {base: i for i, base in enumerate("ACGT")}
_NT["U"] = 3
_MAX64 = (1 << 64) - 1
_MIN_COMPLEXITY_RUN = 32


class Minimizer(NamedTuple):
    key: int
    start: int
    end: int
    strand: int


def _check_sketch_parameters(k, w):
    if isinstance(k, bool) or not isinstance(k, int) or not 1 <= k <= 28:
        raise ValueError("k must be an integer from 1 to 28")
    if isinstance(w, bool) or not isinstance(w, int) or not 1 <= w <= 255:
        raise ValueError("w must be an integer from 1 to 255")


def _hash_kmer(key, mask):
    key = (~key + (key << 21)) & mask
    key ^= key >> 24
    key = (key * 265) & mask
    key ^= key >> 14
    key = (key * 21) & mask
    key ^= key >> 28
    return (key + (key << 31)) & mask


def minimap2_minimizers(sequence: str, k: int = 15, w: int = 10) -> Iterator[Minimizer]:
    """Yield ordinary minimap2 minimizers, including tied positions.

    Coordinates are zero-based, half-open on the original sequence. Lowercase
    is unmasked, U is T, and ambiguous bases break seeds. Reverse-complement
    copies share keys. This follows minimap2 mm_sketch, including its terminal
    partial-window and symmetric-k-mer behavior; HPC (-H) is not supported.

    Adapted from minimap2 sketch.c, copyright (c) 2018- Dana-Farber Cancer
    Institute, (c) 2017-2018 Broad Institute, Inc., under the MIT license reproduced
    in docs/reference_hardness.md.
    """
    _check_sketch_parameters(k, w)
    mask = (1 << (2 * k)) - 1
    reverse_shift = 2 * (k - 1)
    empty = (_MAX64, -1, 0)
    ring = [empty] * w
    minimum = empty
    ring_pos = minimum_pos = valid = forward = reverse = 0

    def unpack(item):
        value, last, strand = item
        return Minimizer(value >> 8, last + 1 - (value & 255), last + 1, strand)

    for pos, base in enumerate(sequence.upper()):
        code = _NT.get(base, 4)
        item = empty
        if code < 4:
            span = min(valid + 1, k)
            forward = ((forward << 2) | code) & mask
            reverse = (reverse >> 2) | ((3 ^ code) << reverse_shift)
            if forward == reverse:
                continue
            strand = int(forward > reverse)
            valid += 1
            if valid >= k:
                key = _hash_kmer(reverse if strand else forward, mask)
                item = ((key << 8) | span, pos, strand)
        else:
            valid = 0
        ring[ring_pos] = item
        if valid == w + k - 1 and minimum[0] != _MAX64:
            for j in (*range(ring_pos + 1, w), *range(ring_pos)):
                if ring[j][0] == minimum[0] and ring[j][1:] != minimum[1:]:
                    yield unpack(ring[j])
        if item[0] <= minimum[0]:
            if valid >= w + k and minimum[0] != _MAX64:
                yield unpack(minimum)
            minimum, minimum_pos = item, ring_pos
        elif ring_pos == minimum_pos:
            if valid >= w + k - 1 and minimum[0] != _MAX64:
                yield unpack(minimum)
            minimum = empty
            for j in (*range(ring_pos + 1, w), *range(ring_pos + 1)):
                if minimum[0] >= ring[j][0]:
                    minimum, minimum_pos = ring[j], j
            if valid >= w + k - 1 and minimum[0] != _MAX64:
                for j in (*range(ring_pos + 1, w), *range(ring_pos + 1)):
                    if ring[j][0] == minimum[0] and ring[j][1:] != minimum[1:]:
                        yield unpack(ring[j])
        ring_pos = (ring_pos + 1) % w
    if minimum[0] != _MAX64:
        yield unpack(minimum)


def _read_fasta(path):
    """Read one record at a time; fail on malformed or duplicate records."""
    with open(path, "rb") as probe:
        compressed = probe.read(2) == b"\x1f\x8b"
    opener = gzip.open if compressed else open
    seen = set()
    name = None
    parts = []

    def record():
        sequence = "".join(parts).upper()
        if not sequence:
            raise ValueError(f"empty FASTA record: {name}")
        bad = _INVALID_BASE.search(sequence)
        if bad:
            raise ValueError(f"invalid FASTA base {bad[0]!r} in {name} at {bad.start()}")
        return name, sequence

    with opener(path, "rt", encoding="ascii") as handle:
        for line_number, line in enumerate(handle, 1):
            line = line.strip()
            if not line:
                continue
            if line.startswith(">"):
                if name is not None:
                    yield record()
                fields = line[1:].split()
                if not fields:
                    raise ValueError(f"missing FASTA name at line {line_number}")
                name = fields[0]
                if name in seen:
                    raise ValueError(f"duplicate FASTA name: {name}")
                seen.add(name)
                parts = []
            else:
                if name is None:
                    raise ValueError(f"sequence before FASTA header at line {line_number}")
                parts.append("".join(line.split()))
    if name is None:
        raise ValueError("FASTA contains no records")
    yield record()


def _low_complexity_tracts(sequence, start, end, window, entropy_cutoff, length_scale):
    """Yield unions of low-entropy windows inside one uninterrupted DNA run.

    Trinucleotide entropy uses all overlapping words, with a window step of
    one base. Runs of 32..window-1 bp use their whole length as one window;
    shorter runs are unassessed. Tract severity averages qualifying windows,
    and tract risk also accounts for the length of their union.
    """
    if end - start < _MIN_COMPLEXITY_RUN:
        return
    window = min(window, end - start)
    n = window - 2
    terms = [c * math.log2(c) if c else 0.0 for c in range(n + 1)]
    counts = Counter(sequence[i:i + 3] for i in range(start, start + n))
    subtotal = math.fsum(terms[c] for c in counts.values())
    tract = None

    def finish():
        left, right, severity_sum, count, minimum_entropy = tract
        size = right - left
        severity = severity_sum / count
        effective_length = severity * size
        return {"start": left, "end": right, "length_bp": size,
                "mean_severity": severity, "min_entropy_bits": minimum_entropy,
                "qualifying_windows": count, "window_bp": window,
                "severity_weighted_length_bp": effective_length,
                "risk": effective_length / (effective_length + length_scale)}

    for pos in range(start, end - window + 1):
        if pos != start:
            old = sequence[pos - 1:pos + 2]
            new = sequence[pos + window - 3:pos + window]
            count = counts[old]
            subtotal += terms[count - 1] - terms[count]
            counts[old] -= 1
            count = counts[new]
            subtotal += terms[count + 1] - terms[count]
            counts[new] += 1
            if (pos - start) % 4096 == 0:
                subtotal = math.fsum(terms[c] for c in counts.values())
        entropy = max(0.0, min(6.0, math.log2(n) - subtotal / n))
        # A tolerance prevents roundoff from creating boundary-only tracts.
        if entropy >= entropy_cutoff - 1e-12:
            continue
        severity = (entropy_cutoff - entropy) / entropy_cutoff
        if tract is not None and pos > tract[1]:
            yield finish()
            tract = None
        if tract is None:
            tract = [pos, pos + window, severity, 1, entropy]
        else:
            tract[1] = pos + window
            tract[2] += severity
            tract[3] += 1
            tract[4] = min(tract[4], entropy)
    if tract is not None:
        yield finish()


class _Measurements:
    def __init__(self):
        self.length = 0
        self.callable_bases = 0
        self.copies = Counter()
        self.gaps = Counter()
        self.gap_penalty_bases = 0.0
        self.seedless_bases = 0
        self.complexity_evaluated_bases = 0
        self.low_complexity_bases = 0
        self.low_complexity_penalty_bases = 0.0
        self.low_complexity_peak_risk = 0.0
        self.low_complexity_tract_count = 0
        self.longest_low_complexity_tract = 0

    def add(self, other):
        self.length += other.length
        self.callable_bases += other.callable_bases
        self.copies.update(other.copies)
        self.gaps.update(other.gaps)
        self.gap_penalty_bases += other.gap_penalty_bases
        self.seedless_bases += other.seedless_bases
        self.complexity_evaluated_bases += other.complexity_evaluated_bases
        self.low_complexity_bases += other.low_complexity_bases
        self.low_complexity_penalty_bases += other.low_complexity_penalty_bases
        self.low_complexity_peak_risk = max(self.low_complexity_peak_risk, other.low_complexity_peak_risk)
        self.low_complexity_tract_count += other.low_complexity_tract_count
        self.longest_low_complexity_tract = max(self.longest_low_complexity_tract, other.longest_low_complexity_tract)

    def report(self, long_gap, high_copy, seed_weight, local_weight, hotspot_weight):
        n = sum(self.copies.values())
        unique = self.copies[1]
        repeated = n - unique
        callable_bases = self.callable_bases
        ambiguity = sum(count * (1 - 1 / copies)
                        for copies, count in self.copies.items()) / n if n else None
        extent = self.gap_penalty_bases / callable_bases if callable_bases else None
        seed_score = (100 * (seed_weight * ambiguity + (1 - seed_weight) * extent)
                      if n and callable_bases else None)
        mean_risk = self.low_complexity_penalty_bases / callable_bases if callable_bases else None
        peak_risk = self.low_complexity_peak_risk if callable_bases else None
        local_risk = ((1 - hotspot_weight) * mean_risk + hotspot_weight * peak_risk
                      if callable_bases else None)
        hardness = ((1 - local_weight) * seed_score + 100 * local_weight * local_risk
                    if seed_score is not None else None)
        gap_bases = sum(size * count for size, count in self.gaps.items())
        cumulative = n50 = 0
        for size, count in sorted(self.gaps.items(), reverse=True):
            cumulative += size * count
            if cumulative * 2 >= gap_bases:
                n50 = size
                break
        long_bases = sum(size * count for size, count in self.gaps.items() if size >= long_gap)
        return {
            "status": "ok" if n else ("no_seeds" if callable_bases else "no_callable_bases"),
            "hardness_score": hardness,
            "seed_hardness_score": seed_score,
            "local_complexity_risk": local_risk,
            "low_complexity_mean_risk": mean_risk,
            "low_complexity_peak_risk": peak_risk,
            "low_complexity_base_fraction": self.low_complexity_bases / callable_bases if callable_bases else None,
            "low_complexity_tract_count": self.low_complexity_tract_count,
            "longest_low_complexity_tract_bp": self.longest_low_complexity_tract,
            "complexity_evaluated_bases": self.complexity_evaluated_bases,
            "complexity_evaluated_fraction": self.complexity_evaluated_bases / callable_bases if callable_bases else None,
            "length_bp": self.length,
            "callable_bases": callable_bases,
            "ambiguous_bases": self.length - callable_bases,
            "ambiguous_base_fraction": (self.length - callable_bases) / self.length,
            "minimizer_occurrences": n,
            "unique_minimizer_occurrences": unique,
            "repetitive_minimizer_fraction": repeated / n if n else None,
            "high_copy_minimizer_fraction": (sum(count for copies, count in self.copies.items()
                                                  if copies >= high_copy) / n if n else None),
            "mean_seed_copy_number": (sum(copies * count for copies, count in self.copies.items())
                                      / n if n else None),
            "mean_log2_seed_copy_number": (sum(math.log2(copies) * count
                                               for copies, count in self.copies.items())
                                           / n if n else None),
            "max_seed_copy_number": max(self.copies, default=0),
            "unique_minimizers_per_kb": unique * 1000 / callable_bases if callable_bases else None,
            "seed_ambiguity": ambiguity,
            "anchor_gap_burden": extent,
            "unique_anchor_gap_bases": gap_bases,
            "longest_unique_anchor_gap_bp": max(self.gaps, default=0),
            "unique_anchor_gap_n50_bp": n50,
            "long_unique_anchor_gap_count": sum(count for size, count in self.gaps.items()
                                                  if size >= long_gap),
            "long_unique_anchor_gap_fraction": long_bases / callable_bases if callable_bases else None,
            "bases_in_seedless_acgt_runs": self.seedless_bases,
            "seedless_acgt_run_fraction": self.seedless_bases / callable_bases if callable_bases else None,
        }


def _measure_sequence(name, sequence, counts, k, w, length_scale, long_gap, intervals,
                      complexity_window, entropy_cutoff, complexity_tracts):
    stats = _Measurements()
    sequence = sequence.replace("U", "T")
    stats.length = len(sequence)
    seeds = iter(minimap2_minimizers(sequence, k, w))
    seed = next(seeds, None)
    for run in _DNA_RUN.finditer(sequence):
        start, end = run.span()
        stats.callable_bases += end - start
        if end - start >= _MIN_COMPLEXITY_RUN:
            stats.complexity_evaluated_bases += end - start
        for tract in _low_complexity_tracts(sequence, start, end, complexity_window, entropy_cutoff, length_scale):
            stats.low_complexity_bases += tract["length_bp"]
            stats.low_complexity_penalty_bases += tract["length_bp"] * tract["risk"]
            stats.low_complexity_peak_risk = max(stats.low_complexity_peak_risk, tract["risk"])
            stats.low_complexity_tract_count += 1
            stats.longest_low_complexity_tract = max(stats.longest_low_complexity_tract, tract["length_bp"])
            complexity_tracts.append({"contig": name, **tract})
        cursor = start
        run_seed_count = 0

        def add_gap(left, right):
            if right <= left:
                return
            size = right - left
            stats.gaps[size] += 1
            stats.gap_penalty_bases += size * (size / (size + length_scale))
            if size >= long_gap:
                intervals.append({"contig": name, "start": left, "end": right,
                                  "length_bp": size})

        while seed is not None and seed.start < end:
            if seed.start < start or seed.end > end:
                raise RuntimeError("minimizer crosses an ambiguous base")
            copies = counts[seed.key]
            stats.copies[copies] += 1
            run_seed_count += 1
            if copies == 1:
                add_gap(cursor, seed.start)
                cursor = max(cursor, seed.end)
            seed = next(seeds, None)
        add_gap(cursor, end)
        if run_seed_count == 0:
            stats.seedless_bases += end - start
    if seed is not None:
        raise RuntimeError("minimizer outside callable sequence")
    return stats


def _quantify_hardness(
    read_records, source, *, k=15, w=10, length_scale=1000, long_gap=1000,
    high_copy=50, seed_weight=0.5, local_weight=0.7, hotspot_weight=0.75,
    complexity_window=200, entropy_cutoff=4.0,
):
    """Return overall/per-contig hardness plus anchor gaps and complexity tracts.

    All FASTA records form ONE minimap2 reference: seed copy numbers include
    both within-record and between-record repeats, on either strand. Run each
    file separately when references will be indexed separately.

    R = mean(1 - 1/c) over minimizer occurrences, c = reference copy count.
    G = sum(g * g/(g+length_scale)) / callable_bases, where g are lengths
        without a globally unique minimizer's span, split at ambiguous bases.
    S = seed_weight * R + (1-seed_weight) * G (the original seed score / 100).
    Low-complexity windows have trinucleotide entropy E < entropy_cutoff and
    severity (entropy_cutoff-E)/entropy_cutoff. For each merged tract of size
    t, effective_length = mean window severity * t, and
    risk = effective_length/(effective_length+length_scale).
    C = (1-hotspot_weight)*base-weighted mean tract risk + hotspot_weight*max risk.
    H = 100 * ((1-local_weight)*S + local_weight*C).

    Defaults give 70% weight to local complexity, with 75% of that component
    from the hardest tract so unique flanks do not dilute every component.
    This emphasizes local coordinate consistency. Set local_weight=0 to
    recover the original score, also returned as seed_hardness_score.

    ``length_scale`` is the gap length or severity-weighted complexity tract
    length receiving half the maximum penalty;
    it is a comparison scale, not an inferred read length. ``long_gap`` only
    controls reported intervals and long-gap metrics. ``high_copy`` is a
    diagnostic cutoff, not minimap2's -f threshold. Scores without seeds are
    None. Ns are reported separately and break both gaps and complexity tracts.
    Complexity windows default to 200 bp; shorter DNA runs down to 32 bp use
    their whole length. Runs below 32 bp are unassessed for complexity.

    Uses two FASTA passes, O(sequence length) time for fixed w and memory for
    distinct seed keys plus one contig, a complexity window and reported intervals. Pure Python;
    intended for local reference FASTAs. Large genomes may require substantial
    time and RAM. Input can be plain or gzipped FASTA. No external dependencies.
    """
    _check_sketch_parameters(k, w)
    for name, value in (("length_scale", length_scale), ("long_gap", long_gap),
                        ("high_copy", high_copy), ("complexity_window", complexity_window)):
        if isinstance(value, bool) or not isinstance(value, int) or value <= 0:
            raise ValueError(f"{name} must be a positive integer")
    if complexity_window < _MIN_COMPLEXITY_RUN:
        raise ValueError(f"complexity_window must be at least {_MIN_COMPLEXITY_RUN}")
    for name, value in (("seed_weight", seed_weight), ("local_weight", local_weight), ("hotspot_weight", hotspot_weight)):
        if isinstance(value, bool) or not isinstance(value, (int, float)) or not 0 <= value <= 1:
            raise ValueError(f"{name} must be between 0 and 1")
    if isinstance(entropy_cutoff, bool) or not isinstance(entropy_cutoff, (int, float)) or not 0 < entropy_cutoff <= 6:
        raise ValueError("entropy_cutoff must be greater than 0 and at most 6 bits")
    counts = Counter()
    signatures = []
    for name, sequence in read_records():
        counts.update(seed.key for seed in minimap2_minimizers(sequence, k, w))
        signatures.append((name, hashlib.sha256(sequence.encode("ascii")).digest()))
    total = _Measurements()
    contigs = []
    intervals = []
    complexity_tracts = []
    for index, (name, sequence) in enumerate(read_records()):
        signature = (name, hashlib.sha256(sequence.encode("ascii")).digest())
        if index >= len(signatures) or signature != signatures[index]:
            raise ValueError("FASTA changed between analysis passes")
        stats = _measure_sequence(name, sequence, counts, k, w, length_scale, long_gap, intervals,
                                  complexity_window, entropy_cutoff, complexity_tracts)
        total.add(stats)
        contigs.append({"name": name, **stats.report(long_gap, high_copy, seed_weight, local_weight, hotspot_weight)})
    if len(contigs) != len(signatures):
        raise ValueError("FASTA changed between analysis passes")
    summary = total.report(long_gap, high_copy, seed_weight, local_weight, hotspot_weight)
    summary["distinct_minimizer_keys"] = len(counts)
    summary["contig_count"] = len(contigs)
    warnings = []
    if summary["ambiguous_bases"]:
        warnings.append("Ambiguous bases are excluded from gap scoring; inspect ambiguous_base_fraction.")
    if total.seedless_bases:
        warnings.append("Some ACGT runs have no seeds; inspect seedless_acgt_run_fraction.")
    if total.complexity_evaluated_bases < total.callable_bases:
        warnings.append("DNA runs below 32 bp are unassessed for complexity; inspect complexity_evaluated_fraction.")
    return {
        "schema_version": 2,
        "score_model": "local-complexity-v2",
        "reference": source,
        "parameters": {"k": k, "w": w, "homopolymer_compressed": False,
                       "length_scale": length_scale, "long_gap": long_gap,
                       "high_copy": high_copy, "seed_weight": seed_weight,
                       "local_weight": local_weight, "hotspot_weight": hotspot_weight,
                       "complexity_window": complexity_window, "entropy_cutoff": entropy_cutoff,
                       "complexity_word_length": 3, "complexity_min_run_bp": _MIN_COMPLEXITY_RUN},
        "summary": summary,
        "contigs": contigs,
        "long_unique_anchor_gaps": intervals,
        "low_complexity_tracts": complexity_tracts,
        "warnings": warnings,
    }


def quantify_reference_hardness(
    fasta_path, *, k=15, w=10, length_scale=1000, long_gap=1000,
    high_copy=50, seed_weight=0.5, local_weight=0.7, hotspot_weight=0.75,
    complexity_window=200, entropy_cutoff=4.0,
):
    """Score FASTA records together using the local-complexity hardness model."""
    path = Path(fasta_path)
    return _quantify_hardness(
        lambda: _read_fasta(path), str(path.resolve()), k=k, w=w,
        length_scale=length_scale, long_gap=long_gap, high_copy=high_copy,
        seed_weight=seed_weight, local_weight=local_weight,
        hotspot_weight=hotspot_weight, complexity_window=complexity_window,
        entropy_cutoff=entropy_cutoff,
    )


def quantify_sequence_hardness(sequence: str, *, name="query", **kwargs):
    """Score one in-memory sequence identically to a single-record FASTA.

    This avoids a temporary FASTA for every realignment window. No copies
    from other queries or the reference enter this sequence's seed counts.
    """
    sequence = sequence.upper()
    if not sequence:
        raise ValueError(f"empty sequence: {name}")
    bad = _INVALID_BASE.search(sequence)
    if bad:
        raise ValueError(f"invalid sequence base {bad[0]!r} in {name} at {bad.start()}")
    return _quantify_hardness(lambda: iter(((name, sequence),)), None, **kwargs)


def _overall_score(report):
    score = report["summary"]["hardness_score"]
    if score is None:
        raise ValueError(f"cannot score reference: {report['summary']['status']}")
    return score


def reference_hardness_score(fasta_path, **kwargs) -> float:
    """Return one overall 0–100 hardness float for all contigs in a FASTA.

    Keyword arguments match quantify_reference_hardness. Raise ValueError if
    the reference has no minimizer seeds and cannot be scored.
    """
    return _overall_score(quantify_reference_hardness(fasta_path, **kwargs))


def _same_file(left, right):
    return left.resolve() == right.resolve() or (left.exists() and right.exists() and left.samefile(right))


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("reference", type=Path, help="reference FASTA (plain or gzip)")
    parser.add_argument("-o", "--output", type=Path, help="output file; default: stdout")
    parser.add_argument("--format", choices=("score", "json"),
                        help="score: one overall float; json: detailed report "
                             "(default: score, or json for .json output)")
    parser.add_argument("--bed", type=Path, help="write long unique-anchor gaps as BED4")
    parser.add_argument("-k", type=int, default=15, help="minimizer k-mer size (default: 15)")
    parser.add_argument("-w", type=int, default=10, help="minimizer window size (default: 10)")
    parser.add_argument("--length-scale", type=int, default=1000,
                        help="half-penalty scale for gaps and severity-weighted complexity tract lengths (bp; default: 1000)")
    parser.add_argument("--long-gap", type=int, default=1000,
                        help="minimum reported gap length in bp (default: 1000)")
    parser.add_argument("--high-copy", type=int, default=50,
                        help="diagnostic high-copy cutoff (default: 50; not minimap2 -f)")
    parser.add_argument("--seed-weight", type=float, default=0.5,
                        help="seed ambiguity versus gap burden within the seed component (default: 0.5)")
    parser.add_argument("--local-weight", type=float, default=0.7,
                        help="weight of local complexity versus the original seed score (default: 0.7; 0 restores original)")
    parser.add_argument("--hotspot-weight", type=float, default=0.75,
                        help="hardest tract versus mean burden within local complexity (default: 0.75)")
    parser.add_argument("--complexity-window", type=int, default=200,
                        help="trinucleotide entropy window in bp (default: 200; minimum: 32)")
    parser.add_argument("--entropy-cutoff", type=float, default=4.0,
                        help="trinucleotide entropy below which windows are low complexity (bits; default: 4)")
    args = parser.parse_args(argv)
    try:
        outputs = [p for p in (args.output, args.bed) if p is not None]
        if (any(_same_file(args.reference, output) for output in outputs)
                or (len(outputs) == 2 and _same_file(*outputs))):
            raise ValueError("reference, report output and BED output must be different files")
        output_format = args.format or ("json" if args.output and args.output.suffix.lower() == ".json" else "score")
        report = quantify_reference_hardness(
            args.reference, k=args.k, w=args.w, length_scale=args.length_scale,
            long_gap=args.long_gap, high_copy=args.high_copy, seed_weight=args.seed_weight,
            local_weight=args.local_weight, hotspot_weight=args.hotspot_weight,
            complexity_window=args.complexity_window, entropy_cutoff=args.entropy_cutoff,
        )
        rendered = (str(_overall_score(report)) + "\n" if output_format == "score"
                    else json.dumps(report, indent=2, allow_nan=False) + "\n")
        if args.bed:
            with args.bed.open("w") as handle:
                for gap in report["long_unique_anchor_gaps"]:
                    handle.write(f"{gap['contig']}\t{gap['start']}\t{gap['end']}\tunique_anchor_gap\n")
        if args.output:
            args.output.write_text(rendered)
        else:
            print(rendered, end="")
    except (OSError, ValueError) as error:
        parser.error(str(error))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
