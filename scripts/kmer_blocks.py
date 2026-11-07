#!/usr/bin/env python3
"""Cut a sequence into repeat-unit copies ("k-mer blocks") with a 31-mer novelty score.

Only canonical 31-mers that occur fewer than --max-occurrence (255) times are
used. With --valid-from, that set is taken from another sequence (the
reference), and the input is cut using only those k-mers; otherwise the input's
own counts decide. Walking the used k-mers in order, a k-mer not yet seen in
the current copy adds 1 (the score is capped at --cap, 100) and a repeated one
subtracts 1; when the score drops below 0, a new copy starts right after the
position of the last maximum.

Records of the input are concatenated. Output records are named
LABEL_bI:START-END|>all (offsets on the concatenated input), the cells format
BlockMatchDP reads.

    python tools/kmer_blocks.py -i ref.fa -l r -o r.kcells.fa
    python tools/kmer_blocks.py -i hap.fa -l q --valid-from ref.fa -o q.kcells.fa
"""
import argparse
import array
import os

try:
    import numpy as np
except ImportError as _error:
    raise ImportError(
        "kmer_blocks.py requires numpy; install it with `python install.py` "
        "or `python -m pip install numpy`"
    ) from _error

K = 31
CODES = bytes("ACGTacgt".find(chr(i)) % 4 if chr(i) in "ACGTacgt" else 4 for i in range(256))


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


def kmer_arrays(seq):
    """(starts, canonical 31-mers) of every 31-mer without an ambiguous base."""
    codes = np.frombuffer(seq.encode().translate(CODES), dtype=np.uint8)
    m = len(codes) - K + 1
    if m <= 0:
        return np.zeros(0, dtype=np.int64), np.zeros(0, dtype=np.uint64)
    bad = np.concatenate(([0], np.cumsum(codes > 3, dtype=np.int64)))
    keep = bad[K:] == bad[:m]
    del bad
    c = codes & np.uint8(3)
    rc = np.uint8(3) - c
    f = np.zeros(m, dtype=np.uint64)
    r = np.zeros(m, dtype=np.uint64)
    tmp = np.empty(m, dtype=np.uint64)
    for j in range(K):
        np.left_shift(f, np.uint64(2), out=f)
        np.bitwise_or(f, c[j:j + m], out=f)
        tmp[:] = rc[j:j + m]
        np.left_shift(tmp, np.uint64(2 * j), out=tmp)
        np.bitwise_or(r, tmp, out=r)
    del c, rc, tmp
    np.minimum(f, r, out=f)
    return np.flatnonzero(keep), f[keep]


def valid_kmer_array(seq, max_occ=255):
    """valid_kmers as a sorted numpy array."""
    kmers, counts = np.unique(kmer_arrays(seq)[1], return_counts=True)
    return kmers[counts < max_occ]


def valid_kmers(seq, max_occ=255):
    return set(valid_kmer_array(seq, max_occ).tolist())


def used_kmers(seq, valid, max_occ):
    """(starts, ids, distinct) of the 31-mers segment() walks: equal k-mers share
    an id in range(distinct). valid: a set, a numpy array (valid_kmer_array) or None."""
    starts, kmers = kmer_arrays(seq)
    if valid is None:
        _, ids, counts = np.unique(kmers, return_inverse=True, return_counts=True)
        keep, distinct = counts[ids] < max_occ, len(counts)
    else:
        if not isinstance(valid, np.ndarray):
            valid = np.sort(np.fromiter(valid, dtype=np.uint64, count=len(valid)))
        # valid is sorted and unique: a k-mer's index in it is its id
        ids = np.searchsorted(valid, kmers)
        keep = (valid[np.minimum(ids, len(valid) - 1)] == kmers if len(valid)
                else np.zeros(len(kmers), dtype=bool))
        distinct = len(valid)
    del kmers
    return starts[keep], array.array("i", ids[keep].astype(np.int32).tobytes()), distinct


def segment(seq, valid=None, max_occ=255, cap=100):
    """Block boundaries [(start, end)] covering seq."""
    starts, km, distinct = used_kmers(seq, valid, max_occ)
    # stamp[k] == epoch: k-mer id k was seen in the copy being walked
    stamp = array.array("i", [-1]) * distinct
    bounds = [0]
    i, n, epoch = 0, len(km), -1
    while i < n:
        epoch += 1
        score, best, best_i, j = 0, 0, i, i
        while j < n:
            k = km[j]
            if stamp[k] == epoch:
                score -= 1
            else:
                stamp[k] = epoch
                score = min(cap, score + 1)
            if score >= best:
                best, best_i = score, j
            if score < 0:
                break
            j += 1
        if j >= n:
            break
        cut = best_i + 1 if best_i + 1 > i else j     # new copy starts after the last maximum
        if cut <= i:
            cut = i + 1
        bounds.append(int(starts[cut]) if cut < n else len(seq))
        i = cut
    bounds.append(len(seq))
    return [(a, b) for a, b in zip(bounds, bounds[1:]) if b > a]


def write_cells(path, label, seq, blocks):
    with open(path + ".tmp", "w") as out:
        for i, (a, b) in enumerate(blocks, 1):
            out.write(f">{label}_b{i}:{a}-{b}|>all\n{seq[a:b]}\n")
    os.replace(path + ".tmp", path)


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("-i", "--input", required=True)
    parser.add_argument("-l", "--label", required=True)
    parser.add_argument("-o", "--output", required=True)
    parser.add_argument("--valid-from", help="take the usable k-mers from this FASTA")
    parser.add_argument("--max-occurrence", type=int, default=255)
    parser.add_argument("--cap", type=int, default=100)
    args = parser.parse_args(argv)
    seq = "".join(s for _, s in records(args.input))
    valid = None
    if args.valid_from:
        valid = valid_kmers("".join(s for _, s in records(args.valid_from)), args.max_occurrence)
    blocks = segment(seq, valid, args.max_occurrence, args.cap)
    write_cells(args.output, args.label, seq, blocks)
    print(f"[kmer_blocks] {len(blocks)} blocks -> {args.output}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
