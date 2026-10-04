"""Which zero-length nested rows on an allele's query edge belong to it.

A nested deletion has a zero-length QUERYCOORD (a point). When one sample
carries several alleles of one path back to back (a repeated row, copies of
one duplication template), a point shared by two alleles cannot say which
one it belongs to. The merged VCF is read with this rule, the same in
gfa_junctions.py, gfa_sample_gaf.py and tools/check_merge_lossless.py:

* the allele's other nested rows are applied first;
* then deletions at its first query point, from path offset 0 on, and
  deletions at its last query point, from the path end back, each only if
  it continues that chain, overlaps no deletion already applied, and does
  not make the allele shorter than its query span.

Standard library only.
"""


def edge_owned(starts, ends, kept_spans, current, span, length):
    """Keys of the edge deletions this allele owns.

    ``starts`` / ``ends``: ``[(first, last, key)]`` deletions of path
    ``[first, last)`` at the allele's first / last query point. ``kept_spans``:
    ``[(first, last)]`` path intervals its other nested rows replace.
    ``current``: its length with those rows applied; ``span``: its query
    length; ``length``: the path length.
    """
    owned, taken = [], list(kept_spans)

    def fits(first, last):
        return (current - (last - first) >= span
                and all(last <= a or b <= first for a, b in taken if a < b))

    covered = 0
    for first, last, key in sorted(starts, key=lambda value: (value[0], value[1])):
        if first == covered and fits(first, last):
            owned.append(key)
            taken.append((first, last))
            current -= last - first
            covered = last
    covered = length
    for first, last, key in sorted(ends, key=lambda value: (-value[1], -value[0])):
        if last == covered and fits(first, last):
            owned.append(key)
            taken.append((first, last))
            current -= last - first
            covered = first
    return owned
