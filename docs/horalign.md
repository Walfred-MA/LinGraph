# horalign: alignment across a tandem repeat array

`horalign` aligns one query sequence to one reference sequence through a large
tandem repeat, such as a centromeric HOR array, and writes one graph CIGAR.

```bash
python scripts/LinGraph.py horalign -r reference.fa -q query.fa -o result.gcigar -t 16
```

Both FASTAs must hold exactly one sequence, on the same strand; the ends of
both are treated as anchored. BlockMatchDP must be built
(`python scripts/install.py`) or given with `--blockmatchdp`.

## Steps and threads

| Step | Threads |
| --- | --- |
| k-mer blocks (`kmer_blocks.py`): repeat-unit copies from 31-mers seen < 255 times in the reference; reference and query cut at the same time | 2 processes |
| BlockMatchDP k-mer matching (`--pair-score log --gap-factor 0.1`) | `-t` |
| Per-block base alignment (each paired copy with ±15 kb flanks, minimap2) | `-t` worker processes |
| Gap relocation (every query gap mapped to the whole reference) | minimap2 `-t` |
| Merge, trim, graph CIGAR | 1 |

Thread count does not change the result.

## Output

One line: the main path `>REF:OPS` walks the reference end to end (`=`, `X`
with reference then query bases, `I` with the inserted query bases, `D`). A
relocated piece is its own segment at its query position, `>REF:STARTH…ENDH`
(forward) or `<REF:ENDH…STARTH` (reverse strand; H counts reference outside the
piece in its own direction); the next `>OPS` resumes the main path. Relocations
are query gaps whose clear best minimap2 hit lies elsewhere, plus aligned
copies that would run the main path backwards (e.g. extra tandem copies).

| Option | Effect |
| --- | --- |
| `--keep-all` | keep intermediate files in `OUTPUT.horalign/` (always kept after a failure) |
| `--no-relocation` | do not map gaps; every query gap is inserted sequence |

## Figure

```bash
python scripts/horalign_plot.py -g child.gcigar -g mother.gcigar \
  --label "HG002_h2 h2tg000007l" --label "HG004_h1 (mother) h1tg000007l" \
  --ref-label "CHM13 chr6" --ref-start 58270000 \
  --query-start 58460000 --query-start 110590000 --query-strand + --query-strand - \
  -o HG002_h2.horalign
```

Reference on top, each query below (the first is the focal one): aligned
blocks, insertions (relocations included) and deletions of at least
`--min-size` (50) bp, as pentagons proportional to size. Rows after the first
report events shared exactly with the first query. Writes PDF (editable text)
and PNG; needs matplotlib.
