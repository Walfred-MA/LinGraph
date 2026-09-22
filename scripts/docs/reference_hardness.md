# Reference alignment hardness for default minimap2

`reference_hardness.py` estimates how difficult it is to anchor alignments
consistently in a reference FASTA. Version 2 emphasizes local low-complexity
severity and tract length, giving the hardest tract substantial weight. It
retains **minimap2 default seeding: k=15, w=10, ordinary minimizers, no assembly
preset**. The command prints one overall float; optional JSON provides its
components, per-contig measurements, anchor gaps, and low-complexity tracts.
The scorer itself is read-only. `graphcigartoref.py` uses it to gate optional
whole-pair tandem realignment, as described below.

**The default numeric scale has changed.** Rescore all comparison files with
the same settings. `--local-weight 0` reproduces the original seed-only score;
JSON also retains that value as `seed_hardness_score`. The JSON schema version
is now 2 and its `score_model` is `local-complexity-v2`.

## Python function

```python
from reference_hardness import quantify_reference_hardness

report = quantify_reference_hardness("reference.fa")
summary = report["summary"]
print(summary["hardness_score"])                 # 0–100; higher is harder
print(summary["repetitive_minimizer_fraction"])
print(summary["longest_unique_anchor_gap_bp"])
print(summary["long_unique_anchor_gap_fraction"])
```

The implementation uses only the Python standard library. Plain and gzip
FASTA are accepted. All records in a FASTA are treated as **one indexed
reference**, so copy counts include repeats within a contig, between contigs,
and on the reverse-complement strand. Analyze files separately when minimap2
will index them separately. Lowercase bases are not masked. Ns and other IUPAC
ambiguity codes break seeds; U is treated as T, as in minimap2.

## Realignment policy

`graphcigartoref.py` calls `quantify_sequence_hardness(sequence)` with the
same defaults and implementation as the FASTA function for the early whole-pair
swap, avoiding a temporary FASTA. Hardness uses query sequence alone, independently of the
reference, other samples, and soft-mask case. The default model combines **70% local complexity / 30% seed-anchor
hardness** (`local_weight=0.7`).

- Early whole-pair swap retains the repeated-traversal tandem trigger. It
  runs default minimap2/mappy only when **overall query hardness <10**.
  This restores the policy before the sliced-region change: score the whole
  query provided to the early pairwise stage, including comparison flanks,
  without applying an owned/core subinterval to hardness. It does not score
  the entire assembly chromosome. Exactly 10 and unscorable queries are rejected.
  It accepts only a strict improvement under `4*matches - 8*mismatches`, with
  gap cost `24 + 8*length` per contiguous I/D run.
- SV-window realignment does **not** use a hardness gate. Its existing window
  selection and size limits still apply. Replacement requires a strict score improvement, with
  match +4, mismatch -8, and gap cost
  `24 + 8*length + 160*log2(masked_bases + 1)`. Masked bases are lowercase
  bases **inside that gap**: query bases for I, reference bases for D.
  The +1 defines the zero-masked case; there is no neighborhood denominator.
  Sequence offsets and adjacent gap runs are included in replacement scoring.
- Adjacent D/I sequence correction remains independent of the hardness gate
  and keeps its existing sequence-support checks and affine score:
  match +1, mismatch -4, gap `4 + (length - 1)`. Its existing short-junction
  safeguard is unchanged. Composite windows with aligned sequence between
  gaps follow the SV-window policy and its masked score, rather than the
  adjacent D/I affine score.

Both replacement policies resolve matches/mismatches from actual sequences
and preserve complete reference/query spans. Equal-score replacements are
rejected. Existing equivalent-gap normalization, trimming and insertion
encoding are separate from these optional alignment decisions. An uppercase
FASTA still receives complexity-based hardness, but supplies no masked bases
for the SV gap penalty. These thresholds are policy choices, not validated
guarantees of improved trio concordance.

## Overall hardness float on stdout

```bash
python3 reference_hardness.py reference.fa
```

The default command prints exactly one floating-point number followed by a
newline, with no header or explanatory text. Higher values (0–100) mean harder.
All contigs in the FASTA contribute to this one score, using pooled seed
occurrences, callable bases, and the hardest low-complexity tract across all
contigs. Cross-contig seed copies count as repeats. It is not an arithmetic
mean of the contig scores. Default minimap2 seeding remains k=15,
w=10, with no assembly preset.

The equivalent Python function returns a float:

```python
from reference_hardness import reference_hardness_score

print(reference_hardness_score("reference.fa"))
```

If the FASTA has no usable minimizers, the CLI exits nonzero with an explanation
on stderr and no value on stdout; the function raises `ValueError`.

## Detailed JSON and BED reports

```bash
python3 reference_hardness.py reference.fa \
  -o hardness.json --bed long_anchor_gaps.bed
```

A `.json` output extension selects detailed JSON automatically. To print
detailed JSON on stdout, use `--format json`. The default output is a scalar
score. BED coordinates are zero-based,
half-open. `--long-gap 1000` controls which gaps are included in BED and the
long-gap statistics; smaller gaps still contribute to the score. JSON includes
the same intervals under `long_unique_anchor_gaps`. Low-complexity tract
coordinates and severities are separately available under
`low_complexity_tracts` in JSON; the existing BED option still reports anchor
gaps only. All coordinates are relative to each FASTA record.

Use `-k` and `-w` only if the mapping invocation overrides default seeding.
HPC (`-H`) is not supported. If using a prebuilt minimap2 index, match its
actual seeding settings. The reference copy counts are computed before any
minimap2 frequency filtering or seed rescue.

## Score definition

The original seed-based component and a new local-complexity component both
lie between zero and one. All weights are configurable heuristic choices;
no alignment-error probabilities or easy/hard thresholds are implied.

**Seed ambiguity R:** let `c_i` be the number of indexed occurrences of the
canonical minimizer at occurrence `i`, and `M` the total occurrence count:

```text
R = sum_i (1 - 1/c_i) / M
```

Unique seeds contribute 0, two-copy seeds 0.5, and highly repeated seeds
approach 1. This is an occurrence-weighted measure: abundant repetitive seeds
receive their full weight. Also inspect `repetitive_minimizer_fraction`
(`c_i > 1`), `mean_seed_copy_number`, `mean_log2_seed_copy_number`, and
`max_seed_copy_number`. The unbounded copy-number metrics help distinguish
very repetitive references whose bounded scores are both near 100.

**Anchor-gap burden G:** a unique anchor is the full k-mer interval of a seed
that occurs exactly once in the entire reference. Merge overlapping unique
anchor spans, then measure their complementary gaps within each uninterrupted
ACGT/U run. Gaps at contig ends are included, and gaps never cross ambiguous
bases or contig boundaries. Let `g_j` be these lengths and `B` the number of
ACGT/U bases:

```text
G = sum_j g_j * g_j / (g_j + L) / B
```

`L = length_scale` (default 1,000 bp) is the gap length receiving half the
maximum per-base penalty. A 100 bp gap contributes about 9 penalty bases with
this default; a 10 kb gap contributes about 9,091. The penalty keeps increasing
smoothly beyond L. Thus, one long gap contributes more than the
same number of bases split among short gaps. This captures the idea that long
repetitive stretches are more difficult to anchor consistently.

`length_scale` is a tunable comparison scale, **not an inferred query length**.
Use the same scale across references. For a different relevant span, for
example 10 kb, pass `length_scale=10000` or `--length-scale 10000`.

The retained seed component is:

```text
S = alpha * R + (1 - alpha) * G
alpha = seed_weight = 0.5
seed_hardness_score = 100 * S
```

### Local low-complexity severity and tract length

Scan every 200 bp window within uninterrupted ACGT/U runs, using all 198
overlapping trinucleotides. Their empirical Shannon entropy is:

```text
E = -sum(p_word * log2(p_word))
window severity = max(0, (E_cutoff - E) / E_cutoff)
E_cutoff = entropy_cutoff = 4 bits
```

There are 64 possible trinucleotides, with maximum entropy 6 bits. A window
below the cutoff has positive severity; homopolymers have severity 1.
Interrupted C/T-rich repeats can also have substantial severity even if their
individual 15-mers differ. The cutoff is a diagnostic choice, not a calibrated
alignment failure threshold.

Merge the full spans of overlapping or touching positive-severity windows
into tracts. Each tract includes the windows' boundary context; its coordinates
are not exact biological repeat boundaries. For tract j:

```text
s_j = mean severity of its qualifying windows
t_j = tract length in bp
e_j = s_j * t_j                       # severity-weighted length
q_j = e_j / (e_j + L)                 # tract risk
L = length_scale = 1000 bp
```

Both lower complexity and longer tracts increase risk. Length is scaled by
severity before applying the saturating transform: a long moderately simple
tract can be difficult, even if it is not a perfect homopolymer. L corresponds
to risk 0.5 at a severity-weighted length of 1,000 bp.

DNA runs of 32–199 bp use their whole length as one window. Runs below 32 bp
are unassessed for complexity and reported in `complexity_evaluated_fraction`.
No window, tract, or anchor gap crosses an ambiguous base or contig boundary.
Lowercase is treated like uppercase, and U like T.

Combine a whole-reference burden with the hardest tract:

```text
C_mean = sum(t_j * q_j) / B            # B = total ACGT/U bases
C_peak = max(q_j), or 0 if no tracts
C = (1 - beta) * C_mean + beta * C_peak
beta = hotspot_weight = 0.75
```

The maximum is taken across all contigs, without averaging contig maxima.
Adding easy flanks lowers the mean burden but leaves the hardest-tract
component intact when that tract is unchanged.

### Final scalar

```text
H = 100 * ((1 - lambda) * S + lambda * C)
lambda = local_weight = 0.7
```

Thus default weights are 15% seed ambiguity, 15% anchor-gap burden, 17.5% mean
local-complexity burden, and **52.5% hardest-tract risk**. This prioritizes local
coordinate difficulty. Whole duplicated but compositionally varied sequences
receive less weight than before; use `seed_hardness_score` when global mapping
placement ambiguity is the main concern. A long homopolymer approaches 100,
and ordinary random sequence is near zero.

To recover the exact previous scalar:

```bash
python3 reference_hardness.py reference.fa --local-weight 0
```

Use `--local-weight`, `--hotspot-weight`, `--complexity-window`,
`--entropy-cutoff`, and `--length-scale` to change the tradeoffs. Compare files
using identical settings. The default weights were chosen to emphasize local
repeat difficulty, not fitted or calibrated against a genome-wide benchmark.

Useful measurements returned alongside H:

| Field | Meaning |
| --- | --- |
| `seed_hardness_score` | Original seed-only 0–100 score for backward comparison |
| `local_complexity_risk` | Weighted mean and hardest-tract risk, between 0 and 1 |
| `low_complexity_mean_risk` | Base-weighted risk burden over all callable bases |
| `low_complexity_peak_risk` | Largest risk of any low-complexity tract |
| `low_complexity_base_fraction` | Fraction of callable bases covered by low-complexity tracts |
| `longest_low_complexity_tract_bp` | Longest merged low-complexity span |
| `complexity_evaluated_fraction` | Callable-base fraction in runs at least 32 bp long |
| `repetitive_minimizer_fraction` | Fraction of seed occurrences with more than one reference copy |
| `unique_minimizers_per_kb` | Unique seed density per 1,000 ACGT/U bases |
| `longest_unique_anchor_gap_bp` | Largest interval without coverage by a globally unique seed |
| `unique_anchor_gap_n50_bp` | Gap length at which descending gaps cover half of all gap bases |
| `long_unique_anchor_gap_fraction` | Fraction of ACGT/U bases in gaps at least `long_gap` bp |
| `high_copy_minimizer_fraction` | Fraction of seed occurrences with at least `high_copy` copies (default 50) |
| `ambiguous_base_fraction` | Fraction excluded because the bases are not ACGT/U |
| `seedless_acgt_run_fraction` | Fraction of callable bases in uninterrupted runs producing no seeds |

`high_copy` is a descriptive cutoff, **not** minimap2's `-f` threshold.
References/contigs without seeds return `hardness_score: null` and a status of
`no_seeds` or `no_callable_bases`, rather than appearing easy. Inspect the
missing-base and seedless-run fractions even when a pooled score is available.

## Interpretation and limitations

- **These are unique-anchor gaps, not exact repeat boundaries or repeat-unit
  sizes.** They approximate the extent of poorly anchored regions. Seedless
  sequence can also create a gap. Short tandem motifs repeated for many kb
  produce a long gap; the motif length itself is not measured.
- Unique exact seeds may distinguish diverged repeat copies, but errors or
  variants in the query can remove those seeds. Conversely, a combination of
  individually repeated seeds can form a unique chain. This score models
  neither effect, so gap lengths can overestimate or underestimate alignment
  difficulty. Long near-identical duplications with scattered distinguishing
  seeds can be underestimated.
- Query length, divergence, sequencing error, chaining, frequency filtering,
  seed rescue, and base-level alignment all influence actual coordinates and
  runtime. Base-level indel placement can remain ambiguous even with strong
  flanking anchors. Reference-only scoring cannot establish coordinate
  consistency; validate it with representative queries before using cutoffs.
- The hardest-tract term intentionally lets a small difficult locus influence
  a large file. This is not the expected error rate of a random base or query.
  Inspect the tract locations and per-contig diagnostics when interpreting it.
- Trinucleotide entropy does not identify every repeat. Long, varied repeat
  units may remain above the cutoff; compositionally biased sequence may be
  flagged even without a tandem-repeat structure. The seed component retains
  some sensitivity to those other repeats. Complexity alone cannot prove that
  an alignment or SV coordinate is correct.
- Ambiguous bases are excluded from G, not labeled as repeats. A mostly-N
  reference can therefore have a low score on its few remaining callable
  bases; `ambiguous_base_fraction` is essential context.
- This pure-Python implementation makes two FASTA passes and retains distinct
  seed counts plus one contig, a complexity window and reported intervals. It is intended for local
  reference FASTAs. Whole-genome inputs may take substantial time and memory.

## Validation and provenance

```bash
python3 -m unittest discover -s tests -p test_reference_hardness.py -v
```

The tests include random sequence, tandem repeats, long versus fragmented
duplications, interrupted C/T-rich sequence, complexity severity and length,
hard-tract persistence after adding easy flanks, independent direct entropy
calculations, reverse complements, Ns, short sequences, malformed FASTA,
gzip input, scalar stdout, JSON/BED output, input protection, and pooled statistics. When minimap2 is installed,
an independent test builds native `.mmi` indexes and compares **every seed
key, coordinate, strand, and multiplicity** against this implementation for
six k/w combinations, including defaults, even k, ties, and terminal partial
windows. This was verified with minimap2 2.31-r1302. The index decoder is only
used by the test; production scoring does not parse `.mmi` files.

The ordinary minimizer algorithm is adapted from Heng Li's minimap2
[`sketch.c`](https://github.com/lh3/minimap2/blob/v2.31/sketch.c).
Default seeding settings and frequency filtering are described in the
[minimap2 manual](https://lh3.github.io/minimap2/minimap2.html) and
[`options.c`](https://github.com/lh3/minimap2/blob/v2.31/options.c).

### Minimap2 license for the adapted sketch implementation

The MIT License

Copyright (c) 2018-     Dana-Farber Cancer Institute
              2017-2018 Broad Institute, Inc.

Permission is hereby granted, free of charge, to any person obtaining
a copy of this software and associated documentation files (the
"Software"), to deal in the Software without restriction, including
without limitation the rights to use, copy, modify, merge, publish,
distribute, sublicense, and/or sell copies of the Software, and to
permit persons to whom the Software is furnished to do so, subject to
the following conditions:

The above copyright notice and this permission notice shall be
included in all copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,
EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF
MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND
NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS
BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN
ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN
CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
SOFTWARE.
