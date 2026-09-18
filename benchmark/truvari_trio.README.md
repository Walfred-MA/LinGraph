# truvari_trio.sh

Benchmark a diploid trio with Truvari while keeping all six haplotype callsets
separate. The script compares each child haplotype against each of the four
parental haplotypes, for **eight pairwise benchmarks**.

## Requirements

- Bash and Unix command-line utilities, including `awk`, `sort`, `grep`, `sed`,
  `comm`, `gzip`, and `mktemp`.
- `truvari` with the `bench` and `consistency` subcommands.
- HTSlib's `bgzip` and `tabix` commands.

All commands must be on `PATH`. The wrapper does not call `bcftools` or
`bedtools`. Record the installed Truvari version with `truvari --version`;
matching and filtering use that version's defaults except for `--refdist`.

## Inputs

Supply exactly two child, two maternal, and two paternal VCFs. Each input is
one haplotype, uses the same reference and contig names, and can be plain `.vcf`
or gzip-compressed `.vcf.gz`.

Every parental VCF must contain nonempty `##referenceCoverage` intervals for
exactly one sample, for example:

```text
##referenceCoverageCoordinateSystem=0-based-half-open
##referenceCoverage=<Sample="mother_h1",Chrom="chr1",Start=0,End=100000>
```

The script intersects coverage across **all four parental haplotypes**, then
passes that intersection to every `truvari bench` run as `--includebed`.
Missing coverage or an empty intersection stops the benchmark. This wrapper
uses the intervals without edge trimming. Child coverage headers are not
required for this filtering.

Inputs must contain alleles supported by the installed Truvari version. The
wrapper adds missing contig headers, sorts, compresses, and indexes temporary
copies; it does not decode graph-encoded alleles or normalize REF/ALT sequences.

## Run

From the repository root:

```bash
bash benchmark/truvari_trio.sh \
  --child child_h1.vcf child_h2.vcf \
  --mother mother_h1.vcf mother_h2.vcf \
  --father father_h1.vcf father_h2.vcf \
  --distance 500 --jobs 8 \
  --fp child.fp.vcf \
  > trio.truvari.tsv 2> trio.truvari.log
```

For LinGraph output, replace each filename with
`cohort_calls/samples/NAME/NAME.vcf`, retaining the required coverage metadata.

| Option | Default | Purpose |
| --- | --- | --- |
| `--distance BP` | `500` | Pass breakpoint distance to Truvari as `--refdist`; `0` removes positional tolerance while retaining other Truvari rules. |
| `--jobs N` | `8` | Maximum number of pairwise benchmarks running concurrently. |
| `--fp FILE` | No VCF output | Write child calls without a match in any parental haplotype. |
| `--keep-temp` | Off on success | Retain prepared VCFs, pairwise results, the coverage BED, and logs. |

## Results

The child VCF is the **base** in every pairwise benchmark; the parent is the
comparison callset. The script uses the child-side `fn.vcf.gz` outputs to find
child calls absent from both maternal, both paternal, or all four parental
haplotypes. It reports each child haplotype separately and also sums their
counts, retaining haplotype copies as separate observations.

| Output | Contents |
| --- | --- |
| `trio.truvari.tsv` | Child totals, parental matches, `trio_consistency`, and `error_rate`, captured from standard output. |
| `trio.truvari.log` | Progress and diagnostic messages, captured from standard error. |
| `child.fp.h1.vcf` | Unmatched child h1 calls, preserving the original child header and sample columns. |
| `child.fp.h2.vcf` | Unmatched child h2 calls, preserving the original child header and sample columns. |
| `child.fp.vcf` | Combined sites-only output with `INFO/TRIO_HAP=h1` or `h2`. |

The last three files are written only with `--fp`. A `.vcf.gz` FP path produces
BGZF-compressed, tabix-indexed outputs. Existing output paths are overwritten.

`trio_consistency` is the fraction of benchmarked child calls found in at least
one parent. `error_rate` is the fraction absent from all four parental
haplotypes. The denominator is Truvari's retained child calls, calculated from
`tp-base.vcf.gz` plus `fn.vcf.gz`, rather than every input row. Zero denominators
are reported as `NA`. See [interpretation notes](README.md#inputs-and-interpretation).

Temporary data are created under `truvari_trio_tmp.*` in the current directory.
They are removed after a successful run unless `--keep-temp` is set, and retained
on failure for inspection. The shared coverage BED is
`parents.all4.covered.bed` inside that directory.
