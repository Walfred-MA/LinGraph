# Trio benchmarks

These scripts measure how often a child's variant calls are also found in at
least one parent. Run the examples from the LinGraph repository root.

| Script | Method | Usage guide |
| --- | --- | --- |
| [QuickTriocheck.py](QuickTriocheck.py) | Match variants by type, position, and size; supports separate haplotype VCFs or one multi-sample VCF. Uses Python's standard library. | [QuickTriocheck README](QuickTriocheck.README.md) |
| [truvari_trio.sh](truvari_trio.sh) | Run eight Truvari comparisons between two child and four parental haplotypes. | [Truvari trio README](truvari_trio.README.md) |

## Inputs and interpretation

Use child and parental calls made against the **same reference**, with matching
contig names. Separate-file benchmarks require the original LinGraph
`##referenceCoverage` headers for parental coverage filtering. Keep those
headers: `tools/grvcf_to_vcf.py` currently removes them from standard VCF output.

QuickTriocheck trims **10,000 bp from each end of every coverage interval** by
default in separate-file mode; the Truvari wrapper does not trim these edges.
Their matching rules and default size filters also differ. Record the options
and Truvari version when comparing results.

A child-only call has no parental match under the selected rules. Although the
scripts label these calls `FP` or report an `error_rate`, they are not confirmed
false positives: true de novo variation and missed parental calls can also
produce this result. The reported consistency measures parental support, not
full Mendelian genotype or transmission consistency.
