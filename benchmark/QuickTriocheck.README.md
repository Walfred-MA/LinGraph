# QuickTriocheck.py

Check whether child variants have a match in either parent. The script supports
standard `GT`, LinGraph's separate FORMAT fields, and older packed `HSV` fields.
It matches SVs by coordinates, type, and size rather than sequence alignment.

## Requirements

- Python 3.9 or later; no third-party Python packages.
- Plain `.vcf` or gzip-compressed `.vcf.gz` files against the same reference.
- For separate haplotype VCFs, retain `##referenceCoverage` headers. With the
  default edge trimming, these are required in **every child and parent VCF**.

No FASTA, VCF index, or Truvari installation is needed.

## Separate haplotype VCFs

Run from the repository root. Supply one or more VCFs per role; a diploid trio
normally has two haplotype VCFs per person:

```bash
python3 benchmark/QuickTriocheck.py \
  --child child_h1.vcf child_h2.vcf \
  --mother mother_h1.vcf mother_h2.vcf \
  --father father_h1.vcf father_h2.vcf \
  --fp child.fp.vcf --tp child.tp.vcf \
  > trio.quick.tsv 2> trio.quick.log
```

For LinGraph output, replace each filename with
`cohort_calls/samples/NAME/NAME.vcf`.

Coverage intervals are zero-based and half-open. For example:

```text
##referenceCoverage=<Sample="mother_h1",Chrom="chr1",Start=0,End=100000>
```

The script trims 10,000 bp from both ends of each raw interval before merging
coverage. In this example, only `[10000, 90000)` remains. Calls in the trimmed
edges are excluded. A child call is also excluded if **any supplied parental
haplotype** lacks coverage at its reference anchor. These exclusions are
reported separately and are removed from the consistency denominator.

Use `--trim-edges 0` to disable edge trimming. Parental coverage headers are
still required; child coverage headers are then optional.

## One multi-sample VCF

With `-i`, `--child`, `--mother`, and `--father` select VCF sample columns:

```bash
python3 benchmark/QuickTriocheck.py \
  -i cohort_calls/cohort.sv.vcf \
  --child HG002 --mother HG004 --father HG003 --prefix-match \
  --fp child.fp.vcf --tp child.tp.vcf \
  > trio.quick.tsv 2> trio.quick.log
```

`--prefix-match` selects haplotype columns such as `HG002_h1` and `HG002_h2`.
Without it, each selector must match a column name exactly. The three roles
must not overlap. This mode excludes child calls with missing genotypes in any
selected parental column; it does not apply coverage-edge trimming.

## Matching options

| Option | Behavior |
| --- | --- |
| `--min-size 50` | Default SV-only mode includes sizes **strictly greater than** 50 bp. |
| `--slop 500` | Default maximum difference for both SV start and end positions. |
| `--size-ratio 0.7` | Require `min(child_size, parent_size) / max(child_size, parent_size) >= 0.7`. |
| `--all` | Include SNPs and small variants; SNP/INDEL matching is exact. |
| `--trim-edges BP` | Change the coverage-edge trim in separate-file mode; default 10,000 bp. |
| `--original-distance BP` | Match SVs using original per-haplotype reference intervals reconstructed from LinGraph annotations. `0` requires overlap; a positive value allows that many bases between intervals. Replaces `--slop`. |
| `--original-distance exact` | Require identical reconstructed reference start and end coordinates. |
| `--VarINS` | Include child alleles labeled as insertion-internal (`INS_`) events, which are excluded by default. |

Original-coordinate matching requires the relevant LinGraph per-observation
annotations; use the original grVCF. Records without recoverable coordinates
are excluded and counted in the report.

## Outputs

- **Standard output:** tab-separated settings, coverage exclusions, parent-match
  counts, fractions, and a breakdown by variant type. Redirect it to a file.
- **Standard error:** script version and input-loading information.
- **`--fp FILE`:** child calls without a match in either parent, annotated with
  `TRIOCHECK_FP` fields.
- **`--tp FILE`:** child calls with a match in at least one parent, annotated with
  parental-match information.

In separate-file mode, outputs are split by child source filename. For the
example above, the FP files are `child.fp.child_h1.vcf` and
`child.fp.child_h2.vcf`; TP files follow the same naming rule. There is no
combined FP/TP file in this mode. Each output preserves its source VCF header
and rows; annotations identify the relevant ALT alleles in multiallelic rows.
A `.gz` output uses ordinary gzip and is not automatically indexed. Existing
output paths are overwritten.

The main fraction is `child_found_in_either_parent_fraction`; its complement
is `child_not_found_in_parents_fraction`. Both use the retained child calls
after filtering. See [interpretation notes](README.md#inputs-and-interpretation).

```bash
python3 benchmark/QuickTriocheck.py --help
```
