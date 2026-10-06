# LinGraph utilities

Trio benchmarking scripts are in [`benchmark/`](../benchmark/README.md), with
separate guides for QuickTriocheck and the Truvari wrapper.

## Prepare assembly FASTAs

Use `prepare_assemblies.py` for soft masking, optional contig-name fixing,
and indexing. It writes prepared copies in a separate output directory:

```bash
python tools/prepare_assemblies.py \
  -i assembly.fa --name HG002_h1 \
  -O prepared_assemblies --contignamefix
```

For name fixing alone, the `-n` value can use either assembly-list notation
(`HG002_h1`) or FASTA-prefix notation (`HG002#1`). The output defaults to a
new sibling FASTA, so the source is preserved:

```bash
python tools/namecontigsfix.py \
  -i assembly.fa -n HG002_h1
samtools faidx assembly.namefixed.fa
```

Use `assembly.namefixed.fa` in the assembly list. Pass `-o output.fa` to
choose another output path. The resulting contig IDs begin with
`HG002#1#` (for example, `HG002#1#chr1`).

Keep native contig names for `CHM13_h1` and `HG38_h1`; omit name fixing for
those references. LinGraph and the standalone sample caller check inputs
without modifying assemblies or indexes and stop with a preparation command
when the checks fail. See [assembly preparation](assembly_preparation.md)
for batch preparation, masking options, and dependencies.

## Convert grVCF to standard VCF

`tools/grvcf_to_vcf.py` exports a plain VCF in one streaming pass, without
reading a reference FASTA. Convert one file at a time:

```bash
python3 tools/grvcf_to_vcf.py \
  -i cohort_calls/cohort.sv.vcf \
  -o cohort_calls/cohort.sv.standard.vcf
```

Repeat for the indel and SNP files. Plain and gzip-compressed inputs are
accepted. The converter:

- Removes nested records on `INS_`, `DEL_`, `SUB_`, or `DUP_` parents, and
  records at `POS=0`, which lack a VCF anchor base.
- Converts `<INS>` with plain `INFO/SEQ` into explicit bases (`ALT=REF+SEQ`).
  Other symbolic alleles, including insertions with graph-encoded sequence,
  remain symbolic; it does not decode them during this export.
- Preserves record order, IDs, REF, QUAL, FILTER, sample names, and GT values.
  INFO retains only `SVTYPE`, `END`, `SVLEN`, and `NSUP`; FORMAT retains GT.
- Does not validate against a reference, normalize or sort alleles, calculate
  allele frequencies, or produce an unplaced-record sidecar. Non-nested
  alternative-locus coordinates remain in their original coordinate system.

This is a representative-allele export, **not a lossless assembly export**.
Keep the original grVCF and its dependencies for nested variation and graph
construction. Check the reported `written`, `nested_removed`, and
`position_zero_removed` counts before using the output downstream.

Use `.vcf.gz` output for BGZF compression (requires `bgzip`); the converter
does not create an index. `--force` replaces an existing output. The legacy
`-r` option is accepted but ignored; `-a` and `--normalize` are not supported.
Run `python3 tools/grvcf_to_vcf.py --help` for supported options.

## Rename reference chromosomes in a VCF

`tools/vcf_chromfix.py` renames reference contigs in LinGraph grVCFs and
plain VCFs. Choose exactly one mode:

```bash
# CHM13 RefSeq accessions NC_060925.1..NC_060948.1 -> chr1..chr22, chrX, chrY
python3 tools/vcf_chromfix.py -i calls.vcf.gz -o calls.chr.vcf.gz --CHM13fix
# Drop everything up to the last '#': HG19#1#chr1 -> chr1
python3 tools/vcf_chromfix.py -i calls.vcf -o calls.chr.vcf --noprefix
# Drop only this leading text; other names stay unchanged
python3 tools/vcf_chromfix.py -i calls.vcf -o calls.chr.vcf --noprefix 'HG19#1#'
# Two-column table (tab or whitespace): column 1 -> column 2
python3 tools/vcf_chromfix.py -i calls.vcf -o calls.chr.vcf --fixtable names.tsv
```

Plain and gzip-compressed inputs are accepted. An output ending in `.gz` is
gzip, not BGZF; run `bgzip` on a plain output if it needs a tabix index.
Names that the chosen mode does not match stay unchanged. The tool:

- Renames the CHROM column and `##contig` IDs.
- Applies the same renames to the reference-contig fields of LinGraph headers:
  `##alternativeLocus` ID, `##referenceCoverage` Chrom, and the
  `##pseudoLinearMapping` Reference path. In these headers, it renames only
  names declared by `##contig`, so sample source paths (`Query`,
  `SourceContig`, and alt-layer `Reference`) stay unchanged.
- Leaves nested rows (CHROM = parent ID `INS_`/`DEL_`/`SUB_`/`DUP_`) and
  record IDs unchanged, so nested rows still point to their parents.
- Drops a `##contig` line whose renamed ID repeats an earlier one, with a
  warning.

## Count INS/DEL SVs

`tools/count_sv.sh` counts insertions and deletions of at least 50 bp
(`-m MIN` changes the threshold) in plain or gzip-compressed VCFs:

```bash
bash tools/count_sv.sh NA19240_h1/NA19240_h1.vcf NA19240_h2/NA19240_h2.vcf
```

Columns are exclusive: `INS_only` and `DEL_only` have one part at or above
the threshold, `both` has an inserted and a deleted part at or above it, and
`total = INS_only + DEL_only + 2*both`. Deleted size is
`|SVLEN|` for DEL and `END - POS` otherwise; inserted size is `SVLEN` for INS
and the length of a plain `INFO/SEQ` (or `SVINSSEQ`) otherwise. Records
without `SVLEN` use explicit REF/ALT lengths. Nested grVCF rows are skipped.
