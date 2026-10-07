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

## Merge individual grVCFs

`tools/merge_grvcfs.py` runs only LinGraph's merge step on per-sample grVCFs
and writes `cohort.sv.vcf`, `cohort.indel.vcf`, `cohort.snp.vcf` and
`cohort.samples.headers.gz` to the output folder:

```bash
python tools/merge_grvcfs.py -I vcfs.list -O merged -t 16
```

The merge runs as stages, each its own process: SV scan, SV merge per
chromosome, SV concat, SNP prepare, SNP merge per chromosome, SNP concat and
publish. Without `--slurm` they run locally one after another (chromosomes one
at a time, `-t` workers). With `--slurm` each stage is a Slurm job: run-once
stages get `-t` CPUs and `--slurm-memory` (default 64G below 100 input grVCFs,
128G from 100 on), and chromosome stages run as up to `--slurm-jobs` (default
20) jobs at once, sized like the pipeline (SV: min(16, `-t`) CPUs and 2G per CPU,
or 64G from 100 input grVCFs on, chr1 doubled; SNP: 32 CPUs, 64G). Finished stages and chromosomes are recorded,
so repeating a failed command reruns only what did not finish. Both modes write
the same files as the single-process merge.

```bash
python tools/merge_grvcfs.py -I vcfs.list -O merged -t 32 --slurm --slurm-jobs 20 \
    --slurm-account ACCOUNT --slurm-partition PARTITION
```

By default it writes the CIGAR version: each sample keeps its own breakpoint
(`TEMPLATEOFFSET`), size and alignment to its row's representative. With
`--exact query_paths.txt --reference-fasta reference.fa` members are realigned
against their assemblies and moved onto the representative's breakpoint;
their other bases become nested, `_F` and `_S` rows. A
`##graphvcfmergeVersion=cigar|exact` header line records which one a file is.

`cohort.samples.headers.gz` keeps every sample's `##referenceCoverage` and
`##pseudoLinearMapping` lines, which merged headers drop. The converter below
needs it to give split samples their coverage.

## Split merged grVCFs, or convert exact to CIGAR

`tools/convert_merged_grvcf.py split` turns merged grVCFs of either version
back into one grVCF per sample; `to-cigar` turns an exact merge into a CIGAR
merge (split, then merge again without realignment). Give all three merged
files of the cohort and the merge's reference FASTA:

```bash
python tools/convert_merged_grvcf.py split \
  -v merged/cohort.sv.vcf merged/cohort.indel.vcf merged/cohort.snp.vcf \
  -r reference.fa -o individual
python tools/convert_merged_grvcf.py to-cigar \
  -v merged/cohort.sv.vcf merged/cohort.indel.vcf merged/cohort.snp.vcf \
  -r reference.fa -O merged_cigar -t 16
```

To add samples to a cohort: split it (or convert it to CIGAR first), then
merge the split files and the new samples' grVCFs with `merge_grvcfs.py`.

- Each call is rebuilt from its row, the nested rows the sample carries and
  its insertion SNPs. From an exact merge, a row, its `_F` pieces and the
  sample's SNPs between them become one call again; `_S` rows stay separate.
- Rebuilt samples reconstruct their assemblies exactly
  (`scripts/tools/check_vcf_lossless.py`). Split rows keep their merged row
  IDs, so merging them again keeps the row names.
- An exact merge keeps no sample breakpoint inside repeated bases. A rebuilt
  insertion takes the equivalent placement that best matches its row's
  representative; a deletion may come back at another equivalent position.
- Files written before the version line existed need `--input-version`.
- Without `cohort.samples.headers.gz` (pass another path with
  `--sample-headers`) the tool warns: the split files then have no coverage,
  so re-merged samples are `.` instead of `0` outside their own calls, and an
  `--exact` merge cannot place them.

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

- Removes nested records on `I_`, `D_`, `SUB_`, or `DUP_` parents (`INS_`/`DEL_`
  in older runs), and
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
- Leaves nested rows (CHROM = parent ID `I_`/`D_`/`SUB_`/`DUP_`, or
  `INS_`/`DEL_` in older runs) and
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
