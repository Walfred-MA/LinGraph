# LinGraph utilities

Trio benchmarking scripts are in [`benchmark/`](../benchmark/README.md), with
separate guides for QuickTriocheck and the Truvari wrapper.

## Rerun every saved SNP chromosome on Slurm

`rerun_snp_chroms.py` reads the existing compact SNP manifest and submits one
separate `sbatch` job for every chromosome, including chromosomes whose parts
already exist. Defaults are 32 CPUs, 64G memory and five hours per job, with at
most 40 jobs running concurrently. Jobs beyond that limit depend on the earlier
job in their slot, using `afterany` so another chromosome's failure does not
block them. All jobs are submitted before the launcher exits; it need not stay
running. No job array or shell heredoc is needed.

After syncing the current scripts and stopping the previous controller and SNP
jobs that use these shards:

```bash
python -m pip install 'polars>=1.0'
python tools/rerun_snp_chroms.py \
  --work cohort_calls/tmp/unfinished-graphcigar-8horly9k \
  --account mchaisso_100 --partition qcb
```

The launcher retains source shards, complete parts and lock files. Existing
parts are replaced atomically when their new merge finishes. It refuses an
active chromosome lock instead of bypassing it. The same Python interpreter
runs the workers. Job IDs and commands are saved in `jobs.jsonl` in a new
`snp-rerun-*` directory under the saved workspace, alongside each job's logs.
Submission success does not mean the merges have finished; run SNP concat
after every job succeeds. This command schedules neither SV jobs nor nested
SNP extraction and does not clean up the saved workspace.

## Recover insertion SNPs from a merged SV VCF

`extract_insertion_snps.py` recovers SNPs using the saved per-observation
`EXTENDGRAPHCIGAR`, `ASSEMBLYCONTIG`, and `QUERYCOORD` fields. It does not
repeat SV clustering or alignment and does not need the old shard directory.
For an unaligned insertion, its VCF ID becomes the SNP CHROM. When a portion
aligns to another insertion, the extractor follows the saved named chunks in
`INFO/SEQ` (or a separate INFO CIGAR) to that destination's coordinates. Named
targets already present in a sample's FORMAT CIGAR are also respected. Chains
are followed recursively; reverse mappings complement ALT bases and reverse
reference coordinates. Raw flanks and query-only I/S bases remain on their
source insertion at their original positions because no destination base
exists for them. This includes recursive insertion rows present in the input.
Provide the merged indel VCF as a second input to include below-cutoff insertion
representatives and any coordinate targets defined there.

Calls are compared to the **final destination REF**. A source mismatch that
becomes reference there is omitted; a source self-match that differs from the
destination is called. Reference coverage is projected too. Observations from
different insertions are merged and deduplicated at the destination before
writing the VCF, retaining sample assembly coordinates and provenance. The
extractor uses saved mappings; it does not infer missing mappings by alignment.

Most merged CIGARs contain bare `X` operations, which identify mismatch
positions but do **not** contain the alternate bases. Supply the original
sample assemblies using a two-column, tab-separated file without a header:

```text
HG02698_h1	/path/to/HG02698_h1.fa
HG02698_h2	/path/to/HG02698_h2.fa
```

Sample names must match the VCF columns exactly, and contig names must match
`ASSEMBLYCONTIG`. Relative FASTA paths are resolved against the TSV directory.
Use uncompressed FASTAs with adjacent `.fai` indexes (`samtools faidx FILE.fa`).
The script reads only the necessary query intervals and reverse-complements
negative-strand observations. It retains the pipeline's traversal-boundary
QUERYCOORD and LABEL_H conventions, contributor labels, sample order, and
reference/missing genotype behavior. `TEMPLATEOFFSET` describes the SV's
original breakpoint; it does not shift positions within the insertion.

```bash
python tools/extract_insertion_snps.py \
  -i cohort_calls/cohort.sv.vcf cohort_calls/cohort.indel.vcf \
  --assemblies assemblies.tsv \
  -o cohort_calls/insertion.snp.vcf -t 8
```

Omit `cohort.indel.vcf` if only SV insertion representatives are wanted. The
result is a complete SNP VCF, sorted by destination contig, position, REF and ALT.
It can be combined with the existing SNP VCF later, but the insertion contig
headers must be added and sample columns must have the same order; do not
concatenate two complete VCF headers. Alternatively, generate a combined file
directly with `--append-to`:

```bash
python tools/extract_insertion_snps.py \
  -i cohort_calls/cohort.sv.vcf cohort_calls/cohort.indel.vcf \
  --assemblies assemblies.tsv --append-to cohort_calls/cohort.snp.vcf \
  -o cohort_calls/cohort.with_insertions.snp.vcf -t 8
```

This copies the existing SNP records first, merges contig headers, and then
appends recovered insertion SNPs. Inputs are retained. Differing sample order,
conflicting contig lengths, and insertion contigs already present in the SNP
records cause an error to prevent an incorrect append.

Plain and gzip VCF inputs/outputs are supported. Each worker handles bounded
batches and one source insertion at a time, saving projected observations as
21-byte records with separate metadata. Final merging loads one destination
locus at a time. A temporary disk index resolves encoded `INFO/SEQ` dependencies
and coordinate mappings, including references to later VCF rows. For external
graph-sequence dependencies, provide their indexed FASTAs with repeated
`-a FILE.fa` options. Use `--tmpdir DIR` to place temporary data on scratch
storage; compressed inputs are decompressed there once. Temporary files are
removed on success or a handled failure. Existing output requires `--force`
and is replaced only after successful extraction. Missing bases, missing
assembly contigs, inconsistent alignments and duplicate insertion IDs fail
explicitly instead of producing incomplete SNP calls. `--assemblies` may be
omitted only when all required mismatch bases are embedded in the VCF.

## Prepare assembly FASTAs

Use `prepare_assemblies.py` for soft masking, optional contig-name fixing,
and indexing. It writes prepared copies in a separate output directory:

```bash
python tools/prepare_assemblies.py \
  -i assembly.fa --name HG002_h1 \
  -O prepared_assemblies --contignamefix
```

For name fixing alone, use `namecontigsfix.py` with a separate output:

```bash
python tools/namecontigsfix.py \
  -i assembly.fa -n 'HG002#1' -o assembly.namefixed.fa
samtools faidx assembly.namefixed.fa
```

Keep native contig names for `CHM13_h1` and `HG38_h1`; omit name fixing for
those references. LinGraph and the standalone sample caller check inputs
without modifying assemblies or indexes and stop with a preparation command
when the checks fail. See [assembly preparation](assembly_preparation.md)
for batch preparation, masking options, and dependencies.

## Convert grVCF to standard VCF

```bash
python tools/grvcf_to_vcf.py \
  -i calls/cohort.sv.vcf \
  -r reference.fa \
  -o calls/cohort.sv.standard.vcf
```

Use the same reference FASTA used for calling, with its adjacent `reference.fa.fai`.
Individual VCFs and merged `.snp.vcf`, `.sv.vcf`, or `.indel.vcf` files are
accepted, including older combined `.all.vcf` files. Input can be plain VCF or
gzip-compressed VCF. Convert one file per invocation; use separate output names
for the SV, indel, and SNP files produced by a current `--all` run.

The converter:

- Writes explicit REF/ALT bases for SNPs, insertions, deletions, and replacements.
  A grVCF `<SUB>` becomes one replacement allele.
- Uses the row's representative sequence and displayed breakpoint. Each
  carrier's GT continues to indicate support for that merged representative;
  the per-observation alleles and `TEMPLATEOFFSET` values are not expanded.
- Saves rows whose CHROM is absent from the reference FASTA to
  `OUTPUT.unplaced.gr.vcf`; use `--unplaced FILE` to choose another path.
  Auxiliary FASTAs supply decoding sequences, not new output reference contigs.
- Preserves sample names, genotype ploidy, phasing, missing calls, IDs, QUAL,
  and FILTER. Old packed `HSV` and `HSNP` sample fields are accepted.
- Replaces graph-specific INFO/FORMAT annotations with standard `GT`, `AC`,
  `AN`, `AF`, and `NS`, computing the counts from the output genotypes.
- Checks REF against the FASTA, supplies padding bases for indels, and sorts
  output by FASTA contig order and POS. A grVCF insertion or deletion at
  boundary zero receives right padding at VCF position 1.
- Decodes graph-encoded `INFO/SEQ` or `INFO/EXTENDGRAPHCIGAR` using literal
  query bases, named FASTA targets, and other grVCF allele IDs. Forward
  references to later rows are supported. Missing sequence, ambiguous `M`
  operations, conflicting IDs, and reference mismatches cause an error.

For encoded sequences referencing additional graph or alternative FASTAs,
provide each catalog with `-a`. Each needs an adjacent `.fai`:

```bash
python tools/grvcf_to_vcf.py \
  -i calls/cohort.sv.vcf -r reference.fa \
  -a savegraph/summary/alternatives.fasta \
  -o calls/cohort.sv.standard.vcf
```

Plain `.vcf` output uses Python and the repository's existing sequence helpers.
Use `--normalize` for left alignment with `bcftools norm`. Compressed output
uses BGZF and includes a CSI index; it also requires `bcftools`:

```bash
python tools/grvcf_to_vcf.py \
  -i calls/cohort.sv.vcf -r reference.fa \
  -o calls/cohort.sv.standard.vcf.gz --normalize
```

Temporary records and resolved graph sequences are stored on disk beside the
output, limiting memory use to individual records and sequences. Existing
output requires `--force`. Input files are never overwritten; final output
is installed after all conversion and optional normalization steps succeed.

The output follows the
[VCF 4.3 allele and genotype conventions](https://samtools.github.io/hts-specs/VCFv4.3.pdf).
