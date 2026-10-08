# NA19240 trio: SV calling with Win50KGraph

This tutorial calls SVs in the YRI trio assembled by HGSVC3: NA19240 (child),
NA19238 (mother) and NA19239 (father), two haplotypes each. Every haplotype is
called independently (`individual` mode) against the precomputed Win50KGraph on
CHM13, starting from the raw assemblies (lazy preparation), and the calls are
then checked for trio consistency, coverage and SV counts.

To build a new pangenome graph from the same trio instead, see
[NA19240trio_graph](../NA19240trio_graph/README.md).

Run every command from one demo folder. Plan for about **300 GB** of disk.
Long commands write their progress to `logs/`; follow one with
`tail -f logs/NAME.log`. Steps marked *(hours)* belong on a compute node.

```bash
mkdir NA19240trio_SVcall && cd NA19240trio_SVcall
mkdir logs
```

## 1. Install LinGraph

```bash
git clone https://github.com/Walfred-MA/LinGraph.git
cd LinGraph
python install.py --conda-env LinGraph -j 16 > ../logs/install.log 2>&1
conda activate LinGraph
python install.py --check-only
cd ..
```

`install.py` creates the `LinGraph` conda environment and compiles the native
tools. Activate the environment in every new shell.

## 2. Prepare the CHM13 reference

```bash
mkdir references && cd references
curl -L -O https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/009/914/755/GCF_009914755.1_T2T-CHM13v2.0/GCF_009914755.1_T2T-CHM13v2.0_genomic.fna.gz
python ../LinGraph/tools/prepare_assemblies.py \
  -i GCF_009914755.1_T2T-CHM13v2.0_genomic.fna.gz --name CHM13_h1 \
  -O chm13 --remask --threads 16 > ../logs/prepare_chm13.log 2>&1
cd ..
```

This writes `references/chm13/GCF_009914755.1_T2T-CHM13v2.0_genomic.fna` with
its index, soft-masked with WindowMasker as used to build Win50KGraph. CHM13
keeps its native contig names (`NC_060925.1` ...).

## 3. Download and reconstruct Win50KGraph *(hours)*

```bash
curl -L -C - https://ndownloader.figshare.com/files/69451449 -o Win50KGraph.tar.gz
tar -xzf Win50KGraph.tar.gz
python LinGraph/scripts/reconstruct_local_graph_folders.py \
  -i cohort_minsetref_v3/summary \
  -r references/chm13/GCF_009914755.1_T2T-CHM13v2.0_genomic.fna \
  --reference-haplotype CHM13_h1 \
  -o Win50KGraph -j 32 --resume > logs/reconstruct.log 2>&1
```

Keep `cohort_minsetref_v3/`: `Win50KGraph/summary` and
`Win50KGraph/references` link back to it. Repeat the command to resume.

## 4. Download the trio assemblies

```bash
mkdir assemblies && cd assemblies
URL=https://ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/HGSVC3/working/20240201_verkko_batch3/assemblies
for sample in NA19240 NA19238 NA19239; do
  for hap in 1 2; do
    curl -L -C - -O $URL/$sample/$sample.vrk-ps-sseq.asm-hap$hap.fasta.gz
    gzip -d $sample.vrk-ps-sseq.asm-hap$hap.fasta.gz
  done
done
cat > assemblies.list <<'EOF'
NA19240_h1 NA19240.vrk-ps-sseq.asm-hap1.fasta
NA19240_h2 NA19240.vrk-ps-sseq.asm-hap2.fasta
NA19238_h1 NA19238.vrk-ps-sseq.asm-hap1.fasta
NA19238_h2 NA19238.vrk-ps-sseq.asm-hap2.fasta
NA19239_h1 NA19239.vrk-ps-sseq.asm-hap1.fasta
NA19239_h2 NA19239.vrk-ps-sseq.asm-hap2.fasta
EOF
cd ..
```

LinGraph reads uncompressed FASTAs only, so the downloads are decompressed.
The assemblies are used as downloaded: their contigs are named like
`haplotype1-0000001`, they are not soft-masked and have no index.

## 5. Call SVs (lazy preparation) *(hours)*

```bash
python LinGraph/scripts/LinGraph.py individual \
  -I assemblies/assemblies.list -G Win50KGraph \
  -r references/chm13/GCF_009914755.1_T2T-CHM13v2.0_genomic.fna --reference-name CHM13_h1 \
  --reference-caches Win50KGraph/references/CHM13_h1_rig \
  -O calls_chm13 -t 64 --force-prepare > logs/individual_chm13.log 2>&1
```

On a SLURM cluster, add `--slurm --slurm-account ACCOUNT --slurm-partition
PARTITION`: each haplotype is called in its own job, all six at once (about one
hour each with 64 CPUs). Keep the launcher running, e.g. with `nohup ... &`.

**Lazy mode.** LinGraph requires every contig to carry its haplotype prefix
(`NA19240#1#haplotype1-0000001`), so contig names stay unique across
assemblies and every call can be traced to its original contig, and it
requires soft-masked assemblies, which help the aligners. Without
`--force-prepare` it stops and prints a preparation command for each
assembly. With it, LinGraph prepares fixed copies (prefixed, masked, indexed)
in `calls_chm13/PreparedAssemblies/` and deletes them after the run; the
original files are never modified. The copies need extra temporary disk, one
full copy of every assembly, so use lazy mode for a few samples and prepare
large cohorts beforehand with `tools/prepare_assemblies.py` (see
[NA19240trio_graph](../NA19240trio_graph/README.md), step 3).

`--reference-caches` reuses the CHM13 alignment shipped with Win50KGraph. Each
haplotype's calls are in `calls_chm13/samples/NAME/NAME.vcf`, with coverage
reports beside them.

## 6. Benchmark

### Trio consistency (Truvari)

`truvari_trio.sh` runs eight Truvari comparisons, each child haplotype against
the four parental haplotypes, inside the reference intervals covered by all
parents, and reports the child SVs found in neither parent. It is run at
`--distance 0` (identical breakpoints) and at Truvari's default 500 bp:

```bash
conda create -y -n truvari -c conda-forge -c bioconda truvari htslib > logs/truvari_env.log 2>&1
conda activate truvari
CALLS=$PWD/calls_chm13/samples; mkdir -p calls_chm13_bench
for d in 0 500; do
  bash LinGraph/benchmark/truvari_trio.sh \
    --child  $CALLS/NA19240_h1/NA19240_h1.vcf $CALLS/NA19240_h2/NA19240_h2.vcf \
    --mother $CALLS/NA19238_h1/NA19238_h1.vcf $CALLS/NA19238_h2/NA19238_h2.vcf \
    --father $CALLS/NA19239_h1/NA19239_h1.vcf $CALLS/NA19239_h2/NA19239_h2.vcf \
    --distance $d --jobs 8 --fp calls_chm13_bench/NA19240.fp.d$d.vcf \
    > calls_chm13_bench/trio.truvari.d$d.txt 2> logs/truvari_d$d.log
done
conda activate LinGraph
```

Results (CHM13):

| | h1, 0 bp | h2, 0 bp | Both, 0 bp | h1, default | h2, default | Both, default |
| --- | ---: | ---: | ---: | ---: | ---: | ---: |
| Child SVs | 14,977 | 15,482 | 30,459 | 14,977 | 15,482 | 30,459 |
| Found in mother | 10,563 | 12,151 | 22,714 | 11,020 | 12,475 | 23,495 |
| Found in father | 11,838 | 10,823 | 22,661 | 12,165 | 11,290 | 23,455 |
| Found in both parents | 7,575 | 7,613 | 15,188 | 8,302 | 8,358 | 16,660 |
| Found in either parent | 14,826 | 15,361 | 30,187 | 14,883 | 15,407 | 30,290 |
| Found in neither parent | 151 | 121 | 272 | 94 | 75 | 169 |
| Error rate | 1.01% | 0.78% | 0.89% | 0.63% | 0.48% | 0.55% |
| Trio consistency | 98.99% | 99.22% | 99.11% | 99.37% | 99.52% | 99.45% |

A child SV found in neither parent is not necessarily an error: de novo
variants and SVs missed in a parent land there too. The calls without a
parental match are in `calls_chm13_bench/NA19240.fp.d*.vcf`.

### Coverage

Each haplotype's `NAME.coverage.summary.tsv` reports the reference and
assembly bases its calls do not cover:

```bash
cat calls_chm13/samples/NA19240_h1/NA19240_h1.coverage.summary.tsv \
    calls_chm13/samples/NA19240_h2/NA19240_h2.coverage.summary.tsv
```

| Scope | Category | NA19240_h1 total (Mb) | NA19240_h1 unmasked (Mb) | NA19240_h2 total (Mb) | NA19240_h2 unmasked (Mb) |
| --- | --- | ---: | ---: | ---: | ---: |
| Reference | Not covered, total | 325.04 | 61.92 | 274.50 | 29.74 |
| Reference | chrY | 61.19 | 14.30 | 61.05 | 14.22 |
| Reference | Other | 263.85 | 47.61 | 213.45 | 15.52 |
| Assembly | Not covered, total | 167.02 | 7.38 | 164.91 | 6.68 |
| Assembly | N gaps | 1.68 | 0.05 | 2.66 | 0.05 |
| Assembly | Scaffold edges | 24.08 | 3.24 | 28.85 | 3.34 |
| Assembly | Other | 141.26 | 4.08 | 133.40 | 3.29 |

NA19240 is female, so all of chrY is uncovered. Most other uncovered bases
are soft-masked repeats (the difference between the two columns). Satellites
are counted separately only when calling with `--satellite-bed`, for example
the CHM13 constitutive satellite BED. `NAME.coverage.missing_reference.bed`
and `NAME.coverage.missing_query.bed` list the uncovered intervals.

### SV counts

`count_sv.sh` counts SVs of at least 50 bp in each VCF:

```bash
bash LinGraph/tools/count_sv.sh \
  calls_chm13/samples/NA19240_h1/NA19240_h1.vcf calls_chm13/samples/NA19240_h2/NA19240_h2.vcf
```

| Haplotype | Insertions only | Deletions only | Both | Total |
| --- | ---: | ---: | ---: | ---: |
| NA19240_h1 | 7,499 | 7,265 | 1,505 | 17,774 |
| NA19240_h2 | 7,858 | 7,500 | 1,569 | 18,496 |

*Both* counts replacements whose inserted and deleted parts are each at least
50 bp; *Total* counts them twice (insertions only + deletions only + 2 ×
both). The Truvari totals above are lower because Truvari is restricted to
the reference intervals covered by all four parental haplotypes.
