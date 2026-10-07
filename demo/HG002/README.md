# Getting started: HG002 on CHM13 and HG19

This tutorial calls SVs on both haplotypes of the GIAB/T2T HG002 Q100
assembly with the precomputed Win50KGraph, on two references:

- **CHM13**, with the reference cache shipped in Win50KGraph.
- **HG19** (UCSC GRCh37), building its reference cache first and reusing the
  CHM13 run's graph alignments.

It then renames the HG19 contigs back to `chr1` … `chrY` and benchmarks the
calls against GIAB: NIST SV v0.6 Tier 1 on HG19 and the Q100 v5.0q structural
variants on CHM13.

Run every command from one demo folder. Plan for about **200 GB** of disk.
Long commands write their progress to `logs/`; follow one with
`tail -f logs/NAME.log`. Steps that take hours are marked; run them on a
compute node, for example by prefixing the command with `nohup` and
appending `&`.

```bash
mkdir hg002_demo && cd hg002_demo
mkdir logs
```

## 1. Install LinGraph

```bash
git clone https://github.com/Walfred-MA/LinGraph.git
cd LinGraph
python3 install.py --conda-env LinGraph -j 16 > ../logs/install.log 2>&1
conda activate LinGraph
python3 install.py --check-only
cd ..
```

`install.py` creates the `LinGraph` conda environment (using mamba when
available), installs the dependencies and compiles the native tools into
`LinGraph/scripts/`. Activate the environment in every new shell.

## 2. Prepare the references

LinGraph needs references as uncompressed, indexed FASTAs.

- **CHM13**: NCBI RefSeq T2T-CHM13v2.0, 24 chromosomes (`NC_060925.1` …),
  the sequences of the Win50KGraph CHM13 cache.
- **HG19**: UCSC hg19, main chromosomes only (chr1–22, X, Y). A reference
  other than CHM13 or GRCh38 uses the `NAME#1#` contig prefix, so
  `--contignamefix` renames `chr1` to `HG19#1#chr1`.

```bash
mkdir references && cd references

curl -L -O https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/009/914/755/GCF_009914755.1_T2T-CHM13v2.0/GCF_009914755.1_T2T-CHM13v2.0_genomic.fna.gz
python3 ../LinGraph/tools/prepare_assemblies.py \
  -i GCF_009914755.1_T2T-CHM13v2.0_genomic.fna.gz --name CHM13_h1 \
  -O chm13 --remask --threads 16 > ../logs/prepare_chm13.log 2>&1

curl -L -O https://hgdownload.soe.ucsc.edu/goldenPath/hg19/bigZips/hg19.fa.gz
gzip -dc hg19.fa.gz > hg19.fa
samtools faidx hg19.fa
samtools faidx hg19.fa chr{1..22} chrX chrY > hg19_main.fa
python3 ../LinGraph/tools/prepare_assemblies.py \
  -i hg19_main.fa --name HG19_h1 -O hg19 \
  --contignamefix --remask --threads 16 > ../logs/prepare_hg19.log 2>&1

cd ..
```

`prepare_assemblies.py` soft-masks with WindowMasker (`--remask` replaces the
downloads' own masking with WindowMasker's, as used to build Win50KGraph),
writes an uncompressed FASTA and runs `samtools faidx`. While masking hg19,
WindowMasker prints many `ambiguous nucleotides` warnings for contigs that
start with long runs of `N`; they are harmless. The results are
`references/chm13/GCF_009914755.1_T2T-CHM13v2.0_genomic.fna` and
`references/hg19/hg19_main.fa`. Keep `references/hg19.fa` (plain `chr` names)
for benchmarking in step 8.

## 3. Download and reconstruct Win50KGraph

```bash
curl -L -C - https://ndownloader.figshare.com/files/69451449 -o Win50KGraph.tar.gz
tar -xzf Win50KGraph.tar.gz

python3 LinGraph/scripts/reconstruct_local_graph_folders.py \
  -i cohort_minsetref_v3/summary \
  -r references/chm13/GCF_009914755.1_T2T-CHM13v2.0_genomic.fna \
  --reference-haplotype CHM13_h1 \
  -o Win50KGraph -j 32 --resume > logs/reconstruct.log 2>&1
```

Reconstruction (hours) rebuilds every local graph into `Win50KGraph/`. Keep
the extracted `cohort_minsetref_v3/` folder: `Win50KGraph/summary` and
`Win50KGraph/references` link back to it. Repeat the command to resume.

## 4. Prepare the HG002 assemblies

Download the HG002 Q100 v1.2 haplotypes. `pat` holds chr1–22 and chrY, `mat`
chr1–22 and chrX. This tutorial uses h1 = paternal, h2 = maternal:

```bash
mkdir assemblies && cd assemblies
for h in pat mat; do
  curl -L -C - -O https://human-pangenomics.s3.amazonaws.com/T2T/HG002/assemblies/hg002v1.2.$h.fasta.gz
done

cat > raw_assemblies.list <<'EOF'
HG002_h1 hg002v1.2.pat.fasta.gz
HG002_h2 hg002v1.2.mat.fasta.gz
EOF

python3 ../LinGraph/tools/prepare_assemblies.py \
  -q raw_assemblies.list -O prepared --contignamefix -j 2 --threads 8 \
  > ../logs/prepare_assemblies.log 2>&1
cd ..
```

The Q100 assembly is not soft-masked, so preparation masks it with
WindowMasker and renames each contig with its haplotype prefix, e.g.
`chr1_PATERNAL` to `HG002#1#chr1_PATERNAL`. The calling input list is
`assemblies/prepared/query_paths.prepared.txt`.

## 5. Call SVs on CHM13

Call both haplotypes against CHM13 (hours):

```bash
python3 LinGraph/scripts/LinGraph.py singular \
  -I assemblies/prepared/query_paths.prepared.txt \
  -G Win50KGraph \
  -r references/chm13/GCF_009914755.1_T2T-CHM13v2.0_genomic.fna --reference-name CHM13_h1 \
  --reference-caches Win50KGraph/references/CHM13_h1_rig \
  -O HG002_calls_chm13 -t 64 > logs/singular_chm13.log 2>&1
```

Each haplotype is called on its own into
`HG002_calls_chm13/samples/HG002_hN/HG002_hN.vcf`, with a coverage report
`HG002_hN.coverage.summary.tsv` next to it. `--reference-caches` reuses the
CHM13 alignment shipped with the package. Calls are on CHM13's NCBI names
(`NC_060925.1` …), the names the CHM13 truth set gets in step 8.

## 6. Call SVs on HG19

Win50KGraph ships caches for CHM13 and GRCh38 only. Build the HG19 cache once
(hours, like calling one haplotype); it is written to
`Win50KGraph/references/HG19_h1_rig`, and no sample is called. Repeat the
command to resume:

```bash
python3 LinGraph/scripts/LinGraph.py singular \
  -G Win50KGraph -r references/hg19/hg19_main.fa --reference-name HG19_h1 \
  --reference-only -O hg19_cache -t 64 > logs/hg19_cache.log 2>&1
```

Then call both haplotypes on HG19. The graph alignment of each haplotype does
not depend on the reference, so `--reuse-alignments` links it from the CHM13
run and only the reference-specific stages run:

```bash
python3 LinGraph/scripts/LinGraph.py singular \
  -I assemblies/prepared/query_paths.prepared.txt \
  -G Win50KGraph \
  -r references/hg19/hg19_main.fa --reference-name HG19_h1 \
  --reuse-alignments HG002_calls_chm13 \
  -O HG002_calls_hg19 -t 64 > logs/singular_hg19.log 2>&1
```

The HG19 cache is found in `Win50KGraph/references/` automatically. Without
`--reuse-alignments`, the haplotypes are aligned to the graph again (hours).

## 7. Rename the HG19 contigs to chr1 … chrY

HG19 calls are on the prepared names (`HG19#1#chr1`). GIAB uses plain names,
so rename the finished VCFs with `tools/vcf_chromfix.py`:

```bash
for h in h1 h2; do
  python3 LinGraph/tools/vcf_chromfix.py \
    -i HG002_calls_hg19/samples/HG002_$h/HG002_$h.vcf \
    -o HG002_calls_hg19/samples/HG002_$h/HG002_$h.chr.vcf --noprefix
done
```

`--noprefix` drops everything up to the last `#` in reference names, in
records and in LinGraph's headers. Keep the original files as input to
LinGraph's own tools; they match the prepared FASTA names.

## 8. Benchmark against GIAB

The benchmark needs truvari and bcftools. Install them in their own
environment:

```bash
conda create -y -n truvari -c conda-forge -c bioconda truvari bcftools htslib > logs/truvari_env.log 2>&1
conda activate truvari
```

`tools/BCFtoolMergeTwoHaplotypes.sh` turns the two haplotype VCFs into one
site-only VCF, `HG002.LinGraph.vcf.gz`: each is sorted and stripped of
genotypes, then they are combined with exact duplicates removed. It writes
into the current folder, so each reference gets its own benchmark folder.

### HG19: GIAB NIST SV v0.6 Tier 1

Download the GIAB SV call set and its confident regions. They use `1` … `22`,
`X`, `Y`; rename them to `chr1` …:

```bash
mkdir bench_hg19 && cd bench_hg19
GIAB=https://ftp-trace.ncbi.nlm.nih.gov/ReferenceSamples/giab/release/AshkenazimTrio/HG002_NA24385_son/NIST_SV_v0.6
wget $GIAB/HG002_SVs_Tier1_v0.6.vcf.gz
wget $GIAB/HG002_SVs_Tier1_v0.6.vcf.gz.tbi
wget $GIAB/HG002_SVs_Tier1_v0.6.bed

for i in $(seq 1 22); do echo -e "$i\tchr$i"; done > giab_to_chr.txt
echo -e "X\tchrX" >> giab_to_chr.txt
echo -e "Y\tchrY" >> giab_to_chr.txt

bcftools annotate --rename-chrs giab_to_chr.txt \
  HG002_SVs_Tier1_v0.6.vcf.gz -Oz -o HG002_SVs_Tier1_v0.6.hg19.vcf.gz
tabix -f -p vcf HG002_SVs_Tier1_v0.6.hg19.vcf.gz
awk 'BEGIN{OFS="\t"} {$1="chr"$1; print}' HG002_SVs_Tier1_v0.6.bed > HG002_SVs_Tier1_v0.6.hg19.bed
```

Combine the haplotypes and benchmark:

```bash
bash ../LinGraph/tools/BCFtoolMergeTwoHaplotypes.sh \
  ../HG002_calls_hg19/samples/HG002_h1/HG002_h1.chr.vcf \
  ../HG002_calls_hg19/samples/HG002_h2/HG002_h2.chr.vcf

truvari bench \
  -b HG002_SVs_Tier1_v0.6.hg19.vcf.gz \
  -c HG002.LinGraph.vcf.gz \
  -f ../references/hg19.fa \
  --includebed HG002_SVs_Tier1_v0.6.hg19.bed \
  --passonly -r 2000 -C 5000 -p 0 -P 0.7 -O 0 --pick multi \
  -o truvari_HG002_hg19

cat truvari_HG002_hg19/summary.json
cd ..
```

### CHM13: GIAB Q100 v5.0q, non-difficult regions

Download the Q100 v5.0q structural-variant benchmark on CHM13v2.0, and rename
its chromosomes to CHM13's NCBI names to match the calls. Check its names
first with `bcftools view -h HG002_CHM13v2.0_v5.0q_stvar.vcf.gz | grep -m3 '^##contig'`;
the table below maps `chr1` … `chrY`:

```bash
mkdir bench_chm13 && cd bench_chm13
Q100=https://ftp-trace.ncbi.nlm.nih.gov/ReferenceSamples/giab/release/AshkenazimTrio/HG002_NA24385_son/v5.0q
wget $Q100/HG002_CHM13v2.0_v5.0q_stvar.vcf.gz
wget $Q100/HG002_CHM13v2.0_v5.0q_stvar.vcf.gz.tbi

for i in $(seq 1 22); do echo -e "chr$i\tNC_0$((60924 + i)).1"; done > chr_to_nc.txt
echo -e "chrX\tNC_060947.1" >> chr_to_nc.txt
echo -e "chrY\tNC_060948.1" >> chr_to_nc.txt

bcftools annotate --rename-chrs chr_to_nc.txt \
  HG002_CHM13v2.0_v5.0q_stvar.vcf.gz -Oz -o HG002_CHM13v2.0_v5.0q_stvar.NC.vcf.gz
tabix -f -p vcf HG002_CHM13v2.0_v5.0q_stvar.NC.vcf.gz
```

The non-difficult regions on CHM13, already with NCBI names, ship with
LinGraph:

```bash
tar -xJf ../LinGraph/benchmark/HG002_Q100.nondifficult.NC.bed.tar.xz
```

Combine the haplotypes and benchmark:

```bash
bash ../LinGraph/tools/BCFtoolMergeTwoHaplotypes.sh \
  ../HG002_calls_chm13/samples/HG002_h1/HG002_h1.vcf \
  ../HG002_calls_chm13/samples/HG002_h2/HG002_h2.vcf

truvari bench \
  -b HG002_CHM13v2.0_v5.0q_stvar.NC.vcf.gz \
  -c HG002.LinGraph.vcf.gz \
  -f ../references/chm13/GCF_009914755.1_T2T-CHM13v2.0_genomic.fna \
  --includebed HG002_Q100.nondifficult.NC.bed \
  -r 2000 -C 5000 --pick multi \
  -o truvari_HG002_nondifficult

cat truvari_HG002_nondifficult/summary.json
cd ..
conda activate LinGraph
```

`benchmark/BenchGIAB.sh` runs the same CHM13 benchmark and also reports it
without GIAB's most extreme truth variants.

## 9. Count SVs and check coverage

Count SVs of at least 50 bp per haplotype:

```bash
for ref in hg19 chm13; do
  for h in h1 h2; do
    bash LinGraph/tools/count_sv.sh HG002_calls_$ref/samples/HG002_$h/HG002_$h.vcf
  done
done
```

`total` counts insertions plus deletions; a call with both an inserted and a
deleted part of at least 50 bp counts twice (`total = INS_only + DEL_only +
2*both`).

Each haplotype's coverage report lists what the calls do not cover:

```bash
cat HG002_calls_hg19/samples/HG002_h1/HG002_h1.coverage.summary.tsv
```

`reference total_unannotated` is the reference sequence this haplotype has
no calls on; `query total_unannotated` is the assembly sequence placed on no
reference interval (N gaps, scaffold edges, sequence that did not align).
The `unmasked_bases` column leaves out soft-masked repeats, mostly
constitutive heterochromatin. Below, the reference side uses
`reference total_unannotated`, and the assembly side `query other`, which
also leaves out N gaps and scaffold edges.

## 10. Results

Results from running this tutorial:

| Reference | Haplotype | Precision vs GIAB | Recall vs GIAB | Unannotated assembly bases (excl. constitutive heterochromatin) | Unannotated reference bases (excl. constitutive heterochromatin) | SVs ≥ 50 bp |
| --- | --- | --- | --- | --- | --- | --- |
| HG19 | h1 (pat) | 0.938 | 0.978 | 13,152,202 | 99,590,004 | 17,176 |
| HG19 | h2 (mat) | 0.938 | 0.978 | 12,245,193 | 25,756,452 | 17,906 |
| CHM13 | h1 (pat) | 0.986 | 0.966 | 4,303,423 | 126,791,810 | 14,230 |
| CHM13 | h2 (mat) | 0.986 | 0.966 | 4,207,899 | 52,815,189 | 14,727 |

- **Precision and recall** are for the diploid call set (h1 and h2
  combined), so both haplotypes share them.
  - HG19 against NIST SV v0.6 Tier 1: 9,432 truth SVs matched, 653 false
    positives, 209 false negatives (F1 0.958).
  - CHM13 against Q100 v5.0q in non-difficult regions: 2,906 truth SVs
    matched, 39 false positives, 103 false negatives (F1 0.976).
- **GIAB is not error-free.** The v0.6 HG19 set is reported to contain about
  5% errors, so the HG19 numbers are a lower bound on accuracy: some "false
  positives" and "false negatives" are errors in the truth set. The newer Q100
  benchmark on CHM13 gives higher agreement.
- **Unannotated bases** are `unmasked_bases` in the coverage reports, which
  leave out constitutive heterochromatin and other soft-masked repeats.
  - Assembly: the `query other` row (e.g. `query other 193 153097044
    4207899 ...` gives 4,207,899), which also leaves out N gaps and scaffold
    edges. About 4 Mb of each haplotype is placed nowhere on CHM13, about
    12–13 Mb on HG19.
  - Reference: the `reference total_unannotated` row. It includes whole
    chromosomes a haplotype lacks: h1 is paternal (chrY, no chrX) and h2
    maternal (chrX, no chrY), which accounts for most of the difference
    between them.
- **SV counts** are higher on HG19 than on CHM13, almost entirely insertions
  (h1: 9,347 insertion-only on HG19 versus 6,136 on CHM13; deletions are
  similar), likely sequence that GRCh37 lacks or represents with a minor
  allele.
