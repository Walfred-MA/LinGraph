# Getting started: HG002 on CHM13, reusing graph alignments

This tutorial calls SVs on both haplotypes of the GIAB/T2T HG002 Q100
assembly with the precomputed Win50KGraph, on CHM13, and benchmarks them
against GIAB's Q100 v5.0q structural variants.

It also shows how to **switch reference without aligning the haplotypes to
the graph again**. The graph alignment of a haplotype does not depend on the
reference, so a call on CHM13 can reuse the alignments of the
[HG19 example](../HG002_HG19/README.md); only the reference-specific stages
run.

- **After the HG19 example:** work in its folder (`hg002_demo`). Steps 1–4
  below are already done there; go to step 5.
- **On its own:** follow steps 1–4, then call in step 5 without
  `--reuse-alignments` (the haplotypes are aligned to the graph, hours).

Plan for about **200 GB** of disk. Long commands write their progress to
`logs/`; follow one with `tail -f logs/NAME.log`. Steps that take hours are
marked; run them on a compute node, for example by prefixing the command with
`nohup` and appending `&`.

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

## 2. Prepare the CHM13 reference

LinGraph needs references as uncompressed, indexed FASTAs. Use NCBI RefSeq
T2T-CHM13v2.0, 24 chromosomes (`NC_060925.1` …), the sequences of the
Win50KGraph CHM13 cache:

```bash
mkdir references && cd references
curl -L -O https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/009/914/755/GCF_009914755.1_T2T-CHM13v2.0/GCF_009914755.1_T2T-CHM13v2.0_genomic.fna.gz
python3 ../LinGraph/tools/prepare_assemblies.py \
  -i GCF_009914755.1_T2T-CHM13v2.0_genomic.fna.gz --name CHM13_h1 \
  -O chm13 --remask --threads 16 > ../logs/prepare_chm13.log 2>&1
cd ..
```

`prepare_assemblies.py` soft-masks with WindowMasker (`--remask` replaces the
download's own masking with WindowMasker's, as used to build Win50KGraph),
writes an uncompressed FASTA and runs `samtools faidx`. CHM13 keeps its NCBI
contig names. The result is
`references/chm13/GCF_009914755.1_T2T-CHM13v2.0_genomic.fna`.

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

## 5. Call SVs on CHM13, reusing the HG19 graph alignments

Calling a haplotype has two parts:

1. **Graph alignment** (`hotspot_search`, `query_align`, `query_summary`,
   `query_blocks`): the haplotype aligned to every local graph. This is the
   hours-long part, and it does not depend on the reference.
2. **Reference stages** (from `local_reference_templates` on): the
   alignments projected onto the chosen reference, then the VCF.

`--reuse-alignments DIR` takes part 1 from an earlier `individual` output made
with the same graph, here the HG19 example's `HG002_calls_hg19`. The files are
hard-linked from `DIR/samples/NAME/`, so it is instant on the same file
system. Check the plan first: the dry run should list `Would link …` lines
for both haplotypes:

```bash
python3 LinGraph/scripts/LinGraph.py individual \
  -I assemblies/prepared/query_paths.prepared.txt \
  -G Win50KGraph \
  -r references/chm13/GCF_009914755.1_T2T-CHM13v2.0_genomic.fna --reference-name CHM13_h1 \
  --reference-caches Win50KGraph/references/CHM13_h1_rig \
  --reuse-alignments HG002_calls_hg19 \
  -O HG002_calls_chm13 -t 64 --dry-run 2>&1 | grep -E "Would link|does not match|not found|ERROR"
```

Then run it:

```bash
python3 LinGraph/scripts/LinGraph.py individual \
  -I assemblies/prepared/query_paths.prepared.txt \
  -G Win50KGraph \
  -r references/chm13/GCF_009914755.1_T2T-CHM13v2.0_genomic.fna --reference-name CHM13_h1 \
  --reference-caches Win50KGraph/references/CHM13_h1_rig \
  --reuse-alignments HG002_calls_hg19 \
  -O HG002_calls_chm13 -t 64 > logs/individual_chm13.log 2>&1
```

- **Reference cache:** `--reference-caches` uses the CHM13 alignment shipped
  with Win50KGraph, so the reference is not aligned either.
- **Safeguards:** a stage is reused only if the earlier run used the same
  alignment settings; otherwise it is computed, and every later stage of that
  haplotype is computed too, never mixed with reused files. The dry run says
  `does not match` when that would happen.
- **Without the HG19 example:** drop `--reuse-alignments HG002_calls_hg19`;
  both haplotypes are then aligned to the graph (hours).

Each haplotype is called into
`HG002_calls_chm13/samples/HG002_hN/HG002_hN.vcf`, with a coverage report
`HG002_hN.coverage.summary.tsv` next to it. Calls are on CHM13's NCBI names
(`NC_060925.1` …), so they need no renaming; the truth set gets the same names
in step 6.

## 6. Benchmark against GIAB Q100 v5.0q

The benchmark needs truvari and bcftools. Install them in their own
environment (skip if you created it in the HG19 example):

```bash
conda create -y -n truvari -c conda-forge -c bioconda truvari bcftools htslib > logs/truvari_env.log 2>&1
conda activate truvari
```

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

`tools/BCFtoolMergeTwoHaplotypes.sh` turns the two haplotype VCFs into one
site-only VCF, `HG002.LinGraph.vcf.gz`, in the current folder: each is sorted
and stripped of genotypes, then they are combined with exact duplicates
removed. Combine the haplotypes and benchmark:

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

`benchmark/BenchGIAB.sh` runs the same benchmark and also reports it without
GIAB's most extreme truth variants.

## 7. Count SVs and check coverage

Count SVs of at least 50 bp per haplotype:

```bash
for h in h1 h2; do
  bash LinGraph/tools/count_sv.sh HG002_calls_chm13/samples/HG002_$h/HG002_$h.vcf
done
```

`total` counts insertions plus deletions; a call with both an inserted and a
deleted part of at least 50 bp counts twice (`total = INS_only + DEL_only +
2*both`).

Each haplotype's coverage report lists what the calls do not cover:

```bash
cat HG002_calls_chm13/samples/HG002_h1/HG002_h1.coverage.summary.tsv
```

`reference total_unannotated` is the reference sequence this haplotype has
no calls on; `query total_unannotated` is the assembly sequence placed on no
reference interval (N gaps, scaffold edges, sequence that did not align).
The `unmasked_bases` column leaves out soft-masked repeats, mostly
constitutive heterochromatin. Below, the reference side uses
`reference total_unannotated`, and the assembly side `query other`, which
also leaves out N gaps and scaffold edges.

## 8. Results

Results from running this tutorial:

| Haplotype | Precision vs GIAB | Recall vs GIAB | Unannotated assembly bases (excl. constitutive heterochromatin) | Unannotated reference bases (excl. constitutive heterochromatin) | SVs ≥ 50 bp |
| --- | --- | --- | --- | --- | --- |
| h1 (pat) | 0.986 | 0.966 | 4,303,423 | 126,791,810 | 14,230 |
| h2 (mat) | 0.986 | 0.966 | 4,207,899 | 52,815,189 | 14,727 |

- **Precision and recall** are for the diploid call set (h1 and h2
  combined), so both haplotypes share them: 2,906 Q100 SVs matched in
  non-difficult regions, 39 false positives, 103 false negatives (F1 0.976).
- **GIAB is not error-free.** On HG19, the older NIST SV v0.6 set is reported
  to contain about 5% errors, and agreement with it is lower (precision
  0.938, recall 0.978; see the [HG19 example](../HG002_HG19/README.md)). The
  Q100 benchmark on CHM13 is newer and assembly-based.
- **Unannotated bases** are `unmasked_bases` in the coverage reports, which
  leave out constitutive heterochromatin and other soft-masked repeats.
  - Assembly: the `query other` row (e.g. `query other 193 153097044
    4207899 ...` gives 4,207,899), which also leaves out N gaps and scaffold
    edges. About 4 Mb of each haplotype is placed nowhere on CHM13, against
    12–13 Mb on HG19.
  - Reference: the `reference total_unannotated` row. It includes whole
    chromosomes a haplotype lacks: h1 is paternal (chrY, no chrX) and h2
    maternal (chrX, no chrY), which accounts for most of the difference
    between them.
- **SV counts** are lower than on HG19, almost entirely in insertions (h1:
  6,136 insertion-only on CHM13 versus 9,347 on HG19; deletions are similar).
