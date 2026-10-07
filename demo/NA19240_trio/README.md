# Getting started: NA19240 trio

This tutorial runs LinGraph on the YRI trio assembled by HGSVC3: NA19240
(child), NA19238 (mother) and NA19239 (father), two haplotypes each. It shows
two ways to use LinGraph:

- **SV calling with a precomputed graph** (`singular`): call each haplotype
  against Win50KGraph, on CHM13 or GRCh38, then check trio consistency.
- **Building a new pangenome graph** (`graph`): build a trio graph on CHM13
  with gene-based windows and the bundled alternative sequences, call every
  haplotype, merge the calls on CHM13 and export the graph as GFA.

Run every command from one demo folder. Plan for about **300 GB** of disk.
The tools print detailed progress, so each long command below writes it to a
file in `logs/`; follow one with `tail -f logs/NAME.log`. Steps that take hours
are marked; run them on a compute node, for example by prefixing the command
with `nohup` and appending `&`.

```bash
mkdir lingraph_demo && cd lingraph_demo
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

LinGraph needs references as uncompressed, indexed FASTAs. Use the same
sequences as the Win50KGraph reference caches:

- **CHM13**: NCBI RefSeq T2T-CHM13v2.0, 24 chromosomes (`NC_060925.1` ...).
- **GRCh38**: the 25 main contigs (chr1–22, X, Y, M) of the no-alt analysis
  set. LinGraph calls on chr1–22, X and Y.

```bash
mkdir references && cd references

curl -L -O https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/009/914/755/GCF_009914755.1_T2T-CHM13v2.0/GCF_009914755.1_T2T-CHM13v2.0_genomic.fna.gz
python3 ../LinGraph/tools/prepare_assemblies.py \
  -i GCF_009914755.1_T2T-CHM13v2.0_genomic.fna.gz --name CHM13_h1 \
  -O chm13 --remask --threads 16 > ../logs/prepare_chm13.log 2>&1

curl -L -O https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/000/001/405/GCA_000001405.15_GRCh38/seqs_for_alignment_pipelines.ucsc_ids/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna.gz
gzip -dc GCA_000001405.15_GRCh38_no_alt_analysis_set.fna.gz > GRCh38_no_alt.fna
samtools faidx GRCh38_no_alt.fna
samtools faidx GRCh38_no_alt.fna chr{1..22} chrX chrY chrM > GRCh38_main.fa
python3 ../LinGraph/tools/prepare_assemblies.py \
  -i GRCh38_main.fa --name HG38_h1 -O hg38 --remask --threads 16 > ../logs/prepare_hg38.log 2>&1

cd ..
```

`prepare_assemblies.py` soft-masks with WindowMasker, writes an uncompressed
FASTA and runs `samtools faidx`. Both downloads are already soft-masked;
`--remask` replaces that masking with WindowMasker's, as used to build
Win50KGraph, because alignment seeding and several calling thresholds use the
soft-masking. Reference contigs keep their native names. The results are
`references/chm13/GCF_009914755.1_T2T-CHM13v2.0_genomic.fna` and
`references/hg38/GRCh38_main.fa`.

## 3. Download and reconstruct Win50KGraph

Win50KGraph is a precomputed graph summary (balanced 50 kb windows) for
calling individual samples.

```bash
curl -L -C - https://ndownloader.figshare.com/files/69451449 -o Win50KGraph.tar.gz
tar -xzf Win50KGraph.tar.gz

python3 LinGraph/scripts/reconstruct_local_graph_folders.py \
  -i cohort_minsetref_v3/summary \
  -r references/chm13/GCF_009914755.1_T2T-CHM13v2.0_genomic.fna \
  --reference-haplotype CHM13_h1 \
  -o Win50KGraph -j 32 --resume > logs/reconstruct.log 2>&1
```

Reconstruction (hours) rebuilds every local graph into `Win50KGraph/`. Its
`summary/` and `references/` entries link back to the extracted
`cohort_minsetref_v3/summary`, so keep that folder. Repeat the command to
resume an interrupted reconstruction.

## 4. Prepare the trio assemblies

Download both haplotypes of each sample from HGSVC3:

```bash
mkdir assemblies && cd assemblies
URL=https://ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/HGSVC3/working/20240201_verkko_batch3/assemblies
for sample in NA19240 NA19238 NA19239; do
  for hap in 1 2; do
    curl -L -C - -O $URL/$sample/$sample.vrk-ps-sseq.asm-hap$hap.fasta.gz
  done
done
```

List one assembly per haplotype, named `SAMPLE_hN`:

```bash
cat > raw_assemblies.list <<'EOF'
NA19240_h1 NA19240.vrk-ps-sseq.asm-hap1.fasta.gz
NA19240_h2 NA19240.vrk-ps-sseq.asm-hap2.fasta.gz
NA19238_h1 NA19238.vrk-ps-sseq.asm-hap1.fasta.gz
NA19238_h2 NA19238.vrk-ps-sseq.asm-hap2.fasta.gz
NA19239_h1 NA19239.vrk-ps-sseq.asm-hap1.fasta.gz
NA19239_h2 NA19239.vrk-ps-sseq.asm-hap2.fasta.gz
EOF

python3 ../LinGraph/tools/prepare_assemblies.py \
  -q raw_assemblies.list -O prepared --contignamefix -j 6 --threads 8 \
  > ../logs/prepare_assemblies.log 2>&1
cd ..
```

These assemblies are unmasked and name contigs like `haplotype1-0000001`.
Preparation (hours) masks them and renames each contig with its haplotype
prefix, e.g. `NA19240#1#haplotype1-0000001`. `-j 6` prepares all six at once;
`--threads 8` masks each one in eight parallel chunks. The calling input list
is `assemblies/prepared/query_paths.prepared.txt`.

## 5. Call SVs with Win50KGraph

Call all six haplotypes against CHM13 (hours):

```bash
python3 LinGraph/scripts/LinGraph.py singular \
  -I assemblies/prepared/query_paths.prepared.txt \
  -G Win50KGraph \
  -r references/chm13/GCF_009914755.1_T2T-CHM13v2.0_genomic.fna --reference-name CHM13_h1 \
  --reference-caches Win50KGraph/references/CHM13_h1_rig \
  -O calls_chm13 -t 64 --merge > logs/singular_chm13.log 2>&1
```

The same graph also calls on GRCh38; only the reference and its cache change.
Use a separate output folder:

```bash
python3 LinGraph/scripts/LinGraph.py singular \
  -I assemblies/prepared/query_paths.prepared.txt \
  -G Win50KGraph \
  -r references/hg38/GRCh38_main.fa --reference-name HG38_h1 \
  --reference-caches Win50KGraph/references/HG38_h1_rig \
  -O calls_hg38 -t 64 --merge > logs/singular_hg38.log 2>&1
```

Each haplotype is called independently into `calls_chm13/samples/NAME/NAME.vcf`
with coverage reports; `--merge` also writes the merged trio VCFs
(`calls_chm13/cohort.*.vcf`). `--reference-caches` reuses the alignment of the
reference that ships with the package instead of recomputing it.

## 6. Check trio consistency

Both checks compare each child haplotype with the four parental haplotypes,
inside the reference intervals covered by all parents.

**QuickTriocheck** (Python only):

```bash
python3 LinGraph/benchmark/QuickTriocheck.py \
  --child  calls_chm13/samples/NA19240_h1/NA19240_h1.vcf calls_chm13/samples/NA19240_h2/NA19240_h2.vcf \
  --mother calls_chm13/samples/NA19238_h1/NA19238_h1.vcf calls_chm13/samples/NA19238_h2/NA19238_h2.vcf \
  --father calls_chm13/samples/NA19239_h1/NA19239_h1.vcf calls_chm13/samples/NA19239_h2/NA19239_h2.vcf \
  -q assemblies/prepared/query_paths.prepared.txt \
  --fp calls_chm13/trio.quick.fp.vcf \
  > calls_chm13/trio.quick.tsv 2> logs/quicktriocheck_chm13.log
```

`-q` supplies the child contig lengths used to skip calls near assembly-contig
edges.

**Truvari** (eight pairwise benchmarks; needs its own environment):

```bash
conda create -y -n truvari -c conda-forge -c bioconda truvari htslib > logs/truvari_env.log 2>&1
conda activate truvari
bash LinGraph/benchmark/truvari_trio.sh \
  --child  calls_chm13/samples/NA19240_h1/NA19240_h1.vcf calls_chm13/samples/NA19240_h2/NA19240_h2.vcf \
  --mother calls_chm13/samples/NA19238_h1/NA19238_h1.vcf calls_chm13/samples/NA19238_h2/NA19238_h2.vcf \
  --father calls_chm13/samples/NA19239_h1/NA19239_h1.vcf calls_chm13/samples/NA19239_h2/NA19239_h2.vcf \
  --distance 500 --jobs 8 --fp calls_chm13/trio.truvari.fp.vcf \
  > calls_chm13/trio.truvari.tsv 2> logs/truvari_chm13.log
conda activate LinGraph
```

Replace `calls_chm13` with `calls_hg38` (or `trio_calls` from step 7) to check
the other call sets. A child call without a parental match is not necessarily
an error: de novo variants and calls missed in a parent also land there.

## 7. Build a trio pangenome graph on CHM13

List CHM13 first, then the six prepared haplotypes:

```bash
printf 'CHM13_h1\t%s\n' "$PWD/references/chm13/GCF_009914755.1_T2T-CHM13v2.0_genomic.fna" > cohort.list
cat assemblies/prepared/query_paths.prepared.txt >> cohort.list
```

Build the graph, call every haplotype, merge on CHM13 and export the GFA. On a
SLURM cluster, graph stages run as jobs; keep this launcher running:

```bash
python3 LinGraph/scripts/LinGraph.py graph \
  -I cohort.list -G trio_graph -O trio_calls \
  -b LinGraph/windowprofs/geneblocks.bed --bed-grouped \
  -r CHM13_h1 --alternative LinGraph/data/alternatives.fa \
  --exact --make-graph -t 32 \
  --slurm --slurm-jobs 50 --slurm-account YOUR_ACCOUNT --slurm-partition YOUR_PARTITION \
  > logs/graph_trio.log 2>&1
```

Without SLURM, drop the last line and set `-t` to the local CPU count.

- `-b LinGraph/windowprofs/geneblocks.bed --bed-grouped`: gene-based windows.
- `--alternative LinGraph/data/alternatives.fa`: the recommended alternative
  sequences for CHM13, imported into the graph.
- `--exact`: merges the cohort calls and realigns merged SVs to the assemblies.
- `--make-graph`: exports `trio_calls/cohort.gfa`.

`trio_graph/` holds the graph cache (resume by repeating the command);
`trio_calls/` holds per-haplotype VCFs in `samples/`, the merged
`cohort.sv.vcf`, `cohort.indel.vcf`, `cohort.snp.vcf` and `cohort.gfa`. Check
these calls with step 6, using `trio_calls` as the folder.
