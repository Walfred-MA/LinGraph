# NA19240 trio: building a pangenome graph

This tutorial builds a pangenome graph from the YRI trio assembled by HGSVC3:
NA19240 (child), NA19238 (mother) and NA19239 (father), two haplotypes each,
on CHM13 with gene-based windows and the bundled alternative sequences. It
calls every haplotype, merges the calls, exports the graph as GFA with one GAF
per haplotype, and then checks that the VCFs, the GFA and the GAFs are
lossless and that the graph's variants are consistent within the trio.

To call SVs with the precomputed Win50KGraph instead, see
[NA19240trio_SVcall](../NA19240trio_SVcall/README.md).

Run every command from one demo folder. Plan for about **300 GB** of disk.
Long commands write their progress to `logs/`; follow one with
`tail -f logs/NAME.log`. Steps marked *(hours)* belong on a compute node.

```bash
mkdir NA19240trio_graph && cd NA19240trio_graph
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

## 2. Prepare the CHM13 reference

```bash
mkdir references && cd references
curl -L -O https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/009/914/755/GCF_009914755.1_T2T-CHM13v2.0/GCF_009914755.1_T2T-CHM13v2.0_genomic.fna.gz
python ../LinGraph/tools/prepare_assemblies.py \
  -i GCF_009914755.1_T2T-CHM13v2.0_genomic.fna.gz --name CHM13_h1 \
  -O chm13 --remask --threads 16 > ../logs/prepare_chm13.log 2>&1
cd ..
```

## 3. Download and prepare the trio assemblies *(hours)*

```bash
mkdir assemblies && cd assemblies
URL=https://ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/HGSVC3/working/20240201_verkko_batch3/assemblies
for sample in NA19240 NA19238 NA19239; do
  for hap in 1 2; do
    curl -L -C - -O $URL/$sample/$sample.vrk-ps-sseq.asm-hap$hap.fasta.gz
  done
done
cat > raw_assemblies.list <<'EOF'
NA19240_h1 NA19240.vrk-ps-sseq.asm-hap1.fasta.gz
NA19240_h2 NA19240.vrk-ps-sseq.asm-hap2.fasta.gz
NA19238_h1 NA19238.vrk-ps-sseq.asm-hap1.fasta.gz
NA19238_h2 NA19238.vrk-ps-sseq.asm-hap2.fasta.gz
NA19239_h1 NA19239.vrk-ps-sseq.asm-hap1.fasta.gz
NA19239_h2 NA19239.vrk-ps-sseq.asm-hap2.fasta.gz
EOF
python ../LinGraph/tools/prepare_assemblies.py \
  -q raw_assemblies.list -O prepared --contignamefix -j 6 --threads 8 \
  > ../logs/prepare_assemblies.log 2>&1
cd ..
```

Preparation soft-masks each assembly and prefixes its contigs with the
haplotype (`NA19240#1#haplotype1-0000001`), so contig names stay unique
across the cohort and every call can be traced to its original contig. It
reads the compressed downloads directly and writes uncompressed, indexed
FASTAs listed in `assemblies/prepared/query_paths.prepared.txt`. For a cohort,
preparing once beforehand is recommended over LinGraph's lazy
`--force-prepare` mode, which copies every assembly again for each run.

## 4. Build the graph, call, merge and export *(hours)*

List CHM13 first, then the six prepared haplotypes:

```bash
printf 'CHM13_h1\t%s\n' "$PWD/references/chm13/GCF_009914755.1_T2T-CHM13v2.0_genomic.fna" > cohort.list
cat assemblies/prepared/query_paths.prepared.txt >> cohort.list

python LinGraph/scripts/LinGraph.py graph \
  -I cohort.list -G trio_graph -O trio_calls \
  -b LinGraph/windowprofs/geneblocks.bed --bed-grouped \
  -r CHM13_h1 --alternative LinGraph/data/alternatives.fa \
  --exact --make-graph -t 32 \
  --slurm --slurm-jobs 50 --slurm-account ACCOUNT --slurm-partition PARTITION \
  > logs/graph_trio.log 2>&1
```

Without SLURM, drop the last option line and set `-t` to the local CPU count.

- `-b LinGraph/windowprofs/geneblocks.bed --bed-grouped`: gene-based windows.
- `--alternative LinGraph/data/alternatives.fa`: the recommended alternative
  sequences for CHM13.
- `--exact`: merges the cohort calls, realigning merged SVs to the assemblies.
- `--make-graph`: exports the GFA and the per-haplotype GAFs.

`trio_graph/` holds the graph (repeat the command to resume). `trio_calls/`
holds:

- `samples/NAME/NAME.vcf`: each haplotype's calls;
- `cohort.sv.vcf`, `cohort.indel.vcf`, `cohort.snp.vcf`: the merged calls;
- `cohort.gfa`: the graph, with `cohort.gfa.variants.tsv` mapping each merged
  record to its path in the graph;
- `gaf/batch_NNN/NAME.gaf`: each haplotype's walk through the graph, and
  `gaf/added_links.gfa`: links the walks add (the final graph is `cohort.gfa`
  plus these links).

## 5. Benchmark

### Losslessness

The checks rebuild each haplotype's assembly from the calls and compare it
base by base, over every mapped region (`##pseudoLinearMapping`) of the
haplotype: first from its own VCF, then from the merged VCFs alone, then by
spelling its GAF records through the GFA.

```bash
O=trio_calls; C=trio_calls/lossless_check; T=LinGraph/tools; mkdir -p $C/vcf
REF=references/chm13/GCF_009914755.1_T2T-CHM13v2.0_genomic.fna
REFS=(-r $REF -r $O/checkpoints/local_reference_templates.fa -r $O/checkpoints/alternative_loci.fa)
VCFS=(); for d in $O/samples/*/; do n=$(basename $d); VCFS+=($d$n.vcf); done

# Per-haplotype VCFs
while read -r name fasta _; do
  python $T/check_vcf_lossless.py -v $O/samples/$name/$name.vcf "${REFS[@]}" -q $fasta \
    -o $C/vcf/$name.vcf_check.tsv > $C/vcf/$name.log 2>&1 &
done < $O/inputs/query_paths.normalized.txt; wait
grep -H -w lossless $C/vcf/*.log

# Merged VCFs, and losses the merge introduced
python $T/check_merge_lossless.py -v $O/cohort.sv.vcf $O/cohort.indel.vcf $O/cohort.snp.vcf \
  -s "${VCFS[@]}" "${REFS[@]}" -q $O/inputs/query_paths.normalized.txt -o $C/merge -t 32
python $T/compare_lossless_reports.py -m $C/merge -c $C/vcf

# GFA and GAFs
cat $O/cohort.gfa $O/gaf/added_links.gfa > $C/cohort.linked.gfa
python $T/check_gaf_lossless.py -g $C/cohort.linked.gfa -a $O/gaf/batch_*/*.gaf \
  -q $O/lingraph/mc-query_paths.list -s "${VCFS[@]}" -o $C/gaf -t 32
rm $C/cohort.linked.gfa
```

Results:

| Haplotype | Own VCF | Merged runs (exact / total) | Mapped bases checked | Merged VCF | Unified | GAF records (exact / total) | GFA/GAF |
| --- | --- | ---: | ---: | --- | --- | ---: | --- |
| CHM13_h1 | lossless | 1,383 / 1,383 | 2,894,884,633 | lossless | yes | 211 / 211 | lossless |
| NA19238_h1 | lossless | 1,432 / 1,432 | 2,850,971,586 | lossless | yes | 516 / 516 | lossless |
| NA19238_h2 | lossless | 1,453 / 1,453 | 2,850,126,679 | lossless | yes | 509 / 509 | lossless |
| NA19239_h1 | lossless | 1,523 / 1,523 | 2,730,526,958 | lossless | yes | 490 / 490 | lossless |
| NA19239_h2 | lossless | 1,518 / 1,518 | 2,854,901,024 | lossless | yes | 520 / 520 | lossless |
| NA19240_h1 | lossless | 1,480 / 1,480 | 2,801,177,944 | lossless | yes | 516 / 516 | lossless |
| NA19240_h2 | lossless | 1,539 / 1,539 | 2,852,542,304 | lossless | yes | 546 / 546 | lossless |

The merge introduced no loss for any haplotype. *Unified* means every carrier
of every merged record sits at that record's breakpoint with its own exact
allele. In the GFA, every GAF record spells its assembly exactly through
linked segments, and every mapped base of every haplotype lies in an exact
GAF record (the same mapped bases as above;
`trio_calls/lossless_check/gaf/summary.tsv`). Losslessness covers the anchored,
mapped assembly bases; unmapped contigs and unanchored insertions are
reported separately and are outside every check.

### Trio consistency of the graph's variants

Each merged record is one variant allele of the graph, carried by the
haplotypes on its line. `trio_line_check.py` checks every record a child
haplotype carries: is it carried by either parent?

```bash
python LinGraph/benchmark/trio_line_check.py \
  trio_calls/cohort.sv.vcf trio_calls/cohort.indel.vcf trio_calls/cohort.snp.vcf \
  --child NA19240 --mother NA19238 --father NA19239 \
  > trio_calls/trio_lines.tsv 2> trio_calls/trio_lines.summary.tsv
cat trio_calls/trio_lines.summary.tsv
```

| Variant | Child haplotype | Records | Found in parents | Child only | Total | Child only (%) |
| --- | --- | ---: | ---: | ---: | ---: | ---: |
| SV deletion | h1 | 23,896 | 8,772 | 162 | 8,934 | 1.81% |
| SV deletion | h2 | 23,896 | 9,001 | 172 | 9,173 | 1.88% |
| SV insertion | h1 | 25,768 | 8,824 | 214 | 9,038 | 2.37% |
| SV insertion | h2 | 25,768 | 9,172 | 245 | 9,417 | 2.60% |
| Small deletion | h1 | 1,028,849 | 392,478 | 8,085 | 400,563 | 2.02% |
| Small deletion | h2 | 1,028,849 | 399,156 | 8,107 | 407,263 | 1.99% |
| Small insertion | h1 | 992,929 | 379,000 | 8,301 | 387,301 | 2.14% |
| Small insertion | h2 | 992,929 | 385,161 | 8,072 | 393,233 | 2.05% |
| SNP | h1 | 7,318,593 | 3,309,027 | 17,184 | 3,326,211 | 0.52% |
| SNP | h2 | 7,318,593 | 3,359,346 | 12,557 | 3,371,903 | 0.37% |
| SV deletion (nested) | h1 | 835 | 197 | 22 | 219 | 10.05% |
| SV deletion (nested) | h2 | 835 | 184 | 9 | 193 | 4.66% |
| SV insertion (nested) | h1 | 152 | 37 | 5 | 42 | 11.90% |
| SV insertion (nested) | h2 | 152 | 40 | 1 | 41 | 2.44% |
| Small deletion (nested) | h1 | 19,546 | 4,107 | 909 | 5,016 | 18.12% |
| Small deletion (nested) | h2 | 19,546 | 3,991 | 755 | 4,746 | 15.91% |
| Small insertion (nested) | h1 | 4,427 | 1,080 | 169 | 1,249 | 13.53% |
| Small insertion (nested) | h2 | 4,427 | 1,008 | 70 | 1,078 | 6.49% |
| SNP (nested) | h1 | 33,392 | 10,120 | 307 | 10,427 | 2.94% |
| SNP (nested) | h2 | 33,392 | 9,807 | 141 | 9,948 | 1.42% |

SVs are at least 50 bp. *Total* counts the child-carried records with a
callable parent; *Child only* are those carried by neither parent.
The check is per record: a parent's allele merged into a neighbouring record
(for example a variant inside an insertion that was merged into another
parent record) counts as child-only, which is why nested records score lower.
HGSVC3's hap1/hap2 come from Strand-seq phasing and are not parent-of-origin
labels, so either parent is accepted for each child haplotype.
`trio_lines.tsv` lists each record with the six genotypes.
