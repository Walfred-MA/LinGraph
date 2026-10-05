# Getting started: NA19240 trio (HGSVC3)

This demo runs LinGraph end to end on the YRI trio from HGSVC3:
NA19240 (child), NA19238 (mother) and NA19239 (father), two haplotypes each.

1. Install LinGraph.
2. Download and prepare CHM13 and GRCh38.
3. Download the Win50KGraph package and reconstruct its local graphs.
4. Download and prepare the six trio assemblies.
5. Call SVs per haplotype with Win50KGraph, on CHM13 and on GRCh38.
6. Benchmark trio consistency.
7. Build a new trio pangenome graph on CHM13 (gene windows plus the bundled
   alternative catalog), call every haplotype, merge on CHM13 and export GFA.

Each step is a script in this folder. Scripts are rerunnable: downloads and
long stages resume where they stopped.

## Setup

Pick an absolute folder with about **300 GB** free, clone LinGraph into it and
run the demo scripts from the clone:

```bash
export DEMO=/absolute/path/lingraph_demo
mkdir -p $DEMO && git clone https://github.com/Walfred-MA/LinGraph.git $DEMO/LinGraph
cd $DEMO/LinGraph/demo/NA19240_trio
source 00_env.sh          # repeat (with DEMO exported) in every new shell
```

Long steps are run with `nohup` on a compute node and log to `$LOGS`
(`$DEMO/logs`). `THREADS` (default 64) sets the CPU budget of steps 3 and 5.

## Steps

| Step | Command | Output |
| --- | --- | --- |
| 1 Install | `bash 01_install.sh` | `$DEMO/LinGraph`, conda envs `LinGraph` and `truvari` |
| 2 References | `bash 02_references.sh` | `$CHM13`, `$HG38` (+ `.fai`) |
| 3 Graph package | `bash 03_graph_package.sh` | `$DEMO/Win50KGraph` |
| 4 Assemblies | `bash 04_assemblies.sh` | `$ASM/prepared/query_paths.prepared.txt` |
| 5 SV calling | `bash 05_singular.sh CHM13` and `bash 05_singular.sh HG38` | `$DEMO/calls_win50k_chm13`, `$DEMO/calls_win50k_hg38` |
| 6 Benchmark | `bash 06_benchmark.sh $DEMO/calls_win50k_chm13` (likewise for the others) | `OUTPUT/benchmark/trio.*.tsv` |
| 7 Graph construction | `bash 07_graph_build.sh` | `$DEMO/trio_graph`, `$DEMO/trio_calls` (`cohort.*.vcf`, `cohort.gfa`) |

Steps 2–6 and 7 run inside the LinGraph environment (`conda activate
LinGraph`). Steps 2, 3 and 4 are independent of each other once step 1 is
done, except that step 3 needs `$CHM13` from step 2.

On the cluster, for example:

```bash
bash 01_install.sh > $LOGS/01.log 2>&1
conda activate LinGraph
nohup bash 02_references.sh > $LOGS/02.log 2>&1 &
nohup bash 04_assemblies.sh > $LOGS/04.log 2>&1 &
# after step 2 finishes:
nohup bash 03_graph_package.sh > $LOGS/03.log 2>&1 &
# after steps 3 and 4:
nohup bash 05_singular.sh CHM13 > $LOGS/05_chm13.log 2>&1 &
nohup bash 07_graph_build.sh > $LOGS/07.log 2>&1 &
```

## Notes

- **Reference choice.** Win50KGraph ships alignment caches for both CHM13
  (`CHM13_h1`) and GRCh38 (`HG38_h1`); use either with the same graph. The
  GRCh38 FASTA is the 25 main contigs of the no-alt analysis set, matching the
  shipped cache. LinGraph calls on chr1–22, X and Y. Keep one output folder per
  reference.
- **Preparation.** The HGSVC3 assemblies are unmasked and use contig names like
  `haplotype1-0000001`. Step 4 soft-masks them with WindowMasker and renames
  contigs to `NA19240#1#haplotype1-0000001`. References keep their native names.
- **Benchmarks.** Both benchmarks read the per-haplotype VCFs in
  `OUTPUT/samples/NAME/NAME.vcf` directly and restrict to the reference
  intervals covered by all four parental haplotypes. Child calls without a
  parental match are not confirmed errors; de novo variants and missed parental
  calls also count there.
- **Graph construction** submits its stages as SLURM jobs (`SLURM_ACCOUNT`,
  `SLURM_PARTITION` in `00_env.sh`); keep the launcher running.
