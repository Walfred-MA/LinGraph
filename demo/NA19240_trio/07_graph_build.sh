#!/usr/bin/env bash
# Step 7: build a new trio pangenome graph on CHM13 with gene-based windows and
# the bundled alternative catalog, call every haplotype, merge on CHM13 and
# export the graph (cohort.gfa).
#   graph cache:  $DEMO/trio_graph    calls + GFA: $DEMO/trio_calls
# Graph stages run as SLURM jobs; keep this launcher running (nohup).
# Rerun the same command to resume. Run inside the LinGraph env.
set -euo pipefail
source "$(dirname "$0")/00_env.sh"

mkdir -p "$DEMO/graph_inputs"
cohort=$DEMO/graph_inputs/cohort.list
{
    printf 'CHM13_h1\t%s\n' "$CHM13"
    cat "$QUERIES"
} > "$cohort"
cat "$cohort"

cd "$LG"
python3 scripts/LinGraph.py graph \
    -I "$cohort" -G "$DEMO/trio_graph" -O "$DEMO/trio_calls" \
    -b windowprofs/geneblocks.bed --bed-grouped \
    -r CHM13_h1 --alternative data/alternatives.fa \
    --exact --mc-graph -t 32 \
    --slurm --slurm-jobs 50 \
    --slurm-account "$SLURM_ACCOUNT" --slurm-partition "$SLURM_PARTITION"

ls "$DEMO/trio_calls"/cohort.*.vcf "$DEMO/trio_calls/cohort.gfa" "$DEMO/trio_calls"/samples/*/*.vcf
