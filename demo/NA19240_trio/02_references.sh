#!/usr/bin/env bash
# Step 2: download CHM13 and GRCh38, keep the contigs the Win50KGraph caches
# were built with, soft-mask with WindowMasker and index.
#   CHM13: NCBI RefSeq GCF_009914755.1, 24 chromosomes (NC_060925.1 ...).
#   GRCh38: chr1-22, X, Y, M of the no-alt analysis set (25 contigs). LinGraph
#           calls on chr1-22, X and Y.
# Run inside the LinGraph env.
set -euo pipefail
source "$(dirname "$0")/00_env.sh"

RAW=$REFS/raw
mkdir -p "$RAW"
cd "$RAW"

NCBI=https://ftp.ncbi.nlm.nih.gov/genomes/all
fetch "$NCBI/GCF/009/914/755/GCF_009914755.1_T2T-CHM13v2.0/GCF_009914755.1_T2T-CHM13v2.0_genomic.fna.gz" \
    GCF_009914755.1_T2T-CHM13v2.0_genomic.fna.gz
fetch "$NCBI/GCA/000/001/405/GCA_000001405.15_GRCh38/seqs_for_alignment_pipelines.ucsc_ids/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna.gz" \
    GCA_000001405.15_GRCh38_no_alt_analysis_set.fna.gz

# GRCh38: keep the 25 main contigs (the analysis-set download is plain gzip).
if [ ! -s GRCh38_no_alt_main.fa.fai ]; then
    gzip -dc GCA_000001405.15_GRCh38_no_alt_analysis_set.fna.gz > GRCh38_no_alt.fna
    samtools faidx GRCh38_no_alt.fna
    samtools faidx GRCh38_no_alt.fna $(printf 'chr%s ' $(seq 1 22) X Y M) > GRCh38_no_alt_main.fa
    samtools faidx GRCh38_no_alt_main.fa
    rm -f GRCh38_no_alt.fna GRCh38_no_alt.fna.fai
fi

# Native reference names stay unchanged (no --contignamefix). --remask applies
# WindowMasker even though the downloads already carry their own soft-masking.
python3 "$LG/tools/prepare_assemblies.py" \
    -i "$RAW/GCF_009914755.1_T2T-CHM13v2.0_genomic.fna.gz" --name CHM13_h1 \
    -O "$REFS/chm13" --remask --threads 16
python3 "$LG/tools/prepare_assemblies.py" \
    -i "$RAW/GRCh38_no_alt_main.fa" --name HG38_h1 \
    -O "$REFS/hg38" --remask --threads 16

ls -l "$CHM13" "$CHM13.fai" "$HG38" "$HG38.fai"
