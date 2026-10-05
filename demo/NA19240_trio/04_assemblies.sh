#!/usr/bin/env bash
# Step 4: download the HGSVC3 Verkko assemblies of the YRI trio and prepare them
# (WindowMasker soft-masking, SAMPLE#HAP#contig names, faidx).
#   NA19240 child, NA19238 mother, NA19239 father; hap1/hap2 each.
# The downloads are unmasked and use contig names like haplotype1-0000001.
# Run inside the LinGraph env.
set -euo pipefail
source "$(dirname "$0")/00_env.sh"

BASE=https://ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/HGSVC3/working/20240201_verkko_batch3/assemblies
mkdir -p "$ASM/raw"
: > "$ASM/raw/raw_assemblies.list"
for sample in NA19240 NA19238 NA19239; do
    for hap in 1 2; do
        file=$sample.vrk-ps-sseq.asm-hap$hap.fasta.gz
        fetch "$BASE/$sample/$file" "$ASM/raw/$file"
        printf '%s_h%s\t%s\n' "$sample" "$hap" "$ASM/raw/$file" >> "$ASM/raw/raw_assemblies.list"
    done
done
cat "$ASM/raw/raw_assemblies.list"

# 6 assemblies at once, each masked with 8 parallel WindowMasker chunks.
python3 "$LG/tools/prepare_assemblies.py" \
    -q "$ASM/raw/raw_assemblies.list" -O "$ASM/prepared" \
    --contignamefix -j 6 --threads 8

cat "$QUERIES"
