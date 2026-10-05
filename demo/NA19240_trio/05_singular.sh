#!/usr/bin/env bash
# Step 5: call SVs for all six haplotypes against the Win50KGraph cache.
#   usage: 05_singular.sh CHM13   -> $DEMO/calls_win50k_chm13
#          05_singular.sh HG38    -> $DEMO/calls_win50k_hg38 (chr1-22, X, Y)
# Each haplotype is called independently into OUTPUT/samples/NAME/NAME.vcf;
# --merge also writes the merged trio VCFs. Run inside the LinGraph env.
set -euo pipefail
source "$(dirname "$0")/00_env.sh"

case "${1:-CHM13}" in
    CHM13) fasta=$CHM13; name=CHM13_h1; out=$DEMO/calls_win50k_chm13 ;;
    HG38)  fasta=$HG38;  name=HG38_h1;  out=$DEMO/calls_win50k_hg38 ;;
    *) echo "usage: $0 CHM13|HG38" >&2; exit 2 ;;
esac

cd "$LG"
# The shipped reference cache is reused directly (no reference re-alignment).
python3 scripts/LinGraph.py singular \
    -I "$QUERIES" -G "$WIN50K" \
    -r "$fasta" --reference-name "$name" \
    --reference-caches "$WIN50K/references/${name}_rig" \
    -O "$out" -t "$THREADS" --merge

ls "$out"/samples/*/*.vcf
