#!/usr/bin/env bash
# Step 3: download the Win50KGraph summary package (~2.4 GB) and reconstruct
# all local graphs from it with the CHM13 reference. Run inside the LinGraph env.
# Rerunning resumes the download and the reconstruction.
set -euo pipefail
source "$(dirname "$0")/00_env.sh"

cd "$DEMO"
fetch https://ndownloader.figshare.com/files/69451449 Win50KGraph.tar.gz
if [ ! -s "$PACKAGE/local_graphs.tsv" ]; then
    tar -xzf Win50KGraph.tar.gz          # -> cohort_minsetref_v3/summary/
fi

# Keep $PACKAGE in place afterwards: Win50KGraph/summary and
# Win50KGraph/references are links into it.
python3 "$LG/scripts/reconstruct_local_graph_folders.py" \
    -i "$PACKAGE" -r "$CHM13" --reference-haplotype CHM13_h1 \
    -o "$WIN50K" -j "$THREADS" --resume

ls -d "$WIN50K/references/CHM13_h1_rig" "$WIN50K/references/HG38_h1_rig"
