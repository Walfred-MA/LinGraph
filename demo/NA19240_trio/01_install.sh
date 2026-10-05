#!/usr/bin/env bash
# Step 1: get the LinGraph code and build the LinGraph conda environment.
# Needs conda (mamba is used automatically when present), make and a C++17 compiler.
set -euo pipefail
source "$(dirname "$0")/00_env.sh"

mkdir -p "$DEMO"
if [ -d "$LG/.git" ]; then
    git -C "$LG" pull --ff-only
else
    git clone https://github.com/Walfred-MA/LinGraph.git "$LG"
fi

cd "$LG"
# Creates/updates the LinGraph env, installs dependencies, compiles the C++ tools.
python3 install.py --conda-env LinGraph -j 16
echo "Installed. In every new shell: conda activate LinGraph"

# Separate environment for the optional Truvari benchmark (step 6).
if ! conda env list | grep -qE '^truvari\s'; then
    if command -v mamba >/dev/null; then
        mamba create -y -n truvari -c conda-forge -c bioconda truvari htslib
    else
        conda create -y -n truvari -c conda-forge -c bioconda truvari htslib
    fi
fi
