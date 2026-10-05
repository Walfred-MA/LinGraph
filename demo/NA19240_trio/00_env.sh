# Source in every shell before running a step:
#   export DEMO=/absolute/path/for/the/demo
#   source demo/NA19240_trio/00_env.sh
# Every path below lives under $DEMO.

if [ -z "${DEMO:-}" ]; then
    echo "Set DEMO to an absolute folder first, e.g. export DEMO=\$PWD/lingraph_demo" >&2
    return 1 2>/dev/null || exit 1
fi
case "$DEMO" in
    /*) ;;
    *) echo "DEMO must be an absolute path: $DEMO" >&2; return 1 2>/dev/null || exit 1 ;;
esac

export LG=$DEMO/LinGraph                      # GitHub checkout
export LOGS=$DEMO/logs
export REFS=$DEMO/references
export ASM=$DEMO/assemblies
export PACKAGE=$DEMO/cohort_minsetref_v3/summary
export WIN50K=$DEMO/Win50KGraph               # reconstructed graph cache

# Prepared (WindowMasker-masked, indexed) references. Names and lengths match
# the reference caches shipped in Win50KGraph.
export CHM13=$REFS/chm13/GCF_009914755.1_T2T-CHM13v2.0_genomic.fna
export HG38=$REFS/hg38/GRCh38_no_alt_main.fa

export QUERIES=$ASM/prepared/query_paths.prepared.txt

# CPU budget per step; adjust to the node you run on.
export THREADS=${THREADS:-64}

# SLURM settings for the graph-construction step (step 7).
export SLURM_ACCOUNT=${SLURM_ACCOUNT:-mchaisso_100}
export SLURM_PARTITION=${SLURM_PARTITION:-qcb}

# fetch URL FILE: resume a partial download, skip a finished one. -f keeps HTTP
# error pages out of FILE; FILE.done marks a completed download.
fetch() {
    if [ -e "$2.done" ]; then
        echo "already downloaded: $2"
        return 0
    fi
    curl -fL --retry 5 -C - -o "$2" "$1" && touch "$2.done"
}

mkdir -p "$LOGS"
if [ -d "$LG/scripts" ]; then
    export PATH="$LG/scripts:$PATH"
fi
