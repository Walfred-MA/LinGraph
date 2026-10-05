#!/usr/bin/env bash
# Step 6: trio consistency of the per-haplotype calls in one output folder.
#   usage: 06_benchmark.sh OUTPUT_DIR      (e.g. $DEMO/calls_win50k_chm13)
# Writes OUTPUT_DIR/benchmark/trio.quick.tsv and, when the truvari env exists,
# OUTPUT_DIR/benchmark/trio.truvari.tsv. Both read LinGraph's per-haplotype
# VCFs directly, using their ##referenceCoverage headers.
set -euo pipefail
source "$(dirname "$0")/00_env.sh"

out=$(cd "${1:?usage: $0 OUTPUT_DIR}" && pwd)
bench=$out/benchmark
mkdir -p "$bench"
vcf() { echo "$out/samples/$1/$1.vcf"; }
CHILD=("$(vcf NA19240_h1)" "$(vcf NA19240_h2)")
MOTHER=("$(vcf NA19238_h1)" "$(vcf NA19238_h2)")
FATHER=("$(vcf NA19239_h1)" "$(vcf NA19239_h2)")
ls -l "${CHILD[@]}" "${MOTHER[@]}" "${FATHER[@]}"

# QuickTriocheck: SVs >= 50 bp matched by position, type and size. -q gives the
# child contig lengths used to trim 10 kb assembly-contig edges.
python3 "$LG/benchmark/QuickTriocheck.py" \
    --child "${CHILD[@]}" --mother "${MOTHER[@]}" --father "${FATHER[@]}" \
    -q "$QUERIES" --fp "$bench/quick.fp.vcf" --tp "$bench/quick.tp.vcf" \
    > "$bench/trio.quick.tsv" 2> "$bench/trio.quick.log"
echo "QuickTriocheck: $bench/trio.quick.tsv"

# Truvari: eight child-vs-parent benchmarks, 500 bp breakpoint tolerance.
TRUVARI_BIN=$(conda info --base)/envs/truvari/bin
if [ -x "$TRUVARI_BIN/truvari" ]; then
    (cd "$bench" && PATH="$TRUVARI_BIN:$PATH" bash "$LG/benchmark/truvari_trio.sh" \
        --child "${CHILD[@]}" --mother "${MOTHER[@]}" --father "${FATHER[@]}" \
        --distance 500 --jobs 8 --fp "$bench/truvari.fp.vcf" \
        > "$bench/trio.truvari.tsv" 2> "$bench/trio.truvari.log")
    echo "Truvari: $bench/trio.truvari.tsv"
else
    echo "truvari env not found ($TRUVARI_BIN); skipped the Truvari benchmark" >&2
fi
