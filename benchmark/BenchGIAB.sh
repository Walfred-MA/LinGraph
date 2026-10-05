#!/usr/bin/env bash
set -euo pipefail

# ================================================================
# BenchGIAB.sh
#
# Usage:
#
#   ./BenchGIAB.sh \
#       --ref GCF_009914755.1_T2T-CHM13v2.0_genomic.fa \
#       --truth HG002_CHM13v2.0_v5.0q_stvar.NC.vcf.gz \
#       --bed HG002_Q100.nondifficult.NC.bed \
#       calls/samples/HG002_h1/HG002_h1.vcf \
#       calls/samples/HG002_h2/HG002_h2.vcf \
#       Out/ > report.txt
#
# The BED ships compressed: tar -xJf HG002_Q100.nondifficult.NC.bed.tar.xz
# ================================================================

usage() {
    echo "Usage: $0 --ref REF.fa --truth TRUTH.vcf.gz --bed REGIONS.bed <hap1.vcf> <hap2.vcf> <outdir>" >&2
    exit 1
}

REF=""
TRUTH=""
BED=""
POSITIONAL=()
while [[ $# -gt 0 ]]; do
    case "$1" in
        --ref)   [[ $# -ge 2 ]] || usage; REF="$2";   shift 2 ;;
        --truth) [[ $# -ge 2 ]] || usage; TRUTH="$2"; shift 2 ;;
        --bed)   [[ $# -ge 2 ]] || usage; BED="$2";   shift 2 ;;
        -h|--help) usage ;;
        -*) echo "Unknown option: $1" >&2; usage ;;
        *) POSITIONAL+=("$1"); shift ;;
    esac
done

if [[ ${#POSITIONAL[@]} -ne 3 || -z "$REF" || -z "$TRUTH" || -z "$BED" ]]; then
    usage
fi

H1="${POSITIONAL[0]}"
H2="${POSITIONAL[1]}"
OUTDIR="${POSITIONAL[2]%/}"

mkdir -p "$OUTDIR"

H1_SITES="$OUTDIR/HG002_h1.sites.vcf.gz"
H2_SITES="$OUTDIR/HG002_h2.sites.vcf.gz"

MERGED="$OUTDIR/HG002.LinGraph.vcf.gz"

BENCH="$OUTDIR/truvari_HG002_nondifficult"


# ================================================================
# Check inputs
# ================================================================

for F in "$H1" "$H2" "$TRUTH" "$REF" "$BED"; do
    if [[ ! -e "$F" ]]; then
        echo "ERROR: missing file: $F" >&2
        exit 1
    fi
done


# ================================================================
# 1. Strip genotypes and sort haplotype VCFs
# ================================================================

echo "[1/5] Preparing haplotype VCFs..." >&2

bcftools sort "$H1" -Ou | \
    bcftools view -G -Oz -o "$H1_SITES"

tabix -f -p vcf "$H1_SITES"


bcftools sort "$H2" -Ou | \
    bcftools view -G -Oz -o "$H2_SITES"

tabix -f -p vcf "$H2_SITES"


# ================================================================
# 2. Merge haplotypes
# ================================================================

echo "[2/5] Merging haplotypes..." >&2

bcftools concat \
    -a \
    -d exact \
    "$H1_SITES" \
    "$H2_SITES" \
    -Ou | \
bcftools sort -Oz -o "$MERGED"

tabix -f -p vcf "$MERGED"


# ================================================================
# 3. Run Truvari
# ================================================================

echo "[3/5] Running Truvari..." >&2

rm -rf "$BENCH"

truvari bench \
    -b "$TRUTH" \
    -c "$MERGED" \
    -f "$REF" \
    -r 2000 \
    -C 5000 \
    --pick multi \
    --includebed "$BED" \
    -o "$BENCH" \
    >&2


# ================================================================
# Truvari versions use either tp-comp or older tp-call naming
# ================================================================

TP_BASE="$BENCH/tp-base.vcf.gz"

if [[ -f "$BENCH/tp-comp.vcf.gz" ]]; then
    TP_CALL="$BENCH/tp-comp.vcf.gz"
elif [[ -f "$BENCH/tp-call.vcf.gz" ]]; then
    TP_CALL="$BENCH/tp-call.vcf.gz"
else
    echo "ERROR: cannot find tp-comp.vcf.gz or tp-call.vcf.gz" >&2
    exit 1
fi

FN="$BENCH/fn.vcf.gz"
FP="$BENCH/fp.vcf.gz"


# ================================================================
# 4. Remove extreme GIAB variants
#
# IMPORTANT:
#
# REMAP and LCR come from GIAB truth.
#
# Therefore:
#
#   tp-base  -> has REMAP/LCR
#   fn       -> has REMAP/LCR
#
#   tp-comp  -> LinGraph, does NOT have REMAP/LCR
#   fp       -> LinGraph, does NOT have REMAP/LCR
#
# For tp-comp we propagate exclusions through MatchId.
# FP is left completely unchanged.
# ================================================================

echo "[4/5] Removing extreme GIAB variants..." >&2

FILTER='INFO/REMAP="partial" || INFO/LCR>=0.99'


TP_BASE_NOEXT="$BENCH/tp-base.no_extreme.vcf.gz"
TP_BASE_EXTREME="$BENCH/tp-base.extreme.vcf.gz"

FN_NOEXT="$BENCH/fn.no_extreme.vcf.gz"

TP_CALL_NOEXT="$BENCH/tp-comp.no_extreme.vcf.gz"

FP_NOEXT="$BENCH/fp.no_extreme.vcf.gz"

EXCLUDED_MATCHIDS="$BENCH/extreme.matchids.txt"


# ----------------------------------------------------------------
# Filter TP-base
# ----------------------------------------------------------------

bcftools view \
    -e "$FILTER" \
    "$TP_BASE" \
    -Oz -o "$TP_BASE_NOEXT"

tabix -f -p vcf "$TP_BASE_NOEXT"


# Save the excluded TP-base variants as well
bcftools view \
    -i "$FILTER" \
    "$TP_BASE" \
    -Oz -o "$TP_BASE_EXTREME"

tabix -f -p vcf "$TP_BASE_EXTREME"


# ----------------------------------------------------------------
# Filter FN
# ----------------------------------------------------------------

bcftools view \
    -e "$FILTER" \
    "$FN" \
    -Oz -o "$FN_NOEXT"

tabix -f -p vcf "$FN_NOEXT"


# ----------------------------------------------------------------
# Extract MatchIds of excluded truth TPs
#
# Truvari places MatchId in both sides of a matched pair.
# ----------------------------------------------------------------

bcftools query \
    -f '%INFO/MatchId\n' \
    "$TP_BASE_EXTREME" | \
awk '$0 != "." && $0 != ""' | \
sort -u > "$EXCLUDED_MATCHIDS"


# ----------------------------------------------------------------
# Remove corresponding tp-comp / tp-call records
# ----------------------------------------------------------------

if [[ -s "$EXCLUDED_MATCHIDS" ]]; then

    {
        bcftools view -h "$TP_CALL"

        bcftools view -H "$TP_CALL" | \
        awk -v idfile="$EXCLUDED_MATCHIDS" '
        BEGIN {
            FS=OFS="\t"

            while ((getline line < idfile) > 0) {
                drop[line] = 1
            }

            close(idfile)
        }

        {
            matchid=""

            n=split($8, info, ";")

            for (i=1; i<=n; i++) {
                if (info[i] ~ /^MatchId=/) {
                    matchid=substr(info[i], 9)
                    break
                }
            }

            if (!(matchid in drop))
                print
        }
        '

    } | bcftools view -Oz -o "$TP_CALL_NOEXT"

else

    # No extreme TP-base variants
    bcftools view \
        "$TP_CALL" \
        -Oz -o "$TP_CALL_NOEXT"

fi

tabix -f -p vcf "$TP_CALL_NOEXT"


# ----------------------------------------------------------------
# FP has no GIAB REMAP/LCR annotation.
#
# Do NOT filter it.
# Make a copy only so naming is consistent.
# ----------------------------------------------------------------

bcftools view \
    "$FP" \
    -Oz -o "$FP_NOEXT"

tabix -f -p vcf "$FP_NOEXT"


# ================================================================
# 5. Count and calculate benchmark
# ================================================================

echo "[5/5] Calculating metrics..." >&2


count_vcf()
{
    bcftools view -H "$1" | \
        awk 'END {print NR+0}'
}


calc_precision()
{
    awk \
        -v tp="$1" \
        -v fp="$2" \
        'BEGIN {
            if (tp+fp == 0)
                print "NA"
            else
                printf "%.6f", tp/(tp+fp)
        }'
}


calc_recall()
{
    awk \
        -v tp="$1" \
        -v fn="$2" \
        'BEGIN {
            if (tp+fn == 0)
                print "NA"
            else
                printf "%.6f", tp/(tp+fn)
        }'
}


calc_f1()
{
    awk \
        -v p="$1" \
        -v r="$2" \
        'BEGIN {
            if (p == "NA" || r == "NA" || p+r == 0)
                print "NA"
            else
                printf "%.6f", 2*p*r/(p+r)
        }'
}


# ================================================================
# Raw counts
# ================================================================

TP_BASE_RAW=$(count_vcf "$TP_BASE")
TP_CALL_RAW=$(count_vcf "$TP_CALL")

FN_RAW=$(count_vcf "$FN")
FP_RAW=$(count_vcf "$FP")


# ================================================================
# Filtered counts
# ================================================================

TP_BASE_FILTERED=$(count_vcf "$TP_BASE_NOEXT")
TP_CALL_FILTERED=$(count_vcf "$TP_CALL_NOEXT")

FN_FILTERED=$(count_vcf "$FN_NOEXT")

# FP is intentionally unchanged
FP_FILTERED="$FP_RAW"


TP_BASE_REMOVED=$((TP_BASE_RAW - TP_BASE_FILTERED))
TP_CALL_REMOVED=$((TP_CALL_RAW - TP_CALL_FILTERED))

FN_REMOVED=$((FN_RAW - FN_FILTERED))

# By definition here:
FP_REMOVED=0


# ================================================================
# Raw Truvari metrics
# ================================================================

PRECISION_RAW=$(calc_precision \
    "$TP_CALL_RAW" \
    "$FP_RAW")

RECALL_RAW=$(calc_recall \
    "$TP_BASE_RAW" \
    "$FN_RAW")

F1_RAW=$(calc_f1 \
    "$PRECISION_RAW" \
    "$RECALL_RAW")


# ================================================================
# Recalculated metrics after removing extreme GIAB truth variants
#
# precision:
#
#   filtered call-side TP / (filtered call-side TP + original FP)
#
# recall:
#
#   filtered truth-side TP / (filtered truth-side TP + filtered FN)
# ================================================================

PRECISION=$(calc_precision \
    "$TP_CALL_FILTERED" \
    "$FP_FILTERED")

RECALL=$(calc_recall \
    "$TP_BASE_FILTERED" \
    "$FN_FILTERED")

F1=$(calc_f1 \
    "$PRECISION" \
    "$RECALL")


# ================================================================
# Report
#
# Everything below goes to stdout -> report.txt
# ================================================================

echo "HG002 LinGraph vs GIAB"
echo "======================"
echo

echo "Input:"
echo "  haplotype 1:       $H1"
echo "  haplotype 2:       $H2"
echo "  merged VCF:        $MERGED"
echo

echo "Benchmark:"
echo "  GIAB truth:        $TRUTH"
echo "  GIAB regions:      $BED"
echo "  Truvari refdist:   2000"
echo

echo "Extreme GIAB filter:"
echo '  REMAP = partial'
echo "  OR"
echo "  LCR >= 0.99"
echo

echo "Raw Truvari counts:"
echo "  TP-base:           $TP_BASE_RAW"
echo "  TP-call:           $TP_CALL_RAW"
echo "  FN:                $FN_RAW"
echo "  FP:                $FP_RAW"
echo

echo "Removed as extreme:"
echo "  TP-base removed:   $TP_BASE_REMOVED"
echo "  TP-call removed:   $TP_CALL_REMOVED"
echo "  FN removed:        $FN_REMOVED"
echo "  FP removed:        $FP_REMOVED"
echo

echo "No-extreme counts:"
echo "  TP-base:           $TP_BASE_FILTERED"
echo "  TP-call:           $TP_CALL_FILTERED"
echo "  FN:                $FN_FILTERED"
echo "  FP:                $FP_FILTERED"
echo

echo "Raw benchmark:"
echo "  Precision:         $PRECISION_RAW"
echo "  Recall:            $RECALL_RAW"
echo "  F1:                $F1_RAW"
echo

echo "No-extreme benchmark:"
echo "  Precision:         $PRECISION"
echo "  Recall:            $RECALL"
echo "  F1:                $F1"
