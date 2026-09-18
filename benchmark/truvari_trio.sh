#!/usr/bin/env bash
set -euo pipefail

# truvari_trio.sh (parent-coverage aware)
#
# Haplotype-resolved trio benchmarking with Truvari.
#
# - Keeps child/mother/father haplotypes separate.
# - Runs 8 independent pairwise Truvari benchmarks:
#     child_h1 vs mother_h1, mother_h2, father_h1, father_h2
#     child_h2 vs mother_h1, mother_h2, father_h1, father_h2
# - Runs pairwise benchmarks in parallel (--jobs, default 8).
# - Uses Truvari defaults except --refdist, controlled by --distance.
# - Restricts every benchmark to the intersection of the ##referenceCoverage
#   intervals from all four parental haplotype VCFs via Truvari --includebed.
# - Uses truvari consistency on child-side fn.vcf.gz files to identify
#   child SVs absent from both maternal and both paternal haplotypes.
# - Reports h1 and h2 errors separately and combined.
# - Writes separate h1/h2 FP VCFs plus a combined sites-only FP VCF.
# - Does not require the external bcftools executable.
# - Repairs missing ##contig header entries and sorts/bgzip/tabix indexes inputs.

usage() {
    cat >&2 <<'USAGE'
Usage:
  bash benchmark/truvari_trio.sh \
    --child  child_h1.vcf child_h2.vcf \
    --mother mother_h1.vcf mother_h2.vcf \
    --father father_h1.vcf father_h2.vcf \
    [--distance 500] \
    [--jobs 8] \
    [--fp child.fp.vcf] \
    [--keep-temp]

Example:
  bash benchmark/truvari_trio.sh \
    --child \
      cohort_calls/samples/NA19240_h1/NA19240_h1.vcf \
      cohort_calls/samples/NA19240_h2/NA19240_h2.vcf \
    --mother \
      cohort_calls/samples/NA19238_h1/NA19238_h1.vcf \
      cohort_calls/samples/NA19238_h2/NA19238_h2.vcf \
    --father \
      cohort_calls/samples/NA19239_h1/NA19239_h1.vcf \
      cohort_calls/samples/NA19239_h2/NA19239_h2.vcf \
    --distance 0 \
    --jobs 8 \
    --fp NA19240.fp1.vcf \
    > NA19240.triocheck01.txt \
    2> NA19240.triocheck01.err &

With --fp NA19240.fp1.vcf, the script writes:
  NA19240.fp1.h1.vcf   child h1 SVs absent from all four parental haplotypes
  NA19240.fp1.h2.vcf   child h2 SVs absent from all four parental haplotypes
  NA19240.fp1.vcf      combined sites-only FP VCF with INFO/TRIO_HAP=h1|h2

The script always reads the 0-based, half-open ##referenceCoverage records
from mother_h1, mother_h2, father_h1, and father_h2, intersects their covered
intervals, and passes the resulting BED to every truvari bench command. Use
--keep-temp to retain the generated BED.
USAGE
    exit 1
}

CHILD=()
MOTHER=()
FATHER=()
DIST=500
JOBS=8
FP=""
KEEP_TEMP=0
INCLUDE_BED_PATH=""

while [[ $# -gt 0 ]]; do
    case "$1" in
        --child)
            [[ $# -ge 3 ]] || usage
            CHILD=("$2" "$3")
            shift 3
            ;;
        --mother)
            [[ $# -ge 3 ]] || usage
            MOTHER=("$2" "$3")
            shift 3
            ;;
        --father)
            [[ $# -ge 3 ]] || usage
            FATHER=("$2" "$3")
            shift 3
            ;;
        --distance)
            [[ $# -ge 2 ]] || usage
            DIST="$2"
            shift 2
            ;;
        --jobs|-j)
            [[ $# -ge 2 ]] || usage
            JOBS="$2"
            shift 2
            ;;
        --fp)
            [[ $# -ge 2 ]] || usage
            FP="$2"
            shift 2
            ;;
        --keep-temp)
            KEEP_TEMP=1
            shift
            ;;
        -h|--help)
            usage
            ;;
        *)
            echo "ERROR: unknown option: $1" >&2
            usage
            ;;
    esac
done

[[ ${#CHILD[@]} -eq 2 ]] || usage
[[ ${#MOTHER[@]} -eq 2 ]] || usage
[[ ${#FATHER[@]} -eq 2 ]] || usage

if ! [[ "$JOBS" =~ ^[1-9][0-9]*$ ]]; then
    echo "ERROR: --jobs must be a positive integer" >&2
    exit 1
fi

if ! [[ "$DIST" =~ ^[0-9]+$ ]]; then
    echo "ERROR: --distance must be a non-negative integer" >&2
    exit 1
fi

for f in "${CHILD[@]}" "${MOTHER[@]}" "${FATHER[@]}"; do
    if [[ ! -r "$f" ]]; then
        echo "ERROR: cannot read input VCF: $f" >&2
        exit 1
    fi
done

for cmd in truvari bgzip tabix gzip awk sort grep sed wc; do
    command -v "$cmd" >/dev/null 2>&1 || {
        echo "ERROR: required command not found: $cmd" >&2
        exit 1
    }
done

TMP=$(mktemp -d ./truvari_trio_tmp.XXXXXX)
SUCCESS=0

cleanup() {
    local rc=$?
    if [[ "$KEEP_TEMP" -eq 1 || "$SUCCESS" -ne 1 || "$rc" -ne 0 ]]; then
        echo "[truvari_trio] Temporary files kept in: $TMP" >&2
    else
        rm -rf "$TMP"
    fi
}
trap cleanup EXIT

echo "[truvari_trio] Truvari executable: $(command -v truvari)" >&2
echo "[truvari_trio] temp directory: $TMP" >&2
echo "[truvari_trio] refdist: $DIST" >&2
echo "[truvari_trio] parallel bench jobs: $JOBS" >&2

# Stream a plain or gzip-compressed VCF to stdout.
stream_vcf() {
    local input="$1"
    if [[ "$input" == *.gz ]]; then
        gzip -cd -- "$input"
    else
        cat -- "$input"
    fi
}

# Prepare a VCF for Truvari without calling the external bcftools executable.
# This also adds missing ##contig lines for contigs seen in the body.
prep_vcf() {
    local input="$1"
    local output="$2"
    local label="$3"

    local raw="$TMP/${label}.raw.vcf"
    local meta="$TMP/${label}.meta.txt"
    local chrom_header="$TMP/${label}.chrom_header.txt"
    local body="$TMP/${label}.body.vcf"
    local body_sorted="$TMP/${label}.body.sorted.vcf"
    local used_contigs="$TMP/${label}.used_contigs.txt"
    local existing_contigs="$TMP/${label}.existing_contigs.txt"
    local missing_contigs="$TMP/${label}.missing_contigs.txt"
    local fixed="$TMP/${label}.fixed.vcf"

    echo "[truvari_trio] preparing $label: $input" >&2

    stream_vcf "$input" > "$raw"

    grep '^##' "$raw" > "$meta" || true
    grep '^#CHROM' "$raw" > "$chrom_header" || true

    if [[ ! -s "$chrom_header" ]]; then
        echo "ERROR: no #CHROM header line found in $input" >&2
        return 1
    fi

    awk -F'\t' '!/^#/ {print}' "$raw" > "$body"

    awk -F'\t' 'NF > 1 && !/^#/ {print $1}' "$body" \
        | LC_ALL=C sort -u > "$used_contigs"

    grep '^##contig=<ID=' "$meta" \
        | sed -E 's/^##contig=<ID=([^,>]+).*/\1/' \
        | LC_ALL=C sort -u > "$existing_contigs" || true

    comm -23 "$used_contigs" "$existing_contigs" > "$missing_contigs" || true

    if [[ -s "$missing_contigs" ]]; then
        local nmiss
        nmiss=$(wc -l < "$missing_contigs" | awk '{print $1}')
        echo "[truvari_trio] $label: adding $nmiss missing ##contig header entries" >&2
    fi

    {
        cat "$meta"
        awk '{print "##contig=<ID=" $0 ">"}' "$missing_contigs"
        cat "$chrom_header"
    } > "$fixed"

    # Truvari/tabix require coordinate-sorted input. Lexical contig order is fine
    # as long as each contig is contiguous and positions increase within contig.
    LC_ALL=C sort -t $'\t' -k1,1 -k2,2n -k3,3 "$body" > "$body_sorted"
    cat "$body_sorted" >> "$fixed"

    bgzip -c "$fixed" > "$output"
    tabix -f -p vcf "$output"

    rm -f "$raw" "$meta" "$chrom_header" "$body" "$body_sorted" \
          "$used_contigs" "$existing_contigs" "$missing_contigs" "$fixed"
}

# Prepare all six haplotypes.
prep_vcf "${CHILD[0]}" "$TMP/ch1.vcf.gz" "ch1"
prep_vcf "${CHILD[1]}" "$TMP/ch2.vcf.gz" "ch2"
prep_vcf "${MOTHER[0]}" "$TMP/m1.vcf.gz" "m1"
prep_vcf "${MOTHER[1]}" "$TMP/m2.vcf.gz" "m2"
prep_vcf "${FATHER[0]}" "$TMP/f1.vcf.gz" "f1"
prep_vcf "${FATHER[1]}" "$TMP/f2.vcf.gz" "f2"

# Extract one haplotype's ##referenceCoverage metadata as a raw BED. The
# generating pipeline declares these coordinates as 0-based, half-open, which
# is already BED's coordinate system. A VCF with multiple coverage samples is
# rejected because this script requires one haplotype per input VCF.
extract_reference_coverage_bed() {
    local input="$1"
    local output="$2"
    local label="$3"

    if ! tabix -H "$input" | awk -v OFS='\t' -v label="$label" '
        function meta_value(record, key, marker, rest, pos) {
            marker=key "=\""
            pos=index(record, marker)
            if (pos) {
                rest=substr(record, pos + length(marker))
                pos=index(rest, "\"")
                return pos ? substr(rest, 1, pos - 1) : ""
            }

            marker=key "="
            pos=index(record, marker)
            if (!pos)
                return ""
            rest=substr(record, pos + length(marker))
            pos=index(rest, ",")
            if (!pos)
                pos=index(rest, ">")
            return pos ? substr(rest, 1, pos - 1) : rest
        }

        /^##referenceCoverageCoordinateSystem=/ {
            coordinate_system=$0
            sub(/^##referenceCoverageCoordinateSystem=/, "", coordinate_system)
            sub(/\r$/, "", coordinate_system)
            next
        }

        /^##referenceCoverage=</ {
            sample=meta_value($0, "Sample")
            chrom=meta_value($0, "Chrom")
            start=meta_value($0, "Start")
            end=meta_value($0, "End")

            if (sample == "" || chrom == "" || start !~ /^[0-9]+$/ ||
                    end !~ /^[0-9]+$/) {
                print "ERROR: malformed ##referenceCoverage header in " label ": " $0 > "/dev/stderr"
                bad=1
                next
            }

            samples[sample]=1
            if ((end + 0) > (start + 0)) {
                print chrom, start + 0, end + 0
                interval_count++
            }
            next
        }

        /^#CHROM/ {next}

        END {
            if (coordinate_system != "" &&
                    coordinate_system != "0-based-half-open") {
                print "ERROR: unsupported ##referenceCoverageCoordinateSystem in " label ": " coordinate_system > "/dev/stderr"
                bad=1
            }
            if (coordinate_system == "") {
                print "[truvari_trio] warning: " label " has no ##referenceCoverageCoordinateSystem; treating Start/End as 0-based, half-open" > "/dev/stderr"
            }

            sample_count=0
            for (sample in samples)
                sample_count++
            if (sample_count != 1) {
                print "ERROR: " label " must contain ##referenceCoverage records for exactly one sample; found " sample_count > "/dev/stderr"
                bad=1
            }
            if (interval_count == 0) {
                print "ERROR: no non-empty ##referenceCoverage intervals found in " label > "/dev/stderr"
                bad=1
            }
            if (bad)
                exit 1
        }
    ' > "$output"; then
        echo "ERROR: failed to extract ##referenceCoverage from $label" >&2
        rm -f "$output"
        return 1
    fi
}

# Sort and merge overlapping or directly adjacent BED intervals.
normalize_bed() {
    local input="$1"
    local output="$2"
    local sorted="$output.sorted"

    LC_ALL=C sort -t $'\t' -k1,1 -k2,2n -k3,3n "$input" > "$sorted"
    awk -F'\t' -v OFS='\t' '
        function emit() {
            if (have)
                print chrom, start, end
        }
        {
            if (!have) {
                chrom=$1
                start=$2
                end=$3
                have=1
            } else if ($1 == chrom && $2 <= end) {
                if ($3 > end)
                    end=$3
            } else {
                emit()
                chrom=$1
                start=$2
                end=$3
            }
        }
        END {emit()}
    ' "$sorted" > "$output"
    rm -f "$sorted"
}

# Intersect two normalized BED files. Intervals from the first file are kept in
# per-contig arrays; the cursor advances monotonically because the second file
# is also sorted. This avoids requiring bedtools.
intersect_beds() {
    local left="$1"
    local right="$2"
    local output="$3"

    awk -F'\t' -v OFS='\t' '
        NR == FNR {
            count[$1]++
            key=$1 SUBSEP count[$1]
            left_start[key]=$2
            left_end[key]=$3
            next
        }
        {
            chrom=$1
            right_start=$2
            right_end=$3
            i=(chrom in cursor) ? cursor[chrom] : 1

            while (i <= count[chrom] &&
                    left_end[chrom SUBSEP i] <= right_start)
                i++
            cursor[chrom]=i

            for (j=i; j <= count[chrom] &&
                    left_start[chrom SUBSEP j] < right_end; j++) {
                start=left_start[chrom SUBSEP j] > right_start \
                    ? left_start[chrom SUBSEP j] : right_start
                end=left_end[chrom SUBSEP j] < right_end \
                    ? left_end[chrom SUBSEP j] : right_end
                if (end > start)
                    print chrom, start, end
            }
        }
    ' "$left" "$right" > "$output"
}

build_parent_include_bed() {
    local label raw normalized
    for label in m1 m2 f1 f2; do
        raw="$TMP/${label}.coverage.raw.bed"
        normalized="$TMP/${label}.coverage.bed"
        extract_reference_coverage_bed "$TMP/${label}.vcf.gz" "$raw" "$label"
        normalize_bed "$raw" "$normalized"
        rm -f "$raw"
    done

    local mother="$TMP/parents.m1_m2.bed"
    local mother_father1="$TMP/parents.m1_m2_f1.bed"
    INCLUDE_BED_PATH="$TMP/parents.all4.covered.bed"

    intersect_beds "$TMP/m1.coverage.bed" "$TMP/m2.coverage.bed" "$mother"
    if [[ ! -s "$mother" ]]; then
        echo "ERROR: mother_h1 and mother_h2 have no shared referenceCoverage region" >&2
        return 1
    fi

    intersect_beds "$mother" "$TMP/f1.coverage.bed" "$mother_father1"
    if [[ ! -s "$mother_father1" ]]; then
        echo "ERROR: the maternal intersection and father_h1 have no shared referenceCoverage region" >&2
        return 1
    fi

    intersect_beds "$mother_father1" "$TMP/f2.coverage.bed" "$INCLUDE_BED_PATH"
    if [[ ! -s "$INCLUDE_BED_PATH" ]]; then
        echo "ERROR: the four parental haplotypes have no shared referenceCoverage region" >&2
        return 1
    fi

    local interval_count covered_bases
    interval_count=$(wc -l < "$INCLUDE_BED_PATH" | awk '{print $1}')
    covered_bases=$(awk -F'\t' '{n += $3 - $2} END {print n + 0}' "$INCLUDE_BED_PATH")
    echo "[truvari_trio] all-four-parent include BED: $INCLUDE_BED_PATH" >&2
    echo "[truvari_trio] include BED: $interval_count intervals, $covered_bases covered bases" >&2
}

build_parent_include_bed

# One pairwise Truvari benchmark. stdout/stderr go into a per-job log so that
# eight parallel jobs do not interleave their output.
run_bench() {
    local child="$1"
    local parent="$2"
    local out="$3"
    local name="$4"
    local log="$TMP/${name}.bench.log"
    local -a bench_args=(
        truvari bench
        -b "$child"
        -c "$parent"
        --refdist "$DIST"
    )

    if [[ -n "$INCLUDE_BED_PATH" ]]; then
        bench_args+=(--includebed "$INCLUDE_BED_PATH")
    fi
    bench_args+=(-o "$out")

    "${bench_args[@]}" > "$log" 2>&1
}

# Run the eight independent comparisons in parallel, in batches of --jobs.
run_all_benches() {
    local -a specs=(
        "ch1 m1"
        "ch1 m2"
        "ch1 f1"
        "ch1 f2"
        "ch2 m1"
        "ch2 m2"
        "ch2 f1"
        "ch2 f2"
    )

    local -a pids=()
    local -a names=()
    local spec child parent name out pid
    local failed=0

    wait_batch() {
        local i rc=0
        for i in "${!pids[@]}"; do
            if wait "${pids[$i]}"; then
                echo "[truvari_trio] completed ${names[$i]}" >&2
            else
                echo "ERROR: Truvari bench failed: ${names[$i]}" >&2
                echo "----- ${names[$i]} log -----" >&2
                cat "$TMP/${names[$i]}.bench.log" >&2 || true
                echo "----- end log -----" >&2
                rc=1
            fi
        done
        pids=()
        names=()
        return "$rc"
    }

    for spec in "${specs[@]}"; do
        child=${spec%% *}
        parent=${spec##* }
        name="${child}_vs_${parent}"
        out="$TMP/$name"

        echo "[truvari_trio] launching $name" >&2
        run_bench "$TMP/${child}.vcf.gz" "$TMP/${parent}.vcf.gz" "$out" "$name" &
        pid=$!
        pids+=("$pid")
        names+=("$name")

        if [[ ${#pids[@]} -ge "$JOBS" ]]; then
            if ! wait_batch; then
                failed=1
            fi
        fi
    done

    if [[ ${#pids[@]} -gt 0 ]]; then
        if ! wait_batch; then
            failed=1
        fi
    fi

    if [[ "$failed" -ne 0 ]]; then
        echo "ERROR: one or more Truvari bench jobs failed" >&2
        return 1
    fi
}

run_all_benches

# Count non-header VCF records, plain or gzip-compressed.
count_vcf() {
    local input="$1"
    if [[ "$input" == *.gz ]]; then
        gzip -cd -- "$input" | awk '!/^#/ {n++} END {print n+0}'
    else
        awk '!/^#/ {n++} END {print n+0}' "$input"
    fi
}

# Count rows in a truvari consistency --output TSV whose COUNT field equals n.
# TSV columns are: CHROM POS ID REF ALT FLAG COUNT
count_consistency_n() {
    local input="$1"
    local n="$2"
    awk -F'\t' -v n="$n" '!/^#/ && $7 == n {c++} END {print c+0}' "$input"
}

# Analyze one child haplotype.
analyze_hap() {
    local h="$1"
    local M1="$TMP/ch${h}_vs_m1"
    local M2="$TMP/ch${h}_vs_m2"
    local F1="$TMP/ch${h}_vs_f1"
    local F2="$TMP/ch${h}_vs_f2"

    local mother_tsv="$TMP/ch${h}.mother_missing.tsv"
    local father_tsv="$TMP/ch${h}.father_missing.tsv"
    local all_tsv="$TMP/ch${h}.all_missing.tsv"

    # These FN files all preserve the original child representation. Therefore,
    # exact-key consistency identifies child calls missed by both/all comparisons.
    truvari consistency -o "$mother_tsv" \
        "$M1/fn.vcf.gz" \
        "$M2/fn.vcf.gz" \
        > "$TMP/ch${h}.mother_missing.report.txt"

    truvari consistency -o "$father_tsv" \
        "$F1/fn.vcf.gz" \
        "$F2/fn.vcf.gz" \
        > "$TMP/ch${h}.father_missing.report.txt"

    truvari consistency -o "$all_tsv" \
        "$M1/fn.vcf.gz" \
        "$M2/fn.vcf.gz" \
        "$F1/fn.vcf.gz" \
        "$F2/fn.vcf.gz" \
        > "$TMP/ch${h}.all_missing.report.txt"

    local tp_m1 tp_m2 tp_f1 tp_f2 fn total
    local miss_m miss_f miss_all found_m found_f found_both found_any
    local consistency error_rate

    tp_m1=$(count_vcf "$M1/tp-base.vcf.gz")
    tp_m2=$(count_vcf "$M2/tp-base.vcf.gz")
    tp_f1=$(count_vcf "$F1/tp-base.vcf.gz")
    tp_f2=$(count_vcf "$F2/tp-base.vcf.gz")

    # Benchmarkable child total: TP-base + FN from any one comparison.
    # The baseline is the same child haplotype in all four comparisons.
    fn=$(count_vcf "$M1/fn.vcf.gz")
    total=$((tp_m1 + fn))

    # Missing from both maternal / both paternal / all four parental haplotypes.
    miss_m=$(count_consistency_n "$mother_tsv" 2)
    miss_f=$(count_consistency_n "$father_tsv" 2)
    miss_all=$(count_consistency_n "$all_tsv" 4)

    found_m=$((total - miss_m))
    found_f=$((total - miss_f))
    found_any=$((total - miss_all))
    found_both=$((total - miss_m - miss_f + miss_all))

    consistency=$(awk -v x="$found_any" -v n="$total" \
        'BEGIN {if(n) printf "%.6f", x/n; else print "NA"}')
    error_rate=$(awk -v x="$miss_all" -v n="$total" \
        'BEGIN {if(n) printf "%.6f", x/n; else print "NA"}')

    echo -e "child_h${h}_total\t$total"
    echo -e "child_h${h}_matched_mother_h1\t$tp_m1"
    echo -e "child_h${h}_matched_mother_h2\t$tp_m2"
    echo -e "child_h${h}_matched_father_h1\t$tp_f1"
    echo -e "child_h${h}_matched_father_h2\t$tp_f2"
    echo -e "child_h${h}_in_mother\t$found_m"
    echo -e "child_h${h}_in_father\t$found_f"
    echo -e "child_h${h}_in_both_parents\t$found_both"
    echo -e "child_h${h}_in_either_parent\t$found_any"
    echo -e "child_h${h}_not_in_parents\t$miss_all"
    echo -e "child_h${h}_error_rate\t$error_rate"
    echo -e "child_h${h}_trio_consistency\t$consistency"

    printf '%s\t%s\t%s\t%s\t%s\t%s\n' \
        "$total" "$found_m" "$found_f" "$found_both" "$found_any" "$miss_all" \
        > "$TMP/ch${h}.numbers"
}

analyze_hap 1
echo
analyze_hap 2

# Combined h1+h2 summary. Haplotype calls remain separate; this is just a sum.
IFS=$'\t' read -r t1 m1 f1 b1 a1 x1 < "$TMP/ch1.numbers"
IFS=$'\t' read -r t2 m2 f2 b2 a2 x2 < "$TMP/ch2.numbers"

TOTAL=$((t1 + t2))
MOTHER_FOUND=$((m1 + m2))
FATHER_FOUND=$((f1 + f2))
BOTH_FOUND=$((b1 + b2))
ANY_FOUND=$((a1 + a2))
MISSING=$((x1 + x2))

FRAC=$(awk -v x="$ANY_FOUND" -v n="$TOTAL" \
    'BEGIN {if(n) printf "%.6f", x/n; else print "NA"}')
ERRFRAC=$(awk -v x="$MISSING" -v n="$TOTAL" \
    'BEGIN {if(n) printf "%.6f", x/n; else print "NA"}')

echo
echo -e "total_child\t$TOTAL"
echo -e "child_in_mother\t$MOTHER_FOUND"
echo -e "child_in_father\t$FATHER_FOUND"
echo -e "child_in_both_parents\t$BOTH_FOUND"
echo -e "child_in_either_parent\t$ANY_FOUND"
echo -e "child_not_in_parents\t$MISSING"
echo -e "error_rate\t$ERRFRAC"
echo -e "trio_consistency\t$FRAC"

# Insert .h1/.h2 before .vcf or .vcf.gz.
fp_hap_name() {
    local fp="$1"
    local hap="$2"
    case "$fp" in
        *.vcf.gz) printf '%s.%s.vcf.gz\n' "${fp%.vcf.gz}" "$hap" ;;
        *.vcf)    printf '%s.%s.vcf\n'    "${fp%.vcf}" "$hap" ;;
        *)        printf '%s.%s.vcf\n'    "$fp" "$hap" ;;
    esac
}

# Write a VCF (plain or bgzip) from a temporary plain VCF.
write_vcf_output() {
    local src="$1"
    local dst="$2"
    mkdir -p "$(dirname "$dst")"

    if [[ "$dst" == *.gz ]]; then
        bgzip -c "$src" > "$dst"
        tabix -f -p vcf "$dst"
    else
        cp "$src" "$dst"
    fi
}

# Write one child haplotype's FP/error VCF, preserving the full child header and
# full original child records. Error = absent from all four parental haplotypes.
write_hap_fp() {
    local h="$1"
    local dst="$2"
    local child="$TMP/ch${h}.vcf.gz"
    local cons="$TMP/ch${h}.all_missing.tsv"
    local outtmp="$TMP/ch${h}.fp.vcf"

    gzip -cd -- "$child" | grep '^#' > "$outtmp"

    awk -F'\t' -v OFS='\t' '
        NR==FNR {
            if ($0 !~ /^#/ && $7 == 4) {
                key=$1 FS $2 FS $3 FS $4 FS $5
                wanted[key]=1
            }
            next
        }
        /^#/ {next}
        {
            key=$1 FS $2 FS $3 FS $4 FS $5
            if (key in wanted)
                print
        }
    ' "$cons" <(gzip -cd -- "$child") >> "$outtmp"

    write_vcf_output "$outtmp" "$dst"
}

# Write a combined sites-only FP VCF with TRIO_HAP indicating h1 or h2.
write_combined_fp() {
    local dst="$1"
    local raw="$TMP/fp.combined.raw.vcf"
    local header="$TMP/fp.combined.header.vcf"
    local body="$TMP/fp.combined.body.vcf"
    local sorted_body="$TMP/fp.combined.body.sorted.vcf"
    local final="$TMP/fp.combined.vcf"

    gzip -cd -- "$TMP/ch1.vcf.gz" | grep '^##' > "$header"
    echo '##INFO=<ID=TRIO_HAP,Number=1,Type=String,Description="Child haplotype source (h1 or h2)">' >> "$header"
    printf '#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n' >> "$header"

    : > "$body"

    emit_sites() {
        local h="$1"
        local child="$TMP/ch${h}.vcf.gz"
        local cons="$TMP/ch${h}.all_missing.tsv"

        awk -F'\t' -v OFS='\t' -v hap="h${h}" '
            NR==FNR {
                if ($0 !~ /^#/ && $7 == 4) {
                    key=$1 FS $2 FS $3 FS $4 FS $5
                    wanted[key]=1
                }
                next
            }
            /^#/ {next}
            {
                key=$1 FS $2 FS $3 FS $4 FS $5
                if (key in wanted) {
                    info=$8
                    if (info=="." || info=="")
                        info="TRIO_HAP=" hap
                    else
                        info=info ";TRIO_HAP=" hap
                    print $1,$2,$3,$4,$5,$6,$7,info
                }
            }
        ' "$cons" <(gzip -cd -- "$child")
    }

    emit_sites 1 >> "$body"
    emit_sites 2 >> "$body"

    LC_ALL=C sort -t $'\t' -k1,1 -k2,2n -k3,3 "$body" > "$sorted_body"
    cat "$header" "$sorted_body" > "$final"

    write_vcf_output "$final" "$dst"
}

if [[ -n "$FP" ]]; then
    FP_H1=$(fp_hap_name "$FP" "h1")
    FP_H2=$(fp_hap_name "$FP" "h2")

    write_hap_fp 1 "$FP_H1"
    write_hap_fp 2 "$FP_H2"
    write_combined_fp "$FP"

    echo "[truvari_trio] child h1 FP VCF: $FP_H1" >&2
    echo "[truvari_trio] child h2 FP VCF: $FP_H2" >&2
    echo "[truvari_trio] combined FP VCF: $FP" >&2
fi

SUCCESS=1
exit 0
