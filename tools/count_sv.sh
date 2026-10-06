#!/usr/bin/env bash
# Count INS and DEL SVs (>= MIN bp, default 50) in VCFs, plain or .gz.
#
# Columns are exclusive: INS_only / DEL_only have one part >= MIN, both has
# the inserted and the deleted part >= MIN; total = INS_only + DEL_only + 2*both.
#   deleted bp:  DEL -> |SVLEN|; others -> END - POS
#   inserted bp: INS -> SVLEN; others -> length of plain INFO/SEQ or SVINSSEQ
#   no SVLEN:    explicit REF/ALT -> deleted len(REF)-1, inserted len(ALT)-1
# Nested grVCF rows (CHROM = parent ID INS_/DEL_/SUB_/DUP_...) are skipped.
#
# Usage: bash tools/count_sv.sh [-m MIN] file.vcf[.gz] ...
set -euo pipefail

MIN=50
if [[ "${1:-}" == "-m" ]]; then
    MIN="$2"
    shift 2
fi
[[ $# -ge 1 ]] || { sed -n '2,11p' "$0" >&2; exit 1; }

printf 'file\tINS_only\tDEL_only\tboth\ttotal\n'
for f in "$@"; do
    if [[ "$f" == *.gz ]]; then gzip -cd -- "$f"; else cat -- "$f"; fi |
    awk -F'\t' -v OFS='\t' -v min="$MIN" -v name="$f" '
        /^#/ || $1 ~ /^(INS|DEL|SUB|DUP)_/ {next}
        {
            delete v
            n = split($8, a, ";")
            for (i = 1; i <= n; i++) {
                p = index(a[i], "=")
                if (p) v[substr(a[i], 1, p - 1)] = substr(a[i], p + 1)
            }
            t = v["SVTYPE"]
            seq = ("SEQ" in v) ? v["SEQ"] : v["SVINSSEQ"]
            if (v["SVLEN"] == "" || v["SVLEN"] == ".") {
                if ($5 ~ /^</ || $5 ~ /,/) next
                del = length($4) - 1
                ins = length($5) - 1
            } else {
                svlen = v["SVLEN"] + 0
                del = (t == "DEL") ? -svlen : ((v["END"] != "") ? v["END"] - $2 : 0)
                if (t == "INS") ins = svlen
                else if (seq ~ /^[ACGTNacgtn]+$/) ins = length(seq)
                else ins = 0
            }
            isins = (ins >= min); isdel = (del >= min)
            nins += (isins && !isdel); ndel += (isdel && !isins); nboth += (isins && isdel)
        }
        END {print name, nins + 0, ndel + 0, nboth + 0, nins + ndel + 2 * nboth}
    '
done
