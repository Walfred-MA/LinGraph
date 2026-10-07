#!/usr/bin/env bash
# Merge the two haplotype VCFs of one sample into one site-only VCF for
# benchmarking (e.g. truvari against GIAB).
#
# Usage (run in the output folder):
#   bash tools/BCFtoolMergeTwoHaplotypes.sh HG002_h1.vcf HG002_h2.vcf
#
# Writes HG002_h1.sites.vcf.gz, HG002_h2.sites.vcf.gz and HG002.LinGraph.vcf.gz
# (with .tbi) in the current folder: each haplotype sorted without genotypes,
# then their union with exact duplicates removed. Needs bcftools and tabix.
set -euo pipefail

if [[ $# -ne 2 ]]; then
    sed -n '2,10p' "$0" >&2
    exit 1
fi

H1=$1
H2=$2

# h1
bcftools sort "$H1" -Ou | \
    bcftools view -G -Oz -o HG002_h1.sites.vcf.gz
tabix -f -p vcf HG002_h1.sites.vcf.gz

# h2
bcftools sort "$H2" -Ou | \
    bcftools view -G -Oz -o HG002_h2.sites.vcf.gz
tabix -f -p vcf HG002_h2.sites.vcf.gz

# Union h1 + h2, removing exact duplicates
bcftools concat \
    -a \
    -d exact \
    HG002_h1.sites.vcf.gz \
    HG002_h2.sites.vcf.gz \
    -Ou | \
bcftools sort -Oz -o HG002.LinGraph.vcf.gz

tabix -f -p vcf HG002.LinGraph.vcf.gz
