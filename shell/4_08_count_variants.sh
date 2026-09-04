#!/bin/bash

set -e
set -o pipefail

DIR_SNP="../vcf/bial_snp/snp"
DIR_INDEL="../vcf/bial_indel/indel"
OUTPUT_DIR="../stats/count_variants"

mkdir -p "$OUTPUT_DIR"

a_snp=0
x_snp=0
for f in "$DIR_SNP"/*.vcf.gz; do
    if [[ "$(basename "$f")" == "NC_041774.1.vcf.gz" ]]; then
        x_n=$(bcftools view -H "$f" | wc -l)
        a_n=0
    else
        a_n=$(bcftools view -H "$f" | wc -l)
        x_n=0
    fi
    a_snp=$((a_snp + a_n))
    x_snp=$((x_snp + x_n))
done
echo "SNP: $a_snp" > "$OUTPUT_DIR/count_variants.txt"
echo "SNP: $x_snp" >> "$OUTPUT_DIR/count_variants.txt"

a_indel=0
x_indel=0
for f in "$DIR_INDEL"/*.vcf.gz; do
    if [[ "$(basename "$f")" == "NC_041774.1.vcf.gz" ]]; then
        x_n=$(bcftools view -H "$f" | wc -l)
        a_n=0
    else
        a_n=$(bcftools view -H "$f" | wc -l)
        x_n=0
    fi
    a_indel=$((a_indel + a_n))
    x_indel=$((x_indel + x_n))
done
echo "INDEL: $a_indel" >> "$OUTPUT_DIR/count_variants.txt"
echo "INDEL: $x_indel" >> "$OUTPUT_DIR/count_variants.txt"
