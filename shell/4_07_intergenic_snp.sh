#!/bin/bash
set -euo pipefail

# Define variables
INPUT_DIR="../vcf/derived/all"
OUTPUT_DIR="../vcf/derived_intergenic"
CHROM_LIST="../list/chromosome.txt"
GFF_FILE="../refseq/GCF_003339765.1_Mmul_10_genomic.gff.gz"
REF_FAI="../refseq/GCF_003339765.1_Mmul_10_genomic.fna.fai"
PARALLEL_JOBS=21

mkdir -p "${OUTPUT_DIR}/mask" "${OUTPUT_DIR}/vcf"

# Make mask for intergenic regions
zcat "${GFF_FILE}" | \
    awk -v OFS='\t' '$3 == "gene" {print $1, $4-1, $5}' | \
    sort -k1,1 -k2,2n > "${OUTPUT_DIR}/mask/genes.bed"

cut -f1,2 "${REF_FAI}" | sort -k1,1 > "${OUTPUT_DIR}/mask/genome.file"

bedtools complement -i "${OUTPUT_DIR}/mask/genes.bed" -g "${OUTPUT_DIR}/mask/genome.file" \
    > "${OUTPUT_DIR}/mask/intergenic.bed"

extract_intergenic_snps() {
    local chrom=$1
    echo "Processing ${chrom}..."
    
    local input_vcf="${INPUT_DIR}/${chrom}.vcf.gz"
    local output_vcf="${OUTPUT_DIR}/vcf/${chrom}.vcf.gz"
    local log_file="${OUTPUT_DIR}/vcf/${chrom}.log"
    
    bcftools view --threads 3 -R "${OUTPUT_DIR}/mask/intergenic.bed" \
        -Oz -o "${output_vcf}" \
        "${input_vcf}" 2> "${log_file}"
        
    bcftools index --threads 3 -f -t "${output_vcf}"
}

export -f extract_intergenic_snps
export INPUT_DIR OUTPUT_DIR

parallel -j "${PARALLEL_JOBS}" extract_intergenic_snps ::: \
  $(awk 'NR>1 {print $7}' "${CHROM_LIST}" | grep -Ev "NC_005943.1|NC_027914.1")