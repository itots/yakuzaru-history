#!/bin/bash

set -e
set -o pipefail

# Define variables
INPUT_DIR="../vcf/bial_snp/snp"
GLOBAL_MASK="../mask/merged/autosome.bed.gz"
OUTPUT_DIR="../stats/data/smc++"
MASTER_LIST="../list/master/master_kinship_filtered.txt"
YAKU_LIST="../list/yaku.txt"
REF_FAI="../refseq/GCF_003339765.1_Mmul_10_genomic.fna.fai"
CHROM_LIST="../list/chromosome.txt"
PARALLEL_JOBS=20

mkdir -p "${OUTPUT_DIR}/refseq" "${OUTPUT_DIR}/inds_intersect" "${OUTPUT_DIR}/inds_global" \
    "${OUTPUT_DIR}/negative_mask" "${OUTPUT_DIR}/smc++_input" "${OUTPUT_DIR}/vcf"

cut -f1,2 "${REF_FAI}" > "${OUTPUT_DIR}/refseq/sizes.txt"
awk -F'\t' 'BEGIN{OFS="\t"} NR>1 && $9=="Yakushima" {print $1}' "${MASTER_LIST}" > "${YAKU_LIST}"

process_chrom(){
    local chrom=$1
    echo "Processing ${chrom}"

    bcftools view --threads 3 -S "${YAKU_LIST}" \
        "${INPUT_DIR}/${chrom}.vcf.gz" | \
        bcftools +fill-tags - -- -t F_MISSING | \
        bcftools view -i 'F_MISSING=0' -Oz -o "${OUTPUT_DIR}/vcf/${chrom}.vcf.gz" \
        &> "${OUTPUT_DIR}/vcf/${chrom}.log"
    bcftools index -f -t "${OUTPUT_DIR}/vcf/${chrom}.vcf.gz"

    bcftools query -f '%CHROM\t%POS0\t%POS\n' "${OUTPUT_DIR}/vcf/${chrom}.vcf.gz" | \
        bedtools merge | \
        bgzip > "${OUTPUT_DIR}/inds_intersect/${chrom}.bed.gz"
    tabix -f -p bed "${OUTPUT_DIR}/inds_intersect/${chrom}.bed.gz"

    # Intersect with global mask
    bedtools intersect -a "${OUTPUT_DIR}/inds_intersect/${chrom}.bed.gz" -b "${GLOBAL_MASK}" | \
        bedtools merge | \
        bgzip > "${OUTPUT_DIR}/inds_global/${chrom}.bed.gz"
    tabix -f -p bed "${OUTPUT_DIR}/inds_global/${chrom}.bed.gz"

    local chrom_sizes="${OUTPUT_DIR}/refseq/${chrom}.sizes"
    awk -v c="${chrom}" '$1==c' "${OUTPUT_DIR}/refseq/sizes.txt" > "${chrom_sizes}"

    bedtools complement -i "${OUTPUT_DIR}/inds_global/${chrom}.bed.gz" -g "${chrom_sizes}" | \
        bgzip > "${OUTPUT_DIR}/negative_mask/${chrom}.bed.gz"
    tabix -f -p bed "${OUTPUT_DIR}/negative_mask/${chrom}.bed.gz"

    local list_samples
    list_samples=$(paste -sd, "${YAKU_LIST}")

    while read -r sample; do
        docker run --rm \
            --user "$(id -u):$(id -g)" \
            --entrypoint bash \
            -v "$(realpath "${OUTPUT_DIR}")":/mnt/smc_out \
            -w /mnt/smc_out \
            terhorst/smcpp:latest -c \
            "smc++ vcf2smc -d ${sample} ${sample} \
            -m negative_mask/${chrom}.bed.gz \
            vcf/${chrom}.vcf.gz \
            smc++_input/${chrom}_${sample}.smc.gz \
            ${chrom} Yakushima:${list_samples} \
            2> smc++_input/${chrom}_${sample}.log"
        if [[ $? -ne 0 ]]; then
            echo "Error: smc++ vcf2smc failed for ${chrom} ${sample}" >&2
        fi
    done < "${YAKU_LIST}"
}

export -f process_chrom
export INPUT_DIR GLOBAL_MASK OUTPUT_DIR YAKU_LIST REF_FAI CHROM_LIST

parallel -j"${PARALLEL_JOBS}" process_chrom ::: $(awk 'NR>1 {print $7}' "${CHROM_LIST}" | grep -Ev "NC_005943.1|NC_041774.1|NC_027914.1")
