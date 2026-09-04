#!/bin/bash

set -e
set -o pipefail

# Define variables
INPUT_DIR="../vcf/derived"
OUTPUT_DIR="../stats/data/genetic_load"
CHROM_LIST="../list/chromosome.txt"
MASTER_LIST="../list/master/master_kinship_filtered.txt"
SAMPLE_LIST="../list/pop_cluster.txt"
PARALLEL_JOBS=20

awk -F'\t' 'BEGIN{OFS="\t"} NR>1 {print $1, $9, $10}' "${MASTER_LIST}" > "${SAMPLE_LIST}"

mkdir -p "${OUTPUT_DIR}/sample_list"

awk '$3 == "East" {print $1}' "${SAMPLE_LIST}" > "${OUTPUT_DIR}/sample_list/East.txt"
awk '$3 == "West" {print $1}' "${SAMPLE_LIST}" > "${OUTPUT_DIR}/sample_list/West.txt"
awk '$3 == "Yakushima" {print $1}' "${SAMPLE_LIST}" > "${OUTPUT_DIR}/sample_list/Yakushima.txt"
awk '$3 == "Yakushima" && $1 !~ /^Mfy/ {print $1}' "${SAMPLE_LIST}" > "${OUTPUT_DIR}/sample_list/Yakushima_wild.txt"

process_chrom_type_effect() {
    local chrom=$1
    local type=$2
    local effect=$3

    if [ "${type}" == "SNP" ]; then
        local input_dir="${INPUT_DIR}/snp"
    elif [ "${type}" == "INDEL" ]; then
        local input_dir="${INPUT_DIR}/indel"
    else
        echo "Unknown type: ${type}"
        exit 1
    fi

    if [ "${effect}" == "synonymous" ]; then
        local filter="(ANN[0].EFFECT has 'synonymous_variant') && !(exists ANN[0].ERRORS)"
    elif [ "${effect}" == "missense" ]; then
        local filter="(ANN[0].EFFECT has 'missense_variant') && !(exists ANN[0].ERRORS)"
    elif [ "${effect}" == "lof" ]; then
        local filter="(exists LOF[0].PERC) && !(exists ANN[0].ERRORS)"
    elif [ "${effect}" == "intergenic" ]; then
        local filter="ANN[0].EFFECT has 'intergenic_region' && !(exists ANN[0].ERRORS)"
    else
        echo "Unknown effect: ${effect}"
        exit 1
    fi

    echo "Processing ${chrom} ${type} ${effect}"

    mkdir -p "${OUTPUT_DIR}/all/${effect}/${type}"
    mkdir -p "${OUTPUT_DIR}/East/${effect}/${type}"
    mkdir -p "${OUTPUT_DIR}/West/${effect}/${type}"
    mkdir -p "${OUTPUT_DIR}/Yakushima/${effect}/${type}"
    mkdir -p "${OUTPUT_DIR}/Yakushima_wild/${effect}/${type}"

    # Filter
    SnpSift filter "${filter}" \
        "${input_dir}/${chrom}.vcf.gz" \
        2> "${OUTPUT_DIR}/all/${effect}/${type}/${chrom}.log" | \
        bcftools +fill-tags -Oz -o "${OUTPUT_DIR}/all/${effect}/${type}/${chrom}.vcf.gz" \
        2>> "${OUTPUT_DIR}/all/${effect}/${type}/${chrom}.log"
    bcftools index "${OUTPUT_DIR}/all/${effect}/${type}/${chrom}.vcf.gz"
    
    # Create population-specific VCFs (recalculate AC/AN after subsetting)
    bcftools view -S "${OUTPUT_DIR}/sample_list/East.txt" -Ou \
        "${OUTPUT_DIR}/all/${effect}/${type}/${chrom}.vcf.gz" \
    | bcftools +fill-tags -Oz -o "${OUTPUT_DIR}/East/${effect}/${type}/${chrom}.vcf.gz" -- -t AC,AN,AF,NS,AC_Het,AC_Hom
    bcftools index -t "${OUTPUT_DIR}/East/${effect}/${type}/${chrom}.vcf.gz"

    bcftools view -S "${OUTPUT_DIR}/sample_list/West.txt" -Ou \
        "${OUTPUT_DIR}/all/${effect}/${type}/${chrom}.vcf.gz" \
    | bcftools +fill-tags -Oz -o "${OUTPUT_DIR}/West/${effect}/${type}/${chrom}.vcf.gz" -- -t AC,AN,AF,NS,AC_Het,AC_Hom
    bcftools index -t "${OUTPUT_DIR}/West/${effect}/${type}/${chrom}.vcf.gz"

    bcftools view -S "${OUTPUT_DIR}/sample_list/Yakushima.txt" -Ou \
        "${OUTPUT_DIR}/all/${effect}/${type}/${chrom}.vcf.gz" \
    | bcftools +fill-tags -Oz -o "${OUTPUT_DIR}/Yakushima/${effect}/${type}/${chrom}.vcf.gz" -- -t AC,AN,AF,NS,AC_Het,AC_Hom
    bcftools index -t "${OUTPUT_DIR}/Yakushima/${effect}/${type}/${chrom}.vcf.gz"

    bcftools view -S "${OUTPUT_DIR}/sample_list/Yakushima_wild.txt" -Ou \
        "${OUTPUT_DIR}/all/${effect}/${type}/${chrom}.vcf.gz" \
    | bcftools +fill-tags -Oz -o "${OUTPUT_DIR}/Yakushima_wild/${effect}/${type}/${chrom}.vcf.gz" -- -t AC,AN,AF,NS,AC_Het,AC_Hom
    bcftools index -t "${OUTPUT_DIR}/Yakushima_wild/${effect}/${type}/${chrom}.vcf.gz"
}

export -f process_chrom_type_effect
export INPUT_DIR OUTPUT_DIR

parallel -j"${PARALLEL_JOBS}" process_chrom_type_effect ::: \
    $(awk 'NR>1 {print $7}' "${CHROM_LIST}" | grep -Ev "NC_005943.1|NC_041774.1|NC_027914.1") ::: \
    "SNP" "INDEL" ::: \
    "synonymous" "missense" "lof" "intergenic"
