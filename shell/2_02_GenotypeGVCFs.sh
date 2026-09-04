#!/bin/bash

set -e
set -o pipefail

# Define variables
INPUT_DIR="../vcf/combined_gvcf"
OUTPUT_DIR="../vcf/joint"
REFSEQ="../refseq/GCF_003339765.1_Mmul_10_genomic.fna"
CHROM_LIST="../list/chromosome.txt"
PARALLEL_JOBS=24
MEMORY="20G"

# Create an output directory for the joint VCF files
mkdir -p "${OUTPUT_DIR}"

genotyping_gvcfs() {
    local chrom=$1
    local sex=${2:-}  # Default to empty (for autosomes and MT)
    local ploidy=2    # Default ploidy for autosomes

    # Determine ploidy based on chromosome type
    case "$chrom" in
        "NC_041774.1") # ChrX
            if [ "$sex" = "male" ]; then
                ploidy=1
            fi
            ;;
        "NC_027914.1" | "NC_005943.1") # ChrY or ChrM
            ploidy=1
            ;;
    esac

    echo "Genotyping: ${chrom} ${sex}"

    # Run GenotypeGVCFs
    if ! gatk --java-options "-Xmx${MEMORY}" GenotypeGVCFs \
        -R "${REFSEQ}" \
        -V "${INPUT_DIR}/${chrom}${sex:+_${sex}}.g.vcf.gz" \
        -O "${OUTPUT_DIR}/${chrom}${sex:+_${sex}}.vcf.gz" \
        -L "${chrom}" \
        -ploidy "${ploidy}" \
        -all-sites &> "${OUTPUT_DIR}/${chrom}${sex:+_${sex}}.log"
    then
        echo "Error: GenotypeGVCFs failed for ${chrom} ${sex}" >&2
        exit 1
    fi    
}

export -f genotyping_gvcfs
export REFSEQ INPUT_DIR OUTPUT_DIR MEMORY

{
    # Process autosomes and mitochondrial DNA (excluding X and Y chromosomes)
    awk 'NR>1 {print $7}' "${CHROM_LIST}" | grep -Ev "NC_041774.1|NC_027914.1"
    
    # Process X chromosome separately for males and females
    echo "NC_041774.1 male"
    echo "NC_041774.1 female"
    
    # Process Y chromosome (only for males)
    echo "NC_027914.1 male"
} | parallel --halt now,fail=1 --colsep ' ' -j"${PARALLEL_JOBS}" genotyping_gvcfs {1} {2}
