#!/bin/bash

set -e
set -o pipefail

# samtools faidx ../refseq/GCF_003339765.1_Mmul_10_genomic.fna

# Create a sequence dictionary for the reference genome
#gatk CreateSequenceDictionary \
#    -R ../refseq/GCF_003339765.1_Mmul_10_genomic.fna \
#    -O ../refseq/GCF_003339765.1_Mmul_10_genomic.dict

# Define directories
SAMPLE_LIST="../list/sample_sex_raw.txt"
CHROM_LIST="../list/chromosome.txt"
INPUT_DIR="../bam/mark_dp"
OUTPUT_DIR="../vcf/gvcf"
REFSEQ="../refseq/GCF_003339765.1_Mmul_10_genomic.fna"
PARALLEL_JOBS=25
MEMORY="40G"

# Prepare output directories
awk '{print $1}' "${SAMPLE_LIST}" | sort -u | xargs -I{} mkdir -p "${OUTPUT_DIR}"/{}

# Function for variant calling for a sample-chromosome
haplotype_calling() {
    local sample_id=$1  # Sample ID
    local sex=$2        # Sex
    local chr=$3        # Chromosome name
    local ploidy=2      # Default ploidy for autosomes

    # Determine ploidy based on chromosome type
    case "$chr" in
        "NC_041774.1") # ChrX
            if [ "$sex" = "male" ]; then
                ploidy=1
            fi
            ;;
        "NC_027914.1" | "NC_005943.1") # ChrY or ChrM
            ploidy=1
            ;;
    esac

    echo -e "${sample_id}\t${sex}\t${chr}\t${ploidy}"

    if ! gatk --java-options "-Xmx${MEMORY}" HaplotypeCaller \
        -R "${REFSEQ}" \
        -I "${INPUT_DIR}/${sample_id}.bam" \
        -O "${OUTPUT_DIR}/${sample_id}/${chr}.g.vcf.gz" \
        -ploidy "${ploidy}" \
        -L "${chr}" \
        --emit-ref-confidence BP_RESOLUTION \
        --output-mode EMIT_ALL_ACTIVE_SITES \
        &> "${OUTPUT_DIR}/${sample_id}/${chr}.log"
    then
        echo "Error: HaplotypeCaller failed for ${sample_id} ${chr}" >&2
        exit 1
    fi
}
export -f haplotype_calling
export REFSEQ INPUT_DIR OUTPUT_DIR MEMORY

# Parallel execution {1}:sample_id {2}:sex {3}:chromosome
parallel --halt now,fail=1 -j"${PARALLEL_JOBS}" --colsep '\t' -a "${SAMPLE_LIST}" haplotype_calling {1} {2} {3} ::: $(awk 'NR>1 {print $7}' "${CHROM_LIST}")

# Check error logs 
grep -rn "Exception" "${OUTPUT_DIR}" --include="*.log" > "${OUTPUT_DIR}/error.log" || true
