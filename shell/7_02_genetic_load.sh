#!/bin/bash

set -euo pipefail

# Define variables
INPUT_DIR="../stats/data/genetic_load"
OUTPUT_DIR="../stats/genetic_load"
CHROM_LIST="../list/chromosome.txt"
PARALLEL_JOBS=20

# Function to process each chromosome, variant type, and effect
process_chrom_type_effect() {
    local chrom=$1
    local type=$2
    local effect=$3

    # Individual-level statistics
    mkdir -p "${OUTPUT_DIR}/individual/${effect}/${type}"

    local samples
    mapfile -t samples < <(bcftools query -l "${INPUT_DIR}/all/${effect}/${type}/${chrom}.vcf.gz")

    bcftools query -f '%CHROM\t%POS[\t%GT]\n' \
        "${INPUT_DIR}/all/${effect}/${type}/${chrom}.vcf.gz" | \
        mawk -v names="$(IFS=,; echo "${samples[*]}")" '
        BEGIN {
            split(names, arr_names, ",")
        }
        {
            for(i=3; i<=NF; i++) {
                split($i, arr_alleles, /[\/|]/)
                a1 = arr_alleles[1]; a2 = arr_alleles[2];

                all_sites[i-2]++;

                if ((a1 ~ /^[01]$/) && (a2 ~ /^[01]$/)) {
                    called_sites[i-2]++;

                    if ((a1=="1" && a2=="0") || (a1=="0" && a2=="1")) {
                        heterozygotes[i-2]++;
                        alleles[i-2]++;     
                        sites[i-2]++;       
                    }
 
                    else if (a1=="1" && a2=="1") {
                        homozygotes[i-2]++;
                        alleles[i-2]+=2; 
                        sites[i-2]++;
                    }
                }

            }
        }
        END {
            for(i=1; i<=length(arr_names); i++) {
                print arr_names[i], "heterozygotes", heterozygotes[i]+0
                print arr_names[i], "homozygotes",  homozygotes[i]+0
                print arr_names[i], "alleles",      alleles[i]+0
                print arr_names[i], "sites",        sites[i]+0
                print arr_names[i], "called_sites", called_sites[i]+0
                print arr_names[i], "all_sites",    all_sites[i]+0
            }
        }' > "${OUTPUT_DIR}/individual/${effect}/${type}/${chrom}.txt"

    # Population-level statistics
    for pop in East West Yakushima Yakushima_wild; do
        mkdir -p "${OUTPUT_DIR}/population/${pop}/${effect}/${type}"
        bcftools query -f '%CHROM\t%POS\t%AC\t%AN\n' \
            "${INPUT_DIR}/${pop}/${effect}/${type}/${chrom}.vcf.gz" \
            > "${OUTPUT_DIR}/population/${pop}/${effect}/${type}/${chrom}.txt"
    done
}

export -f process_chrom_type_effect
export INPUT_DIR OUTPUT_DIR

parallel -j"${PARALLEL_JOBS}" process_chrom_type_effect ::: \
    $(awk 'NR>1 {print $7}' "${CHROM_LIST}" | grep -Ev "NC_005943.1|NC_041774.1|NC_027914.1") ::: \
    "SNP" "INDEL" ::: \
    "synonymous" "missense" "lof" "intergenic"
