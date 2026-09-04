#!/bin/bash

set -e
set -o pipefail

INPUT_VCF_DIR="../stats/data/msmc/vcf_msmc"
INPUT_MASK_DIR="../stats/data/msmc/mask_ind"
OUTPUT_DIR="../stats/data/msmc/msmc_input"
GLOBAL_MASK="../mask/merged/autosome.bed.gz"
MASTER_LIST="../list/master/master_kinship_filtered.txt"
SAMPLE_LIST="../list/jm_msmc.txt"
CHROM_LIST="../list/chromosome.txt"
MSMC_TOOL="../python/msmc-tools/generate_multihetsep.py"
PARALLEL_JOBS=20

# Prepare sample list
# Exclude Okazaki and Takasakiyama beceause N=2 samples are not available for these populations
# For Yakushima, only wild samples are included
awk -F'\t' '
    NR > 1 && $7 == "fuscata" && $9 !~ /^(Takasakiyama|Okazaki)$/ && $11 == "wild" { print $1, $9 }
' OFS='\t' "${MASTER_LIST}" > "${SAMPLE_LIST}"

# Function to get N samples
get_n_samples() {
    local pop=$1
    local chrom=$2
    local n=$3
    awk -v VCF_DIR="${INPUT_VCF_DIR}" -v MASK_DIR="${INPUT_MASK_DIR}" -v chrom="${chrom}" -v pop="${pop}" -v n="$n" '
    BEGIN { count = 0 }
    $2==pop && count<n {
        count++
        samples = (samples ? samples " " : "") VCF_DIR "/" chrom "/" $1 ".vcf.gz"
        masks = (masks ? masks " --mask=" : "--mask=") MASK_DIR "/" chrom "/" $1 ".bed.gz"
        names = (names ? names " " : "") $1
    }
    END {
        if (count == n) {
            print samples "\t" masks "\t" names
        } else {
            print "\t\t"
        }
    }
    ' "${SAMPLE_LIST}"
}

process_chrom() {
    local chrom=$1

    # Two samples for each population
    for pop in $(awk '{print $2}' "${SAMPLE_LIST}" | sort | uniq); do
        IFS=$'\t' read -r samples masks names < <(get_n_samples "${pop}" "${chrom}" 2)
        [ -z "$samples" ] && continue
        outdir="${OUTPUT_DIR}/pop/${pop}_2samples"
        mkdir -p "${outdir}"
        echo "${pop}:${names}" > "${outdir}/${chrom}.log"
        python "${MSMC_TOOL}" \
            --mask="${GLOBAL_MASK}" \
            ${masks} \
            ${samples} \
            > "${outdir}/${chrom}.msmc" \
            2>> "${outdir}/${chrom}.log"
    done

    # Yakushima - Other population pair
    IFS=$'\t' read -r yaku_samples yaku_masks yaku_names < <(get_n_samples "Yakushima" "${chrom}" 2)
    if [ -n "$yaku_samples" ]; then
        for pop in $(awk '$2!="Yakushima"{print $2}' "${SAMPLE_LIST}" | sort | uniq); do
            IFS=$'\t' read -r samples masks names < <(get_n_samples "${pop}" "${chrom}" 2)
            [ -z "$samples" ] && continue
            outdir="${OUTPUT_DIR}/div/Yakushima_${pop}"
            mkdir -p "${outdir}"
            {
                echo "Yakushima:${yaku_names}"
                echo "${pop}:${names}"
            } > "${outdir}/${chrom}.log"
            python "${MSMC_TOOL}" \
                --mask="${GLOBAL_MASK}" \
                ${yaku_masks} ${masks} \
                ${yaku_samples} ${samples} \
                > "${outdir}/${chrom}.msmc" \
                2>> "${outdir}/${chrom}.log"
        done
    fi

    # Western and Eastern population pair
    for pair in "Minoh Hagachizaki" "Minoh Shimokita" "Kojima Hagachizaki" "Kojima Shimokita"; do
        read -r pop1 pop2 <<< "$pair"
        
        IFS=$'\t' read -r samples1 masks1 names1 < <(get_n_samples "${pop1}" "${chrom}" 2)
        IFS=$'\t' read -r samples2 masks2 names2 < <(get_n_samples "${pop2}" "${chrom}" 2)
        
        [[ -z "$samples1" || -z "$samples2" ]] && continue
        
        outdir="${OUTPUT_DIR}/div/${pop1}_${pop2}"
        mkdir -p "${outdir}"
        {
            echo "${pop1}:${names1}"
            echo "${pop2}:${names2}"
        } > "${outdir}/${chrom}.log"
        python "${MSMC_TOOL}" \
            --mask="${GLOBAL_MASK}" \
            ${masks1} ${masks2} \
            ${samples1} ${samples2} \
            > "${outdir}/${chrom}.msmc" \
            2>> "${outdir}/${chrom}.log"
    done

    # Yakushima 5 samples
    IFS=$'\t' read -r samples masks names < <(get_n_samples "Yakushima" "${chrom}" 5)
    if [ -n "$samples" ]; then
        outdir="${OUTPUT_DIR}/pop/Yakushima_5samples"
        mkdir -p "${outdir}"
        echo "Yakushima:${names}" > "${outdir}/${chrom}.log"
        python "${MSMC_TOOL}" \
            --mask="${GLOBAL_MASK}" \
            ${masks} \
            ${samples} \
            > "${outdir}/${chrom}.msmc" \
            2>> "${outdir}/${chrom}.log"
    fi
}

export -f process_chrom
export -f get_n_samples
export INPUT_VCF_DIR INPUT_MASK_DIR SAMPLE_LIST GLOBAL_MASK OUTPUT_DIR MSMC_TOOL

parallel --halt soon,fail=1 -j"${PARALLEL_JOBS}" process_chrom {} ::: $(awk 'NR>1 {print $7}' "${CHROM_LIST}" | grep -Ev "NC_005943.1|NC_041774.1|NC_027914.1")
