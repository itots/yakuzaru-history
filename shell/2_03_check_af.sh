#!/bin/bash

set -e
set -o pipefail

# --- Settings ---
OUT_DIR="../stats/af_dist"
INPUT_DIR="../vcf/joint"
PARALLEL_JOBS=20
BIN_NUM=9

MIN_DP=20
MIN_QUAL=30

DIST_FILE="${OUT_DIR}/summary_af_dist.tsv"
RAW_DIR="${OUT_DIR}/raw_counts"
TMP_DIR="${RAW_DIR}/tmp"

# Get sample list
FIRST_VCF="${INPUT_DIR}/NC_041754.1.vcf.gz"
if [ ! -f "$FIRST_VCF" ]; then echo "Error: First VCF not found."; exit 1; fi

SAMPLE_LIST=$(bcftools query -l "$FIRST_VCF")
SAMPLES_ARRAY=($SAMPLE_LIST)
NUM_SAMPLES=${#SAMPLES_ARRAY[@]}

# Chromosome list
CHROMS=$(seq -f "NC_0417%g.1" 54 73)

mkdir -p "${RAW_DIR}"
for SAMPLE in "${SAMPLES_ARRAY[@]}"; do
    mkdir -p "${RAW_DIR}/${SAMPLE}"
done

mkdir -p "${TMP_DIR}"

# Function: Process one chromosome
calc_af_chrom() {
    local CHROM=$1
    local VCF="${INPUT_DIR}/${CHROM}.vcf.gz"
    
    if [ ! -f "$VCF" ]; then echo "Warning: $VCF not found"; return 0; fi
    
    echo "Processing ${CHROM}..."

    # Extract GT and AD    
    bcftools query \
        -i "TYPE=\"snp\" && QUAL>=${MIN_QUAL}" \
        -f '[%GT:%AD\t]\n' \
        "$VCF" | \
    mawk -v bins="$BIN_NUM" \
        -v tmp_dir="${TMP_DIR}" \
        -v chrom="${CHROM}" \
        -v num_samples="${NUM_SAMPLES}" \
        -v min_dp="${MIN_DP}" '
    BEGIN {
        for(s=1; s<=num_samples; s++) {
            total_het_sites[s] = 0;
            for(b=0; b<bins; b++) counts[s,b] = 0;
        }
    }
    {
        for(s=1; s<=NF; s++) {
            n = split($s, data, ":");
            if (n < 2) continue; 
            
            gt = data[1];
            
            if (gt == "0/1" || gt == "0|1" || gt == "1/0" || gt == "1|0") {
            
                split(data[2], ad, ",");
                ref = ad[1]; 
                alt = ad[2]; 
                
                if(ref != "." && alt != ".") {
                    depth = ref + alt;
                    
                    # Filter: Depth Check
                    if(depth >= min_dp) {
                         af = alt / depth;
                         bin = int(af * bins);
                         if(bin >= bins) bin = bins - 1;
                         counts[s,bin]++;
                         total_het_sites[s]++;
                    }
                }
            }
        }
    }
    END {
        for(s=1; s<=num_samples; s++) {
            outfile = tmp_dir "/" chrom "_sample_" s ".tmp"
            printf "%d", total_het_sites[s] > outfile;
            for(b=0; b<bins; b++) printf "\t%d", counts[s,b] > outfile;
            printf "\n" > outfile;
            close(outfile);
        }
    }'
}

export -f calc_af_chrom
export INPUT_DIR RAW_DIR TMP_DIR BIN_NUM NUM_SAMPLES SAMPLES_ARRAY MIN_DP MIN_QUAL

# Parallel execution
parallel --halt now,fail=1 -j "$PARALLEL_JOBS" calc_af_chrom ::: $CHROMS

# --- Post-processing (Rename) ---
echo "Renaming and organizing files..."
idx=1
for SAMPLE in "${SAMPLES_ARRAY[@]}"; do
    for CHROM in $CHROMS; do
        TMP_FILE="${TMP_DIR}/${CHROM}_sample_${idx}.tmp"
        TARGET_FILE="${RAW_DIR}/${SAMPLE}/${CHROM}.tsv"
        [ -f "$TMP_FILE" ] && mv "$TMP_FILE" "$TARGET_FILE"
    done
    ((idx++))
done

# --- Aggregation ---
echo "Aggregating results..."
echo -ne "Sample\tTotal_Het_Sites" > "$DIST_FILE"
for ((i=0; i<BIN_NUM; i++)); do 
    lower=$(awk "BEGIN {printf \"%.3f\", $i / $BIN_NUM}")
    upper=$(awk "BEGIN {printf \"%.3f\", ($i+1) / $BIN_NUM}")
    echo -ne "\tCount_${lower}-${upper}"
done >> "$DIST_FILE"
for ((i=0; i<BIN_NUM; i++)); do 
    lower=$(awk "BEGIN {printf \"%.3f\", $i / $BIN_NUM}")
    upper=$(awk "BEGIN {printf \"%.3f\", ($i+1) / $BIN_NUM}")
    echo -ne "\tFreq_${lower}-${upper}"
done >> "$DIST_FILE"
echo "" >> "$DIST_FILE"

for SAMPLE in "${SAMPLES_ARRAY[@]}"; do
    awk -v sample="$SAMPLE" -v bins="$BIN_NUM" '
    BEGIN { total=0; for(i=0; i<bins; i++) c[i]=0; }
    {
        total += $1;
        for(i=0; i<bins; i++) c[i] += $(i+2);
    }
    END {
        printf "%s\t%d", sample, total;
        for(i=0; i<bins; i++) printf "\t%d", c[i];
        for(i=0; i<bins; i++) {
            if(total>0) printf "\t%.4f", c[i]/total;
            else printf "\t0.0000";
        }
        printf "\n";
    }' "${RAW_DIR}/${SAMPLE}"/*.tsv >> "$DIST_FILE"
done

rm -rf "${TMP_DIR}"

