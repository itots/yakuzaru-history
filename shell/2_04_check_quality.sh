#!/bin/bash

set -e
set -o pipefail

# --- Settings ---
INPUT_DIR="../vcf/joint"
OUT_DIR="../stats/sample_quality"
PARALLEL_JOBS=20

mkdir -p "${OUT_DIR}"
STATS_DIR="${OUT_DIR}/chrom_stats/autosomes"
mkdir -p "${STATS_DIR}"

OUT_SUMMARY="${OUT_DIR}/summary.tsv"

# Statistics calculation function
run_stats() {
    local vcf_file=$1
    local out_dir=$2
    local filename=$(basename "${vcf_file}")

    echo "Processing ${filename}"
    
    # Process only autosomes (Variant sites only)
    bcftools view -v snps,indels "${vcf_file}" | \
    bcftools stats -s - > "${out_dir}/${filename}.stats"
}
export -f run_stats

echo "Starting Quality Check..."

# Find autosomal VCFs and run stats in parallel
find "${INPUT_DIR}" -name "NC_*.vcf.gz" | \
    grep -Ev "NC_005943.1|NC_041774.1|NC_027914.1" | \
    sort | \
    parallel -j "${PARALLEL_JOBS}" run_stats {} "${STATS_DIR}"

# --- Aggregation ---
echo "Aggregating statistics..."

# Header
echo -e "Sample\tnVariants\tnNonRefHom\tnHets\tnMissing\tMissing_Rate(%)\tTsTv_Ratio\tWeighted_Avg_Depth\tnSingleton" > "${OUT_SUMMARY}"

grep -h "^PSC" "${STATS_DIR}"/*.stats | mawk '
BEGIN { OFS="\t" }
{
    # PSC Format:
    # [3]sample, [4]nRefHom, [5]nNonRefHom, [6]nHets, [7]nTs, [8]nTv, [10]avgDepth, [11]nSingleton, [14]nMissing
    
    sample = $3
    
    nSitesInChrom = $4 + $5 + $6 + $14 
    
    sum_NonRefHom[sample] += $5
    sum_Hets[sample] += $6
    sum_Ts[sample] += $7
    sum_Tv[sample] += $8
    sum_Missing[sample] += $14
    sum_Singleton[sample] += $11
    
    # Weighted Depth Calculation: Sum(AvgDepth * NumSites)
    sum_TotalDepthWeight[sample] += ($10 * nSitesInChrom)
    sum_TotalSites[sample] += nSitesInChrom
}
END {
    for (s in sum_TotalSites) {
        total = sum_TotalSites[s]
        
        # Missing Rate
        if (total > 0) rate = (sum_Missing[s] / total) * 100
        else rate = 0
        
        # Ts/Tv Ratio
        if (sum_Tv[s] > 0) tstv = sum_Ts[s] / sum_Tv[s]
        else tstv = 0
        
        # Weighted Average Depth
        if (total > 0) w_depth = sum_TotalDepthWeight[s] / total
        else w_depth = 0
        
        printf "%s\t%d\t%d\t%d\t%d\t%.4f\t%.4f\t%.2f\t%d\n", \
            s, total, sum_NonRefHom[s], sum_Hets[s], sum_Missing[s], rate, tstv, w_depth, sum_Singleton[s]
    }
}' | sort -k1,1 >> "${OUT_SUMMARY}"

echo "Done. Summary saved to ${OUT_SUMMARY}"

