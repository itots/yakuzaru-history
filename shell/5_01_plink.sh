#!/bin/bash

set -e 
set -o pipefail

# Define variables
INPUT_VCF="../vcf/bial_snp/snp_chr_concat/concatenated.vcf.gz"
OUTPUT_DIR="../stats/data/plink"
CHROM_LIST="../list/chromosome.txt"
JM_LIST="../list/jm.txt"
SAMPLE_SEX="../list/sample_sex_qc_filtered.txt"
CUTOFF_VALUE=0.0884 # remove second degree relations

# Create the output directory if it doesn't exist
mkdir -p "${OUTPUT_DIR}/info"

# Make the list of sex
awk '{
    sex = ($2 == "male") ? 1 : ($2 == "female") ? 2 : 0;
    print "0", $1, sex;
}' "${SAMPLE_SEX}" > "${OUTPUT_DIR}/info/sex.txt"

# Make JM keep file (FID IID)
awk '{print "0", $1}' "${JM_LIST}" > "${OUTPUT_DIR}/info/jm_keep.txt"

mkdir -p "${OUTPUT_DIR}/input"
mkdir -p "${OUTPUT_DIR}/out"

# Convert VCF to PLINK format
plink \
    --vcf "${INPUT_VCF}" \
    --const-fid 0 \
    --set-missing-var-ids @:# \
    --chr-set 20 --allow-extra-chr \
    --make-bed \
    --update-sex "${OUTPUT_DIR}/info/sex.txt" \
    --out "${OUTPUT_DIR}/input/concatenated" \
    &> "${OUTPUT_DIR}/input/concatenated.log"

# Detect and filter related individuals
plink2 \
    --bfile "${OUTPUT_DIR}/input/concatenated" \
    --chr-set 20 --allow-extra-chr \
    --make-king \
    --make-king-table \
    --king-cutoff "${CUTOFF_VALUE}" \
    --out "${OUTPUT_DIR}/out/king"

plink \
    --bfile "${OUTPUT_DIR}/input/concatenated" \
    --chr-set 20 --allow-extra-chr \
    --remove "${OUTPUT_DIR}/out/king.king.cutoff.out.id" \
    --make-bed \
    --out "${OUTPUT_DIR}/out/king_pruned" \
    &> "${OUTPUT_DIR}/out/king_pruned.log"

# Extract Japanese macaques
plink \
    --bfile "${OUTPUT_DIR}/out/king_pruned" \
    --chr-set 20 --allow-extra-chr \
    --keep "${OUTPUT_DIR}/info/jm_keep.txt" \
    --make-bed \
    --out "${OUTPUT_DIR}/out/jm" \
    &> "${OUTPUT_DIR}/out/jm.log"

# LD pruning
plink \
    --bfile "${OUTPUT_DIR}/out/jm" \
    --chr-set 20 --allow-extra-chr \
    --indep-pairwise 50 5 0.2 \
    --out "${OUTPUT_DIR}/out/jm_ld_pruning" \
    &> "${OUTPUT_DIR}/out/jm_ld_pruning.log"

plink \
    --bfile "${OUTPUT_DIR}/out/jm" \
    --chr-set 20 --allow-extra-chr \
    --extract "${OUTPUT_DIR}/out/jm_ld_pruning.prune.in" \
    --make-bed \
    --out "${OUTPUT_DIR}/out/jm_ld_pruned" \
    &> "${OUTPUT_DIR}/out/jm_ld_pruned.log"

# Bp space pruning (all samples)
plink \
    --bfile "${OUTPUT_DIR}/out/king_pruned" \
    --chr-set 20 --allow-extra-chr \
    --geno 0 \
    --make-bed \
    --out "${OUTPUT_DIR}/out/all_geno0" \
    &> "${OUTPUT_DIR}/out/all_geno0.log"

plink \
    --bfile "${OUTPUT_DIR}/out/all_geno0" \
    --chr-set 20 --allow-extra-chr \
    --bp-space 10000 \
    --make-bed \
    --out "${OUTPUT_DIR}/out/all_bp_space" \
    &> "${OUTPUT_DIR}/out/all_bp_space.log"

# IBS distance matrix (ALL samples)
plink \
    --bfile "${OUTPUT_DIR}/out/all_bp_space" \
    --chr-set 20 --allow-extra-chr \
    --autosome \
    --distance 1-ibs \
    --out "${OUTPUT_DIR}/out/ibs" \
    &> "${OUTPUT_DIR}/out/ibs.log"
COM
# PCA (JM only)
plink \
    --bfile "${OUTPUT_DIR}/out/jm_ld_pruned" \
    --chr-set 20 --allow-extra-chr \
    --autosome \
    --pca \
    --out "${OUTPUT_DIR}/out/pca" \
    &> "${OUTPUT_DIR}/out/pca.log"

