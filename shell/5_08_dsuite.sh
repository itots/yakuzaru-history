#!/bin/bash

set -e 
set -o pipefail

# Define variables
INPUT_DIR="../vcf/bial_snp/snp"
OUTPUT_DIR="../stats/dsuite"
SAMPLE_LIST_TREEMIX="../list/sample_treemix.txt"
SAMPLE_LIST_DSUITE="../list/sample_dsuite.txt"
TREE_TREEMIX="../stats/treemix/treemix.0.treeout.gz"
TREE_DSUITE="../stats/data/dsuite/tree.nwk"

mkdir -p "${TREE_DSUITE%/*}"

# Convert to absolute paths
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
to_abs() {
  local p="$1"
  [[ "$p" = /* ]] && { echo "$p"; return; }
  echo "$(cd "$SCRIPT_DIR" && cd "$(dirname "$p")" && pwd)/$(basename "$p")"
}

# Directory setup
INPUT_DIR="$(to_abs "$INPUT_DIR")"
OUTPUT_DIR="$(to_abs "$OUTPUT_DIR")"
SAMPLE_LIST_DSUITE="$(to_abs "$SAMPLE_LIST_DSUITE")"
TREE_DSUITE="$(to_abs "$TREE_DSUITE")"

mkdir -p "$(dirname "$SAMPLE_LIST_DSUITE")"
mkdir -p "$(dirname "$TREE_DSUITE")"

for dir in a x a-x; do
  mkdir -p "${OUTPUT_DIR}/${dir}"
done

# Prepare sample lists for Dsuite
awk -F'\t' 'BEGIN {OFS="\t"} {
    pop_name = ($2 == "fas_Indonesia") ? "Outgroup" : $2;
    print $1, pop_name;    
} ' "${SAMPLE_LIST_TREEMIX}" > "${SAMPLE_LIST_DSUITE}"

zcat "${TREE_TREEMIX}" | head -n 1 | \
    sed -e 's/fas_Indonesia/Outgroup/g' -e 's/:[0-9.eE+-]*//g' | \
    nw_reroot - "Outgroup" > "${TREE_DSUITE}"

for dir in a x; do
   cp "${SAMPLE_LIST_DSUITE}" "${OUTPUT_DIR}/${dir}/sample_dsuite.txt"
done


# Autosome VCFs
autosome_vcfs=()
for file in "${INPUT_DIR}"/*.vcf.gz; do
    if [[ "$(basename "$file")" != "NC_041774.1.vcf.gz" ]]; then
        autosome_vcfs+=("$file")
    fi
done

# Autosome
pushd "${OUTPUT_DIR}/a" >/dev/null

DtriosParallel \
    -n "a" \
    -t "${TREE_DSUITE}" \
    "sample_dsuite.txt" \
    "${autosome_vcfs[@]}" \
    2> "DtriosParallel.log"

Dsuite Fbranch \
    "${TREE_DSUITE}" \
    "DTparallel_sample_dsuite_a_combined_tree.txt" \
    > "Fbranch.txt" \
    2> "Fbranch.log"

dtools_ed.py \
    "Fbranch.txt" \
    "${TREE_DSUITE}" \
    2> "dtools.log"

popd >/dev/null

# X chromosome
pushd "${OUTPUT_DIR}/x" >/dev/null

Dsuite Dtrios \
    -n "x" \
    -t "${TREE_DSUITE}" \
    "${INPUT_DIR}/NC_041774.1.vcf.gz" \
    "sample_dsuite.txt" \
    2> "Dtrios.log"

Dsuite Fbranch \
    "${TREE_DSUITE}" \
    "sample_dsuite_x_tree.txt" \
    > "Fbranch.txt" \
    2> "Fbranch.log"

dtools_ed.py \
    "Fbranch.txt" \
    "${TREE_DSUITE}" \
    2> "dtools.log"

popd >/dev/null

# X-A difference

paste "${OUTPUT_DIR}/a/Fbranch.txt" "${OUTPUT_DIR}/x/Fbranch.txt" | awk '
    NR==1 {
        printf $1; for(i=2;i<=NF/2;i++) printf "\t"$i; print "";
        next
    }
    NR>=2 {
        printf "%s\t%s", $1, $2;
        for(i=3;i<=NF/2;i++) {
            a = $(i);
            x = $(i+NF/2);
            if(a=="nan" || x=="nan") {
                printf "\t%s", "nan";
            } else {
                printf "\t%f", a-x;
            }
        }
        print "";
    }
' > "${OUTPUT_DIR}/a-x/Fbranch.txt"

pushd "${OUTPUT_DIR}/a-x" >/dev/null
dtools_ed.py \
    "Fbranch.txt" \
    "${TREE_DSUITE}" \
    2> "Fbranch.log"

popd >/dev/null
