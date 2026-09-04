#!/bin/bash

set -e 
set -o pipefail

# Define variables
INPUT_DIR="../stats/data/msmc/msmc_input"
OUTPUT_DIR="../stats/msmc"
MSMC_TOOL="../python/msmc-tools/combineCrossCoal.py"
PARALLEL_JOBS=6
THREADS=10

mkdir -p "${OUTPUT_DIR}/pop"
mkdir -p "${OUTPUT_DIR}/div"
mkdir -p "${OUTPUT_DIR}/combined_cc"

run_pop () {
    local run_name=$1
    local time_segments=$2

    echo "Processing population: ${run_name}"

    msmc2_Linux -t "${THREADS}" -o "${OUTPUT_DIR}/pop/${run_name}" \
        -s \
        "${INPUT_DIR}/pop/${run_name}/"*.msmc \
        -p "${time_segments}"
}

run_div () {
    local run_name=$1
    local time_segments=$2

    echo "Processing divergence: ${run_name}"

    msmc2_Linux -t "${THREADS}" -o "${OUTPUT_DIR}/div/${run_name}" \
        -s \
        -I 0-4,0-5,0-6,0-7,1-4,1-5,1-6,1-7,2-4,2-5,2-6,2-7,3-4,3-5,3-6,3-7 \
        "${INPUT_DIR}/div/${run_name}/"*.msmc \
        -p "${time_segments}"
}

export -f run_pop
export -f run_div
export INPUT_DIR OUTPUT_DIR THREADS

TIME_SEGMENTS="30*1+1*2+1*4+1*6"

export TIME_SEGMENTS

run_pop_wrapper () {
    local run_name=$1
    run_pop "$run_name" "${TIME_SEGMENTS}"
}

run_div_wrapper () {
    local run_name=$1
    run_div "$run_name" "${TIME_SEGMENTS}"
}

export -f run_pop_wrapper
export -f run_div_wrapper

# Population runs (excluding special cases)
POP_DIRS=$(find "${INPUT_DIR}/pop" -mindepth 1 -maxdepth 1 -type d -exec basename {} \; \
    | grep -Ev "Yakushima_5samples" | sort)

parallel --halt soon,fail=1 -j"${PARALLEL_JOBS}" run_pop_wrapper {} ::: $POP_DIRS

# Divergence runs
DIV_DIRS=$(find "${INPUT_DIR}/div" -mindepth 1 -maxdepth 1 -type d -exec basename {} \; \
    | sort)

parallel --halt soon,fail=1 -j"${PARALLEL_JOBS}" run_div_wrapper {} ::: $DIV_DIRS

# Special case: Yakushima 5 samples
THREADS=60
export THREADS

run_pop "Yakushima_5samples" "35*1+1*2+1*4+1*6"

# --- Combine cross-coalescence ---
combine_cc () {
    local pair=$1

    local log_file
    log_file=$(ls "${INPUT_DIR}/div/${pair}/"*.log | head -1)

    local pop1 pop2
    pop1=$(sed -n '1p' "$log_file" | cut -d: -f1)
    pop2=$(sed -n '2p' "$log_file" | cut -d: -f1)

    local pop1_file="${OUTPUT_DIR}/pop/${pop1}_2samples.final.txt"
    local pop2_file="${OUTPUT_DIR}/pop/${pop2}_2samples.final.txt"
    local div_file="${OUTPUT_DIR}/div/${pair}.final.txt"

    for f in "$div_file" "$pop1_file" "$pop2_file"; do
        if [ ! -f "$f" ]; then
            echo "Warning: file not found: ${f}" >&2
            return 0
        fi
    done

    echo "Combining cross-coalescence: ${pair} (${pop1} x ${pop2})"

    python "${MSMC_TOOL}" \
        "${div_file}" \
        "${pop1_file}" \
        "${pop2_file}" \
        > "${OUTPUT_DIR}/combined_cc/${pair}.txt" \
        2> "${OUTPUT_DIR}/combined_cc/${pair}.log"
}

export -f combine_cc
export MSMC_TOOL

find "${INPUT_DIR}/div" -mindepth 1 -maxdepth 1 -type d -exec basename {} \; | sort | \
    parallel --halt soon,fail=1 -j"${PARALLEL_JOBS}" combine_cc {}
