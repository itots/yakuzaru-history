#!/bin/bash

set -e
set -o pipefail

# Define variables
THREADS=60
MUT_RATE="7.7e-9"
GEN_TIME="10.4"

SCRIPT_DIR=$(cd "$(dirname "$0")" && pwd)
PROJECT_ROOT=$(cd "${SCRIPT_DIR}/.." && pwd)
CONTAINER_ROOT="/mnt"

HOST_INPUT_DIR="${PROJECT_ROOT}/stats/data/smc++/smc++_input"
HOST_OUTPUT_DIR="${PROJECT_ROOT}/stats/smc++"

CONTAINER_INPUT_DIR="${CONTAINER_ROOT}/stats/data/smc++/smc++_input"
CONTAINER_OUTPUT_DIR="${CONTAINER_ROOT}/stats/smc++"

mkdir -p "${HOST_OUTPUT_DIR}/estimate"

# Check if input directory exists and contains files
if [ ! -d "${HOST_INPUT_DIR}" ] || ! ls "${HOST_INPUT_DIR}"/*.smc.gz >/dev/null 2>&1; then
    echo "Error: No .smc.gz files found in ${HOST_INPUT_DIR}" >&2
    exit 1
fi

echo "Input files verified."

CONTAINER_FILES="${CONTAINER_INPUT_DIR}/*.smc.gz"

# Docker Execution
echo "Starting SMC++ estimate..."

docker run --rm \
    --user "$(id -u):$(id -g)" \
    -v "${PROJECT_ROOT}:${CONTAINER_ROOT}" \
    --entrypoint sh \
    terhorst/smcpp:latest -c "
        smc++ estimate --cores ${THREADS} \
            --timepoints 10 10000 \
            --regularization-penalty 4 \
            -o ${CONTAINER_OUTPUT_DIR}/estimate \
            ${MUT_RATE} \
            ${CONTAINER_FILES} \
            2> ${CONTAINER_OUTPUT_DIR}/estimate/estimate.log \
        && \
        smc++ plot -c -g ${GEN_TIME} \
            ${CONTAINER_OUTPUT_DIR}/estimate/plot.pdf \
            ${CONTAINER_OUTPUT_DIR}/estimate/model.final.json \
            2> ${CONTAINER_OUTPUT_DIR}/estimate/plot.log
    "
