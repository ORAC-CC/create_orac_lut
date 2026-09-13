#!/bin/bash

# Run once to set-up 
# 1) mie dlm code

set -euo pipefail

MAKE_ORAC_LUT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
MIE_DLM_DIR="${MAKE_ORAC_LUT_DIR}/mie/DLM"

echo "========================================"
echo "Building Mie DLM"
echo "MAKE_ORAC_LUT_DIR = ${MAKE_ORAC_LUT_DIR}"
echo "MIE_DLM_DIR       = ${MIE_DLM_DIR}"
echo "========================================"

cd "${MIE_DLM_DIR}"

#
# Load Intel only if this DLM build really uses ifort/ifx.
# Leave commented until confirmed.
#
# module load intel-compilers/2022

make clean
make

echo "========================================"
echo "Built:"
ls -l mie_dlm_single.so
echo "========================================"

Need to setup disort flm as well