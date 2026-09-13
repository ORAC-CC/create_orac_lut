#!/bin/bash

#
# setup_oraclut_env.sh
#
# Prepare a clean standalone ORACLUT IDL environment.
#

set -euo pipefail

#
# Directory containing this script
#
MAKE_ORAC_LUT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
export MAKE_ORAC_LUT_DIR

#These directorys are usally subdirectories of  MAKE_ORAC_LUT_DIR but in principle they could be anywhere
export ORAC_LUT_INPUT_ROOT_DIR="${MAKE_ORAC_LUT_DIR}/input_files"
export ORAC_LUT_OUTPUT_ROOT_DIR="${MAKE_ORAC_LUT_DIR}/luts"

#Solar File
export ORAC_LUT_SOLAR_FILE="${ORAC_LUT_INPUT_ROOT_DIR}/sun/Gueymard2018.sssi"

#Where T-Matrix files live
export ORAC_LUT_TMATRIX_DIR='/network/aopp/matin/eodg/shared/dubovik_tmatrix'

#
# Define a minimal IDL search path.
#
export IDL_PATH="${MAKE_ORAC_LUT_DIR}:${MAKE_ORAC_LUT_DIR}/input_files:${MAKE_ORAC_LUT_DIR}/dubovik:<IDL_DEFAULT>"

#
# Define DLM search path.
#
export IDL_DLM_PATH="${MAKE_ORAC_LUT_DIR}/mie/dlm-code:${MAKE_ORAC_LUT_DIR}/disort2/src:<IDL_DEFAULT>"

#
# Diagnostic output
#
echo "========================================"
echo "ORACLUT environment configured"
echo "MAKE_ORAC_LUT_DIR         = ${MAKE_ORAC_LUT_DIR}"
echo "ORAC_LUT_INPUT_ROOT_DIR   = ${ORAC_LUT_INPUT_ROOT_DIR}"
echo "ORAC_LUT_OUTPUT_ROOT_DIR  = ${ORAC_LUT_OUTPUT_ROOT_DIR}"
echo "ORAC_LUT_TMATRIX_DIR      = ${ORAC_LUT_TMATRIX_DIR}"
echo "IDL_PATH                  = ${IDL_PATH}"
echo "IDL_DLM_PATH              = ${IDL_DLM_PATH}"
echo "========================================"

module load intel-compilers/2022