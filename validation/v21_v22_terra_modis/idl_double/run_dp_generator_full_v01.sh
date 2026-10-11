#!/bin/bash
#
# run_dp_generator_full.sh
#
# PREPARED, NOT YET AUTHORISED TO RUN.  Generate complete Terra MODIS cloud
# LUTs with the unmodified V21 IDL generator and the double-precision DISORT
# backend (idl/generator/call_disort.pro -> DISORT_DP, NSTR 60 as in V21).
#
# The arguments are those of the V21 production create_orac_lut/terra_modis_run
# except the version label, which defaults to 90 (diagnostic double-precision
# reference; never 21).  Output goes only to a new
# results/full_<timestamp>/ directory here; existing LUTs are untouched.
# Expect roughly 3-4 hours per product (see REPORT_idl_double_precision.md).
#
# Usage: bash run_dp_generator_full.sh [pold|p240|pagg ...]   (default: all three)
#        DP_VERSION=90 bash run_dp_generator_full.sh p240

set -u -o pipefail

HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
ROOT="$(cd "$HERE/../../.." && pwd)"
LEGACY="$ROOT/create_orac_lut"
BUILD="$HERE/build"
VERSION="${DP_VERSION:-90}"
TMATRIX=/network/group/aopp/eodg/shared/dubovik_tmatrix/
PRODUCTS=("$@")
[ ${#PRODUCTS[@]} -eq 0 ] && PRODUCTS=(pold p240 pagg)
[ "$VERSION" = 21 ] && { echo "refusing to label double-precision output as version 21" >&2; exit 1; }

declare -A MMFILE=([pold]=liquid-water_old.mm [p240]=liquid-water_240.mm [pagg]=water-ice_agg.mm)
declare -A LUTFILE=([pold]=liquid-water-cloud.lut [p240]=liquid-water-cloud.lut [pagg]=ice-cloud.lut)

RUN_DIR="$HERE/results/full_$(date +%Y%m%dT%H%M%S)"
LOG="$RUN_DIR/run.log"
fail() { echo "DOUBLE-PRECISION GENERATOR FAILED: $*" >&2; exit 1; }

for p in "${PRODUCTS[@]}"; do [ -n "${MMFILE[$p]:-}" ] || fail "unknown product $p"; done
diff -q <(sed -n 11,92p "$LEGACY/call_disort.pro") \
        <(awk '/^pro call_disort/{p=1} p' "$HERE/idl/generator/call_disort.pro" | head -82) > /dev/null \
   || fail "generator shadow no longer contains the V21 call_disort body verbatim"
grep -q "$(/usr/bin/sha256sum "$BUILD/DISORTtoIDL_DP.so" | cut -d' ' -f1)" "$BUILD/BUILD_INFO.txt" \
   || fail "DISORTtoIDL_DP.so does not match BUILD_INFO.txt"
mkdir -p "$RUN_DIR/luts" "$RUN_DIR/input_files" || fail "cannot create $RUN_DIR"
exec > >(tee -a "$LOG") 2>&1

if ! type module > /dev/null 2>&1; then source /etc/profile.d/modules.sh; fi
module load intel-compilers/2022 || fail "module load intel-compilers/2022"
module load idl/890 || fail "module load idl/890"
for entry in "$LEGACY"/input_files/*; do ln -s "$entry" "$RUN_DIR/input_files/$(basename "$entry")"; done
export IDL_DLM_PATH="$BUILD:<IDL_DEFAULT>:+$LEGACY"
export IDL_PATH="$HERE/idl/generator:<IDL_DEFAULT>:+$LEGACY"
unset CREATE_ORAC_CLOUD_LUT_DRIVER IDL_STARTUP ORACLUT_DP_NSTR ORACLUT_DP_CAPTURE_DIR
{
   echo "git_head $(git -C "$ROOT" rev-parse HEAD)"
   echo "version $VERSION products ${PRODUCTS[*]}"
   /usr/bin/sha256sum "$BUILD/DISORTtoIDL_DP.so" "$HERE/idl/generator/call_disort.pro" "$LEGACY"/*.pro
   cat "$BUILD/BUILD_INFO.txt"
} > "$RUN_DIR/provenance.txt"

cd "$RUN_DIR" || fail "cannot cd"
for p in "${PRODUCTS[@]}"; do
   echo "=== $p: ${MMFILE[$p]} (version $VERSION) started $(date --iso-8601=seconds)"
   start=$(date +%s)
   idl -e "create_orac_cloud_lut_wrapper, driver='input_files/driver/terra_modis_cloud.driver', srf_quad=1, mmfile='${MMFILE[$p]}', lutfile='${LUTFILE[$p]}', tmatrix_path='$TMATRIX', version=$VERSION" \
      < /dev/null 2>&1 | tee "$RUN_DIR/idl_$p.txt"
   status=${PIPESTATUS[0]}
   echo "$p elapsed_seconds $(( $(date +%s) - start ))" | tee -a "$RUN_DIR/elapsed.txt"
   [ "$status" -eq 0 ] || fail "$p: IDL exit status $status"
   grep -q -E '% Execution halted|^% Error|failed with code' "$RUN_DIR/idl_$p.txt" && fail "$p: IDL reported an error"
   grep -q 'DISORTTOIDL_DP' "$RUN_DIR/idl_$p.txt" || fail "$p: DISORTtoIDL_DP was not loaded"
   ls "$RUN_DIR"/luts/terra_modis_cloud/*_p${p#p}_v$(printf %02d "$VERSION").nc > /dev/null 2>&1 \
      || fail "$p: expected LUT not written"
done
find "$RUN_DIR/luts" -name '*.nc' | sort | xargs /usr/bin/sha256sum > "$RUN_DIR/outputs.sha256"
echo "DOUBLE-PRECISION GENERATOR SUCCEEDED: $RUN_DIR"
