#!/bin/bash
#
# run_idl_double_validation.sh
#
# Execute the double-precision DISORT validation through IDL:
#
#   1. Interface test of DISORT_DP (DISORTtoIDL_DP) inside IDL, using the V21
#      p240 channel-6 capture.  Nothing else runs unless it passes.
#   2. For each case A-D, the unmodified V21 generator with the validation
#      shadow of call_disort: the historical single-precision results are
#      used for the LUT (a compact V21 recomputation), and every DISORT call is
#      also solved by DISORT_DP at NSTR = 48, 60, 80, 100 (+ a repeat at 60)
#      with inputs and outputs captured.
#   3. For each case, the same generator with the double-precision backend
#      shadow (DISORT_DP at NSTR 60) writing a small end-to-end LUT labelled
#      version 90 (diagnostic; never a V21 label).
#
# Every run writes into a new results/run_<timestamp>/ directory under this
# validation directory; results/latest points to the last successful run.
# Existing LUTs and earlier results are never touched.
#
# Usage: bash run_idl_double_validation.sh

set -u -o pipefail

HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
ROOT="$(cd "$HERE/../../.." && pwd)"
LEGACY="$ROOT/create_orac_lut"
BUILD="$HERE/build"
RESULTS="$HERE/results"
CAPTURE_V21="$HERE/../disort_boundary/results/v21_capture_20261010T175154"
NSTR_LIST="48,60,80,100"
CASES=(A B C D)

STAMP="$(date +%Y%m%dT%H%M%S)"
RUN_DIR="$RESULTS/run_$STAMP"
LOG="$RUN_DIR/run.log"

fail() {
   echo "" >&2
   echo "DOUBLE-PRECISION IDL VALIDATION FAILED: $*" >&2
   [ -f "$LOG" ] && echo "Log: $LOG" >&2
   exit 1
}

# ---------------------------------------------------------------------------
# Pre-flight (no IDL)
# ---------------------------------------------------------------------------
for f in "$BUILD/DISORTtoIDL_DP.so" "$BUILD/DISORTtoIDL_DP.dlm" "$BUILD/BUILD_INFO.txt" \
         "$HERE/idl/dp_interface_test.pro" "$HERE/idl/validation/call_disort.pro" \
         "$HERE/idl/generator/call_disort.pro" "$LEGACY/call_disort.pro" \
         "$LEGACY/create_orac_cloud_lut.pro" "$LEGACY/disort2/src/DISORTtoIDL.so" \
         "$CAPTURE_V21/v21_disort_boundary_input.sav" "$CAPTURE_V21/v21_disort_boundary_output.sav"; do
   [ -s "$f" ] || fail "required file missing or empty: $f"
done
for c in "${CASES[@]}"; do
   [ -s "$HERE/cases/case_$c.lut" ] && [ -s "$HERE/cases/case_$c.driver" ] || fail "case $c inputs missing"
done
# The shadows must contain the V21 call_disort.pro lines 11-92 verbatim.
for shadow in validation generator; do
   diff -q <(sed -n 11,92p "$LEGACY/call_disort.pro") \
           <(awk '/^pro call_disort/{p=1} p' "$HERE/idl/$shadow/call_disort.pro" | head -82) > /dev/null \
      || fail "$shadow shadow no longer contains the V21 call_disort body verbatim"
done
# The DLM must be the one recorded in BUILD_INFO.txt.
grep -q "$(sha256sum "$BUILD/DISORTtoIDL_DP.so" | cut -d' ' -f1)" "$BUILD/BUILD_INFO.txt" \
   || fail "DISORTtoIDL_DP.so does not match BUILD_INFO.txt; rebuild with build_dlm.sh"
mkdir -p "$RESULTS"
[ -e "$RUN_DIR" ] && fail "run directory exists: $RUN_DIR"
mkdir "$RUN_DIR" || fail "cannot create $RUN_DIR"
exec > >(tee -a "$LOG") 2>&1

echo "ORACLUT double-precision DISORT validation through IDL"
echo "started : $(date --iso-8601=seconds)"
echo "host    : $(hostname)"
echo "run dir : $RUN_DIR"

if ! type module > /dev/null 2>&1; then
   [ -f /etc/profile.d/modules.sh ] && source /etc/profile.d/modules.sh
fi
type module > /dev/null 2>&1 || fail "environment modules are not available"
module load intel-compilers/2022 || fail "module load intel-compilers/2022 failed"
module load idl/890 || fail "module load idl/890 failed"
command -v idl > /dev/null || fail "idl not found after module load"
echo "idl     : $(command -v idl)"
echo "modules : $(module list -t 2>&1 | tr '\n' ' ')"
for so in "$BUILD/DISORTtoIDL_DP.so" "$LEGACY/disort2/src/DISORTtoIDL.so" "$LEGACY/mie/dlm-code/mie_dlm_single.so"; do
   missing=$(ldd "$so" | grep 'not found')
   [ -z "$missing" ] || fail "unresolved libraries for $so: $missing"
done

export IDL_DLM_PATH="$BUILD:<IDL_DEFAULT>:+$LEGACY"
unset CREATE_ORAC_CLOUD_LUT_DRIVER IDL_STARTUP ORACLUT_DP_CAPTURE_DIR ORACLUT_DP_NSTR_LIST ORACLUT_DP_NSTR
echo "IDL_DLM_PATH: $IDL_DLM_PATH"

{
   echo "git_head $(git -C "$ROOT" rev-parse HEAD 2>/dev/null)"
   echo "git_dirty_legacy $(git -C "$ROOT" status --porcelain --untracked-files=no -- create_orac_lut | tr '\n' ';')"
   echo "nstr_list $NSTR_LIST"
   sha256sum "$BUILD/DISORTtoIDL_DP.so" "$BUILD/DISORTtoIDL_DP.dlm" "$LEGACY/disort2/src/DISORTtoIDL.so" \
      "$HERE/idl/dp_interface_test.pro" "$HERE/idl/validation/call_disort.pro" "$HERE/idl/generator/call_disort.pro" \
      "$LEGACY/call_disort.pro" "$LEGACY/create_orac_cloud_lut.pro" "$LEGACY/create_orac_cloud_lut_wrapper.pro" \
      "$LEGACY/setup_disort.pro" "$LEGACY/generate_scattering_properties.pro" "$HERE"/cases/case_* \
      "$LEGACY/input_files/inst/terra_modis_v1.inst" "$LEGACY/input_files/microphysics/liquid-water_240.mm"
   cat "$BUILD/BUILD_INFO.txt"
} > "$RUN_DIR/provenance.txt"

stage_inputs() {   # $1 = run directory, $2 = case
   mkdir -p "$1/input_files/lut" "$1/luts"
   for entry in "$LEGACY"/input_files/*; do
      name="$(basename "$entry")"
      [ "$name" = lut ] || [ "$name" = driver ] && continue
      ln -s "$entry" "$1/input_files/$name"
   done
   cp -p "$HERE/cases/case_$2.lut" "$1/input_files/lut/"
   cp -p "$HERE/cases/case_$2.driver" "$1/"
}

expected_calls() {   # diffuse + direct per SZA (+ emission for thermal channels)
   local lut="$HERE/cases/case_$1.lut" channel ntau nre nsza thermal
   channel=$(sed -n 's/^channelid=\[\([0-9]*\)\]/\1/p' "$HERE/cases/case_$1.driver")
   ntau=$(sed -n 1p "$lut" | awk '{print $1}'); nre=$(sed -n 3p "$lut" | awk '{print $1}'); nsza=$(sed -n 5p "$lut" | awk '{print $1}')
   thermal=0; [ "$channel" -ge 20 ] && thermal=1
   echo $(( ntau * nre * (1 + nsza + thermal) ))
}

# ---------------------------------------------------------------------------
# 1. Interface test
# ---------------------------------------------------------------------------
echo ""
echo "=== Interface test ==="
mkdir -p "$RUN_DIR/interface"
cd "$RUN_DIR/interface" || fail "cannot cd"
IDL_PATH="$HERE/idl:<IDL_DEFAULT>:+$LEGACY" idl -e \
   "dp_interface_test, '$CAPTURE_V21/v21_disort_boundary_input.sav', '$CAPTURE_V21/v21_disort_boundary_output.sav', '$RUN_DIR/interface/dp_interface_test.sav'" \
   < /dev/null 2>&1 | tee "$RUN_DIR/interface/idl_output.txt"
grep -q 'DP INTERFACE TEST PASSED' "$RUN_DIR/interface/idl_output.txt" || fail "interface test did not pass"

# ---------------------------------------------------------------------------
# 2./3. Cases
# ---------------------------------------------------------------------------
for c in "${CASES[@]}"; do
   echo ""
   echo "=== Case $c: validation (V21 single precision used, DISORT_DP captured) ==="
   dir="$RUN_DIR/case_$c/validation"
   stage_inputs "$dir" "$c"
   mkdir -p "$dir/captures"
   cd "$dir" || fail "cannot cd"
   start=$(date +%s)
   ORACLUT_DP_CAPTURE_DIR="$dir/captures" ORACLUT_DP_NSTR_LIST="$NSTR_LIST" \
   IDL_PATH="$HERE/idl/validation:<IDL_DEFAULT>:+$LEGACY" \
      idl -e "create_orac_cloud_lut_wrapper, driver='case_$c.driver'" < /dev/null 2>&1 | tee "$dir/idl_output.txt"
   status=${PIPESTATUS[0]}
   echo "elapsed_seconds $(( $(date +%s) - start ))" > "$dir/elapsed.txt"
   [ "$status" -eq 0 ] || fail "case $c validation: IDL exit status $status"
   grep -q -E '% Execution halted|^% Error|failed with code' "$dir/idl_output.txt" && fail "case $c validation: IDL reported an error"
   [ -s "$dir/luts/case_$c/timestamp.txt" ] || fail "case $c validation: generator did not complete"
   ls "$dir"/luts/case_$c/*_v21.nc > /dev/null 2>&1 || fail "case $c validation: no LUT written"
   ncalls=$(ls "$dir"/captures/call_*.sav 2>/dev/null | wc -l)
   expected=$(expected_calls "$c")
   [ "$ncalls" -eq "$expected" ] || fail "case $c validation: $ncalls captured calls, $expected expected"
   echo "case $c validation: $ncalls DISORT calls captured"

   echo ""
   echo "=== Case $c: double-precision generator (version 90, NSTR 60) ==="
   dir="$RUN_DIR/case_$c/dp_generator"
   stage_inputs "$dir" "$c"
   cd "$dir" || fail "cannot cd"
   start=$(date +%s)
   IDL_PATH="$HERE/idl/generator:<IDL_DEFAULT>:+$LEGACY" \
      idl -e "create_orac_cloud_lut_wrapper, driver='case_$c.driver', version=90" < /dev/null 2>&1 | tee "$dir/idl_output.txt"
   status=${PIPESTATUS[0]}
   echo "elapsed_seconds $(( $(date +%s) - start ))" > "$dir/elapsed.txt"
   [ "$status" -eq 0 ] || fail "case $c generator: IDL exit status $status"
   grep -q -E '% Execution halted|^% Error|failed with code' "$dir/idl_output.txt" && fail "case $c generator: IDL reported an error"
   [ -s "$dir/luts/case_$c/timestamp.txt" ] || fail "case $c generator: generator did not complete"
   ls "$dir"/luts/case_$c/*_v90.nc > /dev/null 2>&1 || fail "case $c generator: no version-90 LUT written"
   grep -q 'DISORT_DP\|DISORTTOIDL_DP' "$dir/idl_output.txt" || fail "case $c generator: DISORTtoIDL_DP was not loaded"
done

find "$RUN_DIR" -name '*.nc' -o -name '*.sav' | sort | xargs sha256sum > "$RUN_DIR/outputs.sha256"
ln -sfn "$(basename "$RUN_DIR")" "$RESULTS/latest"
echo ""
echo "DOUBLE-PRECISION IDL VALIDATION RUN SUCCEEDED"
echo "results : $RUN_DIR"
echo "log     : $LOG"
exit 0
