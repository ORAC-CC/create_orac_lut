#!/bin/bash
#
# run_v21_capture.sh
#
# Capture the V21 IDL DISORT boundary (inputs and outputs) for the selected
# Terra MODIS p240 direct-beam call: channel 6, liquid-water_240.mm,
# optical depth 256, effective radius 29 um, solar zenith 80 deg, with the
# full 10 satellite-zenith x 11 relative-azimuth viewing grid of that call.
#
# The legacy V21 IDL sources are used unmodified.  Only call_disort is shadowed
# by capture_call_disort.pro, whose numerical statements are identical to
# create_orac_lut/call_disort.pro; it saves the DISORT arguments before and the
# results after the selected call, then exits before any LUT can be written.
#
# Every run writes into a new results/v21_capture_<timestamp>/ directory, so
# earlier captures are never overwritten.  results/v21_capture_latest points to
# the most recent successful capture.
#
# Usage: bash run_v21_capture.sh

set -u -o pipefail

CASE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
ROOT="$(cd "$CASE/../../.." && pwd)"
LEGACY="$ROOT/create_orac_lut"
RESULTS="$CASE/results"
SHADOW_SOURCE="$CASE/capture_call_disort.pro"
DRIVER_SOURCE="$CASE/p240_ch06_node.driver"
INPUTS="$CASE/run/input_files"

STAMP="$(date +%Y%m%dT%H%M%S)"
RUN_DIR="$RESULTS/v21_capture_$STAMP"
LOG="$RUN_DIR/v21_capture.log"
INPUT_SAV="$RUN_DIR/v21_disort_boundary_input.sav"
OUTPUT_SAV="$RUN_DIR/v21_disort_boundary_output.sav"

fail() {
   echo "" >&2
   echo "V21 DISORT CAPTURE FAILED: $*" >&2
   [ -f "$LOG" ] && echo "Diagnostic log: $LOG" >&2
   exit 1
}

# ---------------------------------------------------------------------------
# Pre-flight checks (no IDL yet)
# ---------------------------------------------------------------------------
for f in "$SHADOW_SOURCE" "$DRIVER_SOURCE" \
         "$LEGACY/call_disort.pro" "$LEGACY/create_orac_cloud_lut.pro" \
         "$LEGACY/create_orac_cloud_lut_wrapper.pro" \
         "$LEGACY/disort2/src/DISORTtoIDL.dlm" "$LEGACY/disort2/src/DISORTtoIDL.so" \
         "$LEGACY/mie/dlm-code/mie_dlm_single.dlm" "$LEGACY/mie/dlm-code/mie_dlm_single.so" \
         "$INPUTS/inst/terra_modis_v1.inst" "$INPUTS/microphysics/liquid-water_240.mm" \
         "$INPUTS/lut/compact_p240_node.lut"; do
   [ -s "$f" ] || fail "required file missing or empty: $f"
done
[ -d "$RESULTS" ] || fail "results directory missing: $RESULTS"
[ -e "$RUN_DIR" ] && fail "run directory already exists: $RUN_DIR"

mkdir "$RUN_DIR" || fail "cannot create $RUN_DIR"
exec > >(tee -a "$LOG") 2>&1

echo "ORACLUT V21 DISORT boundary capture"
echo "started     : $(date --iso-8601=seconds)"
echo "host        : $(hostname)"
echo "run dir     : $RUN_DIR"

# ---------------------------------------------------------------------------
# Modules: the same Intel runtime as the V21 production run (terra_modis_run)
# and IDL 8.9.  DISORTtoIDL.so needs libifcore etc. from intel-compilers/2022.
# ---------------------------------------------------------------------------
if ! type module > /dev/null 2>&1; then
   [ -f /etc/profile.d/modules.sh ] && source /etc/profile.d/modules.sh
fi
type module > /dev/null 2>&1 || fail "environment modules are not available"
module load intel-compilers/2022 || fail "module load intel-compilers/2022 failed"
module load idl/890 || fail "module load idl/890 failed"
command -v idl > /dev/null || fail "idl not found after module load"
echo "idl         : $(command -v idl)"
echo "modules     : $(module list -t 2>&1 | tr '\n' ' ')"

missing=$(ldd "$LEGACY/disort2/src/DISORTtoIDL.so" "$LEGACY/mie/dlm-code/mie_dlm_single.so" | grep 'not found')
[ -z "$missing" ] || fail "unresolved DLM libraries: $missing"

# ---------------------------------------------------------------------------
# Run directory: the wrapper resolves input_files/ and luts/<driver name>/
# relative to the working directory, so everything IDL writes stays here.
# ---------------------------------------------------------------------------
mkdir "$RUN_DIR/idl_shadow" || fail "cannot create idl_shadow"
cp -p "$SHADOW_SOURCE" "$RUN_DIR/idl_shadow/call_disort.pro" || fail "cannot stage call_disort shadow"
cp -p "$DRIVER_SOURCE" "$RUN_DIR/p240_ch06_node.driver" || fail "cannot stage driver"
ln -s "$INPUTS" "$RUN_DIR/input_files" || fail "cannot link input_files"

# Same search paths as the V21 run (create_orac_lut/initpath.bash) with only
# the call_disort shadow placed in front.
export IDL_PATH="$RUN_DIR/idl_shadow:<IDL_DEFAULT>:+$LEGACY"
export IDL_DLM_PATH="<IDL_DEFAULT>:+$LEGACY"
export ORACLUT_V21_DISORT_INPUT="$INPUT_SAV"
export ORACLUT_V21_DISORT_OUTPUT="$OUTPUT_SAV"
unset CREATE_ORAC_CLOUD_LUT_DRIVER

echo "IDL_PATH    : $IDL_PATH"
echo "IDL_DLM_PATH: $IDL_DLM_PATH"
echo "IDL_STARTUP : ${IDL_STARTUP:-<unset>}"

# Provenance of the sources actually executed
{
   echo "git_head $(git -C "$ROOT" rev-parse HEAD 2>/dev/null)"
   echo "git_dirty_legacy $(git -C "$ROOT" status --porcelain --untracked-files=no -- create_orac_lut | tr '\n' ';')"
   sha256sum "$SHADOW_SOURCE" "$RUN_DIR/idl_shadow/call_disort.pro" "$DRIVER_SOURCE" \
      "$LEGACY/call_disort.pro" "$LEGACY/create_orac_cloud_lut.pro" \
      "$LEGACY/create_orac_cloud_lut_wrapper.pro" "$LEGACY/setup_disort.pro" \
      "$LEGACY/generate_scattering_properties.pro" \
      "$LEGACY/disort2/src/DISORTtoIDL.so" "$LEGACY/mie/dlm-code/mie_dlm_single.so" \
      "$INPUTS/lut/compact_p240_node.lut" "$INPUTS/inst/terra_modis_v1.inst" \
      "$INPUTS/microphysics/liquid-water_240.mm"
} > "$RUN_DIR/provenance.txt"

# ---------------------------------------------------------------------------
# Run IDL
# ---------------------------------------------------------------------------
echo ""
echo "Running V21 IDL (this computes Mie properties for one radius first) ..."
cd "$RUN_DIR" || fail "cannot cd to $RUN_DIR"
idl -e "create_orac_cloud_lut_wrapper, driver='p240_ch06_node.driver'" < /dev/null 2>&1 | tee "$RUN_DIR/idl_output.txt"
status=${PIPESTATUS[0]}
echo "idl exit status: $status"

# ---------------------------------------------------------------------------
# Verify
# ---------------------------------------------------------------------------
[ $status -eq 0 ] || fail "IDL exited with status $status"
[ -s "$INPUT_SAV" ] || fail "input capture missing or empty: $INPUT_SAV"
[ -s "$OUTPUT_SAV" ] || fail "output capture missing or empty: $OUTPUT_SAV"
grep -q 'ORACLUT V21 DISORT boundary capture complete' "$RUN_DIR/idl_output.txt" || fail "capture completion message not found in log"
written=$(find "$RUN_DIR/luts" -name '*.nc' 2>/dev/null)
[ -z "$written" ] || fail "a LUT was written, so the selected call was not intercepted: $written"

sha256sum "$INPUT_SAV" "$OUTPUT_SAV" >> "$RUN_DIR/provenance.txt"
ln -sfn "$(basename "$RUN_DIR")" "$RESULTS/v21_capture_latest"

echo ""
echo "V21 DISORT CAPTURE SUCCEEDED"
echo "input  : $INPUT_SAV ($(stat -c %s "$INPUT_SAV") bytes)"
echo "output : $OUTPUT_SAV ($(stat -c %s "$OUTPUT_SAV") bytes)"
echo "log    : $LOG"
echo "latest : $RESULTS/v21_capture_latest"
exit 0
