#!/usr/bin/env bash
# Generate a legacy IDL reference LUT for one compact validation case, capture
# it with a manifest, and compare it with the matching Python product.
#
# Usage (interactive AOPP shell with an IDL licence; tiny cases only):
#
#   scripts/generate_legacy_reference.sh \
#       --formulation aerosol \
#       --driver validation/reference_cases/meteosat10_seviri_aerosol_a79_test_ch09.driver \
#       --case-id aerosol_ir \
#       --python-product validation/generated/aerosol_ir_comparison/python \
#       --reference-dir validation/generated/aerosol_ir_comparison/legacy \
#       --comparison-dir validation/generated/aerosol_ir_comparison
#
# Options:
#   --formulation cloud|aerosol   selects create_orac_cloud_lut_wrapper or
#                                 create_orac_aerosol_lut_wrapper (required)
#   --driver FILE                 legacy driver file (required)
#   --case-id ID                  short label used in file names (required)
#   --python-product PATH         Python NetCDF, or a directory holding exactly one
#   --reference-dir DIR           where the captured reference and logs go
#                                 (default validation/generated/<case-id>_legacy_reference)
#   --comparison-dir DIR          where comparison.json/.md go
#                                 (default validation/generated/<case-id>_comparison)
#   --repeat | --no-repeat        force or suppress the single repeatability run
#                                 (default: repeat only after a DISORT near-singular warning)
#   --no-compare                  stop after capturing the reference
#
# The scientific settings (channels, gas, srf_quad, version, atmosphere) are read
# from the driver file and reported, never re-derived here.
#
# ---------------------------------------------------------------------------
# WHICH LEGACY SOURCE, AND WHY
# ---------------------------------------------------------------------------
# The V2.1 tree imported under create_orac_lut/ contains no cloud or aerosol
# entry points: create_orac_{cloud,aerosol}_lut.pro and their wrappers were
# deleted in the V2.1 tip commit b3c8a91. The only executable path is the
# original working tree, extracted on scratch from the verified preservation
# archive. It must be used as a coherent whole: V2.1 renamed the
# generate_scattering_properties tmatrix keyword (an immediate IDL error), swaps
# the meaning of srf_quad=0 and 1 in load_srfstrarr, and alters the
# surface-pressure grid in load_lutstr, so mixing generations changes the
# science. IDL searches the working directory before IDL_PATH, so the
# generator is run from a clean directory holding no .pro files, and the source
# file of every resolved routine is checked before any calculation starts.
# The extracted tree is read, never written, and never copied into /home.
# ---------------------------------------------------------------------------

set -euo pipefail

REPO_ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd -P)"
LEGACY_ROOT="${REPO_ROOT}/create_orac_lut"

# THE LEGACY SOURCE ROOT. Defined once; every other reference derives from it.
#   create_orac_aerosol_lut.pro          39296 bytes, 2025-11-14 15:40:32 UTC
#   create_orac_aerosol_lut_wrapper.pro  10223 bytes, 2025-11-14 15:40:30 UTC
LEGACY_SRC_DEFAULT=/network/scratch/grainger/oraclut-reconciliation/legacy-reference-source/current_local_create_orac_lut
LEGACY_SRC="${ORACLUT_LEGACY_SRC:-${LEGACY_SRC_DEFAULT}}"
PYTHON="${ORACLUT_PYTHON:-/home/g/grainger/miniforge3/envs/science/bin/python}"

# Legacy .pro files that must be present, and the routines the provenance gate
# resolves. lut_quadrature lives in load_lutstr.pro, read_srfstr and
# bbconstants in input_files/, so those are checked by resolution.
LEGACY_SOURCE_FILES=(
    create_orac_aerosol_lut.pro create_orac_aerosol_lut_wrapper.pro
    create_orac_cloud_lut.pro create_orac_cloud_lut_wrapper.pro
    makerunfile_v2.pro generate_scattering_properties.pro load_inststr.pro
    load_lutstr.pro load_mmdat.pro load_atmstr.pro load_gasstr.pro
    load_srfstrarr.pro load_solar_spectrum.pro setup_disort.pro call_disort.pro
    write_v2_lut.pro create_bwgp.pro read_ri.pro segment.pro legpexp.pro
    quadrature.pro shift_quadrature.pro mie_size_dist_new.pro
    integrate_trapeziodal.pro loadxy.pro loadxyz.pro create_range.pro
    input_files/read_srfstr.pro input_files/bbconstants.pro
)
COHERENT_REQUIRED=(
    create_orac_aerosol_lut create_orac_aerosol_lut_wrapper
    create_orac_cloud_lut create_orac_cloud_lut_wrapper
    generate_scattering_properties load_srfstrarr load_lutstr load_inststr
    load_atmstr load_mmdat load_gasstr load_solar_spectrum setup_disort
    call_disort write_v2_lut create_bwgp read_ri lut_quadrature quadrature
    shift_quadrature mie_size_dist_new integrate_trapeziodal segment legpexp
    read_srfstr bbconstants
)

die() { echo "generate_legacy_reference.sh: $*" >&2; exit 2; }

FORMULATION=""; DRIVER=""; CASE_ID=""; PYTHON_PRODUCT=""; REFERENCE_DIR=""; COMPARISON_DIR=""
REPEAT_MODE=auto; COMPARE=1
while [[ "$#" -gt 0 ]]; do
    case "$1" in
        --formulation) FORMULATION="$2"; shift ;;
        --driver) DRIVER="$2"; shift ;;
        --case-id) CASE_ID="$2"; shift ;;
        --python-product) PYTHON_PRODUCT="$2"; shift ;;
        --reference-dir) REFERENCE_DIR="$2"; shift ;;
        --comparison-dir) COMPARISON_DIR="$2"; shift ;;
        --repeat) REPEAT_MODE=always ;;
        --no-repeat) REPEAT_MODE=never ;;
        --no-compare) COMPARE=0 ;;
        -h|--help) sed -n '2,30p' "${BASH_SOURCE[0]}" | sed 's/^# \{0,1\}//'; exit 0 ;;
        *) die "unknown option: $1" ;;
    esac
    shift
done
[[ -n "${FORMULATION}" && -n "${DRIVER}" && -n "${CASE_ID}" ]] \
    || die "--formulation, --driver and --case-id are required (see --help)"
case "${FORMULATION}" in
    cloud)   WRAPPER=create_orac_cloud_lut_wrapper;   DRIVER_ENV=CREATE_ORAC_CLOUD_LUT_DRIVER ;;
    aerosol) WRAPPER=create_orac_aerosol_lut_wrapper; DRIVER_ENV=CREATE_ORAC_AEROSOL_LUT_DRIVER ;;
    *) die "--formulation must be cloud or aerosol" ;;
esac
case "${REPO_ROOT}" in
    /home/g/grainger/project-oraclut) ;;
    *) die "refusing unexpected repository root: ${REPO_ROOT}" ;;
esac

absolute() { case "$1" in /*) printf '%s\n' "$1" ;; *) printf '%s\n' "${REPO_ROOT}/$1" ;; esac; }
DRIVER="$(absolute "${DRIVER}")"
REFERENCE_DIR="$(absolute "${REFERENCE_DIR:-validation/generated/${CASE_ID}_legacy_reference}")"
COMPARISON_DIR="$(absolute "${COMPARISON_DIR:-validation/generated/${CASE_ID}_comparison}")"
case "${REFERENCE_DIR}" in "${REPO_ROOT}"/*) ;; *) die "reference dir must be inside the repository" ;; esac
case "${COMPARISON_DIR}" in "${REPO_ROOT}"/*) ;; *) die "comparison dir must be inside the repository" ;; esac

# The wrapper derives out_path as luts/<driver stem>, relative to its cwd.
DRIVER_STEM="$(basename "${DRIVER%.*}")"
RUN_DIR="${REFERENCE_DIR}/run"
CANONICAL_DIR="${RUN_DIR}/luts/${DRIVER_STEM}"

require_readable() { [[ -r "$1" ]] || die "preflight failed, missing or unreadable $2: $1"; }

# ---------------------------------------------------------------------------
# Read the case from the driver (five positional values, then key=value).
# ---------------------------------------------------------------------------
require_readable "${DRIVER}" "legacy driver"
mapfile -t DRIVER_LINES < <(sed 's/#.*//' "${DRIVER}" | sed 's/^[[:space:]]*//;s/[[:space:]]*$//' | grep -v '^$')
INST_FILE="${DRIVER_LINES[1]}"; MM_FILE="${DRIVER_LINES[2]}"; LUT_FILE="${DRIVER_LINES[3]}"; ATMOSPHERE="${DRIVER_LINES[4]}"
driver_option() { printf '%s\n' "${DRIVER_LINES[@]}" | sed -n "s/^$1[[:space:]]*=[[:space:]]*//Ip" | head -1; }
CHANNELS="$(driver_option channelid | tr -d '[] ')"; CHANNELS="${CHANNELS:-all}"
SRF_QUAD="$(driver_option srf_quad)"; SRF_QUAD="${SRF_QUAD:-1}"
GAS_FLAG="$(driver_option gas)"; GAS=$([[ "${GAS_FLAG:-0}" == "1" ]] && echo enabled || echo disabled)
NO_RAYLEIGH="$(driver_option no_rayleigh)"; RAYLEIGH=$([[ "${NO_RAYLEIGH:-0}" == "1" ]] && echo disabled || echo enabled)
VERSION="$(driver_option version)"; VERSION="${VERSION:-21}"
[[ "${DRIVER_LINES[0]}" == "input_files" ]] || die "driver input root must be 'input_files' (relative to the run directory)"

# ---------------------------------------------------------------------------
# Preflight
# ---------------------------------------------------------------------------
require_readable "${LEGACY_SRC}" "legacy source tree"
for source_file in "${LEGACY_SOURCE_FILES[@]}"; do
    require_readable "${LEGACY_SRC}/${source_file}" "legacy source file"
done
echo "Preflight: ${#LEGACY_SOURCE_FILES[@]} legacy source files present in ${LEGACY_SRC}"
require_readable "${LEGACY_ROOT}/input_files/inst/${INST_FILE}" "instrument definition"
require_readable "${LEGACY_ROOT}/input_files/microphysics/${MM_FILE}" "microphysics"
require_readable "${LEGACY_ROOT}/input_files/lut/${LUT_FILE}" "LUT grid"
require_readable "${LEGACY_ROOT}/input_files/sun/Gueymard2018.sssi" "solar spectrum"
require_readable "${REPO_ROOT}/mie/dlm-code/mie_dlm_single.dlm" "Mie DLM declaration"
require_readable "${REPO_ROOT}/mie/dlm-code/mie_dlm_single.so" "Mie DLM shared object"
require_readable "${LEGACY_ROOT}/disort2/src/DISORTtoIDL.dlm" "production DISORT DLM declaration"
require_readable "${LEGACY_ROOT}/disort2/src/DISORTtoIDL.so" "production DISORT shared object"
require_readable "${REPO_ROOT}/validation/diagnostics/aerosol_legacy_provenance.pro" "provenance diagnostic"
PLATFORM="$(sed -n 's/^[[:space:]]*platform[[:space:]]*=[[:space:]]*//Ip' "${LEGACY_ROOT}/input_files/inst/${INST_FILE}" | head -1)"
INSTRUMENT="$(sed -n 's/^[[:space:]]*instrument[[:space:]]*=[[:space:]]*//Ip' "${LEGACY_ROOT}/input_files/inst/${INST_FILE}" | head -1)"
if [[ "${CHANNELS}" != "all" ]]; then
    for channel in ${CHANNELS//,/ }; do
        srf="$(sed -n "s/^[[:space:]]*srf[[:space:]]*\[[[:space:]]*${channel}[[:space:]]*\][[:space:]]*=[[:space:]]*//Ip" "${LEGACY_ROOT}/input_files/inst/${INST_FILE}" | head -1)"
        [[ -n "${srf}" ]] || die "channel ${channel} has no SRF entry in ${INST_FILE}"
        require_readable "${LEGACY_ROOT}/input_files/srf/${srf}" "SRF for channel ${channel}"
        if [[ "${GAS}" == "enabled" ]]; then
            require_readable "${LEGACY_ROOT}/input_files/gas/ModtranGasOpd_A${ATMOSPHERE}_${PLATFORM}_${INSTRUMENT}_ch$(printf '%02d' "${channel}").gas" "gas optical depths for channel ${channel}"
        fi
    done
fi

if compgen -G "${REFERENCE_DIR}/*_legacy.nc" >/dev/null; then
    die "a captured reference already exists in ${REFERENCE_DIR} and is never overwritten automatically; move it aside first"
fi
type module >/dev/null 2>&1 || die "the AOPP 'module' command is unavailable; run this in an interactive AOPP shell"

mkdir -p "${REFERENCE_DIR}" "${RUN_DIR}/luts"
ln -sfn "${LEGACY_ROOT}/input_files" "${RUN_DIR}/input_files"
compgen -G "${RUN_DIR}/*.pro" >/dev/null && die "the run directory must contain no .pro files: ${RUN_DIR}"
STDOUT_LOG="${REFERENCE_DIR}/stdout.log"; STDERR_LOG="${REFERENCE_DIR}/stderr.log"
PROVENANCE_LOG="${REFERENCE_DIR}/provenance.log"; MANIFEST="${REFERENCE_DIR}/capture_manifest.json"

# Fingerprint any pre-existing canonical product.
BEFORE_STATE=absent; BEFORE_SHA=""
if compgen -G "${CANONICAL_DIR}/*.nc" >/dev/null; then
    BEFORE_STATE=present
    BEFORE_SHA="$(sha256sum "${CANONICAL_DIR}"/*.nc | cut -d' ' -f1 | sort | tr '\n' ' ')"
fi

# ---------------------------------------------------------------------------
# Environment
# ---------------------------------------------------------------------------
module load intel-compilers/2022
module load idl/890
export IDL_PATH="<IDL_DEFAULT_PATH>:+${LEGACY_SRC}:${REPO_ROOT}/validation/diagnostics"
export IDL_DLM_PATH="<IDL_DEFAULT_DLM>:${REPO_ROOT}/mie/dlm-code:${LEGACY_ROOT}/disort2/src"
unset IDL_STARTUP
export "${DRIVER_ENV}=${DRIVER}"
command -v idl >/dev/null 2>&1 || die "no idl executable on PATH after 'module load idl/890'"
IDL_VERSION="$(idl -e "print, 'IDL ' + !version.release" 2>/dev/null | tail -1 | tr -d '\r')"

echo "============================================================"
echo "LEGACY IDL REFERENCE: ${CASE_ID} (${FORMULATION})"
echo "============================================================"
echo "Host:              $(hostname)"
echo "Date:              $(date -Is)"
echo "IDL version:       ${IDL_VERSION}"
echo "IDL_PATH:          ${IDL_PATH}"
echo "IDL_DLM_PATH:      ${IDL_DLM_PATH}"
echo "Legacy source:     ${LEGACY_SRC}  (read-only)"
echo "Run directory:     ${RUN_DIR}  (no .pro files; input_files symlinked)"
echo "Entry point:       ${WRAPPER} via \$${DRIVER_ENV}"
echo "Driver:            ${DRIVER}"
echo "Instrument:        ${LEGACY_ROOT}/input_files/inst/${INST_FILE}  (${PLATFORM}/${INSTRUMENT})"
echo "Microphysics:      ${LEGACY_ROOT}/input_files/microphysics/${MM_FILE}"
echo "LUT grid:          ${LEGACY_ROOT}/input_files/lut/${LUT_FILE}"
echo "Channels:          ${CHANNELS}"
echo "Atmosphere:        MODTRAN code ${ATMOSPHERE}"
echo "Gas:               ${GAS};  Rayleigh: ${RAYLEIGH};  SRF quadrature: ${SRF_QUAD};  version: ${VERSION}"
echo "Expected output:   ${CANONICAL_DIR}/<legacy name>.nc  (pre-run: ${BEFORE_STATE})"
echo "============================================================"

# ---------------------------------------------------------------------------
# Provenance gate
# ---------------------------------------------------------------------------
( cd "${RUN_DIR}" && idl -e "aerosol_legacy_provenance" ) > "${PROVENANCE_LOG}" 2>&1 \
    || die "the provenance diagnostic failed; see ${PROVENANCE_LOG}"
MIXED_SOURCE_ADVICE="Mixing source generations changes the science: the active V2.1 tree swaps the
meaning of srf_quad=0 and srf_quad=1 in load_srfstrarr, renames the
generate_scattering_properties tmatrix keyword, and alters the surface-pressure
grid in load_lutstr. Check IDL_PATH and that the run directory holds no .pro files."
for routine in "${COHERENT_REQUIRED[@]}"; do
    resolved="$(sed -n "s|^PROVENANCE ${routine} ||p" "${PROVENANCE_LOG}" | head -1)"
    [[ -n "${resolved}" && "${resolved}" != "UNRESOLVED" ]] || die "routine ${routine} did not resolve; see ${PROVENANCE_LOG}"
    case "${resolved}" in
        "${LEGACY_SRC}"/*) ;;
        *) die "routine ${routine} resolved from ${resolved}, not from the selected legacy source tree ${LEGACY_SRC}.
${MIXED_SOURCE_ADVICE}" ;;
    esac
done
CONTAMINATED="$(sed -n 's/^PROVENANCE //p' "${PROVENANCE_LOG}" \
    | awk -v v21="${LEGACY_ROOT}/" 'index($2, v21) == 1 {print "  " $1 " <- " $2}')"
[[ -z "${CONTAMINATED}" ]] || die "these routines resolved from the active V2.1 tree ${LEGACY_ROOT}:
${CONTAMINATED}
${MIXED_SOURCE_ADVICE}"
echo "Routine provenance verified: ${#COHERENT_REQUIRED[@]} required routines resolve from ${LEGACY_SRC};"
echo "                             no routine resolves from the active V2.1 tree ${LEGACY_ROOT}"

# ---------------------------------------------------------------------------
# Run
# ---------------------------------------------------------------------------
run_legacy() {
    local label="$1" out_log="$2" err_log="$3"
    (
        cd "${RUN_DIR}"
        echo "--- ${label}: ${WRAPPER}"
        idl -e "${WRAPPER}"
    ) > >(tee "${out_log}") 2> >(tee "${err_log}" >&2)
}
START_EPOCH="$(date +%s)"
set +e; run_legacy "run 1" "${STDOUT_LOG}" "${STDERR_LOG}"; IDL_STATUS=$?; set -e
END_EPOCH="$(date +%s)"
echo "Legacy run exit status: ${IDL_STATUS}; elapsed $((END_EPOCH - START_EPOCH)) s"

mapfile -t PRODUCTS < <(find "${CANONICAL_DIR}" -maxdepth 1 -type f -name '*.nc' 2>/dev/null | sort)
[[ "${#PRODUCTS[@]}" -eq 1 ]] || die "expected exactly one NetCDF in ${CANONICAL_DIR}, found ${#PRODUCTS[@]}; see ${STDOUT_LOG} and ${STDERR_LOG}"
CANONICAL_NC="${PRODUCTS[0]}"
AFTER_SHA="$(sha256sum "${CANONICAL_NC}" | cut -d' ' -f1)"
if [[ "${BEFORE_STATE}" == "present" && "${BEFORE_SHA}" == *"${AFTER_SHA}"* ]]; then
    die "the canonical product is unchanged (sha256 ${AFTER_SHA}); the legacy run did not write it"
fi
REFERENCE_NC="${REFERENCE_DIR}/$(basename "${CANONICAL_NC%.nc}")_legacy.nc"
cp -- "${CANONICAL_NC}" "${REFERENCE_NC}"
REFERENCE_SHA="$(sha256sum "${REFERENCE_NC}" | cut -d' ' -f1)"; REFERENCE_SIZE="$(stat -c %s "${REFERENCE_NC}")"
[[ "${REFERENCE_SHA}" == "${AFTER_SHA}" ]] || die "copy verification failed for ${REFERENCE_NC}"
[[ -f "${CANONICAL_DIR}/scatfile.sav" ]] && cp -- "${CANONICAL_DIR}/scatfile.sav" "${REFERENCE_DIR}/scatfile.sav"

# ---------------------------------------------------------------------------
# Warnings and the single repeatability run
# ---------------------------------------------------------------------------
SINGULAR_COUNT="$(cat "${STDOUT_LOG}" "${STDERR_LOG}" | grep -c "SGECO says matrix near singular" || true)"; SINGULAR_COUNT="${SINGULAR_COUNT:-0}"
UNDERFLOW_COUNT="$(cat "${STDOUT_LOG}" "${STDERR_LOG}" | grep -c "Floating underflow" || true)"; UNDERFLOW_COUNT="${UNDERFLOW_COUNT:-0}"
REPEAT_JSON="${REFERENCE_DIR}/repeatability.json"
DO_REPEAT=0
case "${REPEAT_MODE}" in always) DO_REPEAT=1 ;; never) DO_REPEAT=0 ;; auto) [[ "${SINGULAR_COUNT}" -gt 0 ]] && DO_REPEAT=1 ;; esac
if [[ "${DO_REPEAT}" -eq 1 ]]; then
    echo "Legacy repeatability check: one repeat run in a fresh IDL process"
    mv -- "${CANONICAL_DIR}" "${CANONICAL_DIR}.run1"
    run_legacy "run 2" "${REFERENCE_DIR}/stdout_run2.log" "${REFERENCE_DIR}/stderr_run2.log" || true
    REPEAT_NC="$(find "${CANONICAL_DIR}" -maxdepth 1 -type f -name '*.nc' 2>/dev/null | head -1)"
    if [[ -n "${REPEAT_NC}" ]]; then
        cp -- "${REPEAT_NC}" "${REFERENCE_DIR}/$(basename "${REPEAT_NC%.nc}")_legacy_run2.nc"
        REPEAT_SHA="$(sha256sum "${REPEAT_NC}" | cut -d' ' -f1)"
        if [[ "${REPEAT_SHA}" == "${REFERENCE_SHA}" ]]; then REPEAT_VERDICT="bitwise identical"
        else
            REPEAT_VERDICT="differs bitwise; see repeatability_comparison"
            PYTHONPATH="${REPO_ROOT}/src${PYTHONPATH:+:${PYTHONPATH}}" "${PYTHON}" -m oraclut.validation.product_comparison \
                "${REFERENCE_NC}" "${REFERENCE_DIR}/$(basename "${REPEAT_NC%.nc}")_legacy_run2.nc" \
                --formulation "${FORMULATION}" --out-dir "${REFERENCE_DIR}/repeatability_comparison" || true
        fi
    else REPEAT_SHA=""; REPEAT_VERDICT="repeat run produced no product"; fi
    printf '{\n  "performed": true,\n  "trigger": "%s",\n  "near_singular_warnings_run1": %s,\n  "run1_sha256": "%s",\n  "run2_sha256": "%s",\n  "verdict": "%s"\n}\n' \
        "${REPEAT_MODE}" "${SINGULAR_COUNT}" "${REFERENCE_SHA}" "${REPEAT_SHA}" "${REPEAT_VERDICT}" > "${REPEAT_JSON}"
else
    printf '{\n  "performed": false,\n  "trigger": "%s",\n  "near_singular_warnings_run1": %s,\n  "verdict": "not performed; no evidence of numerical sensitivity in run 1"\n}\n' \
        "${REPEAT_MODE}" "${SINGULAR_COUNT}" > "${REPEAT_JSON}"
fi

# ---------------------------------------------------------------------------
# Manifest
# ---------------------------------------------------------------------------
PROVENANCE="$(sed -n 's/^PROVENANCE //p' "${PROVENANCE_LOG}" | awk '{printf "%s    \"%s\": \"%s\"", (NR>1 ? ",\n" : ""), $1, $2}')"
CHANNEL_JSON="$([[ "${CHANNELS}" == "all" ]] && echo '"all"' || echo "[${CHANNELS}]")"
cat > "${MANIFEST}" <<EOF
{
  "case_id": "${CASE_ID}",
  "formulation": "${FORMULATION}",
  "timestamp_utc": "$(date -u +%Y-%m-%dT%H:%M:%SZ)",
  "hostname": "$(hostname)",
  "idl_version": "${IDL_VERSION}",
  "legacy_source_tree": "${LEGACY_SRC}",
  "run_directory": "${RUN_DIR}",
  "idl_path": "${IDL_PATH}",
  "idl_dlm_path": "${IDL_DLM_PATH}",
  "entry_point": "${WRAPPER}",
  "driver": "${DRIVER}",
  "instrument_file": "${LEGACY_ROOT}/input_files/inst/${INST_FILE}",
  "microphysics_file": "${LEGACY_ROOT}/input_files/microphysics/${MM_FILE}",
  "lut_definition": "${LEGACY_ROOT}/input_files/lut/${LUT_FILE}",
  "channel_selection": ${CHANNEL_JSON},
  "atmosphere": ${ATMOSPHERE},
  "gas": $([[ "${GAS}" == enabled ]] && echo true || echo false),
  "rayleigh": $([[ "${RAYLEIGH}" == enabled ]] && echo true || echo false),
  "srf_quadrature": ${SRF_QUAD},
  "disort_streams": 60,
  "phase_order": 1000,
  "version": ${VERSION},
  "canonical_output_path": "${CANONICAL_NC}",
  "canonical_state_before_run": "${BEFORE_STATE}",
  "captured_reference_path": "${REFERENCE_NC}",
  "sha256": "${REFERENCE_SHA}",
  "file_size_bytes": ${REFERENCE_SIZE},
  "idl_exit_status": ${IDL_STATUS},
  "elapsed_seconds": $((END_EPOCH - START_EPOCH)),
  "near_singular_warnings": ${SINGULAR_COUNT},
  "floating_underflow_warnings": ${UNDERFLOW_COUNT},
  "stdout_log": "${STDOUT_LOG}",
  "stderr_log": "${STDERR_LOG}",
  "provenance_log": "${PROVENANCE_LOG}",
  "routine_provenance": {
${PROVENANCE}
  }
}
EOF
echo "Captured legacy reference: ${REFERENCE_NC}  sha256 ${REFERENCE_SHA}  (${REFERENCE_SIZE} bytes)"

# ---------------------------------------------------------------------------
# Compare
# ---------------------------------------------------------------------------
if [[ "${COMPARE}" -eq 1 ]]; then
    [[ -n "${PYTHON_PRODUCT}" ]] || die "--python-product is required unless --no-compare is given"
    PYTHON_PRODUCT="$(absolute "${PYTHON_PRODUCT}")"
    if [[ -d "${PYTHON_PRODUCT}" ]]; then
        mapfile -t PY_NC < <(find "${PYTHON_PRODUCT}" -maxdepth 1 -type f -name '*.nc' | sort)
        [[ "${#PY_NC[@]}" -eq 1 ]] || die "expected exactly one Python NetCDF in ${PYTHON_PRODUCT}, found ${#PY_NC[@]}"
        PYTHON_PRODUCT="${PY_NC[0]}"
    fi
    require_readable "${PYTHON_PRODUCT}" "Python product"
    PYTHONPATH="${REPO_ROOT}/src${PYTHONPATH:+:${PYTHONPATH}}" "${PYTHON}" -m oraclut.validation.product_comparison \
        "${REFERENCE_NC}" "${PYTHON_PRODUCT}" --formulation "${FORMULATION}" \
        --legacy-log "${STDOUT_LOG}" --legacy-stderr "${STDERR_LOG}" \
        --repeatability-json "${REPEAT_JSON}" --out-dir "${COMPARISON_DIR}" \
        $([[ "${CHANNELS}" != "all" ]] && echo "--expected-channels ${CHANNELS}")
    echo "Comparison written to ${COMPARISON_DIR}"
fi
