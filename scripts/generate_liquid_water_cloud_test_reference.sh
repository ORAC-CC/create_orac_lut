#!/usr/bin/env bash
set -euo pipefail

# The outer runner may be launched from any directory. The IDL process itself
# must run from LEGACY_ROOT because the historical source uses relative paths
# such as input_files/... and luts/....
REPO_ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd -P)"
LEGACY_ROOT="${REPO_ROOT}/create_orac_lut"
CASE_ID="meteosat10_seviri_liquid_water_stg_cloud_test"
OUTPUT_BASE="${REPO_ROOT}/validation/generated"
FINAL_ROOT="${OUTPUT_BASE}/${CASE_ID}"
STAGING_ROOT="${OUTPUT_BASE}/.${CASE_ID}.tmp"
DRIVER="${LEGACY_ROOT}/input_files/driver/meteosat-10_seviri_cloud.driver"
TEST_LUT="${LEGACY_ROOT}/input_files/lut/liquid-water-cloud_test.lut"
MODEL="${LEGACY_ROOT}/input_files/microphysics/liquid-water_stg.mm"
INSTRUMENT="${LEGACY_ROOT}/input_files/inst/meteosat-10_seviri_v1.inst"
SRF="${LEGACY_ROOT}/input_files/srf/rtcoef_msg_3_seviri_srf_ch01.txt"
RI="${LEGACY_ROOT}/input_files/ri/H2O_Segelstein_1981.ri"
ATMOSPHERE="${LEGACY_ROOT}/input_files/atm/mls.atm"
SOLAR="${LEGACY_ROOT}/input_files/sun/Gueymard2018.sssi"
MIE_DLM="${REPO_ROOT}/mie/dlm-code"
DISORT_DLM="${LEGACY_ROOT}/disort2/src"
CANONICAL_OUTPUT_DIR="${LEGACY_ROOT}/luts/meteosat-10_seviri_cloud"
SOURCE_NC="${CANONICAL_OUTPUT_DIR}/meteosat-10_seviri_m_liquid-water_a01_pstg_v21.nc"
CAPTURE_FILENAME="${CASE_ID}_legacy_reference_v21.nc"

case "${REPO_ROOT}" in
    /home/g/grainger/project-oraclut) ;;
    *) echo "Refusing unexpected repository root: ${REPO_ROOT}" >&2; exit 2 ;;
esac
case "${OUTPUT_BASE}" in
    "${REPO_ROOT}"/*) ;;
    *) echo "Refusing output outside repository: ${OUTPUT_BASE}" >&2; exit 2 ;;
esac

require_readable() {
    if [[ ! -r "$1" ]]; then
        echo "Preflight failed: missing or unreadable $2: $1" >&2
        exit 2
    fi
}

echo "Preflight: current legacy cloud/STG case"
require_readable "${DRIVER}" "Meteosat-10 cloud driver"
require_readable "${INSTRUMENT}" "Meteosat-10 SEVIRI instrument"
require_readable "${TEST_LUT}" "liquid-water-cloud_test LUT"
require_readable "${MODEL}" "STG microphysics"
require_readable "${SRF}" "SEVIRI channel-1 SRF"
require_readable "${RI}" "liquid-water refractive index"
require_readable "${ATMOSPHERE}" "atmosphere code 2"
require_readable "${SOLAR}" "solar spectrum"
require_readable "${MIE_DLM}/mie_dlm_single.dlm" "Mie DLM declaration"
require_readable "${MIE_DLM}/mie_dlm_single.so" "Mie DLM shared object"
require_readable "${DISORT_DLM}/DISORTtoIDL.dlm" "production DISORT DLM declaration"
require_readable "${DISORT_DLM}/DISORTtoIDL.so" "production DISORT shared object"

if ! grep -qiE '^platform[[:space:]]*=[[:space:]]*meteosat-10[[:space:]]*$' "${INSTRUMENT}"; then
    echo "Preflight failed: instrument is not Meteosat-10" >&2; exit 2
fi
if ! grep -qiE '^instrument[[:space:]]*=[[:space:]]*seviri[[:space:]]*$' "${INSTRUMENT}"; then
    echo "Preflight failed: instrument is not SEVIRI" >&2; exit 2
fi
if ! grep -qE '^2[[:space:]]*(#.*)?$' "${DRIVER}"; then
    echo "Preflight failed: driver does not select atmosphere 2" >&2; exit 2
fi

echo "  legacy cwd: ${LEGACY_ROOT}"
echo "  driver: ${DRIVER}"
echo "  model override: liquid-water_stg.mm"
echo "  LUT override: liquid-water-cloud_test.lut"
echo "  channel override: [1]"
echo "  SRF mode: 1; Rayleigh: enabled; gas: disabled; DISORT streams: 60"
echo "  canonical source: ${SOURCE_NC}"
echo "  validation capture: ${FINAL_ROOT}/${CAPTURE_FILENAME}"

if [[ -e "${FINAL_ROOT}" || -L "${FINAL_ROOT}" ]]; then
    if [[ -f "${FINAL_ROOT}/${CAPTURE_FILENAME}" ]]; then
        echo "Refusing to overwrite completed captured reference: ${FINAL_ROOT}" >&2
        exit 2
    fi
    echo "Removing incomplete prior capture: ${FINAL_ROOT}"
    rm -rf -- "${FINAL_ROOT}"
fi
rm -rf -- "${STAGING_ROOT}"
mkdir -p "${STAGING_ROOT}"

promoted=0
cleanup() {
    status=$?
    if (( status != 0 && promoted == 0 )); then
        rm -rf -- "${STAGING_ROOT}"
    fi
    exit "${status}"
}
trap cleanup EXIT

PYTHONPATH="${REPO_ROOT}/src${PYTHONPATH:+:${PYTHONPATH}}" \
python3 -m oraclut.validation.test_case \
    --repository-root "${REPO_ROOT}" \
    --lut-file "${TEST_LUT}" \
    --channels 1 \
    --output-manifest "${STAGING_ROOT}/test_case_manifest.json"

PYTHONPATH="${REPO_ROOT}/src${PYTHONPATH:+:${PYTHONPATH}}" \
python3 -m oraclut.validation.capture fingerprint "${SOURCE_NC}" \
    --output "${STAGING_ROOT}/source_before.json"

(
    source "${REPO_ROOT}/scripts/legacy_idl_env.sh"
    export IDL_PATH="${IDL_PATH}:${REPO_ROOT}/scripts"
    export CREATE_ORAC_CLOUD_LUT_DRIVER="${DRIVER}"
    cd "${LEGACY_ROOT}"
    [[ "${PWD}" == "${LEGACY_ROOT}" ]] || { echo "Could not establish legacy cwd" >&2; exit 2; }

    echo "Preflight: resolving required current cloud-path routines"
    idl -e "resolve_routine, 'create_orac_cloud_lut_wrapper'; resolve_routine, 'create_orac_cloud_lut'; resolve_routine, 'generate_scattering_properties'; resolve_routine, 'create_bwgp'; resolve_routine, 'load_inststr'; resolve_routine, 'load_lutstr'; resolve_routine, 'load_srfstrarr'; resolve_routine, 'load_mmdat'; resolve_routine, 'load_atmstr'; resolve_routine, 'setup_disort'; resolve_routine, 'call_disort'; resolve_routine, 'write_v2_lut'; print, 'Legacy routine preflight passed'; exit"

    start_time="$(date +%s)"
    echo "Running current cloud wrapper from ${PWD}"
    idl -e "create_orac_cloud_lut_wrapper,channelid=[1],srf_quad=1,mmfile='liquid-water_stg.mm',lutfile='liquid-water-cloud_test.lut',version=21"
    end_time="$(date +%s)"
    echo "Elapsed seconds: $((end_time - start_time))"
)

PYTHONPATH="${REPO_ROOT}/src${PYTHONPATH:+:${PYTHONPATH}}" \
python3 -m oraclut.validation.capture capture "${SOURCE_NC}" \
    "${STAGING_ROOT}/${CAPTURE_FILENAME}" \
    --before "${STAGING_ROOT}/source_before.json" \
    --manifest "${STAGING_ROOT}/capture_manifest.json" \
    --configuration "${STAGING_ROOT}/test_case_manifest.json" \
    --idl-version "IDL 8.9.0" \
    --mie-dlm-path "${MIE_DLM}" \
    --disort-dlm-path "${DISORT_DLM}"

mv -- "${STAGING_ROOT}" "${FINAL_ROOT}"
promoted=1
trap - EXIT
echo "Captured verified legacy reference in ${FINAL_ROOT}"
find "${FINAL_ROOT}" -maxdepth 2 -type f -printf '  %p\n' | sort
