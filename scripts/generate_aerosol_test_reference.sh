#!/usr/bin/env bash
set -euo pipefail

# Run the current legacy aerosol/distributed-profile formulation on the compact
# aerosol_test.lut grid. The IDL process must run from create_orac_lut because
# the historical source uses input_files/... relative paths.
REPO_ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd -P)"
LEGACY_ROOT="${REPO_ROOT}/create_orac_lut"
CASE_ID="meteosat10_seviri_aerosol_a79_test"
DRIVER="${REPO_ROOT}/validation/reference_cases/${CASE_ID}.driver"
LEGACY_OUTPUT="${LEGACY_ROOT}/luts/${CASE_ID}"
CAPTURE_DIR="${REPO_ROOT}/validation/generated/${CASE_ID}"

case "${REPO_ROOT}" in
    /home/g/grainger/project-oraclut) ;;
    *) echo "Refusing unexpected repository root: ${REPO_ROOT}" >&2; exit 2 ;;
esac

require_readable() {
    if [[ ! -r "$1" ]]; then
        echo "Preflight failed: missing or unreadable $2: $1" >&2
        exit 2
    fi
}

require_readable "${DRIVER}" "compact aerosol driver"
require_readable "${LEGACY_ROOT}/input_files/inst/meteosat-10_seviri_v1.inst" "SEVIRI instrument"
require_readable "${LEGACY_ROOT}/input_files/microphysics/aerosol_a79.mm" "a79 microphysics"
require_readable "${LEGACY_ROOT}/input_files/lut/aerosol_test.lut" "aerosol test grid"
require_readable "${REPO_ROOT}/mie/dlm-code/mie_dlm_single.dlm" "Mie DLM declaration"
require_readable "${REPO_ROOT}/mie/dlm-code/mie_dlm_single.so" "Mie DLM shared object"
require_readable "${REPO_ROOT}/create_orac_lut/disort2/src/DISORTtoIDL.dlm" "production DISORT DLM declaration"
require_readable "${REPO_ROOT}/create_orac_lut/disort2/src/DISORTtoIDL.so" "production DISORT DLM shared object"

if [[ -e "${CAPTURE_DIR}" || -e "${LEGACY_OUTPUT}" ]]; then
    echo "Refusing to overwrite an existing aerosol reference/output directory" >&2
    echo "  ${CAPTURE_DIR}" >&2
    echo "  ${LEGACY_OUTPUT}" >&2
    exit 2
fi

echo "Compact legacy aerosol reference preflight"
echo "  driver: ${DRIVER}"
echo "  model: aerosol_a79.mm"
echo "  LUT: aerosol_test.lut"
echo "  channel: 1; atmosphere: 2; SRF mode: 1; Rayleigh: enabled; gas: enabled"
echo "  DLMs: ${REPO_ROOT}/mie/dlm-code and ${REPO_ROOT}/create_orac_lut/disort2/src"

(
    source "${REPO_ROOT}/scripts/legacy_idl_env.sh"
    export IDL_PATH="${IDL_PATH}:${REPO_ROOT}/scripts"
    export CREATE_ORAC_AEROSOL_LUT_DRIVER="${DRIVER}"
    cd "${LEGACY_ROOT}"
    [[ "${PWD}" == "${LEGACY_ROOT}" ]] || { echo "Could not establish legacy cwd" >&2; exit 2; }

    echo "Preflight: resolving current aerosol-path routines"
    idl -e "resolve_routine, 'create_orac_aerosol_lut_wrapper'; resolve_routine, 'create_orac_aerosol_lut'; resolve_routine, 'generate_scattering_properties'; resolve_routine, 'create_bwgp'; resolve_routine, 'load_inststr'; resolve_routine, 'load_lutstr'; resolve_routine, 'load_srfstrarr'; resolve_routine, 'load_mmdat'; resolve_routine, 'load_atmstr'; resolve_routine, 'setup_disort'; resolve_routine, 'call_disort'; resolve_routine, 'write_v2_lut'; print, 'Legacy aerosol routine preflight passed'; exit"
    idl -e "create_orac_aerosol_lut_wrapper"
)

count="$(find "${LEGACY_OUTPUT}" -maxdepth 1 -type f -name '*.nc' | wc -l)"
if [[ "${count}" != "1" ]]; then
    echo "Expected exactly one generated aerosol NetCDF file, found ${count}" >&2
    exit 2
fi
mkdir -p "${CAPTURE_DIR}"
source_file="$(find "${LEGACY_OUTPUT}" -maxdepth 1 -type f -name '*.nc' -print)"
cp -- "${source_file}" "${CAPTURE_DIR}/${CASE_ID}_legacy_reference_v21.nc"
printf '%s\n' "${DRIVER}" > "${CAPTURE_DIR}/source_driver.txt"
echo "Captured compact legacy aerosol reference in ${CAPTURE_DIR}"
