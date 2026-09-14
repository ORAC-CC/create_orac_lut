#!/usr/bin/env bash
# Generate the legacy IDL aerosol reference for the established tiny channel-1
# case (Meteosat-10 SEVIRI, aerosol_a79, aerosol_test.lut, channel 1,
# atmosphere 2, srf_quad 1, gas and Rayleigh on) and compare it with the
# Python product.  This is a thin front end to the general provenance-gated
# runner scripts/generate_legacy_reference.sh, which holds all the logic and
# the single definition of the legacy source root.
#
#     cd /home/g/grainger/project-oraclut
#     scripts/generate_aerosol_legacy_reference.sh [--repeat|--no-repeat] [--no-compare]

set -euo pipefail
REPO_ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd -P)"
exec "${REPO_ROOT}/scripts/generate_legacy_reference.sh" \
    --formulation aerosol \
    --driver validation/reference_cases/meteosat10_seviri_aerosol_a79_test.driver \
    --case-id aerosol_visible \
    --python-product create_orac_lut/luts/meteosat-10_seviri_aerosol/meteosat-10_seviri_m_aerosol_a12_pa79_v21.nc \
    --reference-dir validation/generated/aerosol_legacy_reference \
    --comparison-dir validation/generated/aerosol_comparison \
    "$@"
