#!/usr/bin/env bash
# Submit the twenty Sentinel-3A/B SLSTR V23 cloud reference LUTs.
#
# Usage, from an interactive AOPP host:
#   scripts/submit_v23_slstr_cloud_luts.sh             submit every incomplete product
#   scripts/submit_v23_slstr_cloud_luts.sh --dry-run   report what would be submitted
#
# The user's SBATCH_PARTITION environment selects the eligible partitions.
# This script deliberately does not pass --partition.

set -euo pipefail

REPO_ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd -P)"
cd "$REPO_ROOT"

RUN_DIR="runs"
SUBMIT="scripts/submit_oraclut_lut.sh"
LUT_DIR="/network/group/aopp/eodg/RGG004_GRAINGER_ORACFILE/ORAC_LUTS"
VERSION="23"
SLURM_OPTIONS=(--time=72:00:00 --mem=24G)

platforms=(sentinel-3a sentinel-3b)
liquid_water=(liquid-water_old liquid-water_stg liquid-water_240 liquid-water_253 liquid-water_263 liquid-water_273)
ice=(water-ice_sph water-ice_agg water-ice_ghm water-ice_src)

die() { echo "submit_v23_slstr_cloud_luts.sh: $*" >&2; exit 2; }

dry_run=0
case "${1:-}" in
    "") ;;
    --dry-run) dry_run=1 ;;
    -h|--help) sed -n '2,11p' "${BASH_SOURCE[0]}" | sed 's/^# \\{0,1\\}//'; exit 0 ;;
    *) die "unknown argument: $1 (use --dry-run or no argument)" ;;
esac
[[ "$#" -le 1 ]] || die "at most one argument (--dry-run) is accepted"
[[ -n "${SBATCH_PARTITION:-}" ]] || die "SBATCH_PARTITION is not set"
[[ -d "$LUT_DIR" ]] || die "permanent LUT directory not found: $LUT_DIR"
[[ -w "$LUT_DIR" ]] || die "permanent LUT directory is not writable: $LUT_DIR"

run_value() {
    local name="$1" run="$2"
    sed -n "s/^[[:space:]]*${name}[[:space:]]*=[[:space:]]*//p" "$run" |
        sed "s/[[:space:]]*;.*$//; s/^'//; s/'$//"
}

check_run() {
    local run="$1" platform="$2" model="$3" lutfile="$4"
    [[ -f "$run" ]] || die "run file not found: $run"
    [[ "$(run_value platform "$run")" == "$platform" ]] || die "$run: unexpected platform"
    [[ "$(run_value instrument "$run")" == "slstr" ]] || die "$run: instrument is not slstr"
    [[ "$(run_value forward_model "$run")" == "cloud" ]] || die "$run: forward_model is not cloud"
    [[ "$(run_value instfile "$run")" == "${platform}_slstr_v1.inst" ]] || die "$run: unexpected instfile"
    [[ "$(run_value mmfile "$run")" == "${model}.mm" ]] || die "$run: unexpected mmfile"
    [[ "$(run_value lutfile "$run")" == "$lutfile" ]] || die "$run: unexpected Grid B file"
    [[ "$(run_value atmospheres "$run")" == "2" ]] || die "$run: atmosphere is not 2"
    [[ "$(run_value srf_quad "$run")" == "1" ]] || die "$run: srf_quad is not 1"
    [[ "$(run_value gas "$run")" == "0" ]] || die "$run: gas is not 0"
    [[ "$(run_value no_rayleigh "$run")" == "0" ]] || die "$run: no_rayleigh is not 0"
    [[ "$(run_value nstreams "$run")" == "60" ]] || die "$run: nstreams is not 60"
    [[ "$(run_value nmom "$run")" == "1000" ]] || die "$run: nmom is not 1000"
    [[ "$(run_value version "$run")" == "$VERSION" ]] || die "$run: version is not $VERSION"
    [[ "$(run_value out_path "$run")" == "$LUT_DIR" ]] || die "$run: unexpected out_path"
}

submitted=0
complete=0
for platform in "${platforms[@]}"; do
    for model in "${liquid_water[@]}" "${ice[@]}"; do
        case "$model" in
            liquid-water_*) lutfile="liquid-water-cloud-grid-b.lut"; family="water" ;;
            water-ice_*)    lutfile="ice-cloud-grid-b.lut";          family="ice" ;;
        esac
        substance="${model%%_*}"
        shortname="${model#*_}"
        run="${RUN_DIR}/${platform}_slstr_cloud_${model}_v${VERSION}.run"
        product="${LUT_DIR}/${platform}_slstr_m_${substance}_a01_p${shortname}_v${VERSION}.nc"
        check_run "$run" "$platform" "$model" "$lutfile"

        if [[ -s "$product" ]]; then
            echo "COMPLETE  ${platform} ${model}: ${product}"
            complete=$((complete + 1))
            continue
        fi

        job_name="v${VERSION}_${platform}_${family}_${shortname}"
        if [[ "$dry_run" -eq 1 ]]; then
            echo "WOULD SUBMIT  ${platform} ${model}: ${run} -> ${product}"
            echo "    $SUBMIT ${SLURM_OPTIONS[*]} --job-name=${job_name} ${run}"
        else
            echo "SUBMIT  ${platform} ${model}: ${run} -> ${product}"
            "$SUBMIT" "${SLURM_OPTIONS[@]}" --job-name="$job_name" "$run"
        fi
        submitted=$((submitted + 1))
    done
done

if [[ "$dry_run" -eq 1 ]]; then
    echo "Dry run: ${complete} complete, ${submitted} would be submitted; nothing was submitted."
else
    echo "${complete} complete, ${submitted} submitted."
fi

