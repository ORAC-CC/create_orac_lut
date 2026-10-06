#!/usr/bin/env bash
# Submit the forty V25 cloud LUT production calculations.
#
# Usage, from an interactive AOPP host (atmlxint5 is the usual one):
#   scripts/submit_v25_cloud_luts.sh [modis|slstr|all]             submit every incomplete product
#   scripts/submit_v25_cloud_luts.sh --dry-run [modis|slstr|all]   report what would be submitted; submit nothing
#
# V25 is the V24 product set (same platforms, models, channels, Grid B
# sampling and numerics) computed with the adiabatic cloud temperature profile
# recorded in create_orac_luts.py ("Versions and numerical changes" 5):
# every run file sets cloud_temperature_profile = 'adiabatic' and version = 25
# (runs/V25_PRODUCTION.md).  Products: Aqua and Terra MODIS and
# Sentinel-3A/B SLSTR (dual view), each with six liquid-water models (liquid
# water Grid B) and four ice models (ice Grid B).  The run files are
#   runs/<platform>_<instrument>_cloud_<model>_v25.run
# and each writes its NetCDF product directly to LUT_DIR as
#   <platform>_<instrument>_m_<substance>_a01_p<shortname>_v25.nc
#
# Provenance: nothing is submitted unless the production source, inputs, run
# files and batch scripts are unmodified in the working tree and HEAD is
# contained in origin/main (fetched first).  The jobs run this shared working
# tree, so it must not change until they have started; each job log records
# the revision it ran (run_oraclut_lut.slurm).  Every submission is appended
# to MANIFEST (time, host, revision, run file, product, SLURM job ID).
#
# Restart: a product that already exists and is non-empty in LUT_DIR is
# reported as complete and not submitted again; the generator also refuses to
# overwrite an existing product.  V23 and V24 products are never touched.
#
# Each job is submitted through scripts/submit_oraclut_lut.sh, which forwards
# to atmlxint7; SBATCH_PARTITION selects the eligible partitions.

set -euo pipefail

REPO_ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd -P)"
cd "$REPO_ROOT"

RUN_DIR="runs"
SUBMIT="scripts/submit_oraclut_lut.sh"
LUT_DIR="/network/group/aopp/eodg/RGG004_GRAINGER_ORACFILE/ORAC_LUTS"
MANIFEST="validation/v25/production_submissions.tsv"
VERSION="25"
SLURM_OPTIONS=(--time=72:00:00 --mem=24G)
# Paths whose content defines the products.
PRODUCTION_PATHS=(create_orac_luts.py src mie create_orac_lut runs scripts/run_oraclut_lut.slurm
                  scripts/submit_oraclut_lut.sh scripts/submit_v25_cloud_luts.sh)
# ssh to atmlxint7 and git fetch need the caller's ssh agent; a Kerberos-free
# PATH keeps the wrapper usable from a shell whose ticket has expired.

liquid_water=(liquid-water_old liquid-water_stg liquid-water_240 liquid-water_253 liquid-water_263 liquid-water_273)
ice=(water-ice_sph water-ice_agg water-ice_ghm water-ice_src)
baum=(water-ice_agg water-ice_ghm water-ice_src)

die() { echo "submit_v25_cloud_luts.sh: $*" >&2; exit 2; }

dry_run=0
selection="all"
for argument in "$@"; do
    case "$argument" in
        --dry-run) dry_run=1 ;;
        modis|slstr|all) selection="$argument" ;;
        -h|--help) sed -n '2,30p' "${BASH_SOURCE[0]}" | sed 's/^# \{0,1\}//'; exit 0 ;;
        *) die "unknown argument: $argument" ;;
    esac
done
case "$selection" in
    modis) instruments=(modis) ;;
    slstr) instruments=(slstr) ;;
    all)   instruments=(modis slstr) ;;
esac

[[ -n "${SBATCH_PARTITION:-}" ]] || die "SBATCH_PARTITION is not set"
[[ -d "$LUT_DIR" ]] || die "permanent LUT directory not found: $LUT_DIR"
[[ -w "$LUT_DIR" ]] || die "permanent LUT directory is not writable: $LUT_DIR"

# ---------------------------------------------------------------------------
# Provenance: a clean, pushed revision.
# ---------------------------------------------------------------------------
# A dry run reports a failed check and carries on; a submission stops.
provenance() { if [[ "$dry_run" -eq 1 ]]; then echo "DRY RUN WARNING: $*" >&2; else die "$*"; fi; }
revision="$(git rev-parse HEAD)"
dirty="$(git status --porcelain --untracked-files=all -- "${PRODUCTION_PATHS[@]}" | grep -v -E '^\?\? runs/.*_v23\.run$|^\?\? scripts/submit_v23_slstr_cloud_luts\.sh$' || true)"
if [[ -n "$dirty" ]]; then
    echo "$dirty" >&2
    provenance "production paths differ from HEAD $revision; commit (and push) them first"
elif ! git fetch --quiet origin main; then
    provenance "git fetch origin main failed; cannot confirm that $revision is on GitHub"
elif ! git merge-base --is-ancestor "$revision" origin/main; then
    provenance "HEAD $revision is not contained in origin/main; push it first"
else
    echo "Revision: $revision (clean production paths; contained in origin/main)"
fi

# Value of one "name = 'value'" setting in a run file, quotes removed.
run_value() {
    local name="$1" run="$2"
    sed -n "s/^[[:space:]]*${name}[[:space:]]*=[[:space:]]*//p" "$run" | sed "s/[[:space:]]*;.*$//; s/^'//; s/'$//"
}

is_baum() {
    local candidate
    for candidate in "${baum[@]}"; do [[ "$1" == "$candidate" ]] && return 0; done
    return 1
}

# Stop before submitting anything if a run file does not say what this script
# expects: the expected filename below is only valid for these settings.
check_run() {
    local run="$1" platform="$2" instrument="$3" model="$4" lutfile="$5" channels="$6"
    [[ -f "$run" ]] || die "run file not found: $run"
    [[ "$(run_value platform "$run")" == "$platform" ]] || die "$run: platform is not '$platform'"
    [[ "$(run_value instrument "$run")" == "$instrument" ]] || die "$run: instrument is not '$instrument'"
    [[ "$(run_value forward_model "$run")" == "cloud" ]] || die "$run: forward_model is not 'cloud'"
    [[ "$(run_value instfile "$run")" == "${platform}_${instrument}_v1.inst" ]] || die "$run: unexpected instfile"
    [[ "$(run_value mmfile "$run")" == "${model}.mm" ]] || die "$run: mmfile is not ${model}.mm"
    [[ "$(run_value lutfile "$run")" == "$lutfile" ]] || die "$run: lutfile is not $lutfile"
    [[ "$(run_value atmospheres "$run")" == "2" ]] || die "$run: atmospheres is not 2"
    [[ "$(run_value channelid "$run")" == "$channels" ]] || die "$run: unexpected channelid"
    [[ "$(run_value srf_quad "$run")" == "1" ]] || die "$run: srf_quad is not 1 (filename 'm')"
    [[ "$(run_value gas "$run")" == "0" ]] || die "$run: gas is not 0 (filename 'a01')"
    [[ "$(run_value no_rayleigh "$run")" == "0" ]] || die "$run: no_rayleigh is not 0 (filename 'a01')"
    [[ "$(run_value nstreams "$run")" == "60" ]] || die "$run: nstreams is not 60"
    if is_baum "$model"; then
        [[ "$(run_value nmom "$run")" == "1000" ]] || die "$run: a Baum class needs nmom = 1000"
    else
        [[ -z "$(run_value nmom "$run")" ]] || die "$run: nmom must not be set for a Mie class"
    fi
    [[ "$(run_value version "$run")" == "$VERSION" ]] || die "$run: version is not $VERSION"
    [[ "$(run_value cloud_temperature_profile "$run")" == "adiabatic" ]] || die "$run: cloud_temperature_profile is not 'adiabatic' (V25)"
    [[ "$(run_value out_path "$run")" == "$LUT_DIR" ]] || die "$run: out_path is not $LUT_DIR"
}

modis_channels="[$(seq -s ', ' 1 36)]"
slstr_channels="[$(seq -s ', ' 1 9)]"

mkdir -p "$(dirname "$MANIFEST")"
if [[ "$dry_run" -eq 0 && ! -s "$MANIFEST" ]]; then
    printf 'submitted_utc\thost\trevision\tjob_id\tjob_name\trun_file\tproduct\n' > "$MANIFEST"
fi

# Check every run file before anything is submitted.
declare -a queue=()
complete=0
for instrument in "${instruments[@]}"; do
    case "$instrument" in
        modis) platforms=(aqua terra); channels="$modis_channels" ;;
        slstr) platforms=(sentinel-3a sentinel-3b); channels="$slstr_channels" ;;
    esac
    for platform in "${platforms[@]}"; do
        for model in "${liquid_water[@]}" "${ice[@]}"; do
            case "$model" in
                liquid-water_*) lutfile="liquid-water-cloud-grid-b.lut"; family="water" ;;
                water-ice_*)    lutfile="ice-cloud-grid-b.lut";          family="ice" ;;
            esac
            substance="${model%%_*}"
            shortname="${model#*_}"
            run="${RUN_DIR}/${platform}_${instrument}_cloud_${model}_v${VERSION}.run"
            product="${LUT_DIR}/${platform}_${instrument}_m_${substance}_a01_p${shortname}_v${VERSION}.nc"
            check_run "$run" "$platform" "$instrument" "$model" "$lutfile" "$channels"
            if [[ -s "$product" ]]; then
                echo "COMPLETE  ${platform} ${instrument} ${model}: ${product}"
                complete=$((complete + 1))
                continue
            fi
            queue+=("${run}|${product}|v${VERSION}_${platform}_${family}_${shortname}")
        done
    done
done

submitted=0
for entry in "${queue[@]+"${queue[@]}"}"; do
    IFS='|' read -r run product job_name <<< "$entry"
    if [[ "$dry_run" -eq 1 ]]; then
        echo "WOULD SUBMIT  ${run} -> ${product}"
        echo "    $SUBMIT ${SLURM_OPTIONS[*]} --job-name=${job_name} ${run}"
    else
        echo "SUBMIT  ${run} -> ${product}"
        output="$("$SUBMIT" "${SLURM_OPTIONS[@]}" --job-name="$job_name" "$run" 2>&1)" || { echo "$output" >&2; die "submission failed for $run"; }
        echo "$output"
        job_id="$(sed -n 's/^Submitted batch job \([0-9][0-9]*\).*$/\1/p' <<< "$output" | tail -1)"
        [[ -n "$job_id" ]] || die "no SLURM job ID in the submission output for $run"
        printf '%s\t%s\t%s\t%s\t%s\t%s\t%s\n' "$(date -u +%Y-%m-%dT%H:%M:%SZ)" "$(hostname -s)" "$revision" \
            "$job_id" "$job_name" "$run" "$product" >> "$MANIFEST"
    fi
    submitted=$((submitted + 1))
done

if [[ "$dry_run" -eq 1 ]]; then
    echo "Dry run: ${complete} complete, ${submitted} would be submitted; nothing was submitted."
else
    echo "${complete} complete, ${submitted} submitted; recorded in ${MANIFEST}."
fi
