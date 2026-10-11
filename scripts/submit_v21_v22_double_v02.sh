#!/bin/bash
#
# submit_v21_v22_double.sh
#
# Submit the double-precision V21 (IDL) and V22 (Python) Terra MODIS LUTs:
#
#   1. v21dp_validate     IDL integration validation (four compact cases,
#                         double and single precision; ~15 min)
#   2. v21dp_pold/p240/pagg   complete V21 products, each starting only if
#                         job 1 succeeded (afterok; cancelled otherwise)
#   3. v22dp_pold/p240/pagg   complete V22 products (release manifest
#                         releases/lut_v22_double.release), independent
#
# Products go only to the protected staging directories
#   staging/lut_v21_double/   staging/lut_v22_double/
# never to the shared ORAC LUT directory.  A product already staged (and, for
# V21, verified) is not submitted again, so the command is safe to repeat.
#
# Run from any of atmlxint1/2/3/5/6 (it re-runs itself on atmlxint7 over ssh)
# or on atmlxint7.  Partition priority-eodg (atmnode009-010), one CPU each.
# Job IDs are recorded in staging/submissions.tsv; follow progress with
#    bash scripts/monitor_v21_v22_double.sh
#
# Usage: bash scripts/submit_v21_v22_double.sh [--v21-only]
#   --v21-only   submit only the IDL validation and V21 product jobs (e.g. after
#                a V21 failure while the V22 jobs are still running)

set -euo pipefail
export PATH=/usr/bin:/bin:$PATH

REPO="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd -P)"
cd "$REPO"
SUBMIT_HOST=atmlxint7
DISPATCH_HOSTS=" atmlxint1 atmlxint2 atmlxint3 atmlxint5 atmlxint6 "
PYTHON=/home/g/grainger/miniforge3/envs/science/bin/python
STAGE21=staging/lut_v21_double
STAGE22=staging/lut_v22_double
RECORD=staging/submissions.tsv
PRODUCTS=(pold p240 pagg)
declare -A V21FILE=([pold]=terra_modis_m_liquid-water_a01_pold_v21.nc [p240]=terra_modis_m_liquid-water_a01_p240_v21.nc [pagg]=terra_modis_m_water-ice_a01_pagg_v21.nc)
declare -A V22FILE=([pold]=terra_modis_m_liquid-water_a01_pold_v22.nc [p240]=terra_modis_m_liquid-water_a01_p240_v22.nc [pagg]=terra_modis_m_water-ice_a01_pagg_v22.nc)
declare -A V22RUN=([pold]=runs/terra_modis_cloud_liquid-water_old_v22_double.run [p240]=runs/terra_modis_cloud_liquid-water_240_v22_double.run [pagg]=runs/terra_modis_cloud_water-ice_agg_v22_double.run)

die() { echo "submit_v21_v22_double.sh: $*" >&2; exit 2; }
V21_ONLY=0
case "${1:-}" in "") ;; --v21-only) V21_ONLY=1 ;; *) die "unknown option $1" ;; esac

host="$(hostname -s 2>/dev/null || hostname)"; host="${host%%.*}"
if [[ "$host" != "$SUBMIT_HOST" ]]; then
   [[ "${ORACLUT_DISPATCHED:-0}" == 1 ]] && die "dispatched invocation landed on $host"
   [[ "$DISPATCH_HOSTS" == *" $host "* ]] || die "run this from atmlxint1/2/3/5/6 or $SUBMIT_HOST, not $host"
   echo "Dispatching to $SUBMIT_HOST over ssh (repository is on the shared filesystem)"
   ssh -o BatchMode=yes -o ConnectTimeout=20 "$SUBMIT_HOST" \
      bash -lc "$(printf '%q' "cd $(printf '%q' "$REPO") && ORACLUT_DISPATCHED=1 exec bash scripts/submit_v21_v22_double.sh ${1:-}")" \
      || die "ssh to $SUBMIT_HOST failed (an expired Kerberos ticket? run kinit and repeat)"
   exit 0
fi

# ---------------------------------------------------------------------------
# Pre-flight on atmlxint7: every release file verified and committed.
# ---------------------------------------------------------------------------
"$PYTHON" scripts/check_release_manifest.py releases/idl_lut_v21_double.release releases/lut_v22_double.release \
   || die "release files differ from their manifests"
ORAC_LUT_RELEASE_MANIFEST=lut_v22_double.release "$PYTHON" create_orac_luts_release_v03.py --check 22 > /dev/null \
   || die "the V22 double-precision release does not verify"
uncommitted=$(git status --porcelain --untracked-files=all -- releases create_orac_lut/call_disort.pro create_orac_lut/disort2_dp \
              create_orac_luts.py create_orac_luts_v03.py create_orac_luts_release_v03.py src/oraclut CODE_VERSION \
              scripts/run_idl_v21_lut.sh scripts/run_idl_v21_lut.slurm scripts/run_oraclut_lut.slurm runs/*_v22_double.run \
              | grep -v -E '__pycache__' || true)
[[ -z "$uncommitted" ]] || die "release files are not committed (the products must record a committed source state):
$uncommitted"

mkdir -p validation/slurm "$STAGE21" "$STAGE22"
[[ -f "$RECORD" ]] || printf 'time\thost\tcommit\tjob\tslurm_id\tdetail\n' > "$RECORD"
commit=$(git rev-parse HEAD)
unset TMPDIR TMP TEMP
note() { printf '%s\t%s\t%s\t%s\t%s\t%s\n' "$(date --iso-8601=seconds)" "$host" "$commit" "$1" "$2" "$3" >> "$RECORD"; echo "  $1: $2 $3"; }

echo "Submitting double-precision V21/V22 Terra MODIS LUT jobs from $host at $commit"
stamp=$(date +%Y%m%dT%H%M%S)
validate_id=$(sbatch --parsable --job-name=v21dp_validate --time=01:30:00 \
              scripts/run_idl_v21_lut.slurm validate "$STAGE21/integration_validation_$stamp")
note v21dp_validate "$validate_id" "$STAGE21/integration_validation_$stamp"

for p in "${PRODUCTS[@]}"; do
   if [[ -f "$STAGE21/${V21FILE[$p]}.verification.json" ]] && grep -q '"passed": true' "$STAGE21/${V21FILE[$p]}.verification.json"; then
      note "v21dp_$p" skipped "already staged and verified"
   else
      id=$(sbatch --parsable --job-name="v21dp_$p" --dependency="afterok:$validate_id" --kill-on-invalid-dep=yes \
           scripts/run_idl_v21_lut.slurm product "$p" "$STAGE21")
      note "v21dp_$p" "$id" "after $validate_id"
   fi
done

for p in "${PRODUCTS[@]}"; do
   if [[ "$V21_ONLY" == 1 ]]; then
      note "v22dp_$p" skipped "--v21-only"
   elif [[ -f "$STAGE22/${V22FILE[$p]}" ]]; then
      note "v22dp_$p" skipped "already staged"
   else
      id=$(ORAC_LUT_RELEASE_MANIFEST=lut_v22_double.release ORAC_LUT_DISORT_PRECISION=double \
           sbatch --parsable --export=ALL --job-name="v22dp_$p" --partition=priority-eodg --time=14:00:00 --mem=8G \
           scripts/run_oraclut_lut.slurm "${V22RUN[$p]}")
      note "v22dp_$p" "$id" "${V22RUN[$p]}"
   fi
done
echo "Submitted.  Follow with: bash scripts/monitor_v21_v22_double.sh"
