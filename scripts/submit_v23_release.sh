#!/bin/bash
#
# submit_v23_release.sh
#
# Submit the three V23 Terra MODIS LUTs (lut-code-v22.4 code, Grid B)
# (releases/lut_v23.release; launcher create_orac_luts_release_v04.py):
#
#   v23_pold   runs/terra_modis_cloud_liquid-water_old_v23.run
#   v23_p240   runs/terra_modis_cloud_liquid-water_240_v23.run
#   v23_pagg   runs/terra_modis_cloud_water-ice_agg_v23.run
#
# Products go only to the protected staging directory staging/lut_v23/,
# never to the shared ORAC LUT directory.  A product already staged is not
# submitted again, so the command is safe to repeat.
#
# Run from any of atmlxint1/2/3/5/6 (it re-runs itself on atmlxint7 over ssh)
# or on atmlxint7.  Partition priority-eodg (atmnode009-010), one CPU, 24 GB,
# 48 h each.  Job IDs are recorded in staging/submissions.tsv; follow with
#    bash scripts/monitor_v23_release.sh [--compare]
#
# Usage: bash scripts/submit_v23_release.sh
set -euo pipefail
export PATH=/usr/bin:/bin:$PATH
REPO="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd -P)"
cd "$REPO"
SUBMIT_HOST=atmlxint7
DISPATCH_HOSTS=" atmlxint1 atmlxint2 atmlxint3 atmlxint5 atmlxint6 "
PYTHON=/home/g/grainger/miniforge3/envs/science/bin/python
STAGE=staging/lut_v23
RECORD=staging/submissions.tsv
PRODUCTS=(pold p240 pagg)
declare -A FILE=([pold]=terra_modis_m_liquid-water_a01_pold_v23.nc [p240]=terra_modis_m_liquid-water_a01_p240_v23.nc [pagg]=terra_modis_m_water-ice_a01_pagg_v23.nc)
declare -A RUN=([pold]=runs/terra_modis_cloud_liquid-water_old_v23.run [p240]=runs/terra_modis_cloud_liquid-water_240_v23.run [pagg]=runs/terra_modis_cloud_water-ice_agg_v23.run)
die() { echo "submit_v23_release.sh: $*" >&2; exit 2; }
[[ $# -eq 0 ]] || die "takes no arguments"

host="$(hostname -s 2>/dev/null || hostname)"; host="${host%%.*}"
if [[ "$host" != "$SUBMIT_HOST" ]]; then
   [[ "${ORACLUT_DISPATCHED:-0}" == 1 ]] && die "dispatched invocation landed on $host"
   [[ "$DISPATCH_HOSTS" == *" $host "* ]] || die "run this from atmlxint1/2/3/5/6 or $SUBMIT_HOST, not $host"
   echo "Dispatching to $SUBMIT_HOST over ssh (repository is on the shared filesystem)"
   ssh -o BatchMode=yes -o ConnectTimeout=20 "$SUBMIT_HOST" \
      bash -lc "$(printf '%q' "cd $(printf '%q' "$REPO") && ORACLUT_DISPATCHED=1 exec bash scripts/submit_v23_release.sh")" \
      || die "ssh to $SUBMIT_HOST failed (an expired Kerberos ticket? run kinit and repeat)"
   exit 0
fi

# Pre-flight on atmlxint7: every release file verified and committed.
ORAC_LUT_RELEASE_MANIFEST=lut_v23.release "$PYTHON" create_orac_luts_release_v04.py --check 23 > /dev/null \
   || die "the V23 release does not verify"
"$PYTHON" scripts/check_release_manifest.py releases/idl_lut_v21_single.release > /dev/null \
   || die "the V21 reference release files differ from their manifest"
uncommitted=$(git status --porcelain --untracked-files=all -- releases CODE_VERSION create_orac_luts.py create_orac_luts_release_v04.py \
              src/oraclut create_orac_lut/disort2/src scripts/run_v23_release.slurm runs/*_v23.run | grep -v -E '__pycache__' || true)
[[ -z "$uncommitted" ]] || die "release files are not committed (the products must record a committed source state):
$uncommitted"
mkdir -p validation/slurm "$STAGE"
[[ -f "$RECORD" ]] || printf 'time\thost\tcommit\tjob\tslurm_id\tdetail\n' > "$RECORD"
commit=$(git rev-parse HEAD)
unset TMPDIR TMP TEMP
note() { printf '%s\t%s\t%s\t%s\t%s\t%s\n' "$(date --iso-8601=seconds)" "$host" "$commit" "$1" "$2" "$3" >> "$RECORD"; echo "  $1: $2 $3"; }
echo "Submitting V23 Terra MODIS LUT jobs from $host at $commit"
for p in "${PRODUCTS[@]}"; do
   if [[ -f "$STAGE/${FILE[$p]}" ]]; then
      note "v23_$p" skipped "already staged"
   else
      id=$(sbatch --parsable --job-name="v23_$p" scripts/run_v23_release.slurm "${RUN[$p]}")
      note "v23_$p" "$id" "${RUN[$p]}"
   fi
done
echo "Submitted.  Follow with: bash scripts/monitor_v23_release.sh"
