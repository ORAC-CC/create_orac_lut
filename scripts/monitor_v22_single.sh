#!/bin/bash
#
# monitor_v22_single.sh
#
# Report the corrected single-precision V22 LUT jobs (v22sp_*) recorded in
# staging/submissions.tsv and the products staged in staging/lut_v22_single-r2/.
# With --compare, compare each staged product with the production V21 IDL
# product (validation/v21_v22_terra_modis/v21_v22_single/compare_full_luts.py).
# Read-only apart from the comparison results.
#
# Usage: bash scripts/monitor_v22_single.sh [--compare]

set -u
export PATH=/usr/bin:/bin:$PATH
REPO="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd -P)"
cd "$REPO"
PYTHON=/home/g/grainger/miniforge3/envs/science/bin/python

if [[ -f staging/submissions.tsv ]]; then
   echo "== Submissions (v22sp_*)"
   awk -F'\t' 'NR == 1 || $4 ~ /^v22sp_/' staging/submissions.tsv | column -t -s $'\t'
   ids=$(awk -F'\t' '$4 ~ /^v22sp_/ && $5 ~ /^[0-9]+$/ {print $5}' staging/submissions.tsv | paste -sd, -)
   if [[ -n "$ids" ]] && command -v sacct > /dev/null; then
      echo; echo "== SLURM state"
      sacct -j "$ids" --format=JobID,JobName%16,State,Elapsed,NodeList,ExitCode -X 2>/dev/null
   fi
   for id in ${ids//,/ }; do
      log=$(ls validation/slurm/v22sp_*_"$id".out 2>/dev/null | head -1)
      [[ -n "$log" ]] && echo "   $log: $(grep -c '^RT calculation' "$log") DISORT states, last: $(grep '^RT calculation' "$log" | tail -1)"
   done
fi

echo; echo "== Staged products"
for f in staging/lut_v22_single-r2/terra_modis_m_*_v22.nc; do
   [[ -e "$f" ]] || continue
   echo "$(/usr/bin/sha256sum "$f" | cut -c1-16)  $(stat -c %s "$f") bytes  $f"
done

if [[ "${1:-}" == "--compare" ]]; then
   "$PYTHON" validation/v21_v22_terra_modis/v21_v22_single/compare_full_luts.py
fi
