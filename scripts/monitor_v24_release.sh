#!/bin/bash
#
# monitor_v24_release.sh
#
# Report the V24 LUT jobs (v24_*) recorded in
# staging/submissions.tsv and the products staged in staging/lut_v24/.
# With --compare, compare each staged product with the V23 product on the same
# Grid B (validation/v24/scattering_refinements/compare_full_products.py).
# Read-only apart from the comparison results.
#
# Usage: bash scripts/monitor_v24_release.sh [--compare]

set -u
export PATH=/usr/bin:/bin:$PATH
REPO="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd -P)"
cd "$REPO"
PYTHON=/home/g/grainger/miniforge3/envs/science/bin/python

if [[ -f staging/submissions.tsv ]]; then
   echo "== Submissions (v24_*)"
   awk -F'\t' 'NR == 1 || $4 ~ /^v24_/' staging/submissions.tsv | column -t -s $'\t'
   ids=$(awk -F'\t' '$4 ~ /^v24_/ && $5 ~ /^[0-9]+$/ {print $5}' staging/submissions.tsv | paste -sd, -)
   if [[ -n "$ids" ]] && command -v sacct > /dev/null; then
      echo; echo "== SLURM state"
      sacct -j "$ids" --format=JobID,JobName%16,State,Elapsed,NodeList,ExitCode -X 2>/dev/null
   fi
   for id in ${ids//,/ }; do
      log=$(ls validation/slurm/v24_*_"$id".out 2>/dev/null | head -1)
      [[ -n "$log" ]] && echo "   $log: $(grep -c '^RT calculation' "$log") DISORT states, last: $(grep '^RT calculation' "$log" | tail -1)"
   done
fi

echo; echo "== Staged products"
for f in staging/lut_v24/terra_modis_m_*_v24.nc; do
   [[ -e "$f" ]] || continue
   echo "$(/usr/bin/sha256sum "$f" | cut -c1-16)  $(stat -c %s "$f") bytes  $f"
done

if [[ "${1:-}" == "--compare" ]]; then
   "$PYTHON" validation/v24/scattering_refinements/compare_full_products.py
fi
