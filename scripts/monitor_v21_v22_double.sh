#!/bin/bash
#
# monitor_v21_v22_double.sh
#
# Report the double-precision V21/V22 LUT jobs recorded in
# staging/submissions.tsv (squeue/sacct where available) and the state of the
# six staged products.  With --verify, also verify each staged V22 product
# against the historical V22 product (scripts/verify_staged_lut.py).
# Read-only apart from the optional V22 verification records.
#
# Usage: bash scripts/monitor_v21_v22_double.sh [--verify]

set -u
export PATH=/usr/bin:/bin:$PATH
REPO="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd -P)"
cd "$REPO"
PYTHON=/home/g/grainger/miniforge3/envs/science/bin/python
SHARED=/network/group/aopp/eodg/RGG004_GRAINGER_ORACFILE/ORAC_LUTS

if [[ -f staging/submissions.tsv ]]; then
   echo "== Submissions (staging/submissions.tsv)"
   column -t -s $'\t' staging/submissions.tsv | tail -20
   ids=$(awk -F'\t' 'NR > 1 && $5 ~ /^[0-9]+$/ {print $5}' staging/submissions.tsv | paste -sd, -)
   if [[ -n "$ids" ]] && command -v sacct > /dev/null; then
      echo; echo "== SLURM state"
      sacct -j "$ids" --format=JobID,JobName%16,State,Elapsed,NodeList,ExitCode -X 2>/dev/null
   fi
fi

echo; echo "== Staged products"
for f in staging/lut_v21_double/terra_modis_m_*_v21.nc staging/lut_v22_double/terra_modis_m_*_v22.nc; do
   [[ -e "$f" ]] || continue
   verdict="not verified"
   [[ -f "$f.verification.json" ]] && verdict=$(grep -q '"passed": true' "$f.verification.json" && echo verified || echo FAILED)
   echo "$(/usr/bin/sha256sum "$f" | cut -c1-16)  $(stat -c %s "$f") bytes  $verdict  $f"
done

if [[ "${1:-}" == "--verify" ]]; then
   for f in staging/lut_v22_double/terra_modis_m_*_v22.nc; do
      [[ -e "$f" ]] || continue
      "$PYTHON" scripts/verify_staged_lut.py "$f" "$SHARED/$(basename "$f")" --precision double --record "$f.verification.json"
   done
fi
