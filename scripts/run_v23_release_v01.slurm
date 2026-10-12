#!/bin/bash
#SBATCH --job-name=v23
#SBATCH --output=validation/slurm/%x_%j.out
#SBATCH --error=validation/slurm/%x_%j.err
#SBATCH --partition=priority-eodg
#SBATCH --time=48:00:00
#SBATCH --cpus-per-task=1
#SBATCH --mem=24G
#
# SLURM batch script for one V23 Terra MODIS LUT (releases/lut_v23.release:
# the unchanged lut-code-v22.4 scientific code, original Intel single-precision
# DISORT, with the Grid B definitions).  Submitted by scripts/submit_v23_release.sh from
# atmlxint7, which creates the log directory first and drops the login shell's
# TMPDIR.
#
#    sbatch scripts/run_v23_release.slurm runs/<product>_v23.run
#
# The release environment is set here, not inherited.  Run time: the
# definitive V22 products took 4.5-4.6 h on the coarse grids; the V23 smoke
# benchmark (validation/v23_terra_modis/smoke_results.json) scales this to
# about 12 h (liquid water) and 17-26 h (ice) on Grid B, so 48 h leaves a
# margin of about two.  24 GB as for the earlier Grid B productions.  The
# launcher's exit status is the job's, so a failed calculation is never
# reported COMPLETED.

set -euo pipefail

REPO=/home/g/grainger/project-oraclut
PYTHON=/home/g/grainger/miniforge3/envs/science/bin/python
SCRATCH_TMP_ROOT=/network/scratch/grainger/oraclut-tmp
cd "$REPO"

JOB_TMPDIR="${SCRATCH_TMP_ROOT}/job_${SLURM_JOB_ID:-nojob_$$}"
cleanup_tmpdir() {
   if [[ -n "${JOB_TMPDIR:-}" && -d "$JOB_TMPDIR" ]]; then
      case "$JOB_TMPDIR" in "${SCRATCH_TMP_ROOT}"/job_*) rm -rf -- "$JOB_TMPDIR" ;; esac
   fi
}
trap cleanup_tmpdir EXIT
mkdir -p "$JOB_TMPDIR"
export TMPDIR="$JOB_TMPDIR" TMP="$JOB_TMPDIR" TEMP="$JOB_TMPDIR"
export PYTHONPATH="$REPO/src"
export PYTHONUNBUFFERED=1
export ORAC_LUT_RELEASE_MANIFEST=lut_v23.release
export ORAC_LUT_DISORT_PRECISION=single
export ORAC_LUT_DISORT_SP_KERNEL=intel_dlm
unset ORAC_LUT_INTEL_RUNTIME

[[ $# -eq 1 && "$1" == runs/*_v23.run ]] || { echo "usage: sbatch scripts/run_v23_release.slurm runs/<product>_v23.run" >&2; exit 2; }

echo "============================================================"
echo "ORAC V23 LUT SLURM JOB (lut-code-v22.4 scientific code, Grid B, original Intel DISORT)"
echo "Job ID:      ${SLURM_JOB_ID:-unknown}  name: ${SLURM_JOB_NAME:-unknown}"
echo "Node:        $(hostname)"
echo "Start time:  $(date)"
echo "Git:         $(git -C "$REPO" rev-parse HEAD 2>/dev/null || echo unknown)"
echo "Run file:    $1"
echo "Manifest:    $ORAC_LUT_RELEASE_MANIFEST (DISORT $ORAC_LUT_DISORT_PRECISION, kernel $ORAC_LUT_DISORT_SP_KERNEL)"
echo "============================================================"

"$PYTHON" create_orac_luts_release_v04.py "$1"

echo "============================================================"
echo "ORAC V23 LUT SLURM JOB COMPLETED"
echo "End time: $(date)"
echo "============================================================"
