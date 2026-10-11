#!/bin/bash
#SBATCH --job-name=v22sp
#SBATCH --output=validation/slurm/%x_%j.out
#SBATCH --error=validation/slurm/%x_%j.err
#SBATCH --partition=priority-eodg
#SBATCH --time=24:00:00
#SBATCH --cpus-per-task=1
#SBATCH --mem=8G
#
# SLURM batch script for one corrected single-precision V22 Terra MODIS LUT
# (releases/lut_v22_single.release: IDL-compatible arithmetic, original Intel
# DISORT kernel).  Submitted by scripts/submit_v22_single.sh from atmlxint7,
# which creates the log directory first and drops the login shell's TMPDIR.
#
#    sbatch scripts/run_v22_single.slurm runs/<product>_v22_single.run
#
# The release environment is set here, not inherited.  Run time: the
# single-precision V22 products of lut-code-v22.1 took 6.7-7.1 h each on
# atmnode009-010; the corrected release uses the same solver workload, so 24 h
# leaves a threefold margin.  The launcher's exit status is the job's, so a
# failed calculation is never reported COMPLETED.

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
export ORAC_LUT_RELEASE_MANIFEST=lut_v22_single.release
export ORAC_LUT_DISORT_PRECISION=single
export ORAC_LUT_DISORT_SP_KERNEL=intel_dlm
unset ORAC_LUT_INTEL_RUNTIME

[[ $# -eq 1 && "$1" == runs/*_v22_single.run ]] || { echo "usage: sbatch scripts/run_v22_single.slurm runs/<product>_v22_single.run" >&2; exit 2; }

echo "============================================================"
echo "ORAC V22 LUT SLURM JOB (corrected single precision, original Intel DISORT)"
echo "Job ID:      ${SLURM_JOB_ID:-unknown}  name: ${SLURM_JOB_NAME:-unknown}"
echo "Node:        $(hostname)"
echo "Start time:  $(date)"
echo "Git:         $(git -C "$REPO" rev-parse HEAD 2>/dev/null || echo unknown)"
echo "Run file:    $1"
echo "Manifest:    $ORAC_LUT_RELEASE_MANIFEST (DISORT $ORAC_LUT_DISORT_PRECISION, kernel $ORAC_LUT_DISORT_SP_KERNEL)"
echo "============================================================"

"$PYTHON" create_orac_luts_release_v04.py "$1"

echo "============================================================"
echo "ORAC V22 LUT SLURM JOB COMPLETED"
echo "End time: $(date)"
echo "============================================================"
