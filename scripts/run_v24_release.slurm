#!/bin/bash
#SBATCH --job-name=v24
#SBATCH --output=validation/slurm/%x_%j.out
#SBATCH --error=validation/slurm/%x_%j.err
#SBATCH --partition=priority-eodg
#SBATCH --time=72:00:00
#SBATCH --cpus-per-task=1
#SBATCH --mem=24G
#
# SLURM batch script for one V24 Terra MODIS LUT (releases/lut_v24.release:
# the lut-code-v24.1 scientific code = the V23 code with the three V24
# refinements of the Mie-class scattering properties, original Intel
# single-precision DISORT, NSTR 60, the Grid B definitions of V23).  Submitted
# by scripts/submit_v24_release.sh from
# atmlxint7, which creates the log directory first and drops the login shell's
# TMPDIR.
#
#    sbatch scripts/run_v24_release.slurm runs/<product>_v24.run
#
# The release environment is set here, not inherited.  Run time
# (validation/v24/REPORT_v24_scattering_refinements.md, section 8;
# validation/v24/scattering_refinements/results/runtime_estimate.json): the
# V23 jobs took 17.1-18.7 h on atmnode009-010; the V24 liquid-water products
# add about 15 h of refined Mie calculation (measured at five wavelengths,
# interpolated over the 37) and about 3% DISORT time (up to 25% per call for
# phase functions with more than 1000 moments), about 32-34 h in all; the Baum
# pagg product is computed as in V23 (about 19 h).  72 h leaves a margin of
# about two.  24 GB as for V23.  The launcher's exit status is
# the job's, so a failed calculation is never reported COMPLETED.

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
export ORAC_LUT_RELEASE_MANIFEST=lut_v24.release
export ORAC_LUT_DISORT_PRECISION=single
export ORAC_LUT_DISORT_SP_KERNEL=intel_dlm
unset ORAC_LUT_INTEL_RUNTIME

[[ $# -eq 1 && "$1" == runs/*_v24.run ]] || { echo "usage: sbatch scripts/run_v24_release.slurm runs/<product>_v24.run" >&2; exit 2; }

echo "============================================================"
echo "ORAC V24 LUT SLURM JOB (lut-code-v24.1: V23 code + V24 scattering refinements, Grid B, original Intel DISORT)"
echo "Job ID:      ${SLURM_JOB_ID:-unknown}  name: ${SLURM_JOB_NAME:-unknown}"
echo "Node:        $(hostname)"
echo "Start time:  $(date)"
echo "Git:         $(git -C "$REPO" rev-parse HEAD 2>/dev/null || echo unknown)"
echo "Run file:    $1"
echo "Manifest:    $ORAC_LUT_RELEASE_MANIFEST (DISORT $ORAC_LUT_DISORT_PRECISION, kernel $ORAC_LUT_DISORT_SP_KERNEL)"
echo "============================================================"

"$PYTHON" create_orac_luts_release_v04.py "$1"

echo "============================================================"
echo "ORAC V24 LUT SLURM JOB COMPLETED"
echo "End time: $(date)"
echo "============================================================"
