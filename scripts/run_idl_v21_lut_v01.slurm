#!/bin/bash
#SBATCH --job-name=v21dp
#SBATCH --output=validation/slurm/%x_%j.out
#SBATCH --error=validation/slurm/%x_%j.err
#SBATCH --partition=priority-eodg
#SBATCH --time=10:00:00
#SBATCH --cpus-per-task=1
#SBATCH --mem=8G
#
# SLURM batch script for the double-precision V21 IDL LUTs.  Submitted by
# scripts/submit_v21_v22_double.sh from atmlxint7 (which creates the log
# directory first and drops the login shell's TMPDIR).  Arguments are passed
# unchanged to scripts/run_idl_v21_lut.sh:
#    validate <out_dir>
#    product <pold|p240|pagg> <stage_dir>
# A failure of the worker is the job's exit status, so a failed calculation
# is never reported COMPLETED.

set -euo pipefail

REPO=/home/g/grainger/project-oraclut
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

echo "============================================================"
echo "ORAC V21 IDL LUT SLURM JOB (double-precision DISORT)"
echo "Job ID:      ${SLURM_JOB_ID:-unknown}  name: ${SLURM_JOB_NAME:-unknown}"
echo "Node:        $(hostname)"
echo "Start:       $(date --iso-8601=seconds)"
echo "Git:         $(git -C "$REPO" rev-parse HEAD 2>/dev/null || echo unknown)"
echo "Arguments:   $*"
echo "============================================================"

bash "$REPO/scripts/run_idl_v21_lut.sh" "$@"

echo "============================================================"
echo "ORAC V21 IDL LUT SLURM JOB COMPLETED $(date --iso-8601=seconds)"
echo "============================================================"
