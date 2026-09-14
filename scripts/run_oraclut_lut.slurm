#!/bin/bash
#SBATCH --job-name=oraclut_lut
#SBATCH --output=validation/slurm/%x_%j.out
#SBATCH --error=validation/slurm/%x_%j.err
#SBATCH --time=04:00:00
#SBATCH --cpus-per-task=1
#SBATCH --mem=8G
#
# Batch script for the Python ORAC LUT generator.
#
# Submit through scripts/submit_oraclut_lut.sh from atmlxint7; it creates the
# log directory named above before sbatch opens stdout/stderr (a mkdir inside
# this script would be too late) and removes the login shell's per-user
# TMPDIR from the submitted environment so slurmstepd does not try to create
# it on the compute node.  Resources may be overridden on the sbatch command
# line, e.g. --time=12:00:00 --mem=16G, which take precedence over the
# #SBATCH defaults above.
#
# All arguments are passed unchanged to python -m oraclut.generate, e.g.
#   --config configs/meteosat-10_seviri_cloud_liquid-water_stg_test.driver

set -euo pipefail

REPO=/home/g/grainger/project-oraclut
PYTHON=/home/g/grainger/miniforge3/envs/science/bin/python
SCRATCH_TMP_ROOT=/network/scratch/grainger/oraclut-tmp

cd "$REPO"
export PYTHONPATH="$REPO/src"

# ---------------------------------------------------------------------------
# Job-specific temporary directory on scratch (not /tmp on the compute node).
# ---------------------------------------------------------------------------
JOB_TMPDIR="${SCRATCH_TMP_ROOT}/job_${SLURM_JOB_ID:-nojob_$$}"
cleanup_tmpdir() {
    # Only ever remove the directory this job created.
    if [[ -n "${JOB_TMPDIR:-}" && -d "$JOB_TMPDIR" ]]; then
        case "$JOB_TMPDIR" in
            "${SCRATCH_TMP_ROOT}"/job_*) rm -rf -- "$JOB_TMPDIR" ;;
        esac
    fi
}
trap cleanup_tmpdir EXIT
mkdir -p "$JOB_TMPDIR"
if [[ ! -w "$JOB_TMPDIR" ]]; then
    echo "ERROR: job temporary directory is not writable: $JOB_TMPDIR" >&2
    exit 1
fi
export TMPDIR="$JOB_TMPDIR"
export TMP="$JOB_TMPDIR"
export TEMP="$JOB_TMPDIR"

echo "============================================================"
echo "ORAC LUT SLURM JOB"
echo "============================================================"
echo "Job ID:        ${SLURM_JOB_ID:-unknown}"
echo "Job name:      ${SLURM_JOB_NAME:-unknown}"
echo "Node:          $(hostname)"
echo "Start time:    $(date)"
echo "Repository:    $REPO"
echo "Python:        $PYTHON"
echo "Python version: $("$PYTHON" --version 2>&1)"
echo "TMPDIR:        $TMPDIR"
echo "CPUs:          ${SLURM_CPUS_PER_TASK:-unknown}; memory: ${SLURM_MEM_PER_NODE:-unknown} MB; time limit: see squeue -j ${SLURM_JOB_ID:-}"
echo "============================================================"

if [ "$#" -eq 0 ]; then
    echo "ERROR: no ORAC LUT arguments supplied." >&2
    echo >&2
    echo "Examples:" >&2
    echo "  scripts/submit_oraclut_lut.sh --config configs/meteosat-10_seviri_cloud_liquid-water_stg_test.driver" >&2
    echo "  scripts/submit_oraclut_lut.sh --platform meteosat-10 --instrument seviri --microphysics liquid-water_stg.mm --test --channels 1" >&2
    exit 1
fi

echo "Arguments:"
printf '  %q\n' "$@"
echo "============================================================"

"$PYTHON" -m oraclut.generate "$@"

echo "============================================================"
echo "ORAC LUT SLURM JOB COMPLETED"
echo "End time: $(date)"
echo "============================================================"
