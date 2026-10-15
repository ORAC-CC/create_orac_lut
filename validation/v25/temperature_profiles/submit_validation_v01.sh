#!/bin/bash
# Submit the compact V24/V25 cloud-temperature validation (read-only for the
# code and the products) as single-CPU SLURM jobs on priority-eodg, from
# atmlxint7 (dispatched over ssh when run elsewhere).
#
#    bash validation/v25/temperature_profiles/submit_validation.sh [task-prefix]
#
# Outputs: validation/v25/temperature_profiles/work/ (NetCDF, git-ignored);
# logs: validation/slurm/v25_validation/.  Each V25 job prints the SHA-256 of
# the V25 source files it loads (compact_lut.py); the v24 jobs load the sealed
# releases/lut_v24.release files.  A failed task exits non-zero.
set -euo pipefail
REPO=/home/g/grainger/project-oraclut
if [[ "$(hostname -s)" != atmlxint7 ]]; then
   exec ssh atmlxint7 "bash $REPO/validation/v25/temperature_profiles/submit_validation.sh ${1:-}"
fi
cd "$REPO"
HERE=validation/v25/temperature_profiles
LOGS=validation/slurm/v25_validation
mkdir -p "$LOGS" "$HERE/work"
PREFIX="${1:-}"
PY=/home/g/grainger/miniforge3/envs/science/bin/python
ENVS="export ORAC_LUT_DISORT_PRECISION=single ORAC_LUT_DISORT_SP_KERNEL=intel_dlm OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 PYTHONUNBUFFERED=1; unset ORAC_LUT_INTEL_RUNTIME; export TMPDIR=/network/scratch/grainger/oraclut-tmp/v25val_\$SLURM_JOB_ID; mkdir -p \$TMPDIR; trap 'rm -rf \$TMPDIR' EXIT"
CH="1 6 7 20 31 32"

tasks=()
for model in liquid-water_240 water-ice_agg; do
   for variant in v24 isothermal profile constant mirror ttop240 layers20 layers46; do
      tasks+=("${model}_${variant}|$PY $HERE/compact_lut.py $model $variant $HERE/work/${model}_${variant} $CH")
   done
done
for variant in v24 profile; do
   tasks+=("liquid-water_old_${variant}|$PY $HERE/compact_lut.py liquid-water_old $variant $HERE/work/liquid-water_old_${variant} $CH")
done

for task in "${tasks[@]}"; do
   name="${task%%|*}"; command="${task#*|}"
   [[ -n "$PREFIX" && "$name" != "$PREFIX"* ]] && continue
   job=$(sbatch --parsable --job-name="v25val_$name" --partition=priority-eodg --cpus-per-task=1 --mem=3G \
         --time=6:00:00 --output="$LOGS/%x_%j.out" --error="$LOGS/%x_%j.err" \
         --wrap="cd $REPO && $ENVS && echo \"host \$(hostname) git \$(git rev-parse --short HEAD)\" && $command")
   printf '%s\t%s\t%s\t%s\n' "$(date -Iseconds)" "$job" "$name" "$command" >> "$LOGS/submissions.tsv"
   echo "$job $name"
done
