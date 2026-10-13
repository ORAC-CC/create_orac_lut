#!/bin/bash
# Submit the compact V24 scattering-refinement validation (read-only for the
# code and the products) as single-CPU SLURM jobs on priority-eodg, from
# atmlxint7 (dispatched over ssh when run elsewhere).
#
#    bash validation/v24/scattering_refinements/submit_validation.sh [task-prefix]
#
# Outputs: validation/v24/scattering_refinements/work/ (NetCDF/npz, git-ignored)
# and results/ (JSON); logs: validation/slurm/v24_validation/.  Each job prints
# the SHA-256 of the V24 source files it loads (codebase.py); the V23 jobs load
# the sealed releases/lut_v23.release files.  A failed task exits non-zero.
set -euo pipefail
REPO=/home/g/grainger/project-oraclut
if [[ "$(hostname -s)" != atmlxint7 ]]; then
   exec ssh atmlxint7 "bash $REPO/validation/v24/scattering_refinements/submit_validation.sh ${1:-}"
fi
cd "$REPO"
HERE=validation/v24/scattering_refinements
LOGS=validation/slurm/v24_validation
mkdir -p "$LOGS" "$HERE/work" "$HERE/results"
PREFIX="${1:-}"
PY=/home/g/grainger/miniforge3/envs/science/bin/python
ENVS="export ORAC_LUT_DISORT_PRECISION=single ORAC_LUT_DISORT_SP_KERNEL=intel_dlm OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 PYTHONUNBUFFERED=1; unset ORAC_LUT_INTEL_RUNTIME; export TMPDIR=/network/scratch/grainger/oraclut-tmp/v24val_\$SLURM_JOB_ID; mkdir -p \$TMPDIR; trap 'rm -rf \$TMPDIR' EXIT"
CH="1 2 3 6 7 20 31 32"

tasks=()
for code in v23 v24; do
   for model in liquid-water_240 liquid-water_old; do
      tasks+=("optics_${code}_${model}|$PY $HERE/optics.py $code $model $HERE/work/optics_${code}_${model}.npz --channels $CH --efr 1 5 11 21 39")
      tasks+=("compact_${code}_${model}|$PY $HERE/compact_lut.py $code $model $HERE/work/compact_${code}_${model} $CH")
   done
   tasks+=("optics_${code}_water-ice_sph|$PY $HERE/optics.py $code water-ice_sph $HERE/work/optics_${code}_water-ice_sph.npz --channels $CH --efr 5 41 93")
   tasks+=("compact_${code}_water-ice_agg|$PY $HERE/compact_lut.py $code water-ice_agg $HERE/work/compact_${code}_water-ice_agg $CH")
done
for case in liquid-water_240:0.469:5 liquid-water_240:0.646:5 liquid-water_240:0.646:21 liquid-water_240:0.646:39 \
            liquid-water_240:0.858:11 liquid-water_240:2.13:21 liquid-water_240:3.79:1 liquid-water_240:3.79:21 \
            liquid-water_240:11.03:1 liquid-water_240:11.03:21 water-ice_sph:0.646:9 water-ice_sph:0.646:41 \
            water-ice_sph:2.13:41 water-ice_sph:11.03:93; do
   tasks+=("convergence_${case//:/_}|$PY $HERE/convergence.py $case $HERE/results/convergence_${case//:/_}.json")
done
for case in "0.646 39" "0.646 21" "0.469 21" "0.858 11" "3.79 21"; do
   set -- $case
   tasks+=("nmom_liquid-water_240_$1_$2|$PY $HERE/disort_nmom.py liquid-water_240 $1 $2 $HERE/results/disort_nmom_liquid-water_240_$1_$2.json")
done
for wl in 0.469 0.646 1.64 3.79 11.03; do
   tasks+=("timing_liquid-water_240_$wl|$PY $HERE/timing.py liquid-water_240 $wl $HERE/results/timing_liquid-water_240_$wl.json")
done

for task in "${tasks[@]}"; do
   name="${task%%|*}"; command="${task#*|}"
   [[ -n "$PREFIX" && "$name" != "$PREFIX"* ]] && continue
   job=$(sbatch --parsable --job-name="v24val_$name" --partition=priority-eodg --cpus-per-task=1 --mem=3G \
         --time=36:00:00 --output="$LOGS/%x_%j.out" --error="$LOGS/%x_%j.err" \
         --wrap="cd $REPO && $ENVS && echo \"host \$(hostname) git \$(git rev-parse --short HEAD)\" && $command")
   printf '%s\t%s\t%s\t%s\n' "$(date -Iseconds)" "$job" "$name" "$command" >> "$LOGS/submissions.tsv"
   echo "$job $name"
done
