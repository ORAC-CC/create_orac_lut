#!/usr/bin/env bash
# Submit a Python ORAC LUT generation job to SLURM.
#
# Usage (on atmlxint7 only):
#   scripts/submit_oraclut_lut.sh --config configs/<configuration>.driver
#   scripts/submit_oraclut_lut.sh [sbatch options] -- <oraclut.generate options>
#
# Arguments beginning with --time=, --mem=, --cpus-per-task=, --job-name= or
# --partition= are passed to sbatch; everything else goes to
# python -m oraclut.generate.  A literal "--" ends the sbatch options.
#
# The wrapper does not ssh anywhere: log into atmlxint7 first.

set -euo pipefail

REPO_ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd -P)"
BATCH_SCRIPT="scripts/run_oraclut_lut.slurm"
LOG_DIR="validation/slurm"
SUBMIT_HOST="atmlxint7"

# ORACLUT_SBATCH may name an alternative sbatch command for shell-level tests
# of this wrapper (for example a script that echoes its arguments); it is not
# used operationally.
SBATCH="${ORACLUT_SBATCH:-sbatch}"

die() { echo "submit_oraclut_lut.sh: $*" >&2; exit 2; }

usage() {
    sed -n '2,12p' "${BASH_SOURCE[0]}" | sed 's/^# \{0,1\}//'
}

if [[ "$#" -eq 0 || "${1:-}" == "-h" || "${1:-}" == "--help" ]]; then
    usage
    exit 0
fi

# 1. Host check: SLURM jobs may only be submitted from atmlxint7 (AGENTS.md).
host="$(hostname -s 2>/dev/null || hostname)"
if [[ "${host%%.*}" != "$SUBMIT_HOST" ]]; then
    die "SLURM submission is only permitted from ${SUBMIT_HOST}; this host is '${host}'. Log into ${SUBMIT_HOST} and re-run the same command."
fi

# 2. Work from the repository root so relative #SBATCH log paths resolve.
cd "$REPO_ROOT"

# 3. Split sbatch overrides from generator arguments.
sbatch_args=()
gen_args=()
while [[ "$#" -gt 0 ]]; do
    case "$1" in
        --)
            shift
            gen_args+=("$@")
            break
            ;;
        --time=*|--mem=*|--cpus-per-task=*|--job-name=*|--partition=*|--account=*|--qos=*)
            sbatch_args+=("$1")
            ;;
        --time|--mem|--cpus-per-task|--job-name|--partition|--account|--qos)
            [[ "$#" -ge 2 ]] || die "$1 requires a value"
            sbatch_args+=("$1" "$2")
            shift
            ;;
        *)
            gen_args+=("$1")
            ;;
    esac
    shift
done
[[ "${#gen_args[@]}" -gt 0 ]] || die "no generator arguments given (e.g. --config configs/<configuration>.driver)"

# 4. Verify any configuration file named with --config exists.
for ((i = 0; i < ${#gen_args[@]}; i++)); do
    case "${gen_args[i]}" in
        --config=*) config="${gen_args[i]#--config=}" ;;
        --config)
            (( i + 1 < ${#gen_args[@]} )) || die "--config requires a file name"
            config="${gen_args[i+1]}"
            ;;
        *) continue ;;
    esac
    [[ -f "$config" ]] || die "configuration file not found: $config (relative to $REPO_ROOT)"
done

[[ -f "$BATCH_SCRIPT" ]] || die "batch script missing: $BATCH_SCRIPT"

# 5. Create the log directory BEFORE sbatch: SLURM opens the #SBATCH
#    --output/--error files before the batch script runs.
mkdir -p "$LOG_DIR"
[[ -w "$LOG_DIR" ]] || die "log directory is not writable: $LOG_DIR"

# 6. The login shell exports a per-user TMPDIR (/tmp/user/<uid>) that does
#    not exist on compute nodes; slurmstepd would try to create it before the
#    batch script starts and emit "Unable to create TMPDIR ... Permission
#    denied".  Drop it from the submitted environment; the batch script sets
#    its own job-specific TMPDIR on scratch.
unset TMPDIR TMP TEMP

echo "Submitting ORAC LUT job from $host"
echo "  repository: $REPO_ROOT"
echo "  batch script: $BATCH_SCRIPT"
echo "  logs: $LOG_DIR/oraclut_lut_<jobid>.out|.err"
if [[ "${#sbatch_args[@]}" -gt 0 ]]; then
    printf '  sbatch override: %q\n' "${sbatch_args[@]}"
fi
printf '  generator arguments: %q\n' "${gen_args[*]}"

# 7. Submit and show sbatch's own output (job ID).
"$SBATCH" "${sbatch_args[@]+"${sbatch_args[@]}"}" "$BATCH_SCRIPT" "${gen_args[@]}"
