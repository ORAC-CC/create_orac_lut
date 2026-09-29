#!/usr/bin/env bash
# Submit a Python ORAC LUT generation job to SLURM.
#
# Usage, from a normal interactive AOPP host (atmlxint5 is the usual one):
#   scripts/submit_oraclut_lut.sh runs/<run-file>.run [--time=.. --mem=.. --cpus-per-task=..]
#   scripts/submit_oraclut_lut.sh --config configs/<configuration>.driver [--time=.. --mem=.. --cpus-per-task=..]
#   scripts/submit_oraclut_lut.sh [sbatch options] -- <oraclut.generate options>
#
# A runs/<name>.run file is executed by python create_orac_luts.py (the
# IDL-structured path); --config and other options go to python -m
# oraclut.generate (the current path).  Arguments beginning with --time=,
# --mem=, --cpus-per-task=, --job-name=, --partition=, --account= or --qos= are
# passed to sbatch.  A literal "--" ends the sbatch options.
#
# Where the work happens:
#   atmlxint7           validates, creates the log directory, drops the login
#                       shell's TMPDIR and calls sbatch directly.
#   atmlxint1/2/3/5/6   validate locally, then ssh to atmlxint7 and run this
#                       same script there with the arguments preserved exactly;
#                       the calculation runs on a SLURM compute node, never on
#                       an interactive host.
#   anywhere else       refused (atmlxint4 is reserved; compute nodes never submit).
#
# If ssh fails because your Kerberos ticket has expired, run `kinit` and repeat
# the same command.  Nothing is copied between hosts: the repository is on the
# shared filesystem.

set -euo pipefail

REPO_ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd -P)"
BATCH_SCRIPT="scripts/run_oraclut_lut.slurm"
LOG_DIR="validation/slurm"
SUBMIT_HOST="atmlxint7"
PYTHON="/home/g/grainger/miniforge3/envs/science/bin/python"
# Interactive hosts permitted for ORAC/POM work by AGENTS.md; they dispatch to
# SUBMIT_HOST over ssh.  atmlxint4 is deliberately absent.
DISPATCH_HOSTS=(atmlxint1 atmlxint2 atmlxint3 atmlxint5 atmlxint6)

# Test hooks: alternative sbatch/ssh commands for shell-level tests of this
# wrapper (scripts that echo their arguments).  Not used operationally.
SBATCH="${ORACLUT_SBATCH:-sbatch}"
SSH="${ORACLUT_SSH:-ssh}"

die() { echo "submit_oraclut_lut.sh: $*" >&2; exit 2; }

usage() { sed -n '2,26p' "${BASH_SOURCE[0]}" | sed 's/^# \{0,1\}//'; }

if [[ "$#" -eq 0 || "${1:-}" == "-h" || "${1:-}" == "--help" ]]; then
    usage
    exit 0
fi

ORIGINAL_ARGS=("$@")

# ---------------------------------------------------------------------------
# 1. Validate locally (on every host) so obvious mistakes fail before any hop.
# ---------------------------------------------------------------------------
cd "$REPO_ROOT"

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
[[ "${#gen_args[@]}" -gt 0 ]] || die "no generator arguments given (e.g. runs/<run-file>.run or --config configs/<configuration>.driver)"

# A concrete run file or the suffixless maker specification must exist and be
# the only generator argument.
for argument in "${gen_args[@]}"; do
    case "$argument" in
        makerunfile)
            [[ -f "$argument" ]] || die "master specification not found: $argument (relative to $REPO_ROOT)"
            [[ "${#gen_args[@]}" -eq 1 ]] || die "the master specification takes no further generator arguments: ${gen_args[*]}"
            ;;
        *.run)
            [[ -f "$argument" ]] || die "run file not found: $argument (relative to $REPO_ROOT)"
            [[ "${#gen_args[@]}" -eq 1 ]] || die "a run file takes no further generator arguments: ${gen_args[*]}"
            ;;
    esac
done

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

# A maker-style master run expands to one validated SLURM array. The manifest
# freezes the exact concrete configurations checked before submission.
master_run=""
if [[ "${#gen_args[@]}" -eq 1 && "$(basename -- "${gen_args[0]}")" == "makerunfile" ]] && grep -qE '^[[:space:]]*mmfiles[[:space:]]*=' "${gen_args[0]}"; then
    master_run="${gen_args[0]}"
fi

# ---------------------------------------------------------------------------
# 2. Decide by host.  ORACLUT_DISPATCHED=1 marks an invocation that arrived by
#    ssh from a dispatch host; such an invocation must find itself on
#    SUBMIT_HOST and go straight to sbatch, never dispatch again.
# ---------------------------------------------------------------------------
host="$(hostname -s 2>/dev/null || hostname)"
host="${host%%.*}"

is_dispatch_host() {
    local candidate
    for candidate in "${DISPATCH_HOSTS[@]}"; do
        [[ "$1" == "$candidate" ]] && return 0
    done
    return 1
}

if [[ "$host" != "$SUBMIT_HOST" ]]; then
    if [[ "${ORACLUT_DISPATCHED:-0}" == "1" ]]; then
        die "dispatched invocation landed on '${host}', not ${SUBMIT_HOST}; refusing to dispatch again or to run anything here"
    fi
    is_dispatch_host "$host" \
        || die "SLURM submission is only permitted from ${SUBMIT_HOST} (directly) or from ${DISPATCH_HOSTS[*]} (dispatched over ssh); this host is '${host}'"

    # Quote every original argument for the remote bash so spaces and shell
    # metacharacters survive exactly; no unquoted command string is built.
    remote_command="cd $(printf '%q' "$REPO_ROOT") && ORACLUT_DISPATCHED=1 exec scripts/submit_oraclut_lut.sh"
    for argument in "${ORIGINAL_ARGS[@]}"; do
        remote_command+=" $(printf '%q' "$argument")"
    done

    echo "Dispatching ORAC LUT submission from ${host} to ${SUBMIT_HOST} over ssh"
    echo "  repository: $REPO_ROOT (shared filesystem; nothing is copied)"
    printf '  arguments: %q\n' "${ORIGINAL_ARGS[*]}"

    # BatchMode: never prompt for a password; an expired Kerberos ticket fails
    # immediately and is reported below instead.
    ssh_stderr="$(mktemp)"
    trap 'rm -f "$ssh_stderr"' EXIT
    set +e
    "$SSH" -o BatchMode=yes -o ConnectTimeout=20 "$SUBMIT_HOST" bash -lc "$(printf '%q' "$remote_command")" 2> >(tee "$ssh_stderr" >&2)
    status=$?
    set -e
    if [[ "$status" -eq 255 ]] || grep -qiE 'permission denied|gssapi|kerberos|credentials|no ticket|publickey' "$ssh_stderr"; then
        echo >&2
        echo "submit_oraclut_lut.sh: ssh to ${SUBMIT_HOST} failed (exit ${status}). This is usually an expired Kerberos ticket." >&2
        echo "  Run:    kinit" >&2
        printf '  Retry:  scripts/submit_oraclut_lut.sh %s\n' "$(printf '%q ' "${ORIGINAL_ARGS[@]}")" >&2
        exit 3
    fi
    exit "$status"
fi

# ---------------------------------------------------------------------------
# 3. On atmlxint7: direct sbatch (the only place sbatch is ever called).
# ---------------------------------------------------------------------------
# Create the log directory BEFORE sbatch: SLURM opens the #SBATCH
# --output/--error files before the batch script runs.
mkdir -p "$LOG_DIR"
[[ -w "$LOG_DIR" ]] || die "log directory is not writable: $LOG_DIR"

# The login shell exports a per-user TMPDIR (/tmp/user/<uid>) that does not
# exist on compute nodes; slurmstepd would try to create it before the batch
# script starts and emit "Unable to create TMPDIR ... Permission denied".
# Drop it from the submitted environment; the batch script sets its own
# job-specific TMPDIR on scratch.
unset TMPDIR TMP TEMP

echo "Submitting ORAC LUT job from $host"
echo "  repository: $REPO_ROOT"
echo "  batch script: $BATCH_SCRIPT"
echo "  logs: $LOG_DIR/oraclut_lut_<arrayjob>_<task>.out|.err"
if [[ "${#sbatch_args[@]}" -gt 0 ]]; then
    printf '  sbatch override: %q\n' "${sbatch_args[@]}"
fi
printf '  generator arguments: %q\n' "${gen_args[*]}"

if [[ -n "$master_run" ]]; then
    [[ "${#gen_args[@]}" -eq 1 ]] || die "a master run takes no further generator arguments"
    [[ -f "$PYTHON" ]] || die "configured Python interpreter not found: $PYTHON"
    export PYTHONPATH="$REPO_ROOT/src"
    "$PYTHON" -m oraclut.master_run --preflight "$master_run"
    count="$($PYTHON -m oraclut.master_run --count "$master_run")"
    [[ "$count" =~ ^[1-9][0-9]*$ ]] || die "master expansion returned invalid count: $count"
    manifest="$LOG_DIR/manifests/oraclut_master_$$.json"
    "$PYTHON" -m oraclut.master_run --write-manifest "$manifest" "$master_run"
    echo "Array tasks: 0-$((count - 1))"
    echo "Manifest: $manifest"
    "$SBATCH" "${sbatch_args[@]+${sbatch_args[@]}}" --array="0-$((count - 1))" "$BATCH_SCRIPT" --master-manifest "$manifest"
else
    "$SBATCH" "${sbatch_args[@]+${sbatch_args[@]}}" "$BATCH_SCRIPT" "${gen_args[@]}"
fi
