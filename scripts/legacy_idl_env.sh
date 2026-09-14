#!/usr/bin/env bash
# Establish the verified legacy ORAC IDL environment for the current shell.
#
# This file is intentionally meant to be sourced by a short-lived developer
# shell (the reference runner sources it inside a subshell).  It does not edit
# shell startup files or any persistent module configuration.

if ! type module >/dev/null 2>&1; then
    echo "legacy_idl_env.sh: the AOPP module command is unavailable" >&2
    return 1 2>/dev/null || exit 1
fi

module load intel-compilers/2022
module load idl/890

oraclut_root="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"

# IDL source search is recursive; DLM search is deliberately explicit and
# non-recursive so the DISORT4 implementation cannot be exposed accidentally.
export IDL_PATH="<IDL_DEFAULT_PATH>:+${oraclut_root}/"
export IDL_DLM_PATH="<IDL_DEFAULT_DLM>:${oraclut_root}/mie/dlm-code:${oraclut_root}/create_orac_lut/disort2/src"
unset IDL_STARTUP

echo "Legacy ORAC environment established for this shell"
echo "  IDL: $(command -v idl || echo unavailable)"
echo "  IDL_PATH: ${IDL_PATH}"
echo "  IDL_DLM_PATH: ${IDL_DLM_PATH}"
