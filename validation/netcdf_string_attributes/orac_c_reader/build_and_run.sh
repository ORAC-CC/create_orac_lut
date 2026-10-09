#!/usr/bin/env bash
# Read LUT axis "spacing" attributes with ORAC's own, unmodified
# common/nc_get_string_att.c (validation only).
#
#   validation/netcdf_string_attributes/orac_c_reader/build_and_run.sh ORAC_CHECKOUT FILE...
#
# ORAC_CHECKOUT is a read-only clone of https://github.com/ORAC-CC/orac (the
# campaign used master at eb9a1233c1136ef95fa0c37f650fe519e9134ed4); its C
# file is compiled from where it is, not copied into this repository.  The
# NetCDF headers and library come from the project Conda environment.  The
# build happens in a scratch directory (home-filesystem binaries can refuse to
# execute once a Kerberos ticket has expired).  Exit status: the number of
# attribute reads that failed.
set -euo pipefail
HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd -P)"
ENV="${ORACLUT_CONDA_ENV:-/home/g/grainger/miniforge3/envs/science}"
ORAC="${1:?ORAC checkout path}"; shift
BUILD="${TMPDIR:-/tmp}/orac_c_reader.$$"
mkdir -p "$BUILD"
trap 'rm -rf "$BUILD"' EXIT
env -i PATH=/usr/bin:/bin HOME="$BUILD" TMPDIR="$BUILD" /usr/bin/gcc -O0 -o "$BUILD/orac_get_string_att" \
    "$HERE/main.c" "$ORAC/common/nc_get_string_att.c" -I"$ENV/include" -L"$ENV/lib" -lnetcdf -Wl,-rpath,"$ENV/lib"
echo "ORAC source: $ORAC ($(git -C "$ORAC" rev-parse HEAD 2>/dev/null || echo 'revision unknown'))"
failures=0
for file in "$@"; do
    for axis in optical_depth effective_radius satellite_zenith solar_zenith relative_azimuth; do
        if out="$("$BUILD/orac_get_string_att" "$file" "$axis" spacing 2>&1)"; then
            printf '%s\t%s\tOK\t%s\n' "$(basename "$file")" "$axis" "$out"
        else
            printf '%s\t%s\tFAIL\t%s\n' "$(basename "$file")" "$axis" "$out"
            failures=$((failures + 1))
        fi
    done
done
echo "failed reads: $failures"
exit "$failures"
