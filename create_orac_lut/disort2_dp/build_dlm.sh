#!/bin/bash
#
# build_dlm.sh
#
# Build the double-precision DISORT used by the ORAC LUT generators:
#
#   DISORTtoIDL_DP.so  IDL DLM (procedure DISORT_DP, function DISORT_DP_INFO);
#                      the default solver of create_orac_lut/call_disort.pro
#   libdisort_dp.so    the same kernel objects as a plain shared library
#                      exporting disort_, loaded by the Python V22 generator
#                      (src/oraclut/radiative_transfer/legacy_disort.py)
#
# The kernel is the unmodified production DISORT source of
# create_orac_lut/disort2/src compiled with default REAL promoted to 8 bytes.
# The historical single-precision DISORTtoIDL.so is not touched.
# Compiler, flags and sources are exactly those validated in
# validation/v21_v22_terra_modis/REPORT_idl_double_precision.md; the only
# additions are the RPATH to the Intel runtime (so the libraries load without
# the module) and the in-place output directory.
#
# Usage: bash create_orac_lut/disort2_dp/build_dlm.sh   (writes only this directory and ./build)

set -euo pipefail

HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
SRC="$(cd "$HERE/../disort2/src" && pwd)"
OBJ="$HERE/build"
IDL_INCLUDE=/network/software/linux-x86_64/idl/890/idl89/external/include

if ! type module > /dev/null 2>&1; then source /etc/profile.d/modules.sh; fi
module load intel-compilers/2022
IFORT="$(command -v ifort)"
IFORT_DIR="$(dirname "$IFORT")"
INTEL_LIB="$(cd "$IFORT_DIR/../../compiler/lib/intel64_lin" && pwd)"
# Minimal environment: the Kerberos-protected home must not be on PATH.
RUN="env -i PATH=$IFORT_DIR:/usr/bin:/bin HOME=/tmp LD_LIBRARY_PATH=${LD_LIBRARY_PATH:-} INTEL_LICENSE_FILE=${INTEL_LICENSE_FILE:-}"

FFLAGS="-O2 -fp-model precise -fpic -real-size 64 -double-size 64 -integer-size 32"
CFLAGS="-O2 -fPIC -std=gnu99 -Wall -Wextra -Wno-cast-function-type"
LDFLAGS="-shared -fpic -Wl,-rpath,$INTEL_LIB"
FORTRAN_SOURCES=("$SRC/DISORTfunctions.f" "$SRC/PRTFIN.f" "$SRC/disort/ErrPack.f" "$SRC/disort/LINPAK.f"
                 "$SRC/disort/RDI1MACH.f" "$HERE/dp_support.f")

rm -rf "$OBJ"
mkdir -p "$OBJ"
IFORT_VERSION="$($RUN "$IFORT" --version 2>&1 | head -1)"
GCC_VERSION="$(/usr/bin/gcc --version | head -1)"
INFO="ifort=[$IFORT_VERSION] FFLAGS=[$FFLAGS] gcc=[$GCC_VERSION] CFLAGS=[$CFLAGS]"

objects=()
for source in "${FORTRAN_SOURCES[@]}"; do
   object="$OBJ/$(basename "${source%.*}").o"
   $RUN "$IFORT" $FFLAGS -c "$source" -o "$object" 2> "$object.log" || { cat "$object.log"; exit 1; }
   objects+=("$object")
done
$RUN /usr/bin/gcc $CFLAGS -I"$IDL_INCLUDE" -DDP_BUILD_INFO="\"$INFO\"" -c "$HERE/DISORTtoIDL_DP.c" -o "$OBJ/DISORTtoIDL_DP.o"

$RUN "$IFORT" $LDFLAGS -Wl,--version-script="$HERE/exports.map" -Wl,-Bsymbolic -Wl,--no-undefined-version \
   -o "$HERE/DISORTtoIDL_DP.so" "$OBJ/DISORTtoIDL_DP.o" "${objects[@]}"
$RUN "$IFORT" $LDFLAGS -o "$HERE/libdisort_dp.so" "${objects[@]}"

exported="$(/usr/bin/nm -D --defined-only "$HERE/DISORTtoIDL_DP.so" | awk '{print $3}' | sort | tr '\n' ' ')"
case " $exported " in
   *" disort_ "*|*" getmom_ "*|*" r1mach_ "*) echo "kernel symbols leaked from DLM: $exported" >&2; exit 1 ;;
esac
case " $exported " in *" IDL_Load "*) ;; *) echo "IDL_Load not exported" >&2; exit 1 ;; esac
/usr/bin/nm -D --defined-only "$HERE/libdisort_dp.so" | grep -q ' T disort_$' || { echo "disort_ not exported by libdisort_dp.so" >&2; exit 1; }

{
   echo "date: $(date --iso-8601=seconds)"
   echo "host: $(hostname)"
   echo "ifort: $IFORT_VERSION ($IFORT)"
   echo "gcc: $GCC_VERSION"
   echo "FFLAGS: $FFLAGS"
   echo "CFLAGS: $CFLAGS"
   echo "LDFLAGS: $LDFLAGS"
   echo "DLM exported symbols: $exported"
   echo "sources:"
   /usr/bin/sha256sum "${FORTRAN_SOURCES[@]}" "$HERE/DISORTtoIDL_DP.c" "$HERE/exports.map" "$HERE/DISORTtoIDL_DP.dlm"
   echo "kernel objects:"
   /usr/bin/sha256sum "${objects[@]}"
   echo "products:"
   /usr/bin/sha256sum "$HERE/DISORTtoIDL_DP.so" "$HERE/libdisort_dp.so"
} > "$HERE/BUILD_INFO.txt"
sed "s|$HERE/||; s|$SRC/|disort2/src/|" -i "$HERE/BUILD_INFO.txt"
cat "$HERE/BUILD_INFO.txt"
