#!/bin/bash
#
# build_dlm.sh
#
# Build the diagnostic double-precision DISORT DLM (DISORTtoIDL_DP) from the
# unmodified production DISORT sources, plus a plain shared library of the
# same objects (libdisort_dp_ifort.so, exporting disort_) for non-IDL tests.
#
# Precision: default REAL promoted to 8 bytes (-real-size 64); DOUBLE
# PRECISION stays 8 bytes (-double-size 64); INTEGER and LOGICAL stay 4
# bytes.  -fp-model precise forbids value-unsafe optimisation.  The
# historical DISORTtoIDL.so is not touched.
#
# Usage: bash build_dlm.sh        (writes only to ./build)

set -euo pipefail

HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
ROOT="$(cd "$HERE/../../.." && pwd)"
SRC="$ROOT/create_orac_lut/disort2/src"
LOCAL="$HERE/src"
BUILD="$HERE/build"
IDL_INCLUDE=/network/software/linux-x86_64/idl/890/idl89/external/include

if ! type module > /dev/null 2>&1; then source /etc/profile.d/modules.sh; fi
module load intel-compilers/2022
IFORT="$(command -v ifort)"
IFORT_DIR="$(dirname "$IFORT")"
# Minimal environment: the Kerberos-protected home must not be on PATH.
RUN="env -i PATH=$IFORT_DIR:/usr/bin:/bin HOME=/tmp LD_LIBRARY_PATH=${LD_LIBRARY_PATH:-} INTEL_LICENSE_FILE=${INTEL_LICENSE_FILE:-}"

FFLAGS="-O2 -fp-model precise -fpic -real-size 64 -double-size 64 -integer-size 32"
CFLAGS="-O2 -fPIC -std=gnu99 -Wall -Wextra -Wno-cast-function-type"
FORTRAN_SOURCES=("$SRC/DISORTfunctions.f" "$SRC/PRTFIN.f" "$SRC/disort/ErrPack.f" "$SRC/disort/LINPAK.f"
                 "$SRC/disort/RDI1MACH.f" "$LOCAL/dp_support.f")

rm -rf "$BUILD"
mkdir -p "$BUILD"
IFORT_VERSION="$($RUN "$IFORT" --version 2>&1 | head -1)"
GCC_VERSION="$(/usr/bin/gcc --version | head -1)"
INFO="ifort=[$IFORT_VERSION] FFLAGS=[$FFLAGS] gcc=[$GCC_VERSION] CFLAGS=[$CFLAGS]"

objects=()
for source in "${FORTRAN_SOURCES[@]}"; do
   object="$BUILD/$(basename "${source%.*}").o"
   $RUN "$IFORT" $FFLAGS -c "$source" -o "$object" 2> "$object.log" || { cat "$object.log"; exit 1; }
   objects+=("$object")
done
$RUN /usr/bin/gcc $CFLAGS -I"$IDL_INCLUDE" -DDP_BUILD_INFO="\"$INFO\"" -c "$LOCAL/DISORTtoIDL_DP.c" -o "$BUILD/DISORTtoIDL_DP.o"

$RUN "$IFORT" -shared -fpic -Wl,--version-script="$LOCAL/exports.map" -Wl,-Bsymbolic -Wl,--no-undefined-version \
   -o "$BUILD/DISORTtoIDL_DP.so" "$BUILD/DISORTtoIDL_DP.o" "${objects[@]}"
$RUN "$IFORT" -shared -fpic -o "$BUILD/libdisort_dp_ifort.so" "${objects[@]}"
cp "$LOCAL/DISORTtoIDL_DP.dlm" "$BUILD/DISORTtoIDL_DP.dlm"

# The DLM must export IDL_Load and nothing from the DISORT kernel.
exported="$(/usr/bin/nm -D --defined-only "$BUILD/DISORTtoIDL_DP.so" | awk '{print $3}' | sort | tr '\n' ' ')"
case " $exported " in
   *" disort_ "*|*" getmom_ "*|*" r1mach_ "*) echo "kernel symbols leaked from DLM: $exported" >&2; exit 1 ;;
esac
case " $exported " in *" IDL_Load "*) ;; *) echo "IDL_Load not exported" >&2; exit 1 ;; esac

{
   echo "date: $(date --iso-8601=seconds)"
   echo "host: $(hostname)"
   echo "ifort: $IFORT_VERSION ($IFORT)"
   echo "gcc: $GCC_VERSION"
   echo "FFLAGS: $FFLAGS"
   echo "CFLAGS: $CFLAGS"
   echo "DLM exported symbols: $exported"
   echo "git HEAD: $(git -C "$ROOT" rev-parse HEAD)"
   echo "sources:"
   sha256sum "${FORTRAN_SOURCES[@]}" "$LOCAL/DISORTtoIDL_DP.c" "$LOCAL/exports.map" "$LOCAL/DISORTtoIDL_DP.dlm"
   echo "products:"
   sha256sum "$BUILD/DISORTtoIDL_DP.so" "$BUILD/libdisort_dp_ifort.so" "$BUILD/DISORTtoIDL_DP.dlm"
   echo "dynamic dependencies of DISORTtoIDL_DP.so:"
   ldd "$BUILD/DISORTtoIDL_DP.so" | sed 's/ (0x[0-9a-f]*)//'
} > "$BUILD/BUILD_INFO.txt"
cat "$BUILD/BUILD_INFO.txt"
