#!/bin/bash
#
# build_reference.sh
#
# Independent-compiler references for the double-precision DLM: the same
# unmodified production DISORT sources built with gfortran
#   double: -fdefault-real-8 -fdefault-double-8      (REAL = DOUBLE = 8 bytes)
#   quad:   -freal-4-real-16 -freal-8-real-16        (REAL = DOUBLE = 16 bytes)
# each linked to disort_ref_driver.f90.  Writes only to ./build.
#
# Usage: bash build_reference.sh

set -euo pipefail

HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
ROOT="$(cd "$HERE/../../.." && pwd)"
SRC="$ROOT/create_orac_lut/disort2/src"
LOCAL="$HERE/src"
BUILD="$HERE/build"
RUN="env -i PATH=/usr/bin:/bin HOME=/tmp"
BASE="-O2 -std=legacy -ffixed-line-length-0"
SOURCES=("$SRC/DISORTfunctions.f" "$SRC/PRTFIN.f" "$SRC/disort/ErrPack.f" "$SRC/disort/LINPAK.f"
         "$SRC/disort/RDI1MACH.f" "$LOCAL/dp_support.f")
mkdir -p "$BUILD"

for variant in double quad; do
   case $variant in
      double) PRECISION="-fdefault-real-8 -fdefault-double-8" ;;
      quad)   PRECISION="-freal-4-real-16 -freal-8-real-16" ;;
   esac
   out="$BUILD/ref_$variant"
   rm -rf "$out"
   mkdir -p "$out"
   objects=()
   for source in "${SOURCES[@]}"; do
      object="$out/$(basename "${source%.*}").o"
      $RUN /usr/bin/gfortran $BASE $PRECISION -c "$source" -o "$object" 2> "$object.log" || { cat "$object.log"; exit 1; }
      objects+=("$object")
   done
   $RUN /usr/bin/gfortran -O2 $PRECISION -c "$LOCAL/disort_ref_driver.f90" -o "$out/driver.o" -J "$out"
   $RUN /usr/bin/gfortran -o "$out/disort_ref_driver" "$out/driver.o" "${objects[@]}"
   {
      echo "variant: $variant"
      echo "gfortran: $(/usr/bin/gfortran --version | head -1)"
      echo "flags: $BASE $PRECISION"
      sha256sum "${SOURCES[@]}" "$LOCAL/disort_ref_driver.f90" "$out/disort_ref_driver"
   } > "$out/BUILD_INFO.txt"
   cat "$out/BUILD_INFO.txt"
done
