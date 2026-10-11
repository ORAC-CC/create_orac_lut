#!/bin/bash
#
# run_idl_v21_lut.sh
#
# Worker for the double-precision V21 IDL LUTs (normally run by
# scripts/run_idl_v21_lut.slurm; see scripts/submit_v21_v22_double.sh).
#
#   run_idl_v21_lut.sh validate <out_dir>
#       Run the four compact Terra MODIS p240 validation cases
#       (validation/v21_v22_terra_modis/idl_double/cases) through the
#       production IDL path (IDL_PATH = IDL_DLM_PATH = +create_orac_lut,
#       unmodified create_orac_cloud_lut, call_disort.pro with precision
#       selection), once with ORAC_LUT_DISORT_PRECISION=double and once with
#       =single, and compare with the validated results of
#       idl_double/results/run_20261010T214143:
#         double -> the DISORT_DP generator LUTs, within 1 float32 ulp
#         single -> the historical compact V21 LUTs, bitwise
#       Exit status 0 only if all eight comparisons pass.
#
#   run_idl_v21_lut.sh product <pold|p240|pagg> <stage_dir>
#       Generate one complete Terra MODIS V21 product exactly as
#       create_orac_lut/terra_modis_run does (srf_quad=1, original .mm, LUT
#       grid, driver and T-matrix path, version=21) with double-precision
#       DISORT, in a private work directory; verify it against the historical
#       V21 product; write <product>.provenance.json; publish a read-only
#       copy into <stage_dir>.  A product already staged and verified is not
#       recomputed (safe restart).
#
# Writes only below the repository.  Never touches the shared LUT directory.

set -u -o pipefail

REPO="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
LEGACY="$REPO/create_orac_lut"
PYTHON=/home/g/grainger/miniforge3/envs/science/bin/python
REFERENCE="$REPO/validation/v21_v22_terra_modis/idl_double/results/run_20261010T214143"
CASES="$REPO/validation/v21_v22_terra_modis/idl_double/cases"
SHARED=/network/group/aopp/eodg/RGG004_GRAINGER_ORACFILE/ORAC_LUTS
TMATRIX=/network/group/aopp/eodg/shared/dubovik_tmatrix/
MANIFEST="$REPO/releases/idl_lut_v21_double.release"

fail() { echo "run_idl_v21_lut.sh: FAILED: $*" >&2; exit 1; }

setup_environment() {
   if ! type module > /dev/null 2>&1; then source /etc/profile.d/modules.sh; fi
   module load intel-compilers/2022 || fail "module load intel-compilers/2022"
   module load idl/890 || fail "module load idl/890"
   command -v idl > /dev/null || fail "idl not available on $(hostname)"
   export IDL_PATH="<IDL_DEFAULT>:+$LEGACY"
   export IDL_DLM_PATH="<IDL_DEFAULT>:+$LEGACY"
   unset CREATE_ORAC_CLOUD_LUT_DRIVER IDL_STARTUP
   for so in "$LEGACY/disort2_dp/DISORTtoIDL_DP.so" "$LEGACY/disort2/src/DISORTtoIDL.so" "$LEGACY/mie/dlm-code/mie_dlm_single.so"; do
      missing=$(ldd "$so" | grep 'not found')
      [ -z "$missing" ] || fail "unresolved libraries for $so: $missing"
   done
   "$PYTHON" "$REPO/scripts/check_release_manifest.py" "$MANIFEST" > /dev/null || fail "release files differ from $MANIFEST"
   echo "host $(hostname); idl $(command -v idl); modules: $(module list -t 2>&1 | tr '\n' ' ')"
}

stage_inputs() {   # $1 run directory; $2 optional compact case letter
   mkdir -p "$1/luts" "$1/input_files"
   for entry in "$LEGACY"/input_files/*; do
      name="$(basename "$entry")"
      if [ -n "${2:-}" ] && [ "$name" = lut ]; then continue; fi
      ln -s "$entry" "$1/input_files/$name"
   done
   if [ -n "${2:-}" ]; then
      mkdir -p "$1/input_files/lut"
      cp -p "$CASES/case_$2.lut" "$1/input_files/lut/"
      cp -p "$CASES/case_$2.driver" "$1/"
   fi
}

run_idl() {   # $1 directory, $2 IDL statement, $3 log
   # Resolve both paths before changing directory: a relative log path would
   # otherwise be opened relative to the run directory (job 494791).
   local dir log
   dir="$(realpath -m -- "$1")"
   log="$(realpath -m -- "$3")"
   set -- "$dir" "$2" "$log"
   ( cd "$1" && idl -e "$2" < /dev/null > "$3" 2>&1 )
   status=$?
   cat "$3"
   [ "$status" -eq 0 ] || fail "IDL exit status $status ($3)"
   grep -q -E '% Execution halted|^% Error|failed with code' "$3" && fail "IDL reported an error ($3)"
   grep -q 'ORAC LUT generation completed' "$3" || fail "generation did not complete ($3)"
   return 0
}

validate() {
   out="$(realpath -m -- "$1")"
   mkdir -p "$out" || fail "cannot create $out"
   setup_environment
   ok=1
   for c in A B C D; do
      for precision in double single; do
         dir="$out/case_${c}_$precision"
         rm -rf "$dir"
         stage_inputs "$dir" "$c"
         echo "=== case $c, ORAC_LUT_DISORT_PRECISION=$precision"
         ORAC_LUT_DISORT_PRECISION=$precision run_idl "$dir" "create_orac_cloud_lut_wrapper, driver='case_$c.driver'" "$dir/idl.log"
         grep -q "call_disort: DISORT precision = $precision" "$dir/idl.log" || fail "case $c: precision $precision not reported"
         new=$(ls "$dir"/luts/case_$c/*_v21.nc) || fail "case $c $precision: no LUT"
         if [ "$precision" = double ]; then
            grep -q 'Loaded DLM: DISORTTOIDL_DP' "$dir/idl.log" || fail "case $c: DISORTtoIDL_DP not loaded"
            reference=$(ls "$REFERENCE"/case_$c/dp_generator/luts/case_$c/*_v90.nc)
            tolerance=1
         else
            reference=$(ls "$REFERENCE"/case_$c/validation/luts/case_$c/*_v21.nc)
            tolerance=0
         fi
         "$PYTHON" "$REPO/scripts/compare_lut_numeric.py" "$new" "$reference" --max-ulp $tolerance \
            --record "$dir/comparison.json" || ok=0
      done
   done
   [ "$ok" -eq 1 ] || fail "integration validation: at least one comparison failed (see $out/*/comparison.json)"
   echo "V21 IDL INTEGRATION VALIDATION PASSED: double within 1 float32 ulp of the validated DISORT_DP results, single bitwise equal to historical V21"
}

product() {
   name="$1"
   stage="$(realpath -m -- "$2")"
   case "$name" in
      pold) mm=liquid-water_old.mm; lut=liquid-water-cloud.lut; file=terra_modis_m_liquid-water_a01_pold_v21.nc ;;
      p240) mm=liquid-water_240.mm; lut=liquid-water-cloud.lut; file=terra_modis_m_liquid-water_a01_p240_v21.nc ;;
      pagg) mm=water-ice_agg.mm; lut=ice-cloud.lut; file=terra_modis_m_water-ice_a01_pagg_v21.nc ;;
      *) fail "unknown product $name" ;;
   esac
   mkdir -p "$stage" || fail "cannot create $stage"
   if [ -f "$stage/$file" ] && [ -f "$stage/$file.verification.json" ] && grep -q '"passed": true' "$stage/$file.verification.json"; then
      echo "$stage/$file already staged and verified; nothing to do"
      return 0
   fi
   [ -e "$stage/$file" ] && fail "$stage/$file exists without a passing verification record; inspect it and remove it deliberately"
   setup_environment
   work="$stage/work/$name"
   rm -rf "$work"
   stage_inputs "$work"
   started=$(date --iso-8601=seconds)
   start_s=$(date +%s)
   echo "=== $name: $mm, double-precision DISORT, version 21; work $work"
   ORAC_LUT_DISORT_PRECISION=double run_idl "$work" \
      "create_orac_cloud_lut_wrapper, driver='input_files/driver/terra_modis_cloud.driver', srf_quad=1, mmfile='$mm', lutfile='$lut', tmatrix_path='$TMATRIX', version=21" \
      "$work/idl.log"
   grep -q 'call_disort: DISORT precision = double' "$work/idl.log" || fail "$name: double precision not reported"
   grep -q 'Loaded DLM: DISORTTOIDL_DP' "$work/idl.log" || fail "$name: DISORTtoIDL_DP not loaded"
   made="$work/luts/terra_modis_cloud/$file"
   [ -s "$made" ] || fail "$name: expected product $made not written"
   elapsed=$(( $(date +%s) - start_s ))
   "$PYTHON" - "$made" "$MANIFEST" "$work/idl.log" "$started" "$elapsed" "$REPO" <<'EOF' || fail "$name: provenance"
import hashlib, json, os, platform, subprocess, sys
made, manifest, log, started, elapsed, repo = sys.argv[1:]
def sha(p):
   h = hashlib.sha256()
   with open(p, "rb") as f:
      for b in iter(lambda: f.read(1 << 20), b""):
         h.update(b)
   return h.hexdigest()
git = lambda *a: subprocess.run(["/usr/bin/git", "-C", repo, *a], capture_output=True, text=True).stdout.strip()
record = {
   "provenance_schema": "oraclut-idl-lut-provenance-1",
   "product": os.path.basename(made), "product_sha256": sha(made),
   "lut_release": "V21", "generator": "legacy ORAC LUT V2.1 IDL generator (create_orac_lut)",
   "disort_precision": "double", "nstreams": 60, "nmom": 1000,
   "release_manifest": {"path": os.path.relpath(manifest, repo), "sha256": sha(manifest)},
   "git_commit": git("rev-parse", "HEAD"), "git_status_of_release_paths": git("status", "--porcelain", "--untracked-files=no", "--", "create_orac_lut", "releases", "scripts"),
   "started": started, "elapsed_seconds": int(elapsed), "host": platform.node(),
   "slurm_job_id": os.environ.get("SLURM_JOB_ID"), "idl_log_header": open(log).read().splitlines()[:4],
   "idl_command": "create_orac_cloud_lut_wrapper, driver='input_files/driver/terra_modis_cloud.driver', srf_quad=1, mmfile=..., lutfile=..., tmatrix_path=..., version=21",
   "environment": {k: os.environ.get(k) for k in ("ORAC_LUT_DISORT_PRECISION", "IDL_PATH", "IDL_DLM_PATH", "LOADEDMODULES")},
}
open(made[:-3] + ".provenance.json", "w").write(json.dumps(record, indent=1) + "\n")
EOF
   "$PYTHON" "$REPO/scripts/verify_staged_lut.py" "$made" "$SHARED/$file" --precision double \
      --record "$work/verification.json" || fail "$name: verification against the historical V21 product failed"
   cp -p "$made" "$stage/.$file.tmp" && mv "$stage/.$file.tmp" "$stage/$file" || fail "$name: publish"
   cp -p "${made%.nc}.provenance.json" "$stage/${file%.nc}.provenance.json"
   cp -p "$work/verification.json" "$stage/$file.verification.json"
   cp -p "$work/idl.log" "$stage/${file%.nc}.idl.log"
   chmod a-w "$stage/$file" "$stage/${file%.nc}.provenance.json" "$stage/$file.verification.json"
   echo "STAGED $stage/$file ($(/usr/bin/sha256sum "$stage/$file" | cut -d' ' -f1)) in ${elapsed} s"
}

case "${1:-}" in
   validate) [ $# -eq 2 ] || fail "usage: validate <out_dir>"; validate "$2" ;;
   product) [ $# -eq 3 ] || fail "usage: product <pold|p240|pagg> <stage_dir>"; product "$2" "$3" ;;
   *) sed -n '2,32p' "$0"; exit 2 ;;
esac
