"""Write the source/binary/input manifest of an IDL ORAC LUT release (record of how it was made).

   python scripts/make_idl_release_manifest.py single releases/idl_lut_v21_single.release
   python scripts/make_idl_release_manifest.py double releases/idl_lut_v21_double.release

An IDL release is run by the legacy generator under create_orac_lut/ with
IDL_PATH = +create_orac_lut and IDL_DLM_PATH = +create_orac_lut
(create_orac_lut/initpath.bash).  Because the search path is recursive,
every .pro file in the tree is part of the source set.  The manifest lists,
with SHA-256:

  idl-source   every .pro file IDL can resolve (for the historical single-
               precision release, call_disort.pro is its preserved original
               call_disort_v01.pro)
  run          the V21 run script and IDL path set-up
  dlm          DLM descriptors and the compiled DLM libraries IDL loads
  dlm-source   the sources those DLMs were built from
  input        the exact input files the three Terra MODIS products read
  runtime      the Intel Fortran runtime libraries the DLMs link
  lut          (single only) the historical products in the shared directory

The format follows releases/lut_v22.release (oraclut-release-manifest-1
columns) under its own schema name, oraclut-idl-release-manifest-1, so the
Python release launcher (which scans releases/lut_v*.release) never reads it.
The script refuses to overwrite an existing manifest: a manifest is sealed
once written.  Check with  python scripts/check_release_manifest.py <manifest>.
"""

import datetime
import hashlib
import subprocess
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
LEGACY = ROOT / "create_orac_lut"
INTEL_LIB = Path("/network/software/linux-x86_64/compilers/intel/2022/compiler/2022.1.0/linux/compiler/lib/intel64_lin")
SHARED_LUTS = Path("/network/group/aopp/eodg/RGG004_GRAINGER_ORACFILE/ORAC_LUTS")
HISTORICAL_LUTS = ("terra_modis_m_liquid-water_a01_pold_v21.nc", "terra_modis_m_liquid-water_a01_p240_v21.nc",
                   "terra_modis_m_water-ice_a01_pagg_v21.nc")

INPUTS = (["input_files/inst/terra_modis_v1.inst", "input_files/sun/Gueymard2018.sssi", "input_files/atm/mls.atm",
           "input_files/microphysics/liquid-water_old.mm", "input_files/microphysics/liquid-water_240.mm",
           "input_files/microphysics/water-ice_agg.mm", "input_files/ri/oldwater.ri",
           "input_files/ri/H2O_240K_Rowe_2020.ri",
           "input_files/baum/AggregateSolidColumns_SeverelyRough_AllWavelengths_FullPhaseMatrix.nc",
           "input_files/lut/liquid-water-cloud.lut", "input_files/lut/ice-cloud.lut",
           "input_files/driver/terra_modis_cloud.driver"]
          + [f"input_files/srf/rtcoef_eos_1_modis-shifted_srf_ch{c:02d}.txt" for c in range(1, 37)])

DLM_SINGLE = ["disort2/src/DISORTtoIDL.dlm", "disort2/src/DISORTtoIDL.so", "mie/dlm-code/mie_dlm_single.dlm",
              "mie/dlm-code/mie_dlm_single.so"]
DLM_SOURCE_SINGLE = ["disort2/src/" + n for n in ("DISORTfunctions.f", "BDREF.f", "PRTFIN.f", "DISORTwrapper.c",
                                                  "DISORTtoIDL.c", "DISORTtoIDL.h", "DISORTfunctions.h",
                                                  "ambrals-fortran.c", "ambralsfor.c", "ambralsfor.h", "Makefile",
                                                  "make.inc", "disort/ErrPack.f", "disort/LINPAK.f",
                                                  "disort/RDI1MACH.f")]
DLM_DOUBLE = ["disort2_dp/DISORTtoIDL_DP.dlm", "disort2_dp/DISORTtoIDL_DP.so", "disort2_dp/libdisort_dp.so"]
DLM_SOURCE_DOUBLE = ["disort2_dp/" + n for n in ("DISORTtoIDL_DP.c", "dp_support.f", "exports.map", "build_dlm.sh",
                                                 "BUILD_INFO.txt")]
RUNTIME = ["libifcoremt.so.5", "libifport.so.5", "libimf.so", "libsvml.so", "libintlc.so.5", "libirc.so",
           "libifcore.so.5"]


def sha256(path):
   digest = hashlib.sha256()
   with open(path, "rb") as handle:
      for block in iter(lambda: handle.read(1 << 20), b""):
         digest.update(block)
   return digest.hexdigest()


def entry(role, logical, version, path):
   path = Path(path)
   shown = str(path.relative_to(ROOT)) if path.is_relative_to(ROOT) else str(path)
   return f"{role:<11s} {logical:<44s} {version:<7s} {sha256(path)}  {shown}"


def main(precision, output):
   output = ROOT / output if not Path(output).is_absolute() else Path(output)
   if output.exists():
      sys.exit(f"refusing to overwrite sealed manifest {output}")
   if precision not in ("single", "double"):
      sys.exit("precision must be single or double")
   lines = []
   pro_files = sorted(p for p in LEGACY.rglob("*.pro")
                      if "luts" not in p.relative_to(LEGACY).parts and not p.stem[-4:-2] == "_v")
   for path in pro_files:
      logical = str(path.relative_to(LEGACY))
      if logical == "call_disort.pro":
         if precision == "single":
            lines.append(entry("idl-source", logical, "v01", LEGACY / "call_disort_v01.pro"))
         else:
            versions = sorted(LEGACY.glob("call_disort_v[0-9][0-9].pro"))
            latest = versions[-1]
            if sha256(latest) != sha256(path):
               sys.exit(f"call_disort.pro differs from its latest version {latest.name}: version it first")
            lines.append(entry("idl-source", logical, latest.stem[-3:], latest))
      else:
         lines.append(entry("idl-source", logical, "pinned", path))
   for name in ("terra_modis_run", "initpath.bash"):
      lines.append(entry("run", name, "pinned", LEGACY / name))
   for name in DLM_SINGLE + (DLM_DOUBLE if precision == "double" else []):
      lines.append(entry("dlm", name, "pinned", LEGACY / name))
   for name in DLM_SOURCE_SINGLE + (DLM_SOURCE_DOUBLE if precision == "double" else []):
      lines.append(entry("dlm-source", name, "pinned", LEGACY / name))
   for name in INPUTS:
      lines.append(entry("input", name, "pinned", LEGACY / name))
   for name in RUNTIME:
      lines.append(entry("runtime", name, "pinned", INTEL_LIB / name))
   if precision == "single":
      for name in HISTORICAL_LUTS:
         lines.append(entry("lut", name, "product", SHARED_LUTS / name))

   idl_version = "IDL 8.9.0 (linux x86_64 m64), module idl/890"
   ifort = "intel-compilers/2022 (ifort 2021.6.0 runtime)"
   if precision == "single":
      description = """\
# IDL LUT release V21, historical single precision: the legacy ORAC LUT V2.1
# IDL generator as it produced the Terra MODIS V21 products of 2026-10-09/10
# (pold, p240, pagg) with the historical single-precision DISORT DLM
# (create_orac_lut/disort2/src/DISORTtoIDL.so, preserved byte for byte as
# DISORTtoIDL_v01.so).  call_disort.pro is listed as its preserved original,
# call_disort_v01.pro.  The current call_disort.pro reproduces it exactly when
# ORAC_LUT_DISORT_PRECISION=single (README in create_orac_lut/disort2_dp/).
#
# Reproduction evidence: compact recomputations of four Terra MODIS p240
# nodes with these sources equal the production V21 products bitwise in all
# nine RT coefficients (validation/v21_v22_terra_modis/idl_double/results/
# run_20261010T214143; REPORT_idl_double_precision.md).
#
# The lut entries are the historical products in the shared ORAC LUT
# directory; they are archived there and never overwritten."""
   else:
      description = """\
# IDL LUT release V21, double precision (new reference candidate): the same
# legacy ORAC LUT V2.1 IDL generator and inputs as the historical V21 release,
# with DISORT solved by the validated double-precision DLM
# create_orac_lut/disort2_dp/DISORTtoIDL_DP.so (default of call_disort.pro;
# ORAC_LUT_DISORT_PRECISION=double).  NSTR 60, NMOM 1000, original LUT grids,
# float32 NetCDF storage, LUT version 21.  Validation:
# validation/v21_v22_terra_modis/REPORT_idl_double_precision.md.
# Products made from this manifest are recorded, with SHA-256, in the
# companion file of the same name with suffix .products."""
   header = f"""\
# ORAC IDL LUT release manifest
#
{description}
#
# LUT release numbers (V21, V22, ...) and source-file versions (_v01, _v02, ...)
# are independent: each file's version advances only when that file changes.
# "pinned" files have no versioned name; their content is fixed by SHA-256.
#
# Columns: role  logical-name  file-version  sha256  path (relative to the
# repository, or absolute outside it)

manifest_schema = oraclut-idl-release-manifest-1
lut_release = 21
disort_precision = {precision}
nstreams = 60
nmom = 1000
idl = {idl_version}
fortran_runtime = {ifort}
created = {datetime.date.today().isoformat()}
git_head_at_creation = {subprocess.run(['/usr/bin/git', '-C', str(ROOT), 'rev-parse', 'HEAD'], capture_output=True, text=True).stdout.strip()}

"""
   output.write_text(header + "\n".join(lines) + "\n")
   print(f"wrote {output} ({len(lines)} entries)")


if __name__ == "__main__":
   if len(sys.argv) != 3:
      sys.exit(__doc__)
   main(sys.argv[1], sys.argv[2])
