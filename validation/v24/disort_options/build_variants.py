"""Build experimental copies of the production DISORT kernel (validation only).

The production kernel is create_orac_lut/disort2/src/DISORTfunctions.f
(DISORT 2.0 beta, localised for ORAC) behind the C adapter
src/oraclut/radiative_transfer/legacy_disort.c, compiled by
oraclut.radiative_transfer.legacy_disort into build/oraclut/.  Nothing there
is touched.  This script copies the sources into
validation/tmp/disort_options/build/<variant>/ and changes only:

   * PARAMETER MXCMU = 100 -> 256 (so that NSTR up to 256 can be run; array
     bounds only; variants *100 keep 100, for timing comparable with
     production, because DISORT zeroes its MXCMU-sized arrays on every call);
   * DELTAM = .TRUE. and CORINT = .TRUE. (hard-wired in DISORT 2.0) -> values
     taken from a COMMON block that the adapter sets (default .TRUE.);
   * the adapter's fixed ACCUR = 1e-8 -> a value the caller sets (default 1e-8);

and, for variant "double", compiles every Fortran unit with
-fdefault-real-8 -fdefault-double-8 behind an adapter whose REAL arguments
are double.  Compiler flags otherwise follow legacy_disort._build_library
(gfortran -O2 -fPIC -std=legacy -ffixed-line-length-0; gcc -O2 -fPIC).

   env -i PATH=/usr/bin:/bin HOME=/tmp python validation/v24/disort_options/build_variants.py

(the minimal environment avoids executing anything from the Kerberos-
protected home directory).
"""

import re
import subprocess
from pathlib import Path

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[2]
SOURCE = ROOT / "create_orac_lut" / "disort2" / "src"
BUILD = ROOT / "validation" / "tmp" / "disort_options" / "build"
MXCMU = 256

ADAPTER = r"""/* Validation copy of src/oraclut/radiative_transfer/legacy_disort.c with ACCUR,
 * DELTAM and CORINT settable (defaults = production).  REAL_T is float, or
 * double for the double-precision build. */
#include <stdint.h>
typedef REAL_T real_t;
extern int32_t disort_();
extern struct { int32_t deltam; int32_t corint; } oraexp_;
static real_t exp_accur = (real_t)1.0e-8;

void oraclut_set_options(double accur, int32_t deltam, int32_t corint) {
    exp_accur = (real_t)accur;
    oraexp_.deltam = deltam;
    oraexp_.corint = corint;
}

int oraclut_disort(int32_t nlyr, real_t *dtauc, real_t *ssalb, int32_t nmom,
                   real_t *pmom, real_t *temper, real_t wvnlo, real_t wvnhi,
                   int32_t plank, int32_t ntau, real_t *utau, int32_t nstr,
                   int32_t numu, real_t *umu, int32_t nphi, real_t *phi,
                   real_t fbeam, real_t umu0, real_t fisot,
                   real_t *rfldir, real_t *rfldn, real_t *flup, real_t *dfdt,
                   real_t *uavg, real_t *uu, real_t *albmed, real_t *trnmed) {
    int32_t usrtau = 1, usrang = 1, ibcnd = 0, lamber = 1, onlyfl = 0;
    real_t albedo = 0, btemp = 0, ttemp = 0, temis = 0, phi0 = 0;
    real_t accur = exp_accur;
    int32_t prnt[5] = {0, 0, 0, 0, 0};
    int32_t maxcly = nlyr, maxulv = ntau, maxumu = numu, maxphi = nphi, maxmom = nmom;
    return (int)disort_(&nlyr, dtauc, ssalb, &nmom, pmom, temper,
                        &wvnlo, &wvnhi, &usrtau, &ntau, utau, &nstr,
                        &usrang, &numu, umu, &nphi, phi, &ibcnd, &fbeam,
                        &umu0, &phi0, &fisot, &lamber, &albedo, &btemp,
                        &ttemp, &temis, &plank, &onlyfl, &accur, prnt,
                        &maxcly, &maxulv, &maxumu, &maxphi, &maxmom,
                        rfldir, rfldn, flup, dfdt, uavg, uu, albmed,
                        trnmed);
}
"""


def patched_disort(mxcmu):
   text = (SOURCE / "DISORTfunctions.f").read_text()
   old_parameter = "      PARAMETER ( MXCLY = 46, MXULV = 2, MXCMU = 100, MXUMU = 180,"
   assert text.count(old_parameter) == 1
   text = text.replace(old_parameter, f"      PARAMETER ( MXCLY = 46, MXULV = 2, MXCMU = {mxcmu}, MXUMU = 180,")
   save = "      SAVE      DITHER, PASS1, PI, RPD, SQT\n"
   assert text.count(save) == 1
   text = text.replace(save, "      LOGICAL   EXPDLT, EXPCOR\n      COMMON / ORAEXP / EXPDLT, EXPCOR\n" + save)
   # the main-program assignments (6-space indent); SLFTST's own (deeper
   # indented) self-test settings are left alone
   for flag, source in (("DELTAM", "EXPDLT"), ("CORINT", "EXPCOR")):
      pattern = re.compile(rf"^      {flag} = \.TRUE\.$", re.MULTILINE)
      assert len(pattern.findall(text)) == 1, flag
      text = pattern.sub(f"      {flag} = {source}", text)
   return text


def build(variant):
   """variant: single | double (MXCMU 256) or single100 | double100 (production MXCMU = 100, for timing)."""

   out = BUILD / variant
   mxcmu = 100 if variant.endswith("100") else MXCMU
   out.mkdir(parents=True, exist_ok=True)
   (out / "DISORTfunctions_exp.f").write_text(patched_disort(mxcmu))
   real = "double" if variant.startswith("double") else "float"
   (out / "adapter.c").write_text(ADAPTER.replace("REAL_T", real))
   fflags = ["-O2", "-fPIC", "-std=legacy", "-ffixed-line-length-0"]
   if variant.startswith("double"):
      fflags += ["-fdefault-real-8", "-fdefault-double-8"]
   objects = []
   for source in (out / "DISORTfunctions_exp.f", SOURCE / "BDREF.f", SOURCE / "PRTFIN.f", SOURCE / "disort" / "ErrPack.f",
                  SOURCE / "disort" / "LINPAK.f", SOURCE / "disort" / "RDI1MACH.f"):
      obj = out / (source.stem.lower() + ".o")
      subprocess.run(["gfortran", *fflags, "-c", str(source), "-o", str(obj)], check=True)
      objects.append(obj)
   subprocess.run(["gcc", "-O2", "-fPIC", "-std=c11", "-c", str(out / "adapter.c"), "-o", str(out / "adapter.o")], check=True)
   for name, source in (("ambrals_fortran.o", "ambrals-fortran.c"), ("ambrals_c.o", "ambralsfor.c")):
      subprocess.run(["gcc", "-O2", "-fPIC", "-D_GNU_SOURCE", "-std=c11", "-c", str(SOURCE / source), "-o", str(out / name)],
                     check=True)
   library = out / f"libdisort_{variant}.so"
   subprocess.run(["gfortran", "-shared", "-o", str(library), str(out / "adapter.o"), str(out / "ambrals_fortran.o"),
                   str(out / "ambrals_c.o"), *(str(o) for o in objects)], check=True)
   print("built", library)


if __name__ == "__main__":
   for variant in ("single", "double", "single100", "double100"):
      build(variant)
