"""Extract the 90 captured V21 IDL single-precision DISORT calls into a compact reference file.

Source: the IDL validation run idl_double/results/run_20261010T214143,
case_{A,B,C,D}/validation/captures/call_*.sav (pre-call copies of every
argument and the outputs of the original DISORTtoIDL.so, captured inside the
V21 IDL generator; see idl_double/README and REPORT_v21_v22_single_reproduction.md).

Every value is kept exactly as IDL stored it (float32 / int); arrays are cut to
the layers actually passed (NLYR1) and put in Fortran dimension order:
PMOM (NMOM + 1, NLYR1), UU (NUMU, NTAU, NPHI).  TEMPER is NLYR1 + 1 values.

The IDL scattering-property file (scatfile.sav) and the compact V21 LUT of
each case are copied unchanged to reference/idl_compact/case_<X>/.

Usage: python build_v21_reference.py
Writes reference/v21_idl_single_calls.npz, reference/idl_compact/ and
reference/v21_idl_single_calls.json (provenance: sha256 of every source file).
"""

import hashlib
import json
from pathlib import Path
import shutil

import numpy as np
from scipy.io import readsav

HERE = Path(__file__).resolve().parent
RUN = HERE.parent / "idl_double" / "results" / "run_20261010T214143"
OUT_DIR = HERE / "reference"
CASES = ("A", "B", "C", "D")
SCALARS = ("fbeam", "umu0", "phi0", "fisot", "accur", "wvnmlo", "wvnmhi", "btemp", "ttemp", "temis", "albedo")
INTEGERS = ("nmom", "ntau", "numu", "nphi", "nstr_v21", "usrtau", "usrang", "ibcnd", "lamber", "plank", "onlyfl")
OUTPUTS = ("rfldir", "rfldn", "flup", "dfdt", "uavg", "uu", "albmed", "trnmed")


def main():
   OUT_DIR.mkdir(exist_ok=True)
   arrays = {}
   provenance = {"run": str(RUN.relative_to(HERE.parents[2])), "calls": []}
   index = 0
   for case in CASES:
      for path in sorted((RUN / f"case_{case}" / "validation" / "captures").glob("call_*.sav")):
         c = readsav(str(path))
         n = int(c["nlyr1"])
         key = f"c{index:03d}_"
         arrays[key + "case"] = np.array(case)
         arrays[key + "nlyr"] = np.int32(n)
         for name in INTEGERS:
            arrays[key + name] = np.int32(c[name])
         for name in SCALARS:
            arrays[key + name] = np.float32(c[name])
         arrays[key + "dtauc"] = np.asarray(c["dtauc1_in"], np.float32)[:n]
         arrays[key + "ssalb"] = np.asarray(c["ssalb1_in"], np.float32)[:n]
         arrays[key + "pmom"] = np.asfortranarray(np.atleast_2d(np.asarray(c["pmom1_in"], np.float32))[:n].T)
         temper = np.zeros(n + 1, np.float32)
         given = np.asarray(c["temper_in"], np.float32).ravel()
         temper[:min(n + 1, given.size)] = given[:n + 1]
         arrays[key + "temper"] = temper
         for name in ("utau", "umu", "phi"):
            arrays[key + name] = np.asarray(c[name + "_in"], np.float32)
         for name in OUTPUTS:
            value = np.asarray(c[name], np.float32)
            arrays[key + "out_" + name] = np.asfortranarray(value.T) if name == "uu" else value
         provenance["calls"].append({"index": index, "case": case, "capture": path.name,
                                     "sha256": hashlib.sha256(path.read_bytes()).hexdigest()})
         index += 1
   arrays["ncalls"] = np.int32(index)
   provenance["compact"] = []
   for case in CASES:
      source = RUN / f"case_{case}" / "validation" / "luts" / f"case_{case}"
      target = OUT_DIR / "idl_compact" / f"case_{case}"
      target.mkdir(parents=True, exist_ok=True)
      for name in ("scatfile.sav", "terra_modis_m_liquid-water_a01_p240_v21.nc"):
         shutil.copy2(source / name, target / name)
         provenance["compact"].append({"case": case, "file": name, "source": str((source / name).relative_to(HERE.parents[2])),
                                       "sha256": hashlib.sha256((source / name).read_bytes()).hexdigest()})
   np.savez_compressed(OUT_DIR / "v21_idl_single_calls.npz", **arrays)
   (OUT_DIR / "v21_idl_single_calls.json").write_text(json.dumps(provenance, indent=1) + "\n")
   print(f"{index} calls written to {OUT_DIR / 'v21_idl_single_calls.npz'}")


if __name__ == "__main__":
   main()
