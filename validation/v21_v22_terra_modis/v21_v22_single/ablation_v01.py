"""One-at-a-time ablation of the scattering_replica options (Stage 4).

For each option, start from the IDL arithmetic and switch that option alone
to the production-Python value (and conversely, start from the Python
arithmetic and switch that option alone to the IDL value); count the
scatfile elements that differ bitwise from the IDL scatfile.sav.

Usage: python ablation.py   -> results/scattering_ablation.json
"""

import json

import scattering_replica as sr

OPTIONS = ["ri", "wn", "literals", "quad_t", "quad_g", "w1p", "math", "mie", "store"]


def total(rows):
   return sum(v["differing"] for v in rows.values())


def main():
   result = {}
   for case in sr.CASES:
      ref = sr.idl_scatfile(case)
      result[case] = {}
      for option in OPTIONS:
         opt = dict(sr.IDL, **{option: sr.PYTHON[option]})
         idl_minus = sr.compare(sr.scattering(case, opt)[0], ref)
         opt = dict(sr.PYTHON, **{option: sr.IDL[option]})
         python_plus = sr.compare(sr.scattering(case, opt)[0], ref)
         result[case][option] = {"idl_with_python_option": idl_minus, "python_with_idl_option": python_plus}
         print(case, option, "IDL-with-python-option differing:", total(idl_minus), " python-with-IDL-option differing:", total(python_plus), flush=True)
      # Mie memory layout of the Intel kernel (IDL: static Dqv 64-aligned, IDL temporary F11 arrays)
      for dqv_offset, f_offset in ((0, 0), (8, 8), (8, 0)):
         rows = sr.compare(sr.scattering(case, dict(sr.IDL, mie_dqv_offset=dqv_offset, mie_f_offset=f_offset))[0], ref)
         result[case][f"mie_layout_{dqv_offset}_{f_offset}"] = rows
         print(case, "Intel Mie layout dqv", dqv_offset, "f11", f_offset, "differing:", total(rows), flush=True)
   (sr.HERE / "results" / "scattering_ablation.json").write_text(json.dumps(result, indent=1) + "\n")


if __name__ == "__main__":
   main()
