"""Verify a newly generated ORAC LUT against the historical product it replaces.

   python scripts/verify_staged_lut.py NEW.nc HISTORICAL.nc [--precision double] [--record OUT.json]

Checks, failing (exit 1) on any problem:
  * NEW opens as NetCDF and has the same dimensions (names and lengths);
  * the same variables, with the same data types and dimension order;
  * identical coordinate values (optical_depth, effective_radius, solar and
    satellite zenith, relative azimuth, channel identifiers);
  * every RT operator and single-scattering property is finite and not
    identically zero;
  * every text attribute is NC_STRING (the ORAC retrieval reads attributes
    with nc_get_att_string; LUT_FORMAT_CONTRACT.md);
  * a provenance sidecar NEW.provenance.json exists, names the expected
    DISORT precision and records the SHA-256 of NEW.
Writes a JSON verification record when --record is given.  Read-only.
"""

import argparse
import hashlib
import json
import sys
from pathlib import Path

import numpy as np
from netCDF4 import Dataset

COORDINATES = ("optical_depth", "effective_radius", "solar_zenith", "satellite_zenith", "relative_azimuth",
               "channel_id", "solar_channel_id", "thermal_channel_id", "mixed_channel_id")
SCIENCE = ("R_0v", "R_0d", "T_0d", "T_00", "R_dv", "T_dv", "R_dd", "T_dd", "E_md", "extinction_coefficient",
           "extinction_coefficient_ratio", "single_scatter_albedo", "asymmetry_parameter")


def sha256(path):
   digest = hashlib.sha256()
   with open(path, "rb") as handle:
      for block in iter(lambda: handle.read(1 << 20), b""):
         digest.update(block)
   return digest.hexdigest()


def verify(new, old, precision):
   problems, record = [], {"product": str(new), "historical": str(old), "sha256": sha256(new)}
   with Dataset(new) as a, Dataset(old) as b:
      dims_a = {k: len(v) for k, v in a.dimensions.items()}
      dims_b = {k: len(v) for k, v in b.dimensions.items()}
      record["dimensions"] = dims_a
      if dims_a != dims_b:
         problems.append(f"dimensions differ: {dims_a} vs {dims_b}")
      vars_a, vars_b = set(a.variables), set(b.variables)
      if vars_a != vars_b:
         problems.append(f"variables differ: only new {sorted(vars_a - vars_b)}, only historical {sorted(vars_b - vars_a)}")
      for name in sorted(vars_a & vars_b):
         if a[name].dtype != b[name].dtype or a[name].dimensions != b[name].dimensions:
            problems.append(f"{name}: dtype/dimensions {a[name].dtype}{a[name].dimensions} vs {b[name].dtype}{b[name].dimensions}")
      for name in COORDINATES:
         if name in a.variables and name in b.variables:
            if not np.array_equal(np.ma.getdata(a[name][:]), np.ma.getdata(b[name][:])):
               problems.append(f"coordinate {name} differs")
      finite = {}
      for name in SCIENCE:
         if name not in a.variables:
            continue
         # Raw values: netCDF4's automatic masking would also hide legitimate
         # values outside valid_range (e.g. slightly negative transmittances).
         variable = a[name]
         variable.set_auto_maskandscale(False)
         values = np.asarray(variable[:], dtype=np.float64)
         fill = getattr(variable, "_FillValue", None)
         if fill is not None:
            values = np.where(values == float(fill), np.nan, values)
         finite[name] = {"finite_fraction": float(np.mean(np.isfinite(values))),
                         "min": float(np.nanmin(values)), "max": float(np.nanmax(values))}
         if not np.all(np.isfinite(values)):
            problems.append(f"{name}: {int(np.sum(~np.isfinite(values)))} non-finite or fill values")
         if np.nanmax(np.abs(values)) == 0.0:
            problems.append(f"{name}: identically zero")
      record["science_variables"] = finite
      record["global_attributes"] = {k: str(a.getncattr(k))[:200] for k in a.ncattrs()}
   string_check = Path(__file__).resolve().parents[1] / "src" / "oraclut" / "repair_string_attributes.py"
   if string_check.is_file():
      import subprocess
      completed = subprocess.run([sys.executable, "-m", "oraclut.repair_string_attributes", "--check", str(new)],
                                 capture_output=True, text=True, cwd=string_check.parents[2],
                                 env={"PYTHONPATH": str(string_check.parents[1]), "PATH": "/usr/bin:/bin"})
      record["string_attribute_check"] = (completed.stdout + completed.stderr).strip()[-500:]
      if completed.returncode != 0:
         problems.append("text attributes are not all NC_STRING (repair_string_attributes --check failed)")
   sidecar = Path(str(new)[:-3] + ".provenance.json")
   if not sidecar.is_file():
      problems.append(f"provenance sidecar missing: {sidecar}")
   else:
      provenance = json.loads(sidecar.read_text())
      text = json.dumps(provenance)
      if precision and f'"disort_precision": "{precision}"' not in text:
         problems.append(f"provenance does not record disort_precision = {precision}")
      if record["sha256"] not in text and provenance.get("product_sha256") not in (None, record["sha256"]):
         problems.append("provenance records a different product SHA-256")
   record["problems"] = problems
   record["passed"] = not problems
   return record


def main():
   parser = argparse.ArgumentParser(description=__doc__.split("\n")[0])
   parser.add_argument("new", type=Path)
   parser.add_argument("historical", type=Path)
   parser.add_argument("--precision")
   parser.add_argument("--record", type=Path)
   args = parser.parse_args()
   record = verify(args.new, args.historical, args.precision)
   if args.record:
      args.record.write_text(json.dumps(record, indent=1) + "\n")
   for problem in record["problems"]:
      print("PROBLEM", problem)
   print(("VERIFIED " if record["passed"] else "FAILED ") + str(args.new))
   return 0 if record["passed"] else 1


if __name__ == "__main__":
   raise SystemExit(main())
