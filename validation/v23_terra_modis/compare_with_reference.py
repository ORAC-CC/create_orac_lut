"""Validate the V23 Terra MODIS products and compare them with V22 at exactly shared grid nodes.

V23 (releases/lut_v23.release) is the unchanged lut-code-v22.4 scientific code
with the Grid B definitions.  For each product this checks:

* structure: the V22 variable set, dimensions equal to the Grid B definition
  (create_orac_lut/input_files/lut/<grid>_v02.lut through the production
  load_lutstr), coordinates equal to those definitions, 36 channels, no
  NaN/Inf/fill in any floating variable;
* provenance: lut_release V23, the V23 manifest, the sealed V22 driver and the
  CODE_VERSION code release lut-code-v22.4;
* V23 against the definitive V22 product at the nodes whose coordinates are
  exactly shared by the two grids (no interpolation): bitwise counts and
  maximum differences per RT operator.

Quantities that do not depend on the user viewing angles (fluxes: R_dd, T_dd,
R_0d, T_0d, T_00) are expected to agree bitwise.  Intensities at user angles
(R_0v, R_dv, T_dv, E_md) can differ: the original DISORT kernel's intensity
arithmetic depends on the number of user polar angles (README.md).

Usage: python validation/v23_terra_modis/compare_with_reference.py [candidate_dir]
   (default staging/lut_v23)   -> validation/v23_terra_modis/results/<product>.json
"""

import json
from pathlib import Path
import sys

import numpy as np
from netCDF4 import Dataset

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
sys.path.insert(0, str(ROOT / "src"))
from oraclut.idl_mirror import load_inststr, load_lutstr   # noqa: E402

SHARED = Path("/network/group/aopp/eodg/RGG004_GRAINGER_ORACFILE/ORAC_LUTS")
PRODUCTS = {"pold": ("terra_modis_m_liquid-water_a01_pold", "liquid-water-cloud_v02.lut"),
            "p240": ("terra_modis_m_liquid-water_a01_p240", "liquid-water-cloud_v02.lut"),
            "pagg": ("terra_modis_m_water-ice_a01_pagg", "ice-cloud_v02.lut")}
RT = ("R_0v", "R_0d", "T_0d", "T_00", "R_dv", "T_dv", "R_dd", "T_dd", "E_md")
COORDS = {"optical_depth": "opd", "effective_radius": "efr", "solar_zenith": "soz", "satellite_zenith": "saz",
          "relative_azimuth": "raa"}
FILL = np.float32(9.96921e36)


CHANNEL_IDS = {"channels": "channel_id", "solar_channels": "solar_channel_id", "mixed_channels": "mixed_channel_id",
               "thermal_channels": "thermal_channel_id"}


def coordinate(dataset, dim):
   """Coordinate values of a dimension (instrument channel numbers for the channel dimensions)."""
   if dim in CHANNEL_IDS and CHANNEL_IDS[dim] in dataset.variables:
      return np.asarray(dataset[CHANNEL_IDS[dim]][:], np.float64)
   return np.asarray(dataset[dim][:], np.float64) if dim in dataset.variables else np.arange(len(dataset.dimensions[dim]))


def structure(path, grid):
   inputs = ROOT / "create_orac_lut" / "input_files"
   inststr = load_inststr(inputs / "inst" / "terra_modis_v1.inst", requestedchannelid=list(range(1, 37)))
   lutstr = load_lutstr(inputs / "lut" / grid, inststr.max_sat_zenith)
   result = {}
   with Dataset(path) as d:
      d.set_auto_maskandscale(False)
      for dim, field in COORDS.items():
         expected = np.asarray(getattr(lutstr, field), np.float32)
         values = np.asarray(d[dim][:])
         result[dim] = {"n": int(values.size), "expected_n": int(expected.size),
                        "bitwise": bool(values.size == expected.size and np.array_equal(values.view(np.uint32), expected.view(np.uint32))),
                        "equal": bool(values.size == expected.size and np.allclose(values, expected, rtol=1e-6, atol=1e-12))}
      result["channels"] = int(len(d.dimensions["channels"]))
      nonfinite, fill = {}, {}
      for name, var in d.variables.items():
         if np.issubdtype(var.dtype, np.floating):
            a = np.asarray(var[:])
            if np.count_nonzero(~np.isfinite(a)):
               nonfinite[name] = int(np.count_nonzero(~np.isfinite(a)))
            if np.count_nonzero(a == FILL):
               fill[name] = int(np.count_nonzero(a == FILL))
      result["nan_inf"] = nonfinite                      # must be empty
      result["fill"] = fill                              # schema fill (e.g. old* for non-solar channels): compare with V22
      result["variables"] = list(d.variables)
      provenance = json.loads(d.getncattr("provenance_json"))
      header = provenance["release_manifest"]["header"]
      result["provenance"] = {"lut_release": d.getncattr("lut_release"),
                              "manifest": provenance["release_manifest"]["path"],
                              "manifest_sha256": provenance["release_manifest"]["sha256"],
                              "scientific_code": header.get("scientific_code"), "grid": header.get("grid"),
                              "driver": provenance["driver"]["path"], "git_commit": provenance["git_at_start"].get("commit")}
   return result


def shared_nodes(v22_path, v23_path):
   rows = {}
   with Dataset(v22_path) as r, Dataset(v23_path) as c:
      r.set_auto_maskandscale(False)
      c.set_auto_maskandscale(False)
      for name in RT:
         if name not in r.variables:
            continue
         a, b = r[name], c[name]
         ia, ib = [], []
         for dim in a.dimensions:
            x, y = coordinate(r, dim), coordinate(c, dim)
            common = [v for v in x if np.any(np.isclose(y, v, rtol=1e-6, atol=1e-9))]
            ia.append([int(np.argmin(abs(x - v))) for v in common])
            ib.append([int(np.argmin(abs(y - v))) for v in common])
         va = np.asarray(a[:])[np.ix_(*ia)]
         vb = np.asarray(b[:])[np.ix_(*ib)]
         diff = va.view(np.uint32) != vb.view(np.uint32)
         absd = np.abs(va.astype(np.float64) - vb)
         rows[name] = {"shared_elements": int(va.size), "bitwise": int(va.size - diff.sum()), "differing": int(diff.sum()),
                       "max_abs": float(absd.max()) if va.size else 0.0,
                       "p99_abs": float(np.percentile(absd, 99)) if va.size else 0.0,
                       "shared_axes": {dim: len(i) for dim, i in zip(a.dimensions, ia)}}
   return rows


def main():
   candidates = Path(sys.argv[1]) if len(sys.argv) > 1 else ROOT / "staging" / "lut_v23"
   out = HERE / "results"
   out.mkdir(exist_ok=True)
   for product, (stem, grid) in PRODUCTS.items():
      v23 = candidates / f"{stem}_v23.nc"
      if not v23.is_file():
         print(product, "not staged")
         continue
      with Dataset(SHARED / f"{stem}_v22.nc") as r, Dataset(v23) as c:
         same_variables = list(r.variables) == list(c.variables)
         r.set_auto_maskandscale(False)
         v22_fill = {n: int(np.count_nonzero(np.asarray(r[n][:]) == FILL)) for n in r.variables if np.issubdtype(r[n].dtype, np.floating)}
         v22_fill = {n: c for n, c in v22_fill.items() if c}
      record = {"product": str(v23), "v22": str(SHARED / f"{stem}_v22.nc"), "same_variable_set_as_v22": same_variables,
                "structure": structure(v23, grid), "v22_fill": v22_fill, "shared_nodes": shared_nodes(SHARED / f"{stem}_v22.nc", v23)}
      record["fill_as_v22"] = record["structure"]["fill"] == v22_fill
      (out / f"{product}.json").write_text(json.dumps(record, indent=1, default=str) + "\n")
      s = record["structure"]
      print(product, "variables as V22:", same_variables, "| coords bitwise:", {k: v["bitwise"] for k, v in s.items() if isinstance(v, dict) and "bitwise" in v},
            "| channels", s["channels"], "| NaN/Inf", s["nan_inf"] or "none", "| fill as V22", record["fill_as_v22"], "| release", s["provenance"]["lut_release"], s["provenance"]["scientific_code"], s["provenance"]["driver"])
      for name, row in record["shared_nodes"].items():
         print(f"   {name:5s} shared {row['shared_elements']:7d} bitwise {row['bitwise']:7d} max|d| {row['max_abs']:.3g}")


if __name__ == "__main__":
   main()
