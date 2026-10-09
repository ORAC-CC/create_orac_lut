"""Cost of the adaptive Legendre expansion (Phase 3 benchmark).

Validation-only code.  Times generate_scattering_properties (single-thread
CPU) for liquid-water stg on the V23 Grid B liquid effective radii (1-39 um)
at representative MODIS channels, for

   fixed      source revision c6ad545: fixed nmom = 1000 (run from the export
              made by compact_lut_legendre.py in validation/tmp/rev_c6ad545)
   adaptive   the current source

Both use the production 3.5 x effective-radius liquid integration.  Each
variant runs in its own process:

   python validation/legendre_expansion/benchmark_legendre.py [fixed|adaptive]
   python validation/legendre_expansion/benchmark_legendre.py          (both, then summary)
"""

import csv
import json
import os
import subprocess
import sys
import time
from pathlib import Path

for variable in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS"):
   os.environ[variable] = "1"

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
OLD_TREE = ROOT / "validation" / "tmp" / "rev_c6ad545"
CHANNELS = [3, 1, 6, 7, 20, 31]          # 0.469, 0.646, 1.63, 2.13, 3.79 and 11.03 um
RESULTS = HERE / "results"


def run(variant):
   root = OLD_TREE if variant == "fixed" else ROOT
   sys.path.insert(0, str(root / "src"))
   import numpy as np
   from oraclut.idl_mirror import generate_scattering_properties, load_inststr, load_lutstr, load_mmdat, load_srfstrarr

   inputs = ROOT / "create_orac_lut" / "input_files"
   rows = []
   # One call for all channels, as in a LUT run: the 0.55 um reference (whose
   # moments are not passed to DISORT) is calculated once, not per channel.
   inststr = load_inststr(inputs / "inst" / "terra_modis_v1.inst", requestedchannelid=CHANNELS)
   lutstr = load_lutstr(inputs / "lut" / "liquid-water-cloud-grid-b.lut", inststr.max_sat_zenith)
   srfstrarr, nwvl_max = load_srfstrarr(inststr, inputs / "sun" / "Gueymard2018.sssi", 1, inputs)
   mmstr = load_mmdat(inputs / "microphysics" / "liquid-water_stg.mm", inputs)
   started = time.process_time()
   result = generate_scattering_properties(srfstrarr, nwvl_max, inststr, mmstr, lutstr, 1000 if variant == "fixed" else None)
   cpu = time.process_time() - started
   lengths = np.atleast_1d(result[0])
   for l, channel in enumerate(inststr.channelid):
      per_channel = np.ravel(lengths[:, l, :]) if lengths.ndim == 3 else lengths
      rows.append({"variant": variant, "channel": int(channel), "wavelength_um": float(srfstrarr[l].wvl[0]),
                   "mean_moments": float(np.mean(per_channel)), "max_moments": int(np.max(per_channel))})
   rows.append({"variant": variant, "channel": "all (incl. 0.55 um reference)", "wavelength_um": None,
                "cpu_s": cpu, "mean_moments": float(np.mean(lengths)), "max_moments": int(np.max(lengths))})
   print(rows[-1], flush=True)
   (RESULTS / f"benchmark_legendre_{variant}.json").write_text(json.dumps(rows, indent=2) + "\n")


def main():
   RESULTS.mkdir(parents=True, exist_ok=True)
   if len(sys.argv) > 1:
      run(sys.argv[1])
      return
   for variant in ("fixed", "adaptive"):
      subprocess.run([sys.executable, __file__, variant], check=True)
   rows = []
   for variant in ("fixed", "adaptive"):
      rows += json.loads((RESULTS / f"benchmark_legendre_{variant}.json").read_text())
   with open(RESULTS / "benchmark_legendre.csv", "w", newline="") as handle:
      writer = csv.DictWriter(handle, fieldnames=["variant", "channel", "wavelength_um", "cpu_s", "mean_moments", "max_moments"])
      writer.writeheader()
      writer.writerows(rows)
   total = {v: sum(r.get("cpu_s", 0.0) for r in rows if r["variant"] == v) for v in ("fixed", "adaptive")}
   print("total scattering CPU s:", total, "ratio adaptive/fixed:", total["adaptive"] / total["fixed"])


if __name__ == "__main__":
   main()
