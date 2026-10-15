"""Record the inputs of the V25 emission calls of a short Terra MODIS LUT run (read-only for the code).

   ORAC_LUT_DISORT_PRECISION=single ORAC_LUT_DISORT_SP_KERNEL=intel_dlm \
   python validation/v25/temperature_profiles/emission_inputs.py <model> <r_e> <out.npz> <channel> ...

Runs create_orac_cloud_lut (V25 working tree, production profile of the
substance) on Grid B's 24 optical depths at one effective radius and stores,
for every emission call, the in-cloud layer (optical depth, single-scattering
albedo, phase moments), the node's 0.55 um optical depth and the channel's
DISORT window: the realistic cases of layering.py.
"""

import importlib
import sys
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[2]
sys.path.insert(0, str(ROOT / "src"))
sys.path.insert(0, str(ROOT))
driver = importlib.import_module("create_orac_luts")
model, re, out, channels = sys.argv[1], float(sys.argv[2]), Path(sys.argv[3]).resolve(), [int(c) for c in sys.argv[4:]]
substance = model.split("_")[0]
profile = {"liquid-water": "wet_adiabat", "water-ice": "cirrostratus"}[substance]
grid = out.with_suffix(".lut")
grid.write_text(f"""# one-radius Grid B subset for emission_inputs.py
24 uneven_logarithmic # optical depth
1e-10 0.0078125 0.015625 0.03125 0.0625 0.125 0.25 0.5 1.0 2.0 2.82842712474619 4.0 5.65685424949238 8.0 11.3137084989848 16.0 22.6274169979695 32.0 45.254833995939 64.0 90.509667991878 128.0 181.019335983756 256.0
1 uneven_linear # physical effective radius, microns
{re}
1 uneven_linear # solar zenith, degrees
0.0
2 uneven_linear # satellite zenith, degrees
0.0 60.0
1 linear # relative azimuth, degrees, 0-180
0 0
""")
records = []
emission_layers = driver.emission_layers


def recording_layers(emtau, emssa, empmo, tau_055, cloudprof, temp, wnlo, wnhi):
   records.append({"emtau": np.asarray(emtau).copy(), "emssa": np.asarray(emssa).copy(), "empmo": np.asarray(empmo).copy(),
                   "tau_055": float(tau_055), "wnlo": float(wnlo), "wnhi": float(wnhi)})
   return emission_layers(emtau, emssa, empmo, tau_055, cloudprof, temp, wnlo, wnhi)


driver.emission_layers = recording_layers
status = driver.create_orac_cloud_lut(ROOT / "create_orac_lut" / "input_files", "terra_modis_v1.inst", f"{model}.mm", str(grid),
                                      out.parent / out.stem, 2, channelid=channels, gas=0, no_rayleigh=0, srf_quad=1, version=25,
                                      nstreams=60, nmom=1000, cloud_vertical_profile=profile)
if status != 0:
   raise SystemExit(f"create_orac_cloud_lut returned {status}")
# the calls run channel by channel (one SRF point each with srf_quad 1), then optical depth
if len(records) != 24 * len(channels):
   raise SystemExit(f"{len(records)} emission calls, expected {24 * len(channels)}")
nmom = max(r["empmo"].shape[0] for r in records)
moments = np.zeros((len(records), nmom), dtype=np.float32)
for k, r in enumerate(records):
   moments[k, :r["empmo"].shape[0]] = r["empmo"][:, 0]
np.savez(out, model=model, profile=profile, effective_radius=re, channel=np.repeat(channels, 24),
         tau_055=[r["tau_055"] for r in records], wnlo=[r["wnlo"] for r in records], wnhi=[r["wnhi"] for r in records],
         emtau=[r["emtau"][0] for r in records], emssa=[r["emssa"][0] for r in records],
         nmom=[r["empmo"].shape[0] for r in records], empmo=moments)
print(f"{len(records)} emission calls recorded -> {out}")
