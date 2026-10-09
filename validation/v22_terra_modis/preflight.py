"""Read-only preflight of the V22 Terra MODIS run files (no Mie or DISORT work).

    python validation/v22_terra_modis/preflight.py runs/<a>.run [...]

Parses each run file with the production reader (create_orac_luts.read_runfile),
checks the LUT version against CODE_VERSION, resolves the instrument, LUT grid
and microphysics with the production loaders, and prints the resolved settings
and the expected product and provenance paths, refusing if either exists.
"""

import sys
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT)); sys.path.insert(0, str(ROOT / "src"))
import create_orac_luts as c                                     # noqa: E402
from oraclut.idl_mirror import load_inststr, load_lutstr, load_mmdat   # noqa: E402
from oraclut.version import check_lut_version, summarise         # noqa: E402

summary = summarise(ROOT)
print(summary.banner())
status = 0
for run in sys.argv[1:]:
    s = c.read_runfile(run)
    check_lut_version(s["version"], summary.code, source=run)
    in_path = ROOT / s["in_path"]
    inst = load_inststr(in_path / "inst" / s["instfile"])
    lut = load_lutstr(in_path / "lut" / s["lutfile"], inst.max_sat_zenith)
    mm = load_mmdat(in_path / "microphysics" / s["mmfile"], in_path)
    bad = [ch for ch in s["channelid"] if ch not in [int(x) for x in inst.channelid]]
    name = f"{inst.platform}_{inst.instrument}_m_{mm.substance.lower()}_a01_p{mm.shortname.lower()}_v{int(s['version']):02d}"
    out = Path(s["out_path"])
    print(f"\n=== {run}")
    for k in ("platform", "instrument", "forward_model", "instfile", "mmfile", "lutfile", "atmospheres", "gas",
              "no_rayleigh", "srf_quad", "nstreams", "nmom", "version", "reuse_scat", "scat_only"):
        print(f"  {k:<14} {s[k]}")
    print(f"  channels       {len(s['channelid'])}: {s['channelid'][0]}..{s['channelid'][-1]}  invalid: {bad}")
    print(f"  substance      {mm.substance} / {mm.shortname}; components {mm.comptype} {mm.compname}")
    for ax, lab in (("opd", "optical depth"), ("efr", "effective radius"), ("soz", "solar zenith"),
                    ("saz", "satellite zenith"), ("raa", "relative azimuth")):
        v = np.asarray(getattr(lut, ax))
        print(f"  {lab:<17} n={v.size:<3} {getattr(lut, ax + '_spacing'):<19} {np.array2string(v, precision=5, max_line_width=200)}")
    for suffix in (".nc", ".provenance.json"):
        p = out / (name + suffix)
        state = "EXISTS - REFUSE" if p.exists() else "free"
        status |= p.exists()
        print(f"  output         {p}  [{state}]")
    status |= bool(bad)
sys.exit(status)
