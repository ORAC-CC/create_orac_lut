"""Compare the legacy and Python per-layer DISORT inputs for one aerosol state.

Reads ``legacy_state.sav`` (IDL SAVE written by aerosol_layer_capture.pro) and
``python_state.npz`` (written by capture_python_state.py), converts the legacy
arrays to ``legacy_state.npz``, compares every common array strictly, and, as a
controlled check, runs the production DISORT kernel from Python on BOTH sets of
captured inputs so the effect of the input assembly can be separated from the
kernel itself.  Diagnostic only: nothing in ``src/`` is modified.
"""

from __future__ import annotations

import json
from pathlib import Path

import numpy as np
from scipy.io import readsav

from oraclut.radiative_transfer import call_disort

HERE = Path(__file__).resolve().parent
INSPECT_MOMENTS = (0, 1, 2, 10, 100, 500, 999)


def stats(name: str, legacy: np.ndarray, python: np.ndarray) -> dict:
    legacy = np.asarray(legacy, dtype=np.float64)
    python = np.asarray(python, dtype=np.float64)
    entry = {"shape_legacy": list(legacy.shape), "shape_python": list(python.shape)}
    if legacy.shape != python.shape:
        entry["status"] = "shape mismatch"
        return entry
    diff = python - legacy
    absdiff = np.abs(diff)
    index = int(np.argmax(absdiff))
    floor = np.abs(legacy) >= 1e-12
    entry.update({
        "status": "ok",
        "max_abs_difference": float(absdiff.max()),
        "rms_difference": float(np.sqrt(np.mean(diff ** 2))),
        "mean_abs_difference": float(absdiff.mean()),
        "max_rel_difference_where_meaningful": float(np.max(absdiff[floor] / np.abs(legacy[floor]))) if floor.any() else None,
        "index_of_max_abs": list(int(i) for i in np.unravel_index(index, legacy.shape)),
        "legacy_at_max": float(legacy.flat[index]),
        "python_at_max": float(python.flat[index]),
        "elements_differing_above_1e-6": int(np.count_nonzero(absdiff > 1e-6)),
        "bitwise_identical": bool(np.array_equal(legacy.astype(np.float32), python.astype(np.float32))),
    })
    return entry


def main() -> int:
    sav = readsav(HERE / "legacy_state.sav")
    legacy = {
        "level_height_km": sav["level_height"], "level_pressure_hpa": sav["level_pressure"],
        "layer_height_km": sav["layer_height"], "profile_height_km": sav["profile_height"],
        "profile_relative_amount": sav["profile_rext"], "relative_tau_raw": sav["scatreltau_raw"],
        "relative_tau": sav["scatreltau"], "gas_level": sav["gas_level"], "tau_gas": sav["tau_gas"],
        "column_tau_rayleigh": np.float32(sav["column_tau_rayleigh"]), "rayleigh_level": sav["rayleigh_level"],
        "tau_rayleigh": sav["tau_rayleigh"], "tauscat": sav["tauscat"],
        "tau_aerosol_scattering": sav["tau_aerosol_scattering"], "tau_aerosol_absorption": sav["tau_aerosol_absorption"],
        "dtau": sav["dtau"], "ssa": sav["ssa"],
        # IDL PMom[NMom, NLayers] arrives transposed as (NLayers, NMom); restore (moments, layers).
        "pmom": np.ascontiguousarray(np.asarray(sav["pmom"]).T),
        "molecular_moments": sav["molecular_moments"], "particle_moments": sav["particle_moments"],
    }
    legacy_state = {
        "channel": int(sav["state_channel"]), "effective_radius_um": float(sav["state_efr"]),
        "optical_depth": float(sav["state_opd"]), "surface_pressure_hpa": float(sav["state_prs"]),
        "wavelength_um": float(sav["state_wavelength"]), "particle_ssa": float(sav["particle_ssa"]),
        "extinction_ratio": float(sav["particle_bextrat"]), "particle_asymmetry": float(sav["particle_g"]),
        "n_layers": int(sav["nlayers"]), "n_moments": int(sav["nmom"]),
    }
    np.savez(HERE / "legacy_state.npz", state=json.dumps(legacy_state), **legacy)

    py = np.load(HERE / "python_state.npz")
    python_state = json.loads(str(py["state"]))
    python_state["n_layers"] = int(py["dtau"].size)
    python_state["n_moments"] = int(py["pmom"].shape[0])

    order = [
        "level_height_km", "level_pressure_hpa", "layer_height_km", "profile_height_km",
        "profile_relative_amount", "relative_tau_raw", "relative_tau", "gas_level", "tau_gas",
        "column_tau_rayleigh", "rayleigh_level", "tau_rayleigh", "tauscat",
        "tau_aerosol_scattering", "tau_aerosol_absorption", "dtau", "ssa", "pmom",
        "molecular_moments", "particle_moments",
    ]
    comparison = {name: stats(name, legacy[name], py[name]) for name in order}

    # First array in assembly order that differs materially (> 1e-6 absolute).
    first = next((name for name in order if comparison[name].get("elements_differing_above_1e-6", 1)), None)

    pmom_l, pmom_p = legacy["pmom"], np.asarray(py["pmom"])
    moments = {}
    for m in INSPECT_MOMENTS:
        if m < pmom_l.shape[0]:
            moments[str(m)] = stats(f"pmom[{m}]", pmom_l[m], pmom_p[m])
    layers_with_aerosol = [int(i) for i in np.flatnonzero(np.asarray(legacy["tauscat"]) > 0)]

    # Layer table for the aerosol-bearing layers.
    rows = []
    for i in range(legacy["dtau"].size):
        rows.append({
            "layer": i, "height_km": float(legacy["layer_height_km"][i]),
            "relative_tau_legacy": float(legacy["relative_tau"][i]), "relative_tau_python": float(py["relative_tau"][i]),
            "tauscat_legacy": float(legacy["tauscat"][i]), "tauscat_python": float(py["tauscat"][i]),
            "tau_rayleigh": float(legacy["tau_rayleigh"][i]), "tau_gas": float(legacy["tau_gas"][i]),
            "dtau_legacy": float(legacy["dtau"][i]), "dtau_python": float(py["dtau"][i]),
            "ssa_legacy": float(legacy["ssa"][i]), "ssa_python": float(py["ssa"][i]),
            "pmom1_legacy": float(pmom_l[1, i]), "pmom1_python": float(pmom_p[1, i]),
        })

    # Controlled DISORT check: identical kernel and call settings as the pipeline's
    # diffuse calculation (FBeam 0, FIsot 100, umu0 = cos 50 deg), R_dd = flup / (100 pi).
    saz = np.asarray([0.0, 90.0], dtype=np.float32); raa = np.asarray([0.0, 180.0], dtype=np.float32)
    values = saz.copy(); values[values == 90.0] = 89.99
    down = -np.cos(np.deg2rad(values)).astype(np.float32); down[np.abs(down) == 1.0] *= np.float32(0.99999)
    umu = np.concatenate((down, -down[::-1])).astype(np.float32)
    umu0 = np.float32(np.cos(np.deg2rad(np.float32(50.0))))

    def r_dd(arrays):
        dtau = np.asarray(arrays["dtau"], np.float32); ssa = np.asarray(arrays["ssa"], np.float32)
        pmom = np.asfortranarray(np.asarray(arrays["pmom"], np.float32))
        utau = np.asarray([0.0, np.sum(dtau, dtype=np.float32)], np.float32)
        out = call_disort(dtau, ssa, pmom, utau, umu, raa, 0.0, umu0, 100.0, nstreams=60)
        return float(out["flup"][0] / (100.0 * np.pi)), float(out["rfldn"][1] / (100.0 * np.pi))

    hybrid = dict(py)
    # Python inputs but with the legacy profile weighting substituted (diagnostic only).
    w = python_state["particle_ssa"]; ratio = python_state["extinction_ratio"]
    tauscat_h = np.float32(python_state["optical_depth"]) * np.asarray(legacy["relative_tau"], np.float32) * np.float32(ratio)
    dtau_h = (np.asarray(py["tau_gas"]) + np.asarray(py["tau_rayleigh"]) + tauscat_h).astype(np.float32)
    ssa_h = np.divide(np.asarray(py["tau_rayleigh"]) + w * tauscat_h, dtau_h, out=np.zeros_like(dtau_h), where=dtau_h != 0).astype(np.float32)
    ssa_h[ssa_h > 1.0] = 1.0
    mol = np.asarray(py["molecular_moments"]); part = np.asarray(py["particle_moments"])
    pmom_h = np.empty_like(pmom_p)
    for i in range(dtau_h.size):
        st = np.asarray(py["tau_rayleigh"])[i] + w * tauscat_h[i]
        pmom_h[:, i] = mol if tauscat_h[i] <= 0 or st == 0 else (mol * np.asarray(py["tau_rayleigh"])[i] + part * w * tauscat_h[i]) / st
    pmom_h[pmom_h > 1.0] = 1.0
    hybrid = {"dtau": dtau_h, "ssa": ssa_h, "pmom": pmom_h}

    disort_check = {
        "R_dd_T_dd_from_legacy_arrays": r_dd(legacy),
        "R_dd_T_dd_from_python_arrays": r_dd(py),
        "R_dd_T_dd_from_python_arrays_with_legacy_profile_weights": r_dd(hybrid),
        "legacy_product_R_dd": 0.038124, "python_product_R_dd": 0.037788,
        "note": "production DISORT kernel invoked from Python on each captured input set; R_dd = FLUP/(FIsot*pi)",
    }

    report = {
        "state": {"legacy": legacy_state, "python": python_state},
        "arrays": comparison, "pmom_moments": moments,
        "layers_with_aerosol": layers_with_aerosol, "layer_table": rows,
        "first_material_difference": first,
        "disort_check": disort_check,
    }
    (HERE / "comparison.json").write_text(json.dumps(report, indent=2) + "\n")
    print(json.dumps({"first_material_difference": first, "disort_check": disort_check}, indent=2))
    for name in order:
        c = comparison[name]
        print(f"{name:28s} {c.get('status'):>8s} max|d|={c.get('max_abs_difference', float('nan')):.3e} "
              f"rms={c.get('rms_difference', float('nan')):.3e} n>1e-6={c.get('elements_differing_above_1e-6')} bitwise={c.get('bitwise_identical')}")
    print("aerosol layers:", layers_with_aerosol)
    for row in rows:
        if row["layer"] in layers_with_aerosol or row["layer"] in (layers_with_aerosol[-1] + 1 if layers_with_aerosol else 0,):
            print(f"  layer {row['layer']:2d} h={row['height_km']:5.2f} km  relTau L/P {row['relative_tau_legacy']:.5f}/{row['relative_tau_python']:.5f}  "
                  f"tauscat L/P {row['tauscat_legacy']:.5f}/{row['tauscat_python']:.5f}  ray {row['tau_rayleigh']:.5f} gas {row['tau_gas']:.5f}  "
                  f"ssa L/P {row['ssa_legacy']:.5f}/{row['ssa_python']:.5f}  pmom1 L/P {row['pmom1_legacy']:.5f}/{row['pmom1_python']:.5f}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
