"""Write a machine-readable forensic summary of the full R_0v residual."""

from __future__ import annotations

import json
from pathlib import Path

import numpy as np
import xarray as xr


ROOT = Path(__file__).resolve().parents[2]
COMPARISON = ROOT / "validation" / "generated" / "python_full_comparison.json"


def stats(values: np.ndarray) -> dict[str, float | int]:
    values = np.asarray(values, dtype=np.float64).ravel()
    return {
        "count": int(values.size),
        "max_abs": float(np.max(values)),
        "rms_abs": float(np.sqrt(np.mean(values * values))),
        "p99_abs": float(np.quantile(values, 0.99)),
    }


def main() -> None:
    comparison = json.loads(COMPARISON.read_text())
    reference = xr.open_dataset(ROOT / comparison["reference"] if not Path(comparison["reference"]).is_absolute() else comparison["reference"])
    candidate = xr.open_dataset(ROOT / comparison["candidate"] if not Path(comparison["candidate"]).is_absolute() else comparison["candidate"])
    reference_da = reference["R_0v"]
    candidate_da = candidate["R_0v"].transpose(*reference_da.dims)
    ref = np.asarray(reference_da.values, dtype=np.float64)
    py = np.asarray(candidate_da.values, dtype=np.float64)
    absolute = np.abs(ref - py)
    dims = reference_da.dims
    shape = absolute.shape
    axis = {name: np.arange(size) for name, size in zip(dims, shape)}
    sensitive = (
        (axis["solar_zenith"][None, None, :, None, None, None] == 5)
        & (axis["optical_depth"][None, None, None, :, None, None] == 9)
        & (axis["effective_radius"][None, None, None, None, :, None] == 11)
        & (axis["solar_channels"][None, None, None, None, None, :] == 0)
    )
    sensitive = np.broadcast_to(sensitive, shape)
    top_indices = np.argpartition(absolute.ravel(), -50)[-50:]
    top_indices = top_indices[np.argsort(absolute.ravel()[top_indices])[::-1]]
    top = []
    for flat_index in top_indices:
        index = tuple(int(value) for value in np.unravel_index(flat_index, shape))
        top.append({
            "index": list(index),
            "coordinates": {name: float(reference[name].values[i]) for name, i in zip(dims, index)},
            "reference": float(ref[index]),
            "python": float(py[index]),
            "absolute_difference": float(absolute[index]),
            "relative_difference": float(absolute[index] / max(abs(float(ref[index])), 1.0e-12)),
        })
    summary = {
        "variable": "R_0v",
        "dimensions": list(dims),
        "shape": list(shape),
        "all_points": stats(absolute),
        "suspected_sensitive_cluster": stats(absolute[sensitive]),
        "outside_suspected_sensitive_cluster": stats(absolute[~sensitive]),
        "sensitive_cluster_definition": {
            "solar_channel_index": 0,
            "solar_zenith_index": 5,
            "optical_depth_index": 9,
            "effective_radius_index": 11,
            "warning_status": "targeted Python call emitted UPBEAM--SGECO; IDL license unavailable",
        },
        "top_50": top,
    }
    output = ROOT / "validation" / "diagnostics" / "full_r0v_residual_summary.json"
    output.write_text(json.dumps(summary, indent=2) + "\n")
    print(json.dumps(summary, indent=2))


if __name__ == "__main__":
    main()
