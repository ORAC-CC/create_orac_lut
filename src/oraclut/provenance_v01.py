"""Provenance sidecar records for generated ORAC LUT products.

Every LUT written by the generator is accompanied by
``<product>.provenance.json`` (``x_v25.nc`` -> ``x_v25.provenance.json``)
recording the LUT product version, the source release and exact Git commit,
the working-tree status, the UTC generation time, the instrument and particle
model, the identity (path and SHA-256) of every configuration and input
definition file, the scientific grid, the numerical-integration and DISORT
settings, the other scientifically significant run options, the software
versions and the SHA-256 of the product itself.  The NetCDF product is not
changed in any way by this record.

Historical products (V21-V25) have no sidecar; nothing is fabricated for them.
"""

from __future__ import annotations

import datetime
import json
import os
from pathlib import Path
from typing import Any, Mapping

import numpy as np

from .version import REPO_ROOT, sha256_file, summarise

PROVENANCE_FORMAT = 1


def provenance_path(lut_path: str | os.PathLike) -> Path:
    """``x_v25.nc`` -> ``x_v25.provenance.json`` next to the product."""

    lut_path = Path(lut_path)
    return lut_path.with_name(lut_path.stem + ".provenance.json")


def file_identity(path: str | os.PathLike | None, root: str | os.PathLike = REPO_ROOT) -> dict | None:
    """Path (relative to the repository when inside it), size and SHA-256 of a file."""

    if path is None:
        return None
    path = Path(path)
    if not path.is_file():
        return {"path": str(path), "exists": False}
    try:
        shown = str(path.resolve().relative_to(Path(root).resolve()))
    except ValueError:
        shown = str(path)
    return {"path": shown, "size_bytes": path.stat().st_size, "sha256": sha256_file(path)}


def _jsonable(value: Any) -> Any:
    if isinstance(value, Path):
        return str(value)
    if isinstance(value, np.generic):
        return value.item()
    if isinstance(value, np.ndarray):
        return value.tolist()
    if isinstance(value, Mapping):
        return {str(k): _jsonable(v) for k, v in value.items()}
    if isinstance(value, (list, tuple)):
        return [_jsonable(v) for v in value]
    return value


def build_record(
    lut_path: str | os.PathLike,
    *,
    generator: str,
    lut_version: int,
    configuration: Mapping[str, Any],
    configuration_file: str | os.PathLike | None,
    input_files: Mapping[str, str | os.PathLike | None],
    grid: Mapping[str, Any],
    numerics: Mapping[str, Any],
    radiative_transfer: Mapping[str, Any],
    extra: Mapping[str, Any] | None = None,
    root: str | os.PathLike = REPO_ROOT,
) -> dict:
    """Assemble the provenance record for a product that has been written."""

    lut_path = Path(lut_path)
    summary = summarise(root)
    stat = lut_path.stat()
    record = {
        "provenance_format": PROVENANCE_FORMAT,
        "generated_utc": datetime.datetime.now(datetime.timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ"),
        "generator": generator,
        "lut_version": int(lut_version),
        "output": {
            "file": lut_path.name,
            "directory": str(lut_path.parent),
            "size_bytes": stat.st_size,
            "sha256": sha256_file(lut_path),
        },
        "source": summary.record(),
        "configuration_file": None if configuration_file is None else {
            **(file_identity(configuration_file, root) or {}),
            "text": Path(configuration_file).read_text() if Path(configuration_file).is_file() else None,
        },
        "configuration": _jsonable(configuration),
        "input_files": {name: file_identity(path, root) for name, path in input_files.items()},
        "grid": _jsonable(grid),
        "numerics": _jsonable(numerics),
        "radiative_transfer": _jsonable(radiative_transfer),
        "job": {
            "slurm_job_id": os.environ.get("SLURM_JOB_ID"),
            "slurm_array_task_id": os.environ.get("SLURM_ARRAY_TASK_ID"),
            "slurm_job_name": os.environ.get("SLURM_JOB_NAME"),
            "user": os.environ.get("USER"),
            "working_directory": os.getcwd(),
        },
    }
    if extra:
        record["extra"] = _jsonable(extra)
    return record


def write_provenance(lut_path: str | os.PathLike, record: Mapping[str, Any]) -> Path:
    """Write ``record`` atomically as the product's sidecar and return its path."""

    target = provenance_path(lut_path)
    temporary = target.with_name(target.name + ".tmp")
    with open(temporary, "w") as handle:
        json.dump(record, handle, indent=1, sort_keys=True)
        handle.write("\n")
        handle.flush()
        os.fsync(handle.fileno())
    os.replace(temporary, target)
    return target


def read_provenance(path: str | os.PathLike) -> dict:
    return json.loads(Path(path).read_text())
