"""Python access to the preserved ORAC Mie scientific kernel.

The historical IDL DLM is only an IDL registration layer around ``mieint``.
This module builds a repository-local shared library containing the unchanged
Fortran kernel and a small C ABI adapter.  No IDL runtime is loaded.
"""

from __future__ import annotations

from pathlib import Path
import ctypes
import os
import shlex
import subprocess
import threading

import numpy as np


class MieKernelError(RuntimeError):
    """The preserved Mie kernel rejected an input or overflowed."""


_BUILD_LOCK = threading.Lock()
_LIBRARY: ctypes.CDLL | None = None


def _repository_root() -> Path:
    return Path(__file__).resolve().parents[3]


def _library_path() -> Path:
    compiler = Path(os.environ.get("ORACLUT_FORTRAN", "gfortran")).name
    tag = os.environ.get("ORACLUT_FORTRAN_TAG", compiler)
    return _repository_root() / "build" / "oraclut" / f"liboraclut_mie_{tag}.so"


def _build_library(path: Path) -> None:
    root = _repository_root()
    source = root / "mie" / "dlm-code" / "mieint.f"
    adapter = root / "src" / "oraclut" / "optics" / "legacy_mie.c"
    if not source.is_file() or not adapter.is_file():
        raise MieKernelError("preserved Mie kernel sources are unavailable")
    path.parent.mkdir(parents=True, exist_ok=True)
    compiler = os.environ.get("ORACLUT_FORTRAN", "gfortran")
    if Path(compiler).name in {"ifort", "ifx"}:
        fortran_flags = ["-O3", "-fPIC", "-extend_source"]
    else:
        fortran_flags = ["-O3", "-fPIC", "-ffixed-line-length-0"]
    fortran_object = path.with_name("mieint.o")
    adapter_object = path.with_name("legacy_mie_adapter.o")
    subprocess.run(
        [compiler, *fortran_flags, *shlex.split(os.environ.get("ORACLUT_FORTRAN_FLAGS", "")), "-c",
         str(source), "-o", str(fortran_object)],
        check=True,
    )
    subprocess.run(
        ["gcc", "-O3", "-fPIC", "-std=c11", "-c", str(adapter), "-o", str(adapter_object)],
        check=True,
    )
    subprocess.run(
        [compiler, *shlex.split(os.environ.get("ORACLUT_FORTRAN_FLAGS", "")), "-shared", "-o", str(path), str(adapter_object),
         str(fortran_object)],
        check=True,
    )


def _load_library() -> ctypes.CDLL:
    global _LIBRARY
    if _LIBRARY is not None:
        return _LIBRARY
    with _BUILD_LOCK:
        if _LIBRARY is None:
            path = _library_path()
            if not path.is_file():
                _build_library(path)
            library = ctypes.CDLL(str(path))
            library.oraclut_mie_batch.argtypes = [
                ctypes.c_int32,
                np.ctypeslib.ndpointer(np.float64, flags="C_CONTIGUOUS"),
                ctypes.c_double, ctypes.c_double, ctypes.c_int32,
                np.ctypeslib.ndpointer(np.float64, flags="C_CONTIGUOUS"),
                np.ctypeslib.ndpointer(np.float64, flags="C_CONTIGUOUS"),
                np.ctypeslib.ndpointer(np.float64, flags="C_CONTIGUOUS"),
                np.ctypeslib.ndpointer(np.float64, flags="C_CONTIGUOUS"),
                np.ctypeslib.ndpointer(np.float64, flags="C_CONTIGUOUS"),
                np.ctypeslib.ndpointer(np.float64, flags="C_CONTIGUOUS"),
                np.ctypeslib.ndpointer(np.float64, flags="C_CONTIGUOUS"),
                np.ctypeslib.ndpointer(np.float64, flags="C_CONTIGUOUS"),
                np.ctypeslib.ndpointer(np.float64, flags="C_CONTIGUOUS"),
            ]
            library.oraclut_mie_batch.restype = ctypes.c_int
            _LIBRARY = library
    return _LIBRARY


def mie_single_batch(
    size_parameters: np.ndarray,
    refractive_index: complex,
    scattering_cosines: np.ndarray,
) -> dict[str, np.ndarray]:
    """Evaluate the preserved ``mieint`` kernel for a size-parameter array."""

    dx = np.ascontiguousarray(np.asarray(size_parameters, dtype=np.float64).ravel())
    dqv = np.ascontiguousarray(np.asarray(scattering_cosines, dtype=np.float64).ravel())
    if dx.size == 0 or dqv.size == 0:
        raise ValueError("size_parameters and scattering_cosines must be non-empty")
    outputs = [np.empty(dx.size, dtype=np.float64) for _ in range(4)]
    phases = [np.empty((dx.size, dqv.size), dtype=np.float64) for _ in range(4)]
    status = _load_library().oraclut_mie_batch(
        ctypes.c_int32(dx.size), dx,
        ctypes.c_double(float(np.real(refractive_index))),
        ctypes.c_double(float(np.imag(refractive_index))),
        ctypes.c_int32(dqv.size), dqv,
        *outputs, *phases,
    )
    if status:
        raise MieKernelError(f"mieint returned error code {status}")
    return {
        "qext": outputs[0], "qsca": outputs[1], "qbsc": outputs[2], "g": outputs[3],
        "f11": phases[0], "f33": phases[1], "f12": phases[2], "f34": phases[3],
    }
