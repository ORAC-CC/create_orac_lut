"""Python boundary for the current production ORAC DISORT implementation."""

from __future__ import annotations

from pathlib import Path
import ctypes
import os
import shlex
import subprocess
import threading

import numpy as np


class DisortError(RuntimeError):
    """The preserved DISORT kernel reported an error."""


_BUILD_LOCK = threading.Lock()
_LIBRARY: ctypes.CDLL | None = None


def _root() -> Path:
    return Path(__file__).resolve().parents[3]


def _library_path() -> Path:
    compiler = Path(os.environ.get("ORACLUT_FORTRAN", "gfortran")).name
    tag = os.environ.get("ORACLUT_FORTRAN_TAG", compiler)
    return _root() / "build" / "oraclut" / f"liboraclut_disort_{tag}.so"


def _build_library(path: Path) -> None:
    root = _root()
    source_root = root / "create_orac_lut" / "disort2" / "src"
    adapter = root / "src" / "oraclut" / "radiative_transfer" / "legacy_disort.c"
    path.parent.mkdir(parents=True, exist_ok=True)
    compiler = os.environ.get("ORACLUT_FORTRAN", "gfortran")
    if Path(compiler).name in {"ifort", "ifx"}:
        # Match the preserved DISORT Makefile's Intel target.  The legacy
        # makefile passes this compatibility define to every Fortran unit.
        fortran_flags = ["-O3", "-fPIC", "-extend_source", "-D__GFORTRAN__"]
    else:
        fortran_flags = ["-O2", "-fPIC", "-std=legacy", "-ffixed-line-length-0"]
    objects = []
    for source in (
        source_root / "DISORTfunctions.f",
        source_root / "BDREF.f",
        source_root / "PRTFIN.f",
        source_root / "disort" / "ErrPack.f",
        source_root / "disort" / "LINPAK.f",
        source_root / "disort" / "RDI1MACH.f",
    ):
        obj = path.with_name(source.stem.lower() + ".o")
        subprocess.run(
            [compiler, *fortran_flags, *shlex.split(os.environ.get("ORACLUT_FORTRAN_FLAGS", "")),
             "-c", str(source), "-o", str(obj)],
            check=True,
        )
        objects.append(obj)
    adapter_object = path.with_name("legacy_disort_adapter.o")
    ambrals_object = path.with_name("ambrals_fortran.o")
    ambrals_c_object = path.with_name("ambrals_c.o")
    subprocess.run(
        ["gcc", "-O2", "-fPIC", "-std=c11", "-c", str(adapter), "-o", str(adapter_object)],
        check=True,
    )
    subprocess.run(
        ["gcc", "-O2", "-fPIC", "-D_GNU_SOURCE", "-std=c11", "-c",
         str(source_root / "ambrals-fortran.c"), "-o", str(ambrals_object)],
        check=True,
    )
    subprocess.run(
        ["gcc", "-O2", "-fPIC", "-D_GNU_SOURCE", "-std=c11", "-c",
         str(source_root / "ambralsfor.c"), "-o", str(ambrals_c_object)],
        check=True,
    )
    subprocess.run(
        [compiler, *shlex.split(os.environ.get("ORACLUT_FORTRAN_FLAGS", "")), "-shared", "-o", str(path), str(adapter_object),
         str(ambrals_object), str(ambrals_c_object), *(str(obj) for obj in objects)],
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
            f32 = np.ctypeslib.ndpointer(np.float32, flags="C_CONTIGUOUS")
            f32f = np.ctypeslib.ndpointer(np.float32, flags="F_CONTIGUOUS")
            i32 = ctypes.c_int32
            library.oraclut_disort.argtypes = [
                i32, f32, f32, i32, f32f, f32, ctypes.c_float, ctypes.c_float,
                i32, i32, f32, i32, i32, f32, i32, f32,
                ctypes.c_float, ctypes.c_float, ctypes.c_float,
                f32, f32, f32, f32, f32, f32f, f32, f32,
            ]
            library.oraclut_disort.restype = ctypes.c_int
            library.oraclut_getmom.argtypes = [i32, ctypes.c_float, i32, f32]
            library.oraclut_getmom.restype = ctypes.c_int
            library.oraclut_plkavg.argtypes = [ctypes.c_float, ctypes.c_float, ctypes.c_float]
            library.oraclut_plkavg.restype = ctypes.c_float
            _LIBRARY = library
    return _LIBRARY


def getmom(iphas: int, asymmetry: float, highest_order: int) -> np.ndarray:
    """Return the exact production DISORT ``GETMOM`` coefficients."""

    if highest_order < 2:
        raise ValueError("highest_order must be at least 2")
    values = np.empty(highest_order + 1, dtype=np.float32)
    status = _load_library().oraclut_getmom(
        ctypes.c_int32(iphas), ctypes.c_float(asymmetry),
        ctypes.c_int32(highest_order), values,
    )
    if status:
        raise DisortError(f"GETMOM returned error code {status}")
    return values


def plkavg(wavenumber_low: float, wavenumber_high: float, temperature: float) -> float:
    """Return the exact production DISORT ``PLKAVG`` value."""

    return float(_load_library().oraclut_plkavg(
        ctypes.c_float(wavenumber_low), ctypes.c_float(wavenumber_high),
        ctypes.c_float(temperature),
    ))


def disort(
    dtauc: np.ndarray,
    single_scatter_albedo: np.ndarray,
    phase_moments: np.ndarray,
    utau: np.ndarray,
    umu: np.ndarray,
    phi: np.ndarray,
    fbeam: float,
    umu0: float,
    fisot: float,
    *,
    nstreams: int = 60,
    plank: bool = False,
    wavenumber_low: float = 0.0,
    wavenumber_high: float = 0.0,
    temperature: np.ndarray | None = None,
) -> dict[str, np.ndarray]:
    """Call the production kernel with the same argument contract as IDL."""

    dtauc = np.ascontiguousarray(np.asarray(dtauc, dtype=np.float32).ravel())
    ssalb = np.ascontiguousarray(np.asarray(single_scatter_albedo, dtype=np.float32).ravel())
    if dtauc.shape != ssalb.shape:
        raise ValueError("dtauc and single_scatter_albedo must have equal shapes")
    nlyr = dtauc.size
    phase = np.asfortranarray(np.asarray(phase_moments, dtype=np.float32))
    if phase.shape != (phase.shape[0], nlyr):
        raise ValueError("phase_moments must have shape (highest_order + 1, nlayer)")
    highest_order = phase.shape[0] - 1
    utau = np.ascontiguousarray(np.asarray(utau, dtype=np.float32).ravel())
    umu = np.ascontiguousarray(np.asarray(umu, dtype=np.float32).ravel())
    phi = np.ascontiguousarray(np.asarray(phi, dtype=np.float32).ravel())
    if utau.size != 2 or umu.size == 0 or phi.size == 0:
        raise ValueError("production ORAC calls require two utau values and non-empty angles")
    temper = np.ascontiguousarray(
        np.zeros(nlyr, dtype=np.float32) if temperature is None
        else np.asarray(temperature, dtype=np.float32).ravel()
    )
    expected_temperature_sizes = {nlyr + 1} if plank else {nlyr}
    if temper.size not in expected_temperature_sizes:
        raise ValueError(
            "temperature must have one value per layer, or one extra boundary "
            "value for a Planck call"
        )
    ntau = utau.size
    numu = umu.size
    nphi = phi.size
    rfldir = np.zeros(ntau, dtype=np.float32)
    rfldn = np.zeros(ntau, dtype=np.float32)
    flup = np.zeros(ntau, dtype=np.float32)
    dfdt = np.zeros(ntau, dtype=np.float32)
    uavg = np.zeros(ntau, dtype=np.float32)
    uu = np.zeros((numu, ntau, nphi), dtype=np.float32, order="F")
    albmed = np.zeros(numu, dtype=np.float32)
    trnmed = np.zeros(numu, dtype=np.float32)
    status = _load_library().oraclut_disort(
        ctypes.c_int32(nlyr), dtauc, ssalb, ctypes.c_int32(highest_order), phase,
        temper, ctypes.c_float(wavenumber_low), ctypes.c_float(wavenumber_high),
        ctypes.c_int32(int(plank)), ctypes.c_int32(ntau), utau,
        ctypes.c_int32(nstreams), ctypes.c_int32(numu), umu, ctypes.c_int32(nphi), phi,
        ctypes.c_float(fbeam), ctypes.c_float(umu0), ctypes.c_float(fisot),
        rfldir, rfldn, flup, dfdt, uavg, uu, albmed, trnmed,
    )
    if status:
        raise DisortError(f"DISORT returned error code {status}")
    return {
        "rfldir": rfldir, "rfldn": rfldn, "flup": flup, "dfdt": dfdt,
        "uavg": uavg, "uu": uu, "albmed": albmed, "trnmed": trnmed,
    }


def call_disort(
    dtauc: np.ndarray,
    single_scatter_albedo: np.ndarray,
    phase_moments: np.ndarray,
    utau: np.ndarray,
    umu: np.ndarray,
    phi: np.ndarray,
    fbeam: float,
    umu0: float,
    fisot: float,
    *,
    nstreams: int = 60,
    plank: bool = False,
    wavenumber_low: float = 0.0,
    wavenumber_high: float = 0.0,
    temperature: np.ndarray | None = None,
) -> dict[str, np.ndarray]:
    """Reproduce ``call_disort.pro`` layer compression before the kernel call."""

    if plank:
        return disort(
            dtauc, single_scatter_albedo, phase_moments, utau, umu, phi,
            fbeam, umu0, fisot, nstreams=nstreams, plank=True,
            wavenumber_low=wavenumber_low, wavenumber_high=wavenumber_high,
            temperature=temperature,
        )
    dtauc = np.asarray(dtauc, dtype=np.float32).ravel()
    ssalb = np.asarray(single_scatter_albedo, dtype=np.float32).ravel()
    moments = np.asarray(phase_moments, dtype=np.float32)
    compressed_dtau = [dtauc[0]]
    compressed_ssalb = [ssalb[0]]
    compressed_moments = [moments[:, 0]]
    for layer in range(1, dtauc.size):
        if ssalb[layer] == compressed_ssalb[-1] and np.array_equal(
            moments[:, layer], compressed_moments[-1]
        ):
            compressed_dtau[-1] = compressed_dtau[-1] + dtauc[layer]
        else:
            compressed_dtau.append(dtauc[layer])
            compressed_ssalb.append(ssalb[layer])
            compressed_moments.append(moments[:, layer])
    return disort(
        np.asarray(compressed_dtau, dtype=np.float32),
        np.asarray(compressed_ssalb, dtype=np.float32),
        np.asfortranarray(np.column_stack(compressed_moments)),
        utau, umu, phi, fbeam, umu0, fisot, nstreams=nstreams,
    )
