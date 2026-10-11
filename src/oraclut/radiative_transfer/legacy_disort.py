"""Python boundary for the production ORAC DISORT implementation.

Two precisions of the same unmodified DISORT source are available, selected
by the environment variable ORAC_LUT_DISORT_PRECISION exactly as in the IDL
generator's call_disort.pro:

* ``double`` (default): create_orac_lut/disort2_dp/libdisort_dp.so, the
  kernel objects of the IDL DLM DISORTtoIDL_DP (REAL promoted to 8 bytes;
  built by create_orac_lut/disort2_dp/build_dlm.sh);
* ``single``: the historical single-precision kernel built here from
  create_orac_lut/disort2/src (the V22 configuration up to lut-code-v22.1).

GETMOM and PLKAVG always use the single-precision library, as the IDL
generator calls them through the historical DISORTtoIDL DLM.
"""

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


def _disort_single(
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
    """Historical single-precision kernel (the IDL DISORT DLM equivalent)."""

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


DOUBLE_LIBRARY = _root() / "create_orac_lut" / "disort2_dp" / "libdisort_dp.so"
_DOUBLE: ctypes.CDLL | None = None


def disort_precision() -> str:
    """The selected DISORT precision: ``double`` (default) or ``single``."""

    precision = os.environ.get("ORAC_LUT_DISORT_PRECISION", "").strip().lower() or "double"
    if precision not in {"double", "single"}:
        raise ValueError(f"ORAC_LUT_DISORT_PRECISION must be single or double, not {precision!r}")
    return precision


def _load_double_library() -> ctypes.CDLL:
    global _DOUBLE
    if _DOUBLE is not None:
        return _DOUBLE
    with _BUILD_LOCK:
        if _DOUBLE is None:
            if not DOUBLE_LIBRARY.is_file():
                raise DisortError(f"{DOUBLE_LIBRARY} is missing; build it with create_orac_lut/disort2_dp/build_dlm.sh")
            library = ctypes.CDLL(str(DOUBLE_LIBRARY))
            library.disort_.restype = ctypes.c_int32
            _DOUBLE = library
    return _DOUBLE


def _disort_double(
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
    """Double-precision kernel with the argument handling of the IDL DISORT_DP DLM.

    As in IDL, every input is first held at float32 (the IDL variables are
    FLOAT) and then widened exactly to float64; the fixed arguments are those
    of setup_disort (ACCUR = float32 1e-8); MAXCLY etc. equal the actual
    sizes; TEMPER is zero-padded to NLYR + 1 values; outputs are float64.
    """

    f32 = np.float32
    dtauc = np.asarray(dtauc, dtype=f32).ravel()
    ssalb = np.asarray(single_scatter_albedo, dtype=f32).ravel()
    if dtauc.shape != ssalb.shape:
        raise ValueError("dtauc and single_scatter_albedo must have equal shapes")
    nlyr = dtauc.size
    phase = np.asarray(phase_moments, dtype=f32)
    if phase.ndim != 2 or phase.shape[1] != nlyr:
        raise ValueError("phase_moments must have shape (highest_order + 1, nlayer)")
    highest_order = phase.shape[0] - 1
    utau = np.asarray(utau, dtype=f32).ravel()
    umu = np.asarray(umu, dtype=f32).ravel()
    phi = np.asarray(phi, dtype=f32).ravel()
    if utau.size != 2 or umu.size == 0 or phi.size == 0:
        raise ValueError("production ORAC calls require two utau values and non-empty angles")
    given = np.zeros(nlyr, dtype=f32) if temperature is None else np.asarray(temperature, dtype=f32).ravel()
    if given.size not in ({nlyr + 1} if plank else {nlyr}):
        raise ValueError("temperature must have one value per layer, or one extra boundary value for a Planck call")
    if np.any(ssalb > 1.0) or np.any(ssalb < 0.0):
        raise ValueError("single scattering albedo outside [0, 1]")
    temper = np.zeros(nlyr + 1, dtype=np.float64)
    temper[:min(given.size, nlyr + 1)] = given[:nlyr + 1]
    d = lambda a: np.ascontiguousarray(a, dtype=np.float64)
    dtauc, ssalb, utau, umu, phi = d(dtauc), d(ssalb), d(utau), d(umu), d(phi)
    phase = np.asfortranarray(phase, dtype=np.float64)
    ntau, numu, nphi = utau.size, umu.size, phi.size
    out = {k: np.zeros(ntau) for k in ("rfldir", "rfldn", "flup", "dfdt", "uavg")}
    out["uu"] = np.zeros((numu, ntau, nphi), order="F")
    out["albmed"] = np.zeros(numu)
    out["trnmed"] = np.zeros(numu)
    i32 = lambda v: ctypes.byref(ctypes.c_int32(int(v)))
    real = lambda v: ctypes.byref(ctypes.c_double(float(f32(v))))
    ptr = lambda a: a.ctypes.data_as(ctypes.c_void_p)
    prnt = np.zeros(5, dtype=np.int32)
    status = _load_double_library().disort_(
        i32(nlyr), ptr(dtauc), ptr(ssalb), i32(highest_order), ptr(phase), ptr(temper),
        real(wavenumber_low), real(wavenumber_high), i32(1), i32(ntau), ptr(utau), i32(nstreams),
        i32(1), i32(numu), ptr(umu), i32(nphi), ptr(phi), i32(0), real(fbeam), real(umu0),
        real(0.0), real(fisot), i32(1), real(0.0), real(0.0), real(0.0), real(0.0),
        i32(int(bool(plank))), i32(0), real(1.0e-8), ptr(prnt), i32(nlyr), i32(ntau), i32(numu),
        i32(nphi), i32(highest_order),
        ptr(out["rfldir"]), ptr(out["rfldn"]), ptr(out["flup"]), ptr(out["dfdt"]), ptr(out["uavg"]),
        ptr(out["uu"]), ptr(out["albmed"]), ptr(out["trnmed"]),
    )
    if status:
        raise DisortError(f"DISORT (double precision) returned error code {status}")
    return out


def disort(*args, **kwargs) -> dict[str, np.ndarray]:
    """Call DISORT at the selected precision (ORAC_LUT_DISORT_PRECISION, default double)."""

    if disort_precision() == "double":
        return _disort_double(*args, **kwargs)
    return _disort_single(*args, **kwargs)


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
