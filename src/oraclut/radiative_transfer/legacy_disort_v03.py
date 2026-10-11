"""Python boundary for the production ORAC DISORT implementation.

Two precisions of the same unmodified DISORT source are available, selected
by the environment variable ORAC_LUT_DISORT_PRECISION exactly as in the IDL
generator's call_disort.pro:

* ``double`` (default): create_orac_lut/disort2_dp/libdisort_dp.so, the
  kernel objects of the IDL DLM DISORTtoIDL_DP (REAL promoted to 8 bytes;
  built by create_orac_lut/disort2_dp/build_dlm.sh);
* ``single``: the historical single-precision DISORT.  Two compiled kernels
  of the same source exist, selected by ORAC_LUT_DISORT_SP_KERNEL:

  - ``intel_dlm`` (preferred; used when the variable is unset): the original
    V21 IDL DLM create_orac_lut/disort2/src/DISORTtoIDL.so (ifort -O3),
    called unmodified with IDL's argument memory layout, which reproduces
    the V21 IDL DISORT calls bitwise;
  - ``gfortran``: the kernel built here from create_orac_lut/disort2/src with
    gfortran -O2 (the V22 configuration up to lut-code-v22.1), kept as a
    compatibility backend.  For identical inputs it differs from the Intel
    kernel.

GETMOM and PLKAVG are single-precision routines of the same libraries.  With
single precision and the Intel kernel they come from DISORTtoIDL.so, as in the
IDL generator; otherwise (gfortran kernel, or double precision, unchanged
from lut-code-v22.1 and the double-precision development) from the gfortran
library.
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
    if _intel_single_selected():
        return _getmom_intel(iphas, asymmetry, highest_order)
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

    if _intel_single_selected():
        return _plkavg_intel(wavenumber_low, wavenumber_high, temperature)
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


# ---------------------------------------------------------------------------
# The original V21 Intel kernel (create_orac_lut/disort2/src/DISORTtoIDL.so)
# ---------------------------------------------------------------------------
#
# The IDL DLM is an ifort -O3 build whose vectorised loops peel a number of
# elements that depends on the alignment of the array arguments (in
# particular UMU), which changes the summation order.  IDL array data starts
# IDL_OFFSET bytes past an IDL_ALIGNMENT-byte boundary; every array argument
# is therefore copied into a buffer placed in the same way.  With this layout
# the kernel reproduces all 90 captured V21 IDL calls bitwise (offsets 8, 24,
# 40 and 56 modulo 64 all do; 0, 4 and 12 do not:
# validation/v21_v22_terra_modis/v21_v22_single/REPORT_v21_v22_single_reproduction.md).
#
# The DLM's Fortran entry points are called directly; its C wrapper only
# forwards IDL pointers (DISORTwrapper.c) and its IDL_* imports are never
# used, so the library is opened with lazy binding.  It needs the Intel
# 2022 Fortran run-time libraries (module intel-compilers/2022, which the
# V21 IDL runs load); they are opened from an explicit directory so that the
# run time actually used is known.

INTEL_DLM_LIBRARY = _root() / "create_orac_lut" / "disort2" / "src" / "DISORTtoIDL.so"
INTEL_RUNTIME_DIRECTORY = Path(
    "/network/software/linux-x86_64/compilers/intel/2022/compiler/2022.1.0/linux/compiler/lib/intel64_lin"
)
# Dependency order: each library's own NEEDED entries are opened before it.
INTEL_RUNTIME_LIBRARIES = (
    "libintlc.so.5", "libirc.so", "libsvml.so", "libimf.so",
    "libifcore.so.5", "libifport.so.5", "libifcoremt.so.5",
)
IDL_ALIGNMENT = 16
IDL_OFFSET = 8
SINGLE_KERNELS = ("intel_dlm", "gfortran")
_INTEL: ctypes.CDLL | None = None


def single_kernel() -> str:
    """The single-precision kernel: ``intel_dlm`` (preferred, default) or ``gfortran``."""

    kernel = os.environ.get("ORAC_LUT_DISORT_SP_KERNEL", "").strip().lower() or "intel_dlm"
    if kernel not in SINGLE_KERNELS:
        raise ValueError(f"ORAC_LUT_DISORT_SP_KERNEL must be one of {SINGLE_KERNELS}, not {kernel!r}")
    return kernel


def intel_runtime_directory() -> Path:
    """Directory of the Intel run-time libraries (ORAC_LUT_INTEL_RUNTIME overrides the 2022.1.0 default)."""

    override = os.environ.get("ORAC_LUT_INTEL_RUNTIME", "").strip()
    return Path(override) if override else INTEL_RUNTIME_DIRECTORY


def single_library() -> Path:
    """Path of the single-precision kernel library selected by ORAC_LUT_DISORT_SP_KERNEL."""

    return INTEL_DLM_LIBRARY if single_kernel() == "intel_dlm" else _library_path()


def _intel_single_selected() -> bool:
    return disort_precision() == "single" and single_kernel() == "intel_dlm"


def _load_intel_library() -> ctypes.CDLL:
    global _INTEL
    if _INTEL is not None:
        return _INTEL
    with _BUILD_LOCK:
        if _INTEL is None:
            if not INTEL_DLM_LIBRARY.is_file():
                raise DisortError(f"original Intel DISORT kernel {INTEL_DLM_LIBRARY} is missing")
            runtime = intel_runtime_directory()
            for name in INTEL_RUNTIME_LIBRARIES:
                if not (runtime / name).is_file():
                    raise DisortError(
                        f"Intel run-time library {runtime / name} is missing (needed by {INTEL_DLM_LIBRARY.name}; "
                        "set ORAC_LUT_INTEL_RUNTIME to the intel-compilers/2022 lib/intel64_lin directory)"
                    )
                # RTLD_GLOBAL: libifport uses libifcore symbols without
                # listing it as NEEDED (resolved in the DLM's load scope when
                # the run time is found through LD_LIBRARY_PATH).  Standard
                # math names still resolve to glibc libm, which the Python
                # executable loads first, as the IDL executable does.
                ctypes.CDLL(str(runtime / name), mode=os.RTLD_LAZY | os.RTLD_GLOBAL)
            library = ctypes.CDLL(str(INTEL_DLM_LIBRARY), mode=os.RTLD_LAZY | os.RTLD_LOCAL)
            _check_intel_runtime(runtime)
            library.disort_.restype = ctypes.c_int32
            library.getmom_.restype = ctypes.c_int32
            library.plkavg_.restype = ctypes.c_float
            _INTEL = library
    return _INTEL


def _check_intel_runtime(runtime: Path) -> None:
    """Fail if any Intel run-time library mapped into the process is not from ``runtime``."""

    stems = {name.split(".so")[0] for name in INTEL_RUNTIME_LIBRARIES}
    with open("/proc/self/maps") as maps:
        mapped = {line.split()[-1] for line in maps if line.rstrip().endswith((".so", ".so.5"))}
    for path in mapped:
        if Path(path).name.split(".so")[0] in stems and Path(path).resolve().parent != runtime.resolve():
            raise DisortError(f"Intel run-time library {path} is loaded instead of the one in {runtime}")


def idl_placed(values, order: str = "C") -> np.ndarray:
    """Return a float32 copy of ``values`` whose data start IDL_OFFSET bytes past an IDL_ALIGNMENT boundary.

    The copy owns a private buffer (kept alive by the returned array), so the
    kernel may modify it in place without touching the caller's array.
    """

    values = np.asarray(values, dtype=np.float32)
    buffer = np.zeros(values.nbytes + 2 * IDL_ALIGNMENT, dtype=np.uint8)
    start = (IDL_OFFSET - buffer.ctypes.data) % IDL_ALIGNMENT
    array = np.ndarray(values.shape, dtype=np.float32, buffer=buffer, offset=start, order=order)
    array[...] = values
    if array.ctypes.data % IDL_ALIGNMENT != IDL_OFFSET:
        raise AssertionError("IDL array placement failed")
    return array


def _disort_single_intel(
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
    """The original V21 Intel kernel, called as the IDL DLM wrapper calls it.

    Arguments are those of DISORTwrapper.c: float32 REAL arrays and scalars,
    int32 INTEGER and LOGICAL (0/1) values, PMOM(0:NMOM, NLYR) and
    UU(NUMU, NTAU, NPHI) in Fortran order, MAXCLY/MAXULV/MAXUMU/MAXPHI/MAXMOM
    equal to the actual sizes, the fixed setup_disort settings (USRTAU,
    USRANG, LAMBER = 1; IBCND, ONLYFL = 0; PHI0, ALBEDO, BTEMP, TTEMP,
    TEMIS = 0; ACCUR = 1e-8; PRNT all false) and TEMPER(0:NLYR).  Every array,
    outputs included, is a fresh copy placed as IDL places array data;
    DISORT's in-place change of SSALB = 1 (to 1 - DITHER) therefore stays
    inside the copy, as it stays inside call_disort.pro's copy in IDL.
    Scalars are fresh variables for each call: DISORT's first-call self-test
    temporarily overwrites its scalar arguments.
    """

    f32 = np.float32
    dtauc_in = np.asarray(dtauc, dtype=f32).ravel()
    ssalb_in = np.asarray(single_scatter_albedo, dtype=f32).ravel()
    if dtauc_in.shape != ssalb_in.shape:
        raise ValueError("dtauc and single_scatter_albedo must have equal shapes")
    nlyr = dtauc_in.size
    phase_in = np.asarray(phase_moments, dtype=f32)
    if phase_in.ndim != 2 or phase_in.shape[1] != nlyr:
        raise ValueError("phase_moments must have shape (highest_order + 1, nlayer)")
    highest_order = phase_in.shape[0] - 1
    utau_in = np.asarray(utau, dtype=f32).ravel()
    umu_in = np.asarray(umu, dtype=f32).ravel()
    phi_in = np.asarray(phi, dtype=f32).ravel()
    if utau_in.size != 2 or umu_in.size == 0 or phi_in.size == 0:
        raise ValueError("production ORAC calls require two utau values and non-empty angles")
    given = np.zeros(nlyr, dtype=f32) if temperature is None else np.asarray(temperature, dtype=f32).ravel()
    if given.size not in ({nlyr + 1} if plank else {nlyr}):
        raise ValueError("temperature must have one value per layer, or one extra boundary value for a Planck call")
    temper_in = np.zeros(nlyr + 1, dtype=f32)
    temper_in[:given.size] = given

    dtauc_a = idl_placed(dtauc_in)
    ssalb_a = idl_placed(ssalb_in)
    pmom_a = idl_placed(phase_in, order="F")
    temper_a = idl_placed(temper_in)
    utau_a = idl_placed(utau_in)
    umu_a = idl_placed(umu_in)
    phi_a = idl_placed(phi_in)
    ntau, numu, nphi = utau_a.size, umu_a.size, phi_a.size
    out = {name: idl_placed(np.zeros(ntau, dtype=f32)) for name in ("rfldir", "rfldn", "flup", "dfdt", "uavg")}
    out["uu"] = idl_placed(np.zeros((numu, ntau, nphi), dtype=f32), order="F")
    out["albmed"] = idl_placed(np.zeros(numu, dtype=f32))
    out["trnmed"] = idl_placed(np.zeros(numu, dtype=f32))
    prnt = np.zeros(5, dtype=np.int32)

    integer = ctypes.c_int32
    real = lambda value: ctypes.c_float(float(f32(value)))
    nlyr_v, nmom_v, ntau_v, nstr_v = integer(nlyr), integer(highest_order), integer(ntau), integer(int(nstreams))
    numu_v, nphi_v = integer(numu), integer(nphi)
    usrtau_v, usrang_v, ibcnd_v, lamber_v = integer(1), integer(1), integer(0), integer(1)
    plank_v, onlyfl_v = integer(int(bool(plank))), integer(0)
    maxcly_v, maxulv_v, maxumu_v, maxphi_v, maxmom_v = integer(nlyr), integer(ntau), integer(numu), integer(nphi), integer(highest_order)
    wvnmlo_v, wvnmhi_v = real(wavenumber_low), real(wavenumber_high)
    fbeam_v, umu0_v, phi0_v, fisot_v = real(fbeam), real(umu0), real(0.0), real(fisot)
    albedo_v, btemp_v, ttemp_v, temis_v, accur_v = real(0.0), real(0.0), real(0.0), real(0.0), real(1.0e-8)
    ref = ctypes.byref
    ptr = lambda array: array.ctypes.data_as(ctypes.c_void_p)
    status = _load_intel_library().disort_(
        ref(nlyr_v), ptr(dtauc_a), ptr(ssalb_a), ref(nmom_v), ptr(pmom_a), ptr(temper_a),
        ref(wvnmlo_v), ref(wvnmhi_v), ref(usrtau_v), ref(ntau_v), ptr(utau_a), ref(nstr_v),
        ref(usrang_v), ref(numu_v), ptr(umu_a), ref(nphi_v), ptr(phi_a), ref(ibcnd_v),
        ref(fbeam_v), ref(umu0_v), ref(phi0_v), ref(fisot_v), ref(lamber_v), ref(albedo_v),
        ref(btemp_v), ref(ttemp_v), ref(temis_v), ref(plank_v), ref(onlyfl_v), ref(accur_v),
        ptr(prnt), ref(maxcly_v), ref(maxulv_v), ref(maxumu_v), ref(maxphi_v), ref(maxmom_v),
        ptr(out["rfldir"]), ptr(out["rfldn"]), ptr(out["flup"]), ptr(out["dfdt"]), ptr(out["uavg"]),
        ptr(out["uu"]), ptr(out["albmed"]), ptr(out["trnmed"]),
    )
    if status:
        raise DisortError(f"DISORT (original Intel kernel) returned error code {status}")
    return out


def _getmom_intel(iphas: int, asymmetry: float, highest_order: int) -> np.ndarray:
    """GETMOM of the original Intel kernel, called as DISORTwrapper.c calls it (float32 GG)."""

    if abs(np.float32(asymmetry)) >= 1.0:
        raise ValueError("abs GG must be less than 1.0")
    values = idl_placed(np.zeros(highest_order + 1, dtype=np.float32))
    iphas_v, gg_v, nmom_v = ctypes.c_int32(int(iphas)), ctypes.c_float(float(np.float32(asymmetry))), ctypes.c_int32(int(highest_order))
    status = _load_intel_library().getmom_(
        ctypes.byref(iphas_v), ctypes.byref(gg_v), ctypes.byref(nmom_v), values.ctypes.data_as(ctypes.c_void_p)
    )
    if status:
        raise DisortError(f"GETMOM (original Intel kernel) returned error code {status}")
    return np.array(values)


def _plkavg_intel(wavenumber_low: float, wavenumber_high: float, temperature: float) -> float:
    """PLKAVG of the original Intel kernel (arguments converted to float32, as IDL_CvtFlt does)."""

    wnlo_v = ctypes.c_float(float(np.float32(wavenumber_low)))
    wnhi_v = ctypes.c_float(float(np.float32(wavenumber_high)))
    temp_v = ctypes.c_float(float(np.float32(temperature)))
    return float(_load_intel_library().plkavg_(ctypes.byref(wnlo_v), ctypes.byref(wnhi_v), ctypes.byref(temp_v)))


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
    """Call DISORT at the selected precision (ORAC_LUT_DISORT_PRECISION, default double).

    Single precision uses the kernel selected by ORAC_LUT_DISORT_SP_KERNEL
    (``intel_dlm``, the original V21 kernel, unless ``gfortran`` is chosen).
    """

    if disort_precision() == "double":
        return _disort_double(*args, **kwargs)
    if single_kernel() == "intel_dlm":
        return _disort_single_intel(*args, **kwargs)
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
