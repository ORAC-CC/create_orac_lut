"""Regression tests: reproduce the captured V21 IDL single-precision DISORT calls bitwise.

Reference: tests/data/v21_idl_single/v21_idl_single_calls.npz,
the 90 DISORT calls captured inside the V21 IDL generator (four compact Terra
MODIS p240 cases; diffuse, direct-beam and emission calls) with the outputs of
the original Intel DLM create_orac_lut/disort2/src/DISORTtoIDL.so.

Single precision is selected explicitly in every test (the production default
is double).
"""

from pathlib import Path
import os

import numpy as np
import pytest

from oraclut.radiative_transfer import legacy_disort

ROOT = Path(__file__).resolve().parents[1]
REFERENCE = ROOT / "tests/data/v21_idl_single/v21_idl_single_calls.npz"
OUTPUTS = ("rfldir", "rfldn", "flup", "dfdt", "uavg", "uu", "albmed", "trnmed")

pytestmark = pytest.mark.skipif(not REFERENCE.is_file(), reason="V21 IDL capture reference not present")


def _intel_available():
   return legacy_disort.INTEL_DLM_LIBRARY.is_file() and legacy_disort.intel_runtime_directory().is_dir()


needs_intel = pytest.mark.skipif(not _intel_available(), reason="original Intel DLM or Intel 2022 runtime not available")


@pytest.fixture(scope="module")
def reference():
   return np.load(REFERENCE)


@pytest.fixture
def single_intel(monkeypatch):
   monkeypatch.setenv("ORAC_LUT_DISORT_PRECISION", "single")
   monkeypatch.setenv("ORAC_LUT_DISORT_SP_KERNEL", "intel_dlm")


def _call(reference, i):
   k = f"c{i:03d}_"
   n = int(reference[k + "nlyr"])
   plank = bool(reference[k + "plank"])
   return dict(
      dtauc=reference[k + "dtauc"].copy(), single_scatter_albedo=reference[k + "ssalb"].copy(),
      phase_moments=np.asfortranarray(reference[k + "pmom"]), utau=reference[k + "utau"].copy(),
      umu=reference[k + "umu"].copy(), phi=reference[k + "phi"].copy(), fbeam=float(reference[k + "fbeam"]),
      umu0=float(reference[k + "umu0"]), fisot=float(reference[k + "fisot"]), nstreams=int(reference[k + "nstr_v21"]),
      plank=plank, wavenumber_low=float(reference[k + "wvnmlo"]), wavenumber_high=float(reference[k + "wvnmhi"]),
      temperature=reference[k + "temper"].copy() if plank else reference[k + "temper"][:n].copy())


def _bitwise(result, reference, i):
   k = f"c{i:03d}_"
   return [name for name in OUTPUTS
           if not np.array_equal(np.asarray(result[name], np.float32).view(np.uint32), reference[k + "out_" + name].view(np.uint32))]


def _kind(reference, i):
   k = f"c{i:03d}_"
   return "emission" if reference[k + "plank"] else ("direct" if reference[k + "fbeam"] > 0 else "diffuse")


def test_reference_configuration(reference):
   n = int(reference["ncalls"])
   assert n == 90
   cases = [str(reference[f"c{i:03d}_case"]) for i in range(n)]
   assert [cases.count(c) for c in "ABCD"] == [27, 27, 18, 18]
   kinds = [_kind(reference, i) for i in range(n)]
   assert kinds.count("emission") == 9 and kinds.count("direct") == 45 and kinds.count("diffuse") == 36
   for i in range(n):
      k = f"c{i:03d}_"
      # The fixed V21 configuration: 60 streams, 1000 moments, float32 inputs.
      assert int(reference[k + "nstr_v21"]) == 60 and int(reference[k + "nmom"]) == 999
      assert reference[k + "pmom"].shape == (1000, int(reference[k + "nlyr"]))
      assert reference[k + "dtauc"].dtype == np.float32 and reference[k + "pmom"].dtype == np.float32
      assert (int(reference[k + "usrtau"]), int(reference[k + "usrang"]), int(reference[k + "ibcnd"]),
              int(reference[k + "lamber"]), int(reference[k + "onlyfl"])) == (1, 1, 0, 1, 0)
      assert float(reference[k + "accur"]) == float(np.float32(1e-8))
      assert float(reference[k + "phi0"]) == float(reference[k + "albedo"]) == float(reference[k + "btemp"]) == 0.0
      assert float(reference[k + "ttemp"]) == float(reference[k + "temis"]) == 0.0
      if not reference[k + "plank"]:
         assert not np.any(reference[k + "temper"])


@pytest.mark.parametrize("size", [1, 2, 3, 7, 11, 20, 45, 1000, 2001])
def test_idl_placement(size):
   for order in ("C", "F"):
      values = np.arange(size, dtype=np.float32).reshape(size, 1)
      placed = legacy_disort.idl_placed(values, order)
      assert placed.ctypes.data % legacy_disort.IDL_ALIGNMENT == legacy_disort.IDL_OFFSET
      assert placed.dtype == np.float32 and placed.shape == values.shape and np.array_equal(placed, values)
      assert placed.flags["F_CONTIGUOUS" if order == "F" else "C_CONTIGUOUS"]


def test_single_precision_selection_is_explicit(monkeypatch):
   monkeypatch.setenv("ORAC_LUT_DISORT_PRECISION", "single")
   monkeypatch.setenv("ORAC_LUT_DISORT_SP_KERNEL", "gfortran")
   assert legacy_disort.single_kernel() == "gfortran"
   monkeypatch.setenv("ORAC_LUT_DISORT_SP_KERNEL", "intel_dlm")
   assert legacy_disort.single_kernel() == "intel_dlm"
   monkeypatch.setenv("ORAC_LUT_DISORT_SP_KERNEL", "ifort")
   with pytest.raises(ValueError):
      legacy_disort.single_kernel()


@needs_intel
def test_intel_kernel_reproduces_all_v21_calls(reference, single_intel):
   failures = {}
   for i in range(int(reference["ncalls"])):
      differing = _bitwise(legacy_disort.disort(**_call(reference, i)), reference, i)
      if differing:
         failures[i] = differing
   assert not failures, f"calls not bitwise equal to V21 IDL: {failures}"


@needs_intel
def test_intel_kernel_repeated_and_interleaved_calls(reference, single_intel):
   # The failing p240 node (case A capture 24: ch6, tau 256, r_e 29, SZA 80),
   # an emission call and a single-layer call, repeated in different orders so
   # that the dimensions change between calls.
   kinds = [_kind(reference, i) for i in range(int(reference["ncalls"]))]
   emission = kinds.index("emission")
   single_layer = next(i for i in range(int(reference["ncalls"])) if int(reference[f"c{i:03d}_nlyr"]) == 1)
   for i in [23, emission, 23, single_layer, emission, 23, 0, single_layer]:
      assert not _bitwise(legacy_disort.disort(**_call(reference, i)), reference, i), i


@needs_intel
def test_intel_kernel_modifies_arrays_in_place_like_the_idl_dlm(reference, single_intel):
   # The DLM works on the IDL variables' data: DISORT sets UTAU within 1e-4
   # of the column optical depth to TAUC(NLYR), and call_disort.pro passes
   # UTAU through, so IDL's direct-beam call 10 (case A) starts from the UTAU
   # left by diffuse call 9.  SSALB = 1 becomes 1 - DITHER in place.
   assert _kind(reference, 9) == "diffuse" and _kind(reference, 10) == "direct"
   args = _call(reference, 9)
   assert not np.array_equal(args["utau"], reference["c010_utau"])
   legacy_disort.disort(**args)
   assert np.array_equal(args["utau"].view(np.uint32), reference["c010_utau"].view(np.uint32))
   conservative = next(i for i in range(int(reference["ncalls"])) if np.any(reference[f"c{i:03d}_ssalb"] == 1.0))
   args = _call(reference, conservative)
   legacy_disort.disort(**args)
   ones = reference[f"c{conservative:03d}_ssalb"] == 1.0
   assert np.all(args["single_scatter_albedo"][ones] < 1.0) and np.all(args["single_scatter_albedo"][ones] > 0.9999)


def test_call_disort_mirror_copies_layer_arrays_like_idl(reference, single_intel):
   # call_disort.pro copies DTAUC, SSALB and PMOM (compressed, or dtauc1 =
   # dtauc for a Planck call), so the caller's SSALB is never dithered.
   if not _intel_available():
      pytest.skip("original Intel DLM or Intel 2022 runtime not available")
   from types import SimpleNamespace
   from oraclut.idl_mirror.call_disort import call_disort
   conservative = next(i for i in range(int(reference["ncalls"])) if np.any(reference[f"c{i:03d}_ssalb"] == 1.0))
   args = _call(reference, conservative)
   nlyr = args["dtauc"].size
   d = SimpleNamespace(nlyr=nlyr, nmom=999, ntau=2, numu=args["umu"].size, nphi=args["phi"].size, nstr=60)
   ssalb = args["single_scatter_albedo"].copy()
   call_disort(d, args["dtauc"], ssalb, args["phase_moments"], args["utau"], args["umu"], args["phi"],
               args["fbeam"], args["umu0"], args["fisot"])
   assert np.array_equal(ssalb, args["single_scatter_albedo"])


def test_gfortran_kernel_remains_selectable(reference, monkeypatch):
   monkeypatch.setenv("ORAC_LUT_DISORT_PRECISION", "single")
   monkeypatch.setenv("ORAC_LUT_DISORT_SP_KERNEL", "gfortran")
   if not legacy_disort._library_path().is_file():
      pytest.skip("gfortran DISORT library not built")
   for i in (0, 23):
      args = _call(reference, i)
      direct = legacy_disort._disort_single(**_call(reference, i))
      selected = legacy_disort.disort(**args)
      for name in OUTPUTS:
         assert np.array_equal(selected[name], direct[name])


# ---------------------------------------------------------------------------
# Viewing geometry: IDL cos(angle*!dtor) (create_orac_cloud_lut.pro)
# ---------------------------------------------------------------------------

def _case_geometry(case):
   from oraclut.idl_mirror import load_inststr, load_lutstr
   inststr = load_inststr(ROOT / "create_orac_lut/input_files/inst/terra_modis_v1.inst")
   return load_lutstr(ROOT / f"tests/data/v21_idl_single/cases/case_{case}.lut", inststr.max_sat_zenith)


def _idl_umu(saz):
   # create_orac_luts.py, as create_orac_cloud_lut.pro
   from oraclut.idl_math import idl_cos_deg
   f32 = np.float32
   n = saz.size
   tmpsat = saz.copy()
   tmpsat[tmpsat == 90.0] = f32(89.99)
   umu = np.zeros(2 * n, f32)
   umu[0:n] = -f32(1.0) * idl_cos_deg(tmpsat)
   umu[n:2 * n] = -f32(1.0) * umu[n - 1::-1]
   bd = np.abs(umu) == 1.0
   umu[bd] = f32(0.99999) * umu[bd] / np.abs(umu[bd])
   return umu


def test_idl_dtor_is_the_float_system_variable():
   from oraclut.idl_math import IDL_DTOR
   assert IDL_DTOR.dtype == np.float32 and int(IDL_DTOR.view(np.uint32)) == 0x3C8EFA35


@pytest.mark.parametrize("case", "ABCD")
def test_umu_and_umu0_reproduce_v21_captures(reference, case):
   from oraclut.idl_math import idl_cos_deg
   lutstr = _case_geometry(case)
   umu = _idl_umu(lutstr.saz)
   assert lutstr.saz[0] == 0.0 and umu[0] == -np.float32(0.99999) and umu[-1] == np.float32(0.99999)   # nadir rule
   assert lutstr.saz[-1] == np.float32(65.0)                                                # instrument maximum
   nominal = idl_cos_deg(np.float32(50.0))
   calls = [i for i in range(int(reference["ncalls"])) if str(reference[f"c{i:03d}_case"]) == case]
   for i in calls:
      k = f"c{i:03d}_"
      assert np.array_equal(umu.view(np.uint32), reference[k + "umu"].view(np.uint32)), i
      if reference[k + "fbeam"] == 0.0:
         assert nominal.view(np.uint32) == reference[k + "umu0"].view(np.uint32), i
      else:
         assert reference[k + "umu0"] in idl_cos_deg(lutstr.soz), i


def test_idl_cos_matches_float32_rounded_double_cosine_on_production_grids():
   # glibc cosf and the correctly rounded cosine agree on every production
   # angle, so the choice between them does not affect the LUT grids.
   from oraclut.idl_math import IDL_DTOR, idl_cos_deg
   from oraclut.idl_mirror import load_inststr, load_lutstr
   inststr = load_inststr(ROOT / "create_orac_lut/input_files/inst/terra_modis_v1.inst")
   for name in ("liquid-water-cloud.lut", "ice-cloud.lut"):
      lutstr = load_lutstr(ROOT / "create_orac_lut/input_files/lut" / name, inststr.max_sat_zenith)
      for angles in (lutstr.saz, lutstr.soz, np.float32([50.0])):
         angles = np.where(angles == 90.0, np.float32(89.99), angles).astype(np.float32)
         rounded = np.cos((angles * IDL_DTOR).astype(np.float64)).astype(np.float32)
         assert np.array_equal(idl_cos_deg(angles), rounded), name


# ---------------------------------------------------------------------------
# IDL elementary arithmetic (oraclut.idl_math)
# ---------------------------------------------------------------------------

def test_idl_integer_power_is_repeated_squaring():
   from oraclut.idl_math import idl_ipow
   x = np.float32([0.4123456, 2.1305, 11.03])
   assert np.array_equal(idl_ipow(x, 2), x * x)
   assert np.array_equal(idl_ipow(x, 4), (x * x) * (x * x))
   assert idl_ipow(x, 4).dtype == np.float32


def test_idl_exp_uses_the_c_library():
   # IDL FLOAT EXP is the DOUBLE exp rounded to FLOAT (reproduces the V21 T_dv)
   import math
   from oraclut.idl_math import idl_exp
   x = np.linspace(-20, 5, 2001, dtype=np.float32)
   assert np.array_equal(idl_exp(x), np.array([math.exp(float(v)) for v in x]).astype(np.float32))
   y = np.linspace(-20, 5, 2001)
   assert np.array_equal(idl_exp(y), np.array([math.exp(v) for v in y]))
   with pytest.raises(TypeError):
      idl_exp(np.arange(3))


# ---------------------------------------------------------------------------
# End-to-end compact case: the corrected generator against the V21 IDL products
# ---------------------------------------------------------------------------

COMPACT = ROOT / "tests/data/v21_idl_single/idl_compact"
SINGLE_CHANNEL_FILL = {"extinction_coefficient", "extinction_coefficient_ratio", "single_scatter_albedo", "asymmetry_parameter"}


@pytest.mark.slow
@needs_intel
def test_compact_case_a_reproduces_v21_idl_bitwise(tmp_path, monkeypatch):
   from netCDF4 import Dataset
   from scipy.io import readsav
   import create_orac_luts
   monkeypatch.setenv("ORAC_LUT_DISORT_PRECISION", "single")
   monkeypatch.setenv("ORAC_LUT_DISORT_SP_KERNEL", "intel_dlm")
   in_path = tmp_path / "input_files"
   (in_path / "lut").mkdir(parents=True)
   for entry in (ROOT / "create_orac_lut" / "input_files").iterdir():
      if entry.name != "lut":
         (in_path / entry.name).symlink_to(entry)
   (in_path / "lut" / "case_A.lut").write_bytes((ROOT / "tests/data/v21_idl_single/cases/case_A.lut").read_bytes())
   create_orac_luts.create_orac_cloud_lut(in_path, "terra_modis_v1.inst", "liquid-water_240.mm", "case_A.lut", tmp_path / "luts", "2",
                                          channelid=[6], gas=0, no_rayleigh=0, srf_quad=1, reuse_scat=0, scat_only=0,
                                          tmatrix_path=None, version=22, driver=None, nstreams=60, nmom=1000, work_path=tmp_path / "work")
   # scattering properties: bitwise against the IDL scatfile.sav
   idl = readsav(str(COMPACT / "case_A" / "scatfile.sav"))
   python = np.load(tmp_path / "work" / "scatfile.npz")
   for name in ("bext550", "w550", "g550", "phs550", "amom550", "bextrat", "bext", "w", "g", "vavg", "phs", "amom"):
      assert np.array_equal(np.asarray(python[name]).T.ravel().view(np.uint32), np.asarray(idl[name], np.float32).ravel().view(np.uint32)), name
   # every NetCDF variable bitwise, except the single-channel IDL fill artefact
   with Dataset(COMPACT / "case_A" / "terra_modis_m_liquid-water_a01_p240_v21.nc") as r, \
        Dataset(tmp_path / "luts" / "terra_modis_m_liquid-water_a01_p240_v22.nc") as c:
      for d in (r, c):
         d.set_auto_maskandscale(False)
         d.set_auto_chartostring(False)
      assert list(r.variables) == list(c.variables)
      for name in r.variables:
         a, b = np.asarray(r[name][:]), np.asarray(c[name][:])
         if name in SINGLE_CHANNEL_FILL:
            assert np.all(a == np.float32(9.96921e36)) and a.shape == b.shape     # IDL left NetCDF fill for one channel
            continue
         assert a.dtype == b.dtype and a.shape == b.shape and a.tobytes() == b.tobytes(), name


V21_LUTS = Path("/network/group/aopp/eodg/RGG004_GRAINGER_ORACFILE/ORAC_LUTS")


@pytest.mark.skipif(not (V21_LUTS / "terra_modis_m_liquid-water_a01_p240_v21.nc").is_file(), reason="production V21 LUTs not available")
def test_solar_constant_f0_reproduces_the_production_v21_luts():
   from netCDF4 import Dataset
   from oraclut.idl_mirror import load_inststr, load_srfstrarr
   inputs = ROOT / "create_orac_lut" / "input_files"
   inststr = load_inststr(inputs / "inst" / "terra_modis_v1.inst", requestedchannelid=list(range(1, 37)))
   srfstrarr, _ = load_srfstrarr(inststr, inputs / "sun" / "Gueymard2018.sssi", 1, inputs)
   f0 = np.array([s.f0 for s in srfstrarr], np.float32)[np.asarray(inststr.solar_channel_flag, bool)]
   for product in ("liquid-water_a01_p240", "liquid-water_a01_pold", "water-ice_a01_pagg"):
      with Dataset(V21_LUTS / f"terra_modis_m_{product}_v21.nc") as d:
         d.set_auto_maskandscale(False)
         assert np.array_equal(np.asarray(d["F0"][:]).view(np.uint32), f0.view(np.uint32)), product
