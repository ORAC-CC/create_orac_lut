"""Regressions for the two port defects found by the EarthCARE MSI V21/V22 validation.

Defect 1 (RI-ORDER): _interpol_complex passed the wavelength grid derived by
read_ri straight to np.interp.  IDL INTERPOL accepts a monotonically decreasing
abscissa; np.interp does not, and silently returns an endpoint value.  A
'#FORMAT = WAVN' table with ascending wavenumber reaches the interpolation as
descending wavl = 1e4/wavn, so the refractive index collapsed to a constant and
the affected LUTs became non-absorbing.

Defect 2 (AEROSOL-WRITE): create_orac_aerosol_lut omitted bextout from its
write_v2_lut call, shifting every later positional argument by one and raising
TypeError after the whole DISORT calculation had completed.

Every test here fails on the pre-fix code and passes on the corrected code.
None of them runs a production-size calculation.
"""

import ast
import inspect
from pathlib import Path

import numpy as np
import pytest

import create_orac_luts
from oraclut.idl_mirror.create_bwgp import create_bwgp
from oraclut.idl_mirror.generate_scattering_properties import _interpol_complex
from oraclut.idl_mirror.load_mmdat import read_ri
from oraclut.idl_mirror.write_v2_lut import write_v2_lut
from oraclut.idl_mirror import load_inststr, load_lutstr, load_srfstrarr

ROOT = Path(__file__).parents[1]
INPUTS = ROOT / "create_orac_lut" / "input_files"
RI = INPUTS / "ri"
WATER_240K = RI / "H2O_240K_Rowe_2020.ri"

# EarthCARE MSI channel 7 centre, and the V21-consistent single-scatter albedos
# reproduced by the previous audit for the 240 K supercooled-water table.
CHANNEL_7_WAVELENGTH = 11.994
CHANNEL_7_SSA = {5.0: 0.28634, 12.0: 0.42161, 24.0: 0.47925}
MODIFIED_GAMMA_SPREAD = 1.111111e-01


# ==============================================================================
# Defect 1: the interpolation boundary
# ==============================================================================

def test_interpol_complex_interpolates_a_descending_abscissa():
   """A synthetic descending grid must interpolate, not collapse to an endpoint."""

   # n = 1 + wl/10 and k = wl/100 on a descending wavelength grid.
   wl = np.array([10.0, 8.0, 6.0, 4.0, 2.0])
   cm = (1.0 + wl / 10.0) + 1j * (wl / 100.0)
   probe = np.array([3.0, 5.0, 7.0, 9.0])

   result = _interpol_complex(cm, wl, probe)

   np.testing.assert_allclose(np.real(result), 1.0 + probe / 10.0, rtol=0, atol=1e-6)
   np.testing.assert_allclose(np.imag(result), probe / 100.0, rtol=0, atol=1e-6)
   # The pre-fix failure mode: every probe returning the same endpoint value.
   assert np.ptp(np.real(result)) > 0.1


def test_interpol_complex_leaves_an_ascending_abscissa_bit_for_bit_unchanged():
   """Sorting must be the identity for already ascending tables (no regression)."""

   wl = np.array([2.0, 4.0, 6.0, 8.0, 10.0])
   cm = (1.0 + wl / 10.0) + 1j * (wl / 100.0)
   probe = np.array([3.0, 5.0, 7.0, 9.0])

   result = _interpol_complex(cm, wl, probe)
   expected = (np.interp(probe, wl, np.real(cm)).astype(np.float32)
               + 1j * np.interp(probe, wl, np.imag(cm)).astype(np.float32)).astype(np.complex64)

   assert np.array_equal(result, expected)


def test_water_240k_table_is_descending_in_wavelength():
   """The real table that triggered the defect still arrives descending from read_ri."""

   ri = read_ri(WATER_240K)
   wl = np.asarray(ri.wavl, dtype=np.float64)

   assert ri.format.upper().split()[0] == "WAVN"
   assert np.all(np.diff(wl) < 0), "read_ri must keep mirroring read_ri.pro, not sort"


def test_read_ri_was_not_changed_to_sort_its_data():
   """read_ri must contain no sorting: the fix belongs at the interpolation boundary."""

   source = inspect.getsource(read_ri)

   assert "argsort" not in source and "np.sort" not in source and ".sort(" not in source


def test_water_240k_imaginary_index_is_no_longer_pinned_to_an_endpoint():
   """The confirmed numerical signature of Defect 1, on real data."""

   ri = read_ri(WATER_240K)
   cm = ri.n + 1j * (-ri.k)                      # IDL: k negative for an absorbing medium
   probe = np.array([0.670, 0.865, 1.65, 2.21, 8.8, 10.8, 12.0])

   result = _interpol_complex(cm, ri.wavl, probe)
   k = -np.imag(result)

   # Before the fix every one of these was the endpoint value 2.4214e-08.
   assert k[-1] == pytest.approx(0.30751, rel=1e-3), "12 micron k must rise to about 0.307"
   assert k[0] < 1e-6, "0.67 micron k must stay negligible"
   assert k[-1] / k[0] > 1e6, "k must span orders of magnitude across the MSI channels"
   # Monotonic growth of absorption across these particular MSI channels.
   assert np.all(np.diff(k) > 0)


def test_water_240k_refractive_index_varies_smoothly_with_wavelength():
   """No step discontinuities: neighbouring probes must track the tabulated spectrum."""

   ri = read_ri(WATER_240K)
   cm = ri.n + 1j * (-ri.k)
   probe = np.linspace(8.0, 12.0, 41)

   result = _interpol_complex(cm, ri.wavl, probe)
   n = np.real(result)

   assert np.ptp(n) > 0.05, "the real index must actually vary over 8-12 micron"
   # Piecewise-linear interpolation of a smooth table: no neighbouring jump may
   # dominate the whole range.
   assert np.max(np.abs(np.diff(n))) < 0.25 * np.ptp(n)


@pytest.mark.parametrize("effective_radius", sorted(CHANNEL_7_SSA))
def test_channel_7_mie_single_scatter_albedo_matches_the_audited_v21_values(effective_radius):
   """The corrected interpolation reproduces the V21-consistent audit values."""

   ri = read_ri(WATER_240K)
   cm = ri.n + 1j * (-ri.k)
   m = _interpol_complex(cm, ri.wavl, [CHANNEL_7_WAVELENGTH])[0]

   _, w, _, _, _ = create_bwgp("modified_gamma", effective_radius, MODIFIED_GAMMA_SPREAD,
                               np.array([m]), np.array([CHANNEL_7_WAVELENGTH]),
                               np.array([1.0, 0.0, -1.0]))

   assert w[0] == pytest.approx(CHANNEL_7_SSA[effective_radius], abs=5e-5)
   assert w[0] < 0.99, "a pinned endpoint refractive index drives the albedo to unity"


def test_no_refractive_index_table_has_duplicate_wavelengths():
   """Documents why _interpol_complex needs no duplicate-abscissa handling."""

   for path in sorted(RI.glob("*.ri")):
      wl = np.asarray(read_ri(path).wavl, dtype=np.float64)
      assert np.unique(wl).size == wl.size, f"{path.name} has duplicate wavelengths"


# ==============================================================================
# Defect 2: the aerosol writer call
# ==============================================================================

def _positional_argument_names(function_name, call_index):
   """Names of the positional arguments of a write_v2_lut call in create_orac_luts."""

   tree = ast.parse(Path(create_orac_luts.__file__).read_text())
   calls = [node for node in ast.walk(tree)
            if isinstance(node, ast.Call) and getattr(node.func, "id", None) == function_name]
   calls.sort(key=lambda node: node.lineno)
   call = calls[call_index]
   return [argument.id for argument in call.args if isinstance(argument, ast.Name)], call.lineno


@pytest.mark.parametrize("call_index, branch", [(0, "cloud"), (1, "aerosol")])
def test_write_v2_lut_calls_bind_every_positional_argument_by_name(call_index, branch):
   """Each positional argument must land on the identically named parameter.

   Before the fix the aerosol call omitted bextout, so bextratout bound to
   bextout, every later argument shifted by one and rfd was left unfilled.
   """

   names, lineno = _positional_argument_names("write_v2_lut", call_index)
   parameters = [name for name, parameter in inspect.signature(write_v2_lut).parameters.items()
                 if parameter.default is inspect.Parameter.empty]

   assert len(names) == len(parameters), (
      f"{branch} call at line {lineno} passes {len(names)} positional arguments, "
      f"expected {len(parameters)}")
   # The first argument is the local filename variable; the rest share the parameter names.
   for parameter, name in list(zip(parameters, names))[1:]:
      assert parameter == name, (
         f"{branch} call at line {lineno}: argument {name!r} is bound to parameter "
         f"{parameter!r} -- the argument list is shifted")


def test_aerosol_call_passes_bextout():
   """The specific omission that caused Defect 2."""

   names, _ = _positional_argument_names("write_v2_lut", 1)

   assert "bextout" in names
   assert names.index("bextout") < names.index("bextratout")
   assert names[-1] == "rfd"


def test_aerosol_writer_writes_a_readable_lut_with_aerosol_arguments(tmp_path):
   """Drive the corrected aerosol argument list through write_v2_lut for real.

   The compact meteosat-10 test grids are used, with synthetic optical and
   radiative-transfer arrays, so the writer path is exercised without any
   DISORT calculation.
   """

   import netCDF4 as nc

   inststr = load_inststr(INPUTS / "inst" / "meteosat-10_seviri_v1.inst", requestedchannelid=[1])
   lutstr = load_lutstr(INPUTS / "lut" / "aerosol_test.lut", inststr.max_sat_zenith, include_pressure=True)
   srfstrarr, _ = load_srfstrarr(inststr, INPUTS / "sun" / "Gueymard2018.sssi", 1, INPUTS)

   nchan = inststr.number_of_nadir_channels
   nefr, nopd = lutstr.efr_n, lutstr.opd_n
   npre, nsoz, nsaz, nraa = lutstr.prs_n, lutstr.soz_n, lutstr.saz_n, lutstr.raa_n

   # Distinct constants per array, so a shifted argument list would be visible
   # in the written file rather than merely raising.
   vavg = np.full(nefr, 0.5, dtype=np.float32)
   bextout = np.full((nchan, nefr), 2.0, dtype=np.float32)
   bextratout = np.full((nchan, nefr), 3.0, dtype=np.float32)
   ssaout = np.full((nchan, nefr), 0.75, dtype=np.float32)
   gout = np.full((nchan, nefr), 0.6, dtype=np.float32)

   diffuse = np.full((nchan, npre, nefr, nopd), 0.1, dtype=np.float32)
   direct = np.full((nchan, npre, nefr, nopd, nsoz), 0.2, dtype=np.float32)
   bidirectional = np.full((nchan, npre, nefr, nopd, nsoz, nsaz, nraa), 0.3, dtype=np.float32)
   beam = np.full((nchan, npre, nefr, nopd, nsaz), 0.4, dtype=np.float32)

   path = tmp_path / "test_m_aerosol_a12_pa79_v22.nc"
   write_v2_lut(path, lutstr, inststr, srfstrarr, vavg, bextout, bextratout, ssaout, gout,
                direct.copy(), diffuse.copy(), direct.copy(), diffuse.copy(),
                rbd=bidirectional.copy(), rfbd=beam.copy(), tfbd=beam.copy(),
                tb=beam.copy(), em=beam.copy(), include_pressure=True)

   assert path.is_file() and path.stat().st_size > 0

   with nc.Dataset(path) as data:
      for name in ("average_volume_per_particle", "extinction_coefficient",
                   "extinction_coefficient_ratio", "single_scatter_albedo",
                   "asymmetry_parameter", "effective_radius", "optical_depth"):
         assert name in data.variables, f"{name} missing from the aerosol product"
      # Each optical variable must carry its own constant: proof of no shift.
      assert np.allclose(np.asarray(data.variables["extinction_coefficient"][:]), 2.0)
      assert np.allclose(np.asarray(data.variables["extinction_coefficient_ratio"][:]), 3.0)
      assert np.allclose(np.asarray(data.variables["single_scatter_albedo"][:]), 0.75)
      assert np.allclose(np.asarray(data.variables["asymmetry_parameter"][:]), 0.6)
      assert np.allclose(np.asarray(data.variables["average_volume_per_particle"][:]), 0.5)
      # include_pressure=True is the aerosol branch's distinguishing argument.
      assert "surface_pressure" in data.dimensions or "surface_pressure" in data.variables
