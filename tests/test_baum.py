"""Targeted tests for the historical Baum ice optical-property pathway."""

from pathlib import Path

import numpy as np
import pytest

from oraclut.idl_mirror import BaumTable, generate_scattering_properties, load_inststr, load_lutstr, load_mmdat, load_srfstrarr


ROOT = Path(__file__).parents[1]
INPUTS = ROOT / "create_orac_lut" / "input_files"
BAUM_FILES = {
   "agg": INPUTS / "baum" / "AggregateSolidColumns_SeverelyRough_AllWavelengths_FullPhaseMatrix.nc",
   "ghm": INPUTS / "baum" / "GeneralHabitMixture_SeverelyRough_AllWavelengths_FullPhaseMatrix.nc",
   "src": INPUTS / "baum" / "SolidColumns_SeverelyRough_AllWavelengths_FullPhaseMatrix.nc",
}


@pytest.mark.parametrize("shortname", ["agg", "ghm", "src"])
def test_baum_table_loads_and_interpolates_each_habit(shortname):
   table = BaumTable(BAUM_FILES[shortname])
   bext, albedo, asymmetry, phase = table.interpolate(
      wavelengths=[0.55, 10.76829], effective_radii=[1.0, 5.0, 93.0], phase_angles=[0.0, 90.0, 180.0]
   )
   assert table.effective_radii.tolist() == list(np.arange(5.0, 60.0 + 2.5, 2.5))
   assert bext.shape == albedo.shape == asymmetry.shape == (2, 3)
   assert phase.shape == (3, 2, 3)
   assert np.all(np.isfinite(bext)) and np.all(np.isfinite(albedo))
   assert np.all(np.isfinite(asymmetry)) and np.all(np.isfinite(phase))
   assert np.all((albedo >= 0.0) & (albedo <= 1.0))


def test_baum_habits_are_distinct():
   values = [BaumTable(path).interpolate([0.55], [20.0], [90.0])[0][0, 0] for path in BAUM_FILES.values()]
   assert len(set(np.round(values, 7))) == 3


@pytest.mark.parametrize("shortname", ["agg", "ghm", "src"])
def test_baum_habit_reaches_scattering_interface(shortname):
   mmstr = load_mmdat(INPUTS / "microphysics" / f"water-ice_{shortname}.mm", INPUTS)
   assert mmstr.comptype == ["baum"]
   inststr = load_inststr(INPUTS / "inst" / "earthcare_msi_v1.inst", requestedchannelid=[1])
   lutstr = load_lutstr(INPUTS / "lut" / "ice-cloud_test.lut", inststr.max_sat_zenith)
   srfstrarr, nwvl_max = load_srfstrarr(inststr, INPUTS / "sun" / "Gueymard2018.sssi", 1, INPUTS)
   result = generate_scattering_properties(srfstrarr, nwvl_max, inststr, mmstr, lutstr, 24)
   _, bext550, w550, g550, _, _, _, bext, w, g, _, phs, _ = result
   assert bext550.shape == w550.shape == g550.shape == (lutstr.efr_n,)
   assert bext.shape == w.shape == g.shape == (1, 1, lutstr.efr_n)
   assert phs.shape == (24, 1, 1, lutstr.efr_n)
   assert np.all(np.isfinite(bext)) and np.all(np.isfinite(w))
   assert np.all(np.isfinite(g)) and np.all(np.isfinite(phs))
   _, weights = np.polynomial.legendre.leggauss(24)
   assert np.allclose(np.sum(phs * weights[:, None, None, None], axis=0) / 2.0, 1.0, atol=2e-5)
