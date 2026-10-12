from pathlib import Path

import numpy as np

from oraclut.config import (
    read_atmosphere, read_instrument, read_lut_grid, read_microphysics,
    read_srf, read_solar_spectrum, reference_configuration,
)


ROOT = Path(__file__).parents[1]


def test_reference_configuration_uses_current_inputs():
    config = reference_configuration(ROOT)
    assert (config.lut_level, config.revision) == (2, 21)
    assert (config.platform, config.instrument) == ("meteosat-10", "seviri")
    assert config.forward_model == "cloud"
    assert config.channels == tuple(range(1, 12))
    assert config.atmosphere_code == 2
    assert config.gas is False


def test_current_input_readers_match_reference_case():
    config = reference_configuration(ROOT)
    instrument = read_instrument(config.instrument_file)
    grid = read_lut_grid(config.lut_file)
    model = read_microphysics(config.microphysics_file)
    atmosphere = read_atmosphere(config.atmosphere_file, config.atmosphere_code)
    assert instrument.srf_files[1].endswith("ch01.txt")
    np.testing.assert_array_equal(grid.effective_radius, np.arange(1, 40, 2, dtype="f4"))
    assert (model.size_distribution, model.size_parameters) == ("modified_gamma", (12.0, 0.1111111))
    assert model.scattering_code == "mie"
    assert atmosphere.height_km.size == 46
    assert atmosphere.height_km[0] == 100.0
    assert atmosphere.height_km[-1] == 0.0
    srf = read_srf(ROOT / "create_orac_lut/input_files/srf" / instrument.srf_files[1])
    solar = read_solar_spectrum(config.solar_spectrum_file)
    assert srf.coordinate_units == "cm^-1" and srf.coordinate.size == 51
    assert solar.coordinate_units == "microns" and solar.coordinate.size > 1000
