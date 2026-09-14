import inspect

from oraclut.optics import legacy_stg_optics


def test_legacy_optics_has_a_real_legacy_mie_backend():
    parameters = inspect.signature(legacy_stg_optics).parameters
    assert {"effective_radii_microns", "wavelengths_microns", "refractive_indices"} <= set(parameters)
