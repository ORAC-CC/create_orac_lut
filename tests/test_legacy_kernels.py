import numpy as np

from oraclut.optics.legacy import legendre_moments
from oraclut.optics.legacy_mie import mie_single_batch
from oraclut.pipeline import read_channels
from oraclut.config import read_instrument
from oraclut.radiative_transfer import call_disort, getmom, plkavg


ROOT = __import__("pathlib").Path(__file__).parents[1]


def test_legendre_moments_preserve_constant_phase_normalisation():
    qv, qw = np.polynomial.legendre.leggauss(8)
    coefficients = legendre_moments(qv, qw, np.ones(8))
    np.testing.assert_allclose(coefficients[:3], [1.0, 0.0, 0.0], atol=1e-12)


def test_preserved_mie_kernel_returns_finite_batch():
    result = mie_single_batch(
        np.asarray([1.0, 2.0], dtype=np.float64),
        1.335 - 2.4e-9j,
        np.linspace(-0.9, 0.9, 8, dtype=np.float64),
    )
    assert result["qext"].shape == (2,)
    assert np.isfinite(result["qext"]).all()
    assert np.isfinite(result["f11"]).all()


def test_production_disort_kernel_smoke():
    phase = getmom(2, 0.0, 60)[:, None]
    result = call_disort(
        np.asarray([0.1], dtype=np.float32),
        np.asarray([1.0], dtype=np.float32),
        np.asfortranarray(phase),
        np.asarray([0.0, 0.1], dtype=np.float32),
        np.asarray([-0.9, 0.9], dtype=np.float32),
        np.asarray([0.0, 180.0], dtype=np.float32),
        0.0, 0.5, 100.0,
        nstreams=60,
    )
    assert result["uu"].shape == (2, 2, 2)
    assert np.isfinite(result["uu"]).all()
    assert plkavg(1000.0, 1010.0, 250.0) > 0.0


def test_channel_reader_matches_current_seviri_channel_one():
    instrument = read_instrument(ROOT / "create_orac_lut/input_files/inst/meteosat-10_seviri_v1.inst")
    channel = read_channels(ROOT, instrument, (1,), srf_quad=1)[0]
    assert channel.solar
    assert not channel.thermal
    np.testing.assert_allclose(channel.wavelength_microns, 0.6381768, rtol=0, atol=1e-7)
    np.testing.assert_allclose(channel.f0, 513.2258, rtol=0, atol=1e-4)
