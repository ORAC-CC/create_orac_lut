"""IDL 8.9 reference and V21 integration tests for the TOTAL mirror."""

from pathlib import Path

import numpy as np

from oraclut.idl_mirror import idl_total, load_inststr, load_srfstrarr


ROOT = Path(__file__).resolve().parents[1]
INPUTS = ROOT / "create_orac_lut" / "input_files"

V21_TERRA_MODIS_CENTRAL_WAVELENGTH = np.asarray([
   0.6462792754173279, 0.8564648628234863, 0.46574240922927856, 0.5537328124046326,
   1.242307424545288, 1.6293855905532837, 2.1142096519470215, 0.41185253858566284,
   0.44211018085479736, 0.4869867265224457, 0.5297242999076843, 0.5468820929527283,
   0.6657798290252686, 0.6769617199897766, 0.7466169595718384, 0.8663466572761536,
   0.9042461514472961, 0.9356905817985535, 0.9362413883209229, 3.7882111072540283,
   3.9921560287475586, 3.9719860553741455, 4.056727886199951, 4.473208427429199,
   4.545375347137451, 1.3824920654296875, 6.752601146697998, 7.332195281982422,
   8.528982162475586, 9.724987983703613, 11.018420219421387, 12.032552719116211,
   13.36523723602295, 13.668963432312012, 13.896821022033691, 14.17526912689209,
], dtype=np.float32)

V21_TERRA_MODIS_CENTRAL_WAVENUMBER = np.asarray([
   15497.0302734375, 11683.8662109375, 21481.541015625, 18065.705078125,
   8050.8759765625, 6138.048828125, 4731.1455078125, 24290.693359375,
   22622.537109375, 20537.48828125, 18880.806640625, 18287.693359375,
   15021.2646484375, 14773.40234375, 13394.6865234375, 11544.0322265625,
   11064.275390625, 10688.123046875, 10690.1875, 2641.8388671875,
   2505.31298828125, 2518.0546875, 2465.431396484375, 2235.85009765625,
   2200.348388671875, 7235.61767578125, 1481.9339599609375, 1364.814453125,
   1173.2637939453125, 1028.7117919921875, 908.2734375, 831.523193359375,
   748.34716796875, 731.7239379882812, 719.7265625, 705.5760498046875,
], dtype=np.float32)


def test_idl_89_reference_precision_and_dtype():
   float_values = np.asarray([1.0e8, 1.0, 1.0, -1.0e8, 3.0], dtype=np.float32)
   double_values = np.asarray([1.0e16, 1.0, 1.0, -1.0e16, 3.0], dtype=np.float64)
   integer_values = np.asarray([1, 2, 3, 4], dtype=np.int16)
   float_total = idl_total(float_values)
   promoted_total = idl_total(float_values, double=True)
   double_total = idl_total(double_values)
   integer_total = idl_total(integer_values)
   assert float_total.dtype == np.dtype(np.float32)
   assert promoted_total.dtype == np.dtype(np.float64)
   assert double_total.dtype == np.dtype(np.float64)
   assert integer_total.dtype == np.dtype(np.float32)
   assert float_total == np.float32(3.0)
   assert promoted_total == np.float64(5.0)
   assert double_total == np.float64(3.0)
   assert integer_total == np.float32(10.0)


def test_idl_89_reference_dimensions_cumulative_and_order():
   values = np.asarray([1.0e8, 1.0, 1.0, -1.0e8, 3.0, 4.0], dtype=np.float32).reshape((2, 3), order="F")
   expected_cumulative = np.asarray(
      [1.0e8, 1.0, 1.0e8, -1.0e8, 1.0e8, -1.0e8], dtype=np.float32
   ).reshape((2, 3), order="F")
   assert idl_total(values) == np.float32(7.0)
   np.testing.assert_array_equal(idl_total(values, dimension=1), np.asarray([1.0e8, -1.0e8, 7.0], np.float32))
   np.testing.assert_array_equal(idl_total(values, dimension=2), np.asarray([1.0e8, -1.0e8], np.float32))
   np.testing.assert_array_equal(idl_total(values, dimension=2, cumulative=True), expected_cumulative)
   np.testing.assert_array_equal(
      idl_total(np.asarray([1.0e8, 1.0, 1.0, -1.0e8, 3.0], np.float32), cumulative=True),
      np.asarray([1.0e8, 1.0e8, 1.0e8, 0.0, 3.0], np.float32),
   )


def test_idl_total_rejects_unsupported_inputs_and_dimensions():
   with np.testing.assert_raises(TypeError):
      idl_total(np.asarray(["not", "numeric"]))
   with np.testing.assert_raises(ValueError):
      idl_total(np.asarray([], dtype=np.float32))
   with np.testing.assert_raises(ValueError):
      idl_total(np.ones((2, 2), dtype=np.float32), dimension=3)
   with np.testing.assert_raises(TypeError):
      idl_total(np.ones((2, 2), dtype=np.float32), dimension=1.0)


def test_terra_modis_centres_are_bitwise_equal_to_v21():
   instrument = load_inststr(INPUTS / "inst" / "terra_modis_v1.inst")
   spectral_response, _ = load_srfstrarr(
      instrument, INPUTS / "sun" / "Gueymard2018.sssi", 1, INPUTS
   )
   wavelength = np.asarray([item.wvl_centre for item in spectral_response], dtype=np.float32)
   wavenumber = np.asarray([item.wvn_centre for item in spectral_response], dtype=np.float32)
   assert np.array_equal(wavelength, V21_TERRA_MODIS_CENTRAL_WAVELENGTH)
   assert np.array_equal(wavenumber, V21_TERRA_MODIS_CENTRAL_WAVENUMBER)


def test_v22_total_call_sites_do_not_use_numpy_sum():
   paths = [
      ROOT / "create_orac_luts.py",
      ROOT / "src" / "oraclut" / "idl_mirror" / "load_srfstrarr.py",
      ROOT / "src" / "oraclut" / "idl_mirror" / "create_bwgp.py",
      ROOT / "src" / "oraclut" / "idl_mirror" / "generate_scattering_properties.py",
      ROOT / "src" / "oraclut" / "idl_mirror" / "dubovik.py",
      ROOT / "src" / "oraclut" / "idl_mirror" / "write_v2_lut.py",
      ROOT / "src" / "oraclut" / "pipeline.py",
      ROOT / "src" / "oraclut" / "optics" / "legacy.py",
   ]
   for path in paths:
      assert "np.sum(" not in path.read_text(), path
