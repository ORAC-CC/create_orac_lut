"""Focused tests for the two remaining V22 campaign ports.

These tests exercise configuration/array transformations only.  They do not
invoke Mie, T-matrix, DISORT, LUT writing, or SLURM.
"""

from pathlib import Path
from types import SimpleNamespace

import numpy as np
import pytest

import create_orac_luts
from oraclut.idl_mirror.load_srfstrarr import (
   _qm2_segmentation_input,
   _qm2_uses_idl_tabulated_sort,
   read_srfstr,
)
from oraclut.idl_mirror.load_inststr import load_inststr
from oraclut.idl_mirror.segment import segment
from oraclut.idl_mirror.write_v2_lut import _channel_optical_property, _solar_uncertainty_variables


ROOT = Path(__file__).parents[1]
INPUTS = ROOT / "create_orac_lut" / "input_files"
SRFS = INPUTS / "srf"


def _dual_view_instrument(view=7):
   instrument = SimpleNamespace(
      view=view,
      number_of_nadir_channels=2,
      number_of_channels=2 if view == 0 else 4,
      channelid=np.asarray([2, 6], dtype=np.int16),
      srf_file=["ch01", "ch06"],
   )
   for name in (
      "solar_channel_flag", "mixed_channel_flag", "thermal_channel_flag",
      "oldf0", "oldf1", "oldnefr", "oldwvn", "oldb1", "oldb2",
      "oldt1", "oldt2", "oldnebt", "snr", "rgu", "rou", "rua",
      "rub", "ruc", "refbt", "nedt",
   ):
      setattr(instrument, name, np.asarray([10, 20], dtype=np.float32))
   return instrument


def test_earthcare_qm2_srf_ordering_keeps_raw_samples_and_sorts_each_integral():
   instrument = SimpleNamespace(platform="earthcare", instrument="msi")
   source = read_srfstr(SRFS / "rtcoef_earthcare_1_msi_srf_ch01.txt")
   assert np.any(np.diff(source.wvl) < 0.0)

   wvl, response = _qm2_segmentation_input(source, instrument)
   np.testing.assert_array_equal(wvl, source.wvl)
   np.testing.assert_array_equal(response, source.srf)

   segmented_wvl, segmented_response = segment(
      wvl, response, np.float32(0.001), minn=12,
      sort=_qm2_uses_idl_tabulated_sort(source, instrument),
   )
   assert segmented_wvl.size == 21
   assert segmented_response.size == segmented_wvl.size
   # IDL retains the local reversals in the segmented output; only each
   # INT_TABULATED call receives a stable sorted copy for integration.
   assert np.any(np.diff(segmented_wvl) < 0.0)
   assert segmented_wvl[0] == source.wvl[0]
   assert segmented_wvl[-1] == source.wvl[-1]
   assert np.all(np.isfinite(segmented_response))


def test_earthcare_descending_qm2_srf_keeps_legacy_descending_orientation():
   instrument = SimpleNamespace(platform="earthcare", instrument="msi")
   source = read_srfstr(SRFS / "rtcoef_earthcare_1_msi_srf_ch06.txt")
   wvl, response = _qm2_segmentation_input(source, instrument)
   np.testing.assert_array_equal(wvl, source.wvl)
   np.testing.assert_array_equal(response, source.srf)
   segmented_wvl, _ = segment(wvl, response, np.float32(0.001), minn=12)
   assert segmented_wvl.size >= 12
   assert np.all(np.diff(segmented_wvl) < 0.0)


def test_earthcare_qm2_load_matches_idl_logged_point_counts():
   instrument = load_inststr(INPUTS / "inst" / "earthcare_msi_v1.inst", requestedchannelid=[1, 6])
   load_srfstrarr = __import__(
      "oraclut.idl_mirror.load_srfstrarr", fromlist=["load_srfstrarr"]
   ).load_srfstrarr
   segmented, nwvl_max = load_srfstrarr(
      instrument, INPUTS / "sun" / "Gueymard2018.sssi", 2, INPUTS
   )
   assert nwvl_max == 21
   assert [srf.nwvl for srf in segmented] == [21, 19]


def test_non_earthcare_nonmonotonic_srf_is_not_silently_sorted():
   instrument = SimpleNamespace(platform="seviri", instrument="seviri")
   source = SimpleNamespace(
      wvl=np.asarray([1.0, 2.0, 1.5, 3.0], dtype=np.float32),
      srf=np.asarray([0.0, 1.0, 0.5, 0.0], dtype=np.float32),
   )
   wvl, response = _qm2_segmentation_input(source, instrument)
   np.testing.assert_array_equal(wvl, source.wvl)
   np.testing.assert_array_equal(response, source.srf)
   with pytest.raises(ValueError, match="strictly ascending"):
      segment(wvl, response, np.float32(0.001), minn=12)


def test_dual_view_replication_matches_legacy_idl_channel_and_axis_rules():
   instrument = _dual_view_instrument(view=7)
   srf = [SimpleNamespace(nwvl=2), SimpleNamespace(nwvl=3)]
   original_srf = srf[:]
   rt_arrays = tuple(np.arange(2 * 3, dtype=np.float32).reshape(2, 3) + i for i in range(9))
   optical_arrays = tuple(np.arange(2 * 4, dtype=np.float32).reshape(2, 4) + i for i in range(4))

   result = create_orac_luts._replicate_dual_view(instrument, srf, rt_arrays, optical_arrays)
   instrument, srf, rt_arrays, optical_arrays = result

   assert instrument.number_of_channels == 4
   assert instrument.channelid.tolist() == [2, 6, 9, 13]
   assert srf == [original_srf[0], original_srf[1], original_srf[0], original_srf[1]]
   for values, original in zip(rt_arrays, tuple(np.arange(2 * 3, dtype=np.float32).reshape(2, 3) + i for i in range(9))):
      np.testing.assert_array_equal(values, np.concatenate((original, original), axis=0))
   for values, original in zip(optical_arrays, tuple(np.arange(2 * 4, dtype=np.float32).reshape(2, 4) + i for i in range(4))):
      np.testing.assert_array_equal(values, original)
   for name in ("solar_channel_flag", "thermal_channel_flag", "oldf0", "nedt"):
      assert getattr(instrument, name).tolist() == [10, 20, 10, 20]


def test_dual_view_working_operator_shapes_reach_four_channels_once():
   instrument = _dual_view_instrument(view=7)
   # Aerosol T_dv working layout includes pressure; cloud does not.
   aerosol_td = np.zeros((2, 3, 2, 2, 4), dtype=np.float32)
   cloud_td = np.zeros((2, 2, 3, 4), dtype=np.float32)
   aerosol_result = create_orac_luts._replicate_dual_view(
      instrument, [], (aerosol_td,), (np.zeros((2, 2), dtype=np.float32),)
   )
   cloud_result = create_orac_luts._replicate_dual_view(
      _dual_view_instrument(view=9), [], (cloud_td,), (np.zeros((2, 2), dtype=np.float32),)
   )
   assert aerosol_result[2][0].shape == (4, 3, 2, 2, 4)
   assert cloud_result[2][0].shape == (4, 2, 3, 4)
   assert aerosol_result[3][0].shape == (2, 2)
   assert cloud_result[3][0].shape == (2, 2)


@pytest.mark.parametrize(
   ("instrument_file", "requested", "expected"),
   [
      ("envisat_aatsr_v1.inst", [2, 6], [2, 6, 9, 13]),
      ("sentinel-3a_slstr_v1.inst", [2, 8], [2, 8, 11, 17]),
      ("sentinel-3b_slstr_v1.inst", [2, 8], [2, 8, 11, 17]),
   ],
)
def test_actual_dual_view_instrument_definitions_replicate_once(instrument_file, requested, expected):
   instrument = load_inststr(INPUTS / "inst" / instrument_file, requestedchannelid=requested)
   n = instrument.number_of_nadir_channels
   rt = np.zeros((n, 2, 3, 4), dtype=np.float32)
   optical = np.zeros((n, 2), dtype=np.float32)
   instrument, _, rt_arrays, optical_arrays = create_orac_luts._replicate_dual_view(
      instrument, [], (rt,), (optical,)
   )
   assert instrument.number_of_nadir_channels == 2
   assert instrument.number_of_channels == 4
   assert instrument.channelid.tolist() == expected
   assert rt_arrays[0].shape == (4, 2, 3, 4)
   assert optical_arrays[0].shape == (2, 2)


def test_dual_view_replication_is_a_noop_for_nadir_only_instruments():
   instrument = _dual_view_instrument(view=0)
   srf = [object(), object()]
   rt_arrays = (np.ones((2, 3)),)
   optical_arrays = (np.ones((2, 4)),)
   result = create_orac_luts._replicate_dual_view(instrument, srf, rt_arrays, optical_arrays)
   assert result[0] is instrument
   assert result[1] is srf
   assert result[2][0] is rt_arrays[0]
   assert result[3][0] is optical_arrays[0]


def test_dual_view_writer_preserves_idl_optical_fill_slots_and_calibration_metadata():
   instrument = load_inststr(INPUTS / "inst" / "envisat_aatsr_v1.inst", requestedchannelid=[2, 6])
   values = np.asarray([[1.0, 2.0], [3.0, 4.0]], dtype=np.float32)
   instrument, _, _, _ = create_orac_luts._replicate_dual_view(instrument, [], (values,), (values,))

   padded = _channel_optical_property(values, instrument)
   assert padded.shape == (4, 2)
   np.testing.assert_array_equal(padded[:2], values)
   np.testing.assert_array_equal(padded[2:], np.full((2, 2), np.float32(9.96921e36)))

   definitions = _solar_uncertainty_variables(instrument, np.flatnonzero(instrument.solar_channel_flag))
   assert [definition[0] for definition in definitions] == ["rua", "rub", "ruc"]


def test_single_view_writer_retains_snr_metadata():
   instrument = load_inststr(INPUTS / "inst" / "meteosat-10_seviri_v1.inst", requestedchannelid=[1])
   definitions = _solar_uncertainty_variables(instrument, np.flatnonzero(instrument.solar_channel_flag))
   assert [definition[0] for definition in definitions] == ["snr"]
