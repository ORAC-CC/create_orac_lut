"""Concurrency and publication tests for the IDL-structured LUT entry point."""

import numpy as np
from netCDF4 import Dataset

from create_orac_luts import _publish_lut, _scattering_cache_path, private_workdir


def test_private_workdirs_isolate_actual_intermediate_files(tmp_path):
   output = tmp_path / "luts"
   output.mkdir()
   with private_workdir() as first:
      with private_workdir() as second:
         assert first != second
         first_scat = _scattering_cache_path(output, first, reuse_scat=0)
         second_scat = _scattering_cache_path(output, second, reuse_scat=0)
         assert first_scat != second_scat
         np.savez(first_scat, owner=np.asarray([1]))
         np.savez(second_scat, owner=np.asarray([2]))
         (first / "timestamp.txt").write_text("first\n")
         (second / "timestamp.txt").write_text("second\n")
         assert np.load(first_scat)["owner"].tolist() == [1]
         assert np.load(second_scat)["owner"].tolist() == [2]
         assert (first / "timestamp.txt").read_text() == "first\n"
         assert (second / "timestamp.txt").read_text() == "second\n"
      second_path = second
   first_path = first
   assert not first_path.exists()
   assert not second_path.exists()
   assert list(output.iterdir()) == []


def test_reuse_scat_reads_only_the_explicit_persistent_cache(tmp_path):
   output = tmp_path / "luts"
   work = tmp_path / "private"
   assert _scattering_cache_path(output, work, reuse_scat=0) == work / "scatfile.npz"
   assert _scattering_cache_path(output, work, reuse_scat=1) == output / "scatfile.npz"


def test_publication_validates_and_does_not_overwrite_existing_product(tmp_path):
   output = tmp_path / "luts"
   output.mkdir()
   private = tmp_path / "private.nc"
   final = output / "product.nc"
   with Dataset(private, "w") as dataset:
      dataset.createDimension("x", 1)
      dataset.createVariable("value", "f4", ("x",))[:] = [1.0]
   _publish_lut(private, final)
   assert final.is_file() and not private.exists()

   private_again = tmp_path / "private-again.nc"
   with Dataset(private_again, "w") as dataset:
      dataset.createDimension("x", 1)
      dataset.createVariable("value", "f4", ("x",))[:] = [2.0]
   try:
      _publish_lut(private_again, final)
   except FileExistsError:
      pass
   else:
      raise AssertionError("publication unexpectedly overwrote an existing LUT")
