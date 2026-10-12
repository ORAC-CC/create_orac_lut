"""DISORT precision selection after the V21/V22 finalisation.

Single precision (the original V21 Intel kernel) is the default and the
reference; double precision (the suspended experimental development) is
available only by explicit selection.  The IDL call_disort.pro is the
original V21 file.
"""

import os
import subprocess
import sys
from pathlib import Path

import numpy as np
import pytest

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))
from oraclut.radiative_transfer import legacy_disort as L   # noqa: E402

PYTHON = sys.executable


def _layered_case():
   """A small three-layer Rayleigh/cloud/Rayleigh column like the compressed ORAC calls."""
   nmom = 99
   pmom = np.zeros((nmom + 1, 3), np.float32)
   pmom[0, :] = 1.0
   pmom[2, 0] = pmom[2, 2] = 0.1
   g = 0.85
   pmom[:, 1] = g ** np.arange(nmom + 1)
   dtauc = np.array([0.01, 8.0, 0.02], np.float32)
   ssalb = np.array([1.0, 0.999, 1.0], np.float32)
   utau = np.array([0.0, dtauc.sum()], np.float32)
   umu = np.array([-0.9, -0.5, 0.5, 0.9], np.float32)
   phi = np.array([0.0, 90.0, 180.0], np.float32)
   return dtauc, ssalb, pmom, utau, umu, phi


def test_precision_selection(monkeypatch):
   monkeypatch.delenv("ORAC_LUT_DISORT_PRECISION", raising=False)
   monkeypatch.delenv("ORAC_LUT_DISORT_SP_KERNEL", raising=False)
   assert L.disort_precision() == "single" and L.single_kernel() == "intel_dlm"
   monkeypatch.setenv("ORAC_LUT_DISORT_PRECISION", "double")
   assert L.disort_precision() == "double"
   monkeypatch.setenv("ORAC_LUT_DISORT_PRECISION", "quad")
   with pytest.raises(ValueError):
      L.disort_precision()


@pytest.mark.skipif(not L.DOUBLE_LIBRARY.is_file(), reason="double-precision DISORT library not built")
def test_double_and_single_paths(monkeypatch):
   # The historical single-precision kernel overwrites SSALB = 1 in the caller's
   # array (DISORT's dither); production passes freshly built arrays, and so does
   # this test.  The double-precision path must leave its inputs untouched.
   args = _layered_case()
   monkeypatch.setenv("ORAC_LUT_DISORT_PRECISION", "double")
   d = L.disort(*args, 100.0, 0.6, 0.0, nstreams=16)
   assert np.array_equal(args[1], _layered_case()[1]), "double-precision path modified SSALB"
   monkeypatch.setenv("ORAC_LUT_DISORT_PRECISION", "single")
   s = L.disort(*_layered_case(), 100.0, 0.6, 0.0, nstreams=16)
   assert d["uu"].dtype == np.float64 and s["uu"].dtype == np.float32
   # Downward intensity at the top is zero for a pure beam source (exact boundary condition).
   assert np.max(np.abs(d["uu"][:2, 0, :])) == 0.0
   assert abs(d["rfldir"][0] - 60.0) < 1e-5
   # This well-conditioned case: the two precisions agree to single-precision accuracy.
   assert np.max(np.abs(d["uu"] - s["uu"])) < 1e-3
   # Repeat calls are bitwise identical.
   monkeypatch.setenv("ORAC_LUT_DISORT_PRECISION", "double")
   assert np.array_equal(L.disort(*_layered_case(), 100.0, 0.6, 0.0, nstreams=16)["uu"], d["uu"])


def test_call_disort_pro_is_the_original_v21_file():
   # The V21 IDL reference: call_disort.pro is the original (single-precision)
   # file; the experimental precision-selecting version is call_disort_v02.pro.
   assert (ROOT / "create_orac_lut/call_disort.pro").read_bytes() == (ROOT / "create_orac_lut/call_disort_v01.pro").read_bytes()
   assert "DISORT_DP" in (ROOT / "create_orac_lut/call_disort_v02.pro").read_text()


def test_launcher_requires_manifest_choice_for_v22():
   both = subprocess.run([PYTHON, "create_orac_luts_release_v03.py", "--check", "22"], cwd=ROOT,
                         capture_output=True, text=True, env={**os.environ, "ORAC_LUT_RELEASE_MANIFEST": ""})
   assert both.returncode != 0 and "ORAC_LUT_RELEASE_MANIFEST" in both.stderr
   chosen = subprocess.run([PYTHON, "create_orac_luts_release_v03.py", "--check", "22"], cwd=ROOT, capture_output=True,
                           text=True, env={**os.environ, "ORAC_LUT_RELEASE_MANIFEST": "lut_v22_double.release"})
   assert chosen.returncode == 0, chosen.stdout + chosen.stderr


def test_double_release_refuses_single_precision(monkeypatch):
   sys.path.insert(0, str(ROOT))
   import importlib
   launcher = importlib.import_module("create_orac_luts_release_v03")
   monkeypatch.setenv("ORAC_LUT_RELEASE_MANIFEST", "lut_v22_double.release")
   manifest = launcher.find_release_manifest(22)
   monkeypatch.setenv("ORAC_LUT_DISORT_PRECISION", "single")
   monkeypatch.setitem(sys.modules, "oraclut.optics.legacy_mie", sys.modules.get("oraclut.optics.legacy_mie") or importlib.import_module("oraclut.optics.legacy_mie"))
   monkeypatch.setitem(sys.modules, "oraclut.radiative_transfer.legacy_disort", L)
   with pytest.raises(launcher.ReleaseRefusal):
      launcher.verify_libraries(manifest)


@pytest.mark.parametrize("manifest", ["releases/idl_lut_v21_single.release", "releases/idl_lut_v21_double.release",
                                      "releases/lut_v22_single-r2.release", "releases/lut_v22_single-r2.dependencies"])
def test_release_manifests_verify(manifest):
   shared = Path("/network/group/aopp/eodg/RGG004_GRAINGER_ORACFILE/ORAC_LUTS")
   if not shared.is_dir() or not (ROOT / "create_orac_lut/disort2_dp/DISORTtoIDL_DP.so").is_file():
      pytest.skip("shared LUT directory or DP build absent")
   result = subprocess.run([PYTHON, "scripts/check_release_manifest.py", manifest], cwd=ROOT, capture_output=True, text=True)
   assert result.returncode == 0, result.stdout
