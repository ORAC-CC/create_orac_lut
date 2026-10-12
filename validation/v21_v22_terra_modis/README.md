# V21 / V22 Terra MODIS comparison tools

The V21/V22 translation is complete (`../../V21_V22_REFERENCE.md`).  This
directory keeps only the tools still in use:

* `compare_luts.py`: variable-by-variable comparison of two ORAC LUTs;
* `radiance_impact.py`, `operational_domain.py`: ORAC forward-model and
  operational-domain analysis of a LUT pair;
* `v21_v22_single/compare_full_luts.py`: the definitive V21/V22 product check
  (results in `v21_v22_single/results/full_luts/`);
* `v21_v22_single/launcher_smoke.run`: compact smoke test of the V22 release
  launcher (`tests/test_single_precision_release.py`).

The full investigation (reports, diagnostics, IDL and Python captures,
intermediate results) is preserved in Git at commit `f9b0b3b` and tag
`lut-code-v22.4`.
