# Release index

| LUT release | Status | Manifest | Generator | Tag |
|---|---|---|---|---|
| **V21** | **definitive IDL reference** | `idl_lut_v21_single.release` | IDL, `create_orac_lut/terra_modis_run` | `idl-lut-v21` |
| **V22** | **definitive Python reference** (bitwise = V21) | `lut_v22_single-r2.release` | `create_orac_luts_release_v04.py` | `lut-code-v22.4` |
| V22 | superseded: historical gfortran V22 | `lut_v22.release` | `create_orac_luts_release_v03.py` | `lut-code-v22.1` |
| V22 | superseded: corrected V22, revision 1 | `lut_v22_single.release` | `create_orac_luts_release_v04.py` | commit `e4036d0` |
| V22 | experimental: double-precision DISORT | `lut_v22_double.release` | `create_orac_luts_release_v03.py` | commit `37cf766` |
| V21 | experimental: double-precision IDL DISORT | `idl_lut_v21_double.release` | IDL (`call_disort_v02.pro`) | — |

Supporting records of the definitive V22 release:

* `lut_v22_single-r2.dependencies`: every input file and external shared object
  of the three products, with SHA-256 (verify with
  `python scripts/check_release_manifest.py releases/lut_v22_single-r2.dependencies`);
* `lut_v22_single-r2.environment.txt`: the conda environment (`conda list --explicit`);
* `products/lut_v22/`: the products' provenance records.

Manifests are sealed: never edit them.  Several manifests define V22, so the
launcher needs `ORAC_LUT_RELEASE_MANIFEST=<file>`.  See `../V21_V22_REFERENCE.md`.
