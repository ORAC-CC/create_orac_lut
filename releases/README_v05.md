# Release index

| LUT release | Status | Manifest | Generator | Tag |
|---|---|---|---|---|
| **V21** | **definitive IDL reference** | `idl_lut_v21_single.release` | IDL, `create_orac_lut/terra_modis_run` | `idl-lut-v21` |
| **V22** | **definitive Python reference** (bitwise = V21) | `lut_v22_single-r2.release` | `create_orac_luts_release_v04.py` | `lut-code-v22.4` |
| **V23** | **grid-only successor of V22**: `lut-code-v22.4` code, Grid B (`../V23_RELEASE.md`) | `lut_v23.release` | `create_orac_luts_release_v04.py` | (tag after validation) |
| **V24** | **V23 with refined Mie scattering properties**: `lut-code-v24.1` (3.5 r_e liquid limit, nested radius refinement, adaptive Legendre expansion, NMOM = max(1000, L)), Grid B (`../V24_RELEASE.md`) | `lut_v24.release` | `create_orac_luts_release_v04.py` | (tag after validation) |
| **V25** | **V24 with vertically varying cloud temperature in the thermal emission**: `lut-code-v25.2` (cirrostratus profile for ice, wet adiabat for liquid water, 250 K cloud top; geometric emission sub-layers, isothermal cells for thin clouds), Grid B (`../V25_RELEASE.md`) | `lut_v25.release` | `create_orac_luts_release_v04.py` | (naming review) |
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

`lut_v23.release`, `lut_v24.release` and `lut_v25.release` pin their input files as well (entries
`input-*`), including the Grid B definitions; `lut_v23.dependencies` lists the
V23 products' inputs and external shared objects (written by
`scripts/make_release_dependencies.py` from the products' provenance), and
`lut_v24.dependencies` those of the V24 products (2026-10-15).  Manifests are sealed: never edit them.  Several manifests define V22, so the
launcher needs `ORAC_LUT_RELEASE_MANIFEST=<file>`.  See `../V21_V22_REFERENCE.md`.
