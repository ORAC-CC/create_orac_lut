# Double-precision V21/V22 numerical baseline

Status (2026-10-11): **implementation complete and validated on compact
cases; complete Terra MODIS products and the V21/V22 equivalence benchmark
pending** (section 6).  The new double-precision V21 IDL tables are a
*candidate* reference: nothing in the shared ORAC LUT directory has been
replaced, and the historical single-precision V21 products remain the
archived reference until the owner approves the replacement.

## 1. Why

The historical single-precision DISORT 2.1 kernel is ill-conditioned for some
direct-beam calls (UPBEAM/SGECO "matrix near singular"): at those calls it
returns radiances wrong by O(1), sometimes negative, with reflected flux up to
132 % of incident.  They occur inside ORAC's operational angular domain and at
the 80 degree solar-zenith nodes that ORAC interpolates for 70 <= SZA < 75.
Double precision removes them; at the original 60 streams its arithmetic
error is ~3e-8 against a REAL*16 build and its stream-count error is below
0.08 sigma of the MODIS measurement uncertainty.  Evidence:
`validation/v21_v22_terra_modis/disort_boundary/REPORT_disort_boundary.md`,
`.../REPORT_operational_angular_domain.md`, `.../REPORT_idl_double_precision.md`.

## 2. Architecture

```
create_orac_lut/disort2/src/          unmodified DISORT 2.1 source (shared by everything)
   DISORTtoIDL.so  (historical)       single-precision IDL DLM, preserved byte for byte
create_orac_lut/disort2_dp/           double-precision build of the same source
   DISORTtoIDL_DP.so                  IDL DLM: DISORT_DP (same 44 arguments as DISORT)
   libdisort_dp.so                    same kernel objects for the Python generator
create_orac_lut/call_disort.pro       one IDL code path; selects the solver
src/oraclut/radiative_transfer/legacy_disort.py   one Python code path; selects the solver
```

`ORAC_LUT_DISORT_PRECISION` selects the solver in both generators:
`double` (default when unset) or `single` (the historical kernel).  Everything
upstream of the solver — Mie scattering, phase moments (NMOM 1000), spectral
and atmospheric processing, layer compression, LUT grids, NetCDF schema and
float32 storage — is unchanged.  The stream count stays 60.  ORAC holds LUT
coefficients as default REAL (4 bytes; ORAC `common/common_constants.F90`
`sreal = 4`), so float64 storage would gain nothing; the solver returns float64
and each generator rounds once into its float32 arrays, as IDL does.

| Product | Generator | Precision | Release definition |
|---|---|---|---|
| V21 historical | IDL, `call_disort_v01.pro` path | single | `releases/idl_lut_v21_single.release` |
| **V21 (new)** | IDL, `call_disort.pro` (`call_disort_v02.pro`) | **double** (default) | `releases/idl_lut_v21_double.release` |
| V22 historical | Python, `lut-code-v22.1` files (`*_v01.py`) | single | `releases/lut_v22.release` (sealed) |
| **V22 (new)** | Python, active files incl. IDL-compatible `TOTAL` | **double** (default) | `releases/lut_v22_double.release`; `CODE_VERSION` `lut-code-v22.2` |

## 3. Source revisions (SourceVault)

Changed production files, each preserved with SourceVault before and after
the change (file versions are independent of LUT versions):

| File | Historical | New |
|---|---|---|
| `create_orac_lut/call_disort.pro` | `call_disort_v01.pro` | `call_disort_v02.pro` (precision selection) |
| `src/oraclut/radiative_transfer/legacy_disort.py` | `legacy_disort_v01.py` | `legacy_disort_v02.py` (double-precision path, selector) |
| `create_orac_luts.py` | `create_orac_luts_v02.py` | `create_orac_luts_v03.py` (provenance records the precision) |
| `CODE_VERSION` | `CODE_VERSION_v01` (`lut-code-v22.1`) | `CODE_VERSION_v02` (`lut-code-v22.2`, `v22-double`) |
| release launcher | `create_orac_luts_release_v01.py` | `create_orac_luts_release_v03.py` (manifest selection, double library, sealed CODE_VERSION served from its SourceVault copy) |
| `scripts/run_oraclut_lut.slurm`, `scripts/submit_oraclut_lut.sh` | `_v01` | launcher v03, precision/manifest in the job log |

New files: `create_orac_lut/disort2_dp/*` (`_v01` copies of the sources),
`scripts/make_idl_release_manifest.py`, `scripts/make_lut_v22_double_manifest.py`,
`scripts/check_release_manifest.py`, `scripts/verify_staged_lut.py`,
`scripts/compare_lut_numeric.py`, `scripts/run_idl_v21_lut.sh`,
`scripts/run_idl_v21_lut.slurm`, `scripts/submit_v21_v22_double.sh`,
`scripts/monitor_v21_v22_double.sh`, `runs/*_v22_double.run`,
`tests/test_disort_precision.py`.  The historical compiled DLM is also kept as
`create_orac_lut/disort2/src/DISORTtoIDL_v01.so`.

Exact binaries: `create_orac_lut/disort2_dp/BUILD_INFO.txt` (compiler, flags,
source, object and library SHA-256); every file, input, runtime library and
historical product is listed with SHA-256 in the release manifests.  Verify:

```shell
python scripts/check_release_manifest.py releases/idl_lut_v21_single.release releases/idl_lut_v21_double.release
ORAC_LUT_RELEASE_MANIFEST=lut_v22_double.release python create_orac_luts_release_v03.py --check 22
ORAC_LUT_RELEASE_MANIFEST=lut_v22.release        python create_orac_luts_release_v03.py --check 22
```

## 4. Reproducing the historical single-precision products

**V21 (IDL).**  `ORAC_LUT_DISORT_PRECISION=single` makes `call_disort.pro`
execute the original `DISORT` statement with the historical DLM; everything
before it is the original code.  The four compact validation cases reproduce
the historical compact V21 LUTs bitwise in this mode (integration job, section
6), and those compact LUTs equal the production V21 products bitwise in all
RT coefficients.  For the complete products:

```shell
cd create_orac_lut && ORAC_LUT_DISORT_PRECISION=single bash terra_modis_run
```

(run where IDL is permitted; `terra_modis_run` and `initpath.bash` are
preserved as `terra_modis_run_v01` and `initpath_v01.bash`).  The exact
historical source set is `releases/idl_lut_v21_single.release`.

**V22 (Python).**  The sealed historical release still runs from the working
tree:

```shell
ORAC_LUT_RELEASE_MANIFEST=lut_v22.release python create_orac_luts_release_v03.py runs/<product>_v22.run
```

(the launcher serves the sealed `CODE_VERSION` from `CODE_VERSION_v01`); or
check out tag `lut-code-v22.1` in a temporary worktree.

## 5. Producing the double-precision products

```shell
bash scripts/submit_v21_v22_double.sh          # from atmlxint1/2/3/5/6 or atmlxint7
bash scripts/monitor_v21_v22_double.sh [--verify]
```

The script refuses unless every release file verifies against its manifest
and is committed.  It submits to SLURM partition `priority-eodg` from
atmlxint7: the IDL integration validation (four compact cases, double and
single precision), then the three V21 products (only if the validation
passed), and independently the three V22 products.  Products go only to the
protected staging directories `staging/lut_v21_double/` and
`staging/lut_v22_double/` (git-ignored), each with a provenance record; V21
products are verified against the historical V21 products (schema,
dimensions, coordinates, finiteness, NC_STRING attributes) and made
read-only.  Repeating the command skips staged products.

## 6. Validation and the equivalence benchmark

Done (compact Terra MODIS p240 cases A-D: SZA 80 failure, channel 20, a V21
negative reflectance, an ordinary node):

* IDL `DISORT_DP` removes the single-precision failures (17-70 sigma), agrees
  with a REAL*16 build to 3e-8 and converges with stream count to 0.08 sigma at
  NSTR 60 (`REPORT_idl_double_precision.md`).
* The production double-precision kernel objects are bitwise identical to the
  validated build.
* Python `_disort_double` reproduces IDL `DISORT_DP` bitwise for all 90
  captured IDL calls.
* The complete Python generator with double precision agrees with the IDL
  double-precision compact LUTs to <= 2.8e-6 in every RT coefficient
  (relative <= 6e-5 at the smallest reflectances, otherwise <= 2e-6); the
  remaining differences come from 1-3 ulp input differences upstream of
  DISORT (`validation/v21_v22_terra_modis/v22_double/`).
* The double-precision V22 release produces, through the launcher, a product
  bitwise identical to the direct generator run, with the precision recorded.

Pending (after the SLURM jobs): the IDL integration job (double within 1
float32 ulp of the validated results, single bitwise historical); the six
complete products; and the benchmark
`validation/v21_v22_terra_modis/v21_v22_double/compare_double_pairs.py`, which
compares each V21/V22 pair at every node and, through ORAC's own
interpolation (linear in angle, bicubic in log10 tau and r_e) at random
operational states, the forward model (solar reflectance, thermal and mixed
brightness temperature) and both Jacobians against the production MODIS
uncertainty.  Acceptance, fixed in advance: every difference <= 0.1 sigma
(Jacobians scaled by steps of 0.1 in log10 tau and 1 um in r_e).

## 7. Archival and clean-up policy

Kept permanently: SourceVault copies of every production source file; the
historical V21/V22 release manifests and products (shared directory, listed
with SHA-256); the historical DLM; the validated compact test cases and their
reference results under `validation/v21_v22_terra_modis/idl_double/`; the
final reports; the comparison and diagnostic utilities; scientific inputs and
reference papers.  The single-precision comparison reports are marked
historical, not deleted.  Generated products, scratch builds, job logs and
caches stay out of Git (`.gitignore`).  Only unambiguously disposable
generated material is deleted; anything else is inventoried for the owner.

Promotion of the new V21 tables into the shared ORAC LUT directory, and any
replacement of existing products there, requires the owner's explicit
approval.
