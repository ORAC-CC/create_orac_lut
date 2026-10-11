# IDL `TOTAL` compatibility in the ORAC Python LUT generator

## Result

The Python `idl_total` function reproduces the `TOTAL` behavior required by
the reachable V21 LUT calculation: sequential accumulation, IDL accumulation
and result precision, one-based dimensions, first-index-fastest flattening,
`/CUMULATIVE`, and `/DOUBLE`.  All translated `TOTAL` operations in both the
active `create_orac_luts.py`/`idl_mirror` path and the older package pipeline
now use it.  Reductions translating other IDL operations were not changed.

The sealed historical `lut-code-v22.1` `_v01` files and V22 release manifest
remain immutable and therefore still contain the original reductions.  This
change is the active working-source `_v02` successor; a future release must
explicitly select these numbered revisions before a release launcher will use
them.  No existing V22 product or sealed source was rewritten.

All 36 Terra MODIS centre wavelengths and all 36 centre wavenumbers calculated
from the original SRFs are bitwise identical to the values stored in the V21
IDL product.  Before this change, 34 wavelengths and 28 wavenumbers differed.
This validation does not regenerate a LUT, propagate the change through
radiative transfer, or establish complete V21/V22 equivalence.

## IDL semantics established

The primary semantic source is the official [NV5 IDL `TOTAL`
documentation](https://www.nv5geospatialsoftware.com/docs/total.html).  It
specifies that dimensions are one-based; absent/zero dimension reduces the
whole array; default accumulation/result are FLOAT except for DOUBLE,
DCOMPLEX, COMPLEX and 64-bit integer inputs; `/DOUBLE` forces double
accumulation; and `/CUMULATIVE` preserves the input shape.

Those rules were checked directly with the installed IDL 8.9.0 using
[idl_total_reference.pro](../diagnostics/idl_total_reference.pro).  The exact
[captured results](idl_89_reference.txt) cover float32, float64, int16,
one-dimensional and two-dimensional arrays, both relevant dimensions,
cumulative output, `/DOUBLE`, result dtypes and a cancellation-sensitive
sequence.  In particular:

- `TOTAL([1e8, 1, 1, -1e8, 3])` is FLOAT `3`;
- the same call with `/DOUBLE` is DOUBLE `5`;
- `TOTAL([1d16, 1d0, 1d0, -1d16, 3d0])` is DOUBLE `3`.

The installed IDL had `!CPU.TPOOL_MIN_ELTS=100000`.  Every production `TOTAL`
operand in this call graph is smaller: the Mie radius integral is bounded below
25,000 points by its size-parameter limit, `legpexp` permits at most 20,000,
and atmospheric, SRF, angular and component vectors are much smaller.  Thus
the V21 calls use the deterministic non-threaded accumulation order.  The
reference diagnostic explicitly used `/TPOOL_NOTHREAD`.

## Python interface

The reusable implementation is [src/oraclut/idl_total.py](../../src/oraclut/idl_total.py):

```python
idl_total(array, dimension=None, *, cumulative=False, double=False)
```

- `dimension=None` or `0` flattens in Fortran/IDL order and sums sequentially.
- A positive dimension is one-based and removes that axis unless
  `cumulative=True`.
- Sequential accumulation uses `numpy.cumsum`; it deliberately does not use
  `numpy.sum`, whose pairwise reduction tree changes rounding.
- Default real accumulation is float32 except for float64 and 64-bit integer
  inputs; `/DOUBLE` uses float64.  Complex inputs retain/promote their complex
  precision according to the documented rule.
- Empty and non-numeric inputs and invalid dimensions fail explicitly.

The import at `oraclut.idl_mirror.total` re-exports the common implementation,
and `oraclut.idl_mirror.idl_total` is the public mirror-package convenience
name.

Unsupported because the V21 LUT calculation never requests them:
`/INTEGER`, `/NAN`, `/PRESERVE_TYPE`, a numeric `DOUBLE=0` keyword value, and
emulation of IDL's thread-pool-dependent alternate reduction orders.  No
production call uses a dimension or `/CUMULATIVE`; those forms are implemented
and IDL-validated to make ordering and reuse explicit.

## V21 `TOTAL` inventory

All reachable calls have no Dimension argument and do not use `/CUMULATIVE`.
The only keyword in production is `/DOUBLE` in `create_range.pro:62-63`.
Shapes below are the one-dimensional slices actually passed to `TOTAL`, even
where their parent variables are multidimensional.

| V21 source | Argument(s), type and shape | Keywords | Python counterpart |
|---|---|---|---|
| `integrate_trapeziodal.pro:5` | `dx*ym`, FLOAT `[n_srf-1]` | none | `idl_mirror/load_srfstrarr.py:36` |
| `input_files/simpson_integral.pro:4` | `x*y`, `[n_srf]`; FLOAT for the denominator and DOUBLE when IDL solar interpolation promotes the numerator | none | `idl_mirror/load_srfstrarr.py:46`; existing Python input/interpolation precision is otherwise unchanged |
| `create_orac_cloud_lut.pro:351,452,568,660` | FLOAT `scatreltau[nlayer]`, `dtau[nlayer]`, `emTau[ncloudlayer]`, `srf.val[nwvl_max]` | none | `create_orac_luts.py:519,632,716,759` and `pipeline.py:191,318,361` |
| `create_orac_aerosol_lut.pro:343,480,591,658` | same types/shapes as cloud | none | `create_orac_luts.py:960,1062,1130,1168` and `pipeline.py:426,456,496` |
| `create_range.pro:62-63` | component expressions `[ncomp]`; inputs originate as FLOAT | `/DOUBLE` | `idl_mirror/generate_scattering_properties.py:130-131`, `double=True` |
| `create_range.pro:107,110` | FLOAT component ratios for one equal-size group `[ngroup]` | none | `generate_scattering_properties.py:119` |
| `create_range.pro:150-151,163-164,181-182,199-200,209-210` | FLOAT effective-radius diagnostic numerators/denominators `[nmode]` | none | no Python calculation: `type_re` is a verbose-only diagnostic never used in any output |
| `mie_size_dist_new.pro:242-244` | DOUBLE radius-integral terms `[nsize]` | none | `idl_mirror/create_bwgp.py:135-137` and `optics/legacy.py:88-90` |
| `mie_size_dist_new.pro:246` | DOUBLE backscatter terms `[nsize]` | none | no counterpart: optional `Bbac` is never requested by the LUT caller |
| `mie_size_dist_new.pro:250-253` | DOUBLE phase-matrix terms `[nsize]`, once for each Stokes component and angle | none | scalar LUT uses only F11: `create_bwgp.py:142`, dimension 1; `optics/legacy.py:91-93` |
| `mie_size_dist_new.pro:257-260` | DOUBLE distribution integrals `[nsize]` | none | requested `Vavg`: `create_bwgp.py:145`, `optics/legacy.py:94`; optional `Gavg`, `Ravg`, `RVW` are not requested |
| `legpexp.pro:30-31,43` | DOUBLE angular terms `[nmom]` | none | `create_bwgp.py:162,165,173` and `optics/legacy.py:112,115,120` |
| `generate_scattering_properties.pro:217,233,248,268` | phase function times DOUBLE quadrature weights `[nmom]` | none | vectorized over remaining axes at `generate_scattering_properties.py:292,304`, dimension 1 |
| `generate_scattering_properties.pro:313,315,317,319,321` | DOUBLE component mixtures `[ncomp]` at 0.55 um | none | `generate_scattering_properties.py:387-393`; phase vectorization uses dimension 2 |
| `generate_scattering_properties.pro:330,332,334,336,338` | DOUBLE component mixtures `[ncomp]` for each channel/SRF point | none | `generate_scattering_properties.py:399-405`; phase vectorization uses dimension 2 |
| `dubovik/dubovik_lognormal.pro:135-136,142-143,156,176,179` | DOUBLE radius or angle vectors `[nradius]`, `[ntheta]` or `[ndqv-1]` | none | `idl_mirror/dubovik.py:346-365` |
| `dubovik/dubovik_lognormal_multiple_eps.pro:21-23,65-66,69,74` | DOUBLE aspect-ratio mixtures `[neps]` | none | `dubovik.py:378,386,398-405`; phase vectorization uses dimension 1 |
| `read_baran.pro:125` | FLOAT size-bin volume/area terms `[n_Dm_baran]` | none | no counterpart: the Baran reader remains an unavailable V22 particle-property path |
| `write_v2_lut.pro:10,14,18` | integer/byte channel flags `[nchannel]`, result FLOAT | none | `idl_mirror/write_v2_lut.py:85,89,93` |

The old `generate_scattering_properties_simple.pro`, `disort2/ExampleDISORT.pro`,
`input_files/makebbconstantscode.pro`, attic sources and comments containing
the English word "total" are not in the operational LUT call graph.  They were
not treated as V22 call sites.

## Reductions deliberately not replaced

The two remaining `np.sum` calls in `idl_mirror/segment.py` implement the IDL
`INT_TABULATED` algorithm and a validated segmented-SRF integration helper.
They do not translate `TOTAL`, so replacing them would conflate different IDL
semantics.  Ordinary validation/statistics reductions elsewhere are likewise
not scientific translations of V21 `TOTAL` and remain NumPy operations.

## Validation

[tests/test_idl_total.py](../../tests/test_idl_total.py) checks the captured IDL
8.9 results, output dtypes, dimensions, cumulative operation, flattening order,
invalid inputs, all 36 Terra MODIS centres and a source audit preventing the
translated call sites from returning to `np.sum`.

The integration calculation uses `terra_modis_v1.inst`, all original Terra
MODIS SRFs and `Gueymard2018.sssi`, with `srf_quad=1`.  Its bitwise V21 target
is the stored coordinate vectors from
`terra_modis_m_liquid-water_a01_pold_v21.nc` (SHA-256
`d65a05113d7ca9675154639648e0e9f011f75f21b53d6ce71abd087f94feeb45`).
Results:

| Quantity | Elements | Bitwise differences after change | Maximum absolute difference |
|---|---:|---:|---:|
| central wavelength | 36 | 0 | 0 |
| central wavenumber | 36 | 0 | 0 |
| `F0` (additional available intermediate) | 26 | 20 | 2.44140625e-4 |

`F0` improves from 22 to 20 differing values but is not expected to become
identical from summation order alone: the prior investigation established that
the Python path also differs in solar-table input precision and interpolation
arithmetic.  Those issues were not changed here.

Relevant integration tests pass after making the older validation helpers use
the same `idl_total` implementation.  One unrelated existing test remains
unrunnable because `docs/idl_python_port_status.md` is absent, and one legacy
runner test is independently unrunnable because its reference-case driver is
absent.  Neither missing file was invented or repaired as part of this task.
The applicable source-versioning subset passes 45 tests; release-error tests
whose outcome depends on unrelated manifest-message wording or the existing
dirty/uncommitted Git state were deselected.  Its transition assertions were
updated so sealed V22 remains checked against `_v01` and `lut-code-v22.1`,
while the active working source is checked against its latest SourceVault
revision rather than being required to remain identical to V22 forever.

## Scope conclusion

Summation-order differences have been eliminated as a source of disagreement
for every represented V21 `TOTAL` operation.  No relevant summation
compatibility issue remains within the supported interface.  The unavailable
Baran path and optional unused Mie outputs have no V22 computation to update.
This result makes no claim about DISORT or complete LUT equivalence.
