# Terra MODIS V21/V22 LUT differences inside ORAC's operational angular domain

> **Historical (single-precision V21/V22).**  This report analyses the V21 IDL
> and V22 Python products made with the historical single-precision DISORT.
> Since 2026-10-10 both generators use double-precision DISORT by default;
> see `DOUBLE_PRECISION_BASELINE.md` (repository root) for the current
> baseline.  The findings below remain the record of the single-precision
> products, which stay archived.

## Summary

Classification: **2. A small number of isolated significant discrepancies**,
but in `p240` some of them are operationally severe.

- Restricting to the domain ORAC actually uses (day pixels, solar zenith
  < 75 degrees) does **not** remove the worst `p240` failure.  ORAC
  interpolates linearly between the 70 and 80 degree LUT nodes for
  70 <= SZA < 75, giving the failed 80 degree call up to (but not reaching)
  half weight.  At SZA -> 75 the interpolated `R_0v` error is 0.271 and the
  black-surface measurement error 147 sigma (848 sigma at the raw node).
- Independently of the 80 degree node, both LUT pairs contain isolated failed
  direct-beam calls **inside** the domain: in `p240`, mixed channel 20 at
  SZA 10 and 30 degrees (up to 15.4 K, 22 sigma, in brightness temperature)
  and solar channels 6, 7 and 12 at SZA 30-70 degrees (up to 40 sigma); in
  `pold`, channels 13 and 16 at SZA 40 and 60 degrees (5 sigma); in `pagg`,
  channel 11 at SZA 70 degrees (3 sigma).
- Outside those calls the differences are widespread but small: the median
  in-domain call maximum of `|dR_0v|` is 1.2-1.5e-5 and 99.9 % of calls are
  below 1e-3.  Thermal-only channels never exceed 0.002 sigma.  Fewer than
  0.03 % of solar and mixed measurement cases exceed 1 sigma.
- V21 is not uniformly the correct side: in two of the largest in-domain
  `p240` failures V21 has negative direct-beam reflectance and V22 does not.
- The structure (isolated calls, 10^2-10^4 times their neighbours) matches the
  single-precision DISORT instability confirmed for the 80 degree call
  ([DISORT boundary report](disort_boundary/REPORT_disort_boundary.md)).
  Double-precision testing of the candidate calls listed in section 6 is
  justified; it was not performed here.

Everything here is a reanalysis of the existing LUTs; no LUT was regenerated
and no production code, IDL source or DISORT library was changed.

## 1. The operational limits in ORAC

Source: ORAC-CC/orac `master` at `eb9a1233c1136ef95fa0c37f650fe519e9134ed4`,
cloned read-only to the session scratch area on 2026-10-10; the tree was
clean and `git ls-remote` confirmed upstream `master` at the same commit.
Nothing was copied into this repository.

| Item | ORAC behaviour | Source |
|---|---|---|
| Solar zenith limit | `Ctrl%MaxSolZen = 75.` default | `src/read_driver.F90:614` |
| Day / twilight | day iff `SolZenMin <= Sol < MaxSolZen` (**exclusive** at 75); twilight `MaxSolZen <= Sol <= Sunset` | `src/determine_illum.F90:57-64` |
| Use of solar LUT terms | only day pixels run the solar forward model; twilight uses pure thermal channels only, night all thermal channels (mixed channels without their solar part) | `src/fm.F90` header, lines 9-37 |
| Satellite zenith limit | `Ctrl%MaxSatZen = 90.` default; pixels rejected if `SatZen > MaxSatZen` (inclusive) | `src/read_driver.F90:615`, `src/get_geometry.F90:118-131` |
| Configurability | both limits are fields of `Ctrl_t`, and the driver parser is generated from the whole type, so a driver may override them | `src/Makefile:213-223`, `src/read_driver.F90:1234` |
| Angular bracketing | lower node `locate(..., bounded=.true.)`, upper node = lower + 1; weights linear in degrees | `src/set_gzero.F90:95-171`, `src/locate.F90` |
| Angular interpolation | trilinear in satellite zenith, solar zenith, relative azimuth, applied to the LUT coefficients before the forward model; bicubic in optical depth and effective radius | `src/int_lut_tausatsolazire.F90:117-170`, `src/int_lut_routines.F90:97-231` (`LUTIntSelm` default bicubic, `read_driver.F90:624`) |
| Transmission along the view path | `T_00` and `T_0d` are also interpolated at the satellite zenith on the solar-zenith grid (`iSaZSoZ0`) | `src/set_crp_solar.F90:146-158`, `src/int_lut_tausatreonsol.F90` |

So the 75 degree solar limit is a code default and is exclusive.  The 75 degree
satellite limit is an operational choice, not the code default (90); for Terra
MODIS it has no effect, because the LUT satellite-zenith grid ends at the
instrument's 65 degrees and the satellite-zenith screen is applied to the
measured angle.

## 2. LUT nodes actually used

Terra MODIS grids (identical in all six files): solar zenith 0-90 in 10 degree
steps; satellite zenith 0-65 in 10 equal steps; relative azimuth 0-180 in
18 degree steps.

| Solar zenith node | Class | Largest ORAC weight | Used by |
|---|---|---|---|
| 0-70 | in domain | 1 | all day pixels |
| 80 | interpolation support | < 0.5 (supremum 0.5 at SZA -> 75) | day pixels with 70 <= SZA < 75 only |
| 90 | cannot influence a retrieval | 0 | — |

All satellite-zenith and relative-azimuth nodes are in the domain.  For the
view-path transmission (`T_00`, `T_0d` on the solar grid at the satellite
zenith) MODIS angles up to about 65 degrees use solar-grid nodes up to 70
degrees.  Optical depth and effective radius are retained in full, including
the bicubic neighbours.

Because ORAC interpolates coefficients linearly, the error at an operational
geometry with 70 <= SZA < 75 is `(1 - w) dF(70) + w dF(80)`, `w < 0.5`, and
its extreme values lie at the interval ends.  The analysis therefore uses the
nodes 0-70 plus interpolated slices at SZA 72.5 (w = 0.25) and 75- (w = 0.5,
an upper bound).

## 3. Coefficient differences in the operational domain

Maximum absolute difference, RMS, and number of elements over 1e-3 / 1e-2
(operational set; for SZA-dependent variables this includes the 72.5 and 75-
slices).  Maximum relative difference uses elements with |V21| >= 0.01.

| Case | Variable | Max abs | RMS | > 1e-3 | > 1e-2 | Max rel | Location of max abs |
|---|---|---:|---:|---:|---:|---:|---|
| pold | `R_0v` | 0.0372 | 4.8e-5 | 432 | 26 | 0.119 | ch13, τ 8, r_e 5, SZA 40, VZA 43.3, RAA 180 |
| pold | `R_0d` | 0.0034 | 7.1e-5 | 17 | 0 | 0.027 | ch1, τ 256, r_e 19, SZA 60 |
| pold | `T_0d` | 0.0013 | 4.4e-5 | 1 | 0 | 0.035 | ch13, τ 0.016, r_e 39, SZA 0 |
| p240 | `R_0v` | 0.271 | 6.1e-4 | 4416 | 718 | 3.96 | ch6, τ 256, r_e 29, SZA 75-, VZA 0, RAA 180 |
| p240 | `R_0d` | 0.601 | 2.5e-3 | 26 | 6 | 1.23 | ch6, τ 256, r_e 29, SZA 75- |
| p240 | `T_0d` | 0.0083 | 5.9e-5 | 14 | 0 | 0.096 | ch9, τ 0.125, r_e 21, SZA 30 |
| pagg | `R_0v` | 0.0052 | 1.7e-5 | 381 | 0 | 0.068 | ch9, τ 16, r_e 93, SZA 50 |
| pagg | `R_0d` | 0.0129 | 8.0e-5 | 12 | 1 | 0.084 | ch11, τ 0.016, r_e 33, SZA 70 |
| pagg | `T_0d` | 0.0529 | 2.2e-4 | 6 | 3 | 0.170 | ch11, τ 0.016, r_e 33, SZA 70 |

`T_00` (direct transmission) differs by at most 4e-6.  The diffuse and
thermal operators `T_dv`, `R_dv`, `T_dd`, `R_dd` differ by at most 1.1e-3 and
`E_md` by at most 7e-6 in all cases.  For comparison, the raw 80 degree node
reaches `R_0v` 0.541 and `R_0d` 1.20 in `p240`; in `pold` and `pagg` the 80
degree node is unremarkable (max `R_0v` 2.7e-3 and 2.0e-4).

Coordinate profiles (maximum and RMS difference per coordinate value along
every dimension) are in `results/operational_domain/<case>_coordinate_dependence.csv`.
They show no systematic trend with optical depth, effective radius, channel,
solar or satellite zenith.  The maxima are set by individual calls.  The `R_0v`
RMS per coordinate value is 1e-5 to 2e-4 in `pold` and 1e-5 to 3e-5 in
`pagg`; in `p240` it rises to 1e-3 to 3e-3 only at the coordinate values that
contain the failed calls (r_e 29 and 11-19, τ 1-256, SZA 75-).

## 4. Direct-beam calls

A call is one (solar channel, τ, r_e, SZA node) DISORT direct-beam solve,
covering all 10 x 11 viewing directions.  A call is **isolated significant**
when its largest `|dR_0v|` exceeds 1 sigma of the production solar
uncertainty (`R >= 0.01`) and is at least 10 times the median of its six
neighbours in τ, r_e and SZA.  Mixed channels carry no SNR in the LUT, so for
them only the production calibration terms (2.2 % of R) are used at this
level; their measurement significance is assessed in brightness temperature
in section 5.

| Case | In-domain calls | Isolated significant (in domain / 80 deg support) | Not isolated, significant | Median / p99.9 / max in-domain call max |
|---|---:|---|---:|---|
| pold | 70,720 | 15 / 9 | 1 / 0 | 1.5e-5 / 4.2e-4 / 0.037 |
| p240 | 70,720 | 46 / 24 | 0 / 3 | 1.4e-5 / 9.9e-4 / 0.086 |
| pagg | 84,864 | 11 / 0 | 0 / 0 | 1.2e-5 / 2.4e-4 / 0.0052 |

Significant calls are 0.01-0.07 % of all calls and are isolated by factors of
10^2-10^4; there is no population of moderately elevated calls between them
and the floating-point background.  See
[the call-difference plot](../../plots/terra_modis_v21_v22_operational_call_differences.png).

## 5. Operational measurement impact

Representative, not a full retrieval.  The forward model and uncertainties are
those of [`radiance_impact.py`](radiance_impact.py) (validated mirror of the
McGarragh et al. (2018) formulation in `validation/v24/disort_options/orac_fm.py`,
production MODIS uncertainty model): solar black surface at every LUT view
geometry; Lambertian surfaces of albedo 0.05, 0.2 and 0.8 at nadir view;
thermal and mixed channels with the repository mid-latitude-summer atmosphere.
Coefficients are interpolated to the operational solar-zenith set first, as
ORAC does.  With the full grid and identity weights the new code reproduces
the stored `radiance_impact.py` results for all 870 rows exactly.

| Case | Measurement | Cases | > 1 sigma | Max sigma | Max abs | Worst case |
|---|---|---:|---:|---:|---:|---|
| pold | solar, black surface | 6.18e6 | 68 (1.1e-5) | 5.3 | 0.037 | ch16, τ 2, r_e 11, SZA 60, VZA 36, RAA 126 |
| pold | solar, Lambertian (worst albedo) | 6.8e4 | 0 | 0.27 | 1.5e-3 | ch6, SZA 70 |
| pold | mixed BT | 2.24e6 | 0 | 0.33 | 0.28 K | ch20, τ 0.016, r_e 29, SZA 30 |
| pold | thermal BT | 3.4e4 | 0 | 0.002 | 9e-4 K | ch29 |
| p240 | solar, black surface | 6.18e6 | 1192 (1.9e-4) | 147 | 0.271 | ch6, τ 256, r_e 29, SZA 75- |
| p240 | solar, Lambertian (worst albedo) | 6.8e4 | 12 (1.8e-4) | 146 | 0.271 | ch6, τ 256, r_e 29, SZA 75- |
| p240 | mixed BT | 2.24e6 | 505 (2.3e-4) | 22.3 | 15.4 K | ch20, τ 1, r_e 11, SZA 30, VZA 36, RAA 126 |
| p240 | thermal BT | 3.4e4 | 0 | 0.002 | 9e-4 K | ch28 |
| pagg | solar, black surface | 7.70e6 | 215 (2.8e-5) | 3.1 | 5.2e-3 | ch11, τ 0.016, r_e 33, SZA 70 |
| pagg | solar, Lambertian (albedo 0.8) | 8.2e4 | 3 (3.7e-5) | 2.0 | 0.014 | ch11, τ 0.016, r_e 33, SZA 70 |
| pagg | mixed BT | 2.69e6 | 31 (1.2e-5) | 2.1 | 2.8 K | ch20, τ 0.016, r_e 65, SZA 50 |
| pagg | thermal BT | 4.1e4 | 0 | 0.002 | 1e-3 K | ch29 |

Split by geometry, exact in-domain nodes (SZA <= 70) versus interpolated
70 < SZA < 75:

| Case | Solar max sigma, SZA <= 70 / 72.5-75- | Mixed BT max sigma, SZA <= 70 / 72.5-75- |
|---|---|---|
| pold | 5.3 / 1.3 | 0.33 / 0.14 |
| p240 | 40.2 / 147 | 22.3 / 3.5 |
| pagg | 3.1 / 2.4 | 2.1 / 0.17 |

The per-channel comparison with the full LUT grid is in
[the measurement-impact plot](../../plots/terra_modis_v21_v22_operational_measurement_impact.png):
restricting to the operational domain lowers the `p240` channel-6 maximum
from 848 to 147 sigma and leaves every other channel essentially unchanged,
because the other significant calls are inside the domain.

A large coefficient difference does not always become a large measurement
difference: the `p240` `R_0d` error of 0.60 enters only the two-layer
approach (`Rfbd` in `set_crp_solar.F90:129`), and the thermal channels are
insensitive to all differences found.

## 6. Is double precision needed?

Within the operational domain the single-precision DISORT products show
isolated, call-shaped failures with the signature already confirmed for the
80 degree `p240` call: one solve differs from both its neighbours and the
other LUT version by orders of magnitude, sometimes with negative
reflectance.  They are rare but not negligible: in `p240` they reach 40 sigma
(solar) and 22 sigma (mixed BT) at ordinary geometries, and 147 sigma near
SZA 75 through the support node.  Neither V21 nor V22 is reliable at these
calls; which, if either, is correct is unknown without an independent
calculation.

Candidate calls for independent double-precision testing (all MODIS Terra,
solar zenith node given):

| Case | Channel | τ | r_e | SZA | Max dR_0v | Note |
|---|---:|---:|---:|---:|---:|---|
| p240 | 6 | 256, 128, 64 | 29 | 80 | 0.541, 0.131, 0.043 | support node; 75- interpolation affected |
| p240 | 20 | 16, 8 | 19 | 10 | 0.086, 0.053 | mixed |
| p240 | 20 | 1, 2, 0.5, 4 | 11 | 30 | 0.059, 0.043, 0.004, 0.009 | V21 negative at τ 1 |
| p240 | 7 | 4 | 17 | 50 | 0.027 | |
| p240 | 12 | 2 | 7 | 70 | 0.026 | |
| p240 | 6 | 0.125 | 21 | 30 | 0.012 | V21 negative |
| p240 | 26 | 0.25 | 15 | 30 | 0.003 | 12 sigma |
| pold | 13 | 8 | 5 | 40 | 0.037 | |
| pold | 16 | 2 | 11 | 60 | 0.014 | |
| pagg | 11 | 0.016 | 33 | 70 | 0.003 (`T_0d` 0.053) | |

The complete lists are `significant_calls` in
`results/operational_domain/operational_domain_summary.json` and the ranked
`results/operational_domain/<case>_direct_beam_calls.csv` (top 1000).

## 7. Recommended next steps

Follow-up (2026-10-10): step 1 has been carried out for representative nodes
(the 80 degree `p240` channel-6 call, the channel-20 and channel-6 SZA-30
failures, and an ordinary node) through IDL with a double-precision DISORT
DLM; see [the double-precision report](REPORT_idl_double_precision.md).

1. Recompute the candidate calls above with the existing double-precision
   validation kernel (and single precision at a different stream count), as
   done for the 80 degree call, to establish which version, if either, is
   correct at each.
2. Because ORAC uses the 80 degree node for 70 <= SZA < 75, keep checking it
   for any LUT intended for operational use; the 75 degree limit does not make
   it irrelevant.
3. Extend the comparison to the retrieval Jacobians (bicubic `dF/dlog τ` and
   `dF/dr_e`): an isolated failed call perturbs the gradients over its 4 x 4
   bicubic neighbourhood.  This was not assessed here.
4. If double precision removes these failures, as it did for the 80 degree
   call, the case for a double-precision DISORT production build rests on
   operational-domain evidence, not only on the 80 degree node.

## 8. Limitations

- Representative measurement propagation, not a full retrieval; the
  Lambertian cases are at nadir view only; one atmosphere.
- The call-isolation criterion (>= 1 sigma and >= 10 x neighbour median) is a
  diagnostic choice; partial-call failures affecting few views can be missed
  by it, but are still included in the measurement statistics.
- Terra MODIS only; other instruments (with 90 degree satellite-zenith LUT
  grids) would have satellite-zenith support nodes too.

## 9. Reproduction

```shell
python validation/v21_v22_terra_modis/operational_domain.py --validate
python -m pytest tests/test_operational_domain.py
```

`--validate` checks the full-grid reproduction of `radiance_impact.py`.
Outputs: `results/operational_domain/` and the two plots in `plots/`.
