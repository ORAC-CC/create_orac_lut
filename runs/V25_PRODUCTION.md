# V25 ORAC cloud LUT production

## What V25 is

V25 is the V24 cloud LUT product set recomputed with a vertically varying
cloud temperature in the thermal-emission calculation:

- **ice LUTs:** the supplied vertically inhomogeneous cirrostratus profile of
  P. Watts (OCA / EUMETSAT), `cloud_vertical_profile = 'cirrostratus'`;
- **liquid-water LUTs:** the saturated liquid-water (wet) adiabat from the
  reference cloud top over the path length z = τ₀.₅₅ / 20 km⁻¹,
  `cloud_vertical_profile = 'wet_adiabat'`.

Both start from the reference cloud-top temperature 240 K and normalise
`E_md` by B(240 K). The Grid B optical-depth endpoint τ₀.₅₅ = 256 is retained
(accepted by the user).

**Unchanged from V24:**
- platforms, instruments, channels and microphysical models;
- forward model, atmosphere, gas, Rayleigh and DISORT settings (NSTR 60);
- the Mie / Legendre / radius-grid numerics of V24;
- the Grid B product sampling: `liquid-water-cloud-grid-b.lut` and
  `ice-cloud-grid-b.lut`;
- every LUT dimension and variable. The reflection and transmission
  operators are the same calculation as V24; only the emissivity `E_md` of
  thermal and mixed channels changes.

**The supplied profile** (`references/data/ocalut_cloudprofile_Cirrostratus.dat`,
md5 `7d7b472f3ff9681901c03afe85c3202b`):

- origin `/tcenas/home/pwatts/ctpo2/METimCTP_DELIVERY_2/data/cloudprofiles.hdf`,
  cloud type Cirrostratus, supplied by Alessio in connection with Phil
  Watts' vertically inhomogeneous LUT tests for OCA / EUMETSAT;
- 13 total cloud optical depths COT = 0.0625, 0.125, …, 256 (doubling), 99
  rows each: `% Extinction`, `layer dCOT`, `Cumul COT`, `Km-below CT`,
  `Temp above CT`;
- established numerically (`validation/v25`): the extinction shape is the
  same for all 13 COTs and sums to 100 %; `dCOT_i = COT × f_i / 100`; the rows
  are equally spaced in depth from the cloud top to the base; `ΔT = 8 K/km × z`
  in every row; the cloud depth H(COT) rises from 2.0 km (COT 0.0625) to
  11.0 km (COT 128 and 256) — the cloud does not deepen without limit.

**The V25 ice treatment** (`create_orac_luts.py`, "Versions and numerical
changes" 5; `src/oraclut/cloud_temperature.py`):

| Item | V24 | V25 (ice LUTs) |
|---|---|---|
| Cloud temperature in the emission DISORT call | one value for every in-cloud layer (the IDL's 250 K; the value does not affect `E_md`) | the cloud top at the reference temperature T_top = 240 K; T = 240 K + ΔT(F; COT) at the boundaries of equal-optical-depth emission layers, with ΔT from the supplied profile |
| Vertical structure | homogeneous | dτ_i = COT f_i, z_i = H(COT) i/98, ΔT_i = 8 K/km z_i (the supplied model); the relation between the cumulative optical-depth fraction F and ΔT follows from the rows (midpoint rule) |
| Interpolation to the LUT optical-depth grid | – | H(COT) linear in COT between the 13 supplied values (the data are exactly linear in COT from 0.0625 to 1); held at 2 km below COT 0.0625; not extrapolated above 256 (the Grid B maximum); the extinction shape needs no interpolation; every supplied profile is reproduced exactly at its own COT |
| Emission layering | the one in-cloud atmosphere layer | 40 equal-optical-depth sub-layers (DISORT accepts 46), boundaries at F = k/40: total optical depth and the cumulative distribution exact at every boundary; chosen from the layer-count convergence |
| `E_md` normalisation | B(T) of the isothermal cloud | B(240 K); `E_md` may exceed 1 |
| LUT dimensions | – | unchanged: no cloud-temperature dimension |

**The V25 liquid-water treatment** (`WetAdiabatProfile`):

| Item | V25 (liquid-water LUTs) |
|---|---|
| Reference cloud-top state | T_top = 240 K, p_top = 628 hPa (a reference convention, not an altitude) |
| Thermodynamic path length | dz = dτ₀.₅₅ / β with β_ext,055 = 20 km⁻¹, i.e. z = τ₀.₅₅ / 20 km⁻¹; no limit, no H / H_max, no altitude or surface constraint |
| Temperature | the saturated liquid-water (wet) adiabat throughout (Murphy and Koop 2005 vapour pressure over supercooled water, Rogers and Yau 1989 latent heat, pseudo-adiabatic lapse rate), followed continuously with no freezing transition, phase switch or ceiling; the pressure is evolved hydrostatically along the path in the adiabat's own temperature; 4th-order Runge–Kutta, 10 m steps |
| Production grid extreme | τ₀.₅₅ = 256 → z = 12.8 km, 317.5 K, 2957 hPa (finite, monotonic, converged) |
| Emission layering | 40 equal-optical-depth sub-layers, as for ice |

Both treatments prescribe temperature structure, not microphysics: every
sub-layer has the LUT particle model's single-scattering albedo and phase
function, so the diffuse and direct-beam calls (one homogeneous cloud layer)
are unchanged and the reflection and transmission operators are bitwise
V24. The reference temperature 240 K is a convention of the LUT: both
treatments give a temperature **departure** from the cloud top, retrieved
clouds have other top temperatures, and the dependence of the normalised
`E_md` on the reference temperature (quantified in the report, § 8) is an
accepted limitation of V25.

**History.** Two earlier V25 experiments are superseded and recorded in the
report: (1) a constant extinction coefficient β with a depth cap H_max and a
saturated adiabat for both phases (revision `7f6b4e2`; its 40 production jobs
488174–488213 were cancelled by the user); (2) the uncapped z = τ/β adiabat
for both phases (revision `13b9df3`), whose ice branch gave pathological
reference states for thick ice cloud (1630 K at COT 256) and whose liquid
branch is the definitive liquid treatment above. The ice formulation was then
replaced by the supplied cirrostratus profile (`5722280`).

The selection is the run-file setting `cloud_vertical_profile`:
`'cirrostratus'` in the 16 ice run files, `'wet_adiabat'` in the 24
liquid-water run files. A run file without it, or with `'isothermal'`,
reproduces the V24 calculation exactly; the generator refuses a backend for
the other phase (cirrostratus + liquid water, wet_adiabat + water ice).

## Products

There are 40 configurations, as in V24:
- Aqua and Terra MODIS (channels 1–36), and Sentinel-3A and -3B SLSTR (channels
  1–9, dual view);
- each with six liquid-water models (240, 253, 263, 273, old, stg) and four ice
  models (sph, agg, ghm, src).

Run files are `runs/<platform>_<instrument>_cloud_<model>_v25.run` (the V24
run files with `version = 25` and `cloud_vertical_profile`; 16 ice, 24
liquid-water products), and products are written to

    /network/group/aopp/eodg/RGG004_GRAINGER_ORACFILE/ORAC_LUTS/<platform>_<instrument>_m_<substance>_a01_p<shortname>_v25.nc

next to, and never over, the `_v23.nc` and `_v24.nc` products. The generator
refuses to overwrite an existing product. V25 files carry global attributes
recording the treatment (`cloud_vertical_profile`, `cloud_top_temperature_K`,
`cloud_emission_layers`; ice: `cloud_vertical_profile_file`, `_md5`,
`_origin`, `_type`; liquid water: `cloud_top_reference_pressure_hPa`,
`cloud_extinction_coefficient_055um_per_km`), and the `valid_range` of
`E_md` is no longer 0–1.

## Submission and provenance

    scripts/submit_v25_cloud_luts.sh --dry-run      # check all 40 run files; submit nothing
    scripts/submit_v25_cloud_luts.sh [modis|slstr|all] [ice|liquid|both]   # submit every product not yet present

The script submits nothing unless the production paths are unmodified and
HEAD is contained in `origin/main` (it fetches first). It forwards each job
to atmlxint7 over ssh (`scripts/submit_oraclut_lut.sh`), which needs the
caller's ssh agent (`SSH_AUTH_SOCK`) or a valid Kerberos ticket. Resources
are those of V24 (`--time=72:00:00 --mem=24G`).

Each submission is appended to `validation/v25/production_submissions.tsv`
(time, host, Git revision, SLURM job ID, run file, product); the rows of the
cancelled jobs 488174–488213 remain there for provenance. Each job log
(`validation/slurm/v25_*.out`) records the revision and the state of the
source tree when the job started. The jobs execute this shared working tree:
**do not change production source until every job has started.**

## Validation

`validation/v25/results/REPORT_v25_cloud_temperature_profile.md` records the
supplied data and its interpretation, the tests, the compact-grid V24/V25
comparison, the forward-model (top-of-atmosphere) differences, the
Nakajima–King comparisons, the layer convergence, the runtime and the
production status. `validation/v25/check_products.py` checks finished
products; `validation/v25/monitor_jobs.sh` shows job states.
