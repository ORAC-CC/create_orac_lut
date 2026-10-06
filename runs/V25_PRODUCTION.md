# V25 ORAC cloud LUT production

## What V25 is

V25 is the V24 cloud LUT product set recomputed with a vertical cloud
temperature profile in the thermal-emission calculation.

**Unchanged from V24:**
- platforms, instruments, channels and microphysical models;
- forward model, atmosphere, gas, Rayleigh and DISORT settings (NSTR 60);
- the Mie / Legendre / radius-grid numerics of V24;
- the Grid B product sampling: `liquid-water-cloud-grid-b.lut` and
  `ice-cloud-grid-b.lut`;
- every LUT dimension and variable.  The reflection and transmission
  operators are the same calculation as V24; only the emissivity `E_md` of
  thermal and mixed channels changes.

**The V25 scientific change** (`create_orac_luts.py`, "Versions and numerical
changes" 5; `src/oraclut/cloud_temperature.py`):

| Item | V24 | V25 |
|---|---|---|
| Cloud temperature in the emission DISORT call | one value for every in-cloud layer (the IDL's 250 K; the value does not affect `E_md`) | the cloud top at the fixed reference state T_top = 240 K, p_top = 628 hPa; the temperature rises downward along the phase-pure saturated adiabat of the substance |
| Liquid-water LUTs | isothermal | saturated liquid-water adiabat (Murphy and Koop 2005 vapour pressure over supercooled water, L_v of Rogers and Yau 1989) |
| Ice LUTs | isothermal | saturated ice adiabat (Murphy and Koop 2005 vapour pressure over ice and latent heat of sublimation) |
| Optical depth to the thermodynamic path length | – | dz = dτ₀.₅₅ / β, i.e. z(x) = x / β for the cumulative 0.55 µm optical depth x below the cloud top; z is a path-length coordinate along the reference adiabat, not an altitude, not limited, not a LUT variable; z = 0 for τ₀.₅₅ = 0 |
| β (representative volume extinction coefficient at 0.55 µm) | – | 20 km⁻¹ (liquid water), 1 km⁻¹ (ice) |
| Pressure along the adiabat | – | evolved hydrostatically with the adiabat's own temperature, dp/dz = p g / (R_d T), from p_top = 628 hPa; the LUT cloud's position in the atmosphere does not enter |
| Emission layering | the one in-cloud atmosphere layer | 12 equal-optical-depth sub-layers of each in-cloud layer, with the adiabat temperature at every boundary |
| `E_md` normalisation | B(T) of the isothermal cloud | B(240 K); `E_md` may exceed 1 |
| LUT dimensions | – | unchanged: no cloud-temperature, cloud-base-temperature or effective-temperature dimension |

The extinction coefficients are deliberate model constants of V25.  They are
not the per-particle extinction cross-sections of the microphysics, nor a
normalised size-distribution extinction: the LUT fixes no particle number
concentration, so the microphysics does not define a physical thickness.  The
physical cloud structure is defined from the 0.55 µm optical depth (the LUT
coordinate) and is the same for every spectral channel.

The adiabat of the LUT's phase is followed throughout, whatever temperature
is reached: no melting, no 273.15 K ceiling, no phase switch. The reference
state (240 K, 628 hPa) is a convention of the LUT; retrieved clouds have other
top temperatures, and the temperature contrast within the cloud is that of the
reference adiabat. Extreme reference states at the deepest optical-depth
nodes (τ₀.₅₅ = 256: z = 12.8 km for liquid water, 256 km for ice) are a
limitation of this simplified V25 treatment and are documented, not limited;
the dependence on the actual cloud-top temperature and mixed-phase physics
are deferred.

The selection is the run-file setting `cloud_temperature_profile = 'adiabatic'`.
A run file without it, or with `'isothermal'`, reproduces the V24 calculation
exactly (the V24 run files are unchanged and still produce V24).

## Products

There are 40 products, as in V24:
- Aqua and Terra MODIS (channels 1–36), and Sentinel-3A and -3B SLSTR (channels
  1–9, dual view);
- each with six liquid-water models (240, 253, 263, 273, old, stg) and four ice
  models (sph, agg, ghm, src).

Run files are `runs/<platform>_<instrument>_cloud_<model>_v25.run` (the V24
run files with `version = 25` and `cloud_temperature_profile = 'adiabatic'`),
and products are written to

    /network/group/aopp/eodg/RGG004_GRAINGER_ORACFILE/ORAC_LUTS/<platform>_<instrument>_m_<substance>_a01_p<shortname>_v25.nc

next to, and never over, the `_v23.nc` and `_v24.nc` products.  The generator
refuses to overwrite an existing product.  V25 files carry global attributes
(`cloud_temperature_profile`, `cloud_top_temperature_K`, `cloud_adiabat`,
`cloud_top_reference_pressure_hPa`, `cloud_extinction_coefficient_055um_per_km`,
`cloud_emission_sublayers`)
recording the treatment, and the `valid_range` of `E_md` is no longer 0–1.

## Submission and provenance

    scripts/submit_v25_cloud_luts.sh --dry-run      # check all 40 run files; submit nothing
    scripts/submit_v25_cloud_luts.sh                # submit every product not yet present

The script submits nothing unless:
- the production paths are unmodified, and
- HEAD is contained in `origin/main` (it fetches first).

It forwards each job to atmlxint7 over ssh (`scripts/submit_oraclut_lut.sh`),
which needs the caller's ssh agent (`SSH_AUTH_SOCK`) or a valid Kerberos
ticket.  Resources are those of V24 (`--time=72:00:00 --mem=24G`; the
thermodynamic profile costs milliseconds and the sub-layered emission calls
are a negligible share of the DISORT time).

Each submission is appended to `validation/v25/production_submissions.tsv`,
which records time, host, Git revision, SLURM job ID, run file and product.
Each job log (`validation/slurm/v25_*.out`) records the revision and the state
of the source tree when the job started.

The jobs execute this shared working tree.  **Do not change production source
until every job has started.**

## Validation

`validation/v25/results/REPORT_v25_cloud_temperature_profile.md` records the
V24 trace, the tests, the compact-grid V24/V25 comparison, the forward-model
(top-of-atmosphere) differences, the Nakajima–King comparisons, the runtime
and the production runs.  `validation/v25/check_products.py` checks the
finished products.
