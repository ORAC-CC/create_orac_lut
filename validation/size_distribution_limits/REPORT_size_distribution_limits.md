# ORAC particle-size integration-limit investigation

## Executive conclusions

The premise is only partly correct. The fixed **0.001–100 µm radius** domain is
used for the modified-gamma liquid-cloud Mie path. It is passed for all Mie
calls, but the legacy routine deliberately ignores those two values for a
log-normal aerosol and derives mode-dependent bounds instead.

For the repository's liquid-water distribution (effective variance
0.1111111), a common-grid tail experiment finds that an upper limit of about
**3 times effective radius** is sufficient for the cases tested under the
high-fidelity diagnostic band defined below. This is not equivalent to lowering
100 µm in every case:

| effective radius (µm) | candidate upper limit, 3 re (µm) | radius-point change at 0.55 µm versus 100 µm |
|---:|---:|---:|
| 5 | 15 | about 85% fewer |
| 10 | 30 | about 70% fewer |
| 15 | 45 | about 55% fewer |
| 20 | 60 | about 40% fewer |
| 30 | 90 | about 10% fewer |
| 40 | 120 | about 20% more |

The current 100 µm limit is therefore unnecessarily wide for the smaller
cloud radii, appropriate to moderate accuracy for 30–40 µm, but slightly too
narrow for strict convergence of the 40 µm case. Against a 200 µm common-grid
reference, the 40 µm case at 100 µm has worst errors of 4.11e-4 in extinction,
3.99e-4 in phase-function L1 distance, and 2.11e-4 in moments 0–127 across the
six wavelengths. At 120 µm those worst errors fall to 1.82e-5, 2.17e-5, and
1.04e-5 respectively.

The 0.001 µm lower cloud limit is physically much lower than necessary: the
modified-gamma number density is proportional to r^6 near zero, and removing
everything below 0.05 µm made no measurable physical difference on the dense
common grid. It should nevertheless be retained in production for now.
Changing a bound moves every point of the current linear trapezoidal grid;
with `xres=0.4`, this created differences as large as 0.9% in extinction even
when the removed lower tail was physically zero. Raising the lower bound saves
essentially no points compared with reducing the upper bound.

The existing log-normal aerosol limits are scientifically and numerically
preferable to a fixed 0.001–100 µm domain. The three representative production
modes used here have legacy bounds:

| mode | mode radius / spread | legacy bounds (µm) |
|---|---|---|
| small (`aerosol_a79.mm`) | 0.070 / 1.700 | 0.01358–1.4431 |
| accumulation (`volcanic-ash_ctn.mm`) | 0.217 / 1.770 | 0.03717–5.0677 |
| coarse (`aerosol_a75.mm`) | 0.788 / 1.822 | 0.12342–20.1252 |

Forcing the small mode onto 0.001–100 µm is not a safer reference: the broad
linear grid under-resolves the distribution, producing errors up to 88% in
extinction relative to its operational adaptive calculation. The aerosol
formula should not be replaced by fixed bounds.

**Recommendation:** do not change production from this study alone. Prototype
an adaptive modified-gamma upper bound of `3 * effective_radius` in a separate
branch, retain the 0.001 µm lower bound, and couple that change to a radius-grid
convergence improvement or an invariant-node integration design. Validate the
result through full spectral-response integration and compact LUT/retrieval
comparisons before approval. Retain the current log-normal aerosol bounds.

## 1. Current implementation traced

### Active Python/IDL-mirror LUT path

`create_orac_luts.py` calls `generate_scattering_properties` from both its cloud
and aerosol drivers (calls at lines 463–464 and 896–897). That routine obtains
the LUT-dependent component radius and calls `create_bwgp` at 0.55 µm and at
each channel wavelength (`src/oraclut/idl_mirror/generate_scattering_properties.py`,
lines 310–345).

`src/oraclut/idl_mirror/create_bwgp.py`, lines 228–234, calls
`mie_size_dist_new` with `[rm, s, 0.001, 100.0]` and `xres=0.4`. Its legacy IDL
source does the same at `create_orac_lut/create_bwgp.pro`, line 86.

The important distribution branch is in both:

- `src/oraclut/idl_mirror/create_bwgp.py`, lines 80–103;
- `create_orac_lut/mie_size_dist_new.pro`, lines 143–177.

For `modified_gamma`, elements 2 and 3 of the parameter array become the lower
and upper radius limits. For `log_normal`, they are ignored and the code uses

```text
tq = GAUSS_CVF(0.999) = -3.090232...
lower = mode_radius * spread**tq
upper = 4 * mode_radius * spread**(-tq)
```

The separate validated pipeline implements the same behavior in
`src/oraclut/optics/legacy.py`, lines 36–93. Its caller is
`src/oraclut/pipeline.py`, lines 900–921. Thus both current Python paths agree
with the authoritative IDL semantics.

### Radius grid and quadrature

The number of radius points is

```text
max(int(2*pi*(upper-lower)/(wavelength*xres)), 200),  xres = 0.4
```

(`src/oraclut/idl_mirror/create_bwgp.py`, lines 97–103; IDL lines 167–177).
The integration is a uniformly spaced, linear-radius trapezoidal rule. The IDL
trapezoidal weights are defined in `create_orac_lut/quadrature.pro`, lines
75–100, then mapped from [-1,1] to the radius interval by
`create_orac_lut/shift_quadrature.pro`, lines 7–12.

The comment in the authoritative IDL says a size-parameter step of 0.1 is
required for accurate calculation, while the production wrapper explicitly
chooses 0.4 for speed. Endpoint changes therefore alter both truncation and the
sample locations of an already coarsened oscillatory Mie integral.

### Distribution and bulk averaging

`liquid-water_stg.mm`, lines 11–20, selects a single spherical Mie component
with a Hansen–Travis modified-gamma distribution and effective variance
0.1111111. For each LUT effective radius, `generate_scattering_properties`
passes that radius as the distribution parameter.

The IDL distribution equations are at `mie_size_dist_new.pro`, lines 179–191.
Particle area weighting, extinction, scattering, single-scattering albedo,
asymmetry, and scattering-weighted phase matrix are calculated at lines
239–254. The Python mirror is at `src/oraclut/idl_mirror/create_bwgp.py`, lines
105–145. The validation code uses those same equations and calls the preserved,
unchanged Fortran Mie kernel through `mie_single_batch`.

## 2. Study design

The main survey contains:

- liquid clouds at effective radii 5, 10, 15, 20, 30, and 40 µm;
- actual repository log-normal small, accumulation, and coarse aerosol modes;
- 0.55 µm, MODIS 0.646 and 1.628 µm, SLSTR 0.868 µm, MODIS 3.780 µm, and
  MODIS 11.026 µm;
- requested upper limits 5, 10, 20, 40, 60, and 100 µm;
- requested lower limits 0.0001, 0.001, 0.005, 0.01, and 0.05 µm;
- extinction, scattering, single-scattering albedo, asymmetry, phase function,
  and ORAC-convention Legendre moments 0–127;
- a targeted 1000-moment calculation for the 40 µm cloud and coarse aerosol at
  0.55 µm.

Two complementary calculations are reported:

1. `convergence.csv` recomputes the production trapezoidal grid for every
   domain. It measures the result that a literal production bound change would
   produce, including shifted-node numerical effects.
2. `fixed_grid_tail_convergence.csv` masks a shared, denser `xres=0.1` radius
   grid. It isolates the physical contribution of the omitted tail without
   moving all other Mie samples.

The fixed cloud reference was extended to 200 µm in
`cloud_extended_upper_convergence.csv`. Results at 150 and 200 µm differ by at
most 1.18e-7 in extinction, 1.65e-7 in phase L1, and 7.08e-8 in moments, so
200 µm is adequate for this case matrix.

### Diagnostic bands

These are reporting bands, not approved ORAC error budgets. They make the
effect of a threshold choice explicit:

| band | relative extinction/scattering | absolute SSA/g | phase L1 | maximum absolute moment error | interpretation |
|---|---:|---:|---:|---:|---|
| high fidelity | 1e-4 | 2e-5 | 1e-3 | 1e-4 | numerical-change target; optical-depth scaling below 0.01% |
| routine validation | 1e-3 | 2e-4 | 5e-3 | 5e-4 | optical-property differences likely visible only in stringent comparisons |
| screening | 1e-2 | 2e-3 | 2e-2 | 1e-2 | useful for rejecting bad limits, not for approving production |

The 3-re cloud bound meets the high-fidelity band across the survey: worst
errors are 4.24e-5 extinction, 7.21e-5 scattering, 8.58e-6 SSA, 1.03e-5 g,
8.43e-5 phase L1, and 3.98e-5 in moments. A 2.5-re bound does not meet the
routine band for every quantity (worst scattering 1.39e-3 and moment error
6.74e-4). Thus the recommendation changes if a looser, quantity-specific error
budget is adopted; 3 re is the first tested factor that clears all listed
high-fidelity diagnostics.

## 3. Convergence findings

### Cloud upper tail

On the fixed dense grid, the worst errors over all radii and wavelengths are:

| upper (µm), versus 200 µm | relative extinction | phase L1 | max moment 0–127 |
|---:|---:|---:|---:|
| 60 | 7.85e-2 | 5.71e-2 | 2.56e-2 |
| 80 | 7.00e-3 | 6.34e-3 | 2.84e-3 |
| 100 | 4.11e-4 | 3.99e-4 | 2.11e-4 |
| 120 | 1.82e-5 | 2.17e-5 | 1.04e-5 |
| 150 | 1.18e-7 | 1.65e-7 | 7.08e-8 |

These worst values are controlled by the 40 µm distribution. For 30 µm,
100 µm gives only 2.02e-6 extinction error; for radii at or below 20 µm the
100 µm tail is negligible.

The cumulative tables show why a universal 60 µm bound fails: for 40 µm
clouds, radii containing 99%, 99.9%, and 99.99% of the extinction contribution
within 0.0001–100 µm are approximately 77.0, 91.6, and 98.5 µm. This is a clear
case in which a small number tail controls an area-weighted optical quantity.

### Lower tail

For modified-gamma clouds, the physical contribution below 0.05 µm is below
the resolution of the dense-grid experiment. The nonzero operational errors
when changing 0.001 to another lower bound are quadrature-node effects, not a
material small-particle tail. Keeping 0.001 avoids those changes at no
meaningful computational cost.

For aerosols, adding particles below the existing adaptive lower bounds changes
extinction by at most 1.5e-6 in this survey. Raising the lower limit to 0.05 µm,
however, changes extinction by as much as 1.18e-2 because it cuts into the
small/accumulation distributions. The current mode-dependent lower formula is
therefore justified.

### Aerosol upper tail

The existing small and accumulation upper bounds are already about 1.44 and
5.07 µm. For the coarse mode, truncating its 20.125 µm reference at 10 µm gives
1.27% extinction error, 1.33e-2 phase L1 error, and 6.24e-3 moment error.
At 20 µm those fall to 3.61e-6, 2.23e-6, and 1.00e-6. The present adaptive
formula covers all three representative modes without paying for 100 µm.

### Phase functions and high Legendre moments

The phase and moment tests prevent choosing a limit from extinction alone. For
the cloud 3-re study, phase and moments converge at the same factor as the bulk
properties, while 2.5 re misses the chosen bands.

The 1000-moment 0.55 µm stress test also shows that the tail can affect high
orders more strongly than the first 128 moments. For the 40 µm cloud at 60 µm,
maximum moment errors are 3.59e-3 (orders 0–127), 1.06e-2 (128–511), and
1.06e-2 (512–999). At 100 µm they are 1.18e-4, 8.96e-5, and 1.46e-5; at
120 µm all three bands are below 7.43e-6. For the coarse aerosol, 20 µm versus
20.125 µm gives maximum errors of 2.08e-7, 1.67e-7, and 2.45e-16 in those
bands.

## 4. Computational impact

At wavelengths where the 200-point floor does not apply, point count is nearly
proportional to `(upper-lower)/wavelength`. Measured median Mie speed-ups for
cloud domains of 5, 10, 20, 40, and 60 µm relative to 100 µm were approximately
79, 49, 20, 5.7, and 2.7 respectively over the survey. These are kernel timings
on one shared system and should be treated as indicative; Mie cost also grows
with size parameter, so speed-up can exceed the radius-point ratio.

At long wavelengths or narrow aerosol domains, the 200-point floor limits or
eliminates the saving. Conversely, expanding a small aerosol from its adaptive
domain to 100 µm was slower and less accurate because it spread too few linear
nodes over the part of the domain containing particles.

An adaptive 3-re cloud bound would concentrate the saving on the numerous
small-radius calculations, while making the 40 µm endpoint more expensive. A
complete production saving depends on the LUT effective-radius grid and channel
mix and was not extrapolated into a claimed wall-clock percentage.

## 5. Failure regimes and limitations

- The cloud recommendation applies to the repository's effective variance
  0.1111111. A broader modified-gamma distribution needs a larger multiple of
  effective radius; the formula must depend on variance before general use.
- Log-normal limits depend on both mode radius and spread. A fixed multiple of
  effective radius was not substituted for the established quantile formula.
- Wavelength changes Mie oscillations, absorption, and the production point
  count. The distribution controls the broad domain, but convergence must be
  checked over the instrument spectrum and the optical quantity of interest.
- The aerosol cases are spherical Mie representatives. The Dubovik T-matrix
  branch below 6 µm uses its tabulated integration machinery, not this radius
  domain; at thermal wavelengths the wrapper falls back to Mie.
- Mixtures were not recombined into complete aerosol classes. This isolates
  integration behavior but is not a retrieval-level validation.
- The main survey uses 128 phase angles/moments, supplemented by two targeted
  1000-moment visible cases. It does not prove all 1000 moments for every
  wavelength and radius.
- No spectral-response convolution, DISORT calculation, LUT generation, or
  retrieval inversion was run. Translation of optical errors into radiance or
  retrieved-state error remains future work.
- The current `xres=0.4` linear trapezoid produces endpoint-sensitive numerical
  differences. An adaptive limit should not be deployed without addressing or
  explicitly validating this interaction.

## 6. Reproduction and outputs

From the repository root:

```bash
MPLCONFIGDIR=validation/tmp/matplotlib-cache \
PYTHONPATH=src /home/g/grainger/miniforge3/envs/science/bin/python \
  validation/size_distribution_limits/run_study.py --phase-order 128

PYTHONPATH=src /home/g/grainger/miniforge3/envs/science/bin/python \
  validation/size_distribution_limits/run_high_order_check.py
```

Numerical results and plots are under
`validation/size_distribution_limits/results/`. Important tables are:

- `convergence.csv`: literal production-grid endpoint changes and timings;
- `fixed_grid_tail_convergence.csv`: physical tail differences on common nodes;
- `cloud_extended_upper_convergence.csv`: fixed upper limits versus 200 µm;
- `cloud_scaled_upper_convergence.csv`: 2–4 re limits versus 200 µm;
- `high_order_1000_moment_check.csv`: targeted production-order stress checks;
- `cumulative_quantiles.csv`: number/extinction/scattering contribution radii;
- `kernel_timing.csv`: measured Mie timings and point counts.

The cumulative PNG files plot the distribution, cumulative number/extinction/
scattering, and selected cumulative ORAC moments for every studied
distribution. `upper_limit_convergence.png` and
`lower_limit_convergence.png` summarize operational-grid error envelopes and
timings.
