# Legacy reference-case candidates

## Selection criteria

The first reference should be an actually generated product with all inputs
still available, a compact enough grid to inspect, and enough scientific scope
to expose ordering, SRF, particle optics, and RT conventions. The candidate
must also be identified by the exact wrapper/generator family and command-line
overrides; the filename version is insufficient by itself.

## Recommended sequence

### Candidate A — first principal case: Meteosat-10 SEVIRI cloud, liquid-water `stg`

Reference file:

```text
create_orac_lut/luts/meteosat-10_seviri_cloud/
  meteosat-10_seviri_m_liquid-water_a01_pstg_v21.nc
```

Associated evidence:

- copied output driver: `luts/meteosat-10_seviri_cloud/meteosat-10_seviri_cloud.driver`;
- product provenance marker: `luts/meteosat-10_seviri_cloud/git_revision.txt`
  records `c182ef6`;
- the current `create_orac_lut/meteosat-10_seviri_run` is aerosol-only, so the
  exact historical cloud shell command is not retained as a standalone current
  run file;
- driver source: `input_files/driver/meteosat-10_seviri_cloud.driver`;
- instrument: `input_files/inst/meteosat-10_seviri_v1.inst`;
- model: `input_files/microphysics/liquid-water_stg.mm`;
- LUT grid: `input_files/lut/liquid-water-cloud.lut`;
- atmosphere: selector `2`, midlatitude summer profile;
- SRFs and solar spectrum under `input_files/srf/` and `input_files/sun/`.

Observed NetCDF structure:

```text
optical_depth       17
effective_radius    20
satellite_zenith    10
solar_zenith        10
relative_azimuth    11
channels            11
solar_channels       4
mixed_channels       1
thermal_channels     8
surface_pressure     absent
```

Why it is the best first case:

- it is a V2 product with a complete channel set;
- it exercises solar, mixed, and thermal output branches;
- it uses a single, simple Mie liquid-water component, making the first
  particle-property comparison interpretable;
- it exercises SRF/channel metadata across visible and infrared channels;
- it has no pressure dimension, avoiding the current aerosol split’s additional
  structural branch during the first downstream validation.

Its limitation is equally important: it does not validate the current
pressure-aware aerosol path or gas-enabled aerosol behavior.

The retained product is strong output evidence, but the exact historical shell
invocation is not fully recoverable: the current Meteosat-10 run file is
aerosol-only, and the product provenance marker records the older Git revision
`c182ef6`. It should therefore be treated as the first numerical reference,
with the wrapper/command provenance frozen from the copied driver and source
inspection rather than claimed as a definitively retained `makerunfile_v2.pro`
run.

### Candidate B — second principal case: Meteosat-10 SEVIRI aerosol, `a79`

Reference file:

```text
create_orac_lut/luts/meteosat-10_seviri_aerosol/
  meteosat-10_seviri_m_aerosol_a12_pa79_v21.nc
```

Associated evidence:

- run script: `create_orac_lut/meteosat-10_seviri_run`;
- driver: `input_files/driver/meteosat-10_seviri_aerosol.driver`;
- instrument: `input_files/inst/meteosat-10_seviri_v1.inst`;
- model: `input_files/microphysics/aerosol_a79.mm`;
- LUT grid: `input_files/lut/aerosol.lut`;
- atmosphere: selector `2` in the driver and `atmospheres=2` in the run command;
- gas: `gas=1` in the run command;
- output pressure: one `surface_pressure` value in the observed product;
- selected channels: `[1, 2, 3]` in the driver.

Observed NetCDF structure:

```text
optical_depth       20
effective_radius    20
satellite_zenith    10
solar_zenith        10
relative_azimuth    11
surface_pressure     1
channels             3
solar_channels       3
```

The file is 5,573,699 bytes on disk and has a filesystem timestamp of
2025-11-15 16:04:23 UTC. Its output directory also contains `scatfile.sav`,
`timestamp.txt`, a copied driver, and a `git_revision.txt` marker.

`aerosol_a79.mm` is a two-component Mie model: `waf` and `saf`, with
lognormal size parameters and mixing ratios `0.625` and `0.375`. This is a
good second case for checking multi-component mixing, gas/Rayleigh composition,
pressure indexing, and the aerosol split output schema. It is not a full
infrared aerosol case because the driver selects only channels 1–3, although
the current aerosol top-level code contains thermal-channel handling.

### Candidate C — compact diagnostic: Himawari-8 AHI aerosol, `a70`, V14

Reference file:

```text
create_orac_lut/luts/himawari-8_ahi_aerosol/
  himawari-8_ahi_m_aerosol_a12_pa70_v14.nc
```

Observed structure is deliberately small: 2 optical-depth points, 2
effective-radius points, 2 solar and satellite zenith points, 2 relative
azimuth points, 3 pressure points, and 6 solar channels. It is useful for
debugging NetCDF indexing and pressure dimensions, but it is V14 rather than
V2 and should not be the primary modernisation baseline.

The file is 48,997 bytes on disk and has a filesystem timestamp of
2023-08-26 21:17:47 UTC. Its directory includes the copied driver,
`scatfile.sav`, `timestamp.txt`, and a `git_revision.txt` marker.

The `himawari-8_ahi_run` script is also a useful source-level example because it
calls `create_orac_lut_wrapper_v2` and contains the historical `no_rayleigh` /
`reuse_scat` second-pass commands. No complete V2 Himawari product was found
that would make it preferable to Candidate A or B.

### Candidate D — legacy/compatibility diagnostic: Envisat AATSR

Reference file:

```text
create_orac_lut/luts/envisat_aatsr_cloud/
  envisat_aatsr_m_liquid-water_a01_pold_v50.nc
```

This is a small older product and exercises AATSR’s dual-view configuration,
but it uses an older wrapper/version path (`v50`) and is not the recommended
first V2 case. Keep it for later compatibility tests after the V2 schema is
matched.

The file is 59,892 bytes on disk and has a filesystem timestamp of
2023-09-24 08:19:08 UTC. It is useful as a compact compatibility artifact, not
as evidence of the later V2 generation path.

## Reference artifact inventory

The two full V2 principal products have the following read-only inventory:

| product | format | size | filesystem timestamp | provenance marker |
|---|---|---:|---|---|
| Meteosat-10 cloud `stg` | NetCDF4 | 6,640,898 bytes | 2025-11-13 09:47:15 UTC | `c182ef6` |
| Meteosat-10 aerosol `a79` | NetCDF4 | 5,573,699 bytes | 2025-11-15 16:04:23 UTC | `c182ef6` |

Both directories contain a copied driver, `scatfile.sav`, `timestamp.txt`, and
`git_revision.txt`. The product metadata and stored schema are therefore
usable as references, while the exact historical generator command remains a
provenance question.

## Volcanic ash status

The repository contains volcanic-ash microphysics and RI inputs, including
`input_files/microphysics/volcanic-ash_ey1.mm` and `input_files/ri/eyj2010.ri`,
as well as ash LUT definitions. A volcanic-ash driver was not found in the
current driver inventory. The authorised external archive corrects the earlier
local-inventory conclusion: it contains current V2 Meteosat-12 FCI products
`meteosat-12_fci_m_volcanic-ash_a01_pctn_v21.nc`,
`meteosat-12_fci_m_volcanic-ash_a01_pey1_v21.nc`, and
`meteosat-12_fci_m_volcanic-ash_a01_pv81_v21.nc`. Ash remains a follow-on
reproduction target, not the first numerical baseline.

## Validation order for the selected cases

For Candidate A, compare in this order:

1. NetCDF dimensions and coordinate values;
2. instrument metadata, flags, SRF filenames, and solar/thermal channel lists;
3. microphysical arrays (`average_volume_per_particle`, extinction, extinction
   ratio, single-scattering albedo, asymmetry);
4. one monochromatic/SRF point before spectral integration;
5. one DISORT operator at representative optical depth, size, and geometry;
6. full SRF-normalised operators and percentiles/maxima of residuals.

Then repeat the same sequence for Candidate B, adding pressure and gas-level
comparisons. Do not use arbitrary tolerances to hide a first-intermediate
discrepancy; record whether it arises from a grid, units, phase normalisation,
Rayleigh/gas composition, DISORT input, SRF weighting, or NetCDF ordering.

## Provenance cautions

- Existing products have different dates, versions, output dimensions, and
  sometimes different driver-era conventions.
- The output directory is keyed by driver filename, while model and LUT names
  are command-line overrides.
- The nested `create_orac_lut` Git worktree has many pre-existing modifications,
  deletions, and untracked files. These candidate documents do not assert that
  every current file was part of the historical generation environment.
- The source named `create_orac_lut_v2.pro` currently includes two helper files
  deleted from the working tree. Treat the committed `e5758f6` snapshot and the
  existing NetCDF products as separate evidence layers until an executable
  legacy environment is recovered.
