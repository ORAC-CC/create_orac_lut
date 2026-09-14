# Legacy ORAC LUT inputs and outputs

## Configuration layers

The legacy calculation has four configuration layers:

| layer | examples | controls |
|---|---|---|
| run generator | `makerunfile_v2.pro`, `makerunfile.pro` | platform/instrument, model list, forward-model label, version, command-line overrides |
| driver | `input_files/driver/*.driver` | input root, instrument file, default model/LUT file, atmosphere code, optional channels/flags |
| instrument | `input_files/inst/*.inst` | channel IDs, SRFs, solar/thermal/mixed flags, calibration/noise metadata, viewing limit |
| science grids/profiles | `input_files/microphysics/*.mm`, `input_files/lut/*.lut`, atmosphere/gas/SRF/RI/optics trees | particle model, size/profile, optical-depth and geometry coordinates, atmospheric and spectral data |

The driver is text-parsed by the wrapper. The first five non-comment logical
values are the input root, instrument filename, microphysical filename, LUT
filename, and atmosphere selector. Optional assignments follow. Generated IDL
commands override several of those values.

## Instrument inputs

`load_inststr.pro` reads an instrument definition from
`input_files/inst/`. It supplies at least:

- platform, instrument, and instrument version;
- available channel IDs and selected `channelid` subset;
- SRF filenames;
- central wavelength/wavenumber and calibration/noise fields;
- solar, thermal, and mixed-channel flags;
- maximum satellite zenith angle and view configuration.

Examples inspected:

- Meteosat-10 SEVIRI: 11 channels, channels 1–4 solar, 4 mixed, 5–11
  thermal, maximum satellite zenith 90 degrees;
- Envisat AATSR: 7 channels, dual-view configuration, solar channels 1–5,
  thermal channels 5–7, maximum satellite zenith 56 degrees;
- EarthCARE MSI: 7 channels, solar channels 1–4 and thermal channels 5–7,
  maximum satellite zenith 20 degrees.

Dual-view instruments cause channel/SRF state to be duplicated for the forward
view. This must be retained in validation because it changes output shape and
view semantics.

## LUT grid inputs

`load_lutstr.pro` reads the optical-depth, effective-radius, solar-zenith,
satellite-zenith, and relative-azimuth coordinate pairs. The current split
aerosol path can also read a surface-pressure coordinate. The parser expands
the ranges according to the specified grid conventions; do not replace these
with guessed modern grids.

Observed definitions include:

| file | grid evidence |
|---|---|
| `input_files/lut/liquid-water-cloud.lut` | 17 optical depths from approximately `1e-10` to `256`, effective-radius range `1`–`39` with four entries, 10 solar zenith, 10 satellite zenith, 11 relative azimuth |
| `input_files/lut/ash-cloud.lut` | 17 optical depths, 13 explicit effective radii including `0.1`–`15`, 10/10/11 angular grid |
| `input_files/lut/aerosol.lut` | 20 uneven optical depths approximately `0.0078125`–`5.656854`, 20 effective-radius values approximately `0.01`–`10`, 10/10/11 angular grid, and a pressure entry in the current definition |

The exact values and ordering in the files are the validation truth. The
operators use IDL array ordering inherited from these loaders; dimensions in
NetCDF are not evidence that the internal loop order may be changed freely.

## Microphysical model inputs

`load_mmdat.pro` parses `input_files/microphysics/*.mm`. A model describes a
substance and short name, a vertical effective-radius/relative-profile table,
one or more components, size-distribution parameters, scattering code, mixing
ratio, refractive-index file, and optionally aspect-ratio data.

Representative models:

| model | particle-property path |
|---|---|
| `liquid-water_stg.mm` | one user component, modified-gamma size distribution, Mie, `H2O_Segelstein_1981.ri`; profile concentrated at the 3.5 km entry |
| `aerosol_a79.mm` | two user components (`waf`/`saf`), lognormal radius parameters `0.070`, `1.700`, Mie, mixing ratios `0.625`/`0.375` |
| `aerosol_a70.mm` | two Mie components plus one T-matrix `mdc` component with aspect-ratio data and corresponding refractive index |
| `volcanic-ash_ey1.mm` | user `ash_eyja`, lognormal parameters approximately `0.217`, `1.77`, Mie, `eyj2010.ri` |
| `water-ice_agg.mm` | Baum aggregate-solid-column optical properties from a NetCDF full-phase-matrix file |

The parser resolves refractive-index names below `input_files/ri/`. The
available tree includes water, ice, sulphuric-acid, biomass, aerosol, and
volcanic-ash RI files. The model name does not itself select the aerosol or
cloud top-level equations.

## Particle optical-property boundary

`generate_scattering_properties.pro` turns the microphysical description into
the downstream RT representation. The important quantities are:

- volume/average-particle quantities over effective radius;
- volume extinction at each instrument spectral point;
- extinction ratio relative to the 550 nm reference;
- single-scattering albedo;
- asymmetry parameter;
- phase function and Legendre moments.

The legacy component backends include user Mie, user T-matrix, Baran, and Baum
optical properties. Mie/T-matrix calculations may use external compiled code;
Baran/Baum readers load supplied optical-property datasets. A future POM
adapter should produce these scientific quantities and their metadata without
coupling POM internals to the RT solver.

## Spectral inputs and conventions

`load_srfstrarr.pro` resolves SRFs below `input_files/srf/` and loads
`input_files/sun/Gueymard2018.sssi`. It supports:

- `srf_quad=0`: raw SRF behavior;
- `srf_quad=1`: effective central-wavelength/monochromatic behavior;
- `srf_quad=2`: SRF quadrature generated by IDL segmentation.

The inspected V2 generator activates only `srf_quad=1`, despite containing a
band branch. The loader computes solar in-band quantities with Simpson
integration, using wavelength or wavenumber treatment according to spectral
region. Thermal channels use blackbody constants from `bbconstants.pro` and
Planck averaging in the RT path. This is why a center-wavelength port must not
be generalized to all instruments without first matching the selected legacy
path.

## Atmospheric and gas inputs

`load_atmstr.pro` supports the listed simple three-column midlatitude profile
and MODTRAN-style profiles. It retains the relevant atmospheric heights and
constructs RT layers from adjacent profile levels. The atmosphere selector is
the fifth mandatory driver value.

Gas tables follow the naming pattern:

```text
ModtranGasOpd_A<atmospheres>_<platform>_<instrument>_ch<two-digit>.gas
```

`load_gasstr.pro` reads per-level gas optical depth and checks that gas heights
match the selected atmosphere. Current aerosol run scripts explicitly pass
`atmospheres=2,gas=1`; current cloud scripts generally do not pass gas.

Rayleigh is combined in the RT layer state. Historical V2 uses a surface
pressure reference of 1013 hPa; the current aerosol split implementation adds
pressure-dependent output and scales the Rayleigh term by the LUT pressure.

## Radiative-transfer inputs and operators

`setup_disort.pro` configures the DISORT common block with 60 streams, the
atmospheric layer count, angular counts, and phase-moment count. It sets two
user optical-depth boundaries, user angles, Lambertian surface defaults, and
the legacy print/accuracy settings. `call_disort.pro` prepares nonthermal or
Planck thermal calls, compresses identical adjacent layer properties, and
invokes the DISORT IDL DLM.

The V2 operator arrays are channel-first internally. Without pressure their
scientific shapes are equivalent to:

```text
T_dv  [channel, effective_radius, optical_depth, satellite_zenith]
T_dd  [channel, effective_radius, optical_depth]
R_dv  [channel, effective_radius, optical_depth, satellite_zenith]
R_dd  [channel, effective_radius, optical_depth]
R_0v  [solar_channel, effective_radius, optical_depth,
       solar_zenith, satellite_zenith, relative_azimuth]
R_0d  [solar_channel, effective_radius, optical_depth, solar_zenith]
T_0d  [solar_channel, effective_radius, optical_depth, solar_zenith]
T_00  [solar_channel, effective_radius, optical_depth, solar_zenith]
E_md  [thermal_channel, effective_radius, optical_depth, satellite_zenith]
```

The pressure-aware writer inserts `surface_pressure` between the
`effective_radius` and `channels` dimensions. Existing NetCDF declarations in `write_v2_lut.pro` are the exact
schema evidence. The corresponding meanings are diffuse/direct transmission,
diffuse/direct reflection, bidirectional reflectance, diffuse reflectance of a
direct beam, direct/diffuse transmission, and diffuse emissivity.

Solar operators are normalised by SRF weights and include the legacy direct-beam
handling. Thermal operators are produced for thermal channels; mixed channels
are represented in instrument metadata and the relevant solar/thermal arrays.

## Output files and provenance

The wrapper creates `luts/<driver stem>/`. `write_v2_lut` creates a NetCDF4
file with dimensions for the grids, channel metadata, microphysical optical
properties, and RT operators. A typical filename is:

```text
<platform>_<instrument>_<m|b>_<substance>_a<atmosphere-code>_p<shortname>_v<two-digit-version>.nc
```

Examples observed in the repository:

- `luts/meteosat-10_seviri_cloud/meteosat-10_seviri_m_liquid-water_a01_pstg_v21.nc`
- `luts/meteosat-10_seviri_aerosol/meteosat-10_seviri_m_aerosol_a12_pa79_v21.nc`
- `luts/envisat_aatsr_cloud/envisat_aatsr_m_liquid-water_a01_pold_v50.nc`
- `luts/himawari-8_ahi_aerosol/himawari-8_ahi_m_aerosol_a12_pa70_v14.nc`

The V2 generator also writes/copies the driver and saves scattering state. The
copied driver is necessary provenance, but it may not record command-line
overrides such as `mmfile`, `lutfile`, `gas`, `version`, or `tmatrix_path`.

## External/runtime dependencies

Required or expected dependencies include:

- IDL runtime and IDL NetCDF routines;
- the DISORT DLM/Fortran implementation;
- external Mie and possibly T-matrix code;
- the Dubovik T-matrix optical-property directory;
- supplied RI, Baran/Baum, SRF, solar, atmosphere, and gas files;
- host-specific IDL/DLM paths and, for newer scripts, Intel compiler modules.

The verified user-terminal environment loads `intel-compilers/2022` and
`idl/890`, with IDL 8.9.0 and the Oxford licence. The current production
cloud path uses `setup_disort`/`call_disort` and the repository
`create_orac_lut/disort2/src` DLM; DISORT4 is not substituted. Codex's
isolated shell may locate IDL but cannot currently acquire the licence, so the
explicit fast test runner is intended for the user's already working terminal.
