# Current V2 reference case

## Selected product

The initial Python reproduction target is:

```text
/network/group/aopp/eodg/RGG004_GRAINGER_ORACFILE/ORAC_LUTS/
  meteosat-10_seviri_m_liquid-water_a01_pstg_v21.nc
```

This is the user-specified current V2 archive product. It is NetCDF4,
6,640,898 bytes, and was readable at investigation time. Its SHA-256 is
`28b870b8f86dc6af9e591053313267111cf7367d7a851d271f9e274080a2228c`.
The repository-local fingerprint is at
`validation/reference_cases/meteosat10_seviri_liquid_water_stg_v21/fingerprint.json`.

## Current configuration

The product metadata and archived driver establish:

- platform: `meteosat-10`;
- instrument: `seviri`, instrument file `meteosat-10_seviri_v1.inst`;
- particle model: `liquid-water_stg.mm`;
- LUT definition: `liquid-water-cloud.lut`;
- atmosphere: code `2` (midlatitude summer / MODTRAN model);
- channels: 1–11;
- V2 output version: `21`;
- monochromatic/effective-centre mode: `srf_quad=1`;
- Rayleigh: enabled by default; no `no_rayleigh` override is recorded;
- gas: not enabled in the archived cloud driver or current cloud command;
- DISORT: 60 streams;
- output: `m_liquid-water_a01_pstg_v21.nc` in the cloud-driver output family.

The coordinate grids are 17 optical depths, 20 effective radii, 10 satellite
zenith angles, 10 solar zenith angles, and 11 relative azimuth angles. There is
no pressure dimension. The product has 11 channels, 4 solar channels, 1 mixed
channel, and 8 thermal channels.

## Current source path

For the current source tree’s cloud-family implementation, the applicable path
is:

```text
archived/current cloud driver
  -> create_orac_cloud_lut_wrapper
  -> create_orac_cloud_lut
  -> load_inststr/load_lutstr/load_srfstrarr/load_mmdat/load_atmstr
  -> generate_scattering_properties
  -> setup_disort/call_disort
  -> write_v2_lut
```

The retained current `meteosat-10_seviri_run` file is aerosol-only, so it does
not provide the exact cloud shell command used to create this archive product.
The archive contains the copied cloud driver and a `git_revision.txt` marker
recording revision `c182ef67d04967a8afeaeecf4a0157e9ee3b4850`. The reader and
reference description therefore capture the product and current writer schema
without inventing the missing shell provenance.

The current source also contains `makerunfile_v2.pro`, whose V2-era generated
files call `create_orac_lut_wrapper_v2`; that is a distinct source path and is
documented for comparison in `docs/legacy_lut_workflow.md`. The reproduction
target is the current cloud path that matches this product, not an older
version selected solely because its filename contains `v21`.

## Why this is first

It is a complete current V2 product with all SEVIRI channel classes and an
existing reference file. Its liquid-water model is a single modified-gamma Mie
component, making it a clean first test of grid handling, optical-property
loading, spectral metadata, DISORT operator ordering, and NetCDF I/O before
implementing POM or a DISORT backend.

Current V2 volcanic-ash products do exist in the external archive, including
Meteosat-12 FCI products with `pctn`, `pey1`, and `pv81` models. They are
follow-on cases because their particle optics are more complex and are not
needed to validate the first reader/data-layout milestone.

## Remaining uncertainties

- The exact historical/current cloud run script that created this product is
  not retained in the repository; the archived driver and output metadata are
  the available evidence.
- The archive marker points to an older Git revision than the current dirty
  working tree, so code provenance must be recorded separately from the current
  source path.
- Numerical reproduction of the scientific arrays remains future work; this
  task establishes only reading, structure, and reference fingerprints.

The fast development reference is defined separately in
`docs/legacy_test_lut_cases.md`. It uses the existing
`liquid-water-cloud_test.lut` grid with the current STG/cloud/Mie/DISORT path;
it is not a replacement for this full production reference.

## Python pipeline status

| stage | status |
|---|---|
| explicit reference configuration | implemented and validated against repository inputs |
| driver, instrument, LUT-grid, atmosphere, SRF, solar-spectrum, microphysics, and RI readers | implemented and validated with focused tests |
| V2 NetCDF reader | implemented and validated against the external reference |
| V2 NetCDF writer | implemented for supplied arrays; tested with a tiny synthetic file, not yet validated against a generated scientific product |
| legacy-equivalent Mie/scattering calculation | blocked at the current `create_bwgp` Mie/DLM boundary; Python interface defined, no substitute backend used |
| DISORT radiative transfer, SRF-integrated operator generation, and complete LUT generation | not yet implemented |
| POM backend | deliberately excluded from this phase |
