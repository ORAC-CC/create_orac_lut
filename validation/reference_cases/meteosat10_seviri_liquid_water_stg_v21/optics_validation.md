# STG optics validation status

## Target

The target is the current liquid-water STG particle-optics calculation used by
the selected Meteosat-10 SEVIRI V2 LUT at revision v21. The authoritative
source files are:

- `create_orac_lut/input_files/microphysics/liquid-water_stg.mm`;
- `create_orac_lut/input_files/ri/H2O_Segelstein_1981.ri`;
- `create_orac_lut/create_bwgp.pro`;
- `mie/mie_size_dist_new.pro`;
- `create_orac_lut/generate_scattering_properties.pro`;
- `create_orac_lut/input_files/lut/liquid-water-cloud.lut`;
- SEVIRI SRFs under `create_orac_lut/input_files/srf/`.

## Effective calculation

For each effective-radius and selected spectral point, the current code calls
`create_bwgp`. The STG model uses a single modified-gamma Mie component with
`Rm=12`, `S=0.1111111`, radius bounds `0.001`–`100.0` microns, RI from
`H2O_Segelstein_1981.ri`, and `xres=0.4` in the Mie size-distribution call.
The reference wavelength is 0.55 microns. The angular phase function is
returned as F11 on the `Dqv` cosine grid and is expanded into Legendre moments
by `legpexp` for downstream DISORT.

The Python readers currently recover and test the model parameters, RI table,
LUT effective-radius grid, atmosphere, SRF files, and solar spectrum. The
Python optical-property interface is defined in
`src/oraclut/optics/legacy.py`.

## DLM availability

The repository contains `mie/dlm-code/mie_dlm_single.so` and its source, but it
is an IDL DLM with unresolved IDL runtime symbols and unavailable Intel runtime
dependencies in the current environment. A read-only `ctypes` probe fails at
load time on `libifcore.so.5`. No standalone legacy Mie execution was possible,
and no external optical-property intermediate file was identified.

Accordingly, no floating-point comparison of extinction, SSA, phase functions,
or Legendre moments is claimed yet. The exact remaining boundary is loading or
calling the current ORAC Mie/DLM implementation with its required IDL/compiler
runtime. The existing POM-callable Mie implementation was not needed and was
not used in the production path or validation.

## Validation matrix

| quantity | expected source | status |
|---|---|---|
| wavelength/RI grid | `read_ri.pro`, SRF loader | input parsing implemented; numerical optics blocked |
| effective-radius grid | LUT definition | parsed and tested: 20 values, 1–39 microns |
| size distribution | `mie_size_dist_new.pro` | equations and parameters documented; Mie call blocked |
| extinction/scattering/SSA | Mie DLM | blocked |
| phase function F11 | Mie DLM + size weighting | blocked |
| Legendre moments | `legpexp` | blocked until phase function is available |
| reference-wavelength outputs | `create_bwgp` at 0.55 microns | blocked |

No DISORT or complete LUT generation was attempted.
