# Nakajima–King LUT comparison

`compare_nakajima_king.py` (this directory) compares the bidirectional reflectance (`R_0v`) in
two compatible ORAC V2 NetCDF LUTs. It produces four panels: LUT A, LUT B, a
solid/dashed overlay, and pointwise displacement vectors from A to B.

```bash
python validation/nakajima_king/compare_nakajima_king.py LUT_A.nc LUT_B.nc \
  --x-channel 1 --y-channel 3 \
  --solar-zenith 30 --satellite-zenith 0 --relative-azimuth 0 \
  --output validation/nakajima_king/results/example.png
```

Channel arguments are ORAC channel IDs. They are required when a LUT contains
more than two solar channels because no channel pair is universally valid for
a Nakajima–King retrieval. If exactly two solar channels exist, that unique
pair is selected automatically. Geometry arguments must match stored coordinate
values exactly; if omitted, the first value on each geometry axis is used.
Aerosol LUTs with a pressure dimension additionally accept `--surface-pressure`.

Blue curves hold optical depth constant and vary effective radius; their `τ`
labels therefore identify optical-depth isolines. Red curves hold effective
radius constant and vary optical depth; their `rₑ` labels give radius in
microns. Representative labels span the stored grids without labelling every
curve. In the overlay, LUT A is solid and LUT B is dashed. The final panel
compares corresponding `(τ, rₑ)` points, using colour for
`|ΔR| = sqrt(ΔR_x² + ΔR_y²)` and arrows for the displacement from A to B.

The two files must have identical instrument metadata, channels,
optical-depth/effective-radius grids, geometry grids, pressure grids (if any),
and `R_0v` dimensions. The tool deliberately does not interpolate. Current V2
files contain one microphysical phase/model per product rather than a phase
dimension, so phase is selected by choosing the corresponding water or ice LUT
files. Stored central wavelengths may differ by at most 0.0001 µm to accommodate
historical recomputation from an otherwise unchanged SRF; larger spectral
differences are rejected.

Example comparing existing EarthCARE MSI V21 and V22 products:

```bash
python validation/nakajima_king/compare_nakajima_king.py \
  /network/group/aopp/eodg/RGG004_GRAINGER_ORACFILE/ORAC_LUTS/earthcare_msi_m_liquid-water_a01_pstg_v21.nc \
  /network/group/aopp/eodg/RGG004_GRAINGER_ORACFILE/ORAC_LUTS/earthcare_msi_m_liquid-water_a01_pstg_v22.nc \
  --x-channel 1 --y-channel 3 \
  --output validation/nakajima_king/results/earthcare_msi_liquid_water_stg_v21_vs_v22.png
```
