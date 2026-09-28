# Grid B — current recommended LUT grid

Grid B is the definitive LUT sampling specification for new ORAC LUT generation. “Grid B” names the sampling grid only; it is independent of ORAC software, processing, generator, and NetCDF product version numbers.

Use these operational definitions:

| Particle family | Operational LUT definition |
|---|---|
| Liquid-water cloud | `create_orac_lut/input_files/lut/liquid-water-cloud-grid-b.lut` |
| Ice cloud | `create_orac_lut/input_files/lut/ice-cloud-grid-b.lut` |
| Generic tropospheric PA76 aerosol | `create_orac_lut/input_files/lut/aerosol-pa76-grid-b.lut` |

Grid B arose from the value-and-Jacobian convergence study. The three families retain their scientifically distinct tau, effective-radius, solar-zenith, satellite-zenith, and relative-azimuth node sets. All use one explicit satellite-zenith grid covering 0–75°, so smaller-angle instruments use the accessible portion and geostationary instruments are not stretched to their advertised 90° limit.

The physical evidence and exact source arrays remain in the Stage 18 JSON results under `validation/v23/results/`. The concise exact-grid summary and parser/preflight evidence are in `validation/v23/results/GRID_B.md` and `validation/v23/results/grid_b_manifest.json`.

Specialised aerosols, including volcanic ash, Saharan dust, and biomass-burning aerosol, are outside the current Grid B scope and may receive separate future grid studies.
