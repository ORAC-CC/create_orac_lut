# Adaptive Legendre expansion of Mie phase functions

Validation material for the adaptive Legendre expansion
(`src/oraclut/idl_mirror/legendre_expansion.py`, used by
`generate_scattering_properties`). Methodology and results are summarised in
`validation/REPORT_lut_numerics_development.md`.

| File | Purpose |
|---|---|
| `REPORT_legendre_expansion_length.md` | Investigation record from before the implementation (Grainger 1990 King criterion, quadrature requirements, why NMom = 1000 was inadequate) |
| `adaptive_legendre_validation.py` | Coefficient-level validation: production expansion against a long reference and against the directly calculated averaged phase function at about 1400 dense angles |
| `compact_lut_legendre.py` | Compact MODIS, dual-view SLSTR and SEVIRI aerosol LUTs: adaptive against long reference and against fixed nmom = 1000 (revision c6ad545) |
| `compact_aerosol_grid.lut` | Compact aerosol grid used by `compact_lut_legendre.py` (not a production grid) |
| `benchmark_legendre.py` | Single-thread scattering-step cost, fixed against adaptive, over the Grid B liquid radii |
| `results/` | CSV and JSON results of the three scripts |

From the repository root (all of these exceed 30 minutes in total; use `nice -n 19`):

```bash
PYTHONPATH=src nice -n 19 /home/g/grainger/miniforge3/envs/science/bin/python \
  validation/legendre_expansion/adaptive_legendre_validation.py
OMP_NUM_THREADS=1 nice -n 19 /home/g/grainger/miniforge3/envs/science/bin/python \
  validation/legendre_expansion/compact_lut_legendre.py all
nice -n 19 /home/g/grainger/miniforge3/envs/science/bin/python \
  validation/legendre_expansion/benchmark_legendre.py
```

`compact_lut_legendre.py` and `benchmark_legendre.py` export revision `c6ad545`
into `validation/tmp/rev_c6ad545` (with links to `build/`, `create_orac_lut/` and
`mie/`) to run the fixed-1000 expansion. Products and logs go under the ignored
`validation/tmp/`.
