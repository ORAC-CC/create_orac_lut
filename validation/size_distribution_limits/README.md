# Particle-size integration-limit validation

`run_study.py` is a standalone validation driver for the radius-domain
integration used by ORAC's legacy-equivalent Mie path. It calls the preserved
Mie kernel but does not alter production code.

From the repository root:

```bash
MPLCONFIGDIR=validation/tmp/matplotlib-cache \
PYTHONPATH=src /home/g/grainger/miniforge3/envs/science/bin/python \
  validation/size_distribution_limits/run_study.py --phase-order 128
```

The default study covers the six requested liquid-cloud effective radii,
three established repository aerosol modes, six MODIS/SLSTR-relevant
wavelengths, all requested upper and lower limits, and each aerosol mode's
actual legacy adaptive bounds. It records both literal recomputed production
grids and denser common-grid tail differences, because changing a linear-grid
endpoint also moves every quadrature node. Results are written beneath
`results/` as CSV, JSON, and PNG files.

The 128-angle default resolves ORAC moments 0–127. It is a convergence study,
not a replacement for the production default of 1000 moments. A targeted
higher-order run can be made with `--phase-order 1000`, but it is substantially
more expensive and should first be restricted in the source case matrix or run
using the repository's execution policy.

The supplied targeted production-order check is:

```bash
PYTHONPATH=src /home/g/grainger/miniforge3/envs/science/bin/python \
  validation/size_distribution_limits/run_high_order_check.py
```

Scientific interpretation, source tracing, thresholds, and limitations are in
`REPORT_size_distribution_limits.md`.  That report is the exploratory record
from before the decision.

## Validation of the adopted limit (2026-10)

Production now integrates liquid-water modified-gamma distributions to the
first node of the legacy radius lattice at or beyond 3.5 x effective radius
(commit `c6ad545`). The supporting evidence is summarised in
`validation/REPORT_lut_numerics_development.md` and produced by:

| File | Purpose |
|---|---|
| `quadrature_endpoint_diagnostic.py` | Separates the removed tail from node-placement error. It compares reference, naive candidate, endpoint-jittered and lattice integrations against a converged truth (xres 0.05, 0.001–150 µm). `--factor` sets the multiple of r_e. |
| `compact_lut_comparison.py` | Compact MODIS and dual-view SLSTR LUTs through `create_orac_cloud_lut`. Variants: reference (0.001–100 µm), naive candidate, lattice candidate (run-time substitution) and adopted (production code). Includes retrieval-level interpolation and Jacobian differences normalised by instrument noise. |
| `compact_liquid_grid.lut` | Compact grid used by the comparisons (not a production grid). |
| `benchmark_integration_limits.py` | Single-thread Mie-step cost over the Grid B liquid radii. |

Results are in `results/quadrature_endpoint*/`, `results/compact_lut/<variant>/`
and `results/integration_limit_benchmark.csv`. The `candidate` and `lattice`
compact results and `results/quadrature_endpoint/` use the exploratory factor
3. `results/quadrature_endpoint_k3.5/`, `results/quadrature_endpoint_small_radii_k3.5/`
and `results/compact_lut/adopted/` use the adopted 3.5.

```bash
PYTHONPATH=src nice -n 19 /home/g/grainger/miniforge3/envs/science/bin/python \
  validation/size_distribution_limits/quadrature_endpoint_diagnostic.py --factor 3.5 \
  --results validation/size_distribution_limits/results/quadrature_endpoint_k3.5
nice -n 19 /home/g/grainger/miniforge3/envs/science/bin/python \
  validation/size_distribution_limits/compact_lut_comparison.py generate modis adopted
nice -n 19 /home/g/grainger/miniforge3/envs/science/bin/python \
  validation/size_distribution_limits/compact_lut_comparison.py compare adopted
PYTHONPATH=src /home/g/grainger/miniforge3/envs/science/bin/python \
  validation/size_distribution_limits/benchmark_integration_limits.py
```

Since `c6ad545` the `reference`, `candidate` and `lattice` variants of
`compact_lut_comparison.py` re-impose the legacy 0.001–100 µm limit first, so
they still start from the V23 integration.
