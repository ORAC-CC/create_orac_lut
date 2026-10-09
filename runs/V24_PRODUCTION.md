# V24 ORAC cloud LUT production

## What V24 is

V24 is the V23 cloud LUT product set, recomputed with the revised Mie
numerics.

**Unchanged from V23:**
- platforms, instruments, channels and microphysical models;
- forward model, atmosphere, gas and Rayleigh settings;
- the Grid B product sampling: `liquid-water-cloud-grid-b.lut` and
  `ice-cloud-grid-b.lut`.

**Numerical changes since the source that produced V23 (`9d663e9`):**

| Commit | Change | Applies to |
|---|---|---|
| `c6ad545` | Size integration stops at the legacy-lattice node ≥ 3.5 r_e | Liquid-water modified gamma |
| `dbcc42c` | Adaptive Legendre expansion: Nq, L and the DISORT moments come from each averaged phase function; `nmom` is obsolete | Mie classes |
| `96bebbe` | Nested refinement of the radius trapezoid to a size-parameter step ≤ 0.025 | Liquid water |
| `96bebbe` | The same refinement, to a step ≤ 0.05 | Ice spheres |

These are described in `create_orac_luts.py` ("Versions and numerical
changes") and validated in `validation/REPORT_lut_numerics_development.md`.
The Baum ice classes (agg, ghm, src) are computed exactly as in V23
(`nmom = 1000`).

## Products

There are 40 products:
- Aqua and Terra MODIS (channels 1–36), and Sentinel-3A and -3B SLSTR (channels
  1–9, dual view);
- each with six liquid-water models (240, 253, 263, 273, old, stg) and four ice
  models (sph, agg, ghm, src).

Run files are `runs/<platform>_<instrument>_cloud_<model>_v24.run`, and products
are written to

    /network/group/aopp/eodg/RGG004_GRAINGER_ORACFILE/ORAC_LUTS/<platform>_<instrument>_m_<substance>_a01_p<shortname>_v24.nc

next to, and never over, the `_v23.nc` products. The generator refuses to
overwrite an existing product.

## Submission and provenance

    scripts/submit_v24_cloud_luts.sh --dry-run      # check all 40 run files; submit nothing
    scripts/submit_v24_cloud_luts.sh                # submit every product not yet present

The script submits nothing unless:
- the production paths are unmodified, and
- HEAD is contained in `origin/main` (it fetches first).

Each submission is appended to `validation/v24/production_submissions.tsv`,
which records time, host, Git revision, SLURM job ID, run file and product.
Each job log (`validation/slurm/v24_*.out`) records the revision and the state
of the source tree when the job started.

The jobs execute this shared working tree. **Do not change production source
until every job has started.**

## NetCDF text-attribute repair (2026-10-09)

The V24 products as written stored every text attribute as NC_CHAR, which
ORAC's `nc_get_att_string` read of the axis `spacing` attributes rejects.  All
40 V24 products were repaired, metadata only, with
`python -m oraclut.repair_string_attributes` (data verified byte-identical;
originals retained as hard links in
`ORAC_LUTS/originals_before_nc_string_repair_20261009/`).  See
`LUT_FORMAT_CONTRACT.md` and
`validation/netcdf_string_attributes/REPORT_netcdf_string_attribute_repair.md`.
