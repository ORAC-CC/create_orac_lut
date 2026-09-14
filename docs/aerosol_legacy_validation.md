# Legacy aerosol reference and comparison

The Python aerosol generator produces the tiny Meteosat-10 SEVIRI `aerosol_a79`
LUT, but until an independently generated legacy IDL product exists there is
nothing to compare it against. This note describes how that reference is
produced and how the comparison is read.

The cloud formulation is already validated this way: the captured legacy
reference and the Python product agree to 1.3e-4 absolute on the
radiative-transfer operators, with identical structure.

## The case

| Setting | Value |
| --- | --- |
| platform / instrument | meteosat-10 / seviri (`meteosat-10_seviri_v1.inst`) |
| microphysics | `aerosol_a79.mm` (two components, `waf` 0.625 and `saf` 0.375, both Mie) |
| LUT grid | `aerosol_test.lut` |
| channel | 1 |
| atmosphere | MODTRAN code 2 (midlatitude summer) |
| SRF quadrature | 1 |
| gas / Rayleigh | enabled / enabled |
| DISORT streams / phase order | 60 / 1000 |
| grid | optical depth 2, effective radius 2, solar zenith 2, satellite zenith 2, relative azimuth 2, surface pressure 3 = 96 states |

Python product (already generated, not to be overwritten):
`create_orac_lut/luts/meteosat-10_seviri_aerosol/meteosat-10_seviri_m_aerosol_a12_pa79_v21.nc`.

## Which legacy source, and why it matters

The V2.1 tree imported under `create_orac_lut/` contains **no** aerosol entry
point: `create_orac_aerosol_lut.pro` and its wrapper were deleted in the V2.1
tip `b3c8a91`, and V2.1's `makerunfile.pro` calls a procedure
`create_orac_lut_aerosol` that exists nowhere in that tree. The only executable
aerosol path is the original working tree, extracted on scratch from the
verified archive
`…/oraclut-reconciliation/preservation/v21_preserved_current_local_create_orac_lut.tar.gz`:

```text
/network/scratch/grainger/oraclut-reconciliation/legacy-reference-source/current_local_create_orac_lut
```

Its aerosol entry points are the original working files, confirmed by
filesystem metadata:

| File | Size | Modified (UTC) |
| --- | --- | --- |
| `create_orac_aerosol_lut.pro` | 39,296 bytes | 2025-11-14 15:40:32 |
| `create_orac_aerosol_lut_wrapper.pro` | 10,223 bytes | 2025-11-14 15:40:30 |

The runner sets this root once, as `LEGACY_SRC_DEFAULT`; `ORACLUT_LEGACY_SRC`
overrides it. The tree stays on scratch — it is read as the internally
consistent legacy reference implementation and never copied into the
repository or modified. (An earlier copy under
`preservation/git_preservation/…/build/` is stale and no longer complete; the
runner's preflight rejects it.)

Those entry points cannot be run on top of V2.1 infrastructure. Source
comparison shows five routines in the aerosol call graph differ between the two
generations, and the differences are scientific, not cosmetic:

| Routine | Difference |
| --- | --- |
| `generate_scattering_properties` | V2.1 renames the `tmatrix_path` keyword to `tmatrix_dir` and drops the `force_n`/`force_k` blocks |
| `load_srfstrarr` | **V2.1 swaps the meaning of `srf_quad=0` and `srf_quad=1`**: under V2.1 `srf_quad=1` integrates over the SRF, under the working tree it is monochromatic at the effective channel centre — while the output filename still records `m` for mono |
| `load_lutstr` | V2.1 adds `Invalid_N_FLAG` to the surface-pressure `lut_quadrature` call ("latent bug corrected") — directly in the aerosol pressure path |
| `load_inststr` | V2.1 corrects `Count Eq -1` to `Count Eq 0` and the channel-id string/integer comparison |
| `load_atmstr` | V2.1 selects the old profile by `Atmospheres Eq 'midsatm.dat'`, the working tree by `Atmospheres Eq 0` |

`segment.pro`, needed for band mode, is absent from V2.1 altogether;
`create_bwgp.pro` and `read_ri.pro` also differ.

A trial run confirmed this empirically. IDL searches the **current working
directory before `IDL_PATH`**, so running from `create_orac_lut/` resolved the
V2.1 copies regardless of path order, and the run halted immediately:

```text
% Keyword TMATRIX_PATH not allowed in call to: GENERATE_SCATTERING_PROPERTIES
% Execution halted at: CREATE_ORAC_AEROSOL_LUT  352
```

The runner therefore executes IDL from a clean directory containing no `.pro`
files — only a symlink to the repository input hierarchy and a real `luts/`
output directory — with the extracted tree on `IDL_PATH`. Before any
calculation it runs `validation/diagnostics/aerosol_legacy_provenance.pro`,
records the source file of every resolved routine in the manifest, and
**refuses to continue** unless every required routine resolves from the one
selected source tree *and* no routine at all resolves from the active V2.1
tree. Resolve-only verification against the extracted tree confirms all 24
routines come from it, none from V2.1, none unresolved.

The preflight additionally requires all 27 legacy source files the aerosol path
needs — including `makerunfile_v2.pro`, the configuration mechanism the Python
generator replaces — so an incomplete or wrong tree fails before IDL starts.

This is also the generation of source that produced the validated cloud
reference and whose conventions the Python implementation reproduces, so the
comparison is like-for-like. The preserved source is only ever read.

## General runner

`scripts/generate_aerosol_legacy_reference.sh` is now a thin front end to
`scripts/generate_legacy_reference.sh`, which handles both formulations
(`--formulation cloud|aerosol`), any driver, explicit reference/comparison
directories and the same provenance gate. It is the runner used for the
visible + infrared matrix in `docs/validation_matrix.md`.

## Phase A: produce the reference (licensed IDL shell, user-run)

```shell
cd /home/g/grainger/project-oraclut
scripts/generate_aerosol_legacy_reference.sh
```

The runner loads `intel-compilers/2022` and `idl/890`, establishes the known
working environment (Mie and production DISORT2 DLMs only, `IDL_STARTUP`
unset), prints host, IDL version, both paths, every input file, the channel,
atmosphere, gas/Rayleigh state, SRF quadrature and the expected output, then:

1. fingerprints any pre-existing canonical product (SHA-256 and mtime);
2. runs `create_orac_aerosol_lut_wrapper` once;
3. verifies the canonical NetCDF exists and is genuinely new or changed;
4. copies it to
   `validation/generated/aerosol_legacy_reference/meteosat10_seviri_m_aerosol_a12_pa79_v21_legacy.nc`;
5. writes `capture_manifest.json` (timestamp, host, IDL version, legacy source
   tree, every input, all settings, source and captured paths, SHA-256, size,
   routine provenance), plus `stdout.log`, `stderr.log` and `provenance.log`;
6. performs **one** repeatability run in a fresh IDL process if the first run
   emitted a DISORT near-singular warning (`--repeat` forces it, `--no-repeat`
   suppresses it), comparing the two legacy products bitwise;
7. runs the Phase B comparison automatically.

An existing captured reference is never overwritten: move it aside first.

## Phase B: compare

Run automatically by the runner, or by hand:

```shell
PYTHONPATH=src /home/g/grainger/miniforge3/envs/science/bin/python \
    -m oraclut.validation.product_comparison \
    validation/generated/aerosol_legacy_reference/meteosat10_seviri_m_aerosol_a12_pa79_v21_legacy.nc \
    create_orac_lut/luts/meteosat-10_seviri_aerosol/meteosat-10_seviri_m_aerosol_a12_pa79_v21.nc \
    --formulation aerosol \
    --legacy-log validation/generated/aerosol_legacy_reference/stdout.log \
    --out-dir validation/generated/aerosol_comparison
```

It writes `comparison.json` and `comparison.md` and modifies neither product.

**Structural**: dimensions and their order, variable names, shapes, dtypes,
attributes, coordinate values, and fill/finite counts per variable. Result
`PASS` or `FAIL`.

**Numerical**: for every common variable — maximum absolute difference, RMS,
mean absolute difference, maximum relative difference, differing and compared
element counts, nonfinite counts, and the index of the largest discrepancy.
Relative differences are meaningless where the reference is near zero, so the
absolute columns are reported first.

**Aerosol-specific**: surface-pressure coordinate values (950, 1013, 1050), the
position of `surface_pressure` in the operator dimension order, the 96-state
count, channel selection, wavelength, and the configuration record (two
components in order `waf`, `saf`, gas, Rayleigh, atmosphere 2, `srf_quad=1`).
Component structure is verified from the microphysics input, because the V2
product stores only the mixed optical properties.

**Warnings**: `UPBEAM--SGECO says matrix near singular` occurrences are counted
and quoted. A warning accompanied by a complete product is recorded, not
treated as failure; a missing or nonfinite product is.

### Classification

`GREEN` / `AMBER` / `RED` is reported separately from the structural result and
only after the residuals themselves. The bands are stated in
`NUMERICAL_BANDS` and are specific to this comparison — calibrated on the
validated cloud case, where residuals peak near 1.3e-4 on operators of order
one — not a universal scientific tolerance:

- `GREEN`: worst operator residual ≤ 1e-3, no nonfinite or shape problems;
- `AMBER`: worst operator residual ≤ 1e-2;
- `RED`: worse than that, or a shape/order mismatch, nonfinite candidate
  values, or a systematic pattern (more than half the elements differing by
  more than 1e-3).

Fill-value differences are reported in their own section and never drive the
classification: the legacy writer is known to emit fill values for
singleton-channel optics variables where Python writes finite ones.

If legacy repeatability turns out not to be bitwise, the run-to-run envelope is
recorded and bitwise Python equality is not expected.

## If residuals appear

Diagnose the category before changing anything — and change nothing in the
scientific code as part of this validation: input parsing, microphysics,
refractive index, optical properties, pressure profile, SRF treatment, DISORT
inputs, DISORT numerical sensitivity, or NetCDF writing. The first useful
discriminator is whether the single-scattering optics agree: if they do, the
difference is downstream of Mie.
