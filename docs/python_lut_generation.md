# How do I generate an ORAC LUT in Python?

`python -m oraclut.generate` is the single entry point for both legacy
formulations. It replaces the practical purpose of `makerunfile_v2.pro` and
its generated IDL run scripts: you edit one configuration file, and everything
else — instrument definition, SRFs, refractive indices, output naming — is
resolved from the normal `create_orac_lut/input_files/` hierarchy.

## Normal workflow

Four steps. You should not need to edit Python source, the SLURM batch script,
or long command lines.

**1. Copy the template and edit it.**

```shell
cd ~/project-oraclut
cp configs/template.driver configs/my_new_lut.driver
$EDITOR configs/my_new_lut.driver
```

`configs/template.driver` is commented setting by setting. In practice you
edit five things: the instrument definition, the microphysical model, the LUT
grid, the atmosphere code, and the channel list.

**2. Dry-run it (fast, interactive, no calculation).**

```shell
PYTHONPATH=src /home/g/grainger/miniforge3/envs/science/bin/python \
    -m oraclut.generate --config configs/my_new_lut.driver --dry-run
```

Read the printed summary: dimensions, grid-point counts, estimated DISORT
workload, and the exact output file — and whether it already exists.

**3. Log into the submission node.** SLURM submission is only from
`atmlxint7`.

```shell
ssh atmlxint7
cd ~/project-oraclut
```

**4. Submit.**

```shell
scripts/submit_oraclut_lut.sh --config configs/my_new_lut.driver
```

Then `squeue -u "$USER"`, and read `validation/slurm/oraclut_lut_<jobid>.out`.
The full SLURM details are in "Recommended route" below.

For interactive one-off commands the same options work without a
configuration file:

```shell
cd /home/g/grainger/project-oraclut
export PYTHONPATH=src
PYTHON=/home/g/grainger/miniforge3/envs/science/bin/python
```

## Shipped configurations

| File | Formulation | Grid | Size |
| --- | --- | --- | --- |
| `configs/template.driver` | cloud | `liquid-water-cloud_test.lut` | copy this and edit it |
| `configs/meteosat-10_seviri_cloud_liquid-water_stg_test.driver` | cloud | `liquid-water-cloud_test.lut` | 48 grid points, ≈ 30 s |
| `configs/meteosat-10_seviri_aerosol_a79_test.driver` | aerosol | `aerosol_test.lut` | 96 grid points, ≈ 15 s |
| `configs/meteosat-10_seviri_cloud_liquid-water_stg.driver` | cloud | `liquid-water-cloud.lut` | 374,000 grid points × 11 channels — SLURM |
| `configs/meteosat-10_seviri_aerosol_a79.driver` | aerosol | `aerosol.lut` | 440,000 grid points × 11 channels — SLURM |

The two `_test` configurations are the validated development cases; the two
production configurations reproduce the archived V2 product settings.

## Choosing the LUT grid: test versus full

The LUT grid file decides the size of the calculation, so **it is always named
explicitly and there is no hidden default**. In a configuration it is the
fourth positional value; on the command line it is `--lut`. Nothing selects a
production grid implicitly:

```shell
# error: names both choices instead of silently running a large grid
$PYTHON -m oraclut.generate --platform meteosat-10 --instrument seviri \
    --microphysics liquid-water_stg.mm --dry-run
# -> No LUT grid selected, and there is no default: ... use
#    --lut liquid-water-cloud.lut for the full production grid or --test for
#    the compact liquid-water-cloud_test.lut development grid
```

- **Development/test grid:** name `<family>_test.lut` in the configuration, or
  pass `--test` on the command line (the equivalent of `Test = 1` in
  `makerunfile_v2.pro`). Tens of grid points; seconds.
- **Full production grid:** name `<family>.lut`. Hundreds of thousands of grid
  points; hours, so submit it through SLURM.

The family follows the material of the microphysics name exactly as
`makerunfile_v2.pro` did: `liquid-water_*` → `liquid-water-cloud`,
`water-ice_*` → `ice-cloud`, `volcanic-ash_*` → `ash-plume`, `biomass_*` →
`biomass-plume`, `sulphuric-acid_*` → `sulphuric-acid-cloud`, `aerosol_*` →
`aerosol`. `--test` cannot be combined with an explicit grid. The dry-run
summary labels the resolved grid `[full grid]` or
`[compact development/test grid]`, so which one will be generated is never in
doubt.

As in the legacy code the output file name records the instrument, particle
model, atmosphere and revision but **not** the grid, so a test run and a full
run of the same case write the same file name. Set `output=` in the
configuration for experiments, or pass `--overwrite` deliberately; the
dry-run summary always reports whether the output already exists.

## Simplest cloud command

Discrete particle layer (legacy "cloud" formulation), Meteosat-10 SEVIRI,
liquid-water STG, compact test grid, channel 1:

```shell
$PYTHON -m oraclut.generate --config configs/meteosat-10_seviri_cloud_liquid-water_stg_test.driver
# or, without a configuration file:
$PYTHON -m oraclut.generate --platform meteosat-10 --instrument seviri \
  --microphysics liquid-water_stg.mm --test --channels 1
```

Output: `create_orac_lut/luts/meteosat-10_seviri_cloud/meteosat-10_seviri_m_liquid-water_a01_pstg_v21.nc`
(≈ 30 s on an AOPP interactive node).

## Simplest aerosol command

Particles distributed through the column with a surface-pressure dimension
(legacy "aerosol" formulation), Meteosat-10 SEVIRI, `aerosol_a79`, compact
test grid, channel 1, with gas absorption as in the operational run file:

```shell
$PYTHON -m oraclut.generate --config configs/meteosat-10_seviri_aerosol_a79_test.driver
# or, without a configuration file:
$PYTHON -m oraclut.generate --platform meteosat-10 --instrument seviri \
  --microphysics aerosol_a79.mm --test --gas --channels 1
```

Output: `create_orac_lut/luts/meteosat-10_seviri_aerosol/meteosat-10_seviri_m_aerosol_a12_pa79_v21.nc`
(≈ 15 s).

## Dry run

`--dry-run` resolves and prints the complete configuration without running Mie
or DISORT, and never writes anything:

```shell
$PYTHON -m oraclut.generate --config configs/meteosat-10_seviri_cloud_liquid-water_stg_test.driver --dry-run
```

```text
Resolved ORAC LUT configuration
  configuration file:  configs/meteosat-10_seviri_cloud_liquid-water_stg_test.driver
  forward model:       cloud (discrete particle layer)
  platform/instrument: meteosat-10 / seviri (instrument version unknown)
  instrument file:     create_orac_lut/input_files/inst/meteosat-10_seviri_v1.inst
  microphysics:        create_orac_lut/input_files/microphysics/liquid-water_stg.mm
                       liquid-water / stg, 1 component(s), scattering code mie
  LUT grid definition: create_orac_lut/input_files/lut/liquid-water-cloud_test.lut  [compact development/test grid]
  channels:            1 of 11 available: 1
                       solar [1]; thermal none
  atmosphere:          2 (midlatitude summer)
  Rayleigh scattering: on; gas absorption: off
  SRF treatment:       srf_quad=1 (monochromatic at the effective channel centre)
  DISORT streams:      60; phase-function order: 1000
  product:             LUT level 2, revision 21
  LUT dimensions:
      channels             1
      optical_depth        2
      effective_radius     3
      solar_zenith         2
      satellite_zenith     2
      relative_azimuth     2
      surface_pressure     - (cloud formulation: no pressure dimension)
      grid points/channel  48
      grid points total    48
  estimated workload:  18 DISORT states before angular/direct calls
  output:              create_orac_lut/luts/meteosat-10_seviri_cloud/meteosat-10_seviri_m_liquid-water_a01_pstg_v21.nc
  output status:       does not exist yet
dry-run: no Mie or DISORT calculation was executed
```

For the aerosol formulation the `surface_pressure` line carries its real size
(3 for `aerosol_test.lut`, 1 for the production `aerosol.lut`). A full grid is
equally cheap to dry-run and shows what a real run would cost, for example
`liquid-water-cloud.lut` with all 11 SEVIRI channels: 374,000 grid points per
channel, 4,114,000 in total, 153,340 DISORT states. If the output file already
exists the summary says so and a warning is printed, but the dry run still
succeeds.

Every missing or misspelled input fails immediately with a message naming
nearby candidates; nothing is silently substituted.

## Configuration-file example

`--config` accepts a readable driver in the legacy driver-file format: five
positional values (input root, instrument file, microphysics file, LUT file,
MODTRAN atmosphere code) followed by `key=value` options using the legacy
wrapper keyword names. `configs/template.driver` documents every setting; the
shipped examples are listed above. For instance:

```text
# configs/meteosat-10_seviri_aerosol_a79_test.driver
input_files                     # root input hierarchy (create_orac_lut/input_files)
meteosat-10_seviri_v1.inst      # instrument definition (input_files/inst)
aerosol_a79.mm                  # microphysical model (input_files/microphysics)
aerosol_test.lut                # LUT definition (input_files/lut)
2                               # MODTRAN atmosphere code (Midlatitude Summer)
forward_model=aerosol
channelid=[1]
srf_quad=1
gas=1
version=21
```

```shell
$PYTHON -m oraclut.generate --config configs/meteosat-10_seviri_aerosol_a79_test.driver
```

Settings combine as *defaults < configuration file < explicit command line*, so
`--channels 2,3 --no-gas --version 22` on the command line override the file.
`default.mm` / `default.lut` placeholders (as in the legacy per-instrument
drivers) leave that value to the command line; a configuration that leaves the
LUT grid as `default.lut` and is then run without `--lut` or `--test` fails
rather than choosing a grid for you.

Recognised keys: `channelid`, `srf_quad`, `gas`, `no_rayleigh`/`rayleigh`,
`version`, `forward_model`, `platform`, `instrument`, `output`, `streams`,
`phase_order`. Legacy science keywords that are not implemented (`force_n`,
`force_k`, `n_theta`, `mie`, `srfdat`, `tmatrix_path`, `opt_prop_luts`) are
rejected explicitly; workflow-only keywords (`no_screen`, `reuse_scat`,
`scat_only`) are accepted and ignored with a note.

## Relationship to makerunfile_v2

| `makerunfile_v2.pro` setting | Generated legacy argument | Python equivalent |
| --- | --- | --- |
| `platform`, `instrument` | driver `input_files/driver/<platform>_<instrument>_<FM>.driver`, `instfile` from driver | `--platform`, `--instrument` → `input_files/inst/<platform>_<instrument>_v*.inst` (or `--instrument-file`) |
| `MM` entry, e.g. `'liquid-water_stg'` | `mmfile='liquid-water_stg.mm'` | `--microphysics liquid-water_stg.mm` |
| material (text before `_`) → `FM` | `create_orac_cloud_lut_wrapper` or `create_orac_aerosol_lut_wrapper` | `--forward-model` (default derived from the material with the same table: aerosol→aerosol, others→cloud) |
| material → `lutfile` family | `lutfile='liquid-water-cloud.lut'`, `'aerosol.lut'`, `'ice-cloud.lut'`, … | the fourth configuration value, or `--lut`; always explicit, never defaulted |
| `Test = 1` | `lutfile='..._test.lut'` | `--test` (the only shortcut, and it selects the *test* grid) |
| `srf_quad` loop (1 or 2) | `srf_quad=1` (monochromatic) / `2` (band) | `--srf-quad 1|2` |
| driver `channelid=[...]` | `channelID` keyword | `--channels 1,2,3` or `all` (default all) |
| driver atmosphere line / `atmospheres=2` | `atmospheres` argument | `--atmosphere 0..6` (default 2) |
| `gas=1` in run file | `gas` keyword | `--gas` / `--no-gas` |
| `no_rayleigh=1` second call for liquid water | `no_rayleigh` keyword | `--no-rayleigh` (run a second command) |
| `Versions = '21'` | `version=21` | `--version 21` |
| `out_path = 'luts/<driver stem>'` + legacy filename | `<plat>_<inst>_<m|b>_<substance>_a<code>_p<shortname>_v<vv>.nc` | same default directory and filename; `--output` overrides (file or directory) |
| `reuse_scat`, `tmatrix_path`, `source initpath.bash`, `module load` | shell/IDL environment | not needed; no IDL, no shell setup |

The atmospheric-model code in the filename follows the legacy rule: `1X` with
gas absorption for MODTRAN atmosphere `X`, `00` without Rayleigh, `01` with
Rayleigh and no gas.

## Input-file hierarchy

All scientific inputs are the legacy files under `create_orac_lut/input_files/`:

| Directory | Contents |
| --- | --- |
| `inst/` | instrument definitions `<platform>_<instrument>_v1.inst` (24 instruments) |
| `srf/` | spectral response functions referenced by the `.inst` files |
| `lut/` | LUT grid definitions (`*.lut`, `*_test.lut` compact grids) |
| `microphysics/` | particle models (`*.mm`) |
| `ri/` | refractive indices referenced by the `.mm` components |
| `atm/` | MODTRAN atmosphere profiles (codes 1–6) and `midsatm.dat` (code 0) |
| `gas/` | MODTRAN gas optical-depth profiles per instrument channel and atmosphere |
| `sun/` | solar spectrum (`Gueymard2018.sssi`) |
| `baum/` | Baum ice optical-property tables (ice-cloud cases; not yet ported) |
| `driver/` | legacy per-instrument driver files (reference only) |

Python configurations live separately in `configs/`; generated products go to
`create_orac_lut/luts/` (Git-ignored) or wherever `--output` points inside the
repository.

## Output location

Default: `create_orac_lut/luts/<platform>_<instrument>_<formulation>/<legacy filename>`.
`--output` may be a `.nc` path or a directory (which receives the legacy-named
file). Existing files are never overwritten without `--overwrite`. Outputs are
NetCDF V2 layout, product revision from `--version` (default 21).

## Recommended route: SLURM submission from atmlxint7

Substantial ORAC LUT generation (anything beyond the compact `*_test.lut`
grids) should run under SLURM. **SLURM submission is only from `atmlxint7`.**
The repository provides a submission wrapper that does the fiddly parts —
host check, log directory, submission environment — so the normal operational
command is just a configuration file:

```shell
# 1. log into the submission node (the wrapper never does this for you)
ssh atmlxint7
# 2. go to the repository
cd ~/project-oraclut
# 3. choose a configuration in configs/, or copy configs/template.driver and edit it
# 4. optionally dry-run interactively to see the resolved inputs and RT-state count
PYTHONPATH=src /home/g/grainger/miniforge3/envs/science/bin/python -m oraclut.generate \
    --config configs/meteosat-10_seviri_cloud_liquid-water_stg_test.driver --dry-run
# 5. submit
scripts/submit_oraclut_lut.sh --config configs/meteosat-10_seviri_cloud_liquid-water_stg_test.driver
# 6. watch the queue
squeue -u "$USER"
# 7. inspect the logs (job ID from the sbatch message)
less validation/slurm/oraclut_lut_<jobid>.out
cat  validation/slurm/oraclut_lut_<jobid>.err
# 8. the product is where the log's "output:" line says, by default
ls create_orac_lut/luts/meteosat-10_seviri_cloud/
```

Cloud (discrete particle layer) example:

```shell
scripts/submit_oraclut_lut.sh \
    --config configs/meteosat-10_seviri_cloud_liquid-water_stg_test.driver
```

Aerosol (distributed profile with surface-pressure dimension) example:

```shell
scripts/submit_oraclut_lut.sh \
    --config configs/meteosat-10_seviri_aerosol_a79_test.driver
```

Both examples have been run successfully through SLURM (jobs 435685 and
435686, node `atmnode014`), producing
`create_orac_lut/luts/meteosat-10_seviri_cloud/meteosat-10_seviri_m_liquid-water_a01_pstg_v21.nc`
and
`create_orac_lut/luts/meteosat-10_seviri_aerosol/meteosat-10_seviri_m_aerosol_a12_pa79_v21.nc`
(surface-pressure dimension 3).

What `scripts/submit_oraclut_lut.sh` does:

1. refuses to run unless the host is `atmlxint7` (clear message, no ssh);
2. changes to the repository root;
3. creates `validation/slurm/` **before** calling `sbatch` — SLURM opens the
   `#SBATCH --output/--error` files before the batch script runs, so a
   `mkdir` inside the script is too late (this is why job 435683 failed);
4. checks that any `--config` file exists;
5. removes the login shell's per-user `TMPDIR` (`/tmp/user/<uid>`) from the
   submitted environment (see below);
6. runs `sbatch scripts/run_oraclut_lut.slurm <generator arguments>` and
   prints sbatch's "Submitted batch job N" line.

Direct generator options are still accepted for testing and unusual runs,
e.g. `scripts/submit_oraclut_lut.sh --platform meteosat-10 --instrument seviri
--microphysics liquid-water_stg.mm --test --channels 1`. No `PYTHONPATH`,
Python executable, Mie/DISORT paths or output naming need to be given; they
live in `scripts/run_oraclut_lut.slurm` and the configuration.

### SLURM resources and overrides

`scripts/run_oraclut_lut.slurm` requests conservative defaults on the
cluster's default/shared partition (no partition or account is named):

```text
--time=04:00:00   --cpus-per-task=1   --mem=8G
```

Override them on the wrapper command line; these options are forwarded to
`sbatch` and take precedence over the `#SBATCH` defaults:

```shell
scripts/submit_oraclut_lut.sh --time=12:00:00 --mem=16G \
    --config configs/meteosat-10_seviri_cloud_liquid-water_stg.driver
```

Recognised sbatch overrides: `--time`, `--mem`, `--cpus-per-task`,
`--job-name`, `--partition`, `--account`, `--qos` (in `--opt=value` or
`--opt value` form). Anything else goes to the generator; a literal `--`
ends the sbatch options. The calculation is single-threaded; do not raise
`--cpus-per-task` expecting a speed-up.

### Temporary files

The batch script creates `/network/scratch/grainger/oraclut-tmp/job_<jobid>/`
after the job starts, verifies it is writable, exports `TMPDIR`/`TMP`/`TEMP`
to it, and removes it on exit (normal or failure) via an `EXIT` trap. Nothing
else under scratch is touched.

The earlier log message

```text
slurmstepd-atmnode014: error: Unable to create TMPDIR [/tmp/user/27004]: Permission denied
slurmstepd-atmnode014: error: Setting TMPDIR to /tmp
```

is emitted by `slurmstepd` **before the batch script starts**, because the
login shell's `TMPDIR=/tmp/user/<uid>` is forwarded by `sbatch` and cannot be
created on the compute node. It therefore cannot be suppressed from inside
the script; the wrapper prevents it by unsetting `TMPDIR` in the submission
environment. If you call `sbatch` directly instead of using the wrapper, run
`unset TMPDIR` first or expect the harmless warning.

## Runtime considerations

Follow the AOPP execution policy in `AGENTS.md`. The compact `*_test.lut`
cases take 15–30 s and may run interactively on `atmlxint1/2/3/5/6` (never
`atmlxint4`). Production grids are roughly four orders of magnitude larger
(`--dry-run` prints the RT-state estimate: ≈150 000–180 000 states for the
Meteosat-10 SEVIRI production grids, versus 18–36 for the test grids). Expect
hours: submit them through SLURM as above, and raise `--time` accordingly.
Only if a run must be interactive and may exceed about 30 minutes, use
`nice -n 19`. Always dry-run first and scale from a measured test.

## Validation status

Channel-1 (visible) agreement alone is not sufficient evidence that a
formulation reproduces the legacy radiative transfer, because a solar channel
never exercises the thermal emission route. Validation therefore uses a
visible + infrared matrix — SEVIRI channel 1 (0.64 µm, solar) and channel 9
(10.8 µm, thermal-only) — for both formulations, reported separately, with
legacy references produced from the coherent extracted working tree. Full
detail, residual tables and reproduction commands: `docs/validation_matrix.md`.

| Formulation | Visible (ch 1) | Infrared (ch 9) | Visible + IR (ch 1, 9) | Status |
| --- | --- | --- | --- | --- |
| cloud | GREEN, 1.3e-4 | GREEN, 6.6e-5 (`E_md` 1.5e-6) | GREEN | **validated on this matrix** |
| aerosol | GREEN, 9.5e-5 (after correcting the profile-interpolation boundary convention found by the per-layer diagnostic) | GREEN, `E_md` 6.0e-8 (after correcting a ×100 emissivity scaling found by this matrix) | GREEN | **validated on this matrix** |

Both aerosol findings — a factor-100 scaling of `E_md` in the aerosol emission
block, and `np.interp` clamping the vertical profile where IDL `INTERPOL`
extrapolates linearly — were found by this matrix, corrected with one-line and
one-function changes, and are documented with before/after records in
`docs/validation_matrix.md`. The aerosol NetCDF layout now matches legacy completely — dimension
declaration order and the full `surface_pressure` coordinate metadata
including `spacing`, preserved from the LUT definition — with structural PASS
and no attribute differences on all aerosol cases, and scientific arrays
unchanged; see `docs/validation_matrix.md`, Finding 3.

The legacy reference for the singleton-channel cases carries fill values
(9.97e36) in the optics variables — a known writer artefact — where Python
writes finite values; this is reported, not equated.

## Current limitations

- Only single-view instruments; dual-view (`view > 0`) fails explicitly.
- Optics: single-component modified-gamma (STG) and shared-mode log-normal
  Mie component mixtures (e.g. `a79`). T-matrix (`tmatrix_path`), Baum/Baran
  ice tables, `force_n`/`force_k`, `n_theta`, `srfdat` and `opt_prop_luts`
  are not ported and are rejected explicitly.
- Gas absorption is implemented for the aerosol formulation only.
- The two shipped production configurations have been resolved and costed by
  dry run, but not yet generated; only the compact `_test` cases have been run
  end to end and validated.
- `srf_quad=2` (band integration) is accepted but has not been validated
  against a legacy band-mode reference.
- LUT level 2 / V2 NetCDF layout only.
