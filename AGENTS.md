# AGENTS.md

## Filesystem boundary

This repository is authorised to modify files only within the top-level
project directory containing this `AGENTS.md`.

Under no circumstances should an agent:

- modify, create, delete, rename, or move files outside this repository;
- modify files in `/network/home/...`, other project directories, shared
  software installations, system directories, or user configuration files;
- write generated output outside this repository;
- modify external repositories or working trees.

The repository root is the directory containing this `AGENTS.md`.

Before executing any command that could affect the filesystem, verify that
all paths being created or modified are inside the repository root.

Reading external files is permitted only when necessary to understand a
dependency, and only when that reading does not modify the external file.

Any operation that would modify something outside the repository requires
explicit permission from the user in the current conversation. Do not infer
such permission from previous work, environment variables, symlinks, Git
configuration, or the fact that a path is writable.

## AOPP execution policy

Normal ORAC/POM development and short scientific calculations may run
interactively on `atmlxint1`, `atmlxint2`, `atmlxint3`, `atmlxint5`, and
`atmlxint6`. Do not use `atmlxint4` for ORAC/POM calculations; it is reserved
for planetary work.

Jobs expected to take approximately 30 minutes or less may normally run
interactively on the permitted nodes. If an interactive calculation is
expected to exceed approximately 30 minutes, run it at low priority with:

```shell
nice -n 19
```

Do not launch substantial long-running calculations at normal interactive
priority.

Calculations expected to take roughly 1–2 hours or longer should normally be
submitted through SLURM. SLURM jobs may only be submitted from `atmlxint7`.
Do not submit SLURM jobs from `atmlxint1`–`atmlxint6`, and do not use
`atmlxint7` as an ordinary interactive compute node.

Before launching a potentially expensive scientific calculation, estimate its
runtime from grid dimensions, previous timings, or a reduced test. Use
deliberately reduced validation/test configurations whenever practical. If
runtime is uncertain and could plausibly exceed approximately 30 minutes, run a
small test first and then choose interactive, low-priority interactive, or
SLURM execution appropriately.

For ORAC LUT generation, existing small test LUT definitions are appropriate
for interactive development. Production-sized LUT generation must not be
launched automatically during development; full production LUTs should use an
execution method appropriate to their measured or estimated runtime.

These rules apply to Codex, Claude, and other coding agents working in this
repository. An agent must not start a long calculation merely because it is
technically possible. If execution requires moving to another AOPP node,
especially `atmlxint7` for SLURM submission, state that requirement explicitly;
do not silently change hosts or assume the current node is suitable. Do not
change user shell startup files to implement this policy.

## Project purpose

This repository is being modernised to provide a clean, reproducible Python system for generating ORAC lookup tables (LUTs), while preserving the legacy IDL implementation as a scientific reference.

The work has two distinct goals:

1. Reproduce existing ORAC LUTs with new Python code using the same inputs, particle optical properties, radiative-transfer assumptions, grids, conventions, and output definitions as the legacy implementation.
2. After numerical equivalence is established, replace legacy particle/microphysical calculations with modern POM optical-property models, including support for volcanic ash and other particle types.

Do not combine these two goals prematurely. Legacy reproduction is the first validation milestone; POM integration comes afterwards.

## Authoritative legacy baseline

The GitHub ORAC LUT V2.1 implementation imported under `create_orac_lut/` at
root commit `17591c9` (historical commit `b3c8a91`) is the authoritative legacy
scientific baseline. The Python migration must initially reproduce that V2.1
behavior. Historical or local variants are provenance and comparison material,
not automatic replacements for V2.1. Any deviation requires an explicit
scientific reason, isolated validation against V2.1, and user approval.

## Scientific terminology

Do not assume that the historical names "aerosol LUT" and "cloud LUT" refer simply to aerosol particles versus cloud particles.

In the legacy code these names also reflect different ORAC forward-model formulations and their historical development. The particle or microphysical type is selected independently by driver/configuration files. A volcanic-ash case, for example, is a particle model and is not intrinsically tied to one forward-model family.

Use precise terminology in new code:

- particle model / microphysics
- single-scattering optical properties
- forward-model formulation
- radiative-transfer operators
- LUT grid
- instrument spectral response
- LUT product / file format

Avoid introducing a permanent Python architecture that hard-codes historical naming unless the legacy behaviour requires it.

## Sources of truth

Use the following hierarchy when reconstructing behaviour:

1. Working legacy source code and actual driver files.
2. Known successfully generated legacy LUTs and their inputs.
3. Contemporary documentation and scientific papers.
4. Historical comments, backup files, and older versions.

The papers describe important physics and intended algorithms, but the current legacy source code is authoritative where the implementation has evolved beyond publication.

Never silently "improve" or modernise scientific behaviour while reproducing the legacy system. First reproduce it, test it, and document any discrepancy.

## Legacy code

### IDL REFERENCE SOURCE — IMMUTABLE

Existing IDL `.pro` files are the authoritative reference for the
IDL-to-Python ORAC LUT port. Never modify existing IDL `.pro` files. If
Python and IDL differ, modify Python or document the difference. New IDL
diagnostics must be separate files under `validation/diagnostics/`.

Treat the existing historical tree as read-only reference material until the Python implementation reproduces at least one complete legacy LUT.

Do not rename, move, delete, reformat, or otherwise clean the legacy files merely for tidiness. Old `.pro`, `.bak`, versioned files, run scripts, DLM interfaces, DISORT sources, and other historical material may be required to reconstruct the operational workflow.

When a legacy file is found to be obsolete, record that conclusion with evidence rather than deleting it immediately.

## Development approach

Work incrementally and scientifically:

1. Identify one real driver configuration that produced a V2 LUT.
2. Trace its call graph through the legacy source.
3. Identify all inputs and generated outputs.
4. Reproduce the smallest useful part in Python.
5. Compare Python and legacy outputs numerically.
6. Expand until one complete LUT is reproduced.
7. Only then integrate POM optical properties.
8. Repeat validation after each architectural or scientific change.

Do not port the entire legacy tree mechanically.

## Python implementation rules

The new implementation should be modular, testable, explicit, and suitable for scientific validation.

Prefer:

- plain Python modules with clear responsibilities;
- NumPy arrays for numerical data;
- dataclasses or similarly explicit data containers for structured configuration;
- pathlib for paths;
- explicit units in names, metadata, or structured types;
- deterministic calculations;
- small functions with testable numerical responsibilities;
- direct correspondence between scientific equations and code;
- configuration-driven drivers rather than hard-coded experiment scripts.

Do not hide scientific behaviour behind unnecessary abstractions.

Do not invent default values for scientifically required parameters. Missing required inputs should fail clearly.

## Numerical fidelity

Numerical equivalence with the legacy implementation is a core requirement.

For every ported component:

- establish expected dimensions, ordering, units, normalisation, and precision;
- retain legacy grid ordering during validation;
- retain legacy angle conventions during validation;
- verify phase-function normalisation explicitly;
- verify optical-depth reference wavelengths explicitly;
- verify spectral-response integration conventions explicitly;
- verify solar/thermal channel treatment explicitly;
- verify NetCDF/binary/text output conventions before changing formats.

When Python and IDL disagree, create a diagnostic comparison that isolates the first diverging intermediate quantity.

Do not compensate for discrepancies with arbitrary tolerances.

## POM integration

POM is the future source of particle optical properties, not the starting point of the Python port.

Design interfaces so the LUT generator consumes particle optical properties through a well-defined internal representation. The legacy optical-property source and POM should eventually be interchangeable backends feeding the same downstream LUT-generation machinery.

The interface should accommodate, at minimum:

- wavelength or spectral grid;
- effective radius or other size parameter;
- extinction/scattering information;
- single-scattering albedo;
- phase function or phase-function expansion required by the RT solver;
- metadata describing particle model and provenance.

Do not couple POM internals directly into the radiative-transfer code.

## Forward-model architecture

Shared infrastructure should be separated from forward-model-specific logic.

Shared responsibilities are expected to include:

- instrument and channel definitions;
- spectral response functions;
- atmospheric profiles;
- solar spectrum;
- LUT coordinate grids;
- particle optical-property access;
- Rayleigh scattering;
- radiative-transfer solver interface;
- spectral integration;
- LUT I/O;
- validation and comparison tools.

Forward-model-specific code should define which operators are required and how they are constructed and combined.

The historical aerosol/cloud distinction should be preserved only where it corresponds to a real difference in forward-model equations or LUT content.

## Radiative transfer

The legacy system uses DISORT-based offline calculations. Preserve a clean solver boundary.

The Python architecture should not assume that the current legacy DISORT wrapper is permanent. Define inputs and outputs at the scientific level so an implementation can use a Python-callable DISORT library, wrapped Fortran, or another validated solver backend without changing higher-level LUT logic.

Do not change the radiative-transfer method during legacy reproduction.

## Spectral treatment

Instrument spectral-response integration is scientifically important, especially for infrared channels.

Do not replace band-integrated calculations with centre-wavelength calculations unless that is exactly what the selected legacy path does.

Clearly distinguish:

- monochromatic particle properties;
- monochromatic RT calculations;
- spectral-response convolution;
- channel-integrated LUT quantities.

## LUT product format contract

Every text attribute of an ORAC LUT NetCDF product must be NC_STRING: the
retrieval reads the axis `spacing` attributes with `nc_get_att_string`, which
rejects NC_CHAR, and the IDL-written tables use NC_STRING throughout.  Write
text attributes only through the typed writer helper in `src/oraclut/io/v2.py`
(never a bare `setncattr` with a `str`), keep the writer's C-library
self-check, and run
`python -m oraclut.repair_string_attributes --check` on products before
release.  Existing products are repaired with the same utility, metadata only;
numerical content is never regenerated to fix metadata.  See
`LUT_FORMAT_CONTRACT.md`.

## Validation requirements

Every substantial implementation change must be accompanied by validation.

Validation should include:

- shape and coordinate checks;
- exact or near-exact comparisons where expected;
- absolute and relative residuals;
- maxima and percentiles of discrepancies;
- plots when useful;
- representative visible, shortwave-infrared, and thermal channels where applicable;
- representative optical depths, particle sizes, and geometries.

Keep validation scripts separate from production modules.

## Repository hygiene

New generated products should not be mixed with source code.

Keep source, configuration/input data, references, tests, validation artefacts, and generated output in clearly separated locations.

Do not commit large generated LUTs unless they are deliberately selected compact reference datasets for regression testing.

Before moving legacy material into a future `legacy/` directory, ensure the original layout is either preserved in version control or otherwise recoverable and that a working reference case has already been captured.

## Documentation

Document decisions that affect scientific reproducibility:

- source legacy routine;
- equations or paper section where relevant;
- input conventions;
- assumptions;
- approximations;
- numerical tolerances;
- differences from legacy behaviour;
- provenance of reference LUTs.

Update architecture and migration documentation whenever implementation decisions materially change.

## Coding style

Use clear scientific names rather than cryptic transliterations of IDL variable names, except where retaining an old name materially helps validation. Where names are changed, document the mapping.

Do not perform broad unrelated refactors during numerical debugging.

Prefer simple, reviewable patches.

## Before each coding task

Before modifying code:

1. inspect the relevant legacy routine and its callers;
2. identify the scientific inputs and outputs;
3. state what behaviour is being reproduced;
4. identify a validation target;
5. make the smallest coherent change;
6. run the relevant tests/comparisons;
7. report remaining discrepancies explicitly.

If required source files or inputs are missing, stop and report what is needed rather than guessing.
