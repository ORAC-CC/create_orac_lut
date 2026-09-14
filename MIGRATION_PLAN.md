# MIGRATION_PLAN.md

## Objective

Replace the historical ORAC LUT-generation workflow with a clean Python implementation without losing scientific reproducibility.

## Phase 0 - Preserve and inventory

Status: the selected Meteosat-10 SEVIRI liquid-water STG cloud implementation
and full V2 output path are complete. Tiny structure is exact and its
radiative operators agree at single-precision scale, but the captured tiny
file has fill values where Python has finite microphysical optics. The full
case has a documented preserved-DISORT near-singular residual. These issues
must be resolved before claiming scientific equivalence. Broader particle and
forward-model coverage remains future work.

- Preserve the historical working tree unchanged.
- Record Git status, but do not assume Git contains every important file.
- Treat actual filesystem contents as potentially authoritative.
- Identify active driver files, top-level LUT routines, I/O routines, optical-property generators, and DISORT interfaces.
- Separate clearly documented facts from inferred workflow.

Deliverable: a documented map of the legacy system.

The initial forensic map confirms the operational pattern
`makerunfile_v2.pro -> <platform>_<instrument>_run -> IDL wrapper -> LUT
generator`, but it also finds later split aerosol/cloud wrappers and an
incomplete V2.1 redesign. The exact wrapper and command-line overrides must be
recorded with each reference case; `v21` in an output filename is not by itself
enough to identify the call path. See `docs/legacy_lut_workflow.md` and
`docs/legacy_lut_call_graph.md`.

## Phase 1 - Identify a real V2 reference case

Find one LUT that was genuinely generated successfully in the later operational workflow.

Determine:

- driver/run file;
- whether it invokes original, V2, wrapper, aerosol-named, or cloud-named routines;
- particle/microphysical type;
- instrument and channels;
- input files;
- forward-model formulation;
- LUT dimensions and ordering;
- output files;
- code revision if recoverable.

The phrase "V2 LUT" must not be assumed to imply that `create_orac_lut_v2.pro` alone generated it. Establish the actual call path.

Deliverable: one frozen legacy reference case.

Current candidates are ranked in `docs/reference_case_candidates.md`. The
recommended first case is the existing Meteosat-10 SEVIRI liquid-water `stg`
V21 product, followed by the pressure-aware Meteosat-10 aerosol `a79` product.

## Phase 2 - Reconstruct the legacy call graph

Trace the selected driver through the exact routines it calls.

For each routine record:

- caller;
- inputs;
- outputs;
- side effects/files;
- units;
- dimensions;
- numerical conventions;
- external dependencies.

Focus first on the selected reference case rather than the full historical codebase.

Deliverable: a call/data-flow document sufficient to implement the case independently.

## Phase 3 - Establish Python package skeleton

Create only the modules required by the selected case.

Suggested initial responsibilities:

- configuration;
- instrument/SRF reader;
- atmosphere reader;
- LUT grid;
- particle-optics reader/interface;
- radiative-transfer interface;
- LUT assembly;
- output writer;
- validation tools.

Avoid creating placeholder abstractions for parts of the legacy system that have not yet been studied.

Deliverable: importable Python package with tests for parsing and data structures.

## Phase 4 - Reproduce intermediate quantities

Port the legacy workflow in scientifically meaningful stages.

Likely comparison points include:

- input grids;
- spectral sampling;
- particle optical properties;
- optical-depth scaling;
- Rayleigh quantities;
- phase-function representation;
- DISORT layer properties;
- monochromatic RT operators;
- SRF-integrated channel operators;
- final LUT arrays.

At every stage compare Python against IDL before proceeding.

Deliverable: quantitative comparison report with no unexplained divergence in completed stages.

The selected complete path is implemented in `src/oraclut/pipeline.py`.
`validation/generated/python_full_comparison.json` records the full comparison
against the read-only current production reference; difficult DISORT
near-singular cases are validated with the temporary Intel 2022 build.

## Phase 5 - Reproduce one complete legacy LUT

Generate the same selected LUT with Python.

Compare:

- coordinates;
- dimensions;
- metadata;
- each stored operator/quantity;
- numerical residuals.

Define acceptance criteria from the observed numerical precision and algorithmic equivalence rather than using arbitrary broad tolerances.

Deliverable: validated Python reproduction of one complete LUT. The captured
tiny V2 reference and the current full Meteosat-10 production reference are
both generated and compared by the Python pipeline; structural identity,
residuals, compiler choice, and the remaining full-case blocker are documented
in `docs/python_legacy_reproduction.md`.

## Phase 6 - Broaden legacy coverage

Add representative cases spanning:

- both historical forward-model families if still scientifically relevant;
- visible channels;
- shortwave infrared;
- thermal infrared;
- different instruments;
- different particle models.

Do not duplicate code merely to reproduce historical file naming.

Deliverable: regression suite covering the active legacy use cases.

## Phase 7 - Introduce POM

Implement a POM-backed optical-property provider conforming to the common optical-property interface.

Initially keep:

- the same LUT grids;
- the same RT solution;
- the same instrument treatment;
- the same output definitions.

Change only the particle optical properties.

Compare POM-backed output first for a particle case where legacy and POM physics should agree closely, then progress to nonspherical volcanic ash and other new models.

Deliverable: scientifically attributable difference between legacy and POM LUTs.

## Phase 8 - Repository restructuring

Only after working Python reference cases exist:

- move preserved historical material under `legacy/`;
- move papers and documentation under `references/`;
- place reusable input datasets under `inputs/`;
- keep generated products under `output/`;
- maintain compact validation reference datasets under `validation/reference_luts/`;
- update all paths and documentation.

Do not delete history simply because it is no longer part of the production path.

## Phase 9 - Modern scientific evaluation

Once the implementation is stable, assess whether the work supports a new publication.

Potential scientific/technical contributions include:

- reproducible Python replacement of the historical LUT generator;
- documented ORAC forward-model operator architecture;
- modern spectral treatment across solar and thermal channels;
- POM integration;
- modern nonspherical volcanic-ash optics;
- quantitative legacy/new-model comparisons;
- impact on ORAC retrieval behaviour.

Publication is a later outcome, not a constraint on the initial port.
