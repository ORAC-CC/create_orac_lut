# ARCHITECTURE.md

## Purpose

The target system is a Python implementation of ORAC LUT generation that can first reproduce legacy LUTs and later use POM as the particle-optics engine.

The architecture must separate particle microphysics from the ORAC forward-model formulation.

## Conceptual pipeline

```text
driver / configuration
        |
        v
instrument + SRF + atmosphere + LUT grid
        |
        v
particle optical properties
        |
        v
radiative-transfer operator generation
        |
        v
spectral integration / channel treatment
        |
        v
ORAC forward-model-specific LUT assembly
        |
        v
LUT output + metadata
```

During legacy reproduction, particle optical properties come from the same legacy source used by the IDL implementation.

The selected legacy cloud path is now executable from Python. The current
implementation in `src/oraclut/pipeline.py` reads the existing ORAC text
inputs, calls the preserved Mie and production DISORT Fortran kernels through
repository-local C adapters, ports the cloud `setup_disort`/`call_disort`
orchestration, and assembles the observed V2 revision-v21 NetCDF layout. IDL
is retained as the executable reference, not as a runtime dependency of this
Python generator. Details and commands are in
`docs/python_legacy_reproduction.md`.

The user-facing entry point is `python -m oraclut.generate`. It resolves the
original `.inst`, `.lut`, and `.mm` files, then dispatches to the two explicit
legacy orchestration paths. The cloud path keeps the discrete-layer operator
arrays; the aerosol path adds model-profile layer allocation, gas/Rayleigh
composition, pressure-aware operator dimensions, and the 270 K legacy thermal
branch. Both paths share channel/SRF handling, legacy-equivalent Mie optics,
the production DISORT adapter, and the V2 writer. The aerosol implementation
currently supports the shared-mode log-normal Mie component family used by the
selected operational `aerosol_a79.mm` case; unsupported multi-mode or
non-Mie components fail explicitly rather than silently changing physics.

The public names `cloud` and `aerosol` are historical ORAC forward-model names,
not particle categories. In the current source, `cloud` corresponds
approximately to particles confined to a discrete atmospheric layer, while
`aerosol` distributes the particle optical depth through a vertical profile
and adds the pressure-aware output dimension. A particle model, such as
volcanic ash or a liquid-water model, is selected independently of either
formulation. The Python CLI retains these names during reproduction so its
configuration maps directly to `create_orac_cloud_lut` and
`create_orac_aerosol_lut`. Final terminology should be reconsidered only after
both formulations have been compared side by side; no rename is performed in
the reproduction phase.

The observed legacy entry path is not a direct call from a platform
configuration into the scientific routine. `makerunfile_v2.pro` writes a
platform/instrument shell file, which sources host-specific IDL paths, selects
a driver through an environment variable, and invokes an IDL wrapper. The
wrapper parses the driver and applies explicit command-line overrides before
calling the generator. Later current files use separate aerosol and cloud
wrappers, while the older V2-era path uses `create_orac_lut_wrapper_v2`; these
must remain separate during reproduction. The forensic details are in
`docs/legacy_lut_workflow.md`.

After validation, POM becomes an alternative particle-optics backend:

```text
legacy particle optics ----                            >---- common optical-property interface ----> LUT generation
POM -----------------------/
```

## Proposed repository layout

The final layout should be approached incrementally rather than imposed before the legacy workflow is understood.

```text
project-oraclut/
|
|-- AGENTS.md
|-- CLAUDE.md
|-- ARCHITECTURE.md
|-- MIGRATION_PLAN.md
|-- VALIDATION.md
|
|-- src/
|   `-- oraclut/
|       |-- config/
|       |-- instruments/
|       |-- atmosphere/
|       |-- microphysics/
|       |-- optics/
|       |-- radiative_transfer/
|       |-- forward_models/
|       |-- lut/
|       `-- io/
|
|-- inputs/
|   |-- instruments/
|   |-- atmospheres/
|   |-- microphysics/
|   |-- spectral_response/
|   |-- lut_grids/
|   `-- solar/
|
|-- references/
|   |-- papers/
|   `-- documentation/
|
|-- validation/
|   |-- legacy_cases/
|   |-- reference_luts/
|   |-- comparisons/
|   `-- plots/
|
|-- tests/
|-- output/
`-- legacy/          # only after safe migration of the historical tree
```

The existing historical tree should initially remain where it is. Do not move it to `legacy/` until a working reference case has been captured and the move is demonstrably safe.

## Main software boundaries

### Configuration

Represents what legacy run/driver files currently specify.

Configuration should identify, explicitly:

- instrument;
- channels;
- spectral response;
- atmosphere;
- solar spectrum where required;
- LUT coordinate ranges/grids;
- particle/microphysical model;
- forward-model formulation;
- optical-property source/backend;
- output format/location;
- solver settings.

Configuration parsing must remain separate from scientific calculation.

### Instruments

Responsible for channel metadata and spectral response functions.

No microphysical assumptions belong here.

### Atmosphere

Responsible for atmospheric profiles and molecular/Rayleigh quantities needed by the selected forward model.

### Microphysics and optics

`microphysics` describes particle populations and size parameterisation.

`optics` supplies the resulting ensemble optical properties required downstream.

This separation is important because the long-term POM interface belongs at the optical-property boundary rather than inside the LUT or DISORT logic.

### Radiative transfer

Defines a solver-independent scientific API for generating required reflectance, transmission, emissivity, or other operators.

The initial backend must reproduce the legacy DISORT configuration.

### Forward models

Contains only genuine differences in ORAC forward-model formulation.

Do not equate module names automatically with particle type.

A provisional structure may contain historical implementations such as:

```text
forward_models/
   legacy_aerosol.py
   cloud.py
```

but these names should be revisited after the legacy call graph has been reconstructed.

### LUT assembly

Responsible for:

- coordinate definitions;
- dimension ordering;
- operator packaging;
- metadata;
- final output representation.

It should not calculate particle scattering internally.

### I/O

Legacy input and output formats should be isolated behind explicit readers/writers.

This allows scientific reproduction to be validated before deciding whether the long-term storage format should change.

## Design principles

### 1. Reproduce before redesigning physics

Python may improve software structure immediately, but scientific behaviour should remain unchanged until the legacy result is reproduced.

### 2. One internal representation per scientific concept

Avoid allowing each historical code path to invent its own representation of wavelength, optical depth, effective radius, phase function, geometry, or channel data.

### 3. Explicit dimensions and units

Array meaning should be obvious from interfaces and metadata.

### 4. Solver independence

Higher-level LUT code should depend on RT operators, not on the mechanics of an IDL DLM or a particular Fortran wrapper.

### 5. POM independence

Higher-level LUT code should consume optical properties, not POM-specific objects.

### 6. Configuration rather than script proliferation

Historical run files encode configuration by executable code. The new design should move stable scientific choices into declarative configuration where practical, while keeping calculations in Python modules.

### 7. Reproducible provenance

Every generated LUT should be traceable to:

- code version;
- configuration;
- particle model;
- optical-property backend;
- RT backend and settings;
- instrument/SRF data;
- atmospheric data;
- generation date;
- reference wavelength conventions.

## Papers currently supplied as architectural/scientific references

The following uploaded papers should eventually be placed under `references/papers/` using the user's authoritative filenames:

- `2009Thomas2.pdf`
- `2018Sus1.pdf`
- `2018McGarragh1.pdf`

These are historical scientific references, not necessarily exact specifications of the present legacy code. The 2018 McGarragh paper is particularly relevant to the cloud/CC4CL forward model and its offline LUT operators.
