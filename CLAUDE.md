# CLAUDE.md

This repository is a scientific modernisation of the ORAC lookup-table generator.

Read `AGENTS.md` before making changes. Treat it as the primary repository-wide development policy.

Also read:

- `ARCHITECTURE.md` for the intended Python architecture;
- `MIGRATION_PLAN.md` for the staged reconstruction strategy;
- `VALIDATION.md` for numerical acceptance and comparison rules.

Key constraints:

- The legacy IDL implementation is a scientific reference and must not be casually reorganised or "cleaned".
- First reproduce an existing V2 ORAC LUT in Python using legacy inputs and optical properties.
- Integrate POM only after legacy numerical reproduction is demonstrated.
- "Aerosol LUT" and "cloud LUT" are historical forward-model labels, not simply particle categories.
- Driver/configuration files select the particle/microphysical model independently.
- Preserve units, dimensions, grid ordering, angle conventions, spectral-response treatment, optical-depth reference wavelengths, phase-function normalisation, and solar/thermal behaviour during validation.
- Do not invent missing scientific parameters or silently substitute modern formulations.
- When results disagree, isolate the first differing intermediate quantity and diagnose it quantitatively.
- Keep production code, inputs/configuration, references, validation code, tests, and generated output separate.
- Prefer small, reviewable, scientifically testable changes.

Do not duplicate repository policy in this file. If a rule needs to change, update `AGENTS.md` and the relevant design document so Codex and Claude remain aligned.
