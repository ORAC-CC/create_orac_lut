# Repository development workflow review

> Moved on 2026-10-02 from `validation/` to `validation/repository_workflow/`
> (validation-directory tidy, see
> `validation/REPORT_lut_numerics_development.md`). File references below
> describe the layout at the time of the review: `validation/README_nakajima_king.md`
> is now `validation/nakajima_king/README.md`, and
> `validation/REPORT_legendre_expansion_length.md` is now
> `validation/legendre_expansion/REPORT_legendre_expansion_length.md`.

Date: 2026-10-02

This is a workflow design review only. No scientific code, run file, test,
input, production safeguard, Git branch, worktree, or external directory was
changed or created.

## 1. Current repository situation

The repository at `/home/g/grainger/project-oraclut` is already the natural
single scientific-development laboratory. Its main areas have distinct roles:

- `src/oraclut/` and `create_orac_luts.py`: Python production implementation;
- `create_orac_lut/` and `mie/`: authoritative legacy/reference sources and
  preserved scientific kernels;
- `create_orac_lut/input_files/` and `configs/`: scientific configuration and
  input data;
- `runs/`: executable LUT configurations;
- `scripts/`: operational and submission helpers;
- `tests/`: automated regression tests;
- `validation/`: local diagnostics, reports, plots, numerical tables, SLURM
  logs, and generated comparison products;
- `build/`, root `luts/`, and `create_orac_lut/luts/`: generated or compiled
  products rather than source.

At inspection, `HEAD` was `9d663e9` on `main`, exactly tracking
`origin/main`. A legacy `master` branch exists at `c182ef6`, tracking
`origin/master`; it is not a parallel active development line. Tags `v0.0.0`
and `v0.0.1` both exist. The GitHub remote is
`ORAC-CC/create_orac_lut.git`.

The working tree was already deliberately non-clean. Pre-existing changes were:

- modified `AGENTS.md` (the added immutable-IDL-reference rule);
- modified `tests/test_baum.py`;
- untracked `documents/Grainger 1990.pdf`;
- 40 untracked MODIS/SLSTR V23 run files;
- untracked `scripts/submit_v23_slstr_cloud_luts.sh`;
- untracked `tests/test_nakajima_king.py`.

The tracked run set contains `runs/template.run`, three compact cloud test
configurations, three compact aerosol test configurations, and one EarthCARE
aerosol case. The 40 local V23 files follow a consistent
`<platform>_<instrument>_<forward-model>_<microphysics>_v23.run` convention.
Production run files explicitly identify version, channels, grid,
microphysics, numerical settings, and output directory. Submission scripts add
preflight checks, skip complete products, support dry runs, and rely on the
generator's no-overwrite protection.

Generated data are mostly separated from source. `.gitignore` excludes Python
and compiler products, root and legacy LUT output directories, runtime logs,
large Baum staging data, and named validation scratch locations. The current
`validation/` tree contains reports and reusable scripts alongside generated
plots/tables, caches, and about 41 MB of SLURM evidence and other artifacts.
There is no dedicated `development/`, `experiments/`, or repository-local
general scratch directory at present, although `validation/tmp/` and
`validation/scratch/` are reserved in `.gitignore`.

One important inconsistency exists: the final `/validation/` rule ignores the
entire validation tree, including reproducible scripts and reports. Likewise,
`/docs/`, `ARCHITECTURE.md`, `MIGRATION_PLAN.md`, and `VALIDATION.md` are
ignored. `CLAUDE.md` refers readers to those three root design documents, but
they are not currently present. This makes local validation convenient but
prevents GitHub from serving as a complete reproducibility archive unless
files are force-added.

## 2. Recommended development model

Use one authoritative working directory and a staged promotion model:

1. **Investigate in an isolated laboratory area inside this repository.** Put
   disposable scripts, prototypes, reduced configurations, notebooks, timing
   logs, and outputs under `development/<topic>/`. Treat the whole area as
   ignored and removable. It must not be on production import paths.
2. **Evaluate scientifically before changing production.** Reproduce the
   baseline, state the hypothesis and acceptance measures, run reduced tests,
   compare intermediate quantities, and retain enough local evidence to make a
   decision. A failed idea may simply be deleted; Git is not its laboratory
   notebook.
3. **Promote accepted evidence to validation.** Move or rewrite the smallest
   self-contained diagnostic into `validation/<topic>/`, with a short report,
   exact invocation, assumptions, tolerances, and compact results needed to
   reproduce the conclusion. Do not promote caches or bulk output.
4. **Implement an accepted decision in production paths.** Only after the
   scientific choice is explicit should the coherent change be made in
   `src/`, supported legacy-facing Python, tests, inputs, or operational
   configuration. Preserve the legacy source and validate the production path.
5. **Commit the definitive state.** A commit records the accepted code,
   configuration, focused tests, reproducible validation, and documentation as
   one reviewable scientific/software decision. Push that authoritative state
   to GitHub.

This model allows a dirty scientific workbench without confusing it with an
authoritative state. Before promotion or production use, changes must be
classified explicitly rather than assuming that every uncommitted file is
either wrong or ready.

### Recommended locations

| Material | Recommended location | Git treatment |
|---|---|---|
| Disposable experimental code and reduced run configurations | `development/<topic>/` | Ignore; never import from production entry points |
| Temporary numerical output, profiles, compiler products, and timing logs | `development/<topic>/output/` | Ignore and delete when no longer useful |
| Scratch generated by an accepted validation script | `validation/scratch/<topic>/` or `validation/tmp/<topic>/` | Ignore |
| Reproducible validation script, report, compact table/figure, or deliberately selected reference | `validation/<topic>/` | Track after scientific acceptance |
| Generated validation collections and caches | `validation/generated/<topic>/`, cache directories | Ignore |
| Production Python | `src/oraclut/` and established entry points | Track only accepted coherent changes |
| Operational run files | `runs/` | Reserve for reviewed, run-ready configurations |
| Large LUT products | established output directories, separate from source | Do not track unless explicitly selected as a compact regression reference |

If this model is adopted, make a separate, reviewed hygiene change that adds
`/development/` to `.gitignore` and replaces the blanket `/validation/` ignore
with specific generated areas such as `/validation/generated/`,
`/validation/scratch/`, `/validation/tmp/`, `/validation/slurm/`, caches, and
known large binary patterns. Do not perform that change opportunistically
during a scientific investigation.

Experimental code must not be placed in `src/oraclut/`, production scripts, or
`runs/` merely to make it easy to execute. Invoke it by explicit path, keep any
experimental dependencies local to its topic directory, and do not add
`development/` to `PYTHONPATH`. When an experiment needs a production routine,
import the authoritative routine from `src/`; do not copy it and allow the copy
to drift.

## 3. Git and GitHub policy

- Git tracks authoritative scientific states, not every command or exploratory
  edit.
- GitHub is the reproducibility archive for validated source, accepted
  configuration, compact reference evidence, and the reasoning behind
  scientific decisions.
- Experimental intermediate states do not require commits. Failed experiments,
  temporary scripts, caches, scheduler logs, and generated debris should not
  enter history.
- A commit represents one validated scientific/software decision. It should
  state the baseline, scientific reason, affected behavior, validation target,
  tests, tolerances, and any intentional legacy deviation.
- Large generated LUTs do not enter Git unless there is an explicit,
  documented regression need and their size is deliberately acceptable.
- `main` is the continuing authoritative line. Do not create sibling repository
  copies, duplicate clones as development environments, permanent experimental
  branches, or multiple maintained ORAC variants. Do not create worktrees
  unless the user explicitly authorises an exceptional need.
- Tags may identify validated releases or campaign baselines; they should not
  label routine experiments.

The normal sequence is therefore: local investigation, scientific decision,
focused production change, validation, commit, and push. Cleanup of the ignored
laboratory can occur before or after the commit because it is not part of the
authoritative state.

## 4. Production safety rules

Before a production calculation:

1. Record `git rev-parse HEAD` and inspect `git status --short`. A completely
   clean tree is ideal, but where unrelated laboratory work is intentionally
   present, require production-critical paths to match the chosen authoritative
   commit and enumerate any accepted run-file changes explicitly.
2. Use only reviewed run files from `runs/`; keep experimental run files under
   `development/<topic>/`. Confirm platform, instrument, microphysics, grid,
   channels, version, numerical settings, output path, and expected product
   name.
3. Use the existing read/preflight or `--dry-run` path before submission. Check
   currently running jobs and submit only the intended products; never use a
   broad product-set script for a selective restart.
4. Preserve existing submission safeguards, SLURM host policy, resource
   declarations, completion checks, and no-overwrite behavior. Do not point
   operational scripts at experimental modules.
5. Keep LUT products and scheduler logs separate from source. Record a manifest
   containing the commit, run file, input versions/checksums where practical,
   command, job identifier, and final product path.
6. Promote a generated file into Git only as an explicitly selected compact
   regression reference, with provenance and a documented reason.

This preserves released-version reproducibility while allowing unrelated
experiments to remain in the same physical project directory.

## 5. Documentation review

Reviewed Markdown comprised `AGENTS.md`, `CLAUDE.md`,
`validation/README_nakajima_king.md`,
`validation/REPORT_legendre_expansion_length.md`, and the two Markdown files in
`validation/size_distribution_limits/`. `.gitignore`, run templates,
representative V23 runs, submission scripts, and the validation directory were
also inspected for their implied workflow.

No Markdown recommends sibling repository copies, parallel maintained
repositories, long-lived feature branches, automatic merging, or committing
every experiment. `AGENTS.md` already establishes the repository boundary,
incremental scientific validation, source/output separation, and compact
reviewable changes. The validation documents use a sensible convention of a
focused script plus a short scientific report.

Three documentation issues should be resolved in a future dedicated policy
change:

1. narrow the blanket `/validation/` ignore so accepted validation can be
   archived normally;
2. either add and track the design documents referenced by `CLAUDE.md`, or
   remove those stale references;
3. clarify the sentence in `REPORT_legendre_expansion_length.md` saying that
   production Python routines were used "outside the repository", because it
   can be read as endorsing an external code copy even if it intended to mean
   that the calculations were run without changing repository source.

No existing documentation was edited in this review. The only Markdown file
created is this requested report.

## 6. Proposed future AGENTS.md wording

The following can be added after a future policy decision; it was not added in
this review:

> ## Repository-contained scientific development
>
> `/home/g/grainger/project-oraclut` is the single authoritative ORAC LUT
> working directory. Never create a sibling repository, duplicate project copy,
> or alternate authoritative tree. Never create a Git worktree or development
> branch unless the user explicitly authorises it for a specific task.
>
> Keep disposable investigations under the repository's ignored
> `development/<topic>/` area, with temporary output beneath the same topic or
> an approved ignored validation scratch directory. Keep this area off
> production import paths. Do not place experimental run files in `runs/` and
> do not modify production code during an investigation unless the user
> explicitly requests implementation.
>
> When an investigation produces an accepted scientific result, preserve the
> smallest reproducible script, report, compact evidence, and provenance under
> `validation/<topic>/`; then implement and test the coherent production change.
> Git commits represent validated scientific/software decisions. Do not commit
> failed prototypes, temporary scripts, caches, scheduler logs, or generated
> debris, and do not commit large LUT products without an explicit documented
> reason.
>
> Before operational LUT generation, inspect Git status, identify the exact
> authoritative commit and reviewed run file, ensure production-critical files
> contain no experimental changes, use existing preflight/dry-run safeguards,
> and keep generated products separate from source.

## 7. Alternatives considered and rejected

### A. Experiments directly throughout the working tree

This has the lowest startup cost but the highest ambiguity. Untracked scripts
can be mistaken for supported tools, experimental imports can leak into
production, and cleanup becomes risky in an already active tree. Direct edits
are appropriate only after a result is accepted and is being implemented as a
coherent production change. They are not the default laboratory model.

### B. Dedicated ignored area inside the repository

This is the recommended default, specifically `development/<topic>/`. It keeps
the work physically inside the authorised project, makes disposability obvious,
avoids Git-history noise, and gives every experiment a bounded cleanup target.
Its risks are unmanaged accumulation and poor reproducibility; topic naming,
periodic deletion, and promotion of accepted evidence to `validation/` address
those risks.

### C. `/tmp` as the only experiment area

Rejected as the normal workflow. It lies outside the repository boundary,
loses provenance easily, may be cleaned without warning, and encourages scripts
to depend on unrecorded paths. System-created runtime temporaries used by
established code are a separate implementation detail, but agents should not
build the scientific workflow around `/tmp` without explicit permission.

### D. Combined model

Recommended: ignored `development/<topic>/` for exploration, tracked
`validation/<topic>/` for accepted reproducibility evidence, established
production paths only for accepted implementation, and Git/GitHub only for the
validated authoritative state. This provides isolation without another
repository, permanent branch, worktree, or version line.
