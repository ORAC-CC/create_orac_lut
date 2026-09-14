# Current legacy STG particle optics

## Boundary and call chain

For the selected liquid-water STG model, the current source path is:

```text
create_orac_cloud_lut
  -> generate_scattering_properties
       -> create_bwgp
            -> mie_size_dist_new(..., /dlm, xres=0.4, ...)
                 -> Mie_dlm_single (IDL DLM)
                      -> mieint_ (Fortran)
       -> legpexp
  -> current RT layer construction
```

`create_bwgp.pro` is the common wrapper for Mie and Dubovik T-matrix
calculations. STG selects Mie through its `*scattering code` value. T-matrix
branches are not used for this model.

## Exact `create_bwgp` behaviour

Inputs are distribution name, mode/effective-radius parameters `Rm` and `S`,
complex RI array, wavelength array `wl`, and `Dqv`, the cosine scattering-angle
grid. `wl` and RI must have equal element counts. The routine computes
`wn = 1 / wl`, allocates single-precision output arrays, and skips zero
wavelength padding values.

For Mie it calls:

```text
mie_size_dist_new(distname, 1.0, [Rm, S, 0.001, 100.0],
                  wn[i], RI[i], Dqv=Dqv, /dlm, xres=0.4,
                  Bexttmp, Bscatmp, wtmp, gtmp, SPM, Vavg=Vavg)
```

The Mie distribution routine uses double precision internally. It constructs
the size-parameter quadrature with trapezium quadrature, defaulting to
`xres=0.1` unless overridden; `create_bwgp` overrides this to `0.4`. It limits
the number of points to at least 200 and bounds the modified-gamma integration
radius using `[0.001, 100.0]`. The size parameter is
`Dx = 2*pi*R*wavenumber`.

For STG, the current model file supplies:

- distribution: `modified_gamma`;
- parameters: `Rm=12`, `S=0.1111111` at the effective-radius grid’s STG point;
- integration radius bounds: `0.001` to `100.0` microns;
- scattering code: `mie`;
- mixing ratio: `1.0`;
- RI file: `H2O_Segelstein_1981.ri`.

The modified-gamma calculation derives `alpha=(1-3S)/S`, `b=1/(Rm*S)`,
normalises the size weights, and integrates cross-section-weighted quantities:

```text
W1PA = W1P * pi * R^2
W1PV = W1PA * 4/3 * R
Bext = sum(W1PA * Qext)
Bsca = sum(W1PA * Qsca)
w    = Bsca / Bext
g    = sum(W1PA * g_single * Qsca) / Bsca
Vavg = sum(W1PV)
```

When `Dqv` is supplied, the returned phase matrix is averaged with the same
scattering weighting and divided by `Bsca`. `Phi` receives the F11 phase
function component (`SPM[0,*]`) in the current ORAC call. The downstream
`generate_scattering_properties.pro` expands it with `legpexp` on its angular
quadrature and uses the resulting Legendre moments for DISORT.

The reference wavelength call is explicitly `wl=0.55`; channel calls use the
SRF-selected wavelength arrays. `read_ri.pro` stores wavelength and wavenumber
forms and the current source interpolates RI to the selected spectral points
before calling `create_bwgp`. No POM optics are involved.

## Mie/DLM interface

The repository contains `mie/dlm-code/mie_dlm_single.so`, built from
`mie_dlm_single.c`, `mieint.f`, and `mieintnocmplx.f`. The IDL-visible entry
point is `MIE_DLM_SINGLE`, registered by `IDL_Load`. Its IDL calling convention
is:

```text
MIE_DLM_SINGLE, Dx, Cm, Dqv=Dqv,
                DQext, DQsca, DQbsc, Dg, S1, S2, F11, F33, F12, F34
```

`Dx` and `Dqv` are double arrays; `Cm` is a scalar double complex. For a supplied
angle grid, amplitude/phase outputs are arrays with the C wrapper’s declared
shape `[n_angles, n_sizes]` (IDL’s array ordering must be checked at the ABI
boundary). The exported Fortran symbol is `mieint_`, whose source signature is
declared in `mie_dlm_single.h`. The Fortran routine returns Qext, Qsca, Qbsc,
asymmetry, complex S1/S2, and phase-matrix elements.

The existing `.so` is an IDL plugin, not a standalone Python shared library. It
exports `mieint_`, but it also has unresolved IDL runtime symbols and links
against unavailable Intel runtime libraries (`libifcore.so.5`, `libimf.so`,
`libirc.so`, `libiomp5.so`, and related libraries). A direct `ctypes.CDLL`
load therefore fails in the current environment. No rebuild or external
modification was attempted.

The least invasive future options are to run the existing DLM inside a valid
IDL environment or build a repository-local, separately named wrapper around
the existing Fortran routine after its ABI and compiler-runtime requirements
are agreed. Recompiling the shared legacy DLM is not part of this task.

## Python interface and status

`src/oraclut/optics/legacy.py` defines `OpticalProperties`, the common
downstream representation containing wavelength, effective radius, extinction,
SSA, asymmetry, phase moments, and the 0.55-micron reference wavelength. Its
`legacy_stg_optics` entry point deliberately raises at the unresolved Mie/DLM
boundary; it does not call the user’s POM Mie implementation or invent values.

The surrounding STG configuration and input readers are implemented and tested.
The Mie numerical call, phase/moment production, and numerical comparison are
blocked by the unavailable IDL/Intel-runtime DLM loading path. No trustworthy
standalone STG optical-property array was found in the archive; the full V2
reference LUT remains a downstream structural/scientific target, not a source
for isolating the intermediate Mie arrays.

## Future POM boundary

POM must eventually provide the same `OpticalProperties` quantities as a
separate backend. It must not be imported by the legacy backend, and any
comparison must be labelled diagnostic until the current Mie/DLM result has
been reproduced.
