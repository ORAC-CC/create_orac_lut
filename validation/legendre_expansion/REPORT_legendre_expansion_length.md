# Legendre expansion length for ORAC particle phase functions

> Moved on 2026-10-02 from `validation/` to `validation/legendre_expansion/`.
> This is the investigation record as it stood before the Phase 3
> implementation; see `validation/REPORT_lut_numerics_development.md`.

Date: 2026-10-01, third revision. This revision tests the termination method of
Grainger (1990), §4.4–4.5, and replaces the earlier versions of this report. Their
results are summarised in §0.

No production code, run file or LUT was changed. All calculations used the
production Python routines outside the repository.

Thesis source: `documents/Grainger 1990.pdf`. Thesis pp. 53–57 are PDF pp. 70–73.

## 0. What was already established

- `nmom` (1000 in every V23 run file and as the function default) sets three
  things at once: the number of Gauss–Legendre nodes N_q, the number of stored
  coefficients, and the moment count DISORT receives.
- At short wavelengths with large particles, N_q = 1000 aliases even the low
  moments. TOA radiances are wrong by up to 13% (liquid, r_e 39 µm) and 68% (ice
  spheres, r_e 93 µm) at 0.645 µm, and by more at 0.41 µm.
- DISORT's Nakajima–Tanaka correction sums all supplied moments, so truncating an
  accurate series early causes Gibbs ringing at the ORAC scattering angles
  (30–180°).
- **Correction:** the earlier versions attributed the coefficient "noise floor"
  (|ω_l| ≈ 10⁻⁴) to the Mie calculation. It actually comes from errors in numpy's
  `leggauss` weights at N = 8000 (§3). It affected only the high-order
  *references*, and none of the earlier conclusions change.

## 1. The thesis method

These quotations are from the rendered pages, in the thesis notation.

**Expansion.** p(μ) = Σ_{l=0}^{L} ω_l P_l(μ), with ω_l = (2l+1)/2 ∫₋₁¹ p P_l dμ
[4.4.1–2]. Then χ_l = ω_l/(2l+1) [4.4.4], and g = χ₁ = ω₁/3.

**Termination (p. 56).**
- Wiscombe: "a rough value for L such that χ_l < 10⁻⁶ for l ≥ L".
- King (1983): "the (stricter) arbitrary criterion for L such that ω_l < 10⁻⁹ for
  l ≥ L". The thesis then states: "The series was terminated using King's
  criterion."

The tested quantity is therefore **ω_l, the coefficient including the factor
(2l+1)**, not χ_l. It must hold for **all l ≥ L**, not at a single order. The
thesis writes ω_l, not |ω_l|; this report uses the magnitude, since a negative
coefficient cannot by itself signal convergence.

**Accuracy tests (p. 56).**
- "it was ensured that ω₀ was equal to unity and ω₁/3 was equal to g, calculated
  directly from Mie theory";
- "the reconstruction of phase function from the Legendre series was always
  accurate to 6 significant figures";
- following Wiscombe, the χ_l "are positive and decrease monotonically. When these
  rules are violated it is prima facie evidence that the moments in question are
  inaccurate."

**Quadrature (pp. 54–55).**
- "If the phase function is a polynomial of degree L the integrand … is at most a
  polynomial of degree 2L. As Gauss–Legendre quadrature of order N is exact for
  polynomials of degree less than 2N, the order of the quadrature formula must
  exceed the number of significant Legendre coefficients."
- The thesis itself used Lobatto quadrature in θ on [0°, 180°] [4.4.5], for which
  "N < L succeeds" (Table 4.4.1), and for which ω₀ "is an excellent error
  monitor".
- A single drop needs "2x + a few" terms; a size distribution is much smoother.

## 2. ORAC implementation compared with the thesis

**Phase function.** Mie F11 averaged over the modified-gamma size distribution
(`src/oraclut/idl_mirror/create_bwgp.py:67-145,232`; IDL
`create_orac_lut/create_bwgp.pro:86`). It is sampled on N_q = `nmom` Gauss nodes in
μ, numpy `leggauss`
(`src/oraclut/idl_mirror/generate_scattering_properties.py:262-265`).

**Coefficients.** `legpexp` computes lc_n = (2n+1)/2 Σ w p P_n
(`create_bwgp.py:148-179`; IDL `create_orac_lut/legpexp.pro`). That is exactly the
thesis ω_n, with ω₀ = 1. The code then forms χ_n = lc_n/(2n+1), the DISORT `PMOM`
(`generate_scattering_properties.py:372,396,409`). DISORT receives `NMOM = nmom−1`
(`setup_disort.py:25,34`) and resets `PMOM(0)=1`
(`create_orac_lut/disort2/src/DISORTfunctions.f:2824`).

**`legpexp` stop (`create_bwgp.py:173`) against King's criterion.**

| | King (thesis) | `legpexp` |
|---|---|---|
| Normalisation | ω_l | ω_l (identical) |
| Threshold | 10⁻⁹ | 10⁻⁹ (identical; the IDL header comment's "E-5" is stale) |
| Test | ω_l < 10⁻⁹ for **all** l ≥ L | \|ω_n\| < 10⁻⁹ at the **first single** n |
| Can it exceed N_q? | yes (L is chosen first) | no (loop bounded by N_q) |
| Does it set L in ORAC? | – | **no**: `inlc` is discarded and DISORT always gets `nmom` moments |

The routine (described as converted from Grainger's Fortran `Alegpexp`) therefore
contains King's threshold. It applies the threshold as a single-coefficient test
and has no effect on the result.

## 3. Numerical test

**Cases.**
- Production `water-ice_sph`: ice spheres, Warren 2008 refractive index.
- Modified gamma with v_e = 0.111, integrated over 0.001–100 µm, so
  x_max = 2π·100/λ.
- References: N_q = 8000 nodes, 8000 moments.

**Gauss weights.** Newton-refined weights had to replace numpy `leggauss`.
numpy's weights have relative errors of 2.8e-6 at N = 8000 (2.3e-7 at 3000,
1.5e-8 at 2000, 8e-9 at 1000), largest at the node nearest θ = 0. Multiplied by
the 6.7×10⁵ forward peak, they produce a constant spurious tail |ω_l| ≈ 6×10⁻⁴.

**Floor with refined weights.** Re-expanding an exactly band-limited polynomial of
degree 2·NSTOP(x_max) = 3148 still leaves |ω_l| ≈ 1–3×10⁻⁸ out to l ≈ 5200.
This is the double-precision floor of the expansion for this dynamic range,
p(0)/p_min ≈ 5×10⁷.

**Table 1. King's criterion applied to accurate coefficients.**

| Case | x_max | Last l with \|ω_l\| ≥ 10⁻³ / 10⁻⁶ | **King L** (\|ω_l\| < 10⁻⁹ for all l ≥ L) | ω₀−1 | ω₁/3 − g | Reconstruction max rel. err. (all angles) | Wiscombe L (χ < 10⁻⁶) |
|---|---|---|---|---|---|---|---|
| A ice, r_e 93, 0.412 µm (worst) | 1526 | 3061 / 3097 | **5223 computed; true value 3097–3148** | 3e-12 | 3e-12 | 8e-10 (9 s.f.) | 3047: 5e-3 (2.3 s.f.) |
| B ice, r_e 41, 0.412 µm | 1526 | 2995 / 3068 | **3092** | 2e-12 | 2e-12 | 3e-8 (7.5 s.f.) | 2929: 5e-3 |
| F ice, r_e 93, 0.6455 µm | 973 | 1956 / 1981 | **1995** | 2e-12 | 2e-12 | 7e-8 (7.2 s.f.) | 1949: 4e-3 |

Case A's "true value" range follows because the tail from 3148 to 5223 is the
numerical floor described above (≈ 2×10⁻⁸ > 10⁻⁹), not phase-function content.

**Observations.**
- **King's criterion gives physically sensible limits.** L ≈ 2·x_max + a few (the
  thesis's own "2x + a few", with x set by the largest radius integrated). L scales
  as 1/λ. It barely depends on r_e because ORAC always integrates to 100 µm.
- **King-limited series pass the thesis's 6-significant-figure reconstruction
  test; Wiscombe's χ < 10⁻⁶ rule fails it** (about 2 s.f.).
- **King's L converges ORAC observables.** In the earlier DISORT test, L = 3250
  converged TOA reflectance for A to 0.007% and L = 3000 left 9.8%; F converged
  at L = 2000.
- **10⁻⁹ is at the limit of double-precision expansion for the worst case.** For
  A it is not resolvable, and the computed L is conservative (about +70%) but
  harmless.

**Thesis accuracy tests in the ORAC (Gauss-in-μ) setting.**
- ω₀ = 1 and ω₁/3 = g detect only missing forward-peak energy. Production
  N_q = 1000 for A gives ω₀ = 0.794. Aliasing of the higher moments goes
  undetected: with N_q = 3000 (refined weights), ω₀−1 = 3×10⁻¹¹ while the moments
  are wrong by up to 0.15. The thesis's claim that ω₀ is an excellent monitor
  applies to its Lobatto-in-θ quadrature and does not transfer.
- The positivity/monotonicity rule is violated by accurate moments: χ_l rises at
  27–50 orders (l ≈ 9–20, and in the last ≈ 20 orders before the cut-off), and the
  sign alternates near the cut-off. As stated, it would wrongly flag these
  functions.
- **Only the reconstruction test is a reliable check here.**

## 4. Quadrature order needed for L moments

With Gauss–Legendre in μ, moment l is exact when 2N_q − 1 ≥ D + l, where
D ≈ 2·NSTOP(x_max) is the degree of p. The first L moments therefore need

  **N_q ≥ (D + L + 1)/2**, and for L = L_King ≈ D this means **N_q > L**,

as the thesis states. Verified for A:

| N_q | Moments exact (\|Δχ\| ≤ 1e-6) up to l | Predicted (2N_q − 1 − D) |
|---|---|---|
| 2000 | 949 | ≈ 935–950 |
| 3000 | 2954 | ≈ 2850–2950 |
| 3100 | all l < 3100 (\|Δχ\| ≤ 3e-7) | – |

N_q and L are therefore distinct quantities. The current code makes them equal
only by construction (`legpexp` returns N_q coefficients). In this scheme N_q
must slightly exceed the required L, so separating them saves little.

The thesis's Lobatto-in-θ quadrature [4.4.5] can succeed with N_q < L, and could
in principle reduce the Mie cost. This was not tested.

## 5. Answers

1. **Thesis method.** King's criterion: terminate at L such that ω_l < 10⁻⁹ for
   all l ≥ L, with ω_l = (2l+1)χ_l. It is checked by ω₀ = 1, ω₁/3 = g and
   reconstruction to 6 significant figures, with Lobatto quadrature in θ.
2. **Direct applicability.** Yes for the criterion itself: the current ORAC ω_l
   use the same normalisation, so they can be tested directly. Two conditions
   apply:
   - the coefficients must be computed with N_q > L and accurate Gauss weights;
   - the ω₀/g and monotonicity tests do not transfer to ORAC's Gauss-in-μ
     scheme, so reconstruction against directly computed P(Θ) must be the
     accuracy check.
3. **Sensible limits.** Yes: L ≈ 3100 (0.412 µm) and ≈ 2000 (0.645 µm), equal to
   2·x_max + a few, and these converge DISORT TOA radiances. For the worst case,
   10⁻⁹ sits at the double-precision floor and over-estimates L harmlessly.
4. **Quadrature.** N_q ≥ (D + L + 1)/2, i.e. N_q slightly greater than L, with
   accurate weights. The current N_q = 1000 is below the D required whenever
   x_max ≳ 500.
5. **NMom = 1000.** Not scientifically justified. The limit should follow from
   the phase function. King's criterion supplies L, and the size parameter
   (D ≈ 2·NSTOP(x_max)) tells you in advance how large N_q must be to compute
   those moments.
6. **For a future production task** (not done here):
   - choose N_q per phase function from x_max, so that N_q > D;
   - compute ω_l with accurate Gauss weights (refine or replace numpy `leggauss`
     above about 2000 nodes);
   - terminate with King's criterion as a sustained test, "for all l ≥ L", not a
     single-coefficient test;
   - pass that L to DISORT instead of discarding `inlc`;
   - verify by reconstruction against direct Mie P(Θ) over 30–180°;
   - lift the 2000-node cap (`create_bwgp.py:45`; `quadrature.pro:106`);
   - check the Baum (non-spherical) models separately before relying on the
     criterion for them;
   - validate any change as a deliberate departure from the legacy reference.
