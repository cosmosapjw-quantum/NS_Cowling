# Conditional compensated Abel–Liouville concomitant for closure mismatch and sharp interfaces

Date: 2026-11-11  
Research id: COWLING-COMPENSATED-ABEL-INTERFACE-v1  
Scope: specified analytic 2×2 relativistic Cowling radial system only; synthetic validation; no source EOS or stellar-mode claim.

## 1. Claim and scope gate

Let Y(r)=[y1(r),y2(r)] be a 2×2 fundamental matrix whose columns ya=(Wa,Va)^T solve the same first-order system

\[
Y'=M Y,\qquad
M=\begin{pmatrix}
\nu'/\chi_{\rm ad} & A/\chi_{\rm ad}-\ell(\ell+1)e^\lambda\\
-e^\lambda/r^2 & 2\nu'
\end{pmatrix},\qquad
A=(\omega r/c)^2e^{\lambda-2\nu}.
\]

Here prime means d/dr, h=epsilon+P>0, and chi_ad is the derivative entering this perturbation matrix. The background is assumed to obey the cold equilibrium relation

\[
P'=-h\nu',\qquad
\frac{d\epsilon}{dP}=\frac1{\chi_{\rm eq}},\qquad
h'=-\left(1+\frac1{\chi_{\rm eq}}\right)h\nu'.
\]

The exact claim is conditional on this matrix and background first law. It does not assert that a generic nonbarotropic perturbation system is obtained by replacing only chi_eq with chi_ad. Composition-gradient terms, extra variables, dissipation, or a different perturbation system require a fresh trace calculation.

## 2. Smooth-phase derivation

Define the bilinear determinant D=det Y=W1 V2−V1 W2, with no complex conjugation. Abel–Liouville gives

\[
D'=(\operatorname{tr}M)D
=\nu'\left(2+\frac1{\chi_{\rm ad}}\right)D.
\]

The background first law gives

\[
\frac{d}{dr}\ln(he^{-\nu})
=\frac{h'}h-\nu'
=-\nu'\left(2+\frac1{\chi_{\rm eq}}\right).
\]

Therefore,

\[
\Omega:=he^{-\nu}D,\qquad
\Omega'=\delta(r)\Omega,\qquad
\boxed{\delta(r)=\nu'(r)\left(\frac1{\chi_{\rm ad}(r)}-\frac1{\chi_{\rm eq}(r)}\right)}.
\]

This multiplicative identity remains valid when D=0; only log-derivative or ratio diagnostics become undefined. The dimensions are consistent: chi is dimensionless, nu' has dimension inverse length, and delta dr is dimensionless.

For a fixed reference radius r0 within the same smooth phase, define

\[
I(r)=\int_{r_0}^{r}\delta(s)\,ds,\qquad C(r)=e^{-I(r)}.
\]

Then C(r) Omega(r) is constant on that smooth phase. In the equal-closure limit chi_ad=chi_eq, delta=0, C=1, and the prior barotropic concomitant is recovered.

## 3. Interface jump and global piecewise invariant

At an interface rk, use explicitly one-sided values and the stated inner-to-outer orientation

\[
y_o=T_k y_i,\quad h_i,h_o>0,\quad \nu_i,\nu_o.
\]

Since det(Yo)=det(Tk)det(Yi),

\[
\boxed{J_k:=\frac{\Omega_o}{\Omega_i}
=\frac{h_o}{h_i}e^{-(\nu_o-\nu_i)}\det T_k}.
\]

For all interfaces crossed outward from the fixed reference point, ordered by increasing radius, define

\[
\boxed{
\mathcal K(r)=
\Omega(r)\exp\left[-\int_{r_0}^{r}\delta(s)\,ds\right]
\prod_{r_k\in(r_0,r)}J_k^{-1}
}.
\]

It is constant on each smooth segment because the exponential removes the bulk drift. It is continuous across each interface because the factor Jk^-1 removes the exact one-sided determinant jump. Each crossed interface must appear exactly once. For inward integration, use reciprocal jump factors and reverse the ordered path consistently.

For the current slow/rapid matrices, det(Tk)=hi/ho. If the lapse is continuous, Jk=1; the old result is recovered across those interfaces even when a closure mismatch exists in the smooth bulk, after applying the bulk correction. A lapse jump must not be silently discarded.

The interface rule is purely algebraic for a supplied invertible 2×2 transfer matrix. It does not derive a new finite-reaction or dissipative boundary law. If a transfer matrix is singular, the inverse-factor formulation is undefined and the calculation must stop or be reformulated as a higher-dimensional matching problem.

## 4. Numerical representation and strict output contract

The implementation separates four operations:

1. closure_drift_density: evaluate delta(r).
2. integrated_closure_drift: quadrature of delta on one smooth segment, with strictly increasing radii.
3. interface_multiplier: calculate Jk from one-sided h, nu, and det(Tk).
4. compensate_concomitant: calculate K=Omega exp(-I)/product(Jk).

All input arrays must be finite; both sound-speed closure quantities and enthalpy densities must be strictly positive. Interface determinants and multipliers must be nonzero. Complex determinant phases are retained, including small imaginary parts near the real-frequency axis. JSON does not support complex numbers: serialize as an explicit object { "re": ..., "im": ... } and use strict JSON (allow_nan=false). Do not serialize NaN/Infinity and do not replace unavailable values with zero.

The numerical quadrature error for I, ODE integration error for D, one-sided interface evaluation error, and floating-point product/overflow error are distinct. Store them separately. For many interfaces, accumulate logarithmic magnitudes and complex phases rather than multiplying huge/small factors naively; never discard phase. At a zero determinant, the multiplicative identity remains meaningful but log-derivative and observed-ratio diagnostics are undefined.

## 5. Synthetic verification performed

The implementation was tested on constructed coefficient matrices with the specified trace, not on a stellar model.

- 12 deterministic tests pass.
- Constant-coefficient closure-drift quadrature agrees with the analytic integral to absolute error 1.734723475976807e-18.
- A two-segment matrix-exponential test crosses a nontrivial transfer matrix and a lapse discontinuity; the compensated concomitant agrees with its initial value to relative tolerance 2e-12.
- DOP853 independently reproduces the predicted nonbarotropic drift and its compensated constant with relative tolerance 2e-9.
- Multiple arbitrary non-singular interface factors telescope.
- A test first failed because an implementation discarded an imaginary component of order 1e-18; after type-preserving correction, the regression passes.
- Wolfram symbolic reductions returned zero for the bulk compensation residual and for the slow/rapid transfer determinant residuals. With the current det(T)=hi/ho and continuous lapse, the standard interface multiplier residual is also zero.

The actual synthetic run and its exact command are in COMPENSATED_ABEL_RESULTS.json and COMPENSATED_ABEL_VALIDATION_REPORT.md.

## 6. Independent diagnostic integrity finding

The current diagnostics.py on the audited PR #6 snapshot accepts complex phase profiles and uses complex overlaps, but returns interface jumps via float(real(outer-inner)). This keeps only the real projection. A direct value-level example Delta W=-0.75-1.25i is emitted as -0.75, losing an imaginary component of magnitude 1.25.

This is a diagnostic serialization defect, not evidence of a wrong mode equation. The local agent must add a test first and preserve backward compatibility by retaining the legacy real-projection key while adding explicit real, imaginary, and absolute-magnitude fields for both W and V. No phase should be assigned physical meaning without a stated global phase convention. Strict JSON serialization must be tested.

## 7. Literature boundary

Rau & Salaben (2023), arXiv:2309.08540, study radial stability of hybrid stars with first-order transitions while allowing the adiabatic sound speed to differ from the equilibrium sound speed. They also derive a reactive phase-boundary condition, with a stated nonbarotropic scope. This supports treating the two sound-speed closures as distinct physical inputs, but it does not establish the present nonradial Cowling trace formula or license importing their reactive boundary law into this 2×2 problem.

Dong, Kuan & Yunes (2026), arXiv:2609.33002, study complex quadrupolar f- and p1-spectra in full general relativity across sharp interfaces and compare slow/rapid conversion. Their equations and metric perturbations differ from the present Cowling system. Their result is contextual evidence that interface-sensitive spectral information can be important, not a proof of this determinant identity.

Primary records:
- Rau & Salaben: https://arxiv.org/abs/2309.08540
- Dong, Kuan & Yunes: https://arxiv.org/abs/2609.33002
- Tonetto & Lugones (2020), discontinuity gravity modes: https://doi.org/10.1103/PhysRevD.101.123029

No inspected source is cited as proving K. The boxed expressions above are derived here from the explicitly stated ODE trace, background first law, and one-sided transfer determinant.

## 8. Evidence classification and non-claims

- DERIVED: the smooth drift and piecewise compensated identity under the stated ODE/background hypotheses.
- SYMBOLICALLY_CHECKED: the specified drift cancellation and current slow/rapid determinant formulas.
- NUMERICALLY_CHECKED: the finite synthetic matrix-exponential and DOP853 assertions recorded above.
- NOT IMPLEMENTATION_VERIFIED: native Fortran q-to-r equivalence and native profile export.
- UNRESOLVED: whether a general nonbarotropic physical perturbation system retains this exact two-variable matrix after composition terms are included.
- NOT CLAIMED: canonical covariant symplectic current, finite-rate/dissipative interface physics, actual EOS results, actual four fixture-root profiles, or physical f/p1/g labels.

The existing C7X D3–D6 source blocker remains unchanged. No raw EOS, pair manifest, production adapter, source-authority code, TOV result, or source-backed Cowling result is touched.
