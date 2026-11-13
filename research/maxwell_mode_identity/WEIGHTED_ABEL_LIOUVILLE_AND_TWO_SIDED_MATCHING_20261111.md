# Weighted Abel-Liouville invariant and two-sided matching for sharp-Maxwell Cowling modes

Date: 2026-11-11  
Status: DERIVED_CONDITIONALLY_AND_NUMERICALLY_CHECKED_IN_SYNTHETIC_RANGE  
Scope: analytic synthetic cold-barotropic phases only. This does not admit any blocked CompOSE crust EOS.

## Result

For two solutions of the same relativistic Cowling polar-mode first-order system on a smooth cold barotropic phase, define

\[
D(r)=W_1V_2-V_1W_2.
\]

Under the conventions and hypotheses below,

\[
\boxed{\Omega(r)=h(r)e^{-\nu(r)}D(r)=\mathrm{constant}},\qquad h=\epsilon+P.
\]

At a sharp Maxwell interface, with state map inner (high density) to outer (low density), \(y_o=T_{i\to o}y_i\),

\[
\frac{\Omega_o}{\Omega_i}
=\frac{h_o}{h_i}e^{-(\nu_o-\nu_i)}\det T_{i\to o}.
\]

For the current no-thin-shell Maxwell interface, \(\nu_o=\nu_i\), and both slow and rapid transfer maps have \(\det T_{i\to o}=h_i/h_o\). Hence the weighted determinant is continuous across either interface. The determinant ratios telescope through a finite chain of such interfaces.

This is an algebraic Abel-Liouville determinant identity for the reduced two-dimensional ODE. It is not a derivation of a canonical covariant symplectic current or a physical energy flux.

## Conventions and assumptions

Use the static spherical metric

\[
ds^2=-e^{2\nu(r)}c^2dt^2+e^{2\lambda(r)}dr^2+r^2d\Omega^2
\]

with signature (-+++). For fixed \(\omega\) and \(\ell\), the project [W,V] radial system is

\[
y'=M(r;\omega,\ell)y,\quad y=(W,V)^T,
\]

\[
M=\begin{pmatrix}
\nu'/\chi & A/\chi-\ell(\ell+1)e^\lambda\\
-e^\lambda/r^2 & 2\nu'
\end{pmatrix},\qquad
A=(\omega r/c)^2e^{\lambda-2\nu}.
\]

Here \(W\) has units m^3, \(V\) m^2, \(\nu,\lambda,A,\chi\) are dimensionless, \(\nu'\) is m^-1, \(h\) is J m^-3, and \(\chi=dP/d\epsilon>0\). The background must satisfy the cold one-parameter First Law and TOV relation

\[
d\epsilon/dP=1/\chi,\qquad P'=-h\nu'.
\]

The two determinant columns must solve the same coefficient matrix: same background, frequency, \(\ell\), sound-speed closure, variable convention and interface prescription. Complex solutions and frequencies are allowed. The determinant is bilinear without complex conjugation.

## Smooth-phase derivation

Differentiating \(D=W_1V_2-V_1W_2\), the off-diagonal terms cancel and

\[
D'=\operatorname{tr}(M)D
=\nu'(2+1/\chi)D.
\]

Cold First Law plus TOV gives

\[
h'=(d\epsilon/dP+1)P'
=-h\nu'(1+1/\chi).
\]

Therefore

\[
\frac{d}{dr}\log(he^{-\nu})
=\frac{h'}h-\nu'
=-\nu'(2+1/\chi),
\]

so \((he^{-\nu}D)'=0\). This is Liouville's formula: for a fundamental matrix F, \((\det F)'=(\operatorname{tr}M)\det F\).

The result is independent of the shared off-diagonal matrix entries, so the \(\ell(\ell+1)\) factor and shared A(r) do not enter the trace. It is not a cross-frequency invariant: different frequencies give different coefficient matrices and off-diagonal cancellation need not hold.

## Interface jump law

For any linear interface map,

\[
D_o=\det(T_{i\to o})D_i,\qquad
\Omega_o/\Omega_i=(h_o/h_i)e^{-(\nu_o-\nu_i)}\det T_{i\to o}.
\]

With \(\alpha=h_i/h_o\) and \(Q=AV+\nu'W\), the current maps are

\[
T_s=\begin{pmatrix}1&0\\(\alpha-1)\nu'/A&\alpha\end{pmatrix},
\qquad
T_r=\begin{pmatrix}\alpha&(\alpha-1)A/\nu'\\0&1\end{pmatrix}.
\]

Both determinants equal \(\alpha\) and both reduce to identity at zero jump. In two dimensions, for \(J=[[0,1],[-1,0]]\), both satisfy \(T^TJT=\alpha J\). If the implementation uses the reverse transfer orientation, the required determinant is the reciprocal. For an interface with a lapse jump, retain the factor \(e^{-(\nu_o-\nu_i)}\); do not assume it away.

## Two-sided matching consequence

The current ell=2 finite-cutoff surface condition uses

\[
b_1=R\nu'(R),\quad b_2=(\omega R/c)^2e^{\lambda(R)-2\nu(R)},
\]

\[
B_s=b_1W(R)/R^3+b_2V(R)/R^2.
\]

A vector spanning the surface boundary line in [W,V] coordinates is, up to nonzero scale,

\[
y_s(R)=(-b_2/R^2,\ b_1/R^3).
\]

Then \(\det[y_c(R),y_s(R)]=B_s\), exactly the unscaled surface boundary linear form. The current normalized real-frequency residual divides this determinant by positive norms; its zero set is unchanged.

Integrate the regular center solution outward and the surface-line solution inward, with the same ODE and same interface transfers. At any admissible matching radius,

\[
E(\omega;r_m)=h(r_m)e^{-\nu(r_m)}
\det[y_c(r_m;\omega),y_s(r_m;\omega)]
\]

is independent of r_m. Its zero condition agrees with the current one-sided surface shooting condition. This is a useful two-sided cross-check. Do not call the normalized angle a holomorphic Evans function: the norm introduces non-holomorphic dependence for complex states. The raw weighted determinant is the relevant scalar for complex-frequency analysis, subject to endpoint-basis and analyticity checks.

Use the same finite pressure cutoff in both directions. The current endpoint is not automatically the exact zero-pressure surface.

## Nonbarotropic and finite-rate caveats

If the background uses an equilibrium derivative \(\chi_{\rm eq}\) but the perturbation matrix uses a frozen-composition derivative \(\chi_{\rm ad}\), then before additional composition terms,

\[
\Omega'/\Omega=\nu'(1/\chi_{\rm ad}-1/\chi_{\rm eq}).
\]

In a general system the Abel factor is \(\exp[-\int\operatorname{tr}M\,dr]\); its connection to h must be re-derived. A finite-rate interface represented by a 2x2 matrix T(ω) preserves this weighted quantity only when \(\det T=(h_i/h_o)e^{\nu_o-\nu_i}\). A reaction variable, dissipation or extra interface degree of freedom may require a larger system.

## Singular endpoints and conditioning

The W,V system is singular at r=0. Apply the identity on open regular intervals, starting from finite r0 with the correct regular expansion. Near a surface where h tends to zero, the unweighted determinant can become ill-conditioned; use a one-sided sub-surface point or regularized variables.

A two-sided synthetic integration showed maximum weighted-determinant spread 3.5761e-4 with DOP853 unlimited maximum step, but 4.6733e-12 with max_step=100 m; a 25 m cap gave about 1.5e-15. This is a test-specific result, not a universal stellar step setting. Step-size convergence and QR/log-amplitude rescaling with exact determinant-scale accounting are mandatory for local work.

Defect interpretation:
- smooth drift: check diagonal signs, radial/q coordinate chain rule, First Law, chi, lapse factor and ODE error;
- localized interface jump: check transfer orientation, determinant, one-sided h and lapse continuity;
- endpoint-only defect: check singular variables and conditioning;
- finite-rate residual after numerical checks: possible physical interface source, requiring a larger-state derivation.

## Numerical characteristic scalar and scale accounting

The raw weighted determinant is the exact matching scalar and retains holomorphic dependence on complex frequency when the ODE coefficients and initial/terminal basis choices are holomorphic. Its magnitude can nevertheless be numerically unusable when the two columns grow or decay at very different rates. The Euclidean-normalized determinant angle is a real-frequency diagnostic only: conjugation in the norm generally destroys holomorphic dependence.

For robust complex-frequency work, use a multiple-shooting matching matrix or a scaled/exterior determinant. If Y=Q R is used, reconstruct det(Y)=det(Q)det(R), including determinant phase and logarithmic amplitude. Record all nonzero column/row scale factors; Q-only tracking loses determinant information. A scale choice based on a non-holomorphic max-component pivot is not an analytic Evans function.

The center is a singular endpoint. Use a verified regular center series at finite r_c, whose leading relation is W+ell*r_c*V=0; do not evaluate the singular RHS at r=0. The current outer condition is imposed at the finite-pressure cutoff and must not be silently replaced by a true zero-pressure surface.

At a mode root the inner concomitant can be zero. Its ratio to the outer concomitant is then undefined. Return a null value and an explicit UNDEFINED_ZERO_DENOMINATOR status, not NaN or infinity. Complex numbers require explicit {re,im} serialization; strict JSON checks must reject non-finite values.

## Numerical evidence from this increment

42 deterministic pytest cases pass (32 test functions: the earlier interface/profile-helper suite contributes 22 collected cases, and the weighted-invariant file contributes 14 test functions expanding to 20 collected cases). Reported checks:
- smooth-phase maximum relative drift: 7.9072e-11, 201 sampled radii, DOP853;
- smooth endpoint relative drift: 4.7232e-11;
- both interface maps preserve the determinant for complex columns;
- multi-interface factors telescope;
- ell=1,2,3,5 exact matrix-exponential trace checks;
- nonbarotropic chi mismatch reproduces the predicted drift;
- a metric-potential jump is detected;
- two-sided matching across five radii: 4.6733e-12 spread at max_step=100 m;
- unlimited-step negative control: 3.5761e-4 spread.

Commands executed:
- PYTHONPATH=<package-root> python -m pytest -q mode_identity/tests
- python -m compileall -q mode_identity
- JSON schema and DAG parse checks.

The machine-readable computation receipt is in computation_evidence/weighted_wronskian_20261111/computation_manifest_v2.json. That run used an exported artifact directory, not a Git checkout, so its manifest explicitly says the source commit could not be authenticated in that local run. GitHub branch/source identity is recorded separately.

## Literature boundary

Finn (1988), Relativistic stellar pulsations in the Cowling approximation, MNRAS 232, 259, DOI 10.1093/mnras/232.2.259, supports the near-Sturm-Liouville character of the relativistic Cowling problem under its assumptions.

Tonetto & Lugones (2020), Discontinuity gravity modes in hybrid stars, Phys. Rev. D 101, 123029, DOI 10.1103/PhysRevD.101.123029, derives separate slow and rapid sharp-interface conditions.

Dong, Kuan & Yunes (2026), Hybrid-star asteroseismology across sharp phase interfaces, arXiv:2609.33002, is a full-GR contextual source for two-sided basis matching and interface-local tangential displacement. Its formulation is not identical to this Cowling system.

Townsend & Teitler (2013), GYRE: An open-source stellar oscillation code based on a new Magnus Multiple Shooting Scheme, MNRAS, DOI 10.1093/MNRAS/STT1533, provides an established stellar-oscillation example of multiple shooting and discriminant-based eigenfrequency finding. It motivates the numerical architecture, but does not validate this solver's interface matching or this particular weighted concomitant.

Townsend, Kuenzi & Christensen-Dalsgaard (2025), Tools for Characterizing the Numerical Error of Stellar Oscillation Codes, ApJ, DOI 10.3847/1538-4357/addf32, motivates explicit error characterization beyond a single tolerance setting. It does not provide a universal max_step or error threshold for the current Cowling fixture.

No inspected source explicitly states the precise Omega identity derived here. Treat it as an independent derivation for the stated equations; do not attribute it to these papers.

## Final disposition

- Weighted identity: DERIVED_CONDITIONALLY_AND_NUMERICALLY_CHECKED_IN_SYNTHETIC_RANGE.
- Two-sided matching: DERIVED_AND_SYNTHETICALLY_CHECKED; actual fixture roots not yet processed.
- Canonical covariant symplectic current: UNRESOLVED_NOT_CLAIMED.
- Physical mode labels: UNRESOLVED.
- C7X source blocker: unchanged.
