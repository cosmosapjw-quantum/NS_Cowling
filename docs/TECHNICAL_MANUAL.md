# Technical manual: compose_eos 3.3.0 + relativistic Cowling solver

## Architecture

The integrated project contains four linked but separately versioned layers.

1. EOS generator: compose_eos 3.3.0 parses CompOSE data, checks composition and cold thermodynamics, reconstructs master-function gaps and joins crust/core branches. Its raw inputs, normalized tables, matching metadata, 390-case matrix and audit evidence are preserved under eos_generator.
2. Intake: the portable pipeline reads only successful matched EOSs using the bundled compose_eos API. It checks the cold first-law identity and strict monotonicity before exporting source-knot CSV files in MeV fm^-3 and SI columns.
3. Default solver EOS provider (v3): monotone PCHIP of g=epsilon-P versus ln(P), exact preservation of source P and epsilon knots, with n_B reconstructed by the cold first law and mu_B from Euler closure. The optional q-linear P(mu_B) adapter is a v2 legacy comparison method only.
4. Stellar perturbations: the supplied Fortran implementation integrates TOV and ell=2 Cowling fluid equations, locates f and p1 roots by two-sided shooting/matching, and identifies modes by eigenfunction node counts. The old 234-EOS spectrum is a v2 q-linear cohort; one representative direct-provider v3 spectrum is separately recorded.

## Physical conventions and units

Metric signature is (-,+,+,+). The background metric is

\[
ds^2=-e^{2\nu(r)}c^2dt^2+e^{2\lambda(r)}dr^2+r^2d\Omega^2.
\]

The Fortran stellar kernel uses SI: r in m, m(r) in kg, P and epsilon in J m^-3, baryon density in m^-3, chemical potential in J, and omega in s^-1. The CompOSE preparation and table-adapter stages retain a companion representation in MeV, fm^-3 and MeV fm^-3. Conversion factors are 1 MeV=1.602176634e-13 J, 1 fm^-3=1e45 m^-3 and 1 MeV fm^-3=1.602176634e32 J m^-3.

The static TOV equations are

\[
m'=4\pi r^2\epsilon/c^2,\qquad
P'=-(\epsilon+P)\nu',
\]
\[
\nu'=\frac{G[m+4\pi r^3P/c^2]}
{c^2r^2[1-2Gm/(rc^2)]}.
\]

The lapse is normalized to the Schwarzschild exterior at the finite pressure cutoff. This sets the frequency time coordinate for the truncated model.

## Polar-fluid Cowling equations

The perturbations have time dependence exp(i omega t), and

\[
\xi^r=e^{-\lambda}\frac{W(r)}{r^2}Y_{\ell m},\qquad
\xi^a=-\frac{V(r)}{r^2}D^aY_{\ell m}.
\]

The adopted model neglects all Eulerian metric perturbations. For ell=2, define

\[
Q=\left(\frac{\omega r}{c}\right)^2e^{\lambda-2\nu}V+\nu'W.
\]

The radial equations are

\[
\frac{dW}{dr}=\frac{Q}{\chi}-6e^\lambda V,\qquad
\frac{dV}{dr}=2\nu'V-\frac{e^\lambda}{r^2}W,
\]
where chi=dP/d epsilon is the equilibrium beta-equilibrated barotropic response. The regular central solution behaves as W=A r^3+O(r^5) and V=-A r^2/2+O(r^4). The Lagrangian pressure perturbation is Delta P=-(epsilon+P)e^-lambda Q/r^2. Since the tabulated EOS stops at positive pressure, the surface residual is a finite-cutoff diagnostic, not an exact vacuum boundary.

The calculation uses two-sided integration and a normalized matching determinant. Root location and mode identity are separate gates. The f mode is selected by zero interior nodes in W and V; p1 by one interior node in each. The search interval is finite and does not prove completeness of the entire spectrum.

## EOS construction and interpolation contract

The upstream master-function reconstruction and the q-linear solver adapter are different mathematical layers. The upstream EOS generator defines thermodynamic knots and matching bridges. The solver's native representation is a cubic Hermite P(mu_B), with P'=n_B and epsilon=mu_B n_B-P. Its cellwise stability and causality gate is applied analytically throughout each interval. The raw knot first-law identity cannot guarantee that a selected interpolating polynomial has no interior overshoot.

The adapter selects q=n_B/mu_B piecewise linear. It integrates n_B=mu_B q exactly to pressure and derives energy by the Euler relation. The resulting P(mu_B) is exactly cubic per cell and matches the Fortran solver's Hermite interface. The formula c_s^2=q/(q+mu_B q') proves 0<c_s^2<=1 whenever q>0 and q'>=0. Full details and cohort-level changes are in EOS_TO_COWLING_PIPELINE.md and the adapter receipts.

## Mode results and numerical checks

The completed cohort comprises 234 unique EOSs: 78 each for BBB_BHF2, RG_SLY4 and GM_GM1. The target mass is 2.644606257493485e30 kg, or 1.3300106263 nominal solar masses. All retained f roots have node counts (N_W,N_V)=(0,0), and all retained p1 roots have (1,1). The maximum relative target-mass residual is 9.9143e-13. The maximum absolute normalized matching residuals are 4.5603e-8 for f and 3.0892e-10 for p1.

The three-family median frequencies are:
- BBB_BHF2: f=2.5051 kHz, p1=7.0083 kHz.
- RG_SLY4: f=2.3698 kHz, p1=6.7177 kHz.
- GM_GM1: f=1.9202 kHz, p1=6.0090 kHz.

Representative loose-to-tight integration changes reach 9.35e-8 Hz for f and 1.71e-6 Hz for p1. These are representative numerical convergence diagnostics, not global error bounds or estimates of physical accuracy.

## Reuse, reruns and provenance

The integrated release contains source-knot and adapter tables, the original seminar package, the compose_eos 3.3.0 generator release, all mode result tables, the REVTeX paper, source references, test fixtures and audit records. Every adapter receipt links a source hash to an adapter hash. The verifier checks the saved 234-row spectrum and its node identities without running TOV or Cowling integrations. No automatic calculation is triggered by installation or by the verification command.

A new physical closure, a different interpolation representation, frozen-composition response, full-GR perturbations, a sharp interface, or a zero-pressure surface requires a new output directory and a distinct model identifier. Do not mix results from such a model with this beta-equilibrated smooth-barotrope v1 cohort.
