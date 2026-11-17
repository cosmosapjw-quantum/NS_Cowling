# Sharp Maxwell transition in TOV and Cowling: phase-aware reference path v1

## 1. Scope and status

This increment connects the existing native one-conserved-charge Maxwell EOS constructor to a phase-aware stellar-structure and Cowling reference integrator. It is deliberately separate from the legacy smooth-table Fortran TOV/Cowling path: a phase-labelled two-sided EOS table is never passed through a single-valued interpolator across the density jump.

The EOS construction remains native Fortran (`construct_maxwell_transition`). The new stellar-structure and oscillation integration is implemented in Python/SciPy as an independently coded reference path. The Fortran code also contains independent slow- and rapid-conversion interface-transfer routines, exercised through the compiled `case_driver`. The existing production Fortran `solve_tov_star` and `solve_cowling_spectrum` still accept smooth single-valued EOS tables only. Therefore this is **not yet a native Fortran sharp-interface solver**, nor a complete stable-branch implementation.

Supported model: cold, beta-equilibrated, one globally conserved baryon charge, caller-confirmed local electric neutrality in each phase, spherical nonrotating perfect-fluid TOV background, barotropic quadrupolar (`ell=2`) Cowling perturbations, one sharp Maxwell transition, and finite positive surface-pressure cutoff. Both slow and rapid phase-conversion limits are supported at the perturbative interface. A microscopic finite-rate conversion law is not represented.

## 2. Maxwell thermodynamic contract

The existing native builder provides two phase-labelled branches. `Iphase=1` denotes the low-density phase and `Iphase=2` the high-density phase. At the transition,

\[
 P_- = P_+ = P_t,\qquad
 \mu_{B,-}=\mu_{B,+}=\mu_t,
\]

while

\[
 n_+>n_-,\qquad \epsilon_+>\epsilon_-.
\]

At zero temperature the latent-energy identity is

\[
 \boxed{\epsilon_+-\epsilon_-=\mu_t(n_+-n_-).}
\]

Here `-` and `+` denote low- and high-density phases, not the inner/outer sides of a star. In a star with a high-density core, the inner side is the `+` phase and the outer side is the `-` phase.

`MaxwellEOS` splits the labelled table into two separate `DirectBarotrope` providers. It validates transition-pressure equality, chemical-potential continuity, positive density/energy jumps, and the latent-energy identity. No PCHIP or derivative is formed across the jump. The provider's two phase-local first-law residuals are retained in `interface_thermo`.

The caller must set `assume_local_charge_neutrality=True` when invoking the EOS constructor. A composition-free barotrope does not prove local charge neutrality. Composition-bearing EOS tables and multi-charge Gibbs constructions are outside this path.

## 3. Phase-segmented TOV integration

The independent TOV solver uses

\[
 \frac{dm}{dr}=4\pi r^2\frac{\epsilon}{c^2},\qquad
 \frac{dP}{dr}=-(\epsilon+P)\nu',
\]

\[
 \nu'=\frac{G}{c^2}\frac{m/r^2+4\pi rP/c^2}
 {1-2Gm/(rc^2)}.
\]

It integrates in the monotone coordinate

\[
 q=\ln(P_{\rm ref}/P),\qquad
 \frac{dr}{dq}=-\frac{P}{dP/dr},
\]

rather than allowing an adaptive Runge--Kutta stage to query an EOS on the wrong side of a pressure event. The high-density segment is integrated from the centre to the exact endpoint `P=Pt`; the integration is then restarted from the same `(r,m,nu)` with the low-density provider and continued to the finite surface cutoff. The phase transition is not smoothed or sampled through.

At the interface the implementation checks continuity of

\[
 r,\quad m,\quad \nu,\quad \lambda=-\tfrac12\ln(1-2Gm/rc^2),\quad \nu'.
\]

The thermodynamic values `epsilon`, `n_B`, and generally `chi=c_s^2/c^2` are one-sided and may jump. `P` and `mu_B` are continuous. The TOV pressure derivative is one-sided because it depends on `epsilon+P`.

The Seidov onset ratio is recorded as a rapid-conversion stability diagnostic,

\[
 \frac{\Delta\epsilon}{\epsilon_-},\qquad
 \left(\frac{\Delta\epsilon}{\epsilon_-}\right)_{\rm crit}
 =\frac12+\frac{3P_t}{2\epsilon_-}.
\]

This diagnostic does not replace a full maximum-mass and stable-branch search. No global stability claim is made for an integrated configuration solely because the ODE solver returned a solution.

## 4. Cowling system and interface conditions

The implementation uses the same displacement convention as the smooth Cowling solver,

\[
 \xi^r=e^{-\lambda}\frac{W}{r^2}Y_{2m}e^{i\omega t},\qquad
 \xi^\theta=-\frac{V}{r^2}\partial_\theta Y_{2m}e^{i\omega t},\qquad
 \xi^\phi=-\frac{V}{r^2\sin^2\theta}\partial_\phi Y_{2m}e^{i\omega t}.
\]

Define

\[
 A=\left(\frac{\omega r}{c}\right)^2e^{\lambda-2\nu},\qquad
 Q=AV+\nu'W,\qquad h=\epsilon+P.
\]

The first-order Cowling equations are

\[
 \frac{dW}{dr}=\frac{Q}{\chi}-6e^\lambda V,\qquad
 \frac{dV}{dr}=2\nu'V-e^\lambda\frac{W}{r^2}.
\]

The Lagrangian pressure perturbation is

\[
 \Delta P=-h e^{-\lambda}\frac{Q}{r^2}.
\]

At a material phase interface, pressure continuity requires `[Delta P]=0`. Since `r`, `lambda`, and `nu'` are continuous,

\[
 Q_{\rm out}=\frac{h_{\rm in}}{h_{\rm out}}Q_{\rm in}.
\]

The phase-conversion law supplies the second condition.

**Slow conversion** (interface matter does not cross the phase boundary):

\[
 [\xi^r]=0\quad\Longleftrightarrow\quad W_{\rm out}=W_{\rm in}.
\]

**Rapid conversion** (phase conversion keeps up with the perturbation):

\[
 \left[\xi^r-\frac{\Delta P}{P'_0}\right]=0,\qquad
 P'_0=-h\nu',
\]

which in this convention becomes

\[
 W_{\rm out}-\frac{Q_{\rm out}}{\nu'}
 =W_{\rm in}-\frac{Q_{\rm in}}{\nu'}.
\]

After either condition determines `W_out`, the second displacement variable follows from

\[
 V_{\rm out}=\frac{Q_{\rm out}-\nu'W_{\rm out}}{A}.
\]

Both Python and the independent Fortran `ns_interface_transfer` module implement these limits. The Fortran `material_interface_transfer` is the slow limit; `rapid_interface_transfer` is separate and must not be conflated with it.

The Cowling integration is split at the exact interface: integrate the high-density core to `r_t`, apply the selected transfer map, then integrate the low-density envelope to the surface. The existing finite-cutoff surface residual is reused. A frequency scan returns **unlabelled roots**; node counts are sampled diagnostics only. General mode-family tracking through strong transitions, including possible interface modes, is not yet claimed.

## 5. Independent validation evidence

The focused test file is `python_frontend/tests/test_maxwell_star.py`. It uses analytic cold phases

\[
 P=Kn_B^2,\qquad \epsilon=m n_B+P,\qquad
 \mu_B=m+2Kn_B,
\]

for which the Maxwell point can be constructed with known equal pressure and chemical potential and a prescribed density jump.

The following gates are covered:

1. Separate one-sided EOS providers and `P_t`, `mu_t`, `Delta n`, `Delta epsilon` continuity/jump checks.
2. The zero-temperature latent-energy identity and phase-local first-law residuals.
3. TOV restart exactly at `P_t`, with continuous `r,m,nu,lambda,nu'` and the expected high-density-core/low-density-envelope energy and density ordering.
4. Slow and rapid interface transfer conditions, each checking `[Delta P]=0` and its own displacement condition.
5. Independent Fortran-versus-Python transfer-map comparison for both conversion limits using the compiled native `case_driver`.
6. TOV convergence when tightening SciPy tolerances and the logarithmic-pressure step cap.
7. Smooth-interface limit: a nonzero but small density jump is used to compare phase-aware SciPy TOV and Cowling roots with the independent smooth-EOS native Fortran solver.

For the smooth-limit fixture, `Delta n/n_-` is approximately `5.0e-5`; the Seidov onset diagnostic is on the stable side. The phase-aware SciPy TOV result differs from the native Fortran smooth result by approximately `9.8e-6` in radius and `3.5e-5` in mass. At the same target mass, the native Fortran smooth solver returned `f=45896.5117 Hz` and `p1=65878.6305 Hz`. The phase-aware roots differ from those values by at most `4.5e-5` in relative frequency across the tested slow/rapid limits. The slow/rapid roots also converge toward each other as the density jump tends to zero. These are analytic-fixture cross-code checks, not realistic neutron-star predictions.

The stronger-jump fixture is used only to exercise a visibly discontinuous interface and both transfer maps. Its Seidov onset ratio is above the rapid-conversion threshold; it is not presented as a stable hybrid-star configuration.

## 6. Explicit non-claims and next gates

- The SciPy phase-aware TOV/Cowling path is an independent reference implementation, not the existing native Fortran production solver.
- No full stable-branch or `M_max` search is implemented for phase-transition stars.
- No Seidov-based branch filtering, twin-star sequence construction, or radial-mode stability analysis is implemented.
- No generalized mode tracking/overlap continuation is implemented; sign-change roots remain unlabelled outside the smooth-limit benchmark.
- Only one sharp transition and one globally conserved baryon charge are supported.
- No Gibbs mixed phase, multiple transitions, surface tension, finite conversion rate, composition perturbations, elasticity, rotation, magnetic fields, superfluidity, or metric perturbations are included.
- No GW damping times or full-GR polar-mode spectrum are computed.

The next implementation gate is a native Fortran phase-aware structure provider and segmented TOV/Cowling path that consumes the two one-sided EOS branches and the explicit interface event. That work must be cross-checked against the Python reference path before the Python implementation can be treated as a production replacement.
