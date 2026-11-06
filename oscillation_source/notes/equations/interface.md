# Material-interface Cowling contract ()

## Scope

This contract is for a declared, surface-layer-free **slow/material interface** at a finite radius
`r_t > 0` and a strictly positive trial angular frequency `omega > 0`.  It is not a finite-width
mixed phase, rapid conversion law, elastic shell, or zero-frequency limit.

The equilibrium transition has a common pressure `P_t` and one-sided total energy densities
`epsilon_in` and `epsilon_out`.  The background metric variables and `nu'` are evaluated at the
same areal radius.  The perturbation variables use the  convention `[W]=L^3`, `[V]=L^2`.

## Transfer derived from the physical conditions

The adopted material-interface conditions are

\[
[W]=0,\qquad [\Delta P]=0.
\]

With

\[
a=\left(\frac{\omega r_t}{c}\right)^2 e^{\lambda_t-2\nu_t},\qquad
Q=aV+\nu'_t W,
\]

and

\[
\Delta P=-\frac{(\epsilon+P_t)e^{-\lambda_t}}{r_t^2}Q,
\]

continuity of `Delta P` gives

\[
(\epsilon_{\rm in}+P_t)(aV_{\rm in}+\nu'_tW)
=(\epsilon_{\rm out}+P_t)(aV_{\rm out}+\nu'_tW).
\]

Therefore

\[
W_{\rm out}=W_{\rm in},\qquad
V_{\rm out}=qV_{\rm in}+(q-1)\frac{\nu'_tW}{a},\qquad
q=\frac{\epsilon_{\rm in}+P_t}{\epsilon_{\rm out}+P_t}.
\]

The code derives and tests this enthalpy-density ratio.  It does **not** enforce `V_out=V_in`
when the energy density jumps.  It also does not divide by `a` at `omega=0`; that case is outside
this  contract.

### Source comparison caution

Ranea-Sandoval et al., *Oscillation modes of hybrid stars within the relativistic Cowling
approximation*, arXiv:1807.02166, discuss the sharp-interface Cowling problem.  The printed PDF
form of their Eq. (3.6) contains a denominator with `epsilon_+ - P`.  Applying their stated
continuity of the Lagrangian pressure perturbation to the same variable definitions instead gives
the enthalpy factor `epsilon_+ + P` used above.   therefore treats the interface conditions and
the direct derivation above as the 기준; it does not claim that the authors' implementation
used the printed sign.

## Event segmentation

Material jumps are declared by transition pressure.  Since / integrate outward in

\[
q_P=\ln(P_{\rm start}/P),
\]

a transition pressure maps to one exact stop coordinate.  Events must be strictly ordered in
outward direction.  A material event changes the allowed EOS side and applies the transfer above.
A `chi`-only knot changes the derivative law but leaves `(W,V)` continuous.

The integration contract is:

1. integrate only to the next declared event coordinate;
2. preserve the one-sided state on both sides at the same `r_t`;
3. apply the transfer at a material interface;
4. restart from the outer state;
5. never let a single RK trial stage bridge an undeclared discontinuity.

The  `segment_restart_probe` is deliberately a minimal infrastructure test using a smooth zero
RHS between events.  It proves endpoint landing, restart, side preservation, and transfer wiring.
It is not a calibrated discontinuous-star spectrum.

## Restricted mode identity

A sharp density jump may introduce an interface/gravity-mode candidate.  Consequently a root is
not labelled `f` or `p1` merely because of frequency ordering.   accepts an `f`/`p1` family label
only when a root and transfer are valid, a smooth-reference family is supplied, and a declared
overlap diagnostic passes its threshold.  A candidate interface-g mode remains `unresolved`.

The overlap passed to this helper is a diagnostic produced elsewhere;  does not create a new
physical mode-energy inner product.
