# Scientific contracts for the remediation

## 1. Cold one-charge barotrope

Assume (T=0), a single conserved baryon charge, a differentiable homogeneous phase and units with (c=1) for the sound-speed expression. For a thermodynamically consistent one-charge branch,

[
n_B(\mu_B)=\frac{dP}{d\mu_B},\qquad
\epsilon=\mu_B n_B-P.
]

Differentiating the Euler identity and using (P'=n_B),

[
\frac{d\epsilon}{d\mu_B}
=n_B+\mu_B\frac{dn_B}{d\mu_B}-\frac{dP}{d\mu_B}
=\mu_B\frac{dn_B}{d\mu_B}.
]

Consequently,

[
c_s^2=\frac{dP}{d\epsilon}
=\frac{n_B}{\mu_B\,dn_B/d\mu_B}.
]

For positive (n_B,\mu_B), mechanical stability requires (dn_B/d\mu_B>0), and causality requires

[
0<c_s^2\leq1
\quad\Longleftrightarrow\quad
\mu_B\frac{dn_B}{d\mu_B}\geq n_B.
]

Units in the repository: (n_B) in fm(^{-3}), (P,\epsilon) in MeV fm(^{-3}), and (\mu_B) in MeV. Thus (dP/d\mu_B) has the correct units fm(^{-3}), and (c_s^2) is dimensionless.

For monotone (n_B(\mu_B)), an interval also obeys the necessary secant bound

[
n_i\leq
\frac{P_{i+1}-P_i}{\mu_{i+1}-\mu_i}
\leq n_{i+1}.
]

A violation already present in the raw knots cannot be repaired by an interpolant that preserves those same endpoint values and slopes. It is a source-table consistency issue, not merely a spline overshoot.

## 2. Maxwell interface

For two locally neutral phases at (T=0), a sharp Maxwell transition requires

[
P_L(\mu_t)=P_H(\mu_t)=P_t,
\qquad \mu_L=\mu_H=\mu_t,
]

with (n_H>n_L>0) and

[
\Delta\epsilon=\epsilon_H-\epsilon_L
=\mu_t(n_H-n_L)>0.
]

The endpoint jump identity does not establish the differential first law inside either phase. Check both separately.

The Seidov condition is a small-core onset criterion; it is not a substitute for a full radial-mode eigenvalue calculation.

## 3. CompOSE nuclear-group composition

The CompOSE Reference Manual defines particle number fractions as

[
Y_i=\frac{n_i}{n_B}.
]

For a nuclear group (mathcal M_I),

[
A_I^{\rm av}
=\frac{\sum_{j\in\mathcal M_I}A_jY_j}
       {\sum_{j\in\mathcal M_I}Y_j},
\qquad
Y_I=\sum_{j\in\mathcal M_I}Y_j.
]

Therefore the group baryon fraction is

[
A_I^{\rm av}Y_I
=\sum_{j\in\mathcal M_I}A_jY_j.
]

The group fraction (Y_I) itself is a number fraction, not the baryonic mass fraction. Provider-specific nonstandard encodings must be named and independently reconciled; do not silently generalize one provider's convention to all CompOSE tables.

## 4. Evidence semantics

Keep the following claims distinct:

- **Native unit test:** a Fortran routine was called with a test fixture.
- **Native CLI integration:** the compiled Fortran executable was invoked through its public command.
- **Native raw-source end-to-end:** actual CompOSE source files were consumed by the native command and the output was checked.
- **Python mirror:** Python independently recomputed the intended equations; this is not native Fortran execution.
- **Read-back row-count match:** same number of rows; insufficient by itself.
- **Read-back numerical identity:** all relevant physical arrays and metadata satisfy a declared residual budget.
- **Composition preservation:** all source composition values, phase labels, presence masks and applicability masks survive the transformation and native round-trip.

A status is PASS only at the evidence level actually exercised.

## 5. M7 CSS: one speed with a crust-side Maxwell jump

The current native M7 `css` method uses one constant (c_s^2\in(0,1]) over the bridge. It has a sharp first-order crust-to-bridge transition with common (P_t,\mu_t), positive (\Delta n) and (\Delta\epsilon=\mu_t\Delta n), then a continuous state-variable match (n_B,P,\epsilon,\mu_B) to the selected core endpoint. The exact scalar matching equations and representative source-row calculation are in `docs/CSS_SINGLE_SPEED_CRUST_JUMP_CONTRACT_20261115.md`.

The core-side state is continuous, but the core source derivative (c_s^2) need not equal the constant bridge value; therefore do not claim C1. Seidov is reported only as a small-core onset criterion, not a complete radial-mode or Cowling phase-interface stability result. Composition is explicitly dropped; C3/C4 are NOT_AVAILABLE, not PASS.

The 2026-11-14 two-segment smooth CSS artifact and its C⁰/C¹ validation are historical evidence for the previous method. Preserve that artifact unchanged and do not treat it as validation of the current single-speed/jump contract.

## 6. References

1. CompOSE Reference Manual, Version 3.01, especially particle fractions and nuclear-group composition in Eqs. (3.10)–(3.12), (4.20)–(4.22), and thermodynamic consistency in §4.1.2: https://compose.obspm.fr/download/pdf/manual_v3.00.pdf
2. Lindblom, “Spectral representations of neutron-star equations of state,” *Physical Review D* (2010), SciSpace record: https://scispace.com/papers/spectral-representations-of-neutron-star-equations-of-state-2ap7rvsxz2
3. Zero-temperature EOS stability/causality/thermodynamic-consistency discussion in arXiv:2208.03086: https://arxiv.org/abs/2208.03086

The algebraic sound-speed identity in §1 was also symbolically checked in Wolfram Language:
(d\epsilon/d\mu_B=\mu_B n'_B) and (c_s^2=n_B/(\mu_B n'_B)), conditional on (P'=n_B).
