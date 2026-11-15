# M7 CSS contract: one constant sound speed, crust-side first-order jump

Date: 2026-11-15  
Branch: remediation/2026-11-15-m5a-maxwell-eos-contracts

## Decision

The native match_crust_core(method="css") contract is changed from a two-segment, two-speed, continuous-density bridge to:

1. one high-density CSS phase with a single constant \(c_s^2\);
2. a first-order crust-to-CSS transition at \((P_t,\mu_t)\), with \(n_b>n_c\);
3. continuous matching of \(n_B,P,\epsilon,\mu_B\) from the CSS bridge to the selected core endpoint.

The other four smooth matching methods are unchanged. The separate public method="maxwell" path remains fail-closed because composition-aware native sharp-Maxwell integration is a different workstream.

The CSS output remains thermo-only: native matching still requires the caller to set drop_composition=True. The new transition fields are reported explicitly. No composition-preservation claim is made.

## Thermodynamic construction

Assumptions: \(T=0\), one conserved baryon charge, beta-equilibrated locally neutral barotropes, a constant and causal \(0<c_s^2\leq1\) bridge, and the Euler identity \(\epsilon=\mu_Bn_B-P\). Units in code are \(n_B\,[\mathrm{fm}^{-3}]\), \(P,\epsilon\,[\mathrm{MeV\,fm}^{-3}]\), \(\mu_B\,[\mathrm{MeV}]\).

Let the crust-side transition state be \((n_c,P_t,\mu_t)\), and the chosen core endpoint be \((n_2,P_2,\mu_2)\), with \(P_2>P_t\) and \(\mu_2>\mu_t\).

For a single CSS branch,
\[
n_B(\mu)=n_b\left(\frac{\mu}{\mu_t}\right)^{1/c_s^2},
\]
\[
P(\mu)=P_t+n_b\mu_t\frac{c_s^2}{1+c_s^2}
\left[
\left(\frac{\mu}{\mu_t}\right)^{1+1/c_s^2}-1
\right],
\qquad
\epsilon(\mu)=\mu n_B(\mu)-P(\mu).
\]

Enforce the core endpoint density:
\[
n_b=n_2\left(\frac{\mu_t}{\mu_2}\right)^{1/c_s^2}.
\]

The remaining pressure condition becomes a scalar root problem:
\[
P_2-P_t=n_2\mu_t\frac{c_s^2}{1+c_s^2}
\left[
\frac{\mu_2}{\mu_t}
-\left(\frac{\mu_t}{\mu_2}\right)^{1/c_s^2}
\right].
\]
The Fortran implementation solves this equation by bisection on \(c_s^2\in(0,1]\). If no causal root exists, or if the inferred \(n_b\leq n_c\), the CSS candidate is rejected. No endpoint or source knot is changed to force a root.

At the crust-to-CSS interface:
\[
P_{\rm crust}=P_{\rm CSS}=P_t,\qquad
\mu_{\rm crust}=\mu_{\rm CSS}=\mu_t,
\]
\[
\Delta n=n_b-n_c>0,\qquad
\Delta\epsilon=\mu_t\Delta n>0.
\]
This is a Maxwell-type first-order transition at the crust side.

At the CSS-to-core join, \(n_B,P,\mu_B,\epsilon\) are matched continuously by construction. The CSS branch has a single \(c_s^2\), but the source core endpoint's own derivative \(c_s^2\) need not equal it. The report therefore exposes both source endpoint sound speeds and the bridge sound speed; no C1 derivative-continuity claim is made.

## Seidov scope

For this sharp density jump, the small-core onset threshold is reported as
\[
\Delta\epsilon_{\rm crit}=\frac12\epsilon_c+\frac32P_t,
\qquad
\mathcal S=\frac{\Delta\epsilon}{\Delta\epsilon_{\rm crit}}.
\]
The Boolean seidov_small_core_stable records \(\Delta\epsilon<\Delta\epsilon_{\rm crit}\). This is only the Seidov small-core onset criterion. It is not a complete radial-mode or phase-interface Cowling stability calculation.

## Report and row layout

The output density table contains, in order:

- retained crust rows through the crust-side transition row;
- one adjacent CSS-side row at the same \(P_t,\mu_t\) but larger \(n_B,\epsilon\);
- \(n_{\rm bridge}\) interior CSS samples;
- the original core rows starting at the selected core endpoint, which is the continuous CSS/core join.

The native JSON report records:
css_single_sound_speed,
css_transition_nb_low,
css_transition_nb_high,
css_transition_delta_nb,
css_transition_delta_eps,
css_transition_epsilon_low,
css_transition_epsilon_high,
css_transition_pressure,
css_transition_mu,
css_transition_low_index,
css_transition_high_index,
css_core_bridge_index,
seidov_critical_delta_eps,
seidov_ratio, and seidov_small_core_stable.

Indices are zero-based in the returned output table. The first-order transition rows must be adjacent.

## Independently recomputed example

Using the frozen source rows:
- SPG_M1 crust row 364 (zero-based): eos.nb blob 671e4c2b8fe5e77539cb4854586c153f81be3f31, eos.thermo blob 2750bd1ba015b821061aa11cd18e5bca2d20405a;
- BBB_BHF2 core row 69 (zero-based): eos.nb blob 60725116001613f77384bf97e7cc4889af59145d, eos.thermo blob 806e2dffd0cfaf210cfb99604e72aa4f2731a310;
- crust neutron rest-mass gauge shifted by \(939.56533-938.9=0.66533\) MeV as the established representation-only transformation.

The endpoint values are:
\[
(n_c,P_t,\mu_t)=(0.07413102413,\;0.481250995766437,\;954.7895798154601),
\]
\[
(n_2,P_2,\mu_2)=(0.16,\;1.5166868679163006,\;963.9481015750811).
\]
The scalar solution is
\[
c_s^2=0.0128776977924919,\qquad
n_b=0.0762376403934158,
\]
\[
\Delta n=0.00210661626341577,\qquad
\Delta\epsilon=2.01137525697916\;\mathrm{MeV\,fm}^{-3}.
\]
This is an independent calculation from source knots, not yet an implementation-verified M7 output or an EOS physical-admissibility verdict.

## Literature context

The single-speed-plus-first-order-jump architecture is used in the generic hybrid-star CSS parameterization of Alford, Han & Prakash, where the transition pressure, energy-density discontinuity and high-density sound speed are independent control parameters: https://arxiv.org/abs/1302.4732.

Microscopic crust-to-uniform-core calculations also report weak first-order transitions with small density jumps, e.g. Douchin & Haensel (2001), https://arxiv.org/abs/astro-ph/0111092, and the crust-core discussion in Haensel et al. (2013), https://doi.org/10.1088/0004-637X/773/1/11.

This is not the only literature choice. Servignat et al. (2024) present bridge parameterizations designed to ensure continuity of internal energy, pressure and sound speed when stitching a crust to a high-density fit: https://doi.org/10.1103/PhysRevD.109.103022. That is a different interface contract, not evidence that every crust-core transition must be continuous.

## Verification still required

- Native Fortran build on this exact branch HEAD.
- Native match-core ... css ... integration test on the staged GPPVA_FSU2H/GM_GM1 fixture.
- Verify adjacent transition rows, exact jump identity, one bridge speed, and continuous core endpoint.
- CompOSE write/read-back comparison of \(n_B,P,\epsilon,\mu_B\) arrays and transition-row retention.
- Verify TOV/Cowling consumers do not silently treat the jump as a smooth interval. The existing smooth Cowling solver must not be assumed to implement phase-interface jump conditions.
