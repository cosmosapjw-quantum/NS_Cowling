# Changelog v4

- Added a standalone Fortran CompOSE 1-D cold table reader, EOSData serializer, CompOSE writer and C1-C11 consistency verifier.
- Added a Python package frontend which calls the Fortran EOS executable; no silent fallback to the historic Python numerical implementation.
- Added Python wrappers for the existing Fortran TOV and Cowling solvers, including an automatic gfortran source build cache.
- Added a feature-parity matrix that explicitly identifies unported legacy generator/matching/literature/phase-transition capabilities.
- Kept the previous Python compose_eos package as an offline differential reference.
- Kept the v3 direct P-epsilon representation as the default stellar EOS provider.


## v4 follow-up: native master-function density-gap fill

- Added a Fortran PCHIP/First-Law density-gap reconstruction kernel with an integral correction and a causal cubic-Hermite fallback.
- Exposed `fill_density_gap` through the Python frontend and native CLI; cross-phase gaps require an explicit `gap_phase`.
- Added synthetic EOS tests, a real GPPVA source-table test, and a differential comparison against the historic Python implementation.

- Hardened the cubic-Hermite fallback gate: reset candidate admissibility independently of the failed PCHIP candidate and evaluate analytic F'' so dP/dn_B=F'+F'' is checked rather than substituting a duplicate F' term.
- Added regression coverage for both acceptance of a physical Hermite fallback and rejection of a positive-Q1 fallback with negative dP/dn_B.
- Fixed Hermite fallback sampling after early PCHIP rejection: fallback abscissae are now recomputed from the local gap coordinate rather than reusing the potentially partially initialized dense-validation array.


## v4 follow-up: native smooth crust-core matching

- Added six strict-thermo-only Fortran constructions: near-pressure plateau compatibility bridge, two-segment CSS, cubic/quintic pressure--chemical-potential bridges, and cubic/quintic free-energy bridges.
- Added a Python `match_crust_core` wrapper preserving the legacy `(EOSData, info)` return contract and returning the native construction report plus the independent C1--C11 consistency record.
- Added explicit state gates for cold T=0, beta-equilibrium metadata, zero entropy, Q6=Q7, endpoint sound speed, positive monotonic thermodynamics, subluminal bridge sound speed and First-Law closure.
- Composition-bearing inputs now require explicit `drop_composition=True`; composition-aware reconstruction and `literature_beta_eq` remain unported.
- Added differential regression tests against the historic Python package and a self-contained implementation/audit note. Final differential CI status is recorded against the exact branch HEAD.


## v4 follow-up: native one-charge Maxwell transition

- Added a separate Fortran Maxwell construction for two supplied cold, beta-equilibrated, single-phase barotropes, with local charge neutrality explicitly confirmed by the caller.
- Enforced (dP/d\mu_B=n_B) and \(\epsilon=\mu_B n_B-P\) through a cubic-Hermite (P(\mu_B)) representation; mechanical stability and causality are checked over each full interpolant interval.
- Isolated all pressure-equality roots of the piecewise-cubic pressure difference by its stationary points; the v1 policy accepts exactly one isolated coexistence point.
- Wrote two equal-pressure/equal-chemical-potential transition rows with a finite density and energy-density jump, and checked \(\Delta\epsilon=\mu_t\Delta n_B\).
- Masked derivative stencils at every phase-label boundary and blocked complete CompOSE export when phase composition is absent.
- Gibbs, composition reconstruction, Seidov threshold evaluation and phase-interface TOV stability remain unported or unverified; see the self-contained Maxwell note and audit record.
