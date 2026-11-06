# v3 direct-provider research and adversarial audit

## Review sequence

1. Source discovery: verified compose_eos 3.3.0 source-knot definitions, normalized-table columns, existing EOS callback, TOV state access and Cowling solver call graph. The previous provider's cellwise P(mu_B) representation was the source of an unnecessary table transformation.
2. Step-back and alternative construction: compared free-energy Hermite, energy-density Hermite, monotone q-area constraints, and monotone g(ln P). Rejected simple F and epsilon Hermite variants after full-cohort cell checks; rejected exact monotone-q area matching where the endpoint integral feasibility condition failed.
3. Constructive derivation and implementation: derived the g(ln P) closure, first-law density quadrature and causal sound-speed identity; replaced the Fortran provider and refactored cell shooting to log-pressure coordinates. A null controlled perturbation now passes source data through unchanged.
4. Claim verification and adversarial gate: added independent source first-law admission before continuous interpolation admission. The deliberately inconsistent compose_134 fixture was rejected, while all 234 production EOS tables passed the source gate and all-cell PCHIP audit.
5. Numerical and deployment review: ran Python tests, manually compiled and linked the complete Fortran workspace, audited the full Fortran table cohort, recomputed one representative v3 f/p1 spectrum, checked node identities, and compared with the v2 q-linear result.

No native parallel subagent facility was invoked. No simulated delegation or synthetic reviewer result is claimed.

## Principal scientific decision

The source generator already produces first-law-consistent knots. v3 does not repair its EOS and does not create another pressure-energy table. It preserves source P and epsilon knots, imposes a monotone causal continuous g(ln P) interpolant, and reconstructs n_B and mu_B from the cold first law. The density knot difference is reported, not hidden.

## Verification summary

- Python pipeline: 8 tests passed.
- Fortran: 62 modules compiled and both main and cell-shooting executables linked with gfortran Fortran 2018 and runtime checking.
- Fortran source table admission: 234/234; 109,776 cells; zero rejected cells; maximum raw first-law closure 6.1990e-16.
- Negative fixture: raw closure 1.1327e-4; rejected.
- Representative direct-provider spectrum: M residual 1.54e-13; f/p1 matching residuals 4.16e-9 and 5.61e-12; expected node counts preserved.
- FPM 0.13.0 was unavailable; official FPM harness verification is not claimed.
