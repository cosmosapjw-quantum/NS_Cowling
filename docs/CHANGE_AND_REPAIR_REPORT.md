# v3 change and repair report

1. Replaced the default tabular EOS provider's cubic-Hermite (P(\mu_B)) representation with monotone PCHIP (g=\epsilon-P) versus (\ln P).
2. Kept the original pressure and energy-density knots exactly; reconstructed (n_B) by cold first-law quadrature and (\mu_B) by Euler closure.
3. Added independent source knot first-law admission at (2\times10^{-12}), separate from full-cell PCHIP monotonicity/causality checks.
4. Reparameterized cell-by-cell shooting to use log-pressure cells and the new EOS callback. TOV and Cowling differential equations and declared surface/central conditions were not modified.
5. Changed zero-amplitude `cell_solve` to pass the raw normalized table through untouched instead of recalculating energy density through a null perturbation.
6. Added Python `DirectBarotrope`, a CLI audit command, per-table receipts, and a standalone Fortran cohort audit script.
7. Retained v2 q-linear CSVs and the old 234-row mode cohort as explicitly labelled legacy comparison artifacts.
8. Added source-level thermodynamic rejection, causal PCHIP checks, pressure-relative inversion tolerance, and 16-point Gauss-Legendre first-law quadrature.

No original EOS CSV, CompOSE source member, or legacy v2 archive was intentionally overwritten. The method does not claim that all four source knot columns are simultaneously preserved by the continuous interpolant.
