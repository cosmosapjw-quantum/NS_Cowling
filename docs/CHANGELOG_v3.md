# Changelog v3

## Scientific and numerical changes

- Default table EOS provider changed from cubic-Hermite (P(\mu_B)) to monotone PCHIP (g=\epsilon-P) versus (\ln P).
- Original (P,\epsilon) knots are retained; (n_B,\mu_B) follow from the cold first law and Euler identity.
- Source knot first-law gate, positive/monotone pressure-energy checks, analytic cell derivative/causality checks, and no-extrapolation behavior added.
- First-law density integration uses 16-point Gauss-Legendre quadrature per native log-pressure cell.
- The legacy cell-shooting path now queries the direct EOS provider in log-pressure coordinates.
- A zero-amplitude controlled perturbation uses the source table as a direct passthrough.
- v2 q-linear tables/results remain available as legacy comparison artifacts.

## Validation and evidence

- 8 Python pipeline tests passed.
- 62 Fortran modules and both executables compiled/linked manually with gfortran Fortran 2018; FPM 0.13.0 was unavailable.
- All 234 source EOS tables admitted by the direct Fortran provider; 109,776 intervals, none rejected.
- An intentionally inconsistent table was rejected by the raw first-law gate.
- One representative direct-provider Cowling f/p1 spectrum recalculated and matched against its v2 q-linear counterpart.
