# Sharp Maxwell TOV/Cowling reference validation v1

Date: 2026-11-09

## Decision and scope

**Bounded independent reference validation: PASS. Native Fortran sharp-interface production solver: NOT IMPLEMENTED. Stable-branch or astrophysical population claim: NOT MADE.**

This increment separates the true first-order Maxwell transition from the historical smooth maxwell crust-core matching compatibility token. The existing native EOS constructor supplies a phase-labelled two-phase barotrope with a finite density and energy-density jump. The independent stellar reference path is implemented in Python/SciPy and consumes the phase-labelled representation without smoothing across the interface.

The patch is preserved at research_artifacts/maxwell/NS_Cowling_V5_MAXWELL_V1.patch. Its base is the compact v5 source snapshot, not the full GitHub branch; see the handoff contract in GITHUB_MAXWELL_HANDOFF_20261109.md.

## Thermodynamic interface

For a single globally conserved baryon charge, with local electric neutrality imposed independently in each supplied phase, the coexistence point obeys

\[
P_-(\mu_t)=P_+(\mu_t)=P_t,
\qquad
\mu_{B,-}=\mu_{B,+}=\mu_t,
\]

\[
n_B=\frac{dP}{d\mu_B},
\qquad
\epsilon=\mu_B n_B-P.
\]

Consequently,

\[
\Delta\epsilon=\mu_t\Delta n_B,
\]

with the jump defined high-density minus low-density. Both phase providers retain their own thermodynamic interpolants. No interpolation stencil crosses the discontinuity.

## Stellar integration and phase conversion

The phase-aware TOV integration uses logarithmic pressure as the independent variable and terminates the inner high-density integration at the exact transition pressure. It restarts the outer low-density integration at that same interface state. The interface is therefore an integration event, not a sampled EOS knot to be crossed by interpolation.

For the Cowling perturbation, slow and rapid conversion are distinct boundary-value problems:

\[
\begin{aligned}
\text{slow:}\quad &[ \xi^r ]=0,\qquad [\Delta P]=0,\\
\text{rapid:}\quad &[ \Delta P]=0,\qquad
\left[\xi^r-\frac{\Delta P}{P'_0}\right]=0.
\end{aligned}
\]

The independent Fortran/Python transfer comparison covers both maps. The mode-frequency scanner returns unlabelled roots in general; the small-jump smooth-limit fixture is the only case where the first two roots are cross-matched against the already identified smooth f and p1 roots. General mode-family continuation and interface-mode identification remain unresolved.

## Validation evidence

The machine-readable result embedded in the Maxwell-enabled compact package is validation/MAXWELL_REFERENCE_VALIDATION_20261109.json, schema maxwell-interface-validation-v1, status PASS_BOUNDED_REFERENCE_VALIDATION.

### Independent interface-transfer comparison

| Conversion | W relative difference | V relative difference | Python interface residual | Fortran/Python Delta-P jump |
|---|---:|---:|---:|---:|
| Slow | 0 | 2.12867e-16 | 0 | 0 |
| Rapid | 0 | 2.12874e-16 | 4.16056e-16 | 0 |

### Small nonzero-jump smooth-limit benchmark

The analytic fixture has \(\Delta n/n_-=4.99988\times10^{-5}\), below the rapid-onset Seidov diagnostic threshold (ratio 5.00520e-5 versus critical ratio 0.501596). It approaches the smooth EOS limit while retaining a nonzero discontinuity.

| Quantity | Phase-aware reference | Smooth native Fortran | Relative difference |
|---|---:|---:|---:|
| TOV mass (kg) | 3.08386473749e27 | 3.08397366173e27 | -3.53194e-5 |
| TOV radius (m) | 903.161775722 | 903.170582394 | -9.75084e-6 |
| Slow f-like root (Hz) | 45896.5548689 | 45896.5117450 | +9.39590e-7 |
| Slow p1-like root (Hz) | 65878.5486772 | 65878.6305003 | -1.24203e-6 |
| Rapid first root (Hz) | 45894.7915057 | 45896.5117450 | -3.74808e-5 |
| Rapid second root (Hz) | 65876.9809831 | 65878.6305003 | -2.50387e-5 |

The TOV tolerance refinement gives relative changes 1.09383e-9 in mass and 1.69247e-10 in radius. The phase-local first-law residuals are 2.20e-16 and 2.36e-16. The Maxwell identity residual is 1.25666e-10.

### Strong-jump boundary-condition fixture

A separate unit-scale fixture uses n_low=1, n_high=2 fm^-3, mu_t approximately 941 MeV and P_t approximately 1 MeV fm^-3. The Maxwell identity residual is 2.66880e-13. At a representative 1800 Hz interface state, the rapid condition residual is 2.46e-15 and the slow condition residual is zero under their respective maps.

This fixture has Seidov ratio 1.00106, above its critical value 0.501596. It is therefore used only to exercise a large discontinuity and distinguish the interface maps; it is **not** admitted as a stable hybrid-star prediction. The TOV integrator completes the reference integration, but that does not establish radial or nonradial stability.

## Test record and limits

- Focused Maxwell test: python_frontend/tests/test_maxwell_star.py, **7 passed**.
- Independent validation script: python scripts/validate_maxwell_interface.py, status PASS_BOUNDED_REFERENCE_VALIDATION.
- Python compile check: python -m compileall -q python_frontend/compose_cowling python_frontend/compose_eos scripts/validate_maxwell_interface.py, passed.
- Earlier compact regression run: 76 passed, 8 documented corpus/legacy differential tests deselected.
- A later broader rerun exceeded its 120 s execution bound after 19 tests; this rerun is recorded as **timed out**, not passed.
- No GitHub Actions workflow has run on this increment yet.

The phase-aware sharp-interface TOV/Cowling solver is Python/SciPy reference code. The Fortran code independently implements and tests the slow/rapid interface-transfer algebra, but native Fortran solve_tov_star and solve_cowling_spectrum remain smooth-EOS-only. No realistic CompOSE hybrid-EOS cohort, stable branch, maximum mass, twin-star sequence, full-GR polar perturbation, gravitational-wave damping, or generalized mode-family continuation is claimed.

## Provenance

- Smooth solver baseline ZIP SHA-256: 6d2ba545101df4cdc5bb5aed187b640ec84e528bf2005db764792bc4e19d9885.
- Maxwell-enabled compact source ZIP SHA-256: e7095d95b58022f95c1ee8bbc0692fb746eab1c676737388e58dae8ef82366b1.
- Patch SHA-256: de9b51aacbf8d022ffe5dbd0787b7f38d82f6883597f4b15963381530a0f05c5.

These are local package and test artifacts. Their presence in this branch as a patch and documentation is not equivalent to integrating the patch into the native production solver.
