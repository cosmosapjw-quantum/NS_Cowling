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
- Historical note: no GitHub Actions workflow had run when this initial bounded validation note was first prepared. The later exact-PR-head run is recorded below.

The phase-aware sharp-interface TOV/Cowling solver is Python/SciPy reference code. The Fortran code independently implements and tests the slow/rapid interface-transfer algebra, but native Fortran solve_tov_star and solve_cowling_spectrum remain smooth-EOS-only. No realistic CompOSE hybrid-EOS cohort, stable branch, maximum mass, twin-star sequence, full-GR polar perturbation, gravitational-wave damping, or generalized mode-family continuation is claimed.

## Provenance

- Smooth solver baseline ZIP SHA-256: 6d2ba545101df4cdc5bb5aed187b640ec84e528bf2005db764792bc4e19d9885.
- Maxwell-enabled compact source ZIP SHA-256: e7095d95b58022f95c1ee8bbc0692fb746eab1c676737388e58dae8ef82366b1.
- Patch SHA-256: de9b51aacbf8d022ffe5dbd0787b7f38d82f6883597f4b15963381530a0f05c5.

These are local package and test artifacts. Their presence in this branch as a patch and documentation is not equivalent to integrating the patch into the native production solver.


## Final GitHub exact-head verification (2026-11-09)

The source patch is now applied to `research/maxwell-reference-validation-v1` at commit `6b23d024e41661d33b83d6202038db6da918821d`. GitHub Actions [run 42579257980](https://github.com/cosmosapjw-quantum/NS_Cowling/actions/runs/42579257980) completed successfully on the pre-integration head and produced the source commit only after all gates passed:

- Patch applicability and application: PASS.
- Native Fortran EOS build: PASS.
- Focused Maxwell interface tests: **7 passed**.
- Independent Maxwell reference script: `PASS_BOUNDED_REFERENCE_VALIDATION`.
- Full `python_frontend/tests pipeline/tests`: **74 passed**.
- Evidence artifact ID: `13657075483`; ZIP SHA-256: `93347a170c238458fa937f2359e9ae4a8737098c5f9ec9dbdb09a26351f9946b`.

This does not change the PR #2 implementation boundary: its sharp-interface phase-segmented TOV/Cowling path is Python/SciPy, while its Fortran contribution provides the slow/rapid interface-transfer algebra and direct tests. The separate native phase-segmented TOV/Cowling integration is tracked on research/maxwell-native-fortran-v1. The analytic strong-jump fixture remains a boundary-condition test and is not a stable hybrid-star prediction.

A follow-up exact-PR-head Actions run, [42579770412](https://github.com/cosmosapjw-quantum/NS_Cowling/actions/runs/42579770412), checked PR #2 head 2b7a00984dfa215bcca48c7d7adb5f3619013502 and completed with conclusion success. Its exact-head artifact is 13656363884 (SHA-256 0cfba640210b84d3141213b820f5c335d67ea1abef4a9e39b204d4a41b82e4f6).


## Native Fortran follow-up (2026-11-09)

A separate native phase-segmented Fortran TOV/Cowling path has since been added on `research/maxwell-native-fortran-v1` and is under draft [PR #3](https://github.com/cosmosapjw-quantum/NS_Cowling/pull/3). Exact head `98d2f9ce6033015109286e5713639c6626a6f7fb` passed both push and PR GitHub Actions runs `42581655510` and `42581660400`; its analytic fixture passed 18 checks and the production case driver built. This is a new sibling implementation, not a modification of `solve_tov_star` or `solve_cowling_spectrum` in this branch. It still lacks an executed realistic two-table EOS case and a complete mode spectrum. The earlier PR #2 implementation boundary remains unchanged.