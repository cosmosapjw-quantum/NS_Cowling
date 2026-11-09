# Adversarial audit: sharp Maxwell TOV/Cowling reference path (2026-11-09)

## Decision

**Bounded reference-path validation: PASS. Native sharp-interface production solver: NOT IMPLEMENTED. Stable-branch/astrophysical validation: NOT CLAIMED.**

The existing native EOS Maxwell constructor was connected to a two-phase TOV and Cowling integration path in Python/SciPy. The implementation explicitly restarts at the exact pressure transition and never interpolates the density jump as a smooth table. Independent Fortran slow/rapid transfer routines are compiled and compared against the Python transfer map.

## Findings and repairs

1. **Wrong side convention in an initial assertion.** The inner core is the high-density phase; the outer envelope is the low-density phase. A test initially asserted the reverse ordering. The data and test convention were corrected; the reported positive jump is `inner - outer`, consistent with high-density minus low-density thermodynamics.
2. **Pressure-event stage leakage.** A radial integration could evaluate a stage infinitesimally outside a phase's pressure domain. The TOV solver now integrates in `q=ln(P_ref/P)` and restarts at exact `q_t` and `q_surface` endpoints. Query values are clipped only for floating-point endpoint roundoff, not for physical extrapolation.
3. **Unused step controls.** A radial `max_step_m` argument was not a meaningful cap in the logarithmic-pressure independent variable. It was replaced with the explicit `max_step_log_pressure`, which is now passed to both phase integrations. The radius limit is an event, and hitting it before `P_t` or the surface raises an error.
4. **Rapid/slow conversion distinction.** The old Fortran material transfer encoded `[xi^r]=[Delta P]=0`, i.e. the slow limit. A separate `rapid_interface_transfer` was added using `[Delta P]=0` and `[xi^r-Delta P/P'_0]=0`. The Python and Fortran maps are compared at the same interface state.
5. **Numerical root scope.** Mode scans return unlabelled real-frequency roots and sampled node counts. They do not infer `f`, `p1`, or interface-mode identities generically. The smooth-limit test identifies the first two roots against the established smooth Fortran spectrum only for that analytic fixture.

## Verification evidence

Focused command:

```bash
PYTHONPATH=python_frontend python -m pytest -q python_frontend/tests/test_maxwell_star.py
```

Verified result: **7 passed**. This includes the Fortran/Python transfer comparison and the smooth-limit TOV plus Cowling frequency benchmark.

The smooth-limit fixture uses a nonzero density jump with `Delta n/n_- ~= 5e-5` and a Seidov onset ratio below the rapid-conversion threshold. Independent Fortran smooth-EOS reference values at the same target mass are:

- TOV: `M=3.0839736617e27 kg`, `R=903.1705824 m`.
- Cowling: `f=45896.5117446 Hz`, `p1=65878.6304997 Hz`.

The phase-aware Python result differs from the smooth Fortran TOV values by `Delta R/R = -9.75084e-6` and `Delta M/M = -3.53194e-5`. The slow-conversion phase-aware root frequencies differ from the smooth Fortran `f/p1` values by `+9.40e-7` and `-1.24e-6`; the rapid-conversion differences are `-3.74808e-5` and `-2.50387e-5`. Slow-versus-rapid frequency differences are `3.84219e-5` and `2.37973e-5` for the two cross-matched roots.

The independent Fortran and Python transfer maps agree with zero relative difference in `W_out` and a relative difference of `2.13e-16` in `V_out` for both conversion limits; the Fortran `Delta P` jump diagnostic is zero in both. The focused test also checks phase-local first-law residuals around `2.2e-16`, the latent-energy identity (relative residual `1.26e-10` for the small-jump fixture), exact interface restart, continuity of the metric/TOV state, and tightening-tolerance convergence. The reproduced machine-readable values are in `validation/MAXWELL_REFERENCE_VALIDATION_20261109.json`; rerun with `PYTHONPATH=python_frontend python scripts/validate_maxwell_interface.py`.

## Known unverified areas

- The compact-package regression command completed with **76 passed, 8 deselected**. The eight exclusions are the documented corpus/cohort integrity and legacy differential tests that require external/full-corpus assets.
- No GitHub CI has run on this new code.
- Native Fortran `solve_tov_star` and `solve_cowling_spectrum` remain smooth-EOS-only; the sharp-interface stellar integrator is Python/SciPy reference code.
- No realistic CompOSE hybrid EOS cohort, stable branch, maximum mass, Seidov branch filtering, or twin-star sequence has been calculated.
- No generalized mode-family continuation is implemented. Interface modes and possible root crossings must be treated as unresolved.
- No full-GR interface perturbation calculation or GW damping comparison is included.

## Provenance

The solver source is the compact v5 package, pinned before these edits by `NS_Cowling_Fortran_Python_v5_20261109.zip` SHA-256 `6d2ba545101df4cdc5bb5aed187b640ec84e528bf2005db764792bc4e19d9885`. The original solver ZIP and the earlier 234-EOS canonical results were not modified. This audit describes local source/test evidence only; it is not evidence of a GitHub push, CI run, PR merge, or a publication result.
