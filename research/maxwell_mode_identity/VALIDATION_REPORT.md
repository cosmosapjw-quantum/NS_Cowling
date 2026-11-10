# Validation report

Date: 2026-11-11. Local package command: `python -m pytest -q mode_identity/tests`. Result: **22 passed**.

The tests cover:
- slow/rapid transfer pressure continuity and their distinct kinematic conditions;
- determinant alpha and the alpha=1 identity limit;
- rejection of invalid coefficients, states and conversion labels;
- positive norm, quadratic amplitude scaling and invariant trapping fractions;
- one-sided interface polarization with the physical r_t factor;
- symmetric, amplitude/global-sign invariant profile overlaps;
- retention of relative core-envelope phase/sign information in the global overlap;
- refusing a real sign-reversal label for genuinely complex profiles;
- global assignment of spectral branches and ambiguity refusal;
- prohibition on mixing slow and rapid profiles in one continuation call.

**Verification boundary:** these are deterministic synthetic helper tests and symbolic interface-map checks. No production Fortran source was changed or rebuilt. The helper was not run on the four accepted native fixture-root profiles because the exact emitted profiles are not present in the local working directory. Therefore no physical mode-family classification or new stellar spectrum is claimed.

Known limitations: the norm is not a proof of cross-background canonical orthogonality; the overlap is a chosen continuation similarity; node counting continues to use the existing project Python implementation and requires its own mesh/tolerance convergence; the assignment is not predictor-corrector continuation and does not transport a degenerate subspace; the outer profile endpoint follows the existing finite-pressure cutoff and is not automatically the exact zero-pressure surface.
