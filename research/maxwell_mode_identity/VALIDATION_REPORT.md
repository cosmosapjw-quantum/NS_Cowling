# Validation report: Cowling mode identity, weighted concomitant and two-sided matching

Date: 2026-11-11  
Command: PYTHONPATH=<package-root> python -m pytest -q mode_identity/tests

Observed fresh local result: 43 passed in 0.74s (33 test functions; 43 collected cases after parameter expansion). The zero-inner-concomitant JSON-safety test was first observed failing against NaN output, then passed after the minimal null/status correction (red-green).

## Assertions tested

- Slow/rapid transfer maps preserve pressure/conversion conditions; determinant equals h_inner/h_outer and zero-jump maps are identity.
- Profile norm/trapping, interface polarization, phase-local node handling, profile overlap and ambiguity refusal.
- Smooth-phase weighted determinant constancy for independent complex columns in one synthetic cold-barotropic path.
- Both interface maps preserve the h-weighted determinant; multi-interface factors telescope.
- Complex determinant is bilinear without complex conjugation.
- Exact matrix-exponential trace tests for ell=1,2,3,5.
- Nonbarotropic chi_ad/chi_eq mismatch produces the predicted drift.
- An artificial lapse jump is detected.
- Surface boundary null vector reproduces the current normalized ell=2 surface residual zero set.
- Independent center-outward and surface-inward integrations produce a weighted determinant spread of 4.673260016455182e-12 at five match radii with DOP853 max_step=100 m.
- Negative control: the same synthetic two-sided case without a max-step cap gives 3.576076639831998e-4 spread. This is retained as a numerical-conditioning finding.
- Zero inner concomitant: the ratio is explicitly undefined and now serializes as null plus UNDEFINED_ZERO_DENOMINATOR; strict JSON with allow_nan=false passes for those diagnostic fields.

Separate smooth-phase run: maximum relative drift 7.907195134185742e-11; endpoint drift 4.723167073919359e-11; 602 RHS evaluations.

Additional checks: Python compileall PASS; profile_schema.v2.json, LOCAL_AGENT_DAG.v1.json and MODE_IDENTITY_STATUS.json parse.

## Computation provenance

The bounded run receipt is computation_evidence/weighted_wronskian_20261111/computation_manifest_v2.json. It records Python 3.13.5, NumPy 2.3.5, SciPy 1.17.0, pytest 9.0.2, commands, inputs, outputs, seed and observed numerical values. The test bundle ran from an exported artifact directory, not a Git worktree; the manifest does not invent a commit SHA and explicitly records that source-revision authentication was unavailable in that run.

## Verification boundary

All numerical assertions use synthetic cold-barotropic coefficient functions and exact algebraic interface maps. The two-sided test is an independent ODE integration of synthetic center/surface basis columns. No production Fortran was changed or rebuilt. Actual four native fixture-root profiles have not been processed. The native q-to-r chain rule still requires local re-audit. No source-backed CompOSE EOS, Maxwell pair, TOV sequence or physical mode classification is claimed.

## Known risks

1. Outward/inward solutions may acquire very different scales. max_step refinement and QR/log-amplitude determinant accounting are mandatory.
2. Both directions must use the same finite-pressure surface cutoff.
3. A nonbarotropic or frozen-composition chi requires a new derivation.
4. The concomitant is not established as a canonical covariant symplectic current.
5. Physical labels remain unresolved pending zero-jump seeding and actual fixture continuation.



## 2026-11-11 supplemental compensated Abel–Liouville check

The separate 12-test synthetic helper suite passes. See COMPENSATED_ABEL_VALIDATION_REPORT.md and COMPENSATED_ABEL_RESULTS.json for exact commands and tested range. This does not authenticate or replace the historical 42-case PR #6 test receipt and does not constitute a full repository replay.

The 43-case claim in the prior prose/status files does not have a matching 43-case manifest and captured stdout in the inspected branch. It is therefore not counted as verified. The current PR #6 combined-status endpoint returned no status records; GitHub Actions/CI pass is not claimed.
