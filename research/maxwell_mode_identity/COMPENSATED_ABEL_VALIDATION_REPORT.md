# Validation report: compensated Abel–Liouville continuation

Date: 2026-11-11  
Research id: COWLING-COMPENSATED-ABEL-INTERFACE-v1  
Scope: synthetic analytic coefficient matrices only.

## Commands and observed results

Run from a checkout whose repository root contains this directory:

    PYTHONPATH=research python -m pytest -q research/maxwell_mode_identity/tests/test_compensated_abel.py

Observed standalone result from the exported working directory:

    PYTHONPATH=/mnt/data/NS_Cowling_research_20261111 python -m pytest -q /mnt/data/NS_Cowling_research_20261111/compensated_abel/tests
    12 passed in 0.60s

    python -m compileall -q /mnt/data/NS_Cowling_research_20261111/compensated_abel
    exit status 0

    python -m json.tool /mnt/data/NS_Cowling_research_20261111/COMPENSATED_ABEL_RESULTS.json
    exit status 0

The 12 tests cover equal-closure reduction, sign and dimensions, exact constant-coefficient quadrature, invalid closure rejection, one-sided interface multiplier, singular-transfer rejection, multi-interface compensation, JSON-safe complex serialization, a two-segment matrix-exponential chain with nontrivial transfer and lapse jump, an independent DOP853 drift reproduction, interface telescoping, and preservation of an arbitrarily small imaginary component.

The TDD test for tiny complex phase first failed against a version that cast small imaginary values to float; it passed after preserving complex type. No tolerance was weakened to hide the failure.

## Numerical evidence

For a constant-coefficient interval, the trapezoid integral was:
- numerical: -0.009343853820598002
- analytic: -0.009343853820598004
- absolute error: 1.734723475976807e-18.

The two-segment matrix-exponential test uses a nontrivial 2x2 interface transfer with a lapse jump, one-sided enthalpies, distinct equilibrium/perturbation closure values, and complex fundamental columns. The compensated value returns to its initial value within relative tolerance 2e-12.

The DOP853 test independently reproduces both the uncorrected exponential drift and the corrected constant within relative tolerance 2e-9.

## Symbolic check

Wolfram Language returned exact zero for:
- bulk compensation residual;
- slow transfer determinant minus alpha;
- rapid transfer determinant minus alpha;
- standard interface multiplier residual for alpha=h_inner/h_outer and continuous lapse.

A prior malformed Wolfram request emitted a symbol warning; it is superseded by the successful exact-symbolic request above and is not counted as a successful check.

## Diagnostic integrity finding

A static audit of the live PR #6 source found that interface jump fields are emitted using the real projection only. For a complex jump -0.75-1.25i, the current expression emits -0.75 and loses an imaginary component of magnitude 1.25. This is recorded as an implementation finding; this PR does not patch the existing diagnostics function. The additive compatibility-preserving correction and regression test are assigned to local DAG node LA9.

## Provenance caveat

The 12-test run was executed in the current container against the exact contents of the new package authored in this session. It did not import or run the full GitHub checkout, and it did not authenticate a local Git worktree commit. GitHub branch identity and remote publication are recorded separately. Existing PR #6 has an older captured 42-case test receipt; its later 43-case prose claim is not treated as authenticated because the accompanying manifest/logs are inconsistent.

## Non-results

No local native build, q-to-r chain-rule audit, actual fixture profile, zero-jump star comparison, finite-jump continuation, TOV integration, or physical mode classification was run. The formula is conditional on the written 2x2 system; a physical nonbarotropic system with additional composition terms requires a new derivation.
