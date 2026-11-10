# Sharp-Maxwell Cowling mode identity and two-sided matching research

Date: 2026-11-11  
Status: PASS_SYNTHETIC_RESEARCH_HELPERS_AND_WEIGHTED_INVARIANT_NOT_PHYSICAL_MODE_CLASSIFICATION  
Scope: analytic two-phase Maxwell fixture methods only; independent of blocked C7X CompOSE crust composition admission.

This package has three layers: (1) slow/rapid interface maps and phase-aware mode diagnostics/branch tracking; (2) the weighted Abel-Liouville determinant invariant and its two-sided matching consequence; (3) one consolidated local-agent DAG, development contract, mutation policy, validation protocol and handoff prompt.

## Main result

For two solutions of the same cold-barotropic Cowling system,
\[
\Omega=h e^{-\nu}(W_1V_2-V_1W_2)
\]
is constant in each smooth phase. Across an interface,
\[
\Omega_o/\Omega_i=(h_o/h_i)e^{-(\nu_o-\nu_i)}\det T_{i\to o}.
\]
The current slow and rapid transfer maps have det(T)=h_i/h_o and the background lapse is continuous, so the concomitant is continuous. This is an algebraic Abel-Liouville identity, not a proved canonical covariant symplectic current.

The same invariant gives a two-sided matching scalar: integrate the regular center solution outward and the surface-compatible line inward, then evaluate their weighted determinant at multiple radii. Its zero set must agree with the current surface shooting condition.

## Contents

- interface_transfer.py: exact slow/rapid maps and residual checks.
- diagnostics.py: positive profile norm, trapping/localization, interface polarization, profile overlap and conservative branch assignment.
- weighted_wronskian.py: bilinear weighted determinant, interface jump diagnostic, surface boundary basis and matching determinant.
- profile_schema.v2.json: background-complete phase-resolved interchange contract with one-sided interface records.
- WEIGHTED_ABEL_LIOUVILLE_AND_TWO_SIDED_MATCHING_20261111.md: derivation, caveats, numerical results and literature boundary.
- LOCAL_AGENT_DAG_DEVELOPMENT_AND_HANDOFF_20261111.md: single local execution authority, DAG, contracts, mutation gates, TDD, validation and anti-drift audit.
- LOCAL_AGENT_DAG.v1.json: machine-readable DAG and gates.
- HANDOFF_PROMPT_20261111.txt: copyable standalone local-agent prompt.
- MODE_IDENTITY_STATUS.json and VALIDATION_REPORT.md: status and exact synthetic scope.
- computation_evidence/: computation manifest and captured test/environment logs.
- tests/: deterministic regression tests.

## Reproduce

From repository root, execute the following:

    PYTHONPATH=research python -m pytest -q research/maxwell_mode_identity/tests
    python -m compileall -q research/maxwell_mode_identity
    python -m json.tool research/maxwell_mode_identity/profile_schema.v2.json >/dev/null
    python -m json.tool research/maxwell_mode_identity/LOCAL_AGENT_DAG.v1.json >/dev/null

Dependencies: Python 3.10+, NumPy, SciPy, pytest.

## Current synthetic evidence

42 deterministic pytest cases passed (32 test functions after parameter expansion). Smooth-phase maximum relative drift was 7.9072e-11. A two-sided synthetic integration gave 4.6733e-12 spread at five matching radii with 100 m maximum step; the same setup with an unlimited step gave 3.5761e-4. The latter is retained as a negative control and mandates max-step convergence in the local plan.

No actual four native fixture-root profiles have been processed by the new helper. No native Fortran source was changed/rebuilt. Zero-jump star comparison and finite-jump continuation are not run. No physical f, p1 or discontinuity-g label is asserted.

## Mutation boundary

The local-agent contract forbids raw EOS, canonical pair manifests, C7X statuses, source-authority solver semantics and production adapter changes. Additive native profile export is limited to M2. Any RHS, transfer map, root, physical gate or public API change requires explicit approval. Push, PR, merge and Dropbox actions are not authorized by the local handoff.


## Numerical characteristic and JSON policy

The raw weighted determinant carries the exact matching zero set when endpoint bases and the ODE depend analytically on omega. The normalized determinant angle is real-frequency diagnostic only; complex norms are generally non-holomorphic. Complex-frequency calculations must use analytic scaling or multiple shooting and preserve all determinant scale factors. QR requires det(Q)*det(R), not Q alone.

If the inner concomitant is zero or below relative_floor, the observed interface ratio is undefined and is returned as null with UNDEFINED_ZERO_DENOMINATOR. No NaN or Infinity is allowed in JSON receipts. Complex scalar fields must be encoded explicitly as real/imaginary pairs.



## 2026-11-11 continuation: closure-mismatch Abel compensation

The continuation adds a conditional piecewise compensated Abel–Liouville invariant for the exact 2x2 Cowling matrix when the perturbation closure chi_ad differs from the equilibrium derivative chi_eq. Read COMPENSATED_ABEL_LIOUVILLE_NOTE_20261111.md for the derivation and scope limitations. The helper and its 12-test synthetic suite are separate from the native solver.

The single local execution contract is still LOCAL_AGENT_DAG.v1.json and LOCAL_AGENT_DAG_DEVELOPMENT_AND_HANDOFF_20261111.md. Revision 1.1 adds node LA9 for native closure audit, piecewise compensation, complex diagnostic serialization and anti-drift validation. HANDOFF_PROMPT_20261111.txt is the standalone invocation text.

No physical nonbarotropic theorem, native q-to-r verification, actual stellar profile, or f/p1/interface-g classification is claimed. The current PR #6 source also has an identified complex-jump reporting loss; LA9 requires an additive compatibility-preserving regression and output fix.
