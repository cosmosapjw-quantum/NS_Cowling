# Focused EOS integration regression — 2026-11-23

Execution used `COMPOSE_EOS_FORTRAN_EXE=/tmp/ns_eos_reaudit_20261123/build/compose_eos_fortran`, SHA256 `d1600879b3c741bb69f0acd150f05ff129048f244723933b1bb1e79518316588`. OMP/OpenBLAS/MKL/BLIS/vecLib/NumExpr threads were pinned to one. Python performed orchestration and existing test-fixture duties; all new production numerical implementation remains Fortran.

## Final observed result

- Focused pytest: **200 passed, 1 skipped in 19.50 seconds**, exit 0.
- All target source/test files and the pinned binary were byte-identical before and after that run (`FOCUSED_RUN.json`).
- Four existing native scientific programs passed: matching thermodynamic adversarial controls, interval counterexample, physical slice conventions, and matching polynomial extrema controls.
- Two additional native programs passed: quartic polynomial extrema and Maxwell reference-mass normalization. The latter preserves physical energy with maximum observed residual 2.8421709430404007e-14 MeV fm^-3 and rejects a physical +100 MeV/baryon offset as no common chemical-potential support.
- No full TOV campaign was invoked.

`FOCUSED_RUN.json` records the full command, environment, test selection, source/binary SHA256 before and after, and exit status. `FOCUSED_PYTEST.log` and `FOCUSED_PYTEST.xml` contain the final detailed evidence. Native controls have per-build/per-run logs and command/hash receipts under `native_matching_controls` and `extra_native_controls`.

## Initial failures retained and resolved

The initial run was **4 failed / 192 passed / 1 skipped in 19.29 seconds**. Its complete log, XML and hashes are preserved as `INITIAL_FOCUSED_*`, and pre-correction test files as `PRE_test_*.py`. Production gates and tolerance values were not modified in response to these failures.

1. **An overlap fixture incorrectly equated an old endpoint-search failure with physical impossibility.** Its crust/core density domains overlap, but exhaustive native search finds another admissible CSS pair: zero-based i1=72, i2=78 after 336 rejected candidates. The density and energy jumps are positive, sound speed squared is approximately .00215100453, and native serialized C1–C11 verification passes. The regression now asserts the retry, valid endpoint ordering, positive jump, pressure/chemical-potential continuity and full serialized verification. A separate negative control places every core knot below every crust knot, so all candidate density orderings are invalid and exhaustive search must reject.

2. **Two phase-stencil fixtures invalidated only the previously chosen endpoint.** The native retry correctly moves the crust endpoint from 299 to 297 and the core endpoint from 79 to 81. Each new retained-side derivative stencil has six rows of one contiguous phase, without borrowing across the artificial boundary. The corrected tests assert this behavior and the causal source derivative. Two additional negative controls partition the entire crust/core into phases of at most two rows; all possible endpoint derivative stencils must then fail. They reject without writing a combined EOS.

3. **The old composition-absent readback fixture used an acausal full APR source.** Native source verification reports C10 failure, 24 acausal points and maximum c_s^2=1.470861145390916. The revised retained-tail gate correctly prevents every match from inheriting that source tail. This exact APR rejection is retained as a new regression. The positive composition-absent readback control now uses the real causal GM_GM1 thermodynamic files copied byte-for-byte while explicitly omitting only its composition sidecar. SPG_M2→GM_GM1 CSS produces 659 rows, native serialized readback returns 659 rows, and C1–C11 all pass. Both the source verification and complete readback are asserted.

All root-cause receipts and the serialized positive-control source/output are preserved under `rootcause/`. The historical APR table was not modified or truncated.

## Test-file changes authorized after the initial failure report

- `python_frontend/tests/test_native_crust_core_matching.py`
- `python_frontend/tests/test_native_matching_admission.py`
- `python_frontend/tests/test_e390_07_readback.py`

These changes replace mistaken endpoint-specific impossibility assumptions with the intended physical and search contracts, preserve negative acceptance tests, and add all-window-invalid controls. No production source, physical gate, tolerance, source table or catalog was changed by this integration worker.

## Scope and remaining limitation

The one skipped historical M7 test requires `inputs/raw_compose` from an M0 staging step absent in this checkout. Its skip reason is preserved in XML. The real causal composition-absent readback control above executed successfully; the skip is not counted as a pass.

No fixture exhibiting a locally certified first bridge followed by whole-output C11 rejection was observed in this targeted verification. That observation is limited to these tests and does not certify every possible endpoint window or every campaign source. The all-pairs campaign retains its separate whole-output admission gate and failure census.
