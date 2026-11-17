# EOS reconstruction and Maxwell validation Implementation Plan

> **For agentic workers:** execute task-by-task under the existing C4 machine DAG. Do not create a new workflow DAG. Use TDD for every production-code change.

**Goal:** validate cold single-phase EOS reconstruction, preserve the six smooth crust-core matching semantics, and harden the genuine one-charge Maxwell kernel against CompOSE Il-dependent Q3 convention errors.

**Architecture:** native Fortran owns interpolation/reconstruction, polynomial interval extrema, phase validation, coexistence search and output. Python calls the native executable and provides independent end-to-end tests. The existing C4 JSON DAG on research/sharp-maxwell-cohort-admission-v1 remains the sole machine-readable downstream workflow.

**Tech Stack:** Fortran 2018/gfortran, Python >=3.10 (verified 3.13.5), NumPy, SciPy, pytest.

**Spec:** docs/superpowers/specs/2026-11-11-eos-reconstruction-maxwell-design.md

## Global Constraints

- Preserve raw EOS source bytes and source-authority pin a09ec30d74b363652d2d1aefe990557c468679dd.
- Preserve research/fortran-eos-gap-fill-v1 as PR #1 implementation branch; do not merge to main.
- Preserve PR #4 C4 JSON as the only DAG; append nodes, never create a parallel DAG.
- Historical match_crust_core(method=maxwell) remains a smooth plateau and never emits an exact density jump.
- Genuine Maxwell assumes one globally conserved baryon charge and caller-confirmed local charge neutrality in each phase; it is not Gibbs.
- Keep Il=0 and Il=1 Q3 conventions distinct; root search uses Euler mu_E.
- Keep strict reconstruction/source-admission Q3 tolerance 2e-6 distinct from true-Maxwell v1 local tolerance 5e-4.
- Keep C6/C6B/C7 composition/provider blockers; no 65-pair root search or production TOV/Cowling run in this increment.
- Do not insert, modify, reorder or regenerate raw source rows.

## Review Focus

- Il=0 Q3 is charge-shifted while Il=1,Q5≈0 Q3 equals Euler mu. Test both.
- A valid coexistence point must have equal P and Euler mu while density and epsilon jump. Check latent identity.
- Tangential and multiple pressure crossings must not be silently lost or selected.
- Whole-interval causality/stability must be certified by polynomial extrema, not finite sampling.
- Historic maxwell smooth method must keep zero density jump.

---

### Task 1: Il-dependent Q3 convention in genuine Maxwell input validation/output

**Files**
- Modify fortran_eos/src/compose_eos_phase_transition.f90
- Modify python_frontend/tests/test_native_maxwell_transition.py
- Evidence in docs/ADVERSARIAL_AUDIT_MAXWELL_TRANSITION_20261107.md

**Interface:** existing construct_maxwell_transition(low, high, assume_local_charge_neutrality=True, return_report=True) is unchanged.

- [x] Write failing test test_native_maxwell_respects_il0_compose_q3_charge_potential_convention.
- [x] Verify RED: valid Il=0 analytic source phase rejected at Q3/Euler input check.
- [x] Implement minimal fix: recover Euler mu by adding Yq*mn*Q4 to stored Q3 for Il=0; subtract the same term when writing transition and copied source rows.
- [x] Verify GREEN: output Q3 and stored Q3 potential obey Il=0 convention; native C2 and overall checks pass.
- [ ] Commit/update on existing PR #1 branch. Local exported artifact has no .git metadata.

### Task 2: Reconstruction and smooth matching semantics audit

**Files**
- Read: fortran_eos/src/compose_eos_reconstruct.f90
- Read: fortran_eos/src/compose_eos_matching.f90
- Update docs/FORTRAN_EOS_MASTER_RECONSTRUCTION_V1.md
- Update docs/FORTRAN_EOS_CRUST_CORE_MATCHING_V1.md
- Update docs/FORTRAN_EOS_FEATURE_PARITY.md

- [ ] Re-derive F(log n_B), P, epsilon, mu_E and c_s^2 from native source.
- [ ] Verify Il=0/Il=1 Q3 handling in native output and C2 checks.
- [ ] Replay all six smooth methods and assert historic maxwell is smooth_near_pressure_plateau with zero jump.
- [ ] Retain known unported global-spline selection and composition-aware matching as explicit partial scope.

### Task 3: Extend the existing C4 machine DAG

**Files**
- Modify PR #4 branch: data/sharp_maxwell_c4_execution_plan.v1.json
- Update docs/SHARP_MAXWELL_C4_EXECUTION_DAG_V1.md

- [x] Add C4-D3R reconstruction/smooth-matching/true-Maxwell validation node.
- [x] Make Il-aware Q3 relation explicit in C4-D3.
- [x] Make C4-D8 depend on both C4-D7 and C4-D3R.
- [x] Add M0-M5, anti-drift audit and failure taxonomy without weakening the C6B blocker.
- [ ] Execute C4-D3R at the live exact PR #1 head and upload the matching receipts.
- [ ] Keep all 65 pair roots NOT_RUN until source and provider admission gates pass.

### Task 4: Final verification and provenance

**Files**
- Update docs/ADVERSARIAL_AUDIT_MAXWELL_TRANSITION_20261107.md
- Update docs/FORTRAN_EOS_GAP_FILL_TESTS.log
- Update CHANGELOG_v4.md
- Update SHA256SUMS_v4.txt
- Add the EOS design, this plan and EOS_RECONSTRUCTION_MAXWELL_HANDOFF_20261111.txt

- [x] Native build passed with GNU Fortran 14.2.0.
- [x] Full frontend/pipeline suite passed: 68 passed in 13.57s.
- [x] Native polynomial-extrema test passed.
- [ ] Refresh SHA256 manifest after all final documentation and DAG edits; require all entries to pass.
- [ ] Re-fetch PR #1/PR #4 exact heads and CI runs; do not claim CI PASS until the new exact head passes.

## Execution provenance

The local run originated from PR #1 Actions artifact run 42445584718 at exact branch head ca71d8e394f9043d778676c941f88ce01ca7faa9. The exported source tree contains no .git metadata. Therefore the 68-test/build result is an artifact-based local execution, not a committed worktree replay. The RED run was captured before the native change; the final GREEN run included the patch.
