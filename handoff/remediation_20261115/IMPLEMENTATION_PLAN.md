# NS_Cowling EOS Remediation Implementation Plan

> For agentic workers: execute this plan task by task on the remediation branch. Do not merge a pull request.

**Goal:** finish the audit-driven code corrections while preserving scientific invariants and keeping unresolved source/physics work explicitly blocked.

**Architecture:** small fixes stay in the Python frontend and audit scripts; native Fortran contracts are not weakened. Scientifically unresolved interpolation/composition paths remain separate tasks with fail-closed behavior. Every durable result must be linked to the exact source commit, executable and output artifact.

**Tech Stack:** Python 3.10+, NumPy, SciPy, pytest, Fortran 2018, gfortran, GitHub Actions.

**Source snapshot:** `cloud/eos390-fortran-unification-v1@772e85910afb40aea37c292680be24b363586a54`.

## Global constraints

- Keep mass-pair comparison tolerance at exactly (10^{-12}) relative.
- Preserve mechanical stability, causality, beta-equilibrium and first-law gates.
- Do not mutate raw CompOSE EOS inputs.
- Never mark individual M5A parameter cases as native verified from one headline artifact.
- Do not describe thermo-only Maxwell v2 as composition-preserving.
- Do not promote a Python mirror test to native Fortran end-to-end evidence.
- The current M7 CSS has one constant bridge sound speed, a first-order crust-side jump, and continuous state variables at the core join; do not claim C1 or full phase-interface mode stability.
- Do not merge PR #13 or the remediation PR; do not force-push.

## Review focus

1. **Receipt scope:** one native output versus 27 parameter configurations. Test receipt-level evidence separately from per-case evidence.
2. **Round-trip identity:** same row count is insufficient. Compare (n_B,\mu_B,P,\epsilon) and record residuals.
3. **Sparse composition:** presence masks and C3/C4 applicability must survive segmentation/serialization exactly.
4. **Derivative consistency:** endpoint Maxwell identities do not prove branch-interior (n_B=dP/d\mu_B).
5. **Native execution:** distinguish native Fortran unit tests, Python mirror tests, CLI execution and actual raw-source end-to-end output.

---

## Task 1: M5A native evidence scope and read-back identity

**Files:** `scripts/maxwell_audit/m5a_native_readback_verify.py`, `download/m5a_real_css/E390-05A_M5A_REAL_CSS_RECEIPT_20261114.json`, `python_frontend/tests/test_m5a_provenance_correction.py`.

**Implemented changes:**
- Require write, read, row-count match and value match before top-level native verification is true.
- Compare (n_B,\mu_B,P,\epsilon) with maximum relative residual (leq5\times10^{-13}).
- Scope top-level verification to `HEADLINE_BAROTROPE_ONLY`.
- Keep all individual case verification flags false and `NOT_VERIFIED_INDIVIDUALLY`.
- Preserve a list of verified case keys, currently empty.

**Verify:**
```bash
PYTHONPATH=python_frontend COMPOSE_EOS_FORTRAN_ROOT=. COMPOSE_EOS_FORTRAN_EXE=build/compose_eos_fortran/compose_eos_fortran python scripts/maxwell_audit/m5a_native_readback_verify.py
PYTHONPATH=python_frontend COMPOSE_EOS_FORTRAN_ROOT=. COMPOSE_EOS_FORTRAN_EXE=build/compose_eos_fortran/compose_eos_fortran python -m pytest -q python_frontend/tests/test_m5a_provenance_correction.py
```

**Acceptance:** all four primary arrays pass the declared relative tolerance; no unverified parameter case is marked true.

## Task 2: Maxwell v2 composition contract

**Files:** `python_frontend/compose_eos/maxwell_v2.py`, `python_frontend/tests/test_maxwell_v2_composition.py`.

**Implemented changes:**
- Derive repository root from `__file__`, not a developer-specific absolute path.
- Reject `composition_available` or `composition_file_present` input with `COMPOSITION_UNSUPPORTED`.
- Reclassify v2 as Iphase-segmented thermo-only orchestration.
- Do not claim composition is preserved by the successful thermo-only path.

**Verify:**
```PYTHONPATH=python_frontend COMPOSE_EOS_FORTRAN_ROOT=. COMPOSE_EOS_FORTRAN_EXE=build/compose_eos_fortran/compose_eos_fortran python -m pytest -q python_frontend/tests/test_maxwell_v2_composition.py
```

**Acceptance:** composition-bearing source input is explicitly rejected; thermo-only segmentation still reaches the unchanged v1 physical gates.

## Task 3: APR+CSS differential first-law audit

**Files:** `scripts/maxwell_audit/audit_apr_css_maxwell.py`, `python_frontend/tests/test_apr_css_maxwell_audit.py`.

**Implemented changes:**
- Replace `np.gradient(Pchip(mu), mu)` with `PchipInterpolator(mu, P).derivative()(mu_grid)`.
- Make the audit admission verdict require the APR and analytic CSS branch residuals both below the declared (10^{-3}) relative threshold.
- A failed APR source/interpolant check returns `REVIEW_REQUIRED`, even when the endpoint jump identity and TOV/Love example pass.

**Verify:**
```PYTHONPATH=python_frontend COMPOSE_EOS_FORTRAN_ROOT=. COMPOSE_EOS_FORTRAN_EXE=build/compose_eos_fortran/compose_eos_fortran python -m pytest -q python_frontend/tests/test_apr_css_maxwell_audit.py
```

**Not yet implemented:** the production `build_maxwell_eos` still builds independent PCHIP interpolants for pressure, density and energy. Do not assume the audit patch changes that generator. A thermodynamically consistent production interpolant remains an open task.

## Task 4: M7 single-speed CSS with a crust-side first-order jump

**Files:** `fortran_eos/src/compose_eos_matching.f90`, `fortran_eos/app/compose_eos_fortran.f90`, `python_frontend/compose_eos/api.py`, `python_frontend/tests/test_native_crust_core_matching.py`, `scripts/m7/m7_batch_orchestrator.py`.

**Implemented changes:**
- Replace the former two-segment/two-speed CSS construction with one constant (c_s^2\in(0,1]).
- Solve (c_s^2) from the pressure integral between the bridge-side transition state and the selected core endpoint.
- Require a positive crust-side density jump and the exact jump identity (Deltaepsilon=mu_tDelta n).
- Insert two adjacent output rows with equal (P_t,mu_t) and different (n_B,epsilon); then match the bridge continuously to the core state.
- Expose the transition indices, jump values, one bridge speed, and Seidov small-core onset diagnostics through the native JSON and Python report.
- Preserve thermo-only scope and explicit composition dropping. Do not claim C1 at the core-side join if the source derivative differs.

**Verify:**
```bash
PYTHONPATH=python_frontend COMPOSE_EOS_FORTRAN_ROOT=. COMPOSE_EOS_FORTRAN_EXE=build/compose_eos_fortran/compose_eos_fortran python -m pytest -q python_frontend/tests/test_native_crust_core_matching.py
```

**Acceptance:** one (c_s^2) over all bridge samples; adjacent crust/bridge rows with common (P_t,\mu_t), positive (Delta n,Deltaepsilon), and the jump identity; continuous (n_B,P,\epsilon,\mu_B) at the core join; exact report indices; no composition-preservation claim. The prior two-speed CSS artifact remains immutable historical evidence.

---

## Task 4: Native Wigner–Seitz production route

**Files to inspect:** `fortran_eos/src/compose_eos_wigner_seitz.f90`, `fortran_eos/app/compose_eos_fortran.f90`, `fortran_eos/src/compose_eos_io.f90`, `python_frontend/compose_eos/native.py`, `python_frontend/tests/test_wigner_seitz_reconstruction.py`.

**Required work:**
1. Identify exact public routine signatures and `eos_data_t` ownership of the input composition.
2. Design a native CLI/API operation that reads real CompOSE input and returns a provenance-bound reconstructed composition.
3. Preserve original row identity, `Iphase`, masks and provider provenance.
4. Add a test which invokes the native CLI on staged SPG_M1 source data and checks all row-level outputs.
5. Keep synthetic unit tests as unit tests only.

**Stop:** do not label the current module as a full microscopic Thomas–Fermi variational solver; do not claim end-to-end native acceptance without running the native command on the raw source input.

## Task 5: Composition-aware Maxwell v2 (future native implementation)

The current branch deliberately rejects composition-bearing input. A genuine implementation requires an explicit native kernel contract for composition, row-level masks, `Iphase`, phase metadata and round-trip serialization. Do not remove the fail-closed guard until the native path exists and the tests exercise real composition-bearing inputs.

## Task 6: SPG_M1 source convexity

Reconfirm source blob identities and the 147 interval diagnostics. Do not change source data. Seek higher-resolution pressure/chemical-potential/density output or use an independently validated crust model. If a thermodynamic-potential reconstruction is explored, report all source residuals and reject if it exceeds the predeclared source error budget.

## Task 7: Full verification and handoff

Run the exact CI workflow commands from `.github/workflows/fortran-eos-ci.yml`. Create or update a **draft** PR targeting `research/eos390-fortran-unification-v1`. Confirm the workflow run's `head_sha` equals the final branch HEAD. Report pass/skip/fail counts and all remaining blockers. Never merge.

---

## Completion checklist

- [ ] Focused M5A, Maxwell v2 and APR+CSS tests pass.
- [ ] Fortran executable builds with the workflow's gfortran flags.
- [ ] Full `python -m pytest -q python_frontend/tests pipeline/tests` passes.
- [ ] The M5A native read-back values are independently compared, not just counted.
- [ ] No case-level native verification is claimed without a case-specific run.
- [ ] Composition-bearing Maxwell v2 remains fail-closed unless native composition support is actually complete.
- [ ] APR+CSS is not admitted when exact branch differential consistency fails.
- [ ] Wigner–Seitz scope says exactly what was executed.
- [ ] SPG_M1 raw source bytes and all physical gates remain unchanged.
- [ ] Final CI run SHA equals the final branch HEAD.
- [ ] No PR has been merged.
