# COPY-READY PROMPT: NS_Cowling EOS remediation

You are a constrained coding agent working on:

- Repository: `cosmosapjw-quantum/NS_Cowling`
- Branch: `remediation/2026-11-15-m5a-maxwell-eos-contracts`
- Base snapshot: `cloud/eos390-fortran-unification-v1@772e85910afb40aea37c292680be24b363586a54`
- Existing PR #13 is open/draft/unmerged. Do not merge it.
- Execution contract: `handoff/remediation_20261115/execution_manifest.json`
- Implementation plan: `handoff/remediation_20261115/IMPLEMENTATION_PLAN.md`
- Scientific contracts: `handoff/remediation_20261115/SCIENTIFIC_CONTRACTS.md`

## Mission

Finish the audit-driven remediation in the manifest. Some safe patches are already present in the branch; inspect them and do not reimplement them blindly. Run tests, repair reproduced regressions, and complete only work supported by code and numerical evidence. Preserve blockers where a scientific/source decision is required.

## First actions

1. Run `git status --short`, `git branch --show-current`, and `git rev-parse HEAD`.
2. Confirm branch is exactly `remediation/2026-11-15-m5a-maxwell-eos-contracts`.
3. Confirm the ancestry includes `772e85910afb40aea37c292680be24b363586a54`.
4. Read the repository's `AGENTS.md`, `CONTRIBUTING.md`, current workflow and the three handoff files above.
5. Fetch the current remote branch and compare before editing. Never force-push.

## Hard scientific and repository invariants

- Never modify raw CompOSE source EOS files.
- Never weaken the (10^{-12}) relative nucleon-mass gate.
- Never relax mechanical stability, causality, beta-equilibrium or first-law gates.
- Never clip, smooth, rescale or reorder source data to turn a rejected pair into a survivor without a separately approved reconstruction budget and provenance.
- Never claim native verification for all 27 M5A cases from the single headline barotrope.
- Never claim composition is preserved by Maxwell v2 while the native v1 kernel returns thermo-only output.
- Never call the Python Wigner–Seitz mirror a native Fortran end-to-end run.
- Never call the M7 smooth CSS bridge a sharp Maxwell transition. Its audited scope is thermo-only, C⁰, composition dropped, C3/C4 NOT_AVAILABLE.
- Never claim CI PASS without a fresh workflow run whose `head_sha` equals the final branch HEAD.
- Do not merge any PR. Do not force-push.

## Already patched: verify, do not undo

1. M5A native round-trip compares (n_B,\mu_B,P,\epsilon) values with a maximum relative tolerance of (5\times10^{-13}), not merely row count.
2. Receipt-level M5A native verification is scoped to the headline barotrope; individual 27-case native fields are false and marked NOT_VERIFIED_INDIVIDUALLY.
3. Duplicate wrong-name artifact `GPPVA_DD2_m5a_real_css.json` was removed; correct BBB_BHF2 artifact remains.
4. Maxwell v2 rejects `composition_available` and `composition_file_present` input with COMPOSITION_UNSUPPORTED, and no longer reports composition as preserved.
5. Maxwell v2 project root is derived from `__file__`.
6. APR+CSS audit uses `PchipInterpolator.derivative()` and its ADMITTED verdict now requires branch-interior differential first-law residuals below (10^{-3}).

## Execution order

### A. Verify the patches and fix regressions

Run:

```bash
python -m pip install -e '.[test]'
bash fortran_eos/build.sh build/compose_eos_fortran
export PYTHONPATH=python_frontend
export COMPOSE_EOS_FORTRAN_ROOT=.
export COMPOSE_EOS_FORTRAN_EXE=build/compose_eos_fortran/compose_eos_fortran
python -m pytest -q python_frontend/tests/test_m5a_provenance_correction.py
python -m pytest -q python_frontend/tests/test_maxwell_v2_composition.py
python -m pytest -q python_frontend/tests/test_apr_css_maxwell_audit.py
```

Then run the full workflow test command:

```bash
python -m pytest -q python_frontend/tests pipeline/tests
```

Record exact stdout, exit code, skipped count and HEAD SHA. If a failure is unrelated to these patches, identify its root cause and report it; do not hide it or change assertions to green.

### B. Close the APR+CSS production interpolation blocker

The audit script is corrected, but `eos_generator/compose_eos_pkg/compose_eos/phase_transitions.py::build_maxwell_eos` still independently PCHIP-interpolates (P(\mu_B)), (n_B(\mu_B)), and (\epsilon(\mu_B)). This does not guarantee (n_B=dP/d\mu_B).

Required investigation:
1. Reproduce the exact derivative residual on each admitted input phase using the analytic derivative of the actual interpolant.
2. Compare: (i) Hermite (P(\mu_B)) with endpoint slopes (n_B); (ii) shape-preserving (P(\mu_B)) with (n_B=P'); (iii) a thermodynamic-potential reconstruction; (iv) integration of a positive monotone (n_B(\mu_B)).
3. For every method report source-knot residuals in (P,\mu_B,n_B,\epsilon), exact or interval-resolved minimum (n_B), minimum (dn_B/d\mu_B), causality (\mu_B n'_B\geq n_B), and the pressure integration residual.
4. Choose no interpolation method unless it satisfies a predeclared source error budget. If no method does, keep `build_maxwell_eos` admission blocked and mark the artifact REVIEW_REQUIRED.
5. Do not use the Maxwell endpoint identity as a substitute for the branch-interior derivative test.

Do not decide a new error budget from the result you want. If the budget is not specified by source precision/provenance, stop with a precise proposal and evidence.

### C. Complete native Wigner–Seitz integration

The current native module and synthetic unit test are not enough to claim a native production path. Trace `compose_eos_wigner_seitz.f90`, then implement or document a public CLI command that consumes real CompOSE input, calls the native Fortran routine and serializes output with source row keys, Iphase and provenance. Add a SPG_M1 test that invokes that compiled CLI. The Python mirror is a separate test and must be labelled as such.

Do not call the existing algebraic closure reconstruction a full microscopic Thomas–Fermi variational solver unless a native functional minimization and its physical inputs are actually implemented.

### D. Maxwell v2 composition implementation

The current fail-closed guard is correct. Only remove it if all of the following are implemented:
- native operation accepts composition-bearing phases;
- input masks, sparse presence, C3/C4 applicability, Iphase and Nadd semantics have explicit contracts;
- the transition construction defines how composition changes at the phase interface;
- native writer and reader preserve all relevant arrays and metadata;
- tests consume genuine composition-bearing CompOSE input and compare source/output/read-back row-by-row.

Otherwise leave the explicit COMPOSITION_UNSUPPORTED behavior in place.

### E. SPG_M1 Maxwell source inconsistency

Revalidate raw source hashes and rerun the interval audit. The known Iphase-7 segment has 147 intervals and 99 source secants outside the endpoint density bounds. This is not repaired by changing an interpolator while keeping all source knots fixed. Seek higher-resolution source output or another independently admitted crust model. Do not edit raw EOS tables.

### F. Wigner/Maxwell/CSS scope and receipts

Keep the M7 CSS bridge declared as a thermo-only C⁰ bridge; composition dropped and C3/C4 NOT_AVAILABLE. It is not a sharp Maxwell transition and its seam is not C¹. Keep Seidov limited to the small-core onset criterion. Regenerate any audit receipts only after the final code/CI run and stamp them with the exact final Git HEAD, native executable SHA and source data identities.

## TDD and verification rules

- For each new production behavior, add a focused regression test first.
- Observe and record a failing RED run before implementation when the environment permits. In this host session that RED output was not observed; do not invent it.
- Never change a test assertion merely to hide a genuine physical or numerical failure.
- Run the focused tests and then the full workflow suite.
- Run `git diff --check`.
- Verify the final branch diff and all files added/removed.
- Inspect CI job logs, not just the top-level green/red badge.

## Push and final report

Push only to the remediation branch. Create or update a draft PR targeting `research/eos390-fortran-unification-v1`. Never merge.

Return a compact machine-readable completion record with:
`branch`, `head_sha`, `changed_files`, `commit_shas`, `ci_run_id`, `ci_head_sha`, `pytest_passed`, `pytest_skipped`, `pytest_failed`, `fortran_build_exit`, `remaining_blockers`, and `scientific_claims_with_scope`.

A task is complete only when the relevant command was run on the exact final HEAD and its evidence is preserved. If a task is blocked by missing source data, unclear provider semantics, unapproved physical tolerance or unavailable native execution, mark it BLOCKED and state the exact missing item. Do not convert BLOCKED to PASS.
