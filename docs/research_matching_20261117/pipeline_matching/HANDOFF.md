# Pipeline and matching-admission audit, 2026-11-17

Scope: selected C2–C9 public native matching/API/receipt boundaries. Baseline Git HEAD `478c7b5de1514ac512f91156f36551478f36f688`. Production files are frozen for parent integration. No Git mutation, no new native build, no recursive reviewers, one compute thread. Existing native executable `build/spg_finish/compose_eos_fortran`, SHA256 `620164cc7fc0c085b3b3541b9529163c2fa4ebf8dc3aae7c9b2caa52ff94d8cb`.

## Observed defects and repairs

1. **Cold columns lost during independent projection.** Finite Q2, Q5 or Q6 corruption was admitted because rows were projected only to n/P/epsilon/mu. Projection now preserves Q2/Q5/Q6/Q7. Admission applies the native method-specific entropy gate (CSS 1e-10; Maxwell 2e-8), Q5 gate 1e-8 with native Il applicability, Q6−Q7 gate 2e-6, finite temperature |T|≤1e-12, and explicit beta-equilibrium metadata. No physical threshold was loosened.
2. **Wrong report identity and count types.** Error-status or wrong-operation receipts containing otherwise valid fields were admitted; Boolean `pressure_crossings=True` passed as integer one. Status/operation are now necessary checks, and report indices/counts require exact integral types excluding bool. Status remains insufficient by itself.
3. **Native output phase contract lost.** Independent readback admission discarded Iphase. It now checks the explicit native-output layout: CSS thermodynamic-only rows have label 0; native Maxwell output has contiguous 1/2 blocks at the reported boundary. This is not a restriction on generic provider phase IDs. Provider phase 9 is accepted when constant and otherwise valid.
4. **Stellar manifest trusted a status string.** Hash-correct EOS files plus `admission.status=ADMITTED` could enter the stellar path without a pair receipt; changed match-report fields were not bound. The driver now verifies `pair_report_sha256`, completed/survived pair identity, matching manifest/receipt fields, and consistent structured admission checks before the existing EOS-byte checks. These hashes provide within-run binding, not authentication against an actor who rewrites every receipt.
5. **Direct pair source gaps misclassified.** `run_pipeline` previously attempted backend execution for absent source files and could classify a missing input as environment/unclassified. It now returns `NOT_EVALUATED_SOURCE_INCOMPLETE`, `data_provenance`, with the missing inventory, before creating an attempt directory. Existing mocked downstream tests now provide complete mock input inventories.
6. **Declared interpolation support ignored.** Both smooth/CSS matching and sharp `construct_maxwell_transition` now read explicit SUPPORT declarations. For `cross_phase_interpolation=UNAPPROVED`, input must be a constant opaque phase-ID branch. For `cross_segment_interpolation=UNAPPROVED`, the supplied density grid must be one contiguous slice of the original segment's density rows (native I/O numerical tolerance 5e-13). This bounds the native window search and endpoint stencils to the preselected source branch. It is an implementation/input-support limitation, not a physical no-transition result or a certificate of a continuous mesh limit. Inputs without SUPPORT retain historical behavior. Separate outer/lower source segments inherit a parent's SUPPORT when the native segment declaration exists. Edited/deleted SUPPORT is rejected when RESULT.json declares its SHA256. Support content/hash and the selected rows are propagated into API/pipeline evidence.

## Reproduction and verification

Preserved logs in `data/research_matching_20261117/pipeline_matching/`:

- `before.txt`: 15 failing adversarial regressions for cold columns, report identity/counts, warm metadata, phase layout and stellar receipt binding.
- `before_missing_source.txt`: one failing direct-call provenance classification, pytest exit 1.
- `before_support.txt`: two failing support bypasses and two phase-9/no-declaration controls passing.
- `before_support_seal_corrected.txt`: sharp-Maxwell and two sealed-support bypasses fail before correction. Earlier `before_support_seal.txt` is preserved but its sharp-Maxwell test lacked the required caller local-neutrality argument; it is not the authority for that reproduction.
- `verification.txt`: **119 passed, pytest exit 0**, including real native normalized CSS match/write/readback and a retained actual generated TF source rejected at source_support/input_contract. The actual TF declaration is preserved in `source_scope`; the test does not reclassify its sparse sample grid as physical impossibility.

Command (all numerical thread variables 1):

```
OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 \
PYTHONPATH=python_frontend:. \
COMPOSE_EOS_FORTRAN_EXE=$PWD/build/spg_finish/compose_eos_fortran \
python -m pytest -q \
 python_frontend/tests/test_matching_support_20261117.py \
 pipeline/tests/test_matching_pipeline_audit_20261117.py \
 pipeline/tests/test_match_report_admission.py \
 pipeline/tests/test_core_crust_matching_pipeline.py \
 python_frontend/tests/test_cowling_reaudit_contracts.py
```

`git diff --check` returned exit 0 after the production freeze. Existing Maxwell analytic fixtures were corrected to the actual native 1/2 output labels; metadata fixtures now declare their cold temperature and beta flag. Neither change lowers a tested gate. Existing stellar mocks now include the pair receipt contract they model.

## Changed dependencies and remaining scope

Production: `scripts/match_core_crust.py`, `scripts/match_report_admission.py`, `scripts/run_css_survivor_stars.py`, `python_frontend/compose_eos/api.py`.

New regressions: `pipeline/tests/test_matching_pipeline_audit_20261117.py`, `python_frontend/tests/test_matching_support_20261117.py`. Updated fixtures: `pipeline/tests/test_match_report_admission.py`, `pipeline/tests/test_core_crust_matching_pipeline.py`, `python_frontend/tests/test_cowling_reaudit_contracts.py`.

No generation physics, stellar solver equation, native matching polynomial, microscopic model, mass normalization or physical acceptance tolerance was changed. Python here performs explicit independent admission and input/receipt validation; production matching remains native. Parent must run combined fresh-native integration after the native agent's edits. No claim is made that the old executable verifies the new native matcher. Legacy `compose_cowling.intake.prepare_from_matched` is a separate archival classifier/CSV route, not the current CSS census-to-star driver; it was read but not modified in this bounded pass.

Skill context actually read: PhysMath Research Loop 1.3.0 main, POLICY_INTEGRATION.md, SUBAGENT_MODE.md, CODING_LOOP.md; systematic-debugging and verification-before-completion. The parent owns integration/C10; this lane returns the implemented and executed C2–C9 boundary repairs without starting another review cycle.
