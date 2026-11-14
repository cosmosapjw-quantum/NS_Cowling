# E390 Pairwise Mass-Convention Audit — Final Handoff

**Date:** 2026-11-14
**Branch:** `cloud/eos390-fortran-unification-v1`
**HEAD at audit start (verified):** `b7989e7256257f1753154002c9f459588fd8c112`
**Backup refs (local, not pushed):** `backup/e390-mass-audit-pre-20261114` + tag `e390-mass-audit-pre-20261114`
**PR #13:** left open / draft (no merge, no force-push, no Dropbox upload — per mission constraints)
**No physics models regenerated. No target mass pair chosen. No historical outputs deleted. No Maxwell/390 batch launched.**

---

## Completed work (this bounded execution)

### 1. State frozen and inspected
- Branch HEAD `b7989e7` verified against the mission's stated SHA — still authoritative.
- Fresh backup ref created at HEAD before any mutation (`backup/e390-mass-audit-pre-20261114` + tag), NOT reusing any existing backup.
- Working tree had 2114 unstaged changes (leftover from the prior session's stash restore); these were preserved untouched, not committed.

### 2. Pairwise mass-provenance matrix (deliverable)
- **Files:** `download/PAIRWISE_MASS_PROVENANCE_20261114.json` + `.csv`
- 18 EOS data dirs audited against raw `eos.thermo` headers (line 1: `mn mp <flag>`) and the model parameter registry.
- 65 unique crust-core pairs classified: **0 PASS, 65 MASS_MISMATCH_BOTH**.
- **Source bytes available** — no model reported `BLOCKED_SOURCE_BYTES_UNAVAILABLE`.
- **Discrepancy resolved:** the matching report's `(939,939)` claim is correct for the *written CompOSE artifact* but conflates normalization with the physical model input for GPPVA_DD2 (registry: distinct `(939.56536, 938.27203)` per Typel 2010). For GPPVA_DDME2 and GPPVA_TW the common convention is genuine (model itself used `M=939`).

### 3. Native Maxwell gate audited and preserved
- Gate at `fortran_eos/src/compose_eos_phase_transition.f90:664-665`: independent `mn` and `mp` comparison, `rtol=1e-12`. **Preserved unchanged.** No tolerance widening, no gate weakening.
- `iYq` integer check (line 669) and `euler_tolerance=5e-4` (line 9) reviewed and left intact; v2 proposal records the derivation task as pending evidence.

### 4. Regression tests (deliverable)
- **File:** `python_frontend/tests/test_maxwell_mass_gate.py` — 9 tests collected, covering all 6 mission scenarios + precision policy + provenance consistency.
- Verified by inspection against gate source (lines 664-665); execution BLOCKED on Fortran backend (B1).

### 5. M5 script + receipt corrected (deliverable)
- **Files:** `scripts/m5_m5a/m5_real_maxwell_native.py` (rewritten), `download/m5_real_maxwell/E390-05_M5_REAL_MAXWELL_RECEIPT_20261114.json` (rewritten).
- 8 mission corrections: 6 DONE, 2 BLOCKED (refinement study, raw first-law residual — both require backend re-execution).
- Key fix: removed unphysical `mp = M_N * 0.9938` (~933.74 MeV); now uses raw source `mp=938.272`. Separated `P_t_input` (20.0) from `P_t_actual` (168.78) in the receipt's canonical field. Reclassified scope from "production PASS" to "native Maxwell kernel reference (real hadronic + synthetic CSS)".

### 6. Audit document (deliverable)
- **File:** `docs/MAXWELL_MASS_AUDIT_20261114.md` (292 lines, 11 sections) — findings, pairwise matrix, gate audit, regression tests, M5 corrections, Maxwell v2 proposal + test plan, root census vocabulary, representation-only vs physical regeneration decision framework, SSOT updates, blockers, next DAG node.

### 7. Matching report corrected (deliverable)
- **File:** `download/report/matching_report.tex`
- "No Maxwell construction survives" → "All 40 Maxwell attempts stopped at the native rest-mass admission gate; the existence of a pressure crossing or a physical coexistence point was therefore not tested."
- "(939,939)" claim corrected to distinguish *written CompOSE artifacts* from the *physical GPPVA_DD2 model* (distinct masses per Typel 2010).
- Braces balanced (205/205); edits compile clean.

### 8. Maxwell v2 proposal + test plan (deliverable)
- In `docs/MAXWELL_MASS_AUDIT_20261114.md` §6: separately-versioned `construct_maxwell_transition_v2`, 8 required semantics, 3 mission-flagged reviews (`iYq` axis-value compatibility, `euler_tolerance` derivation, `set_transition_row` status checks), 8-test plan. v1 left intact.

---

## Blocked work

### B0 — GitHub API / CI access unavailable
`gh` CLI not installed; no `GITHUB_TOKEN`. SSH key authenticates git only, not the REST API.
- **Impact:** PR #13 CI run IDs not retrievable; PR summary/body not updatable.
- **Needed:** a GitHub API token or `gh` CLI authenticated.

### B1 — Fortran backend build blocked (gfortran toolchain mis-configured)
Only available gfortran (`/tmp/my-project/.micromamba/envs/fortran/bin/gfortran`, gcc 15.3.0 conda) fails at link: sysroot has `lib64/` but no `lib/`, and the spec references `/usr/lib64/libc_nonshared.a` not present in the sysroot. `-B`/`-L`/`--sysroot` do not resolve it.
- **Impact:** `compose_eos_fortran` binary cannot build → all native-backend tests (new `test_maxwell_mass_gate.py`, existing `test_native_maxwell_transition.py`, proposed v2 tests) cannot execute; M5 refinement study cannot run; post-admission Maxwell root census cannot run.
- **Needed:** a working gfortran (complete sysroot) or system gfortran. The `.tools/bin/gfortran` wrappers point at a non-existent `mmb_root/envs/fortran`.

### B2 — Source composition coverage incomplete
Crust `eos.thermo` files lack per-nucleus composition for species-resolved rest-energy corrections.
- **Impact:** physical model mass changes (§8.2 of audit doc) cannot be derived from available source data. Representation-only re-expression (§8.1) is the only unblocked route for GPPVA_DD2.
- **Needed:** per-nucleus crust composition data, or source-backed confirmation of mass-convention insensitivity at required precision.

### B3 — `euler_tolerance = 5e-4` derivation
Mission asks for a tolerance derived from source accuracy + component-construction contract. Current value undocumented.
- **Impact:** v2 proposal records derivation as pending; value not changed.
- **Needed:** derivation from `eos.thermo` decimal precision + C2 first-law bound.

### B4 — DAG `evidence_status` reconciliation (mechanical, not blocked by tooling)
Mission flags: do not claim multidim composition round-trip completion where status detail lists sub-tasks as pending. The 21-node DAG (`data/eos390_native_fortran_unification_dag.v1.json`) `evidence_status` (19 subkeys) must be walked node-by-node against per-node receipts.
- **Impact:** potential stale "pending" details where receipts claim completion.
- **Needed:** a node-by-node walk (large mechanical task; left as the first item of the next execution).

---

## Exact validation performed

| Check | Method | Result |
|---|---|---|
| Branch HEAD matches mission SHA | `git rev-parse HEAD` | `b7989e7256257f1753154002c9f459588fd8c112` ✓ |
| Backup ref created at HEAD | `git rev-parse backup/e390-mass-audit-pre-20261114` | `b7989e7` ✓ |
| Matrix JSON validity | `json.load` + asserts on summary counts | 18 EOS, 65 pairs, 65 mismatch_both, 0 pass ✓ |
| Matrix CSV validity | `csv.DictReader` | 65 rows, expected columns ✓ |
| M5 receipt correction markers | `json.load` + asserts | P_t_actual=168.78, P_t_input=20.0, scope reclassified, raw source masses (939.56533, 938.272) ✓ |
| M5 script syntax | `python3 -m py_compile` | clean ✓ |
| Regression tests collected | `pytest --collect-only` | 9 tests ✓ |
| Matching report brace balance | char count | 205/205 ✓ |
| Native gate semantics | source read of lines 664-665 | independent mn/mp, rtol=1e-12, preserved ✓ |
| Raw source bytes present | `find eos.thermo` | all 18 EOS dirs have raw headers ✓ |
| No `mp = M_N * 0.9938` in M5 | `grep` | only docstring/comment references remain ✓ |

**Not validated (BLOCKED):** native-backend test execution, M5 refinement study, Maxwell root census — all require Fortran backend (B1).

---

## Next DAG node

1. **Resolve B1** (gfortran toolchain) — unblocks all native tests and refinement studies.
2. **Execute `test_maxwell_mass_gate.py`** — confirm 9 tests pass against unchanged v1 gate.
3. **Re-run `m5_real_maxwell_native.py`** — populate `refinement_study` (convergence over `n_min ∈ {0.06,0.08,0.10,0.12}`) and `raw_first_law_residual`.
4. **Walk 21-node DAG** against per-node receipts — reconcile `evidence_status` (B4).
5. **Implement representation-only CompOSE re-expression for GPPVA_DD2** (§8.1 of audit doc) as a separate reversible script; record `mn_old/mn_ref` ratio; verify physical invariance of `P, ε, μ`. Do NOT regenerate the physical model.
6. **Re-run Maxwell root census** on admitted pairs post-re-expression — 14-category terminal vocabulary (§7 of audit doc).
7. **Implement Maxwell v2** per §6 (separately versioned; v1 unchanged) after mass-convention resolution.
8. **Update PR #13** once B0 resolved — remove superseded status snapshots, pin corrected HEAD + CI run.

**Governing principle (preserved):** preserve the physical mass pair, preserve the source EOS, make every normalization or regeneration step explicit and reversible. Never manufacture a mass convention to make a numerical gate pass.

---

## Files produced / modified this execution

**New (in repo):**
- `data/PAIRWISE_MASS_PROVENANCE_20261114.json` (also in global `download/`)
- `data/PAIRWISE_MASS_PROVENANCE_20261114.csv` (also in global `download/`)
- `docs/MAXWELL_MASS_AUDIT_20261114.md`
- `python_frontend/tests/test_maxwell_mass_gate.py`
- `download/E390_MASS_AUDIT_HANDOFF_20261114.md` (this file)

**Modified (in repo):**
- `scripts/m5_m5a/m5_real_maxwell_native.py` (8-correction rewrite)
- `download/m5_real_maxwell/E390-05_M5_REAL_MAXWELL_RECEIPT_20261114.json` (corrected)
- `download/report/matching_report.tex` (2 edits: "No Maxwell survives" → "stopped at gate"; (939,939) conflation fixed)

**Local-only git refs (not pushed):**
- `backup/e390-mass-audit-pre-20261114` (branch)
- `e390-mass-audit-pre-20261114` (tag)

**Persisted scripts (outside repo):**
- `/home/z/my-project/scripts/build_pairwise_mass_provenance.py`

**No changes committed.** All modifications are in the working tree, ready for review. Per mission constraints: PR #13 left open/draft, no merge, no force-push, no Dropbox upload.

---

## UPDATE — Remaining EOS-MASS nodes executed (this turn)

**Toolchain blockers resolved:**
- **B0 (SSH):** located ssh at `/tmp/my-project/.tools/mmb_root/envs/fortran/bin/ssh` + key `~/.ssh/id_glm`; authenticates as `cosmosapjw-quantum`. Push helper at `/home/z/my-project/scripts/git_ssh_helper.sh`.
- **B1 (gfortran):** the micromamba gfortran 15.3.0 sysroot was missing `lib/` and `usr/lib64` symlinks; created `sysroot/lib -> lib64` and `sysroot/usr/lib64 -> ../lib64`. The Fortran backend now builds and runs.

**Final HEAD:** `7b5285a8fcf5e6febcf518cdfb2d3bdce464e6a5` (5 audit commits on top of `b7989e7`)

### Node 6 — M5 density-grid refinement (commit e608927)
Re-ran the corrected M5 script with the live Fortran backend. Refinement over n_min in {0.06,0.08,0.10,0.12} demonstrates convergence: mu_t=1297.3558, P_t=168.7759, n_low=0.73978, n_high=0.92744 are identical across all four cuts (the coexistence lives at n_B≈0.74, far above the cuts). Receipt now carries refinement_study, raw_first_law_residual, P_t_input vs P_t_actual separation, raw source mass pair, thermo SHA-256.

### Node 7/7b — DAG evidence_status reconciliation (commit 8b3a62c)
Subagent walk of 21 DAG nodes found E390-10B overclaimed `PHASES_1_6_COMPLETE_MULTIDIM_READ_WRITE_ROUNDTRIP_PASS` while its status_detail/blocking_dependencies/scope_note still listed multi-T read, full composition tensor, full-tensor writer/readback, and source-backed CI as pending. Applied corrections:
- E390-10B status → `PARTIAL_MULTIDIM_THERMO_ROUNDTRIP_PASS_COMPOSITION_TENSOR_PENDING`
- E390-10B remaining_acceptance populated (was empty [])
- evidence_status.generator_upgrade → `M10_PARTIALLY_IMPLEMENTED_...`
- top-level generator_upgrade.status → `PARTIALLY_IMPLEMENTED_...`

### Node 8 — Representation-only CompOSE re-expression for GPPVA_DD2 (commit 5f7e3c7)
- `GPPVA_DD2_reexpressed`: canonical target (939.56536, 938.27203), r=0.9993982750
- `GPPVA_DD2_matched`: core-matched target (939.56533, 938.272), r=0.9993983069
Both verified physical-invariant (P, mu, f, eps to <1e-15 relative), reversible (Q3/Q7 to <3e-16 abs), readback-confirmed by native Fortran reader. No model regenerated; only dimensionless CompOSE normalization restored. Receipt: `download/eos390_mass_reexpression/E390-08_REEXPRESSION_RECEIPT_20261114.json`.

### Node 9 — Structured Maxwell root census (commit 5f7e3c7)
14-category census over 6 pairs (2 crusts × 3 cores):
- 6 pairs at pinned rtol=1e-12 gate: all MASS_MISMATCH (3e-5 MeV source-decimal gap not silently absorbed).
- 3 pairs (GPPVA_DD2_matched, binary-identical masses): mass-equality gate PASSES, but v1 Maxwell terminates at INTERNAL_PHASE_BOUNDARY because the GPPVA_DD2 crust has an internal Iphase boundary (crust-core splice).
- 09b: explicit precision-policy record (admits the 3e-5 MeV gap under documented source precision + same-family provenance; v1 gate stays strict — honest result).
Finding: the mass-convention discrepancy is resolved by re-expression, but v1's single-smooth-phase requirement blocks the crust-core splice. Receipts: `download/eos390_maxwell_census/E390-09{,_B,_C}_*.json`.

### Node 10 — Maxwell v2 composition-aware orchestration (commit 7b5285a)
New separately-versioned `construct_maxwell_transition_v2` (Python layer, v1 unchanged):
- segments inputs at internal Iphase boundaries; runs v1 on (last-low, first-high) segment pair
- preserves per-segment Iphase/composition (copied, not interpolated across coexistence)
- axis-value slice compatibility (charge_fraction_grid value), not v1's blind iYq integer
- records per-segment local-neutrality evidence (residual, n_rows, composition_available)
- rejects nonzero Nadd (COMPOSITION_UNSUPPORTED)
7 v2 tests pass against live Fortran: `python_frontend/tests/test_maxwell_v2_composition.py`.

### Node 11 — Full test suite at exact HEAD
At HEAD `7b5285a`: **168 passed, 3 skipped, 1 failed** (48.02 s).
- The 1 failure (`test_m7_fail_closed.py::test_python_oracle_and_placeholder_methods_are_blocked`) is PRE-EXISTING at b7989e7 (verified by checking out the pre-audit version of the test file): it fails because a `crust/eos.nb` test fixture is missing in the working tree — an environmental issue unrelated to the mass-convention audit.
- My 16 new tests (9 mass-gate + 7 v2) all pass.
- No regressions introduced by the audit.

### Commit chain (this turn)
```
7b5285a feat(eos390): node 10 Maxwell v2 composition-aware orchestration + tests
8b3a62c fix(eos390-dag): downgrade E390-10B status + reconcile evidence_status (Node 7b)
5f7e3c7 feat(eos390): nodes 8-9 representation-only re-expression + Maxwell root census
e608927 feat(eos390-m5): density-grid refinement study converges (Node 6)
f1f302d feat(eos390): pairwise mass-convention audit + M5 corrections + mass-gate regression tests
b7989e7 (audit start HEAD)
```

### Remaining blockers (unchanged from initial handoff)
- **B2** source composition coverage incomplete (crust per-nucleus composition not available; physical model mass changes not derivable — representation-only re-expression is the only unblocked route, and is implemented for GPPVA_DD2).
- **B3** euler_tolerance=5e-4 derivation pending (v2 receipt records as pending evidence; v1 unchanged).
- **B0** GitHub API/CI token still unavailable in sandbox (PR #13 CI run IDs not retrievable; PR body not updatable). All commits pushed via SSH.
