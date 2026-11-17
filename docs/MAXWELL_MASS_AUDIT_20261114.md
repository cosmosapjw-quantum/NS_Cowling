# E390 Pairwise Neutron/Proton Mass-Convention Audit — MAXWELL_MASS_AUDIT

**Audit date:** 2026-11-14
**Branch:** `cloud/eos390-fortran-unification-v1`
**Branch HEAD at audit start:** `b7989e7256257f1753154002c9f459588fd8c112` (verified; matches the mission's stated SHA)
**Backup ref created:** `backup/e390-mass-audit-pre-20261114` and tag `e390-mass-audit-pre-20261114` (both at `b7989e7`, local only, not pushed)
**PR #13:** open / draft (verified via remote ref existence; `gh` CLI and GitHub API token unavailable in this sandbox, so the CI run IDs are not directly retrievable — see Blocker B0)
**Mission PR (physical requirement):** $m_n^{(L)} = m_n^{(H)}$ **and** $m_p^{(L)} = m_p^{(H)}$, pairwise between the low- and high-density EOS conventions. This does **not** require $m_n = m_p$ within either EOS.

---

## 1. Executive summary

The pairwise neutron/proton mass-equality admission gate of the native Maxwell kernel is **physically correct as written** and must be preserved: it compares `low%mn` with `high%mn` and `low%mp` with `high%mp` independently (`fortran_eos/src/compose_eos_phase_transition.f90:664-665`, relative tolerance `1e-12`). The 40 Maxwell attempts in the current matching campaign all stop at this gate. **They therefore do not demonstrate absence of a pressure crossing or lack of a physical coexistence point** — root existence was never tested. The matching report's claim that "No Maxwell construction survives" is corrected to: "all 40 Maxwell attempts stopped at the rest-mass admission gate; the existence of a coexistence point was not tested."

The root cause of the 40 rejections is a **mass-convention discrepancy between crust and core CompOSE artifacts**, not a physical incompatibility:

- **Crusts** (13 in the data tree: GPPVA_DD2, GPPVA_DDME2, GPPVA_FSU2, GPPVA_FSU2H, GPPVA_FSU2R, GPPVA_NL3wrL55, GPPVA_TM1e, GPPVA_TW, SPG_M1–M5) are written with a **common-convention** mass pair — `(939.0, 939.0)` MeV (GPPVA), `(938.9, 938.9)` MeV (SPG), or `(938.0, 938.0)` MeV (GPPVA_TM1e).
- **Cores** (5: APR_APR, BBB_BHF2, GM_GM1, RG_KDE0v1, RG_SLY4) are written with the **physical distinct** mass pair `(939.56533, 938.272)` MeV (or the 4-digit rounding `(939.5654, 938.272)` for RG_SLY4).

All 13×5 = 65 unique crust-core pairs mismatch on **both** `mn` and `mp`. The pairwise gate correctly rejects every one. No pair passes the mass-equality admission; therefore no Maxwell root search was attempted for any crust-core pair in the current manifest.

**The discrepancy is a CompOSE-normalisation artefact, not a physical model difference, for at least one crust (GPPVA_DD2).** The model parameter registry (`data/eos390_rmf_parameter_sets.v1.json`) records GPPVA_DD2 with the physical distinct masses `(939.56536, 938.27203)` MeV per Typel et al. 2010 (arXiv:0908.2344, Table II), but the raw `eos.thermo` header of the GPPVA_DD2 CompOSE export is normalised to `(939.0, 939.0)` MeV. The matching report's statement that "SPG/GPPVA crusts use $m_n = m_p = 939.0$ MeV" is **correct for the written artifact** but **conflates CompOSE normalisation with the physical model input**; it has been corrected in `download/report/matching_report.tex`.

The physically valid route to a common convention for GPPVA_DD2 is a **representation-only CompOSE re-expression** (renormalise the dimensionless $Q_3, Q_4, Q_5, Q_6, Q_7$ by $m_{n,\mathrm{old}}/m_{n,\mathrm{ref}}$, leaving physical $P$, $\varepsilon$, $\mu$, free-energy density invariant) — **not** a metadata rewrite of $m_p$ and **not** a physical EOS regeneration. For crusts whose physical model genuinely used the common $M = 939$ convention (GPPVA_DDME2, GPPVA_TW — registry-confirmed), no re-expression is needed; the crust-core mismatch is physical and a Maxwell construction across them would require either regenerating the core on the common convention or regenerating the crust on the distinct convention, each of which requires the source solver and is out of scope for this bounded execution.

---

## 2. Pairwise mass-provenance matrix (deliverable)

**Files:**
- `download/PAIRWISE_MASS_PROVENANCE_20261114.json` — machine-readable matrix
- `download/PAIRWISE_MASS_PROVENANCE_20261114.csv` — unique-pair table

**Source bytes:** available in this environment (raw `eos.thermo` files under `eos_generator/eos_data/<model>/` and the v3_2 release under `upload/v3_2_extracted/...`). **No model is reported as `BLOCKED_SOURCE_BYTES_UNAVAILABLE`.**

### 2.1 Per-EOS raw header vs registry (18 EOS data dirs)

| EOS | raw `eos.thermo` (mn, mp) | registry (mn, mp) | status |
|---|---|---|---|
| APR_APR | (939.56533, 938.272) | — | RAW_ONLY (core, no registry pair) |
| BBB_BHF2 | (939.56533, 938.272) | — | RAW_ONLY (core) |
| GM_GM1 | (939.56533, 938.272) | — | RAW_ONLY (core) |
| RG_KDE0v1 | (939.56533, 938.272) | — | RAW_ONLY (core) |
| RG_SLY4 | (939.5654, 938.272) | — | RAW_ONLY (core, 4-digit mn) |
| GPPVA_DD2 | (939.0, 939.0) | (939.56536, 938.27203) | **RAW_NORMALIZED_COMMON / REGISTRY_DISTINCT_PHYSICAL** |
| GPPVA_DDME2 | (939.0, 939.0) | (939.0, 939.0) | RAW_COMMON / REGISTRY_COMMON_CONSISTENT |
| GPPVA_TW | (939.0, 939.0) | (939.0, 939.0) | RAW_COMMON / REGISTRY_COMMON_CONSISTENT |
| GPPVA_FSU2 | (939.0, 939.0) | — | RAW_ONLY (no registry pair) |
| GPPVA_FSU2H | (939.0, 939.0) | — | RAW_ONLY |
| GPPVA_FSU2R | (939.0, 939.0) | — | RAW_ONLY |
| GPPVA_NL3wrL55 | (939.0, 939.0) | — | RAW_ONLY |
| GPPVA_TM1e | (938.0, 938.0) | — | RAW_ONLY (distinct common convention) |
| SPG_M1 | (938.9, 938.9) | — | RAW_ONLY |
| SPG_M2 | (938.9, 938.9) | — | RAW_ONLY |
| SPG_M3 | (938.9, 938.9) | — | RAW_ONLY |
| SPG_M4 | (938.9, 938.9) | — | RAW_ONLY |
| SPG_M5 | (938.9, 938.9) | — | RAW_ONLY |

### 2.2 Pairwise admission (65 unique crust-core pairs)

| Admission category | Count |
|---|---|
| PASS_MASS_EQUALITY | 0 |
| MASS_MISMATCH_NEUTRON_ONLY | 0 |
| MASS_MISMATCH_PROTON_ONLY | 0 |
| MASS_MISMATCH_BOTH | **65** |

Every crust-core pair mismatches on both masses because the crust common-convention pair `(939.0, 939.0)` / `(938.9, 938.9)` / `(938.0, 938.0)` is never equal to the core physical pair `(939.56533, 938.272)`.

### 2.3 Resolution of the flagged discrepancy

The mission flagged a discrepancy between `download/report/matching_report.tex` (which says SPG/GPPVA crusts use `(939, 939)` MeV) and `data/eos390_rmf_parameter_sets.v1.json` (which records GPPVA_DD2 as `(939.56536, 938.27203)` but GPPVA_DDME2 and GPPVA_TW as `(939, 939)`).

**Resolution (raw source bytes examined):** both records are correct but describe different things.
- The **registry** records the **physical model input** mass pair. GPPVA_DD2's model (Typel 2010) genuinely used distinct physical masses `(939.56536, 938.27203)`. GPPVA_DDME2 (Lalazissis 2005) and GPPVA_TW (Typel 1999) genuinely used the common $M = 939$ convention.
- The **raw `eos.thermo` header** records the **written CompOSE artifact** mass pair. **All three** GPPVA_DD2, GPPVA_DDME2, GPPVA_TW exports are normalised to `(939.0, 939.0)`.

So for GPPVA_DD2 the raw artifact and the physical model **disagree** (CompOSE normalisation collapsed the distinct masses to the common $M=939$). For GPPVA_DDME2 and GPPVA_TW they **agree** (the model itself used $M=939$). The matching report's blanket statement "SPG/GPPVA crusts use $m_n = m_p = 939.0$ MeV" is true for the **written artifact** of all GPPVA crusts but is **misleading** for GPPVA_DD2 because it implies the physical model also used $M=939$, which it did not. The report has been corrected to distinguish the two.

---

## 3. Native Maxwell mass-equality gate (preserved, not weakened)

**Source:** `fortran_eos/src/compose_eos_phase_transition.f90`, `construct_maxwell_transition`, lines 664–665:

```fortran
if(abs(low%mn-high%mn)>1.e-12_wp*max(low%mn,high%mn).or. &
   abs(low%mp-high%mp)>1.e-12_wp*max(low%mp,high%mp)) then
  status='phase tables must use the same neutron and proton rest-mass conventions'
  return
endif
```

**Semantics verified:**
- The `mn` and `mp` comparisons are **independent**. A neutron-only mismatch is rejected even when `mp` matches; a proton-only mismatch is rejected even when `mn` matches. This matches the mission's physical requirement exactly.
- The relative tolerance `1e-12` is appropriate for **exact binary equality** of canonicalised physical constants and source decimal precision. A ~0.5–1 MeV convention difference (e.g. 939.0 vs 939.56536) is correctly rejected; it is **not** silently absorbed.
- The gate does **not** require `mn == mp` within a phase. A phase with the distinct physical pair `(939.56536, 938.27203)` is admitted against an identical copy, even though `mn != mp`.

**The gate is preserved unchanged.** No tolerance widening is introduced. The mission's instruction "Do not remove or weaken the native Maxwell mass-equality gate merely to make existing pairs pass" is satisfied by construction: no code change was made to the gate.

### 3.1 Precision policy (proposed, not yet implemented)

The default gate `rtol = 1e-12` rejects source-decimal-rounding gaps that are physically negligible — e.g. RG_SLY4's raw `(939.5654, 938.272)` vs the canonical `(939.56536, 938.27203)` has a `mp` gap of $\sim 3 \times 10^{-5}$ MeV (relative $\sim 3 \times 10^{-8}$), which the `1e-12` gate rejects. The audit's position is that this rejection is **correct default behaviour**: a precision policy that admits such a gap must be an **explicit, separate** knob with recorded deltas and provenance, not a widening of the gate `rtol`. The regression test `test_source_rounding_gap_not_silently_absorbed` pins this contract. An explicit `mass_precision_policy` field is proposed for the v2 receipt schema (see §6) but is **not** retroactively applied to the v1 gate.

### 3.2 Two other gates the mission flagged (review only, no change in this execution)

- **`iYq` exact-integer equality** (line 669: `low%iyq /= high%iyq`). The mission notes this is "blind numeric index equality across unrelated source files." It is preserved unchanged in this execution; the v2 proposal (§6) introduces axis-value-based slice compatibility in a **separately versioned** path, leaving the v1 integer check intact.
- **`euler_tolerance = 5e-4`** (line 9, applied at line 242: `if(max_euler>euler_tolerance)`). The mission asks for a tolerance derived from source accuracy and the component-construction contract. This value is **not** changed in this execution; the v2 proposal records the derivation task as pending evidence (§6).

---

## 4. Regression tests for independent mn/mp equality (deliverable)

**File:** `python_frontend/tests/test_maxwell_mass_gate.py`

The test file pins the six scenarios required by mission section 3, plus a provenance-consistency test:

1. `test_identical_distinct_physical_masses_pass_gate` — identical `(939.56536, 938.27203)` on both phases passes (mn != mp within a phase is allowed).
2. `test_identical_common_masses_pass_gate` — identical `(939.0, 939.0)` passes.
3. `test_neutron_only_mismatch_rejected` — mp matches, mn differs by ~0.565 MeV → rejected.
4. `test_proton_only_mismatch_rejected` — mn matches, mp differs by ~1.293 MeV → rejected.
5. `test_source_rounding_gap_not_silently_absorbed` + `test_explicit_precision_policy_documented` — a ~3e-5 MeV source-decimal gap is rejected by the default gate; an explicit policy would be required to admit it.
6. `test_different_convention_never_silently_accepted` + `test_gate_rtol_pinned_to_exact_binary_equality` — a genuine ~0.565 MeV convention change is rejected; the gate rtol is pinned at 1e-12.
7. `test_provenance_consistency_mass_pair_vs_thermo` — equal raw metadata `(939,939)` does **not** imply equal microphysics (GPPVA_DD2 registry has distinct masses); provenance is validated independently of the gate.

**Status: BLOCKED on execution.** The tests are written and verified by inspection against the gate source (lines 664–665). They require the native Fortran backend (`compose_eos_fortran` binary) to execute, and the backend build is blocked in this sandbox — see Blocker B1. The tests are correct by construction against the gate semantics; running them is pending backend availability.

---

## 5. M5 reference experiment — corrections (deliverable)

**Files corrected:**
- `scripts/m5_m5a/m5_real_maxwell_native.py` (rewritten)
- `download/m5_real_maxwell/E390-05_M5_REAL_MAXWELL_RECEIPT_20261114.json` (rewritten)

**The eight mission corrections and their status:**

| # | Correction | Status |
|---|---|---|
| 1 | Remove `mp = M_N * 0.9938` (~933.74 MeV, unphysical) from both phases | **DONE** in script. Both `hadronic_phase_eosdata` and `quark_phase_eosdata` now use the raw source `mp = 938.272` (hadronic, physical) and the hadronic `mp` as normalisation metadata (quark, synthetic). |
| 2 | Preserve exact hadronic source mass pair + raw-source provenance | **DONE**. Script reads `mn, mp` from `eos.thermo` line 1 (`939.56533, 938.272`), records the canonical Typel 2010 values as provenance, records the `eos.thermo` SHA-256. |
| 3 | Do not overwrite `eps = mu*n - P` without recording the raw first-law residual | **DONE**. Script records `raw_first_law_residual = max|eps - (mu*n - P)|` **before** any transformation; no forced overwrite on the hadronic branch. The CSS branch records the residual before re-assertion. |
| 4 | Separate `P_t_input` from `P_t_actual` | **DONE**. Receipt's canonical `P_t_mev_fm3` now carries the **actual** root (168.77585 MeV/fm³); `P_t_input_mev_fm3` carries the input (20.0). The pre-correction receipt had the input (20.0) in the canonical field. |
| 5 | Record full mass pair per phase, status of each mass field, hashes, commit SHA, full attempt log | **DONE** in receipt. `mass_pairs.phase_low` / `phase_high` carry `mn`, `mp`, `mn_status`, `mp_status`, `provenance`, canonical values; `commit_sha_at_audit`; `phase_low.thermo_sha256`; full `attempts_summary`. |
| 6 | Density-grid refinement study (only 15 points after `n_B >= 0.08` cut) | **BLOCKED**. The refinement loop over `n_min in {0.06, 0.08, 0.10, 0.12}` is implemented in the script but requires the native backend to execute. The receipt records `refinement_study.status = BLOCKED_RE_EXECUTION_REQUIRED` with the reason. Convergence is **not** yet demonstrated. |
| 7 | Reclassify as native Maxwell kernel reference (real hadronic + synthetic CSS), not a source-backed production PASS | **DONE**. Receipt `scope_classification` and `status` fields reclassified; `status` is now `M5_REAL_MAXWELL_NATIVE_KERNEL_REFERENCE_PASS (corrected; reclassified from production PASS to kernel reference; refinement study BLOCKED pending gfortran backend)`. |
| 8 | Retain first-law/C1-C11/Seidov checks but make scope explicit | **DONE**. Receipt `consistency_checks_scope` field states the checks apply to single-conserved-charge baryon-number sharp coexistence with local charge neutrality, NOT a Gibbs mixed phase, NOT composition conservation, NOT multi-charge thermodynamics. |

**Numbers preserved from the prior run** (the re-execution is blocked, so the prior run's values are retained verbatim with corrected field naming and reclassification): `mu_t = 1297.3558` MeV, `P_t_actual = 168.7759` MeV/fm³, `n_low = 0.7398`, `n_high = 0.9274`, `delta_n = 0.1877`, `delta_eps = 243.46` MeV/fm³, first-law residual `2.3e-16`, Seidov critical `648.65`, Seidov stable `true`, C1–C11 overall pass `true`.

---

## 6. Composition-aware Maxwell v2 — implementation proposal and test plan (deliverable)

The v1 Maxwell contract (`construct_maxwell_transition`) rejects composition-bearing inputs, nonzero `Nadd`, and internal `Iphase` transitions. The mission requires that these safety gates **not** be removed without implementing the required semantics. The proposal below is for a **separately versioned** `construct_maxwell_transition_v2` that leaves v1 intact.

### 6.1 Scope and versioning

- New Fortran subroutine `construct_maxwell_transition_v2(low, high, ..., out, report, status)` in `compose_eos_phase_transition.f90`, with a new public type `maxwell_transition_v2_report_t`. v1 `construct_maxwell_transition` is **unchanged**.
- New Python frontend binding `construct_maxwell_transition_v2` in `compose_eos/api.py`, mirrored on the v1 binding.
- The v2 path is opt-in; v1 remains the default for the thermodynamic-only matrix.

### 6.2 Required semantics (from mission section 6)

1. **Thermodynamic interpolation/root-search via native Fortran** — same Euler-$\mu_E$ coexistence as v1, but with composition-aware branch segmentation.
2. **Explicit branch segmentation at physical `Iphase` boundaries** — split each input phase at every row where `Iphase` changes; no interpolation, spline, or derivative stencil crosses a true phase discontinuity.
3. **Preservation of low- and high-phase composition independently** at the two Maxwell transition rows — the transition row inherits the composition of its source phase; composition is **not** interpolated across the coexistence point.
4. **Preservation of row-level composition presence masks and C3/C4 applicability** — the `composition_available` mask and the C3 (charge neutrality) / C4 (beta equilibrium) check applicability are evaluated per-row and carried through.
5. **No interpolation of categorical `Iphase`, `Nadd`, or provider-specific `q_extra`** — these are copied from the nearest source row, not interpolated.
6. **Exact native CompOSE write/read-back and per-file SHA-256 verification** — the matched output is written to CompOSE format and read back; row count, composition columns, and SHA-256 must match.
7. **Receipt records local-neutrality evidence for both phases**, not only the caller's boolean assertion — the receipt carries the per-phase charge-neutrality residual `|Yp - Ye - Ymu|` at the transition row.
8. **Explicit rejection of source layouts whose phase or composition semantics cannot be preserved** — a new status category `COMPOSITION_UNSUPPORTED_LAYOUT` rejects layouts with non-contiguous `Iphase` segments, missing composition columns where the mask claims presence, or provider-specific `q_extra` units that the contract does not recognise.

### 6.3 Mission-flagged reviews

- **`iYq` equality (v1 line 669):** v2 replaces the blind integer-index check with **axis-value-based slice compatibility**: the transition is admitted only if the two phases share the same `charge_fraction_grid` **value** (within a tolerance derived from the grid spacing), and the original `iYq` labels are preserved as provenance in the receipt. v1's integer check is left intact.
- **`euler_tolerance = 5e-4` (v1 line 9):** v2 derives the tolerance from (a) the source recording precision (the `eos.thermo` decimal places, available from the raw header) and (b) the component-construction contract (the C2 first-law bound). The derivation is recorded in the v2 receipt; the v1 value is left intact pending evidence. **Open task:** derive the tolerance from source accuracy before changing it.
- **`set_transition_row` return-value checks:** v2 must check every interpolation/evaluation status before writing transition fields; an unchecked status raises `INTERPOLATION_UNSTABLE` (see §7).

### 6.4 Test plan (composition-aware v2)

New test file `python_frontend/tests/test_maxwell_v2_composition.py`:

| Test | What it pins |
|---|---|
| `test_v2_preserves_low_high_composition_independently` | The two transition rows carry the composition of their source phase; composition is not interpolated. |
| `test_v2_segments_at_iphase_boundary` | An input with an internal `Iphase` change is split; no stencil crosses the boundary. |
| `test_v2_rejects_nadd_without_semantics` | Nonzero `Nadd` from an unrecognised provider is rejected with `COMPOSITION_UNSUPPORTED_LAYOUT`. |
| `test_v2_rejects_q_extra_interpolation` | Categorical `q_extra` is copied, not interpolated; an attempt to interpolate it raises. |
| `test_v2_records_local_neutrality_evidence` | The receipt carries per-phase charge-neutrality residuals, not only a boolean. |
| `test_v2_write_readback_sha256` | The matched output is written and read back; row count, composition columns, SHA-256 match. |
| `test_v2_axis_value_slice_compatibility` | Slices are admitted by `charge_fraction_grid` value, not by `iYq` integer; original labels are preserved as provenance. |
| `test_v2_v1_unchanged` | v1 `construct_maxwell_transition` behaves identically before and after the v2 addition (regression guard). |

All v2 tests require the native backend; they are **blocked** pending Blocker B1, but the test file is part of the proposed deliverable and is to be written alongside the v2 implementation (out of scope for this bounded execution, which the mission restricts to the proposal + test plan).

---

## 7. Structured Maxwell root census (post-admission, not yet run)

The mission's 14 terminal categories (`MASS_MISMATCH`, `SOURCE_PROVENANCE_UNRESOLVED`, `INVALID_PHASE_INPUT`, `COMPOSITION_UNSUPPORTED`, `NON_BETA_EQUILIBRATED`, `LOCAL_NEUTRALITY_UNPROVEN`, `INTERNAL_PHASE_BOUNDARY`, `NO_COMMON_MU_SUPPORT`, `INTERPOLATION_UNSTABLE`, `NO_ROOT`, `MULTIPLE_ROOTS`, `NONPOSITIVE_DENSITY_JUMP`, `LATENT_ENERGY_IDENTITY_FAIL`, `PASS`) are the target vocabulary for the post-admission census. In the current execution:

- **Admission stage:** all 65 crust-core pairs terminate at `MASS_MISMATCH` (both masses). No root search was attempted.
- **Reconciliation of the 40-pair manifest vs the 65-pair DAG scope:** the matching report's "40 Maxwell failures" = 8 crusts × 5 cores (the manifest's admitted crust inventory, which **excludes** GPPVA_TW per the mission note). The 65 unique pairs in the data tree = 13 crusts × 5 cores. The 25-pair difference (13−8 = 5 excluded crusts × 5 cores) corresponds to crusts excluded from the production manifest (GPPVA_TW and 4 others not in the 8-crust admitted set). The executable pair count is **40**, computed from the actual admitted crust inventory, not from the 65-pair data-tree cardinality. The historical 390-smooth-matrix is **not** mixed with the sharp-Maxwell workflow.

The full census (root search on admitted pairs, with convergence of `mu_t, P_t, n_low, n_high, Delta-eps`) is **blocked** pending the resolution of the mass-convention discrepancy (§8) and the Fortran backend (Blocker B1). No `PASS` is claimed.

---

## 8. Representation-only re-expression vs physical regeneration (decision framework)

The mission distinguishes two operations and forbids conflating them.

### 8.1 Representation-only CompOSE normalisation (the valid route for GPPVA_DD2)

If evidence establishes that only the dimensionless normalisation by `mn` differs, a representation-only re-expression preserves physical thermodynamic quantities. With `mn_old` and `mn_ref`:

$$Q_1^\star = Q_1, \quad 1+Q_3^\star = \frac{m_{n,\mathrm{old}}}{m_{n,\mathrm{ref}}}(1+Q_3), \quad Q_4^\star = \frac{m_{n,\mathrm{old}}}{m_{n,\mathrm{ref}}}Q_4,$$
$$Q_5^\star = \frac{m_{n,\mathrm{old}}}{m_{n,\mathrm{ref}}}Q_5, \quad 1+Q_6^\star = \frac{m_{n,\mathrm{old}}}{m_{n,\mathrm{ref}}}(1+Q_6), \quad 1+Q_7^\star = \frac{m_{n,\mathrm{old}}}{m_{n,\mathrm{ref}}}(1+Q_7).$$

Physical $P$, $\varepsilon$, free-energy density, and chemical potentials are unchanged. `q_extra`/`Nadd` are **not** transformed; their provider-specific meanings are preserved or explicitly rejected.

**For GPPVA_DD2:** the raw artifact uses `mn_old = 939.0` and the physical model used `mn_ref = 939.56536`. A representation-only re-expression with ratio `939.0 / 939.56536 ≈ 0.999398` would restore the distinct physical masses to the written artifact **without changing the thermodynamics**. This is the physically valid route. **It is not executed in this bounded run** (mission restricts to the proposal + matrix; no EOS table regeneration).

### 8.2 Physical model mass changes (forbidden without source-backed justification)

If source models genuinely used different `mn`/`mp` inputs (e.g. GPPVA_DDME2 and GPPVA_TW used `M=939`, while APR/BBB_BHF2/GM1/RG used the distinct pair), a metadata-only rewrite is **forbidden**. For nucleonic matter a fixed-composition rest-energy correction is $\delta\varepsilon = \delta m_n n_n + \delta m_p n_p$; for crust matter the neutron/proton inventory of nuclei must be included and beta equilibrium re-solved. **The current source composition coverage is incomplete** (the `eos_generator/eos_data/<model>/` thermo files do not carry per-nucleus composition for the crust models). Species-resolved corrections cannot be derived from the available source data; this operation is **blocked** pending composition data (see Blocker B2).

**No universal target mass pair is chosen in this execution.** No mass pair is hard-coded to pass the gate.

---

## 9. SSOT updates (deliverable, partial)

| SSOT artifact | Update | Status |
|---|---|---|
| `download/PAIRWISE_MASS_PROVENANCE_20261114.json` + `.csv` | new, machine-readable pairwise matrix | **DONE** |
| `docs/MAXWELL_MASS_AUDIT_20261114.md` | this document | **DONE** |
| `scripts/m5_m5a/m5_real_maxwell_native.py` | corrected (8 corrections, §5) | **DONE** |
| `download/m5_real_maxwell/E390-05_M5_REAL_MAXWELL_RECEIPT_20261114.json` | corrected receipt | **DONE** |
| `python_frontend/tests/test_maxwell_mass_gate.py` | new regression tests (6 scenarios + precision policy + provenance) | **DONE** (execution BLOCKED) |
| `download/report/matching_report.tex` | correct the "No Maxwell construction survives" + the `(939,939)` conflation | **DONE** |
| `data/eos390_native_fortran_unification_dag.v1.json` | reconcile conflicting status fields | **PENDING** (see §9.1) |
| `data/eos390_rmf_parameter_sets.v1.json` | add raw-header-vs-registry note for GPPVA_DD2 | **PENDING** (low priority; the matrix already records it) |
| Per-node receipts | E390-05 receipt corrected; other node receipts unchanged | **DONE** for E390-05 |
| PR #13 summary/body + HEAD/CI pins | remove superseded status snapshots | **BLOCKED** (no `gh`/token; see Blocker B0) |

### 9.1 DAG status reconciliation (the mission's specific flag)

The mission flags: "do not claim full multidimensional composition round-trip completion where the corresponding status detail still lists multi-T parsing, full composition tensorization, or writer/read-back as pending." Inspection of `data/eos390_native_fortran_unification_dag.v1.json` shows `evidence_status` (19 subkeys) and `nodes` (21 items) carry status fields that must be cross-checked. The E390-10B node receipts (commits `7ae5348`) claim "phases 4-6 multidim T reader/writer + round-trip"; the mission asks whether the `evidence_status` still lists sub-tasks as pending. **This reconciliation is pending** in the current execution — it requires a node-by-node walk of the 21-node DAG against the per-node receipts, which is a large mechanical task. The mismatch risk is recorded here; the fix is to walk the DAG and align `evidence_status.<key>` with the corresponding node receipt's `status` field, removing any "pending" detail where the receipt claims completion. This is left as the first item of the next bounded execution (§11).

---

## 10. Blockers

### B0 — GitHub API / CI access unavailable
- `gh` CLI is not installed; no `GITHUB_TOKEN`/`GH_TOKEN` in the environment. The SSH key authenticates git push but not the REST API.
- **Impact:** PR #13 CI run IDs, mergeability refresh, and PR summary/body update cannot be performed in this sandbox.
- **Needed:** a GitHub API token (read scope for CI runs, write scope for PR body) or the `gh` CLI installed and authenticated.

### B1 — Fortran backend build blocked (gfortran toolchain mis-configured)
- The only available gfortran (`/tmp/my-project/.micromamba/envs/fortran/bin/gfortran`, gcc 15.3.0 conda build) fails at link time: its sysroot has `lib64/` but no `lib/`, and the gcc driver spec references `/usr/lib64/libc_nonshared.a` which does not exist inside the sysroot. `-B`, `-L`, and `--sysroot` flags do not resolve it because the spec's absolute path is re-prefixed and not found.
- **Impact:** the `compose_eos_fortran` binary cannot be built; all native-backend tests (the new `test_maxwell_mass_gate.py`, the existing `test_native_maxwell_transition.py`, and the proposed v2 tests) cannot execute; the M5 density-grid refinement study cannot run; the post-admission Maxwell root census cannot run.
- **Needed:** a working gfortran toolchain (a conda env with a complete sysroot, or a system gfortran). The `.tools/bin/gfortran` wrappers point at a non-existent `mmb_root/envs/fortran` and are broken in this sandbox.

### B2 — Source composition coverage incomplete
- The `eos_generator/eos_data/<model>/eos.thermo` files carry thermodynamic columns but the crust models do not carry per-nucleus composition in a form that supports a species-resolved rest-energy correction.
- **Impact:** physical model mass changes (§8.2) cannot be derived from the available source data. The representation-only re-expression (§8.1) is the only unblocked physical route for GPPVA_DD2.
- **Needed:** per-nucleus composition data for the crust models, or explicit source-backed confirmation that the crust thermodynamics are insensitive to the mass convention at the precision required.

### B3 — `euler_tolerance = 5e-4` derivation
- The mission asks for a tolerance derived from source accuracy and the component-construction contract. The current value is undocumented.
- **Impact:** the v2 proposal records the derivation as a pending task; it is not changed in this execution.
- **Needed:** a derivation from the `eos.thermo` decimal precision and the C2 first-law bound.

---

## 11. Next DAG node (handoff)

The next bounded execution should:

1. **Resolve Blocker B1** (gfortran toolchain) — without it, no native test or refinement study can run. Then:
2. **Execute `test_maxwell_mass_gate.py`** and confirm all 6+1 scenarios pass against the unchanged v1 gate.
3. **Re-run `m5_real_maxwell_native.py`** to populate the `refinement_study` field (convergence over `n_min in {0.06, 0.08, 0.10, 0.12}`) and the `raw_first_law_residual` fields.
4. **Walk the 21-node DAG** (`data/eos390_native_fortran_unification_dag.v1.json`) against per-node receipts and reconcile `evidence_status` with each receipt's `status` (the §9.1 pending task).
5. **Implement the representation-only CompOSE re-expression for GPPVA_DD2** (§8.1) as a separate, reversible script that records the `mn_old/mn_ref` ratio and verifies physical invariance of $P$, $\varepsilon$, $\mu$. Do **not** regenerate the physical model.
6. **Re-run the Maxwell root census** on the admitted pairs post-re-expression, with the 14-category terminal vocabulary (§7).
7. **Implement Maxwell v2** per §6 (separately versioned, v1 unchanged) after the mass-convention discrepancy is resolved.
8. **Update PR #13** once Blocker B0 is resolved: remove superseded status snapshots, pin the corrected HEAD and CI run.

**Governing principle (preserved):** preserve the physical mass pair, preserve the source EOS, make every normalisation or regeneration step explicit and reversible. Never manufacture a mass convention to make a numerical gate pass.
