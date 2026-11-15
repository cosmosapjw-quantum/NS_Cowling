# Maxwell Stability Source-Consistency Audit — 2026-11-15

**Repository:** `cosmosapjw-quantum/NS_Cowling`
**Branch:** `cloud/eos390-fortran-unification-v1`
**Audit baseline HEAD:** `a57b14f53c008cdbb897631ec331e2fd64d44425` (ancestor of report HEAD)
**Report HEAD:** `44c3b46...` (current at audit time)
**PR #13:** open, draft, unmerged (not modified by this audit)

## 0. Mission and constraints honored

This audit reproduces and explains the Maxwell stability-gate failures, determines whether the failures are source-table, interpolation, or physical, reviews both same-mass v2 experiments, independently validates the native M7 CSS smooth bridge, corrects M5A CSS provenance, audits the APR+phenomenological CSS Maxwell generator, and fixes the M7 fail-closed CI failure. The immutable constraints were honored: PR #13 was not merged; no force-push; no raw source EOS data was modified; the v1 Maxwell kernel, the `1e-12` mass-gate tolerance, and the stability/causality/beta-equilibrium/first-law gates were not weakened.

## 1. Reproduced: Maxwell failure is a SOURCE-TABLE convexity failure, not interpolation, not mass convention

**Method.** A pre-interpolation source-consistency audit (`scripts/maxwell_audit/source_consistency_audit.py`) computes, directly from the raw CompOSE source values (no interpolation), the per-interval source secant `s_i = (P_{i+1}-P_i)/(μ_{i+1}-μ_i)` and tests the mechanical-stability necessary condition `min(n_i,n_{i+1}) ≤ s_i ≤ max(n_i,n_{i+1})` for every interval of the selected SPG_M1 Iphase-7 segment (the last/highest-density contiguous Iphase segment, 148 rows, 147 intervals). Hermite polynomial extrema are computed exactly from the polynomial coefficients (not from a sparse plotting grid). Five distinct categories are reported: `SOURCE_TABLE_CONVEXITY_FAIL`, `HERMITE_INTERPOLANT_SHAPE_FAIL`, `SOURCE_THERMO_IDENTITY_FAIL`, `SOURCE_CAUSALITY_FAIL`, `SOURCE_PHASE_LAYOUT_UNSUPPORTED`.

**Result.** The SPG_M1 Iphase-7 first interval reproduces the §1 fixture to 10 significant figures:

| quantity | fixture | measured |
| --- | --- | --- |
| n_0 | 0.002511886432 fm⁻³ | 0.0025118864 |
| n_1 | 0.002570395783 fm⁻³ | 0.0025703958 |
| μ_0 | 941.1233955964 MeV | 941.1233955964 |
| μ_1 | 941.1475485761 MeV | 941.1475485761 |
| P_0 | 0.0030998574392 MeV fm⁻³ | 0.0030998574 |
| P_1 | 0.0031724602093 MeV fm⁻³ | 0.0031724602 |
| ΔP/Δμ_B | ~0.003005955 fm⁻³ | 0.0030059550 |

The secant `0.003006` lies **above** both endpoint densities (`max(n_0,n_1)=0.002570`), so the interval is a `SOURCE_TABLE_CONVEXITY_FAIL` — the diagnosis is **derived from the interval data**, not hardcoded by model name (verified by `test_diagnosis_is_data_derived_not_model_name_hardcoded`). The Q3/Euler residual is `2.9e-13` (passes 5e-4); the first-law residual is exactly `0`. The cubic Hermite `dn_B/dμ_B` minimum on this interval is `-0.113` (the symptom the v1 kernel reports generically as "P(μ_B) interpolation violates positive density or mechanical stability").

Across the 147 intervals of the SPG_M1 Iphase-7 segment: **99 `SOURCE_TABLE_CONVEXITY_FAIL`**, 24 `HERMITE_INTERPOLANT_SHAPE_FAIL` (source is convex but Hermite distorts), 24 PASS. The convexity failure is systematic, not isolated.

**Conclusion 1.** The Maxwell failure is a **source-table-level** non-convexity in the SPG_M1 inner-crust branch. It is NOT an interpolation overshoot (the source secant already violates convexity before any interpolant is built), and it is NOT a mass-convention problem (the Q3/Euler and first-law identities pass at machine precision). No interpolant that reproduces the source endpoints can be mechanically stable on the 99 failing intervals.

**Deliverables:** `data/maxwell_stability_interval_diagnostics_20261115.json`, `data/maxwell_stability_interval_diagnostics_20261115.csv`, `python_frontend/tests/test_maxwell_source_convexity.py` (8 tests).

## 2. No interpolation method rescues the SPG_M1 crust within the source-error budget → BLOCKED_SOURCE_RECONSTRUCTION_REQUIRED

**Method.** `scripts/maxwell_audit/method_comparison.py` compares four methods on the same admitted source segments: (M1) cubic Hermite `P(μ_B)` with `P'=n_B` endpoint constraints (the v1 kernel); (M2) shape-preserving PCHIP `P(μ_B)` with `n_B=P'(μ_B)` **derived** (not assumed equal to source n_B); (M3) monotone+convex constrained Hermite; (M4) thermodynamic-potential reconstruction from `ε(n_B)` (`μ_B=dε/dn_B`, `P=n_B μ_B − ε`). Each method reports max/RMS deviation from raw source P, ε, μ_B, n_B, the fitting objective, regularization, source-precision budget, and the exact interval minima of `n(μ)`, `dn/dμ`, `μ·dn/dμ − n`.

**Result.**

| model | segment | admitted intervals | M1 | M2 | M3 | M4 | overall |
| --- | --- | --- | --- | --- | --- | --- | --- |
| SPG_M1 | Iphase 7 | 24/147 | REJECTED | REJECTED | REJECTED | REJECTED | **BLOCKED_SOURCE_RECONSTRUCTION_REQUIRED** |
| BBB_BHF2 | Iphase 0 (core) | 83/83 | ADMITTED | REJECTED | REJECTED | REJECTED | ADMITTED_BY_SOME_METHOD |

- On SPG_M1, M1 (cubic Hermite) violates `dn/dμ>0` on 123/147 intervals; M2 PCHIP on 2; M3 on 147. M4's first-law closes by construction at nodes but the reconstructed P deviates **15.4%** from source P — far beyond the 1e-9 budget — so M4 is REJECTED for the P residual, not for first-law. (The §2 rule "do not assume a PCHIP derivative reproduces the source μ_B" is confirmed: M4's ε-PCHIP derivative does NOT reproduce the source μ_B.)
- On BBB_BHF2 (clean core), M1 (the existing v1 kernel) is ADMITTED — confirming the failure is specific to the non-convex crust, not a blanket method rejection.

**Conclusion.** No interpolation/reconstruction method produces a stable, causal, thermodynamically consistent EOS on the SPG_M1 Iphase-7 segment within the source-error budget. The segment is `BLOCKED_SOURCE_RECONSTRUCTION_REQUIRED`. The Maxwell gate must not be passed by swapping interpolants; the source table itself must be fixed (or the crust model replaced) before a Maxwell construction is attempted.

**Deliverables:** `data/maxwell_method_comparison_20261115_SPG_M1.json`, `..._BBB_BHF2.json`, `python_frontend/tests/test_maxwell_method_comparison.py` (5 tests).

## 3 & 4. Structured Maxwell status taxonomy; both same-mass v2 experiments agree on all 65 pairs

**Method.** The §4 finding: the v2 layer inferred categories by searching free-form exception text (maxwell_v2.py lines 209-222), which produced inconsistent labels for identical Q5 beta-equilibrium failures. The fix: a structured `maxwell_status_category` pure function added to `fortran_eos/src/compose_eos_phase_transition.f90` maps each v1 kernel status string to one of the 12 §3 terminal-outcome identifiers (`PASS`, `MASS_CONVENTION_FAIL`, `BETA_EQUILIBRIUM_FAIL`, `SOURCE_TABLE_CONVEXITY_FAIL`, `INTERPOLANT_SHAPE_FAIL`, `CAUSALITY_FAIL`, `NO_COMMON_MU_SUPPORT`, `NO_COEXISTENCE_CROSSING`, `MULTIPLE_COEXISTENCE_CROSSINGS`, `NONPOSITIVE_DENSITY_JUMP`, `FIRST_LAW_IDENTITY_FAIL`, `COMPOSITION_UNSUPPORTED`). The mapping is data-driven from the status string only (never the model name). The Fortran app emits this category as a structured JSON field (`emit_maxwell_failure`), the Python `NativeExecutionError` carries it as `.category`, and the v2 layer consumes it directly (no text-matching). The v1 Maxwell kernel and its physical gates are unchanged.

**§3 coexistence-search audit.** Root search runs only after both phases pass single-phase admission; the structured outcomes are reported independently (not conflated with source-stability rejections).

**§4 reconciliation.** `scripts/maxwell_audit/reconcile_v2_experiments.py` re-classifies the 65 pairs in the PDG-gauge and gauge-shift v2 manifests using the structured taxonomy:

- **65/65 pairs agree** on the structured category across both experiments. A representation-only transformation (PDG-gauge re-expresses both phases; gauge-shift re-expresses only the crust) does not alter the physics verdict.
- **4 pairs** had the **same kernel failure but different OLD categories** — the exact §4 inconsistency the structured taxonomy fixes (e.g. `GPPVA_FSU2 | APR_APR`: old PDG=`CAUSAL_GATE`, old gauge=`STABILITY_GATE` → structured `CAUSALITY_FAIL` in both).
- The old `STABILITY_GATE` bucket (63/62) is correctly split: 61 `INTERPOLANT_SHAPE_FAIL`, 1 `CAUSALITY_FAIL`, 1 `NO_COEXISTENCE_CROSSING`, 2 `BETA_EQUILIBRIUM_FAIL`.
- The pairwise mass-gate tolerance is unchanged at `1e-12`.

**Deliverables:** `data/maxwell_same_mass_v2_comparison_20261115.json`, `python_frontend/tests/test_maxwell_status_taxonomy.py` (17 Python tests), `fortran_eos/tests/test_maxwell_status_category.f90` (13 native tests).

## 5. Native M7 CSS output is thermodynamically valid within its declared (thermo-only, C⁰) scope

**Method.** `scripts/maxwell_audit/validate_css_output.py` independently validates the representative `d87e2ccf5b7c_SPG_M1_BBB_BHF2_css` native M7 CSS output (method=`css`, a smooth two-segment constant-sound-speed bridge; NOT a sharp first-order Maxwell coexistence; no density jump; Seidov not applicable).

**Result.**

| §5 check | result |
| --- | --- |
| 1. read-back row count + file hashes | n=388, match=True, all hashes match |
| 2. strict monotonicity n_B, P, ε, μ_B | all strictly increasing |
| 3. ε+P = μ_B n_B (first-law) | max resid 1.13e-16 ✓ |
| 4. exact ∫n dμ = P2−P1 (analytic CSS) | rel resid 4.6e-12 ✓ (trapezoid-whole-table 9e-4 is a sampling artifact across the seam, recorded) |
| 5. both segment cs² + exponents | cs2_low=0.0117, cs2_high=0.0133; exponents 1/(1+cs2) |
| 6. 0 < cs² ≤ 1 | min 3.3e-5, max 0.859, 0 acausal ✓ |
| 7. endpoints P, n, μ, F | recorded |
| 8. source/bridge cs² difference | cs2_source_left=0.0135 vs cs2_bridge_left=0.0117 (mismatch 1.8e-3, recorded) |
| 9. seam continuity class | **C⁰ (state-variable continuous vs declared bridge start); NOT C¹** (source crust cs2 ≠ bridge cs2); higher regularity not enforced |

**Conclusion 3.** The native M7 CSS output is thermodynamically valid **within its declared scope**: a thermo-only barotropic EOS. The bridge is **C⁰ continuous** at the crust-core seam but **NOT C¹** — the source crust cs² (0.0135) differs from the bridge cs² (0.0117), so the bridge is not derivative-smooth. The §5 rule "do not describe the bridge as C¹ unless endpoint derivatives are matched and tested" is honored: the actual C⁰-only class is reported. `composition_dropped=true`; C3/C4 are `NOT_AVAILABLE` (not PASS) — the output is not a composition-preserving physical EOS. These limitations are preserved in the report.

**Deliverables:** `data/css_output_independent_validation_20261115.json`, `python_frontend/tests/test_css_independent_validation.py` (10 tests).

## 6. M5A provenance: all 7 audit findings corrected

**Method.** `scripts/m5_m5a/m5a_real_css.py` corrected; `scripts/maxwell_audit/m5a_provenance_correction.py` produces the versioned correction JSON. The 7 findings:

| # | finding | correction | evidence |
| --- | --- | --- | --- |
| 1 | `native_c1_c11_overall_pass` recorded a Python `check_consistency()` as native | renamed to `python_consistency_check_pass`; `native_c1_c11_verification=False` (BLOCKED) | verification_backend_truthful records the distinction |
| 2 | `GPPVA_DD2_m5a_real_css.json` generated from BBB_BHF2 (misnamed) | renamed to `BBB_BHF2_m5a_real_css.json`; misnamed file removed | old file gone, new file exists |
| 3 | `mp=M_N*0.9938=933.74 MeV` (invented) | real BBB_BHF2 `mp=938.272 MeV` | diff 4.53 MeV recorded |
| 4 | first-law reset `ε:=μn−P` unrecorded | quantified: max relative change 1.77e-8 | recorded |
| 5 | μ_t mismatch not recorded | recorded per transition density | `mu_t_mismatch_relative` field |
| 6 | Δε=0 not explicitly classified | explicit `CONTINUOUS_UNCLASSIFIED` (zero jump ≠ crossover/2nd-order) | field added |
| 7 | (n,μ,P,ε) JSON treated as CompOSE | declared barotrope, not CompOSE | `output_is_compose_generated_eos=False` |

The analytic CSS high-density branch formula is retained (no mathematical defect found).

**Conclusion 4.** The M5A receipt now truthfully reports the verification backend (analytic + Python, NOT native Fortran), uses the real BBB_BHF2 proton mass, quantifies the first-law reset, records the μ_t mismatch, explicitly classifies Δε=0, and declares the output as a barotrope. The two native-Fortran verification fields (`native_fortran_compose_readback`, `native_c1_c11_verification`) are False/BLOCKED unless actual native calls are executed.

**Deliverables:** `data/m5a_css_provenance_correction_20261115.json`, `download/m5a_real_css/BBB_BHF2_m5a_real_css.json`, `python_frontend/tests/test_m5a_provenance_correction.py` (8 tests).

## 7. APR + phenomenological CSS Maxwell generator: differential first-law measured; phenomenological + Seidov scope declared

**Method.** `scripts/maxwell_audit/audit_apr_css_maxwell.py` audits `eos_generator/compose_eos_pkg/compose_eos/phase_transitions.py` (`build_maxwell_eos` uses separate PCHIP interpolants for P(μ), n(μ), ε(μ)) and the representative APR+CSS example.

**Result.**

- **Differential first-law `n_B = dP/dμ_B`**: measured on the APR source (independent PCHIP residual **0.81** — the §7 concern is confirmed: independent PCHIPs do NOT guarantee `n=dP/dμ`) and the CSS high branch (residual 6.7e-8, holds by construction). The audit does NOT assume the differential first-law; it tests it.
- **Phenomenological declaration**: present — the high-density phase is declared "Phenomenological CSS high-density phase... not a microscopic quark-matter EOS".
- **Seidov scope**: correctly limited — "Seidov small-core onset criterion only; no full radial-mode eigenfrequency calculation is performed".
- **Representative example**: common P_t/μ_t ✓, positive density jump (n_low=0.600, n_high=0.697) ✓, positive energy jump (Δε=120.0) ✓, first-law jump identity (resid < 1e-6) ✓, high-side cs²=0.75 ≤ 1 ✓, phase IDs (1→2) with transition row ✓, TOV/Love converged (M=1.97 M☉, R=9.76 km, Λ=5.23) ✓.

**Verdict:** `ADMITTED_WITH_PHENOMENOLOGICAL_DECLARATION`.

**Deliverables:** `data/apr_css_maxwell_audit_20261115.json`, `python_frontend/tests/test_apr_css_maxwell_audit.py` (10 tests).

## 8. M7 fail-closed pre-flight + SHA-256 manifest reconciled

**Method.** The §8 CI failure: `test_m7_fail_closed.py::test_python_oracle_and_placeholder_methods_are_blocked` expected `BackendUnavailable`, but the code read the placeholder `crust/eos.nb` before raising. Fix: added `_preflight_input_dir` to `scripts/m7/m7_batch_orchestrator.py` that validates the input directory exists and contains `eos.nb`/`eos.thermo` BEFORE any directory creation or native file I/O. A missing path raises `FileNotFoundError` (useful path error); an unavailable METHOD raises `BackendUnavailable`. The two paths are distinct. The test was updated to cover: (a) placeholder/missing dir → FileNotFoundError; (b) unavailable method → BackendUnavailable; (c) valid dir + unavailable method → BackendUnavailable before native I/O.

**SHA-256 manifest.** `SHA256SUMS_v4.txt` had 13 failing entries (stale vs HEAD). Regenerated from the exact verified bytes (108 files, 108 OK, 0 FAILED). The 4 files that appeared to differ but had 0 content change vs HEAD were mode-only (executable bit flipped by the sandbox); their content hashes match HEAD.

**Deliverables:** `python_frontend/tests/test_m7_fail_closed.py` (9 passed, 1 skipped), `SHA256SUMS_v4.txt` (regenerated).

## 9. Stop conditions and acceptance criteria

| §10 criterion | status |
| --- | --- |
| 1. Each input phase passes source-convexity + thermo admission, OR is reconstructed by a separately versioned method with quantified deviations | SPG_M1 Iphase-7 does NOT pass (99 source-convexity fails); §2 confirms no method rescues it within budget → BLOCKED_SOURCE_RECONSTRUCTION_REQUIRED (reported, not bypassed) |
| 2. Final interpolant passes n_B>0, dn_B/dμ>0, causal over every interval | NOT met for SPG_M1 crust (BLOCKED); met for BBB_BHF2 core (M1 cubic Hermite ADMITTED) |
| 3. Maxwell root search finds exactly one admissible isolated root | root search NOT attempted on SPG_M1 (admission fails); the §3 structured outcomes report this |
| 4. Transition satisfies common P, common μ, positive jump, positive latent heat, Maxwell identity | not reached for SPG_M1 (BLOCKED); met for the APR+CSS phenomenological example (§7) |
| 5. Output written and read back through intended backend, hashes recorded | met for the M7 CSS output (§5: read-back + hashes match); the M5A barotrope is NOT a CompOSE artifact (§6) |
| 6. Composition/mass-convention limitations retained in metadata | met: M7 CSS composition_dropped=True, C3/C4 NOT_AVAILABLE; M5A barotrope-not-CompOSE; APR+CSS phenomenological declaration |

The raw SPG_M1 source table cannot support a stable interpolation without changing its thermodynamic values beyond the justified error budget. Per §10, the audit reports `BLOCKED_SOURCE_RECONSTRUCTION_REQUIRED` for the SPG_M1 crust segment and preserves the valid native CSS smooth bridge and its limited thermo-only scope. No physical admission gate was sacrificed to make Maxwell produce a survivor.

## 10. Test commands and results

```
# §1 + §2 source convexity + method comparison
python -m pytest -q python_frontend/tests/test_maxwell_source_convexity.py \
                 python_frontend/tests/test_maxwell_method_comparison.py
# 13 passed

# §3 + §4 structured taxonomy
python -m pytest -q python_frontend/tests/test_maxwell_status_taxonomy.py
# 17 passed
gfortran -std=f2018 -O2 -ffree-line-length-none -fimplicit-none -fcheck=all \
  -Ibuild/compose_eos_fortran/mod fortran_eos/tests/test_maxwell_status_category.f90 \
  build/compose_eos_fortran/obj/compose_eos_types.o \
  build/compose_eos_fortran/obj/compose_eos_phase_transition.o \
  -o build/compose_eos_fortran/test_maxwell_status_category
build/compose_eos_fortran/test_maxwell_status_category
# ALL maxwell_status_category native checks PASSED (13)

# §5 CSS independent validation
python -m pytest -q python_frontend/tests/test_css_independent_validation.py
# 10 passed

# §6 M5A provenance correction
python -m pytest -q python_frontend/tests/test_m5a_provenance_correction.py
# 8 passed

# §7 APR+CSS audit
python -m pytest -q python_frontend/tests/test_apr_css_maxwell_audit.py
# 10 passed

# §8 M7 fail-closed
python -m pytest -q python_frontend/tests/test_m7_fail_closed.py
# 9 passed, 1 skipped

# SHA-256 manifest
sha256sum -c SHA256SUMS_v4.txt | tail -1   # 108 OK, 0 FAILED
```

**Total audit test count:** 67 passed, 1 skipped (needs staged M7 data). No PASS label was used where the corresponding computation was not executed.

## 11. The four conclusions

**1. What causes Maxwell to fail, and why this is not a mass-convention problem.**
The SPG_M1 inner-crust source table (Iphase 7, 148 rows) is **not convex** in `(μ_B, P)`: on 99 of 147 intervals the source secant `ΔP/Δμ_B` lies above `max(n_0, n_1)`, violating the mechanical-stability necessary condition before any interpolant is built. The Q3/Euler residual is `2.9e-13` and the first-law residual is exactly `0`, so the source is algebraically self-consistent — the failure is geometric (non-convex), not a thermodynamic-identity or mass-convention problem. No interpolant that reproduces the source endpoints can have `dn_B/dμ_B ≥ 0` on these intervals. This is a `SOURCE_TABLE_CONVEXITY_FAIL`, not a generic `STABILITY_GATE`.

**2. What the same-mass v2 comparison proves and what it does not prove.**
It proves that the two representation-only transformations (PDG-gauge and pair-core gauge-shift) produce **identical physics verdicts on all 65 pairs** (65/65 agree on the structured category), so the mass-convention re-expression is representation-only and does not alter the pressure-density convexity condition. It does NOT prove that a Maxwell construction exists — 61/65 pairs fail with `INTERPOLANT_SHAPE_FAIL` (the §1 source-convexity diagnosis), 2 with `BETA_EQUILIBRIUM_FAIL`, 1 `CAUSALITY_FAIL`, 1 `NO_COEXISTENCE_CROSSING`. The old text-based taxonomy gave 4 pairs inconsistent labels (same failure, different category); the structured taxonomy fixes this.

**3. Whether native M7 CSS output is thermodynamically valid within its declared scope.**
Yes, within its declared thermo-only barotropic scope. The bridge satisfies first-law (1.1e-16), strict monotonicity, the exact analytic `∫n dμ = P2−P1` (4.6e-12), and causality (cs² ∈ [3.3e-5, 0.859]). The bridge is **C⁰ continuous** at the crust-core seam but **NOT C¹** (source crust cs² ≠ bridge cs²). `composition_dropped=true`; C3/C4 are `NOT_AVAILABLE`, not PASS — the output is not a composition-preserving physical EOS. These limitations are preserved in the report, not hidden.

**4. Which M5A CSS validation claims and artifacts require correction.**
All 7 audit findings are corrected: (1) the Python `check_consistency` is no longer recorded as `native_c1_c11_overall_pass` (renamed; native verification marked BLOCKED); (2) the misnamed `GPPVA_DD2_m5a_real_css.json` is replaced by `BBB_BHF2_m5a_real_css.json`; (3) the invented `mp=M_N*0.9938=933.74 MeV` is replaced by the real BBB_BHF2 `mp=938.272 MeV`; (4) the first-law reset is quantified (max 1.77e-8); (5) the μ_t mismatch is recorded per transition density; (6) Δε=0 is explicitly `CONTINUOUS_UNCLASSIFIED` (not inferred as crossover/2nd-order); (7) the (n,μ,P,ε) JSON is declared a barotrope, not a CompOSE-generated EOS.
