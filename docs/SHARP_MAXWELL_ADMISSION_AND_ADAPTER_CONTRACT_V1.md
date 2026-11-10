# Sharp-Maxwell 65-Pair Cohort: Source Admission and Adapter Contract v1

**Date:** 2026-11-10  
**Research-loop stages:** C2 (source/code/equation mapping), C3 (numerical adapter design)  
**Manifest:** `data/sharp_maxwell_65_cohort_manifest.v1.json`  
**Scope lock:** Do not run the 65 root searches, create EOS tables, or invoke TOV/Cowling until the source-admission and multi-segment envelope blockers below are resolved.

## 1. Decision and cohort identity

The new scientific cohort is `sharp_maxwell`, with 13 crust inputs crossed with 5 core inputs (65 pairs). Its manifest and outputs are independent of the historical six-method `maxwell` rows. Historical products remain read-only provenance; none may be patched, used as an initial guess, or promoted into a sharp result.

The 65-pair source grid is present in the raw CompOSE EOS catalog at `eos_generator/eos_data/` on the PR #3 head `a09ec30d74b363652d2d1aefe990557c468679dd`. The 65 legacy result directories are present on `main` at `matched_eos_all_combinations/matched_eos_390/` (upload commit `261fcdcbeaf6b5d842480e82e249312e8da5532a`). The directory cross product was checked: 5 cores, 13 crusts, 65 unique pair directories, no duplicate pair names.

## 2. Why historical maxwell is not a Maxwell construction

The legacy summary contains 65 rows labelled `method=maxwell` and `spline_mode=maxwell_plateau`. Its two reported interval endpoints do not obey coexistence. Using
`r_x = abs(x_2-x_1)/max(abs(x_1),abs(x_2))`,
the pressure mismatch spans 0.43796–0.86546 (median 0.69700), and the baryon-chemical-potential mismatch spans 0.004081–0.022816 (median 0.009329). Thus none of these endpoint pairs supplies a common `(P_t, mu_t)`. Nonzero reported `delta_epsilon` or `delta_n_B` does not repair that failure.

For a cold Maxwell transition with one globally conserved baryon number and local charge neutrality in each phase, the new output must satisfy
\[
P_L(\mu_t)=P_H(\mu_t)=P_t,\qquad
n_\pm(\mu_t)=\frac{dP_\pm}{d\mu_B},\qquad
\epsilon_\pm=\mu_t n_\pm-P_t,
\]
and therefore
\[
\Delta\epsilon=\epsilon_H-\epsilon_L
=\mu_t(n_H-n_L)=\mu_t\Delta n_B.
\]
The output should carry the two one-sided coexistence rows at exactly common pressure and chemical potential, with distinct density and energy density. A pressure plateau between unequal endpoint pressures or a smooth bridge is not this object.

## 3. Raw-source preflight findings

All values below were read from the raw `eos.thermo` and `eos.compo` files at the pinned PR #3 head. File blob SHAs and row counts are recorded in the manifest. The formula used for the source chemical-potential cross-check was
\[
\mu_B^{(Q_3)}=m_n(1+Q_3),\qquad
\mu_B^{(E)}=\frac{\epsilon+P}{n_B}
=m_n(1+Q_7)+Q_1.
\]
The reported residual is the maximum over knots of
\[
\frac{|\mu_B^{(Q_3)}-\mu_B^{(E)}|}
{\max(|\mu_B^{(Q_3)}|,|\mu_B^{(E)}|)}.
\]
This is explicitly a raw Q3-versus-Euler residual; do not conflate it with a derivative-based C2 residual from another audit.

### 3.1 Crust grid

Every one of the 13 crust files has an `eos.compo` table with matching row count. In all 13, the observed `Iphase` histogram is 217 rows of `0`, followed by the remaining rows of `7`; the change is between source knots 217 and 218. Each composition row has four particle pairs and one nucleus quadruple. These IDs are provider data, not universal physical phase labels. No derivative stencil or interpolation segment may cross this boundary without a documented physical interpretation.

| Crust family | Count | Header masses \((m_n,m_p)\) MeV | Nadd | Max Q3–Euler residual |
| --- | ---: | ---: | ---: | ---: |
| SPG_M1–M5 | 5 | (938.9, 938.9) | 0 | 1.411e-6 |
| GPPVA_DD2 | 1 | (939, 939) | 1 | 8.758e-5 |
| GPPVA_DDME2 | 1 | (939, 939) | 1 | 8.978e-5 |
| GPPVA_FSU2 | 1 | (939, 939) | 1 | 1.199e-6 |
| GPPVA_FSU2H | 1 | (939, 939) | 1 | 1.538e-6 |
| GPPVA_FSU2R | 1 | (939, 939) | 1 | 2.597e-6 |
| GPPVA_NL3wrL55 | 1 | (939, 939) | 1 | 3.214e-6 |
| GPPVA_TM1e | 1 | (938, 938) | 1 | 3.317e-6 |
| GPPVA_TW | 1 | (939, 939) | 1 | 1.476e-6 |

For the eight GPPVA tables, the one additional `eos.thermo` value takes values in 0–3 (model-specific histograms are in the manifest). The additional column is not an ignorable padding field. Its Data Sheet meaning has not been established in this pass; retain every value and block generic interpolation or deletion. Under the proposed strict adapter gate of `2e-6`, five of the eight GPPVA crusts fail the Q3–Euler check (DD2, DDME2, FSU2R, NL3wrL55, TM1e); FSU2, FSU2H, and TW pass this one check. All five SPG tables pass it. Passing this one scalar check is not source admission.

### 3.2 Core grid

| Core | Rows | \((m_n,m_p)\) MeV | eos.compo | Nadd | Max Q3–Euler residual | Max \(|Q_5|\) | Previously reported max \(c_s^2/c^2\) |
| --- | ---: | ---: | ---: | ---: | ---: | ---: | ---: |
| APR_APR | 172 | (939.56533, 938.272) | absent | 0 | 3.62e-16 | 0 | 1.4515 |
| BBB_BHF2 | 84 | (939.56533, 938.272) | absent | 0 | 3.62e-16 | 0 | 0.8440 |
| GM_GM1 | 349 | (939.56533, 938.272) | present, zero pair/quadruple populations in the source parser | 0 | 2.73e-11 | 0 | 0.7562 |
| RG_SLY4 | 197 | (939.5654, 938.272) | present | 1 | 1.128e-4 | 1.383223e-3 | 0.9493 |
| RG_KDE0v1 | 1240 | (939.56533, 938.272) | present | 1 | 1.151e-4 | 6.03452e-4 | 1.1023 |

The quoted maximum sound speeds are carried forward from the preceding source-admission review and were not independently recomputed in this C2/C3 pass. The raw Q3–Euler residuals above were recomputed from the pinned files. APR_APR and RG_KDE0v1 therefore retain a reported acausal high-density tail; any causal-domain restriction must be separately declared and reported alongside full-source results. RG_SLY4 and RG_KDE0v1 also fail the v1 cold beta-equilibrium gate because their source `Q5` is nonzero. The APR and BBB sources have no composition file; GM has a file but the parser detects no nonzero composition populations. These cases must not be conflated.

### 3.3 Rest-mass references

The crust headers use neutron masses 938.9, 939, or 938 MeV, while the cores use approximately 939.565 MeV (and proton mass 938.272 MeV). The raw `Q3` and `Q7` values are normalized by each table's own header mass. The adapter must reconstruct dimensional \(\mu_B\) and \(\epsilon\) using the source header first, then document any common-reference conversion as a thermodynamic transformation. It must not copy a gauge shift from the old smooth matcher or silently replace one source's mass convention with another.

## 4. Current native-kernel and solver contracts

The PR #2 Fortran constructor builds cubic-Hermite \(P(\mu_B)\) branches with \(P'=n_B\), isolates stationary points of the pressure difference, and accepts exactly one crossing. Its static v1 contract is two smooth, cold, beta-equilibrated barotropes with constant `Iphase`, no composition file, and `Nadd=0`. Its internal Q3–Euler tolerance is `5e-4`; that is a kernel validation ceiling, not the scientific cohort admission policy proposed here.

The v1 builder explicitly rejects an input `Iphase` change and rejects nonzero `Nadd`. The strict adapter should use a `2e-6` Q3–Euler gate by default, matching the existing general C2 acceptance level. Any projection of Q3 to Euler must be an explicit repair mode with raw residuals, per-knot deltas, and transition-location sensitivity recorded; no silent use of the kernel's looser tolerance.

The PR #3 production TOV/Cowling path currently accepts two normalized EOS table providers, high-density core and low-density envelope, and one sharp interface. It integrates the core to \(P_t\), restarts the envelope at the same geometry, and checks equality of pressure and chemical potential plus the first-order jump identity. Slow- and rapid-conversion interface transfers are distinct. It does not yet represent the two retained crust source segments plus the core as a multi-segment EOS chain. Therefore, even a successful two-branch root would not by itself authorize a TOV/Cowling run if the adapter had discarded or smoothed the crust's internal source boundary.

A separate implementation issue is confirmed in `oscillation_source/solutions/src/eos/ns_table_eos.f90`: `reconstructed_n_max_rel` is computed after First-Law density reconstruction, but it does not enter `fully_admitted`. The previous direct-provider receipts reported approximately \(7.3\text{–}7.7\times10^{-4}\) drift. This manifest requires that value to be retained and gated, but deliberately does not invent a numerical tolerance before refinement and observable-sensitivity tests establish one.

## 5. Adapter contract

### A. Preserve and segment source data

1. Keep the original CompOSE directory and hashes immutable. Emit a separate normalized adapter artifact with a provenance receipt.
2. Parse `eos.nb`, `eos.thermo`, `eos.t`, `eos.yq`, and `eos.compo` together. Check grid-key coverage and the Nadd column count at every knot.
3. Preserve the provider-defined `Iphase` labels and `Nadd` values. Obtain and record the model Data Sheet interpretation before assigning semantic phase names.
4. Split interpolation domains at every `Iphase` change. Never smooth, differentiate, or construct a Hermite interval across a segment boundary. The 217/218 crust boundary is an observed source boundary, not automatically a Maxwell transition.
5. Preserve particle-pair and nuclear-quadruple content. For the CompOSE nuclear quadruple, retain `Y_I`; convert to mass fraction only with `Y_h=Abar*Y_I` and keep the raw value as provenance.

   **Provider-conformance caveat from C6:** the conversion above is the standard CompOSE interpretation, not an already validated mapping for every source table. For the pinned SPG_M1 and GPPVA_FSU2 rows, blindly applying `Abar*Y_I` produces values above unity, while the provisional direct mapping `Y_h=Y_I` fails native C3/C4 on every `Iphase=7` (inner-crust) row. Their composition semantics remain `BLOCKED_REVIEW_REQUIRED`; see [the C6B composition semantics review](SHARP_MAXWELL_COMPOSITION_SEMANTICS_REVIEW_V1.md). Do not apply either mapping to these providers as an admitted adapter conversion until the provider-level review is complete.

### B. Thermodynamic normalization and admission

1. Compute `P=n_B Q1`, `epsilon=n_B m_n(1+Q7)`, `mu_Q3=m_n(1+Q3)`, and `mu_E=m_n(1+Q7)+Q1` from the raw headers.
2. Require cold `T=0`, entropy `Q2≈0`, beta-equilibrium `Q5≈0`, and `Q6≈Q7` for the strict beta-equilibrated path. A nonzero `Q5` is an admission failure for that path, not something to project away.
3. Apply the adapter-level `max |mu_Q3-mu_E|/scale <= 2e-6` gate by default. A Q3 projection is opt-in and must emit an auditable repair receipt and a before/after sensitivity study.
4. Record `reconstructed_n_max_rel` and make its acceptance criterion a separate gate. Until the criterion is selected and tested, classify affected cases as `DENSITY_DRIFT_REVIEW`, not `FULLY_ADMITTED`.
5. Keep full-source and causal-domain audits separate. No implicit high-density truncation and no causal completion under the original model identity.
6. Require explicit provenance for local charge neutrality in each phase. A composition-free one-dimensional table does not independently prove it.

### C. Sharp root and output

On each pair of smooth, admitted segments, construct monotone \(P(\mu_B)\) with \(P'=n_B\). Isolate every stationary point of \(\Delta P(\mu_B)\), partition into monotone intervals, and count every root. Only `UNIQUE_ROOT` proceeds. Preserve `NO_ROOT`, `MULTIPLE_ROOT`, and input/causality failures as scientific results. The output stores common \((P_t,\mu_t)\), one-sided \(n_\pm,\epsilon_\pm\), \(\Delta n_B\), \(\Delta\epsilon\), Euler residuals, and the residual of \(\Delta\epsilon-\mu_t\Delta n_B\). Output files live only under `sharp_maxwell/`.

### D. TOV/Cowling gate

Before TOV, represent the crust's retained source segments and the sharp core-crust interface without deleting, blending, or crossing source boundaries. This likely needs a multi-segment low-density envelope/provider or a generalized phase-chain TOV background. Then verify the stable stellar branch and Seidov onset diagnostic. Only afterward compute slow- and rapid-conversion Cowling spectra separately.

## 6. Pair-level status and next gates

All 65 pair rows are present in the manifest. Their current status is `BLOCKED_BEFORE_ROOT_SEARCH`, not `NO_ROOT` and not `INPUT_REJECTED`: the root search has not been called. All 65 inherit the crust multi-segment/composition blocker; additional blockers are pair-specific (Nadd semantics, Q3 residual, core Q5, composition mode, acausal tail, and rest-mass reference).

The next allowed stage is C4 source validation: verify all row/key alignments, source-level monotonicity, first-law and composition constraints, segment-boundary behavior, and whole-interval causality; emit row/interval-specific receipts. Then resolve the Data Sheet semantics and decide the multi-segment envelope interface. No root search, TOV, or mode spectrum is part of this checkpoint.

## 7. Verification record

- **Numerically inspected:** 13 crust and 5 core raw CompOSE thermodynamic/composition metadata; file blob SHAs and row counts pinned in the JSON manifest.
- **Cross-product checked:** 65 unique pairs (5 × 13); all pair rows generated in the manifest.
- **Legacy audit:** all 65 historical maxwell rows have endpoint pressure and chemical-potential mismatches as summarized above.
- **Equation check:** \(\Delta\epsilon=\mu_t\Delta n_B\) follows algebraically from \(\epsilon_\pm=\mu_t n_\pm-P_t\); symbolic identity was evaluated.
- **Not done:** no 65-pair sharp root sweep, no new sharp EOS tables, no TOV/Cowling sequence, no mode spectrum, and no claim of full source admission.
- **Execution harness:** no native subagent/delegate tool was available in this session; the C2/C3 work was executed directly and sequentially with source and provenance references.

## References

- [CompOSE Quick Guide for Providers](https://compose.obspm.fr/download/pdf/CompOSE_Quick_Guide_for_Providers.pdf)
- [CompOSE Reference Manual v3.00](https://compose.obspm.fr/download/pdf/manual_v3.00.pdf)
- [PR #2: native sharp-Maxwell constructor](https://github.com/cosmosapjw-quantum/NS_Cowling/pull/2)
- [PR #3: phase-segmented TOV/Cowling](https://github.com/cosmosapjw-quantum/NS_Cowling/pull/3)
