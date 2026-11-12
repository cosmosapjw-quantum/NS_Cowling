# NS_Cowling technical audit and continuation report

**Audit date:** 2026-11-12  
**Disposition:** PARTIAL / NOT RELEASE-READY  
**Repository:** cosmosapjw-quantum/NS_Cowling  
**Cloud execution branch:** cloud/eos390-fortran-unification-v1  
**Current audit baseline:** 0cbe52ac1ace0535dea541beae31a49410abdbe5  
**Sharp-Maxwell reference branch:** research/maxwell-native-fortran-v1  
**Sharp-Maxwell PR:** #3, head 99a0f9c314cab996c3a0bac28299e7b604b6cd01  
**EOS390 PR:** #13, head 0cbe52ac1ace0535dea541beae31a49410abdbe5  

## 1. Executive finding

The repository contains substantial native Fortran EOS functionality and a separately developed sharp-interface Cowling reference path, but the claimed end-to-end regeneration of the 390 crust-core outputs has not been demonstrated. The EOS390 M7 run records 390 task receipts (150 PASS, 30 FAIL, 210 BLOCKED), but source inspection establishes that all 150 PASS records are not production-valid: 30 CSS tasks write the unmodified core table plus a CSS metadata sidecar, and 120 cubic/quintic tasks copy the crust table and mark the result as PASS despite an explicit stub marker. The M7 summary therefore has zero demonstrated complete native Fortran combined EOS outputs. M9's later 26/364 status report is inconsistent with the versioned M7 summary and cannot be reconciled without the original run inputs and logs. Both records are retained as historical evidence; neither is accepted as a release claim.

The current implementation fixes the row-level composition-presence contract and adds fail-closed M7 behavior. At EOS branch head 0cbe52a, the Fortran build passed and the Python/pipeline suite reported 104 passed, 2 skipped, 2 failed. The remaining failures were a blocker-message regex mismatch and a NumPy boolean identity assertion in a regression test; neither is evidence of a physical solver failure, but the whole CI is still FAIL. The follow-up test-only correction is in this audit continuation and requires a fresh exact-head run. Do not describe the EOS390 corpus as regenerated or Dropbox-backed until a fresh source-locked native run completes all read-back and provenance gates.

The separate sharp-Maxwell Cowling branch has bounded analytic/synthetic Fortran validation. It is not yet integrated with the EOS390 batch backend and does not establish realistic source-backed hybrid-star spectra, physical mode identities, or production root convergence.

## 2. Repository state and evidence pins

| Item | Exact state | Interpretation |
| --- | --- | --- |
| PR #13 | open, unmerged; cloud branch head 0cbe52a | EOS390 unification and presence-mask fixes; CI result must be checked at exact head |
| PR #3 | open, unmerged; head 99a0f9c | separate native sharp-Maxwell TOV/Cowling reference path |
| EOS390 top-level DAG | data/eos390_native_fortran_unification_dag.v1.json | physical/source contracts remain a subsystem reference; this report's unified cloud DAG supersedes execution sequencing |
| M7 current summary | download/m7_batch/M7_BATCH_SUMMARY.json | 390 receipts: 150 PASS / 30 FAIL / 210 BLOCKED |
| M9 old report | download/m9/M9_HANDOFF_REPORT.json | historical 5/5 PASS claim is superseded by this adjudication |
| M9 pair manifest | download/m9/M9_PAIR_METHOD_STATUS_MANIFEST.json | 390 rows but 26 PASS / 364 FAIL, inconsistent with M7 summary |
| v3.2 historical outputs | reported 2,109 files | immutable reference; old M9 did not independently prove unchanged state because its check set the result true by construction |

The M7 PASS distribution is 30 each for css, cubic_F, cubic_P, quintic_F and quintic_P. There are no M7 PASS rows for its maxwell method in the recorded summary. The 120 smooth PASS rows are explicit crust-copy stubs. The 30 CSS PASS rows contain eos.css_branch.json and the original core eos.nb/eos.thermo/eos.yq/eos.t; the CSS branch samples are not assembled into a combined crust-core table. Consequently, M7's orchestration/cardinality receipts are evidence of task dispatch and failure retention, not evidence of 150 physical EOS products.

## 3. Method ontology: discontinuity is not interpolation

The historical six-method EOS390 matrix is retained as an immutable identity: maxwell, css, cubic_F, cubic_P, quintic_F, quintic_P. However, the machine DAG records that the existing public token maxwell historically denotes a project-defined smooth near-pressure-plateau Hermite bridge (NPPHB), with continuous density and no latent heat. It is not a thermodynamic sharp Maxwell transition and is not a published construction attribution. The current M7 generator incorrectly labelled this token sharp_first_order and dispatched the Python M5 sharp-Maxwell oracle. The new M7 contract restores the legacy meaning and blocks the method until its native smooth kernel is ported.

A genuine locally neutral, one-conserved-charge sharp Maxwell transition is a separate physical operation. It requires common coexistence pressure and baryon chemical potential, a positive density jump, and the first-law jump identity:
[
P_-(mu_t)=P_+(mu_t)=P_t,quad mu_{B,-}=mu_{B,+}=mu_t,
]
[
n_+>n_->0,qquad Deltaepsilon=mu_tDelta n_B>0.
]
Keep the two one-sided coexistence states adjacent. No spline, finite-difference stencil, interpolation interval, smooth kernel, or root search may cross this interface. If this operation is added to a new task matrix, use an explicit new token such as sharp_maxwell and a new versioned output manifest; do not silently repurpose the historical six-method token or change the historical 390-key identity.

CSS is a parameterized high-density extension. A positive Delta_epsilon creates a sharp first-order jump only when the chosen parameters and first-law conditions are supplied and verified. CSS method names, interpolation method labels, provider Iphase integers, and a smooth curve do not independently classify physical transition order. A zero jump is CONTINUOUS_UNCLASSIFIED unless independent thermodynamic evidence establishes a second-order transition or a crossover.

## 4. EOS backend audit

### 4.1 Fortran/Python boundary

The target architecture remains Fortran numerical backend with a Python frontend restricted to path/argument validation, serialization, native executable invocation, report decoding, hashing and presentation. No production Python numerical fallback is allowed. The native CLI currently exposes read, verify, normalize, write, fill-gap, reconstruct-master, causal-complete, match-core and maxwell-transition operations. The native maxwell-transition constructor is a thermo-only operation and rejects composition-bearing input rather than dropping fields. The CSS segment kernel exists natively, but a full native CSS crust-core join and CLI/read-back path are not established.

### 4.2 Row-level composition presence

The previous global composition_available flag conflated absent records, partial records and explicitly written zero fractions. The new implementation stores per-row masks for particle IDs and the nuclear-group tuple, plus independent C3/C4 applicability. It also prevents the Python frontend from fabricating all-true masks when a parser-originated object has lost its masks. Programmatically constructed EOSData records get explicit initialization under a documented constructor contract; parsed objects keep the exact source masks.

The CompOSE v3.01 manual defines eos.compo rows as grid key, provider-local Iphase, Npairs particle-index/fraction pairs, and Nquad nuclear-group quadruples. Npairs=0 and Nquad=0 are permitted. Only relevant particle fractions need be listed; omitted species are not a license for an internal parser to forget whether a value was actually present. Group tuples contain an identifier, average mass number, average charge number and group number fraction. The manual distinguishes this data structure from provider-defined Iphase labels. See Typel et al., CompOSE Reference Manual v3.01 (2022), arXiv:2203.03209, sections 3.3–3.4 and 4.2.

A sparse-writer defect was found: Nquad=0 was emitted twice on one output path, producing trailing fields on native read-back. The writer was corrected to emit the count exactly once. Duplicate particle/group identifiers in a row are rejected. Explicit zero particle fractions are retained as present values. The current CI test sequence also caught a dissolution edge: C4 is inapplicable once the representative nuclear group dissolves (Abar=0, so Z/A is undefined), while C3's generated baryon inventory remains meaningful. That gate is recorded separately; do not cure it by tolerance relaxation or relabelling absent data.

### 4.3 Provider semantic blocker

SPG/GPPVA final group tuple semantics remain BLOCKED. The available sampled rows are structurally well-formed but the final tuple value does not admit a universal, provider-supported mapping to standard YI across outer and inner crust. Do not guess YI, Xh, Ye, neutron-gas fraction or any other aggregate. Do not relax C3/C4, synthesize Yn/Yp, or drop the composition file. Obtain written definitions for Iphase=0 and Iphase=7, Abar/Zbar, the group fraction and the neutron-gas representation. This is a scientific/data semantic blocker, not a file-width repair.

### 4.4 M2/M2A/M10 and native batch

M2 smooth Hermite interpolation, M2A row-level composition coverage and M10 general multidimensional CompOSE axes/field preservation are explicit dependencies. M2A code changes are in PR #13 and require exact-head CI completion. M2/M10 native kernel ports are not complete. Nadd and q_extra are provider-specific extensions; observed GPPVA_FSU2H q_extra values include categorical pasta-geometry codes and must not be averaged or interpolated as continuous observables. Preserve raw values and provider declarations.

M7 now fails closed for unavailable native methods. Historical maxwell/NPPHB, CSS, cubic_F/P and quintic_F/P remain BLOCKED until their actual native backend and full read-back gates exist. A native sharp Maxwell Fortran kernel is a separate operation and must not be counted as a legacy NPPHB matrix result.

## 5. Sharp-Maxwell Cowling audit

PR #3 is a separate open branch with a native analytic sharp-Maxwell interface implementation and bounded synthetic tests. Its branch has a successful native CI run for the phase-side fixture/import correction. The exact bounded audit in docs/ADVERSARIAL_AUDIT_MAXWELL_TOV_COWLING_20261109.md records synthetic TOV/Cowling tests, slow/rapid transfer comparisons and a generated two-table fixture. It explicitly does not establish realistic source-backed CompOSE hybrid EOS execution, stable-branch survey, positive physical root family continuation, or f/p1/interface-g mode identity.

For the 2x2 Cowling system y' = M y, y=(W,V)^T, same frequency/background/ell/closure, and a smooth cold barotropic phase:
[
mathrm{tr}M=
u'(2+1/chi),quad h' / h=-
u'(1+1/chi),
]
so Liouville's determinant identity gives
[
Omega=h e^{-
u}(W_1V_2-V_1W_2)=mathrm{constant}.
]
Across a linear interface y_o=T_{i\to o}y_i:
[
rac{Omega_o}{Omega_i}
=rac{h_o}{h_i}e^{-(
u_o-
u_i)}det T.
]
This is a reduced-system bilinear Abel-Liouville concomitant, not a claimed canonical covariant symplectic current or physical energy flux. The lapse factor must be retained if nu is discontinuous. The state transfer orientation must be recorded.

The physical slow/material and rapid/equilibrated endpoints are distinct. The LA10 theta-family is an algebraic interpolation of kinematic boundary constraints only; theta is not a microscopic conversion rate, timescale or calibrated kinetic coefficient. There is no reaction flux, affinity, entropy-production law or extra interface degree of freedom in that construction. Its no-thin-shell admissibility remains unresolved. Never call a theta scan a finite-rate kinetics scan.

LA11's finite-interval matching determinant is not automatically a whole-line Evans function. Analyticity of the actual native ODE, endpoint frames, interface maps, gauge/scale factors and a regular Fredholm pencil must be proved before using argument-principle root counts or identifying determinant zero order with physical algebraic multiplicity. QR norms/Gram-Schmidt may introduce complex conjugation and non-holomorphic normalization; determinant phase and every basis/scale cocycle must be retained.

The analytic fixture root values stored in earlier notes are sorted root sets, not physical family labels. Mode identity remains UNRESOLVED until phase-local node counts, eigenfunction/profile overlaps, tangential polarization, interface one-sided values and parameter continuation are checked. Do not identify f, p1 or discontinuity-g modes from frequency ordering alone.

## 6. CI and validation status

The current PR #13 branch CI includes a full Python frontend regression suite and native Fortran build. Recent exact-head attempts have shown:
- build success;
- successive parser-mask failures reduced from 11 to 2 and then the M7 fail-closed TDD tests were added;
- PR13 run 43229325660 failed with 2 tests after the first sparse-writer fix; the failure was one contradictory explicit-zero assertion and one C4 dissolution applicability case;
- commit 1ab9fcf addressed both; its CI surfaced only the Python bool identity assertion in the dissolution test plus the newly added M7 red tests on the subsequent test-only commit;
- commit 0cbe52a contains the M7 implementation and taxonomy corrections; its exact-head CI must be read before this report is promoted to PASS.

The checksum file SHA256SUMS_v4.txt currently contains stale expected hashes for several intentionally edited files. Diagnose and update those entries from the exact committed file bytes; do not delete checksum coverage or call the manifest verified while mismatches remain.

## 7. Release gates and final disposition

No new 390-output regeneration is claimed. M9's earlier 5/5 PASS summary is superseded because C2 read-back did not establish native generation, C4 used a by-construction true flag rather than an independent immutable-baseline comparison, C5 did not establish a Dropbox upload acknowledgement, and the reported M9 pair-status file disagrees with the M7 receipt distribution. Preserve every historical artifact; attach a corrected adjudication receipt instead of rewriting it.

Release requires: (1) source-locked input inventory and paper checksums, (2) provider composition definitions and closure, (3) complete Fortran M2/M2A/M10 kernel and lossless multi-D I/O, (4) a manifest whose historical method identities are exact, (5) native Fortran generation of complete joined tables, (6) independent CompOSE read-back of every PASS output, (7) differential/physical/numerical M8 gates, (8) verified historical-output immutability against a trusted baseline hash map, and (9) actual Dropbox upload acknowledgement for the new versioned output set. The cloud agent may implement/test code and CI now; Dropbox source intake and full 390 numerical generation require the authorized local/cloud filesystem that contains the user's Dropbox mount.

**References**
- Typel et al., CompOSE Reference Manual v3.01 (2022), [arXiv:2203.03209](https://arxiv.org/abs/2203.03209).
- Tonetto & Lugones, Phys. Rev. D 101, 123029 (2020), [arXiv:2003.01259](https://arxiv.org/abs/2003.01259).
- Alford, Han & Prakash, Phys. Rev. D 88, 083013 (2013), [arXiv:1302.4732](https://arxiv.org/abs/1302.4732).
- Dong, Kuan & Yunes, 2026, [arXiv:2609.33002](https://arxiv.org/abs/2609.33002). This is full-GR context, not a drop-in Cowling derivation.
- The NPPHB legacy parameterization is project-defined and is not attributed to those papers.

## 8. Continuation 2026-11-13 (branch cloud/eos05-native-smooth-parity-v1)

Base: `08fe5b814cbe2b3d6e41b9aa44d2e5f687bc1351`. This section is additive; earlier sections are historical evidence.

- **CI at the pinned head.** Runs 43232391185 and 43232373861 at `08fe5b8` succeeded (build, pytest, native extrema/CSS tests, SHA manifest). The pytest count (106 passed, 2 skipped) is from the GitHub job log and was reproduced locally (106 passed, 2 skipped, exit 0; the two skips are the M10 real-input tests because `inputs/` is not staged). The `latest_ci` field of the 2026-11-12 audit receipt is superseded by `data/receipts/ns_cowling/NS_COWLING_AUDIT_ERRATA_20261113.json`.
- **Errata.** The protected M9 adjudication lists `m7_pass_by_method` without `quintic_P` (sum 120 vs 150 PASS). The corrected map, derived from the raw M7 receipts, is css/cubic_F/cubic_P/quintic_F/quintic_P = 30 each. The M7 `maxwell` rows are 30 FAIL and 35 BLOCKED. These remain historical labels: none of the 150 is a valid complete native combined EOS.
- **EOS-05 is not a port task.** The native cubic/quintic F/P bridges and the NPPHB bridge already exist in `compose_eos_matching.f90`.
  - EOS-05a: 104 native-vs-frozen cases (13 crusts x {BBB_BHF2, GM_GM1} x 4 smooth tokens) are `PARITY_MATCH` (worst max relative difference 1.2e-12; branches exact_single_cubic 21, monotone_1knot_cubic 31, exact_single_quintic 9, adaptive_2knot_c2_hermite 43).
  - EOS-05b: the historical `maxwell` key (LEGACY_NPPHB, smooth, no latent heat) is routed to the native kernel under the internal name `legacy_npphb`; 26 pairs `PARITY_MATCH`. Bare `maxwell` stays rejected. `sharp_maxwell` is untouched.
  - **Not a production PASS:** composition is dropped in every case; there is no C1-C11 record, no combined-composition output and no read-back claim. M7 stays fail-closed.
- **Open blocker F-EOS05A-1.** All 52 RG_SLY4 smooth-token cases are rejected by the native guard `beta_equilibrated = max|Q5| < 1e-8` because the RG_SLY4 source has a constant `Q5 = 0.001383223`. RG_KDE0v1 also reads `beta_equilibrated=False`; its matching was not run. 156 of the 390 historical tasks use an RG core. The guard was not relaxed: the meaning of the constant Q5 offset is a source-semantics question.
- **Legacy NPPHB property.** Bridge free energy uses 16-point partial integrals while the first-law root uses a 64-point rule. Against a dense oracle the bridge F error is up to 1.4e-3 of the endpoint dF and the endpoint residual 2.4e-5 of dF; the native first-law residual (same 64-point rule) reads 1.3e-15 and cannot see it. Characterized, not fixed (historical identity).
- **Policy supersession.** `legacy_npphb` follows this cloud DAG, not `data/eos390_native_fortran_unification_dag.v1.json` lines 143/968/1227 (that file is pinned and was not edited). The rejection text of bare `maxwell` still promises a future sharp reuse of the token; this wording drift is logged (F-EOS05B-2).
- **`css` is unchanged** (user decision): historical `css` remains as defined by the cloud DAG and EOS-06 stays BLOCKED. Source evidence that legacy `css` was a continuous two-segment bridge is recorded here as an unresolved ontology question, not a routing change.
- **No Dropbox operation was performed; the 390 matrix was not regenerated.**
