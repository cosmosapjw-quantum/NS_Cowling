# Mass-dependent Cowling research DAG (2026-11-09)

## Authority

- Repository: cosmosapjw-quantum/NS_Cowling
- Research branch: research/mass-dependent-cowling-v1
- Base SHA at branch creation: a09ec30d74b363652d2d1aefe990557c468679dd
- PR #3 remains open/draft/unmerged and is not modified.
- Frozen smooth cohort: 13 crusts × 3 admitted cores × 6 smooth methods = 234 rows / 39 paired core-crust combinations.
- Old spectra have only the fixed target mass 1.3300106263 solar masses; reuse audit says zero modes were recomputed.
- Historical method token maxwell is a smooth bridge, not a sharp density jump.

## Dependency graph

\[
\text{frozen cohort and provenance}
\to TOV\ M(P_c)
\to \text{connected rising branch}
\to \text{adjacent-}P_c\text{ target bracket}
\to \text{native mass-root solve}
\to \text{Cowling roots and W/V identity}
\to \text{six-method paired maps}
\to \text{refinement and benchmarks}
\to \text{publication decision}.
\]

Mechanism and inference extensions branch from the verified mass-dependent spectrum. They are not silently included in the present implementation.

## Node ledger

| ID | Node | Status | Exit evidence |
|---|---|---|---|
| MD-00 | Freeze 234-row cohort, projected EOS tables, release summary and checksums | PASS_READ_ONLY | 13 crusts, 3 core families, 6 smooth methods. All old spectra at 1.3300106 solar masses. |
| MD-01 | Stable-branch extractor and paired estimands | IMPLEMENTED_CI_PENDING | Tests cover first turning point, failed-sample barriers, adjacent target brackets and no-imputation six-method pairing. |
| MD-02 | Exact-backend TOV-only sequence and explicit numerical limits | IMPLEMENTED_CI_PENDING | `solve_cowling_tov_star` calls `cell_solve task=star` with v5 `r0=10 m`, `rtol=1e-9`, `step_cap=0.05`; configurable radius limit defaults to historical 1e6 m and is set to 1e8 m for the exploratory scan. Explicit Pc brackets remain optional and preserve the old default. |
| MD-03 | Six-model high-tail pilot: BBB_BHF2 × GPPVA_FSU2, all six methods and eight masses | RETRY_PENDING | Run 42598463512 failed due a non-comparable TOV backend and the 1e6-m radius cap. Run 42599325312 passed helper tests but stopped at a `solve_model` radius-limit keyword/signature mismatch before native compilation/TOV. The signature correction is pending; baseline parity and branch coverage remain untested with the shared backend. |
| MD-04 | Full 234-EOS × eight-mass execution | NOT_RUN | Manual workflow dispatch with cohort=full; 48 log-pressure TOV samples per EOS. |
| MD-05 | Mass-response and coverage analysis across 39 fixed pairs | NOT_RUN | Complete-pair denominators, paired ranges, per-core/per-crust maps, method sign reversals and peak mass locations. |
| MD-06 | ODE/step/root refinement versus effect size | NOT_RUN | Frequency movement no larger than 0.1 of the smallest reported paired effect. |
| MD-07 | Direct recognized matching-prescription baselines | NOT_STARTED | Simple concatenation as a non-physical control; linear P(epsilon); Canullán-Pascual thermodynamic P(epsilon) interpolation; CUTER-like/documented thermodynamic crust reconstruction. |
| MD-08 | Independent Cowling/full-GR benchmark | NOT_STARTED | Representative low-, reference- and high-mass points; Cowling systematic offsets separate from matching shifts. |
| MD-09 | Eigenfunction/acoustic-cavity mechanism | NOT_STARTED | Matching-region overlap and acoustic travel time; variational kernel only after boundary/matching derivation. |
| MD-10 | Maximum mass, crust thickness/mass and tidal Lambda | PARTIAL_BLOCKED | TOV M,R sequences are in scope. Resolved turning point is only an Mmax candidate; crust event data and Love-number solver are not in this increment. |
| MD-11 | Source-backed CompOSE EOS and NICER/tidal inference | BLOCKED_UPSTREAM | Requires source provenance, thermodynamic admission and structural/tidal observables. |
| MD-12 | Publication claim gate | NOT_READY | JKPS/EPJA after MD-03–08, direct baselines and cohort coverage; A&A/MNRAS after MD-10–11; PRD only with new theory/full-GR result. |

## Analysis contract

For y in {f,p1}:
\[
\Delta_m^{(y)}(M;c,s)=100\left[y_m(M;c,s)/\bar y_6(M;c,s)-1\right].
\]
Only complete six-method cells enter paired summaries. No imputation. Report complete-pair count at every mass. The 234 rows are a finite crossed design, not independent microscopic EOS samples or a posterior.

## Stop and reject conditions

- Reject a target mass not bracketed within the connected rising-M(Pc) branch.
- Do not label a truncated branch endpoint as Mmax.
- A failed TOV sample is a branch barrier; do not bridge it.
- Reject a spectrum row with failed mass closure, failed shooting residual or unresolved/wrong W/V identity.
- Reject incomplete six-method cells from paired estimands; do not impute.
- Do not infer published-method novelty from the six smooth bridges alone; direct matching-prescription baselines are a separate gate.
- A converged Cowling root is not a full-GR QNM and supplies no damping time.
- The legacy maxwell method label is a smooth bridge, not a sharp density-jump construction.
- Do not claim crust mass/thickness, tidal deformability, observational inference or a sensitivity kernel from the current increment.

## Immediate next action

Re-run the six-model pilot after the same-cell TOV-only backend, configurable radius limit and validator correction are pushed. Inspect baseline M/R parity and all sequence samples first. Do not dispatch the full cohort unless the pilot has valid 1.0- and 1.33-solar-mass paired cells and the M/R, mass-residual, shooting-residual and node gates pass. Runs 42598463512 and 42599325312 have zero valid spectra and must not be cited as scientific results.


## Source-identity qualification

MD-02b reproduces the 35 paths and their order from the v5 compile manifest and hashes the current branch source bytes into a runtime build key. The release summary records source_zip_sha256 b2c0dcaf6a5b07956cf8d33c39651d4ade32f44b4c49b33aa08fd84fb458185b, but the Dropbox source archive was not independently downloaded/hashed in this session. Therefore source byte identity versus the archived release is UNVERIFIED; do not call it byte-identical. The next runner manifest records both the reference archive hash and current runtime build key.


## Additional theory/mock node (2026-11-09)

| ID | Node | Status at code-push preparation | Exit evidence |
|---|---|---|---|
| MD-13 | Smooth Cowling forward first variation plus fixed-mass composition | DERIVED_AND_ANALYSIS_MOCK_CHECKED; EXACT-HEAD_CI_PENDING | Matrix tangent S_zeta'=A S_zeta+A_zeta y; simple-root derivative -D_zeta/D_omega; deterministic synthetic ODE root/finite-difference comparison. No physical EOS kernel claim. |
| MD-14 | Whole-domain thermodynamic q-path mock admission | ANALYTIC_EXTREMA_CHECK_ADDED; EXACT-HEAD_CI_PENDING | Polynomial extrema over t in [0,1], q_min=4/5, (dq/dt)_min=137/550; exact pressure-integral endpoint constraint. |

The selected research-loop scope for this increment is R2/R4/R5/R7 plus C5/C6: literature triage, hypothesis/model refinement, analytic derivation, deterministic mock experiment, implementation and exact-head checks. Native subagent dispatch was unavailable and was not run. The next physical acceptance gate remains MD-03, the six-model pilot; MD-04 full-cohort execution is not authorized by a failed or unverified pilot.


## Branch-topology and radial-stability audit (2026-11-09)

| ID | Node | Status | Exit evidence |
|---|---|---|---|
| MD-03 | Six-model 24-point pilot | FAILED_PHYSICAL_GATE | Exact v5 baseline M/R passes at machine precision; first rising-proxy peak is only 0.7432–0.7434 Msun; zero admitted spectra; validator correctly rejects the required 1.0 Msun paired cell. Actions 42610037961, artifact 13668724637. |
| MD-15 | All monotonic M(Pc) candidate segments | IMPLEMENTED_CI_PENDING | New extractor partitions contiguous TOV blocks into rising/falling/flat segments, shares local turning points, preserves failed/missing-index barriers and marks every segment radial_stability_status=unresolved_requires_radial_mode_spectrum. |
| MD-16 | Direct radial eigenvalue branch classification | NOT_IMPLEMENTED_BLOCKS_SECONDARY_BRANCH_USE | Requires a radial Sturm-Liouville solver, explicit Gamma_1 closure, centre/surface boundary conditions, and interface transfer law where a true discontinuity exists. |

The 24-point artifact reveals a common low-compactness first maximum (about 0.743 Msun, R about 1232 km) followed by a large M/R excursion and a later rising interval reaching about 1.915 Msun. This is an observed coarse-grid topology only. It is neither a certified first-branch maximum mass nor proof of stability recovery. A new 96-point scan is staged to test resolution. No full 234-EOS run is authorized.

The radial stability gate is the nodeless fundamental eigenvalue omega_0^2, not the sign of dM/dPc alone. The selected perturbative closure must be declared: Gamma_eq=(epsilon+P)/P*dP/d epsilon for instantaneous barotropic equilibrium, or Gamma_1=(epsilon+P)/P*(partial P/partial epsilon)_(s,Y_i) for frozen composition. These need not coincide. The current finite-positive-pressure surface cutoff also needs a physical continuation or a derived finite-cutoff boundary condition before radial eigenvalues can be trusted.

The detailed equations, pilot table, acceptance gates and references are in docs/MASS_SEQUENCE_BRANCH_STABILITY_AUDIT_20261109.md.


## 96-point pilot execution and 192-point refinement (2026-11-09)

| ID | Node | Status | Exit evidence |
|---|---|---|---|
| MD-15 | Candidate-segment extractor | PASS_IMPLEMENTATION_VERIFIED | Actions 42611801783: mass-sequence tests 10 passed; all six 96-point TOV scans completed with 96/96 valid samples; segment records preserve all turns/gaps and stability=unresolved. |
| MD-17 | Six-model 96-point topology pilot | PASS_EXECUTION_HOLD_SCIENCE | Baseline M/R reproduction passes for all six; 0/48 f/p1 rows and no 1.0 Msun paired cell. First peak shifted materially from the 24-point scan. Artifact 13670407745, SHA-256 ecafc110e8b9d0196dfb4bb341da2ee54273f5d61767cca7c36e502a40737b7d. |
| MD-18 | Six-model 192-point topology refinement | QUEUED_BY_NEXT_PUSH | Compare first-turning-point M/R/Pc and all method-dependent turns against the 96-point values. Do not run spectra or label stability. |

At 96 points, the first sampled maximum ranges from 0.7650195 to 0.7651069 Msun at Pc = 0.226758–0.227326 MeV fm^-3, R = 1879.8–1889.0 km. The maxwell compatibility method, which is a smooth zero-density-jump bridge, exhibits an additional candidate max/min pair inside its low-cs2 bridge. This is an observed sequence feature, not a stable branch label. Full details and the 24-to-96 nonconvergence observation are in docs/MASS_SEQUENCE_BRANCH_STABILITY_AUDIT_20261109.md.


## 192-point pilot and independent TOV benchmark (2026-11-09)

| ID | Node | Status | Exit evidence |
|---|---|---|---|
| MD-18 | Six-model 192-point topology pilot | PASS_EXECUTION_HOLD_SCIENCE | 1152/1152 TOV integrations valid; all six frozen-v5 M/R anchors pass; 0/48 f/p1 rows; no 1.0 Msun paired cell. Artifact 13673520039, SHA-256 e460c295af38a70444e6139fb82aa037b3be4718e29d172db7be31d358abbff5. |
| MD-19 | Same-EOS independent SciPy TOV reference | IMPLEMENTED_CI_PENDING | New DOP853 direct-r SI TOV solver rebuilds native-knot cubic-Hermite P(mu), then compares M/R with v5 cell_solve at five selected Pc values. Acceptance: max relative M/R difference <= 1e-6. |
| MD-20 | 384-point six-model topology refinement | QUEUED_BY_NEXT_PUSH | Resolve remaining first-minimum radius and maxwell-bridge fold location movement. No radial stability labels or spectra. |

The 96-to-192 first-max mass shift is small (about 2.2e-4 relative), but first-max radius moves by about 1.6% and the first-minimum radius by more than a factor of two in the css row. This is not a converged mass-radius topology. The maxwell compatibility bridge's extra fold persists, with its intermediate maximum moving from Pc about 1.1164 to 1.0415 MeV fm^-3 and its intermediate minimum from Pc about 1.8976 to 1.9619 MeV fm^-3. This remains a smooth bridge, not a density-jump Maxwell construction. See docs/MASS_SEQUENCE_BRANCH_STABILITY_AUDIT_20261109.md.


## 384-point pilot and independent TOV diagnostic repair (2026-11-09)

| ID | Node | Status | Evidence |
|---|---|---|---|
| MD-20 | 384-point native six-model TOV topology | PASS_EXECUTION_HOLD_SCIENCE | 2304/2304 TOV integrations valid; all six frozen v5 M/R anchors pass; 0/48 f/p1 roots. First maximum mass is stable to about 2.6e-5 relative from 192 points; first-minimum radius remains unconverged. Actions 42616463263, artifact 13671829050. |
| MD-19a | Initial same-EOS SciPy TOV reference | CONCERN_NUMERICAL_TOLERANCE | 4 native-Hermite EOS tests pass. Five DOP853 vs v5 cases produced max relative mass 3.60164e-7 and radius 1.46572e-6; the 1e-6 radius criterion narrowly failed. |
| MD-19b | Per-case SciPy/native refinement receipt | REPAIR_PUSHED_CI_PENDING | The comparison code now saves per-case coarse/tight SciPy and baseline/tight Fortran results before acceptance. SciPy rtol/max-step and Fortran q-step are refined. 1e-6 remains unchanged. |

The first maximum is nearly converged in mass but not in radius; the following minimum remains badly radius-sensitive. The maxwell-labelled smooth bridge retains its extra local fold pair. These sequence pieces remain kinematic only, and no radial-stability label is assigned.


## Exact-head structural TOV result and reporting repair (2026-11-09)

| ID | Node | Status | Evidence / next condition |
|---|---|---|---|
| MD-19b | Same-EOS independent SciPy TOV acceptance | PASS_NUMERICAL / ACTION_REPORTING_FAILURE | Five-point refined comparison: max relative M=1.6233e-8, R=1.1433e-7, both below unchanged 1e-6. Per-case receipt is preserved in artifact 13673098470. |
| MD-21 | CLI receipt-schema regression repair | IMPLEMENTED_CI_PENDING | The post-compute CLI used obsolete v1 metric keys and raised KeyError; printer now uses the persisted v2 nested schema and a regression test covers the CLI boundary. |
| MD-22 | Six-method mass-pair validator | SCIENCE_HOLD | 0/48 admitted spectra and 0 complete paired cells; first rising-proxy endpoint is about 0.7653 Msun, below the 1.0 Msun target. Do not widen target coverage or extrapolate. |
| MD-23 | Direct radial fundamental-mode stability | NEXT_PHYSICAL_GATE | Implement/verify the nodeless radial eigenvalue omega_0^2 problem with declared perturbative closure and a valid surface condition before allowing modes on secondary rising intervals. |

Run 42617251360 completed 2304/2304 TOV integrations and all frozen-v5 anchor checks. The independent same-interpolation SciPy DOP853 comparison passed its numerical 1e-6 acceptance, but the job status was failure because the CLI summary indexed obsolete keys after writing the passing receipt. The mass-pair validator was skipped by that step failure; even after repairing the reporting bug, it is expected to hold because the selected pilot still has 0/48 modes and no complete 1.0 Msun six-method cell. Overall workflow green is not the exit criterion; physical stability and mode identity remain the controlling gates.


## Surface boundary and radial-stability gate (2026-11-09)

| ID | Node | Status | Evidence / next condition |
|---|---|---|---|
| MD-24 | Shared low-density continuation diagnostic | HISTORICAL_DIAGNOSTIC_ONLY / NO_MODEL_SELECTED | The R13 tangent/curvature tails satisfy a mathematical barotropic closure but are synthetic EOS extensions. They are retained only as historical diagnostics and are excluded from the primary EOS-preserving analysis. |
| MD-25 | Synthetic EOS extension to (P=0) | SUPERSEDED_FOR_PRIMARY_PATH | Do not extend or append EOS rows in the primary study. A physical P=0 continuation is blocked until a sourced low-density crust EOS and a compatible surface boundary are specified. |
| MD-26 | Direct radial fundamental eigenvalue | BLOCKED_ON_PHYSICAL_LOW_DENSITY_EOS_AND_SURFACE_BC | Solve for nodeless (omega_0^2) only after a physical low-density EOS and justified surface boundary condition are available; retain centre regularity, node count and eigenvalue convergence. Initial (Gamma_1=Gamma_{eq}) supports only barotropic-equilibrium stability. |
| MD-27 | Full-star mass-dependent Cowling modes on classified stable segments | BLOCKED_ON_MD-26_AND_PHYSICAL_SURFACE | Require a converged radial-stability label, a physical surface condition, and validated f/p1 root identity before complete-star target-mass scans. |

The surface analysis used the exact six projected EOS CSVs. The first five ((n_b,P,\epsilon,\mu_b)) knots match exactly in all six methods; this permits one shared continuation baseline. However, the local effective polytropic slope changes across those knots, and least-squares exponents depend materially on the fitting window. Hence a broad-window polytrope is not adopted. The C1 tangent polytrope is a starting diagnostic and must be thermodynamically completed and sensitivity-tested before it supplies a radial surface condition.


## Exact-head outcome after CLI repair (2026-11-09)

| ID | Node | Status | Evidence |
|---|---|---|---|
| MD-21 | SciPy TOV receipt-schema CLI repair | PASS_CODE_AND_REGRESSION | Actions 42618432567: 5/5 independent-Hermite/CLI tests pass; CLI prints all v2 receipt groups with no KeyError. |
| MD-22 | Six-method mass-pair validator | SCIENCE_HOLD | Validator intentionally fails: required 1.0 Msun mass has no complete six-method pair; 0/48 roots and 0 complete cells. |
| MD-24 | Shared low-density continuation diagnostic | HISTORICAL_DIAGNOSTIC_ONLY / NO_MODEL_SELECTED | The R13 tangent/curvature tails satisfy a mathematical barotropic closure but are synthetic EOS extensions. They are retained only as historical diagnostics and are excluded from the primary EOS-preserving analysis. |
| MD-25 | Thermodynamic extension to (P=0) | NEXT_EXECUTION_GATE | Enforce first law, (dP=n_b d\mu_b), causality, positivity, and continuation sensitivity before radial surface BC. |
| MD-26 | Direct radial fundamental eigenvalue | BLOCKED_ON_PHYSICAL_LOW_DENSITY_EOS_AND_SURFACE_BC | Solve for nodeless (omega_0^2) only after a physical low-density EOS and justified surface boundary condition are available; retain centre regularity, node count and eigenvalue convergence. Initial (Gamma_1=Gamma_{eq}) supports only barotropic-equilibrium stability. |
| MD-27 | f/p1 on classified stable segments | BLOCKED_ON_MD-26 | No mass-dependent Cowling roots until radial eigenvalue/node and surface-closure gates pass. |

Actions run 42618432567 was executed on code HEAD `2e4dc9982238bd937a92019dd46580b37923998f`. All code and numerical TOV steps passed, including the receipt-schema regression and refined SciPy/native structural comparison. The workflow ended in failure only at the physical 1.0 Msun pairing validator: the first rising-branch proxy ends at approximately 0.7653 Msun, so all 48 requested mass rows are outside it. Artifact 13675450289 has verified SHA-256 `30922c2fafd5260d3c03d6852057530555cd997a64962d5b59d3f8af966af2ee`. A subsequent docs-only commit records the result and the surface-continuation diagnostic; the code tested by the run is unchanged.


## EOS-preserving native-table boundary audit R14 (2026-11-09)

| ID | Node | Status | Evidence / next condition |
|---|---|---|---|
| MD-28 | Native pressure-knot selector and source-integrity guards | PASS_CODE_AND_TESTS | Existing CSV rows are read in source order; pressure must be strictly increasing; cutoff selection returns only exact source-row P, epsilon, n_b and mu_b values; no source file is written. |
| MD-29 | Six-method native finite-pressure TOV cutoff scan | PASS_EXACT_HEAD_48_OF_48 | Actions 42640660436 on exact HEAD `aff4e6ad2789555648f2e564a96bbabd952ee4ac`: 48/48 native TOV integrations, zero failures, all six 323-row source EOS SHA-256 values unchanged, all six row-1 frozen M/R baselines reproduced within 2e-5. |
| MD-30 | Native TOV background-profile preservation | PASS_48_OF_48 | All 48 backend-exported finite-domain background CSVs have verified SHA-256 and row count. They do not extend below the selected positive-pressure endpoint. |
| MD-31 | Physical P=0 completion and full-star radial boundary | BLOCKED_ON_SOURCED_LOW_DENSITY_EOS | No crust EOS, synthetic tail, or vacuum continuation is selected in the primary analysis. A physical P=0 surface and radial boundary condition remain unavailable from the current input tables alone. |
| MD-32 | Finite-domain radial/Cowling diagnostic | DESIGN_DECISION_REQUIRED | May be explored only if explicitly named a truncated-domain eigenproblem with a mathematically justified positive-pressure boundary condition. Do not label its roots as complete-star f/p modes. |

The source cohort is GPPVA_FSU2 x BBB_BHF2 with methods css, cubic_F, cubic_P, maxwell, quintic_F and quintic_P. All six raw projected CSVs have 323 rows and identical first native thermodynamic knots. The minimum pressure is (9.000000079422\times10^{-12}\;\mathrm{MeV\,fm^{-3}}>0), so even the native table-edge endpoint is not the vacuum surface.

The eight selected source rows are 1, 3, 4, 6, 7, 11, 15 and 20 (actual log-pressure offsets 0, 0.311229, 0.467474, 0.779805, 0.936080, 1.560687, 1.994697 and 3.049259 dex). Each method uses its exact source row at the same index and its frozen reference central pressure. No EOS rows are appended, rewritten, reordered or interpolated as new table entries.

At source row 20 relative to row 1, the enclosed mass change is about (-4.54\times10^{-10}) and the radius of that finite-pressure isobar changes by about (-0.2137\%\), approximately (-24.1\,\mathrm{m}). The result says that the mass integral is effectively saturated within this sampled low-pressure part of the supplied tables, while the coordinate location of the selected isobar still moves. It does not estimate the missing outer-crust thickness or the P=0 radius.

The R13 synthetic tail branch is retained for reproducibility but is withdrawn as the primary surface strategy. The legacy synthetic-surface workflow is manual-dispatch only. MD-25 is superseded; MD-26 and MD-27 remain blocked until the physical low-density EOS and boundary issue is resolved. The R14 artifact has ID `13684872040`, size 18,300,545 bytes and SHA-256 `e124339e932636e9e08c92298b2e1269aca3fac4980230083a28a3451421d8f3`.
