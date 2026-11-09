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
