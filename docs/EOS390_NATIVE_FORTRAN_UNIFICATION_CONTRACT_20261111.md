# EOS390 Native Fortran Unification Contract

**Version:** 1.0.0  
**Date:** 2026-11-11  
**Repository:** `cosmosapjw-quantum/NS_Cowling`  
**Work branch:** `research/eos390-fortran-unification-v1`  
**Base:** `research/fortran-eos-gap-fill-v1` at `ecf2271e89b68f0a6ef2e6a9c551deba7862f192`  
**Canonical machine DAG:** `data/eos390_native_fortran_unification_dag.v1.json`  
**Model bibliography:** `data/eos390_source_bibliography.v1.json`

## 1. Goal and authority

The goal is to regenerate the original 390 crust-core combinations from the original CompOSE EOS inputs using the repository's native Fortran calculation path.

There is one scientific-computation authority: **Fortran**. Python is a frontend only. It may validate command arguments and paths, serialize/deserialize data, invoke native executables, decode native reports, and present results. It must not perform production EOS interpolation, root finding, EOS reconstruction, composition solving, crust-core matching, causality calculations, TOV/tidal integration, Cowling integration, or batch numerical orchestration.

The historical Python implementation is preserved as read-only legacy/reference material. The canonical historical source in this repository is `eos_generator/compose_eos_pkg/`; the source Dropbox package is `/ns_eos/EoS_matching_code_package/`. `python_frontend/compose_cowling/maxwell.py` and other Python numerical solver modules are also legacy/reference surfaces until replaced by native Fortran-backed frontend calls. Do not silently keep them as production solvers under a new name.

The previous C4 sharp-Maxwell 65-pair DAG is a separate research task and remains unchanged. This contract governs the 390-combination reconstruction/matching matrix.

## 2. Approaches and selected architecture

Three routes were considered:

1. Retain Python/SciPy as a calculation fallback and gradually add native routines. **Rejected.** It creates two authorities for formulas, root selection, composition and stellar integration; a successful Python run would not prove the native Fortran product.
2. Port only the five smooth methods and leave Maxwell/composition as Python reference. **Rejected.** This does not satisfy the user's target and repeats the existing mismatch.
3. Move the full supported legacy feature surface into the existing Fortran code, keep Python as a thin frontend, and preserve Python only as immutable characterization/differential evidence. **Selected.**

Implement incrementally through the single DAG. Do not create a second EOS390 plan. The current repository already contains a Fortran EOS reader/writer, reconstruction and verification modules, five smooth interpolation kernels, a thermo-only sharp Maxwell kernel, and phase-segmented TOV/Cowling code. The migration must connect these components into one production path instead of adding another numerical package.

## 3. Exact target set and frozen baseline

The target cardinality is:

| Axis | Count | IDs |
|---|---:|---|
| Crust | 13 | SPG_M1, SPG_M2, SPG_M3, SPG_M4, SPG_M5, GPPVA_DD2, GPPVA_DDME2, GPPVA_FSU2, GPPVA_FSU2H, GPPVA_FSU2R, GPPVA_NL3wrL55, GPPVA_TM1e, GPPVA_TW |
| Core | 5 | APR_APR, BBB_BHF2, GM_GM1, RG_SLY4, RG_KDE0v1 |
| Method | 6 | maxwell, css, cubic_F, cubic_P, quintic_F, quintic_P |
| **Total** | **390** | **13 × 5 × 6** |

The frozen v3.2 baseline records 390 generated/read-back pair-method outputs, 234 full-table hard-gate passes, 156 diagnostic-only full-table failures, and a causal-domain pass count of 390. The 156 full-table failures are associated with high-density C10 acausality in APR_APR and RG_KDE0v1. Treat those numbers as historical baseline evidence, not new results from this session. Do not clip, truncate, relabel, or smooth an acausal microscopic tail to make it pass.

Never overwrite:

- `/ns_eos/matched_eos_all_combinations`
- `/ns_eos/Crust_Core_Matching_390_Summary.csv`
- `/EoS_core_crust_release_v3_2_20261106`
- any raw or normalized source EOS directory.

The new output must be written only to a new versioned location, initially:

`/home/cosmosapjw/Dropbox/ns_eos/matched_eos_all_combinations_native_fortran_v1_20261111/`

A retry gets a new attempt identifier or new versioned directory; never overwrite an earlier receipt.

## 4. Method semantics: Maxwell is not a smooth matching method

The legacy Python `method="maxwell"` calls `build_maxwell_plateau_bridge`. It constructs a smooth, near-pressure plateau with continuous baryon density. Its output has zero density jump, even though some metadata called it “1st-Order (Maxwell Plateau)”. That label is physically incorrect. The old implementation must remain in the read-only legacy package for characterization only.

The native smooth matching API now has five methods:

- `css`
- `cubic_F`
- `cubic_P`
- `quintic_F`
- `quintic_P`

The old `maxwell` alias is fail-closed in the native `match_crust_core` route until the composition-aware sharp construction is integrated. This is a temporary safety guard, not the target API.

The final public `method="maxwell"` must call the native sharp first-order construction through the same Fortran-owned matching pipeline. It must not call the old Python plateau or recreate a continuous-density bridge.

### 4.1 Sharp Maxwell conditions

For two cold, beta-equilibrated, locally charge-neutral branches with one globally conserved baryon charge:

\[
P_L(\mu_t)=P_H(\mu_t)=P_t,
\qquad
\mu_{B,L}(\mu_t)=\mu_{B,H}(\mu_t)=\mu_t,
\]

\[
n_L=\frac{dP_L}{d\mu_B},\qquad
n_H=\frac{dP_H}{d\mu_B},\qquad n_H>n_L>0,
\]

\[
\epsilon_i=\mu_t n_i-P_t,
\qquad
\Delta\epsilon=\mu_t\Delta n_B.
\]

This is a local-neutrality, one-conserved-charge Maxwell construction. It is not a Gibbs construction with global charge neutrality and not a finite-width mixed phase. The caller/data contract must justify local neutrality on each branch; it may not be inferred from pressure matching alone.

The resulting EOS has continuous pressure and baryon chemical potential at the interface, but a nonzero density and energy-density jump. The serialized density grid is increasing, with the low-side and high-side coexistence rows adjacent, equal pressure/chemical potential, and distinct densities. There must be no artificial sequence of intermediate densities across the jump.

The Seidov rapid-onset diagnostic must be recorded from the actual jump. It is not a substitute for integrating the complete stellar sequence and does not by itself establish a stable branch, maximum mass, twin-star solution, or a realistic EOS claim.

### 4.2 Branch selection and interpolation

The legacy `find_convex_transition_window` and `get_core_start_index` define important source behavior. Port their scientific meaning to Fortran, then test it against all 65 crust-core model pairs before scaling to 390 method runs.

The low-density provider may contain an outer-crust segment, an inner-crust segment, and provider-specific categorical phase changes. The high-density provider may likewise contain multiple composition labels. Do not search the full tables indiscriminately and accept the first pressure crossing.

- Select a thermodynamically eligible high-density crust branch and a homogeneous/appropriate core branch using explicit criteria.
- Keep source thermodynamic rows and composition rows aligned by the complete CompOSE index triple.
- Do not interpolate across an unresolved source grid gap or a true input thermodynamic discontinuity.
- A change of `Iphase` is categorical composition metadata; it is not by itself proof of a Maxwell density jump. Decide whether a boundary is composition-only, a smooth thermodynamic change, a tabulation gap, or a genuine thermodynamic discontinuity using the physical arrays and model datasheet.
- Piecewise interpolation must respect source phase boundaries. Never interpolate an integer `Iphase` as a continuous quantity.
- Use the existing Hermite P(μ) kernel and stationary-point/root isolation only after the correct provider branches have been selected.
- Require exactly the number of admissible coexistence roots declared by the branch-selection contract. Multiple roots are an explicit ambiguous result, not permission to choose one arbitrarily.

For each interpolated phase interval, check the full represented polynomial interval, not just endpoints or a coarse fixed mesh:

\[
n_B>0,\qquad \frac{dn_B}{d\mu_B}>0,\qquad
\mu_B\frac{dn_B}{d\mu_B}-n_B\ge 0.
\]

The last inequality is the causal bound (c_s^2/c^2\le 1) for the one-charge cold barotrope.

## 5. CompOSE and composition contract

Authority: CompOSE Reference Manual v3.01 (2022), arXiv:2203.03209, plus each EOS's own CompOSE data sheet. The supported legacy feature set is narrower than all CompOSE formats and must remain explicit: cold T=0, one-dimensional fixed-((i_T,i_{Y_q})) slice; species IDs 0, 1, 10, 11; at most one nuclear-group quadruple with group ID 1; and Nadd 0 or 1. Other species, multiple group quadruples, multiple thermodynamic slices, or unsupported Nadd layouts must fail explicitly.

### 5.1 Thermodynamic fields

The native reader/writer owns:

- `eos.nb`, `eos.t`, `eos.yq`, `eos.thermo`, optional `eos.compo`
- full source index identity ((i_T,i_{n_B},i_{Y_q}))
- Q1..Q7, Il, Nadd and q_extra
- nucleon masses, fixed-slice metadata, and row order
- source duplicate-index policy and explicit diagnostics.

Do not use the first line of `eos.yq` as a physical proton fraction. For cold beta-equilibrium tables it may be a coordinate placeholder.

### 5.2 Particle and nuclear-group composition

CompOSE `eos.compo` particle IDs:

| ID | Meaning |
|---:|---|
| 0 | electron |
| 1 | muon |
| 10 | neutron |
| 11 | proton |

The supported nuclear group has ID 1 and is stored as ((I,\bar A,\bar Z,Y_I)), where:

\[
Y_I=\frac{n_I}{n_B},\qquad
Y_h=\bar A Y_I.
\]

The current internal representation uses the baryon fraction (Y_h), not the CompOSE number abundance. Serialization must use (Y_I=Y_h/\bar A), not (Y_h).

Required closure checks (using the supported one-group convention):

\[
Y_n+Y_p+Y_h=1,
\qquad
Y_p+\frac{\bar Z}{\bar A}Y_h=Y_e+Y_\mu.
\]

The legacy writer uses a 2e-3 abundance-closure tolerance; the C3/C4 scientific verification tolerance remains 2e-5 unless a separately justified contract changes it. Do not widen the verification tolerance because a generated table fails.

CompOSE v3.01 permits `Npairs=0` and `Nquad=0`. This is useful for a phase-only `eos.compo` sidecar: preserve `Iphase` without inventing particle fractions or a nuclear group. Therefore keep these states distinct:

- `composition_file_present`: whether the sidecar exists.
- `composition_available`: whether usable composition abundances/group values are present.

The first native patch adds phase-only writer support and baryon/charge closure guards. This code change is only static until the exact branch CI and local round-trip tests run.

### 5.3 Iphase and source phase boundaries

`Iphase` is a categorical, model/provider-specific label. Do not assume that 0, 3, 4, 6, or 7 has a universal physical meaning.

Observed source examples from the normalized tables:

- SPG_M1: 365 rows; `Iphase=0` for 217 rows and `Iphase=7` for 148 rows, switching at row 218. The corresponding density grid has a large source gap from approximately (2.69021\times10^{-4}) to (2.51189\times10^{-3}\,\mathrm{fm}^{-3}). Do not interpolate across that gap as a normal smooth EOS interval.
- GPPVA_FSU2H: 300 rows; `Iphase=0` for 217 rows and `Iphase=7` for 83 rows, switching at row 218. The source table has `Nadd=1`.
- RG_SLY4: 197 rows; `Iphase=4` for 143 rows and `Iphase=3` for 54 rows.
- RG_KDE0v1: 1240 rows; `Iphase=4` for 966 rows and `Iphase=3` for 274 rows.

These labels must be preserved and interpreted from each model's datasheet. They are not automatically mechanical first-order jumps. The composite output must define a deterministic provider-local-to-output phase-ID mapping and store it in a provenance report/phase legend.

For smooth matching, use the legacy-comparable composition bridge only after validating it physically. The legacy construction uses the quintic smootherstep

\[
s(u)=6u^5-15u^4+10u^3
\]

and evolves the nuclear-group baryon fraction and proton fraction from the selected crust boundary toward core composition. Preserve source composition rows outside the bridge exactly. The bridge must solve charge/lepton partition and maintain abundance closure. Preserve provider Q4 in unmodified regions and use the documented lepton/lattice correction only in the bridge. This must be ported to Fortran, not executed by Python.

For a sharp Maxwell transition, **do not use the smooth composition bridge**. Preserve the low-phase and high-phase composition independently at their respective coexistence rows. Evaluate each side from its own provider branch without averaging across the density jump.

### 5.4 Il, Q3, Q4 and Q5

Preserve the CompOSE lepton-inclusion flag and model convention.

For Il=0, the stored Q3 potential is charge-shifted. Recover Euler chemical potential using

\[
\mu_{B,\mathrm{Euler}}=(Q_3+1)m_n+Y_q Q_4 m_n.
\]

For Il=1, the supported cold beta-equilibrium Q5≈0 branch maps Q3 to the Euler baryon chemical potential. Keep source/reconstruction algebraic tolerance and sharp-Maxwell Euler consistency tolerance separate.

### 5.5 Nadd

The supported reader/writer accepts Nadd=0 or 1. Preserve q_extra exactly on untouched rows. The 390 source set mixes Nadd=0 and Nadd=1 (SPG_M1 and GPPVA_FSU2H are confirmed examples). Do not assume two providers' q_extra columns mean the same thing just because both have Nadd=1. If their documented semantics cannot be unified into one output column, record an explicit loss/unsupported status and the reason. Never silently set Nadd=0 or relabel a provider-specific flag.

## 6. Native literature-derived composition

For APR_APR and BBB_BHF2, the source CompOSE core tables may not contain `eos.compo`. The legacy Python path has an explicitly named `literature_beta_eq` approximation, which must be ported as a native Fortran module with provenance. It is not the same as reading an original tabulated composition.

The legacy symmetry-energy approximations recorded in its solver are:

- APR: (S(n)=13u^{2/3}+19.6u^{0.62}\) MeV
- BHF2: (S(n)=13u^{2/3}+17.0u^{0.72}\) MeV
- (u=n_B/0.16\,\mathrm{fm}^{-3})

Solve

\[
4[1-2Y_p(n)]S(n)=\mu_e=\mu_\mu
\]

with a free relativistic electron/muon gas and the correct muon threshold; enforce charge neutrality and (Y_n+Y_p=1) in the homogeneous core. Label every output row/report as literature-derived approximate composition, store model and reference keys, and test it against the legacy solver only as a differential oracle. Do not label the result as provider-tabulated composition.

For DD2, DDME2, TW, FSU2, FSU2H, FSU2R, NL3wrL55, TM1e and the SPG family, use the exact provider table where present. Any further microscopic solver or parameter set must be tied to its original paper and independently checked before use.

## 7. Required native code boundaries

Expected native modules (extend existing modules when possible; do not create duplicate competing authorities):

- `compose_eos_types.f90`: the single canonical EOS data type and typed metadata.
- `compose_eos_io.f90`: strict CompOSE read/write, row identity, composition/Nadd preservation, phase-only rows, transactional output or rollback.
- `compose_eos_reconstruct.f90`: first-law master reconstruction and supported outer-crust corrections.
- a native composition/literature module: particle fractions, nuclear-group conversions, beta equilibrium, electron/muon splitting, model parameter tables and provenance.
- `compose_eos_matching.f90`: five smooth methods plus the public `method=maxwell` route to the true sharp construction.
- `compose_eos_phase_transition.f90`: P(μ) interpolation, root topology, Maxwell coexistence, jump output and composition-aware branch sampling.
- `compose_eos_verify.f90`: C1-C11 and composition/phase-edge verification.
- native TOV/tidal/Cowling modules under `oscillation_source/solutions/src`: reuse and integrate existing Fortran implementations; do not retain Python/SciPy as a production fallback.
- `compose_eos_fortran.f90`: manifest-driven batch command(s) and native reports.

Python modules become thin wrappers. Move or quarantine all numerical Python implementations (including the phase-aware `python_frontend/compose_cowling/maxwell.py` reference path) under legacy or test-only namespace before claiming the frontend-only architecture is complete. A frontend wrapper must not contain a hidden fallback to NumPy/SciPy numerical solvers.

## 8. Numerical and scientific validation

Every generated result must keep separate statuses for:

1. scientific/thermodynamic validity,
2. composition validity,
3. numerical interpolation/root validity,
4. causality domain and full-table causality,
5. TOV/tidal/Cowling validity,
6. source provenance,
7. serialization/readback validity.

Required minimum tests:

- CompOSE round-trip for composition-present, composition-absent, phase-only, Nadd=0, Nadd=1, group-present, group-absent, Iphase changes, and source index offsets.
- Negative tests for unsupported species, Nquad>1, Nadd>1, malformed rows, mismatched index keys, nonfinite/negative fractions, abundance closure, missing sidecars, and stale output sidecars.
- Smooth matching: all five methods; native/legacy differential tests; source-row retention; composition closures; Q4 preservation; first-law residual; interval-wide extrema.
- Maxwell: analytic two-phase fixture; no root, one transverse root, tangent root, multiple roots, near-tangent root, shared-knot root, duplicate-endpoint root, interpolation refinement; equal P and μ; positive Δn; Δε=μΔn; nonzero density jump; composition side preservation; Iphase mapping; Nadd policy.
- Fortran-only TOV/tidal/Cowling: analytic smooth limit; small-jump limit; strong-jump interface fixture; slow and rapid conversion boundary conditions; center/surface conditions; convergence with radial/frequency refinement; compare against independently validated reference results.
- 390 output: exactly 390 unique pair-method records; no duplicate or missing keys; each CompOSE output read back; SHA-256 manifest; C1-C11 receipts; diagnostic-only failures retained and never omitted.

A coarse finite sample cannot prove no narrow negative/unstable/ acausal interval. Use the exact polynomial extrema or stationary-point isolation already used by the native kernels wherever available.

## 9. Mutation gates and permissions

- **M0 Source lock:** read and hash sources; never edit raw EOS inputs.
- **M1 Design lock:** define branch selectors, model metadata, phase mapping and Nadd policy before changing calculation formulas.
- **M2 I/O lock:** schema and round-trip tests before writing matched data.
- **M3 Composition lock:** literature approximations are explicitly named and sourced; no guessed fractions.
- **M4 Smooth-method lock:** all five methods verified with composition and first-law closure.
- **M5 Maxwell lock:** no release of `method=maxwell` until true coexistence, nonzero Δn, composition and provenance pass.
- **M6 Stellar lock:** no production stellar observables until Fortran TOV/tidal/Cowling interfaces pass.
- **M7 Batch lock:** no 390 run before M0-M6 pass.
- **M8 Publication lock:** new versioned output only; readback and checksums required.

Allowed repository mutation is limited to `research/eos390-fortran-unification-v1` and its proposed PR into `research/fortran-eos-gap-fill-v1`. Do not merge or mutate `main` without a separate request. Dropbox output mutation is allowed only in the new versioned output folder and progress snapshot. The following are protected: all raw EOS files, v3.2 archive and generated outputs, the frozen legacy Python code and the old 390 summary.

## 10. Anti-drift and anti-hallucination self-audit

Before every substantive commit, the agent must answer these questions in its receipt:

1. **Authority:** Is each formula owned by exactly one Fortran module? Is there any Python/SciPy production fallback?
2. **Physics:** Does `maxwell` have Δn>0 and satisfy equal P, equal μ_B and Δε=μ_tΔn? If not, the result is not a Maxwell transition.
3. **Composition:** Are Y_I and Y_h kept distinct? Are baryon/charge closures checked? Are Iphase labels provider-specific and mapped? Is Q4/Il=0 treated correctly?
4. **Nadd:** Is q_extra semantics known for both branches? If not, is the drop or block explicit?
5. **Provenance:** Can every output row trace to a raw provider row, a specified interpolation interval, or a named derived solver with paper key?
6. **Numerics:** Were root topology and interval extrema checked? Were tolerances fixed before observing the result?
7. **Failures:** Are FAIL, diagnostic-only, blocked, not-run, no-source and unsupported distinct? Were any failed rows clipped or hidden?
8. **Mutation:** Are raw and v3.2 paths untouched? Is the output versioned?
9. **Evidence:** Are build/test commands, environment, exit codes, logs, source SHA and output hashes attached? Did the agent avoid claiming runtime PASS without execution?
10. **Scope:** Did the agent avoid changing the separate C4 65-pair DAG or introducing a second EOS390 DAG?

If any answer is unknown, record it as unresolved/blocked and do not convert it to PASS. Do not invent test output, agent delegation, paper download, root count, or source-row provenance.

## 11. Failure taxonomy

Use distinct machine-readable statuses:

- `THEORY_CONVENTION_ERROR`
- `NUMERICAL_ALGORITHM_FAILURE`
- `IMPLEMENTATION_DEFECT`
- `COMPOSE_DATA_SEMANTIC_ERROR`
- `SOURCE_PROVENANCE_MISSING`
- `UNSUPPORTED_COMPOSE_LAYOUT`
- `RUNTIME_ENVIRONMENT_BLOCKED`
- `AUTHORIZATION_BLOCKED`
- `NOT_ATTEMPTED`
- `DIAGNOSTIC_ONLY`
- `PASS`

A file existing is not proof of correct physics. A CI pass is not a source-cohort admission. A TOV solution is not proof of a stable physical branch unless the EOS and stability criteria are independently validated.

## 12. Local-agent handoff prompt

The following prompt is intended to be copied verbatim into the local coding agent session.

```text
You are the local implementation and execution agent for cosmosapjw-quantum/NS_Cowling.

PRIMARY GOAL
Regenerate the original 390 crust-core EOS pair-method products with the correct Fortran implementation. All scientific computation must be native Fortran. Python is frontend only. The historical Python implementation is read-only legacy/reference code.

CANONICAL CONTRACTS
- Branch: research/eos390-fortran-unification-v1
- Base: research/fortran-eos-gap-fill-v1
- Machine DAG: data/eos390_native_fortran_unification_dag.v1.json
- Full contract: docs/EOS390_NATIVE_FORTRAN_UNIFICATION_CONTRACT_20261111.md
- Bibliography: data/eos390_source_bibliography.v1.json
- Do not change the separate 65-pair sharp-Maxwell C4 DAG.

LOCAL SOURCE DIRECTORIES (ALL THREE MUST BE READABLE)
- /home/cosmosapjw/Dropbox/ns_eos
- /home/cosmosapjw/Dropbox/ns_oscillation_followup_20261106
- /home/cosmosapjw/Dropbox/EoS_core_crust_release_v3_2_20261106

LOCKED INPUTS AND BASELINES
- Raw EOS source root: /home/cosmosapjw/Dropbox/EoS_core_crust_release_v3_2_20261106/EoS_complete_core_crust_release_v3_2/code_package/eos_data
- Normalized EOS root: /home/cosmosapjw/Dropbox/EoS_core_crust_release_v3_2_20261106/EoS_complete_core_crust_release_v3_2/code_package/eos_data_normalized
- Legacy Python source: /home/cosmosapjw/Dropbox/ns_eos/EoS_matching_code_package
- Historical 390 summary: /home/cosmosapjw/Dropbox/ns_eos/Crust_Core_Matching_390_Summary.csv
- Historical matched outputs: /home/cosmosapjw/Dropbox/ns_eos/matched_eos_all_combinations
- Do not modify these inputs, the v3.2 release, the old summary, or the old matched outputs.

EXACT MATRIX
Crust: SPG_M1, SPG_M2, SPG_M3, SPG_M4, SPG_M5, GPPVA_DD2, GPPVA_DDME2, GPPVA_FSU2, GPPVA_FSU2H, GPPVA_FSU2R, GPPVA_NL3wrL55, GPPVA_TM1e, GPPVA_TW.
Core: APR_APR, BBB_BHF2, GM_GM1, RG_SLY4, RG_KDE0v1.
Methods: maxwell, css, cubic_F, cubic_P, quintic_F, quintic_P.
Expected total: exactly 390.

NON-NEGOTIABLE PHYSICS
- The old Python method=maxwell is a fake smooth pressure plateau with zero density jump. It is never a first-order Maxwell construction.
- Final method=maxwell must satisfy equal pressure and baryon chemical potential, positive nonzero density jump, and Delta epsilon = mu_t Delta n_B. It must not create intermediate densities across the jump.
- Five smooth methods remain css, cubic_F, cubic_P, quintic_F, quintic_P.
- The sharp construction assumes one globally conserved baryon charge with local charge neutrality on each phase; do not call it a Gibbs construction.
- Preserve CompOSE Q1..Q7, Il, Iphase, species ids 0/1/10/11, nuclear-group id 1, Abar/Zbar/Y_I, and Nadd/q_extra semantics. Y_h = Abar*Y_I internally; serialize Y_I = Y_h/Abar. Check baryon and charge closure.
- Iphase is model/provider-specific categorical metadata. Preserve it and produce deterministic provider-to-output phase-ID provenance. Do not interpret label changes as a density jump without checking P, epsilon, mu and n.
- Do not interpolate across unresolved source gaps or real thermodynamic jumps.
- For Il=0, recover Euler mu_B by adding Yq*Q4*m_n to the Q3-derived potential.
- APR_APR and BBB_BHF2 composition may be literature-derived; keep this status explicit and use the exact source approximations/references recorded in the contract.
- Never clip or repair acausal source tails by altering original EOS rows. Keep full-table and causal-domain status separate.

PYTHON BOUNDARY
Python can validate arguments, serialize data, invoke native binaries, decode reports and display results. It must not contain production EOS, composition, interpolation, root-finding, TOV, tidal, Cowling or batch numerical logic. Preserve old Python in the legacy tree for differential checks only.

FIRST ACTIONS
1. cd into the repository and verify git status.
2. git fetch origin and checkout research/eos390-fortran-unification-v1.
3. Verify the exact HEAD against the GitHub handoff and do not force-push.
4. Verify the three Dropbox directories above exist and are readable.
5. Hash all raw and normalized input files and the v3.2 baseline outputs before executing anything.
6. Inspect all 18 model READMEs and eos.pdf datasheets. Update per-model source mapping from official CompOSE pages and original papers. Prefer arXiv, HAL, author repositories, open journal PDFs, CompOSE and Zenodo. If a needed original paper cannot be accessed after these routes, report the exact citation and ask the user for that document; do not invent the source.
7. Build and run the current native tests before further modifications; record exact commands, environment, stdout/stderr and return codes.

BUILD / BASELINE COMMANDS
bash fortran_eos/build.sh build/compose_eos_fortran
COMPOSE_EOS_FORTRAN_ROOT=. COMPOSE_EOS_FORTRAN_EXE=build/compose_eos_fortran/compose_eos_fortran PYTHONPATH=python_frontend python -m pytest -q python_frontend/tests pipeline/tests

IMPLEMENTATION ORDER
Follow the canonical DAG M0 through M9:
- Finish strict CompOSE read/write round-trips, including phase-only eos.compo rows Npairs=0,Nquad=0, full composition closures, Nadd and q_extra rules.
- Port all required composition and literature beta-equilibrium logic from legacy to Fortran.
- Port the five smooth methods with composition preservation.
- Route method=maxwell through the sharp native kernel, selected source branches, composition-aware side sampling and provider-local phase remapping.
- Integrate the native phase-segmented TOV/tidal/Cowling path. Do not use SciPy as production fallback.
- Add a Fortran-owned manifest-driven 390 batch runner.
- Run all exact 390 pair-method tasks, write per-task receipts, perform CompOSE readback, C1-C11, full-table and causal-domain reports, and preserve all diagnostic-only failures.

NEW OUTPUT ONLY
Write generated files to:
 /home/cosmosapjw/Dropbox/ns_eos/matched_eos_all_combinations_native_fortran_v1_20261111/
Never overwrite the original v3.2 outputs or source EOS. Use one directory per crust/core/method, a pair manifest, summary.json/csv, provenance.json, failure_catalog.json and SHA-256 manifests.

FAILURE / STOP RULES
- Never claim PASS without the actual command and observed exit code/output.
- If an input is missing, record the exact path and source model. Do not silently substitute a similarly named table.
- If a paper cannot be retrieved, ask the user for that specific paper.
- If a root is absent, tangent, multiple or ambiguous, record that exact topology; do not choose the visually convenient root.
- If Nadd meanings cannot be unified, stop that combination or use an explicitly recorded loss policy. Do not invent a shared meaning.
- Never alter raw input rows, original outputs, main or the C4 DAG.
- Keep theory, numerical, implementation, data-semantic, provenance and runtime failures separate.
- Do not fabricate subagent receipts, downloads, tests, or outputs.

DELIVERABLES
1. Updated Fortran modules and thin Python frontend.
2. Native 390 pair manifest runner.
3. Per-combination CompOSE output and machine-readable validation receipts.
4. Summary of 390 statuses, preserving FAIL/DIAGNOSTIC_ONLY/BLOCKED.
5. Environment, source and output hashes.
6. Updated DAG status and a concise Dropbox progress snapshot.
7. Exact build/test/batch commands, output paths and any paper requests.

The end goal is the 390 source-backed EOS pair-method regeneration, not merely a passing unit-test suite. Do not claim the goal is complete until all 390 records have explicit outcomes and the generated outputs have been read back and hashed.
```

## 13. Current execution state

Initial safety patch prepared on the new branch:

- The native smooth matcher rejects the historical `maxwell` plateau token.
- The Python frontend rejects the same false alias rather than invoking the plateau.
- The five genuine smooth method tests remain separate.
- The native CompOSE writer now has a phase-only serialization path and composition baryon/charge closure checks; it removes a stale `eos.compo` when the output contract says no sidecar.
- Regression tests cover the temporary Maxwell fail-closed behavior and phase-only Iphase round-trip.

**Not yet verified:** no local Fortran build or pytest execution was possible in this ChatGPT host environment. The branch CI result is the next code-verification evidence. The sharp Maxwell kernel is still thermo-only; composition-aware integration into `match_crust_core`, native literature composition, and the complete 390 run remain pending. Do not describe the initial safety patch as a completed EOS390 implementation.
