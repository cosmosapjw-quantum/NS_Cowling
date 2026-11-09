# NS_Cowling research lineage and handoff receipt

Inventory date: 2026-11-09. Repository: cosmosapjw-quantum/NS_Cowling. Active continuation: research/mass-dependent-cowling-v1. PR #1, #2 and #3 retain their existing states.

## GitHub workstreams

| Workstream | State | Inherited contribution | Not interchangeable with |
|---|---|---|---|
| PR #1: native EOS density-gap fill and smooth crust-core matching | Open, non-draft | Fortran EOS IO, density-gap fill, smooth matching, master-function reconstruction, causality audit and phase-transition construction utilities | Perturbative sharp-interface conversion laws or Cowling spectrum by themselves |
| PR #2: preserve v5 sensitivity results and sharp-Maxwell reference validation | Open, draft, unmerged | Frozen v5 234-row cohort, sharp-interface transfer map, slow/rapid reference tests and adversarial reports | Smooth bridges or fresh mass-dependent spectra |
| PR #3: native Maxwell production TOV/Cowling dispatch | Open, draft, unmerged | Production dispatch for smooth/sharp two-phase backgrounds, slow/rapid roots, limited refinement and independent Python/SciPy references | Historical v5 release backend or 234-row smooth cohort |
| research/mass-dependent-cowling-v1 | Active separate track | Mass-dependent continuation of finite smooth-matching cohort plus fixed-mass response theory/mock experiment | Native sharp first-order Maxwell solver |

The new branch was forked from live native-Maxwell feature head a09ec30d74b363652d2d1aefe990557c468679dd. Fifteen backup/maxwell-native-fortran branches were present in the branch search. No cleanup/destructive mutation was performed.

## Dropbox evidence

| Artifact | Use | Inspection status |
|---|---|---|
| NS_Cowling_REVTeX_v5_20261109.pdf | v5 paper, 234 fixed-mass Cowling roots, EOS projection, node identity, limitations | Text extracted/reviewed |
| ns_oscillation_followup_20261106/PROGRESS.md | 390-to-234 cohort admission, root/mass/surface gates, one-mass scope | Search/folder inventory reviewed |
| FAILURE_ROOT_CAUSE_REPORT_v3_3.pdf | 156 of 390 failures preserve source-core acausality, not matching defects | Text extracted/reviewed |
| NS_Cowling_f_p1_data_tables_v1.zip | Prior numerical tables | Metadata only; not unpacked |
| Integrated_EOS_Cowling_Pipeline_v3.zip | Integrated pipeline | Metadata only; not unpacked |
| NS_Cowling_f_p1_full_release_v1.zip | Public release | Metadata only; 235.86 MB exceeded Dropbox text-read limit |
| NS_Cowling_SOURCE_v5.tar.gz | Compact source distribution | Located by search/listing; not downloaded/hash-verified |
| LECTURE_NOTES.pdf | Earlier solver/equation notes | Not inspected; Dropbox read limit blocked selected large alias |
| FORTIN_CASE2_MATCHING_EVIDENCE.md | Fortin et al. crust-core matching case | Text fetched/reviewed |
| HANDOFF_TOV_COWLING_SOLVER_20260827/02_CONTINUATION_BRIEF.md | Earlier analytic slow-Maxwell project | Text fetched; says physical root refinement NOT_EXECUTED and mode identity NOT_ADMITTED |

Search/metadata discovery is not byte verification. Large ZIP/TAR archives are not source-verified until downloaded, unpacked and checked against hashes/manifests.

## Scientific version boundaries

### v5 smooth cohort
- Starting matrix: 13 crust × 5 candidate cores × 6 methods = 390.
- Admission: 234 retained; 156 source-core acausality cases excluded (78 APR_APR, 78 RG_KDE0v1).
- Crossed cohort: 13 crust ×3 cores (BBB_BHF2, RG_SLY4, GM_GM1) ×6 historical methods =234 rows, 39 fixed pairs.
- Existing spectra are at 1.3300106263 Msun only; spectrum reuse audit says recomputed_modes=0.
- Smooth ell=2 Cowling node identity: f W/V=(0,0), p1=(1,1), with outer-segment checks.
- Historical method label maxwell means a smooth bridge, not a sharp first-order density jump.

### Native sharp-Maxwell track
A true phase boundary with density jump, one-sided EOS values and slow/rapid perturbative conversion. PR #2/#3 belong here. Smooth-cohort W/V node identities do not automatically transfer. August 2026 analytic slow-Maxwell handoff is a distinct early project; its own brief says roots and mode identity were not yet certified.

## v5 source authority

The v5 release compile_manifest.json records 35 source modules under support, methods and calculation, plus calculation/app/cell_solve.f90. The frozen result/input records center_order=2, r0=10 m, rtol=mode_rtol=1e-9, step_cap=0.05, match_index=72, split_index=69, and a finite positive cutoff equal to the first pressure row.

The later solutions refactor is a different implementation and is not the direct v5 comparator. An isolated solve_v5_* frontend now compiles the exact manifest source family. The direct-cell radius limit is an optional namelist argument with historical default 1e6 m; exploratory mass scan explicitly requests 1e8 m and records it. TOV and spectrum use center_order=2.

## Earlier mass-dependent pilot failures preserved

| Actions run | Result | Failure class |
|---|---|---|
| 42598463512 | Failed, zero valid roots | Non-comparable TOV backend (case_driver real_table_star r0=1 m, rtol=1e-8 versus v5 cell_solve r0=10 m, rtol=1e-9, step_cap=0.05) and radius ceiling |
| 42599325312 | Failed before native compilation | Python runner keyword/signature mismatch |
| 42599839757 | Scientific pilot gate failed | Later solutions backend: M/R parity outside gate, historical maxwell-bridge branch and frequency-bracket issues |

These are preserved failure evidence, not physical results. The exact-v5 backend is a repair attempt, not a pass. No new mass-dependent spectrum or Mmax has yet passed.

## Continuation sequence

1. Compile/test exact v5 manifest backend with center_order=2 and configurable radius limit.
2. Rerun six-model pilot BBB_BHF2 × GPPVA_FSU2.
3. Inspect frozen reference TOV M/R parity before reading new roots.
4. Audit every pressure sample and branch label. A single significant decline is a hard unresolved fold candidate; do not bridge a later rise.
5. Require complete six-method pairs at 1.0 and 1.33 Msun, with mass closure, Cowling residual and W/V identity.
6. Only then run full 234×8 and report complete-pair denominators.
7. Keep theorem/mock evidence separate from physical EOS results.

## Local handoff
Read FIXED_MASS_RESPONSE_THEORY_20261109.md, MASS_RESPONSE_MOCK_EXPERIMENT_20261109.md, MASS_DEPENDENT_COWLING_RESEARCH_PLAN_20261109.md, MASS_DEPENDENT_COWLING_DAG_20261109.md and LOCAL_AGENT_RESUME_MASS_DEPENDENT_COWLING_20261109.md.

Native subagent dispatch: NOT_RUN because unavailable in this session. Research used Exa, SciSpace, Parallel, Wolfram, GitHub and Dropbox.


## Additional source-checksum distinction

The v5 release_summary.json records source_zip_sha256:
b2c0dcaf6a5b07956cf8d33c39651d4ade32f44b4c49b33aa08fd84fb458185b

This is a reference checksum embedded in the GitHub release metadata. The Dropbox source archive was not downloaded and independently hashed in this session. The new frontend verifies the 35 manifest paths, their ordering and existence, and derives a runtime build key from the current branch source bytes plus the compile manifest. Its provenance state is therefore MANIFEST_PATH_ORDER_MATCHED / ARCHIVE_BYTE_IDENTITY_UNVERIFIED. Do not describe it as byte-identical to the archived source until the archive is independently retrieved and compared.


## New exact-head research findings (2026-11-09)

### Forward sensitivity and mock verification
- Code/theory commit d1d01e6b114d97cd496640037aafaf7ff7be126e on research/mass-dependent-cowling-v1.
- Fixed-mass IFT identity and a continuous polynomial q=n/mu path gate are documented in FIXED_MASS_RESPONSE_THEORY_20261109.md and MASS_RESPONSE_MOCK_EXPERIMENT_20261109.md.
- Smooth Cowling W,V forward tangent ODE, simple-root derivative, and synthetic DOP853-vs-Brent root refinement are in COWLING_FIRST_VARIATION_RESPONSE_20261109.md.
- Actions 42610037961: mass-sequence helper tests 6 passed; fixed-mass response tests 6 passed; full pilot failed correctly at required 1.0 Msun complete six-method pairing. These mock results are not physical stars or mode predictions.

### 24-point pilot topology
Artifact 13668724637, SHA-256 fb62570658e44f9847d4a08f73174e36f14c091b6960b3b847206232a73b1534.
- Six GPPVA_FSU2 × BBB_BHF2 methods reproduce frozen v5 TOV M/R at reference Pc to near machine precision.
- First sampled M(Pc) peak: 0.743199–0.743409 Msun, radius about 1231–1235 km.
- Subsequent excursion reaches about 0.08 Msun and radii of order 50 km in the coarse scan; the later high-density interval rises to about 1.9147–1.91485 Msun, radius about 9.32 km at the endpoint.
- All six methods show the same qualitative pattern.
- No 0.8–2.0 Msun target lies in the initial rising-proxy branch; zero Cowling spectra are admitted. This is not an Mmax result and not proof of radial stability or its recovery.

### New candidate-segment implementation staged
- extract_candidate_mass_segments divides valid contiguous TOV samples into rising/falling/flat candidate pieces with shared extrema and hard failed/missing-index barriers.
- candidate_mass_segments.csv records candidate pieces, with radial stability explicitly unresolved.
- The original first rising segment is retained only as a conservative mass-root proxy.
- Pilot pressure grid increased from 24 to 96 samples. Exact-head CI and artifact review are pending.
