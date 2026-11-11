# EOS390 Native Fortran Local-Agent Execution DAG and Developer Contract

**Version:** 1.1.0  
**Date:** 2026-11-11  
**Repository:** `cosmosapjw-quantum/NS_Cowling`  
**Branch:** `research/eos390-fortran-unification-v1`  
**Parent before this slice:** `ce24f85cbb4861c28c9f9eea0f3ddea28950fab8`  
**Machine-readable DAG (authority):** `data/eos390_native_fortran_unification_dag.v1.json`  
**Existing 65-pair C4 sharp-Maxwell DAG:** separate and immutable.

## 1. One goal, two execution domains

The end goal is to regenerate the original `/ns_eos` 13-crust x 5-core x 6-method matrix (390 outputs) using the current repository's Fortran code. This is not a request to create a second Python solver. Physical and numerical work belongs in Fortran. Python is restricted to argument/path validation, serialization/deserialization, native executable invocation, report decoding, and presentation. The historical Python EOS package and Python numerical modules are read-only references and characterization sources.

**GitHub/ChatGPT domain:** inspect source, verify literature, design, edit Fortran and contracts, perform source/static checks, add CI tests, commit on the work branch, inspect CI, and preserve evidence.

**Local domain:** verify the actual Dropbox mount, inspect/download source and paper files, run native builds/tests, perform numerical and convergence experiments, run the 390 task manifest, read back each CompOSE output, and upload a new versioned output set to Dropbox. Do not defer work that can be done on the local machine to ChatGPT and do not claim a local run from CI evidence.

## 2. Frozen source matrix and output identity

Crust set (13): SPG_M1, SPG_M2, SPG_M3, SPG_M4, SPG_M5, GPPVA_DD2, GPPVA_DDME2, GPPVA_FSU2, GPPVA_FSU2H, GPPVA_FSU2R, GPPVA_NL3wrL55, GPPVA_TM1e, GPPVA_TW.

Core set (5): APR_APR, BBB_BHF2, GM_GM1, RG_SLY4, RG_KDE0v1.

Historical method axis (6): maxwell, css, cubic_F, cubic_P, quintic_F, quintic_P. Do not silently rename a legacy method while generating a comparison. The future public production meaning of `maxwell` is a sharp Maxwell construction. The old smooth methods remain separately identified; do not use old results as evidence that sharp Maxwell was executed.

The immutable v3.2 benchmark has 390/390 constructed and read back; 234 full-table hard-gate PASS and 156 C10 causality FAIL (all 78 APR_APR and 78 RG_KDE0v1 pair-method outputs). All failed tables are diagnostic-only, not repaired physical EOS. The historical causal-domain subrange may have its own eligibility receipt, but does not change whole-table FAIL.

New local target directory (never overwrite):  
`/home/cosmosapjw/Dropbox/ns_eos/matched_eos_all_combinations_native_fortran_v1_20261111/`

Every attempt must get a unique versioned directory or attempt id. Raw inputs, normalized inputs, historical 390 outputs and v3.2 archives remain unchanged.

## 3. Physical ontology: transition order is not an interpolation method

All calculations use the existing project units: n_B in fm^-3, P and epsilon in MeV fm^-3, mu_B in MeV, and c_s^2/c^2 dimensionless. Keep the (-,+,+,+) metric convention where spacetime equations are involved. In CompOSE cold beta-equilibrium rows, retain the existing Euler/first-law conventions in the canonical contract.

### 3.1 First-order sharp Maxwell interface

For one conserved baryon charge and local electric neutrality in each phase:

```text
P_low(mu_t) = P_high(mu_t) = P_t
mu_B,low = mu_B,high = mu_t
n_low  = dP_low/dmu_B
n_high = dP_high/dmu_B > n_low > 0
epsilon_side = mu_t*n_side - P_t
Delta_epsilon = mu_t*Delta_n_B > 0
```

Pressure and baryon chemical potential are continuous; baryon density and energy density jump. Keep the two coexistence rows adjacent. Never insert a continuous-density interpolation stencil across this interface. This is the locally neutral Maxwell case, not a globally neutral Gibbs mixed phase. Gibbs requires at least (mu_B,mu_Q), charge densities for both phases, and a mixed volume fraction.

### 3.2 Second-order candidate

A continuous density/energy curve is not enough. A second-order classification requires a specified thermodynamic potential and evidence that the potential and first derivatives are continuous while a higher derivative, susceptibility, or order-parameter derivative has a resolved nonanalyticity. Require resolution/refinement robustness and physical provenance. Without these data, use `CONTINUOUS_UNCLASSIFIED`, not `SECOND_ORDER`.

### 3.3 Crossover

Use `CROSSOVER_SMOOTH` only when the relevant thermodynamic potential and derivatives are smooth across a resolved finite-density interval and no thermodynamic nonanalyticity is established. A rapid but smooth change in c_s^2 or a peak in a susceptibility is not alone a first- or second-order transition.

### 3.4 Unresolved interface and phase labels

If Iphase changes, density samples show a gap, or a sharp slope change occurs without enough thermodynamic evidence, preserve the source rows and classify `UNRESOLVED_INTERFACE`. Iphase is a provider-local categorical integer, not a universal physical phase ID. Never sort or linearly interpolate it. When combining providers, emit a deterministic map (provider, source Iphase) -> output Iphase with a phase legend.

A crust-to-medium or medium-to-core interface may carry a density jump. The interface location does not force continuity. The same no-cross-interface rule applies to any intermediate layer. Every interface must be represented in the EOS topology separately from the interpolation law used within each smooth branch.

## 4. CSS: explicit sharp branch versus legacy smooth interpolation

The canonical sharp-interface constant-sound-speed (CSS) construction follows Alford, Han & Prakash (2013), arXiv:1302.4732, and Alford et al. (2015), arXiv:1501.07902. It has an explicit transition pressure, high-side energy-density jump, and constant high-side sound speed. It is not the current v3.2 two-segment continuous bridge.

At a declared interface, take the low-side state (P_t, epsilon_-, n_-), a nonnegative Delta epsilon, and 0<c_s^2<=1. Define

```text
mu_t = (epsilon_- + P_t)/n_-
epsilon_+ = epsilon_- + Delta_epsilon
n_+ = (epsilon_+ + P_t)/mu_t
epsilon(P) = epsilon_+ + (P-P_t)/c_s^2
n(P) = n_+ * [1 + (1+c_s^2)(P-P_t)/
               (c_s^2*(epsilon_+ + P_t))]^(1/(1+c_s^2))
mu_B(P) = (epsilon(P)+P)/n(P)
```

This enforces equal P_t and mu_t across the interface and `Delta_epsilon=mu_t*(n_+-n_-)`. The high-side branch satisfies dP/depsilon=c_s^2 and d epsilon/d n_B=mu_B. A positive jump is a sharp first-order interface. A zero jump is only `CONTINUOUS_UNCLASSIFIED`; it does not prove second-order or crossover. c_s^2=0 is not accepted as a finite-slope CSS branch; values above 1 are rejected.

The interface model is orthogonal to a bridge interpolation choice. CSS may occur above a crust-medium jump or another declared interface; do not force the intermediate layer to be continuous. The native analytic kernel added in this slice is deliberately only a branch-state evaluator. It is not yet wired into the common matcher, CompOSE serializer, TOV/Love/Cowling solver, or 390 runner.

## 5. Name the old false-Maxwell implementation honestly

The old Python `method="maxwell"` function `build_maxwell_plateau_bridge` creates a smooth Hermite pressure bridge. It uses interior knots u=0.22 and 0.78, an inner pressure offset of 0.015*(P2-P1), a 64-point Gauss-Legendre endpoint-integral root solve, and 16-point partial integrals. Its density is continuous and its pressure is not exactly constant; Delta n_B=0. It is neither a sharp Maxwell coexistence construction nor a thermodynamic phase-order classifier.

Literature contains other finite-width replacements/interpolations around a Maxwell crossing, including a polynomial P(mu_B) replacement/interpolation method and the Mixing Interpolation Method (MIM) of A. Ayriyan et al. (2020). Those methods are not algebraically identical to the historical Hermite pressure bridge above. The exact legacy algorithm has no verified source attribution in the inspected repository or searches. Therefore this project assigns the descriptive internal name:

**Legacy Near-Pressure-Plateau Hermite Bridge (NPPHB)**

This is a project-defined descriptive name, not a claim of a published method. Keep the legacy source unchanged. Never expose NPPHB as `maxwell` in the native production API and never include its outputs in the sharp-Maxwell physics validation.

## 6. Generic CompOSE standard-generator upgrade plan

The current project code has a narrow cold, one-dimensional slice reader/writer. The standard generator must be upgraded before admitting all EOS sources:

1. Define a generic in-memory `compose_table_t` with axis arrays and index bounds for T, n_B and Y_q; retain the full (iT,inb,iYq) cell identity. Support 1D, 2D and 3D source tables, plus explicit fixed-(iT,iYq) extraction for EOS390.
2. Parse/write `eos.t`, `eos.nb`, `eos.yq`, `eos.thermo`; model optional `eos.compo`, `eos.micro`, and `eos.mr` as typed sections. Never flatten multidimensional data implicitly.
3. Generalize species and nuclear groups. Do not hard-code only particle IDs {0,1,10,11}, one group, or Nquad<=1 as universal CompOSE restrictions. Retain arbitrary registered IDs and provider declarations.
4. Keep `composition_file_present` distinct from `composition_available`. CompOSE phase-only composition rows with `Npairs=0,Nquad=0` are legal; preserve Iphase without fabricating abundances.
5. Store nuclear group number abundance `Y_I=n_I/n_B`; internal group baryon fraction is `Y_h=Abar*Y_I`. Check `Y_n+Y_p+Y_h=1` and `Y_p+(Zbar/Abar)Y_h=Y_e+Y_mu`. Never blindly reinterpret raw provider abundance columns.
6. Treat Nadd/q_extra as provider-declared semantics. Preserve exact untouched values when meanings agree. Otherwise require a named sidecar or a documented DROP_WITH_RECEIPT operation. Never interpolate an unknown categorical extension.
7. Preserve provider-local Iphase IDs and emit a remapping receipt if outputs combine providers. Never call Iphase=7 "quark" or Iphase=4 "core" without a source definition.
8. Add lossless round-trip fixtures for every source model and independently selected CompOSE examples with each supported grid dimensionality. Unknown mandatory fields fail closed with a specific diagnostic.
9. Only after the generic reader/writer passes its round-trip gate may normalization, literature composition, matching, CSS/Maxwell, or 390 orchestration run on the new canonical object model.

Required output per generated table: source file hashes, provider/model key, exact grid selection, units/schema version, optional-field preservation/loss report, source-local to output-local Iphase map, Nadd policy, solver/code commit, physical gate report, and SHA-256 manifest.

## 7. Canonical machine DAG and mutation gates

The JSON file is the execution authority. The high-level dependencies are:

```text
M0 source lock -> M1 legacy characterization -> M2 native CompOSE round-trip
M2 -> M10 generic multidimensional CompOSE generator
M2 + M10 -> M3 literature composition solvers
M3 + M10 -> M4 five smooth matching models and branch-local composition
M2 + M3 -> M5 sharp Maxwell through the common matcher
M2 -> M5A CSS sharp/high-density branch; M5A is not the common matcher yet
M4 + M5 + M5A -> M6 native TOV/Love/Cowling and interface conditions
M6 -> M7 Fortran manifest-driven 390 tasks -> M8 independent validation -> M9 new output + backup
```

Mutation gates:

- G0: re-read branch, commit SHA, dirty state, and exact tree before changing any file.
- G1: hash the 18 raw EOS folders, 18 normalized EOS folders, original code archives, model PDFs and old 390 outputs. Never mutate them.
- G2: before coding, read the exact module and tests being changed. Do not invent function or type names.
- G3: every physics function has an equation/units/domain contract and an analytic test before connection to the production command.
- G4: build all affected native targets, run the exact test suite and source manifest check. A test not executed is NOT_RUN.
- G5: compare output with independent equations or reference data; no tolerance changes to make a run pass.
- G6: stage outputs under a new attempt path, validate/readback every file, then mark the attempt publishable.
- G7: create a Git commit with parent equal to the inspected branch HEAD and update-ref with expected_sha. Never force-push. Never alter main or the C4 DAG.
- G8: publish a new Dropbox versioned output only after checksums, readback, and failure catalogue are complete.

Stop immediately on an unknown field meaning, paper/source ambiguity, hash mismatch, missing model data, multiple coexistence roots, or disagreement between a stored thermodynamic column and its identity. Record a blocker rather than repairing source data invisibly.

## 8. Anti-drift and hallucination self-audit

Before each commit and at the end of each local run, answer all questions with evidence paths and hashes:

1. What exact commit SHA and tree SHA were read? Did the branch move during this task?
2. Which source files and PDFs were read, with their source paths, revision, DOI/arXiv id, and SHA-256?
3. Is each equation either cited to the primary paper or derived in the report? Are units and assumptions stated?
4. Did I conflate a matching interpolation method, a CSS branch, and a thermodynamic phase transition?
5. Did I ever infer transition order from an interpolation name or Iphase integer?
6. Did any derivative stencil, spline, root bracket, or table iterator cross a declared interface?
7. Did I preserve source data, phase topology, composition, optional CompOSE fields and provider-local Nadd semantics?
8. Did I silently discard a row, species, group, extension, or failing EOS? If yes, STOP and restore/reclassify as blocked.
9. Were all physical calculations done in Fortran? Does Python only invoke native tools and decode reports?
10. What exact shell commands ran, their exit codes, logs, compiler version, floating-point kind, and environment?
11. How many of the 390 unique pair-method keys were attempted, passed, failed, blocked, and not run? Does the key set equal the expected Cartesian product with no duplicates?
12. Can every reported PASS be regenerated from a retained receipt, output hash and log? Otherwise change it to UNVERIFIED.
13. Does any new output overwrite a frozen input or previous attempt? If yes, stop and recover from source backup.
14. Did the final report preserve all failures and uncertainty instead of changing tolerances or redefining the pass criterion?

Required evidence labels: `ESTABLISHED`, `LITERATURE_SUPPORTED`, `DERIVED`, `NUMERICALLY_CHECKED`, `IMPLEMENTATION_VERIFIED`, `CONJECTURAL`, `UNRESOLVED`, `BLOCKED`, `NOT_RUN`. Never use PASS without the exact test/run receipt.

## 9. Local execution order and commands

Start with the exact working directory and verify mounts before copying or running anything:

```bash
cd /path/to/NS_Cowling
git status --short --branch
git rev-parse HEAD
git rev-parse HEAD^{tree}
git diff --exit-code
test -d /home/cosmosapjw/Dropbox/ns_eos
test -d /home/cosmosapjw/Dropbox/EoS_core_crust_release_v3_2_20261106
test -d /home/cosmosapjw/Dropbox/ns_oscillation_followup_20261106
```

Run existing native baseline before feature changes:

```bash
bash fortran_eos/build.sh build/compose_eos_fortran
COMPOSE_EOS_FORTRAN_ROOT=. \
COMPOSE_EOS_FORTRAN_EXE=build/compose_eos_fortran/compose_eos_fortran \
PYTHONPATH=python_frontend python -m pytest -q python_frontend/tests pipeline/tests
gfortran -std=f2018 -O2 -ffree-line-length-none -fimplicit-none -fcheck=all \
  -Ibuild/compose_eos_fortran/mod fortran_eos/tests/test_css_segment.f90 \
  build/compose_eos_fortran/obj/compose_eos_types.o \
  build/compose_eos_fortran/obj/compose_eos_css.o \
  -o build/compose_eos_fortran/test_css_segment
build/compose_eos_fortran/test_css_segment
sha256sum -c SHA256SUMS_v4.txt
```

The current source cohort is not admitted yet. First enumerate all raw and normalized model directories, verify their exact count and names, hash every file, read per-model README/data sheets, and inspect individual source PDFs/references. Retrieve open-source papers from arXiv, journal OA pages or institutional repositories; record URL, DOI/arXiv id, version, downloaded bytes and SHA-256. Ask the researcher only when the model's defining source is paywalled/unavailable and no open equivalent exists.

Do not start 390 regeneration until M0-M8 prerequisites and the exact per-pair manifest are accepted. The runner must create 390 unique keys, write each into a unique staging directory, and read every CompOSE table back with the native parser. It must save all per-pair logs/reports, a machine-readable failure catalog and output hashes. Preserve the v3.2 156 acausal candidates as explicit failures; no clipping or tail removal.

## 10. Required local output and handoff receipt

The local agent must deliver under a new run-id directory:

```text
environment.json
source_manifest.sha256
paper_source_registry.json
build.log
pytest.log
native_unit_tests.log
390_pair_manifest.json
per_pair/<pair_key>/result.json
per_pair/<pair_key>/build.log
per_pair/<pair_key>/readback.json
per_pair/<pair_key>/C1_C11.json
summary.json
failure_catalog.json
phase_topology_and_composition_mapping_receipt.json
legacy_feature_parity_inventory.json
output_manifest.sha256
handoff_result.json
```

`handoff_result.json` must contain: repository, branch, HEAD, parent HEAD, dirty state; source manifest hash; code compiler/version/flags; exact test commands and exit codes; task key-set hash; counts for attempted/PASS/FAIL/BLOCKED/NOT_RUN; all failure IDs; output directory; readback count; output manifest hash; Dropbox upload receipt/revision; explicit statement that raw inputs and old outputs were not modified.

## 11. References and literature status

- Alford, Han & Prakash, *Generic conditions for stable hybrid stars*, Phys. Rev. D 88, 083013 (2013), arXiv:1302.4732. The sharp CSS model uses transition pressure, energy-density jump and high-phase constant sound speed: https://arxiv.org/abs/1302.4732
- Alford et al., *Constraining and applying a generic high-density equation of state*, Phys. Rev. D 92, 083002 (2015), arXiv:1501.07902: https://arxiv.org/abs/1501.07902
- Douchin & Haensel, *A unified equation of state of dense matter and neutron star structure*, A&A 380, 151 (2001), DOI 10.1051/0004-6361:20011402. Their crust-core transition is weakly first order, with a relative density jump of about one percent: https://www.aanda.org/articles/aa/abs/2001/46/aa1755/aa1755.html
- Ayriyan et al., *A mixing interpolation method to mimic pasta phases in compact star matter*, EPJ A 56, 114 (2020), DOI 10.1140/epja/s10050-020-00111-1. This is the distinct MIM scheme, not the legacy Hermite plateau: https://link.springer.com/article/10.1140/epja/s10050-020-00111-1
- Ayriyan, Bastian, Blaschke, Grigorian, Maslov & Voskresensky, Phys. Rev. C 97, 045802 (2018), the parabolic pressure-chemical-potential replacement interpolation discussed by the MIM paper. Do not claim this is identical to NPPHB.
- The historical `build_maxwell_plateau_bridge` inspected in commit `ce24f85cbb4861c28c9f9eea0f3ddea28950fab8` is project-specific in its knot positions, 1.5-percent endpoint pressure difference parameter and Gauss-Legendre integral constraints. No verified published source for these exact implementation choices has been identified.

## 12. Copy-ready local-agent handoff prompt

The following block is the task contract. Execute it literally after checking current GitHub branch state; if a newer commit exists, stop and rebase the instructions onto the newly inspected HEAD rather than overwriting it.

```text
You are the local implementation and execution agent for cosmosapjw-quantum/NS_Cowling.
Your final goal is to regenerate the original Dropbox /ns_eos 390 crust/core/method outputs using the current repository's native Fortran code.

AUTHORITATIVE START:
- Repository: cosmosapjw-quantum/NS_Cowling
- Branch: research/eos390-fortran-unification-v1
- Baseline parent recorded by this handoff: ce24f85cbb4861c28c9f9eea0f3ddea28950fab8
- Canonical machine DAG: data/eos390_native_fortran_unification_dag.v1.json
- Canonical developer contract: docs/EOS390_NATIVE_FORTRAN_LOCAL_AGENT_DAG_20261111.md
- Current implementation slice includes fortran_eos/src/compose_eos_css.f90 and fortran_eos/tests/test_css_segment.f90; inspect actual HEAD and files before assuming they are present.
- Do not modify main, the v3.2 source/archive, /ns_eos historical 390 outputs, raw CompOSE tables, Python legacy, or the separate 65-pair C4 DAG.

NON-NEGOTIABLE ARCHITECTURE:
- All physical/numerical calculations must be in Fortran.
- Python may only validate arguments/paths, serialize/deserialize, invoke native executables, decode reports, and present output.
- Never import Python/SciPy numerical code as a production fallback.
- Preserve failed, blocked, unverified and not-run states. Never promote diagnostic-only results to PASS.
- Never smooth or differentiate across a declared phase interface.

FIRST ACTIONS:
1. Verify git status, branch, HEAD, tree SHA, remote and Dropbox mount paths.
2. Check out the exact current branch; do not assume the recorded baseline SHA is still the remote HEAD.
3. Hash all source, normalized source, source-paper, old 390 and v3.2 files before any mutation.
4. Enumerate the 18 models and per-model README/data sheet/PDF references. Read all EOS-specific paper files and use open-access versions where available. Record URL, DOI/arXiv version and SHA-256.
5. Run baseline Fortran build, current tests, CSS analytic test and sha256sum -c SHA256SUMS_v4.txt.
6. Update only the machine DAG and single developer handoff document when discoveries change scope. Do not make a second plan.

PHYSICAL RULES:
- Sharp first-order Maxwell: common P_t and mu_t, n_high>n_low>0, epsilon=mu_t*n-P, Delta epsilon=mu_t*Delta n_B. State the charge neutrality assumption explicitly (local neutrality for Maxwell; Gibbs/global neutrality needs a separate two-charge solver).
- CSS high branch: epsilon(P)=epsilon_plus+(P-P_t)/cs2; n(P)=n_plus*[1+(1+cs2)(P-P_t)/(cs2*(epsilon_plus+P_t))]**(1/(1+cs2)); mu=(epsilon+P)/n. The interface may be at crust-medium or any source-declared layer boundary.
- For Delta epsilon>0 and valid continuity conditions classify FIRST_ORDER_SHARP.
- For Delta epsilon=0 classify CONTINUOUS_UNCLASSIFIED. Do not guess SECOND_ORDER or CROSSOVER.
- A second-order label needs a specified thermodynamic potential and a resolved higher-derivative/susceptibility nonanalyticity.
- Crossover needs a smooth resolved thermodynamic interval with no established nonanalyticity.
- Preserve provider Iphase values as categorical source-local data; do not infer a universal phase meaning from the integer.
- Keep current legacy css smooth two-segment bridge distinct from the standard sharp CSS branch.
- Name legacy Python build_maxwell_plateau_bridge only as Legacy Near-Pressure-Plateau Hermite Bridge (NPPHB). This is a project-defined descriptive name, not a literature attribution. It is not Maxwell.

GENERIC COMPOSE GENERATOR:
- Implement lossless 1D, 2D and 3D CompOSE axes/cell indices for T, n_B, Y_q.
- Support required eos.t/eos.nb/eos.yq/eos.thermo and typed optional eos.compo/eos.micro/eos.mr.
- Preserve arbitrary species IDs, nuclear group counts/IDs, Nadd extensions and provider provenance. No fixed four-particle/one-group assumptions.
- Y_I=n_I/n_B; Y_h=Abar*Y_I. Validate baryon and charge closure.
- Treat Npairs=0,Nquad=0 as legal phase-only metadata; keep composition_file_present distinct from composition_available.
- Preserve unknown optional fields where possible and fail closed on unsupported mandatory semantics.
- Round-trip all source tables and independent open CompOSE fixtures before applying matching.

IMPLEMENTATION ORDER:
M0 source lock -> M1 legacy characterization -> M2 native I/O -> M10 generic CompOSE generator -> M3 literature composition solvers -> M4 composition-preserving smooth methods -> M5 common sharp Maxwell -> M5A CSS phase segment -> M6 Fortran TOV/Love/Cowling -> M7 exact 390 manifest runner -> M8 independent validation -> M9 new output and backup.
Do not run 390 before M0-M8 prerequisites are admitted.

390 MATRIX:
13 crust models x 5 core models x 6 historical method IDs = 390 exactly.
Every pair-method key must occur once. Keep method identity and physical phase topology separate. Any declared sharp interface is represented by adjacent coexistence states and phase-specific composition; a smooth method is confined to one smooth branch.

ACCEPTANCE:
- Native compiler and flags recorded.
- Build and all tests exit zero.
- Full source hash manifest matches.
- 390 key-set equals Cartesian product exactly, with zero duplicate/missing keys.
- Every output read back by native Fortran.
- All C1-C11 reports and per-pair receipts preserved.
- Whole-table acausality remains FAIL; never clip or trim a source tail to force pass.
- Failed candidate status retained as FAIL or diagnostic-only with reason.
- Output written into a new versioned directory and uploaded only after readback and SHA checks.
- Existing raw/normalized inputs and old 390 outputs remain byte-identical.

ANTI-HALLUCINATION SELF-AUDIT:
Before every commit and final handoff, state exact commit/tree SHA, source hashes, paper identifiers, equations/units, files changed, compiler and commands, exit codes, test counts, all 390 statuses and output hashes. Verify no derivative or interpolator crosses an interface. Check that no transition order was inferred from CSS/cubic/Maxwell names or Iphase labels. If evidence is missing, mark BLOCKED or NOT_RUN; never improvise a citation, value, test, or PASS. No hidden Python fallback. No source data edits. Use a unique staging path; never overwrite earlier output.

REQUIRED OUTPUTS:
environment.json; source_manifest.sha256; paper_source_registry.json; build.log; pytest.log; native_unit_tests.log; 390_pair_manifest.json; per_pair/<key>/result.json; readback.json; C1_C11.json; summary.json; failure_catalog.json; phase_topology_and_composition_mapping_receipt.json; legacy_feature_parity_inventory.json; output_manifest.sha256; handoff_result.json.

STOP CONDITIONS:
- missing/ambiguous paper or data meaning;
- SHA mismatch;
- unknown Nadd/q_extra semantics;
- multiple coexistence roots without specified physical branch rule;
- source or output round-trip loss;
- CI/build/test failures;
- branch HEAD changed unexpectedly;
- any attempt to mutate protected files.
When stopped, preserve the exact blocker and partial evidence. Do not silently change tolerances.

Deliver a concise final handoff with exact commit, test evidence, source admission state, 390 counts, output path, checksums, Dropbox upload receipt, and known remaining gaps.
```


## 13. EOS-specific paper access registry (checked 2026-11-11)

The source-of-record for which paper applies to which EOS is each model's CompOSE `eos.pdf` datasheet, not an assumed family-wide citation. The local agent must archive the exact PDF bytes, DOI/arXiv version, byte count, SHA-256 and source URL in `paper_source_registry.json`.

Open full-text sources found for this matrix include:
- APR: Akmal, Pandharipande & Ravenhall, PRC 58, 1804 (1998), https://arxiv.org/abs/nucl-th/9804027
- BBB_BHF2: Baldo, Bombaci & Burgio, A&A 328, 274 (1997), https://arxiv.org/abs/astro-ph/9707277
- KDE0v1 and SLy4 source discussion: Gulminelli & Raduta (2015), https://arxiv.org/abs/1504.04493
- KDE0v1 parameter paper: Agrawal, Shlomo & Au, PRC 72, 014310 (2005), https://arxiv.org/abs/nucl-th/0505071
- SPG M1-M5: Scurto, Pais & Gulminelli, PRD 109, 103015 (2024), https://arxiv.org/abs/2402.15548
- GPPVA common inner-crust Thomas-Fermi construction: Grill et al., PRC 90, 045803 (2014), https://arxiv.org/abs/1404.2753
- DD2: Typel et al., PRC 81, 015803 (2010), https://arxiv.org/abs/0908.2344
- DD-ME2: Lalazissis et al., PRC 71, 024312 (2005), published version in University of Zagreb repository, https://repozitorij.pmf.unizg.hr/object/pmf:6864/FILE0/download
- FSU2: Chen & Piekarewicz, PRC 90, 044305 (2014), https://arxiv.org/abs/1408.4159
- FSU2R/FSU2H: Negreiros et al., ApJ 863, 104 (2018), https://arxiv.org/abs/1804.00334
- NL3ωρ-L55: Pais & Providência, PRC 94, 015808 (2016), https://arxiv.org/abs/1607.05899
- TM1e: Shen et al., ApJ 891, 148 (2020), https://arxiv.org/abs/2001.10143
- TW-99 open supporting expositions: Avancini et al. (2003), https://arxiv.org/abs/nucl-th/0212080; and https://arxiv.org/abs/0704.0407. These are later TW-99 expositions, **not** the original Typel & Wolter (1999) paper.
- CompOSE reference manual: Typel et al., https://arxiv.org/abs/2203.03209
- Inner-crust CompOSE paper: Malik & Pais, https://arxiv.org/abs/2207.14153
- Tidal density-jump matching reference: Takátsy & Kovács, https://arxiv.org/abs/2007.01139

The GPPVA_TW data sheet found in Dropbox under `/EoS_core_crust_release_v3_2_20261106/.../GPPVA_TW/eos.pdf` is dated 2025-06-26 and explicitly cites three current table-construction papers: C.-J. Xia et al., Commun. Theor. Phys. 74, 095303 (2022); J.-X. Niu et al., arXiv:2506.11492; and C.-J. Xia et al., PRC 105, 045803 (2022), in addition to Typel & Wolter (1999). Do not collapse these into one citation. The data sheet also records seven extra EOS thermo columns (c_e^2,c_s^2,d,Z,A,R_d,R_W); (d) is a categorical pasta geometry code, not a differentiable quantity. The generic generator must preserve these extra fields and their units/semantics.

Exact original PDFs not retrieved here:
- Glendenning & Moszkowski, PRL 67, 2414 (1991), DOI 10.1103/PhysRevLett.67.2414. APS page requires credentials; broad Dropbox and public searches did not find the source paper. Ask the researcher for a PDF if local institutional access fails.
- Typel & Wolter, NPA 656, 331 (1999), DOI 10.1016/S0375-9474(99)00310-3. Later open expositions of TW-99 are available, but not the exact original paper. Ask the researcher for the original if local institutional access fails.
- BPS (1971) and Chabanat et al. (1998) are cited in source data sheets, but their exact article PDFs have not been retrieved. Continue lawful repository/arXiv search; do not substitute an unrelated review silently.

Download policy: public arXiv or institutional PDFs must be downloaded and SHA-256 recorded locally. The connected host in this session could verify their URLs and extracted metadata but could not write the binary PDFs into Dropbox. Do not claim that a PDF has been downloaded merely because its abstract or landing page was fetched.
