# EOS Reconstruction and Maxwell Matching: Local-Agent Contract

**Revision:** 2026-11-11  
**Scope:** return to the pre-existing EOS reconstruction / matching work.  
**Selected coding stages:** C5 implementation and C6 software verification only.  
**Authority:** the existing C4 machine DAG. This file is a development contract, not a second DAG.

## 0. Bounded status

The target is the EOS implementation branch, not the separate LA11 analytic-matching branch. Validate the existing reconstruction, all six smooth crust-core matching methods, and the distinct true first-order Maxwell kernel. Implement only defects established by a deterministic regression.

Observed GitHub identities:
- EOS implementation branch: research/fortran-eos-gap-fill-v1, head ecf2271e89b68f0a6ef2e6a9c551deba7862f192.
- PR #1: open, unmerged; base main at 6eed9c5178afd820c9d814e4138dd0f15d1f5a73.
- PR #1 workflow run 42952472519: completed successfully. Tested PR head ecf2271e89b68f0a6ef2e6a9c551deba7862f192 against the base in merge snapshot 9129d53dfc5b4a216ad53cb118c3fc08ab319572.
- Exact run log: native EOS build successful; frontend/pipeline suite 76 passed in 45.77 s; Fortran polynomial-extrema executable reported PASS; SHA256SUMS_v4.txt verification successful; artifact uploaded.
- The existing Maxwell tests establish an analytic single-crossing fixture and negative controls including non-overlapping chemical-potential domains. Dedicated tests for exact tangency, multiple crossings, and duplicate knot/endpoint roots have not been identified. Those root-topology gates remain open.
- The historical smooth method named maxwell is not the physical first-order Maxwell transition. Keep the APIs distinct.

The physical 65-pair source cohort remains blocked by C6B provider/composition semantics. Do not start C4-D8, C4-D9, C4-D10, any pair root search, TOV sequence, stable-branch scan, or Cowling spectrum.

## 1. Single source of truth and source pins

### 1.1 Canonical execution plan

The only machine-readable DAG is:

- data/sharp_maxwell_c4_execution_plan.v1.json
- plan_id: sharp_maxwell_c4_execution_dag_v1
- branch: research/sharp-maxwell-cohort-admission-v1

Update that DAG in place. Do not create another DAG JSON/YAML, task graph, scheduler manifest, or competing master plan. C4-D3R is the bounded validation node; preserve its dependencies on C4-D0..D7 and its downstream relation to C4-D8..D10.

### 1.2 Observed authority pins

| Role | Ref | SHA | Meaning |
|---|---|---|---|
| EOS code | research/fortran-eos-gap-fill-v1 | ecf2271e89b68f0a6ef2e6a9c551deba7862f192 | observed implementation head |
| PR #1 tested merge snapshot | run 42952472519 | 9129d53dfc5b4a216ad53cb118c3fc08ab319572 | object checked by PR CI |
| PR #1 base | main | 6eed9c5178afd820c9d814e4138dd0f15d1f5a73 | base at time of that run |
| C4 plan branch | research/sharp-maxwell-cohort-admission-v1 | 86e2d96cad43cab9ebee5efd36b83750007421df | parent snapshot before this update |
| C4 source authority | research/maxwell-native-fortran-v1 | a09ec30d74b363652d2d1aefe990557c468679dd | source-authority snapshot in C4 |
| legacy provenance | main | 261fcdcbeaf6b5d842480e82e249312e8da5532a | read-only historical products |

These are observed pins, not permanent assumptions. Re-fetch live refs before every mutation. If any branch moved, re-read changed files and record the new identity. Never assume a Dropbox copy is current.

### 1.3 Required source reads

Production:
- fortran_eos/src/compose_eos_reconstruct.f90 — reconstruct_master and density-gap reconstruction.
- fortran_eos/src/compose_eos_matching.f90 — match_crust_core and six smooth methods.
- fortran_eos/src/compose_eos_phase_transition.f90 — construct_maxwell_transition and coexistence root isolation.
- fortran_eos/src/compose_eos_verify.f90 — C1-C11 checks.
- fortran_eos/src/compose_eos_io.f90 — CompOSE reader/writer. Parsing and field semantics are protected.

Tests:
- python_frontend/tests/test_native_reconstruct_master.py
- python_frontend/tests/test_native_crust_core_matching.py
- python_frontend/tests/test_native_maxwell_transition.py
- fortran_eos/tests/test_polynomial_extrema.f90

Notes:
- docs/FORTRAN_EOS_MASTER_RECONSTRUCTION_V1.md
- docs/FORTRAN_EOS_CRUST_CORE_MATCHING_V1.md
- docs/FORTRAN_EOS_MAXWELL_TRANSITION_V1.md

Read current source and tests before editing. Source controls actual implementation behavior. If documentation and source disagree, stop and record the discrepancy.

## 2. Scientific definitions, units, and thermodynamic contracts

### 2.1 Units and CompOSE variables

Use the current repository convention:
- n_B: fm^-3
- P and epsilon: MeV fm^-3
- neutron/proton rest energies and chemical potentials: MeV
- c_s^2/c^2: dimensionless
- no implicit natural-unit conversion or rest-mass convention changes.

For the supported cold table:
\[
P=n_B Q_1,\qquad
\epsilon=n_B m_n(1+Q_7),\qquad
F=\epsilon/n_B=m_n(1+Q_6).
\]
The cold EOS consistency contract includes Q6 approximately equal to Q7 under the existing tolerance. The Euler baryon chemical potential is
\[
\mu_E=\frac{\epsilon+P}{n_B}.
\]
Do not substitute stored Q3 for mu_E without applying the CompOSE Il convention.

### 2.2 Master-function reconstruction

Let x=ln(n_B), F(x)=epsilon/n_B, and Q1=P/n_B. At zero temperature:
\[
dF=Q_1\,d\ln n_B,\quad F'=Q_1,
\]
\[
P=n_BF',\quad \mu_E=F+F',\quad
\frac{c_s^2}{c^2}=\frac{F'+F''}{F+F'}.
\]

Requirements:
1. Reconstruct F(log(n_B)) through the documented PCHIP path and integral correction.
2. Enforce the interval integral integral(Q1 dln(n_B)) = Delta F.
3. Preserve source rows outside the selected gap.
4. Retain the documented cubic-Hermite fallback and compute derivatives analytically.
5. Check finite values, Q1>0, mu_E>0, positive mechanical derivative and 0<c_s^2/c^2<=1 throughout the gap, not just at knots.
6. Never interpolate, smooth, differentiate or apply a derivative stencil across an Iphase discontinuity.
7. Preserve q_extra unless the documented recompute_extra_cs2 option applies.
8. Require explicit gap_phase for a cross-phase fill.
9. Do not claim full parity with the historical EOS generator. Its global smoothing/candidate-selection and microscopic composition generation remain outside this bounded v1 implementation.

An output formula evaluated using the same code path that created the output is not an independent oracle. Use analytic fixtures and separately implemented integral/derivative checks.

### 2.3 Six historical smooth matching methods

The supported identifiers are exactly:
- maxwell
- css
- cubic_F
- cubic_P
- quintic_F
- quintic_P

All six are smooth, single-valued, thermo-only bridge constructions. Historical match_crust_core(method="maxwell") is a compatibility token for a near-pressure-plateau smooth bridge. It has zero density discontinuity, no latent-energy jump, no sharp interface, is_exact_maxwell=false, and seidov_applicable=false.

Do not turn match_crust_core into a sharp-transition API. Do not rename the historical method in-place. True first-order Maxwell is a distinct operation.

The existing differential test compares the native output with the retained legacy reference and checks matching indices, table arrays, retained source knots, monotonicity, positive Q1, causal c_s^2, first-law residual and Q6=Q7. It explicitly records composition dropping. This is an approximation, not complete composition parity. The test skips if benchmark EOS files are absent; report skips separately.

### 2.4 Genuine one-charge Maxwell transition

construct_maxwell_transition(low, high, local_neutrality_confirmed, ...) consumes two already-computed smooth phase branches. It does not derive hadronic or quark matter. The supported physical construction has one globally conserved charge (baryon number), with local electric neutrality imposed separately on both phases and explicitly confirmed by the caller.

At coexistence:
\[
P_L(\mu_t)=P_H(\mu_t)=P_t,\quad
\mu_E^L(\mu_t)=\mu_E^H(\mu_t)=\mu_t.
\]
Within each smooth phase:
\[
n_B=\frac{dP}{d\mu_E},\quad \epsilon=\mu_E n_B-P.
\]
For a density jump:
\[
n_H>n_L,\quad
\Delta\epsilon=\epsilon_H-\epsilon_L
=\mu_t(n_H-n_L)=\mu_t\Delta n_B.
\]

Output uses two adjacent rows with common P_t and mu_t but different n_B and epsilon. There is no finite-density mixed-phase segment. This is not Gibbs. Do not infer that local neutrality is physically selected in real matter; it is an explicit model assumption.

The kernel's current preconditions include T=0, beta equilibrium, Q2 approximately zero, Q5 approximately zero for the supported neutrino-transparent branch, Q6=Q7, matching mass convention, matching Il/Yq convention and thermo-only input. Composition-bearing input is rejected rather than silently stripped. Preserve the explicit caller confirmation of local neutrality.

### 2.5 CompOSE Q3 and separate tolerances

For the supported cold branches:
\[
m_n(1+Q_3)=
\begin{cases}
\mu_E-Y_qm_nQ_4,& I_l=0,\\
\mu_E,& I_l=1,\ Q_5\simeq0.
\end{cases}
\]
For Il=0, recover the Euler coordinate by adding Yq*mn*Q4. Use Euler mu_E for P(mu_E) interpolation and root finding. Preserve the selected CompOSE Q3 convention on output.

Never conflate:
- strict source/reconstruction Q3-Euler residual threshold: 2e-6;
- true-Maxwell v1 kernel input residual threshold: 5e-4.

Do not widen or replace either threshold without the M3 approval gate.

### 2.6 Hermite interpolation, interval extrema, and root topology

True Maxwell uses cubic Hermite pressure interpolation in Euler mu_E. Each interval matches P at both ends and P'=n_B at both ends. Thus n_B=P' is quadratic and dn_B/dmu_E=P'' is affine.

For the supported positive-mu_E, positive-density regime, validate on every complete interval:
\[
n_B>0,\quad dn_B/d\mu_E>0,\quad
\mu_E\,dn_B/d\mu_E-n_B\ge0.
\]
The last inequality encodes c_s^2/c^2<=1. Diagnostic samples cannot replace exact polynomial extrema.

The current root finder creates a union grid from both phase knot grids. On each union interval, P_L-P_H is cubic and its derivative n_L-n_H is quadratic. The source computes stationary points of this derivative, splits the pressure difference into monotone subintervals, tests endpoints, bisects sign-changing intervals and deduplicates roots before requiring exactly one. Preserve this logic unless a regression demonstrates a defect.

Required deterministic topology regressions:
- disjoint mu support;
- overlapping mu support but no pressure root;
- exactly one transverse root;
- exact tangential/even-multiplicity root at a stationary point;
- multiple distinct roots;
- root at a shared union-grid knot counted once;
- duplicate endpoint/stationary root counted once;
- near-tangent non-root and near-coincident endpoint controls at declared tolerance;
- root location/residual convergence under table refinement.

Do not claim these gates pass from code inspection. Every edge case needs a fixture, an asserted status/root count, and an exact-head executed test. Keep NO_ROOT, MULTIPLE_ROOT, invalid input, acausal input and numerical nonconvergence distinct.

## 3. Immutable source and data semantics

Record exact raw-file hashes, model IDs, row counts, table dimensions, selected T/Yq slices, Nadd, q_extra, Iphase boundaries, masses, Il/iT/iYq, entropy, beta-equilibrium metadata and composition fields.

Hard constraints:
1. Raw CompOSE EOS bytes are immutable.
2. The 52 source files and 65 pair identities are protected.
3. No source row may be added, dropped, reordered or rewritten during admission.
4. Never infer missing composition from a column name or a convenient closure identity.
5. Preserve YI, Abar, Zbar, Yh, Yn, Yp, Ye and Ymu as distinct quantities.
6. Do not impose Yh=YI, Yh=Abar*YI, or a phase-dependent repair without primary evidence and separate authorization.
7. Do not change Iphase, Nadd, q_extra, parser/indexing, masses, units, Q3 convention or tolerances under a generic cleanup task.
8. Never send composition-bearing/Nadd!=0 input to a thermo-only kernel by silently stripping fields.
9. Do not differentiate or interpolate across an Iphase boundary.

C6B is PROVIDER_SEMANTICS_UNRESOLVED__BLOCKED. Provider inquiries remain DRAFT_NOT_SENT. Do not send emails, provider tickets or author contact messages.

## 4. Work sequence in the existing DAG

1. Re-fetch live PR #1, PR #3 source authority, PR #4 and their branch refs.
2. Re-read the canonical C4 JSON, this contract, all source modules, test modules and implementation notes.
3. Verify local source identity. If the source is an archive without Git metadata, record commit identity as unavailable; do not infer it.
4. Reproduce analytic reconstruction and the six smooth matching differential fixtures. Record PASS/FAIL/SKIPPED separately.
5. Reproduce true-Maxwell Il=1 and Il=0 analytic fixtures and every existing negative control.
6. Add missing root-topology regressions. Demonstrate RED on unmodified source first. Only then make a minimal justified correction and demonstrate GREEN.
7. Independently verify whole-interval polynomial extrema and the Maxwell latent-energy identity.
8. Run the full EOS frontend/pipeline suite, native Fortran tests and SHA256 manifest check on the exact candidate commit.
9. Update only C4-D3R in the canonical plan. Keep source admission, C6B, C4-D0..D7, C4-D8..D10 separate.
10. Produce the run receipt and stop. Do not merge either PR or modify main.

## 5. TDD and acceptance

Every production fix must follow:
- RED: add a deterministic test; show failure on unmodified code for the intended reason; record command, assertion, exit code and source SHA.
- GREEN: smallest justified implementation change; record exact diff and changed paths.
- REGRESSION: targeted module, all related tests, full EOS frontend/pipeline suite, native Fortran tests and manifest verification.
- REVIEW: inspect for formula, tolerance, units, data, API and unrelated changes.
- RECEIPT: exact commit, commands, environment, inputs, outputs and hashes in C4-D3R.

A collection failure, skipped test, wrong tree or unrelated runtime error is not valid RED evidence.

| Area | Required evidence | Insufficient evidence |
|---|---|---|
| Reconstruction | analytic first-law fixture, integral constraint, row preservation, fallback, whole-gap admissibility | import success or one pointwise formula |
| Six smooth methods | exact-head legacy differential for all six IDs, first-law and knot checks | assuming all parametrized cases ran without reading collection |
| Historical maxwell | zero density jump, is_exact_maxwell=false, seidov_applicable=false | method name |
| True Maxwell | independent Il=0/1 analytic fixtures; common P/mu_E; positive jump; Delta epsilon=mu_t Delta n; C2/output consistency | pressure equality alone |
| Root topology | explicit no-root, transverse, tangent, multiple, shared-knot and dedup controls | source inspection or one transverse root |
| Numerical stability | exact interval extrema and refinement | coarse sample grid |
| Provenance | live SHAs, command, runtime, input/output hashes, full log, manifest | copied old logs/artifact names |
| Source cohort | C4-D0..D7 and C6B independently cleared | synthetic EOS pass |

Use PASS, PARTIAL, BLOCKED, NOT_RUN, SKIPPED and FAILED as separate states. Every status requires a reason and evidence reference.

## 6. Mutation gates M0-M5

- M0 READ_ONLY: inspect/fetch/hash/build/test. Stop if worktree, branch or authority is ambiguous.
- M1 RESEARCH_ONLY: edit only the existing EOS code/tests/docs and the existing C4 plan branch. No second DAG and no main changes.
- M2 TDD_REQUIRED: every production fix needs recorded RED then GREEN. Preserve stdout/stderr and exit codes.
- M3 HUMAN_APPROVAL_REQUIRED: no autonomous changes to equations, units, masses, Il/Yq interpretation, composition/Iphase/Nadd semantics, parser/indexing, row ordering, tolerances, production APIs, solver callbacks, TOV or Cowling.
- M4 PROTECTED: raw CompOSE, 52 inputs, source-authority SHA, 65 pair IDs, legacy outputs and C6/C6B blocker records.
- M5 SEPARATE_AUTHORIZATION: merge, release, deploy, provider contact, email sending, new services, local-device commands outside the selected worktree, Dropbox upload, or other remote mutations without explicit exact authorization.

Publishing documentation on the existing C4 branch does not authorize merging. CI success does not promote blocked scientific gates.

## 7. Anti-drift and hallucination self-audit

Before editing:
- [ ] Fetch live branch refs and compare every stored pin.
- [ ] Confirm canonical DAG path and plan_id; create no parallel graph.
- [ ] Re-read live Fortran/Python source, not just this document or prior summaries.
- [ ] Keep smooth historical maxwell separate from construct_maxwell_transition.
- [ ] Use Euler mu_E as the root coordinate with correct Il conversion.
- [ ] Keep the 2e-6 and 5e-4 tolerances distinct.
- [ ] Verify raw EOS hashes and row counts unchanged.

Before declaring a test:
- [ ] Record exact argv, cwd, environment, commit/tree and input hashes.
- [ ] Record test collection IDs; do not count skipped tests as passes.
- [ ] Show intended RED reason and GREEN result.
- [ ] Give dimensional meaning and explicit absolute/relative tolerances.
- [ ] Use an independent analytic or numerical oracle for load-bearing values.
- [ ] Test root topology and polynomial extrema across intervals, not from plots.
- [ ] Ensure source/output manifests refer to the same bytes actually used.
- [ ] Check filters, early exits, sorting, fallbacks, interpolation and truncation.

Before publishing:
- [ ] Diff contains only authorized paths.
- [ ] No raw EOS, source pin, physics convention, solver API or tolerance changed without M3/M4 authorization.
- [ ] Preserve failure and blocker records.
- [ ] Do not infer PASS from a plan, paper, compile, or model agreement.
- [ ] Plan, source commit, CI run and artifact refer to the same tested snapshot.
- [ ] Keep PRs unmerged and main unchanged.

Any unchecked item is NOT_VERIFIED. Never fill missing evidence with plausible values.

## 8. Failure taxonomy

Select a primary category:
- THEORY_CONVENTION: First Law, chemical potential, units, phase-equilibrium assumptions or thermodynamic identity.
- NUMERICAL_ALGORITHM: interpolation, root isolation, stationary/tangent roots, extrema, conditioning or tolerance.
- IMPLEMENTATION: Fortran/Python logic, bounds, shapes, serialization, API or lifetime.
- DATA_SEMANTICS: CompOSE meaning, composition, Iphase, Nadd, ordering or provider ambiguity.
- PROVENANCE: stale SHA, source mismatch, hashes, artifact/commit mismatch or incomplete test population.
- RUNTIME_ENVIRONMENT: DNS, missing compiler/dependency, unavailable sandbox, offline local host or ClientError.
- AUTHORIZATION: mutation approval or protected-source boundary.

Never classify runtime launch failure as numerical failure. Never classify a synthetic fixture as physical source admission.

## 9. Exact-head CI evidence and runtime limitations

GitHub Actions run 42952472519 used Python 3.13 and gfortran 15.2.0:
- build: bash fortran_eos/build.sh build/compose_eos_fortran
- tests: python -m pytest -q python_frontend/tests pipeline/tests
- result: 76 passed in 45.77 s
- native polynomial extrema: quartic polynomial extrema: PASS
- manifest: sha256sum -c SHA256SUMS_v4.txt, successful
- artifact ID: 13826230310
- artifact size: 280171082 bytes
- artifact SHA-256: 6b4d884929f6ef8e75057969086e1805f91fa5e28c785f1295161d329db7c5ff

Local-session observations:
- sandbox reported Python 3.13.5, Git 2.47.3, gfortran 14.2.0, CMake 3.31.6;
- direct GitHub clone failed due DNS resolution failure for github.com;
- GitHub Actions artifact download into the conversation succeeded;
- subsequent local container/Python execution returned ClientError, so the artifact was not extracted or executed in that sandbox;
- connected local workstation was offline;
- Dropbox text fetch worked; binary transfer to sandbox was blocked by runtime URL-access policy;
- native subagent delegation tools were unavailable. Native subagent execution is NOT_RUN; this session used serial tools, not simulated agents.

These are runtime limitations, not EOS algorithm failures. A local agent with a working shell/network must re-fetch live GitHub state and execute the commands below. Do not reuse temporary Dropbox download URLs.

## 10. Reproducible command contract

Run from repository root on the exact candidate commit:

~~~bash
set -euo pipefail
git status --short --branch
git rev-parse HEAD
python --version
gfortran --version
git diff --check

bash fortran_eos/build.sh build/compose_eos_fortran

PYTHONPATH=python_frontend \
COMPOSE_EOS_FORTRAN_ROOT=. \
COMPOSE_EOS_FORTRAN_EXE=build/compose_eos_fortran/compose_eos_fortran \
python -m pytest -q python_frontend/tests pipeline/tests

gfortran -std=f2018 -O2 -ffree-line-length-none -fimplicit-none -fcheck=all \
  -Ibuild/compose_eos_fortran/mod \
  fortran_eos/tests/test_polynomial_extrema.f90 \
  build/compose_eos_fortran/obj/compose_eos_types.o \
  build/compose_eos_fortran/obj/compose_eos_reconstruct.o \
  -o build/compose_eos_fortran/test_polynomial_extrema

build/compose_eos_fortran/test_polynomial_extrema
sha256sum -c SHA256SUMS_v4.txt
~~~

Do not claim a command ran unless stdout/stderr and the exit code were returned. Diagnose the exact failure; never substitute a guessed command.

## 11. Required local-agent deliverables

1. eos_reconstruction_validation.json
2. smooth_matching_validation.json
3. true_maxwell_kernel_validation.json
4. run_receipt.json with branch/commit, dirty state, versions, argv, cwd, exit codes, runtime, input/output hashes
5. RED/GREEN logs for each production fix
6. updated C4-D3R result in the existing canonical C4 JSON
7. concise change summary, residual blockers and one next action

These are validation receipts, not additional plans. Each JSON must identify schema version, tested commit, source hashes, test IDs, expected/observed quantities, tolerances, status, evidence paths and non-claims.

## 12. Copyable handoff prompt

~~~text
PROJECT: NS_Cowling
TASK: Continue EOS reconstruction and Maxwell matching validation/implementation.
SELECTED STAGES: C5 then C6 only.
MODE: serial fallback permitted. Native subagents were NOT_RUN previously because no delegation tool was available. Do not simulate subagents.

CANONICAL DAG:
  data/sharp_maxwell_c4_execution_plan.v1.json
  plan_id = sharp_maxwell_c4_execution_dag_v1
  branch = research/sharp-maxwell-cohort-admission-v1

COMPANION CONTRACT:
  docs/EOS_RECONSTRUCTION_MAXWELL_LOCAL_AGENT_CONTRACT_20261111.md

LAST OBSERVED PINS (FETCH LIVE REFS FIRST):
  EOS branch research/fortran-eos-gap-fill-v1
  EOS head ecf2271e89b68f0a6ef2e6a9c551deba7862f192
  C4 plan parent 86e2d96cad43cab9ebee5efd36b83750007421df
  C4 source authority a09ec30d74b363652d2d1aefe990557c468679dd
  PR1 tested merge snapshot 9129d53dfc5b4a216ad53cb118c3fc08ab319572

STARTUP:
- Inspect git status, current branch, HEAD, remotes, compiler/runtime.
- Fetch PR #1, PR #3 source authority, PR #4 and their current refs.
- Read the canonical C4 JSON and this contract.
- Read native reconstruct, smooth matching, true Maxwell and consistency modules.
- Read the three matching/reconstruction test modules and the implementation notes.
- Verify source identity. An archive without Git metadata does not prove its commit.

PHYSICS:
- n_B fm^-3; P and epsilon MeV fm^-3; mn, mp and mu MeV.
- P=n_B*Q1; epsilon=n_B*mn*(1+Q7); F=epsilon/n_B; mu_E=(epsilon+P)/n_B.
- For x=ln(n_B): F'=Q1; P=n_B*F'; mu_E=F+F'; cs2/c2=(F'+F'')/(F+F').
- Preserve integral(Q1 dln(n_B))=Delta F and source rows outside the gap.
- Never interpolate/differentiate across Iphase boundaries.
- Smooth method IDs: maxwell, css, cubic_F, cubic_P, quintic_F, quintic_P.
- Historical match_crust_core(method=maxwell) is smooth only: zero density jump; is_exact_maxwell=false; seidov_applicable=false.
- True first-order transition is separate construct_maxwell_transition on precomputed branches. It assumes one globally conserved baryon charge and requires explicit local-neutrality confirmation.
- At coexistence P_L=P_H=P_t and Euler mu_L=mu_H=mu_t; n=dP/dmu_E; epsilon=mu_E*n-P; Delta epsilon=mu_t*(n_H-n_L), n_H>n_L.
- Root coordinate is Euler mu_E, never raw Q3.
- Il=0: mn*(1+Q3)=mu_E-Yq*mn*Q4; recover mu_E by adding Yq*mn*Q4.
- Il=1: Q3=mu_E only for the supported cold Q5~0 branch.
- Keep strict source/reconstruction tolerance 2e-6 separate from true-Maxwell local-kernel tolerance 5e-4.
- Do not claim local neutrality, Gibbs equivalence, provider composition semantics, Seidov stability or physical cohort admission from a synthetic fixture.

IMPLEMENT:
A. Validate reconstruct_master PCHIP/integral/Hermite paths, source-row retention, phase-boundary guard, q_extra and full-gap admissibility.
B. Replay all six smooth methods against legacy. Confirm historical maxwell has no density jump.
C. Validate true-Maxwell Il=0 and Il=1 fixtures, Q3 conversion, preconditions, latent-energy identity and output consistency.
D. Add deterministic root-topology tests: disjoint support; overlapping but no root; one transverse root; exact tangential/even-multiplicity root; multiple roots; shared-knot root; duplicate endpoint/stationary root; near-tangent non-root; refinement convergence. The current implementation already splits at stationary points; test before changing it.
E. Check whole-interval n>0, dn/dmu>0 and mu*dn/dmu-n>=0 using exact extrema, not point samples.

TDD:
- For each production correction, add a deterministic regression first.
- Show RED on the unmodified source for the intended reason; capture SHA, argv, output and exit code.
- Make the smallest fix; show GREEN; run targeted and complete regression gates.
- A collection failure, skip, wrong tree or unrelated runtime error is not valid RED.

MUTATION:
M0 read-only until source identity is verified.
M1 existing EOS code/tests/docs and existing C4 plan only. No new DAG, no main changes.
M2 RED/GREEN required for every production fix.
M3 human approval for equations, units, masses, Il/Yq, composition/Iphase/Nadd, parser/indexing, row order, tolerance, API, TOV/Cowling changes.
M4 raw EOS bytes, 52 source inputs, source SHA, 65 pair IDs and legacy products immutable.
M5 separate authorization for merge/release/deploy/provider contact/email/new service/local external device/Dropbox mutation.

C6B BLOCKER:
Provider composition semantics remain unresolved. Do not change raw fields or reader/verifier equations. Keep YI, Abar, Zbar, Yh, Yn, Yp, Ye and Ymu distinct. Do not impose Yh=YI or Yh=Abar*YI. Provider queries are DRAFT_NOT_SENT. Do not contact anyone.
Do not start C4-D8/D9/D10, 65-pair roots, TOV or Cowling until C4-D0..D7, C4-D3R and C6B are separately cleared.

COMMANDS (repo root, exact candidate commit):
  bash fortran_eos/build.sh build/compose_eos_fortran
  PYTHONPATH=python_frontend COMPOSE_EOS_FORTRAN_ROOT=. COMPOSE_EOS_FORTRAN_EXE=build/compose_eos_fortran/compose_eos_fortran python -m pytest -q python_frontend/tests pipeline/tests
  gfortran -std=f2018 -O2 -ffree-line-length-none -fimplicit-none -fcheck=all -Ibuild/compose_eos_fortran/mod fortran_eos/tests/test_polynomial_extrema.f90 build/compose_eos_fortran/obj/compose_eos_types.o build/compose_eos_fortran/obj/compose_eos_reconstruct.o -o build/compose_eos_fortran/test_polynomial_extrema
  build/compose_eos_fortran/test_polynomial_extrema
  sha256sum -c SHA256SUMS_v4.txt
Capture exact argv, cwd, versions, source/input/output hashes, full logs and exit codes. Do not claim a run if the runtime did not return.

OUTPUTS:
  eos_reconstruction_validation.json
  smooth_matching_validation.json
  true_maxwell_kernel_validation.json
  run_receipt.json and full logs/hashes
  RED/GREEN regression evidence
  updated C4-D3R in the existing canonical C4 JSON
Use PASS/PARTIAL/BLOCKED/NOT_RUN/SKIPPED/FAILED separately. Unknown is never PASS.

STOP:
- No C4-D8, C4-D9 or C4-D10.
- No 65-pair source root search, TOV sequence or Cowling spectrum.
- Do not merge PR #1/#4; do not modify main.
- Report changed files, exact tested commit, CI/job IDs, test counts including skips, remaining blockers and one next action.
~~~

## 13. Literature anchors

These sources support the physical distinction, not the interpretation of project-specific CompOSE columns:
- L. T. T. Soethe, “Phase Transitions in Neutron Stars” (2019), https://arxiv.org/html/1901.03252
- A. Bhattacharyya, I. N. Mishustin and W. Greiner, “Deconfinement Phase Transition in Compact Stars: Maxwell vs. Gibbs Construction of the Mixed Phase” (2009/2010), https://arxiv.org/abs/0905.0352
- M. G. Alford, S. Han and M. Prakash, “Generic conditions for stable hybrid stars” (2013), https://arxiv.org/abs/1302.4732

Maxwell is a choice of conserved-charge and neutrality assumptions. The Seidov condition is an onset-stability criterion under its stated assumptions, not a substitute for a full stellar stability analysis.

## 14. Reporting contract

Report exact source/plan identities, changed files, mutation gate, argv/cwd/versions/exit codes, pass/fail/skip counts, independently checked identities and residuals, hashes, artifact locations, residual blockers and one next action.

If blocked, identify whether the blocker is theory, numerical algorithm, implementation, data semantics, provenance, runtime or authorization. Never imply that local execution, publication, or Dropbox backup occurred without a successful tool result.
