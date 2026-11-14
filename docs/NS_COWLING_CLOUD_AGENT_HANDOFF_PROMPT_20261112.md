ROLE
You are the cloud coding agent continuing the NS_Cowling research repository. Follow this exact contract and machine DAG. Do not infer project state from chat summaries or old PR descriptions.

REPOSITORY AND LIVE SOURCE PINS
Repository: https://github.com/cosmosapjw-quantum/NS_Cowling
Cloud EOS branch: cloud/eos390-fortran-unification-v1
EOS PR: https://github.com/cosmosapjw-quantum/NS_Cowling/pull/13
Latest source pin at handoff authoring: 0cbe52ac1ace0535dea541beae31a49410abdbe5
Sharp Cowling PR: https://github.com/cosmosapjw-quantum/NS_Cowling/pull/3
Sharp branch head: 99a0f9c314cab996c3a0bac28299e7b604b6cd01
Machine DAG: data/ns_cowling_cloud_agent_dag.v1.json
Detailed contract: docs/NS_COWLING_CLOUD_AGENT_DAG_20261112.md
Audit report: docs/NS_COWLING_AUDIT_REPORT_20261112.md
Re-fetch PR heads, branch refs, tree hashes, workflows and local checkout state before any edit. The listed SHA values are starting pins, not permission to overwrite later changes.

MANDATORY GOAL
The end state remains the original 13 crust x 5 core x 6 historical-method matrix (390 unique keys), regenerated with a Fortran numerical backend and Python as a thin frontend/orchestrator. Separately, validate the sharp-Maxwell relativistic Cowling interface path and its mode identity. The two workstreams share provenance and physics contracts but are not interchangeable.

CRITICAL METHOD IDENTITY
The historical 390 matrix token maxwell means the project-defined smooth NPPHB near-pressure-plateau Hermite bridge. It is not a sharp first-order Maxwell transition. Do not dispatch the historical token to the sharp M5 oracle and do not relabel it sharp_first_order. A genuine sharp Maxwell transition is a separate operation; if added to an expanded matrix, use a new token sharp_maxwell and a new versioned manifest. Do not silently modify the historical 390-key matrix.
cubic_F/P and quintic_F/P are interpolation method identifiers, not physical classifications. Smoothness alone does not prove a crossover.
CSS needs an explicit source-backed P_t, epsilon_low, n_low, Delta_epsilon and cs2/c2, plus native full-table construction and read-back. Never hardcode Delta_epsilon=0.5 epsilon_low or cs2=1/3 as a production prescription.

CURRENT AUDIT FINDING
The historical M7 summary has 390 receipts: 150 PASS, 30 FAIL, 210 BLOCKED. Source code showed all 150 PASS are invalid as complete native combined outputs: 30 CSS rows are core-only EOS plus a sidecar; 120 cubic/quintic rows are crust-copy stubs. M9 reported 26 PASS/364 FAIL, inconsistent with M7. M9's old 5/5 PASS logic used a native reader read-back as if it proved native generation, treated historical unchanged as true by construction, and called local JSON creation a Dropbox upload. Preserve all files but mark the old release claim superseded. A corrected adjudication must be produced from exact raw receipts.

IMPLEMENTATION POLICY
Fortran performs all production numerical EOS operations, interpolation, coexistence root finding, thermodynamic checks, CSS evaluation, TOV/Cowling ODEs, interface maps, mode determinants and roots.
Python is limited to arguments/path checks, strict serializers, native executable invocation, return-code/report decoding, hashes and presentation. Python M5/M5A oracles are characterization tools only, never a silent fallback.
A PASS receipt requires native_fortran backend identity, executable SHA-256, compiler and flags, complete_combined_eos=true, compose_readback_pass=true, physics_contract_id, exact input source hashes, output SHA-256 and C1-C11 report. Missing evidence means BLOCKED/PARTIAL, never PASS.

EOS AND COMPOSITION RULES
CompOSE v3.01 eos.compo rows use grid key, provider-defined Iphase, Npairs particle index/fraction pairs, and Nquad nuclear-group tuples. Npairs=0 and Nquad=0 are legal; relevant species can be omitted. Preserve row/species/group presence separately from values. Explicit zero is present. Do not infer zero from an omitted field unless a source/model contract explicitly derives it.
C3 baryon closure and C4 charge closure are independent per-row checks. CHECKED, PARTIAL, NOT_AVAILABLE and FAIL stay distinct. No dropping composition files, no mapping ambiguous SPG/GPPVA group fields to YI/Xh/Ye/gas fractions, no fabricated particle fractions, no tolerance inflation.
Iphase and Nadd/q_extra are provider-specific categorical/typed metadata, never generic interpolated thermodynamic values. GPPVA_FSU2H q_extra includes categorical pasta geometry and is not Iphase.
M10 must retain all T/nB/Yq axes, exact (iT,inb,iYq) keys, optional files, species/groups, Nadd and metadata. No flattening, nearest-index guessing or silent field loss.

SHARP MAXWELL COWLING
Inner = high-density side, outer = low-density side. Preserve both one-sided interface states. No interpolation across a discontinuity.
For cold same-coefficient 2x2 Cowling states y=(W,V), D=W1 V2 - V1 W2 and trace-first Liouville give Omega=h exp(-nu) D only under the stated barotropic first-law/TOV assumptions and after deriving the actual source RHS and q-to-r chain rule. This is an algebraic Abel-Liouville invariant, not a claimed canonical covariant symplectic current.
Physical slow/material and rapid/equilibrium endpoint maps remain separate. The LA10 theta interpolation is algebraic kinematic boundary interpolation only, not finite-rate kinetics. No flux, affinity, dissipation or calibrated rate has been derived.
LA11 is a finite-interval BVP matching determinant, not automatically a whole-line Evans function. Retain determinant phase, endpoint basis, gauge, column scales and interface chart factors. Do not apply argument principle or multiplicity claims until holomorphy and Fredholm hypotheses are checked on the actual production system.
Do not label roots f/p1/interface-g by frequency ordering. Use phase-local node counts, profiles, overlap/polarization and parameter continuation.

WORK TO PERFORM IN THIS CLOUD TASK
1. Verify exact current refs, check run IDs and audit claims; update run_identity and provenance receipts.
2. Complete E390-02A mask work, sparse eos.compo reader/writer tests and native C3/C4 checks. Ensure explicit-zero/absent semantics are lossless.
3. Complete E390-10 multi-D CompOSE schema/serialization design and Fortran implementation; add lossless round-trip fixtures.
4. Port NPPHB cubic/quintic F/P kernels to Fortran (M2/M10) with frozen reference inputs, analytic derivatives, edge/extrema tests and Python-reference parity receipts.
5. Complete native CSS parameter ingestion and native full crust-core output construction; do not make guessed parameters.
6. Make M7 manifest-driven batch invoke only Fortran. Legacy maxwell maps to native NPPHB. CSS is blocked until source parameters and native joined output exist. Do not mark placeholders PASS.
7. Correct M9 acceptance: cardinality, genuine native output/read-back, SHA hashes, historical immutability against a trusted baseline and actual Dropbox upload acknowledgement are distinct gates. Keep old M9 report unchanged and add a corrected superseding receipt.
8. Integrate PR #3 sharp-Maxwell Cowling work only on a separate integration branch after exact source/semantic diff. Do not merge PRs automatically.
9. Complete phase-resolved one-sided profile export, trace-first invariant from actual RHS, slow/rapid interface checks, finite-interval determinant and mode identity. Synthetic checks remain scoped to synthetic evidence.
10. Produce fresh exact-head CI and failure catalogue. Stop at blockers; no production 390 batch until source/provider semantics and all native kernel gates close.

MUTATION GATES
G0 fetch live refs/tree/status; stop on concurrent drift.
G1 use an isolated branch/worktree; never reset or delete dirty work.
G2 raw EOS, original paper PDFs, historical 390 outputs and v3.2 reference are read-only.
G3 write a focused failing regression test first; observe the intended failure.
G4 one semantic change per patch, no unrelated refactor.
G5 native Fortran compile/unit tests and full Python integration suite.
G6 independent analytic/symbolic checks and exact Python/reference parity where appropriate.
G7 manifest 390/390 unique keys; never collapse missing tasks.
G8 native read-back all claimed outputs and run C1-C11.
G9 compare old data to exact baseline hashes; no "unchanged by construction".
G10 local Dropbox execution and cloud CI are different environments and claims.
G11 no force push, merge, release, raw data edits or unacknowledged Dropbox mutation.

REQUIRED COMMANDS
Build:
  bash fortran_eos/build.sh build/compose_eos_fortran
Full tests:
  PYTHONPATH=python_frontend COMPOSE_EOS_FORTRAN_ROOT=. COMPOSE_EOS_FORTRAN_EXE=build/compose_eos_fortran/compose_eos_fortran python -m pytest -q python_frontend/tests pipeline/tests
Static checks:
  python -m compileall -q python_frontend scripts
  python - <<'PY'
  import ast, json
  from pathlib import Path
  for p in Path("python_frontend").rglob("*.py"):
      ast.parse(p.read_text(encoding="utf-8"), filename=str(p))
  for p in [Path("data/ns_cowling_cloud_agent_dag.v1.json"),
            Path("data/eos390_native_fortran_unification_dag.v1.json")]:
      json.loads(p.read_text(encoding="utf-8"))
  PY
Also run every native Fortran test target in the workflows. Record compiler, flags, stdout/stderr, exit status, collection count and source/output hashes. Do not assume these commands passed until their exact run logs show it.

PROTECTED PATHS
Do not mutate Dropbox/raw data, upload original copyrighted PDFs to GitHub, or overwrite download/m7_batch, download/m8, download/m9, historical v3.2, or previous 390 outputs. Use a new unique versioned output path. Cite papers in comments and bibliography with DOI/arXiv URLs; do not paste copyrighted full text.

ANTI-DRIFT SELF-AUDIT BEFORE EACH CHECKPOINT
Record exact HEAD/tree and changed files; source file blob/hash; input/paper source hashes; method token semantics; physics equations and units; row-level masks; native/Python execution boundary; all test counts and logs; each of 390 task statuses; M9 read-back/immutability/upload receipts; one-sided interface values; actual RHS/chain-rule derivation; determinant basis/scale factors; mode identity evidence; all blocked/failed controls. Unknown or inaccessible is NOT_RUN/BLOCKED, never PASS.

FINAL ACCEPTANCE
Do not claim all 390 regenerated unless exactly 390 unique, correctly interpreted tasks have complete native Fortran outputs, every CompOSE output reads back, C1-C11 pass or have justified scope status, all source and output hashes are recorded, historical data is unchanged against a trusted baseline, and Dropbox confirms the new versioned upload. If a provider semantic blocker remains, preserve it and report exact affected tasks. Keep moving on independent native code work without bypassing the blocker.
