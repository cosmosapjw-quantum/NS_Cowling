> Superseded for current work by [2026-11-16 handoff](CORE_CRUST_REAUDIT_HANDOFF_20261116.md). Historical content below is preserved.

# NS_Cowling: resume after physical normalization and Maxwell/admission repair

Repository: https://github.com/cosmosapjw-quantum/NS_Cowling
Branch: remediation/2026-11-15-m5a-maxwell-eos-contracts
PR:15, keep OPEN/DRAFT/UNMERGED. No force-push or history rewriting.

Read actual remote HEAD first. Baseline for this audit was
`d90d3dca033f759c79ed034b3b8ff939acb43145`; older `6b28fadc...` is historical.
Current authority:
- docs/CORE_CRUST_ADVERSARIAL_REAUDIT_20261115.md
- data/receipts/core_crust_audit_20261115/AUDIT_RECEIPT.json
- docs/MASS_NORMALIZATION_CORRECTION_20261115.md
- data/remediation/MASS_NORMALIZATION_CONTRACT.json
- data/receipts/core_crust_audit_20261115/CENSUS_b11a530.json
- latest PR CI and core-crust-census artifact matching current head/merge SHA.

The user's existing /home/cosmosapjw/Dropbox/NS_Cowling checkout was dirty at
9179edda... with separate generator, native and DAG work. This audit used an
isolated clone. Inspect and preserve that uncommitted work before integrating;
do not overwrite it with this branch. Review overlapping files explicitly.

Completed repairs to preserve:
1. Normalization changes Q representation only. Physical crust F,epsilon,n,P
   and mu_Q stay fixed. Applied dm_mass_gauge_meV/dm_crust=0; reference delta is
   diagnostic only. Do not reintroduce the historical additive energy shift.
2. Native Maxwell exact-midpoint bisection fix and phase-only eos.compo output.
   Iphase preserved, no abundance inventory invented, correct1x1 axes/masks.
3. Python reference Maxwell uses union-knot quartic pressure isolation rather
   than a fixed scan; multiple roots require explicit selection and carry a
   selected-interface scope, not global phase dominance.
4. Strict native T/Yq/Q5 admission and independent row/report validation.
5. Direct pair/census outputs are non-destructive and detailed verifier JSON
   must agree with summary status. CI runs and archives the complete census.

Executed before publication: generator91pass; native/frontend/pipeline357pass,
9skip; new focused native tests29pass. Re-read final CI for authoritative head.
The corrected exploratory census remains80configured/60evaluated/20missing,
10CSS survivors(6 distinct),0Maxwell. Use committed-code receipt/artifact for
exact numbers and output hashes; do not reuse pre-correction numerical outputs.

Next bounded work, in order:
A. Reconcile dirty local work against live branch without raw-source changes.
   Confirm CI head SHA versus tested PR merge SHA; read failures/skips honestly.
B. Replace sampled branchwise first-law/causality acceptance with extrema of
   the declared density/potential interpolant. Test narrow acausal pockets,
   tangency, positive compressibility and one-sided jump boundaries. Secants
   do not certify a higher-order interpolant. Keep physical tolerances fixed.
C. Adjudicate source-backed physical energy conventions for Maxwell inputs.
   Current unequal-mass gate remains:30 attempts are implementation/convention
   blocked, not physics rejections. Explicitly distinguish normalization from
   modifying constituent masses or the Hamiltonian. Preserve data provenance.
D. Resolve source/reconstruction issues: missing GPPVA_DDME2/TW grid files,
   RG beta-equilibrium, duplicate reconstructed SPG_M1/M2/M4, raw SPG convexity,
   APR high-density acausality. Do not manufacture source axes or weaken gates.
E. Implement abundance-aware Maxwell v2 and a generated CSS→phase-segmented
   TOV/Cowling adapter as separate bounded tasks. Current phase-only output
   does not prove abundance preservation or elastic-crust physics.

Verification commands (native and reference packages share a module name;
run their suites in separate processes with the stated PYTHONPATH):
  bash fortran_eos/build.sh build/compose_eos_fortran
  PYTHONPATH=python_frontend COMPOSE_EOS_FORTRAN_EXE=$PWD/build/compose_eos_fortran/compose_eos_fortran python -m pytest -q python_frontend/tests pipeline/tests
  PYTHONPATH=eos_generator:eos_generator/compose_eos_pkg python -m pytest -q eos_generator/tests
  NS_COWLING_MATCH_OUT_ROOT=download/matched_runs/<NEW_UNIQUE_NAME> python scripts/match_core_crust.py
  sha256sum -c SHA256SUMS_v4.txt
  (cd eos_generator && sha256sum -c SHA256SUMS_v3_3.txt)

For each new defect: show the failing counterexample, implement the smallest
correction, execute focused/broader checks, preserve actual receipts and update
hashes. Commit/push bounded changes; confirm exact-head CI. Keep failed physics,
numerical, implementation, environment and unavailable-source cases separate.
Report implementation, numerical, native-readback and scientific scope
separately. No global completion claim from a green workflow.
