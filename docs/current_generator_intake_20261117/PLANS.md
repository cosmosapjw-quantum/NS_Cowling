# Current generator manifest intake plan

## Frozen scope

- Base: `23b1f2099dc15470ac63a100e6b3aae0c397a36b`.
- Preserve the legacy 8 crust x 5 core x 2 method census (80 cells).
- Add a separate, explicit 13 crust x 5 core x 6 method manifest (390 requested cells).
- Permit at most one fresh matching execution: `SPG_M1__GM_GM1__css`, `n_bridge=40`.
- Reuse matching/readback evidence only when identity contracts match; reuse never implies a new attempt or admission.

## Milestones

1. **Contract and baseline — COMPLETE**
   - Remote identities and clean audit worktree verified.
   - Astra review approved manifest-intake-only scope.
2. **Counterexample tests — COMPLETE**
   - Cover tampering, missing/duplicate cells, aggregate count, SUPPORT/diagnostic loss,
     reuse non-promotion, and legacy compatibility.
3. **Bounded implementation — COMPLETE**
   - Manifest builder and separate pipeline intake; no production-physics changes.
4. **Layered validation — COMPLETE**
   - Targeted tests, related regression, one real end-to-end execution, readback and admission.
5. **Independent review and closeout — COMPLETE**
   - Read-only Astra decision: `PASS_SCOPED` for one CSS cell only.
   - Final census: 390 requested, 1 eligible, 1 attempted, 1 admitted.
   - The single execution allowance is consumed; no retry or expansion is authorized.

## Hard stops

- Zero generator or TF optimization runs; no 390-cell physical campaign.
- No energy offsets, tolerance relaxation, phase-label merging/deletion, SUPPORT bypass,
  source substitution, or consumer/stellar-physics edits.
- APR boundary and stellar-consumer work remains
  `HOLD_UNDECLARED_SOURCE_BOUNDARY_LAW` until endpoint, piecewise interpolant, and
  consumer-boundary receipts exist.
- This node cannot establish full EOS390 reproduction, full-crust support, original-author
  calibration, or full APR stellar admission.
