# Retained-phase repair decision — 2026-11-17

Scope: Fortran retained-source phase transport, CompOSE phase-only readback,
receipt binding, and fail-closed CSS stellar intake. Base is
`4a44333377fc6252ad6bcc0c1024614277199f31`.

The fresh read-only reviewer `/root/phase_repair_review` returned
**PASS_SCOPED_REVIEW**, with no remaining blocking findings in the repaired diff.
Requested role: COSTROUTE research-review (`gpt-6-astra/xhigh`). Observed runtime
identity: `NOT_INDEPENDENTLY_OBSERVED`; the request/profile is not runtime proof.

Resolved findings, retained rather than removed from the record:

- **P1:** inferring safe legacy behavior from missing phase fields admitted the
  actual phase-erased APR baseline and constructed a stellar provider. Mandatory
  phase receipts now reject that same baseline in admission and intake. Historical
  files remain readable; their old receipt does not authorize new admission.
- **P2:** integer coercion concealed changed labels: `[0.5, 7.5]` became `[0, 7]`
  while the API claimed preservation. Exact one-dimensional native-range integer
  validation now rejects fractional, nonfinite, boolean, nonvector and out-of-range
  labels before native execution. Original labels are compared without truncation.

The reviewer independently ran 123 targeted tests and 12 support-contract tests,
checked retained native slices, unique bridge labels, phase-only masks and endpoint
stencil confinement. All 27 source hashes, binary hash and five output inventories
in `phase_retention_repair/assembly_003/RECEIPT.json` matched. Intact APR intake
retains boundary pairs (221,222) and (232,233) and raises the typed unsupported
boundary error. The public APR source SUPPORT guard remains unchanged.

Review bindings, also verified by the implementation controller:

```text
production diff SHA256: 48c90c8e5352819d65aaa1b0142661024049ee0a801710f0490294ec559b25da
assembly_003/RECEIPT.json SHA256: 5f61b7babd6d65f358402ccaadcbd941105d16f0dfd03f46be0661d019b8dc97
```

The production diff hash is `git diff --` over the matching module, native app,
EOS API, CSS adapter and match-report admission script (in that path order).
Final full frontend/pipeline regression: **1162 passed, 9 skipped**. The final
data inventory and validation commands are preserved with the report.

Decision: publish the scoped implementation repair under the owner's push
authorization, then verify the containing HEAD's remote CI. This does not admit
full APR matching/stellar evolution, prescribe source boundary interpolation or
interface laws, or establish original-author EOS equivalence. Those remain HOLD
or not established. The next executable step is the authorized commit/push/CI
verification; subsequent physical support work starts from the existing retained
data and boundary-specific evidence, without replaying the full campaign.
