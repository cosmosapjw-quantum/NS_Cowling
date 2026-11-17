# Failure log

Preserve red tests and runtime failures here or in the linked receipt directory. Do not erase
earlier matching failures or repeat generator/TF optimization merely to populate the matrix.

- Initial collection failed because the manifest module did not exist; the red receipt is kept.
- Independent review blocked false requested/attempted semantics, missing cumulative allowance,
  confused launcher/binary identity, and exception-only missing-input handling. Counterexamples
  remain in tests and receipts.
- No physical execution occurred before repair and independent re-review.
- The authorized CSS execution succeeded on its first and only call. No retry was made.
- Remote runs 44014058466 and 44014058906 failed after the first push. The failures were
  classified as a stale tracked digest and a missing CI-only plotting dependency, respectively;
  raw identifiers and the bounded repair are preserved in `receipts/remote_ci_red.json`.
