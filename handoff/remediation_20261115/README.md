# NS_Cowling EOS remediation handoff

**Purpose:** executable, evidence-bound remediation package for the audit findings recorded against PR #13.

- Repository: [cosmosapjw-quantum/NS_Cowling](https://github.com/cosmosapjw-quantum/NS_Cowling)
- Source branch at handoff start: `cloud/eos390-fortran-unification-v1`
- Frozen source HEAD: `772e85910afb40aea37c292680be24b363586a54`
- Remediation branch: `remediation/2026-11-15-m5a-maxwell-eos-contracts`
- Existing PR #13 remains open, draft and unmerged. This package does not merge or rewrite it.

## Package contents

| File | Purpose |
|---|---|
| `execution_manifest.json` | Machine-readable task graph, current patch state, invariants, blockers and acceptance commands |
| `IMPLEMENTATION_PLAN.md` | Ordered task plan with concrete files, interfaces, tests and stop conditions |
| `LOW_PERFORMANCE_AGENT_PROMPT.md` | Copy-ready prompt for a constrained coding agent |
| `SCIENTIFIC_CONTRACTS.md` | Derivations, units, CompOSE composition semantics and physics gates |
| `PATCH_REGISTER.json` | Patch list, file paths, commit IDs and verification state |

## Already patched in this branch

1. M5A native read-back now checks numerical values of (n_B,\mu_B,P,\epsilon), not only the row count.
2. M5A native evidence is scoped to the single headline barotrope. Per-case claims across 27 distinct parameter combinations are explicitly false until each artifact is separately checked.
3. A duplicate M5A JSON with a wrong model name was removed. The removed file was byte-identical (Git blob SHA `664ebbcf22254f3482aef645ecbaeece9c09fa24`) to the correctly named BBB_BHF2 artifact.
4. Maxwell v2 now explicitly rejects composition-bearing inputs because the Fortran v1 kernel is thermo-only. It no longer claims that composition is preserved in successful thermo-only output.
5. Maxwell v2 resolves the repository path from its module location rather than a developer-specific absolute path.
6. APR+CSS audit now evaluates the exact derivative of the PCHIP pressure interpolant and makes the differential first-law check part of the admission verdict.

## Verification state

**Not yet CI-verified.** The latest pre-remediation CI on the frozen source HEAD failed (259 passed, 9 skipped, 1 failed). This branch has not yet produced a GitHub Actions run at the time the package was assembled. The local execution container could not resolve `github.com`, so local clone/pytest execution was blocked. Do not infer PASS from the commits being present.

The test additions were committed before the corresponding implementation changes, but the host could not observe a test-only RED run. The red-observation field is therefore `NOT_OBSERVED`, not PASS.

## Continue here

1. Read `execution_manifest.json`.
2. Follow `IMPLEMENTATION_PLAN.md` in order.
3. Use the prompt in `LOW_PERFORMANCE_AGENT_PROMPT.md` without removing its invariants or stop conditions.
4. Run the GitHub Actions CI on this branch through the draft PR, then fix only evidenced regressions.
5. Keep unresolved physical/source blockers open. Never weaken a gate to turn the status green.

No raw CompOSE source EOS data are to be edited by this package.
