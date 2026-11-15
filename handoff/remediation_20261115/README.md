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
| `../../docs/CSS_SINGLE_SPEED_CRUST_JUMP_CONTRACT_20261115.md` | Single-speed CSS matching equations, jump contract, and verification scope |
| `PATCH_REGISTER.json` | Patch list, file paths, commit IDs and verification state |

## Already patched in this branch

1. M5A native read-back now checks numerical values of (n_B,\mu_B,P,\epsilon), not only the row count.
2. M5A native evidence is scoped to the single headline barotrope. Per-case claims across 27 distinct parameter combinations are explicitly false until each artifact is separately checked.
3. A duplicate M5A JSON with a wrong model name was removed. The removed file was byte-identical (Git blob SHA `664ebbcf22254f3482aef645ecbaeece9c09fa24`) to the correctly named BBB_BHF2 artifact.
4. Maxwell v2 now explicitly rejects composition-bearing inputs because the Fortran v1 kernel is thermo-only. It no longer claims that composition is preserved in successful thermo-only output.
5. Maxwell v2 resolves the repository path from its module location rather than a developer-specific absolute path.
6. APR+CSS audit now evaluates the exact derivative of the PCHIP pressure interpolant and makes the differential first-law check part of the admission verdict.
7. The native M7 `css` matcher is changed to one constant-sound-speed branch with a first-order crust-to-bridge density jump and a continuous bridge-to-core state match. The full contract is in `docs/CSS_SINGLE_SPEED_CRUST_JUMP_CONTRACT_20261115.md`.

## Verification state

**Previous verified CI run (pre single-speed CSS change):** [run 43676024772](https://github.com/cosmosapjw-quantum/NS_Cowling/actions/runs/43676024772), on code-patch HEAD `5d31b41eacf267cfae64e9d60a86a1715a774f14` (see authoritative SHA in `execution_manifest.json`). Fortran build passed; full frontend/direct-provider test command returned **263 passed, 9 skipped, 0 failed**; all workflow steps, including native Fortran unit tests and SHA-256 manifest check, completed successfully.

The two Wigner-Seitz source-data tests and three Cowling E2E tests were skipped because their raw inputs were not staged. Also, some generated core/crust-only component EOS receipts explicitly report C2/C10 physical-validation failures, and GPPVA_TW remains blocked due a provider-source discrepancy. Thus workflow success is not a claim that every generated EOS is physically admissible.

The local execution container could not resolve `github.com`, so local clone/pytest execution was blocked. Focused test-first additions were committed before implementation, but the host could not observe a test-only RED run; the red-observation status remains `NOT_OBSERVED`.

The single-speed CSS implementation initially reached CI run [43677499924](https://github.com/cosmosapjw-quantum/NS_Cowling/actions/runs/43677499924). The Fortran build passed; pytest returned **248 passed, 9 skipped, 15 failed**. Root causes were an endpoint equality check incorrectly applied to the interior-only production grid, an omitted `composition_dropped` report assignment, and stale SHA256 manifest entries. The endpoint check and metadata assignment have since been corrected, and a sharp-seam derivative diagnostic was made not-applicable rather than comparing derivatives across a first-order jump.

A follow-up workflow [43678164119](https://github.com/cosmosapjw-quantum/NS_Cowling/actions/runs/43678164119) was queued for code head `654bdc513b6f3c77f80d828898f81f213b12d02e`, but was **cancelled by the workflow concurrency group** after subsequent handoff/evidence commits changed the PR head. The current final branch snapshot therefore still requires one fresh exact-head workflow run. No CI PASS is claimed for the single-speed CSS implementation yet.

## Continue here

1. Read `execution_manifest.json`.
2. Follow `IMPLEMENTATION_PLAN.md` in order.
3. Use the prompt in `LOW_PERFORMANCE_AGENT_PROMPT.md` without removing its invariants or stop conditions.
4. Run the GitHub Actions CI on this branch through the draft PR, then fix only evidenced regressions.
5. Keep unresolved physical/source blockers open. Never weaken a gate to turn the status green.

The old 2026-11-14 two-segment smooth CSS artifact is historical evidence for the previous implementation and must not be overwritten or treated as validation of the new contract. No raw CompOSE source EOS data are to be edited by this package.
