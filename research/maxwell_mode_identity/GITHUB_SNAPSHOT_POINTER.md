# GitHub snapshot pointer and local archive scope

Repository: https://github.com/cosmosapjw-quantum/NS_Cowling

LA11 branch: `research/maxwell-analytic-matching-la11-v1`
Code/source validation snapshot: `ffc6a75006d63c222871b9c9fc13abdbc699014a` (the helper and test blobs were verified at this commit). Fetch the current PR #10 head live; later metadata commits may advance the branch.
Draft PR #10: https://github.com/cosmosapjw-quantum/NS_Cowling/pull/10
PR base: `research/maxwell-interface-interpolation-la10-v1` at `104b4fa76e191db42e6af26f072533617f45441c`
PR state: OPEN / DRAFT / UNMERGED. CI status query returned no status records, which is NOT a CI PASS.

## Canonical project files on the branch

The sole machine-readable DAG is:
`research/maxwell_mode_identity/LOCAL_AGENT_DAG.v1.json` (contract `NSCOWLING-MODEID-WRONSKIAN-LOCAL-AGENT-v1.4`).

The single detailed development contract is:
`research/maxwell_mode_identity/LOCAL_AGENT_DAG_DEVELOPMENT_AND_HANDOFF_20261111.md`.

The canonical integrated handoff is:
`research/maxwell_mode_identity/HANDOFF_PROMPT_20261111.txt`.

LA11 scientific note, machine-readable results/status, helper, tests and computation receipts are in the same directory. Fetch these files from the branch at the live HEAD; do not rebuild or replace the SSOT from a patch.

## This archive

This archive is a **LA11 evidence and handoff companion**, not a full repository clone. It contains the LA11 audit, exact helper/test bytes, machine-readable LA11 result/status, the LA11 DAG patch for inspection, a standalone local-agent prompt, and computation evidence/receipts. The complete LA0-LA10 context and integrated single DAG/development document remain on the GitHub branch and are fetched live by the handoff prompt.

The local helper and test source bytes were matched to GitHub content blobs:
- helper blob SHA1 `97e9f5c6a134110e92afc68b1b2080dbabbf6756`; SHA256 `087ab499e9466b2008fedca79409e04db6038808bb8d4df03fde9d0f74f4b096`
- test blob SHA1 `8bedc9870cc80ae63b25f03d86454afd59941c9a`; SHA256 `32c31ac29d044c9141077ac5440b60689e6d0d6f11b2b67db291b881bd090e45`

Local receipt: 22 tests passed in 0.21 s; 22 collected in 0.01 s; JSON/AST/compileall and SHA256SUMS verification passed. The local run directory is not a Git worktree, so Git commit identity is taken from the GitHub branch, and the exact source Git blob identities are separately compared.

No native Fortran/q-to-r run, actual endpoint-frame analyticity proof, BVP Fredholm-pencil equivalence, physical multiplicity, stellar root/profile calculation or production integration is claimed.
