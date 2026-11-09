# GitHub handoff: v5 sync and sharp Maxwell reference v1

Date: 2026-11-09

## Repository and branch provenance

- Repository: cosmosapjw-quantum/NS_Cowling
- Existing EOS implementation PR: https://github.com/cosmosapjw-quantum/NS_Cowling/pull/1
- Existing PR branch: research/fortran-eos-gap-fill-v1
- Branch head at handoff: ca71d8e394f9043d778676c941f88ce01ca7faa9
- This artifact branch: research/maxwell-reference-validation-v1
- Artifact branch base: ca71d8e394f9043d778676c941f88ce01ca7faa9

The v5 GitHub overlay describes 33 paths relative to the exact existing PR branch head above. Its expected base tree is b0ee2168c06761c5c5d87fcc4e3ef6892b1b2309. The overlay's old blob identities were prepared against that tree. It has **not** been applied to the repository tree by this artifact-only branch.

## Patch and compact solver baseline

The sharp-interface patch is preserved at research_artifacts/maxwell/NS_Cowling_V5_MAXWELL_V1.patch.

It applies to this exact compact v5 solver ZIP, not directly to the full repository branch:

- Baseline solver ZIP SHA-256: 6d2ba545101df4cdc5bb5aed187b640ec84e528bf2005db764792bc4e19d9885.
- Maxwell-enabled compact source ZIP SHA-256: e7095d95b58022f95c1ee8bbc0692fb746eab1c676737388e58dae8ef82366b1.
- Patch SHA-256: de9b51aacbf8d022ffe5dbd0787b7f38d82f6883597f4b15963381530a0f05c5.

The patch was verified with git apply --check against a fresh extraction of the exact compact v5 baseline. Do not assume it applies to the current GitHub branch without first applying the 33-path v5 overlay and checking the actual result. The overlay archive itself is a separately prepared artifact; this branch does not contain the extracted overlay source tree.

## Contents of this branch

- research_artifacts/maxwell/NS_Cowling_V5_MAXWELL_V1.patch: complete 113,357-character unified source patch against the compact v5 baseline.
- docs/MAXWELL_REFERENCE_VALIDATION_V1_20261109.md: equations, phase-interface conditions, independent test values and explicit non-claims.
- docs/V5_PAIRED_MATCHING_SENSITIVITY_20261109.md: current v5 paired method sensitivity and interpretation limits.
- docs/GITHUB_MAXWELL_HANDOFF_20261109.md: exact provenance and application preconditions.

## Status

This PR is a **reproducibility and source-patch handoff**, not a claim that the v5 overlay or the sharp-interface source patch has been applied to the repository's production paths. The branch does not replace the native Fortran smooth-EOS stellar solver with a phase-segmented production backend.

The independent Python/SciPy sharp-interface reference path and Fortran transfer routines have the bounded tests described in the validation note. No GitHub Actions run is claimed for this patch. The broad compact regression rerun timed out; its earlier 76-pass/8-deselected result is recorded separately and is not promoted to a fresh full-suite PASS.

## Next integration gates

1. Apply the v5 overlay only after confirming the exact branch base SHA and all 33 old-blob identities.
2. Apply this sharp-interface patch only to the exact compact v5 snapshot or a separately reviewed port to the current full repository tree.
3. Add the sharp-interface module to the full repository build graph and extend the native Fortran TOV/Cowling production path if that is the desired end state.
4. Run focused tests, all feasible regression tests, manifest checks and exact-head GitHub Actions.
5. Keep true Maxwell discontinuity outputs separate from the old smooth maxwell matching cohort and from all previously released 234-mode results.


## Final GitHub execution update (2026-11-09)

The earlier status in this handoff records the pre-integration point. The sharp-interface source patch has since been applied to this branch and pushed after exact-head validation:

- Final source commit: `6b23d024e41661d33b83d6202038db6da918821d`.
- Follow-up documentation commits record the exact CI evidence; the current PR head is checked in PR metadata.
- GitHub Actions run [42579257980](https://github.com/cosmosapjw-quantum/NS_Cowling/actions/runs/42579257980) completed successfully.
- Native Fortran EOS build passed; focused Maxwell tests **7/7**; independent validation returned `PASS_BOUNDED_REFERENCE_VALIDATION`; full frontend/pipeline suite **74 passed**.
- Source commit was created by the workflow only after all validation stages passed.
- CI artifact ID `13657075483`, SHA-256 `93347a170c238458fa937f2359e9ae4a8737098c5f9ec9dbdb09a26351f9946b`.

The integrated scope is the phase-aware **Python/SciPy reference TOV/Cowling path** and independent Fortran slow/rapid transfer routines. Native Fortran stellar TOV/Cowling production solvers are still smooth-EOS-only. The 33-path v5 source overlay remains unapplied to this branch and is not represented as synced here. No merge has been performed.
