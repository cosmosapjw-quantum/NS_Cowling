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
- Follow-up exact-PR-head Actions run [42579770412](https://github.com/cosmosapjw-quantum/NS_Cowling/actions/runs/42579770412) checked commit `2b7a00984dfa215bcca48c7d7adb5f3619013502` and completed successfully. Its artifact ID is `13656363884`, SHA-256 `0cfba640210b84d3141213b820f5c335d67ea1abef4a9e39b204d4a41b82e4f6`.

The integrated scope is the phase-aware **Python/SciPy reference TOV/Cowling path** and independent Fortran slow/rapid transfer routines. Native Fortran stellar TOV/Cowling production solvers are still smooth-EOS-only. The 33-path v5 source overlay remains unapplied to this branch and is not represented as synced here. No merge has been performed.


## Native Fortran phase-segmented follow-up (2026-11-09)

The native Fortran implementation has now been ported to the full repository tree on sibling branch `research/maxwell-native-fortran-v1`, tracked by draft [PR #3](https://github.com/cosmosapjw-quantum/NS_Cowling/pull/3). Exact verified head: `98d2f9ce6033015109286e5713639c6626a6f7fb`.

- PR Actions run `42581660400`: success; artifact `13656748986`, SHA-256 `083f3eabf9215ac69cde06677f06c0fc58746bb35b24e5559c3a6184b5e7c5f6`.
- Push Actions run `42581655510`: success; artifact `13656434783`, SHA-256 `b6274813105f8e41f9eec269ffbe747ee3c5b594aea3a7bf36c3427a0cfb665e`.
- The exact-head log reports 12 cell checks and 18 Maxwell checks passed; the analytic small-jump TOV limit has `dM/M=9.1343e-8` and `dR/R=2.8981e-8`.
- `docs/ADVERSARIAL_AUDIT_MAXWELL_NATIVE_FORTRAN_20261109.md` records the verified scope and explicit non-claims; `docs/MAXWELL_NATIVE_FORTRAN_DAG_20261109.md` records the next implementation dependencies.

This follow-up does not merge PR #3. The normalized two-table production path has compile/build evidence but still needs a real end-to-end data fixture. The 33-path v5 overlay remains unapplied and must not be inferred to have been synchronized by this Maxwell integration.
## Final native table-provider and CLI verification (2026-11-09)

Exact verified branch head: `d3649906292b0ec8f1caf98961f12436fa00871b`.

- The generated analytic two-phase normalized-table fixture passed 74 checks, including table loading, provider admission, one-sided `P_t` coexistence, cold first-law closure, table-backed native TOV, and both slow/rapid Cowling shots.
- The workflow then executed the actual `case_driver` with `maxwell_table_cowling`; `status=ok`, `result_valid=true`, and `cowling_shot_valid=true`.
- The 1800 Hz surface residual is a trial residual, not a root. No mode identity or complete spectrum is claimed.
- Push Actions run [42582837414](https://github.com/cosmosapjw-quantum/NS_Cowling/actions/runs/42582837414): success, artifact 13657672304, SHA-256 `ec8333d0afc1fb65b57277c669333d403fbbbd6e18a19fb701aa369397bc2515`.
- PR Actions run [42582845328](https://github.com/cosmosapjw-quantum/NS_Cowling/actions/runs/42582845328): success, artifact 13658265995, SHA-256 `5d5a1683e4ba708f246c275e225b462bfc6968f201cf02d4956c9ceecf2cc278`.

The fixture is an analytic synthetic normalized table pair; a realistic source-backed CompOSE hybrid EOS has not yet been executed. PR #3 remains draft, open and unmerged. The 33-path v5 overlay also remains unapplied.

## Native production generic-dispatch increment (2026-11-09)

The earlier handoff statement that native Fortran `solve_tov_star` and `solve_cowling_spectrum` are smooth-only is superseded by the generic-dispatch implementation on `research/maxwell-native-fortran-v1`.

- Code head: `f8678cb13a18f63bad9f8db65c4cd3e66e491ad3`.
- Production façade: `oscillation_source/solutions/src/production/ns_production_solvers.f90`.
- Maxwell spectrum scanner: `oscillation_source/solutions/src/modes/ns_maxwell_mode_spectrum.f90`.
- Exact-head push workflow: [42584910411](https://github.com/cosmosapjw-quantum/NS_Cowling/actions/runs/42584910411), success; artifact 13658767086, SHA-256 `7a4ecb2784e777256735a809c5a331f9cb4f8bef85963d2fe1348045b17b4443`.
- Checks: 12 cell, 22 analytic Maxwell, 74 normalized-table.
- Smooth TOV cases and normalized-table TOV cases route through `solve_tov_star`; smooth spectrum/profile and Maxwell scans route through `solve_cowling_spectrum`. The Maxwell CLI includes `maxwell_table_spectrum`.
- A pre-dispatch snapshot remains at `backup/maxwell-native-fortran-pre-production-dispatch-20261109`.

The present CI exercises the Maxwell no-root scan path only. Positive-root bisection needs a deterministic bracket fixture; mode-family identity remains unimplemented. Keep PR #3 draft/open/unmerged. No local workstation checkout was changed by this ChatGPT-side GitHub operation.


## Exact-head production dispatch and spectrum CLI update (2026-11-09)

The statement that native Fortran production solvers are smooth-only is superseded by the generic production API and tests on the current feature branch.

- Exact code+workflow head exercised by CI: `22b0faa3c69150f3ae4fffcd2a6cae5d8231bb6c`.
- Code-bearing implementation head: `f8678cb13a18f63bad9f8db65c4cd3e66e491ad3`.
- Push [42585358415](https://github.com/cosmosapjw-quantum/NS_Cowling/actions/runs/42585358415): success; artifact 13658722733; SHA-256 `26313f4315af74718f5fa9d1687658ef603c405807b2f7d44d27057fa38242a1`.
- Pull request [42585366571](https://github.com/cosmosapjw-quantum/NS_Cowling/actions/runs/42585366571): success; artifact 13659056245; SHA-256 `12afce4e910c53b9a0dda180826bea2a98a5986562df0866bb7f9a26da219412`.
- Exact-head results: 12 cell, 22 analytic Maxwell and 74 normalized-table checks passed; production driver built; both `maxwell_table_cowling` and `maxwell_table_spectrum` passed.
- The spectrum smoke scanned three frequencies (1799, 1800, 1801 Hz), all successful, 0 brackets and 0 roots. It does not verify a physical eigenfrequency or positive-root bisection.
- The pre-dispatch branch backup remains `backup/maxwell-native-fortran-pre-production-dispatch-20261109`.

The production generic API now dispatches smooth or Maxwell TOV backgrounds and smooth or Maxwell Cowling scans. PR #3 remains draft/open/unmerged. A realistic source-backed CompOSE hybrid EOS remains untested. The positive-root branch and mode-family identity need dedicated validation. The 33-path v5 overlay remains unapplied.


## Native generic production dispatch and positive-root validation update (2026-11-09)

The former smooth-EOS-only production limitation is superseded. Native Fortran now exposes generic smooth/Maxwell overloads for `solve_tov_star` and `solve_cowling_spectrum`, and existing production case adapters route through them.

Latest exact code/test head: `9e92655a8277922f312f4c31df8e331e8c12df1b`.
Latest exact workflow head tested by push and PR: `ec4527f0d518f6ba19f57eb42d8a8eb00b98dfca`.

The generated analytic two-phase fixture executes both slow and rapid production spectrum scans on 41 frequencies from 100 to 20,000 Hz. Each conversion prescription returns two bisection-refined roots with a fresh final Cowling shot:

- Slow: 1764.8989380477 Hz (residual (-6.9575\times10^{-7})); 4937.6449322607 Hz (residual (9.3334\times10^{-8})).
- Rapid: 4937.0162338857 Hz (residual (8.5331\times10^{-9})); 19493.6326018255 Hz (residual (3.8951\times10^{-7})).

All roots are valid; family identity remains unresolved. These are synthetic fixture roots and not identified physical mode frequencies.

Both exact-head events pass:
- Push [42586438396](https://github.com/cosmosapjw-quantum/NS_Cowling/actions/runs/42586438396): artifact 13658709148, SHA-256 `02d1673cd82045aaacbae8125c74c7ce7b8896b7abe6def70e81df07331c3d04`.
- PR [42586446195](https://github.com/cosmosapjw-quantum/NS_Cowling/actions/runs/42586446195): artifact 13659137810, SHA-256 `478fb960cf8a993bb92bf2bbc544f17160b854e0f0f5ac5a28cc0b14a9ccaeaf`.

The strict CI requires two valid roots for both slow and rapid conversion. Root-frequency convergence, independent Python/SciPy cross-validation, phase-aware mode identity, realistic source-backed CompOSE data and stable M-R analysis remain open. PR #3 remains draft/open/unmerged.
