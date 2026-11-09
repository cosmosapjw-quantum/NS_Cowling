# Adversarial audit: native sharp-Maxwell Fortran TOV/Cowling

Date: 2026-11-09

## Authority and exact source identity

- Repository: `cosmosapjw-quantum/NS_Cowling`
- Branch: `research/maxwell-native-fortran-v1`
- Pull request: [#3](https://github.com/cosmosapjw-quantum/NS_Cowling/pull/3), draft, open, not merged
- Historical verified source head: `d3649906292b0ec8f1caf98961f12436fa00871b`
- Latest code-bearing implementation head: `f8678cb13a18f63bad9f8db65c4cd3e66e491ad3`
- Exact code+workflow head exercised by CI: `22b0faa3c69150f3ae4fffcd2a6cae5d8231bb6c`
- TOV module blob: `6e54f73c07ae68c8113ceb3ff31aff2dd12170e6`
- Cowling module blob: `ce4fac829d3a5c211139d5174541a7053793afc6`
- Table-driven case adapter blob: `1b0853565884a43696afc15cfb072e2bde033231`
- Analytic test fixture blob: `7dc7f6f525b34314257f0645fa1e8026ff85237c`
- Normalized-table fixture blob: `e09a6e5cdb0dcd1ca55fcbbbb3ecee6a619f0a51`

## Gate results

| Gate | Result | Evidence / boundary |
|---|---|---|
| gfortran compile of all solutions modules | PASS | GitHub Actions; GNU Fortran 15.2.0 |
| FPM test graph | PASS | `fpm 0.13.0`; runtime checks enabled (`-fcheck=all`) |
| Existing native cell checks | PASS | 12 checks |
| Analytic sharp-Maxwell fixture | PASS | 18 checks |
| Production `case_driver` link | PASS | Exact-head `fpm build` stage |
| TOV small-jump smooth limit | PASS | `dM/M = 9.1343e-8`; `dR/R = 2.8981e-8` |
| Slow/rapid interface shot | PASS, bounded | Both shots valid at one 1800 Hz trial; both `Delta P` conservation checks below `1e-9` |
| Synthetic normalized-table TOV/Cowling integration | PASS | Two generated T14-format phase tables loaded and admitted; exact-Pt coexistence and first-law closure checked; table-backed TOV and slow/rapid Cowling shots passed |
| Production case-driver table CLI | PASS, synthetic fixture | `maxwell_table_cowling` returned `result_valid=true` and `cowling_shot_valid=true`; one 1800 Hz residual only |
| Realistic CompOSE hybrid-EOS execution | NOT RUN | No real source-backed high/low phase pair has been supplied and audited |
| Cowling frequency scan / mode identity | PARTIAL / IDENTITY NOT IMPLEMENTED | Finite supplied-grid slow/rapid scans and bracket refinement are implemented; the current fixture only verifies the no-root scan path. No positive-root branch, f/p1 continuation, or interface-mode identity is verified. |

### Exact-head GitHub Actions

Both events checked the same final source head `d3649906292b0ec8f1caf98961f12436fa00871b`:

- Push run [42582837414](https://github.com/cosmosapjw-quantum/NS_Cowling/actions/runs/42582837414): success. Artifact ID `13657672304`, SHA-256 `ec8333d0afc1fb65b57277c669333d403fbbbd6e18a19fb701aa369397bc2515`.
- Pull-request run [42582845328](https://github.com/cosmosapjw-quantum/NS_Cowling/actions/runs/42582845328): success. Artifact ID `13658265995`, SHA-256 `5d5a1683e4ba708f246c275e225b462bfc6968f201cf02d4956c9ceecf2cc278`.
- FPM 0.13.0 and GNU Fortran 15.2.0; compile/test flags `-std=f2018 -O0 -g -fcheck=all -ffree-line-length-none -fimplicit-none`.
- Native checks: `cell_checks: 12`; `maxwell_checks: 18`; `maxwell_table_checks: 74`.
- Production CLI smoke: `maxwell_table_cowling` emitted `status=ok`, `result_valid=true`, and `cowling_shot_valid=true` for the generated two-table fixture.
- Table-backed fixture values: (M=9.6356877\times10^{-4}M_\odot^{\rm nominal}), (R=642.793\,\mathrm{m}), (r_t=332.172\,\mathrm{m}), Seidov ratio `0.21276596` below critical `0.50159574`.
- At 1800 Hz the Cowling integration returned a surface residual `0.9935236`, which is a non-root trial value, not an eigenfrequency. The slow-conversion interface pressure perturbations agree to about `3.2e-16` relative; 854 profile points were emitted.

Both workflows also ran `fpm build` for the production case driver and uploaded exact-head source identity evidence.

## Scientific and numerical review

1. The high-density callback is used for the core and the low-density callback for the envelope. The TOV integrator terminates the core at exact `P_t` and restarts the envelope using the same `(r_t,m_t,nu_t)` state. Separate background tables prevent interpolation through the energy-density jump.
2. At `P_t`, the implementation checks pressure equality, chemical-potential continuity, positive `Delta epsilon` and `Delta n_B`, and the cold first-law identity `Delta epsilon = mu_t Delta n_B`. It does not solve for phase coexistence from arbitrary EOS tables; the two supplied tables must already encode coexistence.
3. The Seidov ratio is emitted as an onset diagnostic. It is not a radial-stability calculation or a substitute for tracing stable branches.
4. The Cowling path integrates the core and envelope separately and then applies the existing native slow/material or rapid conversion map. Its test checks that the mapped pressure perturbation agrees on both sides to relative tolerance `1e-9`.
5. The small-jump test compares the new phase-segmented TOV result with the existing smooth single-EOS TOV solution. It is a cross-implementation limit check for this analytic fixture, not a realistic hybrid-star prediction.
6. A second test constructs two normalized T14 tables with (g=\epsilon-P) linear in \(\ln P\), runs the actual reader and provider, and feeds those providers into native TOV/Cowling. The matching construction makes PCHIP and the first-law quadrature representation-exact; it is an integration fixture, not observational EOS data.
7. The production case-driver is separately run against those generated files. This verifies NAMELIST parsing, table admission, stellar integration, and a single-frequency Cowling shot in the actual CLI path.

## Defects caught during the integration loop

- First CI attempt: a stale `pi` symbol in the analytic test fixture stopped compilation. It was replaced by the imported `two_pi` constant.
- Next runtime attempt: the new pressure-coordinate state was mistakenly populated in the legacy `(m,P,nu)` order rather than `(r,m,nu)`. The central pressure was separated into its own scalar output.
- The subsequent run exposed an insufficient central-series initial-step bound. The initial logarithmic-pressure step is now bounded by `0.1*(P_c-P_start)/P_c` as well as the phase pressure span.
- The final exact head passed both PR and push workflows. Earlier failed runs are retained in Actions history and are not relabelled as passes.

## Residual risks and explicit non-claims

- The normalized high-/low-density table provider and production CLI are verified only on a generated analytic two-table fixture. A source-backed CompOSE hybrid-EOS pair remains untested.
- No published realistic hybrid-EOS TOV/Cowling benchmark has been reproduced in this PR.
- A bounded frequency-grid scanner and bisection path are implemented, but this exact-head fixture has no positive bracket. Positive-root refinement convergence, `f/p1` continuation, mode identity tracking and interface-mode classification remain unverified.
- No maximum-mass or stable-branch survey, twin-star analysis, radial stability calculation beyond the local Seidov diagnostic, or parameter inference is claimed.
- No metric perturbations beyond Cowling, rotation, magnetic fields, elasticity, superfluidity, composition stratification, or GW damping.
- The finite positive-pressure surface cutoff remains part of the numerical contract.
- This PR is separate from the 234-row v5 smooth-matching spectrum and does not apply the separately prepared 33-path v5 overlay.

## Audit conclusion

`PASS_BOUNDED_NATIVE_SYNTHETIC_TABLE_AND_CLI` for exact source head `d3649906292b0ec8f1caf98961f12436fa00871b`. The callback and normalized-table fixtures pass; the production CLI loads two separate normalized phase tables, integrates the sharp-interface TOV background, and completes both slow/rapid test shots. The 1800 Hz CLI value is a trial residual, not a mode root. This is not evidence for a realistic CompOSE hybrid-star prediction, a complete spectrum, stable-branch classification, or observational inference.

## Production generic dispatch follow-up (2026-11-09)

### Exact source head

- Dispatch implementation head: `f8678cb13a18f63bad9f8db65c4cd3e66e491ad3`
- New module `ns_production_solvers.f90` blob: `3539dac7f6aaefb5af76c493ff26e009df0e7010`
- New module `ns_maxwell_mode_spectrum.f90` blob: `726a03b219827e77412175c68e3d7c28ef32ef84`
- Production Maxwell adapter blob: `6b911462981cdf503efe649232aae002bbded790`
- Generic-dispatch API regression test blob: `456cd365ecfcdbcb1a0f1055195b935255da14c1`

### Verified result

GitHub Actions push run [42584910411](https://github.com/cosmosapjw-quantum/NS_Cowling/actions/runs/42584910411) passed on exact head `f8678cb13a18f63bad9f8db65c4cd3e66e491ad3`; artifact ID `13658767086`, SHA-256 `7a4ecb2784e777256735a809c5a331f9cb4f8bef85963d2fe1348045b17b4443`.

- GNU Fortran 16.2.0, FPM 0.13.0; runtime checks enabled.
- `cell_checks`: 12 checks passed.
- `maxwell_checks`: 22 checks passed, including smooth and Maxwell TOV generic dispatch, smooth-profile Cowling generic dispatch and slow/rapid Maxwell frequency-grid no-root scans.
- `maxwell_table_checks`: 74 checks passed.
- Production case-driver build and previous `maxwell_table_cowling` synthetic table smoke passed.
- The optional smooth `eos_events(:)` contract is preserved through the generic TOV wrapper. A first run exposed this interface omission; the initial failure is retained in Actions history and the follow-up exact-head run passed after the wrapper was corrected.

### Claim boundary

This run verifies compilation and the generic dispatch/no-root scanning path only. The frequency fixture around 1800 Hz contains no candidate root, so the bisection-positive-root branch was not exercised. Even when a bracketed root is found, mode family identity remains `unresolved`. No realistic CompOSE hybrid EOS, physical eigenfrequency, complete spectrum, f/p1 continuation, or stable-branch inference is claimed.


## Final generic production dispatch and CLI-spectrum gate (2026-11-09)

The latest branch head is `22b0faa3c69150f3ae4fffcd2a6cae5d8231bb6c`; its code-bearing implementation parent is `f8678cb13a18f63bad9f8db65c4cd3e66e491ad3`. Both exact-head workflow events passed:

- Push [42585358415](https://github.com/cosmosapjw-quantum/NS_Cowling/actions/runs/42585358415), artifact 13658722733, SHA-256 `26313f4315af74718f5fa9d1687658ef603c405807b2f7d44d27057fa38242a1`.
- Pull request [42585366571](https://github.com/cosmosapjw-quantum/NS_Cowling/actions/runs/42585366571), artifact 13659056245, SHA-256 `12afce4e910c53b9a0dda180826bea2a98a5986562df0866bb7f9a26da219412`.

Runtime-check results: `cell_checks=12`, `maxwell_checks=22`, `maxwell_table_checks=74`; production build passed. The generated analytic two-table fixture passed both `maxwell_table_cowling` and `maxwell_table_spectrum` through the actual `case_driver`.

The spectrum CLI scanned 1799, 1800 and 1801 Hz. All three shots succeeded with residuals approximately 0.989787, 0.993524 and 0.996295. It found zero sign brackets and zero roots. This verifies a successful finite-grid no-root scan only. The positive-root/bisection branch is not exercised. A successful scan is not a complete spectrum and none of these residual samples is an eigenfrequency.

The implementation also preserves the optional smooth `eos_events(:)` argument through the generic TOV wrapper. A compile failure from its initial omission is recorded at earlier commit `f0723ec...`; the corrected exact-head run passed. A later malformed workflow authoring commit `42b5936...` created no jobs; the workflow was restored from the previous clean YAML and both final events passed.

**Current audit state:** `PASS_BOUNDED_NATIVE_GENERIC_DISPATCH_AND_SYNTHETIC_NO_ROOT_CLI`. Still unverified: source-backed realistic CompOSE phases, positive-root refinement, root-frequency convergence, mode-family identity/continuation, M-R stability sequence and astrophysical inference.


## Positive-root scan and strict slow/rapid CI (2026-11-09)

The fixture writer emits separate slow and rapid NML files and scans 41 frequencies over 100-20,000 Hz. Each scan finds exactly two sign brackets; the bisection path then returns two valid roots with a fresh final Cowling shot.

| Conversion | Root | Frequency (Hz) | Final surface residual | Status |
|---|---:|---:|---:|---|
| slow | 1 | 1764.8989380477 | (-6.9575\times10^{-7}) | valid / ok |
| slow | 2 | 4937.6449322607 | (9.3334\times10^{-8}) | valid / ok |
| rapid | 1 | 4937.0162338857 | (8.5331\times10^{-9}) | valid / ok |
| rapid | 2 | 19493.6326018255 | (3.8951\times10^{-7}) | valid / ok |

All root identities remain `unresolved`. The slow and rapid sets are not paired by root order or nearby frequency; a phase-aware identity invariant is not yet available. The values belong to the synthetic normalized two-table fixture only.

The workflow now requires exactly two roots and both root-valid flags for each conversion prescription. Exact tested heads:
- Code/test: `9e92655a8277922f312f4c31df8e331e8c12df1b`.
- Workflow: `ec4527f0d518f6ba19f57eb42d8a8eb00b98dfca`.

Both workflows passed:
- Push [42586438396](https://github.com/cosmosapjw-quantum/NS_Cowling/actions/runs/42586438396), artifact 13658709148, SHA-256 `02d1673cd82045aaacbae8125c74c7ce7b8896b7abe6def70e81df07331c3d04`.
- PR [42586446195](https://github.com/cosmosapjw-quantum/NS_Cowling/actions/runs/42586446195), artifact 13659137810, SHA-256 `478fb960cf8a993bb92bf2bbc544f17160b854e0f0f5ac5a28cc0b14a9ccaeaf`.

Toolchain: GNU Fortran 16.2.0, FPM 0.13.0; runtime-check flags include `-fcheck=all`. Counts: 12 cell, 22 analytic Maxwell and 74 normalized-table checks.

**Updated audit state:** `PASS_BOUNDED_NATIVE_GENERIC_DISPATCH_AND_SYNTHETIC_SLOW_RAPID_POSITIVE_ROOTS`. Remaining: root-frequency convergence, independent Python/SciPy cross-validation, negative table admission, mode identity/continuation, stable M-R sequence and source-backed realistic EOS.
