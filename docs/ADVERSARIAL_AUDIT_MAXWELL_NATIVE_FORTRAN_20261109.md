# Adversarial audit: native sharp-Maxwell Fortran TOV/Cowling

Date: 2026-11-09

## Authority and exact source identity

- Repository: `cosmosapjw-quantum/NS_Cowling`
- Branch: `research/maxwell-native-fortran-v1`
- Pull request: [#3](https://github.com/cosmosapjw-quantum/NS_Cowling/pull/3), draft, open, not merged
- Exact verified head: `98d2f9ce6033015109286e5713639c6626a6f7fb`
- TOV module blob: `6e54f73c07ae68c8113ceb3ff31aff2dd12170e6`
- Cowling module blob: `ce4fac829d3a5c211139d5174541a7053793afc6`
- Table-driven case adapter blob: `1b0853565884a43696afc15cfb072e2bde033231`
- Analytic test fixture blob: `7dc7f6f525b34314257f0645fa1e8026ff85237c`

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
| Realistic normalized EOS table execution | NOT RUN | CLI/provider path compiles, but no real two-table fixture was executed |
| Complete Cowling spectrum / mode identity | NOT IMPLEMENTED | Single-frequency shooting residual only |

### Exact-head GitHub Actions

Both events checked the same commit `98d2f9ce6033015109286e5713639c6626a6f7fb`:

- Pull-request run [42581660400](https://github.com/cosmosapjw-quantum/NS_Cowling/actions/runs/42581660400): success. Artifact ID `13656748986`, SHA-256 `083f3eabf9215ac69cde06677f06c0fc58746bb35b24e5559c3a6184b5e7c5f6`.
- Push run [42581655510](https://github.com/cosmosapjw-quantum/NS_Cowling/actions/runs/42581655510): success. Artifact ID `13656434783`, SHA-256 `b6274813105f8e41f9eec269ffbe747ee3c5b594aea3a7bf36c3427a0cfb665e`.
- Exact-head job log: `cell_checks: 12 checks passed`; `maxwell_checks: 18 checks passed`; `dM/M=9.1343e-8`; `dR/R=2.8981e-8`.
- Compile/test flags: `-std=f2018 -O0 -g -fcheck=all -ffree-line-length-none -fimplicit-none`.

## Scientific and numerical review

1. The high-density callback is used for the core and the low-density callback for the envelope. The TOV integrator terminates the core at exact `P_t` and restarts the envelope using the same `(r_t,m_t,nu_t)` state. Separate background tables prevent interpolation through the energy-density jump.
2. At `P_t`, the implementation checks pressure equality, chemical-potential continuity, positive `Delta epsilon` and `Delta n_B`, and the cold first-law identity `Delta epsilon = mu_t Delta n_B`. It does not solve for phase coexistence from arbitrary EOS tables; the two supplied tables must already encode coexistence.
3. The Seidov ratio is emitted as an onset diagnostic. It is not a radial-stability calculation or a substitute for tracing stable branches.
4. The Cowling path integrates the core and envelope separately and then applies the existing native slow/material or rapid conversion map. Its test checks that the mapped pressure perturbation agrees on both sides to relative tolerance `1e-9`.
5. The small-jump test compares the new phase-segmented TOV result with the existing smooth single-EOS TOV solution. It is a cross-implementation limit check for this analytic fixture, not a realistic hybrid-star prediction.

## Defects caught during the integration loop

- First CI attempt: a stale `pi` symbol in the analytic test fixture stopped compilation. It was replaced by the imported `two_pi` constant.
- Next runtime attempt: the new pressure-coordinate state was mistakenly populated in the legacy `(m,P,nu)` order rather than `(r,m,nu)`. The central pressure was separated into its own scalar output.
- The subsequent run exposed an insufficient central-series initial-step bound. The initial logarithmic-pressure step is now bounded by `0.1*(P_c-P_start)/P_c` as well as the phase pressure span.
- The final exact head passed both PR and push workflows. Earlier failed runs are retained in Actions history and are not relabelled as passes.

## Residual risks and explicit non-claims

- The normalized high-/low-density table provider path and NAMELIST CLI have compile evidence only; an actual two-table run remains open.
- No published realistic hybrid-EOS TOV/Cowling benchmark has been reproduced in this PR.
- No complete frequency root scan, `f/p1` continuation, mode identity tracking, interface-mode classification, or convergence under a frequency-root refinement study is implemented.
- No maximum-mass or stable-branch survey, twin-star analysis, radial stability calculation beyond the local Seidov diagnostic, or parameter inference is claimed.
- No metric perturbations beyond Cowling, rotation, magnetic fields, elasticity, superfluidity, composition stratification, or GW damping.
- The finite positive-pressure surface cutoff remains part of the numerical contract.
- This PR is separate from the 234-row v5 smooth-matching spectrum and does not apply the separately prepared 33-path v5 overlay.

## Audit conclusion

`PASS_BOUNDED_NATIVE_ANALYTIC_FIXTURE` for the exact source head above. The core native phase-segmented TOV/Cowling path compiles, the analytic discontinuity fixture passes, and the production case driver links. The implementation is suitable for the next integration gate, but it is not yet evidence for realistic hybrid-star predictions or observational inference.