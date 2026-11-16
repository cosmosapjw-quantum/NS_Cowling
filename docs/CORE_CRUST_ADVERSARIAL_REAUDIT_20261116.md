# NS_Cowling adversarial re-audit, 2026-11-16

Baseline: `eda4af96743c9af3369d1810ef733badb1709b93` on PR #15.
Implementation: `d21d44a73588365cbf83ed93d3d5c44604cb76cd` plus final Cowling residual gate `42950053088fa80752826f4aae56d0a0b013c0e4`.
Previous `b0a3f10` handoff and PR body were stale: 26 subsequent commits already implemented interval checks, v2–v5 reconstructions and stellar adapters. This audit reads those changes and tests their actual physical contracts.

## Findings and repairs

1. **Full CLD stationarity was not imposed.** The conditional composition solve used bulk `mu_n-mu_p=mu_e`, omitting the composition derivatives of Coulomb and surface energy. SPG_M1 examples had a bulk residual near 1e-10 MeV but full cell derivative 6.60–7.84 MeV; native DD2 examples gave 9.64–11.30 MeV. The corrected condition at fixed f,R,n_I is

   `0 = mu_p^I-mu_n^I+mu_e + 4*pi*e^2*Phi*R^2*n_Ip + 3*sigma'(y)/(R*n_I)`.

   Python's outer objective now also includes `(m_p-m_n)*n_e`; the omitted `m_n*n_B` is constant at fixed n_B. Native already minimized total energy but mislabeled its 942–946 MeV/baryon output as binding energy; it now reports total and baryon-rest-subtracted energy separately. Both implementations record actual beta and optimizer convergence. Native full-coordinate scaled gradients fell from as much as 7.91e-3 to at most 2.53e-6 MeV/baryon in the three independent finite-difference cases.

2. **The declared interface model was not the cited model.** Scurto et al., arXiv:2402.15548v1, Eqs. 24–25, use `Phi=(2-3*f^(1/3)+f)/5` and `sigma=sigma0*(b+16)/(b+y^-3+(1-y)^-3)`, with model-specific fits to nuclear masses. The repository's historical `sqrt(2-3*f)/2` and `sigma0*(y*(1-y))^b` are different. They are retained explicitly as an exploratory proxy functional, with `source_calibrated:false`. No fitted parameters were invented. The source's `b` and the proxy exponent `b` have different meanings.

3. **CLD emission checked rows, not its declared cubic between rows.** With `epsilon=m_n*n+R(n)`, `mu=m_n+R'`, and `c_s^2/c^2=n*R''/mu`, stability requires `R''>0` and causality requires `D=mu-n*R''>=0`. On each cubic interval `R''` is affine and `D'=-n*R'''` has constant sign. Endpoint checks therefore cover the whole declared interval in real arithmetic; implementation evaluates them in float64. Hidden instability now rejects SPG_M5. Missing gas chemical potentials, nonfinite/invalid physical fields, all optimizer box bounds, and stale output directories are also checked.

4. **TOV causal limits used incompatible maps.** `n(epsilon_max)` did not invert the actual central `P(n)`. A synthetic case returned `P_limit=530.1666066815916`, yet `P(n_max)=530.2551879590255`, with RHS `c_s^2/c^2=1.000156758`. Inverting the actual `P(n)` yields 530.1666066815911 and 0.9999999999999989. The central expansion now uses the RHS's `epsilon(P_c)`. The first unresolved interval terminates the admissible connected prefix. Phase-aware TOV now checks every traversed pressure branch, rather than the center alone, using its actual RHS interpolant and the unchanged 1e-5 gate.

5. **CSS was changed by the stellar adapter.** A combined CSS+core PCHIP changed the declared constant sound speed by up to 54.0694% on preserved v5 SPG_M1+BBB_BHF2 input. `CSSCoreBarotrope` now evaluates `epsilon=epsilon_t+(P-P_t)/s` and `n=n_t*((epsilon+P)/(epsilon_t+P_t))^(1/(1+s))` exactly on the CSS segment, with a separate core interpolant and one-sided derivative at the continuous join. The same test now has zero speed error and 3.33e-16 relative energy error. The historical input's slow root shifts from 2545.963441 to 2536.972166 Hz; this is an implementation comparison, not a v6 physical prediction.

6. **Boundary and evidence contracts could report false success.** Fixed fractional phase/index truncation, dimensionful pressure tolerance, admitted interface-pressure roundoff causing domain errors, final-grid endpoint root loss, nested diagnostics stringification, stale or modified census bytes, receipt overwrite, missing-executable census exit 0, and rejected-source markers ignored when old EOS files remain. Input fingerprints and source scope survive failed attempts. Expected physical rejections remain recorded cases; environment/implementation/unclassified execution failures now make the census exit nonzero. Scanned Cowling candidates must additionally satisfy the existing native default normalized residual tolerance 1e-5.

7. **Python matching could repair invalid inputs silently.** A source `Q3 += 0.01` corruption (~9.4 MeV) on retained or discarded rows still produced a 611-row match. Every source row now passes the same 5e-4 Euler gate as native matching before window selection. Il=0 includes the charge term. Physical energy normalization is unchanged.

## Self-adversarial follow-through

The first corrections were themselves attacked. Fractional proxy exponent b=0.59 exposed a composition maximum selected as equilibrium; finite uniform/logarithmic root-candidate searches now compare full conditional energies. Tiny positive CSS speed exposed `inf > inf` silently accepting invalid analytic states; finite-state checks reject it. A discontinuous sign jump exposed frequency-bracket convergence without a residual zero; final-shot residual admission now rejects it. Native energy-spread stopping alone left appreciable coordinate gradients; simplex and energy convergence were tightened. These are bounded numerical checks, not proofs of global minimization or spectral completeness.

## New data and execution scope

`download/reconstruction_v6_proxy/` contains SPG_M1 (16 rows), GPPVA_DD2 (19), DDME2 (14), TW (14). M2/M4 remain placeholders; M3 is rejected for internal optimizer-bound rows; M5 is rejected for interval-13 instability. Earlier raw and reconstruction trees remain unchanged. Generation receipts record the actual precommit dirty tree and exact generating source hashes; those hashes were checked against the committed source.

The current census has 80 configured, 40 evaluated, 40 source-marked attempts, 8 CSS survivors and 0 Maxwell. Twenty Maxwell attempts stop at the existing mass-convention gate; four CSS APR attempts fail native C10; eight CSS RG attempts fail the beta contract. Counts do not establish microscopic source validity. Every emitted crust is explicitly `exploratory_proxy_CLD`.

The matched stellar model has a finite-pressure cutoff without an outer crust. Imposing Delta P=0 there and normalizing the lapse as vacuum defines an artificial truncated model. Its radii and frequencies are not complete-star observables. Sampled mass maxima, Seidov diagnostics and real Cowling roots do not certify radial-mode stability. The frequency scan reports unlabelled roots and its limited search interval.

All eight v6 proxy survivors were executed from clean committed source `4295005`: each produced a 1.4-solar-mass truncated star and one slow plus one rapid root in 800–4000 Hz. All 16 final normalized residuals meet 1e-5. The machine-readable stellar receipt stores executed-source hashes, input hashes, runtime, node diagnostics and limited scan scope.

A final state-maintenance regression also prevents the historical DAG updater from overwriting the newer audit status. Its old behavior was reproduced before the guard was added.

## Validation and reproducibility

CI #612 (implementation HEAD `4295005`, tested merge `31e596acf975527536575d808da31f3b907abd44`) passed **459 frontend/pipeline tests with 9 skips**, **142 generator tests**, every native unit/build step, and both checksum manifests (187 root / 439 generator entries at that commit).

Machine receipt and retained logs: `data/receipts/core_crust_audit_20261116/`.
Run commands from the repository root:

```sh
bash fortran_eos/build.sh build/compose_eos_fortran
PYTHONPATH=python_frontend COMPOSE_EOS_FORTRAN_ROOT=. COMPOSE_EOS_FORTRAN_EXE=build/compose_eos_fortran/compose_eos_fortran python3 -m pytest -q python_frontend/tests pipeline/tests
PYTHONPATH=eos_generator:eos_generator/compose_eos_pkg python3 -m pytest -q eos_generator/tests
COMPOSE_EOS_FORTRAN_EXE=build/compose_eos_fortran/compose_eos_fortran NS_COWLING_MATCH_OUT_ROOT=build/fresh_census python3 scripts/match_core_crust.py
COMPOSE_EOS_FORTRAN_EXE=build/compose_eos_fortran/compose_eos_fortran python3 scripts/run_css_survivor_stars.py build/fresh_census build/fresh_stars.json
```

Use absent/empty output paths. Native polynomial, CSS, DD-RMF, CLD, Wigner–Seitz, TOV–Love and Cowling build checks are the executable commands in `.github/workflows/fortran-eos-ci.yml`. Local runtime: Python 3.12.3, NumPy 2.4.2, SciPy 1.17.0, gfortran 13.3.0. The local generator suite reports an existing Matplotlib Axes3D installation warning; it does not affect these tests. Remaining raw-M0-dependent skips are recorded, not counted as end-to-end validation.

C1 baseline; C2 equation/data tracing; C3 interpolation/variation derivations; C4 independent work allocation; C5 repairs; C6 RED/GREEN regressions; C7 finite-difference and analytic physical checks; C8 runtime and source/output identities; C9 independent code and post-fix counterexamples; C10 this report, committed receipts and handoff. Wolfram independently checked the Maxwell/TOV extremum derivative identities; primary sources were checked through Exa, Parallel and SciSpace.

Verdict: repaired implementation contracts are executable and tested. The scientific project remains PARTIALLY-COMPLETE. Next work is source-calibrated interface functions/parameters, full variational/global-minimum checks, outer-crust construction/splice, composition-aware matching, physical stellar boundary/mode validation, and the remaining EOS390 campaign. Keep PR #15 draft and unmerged.

Primary references: https://arxiv.org/html/2402.15548v1 ; https://arxiv.org/abs/2003.01259 ; https://compose.obspm.fr/download/pdf/manual_v3.00.pdf (cover v3.01).
