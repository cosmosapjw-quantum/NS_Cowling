# Coupled vector-field band solve: isolated performance proposal

Date: 2026-11-16. The isolated experiment below initially left production and the running frozen campaign source untouched. After independent verification, the exact proposal was promoted to production; see the cutover record at the end. Baseline `spatial_tf.py` SHA-256: `d0c027c4c388a9b37cbfae6ecae753354fce0ae83bbfa1867aad88227a47a826`; isolated prototype SHA: `2d1f1edc72b170b01194f2705e62905383a2d4292bd84cc58067bc2cc1067c02`. The proposal is in `data/research_finish_20261116/vector_band_prototype/PROPOSED.patch`.

## Algebraic equivalence

The coupled vector Newton matrix, in grouped unknown order, is

    J = [[J_omega, D], [D, J_rho]],
    J_omega = L + diag(W_omegaomega),
    J_rho   = L + diag(W_rhorho),
    D       = diag(W_omegarho).

Here L=(hbar c)^2 (-nabla_x^2)/R^2 is tridiagonal; the mixed derivative is local in the radial cell. Permute both rows and columns to (omega_0,rho_0,omega_1,rho_1,...). Every mixed derivative now lies at offset +/-1 within its own cell, and each same-field spatial neighbor lies at offset +/-2. The permuted Jacobian therefore has lower and upper bandwidth two. `scipy.linalg.solve_banded((2,2),...)` solves precisely the same Newton equations, with the same signs and nonlinear residual.

The implementation uses all five diagonals explicitly and interleaves/uninterleaves the right-hand side and solution. It makes no symmetry or positive-definiteness assumption: the radial finite-volume matrix is generally nonsymmetric in ordinary coordinates. A general pivoted band solver also handles the retained nonsingular zero-leading-pivot example; an unpivoted block-Thomas shortcut would require extra assumptions. Exactly singular matrices retain `LinAlgError` behavior. The zero-cross production branch, scalar field, nonlinear line search, equations, field residual thresholds, energy, gradient and pressure definitions remain unchanged.

The *linear solve* changes from dense cubic work and quadratic storage to fixed-band linear work/storage. The full solver still has dense-grid Poisson and field-action operations and a constrained density optimizer, so no claim of linear scaling for the full EOS calculation is made.

## Actual profile benchmarks

The source archives are the retained final GPPVA_FSU2R, GPPVA_FSU2H and GPPVA_TM1e cells nearest n_B=0.05 fm^-3. Their retained profiles have N192 radial cells. N384 benchmark profiles are explicitly labelled linear resamples of those profiles; they are not new stationary N384 EOS results. Each model/mesh was evaluated from cold meson fields and from retained fields perturbed by factors (0.9,1.3,0.7) for (phi,omega,rho), giving 12 paired cases. Both full field/functional evaluation and the screened-electron reduced objective with its gradient were timed five times with one OpenBLAS/OMP thread on the currently loaded machine.

| Grid | Full field evaluation median speedup range | Reduced objective/gradient median speedup range |
| --- | ---: | ---: |
| N192 | 2.60–6.93 | 1.42–3.59 |
| N384 | 4.24–18.49 | 3.49–8.33 |

These ratios are local timing observations, not whole-minimization speed predictions. Across all 12 paired cases, maximum omega/rho differences are 4.27e-14/1.47e-14 MeV; phi and its residual are identical. Energy-density and fixed-profile pressure differences are at most 7.11e-15 MeV fm^-3. The reduced gradient differs by at most 3.29e-14 MeV and the objective by 2.85e-13 MeV. Every final field PDE residual remains under the original criterion. The fixed-profile pressure is a diagnostic comparison on the same input, not an assertion of equilibrium for resampled profiles.

Five additional matrix adversaries cover nonsymmetric, indefinite, zero-leading-pivot, zero-cross and exactly singular matrices. Dense and banded solutions agree within the recorded tolerances, and both raise `LinAlgError` for the exactly singular case. Independent review separately assembled and solved six physical Jacobians and verified the index permutation and unchanged scalar/nonlinear/pressure semantics.

`BENCHMARK.json` retains every timing sample, differences, residuals, seed archive hashes and source hashes. `benchmark_vector_fields.py` reproduces the paired experiment. The test file `python_frontend/tests/test_spatial_tf_vector_band.py` adds dense-oracle comparisons across dimensions 1,2,3 and grid sizes 8,17, plus pivot/indefinite/zero-cross/singular cases. Its baseline failure only means the proposed helper is absent; it is not evidence that the old dense physical solution is wrong.

## Actual N384 minimization

`run_fsu2r_cell.py` runs one actual FSU2R n_B=0.05 fm^-3 N384 minimization from the retained N192 profile with the isolated prototype. It preserves the resulting full profile, actual optimizer/KKT/radius/energy/pressure results, wall time and source/input identities in `FSU2R_N384_CELL.json.gz` and `FSU2R_N384_RECEIPT.json`. It then independently evaluates that final profile using the original dense and prototype band field operators. No duplicate dense full minimization is launched. The receipt, rather than the fixed-profile benchmarks, determines whether the proposal is ready for a coordinated production cutover.

The actual isolated run completed in **149.90 s**, with optimizer success after 124 iterations/312 evaluations and `STATIONARY` status. E/A=947.522780879533 MeV, pressure=0.177261626421206 MeV fm^-3, Y_p=0.0506037906076081. The chemical KKT residual is 3.335e-5 MeV, radius stationarity 3.212e-8 MeV, electron PDE residual 4.90e-13 MeV fm^-2, baryon error 8.53e-15 fm^-3 and charge error 3.22e-20 fm^-3. Dense and banded reevaluations of this final profile give identical float64 energy density and pressure, with meson-field differences no larger than 3.20e-14 MeV. Their independently restarted pressure differs from the optimizer's cached-field report by 2.8e-11 MeV fm^-3, consistent with finite field-solver tolerances. All original physical acceptance thresholds are retained.

An isolated process substituted only the saved prototype module for the production import and ran the new 10 helper tests plus the 53 field/electron regressions: **63 passed in 14.42 s**. The benchmark and actual minimization both verified production source hashes remained unchanged during execution. These results support a coordinated production cutover; they do not provide a measured speedup for a complete dense-versus-band minimization, because the expensive dense optimization was deliberately not duplicated.

Independent fresh cold-field reevaluation of the actual N384 result also passes: energy density, pressure and mu_B are bit-equal for dense and band operators; maximum meson and chemical-potential differences are 3.20e-14 and 5.69e-13 MeV. Scalar/vector PDE residuals are 1.083e-11/2.16e-14 in the module's normalized field-residual units. The independent receipt is `data/research_finish_20261116/independent_review/VECTOR_FINAL_PROFILE_REVIEW.json`.

## Validated production cutover

After the isolated full N384 solve and independent dense-operator reevaluation passed, the exact prototype was installed as production `spatial_tf.py`, SHA-256 `2d1f1edc72b170b01194f2705e62905383a2d4292bd84cc58067bc2cc1067c02`. The conservative electron implementation remains SHA-256 `8be3c6756f1b04a38cb59fdeeacb5c0aed5c98d8337c8ed14c64e7091501bdca`. The 112 completed mesh-density nodes from the repaired dense-vector campaign are retained under their original source identities; the remaining 40 N384 nodes are computed with this band implementation. The interrupted dense runs and per-start attempts are retained. These are versioned calculations of the same functional, not a claim that every raw file was generated by identical code bytes.

A fresh current-operator evaluation of 1,254 retained profile attempts (1,241 unique profiles) found maximum differences of 3.41e-13 MeV/baryon in E/A, 2.84e-12 MeV fm^-3 in pressure and 2.90e-11 MeV in mu_B, with no evaluation failures. This fixed-profile comparison establishes the stated numerical reuse scope; it is not a second dense minimization or a proof of identical optimizer trajectories. The final combined regression with the production band solver passed 971 frontend/pipeline tests, with nine skips, and 142 generator tests. Exact commands and before/after source hashes are in `data/research_finish_20261116/integration/REGRESSION.json`.
