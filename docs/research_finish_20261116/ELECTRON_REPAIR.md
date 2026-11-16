# Screened-electron repair: exact trial replay and conservative Newton solve

Date: 2026-11-16. Base commit: `573a31485bc9af30915e01fe74f259e4fa5ac54f`.
Final production module SHA-256: `8be3c6756f1b04a38cb59fdeeacb5c0aed5c98d8337c8ed14c64e7091501bdca`.

The repaired electron reduction completes the two previously failing DD2 rod starts at n_B=0.0625 fm^-3 and six newly exposed TW/DDME2 slab starts at n_B=0.06 fm^-3. All eight complete cell minimizations satisfy the existing `STATIONARY` classification; this is finite-start, finite-grid evidence, not proof of a global minimum or a complete EOS. No RMF coupling, energy definition, electron equation of state, source table, or physical acceptance threshold was changed. The work is implementation-verified on the retained trials and tests; the algebra below is derived independently of agreement with CompOSE.

## Two distinct numerical failures

The baseline module is preserved byte-for-byte as `data/research_finish_20261116/electron_repair/baseline_spatial_tf_equilibrium.py`, SHA `559f62e1fec635e38d1a098d0d1d7831b6e57c83ed4672c0b18a141609c34a79`. `reproduce_failure.py` imports this saved source and replays the original warm and cold DD2 starts. The warm seed is the retained N192 rod at n_B=0.065 fm^-3, scaled to 0.0625; the cold start uses R=30 fm. The seed archive SHA is `5b607c6e0d5945ac01a0534291b0e42f4c56ab2f035b2808005ac38609992aef`.

Both failures occur in intermediate SLSQP proton profiles, not in an accepted EOS state. Their mean proton densities are 2.0475e-17 and 2.4420e-16 fm^-3. An unconstrained Newton update takes the electron chemical energy below m_e c^2; every electron susceptibility becomes zero. Inverting the Neumann Laplacian then attempts to solve its constant null mode. The captured erroneous chemical-potential steps are 0.1394 and 1.8705 MeV, far larger than the near-threshold interval. Increasing the list of damping factors alone does not remove the singularity.

The first repaired solver (saved as `neutral_pinned_spatial_tf_equilibrium.py`, SHA `11af1e007d997d9eea7465c4d42552440544e26f7070b0b3bd411f13990cb3f9`) fixes that defect. A subsequent all-model campaign exposed six slab failures with mean proton densities 0.00177–0.006 fm^-3, where the dense matrix evaluation of the Poisson action prevented further Newton progress near the internal 2e-12 MeV fm^-2 target. The retained residuals are 2.05e-12–4.46e-12. Independent 70-digit evaluation shows that these old iterates really still have residuals of that scale; they cannot simply be reclassified as converged. Evaluating the same operator by conservative face differences enables further Newton updates. The six slab records contain five distinct electron inputs: the cold TW and DDME2 electron inputs coincide because the starting proton profile is model-independent.

## Equations, gauge and constrained Newton step

Use finite-volume weights w_i>0, sum_i w_i=1, and the discrete positive operator L=-nabla^2/R^2 with no-flux boundaries. The grid identity is w^T L=0. V and mu are in MeV; number densities are in fm^-3; kappa=4 pi e^2 has units MeV fm. Define q_i=mu+V_i and

    n_e(q) = [(sqrt((q-m)(q+m)))/(hbar c)]^3 / (3 pi^2), q>m,
    n_e(q) = 0, q<=m,                  m=m_e c^2.

The product form avoids subtracting two nearly equal squared energies. Its active derivative is chi(q)=q sqrt((q-m)(q+m))/(pi^2 (hbar c)^3). The equations are

    L V + kappa (n_e(mu+V)-n_p) = 0,
    C = sum_i w_i (n_e(mu+V)-n_p) = 0.

Pin V_N=0 during Newton and impose the first N-1 Poisson rows plus neutrality. Weighted summation then recovers the omitted last Poisson row, so this is a gauge choice, not a changed boundary condition. At the end subtract <V> and add it to mu, preserving q.

Let A=L_11+kappa diag(chi_1), z=A^-1(-r_1), and u=A^-1(kappa chi_1). Linearizing charge gives

    delta_mu = [-C - sum_{i<N} w_i chi_i z_i]
               / [sum_{i<N} w_i chi_i (1-u_i) + w_N chi_N],
    delta_V_1 = z-u delta_mu,       delta_V_N=0.

The pinned Laplacian block is nonsingular even as chi tends to zero. Its M-matrix structure gives 0<=u_i<1; the denominator is positive when at least one susceptibility is positive. The implementation evaluates the denominator in the displayed form, avoiding subtraction of two large nearly equal susceptibility sums.

At every proposed potential, solve scalar neutrality again before evaluating the line search. If pbar=<n_p>>0, a valid bracket is mu_low=m-max(V), giving no electrons, and mu_high=hypot(hbar c (3 pi^2 pbar)^(1/3),m)-min(V), giving at least pbar electrons in every cell. A small outward numerical margin protects the upper bracket. Scalar neutrality thus prevents the zero-susceptibility state. The exactly zero-charge case returns zero densities/potential explicitly. Invalid shapes, nonfinite/negative proton densities, and nonpositive/nonfinite radii are refused.

An independent variational check uses J(q)=<q L q>/2+kappa<P_e(q)-n_p q>, with dP_e/dq=n_e. Its Hessian is L+kappa diag(chi), positive in the weighted inner product when some electrons are active. Neutrality elimination minimizes J along the constant mode. This supplies an independent consistency check on the Schur complement and gauge construction.

## Conservative evaluation and precision

For face i+1/2 define c_i=-w_i L_{i,i+1}>0 and F_i=c_i(V_i-V_{i+1}). Evaluate (LV)_i=(F_i-F_{i-1})/w_i, with both boundary fluxes zero. This is the same intended finite-volume operator and retains its constant null mode without multiplying and subtracting large diagonal/off-diagonal products. The Newton matrix remains the existing banded finite-volume matrix.

An extra stress check interpolated the six slab records onto N=192,384,768 grids. Float64 potential bookkeeping with face fluxes passed every N192/N384 check but stalled in four N768 records. This is a numerical-resolution issue, not a physical exclusion of those models. Keeping only the internal potential iterates in `np.longdouble` resolves all 18 checks on the tested Linux platform. Banded linear solves and public outputs remain float64. After rounding the public potential and chemical potential, the function recomputes the returned densities and PDE residual and applies the original physical acceptance criteria. The largest rounded-output residual in this stress check is 3.83e-12 MeV fm^-2, below the unchanged 2e-10 criterion. Extended precision availability is platform-dependent; these checks do not prove support for every mesh or every representable input.

The ordinary internal target remains 2e-12. For nearly empty positive proton profiles it is tightened to min(2e-12,kappa pbar*1e-6); an adversarial cached-potential test exposed the need to prevent an absolute residual floor from hiding relative charge-profile error. The final physical thresholds remain PDE<=2e-10 MeV fm^-2 and absolute mean charge<=2e-11 fm^-3. The cell KKT, baryon, radius and optimizer criteria are unchanged.

## Reproducible evidence and scope

The committed evidence directory contains original inputs, full failed iterates, source snapshots, red logs, intermediate proposals, final profiles and receipts. `RESULT.json` indexes final outcomes and `SHA256SUMS.txt` identifies every evidence file. Intermediate `INITIAL_*`, `PINNED_*`, `FLUX_*` and `EXTENDED_*` files are explicitly historical or prototype evidence; only `GREEN_RECEIPT.json`, `REPAIRED_CELLS.json.gz`, and `FINAL_SLAB_CELLS.json.gz` report final production replays.

Final targeted command:

```sh
PYTHONPATH=python_frontend OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 python3 -m pytest -q \
  python_frontend/tests/test_screened_electron_vacuum_trials.py \
  python_frontend/tests/test_spatial_tf_equilibrium.py \
  python_frontend/tests/test_spatial_tf.py
```

Result: **53 passed in 4.43 s**. This includes the actual eight failure records, dilute/stale/gauge cases in all three dimensions, invalid input refusal, independent integrated-charge/flux checks at N192/N384/N768, and existing full-functional gradient and pressure-derivative checks. Historical red runs retain 2 failed original DD2 tests, 6 failed slab tests, and 2 dilute-profile stopping failures before the corresponding repairs.

At N192, final DD2 warm/cold E/A values are 946.4372194893654 and 946.4372194893638 MeV; chemical KKT residuals are 1.81e-4 and 1.14e-4 MeV. TW three-start E/A spread is 8e-13 MeV, and DDME2 spread is 4.2e-12 MeV. The least accurate accepted DDME2 warm start has chemical KKT residual 0.00255 MeV, still under the existing 0.02 MeV criterion; numerical agreement of energies must not be mistaken for a stronger stationarity bound.

Independent review retains a 70-digit scalar-neutrality/flux-integration oracle for the two original near-vacuum inputs and a separate 70-digit dense-versus-flux comparison for the six slab records. The original pinned repair agrees with the near-vacuum density oracle to 1.91e-10 relative or better, and its potential to 4.4e-25 MeV or better. The final-source near-vacuum oracle recheck reproduces those error bounds. A separate 70-digit exact-geometric residual check after re-solving the six slab records with the final source finds maximum PDE residual 2.733e-13 MeV fm^-2 and relative neutrality error 6.205e-16 on the rounded float64 outputs. See `independent_review/ELECTRON_ORACLE.json` and `independent_review/FINAL_DENSE_ELECTRON_ORACLE.json`.

No geometry transition, continuum convergence, source coefficient calibration, low-density crust join, complete model census, core–crust join, or stellar result is admitted by this repair alone. Those claims belong to separately retained campaigns. CompOSE agreement is not used as a solver pass criterion.
