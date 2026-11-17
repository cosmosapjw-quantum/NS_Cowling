# Bounded native TF coexistence completion

Base commit: `478c7b5de1514ac512f91156f36551478f36f688`.

This lane finishes the existing finite FSU2R N256 slab/tube chord check and compares the N128/N256 candidates. The original source EOS, all old failures and the production physics remain unchanged. The live calculation is in `data/research_matching_20261117/tf/reproduction/`; no N256 FSU2R admission is made until its retained failure has been diagnosed and the bounded chord re-solves finish.

## Definitions and protected scope

The native T=0 finite-range RMF TF functional uses MeV and fm, with `hbar*c = 197.3269804 MeV fm`; densities are fm^-3, energy density and pressure MeV fm^-3, and baryon chemical potential MeV. Each branch is beta equilibrated and neutral in its cell average. Charge is allowed to vary spatially inside a cell. Macroscopic mixing neglects the interphase surface contribution per volume. There is one effective conserved baryon density in this construction; no two-charge Gibbs or interphase potential solution is asserted.

For high-density volume fraction chi, fixed mean density gives
`n = (1-chi)n_A + chi*n_B` and
`epsilon_mix = (1-chi)epsilon_A + chi*epsilon_B`.
Stationarity with respect to both endpoint densities and chi gives
`epsilon_A'(n_A) = epsilon_B'(n_B) = mu_t` and
`P_A = P_B = P_t`, where `P = n*mu - epsilon`.
Thus `epsilon_mix(n) = mu_t*n - P_t`. The actual native re-minimization stencil computes derivatives of energy; it is not just a substitution of the Euler pressure identity.

The retained root tolerances are `abs(delta_mu)<2e-6 MeV` and `abs(delta_P)<2e-7 MeV fm^-3`, verified again at half the density step. Native SLSQP uses `ftol=1e-14`, `maxiter=1800`; accepted nonuniform states require optimizer success, chemical KKT residual below 0.02 MeV, radius stationarity below 2e-4 MeV and neither radius nor density at its search limit. No tolerance is widened or failed status relabeled.

## Reproduced IEEE flag

The previous failed run printed `IEEE_INVALID_FLAG` without a location. Diagnostic-only replay clears and reads that flag around each native solve. Both initial slab and tube return `STATIONARY`, optimizer mode 0, no search bounds and `IEEE_INVALID=true`; their state rows are byte identical to the old passing rows.

A separate executable built with `-ffpe-trap=invalid` on the retained first passing state stops at `fortran_eos/vendor/slsqp/slsqp_optmz.f:787`. The vendor code explicitly computes `0/0` to mark unused bound multipliers as NaNs. The stack and expected exit 136 are in `data/research_matching_20261117/tf/ieee_probe/`. This identifies the flag source in the probe; the flag alone is not evidence of a failed physical field evaluation. The vendor source is unchanged.

## Evidence separation

- Source fidelity: all eight producer files in the prior TF execution record have identical SHA256; see `SOURCE_IDENTITY.json`. This preserves the previous implementation's identity and does not newly validate every source-paper prescription.
- Physics support: the single-baryon, separately neutral construction and energy-derivative checks apply only to the specified finite branches and assumptions.
- Numerical support: the root, h/2 stencil and five chord sample densities must be reported with native statuses. An N128/N256 change is an observed mesh sensitivity, not a rigorous continuum error bound.
- CompOSE comparison: no new CompOSE comparison is performed by this lane; source tables and prior comparison evidence stay separate.

The emitted research thermodynamics deliberately omit mixed-phase composition and mu_Q. Presence masks remain false. Production registries and raw nonmonotone source tables are not replaced.

## Reproduction commands

Use a fresh output directory because the chord driver deliberately refuses to replace retained state files. For the full diagnostic replay, compile the existing TF modules using `fortran_eos/tests/blocker_tf_build.sh` in an isolated build directory, then compile `fortran_eos/tests/tf_finish_diagnostic.f90` against those module objects and `slsqp.o`, linking LAPACK/BLAS. Its arguments are `GPPVA_FSU2R 256` and stdin is the preserved `data/research_coexistence_20261116/tf/input_n256/FSU2R.txt`. The original failure is expected to remain a failure; mode and bounds explain why.

The independent finite-chord driver is `fortran_eos/tests/tf_finish_chord.f90`, linked against the same unchanged module objects. Its arguments are `GPPVA_FSU2R 256` and stdin is the preserved `data/research_matching_20261117/tf/continuation/input.txt`. The first line supplies the previous native root `(n_A,n_B,mu_t,P_t)`, followed by two retained native stationary seed profiles. All ten chord states are new native solves. Subsequent sample points use the preceding converged profile on the same branch. This is an explicit initialization change, not a relabeling of the old failed state.

`fortran_eos/tests/tf_finish_diagnostic_build.sh BUILD` builds a separate diagnostic-only copy of the optimizer adapter and `tf_finish_single_state`; canonical optimizer and physics sources are never edited by it. The copy adds terminal SLSQP merit, equality/bound, iteration/reset and step diagnostics. The single-state program accepts `MODEL N TARGET_DENSITY`, reads one retained profile from stdin, and fails closed on a nonstationary or different-geometry result.

Set `OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1` for every native invocation. The Python summary and exporter frontend is `scripts/tf_finish_summarize.py --build BUILD`; it checks finite root/h/2 values, exact root parameter binding, stationary native chord states and positive interior sampled gaps before invoking the existing Fortran chord exporter. The CSV is native-written and native-read back. Python does no EOS calculation.

## First actual optimizer diagnosis

The independently initialized first tube chord state also returned `REQUIRES_REFINEMENT` and stopped with exit 3. Its optimizer mode is 8, meaning a positive directional derivative in the SLSQP line search, with success false and both search-bound flags false. The finite returned tube has chemical KKT residual 2.1382e-5 MeV, radius stationarity 3.0324e-9 MeV, baryon error 7.7161e-15 fm^-3 and charge error 8.0976e-19 fm^-3. The rejection is an optimizer termination condition, not a changed topology or the printed physical stationarity criteria. This attempt is retained in `continuation/`, never overwritten. Internal SLSQP diagnostics and a new native minimization initialized from this returned finite state are recorded separately.

## Mechanism and actual refinement

An instrumented diagnostic-only adapter repeats the alternative first-tube call with exactly the same printed returned state. In the native variables, the normalized baryon equality is `c = sum_i sqrt(w_i)*(x_n,i+x_p,i)-1`; physically it is the cell-average baryon-density error divided by the target density, up to floating-point summation differences. SLSQP reaches its six-reset positive-direction branch with `f=f0`, zero bound violation and `h3=1.0502709812953981e-13`. The stock source's final relaxed convergence test requires `h3 < tol`, where `tol=10*ftol=1e-13`. This equality test is the failing term even though energy change is zero and the TF physical residuals satisfy their criteria. This identifies finite-precision feasibility drift at stringent optimizer termination; it does not establish absence of the physical branch.

A new native minimization from that retained finite state actually returns mode 0, after 41 iterations and 286 evaluations. Its normalized equality residual is `9.5701224722688494e-14`, below the same unchanged threshold, with zero bound violation. The native epsilon change is about `-7.24e-13 MeV fm^-3`. Chemical KKT and radius stationarity remain below their original gates. This is a fresh optimization with reset SLSQP state, the same type of warm refinement already used by the native family search. The failed state is never reclassified based on its small KKT residual. Both complete logs and producer/input/binary receipts are in `mode_probe/` and `restart_probe/`.

The repaired chord lane initializes at the actually stationary low-density slab and newly stationary low-density tube, then walks the fixed five sample densities using native forward continuation. Its raw results are in `continuation_repaired/`; the earlier `continuation/` failure and original old failure remain separate. No production formula, stopping tolerance, global phase-selection rule or presence mask is changed.
