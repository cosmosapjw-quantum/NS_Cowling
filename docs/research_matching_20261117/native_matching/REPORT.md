# Native core–crust matching: bounded adversarial repair

Base: `478c7b5de1514ac512f91156f36551478f36f688` in the shared isolated worktree `/tmp/ns_match_478c7b5_20261117`. Source ownership was limited to `fortran_eos/src/compose_eos_matching.f90`; no phase-transition kernel, shared build, commit or push was changed by this lane. The before snapshot is byte-identical to the base commit, verified in `data/research_matching_20261117/native_matching/freeze.json`.

Selected work: C2–C9. Read/applied PhysMath 1.3 main, POLICY_INTEGRATION, SUBAGENT_MODE, CODING_LOOP; systematic-debugging and verification-before-completion. This was the native implementation/reproduction part of the parent's existing requested adversarial pass. Theory and pipeline were complementary sibling lanes, not recursively delegated reviewers.

## Physical contract and code map

At zero temperature, use total energy per baryon F (MeV), n (fm^-3), P and epsilon (MeV fm^-3), mu (MeV). The dimensionless sound speed is s=(c_s/c)^2=dP/d epsilon. The native Fortran backend implements all physical calculations; the new Python runner only compiles and executes tests.

For canonical bridges x=ln(n/n1), q=P/n, F_x=q, mu=F+q. With each local polynomial parameter t and scale Delta x, mechanical stability and causality require D=q+q_t/Delta x>0 and mu-D>=0. Positivity requires q,F,mu>0.

For pressure-potential bridges mu=mu0+Delta mu*u, n=P_mu. With local scale Delta mu, require n_t>0, mu*n_t−Delta mu*n>=0, P>0 and epsilon=mu*n−P>0. All polynomials are degree at most five. `certify_profile` and `certify_spline` now check these inequalities at both endpoints and all stationary points found by recursively partitioning at derivative roots. Tangent stationary roots are retained. This is floating-point polynomial extrema coverage, not rigorous interval arithmetic.

## Reproduced defects and corrections

1. **Acausality between 401 accepted diagnostic nodes.** The original quintic_P constructor and validator accept a curve with maximum sampled s=0.99999854118234377 while s=1.0000000031243748 at u=0.8919032228381553. Endpoints are n1=1,n2=3,mu1=1,mu2=2,P1=.01,P2=2.1041881294006204; endpoint dn/du values are 2.1043168390302784 and 1.5780789687536598. The new whole-piece extrema gate rejects this curve regardless of output node placement. The original 401-point arrays remain diagnostic summaries; their min/max fields are not claimed to be global extrema.

2. **Admissible source Q3 residual changes the physical bridge.** A +1e-4 relative Q3 perturbation, within the existing provider Euler tolerance 5e-4, changes quintic_F output by 3.2905030689533987e-5 in Q and causes cubic_P/quintic_P endpoint first-law failure. The old cubic_F state is unchanged but its reported core chemical potential is shifted. Window search, all endpoint construction and reporting now use mu=F+P/n after the unchanged source Q3 gate. The controlled perturbation changes no output Q in the new four-method test. The Il=0 charge-shift source control also succeeds. Physical energies still decode with each source's own reference mass; a 939-to-950 MeV crust normalization change preserves output states.

3. **Pressure-work integration violates the interior first law.** The old near-plateau pressure construction uses GL16 globally for partial integrals and GL64 at the full endpoint across piecewise-cubic knots. For the recorded spline, a centered derivative of its antiderivative disagrees with P/n by 0.0018127–0.0036057 at u=.3–.7. Integrating each smooth cubic piece separately using the same GL16 rule reduces those residuals to 1.7e-13–2.3e-11; the finite-difference error at u=.1,.2 remains below 1e-9. This is a correction to the retained internal historical pressure-bridge implementation, not an assertion that it has been newly exposed as a public method.

4. **Stable decreasing or constant P/n is treated as invalid.** For x=ln(n/.04), F=940+10x−x²/2 and q=10−x, D=q−1>0 although q decreases. The old F bridge rejects it; the new canonical polynomial paths accept it subject to the full physical inequalities. For q=10 and F=940+10x, `find_window` divides by zero and the original trap-enabled build signals SIGFPE. The flat-q compatibility case now uses a guarded finite ratio and the same certified polynomial forms. Its compatibility threshold accounts only for cancellation roundoff in F2−F1. All four F/P methods succeed for both controls. Adaptive monotone-profile fallbacks still have explicit shape limitations; this patch does not assert completeness of the heuristic window or bridge search.

5. **Unused zero division in P profiles.** The original `bridge_profile` evaluates canonical derivative expressions with canonical_dx=0 even for the P representation. Those expressions now execute only for canonical F profiles. Final controls run with invalid/zero/overflow traps enabled.

## Actual verification

The frozen native matching source SHA256 is `72ff79e1709f51fa479996fd051a3735ca539fea64e213c28fd3fabc884391f2`. `freeze.json` also binds the unchanged types/phase-transition sources, original source snapshot and full isolated API binary.

- Before fixes, both native regression programs exit 1 against the exact preserved source: the thermodynamic program records 16 failed controls; the interval program records the accepted acausal witness. Commands, source/test hashes and outputs are under `data/research_matching_20261117/native_matching/red_reproduced/`. The separate trap-enabled flat-q reproduction is `red_flat_q.log`.
- After fixes, all three programs exit 0 under `-fcheck=all -ffpe-trap=zero,invalid,overflow`: source-Q3 invariance, rest-mass normalization, Il=0 decoding, decreasing/constant q, interior first-law differentiation, acausal off-grid witness, flat quartic stationary root, two minima, narrow negative pocket, endpoint minimum and nonfinite coefficients. Records are under `.../green/` and `final_controls.log`.
- Full isolated build: `OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 bash fortran_eos/build.sh build/matching_native/api`, exit 0. `api_build.log` preserves the existing linker executable-stack warning in the unrelated SPG native object.
- Existing API regression/admission tests: 26 passed in 5.23 s, including all four original legacy differential comparisons without tolerance changes. Command and exact API binary SHA256 are in `freeze.json`; output is `api_matching_tests.log`.
- Sibling theory review inspected the transformed polynomial coefficients, physical inequalities and recursive extrema isolation and reported no algebraic defect. Its six public-API/source/topology tests passed in 0.26 s against the same frozen binary (hashes verified); the log is `data/research_matching_20261117/theory/theory_six_green.log`.

Reproduce the controls from the repo root:

```sh
OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 python3 fortran_eos/tests/run_matching_native_adversarial.py
```

The runner compiles only types and a test-local exact source copy whose sole modification is the private-to-public access declaration for tested kernels. It asserts reversing that declaration recovers the original source text and records both hashes. Production exports remain unchanged. CI can call this runner directly.

## Scope and remaining obligations

The native matching repairs are implementation-verified and numerically checked for the recorded controls and existing scoped regressions. No CompOSE success-rate threshold or universal matching theorem is claimed. Source support gaps and downstream admission belong to the pipeline lane. Maxwell phase topology belongs to the existing separate phase-transition kernel and theory tests; it was not modified here. Retained source-table segments and the broader final application build are still subject to the parent's integration validation. No further native edits are planned after this freeze unless a new concrete regression is observed.
