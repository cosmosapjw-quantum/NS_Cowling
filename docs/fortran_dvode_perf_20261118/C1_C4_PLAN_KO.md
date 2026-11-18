# C1–C4: native Fortran performance and physics correction contract

## C1 Fixed baseline
- Git source HEAD: `990d07f4f0383f3777f65ed8d3454263b09e14f2`, isolated branch `research/fortran-dvode-perf-accuracy-20261118`.
- PK1 native homogeneous source `data/native_dense_interpolation_20261118/rmf/PK1_CANONICAL/table/compose`; 4097 direct native rows, `P_c=73.52922470749803 MeV fm^-3`.
- Existing Fortran TOV/Love 4097 rows: `M=1.8614087115102154 M_sun`, `R=13.562322006144168 km`, `k2=0.0032945714419691688`, `Lambda=6.424185513339666`, 84,584 radial steps, 3× elapsed wall 0.83–0.89s.
- Existing 1025 rows wall 0.24–0.26s, 257 rows wall 0.08–0.09s. The k2/Lambda are BUGGY physics baselines, not physically approved.

## C2 Concrete formula / implementation map
- `fortran_eos/src/compose_eos_tov.f90::branch_interp` uses linear O(N) search repeatedly.
- `integrate_branch` has a Dormand-Prince error estimate but accepts every step; missing k7 error coefficient. This violates its own advertised adaptive RK45 accuracy.
- In `tov_rhs` the Q expression multiplies by `cs2=dP/dε`; Hinderer requires `dε/dP=1/cs2`.
- In `love_number_k2` an algebraic 1/2 factor appears where the standard Hinderer expression requires 8/5. The existing Python m6 reference repeats both errors and is **not** an independent scientific reference.
- Python `stellar_sharp.py` still uses DOP853 for Cowling. Native-only release cannot claim full Cowling migration in this scope.

## C3 Numerical plan and fidelity gates
1. Keep the same EOS values, units and source phase intervals; make EOS bracket O(log N).
2. Repair Q and Love normalization against independent analytic identity and Fortran native tests.
3. Implement actual error-based step rejection where the native RK integrator is used, and preserve exact Maxwell discontinuity and surface event.
4. Integrate Netlib DVODE as a selectable native Fortran TOV/Love option (nonstiff Adams MF=10 and BDF MF=22), with no Python numerical solver in the production backend. Adaptation must not bypass physical validation.
5. Compare smooth and sharp Maxwell results across tolerances and at matching accuracy. Keep defaults conservative until demonstrated.

## C4 Resources, DAG, permission
- Host CPU=24 logical cores, free available RAM ≈21GiB. No oversubscription or source mutation.
- Sequential file modifications and parallel bounded read-only verification are safe. Netlib DVODE uses COMMON/static state: never run concurrently in threads of one process, only separate processes for independent EOS.
- C5 implementation → C6 build+tests → C7 GR/Hinderer invariants and Maxwell interface → C8 matched-accuracy/real native source wall timings → C9 one adversarial review and targeted correction → C10 decision, provenance, nonforce push, handoff.
- Failures separate physical formula, numerical accuracy, implementation, I/O, and environment. No fraudulent PASS.

### Original methods
Netlib DVODE `https://netlib.org/ode/vode.f` (vendor SHA-256 `76b3858bb1e32f416a29363dcf56e0889ac25aceb5911289f3f25fe275634e95`). ODEPACK overview `https://computing.llnl.gov/projects/odepack`. Hinderer Tidal Love `https://arxiv.org/abs/0711.2420`; Takátsy–Kovács Maxwell correction `https://arxiv.org/abs/2007.01139`.
