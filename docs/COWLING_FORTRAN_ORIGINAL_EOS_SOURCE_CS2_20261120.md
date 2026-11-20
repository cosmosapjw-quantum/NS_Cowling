# NS_Cowling: original-source equilibrium c_s² in native Fortran

Date: 2026-11-20 KST.
Baseline HEAD: `19e5b16cd0c5542e8fcb91f2fab14f36b1e663df`.
Isolated branch: `research/cowling-fortran-exact-source-cs2-20261120`.
PhysMath: bounded Fortran source-law implementation and scientific verification;
the 10-stage track is not implicitly restarted to exhaustion.
Native skill-bound agent dispatch: **NOT_RUN**, because the host offers no
native child spawn/delegate action. Work was serialized in one scientific
implementation session, while independent *numerical processes* were measured.

## 1. Exact referent, physics conventions and scope

Source EOS bundle:
`data/eos_missing_data_recalc_20261119/artifacts/EOS_ALL_REQUESTED_SOURCE_LOCKED_NESTED_20261119.zip`,
SHA-256
`299de43bb4f0a87ebfc4f1139ab3c8b51c344937ec3637cdff644e9ca8ab2c04`.

One nested `05_SPG5_PK1_10_WHOLE_STARS.zip` contains:
- `FRESH_5x2_WHOLE_STARS/SPG_M2/maxwell/REPRESENTATION.json`,
  SHA-256 `59b7cb4671347c33c25133e6497e8be346b4e5b4f871eece9bcaf59155fd8765`.
- `FRESH_5x2_WHOLE_STARS/SPG_M3/css/REPRESENTATION.json`,
  SHA-256 `f23fb7020c9de1b65c99ac37690840f7a4d6611f3372409b015d8c3bd7e32f3b`.

These represent **the exact piecewise mathematical laws of the source
EOS contract**, not an independently reconstructed microscopic EOS.
The contract itself and the frozen source archive were not modified.

Signature (-,+,+,+); r,m in km, n in fm^-3, chemical potential mu and
F=epsilon/n in MeV, pressure/epsilon in MeV fm^-3,
`h=log(mu/mu_rest)`, and mode omega_geom=2*pi*f/(299792.458 km/s).
The oscillation closure remains `Gamma1=Gamma_eq`: **no
frozen-composition gamma1, bulk composition g-modes or full GR damping
validation** is claimed.

## 2. Definition to Fortran computational form

The new independent module
`fortran_eos/src/compose_eos_cowling_source_exact.f90`
implements the original source's four equilibrium derivative laws.
It computes (c_{s,eq}^2/c^2=dP/d\epsilon) rather than interpolating
c_s² from a coarse PCHIP table. For chemical potential mu:

- Surface polytrope `epsilon=m*n+P/(gamma_s-1)`:
  `c_s²/c² = (gamma_s-1)*(mu-m)/mu`.
- Pressure-primary outer polytrope, source row (n_i,P_i,F_i,gamma_i):
  `c_s²/c² = ((gamma_i-1)*(mu-F_i)+P_i/n_i)/mu`.
  This form is stable at gamma near 1 and requires no nonnative
  density/pressure inversion.
- Causal-exponential bridge
  `y(mu)=n/mu=y_L+Delta_y*Phi_k(t)`,
  `t=(mu-mu_L)/(mu_R-mu_L)`,
  `Phi_k(t)=expm1(k*t)/expm1(k)`:
  `c_s²/c²=y/(y+mu*Delta_y*Phi_k'(t)/(mu_R-mu_L))`.
  A direct Fortran small-argument expm1 series preserves cancellation.
- Upper source Hermite cubic in `P(mu)`:
  `n=dP/dmu`, `chi=dn/dmu=d²P/dmu²`,
  `c_s²/c²=n/(mu*chi)`. Every upper interval is evaluated
  analytically using its source endpoints, never a cross-interval
  sound-speed PCHIP approximation.

The Maxwell source has **two distinct finite one-sided states**
at the same pressure and chemical potential. At the coexistence
band the Fortran evaluator selects the low/high adjacent smooth
Hermite interval according to the explicit phase-side argument.
The original source includes a last-bit M2 pressure doublet
difference `-1.734723475976807e-18 MeV/fm³`; the same
preexisting 2e-12 relative *coexistence-only* continuity
criterion is used for this admitted plateau. All other upper
pressure intervals are required to rise strictly.
No physical threshold has been relaxed to hide failure.

Wolfram algebraic simplification independently produced
`{0,0,0}` for the outer, bridge and upper c_s² identities.
The Fortran formulas were also checked against the actual
original Python source evaluator at representative interiors,
both Maxwell limiting sides and near the EOS knots.

## 3. Fortran integration and computational architecture

The module provides C-compatible original EOS setup,
`cowling_native_source_cs2_eval`, phase-local geometry setup,
and native register/unregister of a **Fortran-to-Fortran** source
coefficient callback.

For a selected one-Maxwell `native_source_cs2` run:
1. Python's frontend reads and verifies the original JSON contract,
   passes immutable surface/outer/bridge/upper arrays to Fortran.
2. Fortran loads and validates all analytic source laws exactly once
   per independent matching computation.
3. For each phase, a native process-local cached Hermite PCHIP lookup
   evaluates **only background geometry**, while the fifth
   coefficient `c_s²` is replaced at each DVODE RHS call by
   the exact contract formula at the instantaneous enthalpy.
4. The legacy phasewise DVODE path restarts on the **same preselected
   strong EOS knot boundaries** used by the original Python
   `hybrid_cs2` callback baseline. Maxwell slow/rapid interface
   transfer remains the previously validated Fortran routine.
5. The Fortran callback is unregistered in a `finally` cleanup
   in the frontend and is never shared between OS processes.
   The existing global DVODE COMMON and SAVEd buffers remain
   explicitly **nonreentrant**.

New backend:
```python
solver = NativeSharpCowlingDVODE(
    background, settings,
    samples=4097,
    coefficient_backend="native_source_cs2",
    source_restart_relative_jump=0.01,
)
mode = solver.solve_mode((3710., 3760.), conversion="slow")
```

Python is the frontend for new source-coefficient evaluation: there
are **zero Python RHS EOS callback calls** in this path.
Remaining numerical components that are NOT yet Fortran-only:
EOS source-contract input marshaling, one-sided background/PCHIP
coefficient preparation, identification of selected strong derivative
knots in Python setup, and the existing sharp-mode SciPy Brent
root-polishing outer loop. Root evaluation calls the native Fortran
DVODE. Do not advertise a complete Fortran-only application.
The prior onejump `fortran_phase_pchip` backend remains a separate,
explicit, less source-faithful path.

## 4. Original-source numerical tests

Direct source-law test:
`python_frontend/tests/test_fortran_source_exact_cs2_20261120.py`.

**8 passed** for M2 and M3 combined, including:
- all four analytic sectors and source domain;
- low/high Maxwell one-sided values;
- perturbations to h on either side of sampled source knots;
- fail-closed invalid argument and deliberately corrupted source.

Fortran-only standalone test:
`fortran_eos/tests/cowling_exact_source_cs2_test.f90`.
Its source inputs are immutable exported original EOS parameters plus
an independent original-Python-source value oracle:
`FORTRAN_EOS_SPG_M2_maxwell.dat.gz`,
`FORTRAN_EOS_SPG_M3_css.dat.gz`.

Fortran source evaluator comparisons:
- M2 **178** probes, M3 **177** probes.
- Across each model the measured maximum normalized relative
  discrepancy was `8.53284e-12` on those probes.
- The self-contained standalone Fortran executable's
  memory peak was **5,356 KB** for M2 and **5,096 KB** for M3.
  The Fortran standalone evaluates source thermodynamic states,
  not whole TOV/Cowling end-to-end physics, so compare memory
  only within its own task.

Compilation: changed project Fortran modules and test built with
`gfortran -O0 -g -fcheck=all`. The untouched legacy Netlib
ODEPACK units were linked with `-O2 -std=legacy`, avoiding the
known incompatible whole-vendor F77 DGEFA bounds instrumentation.
This is **not a whole-vendor memory safety guarantee**.

Reproducible build:
`scripts/build_fortran_source_exact_cs2_test_20261120.sh`.
Frontend-only input exporter:
`scripts/export_fortran_exact_source_cs2_20261120.py`.
The zipped inputs are deterministic `gzip -n -9` outputs,
and the reconstructed Fortran scientific result lines matched
the uncompressed originals byte-for-byte.

## 5. Matched real-star source-derivative mode comparison

Two actual source-locked Maxwell/CSS equilibrium barotropes,
the same native TOV background, phasewise PCHIP **geometry**,
DVODE tolerance, strong knot restarts, Maxwell junction, and
Brent outer loop are used in both variants.
Only the (c_s²) RHS evaluator changes.

| Case | Existing `hybrid_cs2` Python callback frequency | `native_source_cs2` frequency | Native minus callback |
|---|---:|---:|---:|
| M2 Maxwell 1.4 Msun rapid p1 | 5143.719666721 Hz | 5143.719608336 Hz | -0.000058385 Hz |
| M3 CSS 0.7 Msun slow p1 | 3736.210673258 Hz | 3736.210698740 Hz | +0.000025482 Hz |

Both return 1/1 phasewise W,V node counts with
`relative_join_mismatch <1e-5`, and retain interface
traction invariants. The native (c_s²) results differ from
the separate higher-accuracy exact-original-DOP853 references
by `-0.000113227` Hz (M2 rapid p1) and
`+0.001803582` Hz (M3 slow p1).
Neither is an absolute numerical precision certificate.

The exact original source derivative is evaluated without Python
round trips, and the strong-knot restart counts remained
**39** for M2 and **40** for M3, identical to the
previous `hybrid_cs2` backend. The M3 slow case is not
admitted via the inaccurate PCHIP-only source derivative:
it is specifically the Fortran source-exact path, with a
join residual of `2.455e-6`, that passes.

New integrated star tests:
`python_frontend/tests/test_cowling_native_source_cs2_modes_20261120.py`,
**4 passed** (including two actual slow/rapid
source-EOS eigenmode parity calculations).
All test oracles are source-locked.

## 6. Performance, memory and process parallelism

Paired tests reuse identical source-model background and
prepared phase geometry tables. Three warm full-mode solve
observations per backend were taken:

| Case | Old callback median | Fortran analytic source median | Matched-accuracy speedup |
|---|---:|---:|---:|
| M2 rapid p1 | 4.15456 s | 0.02966 s | **140.1x** |
| M3 slow p1 | 6.67482 s | 0.04529 s | **147.4x** |

The eliminated *Python source callback* call counters
were **308,464** and **496,928**, respectively, over the
warm-up and repeated calculations logged in the benchmark.
Fortran backend counters are zero Python callbacks.
Existing model background/PCHIP setup, not shown above,
took about 1.42 s for M2 and 0.85 s for M3 for the
new backend, which is excluded from those warm medians.
These are strong **scoped observations**; no universal
performance factor across all source models or phases
is inferred.

The new source coefficient function uses a cached Fortran
phase-grid cell hint with a binary-search fallback.
As before, four original equilibrium laws are evaluated
without cross-derivative interpolation at EOS knots.
The new Fortran source parameter buffers allocate only
process-local O(N_source) arrays.

Independent real-star process runs, including full EOS loading,
native TOV, coefficient preprocessing and native physics:
- 1 worker for M2 + M3: **7.556697 s** total;
- 2 workers in separate OS processes: **4.712607 s** total;
- speedup **1.603x**, while the canonical per-star JSON outputs
  and SHA-256 manifest are byte-identical.

Maximum RSS observed across independent full-front-end
workers was approximately **101--109 MB per process**.
No full-pipeline memory reduction is claimed.
The original callback path vs new source evaluator could
not be treated as equivalent standalone source-law memory
measurements.

One first process benchmark failed at **Python output
serialization before any numerical solver** with
`SyntaxError: keyword argument repeated: join`.
This is classified `IMPLEMENTATION_FRONTEND_CLI_SYNTAX`,
not physical or numerical failure; raw logs are preserved.
After fixing the duplicated metadata key, the fresh
serial and parallel jobs completed and were byte-identical.

Direct performance evidence:
`data/cowling_fortran_source_cs2_20261120/REAL_STAR_PERFORMANCE.json`;
process output directories
`REAL_P1_VERIFIED/` and `REAL_P2_VERIFIED/`.

## 7. Integrated regression and acceptance

Fresh related regression suite across current new source tests,
C7/C8/C9/C10, the prior PCHIP onejump module, native smooth
Cowling and TOV, sharp Maxwell hybrid, direct original-source
callbacks and empirical mode selector:
**95 passed, 0 failed in 179.10 s**.
Raw test output:
`data/cowling_fortran_source_cs2_20261120/INTEGRATED_REGRESSION_FINAL.log`.

Accepted:
- exact *contract* c_s² laws in native Fortran, by source category;
- numerical parity with original Python source evaluator on tested
  support and one-sided limits;
- original-EOS M2 rapid p1 and M3 slow p1 physical-mode parity,
  with retained phasewise junction conditions and node count;
- process isolation and measured speedup on the fixed workloads.

Not accepted:
- arbitrary microphysical EOS with different representations;
- several Maxwell discontinuities simultaneously;
- every possible upper-knot near-luminal numerical condition;
- source-exact sharp-mode root iteration fully native Fortran;
- complete global f/p1 branch continuation and full-GR damping;
- frozen composition, bulk-buoyancy g-mode claims;
- global EOS-wide 140x speedup or guaranteed memory savings;
- whole-Netlib F77 bounds safety.

**Disposition: PASS_SCOPED_ORIGINAL_EOS_FORTRAN_CS2**, limited to
the two source-locked effective barotropes and their tested modes.
The next constructive work is native Fortran phasewise **Brent
polishing with the same original-source c_s²** and full native
strong-knot scheduler, preserving slow/rapid Maxwell semantics.
After that, stronger convergence and independent branch-identity
checks are scientifically meaningful.

References: Tonetto & Lugones (2020),
https://arxiv.org/html/2003.01259v2 ;
thermodynamically consistent crust-core interpolation,
https://arxiv.org/html/2502.02373 ;
original Netlib DVODE source,
https://netlib.org/ode/vode.f .
