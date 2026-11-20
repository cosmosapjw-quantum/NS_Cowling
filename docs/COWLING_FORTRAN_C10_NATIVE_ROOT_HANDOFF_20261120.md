# NS_Cowling C10: Native Fortran Brent root-polishing and bounded handoff

Date: 2026-11-20 KST.
Parent / authoritative C9 branch commit:
`0952a5c8ca34c5d94136a2c89c705e000c21ef3f`.
Current branch: `research/cowling-c10-native-brent-20261120`.
PhysMath stage: C10 release/continuation judgment, with one explicitly
requested final Fortran implementation slice. Native subagent:
**NOT_RUN**, no native spawn/delegation capability exposed. All
source mutation was restricted to this isolated linked worktree.

## Protected definitions

Metric signature (-,+,+,+); geometric radius and mass in km;
frequency in Hz with `omega_geom = 2*pi*f / (299792.458 km/s)`.
The Cowling equations, one-sided Maxwell boundary conditions,
`Gamma1=Gamma_eq` barotropic closure and EOS input bytes are unchanged.
Frozen-composition `Gamma1`, stratification g-modes, full GR polar metric
perturbations, and sharp-interface native bulk/root polishing remain
**not implemented by this smooth path**.

The scalar spectral determinant is

```text
D(f) = (u_w(f) * v_v(f) - u_v(f) * v_w(f)) /
       (norm_2(u(f)) * norm_2(v(f))),
```

where u, v are the two independent center/surface DVODE shooting vectors
at the matching radius. A sign change is a *necessary bracketing signal*
for an isolated simple zero, not an automatically sufficient physical
f/p-mode certification.

## Concrete implementation and change to the default

- Fortran `fortran_eos/src/compose_eos_cowling_dvode.f90` now implements
  `cowling_native_root` and a private native determinant evaluation helper.
  The method uses bracket-preserving Brent-Dekker inverse interpolation
  with bisection fallback. It initializes a fresh Netlib DVODE state at
  every trial frequency, never shares state between simultaneous threads,
  and returns frequency, root residual, native iterations, RHS calls,
  or a typed nonzero status.
- C ABI call requires a strictly positive, ordered bracket, finite
  parameters, positive numerical tolerances, a physically bounded
  matching coordinate and no active sharp-source callback.
  Same-sign endpoint residuals return `status=3`; active sharp callback
  returns `status=2`. On failure the certified root and residual
  outputs are zero. Iteration exhaustion returns `status=32`.
- New **false-PASS protection**: even after the Brent bracket collapses
  to the requested `root_xtol_hz`, accept a physical root candidate only
  when `abs(D(f)) <= 1e-5`. A numerically discontinuous/underresolved
  sign flip that has `abs(D)>1e-5` returns `status=33`. This gate
  is not a physical mode-identity theorem.
- The native Python wrapper
  `python_frontend/compose_eos/stellar_cowling_native.py`
  now defaults to `solve_mode(..., root_backend="fortran")`.
  `root_backend="python_brentq"` is retained strictly as an explicit
  reference path. The Python frontend prepares data/metadata and does
  not execute the *new root-polishing iterations*.
- Native source-callback segregation, sharp Maxwell transfer and
  previous C8 coefficient lookup remain as before. No source EOS,
  Fortran TOV equation or vendor Netlib source was modified.

This is **not** a full Fortran-only scientific pipeline: source EOS
materialization, PCHIP preparation and mode-profile presentation still
have legacy Python code; `stellar_sharp_native.py` has a separate
sharp-interface mode root path.

## Scientific tests: meaningful success and failure

TDD RED: the newly required `root_backend="fortran"` failed in the
original project because no argument or Fortran symbol existed.
This was observed on the real 1.4 Msun source-locked star.

After the fix: **10 focused C10 tests passed**, exercising native
Adams/BDF f and p1, exact mode node counts, source-locked frequency
parity against explicit SciPy Brent, no Python root-callback dependence,
unbracketed input rejection, active sharp callback rejection and
invalid backend rejection.

Independent all-Fortran test
`fortran_eos/tests/cowling_c10_brent_root_test.f90`:

- A synthetic smooth but physically unconstrained coefficient
  profile has `D(1700 Hz)=-0.7283395`, `D(1800 Hz)=+0.7315829`.
  At a misleading sign flip near 1728.17893 Hz, the determinant
  remains about 0.729 in magnitude. The earlier implementation
  had falsely reported an apparent root (recorded in the raw forensic log).
  Final Fortran returns `status=33`, with no certified root.
- A uniform-density Schwarzschild interior (the acausal, analytic
  incompressible control, **not a realistic neutron-star EOS**)
  gives root frequency **67.473154193 Hz (Adams)** and
  **67.473154194 Hz (BDF)** with residuals `1.42e-11` and
  `1.94e-11`, respectively, consistent with the weak-gravity
  Kelvin estimate. This control validates arithmetic and boundary
  handling under its stated assumptions.

The modified project Fortran module and standalone test were compiled
with `gfortran -O0 -g -fcheck=all`, and Netlib ODEPACK was compiled
separately with `-O2 -std=legacy`. This distinction preserves the
known historical DGEFA Fortran-77 bounds-checker incompatibility;
**full-vendor bounds verification remains unresolved**.
Results: `data/cowling_fortran_c10_native_root_20261120/STANDALONE_FORTRAN_TEST.log`.

Post-change integrated tests for C7/C8/C9/C10, TOV, smooth and sharp
Cowling, EOS, Maxwell and hybrid: **61 passed, 0 failed,
0 skipped in 71.41 seconds**. Reproduce from
`data/cowling_fortran_c10_native_root_20261120/REGRESSION_RUN_FINAL.log`.

## Original EOS source and quantitative accuracy

Smooth source:
`data/core_support_20261117/primary_v1/REPRESENTATION.json`,
SHA-256
`591ac27aa99448512e291684c98f8273efc729d17c63be88e331c20e83bd2653`.

Original locked broader EOS ZIP remains unchanged:
`data/eos_missing_data_recalc_20261119/artifacts/EOS_ALL_REQUESTED_SOURCE_LOCKED_NESTED_20261119.zip`,
SHA-256
`299de43bb4f0a87ebfc4f1139ab3c8b51c344937ec3637cdff644e9ca8ab2c04`.

Native TOV for canonical smooth 1.4 Msun, 8,193 requested Cowling
coefficient samples, `rtol=1e-9`, `atol=1e-12`,
`root_xtol_hz=1e-6`; exact physical bracket inputs unchanged.

| Method, mode | Fortran frequency vs explicit SciPy on *same* native EOS, abs difference | Native normalized D |
|---|---:|---:|
| Adams f | 0.000187756 Hz | 2.802e-8 |
| Adams p1 | 0.001360341 Hz | 1.559e-7 |
| BDF f | 0.000023714 Hz | -7.195e-9 |
| BDF p1 | 0.000858057 Hz | 5.408e-8 |

All root candidates in this set retain `W/V` node count 0/0
(f) or 1/1 (p1), and the center/surface join tolerance used in the
focused tests. These comparisons demonstrate *scope-matched mode parity*
within the existing <0.01 Hz coarse project acceptance envelope,
not guaranteed eigenfrequency precision 1e-6 Hz.

A `rtol` sweep (1e-9, 1e-10, 1e-11, reducing `atol` in concert)
showed the Adams p1 discrepancy between the two root algorithms
falling to ~3.47e-6 Hz at 1e-10. But the f-mode difference between
the native EOS and the separate exact-contract DOP853 reference
remained ~0.0077 Hz at the two tightest tolerances.
The latter is **consistent with a background/coefficient discretization
floor**, not solely root-search algorithm error; this limited test does
not uniquely apportion the full error budget.

Numerical evidence: `TIMING.json` and `TOLERANCE_SWEEP.json`.

## Optimization and process parallelism

Warm per-mode wall-time medians from seven *fully profiled*
mode calculations, with Python `brentq` solely as comparison:

| Method/mode | SciPy reference (s) | Fortran Brent (s) | observed speedup |
|---|---:|---:|---:|
| Adams f | 0.01913 | 0.02788 | 0.69x |
| Adams p1 | 0.04203 | 0.04982 | 0.84x |
| BDF f | 0.04441 | 0.03419 | 1.30x |
| BDF p1 | 0.05732 | 0.07084 | 0.81x |

Native root migration is **not uniformly faster**. No blanket
speed improvement is claimed. Different Brent sampling patterns
change the number of expensive DVODE trial shots. The main purpose
of this change is to transfer root numerics to Fortran while
preserving independent reference and explicit error classification.

A separate real-source 4-star workload evaluates 121 frequencies per
star, refines detected brackets in native Fortran and produces
120--240-point native eigenfunction profiles. Central baryon density
inputs are 0.35, 0.45, 0.54, 0.70 fm^-3.
Nine numerical root candidates were computed with zero rejected
candidates in this small set, retaining explicit non-certification
of global mode identity. Process-isolated worker performance:

| Independent processes | Four-star wall time |
|---|---:|
| 1 | 4.706466 s |
| 2 | 2.387437 s |
| 4 | 1.476012 s |

The 1/2/4-process per-star output JSON and checksums are
**byte-identical**; total time speedup is ~3.19x relative to
the serial implementation on this exact workload.
This is process parallelism, **not unsafe OpenMP into DVODE COMMON**.

One observed end-to-end single-star root scan/polish run at
`n_c=0.54 fm^-3` gave maximum RSS 88,028 KB for the Fortran-root
path versus 87,704 KB for the SciPy-root reference; wall times
1.02 s vs 1.04 s. Neither a memory-saving nor cold-start
performance improvement is established. Four independent workers
duplicate per-process memory and frontend setup.

Scripts:
`scripts/run_cowling_c10_native_root_case.py` (thin frontend),
`scripts/bench_cowling_c10_root_parallel.sh` (bounded, independent
worker processes),
`scripts/benchmark_cowling_c10_native_roots.py`
(timed parity benchmark, no new Python numerical solver).
Scientific receipts and raw trials:
`data/cowling_fortran_c10_native_root_20261120/`.

## C10 disposition and exact next action

**SCOPED_RELEASE_CANDIDATE / C10 CLOSED for the smooth Fortran
scan/root milestone**, provided the final commit and remote/local
identity are verified. Do **not** represent this as full Fortran-only
stellar physics, sharp Maxwell mode-root certification, full EOS
matrix validation, or a merged release.

Next bounded scientific implementation target (a new authorized
coding cycle, not automatically started here):
phasewise sharp-interface native root evaluation and original
source EOS/PCHIP coefficient construction moved from Python to
Fortran, then compare mass/mode trajectories with explicit
EOS identity and phase conversion and quantify interpolation error
and mode branch continuation. Preserve slow/rapid interface
conditions, physical source histories and fixed thresholds.

Primary root-finding algorithm source:
[Netlib double-precision zeroin.f](https://netlib.org/fmm/zeroin.f).
Primary physical context:
[Counsell et al., relativistic Cowling g-modes](https://arxiv.org/abs/2409.20178);
[Kunjipurayil et al., f/p mode EOS study](https://scispace.com/papers/impact-of-the-equation-of-state-on-f-and-p-mode-oscillations-15g1elv5).
