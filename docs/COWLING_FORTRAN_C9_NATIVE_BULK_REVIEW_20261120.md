# NS_Cowling C9: independent numerical review and Fortran-native bulk smooth-star scan

Date: 2026-11-20 KST.
Parent scientific / implementation revision:
`3e90be5645b10007c82285d274bb7c0c0dfa6674`.
Working branch:
`research/cowling-fortran-bulk-scan-c9-20261120`.
PhysMath stage: **C9 independent, bounded code/numerics review**.
Native skill-bound subagents: `NOT_RUN` (no native delegate/spawn tool exposed);
actual work and check were serially managed, with process-isolated numerical jobs.

## 1. Protected numerical and physics invariants

- Metric signature (-,+,+,+); background r and m in km, source n_B in fm^-3.
- c=299792.458 km/s, omega_geom=2*pi*frequency_hz/c.
- One-dimensional equilibrium barotropic Gamma1=Gamma_eq. The source barotrope
  does not determine a frozen-composition adiabatic index. No composition g-mode
  or full-GR oscillation/damping claim.
- Original EOS bytes, native TOV, Maxwell slow/rapid interface transfer,
  sharp source callbacks, ODE tolerances and root-acceptance gates untouched.
- Original Netlib DVODE has shared COMMON, and native module tables are SAVEd:
  no concurrent *threads* may enter one native solver image.
- Native bulk scan is strictly **smooth-star only**. A true sharp interface
  still requires the separate original phasewise solver.

## 2. Current-task findings, reproducible repairs, and tests

**C9-F01, fixed:** `NativeCowlingDVODE.scan_brackets` executed every
frequency as a separate Python-level `matching_residual -> _shoot` call.
That maintained correct numerical physics but forced per-frequency Python
dispatch/lock and failed to retain a native batch boundary. A test that
prohibited Python scalar shots reproduced the defect on the original
smooth-source 1.4 Msun star (Adams and BDF).
Repair: new `cowling_native_scan(...)` in
`fortran_eos/src/compose_eos_cowling_dvode.f90`, with a C-compatible
binding and a thin ctypes frontend. Full frequency validation, DVODE
shots and sign-changing brackets are now evaluated inside Fortran. Each
frequency restarts DVODE (ISTATE=1); no unresolved shared worker state.
On any invalid/failed input, status is nonzero and all output arrays
and cumulative RHS counts are set to zero (no partially admitted scan).

**C9-F02, fixed:** the public native scanner did not explicitly reject
an active one-sided sharp-EOS source callback. An injected callback
regression showed the old routine invoked the callback before rejecting
the shot (old status 11, not the required immediate disjoint-path status).
Repair: explicitly return status=2 *before* scanning when
`associated(active_source)` is true. A regression verifies callback
is not invoked and no results are exposed. Existing sharp solvers are
untouched.

These are scoped repairs to an existing numerical path, not a new
general treatment of nonbarotropic fluid oscillations.

## 3. Source-grounded numerical equivalence and performance

Original real smooth EOS:
`data/core_support_20261117/primary_v1/REPRESENTATION.json`,
SHA-256
`591ac27aa99448512e291684c98f8273efc729d17c63be88e331c20e83bd2653`.
One canonical native-TOV star: M=1.4000000333 Msun,
R=11.7520783024 km, native PCHIP coefficient support 2049 requested
samples, rtol=1e-9, atol=1e-12, max_step_h=.002.
A 241-frequency 1800..7200 Hz spectrum scan was run four times for
the original Python loop calling Fortran and four times for one native
Fortran bulk entry. The maximum absolute native normalized shooting
determinant discrepancy was **0.0** on the tested grid.

| Variant | Four wall-time observations (s) | Median (s) |
|---|---|---:|
| Previous Python-dispatched scalar native shots | 0.55143, 0.53436, 0.55737, 0.53735 | 0.544389 |
| New Fortran bulk native shots | 0.52981, 0.52879, 0.52520, 0.52495 | 0.526995 |

Warm scan-only speedup ~1.033x; **no claim of large single-star speedup**.
Background construction takes a separate observed 0.346321 s.
These measurements preceded the extra callback safety check; that check
does not touch the normal active-source-null data path.
Exact timing evidence: `data/cowling_fortran_c9_bulk_scan_20261120/TIMING.json`.

A separate one-shot end-to-end frontend observation at n_c=0.54 fm^-3
recorded maximum RSS 88,332 KB for scalar native shooting and
87,952 KB for bulk native scanning; both process wall times rounded
to 1.24 s. This is consistent with **no demonstrated memory or
cold-start time improvement**; it is neither an uncertainty estimate
nor a general lightweight-computing claim. The process-parallel method
necessarily duplicates Python/Fortran process memory across workers.

### Independent actual-source star batch (Python frontend, Fortran integration)

Four original smooth EOS TOV star input central baryon densities:
0.35, 0.45, 0.54 and 0.70 fm^-3.
Observed gravitational masses:
0.7686503, 1.1295525, 1.4062437 and 1.7517946 Msun.
Each independent worker uses a separate OS process and calls the native
Fortran TOV and native Fortran bulk Cowling scanner for 241 points.
The legacy Python EOS/PCHIP **setup** remains on the frontend.

| Worker processes | Measured batch wall time |
|---|---:|
| 1 | 5.098019 s |
| 2 | 2.616172 s |
| 4 | 1.404418 s |

All four exact per-star JSON science records are byte-identical
between the 1/2/4-worker runs. Canonical comparison:
`ca2d2709aeea3308c79b65c6291d68c20bebdd7c332c955f8ec88a041c427395`
for `SCIENCE.sha256` in each measured run. The relative 4-worker batch
speedup vs serial is 5.098019/1.404418 ~= 3.63x, under this fixed
four-star source/scanning workload, not an arbitrary 390-EOS speedup.
The batch benchmark excludes mode root refinement and g-mode certification.

After the final sharp-callback fail-closed repair, the same four-star batch
was rerun at 241 frequencies each: **4.950767 s with 1 worker** and
**1.330808 s with 4 independent workers** (~3.72x for this single run).
Both final sets were byte-identical to the original 1/2/4-worker canonical
science receipt. This is the final-code performance readback; variance
and larger source cohorts were not assessed.

Producer scripts:
- `scripts/run_cowling_native_bulk_scan_c9.py`: frontend-only input
  validation, native TOV/Cowling dispatch and JSON output.
- `scripts/bench_cowling_real_star_cases_c9.sh`: independent worker
  process scheduling, complete-case checking and canonical SHA verification.
- `data/cowling_fortran_c9_bulk_scan_20261120/DENSITIES.txt`:
  explicit ordered source-locked input sample.

## 4. Standalone independent Fortran test and vendor limitation

`fortran_eos/tests/cowling_c9_bulk_scan_test.f90` checks the native API
directly, with no Python numerical runtime. Adams and BDF bulk results
match separately restarted scalar native shots, including exact NFEV sums
and bracket flags; invalid frequency fails closed.

The initial whole-program `-fcheck=all` attempt hit
`fortran_eos/vendor/odepack/dgefa.f:87` when BDF invokes a legacy
Fortran-77 matrix interface whose formal second dimension is declared as 1.
This toolchain result is **not** evidence that our new native scan itself
is out-of-bounds, nor does it certify vendor memory safety. Do not mutate
the frozen vendor source merely to suppress this observation.

The corrected scoped check compiles the **modified project Fortran
module and its standalone test** with
`gfortran -O0 -g -fcheck=all -ffree-line-length-none`,
then compiles the original vendored ODEPACK linear algebra in
`-O2 -std=legacy` compatibility mode. The linked scoped-bounds test
returns:
`FORTRAN_C9_BULK_SCAN_PASS_ADAMS_BDF_REJECT_BAD_INPUT`.
Thus project-module bounds verification is PASS; full-vendor bounds
verification is **BLOCKED_LEGACY_ABI_CHECKER**, not PASS.

Reference: [original Netlib DVODE source](https://netlib.org/ode/vode.f)
documents /DVOD01/, /DVOD02/ COMMON state and per-problem restart rules.

## 5. Remaining obligations

- Python currently retains legacy EOS/PCHIP coefficient setup, plus
  `scipy.optimize.brentq` for final mode root extraction. This stage
  moves the **scan** computation into Fortran, but does not yet constitute
  a fully Fortran-only numerical application.
- The Fortran bulk API is limited to smooth backgrounds. Sharp Maxwell
  phasewise callback and junction transport remain supported only via
  the existing separate code path, deliberately unchanged.
- A sign-changing determinant brackets a root but does not by itself
  identify f/p1 or prove branch continuity near avoided crossings.
  Existing surface/junction/node gates still apply to actual mode solutions.
- One benchmark run and one original source family cannot justify
  globally predictive speedups, scaling or error bounds.
- Full Netlib BDF bounds instrumentation remains an isolated
  vendor-compatibility item, not a physics failure.

## 6. Fresh post-repair verification

The final regression replay, including TOV, smooth Adams/BDF,
native Maxwell interface tests, sharp hybrid cs² modes, C7 and C8
evidence tests, and the newly added fail-closed callback case,
finished **51 passed, 0 failed, 0 skipped in 70.70 seconds**.
Retained stdout:
`data/cowling_fortran_c9_bulk_scan_20261120/REGRESSION_RUN_FINAL.log`.
The standalone native Fortran Adams/BDF and invalid-input test also
passed after the final modified source was recompiled with the
native-module-only bounds instrumentation. Detailed receipts and
SHA-256 inventories live under
`data/cowling_fortran_c9_bulk_scan_20261120/`.

C9 disposition: **LOOKS_GOOD_SCOPED_WITH_LIMITATIONS**, backed by actual
post-change checks. This is not a newly performed independent native
subagent review; current runtime exposed no native delegation capability.
The natural next stage is C10 final decision/handoff; subsequent targeted
implementation may port smooth-mode root polishing and full EOS
interpolation to Fortran, preserving exact physics and thresholds.
