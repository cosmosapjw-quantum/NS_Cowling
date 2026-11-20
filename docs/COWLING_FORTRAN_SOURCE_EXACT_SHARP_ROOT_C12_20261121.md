# NS_Cowling C12: source-exact sharp Maxwell roots and native knot restarts

Date: 2026-11-21 KST. Scientific/program baseline:
`7b6bebee997e2a3cad069bddcd8580ea36f63839`.
Working branch:
`research/cowling-fortran-source-exact-sharp-root-c12-20261121`.
Active coding track: one bounded implementation-and-verification slice
following original-EOS Fortran c_s² in C11. This is not a claim
that all C1–C10 phases have been independently re-executed.
Native skill-bound subagents: **NOT_RUN**, because the current host
exposes no native delegate/spawn action. All numerical process-parallel
jobs below are OS processes, not subagents.

## 1. Protected physical question and conventions

For a source-locked sharp first-order Maxwell transition,
can the entire f/p1 eigenfrequency *root evaluation*, including
every original-source c_s² evaluation, ODE history restart at
strong intra-phase source derivative knots, and slow/rapid junction,
run inside Fortran while Python remains an input/output frontend?

Metric signature (-,+,+,+); physical radius and geometrized mass in
km, pressure and epsilon in MeV fm^-3, baryon density in fm^-3,
chemical potential in MeV. Dimensionless enthalpy:
`h=ln(mu/mu_surface)` and shoot coordinate
`s=sqrt(h_start-h)`. For frequency in Hz,
`omega_geom=2*pi*f/(c in km/s)`, with
`c=299792.458 km/s`.

The background is a cold equilibrium barotrope and the Cowling
closure remains **Gamma1=Gamma_eq**. No frozen composition, bulk
stratification g-mode, rotating star, full GR metric response or
GW damping claim is supported.

The two one-sided Maxwell phases share pressure and geometry but
not epsilon. For the scaled Cowling variables `(w,v)`, write
`H=epsilon+P`, `a=omega_geom^2 exp(lambda-2Phi)*r`,
`b=Phi_prime`, `T=H(a*v+b*w)`.
At the unique interface, the original implementation's conditions
are:
```text
slow conversion:  [w]=0, [T]=0
rapid conversion: [v]=0, [T]=0
```
Their algebraic identity was checked in Wolfram: all four
traction/kinematic residual expressions simplified to zero.
This preserves the original source semantics of Tonetto and Lugones,
*Phys. Rev. D* 101, 123029 (2020),
https://arxiv.org/abs/2003.01259 .

## 2. Source provenance (immutable)

Frozen full original effective EOS contract archive:
`data/eos_missing_data_recalc_20261119/artifacts/EOS_ALL_REQUESTED_SOURCE_LOCKED_NESTED_20261119.zip`,
SHA-256
`299de43bb4f0a87ebfc4f1139ab3c8b51c344937ec3637cdff644e9ca8ab2c04`.

The studied nested original files have source-contract SHA-256:
- `SPG_M2/maxwell/REPRESENTATION.json`:
  `59b7cb4671347c33c25133e6497e8be346b4e5b4f871eece9bcaf59155fd8765`.
- `SPG_M3/css/REPRESENTATION.json`:
  `f23fb7020c9de1b65c99ac37690840f7a4d6611f3372409b015d8c3bd7e32f3b`.

These describe the exact cold **effective source-contract piece laws**,
not a newly recovered microscopic EOS.

## 3. Concrete Fortran change and numerical derivation

A new `fortran_eos/src/compose_eos_cowling_sharp_source_root.f90`
exports the C-compatible
`cowling_native_sharp_source_brent`.
It evaluates all bracket residuals, moves the bracket through a
safeguarded Brent-Dekker interpolation/bisection cycle, and returns
the root, signed normalized determinant, iteration count, cumulative
Fortran DVODE RHS count, strong-knot count and typed status.

New `cowling_native_source_select_restarts` in
`fortran_eos/src/compose_eos_cowling_source_exact.f90`
constructs source knot locations directly from the original
surface/outer/bridge/upper chemical potential contracts.
At each strictly phase-interior knot h it evaluates
```text
dh=max(2e-11,1e-9*max(h,1e-5))
J(h)=|cs2(h+dh)-cs2(h-dh)|/
      max(|cs2(h+dh)|,|cs2(h-dh)|)
restart iff J(h) >= 0.01, with h +/- dh in the same phase.
```
Thus **Fortran**, not Python, chooses the strong EOS derivative
breaks, constructs the ordered ODE output waypoints and
one-based restart indices, and creates a fresh DVODE history for
each resulting segment. The separate physical Maxwell density
jump is not smoothed or treated as an intra-phase knot: the
original one-sided interface transfer remains separate.

Original deterministic source-knot selection reproduced:
- M2 Maxwell: **39 strong restart locations** in the two phases.
- M3 CSS: **40 strong restart locations** in the two phases.
All selected `s` positions matched the original Python selector
within `2e-12` in the checked support.

Each native determinant evaluation follows exactly:
1. high phase: regular center -> numerical matching radius;
2. low phase: surface DeltaP=0 -> physical Maxwell interface;
3. original native slow/rapid traction/displacement transfer;
4. high phase: interface -> numerical matching radius.

The normalized two-shooting determinant must approach zero.
A collapsed bracket with `abs(D)>1e-5` returns an explicit
root rejection instead of accepting a false sign flip.
An invalid/unbracketed request returns zero certified root values.
The new root aborts before installing any interpolant when an
unrelated Python source callback remains active.

`python_frontend/compose_eos/stellar_sharp_native.py` exposes
`root_backend="fortran_source_exact"` for
`coefficient_backend="native_source_cs2"`.
It only marshals previously validated physical data into native
Fortran and formats the final eigenfunction. The original
`python_brentq` path is kept as the backwards-compatible reference.
No source physics or default behavior was replaced.

## 4. Actual-source frequencies and accuracy

M2 Maxwell: `n_c=0.30385128926615257 fm^-3`,
`M=1.3999999857354968 Msun`, `R=14.420419623651872 km`,
rapid conversion, p1, bracket [5110,5170] Hz.

M3 CSS: `n_c=0.192510767556007 fm^-3`,
`M=0.7000000086373105 Msun`, `R=14.59794179069283 km`,
slow conversion, p1, bracket [3710,3760] Hz.

With source-exact c_s², 4097 requested one-sided geometry samples,
`rtol=1e-9`, `atol=1e-12`, `max_step_h=.002`,
`root_xtol_hz=1e-5`:

| Case | Original native DVODE + Python Brent (Hz) | C12 native source-exact Fortran Brent (Hz) | Difference |
|---|---:|---:|---:|
| M2 rapid p1 (Adams) | 5143.719608336112 | 5143.719608336112 | 0 within run's stored doubles |
| M3 slow p1 (Adams) | 3736.210698739710 | 3736.210698739710 | 0 within run's stored doubles |

The native root profile still satisfies the original 1/1 phasewise
W,V node classification and `relative_join_mismatch<1e-5`.
The stronger independent higher-precision original-source DOP853
reference frequencies are `5143.719721563306` (M2) and
`3736.2088951581254` Hz (M3); the native errors versus these
separate references are about `-0.00011323` and
`+0.00180358` Hz, respectively.
**A root-finder tolerance is not an absolute astrophysical-frequency
error certificate.** Background interpolation, source geometry,
ODE tolerances and branch identification remain independent
error-budget contributors.

A separate pure Fortran executable tested both Adams and BDF:
- M2 rapid p1: **5143.719608336** Hz (Adams),
  **5143.719703958** Hz (BDF).
- M3 slow p1: **3736.210698740** Hz (Adams),
  **3736.210573488** Hz (BDF).

For the **same frozen source and numerical settings**, these
four native roots differ from the corresponding original
native-source + SciPy-Brent reference by at most about 0.00013 Hz.
This figure is **not** the difference from independent DOP853.
The independent source data are exported in
`ORIGINAL_FORTRAN_ROOT_{M2,M3}.dat.gz`; neither test
executable imports Python.

## 5. Implementation tests and failure classification

TDD baseline: the new Fortran source-knot selector ABI was absent,
so the new test failed before implementation. After implementation,
M2 and M3 selector parity tests passed.

Focused C12 tests include:
- M2/M3 real original-source p1 frequency parity vs
  the existing SciPy outer-root path and external references;
- no Python `matching_residual` calls in native root polishing;
- preserved independently reconstructed physical eigenfunction
  node and interface conditions;
- exact native strong-knot counts and positions;
- unsupported coefficient backend and active external source
  callback rejection, with all certified outputs zeroed;
- degenerate and independently preverified nonbracketing
  intervals rejected.

One initial negative test mistakenly assumed 5600–5700 Hz
contained no root for M3. That physical premise was incorrect,
not a root-solver false positive. The test now first verifies the
signs of the independent reference determinant at its chosen
endpoints, and the complete focused test passed. The raw failure
transcript was preserved with its type `TEST_PREMISE`.

Fresh focused suite: **12 passed, 0 failed**.
Fresh integrated source-locked tests across C7–C12 and all
currently used sharp/smooth, hybrid, Maxwell, TOV and
original-source paths: **107 passed, 0 failed in 175.24 s**.
Logs: `FOCUSED_POST_CALLBACK.log` and
`INTEGRATED_REGRESSION_FINAL.log`.

Standalone pure-Fortran code is in
`fortran_eos/tests/cowling_source_exact_sharp_root_c12.f90`.
The source-locked frontend-only input exporter is
`scripts/export_cowling_source_exact_root_c12_20261121.py`.
A reproducible bounds-checking build helper is
`scripts/build_cowling_source_exact_root_c12_20261121.sh`.
The independently compiled original Fortran modules were checked
with `-O0 -g -fcheck=all`. Untouched Netlib Fortran-77 ODEPACK
was linked using `-O2 -std=legacy`: this does not establish
whole-vendor memory safety.
Rebuilding the modified Fortran modules with `-O2` produced
identical printed scientific results for these four modes.

## 6. Optimization, memory and process-isolated parallelism

Seven alternating warm full-mode timings on the same original
background and PCHIP **geometry** (Python frontend time excluded):

| Case | SciPy outer root median | Native Fortran outer root median | Measured gain |
|---|---:|---:|---:|
| M2 rapid p1 | 0.0311073 s | 0.0252700 s | 1.231x |
| M3 slow p1 | 0.0387935 s | 0.0305119 s | 1.271x |

The new Fortran source-root adds this gain to C11's removal
of Python EOS RHS callbacks. The figures are **not**
global 140x/147x plus 23%/27% compound speedup claims.
The source-locked star background preparation itself took
0.406 s (M2) and 0.212 s (M3), and one-sided coefficient
preparation took 1.288 s and 0.730 s, respectively, in the
benchmark run; these are excluded from the warm medians.

The two physical inputs are compressed with deterministic
`gzip -n -9`; each full input contains original source EOS
parameters plus separate high/low phase PCHIP geometry.
The compressed payloads are about **1.5 MB** (M2) and
**1.0 MB** (M3), down from approximately 4.0 and 2.6 MB
respectively. Raw and gzip content SHA-256 are recorded in
`UNCOMPRESSED_INPUTS.sha256` and
`COMPRESSED_INPUTS.sha256`.

The standalone scientific executable's observed maximum RSS
was roughly **7–9 MB per OS process** for original-source
M2/M3 root finding. That excludes Python EOS/TOV/PCHIP
materialization and must **not** be advertised as an
end-to-end memory reduction.

Bounded OS-process throughput test: **20 repeats** of Adams
and BDF root polishing per original star, **80 native root
solutions total**. The independent source and numerical output
were byte-identical across serial and two-process runs:

| Original-source workers | Total wall time |
|---|---:|
| 1 process | 3.773121 s |
| 2 processes | 2.103439 s |

Speedup **1.793x** on this precise workload.
Canonical per-star output manifest SHA-256:
`bcb02126b1a61541bfc453422fe8331a85fac13d9cbe41f1759d12a2034d5862`.
Every repeated native root was byte-identical to the first
within the same process. The strong-knot selector and all
DVODE COMMON/SAVEd arrays are process-local, not
thread-safe OpenMP shared computations.

This is not a production throughput benchmark across a full EOS
ensemble and not a GPU speedup.

## 7. Reproduction

After checkout at this report's eventual final commit
with the same unchanged source EOS zip:

```bash
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1
export PYTHONPATH=python_frontend
python3 -m pytest -q \
  python_frontend/tests/test_cowling_source_exact_native_root_c12_20261121.py \
  --run-slow
python3 scripts/bench_cowling_source_exact_root_c12_20261121.py
# Regenerate source-locked raw Fortran driver inputs:
python3 scripts/export_cowling_source_exact_root_c12_20261121.py
scripts/build_cowling_source_exact_root_c12_20261121.sh \
  /tmp/cowling_c12_replay/native_sharp_root
/tmp/cowling_c12_replay/native_sharp_root \
  data/cowling_fortran_source_exact_root_c12_20261121/ORIGINAL_FORTRAN_ROOT_M2.dat 1
/tmp/cowling_c12_replay/native_sharp_root \
  data/cowling_fortran_source_exact_root_c12_20261121/ORIGINAL_FORTRAN_ROOT_M3.dat 1
```
The checked-in versioned inputs are `.dat.gz`; use
`gzip -dc < file.dat.gz > /tmp/file.dat` before calling the
pure Fortran driver. The exporter recreates fresh `.dat`
copies only inside its declared input output directory.

## 8. Verdict / bounded next action

**PASS_SCOPED_ORIGINAL_EOS_FORTRAN_SHARP_ROOT** for
the listed original-source M2 rapid and M3 slow p1 cases.
No branch/PR merge, global convergence, universal multi-jump EOS
extension, full frozen-composition dynamics or full-GR
polar-mode certification is implied.

What remains in Python: original source EOS JSON input
materialization and validation, TOV frontend dispatch,
phase geometry PCHIP preparation and final mode-profile
marshaling and node counting. C12 moves the *root iterations
and entire selected internal restart scheduler* to Fortran.

Next constructive stage: Fortran-native background geometry
PCHIP coefficient setup from immutable stellar samples, then
cross-sample and cross-EOS convergence and physical branch
continuation on a greater number of original source stars.
Do not infer a universal f/p1 classification from roots
and two pointwise models.
