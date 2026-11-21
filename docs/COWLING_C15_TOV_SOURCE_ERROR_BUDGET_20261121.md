# NS_Cowling C15: Original-source TOV accuracy and Fortran two-phase batching

Date: 2026-11-21 KST.
Parent: 091c2f6f22d032379324aa900c910c6278a68437.
Research branch: research/cowling-c15-native-tov-source-sweep-20261121.
Scope: a bounded Fortran TOV numerical implementation and TOV-to-Cowling sensitivity study, not a new full ten-stage campaign.
Native research subagents: NOT_RUN (no native delegation tool exists); numerical processes are not research subagents.

## 1. Physical question and provenance

C14 held the TOV background fixed and found nonmonotone Cowling numerical-control sensitivity.
C15 studies the previously unseparated background-TOV contribution, keeping the original cold EOS, central density and Cowling solver settings fixed per star.

Frozen original source archive:
data/eos_missing_data_recalc_20261119/artifacts/EOS_ALL_REQUESTED_SOURCE_LOCKED_NESTED_20261119.zip
SHA256 = 299de43bb4f0a87ebfc4f1139ab3c8b51c344937ec3637cdff644e9ca8ab2c04.

Source contracts:
M2 Maxwell SHA256 59b7cb4671347c33c25133e6497e8be346b4e5b4f871eece9bcaf59155fd8765.
M3 CSS SHA256 f23fb7020c9de1b65c99ac37690840f7a4d6611f3372409b015d8c3bd7e32f3b.

M2 n_c = 0.30385128926615257 fm^-3, M ~ 1.4 Msun, R ~ 14.42041962 km,
rapid-conversion p1, bracket 5110..5170 Hz.
M3 n_c = 0.192510767556007 fm^-3, M ~ 0.7 Msun, R ~ 14.59794179 km,
slow-conversion p1, bracket 3710..3760 Hz.

Metric signature (-,+,+,+). Geometrized r,m are in km.
Source P and epsilon from MeV fm^-3 are converted to km^-2; retain
c = 299792.458 km/s when converting f_Hz to omega_geom = 2 pi f/c.
Closure is Gamma1=Gamma_eq. No full-GR polar metric perturbations, frozen-composition buoyancy, damping or generic g-mode claim.

## 2. Equations and Maxwell identity

Dimensionless enthalpy h = ln(mu/mu_surface), coordinate s = sqrt(h_start-h).
The TOV equations in the conventions of the actual project are:

    dr/dh = -r(r-2m)/(m+4*pi*r^3*P)
    dm/dh = 4*pi*r^2*epsilon*(dr/dh)
    dr/ds = 2*s*r*(r-2m)/(m+4*pi*r^3*P)
    dm/ds = 4*pi*r^2*epsilon*(dr/ds)

Regular near-center initialization:

    h_start = h_c-2*pi*(epsilon_c/3+P_c)*r0^2
    r(0)=r0, m(0)=4*pi*epsilon_c*r0^3/3

At the unique sharp Maxwell interface, [P]=[r]=[m]=0,
while [epsilon] is nonzero. The high-phase terminal (r,m) becomes
the low-phase initial value WITHOUT smoothing, reinterpolating across
the jump or retaining the prior phase's DVODE history.

## 3. Fortran implementation and old/new semantic parity

New source: fortran_eos/src/compose_eos_cowling_tov_two_phase_sweep.f90.
New C ABI: tov_native_two_phase_sweep.

Fortran receives ONE immutable P(h),epsilon(h) original-source sample
grid for each phase and multiple (RTOL,ATOL,max_step_s,method)
controls. It calculates the phase pressure and epsilon PCHIP
slopes via the C13 Fortran Fritsch-Butland routine ONCE,
then solves the high/low TOV equations in fresh DVODE histories.
Every control has a separate success/failure status and
RHS evaluation count; failure zeros its profile data and
does not invalidate other controls.

New optional thin frontend: solve_sharp_background_native_sweep,
in python_frontend/compose_eos/stellar_sharp_native.py.
Python reads original EOS source states once, marshals inputs,
and assembles the existing stellar background and final
physical-mode output. The TOV and Cowling differential
equations, PCHIP EOS derivatives, original-source c_s^2
and Brent roots execute in native Fortran.
Remaining Python numerical functionality: the established
_NativePhase uses SciPy PCHIP to query the completed TOV
r,m profile when constructing the later Cowling geometry.
Therefore no claim that the entire application is Fortran-only.

The original public solve_sharp_background_native and
all old default physical conventions are unchanged.
Original M2/M3 source-locked tests for Adams and BDF verified
the new and old background M/R/geometric states,
the slow/rapid p1 node counts and interface acceptance,
and isolated invalid controls (initial focused test:
10 passed in 10.71 s).

The final fresh C7-C15 relevant regression includes
smooth native TOV, the original Maxwell and source-exact
Cowling modes, C9-C14 and new C15 tests:
136 passed, 0 failed in 212.70 seconds.
Raw evidence:
data/cowling_c15_tov_convergence_20261121/INTEGRATED_REGRESSION.log.

## 4. Actual original-source TOV numerical controls

Fixed Cowling solver:
RTOL=1e-9, ATOL=1e-12, max_step_h=.002,
root_xtol_hz=1e-5; 4097 requested Cowling geometry samples.
The original source EOS, TOV central density and
Maxwell semantics are held fixed.

Three independent TOV axes, each tested with Adams and BDF:
- Requested source grid: 1025,2049,4097,8193, with RTOL=3e-10,
  ATOL=3e-11, max_step_h=.004.
- TOV RTOL: 3e-8,1e-9,3e-10,1e-10, with 2049 requested source samples.
- TOV maximum enthalpy step: .008,.004,.002,.001,
  with 2049 requested samples, RTOL=3e-10.

Effective native maximum s-step is
min(smax/200,max_step_h/(2*smax)).
Duplicate baseline controls are evaluated only once:
20 unique TOV configurations per source star, 40 total.
Actual grids exceed requested sample counts because the
original EOS breakpoints and surface support are inserted.

Results from matched ORIGINAL-source backgrounds:
M2: 17 accepted p1, 3 failed TOV profiles, no
separate Cowling-root rejections.
M3: 16 accepted p1, no failed TOV profile,
4 rejected Cowling-root evaluations, all on the
2049-background BDF at RTOL=3e-10 (and its repeated
nominal max-step controls). The strict normalized
matching residual gate was not relaxed. These are
numerical rejection records, not physical nonexistence
proofs.

M3 TOV source-grid dependence of its slow p1 mode:

| TOV sample request | Adams [Hz] | BDF [Hz] |
|---:|---:|---:|
| 1025 | 3736.212646486 | 3736.212509935 |
| 2049 | 3736.210595628 | REJECTED |
| 4097 | 3736.209216741 | 3736.209165186 |
| 8193 | 3736.208845083 | 3736.208944173 |

Accepted frequency spans on the TOV source-grid axis:
M2 rapid p1, Adams 0.000245508 Hz; BDF 0.000327961 Hz.
M3 slow p1, Adams 0.003801403 Hz; BDF 0.003565762 Hz
(the rejected 2049 BDF control is excluded).
These are finite-control ranges, NOT error bars,
certified convergence orders or absolute eigenfrequency
uncertainty bounds. Distinct source-grid requests
do not correspond to uniform grid-halving after the
source EOS knots are inserted.

For exact control statuses, star masses and radii,
interface geometry, join residuals and native RHS counts:
data/cowling_c15_tov_convergence_20261121/TOV_COWLING_SENSITIVITY_SUMMARY.json
and the original per-control M2/M3 JSON science files.
No invented frequency is assigned to rejected cases.

## 5. Native Fortran executable and memory optimization

Original source-locked M2/M3 two-phase P,E samples were
exported once as compressed deterministic inputs:
C15_ORIGINAL_TOV_SOURCE_M2.dat.gz and M3.dat.gz.
Source hashes are documented in ORIGINAL_TOV_INPUT_EXPORT.log.
The exporter does NOT integrate a TOV stellar ODE in Python.

Fortran-only regression:
fortran_eos/tests/cowling_c15_original_two_phase_tov.f90.

It tests original Maxwell pressure continuity and
epsilon jump, exact high/low r,m continuity,
real M/R reference, Adams and BDF results,
per-control invalid RTOL fail-closed, wrong
Maxwell transition rejection and repeat determinism.

Changed project Fortran was checked with -O0 -g -fcheck=all;
an -O2 rebuild emitted identical results. Frozen
Netlib ODEPACK vendor routines compiled -O2 -std=legacy;
the whole vendor library is NOT bounds-safety certified.

A verified implementation defect was repeated creation of
Fortran array temporaries for strided derivative views
into five-component PCHIP slope arrays. The code now
constructs four contiguous derivative channels ONCE
and releases the larger intermediate matrices before
repeated DVODE controls. The corrected results are
byte-identical to pre-optimization results.
Independent bounds-checked TOV standalone RSS:
M2 8908 KB -> 6588 KB; M3 7804 KB -> 6092 KB.
These are scoped single-process measurements, NOT
whole-pipeline memory savings.

## 6. Timing and nonreentrant process parallelism

Four matched original EOS TOV backgrounds,
median of three warm alternating comparisons,
excluding original source ZIP JSON loading but
including original EOS phase sampling and Fortran
TOV integrations:

| Source star | Four repeated old TOV calls | One four-control native TOV batch | Speedup |
|---|---:|---:|---:|
| M2 | 1.60197 s | 0.40368 s | 3.968x |
| M3 | 0.99656 s | 0.24977 s | 3.990x |

This gain comes from source sample/derivative
reuse and does not imply 4x full code speedup.

Two-star 40-control full pipeline, including
source loading, native TOV and admitted Fortran
Cowling p1 calculations:

| OS workers | Total wall time |
|---|---:|
| 1 | 37.319461 s |
| 2 | 20.499410 s |

Speedup approximately 1.8205x.
Canonical star output JSON files were byte-identical
for serial and parallel executions and compared equal
to pre-optimization outputs.
Peak full-front-end worker RSS around 112--124 MB.
DVODE global COMMON workspaces are never shared
across active threads; all parallelism is OS-process
isolation. No GPU acceleration was claimed.

## 7. Validation gates, limitations and next scientific work

- One original Maxwell interface and two source-locked stars.
- Independent baseline parity for original TOV Adams/BDF,
  mass/radius, stellar geometry, physical p1 mode.
- Failed TOV/control rows and rejected root residuals
  preserved without conflating physics and implementation failure.
- Source TOV sample sensitivity is not the same as
  a true uniform-discretization Richardson convergence order.
- Original TOV center truncation r0 and EOS source
  microphysical representation remain fixed.
- The full mode frequency error budget does not have
  independently certified global bounds.

Status:
PASS_SCOPED_ORIGINAL_EOS_FORTRAN_TOV_SWEEP_WITH_GLOBAL_CONVERGENCE_OPEN.

Next bounded work: verify node-identified f/p1
branch continuation over several central densities
with source-lock provenance, correct mode labels and
per-density rejection statuses preserved.

Primary references:
Lindblom and Indik, enthalpy-based structure
https://arxiv.org/abs/1310.0803
Tonetto and Lugones, slow/rapid junctions
https://arxiv.org/abs/2003.01259
Netlib DVODE local-error documentation
https://netlib.org/ode/vode.f
