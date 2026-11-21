# NS_Cowling C16: Original-source f/p1 local central-density branch continuation

**Date:** 2026-11-21 KST.
**Parent baseline:** `8c0f70eed5595b10933763b3917ad4e045b6cf94`.
**Isolated development branch:**
`research/cowling-c16-source-mode-continuation-20261121`.
**Selected work:** one bounded constructive coding / scientific
verification slice following the C15 TOV source-error budget.
**Native research subagents:** `NOT_RUN`: no host native
delegate/spawn capability. Independent Unix processes for numerical
parallelism are not research subagents.

## 1. Protected scientific question and locked sources

Can f and p1 modes be separately identified and followed
across nearby central densities of TWO already source-locked
Maxwell hybrid stars while retaining exact original-source
c_s^2, slow/rapid interface conversion and independent
phasewise node/junction acceptance?

Original physical EOS archive (immutable):

`data/eos_missing_data_recalc_20261119/artifacts/EOS_ALL_REQUESTED_SOURCE_LOCKED_NESTED_20261119.zip`

SHA-256:
`299de43bb4f0a87ebfc4f1139ab3c8b51c344937ec3637cdff644e9ca8ab2c04`.

Nested original EOS contracts:
- M2 `SPG_M2/maxwell/REPRESENTATION.json`:
  `59b7cb4671347c33c25133e6497e8be346b4e5b4f871eece9bcaf59155fd8765`.
- M3 `SPG_M3/css/REPRESENTATION.json`:
  `f23fb7020c9de1b65c99ac37690840f7a4d6611f3372409b015d8c3bd7e32f3b`.

Original center densities:
M2 `n_anchor=0.30385128926615257 fm^-3`;
M3 `n_anchor=0.192510767556007 fm^-3`.
The source contract independently accepts both anchor +/-2%
densities and keeps the original sharp Maxwell transition.

Physical metric signature (-,+,+,+).
Geometrized mass/radius are in km;
the measured frequency is in Hz, with
`omega_geom=2*pi*f/(c=299792.458 km/s)`.
EOS matter is the original equilibrium cold barotrope:
`Gamma1=Gamma_eq`, not a physical frozen-composition
adiabatic index. No bulk g-mode or full-GR polar
metric perturbation is assumed.

At a sharp Maxwell interface, pressure, mass and
background radius are continuous while energy density
is discontinuous. Original C12 Cowling matching transport:
```text
H=epsilon+P, a=omega_geom^2 exp(lambda-2Phi)*r, b=Phi_prime
T=H*(a*v+b*w)
slow conversion: [w]=0, [T]=0
rapid conversion: [v]=0, [T]=0
```
Root admission still requires both the original
`abs(normalized_signed_matching_determinant)<1e-5`
and independently reconstructed
`relative_join_mismatch<1e-5`.
These gates were **not relaxed**.

## 2. Constructive Fortran upgrade and TDD

Previously, Python implemented the frequency scan as
one `matching_residual(f)` call per proposed frequency,
and counted physical eigenfunction nodes with NumPy.
A sign reversal alone was never a proof of a physical mode.

New source files:
- `fortran_eos/src/compose_eos_cowling_source_scan.f90`
  exports `cowling_native_sharp_source_scan`.
  It validates a monotone positive frequency grid,
  original star data, the Maxwell interface and
  source contract, computes the source derivative
  restart list **once** per physical phase, then
  evaluates the entire 97-point shooting scan in
  Fortran. Each frequency gets a typed status,
  a normalized signed determinant and native
  RHS evaluation count. Failed shots report zero
  certified determinant values.
- `fortran_eos/src/compose_eos_cowling_mode_nodes.f90`
  exports `cowling_native_phase_nodes` and
  `cowling_native_select_mode_branch`.
  It counts phase-local radial W and V sign changes,
  discarding amplitudes smaller than 1e-9 of
  the corresponding phase/component peak. It
  never counts across the physical Maxwell jump.
  A candidate is linked to a prior branch only
  when **exactly one** accepted candidate has
  the correct W/V nodal numbers and
  `abs(f_new-f_old)/f_old<0.12`.
  Ambiguity is an explicit unsuccessful status,
  not an invitation to silently pick the
  closest of two otherwise eligible candidates.
- `fortran_eos/src/compose_eos_cowling_mode_overlap.f90`
  exports `cowling_native_mode_shape_overlap`.
  It compares radial W profiles separately in
  high and low Maxwell phases. Within each
  phase, map the stellar **areal radius** data
  (not radial proper distance) into a
  dimensionless local `u=(r-r_phase_min)/
  (r_phase_max-r_phase_min)` and interpolate
  linearly on a fixed 64-cell grid. The signed
  cross integral and component norms give the
  absolute correlation, defined over both
  phases and separately for each phase.

**Important:** This W-only overlap diagnostic is
*not* a conserved/covariant pulsation-mode inner
product, does not identify mode restoring forces,
and cannot certify the absence of an avoided
crossing. Sign ambiguity of an eigenfunction is
removed by the absolute correlation.
It deliberately **never interpolates across
the Maxwell discontinuity**.

The existing
`fortran_eos/src/compose_eos_cowling_sharp_source_root.f90`
exposes its already validated `exact_onejump_det`
internally to the new native scanner without
modifying any determinant or source-EOS formula.

The Python `stellar_cowling_native.py` frontend
registers the exact new Fortran ABIs.
`stellar_sharp_native.NativeSharpCowlingDVODE`
adds `scan_source_brackets_native`, a thin
marshaling interface. For `coefficient_backend=
"native_source_cs2"`, final physical W/V
node counting now uses Fortran. The SciPy root
backend and all legacy source/callback and
PCHIP paths retain their previous explicit opt-in
behavior; no silent source physics change.

Other Python code in
`scripts/run_c16_source_density.py`
loads original EOS input, creates numerical
work packets, performs source/output serialization,
and dispatches the Fortran scan/Brent/overlap
routines. It holds no alternative EOS physical
law, determinant integrator or PCHIP slope code.
The original C15 Python TOV profile frontend
and final mode profile construction remain
a genuine Python dependency. Thus this is
**not** yet a fully standalone Fortran scientific
application for arbitrary EOS input files.

Negative TDD tests were constructed and observed
before implementation: missing native phase
node ABI, native branch-selector ABI and native
shape-overlap ABI. The raw failure transcripts
are versioned, not mistaken for physics failure.

## 3. Measured local density continuation

For each source model, five actual central
densities were computed at multipliers
`0.98,0.99,1.0,1.01,1.02`.

Fixed original TOV settings:
`rtol=3e-10,atol=3e-11,max_step_h=.004`,
requested TOV background samples 2049.
Fixed Cowling settings:
`rtol=1e-9,atol=1e-12,max_step_h=.002,
root_xtol_hz=1e-5`, requested geometry samples
4097, 360 profile points.
Every shooting RHS evaluates source-exact
one-sided c_s^2 via Fortran and the unchanged
strong-knot restart schedule.

The 97-point frequency scan spans
700..3400 Hz and 3450..5500 Hz;
the latter is not an exhaustive spectrum.
In the same physical model, a sign-changing
interval must survive native Brent and Maxwell
mode-profile admission before receiving a
physical f or p1 name. The production scan
found 2 physical candidate roots per M2 star,
3 per M3 star; extra M3 candidates were
recorded and kept distinct from f/p1.

M2: Maxwell **rapid** conversion, original EOS
and five central density points:

| density factor | M/Msun | R [km] | f [Hz] (W,V=0,0) | p1 [Hz] (W,V=1,1) |
|---:|---:|---:|---:|---:|
| 0.98 | 1.363998497 | 14.430613932 | 1829.339786513 | 5112.358769864 |
| 0.99 | 1.382076430 | 14.425608741 | 1835.875796084 | 5128.166222921 |
| 1.00 | 1.399999986 | 14.420419624 | 1842.311659362 | 5143.719673347 |
| 1.01 | 1.417762893 | 14.415041390 | 1848.648051650 | 5159.021026071 |
| 1.02 | 1.435358283 | 14.409470426 | 1854.886023801 | 5174.072718421 |

M3: CSS **slow** conversion, original EOS
and five central density points:

| density factor | M/Msun | R [km] | f [Hz] (W,V=0,0) | p1 [Hz] (W,V=1,1) |
|---:|---:|---:|---:|---:|
| 0.98 | 0.676800649 | 14.614210220 | 1552.199738568 | 3623.515602409 |
| 0.99 | 0.688378167 | 14.605822193 | 1558.241239186 | 3679.892731125 |
| 1.00 | 0.700000009 | 14.597941791 | 1564.244894076 | 3736.210688357 |
| 1.01 | 0.711665227 | 14.590535998 | 1570.211481315 | 3792.454280319 |
| 1.02 | 0.723373250 | 14.583573041 | 1576.141681434 | 3848.607561523 |

The 20 listed f and p1 branch point
admissions retained the phasewise
zero/one node classification, existing
Maxwell junction gate and f<p1 mode order.
Every adjacent pair on either side of its
anchor has exactly one Fortran-approved
nodal successor within a 12% local frequency
step. No source_EOS physical state was
interpolated across the jump.

## 4. Cross-density eigenfunction shape check

The normalized radial W-shape correlations
were measured on the *actual reconstructed
physical f/p1 profiles*, mapping each Maxwell
phase separately into its own normalized
radius coordinate. Overall minimum absolute
correlations among adjacent density pairs:
- M2 f: at least 0.9999970;
- M2 p1: at least 0.9999650;
- M3 f: at least 0.9999652;
- M3 p1: at least 0.9999986.

The M3 p1 inner, high-density phase alone
has correlations approximately
0.9894..0.9929, while its low phase is
near unity. These values support *local*
nodal and shape continuity, not a global
avoided-crossing theorem or a conserved
mode-inner-product proof.

The first pilot additionally found up to
4-node pressure-mode candidates at the
M3 anchor when the frequency band extended
to 6850 Hz. Those modes are not included in
the reported f/p1 branch continuation and
are not relabeled as low-order modes.

## 5. Validation, performance and limitations

The fresh selected C7--C16 regression
produced **152 passed, 0 failed in 206.95 s**;
the additional C16-only suite produced
**16 passed in 6.45 s**.
All initial missing-ABI red test logs were
preserved. The added tests check the native
original-source scan against the prior
independent source-equivalent detector,
native f/p1 node counts, invalid frequency/
phase inputs, ambiguity failure, and
phasewise overlap sign-invariance.

The pure Fortran executable
`fortran_eos/tests/cowling_c16_original_source_f_p1.f90`
uses original M2/M3 C12 archived source
parameters and unmodified phase geometry.
It independently evaluates original-source
f and p1 signed determinants, both native
Brent roots, source-knot restarts, synthetic
phasewise nodes, branch selector and phase
shape overlap, with **NO Python runtime**.
The checked `-O0 -g -fcheck=all` project
build and optimized `-O2` project build
emitted byte-identical scientific output.
The frozen legacy Netlib F77 ODEPACK is
compiled `-O2 -std=legacy`: bounds
checking of those vendor routines is
not claimed.

Standalone peak RSS for original M2:
**8836 KB**, M3 **7216 KB** for
the exact-source scanner+f/p1 root test.
This is a scoped *standalone Fortran*
measurement. It excludes original
EOS ZIP loading, TOV background
materialization and Python frontend memory.
All 10-density main workflow
Fortran DVODE is process-isolated;
no shared-thread OpenMP state.

A matched-physics comparison with the
existing `scan_brackets` Python-per-frequency
dispatch was actually executed by
`scripts/bench_c16_scanner_speed.py`.
Five alternating warm repeats of 32 original
frequency points per source star gave:

| Case | Python-loop median | Native Fortran median | Gain |
|---|---:|---:|---:|
| M2 | 0.100059679 s | 0.079222767 s | 1.263x |
| M3 | 0.084644694 s | 0.061824602 s | 1.369x |

All checked normalized determinant doubles were
**exactly identical** between the original source
Fortran-DVODE Python-dispatch reference and the
new Fortran batch at those 32 points.
The new scanner avoids repeated Python
per-frequency dispatch and amortizes
original chemical-potential strong-knot
selection once per scan, not once per
frequency point. These timings **exclude**
source EOS input, TOV background and
PCHIP geometry preparation. They are not
an implied full-code universal speedup.

In the FULL two-star density-continuation
workload (10 independent original-source
stellar models, all scan/root/profile
operations, source+TOV frontend preparation,
science serialization), the measured
one-process wall was **17.536732 s** and
two independent processes **9.969628 s**,
a scoped **1.759x** throughput gain.
The exact original-source scientific JSON
hashes were identical between 1/2-process
runs and a second repeated two-process run.
Peak full-worker RSS: approximately
**102--108 MB**. There is **no**
evidence here for total Python application
memory reduction from the standalone
Fortran memory measurement.

**Reproducibility sources:**
- `data/cowling_c16_source_continuation_20261121/PILOT_SCAN.json`
- `data/cowling_c16_source_continuation_20261121/REAL_P2_SHAPE/M2.json`
- `data/cowling_c16_source_continuation_20261121/REAL_P2_SHAPE/M3.json`
- The above two `REAL_P2_SHAPE` JSON files are the
  **single canonical science payload** for Git.
  The final serial and independent parallel trial
  copies have identical hashes; the earlier pre-overlap
  trial has its own distinct schema and hash.
  Each original `SCIENTIFIC.sha256`, timing,
  memory and run log remains versioned. The raw
  duplicate JSON files are excluded from Git
  to avoid repeating over 20,000 textual lines,
  but remain intact in the separately SHA-checked
  local Dropbox full replay backup.
- `data/cowling_c16_source_continuation_20261121/SCANNER_PERFORMANCE.json`
  (created by benchmark)
- `data/cowling_c16_source_continuation_20261121/RECEIPT.json`
- Code: `scripts/run_c16_source_density.py`,
  `scripts/bench_c16_source_parallel.sh`,
  `scripts/build_c16_fortran_original_f_p1.sh`.

## 6. Scientific verdict and explicit uncertainties

Local original-source f/p1 node-identified continuation
and Fortran phasewise W-shape diagnostic:
**SCOPED_VERIFIED** for the ten models and
specified original Maxwell slow/rapid conversions.
The mode order did not invert within ±2% of
each source-locked anchor.

Not established:
- full stellar-branch continuation across
  transition onset/termination, turning points
  or any extended M/R sequence;
- absence of avoided crossing outside the
  tested local range;
- rigorous physical (self-adjoint or covariant)
  mode-inner-product continuity, which depends
  on model closure and matching conditions;
- frozen composition bulk g-modes,
  full-GR metric perturbation, damping rates;
- absolute eigenfrequency error bound.
  C15 explicitly observed TOV background
  sample effects (especially M3) that exceed
  some fixed-background Cowling error controls.

**Next bounded physics stage:** increase
central-density coverage toward transition
onset, vary TOV source-grid and ODE settings
at multiple density points, validate
actual physical eigenfunction overlap
definitions and publish an honest
branch-admission/failure map. Do not
treat a sign flip or matching-frequency
smoothness alone as mode identity.

## Primary literature

Tonetto and Lugones, slow/rapid physical
Maxwell interface conditions:
https://arxiv.org/abs/2003.01259

Sotani, low-mass cold neutron-star
f/p1 avoided crossings:
https://scispace.com/papers/gravitational-wave-asteroseismology-for-low-mass-neutron-m7bkvu8bhl

Rodriguez, Ranea-Sandoval, Chirenti et al.,
pitfalls of mode classification around
avoided crossings in evolving
proto-neutron stars (not equivalent to
these cold hybrid stars):
https://scispace.com/papers/three-approaches-for-the-classification-of-protoneutron-star-1e2ilmv6

Netlib original VODE documentation
(local RTOL/ATOL control not global
eigenfrequency error):
https://www.netlib.org/ode/vode.f
