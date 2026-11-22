# NS_Cowling C17: extended original-EOS Maxwell f/p1 branch and small-core onset audit

**Date:** 2026-11-22 KST.
**Base/SSOT commit:** b956e0a21cc31623d5034ab9f0cf9153b46a1d41
**New isolated branch:** research/cowling-c17-extended-mode-branch-20261122
**Scope:** one constructive numerical-physics feature and bounded verification following C16.
No full C1–C10 restart. Native skill-bound subagent dispatch: NOT_RUN (unavailable).
Independent OS numerical worker processes are not agents.

## 1. Physical model and protected authority

The original immutable M2 Maxwell and M3 CSS effective cold EOSs are loaded
without refitting or changing any source knot or transition density.
Frozen archive:
data/eos_missing_data_recalc_20261119/artifacts/EOS_ALL_REQUESTED_SOURCE_LOCKED_NESTED_20261119.zip

SHA-256:
299de43bb4f0a87ebfc4f1139ab3c8b51c344937ec3637cdff644e9ca8ab2c04

Original nested source identities:
M2 Maxwell: 59b7cb4671347c33c25133e6497e8be346b4e5b4f871eece9bcaf59155fd8765
M3 CSS: f23fb7020c9de1b65c99ac37690840f7a4d6611f3372409b015d8c3bd7e32f3b

Signature (-,+,+,+), physical radius and geometrized mass in km,
energy density and pressure from MeV/fm^3 converted into geometrized km^-2
for the stellar background, and
omega_geom = 2*pi*f_Hz/(c = 299792.458 km/s) for fluid shooting.
Cold equilibrium barotrope; Gamma1=Gamma_eq.
No independent frozen-composition adiabatic derivative,
bulk g-mode family, full-GR polar metric perturbation or GW damping.

The original Maxwell pressure/geometry identities and slow/rapid
fluid-interface matching are retained:
H=epsilon+P, a=omega_geom^2 exp(lambda-2Phi)*r, b=Phi_prime,
T=H*(a*v+b*w).
SLOW: [w]=0 and [T]=0. RAPID: [v]=0 and [T]=0.
The native Fortran root must retain the original absolute normalized
matching determinant and relative eigenfunction join limit of 1e-5.
No acceptance threshold was relaxed.

## 2. Original-source phase transition onset: finite derived result

From original one-sided energy densities and pressure at transition,
the standard infinitesimal high-density-core Maxwell density-jump
diagnostic (Lindblom 1998, Phys Rev D 58, 024008) is:

Delta=(epsilon_plus-epsilon_minus)/(epsilon_minus+P_t);
Delta_crit=(epsilon_minus+3*P_t)/(2*(epsilon_minus+P_t)).

Numerically checked against original source EOS:
- M2: P_t=0.014307398278768153 MeV fm^-3,
  epsilon_minus=8.572084498151984,
  epsilon_plus=16.33826505200197,
  Delta=0.9044754359602759,
  Delta_crit=0.5016662875921976,
  Delta-Delta_crit=+0.40280914836807835.
- M3: P_t=0.5657052249953693,
  epsilon_minus=76.94982713409907,
  epsilon_plus=89.34565739599732,
  Delta=0.15991414732822803,
  Delta_crit=0.507297959618915,
  Delta-Delta_crit=-0.34738381229068693.

The exact source transition high-phase density relative to the original
central-density anchor is:
- M2: n_high/n_anchor = 0.05706326070156743;
- M3: n_high/n_anchor = 0.48894522196075574.

These quantities were evaluated in the new Fortran module and
checked against independent analytic arithmetic. The sign tests describe
the conventional infinitesimal-core threshold only, and do NOT
constitute a global radial-stability theorem or guarantee the
slow-conversion dynamical boundary condition. M2 is NOT followed
to its transition onset here: the lowest chosen density factor is
0.495, still far above M2's 0.057 onset. M3's lowest computed
density factor 0.495 lies approximately 1.24% above the high-phase
transition-onset density.

Primary source: Lee Lindblom, Phase Transitions and the Mass-Radius
Curves of Relativistic Stars, https://arxiv.org/abs/gr-qc/9802072 .

## 3. Constructed Fortran numerical features

New Fortran module:
fortran_eos/src/compose_eos_cowling_extended_branch.f90

C ABIs:
- cowling_native_maxwell_onset_threshold:
  typed original-source (epsilon_minus,epsilon_plus,P_t) input;
  returns Delta, Delta_crit, signed margin and above-criterion flag;
  nonphysical, nonfinite input zeros all numerical outputs and rejects.
- cowling_native_select_mode_branch_overlap:
  frequency, W/V radial node counts, prior root admission statuses,
  high-phase and low-phase absolute W-shape correlations, previous
  frequency and bounded maximum relative frequency change.
  Returns ONE successor only if it is the UNIQUE admissible candidate
  with the specified paired node counts, relative frequency step
  at most 12%, and separately high and low phase W-overlap >= 0.95.

Typed status outcomes: 0 one admitted successor, 3 multiple eligible
successors, 4 no paired-node/frequency candidate, 5 such a candidate
exists but fails the phase-overlap criterion, 1 invalid inputs.
No branch is stitched over a rejected predecessor.
No frequency or shape candidate is silently substituted by a
nearest-neighbour guess.

The 0.95 W-shape overlap threshold is an explicitly configured
NUMERICAL classification diagnostic, not a physical conserved
mode inner-product threshold. W alone is not a sufficient
wave-mode basis across an avoided crossing.
The preceding C16 true phasewise W-overlap is computed by the
previous Fortran source and is not modified.

The original source-exact Cowling physics, DVODE, Maxwell transfer,
Fortran root and PCHIP implementations are UNCHANGED.
The Python ABI registration
python_frontend/compose_eos/stellar_cowling_native.py
only adds the two new native signatures.
Original C16 defaults and science artifacts are untouched.

New Python frontend wrappers and evidence-only scripts:
- scripts/run_c17_extended_density.py:
  one original EOS source per OS process, fixed 111-point frequency
  scan for each central density, then full native Fortran root +
  Maxwell join + phasewise node validation, and unique two-sided
  W-overlap branch selection. Python remains responsible for source
  JSON loading, TOV background sampling, and serialized reports.
  The legacy Python PCHIP query of completed TOV profiles is
  still part of the pre-existing frontend. Claiming the
  *entire* scientific application is Fortran-only is false.
- scripts/bench_c17_extended_parallel.sh:
  independent 1- or 2-process execution, with BLAS/OpenMP single-thread
  limits to avoid unsafe calls to the legacy ODEPACK DVODE COMMON arrays.
- scripts/refine_c17_density_midpoints.py:
  one ONE-TIME finite half-step attempt at four interrupted edges,
  with the original Fortran mode and W overlap gates unchanged.
  No recursive refinement, no acceptance-gate inflation.
- scripts/check_c17_tov_grid.py:
  independent original TOV source grid 4097 and 8193 sensitivity
  at four chosen stellar backgrounds; Cowling grid and local root
  settings fixed; failures are typed not imputed.
- scripts/summarize_c17_extended.py:
  frontend-only parsing, aggregation and deterministic CSV/JSON.

## 4. Extended original EOS 40-star finite computation

Original anchor central density:
M2 = 0.30385128926615257 fm^-3
M3 = 0.192510767556007 fm^-3

The density factors per star are:
0.495,0.50,0.52,0.55,0.60,0.65,0.70,0.75,0.80,0.85,
0.90,0.95,0.98,0.99,1.0,1.01,1.02,1.05,1.10,1.15.

Original TOV controls:
RTOL=3e-10, ATOL=3e-11, max_step_h=.004,
requested source samples 2049.
Cowling controls:
RTOL=1e-9, ATOL=1e-12, max_step_h=.002,
root_xtol_hz=1e-5, requested phase geometry samples=4097,
mode profile points=360. Source derivative strong-knot
restart contrast criterion unchanged (0.01).

The wide fixed scan comprises 111 positive frequency values,
500 to 6000 Hz; scanning is NOT exhaustive over every
potential mode interval. In total 40 stellar backgrounds,
4440 original-source Fortran determinant shots.
All scan results are checked and serialized, with raw scans
losslessly compressed as deterministic gzip sidecars
instead of replicating the whole data in primary JSON.
The original EOS physics produces a total of 287 signed-
determinant candidate intervals, each handled by native Fortran
Brent and independent full Maxwell eigenfunction reconstruction.
A sign change alone is not certified as a physical mode.

Accepted local paired f/p1 nodes and mode order:
- M2: full both-mode branch at 16/20 specified stellar densities;
  f tracked 20/20; p1 tracked 16/20.
- M3: both-mode branch at 7/20 specified stellar densities;
  f tracked 16/20; p1 tracked 7/20.

Every connected point was constrained to f:(W,V)=(0,0),
p1:(W,V)=(1,1), original interface continuity,
and unique root/overlap gate.
Every failure is a typed branch break, with no physical
frequency manufactured for an unadmitted segment.

Notable bounded events:
- M2 p1 could initially be continued down to factor .60,
  while factor .55 was blocked by the single-step 12%
  frequency change. The other lower-density factors
  remain unlinked in the canonical 20-factor scan.
- M3 p1 connects from the anchor down to .90; .85
  fails a numerical nodal/frequency predecessor test.
- M3 at factor 1.05 has a candidate p1 whose *full*
  absolute W-overlap relative to 1.02 was approximately
  0.999988734; the HIGH-PHASE-only overlap was
  0.919350741, below the retained 0.95 threshold.
  This triggered a typed PHASE_SHAPE_OVERLAP_BELOW_THRESHOLD.
- M3 at factor 0.495, only about 1.24% above the
  high-phase transition onset, has radius about 42.492 km
  and 63 signed scan candidates. The proliferation is
  an important underresolution/physical-mode-disambiguation
  warning, not proof of a huge number of independently
  certified stellar modes.

## 5. Finite, non-recursive density midpoint repair experiment

Exactly one intermediate density was inserted per pre-selected
failed edge. Original Fortran mode and physical source gates
were unchanged:

| Edge | Midpoint | Result |
|---|---:|---|
| M2 p1, factor .60 to .55 | .575 | TWO-HOP ADMITTED |
| M3 p1, factor .90 to .85 | .875 | ENDPOINT STILL BLOCKED |
| M3 p1, factor 1.02 to 1.05 | 1.035 | TWO-HOP ADMITTED |
| M3 f, factor .60 to .55 | .575 | ENDPOINT STILL BLOCKED |

M2 p1 new midpoint f_p1 = 3257.12399492596 Hz, reached
the endpoint at factor .55 with p1=3038.301675735732 Hz.
M3 p1 midpoint at factor 1.035:
3932.6262054832187 Hz, endpoint 1.05:
4016.313086414399 Hz, admitted in the two shorter steps.

Thus some initial C17 branch breaks were due to
coarse density steps crossing the acceptance gate, not
demonstrated physical mode disappearance. The two failures
remaining after one midpoint are preserved, without
a recursive, unchecked refinement loop.

The canonical original 40-star status map remains the
unmodified directly sampled result. Midpoint continuations
are separately versioned and do not rewrite those observations.

## 6. Independent TOV background grid crosscheck

Source EOS input is locked and the Cowling settings
remain fixed. Four original-star points were
recomputed with the original Fortran TOV frontend
at requested TOV EOS sample grids 4097 and 8193
instead of 2049: M2 at .60 and 1.0, M3 at .90 and 1.0,
with f and p1 mode attempts.

Out of 16 additional original-source Fortran
mode computations, 15 were admitted and ONE was rejected:
M3 at .90, p1, TOV grid 4097.
Selected shifts relative to the originally admitted
2049-source TOV at fixed center/physics:
- M2 .60 p1: -0.00014917 Hz (4097),
  -0.00027727 Hz (8193);
- M3 1.0 p1: -0.00154627 Hz (4097),
  -0.00185171 Hz (8193);
- M3 .90 p1: rejected at 4097;
  -0.00126548 Hz at 8193.

These finite differences are NOT a rigorous total
eigenfrequency error bound, Richardson rate or
proof of global numerical convergence. The rejected
p1 is not replaced by interpolation.
Raw data: TOV_GRID_CROSSCHECK.json.

## 7. Tests, independent Fortran and resource measurements

New focused tests:
python_frontend/tests/test_c17_extended_branch.py
and python_frontend/tests/test_c17_original_science.py.
A failing preimplementation ABI test was preserved
as RED_NATIVE_OVERLAP_SELECTOR.log, then
the scoped Fortran unit and original-source
finite science checks passed (14 passed in 0.44s
in the first source-record validation run).

The final fresh selected C7--C17 integrated
regression, including the C16 physics paths,
TOV/Cowling Fortran backends, Maxwell interfaces,
all prior sharp/smooth reference variants and
new C17 unit/science tests, produced
**166 passed, 0 failed in 211.67 s**.
Original transcript:
data/cowling_c17_extended_branch_20261122/INTEGRATED_REGRESSION.log.

The independent Fortran-only
fortran_eos/tests/cowling_c17_original_onset_branch.f90
checks original M2/M3 phase-onset quantities,
unambiguous selection, true ambiguous double
candidate refusal and insufficient-overlap refusal.
Changed project sources compiled under
-O0 -g -fcheck=all and -O2,
with byte-identical output in this test.
This standalone bounded Fortran-only program used
2488 KB RSS; this does not measure the full
TOV/Cowling shooting memory or certify the
legacy Netlib F77 vendor bounds.

Two actual 40-star runs including source JSON EOS,
TOV/geometry frontend, source-exact Fortran
scan/Brent/Maxwell/nodal/W overlap, and serialization:

| Workers | Wall time | Relative speed |
|---|---:|---:|
| 1 independent OS process | 63.060987 s | 1.00 |
| 2 independent OS processes | 34.338328 s | 1.84 |

The canonical science JSON + compressed source scan
hash manifest was byte-identical between 1- and
2-process runs. Per-worker maximum RSS was
approximately 113–116 MB. This speedup is scoped
to two independent stars, and does not justify
shared-thread use of non-reentrant ODEPACK.

The lightweight layout keeps one canonical
original source science JSON and one compressed
scan array per star in Git; identical repeated
parallel/serial science JSON copies can be
restored from the user's separate full replay
Dropbox archive. Timing and status logs and
canonical hash records remain in Git.
Do not claim a 1.84x speedup for every
single-star mode or a full-pipeline memory
reduction.

## 8. Verdict, nonclaims and next bounded research

**C17 status:** SCOPED_ORIGINAL_EOS_EXTENDED_BRANCH_MAP_WITH_TYPED_BLOCKS.
This is a numerical-physics construction and failure
map, not a general stability or global mode
continuation theorem.

Outstanding:
- M3 high-density-core W shape differs enough
  from a large density step to cross a 0.95
  arbitrary numerical gate; midpoint repair
  works on one selected segment, but there is
  no physical invariant inner-product proof.
- M3 near onset is dominated by many candidate
  intervals and large changes in radius. The
  sub-percent onset region needs a separate
  physically justified small-core asymptotic
  treatment before any stability or avoided-
  crossing claim can be promoted.
- M2 onset is at factor 0.057, far outside the
  actual explored [0.495,1.15] factors.
- C15/C17 demonstrate nontrivial sensitivity to
  TOV background sampling. No absolute
  eigenfrequency error certification yet.
- Full GR polar and g-mode physics are out of scope.

Next bounded candidate stage:
**C18-ONSET-MODE-INNER-PRODUCT-01**:
derive and verify physically meaningful
phase-conversion-dependent mode overlap and
small-core asymptotics; independently assess
f and p1 over denser near-onset equilibrium
stellar sequences, and develop a fixed source
C16/C17 branch identity/rejection audit.
Do not inflate acceptance gates.

## Primary physical references

Lindblom, Phase Transitions and the Mass-Radius
Curves of Relativistic Stars (1998),
https://arxiv.org/abs/gr-qc/9802072

Sotani, Gravitational wave asteroseismology for
low-mass neutron stars (2020), avoided f/p1 modes,
https://arxiv.org/abs/2008.09839

Tonetto and Lugones (2020), slow/rapid Maxwell
perturbation matching:
https://arxiv.org/abs/2003.01259

Original Netlib DVODE:
https://www.netlib.org/ode/vode.f
