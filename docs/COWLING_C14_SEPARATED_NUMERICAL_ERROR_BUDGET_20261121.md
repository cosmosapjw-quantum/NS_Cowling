# NS_Cowling C14: Source-locked, 3-axis Fortran error-budget experiment

Date: 2026-11-21 KST.
Parent baseline: `e0f6209f5c54b08b80959f9bf1dfe5742240bc2f`.
Isolated branch:
`research/cowling-c14-source-error-budget-20261121`.
PhysMath scope: one bounded C8-style *numerical convergence and
implementation* slice following C13, NOT automatic replay of the
entire research-track stages.
Native skill-bound subagents: **NOT_RUN** (host has no native
spawn/delegate tool). The OS process-parallel calculations are
ordinary numerical workers, not research subagents.

## 1. Scientific problem

The original source-locked C13 runs for M2 Maxwell and M3 CSS
produced p1 frequencies that **did not converge monotonically**
as background geometry sample counts increased. Avoid attributing
this instability solely to the PCHIP algorithm, solely to DVODE, or
to physical mode changes without an independent experiment.

In this experiment we hold the same original cold EOS
and fixed stellar TOV profile while varying **one numerical control
family at a time**, separately for Adams and BDF:

1. Phase-local geometry PCHIP *coarsening* (factor 1,2,4,8 from
   the original full sampled source star); preserve exact
   physical Maxwell phase endpoints and source EOS landmarks
   already present in that fixed high-resolution geometry;
2. DVODE local-error relative tolerance RTOL =
   3e-8, 1e-9, 3e-10, 1e-10;
3. Maximum DVODE s-step HMAX relative to the original
   `base_max_ds`: 2,1,1/2,1/4.

Formally the measured eigenfrequency is
`f(N_geom, tau_DVODE, HMAX_s; EOS,TOV,Maxwell,mf)`.
For each axis the reported finite difference changes that axis
relative to a fixed reference control while the other two are
held unchanged. Interactions such as
`delta_geom(tau_a)-delta_geom(tau_b)` have **not**
been measured here. Thus the three spans cannot be added
as independent components of a total error bound.

The third variable, `s=sqrt(h_start-h)`, is dimensionless;
a change in maximum step *s* is not equivalent to a fixed physical
radius step. Absolute DVODE tolerance is fixed at `1e-12`,
and eigenvalue bracket and root stopping tolerance are fixed.
For every root the original source EOS law, equilibrium closure,
strong chemical-potential derivative knot restart criterion and
slow/rapid Maxwell junction are held unchanged.

**Important:** these experiments investigate the *conditional*
frequency mapping f(N_geometry,RTOL,HMAX) at a fixed TOV profile.
They do **not** independently vary TOV discretization,
the original microscopic EOS, source barotropic equilibrium
closure, crust-core interpolation, branch history or full GR
perturbations.

## 2. Physics and provenance

Metric signature (-,+,+,+), `r` and geometrized
stellar mass in km, `n` in fm^-3, `mu` in MeV,
pressure/energy density in MeV/fm^3.
The perturbation angular frequency is
`omega_geom=2*pi*f_Hz/(299792.458 km/s)`;
do not silently set c=1 when relating root frequencies in Hz.

Cold equilibrium barotrope, **Gamma1=Gamma_eq**,
one sharp Maxwell jump. Interface conditions inherited
from C12 in terms of scaled Cowling variables (w,v):

```text
H=epsilon+p; a=omega_geom^2 exp(lambda-2Phi)*r; b=Phi_prime
T=H*(a*v+b*w)

slow: [w]=0 and [T]=0
rapid: [v]=0 and [T]=0
```

Only the original C12 model brackets are used:
- M2 Maxwell, n_c=0.30385128926615257 fm^-3,
  M about 1.4 M_sun, **rapid p1**,
  root bracket [5110,5170] Hz, 39 source derivative restarts;
- M3 CSS, n_c=0.192510767556007 fm^-3,
  M about 0.7 M_sun, **slow p1**,
  root bracket [3710,3760] Hz, 40 source derivative restarts.

Frozen authority:
`data/eos_missing_data_recalc_20261119/artifacts/EOS_ALL_REQUESTED_SOURCE_LOCKED_NESTED_20261119.zip`
SHA-256:
`299de43bb4f0a87ebfc4f1139ab3c8b51c344937ec3637cdff644e9ca8ab2c04`.

C12 immutable source-EOS/geometry input file uncompressed SHA-256:
- M2:
  `eee499eb50de4300f9987e05bf46a4fee3ddcbfb2e9b9d99deb9ac8a423c35e5`
- M3:
  `f6a7abb566ae68b8ac5f838d410771823634a98024007ca02264bd6e44a0d3b1`.

Original physical EOS contract hashes are retained:
M2 `59b7cb4671347c33c25133e6497e8be346b4e5b4f871eece9bcaf59155fd8765`;
M3 `f23fb7020c9de1b65c99ac37690840f7a4d6611f3372409b015d8c3bd7e32f3b`.
No original source or baseline was modified.

## 3. Concrete Fortran implementation

New public C ABI module:
`fortran_eos/src/compose_eos_cowling_error_budget.f90`:
`cowling_native_sharp_source_sweep`.

This routine accepts an array of per-run
`(RTOL,HMAX,mf)` controls, where mf=10 selects Adams and
mf=22 selects BDF. For every control it calls the previous
original source-exact Fortran Brent solver with a fresh DVODE
history at each actual selected source derivative knot,
and produces frequency, normalized signed matching determinant,
iterations, native RHS evaluation count, source restart count
and a **per-run typed status**. Bad controls get zero certified
numerical outputs, without aborting other valid controls.

The standalone Fortran program
`fortran_eos/tests/cowling_c14_error_budget.f90`:
- Reads the unchanged C12 source EOS/stellar background file,
  including the physically distinct high and low phases;
- Computes the phase geometry PCHIP slopes in Fortran;
- At geometry coarsening factors 1,2,4,8, uses original
  **source mu-defined landmarks** and Maxwell endpoints to
  pin immutable physical breakpoints, while subsampling
  the remaining background geometry points;
- Executes the same original source EOS, Maxwell transfer,
  Fortran Brent and DVODE in each controlled shot;
- Leaves invalid/bracket-rejected roots identified by a
  nonzero status and with no fabricated frequency;
- Emits deterministic plain-text scientific rows
  `C14_ROW<AXIS> factor n_hi n_lo mf rtol maxstep_s
  frequency_Hz det iterations RHS_calls knots status`.

The fixed original M2 sample set has 12124 high-phase and
2029 low-phase points; M3 has 6601 and 2710.
After a stride of 8, while retaining source landmarks,
the counts became M2 3904+501 and M3 1749+628.
Because some original source knots remain pinned,
the factors 1,2,4,8 are **not exact global h-grid doubling**
and cannot justify Richardson extrapolation.

The production `NativeSharpCowlingDVODE` has a new
`sweep_native_numerical_controls` **frontend-only helper**
that marshals controls into a single native Fortran batch,
while preserving the original physics and root behaviors:

```python
controls = [
    dict(rtol=1e-9, max_step_s=base_step_s, method="adams"),
    dict(rtol=3e-10, max_step_s=base_step_s, method="bdf"),
]
results = solver.sweep_native_numerical_controls(
    (5110.,5170.), conversion="rapid", controls=controls)
```

The helper returns each result with typed `status`;
failed controls have `frequency_hz=None` and do not masquerade
as certified eigenmodes. It does not call Python for root
iterations, EOS derivatives or original source sound speeds.
Existing `solve_mode` defaults, source-knot interface
and all other root backends are preserved.

## 4. Numerical results and interpretation

The Fortran-only controlled run on each original star
produced 24 valid controls per star (48/48 successful roots),
eight per axis, four Adams and four BDF.
All retained the original source knot count and
`abs(normalized_determinant)<=1e-5`.

The intra-axis **frequency spans** in Hz are:

| Original star | Numerical variable | Adams (MF=10) | BDF (MF=22) |
|---|---|---:|---:|
| M2 rapid p1 | Geometry coarsening | 0.000109569 | 0.000013441 |
| M2 rapid p1 | DVODE RTOL | 0.000055401 | 0.000229606 |
| M2 rapid p1 | DVODE max s-step | 0.000108295 | 0.000070884 |
| M3 slow p1 | Geometry coarsening | 0.000221548 | 0.000041221 |
| M3 slow p1 | DVODE RTOL | 0.000539952 | 0.000512425 |
| M3 slow p1 | DVODE max s-step | 0.000069154 | 0.000025698 |

These are observed **finite-control sensitivity ranges**,
not global truncation errors, 1-sigma uncertainties
or strict upper bounds on the true eigenfrequency.
Strong nonmonotone control dependence occurs for several
axes, especially M3 relative-tolerance and varying
maximum steps. No empirical asymptotic power law was
admitted. Agreement between Adams/BDF at one tight tolerance
is useful evidence but does not prove geometric/TOV
or microscopic EOS accuracy.

At the baseline rtol=1e-9 and full geometry,
the original M3 Adams normalized matching determinant
is about `9.57e-6`, just below the original `1e-5`
acceptance threshold. This demonstrates sensitivity
near a *preexisting* acceptance boundary; the threshold
was NOT relaxed to admit a candidate.

This standalone experiment evaluates roots on fixed
original-p1 **frequency brackets**. It does not
reconstruct mode eigenfunctions for all 48 control settings.

A distinct bounded follow-up physical profile experiment
reconstructed **12 full original-source eigenfunctions**:
M2 rapid p1 and M3 slow p1, each at three RTOL controls
(3e-8,1e-9,1e-10) and both Adams and BDF methods,
at the full original phase geometry and fixed max s-step.
All 12 retained physical node counts **W=1,V=1**
and met the original `relative_join_mismatch<1e-5`.
The M3 Adams/1e-9 case returned
`relative_join_mismatch=9.565274e-6`: marginally
inside the unchanged historical acceptance limit.
See `MODE_ADMISSIBILITY.json`,
`MODE_ADMISSIBILITY_RUN.log`, and
`scripts/check_cowling_c14_mode_admissibility.py`.
This does not certify 48 eigenfunction identities nor
global branch continuation.

Numerical raw source:
`data/cowling_c14_source_error_budget_20261121/OPT_1/{M2,M3}.science`.
Parsed results and per-axis differences:
`C14_ERROR_BUDGET.json`.
The Python summarizer `scripts/summarize_cowling_c14_budget.py`
only parses Fortran output; it does not solve eigenmodes.

## 5. Build, negative tests, reproducibility and memory

A TDD RED record shows the C ABI was absent before
implementation (`RED_MISSING_NATIVE_BATCH.log`).
The new two-source Python frontend ABI tests
`test_cowling_c14_native_error_budget.py`
check valid paired Adams/BDF outcomes, an invalid
RTOL isolated from valid neighboring controls,
no Python root callback, source-knot parity,
and invalid global bracket fail-closed behavior.
The initial native ABI result was **4 passed**.
A separate RED test then demonstrated that the documented
frontend batch method did not yet exist, after which it was
implemented as a thin Fortran call. The fresh scoped C14
suite reached **6 passed** (see `FOCUSED_FINAL.log`).
The Python helper performs only parameter marshaling and
result formatting; it is not an alternative numerical solver.

A fresh comprehensive regression of the relevant C7–C14
smooth and sharp Cowling, source-EOS, Maxwell, PCHIP,
TOV and physical-mode test paths after the final frontend
addition produced **126 passed, 0 failed in 197.37 s**.
See `INTEGRATED_REGRESSION_FINAL.log`.
The earlier 124-pass run was performed before the final
new public frontend method tests and is retained separately
as preliminary historical evidence.

Fortran numerical checks were performed under two build modes:
- modified project source `-O0 -g -fcheck=all`;
- modified project source `-O2`;
both link untouched Netlib ODEPACK compatibility F77
source compiled with `-O2 -std=legacy`.
The results are **byte-identical across optimized,
bounds-checked, serial and process-parallel executions**,
all with the same original C12 inputs.
This does not certify the memory safety of the entire
legacy F77 vendor distribution.

48 source-locked original-star roots, no Python runtime
inside the Fortran executables:

| Fortran execution | Total wall | Relative to optimized serial |
|---|---:|---:|
| Optimized 1 OS process | 1.751043 s | 1.00x |
| Optimized 2 OS processes | 0.898241 s | 1.95x |
| Project-bounds-checked 1 OS process | 2.382852 s | 0.73x |

Observed maximum RSS: about **7.9–10.3 MB per
Fortran numerical process**. This excludes original EOS
materialization, TOV construction and geometry sampling,
which still have Python frontend preparation costs.
The optimized throughput gain is a **single measured
two-star original-data run**, not a portable scaling
guarantee; DVODE COMMON/global coefficient arrays remain
NONREENTRANT and must never be shared by threads.

`scripts/build_cowling_c14_error_budget.sh` rebuilds the
standalone binary.
`scripts/bench_cowling_c14_parallel.sh` reads archived
original C12 .dat.gz source inputs, decompresses outside
the timed run, and executes separate OS processes;
the generated temporary uncompressed copies are not
retained in the Git evidence folder.
`scripts/summarize_cowling_c14_budget.py` validates
scientific row count, source SHA and status, and
produces reproducible JSON.

## 6. Direct reproduction

After checking out the recorded C14 branch with its original
C12 compressed source and stellar geometry files:

```bash
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1
export PYTHONPATH=python_frontend
bash scripts/build_cowling_c14_error_budget.sh /tmp/ns_c14_checked/c14_budget checked
bash scripts/build_cowling_c14_error_budget.sh /tmp/ns_c14_optimized/c14_budget optimized
bash scripts/bench_cowling_c14_parallel.sh \
    /tmp/ns_c14_optimized/c14_budget \
    data/cowling_fortran_source_exact_root_c12_20261121 \
    2 /tmp/ns_c14_parallel_replay
python3 scripts/summarize_cowling_c14_budget.py
python3 scripts/check_cowling_c14_mode_admissibility.py
python3 -m pytest -q python_frontend/tests/test_cowling_c14_native_error_budget.py --run-slow
```

The summarizer reads the committed canonical `OPT_1`
source results; to analyze a fresh batch, compare its
science rows with those fixed records, or supply them
to a copy of the parsing workflow. Two Fortran methods
are calculated on each axis for each original star.
A workstation may show different wall-time metrics,
but its checked physics should remain within the
admitted source-locked test contracts.

## 7. Verdict and next bounded task

**PASS_SCOPED_FORTRAN_SEPARATED_NUMERICAL_CONTROL_STUDY**:
the new Fortran batched solver interface and the 48-point
fixed-source numerical control grid are scientifically
and implementation-verified within their declared scope.
**GLOBAL_EIGENFREQUENCY_CONVERGENCE_UNRESOLVED**:
these variations are not a clean asymptotic plateau,
and no absolute frequency error bound, physical
branch identity across central density, or EOS-wide
stability conclusion is certified.

Recommended next: independently vary fixed stellar TOV
background integration accuracy and geometry sampling,
then carry out original-EOS central-density continuation
with node-identified p1 and f branches, and perform
genuine node-identity/mode-crossing checks.
Stop using f-digit agreement as a substitute for
physical branch validation.

Primary references:
- Netlib original DVODE documentation:
  https://www.netlib.org/ode/vode.f
  (RTOL/ATOL are **local** error controls,
  not global eigenfrequency error certificates);
- Cowling neutron-star f/p1 modes:
  https://scispace.com/papers/impact-of-the-equation-of-state-on-f-and-p-mode-oscillations-15g1elv5
