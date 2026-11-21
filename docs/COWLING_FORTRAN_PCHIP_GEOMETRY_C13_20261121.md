# NS_Cowling C13: Fortran phase-local PCHIP geometry and source-locked convergence

Date 2026-11-21 KST. Parent HEAD:
`a709208588367032cf34b04f502238f773625041`.
Isolated branch:
`research/cowling-native-pchip-geometry-c13-20261121`.
Bounded new coding stage: phase-local stellar background PCHIP
construction/geometry lookup only, plus scientifically necessary
verification. The complete 10-stage research sequence was not
rerun. Native subagents: **NOT_RUN** (no native spawn/delegate action
available); process-level numerical workers are not subagents.

## 1. Physics and unchanged source contract

Metric signature (-,+,+,+), r and geometrized mass in km,
`c=299792.458 km/s`, frequency f in Hz, and
`omega_geom=2*pi*f/c`. Entalphy `h=log(mu/mu_surface)`
and source radial coordinate `s=sqrt(h_start-h)`. All original
source thermodynamic rules and Maxwell one-sided interpolation
boundaries remain unchanged. Oscillation closure is still
**Gamma1=Gamma_eq**, with no frozen-composition buoyancy or
full-GR metric degrees of freedom.

The authority remains the original source-locked ZIP
`data/eos_missing_data_recalc_20261119/artifacts/EOS_ALL_REQUESTED_SOURCE_LOCKED_NESTED_20261119.zip`,
SHA-256
`299de43bb4f0a87ebfc4f1139ab3c8b51c344937ec3637cdff644e9ca8ab2c04`.

Internal source EOS contracts:
M2 Maxwell SHA-256
`59b7cb4671347c33c25133e6497e8be346b4e5b4f871eece9bcaf59155fd8765`;
M3 CSS SHA-256
`f23fb7020c9de1b65c99ac37690840f7a4d6611f3372409b015d8c3bd7e32f3b`.
Neither original source was modified.

Source-exact c_s² already uses original analytic laws in the
Fortran RHS and strong source derivative knot selection is
native Fortran (C11/C12). This C13 change does not alter those
physical derivatives. The original slow/rapid junction,
surface BC and regular center solution are retained.

## 2. Mathematical construction and code change

New `fortran_eos/src/compose_eos_cowling_pchip.f90`
exports two C-compatible entrypoints:
`cowling_native_pchip_slopes(n,s,values,derivatives,status)`
and `cowling_native_pchip_values(n,s,values,derivatives,nq,queries,out,status)`.
All real arrays use IEEE double and Fortran column-major
5xN layout. Inputs and outputs are validated for finiteness,
strict knot ordering, supported range and non-extrapolation.
Invalid calls zero the outputs and return nonzero status.

The slopes use the documented **Fritsch–Butland
shape-preserving PCHIP** derivative. For adjacent secants
`delta_(i-1)`, `delta_i`, set interior derivative zero on
a zero or sign-changing secant. Otherwise

```math
d_i = (w1+w2)/(w1/delta_(i-1) + w2/delta_i),
w1=2*h_i+h_(i-1), w2=h_i+2*h_(i-1),
h_i=s_(i+1)-s_i.
```

The endpoint derivative uses the SciPy-compatible
one-sided three-point estimate
`((2*h0+h1)*delta0-h0*delta1)/(h0+h1)`,
with a sign-preserving zero limiter and a 3*delta0 cap
when adjacent secants change sign.
At arbitrary interior query s, the Fortran code evaluates
the existing cubic Hermite form in its enclosing cell.
A cached cell hint and binary search fallback bound lookup
cost, while rejecting unsupported extrapolation.

The C ABI is registered in
`python_frontend/compose_eos/stellar_cowling_native.py`.
Within `NativeSharpCowlingDVODE`, **only**
`coefficient_backend="native_source_cs2"` now:
- constructs its five-component PCHIP slopes with Fortran,
  without instantiating `scipy.interpolate.PchipInterpolator`;
- retains one phase per side of Maxwell, never interpolates
  across the physical density jump;
- uses the Fortran field interpolator for the profile's
  physical radius; its phase dictionary has `pchip=None`;
- uses the existing native original-EOS chemical-potential
  strong-knot selector in frontend *setup*, so Python no
  longer computes c_s² differences to classify source
  discontinuities in this backend;
- continues calling fully native C12 Brent/DVODE/original
  EOS c_s² for root polishing.

All other coefficient backends (`pchip`, `hybrid_cs2`,
`direct_source`) preserve their existing Python/SciPy
reference behavior and their separate physics contracts.

Still NOT fully Fortran-only: Python prepares/validates EOS
JSON inputs and TOV reference samples, builds the nonuniform
one-sided phase sampling grid and transfers the arrays,
and handles final profile normalization/node classification.
Therefore 'Fortran-only numerical application' must NOT
be claimed. This change specifically removes native-source
PCHIP coefficient generation and profile geometry evaluation
from Python.

Primary numeric source references:
- SciPy PchipInterpolator:
  https://docs.scipy.org/doc/scipy/reference/generated/scipy.interpolate.PchipInterpolator.html
- SciPy derivative and endpoint source:
  https://github.com/scipy/scipy/blob/main/scipy/interpolate/_cubic.py
- Shape-preserving method: Fritsch & Butland, SIAM
  J. Sci. Stat. Comput. 5, 300–304 (1984).

## 3. Independent verification and test TDD evidence

RED before implementation: new public native PCHIP ABI was
undefined in the original shared library, recorded in
`RED_PCHIP_MISSING_ABI.log`. After implementing it:
- 3 synthetic PCHIP shapes: monotone, flat segments,
  oscillatory/nonmonotone, with irregularly spaced knots;
- a separate endpoint flat-to-decreasing source segment test (Fortran SIGN semantics versus SciPy zero-sign semantics), plus invalid non-finite knots/values, repeated abscissa and
  out-of-domain profile query fail closed;
- both original source M2 and M3, both physical phases:
  Fortran slopes and profile values checked against SciPy
  on actual source-locked background tables;
- new native-source solver construction with the module's
  Python PchipInterpolator replaced by a test stub that
  raises if called;
- original p1 eigenfrequencies, node identities and
  junction constraints retained.

Initial C13 focused tests: **12 passed, 0 failed in 8.60s**. An additional flat-to-decreasing endpoint-limiter regression passed. The final focused **C11/C12/C13 suite reached 33 passed, 0 failed in 32.58s**, recorded in `FOCUSED_FINAL.log`.
Original broad integrated regression across C7–C13,
smooth/star native TOV, Maxwell/slow/rapid, hybrid and
direct-source EOS callbacks and C12 native source-exact:
**119 passed, 0 failed in 189.37s**.
Preserved transcripts: `C13_FOCUSED_INITIAL.log`,
`INTEGRATED_REGRESSION.log`.

An additional completely self-contained Fortran-only program
`fortran_eos/tests/cowling_pchip_c13_standalone.f90`
was compiled `-O0 -g -fcheck=all` and `-O2`.
The source-locked C12 original-star PCHIP snapshots were
downsampled *for independent oracle testing only* to
257 knots and 97 evaluation coordinates for each of:
M2 high, M2 low, M3 high and M3 low.
All four standalone source/SciPy-oracle comparisons passed;
maximum normalized slope deviation per snapshot was
at most **2.22e-16**.
The checked and optimized builds emitted identical
scientific output bytes and rejected out-of-support
profile requests and duplicated knots.
Standalone executable observed RSS was 2.3–2.7 MB per
case; this is NOT the whole star solver's memory use.
The checked SciPy oracle arrays were exported only once
and compressed as deterministic gzip fixtures.
`scripts/export_c13_original_pchip_fixtures.py`
recreates them from the prior locked C12 files.

## 4. Real EOS convergence and limits

Keeping unchanged EOS SHA, Fortran source-exact c_s²,
Fortran Brent and strong knot schedule and native TOV,
one-sided geometry table sampling was varied. Requested
`samples` are NOT the actual phase grid size because
source breakpoints and near-boundary points are inserted.

| Requested geometry samples | M2 Maxwell rapid p1 [Hz] | M3 CSS slow p1 [Hz] |
|---|---:|---:|
| 1025 | 5143.719663516 | 3736.210609812 |
| 2049 | 5143.719647087 | 3736.210628591 |
| 4097 | 5143.719711815 | 3736.210706737 |
| 8193 | 5143.719634196 | 3736.210436364 |

All eight mode calculations passed the existing W/V node
identity 1/1, `relative_join_mismatch<1e-5`,
source-callback-count=0 and Fortran restart-count parity.
The M3 join mismatch at 4097 was ~9.57e-6
(just below the preserved gate); 8193 improved to ~1.39e-6.
But the frequencies are **not monotone with sample count**.
This evidence does NOT establish a Richardson extrapolation
law or global error bars. The signed C13 shifts from the
C12 4097 baseline are approximately +0.00010348 Hz for
M2 and +0.00000800 Hz for M3 at that nominal sample size.
No improvement in absolute frequency accuracy is inferred.

Numerical evidence:
`REAL_SOURCE_GEOMETRY_CONVERGENCE.json`, raw
`CONVERGENCE_RUN.log`.

## 5. Runtime, lightweight and parallel results

On source-locked actual background arrays, timing only
the same 5-field PCHIP **slope construction** with seven
matched-in-memory repeats:
- SciPy-to-Fortran isolated speedup ranged about
  **3.21x to 5.23x** over the tested two EOS and
  1025–8193 sample configurations.
- These tiny millisecond or submillisecond numerical
  tasks are not the main total-runtime cost; the
  full frontend phase preparation remains substantially
  dominated by source EOS/TOV sampling.
- Numeric objects are explicitly phase-local
  5xN double arrays, one reusable calculation per phase.
  No new memory-heavy interpolation object is constructed
  in native_source_cs2.

Whole two-original-star workload at 4097 requested samples,
including Python source and phase preparation, native Fortran
TOV/PCHIP/root and result serialization:
- 1 independent OS worker: **6.924944 s**;
- 2 independent OS workers: **4.380652 s**,
  measured ratio **1.581x**.
Maximum per-worker RSS was **~99–105 MB**.
Do not claim total pipeline memory reduction from a
2–3 MB independent kernel-only test.

The first process comparison failed because timing fields
were incorrectly embedded into the science JSON. Scheduler
time is not part of scientific result identity. This is
classified `OUTPUT_SERIALIZATION_NONDETERMINISTIC_METADATA`,
**not** a numerical source/root failure. The original
per-process raw artifacts are retained as evidence. The
script was corrected to send timing only to the worker
logs and store deterministic observables in science JSON.
A clean fresh 1/2-process rerun then matched exact
per-star scientific JSON bytes and SHA-256.

Validated canonical science manifest checksum:
`077efa2fa53d02562022ab91358b0945fd5c63c34b609fb478d12340909a920e`.

## 6. Verdict and next constrained objective

**PASS_SCOPED_FORTRAN_PCHIP_GEOMETRY_WITH_CONVERGENCE_OPEN**.

- Implementation and backward-compatible data contracts:
  tested and verified on the fixed original M2/M3 stars.
- Normalized source-derivative identities, physically
  meaningful junction requirements and EOS provenance:
  inherited from prior C11/C12, regression tested again.
- Per-source PCHIP interpolation is portable to Fortran,
  with SciPy reference parity on nonuniform/flat/oscillatory
  synthetic data and actual star phases.
- Total mode-frequency convergence across finer phase grids
  remains **UNRESOLVED**; nonmonotone differences and
  proximity to join acceptance indicate a larger joint
  geometry/ODE error-budget study.
- No independent study of full f/p1 mode branch continuation,
  multiple Maxwell transitions, frozen-composition g-modes,
  GPU execution, or a complete Fortran-only EOS importer
  was performed in this bounded slice.

Next scientific coding work should fix the remaining
source/TOV sample materialization boundary and vary
`rtol`, `max_step_h`, profile nodes and source
geometry resolution independently while tracking
continuous mode branches across central density.
Do not inflate acceptance gates or declare the
full EOS regime certified on two stars.
