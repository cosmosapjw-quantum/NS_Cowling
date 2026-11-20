# NS_Cowling: phasewise Fortran sharp-Maxwell native-root implementation

Date: 2026-11-20 KST.
Parent commit: `436d7e4213d94696cf07defe075fed778587d12a`.
Branch: `research/cowling-sharp-native-root-phase-pchip-20261120`.
PhysMath selection: new bounded **C5 + C6** implementation/verification
slice; C8 benchmark checks are local to the new code, not a reopened
full research campaign. Native skill-bound subagents: **NOT_RUN**
(no native agent delegation provided in current tool inventory).

## Science and scope

Metric (-,+,+,+); `r` in km, `f` in Hz,
`omega=2*pi*f/(c in km/s)`, `c=299792.458 km/s`.
The equilibrium barotropic closure is still `Gamma1=Gamma_eq`.
This is **not** a frozen-composition/nonbarotropic or full-GR mode solver.

A single Maxwell jump at common geometry and pressure separates the
high-density inner phase from the low-density outer phase.
Each phase has its own strictly monotone PCHIP coefficient table
in `s=sqrt(h_start-h)`, with no interpolation across the jump.

The scaled eigenvector `(w,v)` propagates with the original two-component
relativistic Cowling DVODE ODE. At the sharp interface define
`H=epsilon+p`, `a=omega² exp(lambda-2Phi) r` and `b=Phi'(r)`.
The matching traction is `T=H(a*v+b*w)`:
```text
slow:  w_inner = w_outer       and T_inner=T_outer
rapid: v_inner = v_outer       and T_inner=T_outer
```
For the low->high transfer, the existing Fortran
`cowling_native_interface_transfer` enforces these relations.
A symbolic algebra check (Wolfram evaluator) reduced both transfer
traction residuals and the respective continuous component residuals
to zero. The physical conditions are consistent with Tonetto & Lugones,
*Phys. Rev. D 101, 123029 (2020)*, Sec. III,
https://arxiv.org/html/2003.01259v2 .

## Concrete Fortran upgrade

New `fortran_eos/src/compose_eos_cowling_sharp_root.f90` contains
`cowling_native_sharp_onejump_brent` (C ABI). Every root iteration
stays in Fortran and executes three native DVODE segments in order:

1. high-density phase: regular center -> interior matching radius;
2. low-density phase: surface DeltaP=0 -> Maxwell interface;
3. native slow or rapid transfer, then high-density phase:
   interface -> matching radius.

The normalized two-shooting determinant is polished by a
bracket-preserving Brent-Dekker method with inverse interpolation and
bisection fallback. A collapsed sign-changing bracket is **not admitted**
if `abs(D)>1e-5`. Invalid frequency/bounds, interface geometry,
nonpositive jump, nonfinite arrays and same-sign endpoints return
typed failure with zero certified root output.

The original DVODE numerical integrator, TOV equations, Maxwell
transfer and the legacy smooth Fortran Brent implementation are
unchanged. Only a tiny `cowling_native_source_callback_is_active()`
exposure was added to the original Fortran module so this new
PCHIP-only path can reject active source callbacks **before they are
silently discarded** on a phase setup.

Python `stellar_sharp_native.NativeSharpCowlingDVODE.solve_mode` supports
`root_backend="fortran_phase_pchip"` as an **explicit opt-in**.
`python_brentq` is retained as the default reference for backwards
compatibility, and hybrid/direct-source modes **must not** be silently
relabelled as PCHIP. Exactly one resolved Maxwell interface and no
knot-waypoint special path are required. The final eigenfunction
reconstruction still invokes existing Fortran native shooting with
Python orchestration and checks join, surface, nodes and Maxwell
traction independently.

**Important:** original EOS-based `hybrid_cs2` and
`direct_source` obtain the actual source `c_s²` via the
existing Python callback at every RHS call. This source-exact
path is NOT ported to Fortran in this slice.
Likewise, EOS extraction, one-sided PCHIP coefficient construction,
and final profile presentation still use the legacy Python frontend.

## Original EOS and numerical equivalence

Immutable source ZIP:
`data/eos_missing_data_recalc_20261119/artifacts/EOS_ALL_REQUESTED_SOURCE_LOCKED_NESTED_20261119.zip`
SHA-256:
`299de43bb4f0a87ebfc4f1139ab3c8b51c344937ec3637cdff644e9ca8ab2c04`.
Internal `SPG_M2/maxwell/REPRESENTATION.json` contract SHA-256:
`59b7cb4671347c33c25133e6497e8be346b4e5b4f871eece9bcaf59155fd8765`.

Native TOV star at `n_c=0.30385128926615257 fm^-3`:
M = **1.3999999857 Msun**, R = **14.4204196237 km**.
The nominal `samples=8193` PCHIP request, including actual source
breakpoints and near-surface points, generated phase tables of
**16,055 and 2,194** entries, respectively.
Fortran uses exactly those two one-sided arrays.

| Mode | Conversion | Original PCHIP + SciPy root (Hz) | Native Fortran root (Hz) | W/V nodes |
|---|---|---:|---:|---:|
| f | slow | 1842.341048875 | 1842.341048875 | 0/0 |
| f | rapid | 1842.311838951 | 1842.311838951 | 0/0 |
| p1 | slow | 5147.989194978 | 5147.989194978 | 1/1 |
| p1 | rapid | 5143.704240601 | 5143.704240601 | 1/1 |

Within the specific comparison runs, the reported frequencies
are numerically identical at stored double precision.
Native maximum join residual among these cases:
`1.586e-6` (below existing `1e-5` threshold).
Native interface traction invariants also passed.

This is **PCHIP-to-PCHIP numerical equivalence**, not full
microphysical source equivalence. The high-precision exact-source
DOP853 M2 rapid p1 reference is 5143.719721563 Hz, about 0.01548 Hz
above the native PCHIP frequency. The M3 CSS slow PCHIP case failed
the existing physical root/join acceptance with 2,049--8,193
samples. It remains explicitly **NUMERICALLY_UNRESOLVED**; this
implementation does not override its failed scientific gate.

Raw scoped result:
`data/cowling_sharp_native_root_20261120/PCHIP_ROOT_BENCHMARK.json`.

## Execution and performance measurements

For each of the four original source cases above, one solver and its
one-sided coefficient tables were reused, and seven warm full
profiled root calculations were timed per backend:

| Case | Old Python dispatch (median s) | Fortran native root (median s) | Speedup |
|---|---:|---:|---:|
| slow f | 0.04509 | 0.03832 | 1.18x |
| rapid f | 0.04350 | 0.03735 | 1.16x |
| slow p1 | 0.09083 | 0.08431 | 1.08x |
| rapid p1 | 0.10830 | 0.09978 | 1.09x |

The identical **source TOV and PCHIP setup is excluded from these
warm timings**, so this is not a full end-to-end speedup claim.
For reference, TOV setup measured 0.389 s and coefficient
preparation 1.205 s in this run. BLAS and OpenMP threads were
limited to one.

Fortran-only numerical execution also has a reproducible exported
input contract (18,249 one-sided PCHIP records + metadata):
`ORIGINAL_SOURCE_PCHIP_ONEJUMP.dat.gz`.
Its gzip-compressed size is about 1.8 MB versus about 4.7 MB
uncompressed, with gzip -n used for deterministic headers.
Uncompressed SHA-256:
`796c9f5fc92cee7dca11381db71b5fcbc59fd0e90e3c2ddec581c3523d7bb1fd`.
Compressed SHA-256:
`d4693cd9a986b54bc3a228430e41a4b0dce29ab7cd518c4b5111aa4c6d0391c1`.

The pure Fortran standalone binary calculates **all eight**
Adams/BDF slow/rapid f/p1 roots from the frozen input without
Python. A final independent rebuild under project-module
`-fcheck=all` reproduced the same result bytes. Measured
single-process Fortran-only maximum RSS was **7,664 KB** in
that rebuild, including the numeric input loading, with
~0.90 s elapsed. This excludes source preparation and
should not be compared as identical functionality to the
Python frontend's roughly 106,468 KB maximum RSS,
which includes EOS/TOV/PCHIP preparation and its own run.

The selected-case standalone executable can execute one slow/rapid
f/p1 mode in a separate OS process (Adams and BDF).
Four independent mode jobs showed:

| Worker processes | Four-case wall (s) |
|---|---:|
| 1 | 1.147294 |
| 2 | 0.571572 |
| 4 | 0.397642 |

The four modes share one original star input, but each worker reloads
its own copy. Every result record is byte-identical across the
process configurations, with SHA-256:
`5dba40ae490deeadc88d0959835fef7e08d37b20e41b0c4922fb0e35cfe8aa9c`.
This is ~2.89x measured case-level throughput gain, **not
thread-safe execution of a shared DVODE image**. The Netlib
`/DVOD01/` and `/DVOD02/` COMMON state requires worker
isolation.

## Tests and source/workflow integrity

Initial RED: **8 failed** because the phasewise Fortran root
interface did not exist. GREEN: 8/8, then 10/10 after adding
smooth native setup cache and unresolved M3 controls.
An additional failure injection demonstrated **one direct C ABI
false substitution**: if the source-exact callback was active,
the new PCHIP root silently dropped the pointer and returned
success. Added the Fortran callback-state gate; targeted
red->green is recorded, with the final targeted test passing.

Fresh integrated regression across C7/C8/C9/C10, new phase root,
native TOV, smooth + sharp Cowling, hybrid, exact-source callback,
slow/rapid Maxwell and evidence contracts:
**83 passed, 0 failed in 134.48 seconds**.
Log: `data/cowling_sharp_native_root_20261120/REGRESSION_SUITE.log`.

Standalone source-locked Fortran numerical test:
`fortran_eos/tests/cowling_sharp_onejump_root_test.f90`.
Rebuild tool:
`scripts/build_sharp_native_onejump_test_20261120.sh`.
Project Fortran is compiled with `gfortran -O0 -g -fcheck=all`.
The frozen Netlib Fortran-77 ODEPACK routines are compiled
separately with `-O2 -std=legacy`, avoiding the known
DGEFA old-form-array bounds instrumentation false positive.
**This does not certify memory safety of all vendor code.**

## Reproduce independently

From this branch and the locked source bundle:
```sh
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1
python3 scripts/benchmark_sharp_native_onejump_root_20261120.py
python3 scripts/export_sharp_onejump_pchip_fortran_20261120.py
gzip -dc data/cowling_sharp_native_root_20261120/ORIGINAL_SOURCE_PCHIP_ONEJUMP.dat.gz \
  > /tmp/cowling_onejump_reference.dat
scripts/build_sharp_native_onejump_test_20261120.sh \
  /tmp/ns_cowling_sharp_native/sharp_onejump_test
/tmp/ns_cowling_sharp_native/sharp_onejump_test /tmp/cowling_onejump_reference.dat
scripts/bench_sharp_native_onejump_parallel_20261120.sh \
  /tmp/ns_cowling_sharp_native/sharp_onejump_test \
  /tmp/cowling_onejump_reference.dat 4 \
  /tmp/cowling_onejump_parallel_output
```
The commands above rebuild/execute only a scoped numeric test. The exact
same host timing and scheduler effects are not guaranteed on another run.
For each reproduction choose a new benchmark output directory.

## Disposition and next phase

**PASS_SCOPED_ONEJUMP_PCHIP_NATIVE_ROOT** for the tested M2
slow/rapid f/p1 and independent Fortran replay with fixed EOS identity.
No acceptance for M3 slow PCHIP, source-exact equivalence, arbitrary
multi-jump EOS, global branch continuity, or full release is claimed.

Next meaningful Fortran-first coding slice:
replace Python per-RHS source callback in `hybrid_cs2` by a native
source-derivative evaluator or fixed phase-local source representation
with rigorously justified knot handling, while retaining the
original source `c_s²` and one-sided discontinuity semantics.
Only then attempt to extend native root integration to the M3 CSS
slow branch; run exact-source DOP853 and node/physical transport
comparisons before promoting any new modes.

References:
- Tonetto & Lugones (2020), *Phys. Rev. D* 101, 123029,
  https://arxiv.org/html/2003.01259v2
- Netlib VODE implementation and COMMON/restart semantics,
  https://netlib.org/ode/vode.f
- Netlib safeguarded root solver, https://netlib.org/fmm/zeroin.f
