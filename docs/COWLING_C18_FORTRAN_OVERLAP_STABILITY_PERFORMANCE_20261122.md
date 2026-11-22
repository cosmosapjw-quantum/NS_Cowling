# NS_Cowling C18: native Fortran phasewise shape overlap, stability and throughput

> **C19 integration correction (2026-11-22).** The strict within-phase
> duplicate-radius rejection documented below passed the isolated synthetic
> C18 tests but rejected C17 real original-source eigenfunctions at a
> legitimate double-sampled internal shooting-match point. The C19 integration
> replaces it with a fail-closed one-duplicate continuity test, preserves
> the Maxwell boundary jump and reproduces all 40 original C17 stellar
> cases. See `docs/COWLING_C19_C17_C18_SOURCE_INTEGRATION_20261122.md`.
> Do **not** merge the original C18-only strict-mesh implementation into
> C17 without this correction.

Date: 2026-11-22 (Asia/Seoul). This is one bounded coding/verification
slice, not a replacement for a full EOS/TOV/Cowling physical validation.
Parent: `b956e0a21cc31623d5034ab9f0cf9153b46a1d41`
(`research/cowling-c16-source-mode-continuation-20261121`).
Worktree: `/home/cosmosapjw/ns_cowling_wt/cowling_c18_overlap_optimization_20261122`.
Branch: `research/cowling-c18-fortran-overlap-safe-perf-20261122`.
The pre-existing dirty worktrees, especially ongoing C17 work, are untouched.
Native subagent delegation: NOT_RUN (host provides no native spawn/delegate tools).

## Scope and physical semantics

Only `fortran_eos/src/compose_eos_cowling_mode_overlap.f90`
changes among production source files. The public C-ABI stays identical.
Neither the original EOS source laws nor TOV, DVODE, Cowling RHS,
Maxwell slow/rapid transport, eigenvalue root criteria, or the node gates
is changed. Python remains frontend/test/orchestration.

For each phase p, compare areal-radius W profiles at 65 quadrature
points in `u=(r-r_min)/(r_max-r_min)`, trapezoid weights.
No interpolation crosses the Maxwell interface. A repeated radius
between the two phases is permitted, but duplicate radii within a phase
are now rejected rather than silently interpolated.

The dimensionless absolute overlap is descriptive, not a covariant
oscillation inner product and not a proof of avoided-crossing absence.
In one phase write W_ip=s_ip*Wtilde_ip where s_ip=max(abs(W_ip)).
Then
- A_ip = s_ip^2 * Atilde_ip
- B_ip = s_2p^2 * Btilde_ip
- C_p  = s_1p*s_2p * Ctilde_p

For global mode maxima S_i=max_p(s_ip) and q_ip=s_ip/S_i, the total
correlation can be evaluated without forming any raw quadratic amplitude:
`abs(sum_p(q_1p*q_2p*Ctilde_p)) /
sqrt(sum_p(q_1p**2*Atilde_p)*sum_p(q_2p**2*Btilde_p))`.
The factors S_i cancel. Phase-local correlations are computed separately
from `Atilde,Btilde,Ctilde`. Input amplitudes must be finite and
nonzero within each phase. This normalizes values before subtraction
and multiplication, avoiding overflow of the raw squared norms.

For profiles up to 513 points per phase, a monotone interval cursor
replaces repeated binary searches over the 65 quadrature abscissae.
Longer profiles retain logarithmic binary lookup. No O(N) scratch
arrays are allocated and there is no SAVEd mutable state.
The overlap kernel is reentrant and can be independently processed
in multiple OpenMP workers; ODEPACK/DVODE is not reentrant and MUST
continue to run in isolated processes, not concurrently in threads.

## Reproduced defect and observed verification

Baseline original module SHA256:
`35daea1943e8e91828e868f63df9ce874b9b43d0acaa01d21973ec6fe7be99c5`.
TDD red test: `cowling_c18_overlap_stability.f90` linked with baseline
`-O0 -g -fcheck=all`: exit 14, tiny finite amplitude status=2,
IEEE overflow/underflow flags. Both final checked Fortran (`-O0 -g -fcheck=all -ffpe-trap=invalid,zero,overflow`)
and optimized Fortran (`-O2 -ffpe-trap=invalid,zero,overflow`) runs
passed `C18_OVERLAP_STABILITY_AND_PHASE_MESH_PASS`.
The test includes positive/negative normalization, W x 1e200 and 1e-200,
cross-normalization 1e200 against 1e-200, and an invalid same-phase
duplicate radius, while preserving the legitimate Maxwell duplicate.

Fresh C14-C16 selected Python frontend/native integration tests:
`PYTHONPATH=python_frontend OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1
python3 -m pytest -q
python_frontend/tests/test_c16_source_mode_continuation.py
python_frontend/tests/test_c15_tov_source_convergence.py
python_frontend/tests/test_cowling_c14_native_error_budget.py`
returned 26 passed, 6 skipped in 21.93 s (fresh final run, exit 0).
The C16 standalone Fortran scientific driver built successfully under
both checked (-O0 -g -fcheck=all) and optimized (-O2) modes.
These build checks alone do not certify an executed original-source
standalone M2/M3 driver in this C18 turn.

## Bounded microbenchmark (synthetic profiles; not full solver)

Hardware: AMD Ryzen 9 5900X (12 cores/24 threads);
Ubuntu 24.04.5 LTS, GNU Fortran 13.3.0.
`scripts/bench_c18_fortran_overlap.sh` builds from exact parent
and current source with `-O2 -fopenmp`.
Each experiment evaluates 20,000 two-phase W-overlap comparisons on
360-point synthetic profiles, 5 alternating runs for each backend
and each concurrency setting.

| Workers | C16 median seconds | C18 median seconds | Kernel-only gain |
|---|---:|---:|---:|
| 1 | 0.1387372 | 0.0684650 | 2.0264x |
| 4 | 0.0375729 | 0.0192355 | 1.9533x |

For the same thread count, the original and new implementation
produce identical displayed 15-digit aggregate checksums.
Between 1/4 workers, reduction ordering changes the least significant
digits; not a byte-reproducibility claim. The 4-worker test parallelizes
independent overlap evaluations only. The whole EOS/TOV/Cowling solver
speedup was NOT measured and is NOT claimed.

## Handoff and limitations

New checked test: `fortran_eos/tests/cowling_c18_overlap_stability.f90`.
Bench driver: `fortran_eos/tests/cowling_c18_overlap_benchmark.f90`.
Bench runner: `scripts/bench_c18_fortran_overlap.sh`.
Only Fortran owns mode-shape arithmetic; Python owns frontend plumbing.

Claim status: finite-amplitude/mesh regression IMPLEMENTATION_VERIFIED;
synthetic kernel speed NUMERICALLY_CHECKED; physical spectrum remains
the inherited C16/C15 scope, with no new full-spectrum certification.
The benchmark does not include EOS ZIP loading, TOV source-grid errors,
full profile reconstruction or root solving. The modified build is
not a published/merged release until verified Git remote action.
Next constructive task: once the parallel C17 branch is stable,
integrate and check this overlap kernel against the 20 C16
source-locked f/p1 eigenfunctions across local central densities,
with matched EOS/TOV numerical error budgets and explicit failure map.
