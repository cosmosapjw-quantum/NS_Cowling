# NS_Cowling C8: Fortran DVODE lookup acceleration and process-isolated parallelism

Date: 2026-11-20 KST. Physics baseline: `d6fdd55d271f397dc7fc9164cddc416ed9ecf655`.
Branch: `research/cowling-fortran-c8-cache-process-parallel-20261120`.
Selected PhysMath stage: C8 (target-specific numerical equivalence and measured cost).
Native subagent: NOT_RUN; current runtime offered no native spawn/delegate tool.

## Strict architecture and scientific scope

New production numerical changes are **Fortran only**:
`fortran_eos/src/compose_eos_cowling_dvode.f90`,
`fortran_eos/src/compose_eos_cowling_tov_dvode.f90`.
The new native probe is also Fortran, and a Bash harness safely orchestrates
independent operating-system processes. No new Python numerical solver was added.
Existing Python remains legacy/frontend orchestration and a verification reference;
the historical Python EOS sampling, PCHIP setup and Brent root-search have **not**
yet been ported to Fortran, so total Python elimination is NOT claimed.

Preserved: equilibrium-barotropic Gamma1=Gamma_eq closure, signature (-,+,+,+),
TOV/Cowling differential equations, source EOS, sharp Maxwell slow/rapid maps,
physical source callback semantics and DVODE tolerance/step controls.

DVODE has documented internal COMMON blocks (/DVOD01/, /DVOD02/) and this
project additionally has module-level SAVEd source/state. No OpenMP thread
is allowed to enter the solver concurrently. Independent worker **processes**
may run safely because each has its own module state and DVODE COMMON.

## Optimization

For native PCHIP lookup in each of the two Fortran modules, retain one
process-local last-sample-interval integer (4 bytes on this compiler).
At the next evaluation, check that interval first, then its immediate
neighbors, and fall back to an exact monotone binary search whenever
the proposed interval does not contain the query. Initial setup resets
the hint. Endpoint clamping, cubic Hermite formula, physical sanity
tests, callback mode and phasewise restart are unchanged.

This introduces no new dynamically allocated arrays. It does not imply
lower peak RSS, and 4-process parallel work duplicates process memory.

## Direct native benchmark (synthetic coefficients, NOT a neutron star)

Fortran source:
`fortran_eos/tests/cowling_c8_native_probe.f90`.
Inputs: 8,193 smooth sample points, 2-component native DVODE
shooting system, geometric frequency omega=2*pi*f/c where
`c=299792.458 km/s`. Compiler gfortran 13.3.0 `-O2`; each job
executes an identical deterministic sequence of native frequency shots
with strict local rtol=1e-9/atol=1e-12 and internal step cap 0.005.
`OPENBLAS_NUM_THREADS=OMP_NUM_THREADS=MKL_NUM_THREADS=1`.

One worker, 45,000 native shots:
- Baseline 34,710,585 RHS evaluations, checksum 2.8647555725735863E+05,
  CPU 7.417657 s and wall 7.42 s.
- Optimized exactly the same 34,710,585 RHS and checksum,
  CPU 6.840924 s and wall 6.84 s.
This is one paired timing observation, not an empirical distribution.

### Separate 8-job benchmark, 12,000 shots per independent job

| Binary | Worker processes | Total shots | Wall (s) |
|---|---:|---:|---:|
| Frozen baseline | 1 | 96,000 | 15.801767 |
| Optimized | 1 | 96,000 | 14.688454 |
| Optimized | 2 | 96,000 | 7.466657 |
| Optimized | 4 | 96,000 | 3.728798 |

All 8 deterministic job records are **byte-identical** across
baseline P1 and optimized P1/P2/P4, with canonical result SHA-256
`4a70378cf7ecb812030dffcf01a902df0d32c5c580d9ffa3fc1bbd41fd165879`.
Optimized P1 is ~1.076x baseline P1; optimized P4 is ~3.939x
optimized P1 or ~4.238x baseline P1, for this synthetic batch only.
No cross-worker state was shared.

Fortran development and test runner:
`scripts/bench_cowling_fortran_c8_process_parallel.sh`.
It requires an existing Fortran probe executable, workers, number of
independent jobs, shot count, destination; it launches only independent
processes. It is a **bounded process-parallel benchmark harness**, NOT a
claim of production parallel speed for source-backed stellar EOSs.

### Memory and bounds

At 9,000 shots: baseline maximum RSS 6,120 KB, optimized 6,224 KB.
Single-run OS accounting is noisy; **no measured memory reduction**.
A separate `gfortran -O0 -g -fcheck=all -std=legacy` executable
completed 400 shots, and its canonical native checksum matched the
optimized production build. Legacy vendor warning output under the
production -O2 compiler was retained in `OPT_BUILD.log` (65 warnings);
a `-std=legacy` bounds build completed with no compiler output.

## Actual original-EOS scientific regression

Fresh `pytest --run-slow` on the revised Fortran backends:
**38 passed, 0 failed in 72.24 seconds**.
Modules include native TOV, native Cowling, sharp Maxwell phase
interfaces, hybrid cs² and C7 source-locked equilibrium closure:
`data/cowling_fortran_c8_20261120/REAL_STAR_REGRESSION.log`.

This is a science-parity pass within these test contracts. It is not
a production star-level timing speedup, general frequency convergence
rate, validation of frozen-composition g-modes, or a full EOS matrix.

## Exact provenance, build and reproduction

Original EOS bundle (unchanged):
`data/eos_missing_data_recalc_20261119/artifacts/EOS_ALL_REQUESTED_SOURCE_LOCKED_NESTED_20261119.zip`
SHA-256 `299de43bb4f0a87ebfc4f1139ab3c8b51c344937ec3637cdff644e9ca8ab2c04`.
Baseline Fortran Cowling core SHA-256:
`6f75a2f218b87af835837a2040ff6c9ce20df7707504209ad3a93916e7097c0b`.
Optimized Cowling core SHA-256:
`4f1e948837230e8e5a904569e398898437f9eb8843bd94514a4a7ea95d518c92`.
Optimized TOV core SHA-256:
`ed8858a14910f9e6eb7a8ea3623678c08f70065bdae197570dc851c5b77bb174`.

The original parent revision and synthetic native source allow both
binaries to be rebuilt from the same compiler flags. The reusable
`scripts/build_cowling_fortran_probe_c8.sh` compiles the unchanged
`vendor/odepack/{vode,dgefa,dgesl,dgbfa,dgbsl}.f`
with `-O2 -ffree-line-length-none -lblas`. The explicit `parent`
option extracts the exact baseline Cowling module from
`d6fdd55d`; the `current` option uses this branch.

```sh
scripts/build_cowling_fortran_probe_c8.sh \
    /tmp/ns_cowling_c8_fortran_20261120/c8_baseline parent
scripts/build_cowling_fortran_probe_c8.sh \
    /tmp/ns_cowling_c8_fortran_20261120/c8_optimized current
```

Fresh rebuild parity (job 3, 4,500 shots) completed with
identical numerical result rows and unchanged cumulative DVODE RHS
evaluation count across the two binaries. Run, for example:

```sh
scripts/bench_cowling_fortran_c8_process_parallel.sh \
    /tmp/ns_cowling_c8_fortran_20261120/c8_optimized 4 8 12000 \
    data/cowling_fortran_c8_20261120/reproduction_P4
```

Raw native output/checksums, OS timings, and regression logs reside
under `data/cowling_fortran_c8_20261120/`.

## Admission and next genuine work

Status C8: **PASS_SCOPED_NATIVE_LOOKUP_AND_PROCESS_ISOLATED_BENCHMARK**.
This is a numeric identity and measured cost conclusion for the
specified synthetic load plus science-parity regression on selected
original EOS cases. Keep process-serial DVODE; no advertised general
thread safety. Follow-up C9: independent adversarial review,
matched-accuracy end-to-end stellar timings, actual EOS/stellar
workload process scheduling, and eventually fully-native thermodynamic
EOS / mode-bracket computation without Python numerical work.

Primary reference for DVODE state and error-control semantics:
https://netlib.org/ode/vode.f
