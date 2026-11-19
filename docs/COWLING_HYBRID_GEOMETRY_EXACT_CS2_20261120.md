# NS_Cowling: hybrid geometry plus exact source sound speed

Date: 2026-11-20 KST.
Parent commit: 6253dc40eee188ea207450ae13a5d4588aa6d398.
Branch: research/cowling-hybrid-geometry-cs2-20261120.

Task: continue and improve the original native Fortran Cowling code without
altering source EOS, TOV, slow or rapid phase conversion or the default backend.
Selected scientific coding scope: implementation, regression, and numerical
accuracy/cost validation (C5/C6/C8). Native subagent NOT_RUN because no native
spawn/delegate tool is available in this host; serial fallback actually ran.

## 1. Physical/numerical issue and proposed method

The previously added direct-source backend evaluates all five background
coefficients (r, Phi, Lambda, dPhi/dr, cs2_over_c2) by Python callback on
every Fortran DVODE RHS call. This eliminates the PCHIP smoothing of sharp
internal sound-speed derivative knots, but also recomputes smooth geometry
at large Python callback cost.

Here the new OPTIONAL coefficient_backend='hybrid_cs2' uses:
- phase-local smooth PCHIP interpolation for r, Phi, Lambda, dPhi/dr;
- the original EOS directly at every RHS for cs2/c2 only;
- the existing native restart at genuine source sound-speed jumps;
- the exact same physical Maxwell interface transfer, never smoothing
  across a density jump.

The Cowling equations use metric signature (-,+,+,+), and
omega_geom=2*pi*f/C_KM_S with c explicitly retained.
Source units: n_B in fm^-3; pressure and epsilon in MeV/fm^3;
r in km, frequencies in Hz. No full GR metric perturbation is solved.

The M3 0.7 Msun source EOS has a strong within-phase cs2 discontinuity.
At one source knot, only 5e-11 in s=sqrt(h0-h) away:
original cs2=0.00263274636, whereas PCHIP cs2=0.501315096.
The relative local coefficient error is 189.415.
This is NOT a 189-fold eigenfrequency error; it demonstrates that a
PCHIP interpolation of cs2 may approximate the wrong differential operator.

## 2. Frozen two-star numerical benchmark

Input bundle:
data/eos_missing_data_recalc_20261119/artifacts/
EOS_ALL_REQUESTED_SOURCE_LOCKED_NESTED_20261119.zip
SHA256 299de43bb4f0a87ebfc4f1139ab3c8b51c344937ec3637cdff644e9ca8ab2c04.
Separate independently source-evaluating DOP853 frequencies are the
prior-stage exact-EOS reference values. Same physical contract, central
density, slow/rapid conversion, native TOV source, ODE control and root
bracket for every compared native backend. Only the Cowling coefficient
representation and its sample count changes. Threads fixed to 1.

### M3 CSS, M=0.7 Msun, slow p1 with W/V=(1,1)

Reference DOP853: 3736.208895158 Hz.

| Backend | Absolute df from DOP853 (Hz) | Wall time (s) | Result |
|:--|--:|--:|:--|
| PCHIP 8193 | N/A | 0.88 | root rejected |
| PCHIP 16385 | N/A | 1.41 | root rejected |
| PCHIP 32769 | N/A | 2.54 | root rejected |
| Direct source callback, sentinel 129 | 0.001581 | 10.17 | admitted |
| Hybrid cs2, geometry 2049 | 0.001706 | 5.56 | admitted |
| Hybrid cs2, geometry 4097 | 0.001778 | 6.61 | admitted |

At the chosen 0.005-Hz error criterion the new hybrid 2049 method
was 1.83x faster than directly evaluating all coefficients
on this particular low-mass test. Direct-source remains slightly closer
in frequency to the reference. Neither result proves 0.001-Hz accuracy.

### M2 Maxwell, M=1.4 Msun, rapid p1 with W/V=(1,1)

Reference DOP853: 5143.719721563 Hz.

| Backend | Absolute df from DOP853 (Hz) | Wall time (s) | Result |
|:--|--:|--:|:--|
| PCHIP 8193 | 0.015483 | 1.30 | admitted at 0.05 Hz only |
| PCHIP 16385 | 0.003511 | 1.74 | admitted at 0.005 Hz |
| PCHIP 32769 | 0.002711 | 2.80 | admitted at 0.005 Hz |
| Direct source callback, sentinel 129 | 0.0000206 | 6.21 | admitted |
| Hybrid cs2, geometry 2049 | 0.0000997 | 5.53 | admitted |
| Hybrid cs2, geometry 4097 | 0.0000548 | 5.14 | admitted |

At 0.005-Hz accuracy PCHIP 16385 is cheaper. At 0.001 Hz, hybrid
4097 is the fastest ACCEPTED among the six measured native variants.
The ratio between hybrid and direct source wall times was about 1.21
here, not a general speedup guarantee.

### Conditional numeric method choice (observed only)

| Error threshold | M3 CSS 0.7 slow | M2 Maxwell 1.4 rapid |
|:--|:--|:--|
| 0.05 Hz | hybrid 2049 | PCHIP 8193 |
| 0.005 Hz | hybrid 2049 | PCHIP 16385 |
| 0.001 Hz | NO ACCEPTED CANDIDATE | hybrid 4097 |

This result is explicitly NOT an empirical calibration for untested
EOSs, masses, conversion conditions or avoided crossings.
The true physical source, local pressure/energy and phase boundaries
were never replaced. Mode labels W1/V2 remain conditional.

## 3. Evidence, reproducibility and checks

- Source code: python_frontend/compose_eos/stellar_sharp_native.py.
- Real physics tests: python_frontend/tests/test_sharp_cowling_hybrid_cs2_20261120.py.
- Hash-locked benchmark evidence tests:
  python_frontend/tests/test_sharp_cowling_hybrid_evidence_20261120.py.
- Executable producer:
  scripts/benchmark_cowling_hybrid_cs2_20261120.py.
- Honest per-case method choice:
  scripts/summarize_hybrid_accuracy_20261120.py.
- Numerically measured outputs:
  data/cowling_hybrid_geometry_20261120/RESULT.json,
  DECISIONS.json, BENCHMARK_RUN.log.

Reproduction shell:
  export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 PYTHONPATH=python_frontend
  python3 scripts/benchmark_cowling_hybrid_cs2_20261120.py
  python3 scripts/summarize_hybrid_accuracy_20261120.py
  python3 -m pytest -q python_frontend/tests/test_sharp_cowling_hybrid_cs2_20261120.py python_frontend/tests/test_sharp_cowling_hybrid_evidence_20261120.py --run-slow

The physical basis for sharp slow/rapid matching is
Tonetto & Lugones (2020), Phys Rev D 101, 123029,
https://arxiv.org/abs/2003.01259 .
No error bound for all modes or publication certification is claimed.

Status at report writing: benchmark numerically checked and new real-star
tests run. Final implementation-verification status requires a fresh
relevant regression suite, observed compiler results and source checks.

## 6. Fresh implementation verification

The original native matcher binary was freshly compiled in the isolated
worktree from its source and the original GPPVA, SHEN and nucleardatapy
sources were independently SHA-checked, without Git redistribution.

Actual **post-change** runtime test executions:

- Previously established legacy suite: **172 passed, 3 skipped**
  (verification/LEGACY_REGRESSION_RUN.log, 42.80 seconds).
- Earlier source-callback, certified-mode and new hybrid tests with real
  stars: **30 passed**
  (verification/COWLING_SOURCE_FULL_RUN.log, 166.58 seconds).
- Final hybrid module after adding explicit error/callback reset regression:
  **9 passed** (verification/HYBRID_FINAL_TEST_RUN.log, 24.80 seconds).

The final nine-test run overlaps eight tests with the earlier 30-test suite;
there are therefore **203 distinct passed tests and 3 skipped** across
these logged suites, not 211 distinct tests. No test failed in these final
three executions. Native callback poison protection has an explicit regression.

The Fortran underlying physical equations and the original default
PCHIP/direct_source implementation code are not modified by this commit;
the newly added hybrid is opt-in. The optional original physical source
callback always resets when Fortran returns, including on source failure.
