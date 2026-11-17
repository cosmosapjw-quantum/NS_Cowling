# Public native generation cohort

Status: FINITE_COHORT_NATIVE_VERIFIED; source bytes unchanged during execution: True.

The public `python3 -m compose_eos.cli generate` and `read` commands were executed for all eighteen models. Core requests used 161 density seeds; APR inserted two physical Gibbs endpoints. SPG requests used 17 points with the explicitly labeled AME2408 candidate surface coefficients. Each TF request used nB={0.003,0.03,0.075} fm^-3 at N16. All thirteen crust requests additionally regenerated separate native outer and BSk22 source-fit segments.

Observed totals: **18/18 models**, **44 native tables**, **916 main-model rows**, and **3048 rows including lower segments**. Four frontend child processes ran at most concurrently, each with one OMP/BLAS thread.

| Model | Main rows | Readback tables | Wall seconds |
|---|---:|---:|---:|
| APR_APR | 163 | 1 | 0.330 |
| BBB_BHF2 | 161 | 1 | 0.333 |
| GM_GM1 | 161 | 1 | 0.482 |
| RG_KDE0v1 | 161 | 1 | 0.270 |
| RG_SLY4 | 161 | 1 | 0.319 |
| SPG_M1 | 17 | 3 | 12.353 |
| SPG_M2 | 17 | 3 | 12.439 |
| SPG_M3 | 17 | 3 | 12.526 |
| SPG_M4 | 17 | 3 | 12.523 |
| SPG_M5 | 17 | 3 | 12.161 |
| GPPVA_DD2 | 3 | 3 | 11.958 |
| GPPVA_DDME2 | 3 | 3 | 11.896 |
| GPPVA_TW | 3 | 3 | 11.946 |
| GPPVA_FSU2 | 3 | 3 | 13.019 |
| GPPVA_FSU2R | 3 | 3 | 12.930 |
| GPPVA_FSU2H | 3 | 3 | 13.020 |
| GPPVA_NL3wrL55 | 3 | 3 | 12.963 |
| GPPVA_TM1e | 3 | 3 | 16.382 |

This verifies native production dispatch, the requested finite numerical results, file integrity, composition-aware serialization and native public readback. It does not establish global TF minima, continuum-mesh convergence, exact original SPG calibration, or physical continuity between the separately exported lower segments. Native support metadata continues to govern interpolation. Final CompOSE thermodynamic values were not inputs or unconditional equality gates.

Every model has a preserved `attempt_001` directory with public command, logs, request, native outputs, dispatch receipt and `COHORT_RESULT.json`. `SUMMARY.json` records all exact source and output checksums. The runnable orchestrator is `scripts/verify_native_generator_cohort.py`.
