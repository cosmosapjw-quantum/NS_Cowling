# R15 native row-217 → row-218 boundary diagnostic

**Status:** `PASS_REPRESENTATIVE_BOUNDARY_READ`  
**Source cohort:** frozen projected EOS tables in `oscillation_release_v1/data/eos_projected/`  
**Representative slice:** 13 crust families × one table per family, fixed to `BBB_BHF2 × css`  
**Preflight authority:** [GitHub Actions run 42650145210](https://github.com/cosmosapjw-quantum/NS_Cowling/actions/runs/42650145210)  
**Frozen spectrum-summary SHA-256:** `ce54c111510f32086065f95356903cd6abad675adce12cab1a06b5cc82350965`

## Question and scope

The full-cohort preflight established that all 234 projected tables share the first 217 low-density source rows and that source row 218 is the first row above the operational cutoff (n_b \le 10^{-3}\,\mathrm{fm}^{-3}). This diagnostic reads the exact adjacent native rows 217 and 218 from one representative `BBB_BHF2 × css` projected table for each of the 13 crust families and compares endpoint ratios for (n_b), (P), (epsilon), and (mu_b).

This is a 13-table representative thermodynamic cross-section, not a full 234-table all-column audit. It does not change the EOS files or insert interpolation rows.

## Endpoint ratios

| Crust family | (n_b(218)/n_b(217)) | (P(218)/P(217)) | (epsilon(218)/epsilon(217)) | (Deltamu_b/mu_b) | (Deltaln P/Deltaln n_b) |
|---|---:|---:|---:|---:|---:|
| GPPVA_DD2 | 7.434 | 6.458 | 7.455 | +0.255% | 0.9298 |
| GPPVA_DDME2 | 7.434 | 6.491 | 7.455 | +0.257% | 0.9324 |
| GPPVA_FSU2 | 11.152 | 6.718 | 11.181 | +0.185% | 0.7898 |
| GPPVA_FSU2H | 7.434 | 5.659 | 7.454 | +0.218% | 0.8639 |
| GPPVA_FSU2R | 7.434 | 5.802 | 7.454 | +0.224% | 0.8764 |
| GPPVA_NL3wrL55 | 7.434 | 4.961 | 7.453 | +0.185% | 0.7983 |
| GPPVA_TM1e | 7.434 | 5.948 | 7.455 | +0.231% | 0.8888 |
| GPPVA_TW | 11.152 | 9.216 | 11.185 | +0.267% | 0.9209 |
| SPG_M1 | 9.337 | 7.212 | 9.363 | +0.237% | 0.8844 |
| SPG_M2 | 8.516 | 6.682 | 8.539 | +0.235% | 0.8868 |
| SPG_M3 | 9.777 | 7.538 | 9.805 | +0.239% | 0.8859 |
| SPG_M4 | 9.777 | 7.567 | 9.805 | +0.240% | 0.8876 |
| SPG_M5 | 9.337 | 7.152 | 9.363 | +0.235% | 0.8807 |

Across these representatives, the density ratio is 7.434–11.152, the pressure ratio is 4.961–9.216, and the energy-density ratio is 7.453–11.185. The baryon chemical potential changes by only +0.185–0.267% between the endpoint samples. The ratio ((P/n_b)_{218}/(P/n_b)_{217}) lies between 0.602 and 0.873.

## Interpretation

The data establish a wide interval between adjacent stored samples. They do **not** establish a thermodynamic discontinuity at fixed pressure or a physical interface: both density and pressure change materially between the endpoint samples, and the files provide no intermediate native rows in this interval.

The two-point secant (Delta\ln P/\Delta\ln n_b) is 0.7898–0.9324. It is only a coarse finite-interval secant across a sparse interval. It is not a local adiabatic index, (c_s^2/c^2), or a causality test. No local derivative-sensitive conclusion can be drawn from these two endpoints alone.

The modest endpoint change in (mu_b) suggests that the chemical-potential values are close across this interval, but does not guarantee a well-resolved derivative or interpolation. Since the pressure and energy-density intervals are coarse, the numerical representation used by the TOV and Cowling solvers should be tested directly.

## Next numerical object

The full-cohort solver-native representation audit has now been completed; see `R15_NATIVE_HERMITE_ROW217_218_AUDIT_20261110.md`. All 234 cubic-Hermite cells pass the exact endpoint stability/causality conditions, (\chi=c_s^2/c^2) is monotonically increasing with zero interior stationary points in every cell, and maximum endpoint First-Law energy closure residual is (2.392\\times10^{-16}). This closes the numerical interpolation concern for the current representation on this cell, not the physical low-density EOS provenance question.

The next useful diagnostic is a source/provenance reconstruction of the low-density rows: identify the original CompOSE and BSk22-derived source records, the crust-family mapping, and any projection/normalization transformations that produced rows 1–218. Keep original source hashes and derived projected-table identity separate. Do not create synthetic EOS rows or infer a physical (P=0) surface from this numerical audit.