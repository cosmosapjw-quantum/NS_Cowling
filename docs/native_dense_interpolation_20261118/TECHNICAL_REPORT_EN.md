# Dense native EOS interpolation cross-validation (2026-11-18)

## Scope and evidence identity

Repository: cosmosapjw-quantum/NS_Cowling, native-dense-interpolation branch;
baseline a0fa5fb93851a9647030f3589c24a62ddea527b1.
New directly computed model states are numerical references **within each fixed physical functional**, not the observationally true neutron-star EOS.
No raw source rows, original H branch or Maxwell phase doublets were modified.

Five native Fortran RMF variants have 4,097 density points each; four core models have 1,025 rows and APR has 1,027 with explicit Gibbs phase boundaries: **25,612** direct model states.
The numerical uniform RMF reference range is nB=0.015..1.0 fm^-3. The core source ranges respect individual contract limits, including BBB/BHF2 <=0.8 fm^-3.
The five TF crust cases add a total of 30 **independent, new N128 stationary droplet** holdouts, with geometry and initial branch fixed. Native TF optimizer failure is never promoted to a datum.

## Physical definitions, sampling and admissibility

Cold beta equilibrium: nB=dP/dmuB, epsilon=muB*nB-P, c_s^2/c^2=nB/(muB*dnB/dmuB). The reference is sampled by **direct native calculations** and the interpolation models are fit to a strictly smaller training subset.
Evaluation at held-out native chemical potentials uses numerical EOS n,P,epsilon and finite-difference native c_s^2 as separate observables. A numerically forced Euler identity is necessary for a consistent interpolant but is not proof of true physical response.
Causality is assessed only within source-supported and native-causal regions. APR 6,8,7 are split before any interpolation. APR native phase 7 becomes acausal around nB >=0.869684 fm^-3 in this numerical source; this is not attributed to interpolant error.
In the aggregate, methods with *both* coverage and sampled physical admission >=99.5% are eligible for the ranking. Low-coverage families do not win merely from their successful subset.
All comparison figures deliberately exclude linear, as requested; the original computation's linear diagnostic record remains accessible and is not used in ranking or visualized.

## Method comparison by physics family (linear excluded)

The n RMS is the holdout-count-weighted RMS of fractional baryon-number-density prediction error, displayed in ppm (1 ppm=10^-6). Coverage is the fraction of held-out points admitting finite positive n and energy; Physical is the fraction of native-causal held-out points meeting 0<c_s^2/c^2<=1.

### Uniform RMF

| Method | n RMS [ppm] | Coverage | Physical (native-causal) | Full-family ranking eligible |
|---|---:|---:|---:|---|
| hermite | 0.37925 | 100.00% | 100.00% | yes |
| exponential | 0.32058 | 100.00% | 100.00% | yes |
| rational | 0.30887 | 100.00% | 100.00% | yes |
| bernstein | 0.70618 | 100.00% | 100.00% | yes |
| power | 1.6894 | 11.19% | 11.19% | NO: insufficient admission |
| pchip | 60.722 | 100.00% | 100.00% | yes |
| akima | 32.324 | 100.00% | 100.00% | yes |
| cubic | 396 | 100.00% | 99.84% | yes |

Eligible best by this numerical holdout score: **rational** (0.30887 ppm). No universal optimum or microscopic superiority is claimed.

### Bulk core

| Method | n RMS [ppm] | Coverage | Physical (native-causal) | Full-family ranking eligible |
|---|---:|---:|---:|---|
| hermite | 11.494 | 100.00% | 100.00% | yes |
| exponential | 10.682 | 98.83% | 99.96% | NO: insufficient admission |
| rational | 10.414 | 98.83% | 99.96% | NO: insufficient admission |
| bernstein | 11.806 | 98.80% | 99.93% | NO: insufficient admission |
| power | 5.2018 | 16.68% | 16.87% | NO: insufficient admission |
| pchip | 33.728 | 100.00% | 99.99% | yes |
| akima | 61.758 | 100.00% | 99.99% | yes |
| cubic | 272.17 | 100.00% | 99.71% | yes |

Eligible best by this numerical holdout score: **hermite** (11.494 ppm). No universal optimum or microscopic superiority is claimed.

### Spatial TF N128 (droplet)

| Method | n RMS [ppm] | Coverage | Physical (native-causal) | Full-family ranking eligible |
|---|---:|---:|---:|---|
| hermite | 74.423 | 100.00% | 100.00% | yes |
| exponential | 77.043 | 100.00% | 100.00% | yes |
| rational | 78.35 | 100.00% | 100.00% | yes |
| bernstein | 70.67 | 100.00% | 100.00% | yes |
| power | 41.266 | 63.33% | 63.33% | NO: insufficient admission |
| pchip | 1297.1 | 100.00% | 100.00% | yes |
| akima | 1505.3 | 100.00% | 100.00% | yes |
| cubic | 1795 | 100.00% | 100.00% | yes |

Eligible best by this numerical holdout score: **bernstein** (70.67 ppm). No universal optimum or microscopic superiority is claimed.

## Source-by-source comparisons (15 fixed inputs)

| Family / model | Best eligible method | n RMS [ppm] | PCHIP comparison [ppm] |
|---|---|---:|---:|
| rmf/NL3_PROVIDER_PRECISION | rational | 0.34259 | 65.565 |
| rmf/PK1_CANONICAL | rational | 0.39812 | 60.725 |
| rmf/TM1_CANONICAL | rational | 0.39994 | 59.52 |
| rmf/DDMEX_PROVIDER_ROUNDED_N0 | rational | 0.14693 | 58.765 |
| rmf/DDMEX_ORIGINAL_TABLE_RHO0 | rational | 0.13996 | 58.765 |
| core/APR_APR | hermite | 0.56287 | 54.902 |
| core/BBB_BHF2 | rational | 23.12 | 36.123 |
| core/GM_GM1 | rational | 0.27442 | 21.67 |
| core/RG_KDE0v1 | rational | 0.69167 | 21.095 |
| core/RG_SLY4 | rational | 0.50959 | 21.393 |
| tf/DDME2 | power | 80.535 | 3200.6 |
| tf/DDMEX | rational | 3.1236 | 113.68 |
| tf/NL3 | hermite | 0.30398 | 32.733 |
| tf/NL3wrL55 | bernstein | 1.6726 | 105.2 |
| tf/TM1e | bernstein | 42.383 | 1366.4 |

## Spatial TF source-generation and failures

The requested 57 new TF densities yielded 26 accepted stationary N128 states directly, plus 7 recovered in an explicit second-anchor interpolation-of-initial-profile solve (same native functional and acceptance tolerance). Three of these were at existing training knots, so **30 genuinely withheld native densities** are available.
The other 24 requested densities remain rejected. The native optimizer mode 8 failures are preserved with their raw STATE.json profiles and no threshold relaxation.

| Model | Requested new densities | Initial accepted | Second-anchor recovered | Distinct holdout |
|---|---:|---:|---:|---:|
| DDME2 | 15 | 5 | 2 | 4 |
| DDMEX | 12 | 5 | 2 | 7 |
| NL3 | 8 | 8 | 0 | 8 |
| NL3wrL55 | 12 | 6 | 0 | 6 |
| TM1e | 10 | 2 | 3 | 5 |

TF benchmark on N128 droplets is not a claim about global pasta geometry, elastic crust, low-frequency Cowling spectra, or true continuous phase equilibrium. Individual TF sample sizes range from four to eight and the density-node distribution is nonuniform after native rejection.
The power-law form is particularly vulnerable to its endpoint derivative: some native intervals fail or are not admitted. Do not rank it from subset-only RMS.

## Method assessment and implementation limits

- Exponential and rational n/mu give consistently small interpolation errors on homogeneous RMF and core native branches with fully declared first-law/euler reconstruction.
- Hermite P(mu) is robust with 100% held-out coverage in the studied core models and accurate derivatives, but not a general guarantee that an arbitrary nonunified crust bridge is causal.
- PCHIP and Akima can be robust and shape-preserving in P(mu), but interpolating only pressure means its derivative can disagree with independently computed true density, making it systematically less accurate in this benchmark.
- Cubic spline may overshoot in chemical susceptibility and causality despite smooth pressure. Its errors and physical-admission failures are separately recorded.
- Bernstein slope matching can preserve positive derivative in a fixed causal interval and is competitive in TF; its fitted endpoint derivatives require source differentiation and a properly defined branch.
- APR Gibbs phase doublets must never be smoothed across phase labels, and native high-density acausality cannot be called an interpolant bug.
- BBB BHF2 native runner itself reconstructs published tabulated many-body energies: it is not a fresh solution of the microscopic BHF equations.
- These native holdout comparisons do not repeat Maxwell root discovery using new interpolation families and do not prove global ground-state matching. A dedicated paired-P(mu) root re-search is a separate obligation.

## Figures (linear excluded)

1. figures/RMF_NATIVE_HOLDOUT_NO_LINEAR.png
2. figures/CORE_NATIVE_HOLDOUT_NO_LINEAR.png
3. figures/TF_NATIVE_HOLDOUT_NO_LINEAR.png
4. figures/NATIVE_HOLDOUT_COVERAGE_NO_LINEAR.png
5. figures/NATIVE_MODEL_METHOD_HEATMAP_NO_LINEAR.png

Method coverage, physical admission, and numerical residuals are stored with per-phase metadata in data/native_dense_interpolation_20261118/cv and /tf_holdout.

## Key supporting literature

- Haensel and Potekhin, Analytical representations of unified equations of state of neutron-star matter (2004): https://www.ioffe.ru/astro/Stars/Paper/hp04.pdf
- Fortin et al., Neutron star radii and crusts: uncertainties and unified EOSs (2016): https://doi.org/10.1103/PhysRevC.94.035804
- Ferreira and Providencia, Neutron Star Properties: Quantifying the Effect of the Crust–Core Matching Procedure (2020): https://www.mdpi.com/2218-1997/6/11/220
- Verma and Mallick, EOS interpolation with first-order phase transitions (2025): https://arxiv.org/html/2501.06453

## Next reproducible action

Run a separate original crust/core source Maxwell crossing re-search with each thermodynamically admitted interpolant, holding source phase and energy reference fixed, to assess whether the actual root and density jump change. Only then infer Maxwell-success interpolation dependence.


## Independent physics and implementation verification

For each newly generated native CompOSE row, recomputed the original cold Euler residual (epsilon+P-nB*muB)/max(1,abs(epsilon+P)). The maximum across 25,612 rows was 5.7e-14. Independent finite differences for de/dn-mu and dP/dn-n*dmu/dn are recorded by phase in aggregated/NATIVE_THERMO_AUDIT.json.
APR's native phase 7 exhibited 59 numerically acausal states, with max c_s^2/c^2 1.1685; these are source-domain diagnostics and were not falsely counted as interpolation-only errors.
Independent provenance/readback/phase/holdout/figure audit: 1096 checks, failures 0.
Frontend/pipeline full suite: **1501 passed, 9 skipped, 4 subtests passed, actual exit 0**; 208 Python source files byte-identical before/after suite. The exact commands, environment and original native binary hashes are in aggregated/FULL_SUITE_RECEIPT.json.

All 5 plotted figures use only the eight declared nonlinear methods; raw linear results remain diagnostic evidence, not a plotted or ranked method.
