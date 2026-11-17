# Executed candidate SPG surface calibration

Status: CANDIDATE_CALIBRATION_NOT_ORIGINAL_SCURTO_FIT. Date: 2026-11-16.
This is an executed calibration experiment for the five author-parameterized GDFM candidates. It is not an independently regenerated crust, does not claim the original Scurto surface coefficients, and is not installed as the generator default.

## Defined candidate model

For an isolated spherical nucleus with fixed A,Z, define y=Z/A, R=(3A/(4 pi n))^(1/3), and

E_candidate(A,Z;n) = A epsilon_GDFM(n,y)/n + 4 pi R^2 sigma(y) + 3 e^2 Z^2/(5R),
sigma(y)=sigma0(b+16)/(b+y^-3+(1-y)^-3).

This finite-nucleus functional is a chosen vacuum sharp-sphere construction motivated by [Scurto2024 Eqs.21-25](https://arxiv.org/html/2402.15548v1). The paper does not establish that its actual mass-fitting implementation uses exactly this functional. The construction omits curvature, shell, pairing, Coulomb exchange, electron screening and electronic binding corrections. The Coulomb term uses continuous charge Z^2, not Z(Z-1). e^2=1.43996448 MeV fm is explicit. Every nucleus is allowed to compress by minimizing over n.

The bulk evaluator is the new three-meson author-parameterized GDFM module, with common nucleon rest energy938.9MeV. The predicted binding is B_candidate=938.9A-E_candidate. It is compared with the AME2016 B/A column multiplied by A. For mass weighting, both observation and model are assigned the same reference rest energy N*939.56536+Z*938.27203; therefore their mass difference is exactly B_AME-B_candidate. This convention removes a purely additive common-mass offset from the residual. It is an explicit candidate convention, not a recovered source atomic-to-nuclear mass implementation.

The training sample is the existing measured-mass candidate selection N,Z>=8 and no estimated marker in AME2016 mass excess:2408nuclei, A=16..270 and y=0.307692..0.6. The pinned raw SHA-256 is2167f57a2a98331e4649b2dd2b658a9006ed4fba1975729ebfe52a42b4b9218a. The [Carreau2019 nine magic/semi-magic nucleus fit](https://arxiv.org/html/1902.07032) and Scurto's2408-mass score are distinct statements. This experiment deliberately trains on the entire2408candidate set; it does not infer that this was the original training protocol.

## Objective and executed variants

The primary objective is the mean squared mass residual divided by an uncertainty

u_i=0.02[N_i*939.56536+Z_i*938.27203-B_AME,i] MeV.

This implements a literal mass-based interpretation of Scurto Eq.31's stated2%mass scale under the explicit mass convention above. A separate fit uses u_i=0.02B_AME,i to measure the effect of interpreting that scale as binding energy. These objectives are not silently equated.

For each model the run executed:

- Full2408training with mass weights.
- Full2408training with binding-energy weights.
- Mass-weighted training on1927nuclei, holding out481nuclei with A modulo5=0.
- A density-grid convergence fit with81points versus161points on n in[.04,.24]fm^-3.
- Three optimizer starts for every fit, using log(sigma0) and log(b+16).
- Direct scalar minimization with the actual GDFM evaluator for12deterministically selected nuclei, independent of cached-grid interpolation.

Candidate parameter search bounds are sigma0 in[.01,10]MeV fm^-2 and b+16 in[.01,10^6]. They are computational search bounds, not claimed physical priors. Density minima first use a global search on the declared finite grid and a local quadratic interpolation. The direct minimizations check selected local brackets; neither calculation proves a global minimum outside the density support.

## Results

The primary full-sample mass-weighted fits produced:

| GDFM candidate | sigma0 / MeV fm^-2 | b | Full RMS mass residual / MeV | Held-out RMS mass residual / MeV |
| --- | ---: | ---: | ---: | ---: |
| SPG_M1 | 1.0415982 | 14.08475 | 21.29196 | 21.35963 |
| SPG_M2 | 1.0656369 | 34.85533 | 15.32784 | 15.39203 |
| SPG_M3 | 1.0611256 | 14.90681 | 22.44424 | 22.51266 |
| SPG_M4 | 1.0658723 | 28.96304 | 20.55594 | 20.62393 |
| SPG_M5 | 1.0177641 | 24.75998 | 15.21843 | 15.28477 |

The held-out column comes from the separate1927training fit, not evaluation of the all-sample fitted parameters. The residuals are substantial. A small scaled Eq.31-like score is expected because2%of a nuclear rest mass is a very large uncertainty; it is not evidence of an accurate nuclear mass model.

The primary fits have no active parameter bounds or density-edge minima. Their optimized densities lie between.143789 and.182916fm^-3 across the five models. The local optimizer Jacobian condition numbers are16.0..23.9 in log-parameter coordinates. These are local numerical diagnostics, not statistical confidence or a proof of global identifiability. All four fit variants and every start are preserved in the JSON.

Selected direct minimizations differ from the161-grid interpolated masses by at most1.07e-4MeV over60tested model/nucleus pairs. The81-to161fit comparison is recorded separately. This is evidence of numerical resolution within the chosen model and support, not evidence that its omitted finite-nucleus physics is negligible.

The JSON also records sigma(y) at y=.01,.03,.05,.1,.3,.5 under mass-weighted, binding-weighted and held-out fits. The training data have y>=.307692. Consequently the small-y surface relevant to neutron-rich clusters is an extrapolation of the fixed p=3 law. Weight sensitivity does not span the missing curvature/shell/pairing terms or alternative surface laws, and cannot be quoted as a complete uncertainty band.

## Reproduction and retained evidence

From the repository root, with a fresh output directory:

```sh
curl -fLsS https://www-nds.iaea.org/amdc/ame2016/mass16.txt -o build/candidate_surface_fit/mass16.txt
PYTHONPATH=python_frontend python3 scripts/fit_candidate_spg_surface.py \
  --ame build/candidate_surface_fit/mass16.txt \
  --output build/candidate_surface_fit/replay
PYTHONPATH=python_frontend python3 -m pytest -q \
  python_frontend/tests/test_candidate_surface_fit.py \
  python_frontend/tests/test_gdfm_uniform.py \
  python_frontend/tests/test_scurto_interface_functions.py
```

The reviewed campaign completed in31.86seconds; timings are machine-specific. The combined targeted suite passed52tests. Four invalid-A/Z tests failed before constructor input validation was repaired; their red log is retained. Independent review also found that a cache keyed only by selected source files could silently reuse bulk energies at a changed density grid. The repaired identity includes the actual density and isotope arrays, AME identity, model ID and transitive bulk dependencies. Existing output directories with content are refused to prevent mixing stale results after an interrupted rerun. Partial reports say RUNNING; the final report says COMPLETE.

Retained under data/research_generator_20261116/candidate_surface/:

- candidate_surface_fit.json: objective, assumptions, coefficients, bounds, residual metrics, all optimizer starts, holdout, numerical convergence, source/software hashes and runtime.
- candidate_surface_residuals.csv: all12040primary predictions, observed binding energies and minimizing densities.
- run.log, test_red.log, test_green.log and SHA256SUMS.txt.

Large cached bulk grids and the acquired AME source remain under build/. The code is isolated in candidate_surface_fit.py and its runner; no existing source-table extraction, legacy proxy, source-calibrated admission or crust output is silently replaced.

## Consequence for the next calculation

This experiment establishes that a labeled model-specific calibration is executable with current inputs; exact historical2408-row identity is not a blanket prerequisite to such a candidate experiment. The next constructive branch can use these explicit coefficients in a separately labeled full CLD solver and compare the resulting composition/thermodynamics with the source. Original-source reproduction still requires identifying or testing the actual finite-nucleus functional and fitting protocol, because this fit's omitted terms and large residuals are physical model limitations. GPPVA's finite-range Thomas-Fermi work remains a separate route and needs none of these candidate surface parameters.

Independent read-only review returned PASS_SCOPED for the defined candidate experiment. Fresh uncached calculations over12nuclei/model, including minimum/maximum training proton fractions, agreed with direct minimization over the full[.04,.24]fm^-3 support to1.03636e-4MeV in nuclear mass. The mechanical stationary residual was at most4.59778e-5MeV fm^-3; numerical parameter derivatives agreed with analytic envelope derivatives to9.92736e-8relative. Vacuum cell-to-nucleus conversion agreed to8.88e-16MeV/A. The verbatim review is retained as independent_review.json; it does not validate the source-fit identity or extrapolation below training y=.3076923.
