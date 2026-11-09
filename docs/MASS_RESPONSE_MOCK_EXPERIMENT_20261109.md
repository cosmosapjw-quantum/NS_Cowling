# Fixed-mass response theorem and thermodynamic-path mock experiment

**Status:** deterministic mock, numerically checked in the analysis environment. **All values are synthetic**, not CompOSE inputs or physical neutron-star predictions.

## Purpose and model

The experiment checks (i) the fixed-gravitational-mass derivative identity, (ii) generic fold conditioning and a possible mode-stationarity cancellation, and (iii) cold-barotrope endpoint and causality preservation along a continuous q=n/μ path.

The synthetic TOV surface is
\[
M(x,\lambda)=M_{\max}(\lambda)-\frac{A(\lambda)}{2}[x-x_t(\lambda)]^2,\quad
x=\ln(P_c/P_\star).
\]
Three core labels, thirteen crust labels and six λ samples are crossed with eight target masses. Only the lower-x inverse root is admitted as the rising branch. Synthetic f and p1 functions have direct λ response and x-dependent response. A second profile makes ωx vanish linearly at the fold.

The thermodynamic-path toy uses
\[
\mu/\mu_1=1+0.2t,\quad
q_\lambda(t)=0.8+0.3t+0.1\lambda\,t(1-t)(t-c),
\]
where c is chosen so the weighted pressure perturbation integrates to zero. The perturbation vanishes at both endpoints. With n=μq, dP/dμ=μq and ε=μ²q−P, all six slices share pressure and energy endpoint values.

## Fixed-mass derivative check

At masses 1.0, 1.33 and 1.8 Msun where the target existed, f and p1 derivative from equation (1) was compared with central finite differences in λ, re-solving x for the same M at each endpoint.
- Checks: 234.
- Maximum relative difference: \(4.19\times10^{-7}\).

This verifies internal consistency of the synthetic implementation and chain rule, not a derivative of the physical Cowling residual.

## Fold scaling

A log-log fit versus Mmax−M gave:
| Synthetic response profile | Fitted slope | Interpretation |
|---|---:|---|
| Generic nonzero ωx | -0.595 | Finite-range approach to the -1/2 fold law |
| ωx=O(x−xt) | -0.00524 | Leading fold divergence cancels |

The theorem predicts generic conditioning singularity, not an unconditional divergence of every eigenfrequency.

## Thermodynamic-path check

Six λ slices and 401 points per slice gave:
- qmin = 0.8;
- cs² in [0.34546, 0.42220];
- endpoint pressure spread \(2.78\times10^{-17}\) in dimensionless units;
- endpoint energy spread at roundoff level.

The path remained positive and monotone, with 0<cs²<1 at all sampled points. These are sampled checks; the continuous-function claim follows from the lemma in FIXED_MASS_RESPONSE_THEORY_20261109.md.

## Crossed mock cohort and coverage

The toy cohort has 39 fixed pairs, six method samples and eight masses. Complete-pair coverage is 39/39 through 1.8 Msun and 33/39 at 2.0 Msun. No missing method is imputed.

| M/Msun | Complete pairs | f median [%] | f p90 [%] | f max [%] | p1 median [%] | p1 p90 [%] | p1 max [%] |
|---:|---:|---:|---:|---:|---:|---:|---:|
| 0.8 | 39 | 0.9495 | 1.2516 | 1.3473 | 1.5381 | 1.8937 | 2.0417 |
| 1.0 | 39 | 0.9006 | 1.1650 | 1.2543 | 1.4686 | 1.7762 | 1.9146 |
| 1.2 | 39 | 0.8360 | 1.0627 | 1.1454 | 1.3885 | 1.6405 | 1.7690 |
| 1.33 | 39 | 0.7830 | 0.9830 | 1.0614 | 1.3156 | 1.5370 | 1.6589 |
| 1.4 | 39 | 0.7419 | 0.9338 | 1.0098 | 1.2740 | 1.4827 | 1.5929 |
| 1.6 | 39 | 0.6481 | 0.7533 | 0.8214 | 1.1418 | 1.3110 | 1.3591 |
| 1.8 | 39 | 0.5136 | 0.6093 | 0.6406 | 0.9205 | 1.1096 | 1.1742 |
| 2.0 | 33 | 0.4031 | 1.5652 | 2.4829 | 0.7493 | 1.4853 | 2.4612 |

At 2.0 Msun coverage drops from 39 to 33 while p90/max widen and medians fall. A median-only curve conceals the combination of fold conditioning, pair exclusion and heterogeneous response. The six missing pairs are the toy models whose synthetic Mmax falls below 2 Msun in at least one method, not random missingness.

## Reproduction

Generator: research_artifacts/mass_response/mock_mass_response.py  
Tests: python_frontend/tests/test_fixed_mass_response.py

The generator writes mock_manifest.json, paired_mock_ranges.csv, mass_coverage_summary.csv, fixed_mass_derivative_checks.csv, turning_point_scaling.csv, thermodynamic_q_path_summary.json and thermodynamic_q_path_samples.csv. The manifest status is MOCK_ONLY_NOT_A_PHYSICAL_EOS_OR_SPECTRUM. GitHub Actions is configured to run the tests and upload outputs with the native pilot artifact.
