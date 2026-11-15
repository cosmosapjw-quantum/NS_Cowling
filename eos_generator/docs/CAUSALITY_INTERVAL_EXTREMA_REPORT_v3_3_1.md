# Causality Interval Extrema Report (`compose_eos` 3.3.1)

- Acceptance basis: `exact_interval_extrema` (degree $\le 6$ stationary polynomial on every piecewise-cubic interval)
- Uniform grid sampling (`samples_per_interval=16/32`): retained strictly as non-gating `grid_diagnostic_*` telemetry.

## Normalized Source Cohort Results

| Model | $n$ | PCHIP $c_s^2$ $[\min, \max]$ | Crossing Kind | First Acausal $n_B$ ($\text{fm}^{-3}$) | $p', e' > 0$ | Stationary Roots | Root Cause |
|---|---:|---|---|---:|---|---:|---|
| `APR_APR` | 172 | [0.000453, 1.471111] | `true_crossing` | 0.8732724588425579 | `True` | 19 | `source_EOS_acausality` |
| `BBB_BHF2` | 84 | [0.000453, 0.857718] | `none` | — | `True` | 5 | `no_acausality_detected` |
| `SPG_M1` | 365 | [0.000000, 0.013538] | `none` | — | `False` | 61 | `no_acausality_detected` |
| `GPPVA_DD2` | 283 | [0.000000, 0.011313] | `none` | — | `False` | 72 | `no_acausality_detected` |
| `RG_KDE0v1` | 1240 | [0.000145, 1.107147] | `true_crossing` | 1.4257049400603699 | `True` | 4 | `source_EOS_acausality` |
| `RG_SLY4` | 197 | [0.000396, 0.954303] | `none` | — | `True` | 14 | `no_acausality_detected` |

## Causal Completions & Adversarial Interval Extrema Tests

- `APR_APR (cs2=0.95)`: `constant_sound_speed_extension` (`true_crossing`), n_match=0.873272 fm^-3, source_max_cs2=1.471111, mu_rel_res=6.568e-07
- `RG_KDE0v1 (cs2=0.80)`: unsplit phase boundary rejected (`True` — `causal completion is not phase-aware; split the source at its Iphase phase boundary first`)
- Narrow superluminal pocket: fixed-grid 16-pt max=0.998191 (<1), 32-pt max=0.999779 (<1), exact polynomial max=1.000200 (>1), crossing_kind=`true_crossing`, missed_by_uniform_grids=`True`
- Tangential $c_s^2=1$ contact: crossing_kind=`tangential_contact`, `causal_completion` rejected=`True`
