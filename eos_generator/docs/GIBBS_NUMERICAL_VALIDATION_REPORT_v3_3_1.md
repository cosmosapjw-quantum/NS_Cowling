# Gibbs Numerical Validation Report (`compose_eos` 3.3.1)

## Summary

- Schema: `compose-eos-gibbs-v2`
- Composition status: `thermo_only`
- Normalization masses: `(mn, mp) = (939.0, 939.0) MeV` — _CompOSE normalization reference only; does not represent constituent masses of either underlying phase_
- Single declared barotropic representation: `PchipInterpolator(mu_b_MeV, nb_fm3)` (24 intervals)
- Staged write/readback status: `PASS`

## Thermodynamic Identity & Interval Extrema Metrics

| Metric | Post-Construction | Post-Readback | Tolerance |
|---|---:|---:|---:|
| Max pressure equality residual ($\text{MeV fm}^{-3}$) | 5.684e-14 | — | — |
| Max global charge residual | 1.984e-16 | — | — |
| Max Euler rel. error $|\varepsilon+P-\mu_B n_B|/(\mu_B n_B)$ | 1.554e-16 | 2.197e-16 | 1.0e-06 |
| Max PCHIP pressure reconstruction rel. error | 1.573e-08 | 1.573e-08 | 5.0e-04 |
| RMS PCHIP pressure reconstruction rel. error | 1.023e-08 | 1.023e-08 | — |
| Interval $\min(dn_B/d\mu_B)$ ($\text{fm}^{-3}/\text{MeV}$) | 9.246520e-04 | 9.246520e-04 | $>0$ |
| Interval $\min g = \mu_B dn_B/d\mu_B - n_B$ ($\text{fm}^{-3}$) | 7.559065e-01 | 7.559065e-01 | $\ge -1.0e-11$ (scaled) |
| Interval $c_s^2$ range $[\min, \max]$ | [0.240466, 0.256815] | [0.240466, 0.256815] | $(0, 1]$ |
| Validation pass | `True` | `True` | `True` |

## Negative Controls

- `zero_root_rejected`: pass=`True` — `mu_b_values must be an increasing grid inside both phase surfaces`
- `tangential_contact_rejected`: pass=`True` — `phase charge densities are indistinguishable; mixed fraction is undefined`
- `non_isolated_interval_rejected`: pass=`True` — `Gibbs pressure equality contains a non-isolated root interval at mu_B=1050, mu_Q=[-100,-80]`
- `lepton_double_count_rejected`: pass=`True` — `had_lep already declares lepton_contributions='included'; refusing to add an external LeptonSector on top`
- `out_of_domain_rejected`: pass=`True` — `mu_b_values must be an increasing grid inside both phase surfaces`
- `acausal_barotrope_rejected`: pass=`True` — `Gibbs synthetic_acausal_check fails dP/dmu_B = n_B integral reconstruction (max_abs=7.625e+00, max_rel=7.943e-02)`
