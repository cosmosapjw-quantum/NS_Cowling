# Source Provenance Manifest (`compose_eos` 3.3.1)

- Reference standard: CompOSE Reference Manual v3.01 (July 2022, Sec. 4.2) & Quick Guide for Users (2023)
- Gibbs `ChargeResolvedPhase` required provenance fields: `mu_b`, `mu_q`, `pressure`, `label`, `source_id`, `source_file_path`, `source_sha256`, `model_id`, `phase_identity`, `mu_b_units`, `mu_q_units`, `pressure_units`, `normalization_convention`, `interpolation_method`, `interpolation_order`, `lepton_contributions`, `_spline`

## Tree: `eos_generator/eos_data` (20 models)

| Model | $n$ | $I_l$ | $T$ | $\max|Q_5|$ | Compo | Beta-Eq | C3 | C4 | C10 |
|---|---:|---:|---:|---:|---|---|---|---|---|
| `APR_APR` | 172 | 1 | 0.0 | 0.00e+00 | `False` | `True` | `NOT_AVAILABLE:None` | `NOT_AVAILABLE:None` | `False` |
| `BBB_BHF2` | 84 | 1 | 0.0 | 0.00e+00 | `False` | `True` | `NOT_AVAILABLE:None` | `NOT_AVAILABLE:None` | `True` |
| `GM_GM1` | 349 | 1 | 0.0 | 0.00e+00 | `True` | `False` | `CHECKED:False` | `CHECKED:False` | `True` |
| `GPPVA_DD2` | — | — | — | — | — | — | — | — | `REJECTED: CompOSE nuclear group abundance implies Yh=Abar*YI > 1; check YI semantics` |
| `GPPVA_DD2_matched` | — | — | — | — | — | — | — | — | `REJECTED: CompOSE nuclear group abundance implies Yh=Abar*YI > 1; check YI semantics` |
| `GPPVA_DD2_reexpressed` | — | — | — | — | — | — | — | — | `REJECTED: CompOSE nuclear group abundance implies Yh=Abar*YI > 1; check YI semantics` |
| `GPPVA_DDME2` | — | — | — | — | — | — | — | — | `REJECTED: CompOSE nuclear group abundance implies Yh=Abar*YI > 1; check YI semantics` |
| `GPPVA_FSU2` | — | — | — | — | — | — | — | — | `REJECTED: CompOSE nuclear group abundance implies Yh=Abar*YI > 1; check YI semantics` |
| `GPPVA_FSU2H` | — | — | — | — | — | — | — | — | `REJECTED: CompOSE nuclear group abundance implies Yh=Abar*YI > 1; check YI semantics` |
| `GPPVA_FSU2R` | — | — | — | — | — | — | — | — | `REJECTED: CompOSE nuclear group abundance implies Yh=Abar*YI > 1; check YI semantics` |
| `GPPVA_NL3wrL55` | — | — | — | — | — | — | — | — | `REJECTED: CompOSE nuclear group abundance implies Yh=Abar*YI > 1; check YI semantics` |
| `GPPVA_TM1e` | — | — | — | — | — | — | — | — | `REJECTED: CompOSE nuclear group abundance implies Yh=Abar*YI > 1; check YI semantics` |
| `GPPVA_TW` | — | — | — | — | — | — | — | — | `REJECTED: CompOSE nuclear group abundance implies Yh=Abar*YI > 1; check YI semantics` |
| `RG_KDE0v1` | — | — | — | — | — | — | — | — | `REJECTED: thermo slice index lies outside eos.t/eos.yq grid index ranges` |
| `RG_SLY4` | — | — | — | — | — | — | — | — | `REJECTED: thermo slice index lies outside eos.t/eos.yq grid index ranges` |
| `SPG_M1` | — | — | — | — | — | — | — | — | `REJECTED: CompOSE nuclear group abundance implies Yh=Abar*YI > 1; check YI semantics` |
| `SPG_M2` | — | — | — | — | — | — | — | — | `REJECTED: CompOSE nuclear group abundance implies Yh=Abar*YI > 1; check YI semantics` |
| `SPG_M3` | — | — | — | — | — | — | — | — | `REJECTED: CompOSE nuclear group abundance implies Yh=Abar*YI > 1; check YI semantics` |
| `SPG_M4` | — | — | — | — | — | — | — | — | `REJECTED: CompOSE nuclear group abundance implies Yh=Abar*YI > 1; check YI semantics` |
| `SPG_M5` | — | — | — | — | — | — | — | — | `REJECTED: CompOSE nuclear group abundance implies Yh=Abar*YI > 1; check YI semantics` |

## Tree: `eos_generator/eos_data_normalized` (18 models)

| Model | $n$ | $I_l$ | $T$ | $\max|Q_5|$ | Compo | Beta-Eq | C3 | C4 | C10 |
|---|---:|---:|---:|---:|---|---|---|---|---|
| `APR_APR` | 172 | 1 | 0.0 | 0.00e+00 | `False` | `True` | `NOT_AVAILABLE:None` | `NOT_AVAILABLE:None` | `False` |
| `BBB_BHF2` | 84 | 1 | 0.0 | 0.00e+00 | `False` | `True` | `NOT_AVAILABLE:None` | `NOT_AVAILABLE:None` | `True` |
| `GM_GM1` | 349 | 1 | 0.0 | 0.00e+00 | `False` | `True` | `NOT_AVAILABLE:None` | `NOT_AVAILABLE:None` | `True` |
| `GPPVA_DD2` | 283 | 1 | 0.0 | 0.00e+00 | `True` | `True` | `CHECKED:True` | `CHECKED:True` | `True` |
| `GPPVA_DDME2` | 288 | 1 | 0.0 | 0.00e+00 | `True` | `True` | `CHECKED:True` | `CHECKED:True` | `True` |
| `GPPVA_FSU2` | 269 | 1 | 0.0 | 0.00e+00 | `True` | `True` | `CHECKED:True` | `CHECKED:True` | `True` |
| `GPPVA_FSU2H` | 300 | 1 | 0.0 | 0.00e+00 | `True` | `True` | `CHECKED:True` | `CHECKED:True` | `True` |
| `GPPVA_FSU2R` | 296 | 1 | 0.0 | 0.00e+00 | `True` | `True` | `CHECKED:True` | `CHECKED:True` | `True` |
| `GPPVA_NL3wrL55` | 298 | 1 | 0.0 | 0.00e+00 | `True` | `True` | `CHECKED:True` | `CHECKED:True` | `True` |
| `GPPVA_TM1e` | 303 | 1 | 0.0 | 0.00e+00 | `True` | `True` | `CHECKED:True` | `CHECKED:True` | `True` |
| `GPPVA_TW` | 289 | 1 | 0.0 | 0.00e+00 | `True` | `True` | `CHECKED:True` | `CHECKED:True` | `True` |
| `RG_KDE0v1` | 1240 | 1 | 0.0 | 6.03e-04 | `True` | `False` | `CHECKED:True` | `CHECKED:True` | `False` |
| `RG_SLY4` | 197 | 1 | 0.0 | 1.38e-03 | `True` | `False` | `CHECKED:True` | `CHECKED:True` | `True` |
| `SPG_M1` | 365 | 1 | 0.0 | 0.00e+00 | `True` | `True` | `CHECKED:True` | `CHECKED:True` | `True` |
| `SPG_M2` | 377 | 1 | 0.0 | 0.00e+00 | `True` | `True` | `CHECKED:True` | `CHECKED:True` | `True` |
| `SPG_M3` | 366 | 1 | 0.0 | 0.00e+00 | `True` | `True` | `CHECKED:True` | `CHECKED:True` | `True` |
| `SPG_M4` | 375 | 1 | 0.0 | 0.00e+00 | `True` | `True` | `CHECKED:True` | `CHECKED:True` | `True` |
| `SPG_M5` | 363 | 1 | 0.0 | 0.00e+00 | `True` | `True` | `CHECKED:True` | `CHECKED:True` | `True` |

