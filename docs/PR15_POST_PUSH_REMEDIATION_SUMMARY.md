# PR #15 Post-Push Remediation and Scientific Validation Summary (`compose_eos` 3.3.1)

- Repository: `cosmosapjw-quantum/NS_Cowling`
- Branch: `remediation/2026-11-15-m5a-maxwell-eos-contracts` (PR #15 targeting `research/eos390-fortran-unification-v1`, kept open and unmerged)
- Package version: `compose_eos` `3.3.1`

## 1. CompOSE Semantic Audit (Manual v3.01, July 2022, Sec. 4.2 & Quick Guide 2023)

- `Q1 = P/n_b` [MeV], `Q2 = s` [1/baryon], `Q3 = \mu_b/m_n - 1`, `Q4 = \mu_q/m_n`, `Q5 = \mu_l/m_n`, `Q6 = f/(n_b m_n) - 1`, `Q7 = \varepsilon/(n_b m_n) - 1`. At $T=0$, $f=\varepsilon \implies Q_6=Q_7$, $\mu_B \equiv m_n(1+Q_3)$.
- `I_l = 1` indicates leptons are included in `eos.thermo`; $Q_5 = 0$ indicates $\mu_l = 0$. Neither `Q_5 = 0` nor `I_l = 1` alone establishes global beta equilibrium across a table when cold neutrality fails or tables omit leptons.
- Provider exceptions preserved and enforced:
  - `RG_KDE0v1` / `RG_SLY4`: `I_l = 1`, $T=0$, $\max|Q_5| \in [6.03\times 10^{-4}, 1.38\times 10^{-3}] \implies \text{beta\_equilibrated}=\text{False}$; raw exports also carry non-normalized `eos.t`/`eos.yq` indices and are rejected by strict `read_compose()` until normalized.
  - Raw `SPG_M1..M5` / `GPPVA_*`: store heavy-cluster fraction $Y_h/Y_q$ in `eos.compo` slot 4 instead of CompOSE $Y_I = Y_h/\bar{A}$, triggering strict `Yh=Abar*YI > 1` rejection in `read_compose()` until algebraic closure normalization runs (`eos_data_normalized/*` achieves $C_3, C_4$ pass).
  - Raw `GM_GM1`: `I_l = 1`, $\max|Q_5| = 0$, but row 49 violates charge neutrality ($|\text{res}| = 2.45\times 10^{-1}$) so `beta_equilibrated=False`.
  - `APR_APR` / `BBB_BHF2`: `I_l = 1`, $\max|Q_5| = 0$, `eos.compo` absent (`C3/C4 = NOT_AVAILABLE`).

## 2. Gibbs Mixed-Phase Validation & Provenance (Phases 3 & 5)

- `ChargeResolvedPhase` carries 12 mandatory provenance fields (`source_id`, `source_file_path`, `source_sha256`, `model_id`, `phase_identity`, `mu_b_units`, `mu_q_units`, `pressure_units`, `normalization_convention`, `interpolation_method`, `interpolation_order`, `lepton_contributions`), propagated into `compose-eos-gibbs-v2` sidecars.
- Output `(mn, mp) = (939.0, 939.0) MeV` is explicitly labeled as a CompOSE normalization reference, not the constituent masses of either underlying phase; `composition_status` is explicitly `thermo_only`.
- `validate_gibbs_barotrope` validates a single declared `PchipInterpolator(mu_b_MeV, nb_fm3)` representation (24 intervals) both post-construction and post-readback using exact piecewise antiderivatives and interval extrema:
  - Pressure-equality residual: $\max = 5.684\times 10^{-14}\text{ MeV fm}^{-3}$; global charge residual $\max = 1.984\times 10^{-16}$.
  - Euler residual $|\varepsilon+P-\mu_B n_B|/(\mu_B n_B)$: `1.554e-16` post-construction, `2.197e-16` post-readback (tol `1.0e-06`).
  - PCHIP integrated pressure max rel residual: `1.573e-08` (RMS `1.023e-08`, tol `5.0e-04`).
  - Interval $\min(dn_B/d\mu_B) = 9.246520\times 10^{-4}\text{ fm}^{-3}/\text{MeV} > 0$, $\min g = \min(\mu_B dn_B/d\mu_B - n_B) = 7.559065\times 10^{-1}\text{ fm}^{-3}$, $c_s^2 \in [0.240466, 0.256815]$.
  - All 6 negative controls (`zero_root`, `tangential_contact`, `non_isolated_interval`, `lepton_double_count`, `out_of_domain`, `acausal_barotrope`) rejected.

## 3. Exact Causality Interval Extrema (Phase 4)

- `causality_profile` and `causal_completion` evaluate exact interval extrema of $c_s^2(x) = e^{p-e} p'(x)/e'(x)$ by finding real roots of the degree-$\le 6$ stationary polynomial $(p'-e')p'e' + p''e' - e''p' = 0$ on every piecewise-cubic interval, verifying $p'>0$ and $e'>0$, and classifying crossings into `true_crossing`, `tangential_contact`, or `none`.
- Fixed 16/32-point grids are non-gating `grid_diagnostic_*` fields only.
- On the adversarial narrow superluminal pocket fixture: 16-pt max = `0.998191 < 1`, 32-pt max = `0.999779 < 1`, while exact polynomial extrema detect `max_cs2 = 1.000200 > 1` (`crossing_kind = true_crossing`, `missed_by_fixed_16_and_32_grids = True`).
- Tangential $c_s^2 = 1$ contact (`crossing_kind = tangential_contact`), subluminal `none`, and unsplit multi-phase `RG_KDE0v1` (`Iphase` boundary) are all rejected by `causal_completion`.

## 4. Staged Readback & Native Fortran Crust Closure CLI (Phases 6 & 7)

- Writers stage output to a temporary directory, run native/Python readback (`Q1..Q7`, mass headers, `Iphase`, `Nadd`, `Il/iT/iYq`, composition arrays & masks, `C3/C4`, source/output SHA-256), and promote only after verification (`eos_generator/receipts/readback_receipts_v3_3_1.json`).
- Native Fortran `reconstruct_compo_closure_ws` (`fortran_eos/src/compose_eos_wigner_seitz.f90`) and CLI command `reconstruct-compo-native` (`fortran_eos/app/compose_eos_fortran.f90`) ingest all 15 raw `SPG_M1..M5` and `GPPVA_*` directories, normalize 1-based `eos.t`/`eos.yq` indices, write CompOSE tables, perform independent native `read_compose` readback, and achieve `max_diff_vs_python_normalizer = 0.0` against `normalize_compose_composition.py` across all arrays and flags.

## 5. Explicit Scientific Scope Boundaries Retained

- `reconstruct-compo-native` performs **row-wise algebraic CompOSE C3/C4 abundance closure reconstruction**; self-consistent spatial Thomas-Fermi density-profile solvers and variational pasta phase minimizations are **not** implemented.
- Cold tables only ($T=0$); multi-temperature 3D CompOSE writers remain out of scope.
- Gibbs outputs are `thermo_only` unless composition transfer is separately modeled and validated.
