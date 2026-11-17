# EOS Generator Scope and Explicit Limitations (`compose_eos` 3.3.1)

## Purpose

The `compose_eos` Python package (v3.3.1) and paired native Fortran executable (`compose_eos_fortran`) form a cold, one-dimensional CompOSE EOS reader, thermodynamic auditor, table reconstructor, crust-core matcher, phase-transition builder (`maxwell`, `css`, `gibbs`), and TOV/Love/Cowling helper. They are not a universal multi-dimensional microscopic EOS generator.

A successful parse, reconstruction, match, root search, or serialization check does not establish that the underlying nuclear model is physically valid. Source admission and physical validation remain separate, fail-closed requirements.

---

## CompOSE Specification & Critical Semantic Rules (Phase 2)

Per the **CompOSE Reference Manual v3.01 (July 2022, Section 4.2)** and the **CompOSE Quick Guide for Users (2023)**:

- Line 1 of `eos.thermo` stores the reference nucleon masses $(m_n, m_p)$ in MeV and the lepton flag $I_l \in \{0, 1\}$.
- Each thermodynamic row $(i_T, i_{n_B}, i_{Y_q})$ stores 7 scaled dimensionless/specific quantities $Q_1 \dots Q_7$ plus $N_{\text{add}}$ optional quantities:
  $$Q_1 = \frac{P}{n_B}\;[\text{MeV}],\qquad Q_2 = \frac{s}{n_B}\;[1],\qquad Q_3 = \frac{\mu_B}{m_n} - 1,\qquad Q_4 = \frac{\mu_Q}{m_n},$$
  $$Q_5 = \frac{\mu_L}{m_n},\qquad Q_6 = \frac{F}{n_B m_n} - 1,\qquad Q_7 = \frac{\varepsilon}{n_B m_n} - 1,$$
  so at $T=0$: $Q_6 = Q_7$, $P = n_B Q_1\;[\text{MeV fm}^{-3}]$, $\mu_B = m_n (Q_3 + 1)\;[\text{MeV}]$, $\mu_Q = m_n Q_4$, $\mu_L = m_n Q_5$, and $\varepsilon = n_B m_n (Q_7 + 1)\;[\text{MeV fm}^{-3}]$.
- **$I_l$ (lepton flag):** $I_l = 1$ means leptons ($e^-$, $\mu^-$) are included in $Q_1 \dots Q_7$; $I_l = 0$ means leptons are excluded (hadronic/baryonic sector only).
- **$Y_q$ (`eos.yq`):** Hadronic net electric charge per baryon. Under global charge neutrality with leptons, $Y_q = Y_e + Y_\mu$ (net $e^- + \mu^-$ fraction per baryon).
- **Normalization masses $(m_n, m_p)$:** In CompOSE general-purpose tables these act as unit-scaling reference constants (for example, $m_n = 938.9$ or $939.565\;\text{MeV}$). In `compose_eos` Gibbs output, the declared header pair `(939.0, 939.0) MeV` is explicitly labelled as a **CompOSE normalization reference**, never as the physical constituent masses of both input phases.

### Critical $\beta$-Equilibrium Rule and Provider Exceptions

**Do not treat $Q_5 = 0$ alone or $I_l = 1$ alone as proof of cold $\beta$-equilibrium.**

Cold $\beta$-equilibrium inference in `compose_eos.io.read_compose` (`beta_equilibrated`) requires **all** of:
1. Cold 1D table: $T = 0\;\text{MeV}$ (single nonzero $T>0$ slices raise `NotImplementedError`).
2. Lepton inclusion: $I_l = 1$ (`Il == 1`).
3. Vanishing lepton chemical potential: $\max |Q_5| < 10^{-8}$.
4. When `eos.compo` is present (`composition_available = True`), local charge neutrality holds across all rows:
   $$\max \left| Y_p + \frac{Z_{\text{av}}}{A_{\text{av}}} Y_h - Y_e - Y_\mu \right| \le 2\times 10^{-4}.$$

Observed repository cohort exceptions requiring explicit source-level handling:
- **`RG_KDE0v1`, `RG_SLY4` (`eos_data/`):** Declare $I_l = 1$ and $T = 0$, yet carry nonzero $\max |Q_5| = 6.03\times 10^{-4}$ (`RG_KDE0v1`) and $1.38\times 10^{-3}$ (`RG_SLY4`) $\implies$ `beta_equilibrated = False`. Raw `eos.t` uses bounds `0..0` while `eos.thermo` uses `iT = 1` (repaired in `eos_data_normalized/`).
- **Raw `SPG_M1..M5` and all `GPPVA_*` (`eos_data/`):** Store $Y_h$ or $Y_q$ in slot 4 (`Y_I` column) and omit free neutrons ($Y_n = 0$), breaking CompOSE C3/C4 closure until normalized via `normalize_compose_composition.py` or native Fortran `reconstruct-compo-native` (`eos_data_normalized/` achieves $I_l=1, \max|Q_5|=0, \beta\text{-equilibrated}=\text{True}$, $\max|r_{C3}|\le 1.11\times 10^{-16}$, $\max|r_{C4}|\le 5.00\times 10^{-6}$).
- **Raw `GM_GM1` (`eos_data/`):** Declares $I_l = 1$ and $\max|Q_5|=0$, yet row 49 violates charge neutrality ($\max|r_{C4}| = 2.45\times 10^{-1}$). `GM_GM1` is non-neutral (`beta_equilibrated = False`), and its corrupted `eos.compo` is omitted in `eos_data_normalized/GM_GM1` (`composition_available = False`).
- **`APR_APR`, `BBB_BHF2`:** Cold lepton-bearing tables ($I_l=1, \max|Q_5|=0$) with no `eos.compo` (`composition_available = False`, C3/C4 status `NOT_AVAILABLE`).

---

## Supported Input/Output Contract & Readback Verification (Phase 6)

The current `EOSData`, CompOSE I/O, and staged writer paths support:

- One selected one-dimensional cold ($T=0$) slice in baryon density $n_B$, with a single temperature point (`eos.t`) and a single charge-fraction grid point (`eos.yq`). Multi-$T$ or multi-$Y_q$ tables and nonzero-$T$ single slices are rejected (`NotImplementedError` / `ValueError`).
- CompOSE lepton flag $I_l \in \{0, 1\}$, `Nadd = 0` or `Nadd = 1`.
- Four tracked particle IDs (`10` neutron, `11` proton, `0` electron, `1` muon) plus at most one averaged nuclear group (`qid = 1`, $A_{\text{av}}, Z_{\text{av}}, Y_I$ with internal $Y_h = A_{\text{av}} Y_I$).
- Missing `eos.compo` and per-row presence masks kept distinct from explicit zeros.
- Category phase identifiers (`Iphase`) preserved across matching and normalization boundaries.

Every generator/writer path (`write_compose`, `write_compose_verified`, `write_maxwell_eos`, `write_gibbs_eos`, `causal_complete_eos.py`, `run_all_390_matching_experiments.py`, CLI `rebuild`/`reconstruct`/`match`) writes to a staging directory, runs `verify_compose_readback` (and `read_maxwell_eos` / `read_gibbs_eos` where applicable) checking all numerical arrays, $Q_1 \dots Q_7$, units/mass headers, `Iphase`/`Nadd`, composition arrays/masks, `C3`/`C4` status and pass flags, and source/output SHA-256 digests before atomic promotion. Failed readbacks exit nonzero and preserve a failure receipt.

---

## Exact Causality Interval Analysis (Phase 4)

For smooth spline intervals in $x = \ln n_B$, $p(x) = \ln P(x)$, $e(x) = \ln \varepsilon(x)$:
$$c_s^2(x) = \frac{dP}{d\varepsilon} = \exp\bigl(p(x) - e(x)\bigr)\frac{p'(x)}{e'(x)}.$$
On each piecewise-cubic interval, stationary points satisfy the degree-$\le 6$ polynomial equation:
$$\bigl(p'(t) - e'(t)\bigr)\,p'(t)\,e'(t) + p''(t)\,e'(t) - e''(t)\,p'(t) = 0.$$
Both `causality_profile` and `causal_completion` isolate all real stationary roots in $(0, h_i)$ on every piecewise-cubic span, verify positivity of the quadratic derivatives $p'(t) > 0, e'(t) > 0$ via vertices and endpoints, compute exact interval $[c_{s,\min}^2, c_{s,\max}^2]$, and distinguish `true_crossing` ($c_s^2 = 1$ crossed into superluminal regime) from `tangential_contact` ($c_s^2 = 1$ touched without exceeding 1, rejected by `causal_completion`) and `none`. Uniform 16/32-point grids are emitted only as non-gating `grid_diagnostic_*` telemetry.

---

## Gibbs Mixed-Phase Thermodynamics, Validation & Provenance (Phases 3 & 5)

- **Required inputs:** Two charge-resolved `ChargeResolvedPhase` surfaces $P(\mu_B, \mu_Q)$ with all 12 provenance fields (`source_id`, `source_file_path`, `source_sha256`, `model_id`, `phase_identity`, `mu_b_units`, `mu_q_units`, `pressure_units`, `normalization_convention`, `interpolation_method`, `interpolation_order`, `lepton_contributions`) plus unit checks and explicit rejection of double-counted lepton pressure.
- **Root isolation:** Partitions fixed-$\mu_B$ $\mu_Q$-slices at the union of both bicubic spline knot vectors, solves the quadratic stationary points of $q_H - q_Q$ on every sub-span, isolates roots on monotone intervals via tight `brentq` (`xtol=1e-15`), rejects tangential pressure contacts and non-isolated equality intervals.
- **Barotropic validation (`validate_gibbs_barotrope`):** Enforced both post-construction (in `gibbs_mixed_phase`) and post-readback (in `read_gibbs_eos`) using a single declared `PchipInterpolator(mu_b, nb)` representation:
  $$\frac{dP}{d\mu_B} = n_B,\qquad P_{\text{rec}}(\mu_B) = P(\mu_0) + \int_{\mu_0}^{\mu_B} n_B(u)\,du,\qquad \varepsilon + P = \mu_B n_B,$$
  $$\frac{dn_B}{d\mu_B} > 0,\qquad g(\mu_B) = \mu_B \frac{dn_B}{d\mu_B} - n_B \ge 0 \iff 0 < c_s^2 = \frac{n_B}{\mu_B \, dn_B/d\mu_B} \le 1.$$
  Interval minima/maxima of $dn_B/d\mu_B$, $g$, and $c_s^2$ are evaluated over every PCHIP cubic interval (endpoints, quadratic vertices, and degree-5 stationary roots), not at endpoints only.
- **Scope note:** The generated Gibbs EOS (`schema: "compose-eos-gibbs-v2"`) carries `composition_status = "thermo_only"` and `(939.0, 939.0) MeV` explicitly as CompOSE normalization reference masses; species-resolved mixed-phase composition is not synthesized.

---

## CLD / Thomas-Fermi / Wigner-Seitz Scope (Phase 7)

- **Production native CLI path added (`reconstruct-compo-native`):** Implemented in `fortran_eos/src/compose_eos_wigner_seitz.f90::reconstruct_compo_closure_ws` and exposed in `compose_eos_fortran reconstruct-compo-native RAW_DIR OUT_DIR [RECEIPT_JSON]`. It ingests admitted raw CompOSE directories (`SPG_M1..M5` and all `GPPVA_*`), enforces cold $T=0$ and composition gates, reconstructs $(Y_n, Y_p, Y_h = A_{\text{av}} Y_I)$ via C3/C4 closure matching `normalize_compose_composition.py` to `0.000e+00`, writes normalized CompOSE files, performs native `read_compose` roundtrip verification, and records deterministic FNV-1a hashes.
- **Explicit physics boundary:** `reconstruct_compo_closure_ws` / `normalize_compose_composition.py` perform **row-wise algebraic baryon and charge closure** preserving provider $(A_{\text{av}}, Z_{\text{av}}, Y_e, Y_\mu)$ and `Iphase`; similarly `compose_eos_cld.f90::cld-reconstruct` is a parameterized liquid-drop estimator. Neither module solves a self-consistent spatial Thomas-Fermi density-profile boundary-value problem or variational multi-geometry nuclear pasta minimization.

---

## Remaining Out-of-Scope Extensions

1. Multi-dimensional $(T \times n_B \times Y_q)$ CompOSE tensor reader/writer.
2. Multiple simultaneous nuclear groups ($N_{\text{quad}} > 1$) or arbitrary particle species outside `{10, 11, 0, 1}`.
3. Self-consistent microscopic RMF/EDF/Thomas-Fermi spatial crust solver or variational pasta solver generating new CompOSE tables from coupling constants alone.
4. Species-resolved mixed-phase composition inside Gibbs `compose-eos-gibbs-v2` tables.
5. Interface-aware phase-discontinuity export in smooth core-crust `match --tov` (separate `run_phase_transition_tov.py` handles Maxwell transitions).
