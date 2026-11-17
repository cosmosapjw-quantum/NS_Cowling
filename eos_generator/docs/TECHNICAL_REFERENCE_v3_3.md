# Technical reference - compose_eos v3.3.1

## Zero-temperature barotrope

The cold EOS uses `n_B` in $\text{fm}^{-3}$ and $P, \varepsilon$ in $\text{MeV fm}^{-3}$. Define $F = \varepsilon / n_B$. The package uses

$$P = n_B^2 \frac{dF}{dn_B}, \qquad \mu_B = F + \frac{P}{n_B}, \qquad \varepsilon = n_B F.$$

Hence $d\varepsilon/dn_B = \mu_B$ and $c_s^2 = dP/d\varepsilon$. Thermodynamic stability requires positive compressibility; relativistic causality requires $0 \le c_s^2 \le 1$.

## CompOSE Manual v3.01 Definitions ($Q_1 \dots Q_7$, $I_l$, $Y_q$)

In CompOSE Reference Manual v3.01 (Sec. 4.2): $Q_1 = P/n_B$, $Q_2 = s/n_B$, $Q_3 = \mu_B/m_n - 1$, $Q_4 = \mu_Q/m_n$, $Q_5 = \mu_L/m_n$, $Q_6 = F/(n_B m_n) - 1$, $Q_7 = \varepsilon/(n_B m_n) - 1$, with lepton flag $I_l \in \{0, 1\}$ and hadronic charge fraction $Y_q$. Cold $\beta$-equilibrium requires $T=0$, $I_l=1$, $\max|Q_5| < 10^{-8}$, and (when `eos.compo` is available) local charge neutrality $\max|Y_p + (Z_{\text{av}}/A_{\text{av}}) Y_h - Y_e - Y_\mu| \le 2\times 10^{-4}$.

## Single-speed CSS core-crust bridge

The core-crust matcher `method="css"` uses one constant $s = c_s^2 \in (0, 1]$ over the entire bridge. Let the crust-side endpoint be $(n_c, P_t, \mu_t)$ and the continuous high-density core endpoint be $(n_2, P_2, \mu_2)$. The bridge-side start density is

$$n_b = n_2 (\mu_t / \mu_2)^{1/s}.$$

The pressure integral gives

$$P_2 - P_t = n_2 \mu_2 \frac{s}{1+s} \left[1 - (\mu_t / \mu_2)^{1+1/s}\right].$$

The scalar equation is solved only on $0 < s \le 1$, and a physical solution additionally requires $n_b > n_c$. At the crust–bridge interface $(P, \mu_B)$ are continuous while $(n_B, \varepsilon)$ jump, with

$$\Delta\varepsilon = \mu_t (n_b - n_c).$$

At the bridge--core endpoint (n_B,P,\epsilon,\mu_B) are continuous. This path is thermo-only and does not infer composition. It is not yet exported through the phase-aware TOV sidecar path; the smooth `match --tov` CLI therefore rejects it. Any Seidov result refers only to the small-core onset criterion.

This is distinct from `build_maxwell_css_eos`, which constructs a sharp transition into a phenomenological CSS high-density phase.

## Sharp Maxwell transition

For locally neutral phases H and Q, coexistence is `P_H(mu_t)=P_Q(mu_t)=P_t`. The energy density and baryon density jump while `P` and `mu_B` are continuous. At coexistence, `Delta epsilon=mu_t Delta n_B`.

The small-core Seidov threshold is `Delta epsilon_crit = epsilon_t/2 + 3P_t/2`. `Delta epsilon < Delta epsilon_crit` is the immediate small-core stable side of the Seidov onset criterion. This is not a full radial-mode stability analysis.

Two input paths are available:
1. `build_maxwell_eos(low_phase, high_phase)` searches the common chemical-potential interval of two tabulated, stable, causal, locally neutral 1-D phases.
2. `build_maxwell_css_eos(hadronic, P_t, Delta epsilon, c_s2_high, ...)` constructs a phenomenological CSS high-density phase. This is a model parameterization, not a microscopic quark-matter table.

The final CompOSE table contains two adjacent rows at equal `P_t,mu_t` and different `n_B,epsilon`. `eos.phase_transitions.json` stores and validates the interface indices and physical values.

## Tidal jump at a sharp interface

In geometric units `G=c=1`, for `y=r H'/H`, the pressure-corrected jump is

$$y^+ - y^- = -\frac{4\pi r_t^3 \Delta\varepsilon}{m_t + 4\pi r_t^3 P_t}.$$

This is the finite-pressure correction by Takátsy and Kovács (2020). The interface metadata includes all quantities required to independently recompute `Delta y`.

## Gibbs construction

For each phase supply a `ChargeResolvedPhase` surface $P(\mu_B, \mu_Q)$ with pressure shape `(len(mu_b), len(mu_q))` and provenance metadata (`source_id`, `source_file_path`, `source_sha256`, `model_id`, `phase_identity`, units, `normalization_convention`, interpolation method/order, `lepton_contributions`). The code computes $n_B = \partial P / \partial \mu_B$ and charge density $q = \partial P / \partial \mu_Q$, checks convexity and lepton double-counting, solves phase-pressure equality, and imposes

$$(1 - \chi) q_H + \chi q_Q + q_l = 0.$$

The mixed baryon density is $(1-\chi)n_H + \chi n_Q$; total energy density follows $\varepsilon = \mu_B n_B + \mu_Q q_{\text{total}} - P_{\text{total}}$ with $q_{\text{total}} = 0$. A zero-temperature electron+muon lepton sector is included when not already present on the phase surfaces. The 1-D CompOSE tables in this release do not supply charge-resolved surfaces; do not pass a `read_compose()` EOS as a substitute. Output headers use `(939.0, 939.0) MeV` explicitly as CompOSE normalization reference masses and record `composition_status = "thermo_only"`.

NPZ input format: `mu_b`, `mu_q`, and `pressure` arrays, with pressure shape `(len(mu_b), len(mu_q))`.

## TOV and numerical units

`solve_tov_with_phase_transition` uses geometric variables in km and km^-2. It converts CompOSE MeV fm^-3 pressure/energy density at the interpolation boundary, integrates the high phase to `P_t`, applies the jump, and integrates the low phase to the first tabulated surface pressure. It rejects missed pressure events, non-causal local sound speed, and trapped-surface conditions. Radius is the finite-table surface radius, not an extrapolated zero-pressure surface.

## Causal completion

For constant `s=c_s^2` the continuation is

$$\mu(n) = \mu_t (n / n_t)^s,$$

$$P(n) = P_t + \frac{s}{1+s} (\varepsilon_t + P_t) \left[(n / n_t)^{1+s} - 1\right], \qquad \varepsilon = \mu n - P.$$

This preserves the first law and has `c_s^2=s` on the tail. It changes the input model, is not uniquely implied by the source table, and must be audited with all three derivative schemes and exact spline interval extrema. An output is diagnostic-only unless its package gate, exact interval causality checks, and PCHIP, cubic-spline and centered-FD checks all pass.

## Failure-cause report

Run:

```bash
PYTHONPATH=compose_eos_pkg python scripts/classify_matching_failure_causes.py \\
  --source-root eos_data_normalized --matched-root /path/to/matched_eos_390_v3_2 \\
  --output-json failure_causes.json --output-csv failure_causes.csv
```

The v3.3 report classified all 156 original full-table FAILs as source-core acausality preserved in the matched high-density tail (78 APR_APR, 78 RG_KDE0v1).

## Reproduction commands

```bash
python scripts/normalize_compose_composition.py --input-root eos_data --output-root eos_data_normalized
PYTHONPATH=compose_eos_pkg python -m pytest -q
PYTHONPATH=compose_eos_pkg python -m compileall -q compose_eos_pkg/compose_eos scripts tests
```

Exact Maxwell from tabulated phases:

```bash
PYTHONPATH=compose_eos_pkg python scripts/construct_maxwell_eos.py \\
  --low-phase /path/to/phase_H --high-phase /path/to/phase_Q --out ./maxwell_eos
```

Maxwell with a CSS high-density phase:

```bash
PYTHONPATH=compose_eos_pkg python scripts/construct_css_maxwell_eos.py \\
  --hadronic eos_data_normalized/APR_APR --pressure-transition 123.020157 \\
  --delta-epsilon 120 --cs2-high 0.75 --n-max 2.0 --out ./APR_CSS_Maxwell
```

Phase-aware TOV/Love solution:

```bash
PYTHONPATH=compose_eos_pkg python scripts/run_phase_transition_tov.py \\
  --eos ./APR_CSS_Maxwell --pc 815.642916 --out ./star.json
```

Gibbs from 2-charge surfaces:

```bash
PYTHONPATH=compose_eos_pkg python scripts/construct_gibbs_eos.py \\
  --hadronic-npz hadronic_surface.npz --high-density-npz quark_surface.npz \\
  --mu-b-min 900 --mu-b-max 1500 --out ./gibbs_eos
```

The NPZ files must contain charge-resolved data; do not synthesize them from a 1-D beta-equilibrium barotrope.


## v3.3.1 implementation note: Gibbs pressure-root isolation

For a fixed baryon chemical potential, each input `RectBivariateSpline` pressure surface restricts to a piecewise-cubic function of charge chemical potential. The pressure difference is therefore piecewise cubic on the union of both charge-axis knot vectors, and its derivative is the charge-density difference, piecewise quadratic.

The Gibbs root solver now partitions the charge domain at that union of knots, isolates every stationary point from the quadratic derivative on each span, and searches each resulting monotone subinterval. This detects narrow root pairs and tangential contacts of the declared interpolant that a fixed sign-scan can miss. A non-isolated pressure-equality interval is rejected explicitly. The guarantee is for the declared spline interpolant, not for arbitrary unsampled physics outside the source table.

Gibbs construction also records its output as a `thermo_only` barotrope (not species-resolved mixed-phase composition), enforces post-construction and post-readback barotropic consistency (`validate_gibbs_barotrope` over `PchipInterpolator(mu_b, n_B)` verifying $P_{\text{rec}}$, Euler identity, interval minima of $dn_B/d\mu_B > 0$, causality margin $g = \mu_B dn_B/d\mu_B - n_B \ge 0$, and $0 < c_s^2 \le 1$), and preserves all 12 `ChargeResolvedPhase` provenance fields in `eos.phase_transitions.json`.

## v3.3.1 implementation note: exact causality interval extrema, staged readback & native Fortran closure

- **Exact causality extrema (`causality.py`):** Uses log-spline derivatives $p(x)=\ln P$, $e(x)=\ln \varepsilon$ ($x=\ln n_B$) and isolates all real roots of $(p'-e')p'e' + p''e' - e''p' = 0$ (degree $\le 6$) on every spline interval to evaluate exact $[c_{s,\min}^2, c_{s,\max}^2]$, enforce $p' > 0, e' > 0$, and distinguish `true_crossing`, `tangential_contact`, and `none`.
- **Staged writer readback (`io.py`, `phase_transitions.py`, `cli.py`):** All CompOSE writers stage output, invoke `verify_compose_readback` (and domain readback where applicable), preserve failure receipts on nonzero exit, and promote atomically only after verification passes.
- **Native Fortran composition closure (`compose_eos_wigner_seitz.f90`):** `compose_eos_fortran reconstruct-compo-native` executes `reconstruct_compo_closure_ws` on cold SPG/GPPVA CompOSE tables with native `read_compose` readback and FNV-1a output hashes (algebraic closure only; not a spatial Thomas-Fermi solver).
