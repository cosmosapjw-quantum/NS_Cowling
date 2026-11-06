# Technical reference - compose_eos v3.3.0

## Zero-temperature barotrope

The cold EOS uses `n_B` in fm^-3 and `P,epsilon` in MeV fm^-3. Define `F=epsilon/n_B`. The package uses

$$P=n_B^2\\frac{dF}{dn_B},\\qquad \\mu_B=F+\\frac{P}{n_B},\\qquad \\epsilon=n_BF.$$

Hence `d epsilon/d n_B = mu_B` and `c_s^2=dP/d epsilon`. Thermodynamic stability requires positive compressibility; relativistic causality requires `0 <= c_s^2 <= 1`.

## Sharp Maxwell transition

For locally neutral phases H and Q, coexistence is `P_H(mu_t)=P_Q(mu_t)=P_t`. The energy density and baryon density jump while `P` and `mu_B` are continuous. At coexistence, `Delta epsilon=mu_t Delta n_B`.

The small-core Seidov threshold is `Delta epsilon_crit = epsilon_t/2 + 3P_t/2`. `Delta epsilon < Delta epsilon_crit` is the immediate small-core stable side of the Seidov onset criterion. This is not a full radial-mode stability analysis.

Two input paths are available:
1. `build_maxwell_eos(low_phase, high_phase)` searches the common chemical-potential interval of two tabulated, stable, causal, locally neutral 1-D phases.
2. `build_maxwell_css_eos(hadronic, P_t, Delta epsilon, c_s2_high, ...)` constructs a phenomenological CSS high-density phase. This is a model parameterization, not a microscopic quark-matter table.

The final CompOSE table contains two adjacent rows at equal `P_t,mu_t` and different `n_B,epsilon`. `eos.phase_transitions.json` stores and validates the interface indices and physical values.

## Tidal jump at a sharp interface

In geometric units `G=c=1`, for `y=r H'/H`, the pressure-corrected jump is

$$y^+-y^-=-\\frac{4\\pi r_t^3\\Delta\\epsilon}{m_t+4\\pi r_t^3P_t}.$$

This is the finite-pressure correction by Takátsy and Kovács (2020). The interface metadata includes all quantities required to independently recompute `Delta y`.

## Gibbs construction

For each phase supply a grid `P(mu_B,mu_Q)` with pressure shape `(len(mu_B),len(mu_Q))`. The code computes `n_B=partial P/partial mu_B` and charge density `q=partial P/partial mu_Q`, checks pressure convexity, solves phase-pressure equality, and imposes

$$ (1-\\chi)q_H+\\chi q_Q+q_l=0. $$

The mixed baryon density is `(1-chi)n_H+chi n_Q`; total energy density follows `epsilon=mu_B*n_B+mu_Q*q_total-P_total`, with `q_total=0`. A zero-temperature electron+muon lepton sector is included. The 1-D CompOSE tables in this release do not supply the required charge-resolved surfaces. Do not pass a `read_compose()` EOS as a substitute.

NPZ input format: `mu_b`, `mu_q`, and `pressure` arrays, with pressure shape `(len(mu_b),len(mu_q))`.

## TOV and numerical units

`solve_tov_with_phase_transition` uses geometric variables in km and km^-2. It converts CompOSE MeV fm^-3 pressure/energy density at the interpolation boundary, integrates the high phase to `P_t`, applies the jump, and integrates the low phase to the first tabulated surface pressure. It rejects missed pressure events, non-causal local sound speed, and trapped-surface conditions. Radius is the finite-table surface radius, not an extrapolated zero-pressure surface.

## Causal completion

For constant `s=c_s^2` the continuation is

$$\\mu(n)=\\mu_t(n/n_t)^s,$$

$$P(n)=P_t+\\frac{s}{1+s}(\\epsilon_t+P_t)[(n/n_t)^{1+s}-1],\\qquad \\epsilon=\\mu n-P.$$

This preserves the first law and has `c_s^2=s` on the tail. It changes the input model, is not uniquely implied by the source table, and must be audited with all three derivative schemes. An output is diagnostic-only unless its package gate and PCHIP, cubic-spline and centered-FD checks all pass.

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
