# compose_eos v3.3.1 - CompOSE EOS auditing and neutron-star phase transitions

This package reads/writes cold CompOSE EOS tables, normalizes selected composition inconsistencies into a separate output tree, audits thermodynamic consistency, constructs core-crust bridges, and exposes sharp Maxwell/Gibbs phase-transition and phase-aware TOV tools.

## Input data and normalization

Raw provider tables under `eos_data/` are preserved. Rebuild normalized data separately:

```bash
python scripts/normalize_compose_composition.py --input-root eos_data --output-root eos_data_normalized
```

The normalizer records source/output SHA-256 hashes and closure residuals. Unsupported or unresolved composition rows cause that model's entire `eos.compo` to be omitted from the normalized copy; raw input is untouched. The package currently supports a cold 1-D slice and at most one averaged nuclear group per row.

Thermodynamic phase identifiers (`Iphase`) are retained independently of the optional `eos.compo` abundance file. Derivative stencils and density-gap operations respect `Iphase` even for thermo-only tables.

## Core-crust matching

Supported methods:

- `first_law_constrained_hermite_pressure_bridge`
- `css`
- `cubic_F`, `cubic_P`
- `quintic_F`, `quintic_P`

The `css` method solves for one constant sound speed squared (c_s^2) across the entire bridge. It has a finite density/energy jump at the crust--CSS interface with continuous (P,\mu_B), and continuous (n_B,P,\epsilon,\mu_B) at the CSS--core endpoint. The root solver rejects cases without a causal solution (0<c_s^2\le1) and a positive crust-side density jump. It does not fall back to two independently fitted sound speeds. This construction is thermo-only and does not infer composition from the CSS barotrope.

The Hermite, cubic and quintic methods remain continuous-density interpolants. The token `maxwell` is not a smooth-matcher method; it denotes the separate sharp local-neutrality phase-transition API.

Recreate the 390-way matching matrix with:

```bash
PYTHONPATH=compose_eos_pkg python scripts/run_all_390_matching_experiments.py \
  --eos-root eos_data_normalized --output ./matched_eos_v3_3 --allow-failed
```

The historical matrix recorded 234 PASS / 156 C10 FAIL, with failures attributed to source-core acausality in APR_APR and RG_KDE0v1. Those counts predate this CSS semantic change and must be recomputed before making claims about the updated matrix.

## Sharp Maxwell and CSS high-density phase

Two locally neutral phase barotropes:

```bash
PYTHONPATH=compose_eos_pkg python scripts/construct_maxwell_eos.py \
  --low-phase /path/to/phase_H --high-phase /path/to/phase_Q --out ./maxwell_eos
```

Phenomenological CSS high-density phase:

```bash
PYTHONPATH=compose_eos_pkg python scripts/construct_css_maxwell_eos.py \
  --hadronic eos_data_normalized/APR_APR --pressure-transition 123.020157 \
  --delta-epsilon 120 --cs2-high 0.75 --n-max 2.0 --out ./APR_CSS_Maxwell
```

The latter is a separate sharp Maxwell transition into a phenomenological CSS high-density phase. It is not the continuous CSS--core matching path above.

A Maxwell output has adjacent equal-pressure/equal-chemical-potential rows and a finite density/energy jump, with `eos.phase_transitions.json` provenance. The phase-aware TOV/Love path applies the pressure-corrected tidal jump of Takátsy and Kovács (2020). Seidov is reported only as a small-core onset criterion, not a full radial-mode stability calculation.

The core-crust `match --tov` helper remains smooth-barotrope-only. It now rejects the single-speed CSS density-jump method rather than passing the jump table to the wrong TOV solver. A phase-interface-aware exporter/solver integration remains open.

## Gibbs mixed phase

Gibbs construction requires two charge-resolved pressure surfaces (P(\mu_B,\mu_Q)) and a lepton sector:

```bash
PYTHONPATH=compose_eos_pkg python scripts/construct_gibbs_eos.py \
  --hadronic-npz hadronic_surface.npz --high-density-npz quark_surface.npz \
  --mu-b-min 900 --mu-b-max 1500 --out ./gibbs_eos
```

Each NPZ must contain `mu_b`, `mu_q`, and `pressure` arrays shaped `(len(mu_b),len(mu_q))`. The supplied 1-D beta-equilibrium CompOSE barotropes do not contain these charge-resolved surfaces and cannot be promoted to Gibbs inputs.

## Causality and reconstruction

Causality is sampled between knots within each smooth phase; derivatives are not taken across phase IDs or large density gaps. An optional constant-(c_s^2) completion produces a separate model alternative and never overwrites the source. The source/completion join is kept as an explicit derivative boundary because the sound-speed derivative may kink there.

## Install and test

```bash
python -m pip install -e ./compose_eos_pkg[test]
PYTHONPATH=compose_eos_pkg python -m pytest -q ../tests
PYTHONPATH=compose_eos_pkg python -m compileall -q compose_eos_pkg/compose_eos scripts ../tests
```
