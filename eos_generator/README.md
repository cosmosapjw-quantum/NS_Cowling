# compose_eos v3.3.0 - CompOSE EOS auditing and neutron-star phase transitions

This package reads/writes cold CompOSE EOS tables, reconstructs selected composition inconsistencies into a separate normalized tree, audits thermodynamic consistency, performs smooth core-crust matching, and now exposes explicit sharp phase-transition and phase-aware TOV tools.

## Input data and normalization

Raw provider tables under `eos_data/` are preserved. Rebuild the normalized data as a separate tree:

```bash
python scripts/normalize_compose_composition.py --input-root eos_data --output-root eos_data_normalized
```

The normalizer records source/output SHA-256 hashes and closure residuals. Unsupported or unresolved composition rows cause that model's entire `eos.compo` to be omitted from the normalized copy; raw input is untouched. The current package supports a cold 1-D slice and at most one averaged nuclear group per row.

## Smooth core-crust matching

The six legacy method tokens (`maxwell`, `css`, `cubic_F`, `cubic_P`, `quintic_F`, `quintic_P`) remain smooth, single-valued barotropic bridges. The `maxwell` token is a compatibility label for a smooth near-pressure plateau; it is **not** a sharp Maxwell density-jump construction.

Recreate the 390-way matrix with:

```bash
PYTHONPATH=compose_eos_pkg python scripts/run_all_390_matching_experiments.py \\
  --eos-root eos_data_normalized --output ./matched_eos_v3_3 --allow-failed
```

The source-matrix audit found 234 PASS and 156 C10 FAIL, all due to acausality already present in APR_APR and RG_KDE0v1 source-core tails. Failed tables are diagnostic-only.

## Sharp Maxwell and CSS phase transitions

Use two locally neutral tabulated phase EOS:

```bash
PYTHONPATH=compose_eos_pkg python scripts/construct_maxwell_eos.py \\
  --low-phase /path/to/phase_H --high-phase /path/to/phase_Q --out ./maxwell_eos
```

Or construct a phenomenological CSS high-density phase from one hadronic table:

```bash
PYTHONPATH=compose_eos_pkg python scripts/construct_css_maxwell_eos.py \\
  --hadronic eos_data_normalized/APR_APR --pressure-transition 123.020157 \\
  --delta-epsilon 120 --cs2-high 0.75 --n-max 2.0 --out ./APR_CSS_Maxwell
```

The CSS phase is phenomenological, not microscopic QCD. A Maxwell output contains adjacent equal-pressure/equal-chemical-potential rows with a finite density/energy jump, plus `eos.phase_transitions.json` provenance.

Run phase-aware TOV and the corrected tidal Love-number jump:

```bash
PYTHONPATH=compose_eos_pkg python scripts/run_phase_transition_tov.py \\
  --eos ./APR_CSS_Maxwell --pc 815.642916 --out ./star.json
```

The tidal jump uses Takátsy-Kovács (2020), including the finite transition-pressure term. The reported Seidov flag applies to the small-core onset criterion only, not a full radial stability calculation.

## Gibbs mixed phase

Gibbs construction requires two charge-resolved pressure surfaces `P(mu_B,mu_Q)` and a lepton sector:

```bash
PYTHONPATH=compose_eos_pkg python scripts/construct_gibbs_eos.py \\
  --hadronic-npz hadronic_surface.npz --high-density-npz quark_surface.npz \\
  --mu-b-min 900 --mu-b-max 1500 --out ./gibbs_eos
```

Each NPZ must contain `mu_b`, `mu_q`, and `pressure` arrays with shape `(len(mu_b),len(mu_q))`. The supplied CompOSE files are 1-D beta-equilibrium barotropes; a Gibbs EOS cannot be reconstructed from them, and the package does not pretend otherwise.

## Causality diagnosis and alternatives

Classify the 390 matched tables against the untouched source cores:

```bash
PYTHONPATH=compose_eos_pkg python scripts/classify_matching_failure_causes.py \\
  --source-root eos_data_normalized --matched-root ./matched_eos_v3_2 \\
  --output-json failure_causes.json --output-csv failure_causes.csv
```

An optional constant-sound-speed causal completion writes a new EOS and never overwrites the source:

```bash
PYTHONPATH=compose_eos_pkg python scripts/causal_complete_eos.py \\
  --input eos_data_normalized/APR_APR --out ./APR_APR_causal --cs2 0.95
```

The CLI runs the full consistency gate and PCHIP/cubic-spline/centered-FD causality checks. It tags a result `diagnostic_only` unless all pass. A causal completion is a model alternative, not a unique physical repair.

## Install and test

```bash
python -m pip install -e ./compose_eos_pkg
PYTHONPATH=compose_eos_pkg python -m pytest -q
PYTHONPATH=compose_eos_pkg python -m compileall -q compose_eos_pkg/compose_eos scripts tests
```

See `docs/TECHNICAL_REFERENCE_v3_3.md`, `docs/FAILURE_ROOT_CAUSE_REPORT_v3_3.md`, and `docs/ADVERSARIAL_AUDIT_REPORT_v3_3.md` for derivations, limits, evidence, and the bounded review/fix log.
