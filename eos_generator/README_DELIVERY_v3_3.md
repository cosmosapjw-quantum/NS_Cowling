# compose_eos v3.3.0 delivery

This bundle contains the complete v3.3 source package, raw and normalized cold CompOSE tables, scripts, tests, four reports in Markdown and PDF, machine-readable audit evidence, validation examples, the built wheel, and the preserved v3.2 matched-390 data archive and source-package provenance.

## Main results

- Full 390-way matching matrix independently re-read and re-gated: **234 PASS, 156 SOURCE_CORE_ACAUSALITY_PRESERVED**. The 156 original failures remain FAIL and trace to high-density APR_APR/KDE0v1 source tails, not to the crust-core matching bridge.
- Full suite: **54 passed**. Python compileall: PASS. Wheel build, clean-target install, API imports and both console entrypoints: PASS.
- Sharp Maxwell construction, phenomenological CSS high-density phase, two-charge Gibbs API, pressure-corrected tidal jump, and phase-aware TOV/Love solver are implemented. Gibbs requires 2-D charge-resolved phase surfaces; it is not fabricated from the supplied 1-D beta-equilibrium tables.
- An APR+CSS Maxwell table is written/read and passes the EOS gates. Three TOV integrations and an independent tidal-jump recomputation are recorded. APR 0.95 causal completion is robust under PCHIP/cubic/centered-FD diagnostics. KDE0v1 0.80 completion remains diagnostic-only.

## Contents

- `compose_eos_pkg/`: importable package v3.3.0
- `scripts/`: normalizer, failure classifier, Maxwell/CSS/Gibbs construction, causal completion and phase-aware TOV commands
- `tests/`: full regression and phase-transition tests
- `eos_data/`, `eos_data_normalized/`: raw and normalized EOS source data, kept separate
- `docs/`: four full reports/references in `.md` and `.pdf`, root-cause CSV/JSON, audit trace, tests and wheel logs
- `validation_examples/`: actual APR+CSS Maxwell EOS, TOV outputs, and causal-completion examples
- `wheels/`: `compose_eos-3.3.0` wheel
- `provenance/`: prior v3.2 source package, matched-390 archive and reproducibility/readback records

## Installation

```bash
python -m pip install wheels/compose_eos-3.3.0-py3-none-any.whl
python -m pytest -q
```

See `docs/TECHNICAL_REFERENCE_v3_3.md` for equations, units, contracts, and reproduction commands. See `docs/FAILURE_ROOT_CAUSE_REPORT_v3_3.md` for the 390-case classification.
