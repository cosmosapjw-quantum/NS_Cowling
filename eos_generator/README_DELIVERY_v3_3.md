# compose_eos v3.3.1 delivery

This bundle contains the `compose_eos` v3.3.1 source package, raw and normalized cold CompOSE tables, scripts, tests, validation/technical reports, machine-readable audit evidence, validation examples, the v3.3.1 wheel, and provenance/readback records.

## Main results

- `compose_eos` version aligned to `3.3.1` across package metadata, CLI (`--version`), docs, wheel, and SHA-256 manifests.
- Phase 2/3/4/5/6/7/8/9 remediation completed:
  - CompOSE cold reader rejects single $T>0$ slices (`NotImplementedError`) and infers `beta_equilibrated` strictly from $T=0 \wedge I_l=1 \wedge \max|Q_5|<10^{-8} \wedge \text{charge-neutral}$.
  - Gibbs construction (`compose-eos-gibbs-v2`) enforces exact knot-union root isolation, 12-field `ChargeResolvedPhase` provenance, lepton double-counting checks, `(939.0, 939.0) MeV` normalization reference labeling, `composition_status = "thermo_only"`, and post-construction / post-readback `validate_gibbs_barotrope` over `PchipInterpolator`.
  - Causality (`causality_profile`, `causal_completion`) uses exact piecewise-cubic interval extrema via the degree-$\le 6$ stationary polynomial, distinguishing `true_crossing`, `tangential_contact`, and `none`.
  - All CompOSE writers perform staged serialization + full `verify_compose_readback` before atomic publication.
  - Native Fortran CLI `compose_eos_fortran reconstruct-compo-native` provides raw-SPG/GPPVA algebraic C3/C4 closure with native readback and `0.000e+00` max difference against `normalize_compose_composition.py` (explicitly bounded as algebraic closure, not spatial Thomas-Fermi).

## Contents

- `compose_eos_pkg/`: importable package v3.3.1
- `scripts/`: normalizer, failure classifier, Maxwell/CSS/Gibbs construction, causal completion, phase-aware TOV, and validation report generators
- `tests/`: regression, adversarial, and phase-transition test suite
- `eos_data/`, `eos_data_normalized/`: raw and normalized EOS source data, kept separate
- `docs/`: technical reference, scope/limitations, v3.3.1 Gibbs/causality/provenance reports, historical v3.3.0 audit reports
- `receipts/`: per-operation staged readback receipts (`readback_receipts_v3_3_1.json`)
- `wheels/`: `compose_eos-3.3.1-py3-none-any.whl`
- `provenance/`: v3.2 provenance metadata (note: large external archive `EoS_revised_code_package_upgraded_v3_2.zip` is excluded from git tracking and omitted from `SHA256SUMS_v3_3.txt`)

## Installation

```bash
python -m pip install wheels/compose_eos-3.3.1-py3-none-any.whl
PYTHONPATH=compose_eos_pkg python -m pytest -q tests
```

See `docs/TECHNICAL_REFERENCE_v3_3.md` and `docs/EOS_GENERATOR_SCOPE_AND_LIMITATIONS.md` for equations, CompOSE v3.01 semantics, contracts, and reproduction commands.
