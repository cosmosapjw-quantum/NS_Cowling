# Change and repair report - compose_eos v3.3.0

## New implementation

- Added `phase_transitions.py`: exact Maxwell coexistence from two 1-D locally neutral barotropes; a single-hadronic-table plus phenomenological CSS high-density phase builder; a 2-charge Gibbs mixed-phase constructor; pressure-corrected tidal jump; phase-aware TOV/Love integration; JSON sidecars and readback validation.
- Added `causality.py`: three-scheme source `c_s^2` profile and explicit constant-sound-speed causal completion, with the original input untouched.
- Added scripts: `classify_matching_failure_causes.py`, `construct_maxwell_eos.py`, `construct_css_maxwell_eos.py`, `construct_gibbs_eos.py`, `run_phase_transition_tov.py`, and upgraded `causal_complete_eos.py`.
- Added phase-transition tests for Maxwell equilibrium, density/energy jump, Seidov threshold, Gibbs pressure/global charge equality and Hessian convexity, write/read sidecars, the corrected tidal jump, phase-aware TOV, and causal completion.

## Important physical contracts

- The old `method="maxwell"` matching token remains a smooth bridge for compatibility. The new `build_maxwell_eos` and `build_maxwell_css_eos` are separate sharp-interface constructions.
- Gibbs requires two-dimensional `P(mu_B,mu_Q)` pressure surfaces. The provided CompOSE tables are one-dimensional beta-equilibrium slices, so only analytic synthetic Gibbs surfaces were used for the regression test; no Gibbs EOS was claimed for the provided data.
- CSS high-density matter is phenomenological and has no microscopic particle fractions. The sidecar states that limitation.
- A causal completion changes the source EOS. APR's tested alternative is cross-scheme causal; KDE0v1's tested `c_s^2=0.80` output remains scheme-sensitive and diagnostic-only.

## Source-failure diagnosis

The prior 390 matched EOS outputs were re-read and re-gated. All 156 FAILs are source-core acausality preserved in the matched tail: 78 APR_APR and 78 RG_KDE0v1. There are zero matching/numerical-candidate failures in that matrix. See `FAILURE_ROOT_CAUSE_REPORT_v3_3.md` and its machine-readable JSON/CSV.

## Tests

- Phase-transition focused tests: 7 passed.
- Full regression suite: 54 passed, one non-fatal Matplotlib 3D-projection import warning.
- Python `compileall`, wheel build and isolated import/CLI smoke test: recorded in `DELIVERY_MANIFEST_v3_3.json`.
- Actual APR+CSS Maxwell table: CompOSE write/read and all hard gates passed.
- Actual phase-aware TOV/Love integration: three settings; jump formula independently reconstructed from stored interface fields.

## Non-claims

No claim is made that all 390 original EOS are causal; no claim is made that Gibbs can be reconstructed from the provided one-dimensional data; no claim is made that the CSS phase is microscopic QCD; no 390-model TOV campaign or full radial-mode stability computation was performed.
