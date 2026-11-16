# Observed phase report audit — bounded preview

Verdict: **PASS_SCOPED** for the data-only actual-profile discovery and transition diagnostic arithmetic in the 106-node preview. This is not final campaign admission, an optimizer energy-error bound, a continuum transition, a common-tangent construction, original author calibration, or a global-minimum certificate.

`OBSERVED_PHASE_REVIEW.json` records the exact report, 16 raw-input and 13 scientific source hashes. The original CLI and prior transition-review receipts were not changed. No production module or solver output was edited and no cell solver was run.

The independent runner replayed all 106 canonical selections and checked each record pointer and descriptive morphology against the retained profiles. There were 104 selected nodes and two preserved FSU2H N192 HOLD nodes at 0.095 and 0.100 fm^-3. It independently constructed every adjacent selected-label change, then checked all 30 windows, 60 samples, 71 available mesh branch pairs and seven raw fine-mesh sign brackets. Pair energies were selected afresh from stationary feasible raw attempts by final descriptive morphology. All energy differences, adjacent-mesh changes, 64-epsilon roundoff scales and bracket flags match exactly.

The new DDME2 rod–slab bracket [0.070, 0.0725] fm^-3 is present because discovery follows actual selected profiles instead of historical requested pairs. At N384, slab-minus-rod E/A is +1.9571150380670588e-4 and -1.9890657347332308e-4 MeV/A; the respective measured N192-to-N384 changes are 7.73843567003496e-9 and 3.1551508072880097e-9 MeV/A. The roundoff scales are about 1.346e-11 MeV/A. These finite diagnostics support the recorded sign bracket; KKT and field residuals have not been translated into an energy uncertainty and do not establish continuum branch ordering.

Near-uniform labels and branch disappearance remain descriptive changes without fabricated sign roots. FSU2/TW uniform-versus-near-uniform label differences are retained separately. Synthetic isolated tests establish that a HOLD node blocks adjacent-change construction, repeated identical anchors deduplicate, and conflicting repeated anchors reject. Raw sign brackets require both named branches at both adjacent endpoints; missing branches are not filled by interpolation.

One source-manifest omission was found: `branch_minima` uses `good()` imported from `scripts/refine_gppva_tf.py`, but the first data-only report did not hash that dependency. The geometry owner expanded the manifest to the full canonical dependencies and regenerated the preview without changing scientific inputs or logic. The replay verifies that repaired 13-source manifest.

Reproduce from repository root:

```sh
OPENBLAS_NUM_THREADS=1 PYTHONPATH=python_frontend:scripts python data/research_finish_20261116/independent_review/cli_review/observed_phase_review/review_observed_phases.py
```

The runner accepts a report path as its first argument, but writes to its own audit directory. Preserve this preview receipt before applying it to a later final report. Final 156-node union and final native export were pending at this review and are outside this verdict.
