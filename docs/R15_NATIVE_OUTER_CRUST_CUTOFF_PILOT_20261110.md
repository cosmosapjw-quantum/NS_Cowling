# R15: native outer-crust finite-pressure cutoff pilot

**Status:** `CI_PENDING`  
**Branch:** `research/r15-native-outer-crust-pilot`  
**Base branch:** `research/mass-dependent-cowling-v1` at `f3635a589176075660649a33c5b28ba8fcc25906`

## Research question

How does the finite-pressure truncation response vary across crust families and core families when the low-density part of the supplied projected EOS tables is sampled at native knots, without rewriting EOS data?

This is a TOV-only diagnostic. It is not a P=0 completion, a physical crust-mass estimate, a complete-star radius calculation, or a radial/Cowling eigenproblem.

## Input authority and intake

The frozen v5 spectrum summary contains 234 unique projected EOS records:

- 13 crust families;
- 3 admitted core families: BBB_BHF2, GM_GM1, RG_SLY4;
- 6 historical smooth matching prescriptions;
- 39 complete crust-core pairs, each with all six methods.

The authoritative solver input is `oscillation_release_v1/data/eos_projected/`. The 234 projected tables have the existing `projected_native_cell_audit=PASS` recorded in the frozen spectrum summary and the prior native projected-table audit. The 234 normalized/original tables are retained as provenance inputs; they are not substituted for the projected solver inputs because the separate native-cell audit records that the full normalized tables are not admissible under the frozen native-knot cubic-Hermite (P(\mu_b)) representation.

The summary authority is pinned by SHA-256:

```text
ce54c111510f32086065f95356903cd6abad675adce12cab1a06b5cc82350965
```

No EOS table is rewritten, reordered, resampled, or extended by this experiment.

## Three-pair pilot design

The selected crossed contrasts are:

| Crust | Core | Role |
|---|---|---|
| GPPVA_FSU2 | BBB_BHF2 | R14 stress case; high p1 matching-response tail in the fixed-mass cohort |
| SPG_M2 | BBB_BHF2 | Change crust family at fixed core |
| GPPVA_FSU2 | GM_GM1 | Change core family at fixed crust |

Each pair retains all six historical smooth methods: `css`, `cubic_F`, `cubic_P`, `maxwell`, `quintic_F`, and `quintic_P`. Here `maxwell` is the legacy smooth bridge token, not a sharp first-order Maxwell density jump.

Total design:

[
3\;\text{pairs}\times6\;\text{methods}\times8\;\text{native knots}
=144\;\text{TOV integrations}.
]

This is a deliberately selected diagnostic pilot, not a representative random sample or a confirmatory estimate over the 39-pair cohort.

## Native low-density cutoff definition

For each projected table:

1. Read the existing rows in source order.
2. Require finite positive (P,\epsilon,n_b,\mu_b), and strict monotonicity of (P) and (n_b) in that source order.
3. Define the *operational low-density segment* as the contiguous source prefix with (n_b\le10^{-3}\,\mathrm{fm}^{-3}).
4. Form eight equally spaced targets in (log_{10} n_b), including the first and last density of that prefix.
5. For each target, select the nearest existing row. All EOS values used by the TOV solve, including (P_{\rm cut}), are copied from that row. Targets do not create new EOS knots.

The 18 projected tables inspected for this pilot have the same low-density source mesh through row 217. The eight selected source rows are:

[
(1,9,17,23,52,73,116,217).
]

The first source row has (n_b=4.67962\times10^{-10}\,\mathrm{fm}^{-3}) and (P=9.000000079422\times10^{-12}\,\mathrm{MeV\,fm}^{-3}). Row 217 is the last source row below the operational (10^{-3}\,\mathrm{fm}^{-3}) ceiling; row 218 is above it. In these tables the row-218 density is approximately (2.29\times10^{-3}) to (3.00\times10^{-3}\,\mathrm{fm}^{-3}), depending on crust family.

The (10^{-3}\,\mathrm{fm}^{-3}) threshold is a reproducible selection boundary for this cutoff diagnostic. It is not asserted to be the crust-core transition density, the pasta boundary, or a physical surface.

## Solver and acceptance

The runner uses the frozen v5 native Fortran `cell_solve task=star` backend:

- (r_0=10\,\mathrm{m});
- (r_{\rm limit}=10^8\,\mathrm{m});
- (mathrm{rtol}=10^{-9});
- (mathrm{step\_cap}=0.05);
- central pressure fixed to the corresponding frozen 1.3300106263-(M_\odot) reference row.

Acceptance requires:

1. exactly 18 projected source tables and 144 requested TOV cases;
2. all 18 projected-table SHA-256 values agree with the frozen summary before execution and remain unchanged afterwards;
3. each input retains its source row order and the exact eight selected native source rows above;
4. all 144 TOV integrations return finite positive enclosed mass and isobar radius;
5. all 18 row-1 native table-edge integrations reproduce frozen reference (M,R) within relative tolerance (2\times10^{-5});
6. all 144 backend-exported background profiles pass SHA-256 and row-count checks;
7. no synthetic EOS table is emitted.

The runner writes a JSON receipt, a case-level CSV, a paired six-method cutoff summary, a background-profile directory, and a hash manifest. A failure is preserved by case and stage; missing method rows are not imputed.

## Interpretation

For each cutoff the observables are (M(r_{\rm cut})) and (r_{\rm cut}), with (P(r_{\rm cut})=P_{\rm cut}>0). Changes are reported relative to the first native table-edge isobar of the same method. The six-method summary is a descriptive paired contrast inside each selected pair.

No result from this pilot can establish:
- the outer-crust mass or physical surface radius below the first table row;
- the actual crust thickness or crust-core transition location;
- radial stability or the sign of (\omega_0^2);
- complete-star Cowling f/p modes or a justified surface perturbation boundary condition.

The next full-cohort gate, if the pilot passes, is 234 tables × 8 cutoffs = 1,872 TOV integrations. It remains separate from the physical (P=0) completion and the radial-stability/Cowling gates.
