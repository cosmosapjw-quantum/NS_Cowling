# R15 full-cohort native-density preflight

**Status:** `PASS_R15_FULL_COHORT_NATIVE_DENSITY_PREFLIGHT`  
**GitHub Actions run:** [42650145210](https://github.com/cosmosapjw-quantum/NS_Cowling/actions/runs/42650145210)  
**Branch:** `research/r15-native-outer-crust-pilot`  
**Preflight source commit:** `ddfa64585bab4aeed4a1102c431027b9efcbfb30`  
**Frozen spectrum-summary SHA-256:** `ce54c111510f32086065f95356903cd6abad675adce12cab1a06b5cc82350965`

## Scope

This is a read-only intake and native-knot selection preflight over the full projected solver-input cohort. It does not execute TOV integrations, modify EOS tables, or generate synthetic EOS rows.

The frozen cohort contains 234 projected tables: 13 crust families, 3 admitted core families (BBB_BHF2, GM_GM1, RG_SLY4), 39 complete crust-core pairs, and 6 historical smooth matching prescriptions per pair.

## Acceptance results

- 234/234 projected CSV tables found; inventory agrees exactly with the frozen 234-row spectrum summary.
- All 234 projected EOS SHA-256 values match the frozen summary and remain unchanged after inspection.
- All source `source_row` values are contiguous, one-based, and agree with source-file order.
- All 234 native tables pass the monotonic positive-pressure and baryon-density requirements used by the cutoff selector.
- 1,872/1,872 native cutoff knots selected, exactly eight per table.
- Zero preflight failures; zero TOV integrations; zero source-table modifications; zero synthetic EOS tables.

The workflow completed successfully. The artifact contains the full JSON receipt, a 234-row table audit CSV, a 1,872-row selected-knot CSV, and a SHA-256 manifest.

## Native mesh findings

Every projected table shares the same eight selected source rows:

- zero-based row indices: `[0, 8, 16, 22, 51, 72, 115, 216]`;
- one-based source rows: `(1, 9, 17, 23, 52, 73, 116, 217)`.

The source tables have 323–687 rows, but the contiguous prefix satisfying (n_b \le 10^{-3}\,\mathrm{fm}^{-3}) contains exactly 217 rows in all 234 tables. The maximum adjacent spacing inside this prefix is constant across the cohort:

[
\max_i \Delta\log_{10} n_{b,i}=0.15499445366657\;\mathrm{dex}.
]

The first ineligible row follows immediately after source row 217. The density jump ratio is invariant under core family and matching method within each crust family:

| Crust family | (n_b(218)/n_b(217)) |
|---|---:|
| GPPVA_DD2 | 7.434364 |
| GPPVA_DDME2 | 7.434364 |
| GPPVA_FSU2 | 11.151546 |
| GPPVA_FSU2H | 7.434364 |
| GPPVA_FSU2R | 7.434364 |
| GPPVA_NL3wrL55 | 7.434364 |
| GPPVA_TM1e | 7.434364 |
| GPPVA_TW | 11.151546 |
| SPG_M1 | 9.337139 |
| SPG_M2 | 8.515572 |
| SPG_M3 | 9.777185 |
| SPG_M4 | 9.777185 |
| SPG_M5 | 9.337139 |

The observed jump therefore tracks the crust-family low-density mesh, rather than the selected core family or smooth matching token. This is a useful diagnostic distinction: the cutoff sequence is common and unambiguous, while the first interval beyond the operational low-density prefix is substantially more sparsely sampled and family-dependent.

## Interpretation boundary and next gate

The result establishes a common native-knot selection rule over all 234 projected tables. It does **not** identify the physical crust-core transition, complete an EOS to (P=0), provide a physical surface radius or crust mass/thickness, or validate radial stability or Cowling-mode boundary conditions.

The full-cohort TOV stage is now complete: 39 pairs × 6 methods × 8 cutoffs = 1,872 integrations passed in [run 42650666976](https://github.com/cosmosapjw-quantum/NS_Cowling/actions/runs/42650666976). Independent receipt, source-hash, native-edge baseline, pair-cell, and all 1,872 background-profile hash/row-count checks also passed after evidence repackaging. See [R15_FULL_COHORT_TOV_RESULTS_20261110.md](R15_FULL_COHORT_TOV_RESULTS_20261110.md) for the full results and artifact provenance.

The large source-row-217-to-218 density interval should remain a separate low-density EOS-segment/interpolation audit; it must not be interpreted as a physical layer thickness.

## Artifact provenance

- GitHub artifact ID: `13690310206`
- Artifact SHA-256: `f2d1cbc22d4d1307b1c71ab3c6005181df3d76e9b5baa0d7da79d98a57a5888b`
- Local ZIP SHA-256 independently recomputed: `f2d1cbc22d4d1307b1c71ab3c6005181df3d76e9b5baa0d7da79d98a57a5888b`
- Receipt SHA-256: `2efd9820d95a0d173f83daec92a3d6ca70d709f63fb58cbb21ab57e4dca758a2`
- Table-audit CSV SHA-256: `5b3a367b92c4b42bccb596a8249a74d88a4de2bbf8e8e643cc6c2b29197b298a`
- Selected-knot CSV SHA-256: `708f94ca8a038774a21b1fdff53cb9582beaf330719676ba2abc905ad134c5c9`
