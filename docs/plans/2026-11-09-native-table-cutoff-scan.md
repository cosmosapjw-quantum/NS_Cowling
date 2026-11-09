# Plan: native-table-only low-pressure boundary audit

## Objective

Quantify the sensitivity of fixed-central-pressure TOV integrations to the finite low-pressure endpoint of the supplied projected EOS, without modifying the EOS data.

## Inputs and authorities

- Frozen reference: `oscillation_release_v1/data/spectrum_summary_234.csv`.
- Six source tables: `oscillation_release_v1/data/eos_projected/GPPVA_FSU2__BBB_BHF2__{css,cubic_F,cubic_P,maxwell,quintic_F,quintic_P}.csv`.
- Solver authority: frozen v5 calculation/methods/support `cell_solve` path and native-knot cubic-Hermite P(mu_b) provider.
- Source table identity is the SHA-256 of the exact on-disk CSV bytes. Do not normalize and write back to a source path.

## Method

1. Read each CSV in source order. Require the physical columns `pressure_pa`, `energy_mev_fm3`, `n_b_fm3`, and `mu_b_mev` to be finite and positive.
2. Require strictly increasing pressure in source order. Do not sort or repair the table.
3. Select eight cutoff rows by nearest native pressure knot to offsets (0, 0.25, 0.5, 0.75, 1, 1.5, 2, 3) dex above P_min. Store the selected row index, source_row, exact P, epsilon, n_b, and mu_b.
4. For each method and selected knot, call `solve_v5_tov_star(source_csv, frozen_pc, p_cut_pa=exact_native_pressure)`. No intermediate EOS CSV is produced.
5. Compare the row-1 native table-edge result with the frozen M/R reference (2e-5 relative tolerance). For each higher-pressure cutoff, report enclosed mass and isobar radius relative to the native table-edge result.
6. Recompute each source-file hash after integration and before sealing the artifact. Fail the run if any source differs.
7. Export the frozen backend's native TOV background CSV for every cutoff, and verify the backend SHA-256 receipt and row count. These profiles are finite-domain structure products, not EOS replacements.

## Expected output

- 48 TOV cases (8 native pressure cutoffs x 6 matching methods).
- 48 native background profiles, one per finite-pressure integration endpoint.
- `native_table_cutoff_scan.json`, `native_table_cutoff_scan.csv`, `native_table_cutoff_sensitivity.csv`, background-profile CSVs, and `manifest.json`.
- No EOS tables in the artifact.

## Failure taxonomy

- `NATIVE_TABLE_INTAKE`: missing columns, invalid physical values, or non-monotone source ordering.
- `NATIVE_KNOT_SELECTION`: requested target offsets cannot be represented by distinct in-domain source knots.
- `NATIVE_TOV`: frozen backend fails at one pressure cutoff.
- `BASELINE_REPRODUCTION`: row-1 TOV M/R differs from frozen reference by more than 2e-5.
- `SOURCE_MUTATION`: source CSV hash changes during or after calculation.
- `BACKGROUND_PROFILE`: background export, schema receipt, row count, or SHA-256 validation fails.

## Stop condition

Stop the scan as scientifically failed if any native-edge baseline fails, any source hash changes, or the complete expected case matrix is not present. A failed scan can still preserve partial receipts, but no sensitivity conclusion is accepted.

## Background-profile handoff

The frozen v5 backend exports its native TOV background profile through `background_output_path`. The R14 runner requests one profile for every one of the 48 finite-pressure isobars and verifies the backend-reported SHA-256 and row count. These files are structural background evidence for a subsequent radial-stability/Cowling study, not new EOS tables. The R14 artifact must retain the explicit finite-pressure/truncated-domain interpretation; a profile ending at positive pressure cannot be relabelled as a complete stellar background.


## Exact-head result (2026-11-09)

Implementation and exact-head workflow validation are complete. GitHub Actions run [42640660436](https://github.com/cosmosapjw-quantum/NS_Cowling/actions/runs/42640660436) passed on branch HEAD `aff4e6ad2789555648f2e564a96bbabd952ee4ac`.

- All 48 requested native TOV cases completed; no case failures.
- All six original 323-row EOS source tables retained their exact before/after SHA-256 identities and expected Git blob SHA1s.
- All six source-row-1 TOV baselines reproduced the frozen M/R reference within the unchanged 2e-5 criterion.
- All 48 background profile CSVs were exported by the frozen v5 TOV backend; all backend SHA-256 receipts and row counts match the downloaded artifact.
- No synthetic EOS table was emitted.
- Downloaded artifact: ID `13684872040`, 18,300,545 bytes, SHA-256 `e124339e932636e9e08c92298b2e1269aca3fac4980230083a28a3451421d8f3`.

Across the six methods, moving the integration endpoint from source row 1 to source row 20 changes the enclosed mass by only about (-4.54\times10^{-10}) relative to the native edge, while the isobar radius decreases by about 24.1 m (approximately 0.2137%). This is a within-table truncation result; it is not the physical missing-crust radius or the P=0 surface correction.

## Supersession and next decision

R13's synthetic low-density tails are retained only as historical mathematical diagnostics. They are not used as EOS input for the primary study and are not selected as a physical crust completion. The R13 workflow has been switched to manual dispatch only.

The next decision is whether to (a) obtain and validate a source-backed low-density crust EOS that is compatible with the supplied table, or (b) run a separately labelled finite-domain radial/Cowling problem whose boundary is explicitly placed at a positive-pressure native isobar. Neither route may be described as a complete P=0 stellar model until a justified physical low-density EOS and boundary condition are supplied. The current TOV cutoff scan alone does not admit a radial-stability or full-star Cowling claim.
