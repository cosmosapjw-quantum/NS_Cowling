# R15 full-cohort native low-density TOV cutoff scan

**Status:** `PASS_R15_FULL_COHORT_NATIVE_OUTER_CRUST_CUTOFF_SCAN`  
**TOV GitHub Actions run:** [42650666976](https://github.com/cosmosapjw-quantum/NS_Cowling/actions/runs/42650666976)  
**Evidence repack/verification:** [42651779185](https://github.com/cosmosapjw-quantum/NS_Cowling/actions/runs/42651779185)  
**Branch:** `research/r15-native-outer-crust-pilot`  
**TOV source commit:** `16f44113464987c2043d1b6c70a028f936af16e8`  
**Repack source commit:** `cf84ae994af6a99f422e47e51716ac4f584e37ae`  
**Frozen spectrum summary SHA-256:** `ce54c111510f32086065f95356903cd6abad675adce12cab1a06b5cc82350965`

## 1. Design and acceptance

The full cohort comprises 13 crust families × 3 admitted core families (BBB_BHF2, GM_GM1, RG_SLY4) = 39 crust-core pairs, with 6 historical matching prescriptions and 8 native low-density cutoff knots per projected EOS table. The legacy `maxwell` method token remains a smooth bridge prescription here, not a sharp first-order Maxwell density jump. The design contains:

[
39\times6\times8=1,872\quad\text{TOV integrations}.
]

The solver is the frozen v5 native Fortran `cell_solve task=star` backend with `r0=10 m`, `r_limit=1e8 m`, `rtol=1e-9`, and `step_cap=0.05`. Each table uses its frozen central-pressure reference associated with the 1.3300106263-(M_\odot) cohort. The cutoff sequence is selected at native source rows `(1, 9, 17, 23, 52, 73, 116, 217)`; no EOS row is inserted or modified.

Acceptance passed:
- 234/234 table hashes matched the frozen summary and remained unchanged;
- 1,872/1,872 cases returned status `ok`;
- all 234 native-edge baselines were reproduced;
- all 312 pair×cutoff six-method summary cells were complete;
- all 1,872 background profiles passed hash and row-count checks;
- failures = 0; synthetic EOS tables = 0.

## 2. Independent verification

The TOV artifact itself exceeded the connector's 512-MiB download limit, so it was downloaded and revalidated inside a dedicated GitHub Actions repack workflow. That workflow checked the preflight manifest, full TOV receipt, result CSV, pair summary, and all profile hashes before splitting the profiles into four independently downloadable shards. The split artifacts were then downloaded locally and independently checked.

- Native-edge baseline: 234/234 pass against relative tolerance `2e-5`; maximum observed relative errors are (1.3323\times10^{-15}) in mass and (1.2213\times10^{-15}) in radius.
- Source integrity: 234/234 projected EOS hashes are identical before/after and equal the frozen `projected_eos_sha256` values.
- Case coverage: 1,872 unique cases, 39 pairs, 6 methods, 8 cutoffs; 1,872 statuses are `ok`.
- Pair summaries: 312/312 complete six-method cells.
- Profiles: 1,872 unique files across four shards; 1,872/1,872 SHA-256 matches and 1,872/1,872 data-row counts match. All six downloaded ZIP artifacts pass integrity checks and their local SHA-256 values match GitHub artifact digests.

## 3. Response across the eight native cutoffs

Changes below are means over all 234 table-method cases at the indicated one-based source row. The denominator is the same method/table native-edge isobar at source row 1.

| Source row | Mean enclosed-mass change | Mean isobar-radius change |
|---:|---:|---:|
| 1 | +0% | +0.000000% |
| 9 | -1.06483e-09% | -0.044355% |
| 17 | -1.65007e-08% | -0.145093% |
| 23 | -3.02174e-07% | -0.388336% |
| 52 | -2.68331e-06% | -0.732000% |
| 73 | -3.15924e-05% | -1.384578% |
| 116 | -0.000335402% | -2.447973% |
| 217 | -0.00301608% | -3.934502% |

At source row 217, the mean enclosed-mass change is (-0.003016\%) and the mean isobar-radius change is (-3.934502\%). Across the 234 individual table-method cases, the row-217 mass response ranges from (-0.004714\%) to (-0.001968\%); the radius response ranges from (-4.629406\%) to (-3.472008\%).

## 4. Core-family and matching-method contrasts at source row 217

| Core family | Mean enclosed-mass change | Mean isobar-radius change | Radius range across 13 crust families |
|---|---:|---:|---:|
| BBB_BHF2 | -0.001972% | -3.475964% | [-3.478395, -3.472073]% |
| GM_GM1 | -0.004694% | -4.622658% | [-4.629344, -4.615399]% |
| RG_SLY4 | -0.002382% | -3.704885% | [-3.706966, -3.700646]% |

The core-family contrast is much larger than the spread among the six matching prescriptions for these low-density cutoff diagnostics. At source row 217, the maximum six-method spread over all 39 pair cells is (0.00147311) percentage points in isobar radius and (2.46122\times10^{-6}) percentage points in enclosed mass. The pressure and baryon-density cutoff values have zero six-method spread within every pair×cutoff cell.

The maximum row-217 cutoff has the common baryon density (n_b=2.69021\times10^{-4}\,\mathrm{fm}^{-3}); its pressure spans approximately (4.97534\times10^{-4}) to (4.97984\times10^{-4}\,\mathrm{MeV\,fm}^{-3}) across crust families.

For context, the mean row-217 radius response is about 1.146694 percentage points more negative for GM_GM1 than for BBB_BHF2. This is a descriptive contrast within the selected 1.3300106263-(M_\odot) reference cohort, not a universal causal statement about core EOS families.

## 5. Physical interpretation boundary

Every result is an internal finite-pressure isobar endpoint with (P_{\rm cut}>0): the reported quantities are (M(r_{\rm cut})) and (r_{\rm cut}). Source row 1 is already at positive pressure. Therefore the row-217 radius reduction is **not** a physical surface-radius correction, crust thickness, or full-star radius change. The enclosed-mass reduction is not a crust-mass estimate.

This run does not complete the EOS to (P=0), establish the crust-core transition or pasta boundary, determine radial stability or the sign of (\omega_0^2), or calculate complete-star Cowling f/p modes. Surface completion and finite-domain perturbation boundary conditions remain separate research gates.

The 217→218 density gap remains a separate low-density EOS representation issue. A 13-family representative endpoint diagnostic is recorded in `R15_ROW217_218_BOUNDARY_DIAGNOSTIC_20261110.md`; its two-point secant is not a local adiabatic index or sound speed.

## 6. Evidence artifacts

Original TOV artifact:
- ID `13689941895`, size 752,016,038 bytes;
- SHA-256 `df1ce1e3813128dc1b3419e95c14e05cefc4bad2cb4e68862aa2e48ac527cd64`.

Repacked evidence artifacts:
- Core artifact ID `13689692158`, SHA-256 `a867eb8311d2e34fbbf6a44ef567ad3ee397e862928ce661ab1faf83b6fe3519`.
- Profile shard 1 ID `13689732229`, SHA-256 `8b4bc51e94c61d2dc461f6d8908da246068c3a9480008f4d9f8259556bc6f872`.
- Profile shard 2 ID `13689966902`, SHA-256 `d04611a33e34bf9d3b73d000d39a2e9b52cb48dd5d4a53ae59db538fa2de493a`.
- Profile shard 3 ID `13690606395`, SHA-256 `66763d4f80c47672888154995877d9a85d907633858ac973a7a168c7a7ba1fba`.
- Profile shard 4 ID `13690366950`, SHA-256 `f7dac88ccd341af9f1e6af51b51128d4c7000cf66ade17672c7ce8dbb956cdc9`.
- Index artifact ID `13690696165`, SHA-256 `9cf31775e7f280a25cd910753db4c685148ce50fba492fed5bbf190b38886612`.

The shard manifests cover all 1,872 profile paths, hashes, and row counts. The compact core artifact contains the complete TOV JSON receipt, case CSV, 312-row pair-cutoff summary, TOV manifest, and full-cohort preflight evidence.

## 7. Next scientific gate

The exact solver-native row-217→218 representation audit is now complete across all 234 projected tables. Actions [42654562160](https://github.com/cosmosapjw-quantum/NS_Cowling/actions/runs/42654562160) passed the Fortran-equivalent whole-cell (P''>0) and (\mu P''-P'\\ge0) conditions in 234/234 cells. The analytic stationary equation for `chi=c_s^2/c^2` has zero interior roots in all cells; `chi` is monotonically increasing, with global range (1.830688\\times10^{-4})–(2.950720\\times10^{-3}). The maximum endpoint First-Law energy closure residual is (2.392\\times10^{-16}). No EOS source hashes changed, no synthetic rows were written, and no TOV integrations were run in this audit. Full report: `docs/R15_NATIVE_HERMITE_ROW217_218_AUDIT_20261110.md`.

R16 has now audited the normalized-to-projected transformation over all 234 tables; see `docs/R16_NORMALIZED_TO_PROJECTED_AUDIT_20261110.md`. Its result concerns the thermodynamic representation and convexity of the matched pressure knots, not upstream microphysical completeness. The next task is R17: compare explicitly constrained reconstruction objectives and trace raw secant/Hermite failures to source, phase and matching provenance before any new TOV or mode runs. The sourced (P=0) low-density EOS and radial boundary condition remain unresolved, so radial (\omega_0^2) and complete-star Cowling modes remain blocked..