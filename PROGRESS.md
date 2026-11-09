# Neutron-star nonradial oscillation follow-up: final checkpoint

> Scope note (2026-11-08): this is the historical v2 234-EOS f/p1 Cowling-spectrum checkpoint. COMPLETE applies only to that cohort, not active v4 EOS work. No new stellar-mode integrations are claimed by the v4 EOS changes.


**Status: COMPLETE for the declared calculation scope. No f/p1 integrations were rerun during resumption.**

## Source and cohort
- Input: /home/cosmosapjw/Dropbox/Neutron_Star_Oscillations_Seminar.zip
- Input SHA-256: b2c0dcaf6a5b07956cf8d33c39651d4ade32f44b4c49b33aa08fd84fb458185b
- Matching matrix: 390 candidates = 13 crusts x 5 core candidates x 6 matching prescriptions.
- Admitted: 234 source/join cases (78 each for BBB_BHF2, RG_SLY4, GM_GM1).
- Excluded: 156 source-core acausality cases (78 APR_APR and 78 RG_KDE0v1).

## Calculation
- ell=2 polar-fluid modes in the relativistic Cowling approximation.
- Cold beta-equilibrated barotropic perturbation closure; fixed spherical, nonrotating perfect-fluid background.
- Finite positive EOS pressure cutoff; no extrapolation to P=0.
- Target mass: 2.644606257493485e30 kg = 1.3300106263127078 nominal solar masses.
- Two-sided shooting/matching. f identity: W,V nodes=(0,0); p1 identity: W,V nodes=(1,1).
- All 234 unique EOS results are present and pass the declared node-identity and native projected-EOS cell audits.
- Parallel recovery already completed: 99 targeted reruns (27 previously incomplete + 72 GM_GM1 p1 identity corrections) with 12 workers.
- Representative convergence already completed: 3 EOS x 3 tolerance levels with 9 workers.

## Aggregate numerical diagnostics
- Maximum relative target-mass residual: 9.914291609902648e-13.
- Maximum absolute normalized matching residual: f 4.560284933008418e-08; p1 3.0891943864075344e-10.
- Maximum inward/outward matching-radius mismatch: f 7.275957614183426e-12 m; p1 9.094947017729282e-12 m.
- Maximum relative matching-mass mismatch: 6.386016275503321e-16.
- Density knots preserved by thermodynamic projection to roundoff; maximum reconstructed pressure change 1.2930514536147284e-3 P_max.
- All 234 projected EOS tables pass the native interpolation-cell audit.
- Representative loose-to-tight maximum frequency changes: f 9.35e-8 Hz; p1 1.71e-6 Hz. This is not a global error bound.

## Median frequencies and full ranges
- BBB_BHF2: f 2.5051 kHz [2.4883, 2.5136]; p1 7.0083 kHz [6.6076, 7.1823]; R 11.2684 km [11.2653, 11.2738].
- RG_SLY4: f 2.3698 kHz [2.3574, 2.3803]; p1 6.7177 kHz [6.3471, 6.8906]; R 11.7706 km [11.7670, 11.7745].
- GM_GM1: f 1.9202 kHz [1.8893, 1.9401]; p1 6.0090 kHz [5.6866, 6.0773]; R 13.8079 km [13.7980, 13.8215].

## Scope and caveats
- These are Cowling fluid frequencies, not full-GR quasinormal modes; no metric perturbations, gravitational-wave damping, rotation, magnetic field, elastic crust, superfluidity, frozen-composition response, or sharp phase-interface transfer modes are included.
- The maxwell suffix in the source matching matrix is a legacy smooth-matching label, not a sharp interface in the oscillation equations.
- The pressure/energy reconstruction is an explicit derived interpolation representation; original normalized EOS tables are preserved separately.
- The surface is a finite positive pressure cutoff, not an exact vacuum surface.
- The REVTeX author field remains a placeholder; replace it before submission. Two minor overfull hbox warnings are recorded in the validation report.

## Deliverables
- Release folder: /tmp/ns_oscillation_followup_20261106/DELIVERY_NS_COWLING_F_P1_v1
- Data archive: DELIVERY_NS_COWLING_F_P1_v1/NS_Cowling_f_p1_data_tables_v1.zip
- REVTeX PDF: DELIVERY_NS_COWLING_F_P1_v1/NS_Cowling_f_p1_REVTeX.pdf
- Full archive: /tmp/ns_oscillation_followup_20261106/NS_Cowling_f_p1_full_release_v1.zip
- Full archive SHA-256: 47ed76103ec510262ffa082dc8e78d9990f6511d623f1a43200bf4738ca4c83d
- Dropbox-synced folder: /home/cosmosapjw/Dropbox/NS_Cowling_f_p1_20261106
- Downloads copies: /home/cosmosapjw/Downloads/NS_Cowling_f_p1_data_tables_v1.zip, NS_Cowling_f_p1_REVTeX.pdf, NS_Cowling_f_p1_full_release_v1.zip

## Resume verification
- Original input and provenance copy have matching SHA-256.
- Archive contains 234 original normalized EOS CSVs, 234 projected EOS CSVs, and 234 projection receipts.
- Summary has 234 unique rows; 78 per core family; mode node identities and projected-EOS audit all PASS.
- Data and full-release ZIP tests pass.
- REVTeX PDF has 5 pages; pdfinfo and pdftotext checks pass.
- All 1,014 delivery-manifest entries verify.
- Dropbox local sync client reports Up to date; synced-folder and Downloads key-file checksums pass.
- BibTeX validation: 5 entries, 0 issues. Abstract logic review: 0 issues.

**Next action:** no f/p1 integrations remain for this v1 scope. Full-GR perturbations, frozen-composition closure, sharp-interface modes, or a zero-pressure surface are separate physical models and should use separately versioned outputs.


## R15 full-cohort native low-density TOV cutoff scan (2026-11-10)

**Status: PASS for the declared finite-positive-pressure truncation diagnostic.** This is not a P=0 surface completion, physical crust-mass/thickness estimate, radial-stability result, or complete-star Cowling spectrum.

- Full-cohort native-density preflight: Actions [42650145210](https://github.com/cosmosapjw-quantum/NS_Cowling/actions/runs/42650145210), 234/234 projected tables, 1,872/1,872 selected native knots, all source hashes unchanged.
- Full TOV scan: Actions [42650666976](https://github.com/cosmosapjw-quantum/NS_Cowling/actions/runs/42650666976), source commit `16f44113464987c2043d1b6c70a028f936af16e8`; 39 pairs × 6 methods × 8 cutoffs = 1,872/1,872 TOV cases, all status `ok`.
- Independent source/baseline/case/profile verification: 234/234 source tables unchanged; 234/234 native-edge M/R baselines reproduced (max relative errors 1.33e-15 mass, 1.22e-15 radius); 312/312 pair×cutoff cells complete; 1,872/1,872 background profile SHA-256 and data-row counts match; zero failures; zero synthetic EOS.
- Original full TOV artifact is 752 MB and exceeds the connector download cap. Repack workflow [42651779185](https://github.com/cosmosapjw-quantum/NS_Cowling/actions/runs/42651779185) verified the source artifact and emitted a compact core artifact, four profile shards, and an index. All six ZIPs were downloaded and their digests/integrity checked.
- At the last selected cutoff (source row 217), mean enclosed-mass response is -0.00301608% and mean internal-isobar radius response is -3.934502%, each relative to the same method/table source-row-1 finite-pressure isobar. The response is descriptive within the frozen 1.3300106263-Msun reference cohort; it is not a physical surface/crust correction.
- At row 217 the core-conditioned mean radius responses are BBB_BHF2 -3.475964%, GM_GM1 -4.622658%, RG_SLY4 -3.704885%. Six-method spread is much smaller than these core-family contrasts in this selected diagnostic.
- Row 217→218 remains a sparse, crust-family-dependent native interval. Its representative two-point secant is not a local sound speed or causality test.

Detailed result and artifact provenance: `docs/R15_FULL_COHORT_TOV_RESULTS_20261110.md`. Native mesh intake: `docs/R15_FULL_COHORT_NATIVE_DENSITY_PREFLIGHT_20261110.md`. Boundary diagnostic: `docs/R15_ROW217_218_BOUNDARY_DIAGNOSTIC_20261110.md`.

**Next numerical question:** inspect the exact native cubic-Hermite (P(\mu_b)) representation over row 217→218, including derivative positivity/monotonicity and whole-cell stationary points. Do not append synthetic rows or infer physical crust/surface completion. MD-31 (sourced P=0 EOS/surface) and MD-32 (finite-domain perturbation BC) remain separate blockers; radial (\omega_0^2) and complete-star Cowling modes remain blocked.

## R15 solver-native cubic-Hermite row-217→218 audit (2026-11-10)

**Status: PASS for the currently admitted interpolant over the audited cell; physical outer-crust completeness remains unresolved.**

- Exact implementation inspected: `ns_mu_bridge.f90`, `ns_cell_eos.f90`, and `ns_table_eos.f90`; interpolation identity `cubic_Hermite_P_of_mu_raw_n_v1`.
- Workflow [42654562160](https://github.com/cosmosapjw-quantum/NS_Cowling/actions/runs/42654562160), audit runner HEAD `f4752ce5c6f0b79741dc79a672af0e39516164cf): tests and 234-table row-217→218 audit passed.
- All 234 cells pass the exact Fortran-equivalent endpoint conditions: (P''>0) and (G=\mu P''-P'\ge0) over the complete cubic interval. Independent analytic stationary-root enumeration finds zero interior stationary points of (chi=c_s^2/c^2); (chi) is monotonically increasing in all 234 cells.
- Whole-cell (chi) range: (1.830688\times10^{-4}) to (2.950720\times10^{-3}); equivalent (c_s/c=0.0135303)–0.0543205. Maximum endpoint First-Law energy closure residual is (2.392\times10^{-16}).
- All projected source hashes remain unchanged; no synthetic rows and no TOV integrations were used in this representation audit.
- Evidence report: `docs/R15_NATIVE_HERMITE_ROW217_218_AUDIT_20261110.md`; artifact ID `13689769387`, ZIP SHA-256 `24e0bd37d8acf81f28443c5e5cd6a739b91cd3f5a0c9ca78f6e5b1fbbb6e549f`.

**Interpretation:** this excludes a hidden non-convex, non-monotone-density, or acausal segment in the current solver-native cubic-Hermite representation across rows 217–218. It does not establish that the supplied EOS is a complete physical outer crust or that its low-density source provenance is complete. The next task is R16 source/provenance tracing of CompOSE/BSk22-derived low-density records and their projection into the first 218 native rows. MD-31/32, radial (\omega_0^2), and complete-star Cowling modes remain blocked.

## R16 full-cohort normalized-to-projected EOS and convexity audit (2026-11-10)

**Status: PASS for the declared reconstruction and solver-cell acceptance; not a statement that the projection is uniquely physical.**

- Actions [42666868668](https://github.com/cosmosapjw-quantum/NS_Cowling/actions/runs/42666868668), code HEAD `858c5d5639e0140e53c6f27b6e504eb261c05b01`; 234 normalized source tables and 234 projected tables audited.
- Independent implementation reproduces all 234 projected tables with maximum relative mismatch 0.0; source/projected SHA-256 values remain unchanged. The full audit has 109,776 adjacent intervals.
- Source (q=n_B/\mu_B) is nondecreasing in all 234 tables: PAV pooling count = 0. The projector's substantive change here is pressure reintegration (P(\mu_B)=P_0+\int n_B\,d\mu_B) and First-Law energy reconstruction, not pooling density knots.
- Source-normalized Hermite cells: 24,465/109,776 fail strict convexity; 24,505/109,776 fail stability or causal-margin acceptance, including 40 stable-but-acausal cells. Projected Hermite cells: 0/109,776 rejected.
- In the source-row-217→218 cell, 144/234 source tables fail (all eight GPPVA families, all 18 tables/family); all five SPG families pass. The projected cell passes in 234/234 tables.
- Source pressure secant slopes (s_i=\Delta P_i/\Delta\mu_i) contain 447 downward transitions (relative decrease threshold (10^{-12})) across all 234 tables; projected pressure secants have zero decreases. The worst source decrease is -21.3709% at the GPPVA_TM1e row-290 knot. This is a discrete convexity obstruction: no convex interpolant can pass through all source pressure knots exactly when successive adjacent secants decrease.
- Density localization (geometric endpoint midpoint): raw cell rejection is 14.47% below (10^{-4}\,\mathrm{fm}^{-3}), 26.51% at (10^{-4})–(10^{-3}), 82.87% at (2\times10^{-3})–(10^{-2}), 44.27% at (10^{-2})–(10^{-1}), 2.99% at (10^{-1})–0.5, and 0% above 0.5. There are no cells in the (10^{-3})–(2\times10^{-3}) bin because of the native mesh gap.
- Row-218 pressure is shifted by +11.54–64.08%; the row-217→218 pressure increment by +13.87–85.93%, relative to normalized source knots. No TOV integrations, synthetic rows, or mode calculations were performed in R16.
- Artifact ID `13697021829`, ZIP SHA-256 `ad9d9858e15ccb013a9591ef63d8e36528202a6ef9d237284765f05c5bbcc8e1`; locally re-downloaded archive and all six manifest payload hashes verified.
- Detailed evidence: `docs/R16_NORMALIZED_TO_PROJECTED_AUDIT_20261110.md`.

**Interpretation:** the original normalized knots satisfy the algebraic identity (\epsilon+P=\mu_B n_B) to machine precision, but this alone does not establish differential compatibility (dP/d\mu_B=n_B). The secant obstruction concerns the stored matched/normalized pressure-knot sequence, not proof that the upstream microphysical EOS is physically unstable. The raw CompOSE package files, normalized matched tables and projected solver tables remain separate provenance identities.

**Next:** R17 should compare constrained thermodynamic reconstructions with declared data-fidelity objectives and explicit movable variables. Preserve the current q-preserving exact-integral projector as baseline A; do not invent statistical errors or weights. Trace the source/phase/composition provenance of raw knot failures before launching further TOV or mode sensitivity runs. MD-31/32, radial (\omega_0^2), and complete-star Cowling modes remain blocked.

## R17 constrained joint-q reconstruction Pareto pilot (2026-11-10)

**Status: PASS for the bounded four-table implementation/admissibility pilot; negative/limited scientific outcome for this objective family.**

- Actions [42669846775](https://github.com/cosmosapjw-quantum/NS_Cowling/actions/runs/42669846775), code HEAD `9fde7e2206a162530c03a77098c5fa7fa19eeba7); 4 representative BBB_BHF2 tables × 5 lambda values × 2 quadrature coordinates = 40 candidates, 14,740 knot rows.
- Quadrature measures: `mu` and `log_n`; lambda values `1, 0.1, 0.01, 0.001, 0.0001`. Lambda and coordinate are sensitivity parameters, not statistical weights.
- All 40 candidates have nondecreasing positive q, positive/increasing P and epsilon, zero failed Hermite cells, zero causal-margin rejects, and zero downward candidate pressure-secant changes.
- Lambda=1 reproduces the frozen projected parent EOS to floating-point roundoff; all four source hashes and all four projected-parent hashes remain unchanged.
- Maximum relative KKT residual is 9.9943e-13 against the 1e-12 stopping tolerance; maximum First-Law residual is 2.2011e-16.
- Mean source-relative pressure RMS mismatch changes only from 2.967337% at lambda=1 to 2.966241% for mu weighting and 2.967068% for log_n weighting at lambda=1e-4. Mean density RMS movement is 0.0801% (mu) and 0.0308% (log_n). The row-218 pressure discrepancies remain +64.0802% (GPPVA_DD2), +30.3615% (GPPVA_TM1e), and +15.8415% (SPG_M1) relative to normalized source knots.
- Interpretation: this constrained joint-q fit is feasible but does not materially resolve the source pressure-knot mismatch. Log-density weighting reduces density movement but gives even less pressure-RMS improvement. This is a scoped negative result for the four BBB_BHF2 examples and this objective family, not a universal no-go.
- No TOV integrations, mode calculations, synthetic EOS rows, or primary EOS writes were performed.
- Final reporting workflow [42669846775](https://github.com/cosmosapjw-quantum/NS_Cowling/actions/runs/42669846775), artifact ID `13697898562`, SHA-256 `8a266c0b5c1316391fd6ede3d5ae65904b869a8dab66fcb29230a3ede5d52799`. The downloaded ZIP, archive integrity and all manifest hashes independently passed. The final artifact includes explicit scale-relative KKT residuals.

Detailed evidence: `docs/R17_JOINT_Q_RECONSTRUCTION_PILOT_20261110.md`.

**Next:** R18 should localize the 447 R16 source pressure-secant decreases by exact source row/density, matching metadata/bridge construction, and available phase/composition provenance. Preserve raw CompOSE package, normalized matched table, and projected solver input as distinct identity layers. Do not scale joint-q fitting or rerun TOV until a distinct reconstruction has a source-backed mechanism and a meaningful Pareto advantage. MD-31/32 and full-star radial/Cowling modes remain blocked.
