# E390-03B: native standalone crust-only EOS construction

## Scope and admission

This node constructs the twelve crust sources admitted by E390-01C:

- SPG_M1--SPG_M5: BSk22 outer-crust source plus model-specific compressible liquid-drop (CLD) inner crust, as documented by Scurto, Pais & Gulminelli, *Phys. Rev. D* **109**, 103015 (2024), arXiv:2402.15548.
- GPPVA_DD2, DDME2, FSU2, FSU2H, FSU2R, NL3wrL55 and TM1e: BSk22 outer-crust source plus model-specific Thomas--Fermi inner crust, with the model papers/datasheets registered in `data/eos390_source_bibliography.v1.json`; the common method paper is Grill et al., *Phys. Rev. C* **90**, 045803 (2014), arXiv:1404.2753.

GPPVA_TW is excluded. The current official CompOSE page 212 states BPS below 0.002 fm^-3, while its linked PDF states BSk22. This unresolved source conflict is not resolved by analogy with the other GPPVA tables.

## Exact source support and boundaries

All rows are 1-based positions in the immutable `eos.nb` source grid. The admitted inputs contain:

| Regime | Source rows | Meaning |
|---|---:|---|
| Outer crust | 1--217 | BSk22 outer-crust rows, `Iphase=0` |
| Neutron drip | 217 -> 218 | Physical BSk22 drip boundary near 2.69e-4 fm^-3 |
| Inner crust / cluster regime | 218--N | Model-specific CLD or Thomas--Fermi rows, `Iphase=7` |
| Provider splice | nominal nB=2e-3 fm^-3 | Separate BSk22-to-model-inner-crust source splice, not neutron drip |

The sampled source grid jumps from the row-217 drip marker to the first inner-crust row at or above 2e-3 fm^-3. No row is synthesized to fill the intervening density interval. Each model terminates at its own exact source endpoint, ranging from approximately 0.054 to 0.0977 fm^-3. No homogeneous high-density core tail is appended.

The physical neutron-drip marker is supported by Pearson et al., *MNRAS* **481**, 2994 (2018), Table 7, and the BSk22 CompOSE table. The source adjudication and exact model row counts/hashes remain in the E390-01C receipt.

## Native construction and lossless data contract

The native Fortran `extract-crust` operation reads the original CompOSE source, checks the source point count and phase split, verifies the BSk22 drip marker, then projects the complete source row support without transformation. The existing native CompOSE writer writes the new directory; a fresh native reader immediately reads it back, and the exact round-trip verifier checks thermodynamic arrays, Q-extra, phase, composition values, masks, C3/C4 applicability, axis metadata and row labels.

- `eos.nb`, `eos.t`, `eos.yq`, `eos.thermo` and `eos.compo` are retained and hashed.
- SPG source tables have Nadd=0. GPPVA tables have Nadd=1; their `eos.thermo:q_extra` values are categorical geometry labels (the observed source codes are 0, 1, 2 and 3).
- `eos.compo:Iphase` and `eos.thermo:q_extra` are distinct fields. They are copied verbatim and never merged, interpolated or numerically remapped.
- All raw eos.compo particle/group values, Iphase, presence masks, and the reader's derived C3/C4 applicability/status are compared across native readback. The source file does not serialize a `row_status` column: `EOSData.row_status` is computed by the native reader under the standard CompOSE group mapping. It is therefore recorded as a reader-derived diagnostic, not attributed to the provider.
- The standard CompOSE convention is retained: `YI` is the combined nuclear-group number fraction and `Xh=Abar*YI` is the group baryon fraction. The native reader maps `Yh=Abar*YI`; the writer maps `YI=Yh/Abar`. These transformations are not changed in E390-03B.
- The E390-02 source audit is pinned by Git blob SHA `419c00b3bf75bf13a41a58ce8edd99b91caaafaf`. It found that for `Iphase=0` rows, the raw group field is closure-compatible with `Xh`, not standard `YI`; for `Iphase=7` rows the raw group field is compatible with neither simple assignment using the available n/p/e/mu fractions. The audit therefore leaves provider semantics unresolved and states that normalized closure is a reconstruction, not independent validation of the raw provider data.
- Accordingly, the all-row `row_status=FAIL` values reported for these files are **reader-derived under standard CompOSE interpretation**, not source-stored row labels. The dedicated `extract-crust` operation uses a narrowly scoped native-writer opt-in to serialize finite raw source values and masks without repairing them; the ordinary CompOSE writer remains fail-closed by default. C3/C4 `NOT_AVAILABLE` is not a closure PASS.
- E390-03B does not run a new CLD/Thomas--Fermi minimization or reconstruct group composition. That is reserved for E390-03, which must preserve the raw parent records, declare the provider-specific mapping/reconstruction, and attach row-level provenance. Until that node resolves the composition semantics, the raw provider composition must not be used as validated composition in matched products.
- C1--C11 diagnostics are recorded in full per model, including the applicability/status of C3 and C4. Thermodynamic checks that pass do not override unavailable or unresolved composition semantics. No clipping, smoothing, renormalization or thermodynamic repair is allowed.

The per-model `E390-03B_RECEIPT.json` records source hashes, method provenance, exact support, categorical counts, output hashes and physical diagnostics. The root `crust_coverage_report.v1.json` distinguishes drip, provider splice, inner-crust support and GPPVA pasta geometry. The root `crust_only_manifest.v1.json` reports the admitted 12-model cohort and its total of 3,883 source rows.

## Non-claims and downstream gates

This node does not execute composition solvers, core-crust matching, any of the 390 smooth method pairs, the separate sharp-Maxwell workflow, TOV/Cowling integration or mode identification. Composition semantics remain BLOCKED for matching until E390-03 completes. GPPVA_TW remains blocked. The five core-only tables have independent C2/C10 physical diagnostic failures already recorded in the E390-03A receipt; those failures continue to block the affected core models from matched production use. Construction/readback is not a physical-admissibility certificate.

The next DAG node after successful E390-03B closure is **E390-03: composition solvers/reconstruction**. Matching-method integration is downstream of that node and must be executed as its own later DAG work.
