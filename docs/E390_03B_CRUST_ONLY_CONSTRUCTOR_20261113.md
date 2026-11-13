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
- All row composition values, row statuses, presence masks and C3/C4 applicability are copied. An absent abundance is not replaced by zero.
- The source thermodynamic rows carry their admitted CLD/Thomas--Fermi provenance. E390-03B does not run a new CLD/Thomas--Fermi minimization, and it does not solve or reconstruct composition; that belongs to E390-03.
- C1--C11 physical diagnostics are recorded per model, separately from construction/readback. No clipping, smoothing, renormalization or thermodynamic repair is allowed.

The per-model `E390-03B_RECEIPT.json` records source hashes, method provenance, exact support, categorical counts, output hashes and physical diagnostics. The root `crust_coverage_report.v1.json` distinguishes drip, provider splice, inner-crust support and GPPVA pasta geometry. The root `crust_only_manifest.v1.json` reports the admitted 12-model cohort and its total of 3,883 source rows.

## Non-claims and downstream gates

This node does not execute composition solvers, core-crust matching, any of the 390 smooth method pairs, the separate sharp-Maxwell workflow, TOV/Cowling integration or mode identification. GPPVA_TW remains blocked. The five core-only tables have independent C2/C10 physical diagnostic failures already recorded in the E390-03A receipt; those failures continue to block the affected core models from matched production use. Construction/readback is not a physical-admissibility certificate.

The next DAG node after successful E390-03B closure is **E390-03: composition solvers/reconstruction**. Matching-method integration is downstream of that node and must be executed as its own later DAG work.
