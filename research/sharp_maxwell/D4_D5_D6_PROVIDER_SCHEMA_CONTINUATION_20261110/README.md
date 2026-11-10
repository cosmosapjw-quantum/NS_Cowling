# C7X D4-D6 provider/schema continuation

Date: 2026-11-10. This is an additive, evidence-only continuation of the C7X Sharp-Maxwell research DAG. It does not modify the source-authority branch or promote any EOS source.

## Finding

Thirteen selected `eos.compo` files (4,975 rows) parse structurally with zero malformed rows and zero row-index mismatches. Across the 11 GPPVA/SPG crust files, all 2,387 outer-crust rows show a mass-fraction-like closure diagnostic, while all 1,151 inner-crust rows have the final nuclear-group field exactly equal to `Y_e`. Neither the CompOSE-standard number-fraction closure nor a simple mass-fraction reinterpretation closes the inner-crust rows. CompOSE v3.01 defines the fourth nuclear-group field as the combined group number fraction, not a mass fraction. No inspected provider sheet documents a phase-dependent alternate encoding. D4 remains blocked; no remapping is authorized.

D5 retains the boundary-delta correction: 5,535 topology-preserving cells, 1,153 causal-polynomial violations, 1,110 cells with `dn/dmu_E <= 0`, zero nonpositive interpolated-density cells, 4,382 passing the selected interval tests. This is not a full replay of every retained cell: 15 newly discovered discrete-state boundary cells were independently evaluated and excluded from the prior 5,550-cell classification.

D6 stops before root search: 55/55 pairs have nonempty tabulated chemical-potential support overlap and require rest-mass rebasing; 0 pairs have both sides passing topology-aware D5, and all 55 retain the crust D4 blocker. No pressure equality, density drift, Maxwell root, TOV, or Cowling result was computed.

## Reproduction and limitations

The Dropbox backup ZIP includes the D4/D5 scripts, synthetic tests, receipts, CSV/JSON evidence, and SHA-256 manifest. Synthetic tests pass 11/11. No native Fortran build, local raw-source full D5 replay, provider contact, or D6-D10 run was performed. See `DAG_FINAL_STATUS.md`, `D4_PROVIDER_SEMANTICS_REASSESSMENT.md`, and `D5_TOPOLOGY_AWARE_UPDATE.md` in this directory.