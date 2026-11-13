# Shared background note (DRAFT_NOT_SENT)

We are validating CompOSE crust tables for use in a neutron-star/Cowling-mode study. Our reader follows CompOSE Reference Manual v3.01 (record layout Eq. 4.19; Y_i = n_i/n_b; for a nuclear group Y_I = sum_j Y_j and X_h = Aav*Y_I). Raw files and parsing are verified: the `eos.compo` keys match `eos.thermo` row by row, and every row has four particle pairs (10, 11, 0, 1) and one nuclear-group quadruple (group index 1).

The CompOSE data sheets define Iphase but not the contents of the quadruple's fourth field, and the public papers contain no composition table. We therefore cannot tell whether the stored value is the CompOSE number fraction Y_I, a baryon (mass) fraction, or a provider-specific quantity, nor whether Aav/Zav refer to the nuclear cluster or the whole Wigner-Seitz cell. Our purpose is only to understand the export convention before using the data; we are not proposing that any value is wrong, and we have not modified the files.

Attachments: `evidence_summary_onepage.md`, `sample_rows_<SOURCE>.csv`, `aggregate_statistics.json`. The full row-level audit is reproducible from the pinned blobs and is public in the project repository (NS_Cowling, PR #4, directory `artifacts/sharp_maxwell/C6B_COMPOSITION_SEMANTICS_20261110T100954Z`).
