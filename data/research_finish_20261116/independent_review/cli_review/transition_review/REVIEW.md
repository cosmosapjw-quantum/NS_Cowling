# Independent transition-report semantics review

Verdict: **PASS_SCOPED** for the inspected diagnostic code and four-model preview arithmetic. This is not admission of the final all-model ledger or native exports.

The review covers `summarize_gppva_geometry.py`, `refine_gppva_geometry.py`, `export_gppva_geometry.py`, and `tf_morphology.py`; exact source and preview-input hashes are in `TRANSITION_SEMANTICS_REVIEW.json`. Source files were unchanged during the checks.

The retained preview has four models and five raw fine-grid sign-change brackets. Their differences, measured N192-to-N384 changes and `64*eps*max(1,|f_0|,|f_1|)` scales were independently recomputed from retained values. A synthetic one-ULP sign flip correctly remains roundoff-near-degenerate even when the two meshes give the same value. A missing intermediate branch prevents a sign bracket across that hole. A fine-grid sign smaller than its measured mesh change is flagged accordingly.

The sign flags mean only that a number exceeds the stated roundoff scale and/or the observed two-mesh change. The latter is not an error bound. The KKT/field/radius tolerances have not been converted into an energy-error budget, so neither flag certifies ordering against optimizer error. Even both flags together do not prove a continuum transition, a common tangent or a global minimum. The inspected code and geometry documentation retain those distinctions. The local density-shift estimate `n^2*|delta(Delta f)|/|Delta P|` has the correct envelope-slope units and is explicitly labelled neither an evaluated root nor an error bound.

Synthetic profile checks separately return uniform, near-uniform composition, noncanonical cylindrical shell and noncanonical spherical shell labels. Noncanonical shapes are excluded only from named simple-morphology pairs; their candidate energies are retained in the source records and complete-family selection. Near-uniform composition is not promoted to an independently resolved new phase.

A read-only check of the chosen repeated planar branches at DD2 n=0.07 and TW n=0.075 found that actual grid spacing R/N decreases from about0.1126fm to0.0563fm between N192 and N384. Their small reported mesh changes are therefore not explained by accidentally comparing a half-cell and full-cell at identical spacing in this preview.

The branch-envelope Simpson checks remain finite density-quadrature diagnostics. Their existing text permits truncation and within-label branch switches; they are not first-law admission tests or interpolation authorization. Final ledger completeness, selected-state provenance and native composition are reviewed separately by the parent reviewer. No solver or production file was modified here.
