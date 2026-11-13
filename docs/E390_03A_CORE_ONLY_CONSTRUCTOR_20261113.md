# E390-03A: native standalone core-only EOS construction

**Scope.** This node constructs five standalone core-only CompOSE tables from source rows admitted by E390-01C. It does not run a crust builder, solve composition, perform crust-core matching, or regenerate the 390-pair matrix.

## Source-authoritative row admission

All row ranges are 1-based positions in the immutable eos.nb data grid.

| Model | Source rows | Core-only rows | Output rows | Additional guard |
|---|---:|---:|---:|---|
| APR_APR | 1–172 | 70–172 | 103 | Keep the 0.075959–0.085 fm^-3 source gap |
| BBB_BHF2 | 1–84 | 70–84 | 15 | Keep the 0.075959–0.16 fm^-3 source gap |
| GM_GM1 | 1–349 | 49–349 | 301 | Iphase=3; do not synthesize the 1e-3 fm^-3 cutoff row |
| RG_KDE0v1 | 1–1240 | 967–1240 | 274 | Model-local Iphase=3 homogeneous matter |
| RG_SLY4 | 1–197 | 144–197 | 54 | Model-local Iphase=3 homogeneous matter |

The exact eos.nb Git blob SHA-1 values are pinned in the source registry in api.py and the E390-01C source-adjudication receipt. The Python frontend verifies that source identity before it asks the native backend to read the table. This hash is a Git blob identity, not a substitute for the independent source-archive SHA-256.

## Numerical and data contract

The constructor performs only an exact row projection. No EOS interpolation, density-gap fill, smoothing, thermodynamic re-normalization, phase remapping, composition reconstruction, or extrapolation is allowed.

The thermodynamic representation follows CompOSE:
\[
P=n_B Q_1,\qquad
\mu_B=m_n(1+Q_3),\qquad
F/n_B=m_n(1+Q_6),\qquad
\epsilon=n_Bm_n(1+Q_7).
\]
Here n_B is in fm^-3, Q1 and m_n,m_p are in MeV, and epsilon and P are in MeV fm^-3. Q2, Q4, Q5, all Nadd values, and source mass/lepton metadata are copied unchanged. These identities are the CompOSE representation, not quantities recomputed by this builder.

Iphase remains an EOS-local categorical code. For GM1, KDE0v1 and SLy4, the admitted core rows must carry the source-adjudicated code 3. APR and BBB-BHF2 retain the phase labels actually present in their admitted source rows. No cross-model phase dictionary is inferred.

The following are preserved row by row: particle/group values, Iphase, row_status, all per-species/group presence masks, and independent C3/C4 applicability. A written zero remains present; an absent pair remains absent. The CompOSE eos.t and eos.yq source index bounds are retained when writing the one-point cold slice. The original eos.nb density-row labels are retained by shifting inb_start by first_row-1.

## Implementation

- Native Fortran extract_core_rows validates source point count, exact row range, cold single-slice scope, and model-specific phase guard, then slices all row-bearing arrays/tensors.
- Native Fortran extract-core writes the new CompOSE directory and immediately reads it back. Exact thermodynamic, phase, composition-presence, tensor, axis, and density-row-key parity is required.
- Python build_core_only is limited to source identity checks, paths, native invocation, readback comparisons, SHA-256 receipts, and presentation. There is no Python numerical fallback.
- Output roots are immutable: an existing attempt path is rejected. Use a fresh versioned directory for every new run.
- eos.micro and eos.mr are not handled by the current EOSData reader/writer and are explicitly not copied by this node.

## Verification and evidence semantics

The CI regression covers all five admitted source models, exact output row counts, strict density monotonicity, native readback, and fail-closed source-hash/overwrite controls. Each production execution writes a model receipt and a root manifest with input/output file SHA-256 values.

A successful extraction/readback establishes implementation-verified row preservation for the declared models. It does not establish physical admissibility. C1-C11 diagnostics are recorded separately; any non-PASS first-law, stability, causality, beta-equilibrium or composition-closure result remains a downstream blocker and is never repaired here.

The CI uses temporary test directories. Production component tables still need to be generated into a new versioned project output directory on the authorized source-data workspace and backed up there before E390-03A can be marked fully complete.
