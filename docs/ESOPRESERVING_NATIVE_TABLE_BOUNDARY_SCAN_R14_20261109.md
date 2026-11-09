# R14: EOS-preserving native-table boundary scan

**Status:** implementation staged for exact-head GitHub Actions verification. No physical P=0 surface or crust model is claimed. The scan also exports the native TOV background profile for every finite-pressure endpoint, to preserve provenance for a separately gated radial-stability/Cowling continuation.

## Research decision

The primary analysis will use each supplied projected EOS CSV exactly as stored. No EOS row is appended, removed, reordered, smoothed, clipped, fitted, or rewritten. The old R13 synthetic low-density continuations remain historical diagnostic artifacts only; they are not admissible input for the current primary analysis and are not used to infer a physical crust or P=0 surface.

The continuous TOV provider remains the frozen v5 native-knot cubic-Hermite representation in P(mu_b), using the source P, n_b, and mu_b knot values. This is interpolation between the existing knots, not a replacement EOS table. The input CSV is passed directly to the solver and its SHA-256 is recorded before and after every method scan.

## Frozen cohort and raw endpoint

The cohort remains GPPVA_FSU2 x BBB_BHF2 with six matching methods: `css`, `cubic_F`, `cubic_P`, `maxwell`, `quintic_F`, and `quintic_P`.

Each projected EOS file contains 323 rows. The six tables share the same first native row:

- source row: 1
- P_min = 9.000000079422e-12 MeV fm^-3 = 1.4419589833248072e21 Pa
- epsilon_min = 4.354134734530827e-7 MeV fm^-3
- n_b,min = 4.67962e-10 fm^-3
- mu_b,min = 930.4654511545 MeV

The first 323 rows remain in their supplied order. The pressure selector rejects non-monotone input rather than sorting it. Each cutoff is an exact pressure value copied from an existing source row.

## Cutoff design

Eight target offsets from the native table-edge pressure are used only to select nearby existing knots. The selected source rows are the same across the six files:

| Target offset (dex) | Source row | Actual offset (dex) | P_cut (MeV fm^-3) |
| ---: | ---: | ---: | ---: |
| 0.00 | 1 | 0.000000 | 9.00000007942e-12 |
| 0.25 | 3 | 0.311229 | 1.84277104749e-11 |
| 0.50 | 4 | 0.467474 | 2.64068378497e-11 |
| 0.75 | 6 | 0.779805 | 5.42059837899e-11 |
| 1.00 | 7 | 0.936080 | 7.76822987837e-11 |
| 1.50 | 11 | 1.560687 | 3.27287616375e-10 |
| 2.00 | 15 | 1.994697 | 8.89078117305e-10 |
| 3.00 | 20 | 3.049259 | 1.00809578095e-08 |

At each cutoff the central pressure is fixed to the frozen six-method reference value, and the existing v5 TOV solver is run. This produces 8 x 6 = 48 TOV integrations.

## Interpretation guard

A run stopped at P_cut > P_min ends on an interior isobar. Its reported M(P_cut) is the mass enclosed by that isobar and R(P_cut) is its coordinate radius. These are truncation diagnostics, not a complete neutron-star mass and radius. The row at P_min is the low-pressure edge of the supplied table, but P_min is still strictly positive; therefore even that run is not mathematically a P=0 surface solution.

No crust or vacuum EOS is inferred below P_min. A physical full-surface radius requires a sourced and specified low-density EOS or a separately labelled approximation. That is outside this EOS-preserving run.

## Acceptance gates

1. The six exact source CSVs are read in their original order and their SHA-256 hashes are unchanged after all runs.
2. The selected cutoff pressures are values from source rows, never interpolated targets.
3. The native table-edge (source row 1) TOV runs reproduce the frozen mass and radius to relative tolerance 2e-5 for all six methods.
4. All 48 requested integrations complete with finite positive enclosed mass and isobar radius.
5. The artifact contains JSON/CSV receipts and native TOV background profiles emitted by the frozen v5 solver; no synthetic EOS table is emitted. Every background CSV has a SHA-256 and row-count receipt.
6. Results are described as finite-pressure isobar truncations. No P=0, crust, radial stability, or Cowling-spectrum claim is inferred from this TOV-only scan.

## Research provenance and context

- O2scl's tabulated-EOS TOV example explicitly installs a default low-density crust EOS before reading the core/table data. This is a separate physical completion step, not evidence that a table endpoint is a vacuum surface: [O2sclpy TOV tabulated-EOS example](https://awsteiner.org/code/o2sclpy/ex_tov.html).
- Chamel & Haensel review the distinct outer- and inner-crust compositions and their structural and dynamical roles: [Physics of Neutron Star Crusts](https://link.springer.com/article/10.12942/lrr-2008-10).
- Fortin et al. quantify how non-unified crust-core matching can produce radius and crust-thickness uncertainty. Their quoted scale (up to about 4% in radius and 30% in crust thickness) concerns crust-core matching at much higher density, and must not be transferred directly to the extremely low-pressure table-edge truncation studied here: [arXiv:1604.01944](https://arxiv.org/abs/1604.01944).
- Davis et al. (2024) and related crust-connection studies motivate keeping matching prescription uncertainty distinct from EOS-table interpolation and boundary truncation: [SciSpace record](https://scispace.com/papers/inference-of-neutron-star-properties-with-unified-crust-core-lkf82hwmxw).
- The 2026 thin-crust approximation study compares approximate crust reconstructions with complete TOV solutions and explicitly uses core-crust interface information. It does not justify appending an unsourced pressure tail to the present table: [Köpp et al. (2026)](https://link.springer.com/article/10.1140/epjc/s10052-026-15443-3).

## Next gate after the TOV scan

Only after the 48-case result and source hashes pass will the project decide whether to run a corresponding finite-boundary Cowling sensitivity study. Any such run must be labelled as a truncated-domain eigenproblem until the physical surface boundary condition below the native EOS endpoint is independently specified. The current TOV scan does not by itself establish a full-star oscillation spectrum. Background profiles are exported only as finite-domain integration products; they do not provide the missing physical EOS below the table edge or a justified P=0 oscillation boundary condition.
