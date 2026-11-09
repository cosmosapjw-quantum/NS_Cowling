# Mass-dependent Cowling research DAG (2026-11-09)

## Authority

- Repository: cosmosapjw-quantum/NS_Cowling
- Research branch: research/mass-dependent-cowling-v1
- Base SHA at branch creation: a09ec30d74b363652d2d1aefe990557c468679dd
- PR #3 remains open/draft/unmerged and is not modified.
- Frozen smooth cohort: 13 crusts × 3 admitted cores × 6 smooth methods = 234 rows / 39 paired core-crust combinations.
- Old spectra have only the fixed target mass 1.3300106263 solar masses; reuse audit says zero modes were recomputed.
- Historical method token maxwell is a smooth bridge, not a sharp density jump.

## Dependency graph

\[
\text{frozen cohort and provenance}
\to TOV\ M(P_c)
\to \text{connected rising branch}
\to \text{adjacent-}P_c\text{ target bracket}
\to \text{native mass-root solve}
\to \text{Cowling roots and W/V identity}
\to \text{six-method paired maps}
\to \text{refinement and benchmarks}
\to \text{publication decision}.
\]

Mechanism and inference extensions branch from the verified mass-dependent spectrum. They are not silently included in the present implementation.

## Node ledger

| ID | Node | Status | Exit evidence |
|---|---|---|---|
| MD-00 | Freeze 234-row cohort, projected EOS tables, release summary and checksums | PASS_READ_ONLY | 13 crusts, 3 core families, 6 smooth methods. All old spectra at 1.3300106 solar masses. |
| MD-01 | Stable-branch extractor and paired estimands | IMPLEMENTED_CI_PENDING | Tests cover first turning point, failed-sample barriers, adjacent target brackets and no-imputation six-method pairing. |
| MD-02 | Exact-backend TOV-only sequence and explicit numerical limits | IMPLEMENTED_CI_PENDING | `solve_cowling_tov_star` calls `cell_solve task=star` with v5 `r0=10 m`, `rtol=1e-9`, `step_cap=0.05`; configurable radius limit defaults to historical 1e6 m and is set to 1e8 m for the exploratory scan. Explicit Pc brackets remain optional and preserve the old default. |
| MD-03 | Six-model high-tail pilot: BBB_BHF2 × GPPVA_FSU2, all six methods and eight masses | WAITING_FOR_PUSH_ACTIONS | Pilot validator requires frozen baseline TOV reproduction and complete six-method 1.0- and 1.33-solar-mass pairs; partial valid pairs are rejected. |
| MD-04 | Full 234-EOS × eight-mass execution | NOT_RUN | Manual workflow dispatch with cohort=full; 48 log-pressure TOV samples per EOS. |
| MD-05 | Mass-response and coverage analysis across 39 fixed pairs | NOT_RUN | Complete-pair denominators, paired ranges, per-core/per-crust maps, method sign reversals and peak mass locations. |
| MD-06 | ODE/step/root refinement versus effect size | NOT_RUN | Frequency movement no larger than 0.1 of the smallest reported paired effect. |
| MD-07 | Direct recognized matching-prescription baselines | NOT_STARTED | Simple concatenation as a non-physical control; linear P(epsilon); Canullán-Pascual thermodynamic P(epsilon) interpolation; CUTER-like/documented thermodynamic crust reconstruction. |
| MD-08 | Independent Cowling/full-GR benchmark | NOT_STARTED | Representative low-, reference- and high-mass points; Cowling systematic offsets separate from matching shifts. |
| MD-09 | Eigenfunction/acoustic-cavity mechanism | NOT_STARTED | Matching-region overlap and acoustic travel time; variational kernel only after boundary/matching derivation. |
| MD-10 | Maximum mass, crust thickness/mass and tidal Lambda | PARTIAL_BLOCKED | TOV M,R sequences are in scope. Resolved turning point is only an Mmax candidate; crust event data and Love-number solver are not in this increment. |
| MD-11 | Source-backed CompOSE EOS and NICER/tidal inference | BLOCKED_UPSTREAM | Requires source provenance, thermodynamic admission and structural/tidal observables. |
| MD-12 | Publication claim gate | NOT_READY | JKPS/EPJA after MD-03–08, direct baselines and cohort coverage; A&A/MNRAS after MD-10–11; PRD only with new theory/full-GR result. |

## Analysis contract

For y in {f,p1}:
\[
\Delta_m^{(y)}(M;c,s)=100\left[y_m(M;c,s)/\bar y_6(M;c,s)-1\right].
\]
Only complete six-method cells enter paired summaries. No imputation. Report complete-pair count at every mass. The 234 rows are a finite crossed design, not independent microscopic EOS samples or a posterior.

## Stop and reject conditions

- Reject a target mass not bracketed within the connected rising-M(Pc) branch.
- Do not label a truncated branch endpoint as Mmax.
- A failed TOV sample is a branch barrier; do not bridge it.
- Reject a spectrum row with failed mass closure, failed shooting residual or unresolved/wrong W/V identity.
- Reject incomplete six-method cells from paired estimands; do not impute.
- Do not infer published-method novelty from the six smooth bridges alone; direct matching-prescription baselines are a separate gate.
- A converged Cowling root is not a full-GR QNM and supplies no damping time.
- The legacy maxwell method label is a smooth bridge, not a sharp density-jump construction.
- Do not claim crust mass/thickness, tidal deformability, observational inference or a sensitivity kernel from the current increment.

## Immediate next action

Re-run the six-model pilot after the same-cell TOV-only backend, configurable radius limit and validator correction are pushed. Inspect baseline M/R parity and all sequence samples first. Do not dispatch the full cohort unless the pilot has valid 1.0- and 1.33-solar-mass paired cells and the M/R, mass-residual, shooting-residual and node gates pass. The first run has zero valid spectra and must not be cited as a scientific result.
