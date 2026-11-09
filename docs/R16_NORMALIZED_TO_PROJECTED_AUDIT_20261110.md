# R16 full-cohort normalized-to-projected EOS and convexity audit

**Status:** `PASS_R16_FULL_COHORT_PROJECTION_RECONSTRUCTION` (implementation and artifact acceptance only)  
**GitHub Actions:** [42666868668](https://github.com/cosmosapjw-quantum/NS_Cowling/actions/runs/42666868668)  
**Code HEAD:** `858c5d5639e0140e53c6f27b6e504eb261c05b01`  
**Frozen spectrum-summary SHA-256:** `ce54c111510f32086065f95356903cd6abad675adce12cab1a06b5cc82350965`  
**Artifact ID:** `13697021829`  
**Artifact SHA-256:** `ad9d9858e15ccb013a9591ef63d8e36528202a6ef9d237284765f05c5bbcc8e1`

## 1. Scope and identity boundary

This audit compares the 234 frozen `eos_original_normalized` CSVs with the 234 frozen `eos_projected` solver-input CSVs named in `spectrum_summary_234.csv`. These are the normalized matched EOS source tables and their derived solver tables. They are **not** a claim that the repository's separately stored raw CompOSE package files are byte-identical to the original files recorded in `original_compose_hashes_json`.

The audit did not modify either source directory, create synthetic EOS rows, run TOV, or calculate oscillation modes. It uses the existing six matching-method labels (`css`, `cubic_F`, `cubic_P`, `maxwell`, `quintic_F`, `quintic_P`); the historical `maxwell` label remains a smooth bridge prescription in this cohort, not a sharp first-order Maxwell interface.

## 2. Exact projection reproduced

The frozen projector implements weighted isotonic regression of

\[
q(\mu_B)=\frac{n_B}{\mu_B}
\]

on trapezoidal chemical-potential weights, then sets

\[
n_B(\mu_B)=\mu_B q(\mu_B),\qquad
P(\mu_B)=P_0+\int_{\mu_0}^{\mu_B}n_B(\tilde\mu)\,d\tilde\mu,
\qquad
\epsilon=\mu_B n_B-P.
\]

Between retained knots, `q` is linear and the pressure integral is evaluated exactly. The source data have `q` nondecreasing in all 234 tables, so PAV performs no pooling: 0 pooled blocks, 0 pooled knots, and 0 knot reductions. The main numerical change is therefore not an alteration of `n_B` by isotonic pooling. It is the reconstruction of pressure from the differential relation `dP/dmu_B = n_B`, followed by recomputation of `epsilon` from the First-Law identity.

## 3. Full-cohort verification

| Check | Result |
|---|---:|
| Frozen normalized source tables | 234 / 234 |
| Frozen projected tables | 234 / 234 |
| Independent projector reconstructions matching projected CSVs | 234 / 234 |
| Maximum relative reconstruction mismatch | 0.0 |
| Normalized/projected source hashes unchanged | 234 / 234 each |
| Source algebraic First-Law residual, maximum | `5.09915761985166e-16` |
| Projected algebraic First-Law residual, maximum | `2.1628977905212927e-16` |
| Source `q` nonmonotone tables | 0 / 234 |
| PAV-pooled tables / knots | 0 / 0 |
| Projected Hermite cells failing stability or causality | 0 / 109,776 |
| TOV integrations / synthetic EOS rows in R16 | 0 / 0 |

The artifact manifest's six data-file SHA-256 values were recomputed locally and all matched. The downloaded ZIP SHA-256 matched the GitHub artifact digest, and `unzip -t` reported no archive errors.

## 4. Raw-source versus projected Hermite-cell admissibility

For each adjacent pair of knots, the cubic-Hermite `P(mu_B)` cell is tested using its endpoint second derivatives and causal margins. The production whole-cell conditions are

\[
P''(\mu_B)>0,\qquad G(\mu_B)=\mu_B P''(\mu_B)-P'(\mu_B)\ge0.
\]

Because `P''` is linear and `G' = mu_B P'''` has fixed sign for a cubic cell at positive chemical potential, checking both endpoints is sufficient for these conditions on the complete cell. The projected representation passes this test for all 109,776 cells. The source-normalized representation does not:

| Source-normalized cells | Count | Fraction |
|---|---:|---:|
| `P''` not strictly positive at both endpoints | 24,465 / 109,776 | 22.2863% |
| Stability or causal-margin condition rejected | 24,505 / 109,776 | 22.3227% |
| Stable but causal-margin rejected | 40 / 109,776 | 0.0364% |
| Projected cells rejected | 0 / 109,776 | 0% |

Every one of the 234 source-normalized tables contains at least one nonstable and one rejected cell. In the specific source-row-217→218 cell, 144/234 source tables fail, while 0/234 projected tables fail. All 8 GPPVA crust families fail in all 18 corresponding core/method tables; all 5 SPG crust families pass in all 18 corresponding tables.

The identity `epsilon + P = mu_B n_B` is an algebraic First-Law closure check at the stored knots. It does **not**, by itself, establish the differential compatibility `dP/dmu_B = n_B`; the cell tests and secant analysis probe that separate requirement.

## 5. Discrete pressure-knot convexity obstruction

Define the adjacent pressure secant slopes

\[
s_i=\frac{P_{i+1}-P_i}{\mu_{i+1}-\mu_i}.
\]

For a convex function `P(mu_B)`, the secant slopes over successive adjacent intervals must be nondecreasing. The source-normalized tables contain 447 downward secant transitions (counting relative decreases beyond `1e-12`), distributed over all 234 tables. The projected tables contain zero downward transitions.

| Statistic | Result |
|---|---:|
| Source downward secant transitions | 447 |
| Source tables with at least one transition | 234 / 234 |
| Projected downward secant transitions | 0 |
| Worst source relative secant drop | -21.3709% |
| Transitions with drop at least 1% | 125 |
| Transitions with drop at least 5% | 25 |
| Transitions with drop at least 10% | 18 |

By matching-method label, the 447 transitions are: `maxwell` 177; each of `css`, `cubic_F`, `cubic_P`, `quintic_F`, and `quintic_P` 54. The largest single relative drop is in `GPPVA_TM1e`, at the source-row-290 knot, with geometric-mean density approximately `0.0754851 fm^-3`.

This gives a necessary-condition obstruction: **no convex interpolant can pass through every source-normalized pressure knot exactly when its adjacent secant slopes decrease**. Changing the cubic-Hermite slope limiter alone cannot remove this obstruction while preserving all of those pressure values. At least some thermodynamic coordinates or the pressure-knot values must change in a convex reconstruction.

This is a mathematical statement about the stored normalized pressure knots. It is not a statement that the original microphysical model itself is physically unstable; the origin of the knot mismatch and its relation to upstream CompOSE records remain separate provenance questions.

## 6. Density localization of the source-normalized cell failures

Cells are binned by the geometric-mean density of their two endpoints. The operational density gap from the last low-density source row to the first ineligible row appears as an empty `1e-3–2e-3 fm^-3` bin.

| Geometric midpoint density (`fm^-3`) | Cells | Source nonstable | Source rejected by stability/causality | Projected rejected |
|---|---:|---:|---:|---:|
| `< 1e-4` | 33,462 | 4,842 (14.47%) | 4,842 (14.47%) | 0 |
| `1e-4–1e-3` | 17,316 | 4,590 (26.51%) | 4,590 (26.51%) | 0 |
| `1e-3–2e-3` | 0 | 0 | 0 | 0 |
| `2e-3–1e-2` | 6,516 | 5,400 (82.87%) | 5,400 (82.87%) | 0 |
| `1e-2–1e-1` | 20,745 | 9,147 (44.09%) | 9,183 (44.27%) | 0 |
| `1e-1–5e-1` | 16,371 | 486 (2.97%) | 490 (2.99%) | 0 |
| `>= 5e-1` | 15,366 | 0 | 0 | 0 |

The failures are not confined to the lowest-density prefix. The `2e-3–1e-2` band is particularly problematic for this source-normalized cohort, and substantial failure fractions also occur in `1e-2–1e-1`. This is a representation diagnostic grouped by density, not a phase classifier; source phase labels and exact provenance should be compared before assigning a physical interpretation to each bin.

## 7. Projection response at the sparse row-217→218 boundary

Across the 234 tables, row-218 pressure is changed by `+11.5416%` to `+64.0802%` relative to source-normalized pressure. The row-217→218 pressure increment is increased by `+13.8659%` to `+85.9332%`. This is consistent with the fact that pressure is re-integrated from the preserved `n_B(mu_B)` data rather than copied from source pressure knots. It should not be interpreted as an outer-crust thickness or physical surface correction.

## 8. Interpretation and next gate

**Established by this run:** the frozen projector is reproduced exactly over all 234 inputs; its projected tables satisfy the declared whole-cell Hermite stability/causality checks; the source-normalized pressure knots contain discrete nonconvexity that prevents an exactly pressure-knot-preserving convex interpolant.

**Not established:** that the projection is the unique or physically preferred EOS reconstruction; that the input source knots are microphysically wrong; physical completeness of an outer crust; a `P=0` surface; radial stability; or full-star Cowling modes. R16 did not rerun TOV, compute modes, or change the previous R15 results.

### Proposed R17 objective

Compare explicitly defined thermodynamic reconstructions under a transparent data-fidelity objective and the constraints `P'>0`, `P''>0`, `mu P''-P'>=0`, and `epsilon=mu n-P`. The current q-preserving exact-integration projector should be retained as baseline A. Any alternative that relaxes pressure, density, or chemical-potential knots must declare which variables are allowed to move and how weights are chosen. Do not infer statistical uncertainties or invent weights from the current files; where uncertainty metadata are absent, treat weights as sensitivity parameters and report a Pareto/error budget rather than a single “best” correction. Keep every alternative in a separate derived EOS artifact with a parent-source hash; do not overwrite the frozen primary projected cohort.

Before TOV/mode sensitivity runs, map source-normalized secant violations and Hermite failures against upstream EOS source identities, CompOSE phase/composition records, and the matching-method transformation. Only then select a bounded set of independently justified alternatives for structural and spectral propagation.

MD-31 (source-backed `P=0` crust completion) and MD-32 (a justified finite-domain perturbation boundary) remain separate. Radial `omega_0^2` and complete-star Cowling modes remain blocked.

## 9. Artifact provenance

- GitHub Actions run: [42666868668](https://github.com/cosmosapjw-quantum/NS_Cowling/actions/runs/42666868668)
- Artifact ID: `13697021829`
- ZIP SHA-256: `ad9d9858e15ccb013a9591ef63d8e36528202a6ef9d237284765f05c5bbcc8e1`
- Frozen summary SHA-256: `ce54c111510f32086065f95356903cd6abad675adce12cab1a06b5cc82350965`
- Receipt SHA-256: `7c4fed99c0ed25a976dae170873c05b01fdb9ac64eeaa094ec98121e85d3828a`
- Table CSV SHA-256: `d46698d02af727dce47ab0c947a0f47b30fa4898af0cf2929c039818025f6389`
- Group CSV SHA-256: `c5c0986e697f9a4bd31712b1036d27b5175c1e5ef05314464848421cf0cfb4fd`
- Cell CSV SHA-256: `e4b8bd2eba9e7c9dce2f91129ed946e222d4bd87620d1adeebed89e6f373a739`
- Density-bin CSV SHA-256: `556b8365852b4b7cf07f69b0a6133f5cc9f4f9437c0ed7e6f884f893c57d38e3`
- Convexity-violation CSV SHA-256: `a4947a82c069bc6eb23bca91f695323bdb08c5e6fef7abb5b3d05f845c308de1`
