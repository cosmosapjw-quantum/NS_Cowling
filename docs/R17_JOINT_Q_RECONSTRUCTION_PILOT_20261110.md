# R17 constrained joint-q reconstruction pilot

**Status:** `PASS_R17_JOINT_Q_FIT_PILOT` for code/test execution and candidate admissibility on the declared four-table slice. The scientific outcome is a **negative/limited Pareto result**: this objective family does not materially repair the source pressure mismatch.

- GitHub Actions: [42669846775](https://github.com/cosmosapjw-quantum/NS_Cowling/actions/runs/42669846775)
- Code HEAD: `9fde7e2206a162530c03a77098c5fa7fa19eeba7`
- Artifact ID: `13697898562`
- Artifact SHA-256: `8a266c0b5c1316391fd6ede3d5ae65904b869a8dab66fcb29230a3ede5d52799`
- Frozen spectrum-summary SHA-256: `ce54c111510f32086065f95356903cd6abad675adce12cab1a06b5cc82350965`

## Scope

Four representative tables were used, all with BBB_BHF2 core:

- `GPPVA_DD2__BBB_BHF2__css`
- `GPPVA_TM1e__BBB_BHF2__css`
- `GPPVA_TM1e__BBB_BHF2__maxwell`
- `SPG_M1__BBB_BHF2__css`

The sensitivity grid is `lambda_density=(1, 0.1, 0.01, 0.001, 0.0001)` with two quadrature measures: trapezoidal in `mu_B` and trapezoidal in `log(n_B)`. The result contains 40 candidate EOS reconstructions and 14,740 knot-level records.

The candidate keeps the current thermodynamic representation:
[
q(\mu_B)=n_B/\mu_B\quad\text{piecewise linear},\qquad
n_B=\mu_B q,
]
[
P(\mu_B)=P_0+\int n_B\,d\mu_B,\qquad
\epsilon=\mu_B n_B-P.
]
The objective combines weighted relative pressure and density residuals. `lambda` and the quadrature coordinate are sensitivity choices, not inferred uncertainties or likelihood weights.

## Acceptance and independent verification

- Unit tests: pass.
- 40/40 candidates generated and marked admissible.
- `lambda=1` reproduces the frozen projected parent EOS; maximum relative mismatch is at floating-point roundoff (below (5\times10^{-12})).
- All 4 normalized source hashes and all 4 projected parent hashes remain unchanged.
- All candidates have positive, increasing thermodynamics and nondecreasing (q).
- 0 nonstable cells, 0 cells rejected by the causal-margin test, 0 candidate pressure-secant decreases.
- Maximum algebraic First-Law residual: (2.2011\times10^{-16}).
- Maximum relative projected-gradient KKT residual: (9.9943\times10^{-13}), meeting the (10^{-12}) stopping tolerance.
- TOV integrations = 0; oscillation modes = 0; synthetic or primary EOS files written = 0.

The absolute projected-gradient residual reaches (4.16\times10^{-9}) because the normalized optimization vector can have a large infinity norm. The independently checked stopping condition is the scale-relative ratio (r_{\rm KKT}/(1+\|x\|_\infty)), not the unscaled absolute value.

The downloaded artifact ZIP matches the GitHub digest, passes `unzip -t`, and all receipt/candidate/knot hashes match the internal manifest.

## Pareto/error-budget result

The following are unweighted RMS relative knot residuals against the normalized source tables, averaged over the four selected EOSs:

| Quadrature measure | \(\lambda\) | Density RMS change | Pressure RMS mismatch |
|---|---:|---:|---:|
| \(\mu_B\) | 1 | 0% | 2.967337% |
| \(\mu_B\) | \(10^{-1}\) | 0.067791% | 2.966365% |
| \(\mu_B\) | \(10^{-2}\) | 0.078534% | 2.966258% |
| \(\mu_B\) | \(10^{-3}\) | 0.080018% | 2.966241% |
| \(\mu_B\) | \(10^{-4}\) | 0.080113% | 2.966241% |
| \(\log n_B\) | 1 | 0% | 2.967337% |
| \(\log n_B\) | \(10^{-1}\) | 0.019737% | 2.967216% |
| \(\log n_B\) | \(10^{-2}\) | 0.029267% | 2.967087% |
| \(\log n_B\) | \(10^{-3}\) | 0.030711% | 2.967069% |
| \(\log n_B\) | \(10^{-4}\) | 0.030752% | 2.967068% |

At \(\lambda=10^{-4}\), the \(\mu_B\)-quadrature fit reduces the mean pressure RMS mismatch by only about **0.001095 percentage points** relative to baseline A, while permitting about 0.0801% mean density RMS movement. The \(\log n_B\)-quadrature fit reduces pressure RMS mismatch by only about **0.000268 percentage points**, with about 0.0308% density RMS movement. Log-density weighting therefore changes density less but does not improve pressure fidelity more.

The maximum pointwise pressure mismatch remains effectively unchanged. At source row 218 it is approximately:

- GPPVA_DD2: +64.0802%
- GPPVA_TM1e: +30.3615%
- SPG_M1: +15.8415%

Candidate density at row 218 changes by a negligible amount in both quadrature choices. The global RMS reduction must not be described as repair of the local sparse-cell pressure mismatch.

## Interpretation

**Established for this bounded pilot:** the joint-q solver is feasible, converges to the declared scale-relative KKT tolerance, and returns stable/causal candidate Hermite cells for all 40 cases.

**Scientific result:** this joint-q objective family stays very close to baseline A and produces no material restoration of the source pressure knots. The normalized source (q=n_B/\mu_B) is already monotone, and the pressure-only part of this constrained objective does not pull the fitted (q) far enough to change the row-218 pressure residual. This is a limited/negative result for the chosen objective and variable representation, not proof that every possible constrained reconstruction must behave the same way.

The pilot covers four tables with BBB_BHF2 core only; it is not a 234-table or core-family generalization. No new TOV, mode, or EOS artifacts were produced.

## Next gate: R18 matching/source localization

R16 found 447 raw pressure-secant decreases. Their distribution is not uniform across matching labels:

- `maxwell`: 177
- each of `css`, `cubic_F`, `cubic_P`, `quintic_F`, `quintic_P`: 54

There is also one low-density secant decrease in the (10^{-4})–(10^{-3}\,\mathrm{fm}^{-3}) bin in every one of the 234 source tables (234 total). The remaining violations occur at higher densities.

The next useful study is therefore to localize the 447 events by (i) exact source row and density, (ii) matching metadata/bridge construction, and (iii) source phase/composition indices where those records are available. The raw CompOSE package, normalized matched table and projected solver table must remain distinct provenance layers. Do not make a broad EOS correction or rerun TOV until a source-backed mechanism and a distinct reconstruction with a meaningful Pareto advantage are identified.

MD-31 ((P=0) physical low-density completion) and MD-32 (justified finite-domain perturbation boundary) remain separate. Radial (\omega_0^2) and complete-star Cowling modes remain blocked.
