# v5 paired crust-core matching-method sensitivity

Date: 2026-11-09

## Design

The v5 crossed cohort contains 13 crust candidates, 3 admitted core families and 6 historical smooth matching prescriptions:

\[
13\times3\times6=234
\]

The admissible core families are BBB_BHF2, GM_GM1 and RG_SLY4 (78 rows each). The six methods are CSS, cubic_F, cubic_P, smooth near-pressure-plateau maxwell, quintic_F and quintic_P (39 rows each).

For each fixed core/crust pair, define the six-method mean \(\bar y_6\), then centre each method response:

\[
\delta_m^{(y)}(c,s)=100\left(\frac{y_m(c,s)}{\bar y_6(c,s)}-1\right),
\qquad y\in\{f,p_1\}.
\]

This paired centring does not privilege one interpolation method as the reference. The method contrasts sum to zero within every pair by construction.

## Results at the v5 target mass

The spectrum cohort uses a target mass of approximately 1.33 solar masses, with the relativistic Cowling approximation and the declared finite positive surface-pressure cutoff.

| Diagnostic | f-mode | p1-mode |
|---|---:|---:|
| Median paired six-method range | 0.0274% | 0.1727% |
| 90th percentile paired range | 0.0733% | 0.4599% |
| Maximum paired range | See machine-readable sensitivity result | 2.2958% |

The largest p1 response occurs for BBB_BHF2 + GPPVA_FSU2. This is a tail case, not a representative prescription effect. In particular, the p1 response cannot be described as uniformly sub-percent across the full crossed cohort.

The maxwell token in this table is the historical **smooth near-pressure-plateau bridge**. It has no density discontinuity or latent heat and is not the true first-order Maxwell construction described in MAXWELL_REFERENCE_VALIDATION_V1_20261109.md.

## Factor structure

For the finite balanced cohort, the descriptive crossed-design variance decomposition gives:

- Core-associated share: 99.91% for f and 94.20% for p1.
- Crust-associated share: approximately 5.63% for p1.
- The matching-method main-effect share is approximately 0.0158% for p1.

Paired conditional spans provide a complementary scale comparison. The median spans across the selected pairs are 20.975% (core), 0.569% (crust) and 0.006% (method) for radius; 26.112%, 0.942% and 0.027% for f; and 16.707%, 6.825% and 0.173% for p1.

These are descriptive decompositions of this selected finite cohort, not posterior variance components, population estimates or causal effects. Two candidate cores were excluded by the source-level causality gate; the 234 combinations are a crossed design, not 234 independent microscopic EOS models.

## Interpretation boundary

At the target mass, the typical f-mode response to the smooth matching prescriptions is very small, while p1 is more sensitive and has a non-negligible tail. Core-EOS variation dominates both mode frequencies in this cohort. No mass-dependent conclusion follows because this release fixes the target mass near 1.33 solar masses.

The report and paired-sensitivity calculation are analysis/reporting updates. They do not modify the pinned v5 solver source or the canonical 234-row spectrum.
