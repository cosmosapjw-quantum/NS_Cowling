# R15 solver-native Hermite audit: source row 217 → 218

**Status:** `PASS_R15_NATIVE_HERMITE_ROW217_218_AUDIT`  
**Workflow:** [42654562160](https://github.com/cosmosapjw-quantum/NS_Cowling/actions/runs/42654562160)  
**Source commit:** `f4752ce5c6f0b79741dc79a672af0e39516164cf`  
**Frozen spectrum-summary SHA-256:** `ce54c111510f32086065f95356903cd6abad675adce12cab1a06b5cc82350965`

## Scope

Read-only numerical representation audit of the adjacent native source rows 217 and 218 in all 234 projected EOS tables. The audit reproduces the exact cubic-Hermite (P(\mu_b)) construction in the frozen Fortran source; it neither changes source tables nor runs TOV integrations.

Implementation sources inspected:
- `oscillation_source/methods/src/eos/ns_mu_bridge.f90`, blob `a071611c3d062eb573b130609607cb843d9775b2`;
- `oscillation_source/methods/src/eos/ns_cell_eos.f90`, blob `4bded648f775c30ba7ff65622690181a8f348314`;
- `oscillation_source/methods/src/eos/ns_table_eos.f90`, blob `8a221cb2fbee904c6f7803782ec6b06d80674e58`.

For (h=\mu_1-\mu_0) and (t=(\mu-\mu_0)/h), the cell is
[
P(t)=a_0+a_1t+a_2t^2+a_3t^3,
]
with
[
a_0=P_0,\quad a_1=hn_0,\quad
a_2=3(P_1-P_0)-h(2n_0+n_1),\quad
a_3=2(P_0-P_1)+h(n_0+n_1).
]
The local thermodynamics are
[
n_b=P'(\mu),\quad \epsilon=\mu n_b-P,\quad
\chi=\frac{c_s^2}{c^2}=\frac{n_b}{\mu P''(\mu)}.
]

## Whole-cell extrema

The production constructor performs exact interval tests:
- (P''(\mu)) is linear, so its minimum occurs at an endpoint.
- (G(\mu)=\mu P''-P') has derivative (G'=\mu P'''); since (P''') is constant, (G) is monotone or constant on the positive-(\mu) cell. Its minimum is also at an endpoint.
- Its accepted conditions are (min(P''_0,P''_1)>0) and (min(G_0,G_1)\ge0).

These imply increasing positive (n_b=P'), increasing energy density (\epsilon'=\mu P''>0), and (0<\chi\le1) across the entire interpolant.

Separately, the audit solved the polynomial stationary equation for (\chi(t)), using
[
N(t)=\mu(t)r(t)^2-hq(t)r(t)-\mu(t)q(t)r'(t),
]
where (q=a_1+2a_2t+3a_3t^2), (r=2a_2+6a_3t), and (\mu(t)=\mu_0+ht).

## Results over all 234 cells

<badge color="success">234 / 234 pass</badge>

| Diagnostic | Result |
|---|---:|
| Cells audited | 234 / 234 |
| Fortran-equivalent stable and causal cells | 234 / 234 |
| Positive (P'') over entire interval | 234 / 234 |
| Nonnegative causal margin over entire interval | 234 / 234 |
| Interior stationary points of (\chi(t)) | 0 across 234 cells |
| (\chi(t)) monotonic direction | Increasing in all 234 cells |
| Global (c_s^2/c^2) range | (1.830688\times10^{-4}) to (2.950720\times10^{-3}) |
| Equivalent (c_s/c) range | 0.0135303 to 0.0543205 |
| Minimum endpoint (P'') | (4.46817\times10^{54}\;\mathrm{Pa/J^2}) |
| Minimum causal margin (G) | (6.72345\times10^{44}\;\mathrm{m^{-3}}) |
| Maximum endpoint (|\mu n_b-P-\epsilon|/|\epsilon|) | (2.392\times10^{-16}) |
| Projected EOS source hash changes | 0 |
| Synthetic EOS rows / TOV integrations | 0 / 0 |

An independent post-artifact calculation recomputed the Hermite endpoint coefficients, (P''), and causal margins for all 234 records. It also evaluated (\chi(t)) on a 10,001-point mesh per cell; all 234 were monotonically increasing. The analytic stationary-root enumeration found no interior roots in any cell. Receipt/table/summary hashes in the artifact manifest were independently recomputed and matched.

## Relation to the sparse native gap

The native density ratio remains large, (n_{b,218}/n_{b,217}=7.43436)–11.15155; the pressure ratio is 4.96053–9.21554. The endpoint secant (Delta\ln P/\Delta\ln n_b) remains 0.789830–0.932357, but it is **not** a local adiabatic index or sound speed. The exact solver interpolant instead yields the much smaller whole-cell (c_s^2/c^2) range stated above, without an interior stationary point.

This rules out a hidden non-convex, non-monotone-density, or acausal segment in the currently admitted row-217→218 cubic-Hermite representation. It does not prove that the supplied low-density segment is a complete physical outer-crust EOS.

## Interpretation boundary and next gate

Established: numerical stability/causality of the existing projected (P(\mu_b)) interpolant on this specific cell, for all 234 tables.

Not established: physical outer-crust completeness, crust-core transition location, physical crust mass/thickness, a (P=0) surface, or a validated radial/Cowling perturbation boundary condition. The interpolant is a representation of supplied endpoint data, not new microphysical source data.

The next high-value study is a provenance and microphysics audit of the underlying low-density source segments (including CompOSE/BSk22-derived records) and the mapping/projection that produces rows 217–218. MD-31 ((P=0) physical completion) and MD-32 (finite-domain perturbation boundary) remain separate gates; radial (\omega_0^2) and complete-star Cowling modes remain blocked.

## Artifact provenance

- Artifact ID: `13689769387`
- ZIP SHA-256: `24e0bd37d8acf81f28443c5e5cd6a739b91cd3f5a0c9ca78f6e5b1fbbb6e549f`
- Receipt SHA-256: `80e214a1d8d8bda1a1e550c1cf1dd233b4da936856f1ce57397d1308f2aef98c`
- Table CSV SHA-256: `d5cc49ae219bcebcefff1143b846ecae370d7ae276805d825944de972fa859cc`
- Crust summary CSV SHA-256: `938f723067d30c14f0872672a7afd301b4958284bb581ac6353dbd6f524c298c`
