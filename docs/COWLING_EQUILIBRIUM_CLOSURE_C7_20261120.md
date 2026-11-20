# NS_Cowling C7: physical equilibrium-closure validation

Date: 2026-11-20 KST. Parent commit:
`94fb92abda7cb1c380fd8eaa1c7209748ff67da5`.
Working branch: `research/cowling-equilibrium-closure-audit-20261120`.

## Scientific question and limits

Can the source-locked sharp EOSs, on their smooth segments and Maxwell
coexistence doublets, support the *existing* equilibrium-barotropic
relativistic Cowling operator with independent numerical evidence?
They cannot, without additional fixed-composition thermodynamic derivatives,
certify a frozen-composition perturbation closure or bulk composition g-modes.

Metric (-,+,+,+). Geometric r,m in km; frequency in Hz,
`omega_geom=2*pi*f/c` with `c=299792.458 km/s`; EOS energy density and
pressure in MeV/fm³; dimensionless enthalpy
`h=log(mu/mu_surface)=int dP/(epsilon+P)`.
The existing source/model equation uses
`Gamma1=Gamma_eq=(epsilon+P)*cs_eq²/(P*c²)`; this identity is a
*model assumption*, not a validated frozen-composition material law.

## Explicit source checks, not interpolated across jumps

For each sampled smooth one-sided segment, evaluate the actual
`SharpContractEnthalpyEOS.state_at_enthalpy` at h-dh, h, h+dh.
Select `dh<=2e-6` and at most 10% of the distance to any original
thermodynamic break. Reject samples at a knot. Independently compute:

```text
dP/dh = epsilon+P
dE/dh = (epsilon+P) / (cs_eq²/c²)
cs_fd²/c² = [P(h+dh)-P(h-dh)]/[E(h+dh)-E(h-dh)].
```

Finite-difference acceptance is a scoped relative diagnostic with
predeclared tolerance 5e-5. This is not a continuum error estimate.

For every declared sharp interface, evaluate *both* original one-sided
states at identical transition enthalpy and independently require
`[P]=0`, `[mu]=0`, positive `DeltaE`,
and `DeltaE=mu_t*(n_high-n_low)`. No finite difference spans the jump.

Frozen source identity:
`data/eos_missing_data_recalc_20261119/artifacts/EOS_ALL_REQUESTED_SOURCE_LOCKED_NESTED_20261119.zip`
SHA-256:
`299de43bb4f0a87ebfc4f1139ab3c8b51c344937ec3637cdff644e9ca8ab2c04`.

For both SPG_M2/maxwell and SPG_M3/css at h =
(0.005, 0.025, 0.055, 0.09, 0.14, 0.18):

| Measured on original one-sided source | M2 Maxwell | M3 CSS |
|---|---:|---:|
| maximum relative dP/dh error | 1.807e-7 | 1.807e-7 |
| maximum relative sound-speed derivative error | 1.199e-7 | 1.199e-7 |
| maximum coexistence residual | 8.006e-16 | 8.598e-16 |
| latent heat, MeV/fm³ | 7.76618055 | 12.39583026 |

These are measured finite samples with a single source EOS cohort, not
a proof covering all 390 methods, densities, masses or compositions.
The declared Maxwell limits contain distinct finite density and energy jumps.

## Why the frozen derivative is unidentifiable

Consider a counterexample with baryon density n and composition fraction x:

```text
epsilon(n,x) = m*n + A*n²/2 + kappa*(x-x_beta(n))²/2
x_beta(n) = a*n/n_star,
P(n,x)=n*(partial epsilon/partial n)_x - epsilon.
```

Here m has units MeV, A MeV fm³, kappa MeV/fm³, n and n_star
fm⁻³, and a and x are dimensionless.
Along the beta equilibrium path `x=x_beta(n)`, for every kappa,
`epsilon_beta = m*n + A*n²/2` and `P_beta=A*n²/2`.
The equilibrium sound-speed ratio is
`cs_eq²/c²=A*n/(m+A*n)`, independent of kappa.
However, at fixed material composition
`cs_frozen²/c²=[A*n+kappa*n*(a/n_star)²]/(m+A*n)`.
Thus even an *exact* table P(epsilon) cannot determine the frozen
response when `a != 0`; different physically stable kappa values
can have identical equilibrium barotropes.
The difference is positive for kappa>0 and n>0.

In relativistic Cowling, the Schwarzschild discriminant term
`A_minus=e^(-lambda/2) [epsilon'/(epsilon+P) -
P'/(Gamma1*P)]` vanishes identically for the equilibrium
`Gamma1=Gamma_eq` and smooth barotropic backgrounds.
For frozen composition, this term generally differs from zero.
This does NOT exclude sharp-interface discontinuity g-modes under
slow conversion; those are a different restoration mechanism and
require separate physical mode identification.

## Code upgrade and scoped verification

New public source-physics check:
`compose_eos.stellar_sharp.audit_equilibrium_cowling_source`.
It emits typed PASS/FAIL with real derivative/coexistence residuals;
a synthetic 10%-corrupted cs² is rejected, and an exact interface
knotted derivative query is refused.

Both Python and native sharp Cowling `solve_mode` outputs explicitly
record `perturbation_closure=GAMMA1_EQUALS_GAMMA_EQ`,
`frozen_composition_gamma1=NOT_DERIVABLE_FROM_EQUILIBRIUM_BAROTROPE`,
and `bulk_composition_g_modes_supported=false`.
No TOV equation, Fortran derivative, Maxwell transfer, source EOS, or
shooting eigenvalue operator is changed.

Producer:
`scripts/audit_cowling_equilibrium_closure_c7_20261120.py`.
Result:
`data/cowling_equilibrium_closure_c7_20261120/RESULT.json`.
Tests:
`python_frontend/tests/test_cowling_equilibrium_closure_c7_20261120.py`.

Reproduce with `OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 PYTHONPATH=python_frontend`:

```sh
python3 scripts/audit_cowling_equilibrium_closure_c7_20261120.py
python3 -m pytest -q \
 python_frontend/tests/test_cowling_equilibrium_closure_c7_20261120.py \
 python_frontend/tests/test_sharp_cowling_hybrid_cs2_20261120.py \
 python_frontend/tests/test_sharp_cowling_hybrid_evidence_20261120.py \
 --run-slow --disable-warnings
```

Fresh observed validation (2026-11-20): the new focused tests first failed
as expected (4 missing-function failures, 1 toy control passed; later
2 native-mode output-metadata failures with 5 other tests deselected).
After repair: 7 C7 tests passed, including real M2 rapid and M3 slow
native modes. The combined source/hybrid/Maxwell/selector suite finished
with **31 passed in 117.41 seconds, exit code 0** and no skipped or
failed tests. Raw output:
`data/cowling_equilibrium_closure_c7_20261120/REGRESSION_RUN.log`.
The existing 203-test historical receipt belongs to its original version
and is not represented as freshly rerun here.
Native-subagent dispatch was unavailable (`native_subagent=NOT_RUN`);
this C7 stage used serial local execution.

## Primary theory sources

- Counsell, Gittins, Andersson & Pnigouras (2024),
  *Neutron star g-modes in the relativistic Cowling approximation*,
  https://arxiv.org/abs/2409.20178 (section 2.3 and Gamma1 closure).
- Montefusco, Antonelli & Gulminelli (2025),
  *Frozen and β-equilibrated f and p modes of cold neutron stars*,
  A&A, doi:10.1051/0004-6361/202452727.
- Tonetto & Lugones (2020), Phys. Rev. D 101, 123029,
  https://arxiv.org/abs/2003.01259 (slow/rapid discontinuity conditions).

## Next scientific step

To model frozen-composition f/p1 and bulk composition g-modes,
derive fixed-composition EOS derivatives for the selected source
microphysics, then implement the nonbarotropic Cowling terms and
repeat mode-continuation / avoided-crossing checks at physically
matched mass and EOS. Do not relabel Gamma_eq or infer a fixed-
composition derivative from an equilibrium source alone.
