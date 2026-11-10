# Intermediate-conversion audit for sharp-interface Cowling modes

**Date:** 2026-11-11  
**Disposition:** `ALGEBRAIC_FAMILY_DERIVED_AND_TESTED`; `PHYSICAL_FINITE_RATE_LAW_NOT_ESTABLISHED`.  
**Scope:** research-only, nonradial sharp-Maxwell Cowling interface in the existing analytic two-phase fixture. No CompOSE EOS, production solver, TOV sequence, or actual stellar root is used here.

## 1. Research question and result

Can the slow/material and rapid/equilibrated nonradial Cowling interface conditions be connected by a one-parameter *algebraic* family while retaining the weighted Abel–Liouville interface concomitant?

Yes. Under the explicit interpolation defined below, the transfer matrix has `det(T_theta)=kappa=h_inner/h_outer` for all `theta in [0,1]`. The endpoints reproduce the slow and rapid maps, and the zero-jump limit is the identity.

This is **not** a finite-rate kinetic law. The parameter `theta` is not a reaction rate, conversion timescale, mobility, or calibrated microscopic coefficient. The construction supplies no reaction flux, chemical affinity, interface dissipation, entropy production, or extra interface degree of freedom. Its admissibility as a no-thin-shell junction condition has not been established. Any scan over `theta` must be described as an *algebraic boundary-condition sensitivity scan*, never as a physical conversion-rate scan.

## 2. Conventions and assumptions

Use signature `(-,+,+,+)` and a static spherical background
[
ds^2=-e^{2
u}c^2dt^2+e^{2lambda}dr^2+r^2dOmega^2.
]
At the interface, inner (`i`) is the high-density side and outer (`o`) the low-density side. Define
[
[X]_i^o:=X_o-X_i.
]
Let `h=epsilon+P` and
[
A=(omega r_t/c)^2e^{lambda-2
u},qquad
Q=AV+
u'W,qquad
Delta P=-h e^{-lambda}Q/r_t^2,qquad
kappa=h_i/h_o>0.
]
The interface geometry and the coefficients `A,nu_prime` are common and continuous on the two sides. Pressure continuity gives
[
[Delta P]_i^o=0quadLongleftrightarrowquad Q_o=kappa Q_i.
]
For a cold one-conserved-charge Maxwell transition, chemical-potential equality gives `h=mu_B*n_B`, hence `kappa=n_i/n_o`. This thermodynamic interpretation is conditional on that one-charge cold equilibrium.

The existing nonradial Cowling endpoint constraints in the current `(W,V)` convention are
[
	ext{slow: }[W]_i^o=0,qquad
	ext{rapid: }[W-Q/
u']_i^o=0.
]
Here `W` has units m^3, `V` m^2, `A` and `kappa` are dimensionless, and `nu_prime` has units m^-1.

## 3. Definition and derivation

Define the *algebraic kinematic interpolation*
[
(1-	heta)[W]_i^o+	heta[W-Q/
u']_i^o=0,qquad 0leq	hetaleq1,
	ag{1}
]
while retaining pressure continuity `Q_o=kappa Q_i`. Set `delta=kappa-1`. Equation (1) implies
[
W_o-W_i=	heta(Q_o-Q_i)/
u'
=	hetadelta Q_i/
u'.
]
Since `Q_i=A V_i+
u'W_i`,
[
W_o=[1+	hetadelta]W_i+	hetadelta(A/
u')V_i.
]
Using pressure continuity to solve for `V_o` then gives
[
V_o=(1-	heta)delta(
u'/A)W_i+
[1+(1-	heta)delta]V_i.
]
Therefore
[
oxed{
T_	heta=
egin{pmatrix}
1+	heta(kappa-1)&	heta(kappa-1)A/
u'\
(1-	heta)(kappa-1)
u'/A&1+(1-	heta)(kappa-1)
end{pmatrix}},qquad
inom{W_o}{V_o}=T_	hetainom{W_i}{V_i}.
	ag{2}
]
The off-diagonal units are required: `A/nu_prime` maps V (m^2) to W (m^3), and `nu_prime/A` maps W to V.

## 4. Exact algebraic properties

With `delta=kappa-1`,
[
egin{aligned}
det T_	heta
&=(1+	hetadelta)(1+(1-	heta)delta)
-	heta(1-	heta)delta^2\
&=1+delta=kappa.
end{aligned}
]
The limits are
[
T_0=egin{pmatrix}1&0\delta
u'/A&kappaend{pmatrix}=T_{m slow},
quad
T_1=egin{pmatrix}kappa&delta A/
u'\0&1end{pmatrix}=T_{m rapid},
quad T_	heta|_{kappa=1}=I.
]
For `J=[[0,1],[-1,0]]`,
[
T_	heta^TJT_	heta=(det T_	heta)J=kappa J.
]
In two dimensions this follows identically from the determinant; it is not an independent symplectic constraint.

For two solutions of the same smooth cold-barotropic Cowling system, the previously derived algebraic concomitant is
[
Omega=h e^{-
u}(W_1V_2-V_1W_2).
]
Its interface jump is
[
rac{Omega_o}{Omega_i}
=rac{h_o}{h_i}e^{-(
u_o-
u_i)}det T_	heta
=e^{-(
u_o-
u_i)}.
]
It is continuous when lapse is continuous. If lapse is discontinuous, retain the residual exponential factor; do not assume it cancels. This is a bilinear Abel–Liouville concomitant, not a claimed canonical covariant symplectic current.

**Independent exact verification:** Wolfram Language `Module` evaluation and SymPy exact algebra both returned zero residuals for the endpoint identities and `T^TJT-kappa J`, `det(T_theta)=kappa`, unit weighted jump for continuous lapse, and identity at zero jump.

## 5. Literature boundary

- **Nonradial physical endpoints:** Tonetto & Lugones (2020), *Discontinuity gravity modes in hybrid stars: assessing the role of rapid and slow phase conversions*, Phys. Rev. D 101, 123029, [arXiv:2003.01259](https://arxiv.org/abs/2003.01259), derive slow and rapid junction conditions for nonradial sharp-interface perturbations. The paper supports the endpoint conditions, not a finite-rate law.
- **Radial algebraic interpolation:** Rau & Sedrakian (2023), *Two first-order phase transitions in hybrid compact stars: higher-order multiplet stars, reaction modes and intermediate conversion speeds*, Phys. Rev. D 107, 103042, [arXiv:2212.09828](https://arxiv.org/abs/2212.09828), studies radial mode stability and an algebraic interpolation. It is not a derivation for the nonradial Cowling `(W,V)` system and must not be transplanted by analogy.
- **Full-GR polar algebraic interpolation:** Dong, Kuan & Yunes (2026), *Hybrid-star asteroseismology across sharp phase interfaces: Slow and rapid conversion in the f- and p1-mode spectra*, [arXiv:2609.33002](https://arxiv.org/abs/2609.33002), also uses an algebraic interpolation. The authors explicitly state intermediate values are not derived from a finite-rate microscopic model or calibrated rates; they also state that no-thin-shell admissibility has not been established.

Thus literature supports the sharp-interface physical endpoints and documents algebraic interpolations in radial/full-GR contexts, but it does not establish a calibrated finite-rate model for this Cowling construction.

## 6. Code and mutation contract

The helper `finite_conversion_interpolation.py` is research-only and exports:
- `algebraic_interpolation_map(theta,kappa,A,nu_prime_per_m)`
- `interpolation_constraint_residual(theta,kappa,A,nu_prime_per_m,state_inner,state_outer)`

The returned map records `model_scope=ALGEBRAIC_BOUNDARY_INTERPOLATION_ONLY_NOT_FINITE_RATE_KINETICS`. Inputs must be finite, `theta in [0,1]`, and `kappa,A,nu_prime_per_m > 0`. This helper is not imported by production modules.

Allowed mutations: additive research files under `research/maxwell_mode_identity/**`, documentation, tests and receipts. Forbidden: source-authority Fortran, raw EOS, C7X admission/manifests, production RHS/interface/root algorithms, TOV/Cowling dispatch, or physical rate inference. Any production/API mutation requires a separate written design and explicit human approval.

## 7. Validation

- TDD RED: focused tests failed collection because the helper module was absent.
- GREEN: `PYTHONPATH=/mnt/data/NS_Cowling_research_20261111/finite_rate_audit python -m pytest -q /mnt/data/NS_Cowling_research_20261111/finite_rate_audit/research/maxwell_mode_identity/tests/test_algebraic_interpolation.py` -> **12 passed**.
- Fresh rerun: **12 passed in 0.04 s**.
- `python -m compileall -q research/maxwell_mode_identity`: exit 0.
- AST parsing of all Python files and JSON parsing of the generated status: PASS.
- SymPy exact identities: PASS.
- Wolfram Language exact `Module` check: PASS without evaluator warnings.

The 12 tests cover both endpoints, six values of `theta` for determinant and weighted jump, pressure and mixed-kinematic residuals, zero-jump identity, invalid domain and nonfinite/nonpositive parameter rejection.

## 8. Explicit non-results

Not run in this increment: existing PR #6 full-suite replay; PR #7 compensated-Abel suite replay; native Fortran build/run; native q-to-r RHS equivalence; actual four fixture-root profile processing; two-sided matching on the actual stellar fixture; zero-jump or finite-jump stellar continuation; source-backed CompOSE/TOV/Cowling; physical finite-rate kinetics; physical f/p1/discontinuity-g classification. GitHub Actions status for the new commit is not claimed unless a nonempty status result is fetched.

`LA10=ALGEBRAIC_FAMILY_DERIVED_AND_TESTED__PHYSICAL_FINITE_RATE_LAW_NOT_ESTABLISHED`.
