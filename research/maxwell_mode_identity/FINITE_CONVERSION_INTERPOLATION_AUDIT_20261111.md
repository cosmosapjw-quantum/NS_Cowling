# Intermediate-conversion audit for sharp-interface Cowling modes

**Date:** 2026-11-11  
**Status:** `ALGEBRAIC_FAMILY_DERIVED_AND_TESTED`; `PHYSICAL_FINITE_RATE_LAW_NOT_ESTABLISHED`  
**Scope:** nonradial, even-parity ℓ=2 Cowling variables `(W,V)` on the existing analytic sharp-Maxwell fixture only. No CompOSE source is admitted and no production solver is changed.

## 1. Research question and conclusion

Can one interpolate between the already-derived slow/material and rapid/equilibrated nonradial Cowling interface conditions while preserving the weighted Abel-Liouville interface concomitant?

Yes, at the level of a **defined algebraic interpolation of kinematic constraints**. With `theta ∈ [0,1]`, the transfer matrix derived below has
`det(T_theta)=kappa=h_inner/h_outer` for every `theta`. It reduces exactly to the slow map at `theta=0`, the rapid map at `theta=1`, and the identity when the density/enthalpy jump vanishes.

This is not a finite-rate microphysical law. It contains no reaction flux, chemical affinity, mobility, conversion timescale, dissipative response, entropy-production relation, or additional interface degree of freedom. No physical interpretation of `theta` is justified by the determinant identity. A numerical scan over `theta`, should one later be authorized, must be titled **algebraic boundary-condition sensitivity scan**, not an intermediate conversion-rate calculation.

## 2. Conventions and assumptions

Use the existing static spherical Cowling convention
[
ds^2=-e^{2
u}c^2dt^2+e^{2lambda}dr^2+r^2dOmega^2,qquad (-+++).
]
At the sharp interface `r_t`, inner means the high-density side and outer the low-density side. Define every jump by
[
[X]_i^o := X_o-X_i.
]
The interface is a no-thin-shell Maxwell surface: pressure, radius and metric potentials are continuous; `h=ε+P` may jump. Let
[
A=(omega r_t/c)^2e^{lambda-2
u},quad
Q=AV+
u'W,quad
Delta P=-h e^{-lambda}Q/r_t^2,quad
kappa=h_i/h_o>0.
]
For a cold one-charge Maxwell transition, `h=μ_B n_B` and chemical-potential equality makes `kappa=n_i/n_o`. The coefficients `A` and `nu'` are evaluated on the continuous interface geometry. Here `W` has units m^3, `V` m^2, `A` and `kappa` are dimensionless, and `nu'` is m^-1.

Continuity of Lagrangian pressure gives
[
[Delta P]_i^o=0quadLongleftrightarrowquad Q_o=kappa Q_i.
]
The physical endpoint kinematic constraints are
[
	ext{slow: }[W]_i^o=0,qquad
	ext{rapid: }[W-Q/
u']_i^o=0.
]
These are the nonradial Cowling-variable forms of the sharp-interface slow and rapid conditions; they are not a finite-rate closure.

## 3. Definition of the algebraic interpolation

Define
[
(1-	heta)[W]_i^o+	heta[W-Q/
u']_i^o=0,
qquad 0leq	hetaleq1,
	ag{1}
]
while retaining `Q_o=kappa Q_i`. The parameter `theta` only mixes two already specified endpoint constraints.

Let `delta=kappa-1`. Equation (1) gives
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
Substituting this result into `Q_o=kappa Q_i` gives
[
V_o=(1-	heta)delta(
u'/A)W_i+
[1+(1-	heta)delta]V_i.
]
Thus
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

The off-diagonal units are required: `A/nu'` has units m and maps V (m^2) to W (m^3); `nu'/A` has units m^-1 and maps W to V (m^2).

## 4. Exact algebraic properties

Writing `delta=kappa-1`,
[
egin{aligned}
det T_	heta
&=(1+	hetadelta)(1+(1-	heta)delta)
-	heta(1-	heta)delta^2\
&=1+delta=kappa.
end{aligned}
]
The endpoint and smooth-limit identities are
[
T_0=egin{pmatrix}1&0\delta
u'/A&kappaend{pmatrix}=T_s,quad
T_1=egin{pmatrix}kappa&delta A/
u'\0&1end{pmatrix}=T_r,quad
T_	heta|_{kappa=1}=I.
]
For the alternating matrix `J=[[0,1],[-1,0]]`,
[
T_	heta^TJT_	heta=(det T_	heta)J=kappa J.
]
In dimension two this last identity is an algebraic consequence of the determinant, not an independent symplectic constraint.

The existing weighted concomitant for two columns of the same cold-barotropic Cowling ODE is
[
Omega=h e^{-
u}(W_1V_2-V_1W_2).
]
Its general interface ratio is
[
rac{Omega_o}{Omega_i}
=rac{h_o}{h_i}e^{-(
u_o-
u_i)}det T_	heta
=e^{-(
u_o-
u_i)}.
]
Therefore it is continuous when the background lapse is continuous, for all `theta`. If a lapse jump is present, the residual factor must be retained; no cancellation may be assumed. This is an algebraic Abel-Liouville concomitant, not a canonical covariant symplectic current.

**Independent exact checks performed:** Wolfram Language (`Module` evaluation without warnings) and SymPy both returned `det(T_theta)=kappa`, exact slow/rapid endpoints, `T^T J T-kappa J=0`, unit weighted jump factor for continuous lapse, and identity at zero jump.

## 5. Literature audit and applicability

### 5.1 Nonradial physical endpoint conditions

Tonetto & Lugones (2020), *Discontinuity gravity modes in hybrid stars: assessing the role of rapid and slow phase conversions*, Phys. Rev. D **101**, 123029, [arXiv:2003.01259](https://arxiv.org/abs/2003.01259), derive the nonradial sharp-interface conditions. Their slow case gives continuity of the radial displacement variable `W` and Lagrangian pressure; the rapid case gives continuity of Lagrangian pressure together with the appropriate jump of `W`, equivalently continuity of `V` in their formulation. They find that the discontinuity g-mode has zero frequency in the rapid limit and finite frequency in the slow case. This source supports the two endpoint conditions used above, not a finite-rate law.

### 5.2 Radial algebraic interpolation

Rau & Sedrakian (2023), *Two first-order phase transitions in hybrid compact stars: higher-order multiplet stars, reaction modes and intermediate conversion speeds*, Phys. Rev. D **107**, 103042, [arXiv:2212.09828](https://arxiv.org/abs/2212.09828), study radial oscillation/stability conditions. Their intermediate prescription interpolates radial displacement junction conditions between slow and rapid endpoints. They explain that the chosen intermediate condition does not follow from their general barotropic interface relation for any particular interface function and that additional microphysical input is required for an accurate description. It is not a derivation for the nonradial `(W,V)` Cowling system and must not be transplanted by analogy.

### 5.3 Full-GR polar algebraic interpolation

Dong, Kuan & Yunes (2026), *Hybrid-star asteroseismology across sharp phase interfaces: Slow and rapid conversion in the f- and p1-mode spectra*, [arXiv:2609.33002](https://arxiv.org/abs/2609.33002), compute full-GR polar spectra for cold, nonrotating hybrid stars and compare the physical slow and rapid limits on identical backgrounds. They also use an algebraic interpolation between the endpoint boundary conditions. Crucially, the authors explicitly state that intermediate values are not derived from a finite-rate microscopic conversion model and should not be interpreted as calibrated rates; they also state that no-thin-shell admissibility of the interpolated condition has not been established. Their full-GR spectra motivate studying interface-sensitive p1 behavior, but their numerical results are not predictions for this Cowling fixture.

### 5.4 What is and is not established

Established by primary literature for the present purpose:
- the distinct ideal slow and rapid endpoint conditions for sharp, nonradial phase interfaces;
- conversion-regime sensitivity of hybrid-star oscillations;
- examples of algebraic interpolation used for sensitivity/stability studies, with explicit warnings against interpreting it as a microphysical rate.

Derived here:
- the specific `T_theta` family for the current Cowling `(W,V)` variables under equation (1);
- determinant preservation `det(T_theta)=kappa` and consequent weighted-concomitant continuity for continuous lapse.

Not established:
- a calibrated relation between `theta` and a physical conversion timescale/rate;
- a finite-rate kinetic boundary law for nonradial Cowling oscillations;
- entropy production, interface dissipation, or the no-thin-shell consistency of the interpolated interface;
- any new numerical stellar spectrum or physical mode label.

## 6. Coding and mutation contract

The research helper is `finite_conversion_interpolation.py` and exports:
- `algebraic_interpolation_map(theta,kappa,A,nu_prime_per_m)`;
- `interpolation_constraint_residual(theta,kappa,A,nu_prime_per_m,state_inner,state_outer)`.

The returned object contains the matrix, determinant and the mandatory scope label
`ALGEBRAIC_BOUNDARY_INTERPOLATION_ONLY_NOT_FINITE_RATE_KINETICS`. It validates finite arguments, `theta in [0,1]`, and positive `kappa,A,nu_prime`. The helper is not imported by production solver modules.

Mutation rules:
- allowed: additive files under `research/maxwell_mode_identity/**`, docs, tests, receipts;
- forbidden: raw EOS, manifests, source-authority Fortran, RHS, interface map used by production, root scanner, public production API, TOV/Cowling dispatch;
- any production change requires a separate written design and explicit human approval;
- no merge is implied by publication of this evidence-only branch.

## 7. Validation evidence

TDD sequence:
1. Wrote focused tests first.
2. Observed expected RED: collection failed because `finite_conversion_interpolation` did not exist.
3. Added the smallest research-only module.
4. Fresh GREEN: `PYTHONPATH=. pytest -q research/maxwell_mode_identity/tests/test_algebraic_interpolation.py` -> **12 passed**.

Independent exact calculations:
- SymPy exact matrix identities -> PASS.
- Wolfram Language `Module` evaluation -> PASS, no messages.
- The 12 pytest cases cover both endpoints, six theta values for the determinant and weighted jump, the mixed kinematic and pressure residuals, zero jump, theta-domain rejection, and invalid/nonfinite parameters.

Not run in this increment:
- the existing PR #6 43-case suite;
- native Fortran compilation or execution;
- native/Python q-to-r RHS comparison;
- the four existing fixture-root profiles;
- two-sided matching on the actual fixture;
- zero-jump or finite-jump stellar continuation;
- any source-backed CompOSE EOS/TOV/Cowling calculation.

## 8. Research disposition

`LA9 = ALGEBRAIC_FAMILY_DERIVED_AND_TESTED__PHYSICAL_FINITE_RATE_LAW_NOT_ESTABLISHED`.

An algebraic `theta` sensitivity scan may be considered later as a separately labelled experiment, only after the local two-sided and zero-jump gates. A physical finite-rate project requires a new, source-supported interfacial constitutive model and a fresh derivation of all jump conditions. Until then:
[
oxed{	ext{PHYSICAL FINITE-RATE KINETICS: BLOCKED}}
]
