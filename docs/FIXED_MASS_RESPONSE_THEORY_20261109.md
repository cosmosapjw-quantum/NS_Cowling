# Fixed-mass response of Cowling modes under crust-core EOS changes

**Status:** analytic identity derived; Wolfram symbolic check returned the identity; deterministic mock experiment numerically checked. This is not a new Cowling perturbation equation or a sharp-interface jump law.

## Definitions and assumptions

Let λ parameterize a differentiable family of cold barotropic EOS constructions. Use the dimensionless central coordinate
\[
x=\ln(P_c/P_\star),\qquad P_\star>0
\]
and define TOV mass and a node-identified Cowling eigenfrequency by \(M=M(x,\lambda)\), \(\omega=\omega(x,\lambda)\). The comparison is at fixed gravitational mass \(M_0\), not fixed central pressure.

Assume locally \(M,\omega\in C^1\), the chosen star lies on one specified connected branch, \(M_x\ne0\), and the Cowling root is simple and mode identity can be continued without crossing/exchange.

## Theorem: fixed-mass response

Under these assumptions there is a local differentiable \(x(\lambda)\) satisfying \(M(x(\lambda),\lambda)=M_0\), and
\[
\boxed{\frac{d\omega}{d\lambda}\bigg|_{M_0}
=\omega_\lambda-\omega_x\frac{M_\lambda}{M_x}.}
\tag{1}
\]
All partial derivatives are evaluated at fixed value of the other argument.

**Proof.** Let \(F(x,\lambda)=M(x,\lambda)-M_0\). Since \(F=0\) and \(F_x=M_x\ne0\), the implicit-function theorem supplies the local function \(x(\lambda)\). Differentiate:
\[
M_x\frac{dx}{d\lambda}+M_\lambda=0,\quad
\frac{dx}{d\lambda}\bigg|_{M_0}=-\frac{M_\lambda}{M_x}.
\tag{2}
\]
Applying the chain rule to \(\omega(x(\lambda),\lambda)\) gives (1). \(\square\)

The first term is the direct EOS response at fixed Pc. The second is the TOV reparameterization needed to keep M fixed. This term is generally nonzero when an EOS change alters M at fixed Pc.

### Root-residual form

If the simple mode root is defined by a two-sided matching residual \(D(\omega,x,\lambda)=0\), with \(D_\omega\ne0\), then \(\omega_\lambda=-D_\lambda/D_\omega\) and \(\omega_x=-D_x/D_\omega\). Therefore
\[
\boxed{
\frac{d\omega}{d\lambda}\bigg|_{M_0}
=-\frac{D_\lambda-D_xM_\lambda/M_x}{D_\omega}.
}
\tag{3}
\]
This can be implemented around the native Cowling shooting determinant. Multiplying D by a smooth nonzero normalization leaves the ratio invariant at a simple root. The result also applies to \(f=\omega/(2\pi)\).

## Fold normal form near the first mass maximum

At a nondegenerate first maximum \(x_t(\lambda)\), assume \(M_x=0\), \(M_{xx}=-A(\lambda)<0\). Locally,
\[
M(x,\lambda)=M_{\max}(\lambda)
-\frac{A(\lambda)}{2}(x-x_t(\lambda))^2+O((x-x_t)^3).
\tag{4}
\]
For target mass \(M_0<M_{\max}\), set \(\delta=M_{\max}-M_0>0\). The two roots are
\[
x_\pm=x_t\pm\sqrt{2\delta/A}.
\tag{5}
\]
The initial rising-mass branch is the lower-x root. At fixed M0,
\[
\frac{dx_-}{d\lambda}
=x_t'-\frac{M_{\max}'}{\sqrt{2A\delta}}
+\frac{A'}{2A}\sqrt{\frac{2\delta}{A}}+O(1).
\tag{6}
\]
If \(\omega_x(x_t)\ne0\) and \(M_{\max}'\ne0\), the generic fixed-mass response scales as \(\delta^{-1/2}\). This is not unconditional: if \(\omega_x=O(x-x_t)=O(\sqrt{\delta})\), the leading divergence can cancel. More generally, \(\omega_x=O(\delta^\alpha)\) gives a structural term \(O(\delta^{\alpha-1/2})\).

At the fold \(M_x=0\), so the IFT hypothesis fails. For \(M_0>M_{\max}\), no star on this branch exists; below the maximum, the same M may correspond to two Pc values on the two sides of the fold. A branch label is therefore part of the scientific datum.

**Practical consequence:** report \(M_x\), distance to the first turning point and mass coverage with every fixed-mass response table. Do not interpolate through the maximum or label a truncated sequence endpoint as Mmax. Use fold-aware continuation if the turning point itself is studied.

## Thermodynamically admissible continuous path in q=n/μ

The six historical matching prescriptions are discrete. A derivative with respect to a method index is undefined until a continuous EOS path is declared. For a cold barotrope parameterized by chemical potential μ>0, define \(q(\mu)=n(\mu)/\mu\). Then
\[
dP/d\mu=\mu q,\qquad \epsilon=\mu^2q-P,\qquad
c_s^2=\frac{dP}{d\epsilon}=\frac{q}{q+\mu q'}.
\tag{7}
\]
If q>0 and q′≥0, then 0<cs²≤1.

**Lemma (convex path preservation).** Let qA and qB be C1 on [μ1,μ2], positive and nondecreasing. Suppose they have identical endpoint values and equal weighted pressure integral:
\[
q_A(\mu_i)=q_B(\mu_i),\quad
\int_{\mu_1}^{\mu_2}\mu q_A\,d\mu=
\int_{\mu_1}^{\mu_2}\mu q_B\,d\mu.
\]
Then \(q_\lambda=(1-\lambda)q_A+\lambda q_B\) preserves positivity, monotonicity, endpoint n, endpoint P and endpoint ε when all paths share P(μ1). Every intermediate path has 0<cs²≤1.

**Proof.** Positivity and nonnegative derivative are preserved by convex combinations. Endpoint q is linear, and the pressure increment is the weighted integral of q, also linear in q. Thus endpoint n and P are λ-independent. The cold Euler identity then fixes endpoint ε, and (7) supplies the sound-speed bound. \(\square\)

This does not imply that the six historic bridges already lie on one common path satisfying these constraints. The actual EOS tables must be audited and mapped to such a path first. Otherwise λ describes a synthetic continuation only.

## Numerical estimators and acceptance

For a finite six-method set at fixed pair and mass,
\[
{\cal R}_6^{(y)}=100\frac{\max_m y_m-\min_m y_m}{\bar y_6},\quad
\Delta_m^{(y)}=100(y_m/\bar y_6-1),\quad y\in\{f,p_1\}.
\]
Use complete six-method cells only; no imputation and no interpretation of 234 rows as 234 independent microscopic EOS samples.

For a continuous-path derivative, estimate \(M_\lambda,M_x,\omega_\lambda,\omega_x\) using refined central differences or a verified tangent/adjoint calculation. Independently retune Pc at λ±h to the same M and compare the finite difference of re-solved frequencies with (1). The fixed-Pc difference is a distinct quantity.

Acceptance: same node-identified mode branch; both stars on the same connected rising branch; exact target-mass closure; derivative stable under h-halving; converged mode roots/residuals; reported |Mx|/fold warning; and no derivative claim for discrete method labels without a declared continuous path.

## Scope

This is a sensitivity identity and thermodynamic-path lemma. It does not derive the Cowling equations, determine frozen-composition adiabatic response, establish radial stability, or derive sharp-interface slow/rapid conversion laws.


## First-variation system for the relativistic Cowling equations

For the current smooth ell-pole system, put y=(W,V)^T, chi=dP/d epsilon, L=ell(ell+1), and a=(omega*r/c)^2 exp(Lambda-2 nu). Then y'=A y with A=[[nu'/chi, a/chi-L exp(Lambda)],[-exp(Lambda)/r^2,2 nu']]. For any differentiable parameter zeta, the tangent S_zeta=partial_zeta y at fixed omega obeys S_zeta'=A S_zeta+A_zeta y, with initial data obtained by differentiating the regular centre expansion.

The full A_zeta must include nu'_zeta, Lambda_zeta, nu_zeta and chi_zeta. At fixed background and varying only chi, A_zeta has first row (-nu'*chi_zeta/chi^2, -a*chi_zeta/chi^2) and a zero second row. For the frequency tangent on a fixed background, A_omega has only the upper-right entry 2a/(omega*chi).

If D(omega,zeta)=0 is a simple shooting/surface residual, then d omega/d zeta=-D_zeta/D_omega. At fixed gravitational mass this combines with the fixed-mass identity above to give d omega/d zeta at fixed M = -(D_zeta-D_x*M_zeta/M_x)/D_omega. Derivatives of the centre expansion, surface/event location and interface/transfer conditions are part of the physical calculation; the synthetic coefficient-only test does not provide them. See docs/COWLING_FIRST_VARIATION_RESPONSE_20261109.md for the full derivation and acceptance contract.
