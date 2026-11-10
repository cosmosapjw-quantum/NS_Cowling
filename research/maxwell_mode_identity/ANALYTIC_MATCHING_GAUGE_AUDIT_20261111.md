# LA11: finite-interval analytic characteristic determinant, basis/gauge covariance, and interface chart singularity

Date: 2026-11-11  
Program: NS_COWLING_MODE_IDENTITY  
Parent: research/maxwell-interface-interpolation-la10-v1 at 104b4fa76e191db42e6af26f072533617f45441c  
Evidence: conditional derivation; exact SymPy/Wolfram checks; 22 synthetic tests. Native equations and stellar endpoint frames remain unverified.

## 1. Finite-interval characteristic determinant

This is a finite-interval boundary-value problem, not automatically a whole-line Evans problem. No far-field stable/unstable splitting, essential spectrum, exponential dichotomy, or asymptotic spectral gap is assumed.

On each smooth layer let
\[
Y'=M_j(r,\zeta)Y,\qquad \zeta\in\Omega\subset\mathbb C,
\]
and at each sharp interface let
\[
y(r_k^+,\zeta)=T_k(\zeta)y(r_k^-,\zeta).
\]
The center-admissible and surface-admissible subspaces have constant dimensions \(p\) and \(n-p\), with full-rank frames \(\Phi_c,\Phi_s\) at a common matching radius. Define
\[
D(\zeta;r_m)=\det[\Phi_c(r_m,\zeta)\;\Phi_s(r_m,\zeta)].
\]
The determinant is bilinear, without complex conjugation.

Sufficient conditions for holomorphy:
1. Each \(M_j(r,\zeta)\) is holomorphic in \(\zeta\), belongs to \(L^1\) in \(r\), and has locally uniformly bounded \(L^1\) norm on compact subsets of \(\Omega\).
2. Each graph interface map \(T_k(\zeta)\) is holomorphic and invertible on the claimed domain, or the complete implicit interface relation is used in a square, full-rank matching matrix.
3. Endpoint dimensions are fixed and full-rank holomorphic frames exist.
4. The parameter domain contains no unhandled coefficient, endpoint, or interface singularity.

The layer propagator is the solution of a Volterra integral equation; its ordered-integral series converges locally uniformly under the \(L^1\) bound, hence is holomorphic in \(\zeta\). Finite products with holomorphic interface maps stay holomorphic. Therefore \(D\) is holomorphic under the stated conditions.

The zero-set statement is
\[
D(\zeta;r_m)=0
\Longleftrightarrow
\operatorname{Ran}\Phi_c(r_m,\zeta)\cap
\operatorname{Ran}\Phi_s(r_m,\zeta)\ne\{0\}.
\]
This is the finite-interval eigenvalue condition.

**Multiplicity gate:** holomorphy and zero-set equivalence do not alone prove that the order of a zero equals algebraic multiplicity of the physical operator. That requires the actual BVP to be a regular holomorphic Fredholm pencil of index zero with nonempty resolvent, and the matching matrix to be an analytic invertible reduction of that pencil. Until proved for the actual Cowling BVP, multiplicity is UNVERIFIED.

## 2. Basis covariance

For holomorphic invertible endpoint basis changes \(A_c(\zeta)\) and \(A_s(\zeta)\),
\[
\widetilde\Phi_c=\Phi_cA_c,\qquad
\widetilde\Phi_s=\Phi_sA_s,
\]
so
\[
\boxed{\widetilde D=D\,\det A_c\,\det A_s.}
\]
The multiplier must be a holomorphic unit throughout the domain, including the contour interior. Being nonzero only on the contour is insufficient: a basis multiplier may have an interior zero and add a spurious zero/winding contribution.

## 3. Gauge covariance and trace cocycle

For a piecewise-\(C^1\), invertible state transformation \(y=G(r,\zeta)z\), holomorphic in \(\zeta\),
\[
M_z=G^{-1}M_yG-G^{-1}G',\qquad
T_{k,z}=G_o^{-1}T_{k,y}G_i,
\]
and at a common matching radius
\[
D_y=\det G(r_m,\zeta)D_z.
\]
Liouville's formula gives
\[
\operatorname{tr}M_z=\operatorname{tr}M_y-\partial_r\log\det G.
\]
For the LA9 weight \(w=he^{-\nu}\), the trace-first drift \(\Gamma=\operatorname{tr}M+h'/h-\nu'\) transforms by the same derivative. The interface multiplier transforms as \(J_z=J_y\det G_i/\det G_o\). When all smooth and one-sided interface factors are retained, the global compensated determinant changes only by the reference-frame factor \(\det G(r_0)\). This is algebraic covariance, not a new physical conservation law.

## 4. Matching-radius factors

For both frame blocks propagated through the same smooth layer,
\[
D(\zeta;r_2)=\det U(r_2,r_1;\zeta)D(\zeta;r_1)
=\exp\left(\int_{r_1}^{r_2}\operatorname{tr}M\,dr\right)D(\zeta;r_1).
\]
Crossing an interface consistently gives
\[
D(\zeta;r_k^+)=\det T_k(\zeta)D(\zeta;r_k^-).
\]
The root set and zero orders are independent of match radius when these factors are nonzero; the determinant value and phase are not.

## 5. LA10 transfer pole and implicit interface matching

Write
\[
A(\omega)=\left(\frac{\omega r_t}{c}\right)^2e^{\lambda_{\rm met,t}-2\nu_t},
\quad \kappa=\frac{h_i}{h_o},\quad
Q=AV+\nu'W,\quad [X]_i^o=X_o-X_i.
\]
The LA10 constraints are
\[
Q_o-\kappa Q_i=0,\qquad
(1-\theta)[W]_i^o+\theta[W-Q/\nu']_i^o=0.
\]
For \(0\le\theta<1\), multiplying the second condition by \(\nu'\) gives
\[
(1-\theta)\nu'[W]_i^o-\theta A[V]_i^o=0.
\]
Thus \(C_i y_i+C_o y_o=0\), with
\[
C_i=\begin{pmatrix}-\kappa\nu'&-\kappa A\\-(1-\theta)\nu'&\theta A\end{pmatrix},
\quad
C_o=\begin{pmatrix}\nu'&A\\(1-\theta)\nu'&-\theta A\end{pmatrix}.
\]
Exact determinants are
\[
\det C_i=-\kappa A\nu',\qquad \det C_o=-A\nu'.
\]
For \(A\ne0\), \(T_\theta=-C_o^{-1}C_i\), which reproduces the LA10 transfer and \(\det T_\theta=\kappa\). But
\[
(T_\theta)_{21}=(1-\theta)(\kappa-1)\frac{\nu'}{A}\sim\omega^{-2}.
\]
For \(0\le\theta<1,\kappa\ne1\), this is a graph-chart singularity at \(A=0\), not a divergence of the implicit blocks. At \(A=0\), \(C_o\) loses rank but the combined block matrix \([C_i\;C_o]\) retains row rank two. There is no unique outer state determined by the inner state, so retain an implicit/multiple-shooting interface relation.

For one-dimensional inner and outer endpoint lines spanned by \(u_i,v_o\), define
\[
E_{\rm imp}=\det[C_i u_i,C_o v_o].
\]
For \(A\ne0,\theta<1,\kappa\ne1\),
\[
\boxed{E_{\rm imp}=A\nu'\det[T_\theta u_i,v_o].}
\]
This removes the \(1/A\) transfer pole and remains finite as \(A\to0\) when the endpoint frames and blocks are analytic. This is an algebraic matching extension, not evidence for a physical zero-frequency stellar mode.

At \(\theta=1\), the rapid constraint is equivalent to \([V]=0\) for \(A\ne0\). The original divided equation loses rank at \(A=0\); use the explicitly selected regular extension \([V]=0\), with
\[
C_i^{(r)}=\begin{pmatrix}-\kappa\nu'&-\kappa A\\0&-1\end{pmatrix},
\quad C_o^{(r)}=\begin{pmatrix}\nu'&A\\0&1\end{pmatrix},
\quad T_r=\begin{pmatrix}\kappa&(\kappa-1)A/\nu'\\0&1\end{pmatrix}.
\]
For \(\kappa=1\), there is no enthalpy jump; return the exact identity map before forming rank-degenerate blocks. These are algebraic branch conventions, not finite-rate kinetics.

## 6. Conditioning and analytic normalization

Use dimensionless states \(\widehat W=W/L^3,\widehat V=V/L^2\), and \(x=A/(\nu' L)\):
\[
\widehat T_\theta=
\begin{pmatrix}
1+\theta\delta&\theta\delta x\\
(1-\theta)\delta/x&1+(1-\theta)\delta
\end{pmatrix},\quad \delta=\kappa-1.
\]
For \(\theta<1,\delta\ne0\), the lower-left entry scales as \(x^{-1}\), while the determinant stays \(\kappa\). Hence generically \(\operatorname{cond}_2(\widehat T)=O(|x|^{-2})=O(|\omega|^{-4})\). This is meaningful only in the stated dimensionless state coordinates.

Euclidean norms, Gram-Schmidt, ordinary unitary QR/SVD and polar factors involve conjugation and are generally nonholomorphic in complex \(\omega\). They can be real-frequency/pointwise conditioning diagnostics, but not the analytic characteristic function used for the argument principle. Any QR, exterior-product, continuous-orthogonalization, or Grassmann-chart method must restore its full phase, volume/radial factor, triangular determinant and chart-transition cocycles. For the current two-state system direct bilinear determinants are simplest.

## 7. Numerical evidence

The research-only helper analytic_matching.py and test_analytic_matching.py pass 22 synthetic tests. Coverage includes bilinear matching determinants, endpoint basis changes, gauge trace/interface cocycles, smooth match-radius Liouville factors, LA10 implicit blocks and graph transfer agreement, the \(A\nu'\) scalar prefactor, zero-frequency graph-pole handling, rapid regular extension, zero-jump identity, nonholomorphic norm normalization, interior-zero basis multiplier, dimensionless conditioning and the finite \(A\to0\) implicit scalar limit.

Final run: 22 passed in 0.17 s; collection: 22 in 0.01 s. Environment: Python 3.13.5, NumPy 2.3.5, SciPy 1.17.0, pytest 9.0.2. SymPy 1.14.0 returned exact zero for seven symbolic checks; a clean Wolfram Language Module returned exact zero for ten checks. The condition scan shows cond2 quadrupling when x is halved while det remains 2.8. The one synthetic implicit determinant has a finite A-to-zero limit.

All inputs are synthetic. Native q-to-r, actual endpoint frames, stellar roots, actual BVP-pencil equivalence and production behavior were not tested.

## 8. Literature

- Beyn, Latushkin, Rottmann-Mathes, arXiv:1210.3952: Theorem 3.9 gives a holomorphic determinant for complementary holomorphic subspace frames and zeros exactly at nontrivial intersections. The Fredholm-pencil/partial-multiplicity results motivate the separate multiplicity gate; they do not remove the need to prove analytic equivalence for the actual Cowling BVP. https://arxiv.org/abs/1210.3952
- Humpherys & Zumbrun, arXiv:math/0508020: polar Evans computation retains an analytic scalar volume factor with a continuous-orthogonalization angular equation; exterior-power cost grows combinatorially. https://arxiv.org/abs/math/0508020
- Ledoux, Malham & Thümmler, arXiv:0710.1037: Grassmannian spectral shooting uses coordinate-patch switching and preserves spectral-parameter analyticity under its construction. https://arxiv.org/abs/0710.1037
- Zumbrun, arXiv:0904.0268v2: centered coordinates and continuous-orthogonalization error analysis in asymptotically constant whole-line problems; its gap hypotheses are not imported into this finite-interval problem. https://arxiv.org/abs/0904.0268v2
- Malham & Niesen, arXiv:math/0605581: Magnus/Gauss-Legendre Evans computation in a stiff Schrödinger example; method-specific order-reduction analysis. https://arxiv.org/abs/math/0605581

## 9. Non-claims

No whole-line Evans theorem, essential-spectrum result, canonical covariant symplectic current, physical finite-rate interface law, physical zero-frequency stellar mode, native Fortran/q-to-r validation, actual EOS/TOV/Cowling result, or physical f/p1/interface-g label is claimed. Algebraic multiplicity and argument-principle root counting remain blocked until their analytic pencil and contour conditions are verified for the actual solver.
