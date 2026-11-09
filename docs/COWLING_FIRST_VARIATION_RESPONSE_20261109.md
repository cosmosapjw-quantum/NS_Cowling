# First-variation response system for relativistic Cowling modes

Date: 2026-11-09  
Track: smooth barotropic crust-core matching cohort  
Status: analytic derivation and deterministic synthetic ODE check completed in the analysis environment; exact-head GitHub CI pending

## 1. Scope and claim boundary

This note differentiates the existing smooth, relativistic Cowling fluid system with respect to a model parameter. It is a forward-sensitivity construction for the current two-variable shooting equations, not a new Cowling approximation, a new variational principle, or a physical prediction. The synthetic test below has a prescribed coefficient profile, not a TOV solution or CompOSE EOS.

The distinction matters because relativistic Cowling pulsation equations and nearly Sturm-Liouville formulations are established (Finn 1988). A variational principle for nonradial relativistic stellar pulsations was already derived by Detweiler and Ipser (1973). The present contribution is narrower: an implementation-ready tangent system for the existing W,V equations and its coupling to fixed-gravitational-mass differentiation.

## 2. First-order matrix form

Use the background metric
\[
ds^2=-e^{2\nu(r)}c^2dt^2+e^{2\Lambda(r)}dr^2+r^2d\Omega^2
\]
and \(\chi(r)=dP/d\epsilon=c_s^2\). With harmonic convention \(e^{i\omega t}\), define
\[
a(r,\omega)=\left(\frac{\omega r}{c}\right)^2e^{\Lambda-2\nu},
\qquad L_\ell=\ell(\ell+1),
\qquad Q=\nu'W+aV.
\]
The current smooth-fluid equations can be written as
\[
y'=\mathsf A y,\qquad y=\begin{pmatrix}W\\V\end{pmatrix},
\]
where
\[
\mathsf A=
\begin{pmatrix}
\nu'/\chi & a/\chi-L_\ell e^\Lambda\\
-e^\Lambda/r^2 & 2\nu'
\end{pmatrix}.
\tag{1}
\]
For \(\ell=2\), \(L_\ell=6\). This is an algebraic rearrangement of the current shooting system; it does not alter the EOS closure or boundary conditions.

## 3. Forward tangent equation

Let \(\zeta\) be a differentiable parameter of a smooth EOS/background family and
\[
S_\zeta(r)=\partial_\zeta y(r,\omega,\zeta)\big|_\omega
\]
the state sensitivity at fixed trial frequency. Differentiating (1) gives
\[
S_\zeta'=\mathsf A S_\zeta+\mathsf A_\zeta y.
\tag{2}
\]
The tangent initial data are the derivative of the regular centre expansion, not an arbitrary zero vector. For a fixed regular normalization independent of \(\zeta\), the leading centre data have zero tangent; higher-order coefficient corrections must be differentiated consistently.

At fixed \(r,\omega\), with \(c\) held fixed,
\[
a_\zeta=a(\Lambda_\zeta-2\nu_\zeta)
\]
and
\[
\mathsf A_\zeta=
\begin{pmatrix}
\nu'_\zeta/\chi-\nu'\chi_\zeta/\chi^2 &
a_\zeta/\chi-a\chi_\zeta/\chi^2-L_\ell e^\Lambda\Lambda_\zeta\\
-e^\Lambda\Lambda_\zeta/r^2 & 2\nu'_\zeta
\end{pmatrix}.
\tag{3}
\]
For the frequency tangent, the background is fixed and
\[
a_\omega=\frac{2a}{\omega},\qquad
\mathsf A_\omega=
\begin{pmatrix}
0 & 2a/(\omega\chi)\\
0 & 0
\end{pmatrix}.
\tag{4}
\]
For an EOS-path parameter that changes only the adiabatic coefficient \(\chi\) on a fixed background,
\[
\mathsf A_\zeta=
\begin{pmatrix}
-\nu'\chi_\zeta/\chi^2 & -a\chi_\zeta/\chi^2\\
0&0
\end{pmatrix}.
\tag{5}
\]
Equation (5) is a useful local diagnostic, but a physical EOS variation generally also changes the TOV background \((\nu,\Lambda)\), the surface radius, and possibly the location of EOS events. Those contributions must be retained before interpreting the tangent as a physical sensitivity kernel.

## 4. Root derivative and fixed-mass composition

For a surface residual \(D(\omega,\zeta)=Q(R;\omega,\zeta)\), with \(R\) fixed in the differentiation,
\[
D_\zeta=\nu'_\zeta W+\nu' S_{\zeta,W}
+a_\zeta V+aS_{\zeta,V},
\]
and the analogous \(D_\omega\) contains \(a_\omega V+aS_{\omega,V}+\nu'S_{\omega,W}\). If the surface/event location changes, include the endpoint term \(D_rR_\zeta\) and differentiate the adopted finite-pressure surface condition itself.

For a two-sided shooting determinant at a fixed matching radius,
\[
D=\det(y_-,y_+),\qquad
D_\zeta=\det(S_{\zeta,-},y_+)+\det(y_-,S_{\zeta,+}).
\tag{6}
\]
If the matching radius moves, include its derivative before evaluating the root; for a smooth common coefficient interval the determinant's matching-radius derivative is proportional to \(D\) and therefore vanishes at a simple root, but this cancellation must not be carried across a coefficient jump or a material interface without checking the transfer law.

At a simple root \(D=0\), \(D_\omega\ne0\),
\[
\frac{d\omega}{d\zeta}\bigg|_{\rm branch}=-\frac{D_\zeta}{D_\omega}.
\tag{7}
\]
Combining this with the fixed-mass identity from FIXED_MASS_RESPONSE_THEORY_20261109.md, with \(x=\ln(P_c/P_\star)\), yields
\[
\boxed{
\left.\frac{d\omega}{d\zeta}\right|_{M_0}
=-\frac{D_\zeta-D_xM_\zeta/M_x}{D_\omega}
}
\tag{8}
\]
provided the same connected TOV branch is used and \(M_x\ne0\). The formula is branch-local and becomes ill-conditioned near a turning point.

## 5. Deterministic synthetic ODE experiment

The regression mock integrates (1) with \(\ell=2\), \(c=1\), radius \(R=1\), prescribed \(\nu=r^2/2\), \(\Lambda=0\), and
\[
\chi(r,\zeta)=0.3+\zeta r^2(1-r)^2.
\]
The centre data are \(W(r_0)=r_0^3\), \(V(r_0)=-r_0^2/2\), and the surface residual is \(D=Q(R)\). Both forward tangent systems (2) for \(\zeta\) and \(\omega\) are integrated with DOP853. A root is refined independently at \(\zeta=0,\pm h\) by Brent's method.

Analysis-environment result:
- dimensionless trial root: \(\omega=8.146085425792169\);
- absolute surface residual: \(8.1\times10^{-14}\);
- tangent prediction: \(d\omega/d\zeta=0.5254905905549\);
- central-difference relative errors at \(h=10^{-3},5\times10^{-4},2.5\times10^{-4}\): \(6.1\times10^{-9},1.5\times10^{-9},3.8\times10^{-10}\).

The approximately factor-four reduction under each halving of \(h\) is the expected second-order central-difference trend. The root is deliberately **not** assigned an \(f\)- or \(p_1\)-family identity. All numbers are dimensionless and synthetic.

Reproduction is provided by research_artifacts/mass_response/mock_cowling_variation.py; the acceptance test is in python_frontend/tests/test_fixed_mass_response.py. The exact-head CI run is the authority for implementation verification.

## 6. Required physical extension before a kernel claim

1. Differentiate the actual TOV background along a declared thermodynamically admissible EOS path to obtain \(\nu_\zeta,\Lambda_\zeta,\nu'_\zeta\), and \(\chi_\zeta\).
2. Differentiate the finite-positive-pressure surface condition and any pressure-event location; the present v5 spectrum is not a \(P=0\) vacuum-surface calculation.
3. Check the derivative of centre regularity coefficients, EOS-cell transitions, and the two-sided matching determinant.
4. Compare tangent predictions with independently retuned \(\zeta\pm h\) roots while preserving mode identity and the same rising-mass branch.
5. For a true sharp Maxwell interface, differentiate the slow/rapid transfer maps and include their boundary terms. The smooth W,V node identity and the smooth tangent equation cannot simply be transplanted to that branch.

No scalar K_cs2(r) kernel, full-GR correction, frozen-composition closure, radial-stability result, or physical mode prediction is claimed here.

## References

- L. S. Finn, “Relativistic stellar pulsations in the Cowling approximation,” MNRAS 232 (1988) 259–275, doi:10.1093/mnras/232.2.259.
- S. L. Detweiler and J. R. Ipser, “A variational principle and a stability criterion for the nonradial modes of pulsation of stellar models in general relativity,” ApJ 185 (1973) 685–707, doi:10.1086/152447.
