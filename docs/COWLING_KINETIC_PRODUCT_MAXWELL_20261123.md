# Cowling kinetic product and ideal Maxwell interfaces

Status: derived and literature-supported; selected algebra checked symbolically.
This note does not change C17/C19 branch acceptance, and no new norm diagnostic
is implemented in C20. It complements the native TOV interpolation migration.

## Conventions and source mapping

Baseline code: C19 122c502d1ef5e55fd890f47d86444df725987032.
Relevant definitions: oscillation_source/notes/equations/cowling.md,
oscillation_source/notes/equations/interface.md,
fortran_eos/src/compose_eos_cowling_dvode.f90, and profile assembly in
python_frontend/compose_eos/stellar_sharp_native.py.

Use signature (-,+,+,+),
\(ds^2=-e^{2\Phi}c^2dt^2+e^{2\Lambda}dr^2+r^2d\Omega^2\),
\(\xi^r=e^{-\Lambda}W Y_{\ell m}/r^2\) and
\(\xi^a=-V D^aY_{\ell m}/r^2\), with unit-normalized harmonics.
Physical amplitudes have [W]=length^3 and [V]=length^2.
Native variables are \(w=W/r^{\ell+1},v=V/r^\ell\).

The exported profile divides both physical W and V by the same peak |W|.
Consequently exported W is dimensionless and V_km_inverse has units km^-1.
An absolute energy inferred from these normalized profiles needs a separate
physical amplitude convention. Normalized overlaps are unaffected by a common
amplitude factor.

## Same-background kinetic bilinear

For one static, spherical, perfect-fluid background in the relativistic Cowling
approximation, the positive kinetic bilinear is

\[
K(\eta,\xi)=\int_{\Sigma_t}\frac{\epsilon+P}{c^2}e^{-\Phi}
\gamma_{ij}\eta^{i*}\xi^j\,dV_\gamma .
\]

Substituting the repository displacement convention and integrating the angles,

\[
K_{ab}=\sum_p\int_p \frac{\epsilon+P}{c^2}e^{\Lambda-\Phi}
\left(\frac{W_a^*W_b}{r^2}+\ell(\ell+1)V_a^*V_b\right)\,dr.
\]

The general invariant expression is supported independently by Counsell,
Gittins, Andersson and Pnigouras, "Neutron star g-modes in the relativistic
Cowling approximation", arXiv:2409.20178v1, section2.4 equation72;
MNRAS536,1967–1979(2025), doi:10.1093/mnras/stae2721.
https://arxiv.org/html/2409.20178v1

Their radial W convention differs; their following specialized expression
must not be copied into this code without a variable transformation.
Here the angular integral is
\(W_a^*W_b/r^4+\ell(\ell+1)V_a^*V_b/r^2\) and
\(dV_\gamma=e^\Lambda r^2dr\,d\Omega\), giving the expression above.

For H=epsilon+P>0, K_aa>0 for a nonzero displacement. With SI physical
amplitudes, K has units kg m^2. The instantaneous kinetic energy referred
to asymptotic time is \(T_\infty=K(\dot\xi,\dot\xi)/2\);
a complex harmonic amplitude has cycle average
\(\langle T_\infty\rangle=\omega^2K_{aa}/4\).

For implementation in regular native variables,

\[
K_{ab}=\sum_p\int_p \frac{H}{c^2}e^{\Lambda-\Phi}r^{2\ell}
(w_a^*w_b+\ell(\ell+1)v_a^*v_b)\,dr.
\]

The integrand is O(r^{2ell}) near a regular center. For ell=2 it is O(r^4).
A flat, constant-H analytic check with
W=A r^{ell+1}, V=-A r^ell/ell, ell=2 gives
\(K=H A^2(R^5-r_0^5)/(2c^2)\).

## Ideal slow and rapid interface terms

Assume conservative adiabatic bulk fluid, a material interface without
a massive shell or surface inertia, and ideal real instantaneous slow/rapid
conversion laws. The antisymmetric Green boundary expression is, up to a
continuous common metric/area factor,

\[
B(\eta,\xi)=\eta^{r*}\Delta P_\xi-\xi^r\Delta P_\eta^* .
\]

Replacing Eulerian deltaP by Lagrangian DeltaP does not change B because
the two P' eta* xi terms cancel.

Slow conversion has [xi^r]=[DeltaP]=0, hence [B]=0.
Ideal rapid conversion has [DeltaP]=0 and [xi^r-DeltaP/P']=0, hence

\[
[B]=[1/P'](\Delta P_\eta^*\Delta P_\xi-
                    \Delta P_\xi\Delta P_\eta^*)=0.
\]

These are the junction conditions of Tonetto and Lugones,
Phys.Rev.D101,123029(2020), arXiv:2003.01259, sectionIII,
slow equations32–33, rapid equations34,41 and alternative42:
https://arxiv.org/abs/2003.01259
https://doi.org/10.1103/PhysRevD.101.123029

Using k=omega/c, P'=-H Phi', and this code's traction convention,

\[
\Delta P=-H e^{-\Lambda}
 (k^2r^2e^{\Lambda-2\Phi}V+\Phi'W)/r^2,
\qquad
\xi^r-\Delta P/P'=-k^2e^{-2\Phi}V/\Phi'.
\]

For nonzero frequency and continuous r,Phi,Phi', the rapid condition is
therefore [V]=0 (equivalently [v]=0). Frequency dependence in the W,V transfer
matrix alone does not imply an additional interface kinetic mass term.
Finite-rate conversion, dissipative response or an interface shell changes
the assumptions and requires a new boundary calculation.

This cancellation supports the conservative same-background operator domain,
with a regular center and free surface DeltaP=0. It is not a proof of
the code's discrete self-adjointness, completeness, or exact orthogonality
at a finite center/surface cutoff. Those remain separate numerical questions.

## Cross-background diagnostic contract

For distinct stars A and B, a phasewise coordinate transport u in [0,1]
must be chosen. One useful construction is

\[
q_{A,p}(u)=\sqrt{(H_A/c^2)e^{\Lambda_A-\Phi_A}\,dr_A/du}\,
\begin{pmatrix}W_A/r_A\\\sqrt{\ell(\ell+1)}V_A\end{pmatrix}.
\]

Then the phasewise L2 self-product equals K_AA, while the cross-product
defines a coordinate-transport-dependent diagnostic. Its normalized magnitude
obeys Cauchy–Schwarz. It is not an orthogonality theorem between different
stellar operators.

A future native diagnostic should report same-background f/p1 products,
per-phase kinetic fractions, phasewise mapped products and vanishing-phase
norm statuses. Integrate phases separately; preserve physical Maxwell jumps
and handle C19's within-phase shooting duplicate explicitly. Do not replace
the original high-phase0.95 gate with a globally dominated weighted overlap.
The M3 factor1.05 p1 rejection remains an original result.

## Symbolic and execution status

A native child derivation inspected source equations and primary literature.
Wolfram evaluated the rapid Green cancellation and the rapid displacement
identity above; a separate exact ell=2 flat-background integral residual was0.
The first child CAS invocation emitted an H-symbol warning; the final integral
check was warning-free. The child's WolframContext search failed internally.
Parent independently checked the general kinetic-product source expression.
The mathematical diagnostic is not implemented or numerically evaluated in C20.
