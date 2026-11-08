# Native Fortran master-function reconstruction v1

## 1. Scope

This note documents the native reconstruction operation exposed as \`reconstruct_master\`. It rebuilds a cold, neutrino-transparent, beta-equilibrated, single-phase barotropic EOS from the supplied free-energy-per-baryon master function. The Fortran routine, rather than the Python reference implementation, owns the interpolation, composition update, outer-crust pressure correction and physical admissibility checks.

The operation is deliberately narrower than a complete port of the historical generator. Its smoothing is a local weighted-cubic regularisation followed by PCHIP, not SciPy's global \`UnivariateSpline\`; it does not port the old global candidate-selection or anomaly-repair policy in full. The old Python generator is a reference for physical conventions and differential context, not an authority when its numerical interpolation conflicts with thermodynamic closure.

The source table may carry composition or may be a composition-free barotrope. In the latter case the routine does not invent proton, nuclear, electron or muon fractions. In the former case it applies the limited baryon/charge/lepton reconstruction described below. This routine is not a microscopic nuclear-structure solver and does not calculate the sequence of equilibrium nuclei in the outer crust.

## 2. CompOSE variables and units

The implementation follows the CompOSE thermodynamic table convention (Typel et al. 2022). The baryon density \(n_B\) is in fm\(^{-3}\), while \(P\) and \(\epsilon\) are in MeV fm\(^{-3}\). The \`mn\` field stores the neutron rest energy \(m_n c^2\) in MeV. The table variables are

\[
Q_1=\frac{P}{n_B},\quad
Q_2=\frac{s}{n_B},\quad
Q_3=\frac{\mu_B}{m_n c^2}-1,
\]
\[
Q_4=\frac{\mu_q}{m_n c^2},\quad
Q_5=\frac{\mu_l}{m_n c^2},\quad
Q_6=\frac{f}{n_Bm_n c^2}-1,\quad
Q_7=\frac{\epsilon}{n_Bm_n c^2}-1.
\]

Thus \(Q_1\) has units of MeV, while \(Q_2,\ldots,Q_7\) are dimensionless. We retain \(\hbar c\) explicitly in the lepton and lattice formulae; no undocumented conversion between fm\(^{-1}\) and MeV is made.

For a cold \(T=0\) barotrope, the free-energy and energy-per-baryon variables should agree. The native v1 input contract requires \(T=0\) metadata, a beta-equilibrium flag, negligible entropy \(Q_2\), and a neutrino-transparent lepton-potential column \(Q_5\simeq0\). A finite \(Q_5\) can represent a different lepton/neutrino thermodynamic state; it must not be silently overwritten by a reconstruction that sets \(Q_5=0\).

## 3. Thermodynamic master-function representation

We define the master function

\[
F(x)=\frac{f}{n_B},\qquad x=\ln n_B,
\]

where \(F\) is in MeV. For the supplied CompOSE table,

\[
F_i=(1+Q_{6,i})m_n c^2.
\]

The first-law-consistent reconstruction derives every thermodynamic quantity from this same \(F(x)\):

\[
Q_1(x)=\frac{dF}{dx},\qquad
\mu_H(x)\equiv\frac{d\epsilon}{dn_B}=F(x)+Q_1(x),
\]
\[
P(x)=n_B Q_1(x),\qquad
\epsilon(x)=n_B F(x),\qquad
Q_6(x)=Q_7(x)=\frac{F(x)}{m_n c^2}-1.
\]

Since \(x=\ln n_B\),

\[
\frac{dF}{dx}=n_B\frac{dF}{dn_B}
\]

and therefore

\[
\mu_H=F+n_B\frac{dF}{dn_B}
=\frac{d\epsilon}{dn_B}.
\]

The pressure follows from the zero-temperature Euler relation,

\[
P=n_B\mu_H-\epsilon=n_B Q_1.
\]

Differentiating the pressure with respect to density gives

\[
D\equiv\frac{dP}{dn_B}=Q_1+\frac{dQ_1}{d\ln n_B}
=F'(x)+F''(x).
\]

Consequently, the adiabatic sound-speed ratio is

\[
\boxed{\frac{c_s^2}{c^2}=\frac{dP}{d\epsilon}
=\frac{Q_1+F''}{F+Q_1}=\frac{D}{\mu_H}.}
\]

This representation is the governing definition of the native output. The routine does not interpolate \(P\), \(\epsilon\), and the enthalpy derivative independently. CompOSE Q3 is a separate convention-dependent column. For \(I_l=0\), the repository verifier uses
\[
(1+Q_3)m_nc^2=\mu_H-Y_qQ_4m_nc^2,
\]
where \(Y_q\) is the charge-fraction metadata value. Output Q3 therefore subtracts \(Y_qQ_4\) from \((F+Q_1)/(m_nc^2)-1\) for \(I_l=0\). For the supported \(I_l=1\) branch, \(Q_5\simeq0\) and the verifier's relation reduces to \((1+Q_3)m_nc^2=\mu_H\). In both branches the sound-speed denominator is \(\mu_H=d\epsilon/dn_B\), not the convention-adjusted Q3 quantity.

## 4. Interpolation and regularisation

The interpolation variable is \(x=\ln n_B\). Two paths are provided.

### 4.1 Unsmoothed path

With \`smooth=False\`, the input samples \(F_i\) are interpolated by the native monotonicity-preserving PCHIP routine. The PCHIP first derivative supplies \(Q_1=F'\), and its analytic second derivative supplies \(F''\). The interpolation is local and shape preserving; it avoids unconstrained high-order spline oscillations between monotone source knots.

### 4.2 Smoothed path

With \`smooth=True\`, the source \(F_i\) are regularised by a local weighted cubic least-squares fit centred at every source knot. The local coordinate is

\[
z_j=\frac{x_j-x_i}{h_i},
\]

where \(h_i\) is the largest distance from the centre to the selected local window. The weights are

\[
w_j=\exp(-2z_j^2).
\]

A polynomial \(a_0+a_1z+a_2z^2+a_3z^3\) is fitted by the \(4\times4\) weighted normal equations, and the smoothed value at the centre is \(a_0\). The half-window radius is controlled by \`smooth_sigma\` and the number of source points. The default is 0.05; the accepted interval is \(0<\texttt{smooth\_sigma}\le0.5\).

The smoothed samples are then passed to PCHIP, and \(Q_1\) and \(F''\) are calculated analytically from that interpolant. This is a defined native regularisation scheme. It is **not** claimed to reproduce the global-knot, smoothing-factor or boundary behaviour of SciPy's historical \`UnivariateSpline\`. Users should assess \`smooth=True\` and \`smooth=False\` as two distinct regularisation choices and test convergence as the source density grid is refined.

The local normal equations are solved with partial pivoting. A singular local system or nonfinite result is a hard error; there is no silent fallback to an unrelated Python numerical routine.

## 5. Mechanical stability and causality

For a cold single-component barotrope, the required conditions are

\[
F>0,\quad Q_1=F'>0,\quad
\mu_H=F+Q_1>0,\quad
D=Q_1+F''>0,
\]
\[
\mu_H-D\ge0
\quad\Longleftrightarrow\quad
0<\frac{c_s^2}{c^2}\le1.
\]

On a PCHIP interval, write \(t=(x-x_i)/h\in[0,1]\). The cubic Hermite master function has coefficients

\[
F(t)=A+Bt+Ct^2+Et^3,
\]
\[
A=F_i,\quad B=hm_i,\quad
C=3(F_{i+1}-F_i)-h(2m_i+m_{i+1}),
\]
\[
E=2(F_i-F_{i+1})+h(m_i+m_{i+1}),
\]

where \(m_i,m_{i+1}\) are the PCHIP endpoint slopes with respect to \(x\).

It follows that \(F\) and \(\mu_H-D\) are cubic polynomials in \(t\), whereas \(Q_1\) and \(D\) are quadratic. The native gate therefore finds the endpoint and interior stationary-point values of the relevant polynomials and tests their minima over the complete interval. For example, the minimum of a cubic on \([0,1]\) is obtained by checking both endpoints and all real roots in \((0,1)\) of its quadratic derivative. A fixed sample grid is not used as the proof for the PCHIP core intervals.

The outer-crust correction integrates a separate PCHIP interpolant for \(Q_1\). Since the integral of each cubic PCHIP segment is quartic, the native gate constructs quartic and cubic coefficients and checks both endpoints and all internal stationary points of \(F\), \(Q_1\), \(\mu_H\), \(D\), and \(\mu_H-D\). The derivative cubic is partitioned at its quadratic critical points; tangential roots are retained and sign-changing monotone intervals are bisected. This replaces the former fixed 256-subinterval mesh. It is a floating-point root-isolation check, not a formal interval-arithmetic proof. The report separately records whether the outer-crust correction was applied.

Any failure of positivity, mechanical stability or causality rejects the reconstructed EOS. The code reports the offending interval and the local thermodynamic quantities; it does not clip \(c_s^2\), flatten the pressure, or return a partially admissible table.

## 6. Composition reconstruction

Composition processing is conditional on the presence and provenance of \`eos.compo\`. If the source has no composition, the routine retains its \(Q_4\) values and its composition placeholders; it does not fabricate \(Y_p\), \(Y_e\), \(Y_\mu\), or nuclear mass/charge numbers.

For a composition-bearing table, the source fractions are clipped to their physical intervals and the baryon fractions are reconstructed as

\[
Y_n=\mathrm{clip}(1-Y_p-Y_h,0,1),\qquad
Y_h=\mathrm{clip}(1-Y_p-Y_n,0,1).
\]

For a nuclear group with mean mass and charge numbers \(A\) and \(Z\), define

\[
Y_q=Y_p+\frac{Z}{A}Y_h,
\]

with \(Z/A=0\) where no nuclear group exists. Charge neutrality requires

\[
Y_e+Y_\mu=Y_q.
\]

At \(T=0\), the Fermi momentum and chemical potential of a lepton species \(\ell\) are

\[
k_{F\ell}=(3\pi^2 n_B Y_\ell)^{1/3},
\qquad
\mu_\ell=\sqrt{(m_\ell c^2)^2+(\hbar c k_{F\ell})^2}.
\]

If the electron-only solution has \(\mu_e\le m_\mu c^2\), the native routine sets \(Y_\mu=0\) and \(Y_e=Y_q\). Otherwise, it solves the monotone scalar equation

\[
\mu_e(n_B,Y_q-Y_\mu)-\mu_\mu(n_B,Y_\mu)=0
\]

on \(0\le Y_\mu\le Y_q\) by bisection. It then sets \(Y_e=Y_q-Y_\mu\) and \(Q_4=-\mu_e/(m_n c^2)\). The output sets \(Q_2=Q_5=0\) and \(Q_6=Q_7\).

This solves the electron-muon Fermi-gas partition at fixed supplied proton/nuclear charge fraction. It does not solve the nuclear many-body equations that determine \(Y_p\), \(A\), or \(Z\); those are inherited from the supplied composition source. Accordingly, the operation is not a replacement for an RMF, Skyrme, APR or nuclear statistical-equilibrium model.

## 7. Nuclear dissolution repair

The optional \`fix_dissolution\` operation addresses a specific legacy-data pathology at the last nuclear row before homogeneous matter. It evaluates the charge-closure residual

\[
r_{\rm ch}=
\left|Y_p+\frac{Z}{A}Y_h-(Y_e+Y_\mu)\right|.
\]

It also checks whether the reconstructed electron fraction at the last nuclear row has an anomalous drop to less than one half of the smaller neighbouring homogeneous/crust value. If either the residual exceeds \`dissolution_tol\` (default 0.015) or the drop criterion is triggered, the routine interpolates the total charge fraction in \(\ln n_B\) between the preceding crust row and the next homogeneous row. At the anomalous row it sets the nuclear group variables \(A,Z,Y_h\) and the muon fraction to zero, and assigns the interpolated charge fraction to \(Y_p\). It then recomputes the electron/muon partition.

This is a data-repair heuristic for a narrow class of remnant-nucleus inconsistencies. It is not a dynamical dissolution calculation, a minimisation of nuclear free energy, or a universal crust-core transition prescription. Set \`fix_dissolution=False\` when the source composition is already trusted or when the heuristic is not justified by provenance.

## 8. Outer-crust electron and lattice pressure correction

The outer-crust branch is considered only when composition is available and the low-density rows satisfy the source pattern: a contiguous region with \(A>0\), \(Y_n\simeq0\), and negligible muon abundance. The routine identifies the last such row as the neutron-drip-side endpoint and chooses a join point at least eight source rows above it, or 2.5 per cent of the table size, whichever is larger, limited by the table end.

For the relativistic, zero-temperature electron Fermi gas, let

\[
n_e=n_B Y_e,\qquad
k_{Fe}=(3\pi^2n_e)^{1/3},\qquad
x_e=\frac{\hbar c\,k_{Fe}}{m_ec^2}.
\]

Its pressure is

\[
P_e=\frac{(m_ec^2)^4}{24\pi^2(\hbar c)^3}
\left[
x_e(2x_e^2-3)\sqrt{1+x_e^2}
+3\,\mathrm{asinh}(x_e)
\right].
\]

For small \(x_e\), direct evaluation suffers cancellation. The implementation uses the series

\[
x_e^5\left(\frac85-\frac47x_e^2+\frac13x_e^4+O(x_e^6)\right)
\]

for the bracket when \(x_e<0.05\). The series coefficients were checked symbolically from the closed form.

The repository's outer-crust correction uses the body-centred-cubic Wigner-Seitz Coulomb-lattice approximation

\[
\epsilon_{\rm lat}
=-0.9\left(\frac{4}{9\pi}\right)^{1/3}
\alpha Z^{2/3}\hbar c\,k_{Fe}n_e,
\qquad
P_{\rm lat}=\frac13\epsilon_{\rm lat}.
\]

The imposed outer-crust master derivative is then

\[
Q_1^{\rm outer}=\frac{P_e+P_{\rm lat}}{n_B}.
\]

The code joins these low-density values to the master derivative at the selected high-density join through a shape-preserving PCHIP in \(\ln n_B\). It reconstructs the free energy backwards from the join by integrating \(Q_1\,d\ln n_B\), so the First Law remains closed by construction. The operation does **not** calculate the equilibrium nuclide sequence, nuclear masses, electron-capture thresholds, neutron drip, or the full BPS outer-crust EOS. It uses the supplied \(A,Z,Y_e\) composition and applies the electron-plus-lattice pressure correction only.

The native report returns \`outer_crust_correction_applied\`; this is false when the input does not meet the detection conditions. A request to reconstruct an EOS does not imply that an outer-crust correction was applied.

## 9. Additional thermodynamic columns

The routine preserves \`q_extra\` by default. Provider-specific extra thermo columns do not have a universal meaning, so the code does not silently overwrite them. If the caller explicitly sets \`recompute_extra_cs2=True\`, the input must have exactly one extra column (\`Nadd=1\)); the routine then records

\[
q_{\rm extra}=\frac{c_s^2}{c^2}.
\]

This flag is a semantic assertion by the caller that the single extra column represents sound speed squared. \`Nadd=0\` is accepted when this option is false; other counts are rejected.

## 10. Strict v1 input contract

The current native v1 path requires:
1. At least five finite source density points, with positive and strictly increasing \(n_B\).
2. Finite positive neutron rest energy \(m_nc^2\), supported \(I_l\in\{0,1\}\), and finite charge-fraction metadata.
3. \(|T|\le10^{-12}\) MeV, beta-equilibrium metadata, \(\max|Q_2|\le2\times10^{-8}\), and \(\max|Q_5|\le10^{-8}\).
4. The cold-source identity \(\max|Q_6-Q_7|\le2\times10^{-6}\).
5. A rowwise source Q3/CompOSE Euler identity relative residual no larger than \(2\times10^{-6}\). For \(I_l=0\), this identity includes \(-Y_qQ_4m_nc^2\); for \(I_l=1\), the correction is absent under the Q5 gate.
6. A single-valued smooth master curve. An equal-pressure pair across a change in \`Iphase\` is rejected because a single \(F(\ln n_B)\) interpolant must not cross a first-order density discontinuity.
7. \`Nadd=0\` or \(1\). Recomputing \(c_s^2\) in the extra column requires \`Nadd=1\`.
8. If composition is marked available, the source \`eos.compo\` provenance and all required composition arrays must be present and finite.

The native code validates complete PCHIP polynomial intervals for both the core branch and integrated outer-crust bridge. A rejected input is not silently passed to the legacy Python implementation.

## 11. Python API and native CLI

Python:

    from compose_eos import reconstruct_master

    rebuilt, report = reconstruct_master(
        eos,
        smooth=True,
        smooth_sigma=0.05,
        fix_dissolution=True,
        dissolution_tol=0.015,
        recompute_extra_cs2=False,
        return_report=True,
    )

The API writes a temporary internal table, invokes the Fortran executable, reads the native output, and then applies the independent C1-C6/C10-C11 consistency gate. It returns either \`EOSData\` or \`(EOSData, report)\`. The report records the interpolation choice, First-Law closure, composition provenance, whether composition was reconstructed, and whether the outer-crust correction was actually applied.

Native CLI, after preparing internal CSV/metadata directories:

    compose_eos_fortran reconstruct-master \
        INPUT_INTERNAL OUTPUT_INTERNAL \
        true 0.05 true 0.015 false

The final flag declares whether the single extra thermo column should be recomputed as \(c_s^2/c^2\).

## 12. Regression evidence and explicit limitations

The native test set contains analytic and adversarial fixtures for:
- exact first-law closure from a linear \(F(\ln n_B)\);
- preservation of composition-free \(Q_4\) rather than fabricated fractions;
- unsmoothed PCHIP reconstruction;
- baryon and charge conservation and \(\mu_e=\mu_\mu\) above the muon threshold;
- dissolution of an inconsistent last nuclear row;
- relativistic-electron plus lattice pressure in an outer-crust fixture;
- preservation of an extra column unless its semantics are explicitly declared as \(c_s^2\);
- rejection of a pressure plateau across an \`Iphase\` boundary;
- rejection of invalid smoothing, nonzero entropy and nonzero neutrino/lepton potential;
- rejection of Q6/Q7 and source Q3/Euler mismatches;
- preservation of the distinct Il=0 Q3 charge-potential convention.

The production GPPVA_FSU2H normalized table passes its native CompOSE and C1-C11 intake checks but is **rejected by this v1 master reconstruction** at the interpolation physical-admissibility gate. This is an important negative result: passing the intake consistency checks does not prove that the selected smoothed \(F(\ln n_B)\) PCHIP is mechanically stable and causal everywhere. The Fortran routine reports the interval and the local thermodynamic quantities instead of silently producing a table. The unsmoothed source path remains a separate option; model-specific follow-up should determine whether the source interpolation, local smoothing width, or outer-crust join causes the first failure before any alternative regularisation is accepted.

This operation does not yet port:
- the historical global UnivariateSpline smoothing factor and candidate-selection policy;
- universal nuclear dissolution or crust-core minimisation;
- the equilibrium nuclear species sequence and complete BPS outer-crust generation;
- microscopic RMF, Skyrme, APR or quark-matter solvers;
- Gibbs mixed phases or phase-interface TOV/Love-number integration;
- full radial-mode or stable-branch TOV sequence diagnostics.

## 13. References

1. S. Typel et al. (2022), *CompOSE Reference Manual v3.01*, [arXiv:2203.03209](https://arxiv.org/abs/2203.03209). The source for \(Q_1,\ldots,Q_7\), the separation of thermodynamic and composition tables, and the table metadata.
2. G. Baym, C. Pethick and P. Sutherland (1971), *The Ground State of Matter at High Densities: Equation of State and Stellar Models*, The Astrophysical Journal **170**, 299-317, [doi:10.1086/151216](https://doi.org/10.1086/151216). A foundational calculation of the cold outer-crust equilibrium nuclide sequence with Coulomb-lattice contributions.
3. G. Baym, H. A. Bethe and C. J. Pethick (1971), *Neutron Star Matter*, Nuclear Physics A **175**, 225-271, [doi:10.1016/0375-9474(71)90281-8](https://doi.org/10.1016/0375-9474(71)90281-8). Microscopic/thermodynamic context for dense neutron-star matter and the crust-core regime.
4. L. Suleiman, J. L. Zdunik and P. Haensel (2024), *Layers of electron captures in the crust of accreting neutron stars*, Astronomy & Astrophysics **690**, A301, [doi:10.1051/0004-6361/202348340](https://doi.org/10.1051/0004-6361/202348340). This work gives the zero-temperature degenerate-electron EOS and a linear-mixing-rule lattice pressure for parent/daughter nuclear mixtures. The present implementation uses only the single-species BCC limit with the supplied nuclear charge number; it does not implement the mixture function or electron-capture kinetics.

The thermodynamic identities and interval-extremum criteria in this note are derived from the stated master function and implemented in the native Fortran module. The local smoothing and dissolution rules are implementation choices whose numerical convergence and model dependence should be assessed with each source EOS.
