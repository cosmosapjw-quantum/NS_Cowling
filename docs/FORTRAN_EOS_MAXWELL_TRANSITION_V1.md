# Native Fortran Maxwell transition v1

## 1. Scope and interpretation

We consider two supplied, smooth, cold, beta-equilibrated, one-dimensional barotropic phases. We construct the **Maxwell limit with one globally conserved charge, baryon number, and local electric neutrality imposed separately in each phase**. It is not a Gibbs mixed phase, not a smooth crust-core interpolation, and not the historical maxwell option in the six-method smooth-matching API.

The inputs are two already-computed phase branches. This module does not derive either branch from a microscopic nuclear or quark model. It does not calculate nuclear composition, electron-capture chains, muon thresholds, Coulomb lattice energies, neutron drip, surface tension, pasta geometry, or screening. It does not infer local electric neutrality from a one-dimensional table. The caller must explicitly assert that both supplied phases represent locally charge-neutral matter.

The native kernel lives in fortran_eos/src/compose_eos_phase_transition.f90. The command-line action is maxwell-transition; the Python frontend exposes construct_maxwell_transition. Python only serializes the arrays, invokes the executable, and returns the native report. The thermodynamic interpolation, root isolation, stability checks, causality checks, and construction of the discontinuity are implemented in Fortran.

## 2. CompOSE convention and units

The definitions below follow the updated CompOSE Reference Manual (Typel et al. 2022). The CompOSE thermodynamic table defines

\[
Q_1=\frac{P}{n_B},\qquad
Q_2=\frac{s}{n_B},\qquad
Q_3=\frac{\mu_B}{m_n}-1,
\]
\[
Q_4=\frac{\mu_q}{m_n},\qquad
Q_5=\frac{\mu_l}{m_n},\qquad
Q_6=\frac{f}{n_Bm_n}-1,\qquad
Q_7=\frac{\epsilon}{n_Bm_n}-1.
\]

Here \(n_B\) is in fm\(^{-3}\), \(P\) and \(\epsilon\) are in MeV fm\(^{-3}\), \(m_n\), \(m_p\), \(\mu_B\), \(\mu_q\), and \(\mu_l\) are in MeV, and \(Q_2\) is dimensionless. Thus

\[
P=n_BQ_1,\qquad
\epsilon=n_Bm_n(1+Q_7),\qquad
\mu_B^{(Q_3)}=m_n(1+Q_3).
\]

The Euler chemical potential reconstructed independently from the pressure and energy density is

\[
\mu_B^{(E)}=\frac{\epsilon+P}{n_B}.
\]

For a cold beta-equilibrated table the implementation requires \(T=0\), \(Q_2=0\), \(Q_5=0\), and \(Q_6=Q_7\), within the stated numerical tolerances. It also requires the tabulated \(Q_3\) chemical potential to agree with the Euler value to relative tolerance \(5\times10^{-4}\). We recompute the output chemical potential from the Euler relation rather than silently retaining a discrepant input value.

The sound-speed convention is kept explicit:

\[
\frac{c_s^2}{c^2}=\left.\frac{dP}{d\epsilon}\right|_{T=0}.
\]

No assumption \(c=1\) is needed to interpret this dimensionless ratio.

## 3. Thermodynamic Maxwell conditions

For a system with one globally conserved baryon number, the Maxwell construction requires a common baryon chemical potential and common pressure. The local charge-neutrality condition is distinct from the global-neutrality constraint used in Gibbs constructions (Hempel et al. 2009; Bhattacharyya et al. 2010):

\[
\mu_B^{(L)}(\mu_t)=\mu_B^{(H)}(\mu_t)=\mu_t,
\qquad
P_L(\mu_t)=P_H(\mu_t)=P_t.
\]

The baryon density in either pure phase follows from the thermodynamic derivative

\[
n_B(\mu_B)=\frac{dP}{d\mu_B}.
\]

The zero-temperature Euler relation then fixes the energy density

\[
\epsilon(\mu_B)=\mu_B n_B(\mu_B)-P(\mu_B).
\]

The two phases generally have different densities at coexistence,

\[
n_L(\mu_t)<n_H(\mu_t),
\]

and therefore different energy densities,

\[
\epsilon_L=\mu_t n_L-P_t,\qquad
\epsilon_H=\mu_t n_H-P_t.
\]

Subtracting these equations gives the latent-energy identity

\[
\boxed{\Delta\epsilon=\epsilon_H-\epsilon_L
      =\mu_t\left(n_H-n_L\right)=\mu_t\Delta n_B.}
\]

This identity is not an optional diagnostic fit. It follows from the common \(P_t\) and \(\mu_t\), and the implementation rejects an output whose numerical residual exceeds \(10^{-10}\) in the report normalization.

In a \(P(\epsilon)\) representation the Maxwell construction is a constant-pressure interval with an energy-density jump. In the \(P(\mu_B)\) representation used by this module, the coexistence point is represented by two adjacent rows with identical \(P_t\) and \(\mu_t\), but distinct \(n_L\), \(n_H\), and \(\epsilon_L\), \(\epsilon_H\). There is no finite-density mixed-phase segment between them.

## 4. Thermodynamic interpolation

A tabulated EOS is discrete. To locate a transition between source knots, we therefore require a differentiable interpolant that preserves the First Law. We construct a cubic Hermite interpolant for pressure as a function of baryon chemical potential on each source interval.

For neighbouring points \((\mu_i,P_i,n_i)\) and \((\mu_{i+1},P_{i+1},n_{i+1})\), the interpolant matches

\[
P(\mu_i)=P_i,\quad P(\mu_{i+1})=P_{i+1},\quad
P'(\mu_i)=n_i,\quad P'(\mu_{i+1})=n_{i+1}.
\]

Its derivative is the interpolated baryon density:

\[
n_B(\mu_B)=P'(\mu_B).
\]

Because the derivative of a cubic is quadratic, \(n_B(\mu_B)\) is quadratic on each interval and \(dn_B/d\mu_B\) is affine. Defining

\[
\epsilon(\mu_B)=\mu_B n_B(\mu_B)-P(\mu_B)
\]

enforces the zero-temperature First Law identically within each interpolating segment:

\[
\frac{d\epsilon}{d\mu_B}
=n_B+\mu_B\frac{dn_B}{d\mu_B}-\frac{dP}{d\mu_B}
=\mu_B\frac{dn_B}{d\mu_B}.
\]

Consequently,

\[
\frac{c_s^2}{c^2}
=\frac{dP/d\mu_B}{d\epsilon/d\mu_B}
=\frac{n_B}{\mu_B\,dn_B/d\mu_B}.
\]

This interpolation is a controlled representation of the supplied barotrope, not a new microscopic model. The result depends on the resolution and quality of the source tables, so the interpolation report and the transition location should be subjected to EOS-refinement studies before observational inference.

## 5. Whole-interval stability and causality gate

The input branches must be mechanically stable and causal throughout the interpolated intervals, not merely at tabulation points. For \(\mu_B>0\) and \(n_B>0\), the conditions

\[
\frac{dn_B}{d\mu_B}>0,\qquad
0<\frac{c_s^2}{c^2}\leq1
\]

are equivalent to

\[
\frac{dn_B}{d\mu_B}>0,\qquad
\mu_B\frac{dn_B}{d\mu_B}-n_B\geq0.
\]

Within a Hermite interval, both \(n_B(\mu_B)\) and
\(\mu_B n'_B(\mu_B)-n_B(\mu_B)\) are quadratic polynomials. We evaluate their endpoint values and their analytic interior vertex, where applicable. Since \(n'_B\) is affine, its minimum is attained at an endpoint. The acceptance gate therefore checks the complete polynomial interval, rather than assuming that a fixed set of sample points proves the bounds.

We also evaluate 33 equally spaced points per interval to summarize the minimum and maximum sound speed in the report. These samples are diagnostic only; they are not the causality acceptance test. Pressure and energy density remain positive because \(P'=n_B>0\), \(d\epsilon/d\mu_B=\mu_B n'_B>0\), and their endpoint values are required to be positive.

## 6. Coexistence-root isolation

The two phases are evaluated on the common chemical-potential support. We retain their common (iY_q) and charge-fraction slice to prevent accidental interpolation between different composition-grid states. We form the union of both phase knot grids. On every resulting subinterval, each phase uses a single cubic Hermite segment, so

\[
\Delta P(\mu_B)=P_L(\mu_B)-P_H(\mu_B)
\]

is cubic, and

\[
\frac{d\Delta P}{d\mu_B}=n_L(\mu_B)-n_H(\mu_B)
\]

is quadratic.

We reconstruct this quadratic from its values at the two endpoints and midpoint, solves for every real stationary point inside the interval, and partitions the interval into monotonic pieces. A bracketed bisection is then applied to each sign-changing monotonic piece. Roots coincident with a partition endpoint are also retained, with a tolerance-based duplicate check. Hence, this isolates every root of the represented piecewise-cubic pressure difference, including a tangential root at a stationary point.

At every point we scale the pressure residual by the larger of the two local phase pressures, rather than by the maximum pressure anywhere in the input tables. Version 1 accepts exactly one pressure-equality root in the common domain. It rejects no-root and multiple-root cases instead of selecting a crossing by a hidden heuristic. A root is finally checked against the pressure-equality residual before the output table is allocated.

## 7. Output table and phase boundary

The output contains the low-density source knots below \(\mu_t\), two coexistence rows, and the high-density source knots above \(\mu_t\). The two coexistence rows have

\[
(P,\mu_B)_L=(P,\mu_B)_H=(P_t,\mu_t),
\]

but

\[
n_L<n_H,\qquad
\epsilon_L=\mu_t n_L-P_t,\quad
\epsilon_H=\mu_t n_H-P_t.
\]

The phase labels are set to 1 and 2 on the corresponding branches. The output is ordered by strictly increasing baryon density, has nondecreasing pressure and increasing energy density, and reports the retained row counts for each phase, the density jump, latent energy density, Euler closure residuals, and separate minimum/maximum sound-speed diagnostics for each pure phase.

The C1-C11 verifier now masks derivative stencils adjacent to any Iphase change, whether or not an eos.compo file exists. A finite-difference derivative must not cross the density discontinuity. The algebraic pressure and energy-density checks still see both coexistence rows, and the monotonic-pressure check permits the equal-pressure pair.

## 8. Strict v1 input contract

Both input tables must satisfy all of the following conditions:

1. At least three points per phase, finite arrays, positive density, pressure, and energy density.
2. Strictly increasing density, pressure, energy density, and Euler-derived baryon chemical potential within each phase.
3. Cold \(T=0\) metadata, zero entropy \(Q_2\), beta-equilibrium metadata, \(Q_5\simeq0\), and \(Q_6\simeq Q_7\).
4. Relative agreement between the CompOSE \(Q_3\) chemical potential and the Euler chemical potential at every knot no worse than \(5\times10^{-4}\).
5. No internal Iphase transition in either supplied phase branch.
6. No composition file and no usable composition arrays. V1 is a barotropic transition builder, not a composition reconstruction algorithm.
7. Nadd=0. Provider-specific extra columns have model-dependent meanings and are not interpolated or discarded silently.
8. Matching neutron and proton rest-mass conventions, matching CompOSE lepton-content flags, and matching (iY_q) and charge-fraction-grid metadata.
9. The caller explicitly passes assume_local_charge_neutrality=True.

The final condition is a physical assertion by the caller, not a property that can be independently verified from a composition-free one-dimensional table. A production workflow should retain the source model documentation proving local charge neutrality and record that provenance next to the generated transition.

## 9. API and command line

Python:

    from compose_eos import construct_maxwell_transition

    hybrid, report = construct_maxwell_transition(
        hadronic_phase,
        high_density_phase,
        assume_local_charge_neutrality=True,
        return_report=True,
    )

The explicit flag is required. Omitting it raises ValueError. The report states that composition was not preserved, the Seidov criterion was not evaluated, and phase-interface TOV support is not verified.

Native CLI, after creating the internal CSV/metadata directories:

    compose_eos_fortran maxwell-transition LOW_INTERNAL HIGH_INTERNAL OUTPUT_INTERNAL true

The CLI flag has the same meaning as the Python flag. Passing false is rejected by the Fortran kernel.

The output is an internal phase-labelled barotrope. It is not a complete composition-preserving CompOSE release. Both the Python write_compose wrapper and the native Fortran writer refuse to export a phase-discontinuous table as a complete CompOSE table when its phase composition is absent.

## 10. Analytic regression fixture

The test suite uses two exact cold barotropes:

\[
P_L(n)=n^2,\quad \epsilon_L(n)=939n+n^2,
\]
\[
P_H(n)=\frac14n^2,\quad
\epsilon_H(n)=940n+\frac14n^2.
\]

Their exact coexistence point is

\[
\mu_t=941\ {\rm MeV},\quad
P_t=1\ {\rm MeV\,fm^{-3}},
\]
\[
n_L=1\ {\rm fm^{-3}},\quad n_H=2\ {\rm fm^{-3}},
\]
\[
\epsilon_L=940\ {\rm MeV\,fm^{-3}},\quad
\epsilon_H=1881\ {\rm MeV\,fm^{-3}}.
\]

Thus \(\Delta n=1\ {\rm fm^{-3}}\) and
\(\Delta\epsilon=941\ {\rm MeV\,fm^{-3}}=\mu_t\Delta n\).
The cubic Hermite representation reproduces these quadratic \(P(\mu_B)\) branches exactly, so this fixture checks the construction independently of the legacy Python EOS generator. Additional tests reject composition-bearing input, missing local-neutrality assertion, disjoint chemical-potential domains, acausal input, a \(Q_3\)/Euler mismatch, Python export refusal and native CLI export refusal without composition.

## 11. What v1 does not claim

- **No Gibbs construction.** Gibbs equilibrium with global electric neutrality requires at least the phase charge densities, the independent charge chemical potential, and a volume fraction. These data are not reconstructible from two one-dimensional barotropes alone (Constantinou et al. 2023).
- **No composition reconstruction.** Electron and muon fractions, proton fraction, nuclear groups, and phase-specific charge constraints are absent from a composition-free input.
- **No crust generator.** Outer-crust nuclear species, electron gas, Wigner-Seitz Coulomb lattice, neutron drip and electron-capture sequence are not generated.
- **No microscopic high-density solver.** RMF, Skyrme, APR, or quark-matter model parameters are not inferred or solved by this routine.
- **No surface or finite-size physics.** Interface tension, Coulomb screening, pasta geometries and mixed-phase morphology are not included.
- **No Seidov verdict.** The energy-density jump is reported, but a stability threshold is not evaluated here.
- **No phase-interface TOV verification.** The output is not certified for a stellar sequence that crosses the density discontinuity. The existing stellar solver must be audited for the correct interface treatment and stable-branch selection before this output is used for mass-radius, tidal or oscillation inference.
- **No complete legacy-generator parity.** Full master-function smoothing/candidate selection, composition-aware crust-core matching and literature-specific EOS solvers remain separate tasks.

A successful build and analytic test do not remove these scope limits. They establish only that the implemented one-charge Maxwell operation satisfies its stated input contract and the tested thermodynamic identities.

## 12. References

1. S. Typel, M. Oertel and T. Klähn, *CompOSE: CompStar online supernova equations of state*, Physics of Particles and Nuclei **46** (2015), 633-664, and the updated CompOSE Reference Manual v3.01, [arXiv:2203.03209](https://arxiv.org/abs/2203.03209). The manual defines the eos.thermo quantities \(Q_1\) through \(Q_7\), their units and the eos.compo separation.
2. M. Hempel, G. Pagliara and J. Schaffner-Bielich, *Conditions for phase equilibrium in supernovae, protoneutron, and neutron stars*, Physical Review D, vol. **80**, article 125014, 2009, [doi:10.1103/PhysRevD.80.125014](https://doi.org/10.1103/PhysRevD.80.125014), [arXiv:0907.2680](https://arxiv.org/abs/0907.2680). This work distinguishes globally conserved charges from locally constrained charge fractions and clarifies the assumptions behind Maxwell and mixed-phase constructions.
3. A. Bhattacharyya, I. N. Mishustin and W. Greiner, *Deconfinement phase transition in compact stars: Maxwell vs. Gibbs construction of the mixed phase*, Journal of Physics G, vol. **37**, article 025201, 2010, [doi:10.1088/0954-3899/37/2/025201](https://doi.org/10.1088/0954-3899/37/2/025201).
4. *Framework for phase transitions between the Maxwell and Gibbs constructions*, Physical Review D, vol. **107**, article 074013, 2023, [doi:10.1103/PhysRevD.107.074013](https://doi.org/10.1103/PhysRevD.107.074013). This provides a further discussion of the distinct local and global charge-neutrality limits.

The exact CompOSE conventions above follow the reference manual. The Hermite closure and the whole-interval polynomial inequalities are derived directly from those definitions and are implemented in the native Fortran kernel.
