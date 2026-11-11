# NS_Cowling v4 physical model and Fortran verification guide

Date: 2026-11-08  
Repository: cosmosapjw-quantum/NS_Cowling  
Branch: research/fortran-eos-gap-fill-v1

## Architecture and scope

The repository contains CompOSE intake/verification, native Fortran EOS transformations, a Python process/API frontend and a Fortran TOV/Cowling solver. Python validates and serializes inputs, invokes executables and packages reports. The EOS calculations described here are performed by Fortran.

Four EOS operations are distinct contracts: fill_density_gap inserts a First-Law-consistent density bridge; match_crust_core constructs six smooth barotropic bridges; construct_maxwell_transition constructs a genuine one-charge coexistence point; reconstruct_master rebuilds a cold single-phase barotrope from F=f/n_B. The historic maxwell token in smooth matching is a near-pressure plateau without a density jump, not the genuine Maxwell constructor.

The stellar solver integrates a static spherical TOV background and polar-fluid modes in the Cowling approximation. It does not solve metric perturbations. Successful CompOSE intake does not prove that every interpolation is mechanically stable or causal between knots.

## CompOSE conventions and units

The CompOSE columns are
\[
Q_1=P/n_B,\quad Q_2=s/n_B,\quad Q_3=\mu_B/(m_nc^2)-1,
\]
\[
Q_4=\mu_q/(m_nc^2),\quad Q_5=\mu_l/(m_nc^2),\quad
Q_6=f/(n_Bm_nc^2)-1,\quad Q_7=\epsilon/(n_Bm_nc^2)-1.
\]
The EOS layer uses n_B in fm^{-3}, P, epsilon and f in MeV fm^{-3}, and m_n c^2 and chemical potentials in MeV. Q1 is in MeV; Q2,...,Q7 are dimensionless. See Typel et al., CompOSE Reference Manual v3.01 (2022), https://arxiv.org/abs/2203.03209.

The stellar integration uses SI: r in m, m in kg, pressure and energy density in J m^{-3}, n_B in m^{-3}, chemical potentials in J and angular frequency in s^{-1}. The conversions are 1 MeV=1.602176634e-13 J, 1 fm^{-3}=1e45 m^{-3}, and 1 MeV fm^{-3}=1.602176634e32 J m^{-3}. EOS expressions retain hbar c; natural units are not silently imposed.

## Cold master-function reconstruction

Let x=ln[n_B/(1 fm^{-3})] and F=f/n_B. For cold source input F_i=(1+Q6_i)m_n c^2. The native reconstruction defines
\[
Q_1=F',\qquad P=n_BQ_1,\qquad \epsilon=n_BF.
\]
It follows that
\[
\mu_H\equiv d\epsilon/dn_B=F+F',
\qquad D\equiv dP/dn_B=F'+F'',
\qquad c_s^2/c^2=D/\mu_H.
\]
mu_H is the enthalpy derivative along the barotropic path. It is not always identical to CompOSE Q3.

For Il=0 the repository verifier uses
\[
(1+Q_3)m_nc^2=(\epsilon+P)/n_B-Y_q\mu_q
=\mu_H-Y_qQ_4m_nc^2,
\]
where Y_q is charge_fraction_grid. Thus output Q3 subtracts Y_q Q4 for Il=0. For supported Il=1, source Q5 must be approximately zero and Q3 reduces to mu_H. The sound-speed denominator remains mu_H in either supported branch.

The input contract requires at least five points, allocated/shape-correct finite arrays, positive strictly increasing n_B, finite positive neutron rest energy, Il=0 or 1, finite charge-fraction metadata and a source flag asserting beta equilibrium. It enforces
\[
|T|\le10^{-12}\ {\rm MeV},\quad \max|Q_2|\le2\times10^{-8},
\quad\max|Q_5|\le10^{-8},\quad\max|Q_6-Q_7|\le2\times10^{-6}.
\]
Every source Q3 must satisfy the charge-aware CompOSE Euler identity to relative tolerance 2e-6. For Il=0 the source identity contains -Y_q Q4 m_n c^2. An equal-pressure pair across an Iphase boundary is rejected because a single-valued F(ln n_B) cannot cross a first-order density jump. N_add is 0 or 1; interpreting the extra column as sound speed requires explicit opt-in.

For smooth=False, native PCHIP interpolates F. For smooth=True, local weighted cubic least-squares values are formed at source knots before PCHIP. This is not exact parity with the historic global SciPy UnivariateSpline candidate-selection policy.

The represented interpolant must satisfy F>0, Q1>0, mu_H>0, D>0 and mu_H-D>=0 within the implemented tolerance. These impose positive energy and pressure per baryon, positive enthalpy derivative, mechanical stability and c_s^2/c^2<=1.

## Outer-crust electron-plus-lattice correction

The correction is conditional on supplied composition; it does not generate equilibrium nuclides. At zero temperature,
\[
n_e=Y_en_B,\quad k_{Fe}=(3\pi^2n_e)^{1/3},\quad x_e=\hbar ck_{Fe}/(m_ec^2),
\]
\[
P_e=\frac{(m_ec^2)^4}{24\pi^2(\hbar c)^3}
[x_e(2x_e^2-3)\sqrt{1+x_e^2}+3\,\operatorname{asinh}x_e].
\]
For small x_e the code uses the bracket expansion
\[
x_e^5(8/5-4x_e^2/7+x_e^4/3+O(x_e^6))
\]
to avoid cancellation; the coefficients were symbolically checked.

The lattice term is the single-species BCC Wigner-Seitz approximation
\[
\epsilon_{\rm lat}=-0.9(4/9\pi)^{1/3}\alpha Z^{2/3}\hbar ck_{Fe}n_e,\quad
P_{\rm lat}=\epsilon_{\rm lat}/3.
\]
The code uses supplied Z; it does not minimise over nuclides, calculate nuclear masses or electron-capture thresholds, or reproduce the complete BPS sequence. See Baym, Pethick and Sutherland (1971), https://doi.org/10.1086/151216.

The outer Q1 bridge is PCHIP-interpolated and integrated backwards from the core-side join:
\[
F(x)=F(x_j)-\int_x^{x_j}Q_1(u)\,du.
\]
The core-side Hermite slope is set equal to bridge Q1 at the join.

For t=(x-x_i)/h in [0,1], Q1 is cubic and its integral F is quartic:
\[
Q_1(t)=q_0+q_1t+q_2t^2+q_3t^3,\quad
F(t)=f_0+h(q_0t+q_1t^2/2+q_2t^3/3+q_3t^4/4).
\]
Thus F, mu_H and mu_H-D are quartic; Q1 and D are cubic. Quartic minima occur at endpoints or internal roots of the cubic derivative. The Fortran helper scales the derivative's quadratic coefficients, partitions [0,1] at critical points of the cubic, retains tangential roots, bisects sign-changing monotone intervals and evaluates all stationary points and endpoints. This replaces the former fixed 256-subinterval mesh. It is floating-point root isolation, not formal interval arithmetic. Native tests cover a flat minimum, two separated minima and a narrow negative pocket between all 257 points of the former mesh.

## Composition, dissolution and phase transitions

Composition is processed only with eos.compo provenance and finite required arrays; composition-free input does not receive invented fractions. Supplied nuclei determine Y_q=Y_p+(Z/A)Y_h, with charge neutrality Y_e+Y_mu=Y_q. Lepton partition above muon threshold is solved by scalar bisection. Optional dissolution repair is a bounded legacy-data heuristic, not a dynamical or microscopic free-energy calculation.

The six match_crust_core methods are smooth bridges. The separate Maxwell kernel assumes local charge neutrality in each supplied phase and one globally conserved baryon charge:
\[
P_L(\mu_t)=P_H(\mu_t),\quad \mu_B^L=\mu_B^H=\mu_t,\quad
n_B=dP/d\mu_B,\quad \epsilon=\mu_Bn_B-P.
\]
Therefore Delta epsilon=mu_t Delta n_B. Gibbs mixed phase requires global neutrality, an independent charge chemical potential, phase charge densities and a volume fraction; those are not available from a generic one-dimensional barotrope. Gibbs, Seidov onset stability and phase-interface TOV/Love/mode conditions are not implemented. See https://arxiv.org/abs/0907.2680.

## TOV and Cowling solver

With signature (-+++),
\[
ds^2=-e^{2\nu}c^2dt^2+e^{2\lambda}dr^2+r^2d\Omega^2.
\]
The static background equations are
\[
m'=4\pi r^2\epsilon/c^2,\quad P'=-(\epsilon+P)\nu',
\quad \nu'=\frac{G[m+4\pi r^3P/c^2]}{c^2r^2[1-2Gm/(rc^2)]}.
\]
The model is static, spherical, nonrotating and isotropic. The code integrates in logarithmic pressure; its surface is a finite positive-pressure cutoff with Schwarzschild time normalisation, not an exact vacuum boundary.

For polar perturbations proportional to exp(i omega t), the implemented l=2 system evolves W,V:
\[
Q=(\omega r/c)^2e^{\lambda-2\nu}V+\nu'W,\quad
W'=Q/\chi-6e^\lambda V,\quad V'=2\nu'V-e^\lambda W/r^2,\quad
\chi=dP/d\epsilon.
\]
At the centre W scales as r^3 and V as -r^2/2. The terminal residual uses the Lagrangian pressure perturbation. Cowling neglects Eulerian metric perturbations. These are fluid-mode frequencies of the approximation, not full-GR quasinormal modes or gravitational-wave damping rates. Rotation, magnetism, elastic crust, superfluidity, frozen-composition response and sharp phase-interface transfer conditions are not included.

## Build, tests and evidence

    python -m pip install -e '.[test]'
    bash fortran_eos/build.sh build/compose_eos_fortran
    python -m pytest -q python_frontend/tests pipeline/tests

CI additionally compiles/runs fortran_eos/tests/test_polynomial_extrema.f90, checks sha256sum -c SHA256SUMS_v4.txt, creates a Git archive of the tested source and uploads it as a GitHub Actions artifact.

Literature-supported claims concern CompOSE, BPS and phase-equilibrium conventions. Derived claims are the master-function identities, sound-speed relation, outer-bridge polynomial degrees and Maxwell latent-energy relation. Numerically checked claims include the small-x_e expansion and adversarial polynomial tests. Implementation-verified claims require exact-head successful build, tests, checksum verification and artifact. Historic global-spline parity, universal dissolution, full BPS species generation, microscopic RMF/Skyrme/APR/quark solvers, Gibbs, Seidov, phase-interface stability and full-GR metric perturbations remain unresolved or unimplemented.

## References

1. Typel et al., CompOSE Reference Manual v3.01 (2022), https://arxiv.org/abs/2203.03209.
2. Baym, Pethick and Sutherland, Astrophysical Journal 170 (1971), https://doi.org/10.1086/151216.
3. Hempel et al., phase equilibrium in supernovae/protoneutron stars/neutron stars, https://arxiv.org/abs/0907.2680.
