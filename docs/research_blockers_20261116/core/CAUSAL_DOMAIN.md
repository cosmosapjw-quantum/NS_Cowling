# Native core causal-domain result (2026-11-16)

Baseline: a49aaf9f82ddfda484e0f965aefdfe73aa7918f2. Selected work: R2–R9 and C2–C8, restricted to the core acausality blocker. No core functional, coefficient, phase construction, or raw EOS output was changed.

## Decision and numerical result

APR, KDE0v1, and SLy4 exceed the relativistic causal bound on their unmodified cold beta-equilibrated branches. This is a physical high-density limitation of the declared nonrelativistic functionals, not a finite-table differentiation artifact. An independently implemented analytic Hessian agrees with derivatives of the native equilibrium solver. The diagnosed numerical edge effects (muon onset and BBB interpolation knots) are separate and now explicit diagnostic statuses.

| Native model | Numerical crossing bracket, n in fm^-3 | Suggested displayed crossing |
| --- | --- | --- |
| APR_APR | [0.8689531515166163, 0.8689531516283750] | 0.86895315 |
| RG_KDE0v1 | [1.4109658931940792, 1.4109658933803437] | 1.41096589 |
| RG_SLY4 | [1.2650056713819504, 1.2650056715682150] | 1.26500567 |

These are numerical sign brackets with width at most 2e-10 fm^-3, not certified interval bounds. Fourth-order native derivatives with six relative steps 0.002, 0.001, 0.0005, 0.0002, 0.0001, 0.00005 produce the same brackets. At the crossing the maximum discrepancy from the analytic Hessian is 8.19e-12 in c_s²/c². The operational boundary buffer is 1e-8 fm^-3, larger than observed method/step differences; it is not a statistical uncertainty. Published coefficient rounding and model uncertainty are not included. A coarse scan locates one causal crossing per model; this is not a theorem excluding arbitrarily narrow unsampled features.

BBB and GM1 have no acausal point in the new 41-point interior scans: maximum c_s²/c² is 0.83308744 on n=0.081…1.499 for BBB, and 0.75629727 on n=0.081…0.999 for GM1. This finite-range statement does not prove global causality or the validity of omitted degrees of freedom.

## Definition and derivation

Assume zero temperature, neutrino-free beta equilibrium, positive baryon density, charge neutrality, and homogeneous non-spin-polarized matter. Leptons are free relativistic electrons/muons (GM1 retains its declared npe convention). The APR coexistence calculation assumes bulk global neutrality with a common lepton gas; it has no finite-size Coulomb or surface correction. Retain c and ħ. Numerical nucleon rest energies are m_n c² and m_p c² in MeV; epsilon and P have units MeV fm^-3.

Let N(n,y) be nuclear energy density including rest energy, y=n_p/n, and L(q) the lepton energy minimized over electron/muon partition at total charge density q=ny. Write F(n,y)=N(n,y)+L(ny). With mu=mu_e and susceptibility chi_l=d(n_e+n_mu)/dmu,
L'=mu and L''=1/chi_l, where

chi_l = sum_l mu sqrt(mu²-m_l² c⁴)/(pi² (ħc)³), for mu>m_l c².

Beta equilibrium is F_y=0 and F_yy>0. Implicit differentiation gives y'=-F_ny/F_yy. Therefore epsilon_eq'=F_n=mu_n and

epsilon_eq'' = F_nn - F_ny²/F_yy,
P = n epsilon_eq' - epsilon_eq,
c_eq²/c² = dP/depsilon_eq = n epsilon_eq''/mu_n.

The implemented nuclear-to-total Hessian map is

F_nn=N_nn+y²/chi_l,
F_ny=N_ny+mu+ny/chi_l,
F_yy=N_yy+n²/chi_l.

N_nn has units MeV fm³, so n F_nn/mu_n is dimensionless. The response at fixed proton fraction, with the lepton partition still relaxed, is n F_nn/mu_n and exceeds or equals the equilibrated response. Thus an equilibrated value above one cannot be cured by freezing proton fraction. Which response governs a particular oscillation depends on reaction time scales; the reported threshold is the cold equilibrium EOS derivative.

The analytic implementation builds the published nuclear functional as separable terms R(n)G(y), differentiates R and G through second order, and applies this Schur complement. It does not differentiate a CompOSE table or reuse the native pressure formula. Native finite differences independently evaluate both (dP/dn)/(depsilon/dn) and n(dmu_n/dn)/mu_n. The first-law residual checks depsilon/dn=mu_n.

### APR two-charge Gibbs derivative

The native mixed interval is n=0.20439495489345874…0.23738424893920132 fm^-3. The acausal crossing is far above it, on the pure high-density branch.

For each nuclear phase use conserved densities (n,q=n_p), chemical potentials (mu_B=mu_n,mu_Q=mu_p-mu_n), and the susceptibility K_i=(Hessian_(n,q) N_i)^(-1). With volume fraction x, K=(1-x)K_L+x K_H+diag(0,chi_l). Let delta n=n_H-n_L and delta q=q_H-q_L. Pressure equality gives delta n dmu_B+delta q dmu_Q=0. Differentiating total baryon density and neutrality yields

dmu_B/dn = (delta q)² /
[(delta q)² K_BB - 2 delta q delta n K_BQ + (delta n)² K_QQ].

Hence c_eq²/c²=(n/mu_B)dmu_B/dn. This is a two-charge Gibbs response, not a constant-pressure one-charge Maxwell response. The native module computes it from independent phase Hessians. Twelve checks at three mixed densities and four steps give maximum analytic/finite-difference discrepancy 3.85e-11. Phase endpoints have distinct one-sided derivatives; a centered stencil crossing them is rejected.

## High-density mechanism

For a regular asymptotic branch epsilon=m c² n+A n^gamma, A>0 and gamma>1, the exact expression is

c_s²/c² = gamma(gamma-1) A n^(gamma-1) / [m c²+gamma A n^(gamma-1)],
lim_(n→infinity) c_s²/c² = gamma-1.

Wolfram Language directly verified this identity, the equilibrium Schur complement, the density/composition Hessian transformation, and the two-charge Gibbs tangent formula; algebraic residuals are zero. This is ordinary model analysis, not a new strong no-go theorem.

In both Skyrme models the leading nuclear term is K n^(8/3)[c0 F_(5/3)(y)+c1 F_(8/3)(y)], K=(3/5)(3pi²)^(2/3), F_a=y^a+(1-y)^a. The coefficients (c0,c1) in MeV fm^5 are (28.7316595,29.073189) for KDE0v1 and (32.47299,49.31277) for SLy4. Both are positive; strict convexity gives the unique leading minimum y=1/2. Lower terms scale as n^(2+alpha), n², n^(5/3), and n^(4/3), with alpha<2/3. On the regular equilibrated asymptotic branch y→1/2 and c_eq²/c²→5/3. This analytic mechanism agrees with the published nonrelativistic Skyrme causality analysis [4].

APR requires more care: for fixed nonzero interior y its potential contains 76.4 y(1-y)n^4, while the pure neutron limit contains 384 n^3. The exponential high-density corrections decay. Fixed-interior-composition asymptotics therefore yield 3; the formal constrained high-density minimum eventually goes to the y=0 boundary and yields 2. This limiting discussion is outside the code's APR support (n≤2 fm^-3) and is not used to locate the measured crossing or assert the smooth continuation of its interior beta root. The actual finite-density Hessian calculation establishes the acausal region needed here.

## Source and implementation traceability

1. [Schneider et al., PRC 100, 025803 (2019)](https://arxiv.org/html/1901.09652v1), Eqs. (10)–(14), Table 1, and Appendix C: APR energy functions and all 21 coefficients. Their g_i are energy-density terms; the native implementation stores g_i/n². This factor is explicitly included in the new radial derivatives. Their published simulation phase treatment is not used as evidence that the project's two-charge Gibbs and a Maxwell construction are identical.
2. [Agrawal, Shlomo and Au, PRC 72, 014310 (2005)](https://arxiv.org/html/nucl-th/0505071v1), Section II and Table 3 second column KDE0(v1): exact KDE0v1 coefficients including alpha=0.1673. This column must not be substituted by the KDE0 fit in Table 5.
3. [Douchin and Haensel (2000)](https://arxiv.org/html/astro-ph/0006135v1), Table 1: SLy4 published parameter precision; bulk functional obtained by removing gradient, Coulomb and spin-orbit contributions in uniform matter.
4. [Constantinou and Prakash (2017)](https://arxiv.org/html/1702.06952v1), degenerate-limit discussion and Eq. (28): nonrelativistic high-density acausality and the n^(8/3) mechanism. Its illustrative pure-neutron thresholds are not substituted for this calculation's npe-muon beta-equilibrated thresholds.

Source reads were performed through Parallel on 2026-11-16. Provider acquisition dates and HTML-render dates are not publication dates. Native source coefficients were directly compared with the relevant parameter tables. Source fidelity, implementation agreement, scientific validity and CompOSE-table agreement remain distinct. No CompOSE value was used in any reported derivative or causal crossing. Historical Python roots fall inside the native brackets; this checks independent implementation of the same construction, not independent microphysics. BBB remains a rounded published-energy interpolation, not a new BHF solve.

## Numerical edge effects and interfaces

At APR n=0.117325 fm^-3 a relative step 1e-4 crosses muon onset and a smooth fourth-order stencil produces a 6.0e-7 sound-speed error. The new status rejects this stencil. At relative steps 5e-5 through 6.25e-6 the analytic disagreement is 1.7e-13…7.3e-12. At BBB n=0.40005 an initial stencil crossed a published-energy PCHIP knot; curvature is branch-sided, and a symmetric derivative is not meaningful there. The original failed evidence is retained as initial_samples.csv and run_initial_failed.log. These effects do not explain the high-density causal crossings.

Public additive APIs in compose_eos_core_domain:
- core_sound_speed_fd(model,n,rel,cs,cs_mu,firstlaw,status)
- core_sound_speed_analytic(model,n,cs,frozen,hyy,beta_residual,status), APR/Skyrme pure phases
- core_sound_speed_gibbs(n,cs,status), interior APR mixed phase
- core_causal_bracket(model,left,right,rel,tol,lo,hi,status)

Only status='ok' validates outputs. Failed calls initialize output scalars to zero; these zeros must never be reported as physical sound speeds. At a phase endpoint choose the requested one-sided branch or report unavailable. The symmetric BBB stencil also rejects its source support endpoints and interpolation knots. Finite-difference status can be resolved by reducing the step only when all stencil points lie on a smooth local branch. Analytic pure-phase evaluation can be used at muon onset because the susceptibility is continuous there.

## Validation and reproducibility

Run: bash fortran_eos/tests/core_domain_run.sh. The isolated test build requires only gfortran and the types/native-core/domain modules; it does not modify the shared production build. The actual compiler invocation used -std=f2018 -O2 -Wall -Wextra -ffree-line-length-none -fcheck=all, one process, and OMP/OPENBLAS/MKL thread counts one. Exit status is 0 and FAILURES=0. The new module has no compiler warning; two unchanged native-module warnings concern real equality and an unused variable.

Artifacts: data/research_blockers_20261116/core/CORE_DOMAIN_SUMMARY.json, core_domain_samples.csv (205 interior states), core_domain_crossings.csv (18 converged brackets), core_domain_gibbs.csv, core_domain_muon_convergence.csv, wolfram_verification.json, compile.log, run.log and exit_status.txt. Source signatures are in the JSON report. Existing raw native core source SHA256 remains a5fd20940c286e1390ad95a5542cc1792ad1a9bb5a7b6a3ebea12b99a881d3c6.

The useful domain is the intersection of source support, numerical solvability, thermodynamic stability and causal range. The causal cutoff alone does not certify microscopic validity below it. Above the crossing the unchanged source model remains reproducible, but its use as a causal relativistic hydrodynamic EOS is invalid.

An optional separately named continuation can impose P=P_s+beta(epsilon-epsilon_s), 0≤beta≤1, with epsilon(n)=[(epsilon_s+P_s)/(1+beta)](n/n_s)^(1+beta)+(beta epsilon_s-P_s)/(1+beta). This is a changed physical model with an explicit matching density and uncertainty, not a correction to APR/SLy4/KDE0v1. It was not implemented or silently selected.

## Integration-review boundary correction

The exported diagnostics initially delegated density acceptance to solve_core_state, which does not enforce core_model_bounds. Consequently a direct API caller could evaluate beyond the declared model domain even though the public generator wrapper was protected. This implementation-boundary bug is now repaired in the diagnostic module: FD, analytic and Gibbs validate finite center density against core_model_bounds; FD also validates the whole [n-2h,n+2h] stencil. The bracket routine inherits these checks. Analytic evaluations accept the inclusive model endpoints, while a symmetric derivative cannot extend beyond them. The numerical derivative formulas and underlying EOS remain unchanged.

The dedicated core_domain_bounds.f90 driver was compiled and run once with the same compiler/thread settings: 20 checks, FAILURES=0, exit 0. Tests cover SLy4 n=1.5, GM1 n=1.1, APR n=1.5, lower/upper stencil overruns with valid centers, NaN/infinity, Gibbs bounds, valid inclusive analytic endpoints, a valid Gibbs state, an APR derivative agreement and the unchanged SLy4 crossing. The broad campaign was not repeated. BOUNDARY_FIX.json binds the initial module SHA256 acec428c5283d86d1015ed4f8ecc695f29eb9d87350c6cb2088b6bd3665cc03f to corrected SHA256 683ef46c86fe72a9b535dc22d0511e15c93905b3879311ef8abc3e9658930c50. CORE_DOMAIN_SUMMARY.json preserves initial_validation_artifacts and records the current artifacts and targeted correction separately.
