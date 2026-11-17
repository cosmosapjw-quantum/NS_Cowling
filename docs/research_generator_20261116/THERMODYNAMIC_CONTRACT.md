# Composition and density-domain contract

The research target comprises five core models and thirteen crust models, each regenerated on its own supported density range with explicit composition. The source archive, chosen functional, species, rest masses, phase prescription and numerical representation define comparison identity.

For cold uniform matter write epsilon(n_n,n_p,n_e,n_mu), including all rest energies. Baryon conservation gives n=n_n+n_p and neutrality gives n_p=n_e+n_mu. Interior stationarity gives mu_n-mu_p=mu_e and mu_mu=mu_e when muons are populated. The empty-muon branch satisfies mu_e<=m_mu*c^2. Along a smooth stationary neutral branch,

    d epsilon_eq/dn = mu_n
    P = sum_i n_i mu_i - epsilon_eq = n mu_n - epsilon_eq
    c_eq^2/c^2 = dP/d epsilon_eq = n epsilon_eq''/epsilon_eq'

These identities were independently simplified with Wolfram Language. With composition coordinates x at a local minimum, epsilon_x=0 and positive composition Hessian H=epsilon_xx, differentiating equilibrium yields

    dx/dn = -H^(-1) epsilon_xn
    epsilon_eq'' = epsilon_nn - epsilon_nx H^(-1) epsilon_xn.

The correction is nonnegative when H is positive definite. Frozen-composition and weak-equilibrium sound speeds therefore answer different physical perturbation questions. Muon onset and phase coexistence require piecewise treatment. A causal beta-equilibrium branch alone does not certify spin, finite-wavelength or frozen-composition stability.

Units in generated tables are n in fm^-3, mu and particle rest energies m*c^2 in MeV, epsilon and P in MeV fm^-3. hbar*c is explicitly stored in MeV fm. Report c_s^2/c^2 as dimensionless. The implementation's MeV mass symbols stand for rest energies.

## Revised use of diagnostics

The previous full-table C10 FAIL flag is retained as a diagnostic of the tested full domain. The 2026-11-16 user instruction requires locating supported regions before model-level exclusion. Every original row remains recorded. The implementation now reports causal/stable sampled segments, first violating brackets, composition availability and Euler residuals. Coarse grid spacing is recorded as continuity uncertainty. A sampled secant estimates an interval average and cannot locate an exact continuous crossing. Smooth functional crossings are independently refined with derivative step checks.

An intrinsic high-density model limitation, a numerical differentiation defect, a wrong source convention and a wrong model parameterization have different corrective actions. Original-provider table comparison and independent functional regeneration identify the attribution. Increasing tolerance to reclassify a failure is not part of this change.

## Crust equations that remain to implement

The Scurto source GDFM is a three-meson density-dependent model with x=n/n_sat. Its legacy four-meson surrogate is preserved as a separately identified historical implementation. The author parameter vectors provide the bulk input; the finite-nucleus surface fit and complete cell variation determine the crust. The AME candidate dataset can support an explicitly documented fit experiment, with source-equivalence assessed from the fit definition and resulting residuals. Author-authenticated labels are evidence of fidelity rather than an authorization requirement for implementation.

For the GPPVA Thomas-Fermi family, the same energy functional must include spatial meson and Coulomb gradient terms. At a fixed mean baryon density, solve the coupled meson-field equations, local neutron/proton Thomas-Fermi relations, electron background and beta/charge constraints with zero normal field derivatives at the cell center and boundary. Minimize the cell energy over its size and each allowed geometry. Droplet, rod and slab solutions, followed by inverse geometries where supported, must be compared on the same thermodynamic normalization. The finite-range gradient functional supplies surface energy; importing the Scurto fitted sigma as a replacement would change the target model.

The full crust completion criteria include source bulk-parameter agreement, composition conservation with per-species/group semantics, converged fields and energy, density/grid/multistart checks, outer-crust provenance, and original-table differences on a common density domain. Existing extraction/readback results provide preserved reference data for these tests.
