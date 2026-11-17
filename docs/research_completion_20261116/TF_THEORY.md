# Finite-range RMF Thomas–Fermi derivation and source audit

This note covers independent nucleonic inner-crust reconstruction for GPPVA_DD2, DDME2, TW, FSU2, FSU2R, FSU2H, NL3wrL55 and TM1e. It does not establish an outer-crust mass model, a core–crust join, complete density coverage, or exact identity to the CompOSE producer's unpublished numerical choices. Source table agreement is a diagnostic. The primary admission conditions are the correct functional, thermodynamic and composition constraints, stationary field/profile equations, and demonstrated spatial/radius/density convergence within a declared physical domain.

## Source chain and approximation

Grill et al., *Phys. Rev. C* **90**, 045803 (2014), arXiv:1404.2753, defines the crust TF framework and compares the minimum-energy Wigner–Seitz (WS) geometry with homogeneous matter. Its implementation references Avancini et al. (2008, 2009, 2010) and Maruyama et al. (2005). Avancini 2009, arXiv:0812.3170, Eqs. (53)–(61), explicitly provides the off-shell TF functional and local n,p,e equilibrium. Avancini 2010, arXiv:1010.3644, Eqs. (9)–(11), explicitly gives local electron screening and density-dependent (DD) rearrangement. Maruyama 2005, arXiv:nucl-th/0503027, states smooth WS connections with zero field gradients at cell boundaries.

The 2012 erratum included in Avancini arXiv:1010.3644v3 corrects the CP temperature-dependent surface-tension implementation and its use of global proton fraction. It does not change the zero-temperature TF equations; its only TF tabular corrections are two finite-temperature pressure misprints.

The source TF approximation includes **local electrons**. Uniform electrons appear in a distinct coexisting-phases approximation and in some other authors' TF models; they must not silently replace the screened source functional. The source crust calculation is nucleonic npe matter. FSU2H has hyperonic core applications, but this does not introduce hyperons into its low-density nucleonic crust. A no-muon crust implementation must still report whether any local electron kinetic Fermi energy exceeds the muon rest energy; above that threshold the npe approximation needs extension. No phenomenological surface-tension fit is added to finite-range TF: meson gradients generate the interface contribution.

The nonrelativistic gradient corrections, shell, pairing and exchange effects absent from this source TF approximation are not manufactured here. They limit fidelity to microscopic nuclei, particularly near neutron drip; a low-density TF continuation cannot be relabeled as a calibrated nuclear mass calculation.

## Units, fields and energy

Write H = ħc in MeV fm, masses as rest energies M in MeV, coordinates in fm, and densities in fm^-3. Raw meson fields φ, w and b have units MeV. The scalar convention is M_i* = M_i − g_s(n)φ, and the isovector vertex is g_r(n) τ_3 / 2 with τ_3(p)=+1, τ_3(n)=−1. Define n=n_n+n_p, n_3=n_p−n_n, q=n_p−n_e and V=e A_0 in MeV. The electromagnetic constant is e²=αħc in MeV fm. U(φ) and W(w,b) have units MeV^4.

For spin-two species, k_F=H(3π²n_i)^(1/3) and
ε_kin,i=(π²H³)^−1 ∫_0^kF k² sqrt(k²+M_i*²) dk.
For electrons M_e*=m_e c². The scalar density is ∂ε_kin/∂M*. The off-shell energy density is

ε = < Σ_i ε_kin,i
 + [ U + H²|∇φ|²/2 − W − H²(|∇w|²+|∇b|²)/2 ]/H³
 + g_w w n + g_r b n_3/2
 + qV − |∇V|²/(8πe²) >.                                      (1)

Scalar gradient energy has a positive sign; vector and electromagnetic gradient terms have negative signs before their constraint fields are eliminated. Minimizing every raw field as if all terms were positive is incorrect. On the stationary Poisson solution, the Coulomb contribution becomes E_C=<qV>/2=<|∇V|²>/(8πe²)≥0. On stationary meson solutions, retaining Eq. (1) avoids accidental omission of nonlinear cross terms. In a homogeneous cell, quartic vector monomials contribute three times their W coefficient to the eliminated vector energy, whereas they enter pressure with their original W coefficient.

The field equations follow directly by variation:

−H² ∇²φ + U_φ = H³ g_s n_s,
−H² ∇²w + W_w = H³ g_w n,
−H² ∇²b + W_b = H³ g_r n_3/2,
−∇²V = 4πe² q.                                                 (2)

All terms in each meson equation have units MeV³; Poisson terms have units MeV fm^-2. Derivatives of DD vertices do not enter Eq. (2), because local nucleon densities are held fixed when varying fields.

## Density-dependent rearrangement and composition

For local DD vertices, the common rearrangement potential is

Σ_R = g_w'(n)wn + g_r'(n)b n_3/2 − g_s'(n)φ n_s,                 (3)

where primes have units fm³. The local functional derivatives are

μ_n(r)=E_Fn*+g_w w−g_r b/2+Σ_R,
μ_p(r)=E_Fp*+g_w w+g_r b/2+Σ_R+V,
μ_e(r)=E_Fe−V.                                                  (4)

This is Avancini 2009 Eqs. (59)–(61) and Avancini 2010 Eq. (11), transformed to the explicit raw-field units above. Using the mean baryon density inside g(n) would define a different functional; omitting Σ_R corrupts both equilibrium and the EOS pressure.

At fixed mean baryon density n_B and global neutrality, minimize ε with
<n_n+n_p>=n_B and <n_p−n_e>=0.
The Lagrangian ε−μ_B<n_n+n_p>−μ_Q<n_p−n_e> gives
μ_n=μ_B, μ_p=μ_B+μ_Q, μ_e=−μ_Q, hence μ_n=μ_p+μ_e.
At a vanishing species density the correct condition is the corresponding one-sided KKT inequality, not an equality imposed on an absent species. Local μ_e is spatially constant at equilibrium while its kinetic part E_Fe=μ_e+V varies.

Mean Y_n,Y_p,Y_e and charge residuals are well-defined for every geometry. A finite cluster A or Z for a rod/slab needs a declared longitudinal normalization and cannot be inferred from spherical-cell formulas. For a sphere, total WS-cell counts 4πR³<n_i>/3 are unambiguous but include dripped neutrons; they are not automatically the cluster A. Maruyama Eq. (15) defines a diagnostic radius from proton-density moments averaged along the radial direction, not a unique diffuse cluster/gas decomposition. Report profiles and explicitly named cell counts rather than inventing source-equivalent nuclear abundances.

## Geometry, boundaries and radius variation

For dimensionality d=1,2,3,
<f> = d/R^d ∫_0^R r^(d−1) f(r) dr
and ∇²_d f=r^(1−d) ∂_r[r^(d−1)∂_r f].
Regularity requires f'(0)=0. At the WS edge, impose zero normal derivatives for all meson and electromagnetic fields. Global charge neutrality is the solvability condition of the Neumann Poisson problem; fix the otherwise arbitrary additive potential with <V>=0.

The five standard WS topologies are d=3 droplet/bubble, d=2 rod/tube and d=1 slab. Dense-centre and hollow-centre profile seeds cover both orientations. A topology label describes the resulting profile, not merely its initial seed; a relaxed uniform cell has no droplet identity. All branches must be compared with the homogeneous state at the same mean density, using the same model and particle content.

At fixed dimensionless profiles x=r/R, let G_s,G_w,G_r be the positive magnitudes of the three meson gradient energies. The envelope theorem after field elimination gives the exact scale derivative

dε/dlogR = −2G_s + 2G_w + 2G_r + 2E_C.                           (5)

The meson gradients scale as R^-2, with their respective off-shell signs. The neutral Poisson inverse scales as R², hence E_C scales as R². At an interior stationary radius, G_s−G_w−G_r=E_C. This relation is derived from Eq. (1), not an imported sharp-interface virial rule. A radius at an artificial search bound is not an interior stationary solution. A uniform state has a redundant radius.

For a fully optimized zero-temperature neutral beta-equilibrium branch, P=n_B μ_B−ε. Independently check P=n_B² d(ε/n_B)/dn_B on smooth branches using reoptimized nearby densities. Comparing the two expressions before profile or radius convergence can misdiagnose a numerical error as model physics. At branch crossings, one-sided derivatives/phase treatment are required.

## Eight-model convention audit

| Model | Source convention to preserve | Transformation used by raw-field functional |
|---|---|---|
| DD2 | Typel et al. 2010; distinct proton/neutron masses; native rho vertex without 1/2 | Set g_r=2 g_r,native, including its density derivative; use local n |
| DDME2 | Lalazissis et al. 2005; M=939 MeV; native rho normalization | Same factor-of-two rho transformation |
| TW | Typel–Wolter 1999; M=939 MeV; density-dependent rational/exponential vertices | Same factor-of-two rho transformation |
| FSU2 | Published FSU parameters are g²; scalar nonlinear powers are powers of g_s φ | U=m_s²φ²/2+κ(g_sφ)³/6+λ(g_sφ)⁴/24 |
| FSU2R | Revised source fit has Λ_ω=0.045, distinct from the older fit | Square-root tabulated g²; preserve revised set identity |
| FSU2H | Same nucleonic field form; hyperonic couplings are irrelevant below onset | Nucleonic crust only; no hyperon population inferred |
| NL3wrL55 | Providência et al., arXiv:1203.3052; κ/M=4.384, λ=−173.31 apply to unscaled raw scalar field | U=m_s²φ²/2+κφ³/6+λφ⁴/24; do not multiply these coefficients by powers of g_s |
| TM1e | Official guide uses M*=M+g_sσ; g_2 in fm^-1; M=938 MeV | φ=−σ, so cubic coefficient becomes −g_2 H; quartic unchanged |

Finite source parameter precision can produce small saturation-property offsets. Record them and test convergence; do not refit source constants to improve a table score. Negative scalar quartic coefficients imply that some parameterizations are effective local models, not globally bounded polynomial theories. Solve and report the physically connected positive-effective-mass branch over the declared domain; a spurious high-field runaway is not admissible simply because it has a lower formal polynomial energy.

## Independent validation and remaining scope

Required implementation checks are: field-equation residuals; weighted finite-volume self-adjointness and zero-flux constant null mode; Coulomb positivity and gauge invariance; local energy directional derivatives against all species chemical potentials (including DD terms); analytic log-radius derivative against finite differences; homogeneous reduction for every model; density, charge and lower-bound KKT residuals; multistart/radius-domain and mesh convergence; pressure from independently reoptimized density differences; and uniform-phase comparison. Numerical failure, source/domain limitations, and physically allowed but acausal or unstable continuation are distinct findings.

A finite list of converged cells is not proof of a global variational minimum, a complete density grid, or reproduction of original producer conventions. A model's acausal high-density region is an applicability boundary to identify, not automatic evidence that its independently implemented low-density branch is wrong. The original nucleus dataset, outer-crust convention and source table normalization remain separate evidence obligations.

Reproducible primary-source URLs, acquired SHA-256 values and review evidence are recorded in tf_literature.json and the tf_review directory. Source acquisition is literature-supported evidence; Eqs. (1)–(5) include explicit unit conversion and independently derived consequences. Execution results must be read from the separate numerical review receipt.

## Executed independent review

The retained functional review checks 18 kinetic-integral cases by direct quadrature and 24 nonuniform states (all eight models in all three dimensions). Maximum relative kinetic error was 1.64e-11; the maximum absolute species directional-derivative difference was 6.67e-8 MeV; the log-radius derivative difference was 5.49e-11 MeV fm^-3. Homogeneous reduction, weighted Poisson self-adjointness, constant null mode, Coulomb positivity and the field-energy identity passed. This is a finite-case functional review, not a proof over all inputs.

An independently optimized DDME2 spherical branch at n_B=.03 fm^-3 gave P=.0793529054 and .0793826234 MeV fm^-3 on 32 and 64 radial cells. The energy per baryon changed by .0003173 MeV; the 64-cell maximum chemical KKT residual was .0002504 MeV. Reoptimizing neighboring densities with steps 3e-5 and 1e-5 fm^-3 reproduced the Euler pressure with relative differences 2.56e-7 and 2.77e-8. This checks one branch, without asserting a geometry-global minimum or full-grid convergence.

Grill2014 Table4 at the same density gives pressure lower by .409% and energy density higher by .06381% than this result, after converting its printed values with hbarc. The approximately .60 MeV/baryon energy offset needs a source normalization and physics investigation. It is not silently corrected, used to refit parameters, or converted into a blanket rejection.

Primary links: [Grill2014](https://arxiv.org/html/1404.2753), [Avancini2009](https://arxiv.org/html/0812.3170), [Avancini2010 plus erratum](https://arxiv.org/html/1010.3644v3), [Maruyama2005](https://arxiv.org/html/nucl-th/0503027).
