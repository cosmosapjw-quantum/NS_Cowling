# Reconstructed Negele–Vautherin source models

The three native CompOSE directories are `printed/`, `d03x100/` and `uniform_beta/`. The first two describe compression of a neutral Wigner–Seitz cell of fixed neutron and proton numbers. The third minimizes the homogeneous energy with respect to proton fraction. These are generated comparison inputs. A physical cold-catalyzed inhomogeneous ground state would additionally require orbital and cell-composition minimization.

## Sources and units

Negele and Vautherin (1972), Phys. Rev. C **5**, 1472, Eq. (3.10–11), Table I on p. 1487, supplies the B, C and D coefficients. Negele and Vautherin (1973), Nucl. Phys. A **207**, 298, pp. 304–306 and Tables 1 and 3, supplies the low-density bulk functional, the adoption of the 1972 gradient parameters, and the 11 cells. The generated quantities use fm, MeV and the project's rest masses. Both sources and the symmetric-matter coefficient checks are identifiable from the repository's original source inputs and source audit. [1972 paper](https://doi.org/10.1103/PhysRevC.5.1472) · [1973 paper](https://doi.org/10.1016/0375-9474(73)90349-7)

The printed model has symmetric D03 = D30 = −6.62699×10³. The alternative has D03 = D30 = −6.62699×10⁵. Only this symmetric pair changes. The second assignment recovers the printed Table II Set 1 spin-independent Skyrme coefficients near saturation. Its status is an internally cross-checked conditional reconstruction. No authenticated published correction was assumed.

## Restricted semiclassical profiles

At a Table 3 anchor with cell volume V₀=(N+Z)/n₀ and printed exterior gas density n_g, take N_b=N−n_gV₀ as the neutron excess. The existing generator's Fermi seed sets R_c=[3(N_b+Z)/(4π×0.14 fm⁻³)]^(1/3) and diffuseness 0.65 fm. Write f(r₀)=1/[1+exp((r₀−R_c)/0.65 fm)] divided by its volume average. The trial species shapes at n₀ are ρ*n=(N−N_b+N_b f)/V₀ and ρ*p=Z f/V₀. These are modeling choices. Both species have independent contrasts λ in [0,1], interpolating convexly between their mean and trial shape. Compression keeps r₀/R₀=r/R and scales the densities by n/n₀. The 49-node matching input uses the same anchor N=950, Z=32 throughout and covers 0.0475–0.0789 fm⁻³. At its lower end this cell differs from the Table 3 equilibrium composition at that density.

At every density, the code compares the grid minima on a 9×9 contrast mesh, refines them by a two-dimensional stencil over 17 levels and also follows the previous density's solution. It conserves both particle numbers on the radial mesh. Each species has kinetic density τ=(3/5)(3π²)^(2/3)ρ^(5/3)+|∇ρ|²/(36ρ), with zero spin current. This chosen kinetic closure omits shell structure, pairing, spin-current relaxation and other semiclassical gradient terms. Energies, normalization, neutral-cell Coulomb direct energy, proton Slater exchange, relativistic electrons and variational dilation pressure use the existing DME kernel. Pressure satisfies P=n²∂(E/A)/∂n at fixed cell and dimensionless trial profile. The density compression probe compares this derivative to the analytic virial at every accepted point.

For the charge derivative, extend the two optimized species shapes to real Z at fixed A,V by scaling their amplitudes to N=A−Z and Z. The closure defines μ_L=∂E_cell/∂Z, μ_e,eff=μ_e,free−average(V_C), μ_Q=μ_L−μ_e,eff, and μ_B=(ε+P)/n−Y_e μ_L. Finite differences with δZ=0.001 verify the calculated charge derivative. Both μ_L and the native Q5 are nonzero. A projected barotrope uses the enthalpy per baryon (ε+P)/n; matching its scalar thermodynamics establishes no beta-equilibrium condition.

CompOSE particle ID 0 records the known electron fraction 32/982. The cell includes a neutron gas, and its A and Z do not directly give the nuclear-group abundance. The sidecar therefore omits the nuclear tuple. Numerical sector ID 1 marks resolved density contrasts and ID 2 marks the homogeneous limit of this ansatz. The two sampled intervals are recorded separately, with no inferred phase boundary inside their gap. The separate 11-row anchor files use each source cell once. Pressures are negative at one printed-mode anchor and two alternative-mode anchors. These data remain in the diagnostic CSV and are excluded from matching.

## Beta-relaxed bulk control

Uniform matter has constant densities, zero gradients and spin density, and τ=τ_TF. The B-dependent part proportional to τ−τ_TF vanishes. The printed and revised coefficient models then reduce to exactly the same homogeneous functional. At 257 log-spaced densities within the printed interval 0.000279–0.0789 fm⁻³, the native code solves μ_p+μ_e−μ_n=0 on the neutron-rich branch, checks local charge curvature and probes 257 trial fractions spanning 10⁻¹¹ through 1−10⁻¹⁰ for competing states. The energy includes the Table 1 bulk polynomial, nonrelativistic free Fermi gases with a common neutron kinetic mass, separate rest masses, proton exchange and relativistic electrons at n_e=n_p. No muons enter this interval in the computed model. The pressure comes from the analytic density derivative. Reoptimized finite differences test this pressure and 0<dP/dε<1. The electron, neutron, proton and explicitly vanishing muon fractions are in eos.compo. The zero entry on the single Yq axis denotes the beta-equilibrated path and the composition sidecar gives its density-dependent Y_e.

This full-domain beta solution is homogeneous throughout. Within the source paper, the tabulated clustered ground states generally have lower energies than the homogeneous matter. The source's 11 uniform-energy values were already compared with the same bulk kernel in the earlier 1973-source report, with a maximum difference of 0.01682 MeV per baryon. The derived bulk table and its strict matching campaign provide a separate diagnostic of this approximation.

## Reproduction

From the repository root, supply a Fortran 2008 compiler and linked LAPACK and BLAS. The second run changes the radial mesh while using identical density nodes. `LAPACK_LIBS` can select versioned libraries on installations without development symlinks.

```sh
bash fortran_eos/run-nv-profile.sh /tmp/nv128 128 49
bash fortran_eos/run-nv-profile.sh /tmp/nv256 256 49
python scripts/verify_nv_reconstructed_profiles.py /tmp/nv128 /tmp/nv256
```

The verification record is `validation.json`. The generator checks a native CompOSE readback. The independent script checks energy and pressure units, generalized Euler balance, beta flags, shape sectors, particle accounting, finite differences and radial convergence. The archived matching inputs were frozen before the final readback and dilation instrumentation. A subsequent run of the checked source reproduced every native file either byte for byte or, for the two restricted eos.thermo files, with at most 8.9×10⁻¹⁶ absolute change in a stored dimensionless or MeV-per-baryon field.
