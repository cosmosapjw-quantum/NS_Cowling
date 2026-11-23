# NV original-paper generation: source, equations, implementation, evidence

Date: 2026-11-23. Baseline: 0fbeeda62f078209cad03798c435933c34e2e1a1.
Branch: research/eos-original-papers-20261123.

## Outcome and claim boundary

The supplied original Negele--Vautherin paper closes the provenance of the published ground-state energy fit. The native `NV1973_author_fit` generation path uses its Table 4 coefficients and preserves its Table 3 data. It generates the requested number of EOS rows on the source density support, including pressure, chemical potential and sound-speed derivatives. The fit is mathematically and numerically the same EOS as the already implemented `NV1974_author_fit`; it must not be counted as an additional distinct surviving EOS.

The original paper also supplies a usable nuclear bulk potential in Table 1. This work implements that potential, its gradient and Hessian, and a separately labeled uniform beta-equilibrium reconstruction. That calculation produces new native numerical diagnostics and agrees with the original eleven uniform-gas energies to a maximum absolute discrepancy of 0.0168109599067829 MeV per baryon. This discrepancy is reported rather than hidden or tuned away.

The full spatially inhomogeneous Hartree--Fock calculation has not been rerun. Its B_IJ, C_IJ, D_IJ coefficient tables and numerical spin-orbit constant W are explicitly referred to another original paper; they are absent from the supplied 1973 paper. Access checks found only metadata, not that paper's full text. The implementation therefore preserves the original fit, the uniform reconstruction, and a complete inhomogeneous HF recomputation as three different claims.

## Source authority and immutable identities

Primary paper: J. W. Negele and D. Vautherin, *Neutron star matter at sub-nuclear densities*, Nuclear Physics A207 (1973), 298--320.
DOI: https://doi.org/10.1016/0375-9474(73)90349-7
User-supplied PDF SHA256:
`5611ae66700d983e7f66881839802c2204e16f91cafbf19c670fbd178042bb5b`.

The complete 23-page source text was read. Decisive equations/tables were inspected visually on journal pages 304, 305, 306, 311 and 317. The copyrighted PDF is not added to Git; only bibliographic identity, factual coefficient/table transcription, native code, and evidence are added.

| Object | SHA256 | Status |
|---|---|---|
| New 1973 source facts | 9370f187fc60be9c6e63ce4aba23bc67d81ce703bb11851bed94e7b062968829 | Added |
| New native module | 81e0fcd0cf7976304a275911d70ded17a613b68cc403d3fd33a973f26b1856bf | Added |
| New regression tests | 0dacb4aa8f4c70d714e54624b032f4ba5b5f0e85ab660260915c4c3ea7fb504f | Added |
| Existing NV1974 source | 284ebc3bfaad9f87677df87eb4b6234c5ea74340189e72ee6b4926fce593e408 | Unchanged |
| Existing NV primary module | 672c4bd5f3aa47a22f8faca1446b5861f48ae3ac1c81896a9e5c0c5ef04c096d | Unchanged |
| Historical NUDA 1973 table | 0a00852b60a02f4b6dfdc8f6305c343ca188c41f9cf0f2eb31983294520efe1e | Unchanged |

The source record fixes document hash, mass convention, density support, both Table 4 fits, all eleven Table 3 rows, and the 28 nonredundant Table 1 coefficients. Parser rejection precedes diagnostic creation. It rejects changed source facts, extra/missing rows or columns, NaN/Infinity, oversized exponents, malformed number tokens and trailing noncomment data.

## Original 1973 versus 1974 distinction

At n=0.00577 fm^-3, the original 1973 Table 3 prints 1.996 MeV; the 1974 author account prints 1.966 MeV. Both values are preserved and shown side by side. The existing historical 1973 table is not erroneous on this point.

The ground-state coefficients in 1973 Table 4, p.317, equal the previously verified 1974 coefficients. The new path directly reuses the existing `nv_primary_state` arithmetic after independently validating the 1973 source contract. No interpolation of eleven sparse energy samples is passed off as an original pressure table. Numerical identity was tested at 4097 points including both endpoints, with exact float64 equality of energy offset, energy density, pressure, chemical potential, sound speed and adiabatic index.

The new receipt sets:
- `semantic_equivalent_model: NV1974_author_fit`;
- `distinct_new_physical_eos: false`;
- `full_hartree_fock_recomputed: false`.

## Published fit: computational contract

Let n be baryon number density in fm^-3, x=ln[n/(10^-4 fm^-3)], Q(x)=sum(k=1..7)c_k*x^(k-1), and F=exp[Q] MeV. The published energy offset is e=c_0+F, with:

```
c_ground = [-4, 0.28822899, 0.59150523, 0.090185940,
            -0.11025614, 0.029377479, -0.0032618465, 0.00013543555]
c_uniform = [-4, 1.4821424, -0.040373482, 0.060455728,
             -0.015307639, 0.0034774416, -0.00043627154, 0.000023383473]
```

With D=F Q' and G=F[(Q')^2+Q''], direct zero-temperature thermodynamics gives:
- epsilon = n(m_n c^2+e);
- P = n D;
- mu_b = m_n c^2+e+D;
- dP/dn = D+G;
- c_s^2/c^2 = (D+G)/mu_b;
- gamma = (D+G)/D.

Dimensions: epsilon and P are MeV fm^-3; e, D, G and mu_b are MeV; c_s^2/c^2 is dimensionless. The source fit includes electrons; no electron energy is added again to the exported crust fit.

Support is exactly 2.79e-4 <= n <= 7.89e-2 fm^-3. No extrapolation is made. A logarithmic grid contains both endpoints exactly. The public API accepts 3..100000 rows. The paper writes nucleon rest masses symbolically; the existing project conventions m_n c^2=939.56542052 MeV and m_p c^2=938.27208816 MeV are recorded rather than attributed numerically to the 1973 paper.

The source-author fit retains its source limitations. For 257 rows the minimum mu_b-m_n c^2 is -0.216823102598255 MeV, and seven sampled rows lie below the free-neutron threshold. The maximum energy fit discrepancy from original Table 3 is 0.115411407322675 MeV. These properties are not repaired by changing the source fit and are not a matching acceptance claim.

## Printed Table 1 bulk functional

The source p.305 prints the symmetric polynomial V^NM(n_n,n_p)=sum(i,j=0..6)V_ij n_n^i n_p^j, with redundant off-diagonal values omitted. Each off-diagonal entry is restored symmetrically, not halved. There are 28 printed lower-triangle coefficients. The preceding description of an older finite-nucleus fit uses a polynomial upper index 7; no missing seventh-index coefficients are invented for the printed Table 1 reconstruction.

The native routine uses nested differentiated Horner evaluation to obtain:
V, gradient=(partial V/partial n_n, partial V/partial n_p), and its symmetric 2x2 Hessian.
V is treated as the uniform nuclear potential energy density in the uniform limit of the p.304 functional. V has units MeV fm^-3, its gradient MeV, and Hessian MeV fm^3; V_ij consequently has units MeV fm^[3(i+j-1)].

The 0<=n_n,n_p<=0.15 fm^-3 rectangle is a conservative diagnostic guard, informed by the original 16x16 grid discussion, not a newly demonstrated physical validity domain of the refitted functional. The source describes shifting and weighting the low-density fit; it does not state a uniform rigorous error bound on this rectangle. The EOS provider still uses the much smaller explicitly listed ground-state density support.

Vacuum has V=0 and both first derivatives zero, as enforced by the printed zero coefficients. Exchange symmetry and analytic derivatives are tested against a direct-power real128 Fortran oracle, independently of the Horner implementation.

## Native uniform matter reconstruction

This diagnostic is the constructive partial realization of the original functional that can be implemented without the missing finite-gradient coefficients. It is not an inhomogeneous crust provider.

At fixed n=n_n+n_p, impose n_e=n_p. The source electron approximation is a uniform relativistic free gas; screening and electron Coulomb exchange are neglected, and proton Coulomb exchange uses the Slater term. Pairing is omitted, as in the source. The common nonrelativistic kinetic mass is explicitly chosen as M=m_n; species rest masses remain distinct.

Define K=(hbar c)^2(3pi^2)^(2/3)/(2m_n c^2),
T_q=(3/5)K*n_q^(5/3), and
E_x=-(3/4)(3/pi)^(1/3)alpha*hbar*c*n_p^(4/3).
Then:
epsilon=m_n c^2 n_n+m_p c^2 n_p+T_n+T_p+V+E_x+epsilon_e;
mu_n=m_n c^2+K n_n^(2/3)+V_n;
mu_p=m_p c^2+K n_p^(2/3)+V_p+(4/3)E_x/n_p;
P=(2/3)(T_n+T_p)+n_n V_n+n_p V_p-V+E_x/3+P_e.

The electron energy and pressure use their exact zero-temperature relativistic Fermi integrals, including electron rest energy. A small-momentum series avoids cancellation when hbar*k_F/(m_e*c)<0.01. The uniform energy offset is formed without subtracting two approximately 939 MeV numbers:
e=[(m_p-m_n)c^2 n_p+T_n+T_p+V+E_x+epsilon_e]/n.

Beta equilibrium solves mu_p+mu_e-mu_n=0 by bounded bisection on 0<n_p<n/2. The final stationary point must have a positive derivative of that residual with respect to n_p at fixed n, i.e. a local minimum with respect to composition. It does not certify a global inhomogeneous minimum or absence of other phases. At n_p=0 the curvature is undefined; the type explicitly exposes `beta_curvature_defined=false` rather than interpreting its unused numeric placeholder as a curvature.

Tiny-composition derivatives use powers differentiated before division, avoiding underflow of n_p^2. Electron curvature is calculated from the Fermi momentum rather than subtracting nearly equal mu_e^2 and m_e^2. Native tests include n_p=0, 1e-300, 1e-30 and 1e-12 fm^-3 under invalid/zero/overflow traps.

Published 11-point comparison:
- maximum uniform Table 4 fit discrepancy: 0.0561928763365067 MeV;
- maximum Table 1 reconstruction discrepancy: 0.0168109599067829 MeV.

These are observations, not calibrated acceptance tolerances. Printed coefficient precision, the common kinetic mass choice, other historical constants, and the old index-7 notation are recorded uncertainties; none has been proved to explain the residual discrepancy.

## Implemented interfaces and evidence

Module: `fortran_eos/src/compose_eos_nv_original.f90`.
Public entry:
`generate_nv_original_native(source_path,n_points,eos,output_dir,status)`.
Successful status is lowercase `ok`; failures return precise status strings and never STOP inside the library.

Additional public routines:
`load_nv_original_source`, `nv1973_bulk_potential`,
`nv1973_uniform_composition`, `nv1973_uniform_beta_state`.

Mandatory output diagnostics:
- `nv1973_states.csv`;
- `nv1973_table3_comparison.csv`;
- `nv1973_uniform_reconstruction.csv`;
- `NV_ORIGINAL_GENERATION.json`.

The CompOSE carrier has no composition file because the published smooth crust fit supplies no defensible smooth composition function. The separately reconstructed uniform composition is present only in its explicitly labeled diagnostic CSV.

Validation command:
`python3 -m pytest python_frontend/tests/test_nv_original_native.py -q`

Result: **43 passed in 2.63 s**, with gfortran Fortran-2018, -O2, -fcheck=all, and -ffpe-trap=invalid,zero,overflow. Python compiles, invokes and reads evidence; the independent numerical oracles themselves are Fortran real128. Requested grids 3, 257 and 513 were generated and read through CompOSE I/O.

Native numerical audit measurements:

| Quantity | Maximum error or minimum bound |
|---|---:|
| Bulk value error / (1+abs(reference)) | 2.868035700973644e-15 |
| Bulk gradient scaled error | 2.763003918874832e-14 |
| Bulk Hessian scaled error | 2.485303125174702e-12 |
| Exchange symmetry scaled error | 2.032967587675036e-12 |
| Uniform beta residual, MeV | 1.932676241267473e-12 |
| Uniform Euler residual, MeV fm^-3 | 2.270406085358445e-14 |
| Uniform chemical potential vs real128 five-point derivative, MeV | 1.958433415438776e-13 |
| Minimum sampled composition curvature, MeV fm^3 | 1.315017242215072e4 |

The polynomial oracle covers a 16x16 grid including zero-density endpoints. The thermodynamic oracle covers 129 beta-equilibrium densities across the entire fit support. Direct five-point differentiation of an independently assembled real128 energy evaluates both species derivatives, including electron and exchange terms. Exact published-fit semantic equality is tested at 4097 densities.

Raw measurements and a generated 257-row carrier are under `verification/`. No full-pair matching count is asserted by this module; the parent campaign performs and records those calculations.

## Full HF dependencies and exact next source

The 1973 p.306 explicitly says the Table 1 V coefficients were used with the original B_IJ, C_IJ and D_IJ from reference 21. Page 308 introduces proton spin-orbit energy (W/2) grad(n_n+2n_p) dot J_p but supplies no numerical W. The finite-gradient DME functions/effective mass and spin-orbit coefficient are therefore not inferable from Table 1 or eleven EOS energy samples.

Required next original source:
J. W. Negele and D. Vautherin, *Density-Matrix Expansion for an Effective Nuclear Hamiltonian*, Physical Review C **5**, 1472--1493 (1972).
DOI: https://doi.org/10.1103/PhysRevC.5.1472

Access investigated on 2026-11-23:
- APS original page: authorization required;
- HAL https://hal.science/in2p3-00010577v1 : record explicitly says **No file**, 0 downloads;
- INSPIRE https://inspirehep.net/literature/79401 : metadata and external identifiers, no full-text document;
- OSTI https://www.osti.gov/biblio/4644874 : bibliographic record;
- SciSpace paper record: open_access=false; full text and requested coefficients not returned;
- Exa/Parallel primary-source discovery: records above, no accessible coefficient-bearing original PDF.

This is a specific current data dependency, not a claim that original HF reconstruction is impossible. The minimum next action is to supply or retrieve the 1972 original and inspect its coefficient tables and W. If it in turn requires an external effective-interaction table, that dependency must be identified from the actual paper before further requests.

Source boundary conditions already recovered for subsequent HF work: even-parity neutron full wavefunctions vanish at the cell radius; odd-parity full radial wavefunctions have zero radial derivative. For a reduced radial function u=rR, the latter means u'=u/R, not u'=0. The source resets density and kinetic density in the last 1/k_F layer to the uniform-gas average. Neutron spin-orbit and pairing are omitted; proton spin-orbit is retained. Global configuration searches require independent initial configurations and variation of N,Z and cell radius. These source requirements are documented; they have not been silently replaced by a simpler Skyrme functional.

## Research-loop disposition

The cycle produced source facts, explicit thermodynamic derivation, native evaluation, a nontrivial uniform functional reconstruction, independent numerical validation, hostile-input fixes, and a precise next-source request. The source-supported fit path is ready for the common generator and all-pair matching campaign. Remaining HF coefficient access is isolated from the successfully executed work. No strong nonexistence or infeasibility claim is made.
