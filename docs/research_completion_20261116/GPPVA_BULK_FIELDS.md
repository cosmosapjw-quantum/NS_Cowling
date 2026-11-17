# GPPVA source bulk and field functional

This node implements the nucleonic bulk inputs for DD2, DDME2, TW, FSU2, FSU2R, FSU2H, NL3wrL55 and TM1e. It supplies scalar/vector/isovector potentials, density-dependent local vertices and derivatives, uniform stationary states and composition to the spatial Thomas–Fermi implementation. FSU2H's hyperonic-core applications do not introduce hyperons into its low-density nucleonic crust.

## Source and normalization

The source input ledger is the prior `data/research_generator_20261116/crust_literature.json`, whose SHA-256 is recorded with the executed campaign. Fresh primary-source checks confirmed the [FSU2/R/H numeric table](https://iopscience.iop.org/0004-637X/863/1/104/suppdata/apjaad049t1_ascii.txt), the [revised FSU2R/H definitions](https://arxiv.org/abs/1708.08681), and [TM1e's J=31.38 MeV, L=40 MeV specification](https://arxiv.org/abs/2001.10143). DD parameters and split/common nucleon masses are also independently cross-checked against the existing native backend; their original references are retained in the executable `SOURCES` map.

Use raw meson fields phi, w, r and particle/meson masses in MeV, spatial distances in fm, and densities in fm^-3. The effective masses are M_i*=M_i-g_sigma phi. The isovector vertex is g_rho*tau3/2, with n3=np-nn. Native DD couplings use the no-half vertex, so their stored g_rho values are doubled in this module without changing physics. All three DD models use the local vector baryon density in their rational/exponential vertices. DDME2 and TW have common neutron/proton mass 939 MeV. DD2 retains the source split values 939.56536 and 938.27203 MeV.

Define potentials in MeV^4:

U(phi)=m_sigma^2 phi^2/2 + a3 phi^3/3 + a4 phi^4/4,
W(w,r)=m_omega^2 w^2/2 + m_rho^2 r^2/2 + b4 w^4/4 + q w^2 r^2.

For FSU2/R/H, a3=kappa*g_sigma^3/2, a4=lambda*g_sigma^4/6, b4=zeta*g_omega^4/6, and q=Lambda*g_omega^2*g_rho^2. The revised R/H Lambda is 0.045. The tables list g_i squared, which is converted once to the vertex.

For NL3wrL55, the source scalar coefficients multiply the unscaled field: a3=(4.384*939 MeV)/2 and a4=-173.31/6. Multiplying these numbers again by powers of g_sigma would be a different functional. For TM1e, phi=-sigma_source maps M*=M+g_sigma*sigma_source to the convention here: a3=-g2*hbar_c, a4=g3, b4=c3. Thus the negative source g2 produces a positive cubic coefficient. The common TM1e nucleon mass is 938 MeV.

## Variational equations and thermodynamics

At fixed local particle densities, the off-shell energy density is

E = E_kin + [U-W]/(hbar_c)^3 + g_omega*w*nB + g_rho*r*n3/2
    + [grad(phi)^2-grad(w)^2-grad(r)^2]/(2*hbar_c).

The fields satisfy

[-(hbar_c)^2 Laplacian]phi + U_phi = (hbar_c)^3 g_sigma rho_s,
[-(hbar_c)^2 Laplacian]w + W_w = (hbar_c)^3 g_omega nB,
[-(hbar_c)^2 Laplacian]r + W_r = (hbar_c)^3 g_rho n3/2.

The vector fields are stationary saddle directions of this off-shell functional; minimizing all raw fields as though their gradient terms were positive would change the equations. Eliminating stationary vector fields produces the familiar positive vector energy. In homogeneous matter this is m_omega^2*w^2/2 + m_rho^2*r^2/2 +3*b4*w^4/4 +3*q*w^2*r^2, divided by (hbar_c)^3.

The local density-dependent rearrangement is

Sigma_R = g_omega'(nB)*w*nB + g_rho'(nB)*r*n3/2 - g_sigma'(nB)*phi*rho_s.

It enters both nucleon chemical potentials equally:
mu_n=EFn+g_omega*w-g_rho*r/2+Sigma_R,
mu_p=EFp+g_omega*w+g_rho*r/2+Sigma_R.

At stationary uniform fields, P=sum(P_kin)+(W-U)/(hbar_c)^3+nB*Sigma_R. Independently differentiating energy with respect to nn and np checks these expressions. The tests also differentiate U and W directly and compare the symmetry derivative S=(1/4)d(mu_n-mu_p)/d(delta) against the analytic rho-propagator formula for common-mass models.

## Physical scalar branch

NL3 and several FSU scalar quartics are negative, so the fitted polynomial cannot be globally minimized over arbitrary raw scalar fields. The source mean-field branch is evaluated in the positive-effective-mass sector. The uniform solver brackets scalar sign changes on a declared 65-point scan, requires positive local scalar stiffness, and ranks retained stationary roots within that sector. The executed crust-density scan has exactly one retained stable root at every point and continues from phi=0 in vacuum. Its maximum scalar field is 48.62634 MeV, well inside the positive-mass interval. This is the physical nuclear branch tested here; an unbounded-field global minimum is not an acceptance criterion for these fitted models.

At very low density, beta equilibrium can sit on an absent-species boundary. An initial interior-only root bracket failed for FSU2 and DD2. The repair minimizes on 0<=Yp<=1 and enforces the corresponding one-sided weak-equilibrium inequality. At the proton-only boundary, the baryon multiplier is mu_p+mu_e; mu_n is the neutron insertion threshold. Normal crust-density solutions are interior and agree with the native DD solver.

## Executed reference campaign

The campaign contains 2,400 uniform field states: 60 densities from vacuum through 0.2 fm^-3 at Yp=0,.1,.3,.5,1 for all eight models. It also computes 80 npemu beta states and 12 comparisons to the existing independent native DD backend. Maximum field-equation residual is 2.98e-14 fm^-3, uniform Euler residual 5.69e-14 MeV fm^-3, and absolute native-comparison difference in epsilon 2.16e-12 MeV fm^-3, pressure 9.95e-13 MeV fm^-3, mu_n 1.90e-11 MeV, mu_e 1.66e-10 MeV, and Yp 5.24e-14. Every sampled state has one stable scalar root.

| Model | n_sat / fm^-3 | E/A / MeV | K / MeV | J / MeV | L / MeV |
| --- | ---: | ---: | ---: | ---: | ---: |
| DD2 | .1490763 | -16.02258 | 242.7232 | 31.67001 | 55.03354 |
| DDME2 | .1519472 | -16.14174 | 250.9021 | 32.30346 | 51.25749 |
| TW | .1530121 | -16.24715 | 240.2638 | 32.76652 | 55.31017 |
| FSU2 | .1504019 | -16.26488 | 237.6852 | 37.57314 | 112.68003 |
| FSU2R | .1504024 | -16.26502 | 237.6894 | 30.68461 | 46.85267 |
| FSU2H | .1504026 | -16.26483 | 237.6823 | 30.54120 | 44.46043 |
| NL3wrL55 | .1482075 | -16.24076 | 271.5544 | 31.68596 | 55.30209 |
| TM1e | .1452229 | -16.26312 | 281.1645 | 31.38192 | 40.00475 |

The published FSU table gives common isoscalar values (.1505,-16.28,238) and J/L (37.6,112.8),(30.7,46.9),(30.5,44.5). Our largest corresponding offsets are 0.0000982 fm^-3, 0.0152 MeV, 0.318 MeV, 0.0413 MeV and 0.120 MeV. NL3wrL55 published J/L 31.66/55.23 differ by 0.0260/0.0721 MeV. These small differences are reported rather than removed by fitting parameters to a table; the exact numerical conventions behind all offsets have not been recovered. The source scalar/rho normalization is independently supported by the field variations, symmetry identity and source-specific values. CompOSE curve equality is a separate comparison, not the pass criterion here.

The 45-test suite passed. Reproduction from repo root:

```sh
PYTHONPATH=python_frontend COMPOSE_EOS_FORTRAN_ROOT=. \
COMPOSE_EOS_FORTRAN_EXE=build/compose_eos_fortran/compose_eos_fortran \
python3 scripts/verify_gppva_bulk.py build/gppva_bulk_replay
PYTHONPATH=python_frontend python3 -m pytest -q python_frontend/tests/test_gppva_rmf.py
```

The fresh output directory contains SUMMARY.json and compressed full state arrays. Retained copies, command/environment details, source hashes and test logs are under data/research_completion_20261116/gppva_bulk/. Finite-range spatial fields, cell-size selection and pasta competition are handled by the separate TF node.

The independent source-normalization and all-eight-model spatial functional review is documented in [TF_THEORY.md](TF_THEORY.md), with its numerical receipt in `data/research_completion_20261116/tf_review/functional_review.json`. A directly executable counterexample for the rejected interior-only beta bracket is retained as `gppva_bulk/reproduce_beta_boundary.py` with its JSON output. The remaining DDME2 versus Grill 2014 Table 4 energy offset is a separate source-comparison issue: the published table directly reports energy in fm^-4, so a CompOSE Q7/header mass-reference shift does not explain it. No energy shift or bulk refit was applied.

The screened-electron envelope reduction and its independent red/green, Poisson, charge, gradient and repaired-profile checks are detailed in `data/research_completion_20261116/gppva_bulk/SCREENED_REVIEW.md`. Production physics was not modified by this reviewer.
