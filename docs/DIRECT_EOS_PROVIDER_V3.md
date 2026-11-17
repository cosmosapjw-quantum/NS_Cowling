# Direct source-knot EOS provider for relativistic Cowling modes (v3)

## Scope and decision

Version 3 changes the default EOS provider, not the TOV equations or the relativistic Cowling perturbation equations. The objective is to consume the already-generated, cold, beta-equilibrated EOS table without projecting its pressure or energy-density knots onto a second EOS table. The original normalized CSVs remain immutable inputs. The optional q-linear (P(\mu_B)) adapter from v2 is retained only as a legacy comparison method.

The adopted representation is a monotone PCHIP of (g(P)=\epsilon(P)-P) against (u=\ln(P/P_*)), followed by first-law reconstruction of (n_B) and (\mu_B). Here (P_*>0) is an arbitrary fixed reference pressure (the implementation uses (P_*=1\,\mathrm{Pa})); only differences in (u) enter derivatives. The method preserves every source (P_i,\epsilon_i) knot. It does not claim to preserve the source (n_{B,i},\mu_{B,i}) simultaneously at every knot.

## What the upstream EOS generator already guarantees

The bundled compose_eos 3.3.0 generator uses a master free energy per baryon (F(n_B)), with (x=\ln n_B) and (Q_1=dF/dx):
[
 P=n_B Q_1,\qquad \mu_B=F+Q_1,\qquad \epsilon=n_B F.
]
At the saved knots this gives
[
 \epsilon_i+P_i=\mu_{B,i}n_{B,i},\qquad
 d\epsilon=\mu_B\,dn_B,\qquad dP=n_B\,d\mu_B.
]
This is already thermodynamically consistent source data. A table, however, stores finitely many samples and does not uniquely determine a continuous interpolation between them. The v2 Fortran provider imposed a cubic-Hermite (P(\mu_B)) cell with (P'=n_B); it then derived (\epsilon=\mu_BP'-P) and checked positivity, compressibility and causality throughout each cell. The failed admission was a property of that specific interpolant, not evidence that the generator's source knots violated the first law.

## Constructive derivation of the v3 representation

Define
[
 u=\ln(P/P_*),\qquad g(u)=\epsilon(P)-P.
]
Interpolate (g(u)) with a shape-preserving, monotone piecewise cubic Hermite interpolant through the original values (g_i=\epsilon_i-P_i). Then set
[
 \epsilon(P)=P+g(u).
]
Differentiating with respect to (P),
[
 \frac{d\epsilon}{dP}=1+\frac{1}{P}\frac{dg}{du},
 \qquad
 \chi\equiv c_s^2=\frac{dP}{d\epsilon}
 =\frac{P}{P+g_u}.
]
Thus (P>0) and (g_u\ge0) imply (0<\chi\le1) throughout every cell. This is the continuous mechanical-stability and causality gate. The implementation additionally checks the analytic minimum of the PCHIP derivative polynomial on each cell, not just sampled query points.

For a cold barotrope the first law and Euler identity are
[
 d\epsilon=\mu_B\,dn_B,\qquad \epsilon+P=\mu_Bn_B.
]
Eliminating (\mu_B) and using (u=\ln P) gives
[
 \frac{d\ln n_B}{du}
 =\frac{P+g_u}{2P+g}.
]
The code normalizes (n_B(P_{\min})) to the first source density and integrates this equation cell by cell. It then computes
[
 \mu_B(P)=\frac{\epsilon(P)+P}{n_B(P)}.
]
The integral is evaluated with a 16-point Gauss-Legendre rule per native log-pressure cell; the same quadrature rule is used for partial-cell queries. The Fortran EOS callback returns (P,\epsilon,n_B,\mu_B,\chi) as one thermodynamically coupled state.

## Why this is the least invasive admissible change

The stellar structure and Cowling kernels already request the EOS as a function of pressure and consume the thermodynamic state returned by a callback. Therefore the background TOV equations, the Cowling equations, the central regularity prescription, the surface residual, and the finite-pressure surface treatment do not need a new physical derivation. The required change is localized to the EOS provider and to the optional cell-by-cell shooting coordinate, which now uses (t=(\ln P-\ln P_i)/(\ln P_{i+1}-\ln P_i)).

The irreducible change is that (n_B(P)) and (\mu_B(P)) are reconstructed from the first law, rather than independently interpolated through all their source knots. This is necessary because an arbitrary interpolation of four tabulated quantities generally violates their differential identities between knots. Preserving all four sets of knots while imposing a monotone causal continuous barotrope is not guaranteed by the input data. The v3 choice prioritizes exact preservation of the two variables that define the pressure-energy barotrope used by TOV, then reconstructs the remaining state variables from the cold first law. It is a declared interpolation closure, not a modification of the microscopic EOS model.

## Alternative representations tested and not selected

1. **Cubic Hermite in free energy (F(\ln n_B)).** This passes the original source knots exactly and derives all thermodynamic quantities from one function. However, the straightforward endpoint-slope cubic produced (F_x+F_{xx}\le0) in 24,374 of 109,776 cells. It therefore cannot be used as a universal stable/causal interpolant for this cohort without a constrained higher-order construction.
2. **Cubic Hermite in (\epsilon(n_B)).** The straightforward endpoint-slope construction produced nonpositive (n_B\epsilon''(n_B)) in 24,465 of 109,776 cells. It was rejected for the same cell-interior reason.
3. **Monotone (q(\mu_B)=n_B/\mu_B) with exact pressure-area constraints.** Requiring (q\ge0), (q'\ge0), and exact preservation of endpoint (P,n_B,\mu_B) imposes an integral feasibility condition on every cell. At least 14,763 of 109,776 cells fell outside the allowed endpoint-area range. Thus it cannot preserve every source knot while imposing these constraints on the full cohort.
4. **v2 q-linear adapter.** The source (q_i) sequences were already nondecreasing and weighted PAV pooled zero knots. The nontrivial v2 change was the between-knot linear-(q) closure, followed by pressure reintegration and Euler energy reconstruction. It retained (n_B) knots to roundoff but changed (P) and (\epsilon) knots. It remains available for reproducibility/comparison, but is no longer the default.

## Source admission and physical limits

The v3 provider first requires the source knot identity
[
 \max_i\frac{|\epsilon_i+P_i-\mu_{B,i}n_{B,i}|}
 {\max(|\epsilon_i+P_i|,\mathrm{tiny})}\le2\times10^{-12}.
]
It separately requires strictly increasing source (P_i,\mu_{B,i}), positive (g_i), and strictly increasing (g_i). Only then does it build and audit the monotone PCHIP cells. This distinction prevents an inconsistent input table from being silently made to look thermodynamically consistent by reconstructing (n_B) and (\mu_B). The deliberately inconsistent fixture `compose_134_si.csv` has a raw closure residual (1.1327\times10^{-4}) and is rejected by the new source gate; it is not part of the 234 accepted production EOS cohort.

For the 234 normalized production EOS tables, the raw knot closure maximum is (6.20\times10^{-16}). The continuous PCHIP derivative is nonnegative in all 109,776 cells, with zero rejected cells. The maximum relative difference between first-law-reconstructed and source (n_B) values at the knots is (9.8825\times10^{-4}). This density difference is the explicit cost of preserving (P,\epsilon) exactly while making the continuous closure satisfy the first law. It must be treated as an interpolation systematic, not hidden as numerical roundoff. The source (n_B,\mu_B) columns remain intact in their original CSVs and are available for comparison.

This is a cold, beta-equilibrated, barotropic closure. It does not implement frozen-composition perturbations, superfluid degrees of freedom, phase-conversion interface conditions, rotation, metric perturbations, gravitational-wave damping, or extrapolation to an exact zero-pressure surface. The stellar surface remains the declared finite positive pressure cutoff. It is not a full-GR mode solver.

## Fortran implementation changes

- `ns_table_eos`: replaces the default cubic-Hermite (P(\mu_B)) provider by (g(\ln P)) PCHIP, analytic derivative checks, first-law density quadrature, and a pressure-to-(\mu_B) inverse for compatibility callers.
- `ns_cell_eos` and `ns_cell_shooting`: use the direct pressure provider and log-pressure cell coordinate; the old `mu_bridge_t` helper remains only for other legacy analytic matching tasks.
- `cell_solve`: when the controlled perturbation amplitude is exactly zero, passes the loaded normalized source table through unchanged. It does not recompute (\epsilon=\mu_Bn_B-P) for a null perturbation. Nonzero amplitudes remain explicit, separate sensitivity experiments.
- The TOV/Cowling differential equations and their declared boundary conditions are unchanged.

The pressure-inverse stopping tolerance is relative to the local pressure scale, not `max(P,1)`. The latter is inappropriate for SI pressures much smaller than one pascal and can stop before resolving the requested pressure.

## Validation record for v3

- Python direct-provider tests: 4 passed; full pipeline tests: 8 passed.
- Fortran: 62 source modules compiled and both `case_driver` and `cell_solve` linked with gfortran in Fortran 2018 mode and runtime checks enabled. The repository's FPM harness could not be run because FPM 0.13.0 is not installed; this is an environment limitation, not a passing FPM test.
- Final Fortran table-provider cohort audit with the raw-knot gate enabled: 234/234 tables admitted; 109,776 cells; zero rejected cells; zero source first-law failures. Maximum raw source closure is 6.1990048995e-16, below the 2e-12 gate. The intentionally inconsistent fixture is rejected with a raw closure of 1.1327080889e-4. Full CSV/JSON results are in data/eos_direct_provider_v3/fortran_direct_eos_cohort.csv and .json.
- A direct-provider representative (\ell=2) Cowling spectrum was recomputed at the same target mass for `GPPVA_FSU2H__BBB_BHF2__quintic_P`: (M=2.6446062574938932\times10^{30}\,\mathrm{kg}), (R=11.2249075\,\mathrm{km}), (f=2.517378706\,\mathrm{kHz}), (p_1=7.205657666\,\mathrm{kHz}). The relative target-mass residual is (1.54\times10^{-13}); normalized matching residuals are (4.16\times10^{-9}) and (5.61\times10^{-12}). Eigenfunction node identities are unchanged: (f:(0,0)), (p_1:(1,1)) on the inner/outer matching segments.
- The old 234-row mode cohort remains the v2 q-linear-adapter result set. It has not been silently relabelled as v3.
- **Remediation update (2026-11-14, NSC_EOS390_FULL_REMEDIATION):** the statement above is superseded by measurement. The stored `direct_v3_spectrum_234.csv` carries fully populated `f_v3_hz`/`p1_v3_hz` columns, and an independent recomputation of six rows (the representative above plus a stratified sample across BBB_BHF2/GM_GM1/RG_SLY4 cores and css/cubic/quintic methods) with the native Cowling solver on the immutable normalized CSV inputs reproduces every stored v3 frequency to within 3e-12 relative, with identical eigenfunction node identity (f:(0,0), p1:(1,1)). The 234 v3 columns are therefore real direct-provider computations, not v2 leftovers. Provenance details: `data/remediation/DIRECT_V3_PROVENANCE_MANIFEST.json`. Two classification notes from the same audit: (i) the 39 `maxwell`-labelled rows carry the historical v1 smooth-plateau method id on legacy-matched inputs — they are legacy reference spectra, not sharp-Maxwell constructions; (ii) the 78 RG_SLY4 rows are numerically valid as computed, but RG_SLY4 fails current source admission for smooth matching (nonzero Q5), so those rows are legacy reference, not production-admitted outputs.

## Literature context

The choice follows the general principle in thermodynamically consistent table interpolation that one should interpolate a thermodynamic potential and derive its conjugate quantities consistently, rather than independently interpolating mutually constrained variables. Relevant examples include Swesty (1996), Timmes & Swesty (2000), and Baturin et al. (2019). The present construction is specifically a one-dimensional cold-barotrope closure in (g(\ln P)); it is not the same method as the multi-dimensional Helmholtz free-energy interpolation in those works. Spectral and continuous-sound-speed parameterizations (Lindblom 2010, 2018; Servignat et al. 2024) are useful for flexible EOS representations but solve a different problem when exact source pressure-energy knots must be retained. Mode studies such as Kunjipurayil et al. (2022) and crust-core interpolation analyses motivate treating the interpolation choice as a frequency-systematic source.

## References

- Swesty, F. D. (1996), *Thermodynamically consistent interpolation for equation of state tables*, Journal of Computational Physics **127**, 118–127. DOI: 10.1006/jcph.1996.0162.
- Timmes, F. X. & Swesty, F. D. (2000), *The Accuracy, Consistency, and Speed of an Electron-Positron Equation of State Based on Table Interpolation of the Helmholtz Free Energy*, The Astrophysical Journal Supplement Series **126**, 501–516. DOI: 10.1086/313304.
- Baturin, V. A. et al. (2019), *Interpolation of equation-of-state data*, Astronomy & Astrophysics **626**, A108. DOI: 10.1051/0004-6361/201935669.
- Lindblom, L. (2010), *Spectral representations of neutron-star equations of state*, Physical Review D **82**, 103011. DOI: 10.1103/PhysRevD.82.103011.
- Lindblom, L. (2018), *Causal representations of neutron-star equations of state*, Physical Review D **97**, 123019. DOI: 10.1103/PhysRevD.97.123019.
- Servignat, B. et al. (2024), *One- and two-argument equation of state parametrizations with continuous sound speed for neutron star simulations*, Physical Review D **109**, 103022. DOI: 10.1103/PhysRevD.109.103022.
- Kunjipurayil, A. et al. (2022), *f- and p1-mode oscillations of neutron stars and their dependence on the equation of state*, Physical Review D **106**, 063005. DOI: 10.1103/PhysRevD.106.063005.
