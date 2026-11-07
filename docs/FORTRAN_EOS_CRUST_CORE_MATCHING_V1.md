# Native Fortran smooth crust-core EOS matching v1

**Scope.** Six historical smooth, one-dimensional, cold barotropic matching constructions: maxwell (compatibility token), css, cubic_F, cubic_P, quintic_F and quintic_P.

**Status.** Native code and the Python wrapper are implemented on the feature branch. Matching-specific CI and the final source-package artifact are the release gates recorded separately in the adversarial audit.

**Excluded by design.** Microscopic composition reconstruction, the legacy literature_beta_eq approximation, and genuine first-order Maxwell or Gibbs phase transitions are not implemented by this v1.

## 1. Problem and conventions

The inputs are two CompOSE-derived cold, one-dimensional barotropes, denoted cr (crust) and co (core). The result is a single-valued barotropic table. A selected upper crust interval and lower core interval are replaced by an interpolated bridge; every source row outside the selected interval is retained.

| Quantity | Meaning | Unit |
|---|---|---|
| n_B | baryon number density | fm^-3 |
| P | pressure | MeV fm^-3 |
| epsilon | total energy density | MeV fm^-3 |
| F=epsilon/n_B | energy per baryon, including rest mass | MeV |
| Q1=P/n_B | pressure per baryon | MeV |
| mu_B | baryon chemical potential | MeV |
| c_s^2=dP/d epsilon | squared sound speed in units of c^2 | dimensionless |

The CompOSE thermodynamic columns are interpreted as
\[
 Q_1=P/n_B,\quad F=(Q_6+1)m_n,\quad
 \mu_B=(Q_3+1)m_n,\quad
 \epsilon=n_B(Q_7+1)m_n.
\]
For cold one-component matter,
\[
 \mu_B=\frac{\epsilon+P}{n_B}=F+Q_1,\qquad
 dF=Q_1\,d\ln n_B,\qquad dP=n_B\,d\mu_B.
\]
The construction uses F as its thermodynamic master function. Interpolation does not introduce an additional microscopic degree of freedom.

The repository's EOS tables use m_n and m_p in MeV, density in fm^-3 and pressure/energy density in MeV fm^-3. The interface does not add a new units-conversion layer.

## 2. Input-state contract

The native routine requires finite EOS arrays, at least three positive strictly increasing density points, positive pressure, valid nucleon masses, T=0, beta_equilibrated metadata, zero entropy within tolerance, and agreement of Q6 and Q7. More explicitly:
\[
 |T|\leq10^{-12}\ {\rm MeV},\qquad
 \max_i|Q_{2,i}|\leq10^{-10},\qquad
 \max_i|Q_{6,i}-Q_{7,i}|\leq2\times10^{-6}.
\]
Both input metadata records must also state beta equilibrium.

If either input has an eos.compo file or usable composition, the caller must explicitly pass drop_composition=True. This is not a cosmetic option. The one-dimensional bridge equations determine a barotrope, not the microscopic neutron, proton, electron, muon, nuclear-cluster or phase-index profiles. Native v1 therefore emits a strict thermo-only result with composition_file_present=False, composition_available=False, Nadd=0, zero particle/nuclear arrays and zero phase labels. The output has Q2=Q5=0 and Q6=Q7. The Q4 charge-chemical-potential column is retained on the source portions after mass normalization and smoothly joined through the bridge; it is not the result of a new composition solve.

The Python API refuses composition-bearing inputs unless drop_composition=True. The Fortran routine repeats that check against serialized metadata. This double gate is intentional.

The native implementation does not support core_compo_mode=literature_beta_eq. The historic Python package invokes model-dependent parabolic symmetry-energy fits for selected microscopic EOSs under that option. It is not a general composition reconstruction and is not silently approximated here.

## 3. Transition-window selection

### 3.1 Core starting point

If the core input has usable composition, the search begins at the first homogeneous row with n_B >= 0.03 fm^-3, where homogeneity is identified by Abar=0. Otherwise it begins at the first row with n_B >= 0.045 fm^-3. These are historic screening conventions, not universal crust-core transition densities.

The search scans the upper crust rows in descending density, limited to approximately the last 150 entries, and scans core rows from the selected starting point upward.

### 3.2 Rest-mass gauge alignment

The historic alignment is
\[
 \Delta m_{\rm cr}=m_n^{\rm co}-m_n^{\rm cr},
 \qquad \Delta m_{\rm co}=0.
\]
The crust F and screening chemical potential are shifted by Delta m_cr; the core is unshifted. This aligns rest-mass conventions. It is not a fitted parameter.

For a candidate pair (n_1,P_1,F_1,mu_1^source) and (n_2,P_2,F_2,mu_2^source), define
\[
 \bar n=\frac{P_2-P_1}{\mu_2^{\rm source}-\mu_1^{\rm source}},
 \quad
 \overline Q_1=\frac{F_2-F_1}{\ln(n_2/n_1)},
\]
\[
 \theta_P=\frac{\bar n-n_1}{n_2-n_1},\qquad
 \theta_F=\frac{\overline Q_1-P_1/n_1}{P_2/n_2-P_1/n_1}.
\]
A candidate must meet
\[
 n_2>1.10n_1,\quad P_2>1.10P_1,\quad
 \mu_2^{\rm source}>\mu_1^{\rm source}+0.2\ {\rm MeV},\quad
 F_2>F_1+0.2\ {\rm MeV},
\]
and 0.12 < theta_P, theta_F < 0.88.

These inequalities are a candidate-screening heuristic. They do not prove convexity, stability or causality. Once a pair is selected, the bridge uses the Euler chemical potentials mu_1=F_1+Q_1,1 and mu_2=F_2+Q_1,2, then constructs and validates the complete bridge.

## 4. Endpoint sound speeds

The endpoint values are estimated from one-sided shape-preserving cubic Hermite interpolants of ln P and ln epsilon as functions of x=ln n_B. Six points are used on each side. At the endpoint,
\[
 c_s^2=\frac{P}{\epsilon}
 \frac{d\ln P/dx}{d\ln\epsilon/dx}.
\]
The PCHIP knot derivatives use the Fritsch--Butland weighted harmonic mean internally and the shape-limited one-sided endpoint formula. Non-positive P or epsilon, a non-positive logarithmic energy derivative, or an endpoint sound speed outside (0,1) causes rejection. The code does not clip or repair invalid endpoint derivatives.

This is a numerical derivative of the supplied table, not a microscopic calculation or an uncertainty analysis.

## 5. The six bridge constructions

Let
\[
 x=\ln n_B,\quad \Delta x=\ln(n_2/n_1),\quad
 \Delta\mu=\mu_2-\mu_1,\quad 0\leq u\leq1.
\]
All six constructions are single-valued smooth barotropes with no density jump. They differ in which thermodynamic variable is interpolated and how the first law is enforced.

### 5.1 maxwell: smooth near-pressure plateau compatibility bridge

The historic method name is retained for API compatibility, but this is **not** an exact Maxwell coexistence transition. It uses a cubic Hermite pressure spline P(u), interior knot positions u_a=0.22 and u_b=0.78, and a near-plateau offset
\[
 \delta P_{\rm plat}=0.015(P_2-P_1).
\]
The inner pressures are P_a=P_t-delta P_plat and P_b=P_t+delta P_plat. The plateau reference P_t is determined by a scalar root solve using a 64-point Gauss--Legendre rule on the full interval, matching the historical reference implementation. Each partial free-energy integral F(u) uses a 16-point Gauss--Legendre rule on [0,u]. The root solve enforces
\[
 F_2-F_1=\Delta x\int_0^1
 \frac{P(u)}{\exp(x_1+u\Delta x)}\,du.
\]
Then
\[
 n_B(u)=\exp(x_1+u\Delta x),\quad Q_1=P/n_B,
\]
\[
 F(u)=F_1+\Delta x\int_0^u
 \frac{P(v)}{\exp(x_1+v\Delta x)}\,dv,\quad \mu_B=F+Q_1,
\]
and
\[
 c_s^2(u)=\frac{dP/du}{\Delta x\,n_B(u)\mu_B(u)}.
\]
The metadata records is_exact_maxwell=False, density_jump_fm3=0 and seidov_applicable=False. A physical Maxwell construction would impose coexistence at equal pressure and chemical potential with a density/energy-density discontinuity. This bridge does not represent that physics.

### 5.2 css: two-segment constant-sound-speed bridge

The bridge uses logarithmic interpolation of chemical potential,
\[
 \mu_B(u)=\mu_1(\mu_2/\mu_1)^u.
\]
A scalar parameter s in (0,1) defines
\[
 \mu_m=\mu_1(\mu_2/\mu_1)^{1-s},\qquad
 n_m=n_1+s(n_2-n_1).
\]
The two segment sound speeds are
\[
 c_{s,A}^2=\frac{\ln(\mu_m/\mu_1)}{\ln(n_m/n_1)},\quad
 c_{s,B}^2=\frac{\ln(\mu_2/\mu_m)}{\ln(n_2/n_m)}.
\]
With gamma_A=1/c_s,A^2 and gamma_B=1/c_s,B^2,
\[
 n_B=n_1(\mu_B/\mu_1)^{\gamma_A}
\]
on the first segment and
\[
 n_B=n_m(\mu_B/\mu_m)^{\gamma_B}
\]
on the second. Pressure is integrated analytically from dP=n_B dmu_B. The parameter s is solved so the pressure increment equals P_2-P_1. The method is continuous in P, epsilon and mu_B, but c_s^2 may jump at the internal CSS knot. It is not a first-order phase transition.

### 5.3 cubic_P and quintic_P: pressure--chemical-potential bridges

Here mu_B(u)=mu_1+u Delta mu. The density n_B(u) is a monotone Hermite bridge with endpoint derivatives inferred from the source sound speeds:
\[
 \left.\frac{dn_B}{du}\right|_0=\frac{n_1\Delta\mu}{\mu_1c_{s,1}^2},
 \quad
 \left.\frac{dn_B}{du}\right|_1=\frac{n_2\Delta\mu}{\mu_2c_{s,2}^2}.
\]
Then
\[
 P(u)=P_1+\Delta\mu\int_0^u n_B(v)\,dv,\quad
 Q_1=P/n_B,\quad F=\mu_B-Q_1,
\]
and
\[
 c_s^2(u)=\frac{n_B(u)\Delta\mu}{\mu_B(u)\,dn_B/du}.
\]
This enforces dP=n_B dmu_B by construction.

The cubic method uses a single quadratic density profile when its shape gates pass, and otherwise a one-interior-knot C1 cubic Hermite spline. The quintic method uses a quartic density profile derived from a fifth-order potential when its bounds pass, preserving the endpoint slopes; otherwise it uses an adaptive two-interior-knot cubic Hermite spline for Q1, with harmonic-mean interior slopes and exact endpoint derivatives. The underlying potential is C2; the interpolated derivative profile is C1. The fallback is an explicit numerical branch of the same matching prescription.

### 5.4 cubic_F and quintic_F: free-energy bridges

Here x(u)=x_1+u Delta x and n_B=exp(x). The interpolated variable is Q1(u), constrained by
\[
 F(u)=F_1+\Delta x\int_0^u Q_1(v)\,dv,\qquad
 \int_0^1 Q_1(u)\,du=\frac{F_2-F_1}{\Delta x}.
\]
The derived quantities are
\[
 \mu_B=F+Q_1,\qquad P=n_BQ_1,\qquad
 c_s^2=\frac{Q_1+(dQ_1/du)/\Delta x}{\mu_B}.
\]
The cubic method uses a quadratic Q1 profile when monotonicity/convexity gates pass, otherwise a one-interior-knot cubic Hermite spline. The quintic method derives a quartic Q1 profile from a fifth-order potential and retains exact endpoint derivatives when accepted; otherwise it uses an adaptive two-interior-knot cubic Hermite construction for Q1 with exact endpoint derivatives. The underlying potential is C2 while Q1 is C1.

Endpoint matching alone is insufficient. The dense validation also checks Q1 positivity, the mechanical-stability numerator Q1+dQ1/dx and the causal ratio 0<c_s^2<1.

## 6. Assembly, output columns and diagnostics

The output grid contains: (i) all crust rows through the selected left knot i1; (ii) exactly n_bridge interior points; and (iii) all core rows from the selected right knot i2. Bridge endpoints are not duplicated. Source density knots outside the removed interval are retained.

The output rest-mass convention is that of the core table. The crust F values receive the fixed mass alignment. At every output point,
\[
 Q_3=\frac{F+Q_1}{m_n}-1,\qquad
 Q_6=Q_7=\frac{F}{m_n}-1,\qquad Q_1=P/n_B.
\]
Q4 on source portions is scaled to the output neutron-mass convention. Across the bridge it is joined with the quintic smootherstep s(u)=6u^5-15u^4+10u^3. The composition arrays and phase labels are zeroed because composition has deliberately been omitted.

Every constructed bridge is checked on 401 points, including both endpoints. The check rejects non-finite values, non-positive n_B, P, F or mu_B, non-increasing density, pressure or energy density, or any c_s^2 outside (0,1). The relative endpoint First-Law integral residual must not exceed 5e-9.

The report includes selected indices, endpoints, source/bridge endpoint sound speeds, minimum/maximum bridge sound speed, pressure span, interpolation/fallback mode and First-Law residual. The Python wrapper also evaluates the native C1--C11 consistency report on the assembled table. It returns the full report in info["consistency"]; info["diagnostic_only"] is true when the independent consistency check reports a failed hard gate. A bridge-level pass therefore does not certify every retained source knot.

## 7. Interfaces

### Python API

~~~python
from compose_eos import read_compose, match_crust_core, write_compose

crust = read_compose("eos_data_normalized/GPPVA_FSU2H")
core = read_compose("eos_data_normalized/GM_GM1")

matched, info = match_crust_core(
    crust, core,
    cr_name="GPPVA_FSU2H", co_name="GM_GM1",
    method="quintic_P", n_bridge=40,
    core_compo_mode="strict_thermo_only",
    drop_composition=True,  # required because this crust has eos.compo
)
write_compose("output/matched_compose", matched)
print(info["construction_kind"])
print(info["consistency"]["overall"])
~~~

The API preserves the legacy tuple convention (matched_eos, info). It raises NativeFeatureUnavailable for the unported literature_beta_eq option; invalid controls and missing composition consent raise ValueError. Native physical/input-state rejections are surfaced as NativeExecutionError.

### Native CLI

Build the executable by running the shell command shown below, then invoke the matching operation with the second command:
build/compose_eos_fortran match-core \
  CRUST_INTERNAL_DIR CORE_INTERNAL_DIR OUTPUT_INTERNAL_DIR \
  quintic_P 40 true
~~~

The three directories contain eos_data.csv and eos_data.nml. The output directory must already exist. The final argument is a required explicit boolean, true or false.

## 8. Verification and benchmark provenance

The regression suite compares native results with the historic Python reference implementation for all six tokens on the repository's GPPVA_FSU2H crust and GM_GM1 core. It compares selected transition indices and the assembled n_B, Q1, Q3, Q4, Q6 and Q7 arrays. It also checks output ordering, source-knot retention, composition omission, positivity, monotonicity, bridge causality and the First-Law residual.

The repository contains legacy matching metadata for this pair. Those records select
\[
 n_1=0.086999997\ {\rm fm}^{-3},\quad
 n_2=0.14954666667\ {\rm fm}^{-3},
\]
\[
 P_1=0.2901099811962072\ {\rm MeV\,fm}^{-3},\quad
 P_2=3.8788625430037915\ {\rm MeV\,fm}^{-3},
\]
\[
 \mu_1=953.407676579\ {\rm MeV},\quad
 \mu_2=979.3980753899036\ {\rm MeV}.
\]
The historical records report zero density jump and is_exact_maxwell=false for all six tokens. Their bridge sound-speed ranges lie inside the stable causal interval. These are legacy reference values, not a claim that this prose independently verifies the native arrays. The release audit records the native CI comparison at the exact tested commit.

Negative tests cover missing composition-drop consent, unavailable literature_beta_eq, nonzero entropy, Q6-Q7 mismatch, non-beta-equilibrated input and invalid controls. CI compiles the executable, runs the frontend/direct-provider tests, verifies the tracked SHA-256 manifest, and packages the exact tested source commit as a workflow artifact.

## 9. Limitations and non-equivalences

- The six smooth barotropic constructions are supported only in strict thermo-only mode. Composition reconstruction, charge-neutrality/lepton splitting, phase-label interpolation and literature symmetry-energy composition approximations are not ported.
- The public maxwell token does not impose coexistence at equal pressure and chemical potential with a density jump. It constructs a smooth near-pressure-plateau interpolation. There is no latent-heat discontinuity and the Seidov criterion is inapplicable.
- A genuine Gibbs construction requires a multi-variable thermodynamic surface, for example P(mu_B,mu_Q), and global charge-neutrality constraints. A one-dimensional beta-equilibrated barotrope is insufficient to reconstruct that surface.
- The css method is a two-segment continuous barotrope with piecewise constant sound speed. It is not interchangeable with all literature CSS parameterizations; some include a sharp first-order transition.
- Full master-function smoothing/candidate selection, microscopic RMF/Skyrme/APR solvers, a full causal-profile port, and complete native TOV-sequence/stable-branch diagnostics remain outside this feature.
- A bridge pass does not repair source-table inconsistencies outside the bridge. The C1--C11 report must be inspected, and astrophysical consequences of the selected matching method assessed separately.

## 10. Literature and authoritative resources

1. CompOSE Collaboration, [CompOSE manual and database documentation](https://compose.obspm.fr/manual/). Official definitions for table organization and thermodynamic/composition interfaces.
2. CompOSE Collaboration, [CompOSE software](https://compose.obspm.fr/software/). Official reader/interpolator software and supporting resources.
3. Servignat, Davis, Novak, Oertel and Pons (2023), [One- and two-argument equation-of-state parameterizations for neutron stars](https://arxiv.org/abs/2311.02653). Relevant to thermodynamically controlled EOS parameterizations and crust matching; it does not imply the six historical methods are unique or microscopic.
4. Canullán-Pascual et al. (2025), [Consistent crust-core interpolation and its effect on non-radial neutron star oscillations](https://arxiv.org/abs/2502.02373). A distinct linear P(mu_B)/constant-sound-speed framework; any density jumps in that framework must not be conflated with the smooth bridges implemented here.
5. Lindblom (2018), [Causal representations of neutron-star equations of state](https://scispace.com/papers/causal-representations-of-neutron-star-equations-of-state-27ptzdab5t). Motivation for explicit stability and causality constraints.
6. Swesty (1996), [Thermodynamically consistent interpolation for equation-of-state tables](https://doi.org/10.1006/jcph.1996.0162). Foundational work on interpolation based on thermodynamic potentials and derivatives.
7. Fortin et al. (2016), [Neutron star radii and crusts: uncertainties and unified equations of state](https://journals.aps.org/prc/abstract/10.1103/PhysRevC.94.035804). Relevant to the effect of crust EOS choices on stellar structure.

The historic implementation in this repository is the primary numerical reference for code-to-code differential testing. The literature motivates thermodynamic constraints and provides alternatives; it does not independently certify the historical interpolation heuristics.


## 11. Verified source snapshot

The six native strict-thermo-only matching methods were built and tested on branch `research/fortran-eos-gap-fill-v1`, final documented HEAD `e2b7ef319b9bb73dc85023a757dc1449c8f22e3c`. Both push GitHub Actions run [42378176837](https://github.com/cosmosapjw-quantum/NS_Cowling/actions/runs/42378176837) and pull-request run [42378180891](https://github.com/cosmosapjw-quantum/NS_Cowling/actions/runs/42378180891) passed on that exact commit.

The push run reports a successful native `gfortran` build, 44 frontend/direct-provider tests passed in 23.50 s, 58/58 tracked SHA-256 checks passed, source archive creation and artifact upload. The pull-request run independently reports 44 tests passed in 18.52 s with the same manifest and packaging gates successful.

The exact branch source snapshot is available as [`NS_Cowling-e2b7ef319b9bb73dc85023a757dc1449c8f22e3c`](https://github.com/cosmosapjw-quantum/NS_Cowling/actions/runs/42378176837/artifacts/13576661259), artifact ID `13576661259`, size 198,970,172 bytes, uploaded ZIP SHA-256 `28793bbbfa0708a8b24d0cc2e128b0833e9eea16e3d30ea547172b172b62474b`, retention expiry 2026-12-07. The artifact ZIP digest and the repository's 58-entry source manifest are distinct integrity checks.
