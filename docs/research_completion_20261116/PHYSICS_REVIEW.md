# Physics review — 2026-11-16

Scope: independent theory and finite numerical review of the five core implementations,
the candidate SPG spherical CLD functional, and the new discrete-mass outer-crust module.
Base commit: `c8ff2d6c97a93f1ccf621051861d3091bfe94d81`.
This report does not admit a global phase census, an authenticated Scurto surface fit,
a spatial GPPVA TF solution, or a combined core–crust/stellar result.

## Acceptance semantics

Physical implementation, source-model identity, numerical solution, and agreement with a
provider table are different claims. A mismatch with CompOSE is diagnostic evidence; it
is not by itself a physical-implementation failure. Conversely, removing a table-match
gate does not authenticate an undocumented calibration or restore omitted degrees of freedom.
A model's acausal extrapolation is to be reported by density interval and underlying
representation. It does not invalidate its independently checked lower-density functional.

The independent checks here use source equations, derivatives of the implemented
energy, composition constraints, boundary variations, and controlled alternate fits.
They do not use CompOSE pressure or composition as their pass thresholds.

## Core review

The new `core_regeneration.generate_core` independently evaluates the requested density
grid, with APR phase endpoints inserted. Its normalization is consistent with the native
cold CompOSE contract: Q1=P/nB, Q3=muB/mn−1, Q4=−mue/mn,
Q6=Q7=epsilon/(nB mn)−1; mn denotes a rest energy in MeV.
The charge chemical potential is −mue and the neutrino chemical potential is zero.
Free-neutron/proton and electron/muon fractions obey baryon and charge closure.
APR mixed states have a separate phase code; ordinary finite differences must not
silently cross a first-order boundary.

| Model | Independently executable physical representation | Limitation that must remain explicit |
|---|---|---|
| APR_APR | Published 21-parameter asymmetric APR branches; analytic species derivatives; bulk two-charge Gibbs construction | Bulk Gibbs neglects interface/Coulomb costs; it is an explicit construction choice. Author-functional agreement is not a full microscopic variational many-body calculation. |
| BBB_BHF2 | Published Paris+Urbana SNM/PNM energies, PCHIP energy interpolation, parabolic isospin, cold beta equilibrium including muons | Original BHF kernel is not re-solved. Microscopic source support ends at 0.8 fm^-3; published points up to 1.5 fm^-3 include source extrapolation. Rounded-table interpolation is part of the representation. |
| GM_GM1 | Nonlinear sigma–omega–rho RMF with published rounded couplings and g_rho/2 vertex; positive-effective-mass stationary roots ranked by energy | The declared output is npe. Adding muons is a separately labeled extension. Common nucleon rest energy is 939 MeV. Finite root bracketing is not a proof that no additional branch exists. |
| RG_KDE0v1 | Full asymmetric spin-saturated uniform Skyrme energy, including analytic composition derivatives | Uses published rounded parameters, alpha=0.1673. Uniform matter omits spatial/gradient, spin-orbit, Coulomb and pairing terms. |
| RG_SLY4 | Full asymmetric spin-saturated uniform Skyrme energy with alpha=1/6 | Same uniform-matter scope; high-density stability/causality limits are separate from equation implementation. |

Primary source links already encoded in these modules are
[APR functional](https://arxiv.org/html/1901.09652v2),
[BBB](https://arxiv.org/pdf/astro-ph/9707277),
[GM1](https://escholarship.org/uc/item/6sw424k1),
[KDE0v1](https://arxiv.org/abs/nucl-th/0505071), and
[SLy4](https://arxiv.org/abs/astro-ph/0006135).
This review inspected the code and source-contract mapping; it does not claim to have
independently rerun the original microscopic BHF/APR calculation.

Independent off-equilibrium test: 72 states, nB=0.113,0.237,0.517,0.913 fm^-3
and Yp=0.03,0.17,0.41, APR low/high evaluated separately. Five-point numerical
partial derivatives of energy density with respect to nn and np were compared with
analytic chemical potentials at two steps. Maximum absolute discrepancy was
1.802e-7 MeV; Euler residuals were retained. These are finite implementation checks,
not high-density physical-admissibility proofs. See
`data/research_completion_20261116/theory_review/core_derivative_review.json`.

## SPG constrained variational identities

Use nB, ni, ng in fm^-3; f is the dense-phase volume fraction; R is in fm;
yi and yg are dense/gas proton fractions. Set
ng=(nB−f ni)/(1−f), ne=f ni yi+(1−f)ng yg, and Delta=ni yi−ng yg.
With e2=alpha*hbar*c, the spherical functional is

```text
epsilon = f epsilon_i(ni,yi) + (1−f) epsilon_g(ng,yg) + epsilon_e(ne)
          + 2*pi*e2*f*Phi(f)*R^2*Delta^2 + 3*f*sigma(yi)/R.
```

This is the two-phase structure of
[Scurto et al. 2024, Eqs.21–25](https://arxiv.org/html/2402.15548v1).
In particular, Eq.22 explicitly contains the proton-density difference of the two
phases. Allowing yg>0 is supported by this source functional. The pure-neutron-gas
assumption in [Carreau et al. 2019, Eq.1](https://arxiv.org/html/1902.07032) is a
restricted ansatz and must not silently override the later source's two-phase form.

At fixed nB, direct Wolfram differentiation independently gives exactly the five
implemented gradients in (f,R,ni,yi,yg), the exact nB partial derivative, and

```text
R*d(epsilon)/dR = 2*epsilon_C − epsilon_S,
R_opt^3 = 3*sigma / (4*pi*e2*Phi*Delta^2).
```

The reduced coordinates t=f ni/nB and v=yg/yi preserve the baryon constraint.
The t=1 endpoint is physical absence of gas; its KKT condition is d(E/A)/dt<=0,
not zero derivative. The exact envelope chemical potential is

```text
B = mu_p,g − mu_n,g + mu_e − 4*pi*e2*f*Phi*R^2*Delta/(1−f)
mu_B = mu_n,g + yg*B + t*d(E/A)/dt
P = nB*mu_B − epsilon.
```

For an interior positive-gas stationary point B=0. For a pure-neutron gas yg=0.
At t=1 the correction involving the t derivative is essential; replacing muB by
the vacuum neutron mass alone gives the wrong outer-side pressure.

Both Wolfram exact residuals are saved with the evaluated code:
`cld_symbolic_verification.{wl,json}` and
`two_phase_symbolic_verification.{wl,json}`.
The former checks the original four-coordinate restricted functional; the latter
checks all five full two-phase derivatives and gives
`{{0,0,0,0,0},0,0}`.

Independent numerical check: 20 off-stationary reduced states across all five
SPG models, including t=1,yg=0. Interior gradient errors are below 1e-6 in the
declared coordinates; nB-envelope errors are below 4e-9 MeV. At ng=0 the gas
kinetic contribution scales as (1−t)^(5/3), so the one-sided first derivative
error scales as h^(2/3), rather than h^4. Observed step-halving ratios are
0.6297–0.6305, consistent with 2^(−2/3). Endpoint maximum error is 3.29e-5,
with the smaller-step error about 2.07e-5. This is an endpoint regularity effect,
not an assertion of fourth-order convergence.

Independent solves at M1,nB=0.0745 and M3,nB=0.078 fm^-3 find positive gas
proton fractions 2.2381e-4 and 1.5778e-5 respectively, projected gradients below
1e-12, and baryon/charge closure at printed floating precision.
These specifically support repairing the previously unstable restricted gas ansatz.
Evidence: `two_phase_numerical_review.{json,log}`.

The new Hessian checks are local curvature diagnostics after finite multistart searches.
A positive insertion margin at a pure-gas solution proves local KKT stability only;
an early return without searching a disconnected positive-gas branch cannot establish
a global minimum. No global-minimum or complete-phase-census claim is justified by
this finite review. Spherical geometry is itself an explicit source approximation.

Native composition detail: this repository defines Yh as the **cluster baryon
fraction**, t. The cluster number abundance is t/A and is recorded separately.
Yn=(1−t)(1−yg), Yp=(1−t)yg, Ye=t yi+(1−t)yg, and
A=4*pi*R^3*ni/3, Z=A yi. These satisfy the native C3/C4 constraints.

## Surface-fit source reconstruction and controlled experiment

Carreau2019 §II.2 defines n_sat as the saturation density of **symmetric** matter.
Its Eq.18 explicitly uses r_sat=(3/(4*pi*n_sat))^(1/3) in
E_S=4*pi*r_sat^2*A^(2/3)*sigma(Yp). It does not redefine that symbol to the
asymmetric saturation density before Eq.18. Its stated calibration nuclei are
40Ca,48Ca,48Ni,58Ni,88Sr,90Zr,114Sn,132Sn,208Pb.

The pinned AME2016 file contains **estimated 48Ni**, marked #, at source line442.
The eight other listed nuclei are measured. Consequently the existing all-2408
N,Z>=8 measured-mass selector does not contain the literal nine-nucleus cohort.
Neither a count match to2408 nor substitution of the older cohort authenticates the
Scurto author's selection/calibration implementation.

Scurto2024 Eq.31 explicitly describes N=2408 and 2%-of-mass scoring.
It states that the surface parameters are fitted per bulk model and cites Carreau
and later work; it does not supply enough detail here to infer that its actual
minimizer used precisely the earlier nine nuclei, fixed r_sat surface area, or this
repository's finite-nucleus compression convention. Its cell surface Eq.23 uses the
geometric R. Treating the fixed-area variant as the exact Scurto fit would exceed
the available evidence.

Constructive experiment performed for every SPG model:

- Same author-parameter GDFM bulk; vacuum unscreened Coulomb; no shell, pairing,
  curvature, exchange or atomic-electron binding corrections.
- Per-nucleus density minimization in [0.04,0.24] fm^-3.
- Nine nuclei including the explicitly estimated48Ni, and an eight-measured-nucleus
  alternative with48Ni held out.
- Both fixed symmetric n_sat surface area and variational geometric surface area.
- The same 2%-of-mass weights are held fixed to isolate the cohort/area changes.
  This is not a claim to have recovered Carreau's unspecified fitting weights.

All20 fits converge with no parameter-bound or density-bound hits. Fixed-area
nine-nucleus fits give sigma0=1.0105–1.0686 MeV fm^-2 and b=12.0336–29.6806.
Geometric-area nine-nucleus fits give sigma0=1.0280–1.0852 and b=6.0950–18.3377.
The four variants are separate candidate calibrations; the operational generator's
surface coefficients were not silently replaced by this review.

A161-to321-point density-grid repetition gives maximum relative sigma0 shift
7.70e-8, maximum absolute b shift4.95e-5, and maximum binding-energy prediction
change7.42e-5 MeV. Therefore the observed fit-prescription differences are much
larger than this density discretization effect. At Yp=0.01, one variant differs
from the prior all2408 geometric candidate surface by33.07%; the neutron-rich
extrapolation uncertainty is material even though finite mass fits converge.
There is no source-authenticated low-Yp calibration in this experiment.

Scripts, exact nine AME records, all fitted coefficients/predictions, optimizer
outcomes and grid comparison are under
`data/research_completion_20261116/theory_review/`.
The AME byte identity is
`2167f57a2a98331e4649b2dd2b658a9006ed4fba1975729ebfe52a42b4b9218a`.
This proves input identity, not cohort or scientific model equivalence.

## Outer-crust source check

Reviewed `outer_crust.py` against the locally acquired Pearson2018 AppendixB
math-alttext and its independent arbitrary-precision energy oracle.
The screening b2 numerator must contain
1+0.684 Z^5+0.222 Z^6, divided by 1+0.222 Z^6.
The initially missing final numerator term was identified by the implementing
agent and confirmed in this review; the corrected file agrees with the source.
The screening-pressure derivative correctly includes the density dependence of
the ion plasma term, and the correlation-energy derivative retains the base10
logarithm factor1/ln(10).

Minimizing g=(epsilon+P)/nB at common pressure over the finite mass catalogue
implements the appropriate cold catalysed discrete-nucleus variational problem.
It must not be replaced by minimization at common density. Two-sided layer
transition densities and neutron drip are explicitly separate. Native Yh=1
means all baryons are bound; cluster number abundance is1/A.

For a catalogue isotope without a source radius, setting the positive finite-size
energy to zero provides a lower bound on its fixed-pressure Gibbs energy:
pointwise e0(v)<=e_true(v) implies min_v[e0(v)+Pv]<=min_v[e_true(v)+Pv].
Thus an unknown-radius isotope that already loses with that lower bound cannot
become the winner after a positive radius correction. Refusing a selected
unknown-radius isotope is appropriate. This argument assumes the stated
equilibrium branch/range and does not prove completeness outside the finite mass
catalogue or below the declared computational/fully ionized regime.

The three arbitrary-precision source-energy derivative tests and the other outer
tests are owned by the outer-crust agent; this report does not duplicate their
entire campaign or claim a macroscopic physical surface at P=0.

## Review verdict and remaining obligations

**PASS_SCOPED** for the recorded core functional derivatives, the SPG constrained
two-phase identities and sampled endpoints, and the inspected outer-crust
equation correction/variational semantics. The checks are implementation evidence,
with the density/model/source limitations above.

**Unresolved source identity:** exact SPG finite-nucleus calibration prescription
and its neutron-rich extrapolation. The20 controlled fits at two resolutions are executable
sensitivity evidence, not recovered author coefficients. BBB remains a
published-energy reconstruction, not a rerun BHF calculation.

**Unresolved numerical completeness:** all-density convergence and phase census,
disconnected minima, non-spherical geometries where relevant, spatial TF
completeness for all intended GPPVA models. Their status belongs to the separate
campaign reports, not to this bounded review. No categorical rejection arises
here solely from a provider-table mismatch or a high-density acausal interval.


## Outer–inner reference and coverage check

The provider's
[SPG(M2) datasheet](https://compose.obspm.fr/download/1D/Crust/M2/eos.pdf)
and the
[DD2 datasheet](https://compose.obspm.fr/download/1D/Crust/DD2/eos.pdf)
state a BSk22 lower segment below nB=2e-3 fm^-3. This is well above the
approximately2.7e-4 fm^-3 physical neutron-drip density of the independently
computed BSk22 outer crust. The provider calls this whole lower segment outer
crust, but its post-drip part requires BSk22 inner-crust physics. The existing
archive inventory independently documents post-drip lower-segment composition.
A strictly pre-drip mass-catalogue calculation cannot regenerate that missing
post-drip interval by extrapolation.

There is also an energy-reference issue that cannot be hidden in table concatenation.
The SPG functional uses the declared common rest energy M0=938.9 MeV, whereas
the mass-catalogue outer crust uses physical neutron/proton masses. A common
shift epsilon'=epsilon+C*nB gives muB'=muB+C and P'=P. It changes neither
equilibrium composition nor the density of neutron drip (when its reference
threshold is shifted consistently). It therefore cannot by itself reconcile
different drip densities.

A physical species-dependent rest-energy restoration instead gives

```text
delta_epsilon = delta_Mn*nn_total + delta_Mp*np_total
              = nB*[delta_Mn + (delta_Mp−delta_Mn)*Ycharge],
delta_Mn = Mn_physical−M0, delta_Mp = Mp_physical−M0.
mu_n'−mu_p' = mu_n−mu_p + Mn_physical−Mp_physical.
```

The approximately1.293 MeV beta-equilibrium shift changes the solved composition.
Pressure is unchanged at fixed species densities, but the re-equilibrated pressure
need not be unchanged. This must be re-minimized, not added only to exported
energy columns. If the kinetic masses/couplings are retained, this is an explicit
rest-energy-corrected candidate model, not an assertion that the original author
used distinct kinetic masses.

The existing vacuum binding-energy fit is invariant under this linear restoration
when its convention is restored consistently:
B'=N Mn_physical+Z Mp_physical−E'=A M0−E.
The reduced CLD energy-per-baryon correction is
delta_Mn+(delta_Mp−delta_Mn)*[t yi+(1−t)yg], with corresponding nonzero
composition derivatives. Thus a consistent candidate calculation is feasible;
the need to re-solve it is a constructive task, not a universal no-go statement.

Two scientifically distinct completion routes remain:

1. Implement the source BSk22 post-drip lower-inner segment, then reproduce the
   explicitly declared density seam with separate thermodynamic diagnostics.
2. Construct a consistently rest-energy-restored mass/CLD candidate and determine
   the finite-model Gibbs competition at common pressure. If a crossing is absent,
   report its actual range and energies; do not manufacture a transition by an
   arbitrary pressure or energy splice.

The first is source-prescription regeneration; the second is an independently
specified candidate EOS. Neither follows merely from adding the current
standalone tables. This review has derived the transformation and established
the source lower-segment obligation; it has not claimed an executed complete
outer–inner construction.


## Executed rest-energy-restored seam experiment

The isolated construction described above was subsequently executed for **all five**
SPG models. No production coefficients or modules were modified. In
`run_rest_mass_seam.py`, the explicit linear contribution is included in the
energy **and its analytic composition gradient before minimization**:
delta_Mn=+0.6654 MeV, delta_Mp=−0.6280 MeV, with retained kinetic masses,
couplings and candidate surface coefficients. This changes the static energy
functional; it is not a mere change of plotting reference and is not promoted to
the original covariant GDFM parameterization with distinct kinetic masses.

At each density three initial configurations are minimized and polished against
the corrected analytic gradients. Gas-proton and muon insertion margins,
all three attempts, composition, and chemical potentials are retained.
The physical common-pressure comparison uses the **energy-derived** BSk22
pressure P=nB²*d(E/A)/dnB and chemical potential d(epsilon)/dnB.
The separately fitted source pressure is evaluated only as a sensitivity case.

The BSk22 production helper intentionally supports only its missing lower segment
through0.002 fm^-3. To locate the near-seam crossings for M1/M2, the isolated
experiment explicitly evaluates the published full-crust energy formula through
0.003 fm^-3, well below its source core-transition density. This does not change
the production helper's support guard or claim microscopic ETFSI regeneration.

The following first-order coexistences are found in this explicitly specified
candidate model, with P in MeV fm^-3 and densities in fm^-3:

| SPG model | n_BSk22 | n_restored_SPG | P_coexistence |
|---|---:|---:|---:|
| M1 | 0.002182394404 | 0.002329803032 | 0.002379705815 |
| M2 | 0.002215739676 | 0.002278216072 | 0.002415477273 |
| M3 | 0.001921727471 | 0.002107758444 | 0.002105908968 |
| M4 | 0.001585411343 | 0.001771564514 | 0.001767257277 |
| M5 | 0.001548148583 | 0.001742784680 | 0.001730670998 |

Each crossing has the correct density ordering, n_restored_SPG>n_BSk22.
At lower common pressure BSk22 has lower Gibbs energy; at higher pressure the
restored SPG candidate has lower Gibbs energy. The maximum chemical-potential
matching residual is5.69e-13 MeV. M1/M2 have no crossing within the source's
strictly declared BSk22 n<=0.002 segment; their roots require a modest explicit
extension of the published BSk22 fit. At the artificial common-density0.002
switch, neither P nor mu is equal. The density switch itself is therefore not
the equilibrium coexistence of these candidate representations.

Independent validation of all five crossings:

- Seven density probes between each pair of coexistence endpoints show both
  energies above their common tangent, within1.5e-15 MeV fm^-3 roundoff.
- Equilibrated-energy finite differences at0.001,0.002 and the SPG endpoint,
  with two steps, reproduce muB within4.75e-9 MeV.
- Every tested endpoint has positive reduced local curvature; the smallest
  scaled Hessian eigenvalue is0.417 in the declared coordinate scaling.
- Corrected beta/composition gradient residuals are below2e-12. The uncorrected
  composition gradient cancels the added linear term explicitly.
- Ye changes by0.00705–0.00881 relative to the original common-mass solution
  at the same endpoint density. This confirms actual re-equilibration; simply
  shifting exported energy could not produce this result.

The source-pressure-fit sensitivity also produces crossings, with shifted
locations, but its independent energy/pressure fits do not satisfy the exact
first law and are not used to admit a thermodynamic construction.
These results establish a **constructive candidate outer-inner interface route**;
they are not evidence that Scurto used this rest-energy prescription or that
the fitted source masses/calibration are exact. They do not resolve the separate
pre-drip mass-catalogue outer crust to BSk22 analytic-fit discrepancy at neutron
drip, and no complete hybrid EOS was exported.

Evidence:
`rest_mass_seam_experiment.{json,log}`,
`rest_mass_seam_validation.{json,log}`,
`run_rest_mass_seam.py`, and `validate_rest_mass_seam.py`.
The earlier0.002-limited experiment is retained separately as
`rest_mass_seam_0p002.{json,log}`.

## Final SPG table and outer-layer review

The retained final SPG campaign contains513 rows per model,2565 total.
An independent raw-file audit compares actual native thermodynamic columns with
the retained variational states. Maximum absolute differences are
5.56e-17 MeV fm^-3 in pressure,1.43e-14 MeV fm^-3 in energy density, and
1.14e-13 MeV in chemical potential. Baryon and charge composition residuals
are below1e-12. Seven probes across each of the original same-functional
CLD/uniform local-neutral common tangents put both energies above the tangent,
with minimum−2.05e-13 MeV fm^-3 at numerical endpoint roundoff.
The script now defaults to the committed campaign location
`data/research_completion_20261116/spg`, so this check does not depend on a
temporary build directory. Evidence:
`spg_final_independent_review.{json,log}` and `run_spg_v2_review.py`
(the script name retains its first review target).

The outer-crust agent's adaptive transition repair was also inspected.
When a third nucleus wins at the apparent crossing of two sampled endpoints,
the algorithm splits the interval and repeats full-catalogue fixed-pressure
minimization. Both phase-side densities are kept at the exact common pressure.
This addresses the observed thin-layer counterexample. Equal endpoint isotope
identities could still conceal an unbracketed reentrant layer, so finite pressure
census and resolution checks remain accurately scoped; no universal phase
completeness claim follows.
