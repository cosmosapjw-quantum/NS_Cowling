# Independent uniform Skyrme regeneration

Implemented two cold homogeneous functionals with explicit neutron, proton,
electron and muon fractions. The species closure is solved at every density.
The functional retains the complete asymmetric spin-saturated Skyrme terms,
including the density-dependent rearrangement derivative.

For n=n_n+n_p and tau_q=(3/5)(3*pi^2)^(2/3)n_q^(5/3), the binding functional is

E = hbar^2/(2m) tau
  + t0/4 [(2+x0)n^2-(2x0+1)(n_n^2+n_p^2)]
  + t3 n^alpha/24 [(2+x3)n^2-(2x3+1)(n_n^2+n_p^2)]
  + [t1(2+x1)+t2(2+x2)] n tau/8
  - [t1(2x1+1)-t2(2x2+1)] (n_n tau_n+n_p tau_p)/8.

Analytic mu_n and mu_p include the derivative of n^alpha.
Rest masses are added explicitly. Pressure uses sum(n_q mu_q)-epsilon.
Ideal relativistic leptons use the shared cold_leptons module.
Neutral beta-equilibrium roots on 0<=Yp<=1/2 are scanned and candidate total
energies compared. Spin-polarized/finite-wavelength phases require separate analysis.

## Parameter provenance

- KDE0v1 uses Agrawal, Shlomo and Au, nucl-th/0505071, Table 3,
  https://arxiv.org/abs/nucl-th/0505071 .
  In particular alpha=0.1673, preserving the published value.
- SLy4 uses Douchin and Haensel, astro-ph/0006135, Table 1,
  https://arxiv.org/abs/astro-ph/0006135 .
  The table publishes rounded t0=-2488.91, t1=486.82, t2=-546.39.
- The common fitted kinetic coefficient 20.73553 MeV fm^2 is corroborated
  by the author-maintained HFODD regression output
  https://www.fuw.edu.pl/~dobaczew/hfodd/pb208-a.out.txt .
  Both the coefficient and separate neutron/proton rest masses are serialized
  in each generation receipt. Exact provider parameter/convention equivalence
  remains an open comparison question.

## Executed results

The 38 targeted tests cover analytic chemical potentials against independent
energy differences, published symmetric saturation, lepton quadrature, baryon
and charge closure, first-law derivatives, density-step refinement and explicit
retention of high-density acausal rows.

The generator ran on all original homogeneous CompOSE density rows.
Native composition-preserving export/readback passed for both models.

| Model | Rows | Max relative pressure difference | Max relative total energy difference | Max absolute Yp difference |
|---|---:|---:|---:|---:|
| KDE0v1 | 274 | 0.0561895 | 0.000444813 | 0.000150633 |
| SLy4 | 54 | 0.1308852 | 0.000167586 | 0.000133048 |

The largest pressure differences occur near the first homogeneous rows,
n=0.04582369 and 0.05754399 fm^-3 respectively. Median signed relative
pressure differences are 0.0001453 and 0.0004131. These residuals are retained
per row. Their attribution requires the provider's exact constants, precision,
phase selection, and pressure construction. The source Q5 conventions are
recorded alongside the regenerated beta-equilibrium Q5=0 convention.

The numerically refined uniform-branch c_s^2/c^2=1 crossings are
KDE0v1 n=1.410965888 fm^-3 and SLy4 n=1.265005668 fm^-3.
Four central-difference step refinements give final changes 1.55e-8 and
9.87e-9 fm^-3. These locate the causal endpoint of the specified uniform
functionals. The source SLy4 comparison grid ends at 1.21 fm^-3, below this
crossing. Physical_domain reports retain both admissible sampled segments
and high-density violations.

## Reproduction

From repository root, after building the native executable:

~~~sh
PYTHONPATH=python_frontend python3 -m pytest -q python_frontend/tests/test_skyrme_uniform.py
COMPOSE_EOS_FORTRAN_ROOT=. COMPOSE_EOS_FORTRAN_EXE=build/compose_eos_fortran/compose_eos_fortran python3 scripts/build_skyrme_uniform.py build/new_skyrme_output
~~~

Executed output is under build/skyrme_research/generated_v4. Tests, refined
crossings, and raw logs are under build/skyrme_research. Each model directory
contains CompOSE files and a full GENERATION_RECEIPT.json with source-file hashes,
all composition rows, original-grid comparisons and domain diagnostics.

Retained output archive, raw logs, source hashes and compact summary are in
data/research_generator_20261116/skyrme/. Ten frozen-composition sound-speed
checks also satisfy c_frozen^2 >= c_equilibrium^2 at the sampled states.
For example SLy4 at n=1.0 fm^-3 gives 0.8368021 versus 0.8318619.

## Density-dependent source fidelity

A bounded sensitivity calculation retained every original homogeneous density.
The pressure difference is below 1% at 194/274 KDE0v1 rows, in sampled ranges
0.1137985--0.1314325 and 0.2102486--1.767846 fm^-3. The last KDE point
1.791349 fm^-3 differs by 1.057%. For SLy4, 40/54 rows meet 1%, comprising
the isolated n=0.1 point and the sampled range 0.2089296--1.21 fm^-3.
These are tabulated fidelity windows with gaps preserved.

At the worst low-density points the source/model pressure differences are
5.619% (KDE0v1) and 13.089% (SLy4). Muons are absent at these states.
Changing to a common 939 MeV nucleon rest mass changes model pressure by
approximately +0.787% and +0.780%; changing the kinetic coefficient to
20.734 MeV fm^2 changes it by -0.0199% and -0.0240%. Higher-precision
SLy4 coefficient variants change the latter pressure by -0.0014%.
Using hbar*c=197.32891 for leptons changes pressure by less than 0.0006%.
Substituting the source's rounded composition changes it by less than 0.018%.
These independently varied conventions account for only a small part of the
observed discrepancy. They were sensitivity calculations with no adopted retuning.

A three-point derivative of the source E/A grid reproduces source pressure
with maximum relative differences 0.0283% (KDE0v1) and 0.301% (SLy4).
The corresponding derivative of regenerated E/A on the identical grid
agrees with its analytic pressure within 0.0221% and 0.460%. Ordinary
finite-grid differentiation therefore does not account for the full low-density
difference. The exact provider functional constants/implementation remain the
next input for reconciliation. Source phase=3 and zero nuclear-group abundance
at the first compared rows are confirmed from the archived eos.compo and
provider eos.pdf phase legend.

The primary Gulminelli-Raduta paper explicitly includes neutron/proton rest
energies in beta equilibrium, Eq.37, https://arxiv.org/html/1504.04493 .
Its homogeneous-gas framework and the provider phase legend were checked.
The present sensitivity result supplies an unresolved reproduction discrepancy,
with its density dependence and tested explanations preserved in
data/research_generator_20261116/skyrme/pressure_fidelity.json.
