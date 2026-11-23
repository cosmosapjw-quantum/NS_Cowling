# EOS physical contract, source conventions, and native verifier audit — 2026-11-23

## Scope and evidence status

This is one complementary physics/source audit for the all-core × all-crust × all-methods campaign. The implementation base is commit `e51dadd132125f0c5389da8c92c50b84eeb524e8`, branch `research/eos-core-crust-all-methods-20261123`. The companion matching audit owns bridge implementation, candidate enumeration, and interval certificates. The campaign owns final source acceptance and survivor counts.

The statements below distinguish direct derivation, primary-source support, source-file inspection, numerical comparison, and implementation verification. A source-format conversion, a consistent sampled table, a certified bridge, and a globally admissible stellar EOS are different claims. The present document does not declare final campaign survivors.

## 1. Physical matching contract

### Definitions and assumptions

The relevant problem is a cold one-dimensional barotrope with conserved baryon number. For an electrically neutral, beta-equilibrated table including leptons, use baryon density \(n>0\), pressure \(P\), total energy density \(\epsilon\) including rest energies, and chemical potential \(\mu\). Units are \(n\) in fm\(^{-3}\), \(P,\epsilon\) in MeV fm\(^{-3}\), and \(\mu\) in MeV. The native variable named `cs2` is the dimensionless ratio \(c_s^2/c^2=dP/d\epsilon\).

At zero temperature,
\[
d\epsilon=\mu\,dn,\qquad \epsilon+P=n\mu,\qquad dP=n\,d\mu.
\]
For CompOSE `Il=0` at a fixed charge fraction \(Y_q\), the density-path conjugate is
\[
\mu_{\rm eff}=\mu_B+Y_q\mu_q.
\]
Consequently, raw Q3-derived \(\mu_B\) is not the denominator of \(dP/d\epsilon\) along this fixed-\(Y_q\) path. At `Il=1` in the beta-equilibrated scope, the effective potential is \(\mu_B\).

**Status:** derived thermodynamic identities under the stated assumptions. Multi-charge Gibbs mixed phases and non-equilibrium composition paths require their own potential and constraints; no such extension is inferred from a one-charge Maxwell check.

### Smooth bridges

Set \(x=\ln(n/n_\star)\), \(F=\epsilon/n\), \(q=P/n\), and \(D=q+q_x\). The first law gives
\[
F_x=q,\quad \mu=F+q,\quad \frac{dP}{dn}=D,\quad
\frac{c_s^2}{c^2}=\frac{D}{\mu}.
\]
A smooth bridge must preserve the physical source endpoint values and satisfy positivity, stability and causal bounds throughout its interval: \(n>0\), \(P>0\), \(\epsilon>0\), \(\mu>0\), \(D>0\), and \(\mu-D\ge0\). A nonnegative pressure surface may be handled separately; a negative-pressure source is not accepted by the native verifier.

In a \(P(\mu)\) representation,
\[
n=P'(\mu),\qquad \epsilon=\mu P'(\mu)-P(\mu),\qquad
\frac{c_s^2}{c^2}=\frac{P'}{\mu P''}.
\]
For positive \(\mu\), stability and causality are \(P'>0\), \(P''>0\), and \(\mu P''-P'\ge0\).

For ordered endpoints \(\mu_1<\mu_2\), causality implies
\[
\frac{n_2}{\mu_2}\ge\frac{n_1}{\mu_1},\qquad
\frac{n_1(\mu_2^2-\mu_1^2)}{2\mu_1}
\le P_2-P_1 \le
\frac{n_2(\mu_2^2-\mu_1^2)}{2\mu_2}.
\]
These follow by integrating the monotonicity of \(n/\mu\). They are necessary endpoint conditions, not proof that a specified interpolation family and its endpoint derivative constraints admit a solution. Failure for a finite candidate search is not a general physical impossibility theorem.

**Status:** derived; consistent with the causal matching bounds in Fortin et al., Section III and Eq. (16), [arXiv:1604.01944](https://arxiv.org/pdf/1604.01944), [Phys. Rev. C 94, 035804](https://journals.aps.org/prc/abstract/10.1103/PhysRevC.94.035804).

A bridge certificate does not certify retained source segments. The full matched artifact must also pass its source-domain, thermodynamic, interpolation, and causal checks. No mass-radius observation, maximum-mass threshold, or preferred transition-density band is introduced as an unstated matching criterion.

### First-order interfaces

A genuine one-charge Maxwell interface has common \(P_t\) and \(\mu_t\), positive jumps \(\Delta n>0\) and \(\Delta\epsilon>0\), and
\[
\Delta\epsilon=\mu_t\Delta n.
\]
The discontinuity is represented as an interface, not by differentiating a smooth stencil through it. A constant-pressure artifact without the coexistence identity is not a Maxwell interface. The Seidov condition concerns stability of a stellar configuration with an infinitesimal new-phase core; it is not an EOS thermodynamic matching veto.

A change of CompOSE normalization \(m_n\), with Q3/Q6/Q7 transformed to preserve the same physical \(\epsilon\) and \(\mu\), is a representation change. Adding \(n\,\Delta\mu\) to the physical energy density is a change of EOS. Fortin et al. explicitly discuss the latter correction; it cannot silently replace physical endpoint matching. Chemical-potential inconsistency can affect crust and stellar-radius estimates: [Zdunik, Fortin & Haensel, A&A 599 A119 (2017)](https://www.aanda.org/articles/aa/full_html/2017/03/aa29975-16/aa29975-16.html).

## 2. Constructive interval checks supplied to the matching implementation

### A sampled-causality counterexample

Let \(\mu\in[1,2]\), \(a=1201/800\), \(\delta=10^{-6}\), and
\[
n(\mu)=20\mu-a^3+\delta-3a\mu^2+\mu^3,\qquad
P(\mu)=1+\int_1^\mu n(z)\,dz,\qquad
\epsilon(\mu)=\mu n(\mu)-P(\mu).
\]
Then
\[
G(\mu)=\mu n'(\mu)-n(\mu)=a^3-\delta-3a\mu^2+2\mu^3.
\]
The exact minimum is at \(\mu=a\), with \(G(a)=-10^{-6}\). All 401 uniformly spaced sample points instead have \(G\ge6.033203125\times10^{-6}>0\).

Direct Wolfram evaluation returned:
- \(d\epsilon/d\mu-\mu\,dn/d\mu=0\) exactly.
- Minimum sampled \(G=6.033203125\times10^{-6}\).
- Exact off-grid \(G(a)=-10^{-6}\).
- \(c_s^2(a)/c^2=1.0000000503153099455\).
- Positive \(n\) and \(\epsilon\) on this interval.

**Status:** analytic counterexample, CAS-checked and numerically checked. It demonstrates the incompleteness of a finite sampled causality test; it does not allege failure of the already present exact polynomial-extremum certificate. The counterexample and extrema condition were supplied to the matching implementation owner.

### Pressure cubic satisfying the first law

For \(t\in[0,1]\), \(x=x_0+s\,t\), \(s>0\), and
\[
P(t)=a_0+a_1t+a_2t^2+a_3t^3,
\]
define
\[
A=P+\frac{P_t}{s}+\frac{P_{tt}}{s^2}+\frac{P_{ttt}}{s^3},\qquad
C=\epsilon_0+A(0).
\]
Then
\[
\epsilon(t)=C e^{st}-A(t)
\]
satisfies \(\epsilon(0)=\epsilon_0\) and \(\epsilon_t=s(\epsilon+P)\) exactly. The causal margin is
\[
G=\epsilon+P-\frac{P_t}{s}
 =Ce^{st}-\left(\frac{2P_t}{s}+\frac{P_{tt}}{s^2}
                    +\frac{P_{ttt}}{s^3}\right),
\]
with
\[
G''=Cs^2e^{st}-\frac{12a_3}{s},\qquad G'''=Cs^3e^{st}.
\]
Thus \(G''\) has at most one zero. When \(12a_3/(Cs^3)>0\), its candidate is
\[
t_\star=s^{-1}\ln\!\left(\frac{12a_3}{Cs^3}\right).
\]
If it lies inside the interval, partition there, isolate roots of \(G'\) by monotone bisection, and evaluate \(G\) at every such critical point and both endpoints. The polynomial pressure and stability margins require their own finite extrema checks. Degenerate \(C=0\) or \(a_3=0\), endpoint equalities, and small-\(s\) cancellation must be handled explicitly.

**Status:** directly derived and Wolfram-verified symbolic identities. This supplies the analytic search contract; final native certificate tests belong to the matching audit.

The Wolfram checks used generic formula queries only. No repository source or private tables were uploaded to literature/CAS services.

## 3. Raw NUDA table conventions

Pinned source package:
`nucleardatapy-1.0.2-py3-none-any.whl`,
at `/home/cosmosapjw/ns_cowling_wt/scratch/all_origin_maxwell_css_20261118/nuda/`.
Column indices in this section are one-based.

| Family | Number | \(n\) column | Energy-offset column | Supplied \(P\) column | Additive neutron-rest reference, MeV |
|---|---:|---:|---:|---:|---:|
| GMRS 2022 | 15 | 1 | 8 | 9 | 939.565346 |
| Gogny D1S, D1M, D1M* | 3 | 1 | 5 | 6 | 939.0 |
| BSk22, BSk24, BSk25, BSk26 ETFSI | 4 | 1 | 10 | 8 | 939.5651828 |

The physical conversion is \(F=E_{\rm offset}+m_nc^2\), \(\epsilon=nF\), using the source convention above; supplied pressure is retained. An Euler chemical potential \((\epsilon+P)/n\) derived during conversion is a constructed field, not an independent proof of thermodynamic consistency. The derivative/first-law checks remain necessary. No zero-density anchor or pressure replacement by a spline gradient is justified merely by converting the source format.

### Gogny

The raw headers identify energy per baryon relative to neutron rest mass. The pinned loader uses `mnuc2_approx=939.0`. This value is independently supported by the authors' own convention: Table II of [Mondal et al., arXiv:2003.03338v3](https://arxiv.org/html/2003.03338v3), [Phys. Rev. C 102, 015802](https://journals.aps.org/prc/abstract/10.1103/PhysRevC.102.015802), explicitly specifies a 939 MeV neutron mass; the Fig. 2 caption uses the same rest-energy reference. It should not be replaced with a modern physical neutron mass without transforming the actual EOS consistently.

**Status:** primary-author convention verified; pinned source-column inspection verified.

### GMRS

The raw table headers identify column 8 as \(E/A-m_n\) and column 9 as pressure in MeV fm\(^{-3}\). The pinned loader explicitly adds `nuda.cst.mnc2=939.565346`. The source reference is [Grams et al., EPJA 58:56 (2022)](https://doi.org/10.1140/epja/s10050-022-00706-w), [arXiv:2203.11645](https://arxiv.org/abs/2203.11645). The [CompOSE H2 description](https://compose.obspm.fr/download/1D/NS/Skyrme/H2/eos.pdf) is a primary provider document.

A direct read of [the CompOSE H2 thermo table](https://compose.obspm.fr/download/1D/NS/Skyrme/H2/eos.thermo) on 2026-11-23 returned header masses \(939.56535,\ 938.27201\) MeV and `Il=1`. The 939.565346 MeV pinned constant agrees with the rounded provider neutron mass to its printed precision. The exact extra digits are attributed to the pinned loader, not asserted as an independently quoted paper value. The current provider's complete raw table is not asserted to be byte-identical to every pinned GMRS table.

**Status:** pinned conversion convention verified; rounded primary-provider cross-check verified.

### BSk: column-identity defect and corrected physical interpretation

The upstream loader assigns one-based column 11 to a quantity named total energy per baryon. In these four files, column 11 is mass density in g cm\(^{-3}\). Using it as MeV/baryon is a dimensional error.

The table contains ETF energy in column 9, shell-corrected ETFSI energy per baryon relative to neutron rest mass in column 10, and neutron chemical potential relative to the same rest mass in column 16. The latter supports an independent Euler check:
\[
E_{\rm offset}+P/n-(\mu_n-m_nc^2)\simeq0.
\]
The total-matter interpretation and reference-energy convention agree with Eqs. (3), (4), and (6) in [Pearson et al., corrected article](https://arxiv.org/html/1903.04981), [MNRAS 481, 2994](https://academic.oup.com/mnras/article/481/3/2994/5090419).

The primary CompOSE thermo headers for all four models give
\[
m_nc^2=939.5651828\ {\rm MeV},\quad
m_pc^2=938.2718440\ {\rm MeV},\quad Il=1.
\]
Exact primary table URLs:
- [BSk22 thermo](https://compose.obspm.fr/download/1D/NS/Skyrme/BSK22/eos.thermo), [density grid](https://compose.obspm.fr/download/1D/NS/Skyrme/BSK22/eos.nb).
- [BSk24 thermo](https://compose.obspm.fr/download/1D/NS/Skyrme/BSK24/eos.thermo), [density grid](https://compose.obspm.fr/download/1D/NS/Skyrme/BSK24/eos.nb).
- [BSk25 thermo](https://compose.obspm.fr/download/1D/NS/Skyrme/BSK25/eos.thermo), [density grid](https://compose.obspm.fr/download/1D/NS/Skyrme/BSK25/eos.nb).
- [BSk26 thermo](https://compose.obspm.fr/download/1D/NS/Skyrme/BSK26/eos.thermo), [density grid](https://compose.obspm.fr/download/1D/NS/Skyrme/BSK26/eos.nb).

| Model | Pinned rows | Density range, fm\(^{-3}\) | Maximum independent Euler residual, MeV |
|---|---:|---|---:|
| BSk22 | 83 | 0.00027–0.0635018 | \(3.60615913042\times10^{-5}\) |
| BSk24 | 74 | 0.00027–0.0749994 | \(4.89920244036\times10^{-5}\) |
| BSk25 | 130 | 0.00027–0.0832615 | \(4.92169662873\times10^{-5}\) |
| BSk26 | 81 | 0.00027–0.0786984 | \(4.95415410171\times10^{-5}\) |

These residuals are at the scale of printed source precision. They establish the column interpretation; they do not certify every interpolation interval.

Pinned raw-file SHA256 identities:
- BSk22: `4afadadf9589d612cc6963d375f4bd74b9f736df0b01076dcaebf33e35badb32`.
- BSk24: `77d3032dbdbc053bb3b0f1651861c67a2bfe1cff6549b485eec1c1416d6b137f`.
- BSk25: `aefe67300cd71f8922b95125f054be96c83990cf25b9f6baaa0e294e5da1a0d2`.
- BSk26: `b2f9991c49f4538914f782846d6603921a72ac2352993085da88a66f3b7b14d5`.

### BSk source revisions must remain distinguishable

A same-density comparison against current CompOSE tables found:
- BSk22: all 83 pinned rows agree within \(4.60\times10^{-8}\) MeV in \(F\) and \(4.00\times10^{-8}\) relative in \(P\).
- BSk26: all 81 rows agree within \(4.64\times10^{-8}\) MeV in \(F\) and \(4.47\times10^{-8}\) relative in \(P\).
- BSk24: 56 unchanged rows, 17 numerically changed rows, and one pinned density 0.067894 fm\(^{-3}\) missing from the current table. Maximum changed-row \(|\Delta F|=0.003620028844\) MeV and relative pressure difference 0.003660323495.
- BSk25: 124 unchanged and 6 changed rows. Maximum \(|\Delta F|=0.002290002026\) MeV and relative pressure difference 0.001902529064.

The six changed BSk25 densities are 0.0768845, 0.0774996, 0.0781196, 0.0787445, 0.0800095 and 0.0806496 fm\(^{-3}\). Therefore the pinned NUDA files must be preserved and labelled as their own source revision; a conversion must not silently replace them with the current provider's values.

**Status:** direct source comparison, numerically checked. The discrepancy is provenance evidence, not proof that either full source is globally thermodynamically exact.

The [author's interpolation documentation](https://www.ioffe.ru/astro/NSG/BSk/nseintr.html) and [subroutine documentation](https://www.ioffe.ru/astro/NSG/BSk/nsesubr.html) also warn that independent analytic fits need not be perfectly thermodynamically consistent. Formula availability does not replace joint validation of \(P,\epsilon,n,\mu\).

## 4. Native verifier repair and actual regression evidence

Owned files:
- `fortran_eos/src/compose_eos_verify.f90`.
- `fortran_eos/tests/test_verify_domains_20261123.f90`.
- `python_frontend/tests/test_verify_reaudit_20261123.py` (compilation/execution harness only).

The old verifier silently set its derivative mask to false when \(d\epsilon/dn\le0\) or raw \(\mu_B\le0\). A corrupt source segment could disappear from the assessment while another segment supplied passing points. Native counterexamples demonstrated this: both a negative-energy phase and a positive-energy phase with negative \(d\epsilon/dn\) returned overall PASS. A thermodynamically consistent fixed-\(Y_q\), `Il=0` control instead falsely failed C11 because it differentiated raw \(\mu_B\).

The repaired entry checks allocation, shapes, finite primary arrays including Q2, metadata/tolerances, and a positive strictly increasing density grid before array arithmetic or logarithms. It enforces pointwise pressure/energy/effective-potential domains. Inside a regular resolved phase, invalid derivatives now contribute `DERIVATIVE_DOMAIN` and force overall failure; they cannot disappear into a passing subset.

For `Il=0`, C7/C9, the first-law sound-speed denominator and the monotonic-potential check use \(\mu_B+Y_q\mu_q\). The algebraic Q3 convention check still uses raw \(\mu_B\) in its proper CompOSE identity. Genuine Maxwell and provider phase-boundary stencils remain excluded from smooth finite-difference diagnostics. The one-sided endpoint formulas use three points, so boundary masking also covers an endpoint when an interface lies on interval 2 or \(N-2\).

C7–C9 remain finite-difference diagnostics, not arbitrary new exact-residual gates. Existing JSON fields remain compatible; the report adds `invalid_derivative_points` and a description of the differential potential. The legacy `min_dmu_b_dnb` field now refers to the effective path potential when `Il=0`.

Observed regression controls:
- Corrupt negative energy: previously PASS, now `PHYSICAL_DOMAIN`.
- Corrupt regular-phase negative energy derivative: previously PASS, now `DERIVATIVE_DOMAIN`, with two invalid derivative points explicitly counted.
- Consistent `Il=0`: previously false C11 failure, now PASS with C7/C9 satisfied.
- Valid smooth EOS and labelled/unlabelled Maxwell EOS: PASS before and after.
- Maxwell interface near an endpoint: PASS with all interface-crossing stencils excluded.
- NaN entropy, zero/duplicate density, missing or malformed arrays, negative pressure, nonfinite metadata and invalid tolerance: clean failure, without a bounds/FPE abort in the debug test.

All 16 native assertions passed with GNU Fortran 13.3.0 in both:
```
-O0 -g -Wall -Wextra -fcheck=all -ffpe-trap=invalid,zero,overflow
-O2 -fcheck=all
```

Both runs use `-std=f2008 -ffree-line-length-none`. Logs are retained under `/tmp/ns_eos_reaudit_20261123/verifier/{red,green-debug,green-opt}` pending the campaign's publication bundle. The debug and optimized run logs are byte-identical.

File identities at the verifier freeze:
- Verifier SHA256: `de5692d5f4e4c122f4ff01aa31abc84db7c827a31bd077c136d80738e431eaf1`.
- Native test SHA256: `dab19defb3b16738d28b0c2119ef714237d58501ff2f5cd98c24ca339ae005ae`.
- Each GREEN run log SHA256: `3040b17762669507bb56b86668a94682e0111d4c054510a24c671b9ffd599b27`.

Reproduce the automated harness with:
```
PYTHONPATH=python_frontend python3 -m pytest -q python_frontend/tests/test_verify_reaudit_20261123.py
```

**Status:** observed RED → GREEN, native implementation verified for the stated regression suite. The full rebuild, real-source readback, all-pair/all-method campaign, and published survivor manifest remain the parent campaign's separate evidence. A sampled native verifier does not itself prove continuum causality of every source interpolation or absence of every possible software defect.


### Automated wrapper and CAS provenance

The dedicated pytest wrapper was also actually run: **2 passed in 1.89 s**. Its log is `/tmp/ns_eos_reaudit_20261123/verifier/pytest.log`; SHA256 `202ac1ab497afafb79433c40f1f98a391d7dc3b2c6833cf64e1dcd234f755e03`.

The earlier successful Wolfram calls returned exact zero residuals and the numerical values recorded above. No persistent public computation URL was supplied by the connector; no such URL is invented here. A later attempt to regenerate one combined result returned an internal service error, and is **not** counted as an additional successful recheck. The exact replayable query and that failed retry response are preserved as `reproducible_cas_query.wl` and `cas_permalink_retry_error.json` in the verifier evidence directory. The analytic identities and the explicit rational counterexample above can be checked independently of connector availability.


## 5. Concrete completed-campaign outcomes

The completed Cartesian journal contains 15,190 requested cells: 8,712 admitted outputs, 5,827 cells with no admitted match, and 651 blocked-source cells. There were zero post-bridge gate failures requiring an unsearched alternative endpoint continuation. The two concrete source failures below were investigated with the frozen campaign executable, SHA256
`d1600879b3c741bb69f0acd150f05ff129048f244723933b1bb1e79518316588`.
These investigations made no production-code changes, source-number changes, or tolerance changes.

### 5.1 GPPVA_NL3wrL55 × PK1_CANONICAL: inherited source Q3 residual

Cell `cell_012222`, method Maxwell, selected `GPPVA_NL3wrL55/phase_000` and `PK1_CANONICAL/full_selected`. The native constructor found one pressure crossing and wrote 765 matched rows. The final native consistency check failed **C2 only**; C10 and C11 passed, with zero invalid derivative points.

A fresh reproduction with the same executable returned the same result. A standalone Fortran row inspector located the maximum residual at source-view row 1 and matched-output row 1. Both retain the original provider thermo row with density index 218:
\[
n=0.0020000001\ {\rm fm}^{-3},\quad
P=0.0018700000935\ {\rm MeV\,fm}^{-3},\quad
\epsilon=1.8796032939024001\ {\rm MeV\,fm}^{-3},
\]
\[
\mu_{\rm Q3}=940.7335764477001\ {\rm MeV},\qquad
\mu_{\rm Euler}=940.73659996112\ {\rm MeV}.
\]
Their difference is \(-0.0030235134198619562\) MeV, or a relative residual
\(3.2139851048496636\times10^{-6}\). The original provider columns are Q1 = 0.935, Q3 = 0.0018461943, Q6 = Q7 = 0.00085367408; the discrepancy is already present in that immutable source.

The Maxwell provider-input Euler threshold is \(5\times10^{-4}\); the assembled-EOS C2 threshold is \(2\times10^{-6}\). Therefore constructor success followed by final rejection is consistent with the two explicitly different gates.

**Classification:** `SOURCE_ALGEBRAIC_RESIDUAL_ABOVE_FINAL_TOLERANCE`. This is neither a finite-difference false failure nor an interface-stencil bug. Rejection is retained without replacing Q3 by an Euler projection or relaxing the threshold. It does not establish impossibility of the underlying nuclear model.

Original thermo SHA256:
`b476efb725b4d1b4a1e3bb8357ba5fb8a9c80c812eda776e06c6112f7e23fe41`.

Evidence is under `/tmp/ns_eos_reaudit_20261123/verifier/cell_012222_repro/`: `CASE_ASSESSMENT.json`, `MATCH_STDOUT.json`, `VERIFY.json`, `SOURCE_C2_ROWS.txt`, `MATCHED_C2_ROWS.txt`, and the Fortran inspector. The campaign trace is `campaign_v1/failure_traces/cell_012222.json.gz` under the audit data directory.

### 5.2 User MYN13 → OMEG1: adopted crust-source domain fails before bridge construction

All seven cells `cell_015184`–`cell_015190` were rejected under the adopted original-source domains. Their inputs remain private; only source identities and scalar/count diagnostics are recorded here.

| Source | Raw rows | Adopted source indices, zero-based half-open | Adopted rows | Adopted density domain, fm\(^{-3}\) |
|---|---:|---|---:|---|
| MYN13 RHF-TF crust | 448 | [0, 182) | 182 | \(6.0222\times10^{-11}\) to 0.06757 |
| OMEG1 core | 6,400 | [157, 6400) | 6,243 | 0.079 to 3.2 |

Source identities:
- Original private ZIP: `ca42416dfdab5a352139e14ce8a1b86cb6e2a516fcc9c0e07f48350ccd09271c`.
- MYN13 raw `RHF-TF.data`: `a4618b645b96695fb561d957caaf47bccebbeb3e687e443c3ab07cfa2f800092`.
- OMEG1 raw `EoS.d`: `76cefb7980974e6cdcaaae8ca84b6b5535df49e780fb77a5959a7ca0c511169d`.
- MYN13 adopted internal CSV: `6cdf195b24af699fc597b2915e3d49caf3221937cfcf70e67a0751ee91052a15`.
- OMEG1 adopted internal CSV: `440ebbead2fdcf72f61e3a6db7b5bc38fae45855189e20fc725b6a54a9c737d8`.

The frozen native Maxwell phase diagnostic and a standalone Fortran replay of the native retained-source/order/integral conditions found:

| Diagnostic | MYN13 adopted crust | OMEG1 adopted core |
|---|---:|---:|
| Nonpositive/nonfinite thermodynamic points | 0 | 0 |
| Pressure decreases beyond native tolerance | 0 | 0 |
| Nonincreasing energy-density edges | 0 | 0 |
| Euler-\(\mu\) decreasing edges | **24** | 0 |
| Causal secant failures | 0 | 0 |
| Maximum \(\Delta P/\Delta\epsilon\) | 0.010595647193585347 | 0.43402223923288169 |
| Native interval PASS | 15 / 181 | **6242 / 6242** |
| Native first-law/monotone-\(\mu\) interval failures | 24 | 0 |
| Native interpolation-shape interval failures | 142 | 0 |

The first MYN13 chemical-potential decrease occurs already between adopted rows **1 and 2**, with \(\Delta\mu_E=-0.0028688916867167791\) MeV. It exceeds native floating-point roundoff by many orders of magnitude. Consequently **no crust retained prefix ending at an eligible endpoint \(j\ge3\) passes**. All 6,241 eligible OMEG1 suffix endpoints pass.

Each of the six non-Maxwell methods examined all **1,123,380** source endpoint pairs and constructed zero bridges. The Fortran diagnostic found that, if the failed retained-source condition is ignored for diagnosis only, all those pairs satisfy the necessary endpoint ordering and causal integral bounds. Thus the operative rejection is the MYN13 retained crust source, not an OMEG1 acausal tail, pressure ordering, causal integral mismatch, or a particular Hermite/CSS bridge family. Maxwell independently rejects the low-density input for nonincreasing Euler chemical potential before searching coexistence.

#### Source precision sensitivity, without changing admission

A second Fortran diagnostic reads the raw decimal tokens and their last printed place \(u_z\), rather than inferring precision from binary floating-point representations. For independent bounds \(z\pm\alpha u_z\), it propagates
\[
\mu_{\rm lower}=
\frac{\epsilon-\alpha u_\epsilon+P-\alpha u_P}{n+\alpha u_n},
\qquad
\mu_{\rm upper}=
\frac{\epsilon+\alpha u_\epsilon+P+\alpha u_P}{n-\alpha u_n}.
\]
With \(\alpha=0.5\), corresponding to the usual nearest-rounding envelope, the first decreasing pair has overlapping chemical-potential uncertainty intervals, but **two** later decreases remain incompatible with any nondecreasing sequence inside those independent envelopes. The first such incompatibility is at row 22. With the wider \(\alpha=1\) envelope, no monotone-sequence incompatibility remains.

This demonstrates sensitivity to source precision. It does **not** establish the author's rounding policy or numerical error budget, and an allowed monotone chemical-potential sequence alone does not establish a thermodynamically consistent causal EOS. The wider envelope was not used as an admission tolerance, and no source entry was changed.

**Classification:** `ADOPTED_RAW_CRUST_MU_NONMONOTONIC_SOURCE_DOMAIN_REJECTION`. Keep all seven original-source outcomes rejected. The evidence does not support a general impossibility claim for MYN13 or for an independently reconstructed higher-precision model. Replacing the source by a higher-precision table or an explicitly labelled thermodynamic reconstruction would require a new source identity and separate validation.

Public scalar/count evidence and Fortran diagnostic source are under
`/tmp/ns_eos_reaudit_20261123/verifier/user_pair_repro/`:
`PUBLIC_ASSESSMENT.json`, `DOMAIN_SUMMARY.txt`, `RAW_PRECISION_SUMMARY.txt`, `MYN13_NATIVE_DIAGNOSTIC.json`, `OMEG1_NATIVE_DIAGNOSTIC.json`, `inspect_pair_domains.f90`, and `inspect_raw_precision.f90`.
The full per-interval values remain **private** under
`/tmp/ns_eos_reaudit_20261123/private_intake/user_pair_diagnostics/` and must not be copied to Git or the public evidence bundle.
