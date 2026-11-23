# Myers–Swiatecki 1966: formula generation and expanded BPS search

## Scope and source identity

The new models are BPS_Myers1966_magic14 and BPS_Myers1966_magic20.
Production mass generation, minimization, candidate filtering, electron/lattice
inversion, phase discovery and EOS construction are Fortran. Python only builds,
invokes and reads. This is implementation-verified for the published formula
and the declared finite domain, not a reconstruction of every UCRL-11980 Appendix D
entry or an experimental validation of the formula. Historical BPS1971_sequence
source/code remain unchanged.

Baseline: 0fbeeda62f078209cad03798c435933c34e2e1a1.
Branch: research/eos-original-papers-20261123.

Primary source: W.D. Myers and W.J. Swiatecki, Nuclear masses and deformations,
Nuclear Physics81 (1966)1–60, DOI 10.1016/S0029-5582(66)80001-9.
User attachment myers1966.pdf: 65 PDF pages, 9,460,361 bytes, SHA256
23d5031429ac8d635fca19db993972d1527b3be03a05f6647042faea01078561.

All substantive text was read. Printed pp1–58 contain the argument/appendices;
PDF pp59–61 are three foldout figures, pp62–63 references; final PDF pages have
no additional usable scientific text. Formula/parameter pages12,22,24,34,54 and
all three foldouts were visually inspected. The p42 empirical light-nucleus
table is a distinct alternative and is not merged into this model. No PDF is
added to Git.

The p2 footnote connects the journal paper with UCRL-11980. Printed p52 identifies
Appendix D of ref39 as the approximately8,000-nucleus table; ref39 names the
27May1965 report. Institutional locators:
https://escholarship.org/uc/item/65v5x6mx
https://escholarship.org/content/qt65v5x6mx/qt65v5x6mx.pdf
https://www.osti.gov/biblio/4514594 .
Actual report-PDF retrieval yielded an empty response and then HTTP403.
No report bytes/Appendix D tables were obtained, so no successful-download SHA
or table-reading claim is made. A user-supplied full UCRL-11980, especially
Appendix D, would resolve the exact-table/edition comparison. The journal DOI
does not itself supply that appendix.

## Equations, units and constants

Mass energies are MeV, lengths fm; N=neutron number, Z=proton number, A=N+Z.
The BPS kernel retains hbar*c and physical rest energies in dimensional constants.

From p13 Eq5, with M_(i-1)<=X<=M_i:

    F(X)=q_i(X-M_(i-1))-(3/5)[X^(5/3)-M_(i-1)^(5/3)]
    q_i=(3/5)[M_i^(5/3)-M_(i-1)^(5/3)]/[M_i-M_(i-1)]
    S=C[(F(N)+F(Z))/(A/2)^(2/3)-c A^(1/3)].

M_0=0; closures are2,8,14 OR20,28,50,82,126,184 for both species.
N=258 is additionally included for illustration (pp13,49). The paper explicitly
calls the high superheavy closures hypothetical. The two light-closure choices
define the two new model identities.

Equations6,9,10 (pp21–26) give the atomic mass excess:

    DeltaM=N Delta_mn+Z Delta_mH-c1 A+c2 A^(2/3)
           +c3 Z^2/A^(1/3)-c4 Z^2/A+delta+V(theta)
    c1=a1[1-kappa((N-Z)/A)^2], c2=a2[1-kappa((N-Z)/A)^2]
    alpha0^2=5(a/r0)^2 A^(-2/3), theta=alpha/alpha0
    E=[(2/5)c2 A^(2/3)-(1/5)c3 Z^2/A^(1/3)]alpha0^2
    F=(4/105)[c2 A^(2/3)+c3 Z^2/A^(1/3)]alpha0^3
    V=E theta^2-F theta^3 cos(3gamma)+S exp(-theta^2).

For theta>=0,F>0, gamma=0 minimizes V; the stable deformed solution is prolate.
Pairing p22: delta=-11/sqrt(A) for even-even, +11/sqrt(A) for odd-odd,
and0 for odd A.

| Constant | Published p34 value |
|---|---:|
| a1 | 15.677 MeV |
| a2 | 18.56 MeV |
| c3 | 0.717 MeV |
| c4 | 1.21129 MeV |
| kappa | 1.79 |
| C | 5.8 MeV |
| c | 0.26 |
| a/r0 | 0.27 |
| neutron mass excess | 8.07144 MeV |
| hydrogen mass excess | 7.28899 MeV |
| pairing coefficient, p22 | 11 MeV |

The source CSV accepts exactly these canonical decimal strings/keys; it is not
a fitting interface. NaN and overflow tokens such as1e999 are rejected before
real conversion. The tentative Wigner term/refitted constants on pp40–42 and
the empirical small-X shell prescription on p42 are not included.

## Stable nuclear deformation

The cubic expansion is unbounded below as theta goes to positive infinity.
An unrestricted global minimization would create runaway, nonnuclear masses.
The implementation instead enumerates stable nuclear local minima, compares
their energies and rejects absent bound shapes. It does not calculate a full
fission path or lifetime.

    V'=2theta f(theta)
    f=E-(3/2)F theta-S exp(-theta^2)
    V''=2E-6F theta+S(4theta^2-2)exp(-theta^2).

The sphere is stable for E-S>0. For S<=0, f decreases strictly, so there is
no nonzero stable minimum. For S>0, f' has two zeros only if
3F/(4S)<1/sqrt(2e). Bisect theta exp(-theta^2)=3F/(4S) on either side of
1/sqrt(2); any stable nonzero minimum is the unique upward zero of f between
these critical points. Its curvature is checked and its energy is compared
with a stable sphere. This covers the small first-order coexistence region
discussed in Appendix A rather than assuming E>S always selects the sphere.
For F=0, S>E>0, theta=sqrt(log(S/E)) and V=E[1+log(S/E)].
Negative cohesion and absent stable local shapes are rejected.

## Mass convention and finite candidate domain

To preserve published mass-excess differences while using the existing BPS
free-neutron rest energy, declare

    u_star*c^2=(mn*c^2)_BPS-Delta_mn(paper)
    M_nuc*c^2=A u_star*c^2+DeltaM_atom(paper)-Z me*c^2
    b=[N mn*c^2+Z mp*c^2-M_nuc*c^2]/A.

No separate bound-atomic electron binding correction is applied; any contribution
already encoded in the atomic mass fit remains, as in the historical BPS convention.
u_star is a reference-energy convention, not a measurement of the atomic mass
constant. The rounded historical1mMU=0.93144MeV conversion is not used as the
modern neutron rest-mass anchor.

The finite grid is2<=N<=258,2<=Z<=184; daughters at count1 are also computed.
The Fig19 p50 convention is Sn>=0 and Sp>=-3MeV, NOT Sp>=0. Candidates also
need stable local shapes and0<b<20MeV/baryon (the BPS primitive's domain).
Candidates with unresolved daughter shapes are excluded. Each variant yields
8,528 candidates from21,659 locally stable shapes. This is a finite formula
catalog, not a claim of global nuclear completeness.

## Optimized BPS equilibrium and complete phase sampling

At fixed pressure the electron/lattice pressure depends on n_e,Z, not A.
Invert P(n_e,Z) once per distinct Z, then evaluate all candidate masses:

    mu(A,Z;P)=[M_nuc*c^2+Z mu_e(P,Z)+(4/3)W_L(P,Z)]/A
    n_b=A n_e/Z.

The existing public bps_isotope_state supplies unchanged relativistic electrons,
bcc lattice, pressure, energy and derivatives. Positive compressibility is
checked. A safeguarded Newton solve in log(n_e) is bracketed; the relative
pressure failure threshold is1e-10. Neutron drip is the first global
mu_b=mn*c^2. The lower pressure equals the existing BPS pressure at
n_b=6.295e-12fm^-3 for56Fe; the new model's lower density uses its own nucleus.

For each nucleus dmu/dP=1/n_b and d2mu/dP2<0. On an interval each Gibbs secant
is a lower bound. Tangents to the endpoint-winning phases bound the global
minimum above. Comparing these piecewise-linear bounds at their kink excludes
hidden competing phases. For distinct endpoint winners, the sufficient
ordering n_low(P_right)<n_high(P_left) ensures a unique pair crossing.
Otherwise split the interval. Every crossing is checked against all candidates
and requires positive density/energy jumps. Unresolved intervals fail closed.

These are floating-point a posteriori bounds with a declared5e-11MeV Gibbs
comparison tolerance, not formal interval-arithmetic certificates. Phases whose
advantage lies below that numerical resolution are not claimed resolved.
Every transition emits both one-sided densities and every phase gets an interior
sample. The output row count may therefore exceed the requested grid count.

## Actual results and primary numerical anchors

For requested_points=129:

| Quantity | magic14 | magic20 |
|---|---:|---:|
| Candidate count | 8528 | 8528 |
| Selected phase count | 13 | 11 |
| Output rows | 166 | 160 |
| Lower n_b [fm^-3] | 6.256784860576778e-12 | 6.256784860576778e-12 |
| Drip n_b [fm^-3] | 2.587036333092080e-4 | 2.587036333092080e-4 |
| Drip P [MeV fm^-3] | 4.909067899495889e-4 | 4.909067899495889e-4 |
| Drip mu_e [MeV] | 26.19090383969095 | 26.19090383969095 |

magic14:52Cr,64Ni,66Ni,86Kr,84Se,82Ge,80Zn,78Ni,76Fe,124Mo,122Zr,120Sr,118Kr.
magic20:52Cr,50Ti,48Ca,84Se,82Ge,80Zn,78Ni,76Fe,122Zr,120Sr,118Kr.

Initial52Cr is the pure formula prediction. The paper recommends using measured
masses where possible (pp51–52) and notes light-nucleus discrepancies.
No measured masses are silently substituted.

myers_vs_bps1971_table2.csv compares new bindings with the12 rounded BPS Table2
anchors. For magic14, the final9 entries (84Se through118Kr) agree within
5e-5MeV/baryon, the half-unit printed rounding scale. Differences for
56Fe,62Ni,64Ni are -0.0128545915853,-0.00719214829086,+0.00227055862204
MeV/baryon. For magic20,76Fe differs by+0.00991298232751MeV/baryon.
This external agreement is useful; it does not prove the entire Appendix D
table or every edition/convention identical.

## Validation

A separate Fortran implementation evaluates Eq9 in real128 and locates local
wells by a dense scan/golden-section search, independently of production Eq10
stationary-root classification. For24 nucleus/variant cases:
maximum mass-excess difference4.121147867408581e-13MeV;
maximum theta difference2.442490654175344e-15.

At17 pressures per variant, a separate log-density bisection for every one of
8,528 candidates reproduces the optimized winner. Maximum mu difference
4.547473508864641e-13MeV; maximum density relative difference
1.953992523340276e-13. The reduced N<=184,Z<=126 domain and full domain give
identical winners at those34 pressures. All full-domain phase nuclei are
strictly interior. This tests saturation within the source domain; no
extrapolation beyond its hypothetical closures is certified.

Further tests cover the F=0 analytic limit, cubic-runaway rejection, competing
spherical/deformed minima, absent stable minima, exact shell closures, p50
approximately10.2MeV double-magic shell depression, invalid domain controls,
source corruption including1e999/NaN, output I/O failure, one-sided coexistence
pairs, grid independence16vs513 and byte identity OMP1vs4.

For magic14 initial generation: first-law finite-difference relative residual
6.27e-12, Euler residual3.41e-13MeV, transition delta-mu2.27e-13MeV,
pressure inversion relative residual1.97e-13. Preserved final JSON/logs give
current-build values. Tiny numerical residuals do not measure mass-model error. All selected phase
nuclei are spherical in these two runs; maximum catalog alpha is0.22883
(magic14) and0.31685(magic20). The source cubic expansion is still an
approximate small-deformation model for the wider candidate catalog.

A controlled single-run timing of the34 fixed-pressure queries measured
0.010650215s for cached Newton inversion and2.980707841s for independent
per-isotope bisection (about280 times). This excludes catalog construction
and output I/O and is not an end-to-end project speedup claim.

The parent's Wolfram derivative check returned zero residuals for the two
derivative identities. Its separate Limit expression had a Power::infy warning;
the native checks above are the executable admission basis.

## Interface and reproducibility

    call generate_myers_bps_native(source_path,n_points,eos,output_dir,status,magic_choice)

magic_choice14or20; n_points16..100000; success status=ok.
The caller supplies an existing output directory and exports the EOS carrier
through the shared application. Composition contains A,Z,Ye=Z/A; free n,p,muon
fractions are zero before drip. Labels start201 and remain distinct at transitions.

Outputs: myers_mass_catalog.csv, myers_bps_states.csv,
myers_bps_transitions.csv, MYERS_BPS_GENERATION.json.
JSON includes requested_points,rows,phase_count,magic_choice,n_min,n_max,
pressure_min,pressure_drip,finiteN/Zbounds and residuals.

Run from repository root:

    python3 -m pytest -q python_frontend/tests/test_myers1966_native.py

Fixture flags: gfortran -std=f2018 -O2 -fcheck=all,
-ffpe-trap=invalid,zero,overflow, -fopenmp. Python generates no production
mass or EOS number.

## Skill context used

PhysMath Research Loop1.3.0 at the parent's supplied actual plugin path.
Read SKILL.md, references/POLICY_INTEGRATION.md, SUBAGENT_MODE.md,
THEORY_LOOP.md, CODING_LOOP.md. Scoped task R2–R9 and C3–C8:
source reading, local-well derivation, actual implementation and bounded
independent numerical verification; no recursive reviewer tree.
Policy effects: preserve old models/raw source; stop failed report retrieval
after a diagnostic retry; implement the usable journal formula under explicit
new identities; repair a diagnosed test-case assumption without changing
tolerances; distinguish numerical accuracy from model accuracy.
Parent owns integration, review, full matching campaign and publication.
