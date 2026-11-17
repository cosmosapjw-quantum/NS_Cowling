# Native EOS theory delta — 2026-11-16

Base commit: 4de49430560930be84d9a4147d71c3c9841da2db. This document covers SPG calibration and lower-crust construction; it does not authorize relabeling candidate EOSs as recovered original models. Physics remains a Fortran responsibility; Python may organize and inspect results.

## Result and evidence boundary

**Literature-supported:** the original SPG source fixes the GDFM bulk functional, spherical CLD expression, two-parameter surface form, AME2016 fit cohort size 2408, and the statement “2% of the corresponding mass”. It uses a same-density inhomogeneous/uniform energy comparison.

**New source delta:** Scurto's 2025 doctoral thesis supplies an explicit CLD derivative derivation (Appendix A), confirms the five SPG model identities (Section 3.3), and repeats the original fit statement (Eq. 3.6). It does not supply missing M1–M5 surface-coefficient values or an unambiguous finite-nucleus density/mass-reference fit contract in the sections inspected. A native implementation can be checked against its variational equations without pretending the author calibration has been recovered.

**Derived and symbolically checked:** the equations below give a constructive reduced CLD functional, its derivatives and beta condition, an efficient conditional calibration minimization, and thermodynamically consistent seam identities. They do not close missing empirical input.

**Unresolved:** source-equivalent sigma0,b calibration; original provider's BSk22-to-SPG joining operation; full BSk22 ETFSI reconstruction. No new numerical SPG EOS or physical transition result is claimed by this theory track.

## Sources and bounded delta inspection

Previously preserved evidence was read first:
- docs/research_finish_20261116/CALIBRATION_SOURCE_CONTRACT_KO.md
- docs/research_finish_20261116/SOURCE_COMPLETION_ADDENDUM_KO.md
- docs/research_finish_20261116/crust_assembly/REPORT_KO.md

New primary source:
- Luigi Scurto, *Neutron Stars: General Modelling within a Relativistic Mean-Field Framework*, Universidade de Coimbra doctoral thesis, issued 2025-09-16.
- Handle: https://hdl.handle.net/10316/120319
- PDF: https://estudogeral.uc.pt/retrieve/282004/PhD_thesis_Luigi_Scurto.pdf
- 5,410,325 bytes; SHA256 cf9ab5e4fc36d569975e3328c3d257e9e677e86a8b18f91ee0ef4072efd90929.
- Archive: data/research_native_20261116/theory/PhD_thesis_Luigi_Scurto.pdf.
- Repository metadata states openAccess and CC BY 4.0. The pdftotext derivative is not a substitute for the PDF if formulas disagree.
- Printed pp.22–24: Eqs.2.27–2.44, CLD surface/Coulomb functional and modified phase equilibrium.
- Printed p.51: Section3.3 specifies SPG(M1–M5) identities.
- Printed pp.105–111 (PDF pages119–125): AppendixA, fixed-total-density constrained derivatives.
- Eq.3.6 repeats 2408 masses and 2% mass; no additional calibration definition was recovered.

Primary original paper: https://arxiv.org/html/2402.15548v1 ; DOI https://doi.org/10.1103/PhysRevD.109.103015 .
Primary antecedent: https://arxiv.org/html/1902.07032 ; DOI https://doi.org/10.1140/epja/i2019-12884-1 .
Carreau2019 Eq.18 uses a symmetric-saturation-radius surface area and fits nine magic/semi-magic nuclei. Therefore citing that antecedent does not uniquely specify SPG's later 2408-nucleus fit.
Primary public source archive (already preserved): https://zenodo.org/records/11084025 ; the 16-column parameter file supplies saturation density, three meson masses and twelve coupling coefficients, not surface coefficients.
New conference-slide lead https://indico.cern.ch/event/1373806/contributions/5952951/attachments/2860432/5004440/Scurto.pdf repeated the surface-fit approach and did not recover coefficients.

Exa found and fetched the thesis handle and primary paper. Parallel independently located Carreau2019 and its primary arXiv/EPJA identities (search_f8c65efe7199bd6707ec433252220bb2). SciSpace found the SPG paper but its targeted extraction reported missing context even for equations present in the directly read primary paper; that output is a retrieval limitation, not evidence that those equations are absent.

## Constructive CLD calculation

Use MeV and fm, retaining hc in fermionic momentum integrals; here e2 means the Coulomb constant e^2/(4 pi epsilon0), in MeV fm. Temperature is zero. Let n and p denote total baryon and proton densities, x=n_I and z=p_I dense-phase densities, and f the dense-phase volume fraction. At fixed n,p,

    n_II=(n-f*x)/(1-f),  p_II=(p-f*z)/(1-f),
    y=z/x,  q=z-p_II=(z-p)/(1-f).

The domain includes 0<f<1, x>=z>=0, n_II>=p_II>=0. Source AppendixA.8–A.10 allows a proton-gas degree of freedom. Restricting p_II to zero is a separate active-boundary calculation; one cannot silently substitute interior stationarity equations there.

For droplets, Phi=(2-3*f^(1/3)+f)/5,

    sigma(y)=sigma0*(b+16)/B(y),
    B(y)=b+y^(-3)+(1-y)^(-3),
    epsilon_C=2*pi*f*e2*Phi*R^2*q^2,
    epsilon_S=3*f*sigma/R.

Minimization with respect to R yields

    R^3=3*sigma/(4*pi*e2*Phi*q^2),
    epsilon_S=2*epsilon_C,
    U=epsilon_S+epsilon_C=(3/2)*epsilon_S.

Thus the native minimizer may analytically eliminate R instead of solving a redundant dimension. Units are R[fm], sigma[MeV/fm^2], U[MeV/fm^3]. The elimination applies only for positive sigma,Phi,q^2 and interior positive radius. q=0 or f=0,1 must be handled as uniform/boundary limits; blindly evaluating the formula there is singular.

Define Delta=(1-y)^(-4)-y^(-4), S=epsilon_S after eliminating R. Then

    d sigma/dy = -3*sigma*Delta/B,
    (partial sigma/partial x)_z = 3*sigma*z*Delta/(B*x^2),
    (partial sigma/partial z)_x = -3*sigma*Delta/(B*x).

For U after radius elimination, with total n,p fixed,

    U_x=3*S*z*Delta/(B*x^2),
    U_z=S/[(1-f)*q]-3*S*Delta/(B*x),
    U_f=S*[3/(2*f)+Phi'(f)/(2*Phi)+1/(1-f)].

At fixed internal x,z,f, variation of total p gives U_p=-S/[(1-f)*q], whereas U_n=0. This latter distinction matters for beta equilibrium.

The full reduced energy density is

    F=f*epsilon_I(x,z)+(1-f)*epsilon_II(n_II,p_II)+U+epsilon_e(p).

Write mu_n^j=(partial epsilon_j/partial n_j) at fixed p_j and mu_p^j as the derivative with respect to the physical proton number at fixed neutron number. Then the interior stationary equations are

    0=f*(mu_n^I-mu_n^II)+U_x,
    0=f*[(mu_p^I-mu_n^I)-(mu_p^II-mu_n^II)]+U_z,
    0=P_II-P_I+z*(mu_p^I-mu_p^II)
      +(x-z)*(mu_n^I-mu_n^II)+U_f.

The surface derivatives cancel in the rearranged pressure condition, but not in either chemical-potential condition. For an internally stationary, beta-equilibrated interior point, the envelope theorem gives

    dF_min/dp|n = mu_p^II-mu_n^II+mu_e-S/[(1-f)*q] = 0,
    dF_min/dn|p = mu_n^II,
    P=n*mu_B-F_min,  mu_B=mu_n^II.

Imposing only mu_n^II-mu_p^II=mu_e would drop the Coulomb/surface contribution to the total proton derivative. At p_II=0 or another active inequality, replace these unconstrained equations by KKT derivatives of the chosen independent variables; the interior envelope formula is not a boundary certificate.

Wolfram symbolic evaluation independently returned {0,0,0,0} for the four U derivative residuals above under x>z>p>0, 0<f<1, sigma0,b,Phi>0. It returned {0,1,0,sigma0,0} for radius-virial residual, radius surface/Coulomb ratio, sigma-y derivative residual, sigma(1/2), and the conditional quadratic fit derivative. This is algebra verification, not a numerical convergence or global-minimum claim.

## Finite-nucleus fitting contract and constructive reduction

For a declared candidate fit, write

    E_i=A_i*e_bulk(n_i,y_i)+E_Si+E_Ci+E_Ki,
    R_i=(3*A_i/(4*pi*n_i))^(1/3),
    E_Si=4*pi*R_i^2*sigma(y_i),
    E_Ci=3*e2*Z_i^2/(5*R_i).

E_K denotes an explicitly selected curvature term; original SPG has a two-parameter surface expression, while CUTER's four-coefficient curvature fit is a different calibration. Choose one density rule in metadata, not implicitly:
1. symmetric saturation area (Carreau2019 Eq.18 antecedent);
2. asymmetric bulk zero-pressure density;
3. total finite-nucleus mass minimization.

For the third rule and a density-independent sigma(y), radius scaling gives

    dE_i/dn_i=A_i*de_bulk/dn_i+(-2*E_Si+E_Ci-E_Ki)/(3*n_i),
    P_bulk=n_i*(2*E_Si-E_Ci+E_Ki)/(3*A_i).

Bulk zero pressure is therefore not finite-nucleus stationarity. If extra density-dependent surface/curvature terms are added, their explicit derivatives must also be included.

At fixed b and fixed n_i, no curvature, set

    a_i(b)=4*pi*R_i^2*(b+16)/(b+y_i^(-3)+(1-y_i)^(-3)),
    d_i=E_bulk,i+E_Ci-M_target,i,
    w_i=1/s_i^2.

Then chi2 is quadratic in sigma0 and has the conditional optimum

    sigma0*(b)=-sum_i(w_i*a_i*d_i)/sum_i(w_i*a_i^2).

Apply declared positivity/bound constraints if needed. This eliminates one fit dimension exactly for the stated fixed-density contract. If n_i is reoptimized for every (sigma0,b), the global profiled objective is not this fixed quadratic; the formula is a conditional update only. A deterministic native fit can combine density minimization, bounded scalar b search and stationary checks, preserving the exact declared data, weights and mass convention.

A further source ambiguity must be retained. If the implementation defines B_model=A*M0-E and compares B_model-B_AME while weighting by 0.02*M_AME, its residual is not automatically the literal Eq.31 mass residual:

    E-M_AME = -(B_model-B_AME)+A*M0-N*Mn-Z*Mp.

Equality requires an explicitly justified remapping of the target mass to the same reference. “Mass weighting” alone does not establish mass-residual equivalence. The source text's 2% statement does not by itself resolve mass versus binding-energy convention, density relaxation, curvature, dataset cuts or optimization details. Neither conditional fitting nor output inversion can manufacture that missing provenance.

## Lower crust and seams without arbitrary shifts

**Source-policy clarification:** the common-tangent equations in this section diagnose an independently thermodynamic construction. They are not a mandatory gate for regenerating a provider-prescribed fixed-density attachment. The model-by-model original contract and minimal native construction are recorded in SOURCE_ATTACHMENT_CONTRACT.md; all thirteen pinned PDFs name BSk22 below nominal .002 fm^-3, while the stored source grids omit rows between the last outer point and the first inner point.

For a cold barotrope whose P(n) is supplied independently, the constructive first-law integration is

    epsilon(n)=n*[epsilon(n0)/n0+integral_(n0)^n P(t)/t^2 dt],
    mu_B=(epsilon+P)/n.

The anchor is physical input; an arbitrary anchor alters mu_B. This exactly enforces the first law for the chosen P curve and anchor, but it does not restore the missing BSk22 ETFSI functional or prove its analytic composition is in chemical equilibrium with the integrated thermodynamics. The preserved BSk discrepancy up to about 0.0255 MeV remains relevant.

For two phases represented in a common energy convention, a zero-temperature neutral Maxwell seam needs

    P_L(n_L)=P_H(n_H),  mu_L(n_L)=mu_H(n_H),
    epsilon_H-epsilon_L=mu_t*(n_H-n_L).

A common shift epsilon->epsilon+a*n leaves P unchanged and shifts every mu_B by a; it cannot create a missing crossing. More generally a conserved reference epsilon->epsilon+a*n_B+b*n_Q induces

    delta mu_n=a, delta mu_p=a+b, delta mu_e=-b,

so delta(mu_n-mu_p-mu_e)=0. Adding separate neutron/proton rest-mass corrections to baryons only instead changes beta equilibrium by delta_mn-delta_mp; it is a different physical approximation, already correctly labeled as a separate candidate in previous work.

The measured seam sensitivity can be stated without fitting the mismatch away. With G(P)=mu_H(P)-mu_L(P),

    G'(P)=1/n_H-1/n_L,
    deltaP_t=(delta_mu_L-delta_mu_H)/G'(P_t).

Small numerical root residuals are not physical accuracy when delta_mu branch-model error is nonzero, particularly as n_H approaches n_L. If no crossing is found inside the actual common support interval, extend support through a specified model calculation or retain a gap. A spline bridge or branch-specific energy offset cannot be labeled as a derivation of the original EOS.

A whole-source-consistent lower crust would require either the original BSk22 ETFSI/mass functional and provider join specification, or a clearly new unified lower branch built from each selected EDF with an explicit nuclear mass approximation. The latter is constructive research, but it changes the source model and must have a separate identifier. This is a bounded missing-input result, not a proof that the original construction is impossible.

## Native implementation acceptance and remaining decision

The theory supplies direct tests for the Fortran reduced functional: radius virial identity; analytic versus finite-difference surface derivatives; full constrained beta derivative; active-boundary KKT signs; first-law residual from regenerated samples; charge neutrality; support and phase-boundary metadata. Compare independent physics, source fidelity, numerical support and CompOSE table differences as four separate claims.

No CompOSE discrepancy alone is a failure of the derived functional. Conversely, a close CompOSE match does not supply author sigma0,b provenance or prove a thermodynamically valid lower-crust join.

The present actionable resolution is to make the native candidate solver and its calibration contract explicit, add the recovered source stationary equations to independent verification, and leave exact-author calibration plus full original lower crust unresolved until the missing input is supplied or a new, explicitly named functional is completed.
