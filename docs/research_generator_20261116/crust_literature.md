# Thirteen crust models: source intake and constructive regeneration routes

Date: 2026-11-16. Base: 677585947e4d126bc27cc3d8f30c265c407a7773.
Scope: independent crusts, with composition. Joined EOSs and stellar calculations are outside this node.

All thirteen archived targets contain composition. They are not intrinsically thermo-only. Exact source-table extraction and independent regeneration from the published model are different accomplishments. This intake removes several unnecessary obstacles to the latter: author-origin SPG bulk parameters are publicly available; GPPVA's source Thomas-Fermi calculation does not require Scurto's fitted surface tension; TW's common lower segment is demonstrably identical to DD2/DDME2 in thermodynamics and composition despite conflicting bibliography. No complete new inner-crust solver is claimed by this literature node.

## Evidence and reproducibility

The companion JSON files are:
- data/research_generator_20261116/crust_literature.json: sources, parameters, conventions and remaining dependencies.
- data/research_generator_20261116/spg_author_candidates.json: five candidate author model vectors, exact raw lines, row indices, source hashes and matching diagnostics.
- data/research_generator_20261116/crust_archive_evidence.json: thirteen archive intervals/hashes and six TW comparisons.
- data/research_generator_20261116/reproduce_crust_intake.py: reruns candidate matching and source-archive comparisons.

Executed from repository root:

```sh
python3 data/research_generator_20261116/reproduce_crust_intake.py
```

Result: 447346 author data rows; candidate indices M1=134057, M2=41210, M3=125583, M4=191537, M5=246826; strict printed-halfbin matches for M1/M3/M5; six outer-segment comparisons with zero mismatching rows; thirteen composition-bearing crust archives.

Large author files stay in build/crust_literature. Reacquire them from the URLs and verify SHA-256 values embedded in the script before reproducing:

| File | Bytes | SHA-256 |
| --- | ---: | --- |
| model_NMP.txt | 39925786 | 5cc222abd2d1f29054e9695f7f8d86c67d950949106cf8e9460ef6edf4ed7556 |
| model_par.txt | 59843103 | 43371e043b7cc572ad283032ce1785cae22e608039dd86cbee64f70ee5d70614 |

The [author Zenodo record](https://zenodo.org/records/11084025) contains README.txt, these two files, and an 18.47 GB plain-text eos.txt, not HDF5. README explicitly defines the latter's four columns as density, energy density, pressure and proton fraction; its outer-crust proton fraction is missing and encoded zero. Do not treat that zero as physical composition. README calls the set 447345 EOSs but gives model numbers 000000 through 447345; the acquired parameter/NMP files have 447346 aligned data rows. The very large EOS file was not downloaded.

## Domain first, then physical applicability

These are actual archived intervals, not a proof of universal model validity. All start at 4.67962e-10 fm^-3. The fourth column is the first row with nonzero provider phase ID. It marks a provider segment boundary; it is not the physical neutron-drip density.

| Model | Rows | Maximum archived n_B / fm^-3 | First nonzero provider phase / fm^-3 |
| --- | ---: | ---: | ---: |
| SPG_M1 | 365 | .07413102413 | .002511886432 |
| SPG_M2 | 377 | .08912509381 | .002290867653 |
| SPG_M3 | 366 | .07943282347 | .002630267992 |
| SPG_M4 | 375 | .09772372210 | .002630267992 |
| SPG_M5 | 363 | .07079457844 | .002511886432 |
| GPPVA_DD2 | 283 | .067000002 | .0020000001 |
| GPPVA_DDME2 | 288 | .071999997 | .0020000001 |
| GPPVA_FSU2 | 269 | .054000001 | .003 |
| GPPVA_FSU2H | 300 | .086999997 | .0020000001 |
| GPPVA_FSU2R | 296 | .082999997 | .0020000001 |
| GPPVA_NL3wrL55 | 298 | .082000002 | .0020000001 |
| GPPVA_TM1e | 303 | .089000002 | .0020000001 |
| GPPVA_TW | 289 | .074000001 | .003 |

A high-density uniform-matter acausality crossing does not invalidate a lower-density crust interval. At each point compare the regenerated and source model on the same temperature, composition constraint, phase and density domain. Record separately (i) source-model pathology, (ii) numeric differentiation or convergence error, (iii) convention mismatch, (iv) implementation error, and (v) use outside the source-supported phase/domain. A density discontinuity between equilibrium phases is not a smooth derivative sample: do not run a finite-difference sound-speed stencil across it.

## SPG M1-M5: author parameters recovered, source family corrected

[Scurto, Pais and Gulminelli 2024](https://arxiv.org/html/2402.15548v1), Sec. II.1 and Eq.20, uses sigma, omega and rho only; no delta meson. The common bare nucleon mass is 938.9 MeV. Each coupling has

g_i(n_B) = a_i + (b_i + d_i x^3) exp(-c_i x),  x = n_B / n_sat(model).

The twelve parameters and the three meson masses are in the author dataset. The recovered candidates have meson masses (560, 782.6, 769) MeV. They differ from the old four-meson Goegelein2008 proxy registry, which uses a fixed 0.16 fm^-3 normalization and a 0.178 fm^-3 saturation target. Source-model conclusions cannot be drawn from failures of those proxy entries.

The paper's model selection also corrects the old registry narrative: M1 maximizes posterior likelihood; M2/M3 represent the low/high extremes of the 68% maximum-mass interval; M4/M5 represent the low/high extremes of the 90% R_1.4 interval. M4 is not selected as a direct-Urca onset tail; indeed its source model does not admit direct Urca in the reported stellar sequence.

| Model | Author zero-based candidate index | Saturation density | Matching status |
| --- | ---: | ---: | --- |
| M1 | 134057 | .1572 | Unique within all nine printed halfbins |
| M2 | 41210 | .1551 | Unique within full printed bins; K differs from ordinary rounding |
| M3 | 125583 | .1612 | Unique within all nine printed halfbins |
| M4 | 191537 | .1621 | Unique within full printed bins; K differs from ordinary rounding |
| M5 | 246826 | .1483 | Unique within all nine printed halfbins |

The matching uses the nine Table9 values from n_sat through Q_sym, with halfbin tolerances .0005 fm^-3 and .05 MeV. M2 has author K=223.9837 MeV versus paper223.9; M4 has246.0713 versus246.0. These exceed ordinary halfbin rounding. The discrepancy is retained explicitly rather than concealed by a relaxed tolerance. Independent reproduction of the archived homogeneous P,epsilon,Y_p trajectory is the appropriate next model-identity check.

There is a separate normalization issue: paper Eqs.2,13,14,29 use the isovector vertex g_rho/2, while the author's raw numeric rho vector uses the no-half vertex Gamma_rho. The parallel implementation node has checked the two hypotheses against the author NMP file: using Gamma_rho directly reproduces J and L for all five candidates, while literally halving the raw numbers does not. Keep raw data untouched and record g_rho(paper)=2 Gamma_rho(data). The paper's displayed pressure expression also contains an omega-sign inconsistency relative to Eq.27; derive pressure from the energy and Euler identity including rearrangement instead of copying the inconsistent displayed line. Candidate extraction by itself does not validate a solver.

The root registry and rmf_solver module description have been corrected to identify the existing SPG branch as a legacy surrogate. Every pre-existing numeric registry value was checked unchanged. Existing inverse-fit receipts remain historical evidence; they are not reclassified as source fits.

### Remaining SPG finite-size work is specific

A subsequent constructive node has executed a separately labeled five-model,2408-nucleus candidate surface fit: see [candidate_surface_fit.md](candidate_surface_fit.md) and data/research_generator_20261116/candidate_surface/. It supplies explicit candidate sigma0,b values for a defined vacuum functional, not the original source coefficients. Substantial15-22MeV mass residuals, held-out results and small-y extrapolation limits are retained. Exact historical source membership was not treated as a prerequisite to this candidate experiment.

Eqs.21-25 define a spherical compressible-liquid-drop energy. At fixed total baryon density and charge neutrality, minimize over high-density baryon/proton densities, volume fraction f and droplet radius R, with low-density phase determined consistently by conservation and chemical equilibrium. Include density-dependent coupling rearrangement in both phases and every derivative of sigma(y). The two finite-size terms are

epsilon_C = 2 pi e^2 f Phi(f) R^2 (n_p^I-n_p^II)^2,
epsilon_S = 3 f sigma(y_p^I) / R,

Phi = (2-3 f^(1/3)+f)/5,
sigma(y) = sigma_0 (b+16) / [b+y^-3+(1-y)^-3].

For f>0, positive sigma and nonzero proton contrast, radius stationarity gives epsilon_S=2 epsilon_C and
R^3 = 3 sigma / [4 pi e^2 Phi (n_p^I-n_p^II)^2].
This is a derived elimination useful for implementation; it does not remove the composition derivatives. With e^2 in MeV fm, sigma in MeV fm^-2 and densities in fm^-3, the right-hand side has units fm^3.

The Zenodo parameter file contains no sigma_0 or b. Thus the remaining source-specific input is five pairs of surface-fit coefficients, or enough original fit specification to reproduce them. The already acquired AME2016 raw bytes and the candidate 2408-row selector alone do not determine those ten coefficients. [Carreau et al. 2019 Sec.II.3](https://arxiv.org/html/1902.07032) fits nine named spherical magic/semi-magic nuclei in its cited method, while Scurto Eq.31 scores 2408 masses with a stated 2%-mass weight. The exact relation between the surface-fit sample and this score remains unresolved. Resolve finite-nucleus energy terms, atomic/nuclear mass conversion, density minimization, shell/pairing/curvature treatment and weighting before calling a new optimizer source-equivalent. An explicitly labeled new calibration can be implemented once those choices are made, but cannot be passed off as the original coefficients.

After obtaining coefficients: solve the full CLD stationarity problem, compare independent finite differences/automatic derivatives, continue in density from multiple initial states, compare minima, and measure mesh/density convergence. The spherical restriction is part of this source model; adding pasta would change its approximation rather than complete a missing source feature. Source composition residuals must cover Y_e,Y_mu, dripped n/p and cluster/cell A,Z where defined, not pressure alone.

## Eight GPPVA models: finite-range Thomas-Fermi is the target

[Grill et al. 2014](https://arxiv.org/html/1404.2753) and [Pais and Providencia 2016](https://arxiv.org/html/1607.05899) provide the finite-range RMF Thomas-Fermi construction. Solve coupled meson-field and electrostatic boundary-value equations with local Fermi densities in Wigner-Seitz cells. Minimize total energy over cell size and competing sphere/cylinder/slab geometries; impose symmetry regularity at the origin, boundary conditions compatible with charge-neutral cells, fixed mean n_B and weak equilibrium. Density-dependent models need local rearrangement terms. Surface energy comes from meson gradients and the self-consistent profiles: **no imported Scurto sigma_0,b fit is required**. Earlier broad statements that all GPPVA models await a mass-fitted CLD surface are too restrictive and identify a different approximation.

| Target | Concrete bulk input route | Remaining source-specific item |
| --- | --- | --- |
| DD2 | Typel2010 exact rational/exponential parameters already in registry | Finite-range TF and source profile/composition definitions |
| DDME2 | Lalazissis2005 exact rational/exponential parameters already in registry | Same TF work; use source common mass convention |
| TW | Typel-Wolter1999 parameters already in registry with rho convention conversion | TF work; annotate conflicting lower-segment bibliography |
| FSU2 | Numeric nonlinear constants acquired from Negreiros2018 Table1 | Convert factorial/sign conventions and verify NMPs, then TF |
| FSU2R | Revised numeric constants from Tolos2017v2 Table1 | Revised Lambda=.045, not old.05; NMP/TF verification |
| FSU2H | Revised numeric constants from Tolos2017v2 Table1 | Nucleonic crust, despite hyperonic-core applications; NMP/TF verification |
| NL3wrL55 | Avancini2012 Table1: g_rho=11.2766, Lambda=.03, L=55.23 | Convert unscaled scalar-field convention and verify NMPs, then TF |
| TM1e | Complete TableI constants from author2020 guide | Map scalar sign and rho/2 convention, then TF |

Complete FSU2/R/H, NL3wrL55 and TM1e numbers are recorded in crust_literature.json, separate from executable registry parameters. FSU tables store g squared, not g. Scalar cubic/quartic factorials, omega quartic and mixed omega-rho normalization must be transformed together. For TM1e, the [author guide](https://user.numazu-ct.ac.jp/~sumi/eos/table4/guide_EOS4.pdf) uses M*=M+g_sigma sigma and an explicit rho/2 vertex; changing to M*=M-g_sigma phi requires phi=-sigma and changes the cubic sign, not the quartic sign. g2 has units fm^-1.

[Avancini et al. 2012 Table1](https://arxiv.org/html/1203.3052) resolves the remaining NL3wrL55 bulk input: g_sigma=10.217, g_omega=12.868, g_rho=11.2766, Lambda=.03, masses(508.194,782.501,763)MeV and M=939MeV. The rho coupling is retuned to S=25.68MeV at k_F=1.15fm^-1; Table2 gives J=31.66 and L=55.23MeV. Its scalar terms are -kappa phi^3/3! -lambda phi^4/4!, with kappa/M=4.384 and lambda=-173.31. These are unscaled scalar fields; the same symbols in the FSU2 convention multiply powers of g_sigma phi. Mixing the tables without this distinction is an implementation error.

[Providencia et al. 2019](https://www.frontiersin.org/journals/astronomy-and-space-sciences/articles/10.3389/fspas.2019.00013/full) provides numerical crust reference values and pasta transitions for FSU2/R/H. It reports no nonspherical pasta for FSU2, whereas FSU2R and FSU2H admit rod/slab intervals. The dynamical-spinodal n_t=.084 fm^-3 for NL3 omega-rho6 is not interchangeable with the TF archived cutoff .082. A spinodal is a linear-instability criterion, not a finite-size phase-minimum computation.

## Lower crust segment, BSk22 and TW

The source archives already provide the complete lower-segment thermodynamics and composition needed for replay. The live [TW page](https://compose.obspm.fr/eos/212) says BPS; its [provider PDF](https://compose.obspm.fr/download/1D/Crust/TW/eos.pdf) says BSk22. To test the practical consequence, the first217data rows were compared against DD2 and DDME2 for eos.nb, eos.thermo and eos.compo. After removing file headers and three leading grid-index fields where applicable, all six comparisons were token-identical. This supports extracting TW's actual archived segment with composition; the remaining conflict concerns historical attribution, not data availability or physical validity.

Independent BSk22 outer-crust reconstruction has a concrete public route. The [BSk22 provider datasheet](https://compose.obspm.fr/download/1D/NS/Skyrme/BSK22/eos.pdf) identifies AME2016 experimental atomic masses, theoretical HFB22 masses for unmeasured nuclei and Welker2017 copper-isotope updates. Minimize Gibbs energy per baryon at fixed pressure over nuclei with the correct electron and lattice corrections. AME2016 has been acquired; the exact historical HFB22 mass table and copper overrides still need intake. BRUSLIB's current landing page is now BSkG3, so blindly substituting its current default dataset would change the model.

A public [Chamel2020 outer-crust C solver](https://zenodo.org/records/3719439), CC-BY4.0, was acquired as build/crust_literature/ocrust.c (SHA-256 cfca995b9a39a81d367112e245eada177d054c605519cb70296591385453e924); it was not compiled or run in this node. It exposes electron-exchange, polarization, Wigner-Seitz and Thomas-Fermi polarization options. Its bundled example uses HFB27, not HFB22. It is a reusable algorithm/reference implementation after supplying the right mass input and source corrections. [Ioffe BSk fits](http://www.ioffe.ru/astro/NSG/BSk/nsesubr.html) also provide composition fits, but an analytic fit is not exact independent reproduction of the mass minimization.

For all sources retain the distinction between cluster A,Z and total cell A_cell,Z_cell, between dripped nucleon number fractions and baryon-weighted cluster abundances, and between provider phase IDs and physical pasta labels. In these exports Iphase=0/7 is not the same coding scheme as an extra pasta descriptor0..6. The lower provider segment extends beyond neutron drip; do not rename its entire217rows as strictly pre-drip outer crust.

## Completion conditions

Available now: thirteen composition-bearing source datasets; source-domain evidence; five explicit author bulk candidate vectors; raw-file hashes and reproducible matching; exact FSU2/R/H, NL3wrL55 and TM1e literature constants; an author outer-crust solver route; a direct TW data-consistency check; corrected legacy SPG provenance.

Remaining independent regeneration work: authenticate SPG candidate identities against homogeneous trajectories; close the SPG surface-fit definition and solve full CLD; implement finite-range TF for all eight GPPVA models; reproduce BSk22 lower-segment mass physics with exact historical inputs. These are identifiable constructive tasks, not a blanket rejection of the models. Source-table replay, new explicit approximations, and source-equivalent independent regeneration must each be labeled accurately.
