# Core EOS literature and reconstruction contract, 2026-11-16

Base commit 677585947e4d126bc27cc3d8f30c265c407a7773. This research covers the five intended cold core components and the three existing uniform DD-RMF functionals. Evidence was acquired through Exa, SciSpace, author code pages, CompOSE provider pages and local provider PDFs. The machine-readable source/model registry is data/research_generator_20261116/core_literature.json.

## Model-specific reconstruction

| Model | Primary implementation route | Composition | Source limits and remaining identity work |
|---|---|---|---|
| APR_APR | A18 + delta-v + UIX* APR Hamiltonian. Schneider et al. give all 21 parameters, LDP/HDP functions and effective masses. Author APREOS Fortran is available under GPLv3. | Solve cold neutral npemu beta equilibrium using derivatives of the full Hamiltonian. | Original CompOSE table ends near 1.34 fm^-3. Determine the actual causal crossing and phase coexistence on the chosen implementation. The cold APR provider and finite-T APREOS table have separately documented phase/lepton conventions. |
| BBB_BHF2 | Paris + Urbana TBF. Published SNM/PNM energies support a parabolic-composition reconstruction. Full microscopic recalculation requires a Bethe-Goldstone reaction-matrix kernel. | Paper Eqs.16–22 solve npemu, with neutron-proton mass difference omitted in mu_n-mu_p. Table3 supplies proton-fraction targets. | Paper states microscopic BHF stability up to about 0.8 fm^-3. Higher values use Eq.29 rational extrapolation. Published table ends at 1.5 fm^-3. Exact extrapolation coefficients are not printed in inspected paper. |
| GM_GM1 | Nonlinear sigma-omega-rho mean field using original coupling-ratio parameters. | Provider specifies npe, with proton fractions in core eos.compo. | Provider table ends at 0.999 fm^-3. Original TableII parameters are rounded. The scalar self-interaction mass convention and exact provider kernel require comparison. |
| RG_KDE0v1 | Uniform spin-saturated Skyrme functional, original Table3 KDE0(v1) coefficients. | Cold neutral npemu beta equilibrium. | Provider is EOS92, with homogeneous phase3. Table maximum 1.791349 fm^-3. Causal, mechanical and chemical stability must be located along the actual branch. |
| RG_SLY4 | Uniform spin-saturated SLy4 functional, published coefficients. | Cold neutral npemu beta equilibrium. | Cold provider EOS134 ends at 1.21 fm^-3. EOS123 is a finite-temperature NSE model and has a different target contract. |
| DD2 / DDME2 / TW | Existing density-dependent RMF uniform field solver with rearrangement self-energy. | Preserve actual configured n,p,e,mu species and beta-equilibrium residuals in outputs. | GPPVA crust-only tables provide a low-density benchmark. A high-density uniform comparison requires the corresponding homogeneous provider/model identity. |

Primary sources and provider metadata follow.

- APR original model is Akmal, Pandharipande and Ravenhall, Phys. Rev. C 58, 1804 (1998), DOI 10.1103/PhysRevC.58.1804. The provider identifies A18 + delta-v + UIX* and npemu. [CompOSE APR](https://compose.obspm.fr/eos/68)
- Schneider, Constantinou, Muccioli and Prakash, Phys. Rev. C 100, 025803 (2019), give the implementation equations and parameter Table I. The author distribution provides the APR source, user guide and example build dependencies Fortran/C/C++, HDF5, BLAS and LAPACK. [Paper](https://arxiv.org/abs/1901.09652) · [Author code page](https://stellarcollapse.org/APREOS.html) · [APREOS repository](https://bitbucket.org/andschn/apreos)
- Baldo, Bombaci and Burgio, A&A 328, 274 (1997), Table1 gives SNM/PNM energies and Table3 gives the Paris+TBF beta-equilibrium EOS and proton fractions. These published numbers are transcribed in bbb2_published_tables.json. The CompOSE page labels particles npe while the primary calculation includes muons, so the discrepancy belongs in provenance and comparison results. [Paper](https://arxiv.org/pdf/astro-ph/9707277) · [Provider](https://compose.obspm.fr/eos/69)
- Glendenning and Moszkowski, PRL 67, 2414 (1991), TableII first row supplies the GM1 coupling ratios squared in fm², 11.79, 7.149 and 4.411, with b=0.002947 and c=-0.001070. Later publications sometimes use slightly different rounded constants. The local provider thermo header specifies reference masses mn=939.56533 and mp=938.272 MeV/c². Those header masses alone do not determine every convention of the original nuclear kernel. [Author manuscript](https://escholarship.org/uc/item/6sw424k1) · [Provider](https://compose.obspm.fr/eos/54)
- Agrawal, Shlomo and Au, PRC 72, 014310 (2005), Table3 KDE0(v1) gives t0=-2553.08, t1=411.70, t2=-419.87, t3=14603.6, x0=0.6483, x1=-0.3472, x2=-0.9268, x3=0.9475 and alpha=0.1673. Replacing alpha by 1/6 changes the defined interaction. The fit includes a positive symmetry-energy slope through 3 times saturation density, which is a fit constraint rather than a full beta-matter causal-domain certificate. [Primary paper](https://arxiv.org/abs/nucl-th/0505071)
- SLy4 originates in Chabanat et al., Nucl. Phys. A635, 231 (1998). Douchin and Haensel Table1 lists t0=-2488.91, t1=486.82, t2=-546.39, t3=13777, x0=0.834, x1=-0.344, x2=-1, x3=1.354, alpha=1/6. Author HFODD output corroborates the usual kinetic coefficient hbar²/(2m)=20.73553 MeV fm². These rounding and kinetic choices must be explicit in comparisons. [Interaction DOI](https://doi.org/10.1016/S0375-9474(98)00180-8) · [Primary SLy4 table](https://ar5iv.labs.arxiv.org/html/astro-ph/0006135) · [Author HFODD output](https://www.fuw.edu.pl/~dobaczew/hfodd/pb208-a.out.txt)
- Gulminelli and Raduta describe the parent unified RG framework. Official KDE0v1 and cold SLy4 pages both state cold beta equilibrium. Their local provider PDFs list species indices n=10, p=11, electron=0, muon=1 and homogeneous phase3. A nonzero stored Q5 must be investigated against these conventions and actual fractions before assigning a physical non-equilibrium interpretation. [RG paper](https://arxiv.org/abs/1504.04493) · [KDE0v1 provider](https://compose.obspm.fr/eos/92) · [Cold SLy4 provider](https://compose.obspm.fr/eos/134)
- The density-dependent models have primary specifications in Typel et al., PRC81, 015803 (DD2), Lalazissis et al., PRC71, 024312 (DDME2), and Typel and Wolter, NPA656, 331 (TW). Rearrangement terms must appear in chemical potentials and pressure. [DD2](https://arxiv.org/abs/0908.2344) · [DDME2](https://doi.org/10.1103/PhysRevC.71.024312) · [TW](https://doi.org/10.1016/S0375-9474(99)00310-3)

## Composition and density-domain calculation

For cold neutral homogeneous matter, solve mu_n-mu_p=mu_e and n_p=n_e+n_mu. Populate muons only when mu_e exceeds m_mu c². With hbar*c in MeV fm, the ideal lepton density is

n_l(mu) = [mu²-m_l² c⁴]^(3/2) / [3 pi² (hbar c)³].

At a smooth composition minimum x(n), define epsilon_eq(n)=epsilon(n,x(n)). Then

P = n d epsilon_eq/dn - epsilon_eq,
cs²/c² = (dP/dn)/(d epsilon_eq/dn)
        = n epsilon_eq'' / epsilon_eq'.

For an interior minimum with a nonsingular stable composition Hessian,

epsilon_eq'' = epsilon_nn - epsilon_nx (epsilon_xx)^(-1) epsilon_xn.

This follows by differentiating epsilon_x=0 along the equilibrium branch. It includes the composition response required for beta-equilibrated sound speed. Frozen-composition derivatives address a different perturbation constraint and need their own label. Muon onsets and phase coexistence require piecewise treatment.

A density-domain report should contain the source tabulation bounds, declared approximation support, finite solver support, each connected stable/causal segment, crossing brackets and derivative/refinement evidence. A model can therefore retain a finite usable interval when its upper tail becomes acausal. An unresolved interpolation interval is recorded separately from a demonstrated model-level acausal interval. The applicable branch is evaluated with its actual species, mass convention and energy functional.

No numerical crossing density is adopted from an unrelated APR or Skyrme variant. All published table maxima in the registry denote available data support. Quantitative physical-domain results require the regenerated branch or the source-specific interpolation and pointwise tests.

## Remaining reconstruction work

APR can proceed through the public Hamiltonian equations or an external APREOS adapter. For cold comparison, add muons consistently and reproduce the selected LDP/HDP coexistence prescription. BBB can immediately support a published-table reconstruction with explicit interpolation error and composition checks. A new microscopic BBB calculation remains a larger kernel task. GM1, KDE0v1 and SLy4 have sufficient published interaction data for direct uniform-matter implementations. Exact agreement with a particular CompOSE archive additionally requires the provider constants, species and phase conventions to be resolved through the numerical comparison.

The overall all-model regeneration goal remains open until every intended component has an executable path, retained composition outputs, source comparison and a model-specific usable density domain.


APR follow-up: the official 2025-03-14 datasheet resolves the provider phase treatment as Gibbs; the current downloaded thermo is byte-identical to the local input. Actual low/high, Maxwell and two-charge Gibbs regeneration, author-kernel comparison, causal-domain localization and residual source discrepancies are documented in [apr_uniform_review.md](apr_uniform_review.md). Phase identity is no longer an unresolved input; exact bulk/data convention provenance still is.
