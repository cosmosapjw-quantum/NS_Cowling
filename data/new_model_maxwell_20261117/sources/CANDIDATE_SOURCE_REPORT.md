# Maxwell candidate source and model map — 2026-11-17

## Decision and bounded next implementation

Recommend a first regenerated homogeneous npemu cohort of **NL3, canonical PK1, and DDME-X**, using the precise parameter vectors in `PARAMETER_INPUTS.json` and preserving every named variant. Keep canonical TM1 as a fourth, explicitly source-diagnostic candidate because its public crust has the strongest direct forward-crossing evidence but the provider paper/file identity is inconsistent. PKDD was downloaded as a resolved fallback; its observed roots on available source support are reverse-oriented, so it is lower priority.

These are genuinely new parameterizations relative to existing native core models APR_APR, BBB_BHF2, GM_GM1, RG_KDE0v1, RG_SLY4 and existing GPPVA functional models DD2, DDME2, TW, FSU2, FSU2R, FSU2H, NL3wrL55, TM1e. The eight GPPVA core counterparts are useful role-new controls, but must not be counted as parameterization-new models. Plain NL3 differs from NL3wrL55 in rho coupling and absent omega-rho term. Plain TM1 differs from TM1e in rho coupling and absent omega-rho term. PK1 and DDME-X are new parameter vectors.

Scope executed: PhysMath 1.3.0 C2/R2 source/model/formula mapping, one source-curve reconnaissance pass. Native delegation was the parent dispatch to `/root/f_candidate_sources`; no recursive child. Main skill, POLICY_INTEGRATION, selected C2/R2 definitions and SUBAGENT_MODE were actually read from `c43/physmath-research-loop`. No EOS regeneration, Maxwell construction, native admission, canonical changes, commits or push occurred. All writes are below this `sources/` directory.

## Exact current provider identities

All provider pages state cold T=0, charge-neutral beta-equilibrated nucleonic n,p,e,mu matter, with a unified Thomas-Fermi crust and homogeneous core. The download directory uses **XMSLZ**, while the displayed provider name is **XMLSLZ**.

| Candidate | Official page | Exact archive | Release evidence | Homogeneous source rows, zero based | Source density range of core |
|---|---|---|---|---|---|
| NL3 | https://compose.obspm.fr/eos/284 | https://compose.obspm.fr/download/1D/NS/RMF/XMSLZ_NL3/eos.zip | Data sheet June 26, 2025; page cites NSXM_2025 update | [106,1078) | 0.058–2 fm^-3 |
| PK1 | https://compose.obspm.fr/eos/285 | https://compose.obspm.fr/download/1D/NS/RMF/XMSLZ_PK1/eos.zip | Data sheet June 26, 2025; NSXM_2025 update | [108,1078) | 0.062–2 fm^-3 |
| DDME-X | https://compose.obspm.fr/eos/281 | https://compose.obspm.fr/download/1D/NS/RMF/XMSLZ_DDMEX/eos.zip | Data sheet June 26, 2025; NSXM_2025 update | [116,1078) | 0.078–2 fm^-3 |
| TM1 | https://compose.obspm.fr/eos/287 | https://compose.obspm.fr/download/1D/NS/RMF/XMSLZ_TM1/eos.zip | Data sheet September 12, 2022; current page lacks 2025 update | [108,1078) | 0.062–2 fm^-3 |
| PKDD, reserve | https://compose.obspm.fr/eos/286 | https://compose.obspm.fr/download/1D/NS/RMF/XMSLZ_PKDD/eos.zip | Data sheet June 26, 2025 | [110,1078) | 0.066–2 fm^-3 |

Each candidate directory contains immutable provider_page.html, reference.zip, raw eos.* files, README, extracted eos.pdf.txt, and PROVENANCE.json with file hashes and actual archive members. There are 1078 rows in every downloaded slice. The table index's 2023 creation date is not the current data-release date. The current source bytes are frozen, not assumed identical to the 2022 numerical release.

Archive SHA-256:

- NL3: c23eb9064105619908267e34a2a1caa69de14972e6279f43c01e8ced03687e4f
- PK1: cb97a8b9bc680a5563b4b4b5943bbbda5cedda0331a6bdbe7210b81b8f0d3368
- DDME-X: 9082340b0b8b61a6eb1c4df806b02061e62117012cc750b301b7994a34ad0d35
- TM1: 82bdde688916e3a1180ca0ec9ec57227f395a33ee92e8b1d8516dad19429342d
- PKDD: 4ff5b520f663c482dc3b7de0ed5d6175927b4b01142116dd18b1b2eec2966f98

## Source identity findings that implementation must preserve

**Established by direct primary PDF comparison:** Xia et al., arXiv:2208.12893v1, Table I labels canonical TM1 parameters as PK1 and canonical PK1 parameters as TM1. Long et al., arXiv:nucl-th/0311031, Table I identifies the canonical vectors unambiguously. The 2025 provider paper arXiv:2506.11492v2 Table I repeats the canonical TM1 vector under PK1, and gives positive g3 for PK1r whereas Long Table I has negative g3. This is an actual tabular conflict, not a conjectured parser issue: downloaded PDF text itself was compared. Do not make the 2025 provider parameter table the sole SSOT.

**Observed source metadata:** current PK1 thermo header is mn=939.573, mp=938.280 MeV; current TM1 header is mn=939.5731, mp=938.2796 MeV, although canonical TM1 uses mn=mp=938 MeV. A CompOSE header is a normalization reference; this alone does not prove the thermodynamic curve is the wrong functional. Resolve numerical model identity by forward evaluation of the frozen primary vectors against raw P, epsilon, mu and composition at common homogeneous densities. Preserve both original fit masses and normalization header and use a representation-only rebase when necessary. Do not alter physical energy to improve a fit or create a crossing.

**DDME-X:** the original Taninah et al. arXiv:1910.13007 supplementary Table II was acquired from the observed arXiv source tar's `suppl.pdf`. It confirms ms=547.332728, gs=10.706722, gw=13.338846, gr=3.619020, bs=1.334964, cs=2.067122, co=1.605966, ar=.620220. It reports rho0=.151832 fm^-3, E/A=-16.097256 MeV, K=267.046210 MeV. Xia's model table uses n0=.152 and the provider data sheet reports rounded .152 with slightly different nuclear matter outputs. Keep an explicitly named provider-rounded reference-density variant (n0=.152) separate from an original-table rho0=.151832 interpretation. This difference is a source-identity diagnostic, not a pressure/energy adjustment parameter. The original main text defines x=rho/rho0; actual numerical comparison should resolve the adopted provider reference density.

**NL3 precision:** original Lalazissis et al. 1997 Table 2 is rounded (gs10.217, gw12.868, gr4.474, g2−10.431, g3−28.885, ms508.194). Long Table I independently gives the more precise vector used in Xia's table (gs10.2169, gw12.8675, gr4.4744, g2−10.4307, g3−28.8851, ms508.1941). Use one named vector consistently; the last printed digits are provenance, not freely fitted values. The original NL3 nuclear saturation controls and Long's comparative NL3 outputs are not identical, so use primary NL3 controls rather than treating every secondary comparative number as authoritative.

## Formula to implementation route

Use the newer `python_frontend/compose_eos/gppva_rmf.py` functional as the reusable reference, together with shared exact zero-temperature leptons. Its `model`, `couplings`, `scalar_potential`, `vector_potential`, `rearrangement`, `baryons`, `beta_state` functions map the relevant model. The present model registry only accepts its eight GPPVA models; create a parameter-injected bounded extension or explicit new-model module rather than monkeypatching global dictionaries. The old `rmf_solver.py` nonlinear sigma documentation/code has a questionable convention flipping both cubic and quartic coefficients. It is not the recommended implementation route.

Source convention: m_i*=m_i+g_sigma*sigma and rho vertex Gamma_rho*tau_i, with tau_n=-tau_p=1 in the 2025 paper. Current gppva convention: phi=-sigma, m_i*=m_i-gs*phi and rho vertex gr*tau3/2 with tau3_p=-tau3_n=1. Thus gr=2*Gamma_rho and the rho mean-field sign is reversed; physical products are unchanged. After phi=-sigma,

U(phi)=ms^2 phi^2/2 - (g2*hbar*c) phi^3/3 + g3 phi^4/4.

Therefore a3=−g2*HBARC and a4=g3; **the quartic sign is retained**. W(omega,rho)=mw^2 omega^2/2+mr^2 rho^2/2+c3 omega^4/4 for NL3/TM1/PK1 (no cross term). NL3 c3=0. PK1r requires a rho^4 term absent from current W and is out of this first cohort.

DD-MEX uses g_i(n)=g_i(n0)*a_i[1+b_i(x+d_i)^2]/[1+c_i(x+d_i)^2], i=sigma,omega, and Gamma_rho(n)=Gamma_rho(n0) exp[-a_rho(x−1)]. Include rearrangement Sigma_R=dg_omega/dn*omega*n+dg_rho/dn*rho*(np−nn)/2−dg_sigma/dn*phi*ns. Pressure includes n*Sigma_R; energy does not receive an independently added rearrangement term. These equations are already represented by the newer gppva code for DD2/DDME2/TW.

Solve n_p=n_e+n_mu, mu_n=mu_p+mu_e, and mu_e=mu_mu when muons are active; handle absent muons and baryonic composition endpoints. Epsilon includes baryon/lepton rest energies. Use fm^-3, MeV, MeV fm^-3 and explicit HBARC=197.3269804 MeV fm. Preserve fit masses exactly. The provider paper uses electron .511 and muon105.66 MeV while current shared code uses more precise masses; explicitly declare the implementation convention when assessing numeric identity.

For source-level forward comparison use raw homogeneous phase **0**, not generic phase3. Here Iphase mirrors geometrical dimension d: 3 droplet,2 rod,1 slab,0 uniform,−2 tube,−3 bubble. The old campaign's phase7 inner-crust selector is not portable to XMLSLZ. 2025 tables have Nadd=7: ce^2,cs^2,d,Z,A,Rd,RW; TM1 has Nadd=5: d,Z,A,Rd,RW. Preserve every original extra column and geometry boundary, using a declared lossless sidecar if the current one-extra-column loader needs a thermo-only projection. Do not interpret particle n,p fractions in a WS cell as free gas fractions.

## Actual P(mu) source reconnaissance and ranking evidence

`source_curve_probe.py` was actually run with exit code 0. `SOURCE_CURVE_SCREEN.json` contains all 198 overlapping component comparisons. It only inspects already available tables, with no source thermodynamic modifications, EOS regeneration, root-solving or matched-EOS output. Existing source components reuse exact prepared_v2 admitted intervals; new components are frozen per phase and split at any nonpositive dP/dn, dmu_E/dn, de/dn or superluminal source secant. The latter is a screening condition, not a certified continuous derivative test.

The pressure curves use mu_E=(epsilon+P)/n from physical raw epsilon=n*mn*(1+Q7), P=n*Q1. Raw source Q3 is separately retained. On these full slices max|mu_Q−mu_E| is 0.09375 MeV (NL3), 0.10288 (PK1), 0.11202 (DDME-X), 0.09599 (TM1), 0.10745 (PKDD); the maxima include high densities. This is another reason the source screen is not an acceptance result. Proper generated-source comparison must diagnose source rounding/consistency without silently rewriting raw values.

1. **NL3:** strongest clean first implementation. Its source vector is independently identifiable, it uses the existing sigma cubic/quartic path, and is distinct from NL3wrL55. The published homogeneous core starts at n=.058, mu_E=948.98050 MeV. No forward sign change appears against existing admitted crust components within this published core domain. However its public crust vs native GM1 core has a forward bracket mu_E=[943.4608815,943.5008200] MeV, deltaP=[−6.89945e−5,+2.80369e−4] MeV fm^-3 and midpoint densities ~.03155→.03978 fm^-3. This is cross-representation exploratory evidence only. Regenerating the homogeneous branch to lower n provides the relevant untested domain; do not extrapolate the public branch.
2. **DDME-X:** complementary density-dependent new model, reusing the rearrangement-tested DD functional. Source core starts .078, mu_E=953.74211 MeV, with no forward core crossing in that finite support. Its public crust vs native GM1 has a low-density forward bracket [942.3005827,942.3068567] MeV and delta_n~.01165 fm^-3. Core role implementation is low complexity once the declared n0 variant is frozen; source reproduction must resolve that variant.
3. **Canonical PK1:** reusable nonlinear vector self-coupling path with exact original Table I. The current homogeneous source starts .062, mu_E=950.93838 MeV; no forward sign change in the inspected support. Its chief current obligation is source numerical identity, because provider paper labels conflict. Canonical PK1 forward evaluation is executable without another literature search.
4. **TM1 diagnostic:** strongest same-representation direct prospect: public XMLSLZ_TM1 crust vs public GM_GM1 core has deltaP=Pcore−Pcrust changing from −1.66980e−5 to +2.19771e−3 MeV fm^-3 over mu_E=[945.4422345,945.6445049] MeV. Midpoint densities are ~.04550→.05494 fm^-3, delta_n~.00944. Against native GM1 the forward bracket is [945.2413989,945.2852410] and delta_n~.00948. These are finite raw-curve brackets, not certified roots. Resolve the source TM1 convention before claiming native reproduction of this particularly promising crust.
5. **PKDD reserve:** data and canonical parameters are available. Existing-support crossings with SPG_M4 and GPPVA_NL3wrL55 are reverse-oriented, e.g. public SPG_M4 bracket [953.9756434,954.1610243] MeV with delta_n~−.01855 fm^-3. Do not admit them merely because pressure equality exists.

For a forward sharp transition define deltaP(mu)=Pcore(mu)−Pcrust(mu). The desired local envelope orientation is deltaP<0 below and >0 above, and d(deltaP)/dmu=ncore−ncrust>0. A reverse crossing is a different phase orientation. Public finite-support absence, reverse crossings, source rounding and generator failure must remain distinct outcomes. No global impossibility or observational viability conclusion follows from this screen.

## Reusable verification obligations for parent

Before the matching campaign, verify the new homogeneous solver at fixed physical density against the chosen frozen source variant; check scalar/vector residuals, beta/charge constraints, epsilon+P=n*mu_B, numerical derivatives with refinement, and symmetric nuclear saturation targets. Keep provider equality separate from physical admission. Extend the homogeneous domain explicitly below its public phase onset if useful, freeze stable branches before comparing partners, and keep physical epsilon and mu fixed under representation rebasing. A positive root then needs the existing independent derivative, causality, ordering and envelope checks. The uniform generated cohort does not reproduce the provider's TF crust by itself.

Read-only baseline: `/home/cosmosapjw/ns_cowling_wt/stellar_extension_20261117`, HEAD `300e5a0532dd137cb09738568fca0f35cfcc9282`. Code-source identities and all artifact hashes are in SOURCE_LEDGER.json. Parent can begin the first three-model generation from PARAMETER_INPUTS.json immediately; no additional broad literature loop is needed.
