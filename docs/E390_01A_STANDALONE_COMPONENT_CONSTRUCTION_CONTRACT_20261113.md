# E390-01A: standalone core-only and crust-only EOS construction contract

**Date:** 2026-11-13  
**Node:** E390-01A (M1A)  
**Contract:** E390-01A-STANDALONE-COMPONENT-CONSTRUCTION v1.0.0  
**Status:** specification complete; source-boundary admission remains a separate blocking gate  
**Repository:** cosmosapjw-quantum/NS_Cowling  
**Source registry:** \`data/eos390_component_source_contract.v1.json\`

## 1. Decision and scope

This node fixes the scientific and data contract for standalone components before either component builder or any crust-core matcher is implemented. It does not generate EOS files, execute a matching method, or run the 390-task matrix.

The current "core" source cohort is not five homogeneous-core-only files. CompOSE identifies APR(APR) and BBB(BHF-BBB2) as non-unified cold EOSs with BPS outer crust and SLy4 inner crust already appended; GM(GM1) contains a Douchin-Haensel crust below \(n_B=10^{-3}\,\mathrm{fm}^{-3}\); RG(KDE0v1) and RG(SLY4) are classified as unified EOSs. Treating each file's minimum density as the start of a homogeneous core would therefore import crust matter into a core-only artifact. This is the central reason for separating raw source coverage from physical component support.

The 13 crust inputs are also not interchangeable density tails. The five SPG models use their model-specific crust construction, with a BSk22 low-density segment and CLD-based inner-crust treatment; the eight GPPVA models use BSk22 at low density and self-consistent Thomas-Fermi inner-crust calculations, including nonspherical pasta where specified. Each crust EOS has its own terminal density. That terminal point is table support, not a universal transition density.

## 2. Source-authoritative component matrix

The intervals below are the exact native \`eos.nb\` grid support in the current branch source files. They are not assertions that every row is physically a core row or that the upper crust endpoint is the exact thermodynamic transition. Git blob identifiers are recorded in the JSON registry as Git SHA-1 identities; the authoritative archive/per-file SHA-256 lineage remains the E390-00 source receipt and inventory.

### 2.1 Core-source cohort

| Model | Raw \(n_B\) support [fm\(^{-3}\)] | Points | Source construction already in file | Core admission |
|---|---:|---:|---|---|
| APR_APR | \(7.9241\times10^{-15}\) to 1.34 | 172 | APR core + SLy4 inner crust + BPS outer crust | **Blocked:** exact APR/SLy4 splice locator not yet admitted |
| BBB_BHF2 | \(7.9241\times10^{-15}\) to 1.50 | 84 | BHF core + SLy4 inner crust + BPS outer crust | **Blocked:** exact BHF/SLy4 splice locator not yet admitted |
| GM_GM1 | \(7.9241\times10^{-15}\) to 1.00 | 349 | GM1 core + Douchin-Haensel crust below \(10^{-3}\) | Source states \(n_B=10^{-3}\) as the crust splice |
| RG_KDE0v1 | \(10^{-7}\) to 1.791349 | 1240 | Unified crust-core table | **Blocked:** KDE0v1-specific phase legend/core onset must be extracted |
| RG_SLY4 | \(10^{-7}\) to 1.21 | 197 | Unified SLy4 table with cluster energy functionals | **Blocked:** SLY4-specific phase legend/core onset must be extracted |

CompOSE source pages: [APR(APR)](https://compose.obspm.fr/eos/68), [BBB(BHF-BBB2)](https://compose.obspm.fr/eos/69), [GM(GM1)](https://compose.obspm.fr/eos/54), [RG(SLY4)](https://compose.obspm.fr/eos/134), and the [CompOSE model selection table](https://compose.obspm.fr/table_selection) for the RG(KDE0v1) unified classification.

Consequences:

1. APR and BHF-BBB2 are non-unified inputs whose low-density rows already contain a different crust model. The core builder must discard those rows at the documented model-specific splice. Both CompOSE abstracts state that no composition information is available; a barotropic core may still be constructed, but composition must remain \`NOT_AVAILABLE\` unless a separate, explicitly sourced reconstruction is produced.
2. GM1 has an explicit public CompOSE statement that the Douchin-Haensel crust is added below \(10^{-3}\,\mathrm{fm}^{-3}\), and that composition is available only for the core. This is a model-specific source boundary, not a rule to copy to other models. The source \`eos.compo\` contains distinct \`Iphase\` codes; retain the model dictionary and row masks.
3. KDE0v1 and SLY4 are unified tables. Their low-density cluster matter is part of the input and must not be called homogeneous core. For each model, extract the phase dictionary from that exact EOS page/data sheet and locate the last inhomogeneous/cluster state before admitting the homogeneous branch. The numerical \`Iphase\` codes seen in one model must not be transferred to the other.

The legacy \`n_B >= 0.045 fm^-3\` fallback is rejected as a physical core boundary. It may be retained only as an explicitly named diagnostic/screening window and never as the authority for core-only construction.

### 2.2 Crust-source cohort

Every row below is the archived \`eos.nb\` support, \(n_{\min}=4.67962\times10^{-10}\,\mathrm{fm}^{-3}\). The per-model upper support and point counts are:

| Model | \(n_{B,\max}\) [fm\(^{-3}\)] | Points | Construction family |
|---|---:|---:|---|
| SPG_M1 | 0.07413102413 | 365 | BSk22 low-density segment + model-specific CLD |
| SPG_M2 | 0.08912509381 | 377 | BSk22 low-density segment + model-specific CLD |
| SPG_M3 | 0.07943282347 | 366 | BSk22 low-density segment + model-specific CLD |
| SPG_M4 | 0.09772372210 | 375 | BSk22 low-density segment + model-specific CLD |
| SPG_M5 | 0.07079457844 | 363 | BSk22 low-density segment + model-specific CLD |
| GPPVA_DD2 | 0.067000002 | 283 | BSk22 + DD2 Thomas-Fermi |
| GPPVA_DDME2 | 0.071999997 | 288 | BSk22 + DD-ME2 Thomas-Fermi |
| GPPVA_FSU2 | 0.054000001 | 269 | BSk22 + FSU2 Thomas-Fermi |
| GPPVA_FSU2H | 0.086999997 | 300 | BSk22 + FSU2H Thomas-Fermi |
| GPPVA_FSU2R | 0.082999997 | 296 | BSk22 + FSU2R Thomas-Fermi |
| GPPVA_NL3wrL55 | 0.082000002 | 298 | BSk22 + NL3-\(\omega\rho\), \(L=55\) MeV Thomas-Fermi |
| GPPVA_TM1e | 0.089000002 | 303 | BSk22 + TM1e Thomas-Fermi |
| GPPVA_TW | 0.074000001 | 289 | BSk22 + TW Thomas-Fermi |

CompOSE identifies 13 neutron-star-crust tables, split into five unified-model crust tables and eight Thomas-Fermi crust tables. The DD2, DD-ME2, FSU2H, NL3\(\omega\rho L55\), TM1e, SPG M3 and SPG M5 pages are individually linked in the component registry. The archived table endpoints above are measured from the actual source grid, not inferred from the model name.

For the SPG/GPPVA records, keep all source-supported outer-crust, neutron-drip and inner-crust rows. The BSk22 low-density source splice below \(2.0\times10^{-3}\,\mathrm{fm}^{-3}\) is **not** the neutron-drip density. Neutron drip is a separate physical transition within the low-density sequence and must be retained with its own source/model evidence. Preserve every source-defined cluster and pasta phase through the source-supported endpoint. Do not append a homogeneous high-density core tail to a crust-only artifact.

The GPPVA family is based on Grill et al.'s self-consistent Thomas-Fermi construction; its nonspherical pasta phases must remain represented where the model data sheet specifies them. The SPG family uses model-specific CLD reconstruction associated with Scurto, Pais and Gulminelli and its cited functional inputs. Neither family permits an unverified scalar remapping of the final \`eos.compo\` tuple.

## 3. Thermodynamic and CompOSE data contract

The canonical barotropic variable is the total energy per baryon, including the source's rest-mass convention:
\[
F(n_B)=\frac{\epsilon(n_B)}{n_B}.
\]
For the cold one-parameter beta-equilibrium slice,
\[
P=n_B\frac{dF}{d\ln n_B},\qquad
\mu_B=F+\frac{P}{n_B},
\]
and
\[
\frac{c_s^2}{c^2}
=\frac{dP/d\ln n_B}{n_B\mu_B}.
\]
With \(Q_1=P/n_B\), this becomes
\[
\frac{c_s^2}{c^2}
=\frac{Q_1+dQ_1/d\ln n_B}{F+Q_1}.
\]

CompOSE v3.01 columns must not be misread as direct dimensional quantities. In particular:
- \(Q_1=P/n_B\);
- \(Q_2=s/n_B\);
- \(Q_3=\mu_b/m_n-1\);
- \(Q_4=\mu_q/m_n\);
- \(Q_5=\mu_l/m_n\);
- \(Q_6=f/(n_Bm_n)-1\);
- \(Q_7=e/(n_Bm_n)-1\).

Thus reconstruct
\[
F=m_n(1+Q_6),\quad
\epsilon/n_B=m_n(1+Q_7),\quad
P=n_BQ_1,\quad
\mu_{B,\mathrm{src}}=m_n(1+Q_3).
\]
At \(T=0\), the free-energy and internal-energy channels must agree within the declared tolerance, so \(Q_6=Q_7\) is a test, not an instruction to overwrite one with the other. Independently compare source \(Q_1\) with \(dF/d\ln n_B\), and source \(Q_3\) with the Euler value \(F+P/n_B\). All differences remain visible as residuals; a failed identity is not repaired by replacing provider data.

The \(m_n\), \(m_p\) and \(I_l\) values are read from each source's \`eos.thermo\` header. A universal \(m_n=939\) MeV constant is forbidden: the archived APR, BBB and GM1 tables, for example, record \(m_n=939.56533\) MeV and \(m_p=938.272\) MeV, whereas other provider files may use a different rounded convention. Preserve the source values, units and lepton convention.

Select \(T=0\) by the actual \`eos.t\` value and record \`iT\). Record the exact \`eos.yq\` value and \`iYq\` used for the source beta-equilibrium table. A stored \(Y_q=0\) coordinate is not a species abundance and must not be silently interpreted as a zero proton fraction. The selected slice identity is part of each component's provenance.

Use phase-local interpolation and derivatives only. Never differentiate, smooth, fit, or interpolate across a first-order jump or an inhomogeneous-to-homogeneous phase boundary. A negative compressibility or \(c_s^2<0\) is a stability failure; \(c_s^2/c^2>1\) is a causality failure. Do not clip either result. The existing C1-C11 framework must report per-field/per-row applicability as CHECKED, PARTIAL, NOT_AVAILABLE or FAIL.

## 4. Composition, phases and sparse records

The source \`eos.compo\` record is keyed by \((i_T,i_{n_B},i_{Y_q})\) and includes an integer \`Iphase\`, \`Npairs\`, particle-index/fraction pairs, \`Nquad\` and nuclear-group tuples. The raw row, exact indices and all presence masks must survive normalization.

- A missing species or missing \`eos.compo\` file is \`NOT_AVAILABLE\`, not a zero abundance.
- An explicitly written zero is present and remains distinct from an omitted field.
- \`Iphase\` is categorical. CompOSE v3.01 states that the physical meaning of each numeric phase index is defined on the corresponding EOS page. There is no global map from integer values to "outer crust", "inner crust", "pasta" or "uniform core".
- The source nuclear-group fraction \(Y_I\), average \(A\), and average \(Z\) remain separate. Convert to \(Y_h=A_{\rm av}Y_I\) only after the provider/model-specific semantics and row applicability have been admitted.
- Preserve the exact \`Nadd\` arity and \`q_extra\` type. Unknown/categorical entries are byte-preserved or explicitly marked unsupported; they are never averaged as a scalar.

APR and BBB explicitly lack compositional information. GM1 has core-only composition coverage, so the crust rows are not zero-filled. KDE0v1 and SLY4 contain source composition records with phase codes that require their own data-sheet interpretation. SPG/GPPVA use the already resolved dual track: preserve the provider's raw records under the standard parser, but production composition reconstruction follows the cited CLD/Thomas-Fermi and model-parameter papers. Do not reinterpret the original provider tuple in place.

## 5. Immutable artifact lineage and component interface

The construction pipeline has five distinct artifact classes:

1. **Raw source:** original provider ZIP and exact files. Immutable, hashed and never overwritten.
2. **Normalized source:** lossless CompOSE parse with axes, row indices, masses, column definitions, phase labels and presence masks. No physical field may be silently changed.
3. **Reconstructed source:** a separately versioned derivative for literature-derived composition or missing source fields. Every reconstructed field has its own citation, model, algorithm version and uncertainty/status.
4. **Standalone component:** \`core_only/<model>\` or \`crust_only/<model>\`, with exact lower/upper support and selection rule.
5. **Matched composite:** a later artifact with both component parents, the matching method ID/version and one-sided endpoint states. This node does not create it.

Each output manifest records at least: model ID, component role, raw archive SHA-256, exact source-file SHA-256 manifest reference, source Git blob ID, CompOSE entry/data-sheet, primary-paper references, source \(n_B\) support, selected \(T\)/\(Y_q\) values and indices, \(m_n,m_p,I_l\), Q-column map, phase dictionary, composition masks, \(N_{\rm add}\) schema, admitted component support, boundary method/value, interpolation version, endpoint \(P,\epsilon,\mu_B\), endpoint slope, source gaps, unsupported fields, transformation lineage and output SHA-256.

No core component may be released while its model-specific homogeneous-core boundary remains unresolved. No crust component may omit source-supported outer crust, neutron drip, inner crust or source-labelled cluster/pasta phases. No component may extrapolate beyond source support. A model-specific source gap is recorded as a gap; it is not filled with a guessed polytrope or neighboring model data.

## 6. Required acceptance tests

The native builders, when implemented in their own DAG nodes, must satisfy all of these:

1. Parse each \`eos.nb\` grid and prove exact point count, finite positive density, strict monotonicity and source bounds.
2. Preserve \`eos.t\`, \`eos.yq\` and every \((i_T,i_{n_B},i_{Y_q})\) key. Select \(T=0\) by value, not by assuming a particular index.
3. Read \(m_n,m_p,I_l\) from each \`eos.thermo\` header and preserve each Q-column definition.
4. Check \(Q_6-Q_7\), \(Q_1-dF/d\ln n_B\), and source \(Q_3\) versus Euler \(\mu_B\) independently.
5. Run C1-C11 only where their source fields are applicable. Missing fields remain NOT_AVAILABLE/PARTIAL.
6. Admit a core boundary only with a model-specific cited criterion and a numeric source-supported locator, or fail closed.
7. Preserve every source-supported crust interval, including neutron drip; keep the BSk22 \(2\times10^{-3}\,\mathrm{fm}^{-3}\) splice distinct from neutron drip.
8. Evaluate derivatives within a single physical phase; report instability/causality failures without clipping.
9. Round-trip unknown \`Nadd\`/\`q_extra\` and phase fields losslessly, or explicitly reject unsupported production use.
10. Produce immutable, versioned output plus source/transform/output hashes and a per-model receipt.

No compilation, numerical builder test or 390-task batch is claimed by this design node.

## 7. Next DAG nodes

**Next: E390-01C, model-specific boundary and phase-dictionary admission.** Resolve the exact APR/SLy4 and BHF/SLy4 splice locators, confirm the GM1 \(10^{-3}\,\mathrm{fm}^{-3}\) source-row selection, and extract the KDE0v1/SLy4 model-specific phase dictionaries and homogeneous-core onset criteria. The same pass locks crust-side phase dictionaries/terminal support for all thirteen crust inputs. Any unresolved item remains blocked rather than receiving an approximate density cutoff.

After E390-01C:
- E390-03A implements the native standalone core-only constructor for the five approved source intervals.
- E390-03B implements the native standalone crust-only constructor for the thirteen source-supported crust tables.
- E390-03 ports the literature-backed composition solvers/reconstruction to Fortran.
- E390-04 and E390-05 then build the smooth matching methods and the distinct sharp-Maxwell method on top of validated components.

## 8. Primary source references

- CompOSE reference manual v3.01, Typel et al., EPJA (2022), [DOI](https://doi.org/10.1140/epja/s10050-022-00847-y), [arXiv:2203.03209](https://arxiv.org/abs/2203.03209). The manual defines Q1-Q7, optional composition data and model-specific \(I_{\rm phase}\) semantics.
- APR: Akmal, Pandharipande & Ravenhall, PRC 58, 1804 (1998), [arXiv:nucl-th/9804027](https://arxiv.org/abs/nucl-th/9804027); BPS 1971; Douchin & Haensel 2001.
- BHF-BBB2: Baldo, Bombaci & Burgio, A&A 328, 274 (1997), [arXiv:astro-ph/9707277](https://arxiv.org/abs/astro-ph/9707277); BPS 1971; Douchin & Haensel 2001.
- GM1: Glendenning & Moszkowski, PRL 67, 2414 (1991), plus the specific CompOSE GM(GM1) table page and its crust statement.
- RG family: Gulminelli & Raduta, PRC 92, 055803 (2015), [arXiv:1504.04493](https://arxiv.org/abs/1504.04493); SLy4 from Chabanat et al. (1998); Danielewicz & Lee (2009).
- SPG crust: Scurto, Pais & Gulminelli, PRD 109, 103015 (2024), [arXiv:2402.15548](https://arxiv.org/abs/2402.15548); Gögelein et al., PRC 77, 025802 (2008); BSk22/Pearson et al. (2018).
- GPPVA crust: Grill et al., PRC 90, 045803 (2014), [DOI](https://doi.org/10.1103/PhysRevC.90.045803); model-specific RMF papers for DD2, DD-ME2, FSU2/H/R, NL3\(\omega\rho L55\), TM1e and TW are listed in the source bibliography JSON.

The exact per-model data-sheet URLs, raw \(n_B\) grid endpoints, point counts and Git blob identities are maintained in \`data/eos390_component_source_contract.v1.json\`.