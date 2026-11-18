# Maxwell ∩ CSS with cubic/quintic Hermite: source-origin independent census

## Source-protected result (2026-11-18)

**95 distinct source records, 635 same-origin pairing candidates**.
Tier mapping: **original** = original provider CompOSE plus original nuda raw crust; **corrected** = energy-primary public thermodynamic reconstructions; **reconstructed** = independently generated local code outputs.
A unique role-qualified source ID is retained for each; there is NO original × corrected or provider × reconstructed pairing.
All 635 candidate rows are available in ALL_CANDIDATE_PAIRS.csv, prioritized by independent crust.
Native finite Maxwell and native CSS both succeeded on **7** pairings. Matched output physical claims are scoped to their actually valid finite source ranges.
**Strict unconditional full inner+outer crust EOS passing both: 0.**
**Conditional code-completed inner+outer stellar EOS passing both and M=1.4M⊙ TOV: 5.**

**Matching 방법은 Maxwell와 CSS만 비교한다.** Cubic/quintic은 EOS 입력자료를 P(mu)로 표현하는 보간 차수로만 사용한다. 보간법 자체가 세 번째 matching construction이 아니다. 검증된 7조합의 2×2 자료는 HERMITED_INPUT_2X2/CAMPAIGN.json과 DUAL_FINITE_SUCCESS.csv에서 확인할 수 있다.

## All finite dual-matching successes (7 pairs)

| Source tier | Crust | Core | Same-functional core exists? | Maxwell P_t [MeV fm⁻³] | CSS P_t [MeV fm⁻³] | Cubic+Maxwell / Quintic+Maxwell | Cubic+CSS / Quintic+CSS | Whole status |
|---|---|---|---|---:|---:|---|---|---|
| reconstructed | GPPVA_FSU2 | GM_GM1 | GPPVA(FSU2) unified inner crust-core EOS; corresponding same-RMF-functional core reference (not used unless explicitly registered) | 0.0777479 | 0.252729 | PASS_NATIVE_MAXWELL_MATCHED / PASS_NATIVE_MAXWELL_MATCHED | PASS_DERIVED_INTERPOLATED_SOURCE / PASS_DERIVED_INTERPOLATED_SOURCE | OUTER_ABSENT |
| reconstructed | SPG_M1 | PK1_CANONICAL | SPG(M1) unified NS EOS; corresponding same-GDFM-functional core exists in CompOSE (not used for the matched candidate) | 0.0174442 | 0.512779 | PASS_NATIVE_MAXWELL_MATCHED / PASS_NATIVE_MAXWELL_MATCHED | PASS_DERIVED_INTERPOLATED_SOURCE / PASS_DERIVED_INTERPOLATED_SOURCE | conditional BSk22+TOV |
| reconstructed | SPG_M2 | PK1_CANONICAL | SPG(M2) unified NS EOS; corresponding same-GDFM-functional core exists in CompOSE (not used for the matched candidate) | 0.0143074 | 0.445755 | PASS_NATIVE_MAXWELL_MATCHED / PASS_NATIVE_MAXWELL_MATCHED | PASS_DERIVED_INTERPOLATED_SOURCE / PASS_DERIVED_INTERPOLATED_SOURCE | conditional BSk22+TOV |
| reconstructed | SPG_M3 | PK1_CANONICAL | SPG(M3) unified NS EOS; corresponding same-GDFM-functional core exists in CompOSE (not used for the matched candidate) | 0.0185507 | 0.565705 | PASS_NATIVE_MAXWELL_MATCHED / PASS_NATIVE_MAXWELL_MATCHED | PASS_DERIVED_INTERPOLATED_SOURCE / PASS_DERIVED_INTERPOLATED_SOURCE | conditional BSk22+TOV |
| reconstructed | SPG_M4 | PK1_CANONICAL | SPG(M4) unified NS EOS; corresponding same-GDFM-functional core exists in CompOSE (not used for the matched candidate) | 0.0178707 | 0.580962 | PASS_NATIVE_MAXWELL_MATCHED / PASS_NATIVE_MAXWELL_MATCHED | PASS_DERIVED_INTERPOLATED_SOURCE / PASS_DERIVED_INTERPOLATED_SOURCE | conditional BSk22+TOV |
| reconstructed | SPG_M5 | PK1_CANONICAL | SPG(M5) unified NS EOS; corresponding same-GDFM-functional core exists in CompOSE (not used for the matched candidate) | 0.0187242 | 0.392475 | PASS_NATIVE_MAXWELL_MATCHED / PASS_NATIVE_MAXWELL_MATCHED | PASS_DERIVED_INTERPOLATED_SOURCE / PASS_DERIVED_INTERPOLATED_SOURCE | conditional BSk22+TOV |
| original | XMLSLZ_TM1 | GM_GM1 | XMLSLZ(TM1) original same-provider uniform phase-0 core exists; source-phase controls are distinct from cross-model pairs | 0.140941 | 0.325076 | PASS_NATIVE_MAXWELL_MATCHED / PASS_NATIVE_MAXWELL_MATCHED | PASS_DERIVED_INTERPOLATED_SOURCE / PASS_DERIVED_INTERPOLATED_SOURCE | SOURCE_INNER_GAPS_OR_EXTRA_PHASES_UNVALIDATED |

**The dual-success set has no pair sharing the same nuclear functional**. The 'same-functional core exists' column identifies the external/unified functional controls, not identity with the actual matched core.

**2×2 source-grounded comparisons:** 7 native admissible finite input pairs × (cubic/quintic P(mu) source representation) × (Maxwell/CSS matching) = 28 native matching evaluations. Maxwell was independently re-run by the Fortran transition generator on phase-restricted, Hermite-interpolated inputs (14/14 passed). CSS was likewise re-run with the native CSS generator (14/14 passed). The earlier root-only local Hermite scan remains an independent diagnostic. Original source knots were retained, readback was checked. These results do not prove complete unified microscopic EOSs.

### Exactly two matching methods × two EOS interpolation orders

| Crust/core | Cubic+Maxwell P_t | Quintic+Maxwell P_t | Cubic+CSS P_t | Quintic+CSS P_t |
|---|---:|---:|---:|---:|
| reconstructed / GPPVA_FSU2+GM_GM1 | 0.0777478555 | 0.0777478885 | 0.2527287 | 0.2527287 |
| reconstructed / SPG_M1+PK1_CANONICAL | 0.0174442116 | 0.017444211 | 0.512778572 | 0.512778572 |
| reconstructed / SPG_M2+PK1_CANONICAL | 0.0143073983 | 0.0143073979 | 0.445754543 | 0.445754543 |
| reconstructed / SPG_M3+PK1_CANONICAL | 0.0185507404 | 0.01855074 | 0.565705225 | 0.565705225 |
| reconstructed / SPG_M4+PK1_CANONICAL | 0.0178706915 | 0.0178706905 | 0.580962229 | 0.580962229 |
| reconstructed / SPG_M5+PK1_CANONICAL | 0.0187241843 | 0.0187241837 | 0.392475445 | 0.392475445 |
| original / XMLSLZ_TM1+GM_GM1 | 0.140941859 | 0.140942021 | 0.3250764 | 0.3250764 |

Note: the two Maxwell values in each row are from native Fortran matched EOSs on Hermite-interpolated source representations, not merely a sign-change root. No cubic_F or quintic_F smooth connector is used as a matching method. Conditional TOV below was performed on source-matched constructions with E outer closure, not on all of the newly interpolated 2×2 variants.

## Conditional whole-star verification (not a unified-microscopic claim)

| Reconstructed crust+core | Maxwell R_1.4 [km] | CSS R_1.4 [km] | Same-source E outer and causal lower join? |
|---|---:|---:|---|
| SPG_M1 + PK1_CANONICAL | 14.4204291 | 14.4158109 | Yes, explicit conditional closure |
| SPG_M2 + PK1_CANONICAL | 14.4204219 | 14.4204276 | Yes, explicit conditional closure |
| SPG_M3 + PK1_CANONICAL | 14.4204326 | 14.4156879 | Yes, explicit conditional closure |
| SPG_M4 + PK1_CANONICAL | 14.4204327 | 14.4259301 | Yes, explicit conditional closure |
| SPG_M5 + PK1_CANONICAL | 14.4204328 | 14.4179570 | Yes, explicit conditional closure |

The five native SPG crust sources lack an outer crust. For these conditional whole models we retained their computed inner crust and appended the *previously reconstructed* BSk22 outer/anchored surface plus an explicit causal low-density connector, without altering source inner knots. Their mass was actually integrated to M≈1.4M⊙. This is not proof that the entire star is generated under one microscopic functional and not proof of global TF minimality.
GPPVA_FSU2+GM1 failed the full-star consumer since inherited source phase boundaries cannot be smoothed across; it is retained as a finite-interface success only.
Original provider XMLSLZ(TM1)+GM1 has an outer source but an invalid/unadmitted hole between accepted crust source components, thus full crust is not approved. It is a **finite Maxwell+CSS success, not a whole-neutron-star success**.

## Prioritized independent crust candidates and constraints

The inventory covers **70 crust source entries** with same-functional reference and original low-density support flags.
External Negele–Vautherin (1973) raw table contains baryon density and E/A-m_n but not pressure P or a verified full EOS; original-table Maxwell cannot be calculated from it without explicitly reconstructing pressure and fixing integration constants. It stays BLOCKED, not a negative physics result.
External Gogny D1S/D1M/D1MS and GMRS 2022 contain source n,P,E/A-m_n but do not have a complete verified same-source outer crust down to the surface nor independently verified total-energy μ and exact CompOSE slice. The original raw corpus is preserved. The external paired EOS is **not fabricated or counted**.
The external BPS crust historical reference is outer-only, BBP is inner-only; without a full matched physics and normalization contract neither alone satisfies the user-mandated inner+outer condition.
Previously user-uploaded MYN13_RHF_TF and OMEG1_EoS are inventoried as original ZIP-only sources with exact source-file SHA. Earlier direct ZIP audit found four mathematical crossings but no physically adopted crust/core Maxwell root. They have no sealed original CompOSE source slice on the remote workstation; therefore new cross-combinations are explicitly blocked, not called no-coexistence.
All 13 provider-original SPG/GPPVA crust tables contain an attached BSk22 outer segment and source inner crust, although native original-source admissible inner components can have holes. Their corrected counterparts retained unchanged outer rows but have energy-primary reconstructed inner pressure; source selection/first-law assumptions differ and are never combined with originals.
The five XMLSLZ original provider crust records are kept as original tier, with their own support provenance. Matched cores are original source of the same tier only.

## Scientific decision rules and limitations

1. **오직 Maxwell와 CSS만 matching construction으로 사용한다.** Cubic/quintic Hermite는 이 두 matching에 공급하는 EOS의 P(mu) 입력 보간법이다. 7개 유한 성공 조합에서 (cubic,quintic)×(Maxwell,CSS)의 28개 조건을 별도로 계산했다. CSS는 Hermite spline 자체가 아니다.
2. A finite Maxwell root + CSS PASS is not a complete NS EOS. Requiring an approved outer crust, complete inner phase chain and 1.4M⊙ TOV shrinks the viable set from 7 finite to 5 *conditional* and 0 strict source-native whole-star results.
3. Generic origin-tier mixing prohibited. A newly created thermodynamic interpolation from an original source is a new reconstructed object, never an original-provider value. No pressure/root location/extrapolation is fabricated to rescue a rejected combination.
4. The coarse cubic P(μ) screen identifies positive roots; for new unscreened pairings a screen-negative result is **not an exhaustive nonlinear root nonexistence proof** (tangential root, multiple narrow roots, invalid source branches). An old rigorous native Maxwell 195 pair campaign was reused with exact source identity.
5. Physical assumption: single conserved baryon charge and local charge neutrality, beta-equilibrated cold barotropes; additional pasta phase competition, formation history, surface tension/interface tension, superfluidity, frozen composition, radial stability and full-GR oscillations remain unresolved.

## Evidence and reproducibility

- data/all_origin_maxwell_css_20261118/SOURCE_CATALOG.json (source/provenance/input support tier)
- data/all_origin_maxwell_css_20261118/CAMPAIGN.json (635 pair results and actual CSS/new Maxwell receipts)
- data/all_origin_maxwell_css_20261118/whole_stars/CAMPAIGN.json and */RESULT.json (6 bounded attempts, 5 conditional dual TOV)
- data/maxwell_three_families_20261117/campaign_v2/CAMPAIGN.json (195 prior independent native Maxwell pairings)
- nucleardatapy 1.0.2 original external distribution wheel SHA recorded in SOURCE_CATALOG, downloaded separately under scratch, not altered or silently republished

## Next scientific boundary
To qualify the full original/corrected crust models for a strict whole-star Maxwell+CSS study, source gaps must be resolved by author-supplied physically consistent low-density tables or independently solved source EOS with identical energy/charge conventions. The 5 conditional SPG+PK1 stars are usable phenomenological comparison baselines, not a validation of a shared microscopic crust-core functional.
