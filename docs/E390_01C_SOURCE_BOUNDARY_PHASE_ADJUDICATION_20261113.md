# E390-01C: Source-exact component boundary and phase adjudication

Date: 2026-11-13  
Repository: `cosmosapjw-quantum/NS_Cowling`  
Execution branch: `cloud/eos390-fortran-unification-v1`  
Pre-node HEAD: `357d9a51d65c41cee45dd987c84a346ef0cc3370`  
Pre-node backup: `backup/e390-01c-pre-20261113`  
DAG: `data/eos390_native_fortran_unification_dag.v1.json`, revision 1.4.7 before this update

## Result

**Node outcome: source adjudication complete; GPPVA_TW crust admission remains blocked by an official-source conflict.** All 5 core inputs and 13 crust inputs have explicit table-row boundaries/support summaries. Five core models are admitted for the next core-only constructor. Twelve crust models are admitted for future crust-only construction; GPPVA_TW is explicitly excluded until its outer-crust provider is resolved.

This is a source-adjudication receipt only. No production code, source EOS table, reconstructed EOS, standalone component, or matched EOS was modified or generated.

## Validation of source grids

The 18 tracked `eos.nb` files were read from the execution branch. Each header point count agrees with the number of parsed density entries. All 18 density grids are positive and strictly increasing. Their exact Git blob SHAs and model endpoints are recorded in the companion JSON receipt.

The source archive's historical E390-00 delivery-manifest SHA-256 is `d671310eb630ffd13f1b278425d945741ad7e08ce14c0d9aff42309c707e6a3b`. That prior receipt records archive CRC and inventory checks. The archive bytes were not independently downloaded or rehashed in this session. This turn instead verifies the individual tracked Git blobs and source row contents.

## Core-only source boundaries

| Core source | Model-specific construction evidence | Source rows to retain as homogeneous core | Boundary evidence |
|---|---|---:|---|
| APR_APR | CompOSE APR page 68 identifies APR core, SLy4 inner crust and BPS outer crust | 70–172 | Rows 1–69 are the pre-attached crust segment; row 69 is nB=0.075959 fm^-3; row 70 is nB=0.085 fm^-3 |
| BBB_BHF2 | CompOSE BBB page 69 identifies BHF-BBB2 core, SLy4 inner crust and BPS outer crust | 70–84 | Rows 1–69 are the pre-attached crust segment; row 69 is nB=0.075959 fm^-3; row 70 is nB=0.16 fm^-3 |
| GM_GM1 | CompOSE page 54 states Douchin–Haensel crust below nB=1.0e-3 fm^-3 and composition is available for core only | 49–349 | eos.compo Iphase changes 1→3 at row 49; row 48 nB=0.00025058488505, row 49 nB=0.048; threshold 0.001 is not a grid member |
| RG_KDE0v1 | Its own CompOSE sheet defines Iphase 4 = heavy nuclei present and 3 = homogeneous matter | 967–1240 | Iphase 4 ends at row 966, nB=0.0414954; Iphase 3 starts row 967, nB=0.04582369 |
| RG_SLY4 | Its own CompOSE sheet defines Iphase 4 = heavy nuclei present and 3 = homogeneous matter | 144–197 | Iphase 4 ends at row 143, nB=0.04786301; Iphase 3 starts row 144, nB=0.05248075 |

APR and BBB share the same first 69 density-grid points, but not an identical thermodynamic table. The first 28 thermodynamic rows are identical; the first difference is at row 29. The source grids then diverge at row 70: APR continues at 0.085 fm^-3 while BBB continues at 0.16 fm^-3. Douchin–Haensel (2001) gives the SLy crust edge near 0.076 fm^-3 and a core-side coexistence density about 1.4% higher. The source gap is retained; no density interpolation, extrapolation, or invented matching point is inserted.

For GM1, the live CompOSE page and tracked 349-point `eos.nb` agree. Its linked 2014 PDF data sheet reports 247 points, a metadata discrepancy recorded here rather than used to override the live page and source blob. The documented 0.001 fm^-3 cutoff lies inside the large source-grid interval between rows 48 and 49; component membership is therefore taken from the core-only composition support and model statement, not by manufacturing a grid point at the cutoff.

## Neutron-drip and 13 crust-only inputs

The official CompOSE BSk22 composition table identifies:
- Iphase 1, rows 1–217: outer crust; row 217 has no free-neutron fraction.
- Iphase 2, rows 218–300: inner crust; the free-neutron fraction first appears at row 218 (Y_n=0.021063081).
- Iphase 0, rows 301–489: homogeneous core.

Pearson et al. (2018), Table 7, reports the BSk22 neutron-drip density as approximately 2.69e-4 fm^-3. The archived SPG/GPPVA crust tables contain the same 217-point low-density outer segment, ending at 0.000269021 fm^-3 (row 217). Their source splice to the model-specific inner-crust calculation is separately documented below 0.002 fm^-3. The 0.002 fm^-3 source splice is not neutron drip. The first stored inner-crust row is model-grid dependent and can be above 0.002; no artificial 0.002 row is inserted.

| Crust input | Rows / points | Iphase source rows | First inner row nB (fm^-3) | Terminal row nB (fm^-3) | Outer provider / method |
|---|---:|---|---:|---:|---|
| SPG_M1 | 365 | 0: 1–217; 7: 218–365 | 0.002511886432 | 0.07413102413 | BSk22 / CLD |
| SPG_M2 | 377 | 0: 1–217; 7: 218–377 | 0.002290867653 | 0.08912509381 | BSk22 / CLD |
| SPG_M3 | 366 | 0: 1–217; 7: 218–366 | 0.002630267992 | 0.07943282347 | BSk22 / CLD |
| SPG_M4 | 375 | 0: 1–217; 7: 218–375 | 0.002630267992 | 0.09772372210 | BSk22 / CLD |
| SPG_M5 | 363 | 0: 1–217; 7: 218–363 | 0.002511886432 | 0.07079457844 | BSk22 / CLD |
| GPPVA_DD2 | 283 | 0: 1–217; 7: 218–283 | 0.0020000001 | 0.067000002 | BSk22 / Thomas–Fermi |
| GPPVA_DDME2 | 288 | 0: 1–217; 7: 218–288 | 0.0020000001 | 0.071999997 | BSk22 / Thomas–Fermi |
| GPPVA_FSU2 | 269 | 0: 1–217; 7: 218–269 | 0.003 | 0.054000001 | BSk22 / Thomas–Fermi |
| GPPVA_FSU2H | 300 | 0: 1–217; 7: 218–300 | 0.0020000001 | 0.086999997 | BSk22 / Thomas–Fermi |
| GPPVA_FSU2R | 296 | 0: 1–217; 7: 218–296 | 0.0020000001 | 0.082999997 | BSk22 / Thomas–Fermi |
| GPPVA_NL3wrL55 | 298 | 0: 1–217; 7: 218–298 | 0.0020000001 | 0.082000002 | BSk22 / Thomas–Fermi |
| GPPVA_TM1e | 303 | 0: 1–217; 7: 218–303 | 0.0020000001 | 0.089000002 | BSk22 / Thomas–Fermi |
| GPPVA_TW | 289 | 0: 1–217; 7: 218–289 | 0.003 | 0.074000001 | **BLOCKED:** live page says BPS; linked PDF says BSk22 |

### Separate GPPVA geometry field

For GPPVA, `eos.thermo` contains a last-column categorical geometry field, distinct from `eos.compo:Iphase`. The provider data sheets define: 0 outer crust, 1 droplets, 2 rods, 3 slabs, 4 tubes, 5 bubbles, 6 homogeneous core. The 8 archived crust tables currently use only codes 0–3. This field is categorical and must be copied losslessly; it must never be numerically interpolated. It is not interchangeable with `Iphase=0/7/6`.

### Provider-source exceptions

- GPPVA_FSU2R: current official CompOSE page 203 explicitly names BSk22 below 0.002 fm^-3; the expected PDF endpoint returned 404, so the provider statement is page-backed, while the phase legend is corroborated by other GPPVA data sheets and the observed source fields.
- GPPVA_TW: official live page 212 says BPS below 0.002 fm^-3, whereas the directly linked official PDF says BSk22. Although its low-density row grid coincides with the common 217-row outer segment in the other crust tables, that is not enough to choose one provider. TW remains `BLOCKED_COMPONENT_ADMISSION` pending a corrected data sheet/provider confirmation.

## Admission and next node

- Core-only admission: 5/5 source boundaries adjudicated. Proceed to E390-03A.
- Crust-only admission: 12/13 source boundaries adjudicated. E390-03B must include only the 12 admitted inputs and fail closed for GPPVA_TW until its provider identity is resolved.
- The constructors must preserve all source rows within their admitted ranges, keep thermodynamic/composition presence masks and provenance, and avoid extrapolation across source gaps.
- No core-only EOS, crust-only EOS, composite EOS, matching method, or 390-pair production output is claimed by this node.

## Primary source links

- CompOSE APR: https://compose.obspm.fr/eos/68
- CompOSE BBB/BHF-BBB2: https://compose.obspm.fr/eos/69
- CompOSE GM1: https://compose.obspm.fr/eos/54
- CompOSE KDE0v1 sheet: https://compose.obspm.fr/download//1D/NS/Skyrme/KDE0v1/eos.pdf
- CompOSE SLy4 sheet: https://compose.obspm.fr/download//1D/NS/Skyrme/SLY4/eos.pdf
- CompOSE BSk22: https://compose.obspm.fr/eos/256
- CompOSE SPG M1–M5 sheets: https://compose.obspm.fr/download/1D/Crust/M1/eos.pdf through M5
- CompOSE GPPVA DD2/DDME2/FSU2/FSU2H/NL3wr/TM1e/TW sheets: https://compose.obspm.fr/download/1D/Crust/<model>/eos.pdf
- CompOSE GPPVA FSU2R: https://compose.obspm.fr/eos/203
- Douchin & Haensel (2001): https://arxiv.org/abs/astro-ph/0111092
- Pearson et al. (2018), MNRAS 481, 2994: https://academic.oup.com/mnras/article/481/3/2994/5090419
- CompOSE Reference Manual v3.01: https://arxiv.org/abs/2203.03209

