# C6B provider provenance (retrieved 2026-11-10)

Evidence classes: **[SPEC]** official CompOSE specification; **[DS]** provider data sheet; **[PRIMARY]** primary paper matched to the exact model; **[RELATED]** related, provenance not confirmed; **[UNRESOLVED]**.

## Standard definitions [SPEC]
CompOSE Reference Manual v3.01 (arXiv:2203.03209; PDF header reads "Version 3.01"; the EPJA DOI 10.1140/epja/s10050-022-00847-y was supplied in the handoff and not independently resolved here). Local copy sha256 `2f2248cd193be578e642b844d50bc604fe8fe637b4aa37a46f331063d1a2b485`.
- Eq. 3.10 `Y_i = n_i/n_b`; Eq. 3.11 `sum B_i Y_i = 1`; Eq. 3.12 `X_i = B_i Y_i`.
- Sec. 4.2.3, record layout Eq. 4.19: `iT inb iYq Iphase Npairs (I Y)*Npairs Nquad (I Aav Zav Y_I)*Nquad`. Raw rows have 18 fields = 3+1+1+8+1+4, matching exactly.
- Eqs. 4.20-4.22: `Aav = sum A_j Y_j / sum Y_j`, `Zav` likewise, `Y_I = sum A_j Y_j / Aav = sum_j Y_j` (combined **number** fraction). So `X_h = Aav*Y_I`.
- Sec. 4.2.3: the meaning of Iphase and of the group index is defined per model in the data sheet; the contributor "needs to define the meaning of the index Iphase".
- Sec. 4.2.8/4.2.1: `eos.t` row 1 = N_Tmin, row 2 = N_Tmax, then T values. SPG_M1 (0,0) and GPPVA_FSU2 (1,1) match their `iT` keys.
- Table 3.x particle indices 10/11/0/1 (neutron/proton/electron/muon) were used as in the C5 reader; the table text was not re-extracted here (**[UNRESOLVED]** only in that sense).

## SPG(M1)
- Data sheet **[DS]**, sha256 `491c23709054bf8af3fa78b5f0563fbb5b1c3bc073dfea1747d3f6daa5e488d2` (https://compose.obspm.fr/download/1D/Crust/M1/eos.pdf; sheet date 2024-02-26): BSk22 outer crust below 2e-3 fm^-3, compressible-liquid-drop inner crust, GDFM density-dependent RMF, 4 particle pairs + 1 heavy-nucleus quadruple, Iphase 0 = outer crust, 7 = inner crust. It gives **no** definition of the quadruple's last field, of the group index, or of how a free neutron gas is recorded.
- **[PRIMARY, provenance of the model]** Scurto, Pais, Gulminelli, Phys. Rev. D 109, 103015 (2024), arXiv:2402.15548 (local pdf sha256 `92b1fa0318bd140abf459bfc654cae236c17d503b4af166e64de1e764c081a7a`). Its "Publicly Available Models" section selects five models M1-M5 for CompOSE; M1 = highest absolute likelihood. Its Table IX NMP values contain J=33.9, Ksym=-181.4 and Mmax=2.54, matching the data sheet's M1 entries. In its CLD model each Wigner-Seitz cell has a cluster part I and a **gas part II** (Eq. 21). The paper does **not** define the eos.compo fields (Aav, Zav, Y_I, Yn, Yp), and I found no composition table to compare against.
- [RELATED, not the table's provenance] Scurto, Pais, Gulminelli, PRC 107, 045806 (2023), arXiv:2212.09355 (sha256 `ef96dee0311986a2d78bc207baa3c4f79c56a2d6c1dd98fbfe228a762974d9e6`): CLD pasta with NL3/NL3wrho, not GDFM. Data sheet's "[4] in preparation" is consistent with the PRD paper but that identification is an inference.
- Outer crust: Pearson et al., MNRAS 481, 2994 (2018) (BSk22), https://compose.obspm.fr/eos/210 (not fetched here).

## GPPVA(FSU2)
- Data sheet **[DS]**, sha256 `660554a75fdeab565c447adc662210c56d2f42b5ca1d1defd12a9f457259475d` (https://compose.obspm.fr/download/1D/Crust/FSU2/eos.pdf; sheet date 2023-04-26; model entry https://compose.obspm.fr/eos/209 not separately fetched): Thomas-Fermi inner crust incl. pasta, FSU2, BSk22 outer crust below 2e-3 fm^-3, Iphase 0/7/6, extra thermo column = pasta geometry (0 outer crust, 1 droplets, ...). No definition of the quadruple's last field.
- **[PRIMARY, method]** Grill, Pais, Providencia, Vidana, Avancini, PRC 90, 045803 (2014), arXiv:1404.2753 (sha256 `2fee7d014cb2077f6fba0750730da025bab55363ca8e3c1508d0ccaebde993e0`): Thomas-Fermi WS cells; Table IV gives only the EoS (energy density, pressure). No A, Z, Y_p composition table.
- **[PRIMARY, publication of this EoS per data sheet]** Providencia, Fortin, Pais, Rabhi, Front. Astron. Space Sci. 6, 13 (2019), arXiv:1811.00786 (sha256 `ecef4d252c1357d99f436546fe8429e351f5ff1f6d36116771ce38aaa231bb43`): supplementary material lists only baryon density, energy density, pressure. No composition.
- [RELATED] Chen & Piekarewicz, PRC 90, 044305 (2014) (FSU2 parameter set; not fetched).
- In the raw table every Iphase=7 row has extra geometry code 1 (droplets); every Iphase=0 row has code 0. The geometry code is distinct from Iphase.

## Conclusion on public evidence
Neither source publication nor either data sheet defines the nuclear-group last field, the group index, or whether Aav/Zav are cluster or cell values. No published composition values exist to compare with the rows. **[UNRESOLVED]**: needs the providers.
