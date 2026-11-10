# One-page summary: unresolved composition convention (both sources)

Status: DRAFT_NOT_SENT (attachment for provider inquiries). Prepared 2026-11-10. Origin: NS_Cowling PR #4, stage C6B/C6C.

**Issue.** The CompOSE `eos.compo` files for SPG(M1) and GPPVA(FSU2) (crust products) were read with the CompOSE v3.01 record layout (Eq. 4.19, 18 fields, 4 particle pairs, 1 nuclear-group quadruple). The data sheets define Iphase 0 (outer crust) and 7 (inner crust) but do not define the fourth quadruple field of the nuclear group, nor whether Aav/Zav are cluster or cell quantities, nor how a free neutron gas is recorded. We want to understand the export convention *before* interpreting the data scientifically.

**Equation sheet (what we need to reconcile).**
- CompOSE: Y_i = n_i/n_b; X_i = B_i Y_i; for a nuclear group Y_I = sum_j Y_j and X_h = Aav * Y_I.
- Charge closure: r_C3 = Yp + (Zav/Aav) Yh - Ye - Ymu.
- Baryon closure: r_C4 = Yn + Yp + Yh - 1.
- Reading "Yh = raw fourth field" closes both for outer crust rows; reading "Yh = Aav * Y_I" closes none. In inner-crust rows neither reading closes both.

**Pinned sources** (authority commit `a09ec30d74b363652d2d1aefe990557c468679dd`, path `eos_generator/eos_data/<SRC>/eos.compo`)

| Source | Rows | git blob SHA-1 | SHA-256 |
|---|---:|---|---|
| SPG_M1 | 365 | `ab79da517271bea1f7de37294e87e019de611699` | `8ce704b9af8a3a7a74e6985d9edc5481a8dc7a8b03e50f7b4d39c3d292dfec27` |
| GPPVA_FSU2 | 269 | `fa1ad2b809eb1c95825e2ee6eb2650cf3f248693` | `0250865a4848a855727b18b6d312450a2d8e3d972bb52a811fe14fee5adfc5d8` |

Row n = n-th line of the file; `inb` = n.

**Aggregates (full tables in `aggregate_statistics.json`; every value derived from the pinned rows).**

| | SPG_M1 Iphase=7 (rows 218-365, 148) | GPPVA_FSU2 Iphase=7 (rows 218-269, 52) |
|---|---|---|
| Yn = Yp = 0 | 148/148 | 52/52 |
| Fourth field = Ye + Ymu (1e-9) | 148/148 | 52/52 |
| Aav, Zav | 82.107-452.669, 21.149-27.512 (non-integer) | 740-2143, 13-36 (integers) |
| Zav/Aav vs Ye | differs in 148/148 (0.061-0.297) | equals Ye in 52/52 |
| Yh implied by C3 | 0.0960-0.5433 | 0.99963-1.00038 |
| Yh implied by C4 | 1.0 | 1.0 |
| C3 vs C4 implied-Yh gap, max | 0.904 (inconsistent 148/148) | 3.8e-4 (consistent 52/52) |
| Baryons unaccounted if C3-implied Yh | 0.457-0.904 | within +-4e-4 |
| Max abs residual with raw field as Yh (C3, C4) | 0.0428, 0.978 | 0.0463, 0.989 |

Outer crust (rows 1-217) is identical in both files: fourth field 1.0, Yn=Yp=Ymu=0, Ye=Zav/Aav, both closures satisfied with Yh = raw field (max |r| 5e-6).

GPPVA only: `eos.thermo` carries an extra geometry column (0 for Iphase=0, 1 for all Iphase=7 rows; "droplets" per data sheet), distinct from Iphase.

**Related public context.** The GPPVA(FSU2) unified product (compose.obspm.fr/eos/215, eos.compo sha256 `a57a6f58e9080c2f501cd5b4d9f9f45ded0833df27eb107bacf36fe55dbf5df1`) repeats rows 1-269 identically; its 157 core rows (Iphase=6) have Aav=Zav=Y_I=0 with Yn, Yp nonzero.

**Links.** CompOSE manual v3.01 https://doi.org/10.1140/epja/s10050-022-00847-y (arXiv:2203.03209). SPG data sheet https://compose.obspm.fr/download/1D/Crust/M1/eos.pdf. GPPVA data sheet https://compose.obspm.fr/download/1D/Crust/FSU2/eos.pdf; model page https://compose.obspm.fr/eos/209. Papers: Scurto, Pais, Gulminelli PRD 109, 103015 (2024), arXiv:2402.15548; Grill et al. PRC 90, 045803 (2014), arXiv:1404.2753; Providência et al. Front. Astron. Space Sci. 6, 13 (2019), arXiv:1811.00786.

**What only the provider can resolve:** the meaning of the fourth quadruple field and group index for Iphase 0/7; whether Aav/Zav describe the cluster or the whole cell; where (if anywhere) the SPG gas neutrons are stored; whether Yn=Yp=0 is deliberate; the normalization needed to reconstruct total baryon content. No numerical mapping is proposed or implied by this note.
