# C6B hypothesis adjudication (per source)

All numbers are from `row_audit_*.csv` / `closure_summary.json` (independent parser; C6 residual maxima reproduced exactly: SPG_M1 C3 0.042790704812013594, C4 0.9779100000000001; GPPVA_FSU2 C3 0.04628324324324324, C4 0.98856).

## Facts established from the raw rows
1. Grammar and keys (H4): both files parse with the manual's Eq. 4.19 layout (18 fields/row, Npairs=4, Nquad=1, particle order 10 11 0 1, group id 1); (iT,inb,iYq) equal the eos.thermo keys row by row; inb=1..N consecutive; N equals the eos.nb count; eos.nb strictly increasing; iT=0 (SPG) / 1 (GPPVA) matches eos.t N_Tmin.
2. Outer crust (Iphase=0, rows 1-217): rows are identical in the two sources (shared BSk22 table). Last field = 1.0, Yn=Yp=Ymu=0, Ye = Z/A. Direct Yh=Y_I closes C3/C4; Aav*Y_I (e.g. 56, 128) fails every row. Under Eq. 4.22 a group with Aav=56 holding all baryons would have Y_I = 1/56, so the stored 1.0 behaves as the baryon fraction X_h, not as Y_I.
3. Inner crust (Iphase=7): the last field equals Ye+Ymu to print precision in **200 of 200** rows (SPG 148, GPPVA 52). It is not a group number fraction (that would be ~1/A ~ 1e-3). Yn=Yp=0 in all rows of both tables.
4. GPPVA inner: Aav, Zav are integers; |Ye - Z/A| <= 5e-6; C3-implied and C4-implied Yh agree in 52/52 rows (0.9996-1.0004, gap <= 3.8e-4, consistent with 3-5 digit rounding). The record is self-consistent with X_h = 1 (A, Z as whole-cell totals), but the stored field is not 1.
5. SPG inner: Aav, Zav non-integer; Z/A != Ye. C3-implied Yh in [0.096, 0.543], C4-implied Yh = 1: incompatible in 148/148 rows (gap up to 0.904). The implied non-cluster baryon fraction is 0.457-0.904 and no Yn is recorded.

## Hypotheses
| | SPG_M1 | GPPVA_FSU2 |
|---|---|---|
| H4 parser/index defect | REJECTED (grammar, keys, counts verified) | REJECTED |
| Standard Y_I conformance | NO: outer field behaves as X_h; inner field equals Yq | NO (same) |
| H1 documented provider convention | NOT SUPPORTED: no document defines one | NOT SUPPORTED |
| H2 mis-populated field | SUPPORTED numerically (field == Yq on all inner rows); not provider-confirmed | SUPPORTED numerically; not provider-confirmed |
| H3 incomplete group representation | SUPPORTED numerically (non-cluster baryons 0.46-0.90 unrecorded; paper's CLD cell has a gas part) ; not provider-confirmed | NOT NEEDED: closure holds with cell-total A,Z, X_h=1 |
| H5 verifier contract | Verifier equations are the standard closure; no scalar from the raw field can satisfy both C3 and C4 | Equations close only with Yh=1, which is not obtainable from the raw field |
| Final | PROVIDER_SEMANTICS_UNRESOLVED__BLOCKED | PROVIDER_SEMANTICS_UNRESOLVED__BLOCKED |

## Rejected as candidate global mappings
- Yh = Aav*Y_I: fails C3/C4 on all 634 rows of both sources.
- Yh = Y_I: closes only Iphase=0; fails all Iphase=7 rows (a pure round-trip mapping, not a physical one).
- Any Iphase switch: not adopted; it would be fit to the closure, not documented.

## Diagnostic only (not a proposal)
Yh = 1 closes every GPPVA row and every outer-crust row, but it ignores the stored field, has no provider documentation, and does not close SPG inner rows. No remediation contract is proposed.

## Stop conditions hit
C3 and C4 imply incompatible Yh for every SPG Iphase=7 row; the meaning of the raw last field cannot be established from documentation.

## Minimal provider question
To the submitters named on the data sheets (SPG(M1): Scurto; GPPVA(FSU2): Pais): (1) In eos.compo, what is the fourth value of the nuclear quadruple for Iphase=0 and Iphase=7 rows, and why does it equal Ye+Ymu on every Iphase=7 row? (2) Are Aav, Zav cluster values or whole-cell totals (SPG: non-integer, Z/A != Ye; GPPVA: integers, Z/A = Ye)? (3) For SPG, where is the free-neutron (gas) fraction recorded, given Yn=0 in all rows? (4) Is Yn=Yp=0 intentional?
