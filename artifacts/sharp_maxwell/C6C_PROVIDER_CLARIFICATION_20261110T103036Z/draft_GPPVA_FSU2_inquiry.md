# DRAFT_NOT_SENT — GPPVA(FSU2) technical inquiry

Status: **DRAFT_NOT_SENT**. Do not transmit without project-owner approval.
Proposed To: Helena Pais (`hpais@uc.pt`; VERIFIED on the data sheet and the CFisUC profile). Responsible export person: unconfirmed.

**Subject:** Question on the composition (eos.compo) convention in CompOSE GPPVA(FSU2) crust table

Dear Prof. Pais,

We are using the CompOSE GPPVA(FSU2) crust table (`1D/Crust/FSU2`, https://compose.obspm.fr/eos/209, sheet dated 2023-04-26) and would like to understand the composition convention of `eos.compo` before interpreting it. The data sheet defines Iphase 0/7/6 and the pasta geometry code in `eos.thermo`, but not the fields of the nuclear-group quadruple. We are not suggesting an error; the attached note gives reproducible checks.

Examples (file blob `fa1ad2b809eb1c95825e2ee6eb2650cf3f248693`):
- row 1 (outer crust): `... 10 0.0 11 0.0 0 0.46429 1 0.0 | 1 1 56. 26. 1.00000`
- row 218 (inner crust, geometry code 1 in `eos.thermo`): `... 10 0.0 11 0.0 0 0.04865 1 0.0 | 1 1 740. 36. 0.04865`

In all 52 inner-crust rows Aav and Zav are integers, Ye = Zav/Aav, Yn = Yp = 0, and the fourth value equals Ye + Yμ (not 1.0). The unified product (https://compose.obspm.fr/eos/215) repeats these rows and, for its core rows (Iphase 6), uses Aav = Zav = fourth value = 0 with Yn, Yp nonzero. Hence we would be grateful if you could clarify:

1. What exactly is the fourth value of the quadruple for Iphase = 0 and 7?
2. Is it Y_I per CompOSE Eq. 4.22, the baryon/mass fraction Aav·Y_I, or another quantity?
3. Are Aav and Zav cluster quantities or whole-cell totals?
4. In the inner-crust rows we find integer Aav/Zav, Ye = Zav/Aav and Yn = Yp = 0, while the fourth value equals Ye + Yμ. Is the intended interpretation that the whole cell is described by the nuclear group?
5. If so, why is the fourth value not 1.0 in these rows? What normalization should be applied to reconstruct total baryon content?
6. Is the `eos.thermo` geometry code (1 = droplets in all Iphase=7 rows) independent of Iphase, and does it affect any composition normalization?
7. Could you share the export schema or a worked row-level example with number densities and mass fractions?

If someone else prepared the CompOSE export, we would be grateful if you could forward this. Thank you very much for the data.

Kind regards,
[sender name / affiliation — to be filled by project owner]

Attachments (DRAFT): `evidence_summary_onepage.md`, `sample_rows_GPPVA_FSU2.csv`, `aggregate_statistics.json`.
