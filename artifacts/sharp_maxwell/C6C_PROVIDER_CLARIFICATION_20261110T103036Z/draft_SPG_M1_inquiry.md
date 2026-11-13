# DRAFT_NOT_SENT — SPG(M1) technical inquiry

Status: **DRAFT_NOT_SENT**. Do not transmit without project-owner approval.
Proposed To: Luigi Scurto — route status in `provider_contacts.md` (data sheet `lscurto@student.uc.pt`, possibly stale; CFisUC page lists `luigiscurto.pr@gmail.com`). Proposed Cc: Helena Pais (`hpais@uc.pt`, verified). Responsible export person: unconfirmed.

**Subject:** Question on the composition (eos.compo) convention in CompOSE SPG(M1) crust table

Dear Dr. Scurto,

We are using the CompOSE SPG(M1) crust table (`1D/Crust/M1`, data sheet dated 2024-02-26; matching the M1 model of Scurto, Pais & Gulminelli, PRD 109, 103015) and would like to understand the composition convention of `eos.compo` before interpreting it. The data sheet defines Iphase (0 outer crust, 7 inner crust) but not the nuclear-group fields. The attached note gives details and reproducible row-level checks; we are not suggesting an error, only that we cannot tell the intended convention from the public documents.

Each row has the form `iT inb iYq Iphase 4×(pair) 1 (group, Aav, Zav, fourth)`. Examples (file blob `ab79da517271bea1f7de37294e87e019de611699`):
- row 1 (outer crust): `... 10 0.0 11 0.0 0 0.46429 1 0.0 | 1 1 56.000 26.000 1.00000`
- row 218 (inner crust): `... 10 0.0 11 0.0 0 0.06087 1 0.0 | 1 1 82.107 24.387 0.06087`

In all 148 inner-crust rows the fourth value equals Ye + Yμ, Yn = Yp = 0, Aav and Zav are non-integer, and Zav/Aav differs from Ye. Hence we would be grateful if you could clarify:

1. For Iphase = 0 and 7, what exactly is the fourth value of the quadruple?
2. Is it the CompOSE number fraction Y_I (so X_h = Aav·Y_I), a baryon/mass fraction, a charge-related quantity, or a provider-specific one?
3. Are Aav and Zav averages over the nuclear cluster only, or over the whole Wigner–Seitz cell including the gas?
4. How is the free-neutron gas of the CLD cell (part II in the paper) represented in the released data?
5. Why are all neutron and proton pair entries zero, including in the inner crust?
6. Is Yn = Yp = 0 intentional, or are those fields unused placeholders for this product?
7. Could you point us to the export script, normalization convention or a worked example that produced rows 1 and 218?

If someone else prepared the CompOSE export, we would be grateful if you could forward this. Thank you very much for making the tables available.

Kind regards,
[sender name / affiliation — to be filled by project owner]

Attachments (DRAFT): `evidence_summary_onepage.md`, `sample_rows_SPG_M1.csv`, `aggregate_statistics.json`.
