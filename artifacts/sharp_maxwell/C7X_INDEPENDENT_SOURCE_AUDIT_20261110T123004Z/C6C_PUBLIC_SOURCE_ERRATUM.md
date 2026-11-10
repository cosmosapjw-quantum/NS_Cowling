C6C public-source erratum and additive source findings

Date: 2026-11-10. Research track: C7X_INDEPENDENT_SOURCE_AUDIT. Repository parent: 6d8d55013d829b91760cacc4de7c025a959b1a38. Raw-source authority: a09ec30d74b363652d2d1aefe990557c468679dd.

1. CompOSE standard schema versus provider conformance

CompOSE Reference Manual v3.01 defines particle fractions Y_i=n_i/n_B, sum_i B_i Y_i=1, and X_i=B_i Y_i. A nuclear-group quadruple stores (I,A_av,Z_av,Y_I), where A_av and Z_av are abundance-weighted average mass and charge numbers and Y_I is the combined number fraction. The standard baryon/mass fraction is X_h=A_av Y_I. Mapping I to the nuclear set M_I is model-specific and belongs in the EOS Data Sheet. Documented standard meaning does not establish provider conformance. Neither Y_h=Y_I nor Y_h=A_av Y_I is adopted as a production repair for unresolved providers.

Reference: S. Typel et al., CompOSE reference manual v3.01, EPJ A 58, 221 (2022), DOI 10.1140/epja/s10050-022-00847-y, equations (3.10)-(3.12), (4.19)-(4.23); manual PDF https://compose.obspm.fr/download/pdf/manual_v3.00.pdf.

2. GM_GM1 additive raw-source correction

The earlier manifest narrative described all 349 GM_GM1 eos.compo rows as Npairs=0 and Nquad=0. The pinned raw file instead gives: 349 rows; Iphase=1: 48 rows, Npairs=0; Iphase=3: 301 rows, Npairs=2; Nquad=0 on all rows. Rows 49-349 contain neutron ID 10 and proton ID 11, with Yn+Yp=1 at stored precision. No malformed rows were found. This corrects descriptive prose only; canonical manifest and historical C6B/C6C receipts remain immutable.

3. Remaining crust composition and Nadd

Each of the 11 remaining crusts has four particle pairs and one nuclear-group quadruple per row. Standard-schema diagnostic A_av Y_I lies far above unity in substantial parts of these files. This is a provider-conformance blocker, not permission to remap. All 11 crust sources retain Iphase=0 to Iphase=7 between knots 217 and 218; the jump is not smoothed or reinterpreted as a Maxwell interface. Seven GPPVA crust sources have Nadd=1 and source-dependent additional-column values in 0-3. RG_SLY4 and RG_KDE0v1 also have Nadd=1. Meanings remain unresolved.

4. Correct Zenodo identity

The related unified RMF dataset record is https://zenodo.org/records/11084025, DOI 10.5281/zenodo.11084025, associated with Scurto, Pais and Gulminelli, arXiv:2402.15548. The earlier Zenodo 11077759 citation was unrelated. The large eos.txt was not downloaded. The record endpoint was not independently re-fetched here, so a current complete file listing is not claimed. This dataset does not contain the CompOSE eos.compo nuclear-group quadruples and does not resolve provider semantics.

5. Immutability

No provider inquiry was sent. No Dropbox Paper document was edited. No raw EOS file, C3/C4 equation, tolerance, canonical 65-pair manifest entry, historical output, C6B receipt, or C6C receipt was changed.
